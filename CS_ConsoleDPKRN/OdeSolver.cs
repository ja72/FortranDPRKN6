using System;
using System.Buffers.Binary;

using static System.Math;

namespace JA
{


    public interface IModel
    {
        int Dof { get; }
        void Acceleration(double t, double[] y, double[] dy, ref double[] ddy);
        bool HasExactSolution { get; }
        void Exact(double t, ref double[] yExact, ref double[] dyExact);
        void Update(double h, double t, double[] y, double[] dy, int steps, int rejectedSteps);
        double[] Y_Init { get; set; }
        double[] dY_Init { get; set; }
    }

    public class OdeSolver
    {
        IModel model;
        IIntegrator integrator;
        double time;
        double[] y;
        double[] dy;
        int steps;
        int rejectedSteps;
        double absoluteTolerance;
        double relativeTolerance;
        double lastStepSize;
        const int pollSteps = 400;


        /// <summary>
        /// Initializes a new instance of the <see cref="OdeSolver"/> class.
        /// </summary>
        /// <param name="omega">The angular frequency.</param>
        /// <param name="y_Init">The initial positions.</param>
        /// <param name="dY_Init">The initial velocities.</param>
        public OdeSolver(IModel model, IntegratorType integratorType)
        {
            this.model = model;

            this.y  = new double[model.Dof];
            this.dy = new double[model.Dof];

            switch (integratorType)
            {
                case IntegratorType.Dprkn6:
                    this.integrator = new Dprkn6Integrator(model);
                break;
                default:
                throw new NotSupportedException( $"Integrator type {integratorType} is not supported." );
            }

            this.absoluteTolerance = 1.0e-8;
            this.relativeTolerance = 1.0e-8;

            Reset();
        }

        /// <summary>
        /// Reset the simulation to the initial conditions.
        /// </summary>
        private void Reset()
        {
            this.lastStepSize = 0.0;
            this.steps = 0;
            this.rejectedSteps = 0;
            this.time = 0.0;
            this.model.Y_Init.CopyTo(this.Y, 0);
            this.model.dY_Init.CopyTo(this.dY, 0);
        }

        /// <summary>The current time.</summary>
        public ref double Time => ref time;
        /// <summary>The current position state.</summary>
        public ref double[] Y => ref y;
        /// <summary>The current velocity state.</summary>
        public ref double[] dY => ref dy;
        public double AbsoluteTolerance
        {
            get => absoluteTolerance;
            set => absoluteTolerance = value;
        }
        public double RelativeTolerance
        {
            get => relativeTolerance;
            set => relativeTolerance = value;
        }
        public int Steps => steps;
        public int RejectedSteps => rejectedSteps;

        /// <summary>Runs the simulation for a estimated number of steps.</summary>
        public void Run(double tEnd, int steps)
            => Run(tEnd, tEnd / steps);

        /// <summary>Runs the simulation a up to a fixed end time, with variable steps.</summary>
        public void Run(double tEnd, double h_init)
        {
            const double exp1 = 1.0 / 6.0;
            const double exp2 = 1.0 / 5.0;

            lastStepSize = h_init;

            int n = y.Length;
            double[] yNext;
            double[] dyNext;
            double[] errY;
            double[] errDy;
            double h = h_init;

            steps = 0;
            rejectedSteps = 0;

            model.Update(h, time, y, dy, steps, rejectedSteps);

            while (time < tEnd)
            {
                if (time + h > tEnd)
                    h = tEnd - time;

                integrator.Step(time, h, y, dy, out yNext, out dyNext, out errY, out errDy);

                double totalError = 0.0;
                for (int i = 0; i < n; i++)
                {
                    totalError = Math.Max(
                        totalError,
                        Math.Abs(errY[i]) /
                        (absoluteTolerance + relativeTolerance * Math.Abs(yNext[i])));
                    totalError = Math.Max(
                        totalError,
                        Math.Abs(errDy[i]) /
                        (absoluteTolerance + relativeTolerance * Math.Abs(dyNext[i])));
                }

                if (totalError <= 1.0)
                {
                    time += h;
                    Array.Copy(yNext, y, n);
                    Array.Copy(dyNext, dy, n);
                    steps++;

                    if(steps % pollSteps == 0)
                        model.Update(h, time, y, dy, steps, rejectedSteps);

                    double scale = totalError == 0.0
                        ? 5.0
                        : 0.9 * Math.Pow(1.0 / totalError, exp1);

                    h *= Math.Max(0.2, Math.Min(5.0, scale));
                }
                else
                {
                    rejectedSteps++;
                    double scale = 0.9 * Math.Pow(1.0 / totalError, exp2);

                    h *= Math.Max(0.1, Math.Min(0.9, scale));
                }
                lastStepSize = h_init;
            }
            model.Update(h, time, y, dy, steps, rejectedSteps);
        }

    }

}
