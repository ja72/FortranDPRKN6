using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using static System.Math;

namespace JA
{
    public class HarmonicOscillatorModel : IModel
    {
        const int width = 10;
        double omega;
        double[] y_Init;
        double[] dy_Init;

        public HarmonicOscillatorModel(double omega, double[] y_Init, double[] dY_Init)
        {
            this.omega=omega;
            this.y_Init = y_Init ?? throw new ArgumentNullException(nameof(y_Init));
            this.dy_Init = dY_Init ?? throw new ArgumentNullException(nameof(dY_Init));
        }
        public bool HasExactSolution
        {
            get => true;
        }
        /// <summary>The system degrees of freedom.</summary>
        public int Dof
        {
            get => 1;
        }
        /// <summary>The acceleration function.</summary>
        public void Acceleration(double t, double[] y, double[] dy, ref double[] ddy)
        {
            ddy[0] = -omega * omega * y[0];
        }

        /// <summary>The exact solution function</summary>
        public void Exact(double t, ref double[] yExact, ref double[] dyExact)
        {
            yExact[0] = Y_Init[0] * Cos(omega * t) + dY_Init[0] / omega * Sin(omega * t);
            dyExact[0] = -Y_Init[0] * omega * Sin(omega * t) + dY_Init[0] * Cos(omega * t);
        }
        public static void Headers()
        {
            Console.WriteLine($"{"step",width} {"t",width} {"y[0]",width} {"dy[0]",width} {"yExact[0]",width} {"dyExact[0]",width} {"h",width} {"absError",width} {"relError",width}");
        }
        public void Update(double h, double t, double[] y, double[] dy, int steps, int rejectedSteps)
        {
            double[] yExact = new double[1];
            double[] dyExact = new double[1];
            Exact(t, ref yExact, ref dyExact);
            double absError = Abs(y[0] - yExact[0]);
            double relError = absError / Abs(yExact[0]);
            Console.WriteLine($"{steps,width} {t,width:F6} {y[0],width:F6} {dy[0],width:F6} {yExact[0],width:F6} {dyExact[0],width:F6} {h,width:F6} {absError,width:F6} {relError,width:P5}");
        }

        /// <summary>Runs the harmonic oscillator for a number of cycles.</summary>
        public OdeSolver Run(int targetCycles, int stepsPerCycle)
        {
            double tEnd = targetCycles * (2.0 * PI) / Omega;
            int steps = targetCycles * stepsPerCycle;
            var sim = new OdeSolver(this, IntegratorType.Dprkn6);
            sim.Run(tEnd, steps);
            return sim; 
        }

        /// <summary>The oscillator frequency in (rad/s)</summary>
        public double Omega
        {
            get => omega;
            set => omega=value;
        }
        /// <summary>Tje inital position state.</summary>
        public double[] Y_Init
        {
            get => y_Init; set => y_Init=value;
        }
        /// <summary>The inital velocity state.</summary>
        public double[] dY_Init
        {
            get => dy_Init; set => dy_Init=value;
        }
    }

    internal class Program
    {
        static void Main(string[] args)
        {

            Console.WriteLine("===========================================================");
            Console.WriteLine("Testing mod_ode with Harmonic Oscillator Example - C#      ");
            Console.WriteLine();

            HarmonicOscillatorModel.Headers();
            HarmonicOscillatorModel model = new HarmonicOscillatorModel(1.0, new double[] { 1.0 }, new double[] { 0.0 });
            OdeSolver sim = model.Run(10, 600);

            double[] yExact = new double[model.Dof];
            double[] dyExact = new double[model.Dof];
            model.Exact(sim.Time, ref yExact, ref dyExact);

            var absError = Abs(sim.Y[0] - yExact[0]);
            var relError = absError / Abs(yExact[0]);

            Console.WriteLine("--- Integration Statistics ---");
            Console.WriteLine($"Successful Steps: {sim.Steps}");
            Console.WriteLine($"Rejected Steps:   {sim.RejectedSteps}");
            Console.WriteLine($"Final Position:   {sim.Y[0]:F9}");
            Console.WriteLine($"Analytical Exact: {yExact[0]:F9}");
            Console.WriteLine($"Absolute Error:   {absError:F9}");
            Console.WriteLine($"Relative Error:   {relError:P5}");
        }
    }
}
