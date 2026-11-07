using System;

namespace JA
{


    /// <summary>
    /// Dormand-Prince RKN(6,4,6) integrator for systems of second-order ODEs:
    ///     y'  = dy
    ///     y'' = f(time, y, dy)
    ///
    /// One call to Step advances the solution by one step of length h.
    /// All calculations use IEEE double precision.
    /// </summary>
    public static class Dprkn6
    {
        #region Constants
        private const int Stages = 6;

        // Stage abscissae.
        private static readonly double[] C =
        {
            0.0,
            0.0870311222683072,
            0.2227182283087289,
            0.5103445105261314,
            0.8118021966205844,
            1.0
        };

        // Coefficients are stored in ordinary C# row-major notation.
        // A[j, s] corresponds to the Fortran coefficient a(j+1, s+1).
        private static readonly double[,] A =
        {
            { 0.0, 0.0037872061902041, 0.0094056263590518, 0.0381643445203362, 0.1118128522770744, 0.1172837260538965 },
            { 0.0, 0.0,                0.0154131557026779, 0.0016641849646879, -0.1264875619727402, 0.0461825129618174 },
            { 0.0, 0.0,                0.0,                0.0903332468351608, 0.2798939912071854, 0.1783510580975610 },
            { 0.0, 0.0,                0.0,                0.0,                0.0643264629472304, 0.1189311654388318 },
            { 0.0, 0.0,                0.0,                0.0,                0.0,                0.0392515374478933 },
            { 0.0, 0.0,                0.0,                0.0,                0.0,                0.0 }
        };

        private static readonly double[,] APrime =
        {
            { 0.0, 0.0870311222683072, 0.0526019342795898, 0.1039860667357774, 0.2736183602568603, 0.1147544074211158 },
            { 0.0, 0.0,                0.1701162940291391, -0.1471018599424707, -0.5731997327318057, 0.0891040854495914 },
            { 0.0, 0.0,                0.0,                0.5534603037328247, 0.8524451044437175, 0.2286704179354032 },
            { 0.0, 0.0,                0.0,                0.0,                0.2589384646518123, 0.4627443831776992 },
            { 0.0, 0.0,                0.0,                0.0,                0.0,                0.1047267060161904 },
            { 0.0, 0.0,                0.0,                0.0,                0.0,                0.0 }
        };

        // Sixth-order integration weights.
        private static readonly double[] B =
        {
            0.1172837260538965,
            0.0461825129618174,
            0.1783510580975610,
            0.1189311654388318,
            0.0392515374478933,
            0.0
        };

        private static readonly double[] BPrime =
        {
            0.1147544074211158,
            0.0891040854495914,
            0.2286704179354032,
            0.4627443831776992,
            0.1047267060161904,
            0.0
        };

        // Fourth-order embedded weights for error estimation.
        private static readonly double[] E =
        {
            0.1221711200155609,
            0.0249767280803513,
            0.1764619729909204,
            0.1384074360341775,
            0.0326462725514081,
            0.0053364703275818
        };

        private static readonly double[] EPrime =
        {
            0.1039860667357774,
            0.1118671569424855,
            0.2104595292415174,
            0.4578132009893992,
            0.0833131707567084,
            0.0325608753341121
        }; 
        #endregion

        /// <summary>
        /// Performs one DPRKN6 integration step.
        ///
        /// Inputs:
        ///   f  - acceleration function y'' = f(time, y, dy)
        ///   time  - current independent variable
        ///   h  - step size
        ///   y  - current position/state vector
        ///   dy - current velocity/first-derivative vector
        ///
        /// Outputs:
        ///   yNext  - position at time + h
        ///   dyNext - velocity at time + h
        ///   errY   - embedded position error estimate
        ///   errDy  - embedded velocity error estimate
        /// </summary>
        public static void Step(
            AccelerationFunction f,
            double t,
            double h,
            double[] y,
            double[] dy,
            out double[] yNext,
            out double[] dyNext,
            out double[] errY,
            out double[] errDy)
        {


            if (f == null)
                throw new ArgumentNullException(nameof(f));
            if (y == null)
                throw new ArgumentNullException(nameof(y));
            if (dy == null)
                throw new ArgumentNullException(nameof(dy));
            if (y.Length != dy.Length)
                throw new ArgumentException("y and dy must have the same length.");

            int n = y.Length;

            // k[stage, component] stores the acceleration at each stage.
            double[,] k = new double[Stages, n];

            double[] yStage = new double[n];
            double[] dyStage = new double[n];
            double[] acceleration = new double[n];

            // Stage 1 evaluation (Fortran stage 1).
            f(t, y, dy, ref  acceleration);
            for (int i = 0; i < n; i++)
                k[0, i] = acceleration[i];

            // Stages 2 through 6.
            for (int s = 1; s < Stages; s++)
            {
                double tStage = t + C[s] * h;

                for (int i = 0; i < n; i++)
                {
                    yStage[i] = y[i] + C[s] * h * dy[i];
                    dyStage[i] = dy[i];
                }

                for (int j = 0; j < s; j++)
                {
                    double h2A = h * h * A[j, s];
                    double hAPrime = h * APrime[j, s];

                    for (int i = 0; i < n; i++)
                    {
                        yStage[i] += h2A * k[j, i];
                        dyStage[i] += hAPrime * k[j, i];
                    }
                }

                Array.Clear(acceleration, 0, n);
                f(tStage, yStage, dyStage, ref acceleration);

                for (int i = 0; i < n; i++)
                    k[s, i] = acceleration[i];
            }

            // Compute the final sixth-order accurate update.
            yNext = new double[n];
            dyNext = new double[n];

            for (int i = 0; i < n; i++)
            {
                yNext[i] = y[i] + h * dy[i];
                dyNext[i] = dy[i];
            }

            double h2 = h * h;

            for (int s = 0; s < Stages; s++)
            {
                double h2B = h2 * B[s];
                double hBPrime = h * BPrime[s];

                for (int i = 0; i < n; i++)
                {
                    yNext[i] += h2B * k[s, i];
                    dyNext[i] += hBPrime * k[s, i];
                }
            }

            // Compute local error estimates from the embedded weights.
            errY = new double[n];
            errDy = new double[n];

            for (int s = 0; s < Stages; s++)
            {
                double h2ErrorWeight = h2 * (B[s] - E[s]);
                double hErrorWeight = h * (BPrime[s] - EPrime[s]);

                for (int i = 0; i < n; i++)
                {
                    errY[i] += h2ErrorWeight * k[s, i];
                    errDy[i] += hErrorWeight * k[s, i];
                }
            }
        }

        /// <summary>
        /// Integrates from the current value of <paramref name="t"/> to
        /// <paramref name="tf"/> using adaptive step-size control.
        /// The accepted final state is written back to <paramref name="y"/>
        /// and <paramref name="dy"/>. The step size is updated through
        /// <paramref name="h"/> for possible reuse by the caller.
        /// </summary>
        public static void RunSimulation(
            AccelerationFunction f,
            ref double h,
            double tf,
            ref double t,
            ref double[] y,
            ref double[] dy, 
            out int steps, 
            out int rejectedSteps,
            UpdateFunction update = null)
        {
            if (f == null)
                throw new ArgumentNullException(nameof(f));
            if (y == null)
                throw new ArgumentNullException(nameof(y));
            if (dy == null)
                throw new ArgumentNullException(nameof(dy));
            if (y.Length != dy.Length)
                throw new ArgumentException("y and dy must have the same length.");
            if (h <= 0.0)
                throw new ArgumentOutOfRangeException(nameof(h), "The step size must be positive.");

            const double absoluteTolerance = 1.0e-8;
            const double relativeTolerance = 1.0e-8;
            const double exp1 = 1.0 / 6.0;
            const double exp2 = 1.0 / 5.0;

            int n = y.Length;
            double[] yNext;
            double[] dyNext;
            double[] errY;
            double[] errDy;

            steps = 0;
            rejectedSteps = 0;

            update?.Invoke(h, t, y, dy, steps, rejectedSteps);

            while (t < tf)
            {
                if (t + h > tf)
                    h = tf - t;

                Step(f, t, h, y, dy, out yNext, out dyNext, out errY, out errDy);

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
                    t += h;
                    Array.Copy(yNext, y, n);
                    Array.Copy(dyNext, dy, n);
                    steps++;

                    update?.Invoke(h, t, y, dy, steps, rejectedSteps);

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
            }
        }
    }
}
