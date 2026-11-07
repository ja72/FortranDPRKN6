using System;

namespace JA
{
    public interface IIntegrator
    {

        IModel Model { get; }
        bool IsAdaptive { get; }
        void Step(
            double t,
            double h,
            double[] y,
            double[] dy,
            out double[] yNext,
            out double[] dyNext,
            out double[] errY,
            out double[] errDy);

    }

    /// <summary>
    /// Represents the type of the integrator.
    /// </summary>
    public enum IntegratorType
    {
        Dprkn6
    }

    /// <summary>
    /// Represents the Dormand-Prince RKN(6,4,6) integrator.
    /// </summary>
    public class Dprkn6Integrator : IIntegrator
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

        public Dprkn6Integrator(IModel model)
        {
            this.Model = model;
        }
        public bool IsAdaptive => true;
        public IModel Model { get; }
        public void Step(
            double t,
            double h,
            double[] y,
            double[] dy,
            out double[] yNext,
            out double[] dyNext,
            out double[] errY,
            out double[] errDy)
        {
            if (y == null)
                throw new ArgumentNullException(nameof(y));
            if (dy == null)
                throw new ArgumentNullException(nameof(dy));
            if (y.Length != dy.Length)
                throw new ArgumentException("y and dy must have the same length.");

            int n = y.Length;

            // k[stage, component] stores the temp_ddy at each stage.
            double[,] k = new double[Stages, n];

            double[] yStage = new double[n];
            double[] dyStage = new double[n];
            double[] temp_ddy = new double[n];

            // Stage 1 evaluation (Fortran stage 1).
            Model.Acceleration(t, y, dy, ref  temp_ddy);
            for (int i = 0; i < n; i++)
                k[0, i] = temp_ddy[i];

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

                Array.Clear(temp_ddy, 0, n);
                Model.Acceleration(tStage, yStage, dyStage, ref temp_ddy);

                for (int i = 0; i < n; i++)
                    k[s, i] = temp_ddy[i];
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
    }

}
