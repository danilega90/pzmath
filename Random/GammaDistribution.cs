// eee.Sheffield.PZ.Math
//
// Copyright ?Ping Zou, 2007
// sg71.cherub@gmail.com

using System;
using System.Collections.Generic;
using System.Text;
using System.IO;

namespace eee.Sheffield.PZ.Math
{
    /// <summary>
    /// Gamma distribution Gamma(k, theta)
    /// k, shape parameter
    /// theta, scale parameter
    /// </summary>
    public class GammaDistribution : PZRandomUnivariate
    {
        #region Fields
        private double _k;
        private double _theta;
        private double _alpha;
        private double _beta;
        private ParkMillerUniform _random;
        #endregion

        #region Property
        public double K
        {
            get { return _k; }
            set { _k = value; }
        }
        public double Theta
        {
            get { return _theta; }
            set { _theta = value; }
        }
        public double Alpha { get { return _alpha; } set { _alpha = value; } }
        public double Beta { get { return _beta; } set { _beta = value; } }
        #endregion

        #region Contructor
        /// <summary>
        /// default constructor 
        /// N(0, 1)
        /// </summary>
        public GammaDistribution() : this(1.0, 2.0) { }

        /// <summary>
        /// Gamma(k, theta)
        /// alpha = k, beta = 1 / theta;
        /// </summary>
        public GammaDistribution(double k, double theta)
            : base()
        {
            _k = k;
            _theta = theta;
        }

        public GammaDistribution(long seed)
            : base(seed)
        {
            _k = 1.0;
            _theta = 2.0;
        }
        public GammaDistribution(double k, double theta, long seed)
            : base(seed)
        {
            _k = k;
            _theta = theta;
        }
        #endregion

        #region override base clas method
        /// <summary>
        /// reset seed and work space
        /// </summary>
        protected override void Reset()
        {
            _alpha = _k;
            _beta = 1 / _theta;
            _random = new ParkMillerUniform(_seed);
        } // Reset()

        /// <summary>
        /// sample a random variable
        /// </summary>
        /// <returns></returns>
        public override double Sample()
        {
            return SampleMarsagliaAndTsang(_random, _alpha, _beta);
        } // Sample()

        /// <summary>
        /// evaluate method
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public override double Evaluate(double x)
        {
            return Evaluator(x, _alpha, _beta);
        } // Evaluate()
        #endregion

        #region sample methods
        /// <summary>
        /// http://en.wikipedia.org/wiki/Gamma_distribution
        /// Knuth
        /// </summary>
        /// <returns></returns>
        public static double SampleKnuth(PZRandomUnivariate random, double k, double theta)
        {
            double integralk = System.Math.Floor(k);
            double fractionalk = k - integralk;
            double xi = 0.0;
            
            if (fractionalk > 0)
            {
                // Gamma(delta, 1), delta = fractional of k
                // step 1
                int m = 1;
                double delta = fractionalk;
                double v0 = System.Math.E / (System.Math.E + delta);
                double xim = 0.0;
                double etam = 0.0;

                while (true)
                {
                    // step 2
                    double v2mM1 = random.Sample();
                    double v2m = random.Sample();
                    // step 3
                    if (v2mM1 <= v0)
                    {
                        // step 4
                        xim = System.Math.Pow(v2mM1, 1 / delta);
                        etam = v2m * System.Math.Pow(xim, delta - 1.0);
                    }
                    else
                    {
                        // ste 5
                        xim = 1 - System.Math.Log(v2mM1);
                        etam = v2m * System.Math.Exp(-1.0 * xim);
                    }
                    // step 6
                    if (etam > System.Math.Pow(xim, delta - 1) * System.Math.Exp(-1.0 * xim))
                        m++;
                    // return step 2
                    else
                    {
                        // step 7
                        xi = xim; // xi ~ Gamma(delta, 1)
                        break;
                    }
                }
            }

            // Gamma(n, 1), n = integral of k
            double sum = 0.0;
            double ui;
            for (int i = 0; i < integralk; i++)
            {
                ui = random.Sample();
                sum += System.Math.Log(ui);
            }

            return theta * (xi - sum);
        } // SampleKnuth()

        /// <summary>
        /// New version based on Marsaglia and Tsang, "A Simple Method for
        /// generating gamma variables", ACM Transactions on Mathematical
        /// Software, Vol 26, No 3 (2000), p363-372.
        /// Implemented by J.D.Lamb@btinternet.com, minor modifications for GSL
        /// by Brian Gough
        /// </summary>
        /// <returns></returns>
        public static double SampleMarsagliaAndTsang(PZRandomUnivariate random, double alpha, double beta)
        {
            if (alpha < 1)
            {
                double u = UniformDistribution.Sampler(random);

                return SampleMarsagliaAndTsang(random, 1.0 + alpha, beta) * System.Math.Pow(u, 1.0 / alpha);
            }

            {
                double x, v, u;
                double d = alpha - 1.0 / 3.0;
                double c = (1.0 / 3.0) / System.Math.Sqrt(d);

                while (true)
                {
                    do
                    {
                        x = NormalDistribution.SamplerPolar(random, 0.0, 1.0);
                        v = 1.0 + c * x;
                    }
                    while (v <= 0);

                    v = v * v * v;
                    u = random.Sample();

                    if (u < 1 - 0.0331 * x * x * x * x)
                        break;

                    if (System.Math.Log(u) < 0.5 * x * x + d * (1 - v + System.Math.Log(v)))
                        break;
                }

                return beta * d * v;
            }
        } // SampleMarsagliaAndTsang()

        public static double Evaluator(double x, double alpha, double beta)
        {
            if (x < 0)
            {
                return 0;
            }
            else if (x == 0)
            {
                if (alpha == 1)
                    return 1 / beta;
                else
                    return 0;
            }
            else if (alpha == 1)
            {
                return System.Math.Exp(-x / beta) / beta;
            }
            else
            {
                double p;
                double lngamma = Gamma.LnGamma(alpha);
                p = System.Math.Exp((alpha - 1) * System.Math.Log(x / beta) - x / beta - lngamma) / beta;
                return p;
            }
        } // Evaluator()
        #endregion

        #region Example
        /// <summary>
        /// Gamma distribution example
        /// </summary>
        public static void ExampleSample()
        {
            ParkMillerUniform random = new ParkMillerUniform();
            string folderName = System.Environment.GetFolderPath(System.Environment.SpecialFolder.Desktop) + "\\";
            double k = 2.0;
            double theta = 2.0;
            double alpha = k;
            double beta = 1 / theta;

            // Knuth method
            string fileName1 = folderName + "Gamma distribution Knuth method.txt";
            FileStream fs1 = new FileStream(fileName1, FileMode.Create, FileAccess.Write);
            StreamWriter w1 = new StreamWriter(fs1);

            for (int i = 0; i < 10000; i++)
                w1.WriteLine(GammaDistribution.SampleKnuth(random, k, theta));

            w1.Close();
            fs1.Close();

            // computational cost of Kunth method
            DateTime startTime1 = DateTime.Now;
            double result1 = 0.0;
            for (int i = 0; i < 10000; i++)
                result1 = GammaDistribution.SampleKnuth(random, k, theta);
            DateTime endTime1 = DateTime.Now;
            TimeSpan duration1 = endTime1 - startTime1;
            Console.WriteLine("ChiSquare distribution Knuth method running time (n = 10,000): " + duration1);

            // MT method
            string fileName2 = folderName + "Gamma distribution MT method.txt";
            FileStream fs2 = new FileStream(fileName2, FileMode.Create, FileAccess.Write);
            StreamWriter w2 = new StreamWriter(fs2);

            for (int i = 0; i < 10000; i++)
                w2.WriteLine(GammaDistribution.SampleMarsagliaAndTsang(random, alpha, beta));

            w2.Close();
            fs2.Close();

            // computational cost of the method
            DateTime startTime2 = DateTime.Now;
            double result2 = 0.0;
            for (int i = 0; i < 10000; i++)
                result2 = GammaDistribution.SampleKnuth(random, alpha, beta);
            DateTime endTime2 = DateTime.Now;
            TimeSpan duration2 = endTime2 - startTime2;
            Console.WriteLine("ChiSquare distribution Knuth method running time (n = 10,000): " + duration2);
        } // Example()
        #endregion
    }
}
