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
    /// t distribution, t(nu)
    /// nu: degree of freedom
    /// pdf: t(x, nu) = Gamma((nu + 1) / 2) / Sqrt(nu * pi) / Gamma(nu / 2) * Pow((1 + x * x / nu), -1 * (nu + 1) / 2);
    /// </summary>
    public class StudentTDistribution : PZRandomUnivariate
    {
        #region Fields
        private double _nu;
        private ParkMillerUniform _random;
        #endregion

        #region Property
        private double Nu
        {
            get { return _nu; }
            set { _nu = value; }
        }
        #endregion

        #region Contructor
        /// <summary>
        /// default constructor 
        /// N(0, 1)
        /// </summary>
        public StudentTDistribution() : this(1.0) { }

        /// <summary>
        /// N(mu, squareSigma)
        /// </summary>
        public StudentTDistribution(double nu)
            : base()
        {
            _nu = nu;
        }

        public StudentTDistribution(long seed)
            : base(seed)
        {
            _nu = 1.0;
        }
        public StudentTDistribution(double nu, long seed)
            : base(seed)
        {
            _nu = nu;
        }
        #endregion

        #region override base clas method
        /// <summary>
        /// reset seed and work space
        /// </summary>
        protected override void Reset()
        {
            _random = new ParkMillerUniform(_seed);
        } // Reset()

        /// <summary>
        /// sample a random variable
        /// </summary>
        /// <returns></returns>
        public override double Sample()
        {
            return Sampler(_random, _nu);
        } // Sample()

        /// <summary>
        /// evaluate method
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public override double Evaluate(double x)
        {
            return Evaluator(x, _nu);
        }
        #endregion

        #region sample method

        /* The t-distribution has the form

           p(x) dx = (Gamma((nu + 1)/2)/(sqrt(pi nu) Gamma(nu/2))
           * (1 + (x^2)/nu)^-((nu + 1)/2) dx

           The method used here is the one described in Knuth 
         */
        public static double Sampler(PZRandomUnivariate random, double nu)
        {

            if (nu <= 2)
            {

                double Y1 = NormalDistribution.SamplerPolar(random, 0.0, 1.0);
                double Y2 = ChiSquareDistribution.Sampler(random, nu);

                double t = Y1 / System.Math.Sqrt(Y2 / nu);

                return t;
            }
            else
            {
                double Y1, Y2, Z, t;
                do
                {
                    Y1 = NormalDistribution.SamplerPolar(random, 0.0, 1.0);
                    Y2 = ExponentialDistribution.Sampler(random, 1 / (nu / 2 - 1));

                    Z = Y1 * Y1 / (nu - 2);
                }
                while (1 - Z < 0 || System.Math.Log(-Y2 - Z) > (1 - Z));

                /* Note that there is a typo in Knuth's formula, the line below
                   is taken from the original paper of Marsaglia, Mathematics of
                   Computation, 34 (1980), p 234-256 */

                t = Y1 / System.Math.Sqrt((1 - 2 / nu) * (1 - Z));
                return t;
            }
        }// Sampler()

        public static double Evaluator(double x, double nu)
        {
            double p;

            double lg1 = Gamma.LnGamma(nu / 2);
            double lg2 = Gamma.LnGamma((nu + 1) / 2);

            p = ((System.Math.Exp(lg2 - lg1) / System.Math.Sqrt(PZMath.M_PI * nu))
                 * System.Math.Pow((1 + x * x / nu), -(nu + 1) / 2));
            return p;
        } // Evaluator()
        #endregion

        #region example codes
        /// <summary>
        /// student t distribution example
        /// </summary>
        public static void ExampleSample()
        {
            ParkMillerUniform random = new ParkMillerUniform();
            string folderName = System.Environment.GetFolderPath(System.Environment.SpecialFolder.Desktop) + "\\";
            double nu = 4.0;

            string fileName1 = folderName + "t distribution.txt";
            FileStream fs1 = new FileStream(fileName1, FileMode.Create, FileAccess.Write);
            StreamWriter w1 = new StreamWriter(fs1);

            for (int i = 0; i < 10000; i++)
                w1.WriteLine(StudentTDistribution.Sampler(random, nu));

            w1.Close();
            fs1.Close();

            // computational cost of the method
            DateTime startTime1 = DateTime.Now;
            double result1 = 0.0;
            for (int i = 0; i < 10000; i++)
                result1 = StudentTDistribution.Sampler(random, nu);
            DateTime endTime1 = DateTime.Now;
            TimeSpan duration1 = endTime1 - startTime1;
            Console.WriteLine("t-distibution running time (n = 10,000): " + duration1);
        } // ExampleSample()
        #endregion
    }
}
