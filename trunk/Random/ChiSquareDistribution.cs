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
    /// Chi-square distribution Chi2(k)
    /// k, degree of freedom
    /// pdf: Chi2(x; k) = Pow(x, k / 2 - 1) * Exp(-x / 2) / Pow(2, k / 2) / Gamma(k / 2)
    /// </summary>
    public class ChiSquareDistribution : PZRandomUnivariate
    {
        #region Fields
        private double _k;
        private ParkMillerUniform _random;
        #endregion

        #region Property
        private double K
        {
            get { return _k; }
            set { _k = value; }
        }
        #endregion

        #region Contructor
        /// <summary>
        /// default constructor 
        /// Chi2(1.0)
        /// </summary>
        public ChiSquareDistribution() : this(1.0) { }

        /// <summary>
        /// Chi2(k)
        /// </summary>
        public ChiSquareDistribution(double k)
            : base()
        {
            _k = k;
        }

        public ChiSquareDistribution(long seed)
            : base(seed)
        {
            _k = 1.0;

        }
        public ChiSquareDistribution(double k, long seed)
            : base(seed)
        {
            _k = k;
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
            return Sampler(_random, _k);
        } // Sample()

        /// <summary>
        /// evaluate method
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public override double Evaluate(double x)
        {
            return Evaluator(x, _k);
        }
        #endregion

        #region sample method
        /// <summary>
        /// </summary>
        /// <returns></returns>
        public static double Sampler(PZRandomUnivariate random, double k)
        {
            double chisq = 2 * GammaDistribution.SampleMarsagliaAndTsang(random, k / 2, 1.0);
            return chisq;
        } // Sampler()

        /// <summary>
        /// pdf: Chi2(x; k) = Pow(x, k / 2 - 1) * Exp(-x / 2) / Pow(2, k / 2) / Gamma(k / 2)
        /// </summary>
        /// <param name="x"></param>
        /// <param name="k"></param>
        /// <returns></returns>
        public static double Evaluator(double x, double k)
        {
            if (x <= 0)
            {
                return 0;
            }
            else
            {
                double p;                
                double lngamma = Gamma.LnGamma(k / 2);

                p = System.Math.Exp((k / 2 - 1) * System.Math.Log(x / 2) - x / 2 - lngamma) / 2;
                return p;
            }
        } // Evaluator()
        #endregion

        #region Example
        /// <summary>
        /// ChiSquare distribution example
        /// </summary>
        public static void ExampleSample()
        {
            ParkMillerUniform random = new ParkMillerUniform();
            string folderName = System.Environment.GetFolderPath(System.Environment.SpecialFolder.Desktop) + "\\";
            double k = 4.0;

            string fileName1 = folderName + "ChiSquare distribution.txt";
            FileStream fs1 = new FileStream(fileName1, FileMode.Create, FileAccess.Write);
            StreamWriter w1 = new StreamWriter(fs1);

            for (int i = 0; i < 10000; i++)
                w1.WriteLine(ChiSquareDistribution.Sampler(random, k));

            w1.Close();
            fs1.Close();

            // computational cost of the method
            DateTime startTime1 = DateTime.Now;
            double result1 = 0.0;
            for (int i = 0; i < 10000; i++)
                result1 = ChiSquareDistribution.Sampler(random, k);
            DateTime endTime1 = DateTime.Now;
            TimeSpan duration1 = endTime1 - startTime1;
            Console.WriteLine("ChiSquare distribution running time (n = 10,000): " + duration1);
        } // Example()
        #endregion
    }
}
