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
    /// Poisson Distribution
    /// </summary>
    public class PoissonDistribution : PZRandomUnivariate
    {
        #region Fields
        private double _lamda;
        private ParkMillerUniform random;
        #endregion

        #region Constructor
        /// <summary>
        /// input lamda
        /// </summary>
        /// <param name="l"></param>
        public PoissonDistribution(double lamda) : base()
        {
            if (lamda <= 0)
                throw new ApplicationException("PoissonDistribution::PossionDistribution(), lamda is less or equal to zero!");

            _lamda = lamda;
            //random = new UniformDistribution();
        }

        public PoissonDistribution(double lamda, long seed) : base(seed)
        {
            if (lamda <= 0)
                throw new ApplicationException("PoissonDistribution::PossionDistribution(), lamda is less or equal to zero!");

            _lamda = lamda;
        }
        #endregion

        #region override base class method
        /// <summary>
        /// reset seed and work space
        /// </summary>
        protected override void Reset()
        {
            random = new ParkMillerUniform(_seed);
        } // Reset()

        /// <summary>
        /// sample a random variable
        /// </summary>
        /// <returns></returns>
        public override double Sample()
        {
            return KnuthSampler();
        } // Sample()

        /// <summary>
        /// Evaluate method
        /// </summary>
        /// <param name="k"></param>
        /// <returns></returns>
        public override double Evaluate(double k)
        {
            double P = 1;
            for (int i = 2; i <= k; i++)
                P = P * i;
            double e = System.Math.Exp(-1.0 * _lamda) * System.Math.Pow(_lamda, k);
            return e / P;
        } // Evaluate()
        #endregion

        #region sample methods
        /// <summary>
        /// Knuth sample method
        /// http://en.wikipedia.org/wiki/Poisson_distribution
        /// slow for large lamda value
        /// </summary>
        /// <returns></returns>
        private double KnuthSampler()
        {
            double L = System.Math.Exp(-1.0 * _lamda);
            double k = 0;
            double p = 1;
            do
            {
                k++;
                double u = random.Sample();
                p = p * u;
            }
            while (p >= L);
            return k - 1;
        } // KnuthSampler()
        #endregion

        #region Example
        /// <summary>
        /// Poisson distribution example
        /// </summary>
        public static void Example()
        {
            string filename1 = "v1.txt";
            FileStream fs1 = new FileStream(filename1, FileMode.Create, FileAccess.Write);
            StreamWriter w1 = new StreamWriter(fs1);
            PoissonDistribution random = new PoissonDistribution(8);
            for (int i = 0; i < 10000; i++)
                w1.WriteLine(random.Sample());
            w1.Close();
            fs1.Close();
        } // Example()
        #endregion
    }
}
