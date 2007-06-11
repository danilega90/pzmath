// eee.Sheffield.PZ.Math
//
// Copyright ?Ping Zou, 2007
// sg71.cherub@gmail.com

using System;
using System.Collections.Generic;
using System.Text;

namespace eee.Sheffield.PZ.Math.Random
{
    /// <summary>
    /// Poisson Distribution
    /// </summary>
    public class PoissonDistribution
    {
        #region Fields
        private double lamda;
        private UniformDistribution random;
        #endregion

        #region Constructor
        /// <summary>
        /// input lamda
        /// </summary>
        /// <param name="l"></param>
        public PoissonDistribution(double l)
        {
            if (l <= 0)
                throw new ApplicationException("PoissonDistribution::PossionDistribution(), lamda is less or equal to zero!");

            lamda = l;
            random = new UniformDistribution();
        }
        #endregion

        #region sample and evaluate method
        /// <summary>
        /// sample method, Knuth method
        /// http://en.wikipedia.org/wiki/Poisson_distribution
        /// </summary>
        /// <returns></returns>
        public int Sample()
        {
            double L = System.Math.Exp(-1.0 * lamda);
            int k = 0;
            double p = 1;
            do
            {
                k++;
                double u = random.Sample();
                p = p * u;
            }
            while (p >= L);
            return k - 1;
        } // Sample()

        /// <summary>
        /// Evaluate method
        /// </summary>
        /// <param name="k"></param>
        /// <returns></returns>
        public double Evaluate(int k)
        {
            int P = 1;
            for (int i = 2; i <= k; i++)
                P = P * i;
            double e = System.Math.Exp(-1.0 * lamda) * System.Math.Pow(lamda, (double)k);
            return e / (double)P;
        } // Evaluate()
        #endregion
    }
}
