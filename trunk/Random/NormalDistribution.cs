// eee.Sheffield.PZ.Math
//
// Copyright ?Ping Zou, 2007
// sg71.cherub@gmail.com

using System;
using System.Collections.Generic;
using System.Text;

namespace eee.Sheffield.PZ.Math
{
    public class NormalDistribution
    {
        #region Fields
        private double mu;
        private double squareSigma;
        private PZMath_random_ParkMiller_Normal random;
        #endregion

        #region Property
        private double Mu
        {
            get { return mu; }
            set { mu = value; }
        }
        private double SquareSigma
        {
            get { return squareSigma; }
            set { squareSigma = value; }
        }
        #endregion

        #region Contructor
        /// <summary>
        /// default constructor 
        /// N(0, 1)
        /// </summary>
        public NormalDistribution() : this(0.0, 1.0)
        {
        }

        /// <summary>
        /// N(mu, squareSigma)
        /// </summary>
        public NormalDistribution(double m, double sigma)
        {
            mu = m;
            squareSigma = sigma * sigma;        
            random = new PZMath_random_ParkMiller_Normal();
        }

        public NormalDistribution(long seed)
        {
            mu = 0.0;
            squareSigma = 1.0;
            random = new PZMath_random_ParkMiller_Normal(seed);
        }
        public NormalDistribution(double m, double sigma, long seed)
        {
            mu = m;
            squareSigma = sigma * sigma;
            random = new PZMath_random_ParkMiller_Normal(seed);
        }

        #endregion

        #region sample and evaluate method
        /// <summary>
        /// sample method
        /// </summary>
        /// <returns></returns>
        public double Sample()
        {
            double r = random.NextVariate();
            return r * System.Math.Sqrt(squareSigma) + mu;
        }

        /// <summary>
        /// evaluate method
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public double Evaluate(double x)
        {
            double t = (x - mu) * (x - mu) / -2.0 / squareSigma;
            double e = System.Math.Exp(t);
            return e / System.Math.Sqrt(2.0 * System.Math.PI * squareSigma);
        }
        #endregion
    }
}
