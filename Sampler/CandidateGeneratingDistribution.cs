// eee.Sheffield.PZ.Math
//
// Copyright ?Ping Zou, 2007
// sg71.cherub@gmail.com

using System;
using System.Collections.Generic;
using System.Text;

namespace eee.Sheffield.PZ.Math
{
    /// <summary>
    /// type of candidate generating distribution
    /// </summary>
    public enum CGDISTRIBUTION { UNIFORM, NORMAL }

    /// <summary>
    /// candidate generating distribution
    /// </summary>
    public class CandidateGeneratingDistribution
    {
        #region Fields
        private CGDISTRIBUTION cgdType;
        private double[] delta;         // uniform distribution bounds
        private PZMath_vector mu;       // normal distribution mu
        private PZMath_matrix sigma;    // normal distribution sigma
        private int dimension;
        private List<NormalDistribution> normalDistributions;
        private List<UniformDistribution> uniformDistributions;
        #endregion

        #region constructor
        /// <summary>
        /// unifor distribution constructor
        /// </summary>
        /// <param name="d"></param>
        public CandidateGeneratingDistribution(double[] d)
        {
            cgdType = CGDISTRIBUTION.UNIFORM;
            dimension = d.Length;
            delta = new double[dimension];

            // init uniform distributions
            uniformDistributions = new List<UniformDistribution>(dimension);
            
            for (int i = 0; i < dimension; i++)
            {
                delta[i] = d[i];
                long seed = DateTime.Now.Millisecond
                    + DateTime.Now.Second
                    + DateTime.Now.Minute + i * 2000;
                uniformDistributions.Add(new UniformDistribution(-1.0 * delta[i], delta[i], seed));
            }
        }

        /// <summary>
        /// normal distribution constructor
        /// </summary>
        /// <param name="m"></param>
        /// <param name="s"></param>
        public CandidateGeneratingDistribution(PZMath_vector m, PZMath_matrix s)
        {
            cgdType = CGDISTRIBUTION.NORMAL;
            dimension = m.Size;
            mu = new PZMath_vector(dimension);
            mu.MemCopyFrom(m);
            sigma = new PZMath_matrix(dimension, dimension);
            sigma.MemCopyFrom(s);

            // init normal distributions
            normalDistributions = new List<NormalDistribution>(dimension);

            for (int i = 0; i < dimension; i++)
            {
                long seed = DateTime.Now.Millisecond
                    + DateTime.Now.Second
                    + DateTime.Now.Minute + i * 2000;
                normalDistributions.Add(new NormalDistribution(mu[i], sigma[i, i], seed));
            }
        }
        #endregion

        #region sample and evalate method
        /// <summary>
        /// random walk sample method
        /// </summary>
        /// <returns></returns>
        public double[] RandomWalkSample(double[] x)
        {
            double[] Y = new double[dimension];
            for (int i = 0; i < dimension; i++)
            {
                if (cgdType == CGDISTRIBUTION.NORMAL)
                    Y[i] = x[i] + normalDistributions[i].Sample();
                else
                    Y[i] = x[i] + uniformDistributions[i].Sample();
            }
            return Y;
        } // Sample()

        /// <summary>
        /// 1st order autogreesive sample
        /// B = I (Identity matrix)
        /// </summary>
        /// <param name="x"></param>
        /// <param name="a"></param>
        /// <returns></returns>
        public double[] FirstOrderAutoGressiveSample(double[] x, double[] a)
        {
            double[] Y = new double[dimension];
            for (int i = 0; i < dimension; i++)
            {
                if (cgdType == CGDISTRIBUTION.NORMAL)
                    Y[i] = a[i] - (x[i] - a[i]) + normalDistributions[i].Sample();
                else
                    Y[i] = a[i] - (x[i] - a[i]) + uniformDistributions[i].Sample();
            }
            return Y;
        }
        /// <summary>
        /// evaluate method, q(x, y) is symmetric, so return 1.0
        /// </summary>
        /// <param name="Z"></param>
        /// <returns></returns>
        public double Evaluate(double[] x, double[] y)
        {
            return 1.0;
        }
        #endregion
    }
}
