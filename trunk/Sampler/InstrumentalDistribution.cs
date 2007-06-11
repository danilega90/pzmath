using System;
using System.Collections.Generic;
using System.Text;

namespace eee.Sheffield.PZ.Math
{
    public delegate double instrumentalDistribution(double[] x);

    /// <summary>
    /// instrumental distribution
    /// by default, it is a multivariate normal distribution
    /// </summary>
    public class InstrumentalDistribution
    {
        #region Fields
        private MultivariateNormalDistribution random;
        #endregion

        #region Constructor
        public InstrumentalDistribution(PZMath_vector m, PZMath_matrix s)
        {
            random = new MultivariateNormalDistribution(m, s);
        }
        #endregion

        #region sample and evaluate method
        public double[] Sample()
        {
            return random.Sample();
        }

        public double Evaluate(double[] x)
        {
            return random.Evaluate(x);
        }
        #endregion
    }
}
