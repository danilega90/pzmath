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
    /// multivariate normal distribution
    /// </summary>
    public class MultivariateNormalDistribution
    {
        #region Fields
        private PZMath_vector mu;
        private PZMath_matrix sigma;
        private int dimension;
        private List<ParkMillerNormal> random;
        private double A;
        private PZMath_matrix invSigma;
        private double detSigma;
        #endregion

        #region Property
        private PZMath_vector Mu
        {
            get { return mu; }
            set { mu = value; }
        }
        private PZMath_matrix Sigma
        {
            get { return sigma; }
            set { sigma = value; }
        }
        #endregion

        #region Contructor
        /// <summary>
        /// N(mu, sigma)
        /// </summary>
        /// <param name="l">lower bound</param>
        /// <param name="u">upper bound</param>
        public MultivariateNormalDistribution(PZMath_vector m, PZMath_matrix s)
        {
            dimension = m.Size;
            mu = new PZMath_vector(dimension);
            mu.MemCopyFrom(m);

            sigma = new PZMath_matrix(dimension, dimension);
            sigma.MemCopyFrom(s);

            if (dimension != sigma.RowCount || dimension != sigma.ColumnCount)
                throw new ApplicationException("MultivariateNormalDistribution::MultivariateNormalDistribution(), mu and sigma in different dimension!");

            Init();
        }

        /// <summary>
        /// sigma is a diagonal matrix
        /// </summary>
        /// <param name="m"></param>
        /// <param name="diag"></param>
        public MultivariateNormalDistribution(PZMath_vector m, PZMath_vector diag)
        {
            if (m.Size != diag.Size)
                throw new ApplicationException("MultivariateNormalDistribution::MultivariateNormalDistribution(), mu and sigma in different dimension!");
            dimension = m.Size;

            mu = new PZMath_vector(dimension);
            mu.MemCopyFrom(m);

            sigma = new PZMath_matrix(dimension, dimension);
            for (int i = 0; i < dimension; i ++)
                sigma[i, i] = diag[i];
            
            Init();
        }
        #endregion

        #region init method
        /// <summary>
        /// init random list
        /// calculate detSigma, invSigma, and A
        /// </summary>
        private void Init()
        {
            random = new List<ParkMillerNormal>(dimension);
            for (int i = 0; i < dimension; i++)
                random.Add(new ParkMillerNormal(DateTime.Now.Millisecond
                + DateTime.Now.Second
                + DateTime.Now.Minute + i * 2000));

            detSigma = sigma.LUDet();
            invSigma = sigma.LUInvert();
            A = 1.0 / System.Math.Pow(2 * System.Math.PI, (double)dimension / 2.0)
                / System.Math.Sqrt(detSigma);
        }

        #endregion

        #region sample and evaluate method
        /// <summary>
        /// sample method
        /// </summary>
        /// <returns></returns>
        public double[] Sample()
        {
            if (!sigma.IsDiagonal())
                throw new ApplicationException("MultivariateNormalDistribution::Sample(), Sigma matrix is not diagonal!");

            double[] sample = new double[dimension];
            for (int i = 0; i < dimension; i++)
            {
                double r = random[i].Sample();
                sample[i] = r * System.Math.Sqrt(sigma[i, i]) + mu[i];
            }
            return sample;
        } // Sample()

        /// <summary>
        /// evaluate method
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public double Evaluate(double[] x)
        {
            if (x.Length != dimension)
                throw new ApplicationException("MultivariateNormalDistribution::Evaluate(), Input's dimention does not match!");
            PZMath_vector input = new PZMath_vector(x);
            PZMath_vector diff = input - mu;
            
            PZMath_matrix columnVector = new PZMath_matrix(diff, CBLAS_TRANSPOSE.CblasTrans);
            PZMath_matrix rowVector = new PZMath_matrix(diff);

            PZMath_matrix t = rowVector * invSigma * columnVector;
            double e = System.Math.Exp(-0.5 * t[0, 0]);
            return e * A;
        } // Evaluate()
        #endregion

        #region example codes
        /// <summary>
        /// multivariate normal distribution example
        /// </summary>
        public static void Example()
        {
            // init a bivariate normal distribution
            double[,] D = { { 1.9, 0 }, 
                            { 0, 1.9 } };
            double[] m = { 1, 2 };
            PZMath_matrix sigam1 = new PZMath_matrix(D);
            PZMath_vector mu1 = new PZMath_vector(m);

            MultivariateNormalDistribution BND = new MultivariateNormalDistribution(mu1, sigam1);
            
            // sample
            string filename1 = "v1.txt";
            FileStream fs1 = new FileStream(filename1, FileMode.Create, FileAccess.Write);
            StreamWriter w1 = new StreamWriter(fs1);
            for (int i = 0; i < 10000; i++)
            {
                double[] sample = BND.Sample();
                w1.Write(sample[0]);
                w1.Write("        ");
                w1.WriteLine(sample[1]);
            }
            w1.Close();
            fs1.Close();
            

            // evalualte
            double[] x = { 1.5, 2.5 };
            double value = BND.Evaluate(x);
            // correct answer
            // value = 0.0734
            Console.WriteLine(value);
        }
        #endregion
    }
}
