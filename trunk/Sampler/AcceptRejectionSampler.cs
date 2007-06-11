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
    /// accept rejection sampling method
    /// </summary>
    public class AcceptRejectionSampler
    {
        #region Fields
        private double c;
        public TargetDistribution ARSTargetDistribution;
        public InstrumentalDistribution ARSInstrumentalDistribution;
        private int N;     // number of MC iterations
        private int n0;    // number of burn-in
        private List<double[]> targetSamples;
        private int totalIteration;
        private double averageAcceptProbability;
        #endregion

        #region Properties
        public int TotalIteration
        {
            get { return totalIteration; }
        }
        public double AverageAcceptProbability
        {
            get { return averageAcceptProbability; }
        }
        #endregion

        #region Constructor
        /// <summary>
        /// c, number of iterations, number of burn-in
        /// </summary>
        /// <param name="inputc"></param>
        /// <param name="inputN"></param>
        /// <param name="inputn0"></param>
        public AcceptRejectionSampler(double inputc, int inputN, int inputn0)
        {
            c = inputc;
            N = inputN;
            n0 = inputn0;
            int run = Convert.ToInt32(N - n0);
            if (n0 >= N)
                throw new ApplicationException("AcceptRejectionSampling::AcceptRejectionSampling(), burn-in length is great than sample length!");
            targetSamples = new List<double[]>(run);
        }

        public AcceptRejectionSampler(double inputc, int inputN)
            : this(inputc, inputN, 0)
        {
        }
        #endregion

        #region Accept Rejection Method
        /// <summary>
        /// do accept rejection sampling
        /// </summary>
        /// <returns></returns>
        public List<double[]> Sample()
        {
            if (ARSInstrumentalDistribution == null
                || ARSTargetDistribution == null)
                throw new ApplicationException("AcceptRrejectionSampling::Sample(), Target distribution or Instrumental distribution is not delegated.");
            
            totalIteration = 0;
            averageAcceptProbability = 0;
            // initial work space
            UniformDistribution U = new UniformDistribution();
            double u;
            double[] Z;
            double acceptProbability;
            
            // MC interation
            bool continueIterate = true;
            int hits = 0;
            while (continueIterate)
            {
                totalIteration++;

                // sample Z from instrumental distribution
                Z = ARSInstrumentalDistribution.Sample();

                // sample u from U(0, 1)
                u = U.Sample();

                // accept probability
                double targetValue =ARSTargetDistribution.Evaluate(Z);
                double instrumentalValue = ARSInstrumentalDistribution.Evaluate(Z);
                acceptProbability = targetValue / c / instrumentalValue;

                averageAcceptProbability += acceptProbability;

                if (u <= acceptProbability)
                {
                    hits++;
                    if (hits > n0)
                        targetSamples.Add(Z);
                    if (targetSamples.Count >= N)
                        continueIterate = false;
                }
            }
            averageAcceptProbability /= (double)totalIteration;
            return targetSamples;
        } // Sample()
        #endregion

        #region Example
        /// <summary>
        /// A-R example
        /// S Chib and E Greenberg, Understanding the Metropolis-Hastings Algorithm, The American Statistician, Vol. 49, No. 4. (Nov., 1995), pp. 327-335
        /// samples are tesitfied by SigmaPlot
        /// </summary>
        public static void Example()
        {
            int N = 6000;
            int n0 = 0;
            string filename1 = "v1.txt";
            FileStream fs1 = new FileStream(filename1, FileMode.Create, FileAccess.Write);
            StreamWriter w1 = new StreamWriter(fs1);

            // target distribution
            //double[,] S = {{1, 0.9},{0.9, 1}};
            double[,] S = { { 2, 0 }, { 0, 1 } };
            double[] m = { 1, 2 };
            PZMath_matrix targetSigma = new PZMath_matrix(S);
            PZMath_vector targetMu = new PZMath_vector(m);

            // instrumental distribution
            double[,] D = {{2, 0},{0, 1}};
            PZMath_matrix instrumentalSigma = new PZMath_matrix(D);
            PZMath_vector instrumentalMu = new PZMath_vector(m);

            double c = System.Math.Sqrt(System.Math.Abs(instrumentalSigma.LUDet()) / System.Math.Abs(targetSigma.LUDet()));

            // prepare workspace
            AcceptRejectionSampler ARSample = new AcceptRejectionSampler(c, N);

            // prepare instrumental distribution
            InstrumentalDistribution instrumentalDistribution = new InstrumentalDistribution(instrumentalMu, instrumentalSigma);

            // prepare target distribution
            TargetDistribution targetDistribution = new TargetDistribution();
            targetDistribution.TargetDistributionFunction = new targetDistribution(AcceptRejectionSampler.ExampleTargetDistributionFunction);

            // assign target distribution and instrumental distribution to the work space
            ARSample.ARSInstrumentalDistribution = instrumentalDistribution;
            ARSample.ARSTargetDistribution = targetDistribution;

            // A-R sample
            List<double[]> samples = ARSample.Sample();

            for (int i = 0; i < N - n0; i++)
            {
                w1.Write(samples[i][0]);
                w1.Write("          ");
                w1.WriteLine(samples[i][1]);
            }

            w1.Close();
            fs1.Close();
        } // Example()

        /// <summary>
        /// example target distribution function
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        private static double ExampleTargetDistributionFunction(double[] x)
        {
            //double[,] S = { { 1, 0.9 }, { 0.9, 1 } };
            double[,] S = { { 2, 0 }, { 0, 1 } };
            double[] m = { 1, 2 };
            PZMath_matrix targetSigma = new PZMath_matrix(S);
            PZMath_vector targetMu = new PZMath_vector(m);

            MultivariateNormalDistribution random = new MultivariateNormalDistribution(targetMu, targetSigma);
            return random.Evaluate(x);
        } // ExampleTargetDistributionFunction()
        #endregion
    }
}
