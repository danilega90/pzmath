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
    /// Metropolis Hastings Sampling algorithm
    /// </summary>
    public class MetropolisHastingsSampler
    {
        #region Fields
        public CandidateGeneratingDistribution MHCandidateGeneratingDistribution;
        public TargetDistribution MHTargetDistribution;
        private int N;     // number of MC iterations
        private int n0;    // number of burn-in
        private List<double[]> targetSamples;
        private List<double[]> totalSamples;
        private double averageAcceptProbability;
        private int dimension;
        private int hits;
        #endregion

        #region Properties
        public List<double[]> TotalSamples
        {
            get { return totalSamples; }
        }
        public int Hits
        {
            get { return hits; }
        }
        public double AverageAcceptProbability
        {
            get { return averageAcceptProbability; }
        }
        #endregion

        #region Constructor
        /// <summary>
        /// number of iterations, number of burn-in
        /// </summary>
        /// <param name="inputN"></param>
        /// <param name="inputn0"></param>
        public MetropolisHastingsSampler(int inputN, int inputn0)
        {
            N = inputN;
            n0 = inputn0;
            int run = N - n0;
            if (n0 >= N)
                throw new ApplicationException("MetropolisHastingsSampling::MetropolisHastingsSampling(), burn-in length is great than sample length!");
            targetSamples = new List<double[]>(run);
            totalSamples = new List<double[]>(N);
        }

        public MetropolisHastingsSampler(int inputN)
            : this(inputN, 0)
        {
        }
        #endregion

        #region MH algorithm

        /// <summary>
        /// init MCMC
        /// </summary>
        /// <param name="initX"></param>
        public void Init(double[] initX)
        {            

            dimension = initX.Length;
            for (int i = 0; i < N; i++)
                totalSamples.Add(new double[dimension]);
            for (int i = 0; i < dimension; i++)
                totalSamples[0][i] = initX[i];
        } // Init()

        /// <summary>
        /// accept probability (accept rate)
        /// </summary>
        /// <returns></returns>
        private double AcceptProbability(double[] x, double[] y)
        {
            double alpha;
            alpha = MHTargetDistribution.Evaluate(y) * MHCandidateGeneratingDistribution.Evaluate(y, x)
                / MHTargetDistribution.Evaluate(x) / MHCandidateGeneratingDistribution.Evaluate(x, y);
            return alpha < 1.0 ? alpha : 1.0;
        } // AcceptProbability()

        /// <summary>
        /// do accept rejection sampling
        /// </summary>
        /// <returns></returns>
        public List<double[]> Sample()
        {
            if (MHTargetDistribution == null
                || MHCandidateGeneratingDistribution == null)
                throw new ApplicationException("MetropolisHastingsSampling::Sample(), Candidate generating distribution or Instrumental distribution is not delegated.");
            
            averageAcceptProbability = 0;

            // initial work space
            UniformDistribution U = new UniformDistribution();
            double u;
            double acceptProbability;
            
            // MCMC interation
            int n = 0;
            hits = 0;
            int Nm1 = N - 1;

            double[] mu = {1, 2};
            while (n < N)
            {          
                // sample Y from candidate generating distribution
                // random walk sample
                double[] Y = MHCandidateGeneratingDistribution.RandomWalkSample(totalSamples[n]);
                // 1st order autogressive sample
                //double[] Y = MHCandidateGeneratingDistribution.FirstOrderAutoGressiveSample(totalSamples[n], mu);

                // sample u from U(0, 1)
                u = U.Sample();

                // accept probability
                acceptProbability = AcceptProbability(totalSamples[n], Y);
                averageAcceptProbability += acceptProbability;

                if (u <= acceptProbability)
                {
                    hits++;
                    if (n < Nm1)
                        for (int i = 0; i < dimension; i++)
                            totalSamples[n + 1][i] = Y[i];

                    if (n >= n0)
                        targetSamples.Add(Y);
                }
                else
                {
                    if (n < Nm1)
                        for (int i = 0; i < dimension; i ++)                            
                            totalSamples[n + 1][i] = totalSamples[n][i];
                    if (n >= n0)
                        targetSamples.Add(totalSamples[n]);
                }
                n++;
            }
            averageAcceptProbability /= (double)N;
            return targetSamples;
        } // Sample()
        #endregion

        #region Example
        /// <summary>
        /// M-H example
        /// S Chib and E Greenberg, Understanding the Metropolis-Hastings Algorithm, The American Statistician, Vol. 49, No. 4. (Nov., 1995), pp. 327-335
        /// 1. random walk - uniform
        /// 2. random walk - normal
        /// 3. 1st order autogressive - uniform (change MH.Sample() samle method)
        /// samples are tesitfied by SigmaPlot
        /// </summary>
        public static void Example()
        {
            int N = 10000;
            int n0 = 1000;

            string filename1 = "v1.txt";
            FileStream fs1 = new FileStream(filename1, FileMode.Create, FileAccess.Write);
            StreamWriter w1 = new StreamWriter(fs1);

            // target distribution
            //double[,] S = {{1, 0.9},{0.9, 1}};
            double[,] S = { { 2, 0 }, { 0, 1 } };
            double[] m = { 1, 2 };
            PZMath_matrix targetSigma = new PZMath_matrix(S);
            PZMath_vector targetMu = new PZMath_vector(m);
            
            // candidate generating distribution
            // normal
            double[,] D = {{2.5, 0},{0, 1.5}};
            double[] mm = { 0, 0 };
            PZMath_matrix cgdSigma = new PZMath_matrix(D);
            PZMath_vector cgdMu = new PZMath_vector(mm);
            // uniform
            double[] d = { 0.75, 1 };     // for random walk
            //double[] d = { 1, 1 };      // for 1st order autogressive

            // prepare workspace
            MetropolisHastingsSampler MHSample = new MetropolisHastingsSampler(N, n0);

            // prepare instrumental distribution
            // multivariate normal candidate generating distribution
            //CandidateGeneratingDistribution cgDistribution = new CandidateGeneratingDistribution(cgdMu, cgdSigma);
            // multivariate-t candidate generating distribution
            CandidateGeneratingDistribution cgDistribution = new CandidateGeneratingDistribution(d);

            // prepare target distribution
            TargetDistribution targetDistribution = new TargetDistribution();
            targetDistribution.TargetDistributionFunction = new targetDistribution(MetropolisHastingsSampler.ExampleTargetDistributionFunction);

            // assign target distribution and instrumental distribution to the work space
            MHSample.MHCandidateGeneratingDistribution = cgDistribution;
            MHSample.MHTargetDistribution = targetDistribution;


            // M-H sample
            double[] initX = { 1, 2 };
            MHSample.Init(initX);
            List<double[]> samples = MHSample.Sample();

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
            //double[,] S = {{1, 0.9},{0.9, 1}};
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
