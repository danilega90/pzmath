// eee.Sheffield.PZ.Math
//
// Copyright ?Ping Zou, 2007
// sg71.cherub@gmail.com

using System;
using System.Collections.Generic;
using System.Text;
using System.IO;

using eee.Sheffield.PZ.Imaging;

namespace eee.Sheffield.PZ.Math
{
    /// <summary>
    /// Geyer and Moller algorithm for a marked point process
    /// C. J. Geyer and J. Moller, Simulation and likehood inference for spatial point process,
    /// Scandinavian Journal of Statistics, Series B, 21:359 - 373, 1994
    /// </summary>
    public class GeyerMollerSampler
    {
        #region Fields
        public TargetDistribution GMTargetDistribution;
        private int N;     // number of MC iterations
        private int n0;    // number of burn-in
        private List<List<PZPoint>> targetSamples;
        private List<List<PZPoint>> totalSamples;
        private double sumAcceptProbability;
        private double averageAcceptProbability;
        private int hits;
        // reference Possion process
        private double lamda;
        private double width;
        private double height;
        private double S;   // v(S);
        // debug variables
        private double sumnx;
        private double averagenx;
        private double sumsx;
        private double averagesx;
        #endregion

        #region Properties
        public List<List<PZPoint>> TotalSamples
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
        public GeyerMollerSampler(int inputN, int inputn0, double l, double w, double h)
        {
            N = inputN;
            n0 = inputn0;
            int run = N - n0;
            if (n0 >= N)
                throw new ApplicationException("GeyerMollerSampler::GeyerMollerSampler(), burn-in length is great than sample length!");
            targetSamples = new List<List<PZPoint>>(run);
            totalSamples = new List<List<PZPoint>>(N);
            lamda = l;
            width = w;
            height = h;
        }

        public GeyerMollerSampler(int inputN, double l, double w, double h)
            : this(inputN, 0, l, w, h)
        {
        }
        #endregion

        #region GeyerMoller algorithm
        /// <summary>
        /// init MCMC
        /// </summary>
        public void Init()
        {
            // MCMC starts from empty point configuration
            List<PZPoint> emptyList = new List<PZPoint>(0);            
            totalSamples.Add(emptyList);
            // v(S)
            S = LebesgueMeasure();
        } // Init()

        /// <summary>
        /// Lebesgue measure
        /// </summary>
        /// <returns></returns>
        private double LebesgueMeasure()
        {
            return 2.0 * lamda * width * height;
        } // LebesgueMeasure()

        /// <summary>
        /// birth accept probability (birth accept rate)
        /// </summary>
        /// <returns></returns>
        private double BirthAcceptProbability(List<PZPoint> x, List<PZPoint> y)
        {
            // calculate birth accept probability 
            double R;
            R = GMTargetDistribution.Evaluate(y) * S
                / GMTargetDistribution.Evaluate(x) / y.Count;
            return R < 1.0 ? R : 1.0;
        } // BirthAcceptProbability()

        /// <summary>
        /// death accept probability
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        private double DeathAcceptProbability(List<PZPoint> x, List<PZPoint> y)
        {
            double R;
            R = GMTargetDistribution.Evaluate(y) * x.Count
                / GMTargetDistribution.Evaluate(x) / S;
            return R < 1.0 ? R : 1.0;
        } // DeathAcceptProbability()

        /// <summary>
        /// do accept rejection sampling
        /// </summary>
        /// <returns></returns>
        public List<List<PZPoint>> Sample()
        {
            if (GMTargetDistribution == null)
                throw new ApplicationException("GeyerMollerSampler::Sample(), Candidate generating distribution or Instrumental distribution is not delegated.");

            averageAcceptProbability = 0.0;
            sumAcceptProbability = 0.0;
            averagesx = 0.0;
            averagenx = 0.0;
            sumnx = 0.0;
            sumsx = 0.0;

            // initial work space
            long seed = DateTime.Now.Millisecond;
            UniformDistribution U = new UniformDistribution(seed + 10);
            double u;
            UniformDistribution Ux = new UniformDistribution(seed + 20);
            UniformDistribution Uy = new UniformDistribution(seed + 30);
            double ux, uy;
            UniformDistribution Ub = new UniformDistribution(seed + 40);
            UniformDistribution Uremove = new UniformDistribution(seed + 50);
            double acceptProbability;
            double birthProbability;
            double deathProbability;

            // MCMC interation
            int n = 0;
            hits = 0;
            int Nm1 = N - 1;

            // start from empty chain.
            while (n < N)
            {
                List<PZPoint> X = new List<PZPoint>(totalSamples[n]);
                int xLength = X.Count;
                bool death = false;
                bool reject = false;

                // Birth or Death
                birthProbability = Ub.Sample();
                if (xLength > 0)
                    deathProbability = 1 - birthProbability;
                else
                {
                    birthProbability = 1.0;
                    deathProbability = 0.0;
                }

                if (birthProbability> deathProbability)
                // Birth
                {
                    death = false;
                    // creat a new point u in S
                    ux = Ux.Sample() * width;
                    uy = Uy.Sample() * height;
                    PZPoint pointU = new PZPoint(ux, uy);

                    // propose Y
                    List<PZPoint> Y = new List<PZPoint>(X);
                    Y.Add(pointU);

                    // accept probability
                    acceptProbability = BirthAcceptProbability(X, Y);
                    sumAcceptProbability += acceptProbability;
                    u = U.Sample();
                    if (u <= acceptProbability)
                    // accept Y
                    {
                        reject = false;
                        hits++;
                        if (n < Nm1)
                            // add Y to totalSample;
                            totalSamples.Add(Y);

                        if (n >= n0)
                            // add Y to targetSample
                            targetSamples.Add(Y);
                    }
                    else
                    // reject Y
                    {
                        reject = true;
                        if (n < Nm1)
                            totalSamples.Add(X);
                        if (n >= n0)
                            targetSamples.Add(X);
                    }
                }
                else
                // Death
                {
                    death = true;
                    // uniformly remove v from X
                    int removeIndex = Convert.ToInt32(Uremove.Sample() * xLength) - 1;
                    removeIndex = removeIndex > 0 ? removeIndex : 0;
                    // propose Y
                    List<PZPoint> Y = new List<PZPoint>(X);
                    Y.RemoveAt(removeIndex);

                    // accept probability
                    acceptProbability = DeathAcceptProbability(X, Y);
                    sumAcceptProbability += acceptProbability;
                    u = U.Sample();
                    if (u <= acceptProbability)
                    // accept Y
                    {
                        reject = false;
                        hits++;
                        if (n < Nm1)
                            // add Y to totalSample;
                            totalSamples.Add(Y);

                        if (n >= n0)
                            // add Y to targetSample
                            targetSamples.Add(Y);
                    }
                    else
                    // reject Y
                    {
                        reject = true;
                        if (n < Nm1)
                            totalSamples.Add(X);
                        if (n >= n0)
                            targetSamples.Add(X);
                    }
                }
                n++;

                // debug output
                averageAcceptProbability = sumAcceptProbability / (double) n;
                double r = 0.3;
                int nx = totalSamples[n - 1].Count;
                int sx = Convert.ToInt32(Sx(totalSamples[n - 1], r));
                sumsx += sx;
                sumnx += nx;
                double hx = GMTargetDistribution.Evaluate(totalSamples[n - 1]);

                Console.Write("n = " + n);
                if (death)
                    Console.Write(" Death");
                else
                    Console.Write(" Birth");
                if (reject)
                    Console.Write(" Reject");
                else
                    Console.Write(" Accept");
                Console.Write(" n(x) = " + nx);
                Console.Write(" s(x) = " + sx);
                Console.Write(" h(x) = " + String.Format("{0:f}", hx));
                Console.Write(" accept rate = " + String.Format("{0:f}", acceptProbability));
                Console.WriteLine(" average accept rate " + String.Format("{0:f}", averageAcceptProbability));
            }
            averageAcceptProbability = sumAcceptProbability / (double) N;
            averagesx = sumsx / (double)N;
            averagenx = sumnx / (double)N;
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

            // file IO
            string filename1 = "v1.txt";
            FileStream fs1 = new FileStream(filename1, FileMode.Create, FileAccess.Write);
            StreamWriter w1 = new StreamWriter(fs1);

            // target distribution
            double lamda = 20;
            double r = 0.3;
            double width = 1.0;
            double height = 1.0;
            double theta1 = 2;
            double theta2 = -1.0;
            
            // prepare workspace
            GeyerMollerSampler GMSampler = new GeyerMollerSampler(N, n0, lamda, width, height);

            // prepare target distribution
            TargetDistribution targetDistribution = new TargetDistribution();
            targetDistribution.TargetPointProcessDistributionFunction = new targetPointProcessDistribution(GeyerMollerSampler.ExampleTargetDistributionFunction);

            // assign target distribution and instrumental distribution to the work space
            GMSampler.GMTargetDistribution = targetDistribution;

            // M-H sample
            GMSampler.Init();

            List<List<PZPoint>> samples = GMSampler.Sample();

            // output samples
            int sampleLength = samples.Count;
            for (int i = 0; i < sampleLength; i++)
            {
                int xLength = samples[i].Count;
                for (int j = 0; j < xLength; j++)
                {
                    w1.Write(String.Format("{0:f}",samples[i][j].x) + "  " + String.Format("{0:f}",samples[i][j].y) + " | ");
                }
                w1.WriteLine();
            }

            w1.Close();
            fs1.Close();
        } // Example()

        /// <summary>
        /// example target distribution function
        /// </summary>
        /// <returns></returns>
        private static double ExampleTargetDistributionFunction(List<PZPoint> x)
        {
            // target Strauss Process parameters
            double r = 0.3;
            double theta1 = 0;
            double theta2 = 0;

            double nx = x.Count;
            double sx = Sx(x, r);

            double hx = System.Math.Exp(1.0 * (nx * theta1 + sx * theta2));
            return hx;
        } // ExampleTargetDistributionFunction()

        /// <summary>
        /// Strauss Proces, s(x)
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        private static double Sx(List<PZPoint> x, double r)
        {
            int xLength = x.Count;
            int Sx = 0;
            for (int i = 0; i < xLength; i++)
            {
                for (int j = 0; j < xLength; j++)
                {
                    if (i == j)
                        continue;
                    else
                    {
                        double dist = x[i].Distance(x[j]);
                        if (dist <= r)
                            Sx++;
                    }
                }
            }

            return (double) Sx / 2.0;
        } // Sx()
        #endregion
    }
}
