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
    public class NormalDistribution : PZRandomUnivariate
    {
        #region Fields
        private double _mu;
        private double _squareSigma;
        private double _sigma;
        private ParkMillerUniform _random;
        #endregion

        #region Property
        private double Mu
        {
            get { return _mu; }
            set { _mu = value; }
        }
        private double SquareSigma
        {
            get { return _squareSigma; }
            set { _squareSigma = value; }
        }
        private double Sigma
        {
            get { return _sigma; }
            set { _sigma = value; }
        }
        #endregion

        #region Contructor
        /// <summary>
        /// default constructor 
        /// N(0, 1)
        /// </summary>
        public NormalDistribution() : this(0.0, 1.0) { }

        /// <summary>
        /// N(mu, squareSigma)
        /// </summary>
        public NormalDistribution(double mu, double sigma) : base()
        {
            _mu = mu;
            _sigma = sigma;
            _squareSigma = sigma * sigma;    
            //random = new ParkMillerNormal();
        }

        public NormalDistribution(long seed) : base(seed)
        {
            _mu = 0.0;
            _squareSigma = 1.0;
            _sigma = 1.0;
            //random = new ParkMillerNormal(seed);
        }
        public NormalDistribution(double mu, double sigma, long seed) : base(seed)
        {
            _mu = mu;
            _squareSigma = sigma * sigma;
            _sigma = sigma;
            //random = new ParkMillerNormal(seed);
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
            return SamplerPolar(_random, _mu, _sigma);
        } // Sample()

        /// <summary>
        /// evaluate method
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public override double Evaluate(double x)
        {
            return Evaluator(x, _mu, _sigma);
        } // Evaluate()
        #endregion

        #region sample/evaluate methods

        /* Of the two methods provided below, I think the Polar method is more
         * efficient, but only when you are actually producing two random
         * deviates.  We don't produce two, because then we'd have to save one
         * in a static variable for the next call, and that would screws up
         * re-entrant or threaded code, so we only produce one.  This makes
         * the Ratio method suddenly more appealing.
         *
         * [Added by Charles Karney] We use Leva's implementation of the Ratio
         * method which avoids calling log() nearly all the time and makes the
         * Ratio method faster than the Polar method (when it produces just one
         * result per call).  Timing per call (gcc -O2 on 866MHz Pentium,
         * average over 10^8 calls)
         *
         *   Polar method: 660 ns
         *   Ratio method: 368 ns
         *
         */
        /* Polar (Box-Mueller) method; See Knuth v2, 3rd ed, p122 */
        public static double SamplerPolar(PZRandomUnivariate random, double mu, double sigma)
        {
            double x, y, r2;
            do
            {
                /* choose x,y in uniform square (-1,-1) to (+1,+1) */
                x = -1 + 2 * random.Sample();
                y = -1 + 2 * random.Sample();

                /* see if it is in the unit circle */
                r2 = x * x + y * y;
            }
            while (r2 > 1.0 || r2 == 0);

            /* Box-Muller transform */
            double r = mu + sigma * y * System.Math.Sqrt(-2.0 * System.Math.Log(r2) / r2);
            return r;
        }  // SamplerPolar ()

        /* Ratio method (Kinderman-Monahan); see Knuth v2, 3rd ed, p130.
         * K+M, ACM Trans Math Software 3 (1977) 257-260.
         *
         * [Added by Charles Karney] This is an implementation of Leva's
         * modifications to the original K+M method; see:
         * J. L. Leva, ACM Trans Math Software 18 (1992) 449-453 and 454-455. */
        public static double SamplerRatio(PZRandomUnivariate random, double mu, double sigma)
        {
            double u, v, x, y, Q;
            double s = 0.449871;    /* Constants from Leva */
            double t = -0.386595;
            double a = 0.19600;
            double b = 0.25472;
            double r1 = 0.27597;
            double r2 = 0.27846;

            do  /* This loop is executed 1.369 times on average  */
            {
                /* Generate a point P = (u, v) uniform in a rectangle enclosing
                   the K+M region v^2 <= - 4 u^2 log(u). */

                /* u in (0, 1] to avoid singularity at u = 0 */
                u = 1 - random.Sample();

                /* v is in the asymmetric interval [-0.5, 0.5).  However v = -0.5
                   is rejected in the last part of the while clause.  The
                   resulting normal deviate is strictly symmetric about 0
                   (provided that v is symmetric once v = -0.5 is excluded). */
                v = random.Sample() - 0.5;

                /* Constant 1.7156 > sqrt(8/e) (for accuracy); but not by too
                   much (for efficiency). */
                v *= 1.7156;

                /* Compute Leva's quadratic form Q */
                x = u - s;
                y = System.Math.Abs(v) - t;
                Q = x * x + y * (a * y - b * x);

                /* Accept P if Q < r1 (Leva) */
                /* Reject P if Q > r2 (Leva) */
                /* Accept if v^2 <= -4 u^2 log(u) (K+M) */
                /* This final test is executed 0.012 times on average. */
            }
            while (Q >= r1 && (Q > r2 || v * v > -4 * u * u * System.Math.Log(u)));
            double r = mu + sigma * (v / u);
            return r;       /* Return slope */
        } // SamplerRatio()


        public static double Evaluator(double x, double mu, double sigma)
        {
            double squareSigma = sigma * sigma;
            double t = (x - mu) * (x - mu) / -2.0 / squareSigma;
            double e = System.Math.Exp(t);
            return e / System.Math.Sqrt(2.0 * System.Math.PI * squareSigma);
        } // Evaluator()

        #endregion

        #region Example
        /// <summary>
        /// normal distribution example
        /// </summary>
        public static void ExampleSample()
        {
            ParkMillerUniform random = new ParkMillerUniform();
            string folderName = System.Environment.GetFolderPath(System.Environment.SpecialFolder.Desktop) + "\\";
            double mu = 1.0;
            double sigma = 2.0;

            // example of Polar method
            string fileName1 = folderName + "Polar method.txt";
            FileStream fs1 = new FileStream(fileName1, FileMode.Create, FileAccess.Write);
            StreamWriter w1 = new StreamWriter(fs1);
            for (int i = 0; i < 10000; i ++)
                w1.WriteLine(NormalDistribution.SamplerPolar(random, mu, sigma));
            w1.Close();
            fs1.Close();

            // computational cost of Polar method
            DateTime startTime1 = DateTime.Now;
            double result1 = 0.0;
            for (int i = 0; i < 10000; i++)
                result1 = NormalDistribution.SamplerPolar(random, mu, sigma);
            DateTime endTime1 = DateTime.Now;
            TimeSpan duration1 = endTime1 - startTime1;
            Console.WriteLine("Polar method running time (n = 10,000): " + duration1);

            // example of Ratio method
            string fileName2 = folderName + "Ratio method.txt";
            FileStream fs2 = new FileStream(fileName2, FileMode.Create, FileAccess.Write);
            StreamWriter w2 = new StreamWriter(fs2);
            for (int i = 0; i < 10000; i ++)
                w2.WriteLine(NormalDistribution.SamplerRatio(random, mu, sigma));
            w2.Close();
            fs2.Close();

            // computational cost of Polar method
            DateTime startTime2 = DateTime.Now;
            double result2 = 0.0;
            for (int i = 0; i < 10000; i++)
                result2 = NormalDistribution.SamplerRatio(random, mu, sigma);
            DateTime endTime2 = DateTime.Now;
            TimeSpan duration2 = endTime2 - startTime1;
            Console.WriteLine("Ratio method running time (n = 10,000): " + duration2);

        } // Example()
        #endregion
    }
}
