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
    /// Gamma distribution Gamma(k, theta)
    /// k, shape parameter
    /// theta, scale parameter
    /// </summary>
    public class GammaDistribution : PZRandomUnivariate
    {
        #region Fields
        private double _k;
        private double _theta;
        private ParkMillerUniform random;
        #endregion

        #region Property
        private double K
        {
            get { return _k; }
            set { _k = value; }
        }
        private double Theta
        {
            get { return _theta; }
            set { _theta = value; }
        }
        #endregion

        #region Contructor
        /// <summary>
        /// default constructor 
        /// N(0, 1)
        /// </summary>
        public GammaDistribution() : this(1.0, 2.0) { }

        /// <summary>
        /// N(mu, squareSigma)
        /// </summary>
        public GammaDistribution(double k, double theta)
            : base()
        {
            _k = k;
            _theta = theta;
            //random = new ParkMillerNormal();
        }

        public GammaDistribution(long seed)
            : base(seed)
        {
            _k = 1.0;
            _theta = 2.0;
            //random = new ParkMillerNormal(seed);
        }
        public GammaDistribution(double k, double theta, long seed)
            : base(seed)
        {
            _k = k;
            _theta = theta;
            //random = new ParkMillerNormal(seed);
        }
        #endregion

        #region override base clas method
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
            return Sampler();
        } // Sample()

        /// <summary>
        /// evaluate method
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public override double Evaluate(double x)
        {
            // // TODO: add this function

            throw new Exception("The method or operation is not implemented.");
        }
        #endregion

        #region sample method
        /// <summary>
        /// http://en.wikipedia.org/wiki/Gamma_distribution
        /// </summary>
        /// <returns></returns>
        private double Sampler()
        {
            double integralk = System.Math.Floor(_k);
            double fractionalk = _k - integralk;
            double xi = 0.0;
            
            if (fractionalk > 0)
            {
                // Gamma(delta, 1), delta = fractional of k
                // step 1
                int m = 1;
                double delta = fractionalk;
                double v0 = System.Math.E / (System.Math.E + delta);
                double xim = 0.0;
                double etam = 0.0;

                while (true)
                {
                    // step 2
                    double v2mM1 = random.Sample();
                    double v2m = random.Sample();
                    // step 3
                    if (v2mM1 <= v0)
                    {
                        // step 4
                        xim = System.Math.Pow(v2mM1, 1 / delta);
                        etam = v2m * System.Math.Pow(xim, delta - 1.0);
                    }
                    else
                    {
                        // ste 5
                        xim = 1 - System.Math.Log(v2mM1);
                        etam = v2m * System.Math.Exp(-1.0 * xim);
                    }
                    // step 6
                    if (etam > System.Math.Pow(xim, delta - 1) * System.Math.Exp(-1.0 * xim))
                        m++;
                    // return step 2
                    else
                    {
                        // step 7
                        xi = xim; // xi ~ Gamma(delta, 1)
                        break;
                    }
                }
            }

            // Gamma(n, 1), n = integral of k
            double sum = 0.0;
            double ui;
            for (int i = 0; i < integralk; i++)
            {
                ui = random.Sample();
                sum += System.Math.Log(ui);
            }

            return _theta * (xi - sum);
        } // Sampler()
        #endregion

        #region Example
        /// <summary>
        /// normal distribution example
        /// </summary>
        public static void Example()
        {
            string filename1 = "v1.txt";
            FileStream fs1 = new FileStream(filename1, FileMode.Create, FileAccess.Write);
            StreamWriter w1 = new StreamWriter(fs1);
            GammaDistribution random = new GammaDistribution (3.0, 2.0);
            for (int i = 0; i < 10000; i++)
                w1.WriteLine(random.Sample());
            w1.Close();
            fs1.Close();
        } // Example()
        #endregion
    }
}
