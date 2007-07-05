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
        private ParkMillerNormal random;
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
            _squareSigma = sigma * sigma;        
            //random = new ParkMillerNormal();
        }

        public NormalDistribution(long seed) : base(seed)
        {
            _mu = 0.0;
            _squareSigma = 1.0;
            //random = new ParkMillerNormal(seed);
        }
        public NormalDistribution(double mu, double sigma, long seed) : base(seed)
        {
            _mu = mu;
            _squareSigma = sigma * sigma;
            //random = new ParkMillerNormal(seed);
        }
        #endregion

        #region override base clas method
        /// <summary>
        /// reset seed and work space
        /// </summary>
        protected override void Reset()
        {
            random = new ParkMillerNormal(_seed);
        } // Reset()

        /// <summary>
        /// sample a random variable
        /// </summary>
        /// <returns></returns>
        public override double Sample()
        {
            double r = random.Sample();
            return r * System.Math.Sqrt(_squareSigma) + _mu;
        } // Sample()

        /// <summary>
        /// evaluate method
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public override double Evaluate(double x)
        {
            double t = (x - _mu) * (x - _mu) / -2.0 / _squareSigma;
            double e = System.Math.Exp(t);
            return e / System.Math.Sqrt(2.0 * System.Math.PI * _squareSigma);
        } // Evaluate()
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
            NormalDistribution random = new NormalDistribution(1.0, 4);
            for (int i = 0; i < 10000; i++)
                w1.WriteLine(random.Sample());
            w1.Close();
            fs1.Close();
        } // Example()
        #endregion
    }
}
