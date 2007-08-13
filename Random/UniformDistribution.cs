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
    /// U(lowerBound, upperBound), uniform distribution
    /// </summary>
    public class UniformDistribution : PZRandomUnivariate
    {
        #region Fields
        private double _lowerBound;
        private double _upperBound;
        private ParkMillerUniform _random;
        private double _range;
        #endregion

        #region Property
        private double LowerBound
        {
            get { return _lowerBound; }
            set { _lowerBound = value; }
        }
        private double UpperBound
        {
            get { return _upperBound; }
            set { _upperBound = value; }
        }
        #endregion

        #region Contructor
        /// <summary>
        /// default constructor U(0, 1]
        /// </summary>
        public UniformDistribution() : this(0.0, 1.0) { }

        /// <summary>
        /// U(lowerBound, upperBound]
        /// </summary>
        /// <param name="l">lower bound</param>
        /// <param name="u">upper bound</param>
        public UniformDistribution(double l, double u) : base()
        {
            if (l == u)
                throw new ApplicationException("UniformDistribution::UniformDistribution(), lower bound equals upper bound!");

            _lowerBound = l < u ? l : u;
            _upperBound = u > l ? u : l;
            
            //_random = new ParkMillerUniform();
            _range = _upperBound - _lowerBound;
        }
        public UniformDistribution(long seed) : base(seed)
        {
            _lowerBound = 0.0;
            _upperBound = 1.0;
            //random = new ParkMillerUniform(seed);
            _range = 1.0;
        }

        public UniformDistribution(double l, double u, long seed) : base(seed)
        {
            if (l == u)
                throw new ApplicationException("UniformDistribution::UniformDistribution(), lower bound equals upper bound!");

            _lowerBound = l < u ? l : u;
            _upperBound = u > l ? u : l;
            
            //random = new ParkMillerUniform(seed);
            _range = _upperBound - _lowerBound;
        }
        #endregion

        #region override base class method
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
            return Sampler(_random) * _range + _lowerBound;
        } // Sample()

        public override double Evaluate(double x)
        {
            return 1.0 / _range;
        } // Evaluate()
        #endregion


        #region sample/evaluate methods
        public static double Sampler(PZRandomUnivariate random)
        {
            double r = random.Sample();
            return r;
        } // Sampler()

        /// <summary>
        /// U(0, range], return int value
        /// </summary>
        /// <param name="random"></param>
        /// <param name="range"></param>
        /// <returns></returns>
        public static int SamplerInt(PZRandomUnivariate random, int range)
        {
            return 0;
        }
        #endregion

        #region example codes
        /// <summary>
        /// uniform distribution example
        /// </summary>
        public static void ExampleSample()
        {
            string folderName = System.Environment.GetFolderPath(System.Environment.SpecialFolder.Desktop) + "\\";
            string fileName1 = "uniform distibution.txt";
            double lowerBound = -0.7;
            double upperBound = 0.7;

            FileStream fs1 = new FileStream(fileName1, FileMode.Create, FileAccess.Write);
            StreamWriter w1 = new StreamWriter(fs1);

            UniformDistribution random = new UniformDistribution(lowerBound, upperBound);

            for (int i = 0; i < 10000; i++)
                w1.WriteLine(random.Sample());

            w1.Close();
            fs1.Close();
        } // ExampleSample()
        #endregion
    }
}
