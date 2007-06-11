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
    public class UniformDistribution
    {
        #region Fields
        private double lowerBound;
        private double upperBound;
        private PZMath_random_ParkMiller random;
        private double range;
        #endregion

        #region Property
        private double LowerBound
        {
            get { return lowerBound; }
            set { lowerBound = value; }
        }
        private double UpperBound
        {
            get { return upperBound; }
            set { upperBound = value; }
        }
        #endregion

        #region Contructor
        /// <summary>
        /// default constructor U(0, 1]
        /// </summary>
        public UniformDistribution() : this(0.0, 1.0)
        {
        }

        /// <summary>
        /// U(lowerBound, upperBound]
        /// </summary>
        /// <param name="l">lower bound</param>
        /// <param name="u">upper bound</param>
        public UniformDistribution(double l, double u)
        {
            lowerBound = l < u ? l : u;
            upperBound = u > l ? u : l;
            if (l == u)
                throw new ApplicationException("UniformDistribution::UniformDistribution(), lower bound equals upper bound!");
            random = new PZMath_random_ParkMiller();
            range = upperBound - lowerBound;
        }
        public UniformDistribution(long seed)
        {
            lowerBound = 0.0;
            upperBound = 1.0;
            random = new PZMath_random_ParkMiller(seed);
            range = 1.0;
        }

        public UniformDistribution(double l, double u, long seed)
        {
            lowerBound = l < u ? l : u;
            upperBound = u > l ? u : l;
            if (l == u)
                throw new ApplicationException("UniformDistribution::UniformDistribution(), lower bound equals upper bound!");
            random = new PZMath_random_ParkMiller(seed);
            range = upperBound - lowerBound;
        }
        #endregion

        #region generate random numbers
        public double Sample()
        {
            double r = random.NextVariate();
            return r * range + lowerBound;
        }
        #endregion

        #region example codes
        /// <summary>
        /// uniform distribution example
        /// </summary>
        public static void Example()
        {
            string filename1 = "v1.txt";
            FileStream fs1 = new FileStream(filename1, FileMode.Create, FileAccess.Write);
            StreamWriter w1 = new StreamWriter(fs1);
            UniformDistribution random = new UniformDistribution(-0.7, 0.7);
            for (int i = 0; i < 10000; i++)
                w1.WriteLine(random.Sample());
            w1.Close();
            fs1.Close();
        }
        #endregion
    }
}
