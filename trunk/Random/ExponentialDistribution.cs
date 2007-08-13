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
    /* randist/exponential.c
     * 
     * Copyright (C) 1996, 1997, 1998, 1999, 2000 James Theiler, Brian Gough
     * 
     * This program is free software; you can redistribute it and/or modify
     * it under the terms of the GNU General Public License as published by
     * the Free Software Foundation; either version 2 of the License, or (at
     * your option) any later version.
     * 
     * This program is distributed in the hope that it will be useful, but
     * WITHOUT ANY WARRANTY; without even the implied warranty of
     * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
     * General Public License for more details.
     * 
     * You should have received a copy of the GNU General Public License
     * along with this program; if not, write to the Free Software
     * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
     */
    public class ExponentialDistribution : PZRandomUnivariate
    {
        #region Fields
        private double _mu;
        private ParkMillerUniform _random;
        #endregion

        #region Property
        private double Mu
        {
            get { return _mu; }
            set { _mu = value; }
        }
        #endregion

        #region Contructor
        /// <summary>
        /// default constructor 
        /// Exponential(1.0)
        /// </summary>
        public ExponentialDistribution() : this(1.0) { }

        /// <summary>
        /// Exponential(k)
        /// </summary>
        public ExponentialDistribution(double mu)
            : base()
        {
            _mu = mu;
        }

        public ExponentialDistribution(long seed)
            : base(seed)
        {
            _mu = 1.0;

        }
        public ExponentialDistribution(double mu, long seed)
            : base(seed)
        {
            _mu = mu;
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
            return Sampler(_random, _mu);
        } // Sample()

        /// <summary>
        /// evaluate method
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public override double Evaluate(double x)
        {
            return Evaluator(x, _mu);
        }
        #endregion

        #region sample/evaluate methods
        /* The exponential distribution has the form
        p(x) dx = exp(-x/mu) dx/mu
        for x = 0 ... +infty */
        public static double Sampler(PZRandomUnivariate random, double mu)
        {
            double e = random.Sample();
            double r = -mu * System.Math.Log(e);
            return r;
        } // Sampler()

        public static double Evaluator(double x, double mu)
        {
            if (x < 0)
            {
                return 0;
            }
            else
            {
                double p = System.Math.Exp(-x / mu) / mu;
                return p;
            }            
        } // Evaluator()
        #endregion

        #region example codes
        /// <summary>
        /// exponential distribution example
        /// </summary>
        public static void ExampleSample()
        {
            ParkMillerUniform random = new ParkMillerUniform();
            string folderName = System.Environment.GetFolderPath(System.Environment.SpecialFolder.Desktop) + "\\";
            double mu = 2.0;

            string fileName1 = folderName + "exponential distribution.txt";
            FileStream fs1 = new FileStream(fileName1, FileMode.Create, FileAccess.Write);
            StreamWriter w1 = new StreamWriter(fs1);

            for (int i = 0; i < 10000; i++)
                w1.WriteLine(ExponentialDistribution.Sampler(random, mu));

            w1.Close();
            fs1.Close();

            // computational cost of the method
            DateTime startTime1 = DateTime.Now;
            double result1 = 0.0;
            for (int i = 0; i < 10000; i++)
                result1 = ExponentialDistribution.Sampler(random, mu);
            DateTime endTime1 = DateTime.Now;
            TimeSpan duration1 = endTime1 - startTime1;
            Console.WriteLine("Exponential running time (n = 10,000): " + duration1);
        } // ExampleSample()
        #endregion


    }
}
