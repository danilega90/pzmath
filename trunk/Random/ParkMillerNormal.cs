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
    /// normal random generator, using Park Miller uniform random generator, N(0, 1) only
    /// </summary>
    public class ParkMillerNormal : PZRandomUnivariate
    {
        #region Fields
        private bool bReady;
        private double dStored;
        private ParkMillerUniform _ran0;
        private ParkMillerUniform _ran1;
        #endregion

        #region Constructor
        public ParkMillerNormal() : base()
            //: this(DateTime.Now.Millisecond
            //    + DateTime.Now.Second
            //    + DateTime.Now.Minute)
        {
        }

        public ParkMillerNormal(long seed) : base (seed)
        {
            bReady = false;		// Initially, no stored variate avaiable!
            dStored = -PZMath_machine.PZMath_DBL_MAX;	// Set to arbitrary number
            _ran0 = new ParkMillerUniform(seed);
            _ran1 = new ParkMillerUniform(seed + 100);
        } // PZMath_random_ParkMiller_Normal()
        #endregion
        
        #region override base class methods
        /// <summary>
        /// reset seed and work space
        /// </summary>
        protected override void Reset()
        {
            bReady = false;		// Initially, no stored variate avaiable!
            dStored = -PZMath_machine.PZMath_DBL_MAX;	// Set to arbitrary number
            _ran0 = new ParkMillerUniform(_seed);
            _ran1 = new ParkMillerUniform(_seed + 100);
        } // Reset()

        /// <summary>
        /// sample random variable
        /// </summary>
        /// <returns></returns>
        public override double Sample()
        {
            double dFac, dR_2, v1, v2;

            if (bReady)	// Test if variate is available
            {
                bReady = false;
                return dStored;
            }

            // Generate point within unit circle
            do
            {
                v1 = 2.0 * _ran0.Sample() - 1.0;
                v2 = 2.0 * _ran1.Sample() - 1.0;
                dR_2 = v1 * v1 + v2 * v2;
            }
            while (dR_2 >= 1.0 || dR_2 == 0.0);

            // Generate new pair of variates
            dFac = System.Math.Sqrt(-2.0 * System.Math.Log(dR_2) / dR_2);
            dStored = v1 * dFac;
            bReady = true;

            return v2 * dFac;
        } // Sample()

        /// <summary>
        /// evaluate probability, mu = 0, squaresigma = 1.0
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public override double Evaluate(double x)
        {
            double t = x * x / -2.0;
            double e = System.Math.Exp(t);
            return e / System.Math.Sqrt(2.0 * System.Math.PI);
        } // Evaluate()
        #endregion

        //#region sample method
        //public double NextVariate()
        //{
        //    double dFac, dR_2, v1, v2;

        //    if (bReady)	// Test if variate is available
        //    {
        //        bReady = false;
        //        return dStored;
        //    }

        //    // Generate point within unit circle
        //    do
        //    {
        //        v1 = 2.0 * _ran0.Sample() - 1.0;
        //        v2 = 2.0 * _ran1.Sample() - 1.0;
        //        dR_2 = v1 * v1 + v2 * v2;
        //    }
        //    while (dR_2 >= 1.0 || dR_2 == 0.0);

        //    // Generate new pair of variates
        //    dFac = System.Math.Sqrt(-2.0 * System.Math.Log(dR_2) / dR_2);
        //    dStored = v1 * dFac;
        //    bReady = true;

        //    return v2 * dFac;
        //} // NextVariate()
        //#endregion

        //#region I/O methods
        //public void SaveState(string filename)
        //{
        //    FileStream fs = new FileStream(filename, FileMode.OpenOrCreate, FileAccess.Write);
        //    StreamWriter w = new StreamWriter(fs);
        //    w.WriteLine(bReady);
        //    w.WriteLine(dStored);
        //    Ran0.SaveState(w);
        //    Ran1.SaveState(w);
        //    w.Close();
        //    fs.Close();
        //} // SaveState()

        //public void ReadState(string filename)
        //{
        //    FileStream fs = new FileStream(filename, FileMode.Open, FileAccess.Read);
        //    StreamReader r = new StreamReader(fs);
        //    bReady = Convert.ToBoolean(r.ReadLine());
        //    dStored = Convert.ToDouble(r.ReadLine());
        //    Ran0.ReadState(r);
        //    Ran1.ReadState(r);
        //    r.Close();
        //    fs.Close();

        //} // ReadState()
        //#endregion

        #region Example
        /// <summary>
        /// normal samples example
        /// testify samples distribution in SigmaPlot
        /// </summary>
        public static void Example()
        {
            string filename1 = "v1.txt";
            string filename2 = "v2.txt";

            FileStream fs1 = new FileStream(filename1, FileMode.Create, FileAccess.Write);
            StreamWriter w1 = new StreamWriter(fs1);
            ParkMillerNormal random1 = new ParkMillerNormal();
            for (int i = 0; i < 10000; i++)
                w1.WriteLine(random1.Sample());
            w1.Close();
            fs1.Close();

            FileStream fs2 = new FileStream(filename2, FileMode.Create, FileAccess.Write);
            StreamWriter w2 = new StreamWriter(fs2);
            ParkMillerNormal random2 = new ParkMillerNormal();
            for (int i = 0; i < 10000; i++)
                w2.WriteLine(random2.Sample());
            w2.Close();
            fs2.Close();
        }
        #endregion
    } // ParkMillerNormal
}
