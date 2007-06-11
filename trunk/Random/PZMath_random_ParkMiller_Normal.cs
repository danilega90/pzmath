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
    /// normal random generator, using Park Miller uniform random generator
    /// </summary>
    public class PZMath_random_ParkMiller_Normal
    {
        #region Fields
        private bool bReady;
        private double dStored;
        private PZMath_random_ParkMiller Ran0, Ran1;
        #endregion

        #region Constructor
        public PZMath_random_ParkMiller_Normal()
            : this(DateTime.Now.Millisecond
                + DateTime.Now.Second
                + DateTime.Now.Minute)
        {
        }

        public PZMath_random_ParkMiller_Normal(long seed)
        {
            bReady = false;		// Initially, no stored variate avaiable!
            dStored = -PZMath_machine.PZMath_DBL_MAX;	// Set to arbitrary number
            Ran0 = new PZMath_random_ParkMiller(seed);
            Ran1 = new PZMath_random_ParkMiller(seed + 100);
        } // PZMath_random_ParkMiller_Normal()
        #endregion

        #region sample method
        public double NextVariate()
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
                v1 = 2.0 * Ran0.NextVariate() - 1.0;
                v2 = 2.0 * Ran1.NextVariate() - 1.0;
                dR_2 = v1 * v1 + v2 * v2;
            }
            while (dR_2 >= 1.0 || dR_2 == 0.0);

            // Generate new pair of variates
            dFac = System.Math.Sqrt(-2.0 * System.Math.Log(dR_2) / dR_2);
            dStored = v1 * dFac;
            bReady = true;

            return v2 * dFac;
        } // NextVariate()
        #endregion

        #region I/O methods
        public void SaveState(string filename)
        {
            FileStream fs = new FileStream(filename, FileMode.OpenOrCreate, FileAccess.Write);
            StreamWriter w = new StreamWriter(fs);
            w.WriteLine(bReady);
            w.WriteLine(dStored);
            Ran0.SaveState(w);
            Ran1.SaveState(w);
            w.Close();
            fs.Close();
        } // SaveState()

        public void ReadState(string filename)
        {
            FileStream fs = new FileStream(filename, FileMode.Open, FileAccess.Read);
            StreamReader r = new StreamReader(fs);
            bReady = Convert.ToBoolean(r.ReadLine());
            dStored = Convert.ToDouble(r.ReadLine());
            Ran0.ReadState(r);
            Ran1.ReadState(r);
            r.Close();
            fs.Close();

        } // ReadState()
        #endregion

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
            PZMath_random_ParkMiller_Normal random1 = new PZMath_random_ParkMiller_Normal();
            for (int i = 0; i < 10000; i++)
                w1.WriteLine(random1.NextVariate());
            w1.Close();
            fs1.Close();

            FileStream fs2 = new FileStream(filename2, FileMode.Create, FileAccess.Write);
            StreamWriter w2 = new StreamWriter(fs2);
            PZMath_random_ParkMiller_Normal random2 = new PZMath_random_ParkMiller_Normal();
            for (int i = 0; i < 10000; i++)
                w2.WriteLine(random2.NextVariate());
            w2.Close();
            fs2.Close();
        }
        #endregion
    } // PZMath_random_ParkMiller_Normal
}
