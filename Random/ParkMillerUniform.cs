// PZ Math Random 
// Park Miller random generator
// learn from Peter I Rockett

// -- Class PZMath_random_ParkMiller, PZMath_random_ParkMiller_Normal
// 09.12.2005

using System;
using System.IO;

namespace eee.Sheffield.PZ.Math
{
	/// <summary>
	/// Park Miller uniform random class
    /// U(0, 1) only
	/// </summary>
	public class ParkMillerUniform : PZRandomUnivariate
	{
        #region Fields
        private int IA = 16807;
        private double AM = (1.0 / (double)(PZMath_machine.PZMath_LONG_MAX));
        private int IQ = 127773;
        private int IR = 2836;
        private int TABLE_SIZE = 32;
        private long NDIV = (1 + (PZMath_machine.PZMath_LONG_MAX - 1) / 32);
        private double RNMX = (1.0 - PZMath_machine.PZMath_DBL_EPSILON);

        private long X;
        private long Y;
        private long[] Z;
        #endregion

		#region constructor
        /// <summary>
        /// default constructor, using current time flag as random seed
        /// </summary>
        public ParkMillerUniform() : base()
        {

        } // ParkMillerUniform()

        /// <summary>
        /// input random seed
        /// </summary>
        /// <param name="lnSeed"></param>
		public ParkMillerUniform(long seed) : base(seed)
		{
        } // ParkMillerUniform()
		#endregion
        
        #region override base class methods
        /// <summary>
        /// reset seed and work space
        /// </summary>
        protected override void Reset()
        {
            Z = new long[TABLE_SIZE];
            long j;
            long k;

            // initialize seed
            Y = _seed;

            if (Y == 0)	// Prevent zero value of idum
                Y = 1;

            // Load the shuffle table after 8 warm-ups
            for (j = TABLE_SIZE + 7; j >= 0; j--)
            {
                // Implement multiplicative congruential generator with Schrage's algorithm
                k = Y / IQ;
                Y = IA * (Y - k * IQ) - IR * k;
                if (Y < 0)
                    Y += PZMath_machine.PZMath_LONG_MAX;

                if (j < TABLE_SIZE)
                    Z[j] = Y;
            }

            X = Z[0];
        } // Reset()

        /// <summary>
        /// sample a random variable
        /// </summary>
        /// <returns></returns>
        public override double Sample()
        {
            // Implement multiplicative congruential generator with Schrage's algorithm
            long k = Y / IQ;
            Y = IA * (Y - k * IQ) - IR * k;
            if (Y < 0)
                Y += PZMath_machine.PZMath_LONG_MAX;

            // Perform Bays-Durham shuffle to remove low-order serial correlations
            long j = X / NDIV;
            X = Z[j];
            Z[j] = Y;
            double dTemp = AM * (double)X;
            if (dTemp > RNMX)
                dTemp = RNMX;

            return dTemp;
        } // Sample()

        public override double Evaluate(double x)
        {
            return 1.0;
        } // Evaluate()
        #endregion

        #region Example
        /// <summary>
        /// uniform samples example
        /// testify samples distribution in SigmaPlot
        /// </summary>
        public static void Example()
        {
            string filename1 = "v1.txt";
            string filename2 = "v2.txt";

            FileStream fs1 = new FileStream(filename1, FileMode.Create, FileAccess.Write);
            StreamWriter w1 = new StreamWriter(fs1);
            ParkMillerUniform random1 = new ParkMillerUniform();
            for (int i = 0; i < 10000; i++)
                w1.WriteLine(random1.Sample());
            w1.Close();
            fs1.Close();

            FileStream fs2 = new FileStream(filename2, FileMode.Create, FileAccess.Write);
            StreamWriter w2 = new StreamWriter(fs2);
            ParkMillerUniform random2 = new ParkMillerUniform();
            for (int i = 0; i < 10000; i++)
                w2.WriteLine(random2.Sample());
            w2.Close();
            fs2.Close();
        }
        #endregion       
    } // ParkMillerUniform	
}
