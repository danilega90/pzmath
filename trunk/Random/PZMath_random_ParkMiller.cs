// PZ Math Random 
// Park Miller random generator
// learn from Peter I Rockett

// -- Class PZMath_random_ParkMiller, PZMath_random_ParkMiller_Normal
// 09.12.2005

using System;
using System.IO;

namespace eee.Sheffield.PZ.Math
{
	#region PZMath_random_ParkMiller
	/// <summary>
	/// Park Miller uniform random class
	/// </summary>
	public class PZMath_random_ParkMiller
	{
        private int IA = 16807;
        private double AM = (1.0 / (double)(PZMath_machine.PZMath_LONG_MAX));
        private int IQ = 127773;
        private int IR = 2836;
        private int TABLE_SIZE = 32;
        private long NDIV = (1 + (PZMath_machine.PZMath_LONG_MAX - 1) / 32);
        private double RNMX = (1.0 - PZMath_machine.PZMath_DBL_EPSILON);

		private long X;
		private long Y;
		private long [] Z;
		private long lnSeed = 12345;

		#region constructor
        /// <summary>
        /// default constructor, using current time flag as random seed
        /// </summary>
        public PZMath_random_ParkMiller()
            : this(DateTime.Now.Millisecond
                + DateTime.Now.Second
                + DateTime.Now.Minute)
        {
        }

        /// <summary>
        /// input random seed
        /// </summary>
        /// <param name="lnSeed"></param>
		public PZMath_random_ParkMiller(long lnSeed)
		{
			Z = new long [TABLE_SIZE];
			long j;
			long k;

            // initialize seed
			Y = lnSeed;

			if(Y == 0)	// Prevent zero value of idum
				Y = 1;
	
			// Load the shuffle table after 8 warm-ups
			for(j = TABLE_SIZE + 7; j >= 0; j--)
			{
				// Implement multiplicative congruential generator with Schrage's algorithm
				k = Y / IQ;
				Y = IA * (Y - k * IQ) - IR * k;
				if(Y < 0)
					Y += PZMath_machine.PZMath_LONG_MAX;

				if(j < TABLE_SIZE)
					Z[j] = Y;
			}

			X = Z[0];
		} // PZMath_random_ParkMiller()
		#endregion

		#region method
		public double NextVariate()
		{
			// Implement multiplicative congruential generator with Schrage's algorithm
			long k = Y / IQ;
			Y = IA * (Y - k * IQ) - IR * k;
			if(Y < 0)
				Y += PZMath_machine.PZMath_LONG_MAX;

			// Perform Bays-Durham shuffle to remove low-order serial correlations
			long j = X / NDIV;
			X = Z[j];
			Z[j] = Y;
			double dTemp = AM * (double) X;
			if(dTemp > RNMX)
				dTemp = RNMX;

			return dTemp;
		} // NextVariate()
		#endregion

		#region save state file IO
		public void SaveState(string filename)
		{
			FileStream fs = new FileStream(filename, FileMode.Open, FileAccess.Write);
			StreamWriter w = new StreamWriter(fs);
			w.WriteLine(X);
			w.WriteLine(Y);
			w.WriteLine(TABLE_SIZE);
			for (int i = 0; i < TABLE_SIZE; i++)
				w.WriteLine(Z[i]);
			w.WriteLine(lnSeed);
			w.Close();
			fs.Close();
		} // SaveState()

		public void ReadState(string filename)
		{
			FileStream fs = new FileStream(filename, FileMode.Open, FileAccess.Read);
			StreamReader r = new StreamReader(fs);
			X = Convert.ToInt64 (r.ReadLine());
			Y = Convert.ToInt64 (r.ReadLine());
			TABLE_SIZE = Convert.ToInt32 (r.ReadLine());
			for (int i = 0; i < TABLE_SIZE; i ++)
				Z[i] = Convert.ToInt64 (r.ReadLine());
			lnSeed = Convert.ToInt64 (r.ReadLine());
			r.Close();
			fs.Close();
		} // ReadState()

        public void SaveState(StreamWriter w)
        {
            w.WriteLine(X);
            w.WriteLine(Y);
            w.WriteLine(TABLE_SIZE);
            for (int i = 0; i < TABLE_SIZE; i++)
                w.WriteLine(Z[i]);
            w.WriteLine(lnSeed);
        } // SaveState()

        public void ReadState(StreamReader r)
        {        
            X = Convert.ToInt64(r.ReadLine());
            Y = Convert.ToInt64(r.ReadLine());
            TABLE_SIZE = Convert.ToInt32(r.ReadLine());
            for (int i = 0; i < TABLE_SIZE; i++)
                Z[i] = Convert.ToInt64(r.ReadLine());
            lnSeed = Convert.ToInt64(r.ReadLine());
        } // ReadState()
		#endregion

		#region access state
		public PZMath_random_state GetState()
		{
			PZMath_random_state state = new PZMath_random_state();
			
			state.x = X;
			state.y = Y;
			state.uTableSize = TABLE_SIZE;
			state.z = new long [TABLE_SIZE];
			for(int i = 0; i < TABLE_SIZE; i++)
				state.z[i] = Z[i];
			state.ln_seed = lnSeed;

			return state;
		} // GetState()
		public void SetState(PZMath_random_state state)
		{
			X = state.x;
			Y = state.y;

			if(state.uTableSize != TABLE_SIZE)
				PZMath_errno.ERROR ("PZMath_random_ParkMiller::SetState() : Incompatible table sizes!", PZMath_errno.PZMath_EFAILED);

			for(int i = 0; i < TABLE_SIZE; i++)
				Z[i] = state.z[i];

			lnSeed = state.ln_seed;
		} // SetState()
		# endregion

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
            PZMath_random_ParkMiller random1 = new PZMath_random_ParkMiller();
            for (int i = 0; i < 10000; i++)
                w1.WriteLine(random1.NextVariate());
            w1.Close();
            fs1.Close();

            FileStream fs2 = new FileStream(filename2, FileMode.Create, FileAccess.Write);
            StreamWriter w2 = new StreamWriter(fs2);
            PZMath_random_ParkMiller random2 = new PZMath_random_ParkMiller();
            for (int i = 0; i < 10000; i++)
                w2.WriteLine(random2.NextVariate());
            w2.Close();
            fs2.Close();
        }
        #endregion

    } // PZMath_random_ParkMiller
	
	#endregion
}
