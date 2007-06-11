// PZ Math Random
// random generator basic class
//
// -- Class PZMath_random
// 09.12.2005


using System;

namespace eee.Sheffield.PZ.Math
{
	/// <summary>
	/// basic random class
	/// </summary>
	
	public class PZMath_random
	{
		public PZMath_random() {}

	} // PZMath_random

	public class PZMath_random_state
	{
		public long x, y;
		public int uTableSize;
		public long[] z; // length of TABLE_SIZE
		public long ln_seed;
	} // PZMath_random_state
}
