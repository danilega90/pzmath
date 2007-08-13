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
	
	public abstract class PZRandomUnivariate : IPZRandomUnivariate
    {
        #region Fields
        protected long _seed;
        #endregion

        #region Properties
        public long Seed { get { return _seed; } }
        #endregion

        #region Constructor
        public PZRandomUnivariate()
        {
            _seed = (long)(DateTime.Now.Millisecond
               + DateTime.Now.Second
               + DateTime.Now.Minute);
            Reset();
        } // PZRandomUnivariate()

        public PZRandomUnivariate(long seed)
        {
            _seed = seed;
            Reset();
        }
        #endregion

        #region Interface method
        public void ResetSeed(long seed)
        {
            _seed = seed;
            Reset();
        }
        protected abstract void Reset();
        public abstract double Sample();
        public abstract double Evaluate(double x);
        #endregion
    } // PZMath_random
}
