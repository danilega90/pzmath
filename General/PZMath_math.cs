// PZ Math
// Based on Math.NET Open source
//
// -- basic declarations
//
// -- Class PZMath_f, Class PZMath_fdf
// 23.11.2005


using System;

namespace eee.Sheffield.PZ.Math
{
	/// <summary>
	/// basic math functions
	/// </summary>
	public class PZMath
	{
        #region static constants
        public static double M_E = 2.71828182845904523536028747135;      /* e */
        public static double M_LOG2E = 1.44269504088896340735992468100;      /* log_2 (e) */
        public static double M_LOG10E = 0.43429448190325182765112891892;      /* log_10 (e) */
        public static double M_SQRT2 = 1.41421356237309504880168872421;      /* sqrt(2) */
        public static double M_SQRT1_2 = 0.70710678118654752440084436210;      /* sqrt(1/2) */
        public static double M_SQRT3 = 1.73205080756887729352744634151;      /* sqrt(3) */
        public static double M_PI = 3.14159265358979323846264338328;      /* pi */
        public static double M_PI_2 = 1.57079632679489661923132169164;      /* pi/2 */
        public static double M_PI_4 = 0.78539816339744830961566084582;     /* pi/4 */
        public static double M_SQRTPI = 1.77245385090551602729816748334;      /* sqrt(pi) */
        public static double M_2_SQRTPI = 1.12837916709551257389615890312;      /* 2/sqrt(pi) */
        public static double M_1_PI = 0.31830988618379067153776752675;      /* 1/pi */
        public static double M_2_PI = 0.63661977236758134307553505349;      /* 2/pi */
        public static double M_LN10 = 2.30258509299404568401799145468;      /* ln(10) */
        public static double M_LN2 = 0.69314718055994530941723212146;      /* ln(2) */
        public static double M_LNPI = 1.14472988584940017414342735135;      /* ln(pi) */
        public static double M_EULER = 0.57721566490153286060651209008;      /* Euler constant */ 
        #endregion

        public static bool PZDebug = false;        
		public PZMath() { }

        #region static method

        public static bool IsOdd(int n)
        {
            return System.Math.IEEERemainder(n, 1) != 0.0;
        } // IsOdd()

        public static bool IsEven(int n)
        {
            return !IsOdd(n);
        } // IsEven()
        #endregion


    } // PZMath

	/// <summary>
	/// delegate of basic function
	/// </summary>
	public delegate double PZMath_delegate_f (double x, object parms);
	/// <summary>
	/// delegate of basic derivative function
	/// </summary>
	public delegate double PZMath_delegate_df (double x, object parms);
	/// <summary>
	/// delegeate of basic f/df
	/// </summary>
	public delegate double PZMath_delegate_fdf (double x, object parms, PZMath_delegate_f f, PZMath_delegate_df df);
	
	/// <summary>
	/// basic f class
	/// </summary>
	public class PZMath_f
	{
		public PZMath_f() {}
		public PZMath_delegate_f f;
		private object parms;

		public object Parms
		{
			set {parms = value;}
			get {return parms;}
		}
		/// <summary>
		/// evalue function
		/// </summary>
		/// <param name="x"></param>
		/// <returns></returns>
		public double FN_EVAL(double x)
		{
			return f(x, parms);
		} // PZMath_FN_EVAL();

	} // PZMath_f

	/// <summary>
	/// basic fdf class
	/// </summary>
	public class PZMath_fdf
	{
		public PZMath_fdf() {}

		public PZMath_delegate_f f;
		public PZMath_delegate_df df;
		public PZMath_delegate_fdf fdf;
		private object parms;

		public object Parms
		{
			set {parms = value;}
			get {return parms;}
		}

	} // PZMath_fdf
}
