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
		public static bool PZDebug = false;

		public PZMath() {}

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
