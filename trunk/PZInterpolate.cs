// PZ Interpolate functions
// 13.12.2005

using System;

namespace eee.Sheffield.PZ.Math
{
	
	/// <summary>
	/// delegate of interpolate function
	/// </summary>
	public delegate double PZInterpolateFunction (double r);	// interpolation function
	
	/// <summary>
	/// various interpolate functions
	/// source code from Peter I Rockett
	/// </summary>
	public class PZInterpolate
	{
		public PZInterpolate() {}

		#region interpolate functions

		/// <summary>
		/// Linear interpolate function
		/// Approximation to sinc sampling function - see Keys or Meijering et al
		/// </summary>
		/// <param name="r">distance from the center</param>
		/// <returns></returns>
		public static double LinearInterpolationFn(double r)
		{
			if( r < 0)
				r = -r;

			if(r >= 1.0)
				return 0.0;

			double a10 = -1, a00 = 1;

			return a10*r + a00;
		} // LinearInterpolationFn()

		/// <summary>
		/// Cubic interpolate function
		/// Approximation to sinc sampling function - see Keys or Meijering et al
		/// </summary>
		/// <param name="r">distance from the center</param>
		/// <returns></returns>
		public static double CubicInterpolationFn(double r)
		{
			if( r < 0)
				r = -r;

			if(r >= 2.0)
				return 0.0;

			double A = -0.5;	// alpha
			double a30 = A + 2, a20 = -(A + 3), a00 = 1;
			double a31 = A, a21 = -5*A, a11 = 8*A, a01 = -4*A;

			if(r < 1)
				return r*(r*(a30*r + a20) /*+ a10*/) + a00;
			else
				return r*(r*(a31*r + a21) + a11) + a01;
		} // CubicInterpolationFn()

		/// <summary>
		/// Quintic interpolate function
		/// Approximation to sinc sampling function - see Meijering et al
		/// </summary>
		/// <param name="r">distance from the center</param>
		/// <returns></returns>
		public static double QuinticInterpolationFn(double r)
		{
			if( r < 0)
				r = -r;

			if(r >= 3.0)
				return 0.0;

			double A = 3.0 / 64.0;	// alpha
			double a50 = 10*A - 21.0/16.0, a40 = -18*A + 45.0/16.0, a20 = 8*A -2.5, a00 = 1;
			double a51 = 11*A - 5.0/16.0, a41 = -88*A +45.0/16.0, a31 = 270*A-10, a21 = -392*A + 17.5, a11 = 265*A - 15, a01 = -66*A + 5;
			double a52 = A, a42 = -14*A, a32 = 78*A, a22 = -216*A, a12 = 297*A, a02 = -162*A;

			if(r < 2)
				if(r < 1)
					return r*(r*(r*(r*(a50*r + a40)) + a20)) + a00;

				else	// 1 <= r < 2
					return r*(r*(r*(r*(a51*r + a41) + a31) + a21) + a11) + a01;

			else	// 2 <= r < 3 
				return r*(r*(r*(r*(a52*r + a42) + a32) + a22) + a12) + a02;

		} // QuinticInterpolationFn()
		
		/// <summary>
		/// Septic interpolate function
		/// Approximation to sinc sampling function - see Meijering et al
		/// </summary>
		/// <param name="r">distance from the center</param>
		/// <returns></returns>
		public static double SepticInterpolationFn(double r)
		{
			if( r < 0)
				r = -r;

			if(r >= 4.0)
				return 0.0;

			double A = -71.0 / 83232.0;	// alpha
			double a70 = 245*A + 821.0/1734.0,a60 = -621*A - 1148.0/867.0,a40 = 760*A + 1960.0/867.0,a30 = 0,a20 = -384*A - 1393.0/578.0,a10 = 0,a00 = 1;
			double a71 = 301*A + 1687.0/6936.0,a61 = -3309*A - 2492.0/867.0,a51 = 14952*A + 32683.0/2312.0,a41 = -35640*A - 128695.0/3468.0,a31 = 47880*A + 127575.0/2312.0,a21 = -36000*A - 13006.0/289.0,a11 = 14168*A + 120407.0/6936.0,a01 = -2352*A - 2233.0/1156.0;
			double a72 = 57*A + 35.0/6936.0,a62 = -1083*A - 175.0/1734.0,a52 = 8736*A + 1995.0/2312.0,a42 = -38720*A - 4725.0/1156.0,a32 = 101640*A + 1575.0/136.0,a22 = -157632*A - 5670.0/289.0,a12 = 133336*A + 42525.0/2312.0,a02 = -47280*A - 8505.0/1156.0;
			double a73 = A,a63 = -27.0*A,a53 = 312.0*A,a43 = -2000.0*A,a33 = 7680.0*A,a23 = -17664.0*A,a13 = 22528.0*A,a03 = -12288.0*A;

			if(r < 3)
				{
				if(r < 2)
					{
					if(r < 1)
						return r*(r*(r*(r*(r*(r*(a70*r + a60) /*+ a50*/) + a40) + a30) + a20) + a10) + a00;
					else	// 1 <= r < 2
						return r*(r*(r*(r*(r*(r*(a71*r + a61) + a51) + a41) + a31) + a21) + a11) + a01;
					}
				else	 // 2 <= r < 3
					return r*(r*(r*(r*(r*(r*(a72*r + a62) + a52) + a42) + a32) + a22) + a12) + a02;
				}
			else	// 3 <= r < 4
				return r*(r*(r*(r*(r*(r*(a73*r + a63) + a53) + a43) + a33) + a23) + a13) + a03;
		} // SepticInterpolationFn()

		/// <summary>
		/// Blakman Harris 
		/// Approximation to sinc sampling function by a Blackman-Harris window
		///- see Lehmann et al
		/// the interpolation kernel = 6 by 6;
		/// </summary>
		/// <param name="r">distance from the center</param>
		/// <returns></returns>
		public static double BlackmanHarrisInterpolationFn(double r)
		{
			if( r < 0)
				r = -r;

			if(r >= 3.0)
				return 0.0;

			double omiga0 = 0.42323;
			double omiga1 = 0.49755;
			double omiga2 = 0.07922;
			double pi_local = 3.14159;

            return System.Math.Sin(pi_local * r) / pi_local / r * (omiga0 + omiga1 * System.Math.Cos(2 * pi_local * 2 * r / 6) + omiga2 * System.Math.Cos(2 * pi_local * 4 * r / 6));
		} // BlackmanHarrisInterpolationFn()
		#endregion
	}
}
