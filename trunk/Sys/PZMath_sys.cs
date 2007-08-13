// PZ Math
//
// -- sys methods
//
// -- Class PZMath_sys
// 25.11.2005


using System;

namespace eee.Sheffield.PZ.Math
{
	/// <summary>
	/// sys methods
	/// </summary>
	
	public class PZMath_sys
	{
		public static double Hypot (double x, double y)
		{
            double xabs = System.Math.Abs(x);
            double yabs = System.Math.Abs(y);
			double min, max;
	
			if (xabs < yabs) 
			{
				min = xabs ;
				max = yabs ;
			} 
			else 
			{
				min = yabs ;
				max = xabs ;
			}
	
			if (min == 0) 
			{
				return max ;
			}

			{
				double u = min / max ;
                return max * System.Math.Sqrt(1 + u * u);
			}
		} // Hypot()

		public static bool Finite(double x)
		{
			double y = x - x;
			//bool status = (y == y);
            bool status = true;
			return status;
		} // Finite()

		public static bool Isnan (double x)
		{
			//bool status = (x != x);
            bool status = false;
			return status;
		} // Isnan()

	} // PZMath_sys
}
