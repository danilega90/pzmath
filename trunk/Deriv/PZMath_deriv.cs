// PZ Math deriv
// Based on Math.NET Open source
//
// -- numerical difference
//
// -- Class PZMath_deriv
// GNU Copyright
/* deriv/deriv.c
 * 
 * Copyright (C) 2004 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */
// 23.11.2005


using System;

namespace eee.Sheffield.PZ.Math
{
	/// <summary>
	/// basic f class
	/// </summary>
	public class PZMath_deriv
	{
		public PZMath_deriv()
		{
		}

		public int Central (PZMath_f f, double x, double h, out double result, out double abserr)
		{
			double r_0, round, trunc, error;
			CentralDeriv (f, x, h, out r_0, out round, out trunc);
			error = round + trunc;
			if (round < trunc && (round > 0 && trunc > 0))
			{
				double r_opt, round_opt, trunc_opt, error_opt;

				/* Compute an optimised stepsize to minimize the total error,
				using the scaling of the truncation error (O(h^2)) and
				rounding error (O(1/h)). */

				double h_opt = h * System.Math.Pow (round / (2.0 * trunc), 1.0 / 3.0);
				CentralDeriv (f, x, h_opt, out r_opt, out round_opt, out trunc_opt);
				error_opt = round_opt + trunc_opt;

				/* Check that the new error is smaller, and that the new derivative 
				is consistent with the error bounds of the original estimate. */

                if (error_opt < error && System.Math.Abs(r_opt - r_0) < 4.0 * error)
				{
					r_0 = r_opt;
					error = error_opt;
				}
			}

			result = r_0;
			abserr = error;

			return PZMath_errno.PZMath_SUCCESS; // 1 means sucess;
		}  //PZMath_deriv_central()

		private static void CentralDeriv 
			(PZMath_f f, double x, double h, out double result, out double abserr_round, out double abserr_trunc)
		{
			/* Compute the derivative using the 5-point rule (x-h, x-h/2, x,
			x+h/2, x+h). Note that the central point is not used.  

			Compute the error using the difference between the 5-point and
			the 3-point rule (x-h,x,x+h). Again the central point is not
			used. */
			double fm1 = f.FN_EVAL(x - h);
			double fp1 = f.FN_EVAL(x + h);

			double fmh = f.FN_EVAL (x - h / 2);
			double fph = f.FN_EVAL (x + h / 2);

			double r3 = 0.5 * (fp1 - fm1);
			double r5 = (4.0 / 3.0) * (fph - fmh) - (1.0 / 3.0) * r3;

            double e3 = (System.Math.Abs(fp1) + System.Math.Abs(fm1)) * PZMath_machine.PZMath_DBL_EPSILON;
            double e5 = 2.0 * (System.Math.Abs(fph) + System.Math.Abs(fmh)) * PZMath_machine.PZMath_DBL_EPSILON + e3;

            double dy = System.Math.Max(System.Math.Abs(r3), System.Math.Abs(r5)) * System.Math.Abs(x) * PZMath_machine.PZMath_DBL_EPSILON;

			/* The truncation error in the r5 approximation itself is O(h^4).
			However, for safety, we estimate the error from r5-r3, which is
			O(h^2).  By scaling h we will minimise this estimated error, not
			the actual truncation error in r5. */

			result = r5 / h;
            abserr_trunc = System.Math.Abs((r5 - r3) / h); /* Estimated truncation error O(h^2) */
            abserr_round = System.Math.Abs(e5 / h) + dy;   /* Rounding error (cancellations) */
		} // central_deriv()

	} // PZMath_deriv
}
