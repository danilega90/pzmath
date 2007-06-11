// PZ Math Multifit nonlinear
//
// -- numerical multi non linear regression
//
// -- Class PZMath_multifit_function
// -- Class PZMath_multifit_function_fdf
//
// GNU Copyright
/* multifit_nlin/gsl_multifit_nlin.h
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Brian Gough
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
// 24.11.2005


using System;

namespace eee.Sheffield.PZ.Math
{
	public delegate int multifit_delegate_f (PZMath_vector x, object parms, PZMath_vector f); 
	public delegate int multifit_delegate_df (PZMath_vector x, object parms, PZMath_matrix df); 
	public delegate int multifit_delegate_fdf (PZMath_vector x, object parms, PZMath_vector f, PZMath_matrix df);

	public class PZMath_multifit_function
	{
		public multifit_delegate_f f;
		public int n;		/* number of functions */
		public int p;		/* number of independent variables */
		private object parms;
		public object Parms
		{
			get {return parms;}
			set {parms = value;}
		}
		public PZMath_multifit_function() {}

		public int FN_EVAL (PZMath_vector x, PZMath_vector y)
		{
			return f(x, parms, y);
		} // FN_EVAL()
	} // PZMath_multifit_function

	public class PZMath_multifit_function_fdf
	{
		public multifit_delegate_f f;
		public multifit_delegate_df df;
		public multifit_delegate_fdf fdf;
		public int n;	/* number of functions */
		public int p;	/* number of independent variables */
		public object parms;

		public PZMath_multifit_function_fdf() {}

		public int FN_EVAL_F (PZMath_vector x, PZMath_vector y)
		{
			return f(x, parms, y);
		} // FN_EVAL_F()

		public int FN_EVAL_DF (PZMath_vector x, PZMath_matrix dy)
		{
			return df(x, parms, dy);
		} // FN_EVAL_DF()

		public int FN_EVAL_F_DF (PZMath_vector x, PZMath_vector y, PZMath_matrix dy)
		{
			return fdf(x, parms, y, dy);
		} // FN_EVAL_F_DF()
	} //PZMath_multifit_function_fdf
}
