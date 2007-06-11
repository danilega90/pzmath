// PZ Math nmsimplex, Nelder Mead Simplex

// 19.12.2005


using System;

namespace eee.Sheffield.PZ.Math
{
	/// <summary>
	/// Nelder Mead Simplex method state class
	/// </summary>
	public class nmsimplex_state_t
	{
		public PZMath_matrix x1 = null;	/* simplex corner points */
		public PZMath_vector y1 = null;/* function value at corner points */
		public PZMath_vector ws1 = null; /* workspace 1 for algorithm */
		public PZMath_vector ws2 = null;/* workspace 2 for algorithm */
	} // nmsimplex_state_t

	public class nmsimplex
	{
		#region comments
		/* multimin/simplex.c
 * 
 * Copyright (C) 2002 Tuomo Keskitalo, Ivo Alxneit
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

		/*
		   - Originally written by Tuomo Keskitalo <tuomo.keskitalo@iki.fi>
		   - Corrections to nmsimplex_iterate and other functions 
			 by Ivo Alxneit <ivo.alxneit@psi.ch>
		   - Additional help by Brian Gough <bjg@network-theory.co.uk>
		*/

		/* The Simplex method of Nelder and Mead,
		   also known as the polytope search alogorithm. Ref:
		   Nelder, J.A., Mead, R., Computer Journal 7 (1965) pp. 308-313.

		   This implementation uses n+1 corner points in the simplex.
		*/
		#endregion
		/// <summary>
		/// empty constructor
		/// </summary>
		public nmsimplex() {}

		public static int Alloc (object vstate, int n)
		{
			(vstate as nmsimplex_state_t).x1 = new PZMath_matrix(n + 1, n);
			(vstate as nmsimplex_state_t).y1 = new PZMath_vector(n + 1);
			(vstate as nmsimplex_state_t).ws1 = new PZMath_vector(n);
			(vstate as nmsimplex_state_t).ws2 = new PZMath_vector(n);
			return PZMath_errno.PZMath_SUCCESS;
		} // alloc()

		public static int Init (object vstate, PZMath_multimin_function f, PZMath_vector x, 
			out double size, PZMath_vector step_size)
		{
			int i;
			double val;
			PZMath_vector xtemp = (vstate as nmsimplex_state_t).ws1;

			/* first point is the original x0 */
			val = f.FN_EVAL(x);
			
			(vstate as nmsimplex_state_t).x1.SetRow(0, x);
			(vstate as nmsimplex_state_t).y1[0] = val;

			/* following points are initialized to x0 + step_size */

			for (i = 0; i < x.size; i++)
			{
				xtemp.MemCopyFrom(x);
				val = xtemp[i] + step_size[i];
				xtemp[i] = val;
				val = f.FN_EVAL(xtemp);
				(vstate as nmsimplex_state_t).x1.SetRow(i + 1, xtemp);
				(vstate as nmsimplex_state_t).y1[i + 1] = val;
			}

			/* Initialize simplex size */
			size = Size (vstate as nmsimplex_state_t);

			return PZMath_errno.PZMath_SUCCESS;
		} // Init()

		public static int Iterate (object vstate, PZMath_multimin_function f,
			PZMath_vector x, out double size, out double fval)
		{

			/* Simplex iteration tries to minimize function f value */
			/* Includes corrections from Ivo Alxneit <ivo.alxneit@psi.ch> */

			/* xc and xc2 vectors store tried corner point coordinates */

			PZMath_vector xc = (vstate as nmsimplex_state_t).ws1;
			PZMath_vector xc2 = (vstate as nmsimplex_state_t).ws2;
			PZMath_vector y1 = (vstate as nmsimplex_state_t).y1;
			PZMath_matrix x1 = (vstate as nmsimplex_state_t).x1;

			int n = y1.size;
			int i;
			int hi = 0, s_hi = 0, lo = 0;
			double dhi, ds_hi, dlo;
			int status;
			double val, val2;

			/* get index of highest, second highest and lowest point */

			dhi = ds_hi = dlo = y1[0];

			for (i = 1; i < n; i++)
			{
				val = y1[i];
				if (val < dlo)
				{
					dlo = val;
					lo = i;
				}
				else if (val > dhi)
				{
					ds_hi = dhi;
					s_hi = hi;
					dhi = val;
					hi = i;
				}
				else if (val > ds_hi)
				{
					ds_hi = val;
					s_hi = i;
				}
			}

			/* reflect the highest value */
			val = MoveCorner (-1.0, vstate as nmsimplex_state_t, hi, xc, f);

			if (val < y1[lo])
			{

				/* reflected point becomes lowest point, try expansion */
				val2 = MoveCorner (-2.0, vstate as nmsimplex_state_t, hi, xc2, f);

				if (val2 < y1[lo])
				{
					x1.SetRow(hi, xc2);
					y1[hi] = val2;
				}
				else
				{
					x1.SetRow(hi, xc);
					y1[hi] = val;
				}
			}

				/* reflection does not improve things enough */

			else if (val > y1[s_hi])
			{
				if (val <= y1[hi])
				{

					/* if trial point is better than highest point, replace 
					   highest point */
					x1.SetRow(hi, xc);
					y1[hi] = val;
				}

				/* try one dimensional contraction */
				val2 = MoveCorner (0.5, vstate as nmsimplex_state_t, hi, xc2, f);

				if (val2 <= y1[hi])
				{
					x1.SetRow(hi, xc2);
					y1[hi] = val2;
				}

				else
				{

					/* contract the whole simplex in respect to the best point */
					status = ContractByBest (vstate as nmsimplex_state_t, lo, xc, f);
					if (status != 0)
					{
						PZMath_errno.ERROR("nmsimplex::Iterate(), nmsimplex_contract_by_best failed", PZMath_errno.PZMath_FAILURE);
					}
				}
			}
			else
			{

				/* trial point is better than second highest point. 
				   Replace highest point by it */

				x1.SetRow(hi, xc);
				y1[hi] = val;
			}

			/* return lowest point of simplex as x */

			lo = y1.MinIndex();
			x1.CopyRow2(lo, x);
			fval = y1[lo];

			/* Update simplex size */

			size = Size (vstate as nmsimplex_state_t);

			return PZMath_errno.PZMath_SUCCESS;
		} // Iterate()
		
		#region nmsimplex static methods
		public static double Size (nmsimplex_state_t state)
		{
			/* calculates simplex size as average sum of length of vectors 
			   from simplex center to corner points:     

			   ( sum ( || y - y_middlepoint || ) ) / n 
			 */

			PZMath_vector s = state.ws1;
			PZMath_vector mp = state.ws2;
			PZMath_matrix x1 = state.x1;
			int i;
			double ss = 0.0;

			/* Calculate middle point */
			CalcCenter (state, mp);

			for (i = 0; i < x1.RowCount; i++)
			{
				
				x1.CopyRow2(i, s);
				PZMath_blas.Daxpy(-1.0, mp, s);
				ss += PZMath_blas.Dnrm2(s);
			}

			return ss / (double) (x1.RowCount);
		} // Size()

		public static int CalcCenter (nmsimplex_state_t state, PZMath_vector mp)
		{
			/* calculates the center of the simplex to mp */

			PZMath_matrix x1 = state.x1;

			int i, j;
			double val;

			for (j = 0; j < x1.ColumnCount; j++)
			{
				val = 0.0;
				for (i = 0; i < x1.RowCount; i++)
				{
					val += x1[i, j];
				}
				val /= x1.RowCount;
				mp[j] = val;
			}

			return PZMath_errno.PZMath_SUCCESS;
		} // CalcCenter()

		public static double MoveCorner (double coeff, nmsimplex_state_t state, int corner, 
			PZMath_vector xc, PZMath_multimin_function f)
		{
			/* moves a simplex corner scaled by coeff (negative value represents 
				mirroring by the middle point of the "other" corner points)
				and gives new corner in xc and function value at xc as a 
				return value 
			*/

			PZMath_matrix x1 = state.x1;

			int i, j;
			double newval, mp;

			if (x1.RowCount < 2)
			{
				PZMath_errno.ERROR ("nmsimplex::MoveCorner(), simplex cannot have less than two corners!", PZMath_errno.PZMath_FAILURE);
			}

			for (j = 0; j < x1.ColumnCount; j++)
			{
				mp = 0.0;
				for (i = 0; i < x1.RowCount; i++)
				{
					if (i != corner)
					{
						mp += (x1[i, j]);
					}
				}
				mp /= (double) (x1.RowCount - 1);
				newval = mp - coeff * (mp - x1[corner, j]);
				xc[j] = newval;
			}

			newval = f.FN_EVAL(xc);

			return newval;
		} // MoveCorner()

		public static int ContractByBest (nmsimplex_state_t state, int best, PZMath_vector xc, 
			PZMath_multimin_function f)
		{

			/* Function contracts the simplex in respect to 
			   best valued corner. That is, all corners besides the 
			   best corner are moved. */

			/* the xc vector is simply work space here */

			PZMath_matrix x1 = state.x1;
			PZMath_vector y1 = state.y1;

			int  i, j;
			double newval;

			for (i = 0; i < x1.RowCount; i++)
			{
				if (i != best)
				{
					for (j = 0; j < x1.ColumnCount; j++)
					{
						newval = 0.5 * (x1[i, j] + x1[best, j]);
						x1[i, j] = newval;
					}

					/* evaluate function in the new point */

					x1.CopyRow2(i, xc);
					newval = f.FN_EVAL(xc);
					y1[i] = newval;
				}
			}

			return PZMath_errno.PZMath_SUCCESS;
		} // ContractByBest()
		#endregion
	}
}
