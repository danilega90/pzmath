// PZ Math integration

/* gsl 1.07/integration/gsl_integration.h
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
// 11.02.2006

using System;

namespace eee.Sheffield.PZ.Math
{	
	public delegate void PZMath_integration_rules (PZMath_f f, double a, double b, 
		ref double result, ref double abserr, ref double defabs, ref double resabs);

	/// <summary>
	/// integration
	/// </summary>
	public class PZMath_integration
	{
		public PZMath_integration() {}

		public static int Qags(PZMath_f f, double a, double b,
			double epsabs, double epsrel, int limit,PZMath_integration_workspace workspace,
			ref double result, ref double abserr)
		{
			PZMath_integration_rules q = new PZMath_integration_rules(PZMath_integration.QK21);
			double area, errsum;
			double res_ext, err_ext;
			double result0 = 0.0, abserr0 = 0.0, resabs0 = 0.0, resasc0 = 0.0;
			double tolerance;

			double ertest = 0;
			double error_over_large_intervals = 0;
			double reseps = 0, abseps = 0, correc = 0;
			int ktmin = 0;
			int roundoff_type1 = 0, roundoff_type2 = 0, roundoff_type3 = 0;
			int error_type = 0, error_type2 = 0;

			int iteration = 0;

			int positive_integrand = 0;
			int extrapolate = 0;
			int disallow_extrapolation = 0;

			ExtrapolationTable table = new ExtrapolationTable();

			/* Initialize results */
			workspace.Initialise(a, b);
			
			result = 0;
			abserr = 0;

			if (limit > workspace.limit)
			{
                //PZMath_errno.ERROR ("iteration limit exceeds available workspace", PZMath_errno.PZMath_EINVAL) ;
                //return 0;
                return PZMath_errno.PZMath_EINVAL;
			}

			/* Test on accuracy */

			if (epsabs <= 0 && (epsrel < 50 * PZMath_machine.PZMath_DBL_EPSILON || epsrel < 0.5e-28))	
			{
                //PZMath_errno.ERROR ("tolerance cannot be acheived with given epsabs and epsrel",
                //    PZMath_errno.PZMath_EBADTOL);
                //return 0;
                return PZMath_errno.PZMath_EBADTOL;
			}
			/* Perform the first integration */

			q (f, a, b, ref result0, ref abserr0, ref resabs0, ref resasc0);

			workspace.SetInitialResult (result0, abserr0);

			tolerance = System.Math.Max (epsabs, epsrel * System.Math.Abs (result0));

			if (abserr0 <= 100 * PZMath_machine.PZMath_DBL_EPSILON * resabs0 && abserr0 > tolerance)
			{
				result = result0;
				abserr = abserr0;

                //PZMath_errno.ERROR ("cannot reach tolerance because of roundoff error on first attempt", PZMath_errno.PZMath_EROUND);
                //return 0;
                return PZMath_errno.PZMath_EROUND;
			}
			else if ((abserr0 <= tolerance && abserr0 != resasc0) || abserr0 == 0.0)
			{
				result = result0;
				abserr = abserr0;

				return PZMath_errno.PZMath_SUCCESS;
			}
			else if (limit == 1)
			{
				result = result0;
				abserr = abserr0;

                //PZMath_errno.ERROR ("a maximum of one iteration was insufficient", PZMath_errno.PZMath_EMAXITER);
                //return 0;
                return PZMath_errno.PZMath_EMAXITER;
			}

			/* Initialization */
			
			table.InitialiseTable();
			table.AppendTable(result0);

			area = result0;
			errsum = abserr0;

			res_ext = result0;
			err_ext = PZMath_machine.PZMath_DBL_MAX;

			positive_integrand = TestPositivity (result0, resabs0);
			iteration = 1;

			do 
			{
				int current_level;
				double a1, b1, a2, b2;
				double a_i = 0.0, b_i = 0.0, r_i = 0.0, e_i = 0.0;
				double area1 = 0, area2 = 0, area12 = 0;
				double error1 = 0, error2 = 0, error12 = 0;
				double resasc1 = 0.0, resasc2 = 0.0;
				double resabs1 = 0.0, resabs2 = 0.0;
				double last_e_i;

				/* Bisect the subinterval with the largest error estimate */

				workspace.Retrieve (ref a_i, ref b_i, ref r_i, ref e_i);

				current_level = workspace.level[workspace.i] + 1;

				a1 = a_i;
				b1 = 0.5 * (a_i + b_i);
				a2 = b1;
				b2 = b_i;

				iteration ++;

				q (f, a1, b1, ref area1, ref error1, ref resabs1, ref resasc1);
				q (f, a2, b2, ref area2, ref error2, ref resabs2, ref resasc2);

				area12 = area1 + area2;
				error12 = error1 + error2;
				last_e_i = e_i;

				/* Improve previous approximations to the integral and test for
					accuracy.

					We write these expressions in the same way as the original
					QUADPACK code so that the rounding errors are the same, which
					makes testing easier. */

				errsum = errsum + error12 - e_i;
				area = area + area12 - r_i;

				tolerance = System.Math.Max (epsabs, epsrel * System.Math.Abs (area));

				if (resasc1 != error1 && resasc2 != error2)
				{
					double delta = r_i - area12;

					if (System.Math.Abs (delta) <= 1.0e-5 * System.Math.Abs (area12) && error12 >= 0.99 * e_i)
					{
						if (extrapolate == 0)
							roundoff_type1++;
						else
							roundoff_type2++;					
					}
					if (iteration > 10 && error12 > e_i)
						roundoff_type3++;
				}

				/* Test for roundoff and eventually set error flag */

				if (roundoff_type1 + roundoff_type2 >= 10 || roundoff_type3 >= 20)
					error_type = 2;       /* round off error */

				if (roundoff_type2 >= 5)
					error_type2 = 1;

				/* set error flag in the case of bad integrand behaviour at
					a point of the integration range */

				if (SubintervalTooSmall (a1, a2, b2))
					error_type = 4;

				/* append the newly-created intervals to the list */

				workspace.Update (a1, b1, area1, error1, a2, b2, area2, error2);

				if (errsum <= tolerance)
					goto compute_result;

				if (error_type >= 1)
					break;

				if (iteration >= limit - 1)
				{
					error_type = 1;
					break;
				}

				if (iteration == 2)       /* set up variables on first iteration */
				{
					error_over_large_intervals = errsum;
					ertest = tolerance;
					table.AppendTable (area);
					continue;
				}

				if (disallow_extrapolation >= 1)
					continue;

				error_over_large_intervals += -last_e_i;

				if (current_level < workspace.maximum_level)
					error_over_large_intervals += error12;

				if (extrapolate == 0)
				{
					/* test whether the interval to be bisected next is the
						smallest interval. */

					if (workspace.LargeInterval ())
						continue;

					extrapolate = 1;
					workspace.nrmax = 1;
				}

				if (error_type2 == 0 && error_over_large_intervals > ertest)
				{
					if (workspace.IncreaseNrmax())
						continue;
				}

				/* Perform extrapolation */

				table.AppendTable (area);

				table.Qelg(ref reseps, ref abseps);

				ktmin ++;

				if (ktmin > 5 && err_ext < 0.001 * errsum)
					error_type = 5;

				if (abseps < err_ext)
				{
					ktmin = 0;
					err_ext = abseps;
					res_ext = reseps;
					correc = error_over_large_intervals;
					ertest = System.Math.Max (epsabs, epsrel * System.Math.Abs (reseps));
					if (err_ext <= ertest)
						break;
				}

				/* Prepare bisection of the smallest interval. */

				if (table.n == 1)
					disallow_extrapolation = 1;

				if (error_type == 5)
					break;

				/* work on interval with largest error */

				workspace.ResetNrmax();
				extrapolate = 0;
				error_over_large_intervals = errsum;

			}
			while (iteration < limit);

			result = res_ext;
			abserr = err_ext;

			if (err_ext == PZMath_machine.PZMath_DBL_MAX)
				goto compute_result;

			if (error_type >=1 || error_type2 >= 1)
			{
				if (error_type2 >= 1)
					err_ext += correc;

				if (error_type == 0)
					error_type = 3;

				if (res_ext != 0.0 && area != 0.0)
				{
					if (err_ext / System.Math.Abs (res_ext) > errsum / System.Math.Abs (area))
						goto compute_result;
				}
				else if (err_ext > errsum)
					goto compute_result;
				else if (area == 0.0)
					goto return_error;
			}

			/*  Test on divergence. */

			{
				double max_area = System.Math.Max (System.Math.Abs (res_ext), System.Math.Abs (area));

				if (positive_integrand  == 0 && max_area < 0.01 * resabs0)
					goto return_error;
			}

			{
				double ratio = res_ext / area;

				if (ratio < 0.01 || ratio > 100.0 || errsum > System.Math.Abs (area))
				error_type = 6;
			}

			goto return_error;

compute_result:

			result = workspace.SumResults();
			abserr = errsum;

return_error:

			if (error_type > 2)
				error_type--;

			if (error_type == 0) 
				return PZMath_errno.PZMath_SUCCESS;
			else if (error_type == 1)
			{
				//PZMath_errno.ERROR("number of iterations was insufficient", PZMath_errno.PZMath_EMAXITER);
				//return 0;
                return PZMath_errno.PZMath_EMAXITER;
			}
			else if (error_type == 2)
			{
				//PZMath_errno.ERROR ("cannot reach tolerance because of roundoff error",
				//	PZMath_errno.PZMath_EROUND);
				//return 0;
                return PZMath_errno.PZMath_EROUND;
			}
			else if (error_type == 3)
			{
                //PZMath_errno.ERROR ("bad integrand behavior found in the integration interval",
                //    PZMath_errno.PZMath_ESING);
                //return 0;
                return PZMath_errno.PZMath_ESING;
			}
			else if (error_type == 4)
			{
                //PZMath_errno.ERROR ("roundoff error detected in the extrapolation table",
                //    PZMath_errno.PZMath_EROUND);
                //return 0;
                return PZMath_errno.PZMath_EROUND;
			}
			else if (error_type == 5)
			{
                //PZMath_errno.ERROR ("integral is divergent, or slowly convergent",
                //    PZMath_errno.PZMath_EDIVERGE);
                //return 0;
                return PZMath_errno.PZMath_EDIVERGE;
			}
			else
			{
                //PZMath_errno.ERROR ("could not integrate function", PZMath_errno.PZMath_EFAILED);
                //return 0;
                return PZMath_errno.PZMath_EFAILED;
			}
		} // Qags()
		
		#region integration rules (QK##)
		#region QK21 comments
		/* Gauss quadrature weights and kronrod quadrature abscissae and
   weights as evaluated with 80 decimal digit arithmetic by
   L. W. Fullerton, Bell Labs, Nov. 1981. */
		#endregion
		public static void QK21(PZMath_f f, double a, double b, 
			ref double result, ref double abserr, ref double resabs, ref double resasc)
		{
			double [] xgk =   /* abscissae of the 21-point kronrod rule */
			{
				0.995657163025808080735527280689003,
				0.973906528517171720077964012084452,
				0.930157491355708226001207180059508,
				0.865063366688984510732096688423493,
				0.780817726586416897063717578345042,
				0.679409568299024406234327365114874,
				0.562757134668604683339000099272694,
				0.433395394129247190799265943165784,
				0.294392862701460198131126603103866,
				0.148874338981631210884826001129720,
				0.000000000000000000000000000000000
			};

			/* xgk[1], xgk[3], ... abscissae of the 10-point gauss rule. 
			   xgk[0], xgk[2], ... abscissae to optimally extend the 10-point gauss rule */

			double [] wg =     /* weights of the 10-point gauss rule */
			{
				0.066671344308688137593568809893332,
				0.149451349150580593145776339657697,
				0.219086362515982043995534934228163,
				0.269266719309996355091226921569469,
				0.295524224714752870173892994651338
			};

		    double [] wgk =   /* weights of the 21-point kronrod rule */
			{
				0.011694638867371874278064396062192,
				0.032558162307964727478818972459390,
				0.054755896574351996031381300244580,
				0.075039674810919952767043140916190,
				0.093125454583697605535065465083366,
				0.109387158802297641899210590325805,
				0.123491976262065851077958109831074,
				0.134709217311473325928054001771707,
				0.142775938577060080797094273138717,
				0.147739104901338491374841515972068,
				0.149445554002916905664936468389821
			};

			double [] fv1 = new double [11];
			double [] fv2 = new double [11];
			QK (11, xgk, wg, wgk, fv1, fv2, f, a, b, ref result, ref abserr, ref resabs, ref resasc);
		} // QK21()

		#region QK comments
		#endregion
		static void QK (int n, double [] xgk, double [] wg, double [] wgk, 
			double [] fv1, double [] fv2,
			PZMath_f f, double a, double b,
			ref double result, ref double abserr, ref double resabs, ref double resasc)
		{

			double center = 0.5 * (a + b);
			double half_length = 0.5 * (b - a);
			double abs_half_length = System.Math.Abs (half_length);
			double f_center = f.FN_EVAL(center);

			double result_gauss = 0;
			double result_kronrod = f_center * wgk[n - 1];

			double result_abs = System.Math.Abs (result_kronrod);
			double result_asc = 0;
			double mean = 0, err = 0;

			int j;

			if (n % 2 == 0)
				result_gauss = f_center * wg[n / 2 - 1];

			for (j = 0; j < (n - 1) / 2; j++)
			{
				int jtw = j * 2 + 1;        /* j=1,2,3 jtw=2,4,6 */
				double abscissa = half_length * xgk[jtw];
				double fval1 = f.FN_EVAL (center - abscissa);
				double fval2 = f.FN_EVAL (center + abscissa);
				double fsum = fval1 + fval2;
				fv1[jtw] = fval1;
				fv2[jtw] = fval2;
				result_gauss += wg[j] * fsum;
				result_kronrod += wgk[jtw] * fsum;
				result_abs += wgk[jtw] * (System.Math.Abs (fval1) + System.Math.Abs (fval2));
			}

			for (j = 0; j < n / 2; j++)
			{
				int jtwm1 = j * 2;
				double abscissa = half_length * xgk[jtwm1];
				double fval1 = f.FN_EVAL (center - abscissa);
				double fval2 = f.FN_EVAL (center + abscissa);
				fv1[jtwm1] = fval1;
				fv2[jtwm1] = fval2;
				result_kronrod += wgk[jtwm1] * (fval1 + fval2);
				result_abs += wgk[jtwm1] * (System.Math.Abs (fval1) + System.Math.Abs (fval2));
			}

			mean = result_kronrod * 0.5;

			result_asc = wgk[n - 1] * System.Math.Abs (f_center - mean);

			for (j = 0; j < n - 1; j++)
				result_asc += wgk[j] * (System.Math.Abs (fv1[j] - mean) + System.Math.Abs (fv2[j] - mean));

			/* scale by the width of the integration region */

			err = (result_kronrod - result_gauss) * half_length;

			result_kronrod *= half_length;
			result_abs *= abs_half_length;
			result_asc *= abs_half_length;

			result = result_kronrod;
			resabs = result_abs;
			resasc = result_asc;
			abserr = RescaleError (err, result_abs, result_asc);
		} // QK()

		static double RescaleError(double err, double result_abs, double result_asc)
		{
			err = System.Math.Abs (err) ;

			if (result_asc != 0 && err != 0)
			{
					double scale = System.Math.Pow ((200 * err / result_asc), 1.5) ;      
					if (scale < 1)
						err = result_asc * scale ;
					else 
						err = result_asc ;
			}
			if (result_abs > PZMath_machine.PZMath_DBL_MIN / (50 * PZMath_machine.PZMath_DBL_EPSILON))
			{
				double min_err = 50 * PZMath_machine.PZMath_DBL_EPSILON * result_abs ;
				if (min_err > err) 
					err = min_err ;
			}	  
			return err ;
		} // RescaleError()
		#endregion

		#region util methods
		static int TestPositivity (double result, double resabs)
		{
			if (System.Math.Abs (result) >= (1 - 50 * PZMath_machine.PZMath_DBL_EPSILON) * resabs)
				return 1;
			else
				return 0;
		} // TestPositivity()

		static bool SubintervalTooSmall (double a1, double a2, double b2)
		{
			const double e = PZMath_machine.PZMath_DBL_EPSILON;
			const double u = PZMath_machine.PZMath_DBL_MIN;

			double tmp = (1 + 100 * e) * (System.Math.Abs (a2) + 1000 * u);

			if (System.Math.Abs (a1) <= tmp && System.Math.Abs  (b2) <= tmp)
				return true;
			else 
				return false;
		} // SubintervalTooSmall()
		#endregion
	} // PZMath_integration

	/// <summary>
	/// integration workspace
	/// </summary>
	public class PZMath_integration_workspace
	{
		public int limit;
		public int size;
		public int nrmax;
		public int i;
		public int maximum_level;
		public double [] alist = null;
		public double [] blist = null;
		public double [] rlist = null;
		public double [] elist = null;
		public int [] order = null;
		public int [] level = null;

		public PZMath_integration_workspace() {}
		
		#region methods
		public void Alloc(int n)
		{
			if (n == 0)
			{
				PZMath_errno.ERROR ("workspace length n must be positive integer",
					PZMath_errno.PZMath_EDOM);
			}
			alist = new double[n];
			blist = new double[n];
			rlist = new double[n];
			elist = new double[n];
			order = new int[n];
			level = new int[n];
			size = 0;
			limit = n;
			maximum_level = 0;
		} // Alloc()

		public void Initialise(double a, double b)
		{
			size = 0;
			nrmax = 0;
			i = 0;
			alist[0] = a;
			blist[0] = b;
			rlist[0] = 0.0;
			elist[0] = 0.0;
			order[0] = 0;
			level[0] = 0;
			maximum_level = 0;
		} // Initialise()

		public void SetInitialResult (double result, double error)
		{
			size = 1;
			rlist[0] = result;
			elist[0] = error;
		} // SetInitialResult()

		public void Retrieve (ref double a, ref double b, ref double r, ref double e)
		{
			a = alist[i] ;
			b = blist[i] ;
			r = rlist[i] ;
			e = elist[i] ;
		} // Retrieve()

		public void Update (double a1, double b1, double area1, double error1,
             double a2, double b2, double area2, double error2)
		{
			
			int i_max = i ;
			int i_new = size ;

			int new_level = level[i_max] + 1;

			/* append the newly-created intervals to the list */
			  
			if (error2 > error1)
			{
				alist[i_max] = a2;        /* blist[maxerr] is already == b2 */
				rlist[i_max] = area2;
				elist[i_max] = error2;
				level[i_max] = new_level;
			      
				alist[i_new] = a1;
				blist[i_new] = b1;
				rlist[i_new] = area1;
				elist[i_new] = error1;
				level[i_new] = new_level;
			}
			else
			{
				blist[i_max] = b1;        /* alist[maxerr] is already == a1 */
				rlist[i_max] = area1;
				elist[i_max] = error1;
				level[i_max] = new_level;
			      
				alist[i_new] = a2;
				blist[i_new] = b2;
				rlist[i_new] = area2;
				elist[i_new] = error2;
				level[i_new] = new_level;
			}
			  
			size ++;

			if (new_level > maximum_level)
				maximum_level = new_level;

			Qpsrt () ;
		} // Update()
		
		void Qpsrt()
		{
			int last = size - 1;

			double errmax ;
			double errmin ;
			int ii, k, top;

			int i_nrmax = nrmax;
			int i_maxerr = order[i_nrmax] ;
			  
			/* Check whether the list contains more than two error estimates */

			if (last < 2) 
			{
				order[0] = 0 ;
				order[1] = 1 ;
				i = i_maxerr ;
				return ;
			}

			errmax = elist[i_maxerr] ;

			/* This part of the routine is only executed if, due to a difficult
				integrand, subdivision increased the error estimate. In the normal
				case the insert procedure should start after the nrmax-th largest
				error estimate. */

			while (i_nrmax > 0 && errmax > elist[order[i_nrmax - 1]]) 
			{
				order[i_nrmax] = order[i_nrmax - 1] ;
				i_nrmax-- ;
			} 

			/* Compute the number of elements in the list to be maintained in
				descending order. This number depends on the number of
				subdivisions still allowed. */
			  
			if(last < (limit/2 + 2)) 
				top = last ;
			else
				top = limit - last + 1;
			  
			/* Insert errmax by traversing the list top-down, starting
				comparison from the element elist(order(i_nrmax+1)). */
			  
			ii = i_nrmax + 1 ;
			  
			/* The order of the tests in the following line is important to
				prevent a segmentation fault */

			while (ii < top && errmax < elist[order[ii]])
			{
				order[ii - 1] = order[ii] ;
				ii ++ ;
			}
			  
			order[ii - 1] = i_maxerr ;
			  
			/* Insert errmin by traversing the list bottom-up */
			  
			errmin = elist[last] ;
			  
			k = top - 1 ;
			  
			while (k > ii - 2 && errmin >= elist[order[k]])
			{
				order[k+1] = order[k] ;
				k-- ;
			}
			  
			order[k+1] = last ;

			/* Set i_max and e_max */

			i_maxerr = order[i_nrmax] ;
			  
			i = i_maxerr ;
			nrmax = i_nrmax ;
		} // Qpsrt()

		public bool LargeInterval ()
		{
			if (level[i] < maximum_level)
				return true ;
			else
				return false ;
		} // LargeInterval()

		public bool IncreaseNrmax ()
		{
			int k;
			int id = nrmax;
			int jupbnd;

			int last = size - 1 ;

			if (last > (1 + limit / 2))
			{
				jupbnd = limit + 1 - last;
			}
			else
			{
				jupbnd = last;
			}
  
			for (k = id; k <= jupbnd; k++)
			{
				int i_max = order[nrmax];
      
				i = i_max ;

				if (level[i_max] < maximum_level)
				{
					return true;
				}

				nrmax ++;

			}
			return false;
		} // IncreaseNrmax()

		public void ResetNrmax ()
		{
			nrmax = 0;
			i = order[0] ;
		} // ResetNrmax()

		public double SumResults ()
		{
			int n = size;

			int k;
			double result_sum = 0;

			for (k = 0; k < n; k++)
				result_sum += rlist[k];

			return result_sum;
		} // SumResult()



		#endregion
	} // PZMath_integration_workspace

	/// <summary>
	/// extrapolation table
	/// </summary>
	public class ExtrapolationTable
	{
		public int n;
		public double[] rlist2 = null;
		public int nres;
		public double[] res3la = null;

		public ExtrapolationTable()
		{
			n = 0;
			nres = 0;
			rlist2 = new double [52];
			res3la = new double [3];
		} // ExtrapolationTable()
		
		#region qelq comments
		/* integration/qelg.c
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
		#endregion
		#region qelq methods
		public void InitialiseTable()
		{
			n = 0;
			nres = 0;
		} // InitialiseTable()

		public void InitialiseTable(double y)
		{
			n = 0;
			nres = 0;
			rlist2[0] = y;
		} // InitialiseTable()

		public void AppendTable (double y)
		{
			rlist2[n] = y;
			n ++;
		} // AppendTable()
		
		public void Qelg (ref double result, ref double abserr)
		{
			double [] epstab = rlist2;
			double [] local_res3la = res3la;
			int nn = n - 1;

			double current = epstab[nn];

			double absolute = PZMath_machine.PZMath_DBL_MAX;
			double relative = 5 * PZMath_machine.PZMath_DBL_EPSILON * System.Math.Abs (current);

			int newelm = nn / 2;
			int n_orig = nn;
			int n_final = nn;
			int i;

			int nres_orig = nres;

			result = current;
			abserr = PZMath_machine.PZMath_DBL_MAX;;

			if (nn < 2)
			{
				result = current;
				abserr = System.Math.Max (absolute, relative);
				return;
			}

			epstab[nn + 2] = epstab[nn];
			epstab[nn] = PZMath_machine.PZMath_DBL_MAX;

			for (i = 0; i < newelm; i++)
			{
				double res = epstab[nn - 2 * i + 2];
				double e0 = epstab[nn - 2 * i - 2];
				double e1 = epstab[nn - 2 * i - 1];
				double e2 = res;

				double e1abs = System.Math.Abs (e1);
				double delta2 = e2 - e1;
				double err2 = System.Math.Abs (delta2);
				double tol2 = System.Math.Max (System.Math.Abs (e2), e1abs) * PZMath_machine.PZMath_DBL_EPSILON;
				double delta3 = e1 - e0;
				double err3 = System.Math.Abs (delta3);
				double tol3 = System.Math.Max (e1abs, System.Math.Abs (e0)) * PZMath_machine.PZMath_DBL_EPSILON;

				double e3, delta1, err1, tol1, ss;

				if (err2 <= tol2 && err3 <= tol3)
				{
					/* If e0, e1 and e2 are equal to within machine accuracy,
						convergence is assumed.  */

					result = res;
					absolute = err2 + err3;
					relative = 5 * PZMath_machine.PZMath_DBL_EPSILON * System.Math.Abs (res);
					abserr = System.Math.Max (absolute, relative);
					return;
				}

				e3 = epstab[nn - 2 * i];
				epstab[nn - 2 * i] = e1;
				delta1 = e1 - e3;
				err1 = System.Math.Abs (delta1);
				tol1 = System.Math.Max (e1abs, System.Math.Abs (e3)) * PZMath_machine.PZMath_DBL_EPSILON;

				/* If two elements are very close to each other, omit a part of
					the table by adjusting the value of n */

				if (err1 <= tol1 || err2 <= tol2 || err3 <= tol3)
				{
					n_final = 2 * i;
					break;
				}

				ss = (1 / delta1 + 1 / delta2) - 1 / delta3;

				/* Test to detect irregular behaviour in the table, and
					eventually omit a part of the table by adjusting the value of
					n. */

				if (System.Math.Abs (ss * e1) <= 0.0001)
				{
					n_final = 2 * i;
					break;
				}

				/* Compute a new element and eventually adjust the value of
					result. */

				res = e1 + 1 / ss;
				epstab[nn - 2 * i] = res;

				{
					double error = err2 + System.Math.Abs (res - e2) + err3;

					if (error <= abserr)
					{
						abserr = error;
						result = res;
					}
				}
			}

			/* Shift the table */

			{
				int limexp = 50 - 1;

				if (n_final == limexp)
				{
					n_final = 2 * (limexp / 2);
				}
			}

			if (n_orig % 2 == 1)
			{
				for (i = 0; i <= newelm; i++)
					epstab[1 + i * 2] = epstab[i * 2 + 3];
			}
			else
			{
				for (i = 0; i <= newelm; i++)
					epstab[i * 2] = epstab[i * 2 + 2];
			}

			if (n_orig != n_final)
			{
				for (i = 0; i <= n_final; i++)
					epstab[i] = epstab[n_orig - n_final + i];
			}

			n = n_final + 1;

			if (nres_orig < 3)
			{
				local_res3la[nres_orig] = result;
				abserr = PZMath_machine.PZMath_DBL_MAX;
			}
			else
			{                           /* Compute error estimate */
				abserr = (System.Math.Abs (result - local_res3la[2]) + System.Math.Abs (result - local_res3la[1])
							+ System.Math.Abs (result - local_res3la[0]));

				local_res3la[0] = local_res3la[1];
				local_res3la[1] = local_res3la[2];
				local_res3la[2] = result;
			}

			/* In QUADPACK the variable table->nres is incremented at the top of
				qelg, so it increases on every call. This leads to the array
				res3la being accessed when its elements are still undefined, so I
				have moved the update to this point so that its value more
				useful. */

			nres = nres_orig + 1;  

			abserr = System.Math.Max (abserr, 5 * PZMath_machine.PZMath_DBL_EPSILON * System.Math.Abs (result));
			return;
		} // qelg()
		#endregion

	} // ExtrapolationTable

}
