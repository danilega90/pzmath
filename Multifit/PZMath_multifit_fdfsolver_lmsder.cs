// PZ Math lmder
//
// -- numerical multi non linear regression
//
// -- Class PZMath_multifit_fdfsolver_lmsder

// GNU Copyright
/* multfit/lmder.c
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


using System;

namespace eee.Sheffield.PZ.Math
{
	public class lmder_state_t 
	{
		public int iter = 0;
		public double xnorm = 0.0d;
		public double fnorm = 0.0d;
		public double delta = 0.0d;
		public double par = 0.0d;
		public PZMath_matrix r = null;
		public PZMath_vector tau = null;
		public PZMath_vector diag = null;
		public PZMath_vector qtf = null;
		public PZMath_vector newton = null;
		public PZMath_vector gradient = null;
		public PZMath_vector x_trial = null;
		public PZMath_vector f_trial = null;
		public PZMath_vector df = null;
		public PZMath_vector sdiag = null;
		public PZMath_vector rptdx = null;
		public PZMath_vector w = null;
		public PZMath_vector work1 = null;
		public PZMath_permutation perm = null;

		public lmder_state_t() {}
	} // lmder_state_t

	public class lmsder
	{
		public lmsder() {}

		public static int Init
			(object vstate, PZMath_multifit_function_fdf fdf, PZMath_vector x, PZMath_vector f, PZMath_matrix J, PZMath_vector dx)
		{
			int status = Init_s (vstate, fdf, x, f, J, dx, 0);
			return status ;
		} // init()

		public static int Alloc (object vstate, int n, int p)
		{
			(vstate as lmder_state_t).r = new PZMath_matrix(n, p);
            (vstate as lmder_state_t).tau = new PZMath_vector(System.Math.Min(n, p));
			(vstate as lmder_state_t).diag = new PZMath_vector(p);
			(vstate as lmder_state_t).qtf = new PZMath_vector(n);
			(vstate as lmder_state_t).newton = new PZMath_vector(p);
			(vstate as lmder_state_t).gradient = new PZMath_vector(p);
			(vstate as lmder_state_t).x_trial = new PZMath_vector(p);
			(vstate as lmder_state_t).f_trial = new PZMath_vector(n);
			(vstate as lmder_state_t).df = new PZMath_vector(n);
			(vstate as lmder_state_t).sdiag = new PZMath_vector(p);
			(vstate as lmder_state_t).rptdx = new PZMath_vector(n);
			(vstate as lmder_state_t).w = new PZMath_vector(n);
			(vstate as lmder_state_t).work1 = new PZMath_vector(p);
			(vstate as lmder_state_t).perm = new PZMath_permutation(p);
			return PZMath_errno.PZMath_SUCCESS;
		} // Alloc()

		public static int Iterate (object state, PZMath_multifit_function_fdf fdf, PZMath_vector x, PZMath_vector f, PZMath_matrix J, PZMath_vector dx)
		{
			int status = Iterate_s (state, fdf, x, f, J, dx, 1);
			return status;
		} // iterate()

		public static int Free (object state)
		{
			return PZMath_errno.PZMath_SUCCESS;
		} // free()

		
		#region lmsder static methods
		private static int Init_s 
			(object vstate, PZMath_multifit_function_fdf fdf, PZMath_vector x, PZMath_vector f, PZMath_matrix J, PZMath_vector dx, int scale)
		{
			lmder_state_t state = (vstate as lmder_state_t);

			//PZMath_matrix r = state.r;
			//PZMath_vector tau = state.tau;
			//PZMath_vector diag = state.diag;
			//PZMath_vector work1 = state.work1;
			//PZMath_permutation perm = state.perm;

			int signum;
			fdf.FN_EVAL_F_DF (x, f, J);


			state.par = 0;
			state.iter = 1;
			state.fnorm = Enorm (f);
			dx.SetAll(0.0d);

			/* store column norms in diag */

			if (scale == 0)
			{
				ComputeDiag (J, state.diag);
			}
			else
			{
				state.diag.SetAll(1.0d);
			}

			/* set delta to 100 |D x| or to 100 if |D x| is zero */

			state.xnorm = ScaledEnorm (state.diag, x);
			state.delta = ComputeDelta (state.diag, x);

			/* Factorize J into QR decomposition */

			state.r.MemCopyFrom(J);
			PZMath_linalg.QRPTDecomp (state.r, state.tau, state.perm, out signum, state.work1);

			return PZMath_errno.PZMath_SUCCESS;
		} // init_s() 

		private static int Iterate_s
			(object vstate, PZMath_multifit_function_fdf fdf, PZMath_vector x, PZMath_vector f, PZMath_matrix J, PZMath_vector dx, int scale)
		{
			lmder_state_t state = (vstate as lmder_state_t);
			/*
			PZMath_matrix r = state.r;
			PZMath_vector tau = state.tau;
			PZMath_vector diag = state.diag;
			PZMath_vector qtf = state.qtf;
			PZMath_vector x_trial = state.x_trial;
			PZMath_vector f_trial = state.f_trial;
			PZMath_vector rptdx = state.rptdx;
			PZMath_vector newton = state.newton;
			PZMath_vector gradient = state.gradient;
			PZMath_vector sdiag = state.sdiag;
			PZMath_vector w = state.w;
			PZMath_vector work1 = state.work1;
			PZMath_permutation perm = state.perm;
			*/

			double prered, actred;
			double pnorm, fnorm1, fnorm1p, gnorm;
			double ratio;
			double dirder;

			int iter = 0;

			double p1 = 0.1, p25 = 0.25, p5 = 0.5, p75 = 0.75, p0001 = 0.0001;

			if (state.fnorm == 0.0) 
			{
				return PZMath_errno.PZMath_SUCCESS;
			}

			/* Compute qtf = Q^T f */
			//f.CopyTo(qtf);
			state.qtf.MemCopyFrom(f);
			PZMath_linalg.QR_QTvec(state.r, state.tau, state.qtf);

			/* Compute norm of scaled gradient */

			ComputeGradientDirection (state.r, state.perm, state.qtf, state.diag, state.gradient);

			{ 
				int iamax = PZMath_blas.Idamax(state.gradient);
                gnorm = System.Math.Abs(state.gradient[iamax] / state.fnorm);
			}

				/* Determine the Levenberg-Marquardt parameter */
lm_iteration:
			iter++ ;
			{
				int status = Lmpar (state.r, state.perm, state.qtf, state.diag, 
					state.delta, ref state.par, state.newton, state.gradient, state.sdiag, 
					dx, state.w);
				if (status > 0)
					return status;
			}

			/* Take a trial step */

			dx.Scale(-1.0d); /* reverse the step to go downhill */
			
			ComputeTrialStep (x, dx, state.x_trial);

			pnorm = ScaledEnorm (state.diag, dx);

			if (state.iter == 1)
			{
				if (pnorm < state.delta)
				{
					System.Diagnostics.Debug.WriteLine("set delta = pnorm = " + pnorm);
					state.delta = pnorm;
				}
			}

				/* Evaluate function at x + p */
				/* return immediately if evaluation raised error */
			{
				int status = fdf.FN_EVAL_F (state.x_trial, state.f_trial);
				if (status > 0)
					return status;
			}

				fnorm1 = Enorm (state.f_trial);

				/* Compute the scaled actual reduction */

				actred = ComputeActualReduction (state.fnorm, fnorm1);
				System.Diagnostics.Debug.WriteLine
					("lmiterate: fnorm = " + state.fnorm + " fnorm1 = " + fnorm1 + " actred = " + actred);

				/* Compute rptdx = R P^T dx, noting that |J dx| = |R P^T dx| */
				ComputeRptdx (state.r, state.perm, dx, state.rptdx);

				fnorm1p = Enorm (state.rptdx);

				/* Compute the scaled predicted reduction = |J dx|^2 + 2 par |D dx|^2 */

			{ 
				double t1 = fnorm1p / state.fnorm;
                double t2 = (System.Math.Sqrt(state.par) * pnorm) / state.fnorm;
    
				prered = t1 * t1 + t2 * t2 / p5;
				dirder = -(t1 * t1 + t2 * t2);
			}

				/* compute the ratio of the actual to predicted reduction */
				if (prered > 0)
					ratio = actred / prered;
				else
					ratio = 0;

				System.Diagnostics.Debug.WriteLine
					("lmiterate: prered = " + prered + " dirder = " + dirder + " ratio = " + ratio);
				/* update the step bound */

				if (ratio > p25)
				{
					System.Diagnostics.Debug.WriteLine("ratio > p25");

					if (state.par == 0 || ratio >= p75)
					{
						state.delta = pnorm / p5;
						state.par *= p5;
						System.Diagnostics.Debug.WriteLine
							("updated step bounds: delta = "+ state.delta + " par = " + state.par);
					}
				}
				else
				{
					double temp = (actred >= 0) ? p5 : p5 * dirder / (dirder + p5 * actred);
					System.Diagnostics.Debug.WriteLine("ratio < p25");

					if (p1 * fnorm1 >= state.fnorm || temp < p1 ) 
					{
						temp = p1;
					}

                    state.delta = temp * System.Math.Min(state.delta, pnorm / p1);

					state.par /= temp;
					System.Diagnostics.Debug.WriteLine
						("updated step bounds: delta = " + state.delta + " par = " + state.par);
				}


				/* test for successful iteration, termination and stringent tolerances */

				if (ratio >= p0001)
				{
					x.MemCopyFrom(state.x_trial);
					f.MemCopyFrom(state.f_trial);

					/* return immediately if evaluation raised error */
				{
					int status = fdf.FN_EVAL_DF (state.x_trial, J);
					if (status > 0)
						return status;
				}

				/* wa2_j  = diag_j * x_j */
				state.xnorm = ScaledEnorm(state.diag, x);
				state.fnorm = fnorm1;
				state.iter ++;

				/* Rescale if necessary */

				if (scale > 0)
					UpdateDiag (J, state.diag);

				{
					int signum;
					state.r.MemCopyFrom(J);
					PZMath_linalg.QRPTDecomp(state.r, state.tau, state.perm, out signum, state.work1);
				}
      
					return PZMath_errno.PZMath_SUCCESS;
				}
                else if (System.Math.Abs(actred) <= PZMath_machine.PZMath_DBL_EPSILON 
					&& prered <= PZMath_machine.PZMath_DBL_EPSILON 
					&& p5 * ratio <= 1.0)
				{
					return PZMath_errno.PZMath_ETOLF;
				}
				else if (state.delta <= PZMath_machine.PZMath_DBL_EPSILON * state.xnorm)
				{
					return PZMath_errno.PZMath_ETOLX;
				}
				else if (gnorm <= PZMath_machine.PZMath_DBL_EPSILON )
				{
					return PZMath_errno.PZMath_ETOLG;
				}
				else if (iter < 10)
				{
					/* Repeat inner loop if unsuccessful */
					goto lm_iteration;
				}

				return PZMath_errno.PZMath_CONTINUE;
		} // Iterate_s()
		#endregion

		#region lmder util methods
		private static double Enorm (PZMath_vector f)
		{
			return PZMath_blas.Dnrm2(f);
		} // Enorm()

		private static void ComputeDiag (PZMath_matrix J, PZMath_vector diag)
		{
			int i, j;
			int n = J.RowCount;
			int p = J.ColumnCount;

			for (j = 0; j < p; j++)
			{
				double sum = 0;
				for (i = 0; i < n; i++)
				{
					double Jij = J[i, j];
					sum += Jij * Jij;
				}
				if (sum == 0)
					sum = 1.0;
                diag[j] = System.Math.Sqrt(sum);
			}
		} // ComputeDiag()

		private static double ScaledEnorm (PZMath_vector d, PZMath_vector f)
		{
			double e2 = 0;
			int i;
			int n = f.Size;
			for (i = 0; i < n; i++)
			{
				double fi = f[i];
				double di = d[i];
				double u = di * fi;
				e2 += u * u;
			}
            return System.Math.Sqrt(e2);
		} // ScaledEnorm()

		private static double ComputeDelta (PZMath_vector diag, PZMath_vector x)
		{
			double Dx = ScaledEnorm (diag, x);
			double factor = 100;  /* generally recommended value from MINPACK */

			return (Dx > 0) ? factor * Dx : factor;
		} // ComputeDelta()
		#endregion

		#region lmpar comments
/* multifit/lmpar.c
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
		#region lmpar methods
		private static int Lmpar (PZMath_matrix r, PZMath_permutation perm, PZMath_vector qtf, 
			PZMath_vector diag, double delta, ref double par_inout, PZMath_vector newton, 
			PZMath_vector gradient, PZMath_vector sdiag, PZMath_vector x, PZMath_vector w)
		{
			double dxnorm, gnorm, fp, fp_old, par_lower, par_upper, par_c;
			double par = par_inout;
		
			int iter = 0;
			System.Diagnostics.Debug.WriteLine("PZMath_multifit_fdfsolver_lmsder::lmpar(): ENTERING lmpar\n");
			
			ComputeNewtonDirection (r, perm, qtf, newton);
			System.Diagnostics.Debug.WriteLine("newton = ");
			newton.DebugWriteLine();
			System.Diagnostics.Debug.WriteLine("diag = ");
			diag.DebugWriteLine();
	
			/* Evaluate the function at the origin and test for acceptance of
			the Gauss-Newton direction. */
			dxnorm = ScaledEnorm(diag, newton);
			fp = dxnorm - delta;

			System.Diagnostics.Debug.WriteLine("dxnorm = " + dxnorm + " delta = " + delta + " fp =" + fp);

			if (fp <= 0.1 * delta)
			{
				x.MemCopyFrom(newton);
				System.Diagnostics.Debug.WriteLine("took newton (fp = " + fp + ", delta = " + delta);

				par_inout = 0;

				return PZMath_errno.PZMath_SUCCESS;
			}

			ComputeNewtonBound (r, newton, dxnorm, perm, diag, w);
			{
				double wnorm = Enorm (w);
				double phider = wnorm * wnorm;

				/* w == zero if r rank-deficient, 
				then set lower bound to zero form MINPACK, lmder.f 
				Hans E. Plesser 2002-02-25 (hans.plesser@itf.nlh.no) */
				if ( wnorm > 0 )
					par_lower = fp / (delta * phider);
				else
					par_lower = 0.0;
			}
			System.Diagnostics.Debug.WriteLine("par = " + par);
			System.Diagnostics.Debug.WriteLine("par_lower = " + par_lower);

			ComputeGradientDirection (r, perm, qtf, diag, gradient);
	
			gnorm = Enorm (gradient);
			System.Diagnostics.Debug.WriteLine("gradient = ");
			gradient.DebugWriteLine();
			System.Diagnostics.Debug.WriteLine("gnorm = " + gnorm);

			par_upper =  gnorm / delta;

			if (par_upper == 0)
			{
                par_upper = PZMath_machine.PZMath_DBL_MIN / System.Math.Min(delta, 0.1);
			}
			System.Diagnostics.Debug.WriteLine("par_upper = " + par_upper);

			if (par > par_upper)
			{
				System.Diagnostics.Debug.WriteLine("set par to par_upper");
				par = par_upper;
			}
			else if (par < par_lower)
			{
				System.Diagnostics.Debug.WriteLine("set par to par_lower");
				par = par_lower;
			}
			if (par == 0)
			{
				par = gnorm / dxnorm;
				System.Diagnostics.Debug.WriteLine("set par to gnorm/dxnorm");
			}
			/* Beginning of iteration */

iteration:

			iter ++;
			System.Diagnostics.Debug.WriteLine("lmpar iteration = " + iter);

// check here!
//#ifdef BRIANSFIX
			/* Seems like this is described in the paper but not in the MINPACK code */

			if (par < par_lower || par > par_upper)
                par = System.Math.Max(0.001 * par_upper, System.Math.Sqrt(par_lower * par_upper));

			/* Evaluate the function at the current value of par */

			if (par == 0)
			{
                par = System.Math.Max(0.001 * par_upper, PZMath_machine.PZMath_DBL_MIN);
				System.Diagnostics.Debug.WriteLine("par = 0 set par to = " + par);
			}

			//Compute the least squares solution of [ R P x - Q^T f, sqrt(par) D x] for A = Q R P^T
		
			System.Diagnostics.Debug.WriteLine("Calling qrsolv with par = " + par);
			{
                double sqrt_par = System.Math.Sqrt(par);
				QRSolv (r, perm, sqrt_par, diag, qtf, x, sdiag, w);
			}

			dxnorm = ScaledEnorm (diag, x);

			fp_old = fp;

			fp = dxnorm - delta;
			System.Diagnostics.Debug.WriteLine("After QRSolv dxnorm = " + dxnorm + " delta = " + delta + " fp = " + fp);
			System.Diagnostics.Debug.WriteLine("sdiag = ");
			sdiag.DebugWriteLine();
			System.Diagnostics.Debug.WriteLine("x = ");
			x.DebugWriteLine();
			System.Diagnostics.Debug.WriteLine("r = ");
			r.DebugWriteLine();
			
			/* If the function is small enough, accept the current value of par */

            if (System.Math.Abs(fp) <= 0.1 * delta)
				goto line220;

			if (par_lower == 0 && fp <= fp_old && fp_old < 0)
				goto line220;

			/* Check for maximum number of iterations */

			if (iter == 10)
				goto line220;

			/* Compute the Newton correction */
			ComputeNewtonCorrection (r, sdiag, perm, x, dxnorm, diag, w);
			
			System.Diagnostics.Debug.WriteLine("Newton correction = ");
			w.DebugWriteLine();

			{
				double wnorm = Enorm (w);
				par_c = fp / (delta * wnorm * wnorm);
			}
			
			System.Diagnostics.Debug.WriteLine("fp = " + fp);
			System.Diagnostics.Debug.WriteLine("par_lower = " + par_lower);
			System.Diagnostics.Debug.WriteLine("par_upper = " + par_upper);
			System.Diagnostics.Debug.WriteLine("par_c = " + par_c);
			
			/* Depending on the sign of the function, update par_lower or par_upper */

			if (fp > 0)
			{
				if (par > par_lower)
				{
					par_lower = par;
					System.Diagnostics.Debug.WriteLine("fp > 0 set par_lower = par = " + par);
				}
			}
			else if (fp < 0)
			{
				if (par < par_upper)
				{
					System.Diagnostics.Debug.WriteLine("fp < 0 set par_upper = par = " + par);
					par_upper = par;
				}
			}

			/* Compute an improved estimate for par */

			System.Diagnostics.Debug.WriteLine("improved estimate par = MAX(" + par_lower + ", " + (par + par_c) + ")");

            par = System.Math.Max(par_lower, par + par_c);
			System.Diagnostics.Debug.WriteLine("improved estimate par = " + par);

			goto iteration;

line220:
						
			System.Diagnostics.Debug.WriteLine("LEAVING lmpar, par = " + par);
			par_inout = par;
			return PZMath_errno.PZMath_SUCCESS;
		} // Lmpar()
		#endregion

		#region static method
		private static void ComputeGradientDirection 
			(PZMath_matrix r, PZMath_permutation p, PZMath_vector qtf, PZMath_vector diag, PZMath_vector g)
		{
			int n = r.ColumnCount;
			int i, j;

			for (j = 0; j < n; j++)
			{
				double sum = 0;

				for (i = 0; i <= j; i++)
				{
					sum += r[i, j] * qtf[i];
				}

				{
				int pj = (int) p[j];
				double dpj = diag[pj];
				g[j] = sum / dpj;
				}
			}
		} // ComputeGradientDirection()

		private static void ComputeNewtonDirection 
			(PZMath_matrix r, PZMath_permutation perm, PZMath_vector qtf, PZMath_vector x)
		{
			/* Compute and store in x the Gauss-Newton direction. If the
			Jacobian is rank-deficient then obtain a least squares
			solution. */

			int n = r.ColumnCount;
			int i, j, nsing;

			for (i = 0 ; i < n ; i++)
			{
				double qtfi = qtf[i];
				x[i] = qtfi;
			}

			nsing = CountNsing (r);
			System.Diagnostics.Debug.WriteLine("nsing = " + nsing);
			r.DebugWriteLine();
			x.DebugWriteLine();

			for (i = nsing; i < n; i++)
				x[i] = 0.0;

			if (nsing > 0)
			{			
				for (j = nsing; j > 0 && (j-- > 0);)
				{
					double rjj = r[j, j];
					double temp = x[j] / rjj;
					x[j] = temp;          
					for (i = 0; i < j; i++)
					{
						double rij = r[i, j];
						double xi = x[i];
						x[i] = xi - rij * temp;
					}
				}
			}
			perm.InversePermute(x);
		} // ComputeNewtonDirection()

		private static void ComputeNewtonBound (PZMath_matrix r, PZMath_vector x, double dxnorm,
			PZMath_permutation perm, PZMath_vector diag, PZMath_vector w)
		{
			/* If the jacobian is not rank-deficient then the Newton step
			provides a lower bound for the zero of the function. Otherwise
			set this bound to zero. */

			int n = r.ColumnCount;
			int i, j;

			int nsing = CountNsing (r);

			if (nsing < n)
			{
				w.SetZero();
				return;
			}

			for (i = 0; i < n; i++)
			{
				int pi = (int) perm[i];

				double dpi = diag[pi];
				double xpi = x[pi];
				w[i] = dpi * (dpi * xpi / dxnorm);
			}
			for (j = 0; j < n; j++)
			{
				double sum = 0;

				for (i = 0; i < j; i++)
				{
					sum += r[i, j] * w[i];
				}

				{
					double rjj = r[j, j];
					double wj = w[j];
					w[j] = (wj - sum) / rjj;
				}
			}
		} // ComputeNewtonBound()

		private static int CountNsing (PZMath_matrix r)
		{
			/* Count the number of nonsingular entries. Returns the index of the
			first entry which is singular. */
			int n = r.ColumnCount;
			int i;

			for (i = 0; i < n; i++)
			{
				double rii = r[i, i];
				if (rii == 0)
					break;
			}
			return i;
		} // CountNsing()
		private static void ComputeTrialStep (PZMath_vector x, PZMath_vector dx, 
			PZMath_vector x_trial)
		{
			int i, N = x.Size;

			for (i = 0; i < N; i++)
			{
				double pi = dx[i];
				double xi = x[i];
				x_trial[i] = xi + pi;
			}
		} // ComputeTrialStep()

		private static double ComputeActualReduction (double fnorm, double fnorm1)
		{
			double actred;

			if (0.1 * fnorm1 < fnorm)
			{
				double u = fnorm1 / fnorm;
				actred = 1 - u * u;
			}
			else
			{
				actred = -1;
			}
			return actred;
		} // ComputeActualReduction()

		private static void ComputeRptdx (PZMath_matrix r, PZMath_permutation p,PZMath_vector dx, PZMath_vector rptdx)
		{
			int i, j, N = dx.Size;

			for (i = 0; i < N; i++)
			{
				double sum = 0;

				for (j = i; j < N; j++)
				{
					int pj = (int) p[j];
					sum += r[i, j] * dx[pj];
				}
				rptdx[i] = sum;
			}
		} // ComputeRptdx()

		private static void UpdateDiag (PZMath_matrix J, PZMath_vector diag)
		{
			int i, j, n = diag.Size;

			for (j = 0; j < n; j++)
			{
				double cnorm, diagj, sum = 0;
				for (i = 0; i < n; i++)
				{
					double Jij = J[i, j];
					sum += Jij * Jij;
				}
				if (sum == 0)
					sum = 1.0;

                cnorm = System.Math.Sqrt(sum);
				diagj = diag[j];

				if (cnorm > diagj)
					diag[j] = cnorm;
				
			}
			
		} // UpdateDiag()

		private static void ComputeNewtonCorrection (PZMath_matrix r, PZMath_vector sdiag, 
			PZMath_permutation p, PZMath_vector x,double dxnorm, PZMath_vector diag, PZMath_vector w)
		{
			int n = r.ColumnCount;
			int i, j;

			for (i = 0; i < n; i++)
			{
				int pi = (int) p[i];

				double dpi = diag[pi];
				double xpi = x[pi];

				w[i] = dpi * (dpi * xpi) / dxnorm;
			}

			for (j = 0; j < n; j++)
			{
				double sj = sdiag[j];
				double wj = w[j];

				double tj = wj / sj;

				w[j] = tj;

				for (i = j + 1; i < n; i++)
				{
					double rij = r[i, j];
					double wi = w[i];

					w[i] = wi - rij * tj;
				}
			}
		} // ComputeNewtonCorrection()
		#endregion

		#region qrsolv comments
	/* This function computes the solution to the least squares system

   phi = [ A x =  b , lambda D x = 0 ]^2
    
   where A is an M by N matrix, D is an N by N diagonal matrix, lambda
   is a scalar parameter and b is a vector of length M.

   The function requires the factorization of A into A = Q R P^T,
   where Q is an orthogonal matrix, R is an upper triangular matrix
   with diagonal elements of non-increasing magnitude and P is a
   permuation matrix. The system above is then equivalent to

   [ R z = Q^T b, P^T (lambda D) P z = 0 ]

   where x = P z. If this system does not have full rank then a least
   squares solution is obtained.  On output the function also provides
   an upper triangular matrix S such that

   P^T (A^T A + lambda^2 D^T D) P = S^T S

   Parameters,
   
   r: On input, contains the full upper triangle of R. On output the
   strict lower triangle contains the transpose of the strict upper
   triangle of S, and the diagonal of S is stored in sdiag.  The full
   upper triangle of R is not modified.

   p: the encoded form of the permutation matrix P. column j of P is
   column p[j] of the identity matrix.

   lambda, diag: contains the scalar lambda and the diagonal elements
   of the matrix D

   qtb: contains the product Q^T b

   x: on output contains the least squares solution of the system

   wa: is a workspace of length N

   */
#endregion
		#region qrsolv
		private static int QRSolv (PZMath_matrix r, PZMath_permutation p, double lambda, 
			PZMath_vector diag, PZMath_vector qtb, PZMath_vector x, PZMath_vector sdiag, PZMath_vector wa)
		{
// check here!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!			
			int n = r.ColumnCount;
			int i, j, k, nsing;

			/* Copy r and qtb to preserve input and initialise s. In particular,
			save the diagonal elements of r in x */
			for (j = 0; j < n; j++)
			{
				double rjj = r[j, j];
				double qtbj = qtb[j];

				for (i = j + 1; i < n; i++)
				{
					double rji = r[j, i];
					r[i, j] = rji;
				}
				x[j] = rjj;
				wa[j] = qtbj;
			}

			/* Eliminate the diagonal matrix d using a Givens rotation */

			for (j = 0; j < n; j++)
			{
				double qtbpj;
				int pj = (int) p[j];
				double diagpj = lambda * diag[pj];

				if (diagpj == 0)
					continue;
				sdiag[j] = diagpj;

				for (k = j + 1; k < n; k++)
					sdiag[k] = 0.0d;

				/* The transformations to eliminate the row of d modify only a
				single element of qtb beyond the first n, which is initially
				zero */

				qtbpj = 0;

				for (k = j; k < n; k++)
				{
					/* Determine a Givens rotation which eliminates the
					appropriate element in the current row of d */

					double sine, cosine;

					double wak = wa[k];
					double rkk = r[k, k];
					double sdiagk = sdiag[k];

					if (sdiagk == 0)
						continue;

                    if (System.Math.Abs(rkk) < System.Math.Abs(sdiagk))
					{
						double cotangent = rkk / sdiagk;
                        sine = 0.5 / System.Math.Sqrt(0.25 + 0.25 * cotangent * cotangent);
						cosine = sine * cotangent;
					}
					else
					{
						double tangent = sdiagk / rkk;
                        cosine = 0.5 / System.Math.Sqrt(0.25 + 0.25 * tangent * tangent);
						sine = cosine * tangent;
					}

					/* Compute the modified diagonal element of r and the
						modified element of [qtb,0] */

					{
						double new_rkk = cosine * rkk + sine * sdiagk;
						double new_wak = cosine * wak + sine * qtbpj;
				        
						qtbpj = -sine * wak + cosine * qtbpj;
						r[k, k] = new_rkk;
						wa[k] = new_wak;
					}

					/* Accumulate the transformation in the row of s */

					for (i = k + 1; i < n; i++)
					{
						double rik = r[i, k];

						double sdiagi = sdiag[i];
					
            
						double new_rik = cosine * rik + sine * sdiagi;
						double new_sdiagi = -sine * rik + cosine * sdiagi;
				            
						r[i, k] = new_rik;
						sdiag[i] = new_sdiagi;
					}
				}

				/* Store the corresponding diagonal element of s and restore the
				corresponding diagonal element of r */

				{
					double rjj = r[j, j];
					double xj = x[j];
    
					sdiag[j] = rjj;
					r[j, j] = xj;
				}

			}

			/* Solve the triangular system for z. If the system is singular then
			obtain a least squares solution */

			nsing = n;

			for (j = 0; j < n; j++)
			{
				double sdiagj = sdiag[j];

				if (sdiagj == 0)
				{
					nsing = j;
					break;
				}
			}

			for (j = nsing; j < n; j++)
				wa[j] = 0.0d;

			for (k = 0; k < nsing; k++)
			{
				double sum = 0;

				j = (nsing - 1) - k;

				for (i = j + 1; i < nsing; i++)
				{
					sum += r[i, j] * wa[i];
				}

				{
					double waj = wa[j];
					double sdiagj = sdiag[j];

					wa[j] = (waj - sum) / sdiagj;
				}
			}

			/* Permute the components of z back to the components of x */

			for (j = 0; j < n; j++)
			{
				int pj = (int) p[j];
				double waj = wa[j];

				x[pj] = waj;
			}
			return PZMath_errno.PZMath_SUCCESS;
		} // QRSolv()
		#endregion
	}
}
