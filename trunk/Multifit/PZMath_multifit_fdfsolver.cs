// PZ Math lmder
//
// -- numerical multi non linear regression
//
// -- Class PZMath_multifit_fdfsolver_type and Class PZMath_multifit_fdfsolver

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
// 24.11.2005


using System;

namespace eee.Sheffield.PZ.Math
{
	public delegate int multifit_fdfsolver_alloc (object state, int n, int p);
	public delegate int multifit_fdfsolver_init (object state, PZMath_multifit_function_fdf fdf, PZMath_vector x, PZMath_vector f, PZMath_matrix J, PZMath_vector dx);
	public delegate int multifit_fdfsolver_iterate (object state, PZMath_multifit_function_fdf fdf, PZMath_vector x, PZMath_vector f, PZMath_matrix J, PZMath_vector dx);
	public delegate void multifit_fdfsolver_free (object state);

	/// <summary>
	/// multifit nonline f df solver type
	/// </summary>
	public class PZMath_multifit_fdfsolver_type
	{
		public string name = "";
		public int size = 0;		// it is useless in C#
		public multifit_fdfsolver_alloc alloc = null;
		public multifit_fdfsolver_init init = null;
		public multifit_fdfsolver_iterate iterate = null;
		public multifit_fdfsolver_free free = null;	// it is useless in C#

		public PZMath_multifit_fdfsolver_type() {}

		public void Init (string inputname)
		{
			switch (inputname)
			{
				case "PZMath_multifit_fdfsolver_lmsder":
					// lmsder solver
				{
					name = "lmsder";
					alloc = new multifit_fdfsolver_alloc(lmsder.Alloc);
					init = new multifit_fdfsolver_init(lmsder.Init);
					iterate = new multifit_fdfsolver_iterate(lmsder.Iterate);
					break;
				}
			}
		} // init()

	} // PZMath_multifit_fdfsolver_type

	/// <summary>
	/// multifit nonline f df solver
	/// </summary>
	public class PZMath_multifit_fdfsolver
	{
		public PZMath_multifit_fdfsolver_type type = null;
		public PZMath_multifit_function_fdf fdf = null;
		public PZMath_vector x = null;
		public PZMath_vector f = null;
		public PZMath_matrix J = null;
		public PZMath_vector dx = null;
		public object state = null;
		
		public PZMath_multifit_fdfsolver() {}

		public int Alloc (PZMath_multifit_fdfsolver_type T, int n, int p)
		{
			x = new PZMath_vector(p);
			f = new PZMath_vector(n);
			J = new PZMath_matrix(n, p);
			dx = new PZMath_vector(p);
			type = T;
			state = new lmder_state_t();
			return (this.type.alloc) (this.state, n, p);
		} // Alloc()

		public int Init (PZMath_multifit_function_fdf f, PZMath_vector x)
		{
			if (this.f.Size != f.n)
			{
				PZMath_errno.ERROR("function size does not match solver", PZMath_errno.PZMath_EBADLEN);
			}
			if (this.x.Size != x.Size)
			{
				PZMath_errno.ERROR("vector length does not match solver", PZMath_errno.PZMath_EBADLEN);
			}
			this.fdf = f;
			this.x = x;
			return (this.type.init) (this.state, this.fdf, this.x, this.f, this.J, this.dx);
		} // init()

		public int Iterate ()
		{
			return (this.type.iterate) (this.state, this.fdf, this.x, this.f, this.J, this.dx);
		} // iterate()

		public double FIT(int i)
		{
			if (i < 0 || i >= x.Size)
				PZMath_errno.ERROR("PZMath_multifit_fdfsolver::FIT(), Index Out of Range", PZMath_errno.PZMath_EFAILED);
			return x[i];
		} // FIT()

		public double ERR(int i, PZMath_matrix covar)
		{
			if (i < 0 || i >= covar.RowCount)
				PZMath_errno.ERROR("PZMath_multifit_fdfsolver::ERR(), Index Out of Range", PZMath_errno.PZMath_EFAILED);
            return System.Math.Sqrt(covar[i, i]);
		} // ERR()

		public void PrintState(int iter)
			// screen output state information
		{
			System.Console.Write("iter: " + iter);
			System.Console.Write(" x = ");
			for (int i = 0; i <x.Size; i ++)
				System.Console.Write(x[i] + " ");
			System.Console.Write("|f(x)| = " + PZMath_blas.Dnrm2(f));
			System.Console.WriteLine(" ");
		} // PrintState()

		public int TestDelta (PZMath_vector dx, PZMath_vector x, double epsabs, double epsrel)
		{
			int i;
			int ok = 1;
			int n = x.Size ;

			if (epsrel < 0.0)
			{
				PZMath_errno.ERROR ("relative tolerance is negative", PZMath_errno.PZMath_EBADTOL);
			}

			for (i = 0 ; i < n ; i++)
			{
				double xi = x[i];
				double dxi = dx[i];
                double tolerance = epsabs + epsrel * System.Math.Abs(xi);

                if (System.Math.Abs(dxi) < tolerance)
					ok = 1;
				else
				{
					ok = 0;
					break;
				}
			}

			if (ok == 1)
				return PZMath_errno.PZMath_SUCCESS ;

			return PZMath_errno.PZMath_CONTINUE;
		} // TestDelta()

		public void Covar(double epsrel, PZMath_matrix covar)
			// calculate covar
		{
			double tolr;

			int i, j, k;
			int kmax = 0;
			
			PZMath_matrix r;
			PZMath_vector tau;
			PZMath_vector norm;
			PZMath_permutation perm;

			int m = J.RowCount;
			int n = J.ColumnCount ;
			  
			if (m < n) 
				PZMath_errno.ERROR ("PZMath_multifit_fdfsolver::Covar(), Jacobian be rectangular M x N with M >= N", PZMath_errno.PZMath_EBADLEN);
			if (covar.ColumnCount != covar.RowCount || covar.RowCount != n)
				PZMath_errno.ERROR ("PZMath_multifit_fdfsolver::Covar(), covariance matrix must be square and match second dimension of jacobian", PZMath_errno.PZMath_EBADLEN);
			r = new PZMath_matrix(m, n);
			tau = new PZMath_vector(n);
			perm = new PZMath_permutation(n);
			norm = new PZMath_vector(n);
			  
			{
				int signum = 0;
				r.MemCopyFrom(J);
				PZMath_linalg.QRPTDecomp(r, tau, perm, out signum, norm);
			}
			  
			  
			/* Form the inverse of R in the full upper triangle of R */

            tolr = epsrel * System.Math.Abs(r[0, 0]);

			for (k = 0 ; k < n ; k++)
			{
				double rkk = r[k, k];

                if (System.Math.Abs(rkk) <= tolr)
				{
					break;
				}

				r[k, k] = 1.0 / rkk;

				for (j = 0; j < k ; j++)
				{
					double t = r[j, k] / rkk;
					r[j, k] = 0.0;

					for (i = 0; i <= j; i++)
					{
						double rik = r[i, k];
						double rij = r[i, j];
			            r[i, k] = rik - t * rij;
					}
				}
				kmax = k;
			}

			/* Form the full upper triangle of the inverse of R^T R in the full
				upper triangle of R */

			for (k = 0; k <= kmax ; k++)
			{
				for (j = 0; j < k; j++)
				{
					double rjk = r[j, k];

					for (i = 0; i <= j ; i++)
					{
						double rij = r[i, j];
						double rik = r[i, k];
					
						r[i, j] = rij + rjk * rik;
					}
				}
			      
				{
					double t = r[k, k];

					for (i = 0; i <= k; i++)
					{
						double rik = r[i, k];
	
						r[i, k] = t * rik;
					}
				}
			}

			/* Form the full lower triangle of the covariance matrix in the
				strict lower triangle of R and in w */

			for (j = 0 ; j < n ; j++)
			{
				int pj = (int) perm[j];
			      
				for (i = 0; i <= j; i++)
				{
					int pi = (int) perm[i];

					double rij;

					if (j > kmax)
					{
						r[i, j] = 0.0;
						rij = 0.0 ;
					}
					else 
					{
						rij = r[i, j];
					}

					if (pi > pj)
					{
						r[pi, pj] = rij;
					} 
					else if (pi < pj)
					{
						r[pj, pi] = rij;
					}

				}
			      
				{ 
					double rjj = r[j, j];
					covar[pj, pj] = rjj;
				}
			}

			     
			/* symmetrize the covariance matrix */

			for (j = 0 ; j < n ; j++)
			{
				for (i = 0; i < j ; i++)
				{
					double rji = r[j, i];

					covar[j, i] = rji;
					covar[i, j] = rji;
				}
			}
		} // Covar()

		public double Chi()
			// return Chi value
		{
			return PZMath_blas.Dnrm2(f);
		} // Chi()

	} // PZMath_multifit_fdfsolver
}
