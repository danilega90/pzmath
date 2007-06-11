// PZ Math Multifit linear
//
// -- numerical multi linear regression
//
//
// GNU Copyright
/* multifit/multilinear.c
 * 
 * Copyright (C) 2000 Brian Gough
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
// 14.12.2005


using System;

namespace eee.Sheffield.PZ.Math
{
	public class PZMath_multifit_linear_workspace
	{
		public int n = 0; /* number of observations */
		public int p = 0; /* number of parameters */
		public PZMath_matrix A = null;
		public PZMath_matrix Q = null;
		public PZMath_matrix QSI = null;
		public PZMath_vector S = null;
		public PZMath_vector t = null;
		public PZMath_vector xt = null;
		public PZMath_vector D = null;

		public PZMath_multifit_linear_workspace() {}
	} // PZMath_multifit_linear_workspace

	/// <summary>
	/// multifit linear class
	/// </summary>
	public class PZMath_multifit_linear
	{
		public int n = 0; /* number of observations */
		public int p = 0; /* number of parameters */
		public PZMath_matrix A = null;
		public PZMath_matrix Q = null;
		public PZMath_matrix QSI = null;
		public PZMath_vector S = null;
		public PZMath_vector t = null;
		public PZMath_vector xt = null;
		public PZMath_vector D = null;

		public PZMath_multifit_linear() {}
		
		#region methods

		/// <summary>
		/// alloc parameters
		/// </summary>
		/// <param name="inputn">number of observations</param>
		/// <param name="inputp">number of parameters</param>
		public void Alloc (int inputn, int inputp)
		{
			n = inputn;                     /* number of observations */
			p = inputp;                     /* number of parameters */

			A = new PZMath_matrix(n, p);
			Q = new PZMath_matrix (p, p);
			QSI = new PZMath_matrix (p, p);
			S = new PZMath_vector (p);
			t = new PZMath_vector (n);
			xt = new PZMath_vector (p);	// calloc()
			D = new PZMath_vector (p);	// calloc()
		} // Alloc(int, int)

		public int Wlinear(PZMath_matrix X, PZMath_vector w, PZMath_vector y, PZMath_vector c, PZMath_matrix cov, 
				out double chisq)
		{
			int rank;
			int status  = WlinearSVD (X, w, y, PZMath_machine.PZMath_DBL_EPSILON, out rank, c, cov, out chisq);
			return status;
		} // Wlinear()


		public int WlinearSVD ( PZMath_matrix X, PZMath_vector w, PZMath_vector y, double tol, 
			out int rank, PZMath_vector c, PZMath_matrix cov, out double chisq)
		{
			if (X.RowCount != y.Size)
				PZMath_errno.ERROR
					("PZMath_multifit_linear::WlinearSVD(), number of observations in y does not match rows of matrix X",
					PZMath_errno.PZMath_EBADLEN);
			else if (X.ColumnCount != c.Size)
				PZMath_errno.ERROR ("PZMath_multifit_linear::WlinearSVD(), number of parameters c does not match columns of matrix X",
							PZMath_errno.PZMath_EBADLEN);
			else if (w.Size != y.Size)
				PZMath_errno.ERROR ("PZMath_multifit_linear::WlinearSVD(), number of weights does not match number of observations",
							PZMath_errno.PZMath_EBADLEN);
			else if (cov.RowCount != cov.ColumnCount)
				PZMath_errno.ERROR ("PZMath_multifit_linear::WlinearSVD(), covariance matrix is not square", PZMath_errno.PZMath_ENOTSQR);
			else if (c.Size != cov.RowCount)
				PZMath_errno.ERROR
					("PZMath_multifit_linear::WlinearSVD(), number of parameters does not match size of covariance matrix",
					PZMath_errno.PZMath_EBADLEN);
			else if (X.RowCount != n || X.ColumnCount != p)
				PZMath_errno.ERROR
					("PZMath_multifit_linear::WlinearSVD(), size of workspace does not match size of observation matrix",
					PZMath_errno.PZMath_EBADLEN);
			else
			{
				int nn = X.RowCount;
				int pp = X.ColumnCount;

				int i, j, p_eff;

				/* Scale X,  A = sqrt(w) X */
				A.MemCopyFrom(X);

				for (i = 0; i < nn; i++)
				{
					double wi = w[i];

					if (wi < 0)
						wi = 0;
					{
						PZMath_vector row = A.Row(i);
                        row.Scale(System.Math.Sqrt(wi));
					}
				}
                
				/* Balance the columns of the matrix A */

				PZMath_linalg.BalanceColumns(A, D);

				/* Decompose A into U S Q^T */

				PZMath_linalg.SVDecompMod (A, QSI, Q, S, xt);

				/* Solve sqrt(w) y = A c for c, by first computing t = sqrt(w) y */

				for (i = 0; i < nn; i++)
				{
					double wi = w[i];
					double yi = y[i];
					if (wi < 0)
						wi = 0;
                    t[i] = System.Math.Sqrt(wi) * yi;
				}

				PZMath_blas.Dgemv (1, 1.0, A, t, 0.0, xt);

				/* Scale the matrix Q,  Q' = Q S^-1 */
				QSI.MemCopyFrom(Q);

				{
					double alpha0 = S[0];
					p_eff = 0;
				    
					for (j = 0; j < pp; j++)
					{
						PZMath_vector column = QSI.Column(j);
						double alpha = S[j];

						if (alpha <= tol * alpha0) {
						alpha = 0.0;
						} else {
						alpha = 1.0 / alpha;
						p_eff++;
						}

						column.Scale(alpha);
					}

					rank = p_eff;
				}

				c.SetZero();

				/* Solution */

				PZMath_blas.Dgemv (0, 1.0, QSI, xt, 0.0, c);

				/* Unscale the balancing factors */

				c.Div(D);

				/* Form covariance matrix cov = (Q S^-1) (Q S^-1)^T */

				for (i = 0; i < pp; i++)
				{
					PZMath_vector row_i = QSI.Row(i);
					double d_i = D[i];

					for (j = i; j < pp; j++)
					{
						PZMath_vector row_j =  QSI.Row(j);
						double d_j = D[j];
						double s;

						PZMath_blas.Ddot (row_i, row_j, out s);

						cov[i, j] = s / (d_i * d_j);
						cov[j, i] = s / (d_i * d_j);
					}
				}

				/* Compute chisq, from residual r = y - X c */

				{
					double r2 = 0;

					for (i = 0; i < nn; i++)
					{
						double yi = y[i];
						double wi = w[i];
						PZMath_vector row = X.Row(i);
						double y_est, ri;
						PZMath_blas.Ddot (row, c, out y_est);
						ri = yi - y_est;
						r2 += wi * ri * ri;
					}

					chisq = r2;
				}

				return PZMath_errno.PZMath_SUCCESS;
			} // end of else
			chisq = 0.0;
			rank = 0;
			return 0;
		} // WlinearSVD()

		#endregion
	} // PZMath_multifit_linear

}
