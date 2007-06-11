// eee.Sheffield.PZ.Math
//
// Copyright ?Ping Zou, 2007
// sg71.cherub@gmail.com


using System;

namespace eee.Sheffield.PZ.Math
{
	/// <summary>
	/// linalg methods
	/// </summary>
	
	public class PZMath_linalg
	{
		#region balance comments
		/* linalg/balance.c
		* 
		* Copyright (C) 2001, 2004 Brian Gough
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

		/* Balance a general matrix by scaling the columns
		 *
		 * B =  A D
		 *
		 * where D is a diagonal matrix
		 */
		#endregion
		#region balance methods
		public static int BalanceColumns (PZMath_matrix A, PZMath_vector D)
		{
			int N = A.ColumnCount;
			int j;

			if (D.Size != N)
			{
				PZMath_errno.ERROR ("PZMath_linalg::BalanceColumns(), length of D must match second dimension of A", PZMath_errno.PZMath_EINVAL);
			}
			
			D.SetAll(1.0);

			for (j = 0; j < N; j++)
			{
				PZMath_vector A_j = A.Column(j);
      
				double s = PZMath_blas.Dasum (A_j);
      
				double f = 1.0;

				if (s == 0.0 || !PZMath_sys.Finite(s))
				{
					D[j] = f;
					continue;
				}

				while (s > 1.0)
				{
					s /= 2.0;
					f *= 2.0;
				}
      
				while (s < 0.5)
				{
					s *= 2.0;
					f /= 2.0;
				}
      
				D[j] = f;

				if (f != 1.0)
					PZMath_blas.Dscal(1.0/f, A_j);
			}

			return PZMath_errno.PZMath_SUCCESS;
		} // BalanceColumns()
		#endregion

		#region QRPT comments
/* Factorise a general M x N matrix A into
 *
 *   A P = Q R
 *
 * where Q is orthogonal (M x M) and R is upper triangular (M x N).
 * When A is rank deficient, r = rank(A) < n, then the permutation is
 * used to ensure that the lower n - r rows of R are zero and the first
 * r columns of Q form an orthonormal basis for A.
 *
 * Q is stored as a packed set of Householder transformations in the
 * strict lower triangular part of the input matrix.
 *
 * R is stored in the diagonal and upper triangle of the input matrix.
 *
 * P: column j of P is column k of the identity matrix, where k =
 * permutation->data[j]
 *
 * The full matrix for Q can be obtained as the product
 *
 *       Q = Q_k .. Q_2 Q_1
 *
 * where k = MIN(M,N) and
 *
 *       Q_i = (I - tau_i * v_i * v_i')
 *
 * and where v_i is a Householder vector
 *
 *       v_i = [1, m(i+1,i), m(i+2,i), ... , m(M,i)]
 *
 * This storage scheme is the same as in LAPACK.  See LAPACK's
 * dgeqpf.f for details.
 * 
 */
		#endregion
		#region QRPT methods
		public static int QRPTDecomp 
			(PZMath_matrix A, PZMath_vector tau, PZMath_permutation p, out int signum, PZMath_vector norm)
		{
			int M = A.RowCount;
			int N = A.ColumnCount;
			signum = 0;

            if (tau.Size != System.Math.Min(M, N))
			{
				PZMath_errno.ERROR("size of tau must be MIN(M,N)", PZMath_errno.PZMath_EBADLEN);
			}
			else if (p.Size != N)
			{
				PZMath_errno.ERROR("permutation size must be N", PZMath_errno.PZMath_EBADLEN);
			}
			else if (norm.Size != N)
			{
				PZMath_errno.ERROR("norm size must be N", PZMath_errno.PZMath_EBADLEN);
			}
			else
			{
				int i;

				signum = 1;

				p.Init();  /* set to identity */

				/* Compute column norms and store in workspace */

				for (i = 0; i < N; i++)
				{
					PZMath_vector c = A.Column(i);
					double x = PZMath_blas.Dnrm2(c);
					norm[i] = x;
				}

                for (i = 0; i < System.Math.Min(M, N); i++)
				{
					/* Bring the column of largest norm into the pivot position */

					double max_norm = norm[i];
					int j;
					int kmax = i;

					for (j = i + 1; j < N; j++)
					{
						double x = norm[j];

						if (x > max_norm)
						{
							max_norm = x;
							kmax = j;
						}
					}

					if (kmax != i)
					{
						A.SwapColumns(i, kmax);
						p.Swap(i, kmax);
						norm.Swap(i, kmax);
						signum = - signum;
					}

					/* Compute the Householder transformation to reduce the j-th
					   column of the matrix to a multiple of the j-th unit vector */
				{
					PZMath_vector c_full = A.Column(i);
					PZMath_vector c = c_full.SubVector(i, M - i);

					double tau_i = HouseholderTransform(c);
					tau[i] = tau_i;

					/* Apply the transformation to the remaining columns */

					if (i + 1 < N)
					{
						PZMath_matrix m = A.Submatrix(i, i + 1, M - i, N - (i + 1));
						HouseholderHM (tau_i, c, m);
					}
				}

					/* Update the norms of the remaining columns too */
					if (i + 1 < M) 
					{
						for (j = i + 1; j < N; j++)
						{
							double x = norm[j];

							if (x > 0.0)
							{
								double y = 0;
								double temp= A[i, j] / x;

                                if (System.Math.Abs(temp) >= 1)
									y = 0.0;
								else
                                    y = x * System.Math.Sqrt(1 - temp * temp);
                      
								/* recompute norm to prevent loss of accuracy */
                                if (System.Math.Abs(y / x) < System.Math.Sqrt(20.0) * PZMath_machine.PZMath_SQRT_DBL_EPSILON)
								{
									PZMath_vector c_full = A.Column(j);
									PZMath_vector c = c_full.SubVector(i + 1, M - (i + 1));
									y = PZMath_blas.Dnrm2(c);
								}
								norm[j] = y;
							}
						}
					}
				}
			}
			return PZMath_errno.PZMath_SUCCESS;
		} // QRPTDecomp()
		#endregion

		#region Householder comments
		/* linalg/householder.c
		 * 
		 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004 Gerard Jungman, Brian Gough
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
		#region Householder methods
		public static double HouseholderTransform (PZMath_vector v)
		{
			/* replace v[0:n-1] with a householder vector (v[0:n-1]) and
			   coefficient tau that annihilate v[1:n-1] */
			int n = v.Size;
			if (n == 1)
			{
				return 0.0; /* tau = 0 */
			}
			else
			{ 
				double alpha, beta, tau ;
				PZMath_vector x = v.SubVector(1, n - 1);
				double xnorm = PZMath_blas.Dnrm2(x);
      
				if (xnorm == 0) 
				{
					return 0.0; /* tau = 0 */
				}
      
				alpha = v[0];
				beta = - (alpha >= 0.0 ? +1.0 : -1.0) * PZMath_sys.Hypot (alpha, xnorm) ;
				tau = (beta - alpha) / beta ;
      
				PZMath_blas.Dscal(1.0 / (alpha - beta), x);
				v[0] = beta;

				return tau;
			}
		} //HouseholderTransform()

		public static int HouseholderHM (double tau, PZMath_vector v, PZMath_matrix A)
		{
			/* applies a householder transformation v,tau to matrix m */
			int i, j;
			if (tau == 0.0)
			{
				return PZMath_errno.PZMath_SUCCESS;
			}
			for (j = 0; j < A.ColumnCount; j++)
			{
				/* Compute wj = Akj vk */
				double wj = A[0, j];
				for (i = 1; i < A.RowCount; i++)  /* note, computed for v(0) = 1 above */
				{
					wj += A[i, j] * v[i];
				}
				/* Aij = Aij - tau vi wj */
				/* i = 0 */
			{
				double A0j = A[0, j];
				A[0, j] = A0j - tau *  wj;
			}

				/* i = 1 .. M-1 */


				for (i = 1; i < A.RowCount; i++)
				{
					double Aij = A[i, j];
					double vi = v[i];
					A[i, j] = Aij - tau * vi * wj;
				}
			}
			
			return PZMath_errno.PZMath_SUCCESS;
		} //HouseholderHM()

		public static int HouseholderMH (double tau, PZMath_vector v, PZMath_matrix A)
		{
			/* applies a householder transformation v,tau to matrix m from the
				right hand side in order to zero out rows */

			int i, j;

			if (tau == 0)
				return PZMath_errno.PZMath_SUCCESS;;

			/* A = A - tau w v' */
			for (i = 0; i < A.RowCount; i++)
			{
				double wi = A[i, 0];

				for (j = 1; j < A.ColumnCount; j++)  /* note, computed for v(0) = 1 above */
				{
					wi += A[i, j] * v[j];
				}
			      
				/* j = 0 */
			      
				{
					double Ai0 = A[i, 0];
					A[i, 0] = Ai0 - tau * wi;
				}

				/* j = 1 .. N-1 */
			      
				for (j = 1; j < A.ColumnCount; j++) 
				{
					double vj = v[j];
					double Aij = A[i, j];
					A[i, j] = Aij - tau * wi * vj;
				}
			}

			return PZMath_errno.PZMath_SUCCESS;
		} // HouseholderMH()

		public static int HouseholderHV (double tau, PZMath_vector v, PZMath_vector w)
		{
			/* applies a householder transformation v to vector w */
			int N = v.Size;
		 
			if (tau == 0)
				return PZMath_errno.PZMath_SUCCESS;

		{
			/* compute d = v'w */
			double d0 = w[0];
			double d1, d;

			PZMath_vector v1 = v.SubVector(1, N - 1);
			PZMath_vector w1 = w.SubVector(1, N - 1);
				
			PZMath_blas.Ddot(v1, w1, out d1);
		    
			d = d0 + d1;

			/* compute w = w - tau (v) (v'w) */		  
		{
			double w0 = w[0];
			w[0] = w0 - tau * d;
		}
			PZMath_blas.Daxpy(-tau * d, v1, w1);
		}
			return PZMath_errno.PZMath_SUCCESS;
		} // HouseholderHV()

		public static int HouseholderHM1 (double tau, PZMath_matrix A)
		{
			/* applies a householder transformation v,tau to a matrix being
			   build up from the identity matrix, using the first column of A as
			   a householder vector */

			int i, j;

			if (tau == 0)
			{
				A[0, 0] = 1.0;
      
				for (j = 1; j < A.ColumnCount; j++)
				{
					A[0, j] = 0.0;
				}

				for (i = 1; i < A.RowCount; i++)
				{
					A[i, 0] = 0.0;
				}

				return PZMath_errno.PZMath_SUCCESS;
			}

			/* w = A' v */
			for (j = 1; j < A.ColumnCount; j++)
			{
				double wj = 0.0;   /* A0j * v0 */

				for (i = 1; i < A.RowCount; i++)
				{
					double vi = A[i, 0];
					wj += A[i, j] * vi;
				}

				/* A = A - tau v w' */

				A[0, j] = - tau *  wj;
      
				for (i = 1; i < A.RowCount; i++)
				{
					double vi = A[i, 0];
					double Aij = A[i, j];
					A[i, j] = Aij - tau * vi * wj;
				}
			}

			for (i = 1; i < A.RowCount; i++)
			{
				double vi = A[i, 0];
				A[i, 0] = -tau * vi;
			}

			A[0, 0] = 1.0 - tau;

			return PZMath_errno.PZMath_SUCCESS;
		} // HouseholderHM1()

		#endregion

		#region QR comments
/* linalg/qr.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman, Brian Gough
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

/* Author:  G. Jungman */

		/* Factorise a general M x N matrix A into
		 *  
		 *   A = Q R
		 *
		 * where Q is orthogonal (M x M) and R is upper triangular (M x N).
		 *
		 * Q is stored as a packed set of Householder transformations in the
		 * strict lower triangular part of the input matrix.
		 *
		 * R is stored in the diagonal and upper triangle of the input matrix.
		 *
		 * The full matrix for Q can be obtained as the product
		 *
		 *       Q = Q_k .. Q_2 Q_1
		 *
		 * where k = MIN(M,N) and
		 *
		 *       Q_i = (I - tau_i * v_i * v_i')
		 *
		 * and where v_i is a Householder vector
		 *
		 *       v_i = [1, m(i+1,i), m(i+2,i), ... , m(M,i)]
		 *
		 * This storage scheme is the same as in LAPACK.  */

		#endregion
		#region QR methods
		public static int QR_QTvec (PZMath_matrix QR, PZMath_vector tau, PZMath_vector v)
		{
			int M = QR.RowCount;
			int N = QR.ColumnCount;

            if (tau.Size != System.Math.Min(M, N))
			{	
				PZMath_errno.ERROR("size of tau must be MIN(M,N)", PZMath_errno.PZMath_EBADLEN);
			}
			else if (v.Size != M)
			{
				PZMath_errno.ERROR("vector size must be N", PZMath_errno.PZMath_EBADLEN);
			}
			else
			{
				int i;

				/* compute Q^T v */
                for (i = 0; i < System.Math.Min(M, N); i++)
				{
					PZMath_vector c = QR.Column(i);
					PZMath_vector h = c.SubVector(i, M - i);
					PZMath_vector w = v.SubVector(i, M - i);
					double ti = tau[i];
					
					HouseholderHV(ti, h, w);
				}
				return PZMath_errno.PZMath_SUCCESS;
			}
			return 0;
		} // QR_QTvec()
		#endregion

		#region SVD coments
		/* Factorise a general M x N matrix A into,
		*
		*   A = U D V^T
		*
		* where U is a column-orthogonal M x N matrix (U^T U = I), 
		* D is a diagonal N x N matrix, 
		* and V is an N x N orthogonal matrix (V^T V = V V^T = I)
		*
		* U is stored in the original matrix A, which has the same size
		*
		* V is stored as a separate matrix (not V^T). You must take the
		* transpose to form the product above.
		*
		* The diagonal matrix D is stored in the vector S,  D_ii = S_i
		*/
		#endregion
		#region SVD methods

		/* Modified algorithm which is better for M>>N */

		public static int SVDecompMod (PZMath_matrix A, PZMath_matrix X, PZMath_matrix V, PZMath_vector S, PZMath_vector work)
		{
			int i, j;

			int M = A.RowCount;
			int N = A.ColumnCount;

			if (M < N)
			{
                PZMath_errno.ERROR ("PZMath_linalg::SVDDecompMod(), svd of MxN matrix, M<N, is not implemented", PZMath_errno.PZMath_EUNIMPL);
			}
			else if (V.RowCount != N)
			{
				PZMath_errno.ERROR ("PZMath_linalg::SVDDecompMod(), square matrix V must match second dimension of matrix A",
					PZMath_errno.PZMath_EBADLEN);
			}
			else if (V.RowCount != V.ColumnCount)
			{
				PZMath_errno.ERROR ("PZMath_linalg::SVDDecompMod(), matrix V must be square", PZMath_errno.PZMath_ENOTSQR);
			}
			else if (X.RowCount != N)
			{
				PZMath_errno.ERROR ("PZMath_linalg::SVDDecompMod(), square matrix X must match second dimension of matrix A",
					PZMath_errno.PZMath_EBADLEN);
			}
			else if (X.RowCount != X.ColumnCount)
			{
				PZMath_errno.ERROR ("PZMath_linalg::SVDDecompMod(), matrix X must be square", PZMath_errno.PZMath_ENOTSQR);
			}
			else if (S.Size != N)
			{
				PZMath_errno.ERROR ("PZMath_linalg::SVDDecompMod(), length of vector S must match second dimension of matrix A",
					PZMath_errno.PZMath_EBADLEN);
			}
			else if (work.Size != N)
			{
				PZMath_errno.ERROR ("PZMath_linalg::SVDDecompMod(), length of workspace must match second dimension of matrix A",
					PZMath_errno.PZMath_EBADLEN);
			}

			if (N == 1)
			{
				PZMath_vector column = A.Column(0);
				double norm = PZMath_blas.Dnrm2(column);

				S[0] = norm;
				V[0, 0] = 1.0;
      
				if (norm != 0.0)
				{
					PZMath_blas.Dscal (1.0/norm, column);
				}

				return PZMath_errno.PZMath_SUCCESS;
			}

			/* Convert A into an upper triangular matrix R */
            
            

			for (i = 0; i < N; i++)
			{
				PZMath_vector c = A.Column(i);
				PZMath_vector v = c.SubVector(i, M - i);
				double tau_i = PZMath_linalg.HouseholderTransform(v);

				/* Apply the transformation to the remaining columns */

				if (i + 1 < N)
				{
					PZMath_matrix m =A.Submatrix(i, i + 1, M - i, N - (i + 1));
					PZMath_linalg.HouseholderHM(tau_i, v, m);
				}

				S[i] = tau_i;
			}

            

			/* Copy the upper triangular part of A into X */

			for (i = 0; i < N; i++)
			{
				for (j = 0; j < i; j++)
				{
					X[i, j] = 0.0;
				}

			{
				double Aii = A[i, i];
				X[i, i] = Aii;
			}

				for (j = i + 1; j < N; j++)
				{
					double Aij = A[i, j];
					X[i, j] = Aij;
				}
			}

           
			/* Convert A into an orthogonal matrix L */

			for (j = N; j > 0 && (j -- > 0);)
			{
				/* Householder column transformation to accumulate L */
				double tj = S[j];
				PZMath_matrix m = A.Submatrix(j, j, M - j, N - j);
				PZMath_linalg.HouseholderHM1(tj, m);
			}           

			/* unpack R into X V S */

			PZMath_linalg.SVDecomp(X, V, S, work);            

			/* Multiply L by X, to obtain U = L X, stored in U */

		{
			PZMath_vector sum = work.SubVector(0, N);

			for (i = 0; i < M; i++)
			{
				PZMath_vector L_i = A.Row(i);
				sum.SetZero();

				for (j = 0; j < N; j++)
				{
					double Lij = L_i[j];
					PZMath_vector X_j = X.Row(j);                   
					PZMath_blas.Daxpy(Lij, X_j, sum);                   

				}

				L_i.MemCopyFrom(sum);
			}
		}
			return PZMath_errno.PZMath_SUCCESS;
		} // SVDecompMod()

		public static int SVDecomp (PZMath_matrix A, PZMath_matrix V, PZMath_vector S, PZMath_vector work)
		{
			int a, b, i, j;

			int M = A.RowCount;
			int N = A.ColumnCount;
            int K = System.Math.Min(M, N);

			if (M < N)
			{
				PZMath_errno.ERROR ("PZMath_linalg::SVDDecomp(),svd of MxN matrix, M<N, is not implemented", PZMath_errno.PZMath_EUNIMPL);
			}
			else if (V.RowCount != N)
			{
				PZMath_errno.ERROR ("PZMath_linalg::SVDDecomp(), square matrix V must match second dimension of matrix A",
					PZMath_errno.PZMath_EBADLEN);
			}
			else if (V.RowCount != V.ColumnCount)
			{
				PZMath_errno.ERROR ("PZMath_linalg::SVDDecomp(),matrix V must be square", PZMath_errno.PZMath_ENOTSQR);
			}
			else if (S.Size != N)
			{
				PZMath_errno.ERROR ("PZMath_linalg::SVDDecomp(),length of vector S must match second dimension of matrix A",
					PZMath_errno.PZMath_EBADLEN);
			}
			else if (work.Size != N)
			{
				PZMath_errno.ERROR ("PZMath_linalg::SVDDecomp(),length of workspace must match second dimension of matrix A",
					PZMath_errno.PZMath_EBADLEN);
			}

			/* Handle the case of N = 1 (SVD of a column vector) */

			if (N == 1)
			{
				PZMath_vector column = A.Column(0);
				double norm = PZMath_blas.Dnrm2(column);

				S[0] = norm;
				V[0, 0] = 1.0;
      
				if (norm != 0.0)
				{
					PZMath_blas.Dscal(1.0 / norm, column);
				}

				return PZMath_errno.PZMath_SUCCESS;
			}
  
		{
			PZMath_vector f = work.SubVector(0, K - 1);
    
			/* bidiagonalize matrix A, unpack A into U S V */
			PZMath_linalg.BidiagDecomp (A, S, f);
			PZMath_linalg.BidiagUnpack2 (A, S, f, V);
    
			/* apply reduction steps to B=(S,Sd) */
			ChopSmallElements (S, f);
    
			/* Progressively reduce the matrix until it is diagonal */
    
			b = N - 1;
    
			while (b > 0)
			{
				double fbm1 = f[b - 1];

				if (fbm1 == 0.0 || PZMath_sys.Isnan (fbm1))
				{
					b--;
					continue;
				}
        
				/* Find the largest unreduced block (a,b) starting from b
				   and working backwards */
        
				a = b - 1;
        
				while (a > 0)
				{
					double fam1 = f[a - 1];

					if (fam1 == 0.0 || PZMath_sys.Isnan (fam1))
					{
						break;
					}
            
					a--;
				}
        
			{
				int n_block = b - a + 1;
				PZMath_vector S_block = S.SubVector(a, n_block);
				PZMath_vector f_block = f.SubVector (a, n_block - 1);
          
				PZMath_matrix U_block = A.Submatrix(0, a, A.RowCount, n_block);
				PZMath_matrix V_block = V.Submatrix(0, a, V.RowCount, n_block);
				QRStep (S_block, f_block, U_block, V_block);
          
				/* remove any small off-diagonal elements */
				ChopSmallElements (S_block, f_block);
			}
			}
		}
			/* Make singular values positive by reflections if necessary */
  
			for (j = 0; j < K; j++)
			{
				double Sj = S[j];
      
				if (Sj < 0.0)
				{
					for (i = 0; i < N; i++)
					{
						double Vij = V[i, j];
						V[i, j] = -Vij;
					}
          
					S[j] = -Sj;
				}
			}
  
			/* Sort singular values into decreasing order */
  
			for (i = 0; i < K; i++)
			{
				double S_max = S[i];
				int i_max = i;
      
				for (j = i + 1; j < K; j++)
				{
					double Sj = S[j];
          
					if (Sj > S_max)
					{
						S_max = Sj;
						i_max = j;
					}
				}
      
				if (i_max != i)
				{
					/* swap eigenvalues */
					S.Swap(i, i_max);
          
					/* swap eigenvectors */
					A.SwapColumns(i, i_max);
					V.SwapColumns(i, i_max);
				}
			}
  
			return PZMath_errno.PZMath_SUCCESS;
		} // SVDecomp()

		public static void ChopSmallElements (PZMath_vector d, PZMath_vector f)
		{
			int N = d.Size;
			double d_i = d[0];

			int i;

			for (i = 0; i < N - 1; i++)
			{
				double f_i = f[i];
				double d_ip1 = d[i + 1];

                if (System.Math.Abs(f_i) < PZMath_machine.PZMath_DBL_EPSILON * (System.Math.Abs(d_i) + System.Math.Abs(d_ip1)))
				{
					f[i] = 0.0;
				}
				d_i = d_ip1;
			}

		} // ChopSmallElements()

		public static void QRStep (PZMath_vector d, PZMath_vector f, PZMath_matrix U, PZMath_matrix V)
		{
			int M = U.RowCount;
			int N = V.RowCount;
			int n = d.Size;
			double y, z;
			double ak, bk, zk, ap, bp, aq, bq;
			int i, k;

			if (n == 1)
				return;  /* shouldn't happen */

			/* Compute 2x2 svd directly */

			if (n == 2)
			{		
				SVD2 (d, f, U, V);
				return;
			}

			/* Chase out any zeroes on the diagonal */

			for (i = 0; i < n - 1; i++)
			{
				double d_i = d[i];
      
				if (d_i == 0.0)
				{			
					ChaseOutIntermediateZero (d, f, U, i);
					return;
				}
			}

			/* Chase out any zero at the end of the diagonal */

		{
			double d_nm1 = d[n - 1];

			if (d_nm1 == 0.0) 
			{		
				ChaseOutTrailingZero (d, f, V);
				return;
			}
		}


			/* Apply QR reduction steps to the diagonal and offdiagonal */

		{
			double d0 = d[0];
			double f0 = f[0];
    
			double d1 = d[1];
			double f1 = f[1];
    
		{	
			double mu = TrailingEigenvalue (d, f);
    
			y = d0 * d0 - mu;
			z = d0 * f0;
		}
    
			/* Set up the recurrence for Givens rotations on a bidiagonal matrix */
    
			ak = 0;
			bk = 0;
    
			ap = d0;
			bp = f0;
    
			aq = d1;
			bq = f1;
		}

			for (k = 0; k < n - 1; k++)
			{
				double c, s;
				CreateGivens (y, z, out c, out s);

				/* Compute V <= V G */

				for (i = 0; i < N; i++)
				{
					double Vip = V[i, k];
					double Viq = V[i, k + 1];
					V[i, k] = c * Vip - s * Viq;
					V[i, k + 1] = s * Vip + c * Viq;
				}

				/* compute B <= B G */

			{
				double bk1 = c * bk - s * z;

				double ap1 = c * ap - s * bp;
				double bp1 = s * ap + c * bp;
				double zp1 = -s * aq;

				double aq1 = c * aq;

				if (k > 0)
				{
					f[k - 1] = bk1;
				}

				ak = ap1;
				bk = bp1;
				zk = zp1;

				ap = aq1;

				if (k < n - 2)
				{
					bp = f[k + 1];
				}
				else
				{
					bp = 0.0;
				}

				y = ak;
				z = zk;
			}
				CreateGivens (y, z, out c, out s);

				/* Compute U <= U G */

				for (i = 0; i < M; i++)
				{
					double Uip = U[i, k];
					double Uiq = U[i, k + 1];
					U[i, k] = c * Uip - s * Uiq;
					U[i, k + 1] = s * Uip + c * Uiq;
				}

				/* compute B <= G^T B */

			{
				double ak1 = c * ak - s * zk;
				double bk1 = c * bk - s * ap;
				double zk1 = -s * bp;

				double ap1 = s * bk + c * ap;
				double bp1 = c * bp;

				d[k] = ak1;

				ak = ak1;
				bk = bk1;
				zk = zk1;

				ap = ap1;
				bp = bp1;

				if (k < n - 2)
				{
					aq = d[k + 2];
				}
				else
				{
					aq = 0.0;
				}

				y = bk;
				z = zk;
			}
			}

			f[n - 2] = bk;
			d[n - 1] = ap;
		} // QRStep()
		
		public static void SVD2 (PZMath_vector d, PZMath_vector f, PZMath_matrix U, PZMath_matrix V)
		{
			int i;
			double c, s, a11, a12, a21, a22;

			int M = U.RowCount;
			int N = V.RowCount;

			double d0 = d[0];
			double f0 = f[0];
  
			double d1 = d[1];

			if (d0 == 0.0)
			{
				/* Eliminate off-diagonal element in [0,f0;0,d1] to make [d,0;0,0] */
				CreateGivens (f0, d1, out c, out s);


				/* compute B <= G^T B X,  where X = [0,1;1,0] */

				d[0] = c * f0 - s * d1;
				f[0] = s * f0 + c * d1;
				d[1] = 0.0;
      
				/* Compute U <= U G */

				for (i = 0; i < M; i++)
				{
					double Uip = U[i, 0];
					double Uiq = U[i, 1];
					U[i, 0] = c * Uip - s * Uiq;
					U[i, 1] = s * Uip + c * Uiq;
				}

				/* Compute V <= V X */

				V.SwapColumns(0, 1);

				return;
			}
			else if (d1 == 0.0)
			{
				/* Eliminate off-diagonal element in [d0,f0;0,0] */
				CreateGivens (d0, f0, out c, out s);
				/* compute B <= B G */

				d[0] = d0 * c - f0 * s;
				f[0] = 0.0;

				/* Compute V <= V G */

				for (i = 0; i < N; i++)
				{
					double Vip = V[i, 0];
					double Viq = V[i, 1];
					V[i, 0] = c * Vip - s * Viq;
					V[i, 1] = s * Vip + c * Viq;
				}

				return;
			}
			else
			{
				/* Make columns orthogonal, A = [d0, f0; 0, d1] * G */
				CreateSchur (d0, f0, d1, out c, out s);
      
				/* compute B <= B G */
      
				a11 = c * d0 - s * f0;
				a21 = - s * d1;
      
				a12 = s * d0 + c * f0;
				a22 = c * d1;
      
				/* Compute V <= V G */
      
				for (i = 0; i < N; i++)
				{
					double Vip = V[i, 0];
					double Viq = V[i, 1];
					V[i, 0] = c * Vip - s * Viq;
					V[i, 1] = s * Vip + c * Viq;
				}
      
				/* Eliminate off-diagonal elements, bring column with largest
				   norm to first column */
      
				if (PZMath_sys.Hypot (a11, a21) < PZMath_sys.Hypot (a12,a22))
				{
					double t1, t2;

					/* B <= B X */

					t1 = a11; a11 = a12; a12 = t1;
					t2 = a21; a21 = a22; a22 = t2;

					/* V <= V X */
					V.SwapColumns(0, 1);
				} 
				CreateGivens (a11, a21, out c, out s);
      
				/* compute B <= G^T B */
				
				d[0] = c * a11 - s * a21;
				f[0] = c * a12 - s * a22;
				d[1] = s * a12 + c * a22;
      
				/* Compute U <= U G */
      
				for (i = 0; i < M; i++)
				{
					double Uip = U[i, 0];
					double Uiq = U[i, 1];
					U[i, 0] = c * Uip - s * Uiq;
					U[i, 1] = s * Uip + c * Uiq;
				}

				return;
			}
		} // SVD2()

		public static void ChaseOutIntermediateZero (PZMath_vector d, PZMath_vector f, PZMath_matrix U, int k0)
		{
			int M = U.RowCount;
			int n = d.Size;
			double c, s;
			double x, y;
			int i, k;

			x = f[k0];
			y = d[k0 + 1];

			for (k = k0; k < n - 1; k++)
			{
						
				CreateGivens (y, -x, out c, out s);
      
				/* Compute U <= U G */
      
				for (i = 0; i < M; i++)
				{
					double Uip = U[i, k0];
					double Uiq = U[i, k + 1];
					U[i, k0] = c * Uip - s * Uiq;
					U[i, k + 1] = s * Uip + c * Uiq;
				}
      
				/* compute B <= G^T B */
      
				d[k + 1] = s * x + c * y;

				if (k == k0)
					f[k] = c * x - s * y;

				if (k < n - 2) 
				{
					double z = f[k + 1];
					f[k + 1] = c * z; 

					x = -s * z ;
					y = d[k + 2];
				}
			}
		} // ChaseOutIntermediateZero()

		public static void ChaseOutTrailingZero (PZMath_vector d, PZMath_vector f, PZMath_matrix V)
		{
			int N = V.RowCount;
			int n = d.Size;
			double c, s;
			double x, y;
			int i, k;

			x = d[n - 2];
			y = f[n - 2];

			for (k = n - 1; k > 0 && (k -- > 0);)
			{	
				CreateGivens (x, y, out c, out s);

				/* Compute V <= V G where G = [c, s ; -s, c] */
      
				for (i = 0; i < N; i++)
				{
					double Vip = V[i, k];
					double Viq = V[i, n - 1];
					V[i, k] = c * Vip - s * Viq;
					V[i, n - 1] = s * Vip + c * Viq;
				}

				/* compute B <= B G */
      
				d[k] = c * x - s * y;

				if (k == n - 2)
					f[k] = s * x + c * y;

				if (k > 0) 
				{
					double z = f[k - 1];
					f[k - 1] = c * z; 

					x = d[k - 1];
					y = s * z ;
				}
			}
		} // ChaseOutTrailingZero()

		public static double TrailingEigenvalue (PZMath_vector d, PZMath_vector f)
		{
			int n = d.Size;

			double da = d[n - 2];
			double db = d[n - 1];
			double fa = (n > 2) ? f[n - 3] : 0.0;
			double fb = f[n - 2];

			double ta = da * da + fa * fa;
			double tb = db * db + fb * fb;
			double tab = da * fb;

			double dt = (ta - tb) / 2.0;

			double mu;

			if (dt >= 0)
				{
				mu = tb - (tab * tab) / (dt + PZMath_sys.Hypot (dt, tab));
				}
			else 
				{
				mu = tb + (tab * tab) / ((-dt) + PZMath_sys.Hypot (dt, tab));
				}
			return mu;
		} // TrailingEigenvalue()

		/* Generate a Givens rotation (cos,sin) which takes v=(x,y) to (|v|,0) 
		 * From Golub and Van Loan, "Matrix Computations", Section 5.1.8 */
		public static void CreateGivens (double a, double b, out double c, out double s)
		{
			if (b == 0)
			{
				c = 1;
				s = 0;
			}
            else if (System.Math.Abs(b) > System.Math.Abs(a))
			{
				double t = -a / b;
                double s1 = 1.0 / System.Math.Sqrt(1 + t * t);
				s = s1;
				c = s1 * t;
			}
			else
			{
				double t = -b / a;
                double c1 = 1.0 / System.Math.Sqrt(1 + t * t);
				c = c1;
				s = c1 * t;
			}
		} // CreateGivens()

		public static void CreateSchur (double d0, double f0, double d1, out double c, out double s)
		{
			double apq = 2.0 * d0 * f0;
  
			if (apq != 0.0)
			{
				double t;
				double tau = (f0*f0 + (d1 + d0)*(d1 - d0)) / apq;
      
				if (tau >= 0.0)
				{
					t = 1.0/(tau + PZMath_sys.Hypot (1.0, tau));
				}
				else
				{
					t = -1.0/(-tau + PZMath_sys.Hypot (1.0, tau));
				}

				c = 1.0 / PZMath_sys.Hypot (1.0, t);
				s = t * c;
			}
			else
			{
				c = 1.0;
				s = 0.0;
			}
		} // CreateSchur()


		#endregion

		#region bidiag comments
		/* Factorise a matrix A into
		*
		* A = U B V^T
		*
		* where U and V are orthogonal and B is upper bidiagonal. 
		*
		* On exit, B is stored in the diagonal and first superdiagonal of A.
		*
		* U is stored as a packed set of Householder transformations in the
		* lower triangular part of the input matrix below the diagonal.
		*
		* V is stored as a packed set of Householder transformations in the
		* upper triangular part of the input matrix above the first
		* superdiagonal.
		*
		* The full matrix for U can be obtained as the product
		*
		*       U = U_1 U_2 .. U_N
		*
		* where 
		*
		*       U_i = (I - tau_i * u_i * u_i')
		*
		* and where u_i is a Householder vector
		*
		*       u_i = [0, .. , 0, 1, A(i+1,i), A(i+3,i), .. , A(M,i)]
		*
		* The full matrix for V can be obtained as the product
		*
		*       V = V_1 V_2 .. V_(N-2)
		*
		* where 
		*
		*       V_i = (I - tau_i * v_i * v_i')
		*
		* and where v_i is a Householder vector
		*
		*       v_i = [0, .. , 0, 1, A(i,i+2), A(i,i+3), .. , A(i,N)]
		*
		* See Golub & Van Loan, "Matrix Computations" (3rd ed), Algorithm 5.4.2 
		*
		* Note: this description uses 1-based indices. The code below uses
		* 0-based indices 
		*/
		#endregion
		#region bidiag methods
		public static int BidiagDecomp (PZMath_matrix A, PZMath_vector tau_U, PZMath_vector tau_V)  
		{
			if (A.RowCount < A.ColumnCount)
				PZMath_errno.ERROR ("PZMath_linalg::BidiagDecomp(), bidiagonal decomposition requires M>=N", PZMath_errno.PZMath_EBADLEN);
			else if (tau_U.Size != A.ColumnCount)
				PZMath_errno.ERROR ("PZMath_linalg::BidiagDecomp(),size of tau_U must be N", PZMath_errno.PZMath_EBADLEN);
			else if (tau_V.Size + 1 != A.ColumnCount)
				PZMath_errno.ERROR ("PZMath_linalg::BidiagDecomp(),size of tau_V must be (N - 1)", PZMath_errno.PZMath_EBADLEN);
			else
			{
				int M = A.RowCount;
				int N = A.ColumnCount;
				int i;
  
				for (i = 0 ; i < N; i++)
				{
					/* Apply Householder transformation to current column */
          
				{
					PZMath_vector c = A.Column(i);
					PZMath_vector v = c.SubVector(i, M - i);
					double tau_i = PZMath_linalg.HouseholderTransform(v);
            
					/* Apply the transformation to the remaining columns */
            
					if (i + 1 < N)
					{
						PZMath_matrix m = A.Submatrix(i, i + 1, M - i, N - (i + 1));
						PZMath_linalg.HouseholderHM(tau_i, v, m);
					}

					tau_U[i] = tau_i;            

				}

					/* Apply Householder transformation to current row */
          
					if (i + 1 < N)
					{
						PZMath_vector r = A.Row(i);
						PZMath_vector v = r.SubVector(i + 1, N - (i + 1));
						double tau_i = PZMath_linalg.HouseholderTransform (v);
              
						/* Apply the transformation to the remaining rows */
              
						if (i + 1 < M)
						{
							PZMath_matrix m = A.Submatrix(i+1, i+1, M - (i+1), N - (i+1));
							PZMath_linalg.HouseholderMH(tau_i, v, m);
						}

						tau_V[i] = tau_i;
					}
				}
			}
        
			return PZMath_errno.PZMath_SUCCESS;
		} // BidiagDecomp()
		
		public static int BidiagUnpack2 (PZMath_matrix A, PZMath_vector tau_U, PZMath_vector tau_V, PZMath_matrix V)
		{
			int M = A.RowCount;
			int N = A.ColumnCount;

            int K = System.Math.Min(M, N);

			if (M < N)
			{
				PZMath_errno.ERROR ("PZMath_linalg::BidiagUnpack2(),matrix A must have M >= N", PZMath_errno.PZMath_EBADLEN);
			}
			else if (tau_U.Size != K)
			{
				PZMath_errno.ERROR ("PZMath_linalg::BidiagUnpack2(),size of tau must be MIN(M,N)", PZMath_errno.PZMath_EBADLEN);
			}
			else if (tau_V.Size + 1 != K)
			{
				PZMath_errno.ERROR ("PZMath_linalg::BidiagUnpack2(),size of tau must be MIN(M,N) - 1", PZMath_errno.PZMath_EBADLEN);
			}
			else if (V.RowCount != N || V.ColumnCount != N)
			{
				PZMath_errno.ERROR ("PZMath_linalg::BidiagUnpack2(),size of V must be N x N", PZMath_errno.PZMath_EBADLEN);
			}
			else
			{
				int i, j;

				/* Initialize V to the identity */
				V.SetIdentity();

				for (i = N - 1; i > 0 && (i -- > 0);)
				{
					/* Householder row transformation to accumulate V */
					PZMath_vector r = A.Row(i);
					PZMath_vector h = r.SubVector(i + 1, N - (i+1));
          
					double ti = tau_V[i];
          
					PZMath_matrix m = V.Submatrix(i + 1, i + 1, N-(i+1), N-(i+1));
					PZMath_linalg.HouseholderHM (ti, h, m);
				}

				/* Copy superdiagonal into tau_v */

				for (i = 0; i < N - 1; i++)
				{
					double Aij = A[i, i+1];
					tau_V[i] = Aij;
				}

				/* Allow U to be unpacked into the same memory as A, copy
				   diagonal into tau_U */

				for (j = N; j > 0 && (j -- > 0);)
				{
					/* Householder column transformation to accumulate U */
					double tj = tau_U[j];
					double Ajj = A[j, j];
					PZMath_matrix m = A.Submatrix (j, j, M-j, N-j);

					tau_U[j] = Ajj;
					PZMath_linalg.HouseholderHM1(tj, m);
				}
			}
			return PZMath_errno.PZMath_SUCCESS;;
		} // BidiagUnpack2()
		#endregion

        #region LU comments
        /* linalg/lu.c
         * 
         * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman, Brian Gough
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

        /* Author:  G. Jungman */

        /* Factorise a general N x N matrix A into,
         *
         *   P A = L U
         *
         * where P is a permutation matrix, L is unit lower triangular and U
         * is upper triangular.
         *
         * L is stored in the strict lower triangular part of the input
         * matrix. The diagonal elements of L are unity and are not stored.
         *
         * U is stored in the diagonal and upper triangular part of the
         * input matrix.  
         * 
         * P is stored in the permutation p. Column j of P is column k of the
         * identity matrix, where k = permutation->data[j]
         *
         * signum gives the sign of the permutation, (-1)^n, where n is the
         * number of interchanges in the permutation. 
         *
         * See Golub & Van Loan, Matrix Computations, Algorithm 3.4.1 (Gauss
         * Elimination with Partial Pivoting).
         */
        #endregion
        #region LU methods
        public static int LUDecomp(PZMath_matrix A, PZMath_permutation p, out int signum)
            // These functions factorize the square matrix A into the LU decomposition PA = LU. 
            // On output the diagonal and upper triangular part of the input matrix A contain the matrix U. 
            // The lower triangular part of the input matrix (excluding the diagonal) contains L. 
            // The diagonal elements of L are unity, and are not stored.
            
            // The permutation matrix P is encoded in the permutation p. 
            // The j-th column of the matrix P is given by the k-th column of the identity matrix, 
            // where k = p_j the j-th element of the permutation vector. 
            // The sign of the permutation is given by signum. 
            // It has the value (-1)^n, where n is the number of interchanges in the permutation.
            
            // The algorithm used in the decomposition is Gaussian Elimination with partial pivoting 
            // (Golub & Van Loan, Matrix Computations, Algorithm 3.4.1). 
        {
            if (A.ColumnCount != A.RowCount)
            {
                PZMath_errno.ERROR("PZMath_linalg::LUDecomp(), LU decomposition requires square matrix", PZMath_errno.PZMath_ENOTSQR);
                signum = 1;
                return PZMath_errno.PZMath_ENOTSQR;
            }
            else if (p.Size != A.ColumnCount)
            {
                PZMath_errno.ERROR("PZMath_linalg::LUDecomp(), permutation length must match matrix size", PZMath_errno.PZMath_EBADLEN);
                signum = 1;
                return PZMath_errno.PZMath_EBADLEN;
            }
            else
            {
                int N = A.ColumnCount;
                int i, j, k;

                signum = 1;
                p.Init();

                for (j = 0; j < N - 1; j++)
                {
                    /* Find maximum in the j-th column */
                    double ajj;
                    double max = System.Math.Abs(A[j, j]);
                    int i_pivot = j;

                    for (i = j + 1; i < N; i++)
                    {
                        double aij = System.Math.Abs(A[i, j]);

                        if (aij > max)
                        {
                            max = aij;
                            i_pivot = i;
                        }
                    }

                    if (i_pivot != j)
                    {
                        A.SwapRows(j, i_pivot);
                        p.Swap(j, i_pivot);
                        signum = -(signum);
                    }

                    ajj = A[j, j];

                    if (ajj != 0.0)
                    {
                        for (i = j + 1; i < N; i++)
                        {
                            double aij = A[i, j] / ajj;
                            A[i, j] = aij;

                            for (k = j + 1; k < N; k++)
                            {
                                double aik = A[i, k];
                                double ajk = A[j, k];
                                A[i, k] = aik - aij * ajk;
                            }
                        }
                    }
                }

                return PZMath_errno.PZMath_SUCCESS;
            }
        } // LUDecomp()

        public static int LUInvert (PZMath_matrix LU, PZMath_permutation p, PZMath_matrix inverse)
            // These functions compute the inverse of a matrix A from its LU decomposition (LU,p), 
            // storing the result in the matrix inverse. 
            // The inverse is computed by solving the system A x = b for each column of the identity matrix. 
            // It is preferable to avoid direct use of the inverse whenever possible, 
            // as the linear solver functions can obtain the same result more efficiently and reliably 
            // (consult any introductory textbook on numerical linear algebra for details).
        {
            int i;
            int n = LU.RowCount;
            int status = PZMath_errno.PZMath_SUCCESS;

            inverse.SetIdentity();

            for (i = 0; i < n; i++)
            {
                PZMath_vector c = inverse.Column(i);
                int status_i = LUSVX(LU, p, c);

                if (status_i == PZMath_errno.PZMath_SUCCESS)
                    status = status_i;
            }
            return status;
        } // LUInvert

        public static int LUSVX(PZMath_matrix LU, PZMath_permutation p, PZMath_vector x)
            // These functions solve the square system A x = b in-place using the LU decomposition of A into (LU,p). 
            // On input x should contain the right-hand side b, 
            // which is replaced by the solution on output.
        {
            if (LU.ColumnCount != LU.RowCount)
            {
                PZMath_errno.ERROR("PZMath_linalg::LUSVX(), LU matrix must be square", PZMath_errno.PZMath_ENOTSQR);
            }
            else if (LU.ColumnCount != p.Size)
            {
                PZMath_errno.ERROR("PZMath_linalg::LUSVX(), permutation length must match matrix size", PZMath_errno.PZMath_EBADLEN);
            }
            else if (LU.ColumnCount != x.Size)
            {
                PZMath_errno.ERROR("PZMath_linalg::LUSVX(), matrix size must match solution/rhs size", PZMath_errno.PZMath_EBADLEN);
            }
            else
            {
                /* Apply permutation to RHS */
                p.Permute(x);

                // todo here
                /* Solve for c using forward-substitution, L c = P b */
                PZMath_blas.Dtrsv(CBLAS_UPLO.CblasLower, CBLAS_TRANSPOSE.CblasNoTrans, CBLAS_DIAG.CblasUnit, LU, x);

                /* Perform back-substitution, U x = c */
                PZMath_blas.Dtrsv(CBLAS_UPLO.CblasUpper, CBLAS_TRANSPOSE.CblasNoTrans, CBLAS_DIAG.CblasNonUnit, LU, x);
            }
            return PZMath_errno.PZMath_SUCCESS;
        } // LUSVX()

        public static double LUDet(PZMath_matrix LU, int signum)
            // These functions compute the determinant of a matrix A from its LU decomposition, LU. 
            // The determinant is computed as the product of the diagonal elements of U 
            // and the sign of the row permutation signum.
        {
            int i;
            int n = LU.RowCount;

            double det = (double)signum;

            for (i = 0; i < n; i++)
            {
                det *= LU[i, i];
            }

            return det;
        } // LUDet()
        #endregion

        #region Cholsky comments
        /* Cholesky Decomposition
         *
         * Copyright (C) 2000  Thomas Walter
         *
         * 03 May 2000: Modified for GSL by Brian Gough
         * 29 Jul 2005: Additions by Gerard Jungman
         * Copyright (C) 2000,2001, 2002, 2003, 2005 Brian Gough, Gerard Jungman
         *
         * This is free software; you can redistribute it and/or modify it
         * under the terms of the GNU General Public License as published by the
         * Free Software Foundation; either version 2, or (at your option) any
         * later version.
         *
         * This source is distributed in the hope that it will be useful, but WITHOUT
         * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
         * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
         * for more details.
         */

        /*
         * Cholesky decomposition of a symmetrix positive definite matrix.
         * This is useful to solve the matrix arising in
         *    periodic cubic splines
         *    approximating splines
         *
         * This algorithm does:
         *   A = L * L'
         * with
         *   L  := lower left triangle matrix
         *   L' := the transposed form of L.
         *
         */
        #endregion
        #region Cholsky methods
        private static double QuietSqrt (double x)  
            /* avoids runtime error, for checking matrix for positive definiteness */
        {
            return (x >= 0) ? System.Math.Sqrt(x) : 1;
        }

        public static int CholeskyDecomp(PZMath_matrix A)
        {
            int M = A.RowCount;
            int N = A.ColumnCount;

            if (M != N)
            {
                PZMath_errno.ERROR("PZMath_linalg::CholeskyDecomp()::cholesky decomposition requires square matrix", PZMath_errno.PZMath_ENOTSQR);
            }

            int i, j, k;
            int status = 0;

            /* Do the first 2 rows explicitly.  It is simple, and faster.  And
             * one can return if the matrix has only 1 or 2 rows.  
             */

            double A_00 = A[0, 0];

            double L_00 = QuietSqrt(A_00);

            if (A_00 <= 0)
            {
                status = PZMath_errno.PZMath_EDOM;
            }

            A[0, 0] = L_00;

            if (M > 1)
            {
                double A_10 = A[1, 0];
                double A_11 = A[1, 1];

                double L_10 = A_10 / L_00;
                double diag = A_11 - L_10 * L_10;
                double L_11 = QuietSqrt(diag);

                if (diag <= 0)
                {
                    status = PZMath_errno.PZMath_EDOM;
                }

                A[1, 0] = L_10;
                A[1, 1] = L_11;
            }

            for (k = 2; k < M; k++)
            {
                double A_kk = A[k, k];

                for (i = 0; i < k; i++)
                {
                    double sum = 0;

                    double A_ki = A[k, i];
                    double A_ii = A[i, i];
                    PZMath_vector ci = A.Row(i);
                    PZMath_vector ck = A.Row(k);

                    if (i > 0)
                    {
                        PZMath_vector di = ci.SubVector(0, i);
                        PZMath_vector dk = ck.SubVector(0, i);
                        PZMath_blas.Ddot(di, dk, out sum);
                    }

                    A_ki = (A_ki - sum) / A_ii;
                    A[k, i] = A_ki;
                }

                {
                    PZMath_vector ck = A.Row(k);
                    PZMath_vector dk = ck.SubVector(0, k);
                    double sum = PZMath_blas.Dnrm2(dk);
                    double diag = A_kk - sum * sum;

                    double L_kk = QuietSqrt(diag);

                    if (diag <= 0)
                    {
                        status = PZMath_errno.PZMath_EDOM;
                    }
                    A[k, k] = L_kk;
                }
            }

            /* Now copy the transposed lower triangle to the upper triangle,
             * the diagonal is common.  
             */

            for (i = 1; i < M; i++)
            {
                for (j = 0; j < i; j++)
                {
                    double A_ij = A[i, j];
                    A[j, i] = A_ij;
                }
            }

            if (status == PZMath_errno.PZMath_EDOM)
            {
                PZMath_errno.ERROR("PZMath_linalg::CholeskyDecomp()::matrix must be positive definite", PZMath_errno.PZMath_EDOM);
            }

            return PZMath_errno.PZMath_SUCCESS;
        } // CholeskyDecomp()
        #endregion

    } // PZMath_linalg
}
