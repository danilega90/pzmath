// PZ Math
//
// -- blas  methods
//
// -- Class PZMath_blas
// 24.11.2005


using System;

namespace eee.Sheffield.PZ.Math
{
	/// <summary>
	/// linalg methods
	/// </summary>
	
    public enum CBLAS_ORDER {CblasRowMajor, CblasColMajor};
    public enum CBLAS_TRANSPOSE {CblasNoTrans, CblasTrans, CblasConjTrans};
    public enum CBLAS_UPLO {CblasUpper, CblasLower};
    public enum CBLAS_DIAG {CblasNonUnit, CblasUnit};
    public enum CBLAS_SIDE {CblasLeft, CblasRight};

	public class PZMath_blas
	{
		
		public static double Dnrm2 (PZMath_vector f)
			// return Norm 2 of vector f i.e. Sqrt(Sum(f[i] ^ 2))
		{
			double scale = 0.0;
			double ssq = 1.0;
			int i;

			if (f.Size <= 0 || f.stride <= 0) 
			{
				return 0;
			} 
			else if (f.size == 1) 
			{
                return System.Math.Abs(f[0]);
			}

			for (i = 0; i < f.size; i++) 
			{
				double x = f[i];

				if (x != 0.0) 
				{
                    double ax = System.Math.Abs(x);

					if (scale < ax) 
					{
						ssq = 1.0 + ssq * (scale / ax) * (scale / ax);
						scale = ax;
					} 
					else 
					{
						ssq += (ax / scale) * (ax / scale);
					}
				}
			}

            return scale * System.Math.Sqrt(ssq);
		} // dnrm2()

		public static void Dscal(double alpha, PZMath_vector X)
			// rescale X by the multicative factor alpha, X = alpha * X
		{
			for (int i = 0; i < X.Size; i ++)
				X[i] = alpha * X[i];
		} // dscal()

		public static int Ddot(PZMath_vector X, PZMath_vector Y, out double result)
			// scalar product of vector X and Y, i.e. (X ^ T)  Y = Sum(X[i] * Y[i]);
			// out result
		{
			result = 0.0;
			if (X.Size == Y.Size)
			{
				for (int i = 0; i < X.Size; i ++)
					result += X[i] * Y[i];
				return PZMath_errno.PZMath_SUCCESS;
			}
			else
			{
				PZMath_errno.ERROR("invalid length", PZMath_errno.PZMath_EBADLEN);
			}
			return 0;
		} // Ddot()

		public static int Daxpy (double alpha, PZMath_vector X, PZMath_vector Y)
			// vector sum, i.e. Y = aX + Y
		{
			if (X.Size == Y.Size)
			{
				for (int i = 0; i < X.Size; i ++)
					Y[i] = alpha * X[i] + Y[i];
				return PZMath_errno.PZMath_SUCCESS;
			}
			else
			{
				PZMath_errno.ERROR("invalid length", PZMath_errno.PZMath_EBADLEN);
			}
			return 0;
		} // Daxpy()

		public static int Idamax (PZMath_vector X)
			// return the index of max element of X, in absolut magtitude
		{
			double max = 0.0;
			double t;
			int index = 0;
			for (int i = 0; i < X.Size; i ++)
			{
                t = System.Math.Abs(X[i]);
				if (t > max)
				{
					max = t;
					index = i;
				}
			}
			return index;
		} // Idamax()

		public static double Dasum(PZMath_vector X)
		{
			double r = 0.0;
			int i;
			int N = X.Size;
			for (i = 0; i < N; i++)
                r += System.Math.Abs(X[i]);
			return r;
		} // Dasum()


		public static int Dgemv(int method, double alpha, PZMath_matrix A, PZMath_vector X, double beta, PZMath_vector Y)
			// These functions compute the matrix-vector product and sum 
			// y = \alpha op(A) x + \beta y, 
			// where op(A) = A, A^T, A^H for method = 
			// 0 - CblasNoTrans,
			// 1 - CblasTrans
		{
			int M = A.RowCount;
			int N = A.ColumnCount;

			if ((method == 0 && N == X.Size && M == Y.Size)
				|| (method == 1 && M == X.Size && N == Y.Size))
			{
				if (method == 0)
					// no trans
				{
					for (int i = 0; i < M; i ++)
						// each element of Y
					{
						double temp = 0;
						for (int j = 0; j < N; j ++)
							// each element of A
							temp += A[i, j] * X[j];
						Y[i] = alpha * temp + beta * Y[i];
					}
				}
				else
					// trans
				{
					for (int i = 0; i < N; i ++)
						// each element of Y
					{
						double temp = 0;
						for (int j = 0; j < M; j ++)
							// each element of A
							temp += A[j, i] * X[j];
						Y[i] = alpha * temp + beta * Y[i];
					}
				}
				return PZMath_errno.PZMath_SUCCESS;
			}
			else
			{
				PZMath_errno.ERROR ("PZMath_blas::Dgemv(), invalid length", PZMath_errno.PZMath_EBADLEN);
				return 0;
			}
		} // Dgemv()

        public static int Dtrsv(CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag, PZMath_matrix A, PZMath_vector x)
        {
            int M = A.RowCount;
            int N = A.ColumnCount;

            if (M != N)
            {
                PZMath_errno.ERROR("PZMath_blas::Dtrsv(), matrix must be square", PZMath_errno.PZMath_ENOTSQR);
            }
            else if (N != x.Size)
            {
                PZMath_errno.ERROR("PZMath_blas::Dtrsv(), invalid length", PZMath_errno.PZMath_EBADLEN);
            }


            bool nonunit = (Diag == CBLAS_DIAG.CblasNonUnit);
            int ix, jx;
            int i, j;
            CBLAS_TRANSPOSE Trans = (TransA != CBLAS_TRANSPOSE.CblasConjTrans) ? TransA : CBLAS_TRANSPOSE.CblasTrans;
            CBLAS_ORDER order = CBLAS_ORDER.CblasRowMajor;

            if (N == 0)
                return PZMath_errno.PZMath_SUCCESS;

            /* form  x := inv( A )*x */

            if ((order == CBLAS_ORDER.CblasRowMajor && Trans == CBLAS_TRANSPOSE.CblasNoTrans
                && Uplo == CBLAS_UPLO.CblasUpper)
                || (order == CBLAS_ORDER.CblasColMajor && Trans == CBLAS_TRANSPOSE.CblasTrans
                && Uplo == CBLAS_UPLO.CblasLower))
            {
                /* backsubstitution */
                ix = N - 1;
                if (nonunit)
                    x[ix] = x[ix] / A[N - 1, N - 1];
                ix--;

                for (i = N - 1; i > 0 && (i-- > 0); )
                {
                    double tmp = x[ix];
                    jx = ix + 1;
                    for (j = i + 1; j < N; j++)
                    {
                        double Aij = A[i, j];
                        tmp -= Aij * x[jx];
                        jx++;
                    }

                    if (nonunit)
                    {
                        x[ix] = tmp / A[i, i];
                    }
                    else
                    {
                        x[ix] = tmp;
                    }
                    ix--;
                }
            }
            else if ((order == CBLAS_ORDER.CblasRowMajor && Trans == CBLAS_TRANSPOSE.CblasNoTrans
                && Uplo == CBLAS_UPLO.CblasLower)
                       || (order == CBLAS_ORDER.CblasColMajor && Trans == CBLAS_TRANSPOSE.CblasTrans
                && Uplo == CBLAS_UPLO.CblasUpper))
            {
                /* forward substitution */
                ix = 0;
                if (nonunit)
                {
                    x[ix] = x[ix] / A[0, 0];
                }
                ix++;
                for (i = 1; i < N; i++)
                {
                    double tmp = x[ix];
                    jx = 0;
                    for (j = 0; j < i; j++)
                    {
                        double Aij = A[i, j];
                        tmp -= Aij * x[jx];
                        jx++;
                    }
                    if (nonunit)
                    {
                        x[ix] = tmp / A[i, i];
                    }
                    else
                    {
                        x[ix] = tmp;
                    }
                    ix++;
                }
            }
            else if ((order == CBLAS_ORDER.CblasRowMajor && Trans == CBLAS_TRANSPOSE.CblasTrans
                && Uplo == CBLAS_UPLO.CblasUpper)
                       || (order == CBLAS_ORDER.CblasColMajor && Trans == CBLAS_TRANSPOSE.CblasNoTrans
                && Uplo == CBLAS_UPLO.CblasLower))
            {

                /* form  x := inv( A' )*x */

                /* forward substitution */
                ix = 0;
                if (nonunit)
                {
                    x[ix] = x[ix] / A[0, 0];
                }
                ix++;
                for (i = 1; i < N; i++)
                {
                    double tmp = x[ix];
                    jx = 0;
                    for (j = 0; j < i; j++)
                    {
                        double Aji = A[j, i];
                        tmp -= Aji * x[jx];
                        jx++;
                    }
                    if (nonunit)
                    {
                        x[ix] = tmp / A[i, i];
                    }
                    else
                    {
                        x[ix] = tmp;
                    }
                    ix++;
                }
            }
            else if ((order == CBLAS_ORDER.CblasRowMajor && Trans == CBLAS_TRANSPOSE.CblasTrans
                && Uplo == CBLAS_UPLO.CblasLower)
                       || (order == CBLAS_ORDER.CblasColMajor && Trans == CBLAS_TRANSPOSE.CblasNoTrans
                && Uplo == CBLAS_UPLO.CblasUpper))
            {

                /* backsubstitution */
                ix = N - 1;
                if (nonunit)
                {
                    x[ix] = x[ix] / A[(N - 1), (N - 1)];
                }
                ix--;
                for (i = N - 1; i > 0 && (i-- > 0); )
                {
                    double tmp = x[ix];
                    jx = ix + 1;
                    for (j = i + 1; j < N; j++)
                    {
                        double Aji = A[j, i];
                        tmp -= Aji * x[jx];
                        jx++;
                    }
                    if (nonunit)
                    {
                        x[ix] = tmp / A[i, i];
                    }
                    else
                    {
                        x[ix] = tmp;
                    }
                    ix--;
                }
            }
            else
            {
                PZMath_errno.ERROR("PZMath_blas::Dtrsv(), unrecognized operation");
            }
            return PZMath_errno.PZMath_SUCCESS;
        } // Dtrsv()
	} // PZMath_blas
}
