// eee.Sheffield.PZ.Math
//
// Copyright ?Ping Zou, 2007
// sg71.cherub@gmail.com

using System;
using eee.Sheffield.PZ.Math;

namespace eee.Sheffield.PZ.Math
{
	/// <summary>
	/// test PZMath multifit 
	/// </summary>
	public class TestPZMathMultifit
	{
        public static void TestNonlinearLeastSquareFit()
        {
            PZMath_multifit_fdfsolver_type T = new PZMath_multifit_fdfsolver_type();
            PZMath_multifit_fdfsolver s = new PZMath_multifit_fdfsolver();

            int N = 40;
            int status;
            int i, iter = 0;

            int n = N;
            int p = 3;

            PZMath_matrix covar = new PZMath_matrix(p, p);

            double[] y = new double[N];
            double[] sigma = new double[N];
            double h = 1e-2; // step-size
            TestNonlinearLeastSquareFitParam d = new TestNonlinearLeastSquareFitParam(n, p, y, sigma, h);
            PZMath_multifit_function_fdf f = new PZMath_multifit_function_fdf();

            double[] x_init = new double[3];
            x_init[0] = 1.0;
            x_init[1] = 0.0;
            x_init[2] = 0.0;
            PZMath_vector x = new PZMath_vector(x_init);

            f.f = new multifit_delegate_f(TestPZMathMultifit.expb_f);
            f.df = new multifit_delegate_df(TestPZMathMultifit.expb_df);
            f.fdf = new multifit_delegate_fdf(TestPZMathMultifit.expb_fdf);
            f.n = n;
            f.p = p;
            f.parms = d;

            PZMath_random_ParkMiller_Normal r = new PZMath_random_ParkMiller_Normal();
            // this is the data to be fiitted
            for (i = 0; i < n; i++)
            {
                double t = i;
                double tt = System.Math.Exp(-0.1 * t);
                d.y[i] = 1.0 + 5 * System.Math.Exp(-0.1 * t);// + r.NextVariate();
                d.sigma[i] = 0.1;
                System.Diagnostics.Debug.WriteLine("data: " + i + " " + d.y[i] + " " + d.sigma[i]);
            }

            T.Init("PZMath_multifit_fdfsolver_lmsder");
            s.Alloc(T, n, p);
            s.Init(f, x);

            s.PrintState(iter);

            do
            {
                iter++;
                status = s.Iterate();
                System.Diagnostics.Debug.WriteLine("status = " + status);

                //printf ("status = %s \n", gsl_strerror(status));

                s.PrintState(iter);

                if (status > 0)
                    break;

                status = s.TestDelta(s.dx, s.x, 1e-4, 1e-4);
            }
            while (status == PZMath_errno.PZMath_CONTINUE && iter < 500);

            s.Covar(0.0, covar);
            covar.DebugWriteLine();

            System.Console.WriteLine("A = " + s.FIT(0) + " +/- " + s.ERR(0, covar));
            System.Console.WriteLine("lambda = " + s.FIT(1) + " +/- " + s.ERR(1, covar));
            System.Console.WriteLine("b = " + s.FIT(2) + " +/- " + s.ERR(2, covar));

            double chi = s.Chi();
            System.Console.WriteLine("chisq/dof = " + (System.Math.Pow(chi, 2.0) / (double)(n - p)));
            System.Console.WriteLine("state = " + status);

            //printf ("status = %s \n", gsl_strerror(status));
        } // TestNonlinearLeastSquareFit()

        // fit data from y = exp(x), x :[0 : 2]
        public static void TestLinearLeastSquareFit()
        {
            int i, n;
            n = 19;

            // -- prepare model parameters
            double xi, yi, chisq;
            PZMath_matrix X, cov;
            PZMath_vector y, w, c;


            X = new PZMath_matrix(n, 3);
            y = new PZMath_vector(n);
            w = new PZMath_vector(n);

            c = new PZMath_vector(3);
            cov = new PZMath_matrix(3, 3);

            // -- data to be fitted
            //PZMath_random_ParkMiller_Normal r = new PZMath_random_ParkMiller_Normal();
            for (i = 0; i < n; i++)
            {
                xi = 0.1 + 0.1 * i;
                yi = System.Math.Exp(xi);	// f(x)
                System.Console.WriteLine(xi + " " + yi);

                X[i, 0] = 1.0;
                X[i, 1] = xi;
                X[i, 2] = xi * xi;
                y[i] = yi;
                w[i] = 1.0;
            }
            PZMath_multifit_linear l = new PZMath_multifit_linear();
            l.Alloc(n, 3);
            l.Wlinear(X, w, y, c, cov, out chisq);

            System.Console.WriteLine("# best fit: Y = " + c[0] + " + " + c[1] + " X + " + c[2] + " X ^ 2");
            System.Console.WriteLine("# covariance matrix:");
            cov.ScreenWriteLine();
            System.Console.WriteLine("# chisq = " + chisq);
        } // TestLinearLeastSquareFit()

        #region Nonlinear f/df
        private class TestNonlinearLeastSquareFitParam
        {
            public int n;
            public int p;
            public double[] y;
            public double[] sigma;
            public double h;			// step-size

            public TestNonlinearLeastSquareFitParam(int inputn, int inputp, double[] inputy, double[] inputsigma, double inputh)
            {
                n = inputn;
                p = inputp;
                y = inputy;
                sigma = inputsigma;
                h = inputh;
            }
        } // Param

        public static int expb_f(PZMath_vector x, object param, PZMath_vector f)
        {
            int n = (param as TestNonlinearLeastSquareFitParam).n;
            double[] y = (param as TestNonlinearLeastSquareFitParam).y;
            double[] sigma = (param as TestNonlinearLeastSquareFitParam).sigma;

            double A = x[0];
            double lambda = x[1];
            double b = x[2];

            int i;
            for (i = 0; i < n; i++)
            {
                // Model Yi = A * exp(-lambda * i) + b
                double t = (double)i;
                double Yi = A * System.Math.Exp(-lambda * t) + b;
                f[i] = (Yi - y[i]) / sigma[i];
            }
            return PZMath_errno.PZMath_SUCCESS;
        } // expb_f()

        public static int expb_df(PZMath_vector x, object param, PZMath_matrix J)
        {
            int n = (param as TestNonlinearLeastSquareFitParam).n;
            double[] sigma = (param as TestNonlinearLeastSquareFitParam).sigma;

            double A = x[0];
            double lambda = x[1];

            int i;
            for (i = 0; i < n; i++)
            {
                // Jacobian matrix J(i, j) = dfi / dxj
                // where fi = (Yi - yi) / sigma[i]
                //		 Yi = A * exp(-lambda * i) + b
                // and the xj are the parameters (A, lambda, b)

                double t = (double)i;
                double s = sigma[i];
                double e = System.Math.Exp(-lambda * t);
                J[i, 0] = e / s;
                J[i, 1] = -t * A * e / s;
                J[i, 2] = 1 / s;
            }

            return PZMath_errno.PZMath_SUCCESS;
        } // expb_df()

        public static int expb_fdf(PZMath_vector x, object param, PZMath_vector f, PZMath_matrix J)
        {
            expb_f(x, param, f);
            expb_df(x, param, J);
            return PZMath_errno.PZMath_SUCCESS;
        } // expb_fdf() 
        #endregion

	} // Test
}
