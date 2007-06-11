// eee.Sheffield.PZ.Math
//
// Copyright ?Ping Zou, 2007
// sg71.cherub@gmail.com


using System;

namespace eee.Sheffield.PZ.Math
{
	public class TestPZMathMultimin
	{
		public static void TestMultiminNMSimplex()
		{
			int np = 2;
            TestMultiminNMSimplexParms parms = new TestMultiminNMSimplexParms(1.0, 2.0);
		
			PZMath_multimin_fminimizer_type T = new PZMath_multimin_fminimizer_type();
			T.Init("PZMath_multimin_fminimizer_nmsimplex");

			PZMath_multimin_fminimizer s = new PZMath_multimin_fminimizer();

			/* Initial vertex size vector */
			PZMath_vector ss = new PZMath_vector(np);
			PZMath_vector x;

			PZMath_multimin_function minex_func = new PZMath_multimin_function();

			int iter = 0, i;
			int status;
			double size;

			/* Set all step sizes to 1 */
			ss.SetAll(1.0);

			/* Starting point */
			x = new PZMath_vector(np);
			x[0] = 5.0;
			x[1] = 7.0;

			/* Initialize method and iterate */
			minex_func.f = new multimin_function(my_f);
			minex_func.n = np;
			minex_func.parms = parms;
			s.Alloc(T, np);
			s.Init(minex_func, x, ss);

			do
			{
				iter++;
				status = s.Iterate();
      
				if (status > 0) 
					break;

				size = s.Size();
				status = s.TestSize (size, 1e-2);

				if (status == 0)
				{
					System.Console.WriteLine  ("converged to minimum at");
				}

				System.Console.Write (iter + " ");

				for (i = 0; i < np; i++)
				{
					System.Console.Write (s.x[i] + " ");
				}
				System.Console.WriteLine ("f() = " + s.fval + "  size = " + size);
			}
			while (status == -2 && iter < 100);
			System.Console.ReadLine();
        } // TestMultiminNMSimplex()

		static double my_f (PZMath_vector v, object parms)
		{
			double x, y;
            double dp1 = (parms as TestMultiminNMSimplexParms).dp1;
            double dp2 = (parms as TestMultiminNMSimplexParms).dp2;
	  
			x = v[0];
			y = v[1];
	 
			return 10.0 * (x - dp1) * (x - dp1) + 20.0 * (y - dp2) * (y - dp2) + 30.0; 
		} // my_f()

	}

	class TestMultiminNMSimplexParms
	{
		public double dp1;
		public double dp2;

        public TestMultiminNMSimplexParms(double p1, double p2)
		{
			dp1 = p1;
			dp2 = p2;
		} // Parms()
	} // Parms
		

}
