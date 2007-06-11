// eee.Sheffield.PZ.Math
//
// Copyright ?Ping Zou, 2007
// sg71.cherub@gmail.com

using System;

namespace eee.Sheffield.PZ.Math
{
	public class TestPZMathIntegration
	{
		public static void TestIntegration()
		{

			PZMath_integration_workspace w = new PZMath_integration_workspace();
			w.Alloc(1000);

			double result = 0, error = 0;
			double expected = -4.0;
			double alpha = 1.0;

            TestIntegrationParms parms = new TestIntegrationParms(alpha);

			PZMath_f F = new PZMath_f();
			F.f = new PZMath_delegate_f(f);
			F.Parms = parms;

			PZMath_integration.Qags(F, 0, 1, 0, 1e-7, 1000, w, ref result, ref error); 

			System.Console.WriteLine ("result          = " + result);
			System.Console.WriteLine ("exact result    = " + expected);
			System.Console.WriteLine  ("estimated error = " + error);
			System.Console.WriteLine  ("actual error    = " + (result - expected));
			System.Console.WriteLine  ("intervals = " + w.size);
			System.Console.ReadLine();
			return;
        } // TestIntegration()
		
		static double f (double x, object parms) 
		{
            double alpha = (parms as TestIntegrationParms).alpha;
			double f = System.Math.Log (alpha * x) / System.Math.Sqrt (x);
			return f;
		} // f()
	}

	class TestIntegrationParms
	{
		public double alpha;
        public TestIntegrationParms(double a) 
		{
			alpha = a;
		}
	} // Parms

	
}
