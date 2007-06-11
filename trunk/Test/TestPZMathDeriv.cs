// eee.Sheffield.PZ.Math
//
// Copyright ?Ping Zou, 2007
// sg71.cherub@gmail.com

using System;

namespace eee.Sheffield.PZ.Math
{
	public class TestPZMathDeriv
	{		
		public static void TestDeriv()
		{
			// the way to use PZMath_ 
			
			// prepare parms
            TestDerivParms x = new TestDerivParms();
			x.x = 2;
			
			// prepare PZMath_f
			PZMath_f F = new PZMath_f();
			F.f = new PZMath_delegate_f(f);
			F.Parms = x;
			
			// use PZMath_deriv
			PZMath_deriv deriv = new PZMath_deriv();
			double result;
			double abserr;
			deriv.Central(F, 2.0d, 0.1d, out result, out abserr);

			System.Console.WriteLine("output " + result);
			System.Console.WriteLine("cos " + 2 * System.Math.Cos(2 * 2));
			System.Console.ReadLine();
        } // TestDeriv()

		public static double f (double x, object parms)
			// define your function, according the Parms defination
		{
            return System.Math.Sin((parms as TestDerivParms).x * x);
		}

        private class TestDerivParms
			// define your Parms
		{
			public int x;
		}
	}	
}
