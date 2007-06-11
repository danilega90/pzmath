// PZ Math
// Based on Math.NET Open source
//
// -- error handling and error codes declaration
// GNU Copyright
/* Author:  B. Gough and G. Jungman */
// -- 
// 23.11.2005


using System;

namespace eee.Sheffield.PZ.Math
{
	/// <summary>
	/// basic f class
	/// </summary>
	public sealed class PZMath_machine
	{
		/* magic constants; mostly for the benefit of the implementation */

		/* -*-MACHINE CONSTANTS-*-
		 *
		 * PLATFORM: Whiz-O-Matic 9000
		 * FP_PLATFORM: IEEE-Virtual
		 * HOSTNAME: nnn.lanl.gov
		 * DATE: Fri Nov 20 17:53:26 MST 1998
		 */
		public const double PZMath_DBL_EPSILON =         2.2204460492503131e-16;
		public const double PZMath_SQRT_DBL_EPSILON =    1.4901161193847656e-08;
		public const double PZMath_ROOT3_DBL_EPSILON =   6.0554544523933429e-06;
		public const double PZMath_ROOT4_DBL_EPSILON =   1.2207031250000000e-04;
		public const double PZMath_ROOT5_DBL_EPSILON =   7.4009597974140505e-04;
		public const double PZMath_ROOT6_DBL_EPSILON =   2.4607833005759251e-03;
		public const double PZMath_LOG_DBL_EPSILON =    (-3.6043653389117154e+01);

		public const double PZMath_DBL_MIN =       2.2250738585072014e-308;
		public const double PZMath_SQRT_DBL_MIN =  1.4916681462400413e-154;
		public const double PZMath_ROOT3_DBL_MIN = 2.8126442852362996e-103;
		public const double PZMath_ROOT4_DBL_MIN = 1.2213386697554620e-77;
		public const double PZMath_ROOT5_DBL_MIN = 2.9476022969691763e-62;
		public const double PZMath_ROOT6_DBL_MIN = 5.3034368905798218e-52;
		public const double PZMath_LOG_DBL_MIN = (-7.0839641853226408e+02);

		public const double PZMath_DBL_MAX =       1.7976931348623157e+308;
		public const double PZMath_SQRT_DBL_MAX =  1.3407807929942596e+154;
		public const double PZMath_ROOT3_DBL_MAX =  5.6438030941222897e+102;
		public const double PZMath_ROOT4_DBL_MAX = 1.1579208923731620e+77;
		public const double PZMath_ROOT5_DBL_MAX = 4.4765466227572707e+61;
		public const double PZMath_ROOT6_DBL_MAX = 2.3756689782295612e+51;
		public const double PZMath_LOG_DBL_MAX = 7.0978271289338397e+02;

		public const double PZMath_FLT_EPSILON =       1.1920928955078125e-07;
		public const double PZMath_SQRT_FLT_EPSILON =  3.4526698300124393e-04;
		public const double PZMath_ROOT3_FLT_EPSILON = 4.9215666011518501e-03;
		public const double PZMath_ROOT4_FLT_EPSILON = 1.8581361171917516e-02;
		public const double PZMath_ROOT5_FLT_EPSILON = 4.1234622211652937e-02;
		public const double PZMath_ROOT6_FLT_EPSILON = 7.0153878019335827e-02;
		public const double PZMath_LOG_FLT_EPSILON =  (-1.5942385152878742e+01);

		public const double PZMath_FLT_MIN =      1.1754943508222875e-38;
		public const double PZMath_SQRT_FLT_MIN =  1.0842021724855044e-19;
		public const double PZMath_ROOT3_FLT_MIN = 2.2737367544323241e-13;
		public const double PZMath_ROOT4_FLT_MIN = 3.2927225399135965e-10;
		public const double PZMath_ROOT5_FLT_MIN = 2.5944428542140822e-08;
		public const double PZMath_ROOT6_FLT_MIN = 4.7683715820312542e-07;
		public const double PZMath_LOG_FLT_MIN =  (-8.7336544750553102e+01);

		public const double PZMath_FLT_MAX =       3.4028234663852886e+38;
		public const double PZMath_SQRT_FLT_MAX =  1.8446743523953730e+19;
		public const double PZMath_ROOT3_FLT_MAX = 6.9814635196223242e+12;
		public const double PZMath_ROOT4_FLT_MAX = 4.2949672319999986e+09;
		public const double PZMath_ROOT5_FLT_MAX = 5.0859007855960041e+07;
		public const double PZMath_ROOT6_FLT_MAX = 2.6422459233807749e+06;
		public const double PZMath_LOG_FLT_MAX =   8.8722839052068352e+01;

		public const double PZMath_SFLT_EPSILON =       4.8828125000000000e-04;
		public const double PZMath_SQRT_SFLT_EPSILON =  2.2097086912079612e-02;
		public const double PZMath_ROOT3_SFLT_EPSILON = 7.8745065618429588e-02;
		public const double PZMath_ROOT4_SFLT_EPSILON = 1.4865088937534013e-01;
		public const double PZMath_ROOT5_SFLT_EPSILON = 2.1763764082403100e-01;
		public const double PZMath_ROOT6_SFLT_EPSILON = 2.8061551207734325e-01;
		public const double PZMath_LOG_SFLT_EPSILON =  (-7.6246189861593985e+00);

		/* !MACHINE CONSTANTS! */


		/* a little internal backwards compatibility */
		public const double PZMath_MACH_EPS  =  PZMath_DBL_EPSILON;



		/* Here are the constants related to or derived from
		 * machine constants. These are not to be confused with
		 * the constants that define various precision levels
		 * for the precision/error system.
		 *
		 * This information is determined at configure time
		 * and is platform dependent. Edit at your own risk.
		 *
		 * PLATFORM: WHIZ-O-MATIC
		 * CONFIG-DATE: Thu Nov 19 19:27:18 MST 1998
		 * CONFIG-HOST: nnn.lanl.gov
		 */

		/* machine precision constants */
		/* #define PZMath_MACH_EPS         1.0e-15 */
		public const double PZMath_SQRT_MACH_EPS  =     3.2e-08;
		public const double PZMath_ROOT3_MACH_EPS =     1.0e-05;
		public const double PZMath_ROOT4_MACH_EPS =     0.000178;
		public const double PZMath_ROOT5_MACH_EPS =     0.00100;
		public const double PZMath_ROOT6_MACH_EPS =      0.00316;
		public const double PZMath_LOG_MACH_EPS =        (-34.54);
		
		
		// The limits for integer types are listed in the following table. 
		//These limits are also defined in the standard header file LIMITS.H.
		public const int PZMath_CHAR_BIT = 8; // Number of bits in the smallest variable that is not a bit field. 
		public const int PZMath_SCHAR_MIN = -128 ; // Minimum value for a variable of type signed char. 
		public const int PZMath_SCHAR_MAX = 127 ; // Maximum value for a variable of type signed char. 
		public const int PZMath_UCHAR_MAX = 255; // Maximum value for a variable of type unsigned char.  (0xff) 
		public const int PZMath_CHAR_MIN = -128; // Minimum value for a variable of type char. 
        public const int PZMath_CHAR_MAX = 127; //  Maximum value for a variable of type char. 
        public const int PZMath_MB_LEN_MAX = 2; //  Maximum number of bytes in a multicharacter constant.  
		public const int PZMath_SHRT_MIN = -32768; // Minimum value for a variable of type short. 
		public const int PZMath_SHRT_MAX = 32767; // Maximum value for a variable of type short. 
		public const int PZMath_USHRT_MAX = 65535; //Maximum value for a variable of type unsigned short.  (0xffff) 
		public const int PZMath_INT_MIN = -2147483647;	//Minimum value for a variable of type int. 
		public const int PZMath_INT_MAX = 2147483647; //  Maximum value for a variable of type int. 
		public const uint PZMath_UINT_MAX = 4294967295; //Maximum value for a variable of type unsigned int.  (0xffffffff) 
		public const int PZMath_LONG_MIN = -2147483647; //Minimum value for a variable of type long. 
		public const int PZMath_LONG_MAX  = 2147483647; //Maximum value for a variable of type long.  
		public const uint PZMath_ULONG_MAX = 4294967295; //Maximum value for a variable of type unsigned long.  (0xffffffff) 

	} // PZMath_machine
}
