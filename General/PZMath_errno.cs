// PZ Math
// Based on Math.NET Open source
//
// -- error handling and error codes declaration
// GNU Copyright
/* err/gsl_errno.h
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
* Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/
// -- 
// 23.11.2005


using System;

namespace eee.Sheffield.PZ.Math
{
	/// <summary>
	/// basic f class
	/// </summary>
	public class PZMath_errno
	{
		public const int PZMath_SUCCESS  = 0; 
		public const int PZMath_FAILURE  = -1;
		public const int PZMath_CONTINUE = -2;  /* iteration has not converged */
		public const int PZMath_EDOM     = 1;   /* input domain error, e.g sqrt(-1) */
		public const int PZMath_ERANGE   = 2;   /* output range error, e.g. exp(1e100) */
		public const int PZMath_EFAULT   = 3;   /* invalid pointer */
		public const int PZMath_EINVAL   = 4;   /* invalid argument supplied by user */
		public const int PZMath_EFAILED  = 5;   /* generic failure */
		public const int PZMath_EFACTOR  = 6;   /* factorization failed */
		public const int PZMath_ESANITY  = 7;   /* sanity check failed - shouldn't happen */
		public const int PZMath_ENOMEM   = 8;   /* malloc failed */
		public const int PZMath_EBADFUNC = 9;   /* problem with user-supplied function */
		public const int PZMath_ERUNAWAY = 10;  /* iterative process is out of control */
		public const int PZMath_EMAXITER = 11;  /* exceeded max number of iterations */
		public const int PZMath_EZERODIV = 12;  /* tried to divide by zero */
		public const int PZMath_EBADTOL  = 13;  /* user specified an invalid tolerance */
		public const int PZMath_ETOL     = 14;  /* failed to reach the specified tolerance */
		public const int PZMath_EUNDRFLW = 15;  /* underflow */
		public const int PZMath_EOVRFLW  = 16;  /* overflow  */
		public const int PZMath_ELOSS    = 17;  /* loss of accuracy */
		public const int PZMath_EROUND   = 18;  /* failed because of roundoff error */
		public const int PZMath_EBADLEN  = 19;  /* matrix, vector lengths are not conformant */
		public const int PZMath_ENOTSQR  = 20;  /* matrix not square */
		public const int PZMath_ESING    = 21;  /* apparent singularity detected */
		public const int PZMath_EDIVERGE = 22;  /* integral or series is divergent */
		public const int PZMath_EUNSUP   = 23;  /* requested feature is not supported by the hardware */
		public const int PZMath_EUNIMPL  = 24;  /* requested feature not (yet) implemented */
		public const int PZMath_ECACHE   = 25;  /* cache limit exceeded */
		public const int PZMath_ETABLE   = 26;  /* table limit exceeded */
		public const int PZMath_ENOPROG  = 27;  /* iteration is not making progress towards solution */
		public const int PZMath_ENOPROGJ = 28;  /* jacobian evaluations are not improving the solution */
		public const int PZMath_ETOLF    = 29;  /* cannot reach the specified tolerance in F */
		public const int PZMath_ETOLX    = 30;  /* cannot reach the specified tolerance in X */
		public const int PZMath_ETOLG    = 31;  /* cannot reach the specified tolerance in gradient */
		public const int PZMath_EOF      = 32;  /* end of file */

		public static void ERROR (string reason, int pzmath_errno)
		{
			string error = "Error: " + pzmath_errno + " because " + reason;
			System.Diagnostics.Debug.WriteLine(error);
			throw new ApplicationException(error);
		} // ERROR()
		
		public static void ERROR (string reason)
		{
			string error = reason;
			System.Diagnostics.Debug.WriteLine(error);
			throw new ApplicationException(error);
		} // ERROR()


	} // PZMath_errno
}
