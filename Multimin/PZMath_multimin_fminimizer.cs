// PZ Math multimin
/* multimin/gsl_multimin.h
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Fabrice Rossi
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

/* Modified by Tuomo Keskitalo to include fminimizer and 
   Nelder Mead related lines */

// 19.12.2005


using System;

namespace eee.Sheffield.PZ.Math
{
	public delegate double multimin_function (PZMath_vector x, object parms);

	public class PZMath_multimin_function
	{
		public multimin_function f;
		public int n;
		public object parms;

		public PZMath_multimin_function() {}
		
		public double FN_EVAL(PZMath_vector x)
		{
			return f(x, parms);
		} // FN_EVAL

	} // PZMath_multimin_function
	
	/// <summary>
	/// multimin fminimizer class
	/// </summary>
	public class PZMath_multimin_fminimizer
	{
		public PZMath_multimin_fminimizer_type type = null;
		public PZMath_multimin_function f = null;
		public double fval = 0.0;
		public PZMath_vector x = null;
		public double size = 0.0;
		public object state = null;

		public PZMath_multimin_fminimizer() {}
		
		#region methods
		public void Alloc(PZMath_multimin_fminimizer_type T, int n)
		{
			type = T;
			x = new PZMath_vector(n);
			state = new nmsimplex_state_t();
			T.alloc(state, n);
		} // Alloc()

		public int Init(PZMath_multimin_function f, PZMath_vector x, PZMath_vector step_size)
		{
			this.f = f;
			this.x.MemCopyFrom(x);
			return type.init(state, this.f, this.x, out this.size, step_size);
		} // Init()
		
		public int Iterate()
		{
			return type.iterate(state, this.f, this.x, out this.size, out this.fval);
		} // Iterate()

		public double Size()
		{
			return size;
		} // Size()
		
		public int TestSize (double size, double epsabs)
		{
			if (epsabs < 0.0)
				PZMath_errno.ERROR ("absolute tolerance is negative", PZMath_errno.PZMath_EBADTOL);
			  
			if (size < epsabs)
			{
				return PZMath_errno.PZMath_SUCCESS;
			}

			return PZMath_errno.PZMath_CONTINUE;
		} // TestSize()

		#endregion
	} // PZMath_multimin_fminimizer

	public delegate int multimin_fminimizer_alloc(object state, int n);
	public delegate int multimin_fminimizer_init(object state, PZMath_multimin_function f,
		PZMath_vector x, out double size, PZMath_vector step_size);
	public delegate int multimin_fminimizer_iterate(object state, PZMath_multimin_function f, 
		PZMath_vector x, out double size, out double fval);

	/// <summary>
	/// multimin fminimizer type 
	/// </summary>
	public class PZMath_multimin_fminimizer_type
	{
		public string name;
		public int size;
		public multimin_fminimizer_alloc alloc;
		public multimin_fminimizer_init init;
		public multimin_fminimizer_iterate iterate;

		public void Init(string inputname)
		{
			switch (inputname)
			{
				case "PZMath_multimin_fminimizer_nmsimplex" :
					// Nelder Mead Simplex
				{
					name = "nmsimplex";
					alloc = new multimin_fminimizer_alloc(nmsimplex.Alloc);
					init = new multimin_fminimizer_init(nmsimplex.Init);
					iterate = new multimin_fminimizer_iterate(nmsimplex.Iterate);
					break;
				}
			}
		} // Init()
	} // PZMath_multimin_fminimizer_type



}
