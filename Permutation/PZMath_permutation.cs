// PZ Math permutation
//
// -- permutation container
//
// -- Class PZMath_permutation
//
// 24.11.2005


using System;


namespace eee.Sheffield.PZ.Math
{
	public class PZMath_permutation
		// A permutation p is represented by an array of n integers in the range 0 to n-1, 
		// where each value p_i occurs once and only once
	{
		private int size = 0;

		private PZMath_vector data = null;

		#region constructors
		public PZMath_permutation() {}
		public PZMath_permutation(int s)
		{
			size = s;
			data = new PZMath_vector(s);
		} // PZMath_permutation(int)
		public PZMath_permutation(PZMath_permutation p)
		{
			size = p.Size;
			data = new PZMath_vector(size);
			for (int i = 0; i < size; i ++)
				data[i] = p[i];
		} //PZMath_permutation(PZMath_permutation permutation)
		#endregion

		#region properties
		public int Size
		{
			get {return size;}
		}
		#endregion

		#region access methods
		public double this[int index] 
		{
			get 
			{
				if (index < 0 || index >= size) 
				{
					PZMath_errno.ERROR("Index Out Of Range", PZMath_errno.PZMath_EFAILED);
				}
				return data[index];
			}
			set 
			{
				if (index < 0 || index >= size) 
				{
					PZMath_errno.ERROR("Index Out Of Range", PZMath_errno.PZMath_EFAILED);
				}
				data[index] = value;
			}
		} // [] overriding
		#endregion

		#region methods
		public void Init ()
		{
			// initialize permutation to identity
			// i.e. (0,1,2,...,n-1).
			for (int i = 0; i < size; i++)
			{
				data[i] = (double) i;
			}
		} // Init()

		public void Swap(int index1, int index2)
		{
			if (index1 < 0 || index1 >= size
				|| index2 < 0 || index2 >= size) 
			{
				PZMath_errno.ERROR("Index Out Of Range", PZMath_errno.PZMath_EFAILED);
			}
			double t = data[index1];
			data[index1] = data[index2];
			data[index2] = t;
		} // Swap()

		public PZMath_permutation Copy()
		{
			return new PZMath_permutation(this);
		} // Copy()

		public void InversePermute(PZMath_vector v)
			//An inverse permutation is a permutation in which 
			// each number and the number of the place which it occupies are exchanged. 
			// For example,
			// p1 (3 8 5 10 9 4 6 1 7 2)
			// p2 (8 10 1 6 3 7 9 2 5 4)
			// are pair of inverse permutation
			// This function applies the inverse of the permutation p to the elements of the vector v
		{
			if (v.Size != size)
				PZMath_errno.ERROR("PZMath_permutation::InversePermute(), Sice of permutation and vector are not equal", PZMath_errno.PZMath_EFAILED);
			PZMath_vector temp = new PZMath_vector(v);
			for (int i = 0; i < size; i ++)
			{
				v[(int) data[i]] = temp[i];
			}
		} // InverseVector()

        public void Permute(PZMath_vector v)
        //In-place Permutations 
        //   permute:    OUT[i]       = IN[perm[i]]     i = 0 .. N-1
        //   invpermute: OUT[perm[i]] = IN[i]           i = 0 .. N-1

        //   PERM is an index map, i.e. a vector which contains a permutation of
        //   the integers 0 .. N-1.

        //   From Knuth "Sorting and Searching", Volume 3 (3rd ed), Section 5.2
        //   Exercise 10 (answers), p 617
        {
            if (v.Size != size)
                PZMath_errno.ERROR("PZMath_permutation::Permute(), Sice of permutation and vector are not equal", PZMath_errno.PZMath_EFAILED);
            PZMath_vector temp = new PZMath_vector(v);
            for (int i = 0; i < size; i++)
            {
                v[i] = temp[(int)data[i]];
            }
        }
		#endregion

	} // PZMath_permutation
}
