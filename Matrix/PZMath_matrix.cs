// eee.Sheffield.PZ.Math
//
// Copyright ?Ping Zou, 2007
// sg71.cherub@gmail.com


using System;
using System.IO;
using System.Text.RegularExpressions;
using System.Collections.Generic;

namespace eee.Sheffield.PZ.Math
{
	/// <summary>
	/// matrix class
    /// 1D array stores 2D matrix
    /// M[row_index, col_index] =
    /// offset + (row_index * col) * stride + col_index * stride;
    /// by default, offset = 0; stride = 1;
	/// </summary>	
	public class PZMath_matrix
    {
        #region Fields
        private double[] data = null;	// 1D array contains data
		private int row;	// - size1 in gsl_matrix, height
		private int col;	// - size2 in gsl_matrix, width
		private int offset;	// offset of the start point of data
		private int stride;	
		private int orgrow;	// original row
		private int orgcol; // original col
		private int length;	// length of data = offset + orgrow * orgcol * stride;
        #endregion

        #region constructor
        public PZMath_matrix()
		{
			row = 0;
			col = 0;
			offset = 0;
			stride = 1;
			orgrow = 0;
			orgcol = 0;
			length = 0;
			data = null;
		} // PZMath_matrix()
		
		public PZMath_matrix(int size)
		{
			col = size;
			row = size;
			offset = 0;
			stride = 1;
			orgcol = size;
			orgrow = size;
			length = orgrow * orgcol * stride;
			data = new double[length];
		} // PZMath_matrix(int size)

		public PZMath_matrix(int r, int c)
			// r -- row
			// c -- col
		{
			row = r;
			col = c;
			offset = 0;
			stride = 1;
			orgrow = r;
			orgcol = c;
			length = orgrow * orgcol * stride;
			data = new double[length];
		} //PZMath_matrix(int row, int col)

		public PZMath_matrix(PZMath_matrix m)
			// value copy
			// reference copy is implemented by SubMatrix()
		{
			row = m.row;
			col = m.col;
			offset = m.offset;
			stride = m.stride;
			orgrow = m.orgrow;
			orgcol = m.orgcol;
			length = m.length;
			data = new double[length];
            m.data.CopyTo(data, 0);
		} // PZMath_matrix(PZMath_matrix matrix)

		public PZMath_matrix(double[,] d)
		{
			row = d.GetLength(1);
			col = d.GetLength(0);
			offset = 0;
			stride = 1;
			orgrow = row;
			orgcol = col;
			length = orgrow * orgcol;
			data = new double[length];
			for (int y = 0; y < row; y ++)
				for (int x = 0; x < col; x ++)
					data[y * col + x] = d[x, y];
		} // PZMath_matrix(double[,] input)

        /// <summary>
        /// true - black - 0
        /// false - white - 255
        /// </summary>
        /// <param name="b"></param>
        public PZMath_matrix(bool[,] b)
        {
            row = b.GetLength(1);
            col = b.GetLength(0);
            offset = 0;
            stride = 1;
            orgrow = row;
            orgcol = col;
            length = orgrow * orgcol;
            data = new double[length];
            for (int y = 0; y < row; y++)
            {
                for (int x = 0; x < col; x++)
                {
                    if (b[x, y])
                        data[y * col + x] = 0.0;
                    else
                        data[y * col + x] = 255.0;
                }
            }
        } // PZMath_matrix(bool[,])

        /// <summary>
        /// genearte a 1D matrix from a row vector
        /// </summary>
        /// <param name="v"></param>
        public PZMath_matrix(PZMath_vector v)
            : this(v, CBLAS_TRANSPOSE.CblasNoTrans)
        {
        } // PZMath_matrix (PZMath_vector)

        /// <summary>
        /// generate a 1D matrix from a column vector
        /// </summary>
        /// <param name="v"></param>
        /// <param name="transV"></param>
        public PZMath_matrix(PZMath_vector v, CBLAS_TRANSPOSE transV)
        {
            if (transV == CBLAS_TRANSPOSE.CblasTrans)
            // column vector
            {
                row = v.Size;
                col = 1;
            }
            else
            // row vector
            {
                row = 1;
                col = v.Size;
            }
            offset = 0;
            stride = 1;
            orgrow = row;
            orgcol = col;
            length = orgrow * orgcol;
            data = new double[length];
            for (int i = 0; i < v.Size; i++)
                data[i] = v[i];
        } // PZMath_matrix(PZMath_vector, CBLAS_TRANSPOSE)
		#endregion

		#region properties
		public int RowCount
		{
			get {return row;}
		}

		public int ColumnCount
		{
			get {return col;}
		}

        public double Max
        {
            get
            {
                double max = 0;
                for (int i = 0; i < length; i++)
                    if (data[i] > max)
                        max = data[i];
                return max;
            }
        }

        public double Min
        {
            get
            {
                double min = PZMath_machine.PZMath_DBL_MAX;
                for (int i = 0; i < length; i++)
                    if (data[i] < min)
                        min = data[i];
                return min;
            }
        }

        /// <summary>
        /// is diagonal matrix?
        /// </summary>
        /// <returns></returns>
        public bool IsDiagonal()
        {
            bool isDiagonal = true;
            if (col != row)
                isDiagonal = false;
            for (int i = 0; i < col; i++)
            {
                for (int j = 0; j < row; j++)
                {
                    if (i == j)
                    {
                        if (this[i, j] == 0)
                            isDiagonal = false;
                    }
                    else
                        if (this[i, j] != 0)
                            isDiagonal = false;
                }
            }

            return isDiagonal;
        }
        #endregion

        #region access methods
        public double this[int indexRow, int indexCol] 
		{
			get 
			{
				if (indexRow < 0 || indexRow >= row 
					|| indexCol < 0 || indexCol >= col) 
				{
					PZMath_errno.ERROR("Index Out Of Range", PZMath_errno.PZMath_EFAILED);
				}
				return data[offset + indexRow * orgcol * stride + indexCol * stride];
			}
			set 
			{
				if (indexRow < 0 || indexRow >= row 
					|| indexCol < 0 || indexCol >= col) 
				{
					PZMath_errno.ERROR("Index Out Of Range", PZMath_errno.PZMath_EFAILED);
				}
				data[offset + indexRow * orgcol * stride + indexCol * stride] = value;
			}
		} // [,] overriding

		public PZMath_vector Column (int indexCol)
		{
			if (indexCol < 0 || indexCol > col)
			{
				PZMath_errno.ERROR("Index Out Of Range", PZMath_errno.PZMath_EFAILED);
			}
			PZMath_vector v = new PZMath_vector();
			v.data = data;
			v.offset = offset + indexCol * stride;
			v.size = row;
			v.stride = orgcol * stride;
			v.length = length;
			return v;
		} // Column()

		public PZMath_vector Row (int indexRow)
			// These functions return a vector view of the i-th row of the matrix m. 
			// The data pointer of the new vector is set to null if i is out of range.
		{
			if (indexRow < 0 || indexRow > row)
			{
				PZMath_errno.ERROR("Index Out Of Range", PZMath_errno.PZMath_EFAILED);
			}

			PZMath_vector v = new PZMath_vector();
			v.data = data;
			v.offset = offset + indexRow * orgcol * stride;
			v.size = col;
			v.stride = stride;
			v.length = length;
			return v;
		} // Row()

		public void CopyRow2(int indexRow, PZMath_vector v)
			//This function copies the elements of the i-th row of the matrix m 
			// into the vector v. 
			// The length of the vector must be the same as the length of the row.
		{
			if (indexRow < 0 || indexRow > row)
			{
				PZMath_errno.ERROR("PZMath_matrix::GetRow(), Index Out Of Range", PZMath_errno.PZMath_EFAILED);
			}
			if (v.size != col)
				PZMath_errno.ERROR("PZMath_matrix::GetRow(), vector does not match the row");
			for (int i = 0; i < col; i ++)
				v[i] = this[indexRow, i];
		} // GetRow()



		public void SwapColumns (int indexCol1, int indexCol2)
		{
			if (indexCol1 < 0 || indexCol1 >= col 
				|| indexCol2 < 0 || indexCol2 >= col) 
			{
				PZMath_errno.ERROR("Index Out Of Range", PZMath_errno.PZMath_EFAILED);
			}

			double [] t = new double[row];
			for (int i = 0; i < row; i ++)
			{
				t[i] = this[i, indexCol1];
				this[i, indexCol1] = this[i, indexCol2];
				this[i, indexCol2] = t[i];
			}
		} // SwapColumns()

		public void SwapRows (int indexRow1, int indexRow2)
		{
			if (indexRow1 < 0 || indexRow1 >= row 
				|| indexRow2 < 0 || indexRow2 >= row) 
			{
				PZMath_errno.ERROR("Index Out Of Range", PZMath_errno.PZMath_EFAILED);
			}
			double [] t = new double[col];
			for (int j = 0; j < col; j ++)
			{
				t[j] = this[indexRow1, j];
				this[indexRow1, j] = this[indexRow2, j];
				this[indexRow2, j] = t[j];
			}
		} // SwapRows()

		public PZMath_matrix Submatrix (int offsetrow,  int offsetcol, int newrow,int newcol)
			// sub matrix, starts from [offsetrow, offsetcolumn], with row = newrowcount; col = newcolumncount
		{
			if (offsetrow < 0 || offsetrow + newrow > row
				|| offsetcol< 0 || offsetcol + newcol > col)
				PZMath_errno.ERROR("Index Out Of Range", PZMath_errno.PZMath_EFAILED);

			PZMath_matrix m = new PZMath_matrix();
			m.data = data;
			m.col = newcol;
			m.row = newrow;
			m.offset = offset + offsetrow * orgcol * stride + offsetcol * stride;
			m.stride = stride;
			m.orgcol = orgcol;
			m.orgrow = orgrow;
			m.length = length;
			return m;
		} // Submatrix()

		public void MemCopyFrom(PZMath_matrix m)
			// m is the source matrix
		{
			if (col != m.col || row != m.row)
				PZMath_errno.ERROR("PZMath_matrix::MemCopy(), Dest matrix is not equal to the source matrix.", PZMath_errno.PZMath_EFAILED);
            m.data.CopyTo(data, 0);
		} // MemCopy()
		#endregion

		#region set methods

		/// <summary>
		/// This function sets the elements of the matrix m to the corresponding elements of the identity matrix, 
		/// m(i,j) = \delta(i,j), i.e. 
		/// a unit diagonal with all off-diagonal elements zero. 
		/// This applies to both square and rectangular matrices.
		/// </summary>
		public void SetIdentity()
		{
			for (int i = 0; i < row; i ++)
				for (int j = 0; j < col; j ++)
					if (i == j)
						this[i, j] = 1.0;
		} // SetIdentity()

		public void SetRow(int r, PZMath_vector v)
		{
			if (r < 0 || r >= row)
				PZMath_errno.ERROR("PZMath_matrix::SetRow(), index out of range");
			if (v.size != col)
				PZMath_errno.ERROR("ZPMath_matrix::SetRow(), vector does not match matrix row");
			for (int i = 0; i < col; i ++)
				this[r, i] = v[i];
		} // SetRow()
		public void Setall (double x)
		{
			for (int i = 0; i < row; i ++)
				for (int j = 0; j < col; j ++)
					this[i, j] = x;
		} // Setall()
		#endregion

		#region compute methods
        /// <summary>
        /// normailize to [0, 1]
        /// </summary>
        public void Normalize()
        {
            Normalize(0, 1);
        }
        /// <summary>
        /// normalize to [lowerBound, upperBound]
        /// </summary>
        /// <param name="lowerBound"></param>
        /// <param name="upperBound"></param>
        public void Normalize(double lowerBound, double upperBound)
        {
            double min = this.Min;
            double max = this.Max;
            
            // scale
            double scale = 0.0;
            if (max != min)
                scale = (upperBound - lowerBound) / (max - min);
            else
                scale = 1.0;

            // normalize
            for (int i = 0; i < length; i++)
                data[i] = (data[i] - min) * scale;
        }

        /// <summary>
        ///Expand the matrix using mirror boundary condition
        /// 
        /// for example 
        ///
        /// A = [
        ///     1  2  3  11
        ///     4  5  6  12
        ///     7  8  9  13
        ///     ]
        ///
        /// B = BoundMirrorExpand(A) will yield
        ///
        ///     5  4  5  6  12  6
        ///     2  1  2  3  11  3
        ///     5  4  5  6  12  6 
        ///     8  7  8  9  13  9 
        ///     5  4  5  6  12  6
        ///
        /// Chenyang Xu and Jerry L. Prince, 9/9/1999
        /// http://iacl.ece.jhu.edu/projects/gvf
        /// </summary>
        /// <param name="wk">kernal width</param>
        /// <param name="hk">kernal height</param>
        /// <returns></returns>
        public PZMath_matrix BoundMirrorExpand(int hk, int wk)
        {          
            // half kernal width and height
            int halfWk = (wk - 1) / 2;
            int halfHk = (hk - 1) / 2;
            // original matrix height & width
            int height = this.RowCount;
            int width = this.ColumnCount;

            // expanded matrix height & width
            int heightA2 = height + 2 * halfHk;
            int widthA2 = width + 2 * halfWk;
            PZMath_matrix expand = new PZMath_matrix(heightA2, widthA2);

            int Wm1aHWK = width - 1 + halfWk;
            int HHKm1 = halfHk - 1;
            // upper row
            for (int xe = halfWk, x = 0; xe <= Wm1aHWK; xe++, x++)
                for (int ye = 0, y = halfHk; ye <= HHKm1; ye++, y--)
                    expand[ye, xe] = this[y, x];
            // bottom row
            int Hm1a2HHk = height - 1 + 2 * halfHk;
            for (int xe = halfWk, x = 0; xe <= Wm1aHWK; xe++, x++)
                for (int ye = height + halfHk, y = height - 2; ye <= Hm1a2HHk; ye++, y--)
                    expand[ye, xe] = this[y, x];
            // left column
            int HWKm1 = halfWk - 1;
            int Hm1aHHK = height - 1 + halfHk;
            for (int xe = 0, x = halfWk; xe <= HWKm1; xe++, x--)
                for (int ye = halfHk, y = 0; ye <= Hm1aHHK; ye++, y++)
                    expand[ye, xe] = this[y, x];
            // right column
            int Wm1a2HWK = width - 1 + 2 * halfWk;
            for (int xe = width + halfWk, x = width - 2; xe <= Wm1a2HWK; xe++, x--)
                for (int ye = halfHk, y = 0; ye <= Hm1aHHK; ye++, y++)
                    expand[ye, xe] = this[y, x];
            // center
            for (int xe = halfWk, x = 0; xe <= Wm1aHWK; xe++, x++)
                for (int ye = halfHk, y = 0; ye <= Hm1aHHK; ye++, y++)
                    expand[ye, xe] = this[y, x];
            // left-top corner
            for (int xc = 0, xe = 2 * halfWk; xc <= HWKm1; xc++, xe--)
                for (int y = 0; y <= HHKm1; y++)
                    expand[y, xc] = expand[y, xe];
            // right-top corner
            for (int xc = width + halfWk, xe = width + halfWk - 2; xc <= Wm1a2HWK; xc++, xe--)
                for (int y = 0; y <= HHKm1; y++)
                    expand[y, xc] = expand[y, xe];
            // left-bottom corner
            for (int xc = 0, xe = 2 * halfWk; xc <= HWKm1; xc++, xe--)
                for (int y = height + halfHk; y <= Hm1a2HHk; y++)
                    expand[y, xc] = expand[y, xe];
            // right-bottom corner
            for (int xc = width + halfWk, xe = width + halfWk - 2; xc <= Wm1a2HWK; xc++, xe--)
                for (int y = height + halfHk; y <= Hm1a2HHk; y++)
                    expand[y, xc] = expand[y, xe];

            return expand;
        }

        /// <summary>
        /// % Shrink the matrix to remove the padded mirror boundaries
        ///
        /// for example 
        ///
        /// A = [
        ///     5  4  5  6  12  6
        ///     2  1  2  3  11  3
        ///     5  4  5  6  12  6 
        ///     8  7  8  9  13  9 
        ///     5  4  5  6  12  6
        ///     ]
        /// 
        /// B = BoundMirrorShrink(A) will yield
        ///
        ///     1  2  3  11
        ///     4  5  6  12
        ///     7  8  9  13
        /// Chenyang Xu and Jerry L. Prince, 9/9/1999
        /// http://iacl.ece.jhu.edu/projects/gvf
        /// </summary>
        /// <returns></returns>
        public PZMath_matrix BoundMirrorShrink(int hk, int wk)
        {
            return this.Submatrix((hk - 1) / 2, (wk - 1) / 2, this.RowCount - (hk - 1), this.ColumnCount - (wk - 1));
        }

        /// <summary>
        /// % Ensure mirror boundary condition
        ///
        /// The number of rows and columns of A must be greater than 2
        ///
        /// for example (X means value that is not of interest)
        /// 
        /// A = [
        ///     X  X  X  X  X   X
        ///     X  1  2  3  11  X
        ///     X  4  5  6  12  X 
        ///     X  7  8  9  13  X 
        ///     X  X  X  X  X   X
        ///     ]
        ///
        /// B = BoundMirrorEnsure(A) will yield
        ///
        ///     5  4  5  6  12  6
        ///     2  1  2  3  11  3
        ///     5  4  5  6  12  6 
        ///     8  7  8  9  13  9 
        ///     5  4  5  6  12  6
        ///
        /// Chenyang Xu and Jerry L. Prince, 9/9/1999
        /// http://iacl.ece.jhu.edu/projects/gvf
        /// </summary>
        /// <returns></returns>
        public void BoundMirrorEnsure(int hk, int wk)
        {
            int height = this.RowCount;
            int width = this.ColumnCount;
            int halfHk = (hk - 1) / 2;
            int halfWk = (wk - 1) / 2;

            // upper and bottom row
            int Wm1mHWK = width - 1 - halfWk;
            int HHKm1 = halfHk - 1;
            for (int x = halfWk; x <= Wm1mHWK; x ++)
                for (int you = 0, yiu = 2 * halfHk, yob = height - halfHk, yib = height - halfHk - 2; you <= HHKm1;you ++, yiu --, yob ++, yib --)
                {
                    this[you, x] = this[yiu, x];
                    this[yob, x] = this[yib, x];
                }
            // left and right row
            int HWKm1 = halfWk - 1;
            for (int xol = 0, xil = 2 * halfWk, xor = width - halfWk, xir = width - halfWk - 2; xol <= HWKm1; xol ++, xil --, xor ++, xir --)
                for (int y = 0; y < height; y ++)
                {
                    this[y, xol] = this[y, xil];
                    this[y, xor] = this[y, xir];
                }
        }

        /// <summary>
        ///Expand the matrix by filling the background intensity
        /// 
        /// for example 
        ///
        /// A = [
        ///     255  0  0  0
        ///     255  255  0  255
        ///     0  255  0  255
        ///     ]
        ///
        /// B = BoundFillingBackgroundIntensityExpand(A, 255) will yield
        ///
        ///     255  255 255 255 255 255
        ///     255  255  0  0  0  0 255
        ///     255  255  255  0  255 255
        ///     255  0  255  0  255  255
        ///     255  255  255  255  255  255
        ///
        /// </summary>
        /// <param name="wk">kernal width</param>
        /// <param name="hk">kernal height</param>
        /// <returns></returns>
        public PZMath_matrix BoundFillingBackgroundIntensityExpand(int hk, int wk, double background)
        {
            // half kernal width and height
            int halfWk = (wk - 1) / 2;
            int halfHk = (hk - 1) / 2;
            // original matrix height & width
            int height = this.RowCount;
            int width = this.ColumnCount;

            // expanded matrix height & width
            int heightA2 = height + 2 * halfHk;
            int widthA2 = width + 2 * halfWk;
            PZMath_matrix expand = new PZMath_matrix(heightA2, widthA2);

            int Wm1aHWK = width - 1 + halfWk;
            int HHKm1 = halfHk - 1;
            // upper row
            for (int xe = halfWk, x = 0; xe <= Wm1aHWK; xe++, x++)
                for (int ye = 0, y = halfHk; ye <= HHKm1; ye++, y--)
                    expand[ye, xe] = background;
            // bottom row
            int Hm1a2HHk = height - 1 + 2 * halfHk;
            for (int xe = halfWk, x = 0; xe <= Wm1aHWK; xe++, x++)
                for (int ye = height + halfHk, y = height - 2; ye <= Hm1a2HHk; ye++, y--)
                    expand[ye, xe] = background;
            // left column
            int HWKm1 = halfWk - 1;
            int Hm1aHHK = height - 1 + halfHk;
            for (int xe = 0, x = halfWk; xe <= HWKm1; xe++, x--)
                for (int ye = halfHk, y = 0; ye <= Hm1aHHK; ye++, y++)
                    expand[ye, xe] = background;
            // right column
            int Wm1a2HWK = width - 1 + 2 * halfWk;
            for (int xe = width + halfWk, x = width - 2; xe <= Wm1a2HWK; xe++, x--)
                for (int ye = halfHk, y = 0; ye <= Hm1aHHK; ye++, y++)
                    expand[ye, xe] = background;
            // center
            for (int xe = halfWk, x = 0; xe <= Wm1aHWK; xe++, x++)
                for (int ye = halfHk, y = 0; ye <= Hm1aHHK; ye++, y++)
                    expand[ye, xe] = this[y, x];
            // left-top corner
            for (int xc = 0, xe = 2 * halfWk; xc <= HWKm1; xc++, xe--)
                for (int y = 0; y <= HHKm1; y++)
                    expand[y, xc] = background;
            // right-top corner
            for (int xc = width + halfWk, xe = width + halfWk - 2; xc <= Wm1a2HWK; xc++, xe--)
                for (int y = 0; y <= HHKm1; y++)
                    expand[y, xc] = background;
            // left-bottom corner
            for (int xc = 0, xe = 2 * halfWk; xc <= HWKm1; xc++, xe--)
                for (int y = height + halfHk; y <= Hm1a2HHk; y++)
                    expand[y, xc] = background;
            // right-bottom corner
            for (int xc = width + halfWk, xe = width + halfWk - 2; xc <= Wm1a2HWK; xc++, xe--)
                for (int y = height + halfHk; y <= Hm1a2HHk; y++)
                    expand[y, xc] = background;
            return expand;
        } // BoundFillingBackgroundIntensityExpand()

        /// <summary>
        /// % Shrink the matrix to remove the padded mirror boundaries
        ///
        /// for example 
        ///
        /// A = [
        ///     255  255 255 255 255 255
        ///     255  255  0  0  0  0 255
        ///     255  255  255  0  255 255
        ///     255  0  255  0  255  255
        ///     255  255  255  255  255  255
        ///     ]
        /// 
        /// B = BoundFillingBackgroundIntensityShrink(A) will yield
        ///
        ///     255  0  0  0
        ///     255  255  0  255
        ///     0  255  0  255
        /// </summary>
        /// <returns></returns>
        public PZMath_matrix BoundFillingBackgroundIntensityShrink(int hk, int wk)
        {
            return this.Submatrix((hk - 1) / 2, (wk - 1) / 2, this.RowCount - (hk - 1), this.ColumnCount - (wk - 1));
        }

        /// <summary>
        /// LU decomposition, P A = L U, return LU
        /// </summary>
        /// <param name="p"></param>
        /// <returns></returns>
        public PZMath_matrix LUDecomp(PZMath_permutation p)
        {
            PZMath_matrix LU = new PZMath_matrix(this);
            int signum;
            PZMath_linalg.LUDecomp(LU, p, out signum);
            return LU;           
        } // LUDecomp()

        /// <summary>
        /// LU Invert
        /// </summary>
        /// <returns></returns>
        public PZMath_matrix LUInvert()
        {
            PZMath_permutation p = new PZMath_permutation(row);
            PZMath_matrix inverse = new PZMath_matrix(row, col);
            PZMath_matrix LU = this.LUDecomp(p);
            PZMath_linalg.LUInvert(LU, p, inverse);
            return inverse;
        } // LUInvert()

        /// <summary>
        /// determination from LU decomp
        /// </summary>
        /// <returns></returns>
        public double LUDet()
        {
            PZMath_permutation p = new PZMath_permutation(row);
            PZMath_matrix LU = new PZMath_matrix(this);
            int signum;
            PZMath_linalg.LUDecomp(LU, p, out signum);
            return PZMath_linalg.LUDet(LU, signum);
        } // LUDet()

        /// <summary>
        /// Cholsky decomposition, A = L * L ^ T, return L (lower triangular matrix)
        /// </summary>
        /// <returns></returns>
        public PZMath_matrix CholeskyDecomp()
        {
            PZMath_matrix CD = new PZMath_matrix(this);
            PZMath_linalg.CholeskyDecomp(CD);
            // copy the lower triangle of CD to L
            PZMath_matrix L = new PZMath_matrix(row, col);
            for (int i = 0; i < row; i++)
            {
                for (int j = 0; j <= i; j++)
                {
                    L[i, j] = CD[i, j];
                }
            }
            return L;
        } // CholskyDecomp()

        #region Convolution methods
        /// <summary>
        /// x-direction convolute with an 1D double kernel
        /// </summary>
        /// <param name="kernel">1D double kernel</param>
        /// <returns></returns>
        public PZMath_matrix ConvoluteRow(PZMath_vector kernel)
        {
            return ConvoluteRowMirrorExpand(kernel);
        } // ConvoluteRow()

        /// <summary>
        /// x-direction convolute with an 1D double kernel
        /// mirror expand source matrix before convolution
        /// </summary>
        /// <param name="kernel">1D double kernel</param>
        /// <returns></returns>
        private PZMath_matrix ConvoluteRowMirrorExpand(PZMath_vector kernel)
        {
            // kernel info
            int kernelSize = kernel.Size;
            int halfKernelSize = (kernelSize - 1) / 2;

            // matrix info
            int height = row;
            int width = col;

            // expand src matrix
            PZMath_matrix expandSrc = this.BoundMirrorExpand(kernelSize, kernelSize);
            PZMath_matrix expandDst = new PZMath_matrix(expandSrc);

            // spatial convolute
            int xStart = halfKernelSize;
            int xEnd = width + halfKernelSize;
            int yStart = halfKernelSize;
            int yEnd = height + halfKernelSize;

            for (int x = xStart; x < xEnd; x++)
            {
                for (int y = yStart; y < yEnd; y++)
                {
                    double localSum = 0.0;
                    // inside kernel
                    for (int k = -1 * halfKernelSize; k <= halfKernelSize; k ++)
                    {
                        localSum += expandSrc[y, x + k] * kernel[k + halfKernelSize];                       
                    }
                    expandDst[y, x] = localSum;
                }
            }

            // shrink dst matrix
            PZMath_matrix dstMatrix = expandDst.BoundMirrorShrink(kernelSize, kernelSize);

            return dstMatrix;
        } // ConvoluteRowMirrorExpand()

        /// <summary>
        /// y-direction convolute with an 1D double kernel
        /// </summary>
        /// <param name="kernel">1D double kernel</param>
        /// <returns></returns>
        public PZMath_matrix ConvoluteColumn(PZMath_vector kernel)
        {
            return ConvoluteColumnMirrorExpand(kernel);
        } // ConvoluteColumn()

        /// <summary>
        /// y_direction convolute with an 1D double kernel
        /// mirror expand source matrix before convolution
        /// </summary>
        /// <param name="kernel">1D double kernel</param>
        /// <returns></returns>
        private PZMath_matrix ConvoluteColumnMirrorExpand(PZMath_vector kernel)
        {
            // kernel info
            int kernelSize = kernel.Size;
            int halfKernelSize = (kernelSize - 1) / 2;

            // matrix info
            int height = row;
            int width = col;

            // expand src matrix
            PZMath_matrix expandSrc = this.BoundMirrorExpand(kernelSize, kernelSize);
            PZMath_matrix expandDst = new PZMath_matrix(expandSrc);

            // spatial convolute
            int xStart = halfKernelSize;
            int xEnd = width + halfKernelSize;
            int yStart = halfKernelSize;
            int yEnd = height + halfKernelSize;

            for (int x = xStart; x < xEnd; x++)
            {
                for (int y = yStart; y < yEnd; y++)
                {
                    double localSum = 0.0;
                    // inside kernel
                    for (int k = -1 * halfKernelSize; k <= halfKernelSize; k++)
                    {
                        localSum += expandSrc[y + k, x] * kernel[k + halfKernelSize];
                    }
                    expandDst[y, x] = localSum;
                }
            }

            // shrink dst matrix
            PZMath_matrix dstMatrix = expandDst.BoundMirrorShrink(kernelSize, kernelSize);

            return dstMatrix;
        } // ConvoluteColumnMirrorExpand()

        /// <summary>
        /// kernel spatial convolution
        /// </summary>
        /// <param name="kernel"></param>
        /// <returns></returns>
        public PZMath_matrix Convolute2D(PZMath_matrix kernel)
        {
            return Convolute2DMirrorExpand(kernel);
        } // Convolute2D()
        /// <summary>
        /// kernel spatial convolution, mirror expand before convolution
        /// </summary>
        /// <param name="kernel"></param>
        /// <returns></returns>
        private PZMath_matrix Convolute2DMirrorExpand(PZMath_matrix kernel)
        {
            // kernel info
            int kernelHeight = kernel.RowCount;
            int kernelWidth = kernel.ColumnCount;
            int halfKernelHeight = (kernelHeight - 1) / 2;
            int halfKernelWidth = (kernelHeight - 1) / 2;

            // matrix info
            int height = row;
            int width = col;

            // expand src matrix
            PZMath_matrix expandSrc = this.BoundMirrorExpand(kernelHeight, kernelWidth);
            PZMath_matrix expandDst = new PZMath_matrix(expandSrc);

            // spatial convolute
            //int xStart = halfKernelWidth - 1;
            //int xEnd = width - halfKernelWidth;
            //int yStart = halfKernelHeight - 1;
            //int yEnd = height - halfKernelHeight;
            int xStart = halfKernelWidth;
            int xEnd = width + halfKernelWidth;
            int yStart = halfKernelHeight;
            int yEnd = height + halfKernelHeight;

            for (int x = xStart; x < xEnd; x++)
            {
                for (int y = yStart; y < yEnd; y++)
                {
                    double localSum = 0.0;
                    // inside kernel
                    for (int kx = -1 * halfKernelWidth; kx <= halfKernelWidth; kx++)
                    {
                        for (int ky = -1 * halfKernelHeight; ky <= halfKernelHeight; ky++)
                        {
                            //localSum += expandSrc[y, x] * kernel[ky + halfKernelHeight, kx + halfKernelWidth];
                            localSum += expandSrc[y + ky, x + kx] * kernel[ky + halfKernelHeight, kx + halfKernelWidth];
                        }
                    }
                    expandDst[y, x] = localSum;
                }
            }

            // shrink dst matrix
            PZMath_matrix dstMatrix = expandDst.BoundMirrorShrink(kernelHeight, kernelWidth);
            return dstMatrix;
        } // Convolute2DMirrorExpand()

        /// <summary>
        /// kernel is a bool template
        /// boolean operation:"AND" only the "true" elements
        /// </summary>
        /// <param name="kernel"></param>
        /// <returns></returns>
        public PZMath_matrix AndTrue2D(bool[,] kernel)
        {
            // kernel information
            int hk = kernel.GetLength(1);
            int wk = kernel.GetLength(0);
            int halfHk = (hk - 1) / 2;
            int halfWk = (wk - 1) / 2;

            // matrix info
            int height = row;
            int width = col;

            PZMath_matrix expandSrcMatrix = this.BoundFillingBackgroundIntensityExpand(hk, wk, 255);
            bool[,] expandSrcBoolMatrix = expandSrcMatrix.ConvertToBoolMatrix(125);
            
            // prepare dst bool matrix
            int expandHeight = expandSrcBoolMatrix.GetLength(1);
            int expandWidth = expandSrcBoolMatrix.GetLength(0);
            bool[,] expandDstBoolMatrix = new bool[expandWidth, expandHeight];

            int xStart = halfWk;
            int xEnd = width + halfWk;
            int yStart = halfHk;
            int yEnd = height + halfHk;

            for (int x = xStart; x < xEnd; x++)
            {
                for (int y = yStart; y < yEnd; y++)
                {
                    bool isHit = true;
                    // inside kernel
                    for (int kx = -1 * halfWk; kx <= halfWk; kx++)
                    {
                        for (int ky = -1 * halfHk; ky <= halfHk; ky++)
                        {
                            // kernel, expandSrcBoolMatrix, and expandDstBoolMatrix are bool[,]
                            if (kernel[kx + halfWk, ky + halfHk])
                                // kernel element is true
                            {
                                if (!expandSrcBoolMatrix[x + kx, y + ky])
                                    // image element is false
                                {
                                    isHit = false;
                                    break;
                                }
                            }                            
                        }
                    }
                    expandDstBoolMatrix[x, y] = isHit;
                }
            }

            PZMath_matrix expandDstMatrix = new PZMath_matrix(expandDstBoolMatrix);
            PZMath_matrix dstMatrix = expandDstMatrix.BoundFillingBackgroundIntensityShrink(hk, wk);

            return dstMatrix;
        } // AndTrue2D()       

        /// <summary>
        /// kernel is a bool template
        /// boolean AND operation for both true and false
        /// </summary>
        /// <param name="kernel"></param>
        /// <returns></returns>
        public PZMath_matrix AndTrueFalse2D(bool[,] kernel)
        {
            // kernel information
            int hk = kernel.GetLength(1);
            int wk = kernel.GetLength(0);
            int halfHk = (hk - 1) / 2;
            int halfWk = (wk - 1) / 2;

            // matrix info
            int height = row;
            int width = col;

            PZMath_matrix expandSrcMatrix = this.BoundFillingBackgroundIntensityExpand(hk, wk, 255);       
            bool[,] expandSrcBoolMatrix = expandSrcMatrix.ConvertToBoolMatrix(125);

            // prepare dst bool matrix
            int expandHeight = expandSrcBoolMatrix.GetLength(1);
            int expandWidth = expandSrcBoolMatrix.GetLength(0);
            bool[,] expandDstBoolMatrix = new bool[expandWidth, expandHeight];

            int xStart = halfWk;
            int xEnd = width + halfWk;
            int yStart = halfHk;
            int yEnd = height + halfHk;

            for (int x = xStart; x < xEnd; x++)
            {
                for (int y = yStart; y < yEnd; y++)
                {                    
                    bool isHit = true;
                    // inside kernel
                    for (int kx = -1 * halfWk; kx <= halfWk; kx++)
                    {
                        for (int ky = -1 * halfHk; ky <= halfHk; ky++)
                        {
                            // kernel is a bool[,]
                            // expandSrcBoolMatrix, and expandDstBoolMatrix are bool[,]
                            if ((kernel[kx + halfWk, ky + halfHk] && !expandSrcBoolMatrix[x + kx, y + ky])
                                || (!kernel[kx + halfWk, ky + halfHk] && expandSrcBoolMatrix[x + kx, y + ky]))
                            {
                                isHit = false;
                                break;
                            }
                        }
                    }
                    expandDstBoolMatrix[x, y] = isHit;
                }
            }

            PZMath_matrix expandDstMatrix = new PZMath_matrix(expandDstBoolMatrix);
            PZMath_matrix dstMatrix = expandDstMatrix.BoundFillingBackgroundIntensityShrink(hk, wk);

            return dstMatrix;
        } // AndTrueFalse2D()       

        #endregion
        
        #region multiresolution methods
        /// P. J. Burt, The pyramid as a structure for efficient computation , 
        ///     in A. Rosenfeld (ed), Multiresolution Image Processing and Analysis, Springer-Verlag, Berlin Heidelberg New York, 1984, pp 6 - 35
        /// W.B. Goh and K.Y. Chan, "Shape description using gradient vector field histograms", 
        ///     in Scale Space Methods in Computer Vision, vol. 2695, Lecture Notes in Computer Science, L. D. Griffin and M. Lillholm, Eds.: Springer Berlin / Heidelberg, 2003, pp. 713-728.     

        /// <summary>
        /// sub sampling, locally average by a kernel, l = 0, ..., N, N is top level, i.e. lowest resolution
        /// G_l(x, y) = sum(sum(w(m, n) * G_l-1(2x + m, 2y + n))), m, n, [-k, k]
        /// </summary>
        /// <param name="?"></param>
        /// <returns></returns>
        public PZMath_matrix Reduce(PZMath_matrix kernel)
        {
            return ReduceMirrorExpand(kernel);
        }
        /// <summary>
        /// sub sampling, locally average by a kernel, mirror expand before convolution
        /// </summary>
        /// <param name="kernel"></param>
        /// <returns></returns>
        private PZMath_matrix ReduceMirrorExpand(PZMath_matrix kernel)
        {
            // kernel info
            int kernelHeight = kernel.RowCount;
            int kernelWidth = kernel.ColumnCount;
            int halfKernelHeight = (kernelHeight - 1) / 2;
            int halfKernelWidth = (kernelHeight - 1) / 2;

            // matrix info
            int height = row;
            int width = col;
            int dstHeight = row / 2;
            int dstWidth = col / 2;

            // expand src matrix
            PZMath_matrix expandSrc = this.BoundMirrorExpand(kernelHeight, kernelWidth);
            PZMath_matrix dstMatrix = new PZMath_matrix(dstHeight, dstWidth);

            // local average (spatial convolute)
            int xOffset = halfKernelWidth;
            int yOffset = halfKernelHeight;
            for (int dstx = 0; dstx < dstWidth; dstx++)
            {
                for (int dsty = 0; dsty < dstHeight; dsty++)
                {
                    double localSum = 0.0;
                    // inside kernel
                    for (int kx = -1 * halfKernelWidth; kx <= halfKernelWidth; kx++)
                    {
                        for (int ky = -1 * halfKernelHeight; ky <= halfKernelHeight; ky++)
                        {
                            int srcx = 2 * dstx + kx + xOffset;
                            int srcy = 2 * dsty + ky + yOffset;
                            int kernelx = kx + halfKernelWidth;
                            int kernely = ky + halfKernelHeight;
                            localSum += expandSrc[srcy, srcx] * kernel[kernely, kernelx];
                        }
                    }
                    dstMatrix[dsty, dstx] = localSum;
                }
            }                       
            return dstMatrix;
        }
        
        /// <summary>
        /// interpolating by a kernel
        /// G_l(x, y) = 4 * sum(sum(w(m, n) * G_l+1((x + m) / 2, (y + n) / 2)), m, n, [-k, k]
        /// </summary>
        /// <param name="kernel"></param>
        /// <returns></returns>
        public PZMath_matrix Expand(PZMath_matrix kernel)
        {
            return ExpandMirrorExpand(kernel);
        } // Expand()
        /// <summary>
        /// interpolating by a kernel, mirror expand before convolute
        /// </summary>
        /// <param name="kernel"></param>
        /// <returns></returns>
        private PZMath_matrix ExpandMirrorExpand(PZMath_matrix kernel)
        {
            // kernel info
            int kernelHeight = kernel.RowCount;
            int kernelWidth = kernel.ColumnCount;
            int halfKernelHeight = (kernelHeight - 1) / 2;
            int halfKernelWidth = (kernelHeight - 1) / 2;

            // matrix info
            int height = row;
            int width = col;
            int dstHeight = row * 2;
            int dstWidth = col * 2;

            // expand src matrix
            PZMath_matrix expandSrc = this.BoundMirrorExpand(kernelHeight, kernelWidth);
            PZMath_matrix dstMatrix = new PZMath_matrix(dstHeight, dstWidth);
            
            // local average (spatial convolute)
            int xOffset = halfKernelWidth;
            int yOffset = halfKernelHeight;
            for (int dstx = 0; dstx < dstWidth; dstx++)
            {
                for (int dsty = 0; dsty < dstHeight; dsty++)
                {
                    double localSum = 0.0;
                    double kernelSum = 0.0;
                    // inside kernel
                    for (int kx = -1 * halfKernelWidth; kx <= halfKernelWidth; kx++)
                    {
                        for (int ky = -1 * halfKernelHeight; ky <= halfKernelHeight; ky++)
                        {
                            if (((dstx + kx) % 2 != 0) || ((dsty + ky) % 2 != 0))
                                continue;
                            int srcx = (dstx + kx) / 2 + xOffset;
                            int srcy = (dsty + ky) / 2 + yOffset;
                            int kernelx = kx + halfKernelWidth;
                            int kernely = ky + halfKernelHeight;
                            localSum += expandSrc[srcy, srcx] * kernel[kernely, kernelx];
                            kernelSum += kernel[kernely, kernelx];
                        }
                    }
                    dstMatrix[dsty, dstx] = localSum / kernelSum;
                }
            }                         
            return dstMatrix;
        } // ExpandMirrorExpand()
        #endregion

        #endregion

        #region operator override
        /// <summary>
        /// operator *
        /// </summary>
        /// <param name="leftMatrix"></param>
        /// <param name="rightMatrix"></param>
        /// <returns></returns>
        public static PZMath_matrix operator *(PZMath_matrix leftMatrix, PZMath_matrix rightMatrix)
        {
            // matrix dimension valid check
            if (leftMatrix.ColumnCount != rightMatrix.RowCount)
                PZMath_errno.ERROR("PZMath_matrix::operator *, left matrix and right matrix are match in dimension.");

            int leftColumnCount = leftMatrix.ColumnCount;
            int leftRowCount = leftMatrix.RowCount;
            int rightColumnCount = rightMatrix.ColumnCount;
            int rightRowCount = rightMatrix.RowCount;

            int newRowCount = leftRowCount;
            int newColumnCount = rightColumnCount;

            PZMath_matrix newMatrix = new PZMath_matrix(newRowCount, newColumnCount);

            // multiply operation
            int newCol, newRow;
            for (newCol = 0; newCol < newColumnCount; newCol ++)
            {
                for (newRow = 0; newRow < newRowCount; newRow ++)
                {
                    double temp = 0;
                    for (int i = 0; i < leftColumnCount; i ++)
                        temp += leftMatrix[newRow, i] * rightMatrix[i, newCol];
                    newMatrix[newRow, newCol] = temp;
                }
            }

            return newMatrix;
        } // operator *

        #endregion

        #region output methods
        public void DebugWriteLine()
			// out put matrix data as Debug.WriteLine()
		{
			if (PZMath.PZDebug)
			{
				for (int y = 0; y < row; y ++)
				{
					for (int x = 0; x < col; x ++)
					{
						System.Diagnostics.Debug.Write(String.Format(" {0, -4}", this[y, x]));
					}
					System.Diagnostics.Debug.WriteLine(" ");
				}
			}
		} // DebugWriteLine()

		public void ScreenWriteLine()
			// out put matrix data as Debug.WriteLine()
		{
            for (int y = 0; y < row; y++)
            {
                for (int x = 0; x < col; x++)
				{
					System.Console.Write(String.Format(" {0, -4}", this[y, x]));
				}
				System.Console.WriteLine(" ");
			}
		} // ScreenWriteLine()

        /// <summary>
        /// write data to file
        /// </summary>
        /// <param name="filename"></param>
        public void WriteFile(string filename)
        {
            FileStream fileStream = new FileStream(filename, FileMode.Create, FileAccess.Write);
            StreamWriter writer = new StreamWriter(fileStream);
            for (int y = 0; y < row; y++)
            {
                for (int x = 0; x < col; x++)
                {
                    writer.Write(String.Format(" {0, -20}", this[y, x]));
                }
                writer.WriteLine();
            }
            writer.Flush();
            writer.Close();
            fileStream.Close();
        }

        /// <summary>
        /// read data from file
        /// </summary>
        /// <param name="fileName"></param>
        public void ReadFile(string fileName)
        {
            // prepare stream
            FileStream fileStream = new FileStream(fileName, FileMode.Open, FileAccess.Read);
            StreamReader reader = new StreamReader(fileStream);
            string readin;
            string delimiting = " ";
            long streamPosition = reader.BaseStream.Position;

            // detect number of columns
            int totalColumn = 0;
            readin = reader.ReadLine();
            string[] values = Regex.Split(readin, delimiting);
            for (int i = 0; i < values.Length; i++)
                if (values[i] != "")
                    totalColumn++;

            // reset stream
            reader.DiscardBufferedData();
            reader.BaseStream.Seek(0, SeekOrigin.Begin);
            reader.BaseStream.Seek(streamPosition, SeekOrigin.Current);

            // detect number of rows
            int totalRow = 0;
            while (!reader.EndOfStream)
            {
                bool newline = false;
                readin = reader.ReadLine();
                values = Regex.Split(readin, delimiting);
                for (int i = 0; i < values.Length; i++)
                {
                    if (values[i] != "")
                    {
                        newline = true;
                        break;
                    }
                }
                if (newline)
                    totalRow++;
            }

            // reset stream
            reader.DiscardBufferedData();
            reader.BaseStream.Seek(0, SeekOrigin.Begin);
            reader.BaseStream.Seek(streamPosition, SeekOrigin.Current);

            // prepare data
            row = totalRow;
            col = totalColumn;
            offset = 0;
            stride = 1;
            orgrow = row;
            orgcol = col;
            length = orgrow * orgcol;
            data = new double[length];

            // read data
            int rowIndex = 0;
            while (!reader.EndOfStream)
            {
                bool newline = false;
                // parse data
                int colIndex = 0;
                readin = reader.ReadLine();
                values = Regex.Split(readin, delimiting);
                string[] effectvalues = new string[totalColumn];
                for (int i = 0; i < values.Length; i++)
                {
                    if (values[i] != "")
                    {
                        newline = true;
                        effectvalues[colIndex] = values[i];
                        colIndex++;
                    }
                }
                if (newline)
                {
                    // parse data
                    for (int i = 0; i < effectvalues.Length; i++)
                        data[rowIndex * col + i] = (double)Convert.ChangeType(effectvalues[i], typeof(double));
                }
                rowIndex++;
            }

            // free source
            reader.Close();
            fileStream.Close();
        }

        /// <summary>
        /// convert to bool matrix, 
        /// black - 1 - true - (smaller than threshold)
        /// white - 0 - false - (larger than threshold)
        /// </summary>
        /// <param name="threshold"></param>
        /// <returns></returns>
        public bool[,] ConvertToBoolMatrix(int threshold)
        {
            bool[,] dstBoolMatrix = new bool[col, row];
            for (int y = 0; y < row; y++)
            {
                for (int x = 0; x < col; x++)
                {
                    if (this[y, x] <= threshold)
                        dstBoolMatrix[x, y] = true;
                    else
                        dstBoolMatrix[x, y] = false;
                }
            }

            return dstBoolMatrix;
        } // ConvertToBoolMatrix()
        #endregion

        #region example codes
        /// <summary>
        /// LU decomposition and LU inverse example
        /// </summary>
        public static void LUExample()
        {
            double[,] m = {{2, 3, 3, 5},
                           {6, 6, 8, 9}, 
                           {10, 11, 12, 13},
                           {14, 15, 17, 17}};
            PZMath_matrix matrix = new PZMath_matrix(m);
            matrix.ScreenWriteLine();

            // LU decomposition
            // correct answer:
            // L =
            // 1.0000         0         0         0
            // 0.1429    1.0000         0         0
            // 0.4286   -0.5000    1.0000         0
            // 0.7143    0.3333   -0.3333    1.0000

            // U =
            // 14.0000   15.0000   17.0000   17.0000
            //       0    0.8571    0.5714    2.5714
            //       0         0    1.0000    3.0000
            //       0         0         0    1.0000

            // P =
            // 0     0     0     1
            // 1     0     0     0
            // 0     1     0     0
            // 0     0     1     0
            PZMath_permutation p = new PZMath_permutation(matrix.RowCount);
            PZMath_matrix LU = matrix.LUDecomp(p);
            LU.ScreenWriteLine();
            
            // LU Invert
            // correct answer:
            // inv(m) = 
            // -2.0833    0.6667    3.5000   -2.4167
            // 1.0000   -1.0000   -1.0000    1.0000
            // 1.0000         0   -3.0000    2.0000
            // -0.1667    0.3333    1.0000   -0.8333
            PZMath_matrix inverse = matrix.LUInvert();
            inverse.ScreenWriteLine();
        } // LUExample()

        /// <summary>
        /// operator * example
        /// </summary>
        public static void OperatorMultiplyExample()
        {
            double[,] m1 = {{2, 3, 3, 5},
                           {6, 6, 8, 9}, 
                           {10, 11, 12, 13}};
            double[,] m2 = {{2, 3, 3},
                            {6, 6, 8}, 
                            {10, 11, 12},
                            {14, 15, 17}};
            PZMath_matrix matrix1 = new PZMath_matrix(m1);
            PZMath_matrix matrix2 = new PZMath_matrix(m2);
            // matrix multiply
            // correct answer
            // 122   132   151
            // 254   277   315
            // 388   423   483
            PZMath_matrix newmatrix = matrix1 * matrix2;
            newmatrix.ScreenWriteLine();
        } // OperatorMultiplyExample()

        /// <summary>
        /// LU Det example
        /// </summary>
        public static void LUDetExample()
        {
            double[,] m = {{2, 3, 3, 5},
                           {6, 6, 8, 9}, 
                           {10, 11, 12, 13},
                           {14, 15, 17, 17}};
            PZMath_matrix matrix = new PZMath_matrix(m);
            double det = matrix.LUDet();
            // correct answer:
            // det(m) = -12
            System.Console.WriteLine(det);
        } // LUDetExample()

        /// <summary>
        /// Cholsky decomposition example
        /// </summary>
        public static void CholeskyDecompExample()
        {
            double[,] m = {{1,    1,    1,    1,    1},
                           {1,    2,    3,    4,    5},
                           {1,    3,    6,   10,   15},
                           {1,    4,   10,   20,   35},
                           {1,    5,   15,   35,   70}};
            PZMath_matrix matrix = new PZMath_matrix(m);
            PZMath_matrix L = matrix.CholeskyDecomp();
            // correct answer:
            // L =
            // 1     0     0     0     0
            // 1     1     0     0     0
            // 1     2     1     0     0
            // 1     3     3     1     0
            // 1     4     6     4     1
            L.ScreenWriteLine();
        }
        #endregion
    } // PZMath_matrix
}
