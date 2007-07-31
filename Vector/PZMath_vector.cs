// PZ Math
//
// -- vector, only hold double data
//
// -- Class PZMath_vector
// 24.11.2005


using System;
using System.IO;
using System.Collections;
using System.Collections.Generic;

namespace eee.Sheffield.PZ.Math
{
	/// <summary>
	/// vector class
	/// </summary>
	
	public class PZMath_vector
		// vector[i] = data[offset + i * stride];
		// for default: offset = 0 and stride = 1
	{
		internal double[] data = null;	// 1D array for contain data
		internal int length;		// length of data = offset + size * stride
		internal int size;		// size of vector
		internal int offset;		// offset of the start point of the data 
		internal int stride;		// stride of index increasing. 
        internal double _mean;   // mean of data
        internal double _stddev; // standard deviation

		#region constructor
		public PZMath_vector()
		{
			size = 0;
			length = 0;
			data = null;
			offset = 0;
			stride = 0;
            _mean = 0.0;
            _stddev = 0.0;
		} // PZMath_vector()
		
		public PZMath_vector(int s)
		{
			size = s;
			offset = 0;
			stride = 1;
			length = s;
			data = new double[length];
            _mean = 0.0;
            _stddev = 0.0;
		} // PZMath_vector(int size)

		public PZMath_vector(double [] d)
		{
			size = d.Length;
			offset = 0;
			stride = 1;
			length = size;
			data = new double [length];
            Array.Copy(d, data, length);
            _mean = CalculateMean();
            _stddev = StandardDeviation();
		} // PZMath_vector(double [] data)

		public PZMath_vector(PZMath_vector v)
			// copy constructor
		{
			size = v.Size;
			offset = v.offset;
			stride = v.stride;
			length = v.length;
			data = new double [length];
            Array.Copy(v.data, data, length);
            _mean = v._mean;
            _stddev = v._stddev;
		} // PZMath_vector(PZMath_vector vector)

        public PZMath_vector(ArrayList al)
        {
            size = al.Count;
            offset = 0;
            stride = 1;
            length = size;
            data = new double[length];
            for (int i = 0; i < length; i++)
            {
                data[i] = (double)al[i];
            }
            _mean = CalculateMean();
            _stddev = StandardDeviation();
        }// 

        public PZMath_vector(List<int> list)
        {
            size = list.Count;
            offset = 0;
            stride = 1;
            length = size;
            data = new double[length];
            for (int i = 0; i < length; i++)
                data[i] = (double)list[i];
            _mean = CalculateMean();
            _stddev = StandardDeviation();
        }

        public PZMath_vector(List<double> list)
        {
            size = list.Count;
            offset = 0;
            stride = 1;
            length = size;
            data = new double[length];
            for (int i = 0; i < length; i++)
                data[i] = (double)list[i];
            _mean = CalculateMean();
            _stddev = StandardDeviation();
        }
        #endregion

        #region properties
        public double Mean
        {
            get
            {
                _mean = CalculateMean();
                return _mean;
            }
        }

        public double StdDev
        {
            get
            {
                _stddev = StandardDeviation();
                return _stddev;
            }
        }

		public int Size
		{
			get {return size;}
		}

		public double this[int index] 
		{
			get 
			{
				if (index < 0 || index >= size) 
				{
					PZMath_errno.ERROR("Error: PZMath_vector::[], Index Out Of Range", PZMath_errno.PZMath_EFAILED);
				}
				return data[offset + index * stride];
			}
			set 
			{
				if (index < 0 || index >= size) 
				{
					PZMath_errno.ERROR("Error: PZMath_vector::[], Index Out Of Range", PZMath_errno.PZMath_EFAILED);
				}
				data[offset + index * stride] = value;
			}
		} // [] overriding
		#endregion

		#region set methods
		public void SetAll (double x)
		{
			for (int i = 0; i < size; i ++)
				this[i] = x;
		} // setall()

		public void SetZero()
		{
            Array.Clear(data, 0, length);
		} // SetZero()

		#endregion

		#region methods
		public int MinIndex()
			// This function returns the index of the minimum value in the vector v. 
			// When there are several equal minimum elements then the lowest index is returned.
		{
			double min = this[0];
			int index = 0;
			for (int i = 1; i < size; i ++)
			{
				if (this[i] < min)
				{
					min = this[i];
					index = i;
				}
			}
			return index;
		} // MinIndex

        public int MaxIndex()
            // This function returns the index of the maximum value in the vector v.
            // When there are several equal maximum elements then the lowest index is returned
        {
            double max = this[0];
            int index = 0;
            for (int i = 1; i < size; i++)
            {
                if (this[i] > max)
                {
                    max = this[i];
                    index = i;
                }
            }
            return index;
        } // MaxIndex()

		public void Swap(int index1, int index2)
		{
			if (index1 < 0 || index1 >= size
				|| index2 < 0 || index2 >= size) 
			{
				PZMath_errno.ERROR("Error: PZMath_vector::Swap(), Index Out Of Range", PZMath_errno.PZMath_EFAILED);
			}
			double t = data[offset + index1 * stride];
			data[offset + index1 * stride] = data[offset + index2 * stride];
			data[offset + index2 * stride] = t;
		} // swap()

		public PZMath_vector SubVector (int o, int newsize)
			// subvector, start from offset and have size = newsize
			// o -- offset
			// data is copy as reference;
		{
			if (o < 0 || o + newsize> size)
				PZMath_errno.ERROR("Error: PZMath_vector::SubVector(), Index Out Of Range", PZMath_errno.PZMath_EFAILED);
			PZMath_vector v = new PZMath_vector(this);
			v.size = newsize;
			v.offset = offset + o * stride;
			v.stride = stride;
			v.data = data;
			return v;
		} // SubVector()

		public void MemCopyFrom(PZMath_vector v)
			// v is the source vector
		{            
            if (data == null)
            {
                size = v.Size;
                offset = v.offset;
                stride = v.stride;
                length = v.length;
                if (v.data != null)
                {
                    data = new double[length];
                    Array.Copy(v.data, data, length);
                }
            }
            else
            {
                if (v.data == null)
                {
                    size = v.Size;
                    offset = v.offset;
                    stride = v.stride;
                    length = v.length;
                    data = null;
                }
                else
                {
                    for (int i = 0; i < size; i++)
                        this[i] = v[i];
                }
            }
            _mean = CalculateMean();
            _stddev = StandardDeviation();
		} // MemCopy()

		public bool Equals (PZMath_vector v)
			// Equals operator
		{
			if (size != v.Size)
				return false;
            else
                return Array.Equals(v.data, data);
		} // Equals()
		#endregion

        #region Operator override
        /// <summary>
        /// operator +
        /// </summary>
        /// <param name="v1"></param>
        /// <param name="v2"></param>
        /// <returns></returns>
        public static PZMath_vector operator +(PZMath_vector v1, PZMath_vector v2)
        {
            if (v1.Size != v2.Size)
                PZMath_errno.ERROR("PZMath_vecto::operator +, v1 and v2 are not match in dimension!");
            int length = v1.Size;
            PZMath_vector newVector = new PZMath_vector(length);
            for (int i = 0; i < length; i++)
                newVector[i] = v1[i] + v2[i];
            return newVector;
        } // operator +

        /// <summary>
        /// operator -
        /// </summary>
        /// <param name="v1"></param>
        /// <param name="v2"></param>
        /// <returns></returns>
        public static PZMath_vector operator -(PZMath_vector v1, PZMath_vector v2)
        {
            if (v1.Size != v2.Size)
                PZMath_errno.ERROR("PZMath_vecto::operator +, v1 and v2 are not match in dimension!");
            int length = v1.Size;
            PZMath_vector newVector = new PZMath_vector(length);
            for (int i = 0; i < length; i++)
                newVector[i] = v1[i] - v2[i];

            return newVector;
        } // operator -
        #endregion

        #region linear operators
        public int Div(PZMath_vector v)
        // This function divides the elements of vector a by the elements of vector b, 
        // a'_i = a_i / b_i. 
        // The two vectors must have the same length.
        // this / v
        {
            if (v.Size != size)
            {
                PZMath_errno.ERROR("PZMath_vector::Div(), vectors must have same length", PZMath_errno.PZMath_EBADLEN);
            }
            else
            {
                for (int i = 0; i < size; i++)
                {
                    this[i] /= v[i];
                }
                return PZMath_errno.PZMath_SUCCESS;
            }
            return 0;
        } // Div()

        public void Scale(double x)
        // multiply the element of vector by factor x
        {
            for (int i = 0; i < size; i++)
                this[i] *= x;
        } // Scale()

        /// <summary>
        /// convolve this with v
        /// </summary>
        /// <param name="v"></param>
        /// <returns></returns>
        public PZMath_vector Convolve(PZMath_vector x)
        {
            // 'this' is b;
            int nx = x.length;
            int nb = length;
            int ny = nx + nb - 1;
            PZMath_vector y = new PZMath_vector(ny);
            y.SetZero();

            for (int i = 0; i < nb; i++)
                for (int j = 0; j < nx; j++)
                    y[j + i] += x[j] * data[i];

            return y;
        } // Convolve()
        #endregion

        #region statistic methods
        public double CalculateMean()
        {
            double mean = 0.0;
            for (int i = 0; i < length; i ++)
                mean += data[i];
            mean /= (double)length;
            return mean;
        } // Mean()

        public double SampleStandardDeviation()
        {
            double mean = CalculateMean();
            double sumSquare = 0.0;
            for (int i = 0; i < length; i++)
                sumSquare += data[i] * data[i];
            double sampleStD = System.Math.Sqrt(sumSquare / (double)length - mean * mean);
            return sampleStD;
        } // SampleStandardDeviation()

        public double PopulationStandardDeviation()
        {
            double mean = CalculateMean();
            double populationStD = 0.0;
            for (int i = 0; i < length; i++)
                populationStD += (data[i] - mean) * (data[i] - mean);
            populationStD = System.Math.Sqrt(populationStD / (double)(length - 1));
            return populationStD;
        } // PopulationStandardDeviation()

        public double StandardDeviation()
        {
            return PopulationStandardDeviation();
        } // StandardDeviation()        

        public double SampleSkewness()
        {
            double m2 = CentralMoment(2);
            double m3 = CentralMoment(3);
            double sampleSkewness = m3 / System.Math.Pow(m2, 3.0 / 2.0);
            return sampleSkewness;
        } // SampleSkewness()

        public double PopulationSkewness()
        {
            double sampleSkewness = SampleSkewness();
            double n = (double)length;
            double populationSkewness = System.Math.Sqrt(n * (n - 1)) / (n - 2.0) * sampleSkewness;
            return populationSkewness;
        } // PopulationSkewness()

        // by default, return population statistics
        public double Skewness()
        {
            return PopulationSkewness();
        } // Skewness()

        public double SampleKurtosis()
        {
            double m4 = CentralMoment(4);
            double m2 = CentralMoment(2);
            double sampleKurtosis = m4 / m2 / m2 - 3;
            return sampleKurtosis;
        } // SampleKurtosis()

        public double PopulationKurtosis()
        {
            double m4 = CentralMoment(4);
            double m2 = CentralMoment(2);
            double n = (double)length;
            double populationKurtosis =
                ((n + 1) * m4 - 3 * (n - 1) * m2 * m2) * (n - 1) /
                (n - 2) / (n - 3) / m2 / m2;
            return populationKurtosis;
        } // PopulationKurtosis()

        // by default, return population statistics
        public double Kurtosis()
        {
            return PopulationKurtosis();
        } // Kurtosis()

        public double CentralMoment(int i)
        {
            if (i < 0)
                throw new ApplicationException("PZMath_vector::CentralMoment(), i < 0!");
            
            double mean = this.Mean;

            double centralMoment = 0;
            for (int j = 0; j < length; j++)
                centralMoment += System.Math.Pow((data[j] - mean), i);
            centralMoment /= (double)length;
            return centralMoment;
        } // CentralMoment()
        #endregion

        #region output methods
        public void DebugWriteLine()
			// out put vector data as Debug.WriteLine()
		{
			if (PZMath.PZDebug)
			{
				for (int x = 0; x < size; x ++)
				{
					System.Diagnostics.Debug.Write(String.Format(" {0, -4}", this[x]));
				}
				System.Diagnostics.Debug.WriteLine(" ");
			}
		} // DebugWriteLine()

        public void WriteFile(string fileName)
        {
            FileStream fileStream = new FileStream(fileName, FileMode.Create, FileAccess.Write);
            StreamWriter writer = new StreamWriter(fileStream);
            for (int i = 0; i < length; i++)
                writer.WriteLine(data[i]);
            writer.Flush();
            writer.Close();
            fileStream.Close();
        } // WriteFile
        #endregion
    } // PZMath_vector
}
