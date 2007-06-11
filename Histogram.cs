// eee.Sheffield.PZ.Math
//
// Copyright ?Ping Zou, 2007
// sg71.cherub@gmail.com

using System;

namespace eee.Sheffield.PZ.Math
{    
    /// <summary>
    /// Histogram
    /// </summary>
    public class Histogram
    {
        #region Fields
        private int[] _values;
        private double _mean = 0.0;
        private double _stdDev = 0.0;
        private int _median = 0;
        private int _min;
        private int _max;

        private int _total = 0; 
        #endregion

        #region Properties
        public int[] Values { get { return _values; } }
        public double Mean { get { return _mean; } }
        public double StdDev { get { return _stdDev; } }
        public int Median { get { return _median; } }
        public int Min { get { return _min; } }
        public int Max { get { return _max; } } 
        #endregion

        // Constructor
        #region Constructor
		public Histogram(int[] values)
        {           
            int v;

            // get histogram array
            int length = values.Length;
            _values = new int[length];
            values.CopyTo(_values, 0);            
            
            // calculate mean, min, max
            _max = 0;
            _min = length - 1;
            for (int i = 0; i < length; i ++)
            {
                _total += _values[i];
            }
            int maxV = 0;
            int minV = _total;

            for (int i = 0; i < length; i++)
            {
                v = values[i];
                if (v > maxV)
                    _max = i;
                if (v < minV)
                    _min = i;
                _mean += i * v;
            }
            _mean /= _total;            

            // calculate stadard deviation
            for (int i = 0; i < length; i++)
            {
                v = _values[i];
                _stdDev += (float)System.Math.Pow(i - _mean, 2) * v;
            }
            _stdDev = (float)System.Math.Sqrt(_stdDev / _total);

            // calculate median
            int h = _total / 2;
            v = 0;
            for (_median = 0; _median < length; _median++)
            {
                v += _values[_median];
                if (v >= h)
                    break;
            }
        } // Histogram()
	    #endregion        
    }
}
