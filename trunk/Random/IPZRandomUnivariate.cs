// eee.Sheffield.PZ.Math
//
// Copyright ?Ping Zou, 2007
// sg71.cherub@gmail.com

using System;
using System.Collections.Generic;
using System.Text;

namespace eee.Sheffield.PZ.Math
{
    /// <summary>
    /// interface for univariate random
    /// </summary>
    public interface IPZRandomUnivariate
    {
        void ResetSeed(long seed);
        double Sample();
        double Evaluate(double x);
    }
}
