// eee.Sheffield.PZ.Math
//
// Copyright ?Ping Zou, 2007
// sg71.cherub@gmail.com

using System;
using System.Collections.Generic;
using System.Text;
using eee.Sheffield.PZ.Imaging;

namespace eee.Sheffield.PZ.Math
{
    /// <summary>
    /// target distribution function delegate
    /// </summary>
    /// <param name="x"></param>
    /// <returns></returns>
    public delegate double targetDistribution(double[] x);
    public delegate double targetPointProcessDistribution(List<PZPoint> x);
    public delegate double targetDensityFunction(LineSegmentConfiguration x);
    /// <summary>
    /// target distribution
    /// </summary>
    public class TargetDistribution
    {
        #region Fields
        public targetDistribution TargetDistributionFunction;
        public targetPointProcessDistribution TargetPointProcessDistributionFunction;
        public targetDensityFunction TargetDensityFunction;
        #endregion

        #region Constructor
        public TargetDistribution() { }
        #endregion

        #region Evaluate method
        //public double Evaluate(double[] x)
        //{
        //    if (TargetDistributionFunction == null)
        //        throw new ApplicationException("TargetDistribution::Evaluate(), target distribution function is not delegated.");
        //    return TargetDistributionFunction(x);
        //}

        //public double Evaluate(List<PZPoint> x)
        //{
        //    if (TargetPointProcessDistributionFunction == null)
        //        throw new ApplicationException("TargetDistribution::Evaluate(), target distribution function is not delegated.");
        //    return TargetPointProcessDistributionFunction(x);
        //}

        //public double Evaluate(LineSegmentConfiguration x)
        //{
        //    if (TargetDensityFunction == null)
        //        throw new ApplicationException("TargetDistribution::Evaluate(), target distribution function is not delegated.");
        //    return TargetDensityFunction(x);
        //}

        public double Evaluate(object x)
        {
            bool isSupported = true;

            if (x is double[])
            {
                if (TargetDistributionFunction == null)
                    throw new ApplicationException("TargetDistribution::Evaluate(), target distribution function is not delegated.");
                return TargetDistributionFunction((x as double[]));
            }
            else if (x is List<PZPoint>)
            {
                if (TargetPointProcessDistributionFunction == null)
                    throw new ApplicationException("TargetDistribution::Evaluate(), target distribution function is not delegated.");
                return TargetPointProcessDistributionFunction((x as List<PZPoint>));
            }
            else if (x is LineSegmentConfiguration)
            {
                if (TargetDensityFunction == null)
                    throw new ApplicationException("TargetDistribution::Evaluate(), target distribution function is not delegated.");
                return TargetDensityFunction((x as LineSegmentConfiguration));
            }
            else
            {
                isSupported = false;  
            }

            if (!isSupported)
                throw new ApplicationException("TargetDistribution::Evaluate(), cannot recognize target distribution function");
            return 0.0;
        } // Evaluate()
        #endregion
    }
}
