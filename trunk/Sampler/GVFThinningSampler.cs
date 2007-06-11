// eee.Sheffield.PZ.Math
//
// Copyright ?Ping Zou, 2007
// sg71.cherub@gmail.com

using System;
using System.Collections.Generic;
using System.Text;
using System.IO;

using eee.Sheffield.PZ.Imaging;

namespace eee.Sheffield.PZ.Math
{
    /// <summary>
    /// GVF Thining sampler, RJMCMC
    /// </summary>
    public class GVFThinningSampler
    {
        #region Fields
        // basic RJMCMC 
        private TargetDistribution _GMTargetDistribution;
        private LineSegmentConfiguration _lineSegmentPool;  // line segment pool, data term should have been assigned.               
        private LineSegmentConfiguration _sourceLineSegmentList;    // source line segment, equals to full configuration        
        private int _totalLineSegment;               // total # of line segments
       
        // RJMCMC parameters
        private int _n;     // number of MCMC iterations        
        private int _n0;    // number of burn-in        
        private List<LineSegmentConfiguration> _targetSamples;  // target samples, after burn-in
        private List<LineSegmentConfiguration> _totalSamples; // all samples, including burn-in
        
        // RJMCMC performance variable
        private double _sumAcceptProbability;       
        private double _averageAcceptProbability;        
        private int _hits;
                
        // reference Possion process
        private double _lamda;        
        private double _width;        
        private double _height;        
        private double _vS;   // v(S);
        
        // debug variables
        private double _sumn = 0.0;
        private double _averagen = 0.0;
        private double _sumnf = 0.0;
        private double _averagenf = 0.0;
        private double _sumnfe = 0.0;
        private double _averagenfe = 0.0;
        private double _sumpC = 0.0;
        private double _averagepC = 0.0;
        private double _sumpL = 0.0;
        private double _averagepL = 0.0;
        private double _sumpG = 0.0;
        private double _averagepG = 0.0;
        private double _sumpI = 0.0;
        private double _averagepI = 0.0;
        private double _sumhp = 0.0;
        private double _averagehp = 0.0;
        private double _sumhd = 0.0;
        private double _averagehd = 0.0;
        private double _sumh = 0.0;
        private double _averageh = 0.0;
        private int _indexMax = 0;        
        #endregion

        #region Properties
        public TargetDistribution GMTargetDistribution { get { return _GMTargetDistribution; } set { _GMTargetDistribution = value; } }
        public LineSegmentConfiguration LineSegmentPool { get { return _lineSegmentPool; } set { _lineSegmentPool = value; } }
        public LineSegmentConfiguration SourceLineSegmentList { get { return _sourceLineSegmentList; } set { _sourceLineSegmentList = value; } }        
        public int TotalLineSegment { get { return _totalLineSegment; } set { _totalLineSegment = value; } }
        // RJMCMC parameters
        public int N { get { return _n; } set { _n = value; } }
        public int N0 { get { return _n0; } set { _n0 = value; } }
        public List<LineSegmentConfiguration> TargetSamples { get { return _targetSamples; } set { _targetSamples = value; } }
        public List<LineSegmentConfiguration> TotalSamples { get { return _totalSamples; } set { _totalSamples = value; } }        
        // RJMCMC performance variables
        public double SumAcceptProbability { get { return _sumAcceptProbability; } set { _sumAcceptProbability = value; } }
        public double AverageAcceptProbability { get { return _averageAcceptProbability; } set { _averageAcceptProbability = value; } }
        public int Hits { get { return _hits; } set { _hits = value; } }
        // reference Possion process
        public double Lamda { get { return _lamda; } set { _lamda = value; } }
        public double Width { get { return _width; } set { _width = value; } }
        public double Height { get { return _height; } set { _height = value; } }
        public double VS { get { return _vS; } set { _vS = value; } }
        // debug variables        
        public double Sumn { get { return _sumn; } set { _sumn = value; } }
        public double Averagen { get { return _averagen; } set { _averagen = value; } }
        public double Sumnf { get { return _sumnf; } set { _sumnf = value; } }
        public double Averagenf { get { return _averagenf; } set { _averagenf = value; } }
        public double Sumnfe { get { return _sumnfe; } set { _sumnfe = value; } }        
        public double Averagenfe { get { return _averagenfe; } set { _averagenfe = value; } }        
        public double SumpC { get { return _sumpC; } set { _sumpC = value; } }
        public double AveragepC { get { return _averagepC; } set { _averagepC = value; } }        
        public double SumpL { get { return _sumpL; } set { _sumpL = value; } }
        public double AveragepL { get { return _averagepL; } set { _averagepL = value; } }
        public double SumpG { get { return _sumpG; } set { _sumpG = value; } }
        public double AveragepG { get { return _averagepG; } set { _averagepG = value; } }
        public double SumpI { get { return _sumpI; } set { _sumpI = value; } }
        public double AveragepI { get { return _averagepI; } set { _averagepI = value; } }
        public double Sumhp { get { return _sumhp; } set { _sumhp = value; } }
        public double Averagehp { get { return _averagehp; } set { _averagehp = value; } }
        public double Sumhd { get { return _sumhd; } set { _sumhd = value; } }
        public double Averagehd { get { return _averagehd; } set { _averagehd = value; } }
        public double Sumh { get { return _sumh; } set { _sumh = value; } }
        public double Averageh { get { return _averageh; } set { _averageh = value; } }        
        public int IndexMax { get { return _indexMax; } }       
        #endregion

        #region Constructor
        /// <summary>
        /// number of iterations, number of burn-in
        /// </summary>
        public GVFThinningSampler(int inputN, int inputn0, double lamda , double w, double h)
        {
            _n = inputN;
            _n0 = inputn0;
            int run = _n - _n0;
            if (_n0 >= _n)
                throw new ApplicationException("GVFThinningSampler::GVFThinningSampler(), burn-in length is great than sample length!");
            _targetSamples = new List<LineSegmentConfiguration>(run);
            _totalSamples = new List<LineSegmentConfiguration>(_n);
            _lamda = lamda;
            _width = w;
            _height = h;
        } // GVFThinningSampler()

        public GVFThinningSampler(int inputN, double l, double w, double h)
            : this(inputN, 0, l, w, h)
        {
        } // GVFThinningSampler()
        #endregion

        #region RJMCMC sampler GeyerMoller algorithm
        /// <summary>
        /// init MCMC
        /// data term has been assigned to the pool
        /// </summary>
        public void Init(LineSegmentConfiguration pool)
        {
            // source line segment list
            _sourceLineSegmentList = new LineSegmentConfiguration(pool);
            // line segment pool
            _lineSegmentPool = new LineSegmentConfiguration(pool);            
            _lineSegmentPool.Empty();
            // RJMCMC starts from full line segment configuration
            LineSegmentConfiguration startList = new LineSegmentConfiguration(pool);
            _totalSamples.Add(startList);

            // v(S)
            _vS = LebesgueMeasure(_lamda, _width, _height);
            _totalLineSegment = pool.N;
        } // Init()

        /// <summary>
        /// Lebesgue measure
        /// </summary>
        /// <returns></returns>
        private double LebesgueMeasure(double lamda, double width, double height)
        {
            return lamda * width * height;
        } // LebesgueMeasure()

        /// <summary>
        /// birth accept probability (birth accept rate)
        /// </summary>
        /// <returns></returns>
        private double BirthAcceptProbability(LineSegmentConfiguration x, LineSegmentConfiguration y)
        {
            // calculate birth accept probability 
            double R;
            double evaY = _GMTargetDistribution.Evaluate(y);
            double evaX = _GMTargetDistribution.Evaluate(x);

            R = evaY * _vS
                / evaX / y.N;
            //R = evaX / evaY;
            return R < 1.0 ? R : 1.0;
        } // BirthAcceptProbability()

        /// <summary>
        /// death accept probability
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        private double DeathAcceptProbability(LineSegmentConfiguration x, LineSegmentConfiguration y)
        {
            double R;
            double evaY = _GMTargetDistribution.Evaluate(y);
            double evaX = _GMTargetDistribution.Evaluate(x);

            R = evaY * x.N
                / evaX / _vS;
            //R = evaX / evaY;
            return R < 1.0 ? R : 1.0;
        } // DeathAcceptProbability()

        /// <summary>
        /// do accept rejection sampling
        /// </summary>
        /// <returns></returns>
        public List<LineSegmentConfiguration> Sample()
        {
            if (_GMTargetDistribution == null)
                throw new ApplicationException("GVFThinningSampler::Sample(), Candidate generating distribution or Instrumental distribution is not delegated.");

            // initial work space
            long seed = DateTime.Now.Millisecond;
            UniformDistribution randomSample = new UniformDistribution(seed + 10);  // sample for accept/reject
            double uSample;
            UniformDistribution randomPool = new UniformDistribution(seed + 30); // sample for birth from pool
            double uPool;
            UniformDistribution randomBirth = new UniformDistribution(seed + 40); // birth/death probability
            double uBirth;
            double uDeath;
            UniformDistribution randomRemove = new UniformDistribution(seed + 50);  // sample for remove line segment
            double uRemove;
            double acceptProbability;

            // MCMC interation
            int n = 0;
            _hits = 0;
            int Nm1 = _n - 1;
            double hMax = Double.MinValue;           
            while (n < _n)
            {                
                // get current coinfiguration
                LineSegmentConfiguration X = new LineSegmentConfiguration(_totalSamples[n]);
                
                int xLength = X.N;
                bool death = false;
                bool reject = false;
                // count hits of line segments
                int sourceLength = _sourceLineSegmentList.N;
                for (int s = 0; s < xLength; s++)
                {
                    int lineSegmentIndex = X.Configure[s].Index;
                    _sourceLineSegmentList.Configure[lineSegmentIndex].Hit++;
                    //for (int si = 0; si < sourceLength; si++)
                    //{
                    //    if (_sourceLineSegmentList.Configure[si].Index == lineSegmentIndex)
                    //    {
                    //        _sourceLineSegmentList.Configure[si].Hit++;
                    //        continue;
                    //    }
                    //}
                }

                /// Birth or Death ?
                #region Birth or Death ?
                uBirth = randomBirth.Sample();
                if (xLength > 0 && xLength < _totalLineSegment)
                    uDeath = 1 - uBirth;
                else if (xLength == 0)
                // only birth
                {
                    uBirth = 1.0;
                    uDeath = 0.0;
                }
                else
                // only death
                {
                    uDeath = 1.0;
                    uBirth = 0.0;
                } 
                #endregion

                if (uBirth > uDeath)
                // Birth
                #region Birth
                {
                    death = false;
                    // sample a new line segment from pool
                    uPool = System.Math.Floor(randomPool.Sample() * _lineSegmentPool.N);
                    LineSegment addedLineSegment = _lineSegmentPool.Configure[(int)uPool];

                    // propose Y
                    LineSegmentConfiguration Y = new LineSegmentConfiguration(X);
                    Y.Add(addedLineSegment);

                    // debug
                    #region debug
                    //double temphy = _GMTargetDistribution.Evaluate(Y);
                    //double temphx = _GMTargetDistribution.Evaluate(X);
                    //temphy = _GMTargetDistribution.Evaluate(Y);
                    //temphx = _GMTargetDistribution.Evaluate(X); 
                    #endregion

                    // accept probability
                    acceptProbability = BirthAcceptProbability(X, Y);
                    _sumAcceptProbability += acceptProbability;
                    uSample = randomSample.Sample();
                    if (uSample <= acceptProbability)
                    // accept Y
                    {
                        reject = false;
                        _hits++;
                        if (n < Nm1)
                            // add Y to totalSample;
                            _totalSamples.Add(Y);
                        if (n >= _n0)
                            // add Y to targetSample
                            _targetSamples.Add(Y);
                        // remove l from pool
                        _lineSegmentPool.RemoveAt((int)uPool);
                        //_totalLineSegment--;
                    }
                    else
                    // reject Y
                    {
                        reject = true;
                        if (n < Nm1)
                            _totalSamples.Add(X);
                        if (n >= _n0)
                            _targetSamples.Add(X);
                    }
                } 
                #endregion
                else
                // Death
                #region Death
                {
                    death = true;
                    // uniformly remove v from X
                    uRemove = (int)System.Math.Floor(randomRemove.Sample() * xLength);
                    // propose Y
                    LineSegmentConfiguration Y = new LineSegmentConfiguration(X);
                    LineSegment removedLineSegment = Y.Configure[(int)uRemove];
                    Y.RemoveAt((int)uRemove);

                    // accept probability
                    acceptProbability = DeathAcceptProbability(X, Y);
                    _sumAcceptProbability += acceptProbability;
                    uSample = randomSample.Sample();
                    if (uSample <= acceptProbability)
                    // accept Y
                    {
                        reject = false;
                        _hits++;
                        if (n < Nm1)
                            // add Y to totalSample;
                            _totalSamples.Add(Y);
                        if (n >= _n0)
                            // add Y to targetSample
                            _targetSamples.Add(Y);
                        // add l back to pool
                        _lineSegmentPool.Add(removedLineSegment);
                        //_totalLineSegment++;
                    }
                    else
                    // reject Y
                    {
                        reject = true;
                        if (n < Nm1)
                            _totalSamples.Add(X);
                        if (n >= _n0)
                            _targetSamples.Add(X);
                    }
                } 
                #endregion

                n++;

                // debug output
                #region debug
                _averageAcceptProbability = _sumAcceptProbability / (double)n;
                double h = _GMTargetDistribution.Evaluate(_totalSamples[n - 1]);
                if (h > hMax)
                {
                    hMax = h;
                    _indexMax = n - 1;
                }
                Console.Write("n = " + (n - 1));
                if (death)
                    Console.Write(" Death");
                else
                    Console.Write(" Birth");
                if (reject)
                    Console.Write(" Reject");
                else
                    Console.Write(" Accept");
                Console.Write(" # = " + _totalSamples[n - 1].N);
                Console.Write(" h(x) = " + String.Format("{0:0.0000e+00}", h));
                Console.Write(" accept rate = " + String.Format("{0:0.0000e+00}", acceptProbability));
                Console.WriteLine(" average accept rate " + String.Format("{0:f}", _averageAcceptProbability)); 
                #endregion
            } // end of while

            #region debug
            Console.Write("hx max = " + String.Format("{0:0.0000}", hMax) + " index = " + _indexMax);
            _averageAcceptProbability = _sumAcceptProbability / (double)_n; 
            #endregion

            return _targetSamples;
        } // Sample()
        #endregion

        #region Example
        /// <summary>
        /// </summary>
        public static void Example()
        {
            int N = 10000;
            int n0 = 1000;

            // file IO
            string filename1 = "v1.txt";
            FileStream fs1 = new FileStream(filename1, FileMode.Create, FileAccess.Write);
            StreamWriter w1 = new StreamWriter(fs1);

            // target distribution
            double lamda = 20;
            double r = 0.3;
            double width = 1.0;
            double height = 1.0;
            double theta1 = 2;
            double theta2 = -1.0;

            // prepare workspace
            GeyerMollerSampler GMSampler = new GeyerMollerSampler(N, n0, lamda, width, height);

            // prepare target distribution
            TargetDistribution targetDistribution = new TargetDistribution();
            targetDistribution.TargetDensityFunction = new targetDensityFunction(GVFThinningSampler.ExampleDensityFunction);

            // assign target distribution and instrumental distribution to the work space
            GMSampler.GMTargetDistribution = targetDistribution;

            // M-H sample
            GMSampler.Init();

            List<List<PZPoint>> samples = GMSampler.Sample();

            // output samples
            int sampleLength = samples.Count;
            for (int i = 0; i < sampleLength; i++)
            {
                int xLength = samples[i].Count;
                for (int j = 0; j < xLength; j++)
                {
                    w1.Write(String.Format("{0:f}", samples[i][j].x) + "  " + String.Format("{0:f}", samples[i][j].y) + " | ");
                }
                w1.WriteLine();
            }

            w1.Close();
            fs1.Close();
        } // Example()

        /// <summary>
        /// example target distribution function
        /// </summary>
        /// <returns></returns>
        public static double ExampleDensityFunction(LineSegmentConfiguration x)
        {
            // calculate hp
            x.FreeEndsFreeSegments();
            x.JoinSegments();
            x.CalculatePC();
            x.CalculateSumPC();
            x.CalculateTotalL();
            x.CalculateAverageL();

            // explicitly calculate sum(pL(s)), sum(pG(s)), and sum(pI(s)) for the efficiency sake
            int count = x.Configure.Count;
            double sumPL = 0.0;
            double sumPG = 0.0;
            double sumPI = 0.0;
            for (int s = 0; s < count; s++)
            {
                // calculate sum(pL(s))
                x.Configure[s].CalculatePL(x.TotalLength, x.AverageLength);
                // try different weight
                double pL = x.Configure[s].PL;
                if (pL < 0.0)
                    sumPL += pL * 100;
                else
                    sumPL += pL;
                x.SumPL = sumPL;
                
                // calculate sum(pG(s))
                x.Configure[s].CalculatePG(x.Gt);
                double pG = x.Configure[s].PG;
                if (pG < 0.0)
                    sumPG += pG * 100;
                else
                    sumPG += pG;
                x.SumPG = sumPG;

                // calculate sum(pI(s))
                x.Configure[s].CalculatePI(x.It, x.KI);
                double pI = x.Configure[s].PI;
                if (pI < 0.0)
                    sumPI += pI * 10000;
                else
                    sumPI += pI;
                x.SumPI = sumPI;
            }
            
            double up = x.W1 * x.N + x.W2 * x.NF + x.W3 * x.NFE +
                x.WC * x.SumPC + x.WL * x.SumPL;
            
            // calculate hd           
            double ud = x.WG * x.SumPG + x.WI * x.SumPI;
            
            // h = hp + hd;
            double u = up + ud;
            //return 1 / h;
            // density function U = exp(-h);
            double h = System.Math.Exp(-1.0 * u);            
            return h;
        } // ExampleTargetDistributionFunction()        
        #endregion
    }
}
