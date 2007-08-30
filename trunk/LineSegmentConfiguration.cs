// eee.Sheffield.PZ.Imaging
//
// Copyright ? Ping Zou, 2007
// sg71.cherub@gmail.com

using System;
using System.Collections.Generic;
using System.Text;
using eee.Sheffield.PZ.Math;
using System.Drawing;
using System.Drawing.Imaging;
using System.IO;
using CH.Combinations;

namespace eee.Sheffield.PZ.Imaging
{
    public class LineSegmentConfiguration
    {
        #region Field
        private List<LineSegment> _configuration = null;
        private List<LineSegment> _joinedConfiguration = null;

        // prior model
        private int _n = 0;     // # of current line segments
        private int _nfe = 0;   // # of free ends
        private int _nf = 0;   // # of free line segments
        private List<double> _pC = null;    // pc 
        private double[,] _pCPair = null;   // pc(si, sj) 
        public List<double> _cs = null;     // cs
       
        private double _totalL = 0.0;      // total length
        private double _averageL = 0.0;    // average length
        private double _sumPC = 0.0;
        private double _sumPL = 0.0;
        // data term           
        private double _sumPG = 0.0;
        private double _sumPI = 0.0;

        // data 
        private Bitmap _srcImage = null;
        private PZMath_matrix _gvfU = null;
        private PZMath_matrix _gvfV = null;
        private PZMath_matrix _gvfMagnitude = null;

        // parameters
        // prior model
        private double _w1 = 0.0;
        private double _w2 = 0.0;
        private double _w3 = 0.0;
        private double _wC = 0.0;
        private double _wL = 0.0;
        private double _thetaMax = System.Math.PI / 4.0;    // theta max
        // data term
        private double _wG = 0.0;  // wG
        private double _wI = 0.0;  // wI
        private int _It = 0; // intensity threshold
        private double _kI = 0.0;  // kI
        private double _Gt = 0.0;  // Gs threshold

        private double _Up = 0.0;
        private double _Ud = 0.0;
        private double _U = 0.0;
        private double _H = 0.0;

        #endregion

        #region Property
        public List<LineSegment> Configure { get { return _configuration; } }
        public List<LineSegment> JoinedConfiguration { get { return _joinedConfiguration; } set { _joinedConfiguration = value; } }

        public double[,] PCPair { get { return _pCPair; } set { _pCPair = value; } }
        // parameters       
        // prior model
        public double W1 { set { _w1 = value; } get { return _w1; } }
        public double W2 { set { _w2 = value; } get { return _w2; } }
        public double W3 { set { _w3 = value; } get { return _w3; } }
        public double WC { set { _wC = value; } get { return _wC; } }
        public double WL { set { _wL = value; } get { return _wL; } }
        public double ThetaMax { set { _thetaMax = value; } get { return _thetaMax; } }
        // data term
        public int It { set { _It = value; } get { return _It; } }
        public double KI { set { _kI = value; } get { return _kI; } }
        public double Gt { set { _Gt = value; } get { return _Gt; } }
        public double WG { set { _wG = value; } get { return _wG; } }
        public double WI { set { _wI = value; } get { return _wI; } }

        // results
        // prior model
        public int N { get { return _n; } }
        public int NFE { get { return _nfe; } }
        public int NF { get { return _nf; } }
        public double TotalLength { get { return _totalL; } }
        public double AverageLength { get { return _averageL; } }
        public double SumPC { get { return _sumPC; } set { _sumPC = value; } }
        public double SumPL { get { return _sumPL; } set { _sumPL = value; } }
        // data term
        public double SumPG { get { return _sumPG; } set { _sumPG = value; } }
        public double SumPI { get { return _sumPI; } set { _sumPI = value; } }

        // data
        public Bitmap Image { get { return _srcImage; } }
        public PZMath_matrix GVFU { get { return _gvfU; } }
        public PZMath_matrix GVFV { get { return _gvfV; } }
        public PZMath_matrix GVFMagnitude { get { return _gvfMagnitude; } }

        public double Ud { get { return _Ud; } set { _Ud = value; } }
        public double Up { get { return _Up; } set { _Up = value; } }
        public double U { get { return _U; } }
        public double H { get { return _H; } }
        #endregion

        #region Constructor
        /// <summary>
        /// constructor, input the line segment list
        /// </summary>
        /// <param name="lineSegmentList"></param>
        public LineSegmentConfiguration(List<LineSegment> lineSegmentList)
        {
            _configuration = new List<LineSegment>(lineSegmentList);
            _joinedConfiguration = new List<LineSegment>(lineSegmentList);
            RemoveShortLineSegment();
            _n = _configuration.Count;
        } // LineSegmentConfiguration()

        public LineSegmentConfiguration(LineSegmentConfiguration c)
        {
            _configuration = new List<LineSegment>(c._configuration);
            _joinedConfiguration = new List<LineSegment>(c._joinedConfiguration);

            _n = c._n;
            _nfe = c._nfe;
            _nf = c._nf;
            _pC = c._pC;
            _pCPair = c._pCPair;

            if (c._pC == null)
                _pC = null;
            else
                _pC = new List<double>(c._pC);
            if (c._pCPair == null)
                _pCPair = null;
            else
            {
                _pCPair = new double[_n, _n];
                Array.Copy(c._pCPair, _pCPair, c._pCPair.Length);
            }

            _totalL = c._totalL;
            _averageL = c._averageL;
            _sumPC = c._sumPC;
            _sumPL = c._sumPL;
            _sumPG = c._sumPG;
            _sumPI = c._sumPI;

            _srcImage = c._srcImage;
            _gvfU = c._gvfU;
            _gvfV = c._gvfV;
            _gvfMagnitude = c._gvfMagnitude;

            // prior model
            _w1 = c._w1;
            _w2 = c._w2;
            _w3 = c._w3;
            _wC = c._wC;
            _wL = c._wL;
            _thetaMax = c._thetaMax;
            // data term
            _It = c._It;
            _kI = c._kI;
            _Gt = c._Gt;
            _wG = c._wG;
            _wI = c._wI;

            _Up = c._Up;
            _Ud = c._Ud;
        } // LineSegmentConfiguration()
        #endregion

        #region Add/Remove method

        /// <summary>
        /// remove short (Ls <= 2 ) line segment
        /// </summary>
        public void RemoveShortLineSegment()
        {
            int index = 0;
            int lineSegmentCount = _configuration.Count;
            while (index < lineSegmentCount)
            {
                LineSegment l = (LineSegment)_configuration[index];
                if (l.PointList.Count <= 2)
                {
                    _configuration.RemoveAt(index);
                    lineSegmentCount--;
                }
                else
                {
                    index++;
                }
            }
        } // RemoveShortLineSegment()

        /// <summary>
        /// remove a line segment at its index
        /// </summary>
        /// <param name="index"></param>
        public void RemoveAtIndex(int index)
        {
            if (_n == 0)
                throw new ApplicationException("LineSegmentConfiguration::RemoveAt(), configuration contains no line segment.");
            bool removed = false;
            for (int s = 0; s < _n; s++)
            {
                if (_configuration[s].Index == index)
                {
                    _configuration.RemoveAt(s);
                    removed = true;
                    _n--;
                    break;
                }
            }
            if (!removed)
                throw new ApplicationException("LineSegmentConfiguration::RemoveAt(), the target has been removed or it does not exist.");
        } // RemoveAtIndex()

        public void RemoveAt(int index)
        {
            if (_n == 0)
                throw new ApplicationException("LineSegmentConfiguration::RemoveAt(), configuration contains no line segment.");
            if (index < 0 || index >= _n)
                throw new ApplicationException("LineSegmentConfiguration::RemoveAt(), index out of range.");
            _configuration.RemoveAt(index);
            _n--;
        } // RemoveAt

        /// <summary>
        /// add a line segment
        /// initialize it before adding
        /// </summary>
        /// <param name="l"></param>
        public void Add(LineSegment l)
        {
            // assigne data term
            //l.Initialize(_srcImage, _gvfU, _gvfV, _gvfMagnitude);
            //l.Index = _n;
            // add line segment
            _configuration.Add(l);
            _n++;
        } // Add()

        /// <summary>
        /// memory copy from the other line segment configuration
        /// </summary>
        /// <param name="c"></param>
        public void MemCopyFrom(LineSegmentConfiguration c)
        {
            _configuration = new List<LineSegment>(c._configuration);
            _joinedConfiguration = new List<LineSegment>(c._joinedConfiguration);
            _srcImage = c._srcImage;
            _gvfU = c._gvfU;
            _gvfV = c._gvfV;
            _gvfMagnitude = c._gvfMagnitude;

            if (c._pC == null)
                _pC = null;
            else
                _pC = new List<double>(c._pC);
            if (c._pCPair == null)
                _pCPair = null;
            else
            {
                _pCPair = new double[_n, _n];
                Array.Copy(c._pCPair, _pCPair, c._pCPair.Length);
            }

            _w1 = c._w1;
            _w2 = c._w2;
            _w3 = c._w3;
            _wC = c._wC;
            _wL = c._wL;
            _thetaMax = c._thetaMax;

            _It = c._It;
            _kI = c._kI;
            _Gt = c._Gt;
            _wG = c.WG;
            _wI = c._wI;

            _n = c._n;
            _nfe = c._nfe;
            _nf = c._nf;
            _totalL = c._totalL;
            _averageL = c._averageL;
            _sumPC = c._sumPC;
            _sumPL = c._sumPL;

            _sumPG = c._sumPG;
            _sumPI = c._sumPI;
            _Up = c._Up;
            _Ud = c._Ud;
        } // MemCopyFrom()

        /// <summary>
        /// empty line segments
        /// </summary>
        public void Empty()
        {
            // new fields
            _configuration.Clear();
            _joinedConfiguration.Clear();
            _pC = null;
            _pCPair = null;

            //if (_pC != null)
            //    _pC.Clear();
            //if (_pCPair != null)
            //    Array.Clear(_pCPair, 0, _pCPair.Length);

            _n = 0;
        } // Empty()
        #endregion

        #region prior model methods
        /// <summary>
        /// initialize line segment index
        /// </summary>
        /// <param name="lineSegmentList"></param>
        private void InitLineSegmentIndex(ref List<LineSegment> lineSegmentList)
        {
            int count = lineSegmentList.Count;
            for (int s = 0; s < count; s++)
                lineSegmentList[s].Index = s;
        } // InitLineSegmentIndex()
        public void InitLineSegmentIndex()
        {
            InitLineSegmentIndex(ref _configuration);
        } // InitLineSegmentIndex()

        private List<LineSegment> JoinSegments(List<LineSegment> lineSegmentList)
        {
            bool continueSearch = true;
            List<LineSegment> joinedLineSegmentList = new List<LineSegment>(lineSegmentList);
            if (joinedLineSegmentList.Count == 1)
                return joinedLineSegmentList;
            while (continueSearch)
            {
                int nl = joinedLineSegmentList.Count;
                int i;
                // for each line segment
                for (i = 0; i < nl; i++)
                {
                    LineSegment current = joinedLineSegmentList[i];
                    PZPoint currentStartPoint = current.StartPoint;
                    PZPoint currentEndPoint = current.EndPoint;
                    int nS = 0;
                    int nE = 0;
                    int indexS = 0, indexE = 0, indexRemove = 0;
                    // seach for the other segments
                    for (int j = 0; j < nl; j++)
                    {
                        if (i == j)
                            continue;
                        #region compare start point & end point
                        // compare start point & end point
                        LineSegment search = lineSegmentList[j];
                        if (currentStartPoint.Equals(search.StartPoint))
                        {
                            nS++;
                            indexS = j;
                        }
                        if (currentStartPoint.Equals(search.EndPoint))
                        {
                            nS++;
                            indexS = j;
                        }
                        if (currentEndPoint.Equals(search.EndPoint))
                        {
                            nE++;
                            indexE = j;
                        }
                        if (currentEndPoint.Equals(search.StartPoint))
                        {
                            nE++;
                            indexE = j;
                        }
                        #endregion
                    }
                    // merge or not
                    LineSegment newl;
                    if (nS == 0 && nE == 1)
                    {
                        indexRemove = indexE;
                        newl = current.MergeTo(joinedLineSegmentList[indexE]);
                    }
                    else if (nS == 1 && nE == 0)
                    {
                        indexRemove = indexS;
                        newl = current.MergeTo(joinedLineSegmentList[indexS]);
                    }
                    else
                        newl = null;
                    // generate new line segment list
                    if (newl != null)
                    {
                        // remove the current and the merged line segmente
                        joinedLineSegmentList.RemoveAt(i);
                        if (indexRemove < i)
                            joinedLineSegmentList.RemoveAt(indexRemove);
                        else
                            joinedLineSegmentList.RemoveAt(indexRemove - 1);
                        // add new line segment
                        joinedLineSegmentList.Add(newl);
                        // stop search
                        break;
                    }
                }
                if (i == nl)
                    continueSearch = false;
            }
            return joinedLineSegmentList;
        } // JoinSegments()
        public void JoinSegments()
        {
            _joinedConfiguration = JoinSegments(_configuration);
        } // JoinSegments()

        /// <summary>
        /// return # of free ends and # of free line segments
        /// </summary>
        /// <param name="lineSegmentList"></param>
        private void FreeEndsFreeSegments(List<LineSegment> lineSegmentList, out int fEnds, out int fLineSegments)
        {
            int length = lineSegmentList.Count;
            fEnds = 0;
            fLineSegments = 0;
            if (length == 1)
            {
                fEnds = 2;
                fLineSegments = 1;
            }
            else
            {
                // for each line segment
                for (int i = 0; i < length; i++)
                {
                    // check connectivity with the other line segment
                    LineSegment current = lineSegmentList[i];
                    bool startPointConnected = false;
                    bool endPointConnected = false;
                    for (int j = 0; j < length; j++)
                    {
                        if (i == j)
                            continue;
                        LineSegment search = lineSegmentList[j];
                        if (current.StartPoint.EqualTo(search.StartPoint)
                            || current.StartPoint.EqualTo(search.EndPoint))
                            startPointConnected = true;
                        if (current.EndPoint.EqualTo(search.StartPoint)
                            || current.EndPoint.EqualTo(search.EndPoint))
                            endPointConnected = true;
                    }
                    if (!startPointConnected)
                        fEnds++;
                    if (!endPointConnected)
                        fEnds++;
                    if (!startPointConnected && !endPointConnected)
                        fLineSegments++;
                }
            }
        } // FreeEndsFreeSegments()        
        public void FreeEndsFreeSegments()
        {
            _n = _configuration.Count;
            FreeEndsFreeSegments(_configuration, out _nfe, out _nf);
        } // FreeEndsFreeSegments()

        /// <summary>
        /// calculate pC(si, sj), before join!!!
        /// </summary>
        /// <param name="pc"></param>
        /// <param name="lineSegmentList"></param>
        private List<double> CalculatePC(List<LineSegment> lineSegmentList, double thetaMax, ref double[,] pCPair)
        {
            List<double> pc = new List<double>(0);
            _cs = new List<double>(0);

            //if (pCPair != null)
            //    Array.Clear(pCPair, 0, pCPair.Length);   
            pCPair = null;
            //pCPair = new double[lineSegmentList.Count, lineSegmentList.Count];
            double thetaij = 0;
            double pcsisj = 0.0;
            double pi = System.Math.PI;
            // search connected/joined pair
            int lineSegments = lineSegmentList.Count;
            int lineSegmentsM1 = lineSegments - 1;
            if (lineSegments == 1)
            {
                pc.Add(0.0);
                _cs.Add(0.0);
                //pCPair[0, 0] = 0.0;
            }
            else
            {
                #region search connected/joined pair
                
                // for each line segment
                for (int i = 0; i < lineSegmentsM1; i++)
                {
                    LineSegment si = lineSegmentList[i];
                    // seach for the other, half configuration
                    for (int j = i + 1; j < lineSegments; j++)
                    {
                        //if (i == j)
                        //    continue;
                        LineSegment sj = lineSegmentList[j];
                        // start point connected/joined
                        if (si.StartPoint.EqualTo(sj.StartPoint))
                        {
                            // thetas
                            thetaij = si.Cs.GetAngle(sj.Cs);
                            thetaij = System.Math.Min(thetaij, pi - thetaij);
                            // thetaij
                            //thetaij = si.DirectionList[0].GetAngle(sj.DirectionList[0]);
                            //thetaij = System.Math.Min(thetaij, pi - thetaij);

                            if (thetaij <= thetaMax)
                                pcsisj = -1.0 * Sigma(thetaij, thetaMax);
                            else
                                pcsisj = 1.0;
                            //if (Double.IsNaN(pcsisj))
                            //{
                            //    System.Console.Write("stop here");
                            //    thetaij = si.Cs.GetAngle(sj.Cs);
                            //}
                            pc.Add(pcsisj);
                            _cs.Add(thetaij);
                            //pCPair[i, j] = pcsisj;

                        }
                        if (si.StartPoint.EqualTo(sj.EndPoint))
                        {
                            // thetas
                            thetaij = si.Cs.GetAngle(sj.Cs);
                            thetaij = System.Math.Min(thetaij, pi - thetaij);
                            // thetaij
                            //thetaij = si.DirectionList[0].GetAngle(sj.DirectionList[sj.Length - 1]);
                            //thetaij = System.Math.Min(thetaij, pi - thetaij);

                            if (thetaij <= thetaMax)
                                pcsisj = -1.0 * Sigma(thetaij, thetaMax);
                            else
                                pcsisj = 1.0;
                            //if (Double.IsNaN(pcsisj))
                            //{
                            //    System.Console.Write("stop here");
                            //    thetaij = si.Cs.GetAngle(sj.Cs);
                            //}
                            pc.Add(pcsisj);
                            _cs.Add(thetaij);
                            //pCPair[i, j] = pcsisj;
                        }
                        // end point connected/joined
                        if (si.EndPoint.EqualTo(sj.StartPoint))
                        {
                            // thetas
                            thetaij = si.Cs.GetAngle(sj.Cs);
                            thetaij = System.Math.Min(thetaij, pi - thetaij);
                            // thetaij
                            //thetaij = si.DirectionList[si.Length - 1].GetAngle(sj.DirectionList[0]);
                            //thetaij = System.Math.Min(thetaij, pi - thetaij);
                            if (thetaij <= thetaMax)
                                pcsisj = -1.0 * Sigma(thetaij, thetaMax);
                            else
                                pcsisj = 1.0;
                            //if (Double.IsNaN(pcsisj))
                            //{
                            //    System.Console.Write("stop here");
                            //    thetaij = si.Cs.GetAngle(sj.Cs);
                            //}
                            pc.Add(pcsisj);
                            _cs.Add(thetaij);
                            //pCPair[i, j] = pcsisj;
                        }
                        if (si.EndPoint.EqualTo(sj.EndPoint))
                        {
                            // thetas
                            thetaij = si.Cs.GetAngle(sj.Cs);
                            thetaij = System.Math.Min(thetaij, pi - thetaij);
                            // thetaij
                            //thetaij = si.DirectionList[si.Length - 1].GetAngle(sj.DirectionList[sj.Length - 1]);
                            //thetaij = System.Math.Min(thetaij, pi - thetaij);

                            if (thetaij <= thetaMax)
                                pcsisj = -1.0 * Sigma(thetaij, thetaMax);
                            else
                                pcsisj = 1.0;
                            //if (Double.IsNaN(pcsisj))
                            //{
                            //    System.Console.Write("stop here");
                            //    thetaij = si.Cs.GetAngle(sj.Cs);
                            //}
                            pc.Add(pcsisj);
                            _cs.Add(thetaij);
                            //pCPair[i, j] = pcsisj;
                        }
                    } // for j
                } // for i 
                #endregion
            }
            return pc;
        } // CalculatePC()

        public void CalculatePC()
        {
            _pC = CalculatePC(_configuration, _thetaMax, ref _pCPair);
        } // CalculatePC()

        /// <summary>
        /// calculate sum pC(si, sj)
        /// </summary>
        /// <param name="pc"></param>
        /// <returns></returns>
        private double CalculateSumPC(List<double> pc)
        {
            int length = pc.Count;
            double sum = 0.0;

            //for (int i = 0; i < length; i++)
            //    sum += pc[i];

            // try different weight
            for (int i = 0; i < length; i++)
            {
                if (pc[i] < 0.0)
                    sum += pc[i] * 10;
                else
                    sum += pc[i];
            }
            // pC comes from half search
            return sum;
        } // CalculateSumPC()
        public void CalculateSumPC()
        {
            _sumPC = CalculateSumPC(_pC);
        } // CalculateSumPC()

        /// <summary>
        /// calculate total length, after join!!!
        /// </summary>
        /// <param name="lineSegmentList"></param>
        /// <returns></returns>
        private double CalculateTotalL(List<LineSegment> lineSegmentList)
        {
            double totalL = 0.0;
            int count = lineSegmentList.Count;
            // for each line segment
            for (int s = 0; s < count; s++)
                totalL += (double)lineSegmentList[s].Ls;
            return totalL;
        } // CalculateTotalL()
        public void CalculateTotalL()
        {
            _totalL = CalculateTotalL(_joinedConfiguration);
        } // CalculateTotalL()

        /// <summary>
        /// calculate average length, after join!!!
        /// </summary>
        /// <param name="n"></param>
        /// <param name="totalL"></param>
        /// <returns></returns>
        private double CalculateAverageL(int n, double totalL)
        {
            double averageL = totalL / (double)n;
            return averageL;
        } // CalculateAverageL()
        public void CalculateAverageL()
        {
            _averageL = CalculateAverageL(_n, _totalL);
        } // CalculateAverageL()

        /// <summary>
        /// calculate sum(pL(s)), after join!!!
        /// given total length and average length
        /// </summary>
        /// <param name="totalL"></param>
        /// <param name="averageL"></param>
        /// <returns></returns>
        private double CalculateSumPL(ref List<LineSegment> lineSegmentList, double totalL, double averageL)
        {
            double sumPL = 0.0;
            int count = lineSegmentList.Count;
            // for each line segment
            for (int s = 0; s < count; s++)
            {
                // calculate PL
                lineSegmentList[s].CalculatePL(totalL, averageL);
                // try different weight
                double pL = lineSegmentList[s].PL;
                if (pL < 0.0)
                    sumPL += pL * 100;
                else
                    sumPL += pL;
                //sumPL += lineSegmentList[s].PL;
            }
            return sumPL;
        } // CalculateSumPL()
        public void CalculateSumPL()
        {
            _sumPL = CalculateSumPL(ref _joinedConfiguration, _totalL, _averageL);
        } // CalculateSumPL()

        /// <summary>
        /// calculate Up(S)
        /// </summary>
        private double CalculateUP(
            List<LineSegment> lineSegmentList, out List<LineSegment> joinedLineSegmentList,
            double thetaMax, double w1, double w2, double w3, double wC, double wL,
            out int n, out int nf, out int nfe,
            out double totalL, out double averageL, out double sumPC, out double sumPL,
            ref List<double> pC, ref double[,] pCPair)
        {
            int count = lineSegmentList.Count;

            if (count == 0)
            {
                n = 0;
                nfe = 0;
                nf = 0;
                totalL = 0.0;
                averageL = 0.0;
                sumPC = 0.0;
                sumPL = 0.0;
                pC = new List<double>(0);
                pCPair = new double[1, 1];
                joinedLineSegmentList = new List<LineSegment>(0);
                double hp = 0.0;
                return hp;
            }
            else
            {
                // initialize the prior model of each line segments;
                // so that Cs and Ls are ready                
                // for each line segment
                for (int s = 0; s < count; s++)
                {
                    lineSegmentList[s].InitializePriorModel();
                    //lineSegmentList[s].Index = s;
                }

                // # of segments, free segments, and free end
                n = count;
                FreeEndsFreeSegments(lineSegmentList, out nfe, out nf);
                // calculate pC(si,sj)
                //if (pC != null)
                //    pC.Clear();
                pC = CalculatePC(lineSegmentList, thetaMax, ref pCPair);
                // calculate sum pc
                sumPC = CalculateSumPC(pC);
                // join segment
                joinedLineSegmentList = JoinSegments(lineSegmentList);
                // total length & average length
                totalL = CalculateTotalL(joinedLineSegmentList);
                averageL = CalculateAverageL(n, totalL);
                sumPL = CalculateSumPL(ref joinedLineSegmentList, totalL, averageL);

                // calculate hp(S)
                double up = w1 * n + w2 * nf + w3 * nfe
                    + wC * sumPC + wL * sumPL;
                return up;
            }
        } // CalculateUP()        
        public void CalculateUP()
        {
            _Up = CalculateUP(
                _configuration, out _joinedConfiguration,
                _thetaMax, _w1, _w2, _w3, _wC, _wL,
                out _n, out _nf, out _nfe,
                out _totalL, out _averageL, out _sumPC, out _sumPL,
                ref _pC, ref _pCPair);
        } // CalculateUP()
        #endregion

        #region data term methods
        /// <summary>
        /// assign data term to line segments
        /// and initialize their data term so that Is and Gs are ready
        /// </summary>
        /// <param name="srcImage"></param>
        /// <param name="gvfU"></param>
        /// <param name="gvfV"></param>
        /// <param name="gvfMagnitude"></param>
        private void AssignDataTerm(ref List<LineSegment> lineSegmentList, Bitmap srcImage, PZMath_matrix gvfU, PZMath_matrix gvfV, PZMath_matrix gvfMagnitude)
        {
            int count = lineSegmentList.Count;
            // for each line segment
            for (int s = 0; s < count; s++)
            {
                lineSegmentList[s].InitializePriorModel();
                lineSegmentList[s].InitializeDataTerm(srcImage, gvfU, gvfV, gvfMagnitude);
            }
        } // AssignDataTerm()
        public void AssignDataTerm(Bitmap srcImage, PZMath_matrix gvfU, PZMath_matrix gvfV, PZMath_matrix gvfMagnitude)
        {
            _srcImage = srcImage;
            _gvfU = gvfU;
            _gvfV = gvfV;
            _gvfMagnitude = gvfMagnitude;
            AssignDataTerm(ref _configuration, _srcImage, _gvfU, _gvfV, _gvfMagnitude);
        } // AssignDataTerm()

        /// <summary>
        /// calculate sum(pG(s))
        /// </summary>
        /// <param name="lineSegmentList"></param>
        /// <returns></returns>
        private double CalculateSumPG(ref List<LineSegment> lineSegmentList, double Gt)
        {
            double sumPG = 0.0;
            int count = lineSegmentList.Count;
            // for each line segment
            for (int s = 0; s < count; s++)
            {
                // calculate pG(s)
                lineSegmentList[s].CalculatePG(Gt);
                // try different weight
                double pG = lineSegmentList[s].PG;
                if (pG < 0.0)
                    sumPG += pG * 100;
                else
                    sumPG += pG;

                //sumPG += lineSegmentList[s].PG;
            }

            return sumPG;
        } // CalculateSumPG()
        public void CalculateSumPG()
        {
            _sumPG = CalculateSumPG(ref _configuration, _Gt);
        } // CalculateSumPG()

        /// <summary>
        /// calculate sum(pI(s))
        /// </summary>
        /// <param name="lineSegmentList"></param>
        /// <returns></returns>
        private double CalculateSumPI(ref List<LineSegment> lineSegmentList, double It, double kI)
        {
            double sumPI = 0.0;
            int count = lineSegmentList.Count;
            // for each line segment
            for (int s = 0; s < count; s++)
            {
                // calculate pI(s)
                lineSegmentList[s].CalculatePI(It, kI);
                // try different weight
                double pI = lineSegmentList[s].PI;
                if (pI < 0.0)
                    sumPI += pI * 10000;
                else
                    sumPI += pI;

                //sumPI += lineSegmentList[s].PI;
            }
            return sumPI;
        } // CalculateSumPI()
        public void CalculateSumPI()
        {
            _sumPI = CalculateSumPI(ref _configuration, _It, _kI);
        } // CalculateSumPI()

        /// <summary>
        /// calculate Ud(S)
        /// </summary>        
        /// <returns></returns>
        private double CalculateUD(ref List<LineSegment> lineSegmentList,
            Bitmap srcImage, PZMath_matrix gvfU, PZMath_matrix gvfV, PZMath_matrix gvfMagnitude,
            double wG, double wI,
            double Gt, double It, double kI,
            out double sumPG, out double sumPI)
        {
            // assign data term to line segments
            AssignDataTerm(ref lineSegmentList, srcImage, gvfU, gvfV, gvfMagnitude);
            // calculate sum(pG(s))
            sumPG = CalculateSumPG(ref lineSegmentList, Gt);
            // calculate sum(pI(s))
            sumPI = CalculateSumPI(ref lineSegmentList, It, kI);
            double ud = 0;
            ud = wG * sumPG + wI * sumPI;
            return ud;
        } // CalculateUD()
        public void CalculateUD()
        {
            _Ud = CalculateUD(ref _configuration, _srcImage, _gvfU, _gvfV, _gvfMagnitude, _wG, _wI, _Gt, _It, _kI, out _sumPG, out _sumPI);
        } // CalculateUD()

        /// <summary>
        /// calculate U(S) = Up(S) + Ud(S)
        /// </summary>
        /// <param name="hp"></param>
        /// <param name="hd"></param>
        /// <returns></returns>
        private double CalculateU(double up, double ud)
        {
            double u = up + ud;
            return u;
        }
        public void CalculateU()
        {
            _U = CalculateU(_Up, _Ud);
        } // CalculateU()

        /// <summary>
        /// calculate H(S) = exp(-U(s))
        /// </summary>
        /// <param name="u"></param>
        /// <returns></returns>
        private double CalculateH(double u)
        {
            double h = System.Math.Exp(-1.0 * u);
            return h;
        } // CalculateH()
        public void CalculateH()
        {
            _H = CalculateH(_U);
        } // CalculateH()
        #endregion

        #region util
        /// <summary>
        /// Sigma(x, M)
        /// </summary>
        /// <param name="x"></param>
        /// <param name="M"></param>
        /// <returns></returns>
        private double Sigma(double x, double M)
        {
            return System.Math.Exp(-1.0 * 10 * x / M);
        }

        /// <summary>
        /// clear loops, remain the one with lower Gs value
        /// </summary>
        public void ClearLoops()
        {
            int nM1 = _n - 1;
            List<int> removeIndexList = new List<int>(0);
            // search loops, record the one with higher Gs value
            for (int i = 0; i < nM1; i++)
            // for each line segment
            {
                LineSegment lineSegmenti = (LineSegment)_configuration[i];
                for (int j = i + 1; j < _n; j++)
                {
                    LineSegment lineSegmentj = (LineSegment)_configuration[j];
                    if ((lineSegmenti.StartPoint.EqualTo(lineSegmentj.StartPoint)
                        && lineSegmenti.EndPoint.EqualTo(lineSegmentj.EndPoint))
                        || (lineSegmenti.StartPoint.EqualTo(lineSegmentj.EndPoint)
                        && lineSegmenti.EndPoint.EqualTo(lineSegmentj.StartPoint)))
                    {
                        if (lineSegmenti.Gs > lineSegmentj.Gs)
                            removeIndexList.Add(i);
                        else
                            removeIndexList.Add(j);
                    }
                }
            }

            // remove 
            int removeCount = removeIndexList.Count;
            removeIndexList.Sort();
            for (int i = removeCount - 1; i >= 0; i--)
                Configure.RemoveAt(removeIndexList[i]);
            _n = _configuration.Count;
        } // ClearLoops()

        /// <summary>
        /// clear free line segments
        /// </summary>
        public void ClearFreeLineSegments()
        {
            bool startFree = true;
            bool endFree = true;
            List<int> removeIndex = new List<int>(0);
            List<LineSegment> joinedLineSegmentConfiguration = JoinSegments(_configuration);
            int lineSegmentCount = joinedLineSegmentConfiguration.Count;
            for (int i = 0; i < lineSegmentCount - 1; i++)
            {
                LineSegment lineSegmenti = (LineSegment)joinedLineSegmentConfiguration[i];
                for (int j = i + 1; j < lineSegmentCount; j++)
                {
                    LineSegment lineSegmentj = (LineSegment)joinedLineSegmentConfiguration[j];

                    // search start point
                    if (lineSegmenti.StartPoint.EqualTo(lineSegmentj.StartPoint)
                        || lineSegmenti.StartPoint.EqualTo(lineSegmentj.EndPoint))
                    {
                        startFree = false;
                        break;
                    }
                    // search end point
                    if (lineSegmenti.EndPoint.EqualTo(lineSegmentj.StartPoint)
                        || lineSegmenti.EndPoint.EqualTo(lineSegmentj.EndPoint))
                    {
                        endFree = false;
                        break;
                    }
                }
                if (startFree && endFree)
                    removeIndex.Add(i);
            }

            // remove free line segments
            int removeCount = removeIndex.Count;
            removeIndex.Sort();
            for (int i = removeCount - 1; i >= 0; i--)
                joinedLineSegmentConfiguration.RemoveAt(removeIndex[i]);

            _configuration = joinedLineSegmentConfiguration;
        } // ClearFreeLineSegments
        #endregion

        #region I/O
        /// <summary>
        /// draw line segment
        /// </summary>
        private Bitmap DrawLineSegments(Bitmap srcImage, List<LineSegment> lineSegmentList)
        {
            int width = srcImage.Width;
            int height = srcImage.Height;
            Bitmap dstImage = new Bitmap(width, height, PixelFormat.Format24bppRgb);

            #region set background white
            // lock destination bitmap data
            BitmapData dstData = dstImage.LockBits(
                new Rectangle(0, 0, width, height),
                ImageLockMode.ReadWrite, PixelFormat.Format24bppRgb);

            int stride = dstData.Stride;
            int offset = stride - width * 3;
            unsafe
            {
                byte* dst = (byte*)dstData.Scan0.ToPointer();
                for (int y = 0; y < height; y++)
                {
                    for (int x = 0; x < width; x++, dst += 3)
                    {
                        dst[0] = 255;
                        dst[1] = 255;
                        dst[2] = 255;
                    }
                    dst += offset;
                }
            }
            dstImage.UnlockBits(dstData);
            #endregion

            // draw line segment on dstImg
            int length = lineSegmentList.Count;
            for (int i = 0; i < length; i++)
            {
                LineSegment current = lineSegmentList[i];
                int points = current.Ls;
                for (int j = 0; j < points; j++)
                {
                    PZPoint p = current.PointList[j];
                    int x = (int)p.x;
                    int y = (int)p.y;
                    dstImage.SetPixel(x, y, Color.Blue);
                    dstImage.SetPixel((int)current.StartPoint.x, (int)current.StartPoint.y, Color.Red);
                    dstImage.SetPixel((int)current.EndPoint.x, (int)current.EndPoint.y, Color.Red);
                }
            }
            return dstImage;
        } // DrawLineSegment()
        public Bitmap DrawLineSegments(Bitmap srcImage)
        {
            return DrawLineSegments(srcImage, _configuration);
        } // DrawLineSegments()

        public PZMath_matrix DrawGsMatrix(Bitmap srcImage)
        {
            foreach (LineSegment l in _configuration)
                l.CalculateGs();
            return DrawGsMatrix(srcImage, _configuration);
        }

        private PZMath_matrix DrawGsMatrix(Bitmap srcImage, List<LineSegment> lineSegmentList)
        {
            int width = srcImage.Width;
            int height = srcImage.Height;

            PZMath_matrix dstMatrix = new PZMath_matrix(height, width);
            dstMatrix.Setall(1.0);

            foreach (LineSegment l in lineSegmentList)
            {
                int points = l.Ls;
                for (int j = 0; j < points; j++)
                {
                    PZPoint p = l.PointList[j];
                    int x = (int)p.x;
                    int y = (int)p.y;
                    dstMatrix[y, x] = l.Gs;
                }
            }

            return dstMatrix;
        } // DrawGsMatrix()

        // draw gs threshold line segments
        private Bitmap DrawGsThreshold(Bitmap srcImage, List<LineSegment> lineSegmentList, double gt)
        {
            int width = srcImage.Width;
            int height = srcImage.Height;    
            Bitmap dstImage = new Bitmap(width, height, PixelFormat.Format24bppRgb);

            #region set background white
            // lock destination bitmap data
            BitmapData dstData = dstImage.LockBits(
                new Rectangle(0, 0, width, height),
                ImageLockMode.ReadWrite, PixelFormat.Format24bppRgb);

            int stride = dstData.Stride;
            int offset = stride - width * 3;
            unsafe
            {
                byte* dst = (byte*)dstData.Scan0.ToPointer();
                for (int y = 0; y < height; y++)
                {
                    for (int x = 0; x < width; x++, dst += 3)
                    {
                        dst[0] = 255;
                        dst[1] = 255;
                        dst[2] = 255;
                    }
                    dst += offset;
                }
            }
            dstImage.UnlockBits(dstData);
            #endregion

            // draw line segment on dstImg
            int length = lineSegmentList.Count;
            for (int i = 0; i < length; i++)
            {
                LineSegment current = lineSegmentList[i];
                if (current.Gs > gt)
                    continue;
                int points = current.Ls;
                for (int j = 0; j < points; j++)
                {
                    PZPoint p = current.PointList[j];
                    int x = (int)p.x;
                    int y = (int)p.y;
                    dstImage.SetPixel(x, y, Color.Blue);
                    dstImage.SetPixel((int)current.StartPoint.x, (int)current.StartPoint.y, Color.Red);
                    dstImage.SetPixel((int)current.EndPoint.x, (int)current.EndPoint.y, Color.Red);
                }
            }
            return dstImage;
        }
        public Bitmap DrawGsThreshold(Bitmap srcImage, double gt)
        {
            int count = _configuration.Count;
            for (int s = 0; s < count; s++)
                _configuration[s].CalculateGs();
            return DrawGsThreshold(srcImage, _configuration, gt);
        }

        private Bitmap DrawGsThreshold(Bitmap srcImage, List<LineSegment> lineSegmentList)
        {
            int width = srcImage.Width;
            int height = srcImage.Height;
            Bitmap dstImage = new Bitmap(width, height, PixelFormat.Format24bppRgb);

            #region set background white
            // lock destination bitmap data
            BitmapData dstData = dstImage.LockBits(
                new Rectangle(0, 0, width, height),
                ImageLockMode.ReadWrite, PixelFormat.Format24bppRgb);

            int stride = dstData.Stride;
            int offset = stride - width * 3;
            unsafe
            {
                byte* dst = (byte*)dstData.Scan0.ToPointer();
                for (int y = 0; y < height; y++)
                {
                    for (int x = 0; x < width; x++, dst += 3)
                    {
                        dst[0] = 255;
                        dst[1] = 255;
                        dst[2] = 255;
                    }
                    dst += offset;
                }
            }
            dstImage.UnlockBits(dstData);
            #endregion

            // draw line segment on dstImg
            int length = lineSegmentList.Count;
            for (int i = 0; i < length; i++)
            {
                LineSegment current = lineSegmentList[i];
                int points = current.Ls;
                for (int j = 0; j < points; j++)
                {
                    PZPoint p = current.PointList[j];
                    int x = (int)p.x;
                    int y = (int)p.y;
                    Color lineColour;
                    if (current.Gs < 0.45)
                        lineColour = Color.Black;
                    else if (current.Gs < 0.55)
                        lineColour = Color.Blue;
                    else if (current.Gs < 0.65)
                        lineColour = Color.Green;
                    else
                        lineColour = Color.Red;
                    dstImage.SetPixel(x, y, lineColour);
                    dstImage.SetPixel((int)current.StartPoint.x, (int)current.StartPoint.y, Color.Red);
                    dstImage.SetPixel((int)current.EndPoint.x, (int)current.EndPoint.y, Color.Red);
                }
            }
            return dstImage;
        }
        public Bitmap DrawGsThreshold(Bitmap srcImage)
        {
            int count = _configuration.Count;
            for (int s = 0; s < count; s++)
                _configuration[s].CalculateGs();
            return DrawGsThreshold(srcImage, _configuration);
        }

        // draw Is threshold line segments
        private Bitmap DrawIsThreshold(Bitmap srcImage, List<LineSegment> lineSegmentList, double It)
        {
            int width = srcImage.Width;
            int height = srcImage.Height;
            Bitmap dstImage = new Bitmap(width, height, PixelFormat.Format24bppRgb);

            #region set background white
            // lock destination bitmap data
            BitmapData dstData = dstImage.LockBits(
                new Rectangle(0, 0, width, height),
                ImageLockMode.ReadWrite, PixelFormat.Format24bppRgb);

            int stride = dstData.Stride;
            int offset = stride - width * 3;
            unsafe
            {
                byte* dst = (byte*)dstData.Scan0.ToPointer();
                for (int y = 0; y < height; y++)
                {
                    for (int x = 0; x < width; x++, dst += 3)
                    {
                        dst[0] = 255;
                        dst[1] = 255;
                        dst[2] = 255;
                    }
                    dst += offset;
                }
            }
            dstImage.UnlockBits(dstData);
            #endregion

            // draw line segment on dstImg
            int length = lineSegmentList.Count;
            for (int i = 0; i < length; i++)
            {
                LineSegment current = lineSegmentList[i];
                if (current.Is > It)
                    continue;
                int points = current.Ls;
                for (int j = 0; j < points; j++)
                {
                    PZPoint p = current.PointList[j];
                    int x = (int)p.x;
                    int y = (int)p.y;
                    dstImage.SetPixel(x, y, Color.Blue);
                    dstImage.SetPixel((int)current.StartPoint.x, (int)current.StartPoint.y, Color.Red);
                    dstImage.SetPixel((int)current.EndPoint.x, (int)current.EndPoint.y, Color.Red);
                }
            }
            return dstImage;
        }
        public Bitmap DrawIsThreshold(Bitmap srcImage, double It)
        {
            int count = _configuration.Count;
            for (int s = 0; s < count; s++)
                _configuration[s].CalculateIs();
            return DrawIsThreshold(srcImage, _configuration, It);
        }

        public void WriteLineSegmentsFile(string fileName)
        {
            FileStream fileStream = new FileStream(fileName, FileMode.Create, FileAccess.Write);
            StreamWriter writer = new StreamWriter(fileStream);
            // for each line segment
            for (int i = 0; i < _n; i++)
            {
                LineSegment l = _configuration[i];
                int length = l.Ls;
                // for each point
                for (int j = 0; j < length; j++)
                {
                    writer.WriteLine(String.Format("{0,20:f}{1,20:f}", l.PointList[j].x, l.PointList[j].y));
                }
            }
            writer.Flush();
            writer.Close();
            fileStream.Close();
        }

        public void WriteEndPointsFile(string fileName)
        {
            FileStream fileStream = new FileStream(fileName, FileMode.Create, FileAccess.Write);
            StreamWriter writer = new StreamWriter(fileStream);
            // for each line segment
            for (int i = 0; i < _n; i++)
            {
                LineSegment l = _configuration[i];
                int length = l.Ls;
                // for each end point
                for (int j = 0; j < length; j++)
                {
                    writer.WriteLine(String.Format("{0,20:f}{1,20:f}", l.StartPoint.x, l.StartPoint.y));
                    writer.WriteLine(String.Format("{0,20:f}{1,20:f}", l.EndPoint.x, l.EndPoint.y));
                }
            }
            writer.Flush();
            writer.Close();
            fileStream.Close();
        }

        public void WriteGsFile(string fileName)
        {
            FileStream fileStream = new FileStream(fileName, FileMode.Create, FileAccess.Write);
            StreamWriter writer = new StreamWriter(fileStream);
            // for each line segment
            for (int s = 0; s < _n; s++)
            {
                LineSegment l = _configuration[s];
                int length = l.Ls;
                // for each end point
                for (int i = 0; i < length; i++)
                    writer.WriteLine(String.Format("{0,20:0.0000}{1,20:0.0000}{2,20:0.0000}", l.PointList[i].x, l.PointList[i].y, l.Gs));
            }
            writer.Flush();
            writer.Close();
            fileStream.Close();
        }

        public void WriteGsHistogram(string fileName)
        {
            FileStream fileStream = new FileStream(fileName, FileMode.Create, FileAccess.Write);
            StreamWriter writer = new StreamWriter(fileStream);
            // for each line segment
            foreach (LineSegment l in _configuration)
                writer.WriteLine(l.Gs);
            writer.Flush();
            writer.Close();
            fileStream.Close();
        }

        public void WriteGsListHistogram(string fileName)
        {
            FileStream fileStream = new FileStream(fileName, FileMode.Create, FileAccess.Write);
            StreamWriter writer = new StreamWriter(fileStream);
            // for each line segment
            foreach (LineSegment l in _configuration)
            {
                int i = 1;  // no junctions/ends
                int countM1 = l._GsList.Count - 1;
                for (; i < countM1; i++)
                    writer.WriteLine(l._GsList[i]);
            }
            writer.Flush();
            writer.Close();
            fileStream.Close();
        } // WriteGsListHistogram()

        public void WriteAngleHistogram(string fileName)
        {
            FileStream fileStream = new FileStream(fileName, FileMode.Create, FileAccess.Write);
            StreamWriter writer = new StreamWriter(fileStream);
            // for each line segment
            foreach (LineSegment l in _configuration)
                writer.WriteLine(l._angle);
            writer.Flush();
            writer.Close();
            fileStream.Close();
        } // WriteAngleHistogram()

        public void WriteAngleListHistogram(string fileName)
        {
            FileStream fileStream = new FileStream(fileName, FileMode.Create, FileAccess.Write);
            StreamWriter writer = new StreamWriter(fileStream);
            // for each line segment            
            foreach (LineSegment l in _configuration)
            {
                int i = 1;  // no junctions/ends 
                int countM1 = l._angleList.Count - 1;
                for (; i < countM1; i++)
                    writer.WriteLine(l._angleList[i]);
            }
            writer.Flush();
            writer.Close();
            fileStream.Close();
        } // WriteAngleListHistogram

        public void WritePGFile(string fileName)
        {
            FileStream fileStream = new FileStream(fileName, FileMode.Create, FileAccess.Write);
            StreamWriter writer = new StreamWriter(fileStream);
            // for each line segment
            for (int s = 0; s < _n; s++)
            {
                LineSegment l = _configuration[s];
                int length = l.Ls;
                // for each end point
                for (int i = 0; i < length; i++)
                    writer.WriteLine(String.Format("{0,20:0.0000}{1,20:0.0000}{2,20:0.0000}", l.PointList[i].x, l.PointList[i].y, l.PG));
            }
            writer.Flush();
            writer.Close();
            fileStream.Close();
        }

        public void WriteCsHistogram(string fileName)
        {
            FileStream fileStream = new FileStream(fileName, FileMode.Create, FileAccess.Write);
            StreamWriter writer = new StreamWriter(fileStream);        
            foreach (double cs in _cs)
            {                
                writer.WriteLine(cs);
            }
            writer.Flush();
            writer.Close();
            fileStream.Close();
        } // WriteCsHistogram()

        public void WritePCHistogram(string fileName)
        {
            FileStream fileStream = new FileStream(fileName, FileMode.Create, FileAccess.Write);
            StreamWriter writer = new StreamWriter(fileStream);          
            foreach (double pc in _pC)
            {
                writer.WriteLine(pc);
            }
            writer.Flush();
            writer.Close();
            fileStream.Close();
        } // WritePCHistogram()

        public void WriteDebugDetails(string fileName, FileMode fileMode, FileAccess fileAccess)
        {
            FileStream fileStream = new FileStream(fileName, fileMode, fileAccess);
            StreamWriter writer = new StreamWriter(fileStream);
            writer.WriteLine("-- prior model --");
            // n nf nfe sumPC totalL averageL sumPL, Up
            writer.WriteLine(String.Format("{0,-6}{1,-6}{2,-6}{3, -12}{4, -12}{5, -12}{6, -12}{7,-12}",
                "n", "nf", "nfe", "sumPC", "total L", "ave L", "sumPL", "Up"));
            writer.WriteLine(String.Format("{0,-6:d}{1,-6:d}{2,-6:d}{3,-12:0.00e+000}{4, -12:d}{5, -12:0.00e+000}{6, -12:0.00e+000}{7, -12:0.00e+000}",
                _n, _nf, _nfe, _sumPC, (int)_totalL, _averageL, _sumPL, _Up));
            writer.WriteLine("-- data term --");
            // sumPG sumPI, Ud, U, h
            writer.WriteLine(String.Format("{0,-12}{1,-12}{2,-12}{3, 20}{4, 12}",
                "sumPG", "sumPI", "Ud", "U", "h"));
            writer.WriteLine(String.Format("{0,-12:0.00e+000}{1,-12:0.00e+000}{2,-12:0.00e+000}{3, 20:0.00e+000}{4, 12:0.00e+000}",
                _sumPG, _sumPI, _Ud, _U, _H));
            // line segment details
            writer.WriteLine("-- Line Segment Details --");
            writer.WriteLine(String.Format("{0,-6}{1,-6}{2,-12}{3,-15}{4,-12}{5,-12}{6,-12}{7,-12}",
                "Index", "Ls", "pL", "Cs", "Is", "pI", "Gs", "pG"));
            // for each line segment
            for (int s = 0; s < _n; s++)
            {
                writer.WriteLine(String.Format("{0,-6:d}{1,-6:d}{2,-12:0.00e+000}{3,-15}{4,-12:0.00e+000}{5,-12:0.00e+000}{6,-12:0.00e+000}{7,-12:0.00e+000}",
                    _configuration[s].Index, _configuration[s].Ls, _configuration[s].PL, _configuration[s].Cs.ToString(), _configuration[s].Is,
                    _configuration[s].PI, _configuration[s].Gs, _configuration[s].PG));
            }
            writer.WriteLine();
            writer.Flush();
            writer.Close();
            fileStream.Close();
        } // WriteDebugDetails()
        public void WriteDebugDetails(string fileName)
        {
            WriteDebugDetails(fileName, FileMode.OpenOrCreate, FileAccess.Write);
        } // WriteDebugDetails()

        public void WriteAllDetails(string fileName, FileMode fileMode, FileAccess fileAccess)
        {
            FileStream fileStream = new FileStream(fileName, fileMode, fileAccess);
            StreamWriter writer = new StreamWriter(fileStream);
            writer.WriteLine("-- prior model --");
            // n nf nfe sumPC totalL averageL sumPL, Up
            writer.WriteLine(String.Format("{0,-6}{1,-6}{2,-6}{3, -12}{4, -12}{5, -12}{6, -12}{7,-12}",
                "n", "nf", "nfe", "sumPC", "total L", "ave L", "sumPL", "Up"));
            writer.WriteLine(String.Format("{0,-6:d}{1,-6:d}{2,-6:d}{3,-12:0.00e+000}{4, -12:d}{5, -12:0.00e+000}{6, -12:0.00e+000}{7, -12:0.00e+000}",
                _n, _nf, _nfe, _sumPC, (int)_totalL, _averageL, _sumPL, _Up));
            writer.WriteLine("-- data term --");
            // sumPG sumPI, Ud, U, h
            writer.WriteLine(String.Format("{0,-12}{1,-12}{2,-12}{3, 20}{4, 12}",
                "sumPG", "sumPI", "Ud", "U", "h"));
            writer.WriteLine(String.Format("{0,-12:0.00e+000}{1,-12:0.00e+000}{2,-12:0.00e+000}{3, 20:0.00e+000}{4, 12:0.00e+000}",
                _sumPG, _sumPI, _Ud, _U, _H));
            // line segment details
            writer.WriteLine("-- Line Segment Details --");
            writer.WriteLine(String.Format("{0,-6}{1,-6}{2,-12}{3,-15}{4,-12}{5,-12}{6,-12}{7,-12}",
                "Index", "Ls", "pL", "Cs", "Is", "pI", "Gs", "pG"));
            // for each line segment
            for (int s = 0; s < _n; s++)
            {
                writer.WriteLine(String.Format("{0,-6:d}{1,-6:d}{2,-12:0.00e+000}{3,-15}{4,-12:0.00e+000}{5,-12:0.00e+000}{6,-12:0.00e+000}{7,-12:0.00e+000}",
                    _configuration[s].Index, _configuration[s].Ls, _configuration[s].PL, _configuration[s].Cs.ToString(), _configuration[s].Is,
                    _configuration[s].PI, _configuration[s].Gs, _configuration[s].PG));
            }
            // si, sj, Csi, Csj, pC(si, sj)
            writer.WriteLine(String.Format("{0,-5}{1,-5}{2,-15}{3,-15}{4,-12}",
                "si", "sj", "Csi", "Csj", "pC(si, sj)"));
            for (int i = 0; i < _n; i++)
            {
                for (int j = i; j < _n; j++)
                {
                    if (i == j)
                        continue;
                    if (_pCPair[i, j] != 0.0)
                        writer.WriteLine(String.Format("{0,-5:d}{1,-5:d}{2,-15}{3,-15}{4,-12:0.00e+000}",
                            i, j, _configuration[i].Cs.ToString(), _configuration[j].Cs.ToString(), _pCPair[i, j]));
                }
            }
            writer.WriteLine();
            writer.Flush();
            writer.Close();
            fileStream.Close();
        } // WriteAllDetails()
        public void WriteAllDetails(string fileName)
        {
            WriteAllDetails(fileName, FileMode.OpenOrCreate, FileAccess.Write);
        } // WriteAllDetails()

        /// <summary>
        /// convert line segment configuration to a matrix
        /// white background, black line segments
        /// </summary>
        /// <param name="srcImage"></param>
        /// <returns></returns>
        public PZMath_matrix ConvertToMatrix(Bitmap srcImage)
        {
            int width = srcImage.Width;
            int height = srcImage.Height;
            PZMath_matrix dstMatrix = new PZMath_matrix(height, width);
            dstMatrix.Setall(255);
            for (int i = 0; i < _n; i++)
            {
                int length = _configuration[i].Ls;
                for (int j = 0; j < length; j++)
                {
                    int x = (int)_configuration[i].PointList[j].x;
                    int y = (int)_configuration[i].PointList[j].y;
                    dstMatrix[y, x] = 0;
                }
            }

            return dstMatrix;
        } // PZMath_matrix()

        #endregion

        #region synthetic model study
        public void SynthethicModelStudy(string fileName)
        {
            // start from full configuration
            int fullLength = this.N;
            int[] input = new int[fullLength];
            for (int i = 0; i < fullLength; i++)
                input[i] = i;

            bool first = true;
            // for each possible length combination
            for (int l = 1; l <= fullLength; l++)
            {
                Combinations<int> combinations = new Combinations<int>(input, l);
                foreach (int[] combination in combinations)
                {
                    // copy of the full configuration
                    LineSegmentConfiguration c = new LineSegmentConfiguration(this);

                    // find the complement set
                    List<int> complement = new List<int>(input);
                    for (int i = 0; i < l; i++)
                        complement.Remove(combination[i]);

                    // get the configuration according the combination
                    for (int i = 0; i < fullLength - l; i++)
                        c.RemoveAtIndex(complement[i]);

                    // calculate Up, Ud, U an H of the current configuration
                    c.CalculateUP();
                    c.CalculateUD();
                    c.CalculateU();
                    c.CalculateH();

                    // output line segment details
                    if (first)
                    {
                        first = false;
                        c.WriteDebugDetails(fileName);
                    }
                    else
                    {
                        c.WriteDebugDetails(fileName, FileMode.Append, FileAccess.Write);
                    }
                }
            }
        } // SynthethicModelStudy()
        #endregion

        #region Example
        #endregion
    }
}
