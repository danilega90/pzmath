// pz 12.06.2006 line segment

using System;
using System.Collections.Generic;
using System.Text;
using eee.Sheffield.PZ.Math;
using System.Drawing;
using System.Drawing.Imaging;
using System.Runtime.Serialization.Formatters.Binary;
using System.Runtime.Serialization;
using System.IO;

namespace eee.Sheffield.PZ.Imaging
{
    public enum ConnectionPointType { Start, End };

    /// <summary>
    /// line segment
    /// </summary>
    [Serializable]
    public class LineSegment : ICloneable
    {
        #region Fields
        // basic
        private int _index = 0;
        private PZPoint _startPoint = null;
        private PZPoint _endPoint = null;
        private List<PZPoint> _pointList = null;
        private List<PZDirection> _directionList = null;

        // prior model
        private PZDirection _Cs = null; // average line segment angle
        private int _Ls = 0;    // length
        private double _pL = 0.0;        // pL(s)

        // data term
        private List<int> _intensityList = null;
        private List<double> _gvfUList = null;
        private List<double> _gvfVList = null;
        private List<double> _gvfMagnitudeList = null;

        private double _Gs = 0.0; // Gs
        public double _angle = 0.0; // angle
        public List<double> _GsList;   // Gs List
        public List<double> _angleList;    // angle List

        private double _pG = 0.0; // pG(s)
        private double _Is = 0.0; // average intensity
        private double _pI = 0.0; // pI(s)      

        // MCMC model
        private int _hit = 0;
        private double _hitProbability = 0.0;

        // connection 
        private List<LineSegment> _startPointConnectionList = new List<LineSegment>(0);
        private List<ConnectionPointType> _startPointConnectionTypeList = new List<ConnectionPointType>(0);
        private List<LineSegment> _endPointConnectionList = new List<LineSegment>(0);
        private List<ConnectionPointType> _endPointConnectionTypeList = new List<ConnectionPointType>(0);
        #endregion

        #region Properties
        public int Index { get { return _index; } set { _index = value; } }
        public int Hit { get { return _hit; } set { _hit = value; } }
        public double HitProbability { get { return _hitProbability; } set { _hitProbability = value; } }
        public PZPoint EndPoint { get { return _endPoint; } set { _endPoint = value; } }
        public PZPoint StartPoint { get { return _startPoint; } set { _startPoint = value; } }
        public List<PZPoint> PointList { get { return _pointList; } set { _pointList = value; } }
        public int Ls { get { return _Ls; } set { _Ls = value; } }

        public List<PZDirection> DirectionList { get { return _directionList; } }

        public List<int> IntensityList { get { return _intensityList; } }
        public List<double> GVFUList { get { return _gvfUList; } }
        public List<double> GVFVList { get { return _gvfVList; } }
        public List<double> GVFMagnitudeList { get { return _gvfMagnitudeList; } }

        public PZDirection Cs { set { _Cs = value; } get { return _Cs; } }
        public double Gs { set { _Gs = value; } get { return _Gs; } }
        public double PG { set { _pG = value; } get { return _pG; } }
        public double Is { set { _Is = value; } get { return _Is; } }
        public double PI { set { _pI = value; } get { return _pI; } }
        public double PL { set { _pL = value; } get { return _pL; } }

        public List<LineSegment> StartPointConnectionList { get { return _startPointConnectionList; } }
        public List<LineSegment> EndPointConnectionList { get { return _endPointConnectionList; } }
        public List<ConnectionPointType> StartPointConnectionTypeList { get { return _startPointConnectionTypeList; } }
        public List<ConnectionPointType> EndPointConnectionTypeList { get { return _endPointConnectionTypeList; } }
        #endregion

        #region constructors
        /// <summary>
        /// empty constructor
        /// </summary>
        public LineSegment()
        {
            _pointList = new List<PZPoint>(0);
            _directionList = new List<PZDirection>(0);

            // data term
            _intensityList = new List<int>(0);
            _gvfUList = new List<double>(0);
            _gvfVList = new List<double>(0);
            _gvfMagnitudeList = new List<double>(0);
        } // LineSegment()
        #endregion

        #region Add points
        /// <summary>
        /// add a start/end point
        /// </summary>
        /// <param name="point"></param>
        public void AddEndPoint(PZPoint point)
        {
            if (_startPoint == null)
            {
                _startPoint = point;
                _pointList.Add(point);
                _Ls++;
            }
            else if (_endPoint == null)
            {
                _endPoint = point;
                _pointList.Add(point);
                _Ls++;
            }
            else
                PZMath_errno.ERROR("Illigal Start/End Points");
        } // AddEndPoint()

        /// <summary>
        /// add a point to the point list
        /// </summary>
        /// <param name="point"></param>
        /// <returns></returns>
        public int AddPoint(PZPoint point)
        {
            _pointList.Add(point);
            _Ls++;
            return _Ls;
        } // AddPoint()
        #endregion
       
        #region valid check
        /// <summary>
        /// check line segment connectivity
        /// </summary>
        /// <returns></returns>
        public bool CheckLineSegmentConnectivity(out string errorReport)
        {
            // check two end points
            if (_startPoint == null || _endPoint == null)
            {
                errorReport = "Line segment is not finished";
                return false;
            }
            // check connectivity
            int lengthM1 = _Ls - 1;
            for (int x = 0; x < lengthM1; x++)
            {
                if (!_pointList[x].Is8NeighbourOf(_pointList[x + 1]))
                {
                    errorReport = "Line segment is not connected.";
                    return false;
                }
            }
            errorReport = "pass check.";
            return true;
        } // CheckLineSegmentConnectivity()
        #endregion

        #region Line Segment relationship method
        /// <summary>
        /// merge the current line segment to the other, end to end
        /// only merge point list, so that the merged line segment needs Initialize()
        /// connections are copied. 
        /// </summary>
        /// <returns></returns>
        public LineSegment MergeTo(LineSegment l)
        {
            LineSegment newl = new LineSegment();
            // this is a, and input (l) is b
            int aLength = this.Ls;
            int bLength = l.Ls;
            int newLength = aLength + bLength - 1;
            int[] pointIndex = new int[newLength];
            #region find start and end points
            // find start and end point
            PZPoint aStart = _startPoint;
            PZPoint aEnd = _endPoint;
            PZPoint bStart = l.StartPoint;
            PZPoint bEnd = l.EndPoint;
            if (aStart.Equals(bStart))
            {
                // aEnd -> aStart(bStart) -> bEnd
                for (int i = 0; i < aLength; i++)
                    pointIndex[i] = aLength - i - 1;
                for (int i = 1; i < bLength; i++)
                    pointIndex[i + aLength - 1] = i;
            }
            else if (aStart.Equals(bEnd))
            {
                // aEnd -> aStart(bEnd) -> bStart
                for (int i = 0; i < aLength; i++)
                    pointIndex[i] = aLength - i - 1;
                for (int i = 1; i < bLength; i++)
                    pointIndex[i + aLength - 1] = bLength - i - 1;
            }
            else if (aEnd.Equals(bStart))
            {
                // aStart -> aEnd(bStart) -> bEnd
                for (int i = 0; i < aLength; i++)
                    pointIndex[i] = i;
                for (int i = 1; i < bLength; i++)
                    pointIndex[i + aLength - 1] = i;
            }
            else
            {
                // aStart- -> aEnd(bEnd) -> bStart
                for (int i = 0; i < aLength; i++)
                    pointIndex[i] = i;
                for (int i = 1; i < bLength; i++)
                    pointIndex[i + aLength - 1] = bLength - i - 1;
            }
            #endregion
            // generate the new line segment            
            newl.Ls = newLength;
            // a -> b
            newl.StartPoint = new PZPoint(_pointList[pointIndex[0]]);
            newl.EndPoint = new PZPoint(l.PointList[pointIndex[newLength - 1]]);
            int indexCount = 0;
            for (int i = 0; i < aLength; i++)
            {
                newl.PointList.Add(_pointList[pointIndex[indexCount]]);
                indexCount++;
            }
            for (int i = 1; i < bLength; i++)
            {
                newl.PointList.Add(l.PointList[pointIndex[indexCount]]);
                indexCount++;
            }
            return newl;
        } // MergeTo()

        /// <summary>
        /// is the line segment connected with l?
        /// </summary>
        /// <param name="l"></param>
        /// <returns></returns>
        public bool IsConnectedWith(LineSegment l)
        {
            bool isConnected = true;
            if ((_startPoint.EqualTo(l._startPoint) && _endPoint.EqualTo(l._endPoint))
                || (_startPoint.EqualTo(l._endPoint) && _endPoint.EqualTo(l._startPoint)))
                // a loop
                isConnected = false;
            else
            {
                if (_startPoint.EqualTo(l._startPoint))
                {
                    _startPointConnectionList.Add(l);
                    _startPointConnectionTypeList.Add(ConnectionPointType.Start);
                }
                else if (_startPoint.EqualTo(l._endPoint))
                {
                    _startPointConnectionList.Add(l);
                    _startPointConnectionTypeList.Add(ConnectionPointType.End);
                }
                else if (_endPoint.EqualTo(l._startPoint))
                {
                    _endPointConnectionList.Add(l);
                    _endPointConnectionTypeList.Add(ConnectionPointType.Start);
                }
                else if (_endPoint.EqualTo(l._endPoint))
                {
                    _endPointConnectionList.Add(l);
                    _endPointConnectionTypeList.Add(ConnectionPointType.End);
                }
                else
                    isConnected = false;
            }
            return isConnected;
        } // IsConnectedWith()

        /// <summary>
        /// is the line segment part of l?
        /// </summary>
        /// <param name="l"></param>
        /// <returns></returns>
        public bool IsPartOf(LineSegment l)
        {
            if (_pointList.Count > l.PointList.Count)
                return false;


            foreach (PZPoint p in _pointList)
            {
                bool pIsPart = false;

                foreach (PZPoint q in l.PointList)
                {
                    if (p.EqualTo(q))
                    {
                        pIsPart = true;
                        break;
                    }
                }

                if (!pIsPart)
                {
                    return false;
                }                
            }

            return true;
        } // IsPartOf()
        #endregion

        #region prior model
        /// <summary>
        /// calculate each point direction, central difference
        /// </summary>
        /// <param name="pList"></param>
        /// <param name="dList"></param>
        private void CalculateDirection(List<PZPoint> pList, ref List<PZDirection> dList)
        {
            int length = pList.Count;
            int lengthM1 = length - 1;
            if (dList.Count > 0)
                dList.Clear();
            for (int i = 0; i < length; i++)
            {
                double dx, dy;
                if (i == 0)
                {
                    PZPoint p0 = pList[0];
                    PZPoint p1 = pList[1];
                    dx = p1.x - p0.x;
                    dy = p1.y - p0.y;

                }
                else if (i > 0 && i < lengthM1)
                {
                    PZPoint pM1 = pList[i - 1];
                    PZPoint pA1 = pList[i + 1];
                    dx = (pA1.x - pM1.x) / 2.0;
                    dy = (pA1.y - pM1.y) / 2.0;
                }
                else
                {
                    PZPoint pM1 = pList[i - 1];
                    PZPoint p = pList[i];
                    dx = p.x - pM1.x;
                    dy = p.y - pM1.y;
                }
                dList.Add(new PZDirection(dx, dy));
            }
        } // CalculateDirection()

        /// <summary>
        /// calculate Cs, the average line segment direction
        /// </summary>
        /// <param name="directionList"></param>
        /// <returns></returns>
        private PZDirection CalculateCs(List<PZDirection> directionList)
        {
            double x = 0.0;
            double y = 0.0;
            int length = directionList.Count;
            for (int i = 0; i < length; i++)
            {
                x += directionList[i].x;
                y += directionList[i].y;
            }
            PZDirection thetas = new PZDirection(x / (double)length, y / (double)length);
            return thetas;
        } // CalculateCs()
        public void CalculateCs()
        {
            _Cs = CalculateCs(_directionList);
        } // CalculateCs()

        /// <summary>
        /// calculate pL(s)
        /// </summary>
        /// <param name="pointList">point list</param>
        /// <param name="Lt">total length of the line segment configuration</param>
        /// <param name="Laverage">average length</param>
        /// <returns></returns>
        private double CalculatePL(List<PZPoint> pointList, double Lt, double Laverage)
        {
            double Ls = Convert.ToDouble(pointList.Count);
            double pL = 0.0;

            if (Lt == Laverage)
                pL = -1.0;
            else
            {
                if (Ls < Laverage)
                    pL = 1.0;
                else
                    pL = Sigma(Ls - Laverage, Lt - Laverage) - 1.0;
            }
            return pL;
        } // CalculatePL()
        /// <summary>
        /// calculate pL(s)
        /// </summary>
        /// <param name="Lt">total length of the line segment configuration</param>
        /// <param name="Laverage">average length</param>
        public void CalculatePL(double Lt, double Laverage)
        {
            _pL = CalculatePL(_pointList, Lt, Laverage);
        } // CalculatePL()
        #endregion

        #region data term
        /// <summary>
        /// assign data term
        /// </summary>
        private void AssignDataTerm
            (Bitmap srcImage, PZMath_matrix gvfU, PZMath_matrix gvfV, PZMath_matrix gvfMagnitude,
            List<PZPoint> pointList, List<PZDirection> directionList,
            ref List<int> intensityList, ref List<double> gvfUList, ref List<double> gvfVList, ref List<double> gvfMagnitudeList)
        {
            // clear lists
            intensityList.Clear();
            gvfUList.Clear();
            gvfVList.Clear();
            gvfMagnitudeList.Clear();

            int length = pointList.Count;
            double offset = 1.0;
            // for each point
            for (int i = 0; i < length; i++)
            {
                PZPoint point = pointList[i];
                double x = point.x;
                double y = point.y;
                // assign intensity
                int intensity = Convert.ToInt32(srcImage.GetPixel((int)x, (int)y).R);
                intensityList.Add(intensity);

                // find normal direction
                PZDirection lineDirection = new PZDirection(directionList[i]);
                lineDirection.Rotate(90);
                // offset 1 or 2 pixels
                y += lineDirection.y * offset;
                x += lineDirection.x * offset;

                // asigne GVF info
                gvfUList.Add(gvfU[(int)y, (int)x]);
                gvfVList.Add(gvfV[(int)y, (int)x]);
                gvfMagnitudeList.Add(gvfMagnitude[(int)y, (int)x]);
            }
        } // AssignDataTerm()

        /// <summary>
        /// calculate Gs
        /// </summary>
        /// <param name="directionList"></param>
        /// <param name="gvfUList"></param>
        /// <param name="gvfVList"></param>
        /// <param name="gvfMagnitudeList"></param>
        /// <returns></returns>
        private double CalculateGs(List<PZDirection> directionList,
            List<double> gvfUList, List<double> gvfVList, List<double> gvfMagnitudeList)
        {
            #region debug
            //string fileName = "centerline.txt";
            //FileStream fs = new FileStream(fileName, FileMode.Create, FileAccess.Write);
            //StreamWriter w = new StreamWriter(fs);
            //w.WriteLine(String.Format(
            //    "{0,20}{1,20}{2,20}{3,20}{4,20}{5,20}{6,20}", "x", "y", "gix", "giy", "uix", "uiy", "guis"));
            #endregion

            int length = directionList.Count;
            
            _GsList = new List<double>(length);
            _angleList = new List<double>(length);

            double Gs = 0.0;
            // for each point
            for (int i = 0; i < length; i++)
            {
                // calculate projection guis (gi project to uis)
                // guis = <gi, uis>, (uis is unit vector, <,> is inner product)
                double guis = 0;
                PZDirection uis = directionList[i];
                double gix = gvfUList[i];
                double giy = gvfVList[i];
                double giMagnitude = gvfMagnitudeList[i];
                //gix *= giMagnitude;
                //giy *= giMagnitude;
                guis = System.Math.Abs(gix * uis.x + giy * uis.y);

                _GsList.Add(guis);
                _angleList.Add(System.Math.Acos(guis));

                #region debug
                //w.WriteLine(String.Format(
                //    "{0,20:0.0000}{1,20:0.0000}{2,20:0.0000}{3,20:0.0000}{4,20:0.0000}{5,20:0.0000}{6,20:0.0000}", 
                //    l.PointList[i].x, l.PointList[i].y, gix, giy, uis.x, uis.y, guis));
                #endregion

                Gs += guis;
            }
            Gs /= (double)length;

            _angle = System.Math.Acos(Gs);

            #region debug
            //w.Flush();
            //w.Close();
            //fs.Close();
            #endregion

            return Gs;
        } // CalculateGs()
        public void CalculateGs()
        {
            _Gs = CalculateGs(_directionList, _gvfUList, _gvfVList, _gvfMagnitudeList);
        } // CalculateGs()

        /// <summary>
        /// calculate pG(s)
        /// </summary>
        /// <param name="Gs"></param>
        /// <param name="Gt"></param>
        /// <returns></returns>
        private double CalculatePG(double Gs, double Gt)
        {
            double pG = 1.0;
            if (Gs > Gt)
                pG = 2.0;
            else
                pG = -1.0 * Sigma(Gs, Gt);
            return pG;
        } // CalculatePG()
        public void CalculatePG(double Gt)
        {
            _pG = CalculatePG(_Gs, Gt);
        } // CalculatePG()

        /// <summary>
        /// calculate Is
        /// </summary>
        /// <returns></returns>
        private double CalculateIs(List<int> intensityList)
        {
            int length = intensityList.Count;
            double Is = 0.0;
            // for each point 
            for (int i = 0; i < length; i++)
            {
                Is += (double)intensityList[i];
            }
            Is /= (double)length;
            return Is;
        } // CalculateIs()
        public void CalculateIs()
        {
            _Is = CalculateIs(_intensityList);
        } // CalculateIs()

        /// <summary>
        /// calculate pI(s)
        /// </summary>
        /// <returns></returns>
        private double CalculatePI(double Is, double It, double kI)
        {
            double threshold = It / kI;
            double pI = 0.0;
            if (Is > threshold)
                pI = 1.0;
            else
                pI = -1.0 * Sigma(Is, threshold);
            return pI;
        } // CalculatePI()
        public void CalculatePI(double It, double kI)
        {
            _pI = CalculatePI(_Is, It, kI);
        } // CalculatePI()
        #endregion

        ///initialize after line segment is obtained.
        #region initialize method
        /// <summary>
        /// initialize prior model
        /// </summary>
        public void InitializePriorModel()
        {
            // check line segment validity
            //string errorMessage;
            //if (!CheckLineSegmentConnectivity(out errorMessage))
            //    throw new ApplicationException("LineSegment::InitializePirorModel(), " + errorMessage);

            // clear connections
            //_startPointConnectionList.Clear();
            //_startPointConnectionTypeList.Clear();
            //_endPointConnectionList.Clear();
            //_endPointConnectionTypeList.Clear();

            // calculate direction list
            CalculateDirection(_pointList, ref _directionList);
            // calculate Cs
            _Cs = CalculateCs(_directionList);
            // assure Ls
            _Ls = _pointList.Count;
        } // InitializePriorModel()

        /// <summary>
        /// initialize data term
        /// </summary>
        /// <param name="srcImage"></param>
        /// <param name="gvfU"></param>
        /// <param name="gvfV"></param>
        /// <param name="gvfMagnitude"></param>
        public void InitializeDataTerm(Bitmap srcImage, PZMath_matrix gvfU, PZMath_matrix gvfV, PZMath_matrix gvfMagnitude)
        {
            // assign data term
            AssignDataTerm(srcImage, gvfU, gvfV, gvfMagnitude, _pointList, _directionList, ref _intensityList, ref _gvfUList, ref _gvfVList, ref _gvfMagnitudeList);
            // calculate Is
            _Is = CalculateIs(_intensityList);
            // calculate Gs
            _Gs = CalculateGs(_directionList, _gvfUList, _gvfVList, _gvfMagnitudeList);
        } // InitializeDataTerm()

        /// <summary>
        /// initialize line segment, Cs, Ls, Is and Gs are then calculated
        /// </summary>
        /// <param name="srcImage"></param>
        /// <param name="gvfU"></param>
        /// <param name="gvfV"></param>
        /// <param name="gvfMagnitude"></param>
        public void Initialize(Bitmap srcImage, PZMath_matrix gvfU, PZMath_matrix gvfV, PZMath_matrix gvfMagnitude)
        {
            InitializePriorModel();
            InitializeDataTerm(srcImage, gvfU, gvfV, gvfMagnitude);
        } // Initialize()
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
        #endregion

        public object Clone(bool doDeepCopy)
        {
            if (doDeepCopy)
            {
                BinaryFormatter BF = new BinaryFormatter();
                MemoryStream memStream = new MemoryStream();

                BF.Serialize(memStream, this);
                memStream.Position = 0;

                return (BF.Deserialize(memStream));
            }
            else
            {
                return (this.MemberwiseClone());
            }
        }

        public object Clone()
        {
            return (Clone(false));
        }
        #region I/O
        #endregion
    }
}
