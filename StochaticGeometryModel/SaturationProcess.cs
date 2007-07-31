// eee.Sheffield.PZ.Math
//
// Copyright ?Ping Zou, 2007
// sg71.cherub@gmail.com

using System;
using System.Collections.Generic;
using System.Text;
using System.Drawing;
using eee.Sheffield.PZ.Imaging;
using eee.Sheffield.PZ.Math;
using System.IO;


namespace eee.Sheffield.PZ.Math
{
    /// <summary>
    /// Saturation process
    /// f_theta(x) = 1 / c(theta) * exp(<t(x), theta>)
    /// t(x) = (n(x), u(x))
    /// u(x) = sum_xi(min(c, m_xi(x))
    /// m_xi(x) = sum_eta(1[0, r](d(xi, eta))
    /// theta = (theta1, theta2)
    /// <t(x), theta> = n(x) * theta1 + u(x) * theta2
    /// theta is equal or less than 0
    /// </summary>
    public class SaturationProcess
    {
        #region Fields
        // Strauss process parameters
        private double _theta1;
        private double _theta2;
        private double _r;          // radius of interaction region
        private double _c;
        private double _lamda;
        private double _width;
        private double _height;

        // MCMC parameters
        private int _total = 10000;         // total # of samples
        private int _burnIn = 1000;        // burn-in #
        private int _space = 100;         // samples are recorded once very "space" samples

        // MCMC records
        private List<PZPoint[]> _sampleList = new List<PZPoint[]>(0);    // recorded samples, each sample is a PZPoint array
        private List<double> _birthAcceptRateList = new List<double>(0);    // birth accept rate
        private List<double> _deathAcceptRateList = new List<double>(0);    // death accept rate
        private List<int> _numberPointsList = new List<int>(0);  // number of points, n(x)
        private List<int> _numberClosePointsList = new List<int>(0);    // number of close point, u(x)
        #endregion

        #region Properties
        public double Theta1 { get { return _theta1; } set { _theta1 = value; } }
        public double Theta2 { get { return _theta2; } set { _theta2 = value; } }
        public double R { get { return _r; } set { _r = value; } }
        public double C { get { return _c; } set { _c = value; } }
        public double Lamda { get { return _lamda; } set { _lamda = value; } }
        public double Width { get { return _width; } set { _width = value; } }
        public double Height { get { return _height; } set { _height = value; } }
        public int Total { get { return _total; } set { _total = value; } }
        public int BurnIn { get { return _burnIn; } set { _burnIn = value; } }
        public int Space { get { return _space; } set { _space = value; } }

        public List<double> BirthAcceptRateList { get { return _birthAcceptRateList; } }
        public List<double> DeathAcceptRateList { get { return _deathAcceptRateList; } }
        public List<int> NumberPointsList { get { return _numberPointsList; } }
        public List<int> NumberClosePointsList { get { return _numberClosePointsList; } }
        #endregion

        #region Constructor
        public SaturationProcess(double theta1, double theta2, double r, double c, double lamda, double width, double height)
        {
            _theta1 = theta1;
            _theta2 = theta2;
            _r = r;
            _c = c;
            _lamda = lamda;
            _width = width;
            _height = height;
        } // StraussProcess()
        #endregion

        #region util methods
        /// <summary>
        /// the number of pairs of points in point array that are separated by distance no more than r
        /// </summary>
        /// <param name="pointArray"></param>
        /// <returns></returns>
        private int S(PZPoint[] pointArray, double r)
        {
            int count = pointArray.Length;
            int s = 0;
            for (int i = 0; i < count - 1; i++)
            {
                for (int j = i + 1; j < count; j++)
                {
                    if (pointArray[i].Distance(pointArray[j]) <= r)
                        s++;

                    //PZPoint xi = pointArray[i];
                    //PZPoint eta = pointArray[j];
                    //if (xi.Distance(eta) <= r)
                    //    s++;
                }
            }
            return s;

            //return s / 2;
        } // S()

        /// <summary>
        /// calculate value of unnormalized density
        /// </summary>
        /// <param name="pointArray"></param>
        /// <param name="theta1"></param>
        /// <param name="theta2"></param>
        /// <param name="r"></param>
        /// <returns></returns>
        private double H(PZPoint[] pointArray, double theta1, double theta2, double r, double c)
        {
            int nx = pointArray.Length;
            double mx = 2.0 * (double)S(pointArray, r);
            double ux = mx < c ? mx : c;
            double h = System.Math.Exp(nx * theta1 + ux * theta2);
            return h;
        } // H()

        public double H(PZPoint[] pointArray)
        {
            return H(pointArray, _theta1, _theta2, _r, _c);
        } // H()
        #endregion

        #region stochatic process method
        public void EvolveStraussProcess()
        {
            EvolveStraussProcessStartFromEmpty(_sampleList,
                _theta1, _theta2, _r, _c,
                _lamda, _width, _height,
                _total, _burnIn, _space,
                _birthAcceptRateList, _deathAcceptRateList,
                _numberPointsList, _numberClosePointsList);
        } // EvolveStraussProcess()

        /// <summary>
        /// evolve Strauss process, C. J. Geyer 1999
        /// start from empty point array
        /// </summary>
        private void EvolveStraussProcessStartFromEmpty(List<PZPoint[]> sampleList,
            double theta1, double theta2, double r, double c,
            double lamda, double width, double height,
            int total, int burIn, int space,
            List<double> birthAcceptRateList, List<double> deathAcceptRateList,
            List<int> numberPointsList, List<int> numberClosePointsList)
        {
            ParkMillerUniform uniformRandom = new ParkMillerUniform();
            PZPoint[] oldPointArray = new PZPoint[0];
            PZPoint[] newPointArray;
            double birthRate = 0;
            double deathRate = 0;
            double birthAcceptRate = 0;
            double deathAcceptRate = 0;
            bool acceptNewPointArray = false;

            int spaceCount = 0;

            // MCMC evolving process
            for (int i = 0; i < _total; i++)
            {
                // birth or death
                if (oldPointArray.Length == 0)
                // birth only
                {
                    birthRate = 1.0;
                    deathRate = 0.0;
                }
                else
                {
                    birthRate = uniformRandom.Sample();
                    deathRate = 1 - birthRate;
                }

                int nx = oldPointArray.Length;
                int sx = S(oldPointArray, r);

                if (birthRate >= deathRate)
                // birth
                {
                    // simulate xi distributed proportional to lamda
                    PZPoint xi = new PZPoint(width * uniformRandom.Sample(), height * uniformRandom.Sample());
                    newPointArray = new PZPoint[nx + 1];
                    // get new point array
                    Array.Copy(oldPointArray, newPointArray, nx);
                    newPointArray[nx] = xi;
                    // birth accept rate
                    birthAcceptRate = BirthAcceptRate(oldPointArray, newPointArray, lamda, width, height, theta1, theta2, r, c);
                    birthAcceptRate = birthAcceptRate < 1.0 ? birthAcceptRate : 1.0;
                    acceptNewPointArray = uniformRandom.Sample() < birthAcceptRate;
                }
                else
                // death
                {
                    int removeIndex = (int)System.Math.Floor(nx * uniformRandom.Sample());
                    int nxNew = nx - 1;
                    newPointArray = new PZPoint[nxNew];
                    // get new point array
                    if (nx - 1 > 0)
                    {
                        for (int indexNew = 0, indexOld = 0; indexNew < nxNew; indexNew++, indexOld++)
                        {
                            if (indexOld != removeIndex)
                            {
                                newPointArray[indexNew] = oldPointArray[indexOld];
                            }
                            else
                            {
                                indexNew--;
                            }
                        }
                    }
                    // death accept rate
                    deathAcceptRate = DeathAcceptRate(oldPointArray, newPointArray, lamda, width, height, theta1, theta2, r, c);
                    deathAcceptRate = deathAcceptRate < 1.0 ? deathAcceptRate : 1.0;
                    acceptNewPointArray = uniformRandom.Sample() < deathAcceptRate;
                }

                if (acceptNewPointArray)
                // accept new point array
                {
                    oldPointArray = newPointArray;
                }
                else
                // reject new point array
                {
                }

                // record MCMC                
                if (i > burIn)
                {
                    numberPointsList.Add(nx);
                    numberClosePointsList.Add(sx);
                    if (birthRate > deathRate)
                        birthAcceptRateList.Add(birthAcceptRate);
                    else
                        deathAcceptRateList.Add(deathAcceptRate);

                    if (spaceCount == space)
                    {
                        spaceCount = 0;
                        PZPoint[] recordPointArray = new PZPoint[oldPointArray.Length];
                        Array.Copy(oldPointArray, recordPointArray, oldPointArray.Length);
                        sampleList.Add(recordPointArray);
                    }
                }

                spaceCount++;
            }
        } // EvolveStraussProcessStartFromEmpty()

        /// <summary>
        /// birth accept rate, C. J. Geyer
        /// R = lamda(S) * h(x U xi) / (n(x) + 1) / h(x)
        /// </summary>
        private double BirthAcceptRate(PZPoint[] oldPointArray, PZPoint[] newPointArray,
            double lamda, double width, double height,
            double theta1, double theta2, double r, double c)
        {
            double lamdaS = lamda * width * height;
            double nx = oldPointArray.Length;
            double hxOld = H(oldPointArray, theta1, theta2, r, c);
            double hxNew = H(newPointArray, theta1, theta2, r, c);
            double R = lamdaS * hxNew / (nx + 1.0) / hxOld;
            return R;
        } // BirthAcceptRate()

        /// <summary>
        /// death accept rate, C. J. Geyer
        /// R = n(x) * h(x \ xi) / lamda(s) / h(x)
        /// </summary>
        private double DeathAcceptRate(PZPoint[] oldPointArray, PZPoint[] newPointArray,
            double lamda, double width, double height,
            double theta1, double theta2, double r, double c)
        {
            double lamdaS = lamda * width * height;
            double nx = oldPointArray.Length;
            double hxOld = H(oldPointArray, theta1, theta2, r, c);
            double hxNew = H(newPointArray, theta1, theta2, r, c);
            double R = nx * hxNew / lamdaS / hxOld;
            return R;
        } // DeathAcceptRate()
        #endregion

        #region Example
        public static void Example()
        {
            string fileName;
            FileStream fileStream;
            StreamWriter streamWriter;

            double r = 0.05;
            double theta1 = 4.0;
            double theta2 = 0.4;
            double c = 4.5;
            double lamda = 1;
            double width = 1;
            double height = 1;
            SaturationProcess process = new SaturationProcess(theta1, theta2, r, lamda, width, height, c);
            process.EvolveStraussProcess();

            // output Strauss process record
            // n(x) u(x)
            fileName = "nx ux.txt";
            fileStream = new FileStream(fileName, FileMode.Create, FileAccess.Write);
            streamWriter = new StreamWriter(fileStream);
            int length = process.NumberPointsList.Count;
            for (int i = 0; i < length; i++)
            {
                streamWriter.WriteLine(String.Format("{0, -10:d}{1, -10:d}", (int)process.NumberPointsList[i], (int)process.NumberClosePointsList[i]));
            }
            streamWriter.Flush();
            streamWriter.Close();
            fileStream.Close();

            // birth accept rate
            fileName = "birth accept rate.txt";
            fileStream = new FileStream(fileName, FileMode.Create, FileAccess.Write);
            streamWriter = new StreamWriter(fileStream);
            length = process.BirthAcceptRateList.Count;
            for (int i = 0; i < length; i++)
            {
                streamWriter.WriteLine(String.Format("{0, -10:f}", (double)process.BirthAcceptRateList[i]));
            }
            streamWriter.Flush();
            streamWriter.Close();
            fileStream.Close();

            // death accept rate
            fileName = "death accept rate.txt";
            fileStream = new FileStream(fileName, FileMode.Create, FileAccess.Write);
            streamWriter = new StreamWriter(fileStream);
            length = process.DeathAcceptRateList.Count;
            for (int i = 0; i < length; i++)
            {
                streamWriter.WriteLine(String.Format("{0, -10:f}", (double)process.DeathAcceptRateList[i]));
            }
            streamWriter.Flush();
            streamWriter.Close();
            fileStream.Close();
        } // Example()
        #endregion
    }
}
