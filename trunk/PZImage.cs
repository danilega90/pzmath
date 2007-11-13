// PZ image 
// 13.12.2005
using System;
using System.Runtime.Serialization.Formatters.Binary;
using System.Runtime.Serialization;
using System.IO;


namespace eee.Sheffield.PZ.Imaging
{
	#region PZPoint class
	/// <summary>
	/// point structure and relative operation
	/// </summary>
    [Serializable]
	public class PZPoint : ICloneable
	{
		public double x;
		public double y;
		public double z;

		#region constructors

		/// <summary>
		/// empty constructor
		/// </summary>
		public PZPoint()
		{
			x = 0.0;
			y = 0.0;
			z = 0.0;
		} // PZPoint()

		public PZPoint(PZPoint p)
		{
			x = p.x;
			y = p.y;
			z = p.z;
		} // PZPoint(PZPoint)


		/// <summary>
		/// 2D constructor
		/// </summary>
		/// <param name="inputx"></param>
		/// <param name="inputy"></param>
		public PZPoint(double inputx, double inputy) 
		{
			x = inputx;
			y = inputy;
			z = 0.0;
		} // PZPoint(double, double)

		/// <summary>
		/// 3D constructor
		/// </summary>
		/// <param name="inputx"></param>
		/// <param name="inputy"></param>
		/// <param name="inputz"></param>
		public PZPoint(double inputx, double inputy, double inputz)
		{
			x = inputx;
			y = inputy;
			z = inputz;
		} // PZPoint(double, double, double)
		#endregion

		public double Distance(PZPoint a)
			// distance between a and current point
		{
			double d = (x - a.x) * (x - a.x);
			d += (y - a.y) * (y - a.y);
			d += (z - a.z) * (z - a.z);
            return System.Math.Sqrt(d);
		} // Distance()

		public void MemCopyFrom(PZPoint p)
		{
			x = p.x;
			y = p.y;
			z = p.z;
		} // MemCopyFrom()

		/// <summary>
		/// extend from a given point, direction and distance
		/// </summary>
		/// <param name="p"></param>
		/// <param name="d"></param>
		/// <param name="dist"></param>
		/// <returns></returns>
		public void ExtendFrom(PZPoint p, PZDirection d, double dist)
		{
			x = p.x + dist * d.x;
			y = p.y + dist * d.y;
			z = p.z + dist * d.z;
		} // ExtendFrom()

        public bool Is8NeighbourOf(PZPoint p)
        {
            bool is8Neighbour = false;
            if ((x == p.x - 1 && y == p.y - 1)
                || (x == p.x - 1 && y == p.y)
                || (x == p.x - 1 && y == p.y + 1)
                || (x == p.x && y == p.y - 1)
                || (x == p.x && y == p.y + 1)
                || (x == p.x + 1 && y == p.y - 1)
                || (x == p.x + 1 && y == p.y)
                || (x == p.x + 1 && y == p.y + 1))
                is8Neighbour = true;
            return is8Neighbour;
        } // Is8NeighbourOf()

        public bool Is4NeighbourOf(PZPoint p)
        {
            bool is4Neighbour = false;
            if ((x == p.x - 1 && y == p.y)
                || (x == p.x + 1 && y == p.y)
                || (x == p.x && y == p.y - 1)
                || (x == p.x && y == p.y + 1))
                is4Neighbour = true;

            return is4Neighbour;
        } // Is4NeighbourOf()

		public override string ToString()
		{
			return "(" + String.Format("{0:0.00}",x) + "," + String.Format("{0:0.00}",y) + ")";
		}

        public bool EqualTo(PZPoint p)
        {
            if (x == p.x && y == p.y && z == p.z)
                return true;
            else
                return false;
        }

        #region Clone method
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
        #endregion
    } // PZPoint
	#endregion

	#region PZDirection class
	/// <summary>
	/// direction structure, with unit length, i.e.
	/// x ^ 2 + y ^ 2 + z ^ 2 = 1
	/// </summary>
    [Serializable]
	public class PZDirection : ICloneable
	{
		public double x;
		public double y;
		public double z;

		#region constructors
		/// <summary>
		/// empty constructor
		/// </summary>
		public PZDirection()
		{
			x = 1.0;
			y = 0.0;
			z = 0.0;
		} // PZDirection()

		public PZDirection(PZDirection d)
		{
			x = d.x;
			y = d.y;
			z = d.z;
		} // PZDirection(PZDirection);

		/// <summary>
		/// 2D constructor
		/// </summary>
		/// <param name="inputx"></param>
		/// <param name="inputy"></param>
		public PZDirection(double inputx, double inputy) 
		{
            double length = System.Math.Sqrt(inputx * inputx + inputy * inputy);
			x = inputx / length;
			y = inputy / length;
			z = 0.0;
        } // PZDirection(double, double)

		/// <summary>
		/// 3D constructor
		/// </summary>
		/// <param name="inputx"></param>
		/// <param name="inputy"></param>
		/// <param name="inputz"></param>
		public PZDirection(double inputx, double inputy, double inputz)
		{
            double length = System.Math.Sqrt(inputx * inputx + inputy * inputy + inputz * inputz);
			x = inputx / length;
			y = inputy / length;
			z = inputz / length;
        } // PZDirection(double, double, double)

        /// <summary>
        /// direction of start poing to end point
        /// </summary>
        /// <param name="start"></param>
        /// <param name="end"></param>
        public PZDirection(PZPoint start, PZPoint end)
        {
            double dx = end.x - start.x;
            double dy = end.y - start.y;
            double length = System.Math.Sqrt(dx * dx + dy * dy);
            x = dx / length;
            y = dy / length;
            z = 0.0;
        } // PZDirection(PZPoint, PZPoint)
		#endregion
		
		/// <summary>
		/// get direction from 2 points
		/// </summary>
		/// <param name="p1">start point</param>
		/// <param name="p2">end point</param>
		/// <returns></returns>
		public void GetDirectionFrom2Points(PZPoint p1, PZPoint p2)
		{
			double dx = p2.x - p1.x;
			double dy = p2.y - p1.y;
            double length = System.Math.Sqrt(dx * dx + dy * dy);
			x = dx / length;
			y = dy / length;
		} // GetDirectionFrom2Points

		/// <summary>
		/// get the normal direction
		/// </summary>
		/// <returns></returns>
		public void GetNormalDirection()
		{
			double tempx = y;
			double tempy = -1.0 * x;
			x = tempx;
			y = tempy;
		} // GetNormalDirection()

		/// <summary>
		/// get angle to direction d in radius
		/// </summary>
		/// <param name="d"></param>
		/// <returns></returns>
		public double GetAngle(PZDirection d)
		{
			double cos = (x * d.x + y * d.y) 
				/ System.Math.Sqrt(x * x + y * y) 
				/ System.Math.Sqrt(d.x * d.x + d.y * d.y);
            
            if (cos > 1.0)
                cos = 1.0;
            if (cos < -1.0)
                cos = 1.0;            

            double angle = System.Math.Acos(cos);
			return angle;
		} // GetAngle()

        public double GetAngle()
        {
            PZDirection h = new PZDirection(1, 0);
            double a = this.GetAngle(h);
            return a;
        } // GetAngle()

        /// <summary>
        /// countclockwise for positive theta in degree
        /// </summary>
        /// <param name="theta"></param>
		public void Rotate(double theta)
			// countclockwise 2D
		{
            double theta_rad = theta / 180.0 * System.Math.PI;
            double temp_x = x * System.Math.Cos(theta_rad) - y * System.Math.Sin(theta_rad);
            double temp_y = x * System.Math.Sin(theta_rad) + y * System.Math.Cos(theta_rad);
			x = temp_x;
			y = temp_y;
		} // Rotate()

		public void MemCopyFrom(PZDirection d)
		{
			x = d.x;
			y = d.y;
			z = d.z;
		} // MemCopyFrom()

		public override string ToString()
		{
			return "(" + String.Format("{0:0.00}",x) + "," + String.Format("{0:0.00}",y) + ")";
		} // ToString()


        #region Clone method
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
        #endregion

        #region Example
        public static void Example()
        {
            PZDirection d1 = new PZDirection(1, 1);     // 45d
            PZDirection d2 = new PZDirection(-1, 1);    // 135d
            PZDirection d3 = new PZDirection(-1, -1);   // 225d
            PZDirection d4 = new PZDirection(1, -1);    // -45d or 315d
            System.Console.WriteLine("d1: " + d1.GetAngle());
            System.Console.WriteLine("d2: " + d2.GetAngle());
            System.Console.WriteLine("d3: " + d3.GetAngle());
            System.Console.WriteLine("d4: " + d4.GetAngle());
        } // Example()
        #endregion
    } // PZDirection
	#endregion

    public class PZPointPair
    {
        public PZPoint P1;
        public PZPoint P2;

        public PZPointPair(PZPoint p1, PZPoint p2)
        {
            P1 = new PZPoint(p1);
            P2 = new PZPoint(p2);
        }
    }
}
