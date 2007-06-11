// eee.Sheffield.PZ.Math
//
// Copyright ? Ping Zou, 2007
// sg71.cherub@gmail.com
using System;
using System.Collections.Generic;
using System.Text;

namespace eee.Sheffield.PZ.Math
{
    /// <summary>
    /// generate random index (unrepeated) in range [min, max]
    /// </summary>
    public class RandomIndex
    {
        #region static methods
        /// <summary>
        /// static method
        /// </summary>
        /// <param name="min"></param>
        /// <param name="max"></param>
        /// <returns></returns>
        public static int[] GenerateRandomIndex(int min, int max)
        {
            if (min >= max)
                throw new ApplicationException("RandomIndex::GenerateRandomIndex(), lower bound is greater than or equals to upper bound.");
            int length = max - min + 1;
            int[] randomIndex = new int[length];
            List<int> temp = new List<int>(length);
            for (int i = 0; i < length; i++)
                temp.Add(i);

            // uniform random number
            long seed = DateTime.Now.Second + DateTime.Now.Millisecond;
            UniformDistribution random = new UniformDistribution(seed);
            for (int i = 0; i < length; i++)
            {
                int index = Convert.ToInt32(random.Sample() * (temp.Count - 1)); 
                randomIndex[i] = min + temp[index];
                temp.Remove(temp[index]);
            }
            return randomIndex;
        } // GenerateRandomIndex()

        public static int[] GenerateRandomIndex(int length)
        {
            return RandomIndex.GenerateRandomIndex(0, length - 1);
        }
        #endregion

        #region example
        public static void Example()
        {
            int[] randomIndex = RandomIndex.GenerateRandomIndex(20);
            for (int i = 0; i < 20; i++)
                System.Console.WriteLine(randomIndex[i]);
        }
        #endregion
    }
}
