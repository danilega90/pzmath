// eee.Sheffield.PZ.Math
//
// Copyright ?Ping Zou, 2007
// sg71.cherub@gmail.com

using System;
using System.Collections.Generic;
using System.Text;

namespace eee.Sheffield.PZ.Math
{
    #region GSL comments
    /* specfunc/elementary.c
     * 
     * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
     * 
     * This program is free software; you can redistribute it and/or modify
     * it under the terms of the GNU General Public License as published by
     * the Free Software Foundation; either version 2 of the License, or (at
     * your option) any later version.
     * 
     * This program is distributed in the hope that it will be useful, but
     * WITHOUT ANY WARRANTY; without even the implied warranty of
     * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
     * General Public License for more details.
     * 
     * You should have received a copy of the GNU General Public License
     * along with this program; if not, write to the Free Software
     * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
     */

    /* Author:  G. Jungman */
    #endregion
    public class Elementary
    {
        #region static methods
        public static int MultiplyE(double x, double y, ref SpecialFunctionResult result)
        {
            double ax = System.Math.Abs(x);
            double ay = System.Math.Abs(y);

            if (x == 0.0 || y == 0.0)
            {
                /* It is necessary to eliminate this immediately.
                 */
                result.Val = 0.0;
                result.Err = 0.0;
                return PZMath_errno.PZMath_SUCCESS;
            }
            else if ((ax <= 1.0 && ay >= 1.0) || (ay <= 1.0 && ax >= 1.0))
            {
                /* Straddling 1.0 is always safe.
                 */
                result.Val = x * y;
                result.Err = 2.0 * PZMath_machine.PZMath_DBL_EPSILON * System.Math.Abs(result.Val);
                return PZMath_errno.PZMath_SUCCESS;
            }
            else
            {
                double f = 1.0 - 2.0 * PZMath_machine.PZMath_DBL_EPSILON;
                double min = System.Math.Min(System.Math.Abs(x), System.Math.Abs(y));
                double max = System.Math.Max(System.Math.Abs(x), System.Math.Abs(y));
                if (max < 0.9 * PZMath_machine.PZMath_SQRT_DBL_MAX || min < (f * PZMath_machine.PZMath_DBL_MAX) / max)
                {
                    result.Val = x * y;
                    result.Err = 2.0 * PZMath_machine.PZMath_DBL_EPSILON * System.Math.Abs(result.Val);
                    PZMath_errno.CheckUnderFlow(result);
                    return PZMath_errno.PZMath_SUCCESS;
                }
                else
                {
                    PZMath_errno.OverFlowError(ref result);
                }

            }
            return PZMath_errno.PZMath_SUCCESS;
        } // MultiplyE()


        public static int MultiplyErrE(double x, double dx, double y, double dy, ref SpecialFunctionResult result)
        {
            int status = MultiplyE(x, y, ref result);
            result.Err += System.Math.Abs(dx * y) + System.Math.Abs(dy * x);
            return status;
        } // MultiplyErrE()


        #endregion
    }
}
