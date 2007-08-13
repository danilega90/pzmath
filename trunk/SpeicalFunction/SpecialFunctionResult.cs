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
    /* specfunc/gsl_sf_result.h
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

    public class SpecialFunctionResult
    {
        #region Fields
        private double _val;
        private double _err;
        #endregion

        #region Properties
        public double Val { get { return _val; } set { _val = value; } }
        public double Err { get { return _err; } set { _err = value; } }
        #endregion

        public void Set(double val, double err)
        {
            _val = val;
            _err = err;
        }
    }

    public class SpecialFunctionE10Result
    {
        #region Fields
        private double _val;
        private double _err;
        private int _e10;
        #endregion

        #region Properties
        public double Val { get { return _val; } set { _val = value; } }
        public double Err { get { return _err; } set { _err = value; } }
        public int E10 { get { return _e10; } set { _e10 = value; } }
        #endregion

        public void Set(double val, double err, int e10)
        {
            _val = val;
            _err = err;
            e10 = _e10;
        }
    }


}
