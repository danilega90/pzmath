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
    /* specfunc/System.Math.Exp.c
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

    public class Exp
    {
        #region static methods

        /* Evaluate the continued fraction for exprel.
 * [Abramowitz+Stegun, 4.2.41]
 */
        public static int ExprelNCF(int N, double x, ref SpecialFunctionResult result)
        {
            double RECUR_BIG = PZMath_machine.PZMath_SQRT_DBL_MAX;
            int maxiter = 5000;
            int n = 1;
            double Anm2 = 1.0;
            double Bnm2 = 0.0;
            double Anm1 = 0.0;
            double Bnm1 = 1.0;
            double a1 = 1.0;
            double b1 = 1.0;
            double a2 = -x;
            double b2 = N + 1;
            double an, bn;

            double fn;

            double An = b1 * Anm1 + a1 * Anm2;   /* A1 */
            double Bn = b1 * Bnm1 + a1 * Bnm2;   /* B1 */

            /* One explicit step, before we get to the main pattern. */
            n++;
            Anm2 = Anm1;
            Bnm2 = Bnm1;
            Anm1 = An;
            Bnm1 = Bn;
            An = b2 * Anm1 + a2 * Anm2;   /* A2 */
            Bn = b2 * Bnm1 + a2 * Bnm2;   /* B2 */

            fn = An / Bn;

            while (n < maxiter)
            {
                double old_fn;
                double del;
                n++;
                Anm2 = Anm1;
                Bnm2 = Bnm1;
                Anm1 = An;
                Bnm1 = Bn;
                an = (PZMath.IsOdd(n) ? ((n - 1) / 2) * x : -(N + (n / 2) - 1) * x);
                bn = N + n - 1;
                An = bn * Anm1 + an * Anm2;
                Bn = bn * Bnm1 + an * Bnm2;

                if (System.Math.Abs(An) > RECUR_BIG || System.Math.Abs(Bn) > RECUR_BIG)
                {
                    An /= RECUR_BIG;
                    Bn /= RECUR_BIG;
                    Anm1 /= RECUR_BIG;
                    Bnm1 /= RECUR_BIG;
                    Anm2 /= RECUR_BIG;
                    Bnm2 /= RECUR_BIG;
                }

                old_fn = fn;
                fn = An / Bn;
                del = old_fn / fn;

                if (System.Math.Abs(del - 1.0) < 2.0 * PZMath_machine.PZMath_DBL_EPSILON)
                    break;
            }

            result.Val = fn;
            result.Err = 2.0 * (n + 1.0) * PZMath_machine.PZMath_DBL_EPSILON * System.Math.Abs(fn);

            if (n == maxiter)
                PZMath_errno.ERROR("error", PZMath_errno.PZMath_EMAXITER);

            return PZMath_errno.PZMath_SUCCESS;
        }


        /*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

        public static int ExpE(double x, ref SpecialFunctionResult result)
        {
            if (x > PZMath_machine.PZMath_LOG_DBL_MAX)
            {
                PZMath_errno.OverFlowError(ref result);
            }
            else if (x < PZMath_machine.PZMath_LOG_DBL_MIN)
            {
                PZMath_errno.UnderFlowError(ref result);
            }
            else
            {
                result.Val = System.Math.Exp(x);
                result.Err = 2.0 * PZMath_machine.PZMath_DBL_EPSILON * System.Math.Abs(result.Val);
                return PZMath_errno.PZMath_SUCCESS;
            }
            return PZMath_errno.PZMath_SUCCESS;
        } // ExpE()

        public static int ExpE10E(double x, ref SpecialFunctionE10Result result)
        {
            if (x > int.MaxValue - 1)
            {
                PZMath_errno.OverFlowErrorE10(ref result);
            }
            else if (x < int.MinValue + 1)
            {
                PZMath_errno.UnderFlowErrorE10(ref result);
            }
            else
            {
                int N = (int)System.Math.Floor(x / PZMath.M_LN10);
                result.Val = System.Math.Exp(x - N * PZMath.M_LN10);
                result.Err = 2.0 * (System.Math.Abs(x) + 1.0) * PZMath_machine.PZMath_DBL_EPSILON * System.Math.Abs(result.Val);
                result.E10 = N;
                return PZMath_errno.PZMath_SUCCESS;
            }
            return PZMath_errno.PZMath_SUCCESS;
        } // ExpE10E()


        public static int ExpMultE(double x, double y, ref SpecialFunctionResult result)
        {
            double ay = System.Math.Abs(y);

            if (y == 0.0)
            {
                result.Val = 0.0;
                result.Err = 0.0;
                return PZMath_errno.PZMath_SUCCESS;
            }
            else if ((x < 0.5 * PZMath_machine.PZMath_LOG_DBL_MAX && x > 0.5 * PZMath_machine.PZMath_LOG_DBL_MIN)
                    && (ay < 0.8 * PZMath_machine.PZMath_SQRT_DBL_MAX && ay > 1.2 * PZMath_machine.PZMath_SQRT_DBL_MIN)
              )
            {
                double ex = System.Math.Exp(x);
                result.Val = y * ex;
                result.Err = (2.0 + System.Math.Abs(x)) * PZMath_machine.PZMath_DBL_EPSILON * System.Math.Abs(result.Val);
                return PZMath_errno.PZMath_SUCCESS;
            }
            else
            {
                double ly = System.Math.Log(ay);
                double lnr = x + ly;

                if (lnr > PZMath_machine.PZMath_LOG_DBL_MAX - 0.01)
                {
                    PZMath_errno.OverFlowError(ref result);
                }
                else if (lnr < PZMath_machine.PZMath_LOG_DBL_MIN + 0.01)
                {
                    PZMath_errno.UnderFlowError(ref result);
                }
                else
                {
                    double sy = System.Math.Sign(y);
                    double M = System.Math.Floor(x);
                    double N = System.Math.Floor(ly);
                    double a = x - M;
                    double b = ly - N;
                    double berr = 2.0 * PZMath_machine.PZMath_DBL_EPSILON * (System.Math.Abs(ly) + System.Math.Abs(N));
                    result.Val = sy * System.Math.Exp(M + N) * System.Math.Exp(a + b);
                    result.Err = berr * System.Math.Abs(result.Val);
                    result.Err += 2.0 * PZMath_machine.PZMath_DBL_EPSILON * (M + N + 1.0) * System.Math.Abs(result.Val);
                    return PZMath_errno.PZMath_SUCCESS;
                }
            }
            return PZMath_errno.PZMath_SUCCESS;
        } // ExpMultE()


        public static int ExpMultE10E(double x, double y, ref SpecialFunctionE10Result result)
        {
            double ay = System.Math.Abs(y);

            if (y == 0.0)
            {
                result.Val = 0.0;
                result.Err = 0.0;
                result.E10 = 0;
                return PZMath_errno.PZMath_SUCCESS;
            }
            else if ((x < 0.5 * PZMath_machine.PZMath_LOG_DBL_MAX && x > 0.5 * PZMath_machine.PZMath_LOG_DBL_MIN)
                    && (ay < 0.8 * PZMath_machine.PZMath_SQRT_DBL_MAX && ay > 1.2 * PZMath_machine.PZMath_SQRT_DBL_MIN)
              )
            {
                double ex = System.Math.Exp(x);
                result.Val = y * ex;
                result.Err = (2.0 + System.Math.Abs(x)) * PZMath_machine.PZMath_DBL_EPSILON * System.Math.Abs(result.Val);
                result.E10 = 0;
                return PZMath_errno.PZMath_SUCCESS;
            }
            else
            {
                double ly = System.Math.Log(ay);
                double l10_val = (x + ly) / PZMath.M_LN10;

                if (l10_val > int.MaxValue - 1)
                {
                    PZMath_errno.OverFlowErrorE10(ref result);
                }
                else if (l10_val < int.MinValue + 1)
                {
                    PZMath_errno.UnderFlowErrorE10(ref result);
                }
                else
                {
                    double sy = System.Math.Sign(y);
                    int N = (int)System.Math.Floor(l10_val);
                    double arg_val = (l10_val - N) * PZMath.M_LN10;
                    double arg_err = 2.0 * PZMath_machine.PZMath_DBL_EPSILON * System.Math.Abs(ly);

                    result.Val = sy * System.Math.Exp(arg_val);
                    result.Err = arg_err * System.Math.Abs(result.Val);
                    result.Err += 2.0 * PZMath_machine.PZMath_DBL_EPSILON * System.Math.Abs(result.Val);
                    result.E10 = N;

                    return PZMath_errno.PZMath_SUCCESS;
                }
            }
            return PZMath_errno.PZMath_SUCCESS;
        } // ExpMultE10E()


        public static int ExpMultErrE(double x, double dx, double y, double dy, ref SpecialFunctionResult result)
        {
            double ay = System.Math.Abs(y);

            if (y == 0.0)
            {
                result.Val = 0.0;
                result.Err = System.Math.Abs(dy * System.Math.Exp(x));
                return PZMath_errno.PZMath_SUCCESS;
            }
            else if ((x < 0.5 * PZMath_machine.PZMath_LOG_DBL_MAX && x > 0.5 * PZMath_machine.PZMath_LOG_DBL_MIN)
                    && (ay < 0.8 * PZMath_machine.PZMath_SQRT_DBL_MAX && ay > 1.2 * PZMath_machine.PZMath_SQRT_DBL_MIN)
              )
            {
                double ex = System.Math.Exp(x);
                result.Val = y * ex;
                result.Err = ex * (System.Math.Abs(dy) + System.Math.Abs(y * dx));
                result.Err += 2.0 * PZMath_machine.PZMath_DBL_EPSILON * System.Math.Abs(result.Val);
                return PZMath_errno.PZMath_SUCCESS;
            }
            else
            {
                double ly = System.Math.Log(ay);
                double lnr = x + ly;

                if (lnr > PZMath_machine.PZMath_LOG_DBL_MAX - 0.01)
                {
                    PZMath_errno.OverFlowError(ref result);
                }
                else if (lnr < PZMath_machine.PZMath_LOG_DBL_MIN + 0.01)
                {
                    PZMath_errno.UnderFlowError(ref result);
                }
                else
                {
                    double sy = System.Math.Sign(y);
                    double M = System.Math.Floor(x);
                    double N = System.Math.Floor(ly);
                    double a = x - M;
                    double b = ly - N;
                    double eMN = System.Math.Exp(M + N);
                    double eab = System.Math.Exp(a + b);
                    result.Val = sy * eMN * eab;
                    result.Err = eMN * eab * 2.0 * PZMath_machine.PZMath_DBL_EPSILON;
                    result.Err += eMN * eab * System.Math.Abs(dy / y);
                    result.Err += eMN * eab * System.Math.Abs(dx);
                    return PZMath_errno.PZMath_SUCCESS;
                }
            }
            return PZMath_errno.PZMath_SUCCESS;
        } // ExpMultErrE()

        public static int ExpMultErrE10E(double x, double dx, double y, double dy, ref SpecialFunctionE10Result result)
        {
            double ay = System.Math.Abs(y);

            if (y == 0.0)
            {
                result.Val = 0.0;
                result.Err = System.Math.Abs(dy * System.Math.Exp(x));
                result.E10 = 0;
                return PZMath_errno.PZMath_SUCCESS;
            }
            else if ((x < 0.5 * PZMath_machine.PZMath_LOG_DBL_MAX && x > 0.5 * PZMath_machine.PZMath_LOG_DBL_MIN)
                    && (ay < 0.8 * PZMath_machine.PZMath_SQRT_DBL_MAX && ay > 1.2 * PZMath_machine.PZMath_SQRT_DBL_MIN)
              )
            {
                double ex = System.Math.Exp(x);
                result.Val = y * ex;
                result.Err = ex * (System.Math.Abs(dy) + System.Math.Abs(y * dx));
                result.Err += 2.0 * PZMath_machine.PZMath_DBL_EPSILON * System.Math.Abs(result.Val);
                result.E10 = 0;
                return PZMath_errno.PZMath_SUCCESS;
            }
            else
            {
                double ly = System.Math.Log(ay);
                double l10_val = (x + ly) / PZMath.M_LN10;

                if (l10_val > int.MaxValue - 1)
                {
                    PZMath_errno.OverFlowErrorE10(ref result);
                }
                else if (l10_val < int.MinValue + 1)
                {
                    PZMath_errno.UnderFlowErrorE10(ref result);
                }
                else
                {
                    double sy = System.Math.Sign(y);
                    int N = (int)System.Math.Floor(l10_val);
                    double arg_val = (l10_val - N) * PZMath.M_LN10;
                    double arg_err = dy / System.Math.Abs(y) + dx + 2.0 * PZMath_machine.PZMath_DBL_EPSILON * System.Math.Abs(arg_val);

                    result.Val = sy * System.Math.Exp(arg_val);
                    result.Err = arg_err * System.Math.Abs(result.Val);
                    result.Err += 2.0 * PZMath_machine.PZMath_DBL_EPSILON * System.Math.Abs(result.Val);
                    result.E10 = N;

                    return PZMath_errno.PZMath_SUCCESS;
                }
            }
            return PZMath_errno.PZMath_SUCCESS;
        } // ExpMultErrE10E()


        public static int Expm1E(double x, ref SpecialFunctionResult result)
        {
            double cut = 0.002;

            if (x < PZMath_machine.PZMath_LOG_DBL_MIN)
            {
                result.Val = -1.0;
                result.Err = PZMath_machine.PZMath_DBL_EPSILON;
                return PZMath_errno.PZMath_SUCCESS;
            }
            else if (x < -cut)
            {
                result.Val = System.Math.Exp(x) - 1.0;
                result.Err = 2.0 * PZMath_machine.PZMath_DBL_EPSILON * System.Math.Abs(result.Val);
                return PZMath_errno.PZMath_SUCCESS;
            }
            else if (x < cut)
            {
                result.Val = x * (1.0 + 0.5 * x * (1.0 + x / 3.0 * (1.0 + 0.25 * x * (1.0 + 0.2 * x))));
                result.Err = 2.0 * PZMath_machine.PZMath_DBL_EPSILON * System.Math.Abs(result.Val);
                return PZMath_errno.PZMath_SUCCESS;
            }
            else if (x < PZMath_machine.PZMath_LOG_DBL_MAX)
            {
                result.Val = System.Math.Exp(x) - 1.0;
                result.Err = 2.0 * PZMath_machine.PZMath_DBL_EPSILON * System.Math.Abs(result.Val);
                return PZMath_errno.PZMath_SUCCESS;
            }
            else
            {
                PZMath_errno.OverFlowError(ref result);
            }
            return PZMath_errno.PZMath_SUCCESS;
        } // Expm1E()


        public static int ExprelE(double x, ref SpecialFunctionResult result)
        {
            double cut = 0.002;

            if (x < PZMath_machine.PZMath_LOG_DBL_MIN)
            {
                result.Val = -1.0 / x;
                result.Err = PZMath_machine.PZMath_DBL_EPSILON * System.Math.Abs(result.Val);
                return PZMath_errno.PZMath_SUCCESS;
            }
            else if (x < -cut)
            {
                result.Val = (System.Math.Exp(x) - 1.0) / x;
                result.Err = 2.0 * PZMath_machine.PZMath_DBL_EPSILON * System.Math.Abs(result.Val);
                return PZMath_errno.PZMath_SUCCESS;
            }
            else if (x < cut)
            {
                result.Val = (1.0 + 0.5 * x * (1.0 + x / 3.0 * (1.0 + 0.25 * x * (1.0 + 0.2 * x))));
                result.Err = 2.0 * PZMath_machine.PZMath_DBL_EPSILON * System.Math.Abs(result.Val);
                return PZMath_errno.PZMath_SUCCESS;
            }
            else if (x < PZMath_machine.PZMath_LOG_DBL_MAX)
            {
                result.Val = (System.Math.Exp(x) - 1.0) / x;
                result.Err = 2.0 * PZMath_machine.PZMath_DBL_EPSILON * System.Math.Abs(result.Val);
                return PZMath_errno.PZMath_SUCCESS;
            }
            else
            {
                PZMath_errno.OverFlowError(ref result);
            }
            return PZMath_errno.PZMath_SUCCESS;
        } // ExprelE()


        public static int Exprel2E(double x, ref SpecialFunctionResult result)
        {
            double cut = 0.002;

            if (x < PZMath_machine.PZMath_LOG_DBL_MIN)
            {
                result.Val = -2.0 / x * (1.0 + 1.0 / x);
                result.Err = 2.0 * PZMath_machine.PZMath_DBL_EPSILON * System.Math.Abs(result.Val);
                return PZMath_errno.PZMath_SUCCESS;
            }
            else if (x < -cut)
            {
                result.Val = 2.0 * (System.Math.Exp(x) - 1.0 - x) / (x * x);
                result.Err = 2.0 * PZMath_machine.PZMath_DBL_EPSILON * System.Math.Abs(result.Val);
                return PZMath_errno.PZMath_SUCCESS;
            }
            else if (x < cut)
            {
                result.Val = (1.0 + 1.0 / 3.0 * x * (1.0 + 0.25 * x * (1.0 + 0.2 * x * (1.0 + 1.0 / 6.0 * x))));
                result.Err = 2.0 * PZMath_machine.PZMath_DBL_EPSILON * System.Math.Abs(result.Val);
                return PZMath_errno.PZMath_SUCCESS;
            }
            else if (x < PZMath_machine.PZMath_LOG_DBL_MAX)
            {
                result.Val = 2.0 * (System.Math.Exp(x) - 1.0 - x) / (x * x);
                result.Err = 2.0 * PZMath_machine.PZMath_DBL_EPSILON * System.Math.Abs(result.Val);
                return PZMath_errno.PZMath_SUCCESS;
            }
            else
            {
                PZMath_errno.OverFlowError(ref result);
            }
            return PZMath_errno.PZMath_SUCCESS;

        } // Exprel2E()


        public static int ExprelNE(int N, double x, ref SpecialFunctionResult result)
        {
            if (N < 0)
            {
                PZMath_errno.DomainError(ref result);
            }
            else if (x == 0.0)
            {
                result.Val = 1.0;
                result.Err = 0.0;
                return PZMath_errno.PZMath_SUCCESS;
            }
            else if (System.Math.Abs(x) < PZMath_machine.PZMath_ROOT3_DBL_EPSILON * N)
            {
                result.Val = 1.0 + x / (N + 1) * (1.0 + x / (N + 2));
                result.Err = 2.0 * PZMath_machine.PZMath_DBL_EPSILON;
                return PZMath_errno.PZMath_SUCCESS;
            }
            else if (N == 0)
            {
                return ExpE(x, ref result);
            }
            else if (N == 1)
            {
                return ExprelE(x, ref result);
            }
            else if (N == 2)
            {
                return Exprel2E(x, ref result);
            }
            else
            {
                if (x > N && (-x + N * (1.0 + System.Math.Log(x / N)) < PZMath_machine.PZMath_LOG_DBL_EPSILON))
                {
                    /* x is much larger than n.
                     * Ignore polynomial part, so
                     * exprel_N(x) ~= e^x N!/x^N
                     */
                    SpecialFunctionResult lnf_N = new SpecialFunctionResult();
                    double lnr_val;
                    double lnr_err;
                    double lnterm;
                    Gamma.LnFactE(N, ref lnf_N);
                    lnterm = N * System.Math.Log(x);
                    lnr_val = x + lnf_N.Val - lnterm;
                    lnr_err = PZMath_machine.PZMath_DBL_EPSILON * (System.Math.Abs(x) + System.Math.Abs(lnf_N.Val) + System.Math.Abs(lnterm));
                    lnr_err += lnf_N.Err;
                    return ExpErrE(lnr_val, lnr_err, ref result);
                }
                else if (x > N)
                {
                    /* Write the identity
                     *   exprel_n(x) = e^x n! / x^n (1 - Gamma[n,x]/Gamma[n])
                     * then use the asymptotic expansion
                     * Gamma[n,x] ~ x^(n-1) e^(-x) (1 + (n-1)/x + (n-1)(n-2)/x^2 + ...)
                     */
                    double ln_x = System.Math.Log(x);

                    SpecialFunctionResult lnf_N = new SpecialFunctionResult();
                    double lg_N;
                    double lnpre_val;
                    double lnpre_err;
                    Gamma.LnFactE(N, ref lnf_N);    /* System.Math.Log(N!)       */
                    lg_N = lnf_N.Val - System.Math.Log(N);       /* System.Math.Log(Gamma(N)) */
                    lnpre_val = x + lnf_N.Val - N * ln_x;
                    lnpre_err = PZMath_machine.PZMath_DBL_EPSILON * (System.Math.Abs(x) + System.Math.Abs(lnf_N.Val) + System.Math.Abs(N * ln_x));
                    lnpre_err += lnf_N.Err;
                    if (lnpre_val < PZMath_machine.PZMath_LOG_DBL_MAX - 5.0)
                    {
                        int stat_eG;
                        SpecialFunctionResult bigG_ratio = new SpecialFunctionResult();
                        SpecialFunctionResult pre = new SpecialFunctionResult();
                        int stat_ex = ExpErrE(lnpre_val, lnpre_err, ref pre);
                        double ln_bigG_ratio_pre = -x + (N - 1) * ln_x - lg_N;
                        double bigGsum = 1.0;
                        double term = 1.0;
                        int k;
                        for (k = 1; k < N; k++)
                        {
                            term *= (N - k) / x;
                            bigGsum += term;
                        }
                        stat_eG = ExpMultE(ln_bigG_ratio_pre, bigGsum, ref bigG_ratio);
                        if (stat_eG == PZMath_errno.PZMath_SUCCESS)
                        {
                            result.Val = pre.Val * (1.0 - bigG_ratio.Val);
                            result.Err = pre.Val * (2.0 * PZMath_machine.PZMath_DBL_EPSILON + bigG_ratio.Err);
                            result.Err += pre.Err * System.Math.Abs(1.0 - bigG_ratio.Val);
                            result.Err += 2.0 * PZMath_machine.PZMath_DBL_EPSILON * System.Math.Abs(result.Val);
                            return stat_ex;
                        }
                        else
                        {
                            result.Val = 0.0;
                            result.Err = 0.0;
                            return stat_eG;
                        }
                    }
                    else
                    {
                        PZMath_errno.OverFlowError(ref result);
                    }
                }
                else if (x > -10.0 * N)
                {
                    return ExprelNCF(N, x, ref result);
                }
                else
                {
                    /* x -> -Inf asymptotic:
                     * exprel_n(x) ~ e^x n!/x^n - n/x (1 + (n-1)/x + (n-1)(n-2)/x + ...)
                     *             ~ - n/x (1 + (n-1)/x + (n-1)(n-2)/x + ...)
                     */
                    double sum = 1.0;
                    double term = 1.0;
                    int k;
                    for (k = 1; k < N; k++)
                    {
                        term *= (N - k) / x;
                        sum += term;
                    }
                    result.Val = -N / x * sum;
                    result.Err = 2.0 * PZMath_machine.PZMath_DBL_EPSILON * System.Math.Abs(result.Val);
                    return PZMath_errno.PZMath_SUCCESS;
                }
            }
            return PZMath_errno.PZMath_SUCCESS;
        } // ExprelNE()


        public static int ExpErrE(double x, double dx, ref SpecialFunctionResult result)
        {
            double adx = System.Math.Abs(dx);

            if (x + adx > PZMath_machine.PZMath_LOG_DBL_MAX)
            {
                PZMath_errno.OverFlowError(ref result);
            }
            else if (x - adx < PZMath_machine.PZMath_LOG_DBL_MIN)
            {
                PZMath_errno.UnderFlowError(ref result);
            }
            else
            {
                double ex = System.Math.Exp(x);
                double edx = System.Math.Exp(adx);
                result.Val = ex;
                result.Err = ex * System.Math.Max(PZMath_machine.PZMath_DBL_EPSILON, edx - 1.0 / edx);
                result.Err += 2.0 * PZMath_machine.PZMath_DBL_EPSILON * System.Math.Abs(result.Val);
                return PZMath_errno.PZMath_SUCCESS;
            }
            return PZMath_errno.PZMath_SUCCESS;
        } // ExpErrE()


        public static int ExpErrE10E(double x, double dx, ref SpecialFunctionE10Result result)
        {
            double adx = System.Math.Abs(dx);

            if (x + adx > int.MaxValue - 1)
            {
                PZMath_errno.OverFlowErrorE10(ref result);
            }
            else if (x - adx < int.MinValue + 1)
            {
                PZMath_errno.UnderFlowErrorE10(ref result);
            }
            else
            {
                int N = (int)System.Math.Floor(x / PZMath.M_LN10);
                double ex = System.Math.Exp(x - N * PZMath.M_LN10);
                result.Val = ex;
                result.Err = ex * (2.0 * PZMath_machine.PZMath_DBL_EPSILON * (System.Math.Abs(x) + 1.0) + adx);
                result.E10 = N;
                return PZMath_errno.PZMath_SUCCESS;
            }
            return PZMath_errno.PZMath_SUCCESS;
        } // ExpErrE10E()

        #endregion

    }
}
