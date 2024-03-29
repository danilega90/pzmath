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
    /* specfunc/zeta.c
     * 
     * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004 Gerard Jungman
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

    public class Zeta
    {
        #region static parameters
        public static double LogTwoPi_ = 1.8378770664093454835606594728111235279723;

        /* chebyshev fit for (s(t)-1)Zeta[s(t)]
         * s(t)= (t+1)/2
         * -1 <= t <= 1
         */

        static double[] zeta_xlt1_data = {
          1.48018677156931561235192914649,
          0.25012062539889426471999938167,
          0.00991137502135360774243761467,
         -0.00012084759656676410329833091,
         -4.7585866367662556504652535281e-06,
          2.2229946694466391855561441361e-07,
         -2.2237496498030257121309056582e-09,
         -1.0173226513229028319420799028e-10,
          4.3756643450424558284466248449e-12,
         -6.2229632593100551465504090814e-14,
         -6.6116201003272207115277520305e-16,
          4.9477279533373912324518463830e-17,
         -1.0429819093456189719660003522e-18,
          6.9925216166580021051464412040e-21,
        };

        static ChebyshevSeries zeta_xlt1_cs = new ChebyshevSeries(zeta_xlt1_data, 13, -1, 1, 8);

        /* chebyshev fit for (s(t)-1)Zeta[s(t)]
         * s(t)= (19t+21)/2
         * -1 <= t <= 1
         */
        static double[] zeta_xgt1_data = {
          19.3918515726724119415911269006,
           9.1525329692510756181581271500,
           0.2427897658867379985365270155,
          -0.1339000688262027338316641329,
           0.0577827064065028595578410202,
          -0.0187625983754002298566409700,
           0.0039403014258320354840823803,
          -0.0000581508273158127963598882,
          -0.0003756148907214820704594549,
           0.0001892530548109214349092999,
          -0.0000549032199695513496115090,
           8.7086484008939038610413331863e-6,
           6.4609477924811889068410083425e-7,
          -9.6749773915059089205835337136e-7,
           3.6585400766767257736982342461e-7,
          -8.4592516427275164351876072573e-8,
           9.9956786144497936572288988883e-9,
           1.4260036420951118112457144842e-9,
          -1.1761968823382879195380320948e-9,
           3.7114575899785204664648987295e-10,
          -7.4756855194210961661210215325e-11,
           7.8536934209183700456512982968e-12,
           9.9827182259685539619810406271e-13,
          -7.5276687030192221587850302453e-13,
           2.1955026393964279988917878654e-13,
          -4.1934859852834647427576319246e-14,
           4.6341149635933550715779074274e-15,
           2.3742488509048340106830309402e-16,
          -2.7276516388124786119323824391e-16,
           7.8473570134636044722154797225e-17
        };

        static ChebyshevSeries zeta_xgt1_cs = new ChebyshevSeries(zeta_xgt1_data, 29, -1, 1, 17);


        /* chebyshev fit for Ln[Zeta[s(t)] - 1 - 2^(-s(t))]
         * s(t)= 10 + 5t
         * -1 <= t <= 1; 5 <= s <= 15
         */
        static double[] zetam1_inter_data = {
          -21.7509435653088483422022339374,
          -5.63036877698121782876372020472,
           0.0528041358684229425504861579635,
          -0.0156381809179670789342700883562,
           0.00408218474372355881195080781927,
          -0.0010264867349474874045036628282,
           0.000260469880409886900143834962387,
          -0.0000676175847209968878098566819447,
           0.0000179284472587833525426660171124,
          -4.83238651318556188834107605116e-6,
           1.31913788964999288471371329447e-6,
          -3.63760500656329972578222188542e-7,
           1.01146847513194744989748396574e-7,
          -2.83215225141806501619105289509e-8,
           7.97733710252021423361012829496e-9,
          -2.25850168553956886676250696891e-9,
           6.42269392950164306086395744145e-10,
          -1.83363861846127284505060843614e-10,
           5.25309763895283179960368072104e-11,
          -1.50958687042589821074710575446e-11,
           4.34997545516049244697776942981e-12,
          -1.25597782748190416118082322061e-12,
           3.61280740072222650030134104162e-13,
          -9.66437239205745207188920348801e-14
        };

        static ChebyshevSeries zetam1_inter_cs = new ChebyshevSeries(zetam1_inter_data, 22, -1, 1, 12);



        static int ZETA_POS_TABLE_NMAX = 100;

        static double[] zetam1_pos_int_table = {
         -1.5,                               /* zeta(0) */
          0.0,       /* FIXME: Infinity */   /* zeta(1) - 1 */
          0.644934066848226436472415166646,  /* zeta(2) - 1 */
          0.202056903159594285399738161511,
          0.082323233711138191516003696541,
          0.036927755143369926331365486457,
          0.017343061984449139714517929790,
          0.008349277381922826839797549849,
          0.004077356197944339378685238508,
          0.002008392826082214417852769232,
          0.000994575127818085337145958900,
          0.000494188604119464558702282526,
          0.000246086553308048298637998047,
          0.000122713347578489146751836526,
          0.000061248135058704829258545105,
          0.000030588236307020493551728510,
          0.000015282259408651871732571487,
          7.6371976378997622736002935630e-6,
          3.8172932649998398564616446219e-6,
          1.9082127165539389256569577951e-6,
          9.5396203387279611315203868344e-7,
          4.7693298678780646311671960437e-7,
          2.3845050272773299000364818675e-7,
          1.1921992596531107306778871888e-7,
          5.9608189051259479612440207935e-8,
          2.9803503514652280186063705069e-8,
          1.4901554828365041234658506630e-8,
          7.4507117898354294919810041706e-9,
          3.7253340247884570548192040184e-9,
          1.8626597235130490064039099454e-9,
          9.3132743241966818287176473502e-10,
          4.6566290650337840729892332512e-10,
          2.3283118336765054920014559759e-10,
          1.1641550172700519775929738354e-10,
          5.8207720879027008892436859891e-11,
          2.9103850444970996869294252278e-11,
          1.4551921891041984235929632245e-11,
          7.2759598350574810145208690123e-12,
          3.6379795473786511902372363558e-12,
          1.8189896503070659475848321007e-12,
          9.0949478402638892825331183869e-13,
          4.5474737830421540267991120294e-13,
          2.2737368458246525152268215779e-13,
          1.1368684076802278493491048380e-13,
          5.6843419876275856092771829675e-14,
          2.8421709768893018554550737049e-14,
          1.4210854828031606769834307141e-14,
          7.1054273952108527128773544799e-15,
          3.5527136913371136732984695340e-15,
          1.7763568435791203274733490144e-15,
          8.8817842109308159030960913863e-16,
          4.4408921031438133641977709402e-16,
          2.2204460507980419839993200942e-16,
          1.1102230251410661337205445699e-16,
          5.5511151248454812437237365905e-17,
          2.7755575621361241725816324538e-17,
          1.3877787809725232762839094906e-17,
          6.9388939045441536974460853262e-18,
          3.4694469521659226247442714961e-18,
          1.7347234760475765720489729699e-18,
          8.6736173801199337283420550673e-19,
          4.3368086900206504874970235659e-19,
          2.1684043449972197850139101683e-19,
          1.0842021724942414063012711165e-19,
          5.4210108624566454109187004043e-20,
          2.7105054312234688319546213119e-20,
          1.3552527156101164581485233996e-20,
          6.7762635780451890979952987415e-21,
          3.3881317890207968180857031004e-21,
          1.6940658945097991654064927471e-21,
          8.4703294725469983482469926091e-22,
          4.2351647362728333478622704833e-22,
          2.1175823681361947318442094398e-22,
          1.0587911840680233852265001539e-22,
          5.2939559203398703238139123029e-23,
          2.6469779601698529611341166842e-23,
          1.3234889800848990803094510250e-23,
          6.6174449004244040673552453323e-24,
          3.3087224502121715889469563843e-24,
          1.6543612251060756462299236771e-24,
          8.2718061255303444036711056167e-25,
          4.1359030627651609260093824555e-25,
          2.0679515313825767043959679193e-25,
          1.0339757656912870993284095591e-25,
          5.1698788284564313204101332166e-26,
          2.5849394142282142681277617708e-26,
          1.2924697071141066700381126118e-26,
          6.4623485355705318034380021611e-27,
          3.2311742677852653861348141180e-27,
          1.6155871338926325212060114057e-27,
          8.0779356694631620331587381863e-28,
          4.0389678347315808256222628129e-28,
          2.0194839173657903491587626465e-28,
          1.0097419586828951533619250700e-28,
          5.0487097934144756960847711725e-29,
          2.5243548967072378244674341938e-29,
          1.2621774483536189043753999660e-29,
          6.3108872417680944956826093943e-30,
          3.1554436208840472391098412184e-30,
          1.5777218104420236166444327830e-30,
          7.8886090522101180735205378276e-31
        };


        static int ZETA_NEG_TABLE_NMAX = 99;
        //static int ZETA_NEG_TABLE_SIZE = 50;
        static double[] zeta_neg_int_table = {
         -0.083333333333333333333333333333,     /* zeta(-1) */
          0.008333333333333333333333333333,     /* zeta(-3) */
         -0.003968253968253968253968253968,     /* ...      */
          0.004166666666666666666666666667,
         -0.007575757575757575757575757576,
          0.021092796092796092796092796093,
         -0.083333333333333333333333333333,
          0.44325980392156862745098039216,
         -3.05395433027011974380395433027,
          26.4562121212121212121212121212,
         -281.460144927536231884057971014,
          3607.5105463980463980463980464,
         -54827.583333333333333333333333,
          974936.82385057471264367816092,
         -2.0052695796688078946143462272e+07,
          4.7238486772162990196078431373e+08,
         -1.2635724795916666666666666667e+10,
          3.8087931125245368811553022079e+11,
         -1.2850850499305083333333333333e+13,
          4.8241448354850170371581670362e+14,
         -2.0040310656516252738108421663e+16,
          9.1677436031953307756992753623e+17,
         -4.5979888343656503490437943262e+19,
          2.5180471921451095697089023320e+21,
         -1.5001733492153928733711440151e+23,
          9.6899578874635940656497942895e+24,
         -6.7645882379292820990945242302e+26,
          5.0890659468662289689766332916e+28,
         -4.1147288792557978697665486068e+30,
          3.5666582095375556109684574609e+32,
         -3.3066089876577576725680214670e+34,
          3.2715634236478716264211227016e+36,
         -3.4473782558278053878256455080e+38,
          3.8614279832705258893092720200e+40,
         -4.5892974432454332168863989006e+42,
          5.7775386342770431824884825688e+44,
         -7.6919858759507135167410075972e+46,
          1.0813635449971654696354033351e+49,
         -1.6029364522008965406067102346e+51,
          2.5019479041560462843656661499e+53,
         -4.1067052335810212479752045004e+55,
          7.0798774408494580617452972433e+57,
         -1.2804546887939508790190849756e+60,
          2.4267340392333524078020892067e+62,
         -4.8143218874045769355129570066e+64,
          9.9875574175727530680652777408e+66,
         -2.1645634868435185631335136160e+69,
          4.8962327039620553206849224516e+71,    /* ...        */
         -1.1549023923963519663954271692e+74,    /* zeta(-97)  */
          2.8382249570693706959264156336e+76     /* zeta(-99)  */
        };


        /* coefficients for Maclaurin summation in hzeta()
         * B_{2j}/(2j)!
         */
        static double[] hzeta_c = {
          1.00000000000000000000000000000,
          0.083333333333333333333333333333,
         -0.00138888888888888888888888888889,
          0.000033068783068783068783068783069,
         -8.2671957671957671957671957672e-07,
          2.0876756987868098979210090321e-08,
         -5.2841901386874931848476822022e-10,
          1.3382536530684678832826980975e-11,
         -3.3896802963225828668301953912e-13,
          8.5860620562778445641359054504e-15,
         -2.1748686985580618730415164239e-16,
          5.5090028283602295152026526089e-18,
         -1.3954464685812523340707686264e-19,
          3.5347070396294674716932299778e-21,
         -8.9535174270375468504026113181e-23
        };

        static int ETA_POS_TABLE_NMAX = 100;
        static double[] eta_pos_int_table = {
        0.50000000000000000000000000000,  /* eta(0) */
        PZMath.M_LN2,                            /* eta(1) */
        0.82246703342411321823620758332,  /* ...    */
        0.90154267736969571404980362113,
        0.94703282949724591757650323447,
        0.97211977044690930593565514355,
        0.98555109129743510409843924448,
        0.99259381992283028267042571313,
        0.99623300185264789922728926008,
        0.99809429754160533076778303185,
        0.99903950759827156563922184570,
        0.99951714349806075414409417483,
        0.99975768514385819085317967871,
        0.99987854276326511549217499282,
        0.99993917034597971817095419226,
        0.99996955121309923808263293263,
        0.99998476421490610644168277496,
        0.99999237829204101197693787224,
        0.99999618786961011347968922641,
        0.99999809350817167510685649297,
        0.99999904661158152211505084256,
        0.99999952325821554281631666433,
        0.99999976161323082254789720494,
        0.99999988080131843950322382485,
        0.99999994039889239462836140314,
        0.99999997019885696283441513311,
        0.99999998509923199656878766181,
        0.99999999254955048496351585274,
        0.99999999627475340010872752767,
        0.99999999813736941811218674656,
        0.99999999906868228145397862728,
        0.99999999953434033145421751469,
        0.99999999976716989595149082282,
        0.99999999988358485804603047265,
        0.99999999994179239904531592388,
        0.99999999997089618952980952258,
        0.99999999998544809143388476396,
        0.99999999999272404460658475006,
        0.99999999999636202193316875550,
        0.99999999999818101084320873555,
        0.99999999999909050538047887809,
        0.99999999999954525267653087357,
        0.99999999999977262633369589773,
        0.99999999999988631316532476488,
        0.99999999999994315658215465336,
        0.99999999999997157829090808339,
        0.99999999999998578914539762720,
        0.99999999999999289457268000875,
        0.99999999999999644728633373609,
        0.99999999999999822364316477861,
        0.99999999999999911182158169283,
        0.99999999999999955591079061426,
        0.99999999999999977795539522974,
        0.99999999999999988897769758908,
        0.99999999999999994448884878594,
        0.99999999999999997224442439010,
        0.99999999999999998612221219410,
        0.99999999999999999306110609673,
        0.99999999999999999653055304826,
        0.99999999999999999826527652409,
        0.99999999999999999913263826204,
        0.99999999999999999956631913101,
        0.99999999999999999978315956551,
        0.99999999999999999989157978275,
        0.99999999999999999994578989138,
        0.99999999999999999997289494569,
        0.99999999999999999998644747284,
        0.99999999999999999999322373642,
        0.99999999999999999999661186821,
        0.99999999999999999999830593411,
        0.99999999999999999999915296705,
        0.99999999999999999999957648353,
        0.99999999999999999999978824176,
        0.99999999999999999999989412088,
        0.99999999999999999999994706044,
        0.99999999999999999999997353022,
        0.99999999999999999999998676511,
        0.99999999999999999999999338256,
        0.99999999999999999999999669128,
        0.99999999999999999999999834564,
        0.99999999999999999999999917282,
        0.99999999999999999999999958641,
        0.99999999999999999999999979320,
        0.99999999999999999999999989660,
        0.99999999999999999999999994830,
        0.99999999999999999999999997415,
        0.99999999999999999999999998708,
        0.99999999999999999999999999354,
        0.99999999999999999999999999677,
        0.99999999999999999999999999838,
        0.99999999999999999999999999919,
        0.99999999999999999999999999960,
        0.99999999999999999999999999980,
        0.99999999999999999999999999990,
        0.99999999999999999999999999995,
        0.99999999999999999999999999997,
        0.99999999999999999999999999999,
        0.99999999999999999999999999999,
        1.00000000000000000000000000000,
        1.00000000000000000000000000000,
        1.00000000000000000000000000000,
        };


        static int ETA_NEG_TABLE_NMAX = 99;
        //static int ETA_NEG_TABLE_SIZE = 50;
        static double[] eta_neg_int_table = {
         0.25000000000000000000000000000,   /* eta(-1) */
        -0.12500000000000000000000000000,   /* eta(-3) */
         0.25000000000000000000000000000,   /* ...      */
        -1.06250000000000000000000000000,
         7.75000000000000000000000000000,
        -86.3750000000000000000000000000,
         1365.25000000000000000000000000,
        -29049.0312500000000000000000000,
         800572.750000000000000000000000,
        -2.7741322625000000000000000000e+7,
         1.1805291302500000000000000000e+9,
        -6.0523980051687500000000000000e+10,
         3.6794167785377500000000000000e+12,
        -2.6170760990658387500000000000e+14,
         2.1531418140800295250000000000e+16,
        -2.0288775575173015930156250000e+18,
         2.1708009902623770590275000000e+20,
        -2.6173826968455814932120125000e+22,
         3.5324148876863877826668602500e+24,
        -5.3042033406864906641493838981e+26,
         8.8138218364311576767253114668e+28,
        -1.6128065107490778547354654864e+31,
         3.2355470001722734208527794569e+33,
        -7.0876727476537493198506645215e+35,
         1.6890450341293965779175629389e+38,
        -4.3639690731216831157655651358e+40,
         1.2185998827061261322605065672e+43,
        -3.6670584803153006180101262324e+45,
         1.1859898526302099104271449748e+48,
        -4.1120769493584015047981746438e+50,
         1.5249042436787620309090168687e+53,
        -6.0349693196941307074572991901e+55,
         2.5437161764210695823197691519e+58,
        -1.1396923802632287851130360170e+61,
         5.4180861064753979196802726455e+63,
        -2.7283654799994373847287197104e+66,
         1.4529750514918543238511171663e+69,
        -8.1705519371067450079777183386e+71,
         4.8445781606678367790247757259e+74,
        -3.0246694206649519336179448018e+77,
         1.9858807961690493054169047970e+80,
        -1.3694474620720086994386818232e+83,
         9.9070382984295807826303785989e+85,
        -7.5103780796592645925968460677e+88,
         5.9598418264260880840077992227e+91,
        -4.9455988887500020399263196307e+94,
         4.2873596927020241277675775935e+97,
        -3.8791952037716162900707994047e+100,
         3.6600317773156342245401829308e+103,
        -3.5978775704117283875784869570e+106    /* eta(-99)  */
        };


        #endregion

        #region static methods

        /* assumes s >= 0 and s != 1.0 */

        public static int RiemannZetaSgt0(double s, ref SpecialFunctionResult result)
        {
            if (s < 1.0)
            {
                SpecialFunctionResult c = new SpecialFunctionResult();
                ChebyshevSeries.ChebEvalE(zeta_xlt1_cs, 2.0 * s - 1.0, ref c);
                result.Val = c.Val / (s - 1.0);
                result.Err = c.Err / System.Math.Abs(s - 1.0) + PZMath_machine.PZMath_DBL_EPSILON * System.Math.Abs(result.Val);
                return PZMath_errno.PZMath_SUCCESS;
            }
            else if (s <= 20.0)
            {
                double x = (2.0 * s - 21.0) / 19.0;
                SpecialFunctionResult c = new SpecialFunctionResult();
                ChebyshevSeries.ChebEvalE(zeta_xgt1_cs, x, ref c);
                result.Val = c.Val / (s - 1.0);
                result.Err = c.Err / (s - 1.0) + PZMath_machine.PZMath_DBL_EPSILON * System.Math.Abs(result.Val);
                return PZMath_errno.PZMath_SUCCESS;
            }
            else
            {
                double f2 = 1.0 - System.Math.Pow(2.0, -s);
                double f3 = 1.0 - System.Math.Pow(3.0, -s);
                double f5 = 1.0 - System.Math.Pow(5.0, -s);
                double f7 = 1.0 - System.Math.Pow(7.0, -s);
                result.Val = 1.0 / (f2 * f3 * f5 * f7);
                result.Err = 3.0 * PZMath_machine.PZMath_DBL_EPSILON * System.Math.Abs(result.Val);
                return PZMath_errno.PZMath_SUCCESS;
            }
        } // RiemannZetaSgt0()

        public static int RiemannZeta1msSlt0(double s, ref SpecialFunctionResult result)
        {
            if (s > -19.0)
            {
                double x = (-19 - 2.0 * s) / 19.0;
                SpecialFunctionResult c = new SpecialFunctionResult();
                ChebyshevSeries.ChebEvalE(zeta_xgt1_cs, x, ref c);
                result.Val = c.Val / (-s);
                result.Err = c.Err / (-s) + PZMath_machine.PZMath_DBL_EPSILON * System.Math.Abs(result.Val);
                return PZMath_errno.PZMath_SUCCESS;
            }
            else
            {
                double f2 = 1.0 - System.Math.Pow(2.0, -(1.0 - s));
                double f3 = 1.0 - System.Math.Pow(3.0, -(1.0 - s));
                double f5 = 1.0 - System.Math.Pow(5.0, -(1.0 - s));
                double f7 = 1.0 - System.Math.Pow(7.0, -(1.0 - s));
                result.Val = 1.0 / (f2 * f3 * f5 * f7);
                result.Err = 3.0 * PZMath_machine.PZMath_DBL_EPSILON * System.Math.Abs(result.Val);
                return PZMath_errno.PZMath_SUCCESS;
            }
        } // RiemannZeta1msSlt0()

        /* works for 5 < s < 15*/
        public static int RiemannZetaMinus1IntermediateS(double s, ref SpecialFunctionResult result)
        {
            double t = (s - 10.0) / 5.0;
            SpecialFunctionResult c = new SpecialFunctionResult();
            ChebyshevSeries.ChebEvalE(zetam1_inter_cs, t, ref c);
            result.Val = System.Math.Exp(c.Val) + System.Math.Pow(2.0, -s);
            result.Err = (c.Err + 2.0 * PZMath_machine.PZMath_DBL_EPSILON) * result.Val;
            return PZMath_errno.PZMath_SUCCESS;
        } // RiemannZetaMinus1IntermediateS()

        /* assumes s is large and positive
         * write: zeta(s) - 1 = zeta(s) * (1 - 1/zeta(s))
         * and expand a few terms of the product formula to evaluate 1 - 1/zeta(s)
         *
         * works well for s > 15
         */
        public static int RiemannZetaMinus1LargeS(double s, ref SpecialFunctionResult result)
        {
            double a = System.Math.Pow(2.0, -s);
            double b = System.Math.Pow(3.0, -s);
            double c = System.Math.Pow(5.0, -s);
            double d = System.Math.Pow(7.0, -s);
            double e = System.Math.Pow(11.0, -s);
            double f = System.Math.Pow(13.0, -s);
            double t1 = a + b + c + d + e + f;
            double t2 = a * (b + c + d + e + f) + b * (c + d + e + f) + c * (d + e + f) + d * (e + f) + e * f;
            /*
            double t3 = a*(b*(c+d+e+f) + c*(d+e+f) + d*(e+f) + e*f) +
                        b*(c*(d+e+f) + d*(e+f) + e*f) +
                        c*(d*(e+f) + e*f) +
                        d*e*f;
            double t4 = a*(b*(c*(d + e + f) + d*(e + f) + e*f) + c*(d*(e+f) + e*f) + d*e*f) +
                        b*(c*(d*(e+f) + e*f) + d*e*f) +
                        c*d*e*f;
            double t5 = b*c*d*e*f + a*c*d*e*f+ a*b*d*e*f+ a*b*c*e*f+ a*b*c*d*f+ a*b*c*d*e;
            double t6 = a*b*c*d*e*f;
            */
            double numt = t1 - t2 /* + t3 - t4 + t5 - t6 */;
            double zeta = 1.0 / ((1.0 - a) * (1.0 - b) * (1.0 - c) * (1.0 - d) * (1.0 - e) * (1.0 - f));
            result.Val = numt * zeta;
            result.Err = (15.0 / s + 1.0) * 6.0 * PZMath_machine.PZMath_DBL_EPSILON * result.Val;
            return PZMath_errno.PZMath_SUCCESS;
        } // RiemannZetaMinus1LargeS()

        public static int HZetaE(double s, double q, ref SpecialFunctionResult result)
        {

            if (s <= 1.0 || q <= 0.0)
            {
                PZMath_errno.DomainError(ref result);
            }
            else
            {
                double max_bits = 54.0;
                double ln_term0 = -s * System.Math.Log(q);

                if (ln_term0 < PZMath_machine.PZMath_LOG_DBL_MIN + 1.0)
                {
                    PZMath_errno.UnderFlowError(ref result);
                }
                else if (ln_term0 > PZMath_machine.PZMath_LOG_DBL_MAX - 1.0)
                {
                    PZMath_errno.OverFlowError(ref result);
                }
                else if ((s > max_bits && q < 1.0) || (s > 0.5 * max_bits && q < 0.25))
                {
                    result.Val = System.Math.Pow(q, -s);
                    result.Err = 2.0 * PZMath_machine.PZMath_DBL_EPSILON * System.Math.Abs(result.Val);
                    return PZMath_errno.PZMath_SUCCESS;
                }
                else if (s > 0.5 * max_bits && q < 1.0)
                {
                    double p1 = System.Math.Pow(q, -s);
                    double p2 = System.Math.Pow(q / (1.0 + q), s);
                    double p3 = System.Math.Pow(q / (2.0 + q), s);
                    result.Val = p1 * (1.0 + p2 + p3);
                    result.Err = PZMath_machine.PZMath_DBL_EPSILON * (0.5 * s + 2.0) * System.Math.Abs(result.Val);
                    return PZMath_errno.PZMath_SUCCESS;
                }
                else
                {
                    /* Euler-Maclaurin summation formula 
                     * [Moshier, p. 400, with several typo corrections]
                     */
                    int jmax = 12;
                    int kmax = 10;
                    int j, k;
                    double pmax = System.Math.Pow(kmax + q, -s);
                    double scp = s;
                    double pcp = pmax / (kmax + q);
                    double ans = pmax * ((kmax + q) / (s - 1.0) + 0.5);

                    for (k = 0; k < kmax; k++)
                    {
                        ans += System.Math.Pow(k + q, -s);
                    }

                    for (j = 0; j <= jmax; j++)
                    {
                        double delta = hzeta_c[j + 1] * scp * pcp;
                        ans += delta;
                        if (System.Math.Abs(delta / ans) < 0.5 * PZMath_machine.PZMath_DBL_EPSILON)
                            break;
                        scp *= (s + 2 * j + 1) * (s + 2 * j + 2);
                        pcp /= (kmax + q) * (kmax + q);
                    }

                    result.Val = ans;
                    result.Err = 2.0 * (jmax + 1.0) * PZMath_machine.PZMath_DBL_EPSILON * System.Math.Abs(ans);
                    return PZMath_errno.PZMath_SUCCESS;
                }
            }

            return PZMath_errno.PZMath_SUCCESS;
        } // HZetaE()


        public static int ZetaE(double s, ref SpecialFunctionResult result)
        {

            if (s == 1.0)
            {
                PZMath_errno.DomainError(ref result);
            }
            else if (s >= 0.0)
            {
                return RiemannZetaSgt0(s, ref result);
            }
            else
            {
                /* reflection formula, [Abramowitz+Stegun, 23.2.5] */

                SpecialFunctionResult zeta_one_minus_s = new SpecialFunctionResult();
                int stat_zoms = RiemannZeta1msSlt0(s, ref zeta_one_minus_s);
                double sin_term = (System.Math.IEEERemainder(s, 2.0) == 0.0) ? 0.0 : System.Math.Sin(0.5 * PZMath.M_PI * System.Math.IEEERemainder(s, 4.0)) / PZMath.M_PI;

                if (sin_term == 0.0)
                {
                    result.Val = 0.0;
                    result.Err = 0.0;
                    return PZMath_errno.PZMath_SUCCESS;
                }
                else if (s > -170)
                {
                    /* We have to be careful about losing digits
                     * in calculating System.Math.Pow(2 Pi, s). The gamma
                     * function is fine because we were careful
                     * with that implementation.
                     * We keep an array of (2 Pi)^(10 n).
                     */
                    double[] twopi_pow = { 1.0,
                                     9.589560061550901348e+007,
                                     9.195966217409212684e+015,
                                     8.818527036583869903e+023,
                                     8.456579467173150313e+031,
                                     8.109487671573504384e+039,
                                     7.776641909496069036e+047,
                                     7.457457466828644277e+055,
                                     7.151373628461452286e+063,
                                     6.857852693272229709e+071,
                                     6.576379029540265771e+079,
                                     6.306458169130020789e+087,
                                     6.047615938853066678e+095,
                                     5.799397627482402614e+103,
                                     5.561367186955830005e+111,
                                     5.333106466365131227e+119,
                                     5.114214477385391780e+127,
                                     4.904306689854036836e+135
                                    };
                    int n = (int)System.Math.Floor((-s) / 10.0);
                    double fs = s + 10.0 * n;
                    double p = System.Math.Pow(2.0 * PZMath.M_PI, fs) / twopi_pow[n];

                    SpecialFunctionResult g = new SpecialFunctionResult();
                    int stat_g = Gamma.GammaE(1.0 - s, ref g);
                    result.Val = p * g.Val * sin_term * zeta_one_minus_s.Val;
                    result.Err = System.Math.Abs(p * g.Val * sin_term) * zeta_one_minus_s.Err;
                    result.Err += System.Math.Abs(p * sin_term * zeta_one_minus_s.Val) * g.Err;
                    result.Err += PZMath_machine.PZMath_DBL_EPSILON * (System.Math.Abs(s) + 2.0) * System.Math.Abs(result.Val);
                    return PZMath_errno.ErrorSelect2(stat_g, stat_zoms);
                }
                else
                {
                    /* The actual zeta function may or may not
                     * overflow here. But we have no easy way
                     * to calculate it when the prefactor(s)
                     * overflow. Trying to use System.Math.Log's and exp
                     * is no good because we loose a couple
                     * digits to the exp error amplification.
                     * When we gather a little more patience,
                     * we can implement something here. Until
                     * then just give up.
                     */
                    PZMath_errno.OverFlowError(ref result);
                }
            }
            return PZMath_errno.PZMath_SUCCESS;
        } // ZetaE()


        public static int ZetaIntE(int n, ref SpecialFunctionResult result)
        {

            if (n < 0)
            {
                if (!PZMath.IsOdd(n))
                {
                    result.Val = 0.0; /* exactly zero at even negative integers */
                    result.Err = 0.0;
                    return PZMath_errno.PZMath_SUCCESS;
                }
                else if (n > -ZETA_NEG_TABLE_NMAX)
                {
                    result.Val = zeta_neg_int_table[-(n + 1) / 2];
                    result.Err = 2.0 * PZMath_machine.PZMath_DBL_EPSILON * System.Math.Abs(result.Val);
                    return PZMath_errno.PZMath_SUCCESS;
                }
                else
                {
                    return ZetaE((double)n, ref result);
                }
            }
            else if (n == 1)
            {
                PZMath_errno.DomainError(ref result);
            }
            else if (n <= ZETA_POS_TABLE_NMAX)
            {
                result.Val = 1.0 + zetam1_pos_int_table[n];
                result.Err = 2.0 * PZMath_machine.PZMath_DBL_EPSILON * System.Math.Abs(result.Val);
                return PZMath_errno.PZMath_SUCCESS;
            }
            else
            {
                result.Val = 1.0;
                result.Err = PZMath_machine.PZMath_DBL_EPSILON;
                return PZMath_errno.PZMath_SUCCESS;
            }

            return PZMath_errno.PZMath_SUCCESS;
        } // ZetaIntE()

        public static int Zetam1E(double s, ref SpecialFunctionResult result)
        {
            if (s <= 5.0)
            {
                int stat = ZetaE(s, ref result);
                result.Val = result.Val - 1.0;
                return stat;
            }
            else if (s < 15.0)
            {
                return RiemannZetaMinus1IntermediateS(s, ref result);
            }
            else
            {
                return RiemannZetaMinus1LargeS(s, ref result);
            }
        } // Zetam1E()


        public static int Zetam1IntE(int n, ref SpecialFunctionResult result)
        {
            if (n < 0)
            {
                if (!PZMath.IsOdd(n))
                {
                    result.Val = -1.0; /* at even negative integers zetam1 == -1 since zeta is exactly zero */
                    result.Err = 0.0;
                    return PZMath_errno.PZMath_SUCCESS;
                }
                else if (n > -ZETA_NEG_TABLE_NMAX)
                {
                    result.Val = zeta_neg_int_table[-(n + 1) / 2] - 1.0;
                    result.Err = 2.0 * PZMath_machine.PZMath_DBL_EPSILON * System.Math.Abs(result.Val);
                    return PZMath_errno.PZMath_SUCCESS;
                }
                else
                {
                    /* could use gsl_sf_zetam1_e here but subtracting 1 makes no difference
                       for such large values, so go straight to the result */
                    return ZetaE((double)n, ref result);
                }
            }
            else if (n == 1)
            {
                PZMath_errno.DomainError(ref result);
            }
            else if (n <= ZETA_POS_TABLE_NMAX)
            {
                result.Val = zetam1_pos_int_table[n];
                result.Err = 2.0 * PZMath_machine.PZMath_DBL_EPSILON * System.Math.Abs(result.Val);
                return PZMath_errno.PZMath_SUCCESS;
            }
            else
            {
                return Zetam1E(n, ref result);
            }
            return PZMath_errno.PZMath_SUCCESS;
        } // Zetam1IntE()


        public static int EtaIntE(int n, ref SpecialFunctionResult result)
        {
            if (n > ETA_POS_TABLE_NMAX)
            {
                result.Val = 1.0;
                result.Err = PZMath_machine.PZMath_DBL_EPSILON;
                return PZMath_errno.PZMath_SUCCESS;
            }
            else if (n >= 0)
            {
                result.Val = eta_pos_int_table[n];
                result.Err = 2.0 * PZMath_machine.PZMath_DBL_EPSILON * System.Math.Abs(result.Val);
                return PZMath_errno.PZMath_SUCCESS;
            }
            else
            {
                /* n < 0 */

                if (!PZMath.IsOdd(n))
                {
                    /* exactly zero at even negative integers */
                    result.Val = 0.0;
                    result.Err = 0.0;
                    return PZMath_errno.PZMath_SUCCESS;
                }
                else if (n > -ETA_NEG_TABLE_NMAX)
                {
                    result.Val = eta_neg_int_table[-(n + 1) / 2];
                    result.Err = 2.0 * PZMath_machine.PZMath_DBL_EPSILON * System.Math.Abs(result.Val);
                    return PZMath_errno.PZMath_SUCCESS;
                }
                else
                {
                    SpecialFunctionResult z = new SpecialFunctionResult();
                    SpecialFunctionResult p = new SpecialFunctionResult();
                    int stat_z = ZetaIntE(n, ref z);
                    int stat_p = Exp.ExpE((1.0 - n) * PZMath.M_LN2, ref p);
                    int stat_m = Elementary.MultiplyE(-p.Val, z.Val, ref result);
                    result.Err = System.Math.Abs(p.Err * (PZMath.M_LN2 * (1.0 - n)) * z.Val) + z.Err * System.Math.Abs(p.Val);
                    result.Err += 2.0 * PZMath_machine.PZMath_DBL_EPSILON * System.Math.Abs(result.Val);
                    return PZMath_errno.ErrorSelect3(stat_m, stat_p, stat_z);
                }
            }
        } // EtaIntE()


        public static int EtaE(double s, ref SpecialFunctionResult result)
        {

            if (s > 100.0)
            {
                result.Val = 1.0;
                result.Err = PZMath_machine.PZMath_DBL_EPSILON;
                return PZMath_errno.PZMath_SUCCESS;
            }
            else if (System.Math.Abs(s - 1.0) < 10.0 * PZMath_machine.PZMath_ROOT5_FLT_EPSILON)
            {
                double del = s - 1.0;
                double c0 = PZMath.M_LN2;
                double c1 = PZMath.M_LN2 * (PZMath.M_EULER - 0.5 * PZMath.M_LN2);
                double c2 = -0.0326862962794492996;
                double c3 = 0.0015689917054155150;
                double c4 = 0.00074987242112047532;
                result.Val = c0 + del * (c1 + del * (c2 + del * (c3 + del * c4)));
                result.Err = 2.0 * PZMath_machine.PZMath_DBL_EPSILON * System.Math.Abs(result.Val);
                return PZMath_errno.PZMath_SUCCESS;
            }
            else
            {
                SpecialFunctionResult z = new SpecialFunctionResult();
                SpecialFunctionResult p = new SpecialFunctionResult();
                int stat_z = ZetaE(s, ref z);
                int stat_p = Exp.ExpE((1.0 - s) * PZMath.M_LN2, ref p);
                int stat_m = Elementary.MultiplyE(1.0 - p.Val, z.Val, ref result);
                result.Err = System.Math.Abs(p.Err * (PZMath.M_LN2 * (1.0 - s)) * z.Val) + z.Err * System.Math.Abs(p.Val);
                result.Err += 2.0 * PZMath_machine.PZMath_DBL_EPSILON * System.Math.Abs(result.Val);
                return PZMath_errno.ErrorSelect3(stat_m, stat_p, stat_z);
            }
        } // EtaE()
        #endregion
    }
}
