/*
 *  This file is part of PowellSnakes.
 *
 *  PowellSnakes is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  PowellSnakes is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with PowellSnakes; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

/*
 *  PowellSnakes is being developed by Pedro Carvalho with the scientific advice of:
 *
 *  Dr.		Michael Hobson
 *  Prof.	Anthony Lasenby
 *	Dr.		Graca Rocha
 */

/*
 *  Copyright (C) 2005, 2012, Pedro Carvalho
 *  Author: Pedro Carvalho
 */



//----------------------------------
#ifndef DEFAULTSH
#define DEFAULTSH
//----------------------------------

#define GLOBAL_TWOSTEPSDETECTNOPRIOR_DEF		2
//Use two steps for detection when no prior planes

#define GLOBAL_NONBLINDDETECTION_DEF			0
// Non blind detection

#define GLOBAL_TWOSTEPSDETECTPRIORP_DEF			1
//Use one steps only for detection when no prior planes

#define GLOBAL_CATMERGE_AVGT_DEF				1
//Choose the brightest source always

#define GLOBAL_NOTALIGNEDOBJS_DEF		2
// 2 ; objects are not aligned, 2 use centroid

#define GLOBAL_NPRIORPLANES_DEF			0
// No prior templates

#if	defined(PWS_IPAC) || defined(LFIDPC)
#define GLOBAL_NFREQS_DEF				1
#else
#define GLOBAL_NFREQS_DEF				6
#endif

#if	defined(PWS_IPAC) || defined(LFIDPC)
#define GLOBAL_SZDETECTION_DEF			0
#else
#define GLOBAL_SZDETECTION_DEF			1
#endif

#define GLOBAL_PWSSEARCHPOS_DEF			0
// No ; Use symmetry properties

#define GLOBAL_USE2DFORMULA_DEF			1
// NO ; Search only Scale

#define GLOBAL_NSCALEBINS_SZ_DEF		2048
#define GLOBAL_NSCALEBINS_PS_DEF		1024
// N of source scale bins 

#define GLOBAL_USEBOUNDS_DEF			1
// Yes

#define GLOBAL_ASSESS_TYPE_DEF			1
// 1 -  Bayesian ; variable sigma threshold
// 0 -  Frequentist; fixed sigma threshold

#define GLOBAL_TEMPLATESIGMA_DEF		0.0
// The template noise RMS level; percentage

#define GLOBAL_BOUNDTOL_DEF				0.0005
// Bounding tolerance

#define GLOBAL_SIGMA_THRESHOLD_DEF		4.00
// Frequentist sigma

#define GLOBAL_SSUBLIMITGAUSS_DEF		-0.03
// Gaussianity threshold for source subtraction

#define DEFREJECTLEV_NPMAP				4.0

#define GLOBAL_MELTMAXDISTSZ_DEF		5.0
//SZ min maxdist 5 arcmin

#define GLOBAL_MELTMAXDISTSZMAXLIMIT_DEF	10.0
//SZ max limit maxdist 10.0 arcmin

#define GLOBAL_MELTMAXDISTPS_DEF		1.2
//PS ; times the antenna sigma rounded to nest pixel gives the distance

#define GLOBAL_MELTMAXDISTPSMAXLIMIT_DEF	-1.0
//PS masdist max limit ; Not used

#define GLOBAL_CATALOGUESIGMA_DEF		5.0
// do not filter

#define GLOBAL_APODIZEMAPS_DEF			0

#define GLOBAL_FLUXTHRESHOLD_SZ_DEF		5.0e-4
// 5e-4 arcmin^2

#define GLOBAL_FLUXTHRESHOLD_PS_DEF		200.0
// 200 mJys

#define GLOBAL_OUTPUTCOORDS_DEF				2
// Galactic

#define GLOBAL_HARDCONSTLEVEL_DEF			0
// 1 bit to remove 0 scale
// 0 loose

#define GLOBAL_OUTPUTDEGREES_SZ_DEF			1
// output degrees
#define GLOBAL_OUTPUTDEGREES_PS_DEF			0
// output rads
#define GLOBAL_GALACTICSIGMA_DEF			1.0
// *** changed ***
// = 0.0 -> real plane evaluation 
// != 0.0 -> filtering (likelihood) 

#define GLOBAL_OUTUNITS_DEF					0
// No units translation

#define GLOBAL_GALACTICCUT_PS_DEF			0.0
// These variables are reserve for future use (now always = 0.0)
#define GLOBAL_GALACTICCUT_SZ_DEF			0.0


// Linear filler, with fixed scales 0.0
#define GLOBAL_FILTERINGSCALES_DEF			L"0.0,33.0,0.0,0.80,0.90,1.01,1.14,1.28,1.44,1.61,1.81,2.04,2.29,2.58,2.89,3.25,3.66,4.11,4.62,5.19,5.84,6.56,7.38,8.29,9.32,10.47,11.77,13.23,14.88,16.72,18.79,21.12,23.74,26.69,30.00,33.75,37.97,42.71,48.05"

// Log filler 1.0; the first bin is 2nd element in the list
//#define GLOBAL_FILTERINGSCALES_DEF		L"1.0,0.8,6.0,8.0,9.0,11.0,13.0,14.0,16.0,17.0,19.0,21.0,23.0,25.0,29.0,33.0,38.0,42.0,46.0,50.0,58.0,67.0,83.0"

// Hardcoded scales 2.0; the scales follow (arcmin)
//#define GLOBAL_FILTERINGSCALES_DEF			L"2.0,0.0,0.80,0.90,1.01,1.14,1.28,1.44,1.61,1.81,2.04,2.29,2.58,2.89,3.25,3.66,4.11,4.62,5.19,5.84,6.56,7.38,8.29,9.32,10.47,11.77,13.23,14.88,16.72,18.79,21.12,23.74,26.69,30.00,33.72,37.90,42.60,47.89"
//#define GLOBAL_FILTERINGSCALES_DEF			L"2.0,0.0,0.80,0.90,1.01,1.14,1.28,1.44,1.61,1.81,2.04,2.29,2.58,2.89,3.25,3.66,4.11,4.62,5.19,5.84,6.56,7.38,8.29,9.32,10.47,11.77,13.23,14.88,16.72,18.79,21.12,23.74,26.69,30.00"
// Linear filler 3.0
//#define GLOBAL_FILTERINGSCALES_DEF			L"3.0,6.0,8.0,9.0,11.0,13.0,14.0,16.0,17.0,19.0,21.0,23.0,25.0,29.0,33.0,38.0,42.0,46.0,50.0,58.0,67.0,83.0"


#define GLOBAL_CACHESZ_DEF					1.5

#define GLOBAL_SZPROFILE_DEF				2
// Ngai07 profile


// Mparsec
#define GLOBAL_PRIORRADDISTMIN_DEF					0.001
#define GLOBAL_PRIORRADDISTMAX_DEF					1.000
// Solar masses
#define GLOBAL_PRIORMASSMIN_DEF						1.0
// default -> do not calibrate Cross-PowerSp
// make negative to calibrate
#define GLOBAL_PRIORMASSMAX_DEF						1.0
// default -> Countour plots are Y5r500 <-> ThetaS
// make negative to Beta <-> ThetaS

// redshift
#define GLOBAL_PRIORMAXZ_DEF						3.0
// KeV
#define GLOBAL_PRIORTEMPMIN_DEF						0.01
#define GLOBAL_PRIORTEMPMAX_DEF						20.0

#define GLOBAL_PRIORGASMASSRATIO_DEF				0.10
// Now sigma8 is a flag for the QA
// >= 0.0 do the QA (default)
// < 0.0 write the contours to disc
#define GLOBAL_PRIORSIGMA8_DEF						0.9

#define GLOBAL_USEPRESSSCHT_DEF						2

#define GLOBAL_MN_NLIVEPOINTS_DEF					200
#define GLOBAL_MN_NINDIVSAMPLES_DEF					20
#define GLOBAL_MN_FRACTOLEV_DEF						0.2
#define GLOBAL_MN_XTRAENLFACT_DEF					1.1

#define GLOBAL_OUTPURITY_DEF						10.0

#define GLOBAL_DEFAULTDATASET						L"unknown"

#define GLOBAL_REJECTMASKFNAME_DEF					L"PwS_RejectionMask"

#define GLOBAL_FLUXCALIBCTE_DEF						1.0

#define GLOBAL_PRIORFLUXMIN_PS_DEF					200.0
#define GLOBAL_PRIORFLUXMAX_PS_DEF					1000000.0
#define GLOBAL_PRIORFLUXEXP_PS_DEF					1.2

#define GLOBAL_PRIORFLUXMIN_SZ_DEF					0.0005
#define GLOBAL_PRIORFLUXMAX_SZ_DEF					0.20
#define GLOBAL_PRIORFLUXEXP_SZ_DEF					1.60

#define GLOBAL_SZPROFVARALPHAMIN_DEF				0.70
#define GLOBAL_SZPROFVARALPHAMAX_DEF				1.80
#define GLOBAL_SZPROFVARALPHABIN_DEF				-1
// Profile var *not* active
#define GLOBAL_SZPROFVARBETAMIN_DEF					4.0
#define GLOBAL_SZPROFVARBETAMAX_DEF					7.0
#define GLOBAL_SZPROFVARBETABIN_DEF					-1
// Profile var *not* active

#define GLOBAL_SZPROFVARC500MIN_DEF					1.0
#define GLOBAL_SZPROFVARC500MAX_DEF					1.0
#define GLOBAL_SZPROFVARC500BIN_DEF					1

#define GLOBAL_PRIORSRCSCALEMIN_PS_DEF				0.0001
#define GLOBAL_SRCMAXSCALE_PS_DEF					2.0
#define GLOBAL_SRCMAXSCALE_PSNB_DEF					2.0
#define GLOBAL_PRIORSRCSCALEEXP_PS_DEF				((double)(1.0/(GLOBAL_SRCMAXSCALE_PS_DEF * 2.0)))

#define GLOBAL_PRIORSRCSCALEMIN_SZ_DEF				1.30
#define GLOBAL_PRIORSRCSCALEMAX_SZ_DEF				40.0
#define GLOBAL_PRIORSRCSCALEEXP_SZ_DEF				0.20

#define GLOBAL_SRCMAXSCALE_PS_DEF					2.0
#define GLOBAL_SRCMAXSCALE_PSNB_DEF					2.0

#define GLOBAL_SRCMAXSCALE_SZ_DEF					50.0

#define GLOBAL_OUTPUTRADIUSCAL_SZ_DEF				1.0
#define GLOBAL_OUTPUTRADIUSCAL_PS_DEF				1.0

#define GLOBAL_JFESTIMATTYPE_DEF					1
// Jeffreys' estimator 0 - mode, 1 - mean

#define GLOBAL_PRIORAVSRCSZ_DEF						679.0
#define GLOBAL_PRIORAVSRCPS_DEF						358.0


#endif //DEFAULTSH

