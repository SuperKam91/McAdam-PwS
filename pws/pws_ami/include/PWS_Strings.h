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


//-------------------------------------
#ifndef PWS_STRINGSH
#define PWS_STRINGSH

#include "ZEUS_Strings.h"

//-------------------------------------
// Error mesages

#define PROGRAMSTRINGID								L"\"PwS_Ver_3.41\""
#define USAGEMSG1									L"\n\nPowellSnakes Ver. 3.41\n"

#define USAGEMSG2									L"CamPwSCore64 <ParamsFile> [OutputFile] [IO_flag] [FirstPatchNumber] [LastPatchNumber] [Catalogue0] .. [CatalogueN]\n\n"
#define PARAMSFILENAME								L"Parameters file name => "
#define OUTPUTFILENAME								L"Output file name => "
#define CONTEXTID									L"Context ID => "
#define IOFLAG										L"IO flag => "

#if	defined(WIN32) && !defined(HFIDMC) && !defined(LFIDPC) && !defined(AMI)
#define	PWSLOGGER									L"PWSLOGGER"
#define	PWSMAPCUTOPT								L"PWSMAPOPT"
#else
#define	PWSLOGGER									"PWSLOGGER"
#define	PWSMAPCUTOPT								"PWSMAPOPT"
#endif

#define	FINALCATSUFFIX								L"_cat"
#define FIRSTPATCH									L"First patch => "
#define LASTPATCH									L"Last patch => "
#define LASTPATCHN									L"last"
#define STARTPATCHNUMBER							L"Starting patch => "
#define REPORTINGPATCHNUMBER						L"\nReporting patch => "
#define REPPATCHNUMBERNOTVALID						L"\nThis patch %04d is not valid. Detection will not be attempted."
#define STARTREADINGMAPS							L"Start reading maps."
#define ENDREADINGMAPS								L"Finish reading maps."
#define STARTREADINGBACKGMAPS						L"Start reading background maps."
#define ENDREADINGBACKGMAPS							L"Finish reading backgorund maps."
#define VERYBRIGHTOBJFOUND							L"Very bright object found."
#define FINISHINGPROGRAM							L"Finishing program.\n\n\n"
#define CATALOGUENAME								L"Catalogue name => "
#define GLOBAL_YES									L"Yes"
#define GLOBAL_NO									L"No"
#define GLOBAL_PARAMETERS							L"\nParameter's values\n"
#define GLOBAL_ASSTYPEBAYS							L"Bayesian"
#define GLOBAL_ASSTYPEGLRT							L"GLRT"
#define GLOBAL_MERGTYPEHIGH							L"Highest detection"
#define GLOBAL_MERGTYPEAVER							L"Average"
#define GLOBAL_SZPROFILENFW							L"Generalised NFW"
#define GLOBAL_SZPROFILEBETA						L"Beta"
#define GLOBAL_SZPROFILEOTHER						L"Other"
#define GLOBAL_JFESTIMATORMEAN						L"Mean estimator"
#define GLOBAL_JFESTIMATORMODE						L"Mode estimator"
#define GLOBAL_MASSFUNCTPS							L"Press-Schechter"
#define GLOBAL_MASSFUNCTJK							L"Jenkins"
#define GLOBAL_STARTSAMPLING						L"\n Start sampling ...\n\n"

#ifdef WIN32
#define GLOBAL_FORMATSAMPLPAR						L"%ls -> Mean -> %8.4f, Mode -> %8.4f, StDev -> %8.4f, HPD_l -> %8.4f, HPD_h -> %8.4f"
#else
#define GLOBAL_FORMATSAMPLPAR						L"%S -> Mean -> %8.4f, Mode -> %8.4f, StDev -> %8.4f, HPD_l -> %8.4f, HPD_h -> %8.4f"
#endif

#define GLOBAL_FORMATDEGENPAR						L"Slope (*10^3) -> %8.4f, Correlation -> %8.4f, y0 (*10^3) -> %8.4f, Slope* (*10^3)-> %8.4f, y0* (*10^3)-> %8.4f"


#define GLOBAL_FORMATSAMPLDETECT					L"rho -> %g, Err -> %g, Evals -> %d\n"
#define GLOBAL_FORMATNPSTAT							L"\nSigma -> %f, X -> %.2f, Y -> %.2f"

#define MYERRORLOGLNTXT								L"This patch had an exception "
#define MYERRORSEPARATOR							L" => "

#define ERROR_COD_PWS_OFFSET						100000
#define PATCHNUMBERSTR								L". Patch number -> "
#define SRCINDEXSTR									L"Non-blind detection. Src index -> "



#define ERRCOD_PWS_NOTAINTEGER						ERROR_COD_PWS_OFFSET + 0
#define ERRMSG_PWS_NOTAINTEGER						L"Can't understand this as integer. "

#define ERRCOD_PWS_CANTFINDFILE						ERROR_COD_PWS_OFFSET + 1
#define ERRMSG_PWS_CANTFINDFILE						L"Can't find file. "

#define ERRCOD_PWS_INVCROSSPOWER					ERROR_COD_PWS_OFFSET + 2
#define ERRMSG_PWS_INVCROSSPOWER					L"The cross power matrix is not positive definitive or is singular."

#define ERRCOD_PWS_FREQNOTFOUND						ERROR_COD_PWS_OFFSET + 3
#define ERRMSG_PWS_FREQNOTFOUND						L"Inconsistency error. Frequency not found."

#define ERRCOD_PWS_PARAMOUTOFRANGE					ERROR_COD_PWS_OFFSET + 4
#define ERRMSG_PWS_PARAMOUTOFRANGE					L"Parameter out of range. "

#define ERRCOD_PWS_INVALIDSTRING					ERROR_COD_PWS_OFFSET + 5
#define ERRMSG_PWS_INVALIDSTRING					L"Invalid string"

#define ERRCOD_PWS_NOTINITIALISED					ERROR_COD_PWS_OFFSET + 6
#define ERRMSG_PWS_NOTINITIALISED					L"Object not initialised"

#define ERRCOD_PWS_PARATOOMANYPS					ERROR_COD_PWS_OFFSET + 7
#define ERRMSG_PWS_PARATOOMANYPS					L"Too many bright sources in patch. "

#define ERRCOD_PWS_INVALPRIORVAL					ERROR_COD_PWS_OFFSET + 8
#define ERRMSG_PWS_INVALPRIORVAL					L"Invalid prior value -> "

#define ERRCOD_PWS_INVALPRIORCONG					ERROR_COD_PWS_OFFSET + 9
#define ERRMSG_PWS_INVALPRIORCONG					L"Invalid non-blind configuration. "
#define ERRMSG_PWS_INVALPRIORDETECT					L"Invalid non-blind detection. "

#define ERRCOD_PWS_MPIERROR_						100000


// Zone identifiers
// All this identifiers will be preceeded by the Zone prefix

#define ZONE_FOREGOBJID					L"foreg_obj_id"
#define ZONE_FOREGOBJMFINITFACT			L"foreg_obj_mfilter_init_scale"
#define ZONE_ANTENNAPREFIX				L"antenna_"
#define ZONE_ANTENNAID					L"id_"
#define ZONE_ANTENNAFWHM				L"fwhm_"
#define ZONE_ANTENNAIDGAINFACT			L"gain_"
#define ZONE_ANTENNAIDPIXELNOISE		L"pixel_noise_"
#define ZONE_ANTENNAFREQ				L"freq_"
#define ZONE_AVGNOBJS					L"avg_n_objs"
#define ZONE_PRIORAMPLBETAPARAM			L"prior_ampl_beta"
#define ZONE_PRIORAMPLMINPARAM			L"prior_ampl_min"

#define	PATCH_OBSERVATION_FILE			L"map"
#define	PATCH_BACKCMB_FILE				L"CMB"


#if	defined(HFIDMC) || defined(LFIDPC)
#define	PATCH_OBSERVATION_DIR			L"_FreqMaps"
#define	PATCH_BACKCMB_DIR				L"_BackCMB"
#else
#define	PATCH_OBSERVATION_DIR			L"FreqMaps"
#define	PATCH_BACKCMB_DIR				L"BackCMB"
#endif

#define PHYMATH_REDSHIFT				L"Redshift"
#define PHYMATH_0CURVONLY				L"Sorry! Available for 0 curvature only"
#define OUTPUTFILELEADINGCHAR			L'_'


#endif //PWS_STRINGSH

