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


//---------------------
#ifndef MC_STRINGSH
#define MC_STRINGSH
//---------------------


#define ERROR_COD_MC_MSG							L"Error code No -> "							
#define ERROR_COD_MC_OFFSET							10

#define ERRCOD_MC_INVARGSTR							ERROR_COD_MC_OFFSET + 0
#define ERRMSG_MC_INVARGSTR							L"Invalid command line string."

#define ERRCOD_MC_BADINITFILE						ERROR_COD_MC_OFFSET + 4
#define ERRMSG_MC_BADINITFILE						L"Bad initialization. This operation needs a file and a suffix."

#define ERRCOD_MC_BADINITFREQ						ERROR_COD_MC_OFFSET + 5
#define ERRMSG_MC_BADINITFREQ						L"Bad initialization. This operation needs a frequency value (GHz)."

#define ERRCOD_MC_BADINITHPFILE						ERROR_COD_MC_OFFSET + 6
#define ERRMSG_MC_BADINITHPFILE						L"Bad initialization. This operation needs a Healpix file name."

#define ERRCOD_MC_INVALIDPARAM						ERROR_COD_MC_OFFSET + 7
#define ERRMSG_MC_INVALIDPARAM						L"Invalid parameter value or format."

#define ERRCOD_MC_HPNOTCOORDYS						ERROR_COD_MC_OFFSET + 8
#define ERRMSG_MC_HPNOTCOORDYS						L"This Healpix map does not have a COORDSYS keyword. Assuming values in the Parameter file. "

#define ERRCOD_MC_HPNOTTUNIT1						ERROR_COD_MC_OFFSET + 9
#define ERRMSG_MC_HPNOTTUNIT1						L"This Healpix map does not have an TUNIT1 keyword. Assuming values in the Parameter file. "

#define ERRCOD_MC_HPINVALIDCOORDYS					ERROR_COD_MC_OFFSET + 12
#define ERRMSG_MC_HPINVALIDCOORDYS					L"Invalid COORDSYS keyword. Possible values are <G|E|C>. Assuming values in the Parameter file. "

#define ERRCOD_MC_HPINVALIDTUNIT1					ERROR_COD_MC_OFFSET + 13
#define ERRMSG_MC_HPINVALIDTUNIT1					L"Invalid TUNIT1 keyword. Assuming values in the Parameter file. "

#define ERRCOD_MC_HPNSIDEMISMATCH					ERROR_COD_MC_OFFSET + 14
#define ERRMSG_MC_HPNSIDEMISMATCH					L"This set of maps requires a different Geom NSide. "

#define ERRCOD_MC_INVALIDANT						ERROR_COD_MC_OFFSET + 16
#define ERRMSG_MC_INVALIDANT						L"Invalid antenna. "

#define ERRCOD_MC_INVALIDFREQ						ERROR_COD_MC_OFFSET + 17
#define ERRMSG_MC_INVALIDFREQ						L"Invalid Frequency. "

#define ERRCOD_MC_MAPNOTFOUND						ERROR_COD_MC_OFFSET + 18
#define ERRMSG_MC_MAPNOTFOUND						L"Cannot find this map and this map MUST be available. "


#define	MC_IDHEALPIXSPIN							L"Healpix with spin"
#define	MC_CTELATITUDE								L"Constant latitude"
#define	MC_NONBLIND									L"Non Blind"

#define MC_PARAMSFILENAME							L"Parameters file name => "
#define MC_CONTEXTID								L"Context ID => "
#define MC_ARGSFLAG									L"Command flag => "
#define MC_IOFLAG									L"IO flag => "
#define MC_FIRSTPATCH								L"First patch => "
#define MC_LASTPATCH								L"Last patch => "
#define MC_STR1ARG									L"str 1st arg => "
#define MC_STR2ARG									L"str 2nd arg => "
#define MC_ARGPOINTING								L"Pointing => "
#define MC_PARAMETERS								L"\nParameter's values\n"
#define MC_YES										L"Yes"
#define MC_NO										L"No"

#define	HEALPIX_FREQ_WARNING						L"FREQ Keyword not defined. "
#define	HEALPIX_BEAM_NOSMOOTH						L"Map assumed NOT smoothed by the beam. "
#define	HEALPIX_TUNIT1_WARNING						L"TUNIT1 Keyword not defined. "
#define	HEALPIX_COORDSYS_WARNING					L"COORDSYS Keyword not defined. "
#define	HEALPIXMAP_WARNING							L" . Healpix map -> "
#define	HEALPIX_LOWRES_WARNING						L"This map has a  resolution lower than the nominal value. "
#define HEALPIX_CURRNSIDE							L"Current NSide -> "
#define HEALPIX_NOMINALNSIDE						L" .Nominal NSide -> "
#define WARMSG_MC_BADPIXELS							L"Bad pixels have been left in this patch: "
#define WARMSG_MC_BADPIXELSEND						L"bad pixels"
#define WARMSG_MC_NOPTGSFILES						L"No pointings files found. Is this the right input directory (output in the parameters file) ? -> "
#define WARMSG_MC_COORDSYSMIS						L"Coordsys mismatch. This may cause aliasing. Input coordys : "
#define WARMSG_MC_COORDSYSOUT						L" , Output coordsys : "
#define WARMSG_MC_CREATEDIR							L"Cannot create dir -> "
#define	PROCESSING_FILE								L"Processing file -> "
#define	ALREADYPROC_FILE							L"File has already been processed. "

#define HEALPIX_COORDSYS_MSG						L"COORDSYS = "
#define HEALPIX_COORDS_ECL							L"ECLIPTIC"
#define HEALPIX_COORDS_GAL							L"GALACTIC"
#define HEALPIX_COORDS_EQU							L"EQUATORIAL"
#define HEALPIX_UNITS_MSG							L"UNITS = "
#define HEALPIX_UNITSSCAL_MSG						L"UNITS Scaling = "

#define HEALPIX_UNITS_THERMO						L"Thermodynamic T / K "
#define HEALPIX_UNITS_ANTENN						L"Antenna T / K "
#define HEALPIX_UNITS_BRIGHT						L"Brightness MJys / Sr "

#define HEALPIX_GEOMNSIDE_MSG						L"Geom NSide -> "
#define HEALPIX_MAPNSIDE_MSG						L", Max NSide -> "

#define HEALPIX_BEAM_MSG							L"BEAM FWHM (arcmin) = "
#define HEALPIX_BEAM_CORRECTED_MSG					L"BEAM FWHM ***CORRECTED***(arcmin) = "
#define HEALPIX_ORGBEAM_MSG							L"FITS HEADER BEAM FWHM (arcmin) = "

#define HEALPIX_FREQ_MSG							L"Map Frequency (GHz) = "

#define MC_MASKREADINGBADPIX						L"Reading the bad pixels mask. "
#define MC_MASKREADINGBADPIXSMOOTH					L"Reading the bad pixels transparency mask. "

#define MC_MASKREADINGREJECT						L"Reading/updating the source rejection mask. "

#define MC_MASKCAHNGINGRES							L"Changing resolution of the mask. "
#define MC_MASKFILTERMAP							L"Making the Filtered map. "
#define MC_MASKFILTERMASK1							L"Enlarging the bad pixels mask. "
#define MC_MASKFILTERMASK2							L"Filtering the mask. "

#define MC_MASKFINALMAP								L"Making the final map. "
#define MC_MASKFORWARDTRANS							L"Harmonic forward transform. "
#define MC_MASKBACKTRANS							L"Harmonic inverse transform. "
#define MC_MASKFILTERING							L"Applying the filter. "
#define MC_MASKMASKALERT							L"WARNING: The COORDSYS of the map and the mask MUST be the SAME, usually GALACTIC. "

#ifdef PWS_IPAC
#define MC_MASKREMOVENAME_DEF	L""
#define MC_GALACTICCUT_DEF		0.0
#else
#define MC_MASKREMOVENAME_DEF	L""
#define MC_GALACTICCUT_DEF		5.0
#endif

#define MC_MASKREJECTNAME_DEF	L""
#define MC_NONBLINDPTGSFILE_DEF	L""
#define MC_PSCATFILENAME_DEF	L""

#define MC_MAPCUTTERID_DEF		1
#define MC_NSIDE_DEF			2048 
#define	MC_PATCHSZ_DEF			512 
//#define MC_PBORDER_DEF			-1
#define MC_PBORDER_DEF			0
#define MC_NLINEPTS_DEF			16 
#define MC_EPOCH_DEF			2000.0
#define MC_MASLENLAGTH_DEF		0.95
#define MC_RMSREJECTLEVEL_DEF	4.0	
#define MC_MASKRADIUSRATIO_DEF	-1.0
#define MC_SYNCID_DEF			0

#define MC_UNITS_DEF			1 // Thermo_T
#define MC_MAPTYPE_DEF			0 // Observation
#define MC_COORDSYS_DEF			2
#define MC_CONVFACTOR_DEF		1.0
#define MC_ANTFWHM_DEF			-1.0
#define MC_FREQ_DEF				-1.0 // invalid
#define MC_FILEPROCESSED_DEF	0
#define MC_PERCENTREJECT_DEF	0.50
#define MC_PTGSCOORDSYS_DEF		2
#define MC_MASKFWHM_DEF			66.0
#define MC_THRESSUB_DEF			4.8
#define MC_THRESMASK_DEF		15.0
#define MC_PTGSMAPNAME_DEF		L"PwS_PtgsMap"

#define	HEALPIX_BEAM			"BEAM"
#define	HEALPIX_TUNIT1			"TUNIT1"
#define	HEALPIX_FREQ			"FREQ"
#define HEALPIX_UNIT_THERMO		"CMB"
#define HEALPIX_UNIT_ANTENN		"ANT"
#define HEALPIX_UNIT_BRIGHT		"SR"
#define HEALPIX_UNIT_SCALE		"MK"
#define MC_PTCHMAXDIGNF			L"%04d_"
#define MC_SEPARATORS			L",;"
#define MC_LISTSEP				L","
#define MC_CUTTINGFORMSTR		L"Patch-> %04d, BP-> %02d, Th-> %6.2f, Phi-> %6.2f, Spin-> %+6.2f, Reject-> %6.2f\n"
#define MC_PIXELSFORMSTR		L"Pch-> %04d, Init_badpix -> %7.3f, End_badpix -> %7.3f, IPix-> %5d, Off -> %+7.3e, Rms -> %7.3e\n"
#define MC_REJECTPATCHSTR		L"This patch (Src index)-> %04d, has been rejected\n"
#define MC_IMGDRAWFORMSTR		L"Drawing patch number -> %04d\n"

#define MC_NBLDOUTCUTFORMSTR	L"Source in mask, Patch-> %04d, BP-> %02d, Th-> %6.2f, Phi-> %6.2f, Orginal src index -> %04d\n"

#if defined(HFIDMC) ||  defined(LFDMC)
#define	MT_OBSERVATION_DIR		L"_FreqMaps"
#define	MT_BACKCMB_DIR			L"_BackCMB"
#else
#define	MT_OBSERVATION_DIR		L"FreqMaps"
#define	MT_BACKCMB_DIR			L"BackCMB"
#endif

#define	MT_OBSERVATION_FILE			L"map"
#define	MT_BACKCMB_FILE				L"CMB"

#define MT_FILEPTGSDEFNAME			L"PwS_Pointings"
#define MT_FILEDETECTFNAME			L"PwS_Detections"

#define	MT_OBSERVATION_FILE_NAME		L"Temperature map"
#define MT_BACKGROUND_FILE_NAME			L"Background template "
#define MT_FILENAME_TXT					L"File type is -> "

#define MC_FINISHINGPROGRAM		L"Finishing program.\n\n\n"


#endif //MC_STRINGS

