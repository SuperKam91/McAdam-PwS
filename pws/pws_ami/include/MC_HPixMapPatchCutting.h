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



//---------------------------------------------------------------------------
#ifndef MC_GLOBALSH
#define MC_GLOBALSH

#include <vector>
#include "pointing.h"
#include "ZEUS_InOut.h"
#include "ZEUS_WorkSpace.h"
#include "ZEUS_GlobalValuesStore.h"
#include "ZEUS_PhysicsMath.h"
#include "ZEUS_Strings.h"

#define HEALPIX_ATOM_PREC					float
#define MC_MAXNSIDE							2048						
#define MC_DISPLAYNSIDE						512						

//#define MULTI_THREAD

#ifdef MULTI_THREAD
typedef Zeus::MThGlobalVars					GlobalVars;
#else
typedef Zeus::SThGlobalVars					GlobalVars;
#endif

class GeneralMapCutter;
class ConstGalacticLatCutter;
class HealpixCutter;
class PatchFinder;
class PixExtractProc;

enum  ValidNSide {NSIDE_512=512,NSIDE_1024=1024,NSIDE_2048=2048};

struct PlanckMapsInfoType
{
	int			Freq_;
	double		EffFreq_;
	int			OperatNSide_;
	double		AntFWHM_;
	bool operator<(const PlanckMapsInfoType& sec)
	{return Freq_ < sec.Freq_;}
};

typedef std::vector<PlanckMapsInfoType>	PlanckMapsInfoCollType;


struct PixExtractProcCtorArgsType
{
	struct MapFileInfoType
	{
		int								FreqID_;
		Zeus::MapType					mapT_;
		int								MapCoordSys_;
		Zeus::PlanckUnitsValuesTranform::UnitsType	MapUnits_;
		double							UnitFactor_;
		double							Freq_;
		std::wstring					HP_FName_;
		std::wstring					PS_Catal_FName_;
		double							Ant_FWHM_;
		double							Thresh_sub_;
		double							Thresh_mask_;
		int								NSide_;
		int								Processed_;
	};

	typedef std::vector<MapFileInfoType>	MapFileInfoCollType;

	double							Epoch_;
	double							RmsRejectLevel_;
	double							MaskEnlagThreshold_;
	double							MaskSmoothFWHM_;
	double							MaskRadiusRatio_;
	std::wstring					HP_Dir_;
	std::wstring					Masks_Dir_;
	std::wstring					Masks_DirIn_;
	std::wstring					DataBuffer_;
	std::wstring					DirOut_;
	std::wstring					DirOutMaps_;
	std::wstring					MaskRejectName_;
	std::wstring					MaskRemoveName_;
	int								syncID_;
	int								Data2Buffer_;

	MapFileInfoCollType				FileCollection_;
};

struct ProcArgs
{
	int				ExecId_;
	int				NArgs_;
	int				EnvId;
	std::wstring	ParamFile_;
	wchar_t			flag_;
	pointing		ptg_;
	int				FrstPatch_;
	int				LastPatch_;
	std::wstring	fstString_;
	std::wstring	secString_;
};

#if	defined(WIN32) && !defined(HFIDMC) && !defined(LFIDPC)
void PreProcessArgs(int argc,wchar_t* argv[],ProcArgs& args);
#else
void PreProcessArgs(int argc,char* argv[],ProcArgs& args);
#endif


void				LauchHealpixPatchCutter(int mode, int endID,const std::wstring& fstString,const std::wstring& secString);
void				LauchHealpixPixExtractor(int endID,int OpId,long FPatch,long LPatch);
void				LauchHealpixPatchFinder(int endID,const pointing& ptg);
void				Extract1stMap(int endID);
void				MapCutterPringArgs(const ProcArgs& Args);
int					InitGlobalVars(int endID,const std::wstring& fname);
void				Do_ProcessArgs(std::vector<std::wstring>& ,ProcArgs& args);
void				PrintUsage(void);
void				ProcessEnvVars(ProcArgs& args);
std::wstring		GetFileTypeName(Zeus::MapType mapT);

GeneralMapCutter*	GetMapCutter(int MapCutterID,EnvIDsType InOutEnvID,const Zeus::PatchGeomType::HeaderType& Header,
								 const std::wstring& DirOut,
								 const std::wstring& NonBlindPtgsFile,
								 const std::wstring& MaskRejectFile,
								 const std::wstring& Chi2MaskFile,
								 const std::wstring& DirIn,coordsys PtgsCoordSys,double PercentReject,
								 const std::wstring& DataBuffer,int Object2buffer,double R500Ratio);


#endif // PCUT_GLOBALSH 

