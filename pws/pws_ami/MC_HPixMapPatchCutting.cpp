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

// HPixMapPatchCutting.cpp : Defines the entry point for the console application.
//
#include <iostream>
#include <memory>
#include "ZEUS_Exceptions.h"
#include "ZEUS_PhysicsMath.h"
#include "ZEUS_WorkSpace.h"
#include "FFTW_Traits.h"
#include "ZEUS_Strings.h"
#include "MC_HPixMapPatchCutting.h"
#include "MC_HealpixCutter.h"
#include "MC_PatchFinder.h"
#include "MC_PixExtractProc.h"

#ifdef _DEBUG
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#endif

#ifdef WIN32
//#include <Windows.h>
#include <conio.h>
#undef PWSMPI
#endif  //WIN32

#ifdef PWSMPI
#include <mpi.h>
#endif //PWSMPI

#if	defined(WIN32) && !defined(HFIDMC) && !defined(LFIDPC)
int _tmain(int argc, _TCHAR* argv[])
#else
int main(int argc, char* argv[])
#endif //WIN32
{

	if(argc < 2){
		PrintUsage();
		return -1;
	}

	ProcArgs			args;
	std::wstring		ExcptString;
#ifdef PWSMPI
	int					PwSMPIerr(MPI_SUCCESS);
	int					PwSnumprocs;
#endif //PWSMPI

	try{
#ifdef PWSMPI
		if(
			((PwSMPIerr = MPI_Init(&argc,&argv)) != MPI_SUCCESS) ||
			((PwSMPIerr = MPI_Comm_size(MPI_COMM_WORLD,&PwSnumprocs)) != MPI_SUCCESS)
			)
			throw Zeus::libMPIException(PwSMPIerr,L"MPI Error: Unable to initialise",L"MPI_Init");
#endif //PWSMPI
		PixExtractProc::SetMapStaticInfo();
		fftw::init_fftw();
		Zeus::ConManager::Initialize();
		PreProcessArgs(argc,argv,args);
		(Zeus::ConManager::Instance())->Initialise(argc,argv,args.ExecId_);
	}
#ifdef PWSMPI
 	catch(Zeus::libMPIException& err)
	{
		ExcptString = (std::wstring(PWSMPIINITEXPTMSG) + err.what_Xmsg());
		PwSMPIerr = err.GetErr();
	}
#endif // PWSMPI
 	catch(Zeus::libException& err)
	{
		ExcptString = (std::wstring(ZEUSPWSINITEXPTMSG) + err.what_Xmsg());
	}
	catch(...)
	{
		ExcptString = std::wstring(STDPWSINITEXPTMSG);
	}

	if(!(ExcptString.empty()))
	{
		Zeus::PrintErr2Console(ExcptString);
#ifdef PWSMPI
		if(PwSMPIerr!=MPI_SUCCESS)
			MPI_Abort(MPI_COMM_WORLD, PwSMPIerr);
		else
			MPI_Finalize();
#endif //PWSMPI
#ifdef WIN32
		wprintf(L"\n\nPress any key to terminate the program\n\n");
 		_getch();
#endif //WIN32
 		std::exit(-1);
	}
#ifdef PWSMPI
	int rk;
	MPI_Comm_rank(MPI_COMM_WORLD,&rk);
#endif
	try{
#ifdef PWSMPI
		if(!rk)
		{
#endif
			InitGlobalVars(args.EnvId,args.ParamFile_);
			if(args.flag_ == L'$')
			{
				std::vector<Zeus::VariableIdType>	ids;
				std::vector<Zeus::DBField>			values;

				ids.push_back(Zeus::VariableIdType(std::wstring(MC_COMMANDFLAG),
					Zeus::DBField(std::wstring(L""))));
				ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_FIRSTPATCH),
					Zeus::DBField(static_cast<int>(-1))));
				ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_LASTPATCH),
					Zeus::DBField(static_cast<int>(-1))));

				(GlobalVars::Instance())->GetVar(ids,values);

				std::vector<Zeus::DBField>::iterator v_ptr(values.begin());

				std::wstring	CommandFlagStr(*((*v_ptr++).GetPtr<std::wstring>()));
				
				args.FrstPatch_		= (*v_ptr++).Get<int>();
				args.LastPatch_		= (*v_ptr++).Get<int>();

				if(!(CommandFlagStr.empty()))
				{
					args.flag_	= Zeus::ToCase(CommandFlagStr.at(0),true);
				}
			}

			if( (args.flag_ != L'c') &&
				(args.flag_ != L'p') &&
				(args.flag_ != L'i') &&
				(args.flag_ != L't') &&
				(args.flag_ != L'm') &&
				(args.flag_ != L'd') &&
				(args.flag_ != L's') &&
				(args.flag_ != L'f') &&
				(args.flag_ != L'r') &&
				(args.flag_ != L'h') &&
				(args.flag_ != L'j')
			)
				throw Zeus::libException(ERRCOD_MC_INVARGSTR,ERRMSG_MC_INVARGSTR,L"main");

			MapCutterPringArgs(args);

			switch(args.flag_)
			{
			case L'p' :
				// Makes the pointing files
				LauchHealpixPatchCutter(0,args.EnvId,args.fstString_,args.secString_);
				break;
			case L's' :
				// Creates a Healpix map with the pointings
				LauchHealpixPixExtractor(args.EnvId,5,args.FrstPatch_,args.LastPatch_);
				break;
			case L'd' :
				// Creates a Healpix map with the detections
				LauchHealpixPatchCutter(2,args.EnvId,args.fstString_,args.secString_);
				break;
			case L'f' :
				// Converts a catalogue in the internal format to one of the may output formats
				// When using the parameter file fstString_ is read from nonblindptgs file
				// and secString_ from the rejectmask
				LauchHealpixPatchCutter(3,args.EnvId,args.fstString_,args.secString_);
				break;
			case L'h' :
				// Makes a CHI^2 catalogue
				// When using the parameter file fstString_ is read from nonblindptgs file
				// and secString_ from the rejectmask
				LauchHealpixPatchCutter(4,args.EnvId,args.fstString_,args.secString_);
				break;
			case L'c' :
				// 0, to extract the pixels from the healpix maps
				LauchHealpixPixExtractor(args.EnvId,0,args.FrstPatch_,args.LastPatch_);
				break;
			case L't' :
				// 1, to extract the fst column from the healpix files
				LauchHealpixPixExtractor(args.EnvId,1,args.FrstPatch_,args.LastPatch_);
				break;
			case L'm' :
				// 2, to make a mask
				LauchHealpixPixExtractor(args.EnvId,2,args.FrstPatch_,args.LastPatch_);
				break;
			case L'j' :
				// 3, join multiple masks
				LauchHealpixPixExtractor(args.EnvId,3,args.FrstPatch_,args.LastPatch_);
				break;
			case L'r' :
				// 4, Make maps/masks without point sources
				// New maps will be in the Masks directory
				LauchHealpixPixExtractor(args.EnvId,4,args.FrstPatch_,args.LastPatch_);
				break;
			case L'i' :
				LauchHealpixPatchFinder(args.EnvId,args.ptg_);
				break;
			default:
				PrintUsage();
			}
#ifdef PWSMPI
		}
#endif
	}
#ifdef PWSMPI
	catch(Zeus::libMPIException& err)
	{
		ExcptString = (std::wstring(PWSMPIINITEXPTMSG) + err.what_Xmsg());
		PwSMPIerr = err.GetErr();
	}
#endif //PWSMPI
	catch(Zeus::libException& err)
	{
		ExcptString		= (std::wstring(ZEUSPWSEXCPTMSG) + err.what_Xmsg());
	}
	catch(std::exception& err)
	{
		ExcptString		= (std::wstring(STANDARTEXCPTMSG) + Zeus::Achar2Wstr(err.what()));
	}
	catch(...)
	{
		ExcptString = std::wstring(STANDARTEXCPTMSG);
	}
	
	if(!(ExcptString.empty()))
	{
#ifdef PWSMPI
		if(PwSMPIerr!=MPI_SUCCESS)
			MPI_Abort(MPI_COMM_WORLD, PwSMPIerr);
#endif //PWSMPI
		(Zeus::ConManager::Instance())->PrintErr2Console(ExcptString);
		ExcptString.clear();
	}
	
	(Zeus::ConManager::Instance())->Flush();
	(Zeus::ConManager::Instance())->PrintStr2Console(MC_FINISHINGPROGRAM);
	(Zeus::ConManager::Instance())->CloseStreams();
#ifdef WIN32
	 wprintf(L"\n\nPress any key to terminate the program\n\n");
 	_getch();
#endif //WIN32
#ifdef PWSMPI
	MPI_Finalize();
#endif //PWSMPI
	std::exit(0);
	return 0;
}

int InitGlobalVars(int endID,const std::wstring& fname)
{
	(GlobalVars::Instance())->InitializeFromFile(fname,endID);
	return 1;
}

void PrintUsage(void)
{
#ifdef WIN32
		wprintf(L"\n\nMapCutter ver 3.60\n");
		wprintf(L"CamMapCutter64 <ParamsFile> [c | p | i | f | m | d | s | j | t] [envID] [IO_flag] [args] ...\n\n");
		wprintf(L"\n\nPress any key to terminate the program\n\n");
		_getch();
#else
		printf("\n\nMapCutter ver 3.60\n");
		printf("CamMapCutter64 <ParamsFile> [c | p | i | f | m | d | s | j | t] [endID] [IO_flag] [args] ...\n\n");
#endif //WIN32
}

#if	defined(WIN32) && !defined(HFIDMC) && !defined(LFIDPC)
void PreProcessArgs(int argc,wchar_t* argv[],ProcArgs& args)
#else
void PreProcessArgs(int argc,char* argv[],ProcArgs& args)
#endif
{
	args.NArgs_ = argc;
	std::vector<std::wstring> temp;
	for(int i=1;i<argc;++i)
#if	defined(WIN32) && !defined(HFIDMC) && !defined(LFIDPC)
	{temp.push_back(std::wstring(argv[i]));}
#else
	{temp.push_back(Zeus::Achar2Wstr(argv[i]));}
#endif
	Do_ProcessArgs(temp,args);
}

void	Do_ProcessArgs(std::vector<std::wstring>& strargs,ProcArgs& args)
{

	std::vector<std::wstring>::const_iterator piv(strargs.begin());
	std::vector<std::wstring>::const_iterator const end(strargs.end());
	args.flag_		= L'$';
	args.ptg_.theta	= -1.0;
	args.ptg_.phi	= -1.0;
	args.FrstPatch_ = -1;
	args.LastPatch_ = -1;
	args.ExecId_	= 0x03;
	args.fstString_.clear();
	args.secString_.clear();

#ifdef	HFIDMC
	args.EnvId			= 1;
#elif	LFIDPC
	args.EnvId			= 2;
#elif	MyDEBUG
	args.EnvId			= 0;
#else
	args.EnvId			= 3;
#endif

	args.ParamFile_ = *piv;

	ProcessEnvVars(args);

	if(++piv != end)
	{
		args.flag_		= Zeus::ToCase(piv->at(0),true);
	}

	if( (args.flag_ != L'c') &&
		(args.flag_ != L'p') &&
		(args.flag_ != L'i') &&
		(args.flag_ != L't') &&
		(args.flag_ != L'm') &&
		(args.flag_ != L'd') &&
		(args.flag_ != L's') &&
		(args.flag_ != L'f') &&
		(args.flag_ != L'r') &&
		(args.flag_ != L'h') &&
		(args.flag_ != L'j')
	)
		return;

	if(piv == end) return;
	else 	++piv;	

	if(piv == end) return;

	long	tEnvId;
	double	temp;

	if(!((*piv).empty()) && !(Zeus::get_number(*piv,tEnvId)))
	{args.ExecId_	= (int)tEnvId;}
	
	if(++piv == end) return;


	switch(args.flag_)
	{
	case	L'c':
		if(!((*piv).empty()) && !(Zeus::get_number(*piv,tEnvId)))
		{args.FrstPatch_	= (int)tEnvId;}
		if(++piv == end) return;
		if(!((*piv).empty()) && !(Zeus::get_number(*piv,tEnvId)))
		{args.LastPatch_	= (int)tEnvId;}
		break;
	case	L'i':
		if(!((*piv).empty()) && !(Zeus::get_number(*piv,temp)))
		{args.ptg_.theta	= ((90.0 - temp) / RAD2DEGREE);}
		if(++piv == end) return;
		if(!((*piv).empty()) && !(Zeus::get_number(*piv,temp)))
		{args.ptg_.phi	= (temp / RAD2DEGREE);}
		break;
	case	L'f':
	case	L'd':
	case	L'h':
		args.fstString_	= *piv;
		if(++piv == end) return;
		args.secString_	= *piv;
		break;
	case	L's':
		args.fstString_	= *piv;
		break;
	default:
		break;
	}
}

void	ProcessEnvVars(ProcArgs& args)
{
	std::wstring	buf(Zeus::GetEnvVar(PWSLOGGER));

	if(!buf.empty())
	{
		long	lNumber;

		if(Zeus::get_number(buf,lNumber))
			throw Zeus::libException(ERRCOD_PWS_NOTAINTEGER,ERRMSG_PWS_NOTAINTEGER,buf);
		args.ExecId_ = static_cast<int>(lNumber);
	}

	buf = Zeus::GetEnvVar(PWSMAPCUTOPT);

	if(!buf.empty())
	{
		args.flag_		= Zeus::ToCase(buf.at(0),true);
	}
}

void	LauchHealpixPatchCutter(int mode,int envID,const std::wstring& fstString,const std::wstring& secString)
{
	int									mapcutterid;
	std::vector<Zeus::VariableIdType>	ids;
	std::vector<Zeus::DBField>			values;
	Zeus::PatchGeomType::HeaderType		Header;
	std::wstring						DirOut;
	std::wstring						DirIn;
	std::wstring						MaskRejectFile;
	std::wstring						Chi2MaskFile;
	std::wstring						NonBldPtgsFile;
	int									PtgsCoordSys;
	double								PercentReject;
	double								R500Ratio;
	int									Objects2Buffer;
	std::wstring						DataBuffer;

	ids.push_back(Zeus::VariableIdType(std::wstring(MC_NSIDE),
		Zeus::DBField(static_cast<int>(MC_NSIDE_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(MC_PATCHSZ),
		Zeus::DBField(static_cast<int>(MC_PATCHSZ_DEF))));
	// This is now information for shuffling the patches 1 = shuffle
	ids.push_back(Zeus::VariableIdType(std::wstring(MC_PBORDER),
		Zeus::DBField(static_cast<int>(MC_PBORDER_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(MC_NLINEPTS),
		Zeus::DBField(static_cast<int>(MC_NLINEPTS_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(MC_MAPCUTTERID),
		Zeus::DBField(static_cast<int>(MC_MAPCUTTERID_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(MC_GALACTICCUT),
		Zeus::DBField(static_cast<double>(MC_GALACTICCUT_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_POINTINGS),
		Zeus::DBField()));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_DIRINMASKS),
		Zeus::DBField()));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_DIRBUFERDATA),
		Zeus::DBField(std::wstring(L""))));
	ids.push_back(Zeus::VariableIdType(std::wstring(MC_MASKREJECTNAME),
		Zeus::DBField(std::wstring(MC_MASKREJECTNAME_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(MC_MASKREMOVENAME),
		Zeus::DBField(std::wstring(MC_MASKREMOVENAME_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(MC_NONBLINDPTGSFILE),
		Zeus::DBField(std::wstring(MC_NONBLINDPTGSFILE_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(MC_PTGSCOORDSYS),
		Zeus::DBField(static_cast<int>(MC_PTGSCOORDSYS_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(MC_PERCENTREJECT),
		Zeus::DBField(static_cast<double>(MC_PERCENTREJECT_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_SZVIRIALRATIO),
		Zeus::DBField(static_cast<double>(GLOBAL_SZVIRIALRATIO_DEF))));

	(GlobalVars::Instance())->GetVar(ids,values);

	std::vector<Zeus::DBField>::iterator v_ptr(values.begin());

	(Zeus::ConManager::Instance())->PrintStr2Console(MC_PARAMETERS);

	Header.NSide_					= (*v_ptr++).Get<int>();
	Zeus::DumpScalarVariable(MC_NSIDE,Header.NSide_);
	Header.PtchSz_					= (*v_ptr++).Get<int>();
	if(Header.PtchSz_ < 0)
	{Header.PtchSz_ = 512;}
	Zeus::DumpScalarVariable(MC_PATCHSZ,Header.PtchSz_);
	
	Header.PtchBorder_				= -1;
	// Now Patch border contains Shuffle
	Header.Shuffle_					= ((*v_ptr++).Get<int>()!=0);
	Zeus::DumpScalarVariable(MC_PBORDER,Header.PtchBorder_);

	Header.NPatchesPerMainPix_		= (*v_ptr++).Get<int>();
	Zeus::DumpScalarVariable(MC_NLINEPTS,Header.NPatchesPerMainPix_);
	
	mapcutterid						= (*v_ptr++).Get<int>();
	Header.GalacticCut_				= (*v_ptr++).Get<double>();
	Zeus::DumpScalarVariable(MC_GALACTICCUT,Header.GalacticCut_);

	std::wstring					DummyGcc(*((*v_ptr++).GetPtr<std::wstring>()));
	std::wstring					DummyGcc1(*((*v_ptr++).GetPtr<std::wstring>()));
	std::wstring					DummyGcc2(*((*v_ptr++).GetPtr<std::wstring>()));

	DirOut							= Zeus::CorrectDir(Zeus::MakePath(DummyGcc));
	Zeus::DumpScalarVariable(GLOBAL_POINTINGS,DirOut);
	DirIn							= Zeus::CorrectDir(Zeus::MakePath(DummyGcc1));
	Zeus::DumpScalarVariable(GLOBAL_DIRINMASKS,DirIn);
	if(DummyGcc2.size()<2)
	{
		DataBuffer.clear();
		Objects2Buffer=0;
	}
	else
	{
		DataBuffer = Zeus::CorrectDir(Zeus::MakePath(DummyGcc2),1);
		Objects2Buffer=1;
	}
	Zeus::DumpScalarVariable(GLOBAL_DIRBUFERDATA,DataBuffer);

	MaskRejectFile					= *((*v_ptr++).GetPtr<std::wstring>());
	Zeus::DumpScalarVariable(MC_MASKREJECTNAME,MaskRejectFile);
	Chi2MaskFile					= *((*v_ptr++).GetPtr<std::wstring>());
	Zeus::DumpScalarVariable(MC_MASKREMOVENAME,Chi2MaskFile);
	NonBldPtgsFile					= *((*v_ptr++).GetPtr<std::wstring>());
	Zeus::DumpScalarVariable(MC_NONBLINDPTGSFILE,NonBldPtgsFile);
	PtgsCoordSys					= (*v_ptr++).Get<int>();
	PercentReject					= ((*v_ptr++).Get<double>());
	Zeus::DumpScalarVariable(MC_PERCENTREJECT,PercentReject);
	PercentReject					/= 100.0;
	R500Ratio = ((*v_ptr++).Get<double>());
	Zeus::DumpScalarVariable(GLOBAL_SZVIRIALRATIO, R500Ratio);

	if((PtgsCoordSys<0) || (PtgsCoordSys>2))
		PtgsCoordSys = 2; // default Galactic
	if(mapcutterid<0)
		mapcutterid = 1; // const latitude
	Zeus::DumpScalarVariable(MC_PTGSCOORDSYS,PtgsCoordSys?(PtgsCoordSys==2?std::wstring(L"G"):std::wstring(L"J")):std::wstring(L"E"));
	Zeus::DumpScalarVariable(MC_MAPCUTTERID,mapcutterid?(mapcutterid>=2?std::wstring(MC_NONBLIND):std::wstring(MC_CTELATITUDE)):std::wstring(MC_IDHEALPIXSPIN));
	
	std::auto_ptr<GeneralMapCutter> hp_cutter(GetMapCutter(mapcutterid,static_cast<EnvIDsType>(envID),Header,DirOut,
		NonBldPtgsFile,MaskRejectFile,Chi2MaskFile,DirIn,static_cast<coordsys>(PtgsCoordSys),
		PercentReject,DataBuffer,Objects2Buffer,R500Ratio));

	if(!mode)
	{
		hp_cutter->Initialize();
	}

	switch(mode)
	{
	case	2:
		hp_cutter->CreatHealpixMapWithDetections(fstString.empty()?NonBldPtgsFile:fstString,secString.empty()?MaskRejectFile:secString);
		break;
	case	3:
		hp_cutter->ConvertCatalogue(fstString.empty()?NonBldPtgsFile:fstString,secString.empty()?(Chi2MaskFile.empty()?MaskRejectFile:Chi2MaskFile):secString);
		break;
	case	4:
		hp_cutter->MakeChi2Catalogue(fstString.empty()?NonBldPtgsFile:fstString,secString.empty()?MaskRejectFile:secString);
		break;
	default:
		hp_cutter->MakePointingsFiles();
	}

}

void	LauchHealpixPixExtractor(int envID,int OpId,long FPatch,long LPatch)
{
	int	NFiles;
	std::vector<Zeus::VariableIdType>	ids;
	std::vector<Zeus::DBField>			values;
	PixExtractProcCtorArgsType			temp;
	PixExtractProcCtorArgsType::MapFileInfoType	FileLineTemp;
	std::wstring						fileID;
	wchar_t								buffer[DOUBLETXTMAXSZ];

	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_POINTINGS),
		Zeus::DBField()));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_DIRIN),
		Zeus::DBField()));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_DIRINMAPS),
		Zeus::DBField()));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_DIRINMASKS),
		Zeus::DBField()));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_DIRINMASKSIN),
		Zeus::DBField(std::wstring(L""))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_DIRBUFERDATA),
		Zeus::DBField(std::wstring(L""))));

	ids.push_back(Zeus::VariableIdType(std::wstring(MC_MASKREJECTNAME),
		Zeus::DBField(std::wstring(MC_MASKREJECTNAME_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(MC_MASKREMOVENAME),
		Zeus::DBField(std::wstring(MC_MASKREMOVENAME_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(MC_MASLENLAGTH),
		Zeus::DBField(static_cast<double>(MC_MASLENLAGTH_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(MC_MASKFWHM),
		Zeus::DBField(static_cast<double>(MC_MASKFWHM_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(MC_EPOCH),
		Zeus::DBField(static_cast<double>(MC_EPOCH_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(MC_RMSREJECTLEVEL),
		Zeus::DBField(static_cast<double>(MC_RMSREJECTLEVEL_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(MC_NOFFILES),
		Zeus::DBField()));
	ids.push_back(Zeus::VariableIdType(std::wstring(MC_SYNCID),
		Zeus::DBField(static_cast<int>(MC_SYNCID_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(MC_MASKRADIUSRATIO),
		Zeus::DBField(static_cast<double>(MC_MASKRADIUSRATIO_DEF))));

	(GlobalVars::Instance())->GetVar(ids,values);

	std::vector<Zeus::DBField>::iterator v_ptr(values.begin());

	std::wstring	DummyGcc1(*((*v_ptr++).GetPtr<std::wstring>()));
	std::wstring	DummyGcc2(*((*v_ptr++).GetPtr<std::wstring>()));
	std::wstring	DummyGcc3(*((*v_ptr++).GetPtr<std::wstring>()));
	std::wstring	DummyGcc4(*((*v_ptr++).GetPtr<std::wstring>()));
	std::wstring	DummyGcc5(*((*v_ptr++).GetPtr<std::wstring>()));
	std::wstring	DummyGcc6(*((*v_ptr++).GetPtr<std::wstring>()));

	(Zeus::ConManager::Instance())->PrintStr2Console(MC_PARAMETERS);

	temp.DirOut_			= Zeus::CorrectDir(Zeus::MakePath(DummyGcc1));
	Zeus::DumpScalarVariable(GLOBAL_POINTINGS,temp.DirOut_);
	temp.DirOutMaps_		= Zeus::CorrectDir(Zeus::MakePath(DummyGcc2));
	Zeus::DumpScalarVariable(GLOBAL_DIRIN,temp.DirOutMaps_);
	temp.HP_Dir_			= Zeus::CorrectDir(Zeus::MakePath(DummyGcc3));
	Zeus::DumpScalarVariable(GLOBAL_DIRINMAPS,temp.HP_Dir_);
	temp.Masks_Dir_			= Zeus::CorrectDir(Zeus::MakePath(DummyGcc4));
	Zeus::DumpScalarVariable(GLOBAL_DIRINMASKS,temp.Masks_Dir_);

	if(DummyGcc5.size()<2)
	{temp.Masks_DirIn_ = temp.Masks_Dir_;}
	else
	{temp.Masks_DirIn_ = Zeus::CorrectDir(Zeus::MakePath(DummyGcc5));}
	Zeus::DumpScalarVariable(GLOBAL_DIRINMASKSIN,temp.Masks_DirIn_);

	if(DummyGcc6.size()<2)
	{temp.DataBuffer_.clear();temp.Data2Buffer_ = 0;}
	else
	{temp.DataBuffer_ = Zeus::CorrectDir(Zeus::MakePath(DummyGcc6),1);temp.Data2Buffer_ = 1;}
	Zeus::DumpScalarVariable(GLOBAL_DIRBUFERDATA,temp.DataBuffer_);


	temp.MaskRejectName_	= *((*v_ptr++).GetPtr<std::wstring>());
	Zeus::DumpScalarVariable(MC_MASKREJECTNAME,temp.MaskRejectName_);
	temp.MaskRemoveName_	= *((*v_ptr++).GetPtr<std::wstring>());
	Zeus::DumpScalarVariable(MC_MASKREMOVENAME,temp.MaskRemoveName_);
	temp.MaskEnlagThreshold_= static_cast<float>((*v_ptr++).Get<double>());
	Zeus::DumpScalarVariable(MC_MASLENLAGTH,temp.MaskEnlagThreshold_);
	temp.MaskSmoothFWHM_	= (*v_ptr++).Get<double>();
	Zeus::DumpScalarVariable(MC_MASKFWHM,temp.MaskSmoothFWHM_);
	temp.Epoch_				= (*v_ptr++).Get<double>();
	Zeus::DumpScalarVariable(MC_EPOCH,temp.Epoch_);
	temp.RmsRejectLevel_	= (*v_ptr++).Get<double>();
	Zeus::DumpScalarVariable(MC_RMSREJECTLEVEL,temp.RmsRejectLevel_);
	NFiles					= (*v_ptr++).Get<int>();
	Zeus::DumpScalarVariable(MC_NOFFILES,NFiles);
	temp.syncID_			= (*v_ptr++).Get<int>();
	Zeus::DumpScalarVariable(MC_SYNCID,temp.syncID_);
	temp.MaskRadiusRatio_	= (*v_ptr++).Get<double>();
	Zeus::DumpScalarVariable(MC_MASKRADIUSRATIO,temp.MaskRadiusRatio_);

	Zeus::DBField	tDirMasksIn;
	(GlobalVars::Instance())->GetVar(Zeus::VariableIdType(std::wstring(GLOBAL_DIRINMASKSIN),
		Zeus::DBField(temp.Masks_Dir_)),tDirMasksIn);

	temp.Masks_DirIn_		= Zeus::CorrectDir(Zeus::MakePath((*(tDirMasksIn.GetPtr<std::wstring>()))));
	Zeus::DumpScalarVariable(GLOBAL_DIRINMASKSIN,temp.Masks_DirIn_);

	ids.clear();

	for(int i=0;i!=NFiles;++i)
	{
		fileID = std::wstring(MC_FILEIDSTR);
		
		fileID += Zeus::PutNumber2Txt(i);
		ids.push_back(Zeus::VariableIdType(fileID + std::wstring(MC_FILEFREQIDSTR),
			Zeus::DBField()));
		ids.push_back(Zeus::VariableIdType(fileID + std::wstring(MC_FILEMAPTYPESTR),
			Zeus::DBField(static_cast<int>(MC_MAPTYPE_DEF))));
		ids.push_back(Zeus::VariableIdType(fileID + std::wstring(MC_FILEUNITSSTR),
			Zeus::DBField(static_cast<int>(MC_UNITS_DEF))));
		ids.push_back(Zeus::VariableIdType(fileID + std::wstring(MC_FILECOORDSYSSTR),
			Zeus::DBField(static_cast<int>(MC_COORDSYS_DEF))));
		ids.push_back(Zeus::VariableIdType(fileID + std::wstring(MC_FILEUNITCONVFACTORSTR),
			Zeus::DBField(static_cast<double>(MC_CONVFACTOR_DEF))));
		ids.push_back(Zeus::VariableIdType(fileID + std::wstring(MC_FILEFREQSTR),
			Zeus::DBField(static_cast<double>(MC_FREQ_DEF))));
		ids.push_back(Zeus::VariableIdType(fileID + std::wstring(MC_FILEMAPFILENAMESTR),
			Zeus::DBField()));
		ids.push_back(Zeus::VariableIdType(fileID + std::wstring(MC_FILEPSCATFILENAMESTR),
			Zeus::DBField(std::wstring(MC_PSCATFILENAME_DEF))));
		ids.push_back(Zeus::VariableIdType(fileID + std::wstring(MC_FILEANTFWHMSTR),
			Zeus::DBField(static_cast<double>(MC_ANTFWHM_DEF))));
		ids.push_back(Zeus::VariableIdType(fileID + std::wstring(MC_THRESSUBSTR),
			Zeus::DBField(static_cast<double>(MC_THRESSUB_DEF))));
		ids.push_back(Zeus::VariableIdType(fileID + std::wstring(MC_THRESMASKSTR),
			Zeus::DBField(static_cast<double>(MC_THRESMASK_DEF))));
		ids.push_back(Zeus::VariableIdType(fileID + std::wstring(MC_FILEPROCESSEDSTR),
			Zeus::DBField(static_cast<int>(MC_FILEPROCESSED_DEF))));
	}

	(GlobalVars::Instance())->GetVar(ids,values);

	v_ptr = values.begin();
	temp.FileCollection_.clear();
	
	(Zeus::ConManager::Instance())->PrintStr2Console(L".\n");

	for(int i=0;i!=NFiles;++i)
	{
		fileID = std::wstring(MC_FILEIDSTR);
		fileID += Zeus::PutNumber2Txt(i);
		FileLineTemp.FreqID_		= (*v_ptr++).Get<int>();
		Zeus::DumpScalarVariable(fileID + std::wstring(MC_FILEFREQIDSTR),FileLineTemp.FreqID_);
		FileLineTemp.mapT_			= static_cast<Zeus::MapType>((*v_ptr++).Get<int>());
		Zeus::DumpScalarVariable(fileID + std::wstring(MC_FILEMAPTYPESTR),GetFileTypeName(FileLineTemp.mapT_));
		FileLineTemp.MapUnits_		= static_cast<Zeus::PlanckUnitsValuesTranform::UnitsType>((*v_ptr++).Get<int>());
		Zeus::DumpScalarVariable(fileID + std::wstring(MC_FILEUNITSSTR),FileLineTemp.MapUnits_?(FileLineTemp.MapUnits_==Zeus::PlanckUnitsValuesTranform::BRIGHTNESS?std::wstring(L"BRIGHTNESS"):std::wstring(L"THERMO_T")):std::wstring(L"ANTENNA_T"));
		FileLineTemp.MapCoordSys_	= (*v_ptr++).Get<int>();
		Zeus::DumpScalarVariable(fileID + std::wstring(MC_FILECOORDSYSSTR),FileLineTemp.MapCoordSys_?(FileLineTemp.MapCoordSys_==2?std::wstring(L"G"):std::wstring(L"J")):std::wstring(L"E"));
		FileLineTemp.UnitFactor_	= (*v_ptr++).Get<double>();
		Zeus::DumpScalarVariable(fileID + std::wstring(MC_FILEUNITCONVFACTORSTR),FileLineTemp.UnitFactor_);
		FileLineTemp.Freq_			= (*v_ptr++).Get<double>();
		Zeus::DumpScalarVariable(fileID + std::wstring(MC_FILEFREQSTR),FileLineTemp.Freq_);
		FileLineTemp.HP_FName_		= std::wstring(*((*v_ptr++).GetPtr<std::wstring>()));
		Zeus::DumpScalarVariable(fileID + std::wstring(MC_FILEMAPFILENAMESTR),FileLineTemp.HP_FName_);
		FileLineTemp.PS_Catal_FName_		= std::wstring(*((*v_ptr++).GetPtr<std::wstring>()));
		Zeus::DumpScalarVariable(fileID + std::wstring(MC_FILEPSCATFILENAMESTR),FileLineTemp.PS_Catal_FName_);
		FileLineTemp.Ant_FWHM_		= (*v_ptr++).Get<double>();
		Zeus::DumpScalarVariable(fileID + std::wstring(MC_FILEANTFWHMSTR),FileLineTemp.Ant_FWHM_);
		FileLineTemp.Thresh_sub_		= (*v_ptr++).Get<double>();
		Zeus::DumpScalarVariable(fileID + std::wstring(MC_THRESSUBSTR),FileLineTemp.Thresh_sub_);
		FileLineTemp.Thresh_mask_		= (*v_ptr++).Get<double>();
		Zeus::DumpScalarVariable(fileID + std::wstring(MC_THRESMASKSTR),FileLineTemp.Thresh_mask_);
		FileLineTemp.Processed_		= (*v_ptr++).Get<int>();
		Zeus::DumpScalarVariable(fileID + std::wstring(MC_FILEPROCESSEDSTR),FileLineTemp.Processed_?std::wstring(GLOBAL_YES):std::wstring(GLOBAL_NO));
		temp.FileCollection_.push_back(FileLineTemp);

		(Zeus::ConManager::Instance())->PrintStr2Console(L".\n");
	}

	Zeus::PlaneBoundsType::SetDefaultUseBounds(Zeus::UB_NOUSE);
	Zeus::PlaneBoundsType::SetDefaultBoundingFactor(static_cast<float>(0.0001));

	std::auto_ptr<PixExtractProc> ptPixExtract(new PixExtractProc(static_cast<EnvIDsType>(envID),temp));

	switch(OpId)
	{
	case	1:
		ptPixExtract->Extract1stMap(temp.MaskRejectName_);
		break;
	case	2:
		ptPixExtract->DoPixExtractionAllMaps(2,0,0);
		break;
	case	3:
		ptPixExtract->JoinManyMasks();
		break;
	case	4:
		ptPixExtract->SubSourcesDoProcUpAllMaps();
		break;
	case	5:
		ptPixExtract->CreatHealpixMapWithPtgs(temp.MaskRejectName_);
		break;
	default:
		ptPixExtract->DoPixExtractionAllMaps(0,FPatch,LPatch);
		break;
	}
}

std::wstring	GetFileTypeName(Zeus::MapType mapT)
{
	std::wstring temp;

	switch(mapT)
	{
	case Zeus::MAPTYPE_OBSERVATION:
		temp = std::wstring(MT_OBSERVATION_FILE_NAME);
		break;
	case Zeus::MAPTYPE_BACKGROUND:
		temp = std::wstring(MT_BACKGROUND_FILE_NAME);
		break;
	default:
		temp = std::wstring(L"<INVALID>");
	}

	return temp;
}

void	LauchHealpixPatchFinder(int envID,const pointing& ptg)
{

	if((ptg.theta < 0) || (ptg.theta > PI) || (ptg.phi < 0) || (ptg.phi > PITIMES2))
		throw Zeus::libException(ERRCOD_MC_INVARGSTR,ERRMSG_MC_INVARGSTR,L"LauchHealpixPatchFinder");

	std::vector<Zeus::VariableIdType>	ids;
	std::vector<Zeus::DBField>			values;
	std::wstring						DirOut;

	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_POINTINGS),
		Zeus::DBField()));

	(GlobalVars::Instance())->GetVar(ids,values);

	std::vector<Zeus::DBField>::iterator v_ptr(values.begin());

	(Zeus::ConManager::Instance())->PrintStr2Console(MC_PARAMETERS);

	std::wstring		DummyGcc1(*((*v_ptr++).GetPtr<std::wstring>()));
	DirOut				= Zeus::CorrectDir(Zeus::MakePath(DummyGcc1));
	Zeus::DumpScalarVariable(GLOBAL_POINTINGS,DirOut);
	std::auto_ptr<PatchFinder> ptFinder(new PatchFinder(static_cast<EnvIDsType>(envID),DirOut));
	ptFinder->FindPacthesFromPtg(ptg);
}
//
GeneralMapCutter*	GetMapCutter(int MapCutterID,EnvIDsType InOutEnvID,const Zeus::PatchGeomType::HeaderType& Header,const std::wstring& DirOut,
								 const std::wstring& NonBlindPtgsFile,
								 const std::wstring& MaskRejectFile,const std::wstring& Chi2MaskFile,
								 const std::wstring& DirIn,coordsys PtgsCoordSys,double PercentReject,
								 const std::wstring& DataBuffer, int Object2buffer, double R500Ratio)
{
	switch(MapCutterID)
	{
	case	0:
		return new HealpixCutter(static_cast<EnvIDsType>(InOutEnvID),Header,DirOut,MaskRejectFile,Chi2MaskFile,DirIn,
									PtgsCoordSys,PercentReject,DataBuffer,Object2buffer);
	case	1:
		return new ConstGalacticLatCutter(static_cast<EnvIDsType>(InOutEnvID),Header,DirOut,MaskRejectFile,Chi2MaskFile,DirIn,
											PtgsCoordSys,PercentReject,DataBuffer,Object2buffer);
	default:
		return new NonBlindCutter(static_cast<EnvIDsType>(InOutEnvID), R500Ratio, Header, DirOut, NonBlindPtgsFile, MaskRejectFile, Chi2MaskFile, DirIn,
									PtgsCoordSys,PercentReject,DataBuffer,Object2buffer);
	}
}
//
void				MapCutterPringArgs(const ProcArgs& Args)
{
	std::wstring temp(MC_PARAMSFILENAME);
	(Zeus::ConManager::Instance())->PrintStr2Console(temp + Args.ParamFile_);
	temp = MC_CONTEXTID;
	(Zeus::ConManager::Instance())->PrintStr2Console(temp + Zeus::PutNumber2Txt(Args.EnvId));
	temp = MC_ARGSFLAG;

	wchar_t tbuffer[5];
	tbuffer[0] = Args.flag_;tbuffer[1] = ((wchar_t)0);
	(Zeus::ConManager::Instance())->PrintStr2Console(temp + std::wstring(tbuffer));
	temp = MC_IOFLAG;
	(Zeus::ConManager::Instance())->PrintStr2Console(temp + Zeus::PutNumber2Txt(Args.ExecId_));
	temp = MC_FIRSTPATCH;
	if(Args.FrstPatch_ < 0)
	{
		(Zeus::ConManager::Instance())->PrintStr2Console(temp + Zeus::PutNumber2Txt(0));
	}
	else
	{
		(Zeus::ConManager::Instance())->PrintStr2Console(temp + Zeus::PutNumber2Txt(Args.FrstPatch_));
	}
	temp = MC_LASTPATCH;
	if((Args.LastPatch_ < 0) )
	{
		(Zeus::ConManager::Instance())->PrintStr2Console(temp + std::wstring(LASTPATCHN));
	}
	else
	{
		(Zeus::ConManager::Instance())->PrintStr2Console(temp + Zeus::PutNumber2Txt(Args.LastPatch_));
	}
	if(!(Args.fstString_.empty()))
	{
		temp = MC_STR1ARG;
		(Zeus::ConManager::Instance())->PrintStr2Console(temp + Args.fstString_);
	}
	if(!(Args.secString_.empty()))
	{
		temp = MC_STR2ARG;
		(Zeus::ConManager::Instance())->PrintStr2Console(temp + Args.secString_);
	}
	if((Args.ptg_.theta >= 0.0) && (Args.ptg_.phi >= 0.0))
	{
		temp = MC_ARGPOINTING;
		(Zeus::ConManager::Instance())->PrintStr2Console(temp + (Zeus::PutNumber2Txt(Args.ptg_.phi * RAD2DEGREE)) + std::wstring(L" long ; ") + Zeus::PutNumber2Txt((PIOVER2 - Args.ptg_.theta) * RAD2DEGREE) + std::wstring(L" lat"));
	}
}






