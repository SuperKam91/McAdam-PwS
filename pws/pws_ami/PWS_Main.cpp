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
 *  Copyright (C) 2005, 2009, Pedro Carvalho
 *  Author: Pedro Carvalho
 */


//----------------------------------
#include <time.h>
#include <iostream>
#include <memory>
#include "PWS_Globals.h"
#include "PWS_BackgroundProcessor.h"
#include "PWS_PatchProcessor.h"
#include "PWS_BackgroundProcessor.h"
#include "ZEUS_Priors.h"
#include "ZEUS_Debug.h"

//----------------------------------


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

#ifdef AMI
PatchProcessor*		RetVal(NULL);
#endif

#if	defined(WIN32) && !defined(HFIDMC) && !defined(LFIDPC) && !defined(AMI)
int _tmain(int argc, _TCHAR* argv[])
{
#elif AMI
int AMI_initializePwSlike(const char* ParFName)
{
	char bufferAMI[DOUBLETXTMAXSZ];
	int argc(2);
	char* argv[2] = { NULL, bufferAMI };
	strcpy(bufferAMI, ParFName);
#else
int main(int argc, char* argv[])
{
#endif //WIN32

	if (argc < 2)
	{
		PrintUsage();
		return -1;
	}

	//---------------------------------------------------
	double ZeroPlaneFrequency(1.0e32);
	double ZeroPlaneSZ_spectralSig(1.0e32);

	wchar_t				buffer[DOUBLETXTMAXSZ];
	CommandLineArgs		args;
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
		Zeus::ConManager::Initialize();
		fftw::init_fftw();
		PreProcessArgs(argc, argv, args);
		(Zeus::ConManager::Instance())->Initialise(argc, argv, args.ExecId_);
	}
#ifdef PWSMPI
	catch(Zeus::libMPIException& err)
	{
		ExcptString = (std::wstring(PWSMPIINITEXPTMSG) + err.what_Xmsg());
		PwSMPIerr = err.GetErr();
	}
#endif // PWSMPI
	catch (Zeus::libException& err)
	{
		ExcptString = (std::wstring(ZEUSPWSINITEXPTMSG) + err.what_Xmsg());
	}
	catch (...)
	{
		ExcptString = std::wstring(STDPWSINITEXPTMSG);
	}

	if (!(ExcptString.empty()))
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

	try{
		(PlanckInfo::Instance())->Initialise(args, time(NULL));
		const GlobalScalarVarsType&		GlbVars((PlanckInfo::Instance())->GetGlobalVars());

		PringArgs(args);

		ZeroPlaneFrequency = static_cast<double>(GlbVars.FreqsColl_[0].freq_);
		if ((GlbVars.SZ_ == 1) && (GlbVars.N_ObsPlanes_ == 1) && ((GlbVars.NonBlindDetection_ % 4) == 3) && (GlbVars.AssessmentKind_ == 0))
		{
			ZeroPlaneSZ_spectralSig = static_cast<double>(((PlanckInfo::Instance())->GetPlanckStaticInfoByFreq(GlbVars.FreqsColl_[0].freq_)).SZ_Signal_);
			(PlanckInfo::Instance())->SetSZGainTo1();
		}
	}
#ifdef PWSMPI
	catch(Zeus::libMPIException& err)
	{
		ExcptString = (std::wstring(PWSMPIINITEXPTMSG) + err.what_Xmsg());
		PwSMPIerr = err.GetErr();
	}
#endif //PWSMPI
	catch (Zeus::libException& err)
	{
		ExcptString = (std::wstring(ZEUSPWSEXCPTMSG) + err.what_Xmsg());
	}
	catch (std::exception& err)
	{
		ExcptString = (std::wstring(STANDARTEXCPTMSG) + Zeus::Achar2Wstr(err.what()));
	}
	catch (...)
	{
		ExcptString = std::wstring(UNDEFINEDEXCPTMSG);
	}

	if (!(ExcptString.empty()))
	{
		(Zeus::ConManager::Instance())->PrintErr2Console(ExcptString);
		(Zeus::ConManager::Instance())->CloseStreams();
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

	// AMI successfull exit
#ifdef AMI
	return 0;
#endif

	int							loop(true);

	do
	{
		try{
			if((args.FirstPatchNumber_ <= args.LastPatchNumber_) && (args.OutFName_.empty()))
			{
				loop = false;
				(PlanckInfo::Instance())->errParameterOutOfRange(std::wstring(GLOBAL_INTCATNAME));
			}

			for(;args.FirstPatchNumber_<=args.LastPatchNumber_;++(args.FirstPatchNumber_))
			{
				const Zeus::PatchGeomLineType&	GeomProps(((PlanckInfo::Instance())->GetGeoProps()).Storage_.at(args.FirstPatchNumber_));

				if(!(GeomProps.PatchValid_))
				{
					(Zeus::ConManager::Instance())->PrintStr2Console(std::wstring(STARTPATCHNUMBER) + Zeus::PutNumber2Txt(args.FirstPatchNumber_));
					std::auto_ptr<PatchProcessor> PtProc(new PatchProcessor(args.FirstPatchNumber_));
					PtProc->Initialise();
					PtProc->FindSources(DO_DETECTION);
					(Zeus::ConManager::Instance())->PrintStr2Console(std::wstring(REPORTINGPATCHNUMBER) + Zeus::PutNumber2Txt(args.FirstPatchNumber_));
					PtProc->ReportCurrentSrcs(GeomProps.SrcIndex_);
				}
				else
				{
					wchar_t	buffer[BUFFERMAXCHAR];

					PRINTINTOBUFFERFUNCT (buffer,BUFFERMAXCHAR,REPPATCHNUMBERNOTVALID,
						args.FirstPatchNumber_
					);
					(Zeus::ConManager::Instance())->PrintStr2Console(buffer);

					(PlanckInfo::Instance())->AppendNonValidObjs2File(GeomProps.SrcIndex_,args.FirstPatchNumber_);
				
				}
				(Zeus::ConManager::Instance())->Flush();
			}

			loop = false;

			(PlanckInfo::Instance())->FlushIntermedCat();

			if(!(args.FinalCatFName_.empty()))
			{
#ifdef PWSMPI
				int rk;
				MPI_Comm_rank(MPI_COMM_WORLD,&rk);
				if(!rk)
				{
#endif
				(PlanckInfo::Instance())->DoCompressCatalogue(args,ZeroPlaneFrequency,ZeroPlaneSZ_spectralSig);
#ifdef PWSMPI
				}
#endif
				(Zeus::ConManager::Instance())->Flush();
			}
			break;
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
			ExcptString = (std::wstring(MYERRORLOGLNTXT) + Zeus::PutNumber2Txt(args.FirstPatchNumber_) + 
							std::wstring(MYERRORSEPARATOR) + ExcptString);
			(Zeus::ConManager::Instance())->PrintErr2Console(ExcptString);
			ExcptString.clear();
			if(((PlanckInfo::Instance())->GetGlobalVars()).NonBlindDetection_ && loop)
			{
				(PlanckInfo::Instance())->AppendNonValidObjs2File(((PlanckInfo::Instance())->GetGeoProps()).Storage_.at(args.FirstPatchNumber_).SrcIndex_,args.FirstPatchNumber_);
			}
		}

		(Zeus::ConManager::Instance())->Flush();
		if(loop)
		{++(args.FirstPatchNumber_);}
	}while(loop);

	(Zeus::ConManager::Instance())->PrintStr2Console(FINISHINGPROGRAM);
	(Zeus::ConManager::Instance())->CloseStreams();
#ifdef WIN32
	 wprintf(L"\n\nPress any key to terminate the program\n\n");
 	_getch();
#endif //WIN32
#ifdef PWSMPI
	MPI_Finalize();
#endif //PWSMPI
#ifdef _DEBUG
//_CrtDumpMemoryLeaks();
#endif
	std::exit(0);
	// Never gets to here
	return 0;
}

void PringArgs(const CommandLineArgs&	args)
{
	std::wstring temp(PARAMSFILENAME);
	(Zeus::ConManager::Instance())->PrintStr2Console(temp + args.fname_);
	temp = OUTPUTFILENAME;
	(Zeus::ConManager::Instance())->PrintStr2Console(temp + args.OutFName_);
	temp = CONTEXTID;
	(Zeus::ConManager::Instance())->PrintStr2Console(temp + Zeus::PutNumber2Txt(args.ContextID_));
	temp = IOFLAG;
	(Zeus::ConManager::Instance())->PrintStr2Console(temp + Zeus::PutNumber2Txt(args.ExecId_));
	temp = FIRSTPATCH;
	if(args.FirstPatchNumber_ < 0)
	{
		(Zeus::ConManager::Instance())->PrintStr2Console(temp + Zeus::PutNumber2Txt(0));
	}
	else
	{
		(Zeus::ConManager::Instance())->PrintStr2Console(temp + Zeus::PutNumber2Txt(args.FirstPatchNumber_));
	}
	temp = LASTPATCH;
	if((args.LastPatchNumber_ < 0) )
	{
		(Zeus::ConManager::Instance())->PrintStr2Console(temp + std::wstring(LASTPATCHN));
	}
	else
	{
		(Zeus::ConManager::Instance())->PrintStr2Console(temp + Zeus::PutNumber2Txt(args.LastPatchNumber_));
	}

	if(!(args.FinalCatFName_.empty()))
	{
		temp = CATALOGUENAME;
		(Zeus::ConManager::Instance())->PrintStr2Console(temp + args.FinalCatFName_);
	}

	(Zeus::ConManager::Instance())->PrintStr2Console(L".\n");
}

void PrintUsage(void)
{
	Zeus::PrintErr2Console(USAGEMSG1);
	Zeus::PrintErr2Console(USAGEMSG2);
#ifdef WIN32
	Zeus::PrintErr2Console(PRESSANYKEYENDMSG);
	_getch();
#endif //WIN32
}
//
#if	defined(WIN32) && !defined(HFIDMC) && !defined(LFIDPC) && !defined(AMI)
void PreProcessArgs(int argc, _TCHAR* argv[],CommandLineArgs &args)
#else
void PreProcessArgs(int argc, char* argv[],CommandLineArgs &args)
#endif //WIN32
{
	std::wstring	temp;
	long			lNumber;

#ifdef	HFIDMC
	args.ContextID_			= 1;
#elif	LFIDPC
	args.ContextID_			= 2;
#elif	MyDEBUG
	args.ContextID_			= 0;
#else
	args.ContextID_			= 3;
#endif


	args.Argc_				= argc;
	if(argc >= 2)
	{
		args.fname_ = Zeus::FullTrim(Zeus::CharPtr2String(argv[1]));
	}
	else
	{
		args.fname_ = std::wstring(L"Oooops");
	}
	args.ExecId_			= 0x03;
	args.FirstPatchNumber_	= -1;
	args.LastPatchNumber_	= -1;

	ProcessEnvVariables(args);

	if(argc < 3) return;

	args.OutFName_			= Zeus::CharPtr2String(argv[2]);

	if(argc < 4) return;

	if(Zeus::get_number(temp = Zeus::CharPtr2String(argv[3]),lNumber))
		throw Zeus::libException(ERRCOD_PWS_NOTAINTEGER,ERRMSG_PWS_NOTAINTEGER,temp);

	args.ExecId_			= static_cast<int>(lNumber);

	if(argc < 5) return;

	if(Zeus::get_number(temp = Zeus::CharPtr2String(argv[4]),lNumber))
		throw Zeus::libException(ERRCOD_PWS_NOTAINTEGER,ERRMSG_PWS_NOTAINTEGER,temp);

	args.FirstPatchNumber_ = static_cast<int>(lNumber); 

	if(argc < 6) return;

	if(Zeus::get_number(temp = Zeus::CharPtr2String(argv[5]),lNumber))
		throw Zeus::libException(ERRCOD_PWS_NOTAINTEGER,ERRMSG_PWS_NOTAINTEGER,temp);

	args.LastPatchNumber_ = static_cast<int>(lNumber);

	for(int i = 6;i<argc;++i)
	{
		temp = Zeus::CharPtr2String(argv[i]);

		if(!(temp.empty()) &&  (temp.at(temp.size()-1) == OUTPUTFILELEADINGCHAR))
		{
			args.FinalCatFName_ = temp;
			continue;
		}

		args.CatalogueFNameColl_.push_back(Zeus::CharPtr2String(argv[i]));
	}
}

void	ProcessEnvVariables(CommandLineArgs &args)
{
	std::wstring	buf(Zeus::GetEnvVar(PWSLOGGER));

	if(buf.empty()) return;

	long	lNumber;

	if(Zeus::get_number(buf,lNumber))
		throw Zeus::libException(ERRCOD_PWS_NOTAINTEGER,ERRMSG_PWS_NOTAINTEGER,buf);

	args.ExecId_ = static_cast<int>(lNumber);
}

#ifdef AMI

//
// This main it is just for debugging
/*
int _tmain(int argc, _TCHAR* argv[])
{
	const int NParams(4);
	const int Patch(1);
	const int NoDuocimation(true);
	double In_params[6] = { 50.0, 50.0, 0.003, 10.0, 4.6, 0.81 };
//	double In_params[6] = {50.0, 50.0, 0.15, 39.0};
	double Out_params[6];
	double LikeValue;
	
	AMI_initializePwSlike(Zeus::Wstr2Str(std::wstring(argv[1])).c_str());
	AMI_fetchPatch(&Patch);
	for (int i = 0; i < 1;i++)
		AMI_PwSLikelihood(&LikeValue, &NParams, In_params, Out_params, &NoDuocimation);

	wchar_t	buffer[BUFFERMAXCHAR];

	PRINTINTOBUFFERFUNCT(buffer, BUFFERMAXCHAR, L"Likevalue -> %g, %g, %g, %g, %g, %g, %g\n", LikeValue, Out_params[0], Out_params[1], Out_params[2], Out_params[3], Out_params[4], Out_params[5]);
	(Zeus::ConManager::Instance())->PrintStr2Console(buffer);

	AMI_ReleasePatch();
}
*/
//
int	AMI_fetchPatch(const int *ptrPatchNumber)
{
	int PatchNumber(*ptrPatchNumber);
	std::wstring		ExcptString;

	if (RetVal) delete RetVal;
	RetVal = NULL;

	try{

		const Zeus::PatchGeomLineType&	GeomProps(((PlanckInfo::Instance())->GetGeoProps()).Storage_.at(PatchNumber));

		if (!(GeomProps.PatchValid_))
		{
			(Zeus::ConManager::Instance())->PrintStr2Console(std::wstring(STARTPATCHNUMBER) + Zeus::PutNumber2Txt(PatchNumber));
			RetVal = new PatchProcessor(PatchNumber);
//
			RetVal->Initialise();
			//
			RetVal->FindSources(NO_DETECTION);
			(Zeus::ConManager::Instance())->PrintStr2Console(std::wstring(REPORTINGPATCHNUMBER) + Zeus::PutNumber2Txt(PatchNumber));
		}
		else
		{
			wchar_t	buffer[BUFFERMAXCHAR];

			PRINTINTOBUFFERFUNCT(buffer, BUFFERMAXCHAR, REPPATCHNUMBERNOTVALID,PatchNumber);
			(Zeus::ConManager::Instance())->PrintStr2Console(buffer);
		}
	}
	catch (Zeus::libException& err)
	{
		ExcptString = (std::wstring(ZEUSPWSEXCPTMSG) + err.what_Xmsg());
	}
	catch (std::exception& err)
	{
		ExcptString = (std::wstring(STANDARTEXCPTMSG) + Zeus::Achar2Wstr(err.what()));
	}
	catch (...)
	{
		ExcptString = std::wstring(STANDARTEXCPTMSG);
	}

	if (!(ExcptString.empty()))
	{
		ExcptString = (std::wstring(MYERRORLOGLNTXT) + Zeus::PutNumber2Txt(PatchNumber) + std::wstring(MYERRORSEPARATOR) + ExcptString);
		(Zeus::ConManager::Instance())->PrintErr2Console(ExcptString);
		ExcptString.clear();
		if (RetVal) delete RetVal;
		RetVal = NULL;
	}

	(Zeus::ConManager::Instance())->Flush();
	return -(RetVal == NULL);
}
//
void	AMI_ReleasePatch(void)
{
	if (RetVal) delete RetVal;
	RetVal = NULL;
}
// all doubles
// 
// 0-> Xcoord
// 1-> Ycoord
// 2-> Y5r500
// 3-> Theta_s
// 4-> alpha
// 5-> beta
// 6-> c500
// 7-> gamma


int		AMI_PwSLikelihood(double* LikeValue, const int *ptrNParams, const double *In_params, double *Out_params, const int *ptrNoDuocimation)
{
	if (!RetVal) return -1;
	AmiLikeParams		ParamsIn;
	AmiLikeParams		ParamsOut;
	int NParams(*ptrNParams);

	if (NParams != 4 && NParams != 6)
		return -1;
	switch (NParams)
	{
	case 6 :
		ParamsIn.beta_ = In_params[4];
		ParamsIn.alpha_ = In_params[5];
	case 4 :
		ParamsIn.XCoord_ = In_params[0];
		ParamsIn.YCoord_ = In_params[1];
		ParamsIn.Ytot_ = In_params[2];
		ParamsIn.ThetaS_ = In_params[3];
		break;
	default:
		return -1;
	}
	ParamsIn.NParams_ = ParamsOut.NParams_ = NParams;
	if (RetVal->AmiPwSLikelihood(LikeValue, ParamsIn, ParamsOut, true))
		return -1;

	Out_params[4] = ParamsOut.beta_;
	Out_params[5] = ParamsOut.alpha_;
	Out_params[0] = ParamsOut.XCoord_;
	Out_params[1] = ParamsOut.YCoord_;
	Out_params[2] = ParamsOut.Ytot_;
	Out_params[3] = ParamsOut.ThetaS_;
	return 0;
}

#endif //AMI

