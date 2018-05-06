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


//----------------------------------------

#include <memory>
#include "ZEUS_ObjSZNgai.h"
#include "ZEUS_InOutPipeline.h"
#include "PWS_CatalogueComp.h"
#include "PWS_GlobalInfoStore.h"

#ifdef PWSMPI
#include <mpi.h>
#endif //PWSMPI

//----------------------------------------

void	PlankInfoStore::InitialiseStatic(void)
{
	PlanckStaticInfoAtomType temp;
	StaticInfoColl_.clear();
	temp.Freq_			= 30;
	temp.AntFWHM_		= 32.628;
	temp.OpNSide_		= 512;
#ifdef PLANCKNOMINALCHANNELS
	temp.SZ_Signal_		= -0.144;
#else
	temp.SZ_Signal_		= -0.144; // -0.12264
#endif
	StaticInfoColl_.push_back(temp);
//
	temp.Freq_			= 44;
	temp.AntFWHM_		= 27.958;
	temp.OpNSide_		= 512;
#ifdef PLANCKNOMINALCHANNELS
	temp.SZ_Signal_		= -0.293;
#else
	temp.SZ_Signal_		= -0.293; // -0.294075
#endif
	StaticInfoColl_.push_back(temp);
//
	temp.Freq_			= 70;
	temp.AntFWHM_		= 13.046;
	temp.OpNSide_		= 1024;
#ifdef PLANCKNOMINALCHANNELS
	temp.SZ_Signal_		= -0.633;
#else
	temp.SZ_Signal_		= -0.633; // -0.637099
#endif
	StaticInfoColl_.push_back(temp);
//
	temp.Freq_			= 100;
	temp.AntFWHM_		= 9.880;
	temp.OpNSide_		= 1024;
#ifdef PLANCKNOMINALCHANNELS
	temp.SZ_Signal_		= -0.970;
#else
	temp.SZ_Signal_		= -0.98227;
#endif
	StaticInfoColl_.push_back(temp);
//
	temp.Freq_			= 143;
	temp.AntFWHM_		= 7.180;
	temp.OpNSide_		= 2048;
#ifdef PLANCKNOMINALCHANNELS
	temp.SZ_Signal_		= -1.049;
#else
	temp.SZ_Signal_		= -1.032608;
#endif
	StaticInfoColl_.push_back(temp);
//
	temp.Freq_			= 217;
	temp.AntFWHM_		= 4.870;
	temp.OpNSide_		= 2048;
#ifdef PLANCKNOMINALCHANNELS
	temp.SZ_Signal_		= -0.008;
#else
	temp.SZ_Signal_		= 0.0991965;
#endif
	StaticInfoColl_.push_back(temp);
//
	temp.Freq_			= 353;
	temp.AntFWHM_		= 4.41;
	temp.OpNSide_		= 2048;
#ifdef PLANCKNOMINALCHANNELS
	temp.SZ_Signal_		= 1.725;
#else
	temp.SZ_Signal_		= 1.781322;
#endif
	StaticInfoColl_.push_back(temp);
//
	temp.Freq_			= 545;
	temp.AntFWHM_		= 4.47;
	temp.OpNSide_		= 2048;
#ifdef PLANCKNOMINALCHANNELS
	temp.SZ_Signal_		= 0.908;
#else
	temp.SZ_Signal_		= 0.838389;
#endif
	StaticInfoColl_.push_back(temp);
//
	temp.Freq_			= 857;
	temp.AntFWHM_		= 4.23;
	temp.OpNSide_		= 2048;
#ifdef PLANCKNOMINALCHANNELS
	temp.SZ_Signal_		= 0.064;
#else
	temp.SZ_Signal_		= 0.0589986;
#endif
	StaticInfoColl_.push_back(temp);
//
	temp.Freq_			= 0;
	temp.AntFWHM_		= 0.0;
	temp.OpNSide_		= 0;
	temp.SZ_Signal_		= 0.0;
	StaticInfoColl_.push_back(temp);
}

//
void	PlankInfoStore::ChangeStatic(int freq, double& FWHM)
{
	PlanckStaticInfoAtomType&	FreqInfo(GetPlanckStaticInfoByFreq(freq));
	if(FreqInfo.Freq_ == 0) return;
	if(FWHM < 0.0)
	{
		FWHM = FreqInfo.AntFWHM_;
		return;
	}
	FreqInfo.AntFWHM_ = FWHM;
}
//
void	PlankInfoStore::ProcessXtraArgs(CommandLineArgs& args)
{
	std::vector<std::wstring>	temp;
	std::vector<std::wstring>	temp1;
	std::wstring				temp2;

	Zeus::SplitSSString(GlobalScalarVars_.FinalCatalogueFileName_,std::wstring(L","),temp);

	std::vector<std::wstring>::const_iterator	piv(temp.begin());
	std::vector<std::wstring>::const_iterator	end(temp.end());

	for(;piv != end;++piv)
	{
		if(!(piv->empty()) && (piv->at(piv->size()-1) == OUTPUTFILELEADINGCHAR))
		{
			if(args.FinalCatFName_.empty())
			{
				temp2 = Zeus::FullTrim(*piv);
			}
			continue;
		}
		temp1.push_back(Zeus::FullTrim(*piv));
	}

	if(args.FinalCatFName_.empty() && !(temp2.empty()))
	{
		args.FinalCatFName_ = Zeus::FullTrim(temp2);
	}

	if(!(temp1.empty()))
	{
		args.CatalogueFNameColl_.insert(args.CatalogueFNameColl_.end(),temp1.begin(),temp1.end());
	}
}
//
void	PlankInfoStore::Initialise(CommandLineArgs& args,time_t InitialTime)
{
	InitialTime_ = InitialTime;
	InitialiseStatic();
	ContextID_ = args.ContextID_;
	(GlobalVars::Instance())->InitializeFromFile(args.fname_,args.ContextID_);

	ReadGlobalScalarVars();
#ifdef PWSMPI
	int	PwSMPIerr;
	if((PwSMPIerr=MPI_Comm_rank(MPI_COMM_WORLD,&GlobalScalarVars_.MPI_rank_))!=MPI_SUCCESS)
		throw Zeus::libMPIException(PwSMPIerr,L"MPI Error: cannot get thread rank",L"MPI_Comm_rank");
#endif //PWSMPI

	if((args.FirstPatchNumber_ < 0) && (GlobalScalarVars_.F_Patch_ >= 0))
	{args.FirstPatchNumber_ = GlobalScalarVars_.F_Patch_;}
	if((args.LastPatchNumber_ < 0) && (GlobalScalarVars_.L_Patch_ >= 0))
	{args.LastPatchNumber_ = GlobalScalarVars_.L_Patch_;}
	if(args.FirstPatchNumber_ < 0) args.FirstPatchNumber_ = 0;

	if(args.OutFName_.empty())
	{
		wchar_t			buffer[DOUBLETXTMAXSZ];
		swprintf(buffer,DOUBLETXTMAXSZ,L"%04d",GlobalScalarVars_.Sync_ID_);
		args.OutFName_ = IntermFName_ = (GlobalScalarVars_.IntermediateCatalogueFileName_.substr(0,GlobalScalarVars_.IntermediateCatalogueFileName_.size()-4) + std::wstring(buffer));
	}
	else
	{
		IntermFName_ = args.OutFName_;
	}

	if(!(GlobalScalarVars_.FinalCatalogueFileName_.empty()))
	{
		ProcessXtraArgs(args);
	}

	if(!(args.OutFName_.empty()))
	{args.CatalogueFNameColl_.push_back(args.OutFName_);}

	if(args.CatalogueFNameColl_.size() > 1)
	{
		std::sort(args.CatalogueFNameColl_.begin(),args.CatalogueFNameColl_.end());
		args.CatalogueFNameColl_.erase(std::unique(args.CatalogueFNameColl_.begin(),args.CatalogueFNameColl_.end()),
			args.CatalogueFNameColl_.end());
	}

	ReadGeoPropsHeader();

	if((args.LastPatchNumber_  > (GlobalScalarVars_.TotNPatches_-1)) || (args.LastPatchNumber_ < 0))
		args.LastPatchNumber_ = GlobalScalarVars_.TotNPatches_ - 1;

	ReadMapInfo();

	Zeus::PlaneBoundsType::SetDefaultUseBounds(GlobalScalarVars_.DefaultUseBounds_ ? Zeus::UB_USE : Zeus::UB_NOUSE);
	Zeus::PlaneBoundsType::SetDefaultBoundingFactor(static_cast<float>(GlobalScalarVars_.BoundTol_));

	if(GlobalScalarVars_.MaskFileName_.size() >= 2)
	{
		ReadHealpixMaskFile();
	}

	if(!(IntermFName_.empty()))
	{
		std::wstring	tDir;
		std::wstring	tFName(Zeus::ExtractFileName(IntermFName_,tDir));

		ObjWriter_ = Zeus::GetObjWriterHandler(ContextID_, GlobalScalarVars_.DirIn_ , tFName,TotalN_WrittenObjs_);
		
		if(ObjWriter_->Remove())
		{
			(Zeus::ConManager::Instance())->PrintStr2Console(std::wstring(L"Could not find/delete this object -> ") + ObjWriter_->GetCollID());
		}
		else
		{
			(Zeus::ConManager::Instance())->PrintStr2Console(std::wstring(L"Object successfully deleted -> ") + ObjWriter_->GetCollID());		
		}

		ObjWriter_->Initialize();

#if defined(HFIDMC) && defined(HFIDMC_EXTOBJECTS)
		if(GlobalScalarVars_.Data2Buffer_)
		{
#ifdef WIN32
			std::wstring	CatsExt(L"Cats\\");
#else
			std::wstring	CatsExt(L"Cats/");
#endif
			std::wstring	FullPath(GlobalScalarVars_.DirBuffer_ + CatsExt);
			Zeus::CreateDir(FullPath);

			std::wstring	tDir;
			std::wstring	tFName(Zeus::ExtractFileName(IntermFName_,tDir));

			ObjWriterExternals_ = Zeus::GetObjWriterHandler(1000,FullPath, tFName ,TotalN_WrittenObjsExt_);
			if(ObjWriterExternals_->Remove())
			{
				(Zeus::ConManager::Instance())->PrintStr2Console(std::wstring(L"Could not find/delete this object -> ") + ObjWriterExternals_->GetCollID());
			}
			else
			{
				(Zeus::ConManager::Instance())->PrintStr2Console(std::wstring(L"Object successfully deleted -> ") + ObjWriterExternals_->GetCollID());		
			}

			ObjWriterExternals_->Initialize();
		}
// if contours then creat contours group
		if(	(GlobalScalarVars_.SZ_ == 1) &&  (GlobalScalarVars_.N_ObsPlanes_ != 1) && (GlobalScalarVars_.Jf_Estimator_ == 0) &&
			(GlobalScalarVars_.AssessmentKind_ == 0)
		)
		{
			std::wstring	GrpName(GlobalScalarVars_.DirOut_ + std::wstring(L"_DegCont"));
			std::string		tGrpName(Zeus::Wstr2Str(GrpName));
			int				retries(MAXRETRIES);

#ifndef PWSMPI
_retry:
#endif
			if(PIOCheckGroup(const_cast<char *>(tGrpName.c_str())))
			{
				PIODOUBLE		PV[20];
				PIOErr			MyErr;
#ifdef PWSMPI
				int rk;
				MPI_Comm_rank(MPI_COMM_WORLD,&rk);
				if(!rk)
				{
_retry:
#endif

					for(int i=0;i<20;++i) PV[i]=0.0;
					if(MyErr = PIOCreateIMG2DGrp(const_cast<char *>(tGrpName.c_str()),"TAN","GALACTIC",GlobalScalarVars_.ProfParam_.ContBinsDef_,GlobalScalarVars_.ProfParam_.ContBinsDef_,GlobalScalarVars_.PixSz_ * RAD2DEGREE,GlobalScalarVars_.PixSz_ * RAD2DEGREE,(static_cast<double>(GlobalScalarVars_.ProfParam_.ContBinsDef_)/2.0),(static_cast<double>(GlobalScalarVars_.ProfParam_.ContBinsDef_)/2.0),PV)) 
					{
						if(--retries >= 0)
						{Zeus::MySleep(1);goto _retry;}
						throw Zeus::libException(ERROR_COD_HFIDMCERRCREATFL,std::wstring(ERROR_MSG_HFIDMCERRCREATFL) + std::wstring(L" -> ") + GrpName,L"PlankInfoStore::Initialise");
					}
#ifdef PWSMPI
				}
#endif
			}
#ifdef PWSMPI
			MPI_Barrier(MPI_COMM_WORLD);
#endif
		}
//
#endif
	}
}

void	PlankInfoStore::ReadGlobalScalarVars(void)
{
	std::vector<Zeus::VariableIdType>	ids;
	std::vector<Zeus::DBField>			values;

	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_DIRIN),
		Zeus::DBField()));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_DIRINMAPS),
		Zeus::DBField()));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_POINTINGS),
		Zeus::DBField()));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_DIRINMASKS),
		Zeus::DBField()));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_DIRINMASKSIN),
		Zeus::DBField(std::wstring(L""))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_DIROUTDATA),
		Zeus::DBField(std::wstring(L""))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_DIRBUFERDATA),
		Zeus::DBField(std::wstring(L""))));

	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_MASKFILENAME),
		Zeus::DBField(std::wstring(GLOBAL_REJECTMASKFNAME_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_DATASETNAME),
		Zeus::DBField(std::wstring(GLOBAL_DEFAULTDATASET))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_NFREQS),
		Zeus::DBField(static_cast<int>(GLOBAL_NFREQS_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_FIRSTPATCH),
		Zeus::DBField(static_cast<int>(-1))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_LASTPATCH),
		Zeus::DBField(static_cast<int>(-1))));
	ids.push_back(Zeus::VariableIdType(std::wstring(MC_SYNCID),
		Zeus::DBField(static_cast<int>(0))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_INTCATNAME),
		Zeus::DBField(std::wstring(L""))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_FINALCATNAME),
		Zeus::DBField(std::wstring(L""))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_NPRIORPLANES),
		Zeus::DBField(static_cast<int>(GLOBAL_NPRIORPLANES_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_NOTALIGNEDOBJS),
		Zeus::DBField(static_cast<int>(GLOBAL_NOTALIGNEDOBJS_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_PWSSEARCHPOS),
		Zeus::DBField(static_cast<int>(GLOBAL_PWSSEARCHPOS_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_USE2DFORMULA),
		Zeus::DBField(static_cast<int>(GLOBAL_USE2DFORMULA_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_USEBOUNDS),
		Zeus::DBField(static_cast<int>(GLOBAL_USEBOUNDS_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_ASSESS_TYPE),
		Zeus::DBField(static_cast<int>(GLOBAL_ASSESS_TYPE_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_BOUNDTOL),
		Zeus::DBField(static_cast<double>(GLOBAL_BOUNDTOL_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_SIGMA_THRESHOLD),
		Zeus::DBField(static_cast<double>(GLOBAL_SIGMA_THRESHOLD_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_SSUB_THRESHOLD),
		Zeus::DBField(static_cast<double>(GLOBAL_SSUBLIMITGAUSS_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_SZDETECTION),
		Zeus::DBField(GLOBAL_SZDETECTION_DEF)));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_CATALOGUESIGMA),
		Zeus::DBField(GLOBAL_CATALOGUESIGMA_DEF)));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_OUTPUTCOORDS),
		Zeus::DBField(static_cast<int>(GLOBAL_OUTPUTCOORDS_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_CATMERGE_AVGT),
		Zeus::DBField(static_cast<int>(GLOBAL_CATMERGE_AVGT_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_GALACTICSIGMA),
		Zeus::DBField(static_cast<double>(GLOBAL_GALACTICSIGMA_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_OUTUNITS),
		Zeus::DBField(static_cast<int>(GLOBAL_OUTUNITS_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_FILTERINGSCALES),
		Zeus::DBField(std::wstring(GLOBAL_FILTERINGSCALES_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_SZPROFILE),
		Zeus::DBField(static_cast<int>(GLOBAL_SZPROFILE_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_CACHESZ),
		Zeus::DBField(static_cast<double>(GLOBAL_CACHESZ_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_JFESTIMATTYPE),
		Zeus::DBField(static_cast<int>(GLOBAL_JFESTIMATTYPE_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_NONBLINDDETECTION),
		Zeus::DBField(static_cast<int>(GLOBAL_NONBLINDDETECTION_DEF))));

	(GlobalVars::Instance())->GetVar(ids,values);

	std::vector<Zeus::DBField>::iterator v_ptr(values.begin());

	std::wstring						DummyGcc1(*((*v_ptr++).GetPtr<std::wstring>()));
	std::wstring						DummyGcc2(*((*v_ptr++).GetPtr<std::wstring>()));
	std::wstring						DummyGcc3(*((*v_ptr++).GetPtr<std::wstring>()));
	std::wstring						DummyGcc4(*((*v_ptr++).GetPtr<std::wstring>()));
	std::wstring						DummyGcc5(*((*v_ptr++).GetPtr<std::wstring>()));
	std::wstring						DummyGcc6(*((*v_ptr++).GetPtr<std::wstring>()));
	std::wstring						DummyGcc7(*((*v_ptr++).GetPtr<std::wstring>()));

	(Zeus::ConManager::Instance())->PrintStr2Console(GLOBAL_PARAMETERS);

	GlobalScalarVars_.DirIn_			= Zeus::CorrectDir(Zeus::MakePath(DummyGcc1));
	Zeus::DumpScalarVariable(GLOBAL_DIRIN,GlobalScalarVars_.DirIn_);
	GlobalScalarVars_.DirInMaps_		= Zeus::CorrectDir(Zeus::MakePath(DummyGcc2));
	Zeus::DumpScalarVariable(GLOBAL_DIRINMAPS,GlobalScalarVars_.DirInMaps_);
	GlobalScalarVars_.DirPointings_		= Zeus::CorrectDir(Zeus::MakePath(DummyGcc3));
	Zeus::DumpScalarVariable(GLOBAL_POINTINGS,GlobalScalarVars_.DirPointings_);
	GlobalScalarVars_.DirInMasks_		= Zeus::CorrectDir(Zeus::MakePath(DummyGcc4));
	Zeus::DumpScalarVariable(GLOBAL_DIRINMASKS,GlobalScalarVars_.DirInMasks_);

	if(DummyGcc5.size()<2)
	{GlobalScalarVars_.DirInMasksIn_ = GlobalScalarVars_.DirInMasks_;}
	else
	{GlobalScalarVars_.DirInMasksIn_ = Zeus::CorrectDir(Zeus::MakePath(DummyGcc5));}
	Zeus::DumpScalarVariable(GLOBAL_DIRINMASKSIN,GlobalScalarVars_.DirInMasksIn_);

	if(DummyGcc6.size()<2)
	{GlobalScalarVars_.DirOut_ = GlobalScalarVars_.DirIn_;}
	else
	{GlobalScalarVars_.DirOut_ = Zeus::CorrectDir(Zeus::MakePath(DummyGcc6));}
	Zeus::DumpScalarVariable(GLOBAL_DIROUTDATA,GlobalScalarVars_.DirOut_);

	if(DummyGcc7.size()<2)
	{
		GlobalScalarVars_.DirBuffer_.clear();
		GlobalScalarVars_.Data2Buffer_	= 0;
	}
	else
	{
		GlobalScalarVars_.DirBuffer_ = Zeus::CorrectDir(Zeus::MakePath(DummyGcc7),1);
		GlobalScalarVars_.Data2Buffer_	= 1;	
	}
	Zeus::DumpScalarVariable(GLOBAL_DIRBUFERDATA,GlobalScalarVars_.DirBuffer_);

	GlobalScalarVars_.MaskFileName_		= *((*v_ptr++).GetPtr<std::wstring>());
	Zeus::DumpScalarVariable(GLOBAL_MASKFILENAME,GlobalScalarVars_.MaskFileName_);
	GlobalScalarVars_.OutputDataSetName_ = *((*v_ptr++).GetPtr<std::wstring>());
	Zeus::DumpScalarVariable(GLOBAL_DATASETNAME,GlobalScalarVars_.OutputDataSetName_);
	GlobalScalarVars_.N_ObsPlanes_		= (*v_ptr++).Get<int>();
	Zeus::DumpScalarVariable(GLOBAL_NFREQS,GlobalScalarVars_.N_ObsPlanes_);
	GlobalScalarVars_.F_Patch_			= (*v_ptr++).Get<int>();
	Zeus::DumpScalarVariable(GLOBAL_FIRSTPATCH,GlobalScalarVars_.F_Patch_);
	GlobalScalarVars_.L_Patch_			= (*v_ptr++).Get<int>();
	Zeus::DumpScalarVariable(GLOBAL_LASTPATCH,GlobalScalarVars_.L_Patch_);
	GlobalScalarVars_.Sync_ID_			= (*v_ptr++).Get<int>();
	Zeus::DumpScalarVariable(MC_SYNCID,GlobalScalarVars_.Sync_ID_);
	GlobalScalarVars_.IntermediateCatalogueFileName_	= *((*v_ptr++).GetPtr<std::wstring>());
	Zeus::DumpScalarVariable(GLOBAL_INTCATNAME,GlobalScalarVars_.IntermediateCatalogueFileName_);
	GlobalScalarVars_.FinalCatalogueFileName_			= *((*v_ptr++).GetPtr<std::wstring>());
	Zeus::DumpScalarVariable(GLOBAL_FINALCATNAME,GlobalScalarVars_.FinalCatalogueFileName_);
	GlobalScalarVars_.N_PriorPlanes_	= (*v_ptr++).Get<int>();
	GlobalScalarVars_.NotAlignedObjs_	= (*v_ptr++).Get<int>();
	Zeus::DumpScalarVariable(GLOBAL_NOTALIGNEDOBJS,GlobalScalarVars_.NotAlignedObjs_?std::wstring(GLOBAL_YES):std::wstring(GLOBAL_NO));
	GlobalScalarVars_.SearchPos_		= (*v_ptr++).Get<int>();
	Zeus::DumpScalarVariable(GLOBAL_PWSSEARCHPOS,GlobalScalarVars_.SearchPos_?std::wstring(GLOBAL_YES):std::wstring(GLOBAL_NO));
	GlobalScalarVars_.Use2DFormula_		= (*v_ptr++).Get<int>();
	Zeus::DumpScalarVariable(GLOBAL_USE2DFORMULA,GlobalScalarVars_.Use2DFormula_?std::wstring(GLOBAL_YES):std::wstring(GLOBAL_NO));
	GlobalScalarVars_.DefaultUseBounds_	= (*v_ptr++).Get<int>();
	Zeus::DumpScalarVariable(GLOBAL_USEBOUNDS,GlobalScalarVars_.DefaultUseBounds_?std::wstring(GLOBAL_YES):std::wstring(GLOBAL_NO));
	GlobalScalarVars_.AssessmentKind_	= (*v_ptr++).Get<int>();
	Zeus::DumpScalarVariable(GLOBAL_ASSESS_TYPE,GlobalScalarVars_.AssessmentKind_?std::wstring(GLOBAL_ASSTYPEBAYS):std::wstring(GLOBAL_ASSTYPEGLRT));
	GlobalScalarVars_.BoundTol_			= (*v_ptr++).Get<double>();
	Zeus::DumpScalarVariable(GLOBAL_BOUNDTOL,GlobalScalarVars_.BoundTol_);
	GlobalScalarVars_.NP_Sigma_			= (*v_ptr++).Get<double>();
	Zeus::DumpScalarVariable(GLOBAL_SIGMA_THRESHOLD,GlobalScalarVars_.NP_Sigma_);
	GlobalScalarVars_.SSubGaussLevel_			= (*v_ptr++).Get<double>();
	Zeus::DumpScalarVariable(GLOBAL_SSUB_THRESHOLD,GlobalScalarVars_.SSubGaussLevel_);
	GlobalScalarVars_.SZ_				= (*v_ptr++).Get<int>();
	Zeus::DumpScalarVariable(GLOBAL_SZDETECTION,GlobalScalarVars_.SZ_?std::wstring(GLOBAL_YES):std::wstring(GLOBAL_NO));
	GlobalScalarVars_.OutputSigma_		= (*v_ptr++).Get<double>();
	Zeus::DumpScalarVariable(GLOBAL_CATALOGUESIGMA,GlobalScalarVars_.OutputSigma_);
	GlobalScalarVars_.OutputCoords_		= (*v_ptr++).Get<int>();
	Zeus::DumpScalarVariable(GLOBAL_OUTPUTCOORDS,GlobalScalarVars_.OutputCoords_?(GlobalScalarVars_.OutputCoords_==2?std::wstring(L"G"):std::wstring(L"J")):std::wstring(L"E"));
	GlobalScalarVars_.OutputMergeAvGt_	 = (*v_ptr++).Get<int>();
	Zeus::DumpScalarVariable(GLOBAL_CATMERGE_AVGT,GlobalScalarVars_.OutputMergeAvGt_?std::wstring(GLOBAL_MERGTYPEHIGH):std::wstring(GLOBAL_MERGTYPEAVER));
	GlobalScalarVars_.OutputGalSigma_  = (*v_ptr++).Get<double>();
	Zeus::DumpScalarVariable(GLOBAL_GALACTICSIGMA,GlobalScalarVars_.OutputGalSigma_);
	GlobalScalarVars_.OutputUnits_		= (*v_ptr++).Get<int>();
	Zeus::DumpScalarVariable(GLOBAL_OUTUNITS,GlobalScalarVars_.OutputUnits_);
	std::wstring						ScaleBinsCollStr(*((*v_ptr++).GetPtr<std::wstring>()));
	Zeus::DumpScalarVariable(GLOBAL_FILTERINGSCALES,ScaleBinsCollStr);
	GlobalScalarVars_.ProfParam_.SZ_Profile_	= (*v_ptr++).Get<int>();
	Zeus::DumpScalarVariable(GLOBAL_SZPROFILE,(GlobalScalarVars_.ProfParam_.SZ_Profile_ == 2)?std::wstring(GLOBAL_SZPROFILENFW):((GlobalScalarVars_.ProfParam_.SZ_Profile_==1)?std::wstring(GLOBAL_SZPROFILEBETA):std::wstring(GLOBAL_SZPROFILEOTHER)));
	GlobalScalarVars_.CacheSz_			= (*v_ptr++).Get<double>();
	Zeus::DumpScalarVariable(GLOBAL_CACHESZ,GlobalScalarVars_.CacheSz_);
	GlobalScalarVars_.Jf_Estimator_		= (*v_ptr++).Get<int>();
	Zeus::DumpScalarVariable(GLOBAL_JFESTIMATTYPE,GlobalScalarVars_.Jf_Estimator_?std::wstring(GLOBAL_JFESTIMATORMEAN):std::wstring(GLOBAL_JFESTIMATORMODE));
	GlobalScalarVars_.NonBlindDetection_		= (*v_ptr++).Get<int>();
	Zeus::DumpScalarVariable(GLOBAL_NONBLINDDETECTION,GlobalScalarVars_.NonBlindDetection_?std::wstring(GLOBAL_YES):std::wstring(GLOBAL_NO));
	
	if(GlobalScalarVars_.N_PriorPlanes_)
	{
		GlobalScalarVars_.N_PriorPlanes_ = GlobalScalarVars_.N_ObsPlanes_;
	}
	Zeus::DumpScalarVariable(GLOBAL_NPRIORPLANES,GlobalScalarVars_.N_PriorPlanes_);

	if(!Zeus::GetNumbersHomogenColl(ScaleBinsCollStr,std::wstring(L","),GlobalScalarVars_.ScalesColl_))
	{errInvalidString(ScaleBinsCollStr);}

	GlobalScalarVars_.ContextID_		= ContextID_;

	ids.clear();
	values.clear();

	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_NSCALEBINS),
		Zeus::DBField(static_cast<int>(GlobalScalarVars_.SZ_?GLOBAL_NSCALEBINS_SZ_DEF:GLOBAL_NSCALEBINS_PS_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_SZPROFALPHA),
		Zeus::DBField(static_cast<double>(GLOBAL_SZPROFALPHA_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_SZPROFBETA),
		Zeus::DBField(static_cast<double>(GLOBAL_SZPROFBETA_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_SZPROFGAMMA),
		Zeus::DBField(static_cast<double>(GLOBAL_SZPROFGAMMA_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_SZPROF_C500),
		Zeus::DBField(static_cast<double>(GLOBAL_SZPROF_C500_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_SZVIRIALRATIO),
		Zeus::DBField(static_cast<double>(GLOBAL_SZVIRIALRATIO_DEF))));

	//----------------------------------------------------------------------

	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_OUTPUTRADIUSCAL),
		Zeus::DBField(static_cast<double>(GlobalScalarVars_.SZ_?GLOBAL_OUTPUTRADIUSCAL_SZ_DEF:GLOBAL_OUTPUTRADIUSCAL_PS_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_FLUXCALIBCTE),
		Zeus::DBField(static_cast<double>(GLOBAL_FLUXCALIBCTE_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_SRCMAXSCALE),
		Zeus::DBField(static_cast<double>(GlobalScalarVars_.SZ_?GLOBAL_SRCMAXSCALE_SZ_DEF:(GlobalScalarVars_.NonBlindDetection_?GLOBAL_SRCMAXSCALE_PSNB_DEF:GLOBAL_SRCMAXSCALE_PS_DEF)))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_TWOSTEPSDETECT),
		Zeus::DBField(static_cast<int>(GLOBAL_TWOSTEPSDETECTNOPRIOR_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_MELTMAXDIST),
		Zeus::DBField(static_cast<double>(GlobalScalarVars_.SZ_?GLOBAL_MELTMAXDISTSZ_DEF:GLOBAL_MELTMAXDISTPS_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_MELTMAXDISTSZLIMIT),
		Zeus::DBField(static_cast<double>(GlobalScalarVars_.SZ_?GLOBAL_MELTMAXDISTSZMAXLIMIT_DEF:GLOBAL_MELTMAXDISTPSMAXLIMIT_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_FLUXTHRESHOLD),
		Zeus::DBField(static_cast<double>(GlobalScalarVars_.SZ_?GLOBAL_FLUXTHRESHOLD_SZ_DEF:GLOBAL_FLUXTHRESHOLD_PS_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_OUTPUTLAT),
		Zeus::DBField(static_cast<int>(GLOBAL_HARDCONSTLEVEL_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_OUTPUTDEGREES),
		Zeus::DBField(static_cast<int>(GlobalScalarVars_.SZ_?GLOBAL_OUTPUTDEGREES_SZ_DEF:GLOBAL_OUTPUTDEGREES_PS_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_GALACTICCUT),
		Zeus::DBField(static_cast<double>(GlobalScalarVars_.SZ_?GLOBAL_GALACTICCUT_SZ_DEF:GLOBAL_GALACTICCUT_PS_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_PRIORFLUXMIN),
		Zeus::DBField(static_cast<double>(GlobalScalarVars_.SZ_?GLOBAL_PRIORFLUXMIN_SZ_DEF:GLOBAL_PRIORFLUXMIN_PS_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_PRIORFLUXMAX),
		Zeus::DBField(static_cast<double>(GlobalScalarVars_.SZ_?GLOBAL_PRIORFLUXMAX_SZ_DEF:GLOBAL_PRIORFLUXMAX_PS_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_PRIORRADDISTMIN),
		Zeus::DBField(static_cast<double>(GLOBAL_PRIORRADDISTMIN_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_PRIORRADDISTMAX),
		Zeus::DBField(static_cast<double>(GLOBAL_PRIORRADDISTMAX_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_PRIORMASSMIN),
		Zeus::DBField(static_cast<double>(GLOBAL_PRIORMASSMIN_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_PRIORMASSMAX),
		Zeus::DBField(static_cast<double>(GLOBAL_PRIORMASSMAX_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_PRIORMAXZ),
		Zeus::DBField(static_cast<double>(GLOBAL_PRIORMAXZ_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_PRIORTEMPMIN),
		Zeus::DBField(static_cast<double>(GLOBAL_PRIORTEMPMIN_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_PRIORTEMPMAX),
		Zeus::DBField(static_cast<double>(GLOBAL_PRIORTEMPMAX_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_PRIORGASMASSRATIO),
		Zeus::DBField(static_cast<double>(GLOBAL_PRIORGASMASSRATIO_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_PRIORSIGMA8),
		Zeus::DBField(static_cast<double>(GLOBAL_PRIORSIGMA8_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_USEPRESSSCHT),
		Zeus::DBField(static_cast<int>(GLOBAL_USEPRESSSCHT_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_MN_NLIVEPOINTS),
		Zeus::DBField(static_cast<int>(GLOBAL_MN_NLIVEPOINTS_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_MN_NINDIVSAMPLES),
		Zeus::DBField(static_cast<int>(GLOBAL_MN_NINDIVSAMPLES_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_MN_FRACTOLEV),
		Zeus::DBField(static_cast<double>(GLOBAL_MN_FRACTOLEV_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_MN_XTRAENLFACT),
		Zeus::DBField(static_cast<double>(GLOBAL_MN_XTRAENLFACT_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_PRIORFLUXEXP),
		Zeus::DBField(static_cast<double>(GlobalScalarVars_.SZ_?GLOBAL_PRIORFLUXEXP_SZ_DEF:GLOBAL_PRIORFLUXEXP_PS_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_PRIORSRCSCALEMIN),
		Zeus::DBField(static_cast<double>(GlobalScalarVars_.SZ_?GLOBAL_PRIORSRCSCALEMIN_SZ_DEF:GLOBAL_PRIORSRCSCALEMIN_PS_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_PRIORSRCSCALEMAX),
		Zeus::DBField(static_cast<double>(GlobalScalarVars_.SZ_?GLOBAL_PRIORSRCSCALEMAX_SZ_DEF:(GlobalScalarVars_.NonBlindDetection_?GLOBAL_SRCMAXSCALE_PSNB_DEF:GLOBAL_SRCMAXSCALE_PS_DEF)))));//
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_PRIORSRCSCALEEXP),
		Zeus::DBField(static_cast<double>(GlobalScalarVars_.SZ_?GLOBAL_PRIORSRCSCALEEXP_SZ_DEF:GLOBAL_PRIORSRCSCALEEXP_PS_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_PRIORAVSRCPATCH),
		Zeus::DBField(static_cast<double>(GlobalScalarVars_.SZ_?GLOBAL_PRIORAVSRCSZ_DEF:GLOBAL_PRIORAVSRCPS_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_OUTPURITY),
		Zeus::DBField(static_cast<double>(GLOBAL_OUTPURITY_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_SZPROFVARALPHAMIN),
		Zeus::DBField(static_cast<double>(GLOBAL_SZPROFVARALPHAMIN_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_SZPROFVARALPHAMAX),
		Zeus::DBField(static_cast<double>(GLOBAL_SZPROFVARALPHAMAX_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_SZPROFVARALPHABIN),
		Zeus::DBField(static_cast<int>(GLOBAL_SZPROFVARALPHABIN_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_SZPROFVARBETAMIN),
		Zeus::DBField(static_cast<double>(GLOBAL_SZPROFVARBETAMIN_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_SZPROFVARBETAMAX),
		Zeus::DBField(static_cast<double>(GLOBAL_SZPROFVARBETAMAX_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_SZPROFVARBETABIN),
		Zeus::DBField(static_cast<int>(GLOBAL_SZPROFVARBETABIN_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_SZPROFVARC500MIN),
		Zeus::DBField(static_cast<double>(GLOBAL_SZPROFVARC500MIN_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_SZPROFVARC500MAX),
		Zeus::DBField(static_cast<double>(GLOBAL_SZPROFVARC500MAX_DEF))));
	ids.push_back(Zeus::VariableIdType(std::wstring(GLOBAL_SZPROFVARC500BIN),
		Zeus::DBField(static_cast<int>(GLOBAL_SZPROFVARC500BIN_DEF))));
//
	(GlobalVars::Instance())->GetVar(ids,values);

	v_ptr = values.begin();

	GlobalScalarVars_.N_ScaleBins_		= (*v_ptr++).Get<int>();
	Zeus::DumpScalarVariable(GLOBAL_NSCALEBINS,GlobalScalarVars_.N_ScaleBins_);

	GlobalScalarVars_.ProfParam_.MNFW_alpha_	= (*v_ptr++).Get<double>();
	Zeus::DumpScalarVariable(GLOBAL_SZPROFALPHA,GlobalScalarVars_.ProfParam_.MNFW_alpha_);
	GlobalScalarVars_.ProfParam_.MNFW_beta_	= (*v_ptr++).Get<double>();
	Zeus::DumpScalarVariable(GLOBAL_SZPROFBETA,GlobalScalarVars_.ProfParam_.MNFW_beta_);
	GlobalScalarVars_.ProfParam_.MNFW_gamma_	= (*v_ptr++).Get<double>();
	Zeus::DumpScalarVariable(GLOBAL_SZPROFGAMMA,GlobalScalarVars_.ProfParam_.MNFW_gamma_);
	GlobalScalarVars_.ProfParam_.MNFW_C500_	= (*v_ptr++).Get<double>();
	Zeus::DumpScalarVariable(GLOBAL_SZPROF_C500,GlobalScalarVars_.ProfParam_.MNFW_C500_);

	GlobalScalarVars_.ProfParam_.R500_Ratio_ = ((*v_ptr++).Get<double>());
	GlobalScalarVars_.ProfParam_.VirialRatio_ = GlobalScalarVars_.ProfParam_.R500_Ratio_ * GlobalScalarVars_.ProfParam_.MNFW_C500_;

	Zeus::DumpScalarVariable(GLOBAL_SZVIRIALRATIO,GlobalScalarVars_.ProfParam_.VirialRatio_);
	GlobalScalarVars_.ProfParam_.RadiusCalCte_ = (*v_ptr++).Get<double>();
	Zeus::DumpScalarVariable(GLOBAL_OUTPUTRADIUSCAL,GlobalScalarVars_.ProfParam_.RadiusCalCte_);
	GlobalScalarVars_.ProfParam_.FluxCalibCte_	= (*v_ptr++).Get<double>();
	Zeus::DumpScalarVariable(GLOBAL_FLUXCALIBCTE,GlobalScalarVars_.ProfParam_.FluxCalibCte_);
	GlobalScalarVars_.SrcMaxScale_				= (*v_ptr++).Get<double>();
	Zeus::DumpScalarVariable(GLOBAL_SRCMAXSCALE,GlobalScalarVars_.SrcMaxScale_);
	GlobalScalarVars_.TwoStepsDetection_		= (*v_ptr++).Get<int>();
	Zeus::DumpScalarVariable(GLOBAL_TWOSTEPSDETECT,GlobalScalarVars_.TwoStepsDetection_);
	GlobalScalarVars_.OutputMeltMaxDist_		= (*v_ptr++).Get<double>();
	Zeus::DumpScalarVariable(GLOBAL_MELTMAXDIST,GlobalScalarVars_.OutputMeltMaxDist_);
	GlobalScalarVars_.OutputMeltSZMaxLimit_		= (*v_ptr++).Get<double>();
	Zeus::DumpScalarVariable(GLOBAL_MELTMAXDISTSZLIMIT,GlobalScalarVars_.OutputMeltSZMaxLimit_);
	GlobalScalarVars_.OutputFluxThreshold_		= (*v_ptr++).Get<double>();
	Zeus::DumpScalarVariable(GLOBAL_FLUXTHRESHOLD,GlobalScalarVars_.OutputFluxThreshold_);
	GlobalScalarVars_.OutputLat_				= (*v_ptr++).Get<int>();
	Zeus::DumpScalarVariable(GLOBAL_OUTPUTLAT,GlobalScalarVars_.OutputLat_?std::wstring(GLOBAL_YES):std::wstring(GLOBAL_NO));
	GlobalScalarVars_.OutputDegrees_			= (*v_ptr++).Get<int>();
	Zeus::DumpScalarVariable(GLOBAL_OUTPUTDEGREES,GlobalScalarVars_.OutputDegrees_?std::wstring(GLOBAL_YES):std::wstring(GLOBAL_NO));
	GlobalScalarVars_.OutputGalCut_				= (*v_ptr++).Get<double>();	
	Zeus::DumpScalarVariable(GLOBAL_GALACTICCUT,GlobalScalarVars_.OutputGalCut_);
	GlobalScalarVars_.PriorFluxMin_				= (*v_ptr++).Get<double>();
	Zeus::DumpScalarVariable(GLOBAL_PRIORFLUXMIN,GlobalScalarVars_.PriorFluxMin_);
	GlobalScalarVars_.PriorFluxMax_				= (*v_ptr++).Get<double>();
	Zeus::DumpScalarVariable(GLOBAL_PRIORFLUXMAX,GlobalScalarVars_.PriorFluxMax_);
	GlobalScalarVars_.PriorRadDistMin_			= (*v_ptr++).Get<double>();
	Zeus::DumpScalarVariable(GLOBAL_PRIORRADDISTMIN,GlobalScalarVars_.PriorRadDistMin_);
	GlobalScalarVars_.PriorRadDistMax_			= (*v_ptr++).Get<double>();
	Zeus::DumpScalarVariable(GLOBAL_PRIORRADDISTMAX,GlobalScalarVars_.PriorRadDistMax_);
	GlobalScalarVars_.PriorMassMin_				= (*v_ptr++).Get<double>();
	Zeus::DumpScalarVariable(GLOBAL_PRIORMASSMIN,GlobalScalarVars_.PriorMassMin_);
	GlobalScalarVars_.PriorMassMax_				= (*v_ptr++).Get<double>();
	Zeus::DumpScalarVariable(GLOBAL_PRIORMASSMAX,GlobalScalarVars_.PriorMassMax_);
	GlobalScalarVars_.PriorMaxZ_				= (*v_ptr++).Get<double>();
	Zeus::DumpScalarVariable(GLOBAL_PRIORMAXZ,GlobalScalarVars_.PriorMaxZ_);
	GlobalScalarVars_.PriorTempMin_				= (*v_ptr++).Get<double>();
	Zeus::DumpScalarVariable(GLOBAL_PRIORTEMPMIN,GlobalScalarVars_.PriorTempMin_);
	GlobalScalarVars_.PriorTempMax_				= (*v_ptr++).Get<double>();
	Zeus::DumpScalarVariable(GLOBAL_PRIORTEMPMAX,GlobalScalarVars_.PriorTempMax_);
	GlobalScalarVars_.PriorGassMassRatio_		= (*v_ptr++).Get<double>();
	Zeus::DumpScalarVariable(GLOBAL_PRIORGASMASSRATIO,GlobalScalarVars_.PriorGassMassRatio_);
	GlobalScalarVars_.Sigma8_					= (*v_ptr++).Get<double>();
	Zeus::DumpScalarVariable(GLOBAL_PRIORSIGMA8,GlobalScalarVars_.Sigma8_);
	GlobalScalarVars_.PriorUsePressSchether_	= (*v_ptr++).Get<int>();
	Zeus::DumpScalarVariable(GLOBAL_USEPRESSSCHT,(GlobalScalarVars_.PriorUsePressSchether_==1)?std::wstring(GLOBAL_MASSFUNCTPS):(GlobalScalarVars_.PriorUsePressSchether_==2)?std::wstring(GLOBAL_MASSFUNCTJK):std::wstring(GLOBAL_SZPROFILEOTHER));
	GlobalScalarVars_.MuNe_NLivePoints_			= (*v_ptr++).Get<int>();
	Zeus::DumpScalarVariable(GLOBAL_MN_NLIVEPOINTS,GlobalScalarVars_.MuNe_NLivePoints_);
	GlobalScalarVars_.MuNe_NIndivSamples_		= (*v_ptr++).Get<int>(); 
	Zeus::DumpScalarVariable(GLOBAL_MN_NINDIVSAMPLES,GlobalScalarVars_.MuNe_NIndivSamples_);
	GlobalScalarVars_.MuNe_FractTolEv_			= (*v_ptr++).Get<double>();
	Zeus::DumpScalarVariable(GLOBAL_MN_FRACTOLEV,GlobalScalarVars_.MuNe_FractTolEv_);
	GlobalScalarVars_.MuNe_xtraEnlgFactor_		= (*v_ptr++).Get<double>();
	Zeus::DumpScalarVariable(GLOBAL_MN_XTRAENLFACT,GlobalScalarVars_.MuNe_xtraEnlgFactor_);
	GlobalScalarVars_.PriorFluxExp_				= (*v_ptr++).Get<double>();
	Zeus::DumpScalarVariable(GLOBAL_PRIORFLUXEXP,GlobalScalarVars_.PriorFluxExp_);
	GlobalScalarVars_.PriorSrcMinScale_			= (*v_ptr++).Get<double>();
	Zeus::DumpScalarVariable(GLOBAL_PRIORSRCSCALEMIN,GlobalScalarVars_.PriorSrcMinScale_);
	GlobalScalarVars_.PriorSrcMaxScale_			= (*v_ptr++).Get<double>();
	Zeus::DumpScalarVariable(GLOBAL_PRIORSRCSCALEMAX,GlobalScalarVars_.PriorSrcMaxScale_);
	GlobalScalarVars_.PriorSrcScaleExp_			= (*v_ptr++).Get<double>();
	GlobalScalarVars_.PriorAvObjectsPatch_		= (*v_ptr++).Get<double>();
	Zeus::DumpScalarVariable(GLOBAL_PRIORAVSRCPATCH,GlobalScalarVars_.PriorAvObjectsPatch_);
	GlobalScalarVars_.OutputPurity_				= (*v_ptr++).Get<double>(); 
	Zeus::DumpScalarVariable(GLOBAL_OUTPURITY,GlobalScalarVars_.OutputPurity_);
	GlobalScalarVars_.OutputPurity_ *= 0.6;
	GlobalScalarVars_.ProfParamVar_.alpha_.min_	= (*v_ptr++).Get<double>();
	Zeus::DumpScalarVariable(GLOBAL_SZPROFVARALPHAMIN,GlobalScalarVars_.ProfParamVar_.alpha_.min_);
	GlobalScalarVars_.ProfParamVar_.alpha_.max_	= (*v_ptr++).Get<double>();
	Zeus::DumpScalarVariable(GLOBAL_SZPROFVARALPHAMAX,GlobalScalarVars_.ProfParamVar_.alpha_.max_);
	GlobalScalarVars_.ProfParamVar_.alpha_.nbins_	= (*v_ptr++).Get<int>();
	if(GlobalScalarVars_.ProfParamVar_.alpha_.nbins_>0)
		GlobalScalarVars_.ProfParamVar_.alpha_.nbins_ = Zeus::getLessPow2p1(GlobalScalarVars_.ProfParamVar_.alpha_.nbins_);
	Zeus::DumpScalarVariable(GLOBAL_SZPROFVARALPHABIN,GlobalScalarVars_.ProfParamVar_.alpha_.nbins_);
	GlobalScalarVars_.ProfParamVar_.beta_.min_	= (*v_ptr++).Get<double>();
	Zeus::DumpScalarVariable(GLOBAL_SZPROFVARBETAMIN,GlobalScalarVars_.ProfParamVar_.beta_.min_);
	GlobalScalarVars_.ProfParamVar_.beta_.max_	= (*v_ptr++).Get<double>();
	Zeus::DumpScalarVariable(GLOBAL_SZPROFVARBETAMAX,GlobalScalarVars_.ProfParamVar_.beta_.max_);
	GlobalScalarVars_.ProfParamVar_.beta_.nbins_	= (*v_ptr++).Get<int>();
	if(GlobalScalarVars_.ProfParamVar_.beta_.nbins_>0)
		GlobalScalarVars_.ProfParamVar_.beta_.nbins_ = Zeus::getLessPow2p1(GlobalScalarVars_.ProfParamVar_.beta_.nbins_);
	Zeus::DumpScalarVariable(GLOBAL_SZPROFVARBETABIN,GlobalScalarVars_.ProfParamVar_.beta_.nbins_);

	GlobalScalarVars_.ProfParamVar_.C500_.min_ = (*v_ptr++).Get<double>();
	Zeus::DumpScalarVariable(GLOBAL_SZPROFVARC500MIN, GlobalScalarVars_.ProfParamVar_.C500_.min_);
	GlobalScalarVars_.ProfParamVar_.C500_.max_ = (*v_ptr++).Get<double>();
	Zeus::DumpScalarVariable(GLOBAL_SZPROFVARC500MAX, GlobalScalarVars_.ProfParamVar_.C500_.max_);
	GlobalScalarVars_.ProfParamVar_.C500_.nbins_ = (*v_ptr++).Get<int>();
	if (GlobalScalarVars_.ProfParamVar_.C500_.nbins_>0)
		GlobalScalarVars_.ProfParamVar_.C500_.nbins_ = Zeus::getLessPow2p1(GlobalScalarVars_.ProfParamVar_.C500_.nbins_);
	Zeus::DumpScalarVariable(GLOBAL_SZPROFVARC500BIN, GlobalScalarVars_.ProfParamVar_.C500_.nbins_);

	if(GlobalScalarVars_.PriorMassMax_ < 0.0)
	{
		if(GlobalScalarVars_.ProfParamVar_.IsInit())
		{GlobalScalarVars_.ProfParam_.ContBinsDef_ = GlobalScalarVars_.ProfParamVar_.beta_.nbins_;}
	}

	if(GlobalScalarVars_.TwoStepsDetection_ < 1)
	{GlobalScalarVars_.TwoStepsDetection_ = 1;}
	if(GlobalScalarVars_.PriorSrcMaxScale_ < 0.0)
	{GlobalScalarVars_.PriorSrcMaxScale_ = GlobalScalarVars_.SrcMaxScale_;}
	if(!(GlobalScalarVars_.SZ_))
	{GlobalScalarVars_.PriorSrcScaleExp_ = 1.0/(2.0 * GlobalScalarVars_.PriorSrcMaxScale_);}
	Zeus::DumpScalarVariable(GLOBAL_PRIORSRCSCALEEXP,GlobalScalarVars_.PriorSrcScaleExp_);

	ProcessScalesColl();
	//GlbVars_.ProfParam_.VirialRatio_
	{
		Zeus::GammaFuncts f;
		const double Y_sphR500(f.betainc((3.0-GlobalScalarVars_.ProfParam_.MNFW_gamma_)/GlobalScalarVars_.ProfParam_.MNFW_alpha_,
			(GlobalScalarVars_.ProfParam_.MNFW_beta_-3.0)/GlobalScalarVars_.ProfParam_.MNFW_alpha_,
			std::pow(GlobalScalarVars_.ProfParam_.MNFW_C500_,GlobalScalarVars_.ProfParam_.MNFW_alpha_)/(1.0+std::pow(GlobalScalarVars_.ProfParam_.MNFW_C500_,GlobalScalarVars_.ProfParam_.MNFW_alpha_))));

		const double Y_sph5R500(f.betainc((3.0-GlobalScalarVars_.ProfParam_.MNFW_gamma_)/GlobalScalarVars_.ProfParam_.MNFW_alpha_,
			(GlobalScalarVars_.ProfParam_.MNFW_beta_-3.0)/GlobalScalarVars_.ProfParam_.MNFW_alpha_,
			std::pow(GlobalScalarVars_.ProfParam_.VirialRatio_,GlobalScalarVars_.ProfParam_.MNFW_alpha_)/(1.0+std::pow(GlobalScalarVars_.ProfParam_.VirialRatio_,GlobalScalarVars_.ProfParam_.MNFW_alpha_))));

		const double tBetaCte((GlobalScalarVars_.ProfParam_.MNFW_beta_<=5.05)?5.05:GlobalScalarVars_.ProfParam_.MNFW_beta_);

		double M2(f.betainc((5.0-GlobalScalarVars_.ProfParam_.MNFW_gamma_)/GlobalScalarVars_.ProfParam_.MNFW_alpha_,
			(tBetaCte-5.0)/GlobalScalarVars_.ProfParam_.MNFW_alpha_,
			std::pow(GlobalScalarVars_.ProfParam_.VirialRatio_,GlobalScalarVars_.ProfParam_.MNFW_alpha_)/(1.0+std::pow(GlobalScalarVars_.ProfParam_.VirialRatio_,GlobalScalarVars_.ProfParam_.MNFW_alpha_))));

		M2	*= (f.beta((5.0-GlobalScalarVars_.ProfParam_.MNFW_gamma_)/GlobalScalarVars_.ProfParam_.MNFW_alpha_,
			(tBetaCte-5.0)/GlobalScalarVars_.ProfParam_.MNFW_alpha_));

		double Norm(f.betainc((3.0-GlobalScalarVars_.ProfParam_.MNFW_gamma_)/GlobalScalarVars_.ProfParam_.MNFW_alpha_,
			(tBetaCte-3.0)/GlobalScalarVars_.ProfParam_.MNFW_alpha_,
			std::pow(GlobalScalarVars_.ProfParam_.VirialRatio_,GlobalScalarVars_.ProfParam_.MNFW_alpha_)/(1.0+std::pow(GlobalScalarVars_.ProfParam_.VirialRatio_,GlobalScalarVars_.ProfParam_.MNFW_alpha_))));

		Norm *= (f.beta((3.0-GlobalScalarVars_.ProfParam_.MNFW_gamma_)/GlobalScalarVars_.ProfParam_.MNFW_alpha_,
			(tBetaCte-3.0)/GlobalScalarVars_.ProfParam_.MNFW_alpha_));

		GlobalScalarVars_.ProfParam_.MNFW_Ratio_CY500CYR500_ = Y_sph5R500/Y_sphR500;
		Zeus::DumpScalarVariable(GLOBAL_SZPROFCY500CYR500,GlobalScalarVars_.ProfParam_.MNFW_Ratio_CY500CYR500_);
		GlobalScalarVars_.ProfParam_.MNFW_2ndMoment1D_ = M2/(3.0*Norm); // * 2/(3 * 2)
		Zeus::DumpScalarVariable(GLOBAL_SZPROF2NDMOMENT,GlobalScalarVars_.ProfParam_.MNFW_2ndMoment1D_);
	}
	// Yes, subtract sources
	GlobalScalarVars_.SubSources_ = 1;

	(Zeus::ConManager::Instance())->PrintStr2Console(L".\n");
}

void	PlankInfoStore::ReadGeoPropsHeader(void)
{
	std::wstring			tDir(GlobalScalarVars_.DirPointings_);

//	if(GlobalScalarVars_.Data2Buffer_)
	if(false)
	{
#ifdef WIN32
		std::wstring	GeomExt(L"Geom\\");
#else
		std::wstring	GeomExt(L"Geom/");
#endif
		if(!Zeus::CheckDir(tDir = (GlobalScalarVars_.DirBuffer_  + GeomExt)))
		{
			tDir = Zeus::RemoveInstanceName(GlobalScalarVars_.DirBuffer_) + GeomExt;
		}
	}

/*
GlobalScalarVars_.Data2Buffer_?1000:ContextID_
*/
	std::auto_ptr<Zeus::GenCollReader<Zeus::PatchGeomType> >
		FReader(Zeus::GetPatchGeomFileReaderHandler(Loki::Type2Type<Zeus::PatchGeomType>(),
		ContextID_,
		tDir,std::wstring(PATCH_GEOMPROPSFNAME)));

	FReader->Initialize();
	FReader->Read();
	FReader->Release(PatchGeomInfo_);

	GlobalScalarVars_.PatchBorder_	= (PatchGeomInfo_.Header_.PtchBorder_ >> 1);
	GlobalScalarVars_.PatchSz_		= PatchGeomInfo_.Header_.PtchSz_;
	GlobalScalarVars_.NSize_		= PatchGeomInfo_.Header_.NSide_;
	GlobalScalarVars_.TotNPatches_	= PatchGeomInfo_.Header_.NTotalPatches_;

	const double pixArea(PITIMES4 / static_cast<double>(12 * PatchGeomInfo_.Header_.NSide_ * PatchGeomInfo_.Header_.NSide_));
	GlobalScalarVars_.PixSz_		= std::sqrt(pixArea);
	double t(static_cast<double>(PatchGeomInfo_.Header_.PtchSz_ - (GlobalScalarVars_.PatchBorder_<<1)));
	GlobalScalarVars_.PriorAvObjectsPatch_ *= (t*t*pixArea);

	Zeus::DumpScalarVariable(MAKEWCHARARR(GEOHEADERID_NSIZE),GlobalScalarVars_.NSize_);
	Zeus::DumpScalarVariable(MAKEWCHARARR(GEOHEADERID_PATCHSZ),GlobalScalarVars_.PatchSz_);
	Zeus::DumpScalarVariable(MAKEWCHARARR(GEOHEADERID_PATCHBORDER),GlobalScalarVars_.PatchBorder_);
	Zeus::DumpScalarVariable(MAKEWCHARARR(GEOHEADERID_TOTPATCHES),GlobalScalarVars_.TotNPatches_);
	Zeus::DumpScalarVariable(MAKEWCHARARR(GEOHEADERID_PIXSZ),GlobalScalarVars_.PixSz_);
	Zeus::DumpScalarVariable(MAKEWCHARARR(GEOHEADERID_PRAVOBJINPATCH),GlobalScalarVars_.PriorAvObjectsPatch_);

	(Zeus::ConManager::Instance())->PrintStr2Console(L".\n");
}

void	PlankInfoStore::ReadMapInfo(void)
{
	std::vector<Zeus::VariableIdType>	ids;
	std::vector<Zeus::DBField>			values;
	wchar_t								buffer[DOUBLETXTMAXSZ];
	std::wstring						ID;
	std::wstring						ID1;
	std::wstring						ID2;
	std::wstring						ID3;
	std::wstring						IDFWHM;

	std::wstring						StrTemp;

	for(int i=0;i!=GlobalScalarVars_.N_ObsPlanes_;++i)
	{
		ID					= GLOBAL_OBSPLANEFREQ;
		IDFWHM				= GLOBAL_OBSPLANEFWHM;
		std::wstring		tStr(Zeus::PutNumber2Txt(i));
		ID					+= tStr;
		IDFWHM				+= tStr;		
		ids.push_back(Zeus::VariableIdType(ID,Zeus::DBField()));
		ids.push_back(Zeus::VariableIdType(IDFWHM,Zeus::DBField(static_cast<double>(-1.0))));
	}

	(GlobalVars::Instance())->GetVar(ids,values);

	std::vector<Zeus::DBField>::iterator v_ptr(values.begin());

	double					t1;

	for(int i=0;i!=GlobalScalarVars_.N_ObsPlanes_;++i)
	{

		ID					= GLOBAL_OBSPLANEFREQ;
		IDFWHM				= GLOBAL_OBSPLANEFWHM;
		std::wstring		tStr(Zeus::PutNumber2Txt(i));

		ID		+= tStr;
		IDFWHM	+= tStr;
		ObsFreqsType			t((*v_ptr++).Get<int>());
		GlobalScalarVars_.FreqsColl_.push_back(t);
		t1 = (*v_ptr++).Get<double>();
		ChangeStatic(t.freq_, t1);
		Zeus::DumpScalarVariable(ID,t.sign_);
		Zeus::DumpScalarVariable(IDFWHM,t1);
	}
	
	(Zeus::ConManager::Instance())->PrintStr2Console(L".\n");
}

void		PlankInfoStore::ReadHealpixMaskFile(void)
{
	if(GlobalScalarVars_.MaskFileName_.size() < 2)
		return;

	std::wstring				tDir;
	std::wstring				tFName(Zeus::ExtractFileName(GlobalScalarVars_.MaskFileName_,tDir));
	int							Context;

	if(	GlobalScalarVars_.Data2Buffer_ &&
		tDir.empty() &&
		(GlobalScalarVars_.MaskFileName_ != std::wstring(MC_MASKREJNAME_DEF))
		)
	{
#ifdef WIN32
		std::wstring	MasksExt(L"Masks\\");
#else
		std::wstring	MasksExt(L"Masks/");
#endif
		tDir		= Zeus::RemoveInstanceName(GlobalScalarVars_.DirBuffer_) + MasksExt;
		Context		= 1000;
	}
	else
	{
		Context = GlobalScalarVars_.ContextID_;
		if(tDir.empty())
			tDir = (GlobalScalarVars_.MaskFileName_ != std::wstring(MC_MASKREJNAME_DEF))?GlobalScalarVars_.DirInMasksIn_:GlobalScalarVars_.DirInMasks_;
	}


	std::auto_ptr<Zeus::GenHealpixReader<HEALPIX_ATOM_PREC> >	
		FReader(Zeus::GetGenHealpixReader(Loki::Type2Type<HEALPIX_ATOM_PREC>(),Context,tDir,tFName));

	FReader->Initialize(0);
	FReader->Read();
	FReader->Release(MaskMap_);
}
//
int		PlankInfoStore::RemoveObjectInMask(Zeus::PeakCollType& Peaks)
{

	Zeus::PeakCollType::iterator	newEnd(std::remove_if(Peaks.begin(),Peaks.end(),RemoveObjsInMaskFunctor()));
	if(newEnd != Peaks.end())
	{Peaks.erase(newEnd,Peaks.end());}

	return static_cast<int>(Peaks.size());
}
//
int		PlankInfoStore::AppendDetectObjs2File(int PatchN,Zeus::PeakCollType& Peaks)
{
	if(Peaks.empty())
		return 0;

	TranslatingPix2Pointings(PatchN,Peaks);
	
	if(!(GlobalScalarVars_.NonBlindDetection_))
	{
		RemoveObjectInMask(Peaks);
	}

	TotalN_WrittenObjs_ += ObjWriter_->Write(Peaks);

#if (defined(HFIDMC) || defined(LFIDPC)) && defined(HFIDMC_EXTOBJECTS)
	if(ObjWriterExternals_)
	{
		TotalN_WrittenObjsExt_ += ObjWriterExternals_->Write(Peaks);
	}
#endif

	PrintCurrentStatus(PatchN);

	return static_cast<int>(Peaks.size());
}
//
int		PlankInfoStore::AppendNonValidObjs2File(int SrcIndex,int PatchN)
{
	//Zeus::PatchGeomType
	wchar_t							buffer[BUFFERMAXCHAR];	
	Zeus::PeakCollType				NonValidPeak;
	Zeus::PeakCollType::value_type	temp;

	memset(&temp,0,sizeof(temp));

	temp.DetectID_				= SrcIndex;
	temp.PatchNumber_			= PatchN;
	temp.PK_BayesDetectStat_	= Zeus::PeakCollType::value_type::PK_DET_ERROR;
	temp.JF_lnRho_				= -1.0;
	temp.DetectionSigma_		= -1.0;
	temp.SrcAmplNormalised_		= -1.0;

	NonValidPeak.push_back(temp);

	TotalN_WrittenObjs_ += ObjWriter_->Write(NonValidPeak);

#if (defined(HFIDMC) || defined(LFIDPC)) && defined(HFIDMC_EXTOBJECTS)
	if(ObjWriterExternals_)
	{
		TotalN_WrittenObjsExt_ += ObjWriterExternals_->Write(NonValidPeak);
	}
#endif

	PRINTINTOBUFFERFUNCT
			(buffer,BUFFERMAXCHAR,NONVALIDSTATFORMAT,
			PatchN,
			static_cast<int>(time(NULL) - InitialTime_)
			);
		(Zeus::ConManager::Instance())->PrintStr2Console(buffer);

	return 1;
}
//

std::wstring	PlankInfoStore::GetCurrPatchPtgsOutputName(int patchN,int ptgType) const
{
	std::wstring	ext(ptgType?PHIFILEEXT:THETAFILEEXT);
	
	return  std::wstring(PTGSPREFIX) + Zeus::PutNumber2Txt(GlobalScalarVars_.NSize_) + std::wstring(L"_") + Zeus::PutNumber2Txt(patchN) + ext;
}

void	PlankInfoStore::ReadPtgs(int patchN,Zeus::LArr2D<double>& CoLat,Zeus::LArr2D<double>& Long) const
{
	std::wstring			tDir(GlobalScalarVars_.DirPointings_);

	if(GlobalScalarVars_.Data2Buffer_)
	{
#ifdef WIN32
		std::wstring	PtgsExt(L"Ptgs\\");
#else
		std::wstring	PtgsExt(L"Ptgs/");
#endif
		if(!Zeus::CheckDir(tDir = (GlobalScalarVars_.DirBuffer_  + PtgsExt)))
		{
			tDir = Zeus::RemoveInstanceName(GlobalScalarVars_.DirBuffer_) + PtgsExt;
		}
	}

	std::auto_ptr<Zeus::GenCollReader<Zeus::LArr2D<double> > >
		FReaderColat(Zeus::GetWrkSpFileReaderHandler(Loki::Type2Type<Zeus::LArr2D<double> >(),GlobalScalarVars_.Data2Buffer_?1000:ContextID_,
		GetCurrPatchPtgsOutputName(patchN,0),GlobalScalarVars_.PatchSz_,GlobalScalarVars_.PatchSz_,tDir));
	FReaderColat->Initialize();
	FReaderColat->Read();
	FReaderColat->Release(CoLat);

	std::auto_ptr<Zeus::GenCollReader<Zeus::LArr2D<double> > >
		FReaderLong(Zeus::GetWrkSpFileReaderHandler(Loki::Type2Type<Zeus::LArr2D<double> >(),GlobalScalarVars_.Data2Buffer_?1000:ContextID_,
		GetCurrPatchPtgsOutputName(patchN,1),GlobalScalarVars_.PatchSz_,GlobalScalarVars_.PatchSz_,tDir));
	FReaderLong->Initialize();
	FReaderLong->Read();
	FReaderLong->Release(Long);
}


void	PlankInfoStore::TranslatingPix2Pointings(int PatchN,Zeus::PeakCollType& Peaks)
{
	Zeus::LArr2D<double> CoLat;
	Zeus::LArr2D<double> Long;

	Trafo	TransCoordsEcl(COORDSEPOCH,COORDSEPOCH,Galactic,Ecliptic);
	Trafo	TransCoordsEqu(COORDSEPOCH,COORDSEPOCH,Galactic,Equatorial);

	ReadPtgs(((PlanckInfo::Instance())->GetGeoProps()).Storage_.at(PatchN).SrcIndex_,CoLat,Long);
	Zeus::PeakCollType::iterator				piv(Peaks.begin());
	Zeus::PeakCollType::const_iterator		const end(Peaks.end());
	for(;piv != end;++piv)
	{


		if(GlobalScalarVars_.NotAlignedObjs_)
		{
			double		tYCoord(GlobalScalarVars_.AssessmentKind_?(GlobalScalarVars_.Jf_Estimator_?piv->Pos_.JF_YCoord_.Mean_:piv->Pos_.JF_YCoord_.Mode_):piv->Pos_.YCoord_);
			double		tXCoord(GlobalScalarVars_.AssessmentKind_?(GlobalScalarVars_.Jf_Estimator_?piv->Pos_.JF_XCoord_.Mean_:piv->Pos_.JF_XCoord_.Mode_):piv->Pos_.XCoord_);
			int			YPix(Zeus::toInt(tYCoord + 0.5)),XPix(Zeus::toInt(tXCoord + 0.5));
			double		DeltaY(tYCoord - static_cast<double>(YPix)),DeltaX(tXCoord -static_cast<double>(XPix));
			int			DirY(DeltaY>=0?1:-1),DirX(DeltaX>=0?1:-1);
			vec3		vecBase(pointing(CoLat(YPix,XPix),Long(YPix,XPix)));
			vec3		vec_y(pointing(CoLat(YPix+DirY,XPix),Long(YPix+DirY,XPix)));
			vec3		vec_x(pointing(CoLat(YPix,XPix+DirX),Long(YPix,XPix+DirX)));
			vec3		vec_yx(pointing(CoLat(YPix+DirY,XPix+DirX),Long(YPix+DirY,XPix+DirX)));
			pointing	result(Zeus::BiLinearInterpolation(vecBase,vec_y,vec_x,vec_yx,std::abs(DeltaY),std::abs(DeltaX)));
			piv->GalPt_Colat_	= result.theta;
			piv->GalPt_Long_	= result.phi;
		}
		else
		{
			if(GlobalScalarVars_.AssessmentKind_)
			{
				int tJF_YPix(GlobalScalarVars_.Jf_Estimator_?piv->Pos_.JF_YPix_.Mean_:piv->Pos_.JF_YPix_.Mode_);
				int tJF_XPix(GlobalScalarVars_.Jf_Estimator_?piv->Pos_.JF_XPix_.Mean_:piv->Pos_.JF_XPix_.Mode_);

				piv->GalPt_Colat_	= CoLat(tJF_YPix,tJF_XPix);
				piv->GalPt_Long_	= Long(tJF_YPix,tJF_XPix);		
			}
			else
			{
				piv->GalPt_Colat_	= CoLat(piv->Pos_.YPix_,piv->Pos_.XPix_);
				piv->GalPt_Long_	= Long(piv->Pos_.YPix_,piv->Pos_.XPix_);		
			}
		}
/*
		if(	(MaskMap_.Map().size())	&&
			(MaskMap_[MaskMap_.ang2pix(pointing(piv->GalPt_Colat_, piv->GalPt_Long_))] < 0.5)
			)
		{
			piv->UnderGalMask_ = 1;
		}
		else
		{
			piv->UnderGalMask_ = 0;		
		}
*/
		piv->UnderGalMask_ = 0;

		pointing tempEcl(TransCoordsEcl(pointing(piv->GalPt_Colat_,piv->GalPt_Long_)));
		pointing tempEqu(TransCoordsEqu(pointing(piv->GalPt_Colat_,piv->GalPt_Long_)));

		piv->CoordPt_Lat_	= (PIOVER2 - tempEcl.theta) * RAD2DEGREE;
		piv->CoordPt_Long_  = tempEcl.phi * RAD2DEGREE;

		piv->EquPt_Lat_	= (PIOVER2 - tempEqu.theta) * RAD2DEGREE;
		piv->EquPt_Long_  = tempEqu.phi * RAD2DEGREE;

		++TotalN_Objs_;
		if( (!(GlobalScalarVars_.SZ_) && (piv->SrcFlux_mJys_ >= GlobalScalarVars_.OutputFluxThreshold_)) ||
			((GlobalScalarVars_.SZ_) && (piv->SrcCompt_arcmin2_ >= GlobalScalarVars_.OutputFluxThreshold_))
			)
		{
			++TotalN_BrightObjs_;
		}
	}
}

void	PlankInfoStore::TranslatingPix2PointingsLib(int PatchN,int NPixCoords,unsigned int Xpix[],unsigned int Ypix[],float Colatitude[],float Longitude[])
{
	Zeus::LArr2D<double> CoLat;
	Zeus::LArr2D<double> Long;
	
	ReadPtgs(PatchN,CoLat,Long);

	for(int i=0;i<NPixCoords;++i)
	{
		Colatitude[i]	= static_cast<float>(CoLat(Ypix[i],Xpix[i]));
		Longitude[i]	= static_cast<float>(Long(Ypix[i],Xpix[i]));
	}
}

int		PlankInfoStore::DoCompressCatalogue(CommandLineArgs& args,double frequency,double SZ_spectralCte)
{
	std::auto_ptr<CompressCatalogue> CatCompressor(new CompressCatalogue(ContextID_,args.CatalogueFNameColl_,args.FinalCatFName_,MaskMap_));
	if(!(CatCompressor->Initialize()))
		return 0;
	CatCompressor->Do_CompressCatalogue();
	CatCompressor->SetOutputFields(frequency,SZ_spectralCte);
	int NSrcs(CatCompressor->WriteOutCatalogue());

//
#if defined(HFIDMC) || defined(LFIDPC)
	wchar_t			buffer[DOUBLETXTMAXSZ];
	std::wstring	tName(L"pws_DetectionMutex_");
	swprintf(buffer,DOUBLETXTMAXSZ,L"%d",GlobalScalarVars_.Sync_ID_);
	tName += std::wstring(buffer);
	Healpix_Map<HEALPIX_ATOM_PREC>	temp(2,RING,nside_dummy());
	Zeus::CreatReadableHealpixFile(ContextID_,GlobalScalarVars_.DirInMasks_,tName,Galactic,temp);

#endif
//
	return NSrcs;

}

void	PlankInfoStore::ProcessScalesColl(void)
{
	if(GlobalScalarVars_.ScalesColl_.empty())
		return;

	FilterScalesCollType::iterator			piv(GlobalScalarVars_.ScalesColl_.begin());
	FilterScalesCollType::const_iterator	const end(GlobalScalarVars_.ScalesColl_.end());
	double									t;

	FilterScalesCollType					tList;
	
	GlobalScalarVars_.ScalesMaxNumber_		= -1;
	GlobalScalarVars_.ScalesFillerType_		= Zeus::toInt(*piv); ++piv;

	switch(GlobalScalarVars_.ScalesFillerType_)
	{
//
	case 1:			// Log
		GlobalScalarVars_.ScalesFirstElement_	= *piv; ++piv;
		for(;piv != end;++piv)
		{
			t = ((static_cast<double>(*piv) * static_cast<double>(GlobalScalarVars_.N_ScaleBins_)) / (double)100.0);
			if(t > static_cast<double>(GlobalScalarVars_.N_ScaleBins_))
				t = static_cast<double>(GlobalScalarVars_.N_ScaleBins_);
			tList.push_back(t);
		}			
		break;
//
	case 2:   // Hardcoded
		for(;piv != end;++piv)
		{tList.push_back(*piv);}

		GlobalScalarVars_.N_ScaleBins_ = static_cast<int>(tList.size());
		GlobalScalarVars_.SrcMaxScale_ = tList[tList.size()-1];
		GlobalScalarVars_.PriorSrcMaxScale_ = -1;
		break;
//
	case 3:   // Linear
		for(;piv != end;++piv)
		{
			t = ((static_cast<double>(*piv) * static_cast<double>(GlobalScalarVars_.N_ScaleBins_)) / (double)100.0);
			if(t > static_cast<double>(GlobalScalarVars_.N_ScaleBins_))
				t = static_cast<double>(GlobalScalarVars_.N_ScaleBins_);
			tList.push_back(t);
		}			
		break;
//
	default: // Linear but with fixed scales
		GlobalScalarVars_.ScalesMaxNumber_		= Zeus::toInt(*piv); ++piv;
		for(;piv != end;++piv)
		{tList.push_back(*piv);}			
		break;
	}
	GlobalScalarVars_.ScalesColl_.swap(tList);
}
//
void		PlankInfoStore::ComputeStats1Map(const Zeus::RealPlane<double>& ws,double rejectLevel,double& bias,double& rms,int& PixIncluded,double guess) const
{
	const   Zeus::RealPlane<double>::DataInnerType&	InnerArray(ws.GetInnerData());
	int		NPixIncluded(0);
	double	Bias(0.0);
	double	Variance(0.0);

	Zeus::LArr2D<double>::const_iterator	const pixBeg(InnerArray.begin());
	Zeus::LArr2D<double>::const_iterator	pixPiv(pixBeg);
	Zeus::LArr2D<double>::const_iterator	const pixEnd(pixPiv + InnerArray.getSz());
	guess *= DEFREJECTLEV_NPMAP;
//TODO OpenMP
	for(;pixPiv != pixEnd;++pixPiv)
	{
		if(((guess >= 0.0) && (std::abs(*pixPiv) > guess)) || IsPixInBorder(static_cast<long>(pixPiv - pixBeg)))
			continue;
		Bias		+= *pixPiv;
		Variance	+= (*pixPiv * *pixPiv);
		++NPixIncluded;
	}

	Variance -= ((Bias * Bias) / static_cast<double>(NPixIncluded));
	double tRms(std::sqrt(Variance / static_cast<double>(NPixIncluded)));
	Bias /= static_cast<double>(NPixIncluded);


	int		tNPixIncluded(0);
	double	tBias(0.0);
	double	tVariance(0.0);

	pixPiv  = InnerArray.begin();
//TODO OpenMP
	for(;pixPiv != pixEnd;++pixPiv)
	{
		if((std::abs(*pixPiv - Bias) < (rejectLevel * tRms))  && !(IsPixInBorder(static_cast<long>(pixPiv - pixBeg))))
		{
			++tNPixIncluded;
			tBias			+= *pixPiv;
			tVariance		+= (*pixPiv * *pixPiv);
		}
	}

	tVariance		-= ((tBias * tBias) / static_cast<double>(tNPixIncluded)); 
	rms				=  std::sqrt(tVariance / static_cast<double>(tNPixIncluded));
	bias			=  tBias / static_cast<double>(tNPixIncluded);
	PixIncluded		=  tNPixIncluded;
}
//
void		PlankInfoStore::doTruncateOutliers(Zeus::RealPlane<double>& BackgReal,double bias,double rms, MaskingType mask)
{
	Zeus::RealPlane<double>::DataInnerType&	InnerArray(BackgReal.GetInnerData());

	Zeus::LArr2D<double>::const_iterator	const pixOrg(InnerArray.begin());
	Zeus::LArr2D<double>::iterator			pixPiv(InnerArray.begin());
	Zeus::LArr2D<double>::const_iterator	const pixEnd(pixPiv + InnerArray.getSz());
	int			Metric,YSz,inborders,toff;
	double		SigLevel;

	BackgReal.GetSz(YSz,Metric);

	for(;pixPiv != pixEnd;++pixPiv)
	{
		SigLevel = (*pixPiv - bias) / rms;
		if(SigLevel <= STATREJECTMASKTHRESHOLD)
			continue;
		toff	= (static_cast<int>(pixPiv - pixOrg));
		inborders = IsPixInBorder(toff);
		if(
			(mask != MASK_ALL)
			&& (((mask == MASK_CORE)   && inborders) || ((mask == MASK_BORDER) && !inborders))
			)
			continue;
		*pixPiv	= bias + (STATREJECTMASKTHRESHOLD * rms);
	}
}
//
