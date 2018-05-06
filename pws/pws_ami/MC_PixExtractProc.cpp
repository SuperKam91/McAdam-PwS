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

//--------------------------------

#include "alm.h"
#include "healpix_map.h"
#include "alm_powspec_tools.h"
#include "alm_healpix_tools.h"
#include "healpix_map_fitsio.h"
#include "MC_HPixMapPatchCutting.h"
#include "MC_HealpixCutter.h"
#include "MC_PixExtractProc.h"
#include "MC_AntennaGaussian.h"
#include "MC_AntennaPixWindow.h"
#include "ZEUS_Debug.h"

#define MASKADDNOISELEVEL	3.0

//--------------------------------

PlanckMapsInfoCollType					PixExtractProc::MapStaticInfo_;
PlanckMapsInfoCollType::const_iterator	PixExtractProc::PlanckEnd_;

void			PixExtractProc::SetMapStaticInfo(void)
{
	PlanckMapsInfoType temp;
	MapStaticInfo_.clear();
	temp.Freq_			= 30;
	temp.EffFreq_		= 30.0;
	temp.AntFWHM_		= 32.628;
	temp.OperatNSide_	= 1024; //512
	MapStaticInfo_.push_back(temp);
	temp.Freq_			= 44;
	temp.EffFreq_		= 44.0;
	temp.AntFWHM_		= 27.958;
	temp.OperatNSide_	= 1024; //512
	MapStaticInfo_.push_back(temp);
	temp.Freq_			= 70;
	temp.EffFreq_		= 70.0;
	temp.AntFWHM_		= 13.046;
	temp.OperatNSide_	= 1024; //1024
	MapStaticInfo_.push_back(temp);
	temp.Freq_			= 100;
	temp.EffFreq_		= 100.0;
	temp.AntFWHM_		= 9.88;
	temp.OperatNSide_	= 2048; //1024
	MapStaticInfo_.push_back(temp);
	temp.Freq_			= 143;
	temp.EffFreq_		= 143.0;
	temp.AntFWHM_		= 7.18;
	temp.OperatNSide_	= 2048;
	MapStaticInfo_.push_back(temp);
	temp.Freq_			= 217;
	temp.EffFreq_		= 217.0;
	temp.AntFWHM_		= 4.87;
	temp.OperatNSide_	= 2048;
	MapStaticInfo_.push_back(temp);
	temp.Freq_			= 353;
	temp.EffFreq_		= 353.0;
	temp.AntFWHM_		= 4.650;
	temp.OperatNSide_	= 2048;
	MapStaticInfo_.push_back(temp);
	temp.Freq_			= 545;
	temp.EffFreq_		= 545.0;
	temp.AntFWHM_		= 4.720;
	temp.OperatNSide_	= 2048;
	MapStaticInfo_.push_back(temp);
	temp.Freq_			= 857;
	temp.EffFreq_		= 857.0;
	temp.AntFWHM_		= 4.390;
	temp.OperatNSide_	= 2048;
	MapStaticInfo_.push_back(temp);
	PlanckEnd_ = MapStaticInfo_.end();
}
//
void			PixExtractProc::ReadBadPixMask(const std::wstring& DirName)
{
	std::wstring	ActualFileName;
	(Zeus::ConManager::Instance())->PrintStr2Console(std::wstring(MC_MASKREADINGBADPIX) + std::wstring(MC_MASKMASKALERT) + std::wstring(L"\n\n"));
	
	int fileExits(ReadMaskMap(DirName,MaskRemoveName_,CutOutMask_,ActualFileName));
	
	if((MaskRemoveName_.size() >= 2) && !fileExits)
	{errHealpix(ERRCOD_MC_MAPNOTFOUND,ERRMSG_MC_MAPNOTFOUND,std::wstring(L""),ActualFileName,1);}

	if(CutOutMask_.Map().size() && (CutOutMask_.Scheme() != RING))
	{CutOutMask_.swap_scheme();}	
}
//
void			PixExtractProc::ReadSrcCatalogue(const std::wstring& SrsCatFName,Zeus::NonBlingCatType& SourceCat)
{

	std::wstring	tDir;
	std::wstring	tFName(Zeus::ExtractFileName(SrsCatFName,tDir));
	int				Context;

	if(Data2Buffer_ && tDir.empty())
	{
#ifdef WIN32
		std::wstring	CatsExt(L"Cats\\");
#else
		std::wstring	CatsExt(L"Cats/");
#endif
		tDir = DataBuffer_ + CatsExt;
		Context = 1000;
	}
	else
	{
		Context = InOutEnvID_;
		if(tDir.empty())  tDir = DirOut_;
	}

	Zeus::CatalogueFormatType	tCat;

	std::auto_ptr<Zeus::GenCollReader<Zeus::CatalogueFormatType> >
		FReader(Zeus::GetGenCollFileReaderHandler(Loki::Type2Type<Zeus::CatalogueFormatType>(),Context,tDir , tFName));

	FReader->Initialize();
	FReader->Read();
	FReader->Release(tCat);

	SourceCat.Header_.CoordsType_		= tCat.Header_.CoordsType_;
	SourceCat.Header_.CoordSystem_		= tCat.Header_.CoordSystem_;
	SourceCat.Header_.Epoch_			= tCat.Header_.Epoch_;

	Zeus::CatalogueFormatType::StorageType::const_iterator	OrgPiv(tCat.Storage_.begin());
	Zeus::CatalogueFormatType::StorageType::const_iterator	const OrgEnd(tCat.Storage_.end());
	Zeus::NonBlingCatType::StorageType::value_type			t;

	SourceCat.Storage_.clear();

	for(;OrgPiv != OrgEnd; ++OrgPiv)
	{
		t.ptg_.phi				= OrgPiv->GalLongDegs_;
		t.ptg_.theta			= OrgPiv->GalLatDegs_;
		t.ErrPosition_			= OrgPiv->ErrorBars_.TotalPosErrorBar_;
		t.PredFluxGLRT_			= OrgPiv->FluxComptGLRT_;
		t.PredRadiusGLRT_		= OrgPiv->RadiusGLRT_;
//		t.SNR_					= ((OrgPiv->DetectSigma_ < OrgPiv->NormalAmpl_) ? OrgPiv->DetectSigma_ : OrgPiv->NormalAmpl_) ;
		t.SNR_					= OrgPiv->NormalAmpl_;
		t.PredFluxBay_			= OrgPiv->FluxCompt_;
		t.PredRadiusBay_		= OrgPiv->Radius_;
		t.ErrFlux_				= OrgPiv->ErrorBars_.FluxErrorBar_;
		t.ErrRadius_			= OrgPiv->ErrorBars_.RadiusErrorBar_;

		SourceCat.Storage_.push_back(t);
	}
}
//
void			PixExtractProc::FillInDefaultValues(PixExtractProcCtorArgsType::MapFileInfoCollType& MapColl)
{
	PixExtractProcCtorArgsType::MapFileInfoCollType::iterator		piv(MapColl.begin());
	PixExtractProcCtorArgsType::MapFileInfoCollType::const_iterator	const end(MapColl.end());
	PlanckMapsInfoCollType::const_iterator							StaticInfo;

	for(;piv != end;++piv)
	{
		if((piv->Freq_ < 0.0) || (piv->Ant_FWHM_ < 0.0))
		{
			StaticInfo = GetPlanckData(std::abs(piv->FreqID_) % 1000);
		}
		if(piv->Freq_ < 0.0)
		{
			piv->Freq_ = StaticInfo->EffFreq_;
		}
		if(piv->Ant_FWHM_ < 0.0)
		{
			piv->Ant_FWHM_ = StaticInfo->AntFWHM_;
		}	
	}
}
//
PlanckMapsInfoCollType::const_iterator	PixExtractProc::GetPlanckData(int freq)
{
	PlanckMapsInfoCollType::const_iterator	piv(MapStaticInfo_.begin());
	PlanckMapsInfoCollType::const_iterator	const end(MapStaticInfo_.end());

	for(;piv != end;++piv)
	{if(piv->Freq_ == freq) return piv;}
	errInvalidFreq(freq);
	return end;
}
//
void			PixExtractProc::errHealpixNside(int errCode,wchar_t* msg, int NSideGeo,int NSideMap) const
{
	std::wstring errstring(msg);

	errstring += std::wstring(HEALPIX_GEOMNSIDE_MSG);
	errstring += Zeus::PutNumber2Txt(NSideGeo);
	errstring += std::wstring(HEALPIX_MAPNSIDE_MSG);
	errstring += Zeus::PutNumber2Txt(NSideMap);
	errGeomFile(errCode,errstring.c_str());
}
//
void			PixExtractProc::GetMapUnits(const std::string& unitStr)
{
	std::string tstr(Zeus::ToCase(unitStr,false));
	if(unitStr.empty())
	{
		errHealpix(ERRCOD_MC_HPINVALIDTUNIT1,ERRMSG_MC_HPINVALIDTUNIT1,HP_Dir_,CurMapPiv_->HP_FName_,0);
		return;
	}
	if(tstr.find(HEALPIX_UNIT_BRIGHT)!= std::string::npos)
	{
		CurMapPiv_->MapUnits_ = Zeus::PlanckUnitsValuesTranform::BRIGHTNESS;
		return;
	}
	if(tstr.find(HEALPIX_UNIT_THERMO)!= std::string::npos)
	{
		CurMapPiv_->MapUnits_ = Zeus::PlanckUnitsValuesTranform::THERMO_T;
	}
	else if(tstr.find(HEALPIX_UNIT_ANTENN)!= std::string::npos)
	{CurMapPiv_->MapUnits_ = Zeus::PlanckUnitsValuesTranform::ANTENNA_T;}
	else
	{
		errHealpix(ERRCOD_MC_HPINVALIDTUNIT1,ERRMSG_MC_HPINVALIDTUNIT1,HP_Dir_,CurMapPiv_->HP_FName_,0);
		return;
	}
	return;
}
//
void			PixExtractProc::ProcCoordSys(const Zeus::HealpixHeaderAtomType& key)
{
	if(key.ValueType_ == Zeus::NO_KEY)
	{
		errHealpix(ERRCOD_MC_HPNOTCOORDYS,ERRMSG_MC_HPNOTCOORDYS,HP_Dir_,CurMapPiv_->HP_FName_,0);
	}
	else
	{
		(Zeus::ConManager::Instance())->PrintWarning2Console(Zeus::Achar2Wstr(key.stringValue_.c_str()));
		int t(Zeus::GetCoordSysFromStr(key.stringValue_));
		if(t>=0) CurMapPiv_->MapCoordSys_ = t;
	}
	PrintCoordSys(CurMapPiv_->MapCoordSys_);
	if(Galactic != CurMapPiv_->MapCoordSys_)
	{
		PrintCoordSysMismatch(CurMapPiv_->MapCoordSys_);
		Trafo_ = new Trafo(Epoch_,Epoch_,Galactic,static_cast<coordsys>(CurMapPiv_->MapCoordSys_));
	}
	else{delete Trafo_;Trafo_ = 0;}
}

//
void			PixExtractProc::ProcUnits(const Zeus::HealpixHeaderAtomType& key)
{
	if(key.ValueType_ == Zeus::NO_KEY)
	{
		errHealpix(ERRCOD_MC_HPNOTTUNIT1,ERRMSG_MC_HPNOTTUNIT1,HP_Dir_,CurMapPiv_->HP_FName_,0);
	}
	else
	{
		(Zeus::ConManager::Instance())->PrintWarning2Console(Zeus::Achar2Wstr(key.stringValue_.c_str()));
		GetMapUnits(key.stringValue_);
	}
	
	PrintMapUnits();

	UnitsTransf_.Set(CurMapPiv_->MapUnits_,Zeus::PlanckUnitsValuesTranform::BRIGHTNESS,CurMapPiv_->Freq_,CurMapPiv_->UnitFactor_);

}
//
void			PixExtractProc::ReplaceBadPixSubBias(void)
{
	Zeus::LArr2D<double>::iterator			pixPiv(Pix_.begin());
	Zeus::LArr2D<double>::const_iterator	const pixEnd(Pix_.end());

	for(;pixPiv != pixEnd;++pixPiv)
	{
		if((*pixPiv) < HEALPIX_BADPIX)
		{
			*pixPiv	= 0.0;
		}
		else
		{
			*pixPiv	-= CurrStat_.CurrPatchBias_;
		}
	}
}
//
void			PixExtractProc::AddNoise2Pix(Zeus::LArr2D<double>& Pixels,Zeus::LArr2D<double>& PixMask)
{
	Zeus::LArr2D<double>::iterator			pivPix(Pixels.begin());
	Zeus::LArr2D<double>::const_iterator	const pixEnd(pivPix + Pixels.getSz());
	Zeus::LArr2D<double>::const_iterator	pivMask(PixMask.begin());


	double	NoiseLevel(CurrStat_.CurrPatchRms_ / MASKADDNOISELEVEL);

	for(;pivPix != pixEnd;++pivMask,++pivPix)
	{
		if((*pivMask) < HEALPIX_BADPIX)
		{
			*pivPix = (NoiseLevel * NoiseGen_->RandGauss());
		}
	}
}
//
void			PixExtractProc::NormaliseRMS(Zeus::LArr2D<double>& PixDest)
{
	double	Bias(0.0);
	double	Variance(0.0);
	const	int NPixIncluded(PixDest.getSz());

	Zeus::LArr2D<double>::iterator			pixPiv(PixDest.begin());
	Zeus::LArr2D<double>::const_iterator	const pixEnd(pixPiv + PixDest.getSz());

	for(;pixPiv != pixEnd;++pixPiv)
	{
		Bias		+= *pixPiv;
		Variance	+= (*pixPiv * *pixPiv);
	}
	Variance -= ((Bias * Bias) / static_cast<double>(NPixIncluded));
	double tRms(std::sqrt(Variance / static_cast<double>(NPixIncluded)));
	Bias /= static_cast<double>(NPixIncluded);
	const double NoiseLevel(CurrStat_.CurrPatchRms_/tRms);

	pixPiv = PixDest.begin();
	for(;pixPiv != pixEnd;++pixPiv)
	{
		*pixPiv = ((*pixPiv - Bias) * NoiseLevel) + Bias;
	}
}
//
double			PixExtractProc::TryBestRotation(PatchDirections& BestRotation,const Zeus::LArr2D<double>& Mask)
{
	double			RejectRatio;
	double			MinRejectPercent(1e30);
	int				initialRot(BestRotation),i(BestRotation);
	
	do
	{
		RejectRatio = PixTry1Rotation(static_cast<PatchDirections>(i),Mask);
		if(RejectRatio < MinRejectPercent)
		{
			MinRejectPercent	= RejectRatio;
			BestRotation		= static_cast<PatchDirections>(i);
			if(MinRejectPercent == 0.0)
				return 0.0;
		}
	++i;i %= 4;
	}while(i!= initialRot);

	return MinRejectPercent;
}
//
double			PixExtractProc::PixTry1Rotation(PatchDirections dir,const Zeus::LArr2D<double>& Mask)
{
	const int	CentralPix(PtchSz_ >> 1);
	int			treject(0);
	int			New_j,New_i;

	for(int j= -CentralPix; j<CentralPix; ++j)
	{
		for(int i= -CentralPix; i<CentralPix; ++i)
		{
			if(Mask(j + CentralPix,i + CentralPix) < 0.5)
			{
				switch(dir)
				{
				case	PTCHDIR_NAT:
					New_j = j;
					New_i = i;
					break;
				case	PTCHDIR_ROT90:
					New_j =  i;
					New_i = -j;
					break;
				case	PTCHDIR_ROT180:
					New_j = -j;
					New_i = -i;
					break;
				case	PTCHDIR_ROT270:
					New_j = -i;
					New_i =  j;
					break;
				default:
					New_j = j;
					New_i = i;
					break;
				}

				if(New_j >= CentralPix ) New_j = CentralPix - 1;
				if(New_i >= CentralPix ) New_i = CentralPix - 1;

				New_j += CentralPix;
				New_i += CentralPix;

				if(Mask(New_j,New_i) < 0.5)
				{
					++treject;
				}
			}
		}
	}

	return static_cast<double>(treject) / static_cast<double>(Ptgs_.getSz());
}

//
int				PixExtractProc::CheckBadPix(void)
{
	Zeus::LArr2D<double>::const_iterator	pivPixOrg(Pix_.begin());
	Zeus::LArr2D<double>::const_iterator	const endPixOrg(pivPixOrg + Pix_.getSz());

	for(;pivPixOrg != endPixOrg;++pivPixOrg)
	{
		if(*pivPixOrg == Healpix_undef)
			return 0;
	}
	return 1;
}
//
void			PixExtractProc::MakeMask(Zeus::LArr2D<double>& NewMask)
{
	Zeus::LArr2D<double>::const_iterator	pivPixOrg(Pix_.begin());
	Zeus::LArr2D<double>::const_iterator	const endPixOrg(pivPixOrg + Pix_.getSz());
	Zeus::LArr2D<double>::iterator			pivMask(NewMask.begin());

	for(;pivPixOrg != endPixOrg;++pivPixOrg,++pivMask)
	{
		if(*pivPixOrg == Healpix_undef)
		{*pivMask = 0.0;}
		else{*pivMask = 1.0;}
	}
}
//
double			PixExtractProc::CombinePatchesX(Zeus::LArr2D<double>& NewMap,Zeus::LArr2D<double>& Mask)
{
	Zeus::LArr2D<double>::iterator			pivPixOrg(Pix_.begin());
	Zeus::LArr2D<double>::const_iterator	const endPixOrg(pivPixOrg + Pix_.getSz());
	Zeus::LArr2D<double>::iterator			pivNewMap(NewMap.begin());
	Zeus::LArr2D<double>::iterator			pivMask(Mask.begin());
	int										BadPixels(0);

	for(;pivPixOrg != endPixOrg;++pivPixOrg,++pivMask,++pivNewMap)
	{
		if((*pivPixOrg == Healpix_undef) && (*pivNewMap == Healpix_undef))
		{
			*pivPixOrg	= Healpix_undef;
			*pivMask	= 0.0;
			++BadPixels;
			continue;
		}

		if((*pivNewMap == Healpix_undef) || (*pivMask >= 1.0))
		{
			goto CombinePatchesX_loop;
		}

		if((*pivPixOrg == Healpix_undef) || (*pivMask <= 0.0))
		{
			*pivPixOrg	= *pivNewMap;
			goto CombinePatchesX_loop;
		}

		*pivPixOrg	= ((*pivPixOrg) * (*pivMask) + (*pivNewMap) * (1.0 - *pivMask));

CombinePatchesX_loop:
		*pivMask	= 1.0;
	}

	return static_cast<double>(BadPixels)/static_cast<double>(Pix_.getSz());
}
//
void			PixExtractProc::RenormaliseMask(Zeus::LArr2D<double>& DestPix)
{
	Zeus::LArr2D<double>::iterator			pivPixOrg(DestPix.begin());
	Zeus::LArr2D<double>::const_iterator	const endPixOrg(pivPixOrg + DestPix.getSz());
	double	Average(0.0);
	double	Max(-1.0e38);

	for(;pivPixOrg != endPixOrg;++pivPixOrg)
	{
		if(*pivPixOrg > Max)
			Max = *pivPixOrg;
		Average	+=	*pivPixOrg;
	}
	Average	/=	static_cast<double>(DestPix.getSz());
	//Max		-= Average;
	for(pivPixOrg=DestPix.begin();pivPixOrg != endPixOrg;++pivPixOrg)
	{
		//(*pivPixOrg) -= Average;
		if(*pivPixOrg < 0.0) *pivPixOrg = 0.0;
		*pivPixOrg	/=	Max;	
	}
}
//
double			PixExtractProc::MakeNewPatchByRotation(void)
{
	Zeus::LArr2D<double>	Mask(Pix_.getSz(),Pix_.getPtrMetric());
	Zeus::LArr2D<double>	RotatedPatch(Pix_.getSz(),Pix_.getPtrMetric());
	double					UncoveredArea;

	MakeMask(Mask);
	PatchDirections	CurrDirection;
	for(int	i=PTCHDIR_ROT180;i!=PTCHDIR_ROTEND;++i)
	{
		CurrDirection = static_cast<PatchDirections>(i);
		UncoveredArea = TryBestRotation(CurrDirection,Mask);
		if(CurrDirection==PTCHDIR_NAT)
			return UncoveredArea;

		SmoothMask(Mask);

		PixDoRotation(CurrDirection,RotatedPatch);
		if(CombinePatchesX(RotatedPatch,Mask) == 0.0)
			return 0.0;
	}
	return UncoveredArea;
}
//
void			PixExtractProc::SmoothMask(Zeus::LArr2D<double>& Mask)
{
	Zeus::LArr2D<double>	DestPix;

	AntennaSmooth(Mask,DestPix,CurMapPiv_->Freq_,MaskSmoothFWHM_ * SMOOTHFACTOR);
	RenormaliseMask(DestPix);
	Mask.swap(DestPix);
}
//

//
void	PixExtractProc::DoPixExtraction1Map(long LPatch)
{
	if((LPatch < 0) || (LPatch > (TotalPatches_-1)))
	{
		LPatch = TotalPatches_-1;
	}

	while(CurrPatchNumber_ <= LPatch)
	{
		char	bufferS[BUFFERMAXCHAR];	
		if(!(PatchGeom_.Storage_.at(CurrPatchNumber_).PatchValid_))
		{
			ReadPtgsFile(PatchGeom_.Storage_.at(CurrPatchNumber_).SrcIndex_);
			Ptgs2Pixels(IllPixMask_);
			CurrStat_.Reset();
			if(!CheckBadPix())
			{
				if((CurrStat_.PercentBadpixAfterRot_ =  MakeNewPatchByRotation())!= 0.0)
				{
					Zeus::LArr2D<double>	DestPix;
					Zeus::LArr2D<double>	temp;
					Zeus::LArr2D<double>	Mask(Pix_.getSz(),Pix_.getPtrMetric());

					MakeMask(Mask);
/*
					if(CurrStat_.PercentBadpixAfterRot_ > 10.0)
					{
						sprintf(bufferS,"%s_%d_","BadpixMask",CurrPatchNumber_);
						Zeus::DumpInOut_2d(bufferS,256,256,512,128*512+128,Mask.begin(),1.0);
					}
*/
					SmoothMask(Mask);
					ComputeStats();
					temp.Assign(Pix_);
					ReplaceBadPixSubBias();
					AntennaSmooth(Pix_,DestPix,CurMapPiv_->Freq_,MaskSmoothFWHM_ * SMOOTHFACTOR);
					NormaliseRMS(DestPix);
					AddNoise2Pix(DestPix,temp);
					CombinePatchesX(DestPix,Mask);
				}
				else
				{
					ComputeStats();
					ReplaceBadPixSubBias();			
				}
			}
/*		
			if(CurrStat_.PercentBadpixAfterRot_ > 10.0)
			{
				sprintf(bufferS,"%s_%d_","Patches",CurrPatchNumber_);
				Zeus::DumpInOut_2d(bufferS,256,256,512,128*512+128,Pix_.begin(),1.0);
			}
*/
			if (CurMapPiv_->FreqID_ < 0)
			{
				InvertPatch();
			}
			ComputeStats();		
			WritePixels2File(Pix_);
			ReportPixCurrentParams();
		}
		else
		{
			wchar_t	buffer[BUFFERMAXCHAR];

			PRINTINTOBUFFERFUNCT (buffer,BUFFERMAXCHAR,MC_REJECTPATCHSTR,PatchGeom_.Storage_.at(CurrPatchNumber_).SrcIndex_);
			(Zeus::ConManager::Instance())->PrintStr2Console(buffer);			
		}
		++CurrPatchNumber_;
	}
}
//
void	PixExtractProc::AntennaSmooth(const Zeus::LArr2D<double>& OrgPix,Zeus::LArr2D<double>& DestPix,double Freq,double FWHM)
{
	Zeus::FourierPlane<std::complex<double> > Spectrum(FourierMachine_->GetXComplexBufferSz() - 1,true);
	FourierMachine_->Real2Fourier(OrgPix,Spectrum.GetInnerData());
	AntennasCollType::iterator	ptrAnt(CreateAntennas(MC_antenna::ANTENNA_GAUSSIAN,Freq,FWHM));
	Zeus::MultiplyInPlace(Spectrum,ptrAnt->Antenna_->GetAntennaBufferRef(),Zeus::UB_NOUSE);
	DestPix.Make(OrgPix.getSz(),OrgPix.getPtrMetric());
	FourierMachine_->Fourier2Real(Spectrum.GetInnerData(),DestPix);
	NormalisePix(1.0/static_cast<double>(PtchSz_*PtchSz_),DestPix);
}

void	PixExtractProc::ReadGeoProps(int JustHeader)
{
	wchar_t buffer[DOUBLETXTMAXSZ];

	std::auto_ptr<Zeus::GenCollReader<Zeus::PatchGeomType> >
		FReader(Zeus::GetPatchGeomFileReaderHandler(Loki::Type2Type<Zeus::PatchGeomType>(),InOutEnvID_,DirOut_,std::wstring(PATCH_GEOMPROPSFNAME)));

	Zeus::PatchGeomType::HeaderType * hdPtr(FReader->Initialize());

	TotalPatches_	= hdPtr->NTotalPatches_;
	PtchSz_			= hdPtr->PtchSz_;
	GeomNSide_		= hdPtr->NSide_;
	BorderSz_		= hdPtr->PtchBorder_;
	PixSz_			= std::sqrt(PITIMES4 / (12.0 * static_cast<double>(GeomNSide_) * static_cast<double>(GeomNSide_)));
	
	if(!JustHeader)
	{
		FReader->Read();
		FReader->Release(PatchGeom_);
	}
}

void	PixExtractProc::SetNoiseGenerator(void)
{
	delete NoiseGen_;
	NoiseGen_ = new Zeus::PwSCoreRandGen(-1);
	NoiseGen_->RandDouble();
}

int		PixExtractProc::GetMaxNSide(void)
{
	int MaxNSide(0);
	PlanckMapsInfoCollType::const_iterator	tpiv;
	PixExtractProcCtorArgsType::MapFileInfoCollType::const_iterator
		piv(FileCollection_.begin());
	PixExtractProcCtorArgsType::MapFileInfoCollType::const_iterator
		const end(FileCollection_.end());

	for(;piv != end;++piv)
	{
		if((tpiv = GetPlanckData(std::abs(piv->FreqID_)%1000)) == PlanckEnd_)
			continue;
		if(tpiv->OperatNSide_ > MaxNSide) MaxNSide = tpiv->OperatNSide_;
	}
	return MaxNSide;
}

void	PixExtractProc::Initialize(void)
{
	int MaxNSide;

	GeomNSide_ = 2048;

	if(!OpMode_)
	{
		ReadGeoProps(0);

		if((MaxNSide = GetMaxNSide())> GeomNSide_)
			errHealpixNside(ERRCOD_MC_HPNSIDEMISMATCH,ERRMSG_MC_HPNSIDEMISMATCH,GeomNSide_,MaxNSide);

		Ptgs_.Make(PtchSz_ * PtchSz_,PtchSz_);
		InitFourierMachine();
		SetNoiseGenerator();
	}

	if(OpMode_)
	{
		ReadBadPixMask(Masks_DirIn_);
	}
	else
	{	
		(Zeus::ConManager::Instance())->PrintStr2Console(std::wstring(MC_MASKREADINGBADPIX) + std::wstring(L"\n\n"));

		std::wstring	Dummy0(MC_MASKBADPIX_DEF);
		std::wstring	Dummy1;

		// Pass dir together with the name
		ReadMaskMap(std::wstring(L""),Zeus::CorrectDir(Masks_Dir_,1) + Dummy0,IllPixMask_,Dummy1);
	}
}


void	PixExtractProc::PrintMapUnits(void) const
{
	std::wstring	msg(HEALPIX_UNITS_MSG);
	wchar_t buffer[DOUBLETXTMAXSZ];

	switch(CurMapPiv_->MapUnits_)
	{
	case Zeus::PlanckUnitsValuesTranform::BRIGHTNESS :
		msg += std::wstring(HEALPIX_UNITS_BRIGHT);
		break;
	case Zeus::PlanckUnitsValuesTranform::THERMO_T :
		msg += std::wstring(HEALPIX_UNITS_THERMO);
		break;
	case Zeus::PlanckUnitsValuesTranform::ANTENNA_T :
		msg += std::wstring(HEALPIX_UNITS_ANTENN);
		break;
	}
	
	msg	+= std::wstring(L" ; ");
	msg	+= (std::wstring(HEALPIX_UNITSSCAL_MSG) + Zeus::PutNumber2Txt(CurMapPiv_->UnitFactor_));
	msg	+= std::wstring(L"\n\n");

	(Zeus::ConManager::Instance())->PrintStr2Console(msg);

}

void	PixExtractProc::PrintObsMapLowRes(int CurrNSide,int OperaNSide) const
{
	std::wstring msg(HEALPIX_LOWRES_WARNING);
	msg += std::wstring(HEALPIX_CURRNSIDE);
	msg += Zeus::PutNumber2Txt(CurrNSide);
	msg += std::wstring(HEALPIX_NOMINALNSIDE);
	msg += Zeus::PutNumber2Txt(OperaNSide);
	PrintWarning(msg);
}

void	PixExtractProc::ProcFreq(const Zeus::HealpixHeaderAtomType& key)
{
	if(key.ValueType_ == Zeus::NO_KEY)
	{PrintWarning(HEALPIX_FREQ_WARNING);}
	else
	{
		int freq(Zeus::toInt(std::floor(key.PODValue_.doubleType_) + 0.5));
		if(freq)
			CurMapPiv_->Freq_ = freq;
	}
	PrintFreq();
}

int		PixExtractProc::ReadMaskMap(const std::wstring& dirName,const std::wstring& fname,Healpix_Map<HEALPIX_ATOM_PREC>& Mask,std::wstring& ActualFileName)
{
	if(fname.size() < 2)
		return false;

	std::wstring				tDir;
	std::wstring				ExcptStr;
	std::wstring				tFName(Zeus::ExtractFileName(fname,tDir));
	int							Context;

	if(Data2Buffer_ && tDir.empty())
	{
#ifdef WIN32
		std::wstring	MasksExt(L"Masks\\");
#else
		std::wstring	MasksExt(L"Masks/");
#endif
		tDir	= Zeus::RemoveInstanceName(DataBuffer_) + MasksExt;
		Context		= 1000;
	}
	else
	{
		Context = InOutEnvID_;
		if(tDir.empty())
			tDir = dirName;
	}

	std::auto_ptr<Zeus::GenHealpixReader<HEALPIX_ATOM_PREC> >	
		FReader(Zeus::GetGenHealpixReader(Loki::Type2Type<HEALPIX_ATOM_PREC>(),Context,tDir,tFName));

	ActualFileName = FReader->GetCollID();

	try
	{
		FReader->Initialize(0);
		FReader->Read();
		FReader->Release(Mask);
		return true;
	}
	catch(Zeus::libException& err)
	{
		ExcptStr = err.what_Xmsg();
	}
	catch(...)
	{
		ExcptStr = std::wstring(L"Non PwS exception.");
	}
	if(!(ExcptStr.empty()))
	{
		Healpix_Map<HEALPIX_ATOM_PREC> temp;
		Mask.swap(temp);
		(Zeus::ConManager::Instance())->PrintErr2Console(ExcptStr);
		(Zeus::ConManager::Instance())->PrintErr2Console(std::wstring(L"Cannot find/read file/object -> ") + FReader->GetCollID() + std::wstring(L"\n\n"));
		return false;	
	}
	return false;
}

void	PixExtractProc::OpenHealpixMap(const std::wstring& DirName,const std::wstring& MapName,Healpix_Map<HEALPIX_ATOM_PREC>& map)
{
	Zeus::HealpixHeaderAtomType	Atom;
	Zeus::HealpixHeaderType		head;
	std::wstring				tDir;
	std::wstring				tFName(Zeus::ExtractFileName(MapName,tDir));
	int							Context;

	if(Data2Buffer_ && tDir.empty())
	{
#ifdef WIN32
		std::wstring	MapsExt(L"Maps\\");
#else
		std::wstring	MapsExt(L"Maps/");
#endif
		tDir = DataBuffer_ + MapsExt;
		Context = 1000;
	}
	else
	{
		Context = InOutEnvID_;
		if(tDir.empty())  tDir = DirName;
	}
	std::auto_ptr<Zeus::GenHealpixReader<HEALPIX_ATOM_PREC> >	
		FReader(Zeus::GetGenHealpixReader(Loki::Type2Type<HEALPIX_ATOM_PREC>(),Context,
		tDir,tFName));

	std::string	DummyGcc1;
	double		DummyGcc2(0.0);
	Atom.Set(std::string("COORDSYS"),DummyGcc1);
	head.Columns_.push_back(Atom);
	Atom.Set(std::string(HEALPIX_TUNIT1),DummyGcc1);
	head.Columns_.push_back(Atom);
	Atom.Set(std::string(HEALPIX_BEAM),DummyGcc2);
	head.Columns_.push_back(Atom);
	Atom.Set(std::string(HEALPIX_FREQ),DummyGcc2);
	head.Columns_.push_back(Atom);
	FReader->Initialize(&head);
	CurMapPiv_->NSide_ = head.NSide_;
	PrintFileType(GetFileTypeName(CurMapPiv_->mapT_));
	ProcCoordSys(head.Columns_.at(0));
	ProcUnits(head.Columns_.at(1));
	ProcFreq(head.Columns_.at(3));
	FReader->Read();
	FReader->Release(map);
	if(map.Scheme() != RING)
	{
		map.swap_scheme();
	}

}
PixExtractProc::AntennasCollType::iterator
	PixExtractProc::CreateAntennas(MC_antenna::AntennaIDType Id,double freq,double FWHM)
{
	AntennasCollType::iterator			piv(AntColl_.begin());
	AntennasCollType::const_iterator	const end(AntColl_.end());

	for(;piv != end;++piv)
	{
		if((piv->Id_ == Id) && (piv->Freq_ == freq))
			return piv;
	}
	AntennaInfoType temp;
	temp.Id_		= Id;
	temp.Freq_		= freq;
	temp.Antenna_	= AntennaFactory(Id,freq,FWHM);
	AntColl_.push_back(temp);
	return AntColl_.end() - 1;
}

void	PixExtractProc::ReadPtgsFile1Coord(int coord,const std::wstring& fn)
{
	Zeus::LArr2D<double>	ws;
	std::wstring			tDir(DirOut_);

	if(Data2Buffer_)
	{
#ifdef WIN32
		std::wstring	PtgsExt(L"Ptgs\\");
#else
		std::wstring	PtgsExt(L"Ptgs/");
#endif
		if(!Zeus::CheckDir(tDir = (DataBuffer_  + PtgsExt)))
		{
			tDir = Zeus::RemoveInstanceName(DataBuffer_) + PtgsExt;
		}
	}

	std::auto_ptr<Zeus::GenCollReader<Zeus::LArr2D<double> > >
		FReader(Zeus::GetWrkSpFileReaderHandler(Loki::Type2Type<Zeus::LArr2D<double> >(),
		Data2Buffer_?1000:InOutEnvID_,fn,PtchSz_,PtchSz_,tDir));
	FReader->Initialize();
	FReader->Read();
	FReader->Release(ws);
	PutInPtgs(coord,ws);
}

void	PixExtractProc::PutInPtgs(int coord,const Zeus::LArr2D<double>&	ws)
{
	Zeus::LArr2D<double>::const_iterator			pivIn(ws.begin());
	Zeus::LArr2D<double>::const_iterator	const	pivEnd(pivIn + ws.getSz());

	Zeus::LArr2D<pointing>::iterator			pivOut(Ptgs_.begin());

	if(coord)
	{
		for(;pivIn!=pivEnd;++pivOut,++pivIn)
		{pivOut->phi = *pivIn;}
	}
	else{
		for(;pivIn!=pivEnd;++pivOut,++pivIn)
		{pivOut->theta = *pivIn;}
	}
}
//
void	PixExtractProc::PixDoRotation(PatchDirections dir,Zeus::LArr2D<double>& Pixels)
{
	const int	CentralPix(PtchSz_ >> 1);
	int			New_j,New_i;

	for(int j= -CentralPix; j<CentralPix; ++j)
	{
		for(int i= -CentralPix; i<CentralPix; ++i)
		{
			switch(dir)
			{
			case	PTCHDIR_NAT:
				New_j = j;
				New_i = i;
				break;
			case	PTCHDIR_ROT90:
				New_j =  i;
				New_i = -j;
				break;
			case	PTCHDIR_ROT180:
				New_j = -j;
				New_i = -i;
				break;
			case	PTCHDIR_ROT270:
				New_j = -i;
				New_i =  j;
				break;
			default:
				New_j = j;
				New_i = i;
				break;
			}

			if(New_j >= CentralPix ) New_j = CentralPix - 1;
			if(New_i >= CentralPix ) New_i = CentralPix - 1;

			New_j += CentralPix;
			New_i += CentralPix;

			Pixels(j + CentralPix,i + CentralPix) = Pix_(New_j,New_i);
		}
	}
}
//
void	PixExtractProc::Ptgs2Pixels(const Healpix_Map<HEALPIX_ATOM_PREC>& BadPixMask)
{
	Zeus::LArr2D<pointing>::const_iterator	piv(Ptgs_.begin());
	Zeus::LArr2D<pointing>::const_iterator	const end(piv + Ptgs_.getSz());
	Pix_.Make(PtchSz_ * PtchSz_,PtchSz_);
	Zeus::LArr2D<double>::iterator			pixPiv(Pix_.begin());

	for(;piv != end;++piv,++pixPiv)
	{
		if(BadPixMask.Map().size() && (BadPixMask[BadPixMask.ang2pix(*piv)] < 0.5))
		{
			*pixPiv = Healpix_undef;	
			continue;
		}

		*pixPiv = PixelInterpolation(*piv);
		
		if(*pixPiv > HEALPIX_BADPIX)
		{
			*pixPiv = UnitsTransf_(*pixPiv);
			if(*pixPiv > 1e+10)
			{
				*pixPiv = Healpix_undef;
			}
		}
		else
		{
			*pixPiv = Healpix_undef;
		}
	}
}

void	PixExtractProc::NormalisePix(double factor,Zeus::LArr2D<double>& DestPix)
{
	Zeus::LArr2D<double>::iterator			pixPiv(DestPix.begin());
	Zeus::LArr2D<double>::const_iterator	const pixEnd(pixPiv + DestPix.getSz());

	for(;pixPiv != pixEnd;++pixPiv)
	{*pixPiv *= factor;}
}

void 	PixExtractProc::InvertPatch(void)
{
	Zeus::LArr2D<double>::iterator			pixPiv(Pix_.begin());
	Zeus::LArr2D<double>::const_iterator	const pixEnd(pixPiv + Pix_.getSz());

	for (; pixPiv != pixEnd; ++pixPiv)
	{
		if (*pixPiv > HEALPIX_BADPIX)
		{
			*pixPiv = -(*pixPiv);
		}
	}
}

void 	PixExtractProc::ComputeStats(void)
{
	int		NPixIncluded(0);
	double	Bias(0.0);
	double	Variance(0.0);

	Zeus::LArr2D<double>::const_iterator	pixPiv(Pix_.begin());
	Zeus::LArr2D<double>::const_iterator	const pixEnd(pixPiv + Pix_.getSz());

	for(;pixPiv != pixEnd;++pixPiv)
	{
		if(*pixPiv > HEALPIX_BADPIX)
		{
			++(NPixIncluded);
			Bias		+= *pixPiv;
			Variance	+= (*pixPiv * *pixPiv);
		}
	}
	if(!NPixIncluded)
	{
		CurrStat_.CurrPatchRms_		= 0.0;
		CurrStat_.CurrPatchBias_	= 0.0;
		return;
	}
	Variance -= ((Bias * Bias) / static_cast<double>(NPixIncluded));
	double tRms(std::sqrt(Variance / static_cast<double>(NPixIncluded)));
	Bias /= static_cast<double>(NPixIncluded);
	CurrStat_.NPixIncluded_		= 0;
	CurrStat_.CurrPatchBias_	= CurrStat_.CurrPatchRms_ = 0.0;
	
	pixPiv  = Pix_.begin();
	for(;pixPiv != pixEnd;++pixPiv)
	{
		if((*pixPiv > HEALPIX_BADPIX) && (std::abs(*pixPiv - Bias) < (RmsRejectLevel_ * tRms)))
		{
			++(CurrStat_.NPixIncluded_);
			CurrStat_.CurrPatchBias_	+= *pixPiv;
			CurrStat_.CurrPatchRms_		+= (*pixPiv * *pixPiv);
		}
	}
	if(!(CurrStat_.NPixIncluded_))
	{
		CurrStat_.CurrPatchRms_		= 0.0;
		CurrStat_.CurrPatchBias_	= 0.0;
		return;
	}
	CurrStat_.CurrPatchRms_		-= ((CurrStat_.CurrPatchBias_ * CurrStat_.CurrPatchBias_) / static_cast<double>(CurrStat_.NPixIncluded_)); 
	CurrStat_.CurrPatchRms_		= std::sqrt(CurrStat_.CurrPatchRms_ / static_cast<double>(CurrStat_.NPixIncluded_));
	CurrStat_.CurrPatchBias_	/= static_cast<double>(CurrStat_.NPixIncluded_);
}


double	PixExtractProc::PixelInterpolation(const pointing& ptg0)
{
	fix_arr<int,4>			pixIndex;
	fix_arr<double,4>		wgt;
	fix_arr<double,4>		Values;
	double					TotalWgt(0.0);
	int						IncludedPix(0);

	map_.get_interpol(ptg0,pixIndex,wgt);

	Values[0] = map_[pixIndex[0]];
	Values[1] = map_[pixIndex[1]];
	Values[2] = map_[pixIndex[2]];
	Values[3] = map_[pixIndex[3]];

	for(int i=0;i<4;++i)
	{
		if((Values[i] <= HEALPIX_BADPIX) || (Values[i] >= DVDCORRUPTED_BADPIX))
		{
			wgt[i] = Values[i] = 0.0;
		}
		else
		{++IncludedPix;}
		TotalWgt += wgt[i];
	}

	if((IncludedPix < 3) || (TotalWgt < 0.1)) return Healpix_undef;

	return  (((wgt[0] * Values[0]) + (wgt[1] * Values[1]) + (wgt[2] * Values[2]) + (wgt[3] * Values[3]))/TotalWgt);
}

void	PixExtractProc::InitialiseSrcSub(Zeus::NonBlingCatType& SourceCat)
{
	//ReadGeoProps(1);
	
	// We don't need to work with the actual patch size and other geometric properties
	// when doing the subtraction of sources from the maps
	// So we adjust it to a more convenient value to do the subtraction


	TotalPatches_	= 0;
	PtchSz_			= 128;
	GeomNSide_		= 2048;
	PixSz_			= std::sqrt(PITIMES4 / (12.0 * static_cast<double>(GeomNSide_) * static_cast<double>(GeomNSide_)));

	InitFourierMachine();
}

void	PixExtractProc::MakeSrcsMaskMaps(Healpix_Map<HEALPIX_ATOM_PREC>& TempMask,Healpix_Map<HEALPIX_ATOM_PREC>& SrcsMap,Zeus::NonBlingCatType& SourceCat)
{
	Zeus::NonBlingCatType::StorageType::const_iterator	piv(SourceCat.Storage_.begin());
	Zeus::NonBlingCatType::StorageType::const_iterator	const end(SourceCat.Storage_.end());
	int													Index(0);

	//	Current implementation only handles point sources yet

	ObjFactory::ObjFactoryArgs	tObjArgs;

	tObjArgs.ObjType_	= 0;		// ID for point sources
	// !!! Important -> we are not taking into consideration possible flux calibration constants
	// If the profile needs a calibration constant, that constant must be introduced here
	// This is the consversion for point sources
	
	Zeus::PlanckUnitsValuesTranform	Uniconv(Zeus::PlanckUnitsValuesTranform::BRIGHTNESS,CurMapPiv_->MapUnits_,CurMapPiv_->Freq_,1.0/(static_cast<double>(PtchSz_ * PtchSz_) * PixSz_ * PixSz_ * CurMapPiv_->UnitFactor_ * 1.0e9));
	tObjArgs.PatchSz_	= PtchSz_;
	tObjArgs.PixSz_		= PixSz_ ;
	tObjArgs.SynID_		= syncID_;

	ObjFactory	ObjMaker(&tObjArgs);

	RealPlaneSurfType							Patch(PtchSz_);
	Zeus::FourierPlane<std::complex<double> >	Spectrum(FourierMachine_->GetXComplexBufferSz() - 1,true);
	AntennasCollType::iterator					ptrAnt(CreateAntennas(MC_antenna::ANTENNA_GAUSSIAN,CurMapPiv_->Freq_,CurMapPiv_->Ant_FWHM_));

	for(;piv != end;++piv)
	{
		if(piv->SNR_ < CurMapPiv_->Thresh_sub_)
			continue;

		// Shape creation - inline because of creation and destruction of the temporary planes
		// Not radius calibration constants yet. Calibration constant depends on the selected profile
		
		ObjInfo		OInf(ObjMaker.MakeObjectShape(piv->PredRadiusBay_));
		
		FourierMachine_->Real2Fourier(OInf.surf_.GetInnerData(),Spectrum.GetInnerData());
		Zeus::MultiplyInPlace(Spectrum,ptrAnt->Antenna_->GetAntennaBufferRef(),Zeus::UB_NOUSE);
		FourierMachine_->Fourier2Real(Spectrum.GetInnerData(),Patch.GetInnerData());
		Zeus::MultFunctor<double>	DummyGcc1(Uniconv(piv->PredFluxBay_));
		Patch.Transform(DummyGcc1,0,Zeus::UB_NOUSE);
		// end of shape creation - shape has the right units

		//char		bufferS[BUFFERMAXCHAR];	
		//sprintf(bufferS,"%s_%d_","SourceShape",static_cast<int>(Zeus::toInt(piv->SNR_)));
		//Zeus::DumpInOut_2d(bufferS,256,256,256,0,Patch.GetInnerData().begin(),1.0);

		bool mask(piv->SNR_ >= CurMapPiv_->Thresh_mask_);

		ProjectOntoSphere(mask,Patch,mask?TempMask:SrcsMap,piv->ptg_);

		if(!(Index % 80))
#ifdef WIN32
		{std::wcout << L"\n";}
		std::wcout << (mask?L".":L"*");
#else
		{std::cout << "\n";}
		std::cout << (mask?".":"*");
#endif
		++Index;
	}
#ifdef WIN32
		{std::wcout << L"\n\n";}
#else
		{std::cout << "\n\n";}
#endif
}


void	PixExtractProc::Project1Pix(double Value,const pointing& ptg,Healpix_Map<HEALPIX_ATOM_PREC>& Map)
{
	fix_arr<int,4>			pixIndex;
	fix_arr<double,4>		wgt;

	Map.get_interpol(ptg,pixIndex,wgt);

	for(int i=0;i<4;++i)
	{
		if(Map[pixIndex[i]] > HEALPIX_BADPIX)
			Map[pixIndex[i]] += (wgt[i] * Value);
	}
}

void	PixExtractProc::ProjectOntoSphere(bool mask,RealPlaneSurfType& Patch,Healpix_Map<HEALPIX_ATOM_PREC>& Map,const pointing& posIn)
{
	rotmatrix	rot0;
	pointing	ptg(posIn.theta - PIOVER2,posIn.phi);
	rot0.Make_CPAC_Euler_Matrix(ptg.phi,ptg.theta,0.0);
	const int Limit(PtchSz_ >> 1);
	double	MaskRadius;
	const	double	Sigma((static_cast<double>(EvalMaskRadiusPix(Patch)) * PixSz_ * 2.0) / FWHM2SIGMA);
	
	if(MaskRadiusRatio_ > 0.0)
	{
		MaskRadius = (Sigma * MaskRadiusRatio_);
	}
	else
	{
		MaskRadius = (Sigma * std::sqrt(2.0 * std::log(CurMapPiv_->Thresh_mask_ * MASKINGRMSFACTOR)));		
	}
		
	MaskRadius	*= MaskRadius;
	for(int i= -Limit;i <= Limit; ++i)
	{
		const double ZDisp(i*PixSz_);
		for(int j= -Limit;j < Limit; ++j)
		{
			const double	YDisp(j*PixSz_);
			pointing		tptg(rot0.Transform(vec3(1.0,YDisp,ZDisp)));
			double			Value(Patch(i,j));
			if(mask)
			{
				const double	PixDistSq(YDisp*YDisp+ZDisp*ZDisp);
				if(MaskRadius >= PixDistSq)
				{Value = -1.0;}
				else
				{
					EvalValueForMasking(Value,tptg);
				}
			}
			if(Value != 0.0)
			{
				Project1Pix(Value,tptg,Map);
			}
		}
	}
}

void	PixExtractProc::EvalValueForMasking(double& Value,const pointing& tptg)
{
	std::vector<int>	pixColl;
	
	map_.query_disc_inclusive(tptg,PixSz_ * static_cast<double>(3.0),pixColl);
	std::vector<int>::const_iterator	pixPiv(pixColl.begin());
	std::vector<int>::const_iterator	const pixEnd(pixColl.end());
	int									NPix(0);
	double								SqAcum(0.0);
	double								Acum(0.0),temp,LRms;

	for(;pixPiv != pixEnd;++pixPiv)
	{
		if((temp = map_[*pixPiv]) <= HEALPIX_BADPIX)
			continue;
		++NPix;
		SqAcum	+= (temp*temp);
		Acum	+=	temp;
	}

	if(NPix)
	{
		SqAcum -= ((Acum * Acum) / static_cast<double>(NPix));
		if(SqAcum <= 0.0)
		{LRms = 0.0;}
		else
		{LRms = (std::sqrt(SqAcum / static_cast<double>(NPix)));}
	}
	else
	{LRms = 1.0e32;}

	Value = (((Value * MASKINGRMSFACTOR)  < LRms)?0.0:-1.0);

	return;
}

void	PixExtractProc::SubSourcesDoProcUpAllMaps(void)
{
	CurMapPiv_		= FileCollection_.begin();
	PixExtractProcCtorArgsType::MapFileInfoCollType::const_iterator const end(FileCollection_.end());
	Zeus::NonBlingCatType			SourceCat;
	Healpix_Map<HEALPIX_ATOM_PREC>	TempMask(MC_MAXNSIDE,RING,nside_dummy());
	Healpix_Map<HEALPIX_ATOM_PREC>	SrcsMap(MC_MAXNSIDE,RING,nside_dummy());

	InitialiseSrcSub(SourceCat);

	{
		// Deallocates the memory taken by map
		Healpix_Map<HEALPIX_ATOM_PREC>	Dummy(MC_MAXNSIDE,RING,nside_dummy());
		IllPixMask_.swap(Dummy);
		IllPixMask_.fill(1.0);
	}

	for(;CurMapPiv_ != end;++CurMapPiv_)
	{
		PrintProcessingFile();
		if(CurMapPiv_->Processed_)
		{
			(Zeus::ConManager::Instance())->PrintStr2Console(std::wstring(ALREADYPROC_FILE) + std::wstring(L"\n\n\n"));
			continue;
		}

		ReadSrcCatalogue(CurMapPiv_->PS_Catal_FName_,SourceCat);

		if(SourceCat.Header_.CoordsType_)
		{NormaliseCoords(SourceCat);}
		if(SourceCat.Header_.CoordSystem_ != CurMapPiv_->MapCoordSys_)
		{TranslateCoords(SourceCat,static_cast<coordsys>(CurMapPiv_->MapCoordSys_));}

		OpenHealpixMap(HP_Dir_,CurMapPiv_->HP_FName_,map_);

		TempMask.fill(1.0);
		SrcsMap.fill(0.0);

		// Make the maps with the sources and masks
		MakeSrcsMaskMaps(TempMask,SrcsMap,SourceCat);
		ClampLevelMask(TempMask);		

		std::wstring				tDummy;
		std::wstring				tFName(Zeus::ExtractFileName(CurMapPiv_->HP_FName_,tDummy));

		Zeus::CreatReadableHealpixFile(InOutEnvID_,Masks_Dir_,std::wstring(SRCMAPPREFIX) + tFName,Galactic,SrcsMap);
#if (defined(HFIDMC) || defined(LFIDPC)) && defined(HFIDMC_EXTOBJECTS)
		if(Data2Buffer_)
		{
			Zeus::CreatReadableExternalHealpixFile(DataBuffer_,std::wstring(SRCMAPPREFIX) + tFName,Galactic,SrcsMap);
		}
#endif

		// Write bad pix mask for this frequency

		Zeus::CreatReadableHealpixFile(InOutEnvID_,Masks_Dir_,std::wstring(MASKMAPPREFIX) + tFName,Galactic,TempMask);
#if (defined(HFIDMC) || defined(LFIDPC)) && defined(HFIDMC_EXTOBJECTS)
		if(Data2Buffer_)
		{
			Zeus::CreatReadableExternalHealpixFile(DataBuffer_,std::wstring(MASKMAPPREFIX) + tFName,Galactic,TempMask);
		}
#endif
		// Subtract map src from original map

		// Write new map (with srcs subtracted) for this frequency

		Zeus::CreatReadableHealpixFile(InOutEnvID_,Masks_Dir_,std::wstring(NEWMAPPREFIX) + tFName,Galactic,SubAddMapsInPlace(map_,SrcsMap,0));
#if (defined(HFIDMC) || defined(LFIDPC)) && defined(HFIDMC_EXTOBJECTS)
		if(Data2Buffer_)
		{
			Zeus::CreatReadableExternalHealpixFile(DataBuffer_,std::wstring(NEWMAPPREFIX) + tFName,Galactic,map_);
		}
#endif
		// Merge mask map with temporary mask

		MergeMapsInPlace(TempMask);

		(Zeus::ConManager::Instance())->PrintStr2Console(L"\n\n\n");
	}

	{
		// Deallocates the memory taken by map
		Healpix_Map<HEALPIX_ATOM_PREC>	Dummy;
		Dummy.swap(map_);
	}

	ReadBadPixMask(Masks_Dir_);

	// Merge the badpix mask with the temporary
	// Write new bad pix mask

	std::wstring				tDummy;
	std::wstring				tFName(Zeus::ExtractFileName(MaskRemoveName_,tDummy));

	Zeus::CreatReadableHealpixFile(InOutEnvID_,Masks_Dir_,std::wstring(NEWMAPPREFIX) + tFName,Galactic,MergeMapsInPlace(CutOutMask_));
#if (defined(HFIDMC) || defined(LFIDPC)) && defined(HFIDMC_EXTOBJECTS)
	if(Data2Buffer_)
	{
		Zeus::CreatReadableExternalHealpixFile(DataBuffer_,std::wstring(NEWMAPPREFIX) + tFName,Galactic,IllPixMask_);
	}
#endif
}
//
void	PixExtractProc::NormaliseCoords(Zeus::NonBlingCatType& SrsCat)
{
	Zeus::NonBlingCatType::StorageType::iterator			piv(SrsCat.Storage_.begin());
	Zeus::NonBlingCatType::StorageType::const_iterator	const end(SrsCat.Storage_.end());

	for(;piv != end;++piv)
	{
		if((piv->ptg_.phi   < -1.0e10) || (piv->ptg_.theta < -1.0e10))
			continue;
		piv->ptg_.theta		= ((90.0 - piv->ptg_.theta) / RAD2DEGREE); 
		piv->ptg_.phi		/= RAD2DEGREE;	
	}
}
//
void	PixExtractProc::TranslateCoords(Zeus::NonBlingCatType& SrsCat,const coordsys PtgsCoordSystem)
{
	Trafo	CoordTranslator(SrsCat.Header_.Epoch_,MAPEPOCH,SrsCat.Header_.CoordSystem_,PtgsCoordSystem);
	Zeus::NonBlingCatType::StorageType::iterator			piv(SrsCat.Storage_.begin());
	Zeus::NonBlingCatType::StorageType::const_iterator	const end(SrsCat.Storage_.end());

	for(;piv != end;++piv)
	{
		if((piv->ptg_.phi   < -1.0e10) || (piv->ptg_.theta < -1.0e10))
			continue;
		piv->ptg_ = CoordTranslator(piv->ptg_);
	}
}
//
Healpix_Map<HEALPIX_ATOM_PREC>&	PixExtractProc::SubAddMapsInPlace(Healpix_Map<HEALPIX_ATOM_PREC>& Accumul,Healpix_Map<HEALPIX_ATOM_PREC>& SrcsMap,int AddSub)
{

	if(Accumul.Nside() != SrcsMap.Nside())
	{
		Healpix_Map<HEALPIX_ATOM_PREC>	temp(Accumul.Nside(),RING,nside_dummy());
		temp.Import_degrade(SrcsMap,true);
		SrcsMap.swap(temp);
	}

	const int	NPix(Accumul.Npix());

	for(int i=0;i<NPix;++i)
	{
		if(Accumul[i]  > HEALPIX_BADPIX)
		{
			if(AddSub) {Accumul[i] += SrcsMap[i];}
			else {Accumul[i] -= SrcsMap[i];}
		}
	}
	return Accumul;
}

void	PixExtractProc::DoPixExtractionAllMaps(int OpMode,long FPatch,long LPatch)
{
	OpMode_ = OpMode;
	CurMapPiv_		= FileCollection_.begin();
	PixExtractProcCtorArgsType::MapFileInfoCollType::const_iterator const end(FileCollection_.end());
	Initialize();

	if(FPatch < 0) FPatch = 0;
	if((FPatch > LPatch) && (LPatch >=0))
		goto _go_out;

	for(CurrPatchNumber_ = FPatch;CurMapPiv_ != end;++CurMapPiv_,CurrPatchNumber_ = FPatch)
	{
		PrintProcessingFile();
		if(CurMapPiv_->Processed_)
		{
			(Zeus::ConManager::Instance())->PrintStr2Console(std::wstring(ALREADYPROC_FILE) + std::wstring(L"\n\n\n"));
			continue;
		}

		OpenHealpixMap(HP_Dir_,CurMapPiv_->HP_FName_,map_);

		MaskMap();

		if(!OpMode_)
		{
			PreProcessMap();
			if(map_.Scheme() != NEST)
			{
				map_.swap_scheme();
			}
			CreateDir();
			DoPixExtraction1Map(LPatch);
		}
		(Zeus::ConManager::Instance())->PrintStr2Console(L"\n\n\n");
	}
_go_out:
	if(OpMode_)
	{
		WriteMasks();
	}
#if defined(HFIDMC) || defined(LFIDPC)
	else
	{
		wchar_t			buffer[DOUBLETXTMAXSZ];
		std::wstring	tName(MC_PIPELINESYNC_DEF);
		swprintf(buffer,DOUBLETXTMAXSZ,L"%d",syncID_);
		tName += std::wstring(buffer);

		Healpix_Map<HEALPIX_ATOM_PREC>	temp(2,RING,nside_dummy());
		Zeus::CreatReadableHealpixFile(InOutEnvID_,Masks_Dir_,
			tName,Galactic,temp);
	}
#endif
}

std::wstring	PixExtractProc::GetFileTypeDir(void)
{
	std::wstring temp;
	switch(CurMapPiv_->mapT_)
	{
	case Zeus::MAPTYPE_OBSERVATION:
		temp = std::wstring(MT_OBSERVATION_DIR);
		break;
	case Zeus::MAPTYPE_BACKGROUND:
		temp = std::wstring(MT_BACKCMB_DIR);
		break;
	default:
		return std::wstring(L"<INVALID>");
	}
	return temp;
}


std::wstring	PixExtractProc::GetFileTypePrefix(void) const
{
	std::wstring temp;

	switch(CurMapPiv_->mapT_)
	{
	case Zeus::MAPTYPE_OBSERVATION:
		temp = std::wstring(MT_OBSERVATION_FILE);
		break;
	case Zeus::MAPTYPE_BACKGROUND:
		temp += std::wstring(MT_BACKCMB_FILE);
		break;
	default:
		return std::wstring(L"<INVALID>");
	}
	return temp;
}


void	PixExtractProc::CreateDir(void)
{
	bool result(true);
	DirOutCurrMap_.clear();
	std::wstring	suffix(GetFileTypeDir());

	std::wstring temp(Zeus::CorrectDir(DirOutMaps_ + suffix));

	if(Data2Buffer_)
	{
		result = Zeus::CreateDir(temp = Zeus::CorrectDir(DataBuffer_ + suffix,1));
	}
	else
	{
#if	!defined(HFIDMC) && !defined(LFIDPC)
		result = Zeus::CreateDir(temp);
#endif
	}

	if(result)
	{DirOutCurrMap_ = temp;}
	else
	{PrintDirErr(temp);}
}

void	PixExtractProc::PreProcessMap(void)
{
	PlanckMapsInfoCollType::const_iterator	tpiv;
	if((tpiv = GetPlanckData(std::abs(CurMapPiv_->FreqID_)%1000)) == PlanckEnd_)
		errHealpix(ERRCOD_MC_INVALIDFREQ,ERRMSG_MC_INVALIDFREQ,HP_Dir_,CurMapPiv_->HP_FName_,1);

	if(CurMapPiv_->NSide_ != tpiv->OperatNSide_)
	{
		if(CurMapPiv_->NSide_ < tpiv->OperatNSide_)
		{
			PrintObsMapLowRes(CurMapPiv_->NSide_,tpiv->OperatNSide_);
		}
		else
		{
			//const nside_dummy dummy;
			Healpix_Map<HEALPIX_ATOM_PREC>	temp(tpiv->OperatNSide_,RING,nside_dummy());
			temp.Import_degrade(map_,true);
			map_.swap(temp);
			CurMapPiv_->NSide_ = tpiv->OperatNSide_;
		}
	}
}

//
// ************************************
// ** Main entring point for masking **
// ************************************
//
void	PixExtractProc::MaskMap(void)
{
	Healpix_Map<HEALPIX_ATOM_PREC>	OutOrg;
	Healpix_Map<HEALPIX_ATOM_PREC>	IllPixTemp;

	if(!OpMode_ || !CheckBadPixels()) return;

	if(!(IllPixMask_.Map().size()))
	{
		Healpix_Map<HEALPIX_ATOM_PREC>	temp(MC_MAXNSIDE,RING,nside_dummy());
		temp.fill(1.0);
		IllPixMask_.swap(temp);
	}

	{
		Healpix_Map<HEALPIX_ATOM_PREC>	temp(CurMapPiv_->NSide_,RING,nside_dummy());
		temp.fill(1.0);
		IllPixTemp.swap(temp);
		MarkBadPixels(IllPixTemp);
	}

	if(IllPixTemp.Nside() != MC_MAXNSIDE)
	{
		Healpix_Map<HEALPIX_ATOM_PREC>	temp(MC_MAXNSIDE,RING,nside_dummy());
		temp.Import_upgrade(IllPixTemp);
		IllPixTemp.swap(temp);
	}

	const int	NPix(IllPixMask_.Npix());

	for(int i=0;i<NPix;++i)
	{
		if(IllPixTemp[i] < 0.5)
		{IllPixMask_[i] = 0.0;}
	}
}
//
void	PixExtractProc::WriteMasks(void)
{
	Healpix_Map<HEALPIX_ATOM_PREC>	RejectionMask;

	if(CutOutMask_.Map().size())
	{
		if(CutOutMask_.Nside() != MC_MAXNSIDE)
		{
			Healpix_Map<HEALPIX_ATOM_PREC>	temp(MC_MAXNSIDE,RING,nside_dummy());
			temp.Import_upgrade(CutOutMask_);
			CutOutMask_.swap(temp);
		}

		if(!(IllPixMask_.Map().size()))
		{
			Healpix_Map<HEALPIX_ATOM_PREC>	temp(MC_MAXNSIDE,RING,nside_dummy());
			temp.fill(1.0);
			IllPixMask_.swap(temp);
		}

		const int	tNPix(CutOutMask_.Npix());

		for(int i=0;i<tNPix;++i)
		{
			if(CutOutMask_[i] < 0.5)
			{
				IllPixMask_[i] = 0.0;
			}
		}	

	}

	if((MaskRejectName_.size() >= 2) || IllPixMask_.Map().size())
	{

		if(MaskRejectName_.size() < 2)
		{
			Healpix_Map<HEALPIX_ATOM_PREC>	temp(MC_MAXNSIDE,RING,nside_dummy());
			temp.fill(1.0);
			RejectionMask.swap(temp);
		}
		else
		{
			Healpix_Map<HEALPIX_ATOM_PREC>	temp;

			(Zeus::ConManager::Instance())->PrintStr2Console(std::wstring(MC_MASKREADINGREJECT) + std::wstring(MC_MASKMASKALERT) + std::wstring(L"\n\n"));

			std::wstring	ActualFileName;
			int fileExits(ReadMaskMap(Masks_DirIn_,MaskRejectName_,temp,ActualFileName));
		
			if((MaskRejectName_.size() >= 2) && !fileExits)
			{
				errObjectNotFound(ERRCOD_MC_MAPNOTFOUND,std::wstring(L""),ERRMSG_MC_MAPNOTFOUND,ActualFileName,1);
			}
 
			if(temp.Map().size())
			{
				if(temp.Scheme() != RING)
				{temp.swap_scheme();}

				if(temp.Nside() != MC_MAXNSIDE)
				{
					Healpix_Map<HEALPIX_ATOM_PREC>	temp1(MC_MAXNSIDE,RING,nside_dummy());
					temp1.Import_upgrade(temp);
					RejectionMask.swap(temp1);
				}
				else
				{
					RejectionMask.swap(temp);
				}
			}
			else
			{
				Healpix_Map<HEALPIX_ATOM_PREC>	temp1(MC_MAXNSIDE,RING,nside_dummy());
				temp1.fill(1.0);
				RejectionMask.swap(temp1);
			}
		}
	}

	if(IllPixMask_.Map().size())
	{
		const int	tNPix(IllPixMask_.Npix());
		Healpix_Map<HEALPIX_ATOM_PREC>	OutMask(MC_MAXNSIDE,RING,nside_dummy());

		// Write the smooth illpix mask to the disk

		(Zeus::ConManager::Instance())->PrintStr2Console(MC_MASKFILTERMASK1);
		float MapBias;

		// Enlarge Badpix mask by 2 pixels
		OutMask.fill(1.0);
		for(int i=0;i<tNPix;++i)
		{	
			if(IllPixMask_[i] < 0.5)
			{SetBadPix(i,OutMask);}
		}

		IllPixMask_.swap(OutMask);

		// Write the illpix mask to the disk IllPixMask_
		
		Zeus::CreatReadableHealpixFile(InOutEnvID_,Masks_Dir_,
			std::wstring(MC_MASKBADPIX_DEF),Galactic,IllPixMask_);

#if (defined(HFIDMC) || defined(LFIDPC)) && defined(HFIDMC_EXTOBJECTS)
		if(Data2Buffer_)
		{
			Zeus::CreatReadableExternalHealpixFile(DataBuffer_,std::wstring(MC_MASKBADPIX_DEF),Galactic,IllPixMask_);
		}
#endif

		(Zeus::ConManager::Instance())->PrintStr2Console(MC_MASKFILTERMASK2);
		MapBias = static_cast<float>(IllPixMask_.average());
		FilterMap(IllPixMask_,OutMask,MaskSmoothFWHM_);
		OutMask.Add(MapBias);

		for(int i=0;i<tNPix;++i)
		{
#ifdef MAKEREDUCEDMASKS
			OutMask[i] = ((OutMask[i] < 0.03)? 0.0 : 1.0);
#else
			OutMask[i] = (((OutMask[i] < MaskEnlagThreshold_) || (IllPixMask_[i] < 0.5))? 0.0 : 1.0);
#endif
		}

		for(int i=0;i<tNPix;++i)
		{
			if(OutMask[i] < 0.5)
			{
				RejectionMask[i] = 0.0;
			}
		}	
	}
	else
	{
		Healpix_Map<HEALPIX_ATOM_PREC>	Dummy(MC_MAXNSIDE,RING,nside_dummy());

		Dummy.fill(1.0);

		Zeus::CreatReadableHealpixFile(InOutEnvID_,Masks_Dir_,
			std::wstring(MC_MASKBADPIX_DEF),Galactic,Dummy);

#if (defined(HFIDMC) || defined(LFIDPC)) && defined(HFIDMC_EXTOBJECTS)
		if(Data2Buffer_)
		{
			Zeus::CreatReadableExternalHealpixFile(DataBuffer_,std::wstring(MC_MASKBADPIX_DEF),Galactic,Dummy);
		}
#endif
	}

	if(RejectionMask.Map().size())
	{
		Zeus::CreatReadableHealpixFile(InOutEnvID_,Masks_Dir_,
			std::wstring(MC_MASKREJNAME_DEF),Galactic,RejectionMask);

#if (defined(HFIDMC) || defined(LFIDPC)) && defined(HFIDMC_EXTOBJECTS)
		if(Data2Buffer_)
		{
			Zeus::CreatReadableExternalHealpixFile(DataBuffer_,std::wstring(MC_MASKREJNAME_DEF),Galactic,RejectionMask);
		}
#endif
	}
	else
	{
		Healpix_Map<HEALPIX_ATOM_PREC>	Dummy(MC_MAXNSIDE,RING,nside_dummy());

		Dummy.fill(1.0);

		Zeus::CreatReadableHealpixFile(InOutEnvID_,Masks_Dir_,
			std::wstring(MC_MASKREJNAME_DEF),Galactic,Dummy);
#if (defined(HFIDMC) || defined(LFIDPC)) && defined(HFIDMC_EXTOBJECTS)
		if(Data2Buffer_)
		{
			Zeus::CreatReadableExternalHealpixFile(DataBuffer_,std::wstring(MC_MASKREJNAME_DEF),Galactic,Dummy);
		}
#endif
	}
}
//
void		PixExtractProc::SetBadPix(int Pix,Healpix_Map<HEALPIX_ATOM_PREC>& Map)
{
	std::vector<int>	pixColl;
	
	IllPixMask_.query_disc(IllPixMask_.pix2ang(Pix),PixSz_ * 2.0,pixColl);
	std::vector<int>::const_iterator	pixPiv(pixColl.begin());
	std::vector<int>::const_iterator	const pixEnd(pixColl.end());

	for(;pixPiv != pixEnd;++pixPiv)
	{Map[*pixPiv] = 0.0;}
}
//
int		PixExtractProc::CheckBadPixels(void)
{
	const int	NPix(map_.Npix());
	float		Dummy;

	for(int i=0;i<NPix;++i)
	{
		if((Dummy = map_[i]) <= HEALPIX_BADPIX)
			return 1;
	}

	return 0;
}
//
void	PixExtractProc::MarkBadPixels(Healpix_Map<HEALPIX_ATOM_PREC>& Mask)
{
	const int	NPix(Mask.Npix());

	for(int i=0;i<NPix;++i)
	{
		if(map_[i] <= HEALPIX_BADPIX)
			Mask[i] = 0.0;
	}
}
//
void	PixExtractProc::CombineMaps(const Healpix_Map<HEALPIX_ATOM_PREC>& OutOrg, const Healpix_Map<HEALPIX_ATOM_PREC>& Mask,Healpix_Map<HEALPIX_ATOM_PREC>& Result)
{
	const int	NPix(Mask.Npix());
	float		weight;

	for(int i=0;i<NPix;++i)
	{

		if(Mask[i] > 0.0)
		{
			if(Mask[i] < 1.0)
			{
				weight = Mask[i];
				Result[i] = (Result[i] * weight + OutOrg[i] * (1.0 - weight));
			}
		}
		else
		{
			Result[i] = OutOrg[i];
		}
	}
}
//
void	PixExtractProc::FilterMap(const Healpix_Map<HEALPIX_ATOM_PREC>& InMap,Healpix_Map<HEALPIX_ATOM_PREC>& OutMap,double fwhm)
{
	{
		Healpix_Map<HEALPIX_ATOM_PREC>	temp(InMap);
		OutMap.swap(temp);
	}

	int	maxL(Zeus::toInt((((5.0 * FWHM2SIGMA*RAD2ARCMIN) / fwhm)) + 1.0));

	Alm<xcomplex<float> > alm(maxL,maxL);

	double avg(OutMap.average());
	OutMap.Add(-avg);
	if (OutMap.Scheme()==NEST)
	{
		OutMap.swap_scheme();
	}

	(Zeus::ConManager::Instance())->PrintStr2Console(MC_MASKFORWARDTRANS);

	map2alm_iter(OutMap,alm,0);

	(Zeus::ConManager::Instance())->PrintStr2Console(MC_MASKFILTERING);

	smoothWithGauss(alm, (fwhm / RAD2ARCMIN));

	(Zeus::ConManager::Instance())->PrintStr2Console(MC_MASKBACKTRANS);

	alm2map(alm,OutMap);
}
//
MC_antenna*	PixExtractProc::AntennaFactory(MC_antenna::AntennaIDType Id,double Freq, double FWHM)
{
	MC_antenna* ant(0);

	switch(Id)
	{
	case MC_antenna::ANTENNA_GAUSSIAN:
		ant = new AntennaGaussian(std::sqrt(PITIMES4 / (GeomNSide_ * GeomNSide_ * 12.0)) * RAD2ARCMIN,
			PtchSz_ >> 1,-1.0,Freq,FWHM,FWHM);
		break;
	case MC_antenna::ANTENNA_WINDOW:
		ant = new AntennaPixWindow(std::sqrt(PITIMES4 / (GeomNSide_ * GeomNSide_ * 12.0)) * RAD2ARCMIN, PtchSz_ >> 1,FWHM);
		break;
	}
	ant->Initialise();
	return ant;
}
//
void	PixExtractProc::ReportPixCurrentParams(void) const
{
	wchar_t	buffer[BUFFERMAXCHAR];

	PRINTINTOBUFFERFUNCT (buffer,BUFFERMAXCHAR,MC_PIXELSFORMSTR,
		CurrPatchNumber_,
		PatchGeom_.Storage_[CurrPatchNumber_].InitRejectRatio_  * 100.0,
		CurrStat_.PercentBadpixAfterRot_ * 100.0,
		CurrStat_.NPixIncluded_,
		CurrStat_.CurrPatchBias_,
		CurrStat_.CurrPatchRms_
	);

	(Zeus::ConManager::Instance())->PrintStr2Console(std::wstring(buffer));
}

void	PixExtractProc::PrintBeamSigma(double orgFWHM,double newFWHM) const
{
	std::wstring msg;
	if(orgFWHM > 0.0)
	{
		msg = HEALPIX_ORGBEAM_MSG;
		msg += Zeus::PutNumber2Txt(orgFWHM);
		msg += std::wstring(L"\n");
		(Zeus::ConManager::Instance())->PrintStr2Console(msg);
		if(newFWHM <= 0.0)
		{
			(Zeus::ConManager::Instance())->PrintStr2Console(L"\n");
			return;
		}
		if(newFWHM != orgFWHM)
		{msg = std::wstring(HEALPIX_BEAM_CORRECTED_MSG);}
	}
	else
	{
		if(newFWHM <= 0.0)
		{
			(Zeus::ConManager::Instance())->PrintStr2Console(L"\n");
			return;
		}
		msg = std::wstring(HEALPIX_BEAM_MSG);
	}
	
	msg += Zeus::PutNumber2Txt(newFWHM);
	msg += std::wstring(L"\n\n");
	(Zeus::ConManager::Instance())->PrintStr2Console(msg);
}

void	PixExtractProc::Extract1stMap(const std::wstring& fName)
{
	std::string	cmdprefixR("rm -R ");

#ifdef REMOVEFULLTEMPDIR
	system((cmdprefixR + Zeus::Wstr2Str(fName)).c_str());
#else
	std::string	cmdprefix("rm ");
	std::string	MainPath(Zeus::Wstr2Str(Zeus::CorrectDir(fName,1)));

#ifdef WIN32
		std::string	NonCatsFiles("Cats\\PwSCPAC*");
		std::string	TempMasksFiles("Masks\\PwS_*");
#else
		std::string	NonCatsFiles("Cats/PwSCPAC*");
		std::string	TempMasksFiles("Masks/PwS_*");
#endif

	system((cmdprefix + MainPath + NonCatsFiles).c_str());
	system((cmdprefix + MainPath + TempMasksFiles).c_str());
	system((cmdprefixR + MainPath + std::string("Geom")).c_str());
	system((cmdprefixR + MainPath + std::string("Ptgs")).c_str());
	system((cmdprefixR + MainPath + std::string("_FreqMaps")).c_str());
	system((cmdprefixR + MainPath + std::string("_BackCMB")).c_str());

#endif // REMOVEFULLTEMPDIR 
}

void	PixExtractProc::JoinManyMasks(void)
{
	CurMapPiv_		= FileCollection_.begin();
	PixExtractProcCtorArgsType::MapFileInfoCollType::const_iterator const end(FileCollection_.end());

	if(CurMapPiv_ == end)
		return;

	Healpix_Map<HEALPIX_ATOM_PREC>	Acumulator(MC_MAXNSIDE,RING,nside_dummy());
	Acumulator.fill(1.0);

	for(;CurMapPiv_ != end;++CurMapPiv_)
	{
		std::wstring	dummy(PROCESSING_FILE);
		dummy	+= CurMapPiv_->HP_FName_;
		dummy	+= std::wstring(L"\n");
		(Zeus::ConManager::Instance())->PrintStr2Console(dummy);
		if(CurMapPiv_->Processed_)
		{
			(Zeus::ConManager::Instance())->PrintStr2Console(std::wstring(ALREADYPROC_FILE) + std::wstring(L"\n"));
			continue;
		}
		{
			Healpix_Map<HEALPIX_ATOM_PREC>		map;
			if(CurMapPiv_->HP_FName_.size() < 2)
				continue;
			OpenHealpixMap(Masks_Dir_,CurMapPiv_->HP_FName_,map);

			if(map.Scheme() != RING)
			{map.swap_scheme();}

			if(map.Nside() != MC_MAXNSIDE)
			{
				Healpix_Map<HEALPIX_ATOM_PREC>	temp1(MC_MAXNSIDE,RING,nside_dummy());
				temp1.Import_upgrade(map);
				map.swap(temp1);
			}
			
			const int	tNPix(map.Npix());

			for(int i=0;i<tNPix;++i)
			{
				if(map[i] < 0.5)
				{
					Acumulator[i] = 0.0;
				}
			}	
		}
	}

	std::wstring			tDir;
	std::wstring			tFName(Zeus::ExtractFileName(MaskRejectName_,tDir));
	int						Context;

	if(Data2Buffer_ && tDir.empty())
	{
#ifdef WIN32
		std::wstring	MapsExt(L"Maps\\");
#else
		std::wstring	MapsExt(L"Maps/");
#endif
		tDir = DataBuffer_ + MapsExt;
		Zeus::CreateDir(tDir);
		Context = 1000;
	}
	else
	{
		Context = InOutEnvID_;
		if(tDir.empty())  tDir = Masks_Dir_;
	}

	Zeus::CreatReadableHealpixFile(Context,tDir,tFName,Galactic,Acumulator);
}
//
void	PixExtractProc::ReadPtgsFileIn(int PatchNumber, Zeus::LArr2D<double>& theta,Zeus::LArr2D<double>& phi)
{
	std::wstring			tDir(DirOut_);

	if(Data2Buffer_)
	{
#ifdef WIN32
		std::wstring	PtgsExt(L"Ptgs\\");
#else
		std::wstring	PtgsExt(L"Ptgs/");
#endif
		if(!Zeus::CheckDir(tDir = (DataBuffer_  + PtgsExt)))
		{
			tDir = Zeus::RemoveInstanceName(DataBuffer_) + PtgsExt;
		}
	}

	std::auto_ptr<Zeus::GenCollReader<Zeus::LArr2D<double> > >
		FReader(Zeus::GetWrkSpFileReaderHandler(Loki::Type2Type<Zeus::LArr2D<double> >(),Data2Buffer_?1000:InOutEnvID_,GetCurrPatchPtgsOutputName(PatchNumber,0),
		PtchSz_,PtchSz_,tDir));

	FReader->Initialize();
	FReader->Read();
	FReader->Release(theta);

	std::auto_ptr<Zeus::GenCollReader<Zeus::LArr2D<double> > >
		FReader1(Zeus::GetWrkSpFileReaderHandler(Loki::Type2Type<Zeus::LArr2D<double> >(),Data2Buffer_?1000:InOutEnvID_,GetCurrPatchPtgsOutputName(PatchNumber,1),
		PtchSz_,PtchSz_,tDir));

	FReader1->Initialize();
	FReader1->Read();
	FReader1->Release(phi);
}
//
void	PixExtractProc::CreatHealpixMapWithPtgs(const std::wstring& fName)
{
	Zeus::LArr2D<double>	theta;
	Zeus::LArr2D<double>	phi;
	Healpix_Map<HEALPIX_ATOM_PREC>	HPixMap(MC_MAXNSIDE,RING,nside_dummy());
	HPixMap.fill(0.0);
	int			PatchesProcessed(0);

	ReadGeoProps(0);

	const int	TotNPatches(PatchGeom_.Storage_.size());	

	for(int PatchNumber=0;PatchNumber != TotNPatches;++PatchNumber)
	{
		try
		{
			if(!((PatchGeom_.Storage_.at(PatchNumber)).PatchValid_))
			{
				ReadPtgsFileIn(PatchNumber,theta,phi);
			}
			else continue;
		}
		catch(...)
		{
			break;
		}
		for(int j=BorderSz_;j<(PtchSz_ - BorderSz_);++j)
		{
			for(int i=BorderSz_;i<(PtchSz_ - BorderSz_);++i)
			{
				HPixMap[HPixMap.ang2pix(pointing(theta(j,i),phi(j,i)))] += 1.0;
			}
		}
		if(!(PatchesProcessed % 80))
#ifdef WIN32
		{std::wcout << L"\n";}
		std::wcout << L".";
#else
		{std::cout << "\n";}
		std::cout << ".";
#endif
		++PatchesProcessed;
	}


	if(PatchesProcessed)
	{

		std::wstring			tDir;
		std::wstring			tFName(Zeus::ExtractFileName(fName,tDir));
		int						Context;

		if(Data2Buffer_ && tDir.empty())
		{
	#ifdef WIN32
			std::wstring	MapsExt(L"Maps\\");
	#else
			std::wstring	MapsExt(L"Maps/");
	#endif
			tDir = DataBuffer_ + MapsExt;
			Zeus::CreateDir(tDir);
			Context = 1000;
		}
		else
		{
			Context = InOutEnvID_;
			if(tDir.empty())  tDir = Masks_Dir_;
		}

		Zeus::CreatReadableHealpixFile(Context,tDir,tFName,Galactic,HPixMap);
	}
	else
	{
		(Zeus::ConManager::Instance())->PrintStr2Console(std::wstring(L"Error reading the pointings or no patches to paint \n\n"));
	}
}
//

