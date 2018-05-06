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
#include <memory>
#include "ZEUS_InOut.h"
#include "PWS_BackgroundProcessor.h"
#include "PWS_PatchProcessor.h"
#include "PWS_OddsEval.h"
#include "PWS_BrentMinimizer.h"
#include "PWS_PowellSnakes.h"
#include "PWS_NestSamplCases.h"
#include "ZEUS_Debug.h"
#include "qa_assess_contours.h"

//----------------------------------
//
static const int		YPLOTSZ(256);
#ifdef JUSTCENTRALPIXEL
static const int		POS_LIM(1);
#else
static const int		POS_LIM(4);
#endif
static const double		NSIGMA(6.0);
static const int		NOFPRIORS(4);
static const double		RADIUSTOL(1.0e-8);

std::wstring	PatchProcessor::GetFileTypePrefix(Zeus::MapType mapT) const
{
	std::wstring temp;

	switch(mapT)
	{
	case Zeus::MAPTYPE_OBSERVATION:
		temp = std::wstring(PATCH_OBSERVATION_FILE);
		break;
	case Zeus::MAPTYPE_BACKGROUND:
		temp += std::wstring(PATCH_BACKCMB_FILE);
		break;
	default:
		return std::wstring(L"<INVALID>");
	}
	return temp;
}
//
std::wstring	PatchProcessor::GetFileTypeDir(Zeus::MapType mapT) const
{
	std::wstring temp;
	switch(mapT)
	{
	case Zeus::MAPTYPE_OBSERVATION:
		temp = std::wstring(PATCH_OBSERVATION_DIR);
		break;
	case Zeus::MAPTYPE_BACKGROUND:
		temp = std::wstring(PATCH_BACKCMB_DIR);
		break;
	default:
		return std::wstring(L"<INVALID>");
	}
	return temp;
}

std::wstring	PatchProcessor::GetCurrPatchName(int patchN,int freq,Zeus::MapType mapT,std::wstring& Directory) const
{
	std::wstring			tDir(Zeus::CorrectDir((GlbVars_.DirIn_ +  GetFileTypeDir(mapT))));

	if(GlbVars_.Data2Buffer_)
	{
		tDir = Zeus::CorrectDir(GlbVars_.DirBuffer_ + GetFileTypeDir(mapT),1);
	}

	Directory	= tDir;
	return GetFileTypePrefix(mapT) + Zeus::PutNumber2Txt(GlbVars_.NSize_) + std::wstring(L"_") + Zeus::PutNumber2Txt(freq) + std::wstring(L"_") + Zeus::PutNumber2Txt(patchN);
}

void			PatchProcessor::InitializeFourierFarm(void)
{
	int i;
	for(i = GlbVars_.PatchSz_;i > MINIMUMDUOCIMATION;i >>= 1)
	{
		Zeus::FourierMachine*	f(new Zeus::FourierMachine(i));
		f->Initialize();
		FourierFarm_.push_back(f);
	}
}

void			PatchProcessor::Initialise(void)
{
	InitializeFourierFarm();
	SetNoiseDevice();
	MapsN_ = GlbVars_.N_ObsPlanes_;
	delete Background_;
	Background_ = new BackgroundProcessor(PatchNumber_,*this);
	Background_->Initialize();
	Zone_ = (PlanckInfo::Instance())->GetZone(PatchNumber_);
	Zone_->Initialise();
	delete	OddsEval_;
	OddsEval_	= new OddsEval(*this);
	delete	BrentMin_;
//	NSnakesDim_	= (GlbVars_.Use2DFormula_?NPARAMS_2D:NPARAMS_1D);
	NSnakesDim_	= NPARAMS_1D;
	BrentMin_	= new BrentLineMinimizer(OddsEval_,NSnakesDim_,LMBRENTMAXITER,LMBRENTTOL);
	SetupPowell();
	// 
	GammaFunct_.setbetaAB((3.0-GlbVars_.ProfParam_.MNFW_gamma_)/GlbVars_.ProfParam_.MNFW_alpha_,(GlbVars_.ProfParam_.MNFW_beta_-3.0)/GlbVars_.ProfParam_.MNFW_alpha_);
	y0EvalIntegrator_ =	new Zeus::SimpsonIntegrator<double,PatchProcessor>(2,MF_MAXITER,this,&PatchProcessor::y0Funct,1.0e-3,true);
	Evaly0_denominator_ = y0Cylindircal(GlbVars_.ProfParam_.VirialRatio_) * (GlbVars_.PixSz_ * GlbVars_.PixSz_  * SR2ARCMIN2);
}
//
void			PatchProcessor::MakeWhiteningStructs(double& LikeSurfISNR2Cte)
{
	PatchFilterPropsCache_.clear();
//	ReadObsMaps((GlbVars_.N_PriorPlanes_ && (CurrentIter_ != (GlbVars_.TwoStepsDetection_ - 1)))?Zeus::MAPTYPE_BACKGROUND:Zeus::MAPTYPE_OBSERVATION);
//	ComputeStats();
//	MaskMaps((GlbVars_.SZ_)?MASK_ALL:MASK_BORDER);
	Background_->MakeWhiteningStructs();

	AntennaFourVectors_	= Background_->WhiteFourierMaps(Zone_->GetAntVectorsRef());
	Zeus::PowerFunctorVect<double>	DummyGcc(ApodFunct_.GetPwSpApodCorrect());
	AntennaSFE_			= Transform(Loki::Type2Type<double>(),AntennaFourVectors_,DummyGcc,Zeus::UB_NOUSE);

	*(AntennaSFE_.GetInnerData().begin()) = 0.0;

	//Zeus::DumpInOut_2d("Filter",256,129,514,0,AntennaSFE_.GetInnerData().begin(),1.0);
	AntennaSFE_real_	= Fourier2Real(AntennaSFE_);
	Zeus::MultFunctor<double>	DummyGcc1(1.0 / static_cast<double>(GlbVars_.PatchSz_ * GlbVars_.PatchSz_));
	AntennaSFE_real_.Transform(DummyGcc1,0,Zeus::UB_NOUSE);

	//Zeus::DumpInOut_2d("FilterRealSpace",256,256,256,0,AntennaSFE_real_.GetInnerData().begin(),1.0);

	LikeSurfISNR2Cte = *(AntennaSFE_real_.GetInnerData().begin());

	PatchFilterPropsCache_.insert(PatchFilterPropsCacheType::value_type(ObjFilterParams(),PatchFilterPropsType(LikeSurfISNR2Cte)));
}

void			PatchProcessor::EvalEquivBeamSigmaPix(int& FullWith,int& HalfMax)
{
	const   Zeus::RealPlane<double>::DataInnerType&	InnerArray(AntennaSFE_real_.GetInnerData());

	Zeus::LArr2D<double>::const_iterator	const pix00(InnerArray.begin());
	Zeus::LArr2D<double>::const_iterator	pixPiv(pix00 + (InnerArray.getPtrMetric()>>1) - 1);
	Zeus::LArr2D<double>::const_iterator	const pixEndRow(pixPiv + 1);

	Zeus::LArr2D<double>::const_iterator	const pixEnd(pix00 - 1);

	double val(std::abs(*pix00) * EQUIVBEAMSIGMAFACTOR);

//TODO OpenMP
	for(;(pixPiv != pixEnd) && (std::abs(*pixPiv) < val);--pixPiv)
		;
	++pixPiv;
	FullWith = static_cast<int>(pixPiv - pix00);

	val = std::abs(*pix00) * 0.5;
	pixPiv = pix00;
//TODO OpenMP
	for(;(pixPiv != pixEndRow) && (std::abs(*pixPiv) > val);++pixPiv)
		;
	EqBeamSigmaPix_ = HalfMax = static_cast<int>(pixPiv - pix00);
}

void			PatchProcessor::SetupPowell(void)
{
	MinArrayType		InitDirScales(NSnakesDim_);
	std::vector<int>	AllowSubSpaces(NSnakesDim_);
	delete	PwSnake_;
	PwSnake_	= new PowellSnake(BrentMin_,OddsEval_,NSnakesDim_,PWMAXITER,PWTOL);
	InitDirScales.at(0) = 1.0;
	InitDirScales.at(1) = 1.0;
	InitDirScales.at(2) = 1.0;
	if(NSnakesDim_ > 3) InitDirScales.at(3) = 1.0;

	PwSnake_->SetDirectionScales(InitDirScales);

	if(GlbVars_.SearchPos_)
	{
		AllowSubSpaces.at(0) = 1;
		AllowSubSpaces.at(1) = 1;
	}
	else
	{
		AllowSubSpaces.at(0) = 0;
		AllowSubSpaces.at(1) = 0;
	}
	AllowSubSpaces.at(2) = 1;
	if(NSnakesDim_ > 3) AllowSubSpaces.at(3) = 1;
	PwSnake_->SetAllowedSubSpaces(AllowSubSpaces);
}

void			PatchProcessor::MakeMainSurfaces(double LikeSurfISNR2Cte,double& LikeSurfSigma, int freeMem)
{
	int						NPix;
	double					bias;
	MapsCollFourierType		MainFourVectors;
	MapsCollFourierType		FModeVectors;

	FModeVectors.MakeNewSz(GlbVars_.PatchSz_ >> 1,false);
	MakeFourSurfaces();
	PutMap2FModeVectors(FModeVectors);
	ReleaseFModeMaps();
	MainFourVectors = Background_->WhiteFourierMaps(FModeVectors);
	FModeVectors.Release();
	MapsMainSurfFour_ = Multiply(Loki::Type2Type<std::complex<double> >(),AntennaFourVectors_,MainFourVectors,Zeus::UB_NOUSE);
	MainFourVectors.Release();
	if(GlbVars_.SZ_ && freeMem)
	{
		AntennaFourVectors_.Release();
		Background_->ReleaseWhiteningStructs();
	}
	
	MapsMainSurfFour_.AverageMode1();	

	LikelihoodSurf_ = Fourier2Real(MapsMainSurfFour_);

	//DumpInOut_2d("MapsMainSurfFour_",512,256,257,1,MapsMainSurfFour_.GetInnerData().begin(),1.0,1);

	Zeus::MultFunctor<double>	DummyGcc(ApodFunct_.GetPwSpApodCorrect()/static_cast<double>(GlbVars_.PatchSz_ * GlbVars_.PatchSz_));
	LikelihoodSurf_.Transform(DummyGcc,0,Zeus::UB_NOUSE);


	double SqrtISNR(std::sqrt(LikeSurfISNR2Cte));

	//const std::complex<double> *Dummy(MapsMainSurfFour_.GetInnerData().begin());
	//DumpInOut_2d("MapsMainSurfFour",256,256,257,0,Dummy,1.0 / SqrtISNR,1);
	//exit(-1);
	{
		(PlanckInfo::Instance())->ComputeStats1Map(LikelihoodSurf_,DEFREJECTLEV_NPMAP,bias,LikeSurfSigma,NPix,SqrtISNR);
	}

	double SigIsnrRatio(LikeSurfSigma / SqrtISNR);

	wchar_t	buffer[BUFFERMAXCHAR];	
	PRINTINTOBUFFERFUNCT
		(buffer,BUFFERMAXCHAR,FILTREPFORMAT,
		0.0,
		LikeSurfSigma,
		SqrtISNR,
		SigIsnrRatio
		);
	(Zeus::ConManager::Instance())->PrintStr2Console(buffer);
}

void			PatchProcessor::FindSourcePositions(double LikeSurfISNR2Cte,int StoreNoise)
{
	//TODO
	ObjFilterParams		pLikeParams;
	RealPlaneSurfType	tempLikelihood;
	double				LikeSurfISNR2AuxCte,LikeSurfAuxSigma;
	std::vector<double>	tGauss;
	const	int			ScalesCollSz(GlbVars_.ScalesColl_.size());

	for(int i=0;i < ScalesCollSz;++i)
	{
		tGauss.push_back(MakeLikelihoodAtScale(i,tempLikelihood,LikeSurfISNR2AuxCte,LikeSurfAuxSigma,pLikeParams,StoreNoise));
		ScanLikeSurface(pLikeParams,tempLikelihood,LikeSurfISNR2AuxCte,LikeSurfAuxSigma);
		if(GlbVars_.NotAlignedObjs_)
		{
			EvalPosMaxima(tempLikelihood,Zeus::toInt(GlbVars_.ScalesColl_[i] + 0.5));
		}
	}
	CheckSubSources(tGauss);
	SelectMaxima();
}
//
void			PatchProcessor::CheckSubSources(std::vector<double>& tGauss)
{
	wchar_t		buffer[BUFFERMAXCHAR];
	double t(0.0);
	std::vector<double>::const_iterator	piv(tGauss.begin());
	std::vector<double>::const_iterator	const end(tGauss.end());

	for(;piv != end;++piv)
	{
		if((*piv)<1.0)
		{t += std::log(*piv);}
	}

	t	/= (static_cast<double>(tGauss.size()));

	if(t<GlbVars_.SSubGaussLevel_)
	{(PlanckInfo::Instance())->SetSubtractSources(1);}
	else
	{(PlanckInfo::Instance())->SetSubtractSources(0);}

	PRINTINTOBUFFERFUNCT
		(buffer,BUFFERMAXCHAR,L"Average log Gaussianity -> %5.2g, Subtract sources ? -> %ls\n",t,((t<GlbVars_.SSubGaussLevel_)?L"YES":L"NO"));

	(Zeus::ConManager::Instance())->PrintStr2Console(buffer);
}
//
void			PatchProcessor::FindOptimalScaleAtPixel(MaxLikeSurfMaximum& LikePeak)
{
	const	int		ScalesCollSz(GlbVars_.ScalesColl_.size());	
	int				MaxIndex, MaxRadBin, MaxSign, sign;
	double			MaxLike(-1.0e32),MaxLiketemp,ISNR,MaxISNR;	

	for(int i=0;i<ScalesCollSz;++i)
	{
		const ObjFilterParams	RadBin(Zeus::toInt(GlbVars_.ScalesColl_[i] + 0.5));
		MaxLiketemp = OddsEval_->GetOddsResultRawOptimAmp(LikePeak.Coords_.YPix_, LikePeak.Coords_.XPix_, RadBin, ISNR, sign);
		if(MaxLiketemp > MaxLike)
		{
			MaxSign		= sign;
			MaxISNR		= ISNR;
			MaxLike		= MaxLiketemp;
			MaxIndex	= i;
			MaxRadBin	= RadBin.BinScale_;
		}
	}
	ISNR						= std::abs(ISNR);
	LikePeak.LikeRatio_			= std::abs(MaxLike);
	LikePeak.NormAmpl_			= std::sqrt(2.0 * LikePeak.LikeRatio_);
	LikePeak.Correlation_		= LikePeak.NormAmpl_ * ISNR;
	LikePeak.flux_				= (LikePeak.NormAmpl_ / ISNR) * static_cast<double>(sign);
	LikePeak.FiltParams_.BinScale_ = MaxRadBin;
	LikePeak.GaussianIndex_		= ScalesProps_[MaxIndex].Gaussianity_;
}
//
double			PatchProcessor::MakeLikelihoodAtRealScale(const Zeus::ObjFilterRealParams& params,RealPlaneSurfType& tempLikelihood)
{
	// get the profile
	SrcInfoType						ObjSurf(Zone_->GetObjOfParamsRealValue(params));
	
	// Zeus::ObjFilterRealParams& params
	// -> Fourier domain
	MapsMainSurfFourType			tempFour(Real2Fourier(ObjSurf.surf_));
	
	// Take the power
	Zeus::PowerFunctor<double>		DummyGcc1;
	AntMainSurfType					ObjTemp(Transform(Loki::Type2Type<double>(),tempFour,DummyGcc1,Zeus::UB_NOUSE));
	
	// Multiply by the antenna single frequency equivalent
	MultiplyInPlace(ObjTemp,AntennaSFE_,Zeus::UB_NOUSE);
	
	// Get the sum -> ISNR2
	double ISNR2,dummy;
	ObjTemp.IntegrateSurface(Zeus::SEL_SUM,ISNR2,dummy,Zeus::UB_NOUSE);
	
	//Multiply the object by the SFE of the maps
	MultiplyInPlace(tempFour,MapsMainSurfFour_,Zeus::UB_NOUSE);

	// -> pixel domain
	tempLikelihood = Fourier2Real(tempFour);

	// Normalise
	Zeus::MultFunctor<double>	DummyGcc2(ApodFunct_.GetPwSpApodCorrect()/static_cast<double>(GlbVars_.PatchSz_ * GlbVars_.PatchSz_));
	tempLikelihood.Transform(DummyGcc2,0,Zeus::UB_NOUSE);

	return ISNR2;
}
//
double			PatchProcessor::MakeLikelihoodAtScale(int BinIndex,RealPlaneSurfType& tempLikelihood,
													double&	LikeSurfISNR2AuxCte,double& LikeSurfAuxSigma,
													ObjFilterParams& objParam, int StoreNoise)
{
	wchar_t		buffer[BUFFERMAXCHAR];	
	char		bufferS[BUFFERMAXCHAR];	
	int			ScaleBin;


	if(GlbVars_.ScalesFillerType_)
	{
		ScaleBin	= Zeus::toInt(GlbVars_.ScalesColl_[BinIndex] + 0.5);
		objParam	= ObjFilterParams(ScaleBin);
	}
	else
	{
		objParam = Zone_->TranslateParamsInverse(Zeus::ObjFilterRealParams(GlbVars_.ScalesColl_[BinIndex]));
		ScaleBin	= BinIndex;
	}

	double		dummy;

	SrcInfoType ObjSurf(GlbVars_.ScalesFillerType_?Zone_->GetObjOfRadiusBin(ScaleBin):Zone_->GetObjOfRadiusRealValue(GlbVars_.ScalesColl_[BinIndex]));

	// Get the real space object

	//sprintf(bufferS,"%s_%d_","NgaiProfile",ScaleBin);
	//DumpInOut_2d(bufferS,256,256,512,0,ObjSurf.surf_.GetInnerData().begin(),1.0);

	MapsMainSurfFourType tempFour(Real2Fourier(ObjSurf.surf_));
	Zeus::PowerFunctor<double>	DummyGcc1;
	AntMainSurfType		ObjTemp(Transform(Loki::Type2Type<double>(),tempFour,DummyGcc1,Zeus::UB_NOUSE));

	//AntTemp has now the square of the Object fourier transform

	//sprintf(bufferS,"%s_%d_","Obj2Before",ScaleBin);
	//Zeus::DumpInOut_2d(bufferS,256,129,129,0,ObjTemp.GetInnerData().begin(),1.0);

	MultiplyInPlace(ObjTemp,AntennaSFE_,Zeus::UB_NOUSE);

	//Multilpy by the antenna single frequency equivalent surface
	//sprintf(bufferS,"%s_%d_","Obj2AntennaSFE",ScaleBin);
	//Zeus::DumpInOut_2d(bufferS,256,129,129,0,AntennaSFE_.GetInnerData().begin(),1.0);
	//sprintf(bufferS,"%s_%d_","Obj2Antenna",ScaleBin);
	//DumpInOut_2d(bufferS,256,129,257,0,ObjTemp.GetInnerData().begin(),1.0);

	ObjTemp.IntegrateSurface(Zeus::SEL_SUM,LikeSurfISNR2AuxCte,dummy,Zeus::UB_NOUSE);
	
	// Evaluate the ISNR

	// PatchFilterPropsCache_.insert(PatchFilterPropsCacheType::value_type(objParam,PatchFilterPropsType(LikeSurfISNR2AuxCte)));
	
	//put it in the cache

	MultiplyInPlace(tempFour,MapsMainSurfFour_,Zeus::UB_NOUSE);

	//Multiply the object by the SFE of the maps

	tempLikelihood = Fourier2Real(tempFour);

	//tempLikelihood = Fourier2Real(MapsMainSurfFour_);

	Zeus::MultFunctor<double>	DummyGcc2(ApodFunct_.GetPwSpApodCorrect()/static_cast<double>(GlbVars_.PatchSz_ * GlbVars_.PatchSz_));
	tempLikelihood.Transform(DummyGcc2,0,Zeus::UB_NOUSE);

	//sprintf(bufferS,"%s_%d_","tempLikelihood",ScaleBin);
	//DumpInOut_2d(bufferS,256,256,512,128 * 512 + 128,tempLikelihood.GetInnerData().begin(),1.0);

	// Adjust and normalise

	int	NPix;
	double SqrtISNR(std::sqrt(LikeSurfISNR2AuxCte));
	
	(PlanckInfo::Instance())->ComputeStats1Map(tempLikelihood,DEFREJECTLEV_NPMAP,dummy,LikeSurfAuxSigma,NPix,SqrtISNR);

	if(StoreNoise)
	{
		int dummy1,dummy2;
		const int	XCentralPix((PatchesGeoInfo_.Storage_[PatchNumber_].SrcXCoord_ >= 0)?PatchesGeoInfo_.Storage_[PatchNumber_].SrcXCoord_:(GlbVars_.PatchSz_ >> 1));
		const int	YCentralPix((PatchesGeoInfo_.Storage_[PatchNumber_].SrcYCoord_ >= 0)?PatchesGeoInfo_.Storage_[PatchNumber_].SrcYCoord_:(GlbVars_.PatchSz_ >> 1));

		ScalesProps_[BinIndex].Th1OverSigma_		= SqrtISNR;
		ScalesProps_[BinIndex].LikelihoodSigma_		= LikeSurfAuxSigma;
		ScalesProps_[BinIndex].LikeCentralPix_		= GetMaximumLike(tempLikelihood,YCentralPix,XCentralPix,((GlbVars_.NonBlindDetection_ > 0)?2:0),dummy1,dummy2);
	}
	
	ScalesProps_[BinIndex].SrcScale_	= ObjSurf.RealParams_.RealScale_;
	ScalesProps_[BinIndex].Gaussianity_ = LikeSurfAuxSigma / SqrtISNR;

	PRINTINTOBUFFERFUNCT
		(buffer,BUFFERMAXCHAR,FILTREPFORMAT,
		ObjSurf.RealParams_.RealScale_,
		LikeSurfAuxSigma,
		SqrtISNR,
		ScalesProps_[BinIndex].Gaussianity_
		);

	(Zeus::ConManager::Instance())->PrintStr2Console(buffer);

//	sprintf(bufferS,"%s_%d","Likelihood",ScaleBin);
//	Zeus::DumpInOutFits_2d(bufferS, tempLikelihood.GetInnerData().getPtrMetric(), tempLikelihood.GetInnerData().getPtrMetric(), tempLikelihood.GetInnerData().getPtrMetric(), 0, tempLikelihood.GetInnerData().begin(), 1.0 / SqrtISNR);
	return 	ScalesProps_[BinIndex].Gaussianity_;
}
//
void			PatchProcessor::ReadObsMaps(int dataset)
{
	ObsFreqsCollType::const_iterator	piv(GlbVars_.FreqsColl_.begin());
	ObsFreqsCollType::const_iterator	const end(GlbVars_.FreqsColl_.end());
	ObsPriorMapType						temp;
	double								Dummy;

	ObsPriorMapsColl_.clear();
	(Zeus::ConManager::Instance())->PrintStr2Console(STARTREADINGMAPS);

	for(;piv != end;++piv)
	{
		temp.Freq_ = piv->freq_;
		temp.Order_ = GetFreqOrder(temp.Freq_);
		temp.ws_.MakeNewSz(GlbVars_.PatchSz_,false); //Needs initialization
		Read1ObsMap(PatchNumber_, *piv, static_cast<Zeus::MapType>(dataset), temp.ws_.GetInnerData());
		temp.ws_.RemoveBias(Dummy);

		//char	bufferS[BUFFERMAXCHAR];	
		//sprintf(bufferS,"%s_%d","Map",temp.Freq_);
		//Zeus::DumpInOut_2d(bufferS,512,512,512,0,temp.ws_.GetInnerData().begin(),1.0);

		ObsPriorMapsColl_.push_back(temp);
	}

	std::sort(ObsPriorMapsColl_.begin(),ObsPriorMapsColl_.end());
	AdjustMapsOrderIndex();
	(Zeus::ConManager::Instance())->PrintStr2Console(ENDREADINGMAPS);

}
void			PatchProcessor::Read1ObsMap(int patchN, const ObsFreqsType& freq, Zeus::MapType mapT, Zeus::LArr2D<double>& ws) const
{
	std::wstring	Directory;
	std::wstring	FName(GetCurrPatchName(patchN, freq.sign_, mapT, Directory));

	std::auto_ptr<Zeus::GenCollReader<Zeus::LArr2D<double> > >
		FReader(Zeus::GetWrkSpFileReaderHandler(Loki::Type2Type<Zeus::LArr2D<double> >(),
		GlbVars_.Data2Buffer_?1000:GlbVars_.ContextID_,FName,GlbVars_.PatchSz_,GlbVars_.PatchSz_,Directory));
	FReader->Initialize();
	FReader->Read();
	FReader->Release(ws);
}

void			PatchProcessor::SetNoiseDevice(void)
{
	delete	NoiseGen_;
	NoiseGen_ = new Zeus::PwSCoreRandGen(-1 /* use time to seed the generator*/);
}
//
void 			PatchProcessor::ScanLikeSurface(const ObjFilterParams& FilterParams,const RealPlaneSurfType& surf,double ISNR2,double Sigma)
{
	const unsigned long ptrMetric(surf.GetInnerData().getPtrMetric());
	unsigned long OffSet;
	const double * const LikeSurfOrg_ptr(surf.GetInnerData().begin());
	const double *LikeSurf0_ptr(LikeSurfOrg_ptr);
	const double *LikeSurf1_ptr(LikeSurf0_ptr + 1);
	const double *LikeSurf2_ptr(LikeSurf0_ptr + 2);
	const double *LikeSurf3_ptr(LikeSurf0_ptr + ptrMetric);
	const double *LikeSurf4_ptr(LikeSurf3_ptr + 1);
	const double *LikeSurf5_ptr(LikeSurf3_ptr + 2);
	const double *LikeSurf6_ptr(LikeSurf3_ptr + ptrMetric);
	const double *LikeSurf7_ptr(LikeSurf6_ptr + 1);
	const double *LikeSurf8_ptr(LikeSurf6_ptr + 2);
	const double *LikeSurfEnd_ptr(LikeSurf0_ptr + surf.GetInnerData().getSz());
	const double sqtrISNR(std::sqrt(ISNR2));
	const double ThreeSigmaThreshold(SIGMAPEAKSELECT * ((Sigma < sqtrISNR)?Sigma:sqtrISNR));
	const double GaussianIndex(Sigma/sqtrISNR);

	MaxLikeSurfMaximum	temp;

	temp.FiltParams_ = FilterParams;
	
	for(;LikeSurf8_ptr != LikeSurfEnd_ptr;
		++LikeSurf0_ptr,++LikeSurf1_ptr,++LikeSurf2_ptr,
		++LikeSurf3_ptr,++LikeSurf4_ptr,++LikeSurf5_ptr,
		++LikeSurf6_ptr,++LikeSurf7_ptr,++LikeSurf8_ptr
		)
	{

		if( (*LikeSurf4_ptr > *LikeSurf0_ptr) && (*LikeSurf4_ptr > *LikeSurf1_ptr) && (*LikeSurf4_ptr > *LikeSurf2_ptr) &&
			(*LikeSurf4_ptr > *LikeSurf3_ptr) && (*LikeSurf4_ptr > *LikeSurf5_ptr) &&
			(*LikeSurf4_ptr > *LikeSurf6_ptr) && (*LikeSurf4_ptr > *LikeSurf7_ptr) && (*LikeSurf4_ptr > *LikeSurf8_ptr)
		)
		{
			if(*LikeSurf4_ptr <= ThreeSigmaThreshold)
				continue;

			OffSet = static_cast<unsigned long>(LikeSurf4_ptr - LikeSurfOrg_ptr);

			temp.Coords_.XPix_ = (OffSet % ptrMetric);
			temp.Coords_.YPix_ = (OffSet / ptrMetric);

			if(!CheckLikeMaximum(temp))
				continue;

			temp.Coords_.YCoord_ = static_cast<double>(temp.Coords_.YPix_);
			temp.Coords_.XCoord_ = static_cast<double>(temp.Coords_.XPix_);
			temp.GaussianIndex_ = GaussianIndex;
			temp.Correlation_	= *LikeSurf4_ptr;
			temp.NormAmpl_		= temp.Correlation_ / std::sqrt(std::abs(ISNR2));
			temp.LikeRatio_		= (temp.Correlation_ * temp.Correlation_) / (2.0 * ISNR2);
			if(temp.Correlation_ < 0.0) temp.LikeRatio_ = -temp.LikeRatio_;
			temp.flux_			=  temp.Correlation_ / ISNR2;
			
			MaximaColl_.push_back(temp);
		}
	}
}

void			PatchProcessor::SelectMaxima(void)
{
	std::sort(MaximaColl_.begin(),MaximaColl_.end());
	MarkLowerFluxes();
	MaxLikeSurfMaximaCollType::iterator	newEnd(std::remove_if(MaximaColl_.begin(),MaximaColl_.end(),LowerFluxesRemoveMark));
	if(newEnd != MaximaColl_.end())
	{MaximaColl_.erase(newEnd,MaximaColl_.end());}

	if((CurrentIter_ == (GlbVars_.TwoStepsDetection_ - 1)) && GlbVars_.NonBlindDetection_)
	{
		RemoveAllButCentralPeak();
		MaxLikeSurfMaximaCollType::iterator	newEnd(std::remove_if(MaximaColl_.begin(),MaximaColl_.end(),LowerFluxesRemoveMark));
		if(newEnd != MaximaColl_.end())
		{MaximaColl_.erase(newEnd,MaximaColl_.end());}
	}
}
//
double			PatchProcessor::GetMaximumLike(const RealPlaneSurfType& tempLikelihood,int YCent,int XCent,int Tol,int& YMax,int& XMax)
{
	double LikeMaxValue(Zeus::logZERO);
	int YInfLimt,YSupLimt,XInfLimt,XSupLimt;
	const unsigned long ptrMetric(tempLikelihood.GetInnerData().getPtrMetric());
	long Line;
	const double * const LikeSurfOrg_ptr(tempLikelihood.GetInnerData().begin());

	if((YInfLimt = YCent-Tol)<0) YInfLimt=0;
	if((XInfLimt = XCent-Tol)<0) XInfLimt=0;
	if((YSupLimt = YCent+Tol)>=GlbVars_.PatchSz_) YSupLimt = GlbVars_.PatchSz_-1;
	if((XSupLimt = XCent+Tol)>=GlbVars_.PatchSz_) XSupLimt = GlbVars_.PatchSz_-1;

	for(int j=YInfLimt;j<=YSupLimt;++j)
	{
		Line = (ptrMetric*j);
		for(int i=XInfLimt;i<=XSupLimt;++i)
		{
			if(*(LikeSurfOrg_ptr + (Line + i)) > LikeMaxValue)
			{
				LikeMaxValue = *(LikeSurfOrg_ptr + (Line + i));
				YMax = j;XMax = i;
			}
		}	
	}
	return LikeMaxValue;
}
//
void			PatchProcessor::RemoveAllButCentralPeak(void)
{
	if(MaximaColl_.empty())
		return;

	double	tLikeRatio(-1.0e32);

	int	XCentralPix((PatchesGeoInfo_.Storage_[PatchNumber_].SrcXCoord_ >= 0)?PatchesGeoInfo_.Storage_[PatchNumber_].SrcXCoord_:(GlbVars_.PatchSz_ >> 1));
	int	YCentralPix((PatchesGeoInfo_.Storage_[PatchNumber_].SrcYCoord_ >= 0)?PatchesGeoInfo_.Storage_[PatchNumber_].SrcYCoord_:(GlbVars_.PatchSz_ >> 1));

	MaxLikeSurfMaximaCollType::iterator			piv(MaximaColl_.begin());
	MaxLikeSurfMaximaCollType::iterator			const end(MaximaColl_.end());
	int		LimitPixels(GlbVars_.NonBlindDetection_ >> 4);
	// GlbVars_.ScalesFillerType_==2 -> hardcoded scales
	const int	PriorType((GlbVars_.ScalesFillerType_==2)? 2:(GlbVars_.NonBlindDetection_ % 8));
	
	MaxLikeSurfMaximaCollType::iterator			pivAux(end);
	
	for(;piv != end;++piv)
	{
		if((std::abs(piv->Coords_.YPix_ - YCentralPix) <= LimitPixels)	&&
			(std::abs(piv->Coords_.XPix_ - XCentralPix) <= LimitPixels)
			)
		{
			if(piv->LikeRatio_ > tLikeRatio)
			{
				tLikeRatio = piv->LikeRatio_;
				pivAux	= piv;
			}

		}
	}

	if(pivAux == end)
	{
		MaxLikeSurfMaximaCollType::value_type	temp;
/*
		if(
			!PriorType										 ||
			(GlbVars_.AssessmentKind_	&& ((PriorType == 2) || (PriorType == 6)))
		)

		{
			pivAux = MaximaColl_.begin();
		}
		else
*/
		{
			temp.Coords_.XPix_		= XCentralPix;
			temp.Coords_.YPix_		= YCentralPix;	
			temp.Coords_.XCoord_	= static_cast<double>(XCentralPix);	
			temp.Coords_.YCoord_	= static_cast<double>(YCentralPix);	
			
			if(!(GlbVars_.AssessmentKind_) && (PriorType & 0x02))
			{temp.NonBlindTarget_ = -1;}
			else
			{temp.NonBlindTarget_ = 1;}

			FindOptimalScaleAtPixel(temp);
			MaximaColl_.push_back(temp);
			pivAux = MaximaColl_.end();
			return;
		}
	}
	else
	{
		if(!(GlbVars_.AssessmentKind_) && (PriorType & 0x02))
		{pivAux->NonBlindTarget_ = -1;}
		else
		{pivAux->NonBlindTarget_ = 1;}
		++pivAux;
	}

	for(;pivAux != end;++pivAux)
	{
		pivAux->Removed_ = 1;
	}
}

void			PatchProcessor::MarkLowerFluxes(void)
{
	MaxLikeSurfMaximaCollType::iterator			piv(MaximaColl_.begin());
	MaxLikeSurfMaximaCollType::const_iterator	const end(MaximaColl_.end());
	MaxLikeSurfMaximaCollType::iterator			pivAux;
	for(;piv != end;++piv)
	{
		if(piv->Removed_)
			continue;
		for(pivAux = (piv+1);pivAux != end;++pivAux)
		{
			if((std::abs(pivAux->Coords_.YPix_ - piv->Coords_.YPix_) <= 3)	&&
				(std::abs(pivAux->Coords_.XPix_ - piv->Coords_.XPix_) <= 3)
				)
			{
				pivAux->Removed_ = 1;
			}
		}
	}
}


void			PatchProcessor::EvalPosMaxima(const RealPlaneSurfType& surf,int BinScale)
{
	MaxLikeSurfMaximaCollType::iterator			piv(MaximaColl_.begin());
	MaxLikeSurfMaximaCollType::const_iterator	const end(MaximaColl_.end());

	const double *	const	Org(surf.GetInnerData().begin());
	const unsigned long		ptrMetric(surf.GetInnerData().getPtrMetric());
	const double *			Like_ptr;
	double					val[9];


	for(;piv != end;++piv)
	{
		if((piv->FiltParams_).BinScale_ != BinScale)
			continue;

		Like_ptr = Org + (piv->Coords_.XPix_ + (piv->Coords_.YPix_ * ptrMetric));

		val[0] = *(Like_ptr - (ptrMetric + 1));
		val[1] = *(Like_ptr - ptrMetric);
		val[2] = *(Like_ptr - (ptrMetric - 1));
		val[3] = *(Like_ptr - 1);
		val[4] = *Like_ptr;
		val[5] = *(Like_ptr + 1);
		val[6] = *(Like_ptr + (ptrMetric - 1));
		val[7] = *(Like_ptr + ptrMetric);
		val[8] = *(Like_ptr + ptrMetric + 1);
		switch(GlbVars_.NotAlignedObjs_)
		{
		case	1:
			Centroid(val);
			break;
		case	2:
			MaximaByQuadInterpolation(val);
			break;
		default:
			Centroid(val);
		}
		piv->Coords_.YCoord_ = static_cast<double>(piv->Coords_.YPix_) + val[1];
		piv->Coords_.XCoord_ = static_cast<double>(piv->Coords_.XPix_) + val[0];
	}
}


RealPlaneSurfType	PatchProcessor::MakeWhitenedSource(double flux,double DeltaY, double DeltaX,const Zeus::ObjFilterRealParams& ObjParams,double& CentralPixAmpl)
{
	SrcInfoType				ObjSurf(Zone_->GetSrcObjRaw(ObjParams));
	MapsMainSurfFourType	tempFour(Real2Fourier(ObjSurf.surf_));
	MultiplyInPlace(tempFour,AntennaSFE_,Zeus::UB_NOUSE);
	tempFour.ShiftInPlace(DeltaY,DeltaX);
	RealPlaneSurfType tMapSource(Fourier2Real(tempFour));
	Zeus::MultFunctor<double>	DummyGcc(flux / static_cast<double>(GlbVars_.PatchSz_ * GlbVars_.PatchSz_));
	tMapSource.Transform(DummyGcc,0,Zeus::UB_NOUSE);

	//DumpInOut_2d("WhitenedSource",256,256,256,tMapSource.GetInnerData().begin(),1.0 / LikeSurfSigma_);

	CentralPixAmpl = ObjSurf.CentralPixAmpl_;
	return tMapSource;
}
//
PatchProcessor::~PatchProcessor(void)
{
	delete	y0EvalIntegrator_;
	delete	NestSampler_;
	delete	PwSnake_;
	delete	OddsEval_;
	delete	BrentMin_;
	delete	Background_;
	delete	NoiseGen_;
	FourierFarmType::iterator		piv(FourierFarm_.begin());
	FourierFarmType::const_iterator	const end(FourierFarm_.end());
	for(;piv!=end;++piv)
	{delete *piv;}
}
//
bool	PatchProcessor::Find1PeakOptimalParams(MaxLikeSurfMaximum& LikePeak,Zeus::PeakType &PeakResult)
{
	MinArrayType	InitGuessFinal(NSnakesDim_);

	if(GlbVars_.NotAlignedObjs_)
	{
		InitGuessFinal[0] = LikePeak.Coords_.YCoord_;
		InitGuessFinal[1] = LikePeak.Coords_.XCoord_;
	}
	else
	{
		InitGuessFinal[0] = static_cast<double>(LikePeak.Coords_.YPix_);
		InitGuessFinal[1] = static_cast<double>(LikePeak.Coords_.XPix_);
	}

	InitGuessFinal[2]	=  static_cast<double>(LikePeak.FiltParams_.BinScale_);

	PeakResult.InitSrcFluxEst_		= LikePeak.flux_;
	PeakResult.InitCorrelation_		= LikePeak.Correlation_;
	PeakResult.GaussianIndex_		= LikePeak.GaussianIndex_;

	if(NSnakesDim_ > 3)
	{
		InitGuessFinal[3] = PeakResult.InitSrcFluxEst_;
	}

	if(
		!(PwSnake_->PowellMinimunSearch(InitGuessFinal,PeakResult.Odds_)) ||
		(PeakResult.Odds_ > 0.0)
	)
	{return false;}

	PeakResult.Odds_			= -PeakResult.Odds_;
	PeakResult.Pos_.YCoord_		= InitGuessFinal[0];
	PeakResult.Pos_.XCoord_		= InitGuessFinal[1];
	PeakResult.Pos_.YPix_		= Zeus::toInt(PeakResult.Pos_.YCoord_ + 0.5);
	PeakResult.Pos_.XPix_		= Zeus::toInt(PeakResult.Pos_.XCoord_ + 0.5);
	PeakResult.UnTransParams_	= Zeus::ObjFilterRealParams(InitGuessFinal[2]);
	PeakResult.RealParams_		= Zone_->TranslateParams(PeakResult.UnTransParams_);
	PeakResult.ISNR2_			= GetISNR(PeakResult.UnTransParams_);

	PeakResult.SrcAmplNormalised_ = std::sqrt(2.0 * PeakResult.Odds_);

	if(NSnakesDim_ > 3)
	{
		PeakResult.SrcFlux_ = InitGuessFinal[3];
	}
	else
	{
		PeakResult.SrcFlux_ = PeakResult.SrcAmplNormalised_ / std::sqrt(PeakResult.ISNR2_);
	}

	return true;
}

int		PatchProcessor::EvalExpectNCorrectCtes(double& lnCv,double& lnCu)
{
	if(!DetectedObj_)
	{
		lnCv	= Zeus::logZERO;
		lnCu	= Zeus::logZERO;
		return 1;
	}
	if(DetectedObj_ == 1)
	{
		lnCv	= Zeus::logZERO;
		lnCu	= -GlbVars_.PriorAvObjectsPatch_;	
		return 1;
	}

	const double lnMiu(std::log(GlbVars_.PriorAvObjectsPatch_));

	lnCu = 0.0;
	lnCv = Zeus::logZERO;

	double t;
	for(int k=1;k<DetectedObj_;++k)
	{
		t		= k*lnMiu - Zeus::logfactorial(k);
		lnCu	= Zeus::AddLog(lnCu,t);
		lnCv	= Zeus::AddLog(lnCv,t+std::log(static_cast<double>(k)));
	}
	
	lnCu   -= GlbVars_.PriorAvObjectsPatch_;
	lnCv   -= GlbVars_.PriorAvObjectsPatch_;

	return 1;
}

double		PatchProcessor::Eval_lnModelRatio(Zeus::PeakType & PeakRes)
{
	double	CurrExpectObjPeaks;
	double	lnCv,lnCu;

	EvalExpectNCorrectCtes(lnCv,lnCu);

	double t(Zeus::SubLog(std::log(GlbVars_.PriorAvObjectsPatch_),lnCv) - Zeus::SubLog(0.0,lnCu));
	
	if(t >= 6.90)
		return Zeus::logZERO;

	CurrExpectObjPeaks = std::exp(t);

	CurrExpectObjPeaks -= static_cast<double>(DetectedObj_);
	if(CurrExpectObjPeaks < 0.0)
		return Zeus::logZERO;

	// Correct expected number of peaks
	// ln(4) 4 times more background peaks than true peaks

	return std::log(CurrExpectObjPeaks / ExpectedBackgPeaks_);
}
//
double		PatchProcessor::Eval_ExpectedBackgPeaks(void)
{
	ExpectedBackgPeaks_ = 4.0 * GlbVars_.PriorAvObjectsPatch_;
	return ExpectedBackgPeaks_;
}

int	PatchProcessor::BayesLaplaceAssess(Zeus::PeakType & PeakRes)
{
	wchar_t	buffer[BUFFERMAXCHAR];

	StatProps						Evidence;
	PeakRes.JF_lnRhoTh_				= 0.0 ;
	PeakRes.JF_lnRho_				= Zeus::logZERO;
	PeakRes.PK_BayesDetectStat_		= Zeus::PeakType::PK_DET_ERROR;
	PeakRes.JF_lnEvidence_			= Zeus::logZERO;
	PeakRes.JF_lnEvidenceErrBar_	= Zeus::logZERO;
	PeakRes.JF_lnModelRatio_		= Zeus::logZERO;
	PeakRes.JF_lnFormFactor_		= Zeus::logZERO;
	PeakRes.JF_UsedSigmaThreshold_	= -1.0;

	if(DetectedObj_ > ((8.0 * std::sqrt(GlbVars_.PriorAvObjectsPatch_)) + GlbVars_.PriorAvObjectsPatch_))
		return Zeus::PeakType::PK_DET_ERROR;

	NestSampler_->NextPeak(&PeakRes,PatchesGeoInfo_.Storage_[PatchNumber_]);

	const int EvStat(NestSampler_->GetLogEvidence(Evidence,(NestSampler_->GetEvCorrectingCte() - PENALTY1SOURCE - (THRESHOLD10PERC/2.0))));

	if(EvStat)
	{
		PRINTINTOBUFFERFUNCT
			(buffer,BUFFERMAXCHAR,
			L"Nested sampler did not converge. Patch No -> %4d, X -> %3d , Y -> %3d\n",
			PeakRes.PatchNumber_,PeakRes.Pos_.XPix_,PeakRes.Pos_.YPix_
			);

		(Zeus::ConManager::Instance())->PrintErr2Console(buffer);
	}

	PeakRes.JF_lnEvidence_			= Evidence.Mean_  + NestSampler_->GetEvCorrectingCte();
	PeakRes.JF_lnEvidenceErrBar_	= Evidence.StDev_;

	const double mL(NestSampler_->GetMaxLogLike());

	if(PeakRes.Odds_ != mL)
	{
		PeakRes.Odds_ = mL;
		PeakRes.SrcAmplNormalised_ = std::sqrt(2.0 * PeakRes.Odds_);
		NP_Assessment(PeakRes);
	}

	PeakRes.JF_lnFormFactor_		= PeakRes.JF_lnEvidence_ - PeakRes.Odds_;
	
	if(!(GlbVars_.NonBlindDetection_))
	{
		PeakRes.JF_lnModelRatio_		= -Eval_lnModelRatio(PeakRes);
	}
	else
	{
		PeakRes.JF_lnModelRatio_ = 0.0;
	}

	if(PeakRes.JF_lnModelRatio_ != Zeus::logINF)
	{
		const double t(PeakRes.JF_lnModelRatio_ - PeakRes.JF_lnFormFactor_);
		if(t > 0.0)
			PeakRes.JF_UsedSigmaThreshold_	= std::sqrt(2.0 * t);
		else PeakRes.JF_UsedSigmaThreshold_ = 0.0;
	}
	else
	{
		PeakRes.JF_UsedSigmaThreshold_	= Zeus::logINF;	
	}

	PeakRes.JF_lnRho_				= PeakRes.JF_lnEvidence_ - PeakRes.JF_lnModelRatio_;
	if(PeakRes.JF_lnRho_ > THRESHOLD10PERC)
	{
		++DetectedObj_;
	}

	PeakRes.PK_BayesDetectStat_		= ((!EvStat) ? ((PeakRes.JF_lnRho_>PeakRes.JF_lnRhoTh_)?Zeus::PeakType::PK_DET_OK : Zeus::PeakType::PK_DET_NODET) : Zeus::PeakType::PK_DET_BAYNOCONV);

	if(GlbVars_.NonBlindDetection_)
	{
		PeakRes.JF_lnModelRatio_ = PeakRes.PeakNonBlindInfo_.SNR_;
	}


	PRINTINTOBUFFERFUNCT
		(buffer,BUFFERMAXCHAR,
		GLOBAL_FORMATSAMPLDETECT,
		PeakRes.JF_lnRho_,PeakRes.JF_lnEvidenceErrBar_,NestSampler_->GetNLikeEval()
		);

	(Zeus::ConManager::Instance())->PrintStr2Console(buffer);

	return PeakRes.PK_BayesDetectStat_;
}

void	PatchProcessor::InjectNonBlindParams2Peak(Zeus::PeakType &PeakResult)
{
	PeakResult.PeakNonBlindInfo_.CollLstIndex_	= PatchesGeoInfo_.Storage_[PatchNumber_].BPixel_;
//	PeakResult.PeakNonBlindInfo_.SrcIndex_		= PatchesGeoInfo_.Storage_[PatchNumber_].PatchNumber_;
	PeakResult.PeakNonBlindInfo_.SrcIndex_		= PatchesGeoInfo_.Storage_[PatchNumber_].SrcIndex_;
	PeakResult.PeakNonBlindInfo_.ErrFlux_		= PatchesGeoInfo_.Storage_[PatchNumber_].ErrFlux_;
	PeakResult.PeakNonBlindInfo_.ErrPos_		= PatchesGeoInfo_.Storage_[PatchNumber_].ErrPos_;
	PeakResult.PeakNonBlindInfo_.ErrRadius_		= PatchesGeoInfo_.Storage_[PatchNumber_].ErrRadius_;
	PeakResult.PeakNonBlindInfo_.PredFlux_		= PatchesGeoInfo_.Storage_[PatchNumber_].PredFlux_;
	PeakResult.PeakNonBlindInfo_.PredRadius_	= PatchesGeoInfo_.Storage_[PatchNumber_].PredRadius_;
	PeakResult.PeakNonBlindInfo_.SNR_			= PatchesGeoInfo_.Storage_[PatchNumber_].SNR_;
	PeakResult.PeakNonBlindInfo_.ErrRadiusHigh_	= PatchesGeoInfo_.Storage_[PatchNumber_].ErrRadiusHigh_;
	PeakResult.PeakNonBlindInfo_.ErrRadiusLow_	= PatchesGeoInfo_.Storage_[PatchNumber_].ErrRadiusLow_;
	PeakResult.PeakNonBlindInfo_.ErrFluxHigh_	= PatchesGeoInfo_.Storage_[PatchNumber_].ErrFluxHigh_;
	PeakResult.PeakNonBlindInfo_.ErrFluxLow_	= PatchesGeoInfo_.Storage_[PatchNumber_].ErrFluxLow_;
	PeakResult.PeakNonBlindInfo_.FluxBay_		= PatchesGeoInfo_.Storage_[PatchNumber_].FluxBay_;
	PeakResult.PeakNonBlindInfo_.RadiusBay_		= PatchesGeoInfo_.Storage_[PatchNumber_].RadiusBay_;
	PeakResult.PeakNonBlindInfo_.X0Y0Ptg_		= PatchesGeoInfo_.Storage_[PatchNumber_].X0Y0Ptg_;
	PeakResult.PeakNonBlindInfo_.SrcXCoord_		= PatchesGeoInfo_.Storage_[PatchNumber_].SrcXCoord_;
	PeakResult.PeakNonBlindInfo_.SrcYCoord_		= PatchesGeoInfo_.Storage_[PatchNumber_].SrcYCoord_;
	PeakResult.PeakNonBlindInfo_.SourcePtg_		= (((PeakResult.PeakNonBlindInfo_.SrcXCoord_ >= 0) && (PeakResult.PeakNonBlindInfo_.SrcYCoord_ >= 0))?PatchesGeoInfo_.Storage_[PatchNumber_].SourcePtg_:PatchesGeoInfo_.Storage_[PatchNumber_].X0Y0Ptg_);
	PeakResult.PeakNonBlindInfo_.QAIN_SrcPtg_	= PatchesGeoInfo_.Storage_[PatchNumber_].QAIN_SrcPtg_;
	PeakResult.PeakNonBlindInfo_.QAIN_CyR500	= PatchesGeoInfo_.Storage_[PatchNumber_].QAIN_CyR500;
	PeakResult.PeakNonBlindInfo_.QAIN_T500		= PatchesGeoInfo_.Storage_[PatchNumber_].QAIN_T500;
	PeakResult.PeakNonBlindInfo_.QAIN_detectable= PatchesGeoInfo_.Storage_[PatchNumber_].QAIN_detectable;
}

PP_FindOptParamsResType	PatchProcessor::CharacterAndAssess(int DoDetection)
{
	if(MaximaColl_.empty()) return PP_FOP_END;

	Zeus::PeakType								tPeak;
	MaxLikeSurfMaximaCollType::iterator			piv(MaximaColl_.begin());
	MaxLikeSurfMaximaCollType::const_iterator	const pivEnd(MaximaColl_.end());
	int											NLoopsDebug(0);
	IntCatalFilterFunctor						InCoreFunctor(GlbVars_.PatchBorder_,GlbVars_.PatchSz_);
	wchar_t	buffer[BUFFERMAXCHAR];

	OddsEval_->ClearCache();

	for(;piv!=pivEnd;++piv,++NLoopsDebug)
	{
		tPeak.PK_BayesDetectStat_	= Zeus::PeakType::PK_DET_NODET;

		if(piv->NonBlindTarget_ >= 0)
		{
			tPeak.PK_NP_DetectStat_		= Zeus::PeakType::PK_DET_NODET;
			const bool Find1PeakRes(Find1PeakOptimalParams(*piv,tPeak));
			OddsEval_->ClearCache();
			if(!Find1PeakRes)
				continue;
			if(IsOutsideDetectArea(tPeak))
				continue;
			NP_Assessment(tPeak,1);
			if(piv->NonBlindTarget_ == 1)
			{
				tPeak.PK_NP_DetectStat_ = Zeus::PeakType::PK_DET_OK;
			}
		}
		else
		{
			double	PredRad1(PatchesGeoInfo_.Storage_[PatchNumber_].PredRadius_);
			int		sign;

			if(PredRad1 < 0.0)
				errPriorValues(L"PatchProcessor::CharacterAndAssess",PatchesGeoInfo_.Storage_[PatchNumber_].SrcIndex_,PredRad1);

			tPeak.RealParams_			= Zeus::ObjFilterRealParams(PredRad1);

			if(GlbVars_.ScalesFillerType_ == 2) // If hardcoded scales
			{
//				tPeak.RealParams_			= Zone_->GetNearestScale(tPeak.RealParams_);
				tPeak.RealParams_			= Zone_->GetNearestScaleButFirst(tPeak.RealParams_);			
			}

			tPeak.UnTransParams_		= Zone_->TranslateParamsInverseInterpol(tPeak.RealParams_);

			tPeak.InitSrcFluxEst_		= piv->flux_;
			tPeak.InitCorrelation_		= piv->Correlation_;
			tPeak.GaussianIndex_		= piv->GaussianIndex_;
			tPeak.Pos_.XCoord_			= static_cast<double>((PatchesGeoInfo_.Storage_[PatchNumber_].SrcXCoord_ >= 0)?PatchesGeoInfo_.Storage_[PatchNumber_].SrcXCoord_:(GlbVars_.PatchSz_ >> 1));
			tPeak.Pos_.YCoord_			= static_cast<double>((PatchesGeoInfo_.Storage_[PatchNumber_].SrcYCoord_ >= 0)?PatchesGeoInfo_.Storage_[PatchNumber_].SrcYCoord_:(GlbVars_.PatchSz_ >> 1));
			tPeak.Pos_.YPix_			= Zeus::toInt(tPeak.Pos_.YCoord_ + 0.5);
			tPeak.Pos_.XPix_			= Zeus::toInt(tPeak.Pos_.XCoord_ + 0.5);
			tPeak.Odds_ = GetOddsResultRawOptimAmpRange(tPeak.Pos_.YPix_, tPeak.Pos_.XPix_, 2, tPeak.RealParams_, tPeak.ISNR2_, sign);
			tPeak.SrcAmplNormalised_	= std::sqrt(2.0 * tPeak.Odds_);
			tPeak.DetectionSigma_		= tPeak.SrcAmplNormalised_ / piv->GaussianIndex_;
			tPeak.SrcFlux_				= (tPeak.SrcAmplNormalised_ / std::sqrt(tPeak.ISNR2_)) * static_cast<double>(sign);
			tPeak.UsedSigmaThreshold_	= GlbVars_.NP_Sigma_;
			tPeak.PK_NP_DetectStat_		= Zeus::PeakType::PK_DET_OK;
		}

		if(tPeak.PK_NP_DetectStat_==Zeus::PeakType::PK_DET_OK)
		{
			if(CurrentIter_ == (GlbVars_.TwoStepsDetection_ - 1))
			{
				PRINTINTOBUFFERFUNCT
					(buffer,BUFFERMAXCHAR,
					GLOBAL_FORMATNPSTAT,
					tPeak.SrcAmplNormalised_,tPeak.Pos_.XCoord_,tPeak.Pos_.YCoord_
					);
				(Zeus::ConManager::Instance())->PrintStr2Console(buffer);

				if(!InCoreFunctor(tPeak))
				{
					int stat(Zeus::PeakType::PK_DET_ERROR);
					if(!(GlbVars_.NonBlindDetection_) || (piv->NonBlindTarget_!=0))
					{
						if(GlbVars_.NonBlindDetection_  && DoDetection)
						{
							InjectNonBlindParams2Peak(tPeak);
							if(!(GlbVars_.AssessmentKind_))
							{
								NonBlindMaxLikeReevaluation(*piv,tPeak);
								NP_Assessment(tPeak);
							}
						}
						if(GlbVars_.AssessmentKind_ && DoDetection)
						{
							(Zeus::ConManager::Instance())->PrintStr2Console(GLOBAL_STARTSAMPLING);
							stat = BayesLaplaceAssess(tPeak);

							if((stat == Zeus::PeakType::PK_DET_OK)		||
								(stat == Zeus::PeakType::PK_DET_NODET)	||
								(stat == Zeus::PeakType::PK_DET_BAYNOCONV))
							{
								JF_EvalParams(tPeak);
							}

							NestSampler_->ClearSamples();
						}
						if((!(GlbVars_.AssessmentKind_))  || (stat == Zeus::PeakType::PK_DET_ERROR))
						{
							NP_Eval_EvalParams(tPeak,stat);
						}

						if((tPeak.RealParams_.RealScale_ > 0.001) &&  (stat != Zeus::PeakType::PK_DET_ERROR))
						{
							tPeak.y0_ = y0Cylindircal(std::sqrt((GlbVars_.PixSz_ * GlbVars_.PixSz_  * SR2ARCMIN2) / PI) / tPeak.RealParams_.RealScale_) / Evaly0_denominator_;
						}
						else
						{tPeak.y0_ = -1.0;}					}
					else
					{
						NP_Eval_EvalParams(tPeak,stat);
					}
				}
			}
			
			tPeak.PatchNumber_ = PatchNumber_;
			
			if(!SubtractMaskSource(tPeak))
			{
				tPeak.Masked_ = 1;
				++NBrightPsInPatch_;
				if(NBrightPsInPatch_ > MAXNOFBRIGHTPSPERPATCH)
					err(ERRCOD_PWS_PARATOOMANYPS,ERRMSG_PWS_PARATOOMANYPS,L"PatchProcessor::CharacterAndAssess");
				(Zeus::ConManager::Instance())->PrintStr2Console(VERYBRIGHTOBJFOUND);
				PeakColl_.push_back(tPeak);
				return PP_FOP_BRGHTSRC;
			}
			else
			{
				tPeak.Masked_ = 0;
				PeakColl_.push_back(tPeak);
			}
		}
	}
	return PP_FOP_END;
}
//
int		PatchProcessor::SetupPriors(Zeus::PriorsTypeCollType&	tPriType, int PriorsSel)
{
	tPriType.clear();
	int	PriorType(PriorsSel % 8);
	
	
	switch(PriorType)
	{
	case 1:
		tPriType.push_back(Zeus::Priors::Gaussian);
		tPriType.push_back(Zeus::Priors::Gaussian);
		tPriType.push_back(GlbVars_.SZ_?Zeus::Priors::RadiusNinf:Zeus::Priors::Exponential);
		tPriType.push_back(Zeus::Priors::Uniform);
		break;
	case 2:
		tPriType.push_back(Zeus::Priors::Uniform);
		tPriType.push_back(Zeus::Priors::Uniform);
		tPriType.push_back(Zeus::Priors::Gaussian);
		tPriType.push_back(Zeus::Priors::Uniform);
		break;
	case 3:
		tPriType.push_back(Zeus::Priors::Gaussian);
		tPriType.push_back(Zeus::Priors::Gaussian);
		tPriType.push_back(Zeus::Priors::Gaussian);
		tPriType.push_back(Zeus::Priors::Uniform);
		break;
	case 4:
		tPriType.push_back(Zeus::Priors::Uniform);
		tPriType.push_back(Zeus::Priors::Uniform);
		tPriType.push_back(GlbVars_.SZ_?Zeus::Priors::RadiusNinf:Zeus::Priors::Exponential);
		tPriType.push_back(Zeus::Priors::Gaussian);
	case 5:
		tPriType.push_back(Zeus::Priors::Gaussian);
		tPriType.push_back(Zeus::Priors::Gaussian);
		tPriType.push_back(GlbVars_.SZ_?Zeus::Priors::RadiusNinf:Zeus::Priors::Exponential);
		tPriType.push_back(Zeus::Priors::Gaussian);
		break;
	case 6:
		tPriType.push_back(Zeus::Priors::Uniform);
		tPriType.push_back(Zeus::Priors::Uniform);
		tPriType.push_back(Zeus::Priors::Gaussian);
		tPriType.push_back(Zeus::Priors::Gaussian);
		break;
	case 7:
		tPriType.push_back(Zeus::Priors::Gaussian);
		tPriType.push_back(Zeus::Priors::Gaussian);
		tPriType.push_back(Zeus::Priors::Gaussian);
		tPriType.push_back(Zeus::Priors::Gaussian);
		break;
	default:
		tPriType.push_back(Zeus::Priors::Uniform);
		tPriType.push_back(Zeus::Priors::Uniform);
//
		if(GlbVars_.SZ_)
		{
//  PriorsSel & 0x08
// GlbVars_.Use2DFormula_ == 8
			if(!(PriorsSel & 0x08))
			{
				tPriType.push_back(Zeus::Priors::Marginal);
				tPriType.push_back(Zeus::Priors::Conditional);
//				tPriType.push_back(Zeus::Priors::Exponential);
//				tPriType.push_back(Zeus::Priors::Power);
//				tPriType.push_back(Zeus::Priors::RadiusNinf);
//				tPriType.push_back(Zeus::Priors::Uniform);
			}
			else
			{
				tPriType.push_back(Zeus::Priors::RadiusNinf);
				tPriType.push_back(Zeus::Priors::Uniform);
			}
		}
		else
		{
			tPriType.push_back(Zeus::Priors::Exponential);
			tPriType.push_back(Zeus::Priors::Power);
		}
	}
//
	return PriorType;
}

void	PatchProcessor::FindSources(int DoDetection)
{
	int		NIterations(GlbVars_.TwoStepsDetection_);
	double	LikeSurfISNR2Cte;
	int		StoreNoise(0);
	
	PeakColl_.clear();
	for(;CurrentIter_ < NIterations;++CurrentIter_)
	{
		LikeSurfISNR2Cte = 0.0;
		NBrightPsInPatch_ = 0;
		MakeWhiteningStructs(LikeSurfISNR2Cte);
		if(GlbVars_.AssessmentKind_ && DoDetection && (CurrentIter_ == (NIterations - 1)))
		{
			int BeamFullWith,BeamHalfWith;
			EvalEquivBeamSigmaPix(BeamFullWith,BeamHalfWith);
			Eval_ExpectedBackgPeaks();
			DetectedObj_ = 0;
			delete	NestSampler_;
			Zeus::PriorsTypeCollType	tPriType;
			SetupPriors(tPriType,GlbVars_.NonBlindDetection_);
			NestSampler_	= new NestSamplPwSObserParamNew(tPriType,*OddsEval_,GlbVars_.MuNe_NLivePoints_,GlbVars_.MuNe_NIndivSamples_,
				NMMAXSTEPS,GlbVars_.MuNe_FractTolEv_,GlbVars_.MuNe_xtraEnlgFactor_,BeamFullWith,BeamHalfWith,
				NMSEED);
		}
//
		if(GlbVars_.OutputUnits_ && (CurrentIter_ == (GlbVars_.TwoStepsDetection_ - 1)) && (!StoreNoise))
		{
			double LikeSurfSigma(0.0);
			
			(Zeus::ConManager::Instance())->PrintStr2Console(L"Evaluating SNR ... \n");

			ReadObsMaps(GlbVars_.N_PriorPlanes_?Zeus::MAPTYPE_BACKGROUND:Zeus::MAPTYPE_OBSERVATION);
			ComputeStats();
			MaskMaps(GlbVars_.SZ_?MASK_ALL:MASK_BORDER);
			MakeMainSurfaces(LikeSurfISNR2Cte,LikeSurfSigma,0);
			GetPatchNoiseProperties();

			(Zeus::ConManager::Instance())->PrintStr2Console(L"Evaluating SNR. Done \n");
			StoreNoise = 0;
		}
		else{StoreNoise = 1;}
//
		ReadObsMaps((GlbVars_.N_PriorPlanes_ && (CurrentIter_ != (GlbVars_.TwoStepsDetection_ - 1)))?Zeus::MAPTYPE_BACKGROUND:Zeus::MAPTYPE_OBSERVATION);
		ComputeStats();
		MaskMaps((GlbVars_.SZ_)?MASK_ALL:MASK_BORDER);

		PeakColl_.clear();
		do
		{
			double LikeSurfSigma(0.0);

			MakeMainSurfaces(LikeSurfISNR2Cte,LikeSurfSigma,1);
			MaximaColl_.clear();
			if((CurrentIter_ == (NIterations - 1)) && (!DoDetection))
			{return;}
			FindSourcePositions(LikeSurfISNR2Cte,StoreNoise);

			if(CharacterAndAssess(DoDetection) == PP_FOP_END)
			{
				break;
			}
			else
			{
				StoreNoise = 1;
			}
		}while(true);
		if(PeakColl_.empty() && (!(GlbVars_.NonBlindDetection_)))
		{
			break;
		}
	}

	if(GlbVars_.NonBlindDetection_)
	{
		RemoveSrcButCentral();
	}
	else
	{
		RemoveSrcGreyArea();
	}
}

void	PatchProcessor::RemoveSrcButCentral(void)
{
	Zeus::PeakCollType::const_iterator	pivCurrMax(PeakColl_.end());
	Zeus::PeakCollType::const_iterator	piv(PeakColl_.begin());
	Zeus::PeakCollType::const_iterator	const end(PeakColl_.end());
	const int	XCentralPix((PatchesGeoInfo_.Storage_[PatchNumber_].SrcXCoord_ >= 0)?PatchesGeoInfo_.Storage_[PatchNumber_].SrcXCoord_:(GlbVars_.PatchSz_ >> 1));
	const int	YCentralPix((PatchesGeoInfo_.Storage_[PatchNumber_].SrcYCoord_ >= 0)?PatchesGeoInfo_.Storage_[PatchNumber_].SrcYCoord_:(GlbVars_.PatchSz_ >> 1));

	int			PixCoordX;
	int			PixCoordY;
	int			LimitPixels(GlbVars_.NonBlindDetection_ >> 4);
	double		DetLevel,MaxDetLevel(-1.0e23);

	for(;piv != end;++piv)
	{
/*
		if(GlbVars_.AssessmentKind_)
		{
			PixCoordX	= (GlbVars_.Jf_Estimator_?piv->Pos_.JF_XPix_.Mean_:piv->Pos_.JF_XPix_.Mode_);
			PixCoordY	= (GlbVars_.Jf_Estimator_?piv->Pos_.JF_YPix_.Mean_:piv->Pos_.JF_YPix_.Mode_);
		}
		else
		{
			PixCoordX = piv->Pos_.XPix_;
			PixCoordY = piv->Pos_.YPix_;
		}
*/

		PixCoordX = piv->Pos_.XPix_;
		PixCoordY = piv->Pos_.YPix_;

		if((std::abs(PixCoordY - YCentralPix) <= LimitPixels)	 &&
			(std::abs(PixCoordX - XCentralPix) <= LimitPixels) 
			)
		{
			if(GlbVars_.AssessmentKind_)
			{
				DetLevel	= piv->JF_lnRho_;
			}
			else
			{
				DetLevel	= piv->SrcAmplNormalised_;			
			}
			if(DetLevel > MaxDetLevel)
			{
				MaxDetLevel = DetLevel;
				pivCurrMax	= piv;
			}
		}
	}

	if(pivCurrMax != PeakColl_.end())
	{
		Zeus::PeakType	t(*pivCurrMax);
		PeakColl_.clear();
		PeakColl_.push_back(t);
	}
	else
	{
		PeakColl_.clear();
	}
}

void	PatchProcessor::MaskSource(double flux,double YPos, double XPos,const Zeus::ObjFilterRealParams& ObjParams,double& CentralPixAmpl)
{
	int		YPix(Zeus::toInt(YPos + 0.5)),XPix(Zeus::toInt(XPos + 0.5));
	double	DeltaY,DeltaX;

 	if(GlbVars_.NotAlignedObjs_)
	{
		DeltaY	= (YPos - static_cast<double>(YPix));
		DeltaX  = (XPos - static_cast<double>(XPix));
	}
	else
	{
		DeltaY	= 0.0;
		DeltaX  = 0.0;
	}
	SrcInfoType		ObjSurf(Zone_->GetSrcObjRaw(ObjParams));
	MapsMainSurfFourType tempFour(Real2Fourier(ObjSurf.surf_));
	// Antenna for the 0 map
	// In SZ we don't mask
	MapsMainSurfFourType AntennaSurf(Zone_->GetAntenna(0));
	MultiplyInPlace(tempFour,AntennaSurf,Zeus::UB_NOUSE);
	tempFour.ShiftInPlace(DeltaY,DeltaX);
	RealPlaneSurfType tSrc(Fourier2Real(tempFour));
	Zeus::MultFunctor<double>	DummyGcc(1.0/static_cast<double>(GlbVars_.PatchSz_ * GlbVars_.PatchSz_));
	tSrc.Transform(DummyGcc,0,Zeus::UB_NOUSE);

	int Metric,YSz;
	(ObsPriorMapsColl_[0]).ws_.GetSz(YSz,Metric);

	MaskBordersCollType	MaskBordersXColl;
	MaskBordersCollType	MaskBordersYColl;
	CentralPixAmpl = *(tSrc.GetInnerData().begin());
	// Get X mask
	SrcMasker	DummyGcc1(flux,*(tSrc.GetInnerData().begin()),Metric,MaskBordersXColl,0);
	tSrc.Transform(DummyGcc1,0,Zeus::UB_NOUSE);
	// Get Y mask
	SrcMasker	DummyGcc2(flux,*(tSrc.GetInnerData().begin()),Metric,MaskBordersYColl,1);
	tSrc.Transform(DummyGcc2,1,Zeus::UB_NOUSE);
	// Get X pixels
	GetMaskBorderValues(MaskBordersXColl,(ObsPriorMapsColl_[0]).ws_,YPix,XPix,0);
	// Get Y pixels
	GetMaskBorderValues(MaskBordersYColl,(ObsPriorMapsColl_[0]).ws_,YPix,XPix,1);

	RealPlaneSurfType	SurfaceMask(Metric,false);
	Zeus::SetValueFunctor<double>	DummyGcc3(MASKVALUE);
	SurfaceMask.Transform(DummyGcc3,0,Zeus::UB_NOUSE);
	// Render X mask
	RenderSurfaceMask(MaskBordersXColl,SurfaceMask,Metric,0);
	// Render Y mask
	RenderSurfaceMask(MaskBordersYColl,SurfaceMask,Metric,1);

	static  int DebugCounter(1);
//	char	Debugbuffer[1024];

	//sprintf(Debugbuffer,"MaskSourceBefore_%d_",DebugCounter);
	//Zeus::DumpInOut_2d(Debugbuffer,256,256,256,0,(ObsPriorMapsColl_[0]).ws_.GetInnerData().begin(),1.0);

	PixMasker	DummyGcc4;
	Zeus::Transform(DummyGcc4,(ObsPriorMapsColl_[0]).ws_,YPix,XPix,SurfaceMask,false,Zeus::UB_NOUSE);

	//sprintf(Debugbuffer,"MaskSourceAfter_%d_",DebugCounter);
	//Zeus::DumpInOut_2d(Debugbuffer,256,256,256,0,(ObsPriorMapsColl_[0]).ws_.GetInnerData().begin(),1.0);
}

void	PatchProcessor::GetMaskBorderValues(MaskBordersCollType& MskBColl,const RealPlaneSurfType& surf,int YPix,int XPix,int Swap)
{
	MaskBordersCollType::iterator	piv(MskBColl.begin());
	MaskBordersCollType::iterator	const end(MskBColl.end());
	RealPlaneSurfType::DataInnerType::const_iterator const pivOrg(surf.GetInnerData().begin());
	int	offset,offsetAux,Metric,YSz;

	surf.GetSz(YSz,Metric);
	for(;piv != end;++piv)
	{
		piv->Out_		= 0;
		offset			= YPix + piv->y_ ;
		if((offset < 0) || ((offset >= Metric)))
		{
			piv->PixVal_ = 0.0;
			if(!Swap)
			{
				piv->Out_	 = 1;
			}
			continue;
		}
		offset		*= Metric;
		offsetAux	= XPix + piv->x_;
		if((offsetAux < 0) || ((offsetAux >= Metric)))
		{
			piv->PixVal_ = 0.0;
			if(Swap)
			{
				piv->Out_	 = 1;
			}
			continue;
		}
		offset		+= offsetAux;
		piv->PixVal_ = *(pivOrg + offset);
	}

	std::sort(MskBColl.begin(),MskBColl.end(),Swap?SortMaskBordersY:SortMaskBordersX);
}

void	PatchProcessor::RenderSurfaceMask(const MaskBordersCollType& MskBColl,RealPlaneSurfType& tSurf,int Metric,int Swap)
{
	double	slope;
	int		offset;

	RealPlaneSurfType::DataInnerType::iterator const pivOrg(tSurf.GetInnerData().begin());
	MaskBordersCollType::const_iterator	piv(MskBColl.begin());
	MaskBordersCollType::const_iterator	pivAux,pivAhead;
	double	NewValue,OldValue;
	int		TotOffset;

	MaskBordersCollType::const_iterator	const end(MskBColl.end());
	while(piv != end)
	{
		for(pivAux = piv;(pivAux != end) && (Swap?(pivAux->x_ == piv->x_):(pivAux->y_ == piv->y_));++pivAux)
			;
		pivAhead = pivAux - 1;
		if(piv->Out_ != 1)
		{
			slope = (pivAhead->PixVal_ - piv->PixVal_);
			if(Swap)
			{
				offset = (piv->x_ < 0? piv->x_ + Metric : piv->x_);
				slope /= (static_cast<double>(pivAhead->y_) - static_cast<double>(piv->y_));
			}
			else
			{
				offset = (Metric * (piv->y_ < 0? piv->y_ + Metric : piv->y_));
				slope /= (static_cast<double>(pivAhead->x_) - static_cast<double>(piv->x_));
			}
			for(int i = (Swap?piv->y_:piv->x_);i <= (Swap?pivAhead->y_:pivAhead->x_);++i)
			{
				if(Swap)
				{
					TotOffset = offset + Metric * (i < 0? i + Metric : i);
				}
				else
				{
					TotOffset = offset + (i < 0? i + Metric : i);
				}

				OldValue  = *(pivOrg + TotOffset);
				NewValue  = slope * (Swap?static_cast<double>(i - piv->y_):static_cast<double>(i - piv->x_)) + piv->PixVal_;
				if(OldValue == MASKVALUE)
				{
					*(pivOrg + TotOffset) = NewValue;
				}
				else
				{
					*(pivOrg + TotOffset) = (OldValue + NewValue) / 2.0;
				}
			}
		}
		piv = pivAux;
	}
}
//
int		PatchProcessor::ReportCurrentSrcs(int SrcIdex)
{
	if(GlbVars_.NonBlindDetection_)
	{
		if(PeakColl_.empty())
		{
			(PlanckInfo::Instance())->AppendNonValidObjs2File(SrcIdex,PatchNumber_);
			return 1;
		}
	}

	SetMultiScaleY(PeakColl_);

	if(GlbVars_.OutputUnits_)
	{
		SwapPatchNoiseEstimates();
	}

	return (PlanckInfo::Instance())->AppendDetectObjs2File(PatchNumber_,PeakColl_);
}
//
int		PatchProcessor::SetMultiScaleY(Zeus::PeakCollType& PeakColl)
{
	if(PeakColl.empty()) return 0;

	unsigned int NMaxScales((unsigned int) GlbVars_.ScalesMaxNumber_);
	Zeus::ScaleLikeNoise	t;

	Zeus::PeakCollType::iterator				pivPeak(PeakColl.begin());
	Zeus::PeakCollType::iterator				const endPeak(PeakColl.end());
	ScalesPropsCollType::const_iterator			const end(ScalesProps_.end());

	ScalesPropsCollType::const_iterator			piv;

	for(;pivPeak != endPeak;++pivPeak)
	{
		(pivPeak->ScaleLikeNoise_).clear();
		piv	=	ScalesProps_.begin();

		for(unsigned int i=0;(piv!=end) && (i<NMaxScales);++piv,++i)
		{
			const double Sth2(piv->Th1OverSigma_*piv->Th1OverSigma_);
			const double sigmaTh(1.0/piv->Th1OverSigma_);
			const double sigmaPatch(piv->LikelihoodSigma_/Sth2);
			
			t.Scale_	= piv->SrcScale_;
			t.Like_		= piv->LikeCentralPix_ / Sth2;
//			t.Noise_	= ((sigmaTh>sigmaPatch)?sigmaTh:sigmaPatch);
			t.Noise_	= sigmaTh;
			(pivPeak->ScaleLikeNoise_).push_back(t);
		}
	}
	return ((ScalesProps_.size() > NMaxScales)?NMaxScales:ScalesProps_.size());
}
//
void	PatchProcessor::SwapPatchNoiseEstimates(void)
{

	if(ScalesProps_.empty()) return;

	std::sort(ScalesProps_.begin(),ScalesProps_.end());

	ScalesPropsCollType::value_type		t;
	double								tSigma;
	Zeus::PeakCollType::iterator		piv(PeakColl_.begin());
	Zeus::PeakCollType::const_iterator	const end(PeakColl_.end());

	for(;piv != end;++piv)
	{
		t.SrcScale_ = piv->RealParams_.RealScale_;
		ScalesPropsCollType::iterator	uBound(std::upper_bound(ScalesProps_.begin(),ScalesProps_.end(),t));

		if(uBound == ScalesProps_.begin())
		{
			const double thSigma		= 1.0 / uBound->Th1OverSigma_;
			const double ptchSigma		= uBound->LikelihoodSigma_ / (uBound->Th1OverSigma_ * uBound->Th1OverSigma_);
//			tSigma = ((ptchSigma < thSigma) ? ptchSigma : thSigma);

			tSigma = thSigma;

		}
		else if(uBound == ScalesProps_.end())
		{
			const double thSigma		= 1.0 / (uBound-1)->Th1OverSigma_;
			const double ptchSigma		= (uBound-1)->LikelihoodSigma_ / ((uBound-1)->Th1OverSigma_ * (uBound-1)->Th1OverSigma_);
//			tSigma = ((ptchSigma < thSigma) ? ptchSigma : thSigma); 

			tSigma = thSigma;
		}
		else if(std::abs(uBound->SrcScale_ - (uBound-1)->SrcScale_)< 0.1)
		{
			const double thSigma		= 1.0 / uBound->Th1OverSigma_;
			const double ptchSigma		= uBound->LikelihoodSigma_ / (uBound->Th1OverSigma_ * uBound->Th1OverSigma_);
			const double thSigma1		= 1.0 / (uBound-1)->Th1OverSigma_;
			const double ptchSigma1		= (uBound-1)->LikelihoodSigma_ / ((uBound-1)->Th1OverSigma_ * (uBound-1)->Th1OverSigma_);
//			const double s0((thSigma>ptchSigma)?thSigma:ptchSigma);
//			const double s1((thSigma1>ptchSigma1)?thSigma1:ptchSigma1);
			const double s0(thSigma);
			const double s1(thSigma1);

			tSigma = ((s1>s0)?s1:s0);
		}
		else
		{
			const double thSigma		= 1.0 / uBound->Th1OverSigma_;
			const double ptchSigma		= uBound->LikelihoodSigma_ / (uBound->Th1OverSigma_ * uBound->Th1OverSigma_);
			const double thSigma1		= 1.0 / (uBound-1)->Th1OverSigma_;
			const double ptchSigma1		= (uBound-1)->LikelihoodSigma_ / ((uBound-1)->Th1OverSigma_ * (uBound-1)->Th1OverSigma_);
//			const double s0((thSigma>ptchSigma)?thSigma:ptchSigma);
//			const double s1((thSigma1>ptchSigma1)?thSigma1:ptchSigma1);
			const double s0(thSigma);
			const double s1(thSigma1);

			tSigma = ((s0-s1)/((uBound)->SrcScale_ - (uBound-1)->SrcScale_)) * (t.SrcScale_- (uBound-1)->SrcScale_) + s1;
		}

		piv->DetectionSigma_ = piv->SrcAmplNormalised_ = (piv->SrcFlux_ / tSigma);
	}
}
//
void	PatchProcessor::JF_EvalParams(Zeus::PeakType & PeakRes)
{
	wchar_t						buffer[BUFFERMAXCHAR];
	ParamsStatsCollType			Params;
	LinearCorrParamsType		LineParams;

	NestSampler_->GetParamsStatsEx(Params);

	NestSampler_->GetLinearCorrelatioParams((Params[3].StDev_*Params[3].StDev_)/(Params[2].StDev_*Params[2].StDev_),LineParams);

	PeakRes.Pos_.JF_XCoord_.Mode_		= Params[0].Mode_;
	PeakRes.Pos_.JF_YCoord_.Mode_		= Params[1].Mode_;
	PeakRes.Pos_.JF_XCoord_.Mean_		= Params[0].Mean_;
	PeakRes.Pos_.JF_YCoord_.Mean_		= Params[1].Mean_;
	PeakRes.Pos_.JF_XPix_.Mode_			= Zeus::toInt(PeakRes.Pos_.JF_XCoord_.Mode_ + 0.5);
	PeakRes.Pos_.JF_YPix_.Mode_			= Zeus::toInt(PeakRes.Pos_.JF_YCoord_.Mode_ + 0.5);
	PeakRes.Pos_.JF_XPix_.Mean_			= Zeus::toInt(PeakRes.Pos_.JF_XCoord_.Mean_ + 0.5);
	PeakRes.Pos_.JF_YPix_.Mean_			= Zeus::toInt(PeakRes.Pos_.JF_YCoord_.Mean_ + 0.5);
	PeakRes.JF_Radius_.Mode_			= Params[2].Mode_;
	PeakRes.JF_SrcFlux_.Mode_			= Params[3].Mode_;
	PeakRes.JF_Radius_.Mean_			= Params[2].Mean_;
	PeakRes.JF_SrcFlux_.Mean_			= Params[3].Mean_;
	double	FluxConvCte;
	if(GlbVars_.SZ_)
	{
		PeakRes.SrcFlux_mJys_			= -1.0;
		FluxConvCte	= (0.001 * GlbVars_.PixSz_ * GlbVars_.PixSz_ * SR2ARCMIN2)  * GlbVars_.ProfParam_.FluxCalibCte_;
		PeakRes.SrcCompt_arcmin2_		= (GlbVars_.Jf_Estimator_?PeakRes.JF_SrcFlux_.Mean_:PeakRes.JF_SrcFlux_.Mode_) * FluxConvCte;
	}
	else
	{
		PeakRes.SrcCompt_arcmin2_		= -1.0;
		FluxConvCte = (GlbVars_.PixSz_ * GlbVars_.PixSz_ * 1e9 * GlbVars_.ProfParam_.FluxCalibCte_);
		PeakRes.SrcFlux_mJys_			= (GlbVars_.Jf_Estimator_?PeakRes.JF_SrcFlux_.Mean_:PeakRes.JF_SrcFlux_.Mode_) * FluxConvCte;
	}
	// Error bar is 95% probability ~ 2 sigma for Gaussian distribution
	PeakRes.ErrorBars_.TotalPosErrorBar_ = 0.5 * ERRPOS95CTE * (Params[0].StDev_ + Params[1].StDev_) * GlbVars_.PixSz_ * RAD2ARCMIN;
	// rescaling to 1 sigma is the sense that 2 * value = Error bar for 95%
	PeakRes.ErrorBars_.TotalPosErrorBar_ /= 2.0;
	// Error bars are 1 sigma
	PeakRes.ErrorBars_.RadiusErrorBar_		=  Params[2].StDev_ * GlbVars_.ProfParam_.RadiusCalCte_;
	PeakRes.ErrorBars_.FluxErrorBar_		=  Params[3].StDev_ * FluxConvCte;
	PeakRes.ErrorBars_.LowFluxErrorBar_		=  (((GlbVars_.Jf_Estimator_?Params[3].Mean_:Params[3].Mode_) - Params[3].LowHPD_) * FluxConvCte)/2.0;
	PeakRes.ErrorBars_.HighFluxErrorBar_	=  ((Params[3].HighHPD_ - (GlbVars_.Jf_Estimator_?Params[3].Mean_:Params[3].Mode_)) * FluxConvCte)/2.0;
	PeakRes.ErrorBars_.LowRadiusErrorBar_	=  (((GlbVars_.Jf_Estimator_?Params[2].Mean_:Params[2].Mode_) - Params[2].LowHPD_) * GlbVars_.ProfParam_.RadiusCalCte_)/2.0;
	PeakRes.ErrorBars_.HighRadiusErrorBar_	=  ((Params[2].HighHPD_ - (GlbVars_.Jf_Estimator_?Params[2].Mean_:Params[2].Mode_)) * GlbVars_.ProfParam_.RadiusCalCte_)/2.0;
	
	PeakRes.ErrorBars_.DegenSlope_			= (LineParams.Slope_ * FluxConvCte)/GlbVars_.ProfParam_.RadiusCalCte_;
	PeakRes.ErrorBars_.DegenCorr_			= LineParams.Corr_;
	PeakRes.ErrorBars_.DegenOrd_			= LineParams.Ord_ * FluxConvCte;
	PeakRes.ErrorBars_.DegenOrdYErr_			= LineParams.ErrOrd_ * FluxConvCte;
	PeakRes.ErrorBars_.DegenSlopeYErr_		= (LineParams.ErrSlope_ *  FluxConvCte) / GlbVars_.ProfParam_.RadiusCalCte_;

	PRINTINTOBUFFERFUNCT
		(buffer,BUFFERMAXCHAR,
		GLOBAL_FORMATSAMPLPAR,
		L"x",Params[0].Mean_,Params[0].Mode_,Params[0].StDev_,Params[0].LowHPD_,Params[0].HighHPD_
		);
	(Zeus::ConManager::Instance())->PrintStr2Console(buffer);
	PRINTINTOBUFFERFUNCT
		(buffer,BUFFERMAXCHAR,
		GLOBAL_FORMATSAMPLPAR,
		L"y",Params[1].Mean_,Params[1].Mode_,Params[1].StDev_,Params[1].LowHPD_,Params[1].HighHPD_
		);
	(Zeus::ConManager::Instance())->PrintStr2Console(buffer);
	PRINTINTOBUFFERFUNCT
		(buffer,BUFFERMAXCHAR,
		GLOBAL_FORMATSAMPLPAR,
		L"Y",Params[3].Mean_ * FluxConvCte,Params[3].Mode_ * FluxConvCte,Params[3].StDev_ * FluxConvCte,Params[3].LowHPD_ * FluxConvCte,Params[3].HighHPD_ * FluxConvCte
		);
	(Zeus::ConManager::Instance())->PrintStr2Console(buffer);
	PRINTINTOBUFFERFUNCT
		(buffer,BUFFERMAXCHAR,
		GLOBAL_FORMATSAMPLPAR,
		L"T",Params[2].Mean_* GlbVars_.ProfParam_.RadiusCalCte_,Params[2].Mode_* GlbVars_.ProfParam_.RadiusCalCte_,Params[2].StDev_* GlbVars_.ProfParam_.RadiusCalCte_,Params[2].LowHPD_* GlbVars_.ProfParam_.RadiusCalCte_,Params[2].HighHPD_* GlbVars_.ProfParam_.RadiusCalCte_
		);
	(Zeus::ConManager::Instance())->PrintStr2Console(buffer);

	PRINTINTOBUFFERFUNCT
		(buffer,BUFFERMAXCHAR,
		GLOBAL_FORMATDEGENPAR,
		(PeakRes.ErrorBars_.DegenSlope_*1000.0),PeakRes.ErrorBars_.DegenCorr_,(PeakRes.ErrorBars_.DegenOrd_*1000.0),(PeakRes.ErrorBars_.DegenSlopeYErr_*1000.0),(PeakRes.ErrorBars_.DegenOrdYErr_*1000.0)
		);
	(Zeus::ConManager::Instance())->PrintStr2Console(buffer);

}

void	PatchProcessor::NP_Eval_EvalParams(Zeus::PeakType & PeakRes,int stat)
{
	const double	PixSq(static_cast<double>(GlbVars_.PixSz_ * GlbVars_.PixSz_));
	double			FluxCte;
		
	// sigma = (1.0 / std::sqrt(PeakRes.ISNR2_))
	// this needs to be changed; this formula is not exact
	// check PwSI
	PeakRes.ErrorBars_.FluxErrorBar_	= (1.0 / std::sqrt(PeakRes.ISNR2_));

	if(GlbVars_.SZ_)
	{
		FluxCte = (0.001 * PixSq * SR2ARCMIN2 * GlbVars_.ProfParam_.FluxCalibCte_);
		PeakRes.SrcFlux_mJys_			= -1.0;
		PeakRes.SrcCompt_arcmin2_		= PeakRes.SrcFlux_ * FluxCte;
	}
	else
	{
		FluxCte = (PixSq * 1.0e9 * GlbVars_.ProfParam_.FluxCalibCte_);
		PeakRes.SrcCompt_arcmin2_		= -1.0;
		PeakRes.SrcFlux_mJys_			= PeakRes.SrcFlux_ * FluxCte;
	}

	PeakRes.ErrorBars_.FluxErrorBar_ *= FluxCte;

	// No error bar for position
	PeakRes.ErrorBars_.TotalPosErrorBar_	= Zeus::logZERO;
	// No error bar for radius
	PeakRes.ErrorBars_.RadiusErrorBar_		= Zeus::logZERO;

	PeakRes.Pos_.JF_XCoord_.Mode_	= PeakRes.Pos_.JF_XCoord_.Mean_ = PeakRes.Pos_.XCoord_;
	PeakRes.Pos_.JF_YCoord_.Mode_	= PeakRes.Pos_.JF_YCoord_.Mean_ = PeakRes.Pos_.YCoord_;
	PeakRes.Pos_.JF_XPix_.Mode_		= PeakRes.Pos_.JF_XPix_.Mean_	= PeakRes.Pos_.XPix_;
	PeakRes.Pos_.JF_YPix_.Mode_		= PeakRes.Pos_.JF_YPix_.Mean_	= PeakRes.Pos_.YPix_;
	PeakRes.JF_Radius_.Mode_		= PeakRes.JF_Radius_.Mean_		= PeakRes.RealParams_.RealScale_;
	PeakRes.JF_SrcFlux_.Mode_		= PeakRes.JF_SrcFlux_.Mean_		= PeakRes.SrcFlux_;
	if(stat == Zeus::PeakType::PK_DET_ERROR)
	{
		PeakRes.JF_lnRhoTh_				= Zeus::logZERO;
		PeakRes.JF_lnRho_				= Zeus::logZERO;
		PeakRes.JF_lnEvidence_			= Zeus::logZERO;
		PeakRes.JF_lnEvidenceErrBar_	= Zeus::logZERO;
		PeakRes.JF_lnModelRatio_		= Zeus::logZERO;
		PeakRes.JF_lnFormFactor_		= Zeus::logZERO;
		PeakRes.JF_UsedSigmaThreshold_	= -1.0;
	}
	PeakRes.PK_BayesDetectStat_		= stat;
}


void	PatchProcessor::SubSrcFromMap(MaskMapType& Map)
{
	int Metric,YSz;
	Map.ws_.GetSz(YSz,Metric);

	MaskBordersCollType				MaskBordersXColl;
	MaskBordersCollType				MaskBordersYColl;

	Zeus::RealPlane<double>			tSurface(Metric,true);
	int								YPix,XPix;
	double							DeltaY,DeltaX;
	Zeus::PeakCollType::const_iterator	piv(PeakColl_.begin());
	Zeus::PeakCollType::const_iterator	const end(PeakColl_.end());
	const double					NormCte(static_cast<double>(GlbVars_.PatchSz_ * GlbVars_.PatchSz_));
	RealPlaneSurfType				SurfaceMask;

	//Zeus::DumpInOut_2d("SubSrcFromMapBeforeSubtraction",256,256,256,0,Map.ws_.GetInnerData().begin(),1.0);
/*
TODO
(GlbVars_.SZ_) &&
((DetectionSigma < 10) || ((std::abs((GlbVars_.SrcMaxScale_ - ObjParams.RealScale_)/ GlbVars_.SrcMaxScale_) < 0.02) && (DetectionSigma < 20)))

*/
	for(;piv != end;++piv)
	{
		if(
			(piv->Masked_) ||
			((GlbVars_.SZ_) && ((piv->SrcAmplNormalised_ < SIGMASUBTRACTBACK) || ((std::abs((GlbVars_.SrcMaxScale_ - piv->RealParams_.RealScale_)/ GlbVars_.SrcMaxScale_) < 0.02) && (piv->SrcAmplNormalised_ < 20.0))))
			)
			continue;

		YPix	= Zeus::toInt(piv->Pos_.YCoord_ + 0.5);
		XPix	= Zeus::toInt(piv->Pos_.XCoord_ + 0.5);
		if(GlbVars_.NotAlignedObjs_)
		{DeltaY	= (piv->Pos_.YCoord_ - static_cast<double>(YPix));DeltaX  = (piv->Pos_.XCoord_ - static_cast<double>(XPix));}
		else
		{DeltaY	= 0.0;DeltaX  = 0.0;}

		{
			Zeus::ObjFilterRealParams	UnTransParamsTemp(piv->UnTransParams_);
/*
			double CorrectionTemp(1.0);
			if(piv->SrcAmplNormalised_ > 10.0)
			{
				CorrectionTemp = ((piv->SrcAmplNormalised_ > 12.0)?0.9:1.0-0.05*(piv->SrcAmplNormalised_ - 10.0));
				UnTransParamsTemp.RealScale_ *= CorrectionTemp;
			}
*/
			SrcInfoType					ObjSurf(Zone_->GetSrcObjRaw(UnTransParamsTemp));
			MapsMainSurfFourType		tempFour(Real2Fourier(ObjSurf.surf_));
			int AntennaIndex(GlbVars_.GetObsFreqIndex(Map.Freq_));
			if(AntennaIndex < 0)
				err(ERRCOD_PWS_FREQNOTFOUND,ERRMSG_PWS_FREQNOTFOUND,L"PatchProcessor::SubSrcFromMap");
			MapsMainSurfFourType AntennaSurf(Zone_->GetAntenna(AntennaIndex));
			MultiplyInPlace(tempFour,AntennaSurf,Zeus::UB_NOUSE);
			tempFour.ShiftInPlace(DeltaY,DeltaX);
			RealPlaneSurfType tSrc(Fourier2Real(tempFour));
			double SrcFluxTemp(piv->SrcFlux_);
/*
			if(piv->SrcAmplNormalised_ > 10.0)
			{
				SrcFluxTemp *= (CorrectionTemp*CorrectionTemp);
			}
*/
			Zeus::MultFunctor<double>	DummyGcc1(SrcFluxTemp / NormCte);
			tSrc.Transform(DummyGcc1,0,Zeus::UB_NOUSE);
			SubAddFunctor	DummyGcc2(true);
			Zeus::Transform(DummyGcc2,tSurface,YPix,XPix,tSrc,false,Zeus::UB_NOUSE);
		}
	}

	Zeus::AddInPlace(Map.ws_,tSurface,Zeus::UB_NOUSE);

	//Zeus::DumpInOut_2d("SubSrcFromMapAfterSubtraction",256,256,256,0,Map.ws_.GetInnerData().begin(),1.0);

	piv = PeakColl_.begin();

	for(;piv != end;++piv)
	{
		if(!(piv->Masked_))
			continue;
		YPix	= Zeus::toInt(piv->Pos_.YCoord_ + 0.5);
		XPix	= Zeus::toInt(piv->Pos_.XCoord_ + 0.5);
		if(GlbVars_.NotAlignedObjs_)
		{
			DeltaY	= (piv->Pos_.YCoord_ - static_cast<double>(YPix));
			DeltaX  = (piv->Pos_.XCoord_ - static_cast<double>(XPix));
		}
		else
		{
			DeltaY	= 0.0;
			DeltaX  = 0.0;
		}

		{
			SrcInfoType				ObjSurf(Zone_->GetSrcObjRaw(piv->UnTransParams_));
			MapsMainSurfFourType	tempFour(Real2Fourier(ObjSurf.surf_));
			int AntennaIndex(GlbVars_.GetObsFreqIndex(Map.Freq_));
			MapsMainSurfFourType AntennaSurf(Zone_->GetAntenna(AntennaIndex));
			MultiplyInPlace(tempFour,AntennaSurf,Zeus::UB_NOUSE);
			tempFour.ShiftInPlace(DeltaY,DeltaX);
			RealPlaneSurfType tSrc(Fourier2Real(tempFour));
			{
				Zeus::MultFunctor<double>	DummyGcc1(1.0/NormCte);
				tSrc.Transform(DummyGcc1,0,Zeus::UB_NOUSE);
				// Get X mask
				SrcMasker	DummyGcc2(piv->SrcFlux_,*(tSrc.GetInnerData().begin()),GlbVars_.PatchSz_,MaskBordersXColl,0);
				tSrc.Transform(DummyGcc2,0,Zeus::UB_NOUSE);
				// Get Y mask
				SrcMasker	DummyGcc3(piv->SrcFlux_,*(tSrc.GetInnerData().begin()),GlbVars_.PatchSz_,MaskBordersYColl,1);
				tSrc.Transform(DummyGcc3,1,Zeus::UB_NOUSE);
				// Get X pixels
				GetMaskBorderValues(MaskBordersXColl,Map.ws_,YPix,XPix,0);
				// Get Y pixels
				GetMaskBorderValues(MaskBordersYColl,Map.ws_,YPix,XPix,1);

				SurfaceMask.MakeNewSz(GlbVars_.PatchSz_,false);

				Zeus::SetValueFunctor<double>	DummyGcc4(MASKVALUE);
				SurfaceMask.Transform(DummyGcc4,0,Zeus::UB_NOUSE);
				// Render X mask
				RenderSurfaceMask(MaskBordersXColl,SurfaceMask,GlbVars_.PatchSz_,0);
				// Render Y mask
				RenderSurfaceMask(MaskBordersYColl,SurfaceMask,GlbVars_.PatchSz_,1);
				PixMasker	DummyGcc5;
				Zeus::Transform(DummyGcc5,Map.ws_,YPix,XPix,SurfaceMask,false,Zeus::UB_NOUSE);
			}
		}
	}

	//Zeus::DumpInOut_2d("SubSrcFromMapMaskComplete",256,256,256,0,Map.ws_.GetInnerData().begin(),1.0);
}


void	PatchProcessor::NonBlindMaxLikeReevaluation(MaxLikeSurfMaximum& LikePeak,Zeus::PeakType &PeakResult)
{
	const int	PriorType(GlbVars_.NonBlindDetection_ % 4);
	const int	PriorSZ_NInfo((GlbVars_.NonBlindDetection_ & 0x08) != 0);
	const int	PriorPixTol(GlbVars_.NonBlindDetection_ >> 4);
	const int	XCentralPix((PatchesGeoInfo_.Storage_[PatchNumber_].SrcXCoord_ >= 0)?PatchesGeoInfo_.Storage_[PatchNumber_].SrcXCoord_:(GlbVars_.PatchSz_ >> 1));
	const int	YCentralPix((PatchesGeoInfo_.Storage_[PatchNumber_].SrcYCoord_ >= 0)?PatchesGeoInfo_.Storage_[PatchNumber_].SrcYCoord_:(GlbVars_.PatchSz_ >> 1));
	int			MaxX,MaxY;
	int			sign(0);
	double		Isnr;
	double		logFlux(Zeus::logZERO),logRadius(Zeus::logZERO);

	if(GlbVars_.SearchPos_ && (NSnakesDim_ > 3))
		errPriorConfig(L"PatchProcessor::NonBlindMaxLikeReevaluation",PeakResult.PeakNonBlindInfo_.SrcIndex_,ERRMSG_PWS_INVALPRIORCONG);

	switch(PriorType)
	{
	case	1:
	{
		MinArrayType		InitGuessFinal(NSnakesDim_);
		double				tOdds;

		InitGuessFinal[0] = static_cast<double>(YCentralPix);
		InitGuessFinal[1] = static_cast<double>(XCentralPix);
		InitGuessFinal[2] = static_cast<double>(LikePeak.FiltParams_.BinScale_);

		if(!(PwSnake_->PowellMinimunSearch(InitGuessFinal,PeakResult.Odds_)))
		{
			errPriorConfig(L"PatchProcessor::NonBlindMaxLikeReevaluation",PeakResult.PeakNonBlindInfo_.SrcIndex_,ERRMSG_PWS_INVALPRIORDETECT);
		}
		PeakResult.Odds_			= -PeakResult.Odds_;
		PeakResult.Pos_.YCoord_		= InitGuessFinal[0];
		PeakResult.Pos_.XCoord_		= InitGuessFinal[1];
		PeakResult.UnTransParams_	= Zeus::ObjFilterRealParams(InitGuessFinal[2]);
		PeakResult.RealParams_		= Zone_->TranslateParams(PeakResult.UnTransParams_);
		OddsEval_->GetOddsResultRawOptimAmp(YCentralPix, XCentralPix, ObjFilterParams(Zeus::toInt(InitGuessFinal[2])), PeakResult.ISNR2_, sign);
//		PeakResult.ISNR2_			= GetISNR(PeakResult.UnTransParams_);
	}
		break;
	case	2:
	{
		if(PeakResult.PeakNonBlindInfo_.PredRadius_ < 0.0)
			errPriorValues(L"PatchProcessor::NonBlindMaxLikeReevaluation",PeakResult.PeakNonBlindInfo_.SrcIndex_,PeakResult.PeakNonBlindInfo_.PredRadius_);

		double	tOdds, MaxOdds(-1.0e30), MaxIsnr;
		double	val[9];
		int		tsign(0);
		Zeus::LArr2D<double>	Buffer(((PriorPixTol<<1) + 1)*((PriorPixTol<<1) + 1),((PriorPixTol<<1) + 1));

		Zeus::ObjFilterRealParams in(PeakResult.PeakNonBlindInfo_.PredRadius_);

		for(int j = -PriorPixTol;j <= PriorPixTol;++j)
		{
			for(int i = -PriorPixTol;i <= PriorPixTol;++i)
			{
				tOdds = OddsEval_->GetOddsResultRawOptimAmp(YCentralPix + j,XCentralPix + i,in,Isnr,tsign);

				if(tOdds >= 0.0)
					Buffer(j+PriorPixTol,i+PriorPixTol) = tOdds;
				else
					Buffer(j+PriorPixTol,i+PriorPixTol) = 0.0;

				if(tOdds > MaxOdds)
				{
					MaxY = j+PriorPixTol;MaxX = i+PriorPixTol;
					MaxOdds = tOdds;
					MaxIsnr	= Isnr;
					sign = tsign;
				}
			}
		}

		for(int j=-1;j<=1;++j)
		{
			for(int i=-1;i<=1;++i)
			{
				val[(j+1)*3+(i+1)] = ((((MaxY + j)>=0) && ((MaxX + i)>=0))?Buffer(MaxY + j,MaxX + i):0.0);  
			}
		}

		Centroid(val);

		const	int tXPixTransl(XCentralPix - PriorPixTol);
		const	int tYPixTransl(YCentralPix - PriorPixTol);

		PeakResult.Odds_				= MaxOdds;
		PeakResult.ISNR2_				= MaxIsnr;
		PeakResult.Pos_.YCoord_			= static_cast<double>(MaxY + tYPixTransl) + val[1];
		PeakResult.Pos_.XCoord_			= static_cast<double>(MaxX + tXPixTransl) + val[0];
		PeakResult.RealParams_			= Zeus::ObjFilterRealParams(PeakResult.PeakNonBlindInfo_.PredRadius_);
		PeakResult.UnTransParams_		= Zone_->TranslateParamsInverseInterpol(PeakResult.RealParams_);
	
	}

		break;
	case	3:
		if(	(GlbVars_.SZ_ == 1) && 
			(GlbVars_.N_ObsPlanes_ != 1) &&
			(GlbVars_.Jf_Estimator_ == 0) &&
			(GlbVars_.AssessmentKind_ == 0)
		)
		{
			static int		FileID(0);
			double			FluxBestFit,RadiusBestFit,NSigmaLow,NSigmaHigh;
			Zeus::UniformPrior	*BetaPrior(NULL);
			double MaxBeta;
			double MinBeta;
			double	BetaStep;
			Zeus::GammaFuncts Gammafun;


			Zeus::ObjFilterRealParams	dummy;
			std::vector<LikePixBufferType> pixbufferArr(GlbVars_.ProfParam_.ContBinsDef_);
			for (int i = 0; i < GlbVars_.ProfParam_.ContBinsDef_; ++i)
			{
				pixbufferArr[i].InUse_ = 0;
				pixbufferArr[i].ISNR2_ = 0.0;
				pixbufferArr[i].pixbuffer.Make(((POS_LIM << 1) - 1) * ((POS_LIM << 1) - 1), ((POS_LIM << 1) - 1));
				pixbufferArr[i].pixbuffer.reset(Zeus::logZERO);
			}

			const double	FluxConvCte(1.0/((GlbVars_.PixSz_ * GlbVars_.PixSz_) * (0.001 * SR2ARCMIN2)));
// set priors
			Zeus::PriorsTypeCollType	tPriType;
			SetupPriors(tPriType,GlbVars_.Use2DFormula_);

			FluxBestFit		= ((PeakResult.PeakNonBlindInfo_.FluxBay_ > 0.0) ? PeakResult.PeakNonBlindInfo_.FluxBay_: PeakResult.PeakNonBlindInfo_.PredFlux_);
			RadiusBestFit	= ((PeakResult.PeakNonBlindInfo_.RadiusBay_ > 0.0) ? PeakResult.PeakNonBlindInfo_.RadiusBay_: PeakResult.PeakNonBlindInfo_.PredRadius_);
// set sigmas (range)
			NSigmaHigh		= NSIGMA;
			NSigmaLow		= NSIGMA;
// set amplitudes
			const double fluxLimite(PeakResult.PeakNonBlindInfo_.ErrFluxHigh_ > PeakResult.PeakNonBlindInfo_.ErrFluxLow_?PeakResult.PeakNonBlindInfo_.ErrFluxHigh_:PeakResult.PeakNonBlindInfo_.ErrFluxLow_);
			double	MaxAmpl(FluxBestFit + (fluxLimite * NSigmaHigh));
			double	MinAmpl(FluxBestFit - (fluxLimite * NSigmaLow));
			if(MinAmpl < 0.0)
			{MinAmpl = 0.0;}
#ifdef JUSTCENTRALPIXEL
			MinAmpl = 0.0;
#endif

			const double	deltaAmpl((MaxAmpl - MinAmpl)/static_cast<double>(YPLOTSZ-1));
// set radii
			const double radiiLimite(PeakResult.PeakNonBlindInfo_.ErrRadiusHigh_ > PeakResult.PeakNonBlindInfo_.ErrRadiusLow_?PeakResult.PeakNonBlindInfo_.ErrRadiusHigh_:PeakResult.PeakNonBlindInfo_.ErrRadiusLow_);
			double	MaxRadius(RadiusBestFit + (radiiLimite * NSigmaHigh));
			double	MinRadius(RadiusBestFit - (radiiLimite * NSigmaLow));
			if(MinRadius < 0.0)
			{MinRadius = 0.0;}
#ifdef JUSTCENTRALPIXEL
			MinRadius = 0.0;
#endif
			const double	deltaRadius((MaxRadius - MinRadius)/static_cast<double>(GlbVars_.ProfParam_.ContBinsDef_-1));


 			ClearISNRCache();
// set priors parameters 
// BUG in Gcc http://stackoverflow.com/questions/16143185/passing-const-stdauto-ptr-as-this-argument-of-stdauto-ptr-tpoperator-s
// Need to use regular pointers ... Grrrrrrrrrrrrr

			std::vector<Zeus::Priors*>	PriorsColl(CONTPRIOR_SIZE,static_cast<Zeus::Priors*>(0));
			try{
				ContoursSetPriorsParameters(tPriType, PriorsColl, MinRadius, MaxRadius, MinAmpl, MaxAmpl, PeakResult);
			}
			catch(...){
				for(int i=0;i<CONTPRIOR_SIZE;++i)
				{delete PriorsColl[i];}
				throw;
			}

			(Zeus::ConManager::Instance())->PrintStr2Console(std::wstring(L"Sampling started..."));

			double logPriorTheta_s,ISNR2;
			double C500MinValue(1.0), C500MaxValue(1.0), C500ValueStep(0.0);
			double alphaValueStep(0.0);
			double betaValueStep(0.0);
			double YtotY500ConvCte(SZCATALOGUEDEFAULTVALUE);
			double Y500MinValue(SZCATALOGUEDEFAULTVALUE), Y500MaxValue(SZCATALOGUEDEFAULTVALUE);
			double Y500MinValueStep(0.0);
			//
			int stepBeta, stepAlpha, maxProfK, stepPrint, profAlphaMax, profBetaMax, profC500Max;

			if(GlbVars_.ProfParamVar_.IsInit())
			{
				stepAlpha		= ProfVarGetAlphaStep(PeakResult.PeakNonBlindInfo_.SNR_);
				if(((GlbVars_.ProfParamVar_.alpha_.nbins_-1) / stepAlpha)<2) stepAlpha >>= 1;
				if ((GlbVars_.ProfParamVar_.alpha_.nbins_ - 1) < stepAlpha) stepAlpha = (GlbVars_.ProfParamVar_.alpha_.nbins_ - 1);
				if(stepAlpha<1) stepAlpha=1;

				stepBeta		= ProfVarGetBetaStep(PeakResult.PeakNonBlindInfo_.SNR_);
				if(((GlbVars_.ProfParamVar_.beta_.nbins_-1) / stepBeta)<2) stepBeta >>= 1;
				if ((GlbVars_.ProfParamVar_.beta_.nbins_ - 1) < stepBeta) stepBeta = (GlbVars_.ProfParamVar_.beta_.nbins_ - 1);
				if(stepBeta<1) stepBeta=1;

				profC500Max = ((GlbVars_.ProfParamVar_.C500_.nbins_ <= 1) ? 1 : GlbVars_.ProfParamVar_.C500_.nbins_);
				C500MinValue = GlbVars_.ProfParamVar_.C500_.min_;
				C500MaxValue = GlbVars_.ProfParamVar_.C500_.max_;
				C500ValueStep = (GlbVars_.ProfParamVar_.C500_.max_ - GlbVars_.ProfParamVar_.C500_.min_) / static_cast<double>((GlbVars_.ProfParamVar_.C500_.nbins_<=1) ? 1 : GlbVars_.ProfParamVar_.C500_.nbins_ - 1);
				alphaValueStep = ((GlbVars_.ProfParamVar_.alpha_.max_ - GlbVars_.ProfParamVar_.alpha_.min_) / static_cast<double>((GlbVars_.ProfParamVar_.alpha_.nbins_ <= 1) ? 1 : GlbVars_.ProfParamVar_.alpha_.nbins_ - 1));
				betaValueStep = ((GlbVars_.ProfParamVar_.beta_.max_ - GlbVars_.ProfParamVar_.beta_.min_) / static_cast<double>((GlbVars_.ProfParamVar_.beta_.nbins_ <= 1) ? 1 : GlbVars_.ProfParamVar_.beta_.nbins_ - 1));
				// assumes Ytot is Y @ R500 * R500_Ratio_
				Y500MinValue = MinAmpl / GetYtotConvY500(GlbVars_.ProfParamVar_.alpha_.min_, GlbVars_.ProfParamVar_.beta_.min_, GlbVars_.ProfParam_.MNFW_gamma_, GlbVars_.ProfParamVar_.C500_.min_, GlbVars_.ProfParam_.R500_Ratio_, Gammafun);
				Y500MaxValue = (MaxAmpl / GetYtotConvY500(GlbVars_.ProfParamVar_.alpha_.max_, GlbVars_.ProfParamVar_.beta_.max_, GlbVars_.ProfParam_.MNFW_gamma_, GlbVars_.ProfParamVar_.C500_.max_, GlbVars_.ProfParam_.R500_Ratio_, Gammafun));
				Y500MinValueStep = (Y500MaxValue - Y500MinValue) / static_cast<double>(YPLOTSZ - 1);

				if(GlbVars_.PriorMassMax_ < 0.0)
				{
					stepBeta=1;
					MaxBeta= GlbVars_.ProfParamVar_.beta_.max_;
					MinBeta= GlbVars_.ProfParamVar_.beta_.min_;
					BetaStep= ((GlbVars_.ProfParamVar_.beta_.max_- GlbVars_.ProfParamVar_.beta_.min_) / static_cast<double>(GlbVars_.ProfParamVar_.beta_.nbins_));
					BetaPrior = new Zeus::UniformPrior(GlbVars_.ProfParamVar_.beta_.min_, GlbVars_.ProfParamVar_.beta_.max_);
				}
				maxProfK		= ((profBetaMax=GlbVars_.ProfParamVar_.beta_.nbins_) * (profAlphaMax=GlbVars_.ProfParamVar_.alpha_.nbins_))/(stepBeta*stepAlpha);
				stepPrint		= 10*((maxProfK<=1)?1:(maxProfK>>1));
			}
			else
			{
				stepAlpha = stepBeta = maxProfK = stepPrint = profAlphaMax = profC500Max = profBetaMax = 1;
				C500MaxValue = C500MinValue = 1.0; C500ValueStep = 0.0;
			}

			const double Theta500Step(((MaxRadius*C500MaxValue) - (MinRadius*C500MinValue)) / static_cast<double>(GlbVars_.ProfParam_.ContBinsDef_ - 1));

			Zeus::LArr2D<double>		ContoursBuffer(GlbVars_.ProfParam_.ContBinsDef_ * GlbVars_.ProfParam_.ContBinsDef_,GlbVars_.ProfParam_.ContBinsDef_,Zeus::logZERO);
			/*
			{
				wchar_t buf[BUFFERMAXCHAR];
				PRINTINTOBUFFERFUNCT(buf, BUFFERMAXCHAR, L"\nThetasMin -> %9.3g,ThetasMax -> %9.3g,deltaThetas -> %9.3g,YMin -> %9.3g,YMax -> %9.3g,deltaY -> %9.3g,C500Min -> %9.3g,C500Max -> %9.3g,deltaC500 -> %9.3g\n", MinRadius, MaxRadius, deltaRadius, MinAmpl, MaxAmpl, deltaAmpl, C500MinValue, C500MaxValue, Theta500Step);
				(Zeus::ConManager::Instance())->PrintStr2Console(buf);
			}
			*/
			for(int profAlpha=0;profAlpha<profAlphaMax;profAlpha+= stepAlpha)
			{
				for (int profBeta = 0; profBeta<profBetaMax; profBeta += stepBeta)
				{
					const int CurrProfileIndex((profAlpha*profBetaMax)+profBeta);

					// reset likelihood buffer
					for (int i = 0; i < GlbVars_.ProfParam_.ContBinsDef_; ++i)
					{
						pixbufferArr[i].InUse_ = 0;
						pixbufferArr[i].ISNR2_ = Zeus::logZERO;
						pixbufferArr[i].pixbuffer.reset(Zeus::logZERO);
					}

					for (int profC500 = 0;profC500 < profC500Max;++profC500)
					{
// new cycle
						const double C500piv(C500MinValue + (static_cast<double>(profC500)*C500ValueStep));

						if (GlbVars_.ProfParamVar_.IsInit() && (profC500Max > 1))
						{
							YtotY500ConvCte = GetYtotConvY500(GlbVars_.ProfParamVar_.alpha_.min_ + (static_cast<double>(profAlpha)* alphaValueStep), GlbVars_.ProfParamVar_.beta_.min_ + (static_cast<double>(profBeta)* betaValueStep), GlbVars_.ProfParam_.MNFW_gamma_, C500piv, GlbVars_.ProfParam_.R500_Ratio_ , Gammafun);
						}

						for (int j = 0;j<GlbVars_.ProfParam_.ContBinsDef_;++j)
						{
							if (!((CurrProfileIndex*GlbVars_.ProfParam_.ContBinsDef_ + j) % stepPrint))
							wprintf(L".");

							const double Theta_s(((MinRadius*C500MinValue) + static_cast<double>(j)*Theta500Step) / C500piv);

							if ((Theta_s<(MinRadius - RADIUSTOL)) || (Theta_s>(MaxRadius + RADIUSTOL)))
								continue;

							int ThetaSInd(Zeus::toInt(((Theta_s - MinRadius) / deltaRadius)+0.5));
							 
							if (ThetaSInd < 0)
							{ThetaSInd = 0;}
							else if (ThetaSInd > (GlbVars_.ProfParam_.ContBinsDef_-1))
							{ThetaSInd = GlbVars_.ProfParam_.ContBinsDef_ - 1;}

//							const double Theta_s(MinRadius + ((static_cast<double>(j)* deltaRadius)));

							logPriorTheta_s = PriorsColl[CONTPRIOR_THETA]->GetPDF(Theta_s);
							logPriorTheta_s = ((logPriorTheta_s > 0.0) ? std::log(logPriorTheta_s) : Zeus::logZERO);

							
							if (!(pixbufferArr[ThetaSInd].InUse_))
							{
								RealPlaneSurfType	tempLikelihood;

								ISNR2 = MakeLikelihoodAtRealScale(Zeus::ObjFilterRealParams(Theta_s, GlbVars_.ProfParam_.MNFW_alpha_, GlbVars_.ProfParam_.MNFW_beta_, GlbVars_.ProfParam_.MNFW_gamma_, (GlbVars_.ProfParamVar_.IsInit() && (profC500Max > 1)) ? C500piv : GlbVars_.ProfParam_.MNFW_C500_, CurrProfileIndex), tempLikelihood);
								pixbufferArr[ThetaSInd].ISNR2_ = ISNR2;

								if (Zeus::pwsIsNaN(ISNR2))
								{
									wchar_t buf[BUFFERMAXCHAR];
									PRINTINTOBUFFERFUNCT(buf, BUFFERMAXCHAR, L"\nFilter returned NaN, Profile -> %d, radius -> %9.3f, SrcIndex -> %d, SNR -> %9.3f\n", CurrProfileIndex, Theta_s, PeakResult.PeakNonBlindInfo_.SrcIndex_, PeakResult.PeakNonBlindInfo_.SNR_);
									(Zeus::ConManager::Instance())->PrintStr2Console(buf);
									continue;
								}
								const double*	const	pivOrg((tempLikelihood.GetInnerData()).begin());
								const unsigned long		ArrMetric((tempLikelihood.GetInnerData()).getPtrMetric());
								int YOff, tY, tX;
								for (int l = -POS_LIM + 1; l < POS_LIM; ++l)
								{
									YOff = ArrMetric * (tY = (YCentralPix + l));
									for (int k = -POS_LIM + 1; k < POS_LIM; ++k)
									{
										if (IsOutsideDetectArea(tY, (tX = XCentralPix + k)))
											continue;
										pixbufferArr[ThetaSInd].pixbuffer(l + POS_LIM - 1, k + POS_LIM - 1) = *(pivOrg + YOff + tX);
									}
								}
								pixbufferArr[ThetaSInd].InUse_ = 1;
							}
							else
							{
							//	wchar_t buf[BUFFERMAXCHAR];
							//	PRINTINTOBUFFERFUNCT(buf, BUFFERMAXCHAR, L"\nOOooooooooooooooooooooooooooooo ps !!!!\n");
							//	(Zeus::ConManager::Instance())->PrintStr2Console(buf);
							}
							// Flux loop
							double logPost(Zeus::logZERO);
							for (int i = 0; i<YPLOTSZ; ++i, logPost = Zeus::logZERO)
							{
								double YFlux;
								if (!(GlbVars_.ProfParamVar_.IsInit()) || (profC500Max <= 1))
								{
									YFlux = MinAmpl + (static_cast<double>(i)*deltaAmpl);
								}
								else
								{
									const double Y500(Y500MinValue + static_cast<double>(i)* Y500MinValueStep);
									YFlux = YtotY500ConvCte * Y500;
									if ((YFlux < MinAmpl) || (YFlux > MaxAmpl))
										continue;
								}
								double logPriorYFLux(PriorsColl[CONTPRIOR_Y]->GetPDF(YFlux));
								logPriorYFLux = ((logPriorYFLux > 0.0) ? std::log(logPriorYFLux) : Zeus::logZERO);

								const double tIntFlux(YFlux*FluxConvCte);
								{
									for (int l = -POS_LIM + 1; l < POS_LIM; ++l)
									{
										for (int k = -POS_LIM + 1; k < POS_LIM; ++k)
										{
											logPost = Zeus::AddLog(logPost, tIntFlux * pixbufferArr[ThetaSInd].pixbuffer(l + POS_LIM - 1, k + POS_LIM - 1));
										}
									}
									logPost -= (0.5 * pixbufferArr[ThetaSInd].ISNR2_ * tIntFlux * tIntFlux);
								}
								// i and j swapped 
								// ContoursBuffer(i,j) = (logPost + logPriorYFLux + logPriorTheta_s);
								if (GlbVars_.PriorMassMax_ < 0.0)
								{
									double logPriorBeta(BetaPrior->GetPDF(GlbVars_.ProfParamVar_.beta_.min_ + (static_cast<double>(profBeta)* BetaStep)));
									logPriorBeta = ((logPriorBeta > 0.0) ? std::log(logPriorBeta) : Zeus::logZERO);
									ContoursBuffer(profBeta, j) = Zeus::AddLog(ContoursBuffer(profBeta, j), (logPost + logPriorYFLux + logPriorTheta_s + logPriorBeta));
								}
								else
								{
									ContoursBuffer(i, j) = Zeus::AddLog(ContoursBuffer(i, j), (logPost + logPriorYFLux + logPriorTheta_s));
								}
							}
							// new cycle
						}
					}
				}
			}
//
			for(int i=0;i<CONTPRIOR_SIZE;++i)
			{delete PriorsColl[i];}
			delete BetaPrior;
			Zeus::NormalizeProb(ContoursBuffer.begin(),ContoursBuffer.getSz(),ContoursBuffer.getPtrMetric());
			if(GlbVars_.PriorMassMax_ < 0.0)
			{
				PostProcessContours(PeakResult,ContoursBuffer, MaxBeta, MinBeta, MaxRadius, MinRadius);
			}
			else
			{
				PostProcessContours(PeakResult, ContoursBuffer, (!(GlbVars_.ProfParamVar_.IsInit()) || (profC500Max <= 1)) ? MaxAmpl : Y500MaxValue, (!(GlbVars_.ProfParamVar_.IsInit()) || (profC500Max <= 1)) ? MinAmpl : Y500MinValue, MaxRadius*C500MaxValue, MinRadius*C500MinValue);
			}
			(Zeus::ConManager::Instance())->PrintStr2Console(std::wstring(L"\nSampling ended ..."));
			++FileID;
		}

	default:
		return;
	}

	PeakResult.Pos_.YPix_			= Zeus::toInt(PeakResult.Pos_.YCoord_ + 0.5);
	PeakResult.Pos_.XPix_			= Zeus::toInt(PeakResult.Pos_.XCoord_ + 0.5);
	PeakResult.SrcAmplNormalised_	= std::sqrt(2.0 * PeakResult.Odds_);
	PeakResult.SrcFlux_				= PeakResult.SrcAmplNormalised_ / std::sqrt(PeakResult.ISNR2_);
	if (sign)
	{
		PeakResult.SrcFlux_			*= static_cast<double>(sign);
	}
}
//
void	PatchProcessor::ContoursSetPriorsParameters(Zeus::PriorsTypeCollType& priorsTypes, std::vector<Zeus::Priors*> & PriorsColl, double MinRadius, double MaxRadius, double MinAmpl, double MaxAmpl,const Zeus::PeakType &PeakResult)
{
	if(priorsTypes.size() < CONTPRIOR_SIZE)
		throw Zeus::libException(ERROR_COD_ZEUSNUMWRONGARGS,L"Currently PwS only supports 4 parameters.",this);

	Zeus::MarginalPrior*	tMarginal(0);
	Zeus::ConditionalPrior*	tConditional(0);

	int DummyFull,BeamHalfRadiusPix;
	double CorrectedBeamSigma;

	EvalEquivBeamSigmaPix(DummyFull,BeamHalfRadiusPix);	
	// Priors on position are currently hardcoded to uniform
	// X
	PriorsColl[CONTPRIOR_XPOS]= new Zeus::UniformPrior();
	// Y
	PriorsColl[CONTPRIOR_YPOS]= new Zeus::UniformPrior();

	// Theta
	switch(priorsTypes[CONTPRIOR_THETA])
	{
	case	Zeus::Priors::Gaussian:
			if ((PeakResult.PeakNonBlindInfo_.PredRadius_ < 0.0) || (PeakResult.PeakNonBlindInfo_.ErrRadius_ <= 0.0))
			{
				const double errValue((PeakResult.PeakNonBlindInfo_.PredRadius_ < 0.0) ? PeakResult.PeakNonBlindInfo_.PredRadius_ : PeakResult.PeakNonBlindInfo_.ErrRadius_);
				errPriorValues(L"NestSamplPwSObserParamNew::SetRadiusGaussianPriorParams", PeakResult.PeakNonBlindInfo_.SrcIndex_, errValue);
			}
			PriorsColl[CONTPRIOR_THETA] = new Zeus::GaussianPrior(PeakResult.PeakNonBlindInfo_.RadiusBay_, PeakResult.PeakNonBlindInfo_.ErrRadius_);
			break;
	case	Zeus::Priors::Uniform:
			PriorsColl[CONTPRIOR_THETA]= new Zeus::UniformPrior(MinRadius,MaxRadius);
			break;
	case	Zeus::Priors::Power:
			PriorsColl[CONTPRIOR_THETA]= new Zeus::PowerlawPrior(MinRadius,MaxRadius,GlbVars_.PriorSrcScaleExp_);
			break;
	case	Zeus::Priors::Exponential:
			PriorsColl[CONTPRIOR_THETA]= new Zeus::ExponencialPrior(MinRadius,MaxRadius,GlbVars_.PriorSrcScaleExp_);
			break;
	case	Zeus::Priors::RadiusNinf:
			CorrectedBeamSigma = (GlbVars_.PixSz_ * RAD2ARCMIN *  BeamHalfRadiusPix * 0.8493218) ; //0.8493218 = 2/SQRT(8*LN(2)) 
			CorrectedBeamSigma *= CorrectedBeamSigma;
			CorrectedBeamSigma *= (1.0/GlbVars_.ProfParam_.MNFW_2ndMoment1D_);
			PriorsColl[CONTPRIOR_THETA]= new Zeus::RadiusNinfPrior(std::sqrt(CorrectedBeamSigma),MinRadius,MaxRadius);
			break;
	case	Zeus::Priors::Marginal:
			{
				std::vector<double> tDummy(PriorsThetaMarginal,PriorsThetaMarginal+PriorsThetaMarginalSIZE);
				PriorsColl[CONTPRIOR_THETA] =  tMarginal = new Zeus::MarginalPrior(tDummy);
			}
			break;
	default :
			throw Zeus::libException(ERROR_COD_TYPEMISMATCH,L"A prior not allowed was used.",L"Theta_s");
			break;
	}

	// GlbVars_.Use2DFormula_
	switch(priorsTypes[CONTPRIOR_Y])
	{		
	case	Zeus::Priors::Uniform:
			PriorsColl[CONTPRIOR_Y]= new Zeus::UniformPrior(MinAmpl,MaxAmpl);
//			PriorsColl[CONTPRIOR_Y]= new Zeus::PowerlawPrior(MinAmpl,MaxAmpl,-0.666666); // - 2/3
			break;
	case	Zeus::Priors::Power:
			PriorsColl[CONTPRIOR_Y]= new Zeus::PowerlawPrior(MinAmpl,MaxAmpl,-0.666666); // - 2/3
			break;
	case	Zeus::Priors::Exponential:
			PriorsColl[CONTPRIOR_Y]= new Zeus::ExponencialPrior(MinAmpl,MaxAmpl,GlbVars_.PriorFluxExp_);
			break;
	case	Zeus::Priors::Conditional:
			{
				std::vector<double> tDummy(PriorsYcondTheta,PriorsYcondTheta+PriorsYcondThetaSIZE);
				std::vector<double> tDummyPDf(PriorsCondPDF,PriorsCondPDF+PriorsCondPDFSIZE);
				PriorsColl[CONTPRIOR_Y] =  tConditional = new Zeus::ConditionalPrior(tDummy,tDummyPDf,(PriorsYcondThetaSIZE / PriorsThetaMarginalSIZE));
			}
			break;
	default :
			throw Zeus::libException(ERROR_COD_TYPEMISMATCH,L"A prior not allowed was used.",L"Y5r500");
			break;
	}
//
	if(tMarginal && tConditional)
	{
		tMarginal->SetSlaveObject(tConditional);
	}
}
// Assumes priors on position are ALWAYS flat
double	PatchProcessor::MargPosGridSample(int YCentralPix,int XCentralPix,int PixRange,double thetaS,double FluxIntUnits)
{
	double logPost(Zeus::logZERO);
	Zeus::ObjFilterRealParams in,dummy;

	in.RealScale_ = thetaS;

	for(int l=-PixRange+1;l<PixRange;++l)
	{
		for(int k=-PixRange+1;k<PixRange;++k)
		{logPost = Zeus::AddLog(OddsEval_->GetOddsResultRaw(YCentralPix+l,XCentralPix+k,in,dummy,FluxIntUnits),logPost);}
	}

	return logPost;
}
#ifdef AMI
int	PatchProcessor::AmiPwSLikelihood(double* LikeValue, const AmiLikeParams& ParamsIN, AmiLikeParams& ParamsOut, bool Noduocimation)
{
	const int	XCentralPix((PatchesGeoInfo_.Storage_[PatchNumber_].SrcXCoord_ >= 0) ? PatchesGeoInfo_.Storage_[PatchNumber_].SrcXCoord_ : (GlbVars_.PatchSz_ >> 1));
	const int	YCentralPix((PatchesGeoInfo_.Storage_[PatchNumber_].SrcYCoord_ >= 0) ? PatchesGeoInfo_.Storage_[PatchNumber_].SrcYCoord_ : (GlbVars_.PatchSz_ >> 1));
	const double	FluxConvCte(1.0 / ((GlbVars_.PixSz_ * GlbVars_.PixSz_) * (0.001 * SR2ARCMIN2)));
	const double pixSzArcSec((GlbVars_.PixSz_ * 648000.0) / PI);
	const int		yPixShift(Zeus::toInt((ParamsIN.YCoord_ / pixSzArcSec) + 0.5));
	const int		xPixShift(Zeus::toInt((ParamsIN.XCoord_ / pixSzArcSec) + 0.5));
	const double	yShiftOff((ParamsIN.YCoord_ / pixSzArcSec) - static_cast<double>(yPixShift));
	const double	xShiftOff((ParamsIN.XCoord_ / pixSzArcSec) - static_cast<double>(xPixShift));

	ParamsOut.YCoord_ = (yShiftOff + static_cast<double>(yPixShift)) *pixSzArcSec;
	ParamsOut.XCoord_ = (xShiftOff + static_cast<double>(xPixShift)) *pixSzArcSec;

	int CurrProfileIndex(-1);

	if (GlbVars_.ProfParamVar_.IsInit() && (ParamsIN.NParams_ >= 6))
	{
		const double	AlphaStep((GlbVars_.ProfParamVar_.alpha_.max_ - GlbVars_.ProfParamVar_.alpha_.min_) / static_cast<double>((GlbVars_.ProfParamVar_.alpha_.nbins_ <= 1) ? 1 : GlbVars_.ProfParamVar_.alpha_.nbins_ - 1));
		const double	BetaStep((GlbVars_.ProfParamVar_.beta_.max_ - GlbVars_.ProfParamVar_.beta_.min_) / static_cast<double>((GlbVars_.ProfParamVar_.beta_.nbins_ <= 1) ? 1 : GlbVars_.ProfParamVar_.beta_.nbins_ - 1));

		int		AlphaPix(Zeus::toInt(((ParamsIN.alpha_ - GlbVars_.ProfParamVar_.alpha_.min_) / AlphaStep) + 0.5));
		if (AlphaPix < 0) AlphaPix = 0;
		if (AlphaPix >= GlbVars_.ProfParamVar_.alpha_.nbins_) AlphaPix = (GlbVars_.ProfParamVar_.alpha_.nbins_-1);
		ParamsOut.alpha_ = GlbVars_.ProfParamVar_.alpha_.min_ + (AlphaPix * AlphaStep);
		int		BetaPix(Zeus::toInt(((ParamsIN.beta_ - GlbVars_.ProfParamVar_.beta_.min_) / BetaStep) + 0.5));
		if (BetaPix < 0) BetaPix = 0;
		if (BetaPix >= GlbVars_.ProfParamVar_.beta_.nbins_) BetaPix = (GlbVars_.ProfParamVar_.beta_.nbins_ - 1);

		ParamsOut.beta_ = GlbVars_.ProfParamVar_.beta_.min_ + (BetaPix * BetaStep);

		CurrProfileIndex = ((AlphaPix*((GlbVars_.ProfParamVar_.beta_.nbins_ < 1) ? 1 : GlbVars_.ProfParamVar_.beta_.nbins_)) + BetaPix);
	}
	else
	{
		ParamsOut.alpha_ = GlbVars_.ProfParam_.MNFW_alpha_;
		ParamsOut.beta_ = GlbVars_.ProfParam_.MNFW_beta_;
	}

	Zeus::ObjFilterRealParams in(ParamsOut.ThetaS_= ParamsIN.ThetaS_, GlbVars_.ProfParam_.MNFW_alpha_, GlbVars_.ProfParam_.MNFW_beta_, GlbVars_.ProfParam_.MNFW_gamma_, GlbVars_.ProfParam_.MNFW_C500_, CurrProfileIndex);

	*LikeValue = OddsEval_->GetOddsResultProfVarRaw(YCentralPix + yPixShift, XCentralPix + xPixShift, yShiftOff, xShiftOff, in, Noduocimation, (ParamsOut.Ytot_ = ParamsIN.Ytot_)*FluxConvCte);

	return 0;
}
#endif //AMI
//
void	PatchProcessor::PostProcessContours(Zeus::PeakType &PeakResult,Zeus::LArr2D<double>	&Buffer,double MaxAmpl,double MinAmpl,double MaxRadius,double MinRadius)
{
	if(GlbVars_.Sigma8_ < 0.0)
	{
		Zeus::ContoursOutAux aux(PeakResult.PeakNonBlindInfo_.QAIN_CyR500,PeakResult.PeakNonBlindInfo_.QAIN_T500);
		if((PeakResult.PeakNonBlindInfo_.QAIN_CyR500 < -1.0e+20) || (PeakResult.PeakNonBlindInfo_.QAIN_T500 < -1.0e+20))
		{WriteDegeneracyData2Media(Buffer,PeakResult,MaxRadius,MinRadius,MaxAmpl,MinAmpl,NULL);}
		else
		{WriteDegeneracyData2Media(Buffer,PeakResult,MaxRadius,MinRadius,MaxAmpl,MinAmpl,&aux);}
	}
	else
	{
		const int		nscales(GlbVars_.ProfParam_.ContBinsDef_);
		const double	input_theta_s(PeakResult.PeakNonBlindInfo_.QAIN_T500);
		const double	input_cy5r500(PeakResult.PeakNonBlindInfo_.QAIN_CyR500);

#ifndef MYCONTOURSASS
		QA_CONTOURS_FUNCTION(&nscales,&input_theta_s,&input_cy5r500,&MinRadius,
			&MaxRadius,&MinAmpl,&MaxAmpl,Buffer.begin(),PeakResult.QAResult_+ 4);
#else
		myQA_ContAssess(nscales,input_theta_s,input_cy5r500,MinRadius,
			MaxRadius, MinAmpl,MaxAmpl,Buffer.begin(),PeakResult.QAResult_+ 4,GlbVars_.Jf_Estimator_);
#endif
		PeakResult.QAResult_[0]		= PeakResult.PeakNonBlindInfo_.QAIN_SrcPtg_.theta;
		PeakResult.QAResult_[1]		= PeakResult.PeakNonBlindInfo_.QAIN_SrcPtg_.phi;
		PeakResult.QAResult_[2]		= input_theta_s;
		PeakResult.QAResult_[3]		= input_cy5r500;
		PeakResult.CollListIndex_	= (1.0 * static_cast<double>(PeakResult.PeakNonBlindInfo_.CollLstIndex_));
	}
}
//
void	PatchProcessor::WriteDegeneracyData2Media(const Zeus::LArr2D<double>& data,const Zeus::PeakType& PeakResult, double XMax, double XMin, double YMax, double YMin,Zeus::ContoursOutAux* aux)
{

	std::wstring	ObjectName(MakeDegeneracyPlotsFileName(PeakResult));

	std::auto_ptr<Zeus::GenCollWriter<Zeus::LArr2D<double> > > 
		FWriter0(Zeus::GetGenCollFileWriterHandler(Loki::Type2Type<Zeus::LArr2D<double> >(),GlbVars_.ContextID_,
		ObjectName,
		GlbVars_.DirOut_ + std::wstring(L"_DegCont"),data.getPtrMetric(),GlbVars_.PixSz_,YMax,YMin,XMax,XMin,aux));

	FWriter0->Initialize();
	FWriter0->Write(data);
	FWriter0->Flush();

#if defined(HFIDMC) && defined(HFIDMC_EXTOBJECTS)
	if(GlbVars_.Data2Buffer_)
	{
#ifdef WIN32
		std::wstring	degExt(L"DegCont\\");
#else
		std::wstring	degExt(L"DegCont/");
#endif
		std::wstring	FullPath(GlbVars_.DirBuffer_ + degExt);
		Zeus::CreateDir(FullPath);

		std::auto_ptr<Zeus::GenCollWriter<Zeus::LArr2D<double> > > 
			FWriter1(Zeus::GetGenCollFileWriterHandler(Loki::Type2Type<Zeus::LArr2D<double> >(),1000,
			ObjectName,
			FullPath,GlbVars_.PatchSz_,GlbVars_.PixSz_,YMax,YMin,XMax,XMin,aux));

		FWriter1->Initialize();
		FWriter1->Write(data);
		FWriter1->Flush();
	}
#endif

}
//
double	PatchProcessor::GetOddsResultRawOptimAmpRange(int YPix,int XPix,int range, const Zeus::ObjFilterRealParams& inFilterParams,double& Isnr,int& sign)
{
	Zeus::ObjFilterRealParams outFilterParams;
	double		t;
	double		MaxCorr(Zeus::logZERO);

	ObjFilterParams FilterParams(Zone_->TranslateParamsInverse(inFilterParams));
	SrcInfoType		ObjSurf(Zone_->GetSrcObj(FilterParams,outFilterParams,false));
	Isnr =			GetISNR(FilterParams,ObjSurf);

	for(int j=YPix-range;j<=YPix+range;++j)
	{
		for(int i=XPix-range;i<=XPix+range;++i)
		{
			if(MaxCorr < (t=Zeus::Correlation(LikelihoodSurf_,ObjSurf.surf_,j,i)))
			{MaxCorr=t;}
		}
	}
	sign = (MaxCorr<0.0?-1:1);
	return  (MaxCorr * MaxCorr) / (2.0 * Isnr);
}
//
int		PatchProcessor::RemoveSrcGreyArea(void)
{
	Zeus::PeakCollType::iterator	newEnd(std::remove_if(PeakColl_.begin(),PeakColl_.end(),
		IntCatalFilterFunctor(GlbVars_.PatchBorder_,GlbVars_.PatchSz_)
		));
	if(newEnd != PeakColl_.end())
	{PeakColl_.erase(newEnd,PeakColl_.end());}

	return static_cast<int>(PeakColl_.size());
}
//
int		PatchProcessor::ProfVarGetBetaStep(double snr)
{
	if(snr >= 18.0)	return 1;
	if(snr >= 9.0)  return 2;
	if(snr >= 6.0)	return 4;
	return 8;
}
//
int		PatchProcessor::ProfVarGetAlphaStep(double snr)
{
	if(snr >= 7.0)	return 1;
	return 2;
}

double GetYtotConvY500(double alpha, double beta, double gamma, double c500, double R500Ratio, Zeus::GammaFuncts& Gammafun)
{
	const double Y_sphR500(Gammafun.betainc((3.0 - gamma) / alpha, (beta - 3.0) / alpha, std::pow(c500, alpha) / (1.0 + std::pow(c500, alpha))));
	const double Y_sphTotR500(Gammafun.betainc((3.0 - gamma) / alpha, (beta - 3.0) / alpha, std::pow(R500Ratio * c500, alpha) / (1.0 + std::pow(R500Ratio * c500, alpha))));

	return Y_sphTotR500 / Y_sphR500;
}


