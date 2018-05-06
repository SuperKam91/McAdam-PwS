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


#include "PWS_Globals.h"
#include "PWS_NestSamplCases.h"

const int	NestSamplPwSObserParamNew::MaxDims_ = 4; //Number of dimensions

NestSamplPwSObserParamNew::NestSamplPwSObserParamNew(const Zeus::PriorsTypeCollType& PriorsCollT,OddsEval& LikeEval,
						int NLivePoints,int EvNIndivSamples,int MaxSteps,
						double FractTolEv,double xtraEnlFactor,
						int BeamFullRadiusPix,int BeamHalfRadiusPix,int seed)
		:	NestSamplEllipBound(NLivePoints,EvNIndivSamples,MaxSteps,FractTolEv,xtraEnlFactor,seed),
			LikeEval_(LikeEval),
			BeamFullRadiusPix_(BeamFullRadiusPix),
			BeamHalfRadiusPix_(BeamHalfRadiusPix),
			GlbVars_((PlanckInfo::Instance())->GetGlobalVars()),
			InnerCore_(0.0)
{
	if(PriorsCollT.size() < MaxDims_)
		throw Zeus::libException(ERROR_COD_ZEUSNUMWRONGARGS,ERROR_MSG_ZEUSNUMWRONGARGS,this);
	
	PriorsSlots_.resize(MaxDims_,0);
	PriorsSlotsCollType::iterator		piv(PriorsSlots_.begin());
	PriorsSlotsCollType::const_iterator	const end(PriorsSlots_.end());
	Zeus::PriorsTypeCollType::const_iterator	pivT(PriorsCollT.begin());
	double	CorrectedBeamSigma;
	Zeus::MarginalPrior*	tMarginal(0);
	Zeus::ConditionalPrior*	tConditional(0);


  	for(;piv != end;++piv,++pivT)
	{
		switch(*pivT)
		{
		case	Zeus::Priors::Uniform:
				*piv = new Zeus::UniformPrior();
				break;
		case	Zeus::Priors::Power:
				*piv = new Zeus::PowerlawPrior(GlbVars_.PriorFluxMin_,GlbVars_.PriorFluxMax_,GlbVars_.PriorFluxExp_);
				break;
		case	Zeus::Priors::Exponential:
				*piv = new Zeus::ExponencialPrior(GlbVars_.PriorSrcMinScale_,GlbVars_.PriorSrcMaxScale_,GlbVars_.PriorSrcScaleExp_);
				break;
		case	Zeus::Priors::Gaussian:
				*piv =  new Zeus::GaussianPrior();
				break;
		case	Zeus::Priors::RadiusNinf:
				CorrectedBeamSigma = (GlbVars_.PixSz_ * RAD2ARCMIN *  BeamHalfRadiusPix * 0.8493218) ; //0.8493218 = 2/SQRT(8*LN(2))
				CorrectedBeamSigma *= CorrectedBeamSigma;
				if(GlbVars_.SZ_)
				{CorrectedBeamSigma *= (1.0/GlbVars_.ProfParam_.MNFW_2ndMoment1D_);}
				*piv =  new Zeus::RadiusNinfPrior(std::sqrt(CorrectedBeamSigma),GlbVars_.PriorSrcMinScale_,GlbVars_.PriorSrcMaxScale_);
				break;
		case	Zeus::Priors::Marginal:
				{
					if(tMarginal)
					{
						DisposePriors();
						throw Zeus::libException(ERROR_COD_TYPEMISMATCH,L"Only one marginal prior allowed",this);					
					}
					std::vector<double> tDummy(PriorsThetaMarginal,PriorsThetaMarginal+PriorsThetaMarginalSIZE);
					*piv =  tMarginal = new Zeus::MarginalPrior(tDummy);
				}
				break;
		case	Zeus::Priors::Conditional:
				{
					if(tConditional)
					{
						DisposePriors();
						throw Zeus::libException(ERROR_COD_TYPEMISMATCH,L"Only one conditional prior allowed",this);					
					}
					std::vector<double> tDummy(PriorsYcondTheta,PriorsYcondTheta+PriorsYcondThetaSIZE);
					std::vector<double> tDummyPDf(PriorsCondPDF,PriorsCondPDF+PriorsCondPDFSIZE);

					*piv =  tConditional = new Zeus::ConditionalPrior(tDummy,tDummyPDf,(PriorsYcondThetaSIZE / PriorsThetaMarginalSIZE));
				}
				break;
		default :
				DisposePriors();
				throw Zeus::libException(ERROR_COD_TYPEMISMATCH,L"Prior not allowed",this);
				break;
		}
	}
	
	if(tMarginal && tConditional)
	{
		tMarginal->SetSlaveObject(tConditional);
	}

	FluxConvCte_	= (GlbVars_.ProfParam_.FluxCalibCte_ * GlbVars_.PixSz_ * GlbVars_.PixSz_) * (GlbVars_.SZ_?(0.001 * SR2ARCMIN2):1.0e9);
	FluxConvCte_	=  1.0 / FluxConvCte_;
}

double	NestSamplPwSObserParamNew::do_Reset(int& NPriorDim,int& NLikeDim)
{
	NPriorDim = NLikeDim = MaxDims_;
	return PeakProps_->Odds_;
}

int		NestSamplPwSObserParamNew::do_TranslateSample(MyArrType& LikelihoodSample,const MyArrType& PriorSample)
{
	LikelihoodSample[0] =	PriorsSlots_[0]->GetSample(PriorSample[0]);
	LikelihoodSample[1] =	PriorsSlots_[1]->GetSample(PriorSample[1]);
	LikelihoodSample[2] =	PriorsSlots_[2]->GetSample(PriorSample[2]);
	LikelihoodSample[3] =	(PriorsSlots_[3]->GetSample(PriorSample[3])) * FluxConvCte_;
	return 1;
}

double		NestSamplPwSObserParamNew::do_LikelihoodEval(const MyArrType& LikelihoodSample)
{

	Zeus::ObjFilterRealParams	in(LikelihoodSample[2]);
	Zeus::ObjFilterRealParams	out;

	double						Isnr;
	int							sign;
	const double result(LikeEval_.GetOddsResultRaw(Zeus::toInt(LikelihoodSample[1] + 0.5),Zeus::toInt(LikelihoodSample[0] + 0.5),in,out,LikelihoodSample[3]));
	const double x(static_cast<double>(PeakProps_->Pos_.XCoord_));
	const double y(static_cast<double>(PeakProps_->Pos_.YCoord_));
	const double xdif(LikelihoodSample[0] - x);
	const double ydif(LikelihoodSample[1] - y);
	const double rad2(xdif*xdif+ydif*ydif);
	if(rad2 <= InnerCore_) return result;
	if(result >= PeakProps_->Odds_) return Zeus::logZERO;
	const double result1(LikeEval_.GetOddsResultRawOptimAmp(Zeus::toInt(((LikelihoodSample[1] + y) / 2.0) + 0.5), Zeus::toInt(((LikelihoodSample[0] + x) / 2.0) + 0.5), in, Isnr, sign));
	if(result1 < result) return Zeus::logZERO;
	return result;
}
//
void	NestSamplPwSObserParamNew::SetPosGaussianPriorParams(Zeus::PeakType* peak,const Zeus::PatchGeomLineType& SrcGeom)
{
	Zeus::Priors::PriorParamsType	params;

	const double	XCentralPix(static_cast<double>((SrcGeom.SrcXCoord_ >= 0)?SrcGeom.SrcXCoord_:(GlbVars_.PatchSz_ >> 1)));
	const double	YCentralPix(static_cast<double>((SrcGeom.SrcYCoord_ >= 0)?SrcGeom.SrcYCoord_:(GlbVars_.PatchSz_ >> 1)));

	if(peak->PeakNonBlindInfo_.ErrPos_ <= 0.0)
		errPriorValues(L"NestSamplPwSObserParamNew::SetPosGaussianPriorParams",peak->PeakNonBlindInfo_.SrcIndex_,peak->PeakNonBlindInfo_.ErrPos_);

	const double posUncert(peak->PeakNonBlindInfo_.ErrPos_ / (SQRT2 * GlbVars_.PixSz_   * RAD2ARCMIN));
	params.push_back(XCentralPix);
	params.push_back(posUncert);
	PriorsSlots_[0]->SetParams(params);
	params.at(0) = YCentralPix;
	PriorsSlots_[1]->SetParams(params);
	EvCorrectingCte_	= 0.0;
}
//
void	NestSamplPwSObserParamNew::SetRadiusGaussianPriorParams(Zeus::PeakType* peak,const Zeus::PatchGeomLineType& SrcGeom)
{
	Zeus::Priors::PriorParamsType	params;

	if((peak->PeakNonBlindInfo_.PredRadius_ < 0.0) || (peak->PeakNonBlindInfo_.ErrRadius_ <= 0.0))
	{
		double errValue((peak->PeakNonBlindInfo_.PredRadius_ < 0.0)?peak->PeakNonBlindInfo_.PredRadius_:peak->PeakNonBlindInfo_.ErrRadius_);
		errPriorValues(L"NestSamplPwSObserParamNew::SetRadiusGaussianPriorParams",peak->PeakNonBlindInfo_.SrcIndex_,errValue);
	}
	params.push_back(peak->PeakNonBlindInfo_.PredRadius_);
	params.push_back(peak->PeakNonBlindInfo_.ErrRadius_);
	PriorsSlots_[2]->SetParams(params);
}

//
void	NestSamplPwSObserParamNew::SetFluxGaussianPriorParams(Zeus::PeakType* peak,const Zeus::PatchGeomLineType& SrcGeom)
{
	Zeus::Priors::PriorParamsType	params;

	params.push_back(peak->PeakNonBlindInfo_.PredFlux_);
	params.push_back(peak->PeakNonBlindInfo_.ErrFlux_);

	if((peak->PeakNonBlindInfo_.PredFlux_ < -1.0e30) || (peak->PeakNonBlindInfo_.ErrFlux_ <= 0.0))
	{
		double errValue((peak->PeakNonBlindInfo_.PredFlux_ < 0.0)?peak->PeakNonBlindInfo_.PredFlux_:peak->PeakNonBlindInfo_.ErrFlux_);
		errPriorValues(L"NestSamplPwSObserParamNew::SetRadiusGaussianPriorParams",peak->PeakNonBlindInfo_.SrcIndex_,errValue);
	}

	PriorsSlots_[3]->SetParams(params);
}
//
void	NestSamplPwSObserParamNew::SetFluxUniformPriorParams(Zeus::PeakType* peak,const Zeus::PatchGeomLineType& SrcGeom)
{
	double	PredFluxErrorBar(1.0 / std::sqrt(peak->ISNR2_));
	double	FluxCte(GlbVars_.ProfParam_.FluxCalibCte_ * GlbVars_.PixSz_ * GlbVars_.PixSz_);
	double	Flux(peak->SrcFlux_);

	if(GlbVars_.SZ_)
	{
		FluxCte				*= (0.001 * SR2ARCMIN2);
	}
	else
	{
		FluxCte				*= 1.0e9;
	}
	
	PredFluxErrorBar	*= FluxCte;
	Flux				*= FluxCte;

	double	FluxMin(Flux - (5.0 * PredFluxErrorBar));
	double	FluxMax(Flux + (5.0 * PredFluxErrorBar));

	if(FluxMin < GlbVars_.PriorFluxMin_)
	{
		FluxMin = GlbVars_.PriorFluxMin_;
	}

	if(FluxMax > GlbVars_.PriorFluxMax_)
	{
		FluxMax = GlbVars_.PriorFluxMax_;
	}

	Zeus::Priors::PriorParamsType	params;

	params.push_back(FluxMin);
	params.push_back(FluxMax);
	PriorsSlots_[3]->SetParams(params);
}
//
void	NestSamplPwSObserParamNew::SetPosFlatPriorBounds(Zeus::PeakType* peak,const Zeus::PatchGeomLineType& SrcGeom)
{
	int YMin,YMax,XMin,XMax;

	const double	MaxSzCte(GlbVars_.SZ_?GlbVars_.ProfParam_.VirialRatio_:1.0);
	int	PriorPosRange(BeamFullRadiusPix_ + ((Zeus::toInt(((peak->RealParams_.RealScale_ * MaxSzCte )/ (GlbVars_.PixSz_   * RAD2ARCMIN)) + 1.0))<<1));  
	PriorPosRange	<<= 1;


	double PriorArea(static_cast<double>(GlbVars_.PatchSz_ - (GlbVars_.PatchBorder_<<1)));
	PriorArea	*= PriorArea;
	PriorArea	/= GlbVars_.PriorAvObjectsPatch_;

	if(PriorPosRange > Zeus::toInt(std::sqrt(PriorArea) + 0.5))
	{
		PriorPosRange = Zeus::toInt(std::sqrt(PriorArea) + 0.5);
	}

	PriorPosRange	>>= 1;

// Making life easier for bright PS
	if(!(GlbVars_.SZ_) &&
		(peak->SrcAmplNormalised_ > THRESHOLDFORMASKING)
		)
	{
		if(PriorPosRange > 2)
			PriorPosRange = 2;
	}

	XMin	= peak->Pos_.XPix_ - PriorPosRange;
	if(XMin < GlbVars_.PatchBorder_)	XMin = GlbVars_.PatchBorder_;
	XMax	= peak->Pos_.XPix_ + PriorPosRange;
	if(XMax > (GlbVars_.PatchSz_ - GlbVars_.PatchBorder_))	XMax = (GlbVars_.PatchSz_ - GlbVars_.PatchBorder_);

	YMin	= peak->Pos_.YPix_ - PriorPosRange;
	if(YMin < GlbVars_.PatchBorder_)	YMin = GlbVars_.PatchBorder_;
	YMax	= peak->Pos_.YPix_ + PriorPosRange;
	if(YMax > (GlbVars_.PatchSz_ - GlbVars_.PatchBorder_))	YMax = (GlbVars_.PatchSz_ - GlbVars_.PatchBorder_);

	const double usedPriorArea(static_cast<double>(XMax - XMin) * static_cast<double>(YMax - YMin));

	Zeus::Priors::PriorParamsType	params;

	params.push_back(XMin);params.push_back(XMax);
	PriorsSlots_[0]->SetParams(params);
	params.clear();
	params.push_back(YMin);params.push_back(YMax);
	PriorsSlots_[1]->SetParams(params);

	EvCorrectingCte_	= std::log(usedPriorArea / PriorArea);
}


//
void	NestSamplPwSObserParamNew::NextPeak(Zeus::PeakType* peak,const Zeus::PatchGeomLineType& SrcGeom)
{
	PeakProps_		= peak;

	if(PriorsSlots_[0]->GetPriorType() == Zeus::Priors::Gaussian)
	{
		SetPosGaussianPriorParams(peak,SrcGeom);
	}
	else
	{
		SetPosFlatPriorBounds(peak,SrcGeom);
	}

	if(PriorsSlots_[2]->GetPriorType() == Zeus::Priors::Gaussian)
	{
		SetRadiusGaussianPriorParams(peak,SrcGeom);
	}

	switch(PriorsSlots_[3]->GetPriorType())
	{
	case Zeus::Priors::Gaussian:
		SetFluxGaussianPriorParams(peak,SrcGeom);
		break;
	case Zeus::Priors::Uniform:
		SetFluxUniformPriorParams(peak,SrcGeom);
		break;
	default:
		break;
	};
//

	const double t(static_cast<double>(BeamHalfRadiusPix_));
	const double t1(peak->RealParams_.RealScale_/ (GlbVars_.PixSz_   * RAD2ARCMIN));

	InnerCore_		= (t + t1)*(t + t1);

	NestedSamplerBase::Initialise();
}


//-----------------------------------------------------------------------

//
// ---------------------------------------------------
// ----- Simple example with ellipsoidal bounding  ---
// --------         Implementation          ----------

double	NSigma[DIMS] = {1.0,0.5,0.5,1.0};

double		NestSamplSimplGaussWithBounds::do_Reset(int& NPriorDim,int& NLikeDim)
{
	NPriorDim		= MaxDims_;
	NLikeDim		= MaxDims_;
	return MaxLike_;
}
int			NestSamplSimplGaussWithBounds::do_TranslateSample(MyArrType& LikelihoodSample,const MyArrType& PriorSample)
{
	for(int i=0;i<MaxDims_;++i)
	{
		LikelihoodSample[i] = (PriorSample[i] - 0.5) * PriorMaxRange_;
	}
	return 1;
}
double		NestSamplSimplGaussWithBounds::do_LikelihoodEval(const MyArrType& LikelihoodSample)
{
	double LnLike(0.0);

	for(int i=0;i<MaxDims_;++i)
	{
		LnLike += ((LikelihoodSample[i] * LikelihoodSample[i])/(NSigma[i]*NSigma[i]));
	}

	return  MaxLike_ - (0.5 * LnLike) ;
}

//
// ---------------------------------------------------
// -------- Simple example just for testing ----------
// --------         Implementation          ----------

double				NestSamplSimpleGaussian::do_Initialise(void)
{
	return MaxLike_;
}

int					NestSamplSimpleGaussian::do_GetInitialSamples(int NSamples,LivePointsSetType& Values)
{
	MyArrType	temp(MaxDims_);
	double		total;

	Values.clear();

	for(int i=0;i<NSamples;++i)
	{
		total = 0.0;
		for(int j=0;j<MaxDims_;++j)
		{
			temp[j]		=	(RandomGen_.RandDouble() - 0.5) * PriorMaxRange_;
			total		+=	(temp[j] * temp[j]);
		}

		Values.push_back(LivePointType(temp,temp,-(0.5 * total / (NSigma[0]*NSigma[0])) + MaxLike_));
	}

	return NSamples;
}

LivePointType		NestSamplSimpleGaussian::do_GetNextSample(int& NlikeEval,double MinLikeValue)
{
	MyArrType	temp(MaxDims_);
	NlikeEval = 0;
	const double r1(NSigma[0] * std::sqrt(-2.0*(MinLikeValue - MaxLike_)));
	double		r;
	const double limit(PriorMaxRange_ / 2.0);
	double LikeValue;
	int bad;

	do
	{
		bad = 0;
		r = r1;
		r *= std::pow(RandomGen_.RandDouble(), 1.0 / static_cast<double>(MaxDims_));

		for(int i = 0; i < MaxDims_; ++i )
		{
			temp[i] = RandomGen_.RandGauss();
		}

		double rr(0.0);

		for(int i = 0; i < MaxDims_; ++i)
		{
			rr += (temp[i] * temp[i]);
		}

		// Scale onto new radius
		double s(r / std::sqrt(rr));

		for(int i = 0;i<MaxDims_;++i)
		{
			temp[i] *= s;
			if(std::abs(temp[i]) > limit)
			{
				bad = 1;
				break;
			}
		}
		if(!bad)
			LikeValue = -(0.5 * r * r / (NSigma[0]*NSigma[0])) + MaxLike_;
		++NlikeEval;
	}while(bad);

	return LivePointType(temp,temp,LikeValue);
}




