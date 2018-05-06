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

#include "ZEUS_InOut.h"
#include "ZEUS_Debug.h"
#include "PWS_NestedSampler.h"


void	NestedSamplerBase::Initialise(void)
{
	NTotalLikeEval_		= 0;	
    Glb_lnEv_			= Zeus::logZERO;
	Glb_lnLastLike_		= Zeus::logZERO;
	Glb_lnWg_			= Init_lnWg_;
	NestSamplCnt_.clear();
	NestSamplCnt_.resize(NIndivSamples_);
	LivePointsSet_.clear();
	Samples_.clear();
	MaxLoglike_			= do_Initialise();
}

void	NestedSamplerBase::GetInitialSamples(void)
{
	LivePointsSetType	InitValues;
	NTotalLikeEval_ += do_GetInitialSamples(NLivePoints_,InitValues);
	
	LivePointsSetType::iterator				piv(InitValues.begin());
	LivePointsSetType::const_iterator		const end(InitValues.end());
	for(int i=0;piv != end;++i,++piv)
	{
		if(piv->LogLike_ > MaxLoglike_)
			MaxLoglike_ = piv->LogLike_;
		piv->ext_ = RandomGen_.RandDouble();
		LivePointsSet_.push_back(*piv);
	}
	std::make_heap(LivePointsSet_.begin(),LivePointsSet_.end(),std::greater<LivePointType>());
}

int		NestedSamplerBase::GetLogEvidence(StatProps& Ev, double CorrFactor)
{
	std::wstring	ExcptString;
	GetInitialSamples();
	// Start of the main loop
	int		i(0);
	int		result(1);
	int		secs;
	time_t	InitialTime(time(NULL));

	for(;i != MaxSteps_ ;++i)
	{
		SamplePriorVolume();

		if(result = CheckConvergence())
		{
			result = 1 - result;
			break;
		}

		if((result = GetNext_lnLikeSample()) < 0)
		{break;}

		secs = (static_cast<int>(time(NULL) - InitialTime));

		if(secs > 300) // 5 minutes
		{
			result = -1;
			break;		
		}
		else
		{
			if(secs > 60)
			{
				ComputeEvidence(Ev);
				if((Ev.Mean_ + CorrFactor) <= 0.0)
				{
					result = -1;
					break;		
				}
			}
		}
	}

	if((i == MaxSteps_) || (result > 0))
	{
		result = -1;	
	}

	AddRemainingPoints();
	ComputeEvidence(Ev);

	return result;
}

void	NestedSamplerBase::SamplePriorVolume(void)
{
	const double	lnMinLike(PeekAtLowestLP().LogLike_);
	double			t;
	double			LikeValue;


	LikeValue	= Zeus::AddLog(Glb_lnLastLike_,lnMinLike);
	LikeValue	-= LOG2;
	Glb_lnEv_	= Zeus::AddLog(Glb_lnEv_,Glb_lnEvBits_ = (Glb_lnWg_ + LikeValue));

	for(int i=0;i<NIndivSamples_;++i)
	{
		t = std::log(RandomGen_.RandDouble()) / static_cast<double>(NLivePoints_);
		NestSamplCnt_[i].EvidIndiv_ = Zeus::AddLog(NestSamplCnt_[i].EvidIndiv_,std::log(1.0 - std::exp(t)) + NestSamplCnt_[i].t_Indiv_ + LikeValue);
		NestSamplCnt_[i].t_Indiv_	+= t;		
	}

	Glb_lnWg_		-= Glb_lnPriorVolInc_;
	Glb_lnLastLike_	= lnMinLike;
}

void	NestedSamplerBase::AddRemainingPoints(void)
{
	LivePointsSetType::iterator			piv(LivePointsSet_.begin());
	LivePointsSetType::const_iterator	const end(LivePointsSet_.end());
	const double						lnN(std::log(static_cast<double>(NLivePoints_)));
	const double						pLnWg(Glb_lnWg_ - Init_lnWg_ - lnN);

	for(;piv != end;++piv)
	{
		Glb_lnEv_	= Zeus::AddLog(Glb_lnEv_,piv->LogEvBits_ = (pLnWg + piv->LogLike_));
		for(int i=0;i<NIndivSamples_;++i)
		{
			NestSamplCnt_[i].EvidIndiv_ = Zeus::AddLog(NestSamplCnt_[i].EvidIndiv_, NestSamplCnt_[i].t_Indiv_ - lnN + piv->LogLike_);
		}
		Samples_.push_back(*piv);
	}
	LivePointsSet_.clear();
}

void	NestedSamplerBase::ComputeEvidence(StatProps& Ev)
{
	const int centre(NIndivSamples_ >> 1);
	std::sort(NestSamplCnt_.begin(),NestSamplCnt_.end());
	Ev.Median_	= (NestSamplCnt_[centre].EvidIndiv_ + NestSamplCnt_[centre-1].EvidIndiv_) / 2.0;

	double mean(0.0),var(0.0);
	double wnew,diff;

	for(int i=0;i<NIndivSamples_;++i)
	{
        wnew = 1.0 / (i + 1.0);
		diff = NestSamplCnt_[i].EvidIndiv_ - mean;
        mean += (wnew * diff);
        var = (1.0 - wnew) * (var + wnew * diff * diff);
	}

	Ev.Mean_	= mean;
	Ev.StDev_	= std::sqrt(var);
}

void	NestedSamplerBase::PutLivePointInHeap(const LivePointType& lp)
{
	if(MaxLoglike_ < lp.LogLike_)
	{
		MaxLoglike_ = lp.LogLike_;
	}
	std::pop_heap(LivePointsSet_.begin(),LivePointsSet_.end(),std::greater<LivePointType>());
	LivePointsSetType::iterator endLivPt(LivePointsSet_.end() - 1);
	endLivPt->LogEvBits_ = Glb_lnEvBits_; 
	Samples_.push_back(*endLivPt);
	*endLivPt			= lp;
	endLivPt->ext_		= RandomGen_.RandDouble();
	std::push_heap(LivePointsSet_.begin(),LivePointsSet_.end(),std::greater<LivePointType>());
}
//
void	NestedSamplerBase::GetParamsStatsEx(ParamsStatsCollType& Params)
{
	Params.clear();
	if(Samples_.empty())
		return;

	const int	NParam(Samples_[0].GetNTransDims());

	getValFromLivePoint		ValAccessor;
	getWeightFromLivePoint	WeightAccessor(Glb_lnEv_);

	Params.resize(NParam);
	double t;
	for(int p=0;p<NParam;++p)
	{
		ValAccessor.setIndex(p);
		{
			PwS_HistogramType Hist(Samples_,ValAccessor,WeightAccessor);
			Hist.GetHPDLimits(HPDINTERVAL,Params[p].LowHPD_,Params[p].HighHPD_);
			Params[p].Mean_		= Hist.GetExpectedFunctOverHPD(EXPECTVINTERVAL,value);
			Params[p].Mode_		= Hist.GetExpectedFunctOverHPD(MODEINTERVAL,value);
			t					= Hist.GetExpectedFunctOverHPD(1.0,value);
			Params[p].StDev_	= std::sqrt(Hist.GetExpectedFunctOverHPD(1.0,valueSQ) - (t * t));
		}
	}
}
//
void	NestedSamplerBase::GetLinearCorrelatioParams(double VarRatio,LinearCorrParamsType& Params)
{
	if(Samples_.empty())
		return;

	// 2 -> Radius; 3 -> Flux
	getPairFromLivePoint		ValAccessor(2,3);
	getWeightFromLivePoint		WeightAccessor(Glb_lnEv_);

	PwS_LinearFitType LinearFit(Samples_,ValAccessor,WeightAccessor,VarRatio);

	LinearFit.GetParameters(Params.Ord_,Params.Slope_,Params.Corr_);
	LinearFit.GetParamUncertain(Params.ErrOrd_,Params.ErrSlope_);
}
//
void	NestedSamplerBase::GetParamsStats(ParamsStatsCollType& Params)
{
	Params.clear();
	if(Samples_.empty())
		return;

	//HistDebug				DummyDebug;
	//std::vector<HistDebug>	Debug;

	double			TotProb;
	double			prMode;
	double			Prob;
	const double*	SampleParam;
	const int		NParam(Samples_[0].GetNTransDims());

	Params.resize(NParam);

	LivePointsSetType::const_iterator			pivS;
	LivePointsSetType::const_iterator	const	endS(Samples_.end());

	for(int p=0;p<NParam;++p)
	{
		//Debug.clear();
		prMode					= -1.0;
		Params[p].Mean_			= 0.0;
		Params[p].StDev_		= 0.0;
		TotProb					= 0.0;
		pivS					= Samples_.begin();
		for(;pivS != endS;++pivS)
		{
			SampleParam	= pivS->GetTransDimValues();
			Prob		= std::exp(pivS->LogEvBits_ - Glb_lnEv_);
			if(Prob > prMode)
			{
				prMode				= Prob;
				Params[p].Mode_		= SampleParam[p];
			}
			//DummyDebug.Prob_	= Prob;
			//DummyDebug.x_		= SampleParam[p];
			//Debug.push_back(DummyDebug);
			Params[p].Mean_		+= (SampleParam[p] * Prob);
			Params[p].StDev_	+= (SampleParam[p] * SampleParam[p] * Prob);
			TotProb += Prob;
		}
		//DumpHistogram("Histogram",Debug ,1.0,static_cast<int>(Params[p].Mean_));
		Params[p].Mean_		/= TotProb;
		Params[p].StDev_	/= TotProb;
		Params[p].StDev_ = std::sqrt(Params[p].StDev_ - (Params[p].Mean_ * Params[p].Mean_));
	}
}
