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



//----------------------------------------------------------------
#ifndef PSNAKES_NESTEDSAMPLERH
#define PSNAKES_NESTEDSAMPLERH

#include <time.h>
#include <queue>
#include "ZEUS_PhysicsMath.h"
#include "ZEUS_GaussianRandomGen.h"
#include "ZEUS_WorkSpace.h"

#define	LOG2			0.69314718055994531
#define NMMAXSTEPS		40000
// Use time 
#define	NMSEED			-1
// PwS sampler 
#define MNSAMPLERID		0

#define HPDINTERVAL		0.95
#define EXPECTVINTERVAL	0.95

#define MODEINTERVAL	0.10
// To make the evaluation more stable

inline	double value(double x)
{return  x;}
//
inline	double valueSQ(double x)
{return  x*x;}
//
struct	StatProps
{
	double	Mean_;
	double	StDev_;
	double  Median_;
};
//
struct	ParamsStatProps
{
	double	Mean_;
	double	Mode_;
	double	StDev_;
	double	LowHPD_;
	double	HighHPD_;

	ParamsStatProps(void)
		:Mean_(Zeus::pwsNaN),Mode_(Zeus::pwsNaN),StDev_(Zeus::pwsNaN),
		LowHPD_(Zeus::pwsNaN),HighHPD_(Zeus::pwsNaN)
	{}
};


struct	NestSampleAtom
{
	double	EvidIndiv_;
	double	t_Indiv_;

	inline NestSampleAtom(void)
		: EvidIndiv_(Zeus::logZERO),t_Indiv_(0)
	{}
/*
	inline NestSampleAtom(double EvidIndiv, double t_Indiv)
		: EvidIndiv_(Zeus::logZERO),t_Indiv_(0)
	{}
*/
	inline bool	operator<(const NestSampleAtom& rhs) const
	{
		return EvidIndiv_ < rhs.EvidIndiv_;
	}
};

typedef	std::vector<ParamsStatProps>	ParamsStatsCollType;
typedef	Zeus::VectTriMatxBase<double>	DimCoordsType;

struct LinearCorrParamsType
{
	double	Slope_;
	double	Ord_;
	double	Corr_;
	double	ErrSlope_;
	double	ErrOrd_;

	LinearCorrParamsType(void)
		:Slope_(-Zeus::INF),Ord_(-Zeus::INF),Corr_(-Zeus::INF),ErrSlope_(-Zeus::INF),ErrOrd_(-Zeus::INF)
	{}
};

class	LivePointType
{
public:

	double			LogLike_;
	double			ext_;
	double			LogEvBits_;

	inline LivePointType(void)
		: LogLike_(Zeus::logZERO),ext_(-1.0),LogEvBits_(Zeus::logZERO)
	{}

	inline LivePointType(const Zeus::LArr1D<double>& dataDim,const Zeus::LArr1D<double>& dataTrans,double logLike)
		: DimCoords_(dataDim),TransCoords_(dataTrans),LogLike_(logLike)
	{}	
	inline int	GetNDims(void) const
	{return DimCoords_.GetSz();}
	inline int	GetNTransDims(void) const
	{return TransCoords_.GetSz();}

	inline double* GetDimValues(void)
	{return DimCoords_.GetInnerData().begin();}

	inline const double* GetDimValues(void) const
	{return DimCoords_.GetInnerData().begin();}

	inline double* GetTransDimValues(void)
	{return TransCoords_.GetInnerData().begin();}

	inline const double* GetTransDimValues(void) const
	{return TransCoords_.GetInnerData().begin();}

	inline bool	operator<(const LivePointType& rhs) const
	{return (LogLike_ < rhs.LogLike_) || ((LogLike_ == rhs.LogLike_) && (ext_ < rhs.ext_));}
	inline bool	operator>(const LivePointType& rhs) const
	{return (LogLike_ > rhs.LogLike_) || ((LogLike_ == rhs.LogLike_) && (ext_ > rhs.ext_));}

private:
	DimCoordsType	DimCoords_;
	DimCoordsType	TransCoords_;
};

typedef std::vector<LivePointType>					LivePointsSetType;

class NestedSamplerBase
{
protected:
	typedef std::vector<NestSampleAtom>					NestSamplCntCollType;

//	----------------- Interface ------------------------
	// Must return the Maximum of the likelihood
	virtual double				do_Initialise(void) = 0;
	virtual int					do_GetInitialSamples(int NSamples,LivePointsSetType& Values) = 0;
	virtual LivePointType		do_GetNextSample(int& NlikeEval,double MinLikeValue) = 0;
// -----------------------------------------------------

	inline const LivePointsSetType& GetLP_Coll(void) const
	{
		return LivePointsSet_;
	}

	inline double					Get_lnAvPriorVolume(void) const
	{
		return Glb_lnWg_ - Init_lnWg_;
	}

	inline	const LivePointType&	PeekAtLowestLP(void) const
	{
		return LivePointsSet_[0];
	}

	Zeus::MultNestRandGen			RandomGen_;

public:
	NestedSamplerBase(int NLivePoints, int EvNIndivSamples, int MaxSteps,double FractTolEv,
		int seed)
		: NLivePoints_(NLivePoints),NIndivSamples_(EvNIndivSamples),MaxSteps_(MaxSteps),
		FractTolEv_(std::log(FractTolEv)),
		Init_lnWg_(std::log(1.0 - std::exp(-1.0 / static_cast<double>(NLivePoints_)))),
		Glb_lnPriorVolInc_(1.0 / static_cast<double>(NLivePoints_)),
		RandomGen_(seed)
	{}

	void	Initialise(void);
	int		GetLogEvidence(StatProps& Ev, double CorrFactor);
	void	GetParamsStats(ParamsStatsCollType& Params);
	void	GetParamsStatsEx(ParamsStatsCollType& Params);
	void	GetLinearCorrelatioParams(double VarRatio,LinearCorrParamsType& Params);
//
	inline  int GetNLikeEval(void) const
	{return NTotalLikeEval_;}
	inline	double GetAvLnEvidence(void) const
	{return Glb_lnEv_;}
	inline	double GetMaxLogLike(void) const
	{
		return MaxLoglike_;
	}
	inline	int	ClearSamples(void)
	{
		int t(static_cast<int>(Samples_.size()));
		Samples_.clear();
		return t;
	}
	virtual ~NestedSamplerBase(void)
	{}
private:
	inline NestedSamplerBase(const NestedSamplerBase& );
	inline NestedSamplerBase& operator=(const NestedSamplerBase& rhs);

	void							PutLivePointInHeap(const LivePointType& lp);

	inline int						CheckConvergence(void) const
	{
		const double tW(Glb_lnWg_ - Init_lnWg_);

		if((MaxLoglike_ + tW) < (FractTolEv_ + Glb_lnEv_))
			return 1;
		if(tW < -42.0)
			return -1;
		return 0;
	}

	void							SamplePriorVolume(void);
	void							GetInitialSamples(void);
	void							AddRemainingPoints(void);
	void							ComputeEvidence(StatProps& Ev);

	inline int						GetNext_lnLikeSample(void)
	{
		int		lev;
		time_t	InitialTime(time(NULL));

		const	double LowestLikeValue(PeekAtLowestLP().LogLike_);
		LivePointType		nextSample(do_GetNextSample(lev,LowestLikeValue));

		if((lev < 0) || (static_cast<int>(time(NULL) - InitialTime) > 30)) // 0.5 minutes
		{return -1;}

		NTotalLikeEval_		+= lev;
		while(nextSample.LogLike_ < LowestLikeValue)
		{

			nextSample			= do_GetNextSample(lev,LowestLikeValue);
			
			if((lev < 0) || (static_cast<int>(time(NULL) - InitialTime) > 30)) // 0.5 minutes
			{return -1;}
			
			NTotalLikeEval_		+= lev;
		}
		PutLivePointInHeap(nextSample);
		return 0;
	}


	const int						NLivePoints_;
	const int						NIndivSamples_;
	const int						MaxSteps_;
	const double					FractTolEv_;
	const double					Init_lnWg_;
	const double					Glb_lnPriorVolInc_;
	double							MaxLoglike_;
	double							Glb_lnLastLike_;
	int								NTotalLikeEval_;
    double							Glb_lnEv_;
	double							Glb_lnWg_;
	double							Glb_lnEvBits_;
	NestSamplCntCollType			NestSamplCnt_;
	LivePointsSetType				Samples_;
	LivePointsSetType				LivePointsSet_;

};
//
struct	getValFromLivePoint
{
	inline getValFromLivePoint(int index=0)
		:index_(index)
	{}
//
	inline int setIndex(int NewIndex)
	{
		int t(index_);
		index_ = NewIndex;
		return t;
	}
//
	inline double	operator()(LivePointsSetType::const_iterator p)
	{return (p->GetTransDimValues())[index_];}

	int index_;
};
//
struct	getPairFromLivePoint
{
	inline getPairFromLivePoint(int x=0,int y=1)
		:x_(x),y_(y)
	{}
//
	inline std::pair<int,int> setIndexes(int x,int y)
	{
		int xOld(x_);
		int yOld(y_);
		x_ = x; y_ = y;
		return std::pair<int,int>(xOld,yOld);
	}
//
	inline std::pair<double,double>	operator()(LivePointsSetType::const_iterator p)
	{
		return std::pair<double,double>((p->GetTransDimValues())[x_],(p->GetTransDimValues())[y_]);
	}

	int x_;
	int y_;
};
//
struct	getWeightFromLivePoint
{
	inline getWeightFromLivePoint(double Glb_lnEv)
		:Glb_lnEv_(Glb_lnEv)
	{}
//
	inline double	operator()(LivePointsSetType::const_iterator p)
	{return std::exp(p->LogEvBits_ - Glb_lnEv_);}

	double Glb_lnEv_;
};

typedef Zeus::HistogramMaker<LivePointType,getValFromLivePoint,getWeightFromLivePoint> PwS_HistogramType;

typedef Zeus::LinearRegression<LivePointType,getPairFromLivePoint,getWeightFromLivePoint> PwS_LinearFitType;

//
#endif //PSNAKES_NESTEDSAMPLERH

