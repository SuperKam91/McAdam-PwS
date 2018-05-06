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
#ifndef PSNAKES_NESTSAMPLECASESH
#define PSNAKES_NESTSAMPLECASESH

#include "ZEUS_PhysicsMath.h"
#include "ZEUS_Priors.h"
#include "PWS_OddsEval.h"
#include "PWS_GlobalInfoStore.h"
#include "PWS_NestedSampler.h"
#include "PWS_NestSamplEllipsBound.h"

class	NestSamplPwSObserParamNew : public NestSamplEllipBound
{
	static const int			MaxDims_;

	typedef	std::vector<Zeus::Priors*>	PriorsSlotsCollType;

//  -------------------------------------------
	virtual double				do_Reset(int& NPriorDim,int& NLikeDim);
	virtual int					do_TranslateSample(MyArrType& LikelihoodSample,const MyArrType& PriorSample);
	virtual double				do_LikelihoodEval(const MyArrType& LikelihoodSample);
//  -------------------------------------------

	inline  void					errPriorValues(const wchar_t* FunctName,int PatchNumber,double PriorValue) const
	{
		std::wstring errstring(std::wstring(SRCINDEXSTR) + Zeus::PutNumber2Txt(PatchNumber) + std::wstring(L"\n"));
		errstring	+= (std::wstring(ERRMSG_PWS_INVALPRIORVAL) + Zeus::PutNumber2Txt(PriorValue));
		throw Zeus::libException(ERRCOD_PWS_INVALPRIORVAL,errstring,FunctName);
	}
	
	
	inline void	DisposePriors(void)
	{
		PriorsSlotsCollType::iterator	piv(PriorsSlots_.begin());
		PriorsSlotsCollType::iterator	const end(PriorsSlots_.end());
		for(;piv != end;++piv)
		{delete *piv;}
	}
//
	void	SetPosFlatPriorBounds(Zeus::PeakType* peak,const Zeus::PatchGeomLineType& SrcGeom);
	void	SetPosGaussianPriorParams(Zeus::PeakType* peak,const Zeus::PatchGeomLineType& SrcGeom);
	void	SetRadiusGaussianPriorParams(Zeus::PeakType* peak,const Zeus::PatchGeomLineType& SrcGeom);
	void	SetFluxGaussianPriorParams(Zeus::PeakType* peak,const Zeus::PatchGeomLineType& SrcGeom);
	void	SetFluxUniformPriorParams(Zeus::PeakType* peak,const Zeus::PatchGeomLineType& SrcGeom);

public:
	NestSamplPwSObserParamNew(const Zeus::PriorsTypeCollType& PriorsCollT,OddsEval& LikeEval,int NLivePoints,
								int EvNIndivSamples,
								int MaxSteps,double FractTolEv,double xtraEnlFactor,
								int BeamFullRadiusPix,int BeamHalfRadiusPix,int seed=-1);

	void			NextPeak(Zeus::PeakType* peak,const Zeus::PatchGeomLineType& SrcGeom);
	inline	double	GetEvCorrectingCte(void)
	{
		return EvCorrectingCte_;
	}

	virtual		~NestSamplPwSObserParamNew(void)
	{
		DisposePriors();
	}

private:

	PriorsSlotsCollType				PriorsSlots_;
	const GlobalScalarVarsType&		GlbVars_;
	double							EvCorrectingCte_;
	double							FluxConvCte_;
	int								BeamFullRadiusPix_;
	int								BeamHalfRadiusPix_;
	double							InnerCore_;
	Zeus::PeakType*					PeakProps_;
	OddsEval&						LikeEval_;
};

// --------------------------------------------------------
// - Simple example just for testing ellipsoidal bounding -
// --------         Declaration             ---------------

#define	LN2PI			1.8378771
#define	LN2PIOVER2		0.9189385
#define	DIMS			4
#define PRIORRANGE		10.0

extern double	NSigma[DIMS];

class	NestSamplSimplGaussWithBounds : public NestSamplEllipBound
{
//  -------------------------------------------
	virtual double				do_Reset(int& NPriorDim,int& NLikeDim);
	virtual int					do_TranslateSample(MyArrType& LikelihoodSample,const MyArrType& PriorSample);
	virtual double				do_LikelihoodEval(const MyArrType& LikelihoodSample);
//  -------------------------------------------
public:
	NestSamplSimplGaussWithBounds(int NLivePoints,int EvNIndivSamples,int MaxSteps,double FractTolEv,double xtraEnlFactor,int seed=-1)
		: NestSamplEllipBound(NLivePoints,EvNIndivSamples,MaxSteps,FractTolEv,xtraEnlFactor,seed),MaxDims_(DIMS),
		PriorMaxRange_(PRIORRANGE)
	{
		MaxLike_ = ((static_cast<double>(-MaxDims_)/2.0) * LN2PI);
		for(int i=0;i<MaxDims_;++i)
		{
			MaxLike_ -= std::log(NSigma[i]);
		}
	}
	virtual	~NestSamplSimplGaussWithBounds(void)
	{}
private:
	const int			MaxDims_;
	double				MaxLike_;
	const double		PriorMaxRange_;
};


// ---------------------------------------------------
// -------- Simple example just for testing ----------
// --------         Declaration             ----------


class NestSamplSimpleGaussian : public NestedSamplerBase
{
	typedef	Zeus::LArr1D<double>	MyArrType;
protected:

//	----------------- Interface --------------
	virtual double				do_Initialise(void);
	virtual int					do_GetInitialSamples(int NSamples,LivePointsSetType& Values);
	virtual LivePointType		do_GetNextSample(int& NlikeEval,double MinLikeValue);
// -------------------------------------------

public:
	NestSamplSimpleGaussian(int NLivePoints, int NIndivSamples, int MaxSteps,double FractTolEv,int seed=-1)
		: NestedSamplerBase(NLivePoints,NIndivSamples,MaxSteps,FractTolEv,seed),MaxDims_(DIMS),
			PriorMaxRange_(PRIORRANGE),MaxLike_(static_cast<double>(-MaxDims_) * (LN2PIOVER2 + std::log(NSigma[0])))
	{}
	virtual	~NestSamplSimpleGaussian(void)
	{}
private:
	const int			MaxDims_;
	const double		MaxLike_;
	const double		PriorMaxRange_;
	LivePointType		GetNextSample(int& NlikeEval);
};



#endif


