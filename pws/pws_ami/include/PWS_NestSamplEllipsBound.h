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
#ifndef PSNAKES_NESTSAMPLELLIPSBOUNDH
#define PSNAKES_NESTSAMPLELLIPSBOUNDH

#include "newmatap.h"                // need matrix applications
#include "PWS_NestedSampler.h"

// EXP(-p) = p; p / 2
#define	PERCENTREBOUND		0.2835716451

class NestSamplEllipBound : public NestedSamplerBase
{
//	----------------- Interface --------------
	virtual double				do_Initialise(void);
	virtual int					do_GetInitialSamples(int NSamples,LivePointsSetType& Values);
	virtual LivePointType		do_GetNextSample(int& NlikeEval,double MinLikeValue);
protected:
//  -------------------------------------------
	typedef	Zeus::LArr1D<double>		MyArrType;
	typedef	std::vector<double>			TempArrType;
	typedef	std::vector<TempArrType>	SamplesCollType;


	// Must return the number of prior dimensions and Maximum of the likelihood;
	virtual double				do_Reset(int& NPriorDim,int& NLikeDim)=0;
	// Return 0 if the parameters produce an out-of-range likelihood
	virtual int					do_TranslateSample(MyArrType& LikelihoodSample,const MyArrType& PriorSample)=0;
	virtual double				do_LikelihoodEval(const MyArrType& LikelihoodSample)=0;
//  -------------------------------------------

public:
	NestSamplEllipBound(int NLivePoints,int EvNIndivSamples,int MaxSteps,double FractTolEv,double xtraEnlFactor,int seed)
		:NestedSamplerBase(NLivePoints,EvNIndivSamples,MaxSteps,FractTolEv,seed),
		GammaFn(),NLivePoints_(NLivePoints),xtraEnlFactor_(xtraEnlFactor)
	{}
	virtual	~NestSamplEllipBound(void)
	{}
private:
	inline double	lnEllipsoidVolume(const NEWMAT::DiagonalMatrix& EigValues) const
	{
		double LnV(lnVolumeCte_);
		for(int i=0;i<PriorSampleNDims_;++i)
		{
			LnV += std::log(EigValues.element(i));
		}
		return LnV;
	}

	void			GetSampleFromUnitSphere(TempArrType& Sampl);
	int				Distort2Ellipsoid(const TempArrType& UnitSamples,MyArrType& Samples);
	void			CenterSamples(SamplesCollType& Samples);
	void			LivePointAverageVal(void);
	int				NeedRebound(void);

	// ------------ TODO ---------------------
	int		CheckReboundNeeded(void);
	void	BuidCovMatrix(SamplesCollType& Samples,NEWMAT::SymmetricMatrix& Cov);

	void	BuildEigenValuesVectors(void);
	double	GetEnlargmentFactor(NEWMAT::DiagonalMatrix& EigenValues,SamplesCollType& CentredSamples);

	Zeus::GammaFuncts			GammaFn;
	int							PriorSampleNDims_;
	int							LikeSampleNDims_;
	const int					NLivePoints_;
	const double				xtraEnlFactor_;
	int							NLikeEvals_;
	int							NSuccLikeEvals_;
	NEWMAT::Matrix				EigenVectors_;
	TempArrType					MassCentre_;
	TempArrType					wUnitSphereSample_;
	MyArrType					wSamples_;
	MyArrType					wLikeSample_;
	double						lnAvPriorVolume_;
	double						lnEllipsoidVolume_;
	double						lnVolumeCte_;
};



#endif //PSNAKES_NESTSAMPLELLIPSBOUNDH

