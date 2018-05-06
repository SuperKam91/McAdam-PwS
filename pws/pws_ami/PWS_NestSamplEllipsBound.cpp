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


#include "PWS_NestSamplEllipsBound.h"


double				NestSamplEllipBound::do_Initialise(void)
{
	double MaxLikelihood(do_Reset(PriorSampleNDims_,LikeSampleNDims_));
	wUnitSphereSample_.resize(PriorSampleNDims_);
	wSamples_.Make(PriorSampleNDims_);
	wLikeSample_.Make(LikeSampleNDims_);
	NLikeEvals_		= 0;
	NSuccLikeEvals_ = NLivePoints_;
	lnVolumeCte_	= ((static_cast<double>(PriorSampleNDims_)/2.0)*LNPI) - GammaFn.gammln((static_cast<double>(PriorSampleNDims_)/2.0) + 1.0);
	return MaxLikelihood;
}


void				NestSamplEllipBound::GetSampleFromUnitSphere(TempArrType& Sampl)
{
	double rr;

	do
	{
		rr = 0.0;
		for(int i = 0; i < PriorSampleNDims_; ++i )
		{
			Sampl[i] = RandomGen_.RandGauss();
			rr += (Sampl[i] * Sampl[i]);
		}
	}
	while(rr == 0.0);

	const double	Radius(std::pow(RandomGen_.RandDouble(), 1.0 / static_cast<double>(PriorSampleNDims_)));
	const double	s(Radius/std::sqrt(rr));
	
	for(int i = 0; i < PriorSampleNDims_; ++i )
	{
		Sampl[i] *= s;
	}
}


int					NestSamplEllipBound::Distort2Ellipsoid(const TempArrType& UnitSamples,MyArrType& Samples)
{	
	Samples.reset(0.0);
	for(int j=0;j<PriorSampleNDims_;++j)
	{
		const double Mult(UnitSamples[j]);
		for(int i=0;i<PriorSampleNDims_;++i)
		{
			Samples[i] += (EigenVectors_.element(i,j) * Mult); 
		}
	}

	for(int i=0;i<PriorSampleNDims_;++i)
	{
		Samples[i] += MassCentre_[i];
		if((Samples[i]< 0.0) || (Samples[i]>=1.0))
		{
			return 0;
		}
	}
	return 1;
}

double				NestSamplEllipBound::GetEnlargmentFactor(NEWMAT::DiagonalMatrix& EigenValues,SamplesCollType& CentredSamples)
{
	NEWMAT::ArrayLengthSpecifier	Dummy(PriorSampleNDims_);
	NEWMAT::DiagonalMatrix			InvEigen(Dummy);
	NEWMAT::Matrix					InCorr(PriorSampleNDims_,PriorSampleNDims_);
	double							Accumul;
	double							MaxDist(-1.0e+50);

	for(int i=0;i<PriorSampleNDims_;++i)
	{
		if(EigenValues.element(i) < 1.0e-18)
			EigenValues.element(i) = 1.0e-18;
		InvEigen.element(i) = 1.0 / EigenValues.element(i);
	}

	InCorr = EigenVectors_ * (InvEigen * EigenVectors_.t());

	SamplesCollType::const_iterator	piv(CentredSamples.begin());
	SamplesCollType::const_iterator	const end(CentredSamples.end());

	for(;piv!=end;++piv)
	{
		Accumul = 0.0;
		const TempArrType& Sample(*piv);
		for(int i=0;i<PriorSampleNDims_;++i)
		{
			for(int j=0;j<PriorSampleNDims_;++j)
			{
				Accumul += (InCorr.element(i,j) * Sample[i] * Sample[j]);
			}		
		}
		if(MaxDist < Accumul) MaxDist = Accumul;
	}
	if(MaxDist <= 0.0) MaxDist = 1.0e-18;
	return std::sqrt(MaxDist);
}

void				NestSamplEllipBound::BuildEigenValuesVectors(void)
{
	SamplesCollType					CentredSamples;
	double							DimensionsRatio;
	NEWMAT::ArrayLengthSpecifier	Dummy(PriorSampleNDims_);
	NEWMAT::SymmetricMatrix			Cov(Dummy);
	NEWMAT::DiagonalMatrix			EigenValues;
	const LivePointsSetType&		Lv(GetLP_Coll());

	CentredSamples.clear();
	CentredSamples.resize(Lv.size());

	BuidCovMatrix(CentredSamples,Cov);
	
	NEWMAT::EigenValues(Cov,EigenValues,EigenVectors_);

	const double EnlFactor(GetEnlargmentFactor(EigenValues,CentredSamples));

	for(int i=0;i<PriorSampleNDims_;++i)
	{
		EigenValues.element(i) = std::sqrt(EigenValues.element(i)) * EnlFactor;
	}

	lnAvPriorVolume_	= Get_lnAvPriorVolume();
	lnEllipsoidVolume_	= lnEllipsoidVolume(EigenValues);
	
	if(lnAvPriorVolume_ > lnEllipsoidVolume_)
	{
		DimensionsRatio	= std::exp((lnAvPriorVolume_ - lnEllipsoidVolume_)/static_cast<double>(PriorSampleNDims_)) * xtraEnlFactor_;
	}
	else
	{
		DimensionsRatio	= xtraEnlFactor_;
	}


	for(int j=0;j<PriorSampleNDims_;++j)
	{
		for(int i=0;i<PriorSampleNDims_;++i)
		{
			EigenVectors_.element(i,j) *= (EigenValues.element(j) * DimensionsRatio); 
		}
	}
	NLikeEvals_ = NSuccLikeEvals_ = 0;
}

void				NestSamplEllipBound::BuidCovMatrix(SamplesCollType& CentredSamples,NEWMAT::SymmetricMatrix& Cov)
{
	Cov = 0.0;
	CenterSamples(CentredSamples);
	SamplesCollType::const_iterator			piv(CentredSamples.begin());
	SamplesCollType::const_iterator	const	end(CentredSamples.end());
	double	NFactor(static_cast<double>(CentredSamples.size()));

	for(;piv != end;++piv)
	{
		const TempArrType&	pV(*piv);
		for(int i=0;i<PriorSampleNDims_;++i)
		{
			for(int j=i;j<PriorSampleNDims_;++j)
			{
				Cov.element(i,j) += (pV[i] * pV[j]);
			}
		}
	}

	for(int i=0;i<PriorSampleNDims_;++i)
	{
		for(int j=i;j<PriorSampleNDims_;++j)
		{
			Cov.element(i,j) /= NFactor;
		}
	}
}

void				NestSamplEllipBound::CenterSamples(SamplesCollType& Samples)
{

	const LivePointsSetType&		Lv(GetLP_Coll());

	LivePointAverageVal();

	LivePointsSetType::const_iterator			pivLP(Lv.begin());
	LivePointsSetType::const_iterator	const	endLP(Lv.end());
	SamplesCollType::iterator					pivSampl(Samples.begin());

	for(;pivLP != endLP;++pivLP,++pivSampl)
	{
		pivSampl->resize(PriorSampleNDims_);
		const	double* const	LPValues(pivLP->GetDimValues());
		TempArrType&			SamplVal(*pivSampl);
		for(int i=0;i < PriorSampleNDims_;++i)
		{
			SamplVal[i] = LPValues[i] - MassCentre_[i];
		}
	}
}

void				NestSamplEllipBound::LivePointAverageVal(void)
{
	const LivePointsSetType& Lv(GetLP_Coll());

	double	N(static_cast<double>(Lv.size()));

	MassCentre_.clear();
	MassCentre_.resize(PriorSampleNDims_,0.0);

	LivePointsSetType::const_iterator			piv(Lv.begin());
	LivePointsSetType::const_iterator	const	end(Lv.end());

	for(;piv != end;++piv)
	{
		const double* const v(piv->GetDimValues());
		for(int i=0;i<PriorSampleNDims_;++i)
		{
			MassCentre_[i] += *(v + i);
		}
	}
	for(int j=0;j<PriorSampleNDims_;++j)
	{
		MassCentre_[j] /= N;
	}
}

int					NestSamplEllipBound::do_GetInitialSamples(int NSamples,LivePointsSetType& Values)
{
	MyArrType		PriorSample(PriorSampleNDims_);
	MyArrType		LikeSample(LikeSampleNDims_);
	double			Val;

	Values.clear();
	for(int j=0;j<NSamples;++j)
	{
		do
		{
			do
			{
				for(int i=0;i<PriorSampleNDims_;++i)
				{
					PriorSample[i] = RandomGen_.RandDouble();
				}
			}
			while(!do_TranslateSample(LikeSample,PriorSample));
		}
		while((Val=do_LikelihoodEval(LikeSample)) == Zeus::logZERO);
		Values.push_back(LivePointType(PriorSample,LikeSample,Val));
	}
	return NSamples;
}


LivePointType		NestSamplEllipBound::do_GetNextSample(int& NlikeEval,double MinLikeValue)
{
	double	LikeValue;
	int		NSamples(0);
	NlikeEval = 0;
	if(NeedRebound())
	{
		BuildEigenValuesVectors();
	}

	do
	{
		do
		{
			GetSampleFromUnitSphere(wUnitSphereSample_);
			if(Distort2Ellipsoid(wUnitSphereSample_,wSamples_) && do_TranslateSample(wLikeSample_,wSamples_))
			{break;}
			if(++NSamples > 1000000)
			{
				NlikeEval  = -1;
				return LivePointType(wSamples_,wLikeSample_,Zeus::logZERO);			
			}
		}while(true);

		if(++NlikeEval >= 500000)
		{
			NlikeEval  = -1;
			return LivePointType(wSamples_,wLikeSample_,Zeus::logZERO);
		}
	}while((LikeValue = do_LikelihoodEval(wLikeSample_)) < MinLikeValue);
	NLikeEvals_			+= NlikeEval;
	++NSuccLikeEvals_;
	return LivePointType(wSamples_,wLikeSample_,LikeValue);
}

int					NestSamplEllipBound::NeedRebound(void)
{
	if(NSuccLikeEvals_ >= NLivePoints_)
		return 1;
	if((NLikeEvals_ >= NLivePoints_) && ((static_cast<double>(NSuccLikeEvals_)/static_cast<double>(NLikeEvals_)) < PERCENTREBOUND))
		return 1;
	return 0;
}

