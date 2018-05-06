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


//----------------------------
#include "ZEUS_WorkSpace.h"
#include "MC_Antenna.h"
//----------------------------


void	MC_antenna::Initialise(void)
{

	double	freqY,freqX,temp;
	do_Initialise(AntProps_);
	SamplingFactor_	 = 1.0/(AntProps_.PeriodArc_ * static_cast<double>(AntProps_.NSamples_ << 1));

	if(Zeus::PlaneBoundsType::GetDefaultUseBounds())
	{
		do_GetAntennaBounds(freqY,freqX,Zeus::PlaneBoundsType::GetDefaultBoundingFactor());

		temp = freqY / SamplingFactor_;
		if(temp <= 0.0){ AntProps_.BoundSampleY_ = 0;}
		else
		{
			if(temp >= static_cast<double>(AntProps_.NSamples_-1))
			{AntProps_.BoundSampleY_ = Zeus::PlaneBoundsType::INVALID_BOUND;}
			else{AntProps_.BoundSampleY_ = Zeus::toInt(temp + 1.5);}
		}
		temp = freqX / SamplingFactor_;
		if(temp <= 0.0){ AntProps_.BoundSampleX_ = 0;}
		else
		{
			if(temp >= static_cast<double>(AntProps_.NSamples_-1))
			{AntProps_.BoundSampleX_ = Zeus::PlaneBoundsType::INVALID_BOUND;}
			else{AntProps_.BoundSampleX_ = Zeus::toInt(temp + 1.5);}
		}
	}
	else
	{
		AntProps_.BoundSampleY_ = Zeus::PlaneBoundsType::INVALID_BOUND;
		AntProps_.BoundSampleX_ = Zeus::PlaneBoundsType::INVALID_BOUND;
	}

	double Sum,SumSq;

	Zeus::FourierPlane<double>
		tBuffer(AntProps_.NSamples_,Zeus::PlaneBoundsType(AntProps_.BoundSampleY_,AntProps_.BoundSampleX_) ,true);
	
	FillEvalAntennaValues(tBuffer,Sum,SumSq);

	if(AntProps_.Gain_ > 0.0)
	{
		double temp(AntProps_.Gain_ / Sum);
		AntProps_.IntSq_		=  SumSq * temp* temp;
		Zeus::MultFunctor<double>	DummyGcc(temp);
		tBuffer.Transform(DummyGcc,Zeus::UB_NOUSE);
		AntProps_.Int_			= AntProps_.Gain_;
	}
	else
	{
		AntProps_.IntSq_		= SumSq;
		AntProps_.Int_			= Sum;
		AntProps_.Gain_			= 1.0;
	}

	AntProps_.ZeroFourMode_ = (tBuffer.GetInnerData())(0,0);
	buffer_.Swap(tBuffer);
}

void 	MC_antenna::FillEvalAntennaValues(Zeus::FourierPlane<double>& surf,double& Sum,double& SumSq)
{
	const int tNSamples(AntProps_.NSamples_ + 1);

	for(int i=0;i!= tNSamples;++i)
	{
		for(int j=0;j != tNSamples;++j)
		{(surf.GetInnerData())(i,j) =
		 do_GetAntennaAtFreq(static_cast<double>(i)*SamplingFactor_,static_cast<double>(j)*SamplingFactor_) / AmplNormCte_;}
	}

	const int temp(-tNSamples + 1);
	const int temp1((tNSamples - 1)<<1);
	for(int i= -1;i!=temp;--i)
	{
		for(int j=0; j!= tNSamples;++j)
		{
			(surf.GetInnerData())(i + temp1,j) =
			do_GetAntennaAtFreq(static_cast<double>(i)*SamplingFactor_,static_cast<double>(j)*SamplingFactor_) / AmplNormCte_;
		}
	}

	if(AntProps_.Gain_ < 0.0 )
	{
		Zeus::MultFunctor<double>	DummyGcc(1.0 /(surf.GetInnerData())(0,0));
		surf.Transform(DummyGcc,Zeus::UB_NOUSE);
	}
	surf.IntegrateSurface(Zeus::SEL_SUMSUMSQ,Sum,SumSq);

}

