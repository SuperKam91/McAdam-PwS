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
#include "PWS_Antenna.h"
#include "ZEUS_Debug.h"
//----------------------------


void	Antenna::Initialise(void)
{
	double	freqY,freqX,temp;
	do_Initialise(AntProps_);
	SamplingFactor_	 = 1.0/(AntProps_.PeriodArc_ * static_cast<double>(AntProps_.NSamples_<< 1));

	if(Zeus::PlaneBoundsType::GetDefaultUseBounds())
	{
		do_GetAntennaBounds(freqY,freqX,Zeus::PlaneBoundsType::GetDefaultBoundingFactor());

		temp = freqY / SamplingFactor_;
		if(temp <= 0.0){ AntProps_.BoundSampleY_ = 1;}
		else
		{
			if(temp > static_cast<double>(AntProps_.NSamples_))
			{AntProps_.BoundSampleY_ = AntProps_.NSamples_ + 1;}
			else{AntProps_.BoundSampleY_ = Zeus::toInt(temp + 1.5);}
		}
		temp = freqX / SamplingFactor_;
		if(temp <= 0.0){ AntProps_.BoundSampleX_ = 1;}
		else
		{
			if(temp > static_cast<double>(AntProps_.NSamples_))
			{AntProps_.BoundSampleX_ = AntProps_.NSamples_ + 1;}
			else{AntProps_.BoundSampleX_ = Zeus::toInt(temp + 1.5);}
		}
	}
	else
	{
		AntProps_.BoundSampleY_ = Zeus::PlaneBoundsType::INVALID_BOUND;
		AntProps_.BoundSampleX_ = Zeus::PlaneBoundsType::INVALID_BOUND;
	}

	std::complex<double> Sum,SumSq;

	Zeus::FourierPlane<std::complex<double> >
		tBuffer(AntProps_.NSamples_,Zeus::PlaneBoundsType(AntProps_.BoundSampleY_,AntProps_.BoundSampleX_) ,true);
	
	FillEvalAntennaValues(tBuffer,Sum,SumSq);

	AntProps_.IntSq_		= SumSq.real();
	AntProps_.Int_			= Sum.real();
	Zeus::MultFunctor<std::complex<double> > Dummy(std::complex<double>(AntProps_.Gain_,0.0));
	//G++ BUG; can't resolve function otherwise !! Grrrrrrrr
	tBuffer.Transform(Dummy,Zeus::UB_NOUSE);
	buffer_.Swap(tBuffer);
}

void 	Antenna::FillEvalAntennaValues(Zeus::FourierPlane<std::complex<double> >& surf,
									   std::complex<double>& Sum,std::complex<double>& SumSq)
{
	Zeus::FourierPlane<std::complex<double> >::DataInnerType& DataArray(surf.GetInnerData());
	double		ValueLimit;
	int			EndY,EndX;
	const	int HasBounds(surf.HasBounds());
	int			tY,tX,BoundY,BoundX;
	surf.GetSz(tY,tX);
	surf.GetBounds(BoundY,BoundX);

	if(HasBounds)
	{
		EndY	= BoundY;
		EndX	= BoundX;
		ValueLimit = std::abs(do_GetAntennaAtFreq(0.0,0.0)) * surf.GetBoundFactor();
	}
	else
	{
		EndY	= (tY >> 1) + 1;
		EndX	= tX;
	}

	int			F_XBound(-1),F_YBound(-1);

	for(int i=0;i!= EndY;++i)
	{
		for(int j=0;j != EndX;++j)
		{
			DataArray(i,j) =
			do_GetAntennaAtFreq(static_cast<double>(i)*SamplingFactor_,static_cast<double>(j)*SamplingFactor_);
			if(HasBounds && (std::abs(DataArray(i,j)) > ValueLimit))
			{
				if(i > F_YBound) F_YBound = i;
				if(j > F_XBound) F_XBound = j;
			}
		}
	}

	const int	temp((EndY	!= ((tY >> 1) + 1)) ? -EndY : -EndY + 1);
	const int	temp1(tY);

	for(int i= -1;i!=temp;--i)
	{
		for(int j=0; j!= EndX;++j)
		{
			DataArray(i + temp1,j) = 
			do_GetAntennaAtFreq(static_cast<double>(i)*SamplingFactor_,static_cast<double>(j)*SamplingFactor_);
			if(HasBounds && (std::abs(DataArray(i + temp1,j)) > ValueLimit))
			{
				if((-i) > F_YBound) F_YBound = -i;
				if(j > F_XBound) F_XBound = j;
			}
		}
	}

	if(HasBounds)
	{
		F_YBound	= Zeus::toInt(static_cast<double>(F_YBound) +  1.5);
		F_XBound	= Zeus::toInt(static_cast<double>(F_XBound) +  1.5);

		if(F_YBound > (tY >> 1) + 1)
		{BoundY = (tY >> 1) + 1;}
		else{BoundY = F_YBound;}
		if(F_XBound > tX)
		{BoundX = tX;}
		else{BoundX = F_XBound;}
		AntProps_.BoundSampleY_ = BoundY;
		AntProps_.BoundSampleY_ = BoundX;
		surf.SetBounds(BoundY,BoundX);
	}

	Zeus::MultFunctor<std::complex<double> > Dummy(std::complex<double>(1.0,0.0) / DataArray(0,0));
	//G++ BUG; can't resolve function otherwise !! Grrrrrrrr
	surf.Transform(Dummy,Zeus::UB_NOUSE);
	surf.IntegrateSurface(Zeus::SEL_SUMSUMSQ,Sum,SumSq,Zeus::UB_NOUSE);
}

