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



//---------------------------------------

#include "ZEUS_Object.h"

//---------------------------------------

namespace Zeus
{

void SrcObject::CheckBoundaries(void)
{
	if(Zeus::PlaneBoundsType::GetDefaultUseBounds())
	{
		double YBound,XBound,temp;
		do_GetObjBounds(YBound,XBound);
		temp = YBound / ObjProps_.PeriodArc_;
		if(temp <= 0.0){ ObjProps_.BoundSampleY_ = 1;}
		else
		{
			if(temp > static_cast<double>(ObjProps_.NSamples_ >> 1) )
			{ObjProps_.BoundSampleY_ = (ObjProps_.NSamples_ >> 1) + 1;}
			else{ObjProps_.BoundSampleY_ = Zeus::toInt(temp + 1.5);}
		}
		temp = XBound / ObjProps_.PeriodArc_;
		if(temp <= 0.0){ ObjProps_.BoundSampleX_ = 1;}
		else
		{
			if(temp > static_cast<double>(ObjProps_.NSamples_ >> 1))
			{ObjProps_.BoundSampleX_ = (ObjProps_.NSamples_ >> 1) + 1;}
			else{ObjProps_.BoundSampleX_ = Zeus::toInt(temp + 1.5);}
		}
	}
	else
	{
		ObjProps_.BoundSampleY_ = Zeus::PlaneBoundsType::INVALID_BOUND;
		ObjProps_.BoundSampleX_ = Zeus::PlaneBoundsType::INVALID_BOUND;
	}
}

void SrcObject::FillPixSurface(bool NoDuocimation)
{
	int tDuoNewSz;
	CheckBoundaries();
	if(	!NoDuocimation &&
		(ObjProps_.BoundSampleY_ != Zeus::PlaneBoundsType::INVALID_BOUND) &&
		(ObjProps_.BoundSampleX_ != Zeus::PlaneBoundsType::INVALID_BOUND) )
	{
		tDuoNewSz = EvalDuocimation(ObjProps_.BoundSampleY_,ObjProps_.BoundSampleX_,ObjProps_.NSamples_);
	}
	else{tDuoNewSz = ObjProps_.NSamples_;}

	Zeus::RealPlane<double>
		tBuffer(tDuoNewSz,Zeus::PlaneBoundsType(ObjProps_.BoundSampleY_,ObjProps_.BoundSampleX_) ,true);
	GetValues(NoDuocimation,tBuffer);
	double Sum,SumSq;
	tBuffer.IntegrateSurface(Zeus::SEL_SUMSUMSQ,Sum,SumSq);
	if(Sum != 1.0)
	{
		Zeus::MultFunctor<double>	DummyGcc(1.0/Sum);
		tBuffer.Transform(DummyGcc,0);
	}
	ObjProps_.SumSq_ = SumSq /(Sum * Sum);
	ObjProps_.Mode0_ = (tBuffer.GetInnerData())(0,0);
	buffer_.Swap(tBuffer);
}

void	SrcObject::GetValues(bool NoDuocimation, Zeus::RealPlane<double>& buff)
{
	Zeus::RealPlane<double>::DataInnerType& DataArray(buff.GetInnerData());
	double		ValueLimit;
	int			StartY,StartX,EndY,EndX;
	const	int HasBounds(buff.HasBounds());
	int			tY,tX,BoundY,BoundX;

	buff.GetSz(tY,tX);

	const	int BuffSz(tY);
	const	int tHalf((BuffSz >> 1) + 1);

	if(HasBounds)
	{
		buff.GetBounds(BoundY,BoundX);
		StartY	= ((BoundY != tHalf) ? -BoundY + 1 : -BoundY + 2);
		EndY	= BoundY;
		StartX	= ((BoundX != tHalf) ? -BoundX + 1 : -BoundX + 2);
		EndX	= BoundX;
		ValueLimit = do_GetObjAtCoord(0.0,0.0) * 1.0e-6;
	}
	else
	{
		StartY	= -tHalf + 2;
		EndY	= tHalf;
		StartX	= -tHalf + 2;
		EndX	= tHalf;
	}

	int			i,j,iAddr,F_XBound(-1),F_YBound(-1);
	double		FunctValue;

	for(i=StartY;i!= EndY;++i)
	{
		iAddr = (i<0?i+BuffSz:i);
		for(j=StartX;j != EndX;++j)
		{
			DataArray(iAddr,(j<0?j+BuffSz:j)) = FunctValue =
			do_GetObjAtCoord((static_cast<double>(i) + ObjProps_.YShift_) * ObjProps_.PeriodArc_,(static_cast<double>(j) + ObjProps_.XShift_) * ObjProps_.PeriodArc_);
			if(HasBounds && (FunctValue > ValueLimit))
			{
				if(std::abs(i) > F_YBound) F_YBound = std::abs(i);
				if(std::abs(j) > F_XBound) F_XBound = std::abs(j);
			}
		}
	}

	if(HasBounds)
	{
		F_YBound			= Zeus::toInt(static_cast<double>(F_YBound) + ObjProps_.YShift_ + 1.5);
		F_XBound			= Zeus::toInt(static_cast<double>(F_XBound) + ObjProps_.XShift_ + 1.5);

		if(F_YBound > tHalf)
		{F_YBound = tHalf;}
		if(F_XBound > tHalf )
		{F_XBound = tHalf;}

		ObjProps_.BoundSampleY_ = F_YBound;
		ObjProps_.BoundSampleX_ = F_XBound;
		buff.SetBounds(F_YBound,F_XBound);
/*
		if((NoDuocimation == false) && (F_YBound != tHalf) && (F_XBound != tHalf))
		{
			const int tSamples(EvalDuocimation(F_YBound,F_XBound,BuffSz));
			if(tSamples < BuffSz)
			{
				Zeus::RealPlane<double>
					tBuffer(tSamples,true);
				Zeus::CopyValuesFrom(tBuffer,buff);
				buff.Swap(tBuffer);
				return;
			}
		}
*/
		if((F_YBound < BoundY) || (F_XBound < BoundX))
		{
			buff.CleanUp(Zeus::PlaneBoundsType(BoundY,BoundX));
		}

		//Put new boundaries and clean 
		return;
	}
}

}

