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


//-----------------------------------------
#ifndef	SOURCEOBJECTGAUSSH
#define SOURCEOBJECTGAUSSH

#include "ZEUS_Object.h"

//-----------------------------------------

namespace Zeus
{

class SrcObjGauss: public SrcObject
{
protected:
	inline virtual	void					do_Initialise(const SrcObjectPropsType& props,double boundTol)
	{
		boundTol_	= boundTol;
		filtP_		= props.ObjParams_;
		if((4.0 * props.ObjParams_.RealScale_) <= PeriodArc_) Coef_ = 1.0;
		else Coef_	= -1.0/(2.0 * props.ObjParams_.RealScale_ * props.ObjParams_.RealScale_);
	}
	inline virtual	double			do_GetObjAtCoord(double YCoord,double XCoord) const
	{
		if(Coef_ > 0.0)
		{
			if((std::abs(YCoord) <= PeriodArcOver2_) && (std::abs(XCoord) <= PeriodArcOver2_))
				return 1.0;
			else return 0.0;
		}
		return std::exp(Coef_*(YCoord*YCoord+XCoord*XCoord));
	}
	inline virtual	void			do_GetObjBounds(double& YBound,double& XBound) const
	{
		if(Coef_ > 0.0)
		{
			YBound = XBound = PeriodArc_;
			return;
		}
		YBound = std::sqrt(std::log(boundTol_)/Coef_);
		XBound = YBound;
		YBound += std::abs(YShift_);
		XBound += std::abs(XShift_);
	}
	inline virtual	void			do_ChangeObjParams(const Zeus::ObjFilterRealParams& objP,double YShift,double XShift)
	{
		filtP_	= objP;
		YShift_	= YShift;XShift_ = XShift;
		if((4.0 * objP.RealScale_) <= PeriodArc_) Coef_ = 1.0;
		else Coef_	= -1.0/(2.0 * objP.RealScale_ * objP.RealScale_);
	}
public:
	SrcObjGauss(const std::wstring& name,double PeriodArc,int NSamples)
		:SrcObject(name,PeriodArc,NSamples),PeriodArc_(PeriodArc),PeriodArcOver2_(PeriodArc/2.0),YShift_(0.0),XShift_(0.0)
	{}
	virtual ~SrcObjGauss(void)
	{}
private:
	double					boundTol_;
	double					PeriodArc_;
	double					PeriodArcOver2_;
	double					YShift_;
	double					XShift_;
	double					Coef_;
	Zeus::ObjFilterRealParams		filtP_;
};

} //Namespace Zeus




#endif //SOURCEOBJECTGAUSSH


