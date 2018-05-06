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



//------------------------------------------------
#ifndef SOURCEOBJECTSZH
#define SOURCEOBJECTSZH

//------------------------------------------------

#include "ZEUS_Object.h"

//-----------------------------------------

namespace Zeus
{

class SrcObjSZ: public SrcObject
{
protected:
	inline virtual	void					do_Initialise(const SrcObjectPropsType& props,double boundTol)
	{
		boundTol_	= boundTol;
		Rc_			= props.ObjParams_.RealScale_;
		Rv_			= Rc_ * RcRvRatio_;
	}

	inline virtual	double					do_GetObjAtCoord(double YCoord,double XCoord) const
	{
		YCoord -= YShift_ ; XCoord -= XShift_;

		if(Rc_ < OBJEPS)
		{
			if((std::abs(YCoord) < PeriodArc_) && (std::abs(XCoord) < PeriodArc_))
				return 1.0;
			else return 0.0;
		}
		double d2(YCoord*YCoord+XCoord*XCoord);
		return (Rc_ * Rv_ / (Rv_ - Rc_)) * ((1.0/std::sqrt(Rc_*Rc_ + d2))   - (1.0/std::sqrt(Rv_*Rv_ + d2)));
	}

	inline virtual	void					do_GetObjBounds(double& YBound,double& XBound) const
	{
		XBound = YBound = (Rc_ * std::sqrt(1 - boundTol_ )) / std::sqrt(boundTol_);

		if(YBound < PeriodArc_)
		{
			YBound = XBound = PeriodArc_;
		}

		YBound += std::abs(YShift_);
		XBound += std::abs(XShift_);
	}
	inline virtual	void					do_ChangeObjParams(const Zeus::ObjFilterRealParams& objP,double YShift,double XShift)
	{
		YShift_	= YShift;XShift_ = XShift;
		Rc_			= objP.RealScale_;
		Rv_			= Rc_ * RcRvRatio_;
	}
public:
	SrcObjSZ(const std::wstring& name,double PeriodArc,int NSamples,double RcRvRatio)
		:SrcObject(name,PeriodArc,NSamples),PeriodArc_(PeriodArc),YShift_(0.0),XShift_(0.0),
		RcRvRatio_(RcRvRatio)
	{}
	virtual ~SrcObjSZ(void)
	{}
private:
	double					boundTol_;
	double					PeriodArc_;
	double					YShift_;
	double					XShift_;
	double					Rc_;
	double					Rv_;
	double					RcRvRatio_;
};

} // Namespace Zeus

#endif //SOURCEOBJECTSZH

