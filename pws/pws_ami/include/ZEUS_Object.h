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


//--------------------------------------------

#ifndef	SOURCEOBJECTH
#define SOURCEOBJECTH

#include "ZEUS_WorkSpace.h"
#include "ZEUS_Factory.h"
#include "ZEUS_InOut.h"
//----------------------------------
#define	PROFILEBETABINS					32
#define	PROFILEALPHABINS				4
#define	GNFWALPHAMIN					(1.0510-0.8)
#define	GNFWALPHAMAX					(1.0510+0.8)
#define	GNFWBETAMIN						(5.4905-1.49)
#define	GNFWBETAMAX						(5.4905+1.49)

namespace Zeus
{

//
struct SrcObjectPropsType
{
	double						PeriodArc_;
	int							NSamples_;
	int							BoundSampleY_;
	int							BoundSampleX_;
	double						YShift_;
	double						XShift_;
	Zeus::ObjFilterRealParams	ObjParams_;
	double						SumSq_;
	double						Mode0_;
	SrcObjectPropsType(void)
		:BoundSampleY_(Zeus::PlaneBoundsType::INVALID_BOUND),BoundSampleX_(Zeus::PlaneBoundsType::INVALID_BOUND),
		SumSq_(0.0)
	{}
};
//
class SrcObject: public Zeus::FactoryObjType
{
protected:
	virtual	void					do_Initialise(const SrcObjectPropsType& props,double boundTol) = 0;
	virtual	double					do_GetObjAtCoord(double YCoord,double XCoord) const  = 0;
	virtual	void					do_GetObjBounds(double& YBound,double& XBound) const = 0;
	virtual	void					do_ChangeObjParams(const Zeus::ObjFilterRealParams& objP,double YShift,double XShift) = 0;
public:
	typedef Zeus::RealPlane<double>		BufferType;
	typedef BufferType::DataInnerType	BufferDataType;

	inline virtual void	Initialise(void)
	{
		do_Initialise(ObjProps_,Zeus::PlaneBoundsType::GetDefaultUseBounds()?Zeus::PlaneBoundsType::GetDefaultBoundingFactor():-1.0);
		FillPixSurface(false);
	}

	inline void	ChangeObjParams(const Zeus::ObjFilterRealParams& ObjParams,double YShift,double XShift, bool NoDuocimation)
	{
		ObjProps_.ObjParams_	= ObjParams;
		ObjProps_.YShift_		= YShift;
		ObjProps_.XShift_		= XShift;
		do_ChangeObjParams(ObjParams,YShift,XShift);
		FillPixSurface(NoDuocimation);
	}

	SrcObject(const std::wstring& name,double PeriodArc,int NSamples)
		:Zeus::FactoryObjType(name)
	{
		ObjProps_.PeriodArc_	= PeriodArc;
		ObjProps_.NSamples_		= NSamples;
		ObjProps_.YShift_		= 0.0;
		ObjProps_.XShift_		= 0.0;
	}

	inline double	GetObjAtCoordRaw(double YCoord,double XCoord) const
	{return do_GetObjAtCoord(YCoord, XCoord);}
	inline const SrcObjectPropsType&	GetObjectProp(void) const
	{return ObjProps_;}

	inline const BufferType&			GetObjectBufferRef(void) const
	{return buffer_;}
	inline BufferType&					GetObjectBufferNonConstRef(void)
	{return buffer_;}
	inline BufferType					GetObjectBuffer(void) const
	{return buffer_;}
	inline BufferType					CopyObjectBuffer(void) const
	{return buffer_.Clone();}

	inline void							ReleaseObjectBuffer(void)
	{buffer_.Release();}

	virtual ~SrcObject(void)
	{}
private:

	void						FillPixSurface(bool NoDuocimation);
	void						CheckBoundaries(void);
	void						GetValues(bool NoDuocimation, Zeus::RealPlane<double>& buff);

	inline static int			EvalDuocimation(int BoundSampleY,int BoundSampleX,int NSamples)
	{
		int i;
		--BoundSampleY;--BoundSampleX;
		for(i=NSamples / 2;
			(i >= MINIMUMDUOCIMATION) && (i >= BoundSampleY) && (i >= BoundSampleX);
			i >>= 1) ;
		return i<<2;
	}

	SrcObjectPropsType			ObjProps_;
	BufferType					buffer_;
};

} //Namespace Zeus

#endif //SOURCEOBJECTH

