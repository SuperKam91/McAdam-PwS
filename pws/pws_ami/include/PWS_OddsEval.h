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

#ifndef PWSODDSEVALH
#define PWSODDSEVALH

#define ODDSVAL_NOAMPLITUDE		-1e30

#include "PWS_Globals.h"
#include "PWS_GlobalInfoStore.h"
#include "PWS_PatchProcessor.h"

//--------------------------------------------

class OddsEval
{

public:
	OddsEval(PatchProcessor& PatchProc)
		:GlbVars_((PlanckInfo::Instance())->GetGlobalVars()),PatchProc_(PatchProc),
		OutOfBoundsValue_(ODDSVAL_OUTOFBOUNDS)
	{}

	inline double	OddsValue(int Size,const MinArrayType& Params)
	{
		double Result;
		if(!EvaluateOdds(Size,Params,Result)) 
			return OutOfBoundsValue_;
		if(Result > OutOfBoundsValue_) OutOfBoundsValue_ = Result;
		return  (-Result) ;
	}

	inline double	GetOddsResultRaw(int YPix,int XPix,const Zeus::ObjFilterRealParams& inFilterParams,Zeus::ObjFilterRealParams& outFilterParams,double Amplitude)
	{
		double ISNR2;

		if((inFilterParams.RealScale_ < 0.0) || !CheckBounds(YPix,XPix))
			return Zeus::logZERO;

		const  double SurfCorrelation(PatchProc_.EvalSurfCorrelSample(inFilterParams,outFilterParams,YPix,XPix,ISNR2));

		return Amplitude * (SurfCorrelation - ((Amplitude * ISNR2) / 2.0));
	}

	inline double	GetOddsResultProfVarRaw(int YPix, int XPix, double YOffPix, double XOffPix, const Zeus::ObjFilterRealParams& inFilterParams, bool Noduocimation, double Amplitude)
	{
		double ISNR2;

		if ((inFilterParams.RealScale_ < 0.0) || !CheckBounds(YPix, XPix))
			return Zeus::logZERO;

		const  double SurfCorrelation(PatchProc_.EvalSurfCorrelprofVarSample(inFilterParams, Noduocimation, YPix, XPix, YOffPix, XOffPix, ISNR2));

		return Amplitude * (SurfCorrelation - ((Amplitude * ISNR2) / 2.0));
	}

	inline double	GetOddsResultRawOptimAmp(int YPix,int XPix,const Zeus::ObjFilterRealParams& inFilterParams,double& Isnr,int& sign)
	{
		Zeus::ObjFilterRealParams outFilterParams;

		if((inFilterParams.RealScale_ < 0.0) || !CheckBounds(YPix,XPix))
			return -OutOfBoundsValue_;
		const  double SurfCorrelation(PatchProc_.EvalSurfCorrelSample(inFilterParams,outFilterParams,YPix,XPix,Isnr));
		sign = (SurfCorrelation < 0.0? -1 : 1);
		return  (SurfCorrelation * SurfCorrelation) / (2.0 * Isnr);
	}
//
	inline double	GetOddsResultRawOptimAmp(int YPix, int XPix, const ObjFilterParams& inFilterParams, double& Isnr, int& sign)
	{
		if(!CheckBounds(YPix,XPix))
			return -OutOfBoundsValue_;
		const  double SurfCorrelation(PatchProc_.EvalSurfCorrel(inFilterParams,YPix,XPix,Isnr));
		sign = (SurfCorrelation < 0.0? -1 : 1);
		return  (SurfCorrelation * SurfCorrelation) / (2.0 * Isnr);
	}

	bool			GetOddsResult(const OddEvalBinParams& BinParam,double Amplitude,double& result);

	inline	void	ClearCache(void)
	{
		OddsCache_.clear();
	}
private:
	bool		EvaluateOdds(int Size,const MinArrayType& Params,double& Result);

	inline bool	CheckBounds(int YPix,int XPix) const
	{

		if(	(YPix < GlbVars_.PatchBorder_) || (YPix >= (GlbVars_.PatchSz_ - GlbVars_.PatchBorder_)) ||
			(XPix < GlbVars_.PatchBorder_) || (XPix >= (GlbVars_.PatchSz_ - GlbVars_.PatchBorder_))
			)
			return false;

		return true;
	}

	typedef std::map<OddEvalBinParams, OddsResultType>	OddsCacheType;

	const GlobalScalarVarsType&		GlbVars_;
	PatchProcessor&					PatchProc_;
	double							OutOfBoundsValue_;
	OddsCacheType					OddsCache_;
};

#endif //PWSODDSEVALH

