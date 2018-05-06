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
#include "PWS_PatchProcessor.h"
#include "PWS_OddsEval.h"

//--------------------------------------------

//
bool	OddsEval::EvaluateOdds(int Size,const MinArrayType& Params,double& Result)
{
	OddEvalBinParams	OddsParams(Zeus::toInt(Params[0] + 0.5),Zeus::toInt(Params[1] + 0.5),ObjFilterParams(Zeus::toInt(Params[2])));
	double				Amplitude;
	double				FirstResult;
	double				SecondResult;
	
	if(Size > 3)
	{Amplitude = Params[3];}
	else{Amplitude = ODDSVAL_NOAMPLITUDE;}

	if(!CheckBounds(OddsParams.YPix_,OddsParams.YPix_))
		return false;
	double Delta(Params[2] - static_cast<double>(OddsParams.Params_.BinScale_));
	if(!GetOddsResult(OddsParams,Amplitude,FirstResult))
		return false;
	Result = FirstResult;
	if(Delta == 0.0) return true;
	++(OddsParams.Params_.BinScale_);
	if(!GetOddsResult(OddsParams,Amplitude,SecondResult))
		return false;
	Result = Zeus::LinearInterpolation(FirstResult,SecondResult,Delta);
	return true;
}

bool	OddsEval::GetOddsResult(const OddEvalBinParams& BinParam,double Amplitude,double& OutResult)
{
	double SurfCorrelation,ISNR2,LikeRatio;

	OddEvalBinParams	tBinParam(BinParam);
	tBinParam.Params_ = PatchProc_.CorrectBinValues(tBinParam.Params_);		
	OddsCacheType::const_iterator piv(OddsCache_.find(tBinParam));

	if (piv == OddsCache_.end())
	{
		SurfCorrelation	= PatchProc_.EvalSurfCorrel(tBinParam.Params_,tBinParam.YPix_,tBinParam.XPix_,ISNR2);
		if(SurfCorrelation >= std::sqrt(ISNR2))
		{
			OddsCache_.insert(OddsCacheType::value_type(tBinParam,OddsResultType(SurfCorrelation,ISNR2)));
		}
	}
	else
	{
		SurfCorrelation = (piv->second).Correlation_;
		ISNR2			= (piv->second).ISNR2_;
	}
	if(Amplitude != ODDSVAL_NOAMPLITUDE)
	{
		LikeRatio = Amplitude * (SurfCorrelation - ((Amplitude * ISNR2) / 2.0));
	}
	else
	{
		LikeRatio = ((SurfCorrelation * SurfCorrelation) / (2.0 * ISNR2));
	}

	OutResult = LikeRatio;

	return true;
}

