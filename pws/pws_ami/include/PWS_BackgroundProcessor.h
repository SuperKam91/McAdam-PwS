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
#ifndef BACKGROUNDPROCESSORH
#define BACKGROUNDPROCESSORH

#include "ZEUS_GaussianRandomGen.h"
#include "PWS_Globals.h"
#include "PWS_GlobalInfoStore.h"

//----------------------------

class BackgroundProcessor
{
public:
	 BackgroundProcessor(int PatchNumber,PatchProcessor& PP)
		:PatchNumber_(PatchNumber),GlbVars_((PlanckInfo::Instance())->GetGlobalVars()),
		PatchProcessor_(PP)
	{}
//
	void Initialize(void);
//
	void MakeWhiteningStructs(void);
//
	inline MapsCollFourierType WhiteFourierMaps(const MapsCollFourierType& OrgMaps)
	{
		return Zeus::Multiply(Loki::Type2Type<Zeus::AlgVect<std::complex<double> > >(),WhiteningStruct_,OrgMaps,Zeus::UB_NOUSE);
	}
//
	inline void ReleaseWhiteningStructs(void)
	{
		WhiteningStruct_.Release();
	}
//
	inline MaskMapCollType&	GetBackgMapColl(void)
	{
		return BackgMaps_;
	}
//
private:
//
	void							MakeCxPowerMaps(void);
//
	void							Make1ModeCrossPower(int off,int YSz,int metric);
//
	void							ReadBackgroundMaps(MaskMapCollType&	BackgMaps);
	void							Read1BackgroundMapData(int patchN, const ObsFreqsType& freq, Zeus::MapType mapT, Zeus::LArr2D<double>& ws) const;
//
	inline int						ModeIndexer(int fst,int sec) const
	{return (WhiteModeSz_ * sec -((sec * sec - sec)/2)) + (fst - sec);}
//
	inline void 	AdjustMapsOrderIndex(void)
	{
		BackgFourMapCollType::iterator			piv(BackMapsColl_.begin());
		BackgFourMapCollType::const_iterator	const end(BackMapsColl_.end());
		for(int i=0;piv != end;++piv,++i)
		{
			piv->FreqOrd_ = i;
		}
	}
//
	inline	double	GetAvSamplesPS_smoothing(double FWHMpix)
	{
		// FWHMpix is in pix units 2.91 = 5.0 arcmin / PlanckpixNSide2048 
		double tSamples(PWSAVSAMPLES);
		if(FWHMpix < 2.91) return tSamples;
		return	(tSamples * 2.91) / FWHMpix;
	}
//
	inline void						MakeWhStruct(void)
	{
		int YSz,XSz;
		WhiteningStruct_.GetSz(YSz,XSz);
		const int	EndOff(YSz * XSz);
		int			CurrOff(0);
		for(;CurrOff != EndOff;++CurrOff)
		{
			Make1ModeCrossPower(CurrOff,YSz,XSz);
		}
	}

	const GlobalScalarVarsType&	GlbVars_; 
	PatchProcessor&				PatchProcessor_;
	int							PatchNumber_;
	int							WhiteModeSz_;

	Zeus::LArr2D<double>		CholWrkSpaceRect_;	// buffer for Cholesky inversion
	BackgFourMapCollType		BackMapsColl_;		// Collection for holding the CxPower matrix raw data
	WhiteStructType				WhiteningStruct_;	// Whitening struct
	MaskMapCollType				BackgMaps_;			// Background maps
	Zeus::PwSCoreRandGen		*NoiGen_;
};

template<typename T>
class InversionFunct
{
public:
	InversionFunct(int WhiteModeSz,int PatchNumber,Zeus::LArr2D<T>& WrkSp,const BackgroundProcessor& Background)
		:WhiteModeSz_(WhiteModeSz),PatchNumber_(PatchNumber),WrkSp_(WrkSp),Background_(Background)
	{}
	void operator() (Zeus::AlgTriMatrix<T>& CrossP,int Inbound)
	{
		const int	Metric(CrossP.GetMetric());
		int			i;

		T	Average;
		Zeus::AlgTriMatrix<T>	TriRes(Metric,false);
		Zeus::AlgTriMatrix<T>	CrossPClone(CrossP.Clone());
		for(i=0;i<2;++i)
		{
			Tri2Rect(CrossPClone,WrkSp_);
			Average = 1.0 / EvalAverageEingValue(WrkSp_);
			RectMult(WrkSp_,Average);
			if(CholInve(WrkSp_))
			{
				if(!i)
				{
					PutNotDiagToZero(CrossPClone);
					continue;
				}
				err(ERRCOD_PWS_INVCROSSPOWER,ERRMSG_PWS_INVCROSSPOWER);
			}
			RectMult(WrkSp_,std::sqrt(Average));
			break;
		}
		Rect2Tri(WrkSp_,TriRes);
		CrossP = TriRes;
	}
private:
	typedef Zeus::AlgTriMatrix<T>	ArgType;
	inline  void					err(int errCode,const wchar_t* msg) const
	{
		std::wstring errstring(msg);
		errstring += std::wstring(PATCHNUMBERSTR);
		errstring += Zeus::PutNumber2Txt(PatchNumber_);
		throw Zeus::libException(errCode,errstring,*this);
	}
	void	Tri2Rect(const Zeus::AlgTriMatrix<T>& Tri ,Zeus::LArr2D<T>& Rect);
	void	Rect2Tri(const Zeus::LArr2D<T>& Rect, Zeus::AlgTriMatrix<T>& Tri);
	void	InvDiag(const Zeus::AlgTriMatrix<T>& CrossP,Zeus::AlgTriMatrix<T>& Inv);
	int		CholInve(Zeus::LArr2D<T>& data);
	T		EvalAverageEingValue(Zeus::LArr2D<T>& Rect);
	void	RectMult(Zeus::LArr2D<T>& Rect,const T Cte);
	void	PutNotDiagToZero(Zeus::AlgTriMatrix<T>& mat);

	int								WhiteModeSz_;
	int								PatchNumber_;
	Zeus::LArr2D<T>&				WrkSp_;
	const BackgroundProcessor&		Background_;
};


#endif //BACKGROUNDPROCESSORH

