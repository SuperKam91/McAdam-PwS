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


//-------------------------------------
#ifndef PWS_ZONEH
#define PWS_ZONEH

#include "PWS_Globals.h"
#include "PWS_Antenna.h"
#include "ZEUS_Object.h"

//-------------------------------------

class Zone:public HandleRefCounter
{
public:
	Zone(int ZoneIdx);
	inline void	Initialise(void)
	{
		if(!IsInit_)
		{
			LockerType lock(getLocker());
			if(!IsInit_)
			{
				do_Initialise();
				IsInit_ = 1;
			}
		}
	}
	inline const MapsCollFourierType& GetAntVectorsRef(void) const
	{return FModeVectors_;}
	inline MapsCollFourierType GetAntVectorsCopy(void) const
	{return FModeVectors_.Clone();}

	inline Zeus::FourierPlane<std::complex<double> > GetAntenna(int AntennaIndex)
	{
		int YSz,XSz;
		FModeVectors_.GetSz(YSz,XSz);
		Zeus::ExtractSurfFromVectorsFunct<std::complex<double> > Extractor(FModeVectors_,XSz - 1,0/*100 GHz*/);
		FModeVectors_.Transform(Extractor,Zeus::UB_NOUSE);
		return Extractor.GetSurface();
	}

	inline double			GetAntenna0pix(int AntennaIndex)
	{
		return ((AntennaColl_.at(AntennaIndex)).AntHandle_->GetAntenaProp()).Int_;
	}
//
	SrcInfoType		GetSrcObj(const ObjFilterParams& objParam,Zeus::ObjFilterRealParams& ObjRealParam,bool NoDuocimation);
//
	SrcInfoType		GetSrcObjFromRealParam(const Zeus::ObjFilterRealParams& ObjRealParam,bool NoDuocimation);

	inline SrcInfoType		GetSrcObjFromRealParamNoCache(const Zeus::ObjFilterRealParams& ObjRealParam, bool NoDuocimation)
	{
		LockerType lock(getLocker());

		return MakeNewObjShape(ObjRealParam, 0.0, 0.0, NoDuocimation);
	}

//
	inline double	GetAntennasAverageFWHMpix(void)
	{
		AntennaCollType::const_iterator	piv(AntennaColl_.begin());
		AntennaCollType::const_iterator	const end(AntennaColl_.end());
		double		Weights(0.0);
		double		tFWHM(0.0);

		for(;piv != end;++piv)
		{
			const AntennaPropsType& AntPropRef(piv->AntHandle_->GetAntenaProp());	
			tFWHM		+= AntPropRef.FWHM_;
			Weights		+= std::abs(AntPropRef.Gain_);
		}
		return tFWHM / (Weights * (((AntennaColl_.begin())->AntHandle_->GetAntenaProp()).PeriodArc_));
	}

	long long				PurgeCache(void);

	inline SrcInfoType		GetSrcObjRaw(const Zeus::ObjFilterRealParams& objParam)
	{
		LockerType lock(getLocker());
		return MakeNewObjShape(TranslateParams(objParam),0.0,0.0,true);
	}

	inline SrcInfoType		GetObjOfRadiusBin(int RadiusBin)
	{
		Zeus::ObjFilterRealParams	Dummy;
		ObjFilterParams objParam(RadiusBin);
		return  GetSrcObj(objParam,Dummy,true);
	}
//
	inline SrcInfoType		GetObjOfRadiusRealValue(double Radius)
	{
		return  GetSrcObjFromRealParam(Zeus::ObjFilterRealParams(Radius),true);
	}
//
	inline SrcInfoType		GetObjOfParamsRealValue(const Zeus::ObjFilterRealParams& params)
	{
		return  GetSrcObjFromRealParam(params,true);
	}

//
	inline Zeus::ObjFilterRealParams	TranslateParams(const Zeus::ObjFilterRealParams& objParam) const
	{
		return Zeus::ObjFilterRealParams(ScaleBins_.getValue(objParam.RealScale_));
	}

	inline Zeus::ObjFilterRealParams	GetNearestScale(const Zeus::ObjFilterRealParams& objParam) const
	{
		return Zeus::ObjFilterRealParams(ScaleBins_.GetNearestScale(objParam.RealScale_));
	}

	inline Zeus::ObjFilterRealParams	GetNearestScaleButFirst(const Zeus::ObjFilterRealParams& objParam) const
	{
		return Zeus::ObjFilterRealParams(ScaleBins_.GetNearestScaleButFirst(objParam.RealScale_));
	}

	inline Zeus::ObjFilterRealParams	TranslateParamsInverseInterpol(const Zeus::ObjFilterRealParams& objParam) const
	{
		return Zeus::ObjFilterRealParams(ScaleBins_.GetBinInterpol(objParam.RealScale_));
	}

	inline ObjFilterParams	TranslateParamsInverse(const Zeus::ObjFilterRealParams& objParam) const
	{
		return ObjFilterParams(ScaleBins_.GetBin(objParam.RealScale_));
	}

	inline Zeus::ObjFilterRealParams	TranslateParams(const ObjFilterParams& objParam) const
	{
		return Zeus::ObjFilterRealParams(ScaleBins_.getValue(objParam.BinScale_));
	}
	
	inline ObjFilterParams	CorrectBinValues(const ObjFilterParams& objParam) const
	{
		return ObjFilterParams(ScaleBins_.CorrectBin(objParam.BinScale_));
	}
	
	inline void ClearCache(void)
	{
		SrcCache_.clear();
		CacheSz_	= 0;
		CacheTime_	= 1;
	}


	~Zone(void)
	{
		delete SrcObj_;
	}
private:
	Zone(const Zone&);
	Zone& operator=(const Zone&);

	typedef  Zeus::ObjHandle<Antenna>::Type		AntennaHandleType;

	struct ZoneAntInfoType
	{
		Antenna::BufferDataType::const_iterator	Org_;
		AntennaHandleType						AntHandle_;
	};

	struct	HelperCachePurgeAtomType
	{
		ObjFilterParams		Index_;
		float				priority_;
		int					Sz_;
		bool	operator<(const HelperCachePurgeAtomType& rhs) const
		{return priority_ > rhs.priority_;}
	};

	typedef std::vector<HelperCachePurgeAtomType>	HelperCachePurgeType;		


	typedef  std::vector<ZoneAntInfoType>		AntennaCollType;

	void						do_Initialise(void);
	void						MakeAntennas(void);
	SrcInfoType					MakeNewObjShape(const Zeus::ObjFilterRealParams& objParam,double YShift,double XShift,bool NoDuocimation);

	inline const PlanckStaticInfoAtomType&	GetPlanckStaticInfoByFreq(int freq) const
	{
		PlanckStaticType::const_iterator	piv(StaticInfoColl_.begin());
		PlanckStaticType::const_iterator	const end(StaticInfoColl_.end());
		for(;piv != end;++piv)
		{
			if(piv->Freq_ == freq)
			return *piv;
		}
		return *(end-1);
	}

	inline void							Put1FMode2Vector(int mode)
	{
		Zeus::AlgVect<std::complex<double> > d(GlbVars_.N_ObsPlanes_,false);
		Zeus::AlgVect<std::complex<double> >::DataInnerType::iterator	orgVector(d.GetInnerData().begin());
		AntennaCollType::const_iterator piv(AntennaColl_.begin());
		AntennaCollType::const_iterator const end(AntennaColl_.end());
		for(;piv != end;++piv,++orgVector)
		{
			*orgVector = *(piv->Org_ + mode);
		}
		(FModeVectors_.GetInnerData().begin() + mode)->Swap(d);
	}

	inline void							Put2Vectors(void)
	{
		int YSz,XSz;
		FModeVectors_.GetSz(YSz,XSz);
		const int	EndOff(YSz * XSz);
		int			CurrOff(0);
		for(;CurrOff != EndOff;++CurrOff)
		{Put1FMode2Vector(CurrOff);}
	}

	inline void							ReleaseAntBuffers(void)
	{
		AntennaCollType::iterator piv(AntennaColl_.begin());
		AntennaCollType::const_iterator const end(AntennaColl_.end());
		for(;piv != end;++piv)
		{
			piv->AntHandle_->ReleaseAntennaBuffer();
			piv->Org_ = 0;
		}
	}

	inline long long					UpdateCacheSz(const SrcInfoType& Src)
	{
		int	XSz,YSz;
		Src.surf_.GetSz(YSz,XSz);
		CacheSz_	+= (YSz * XSz * sizeof(double));
		return CacheSz_;
	}

	inline long long					PutIntoCache(const ObjFilterParams& ObjParams,SrcInfoType& Src)
	{
		if(CacheSz_ >= CacheMaxSz_)
		{
			PurgeCache();
		}
		Src.LastHit_	= CacheTime_;
		SrcCache_.insert(SrcObjctCacheType::value_type(ObjParams,Src));
		return UpdateCacheSz(Src);
	}
	
	inline  void		errParOutOfRange(void) const
	{
		throw Zeus::libException(ERRCOD_PWS_PARAMOUTOFRANGE,ERRMSG_PWS_PARAMOUTOFRANGE,"Zone::MakeObject");
	}

	
	Zeus::SrcObject						*MakeObject(int ObjType);

	int									ZoneIdx_;

	long long							CacheTime_;
	long long							CacheSz_;
	const	GlobalScalarVarsType&		GlbVars_;
	const	PlanckStaticType&			StaticInfoColl_;
	baseIntVolatile						IsInit_;
	AntennaCollType						AntennaColl_;
	MapsCollFourierType					FModeVectors_;
	Zeus::ValuesBins<double>			ScaleBins_;
	SrcObjctCacheType					SrcCache_;
	Zeus::SrcObject						*SrcObj_;
	const long long						CacheMaxSz_;

};




#endif //PWS_ZONEH

