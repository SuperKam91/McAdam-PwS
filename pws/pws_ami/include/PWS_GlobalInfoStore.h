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



//---------------------------------
#ifndef PWSGLOBALINFOSTOREH
#define PWSGLOBALINFOSTOREH
#include <time.h>
#include "PWS_Globals.h"
#include "PWS_Zone.h"
#include "ZEUS_InOut.h"

//---------------------------------

#define	CURRENTSTATFORMAT	L"PN => %05d, Obj => %05d, BObj => %05d, WObj => %05d, t => %05d\n"
#define	NONVALIDSTATFORMAT	L"The current patch was not valid => %05d, t => %05d\n"

struct RemoveObjsInMaskFunctor
{
	inline bool operator()(const Zeus::PeakType& peak) const
	{
		if(peak.UnderGalMask_) return 1;
		return 0;
	}
};


class PlankInfoStore
{
public:
	PlankInfoStore(void)
		:TotalN_Objs_(0),TotalN_BrightObjs_(0),TotalN_WrittenObjs_(0),ObjWriter_(0),ObjWriterExternals_(0),
		TotalN_WrittenObjsExt_(0)
	{}
	void									Initialise(CommandLineArgs& args,time_t InitialTime);
	void									ProcessXtraArgs(CommandLineArgs& args);
//
	inline	const GlobalScalarVarsType&		GetGlobalVars(void) const
	{return GlobalScalarVars_;}
//
	inline int	SetScales2HardCoded(int NBins)
	{
		int s(GlobalScalarVars_.ScalesColl_.size());
		
		GlobalScalarVars_.ScalesColl_.resize(NBins);

		for(int i=0;i<NBins;++i)
		{
			GlobalScalarVars_.ScalesColl_[i]=i;
		}
		return s;
	}
//
	inline int	SetSubtractSources(int v)
	{
		int t(GlobalScalarVars_.SubSources_);

		GlobalScalarVars_.SubSources_ = v;
		
		return t;
	}
//
	inline const PlanckStaticType&			GetPlanckStaticInfo(void) const
	{return StaticInfoColl_;}
//
	inline	const Zeus::PatchGeomType&		GetGeoProps(void) const
	{return PatchGeomInfo_;}
//
	inline ZoneHandleType					GetZone(int PatchN)
	{
		int ZIdx(GetZoneIndex(PatchN));
		ZoneCollType::const_iterator	piv(ZoneColl_.begin());
		ZoneCollType::const_iterator	const end(ZoneColl_.end());
		for(;piv !=  end;++piv)
		{
			if(ZIdx == piv->ZoneIdx_)
			{return piv->ZoneHandle_;}
		}
		ZoneHandleType zptr(new Zone(ZIdx));
		ZoneColl_.push_back(ZoneInfoType(ZIdx,zptr));
		return zptr;
	}
//
	inline void ReleaseData(void)
	{
		ZoneColl_.clear();
	}
//
	inline PlanckStaticInfoAtomType&	GetPlanckStaticInfoByFreq(int freq)
	{
		PlanckStaticType::iterator	piv(StaticInfoColl_.begin());
		PlanckStaticType::iterator	const end(StaticInfoColl_.end());
		for(;(piv != end) && (piv->Freq_ != freq);++piv) ;
		if(piv == end) return *(end-1);
		return *piv;
	}
//
	inline void							SetSZGainTo1(void)
	{
		PlanckStaticType::iterator	piv(StaticInfoColl_.begin());
		PlanckStaticType::iterator	const end(StaticInfoColl_.end());
		for(;(piv != end);++piv)
		{
			piv->SZ_Signal_ = (piv->SZ_Signal_ / std::abs(piv->SZ_Signal_));
		}
	}
//
	inline  void		FlushIntermedCat(void)
	{
		wchar_t	buffer[BUFFERMAXCHAR];

		if(ObjWriter_)
		{
			ObjWriter_->Flush();
			PRINTINTOBUFFERFUNCT
				(buffer,BUFFERMAXCHAR,L"Wrote %d lines into -> %ls\n",TotalN_WrittenObjs_,(ObjWriter_->GetCollID()).c_str());
		}
#if (defined(HFIDMC) || defined(LFIDPC)) && defined(HFIDMC_EXTOBJECTS)
		if(ObjWriterExternals_)
		{
			ObjWriterExternals_->Flush();
			PRINTINTOBUFFERFUNCT
				(buffer,BUFFERMAXCHAR,L"Wrote %d lines into -> %ls\n",TotalN_WrittenObjsExt_,(ObjWriterExternals_->GetCollID()).c_str());
		}
#endif
	}
//
	inline  void		errParameterOutOfRange(const std::wstring& ParamName) const
	{
		std::wstring errstring(ERRMSG_PWS_PARAMOUTOFRANGE);
		errstring += std::wstring(L" -> ");
		errstring += ParamName;
		throw Zeus::libException(ERRCOD_PWS_PARAMOUTOFRANGE,errstring,*this);
	}
//
	int										AppendDetectObjs2File(int PatchN,Zeus::PeakCollType& Peaks);
	int										AppendNonValidObjs2File(int SrcIndex,int PatchN);
	int										DoCompressCatalogue(CommandLineArgs& args,double frequency,double SZ_spectralCte);
	void									TranslatingPix2PointingsLib(int PatchN,int NPixCoords,unsigned int Xpix[],unsigned int Ypix[],float Colatitude[],float Longitude[]);
	void									ReadHealpixMaskFile(void);
	int										RemoveObjectInMask(Zeus::PeakCollType& Peaks);
	void									ComputeStats1Map(const Zeus::RealPlane<double>& ws,double rejectLevel,double& bias,double& rms,int& PixIncluded,double guess = -1.0) const;
//
	inline void								TruncateOutliers(Zeus::RealPlane<double>& BackgMap,MaskingType mask)
	{
		double	bias,rms;
		int		PixIncluded;
		wchar_t buffer[BUFFERMAXCHAR];

		ComputeStats1Map(BackgMap,STATREJECTTHRESHOLD,bias,rms,PixIncluded);
		PRINTINTOBUFFERFUNCT
			(buffer,BUFFERMAXCHAR,L"Now masking ----> Pix -> %04d,Bias -> %6.4g,RMS -> %6.4g\n",
			PixIncluded,
			bias,
			rms
			);
		(Zeus::ConManager::Instance())->PrintStr2Console(buffer);

		doTruncateOutliers(BackgMap,bias,rms,mask);
	}
//
	~PlankInfoStore(void)
	{
		delete ObjWriter_;
		delete ObjWriterExternals_;
	}
private:
	struct ZoneInfoType
	{
		int					ZoneIdx_;
		ZoneHandleType		ZoneHandle_;
		ZoneInfoType(int ZoneIdx,ZoneHandleType Z_ptr)
			:ZoneIdx_(ZoneIdx),ZoneHandle_(Z_ptr)
		{}
	};

	typedef  std::vector<ZoneInfoType>	ZoneCollType;

	void				InitialiseStatic(void);
	void				ChangeStatic(int freq, double& FWHM);
	void				ReadGlobalScalarVars(void);
	void				ReadGeoPropsHeader(void);
	void				ReadMapInfo(void);
	void				doTruncateOutliers(Zeus::RealPlane<double>& BackgReal,double bias,double rms,MaskingType mask);
//
	inline bool 		InbordersCoords(int YCoord,int XCoord) const
	{
		if((YCoord < GlobalScalarVars_.PatchBorder_)							||
			(YCoord >= (GlobalScalarVars_.PatchSz_ - GlobalScalarVars_.PatchBorder_))	||
			(XCoord < GlobalScalarVars_.PatchBorder_)						||
			(XCoord >= (GlobalScalarVars_.PatchSz_ - GlobalScalarVars_.PatchBorder_))
			)
			return true;
		return false;
	}
//
	inline bool						IsPixInBorder(long OffSet) const
	{
		const int	XCoord(OffSet % GlobalScalarVars_.PatchSz_);
		const int	YCoord(OffSet / GlobalScalarVars_.PatchSz_);

		return InbordersCoords(YCoord,XCoord);
	}
//	
	void				TranslatingPix2Pointings(int PatchN,Zeus::PeakCollType& Peaks);
	std::wstring		GetCurrPatchPtgsOutputName(int patchN,int ptgType) const;
	void				ReadPtgs(int patchN,Zeus::LArr2D<double>& CoLat,Zeus::LArr2D<double>& Long) const;
	void				ProcessScalesColl(void);
//
	inline	void		PrintCurrentStatus(int PatchN) const
	{
		wchar_t	buffer[BUFFERMAXCHAR];	
		PRINTINTOBUFFERFUNCT
			(buffer,BUFFERMAXCHAR,CURRENTSTATFORMAT,
			PatchN,
			TotalN_Objs_,
			TotalN_BrightObjs_,
			TotalN_WrittenObjs_,
			static_cast<int>(time(NULL) - InitialTime_)
			);
		(Zeus::ConManager::Instance())->PrintStr2Console(buffer);
	}

	inline  void		errCantFindFile(const std::wstring& fname) const
	{
		std::wstring errstring(ERRMSG_PWS_CANTFINDFILE);
		errstring += std::wstring(L" -> ");
		errstring += GlobalScalarVars_.DirPointings_;errstring += fname;
		throw Zeus::libException(ERRCOD_PWS_CANTFINDFILE,errstring,*this);
	}

	inline  void		errInvalidString(const std::wstring& str) const
	{
		std::wstring errstring(ERRMSG_PWS_INVALIDSTRING);
		errstring += std::wstring(L" -> ");errstring += str;
		throw Zeus::libException(ERRCOD_PWS_INVALIDSTRING,errstring,*this);
	}

	inline int GetZoneIndex(int PatchN)
	{return 0;}
	
	int									ContextID_;
	std::wstring						IntermFName_;
	ZoneCollType						ZoneColl_;
	GlobalScalarVarsType				GlobalScalarVars_;
	Zeus::PatchGeomType					PatchGeomInfo_;
	PlanckStaticType					StaticInfoColl_;
	time_t								InitialTime_;
	int									TotalN_Objs_;
	int									TotalN_BrightObjs_;
	int									TotalN_WrittenObjs_;
	int									TotalN_WrittenObjsExt_;
	Healpix_Map<HEALPIX_ATOM_PREC>		MaskMap_;
	Zeus::GenCollWriter<Zeus::PeakCollType>	*ObjWriter_;
	Zeus::GenCollWriter<Zeus::PeakCollType>	*ObjWriterExternals_;
};


#ifdef MULTI_THREAD
typedef Zeus::Multi_SingletonHolder<PlankInfoStore>::Type	PlanckInfo;
#else
typedef Zeus::Single_SingletonHolder<PlankInfoStore>::Type	PlanckInfo;
#endif

// Filler for the class Bins to map a linear scale into the radius scales

template<typename T>
struct HardCodedFiller
{
	const std::vector<double> t;

	inline HardCodedFiller(const std::vector<double>& Scales)
		:t(Scales)
	{
		(PlanckInfo::Instance())->SetScales2HardCoded(Scales.size());
	}

	inline int	GetNBins(void) const
	{return (t.size()-1);}

	inline T	operator()(int i) const
	{
		return t[i];
	}
};

template<typename T>
struct LogFiller
{
	const T		lnMinV_;
	const T		lnStep_;
	const int	Bins_;

	inline LogFiller(T	MajV,T	MinV,int NBins)
		:lnMinV_(std::log(MinV)),lnStep_(std::log(MajV/MinV)/static_cast<T>(NBins-1)),Bins_(NBins)
	{}

	inline int	GetNBins(void) const
	{return Bins_;}

	inline T	operator()(int i) const
	{
		if(!i) return static_cast<T>(0.0);
		return std::exp(lnMinV_ + lnStep_ * static_cast<T>(i-1));
	}
};

template<typename T>
struct LinearFiller
{
	const T			MajV_;
	const int		NBins_;

	inline LinearFiller(T	MajV,int NBins)
		:MajV_(MajV),NBins_(NBins)
	{}

	inline int	GetNBins(void) const
	{return NBins_;}

	inline T	operator()(int i) const
	{
		return (MajV_ / static_cast<T>(NBins_)) * static_cast<T>(i);		
	}
};

#endif //PWSGLOBALINFOSTOREH

