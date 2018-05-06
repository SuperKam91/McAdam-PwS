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
#include "PWS_GlobalInfoStore.h"
#include "PWS_AntGauss.h"
#include "ZEUS_ObjGauss.h"
#include "ZEUS_ObjSZNgai.h"
#include "ZEUS_ObjSZ.h"
#include "ZEUS_Debug.h"
#include "PWS_Zone.h"

//-------------------------------------

Zone::Zone(int ZoneIdx)
	:IsInit_(0),ZoneIdx_(ZoneIdx),GlbVars_((PlanckInfo::Instance())->GetGlobalVars()),
	StaticInfoColl_((PlanckInfo::Instance())->GetPlanckStaticInfo()),SrcObj_(0),CacheTime_(0),
	CacheSz_(0),CacheMaxSz_(static_cast<long long>(GlbVars_.CacheSz_*SRCCACHEMAXSZ))
{
	switch(GlbVars_.ScalesFillerType_)
	{
	case 0:
	case 3:
		{
		LinearFiller<double>	DummyGcc0(GlbVars_.SrcMaxScale_,GlbVars_.N_ScaleBins_);
		ScaleBins_.Init(DummyGcc0);
		}
		break;
	case 1:
		{
		LogFiller<double>	DummyGcc1(GlbVars_.SrcMaxScale_,GlbVars_.ScalesFirstElement_,GlbVars_.N_ScaleBins_);
		ScaleBins_.Init(DummyGcc1);
		}
		break;
	case 2:
		{
		HardCodedFiller<double>	DummyGcc2(GlbVars_.ScalesColl_);
		ScaleBins_.Init(DummyGcc2);
		}
		break;
	}
}

void	Zone::MakeAntennas(void)
{
	ObsFreqsCollType::const_iterator	piv(GlbVars_.FreqsColl_.begin());
	ObsFreqsCollType::const_iterator	const end(GlbVars_.FreqsColl_.end());
	ZoneAntInfoType						tAnt;
	std::wstring						AntName(L"AntennaGaussian");

	for(;piv != end;++piv)
	{
		const PlanckStaticInfoAtomType&	plkFreqInfo(GetPlanckStaticInfoByFreq(piv->freq_));
		if(plkFreqInfo.Freq_ == 0)
			continue;
		tAnt.AntHandle_ = AntennaHandleType(new AntennaGaussian(AntName,GlbVars_.PixSz_ * RAD2ARCMIN,GlbVars_.PatchSz_ >> 1,GlbVars_.SZ_?plkFreqInfo.SZ_Signal_:1.0,plkFreqInfo.Freq_,plkFreqInfo.AntFWHM_,plkFreqInfo.AntFWHM_));
		tAnt.AntHandle_->Initialise();
		tAnt.Org_		= ((tAnt.AntHandle_->GetAntennaBufferRef()).GetInnerData()).begin();
		AntennaColl_.push_back(tAnt);
	}
}

void	Zone::do_Initialise(void)
{
	CacheTime_	= 1;
	CacheSz_	= 0;
	MakeAntennas();
	FModeVectors_.MakeNewSz(GlbVars_.PatchSz_ >> 1,false);
	Put2Vectors();
	ReleaseAntBuffers();
	delete SrcObj_;
	SrcObj_ = MakeObject(GlbVars_.SZ_?GlbVars_.ProfParam_.SZ_Profile_:0);
	SrcInfoType temp;
	temp.surf_.Swap(SrcObj_->GetObjectBufferNonConstRef());
	temp.CentralPixAmpl_	= *(temp.surf_.GetInnerData().begin());
	temp.RealParams_		= Zeus::ObjFilterRealParams();	
	PutIntoCache(ObjFilterParams(),temp);
}

SrcInfoType	Zone::GetSrcObjFromRealParam(const Zeus::ObjFilterRealParams& ObjRealParam,bool NoDuocimation)
{
	LockerType lock(getLocker());

	if(NoDuocimation || !(Zeus::PlaneBoundsType::GetDefaultUseBounds()))
		return MakeNewObjShape(ObjRealParam,0.0,0.0,NoDuocimation);

	ObjFilterParams	tobjParam(TranslateParamsInverse(ObjRealParam));
	
	SrcObjctCacheType::iterator piv(SrcCache_.find(tobjParam));
	if (piv == SrcCache_.end())
	{
		SrcInfoType temp(MakeNewObjShape(ObjRealParam,0.0,0.0,false));
		PutIntoCache(tobjParam,temp);
		return temp;
	}
	piv->second.LastHit_	= CacheTime_;
	return piv->second;
}

SrcInfoType	Zone::GetSrcObj(const ObjFilterParams& objParam,Zeus::ObjFilterRealParams& ObjRealParam,bool NoDuocimation)
{
	// this needs to be optsimized
	// SrcInfoType must build the surface with an antenna clone; for another eon

	LockerType lock(getLocker());
	ObjFilterParams	tobjParam(CorrectBinValues(objParam));
	ObjRealParam	= TranslateParams(tobjParam);

	if(NoDuocimation || !(Zeus::PlaneBoundsType::GetDefaultUseBounds()))
		return MakeNewObjShape(ObjRealParam,0.0,0.0,NoDuocimation);
	SrcObjctCacheType::iterator piv(SrcCache_.find(tobjParam));
	if (piv == SrcCache_.end())
	{
		SrcInfoType temp(MakeNewObjShape(ObjRealParam,0.0,0.0,false));
		PutIntoCache(tobjParam,temp);
		return temp;
	}
	piv->second.LastHit_	= CacheTime_;
	return piv->second;
}

SrcInfoType	Zone::MakeNewObjShape(const Zeus::ObjFilterRealParams& objParam,double YShift,double XShift,bool NoDuocimation)
{
	SrcInfoType temp;
	SrcObj_->ChangeObjParams(objParam,YShift,XShift,NoDuocimation);

	temp.surf_.Swap(SrcObj_->GetObjectBufferNonConstRef());
	temp.CentralPixAmpl_	= *(temp.surf_.GetInnerData().begin());
	temp.RealParams_		= objParam;
	return temp;
}

Zeus::SrcObject	*Zone::MakeObject(int ObjType)
{
	Zeus::SrcObjGauss		*tempG;
	Zeus::SrcObjSZ			*tempSZ;
	Zeus::SrcObjSZNgai		*tempSZNgai;
	int						doProfVar(0);
	switch(ObjType)
	{
	case 0:
		tempG		= new Zeus::SrcObjGauss(L"GaussianObj",GlbVars_.PixSz_ * RAD2ARCMIN,GlbVars_.PatchSz_);
		tempG->Initialise();
		return tempG;
	case 1:
		tempSZ		= new Zeus::SrcObjSZ(L"BetaSZ_Obj",GlbVars_.PixSz_ * RAD2ARCMIN,GlbVars_.PatchSz_,GlbVars_.ProfParam_.VirialRatio_);
		tempSZ->Initialise();
		return tempSZ;
	case 2:
		doProfVar = ((GlbVars_.SZ_ == 1) && (GlbVars_.N_ObsPlanes_ != 1) && (GlbVars_.Jf_Estimator_ == 0) && (GlbVars_.AssessmentKind_ == 0) && GlbVars_.ProfParamVar_.IsInit());
		tempSZNgai = new Zeus::SrcObjSZNgai(L"NgaiSZ_Obj",GlbVars_.PixSz_ * RAD2ARCMIN,GlbVars_.PatchSz_,
			GlbVars_.ProfParam_,GlbVars_.ContextID_,GlbVars_.DirInMasks_,GlbVars_.Sync_ID_,doProfVar?&(GlbVars_.ProfParamVar_):NULL);
		tempSZNgai->Initialise();
		return tempSZNgai;
	default:
		errParOutOfRange();
	}
	return 0;
}

long long	Zone::PurgeCache(void)
{
	const long long						CacheSzLimit(CacheMaxSz_ - (CacheMaxSz_>>2));
	long long							CurrObjSz(0);
	HelperCachePurgeAtomType			Atom;
	HelperCachePurgeType				priorityQueue;
	int									YSz,XSz;
	SrcObjctCacheType::const_iterator	pivMap(SrcCache_.begin());
	SrcObjctCacheType::const_iterator	const endMap(SrcCache_.end());
	for(;pivMap != endMap;++pivMap)
	{
		pivMap->second.surf_.GetSz(YSz,XSz);
		Atom.Sz_		= YSz*XSz;
		Atom.priority_	= static_cast<float>(CacheTime_ - pivMap->second.LastHit_)/static_cast<float>(Atom.Sz_);
		Atom.Index_		= pivMap->first;
		priorityQueue.push_back(Atom);
	}

	std::sort(priorityQueue.begin(),priorityQueue.end());

	HelperCachePurgeType::const_iterator	pivHelper(priorityQueue.begin());
	HelperCachePurgeType::const_iterator	const endHelper(priorityQueue.end());

	for(;(pivHelper != endHelper) && (CacheSz_ > CacheSzLimit);++pivHelper)
	{
		CacheSz_ -= (SrcCache_.erase(pivHelper->Index_) * static_cast<long long>(pivHelper->Sz_ * sizeof(double)));
	}
	return CacheSz_;
}
