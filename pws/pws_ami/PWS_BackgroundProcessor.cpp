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


//----------------------------------
#include <memory>
#include "PWS_BackgroundProcessor.h"
#include "ZEUS_InOut.h"
#include "ZEUS_WorkSpace.h"
#include "ZEUS_Debug.h"
#include "PWS_PatchProcessor.h"

//----------------------------------
//
void	BackgroundProcessor::Initialize(void)
{
	WhiteModeSz_	= 	GlbVars_.N_ObsPlanes_;
	NoiGen_			=   PatchProcessor_.GetNoiseDevice();
}
//
void	BackgroundProcessor::ReadBackgroundMaps(MaskMapCollType& BackgMaps)
{
	ObsFreqsCollType::const_iterator	piv(GlbVars_.FreqsColl_.begin());
	ObsFreqsCollType::const_iterator	const end(GlbVars_.FreqsColl_.end());
	double								Dummy;

	(Zeus::ConManager::Instance())->PrintStr2Console(STARTREADINGBACKGMAPS);

	for(;piv != end;++piv)
	{
		MaskMapType		tMaskBackg(piv->freq_);

		tMaskBackg.ws_.MakeNewSz(GlbVars_.PatchSz_,false); //Needs initialisation
		Read1BackgroundMapData(PatchNumber_, *piv, static_cast<Zeus::MapType>(GlbVars_.N_PriorPlanes_ ? Zeus::MAPTYPE_BACKGROUND : Zeus::MAPTYPE_OBSERVATION), tMaskBackg.ws_.GetInnerData());
		tMaskBackg.ws_.RemoveBias(Dummy);
		BackgMaps.push_back(tMaskBackg);
	}

	(Zeus::ConManager::Instance())->PrintStr2Console(ENDREADINGBACKGMAPS);

}
//
void	BackgroundProcessor::MakeCxPowerMaps(void)
{

	BackMapsColl_.clear();
	BackgMaps_.clear();

	ReadBackgroundMaps(BackgMaps_);

	MaskMapCollType::const_iterator			pivObsFreq(BackgMaps_.begin());
	MaskMapCollType::const_iterator			const endObsFreq(BackgMaps_.end());
	BackgFourMap							tBackgMap;
	MaskMapCollType							WorkBackgMaps;

	for(;pivObsFreq != endObsFreq;++pivObsFreq)
	{
		tBackgMap.Freq_		= pivObsFreq->Freq_;
		tBackgMap.FreqOrd_	= GetFreqOrder(pivObsFreq->Freq_);
		{
			// copy the map
			MaskMapType	Map(tBackgMap.Freq_,pivObsFreq->ws_.Clone());
			// 
			if((PatchProcessor_.DetectedSrsEmpty() || !(GlbVars_.SubSources_)))
			{
				(PlanckInfo::Instance())->TruncateOutliers(Map.ws_,MASK_ALL);
			}
			else
			{
				PatchProcessor_.SubSrcFromMap(Map);
				(PlanckInfo::Instance())->TruncateOutliers(Map.ws_,MASK_BORDER);
			}

			// Apodize maps
			Map.ws_.Transform(PatchProcessor_.GetApodisingDevice(1),0,Zeus::UB_NOUSE);

			tBackgMap.ws_	= PatchProcessor_.Real2Fourier(Map.ws_);
		}
		tBackgMap.Org_	= tBackgMap.ws_.GetInnerData().begin();
		BackMapsColl_.push_back(tBackgMap);
	}
}
//
void	BackgroundProcessor::Read1BackgroundMapData(int patchN, const ObsFreqsType& freq, Zeus::MapType mapT, Zeus::LArr2D<double>& ws) const
{
	std::wstring	Directory;
	std::wstring	FName(PatchProcessor_.GetCurrPatchName(patchN, freq.sign_, mapT, Directory));

	std::auto_ptr<Zeus::GenCollReader<Zeus::LArr2D<double> > >
		FReader(Zeus::GetWrkSpFileReaderHandler(Loki::Type2Type<Zeus::LArr2D<double> >(),
		GlbVars_.Data2Buffer_?1000:GlbVars_.ContextID_,FName,GlbVars_.PatchSz_,GlbVars_.PatchSz_,Directory));
	FReader->Initialize();
	FReader->Read();
	FReader->Release(ws);
}
//
void	BackgroundProcessor::Make1ModeCrossPower(int off,int YSz,int metric)
{
	Zeus::AlgTriMatrix<double> d(WhiteModeSz_,true);
	const double NormCte(static_cast<double>(GlbVars_.PatchSz_ * GlbVars_.PatchSz_ ));
	int	ymode(off / metric);
	if(ymode >  (YSz >> 1)) ymode -= YSz;
	const int	xmode(off % metric);
	double	l(std::sqrt(static_cast<double>(ymode*ymode) + static_cast<double>(xmode*xmode)));
	std::complex<double>		mode0,mode1;
	double		Result;

	l *=  (PITIMES2 / (2.0 * static_cast<double>(metric - 1) * GlbVars_.PixSz_));

	Zeus::AlgTriMatrix<double>::DataInnerType::iterator	orgMatx(d.GetInnerData().begin());

	BackgFourMapCollType::const_iterator	piv(BackMapsColl_.begin());
	BackgFourMapCollType::const_iterator	pivAux2;
	BackgFourMapCollType::const_iterator	const end(BackMapsColl_.end());

	for(;piv != end;++piv)
	{
		mode0 = *(piv->Org_ + off);
		for(pivAux2 = piv;pivAux2 != end;++pivAux2)
		{
			mode1 = *(pivAux2->Org_ + off);
//		Debug code; remove correlation between channels
#ifdef REMOVECORRELATION
			if(pivAux2 == piv)
			{
				Result = Zeus::Cx_XPower(mode1,mode0) / NormCte;
			}
			else
			{Result = 0.0;}
#else
			Result = Zeus::Cx_XPower(mode1,mode0) / NormCte;
#endif
			*(orgMatx + ModeIndexer(pivAux2->FreqOrd_,piv->FreqOrd_)) = Result;
		}
	}

	(WhiteningStruct_.GetInnerData().begin() + off)->Swap(d);
}
//
void 	BackgroundProcessor::MakeWhiteningStructs(void)
{
	MakeCxPowerMaps();
	std::sort(BackMapsColl_.begin(),BackMapsColl_.end());
	AdjustMapsOrderIndex();
	WhiteningStruct_.Release();
	WhiteningStruct_.MakeNewSz(GlbVars_.PatchSz_ >> 1,false);
	CholWrkSpaceRect_.Make(WhiteModeSz_*WhiteModeSz_,WhiteModeSz_,0.0);
	MakeWhStruct();
	BackMapsColl_.clear();
	
	WhiteningStruct_.AzAv(GetAvSamplesPS_smoothing(PatchProcessor_.ZoneGetAverageFWHMpix()),(GlbVars_.PriorMassMin_ < 0.0),Zeus::UB_NOUSE);
	InversionFunct<double>	DummyGcc(WhiteModeSz_,PatchNumber_,CholWrkSpaceRect_,*this);
	WhiteningStruct_.Transform(DummyGcc,Zeus::UB_NOUSE);
}
//
template<typename T>
void  	InversionFunct<T>::PutNotDiagToZero(Zeus::AlgTriMatrix<T>& mat)
{
	Zeus::AlgTriMatrix<double>::DataInnerType::iterator			pivSrc(mat.GetInnerData().begin());
	Zeus::AlgTriMatrix<double>::DataInnerType::const_iterator	const endSrc(pivSrc + mat.GetInnerData().getSz());

	int		Metric(mat.GetMetric());
	int		i(0),j(0),k(0);

	for(;pivSrc != endSrc;++pivSrc,++j)
	{
		if(j == i)
		{
			i += (Metric - k);
			++k;
		}
		else
		{*pivSrc = ((T)0.0);}
	}

}
//
template<typename T>
int  	InversionFunct<T>::CholInve(Zeus::LArr2D<T>& data)
{

	int							i,j,k;
	const int					LinOffset(WhiteModeSz_);
	typename Zeus::LArr2D<T>::iterator	const dataOrg(data.begin());
	Zeus::LArr1D<T>				t_Diagonal(WhiteModeSz_);
	typename Zeus::LArr1D<T>::iterator	const pivDiagonal(t_Diagonal.begin());
	T							t_sum;

	for(i=0;i<WhiteModeSz_;++i)
	{
		for(j=i;j<WhiteModeSz_;++j)
		{
			for(k=i-1,t_sum = *(dataOrg + (i*LinOffset) + j);k>=0;--k)
			{
				t_sum -= (*(dataOrg + (i*LinOffset) + k) * *(dataOrg + (j*LinOffset) + k));
			}
			if(i == j)
			{
				if(t_sum <= 0.0)
				{
					return 1;
				}
				*(pivDiagonal + i) = std::sqrt(t_sum);
			}
			else
			{*(dataOrg + (j*LinOffset) + i) = t_sum / *(pivDiagonal + i);}
		}
	}

	for(i=0;i<WhiteModeSz_;++i)
	{
		*(dataOrg + (i*LinOffset) + i) = ((T)1.0) / *(pivDiagonal + i);
		for(j=i+1,t_sum = ((T)0.0);j<WhiteModeSz_;t_sum = ((T)0.0),++j)
		{
			for(k=i;k<j;++k)
			{t_sum -= (*(dataOrg + (j*LinOffset) + k) * *(dataOrg + (k*LinOffset) + i));}
			*(dataOrg + (j*LinOffset) + i) = t_sum / *(pivDiagonal + j);
		}
	}
	return 0;
}
//
template<typename T>
void   	InversionFunct<T>::InvDiag(const Zeus::AlgTriMatrix<T>& CrossP,Zeus::AlgTriMatrix<T>& Inv)
{
	typename ArgType::DataInnerType::iterator		orgDst(Inv.GetInnerData().begin());
	typename ArgType::DataInnerType::const_iterator	orgSrc(CrossP.GetInnerData().begin());
	typename ArgType::DataInnerType::const_iterator	const endSrc(CrossP.GetInnerData().end());
	for(;orgSrc != endSrc;++orgSrc,++orgDst)
	{
		if((*orgSrc) <= 0.0)
			err(ERRCOD_PWS_INVCROSSPOWER,ERRMSG_PWS_INVCROSSPOWER);
		*orgDst = static_cast<T>((static_cast<T>(1.0) / std::sqrt(*orgSrc)));
	}
}
//
template<typename T>
void InversionFunct<T>::Tri2Rect(const Zeus::AlgTriMatrix<T>& Tri ,Zeus::LArr2D<T>& Rect)
{
	int ArrMetric(static_cast<int>(Rect.getPtrMetric())),i,j;
	Rect.reset();
	typename Zeus::AlgTriMatrix<T>::DataInnerType::const_iterator	pivOrg(Tri.GetInnerData().begin());
	typename Zeus::AlgTriMatrix<T>::DataInnerType::const_iterator	const endOrg(pivOrg + Tri.GetInnerData().getSz());
	typename Zeus::LArr2D<T>::iterator								pivDest(Rect.begin());
	for(i=0,j= 1;pivOrg != endOrg;++i,++pivOrg,++pivDest)
	{
		if(i && !(i % ArrMetric))
		{
			pivDest += j;
			i=j;++j;
		}
		*pivDest = *pivOrg;
	}
}
//
template<typename T>
void InversionFunct<T>::Rect2Tri(const Zeus::LArr2D<T>& Rect, Zeus::AlgTriMatrix<T>& Tri)
{
	int ArrMetric(static_cast<int>(Rect.getPtrMetric())),i,j;

	typename Zeus::AlgTriMatrix<T>::DataInnerType::iterator			pivTri(Tri.GetInnerData().begin());
	typename Zeus::AlgTriMatrix<T>::DataInnerType::iterator			endTri(pivTri + Tri.GetInnerData().getSz());
	typename Zeus::LArr2D<T>::const_iterator						pivRect(Rect.begin());
	for(j=0;j!=ArrMetric;++j){
		for(i=j;i!=ArrMetric;++i)
		{
			*pivTri = *(pivRect + (i*ArrMetric + j));
			++pivTri;
		}
	}
}
//
template<typename T>
T	 InversionFunct<T>::EvalAverageEingValue(Zeus::LArr2D<T>& Rect)
{
	T	Average(0.0);
	int ArrMetric(static_cast<int>(Rect.getPtrMetric())),i;

	typename Zeus::LArr2D<T>::iterator		pivRect(Rect.begin());
	for(i=0;i<ArrMetric;++i,pivRect += (ArrMetric + 1))
	{
		Average += *pivRect;
	}
	
	return Average / static_cast<double>(ArrMetric);
}
//
template<typename T>
void InversionFunct<T>::RectMult(Zeus::LArr2D<T>& Rect,const T Cte)
{
	typename Zeus::LArr2D<T>::iterator			pivRect(Rect.begin());
	typename Zeus::LArr2D<T>::const_iterator		endRect(pivRect + Rect.getSz());

	for(;pivRect != endRect; ++pivRect)
	{*pivRect *= Cte;}
}
//
template void InversionFunct<double>::Tri2Rect(const Zeus::AlgTriMatrix<double>& Tri ,Zeus::LArr2D<double>& Rect);
template void InversionFunct<double>::Rect2Tri(const Zeus::LArr2D<double>& Rect, Zeus::AlgTriMatrix<double>& Tri);
template void InversionFunct<double>::InvDiag(const Zeus::AlgTriMatrix<double>& CrossP,Zeus::AlgTriMatrix<double>& Inv);
template int  InversionFunct<double>::CholInve(Zeus::LArr2D<double>& data);
template double InversionFunct<double>::EvalAverageEingValue(Zeus::LArr2D<double>& Rect);
template void InversionFunct<double>::RectMult(Zeus::LArr2D<double>& Rect,const double Cte);
template void InversionFunct<double>::PutNotDiagToZero(Zeus::AlgTriMatrix<double>& mat);

