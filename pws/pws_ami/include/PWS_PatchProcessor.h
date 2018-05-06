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
#ifndef PATCHPROCESSORH
#define PATCHPROCESSORH

#include "PWS_Globals.h"
#include "ZEUS_Priors.h"
#include "PWS_GlobalInfoStore.h"
#include "PWS_BackgroundProcessor.h"
#include "ZEUS_GaussianRandomGen.h"
#include "ZEUS_Debug.h"

//----------------------------


#define EQUIVBEAMSIGMAFACTOR	0.005
#define THRESHOLD10PERC			2.751535313
#define PENALTY1SOURCE			1.38629
#define FILTREPFORMAT			L"Scale -> %3.4f,Sigma -> %3.4f,Isnr -> %3.4f,ratio -> %3.4f\n"
#define PATCHREPFORMAT			L"Freq. -> %04d,Bias -> %6.4e,RMS -> %6.4e\n"

#define PWSMASKLEN				2
#define MAXNOFBRIGHTPSPERPATCH	100

#define CONTPRIOR_XPOS			0
#define CONTPRIOR_YPOS			1
#define CONTPRIOR_THETA			2
#define CONTPRIOR_Y				3
#define CONTPRIOR_SIZE			4

class	NestSamplPwSObserParamOld;
class	NestSamplPwSObserParamNew;


double GetYtotConvY500(double alpha, double beta, double gamma, double c500, double R500Ratio, Zeus::GammaFuncts& Gammafun);

inline bool	LowerFluxesRemoveMark(const MaxLikeSurfMaximum& peak)
{
	return peak.Removed_;
}

struct LikePixBufferType
{
	int						InUse_;
	double					ISNR2_;
	Zeus::LArr2D<double>	pixbuffer;
};

struct IntCatalFilterFunctor
{
	int PatchBorder_,PatchSz_;

	inline bool operator()(const Zeus::PeakType& peak) const
	{
		if(	(peak.Pos_.YPix_ <	PatchBorder_)				||
			(peak.Pos_.YPix_ >= (PatchSz_ - PatchBorder_))	||
			(peak.Pos_.XPix_ <	PatchBorder_)				||
			(peak.Pos_.XPix_ >= (PatchSz_ - PatchBorder_))
			)
			return true;

		return false;
	}

	IntCatalFilterFunctor(int PatchBorder,int PatchSz)
		:PatchBorder_(PatchBorder<<1),PatchSz_(PatchSz)
	{}
};


template<typename T>
struct ApodiseFunctor
{
	int YSz_,XSz_,Direct_,Border_;

	ApodiseFunctor(int YSz,int XSz,int Border)
		:YSz_(YSz),XSz_(XSz),Direct_(1),Border_(Border>>1)
	{}

	inline void		SetMode(int mode)
	{Direct_ = mode;}
	inline  static T	GetPwSpApodCorrect(void)
	{
		return 0.25; //(0.5 * 0.5)
		//return 1.0;
	}
	inline T		ApodiseFunction(int YCoord,int XCoord) const
	{
		if((YCoord < 0) || (YCoord >= YSz_) || (XCoord < 0) || (XCoord >= XSz_))
			return 0.0;
		return (std::sin((PI * static_cast<T>(XCoord) / static_cast<T>(XSz_))) * std::sin((PI * static_cast<T>(YCoord) / static_cast<T>(YSz_))));
		//return 1.0;
	}

	inline void	operator()(T& x,int offset) const
	{
		int XCoord(offset % XSz_);
		int YCoord(offset / XSz_);
		
		if(Direct_)
		{
			x *= ApodiseFunction(YCoord,XCoord);
		}
		else
		{
			if(
				(YCoord < Border_) || ((YSz_ - YCoord) <= Border_) ||
				(XCoord < Border_) || ((XSz_ - XCoord) <= Border_)
				)
				return;
			x /= (std::sin((PI * static_cast<T>(XCoord) / static_cast<T>(XSz_))) * std::sin((PI * static_cast<T>(YCoord) / static_cast<T>(YSz_))));
			return;
		}
	}
};

class SubAddFunctor
{
public:
	SubAddFunctor(bool SubAdd)
		:SubAdd_(SubAdd)
	{}
	double operator() (double patchVal,double surfVal)
	{
		if(SubAdd_) return (surfVal - patchVal);
		return (surfVal + patchVal);
	}
private:
	const bool	SubAdd_;
};

struct PixMasker
{
public:
	PixMasker(void)
	{}
	inline double operator() (double patchVal,double surfVal)
	{
		return (patchVal == MASKVALUE)? surfVal : patchVal;
	}
};

enum  SrcMaskerState{SMASK_IN = -1,SMASK_OUT = +1};
enum  SrcMaskerTransState{SMASK_OUTIN = -1,SMASK_INOUT = 1};


struct MaskBorders
{
	int						y_;
	int						x_;
	int						Out_;
	double					PixVal_;
	SrcMaskerTransState		Trans_;
	inline MaskBorders(void)
	{}
	inline MaskBorders(int y,int x,SrcMaskerTransState Trans)
		:y_(y),x_(x),Trans_(Trans),Out_(0)
	{}
};

inline bool	SortMaskBordersX(const MaskBorders& lhs,const MaskBorders& rhs)
{
	if(lhs.y_ != rhs.y_) return lhs.y_ < rhs.y_;
	if (static_cast<int>(lhs.Trans_) != static_cast<int>(rhs.Trans_))
		return (static_cast<int>(lhs.Trans_) < static_cast<int>(rhs.Trans_));
	return lhs.x_ < rhs.x_;
}

inline bool	SortMaskBordersY(const MaskBorders& lhs,const MaskBorders& rhs)
{
	if(lhs.x_ != rhs.x_) return lhs.x_ < rhs.x_;
	if (static_cast<int>(lhs.Trans_) != static_cast<int>(rhs.Trans_))
		return (static_cast<int>(lhs.Trans_) < static_cast<int>(rhs.Trans_));
	return lhs.y_ < rhs.y_;
}


typedef std::vector<MaskBorders>	MaskBordersCollType;

class SrcMasker
{
public:
	SrcMasker(double Amplitude,double NormalCte,int metric,MaskBordersCollType& MaskBorders,int swapCoords)
		:Amplitude_(Amplitude),Level_(Amplitude * NormalCte * 3.355e-4),
		Metric_(metric),Previous_(-1000000),MaskBorders_(MaskBorders),SwapCoords_(swapCoords)
	{}
	inline void operator() (double& Val,int InBound)
	{
		int x,y;
		CoordsFromOffset(InBound,y,x);

		if(SwapCoords_)
		{
			if(x != Previous_)
			{
				PreviousState_  = SMASK_OUT;
				Previous_		= x;
			}
		}
		else
		{
			if(y != Previous_)
			{
				PreviousState_  = SMASK_OUT;
				Previous_		= y;
			}
		}

		if((Val*Amplitude_) > Level_)
		{
			if(PreviousState_ == SMASK_IN)
				return;
			if(SwapCoords_)
			{y -= (PWSMASKLEN+1);}
			else
			{x -= (PWSMASKLEN+1);}
			MaskBorders_.push_back(MaskBorders(y,x,SMASK_OUTIN));
			PreviousState_ = SMASK_IN;
		}
		else
		{
			if(PreviousState_ == SMASK_OUT)
				return;
			if(SwapCoords_)
			{y += PWSMASKLEN;}
			else
			{x += PWSMASKLEN;}
			MaskBorders_.push_back(MaskBorders(y,x,SMASK_INOUT));
			PreviousState_ = SMASK_OUT;
		}
	}
private:
	inline void CoordsFromOffset(int offset,int& y,int& x) const
	{
		x	= (offset % Metric_);
		y	= (offset / Metric_);
		if(y > (Metric_ >> 1))
			y -= Metric_;
		if(x > (Metric_ >> 1))
			x -= Metric_;
	}
	SrcMaskerState			PreviousState_;
	const double			Amplitude_;
	const double			Level_;
	const int				Metric_;
	int						Previous_;
	int						SwapCoords_;
	MaskBordersCollType&	MaskBorders_;

};

class PatchProcessor
{
public:
	PatchProcessor(int PatchNumber)
		:PatchNumber_(PatchNumber),GlbVars_((PlanckInfo::Instance())->GetGlobalVars()),PatchesGeoInfo_((PlanckInfo::Instance())->GetGeoProps()),
		StaticInfoColl_((PlanckInfo::Instance())->GetPlanckStaticInfo()),NoiseGen_(0),Background_(0),
		OddsEval_(0),BrentMin_(0),PwSnake_(0),NestSampler_(0),CurrentIter_(0),EqBeamSigmaPix_(-1),
		NBrightPsInPatch_(0),
		ApodFunct_(GlbVars_.PatchSz_,GlbVars_.PatchSz_,GlbVars_.PatchBorder_),y0EvalIntegrator_(0)
	{
		ScalesProps_.resize(GlbVars_.ScalesColl_.size());
	}
	void									Initialise(void);
	void									MakeMainSurfaces(double LikeSurfISNR2Cte,double& LikeSurfSigma, int freeMem);
	void									FindSourcePositions(double LikeSurfISNR2Cte,int StoreNoise);
	inline void								GetPatchNoiseProperties(void)
	{
		ObjFilterParams		objParams;
		RealPlaneSurfType	dummySurface;
		double				Dummy1,Dummy2;
		const	int			ScalesCollSz(GlbVars_.ScalesColl_.size());

		for(int i=0;i < ScalesCollSz;++i)
		{
			MakeLikelihoodAtScale(i,dummySurface,Dummy1,Dummy2,objParams,1);
		}
	}
	void									CheckSubSources(std::vector<double>& tGauss);
	void									FindOptimalScaleAtPixel(MaxLikeSurfMaximum& LikePeak);
	PP_FindOptParamsResType					CharacterAndAssess(int DoDetection);	
	std::wstring							GetCurrPatchName(int patchN,int freq,Zeus::MapType mapT,std::wstring& Directory) const;
	void									FindSources(int DoDetection);
	int										ReportCurrentSrcs(int SrcIdex);
	int										SetMultiScaleY(Zeus::PeakCollType& PeakColl);
	void									SwapPatchNoiseEstimates(void);
	void									SubSrcFromMap(MaskMapType& Map);
//--------------------

	inline  double ZoneGetAverageFWHMpix(void)
	{
		return Zone_->GetAntennasAverageFWHMpix();
	}


	inline  void ZoneClearCache(void)
	{
		Zone_->ClearCache();
	}
	
	inline ObjFilterParams					CorrectBinValues(const ObjFilterParams& objParam) const
	{
		return Zone_->CorrectBinValues(objParam);
	}

	inline	double							EvalSurfCorrel(const ObjFilterParams& FilterParams,int YPix,int XPix,double& ISNR2)
	{
		Zeus::ObjFilterRealParams Dummy;
		SrcInfoType		ObjSurf(Zone_->GetSrcObj(FilterParams,Dummy,false));
		ISNR2 = GetISNR(FilterParams,ObjSurf);

		//DumpInOut_2d("ObjSurf",256,256,256,ObjSurf.surf_.GetInnerData().begin(),1.0);
		//DumpInOut_2d("LikelihoodReal",256,256,256,LikelihoodSurf_.GetInnerData().begin(),1.0);

		return Zeus::Correlation(LikelihoodSurf_,ObjSurf.surf_,YPix,XPix);
	}
//
	inline	double							EvalSurfCorrelSample(const Zeus::ObjFilterRealParams& inFilterParams,Zeus::ObjFilterRealParams&	outFilterParams,int YPix,int XPix,double& ISNR2)
	{
		ObjFilterParams FilterParams(Zone_->TranslateParamsInverse(inFilterParams));
		SrcInfoType		ObjSurf(Zone_->GetSrcObj(FilterParams,outFilterParams,false));
		ISNR2 = GetISNR(FilterParams,ObjSurf);
		return Zeus::Correlation(LikelihoodSurf_,ObjSurf.surf_,YPix,XPix);
	}
//
	inline	double		EvalSurfCorrelprofVarSample(const Zeus::ObjFilterRealParams& inFilterParams, bool noduocimation, int YPix, int XPix, double YOffPix, double XOffPix, double& ISNR2)
	{
		SrcInfoType			ObjSurf(Zone_->GetSrcObjFromRealParamNoCache(inFilterParams, noduocimation));
// shift the object
		MapsMainSurfFourType		tempFour(Real2Fourier(ObjSurf.surf_));
		tempFour.ShiftInPlace(YOffPix, XOffPix);
		RealPlaneSurfType tSrc(Fourier2Real(tempFour));
		Zeus::MultFunctor<double>	DummyGcc(1.0 / static_cast<double>(GlbVars_.PatchSz_ * GlbVars_.PatchSz_));
		tSrc.Transform(DummyGcc, 0, Zeus::UB_NOUSE);
		ISNR2 = ISNRFromSurf(tSrc);
		return Zeus::Correlation(LikelihoodSurf_, tSrc, YPix, XPix);
	}
//
	inline  double							GetISNR(const ObjFilterParams& SurfaceParams,const SrcInfoType& Surface)
	{
		ObjFilterParams tFiltParams(Zone_->CorrectBinValues(SurfaceParams));
		PatchFilterPropsCacheType::const_iterator piv(PatchFilterPropsCache_.find(tFiltParams));
		if (piv == PatchFilterPropsCache_.end())
		{
			double tISNR(ISNRFromSurf(Surface.surf_));
			PatchFilterPropsCache_.insert(PatchFilterPropsCacheType::value_type(tFiltParams,PatchFilterPropsType(tISNR)));
			return tISNR;
		}
		return (piv->second).Isnr2_;
	}
//
	inline double							GetISNR(const Zeus::ObjFilterRealParams& FiltParams)
	{
		SrcInfoType		ObjSurf(Zone_->GetSrcObjRaw(FiltParams));
		return ISNRFromSurf(ObjSurf.surf_);
	}

	inline double							ISNRFromSurf(const RealPlaneSurfType& Surface)
	{
		MapsMainSurfFourType		ObjFour(Real2Fourier(Surface));
		Zeus::PowerFunctor<double>	Dummy;
		AntMainSurfType				ObjFour2(Transform(Loki::Type2Type<double>(),ObjFour,Dummy,Zeus::UB_NOUSE));
		RealPlaneSurfType			ObjReal2(Fourier2Real(ObjFour2));
		int YSz,XSz;
		ObjReal2.GetSz(YSz,XSz);
		return Zeus::DotProduct(AntennaSFE_real_,ObjReal2,Zeus::UB_NOUSE) / static_cast<double>(YSz*XSz);
	}

	inline const PlanckStaticInfoAtomType&	GetPlanckStaticInfoByFreq(int freq) const
	{
		PlanckStaticType::const_iterator	piv(StaticInfoColl_.begin());
		PlanckStaticType::const_iterator	const end(StaticInfoColl_.end());
		for(;(piv != end) && (piv->Freq_ != freq);++piv) ;
		if(piv == end) return *(end-1);
		return *piv;
	}

	template<typename T>
	inline int GetFourMachineIndex(const Zeus::DataPlaneSurface<T>& rp)
	{
		int YSz,XSz;
		rp.GetSz(YSz,XSz);
		int	r(GlbVars_.PatchSz_ / YSz);
		int i,duocim;
		for(duocim = 1,i = 0; duocim < r;++i,duocim <<= 1)  ;
		return i;
	}

	inline Zeus::FourierPlane<std::complex<double> > Real2Fourier(const Zeus::RealPlane<double>& rp)
	{
		int MachineIndex(GetFourMachineIndex(rp));
		Zeus::FourierPlane<std::complex<double> >   temp((FourierFarm_.at(MachineIndex))->GetXComplexBufferSz() - 1);
		(FourierFarm_.at(MachineIndex))->Real2Fourier(rp.GetInnerData(),temp.GetInnerData());
		return temp;
	}

	inline Zeus::RealPlane<double>			Fourier2Real(const Zeus::FourierPlane<std::complex<double> >& rp)
	{
		int MachineIndex(GetFourMachineIndex(rp));
		Zeus::RealPlane<double>   temp((FourierFarm_.at(MachineIndex))->GetRealBufferSz());
		(FourierFarm_.at(MachineIndex))->Fourier2Real(rp.GetInnerData(),temp.GetInnerData());
		return temp;
	}

	inline Zeus::RealPlane<double>			Fourier2Real(const Zeus::FourierPlane<double>& rp)
	{
		int MachineIndex(GetFourMachineIndex(rp));
		Zeus::RealPlane<double>   temp((FourierFarm_[MachineIndex])->GetRealBufferSz());
		(FourierFarm_[MachineIndex])->Fourier2Real(rp.GetInnerData(),temp.GetInnerData());
		return temp;
	}

	inline const ObsPriorMapType&			GetMapRef(int freq,Zeus::MapType MapT)
	{
		ObsPriorMapCollType::iterator mapPiv;
		if((mapPiv = GetMapAtFreq(freq,MapT)) == ObsPriorMapsColl_.end())
		{
			err(ERRCOD_PWS_FREQNOTFOUND,ERRMSG_PWS_FREQNOTFOUND,L"PatchProcessor::GetMapRef");
		}	
		return *mapPiv;
	}

	inline	const ApodiseFunctor<double>& GetApodisingDevice(int mode)
	{
		ApodFunct_.SetMode(mode);
		return ApodFunct_;
	}

	inline Zeus::PwSCoreRandGen	*	GetNoiseDevice(void)
	{
		return NoiseGen_;
	}

	inline bool DetectedSrsEmpty(void) const
	{
		return PeakColl_.empty();
	}

	inline OddsEval*	GetOddsEval(void)
	{
		return OddsEval_;
	}
#ifdef AMI
	int AmiPwSLikelihood(double* LikeValue, const AmiLikeParams& ParamsIN, AmiLikeParams& ParamsOut,bool Noduocimation);
#endif
	~PatchProcessor(void);
private:
//
	inline  void					err(int errCode,const wchar_t* msg,const wchar_t* FunctName) const
	{
		std::wstring errstring(msg);
		errstring += std::wstring(PATCHNUMBERSTR);
		errstring += Zeus::PutNumber2Txt(PatchNumber_);
		throw Zeus::libException(errCode,errstring,FunctName);
	}
//
	inline  void					errPriorValues(const wchar_t* FunctName,int PatchNumber,double PriorValue) const
	{
		std::wstring errstring(std::wstring(SRCINDEXSTR) + Zeus::PutNumber2Txt(PatchNumber) + std::wstring(L"\n"));
		errstring	+= (std::wstring(ERRMSG_PWS_INVALPRIORVAL) + Zeus::PutNumber2Txt(PriorValue));
		throw Zeus::libException(ERRCOD_PWS_INVALPRIORVAL,errstring,FunctName);
	}
//
	inline  void					errPriorConfig(const wchar_t* FunctName,int PatchNumber,const wchar_t* msg) const
	{
		std::wstring errstring(std::wstring(SRCINDEXSTR) + Zeus::PutNumber2Txt(PatchNumber) + std::wstring(L"\n"));
		errstring	+= (std::wstring(msg));
		throw Zeus::libException(ERRCOD_PWS_INVALPRIORCONG,errstring,FunctName);
	}
//
	inline void 	AdjustMapsOrderIndex(void)
	{
		ObsPriorMapCollType::iterator			piv(ObsPriorMapsColl_.begin());
		ObsPriorMapCollType::const_iterator	const end(ObsPriorMapsColl_.end());
		for(int i=0;piv != end;++piv,++i)
		{
			piv->Order_ = i;
		}
	}
//
	typedef std::vector<Zeus::FourierMachine*>	FourierFarmType;

	double						GetOddsResultRawOptimAmpRange(int YPix,int XPix,int range, const Zeus::ObjFilterRealParams& inFilterParams,double& Isnr,int& sign);
	int							RemoveSrcGreyArea(void);
	void						RemoveSrcButCentral(void);
	double						MakeLikelihoodAtScale(int BinIndex,RealPlaneSurfType& tempLikelihooddouble,
													double&	LikeSurfISNR2AuxCte,double& LikeSurfAuxSigma,ObjFilterParams& objParam, int StoreNoise);
	double						MakeLikelihoodAtRealScale(const Zeus::ObjFilterRealParams& params,RealPlaneSurfType& tempLikelihood);
	void						PreprocessMaps(void);
	int							SetupPriors(Zeus::PriorsTypeCollType&	tPriType, int PriorsSel);
	void						MakeWhiteningStructs(double& LikeSurfISNR2Cte);
	void						EvalEquivBeamSigmaPix(int& FullWith,int& HalfMax);
	void						SetupPowell(void);
	void						InitializeFourierFarm(void);
	void						EvalPosMaxima(const RealPlaneSurfType& surf,int BinScale);
	void						SelectMaxima(void);
	void						RemoveAllButCentralPeak(void);
	double						GetMaximumLike(const RealPlaneSurfType& tempLikelihood,int YCent,int XCent,int Tol,int& YMax,int& XMax);
	void 						ScanLikeSurface(const ObjFilterParams& FilterParams,const RealPlaneSurfType& surf,double ISNR2,double Sigma);
	void						SetNoiseDevice(void);
	void						ReadObsMaps(int dataset);
	void						Read1ObsMap(int patchN, const ObsFreqsType& freq, Zeus::MapType mapT, Zeus::LArr2D<double>& ws) const;
	std::wstring				GetFileTypePrefix(Zeus::MapType mapT) const;
	std::wstring				GetFileTypeDir(Zeus::MapType mapT) const;
	void						MarkLowerFluxes(void);
	bool						Find1PeakOptimalParams(MaxLikeSurfMaximum& LikePeak,Zeus::PeakType &PeakResult);
	void						NonBlindMaxLikeReevaluation(MaxLikeSurfMaximum& LikePeak,Zeus::PeakType &PeakResult);
	void						PostProcessContours(Zeus::PeakType &PeakResult,Zeus::LArr2D<double>	&Buffer,double MaxAmpl,double MinAmpl,double MaxRadius,double MinRadius);
	double						MargPosGridSample(int YCentralPix,int XCentralPix,int PixRange,double thetaS,double FluxIntUnits);
	void						ContoursSetPriorsParameters(Zeus::PriorsTypeCollType& priorsTypes, std::vector<Zeus::Priors*> & PriorsColl, double MinRadius, double MaxRadius, double MinAmpl, double MaxAmpl,const Zeus::PeakType &PeakResult);
	void						WriteDegeneracyData2Media(const Zeus::LArr2D<double>& data,const Zeus::PeakType& PeakResult, double XMax, double XMin, double YMax, double YMin,Zeus::ContoursOutAux* aux);
	void						InjectNonBlindParams2Peak(Zeus::PeakType &PeakResult);
	RealPlaneSurfType			MakeWhitenedSource(double Amplitude,double DeltaY, double DeltaX,const Zeus::ObjFilterRealParams& ObjParams,double& CentralPixAmpl);
	void						MaskSource(double Amplitude,double YPos, double XPos,const Zeus::ObjFilterRealParams& ObjParams,double& CentralPixAmpl);
	void						NP_Eval_EvalParams(Zeus::PeakType & PeakRes,int stat);
	void						JF_EvalParams(Zeus::PeakType & PeakRes);
	void						GetMaskBorderValues(MaskBordersCollType& MskBColl,const RealPlaneSurfType& surf,int YPix,int XPix,int Swap);
	void						RenderSurfaceMask(const MaskBordersCollType& MskBColl,RealPlaneSurfType& tSurf,int Metric,int Swap);
	double						Eval_lnModelRatio(Zeus::PeakType & PeakRes);
	double						Eval_ExpectedBackgPeaks(void);
	int							EvalExpectNCorrectCtes(double& lnCv,double& lnCu);
	int							ProfVarGetBetaStep(double snr);
	int							ProfVarGetAlphaStep(double snr);
//-----------------------------------
//
	inline std::wstring			MakeDegeneracyPlotsFileName(const Zeus::PeakType& PeakResult)
	{
		wchar_t			buffer[BUFFERMAXCHAR];

		swprintf(buffer,BUFFERMAXCHAR,L"SZcolat%dlong%dpatchNo%d",Zeus::toInt(PeakResult.PeakNonBlindInfo_.SourcePtg_.theta * RAD2ARCMIN),
			Zeus::toInt(PeakResult.PeakNonBlindInfo_.SourcePtg_.phi * RAD2ARCMIN),PeakResult.PeakNonBlindInfo_.SrcIndex_);

		return std::wstring(buffer);
	}
//
	inline void						ReBuildMapsMainSurfFour(void)
	{
		MapsMainSurfFour_ = Real2Fourier(LikelihoodSurf_);
	}
	
	inline int	NP_Assessment(Zeus::PeakType & PeakRes, int checkRadii=0)
	{
		PeakRes.UsedSigmaThreshold_ = GlbVars_.NP_Sigma_;
		PeakRes.DetectionSigma_		= PeakRes.SrcAmplNormalised_ / PeakRes.GaussianIndex_;

		if(PeakRes.SrcAmplNormalised_ > GlbVars_.NP_Sigma_)
		{
			if(
				(checkRadii) && (GlbVars_.SZ_)
				&& (std::abs((GlbVars_.SrcMaxScale_ - PeakRes.RealParams_.RealScale_)/ GlbVars_.SrcMaxScale_) < 0.02) 
				&& (PeakRes.SrcAmplNormalised_ < 10)
				)
				return (PeakRes.PK_NP_DetectStat_ = Zeus::PeakType::PK_DET_NODET);

			return (PeakRes.PK_NP_DetectStat_ = Zeus::PeakType::PK_DET_OK);
		}
		return (PeakRes.PK_NP_DetectStat_ = Zeus::PeakType::PK_DET_NODET);
	}

	int	BayesLaplaceAssess(Zeus::PeakType & PeakRes);
//
	inline	int	IsOutsideDetectArea(int YPix,int XPix) const
	{
		if(	(YPix < GlbVars_.PatchBorder_) || (YPix > (GlbVars_.PatchSz_ - GlbVars_.PatchBorder_)) ||
			(XPix < GlbVars_.PatchBorder_) || (XPix >= (GlbVars_.PatchSz_ - GlbVars_.PatchBorder_))
			)
		return true;

		return false;
	}
//
	inline	int	IsOutsideDetectArea(const Zeus::PeakType & PeakRes) const
	{return IsOutsideDetectArea(PeakRes.Pos_.YPix_,PeakRes.Pos_.XPix_);}
//
	inline bool WasSrcMasked(int YCoord,int XCoord)
	{
		MaskSrcCoordCoolType::const_iterator	piv(MaskedSrcCoords_.begin());
		MaskSrcCoordCoolType::const_iterator	const end(MaskedSrcCoords_.end());
		for(;piv != end; ++piv)
		{
			if((std::sqrt(static_cast<double>(YCoord - piv->YCoord_)*static_cast<double>(YCoord - piv->YCoord_) + static_cast<double>(XCoord - piv->XCoord_)*static_cast<double>(XCoord - piv->XCoord_)))<=1.414)
				return true;
		}
		return false;
	}

	inline void					SubtractSource(double flux,double YPos, double XPos,const Zeus::ObjFilterRealParams& ObjParams,double DetectionSigma,double& CentralPixAmpl,bool Subtract = true)
	{
		if(DetectionSigma < SIGMASUBTRACT)
		{
			if(GlbVars_.SZ_)
			{
				CentralPixAmpl = -1.0;
				return;
			}
			SrcInfoType	ObjSurf(Zone_->GetSrcObjRaw(ObjParams));
			CentralPixAmpl = ObjSurf.CentralPixAmpl_ * Zone_->GetAntenna0pix(0);
			return;
		}
		else
		{
			if(
				(GlbVars_.SZ_) &&
				((std::abs((GlbVars_.SrcMaxScale_ - ObjParams.RealScale_)/ GlbVars_.SrcMaxScale_) < 0.02) && (DetectionSigma < 20.0))
			)
			{
				CentralPixAmpl = -1.0;
				return;
			}
		}

		int		YPix(Zeus::toInt(YPos + 0.5)),XPix(Zeus::toInt(XPos + 0.5));
		double	DeltaY,DeltaX;

		if(GlbVars_.NotAlignedObjs_)
		{DeltaY	= (YPos - static_cast<double>(YPix));DeltaX  = (XPos - static_cast<double>(XPix));}
		else
		{DeltaY	= 0.0;DeltaX  = 0.0;}

		RealPlaneSurfType	Dummy1(MakeWhitenedSource(flux,DeltaY, DeltaX,ObjParams,CentralPixAmpl));
		SubAddFunctor		Dummy(Subtract);

		//static  int DebugCounter(1);
		//char	Debugbuffer[1024];

		//sprintf(Debugbuffer,"LikelihoodSubtractionBefore_%d_",DebugCounter);
		//Zeus::DumpInOut_2d(Debugbuffer,256,256,256,0,LikelihoodSurf_.GetInnerData().begin(),1.0);

		Zeus::Transform(Dummy,LikelihoodSurf_,YPix,XPix,Dummy1,false,Zeus::UB_NOUSE);

		//sprintf(Debugbuffer,"LikelihoodSubtractionAfter_%d_",DebugCounter);
		//Zeus::DumpInOut_2d(Debugbuffer,256,256,256,0,LikelihoodSurf_.GetInnerData().begin(),1.0);

		if(!(GlbVars_.SZ_))
		{CentralPixAmpl *= Zone_->GetAntenna0pix(0);}
		else{CentralPixAmpl = -1.0;}
	}

	inline int					SubtractMaskSource(Zeus::PeakType &PeakResult)
	{
		if(
			!(GlbVars_.SZ_)	&& (PeakResult.SrcAmplNormalised_ > THRESHOLDFORMASKING) && 
			(!CurrentIter_ || WasSrcMasked(PeakResult.Pos_.YPix_,PeakResult.Pos_.XPix_))
			)
		{   
			// Make a backup of the first map; PS detection only
			// Very serious bug; 2nd iterations makes backup map without the bright sources !!!!!!!
			//
			if(CurrentIter_ != (GlbVars_.TwoStepsDetection_ - 1))
			{
				MaskedSrcCoords_.push_back(MaskSrcCoords(PeakResult.Pos_.YPix_,PeakResult.Pos_.XPix_));
			}
			double DummyMaskSource;
			MaskSource(PeakResult.SrcFlux_,PeakResult.Pos_.YCoord_, PeakResult.Pos_.XCoord_,PeakResult.UnTransParams_,DummyMaskSource);
			return 0; //Mask
		}
		
		double DummySubtractSource;
		SubtractSource(PeakResult.SrcFlux_,PeakResult.Pos_.YCoord_,PeakResult.Pos_.XCoord_,PeakResult.UnTransParams_,PeakResult.SrcAmplNormalised_,DummySubtractSource,true);
/*
		//DEBUG

		MapsMainSurfFourType	tempLike0Scale(Real2Fourier(LikelihoodSurf_));
		SrcInfoType ObjSurf(Zone_->GetObjOfRadiusBin(511));
		MapsMainSurfFourType tempFour(Real2Fourier(ObjSurf.surf_));
		MultiplyInPlace(tempLike0Scale,tempFour,Zeus::UB_NOUSE);
*/
		//char	bufferS[BUFFERMAXCHAR];	
		//sprintf(bufferS,"%s_%d_","LikelihoodSingleScale",511);
		//DumpInOut_2d(bufferS,256,256,256,0,LikelihoodSurf_.GetInnerData().begin(),1.0/(0.2817*static_cast<double>(GlbVars_.PatchSz_ * GlbVars_.PatchSz_)));

		return 1; // Subtract
	}
//
	inline void					Put1FMode2Vector(MapsCollFourierType& FModeVectors,int mode)
	{
		Zeus::AlgVect<std::complex<double> > d(MapsN_,false);
		Zeus::AlgVect<std::complex<double> >::DataInnerType::iterator	orgVector(d.GetInnerData().begin());
		ObsPriorMapCollType::const_iterator piv(ObsPriorMapsColl_.begin());
		ObsPriorMapCollType::const_iterator const end(ObsPriorMapsColl_.end());
		for(;piv != end;++piv,++orgVector)
		{
			*orgVector = *(piv->CxOrg_ + mode);
		}
		(FModeVectors.GetInnerData().begin() + mode)->Swap(d);
	}

	inline void					MakeFourSurfaces(void)
	{

		ObsPriorMapCollType::iterator		piv(ObsPriorMapsColl_.begin());
		ObsPriorMapCollType::iterator const end(ObsPriorMapsColl_.end());
		for(;piv != end;++piv)
		{
			{
				Zeus::FourierPlane<std::complex<double> > Dummy(Real2Fourier(piv->ws_));
				piv->C_ws_.Swap(Dummy);
			}
			piv->CxOrg_ = (piv->C_ws_).GetInnerData().begin();
			*((piv->C_ws_).GetInnerData().begin()) = std::complex<double>(0.0,0.0);
		}
	}

	inline void					PutMap2FModeVectors(MapsCollFourierType& FModeVectors)
	{
		int			YSz,XSz;
		FModeVectors.GetSz(YSz,XSz);
		const int	EndOff(YSz * XSz);
		int			CurrOff(0);
		for(;CurrOff != EndOff;++CurrOff)
		{
			Put1FMode2Vector(FModeVectors,CurrOff);
		}
	}

	inline void					ReleaseFModeMaps(void)
	{
		ObsPriorMapCollType::iterator piv(ObsPriorMapsColl_.begin());
		ObsPriorMapCollType::iterator const end(ObsPriorMapsColl_.end());
		for(;piv != end;++piv)
		{piv->C_ws_.Release();}
	}

	inline void					ComputeStats(void)
	{
		wchar_t	buffer[BUFFERMAXCHAR];	
		ObsPriorMapCollType::iterator piv(ObsPriorMapsColl_.begin());
		ObsPriorMapCollType::iterator const end(ObsPriorMapsColl_.end());
		for(;piv != end;++piv)
		{
			(PlanckInfo::Instance())->ComputeStats1Map(piv->ws_,STATREJECTTHRESHOLD,piv->Bias_,piv->Rms_,piv->NpixIncluded_);
			PRINTINTOBUFFERFUNCT
				(buffer,BUFFERMAXCHAR,PATCHREPFORMAT,
				piv->Freq_,
				piv->Bias_,
				piv->Rms_
				);

		(Zeus::ConManager::Instance())->PrintStr2Console(buffer);

		}
	}

	inline void					MaskMaps(MaskingType mask)
	{
		ObsPriorMapCollType::iterator piv(ObsPriorMapsColl_.begin());
		ObsPriorMapCollType::iterator const end(ObsPriorMapsColl_.end());
		for(;piv != end;++piv)
		{
			(PlanckInfo::Instance())->TruncateOutliers(piv->ws_,mask);

			//char	bufferS[BUFFERMAXCHAR];	
			//sprintf(bufferS,"%s_%d_","MapMasked",piv->Freq_);
			//Zeus::DumpInOut_2d(bufferS,256,256,256,0,piv->ws_.GetInnerData().begin(),1.0);
		}
	}

	inline ObsPriorMapCollType::iterator GetMapAtFreq(int freq,Zeus::MapType MapT)
	{
		ObsPriorMapCollType::iterator		piv(ObsPriorMapsColl_.begin());
		ObsPriorMapCollType::const_iterator	const pivEnd(ObsPriorMapsColl_.end());
		for(;piv!=pivEnd;++piv)
		{
			if(piv->Freq_ == freq)
				return piv;
		}
		return ObsPriorMapsColl_.end();
	}

	inline bool					CheckLikeMaximum(const MaxLikeSurfMaximum& LMax) const
	{
		if( (LMax.Coords_.XPix_ < GlbVars_.PatchBorder_) ||
			(LMax.Coords_.YPix_ < GlbVars_.PatchBorder_) ||
			(LMax.Coords_.XPix_ >= (GlbVars_.PatchSz_ - GlbVars_.PatchBorder_)) ||
			(LMax.Coords_.YPix_ >= (GlbVars_.PatchSz_ - GlbVars_.PatchBorder_))
		) return false;

		return true;
	}

	inline void	MaximaByQuadInterpolation(double val[9]) const
	{
		double hor(Zeus::FindMaxQuadFromSamples(val[3],val[4],val[5]));
		double ver(Zeus::FindMaxQuadFromSamples(val[1],val[4],val[7]));
	/*
		double diagUp(FindMaxQuadFromSamples(val[0],val[4],val[8]));
		double diagDown(FindMaxQuadFromSamples(val[6],val[4],val[2]));

		hor += diagUp	; ver += diagUp;
		hor += diagDown	; ver -= diagDown;
		val[0] = hor/((T)3.0);val[1] = ver/((T)3.0);
	*/
		val[0] = hor;
		val[1] = ver;
	}

	inline void	Centroid(double val[9]) const
	{
		double x(-(val[0])),y(-(val[0])),w(val[0]);
		y -= val[1];			  w += val[1];
		y -= val[2]; x += val[2]; w += val[2];
		
		x -= val[3]; w += val[3];
					 w += val[4];
		x += val[5]; w += val[5];

		y += val[6]; x -= val[6]; w += val[6];
		y += val[7];			  w += val[7];
		y += val[8]; x += val[8]; w += val[8];
		if(w != 0.0)
		{
			val[0] = x / w;
			val[1] = y / w;
		}
		else
		{
			val[0] = val[1] = 0.0;		
		}
	}

	inline	double y0Funct(double x)
	{
		return GammaFunct_.betainc(Evaly0_param_/(std::pow(std::cos(x),GlbVars_.ProfParam_.MNFW_alpha_) + Evaly0_param_)) * std::cos(x);
	}

	inline	double y0Cylindircal(double delta)
	{
		bool ok;
		Evaly0_param_ = std::pow(delta,GlbVars_.ProfParam_.MNFW_alpha_);
		return y0EvalIntegrator_->Integrate(0.0,PIOVER2,ok);
	}
//
	inline	void ClearISNRCache(void)
	{
		PatchFilterPropsCache_.clear();
	}
//------------------------------------
	double							Evaly0_param_;
	double							Evaly0_denominator_;

	const Zeus::PatchGeomType&		PatchesGeoInfo_;
	const GlobalScalarVarsType&		GlbVars_;
	const PlanckStaticType			StaticInfoColl_;
	double							ExpectedBackgPeaks_;
	int								NBrightPsInPatch_;
	int								DetectedObj_;
	int								EqBeamSigmaPix_;
	int								CurrentIter_;
	int								PatchNumber_;
	int								MapsN_;
	int								NSnakesDim_;
	FourierFarmType					FourierFarm_;
	BackgroundProcessor				*Background_;
	ZoneHandleType					Zone_;
	ObsPriorMapCollType				ObsPriorMapsColl_;
	Zeus::PwSCoreRandGen			*NoiseGen_;
	MapsCollFourierType				AntennaFourVectors_;
	AntMainSurfType					AntennaSFE_;
	RealPlaneSurfType				AntennaSFE_real_;				
	MapsMainSurfFourType			MapsMainSurfFour_;
	MaxLikeSurfMaximaCollType		MaximaColl_;
	RealPlaneSurfType				LikelihoodSurf_;
	PatchFilterPropsCacheType		PatchFilterPropsCache_;
	OddsEval						*OddsEval_;
	BrentLineMinimizer				*BrentMin_;
	PowellSnake						*PwSnake_;
	NestSamplPwSObserParamNew		*NestSampler_;
	Zeus::PeakCollType				PeakColl_;
	ApodiseFunctor<double>			ApodFunct_;
	MaskSrcCoordCoolType			MaskedSrcCoords_;
	ScalesPropsCollType				ScalesProps_;
	Zeus::GammaFuncts				GammaFunct_;
	Zeus::SimpsonIntegrator< double,PatchProcessor> *y0EvalIntegrator_;

};

#endif //PATCHPROCESSORH

