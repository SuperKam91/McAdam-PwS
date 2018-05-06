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
#ifndef PWS_GLOBALSH
#define PWS_GLOBALSH

#include "vec3.h"
#include "pointing.h"
#include "trafos.h"
#include "ZEUS_WorkSpace.h"
#include "ZEUS_InOut.h"
#include "ZEUS_PhysicsMath.h"
#include "ZEUS_GlobalValuesStore.h"
#include "ZEUS_FourierMachine.h"
#include "ZEUS_Bins.h"
#include "ZEUS_Object.h" 
#include "PWS_Strings.h"
#include "PWS_Defaults.h"
//-------------------------------------

#define HEALPIX_ATOM_PREC					float

#undef MULTI_THREAD

#define MAPREADERTINY			1e-20

#ifdef MULTI_THREAD
typedef Zeus::MThGlobalVars										GlobalVars;
typedef Zeus::HandleStorageBase<Loki::ObjectLevelLockable>		HandleRefCounter
#else
typedef Zeus::SThGlobalVars										GlobalVars;
typedef Zeus::HandleStorageBase<Loki::SingleThreaded>			HandleRefCounter;
#endif

// 1.5 Gb
#define	SRCCACHEMAXSZ			(static_cast<long long>(1024*1024*1024))
#define COORDSEPOCH				2000.0
#define ODDSVAL_OUTOFBOUNDS		1e4
#define LMBRENTTOL				((double)0.0005)
#define LMBRENTMAXITER			100
#define PWTOL					((double)0.0005)
#define PWMAXITER				100
#define NPARAMS_2D				4 //X,Y,R,A
#define NPARAMS_1D				3 //X,Y,R
#define FAILSBEFORELEAVE		5
#define THRESHOLDFORMASKING		40.0
#define MASKVALUE				-1e20
#define RMS_TEMPLATENOISERATIO	20.0
#define MINIMUMNOISELEVEL		1.0e-6
#define MINIMUNSIGMALEVEL		0.8
#define MAXRADIUSCTE			3.0
#define SIGMAPEAKSELECT			3.0
#define SIGMASUBTRACT			4.3
#define SIGMASUBTRACTBACK		5.5
#define STATREJECTTHRESHOLD		4.0
#define STATREJECTMASKTHRESHOLD	8.0
//#define STATREJECTMASKTHRESHOLD	12.0
#define PSDEFPATCHSIZE			256
#define DO_DETECTION			1
#define NO_DETECTION			0

typedef double							MinArrayAtomType;
typedef std::vector<MinArrayAtomType>	MinArrayType;

class PatchProcessor;
class BackgroundProcessor;
class Zone;
class Antenna;
class BrentLineMinimizer;
class OddsEval;
class PowellSnake;
class CompressCatalogue;
class NestedSamplerBase;
class NestSamplPwSObserParam;

struct CommandLineArgs
{
	int							Argc_;
	int							ContextID_;
	int							ExecId_;
	int							FirstPatchNumber_;
	int							LastPatchNumber_;
	std::wstring				fname_;
	std::wstring				OutFName_;
	std::wstring				FinalCatFName_;
	std::vector<std::wstring>	CatalogueFNameColl_;
};

typedef  Zeus::ObjHandle<Zone>::Type		ZoneHandleType;

enum  ValidNSide {NSIDE_512=512,NSIDE_1024=1024,NSIDE_2048=2048};

struct PlanckStaticInfoAtomType
{
	int			Freq_;
	double		AntFWHM_;
	int			OpNSide_;
	double		SZ_Signal_;
	inline bool operator<(const PlanckStaticInfoAtomType& sec) const
	{return Freq_ < sec.Freq_;}
};

typedef std::vector<PlanckStaticInfoAtomType>	PlanckStaticType;

struct ObsFreqsType
{
	int		sign_;
	int		freq_;

	explicit ObsFreqsType(int freq)
		:freq_(std::abs(freq) % 1000), sign_(freq)
	{}
};

typedef std::vector<ObsFreqsType>		ObsFreqsCollType;

typedef std::vector<double>				FilterScalesCollType;

struct PriorMapsInfoType
{
	int			MapType_;
	int			Freq_;
	double		MaxL_;
	double		PixNoiseRMS_;
};

typedef std::vector<PriorMapsInfoType>		PriorMapsInfoCollType;

typedef Zeus::FourierPlane<Zeus::AlgTriMatrix<double> >	WhiteStructType;

typedef Zeus::FourierPlane<Zeus::AlgVect<std::complex<double> > >	MapsCollFourierType;

typedef Zeus::FourierPlane<double>					AntMainSurfType;
typedef Zeus::FourierPlane<std::complex<double> >	MapsMainSurfFourType;
typedef Zeus::RealPlane<double>						RealPlaneSurfType;

struct AntennaPropsType{
	bool			IsCircSymm_;
	double			PeriodArc_;
	int				NSamples_;
	int				BoundSampleY_;
	int				BoundSampleX_;
	double			Gain_;
	double			Freq_;
	double			IntSq_;
	double			Int_;
	double			FWHM_;
};
//
// holds the binned values of the shape parameters
// These vales are not translated into real values
struct	ObjFilterParams
{
	// could be other parameter;
	// but they must change the way the filtering takes place

	int			BinScale_;
	inline ObjFilterParams(void)	// Default for point like objects
		:BinScale_(0)
	{}
	inline explicit ObjFilterParams(int BinScale)
		:BinScale_(BinScale)
	{}

	inline bool operator<(const ObjFilterParams& rhs) const
	{return BinScale_ < rhs.BinScale_;}
};
//
struct PatchFilterPropsType
{
	double	Isnr2_;
	inline PatchFilterPropsType(void)
		:Isnr2_(0.0)
	{}

	inline explicit PatchFilterPropsType(double Isnr2)
		:Isnr2_(Isnr2)
	{}
};
//
typedef std::map<ObjFilterParams,PatchFilterPropsType>	PatchFilterPropsCacheType;
//
struct SrcInfoType
{
	long long						LastHit_;
	Zeus::ObjFilterRealParams		RealParams_;
	RealPlaneSurfType				surf_;
	RealPlaneSurfType::AtomType		CentralPixAmpl_;
};
//
typedef std::map<ObjFilterParams, SrcInfoType>	SrcObjctCacheType;
//
struct	OddsResultType
{
	double	Correlation_;
	double  ISNR2_;

	inline OddsResultType(void)
		:Correlation_(0.0),ISNR2_(0.0)
	{}

	inline OddsResultType(double Correlation,double ISNR2)
		:Correlation_(Correlation),ISNR2_(ISNR2)
	{}

};
//
enum	PP_FindOptParamsResType		{PP_FOP_END = 0,PP_FOP_BRGHTSRC = 1,PP_FOP_RESCAN = 2};
//
struct MaskSrcCoords
{
	int YCoord_;
	int XCoord_;

	inline MaskSrcCoords(int y,int x)
		:YCoord_(y),XCoord_(x)
	{}
};
//
typedef std::vector<MaskSrcCoords>	MaskSrcCoordCoolType;
//
struct	OddEvalBinParams
{
	int					YPix_;
	int					XPix_;
	ObjFilterParams		Params_;

	OddEvalBinParams(int YPix,int XPix,const ObjFilterParams& Params)
		:YPix_(YPix),XPix_(XPix),Params_(Params)
	{}

	inline bool operator<(const OddEvalBinParams& rhs) const
	{
		if(Params_.BinScale_ != rhs.Params_.BinScale_)
			return Params_.BinScale_ < rhs.Params_.BinScale_;
		if(YPix_ != rhs.YPix_)
			return YPix_ < rhs.YPix_;
		return XPix_ < rhs.XPix_;
	}
};
//
struct MaxLikeSurfMaximum
{
	int						Removed_;
	int						NonBlindTarget_;
	ObjFilterParams			FiltParams_;
	Zeus::ObjPositionParams	Coords_;
	double					LikeRatio_;
	double					NormAmpl_;
	double					Correlation_;
	double					flux_;
	double					GaussianIndex_;
	// = SigmPwS_p / SigmPwS_t = (Sigma_p / Sigma_t^2) / (1/Sigma_t) =  Sigma_p / Sigma_t

	MaxLikeSurfMaximum(void)
		: Removed_(0),NonBlindTarget_(0)
	{}
	inline bool	operator<(const MaxLikeSurfMaximum& rhs) const
	{return LikeRatio_ > rhs.LikeRatio_;}
};

typedef	std::vector<MaxLikeSurfMaximum>	MaxLikeSurfMaximaCollType;

struct ObsPriorMapType
{
	typedef Zeus::RealPlane<double>						RealSurfType;
	typedef Zeus::FourierPlane<std::complex<double> >	ComplexSurfType;

	int				Order_;
	int				Freq_;
	double			Rms_;
	double			Bias_;
	int				NpixIncluded_;
	RealSurfType	ws_;
	ComplexSurfType	C_ws_;
	ComplexSurfType::DataInnerType::const_iterator	CxOrg_;	

	bool operator<(const ObsPriorMapType& rhs) const
	{
		return Order_ < rhs.Order_;
	}
};

typedef std::vector<ObsPriorMapType>	ObsPriorMapCollType;

struct	ScalesPropsType
{
	double	SrcScale_;
	double	LikeCentralPix_;
	double	Th1OverSigma_;
	double	LikelihoodSigma_;
	double	Gaussianity_;
	// = SigmPwS_p / SigmPwS_t = (Sigma_p / Sigma_t^2) / (1/Sigma_t) =  Sigma_p / Sigma_t

//SZCATALOGUEDEFAULTVALUE

	ScalesPropsType(void)
		:SrcScale_(SZCATALOGUEDEFAULTVALUE),LikeCentralPix_(SZCATALOGUEDEFAULTVALUE),
		Th1OverSigma_(SZCATALOGUEDEFAULTVALUE),
		LikelihoodSigma_(SZCATALOGUEDEFAULTVALUE),Gaussianity_(SZCATALOGUEDEFAULTVALUE)
	{}

	bool	operator<(const ScalesPropsType& rhs) const
	{
		return SrcScale_ < rhs.SrcScale_;
	}
};


typedef	std::vector<ScalesPropsType>	ScalesPropsCollType;


struct MaskMapType
{
	typedef Zeus::RealPlane<double> SurfType;
	int				Freq_;
	SurfType		ws_;

	MaskMapType(int Freq)
		:Freq_(Freq)
	{}
	MaskMapType(int Freq,Zeus::RealPlane<double> ws)
		:Freq_(Freq),ws_(ws)
	{}
};

typedef std::vector<MaskMapType>	MaskMapCollType;

struct BackgFourMap
{
	typedef Zeus::FourierPlane<std::complex<double> > SurfType;
	int				FreqOrd_;
	int				Freq_;
	SurfType::DataInnerType::const_iterator		Org_;	
	SurfType									ws_;

	inline bool operator<(const BackgFourMap& rhs) const
	{
		return FreqOrd_ < rhs.FreqOrd_;
	}
};

typedef std::vector<BackgFourMap>	BackgFourMapCollType;

typedef enum{MASK_CORE=0x02,MASK_BORDER=0x01,MASK_ALL=0x03} MaskingType;

struct GlobalScalarVarsType
{
	std::wstring		DirIn_;
	std::wstring		DirBuffer_;
	std::wstring		DirOut_;
	std::wstring		DirInMasks_;
	std::wstring		DirInMaps_;
	std::wstring		DirPointings_;
	std::wstring		DirInMasksIn_;
	std::wstring		MaskFileName_;
	std::wstring		OutputDataSetName_;
	std::wstring		IntermediateCatalogueFileName_;
	std::wstring		FinalCatalogueFileName_;
	int					MPI_rank_;
	int					Data2Buffer_;
	int					TwoStepsDetection_;
	int					NonBlindDetection_;
	int					TotNPatches_;
	int					N_ObsPlanes_;
	int					F_Patch_;
	int					L_Patch_;
	int					Sync_ID_;
	int					N_PriorPlanes_;
	int					PriorFreq_;
	int					DefaultUseBounds_;
	int					PatchSz_;
	int					PatchBorder_;
	int					NotAlignedObjs_;
	int					SearchPos_;
	int					Use2DFormula_;
	int					N_ScaleBins_;
	int					NSize_;
	int					ContextID_;
	int					AssessmentKind_;
	int					SubSources_;
	double				PixSz_;
	double				BoundTol_;
	double				NP_Sigma_;
	double				SSubGaussLevel_;
	int					Jf_Estimator_;
	int					SZ_;
	Zeus::SZPS_ProfParamType		ProfParam_;
	Zeus::SZ_ProfParamRangeType		ProfParamVar_;
	double				CacheSz_;
	double				PriorFluxExp_;
	double				PriorSrcMinScale_;
	double				PriorSrcMaxScale_;
	double				PriorSrcScaleExp_;
	double				PriorAvObjectsPatch_;
	double				SrcMaxScale_;
	double				PriorFluxMin_;
	double				PriorFluxMax_;
	double				PriorMaxZ_;
	double				PriorRadDistMax_;
	double				PriorRadDistMin_;
	double				PriorMassMax_;
	double				PriorMassMin_;
	double				PriorTempMax_;
	double				PriorTempMin_;
	double				PriorGassMassRatio_;
	int					PriorUsePressSchether_;
	double				Sigma8_;
	int					MuNe_NLivePoints_;
	int					MuNe_NIndivSamples_;
	double				MuNe_FractTolEv_;
	double				MuNe_xtraEnlgFactor_;
	int					OutputMergeAvGt_;
	int					OutputCoords_;
	int					OutputLat_;
	int					OutputDegrees_;
	double				OutputSigma_;
	double				OutputMeltMaxDist_;
	double				OutputMeltSZMaxLimit_;
	double				OutputFluxThreshold_;
	double				OutputGalCut_;
	double				OutputGalSigma_;
	int					OutputUnits_;
	double				OutputPurity_;
	// 0 - flux mJys / arcmin^2 (no translation); 1 - brightness MJys/Sr; 2 - Antenna K;3 - Thermo K

	int						ScalesFillerType_;
	int						ScalesMaxNumber_;
	double					ScalesFirstElement_;
	FilterScalesCollType	ScalesColl_;
	ObsFreqsCollType		FreqsColl_;

	inline int	GetObsFreqIndex(int freq) const
	{
		ObsFreqsCollType::const_iterator const Org(FreqsColl_.begin());
		ObsFreqsCollType::const_iterator piv(Org);
		ObsFreqsCollType::const_iterator const end(FreqsColl_.end());
		for(;(piv != end) && (piv->freq_ != freq);++piv)
			;
		if(piv == end)
			return -1;
		return static_cast<int>(piv - Org);
	}
};

#if	defined(WIN32) && !defined(HFIDMC) && !defined(LFIDPC) && !defined(AMI)
void PreProcessArgs(int argc, _TCHAR* argv[],CommandLineArgs& args);
#else
void PreProcessArgs(int argc, char* argv[],CommandLineArgs& args);
#endif
//
inline int	GetFreqOrder(int Freq)
{
	static int def_v(8);
	switch(Freq)
	{
	case 100:
		return 0;
	case 143:
		return 1;
	case 217:
		return 2;
	case 353:
		return 3;
	case 545:
		return 4;
	case 857:
		return 5;
	case 30:
		return 6;
	case 44:
		return 7;
	case 70:
		return 8;
	default:
		++def_v;
		return def_v;
	}
}
//
void	ProcessEnvVariables(CommandLineArgs &args);
void	PrintUsage(void);
void	PringArgs(const CommandLineArgs&	args);

#ifdef AMI

struct AmiLikeParams
{
	int		NParams_;
	double XCoord_;
	double YCoord_;
	double ThetaS_;
	double Ytot_;
	double alpha_;
	double beta_;
	double c500_;
	double gamma_;

	AmiLikeParams(void) :
		NParams_(4),
		XCoord_(SZCATALOGUEDEFAULTVALUE), YCoord_(SZCATALOGUEDEFAULTVALUE), ThetaS_(SZCATALOGUEDEFAULTVALUE), Ytot_(SZCATALOGUEDEFAULTVALUE),
		alpha_(SZCATALOGUEDEFAULTVALUE), beta_(SZCATALOGUEDEFAULTVALUE), c500_(SZCATALOGUEDEFAULTVALUE), gamma_(SZCATALOGUEDEFAULTVALUE)
	{}
	/*
	AmiLikeParams(double XCoord, double  YCoord, double  ThetaS, double  Ytot) :
		XCoord_(XCoord), YCoord_(YCoord), ThetaS_(ThetaS), Ytot_(Ytot),
		alpha_(SZCATALOGUEDEFAULTVALUE), beta_(SZCATALOGUEDEFAULTVALUE), c500_(SZCATALOGUEDEFAULTVALUE), gamma_(SZCATALOGUEDEFAULTVALUE)
	{}
	AmiLikeParams(double XCoord, double  YCoord, double  ThetaS, double  Ytot, double  alpha, double  beta) :
		XCoord_(XCoord), YCoord_(YCoord), ThetaS_(ThetaS), Ytot_(Ytot), alpha_(alpha), beta_(beta),
		c500_(SZCATALOGUEDEFAULTVALUE), gamma_(SZCATALOGUEDEFAULTVALUE)
	{}
	*/
};
//

extern "C" {
	int			AMI_initializePwSlike(const char* ParFName);
	int			AMI_fetchPatch(const int *ptrPatchNumber);
	void		AMI_ReleasePatch(void);
	int			AMI_PwSLikelihood(double* LikeValue,const int *ptrNParams, const double *In_params, double *Out_params,const int *ptrNoDuocimation);
}
#endif
//


#endif //PWS_GLOBALSH

