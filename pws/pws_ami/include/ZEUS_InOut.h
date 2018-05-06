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


#ifndef ZEUS_INOUTMANIP
#define ZEUS_INOUTMANIP

#ifdef WIN32
#include <io.h>
#include <errno.h>
#include <direct.h>
#include <stdio.h>
#include <stdlib.h>
#else
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#endif //WIN32

#include <memory>
#include <complex>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <locale>


#ifdef PWSMPI
#include <mpi.h>
#endif //PWSMPI

#ifdef HFIDMC
#ifdef WIN32
#include "PIOLib.h"
#else
#include "HL2_PIOLIB/PIOLib.h"
#endif // win32
#include "CamPwSCore64_param.h"
#endif //HFIDMC

#ifdef LFIDPC
#include "arr.h"
#include "focalplane_db.h"
#include "cxxutils.h"
#include "paramfile.h"
#include "iohandle_current.h" 
#endif // LFIDPC

#include "trafos.h"
#include "healpix_map.h"
#include "fitshandle.h"
#include "healpix_map_fitsio.h"
#include "ZEUS_General.h"
#include "ZEUS_SingletonSerialize.h"
#include "ZEUS_StorageBaseManip.h"
#include "ZEUS_Exceptions.h"
#include "ZEUS_WorkSpace.h"

#define GLOBAL_GNFWPROFILECOEFS_DEF					L"GNFWProfileCoefs"

#define WARNINGMSG				L"WARNING : "
#define FILESTRBUFFER			8192
#define SCANFMAXARGS			32
#define BUFFERMAXCHAR			4096
#define STRBUFFER				BUFFERMAXCHAR
#define DOUBLETXTMAXSZ			1024
#define SZCATALOGUEDEFAULTVALUE	-1.63750e+30
#define MC_MASKREJNAME_DEF		L"PwS_RejectionMask"
#define MC_BADPIXSMOOTH_DEF		L"PwS_TransparencyMask"
#define MC_MASKBADPIX_DEF		L"PwS_BadPixMask"
#define PATCH_GEOMPROPSFNAME	L"PatchGeomProps"
#define MC_PIPELINESYNC_DEF		L"PwS_SyncPixCut"

#define ICAT_FILEEXT			L".csv"
#define FCAT_FILEEXT			L".csv"
#define BIN_FILEEXT				L".dat"
#define PARAM_FILEEXT			L".par"
#define FITS_FILEEXT			L".fits"
//
#define ORANGE10NCOLUMNS		21
#define QARESULTSSZ				27		
//
enum	EnvIDsType{ENVID_DEFAULT=0,ENVID_HFI=1,ENVID_LFI=2,ENVID_FITS=3};
//
#ifdef WIN32
#define PRINTINTOBUFFERFUNCT	_snwprintf_s
#else
#define PRINTINTOBUFFERFUNCT	swprintf
#endif
//
#ifdef WIN32
#define PRINTINTBUFFERCHAR	_snprintf_s
#else
#define PRINTINTBUFFERCHAR	snprintf
#endif
//
namespace Zeus
{
//
template<typename T>
struct	JF_EstimatorType
{
	T		Mode_;
	T		Mean_;

	JF_EstimatorType(void)
		:Mode_(static_cast<T>(-1.0)),Mean_(static_cast<T>(-1.0))
	{}
};
//
struct ContoursOutAux
{
	double		In_Y_;
	double		In_Theta_;
	ContoursOutAux(void)
		:In_Y_(SZCATALOGUEDEFAULTVALUE),In_Theta_(SZCATALOGUEDEFAULTVALUE)
	{}
	ContoursOutAux(double Y,double Theta)
		:In_Y_(Y),In_Theta_(Theta)
	{}

};
//
struct	ObjPositionParams
{
	int							YPix_;
	int							XPix_;
	double						YCoord_;
	double						XCoord_;
	JF_EstimatorType<int>		JF_YPix_;
	JF_EstimatorType<int>		JF_XPix_;
	JF_EstimatorType<double>	JF_YCoord_;
	JF_EstimatorType<double>	JF_XCoord_;
	ObjPositionParams(void)
	{}
};
//
struct	ObjFilterRealParams
{
	// could be other parameter;
	// but they must change the way the filtering takes place
	double	RealScale_;
	double	a0_,a1_,a2_,a3_;
	int		i0_;
	inline ObjFilterRealParams(void)	// Default for point like objects
		:RealScale_(0), a0_(Zeus::pwsNaN), a1_(Zeus::pwsNaN), a2_(Zeus::pwsNaN), a3_(Zeus::pwsNaN), i0_(-1)
	{}
	inline ObjFilterRealParams(double RealScale, double a0, double a1, double a2, double a3, int i0)	// Default for point like object
		: RealScale_(RealScale), a0_(a0), a1_(a1), a2_(a2), a3_(a3), i0_(i0)
	{}
	inline explicit ObjFilterRealParams(double RealScale)
		:RealScale_(RealScale), a0_(Zeus::pwsNaN), a1_(Zeus::pwsNaN), a2_(Zeus::pwsNaN), a3_(Zeus::pwsNaN), i0_(-1)
	{}
};
//
struct SrcErrorBarsType
{
	// Assumed Gaussian
	// Error bars are 1 sigma
	double FluxErrorBar_;
	double RadiusErrorBar_;
	// Each component of the radius is gaussian
	// But the PDF of the radius itself is NOT gaussian
	double TotalPosErrorBar_;

	// Error bars are HPD ~95% (2 sigma)
	// General distribution
	double LowFluxErrorBar_;
	double LowRadiusErrorBar_;
	double HighFluxErrorBar_;
	double HighRadiusErrorBar_;
//
	double DegenSlope_;
	double DegenOrd_;
	double DegenCorr_;
	// The next fields do *NOT* hold the errors
	// but instead the Slope and Ord considering the errors on X and Y
	double DegenSlopeYErr_;
	double DegenOrdYErr_;
//
	SrcErrorBarsType(void)
		:FluxErrorBar_(SZCATALOGUEDEFAULTVALUE),RadiusErrorBar_(SZCATALOGUEDEFAULTVALUE),
		TotalPosErrorBar_(SZCATALOGUEDEFAULTVALUE),LowFluxErrorBar_(SZCATALOGUEDEFAULTVALUE),
		LowRadiusErrorBar_(SZCATALOGUEDEFAULTVALUE),HighFluxErrorBar_(SZCATALOGUEDEFAULTVALUE),
		HighRadiusErrorBar_(SZCATALOGUEDEFAULTVALUE),
		DegenSlope_(SZCATALOGUEDEFAULTVALUE),DegenOrd_(SZCATALOGUEDEFAULTVALUE),DegenCorr_(SZCATALOGUEDEFAULTVALUE),
		DegenSlopeYErr_(SZCATALOGUEDEFAULTVALUE),DegenOrdYErr_(SZCATALOGUEDEFAULTVALUE)
	{}
};
//
struct	ParamRangeType
{
	double	min_;
	double	max_;
	int		nbins_;

	inline ParamRangeType(void)
		:min_(SZCATALOGUEDEFAULTVALUE),max_(SZCATALOGUEDEFAULTVALUE),nbins_(-1)
	{}
	inline ParamRangeType(double min,double max,int nbins)
		:min_(min),max_(max),nbins_(nbins)
	{}

	inline int IsInit(void) const
	{return !((min_ == SZCATALOGUEDEFAULTVALUE) || (max_ == SZCATALOGUEDEFAULTVALUE) || (nbins_ < 1));}
};
//
struct	SZ_ProfParamRangeType
{
	ParamRangeType		alpha_;
	ParamRangeType		beta_;
	ParamRangeType		gamma_;
	ParamRangeType		C500_;

	SZ_ProfParamRangeType(void)
		:alpha_(), beta_(), gamma_(0.3081, 0.3081, 1), C500_(1.0, 1.0, 1)
	{}
	inline int IsInit(void) const
	{return (alpha_.IsInit() && beta_.IsInit()  && gamma_.IsInit() && C500_.IsInit());}
};
//
struct	SZPS_ProfParamType
{
	int			SZ_Profile_;
	int			ContBinsDef_;
	double		VirialRatio_;
	double		R500_Ratio_;
	double		MNFW_alpha_;
	double		MNFW_beta_;
	double		MNFW_gamma_;
	double		MNFW_C500_;
	double		MNFW_Ratio_CY500CYR500_;
	double		MNFW_2ndMoment1D_;
	double		FluxCalibCte_;
	double		RadiusCalCte_;

	SZPS_ProfParamType(void)
		:SZ_Profile_(2), ContBinsDef_(256), VirialRatio_(SZCATALOGUEDEFAULTVALUE), R500_Ratio_(SZCATALOGUEDEFAULTVALUE), MNFW_alpha_(SZCATALOGUEDEFAULTVALUE), MNFW_beta_(SZCATALOGUEDEFAULTVALUE),
		MNFW_gamma_(SZCATALOGUEDEFAULTVALUE),MNFW_C500_(SZCATALOGUEDEFAULTVALUE),MNFW_Ratio_CY500CYR500_(SZCATALOGUEDEFAULTVALUE),
		MNFW_2ndMoment1D_(SZCATALOGUEDEFAULTVALUE),
		FluxCalibCte_(1.0),RadiusCalCte_(1.0)
	{}
};
//
//
struct	PeakNonBlindInfoType
{
	pointing		X0Y0Ptg_; // centre of the patch
	pointing		SourcePtg_; // source coordinates
	pointing		QAIN_SrcPtg_; // QA true location of the source
	double			PredFlux_;
	double			PredRadius_;
	double			FluxBay_;
	double			RadiusBay_;
	double			SNR_;
	double			ErrRadius_;
	double			ErrFlux_;
	double			ErrPos_;
	double			ErrRadiusHigh_;
	double			ErrRadiusLow_;
	double			ErrFluxHigh_;
	double			ErrFluxLow_;
	double			QAIN_CyR500;
	double			QAIN_T500;
	double			QAIN_detectable;
	int				SrcIndex_;
	int				SrcXCoord_;
	int				SrcYCoord_;
	int				CollLstIndex_;

	PeakNonBlindInfoType(void)
		:X0Y0Ptg_(-1.0,-1.0),SourcePtg_(-1.0,-1.0),QAIN_SrcPtg_(-1.0,-1.0),SrcIndex_(-1),SrcXCoord_(-1),SrcYCoord_(-1),PredFlux_(SZCATALOGUEDEFAULTVALUE),PredRadius_(SZCATALOGUEDEFAULTVALUE),FluxBay_(SZCATALOGUEDEFAULTVALUE),RadiusBay_(SZCATALOGUEDEFAULTVALUE),
		SNR_(SZCATALOGUEDEFAULTVALUE),ErrRadius_(SZCATALOGUEDEFAULTVALUE),
		ErrFlux_(SZCATALOGUEDEFAULTVALUE),ErrPos_(SZCATALOGUEDEFAULTVALUE),ErrRadiusHigh_(SZCATALOGUEDEFAULTVALUE),ErrRadiusLow_(SZCATALOGUEDEFAULTVALUE),
		ErrFluxHigh_(SZCATALOGUEDEFAULTVALUE),ErrFluxLow_(SZCATALOGUEDEFAULTVALUE),QAIN_CyR500(SZCATALOGUEDEFAULTVALUE),QAIN_T500(SZCATALOGUEDEFAULTVALUE),QAIN_detectable(SZCATALOGUEDEFAULTVALUE),
		CollLstIndex_(-1)
	{}
};
//
struct ScaleLikeNoise
{
	double Scale_;
	double Like_;
	double Noise_;

	ScaleLikeNoise(void)
		:Scale_(SZCATALOGUEDEFAULTVALUE),Like_(SZCATALOGUEDEFAULTVALUE),
		Noise_(SZCATALOGUEDEFAULTVALUE)
	{}
};
//
typedef std::vector<ScaleLikeNoise>		ScaleLikeNoiseColl;
//
struct PeakType
{
	enum	PK_DetectStatType	{PK_DET_ERROR=0,PK_DET_OK=1,PK_DET_NODET=-1,PK_DET_BAYNOCONV=2};

	PeakNonBlindInfoType		PeakNonBlindInfo_;
	double						QAResult_[QARESULTSSZ];
	int							DetectID_;
	int							PK_NP_DetectStat_;
	int							PK_BayesDetectStat_;
	double						GaussianIndex_;
	// = SigmPwS_p / SigmPwS_t = (Sigma_p / Sigma_t^2) / (1/Sigma_t) =  Sigma_p / Sigma_t
	double						y0_;
	// DO NOT contain y0; Contains the Cte to convert Ytot into y0
	double						InitSrcFluxEst_;
	double						InitCorrelation_;
	double						CollListIndex_;
	double						SrcFlux_;
	double						SrcFlux_mJys_;
	double						SrcCompt_arcmin2_;
	double						SrcAmplNormalised_;
	double						DetectionSigma_;
	double						UsedSigmaThreshold_;
	double						JF_UsedSigmaThreshold_;
	double						JF_lnRhoTh_;
	double						JF_lnRho_;
	double						JF_lnEvidence_;
	double						JF_lnEvidenceErrBar_;
	// Evidence calibrate with the true priors on position and divided
	// by the MaxLikelihood
	double						JF_lnFormFactor_;
	double						JF_lnModelRatio_;
	JF_EstimatorType<double>	JF_Radius_;
	JF_EstimatorType<double>	JF_SrcFlux_;
	ObjPositionParams			Pos_;
	double						GalPt_Colat_;
	double						GalPt_Long_;
	double						CoordPt_Lat_;
	double						CoordPt_Long_;
	double						EquPt_Lat_;
	double						EquPt_Long_;
	ObjFilterRealParams			UnTransParams_;
	ObjFilterRealParams			RealParams_;
	SrcErrorBarsType			ErrorBars_;
	ScaleLikeNoiseColl			ScaleLikeNoise_;
	double						Odds_;
	double						ISNR2_;
	int							PatchNumber_;
	int							Masked_;
	int							UnderGalMask_;
	int							Hited_;
	PeakType(void)
		:PK_NP_DetectStat_(PK_DET_NODET),PK_BayesDetectStat_(PK_DET_ERROR),GaussianIndex_(1.0),JF_UsedSigmaThreshold_(-1.0),
		JF_lnRhoTh_(Zeus::logZERO),JF_lnRho_(Zeus::logZERO),JF_lnEvidence_(Zeus::logZERO),
		JF_lnEvidenceErrBar_(Zeus::logZERO),JF_lnFormFactor_(Zeus::logZERO),JF_lnModelRatio_(Zeus::logZERO),y0_(Zeus::logZERO),CollListIndex_(-1.0),
		UnderGalMask_(0),Hited_(0)
	{
		for(int i=0;i<QARESULTSSZ;++i)
		{QAResult_[i] = -1.0;}
	}
};
//
inline bool IsSamePeak(const PeakType& fst,const PeakType& sec)
{
	if(
		(fst.PatchNumber_ == sec.PatchNumber_)	&&
		(fst.Pos_.YPix_  == sec.Pos_.YPix_)		&&
		(fst.Pos_.XPix_  == sec.Pos_.XPix_)
	) return true;

	return false;
}
//
typedef	std::vector<PeakType>											PeakCollType;
//
typedef CollectionDiskStorageHelper<PeakCollType>						PeakCollReadbleType;
//
struct CatLineType
{
	int					ID_;
	int					Patch_;
	int					freqId_;
	int					CollLstIndex_;
	double				lnRho_;
	double				NormalAmpl_;
	double				DetectSigma_;
	double				GalLatDegs_;
	double				GalLongDegs_;
	double				PatchGalLatDegs_;
	double				PatchGalLongDegs_;
	double				PatchSpin_;
	double				FluxCompt_;
	double				FluxComptGLRT_;
	double				Radius_;
	double				RadiusGLRT_;
	SrcErrorBarsType	ErrorBars_;
	double				Gaussianity_;
	double				lnEvidence_;
	double				lnPenaltySrc_;
	// Now contains y0
	double				PatchMF_sigmaSqr_;
	double				SZ_ConversionCte_;
	// When not in chi2 mode contains the cte to convert Y5R500 -> y0
	ScaleLikeNoiseColl	ScaleLikeNoise_;
	double				QAResult_[QARESULTSSZ];


	CatLineType(void)
		:ID_(-1),Patch_(-1),freqId_(-1),CollLstIndex_(-1),lnRho_(SZCATALOGUEDEFAULTVALUE),lnEvidence_(SZCATALOGUEDEFAULTVALUE),lnPenaltySrc_(SZCATALOGUEDEFAULTVALUE),
		NormalAmpl_(SZCATALOGUEDEFAULTVALUE),DetectSigma_(SZCATALOGUEDEFAULTVALUE),GalLatDegs_(SZCATALOGUEDEFAULTVALUE),
		PatchGalLatDegs_(SZCATALOGUEDEFAULTVALUE),PatchGalLongDegs_(SZCATALOGUEDEFAULTVALUE),PatchSpin_(SZCATALOGUEDEFAULTVALUE),
		GalLongDegs_(SZCATALOGUEDEFAULTVALUE),FluxCompt_(SZCATALOGUEDEFAULTVALUE),FluxComptGLRT_(SZCATALOGUEDEFAULTVALUE),
		Radius_(SZCATALOGUEDEFAULTVALUE),RadiusGLRT_(SZCATALOGUEDEFAULTVALUE),Gaussianity_(SZCATALOGUEDEFAULTVALUE),
		PatchMF_sigmaSqr_(SZCATALOGUEDEFAULTVALUE),SZ_ConversionCte_(SZCATALOGUEDEFAULTVALUE),ErrorBars_()
	{
		for(int i=0;i<QARESULTSSZ;++i)
		{QAResult_[i] = -1.0;}	
	}
};
// ***********
struct QA_ProfilesLstLineType
{
	double	Alpha_;
	double	Beta_;
	double	Gamma_;
	double	C500_;
	double	ConvYtot2Y500_;

	QA_ProfilesLstLineType(void)
		:Alpha_(SZCATALOGUEDEFAULTVALUE),Beta_(SZCATALOGUEDEFAULTVALUE),Gamma_(SZCATALOGUEDEFAULTVALUE),
		C500_(SZCATALOGUEDEFAULTVALUE),ConvYtot2Y500_(SZCATALOGUEDEFAULTVALUE)
	{}
};
//
struct QA_CltLstLineType
{
	double	IN_theta;
	double	IN_phi;
	double	IN_CyR500;
	double	IN_T500;
	double	IN_M500;
	double	IN_z;
	double	IN_Vrec;
	double	IN_detectable;
	double	Cat_GLAT;
	double	Cat_GLON;
	double	Cat_yc;
	double	Cat_yc_error;
	double	Cat_CY5R500;
	double	Cat_LOWER_ERR_CY5R500F;
	double	Cat_UPPER_ERR_CY5R500F;
	double	Cat_T500;
	double	Cat_LOWER_ERR_T500F;
	double	Cat_UPPER_ERR_T500F;
	double	Cat_SNR;
	double	Cat_Patch_No;
	double	Cat_DEGEN_PARAM;
	double	Cat_SNR_30GHZ;
	double	Cat_SNR_44GHZ;
	double	Cat_SNR_70GHZ;
	double	Cat_SNR_100GHZ;
	double	Cat_SNR_143GHZ;
	double	Cat_SNR_217GHZ;
	double	Cat_SNR_353GHZ;
	double	Cat_SNR_545GHZ;
	double	Cat_SNR_857GHZ;
	double	Cat_FLUX5R500_30GHZ;
	double	Cat_FLUX5R500_44GHZ;
	double	Cat_FLUX5R500_70GHZ;
	double	Cat_FLUX5R500_100GHZ;
	double	Cat_FLUX5R500_143GHZ;
	double	Cat_FLUX5R500_217GHZ;
	double	Cat_FLUX5R500_353GHZ;
	double	Cat_FLUX5R500_545GHZ;
	double	Cat_FLUX5R500_857GHZ;
	double	Cat_CHI2;
	double	Cat_ERR_RAD;
	double	Cat_CY5R500F;
	double	Cat_T500F;
	int		index_;
//
	QA_CltLstLineType(void):
		IN_theta(SZCATALOGUEDEFAULTVALUE),IN_phi(SZCATALOGUEDEFAULTVALUE),IN_CyR500(SZCATALOGUEDEFAULTVALUE),IN_T500(SZCATALOGUEDEFAULTVALUE),IN_M500(SZCATALOGUEDEFAULTVALUE),
		IN_z(SZCATALOGUEDEFAULTVALUE),IN_Vrec(SZCATALOGUEDEFAULTVALUE),IN_detectable(SZCATALOGUEDEFAULTVALUE),
		Cat_GLAT(SZCATALOGUEDEFAULTVALUE),Cat_GLON(SZCATALOGUEDEFAULTVALUE),
		Cat_yc(SZCATALOGUEDEFAULTVALUE),Cat_yc_error(SZCATALOGUEDEFAULTVALUE),Cat_CY5R500(SZCATALOGUEDEFAULTVALUE),Cat_LOWER_ERR_CY5R500F(SZCATALOGUEDEFAULTVALUE),Cat_UPPER_ERR_CY5R500F(SZCATALOGUEDEFAULTVALUE),
		Cat_T500(SZCATALOGUEDEFAULTVALUE),Cat_LOWER_ERR_T500F(SZCATALOGUEDEFAULTVALUE),Cat_UPPER_ERR_T500F(SZCATALOGUEDEFAULTVALUE),Cat_SNR(SZCATALOGUEDEFAULTVALUE),Cat_Patch_No(SZCATALOGUEDEFAULTVALUE),
		Cat_DEGEN_PARAM(SZCATALOGUEDEFAULTVALUE),Cat_SNR_30GHZ(SZCATALOGUEDEFAULTVALUE),Cat_SNR_44GHZ(SZCATALOGUEDEFAULTVALUE),Cat_SNR_70GHZ(SZCATALOGUEDEFAULTVALUE),Cat_SNR_100GHZ(SZCATALOGUEDEFAULTVALUE),
		Cat_SNR_143GHZ(SZCATALOGUEDEFAULTVALUE),Cat_SNR_217GHZ(SZCATALOGUEDEFAULTVALUE),Cat_SNR_353GHZ(SZCATALOGUEDEFAULTVALUE),Cat_SNR_545GHZ(SZCATALOGUEDEFAULTVALUE),Cat_SNR_857GHZ(SZCATALOGUEDEFAULTVALUE),
		Cat_FLUX5R500_30GHZ(SZCATALOGUEDEFAULTVALUE),Cat_FLUX5R500_44GHZ(SZCATALOGUEDEFAULTVALUE),Cat_FLUX5R500_70GHZ(SZCATALOGUEDEFAULTVALUE),Cat_FLUX5R500_100GHZ(SZCATALOGUEDEFAULTVALUE),
		Cat_FLUX5R500_143GHZ(SZCATALOGUEDEFAULTVALUE),Cat_FLUX5R500_217GHZ(SZCATALOGUEDEFAULTVALUE),Cat_FLUX5R500_353GHZ(SZCATALOGUEDEFAULTVALUE),Cat_FLUX5R500_545GHZ(SZCATALOGUEDEFAULTVALUE),Cat_FLUX5R500_857GHZ(SZCATALOGUEDEFAULTVALUE),
		Cat_CHI2(SZCATALOGUEDEFAULTVALUE),Cat_ERR_RAD(SZCATALOGUEDEFAULTVALUE),Cat_CY5R500F(SZCATALOGUEDEFAULTVALUE),Cat_T500F(SZCATALOGUEDEFAULTVALUE),index_(-1)
	{}
};
//
typedef	std::vector<QA_CltLstLineType>			QA_CltLstCollType;
typedef	std::vector<QA_ProfilesLstLineType>		QA_ProfilesLstCollType;

// *********** 
typedef	std::vector<CatLineType>		CatLineCollType;
//
struct PCCFormatHeaderType
{
	static const int			CatNColumns_;
	static const int			CatHeaderNColumns_;
	static const wchar_t*		CatPrintStr_;
	static const wchar_t*		CatHeaderPrintStr_;
	static const wchar_t*		CatReadStr_;
	static const wchar_t*		CatHeaderReadStr_;

	int					CoordsType_;
	coordsys			CoordSystem_;
	double				Epoch_;
	int					NColumns_;
	int					NRows_;
	// new fields
	int					DetectionType_;
	int					Estimator_;
	int					PriorsType_;
	int					CollLstSz_;

	SZPS_ProfParamType	SZ_params_;
//
	PCCFormatHeaderType(void)
		:CoordsType_(1), // Needs to convert to CoLat and rads
		CoordSystem_(Galactic),Epoch_(2000.0),NColumns_(CatNColumns_),NRows_(-1), SZ_params_(),
		DetectionType_(-1),Estimator_(-1),PriorsType_(-1),CollLstSz_(-1)
	{}
};
//
struct OutputExtensionHeaderType
{
	std::wstring	DataSetName_;
//
	OutputExtensionHeaderType(void)
		:DataSetName_(L"Unknown")
	{}
};
//
struct OutputHeaderType
{
	PCCFormatHeaderType			PCCHeader_;
	OutputExtensionHeaderType	ExtHeader_;
//
	OutputHeaderType(void)
	{}
//
	OutputHeaderType(const PCCFormatHeaderType& PCCHeader)
		:PCCHeader_(PCCHeader)
	{}
//
	OutputHeaderType(const PCCFormatHeaderType& PCCHeader,const OutputExtensionHeaderType& ExtHeader)
		:PCCHeader_(PCCHeader),ExtHeader_(ExtHeader)
	{}
};
//
struct ExtFormatLnType
{
	int		NData_;
	double	CHI2_;
	double	Flux5R500_030_;
	double	SNR_030_;
	double	Flux5R500_044_;
	double	SNR_044_;
	double	Flux5R500_070_;
	double	SNR_070_;
	double	Flux5R500_100_;
	double	SNR_100_;
	double	Flux5R500_143_;
	double	SNR_143_;
	double	Flux5R500_217_;
	double	SNR_217_;
	double	Flux5R500_353_;
	double	SNR_353_;
	double	Flux5R500_545_;
	double	SNR_545_;
	double	Flux5R500_857_;
	double	SNR_857_;
	double	SNRF_;

	ExtFormatLnType(void)
		:
		NData_(-1),CHI2_(SZCATALOGUEDEFAULTVALUE),
		Flux5R500_030_(SZCATALOGUEDEFAULTVALUE),SNR_030_(SZCATALOGUEDEFAULTVALUE),
		Flux5R500_044_(SZCATALOGUEDEFAULTVALUE),SNR_044_(SZCATALOGUEDEFAULTVALUE),
		Flux5R500_070_(SZCATALOGUEDEFAULTVALUE),SNR_070_(SZCATALOGUEDEFAULTVALUE),
		Flux5R500_100_(SZCATALOGUEDEFAULTVALUE),SNR_100_(SZCATALOGUEDEFAULTVALUE),
		Flux5R500_143_(SZCATALOGUEDEFAULTVALUE),SNR_143_(SZCATALOGUEDEFAULTVALUE),
		Flux5R500_217_(SZCATALOGUEDEFAULTVALUE),SNR_217_(SZCATALOGUEDEFAULTVALUE),
		Flux5R500_353_(SZCATALOGUEDEFAULTVALUE),SNR_353_(SZCATALOGUEDEFAULTVALUE),
		Flux5R500_545_(SZCATALOGUEDEFAULTVALUE),SNR_545_(SZCATALOGUEDEFAULTVALUE),
		Flux5R500_857_(SZCATALOGUEDEFAULTVALUE),SNR_857_(SZCATALOGUEDEFAULTVALUE),
		SNRF_(SZCATALOGUEDEFAULTVALUE)
	{}
};
//
struct OutputLineType
{
	CatLineType			Cat_;
	ExtFormatLnType		Ext_;
};
//
typedef std::vector<OutputLineType>												OutputLineCollType;
//
typedef CollectionDiskStorageHelper<OutputLineCollType,OutputHeaderType>		OutputFormatType;
//
typedef CollectionDiskStorageHelper<CatLineCollType,PCCFormatHeaderType>		CatalogueFormatType;
//
typedef CollectionDiskStorageHelper<QA_CltLstCollType>							QA_CltLstCatalogueType;
// 
typedef CollectionDiskStorageHelper<QA_ProfilesLstCollType>						QA_ProfilesLstType;
//
struct	NonBlindHeader
{
	int					CoordsType_;
	coordsys			CoordSystem_;
	double				Epoch_;
	SZPS_ProfParamType	SZ_params_;
//
	NonBlindHeader(void)
		:CoordsType_(1),  // Needs to convert to CoLat and rads
		CoordSystem_(Galactic),Epoch_(2000.0)
	{}
};
//
struct PatchGeomHeaderType
{
	int			NSide_;
	int			Shuffle_;
	int			PtchSz_;
	int			PtchBorder_;
	int			NPatchesPerMainPix_;
	int			NTotalPatches_;
	int			CollListSz_;
	double		GalacticCut_;

	PatchGeomHeaderType(void):
		NSide_(-1),Shuffle_(0),PtchSz_(-1),PtchBorder_(-1),NPatchesPerMainPix_(-1),NTotalPatches_(-1),
		CollListSz_(-1),GalacticCut_(-1.0)
	{}
};
//
struct PatchGeomLineType
{
	int				PatchNumber_;
	int				BPixel_;
	double			Spin_;
	pointing		X0Y0Ptg_;	// centre of the patch
	pointing		SourcePtg_; // pointing of the source
	pointing		QAIN_SrcPtg_; // QA true location of the source
	pointing		X0Y0_;
	pointing		XLY0_;
	pointing		X0YL_;
	pointing		XLYL_;
	double			InitRejectRatio_;
	double			FinalRejectRatio_;
	int				PatchValid_;
	double			PredFlux_;
	double			PredRadius_;
	double			SNR_;
	double			ErrRadius_;
	double			ErrFlux_;
	double			ErrPos_;
	int				SrcIndex_;
	int				SrcXCoord_;
	int				SrcYCoord_;
	double			FluxBay_;
	double			RadiusBay_;
	double			ErrRadiusHigh_;
	double			ErrRadiusLow_;
	double			ErrFluxHigh_;
	double			ErrFluxLow_;
	double			QAIN_CyR500;
	double			QAIN_T500;
	double			QAIN_detectable;

	inline bool operator<(const PatchGeomLineType& rhs) const
	{return PatchNumber_ < rhs.PatchNumber_;}

	PatchGeomLineType(void):
		X0Y0Ptg_(-1.0,-1.0),SourcePtg_(-1.0,-1.0),QAIN_SrcPtg_(-1.0,-1.0),SrcXCoord_(-1),SrcYCoord_(-1),QAIN_CyR500(SZCATALOGUEDEFAULTVALUE),
		QAIN_T500(SZCATALOGUEDEFAULTVALUE),QAIN_detectable(1.0)
	{}
};
//
struct NonBlingCatLineType
{
	int					Index_;
	pointing			QAIN_SrcPtg_;
	pointing			ptg_;
	pointing			PatchPtg_;
	double				Spin_;
	double				PredFluxGLRT_;
	double				PredRadiusGLRT_;
	double				PredFluxBay_;
	double				PredRadiusBay_;
	double				SNR_;
	double				ErrFlux_;
	double				ErrRadius_;
	double				ErrPosition_;
	double				ErrRadiusHigh_;
	double				ErrRadiusLow_;
	double				ErrFluxHigh_;
	double				ErrFluxLow_;
	double				QAIN_CyR500;
	double				QAIN_T500;
	double				QAIN_detectable;
	int					QAIN_ProfileIndex_;
	int					flagged_;

	NonBlingCatLineType(void):
		QAIN_SrcPtg_(-2.0e10,-2.0e10),PatchPtg_(-2.0e10,-2.0e10),QAIN_CyR500(SZCATALOGUEDEFAULTVALUE),QAIN_T500(SZCATALOGUEDEFAULTVALUE),
		QAIN_detectable(SZCATALOGUEDEFAULTVALUE),QAIN_ProfileIndex_(-1),flagged_(-1)
	{}
};
//
typedef std::vector<NonBlingCatLineType> NonBlingCatCollType;
// 
typedef CollectionDiskStorageHelper<NonBlingCatCollType, NonBlindHeader>		NonBlingCatType;
//
typedef CollectionDiskStorageHelper<std::vector<PatchGeomLineType>,PatchGeomHeaderType>		PatchGeomType;
//
#ifdef HFIDMC
class	HFIDMC_parContent_helper: public HandleStorageBase<Loki::SingleThreaded>
{
	HFIDMC_parContent_helper(const HFIDMC_parContent_helper&);
	HFIDMC_parContent_helper& operator=(const HFIDMC_parContent_helper&);

	parContent*	hfi_params_;
public:
	HFIDMC_parContent_helper(void)
		:HandleStorageBase<Loki::SingleThreaded>(),hfi_params_(0)
	{}

	inline void	SetParContent(parContent* newPar) //Get the ownership 0f parContent
	{
		if(hfi_params_ && (hfi_params_ != newPar))
		{free_parContent(&hfi_params_);}
		hfi_params_ = newPar;
	}

	inline const parContent* GetParContent(void) const
	{return hfi_params_;}

	~HFIDMC_parContent_helper()
	{
		if(hfi_params_)
		{
		//	free_parContent(&hfi_params_);
		}
	}
};
//

typedef ObjHandle<HFIDMC_parContent_helper>::Type			HFIDMC_parContent_sptr;

//typedef std::auto_ptr<HFIDMC_parContent_helper>					HFIDMC_parContent_sptr;

typedef CollectionDiskStorageHelper<std::map<std::wstring, DBField>, HFIDMC_parContent_sptr > ParamVarsStoreType;
//
#else
typedef CollectionDiskStorageHelper<std::map<std::wstring, DBField> >						ParamVarsStoreType;
#endif
//
enum OBJECTTYPE
{
	OBJECTTYPE_MAP=0,
	OBJECTTYPE_TABLE=1,
	OBJECTTYPE_VECT=2
};
//
enum MapType
{
	MAPTYPE_OBSERVATION=0,
	MAPTYPE_BACKGROUND=1,
	MAPTYPE_INVALID=11
};
//
inline bool IsObservation(int MapType)
{return (MapType == MAPTYPE_OBSERVATION);}
//
inline bool IsBackground(int MapType)
{return (MapType == MAPTYPE_BACKGROUND);}
//
enum ReaderWriterStatType {UNINIT=0,INIT=1,RELEASED=2} ;
//
enum HealpixHeaderAtomInnerType{NO_KEY=-1,INT_TYPE=0,FLOAT_TYPE=2,DOUBLE_TYPE=3,BOOL_TYPE=4,STRING_TYPE=5}	;
//
struct HealpixHeaderAtomType
{

	std::string					Keyword_;
	HealpixHeaderAtomInnerType	ValueType_;
	union{
		int				intType_;
		float			floatType_;
		double			doubleType_;
		bool			boolType_;
	}PODValue_;
	std::string			stringValue_;
//
	HealpixHeaderAtomType(void)
		:ValueType_(NO_KEY)
	{}
	template<typename T>
	HealpixHeaderAtomType(const std::string& key,const T& value)
		:ValueType_(NO_KEY)
	{Set<T>(key,value);}
	HealpixHeaderAtomType(const std::string& key,int Type,void* value)
	{
		Keyword_ = key;
		switch(Type)
		{
		case INT_TYPE:
			PODValue_.intType_		= *(reinterpret_cast<int*>(value));
			ValueType_				= INT_TYPE;
			break;
		case FLOAT_TYPE:
			PODValue_.floatType_	= *(reinterpret_cast<float*>(value));
			ValueType_				= FLOAT_TYPE;
			break;
		case DOUBLE_TYPE:
			PODValue_.doubleType_	= *(reinterpret_cast<double*>(value));
			ValueType_				= DOUBLE_TYPE;
			break;
		case BOOL_TYPE:
			PODValue_.boolType_		= *(reinterpret_cast<bool*>(value));
			ValueType_				= BOOL_TYPE;
			break;
		case STRING_TYPE:
			stringValue_			= std::string(reinterpret_cast<char*>(value));
			ValueType_				= STRING_TYPE;
			break;
		default:
			ErrTypeMismatch();
		}
	}

	inline void Reset(void)
	{ValueType_ = NO_KEY;stringValue_.clear();}
	template<typename T>
	void	Set(const std::string& key,const T& value);
	template<typename T>
	void	Set(const T& value);
	inline	void	ErrTypeMismatch(void) const
	{
		throw Zeus::libException(ERROR_COD_TYPEMISMATCH,ERROR_MSG_TYPEMISMATCH,*this);
	}
};
//
namespace Private
{
template<typename T>
struct HHAtomTypeHelper
{};
template <> struct HHAtomTypeHelper<int>
{
	enum {VALUETYPE = INT_TYPE};
	static int& Storage(HealpixHeaderAtomType& t)
	{return t.PODValue_.intType_;}
};

template <> struct HHAtomTypeHelper<float>
{
	enum {VALUETYPE = FLOAT_TYPE};
	static float& Storage(HealpixHeaderAtomType& t)
	{return t.PODValue_.floatType_;}
};

template <> struct HHAtomTypeHelper<double>
{
	enum {VALUETYPE = DOUBLE_TYPE};
	static double& Storage(HealpixHeaderAtomType& t)
	{return t.PODValue_.doubleType_;}
};

template <> struct HHAtomTypeHelper<bool>
{
	enum {VALUETYPE = BOOL_TYPE};
	static bool& Storage(HealpixHeaderAtomType& t)
	{return t.PODValue_.boolType_;}
};

template <> struct HHAtomTypeHelper<std::string>
{
	enum {VALUETYPE = STRING_TYPE};
	static std::string& Storage(HealpixHeaderAtomType& t)
	{return t.stringValue_;}
};

}
//
struct HealpixHeaderType
{
	typedef std::vector<HealpixHeaderAtomType>	HealpixHeaderItemCollType;
	int							NSide_;
	HealpixHeaderItemCollType	Columns_;
};
//
inline std::wstring CharPtr2String(const wchar_t* str)
{return std::wstring(str);}
//
inline std::wstring CharPtr2String(const char* str)
{return Zeus::Achar2Wstr(str);}
//
class  ConManagerClass
{
public:
	inline void	PrintStr2Console(const std::wstring& str)
	{
		std::wcout << AddMsg2Log(str) << L"\n";
	}
//
	inline void	PrintStr2Console(const wchar_t* str)
	{
		std::wcout << AddMsg2Log(str) << L"\n";
	}
//
	inline void	PrintErr2Console(const wchar_t* str)
	{
		std::wcerr << AddErr2Log(str) << L"\n";
	}
//
	inline void	PrintErr2Console(const std::wstring& str)
	{
		std::wcerr << AddErr2Log(str) << L"\n";
	}
//
	inline void	PrintWarning2Console(const std::wstring& WarmMsg)
	{
		PrintStr2Console(std::wstring(WARNINGMSG) + WarmMsg + std::wstring(L"\n\n"));
	}
//
	inline void	PrintWarning2Console(const wchar_t* WarmMsg)
	{
		PrintStr2Console(std::wstring(WARNINGMSG) + std::wstring(WarmMsg) + std::wstring(L"\n\n"));
	}
//
	inline void	Flush(void)
	{
#ifdef	HFIDMC
#elif	LFIDPC
#else
		if(MsgFile_)
		{
			MsgFile_->flush();
		}
		if(ErrFile_)
		{
			ErrFile_->flush();
		}	
#endif
	}
//
	inline void	CloseStreams(int Sel = 0x03)
	{
#ifdef	HFIDMC
		if(Logger_)
		{
			PIOFreeLogger(&Logger_);
			Logger_ = 0;
		}
		if(Sel & 0x01) MsgFile_ = 0;
		if(Sel & 0x10) ErrFile_ = 0;
#else
		if(MsgFile_ && (Sel & 0x01))
		{
			MsgFile_->flush();
			delete MsgFile_;
			MsgFile_ = 0;
		}

		if(ErrFile_ && (Sel & 0x10))
		{
			ErrFile_->flush();
			delete ErrFile_;
			ErrFile_ = 0;
		}	
#endif
	}
//
#ifdef HFIDMC
	inline PIODBObjectLog		*GetLogger(void)
	{return Logger_;}
#endif
//
#if	defined(WIN32) && !defined(HFIDMC) && !defined(LFIDPC) && !defined(AMI)
	void	Initialise(int argc, _TCHAR* argv[],int ExecId);
#else
	void	Initialise(int argc, char* argv[],int ExecId);
#endif
//
	ConManagerClass(void)
		:MsgFile_(0),ErrFile_(0),execID_(-1)
#ifdef	HFIDMC
	,Logger_(0)
#endif
	{}
	~ConManagerClass(void)
	{
		CloseStreams();
	}
private:
	inline const std::wstring&	AddMsg2Log(const std::wstring& str)
	{
		if(MsgFile_ && !(str.empty()))
		{
			std::wstring	t(GetLogPrefix());
			std::wstring	t1(LeftTrim(str));

			
			if(!(t1.empty()))
			{
				t				+= t1;
#ifdef	HFIDMC
				t				+= std::wstring(L"\n");
				PIOInfoLogger(Logger_,const_cast<char *>(Wstr2Str(t).c_str()));
#else
				*MsgFile_		<< t;
				std::endl(*MsgFile_);
#endif
			}
		}
		return str;
	}
//
	inline const wchar_t*		AddMsg2Log(const wchar_t* str)
	{
		if(MsgFile_ && (wcslen(str)))
		{
			std::wstring	t(GetLogPrefix());
			std::wstring	t1(LeftTrim(std::wstring(str)));

			if(!(t1.empty()))
			{
				t				+= t1;
#ifdef	HFIDMC
				t				+= std::wstring(L"\n");
				PIOInfoLogger(Logger_,const_cast<char *>(Wstr2Str(t).c_str()));
#else
				*MsgFile_		<< t;
				std::endl(*MsgFile_);
#endif
			}
		}
		return str;
	}
//
	inline const std::wstring&	AddErr2Log(const std::wstring& Err)
	{
		if(ErrFile_ && !(Err.empty()))
		{
			std::wstring	t(GetLogPrefix());
			std::wstring	t1(LeftTrim(Err));
			if(!(t1.empty()))
			{
				t				+= t1;

#ifdef	HFIDMC
				t				+= std::wstring(L"\n");
				PIOErrorLogger(Logger_,const_cast<char *>(Wstr2Str(t).c_str()));
#else
				*ErrFile_		<< t;
				std::endl(*ErrFile_);
#endif
			}
		}
		return Err;
	}
//
	inline const wchar_t*	AddErr2Log(const wchar_t* Err)
	{
		if(ErrFile_ && (wcslen(Err)))
		{
			std::wstring	t(GetLogPrefix());
			std::wstring	t1(LeftTrim(std::wstring(Err)));

			if(!(t1.empty()))
			{
				t				+= t1;

#ifdef	HFIDMC
				t				+= std::wstring(L"\n");
				PIOErrorLogger(Logger_,const_cast<char *>(Wstr2Str(t).c_str()));
#else
				*ErrFile_		<< t;
				std::endl(*ErrFile_);
#endif
			}
		}
		return Err;
	}
//
	inline std::string			GetLogFileName(const char* str)
	{
		char	buf[BUFFERMAXCHAR];
		PRINTINTBUFFERCHAR(buf,BUFFERMAXCHAR,"%s_%04d.txt",str,(execID_>>2));
		return std::string(buf);
	}
//
	inline	std::wstring	GetTimeStamp(void)
	{
		time_t ltime;
		time( &ltime );
#ifdef WIN32
		char	buf[BUFFERMAXCHAR];
		if(ctime_s(buf,BUFFERMAXCHAR,&ltime))
		{
			return std::wstring(GLOBAL_TIMEERROR);
		}
#else
		char	*buf(ctime(&ltime));
		if(!buf)
		{
			return std::wstring(GLOBAL_TIMEERROR);
		}
#endif
		buf[strlen(buf) - 1] = ((char)0);
		return CharPtr2String(buf);
	}
//
#ifdef PWSMPI
	inline	std::wstring GetMPIprefix(void) 
	{
		int res;
		MPI_Initialized(&res);
		if(res)
		{
			int rank;
			wchar_t	buffer[BUFFERMAXCHAR];
			MPI_Comm_rank(MPI_COMM_WORLD,&rank);
			PRINTINTOBUFFERFUNCT
				(buffer,BUFFERMAXCHAR,L"-- MPI thread (%d) -- ",rank);
			return std::wstring(buffer);
		}
		else
		{return std::wstring(L"-- MPI not initislised -- ");}
	}
#endif //PWSMPI
//
	inline	std::wstring GetLogPrefix(void)
	{
		std::wstring	t(GetTimeStamp());
#ifdef PWSMPI
		t	+= GetMPIprefix();
#endif //PWSMPI
		t	+= std::wstring(GLOBAL_LOGSEPARATOR);
		return t; 
	}
//
#ifdef	HFIDMC
	PIODBObjectLog		*Logger_;
	PIOLOGGERSTRING		logStr;
	PIOErr				MyErr;
	int					ErrFile_;
	int					MsgFile_;
#else
	std::wofstream		*ErrFile_;
	std::wofstream		*MsgFile_;
#endif
	int					execID_;
//
};

//
template<typename T>
class GenHealpixReader
{
private:
	ReaderWriterStatType		Stat_;
	bool						xcptThrow_;
//
	GenHealpixReader(const GenHealpixReader& rhs);
	GenHealpixReader& operator= (const GenHealpixReader& rhs);
protected:
	typedef T					HealpixAtomType;
	Healpix_Map<T>*				HPixMap_;
//
	inline void		errX(int errCode,wchar_t* msg) const
	{
		std::wstring errstring(msg);

		errstring += std::wstring(L" -> ");
		errstring += do_GetCollID();
		throw Zeus::libException(errCode,errstring,*this);
	}
//
	virtual std::wstring	do_GetCollID(void) const = 0;
	virtual bool			do_Initialize(HealpixHeaderType * hd) = 0;
	virtual void			do_Release(void) = 0;
	virtual bool			do_Read(void) = 0;
public:
	//
	GenHealpixReader(bool xcptThrow=true)
		:Stat_(UNINIT),xcptThrow_(xcptThrow),HPixMap_(0)
	{}
	//
	inline void			GetCollStat(int& st,int& xpt)
	{st	= static_cast<int>(Stat_); xpt	= xcptThrow_;}
	//
	inline bool Initialize(HealpixHeaderType * hd)
	{
		delete HPixMap_;HPixMap_ = 0;
		HPixMap_ = new Healpix_Map<T>();
		bool  res(do_Initialize(hd));
		if(!res)
		{
			if(xcptThrow_) errX(ERROR_COD_ZEUSINITREADER,ERROR_MSG_ZEUSINITREADER);
			return false;
		}
		Stat_ = INIT;
		return true;
	}
	//
	inline std::wstring	GetCollID(void)
	{
		return do_GetCollID();
	}
	//
	void Release(Healpix_Map<T>& data)
	{
		if((Stat_ == UNINIT) || (Stat_ == RELEASED)) return;
		data.swap(*HPixMap_);
		delete HPixMap_; HPixMap_ = 0;
		do_Release();
		Stat_ = RELEASED;
	}
	//
	bool Read(void)
	{
		if((Stat_ == UNINIT) || (Stat_ == RELEASED))
		{
			if(xcptThrow_)
				errX(ERROR_COD_ZEUSREADERNOTINIT,ERROR_MSG_ZEUSREADERNOTINIT);
			return false;
		}
		bool  res(do_Read());
		if(!res)
		{
			delete HPixMap_; HPixMap_ = 0;
			if(xcptThrow_)
				errX(ERROR_COD_MYINOUT_FILEREAD,ERROR_MSG_MYINOUT_FILEREAD);
			return false;
		}
		return true;
	}
	//
	virtual ~GenHealpixReader(void)
	{delete HPixMap_;HPixMap_ = 0;}

};

//
template<typename T>
class GenHealpixWriter
{
private:
	ReaderWriterStatType		Stat_;
	bool						xcptThrow_;

	GenHealpixWriter(const GenHealpixWriter& rhs);
	GenHealpixWriter& operator= (const GenHealpixWriter& rhs);
protected:
	typedef T					HealpixAtomType;
	const std::wstring			HealpixMap_ID_;

	inline void		errX(int errCode,wchar_t* msg) const
	{
		std::wstring errstring(msg);

		errstring += std::wstring(L" -> ");
		errstring += do_GetCollID();
		throw Zeus::libException(errCode,errstring,*this);
	}

	virtual std::wstring		do_GetCollID(void) const = 0;
	virtual bool				do_Initialize(void) = 0;
	virtual int					do_Write(const Healpix_Map<T>& data,HealpixHeaderType * hd) = 0;
public:
	//
	GenHealpixWriter(bool xcptThrow=true)
		:Stat_(UNINIT),xcptThrow_(xcptThrow)
	{}
	//
	inline void			GetCollStat(int& st,int& xpt)
	{st	= static_cast<int>(Stat_); xpt	= xcptThrow_;}
	//
	inline std::wstring	GetCollID(void)
	{
		return do_GetCollID();
	}
	//
	inline bool Initialize(void)
	{
		bool  res(do_Initialize());
		if(!res)
		{
			if(xcptThrow_) errX(ERROR_COD_ZEUSINITWRITER,ERROR_MSG_ZEUSINITWRITER);
			return false;
		}
		Stat_ = INIT;
		return true;
	}
	//
	bool Write(const Healpix_Map<T>& data,HealpixHeaderType * hd)
	{
		if((Stat_ == UNINIT) || (Stat_ == RELEASED))
		{
			if(xcptThrow_) errX(ERROR_COD_ZEUSWRITERNOTINIT,ERROR_MSG_ZEUSWRITERNOTINIT);
			return false;
		}
		int		res(do_Write(data,hd));
		if(res < 0)
		{
			if(xcptThrow_) errX(ERROR_COD_MYINOUT_FILEWRITE,ERROR_MSG_MYINOUT_FILEWRITE);
			return false;
		}
		Stat_ = RELEASED;
		return true;
	}
	//
	virtual ~GenHealpixWriter(void)
	{}

};
//
template<class T>
class GenCollReader
{
private:
	ReaderWriterStatType	Stat_;
	bool					xcptThrow_;
	GenCollReader(const GenCollReader& rhs);
	GenCollReader& operator= (const GenCollReader& rhs);

protected:
	inline void		errX(int errCode,wchar_t* msg) const
	{
		std::wstring errstring(msg);

		errstring += std::wstring(L" -> ");
		errstring += do_GetCollID();
		throw Zeus::libException(errCode,errstring,*this);
	}

	virtual	T&						do_GetData(void) = 0;
	virtual void					do_DisposeData(void) = 0;
	virtual std::wstring			do_GetCollID(void) const = 0;
	virtual typename T::HeaderType *do_Initialize(void) = 0;
	virtual bool					do_Read(void) = 0;
public:
	typedef T		ReaderCollType;
	//
	GenCollReader(bool xcptThrow=true)
		:Stat_(UNINIT),xcptThrow_(xcptThrow)
	{}
	//
	inline void			GetCollStat(int& st,int& xpt)
	{st	= static_cast<int>(Stat_); xpt	= xcptThrow_;}
	//
	inline std::wstring	GetCollID(void)
	{
		return do_GetCollID();
	}
	//
	inline typename T::HeaderType *Initialize(void)
	{
		if(Stat_ == INIT)
		{do_DisposeData(); Stat_ = RELEASED;}

		typename T::HeaderType *res(do_Initialize());
		if(!res)
		{
			if(xcptThrow_)
				errX(ERROR_COD_ZEUSINITREADER,ERROR_MSG_ZEUSINITREADER);
		}
		Stat_ = INIT;
		return res;
	}
	//
	inline void Release(T & data)
	{
		if(Stat_ != INIT) return;
		data.swap(do_GetData());
		do_DisposeData();
		Stat_ = RELEASED;
	}
	//
	inline bool Read(void)
	{
		if(Stat_ != INIT)
		{
			if(xcptThrow_) errX(ERROR_COD_ZEUSREADERNOTINIT,ERROR_MSG_ZEUSREADERNOTINIT);
			return false;
		}

		bool  res(do_Read());
		if(!res)
		{
			if(xcptThrow_) errX(ERROR_COD_MYINOUT_FILEREAD,ERROR_MSG_MYINOUT_FILEREAD);
			return false;
		}
		return true;
	}
	//
	virtual ~GenCollReader(void)
	{}
};
//
template<class T>
class GenCollWriter
{
private:
	ReaderWriterStatType		Stat_;
	bool						xcptThrow_;
	GenCollWriter(const GenCollWriter& rhs);
	GenCollWriter& operator= (const GenCollWriter& rhs);

protected:
	typedef		T				DataType;

	inline void					errX(int errCode,wchar_t* msg) const
	{
		std::wstring errstring(msg);

		errstring += std::wstring(L" -> ");
		errstring += do_GetCollID();
		throw Zeus::libException(errCode,errstring,*this);
	}

	virtual std::wstring		do_GetCollID(void) const = 0;
	virtual void				do_DisposeData(void) = 0;
	virtual bool				do_Initialize(void) = 0;
	virtual int					do_Write(const T & data) = 0;
	virtual int					do_Flush(void) = 0;
	virtual int					do_Remove(void) = 0;
	//
public:
	GenCollWriter( bool xcptThrow=true)
		:Stat_(UNINIT),xcptThrow_(xcptThrow)
	{}
	//
	inline int			Remove(void)
	{
		if(Stat_ !=UNINIT)
			return -1;
		return do_Remove();
	}
	//
	inline bool			Initialize(void)
	{
		if(Stat_ == INIT)
		{do_DisposeData();Stat_ = RELEASED;}
		
		bool  res(do_Initialize());
		if(!res)
		{
			if(xcptThrow_) errX(ERROR_COD_ZEUSINITWRITER,ERROR_MSG_ZEUSINITWRITER);
			return false;
		}
		Stat_ = INIT;
		return true;
	}
	//
	inline std::wstring	GetCollID(void)
	{
		return do_GetCollID();
	}
	//
	inline void			GetCollStat(int& st,int& xpt)
	{st	= static_cast<int>(Stat_); xpt	= xcptThrow_;}
	//
	inline int			Write(const T & data)
	{
		if((Stat_ == UNINIT) || (Stat_ == RELEASED))
		{
			if(xcptThrow_)
				errX(ERROR_COD_ZEUSWRITERNOTINIT,ERROR_MSG_ZEUSWRITERNOTINIT);
			return false;
		}
		int   res(do_Write(data));
		if(res < 0)
		{
			if(xcptThrow_)
				errX(ERROR_COD_MYINOUT_FILEWRITE,ERROR_MSG_MYINOUT_FILEWRITE);
			return -1;
		}
		return res;
	}
	//
	inline int			Flush(void)
	{
		if((Stat_ == UNINIT) || (Stat_ == RELEASED))
		{return -1;}
		int res(do_Flush());
		if(res < 0)
		{
			if(xcptThrow_)
				errX(ERROR_COD_MYINOUT_FILEWRITE,ERROR_MSG_MYINOUT_FILEWRITE);
			return -1;
		}
		do_DisposeData();
		Stat_ = RELEASED;
		return 0;
	}
//
	virtual ~GenCollWriter(void)
	{}
};
//
template<typename T>
inline int		Number2Str(wchar_t* buffer,int sz,const wchar_t* format,T value,int width=4)
{
#ifdef WIN32
	return _snwprintf_s(buffer,sz,sz,format,width,value);
#else
	return swprintf(buffer,sz,format,width,value);
#endif
}

template<typename T>
inline int		Number2Str(wchar_t* buffer,int sz,const wchar_t* format,std::complex<T> value,int width=4)
{
#ifdef WIN32
	return _snwprintf_s(buffer,sz,sz,format,width,value.real(),width,value.imag());
#else
	return swprintf(buffer,sz,format,width,value.real(),width,value.imag());
#endif
}


//
namespace Private
{
template<typename T>
struct NumTxtFormatType
{};

template<>
struct NumTxtFormatType<int>
{
	static  const wchar_t* PrintfFormat;
};

template<>
struct NumTxtFormatType<double>
{
	static  const wchar_t* PrintfFormat;
};

template<>
struct NumTxtFormatType<float>
{
	static  const wchar_t* PrintfFormat;
};

template<>
struct NumTxtFormatType<std::complex<double> >
{
	static  const wchar_t* PrintfFormat;
};

template<>
struct NumTxtFormatType<std::complex<float> >
{
	static  const wchar_t* PrintfFormat;
};


#ifdef HFIDMC
template<typename T>
struct PIOTypeChooserType
{};

template<>
struct PIOTypeChooserType<int>
{static  PIOSTRING PIO_TYPE;};

template<>
struct PIOTypeChooserType<float>
{static  PIOSTRING PIO_TYPE;};

template<>
struct PIOTypeChooserType<double>
{static  PIOSTRING PIO_TYPE;};

template<>
struct PIOTypeChooserType<bool>
{static  PIOSTRING PIO_TYPE;};

template<>
struct PIOTypeChooserType<std::string>
{static  PIOSTRING PIO_TYPE;};

#endif
}
//
#ifdef HFIDMC
inline int GetPIOType(PIOSTRING str)
{
	PIOSTRING	buf[]={	"PIOINT",
						"$$$$$$",
						"PIOFLOAT",
						"PIODOUBLE",
						"PIOBOOL",
						"PIOSTRING"
					 };
	for(int i=0;i<5;++i)
	{
		if(!strcmp(str,buf[i])) return i;
	}
	return -1;
}

union HFIDMC_ReadKeyBufType
{
	PIOINT				intType_;
	PIOFLOAT			floatType_;
	PIODOUBLE			doubleType_;
	PIOINT				boolType_;
	PIOSTRING			strType_;
};

#endif
//
std::wstring& 	MakePath(std::wstring& str);
std::wstring	CorrectDir(const std::wstring& dir,int ForceAppend=0);
int				RemoveFile(OBJECTTYPE ObjT,const std::wstring& Dir,const std::wstring& file);
//
template<typename T>
inline std::wstring	PutNumber2Txt(const T& value,int width=4)
{
	wchar_t buffer[STRBUFFER];

	Number2Str<T>(buffer,DOUBLETXTMAXSZ,Private::NumTxtFormatType<T>::PrintfFormat,value,width);
	return buffer;
}
//
inline void	PrintErr2Console(const wchar_t* str)
{
#ifdef WIN32
	std::wcerr << str << L"\n";
#else
	std::cerr << Zeus::Wstr2Str(std::wstring(str)) << "\n";
#endif //WIN32
}
//
inline std::wstring	RemoveInstanceName(const std::wstring& Dir)
{
	if(Dir.empty()) return std::wstring(L"");
	
	std::wstring::size_type chPtr;

	if((chPtr = Dir.find(L"run0123456789_")) == std::wstring::npos)
		return Dir;

	return Dir.substr(0,chPtr);
}
//
inline std::wstring	ExtractFileName(const std::wstring& in,std::wstring& Dir)
{
	Dir.clear();
	if(in.empty()) return std::wstring(L"");

#ifdef WIN32
	std::wstring	DirSep(L":\\");
#else
	std::wstring	DirSep(L":/");
#endif //WIN32

	std::wstring tin(Zeus::FullTrim(in));
	std::wstring::size_type	SepInd(tin.find_last_of(DirSep));
	if(SepInd == std::wstring::npos)
		return tin;
	Dir = CorrectDir(tin.substr(0,SepInd));
	return tin.substr(SepInd+1);
}
//
inline std::wstring	SplitFileName(const std::wstring& in,std::wstring& ext)
{
	ext.clear();
	if(in.empty()) return std::wstring(L"");

	std::wstring	ExtSep(L".");

	std::wstring::size_type	SepInd(in.find_last_of(ExtSep));
	if(SepInd == std::wstring::npos)
		return in;
	ext = in.substr(SepInd+1);
	return in.substr(0,SepInd);
}
//
inline void	PrintErr2Console(const std::wstring& str)
{
#ifdef WIN32
	std::wcerr << str << L"\n";
#else
	std::cerr << Zeus::Wstr2Str(str) << "\n";
#endif //WIN32
}
//
inline bool CheckFileExists(const std::wstring& fname)
{
#ifdef	WIN32
	if(_waccess_s(fname.c_str(),0) == 0) return true;
#else
	if(access(Zeus::Wstr2Str(fname).c_str(),F_OK) == 0) return true;
#endif //WIN32
	return false;
}
//
inline bool CheckDir(const std::wstring& DirName)
{
#ifdef WIN32
	if(_waccess_s(DirName.c_str(),0) == 0)
		return true;
	return false;
#else
	if(access(Zeus::Wstr2Str(DirName).c_str(),F_OK) == 0)
		return true;
	return false;

#endif //WIN32
}
//
inline bool CreateDir(const std::wstring& DirName)
{
#ifdef WIN32
	if(_waccess_s(DirName.c_str(),0) == 0)
		return true;
	if((_wmkdir(DirName.c_str()) == 0) ||
		(errno == EEXIST)) return true;

	return false;
#else
	if(access(Zeus::Wstr2Str(DirName).c_str(),F_OK) == 0)
		return true;
	if((mkdir(Zeus::Wstr2Str(DirName).c_str(),
		S_IRUSR | S_IWUSR | S_IXUSR | S_IRGRP | S_IWGRP | S_IXGRP | S_IROTH | S_IWOTH | S_IXOTH) == 0) ||
		(errno == EEXIST)) return true;

	return false;

#endif //WIN32
}
//
void	CreatReadableHealpixFile(int EnvId,const std::wstring& DirName,const std::wstring& FName,
								 coordsys CoordSys,const Healpix_Map<float>& HealpixData);
//
inline int	GetCoordSysFromStr(const std::string& str)
{
	switch(ToCase(str.at(0),false)) //Upper case
	{
	case 'G':
		return Galactic;
	case 'C':
		return Equatorial;
	case 'E':
		return Ecliptic;
	default:
		return -1;
	}	
}
//
inline void	CreatReadableExternalHealpixFile(const std::wstring& DirName,const std::wstring& FName,
								 coordsys CoordSys,const Healpix_Map<float>& HealpixData)
{
#ifdef WIN32
	std::wstring	MapsExt(L"Masks\\");
#else
	std::wstring	MapsExt(L"Masks/");
#endif
	std::wstring	tDir(Zeus::RemoveInstanceName(DirName)+MapsExt);
	Zeus::CreateDir(tDir);
	Zeus::CreatReadableHealpixFile(1000,tDir,FName,CoordSys,HealpixData);
}
//
inline int	GetCoordSysFromStr(const std::wstring& str)
{
	return GetCoordSysFromStr(Wstr2Str(str));
}
//
inline	std::string				GetStrFromCoordsys(coordsys coord)
{
	switch(coord)
	{
	case Equatorial:
		return std::string("EQUATORIAL");
	case Ecliptic:
		return std::string("ECLIPTIC");
	default:
		return std::string("GALACTIC");
	}	
}
//
inline	std::wstring			GetWStrFromCoordsys(coordsys coord)
{
	return Achar2Wstr(GetStrFromCoordsys(coord).c_str());
}
//
inline	Healpix_Ordering_Scheme	GetOrderingFromStr(const std::string& str)
{
	return	((str == std::string("RING")) ? RING : NEST);
}
//
inline	Healpix_Ordering_Scheme	GetOrderingFromStr(const std::wstring& str)
{
	return GetOrderingFromStr(Wstr2Str(str));
}
//
inline	std::string				GetStrFromOrdering(Healpix_Ordering_Scheme ord)
{
	return std::string((ord==NEST)?"NEST":"RING");
}
//
inline	std::wstring			GetWStrFromOrdering(Healpix_Ordering_Scheme ord)
{
	return std::wstring((ord==NEST)?L"NEST":L"RING");
}
//
#ifdef MULTI_THREAD
typedef Multi_SingletonHolder<ConManagerClass>::Type		ConManager;
#else
typedef Single_SingletonHolder<ConManagerClass>::Type		ConManager;
#endif
//
template<typename T>
inline void	DumpScalarVariable(const wchar_t* name,const T& var)
{
	wchar_t buff[STRBUFFER];
#ifdef WIN32
	_snwprintf_s(buff,STRBUFFER-1,STRBUFFER-1,L"%-20.20s",name);
#else
	swprintf(buff,STRBUFFER-1,L"%-20.20S",name);
#endif
	std::wstring	t(buff);
	t	+=	L" => ";
	t	+=	Zeus::PutNumber2Txt(var);
	(Zeus::ConManager::Instance())->PrintStr2Console(t);
}
//
template<typename T>
inline void	DumpScalarVariable(const std::wstring& name,const T& var)
{
	wchar_t buff[STRBUFFER];
#ifdef WIN32
	_snwprintf_s(buff,STRBUFFER-1,STRBUFFER-1,L"%-20.20s",name.c_str());
#else
	swprintf(buff,STRBUFFER-1,L"%-20.20S",name.c_str());
#endif
	std::wstring	t(buff);
	t	+=	L" => ";
	t	+=	Zeus::PutNumber2Txt(var);
	(Zeus::ConManager::Instance())->PrintStr2Console(t);
}
//
inline void	DumpScalarVariable(const wchar_t* name,const std::wstring& var)
{
	wchar_t buff[STRBUFFER];
#ifdef WIN32
	_snwprintf_s(buff,STRBUFFER-1,STRBUFFER-1,L"%-20.20s",name);
#else
	swprintf(buff,STRBUFFER-1,L"%-20.20S",name);
#endif
	std::wstring	t(buff);
	t	+=	L" => ";
	t	+=	var;
	(Zeus::ConManager::Instance())->PrintStr2Console(t);
}
//
inline void	DumpScalarVariable(const std::wstring& name,const std::wstring& var)
{
	wchar_t buff[STRBUFFER];
#ifdef WIN32
	_snwprintf_s(buff,STRBUFFER-1,STRBUFFER-1,L"%-20.20s",name.c_str());
#else
	swprintf(buff,STRBUFFER-1,L"%-20.20S",name.c_str());
#endif
	std::wstring	t(buff);
	t	+=	L" => ";
	t	+=	var;
	(Zeus::ConManager::Instance())->PrintStr2Console(t);
}

//
#ifdef HFIDMC
inline	void	HFIDMC_ReportMemoryLeak(PIOErr MyErr,char* formatStr,const char* fieldName)
{
	char buff[STRBUFFER];
	sprintf(buff,formatStr,fieldName,PIOErrMess(MyErr));
	(Zeus::ConManager::Instance())->PrintStr2Console(Achar2Wstr(buff));
}

inline PIOErr	EraseGroup(const std::wstring& GName)
{
	std::string	tStr(Wstr2Str(GName));

	if(PIOCheckGroup(const_cast<char *>(tStr.c_str())))
	{
		return PIODeleteGroup(const_cast<char *>(tStr.c_str()));
	}
	return (PIOErr) 0;
}

#endif

#if	defined(WIN32) && !defined(HFIDMC) && !defined(LFIDPC) && !defined(AMI)
inline	std::wstring	GetEnvVar(wchar_t * varName)
{
	wchar_t		buff[STRBUFFER];
	size_t		pReturnValue;

	if(_wgetenv_s(&pReturnValue,buff,varName) || (pReturnValue == 0))
		return std::wstring();
	
	return std::wstring(buff);
}
#else
inline	std::wstring	GetEnvVar(char * varName)
{
	char *ptr;

	if(ptr = getenv(varName))
	{return Achar2Wstr(ptr);}
	
	return std::wstring();

}
#endif

template<typename T>
GenHealpixWriter<T>
*GetGenHealpixWriter(Loki::Type2Type<T>,int ID,const std::wstring& DirName,const std::wstring& MapName,
					 int NSide,Healpix_Ordering_Scheme Scheme,coordsys CoordSys);
//
template<typename T>
GenHealpixReader<T>
*GetGenHealpixReader(Loki::Type2Type<T>,int ID,const std::wstring& DirName,const std::wstring& MapName);
//
template<typename T>
GenCollReader<LArr2D<T> >
*GetWrkSpFileReaderHandler(Loki::Type2Type<LArr2D<T> >,int ContextID,const std::wstring& fn, int YSz,int Metric,const std::wstring& DirN);
//
template<typename T>
Zeus::GenCollWriter<LArr2D<T> >
*GetGenCollFileWriterHandler(Loki::Type2Type<LArr2D<T> >,int ID,const std::wstring& fname,const std::wstring& DirN,int AxSz,double PixSz,double PatchCx,double PatchCy,double XMax,double XMin,ContoursOutAux* aux=NULL);
//
GenCollWriter<CatalogueFormatType>
*GetCatWriterHandler(int ContextID, const std::wstring& fn,const std::wstring& GrpName,int NColumns);
//
GenCollWriter<OutputFormatType>
*GetFormatCatWriterHandler(int ContextID,int DetectionType,const std::wstring& fn,const std::wstring& GrpName);
//
GenCollReader<PeakCollReadbleType>
*GetObjReaderHandler(int ContextID,const std::wstring& Dir,const std::wstring& fn);
//
GenCollWriter<PeakCollType>
*GetObjWriterHandler(int ContextID,const std::wstring& Dir,const std::wstring& fn,int InitObjID);
//
GenCollReader<PatchGeomType>
*GetPatchGeomFileReaderHandler(Loki::Type2Type<PatchGeomType>,int ContextID,const std::wstring& DirName,const std::wstring& fn);
//
GenCollWriter<PatchGeomType>
*GetGenCollFileWriterHandler(Loki::Type2Type<PatchGeomType>,int ID,const std::wstring& DirName,const std::wstring& fname);
//
GenCollReader<CatalogueFormatType>
*GetGenCollFileReaderHandler(Loki::Type2Type<CatalogueFormatType>,int ID,const std::wstring& Dname,const std::wstring& fname);
//
GenCollReader<QA_CltLstCatalogueType>
*GetQA_ColLstReaderHandler(int ID,const std::wstring& Dname,const std::wstring& fname);
//
GenCollReader<QA_ProfilesLstType>
*GetQA_ProfilesLstReaderHandler(int ID,const std::wstring& Dname,const std::wstring& fname);
//
template<typename T>
GenCollReader<LArr1D<T> >
*GetVectFileHandlerReader(Loki::Type2Type<LArr1D<T> >,int ContextID,const std::wstring& fn, int Sz,const std::wstring& DirN);
//
template<typename T>
Zeus::GenCollWriter<LArr1D<T> >
*GetVectFileHandlerWriter(Loki::Type2Type<LArr1D<T> >,int ID,const std::wstring& fname,const std::wstring& DirN);
//
//

} // end of namespace ZEUS



#endif //ZEUS_INOUTMANIP 
