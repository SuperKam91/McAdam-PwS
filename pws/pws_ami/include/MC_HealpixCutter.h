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
#ifndef MCHEALPIXCUTTER
#define MCHEALPIXCUTTER

#include "ZEUS_WorkSpace.h"
#include "ZEUS_InOut.h"
#include "MC_HPixMapPatchCutting.h"

//---------------------------------

#define MAXBASEPIXEDGE				1.50
#define PANICFACTOR					1.05
#define MC_PATCHSZLARGE				256
#define MC_PATCHSZSMALL				128
#define MC_PATCHSZTINY				64
#define HEALPIX_BADPIX				-1e20
//const double Healpix_undef=-1.6375e30;
#define DVDCORRUPTED_BADPIX			1e20
#define PSDEFPATCHSIZE				256
#define MAPEPOCH					2000.0
#define MAXAGGRRADIUS				0.00145444104332646688102977729334
 //5 arcmin -> RAD

enum BasePixels					{BP_PON0 = 0x00,BP_PON1 = 0x01,BP_PON2 = 0x02,BP_PON3 = 0x03,
								BP_EQ0 = 0x04,BP_EQ1 = 0x05,BP_EQ2 = 0x06,BP_EQ3 = 0x07,
								BP_POS0 = 0x08,BP_POS1 = 0x09,BP_POS2 = 0x0A,BP_POS3 = 0x0B,BP_EMPTY = 0x0C};

enum	BasePixCorners			{BPC_DOWN=0x00,BPC_RIGHT=0x01,BPC_LEFT=0x02,BPC_UP=0x03};
enum	DirBasePixels			{BPD_DOWN=0,BPD_LEFT=2,BPD_UP=4,BPD_RIGHT=6};

enum	PatchDirections			{PTCHDIR_NAT=0,PTCHDIR_ROT180=1,PTCHDIR_ROT90=2,PTCHDIR_ROT270=3,PTCHDIR_ROTEND=4};


inline bool SelectOnCollLst(const Zeus::CatLineCollType::value_type& t)
{return ((t.CollLstIndex_ < 0)?true:false);}

inline bool SortColLstByIndex(const Zeus::CatLineCollType::value_type& first,const Zeus::CatLineCollType::value_type& sec)
{return first.CollLstIndex_ < sec.CollLstIndex_ ;}

inline bool SortPSCatBySNR(const Zeus::NonBlingCatLineType& first,const Zeus::NonBlingCatLineType& sec)
{
	return first.SNR_ < sec.SNR_;
}

inline bool SortChi2ByIdFreq(const Zeus::CatalogueFormatType::StorageType::value_type& first,
							 const Zeus::CatalogueFormatType::StorageType::value_type& sec)
{
	if(first.ID_ != sec.ID_)
		return first.ID_ < sec.ID_;

	return first.freqId_ < sec.freqId_;
}
//
class	GeneralMapCutter
{
	GeneralMapCutter(const GeneralMapCutter&);
	GeneralMapCutter& operator=(const GeneralMapCutter&);

	struct Chi2InfoType
	{
		int					freqId_;
		std::wstring		fileName_;

		Chi2InfoType(int freqId,const std::wstring& fileName)
			:freqId_(freqId),fileName_(fileName)
		{}
	};


	void					ReportCuttingCurrentParams(void) const;

	void 					WritePtgData2Output(const pointing& CentralPixel) const;
	void					WriteGeomData2Output(void) const;
	void					PutCurrGeomData2Place(void);
	std::wstring			GetCurrPatchPtgsOutputName(int PatchNumber,int ptgType) const;
	Zeus::LArr2D<double>&	PutPtgs2Workspace(int ptgAngle,Zeus::LArr2D<double>& ws) const;
	void					SetMap(void);
	void					ReadMaskMap(const std::wstring& fname,Healpix_Map<HEALPIX_ATOM_PREC>& Mask);
	inline int				PtgOutMask(void) const
	{
		pointing	temp(CurrCentralPtg_);

		return CheckCentralPixel(temp);
	}

	void					TranslatingPtgs(void);
	void					AdjustRejectionMask(void);
	void					CreateObservedPixels(void);
	void					ReadPtgsFileIn(int PatchNumber, Zeus::LArr2D<double>& theta,Zeus::LArr2D<double>& phi);
	void					WriteOutCats(const std::wstring& fname,const Zeus::OutputFormatType& cat);
	int						ParseChi2InCatNames(const std::wstring& InCatCollNames,std::vector<Chi2InfoType>& fNames);
	int						ReadChi2CatIn(const std::vector<Chi2InfoType>& fNames,Zeus::CatalogueFormatType& Data);
	void					DoMakeChi2Catalogue(Zeus::CatalogueFormatType& Data,Zeus::OutputFormatType& catOut);
	void					DoConvertCatalogue(const Zeus::CatalogueFormatType& catIn,Zeus::OutputFormatType& catOut);
	void					Chi2OneRow(Zeus::CatalogueFormatType::StorageType::const_iterator& f,
									 Zeus::CatalogueFormatType::StorageType::const_iterator& s,
									 Zeus::OutputFormatType::StorageType::value_type& t
									 );
	void					GetFreqIdFields(int FreqId,double*& ptrSNR,double*& ptrFlux,Zeus::ExtFormatLnType& t);
	void					SetMainCatFields(Zeus::CatalogueFormatType::StorageType::const_iterator& f,Zeus::CatLineType& t, double BayFlux);
	int						GetFreqID(const std::wstring& fname);
protected:
	int									CurrBase_;
	int									CurrPatchNumber_;
	const EnvIDsType					InOutEnvID_;
	const std::wstring					DirOut_;
	const std::wstring					DirIn_;
	const std::wstring					DataBuffer_;
	const std::wstring					MaskRejectFile_;
	const std::wstring					Chi2MaskFile_;
	Zeus::LArr2D<pointing>				Ptgs_;
	pointing							CurrCentralPtg_;
	pointing							CurrSrcPtg_;
	Healpix_Map<HEALPIX_ATOM_PREC>		HPixMap_;
	Healpix_Map<HEALPIX_ATOM_PREC>		HPixMask_;
	Healpix_Map<HEALPIX_ATOM_PREC>		pixCovered_;
	Zeus::PatchGeomType					PatchsGeom_;
	double								PixelSz_;
	double								Spin_;
	double								NonBlind_Snr_;
	double								PredictedRadius_;
	double								PredictedFlux_;
	double								RadiusBay_;
	double								FluxBay_;
	double								ErrRadius_;
	double								ErrFlux_;
	double								ErrPos_;
	double								ErrRadiusHigh_;
	double								ErrRadiusLow_;
	double								ErrFluxHigh_;
	double								ErrFluxLow_;
	double								IniRejectPercent_;
	double								QAIN_CyR500;
	double								QAIN_T500;
	double								QAIN_detectable;
	pointing							QAIN_SrcPtg_;
	int									NonBlindSrcIndex_;
	int									Objects2buffer_;
	const double						PercentReject_;
	const coordsys						PtgsCoordSystem_;

	void								NormaliseCoords(Zeus::NonBlingCatType& SrsCat);
	void								NormalisePatchCoords(Zeus::NonBlingCatType& SrsCat);

	void								TranslateCoords(Zeus::NonBlingCatType& SrsCat,const coordsys	PtgsCoordSystem);
	int									ReadCatalogue(const std::wstring& fname,Zeus::CatalogueFormatType& cat);
	int									ReadCatalogue2NBlind(const std::wstring& SrsCatFName,Zeus::NonBlingCatType& SourceCat);

	virtual void			do_Initialize(void) = 0;
	virtual int				do_SetNewPatchCentralPixel(void) = 0;
	virtual int				CheckCentralPixel(const pointing&	temp) const  = 0;
	virtual	void			CentralPtgOnMask(void);
	virtual	int				CheckMinimumGoodPix(void);
	virtual void			ShufflePatchNumbers(Zeus::PatchGeomType& pg)
	{
		Zeus::PatchGeomType::StorageType::iterator			piv(pg.Storage_.begin());
		Zeus::PatchGeomType::StorageType::const_iterator	const end(pg.Storage_.end());

		for(;piv != end;++piv)
		{piv->SrcIndex_ = piv->PatchNumber_;}
		std::sort(pg.Storage_.begin(),pg.Storage_.end());
	}

	
	inline void				errInvParam(int errCode,const wchar_t* msg,const wchar_t* ParamName) const
	{
		std::wstring errstring(msg);

		errstring += std::wstring(L" -> ");
		errstring += ParamName;
		throw Zeus::libException(errCode,errstring,*this);
	}

	// Default projection is Gnomonic
	virtual void			do_Projection(void);

public:
	GeneralMapCutter(EnvIDsType InOutEnvID,const Zeus::PatchGeomType::HeaderType& Header,const std::wstring& DirOut,
		const std::wstring& MaskRejectFile,const std::wstring& Chi2MaskFile,const std::wstring& DirIn,coordsys PtgsCoordSystem,double PercentReject,
		const std::wstring& DataBuffer,int Object2buffer)
		: InOutEnvID_(InOutEnvID),DirOut_(DirOut),MaskRejectFile_(MaskRejectFile),Chi2MaskFile_(Chi2MaskFile),
		DirIn_(DirIn),PtgsCoordSystem_(PtgsCoordSystem),PercentReject_(PercentReject),
		DataBuffer_(DataBuffer),Objects2buffer_(Object2buffer),
		NonBlind_Snr_(0.0),NonBlindSrcIndex_(0),PredictedRadius_(0.0),PredictedFlux_(0.0),
		ErrRadius_(0.0),ErrFlux_(0.0),ErrPos_(0.0),CurrSrcPtg_(-2.0e10,-2.0e10),QAIN_SrcPtg_(-2.0e10,-2.0e10),QAIN_CyR500(SZCATALOGUEDEFAULTVALUE),
		QAIN_T500(SZCATALOGUEDEFAULTVALUE),QAIN_detectable(1.0),ErrRadiusHigh_(SZCATALOGUEDEFAULTVALUE),ErrRadiusLow_(SZCATALOGUEDEFAULTVALUE),
		ErrFluxHigh_(SZCATALOGUEDEFAULTVALUE),ErrFluxLow_(SZCATALOGUEDEFAULTVALUE),FluxBay_(SZCATALOGUEDEFAULTVALUE),RadiusBay_(SZCATALOGUEDEFAULTVALUE)

	{
		PatchsGeom_.Header_ = Header;
	}

	void			Initialize(void);
	void			MakePointingsFiles(void);
	void			CreatHealpixMapWithDetections(const std::wstring& CatFname,const std::wstring& HealpixFname);
	void			ConvertCatalogue(const std::wstring& InCatName,const std::wstring& OutCatName);
	void			MakeChi2Catalogue(const std::wstring& InCatCollNames,const std::wstring& OutCatName);
	void			RemoveChi2Mask(Zeus::OutputFormatType& cat);
	void			RemoveSuspicious(Zeus::OutputFormatType& cat);
	virtual			~GeneralMapCutter(void) = 0;
};


class	HealpixCutter : public GeneralMapCutter
{
	double						ColumnMetricDist_;
	double						RowMetricDist_;
	int							SqrtNPatches_;
	int							BaseRowPix_;
	int							PatchLine_;
	int							PatchColumn_;
	DirBasePixels				RowDir_;
	DirBasePixels				ColumnDir_;
	int							Healpix00Pix_;

	//
	inline int		PixAt(int HealPixPix,int dist,DirBasePixels dirT) const
	{
		fix_arr<int,8>	res;

		for(int i=0;i<dist;++i)
		{
			HPixMap_.neighbors(HealPixPix,res);
			HealPixPix = res[dirT];
		}
		return HealPixPix;
	}

	//
	int				PixAtDist(int HealPixPix,double dist,DirBasePixels dirT) const;
	//
	int				PixAtBorder(int HealPixPix,DirBasePixels dirT) const;
	//
	inline int		GetNestedCurrBasePix(int pix) const
	{return (pix >> (HPixMap_.Order() << 1));}
	//
	inline int		GetNestedBasePixCorner(BasePixels BasePix,BasePixCorners corner) const
	{
		int t_bp(BasePix),t_cor(corner);
		for(int i=0;i < HPixMap_.Order();++i)
		{t_bp <<= 2;t_bp |= t_cor;}
		return t_bp;
	}

	//
	double			EvalSpin(int CCorner,double dist) const;
	inline double	GetPixDistance(int pix1,int pix2) const
	{
		vec3	t_v0(HPixMap_.pix2ang(pix1));
		vec3	t_v1(HPixMap_.pix2ang(pix2));
		return std::acos(dotprod(t_v0,t_v1)/(t_v0.Length()*t_v1.Length()));
	}

	//
	inline double	Dist2Border(int pix,DirBasePixels dirT) const
	{
		int PixBorder(PixAtBorder(pix,dirT));
		return GetPixDistance(pix,PixBorder);
	}

	//
	void			ReportCuttingCurrentParams(void) const;
	//
	int				SetNewBasePixelState(void);

protected:
	virtual void	do_Initialize(void);
	virtual int		do_SetNewPatchCentralPixel(void);
	virtual int		CheckCentralPixel(const pointing&	temp) const
	{

		const double tLatLimit(PatchsGeom_.Header_.GalacticCut_ / RAD2DEGREE);

		if( 
			((temp.theta < 0.0)	|| (temp.phi < 0.0)) ||
			((temp.theta > ((PIOVER2 - tLatLimit) + 1.0e-8)) &&
			(temp.theta < ((PIOVER2 + tLatLimit) - 1.0e-8)))
			)
			return 0;

		return 1;		
	}

public:
	HealpixCutter(EnvIDsType InOutEnvID,const Zeus::PatchGeomType::HeaderType& Header,const std::wstring& DirOut,
		const std::wstring& MaskRejectFile,const std::wstring& Chi2MaskFile,
		const std::wstring& DirIn,coordsys PtgsCoordSys,double PercentReject,
		const std::wstring& DataBuffer,int Object2buffer)
		: GeneralMapCutter(InOutEnvID,Header,DirOut,MaskRejectFile,Chi2MaskFile,DirIn,PtgsCoordSys,PercentReject,
		DataBuffer,Object2buffer)
	{}
	//
	virtual ~HealpixCutter(void)
	{}
	//
};

class	ConstGalacticLatCutter : public GeneralMapCutter
{
	double			HalfPatchSz_;
	int				NorthSouth_;
	double			CurrCoLat_;
	double			LatIncr_;

	double			CurrLong_;
	double			LongIncr_;

	void			SetRingsProps(void);

protected:
	virtual void	do_Initialize(void);
	virtual int		do_SetNewPatchCentralPixel(void);
	virtual void	ShufflePatchNumbers(Zeus::PatchGeomType& pg);
	virtual int		CheckCentralPixel(const pointing&	temp) const
	{

		const double tLatLimit(PatchsGeom_.Header_.GalacticCut_ / RAD2DEGREE);

		if( 
			((temp.theta < 0.0)	|| (temp.phi < 0.0)) ||
			((temp.theta > ((PIOVER2 - tLatLimit) + 1.0e-8)) &&
			(temp.theta < ((PIOVER2 + tLatLimit) - 1.0e-8)))
			)
			return 0;

		return 1;	
	}

public:
	ConstGalacticLatCutter(EnvIDsType InOutEnvID,const Zeus::PatchGeomType::HeaderType& Header,const std::wstring& DirOut,
		const std::wstring& MaskRejectFile,const std::wstring& Chi2MaskFile,
		const std::wstring& DirIn,coordsys PtgsCoordSys,double PercentReject,
		const std::wstring& DataBuffer,int Object2buffer)
		: GeneralMapCutter(InOutEnvID,Header,DirOut,MaskRejectFile,Chi2MaskFile,DirIn,PtgsCoordSys,
		PercentReject,DataBuffer,Object2buffer)
	{}
	virtual ~ConstGalacticLatCutter(void)
	{}
};

class	NonBlindCutter : public GeneralMapCutter
{
	const std::wstring		NonBlindPtgsFile_;
	Zeus::NonBlingCatType	SrsCat_;
	double					R500Ratio_;
	double					YtotAux_;
	double					tAlpha_;
	Zeus::GammaFuncts		f_;
	Zeus::NonBlingCatType::StorageType::const_iterator	pivSrcCat_;
	Zeus::NonBlingCatType::StorageType::const_iterator	endSrcCat_;
//
	int						ReadQA_ColLst(const std::wstring& fname,Zeus::QA_CltLstCatalogueType& cat);
//
	int						ReadQA_Profiles(const std::wstring& fname,Zeus::QA_ProfilesLstType& cat);
// 
	int						XtendNB_Info(Zeus::QA_CltLstCatalogueType& cat,double AggRadius /*arcmin*/);
//
	int						MatchProfileQA(Zeus::QA_ProfilesLstType& proflst);
//
	int						SingleProfileQA(void);
//
	int						FindClosestDetection2NB(double AggRadius /*arcmin*/,const Zeus::QA_CltLstCatalogueType::StorageType::value_type& t,Zeus::NonBlingCatType::StorageType::iterator& out);
//
	inline int				ReadCatalogueIn(void)
	{
		return ReadCatalogue2NBlind(NonBlindPtgsFile_,SrsCat_);
	}

	inline	double ytotFunct(double x)
	{
		return f_.betainc(YtotAux_ / (std::pow(std::cos(x), tAlpha_) + YtotAux_)) * std::cos(x);
	}

//
protected:
	virtual void			do_Initialize(void);
	virtual int				do_SetNewPatchCentralPixel(void);
	virtual int				CheckCentralPixel(const pointing&	temp) const;
	virtual void			CentralPtgOnMask(void);
	virtual	int				CheckMinimumGoodPix(void);
	virtual void			ShufflePatchNumbers(Zeus::PatchGeomType& pg);


public:
	NonBlindCutter(EnvIDsType InOutEnvID, double R500Ratio, const Zeus::PatchGeomType::HeaderType& Header, const std::wstring& DirOut,
		const std::wstring& NonBlindPtgsFile,const std::wstring& MaskRejectFile,const std::wstring& Chi2MaskFile,
		const std::wstring& DirIn,
		coordsys PtgsCoordSys,double PercentReject,const std::wstring& DataBuffer,int Object2buffer)
		: GeneralMapCutter(InOutEnvID,Header,DirOut,MaskRejectFile,Chi2MaskFile,DirIn,PtgsCoordSys,
							PercentReject,DataBuffer,Object2buffer),
							NonBlindPtgsFile_(NonBlindPtgsFile), R500Ratio_(R500Ratio)
	{}
	virtual ~NonBlindCutter(void)
	{}
};
//
inline bool RemoveSuspiciousFunct(const Zeus::OutputFormatType::StorageType::value_type& p)
{
	if(
#ifdef APPLY_SNRCUT
		(p.Cat_.NormalAmpl_ < APPLY_SNRCUT) ||
#endif
		(p.Ext_.CHI2_ > 16.0) ||
		((p.Cat_.NormalAmpl_ < 15.0) && (((p.Ext_.SNRF_ - p.Cat_.NormalAmpl_)/p.Cat_.NormalAmpl_) >= 1.0))
		)
		return true;

	return false;
}
//
struct Chi2UnderMaskRemFunct
{
	inline bool operator()(const Zeus::OutputFormatType::StorageType::value_type& p) const
	{
		if(HPixMask_[HPixMask_.ang2pix(pointing((90.0-p.Cat_.GalLatDegs_)*Conv2Rads_,p.Cat_.GalLongDegs_*Conv2Rads_))] < 0.5)
			return true;

		return false;
	}
//
	Chi2UnderMaskRemFunct(const Healpix_Map<HEALPIX_ATOM_PREC>& HPixMask)
		:HPixMask_(HPixMask),Conv2Rads_(PI/180.0)
	{}

	const Healpix_Map<HEALPIX_ATOM_PREC>&	HPixMask_;
	const double							Conv2Rads_;
};

inline bool CltLst_NotMatchFunctor(const Zeus::QA_CltLstCatalogueType::StorageType::value_type& v)
{
		if(v.Cat_SNR< 0.5)
			return true; // delete
		return false;
}

inline bool SortNonBlindByColat(const Zeus::NonBlingCatType::StorageType::value_type& first,const Zeus::NonBlingCatType::StorageType::value_type& sec)
{
	return first.ptg_.theta < sec.ptg_.theta;
}

inline bool SortNonBlindByIndex(const Zeus::NonBlingCatType::StorageType::value_type& first,const Zeus::NonBlingCatType::StorageType::value_type& sec)
{
	return first.Index_ < sec.Index_;
}

inline bool NonBlindNotFinalCatFunctor(const Zeus::NonBlingCatType::StorageType::value_type& v)
{
	if(v.flagged_<0)
			return true; // delete
	return false;
}

//void	CreatHealpixFileWithPatches(int NSide,int PatchSz,int NPatches);

#endif //MCHEALPIXCUTTER

