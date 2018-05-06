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



//---------------------------
#ifndef MCPIXEXTRACTPROC
#define MCPIXEXTRACTPROC

#include <memory>
#include "trafos.h"
#include "ZEUS_InOut.h"
#include "ZEUS_FourierMachine.h"
#include "ZEUS_PhysicsMath.h"
#include "ZEUS_GaussianRandomGen.h"
#include "MC_Defaults.h"
#include "MC_Antenna.h"
#include "MC_ObjFactory.h"

//---------------------------

#define	MASKBADPIXGROWFACTOR	0.0
#define	SMOOTHFACTOR			0.25
#define	MASKCLAMPLEVEL			0.75
#define	MASKINGRMSFACTOR		10.0

#define	SRCMAPPREFIX			L"Src_"
#define	MASKMAPPREFIX			L"Mask_"
#define	NEWMAPPREFIX			L"New_"

struct	AntennaInfoType
{
	MC_antenna::AntennaIDType		Id_;
	double							Freq_;
	MC_antenna*						Antenna_;
};

class PixExtractProc
{
	typedef	std::vector<AntennaInfoType>	AntennasCollType;

public:
	static void	SetMapStaticInfo(void);

	struct	StatInfoType
	{
		double		CurrPatchRms_;
		double		CurrPatchBias_;
		double		PercentBadpixAfterRot_;
		int			NPixIncluded_;

		void		Reset(void)
		{CurrPatchRms_ = CurrPatchBias_ = PercentBadpixAfterRot_ = 0.0;NPixIncluded_ = 0;}
	};


	PixExtractProc(EnvIDsType InOutEnvID,PixExtractProcCtorArgsType& args)
		:InOutEnvID_(InOutEnvID),CurrPatchNumber_(0),TotalPatches_(0),Trafo_(0),NoiseGen_(0),FourierMachine_(0),OpMode_(0)
	{
		Epoch_			= args.Epoch_;
		RmsRejectLevel_ = args.RmsRejectLevel_;
		HP_Dir_			= args.HP_Dir_;
		Masks_Dir_		= args.Masks_Dir_;
		Masks_DirIn_	= args.Masks_DirIn_;
		DirOut_			= args.DirOut_;
		DirOutMaps_		= args.DirOutMaps_;
		DataBuffer_		= args.DataBuffer_;
		Data2Buffer_	= args.Data2Buffer_;

		MaskRejectName_	= args.MaskRejectName_;
		MaskRemoveName_	= args.MaskRemoveName_;
		MaskEnlagThreshold_ = args.MaskEnlagThreshold_;
		MaskSmoothFWHM_	= args.MaskSmoothFWHM_;
		MaskRadiusRatio_ = args.MaskRadiusRatio_;
		syncID_			= args.syncID_;
		FileCollection_.swap(args.FileCollection_);
		// 50331648 = 12 * 2048^2
		PixSz_			= std::sqrt(PITIMES4 / 50331648.0);
		FillInDefaultValues(FileCollection_);
	}

	void	DoPixExtractionAllMaps(int OpMode,long FPatch,long LPatch);
	void	Extract1stMap(const std::wstring& fName);
	void	JoinManyMasks(void);
	void	CreatHealpixMapWithPtgs(const std::wstring& fName);
	void	SubSourcesDoProcUpAllMaps(void);

	~PixExtractProc(void)
	{
		delete Trafo_;
		delete NoiseGen_;
		delete FourierMachine_;
		AntennasCollType::iterator			piv(AntColl_.begin());
		AntennasCollType::const_iterator	const end(AntColl_.end());
		for(;piv != end;++piv)
		{delete piv->Antenna_;}
	}
private:
	PixExtractProc(const PixExtractProc&);
	PixExtractProc& operator=(const PixExtractProc&);

	void						NormaliseCoords(Zeus::NonBlingCatType& SrsCat);
	void						TranslateCoords(Zeus::NonBlingCatType& SrsCat,const coordsys	PtgsCoordSystem);

	inline	void				InitFourierMachine(void)
	{
		delete FourierMachine_;
		FourierMachine_ = new Zeus::FourierMachine(PtchSz_);
		FourierMachine_->Initialize();	
	}
//
	void			ReadSrcCatalogue(const std::wstring& SrsCatFName,Zeus::NonBlingCatType& SourceCat);
//
	void			ReadBadPixMask(const std::wstring& DirName);
//
	inline void		ClampLevelMask(Healpix_Map<HEALPIX_ATOM_PREC>& Mask)
	{
		const int	NPix(Mask.Npix());

		for(int i=0;i<NPix;++i)
		{
			if(Mask[i] < MASKCLAMPLEVEL)
			{Mask[i] = 0.0;}
			else
			{Mask[i] = 1.0;}
		}
	}
//
	inline Healpix_Map<HEALPIX_ATOM_PREC>&	MergeMapsInPlace(const Healpix_Map<HEALPIX_ATOM_PREC>& Mask)
	{
		const int	NPix(IllPixMask_.Npix());

		for(int i=0;i<NPix;++i)
		{if(Mask[i] < MASKCLAMPLEVEL) IllPixMask_[i] = 0.0;}
		return IllPixMask_;
	}
//
	inline int			EvalMaskRadiusPix(const RealPlaneSurfType& Patch)
	{

		const double	*piv(Patch.GetInnerData().begin());
		const double	MaxValue(*piv);
		int				Xsz,Ysz,i;

		Patch.GetSz(Ysz,Xsz);
		Xsz >>= 1;

		for(i=0;i<=Xsz;++i)
		{
			if(MaxValue > (2.0 * piv[i]))
				return i;
		}

		return i;
	}
//		
	Healpix_Map<HEALPIX_ATOM_PREC>&	SubAddMapsInPlace(Healpix_Map<HEALPIX_ATOM_PREC>& Accumul,Healpix_Map<HEALPIX_ATOM_PREC>& SrcsMap,int AddSub);
	void					ReadPtgsFileIn(int PatchNumber, Zeus::LArr2D<double>& theta,Zeus::LArr2D<double>& phi);	
	void					MakeSrcsMaskMaps(Healpix_Map<HEALPIX_ATOM_PREC>& TempMask,Healpix_Map<HEALPIX_ATOM_PREC>& SrcsMap,Zeus::NonBlingCatType& SourceCat);
	void					InitialiseSrcSub(Zeus::NonBlingCatType& SourceCat);
	void					FillInDefaultValues(PixExtractProcCtorArgsType::MapFileInfoCollType& MapColl);
	void					Initialize(void);
	void					DoPixExtraction1Map(long LPatch);
	void					NormaliseRMS(Zeus::LArr2D<double>& PixDest);
	void					ReplaceBadPixSubBias(void);
	void					AddNoise2Pix(Zeus::LArr2D<double>& Pixels,Zeus::LArr2D<double>& PixMask);
	void					NormalisePix(double factor,Zeus::LArr2D<double>& DestPix);
	void					OpenHealpixMap(const std::wstring& DirName,const std::wstring& MapName,Healpix_Map<HEALPIX_ATOM_PREC>& map);
	int						ReadMaskMap(const std::wstring& dirName,const std::wstring& fname,Healpix_Map<HEALPIX_ATOM_PREC>& Mask,std::wstring& ActualFileName);
	inline std::wstring		GetCurrPatchPtgsOutputName(int patchN,int ptgType) const
	{
		std::wstring	ext(ptgType?PHIFILEEXT:THETAFILEEXT);

		return std::wstring(PTGSPREFIX) + Zeus::PutNumber2Txt(GeomNSide_) + std::wstring(L"_") + Zeus::PutNumber2Txt(patchN) + ext;
	}

	inline std::wstring		GetCurrPatchOutputName(int patchN) const
	{
		return GetFileTypePrefix() + Zeus::PutNumber2Txt(GeomNSide_) + std::wstring(L"_") + Zeus::PutNumber2Txt(CurMapPiv_->FreqID_) + std::wstring(L"_") + Zeus::PutNumber2Txt(patchN);
	}

	void					ProcCoordSys(const Zeus::HealpixHeaderAtomType& key);
	void					ProcUnits(const Zeus::HealpixHeaderAtomType& key);
	void					ProcFreq(const Zeus::HealpixHeaderAtomType& key);

	inline void				PrintWarning(const std::wstring& msg) const
	{
		std::wstring temp(msg);
		temp += std::wstring(HEALPIXMAP_WARNING) + HP_Dir_ + CurMapPiv_->HP_FName_;
		(Zeus::ConManager::Instance())->PrintStr2Console(temp);
	}
	inline void				PrintProcessingFile(void)
	{
		std::wstring temp(L"\n");
		temp	+= (std::wstring(PROCESSING_FILE) + HP_Dir_ + CurMapPiv_->HP_FName_ + std::wstring(L"\n\n"));
		(Zeus::ConManager::Instance())->PrintStr2Console(temp);
	}

	inline	static			std::wstring GetCoordSysMeg(int coordsys)
	{
		switch(coordsys)
		{
		case Galactic:
			return std::wstring(HEALPIX_COORDS_GAL);
		case Equatorial:
			return std::wstring(HEALPIX_COORDS_EQU);
		case Ecliptic:
			return std::wstring(HEALPIX_COORDS_ECL);
		default:
			return std::wstring();
		}
	}
	inline	void			PrintCoordSys(int coordsys) const
	{
		std::wstring	msg(HEALPIX_COORDSYS_MSG);
		msg	+= GetCoordSysMeg(coordsys);
		msg += std::wstring(L"\n\n");
		(Zeus::ConManager::Instance())->PrintStr2Console(msg);
	}

	inline	void			PrintCoordSysMismatch(int coordsys1) const
	{
		std::wstring	msg(WARMSG_MC_COORDSYSMIS);
		msg	+= GetCoordSysMeg(coordsys1);
		msg += std::wstring(WARMSG_MC_COORDSYSOUT);
		msg	+= GetCoordSysMeg(Galactic);
		PrintWarning(msg);
	}

	inline  void			errGeomFile(int errCode,const wchar_t* msg) const
	{
		std::wstring errstring(msg);

		errstring += std::wstring(L" -> ");
		errstring += DirOut_;
		throw Zeus::libException(errCode,errstring,*this);
	}

	inline static void			errInvalidFreq(int freq)
	{
		std::wstring errstring(ERRMSG_MC_INVALIDFREQ);

		errstring += std::wstring(L" -> ");
		errstring += Zeus::PutNumber2Txt(freq);
		throw Zeus::libException(ERRCOD_MC_INVALIDFREQ,errstring,std::wstring(L"PixExtractProc::errInvalidFreq"));
	}

	inline  void			errObjectNotFound(int errCode,const std::wstring& DirName,const wchar_t* msg,const std::wstring& FileName,int throwEx) const
	{
		std::wstring errstring(ERROR_COD_MC_MSG);
		
		errstring += Zeus::PutNumber2Txt(errCode);
		errstring += std::wstring(L"; ");
		errstring += DirName;errstring += FileName;
		errstring += std::wstring(L"\n");
		errstring += msg;
		if(throwEx)
		{
			throw Zeus::libException(errCode,errstring,*this);
		}
		else
		{
			(Zeus::ConManager::Instance())->PrintStr2Console(errstring + std::wstring(L"\n\n"));
		}
	}

	inline  void			errHealpix(int errCode,const wchar_t* msg,const std::wstring& DirName,const std::wstring& FileName,int throwEx) const
	{
		errObjectNotFound(errCode,DirName,msg,FileName,throwEx);
	}


	void					errHealpixNside(int errCode,wchar_t* msg, int NSideGeo,int NSideMap) const;
	void					GetMapUnits(const std::string& unitStr);
	void					PrintMapUnits(void) const;
	void					ReadPtgsFile1Coord(int coord,const std::wstring& fname);
	inline void				ReadPtgsFile(int patchN)
	{
		//theta first
		ReadPtgsFile1Coord(0,GetCurrPatchPtgsOutputName(patchN,0));
		//phi next
		ReadPtgsFile1Coord(1,GetCurrPatchPtgsOutputName(patchN,1));
	}

	void					ComputeStats(void);
	void					InvertPatch(void);
	void					PutInPtgs(int coord,const Zeus::LArr2D<double>&	ws);
	void					PrintBeamSigma(double orgFWHM,double newFWHM) const;
	void					ReportPixCurrentParams(void) const;
	inline void				PrintFreq(void) const
	{
		std::wstring msg(HEALPIX_FREQ_MSG);

		msg += Zeus::PutNumber2Txt(CurMapPiv_->Freq_);
		(Zeus::ConManager::Instance())->PrintStr2Console(msg + std::wstring(L"\n\n"));
	}

	inline void				PrintDirErr(const std::wstring& temp) const
	{
		std::wstring msg(WARMSG_MC_CREATEDIR);
		msg += temp;
		PrintWarning(msg);
	}

	inline void				PrintFileType(const std::wstring& temp) const
	{
		std::wstring msg(MT_FILENAME_TXT);
		msg += temp;
		(Zeus::ConManager::Instance())->PrintStr2Console(msg + std::wstring(L"\n\n"));
	}
	void					PrintObsMapLowRes(int CurrNSide,int OperaNSide) const;
	double					PixelInterpolation(const pointing& ptg0);
	void					Ptgs2Pixels(const Healpix_Map<HEALPIX_ATOM_PREC>& BadPixMask);
	void					ReadGeoProps(int JustHeader);
	template<typename T>
	void 					WritePixels2File(Zeus::LArr2D<T>& data)
	{

		pointing CentralPtg(PatchGeom_.Storage_.at(CurrPatchNumber_).X0Y0Ptg_);
		std::auto_ptr<Zeus::GenCollWriter<Zeus::LArr2D<T> > > 
			FWriter(Zeus::GetGenCollFileWriterHandler(Loki::Type2Type<Zeus::LArr2D<T> >(),
			Data2Buffer_?1000:InOutEnvID_,GetCurrPatchOutputName(CurrPatchNumber_),DirOutCurrMap_,
			data.getPtrMetric(),PixSz_,CentralPtg.phi,CentralPtg.theta,
			static_cast<double>(data.getPtrMetric()) * PixSz_,0.0));

		FWriter->Initialize();
		FWriter->Write(data);
		FWriter->Flush();
	}
	void					SetNoiseGenerator(void);
	int						GetMaxNSide(void);
	void					AntennaSmooth(const Zeus::LArr2D<double>& OrgPix,Zeus::LArr2D<double>& DestPix,double Freq,double FWHM);
	void					SmoothMask(Zeus::LArr2D<double>& Mask);
	double					CombinePatchesX(Zeus::LArr2D<double>& NewMap,Zeus::LArr2D<double>& Mask);
	void					RenormaliseMask(Zeus::LArr2D<double>& DestPix);
	void					MakeMask(Zeus::LArr2D<double>& NewMask);
	int						CheckBadPix(void);
	double					PixTry1Rotation(PatchDirections dir,const Zeus::LArr2D<double>& Mask);
	double					TryBestRotation(PatchDirections& BestRotation,const Zeus::LArr2D<double>& Mask);
	double					MakeNewPatchByRotation(void);
	void					PreProcessMap(void);
	void					MaskMap(void);
	void					WriteMasks(void);
	int						CheckBadPixels(void);
	void					SetBadPix(int Pix,Healpix_Map<HEALPIX_ATOM_PREC>& Map);
	void					MarkBadPixels(Healpix_Map<HEALPIX_ATOM_PREC>& Mask);
	void					CombineMaps(const Healpix_Map<HEALPIX_ATOM_PREC>& OutOrg , const Healpix_Map<HEALPIX_ATOM_PREC>& Mask,Healpix_Map<HEALPIX_ATOM_PREC>& Result);
	void					FilterMap(const Healpix_Map<HEALPIX_ATOM_PREC>& InMap,Healpix_Map<HEALPIX_ATOM_PREC>& OutMap,double fwhm);
	void					CreateDir(void);
	std::wstring			GetFileTypeDir(void);
	std::wstring			GetFileTypePrefix(void) const;
	void					PixDoRotation(PatchDirections dir,Zeus::LArr2D<double>& Pixels);
	void					ProjectOntoSphere(bool mask,RealPlaneSurfType& Patch,Healpix_Map<HEALPIX_ATOM_PREC>& Map,const pointing& pos);
	void					Project1Pix(double Value,const pointing& ptg,Healpix_Map<HEALPIX_ATOM_PREC>& Map);
	void					EvalValueForMasking(double& Value,const pointing& tptg);

	MC_antenna*				AntennaFactory(MC_antenna::AntennaIDType Id,double Freq, double FWHM);
	AntennasCollType::iterator	CreateAntennas(MC_antenna::AntennaIDType Id,double freq,double FWHM);


	static PlanckMapsInfoCollType::const_iterator	GetPlanckData(int freq);

	//
	EnvIDsType				InOutEnvID_;
	std::wstring			HP_Dir_;
	std::wstring			Masks_Dir_;
	std::wstring			Masks_DirIn_;
	std::wstring			DirOut_;
	std::wstring			DirOutMaps_;
	std::wstring			DataBuffer_;
	std::wstring			DirOutCurrMap_;
	std::wstring			MaskRejectName_;
	std::wstring			MaskRemoveName_;
	int						syncID_;
	int						Data2Buffer_;

	Healpix_Map<HEALPIX_ATOM_PREC>	map_;
	Healpix_Map<HEALPIX_ATOM_PREC>	IllPixMask_;
	Healpix_Map<HEALPIX_ATOM_PREC>	CutOutMask_;


	double					Epoch_;
	double					RmsRejectLevel_;
	double					PixSz_;
	double					MaskRadiusRatio_;
	int						GeomNSide_;
	int						PtchSz_;
	int						BorderSz_;
	Zeus::PatchGeomType		PatchGeom_;
	HEALPIX_ATOM_PREC		MaskEnlagThreshold_;
	double					MaskSmoothFWHM_;
	int						OpMode_;
	PixExtractProcCtorArgsType::MapFileInfoCollType	FileCollection_;
	int								CurrPatchNumber_;
	int								TotalPatches_;
	PixExtractProcCtorArgsType::MapFileInfoCollType::iterator CurMapPiv_;
	StatInfoType					CurrStat_;
	Trafo							*Trafo_;
	Zeus::LArr2D<pointing>			Ptgs_;
	Zeus::LArr2D<double>			Pix_;
	Zeus::PlanckUnitsValuesTranform		UnitsTransf_;
	Zeus::PwSCoreRandGen			*NoiseGen_;
	Zeus::FourierMachine			*FourierMachine_;
	AntennasCollType				AntColl_;
	static PlanckMapsInfoCollType	MapStaticInfo_;
	static PlanckMapsInfoCollType::const_iterator PlanckEnd_;
};


#endif //MCPIXEXTRACTPROC
