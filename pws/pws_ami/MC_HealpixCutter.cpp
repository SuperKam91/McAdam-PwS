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

//-------------------------------
#include <stdlib.h>
#include <stdio.h>
#include <memory>
#include <time.h>
#include "rotmatrix.h"
#include "fitshandle.h"
#include "ZEUS_PhysicsMath.h"
#include "ZEUS_InOutTxtFile.h"
#include "ZEUS_InOutBinFile.h"
#include "MC_Defaults.h"
#include "MC_HealpixCutter.h"
#include "MC_HPixMapPatchCutting.h"

//-------------------------------

int		HealpixCutter::PixAtDist(int HealPixPix,double dist,DirBasePixels dirT) const
{
	int		Org(HealPixPix);
	int		OrgBasePix(GetNestedCurrBasePix(HealPixPix));
	int		HealPixPixOld;
	double	CurrDist;

	fix_arr<int,8>	res;
	do{
		HealPixPixOld = HealPixPix;
		HPixMap_.neighbors(HealPixPix,res);
		HealPixPix = res[dirT];
		if(GetNestedCurrBasePix(HealPixPix) != OrgBasePix)
			return HealPixPixOld;
	}
	while((CurrDist = GetPixDistance(Org,HealPixPix)) < dist);
	if((CurrDist - dist) <= (dist - GetPixDistance(Org,HealPixPixOld)))
		return HealPixPix;
	return HealPixPixOld;
}
//
void	HealpixCutter::do_Initialize(void)
{
/*
	if(!(PatchsGeom_.Header_.NPatchesPerMainPix_))
	{
		SqrtNPatches_		= Zeus::toInt(((MAXBASEPIXEDGE) / ((PatchsGeom_.Header_.PtchSz_ - (PatchsGeom_.Header_.PtchBorder_ << 1)) * PixelSz_))  + 1.0);
		PatchsGeom_.Header_.NPatchesPerMainPix_ = (SqrtNPatches_*SqrtNPatches_);
	}
	else
	{
		SqrtNPatches_ = PatchsGeom_.Header_.NPatchesPerMainPix_;
		PatchsGeom_.Header_.NPatchesPerMainPix_ *= PatchsGeom_.Header_.NPatchesPerMainPix_;
	}
*/
	SqrtNPatches_ = 16;
	CurrBase_	= BP_PON0;
	SetNewBasePixelState();
	CurrCentralPtg_ = HPixMap_.pix2ang(Healpix00Pix_);
}

//
int		HealpixCutter::PixAtBorder(int HealPixPix,DirBasePixels dirT) const
{
	fix_arr<int,8>	res;
	int		BPix(GetNestedCurrBasePix(HealPixPix));
	int		HealPixPixOld;
	
	do{
		HealPixPixOld = HealPixPix;
		HPixMap_.neighbors(HealPixPix,res);
		HealPixPix = res[dirT];
	}
	while(GetNestedCurrBasePix(HealPixPix) == BPix);
	return HealPixPixOld;
}

//
int		HealpixCutter::SetNewBasePixelState(void)
{

	RowDir_		= ((CurrBase_ < BP_POS0)? BPD_RIGHT: BPD_UP);
	ColumnDir_	= ((CurrBase_ < BP_POS0)? BPD_UP: BPD_RIGHT);

	int	BasePix(GetNestedBasePixCorner(static_cast<BasePixels>(CurrBase_),BPC_LEFT));
	PatchLine_			= 0;
	PatchColumn_		= 0;
	RowMetricDist_		= (PANICFACTOR * Dist2Border(BasePix,RowDir_)) / ((double)SqrtNPatches_);
	BaseRowPix_			= PixAtDist(BasePix,RowMetricDist_/2.0,RowDir_);
	ColumnMetricDist_	= (PANICFACTOR * Dist2Border(BaseRowPix_,ColumnDir_)) / ((double)SqrtNPatches_);
	Healpix00Pix_		= PixAtDist(BaseRowPix_,ColumnMetricDist_/2.0,ColumnDir_);
	Spin_				= EvalSpin(Healpix00Pix_,ColumnMetricDist_);
	return Healpix00Pix_;
}
//
int		HealpixCutter::do_SetNewPatchCentralPixel(void)
{
	if(++PatchColumn_ == SqrtNPatches_)
	{
		if(++PatchLine_ == SqrtNPatches_)
		{
			CurrBase_ = static_cast<BasePixels>(CurrBase_ +  1);
			if(CurrBase_ == BP_EMPTY) return (Healpix00Pix_= -1);
			return SetNewBasePixelState();
		}
		PatchColumn_			= 0;
		BaseRowPix_				= PixAtDist(BaseRowPix_,RowMetricDist_,RowDir_);
		ColumnMetricDist_		= (PANICFACTOR * Dist2Border(BaseRowPix_,ColumnDir_)) / ((double)SqrtNPatches_);
		Healpix00Pix_			= PixAtDist(BaseRowPix_,ColumnMetricDist_/2.0,ColumnDir_);
		Spin_					= EvalSpin(Healpix00Pix_,ColumnMetricDist_);
		return Healpix00Pix_;
	}
	
	Healpix00Pix_			= PixAtDist(Healpix00Pix_,ColumnMetricDist_,ColumnDir_);
	Spin_					= EvalSpin(Healpix00Pix_,ColumnMetricDist_);
	CurrCentralPtg_			= HPixMap_.pix2ang(Healpix00Pix_);
	return Healpix00Pix_;
}


//
double	HealpixCutter::EvalSpin(int CCorner,double dist) const
{
	rotmatrix rot;
	DirBasePixels	SpinDirForw;
	DirBasePixels	SpinDirBack;

	if(ColumnDir_ == BPD_UP)
	{SpinDirForw = BPD_UP;SpinDirBack = BPD_DOWN;}
	else
	{SpinDirForw = BPD_RIGHT;SpinDirBack = BPD_LEFT;}

	pointing	pt0(HPixMap_.pix2ang(CCorner));
	pt0.theta	-= PIOVER2;
	rot.Make_CPAC_Euler_Matrix(pt0.phi,pt0.theta,0.0);
	rot.Transpose();
	vec3		BVec(HPixMap_.pix2ang(PixAtDist(CCorner,dist/2.0,SpinDirBack)));
	vec3		FVec(HPixMap_.pix2ang(PixAtDist(CCorner,dist/2.0,SpinDirForw)));
	BVec = rot.Transform(BVec);
	FVec = rot.Transform(FVec);
	FVec = FVec - BVec;
	FVec.x = 0.0;
	FVec.Normalize();
	if(FVec.z >= 0.0) return std::acos(FVec.y);
	return -(std::acos(FVec.y));
}
//
GeneralMapCutter::~GeneralMapCutter(void)
{}
//
void	GeneralMapCutter::PutCurrGeomData2Place(void)
{
	Trafo	CoordTranslator(MAPEPOCH,MAPEPOCH,PtgsCoordSystem_,Galactic);

	Zeus::PatchGeomType::StorageType::value_type	temp;

	temp.BPixel_			= CurrBase_;
	temp.PatchNumber_		= CurrPatchNumber_;
	temp.Spin_				= Spin_;
	temp.InitRejectRatio_	= IniRejectPercent_;
	temp.FinalRejectRatio_  = 0.0;
	temp.PatchValid_		= 0;
	temp.SNR_				= NonBlind_Snr_;
	temp.SrcIndex_			= NonBlindSrcIndex_;
	temp.PredFlux_			= PredictedFlux_;
	temp.PredRadius_		= PredictedRadius_;
	temp.ErrRadius_			= ErrRadius_;
	temp.ErrFlux_			= ErrFlux_;
	temp.ErrPos_			= ErrPos_;
	temp.ErrRadiusHigh_		= ErrRadiusHigh_;
	temp.ErrRadiusLow_		= ErrRadiusLow_;
	temp.ErrFluxHigh_		= ErrFluxHigh_;
	temp.ErrFluxLow_		= ErrFluxLow_;
	temp.FluxBay_			= FluxBay_;
	temp.RadiusBay_			= RadiusBay_;
	temp.X0Y0Ptg_			= CurrCentralPtg_;
	temp.SourcePtg_			= CurrSrcPtg_;
	temp.QAIN_SrcPtg_		= QAIN_SrcPtg_;
	temp.QAIN_CyR500		= QAIN_CyR500;
	temp.QAIN_T500			= QAIN_T500;
	temp.QAIN_detectable	= QAIN_detectable;

	temp.X0Y0_				= Ptgs_(0,0);
	temp.XLY0_				= Ptgs_(0,PatchsGeom_.Header_.PtchSz_-1);
	temp.X0YL_				= Ptgs_(PatchsGeom_.Header_.PtchSz_-1,0);
	temp.XLYL_				= Ptgs_(PatchsGeom_.Header_.PtchSz_-1,PatchsGeom_.Header_.PtchSz_-1);
	
	if(
		((CurrSrcPtg_.theta >= 0.0) && (CurrSrcPtg_.phi >= 0.0)) &&
		((CurrCentralPtg_.theta != CurrSrcPtg_.theta) || (CurrCentralPtg_.phi != CurrSrcPtg_.phi)) &&
		(PatchsGeom_.Header_.NPatchesPerMainPix_ != 0)
		)
	{
		rotmatrix	rot0;
		rotmatrix	rot1;
		pointing	ptg0(CurrCentralPtg_);
		pointing	ptgSrc(CurrSrcPtg_);

		rot1.Make_Axis_Rotation_Transform(vec3(1.0,0.0,0.0),Spin_);
		ptg0.theta -= PIOVER2;
		rot0.Make_CPAC_Euler_Matrix(ptg0.phi,ptg0.theta,0.0);
		rot0	= rot0 * rot1;
		rot0.Transpose();
		ptgSrc.normalize();
		vec3 Vpix(rot0.Transform(static_cast<vec3>(ptgSrc)));
		Vpix *= (1.0/Vpix.x);
		temp.SrcXCoord_  = (PatchsGeom_.Header_.PtchSz_ >> 1) + Zeus::toInt((Vpix.y/PixelSz_) + 0.5);
		temp.SrcYCoord_  = (PatchsGeom_.Header_.PtchSz_ >> 1) + Zeus::toInt((Vpix.z/PixelSz_) + 0.5);
	}
	else
	{
		temp.SrcXCoord_ = - 1;temp.SrcYCoord_ = -1;
	}

	PatchsGeom_.Storage_.push_back(temp);
}
//
void 	GeneralMapCutter::WritePtgData2Output(const pointing& CentralPixel) const
{
	Zeus::LArr2D<double>	ws(PatchsGeom_.Header_.PtchSz_ * PatchsGeom_.Header_.PtchSz_,PatchsGeom_.Header_.PtchSz_);
	std::wstring			tDir(DirOut_);

	if(Objects2buffer_)
	{
#ifdef WIN32
		std::wstring	PtgsExt(L"Ptgs\\");
#else
		std::wstring	PtgsExt(L"Ptgs/");
#endif
//		tDir = Zeus::RemoveInstanceName(DataBuffer_) + PtgsExt;
		tDir = DataBuffer_ + PtgsExt;
		Zeus::CreateDir(tDir);
	}

	//theta first
	std::auto_ptr<Zeus::GenCollWriter<Zeus::LArr2D<double> > > 
		FWriter0(Zeus::GetGenCollFileWriterHandler(Loki::Type2Type<Zeus::LArr2D<double> >(),Objects2buffer_?1000:InOutEnvID_,
		GetCurrPatchPtgsOutputName(CurrPatchNumber_,0),
		tDir,PatchsGeom_.Header_.PtchSz_,PixelSz_,CentralPixel.phi * RAD2DEGREE,(PIOVER2 - CentralPixel.theta) * RAD2DEGREE,
		static_cast<double>(PatchsGeom_.Header_.PtchSz_) * PixelSz_,0.0));
	
	FWriter0->Initialize();
	FWriter0->Write(PutPtgs2Workspace(0,ws));
	FWriter0->Flush();

	//phi second
	std::auto_ptr<Zeus::GenCollWriter<Zeus::LArr2D<double> > > 
		FWriter1(Zeus::GetGenCollFileWriterHandler(Loki::Type2Type<Zeus::LArr2D<double> >(),Objects2buffer_?1000:InOutEnvID_,
		GetCurrPatchPtgsOutputName(CurrPatchNumber_,1),
		tDir,PatchsGeom_.Header_.PtchSz_,PixelSz_,CentralPixel.phi * RAD2DEGREE,(PIOVER2 - CentralPixel.theta) * RAD2DEGREE,
		static_cast<double>(PatchsGeom_.Header_.PtchSz_) * PixelSz_,0.0));

	FWriter1->Initialize();
	FWriter1->Write(PutPtgs2Workspace(1,ws));
	FWriter1->Flush();
}
//
void 	GeneralMapCutter::WriteGeomData2Output(void) const
{

	
	std::auto_ptr<Zeus::GenCollWriter<Zeus::PatchGeomType> > 
		FWriter(Zeus::GetGenCollFileWriterHandler(Loki::Type2Type<Zeus::PatchGeomType>(),InOutEnvID_,DirOut_,std::wstring(PATCH_GEOMPROPSFNAME)));

	if(FWriter->Remove())
	{
		(Zeus::ConManager::Instance())->PrintStr2Console(std::wstring(L"Could not find/delete this object -> ") + FWriter->GetCollID());
	}
	else
	{
		(Zeus::ConManager::Instance())->PrintStr2Console(std::wstring(L"Object successfully deleted -> ") + FWriter->GetCollID());		
	}

	FWriter->Initialize();
	FWriter->Write(PatchsGeom_);
	FWriter->Flush();

#if (defined(HFIDMC) || defined(LFIDPC)) && defined(HFIDMC_EXTOBJECTS)
	if(Objects2buffer_)
	{
#ifdef WIN32
		std::wstring	MasksExt(L"Geom\\");
#else
		std::wstring	MasksExt(L"Geom/");
#endif
		std::wstring	FullPath(DataBuffer_ + MasksExt);

		Zeus::CreateDir(FullPath);
		std::auto_ptr<Zeus::GenCollWriter<Zeus::PatchGeomType> > 
			FWriterExt(Zeus::GetGenCollFileWriterHandler(Loki::Type2Type<Zeus::PatchGeomType>(),1000,FullPath,std::wstring(PATCH_GEOMPROPSFNAME)));

		if(FWriterExt->Remove())
		{
			(Zeus::ConManager::Instance())->PrintStr2Console(std::wstring(L"Could not find/delete this object -> ") + FWriterExt->GetCollID());
		}
		else
		{
			(Zeus::ConManager::Instance())->PrintStr2Console(std::wstring(L"Object successfully deleted -> ") + FWriterExt->GetCollID());		
		}

		FWriterExt->Initialize();
		FWriterExt->Write(PatchsGeom_);
		FWriterExt->Flush();
	}
#endif
}
//
Zeus::LArr2D<double>&	GeneralMapCutter::PutPtgs2Workspace(int ptgAngle,Zeus::LArr2D<double>& ws) const
{
	Zeus::LArr2D<pointing>::const_iterator	pivIn(Ptgs_.begin());
	Zeus::LArr2D<pointing>::const_iterator	const pivEnd(pivIn + Ptgs_.getSz());
	Zeus::LArr2D<double>::iterator	pivOut(ws.begin());
	if(ptgAngle)
	{
		for(;pivIn!=pivEnd;++pivOut,++pivIn)
		{*pivOut = static_cast<double>(pivIn->phi);}
	}
	else{
		for(;pivIn!=pivEnd;++pivOut,++pivIn)
		{*pivOut = static_cast<double>(pivIn->theta);}
	}
	return ws;
}

std::wstring	GeneralMapCutter::GetCurrPatchPtgsOutputName(int PatchNumber,int ptgType) const
{
	std::wstring	ext(ptgType?PHIFILEEXT:THETAFILEEXT);

	return std::wstring(PTGSPREFIX) + Zeus::PutNumber2Txt(PatchsGeom_.Header_.NSide_) + std::wstring(L"_") + Zeus::PutNumber2Txt(PatchNumber) + ext;
}
//
void	GeneralMapCutter::SetMap(void)
{
	if(PatchsGeom_.Header_.NSide_ == -1)
		errInvParam(ERRCOD_MC_INVALIDPARAM,ERRMSG_MC_INVALIDPARAM,L"nside");
	{
		Healpix_Map<HEALPIX_ATOM_PREC>	Dummy(PatchsGeom_.Header_.NSide_,NEST,nside_dummy());
		HPixMap_.swap(Dummy);
	}

	if(PatchsGeom_.Header_.PtchSz_ == -1)
	{
		if(PatchsGeom_.Header_.NSide_  < 2048)
		{
			if(PatchsGeom_.Header_.NSide_  >= 1024) PatchsGeom_.Header_.PtchSz_	 = MC_PATCHSZSMALL;
			else PatchsGeom_.Header_.PtchSz_ = MC_PATCHSZTINY;
		}
		else
		{
			PatchsGeom_.Header_.PtchSz_ = MC_PATCHSZLARGE;
		}
	}

	if(PatchsGeom_.Header_.PtchBorder_ == -1)
	{
		PatchsGeom_.Header_.PtchBorder_	= (PatchsGeom_.Header_.PtchSz_ >> 2);
	}
}
//
void	GeneralMapCutter::ReadMaskMap(const std::wstring& fname,Healpix_Map<HEALPIX_ATOM_PREC>& Mask)
{
	std::wstring				tDir;
	std::wstring				ExcptStr;
	std::wstring				tFName(Zeus::ExtractFileName(fname,tDir));
	int							Context;

	if(Objects2buffer_ && tDir.empty())
	{
#ifdef WIN32
		std::wstring	MasksExt(L"Masks\\");
#else
		std::wstring	MasksExt(L"Masks/");
#endif
		tDir = Zeus::RemoveInstanceName(DataBuffer_) + MasksExt;
		Context = 1000;
	}
	else
	{
		Context = InOutEnvID_;
		if(tDir.empty())  tDir = DirIn_;
	}

	std::auto_ptr<Zeus::GenHealpixReader<HEALPIX_ATOM_PREC> >	
		FReader(Zeus::GetGenHealpixReader(Loki::Type2Type<HEALPIX_ATOM_PREC>(),Context,tDir,tFName));

	try
	{
		FReader->Initialize(0);
		FReader->Read();
		FReader->Release(Mask);
	}
	catch(Zeus::libException& err)
	{
		ExcptStr = err.what_Xmsg();
	}
	catch(...)
	{
		ExcptStr = std::wstring(L"Non PwS exception.");
	}
	if(!(ExcptStr.empty()))
	{
		Healpix_Map<HEALPIX_ATOM_PREC> temp;
		Mask.swap(temp);
		(Zeus::ConManager::Instance())->PrintStr2Console(ExcptStr);
		(Zeus::ConManager::Instance())->PrintStr2Console(std::wstring(L"Cannot find/read file/object -> ") + FReader->GetCollID() + std::wstring(L"\n\n"));
	}
}
//
void	GeneralMapCutter::Initialize(void)
{
	Spin_				= 0.0;

	PatchsGeom_.Storage_.clear();

	SetMap();
	Ptgs_.Make(PatchsGeom_.Header_.PtchSz_*PatchsGeom_.Header_.PtchSz_,PatchsGeom_.Header_.PtchSz_);
	PixelSz_			= std::sqrt(PITIMES4 / HPixMap_.Npix());
	CurrPatchNumber_	= 0;

	Healpix_Map<HEALPIX_ATOM_PREC>	Dummy(MC_MAXNSIDE,RING,nside_dummy());
	pixCovered_.swap(Dummy);
	pixCovered_.fill(0.0);

	do_Initialize(); 
}
//
void	GeneralMapCutter::CentralPtgOnMask(void)
{}
//
void	NonBlindCutter::CentralPtgOnMask(void)
{
	Zeus::PatchGeomType::StorageType::value_type	temp;

	temp.BPixel_			= CurrBase_;
	temp.PatchNumber_		= CurrPatchNumber_;
	temp.Spin_				= -1.0;
	temp.InitRejectRatio_	= -1.0;
	temp.FinalRejectRatio_  = 0.0;
	temp.PatchValid_		= -1;
	temp.SNR_				= -1.0;
	temp.SrcIndex_			= NonBlindSrcIndex_;
	temp.PredFlux_			= -1.0;
	temp.PredRadius_		= -1.0;
	temp.ErrRadius_			= -1.0;
	temp.ErrFlux_			= -1.0;
	temp.ErrPos_			= -1.0;
	temp.X0Y0Ptg_			= CurrCentralPtg_;
	temp.X0Y0_				= CurrCentralPtg_;
	temp.XLY0_				= CurrCentralPtg_;
	temp.X0YL_				= CurrCentralPtg_;
	temp.XLYL_				= CurrCentralPtg_;

	PatchsGeom_.Storage_.push_back(temp);

	wchar_t	buffer[BUFFERMAXCHAR];

	PRINTINTOBUFFERFUNCT (buffer,BUFFERMAXCHAR,MC_NBLDOUTCUTFORMSTR,
		CurrPatchNumber_,
		static_cast<int>(CurrBase_),
		static_cast<double>((PIOVER2 - CurrCentralPtg_.theta)*RAD2DEGREE),
		static_cast<double>(CurrCentralPtg_.phi*RAD2DEGREE),
		static_cast<int>(NonBlindSrcIndex_)
	);

	(Zeus::ConManager::Instance())->PrintStr2Console(buffer);

	++CurrPatchNumber_;

	return;
}
//
void	GeneralMapCutter::MakePointingsFiles(void)
{
	
	(Zeus::ConManager::Instance())->PrintStr2Console(MC_MASKREADINGBADPIX);
	(Zeus::ConManager::Instance())->PrintStr2Console(L"\n\n");

	std::wstring	Dummy0(MC_MASKBADPIX_DEF);

	ReadMaskMap(Zeus::CorrectDir(DirIn_,1) + Dummy0,HPixMask_);
	
	if(HPixMask_.Map().size() && (HPixMask_.Scheme() != RING))
	{
		HPixMask_.swap_scheme();
	}

	do
	{
		if(PtgOutMask())
		{
			do_Projection();
			if(CheckMinimumGoodPix())
			{
				if(PtgsCoordSystem_ != Galactic)
				{TranslatingPtgs();}
				PutCurrGeomData2Place();
				WritePtgData2Output(CurrCentralPtg_);
				CreateObservedPixels();
				ReportCuttingCurrentParams();
				++CurrPatchNumber_;
				continue;
			}
		}
		
		CentralPtgOnMask();
	}
	while(do_SetNewPatchCentralPixel() >= 0);
	
	PatchsGeom_.Header_.NTotalPatches_ = CurrPatchNumber_;
	ShufflePatchNumbers(PatchsGeom_);
	WriteGeomData2Output();
	AdjustRejectionMask();
}
//
void	GeneralMapCutter::AdjustRejectionMask(void)
{
	Healpix_Map<HEALPIX_ATOM_PREC>	rejectMask(MC_MAXNSIDE,RING,nside_dummy());

	const std::wstring	Dummy0(MC_MASKREJNAME_DEF);

	// Read rejection Mask

	ReadMaskMap(Zeus::CorrectDir(DirIn_,1) + Dummy0,rejectMask);
	
	if(!(rejectMask.Map().size()))
		return;

	if(rejectMask.Scheme() != RING)
	{
		rejectMask.swap_scheme();
	}

	// Merge Masks

	const int	tNPix(pixCovered_.Npix());

	for(int i=0;i<tNPix;++i)
	{ 
		if(pixCovered_[i] < 0.5)
		{
			rejectMask[i] = 0.0;
		}
	}	

	// Write Mask

	Zeus::CreatReadableHealpixFile(InOutEnvID_,DirIn_,Dummy0,Galactic,rejectMask);
//
#if (defined(HFIDMC) || defined(LFIDPC)) && defined(HFIDMC_EXTOBJECTS)
	if(Objects2buffer_)
	{
		Zeus::CreatReadableExternalHealpixFile(DataBuffer_,Dummy0,Galactic,rejectMask);
	}
#endif
	Healpix_Map<HEALPIX_ATOM_PREC> Dummy;
	pixCovered_.swap(Dummy);
}
//
void	GeneralMapCutter::CreateObservedPixels(void)
{

	const int	BorderSz(PatchsGeom_.Header_.PtchBorder_ + 1);
	const int	PtchLimit(PatchsGeom_.Header_.PtchSz_ - PatchsGeom_.Header_.PtchBorder_ - 1);
	fix_arr<int,8>	NeighB;

	for(int j=BorderSz;j<PtchLimit;++j)
	{
		for(int i=BorderSz;i<PtchLimit;++i)
		{
			pixCovered_.neighbors(pixCovered_.ang2pix(Ptgs_(j,i)), NeighB);
			for(int k=0;k<8;++k)
			{
				const int p(NeighB[k]);

				if(p > -1) pixCovered_[p] = 1.0;
			}
		}
	}
}
//
void	GeneralMapCutter::ReadPtgsFileIn(int PatchNumber, Zeus::LArr2D<double>& theta,Zeus::LArr2D<double>& phi)
{
	std::wstring			tDir(DirOut_);

	if(Objects2buffer_)
	{
#ifdef WIN32
		std::wstring	PtgsExt(L"Ptgs\\");
#else
		std::wstring	PtgsExt(L"Ptgs/");
#endif
		if(!Zeus::CheckDir(tDir = (DataBuffer_  + PtgsExt)))
		{
			tDir = Zeus::RemoveInstanceName(DataBuffer_) + PtgsExt;
		}
	}

	std::auto_ptr<Zeus::GenCollReader<Zeus::LArr2D<double> > >
		FReader(Zeus::GetWrkSpFileReaderHandler(Loki::Type2Type<Zeus::LArr2D<double> >(),Objects2buffer_?1000:InOutEnvID_,GetCurrPatchPtgsOutputName(PatchNumber,0),
		PatchsGeom_.Header_.PtchSz_,PatchsGeom_.Header_.PtchSz_,tDir));

	FReader->Initialize();
	FReader->Read();
	FReader->Release(theta);

	std::auto_ptr<Zeus::GenCollReader<Zeus::LArr2D<double> > >
		FReader1(Zeus::GetWrkSpFileReaderHandler(Loki::Type2Type<Zeus::LArr2D<double> >(),Objects2buffer_?1000:InOutEnvID_,GetCurrPatchPtgsOutputName(PatchNumber,1),
		PatchsGeom_.Header_.PtchSz_,PatchsGeom_.Header_.PtchSz_,tDir));

	FReader1->Initialize();
	FReader1->Read();
	FReader1->Release(phi);

}
//
void	GeneralMapCutter::TranslatingPtgs(void)
{
	Zeus::LArr2D<pointing>::iterator		piv(Ptgs_.begin());
	Zeus::LArr2D<pointing>::const_iterator	const end(piv + Ptgs_.getSz());

	Trafo	CoordTranslator(MAPEPOCH,MAPEPOCH,PtgsCoordSystem_,Galactic);

	for(;piv != end; ++piv)
	{
		*piv = CoordTranslator(*piv);
	}
}
//
int		NonBlindCutter::CheckMinimumGoodPix(void)
{
	IniRejectPercent_ = 0.0;
	return 1;
}
//
int		GeneralMapCutter::CheckMinimumGoodPix(void)
{
	Zeus::LArr2D<pointing>::const_iterator		piv(Ptgs_.begin());
	Zeus::LArr2D<pointing>::const_iterator	const end(piv + Ptgs_.getSz());
	int treject(0);
	IniRejectPercent_ = 0.0;

	if(!(HPixMask_.Map().size()))
	{
		return 1;
	}

	for(;piv != end; ++piv)
	{
		pointing	temp(*piv);

		if(HPixMask_[HPixMask_.ang2pix(temp)] < 0.5)
		{++treject;}
	}
	IniRejectPercent_ = static_cast<double>(treject) / static_cast<double>(Ptgs_.getSz());
	if(IniRejectPercent_ > PercentReject_)
		return 0;
	return 1;
}
//
void	GeneralMapCutter::do_Projection(void)
{
	rotmatrix	rot0;
	rotmatrix	rot1;
	pointing	ptg0(CurrCentralPtg_);

	rot1.Make_Axis_Rotation_Transform(vec3(1.0,0.0,0.0),Spin_);
	ptg0.theta -= PIOVER2;
	rot0.Make_CPAC_Euler_Matrix(ptg0.phi,ptg0.theta,0.0);
	rot0	= rot0 * rot1;
	
	double	start((static_cast<double>(-PatchsGeom_.Header_.PtchSz_)/2.0)*PixelSz_);
	for(int i=0; i < PatchsGeom_.Header_.PtchSz_;++i)
	{
		for(int j=0; j < PatchsGeom_.Header_.PtchSz_;++j)
		{
			vec3	pnt(1.0,start+j*PixelSz_,start+i*PixelSz_);
			Ptgs_(i,j) = rot0.Transform(pnt);
		}
	}
}
//
void	GeneralMapCutter::ReportCuttingCurrentParams(void) const
{
	wchar_t	buffer[BUFFERMAXCHAR];

	PRINTINTOBUFFERFUNCT (buffer,BUFFERMAXCHAR,MC_CUTTINGFORMSTR,
		CurrPatchNumber_,
		static_cast<int>(CurrBase_),
		static_cast<double>((PIOVER2 - CurrCentralPtg_.theta)*RAD2DEGREE),
		static_cast<double>(CurrCentralPtg_.phi*RAD2DEGREE),
		static_cast<double>(Spin_*RAD2DEGREE),
		static_cast<double>(IniRejectPercent_*100.0)
	);

	(Zeus::ConManager::Instance())->PrintStr2Console(buffer);
}
//
void	ConstGalacticLatCutter::SetRingsProps(void)
{
	const double RingSz(PITIMES2 * std::sin(CurrCoLat_));
	double SzIncrement(RingSz / Zeus::toInt((RingSz/HalfPatchSz_) + 1.5));

	if(RingSz < 1.0e-16)
	{
		LongIncr_	= PITIMES4;
		CurrLong_	= 0.0;
	}
	else
	{
		LongIncr_	= (PITIMES2 / RingSz) * SzIncrement;
		CurrLong_	= LongIncr_ / 2.0;
	}
	NorthSouth_		= 0;
}
//
void	ConstGalacticLatCutter::do_Initialize(void)
{
	//double GalCutTemp(PIOVER2 - (PatchsGeom_.Header_.GalacticCut_ / RAD2DEGREE));
	CurrBase_		= -1;
	HalfPatchSz_	= (PatchsGeom_.Header_.PtchSz_ - (PatchsGeom_.Header_.PtchBorder_ << 1)) * (PixelSz_ / 2.0);
	LatIncr_		= PIOVER2 / static_cast<double>(Zeus::toInt((PIOVER2 / HalfPatchSz_) + 1.5));
	CurrCoLat_		= 0.0;
	SetRingsProps();
	CurrCentralPtg_	 = pointing(CurrCoLat_,CurrLong_);
}
//
int		ConstGalacticLatCutter::do_SetNewPatchCentralPixel(void)
{
	CurrLong_		+= LongIncr_;
	if(CurrLong_ < PITIMES2 )
	{goto terminate_SetNewPatchCentralPixel;}
	if(!NorthSouth_)
	{
		CurrLong_		= ((LongIncr_ > PITIMES2)?0.0:LongIncr_ / 2.0);
		if((PI - (2.0*CurrCoLat_)) < (LatIncr_ / 32.0))
			return -1; // This is the end my friend
		CurrCoLat_		= PI - CurrCoLat_; 
		NorthSouth_		= 1;
		goto terminate_SetNewPatchCentralPixel;
	}

	CurrCoLat_		= (PI - CurrCoLat_) + LatIncr_;
/*
	if(CurrCoLat_ > (PIOVER2 - (PatchsGeom_.Header_.GalacticCut_ / RAD2DEGREE)))
		return -1; // This is the end my friend
*/
	if(CurrCoLat_ > (PIOVER2 + 1.0e-8))
		return -1; // This is the end my friend

	SetRingsProps();
terminate_SetNewPatchCentralPixel:
	CurrCentralPtg_ = pointing(CurrCoLat_,CurrLong_);
	return 1;
}
//
void	ConstGalacticLatCutter::ShufflePatchNumbers(Zeus::PatchGeomType& pg)
{
	Zeus::PatchGeomType::StorageType::iterator			piv(pg.Storage_.begin());
	Zeus::PatchGeomType::StorageType::const_iterator	const end(pg.Storage_.end());

	for(;piv != end;++piv)
	{piv->SrcIndex_ = piv->PatchNumber_;}
#ifndef NOSHUFFLE
	std::vector<int>	t(pg.Storage_.size());

	piv	= pg.Storage_.begin();
	for(int i=0;piv != end;++piv,++i)
	{t[i] = piv->PatchNumber_;}

	srand(unsigned(time(NULL)));
	std::random_shuffle(t.begin(),t.end());

	piv	= pg.Storage_.begin();

	for(int i=0;piv != end;++piv,++i)
	{piv->PatchNumber_ = t[i];}
#endif //NOSHUFFLE
	std::sort(pg.Storage_.begin(),pg.Storage_.end());
}
//
void	NonBlindCutter::ShufflePatchNumbers(Zeus::PatchGeomType& pg)
{
	Zeus::PatchGeomType::StorageType::iterator			piv(pg.Storage_.begin());
	Zeus::PatchGeomType::StorageType::const_iterator	const end(pg.Storage_.end());

	for(;piv != end;++piv)
	{piv->SrcIndex_ = piv->PatchNumber_;}
#ifndef NOSHUFFLE
	if(PatchsGeom_.Header_.Shuffle_)
	{
		std::vector<int>	t(pg.Storage_.size());

		piv	= pg.Storage_.begin();
		for(int i=0;piv != end;++piv,++i)
		{t[i] = piv->PatchNumber_;}

		srand(unsigned(time(NULL)));
		std::random_shuffle(t.begin(),t.end());

		piv	= pg.Storage_.begin();

		for(int i=0;piv != end;++piv,++i)
		{piv->PatchNumber_ = t[i];}
	}
#endif //NOSHUFFLE
	std::sort(pg.Storage_.begin(),pg.Storage_.end());
}
//
void	NonBlindCutter::do_Initialize(void)
{
	if(ReadCatalogueIn()<0)
		errInvParam(ERRCOD_MC_INVALIDPARAM,ERRMSG_MC_INVALIDPARAM,L"The Non-blind catalogue has the wrong format.");
	NormalisePatchCoords(SrsCat_);
	if(SrsCat_.Header_.CoordsType_)
	{NormaliseCoords(SrsCat_);}
	if(SrsCat_.Header_.CoordSystem_ != PtgsCoordSystem_)
	{TranslateCoords(SrsCat_,PtgsCoordSystem_);}

	if(Chi2MaskFile_.size() > 2)
	{
		Zeus::QA_CltLstCatalogueType	QACOlatList;
		Zeus::QA_ProfilesLstType		QAProfList;
//
		if(((PatchsGeom_.Header_.CollListSz_ = ReadQA_ColLst(Chi2MaskFile_,QACOlatList))<0) || (XtendNB_Info(QACOlatList,MAXAGGRRADIUS)<0))
			errInvParam(ERRCOD_MC_INVALIDPARAM,ERRMSG_MC_INVALIDPARAM,L"The QA colated list catalogue has the wrong format or contain non-matched sources.");
//
		if(MaskRejectFile_.size() > 2)
		{
			if(((ReadQA_Profiles(MaskRejectFile_,QAProfList))<0) || (MatchProfileQA(QAProfList)<0))
				errInvParam(ERRCOD_MC_INVALIDPARAM,ERRMSG_MC_INVALIDPARAM,L"The Profiles file has the wrong format or there are non-matched profiles.");
		}
		else
		{SingleProfileQA();}
	}
//
	{
		wchar_t	buffer[BUFFERMAXCHAR];
		PRINTINTOBUFFERFUNCT
			(buffer,BUFFERMAXCHAR,L"Non-blind catalogue has %d lines.\n",
			static_cast<int>(SrsCat_.Storage_.size()));

			(Zeus::ConManager::Instance())->PrintStr2Console(buffer);
	}
//
	if(!(SrsCat_.Storage_.size()))
		errInvParam(ERRCOD_MC_INVALIDPARAM,ERRMSG_MC_INVALIDPARAM,L"The non-blind list of sources is empty. Pipeline must abort.");

	pivSrcCat_ = SrsCat_.Storage_.begin();
	endSrcCat_ = SrsCat_.Storage_.end();

	if(
		(pivSrcCat_->PatchPtg_.phi < 0.0)			||
		(pivSrcCat_->PatchPtg_.theta < 0.0)			||
		(pivSrcCat_->PatchPtg_.phi > PITIMES2)		||
		(pivSrcCat_->PatchPtg_.theta > PI)			||
		(PatchsGeom_.Header_.NPatchesPerMainPix_ == 0)
	)
	{
		CurrCentralPtg_		= pivSrcCat_->ptg_ ;
	}
	else
	{
		CurrCentralPtg_		= pivSrcCat_->PatchPtg_ ;	
	}

	CurrSrcPtg_			= pivSrcCat_->ptg_;
	Spin_				= pivSrcCat_->Spin_;
	NonBlind_Snr_		= pivSrcCat_->SNR_;
	NonBlindSrcIndex_	= pivSrcCat_->Index_;
	PredictedRadius_	= pivSrcCat_->PredRadiusGLRT_;
	PredictedFlux_		= pivSrcCat_->PredFluxGLRT_;
	ErrRadius_			= pivSrcCat_->ErrRadius_;
	ErrFlux_			= pivSrcCat_->ErrFlux_;
	ErrPos_				= pivSrcCat_->ErrPosition_;
	ErrRadiusHigh_		= pivSrcCat_->ErrRadiusHigh_;
	ErrRadiusLow_		= pivSrcCat_->ErrRadiusLow_;
	ErrFluxHigh_		= pivSrcCat_->ErrFluxHigh_;
	ErrFluxLow_			= pivSrcCat_->ErrFluxLow_;
	FluxBay_			= pivSrcCat_->PredFluxBay_;
	RadiusBay_			= pivSrcCat_->PredRadiusBay_;
	QAIN_CyR500			= pivSrcCat_->QAIN_CyR500;
	QAIN_T500			= pivSrcCat_->QAIN_T500;
	QAIN_detectable		= pivSrcCat_->QAIN_detectable;
	QAIN_SrcPtg_		= pivSrcCat_->QAIN_SrcPtg_;
	CurrBase_			= pivSrcCat_->flagged_;

}
//
int		NonBlindCutter::ReadQA_Profiles(const std::wstring& fname,Zeus::QA_ProfilesLstType& cat)
{
	std::wstring				tDir;
	std::wstring				ExcptStr;
	std::wstring				tFName(Zeus::ExtractFileName(fname,tDir));
	int							Context;

	cat.Storage_.clear();

	if(Objects2buffer_ && tDir.empty())
	{
#ifdef WIN32
		std::wstring	CatsExt(L"Cats\\");
#else
		std::wstring	CatsExt(L"Cats/");
#endif
		tDir = DataBuffer_ + CatsExt;
		Context = 1000;
	}
	else
	{
		Context = InOutEnvID_;
		if(tDir.empty())  tDir = DirIn_;
	}

	std::auto_ptr<Zeus::GenCollReader<Zeus::QA_ProfilesLstType> >	
		FReader(Zeus::GetQA_ProfilesLstReaderHandler(Context,tDir,tFName));

	try
	{
		FReader->Initialize();
		FReader->Read();
	}
	catch(Zeus::libException& err)
	{
		ExcptStr = err.what_Xmsg();
	}
	catch(...)
	{
		ExcptStr = std::wstring(L"Non PwS exception.");
	}
	if(!(ExcptStr.empty()))
	{
		(Zeus::ConManager::Instance())->PrintStr2Console(ExcptStr);
		(Zeus::ConManager::Instance())->PrintStr2Console(std::wstring(L"Cannot find/read file/object -> ") + FReader->GetCollID() + std::wstring(L"\n\n"));
		return -1;
	}

	FReader->Release(cat);

	wchar_t	buffer[BUFFERMAXCHAR];
	swprintf(buffer,BUFFERMAXCHAR,L"Read %d lines from catalogue -> %ls\n",(int)cat.Storage_.size(),(FReader->GetCollID()).c_str());
	(Zeus::ConManager::Instance())->PrintStr2Console(buffer);

	return cat.Storage_.size();
}
//
int		NonBlindCutter::ReadQA_ColLst(const std::wstring& fname,Zeus::QA_CltLstCatalogueType& cat)
{
	std::wstring				tDir;
	std::wstring				ExcptStr;
	std::wstring				tFName(Zeus::ExtractFileName(fname,tDir));
	int							Context;

	cat.Storage_.clear();

	if(Objects2buffer_ && tDir.empty())
	{
#ifdef WIN32
		std::wstring	CatsExt(L"Cats\\");
#else
		std::wstring	CatsExt(L"Cats/");
#endif
		tDir = DataBuffer_ + CatsExt;
		Context = 1000;
	}
	else
	{
		Context = InOutEnvID_;
		if(tDir.empty())  tDir = DirIn_;
	}

	std::auto_ptr<Zeus::GenCollReader<Zeus::QA_CltLstCatalogueType> >	
		FReader(Zeus::GetQA_ColLstReaderHandler(Context,tDir,tFName));

	try
	{
		FReader->Initialize();
		FReader->Read();
	}
	catch(Zeus::libException& err)
	{
		ExcptStr = err.what_Xmsg();
	}
	catch(...)
	{
		ExcptStr = std::wstring(L"Non PwS exception.");
	}
	if(!(ExcptStr.empty()))
	{
		(Zeus::ConManager::Instance())->PrintStr2Console(ExcptStr);
		(Zeus::ConManager::Instance())->PrintStr2Console(std::wstring(L"Cannot find/read file/object -> ") + FReader->GetCollID() + std::wstring(L"\n\n"));
		return -1;
	}

	FReader->Release(cat);

	wchar_t	buffer[BUFFERMAXCHAR];
	swprintf(buffer,BUFFERMAXCHAR,L"Read %d lines from catalogue -> %ls\n",(int)cat.Storage_.size(),(FReader->GetCollID()).c_str());
	(Zeus::ConManager::Instance())->PrintStr2Console(buffer);

/*
	(Zeus::ConManager::Instance())->PrintStr2Console(buffer);
	{
		Zeus::QA_CltLstCatalogueType::StorageType::const_iterator	piv(cat.Storage_.begin());
		Zeus::QA_CltLstCatalogueType::StorageType::const_iterator	const end(cat.Storage_.end());
		for(;piv != end;++piv)
		{
			const double LatDeg((PIOVER2-piv->Cat_GLAT)*RAD2DEGREE);
			const double LongDeg((piv->Cat_GLON)*RAD2DEGREE);

			swprintf(buffer,BUFFERMAXCHAR,L"[GLON-> %9.6g, GLAT-> %9.6g, SNR-> %9.6g]\n",
				LongDeg,LatDeg,piv->Cat_SNR);

			(Zeus::ConManager::Instance())->PrintStr2Console(buffer);
		}
	}
*/
	return cat.Storage_.size();
}
//
int		NonBlindCutter::do_SetNewPatchCentralPixel(void)
{
	++pivSrcCat_;
	if(pivSrcCat_ == endSrcCat_)
		return -1; // This is the end my friend

	if(
		(pivSrcCat_->PatchPtg_.phi < 0.0)			||
		(pivSrcCat_->PatchPtg_.theta < 0.0)			||
		(pivSrcCat_->PatchPtg_.phi > PITIMES2)		||
		(pivSrcCat_->PatchPtg_.theta > PI)			||
		(PatchsGeom_.Header_.NPatchesPerMainPix_ == 0)
	)
	{
		CurrCentralPtg_		= pivSrcCat_->ptg_ ;
	}
	else
	{
		CurrCentralPtg_		= pivSrcCat_->PatchPtg_ ;	
	}
	CurrSrcPtg_			= pivSrcCat_->ptg_;
	NonBlind_Snr_		= pivSrcCat_->SNR_;
	NonBlindSrcIndex_	= pivSrcCat_->Index_;
	PredictedRadius_	= pivSrcCat_->PredRadiusGLRT_;
	PredictedFlux_		= pivSrcCat_->PredFluxGLRT_;
	ErrRadius_			= pivSrcCat_->ErrRadius_;
	ErrFlux_			= pivSrcCat_->ErrFlux_;
	ErrPos_				= pivSrcCat_->ErrPosition_;
	ErrRadiusHigh_		= pivSrcCat_->ErrRadiusHigh_;
	ErrRadiusLow_		= pivSrcCat_->ErrRadiusLow_;
	ErrFluxHigh_		= pivSrcCat_->ErrFluxHigh_;
	ErrFluxLow_			= pivSrcCat_->ErrFluxLow_;
	FluxBay_			= pivSrcCat_->PredFluxBay_;
	RadiusBay_			= pivSrcCat_->PredRadiusBay_;
	Spin_				= pivSrcCat_->Spin_;
	QAIN_CyR500			= pivSrcCat_->QAIN_CyR500;
	QAIN_T500			= pivSrcCat_->QAIN_T500;
	QAIN_detectable		= pivSrcCat_->QAIN_detectable;
	QAIN_SrcPtg_		= pivSrcCat_->QAIN_SrcPtg_;
	CurrBase_			= pivSrcCat_->flagged_;

	return 1;
}
//
void	GeneralMapCutter::NormaliseCoords(Zeus::NonBlingCatType& SrsCat)
{
	Zeus::NonBlingCatType::StorageType::iterator			piv(SrsCat.Storage_.begin());
	Zeus::NonBlingCatType::StorageType::const_iterator	const end(SrsCat.Storage_.end());

	for(;piv != end;++piv)
	{
		if((piv->ptg_.phi   < -1.0e10) || (piv->ptg_.theta < -1.0e10))
			continue;
		piv->ptg_.theta		= ((90.0 - piv->ptg_.theta) / RAD2DEGREE); 
		piv->ptg_.phi		/= RAD2DEGREE;	
	}
}
//
void	GeneralMapCutter::NormalisePatchCoords(Zeus::NonBlingCatType& SrsCat)
{
	Zeus::NonBlingCatType::StorageType::iterator			piv(SrsCat.Storage_.begin());
	Zeus::NonBlingCatType::StorageType::const_iterator	const end(SrsCat.Storage_.end());

	for(;piv != end;++piv)
	{
		if((piv->PatchPtg_.phi   < 0.0) || (piv->PatchPtg_.theta < -1.0e10))
			continue;
		piv->PatchPtg_.theta		= ((90.0 - piv->PatchPtg_.theta) / RAD2DEGREE); 
		piv->PatchPtg_.phi		/= RAD2DEGREE;	
	}
}
//
void	GeneralMapCutter::TranslateCoords(Zeus::NonBlingCatType& SrsCat,const coordsys PtgsCoordSystem)
{
	Trafo	CoordTranslator(SrsCat.Header_.Epoch_,MAPEPOCH,SrsCat.Header_.CoordSystem_,PtgsCoordSystem);
	Zeus::NonBlingCatType::StorageType::iterator			piv(SrsCat.Storage_.begin());
	Zeus::NonBlingCatType::StorageType::const_iterator	const end(SrsCat.Storage_.end());

	for(;piv != end;++piv)
	{
		if((piv->ptg_.phi   < -1.0e10) || (piv->ptg_.theta < -1.0e10))
			continue;
		piv->ptg_ = CoordTranslator(piv->ptg_);
	}
}
//
int		GeneralMapCutter::ReadCatalogue2NBlind(const std::wstring& SrsCatFName,Zeus::NonBlingCatType& SourceCat)
{
	Zeus::CatalogueFormatType	tCat;

	if(!ReadCatalogue(SrsCatFName,tCat))
		return -1;

// SZ_params_

	SourceCat.Header_.CoordsType_		= tCat.Header_.CoordsType_;
	SourceCat.Header_.CoordSystem_		= tCat.Header_.CoordSystem_;
	SourceCat.Header_.Epoch_			= tCat.Header_.Epoch_;
	SourceCat.Header_.SZ_params_		= tCat.Header_.SZ_params_;

	Zeus::CatalogueFormatType::StorageType::const_iterator	OrgPiv(tCat.Storage_.begin());
	Zeus::CatalogueFormatType::StorageType::const_iterator	const OrgEnd(tCat.Storage_.end());
	Zeus::NonBlingCatType::StorageType::value_type			t;

	SourceCat.Storage_.clear();

	for(;OrgPiv != OrgEnd; ++OrgPiv)
	{
		t.Index_				= OrgPiv->ID_;
		t.ptg_.phi				= OrgPiv->GalLongDegs_;
		t.ptg_.theta			= OrgPiv->GalLatDegs_;
		t.PatchPtg_.phi			= OrgPiv->PatchGalLongDegs_;
		t.PatchPtg_.theta		= OrgPiv->PatchGalLatDegs_;
		t.ErrPosition_			= OrgPiv->ErrorBars_.TotalPosErrorBar_;
		t.PredFluxGLRT_			= OrgPiv->FluxComptGLRT_;
		t.PredRadiusGLRT_		= OrgPiv->RadiusGLRT_;
//		t.SNR_					= ((OrgPiv->DetectSigma_ < OrgPiv->NormalAmpl_) ? OrgPiv->DetectSigma_ : OrgPiv->NormalAmpl_) ;
		t.SNR_					= OrgPiv->NormalAmpl_;
		t.PredFluxBay_			= OrgPiv->FluxCompt_;
		t.PredRadiusBay_		= OrgPiv->Radius_;
		t.ErrFlux_				= OrgPiv->ErrorBars_.FluxErrorBar_;
		t.ErrRadius_			= OrgPiv->ErrorBars_.RadiusErrorBar_;
		t.ErrRadiusHigh_		= OrgPiv->ErrorBars_.HighRadiusErrorBar_;
		t.ErrRadiusLow_			= OrgPiv->ErrorBars_.LowRadiusErrorBar_;
		t.ErrFluxHigh_			= OrgPiv->ErrorBars_.HighFluxErrorBar_;
		t.ErrFluxLow_			= OrgPiv->ErrorBars_.LowFluxErrorBar_;
		t.PredFluxBay_			= OrgPiv->FluxCompt_;
		t.PredRadiusBay_		= OrgPiv->Radius_;
		t.Spin_					= OrgPiv->PatchSpin_;

		SourceCat.Storage_.push_back(t);
	}

	return SourceCat.Storage_.size();
}
//
int		NonBlindCutter::FindClosestDetection2NB(double AggRadius /*arcmin*/,const Zeus::QA_CltLstCatalogueType::StorageType::value_type& t,Zeus::NonBlingCatType::StorageType::iterator& out)
{
	Zeus::NonBlingCatType::StorageType::value_type temp;

	temp.ptg_.theta = t.IN_theta + AggRadius;
	if(temp.ptg_.theta > PI) temp.ptg_.theta = PI;

	Zeus::NonBlingCatType::StorageType::iterator	uBound(std::upper_bound(SrsCat_.Storage_.begin(),SrsCat_.Storage_.end(),temp,SortNonBlindByColat));

	temp.ptg_.theta = t.IN_theta - AggRadius;
	if(temp.ptg_.theta < 0.0) temp.ptg_.theta = 0.0;

	Zeus::NonBlingCatType::StorageType::iterator	lBound(std::lower_bound(SrsCat_.Storage_.begin(),SrsCat_.Storage_.end(),temp,SortNonBlindByColat));

	double CurrMinDist(PITIMES2);
	double CurrDist;
	for(;lBound != uBound;++lBound)
	{
		if((CurrDist=Zeus::Dist2ptgs(t.IN_theta,t.IN_phi,lBound->ptg_.theta,lBound->ptg_.phi)) > AggRadius)
			continue;
		if(CurrDist < CurrMinDist)
		{
			CurrMinDist=CurrDist;
			out=lBound;
		}
	}
	if(CurrMinDist==PITIMES2)
		return 0;

	return 1;
}
//
int		NonBlindCutter::MatchProfileQA(Zeus::QA_ProfilesLstType& proflst)
{

	Zeus::NonBlingCatType::StorageType::iterator		piv(SrsCat_.Storage_.begin());
	Zeus::NonBlingCatType::StorageType::const_iterator const end(SrsCat_.Storage_.end());

	std::auto_ptr<Zeus::SimpsonIntegrator< double, NonBlindCutter> > yTotEvalIntegrator(new Zeus::SimpsonIntegrator< double, NonBlindCutter>(2, MF_MAXITER, this, &NonBlindCutter::ytotFunct, 5.0e-4, true));

	double Y_sphTot, Y_sphR500;
	double Y_cylTot;

	bool ok;
	for(;piv != end;++piv)
	{
		const int ProfileIndex(piv->QAIN_ProfileIndex_-1);

		if((ProfileIndex<0) || (ProfileIndex>=proflst.Storage_.size()))
			return -1;

//		f_.setbetaAB((3.0 - proflst.Storage_[ProfileIndex].Gamma_) / proflst.Storage_[ProfileIndex].Alpha_, (proflst.Storage_[ProfileIndex].Beta_ - 3.0) / proflst.Storage_[ProfileIndex].Alpha_);
//		tAlpha_ = proflst.Storage_[ProfileIndex].Alpha_;
//		YtotAux_ = std::pow(R500Ratio_*proflst.Storage_[ProfileIndex].C500_, proflst.Storage_[ProfileIndex].Alpha_);
//		Y_cylTot= yTotEvalIntegrator->Integrate(0.0, PIOVER2, ok);

		Y_sphR500 = (f_.betainc((3.0 - proflst.Storage_[ProfileIndex].Gamma_) / proflst.Storage_[ProfileIndex].Alpha_,
			(proflst.Storage_[ProfileIndex].Beta_ - 3.0) / proflst.Storage_[ProfileIndex].Alpha_,
			std::pow(proflst.Storage_[ProfileIndex].C500_, proflst.Storage_[ProfileIndex].Alpha_) / (1.0 + std::pow(proflst.Storage_[ProfileIndex].C500_, proflst.Storage_[ProfileIndex].Alpha_))));

		Y_sphTot = (f_.betainc((3.0 - proflst.Storage_[ProfileIndex].Gamma_) / proflst.Storage_[ProfileIndex].Alpha_,
			(proflst.Storage_[ProfileIndex].Beta_ - 3.0) / proflst.Storage_[ProfileIndex].Alpha_,
			std::pow(R500Ratio_*proflst.Storage_[ProfileIndex].C500_, proflst.Storage_[ProfileIndex].Alpha_) / (1.0 + std::pow(R500Ratio_*proflst.Storage_[ProfileIndex].C500_, proflst.Storage_[ProfileIndex].Alpha_))));

		piv->QAIN_T500 /= proflst.Storage_[ProfileIndex].C500_;
		piv->QAIN_CyR500 *= (Y_sphTot / Y_sphR500);

//		piv->QAIN_CyR500 *= (Y_cylTot / Y_sphR500);
	}
	return SrsCat_.Storage_.size();
}
//
int		NonBlindCutter::SingleProfileQA(void)
{
	Zeus::NonBlingCatType::StorageType::iterator		piv(SrsCat_.Storage_.begin());
	Zeus::NonBlingCatType::StorageType::const_iterator const end(SrsCat_.Storage_.end());
	for(;piv != end;++piv)
	{
		piv->QAIN_T500		/=  SrsCat_.Header_.SZ_params_.MNFW_C500_;
		piv->QAIN_CyR500	*=  SrsCat_.Header_.SZ_params_.MNFW_Ratio_CY500CYR500_;
	}
	return SrsCat_.Storage_.size();
}
//
int		NonBlindCutter::XtendNB_Info(Zeus::QA_CltLstCatalogueType& cat,double AggRadius /*arcmin*/)
{
	{
		Zeus::QA_CltLstCatalogueType::StorageType::iterator			piv(cat.Storage_.begin());
		Zeus::QA_CltLstCatalogueType::StorageType::const_iterator	const end(cat.Storage_.end());
		for(int i=0;piv!=end;++piv,++i)
		{piv->index_ = i;}
	}
	Zeus::QA_CltLstCatalogueType::StorageType::iterator	QA_CltLstNewEnd;
	QA_CltLstNewEnd = std::remove_if(cat.Storage_.begin(),cat.Storage_.end(),CltLst_NotMatchFunctor);
	if(QA_CltLstNewEnd != cat.Storage_.end())
	{cat.Storage_.erase(QA_CltLstNewEnd,cat.Storage_.end());}
	{
		wchar_t	buffer[BUFFERMAXCHAR];
		PRINTINTOBUFFERFUNCT
			(buffer,BUFFERMAXCHAR,L"Collated list has %d non-empty lines.\n",static_cast<int>(cat.Storage_.size()));

		(Zeus::ConManager::Instance())->PrintStr2Console(buffer);

		Zeus::NonBlingCatType::StorageType::iterator piv(SrsCat_.Storage_.begin());
		Zeus::NonBlingCatType::StorageType::const_iterator const end(SrsCat_.Storage_.end());
		for(;piv!=end;++piv)
		{piv->flagged_ = -1;}	
	}

	if(cat.Storage_.size()!=0)
	{
		std::sort(SrsCat_.Storage_.begin(),SrsCat_.Storage_.end(),SortNonBlindByColat);

		Zeus::NonBlingCatType::StorageType::iterator				closest;
		Zeus::QA_CltLstCatalogueType::StorageType::iterator			piv(cat.Storage_.begin());
		Zeus::QA_CltLstCatalogueType::StorageType::const_iterator	const end(cat.Storage_.end());
		for(;piv!=end;++piv)
		{
			if(FindClosestDetection2NB(AggRadius,*piv,closest))
			{
				closest->flagged_			= piv->index_;
				closest->QAIN_CyR500		= piv->IN_CyR500;
				closest->QAIN_T500			= piv->IN_T500;
				closest->QAIN_detectable	= piv->IN_detectable;
				closest->QAIN_ProfileIndex_	= Zeus::toInt(piv->IN_Vrec+0.5);
				closest->QAIN_SrcPtg_		= pointing(piv->IN_theta,piv->IN_phi);
			}
			else
			{
				const double LatDeg((PIOVER2-(*piv).Cat_GLAT)*RAD2DEGREE);
				const double LongDeg((*piv).Cat_GLON*RAD2DEGREE);

				{
					wchar_t	buffer[BUFFERMAXCHAR];
					PRINTINTOBUFFERFUNCT
						(buffer,BUFFERMAXCHAR,L"Warning! I could not find this source [GLON-> %9.6g, GLAT-> %9.6g, SNR-> %9.6g] in my catalogue. This detection won't be processed\n",
						LongDeg,LatDeg,(*piv).Cat_SNR);

						(Zeus::ConManager::Instance())->PrintStr2Console(buffer);
				}
				// if this error is serious enough then return
			}
		}
	}
	Zeus::NonBlingCatType::StorageType::iterator	NonBlindNewEnd;
	NonBlindNewEnd = std::remove_if(SrsCat_.Storage_.begin(),SrsCat_.Storage_.end(),NonBlindNotFinalCatFunctor);
	if(NonBlindNewEnd != SrsCat_.Storage_.end())
	{SrsCat_.Storage_.erase(NonBlindNewEnd,SrsCat_.Storage_.end());}

	std::sort(SrsCat_.Storage_.begin(),SrsCat_.Storage_.end(),SortNonBlindByIndex);
	return SrsCat_.Storage_.size();
}
//
int		GeneralMapCutter::ReadCatalogue(const std::wstring& fname,Zeus::CatalogueFormatType& cat)
{
	std::wstring				tDir;
	std::wstring				ExcptStr;
	std::wstring				tFName(Zeus::ExtractFileName(fname,tDir));
	int							Context;

	cat.Header_ = Zeus::CatalogueFormatType::HeaderType();
	cat.Storage_.clear();

	if(Objects2buffer_ && tDir.empty())
	{
#ifdef WIN32
		std::wstring	CatsExt(L"Cats\\");
#else
		std::wstring	CatsExt(L"Cats/");
#endif
		tDir = DataBuffer_ + CatsExt;
		Context = 1000;
	}
	else
	{
		Context = InOutEnvID_;
		if(tDir.empty())  tDir = DirIn_;
	}

	std::auto_ptr<Zeus::GenCollReader<Zeus::CatalogueFormatType> >	
		FReader(GetGenCollFileReaderHandler(Loki::Type2Type<Zeus::CatalogueFormatType>(),Context,tDir,tFName));

	try
	{
		FReader->Initialize();
		FReader->Read();
	}
	catch(Zeus::libException& err)
	{
		ExcptStr = err.what_Xmsg();
	}
	catch(...)
	{
		ExcptStr = std::wstring(L"Non PwS exception.");
	}
	if(!(ExcptStr.empty()))
	{
		(Zeus::ConManager::Instance())->PrintStr2Console(ExcptStr);
		(Zeus::ConManager::Instance())->PrintStr2Console(std::wstring(L"Cannot find/read file/object -> ") + FReader->GetCollID() + std::wstring(L"\n\n"));
		return 0;
	}

	FReader->Release(cat);

	wchar_t	buffer[BUFFERMAXCHAR];
	swprintf(buffer,BUFFERMAXCHAR,L"Read %d lines from catalogue -> %ls\n",(int)cat.Storage_.size(),(FReader->GetCollID()).c_str());
	(Zeus::ConManager::Instance())->PrintStr2Console(buffer);

	return cat.Storage_.size();
}
//
int		NonBlindCutter::CheckCentralPixel(const pointing&	temp) const
{
	if(
		(temp.phi < 0.0)
		|| (temp.phi >= PITIMES2)
		|| (temp.theta < 0.0)
		|| (temp.theta >= PI)
	)
		return 0;

	return 1;
}
//
void	GeneralMapCutter::CreatHealpixMapWithDetections(const std::wstring& CatFname,const std::wstring& HealpixFname)
{
	Zeus::NonBlingCatType			cat;

#if defined(HFIDMC) || defined(LFIDPC)
	Healpix_Map<HEALPIX_ATOM_PREC>	HPixMap(MC_MAXNSIDE,RING,nside_dummy());
#else
	Healpix_Map<HEALPIX_ATOM_PREC>	HPixMap(MC_DISPLAYNSIDE,RING,nside_dummy());
#endif

	const double PixSz(4.0*(std::sqrt(PITIMES4 / HPixMap.Npix())));

	if(!(ReadCatalogue2NBlind(CatFname,cat)))
		return;

	if(cat.Header_.CoordsType_)
	{NormaliseCoords(cat);}
	if(cat.Header_.CoordSystem_ != PtgsCoordSystem_)
	{TranslateCoords(cat,PtgsCoordSystem_);}

	HPixMap.fill(0.0);

	std::sort(cat.Storage_.begin(),cat.Storage_.end(),SortPSCatBySNR);

	Zeus::NonBlingCatType::StorageType::const_iterator	piv(cat.Storage_.begin());
	Zeus::NonBlingCatType::StorageType::const_iterator	const end(cat.Storage_.end());
	
	std::vector<int>	pixColl;
	for(int i=0;piv != end;++piv)
	{
		if(
			(piv->ptg_.phi   < -1.0e10)		||
			(piv->ptg_.theta < -1.0e10)		||
			(piv->PredRadiusGLRT_ < -1.0e10)
		)
			continue;

		double tRad(5.0 * piv->PredRadiusGLRT_ / RAD2ARCMIN);
		pixColl.clear();
		if(tRad < PixSz) tRad += PixSz;
		HPixMap.query_disc_inclusive(piv->ptg_,tRad,pixColl);
		std::vector<int>::const_iterator	pixPiv(pixColl.begin());
		std::vector<int>::const_iterator	const pixEnd(pixColl.end());

		for(;pixPiv != pixEnd;++pixPiv)
		{
			HPixMap[*pixPiv] = piv->SNR_;
		}
		if(!(i % 80))
#ifdef WIN32
		{wprintf(L"\n");}
		wprintf(L".");
#else
		{printf("\n");}
		printf(".");
#endif
		++i;
	}
#ifdef WIN32
		wprintf(L"\n\n\n");
#else
		printf("\n\n\n");
#endif

	std::wstring			tDir;
	std::wstring			tFName(Zeus::ExtractFileName(HealpixFname,tDir));
	int						Context;


	if(Objects2buffer_ && tDir.empty())
	{
#ifdef WIN32
		std::wstring	MapsExt(L"Maps\\");
#else
		std::wstring	MapsExt(L"Maps/");
#endif
		tDir = DataBuffer_ + MapsExt;
		Zeus::CreateDir(tDir);
		Context = 1000;
	}
	else
	{
		Context = InOutEnvID_;
		if(tDir.empty())  tDir = DirIn_;
	}

	Zeus::CreatReadableHealpixFile(Context,tDir,HealpixFname,Galactic,HPixMap);
}
//
void	GeneralMapCutter::WriteOutCats(const std::wstring& fname,const Zeus::OutputFormatType& cat)
{
	std::wstring				ExcptStr;
	std::wstring				tDir;
	std::wstring				tFName(Zeus::ExtractFileName(fname,tDir));

#if (defined(HFIDMC) || defined(LFIDPC))
	if(tDir.empty())  tDir = DirIn_ + std::wstring(L"_Out");
#else
	if(tDir.empty())  tDir = DirIn_;
#endif

	std::auto_ptr<Zeus::GenCollWriter<Zeus::OutputFormatType> >	
		FWriter(Zeus::GetFormatCatWriterHandler(InOutEnvID_,cat.Header_.PCCHeader_.DetectionType_,tFName,tDir));

	if(FWriter->Remove())
	{
		(Zeus::ConManager::Instance())->PrintStr2Console(std::wstring(L"Could not find/delete this object -> ") + FWriter->GetCollID());
	}
	else
	{
		(Zeus::ConManager::Instance())->PrintStr2Console(std::wstring(L"Object successfully deleted -> ") + FWriter->GetCollID());
	}

	try
	{
		FWriter->Initialize();
		FWriter->Write(cat);
		FWriter->Flush();
	}
	catch(Zeus::libException& err)
	{
		ExcptStr = err.what_Xmsg();
	}
	catch(...)
	{
		ExcptStr = std::wstring(L"Non PwS exception.");
	}
	if(!(ExcptStr.empty()))
	{
		(Zeus::ConManager::Instance())->PrintStr2Console(ExcptStr);
		(Zeus::ConManager::Instance())->PrintStr2Console(std::wstring(L"Cannot write file/object -> ") + FWriter->GetCollID() + std::wstring(L"\n\n"));
		return;
	}
#if (defined(HFIDMC) || defined(LFIDPC)) && defined(HFIDMC_EXTOBJECTS)
	if(Objects2buffer_)
	{
		std::wstring ExcptStr;
#ifdef WIN32
		std::wstring	CatsExt(L"Cats\\");
#else
		std::wstring	CatsExt(L"Cats/");
#endif
		tDir = DataBuffer_ + CatsExt;
		Zeus::CreateDir(tDir);

		std::auto_ptr<Zeus::GenCollWriter<Zeus::OutputFormatType> >	
			FWriter(Zeus::GetFormatCatWriterHandler(3,cat.Header_.PCCHeader_.DetectionType_,tFName,tDir));

		if(FWriter->Remove())
		{
			(Zeus::ConManager::Instance())->PrintStr2Console(std::wstring(L"Could not find/delete this object -> ") + FWriter->GetCollID());
		}
		else
		{
			(Zeus::ConManager::Instance())->PrintStr2Console(std::wstring(L"Object successfully deleted -> ") + FWriter->GetCollID());
		}

		try
		{
			FWriter->Initialize();
			FWriter->Write(cat);
			FWriter->Flush();
		}
		catch(Zeus::libException& err)
		{
			ExcptStr = err.what_Xmsg();
		}
		catch(...)
		{
			ExcptStr = std::wstring(L"Non PwS exception.");
		}
		if(!(ExcptStr.empty()))
		{
			(Zeus::ConManager::Instance())->PrintStr2Console(ExcptStr);
			(Zeus::ConManager::Instance())->PrintStr2Console(std::wstring(L"Cannot write file/object -> ") + FWriter->GetCollID() + std::wstring(L"\n\n"));
			return;
		}
	}
#endif

	wchar_t	buffer[BUFFERMAXCHAR];
	swprintf(buffer,BUFFERMAXCHAR,L"Wrote %d lines into catalogue %ls\n",(int)cat.Storage_.size(),(FWriter->GetCollID()).c_str());
	(Zeus::ConManager::Instance())->PrintStr2Console(buffer);
}
//
void	GeneralMapCutter::DoConvertCatalogue(const Zeus::CatalogueFormatType& catIn,Zeus::OutputFormatType& catOut)
{
// 
	if((catIn.Header_.CollLstSz_>=0) && (catIn.Header_.DetectionType_ == 1000))
	{
		Zeus::CatLineCollType tcatIn(catIn.Storage_);

		Zeus::CatLineCollType::iterator	NewEnd;
		NewEnd = std::remove_if(
			tcatIn.begin(),tcatIn.end(),
			SelectOnCollLst
		);
		if(NewEnd != tcatIn.end())
		{tcatIn.erase(NewEnd,tcatIn.end());}

		std::sort(
			tcatIn.begin(),tcatIn.end(),
			SortColLstByIndex
		);

		Zeus::OutputFormatType::StorageType::value_type	temp;

		catOut.Header_.ExtHeader_ = Zeus::OutputExtensionHeaderType();
		catOut.Header_.PCCHeader_ = catIn.Header_;

		Zeus::CatLineCollType::const_iterator	piv(tcatIn.begin());
		Zeus::CatLineCollType::const_iterator	const end(tcatIn.end());
		int i(0);
		for(;piv!=end;++piv,++i)
		{
			for(;i<piv->CollLstIndex_;++i)
			{catOut.Storage_.push_back(Zeus::OutputFormatType::StorageType::value_type());}
			temp.Cat_ = *piv;
			catOut.Storage_.push_back(temp);
		}
		for(;i<catIn.Header_.CollLstSz_;++i)
		{catOut.Storage_.push_back(Zeus::OutputFormatType::StorageType::value_type());}
	}
	else
	{
		Zeus::OutputFormatType::StorageType::value_type	temp;
		catOut.Header_.ExtHeader_ = Zeus::OutputExtensionHeaderType();
		catOut.Header_.PCCHeader_ = catIn.Header_;

		Zeus::CatalogueFormatType::StorageType::const_iterator	piv(catIn.Storage_.begin());
		Zeus::CatalogueFormatType::StorageType::const_iterator	const end(catIn.Storage_.end());
		for(;piv!=end;++piv)
		{
			temp.Cat_ = *piv;
			catOut.Storage_.push_back(temp);
		}
	}
}
//
void	GeneralMapCutter::ConvertCatalogue(const std::wstring& InCatName,const std::wstring& OutCatName)
{

	Zeus::CatalogueFormatType	InCat;

	if(!ReadCatalogue(InCatName,InCat))
		return;
	
	Zeus::OutputFormatType		DataOut;

	DoConvertCatalogue(InCat,DataOut);

	WriteOutCats(OutCatName,DataOut);
}
//
void	GeneralMapCutter::MakeChi2Catalogue(const std::wstring& InCatCollNames,const std::wstring& OutCatName)
{
	std::vector<Chi2InfoType>			fNames;
	Zeus::CatalogueFormatType			Data;
	Zeus::OutputFormatType				catOut;

	if(ParseChi2InCatNames(InCatCollNames,fNames) && ReadChi2CatIn(fNames,Data))
	{
		// Setup the header for the output catalogue
		DoMakeChi2Catalogue(Data,catOut);
	}

	if(Chi2MaskFile_.size() >= 2)
	{
		ReadMaskMap(Chi2MaskFile_,HPixMask_);
		if(HPixMask_.Map().size()!=0)
			RemoveChi2Mask(catOut);
	}

	RemoveSuspicious(catOut);
	WriteOutCats(OutCatName,catOut);
}
//
void	GeneralMapCutter::RemoveSuspicious(Zeus::OutputFormatType& cat)
{

	Zeus::OutputFormatType::StorageType::iterator	newEnd;
	
	newEnd = std::remove_if(cat.Storage_.begin(),cat.Storage_.end(),RemoveSuspiciousFunct);

	if(newEnd != cat.Storage_.end())
	{cat.Storage_.erase(newEnd,cat.Storage_.end());}
}
//
void	GeneralMapCutter::RemoveChi2Mask(Zeus::OutputFormatType& cat)
{

	Zeus::OutputFormatType::StorageType::iterator	newEnd;
	
	newEnd = std::remove_if(cat.Storage_.begin(),cat.Storage_.end(),Chi2UnderMaskRemFunct(HPixMask_));

	if(newEnd != cat.Storage_.end())
	{cat.Storage_.erase(newEnd,cat.Storage_.end());}
}
//
int		GeneralMapCutter::ParseChi2InCatNames(const std::wstring& InCatCollNames,std::vector<Chi2InfoType>& fNames)
{
	std::vector<std::wstring>	tStrings;
	std::wstring				tStr;

	Zeus::SplitSSString(InCatCollNames,std::wstring(L","),tStrings);

	std::vector<std::wstring>::iterator			piv(tStrings.begin());
	std::vector<std::wstring>::const_iterator	end(tStrings.end());

	if(!tStrings.size())
		return 0;

	for(;piv != end;++piv)
	{
		*piv = Zeus::FullTrim(*piv);
	}
	
	if(tStrings.size() > 1)
	{
		std::sort(tStrings.begin(),tStrings.end());
		tStrings.erase(std::unique(tStrings.begin(),tStrings.end()),tStrings.end());
	}

	piv	= tStrings.begin();
	end = tStrings.end();

	for(;piv != end;++piv)
	{
		fNames.push_back(Chi2InfoType(GetFreqID(*piv),*piv));
	}

	return fNames.size();
}
//
int		GeneralMapCutter::GetFreqID(const std::wstring& fname)
{

	static const wchar_t *	freqsStr[11]	= {L"_SZ_",L"_COMP_",L"_100_",L"_143_",L"_217_",L"_353_",L"_545_",L"_857_",L"_30_",L"_44_",L"_70_"};
	static int				freqs[11]	= {-1,0,1,2,3,4,5,6,7,8,9};

	
	std::wstring	tDir;
	std::wstring	tFName(Zeus::ExtractFileName(fname,tDir));

	for(int i=0;i<11;++i)
	{
		if(tFName.find(freqsStr[i])!= std::wstring::npos)
			return freqs[i];
	}
	return -1000;
}
//
int		GeneralMapCutter::ReadChi2CatIn(const std::vector<Chi2InfoType>& fNames,Zeus::CatalogueFormatType& Data)
{
	Data.Header_ = Zeus::CatalogueFormatType::HeaderType();
	Data.Storage_.clear();

	std::vector<Chi2InfoType>::const_iterator	piv(fNames.begin());
	std::vector<Chi2InfoType>::const_iterator	const end(fNames.end());

	for(;piv != end;++piv)
	{
		if(piv->freqId_ > -2)
		{
			Zeus::CatalogueFormatType tData;
			if(!ReadCatalogue(piv->fileName_,tData))
				continue;
			Data.Header_ = tData.Header_;
			Zeus::CatalogueFormatType::StorageType::iterator		piv2(tData.Storage_.begin());
			Zeus::CatalogueFormatType::StorageType::const_iterator	const end2(tData.Storage_.end());
			for(;piv2 != end2;++piv2)
			{
				piv2->freqId_ = piv->freqId_;
			}
			
			Data.Storage_.insert(Data.Storage_.end(),tData.Storage_.begin(),tData.Storage_.end());
		}
	}

	std::sort(Data.Storage_.begin(),Data.Storage_.end(),SortChi2ByIdFreq);

	return Data.Storage_.size();
}
//
void	GeneralMapCutter::DoMakeChi2Catalogue(Zeus::CatalogueFormatType& Data,Zeus::OutputFormatType& catOut)
{
	catOut.Header_ = Zeus::OutputFormatType::HeaderType();
	catOut.Header_.PCCHeader_ = Data.Header_ ;
	catOut.Storage_.clear();

	if(!Data.Storage_.size())
		return;

	Zeus::CatalogueFormatType::StorageType::const_iterator	piv(Data.Storage_.begin());
	Zeus::CatalogueFormatType::StorageType::const_iterator	pivAux(piv);
	Zeus::CatalogueFormatType::StorageType::const_iterator	const end(Data.Storage_.end());
	int	CurrentID(piv->ID_);

	for(;;)
	{
		{
			Zeus::OutputFormatType::StorageType::value_type temp;
			for(;(pivAux != end) && (CurrentID == pivAux->ID_);++pivAux)
				;
			
			Chi2OneRow(piv,pivAux,temp);

			catOut.Storage_.push_back(temp);
			piv = pivAux;
			if(piv!=end)
			{CurrentID = piv->ID_;}
			else break;
		}
	}
}
//
void	GeneralMapCutter::Chi2OneRow(Zeus::CatalogueFormatType::StorageType::const_iterator& f,
									 Zeus::CatalogueFormatType::StorageType::const_iterator& s,
									 Zeus::OutputFormatType::StorageType::value_type& t
									 )
{
	double	BayFlux;
	double	tFluxGLRT;

	if(f->freqId_ >= 0)
		errInvParam(ERRCOD_MC_INVALIDPARAM,ERRMSG_MC_INVALIDPARAM,L"Chi2, missing main detection field.");
	else
	{
		t.Cat_			= *f;
		BayFlux			= f->FluxCompt_;
		tFluxGLRT		= f->FluxComptGLRT_;
		t.Ext_.SNRF_	= f->NormalAmpl_;
		++f;
	}
// process extra info if available
	if(t.Cat_.FluxComptGLRT_ < -1.0e20)
		return;

	double	Chi2(0.0);
	int		N(-1);
	double	tChi;
	double	*ptrSNR,*ptrFlux;

	for(;f != s;++f)
	{
		if(f->freqId_ == 0)
		{
			tFluxGLRT = f->FluxComptGLRT_;
			SetMainCatFields(f,t.Cat_,BayFlux);
		}
		else
		{
			if((f->SZ_ConversionCte_ > -1.0e20) && (f->FluxComptGLRT_ > -1.0e20) && (f->PatchMF_sigmaSqr_ > -1.0e20))
			{
				tChi	= (tFluxGLRT * f->SZ_ConversionCte_) - f->FluxComptGLRT_;
				tChi	*= (tChi / f->PatchMF_sigmaSqr_);
				Chi2	+= tChi;
				++N;
				// correct for extra information
				GetFreqIdFields(f->freqId_-1,ptrSNR,ptrFlux,t.Ext_);
				*ptrFlux	= f->FluxComptGLRT_;
//				*ptrSNR		= ((f->NormalAmpl_ < f->DetectSigma_) ? f->NormalAmpl_ : f->DetectSigma_);
				*ptrSNR		= f->NormalAmpl_;
			}
		}
	}

	if(N > 0)
	{t.Ext_.CHI2_ = Chi2 / static_cast<double>(N);}
}
//
void	GeneralMapCutter::SetMainCatFields(Zeus::CatalogueFormatType::StorageType::const_iterator& f,Zeus::CatLineType& t, double BayFLux)
{
	if(f->FluxComptGLRT_  > -1.0e20)
	{
		const double	Adjust(BayFLux / f->FluxComptGLRT_);

// Use theoretical sigma
		t.DetectSigma_		= f->NormalAmpl_;
		t.NormalAmpl_		= f->NormalAmpl_;
//		t.NormalAmpl_		= f->NormalAmpl_;
//		t.FluxComptGLRT_	= BayFLux;
		t.FluxComptGLRT_	= f->FluxComptGLRT_;
		t.Gaussianity_		= f->Gaussianity_;
		t.RadiusGLRT_		= f->RadiusGLRT_;
		t.ScaleLikeNoise_	= f->ScaleLikeNoise_;

/*
		Zeus::ScaleLikeNoiseColl::iterator			pivS(t.ScaleLikeNoise_.begin());
		Zeus::ScaleLikeNoiseColl::const_iterator	const endS(t.ScaleLikeNoise_.end());

		for(;pivS != endS; ++pivS)
		{
			pivS->Like_		*= Adjust;
			pivS->Noise_	*= Adjust;
		}
*/
	}
}
//
void	GeneralMapCutter::GetFreqIdFields(int FreqId,double*& ptrSNR,double*& ptrFlux,Zeus::ExtFormatLnType& t)
{
	switch(FreqId)
	{
	case 0:
		ptrFlux	= &(t.Flux5R500_100_);
		ptrSNR	= &(t.SNR_100_);
		break;
	case 1:
		ptrFlux	= &(t.Flux5R500_143_);
		ptrSNR	= &(t.SNR_143_);
		break;
	case 2:
		ptrFlux	= &(t.Flux5R500_217_);
		ptrSNR	= &(t.SNR_217_);
		break;
	case 3:
		ptrFlux	= &(t.Flux5R500_353_);
		ptrSNR	= &(t.SNR_353_);
		break;
	case 4:
		ptrFlux	= &(t.Flux5R500_545_);
		ptrSNR	= &(t.SNR_545_);
		break;
	case 5:
		ptrFlux	= &(t.Flux5R500_857_);
		ptrSNR	= &(t.SNR_857_);
		break;
	case 6:
		ptrFlux	= &(t.Flux5R500_030_);
		ptrSNR	= &(t.SNR_030_);
		break;
	case 7:
		ptrFlux	= &(t.Flux5R500_044_);
		ptrSNR	= &(t.SNR_044_);
		break;
	case 8:	
		ptrFlux	= &(t.Flux5R500_070_);
		ptrSNR	= &(t.SNR_070_);
		break;
	default:
		ptrSNR	= 0;
		ptrFlux = 0;
		break;
	}
}

//
