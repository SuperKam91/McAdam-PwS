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


//------------------------------------
#include "ZEUS_GaussianRandomGen.h"
#include "ZEUS_InOutPipeline.h"
#ifdef PWSMPI
#include <mpi.h>
#endif //PWSMPI

//------------------------------------

namespace	Zeus
{
//
#ifdef	HFIDMC
//
void		CatalogueWriterHFIDMC::CreatStuff(void)
{
	PIOErr			MyErr;
	std::string		tGrpNameStr(Wstr2Str(GrpName_));
	int retries(MAXRETRIES);

_retry:
	if(PIOCheckGroup(const_cast<char *>(tGrpNameStr.c_str())))
	{
		if(MyErr = PIOCreateTAB2DGrp(const_cast<char *>(tGrpNameStr.c_str()),NColumns_)) 
		{
			if(--retries >= 0)
			{
				MySleep(1);
				goto _retry;
			}
			HFIDMC_errIO(ERROR_COD_HFIDMCERRCREATFL,ERROR_MSG_HFIDMCERRCREATFL,GrpName_);
		}
	}

	if (!(TAB2Dgrp_ = PIOOpenTAB2DGrp(const_cast<char *>(tGrpNameStr.c_str()),"w")))
	{HFIDMC_errIO(ERROR_COD_HFIDMCERROPENFL,ERROR_MSG_HFIDMCERROPENFL,GrpName_);}

	if(PIOCheckObject(const_cast<char *>((Wstr2Str(FName_)).c_str()),TAB2Dgrp_))
	{
		if(MyErr = PIOCreateTAB2DObject(const_cast<char *>((Wstr2Str(FName_)).c_str()),"PIODOUBLE",TAB2Dgrp_)) 
		{HFIDMC_errIO_TAB2D(MyErr,ERROR_MSG_HFIDMCCREATEIMGFORMAT,FName_);}
	}

}
//
void		CatWriterOutputHFIDMC::CreatStuff(void)
{
	PIOErr			MyErr;
	std::string		tGrpNameStr(Wstr2Str(GrpName_));
	int retries(MAXRETRIES);

_retry:
	if(PIOCheckGroup(const_cast<char *>(tGrpNameStr.c_str())))
	{
		if(MyErr = PIOCreateTAB2DGrp(const_cast<char *>(tGrpNameStr.c_str()), GetNColumns())) 
		{
			if(--retries >= 0)
			{
				MySleep(1);
				goto _retry;
			}
			HFIDMC_errIO(ERROR_COD_HFIDMCERRCREATFL,ERROR_MSG_HFIDMCERRCREATFL,GrpName_);
		}
	}
	if (!(TAB2Dgrp_ = PIOOpenTAB2DGrp(const_cast<char *>(tGrpNameStr.c_str()),"w")))
	{HFIDMC_errIO(ERROR_COD_HFIDMCERROPENFL,ERROR_MSG_HFIDMCERROPENFL,GrpName_);}
	if(MyErr = PIOCreateTAB2DObject(const_cast<char *>((Wstr2Str(FName_)).c_str()),"PIODOUBLE",TAB2Dgrp_)) 
	{HFIDMC_errIO_TAB2D(MyErr,ERROR_MSG_HFIDMCCREATEIMGFORMAT,FName_);}
}
//
void		WriteObjResultsHFIDMC::CreatStuff(void)
{
	PIOErr			MyErr;
	std::string		tGrpNameStr(Wstr2Str(GrpName_));
	int retries(MAXRETRIES);
	std::auto_ptr<PwSCoreRandGen> NoiseGen_(new Zeus::PwSCoreRandGen(-1));

#ifdef PWSMPI
	int rk;
	MPI_Comm_rank(MPI_COMM_WORLD,&rk);
	if(!rk)
	{
#endif
_retry:
		if(PIOCheckGroup(const_cast<char *>(tGrpNameStr.c_str())))
		{
			MySleep(Zeus::toInt(NoiseGen_->RandDouble()*20.0));
			if(MyErr = PIOCreateTAB2DGrp(const_cast<char *>(tGrpNameStr.c_str()),NColumns_)) 
			{
				if(--retries >= 0)
				{
					MySleep(Zeus::toInt(NoiseGen_->RandDouble()*10.0));
					goto _retry;
				}
				HFIDMC_errIO(ERROR_COD_HFIDMCERRCREATFL,ERROR_MSG_HFIDMCERRCREATFL,GrpName_);
			}
		}
#ifdef PWSMPI
	}
	MPI_Barrier(MPI_COMM_WORLD);
	NoiseGen_->RandInit(rk<<4);
	MySleep(Zeus::toInt(NoiseGen_->RandDouble()*20.0));
#endif
	if (!(TAB2Dgrp_ = PIOOpenTAB2DGrp(const_cast<char *>(tGrpNameStr.c_str()),"w")))
	{HFIDMC_errIO(ERROR_COD_HFIDMCERROPENFL,ERROR_MSG_HFIDMCERROPENFL,GrpName_);}

	if(PIOCheckObject(const_cast<char *>((Wstr2Str(FName_)).c_str()),TAB2Dgrp_))
	{
		if(MyErr = PIOCreateTAB2DObject(const_cast<char *>((Wstr2Str(FName_)).c_str()),"PIODOUBLE",TAB2Dgrp_)) 
		{HFIDMC_errIO_TAB2D(MyErr,ERROR_MSG_HFIDMCCREATEIMGFORMAT,FName_);}
	}

}
//
void		CatalogueWriterHFIDMC::CreateHeader(const CatalogueFormatType & data)
{
	PIOErr			MyErr;
	std::string		tStr(Wstr2Str(FName_));
	int				tCols(static_cast<int>(NColumns_));
	int				tRows(static_cast<int>(data.Storage_.size()));
	
	HFIDMC_WRITEKEYWORDTAB("PIOINT",data.Header_.CoordsType_,HFIDMC_NONBLINDHEADER_DEGREES,"Degrees Rads",const_cast<char *>(tStr.c_str()));

	HFIDMC_WRITEKEYWORDTAB("PIOINT",data.Header_.CoordSystem_,HFIDMC_NONBLINDHEADER_COORDSYS,"Coords Sys",const_cast<char *>(tStr.c_str()));

	HFIDMC_WRITEKEYWORDTAB("PIODOUBLE",data.Header_.Epoch_,HFIDMC_NONBLINDHEADER_EPOCH,"Epoch",const_cast<char *>(tStr.c_str()));

	HFIDMC_WRITEKEYWORDTAB("PIOINT",tCols,HFIDMC_NONBLINDHEADER_CATNCOL,"NColumns",const_cast<char *>(tStr.c_str()));

	HFIDMC_WRITEKEYWORDTAB("PIOINT",tRows,HFIDMC_NONBLINDHEADER_CATNROW,"NRows",const_cast<char *>(tStr.c_str()));

	HFIDMC_WRITEKEYWORDTAB("PIOINT",data.Header_.DetectionType_,HFIDMC_NONBLINDHEADER_DETCTTYPE,"Detection type",const_cast<char *>(tStr.c_str()));

	HFIDMC_WRITEKEYWORDTAB("PIOINT",data.Header_.Estimator_,HFIDMC_NONBLINDHEADER_ESTIMATOR,"Estimator",const_cast<char *>(tStr.c_str()));

	HFIDMC_WRITEKEYWORDTAB("PIOINT",data.Header_.PriorsType_,HFIDMC_NONBLINDHEADER_PRIORTYPE,"Priors type",const_cast<char *>(tStr.c_str()));

	HFIDMC_WRITEKEYWORDTAB("PIOINT",data.Header_.SZ_params_.SZ_Profile_,HFIDMC_NONBLINDHEADER_SZPROF,"SZ profile",const_cast<char *>(tStr.c_str()));

	HFIDMC_WRITEKEYWORDTAB("PIODOUBLE",data.Header_.SZ_params_.VirialRatio_,HFIDMC_NONBLINDHEADER_SZVRATIO,"SZ Virial ratio",const_cast<char *>(tStr.c_str()));

	HFIDMC_WRITEKEYWORDTAB("PIODOUBLE",data.Header_.SZ_params_.MNFW_alpha_,HFIDMC_NONBLINDHEADER_ALPHA,"SZ alpha",const_cast<char *>(tStr.c_str()));

	HFIDMC_WRITEKEYWORDTAB("PIODOUBLE",data.Header_.SZ_params_.MNFW_beta_,HFIDMC_NONBLINDHEADER_BETA,"SZ beta",const_cast<char *>(tStr.c_str()));

	HFIDMC_WRITEKEYWORDTAB("PIODOUBLE",data.Header_.SZ_params_.MNFW_gamma_,HFIDMC_NONBLINDHEADER_GAMMA,"SZ gamma",const_cast<char *>(tStr.c_str()));

	HFIDMC_WRITEKEYWORDTAB("PIODOUBLE",data.Header_.SZ_params_.MNFW_C500_,HFIDMC_NONBLINDHEADER_C500,"SZ C500",const_cast<char *>(tStr.c_str()));

	HFIDMC_WRITEKEYWORDTAB("PIODOUBLE",data.Header_.SZ_params_.MNFW_Ratio_CY500CYR500_,HFIDMC_NONBLINDHEADER_CY500CYR500,"SZ CY500CYR500",const_cast<char *>(tStr.c_str()));

	HFIDMC_WRITEKEYWORDTAB("PIODOUBLE",data.Header_.SZ_params_.FluxCalibCte_,HFIDMC_NONBLINDHEADER_FLUXCAL,"SZ Flux Calib",const_cast<char *>(tStr.c_str()));

	HFIDMC_WRITEKEYWORDTAB("PIODOUBLE",data.Header_.SZ_params_.RadiusCalCte_,HFIDMC_NONBLINDHEADER_RADCAL,"SZ Radius Calib",const_cast<char *>(tStr.c_str()));

	HFIDMC_WRITEKEYWORDTAB("PIOINT",data.Header_.CollLstSz_,HFIDMC_NONBLINDHEADER_COLLSTSZ,"Collated list Sz",const_cast<char *>(tStr.c_str()));
}
//
void		CatWriterOutputHFIDMC::CreateHeaderHelper(const OutputFormatType & data)
{
	PIOErr			MyErr;
	std::string		tStr(Wstr2Str(FName_));
	int				tCols(GetNColumns());
	int				tRows(static_cast<int>(data.Storage_.size()));
	std::string		tStr1(Wstr2Str(data.Header_.ExtHeader_.DataSetName_));
	const char		*ptrStr(tStr1.c_str());

	HFIDMC_WRITEKEYWORDTAB("PIOINT",data.Header_.PCCHeader_.CoordsType_,HFIDMC_NONBLINDHEADER_DEGREES,"Degrees Rads",const_cast<char *>(tStr.c_str()));

	HFIDMC_WRITEKEYWORDTAB("PIOINT",data.Header_.PCCHeader_.CoordSystem_,HFIDMC_NONBLINDHEADER_COORDSYS,"Coords Sys",const_cast<char *>(tStr.c_str()));

	HFIDMC_WRITEKEYWORDTAB("PIODOUBLE",data.Header_.PCCHeader_.Epoch_,HFIDMC_NONBLINDHEADER_EPOCH,"Epoch",const_cast<char *>(tStr.c_str()));

	HFIDMC_WRITEKEYWORDTAB("PIOINT",tCols,HFIDMC_NONBLINDHEADER_CATNCOL,"NColumns",const_cast<char *>(tStr.c_str()));

	HFIDMC_WRITEKEYWORDTAB("PIOINT",tRows,HFIDMC_NONBLINDHEADER_CATNROW,"NRows",const_cast<char *>(tStr.c_str()));

	HFIDMC_WRITEKEYWORDTAB("PIOINT",data.Header_.PCCHeader_.DetectionType_,HFIDMC_NONBLINDHEADER_DETCTTYPE,"Detection type",const_cast<char *>(tStr.c_str()));

	HFIDMC_WRITEKEYWORDTAB("PIOINT",data.Header_.PCCHeader_.Estimator_,HFIDMC_NONBLINDHEADER_ESTIMATOR,"Estimator",const_cast<char *>(tStr.c_str()));

	HFIDMC_WRITEKEYWORDTAB("PIOINT",data.Header_.PCCHeader_.PriorsType_,HFIDMC_NONBLINDHEADER_PRIORTYPE,"Priors type",const_cast<char *>(tStr.c_str()));

	HFIDMC_WRITEKEYWORDTAB("PIOINT",data.Header_.PCCHeader_.SZ_params_.SZ_Profile_,HFIDMC_NONBLINDHEADER_SZPROF,"SZ profile",const_cast<char *>(tStr.c_str()));

	HFIDMC_WRITEKEYWORDTAB("PIODOUBLE",data.Header_.PCCHeader_.SZ_params_.VirialRatio_,HFIDMC_NONBLINDHEADER_SZVRATIO,"SZ Virial ratio",const_cast<char *>(tStr.c_str()));

	HFIDMC_WRITEKEYWORDTAB("PIODOUBLE",data.Header_.PCCHeader_.SZ_params_.MNFW_alpha_,HFIDMC_NONBLINDHEADER_ALPHA,"SZ alpha",const_cast<char *>(tStr.c_str()));

	HFIDMC_WRITEKEYWORDTAB("PIODOUBLE",data.Header_.PCCHeader_.SZ_params_.MNFW_beta_,HFIDMC_NONBLINDHEADER_BETA,"SZ beta",const_cast<char *>(tStr.c_str()));

	HFIDMC_WRITEKEYWORDTAB("PIODOUBLE",data.Header_.PCCHeader_.SZ_params_.MNFW_gamma_,HFIDMC_NONBLINDHEADER_GAMMA,"SZ gamma",const_cast<char *>(tStr.c_str()));

	HFIDMC_WRITEKEYWORDTAB("PIODOUBLE",data.Header_.PCCHeader_.SZ_params_.MNFW_C500_,HFIDMC_NONBLINDHEADER_C500,"SZ C500",const_cast<char *>(tStr.c_str()));

	HFIDMC_WRITEKEYWORDTAB("PIODOUBLE",data.Header_.PCCHeader_.SZ_params_.MNFW_Ratio_CY500CYR500_,HFIDMC_NONBLINDHEADER_CY500CYR500,"SZ CY500CYR500",const_cast<char *>(tStr.c_str()));

	HFIDMC_WRITEKEYWORDTAB("PIODOUBLE",data.Header_.PCCHeader_.SZ_params_.FluxCalibCte_,HFIDMC_NONBLINDHEADER_FLUXCAL,"SZ Flux Calib",const_cast<char *>(tStr.c_str()));

	HFIDMC_WRITEKEYWORDTAB("PIODOUBLE",data.Header_.PCCHeader_.SZ_params_.RadiusCalCte_,HFIDMC_NONBLINDHEADER_RADCAL,"SZ Radius Calib",const_cast<char *>(tStr.c_str()));

	HFIDMC_WRITEKEYWORDTAB("PIOSTRING",ptrStr[0],"DATASETNAME","Data set name",const_cast<char *>(tStr.c_str()));
}
//
int			CatalogueWriterHFIDMC::do_Write(const CatalogueFormatType & dataIN)
{
	CreateHeader(dataIN);

	if(dataIN.Storage_.empty())
		return true;

	PIOLONG					MyErr;
	LArr2D<double>			dataOUT(dataIN.Storage_.size() * NColumns_,NColumns_,SZCATALOGUEDEFAULTVALUE);
	
	Move2Array2D(dataIN,dataOUT);

	PIOSTRING	command;
	
	sprintf(command,"tab=0:%" PRId64 ",0:%" PRId64,(PIOLONG)(NColumns_ - 1),(PIOLONG)(dataIN.Storage_.size() - 1));

	if((MyErr = PIOWriteTAB2DObject(dataOUT.begin(),const_cast<char *>((Wstr2Str(FName_)).c_str()),"PIODOUBLE",command,TAB2Dgrp_)) < 0) 
	{
		HFIDMC_errIO_TAB2D((PIOErr) MyErr,ERROR_MSG_HFIDMCWRITEIMGFORMAT,FName_);
	}

	return true;
}
//
int			CatWriterOutputHFIDMC::do_Write(const OutputFormatType & dataIN)
{
	CreateHeader(dataIN);

	const int NCol(GetNColumns());

	if(dataIN.Storage_.empty())
		return true;

	PIOLONG					MyErr;
	LArr2D<double>			dataOUT(dataIN.Storage_.size() * NCol, NCol,SZCATALOGUEDEFAULTVALUE);
	
	Move2Array2D(dataIN,dataOUT);

	PIOSTRING	command;
	
	sprintf(command,"tab=0:%" PRId64 ",0:%" PRId64,(PIOLONG)(NCol - 1),(PIOLONG)(dataIN.Storage_.size() - 1));

	if((MyErr = PIOWriteTAB2DObject(dataOUT.begin(),const_cast<char *>((Wstr2Str(FName_)).c_str()),"PIODOUBLE",command,TAB2Dgrp_)) < 0) 
	{
		HFIDMC_errIO_TAB2D((PIOErr) MyErr,ERROR_MSG_HFIDMCWRITEIMGFORMAT,FName_);
	}

	return true;
}
//
int			ReadPatchesGeomInfoHFIDMC::ReadHeader(void)
{
	PIOErr			MyErr;
	PIOSTRING		dummyStr;
	PIODOUBLE		dummyF;
	std::string		tStr(Wstr2Str(FName_));

	HFIDMC_READKEYWORDTAB("PIOINT",NRows_,HFIDMC_GEOMHEADERINFO_NROW,const_cast<char *>(tStr.c_str()));
	HFIDMC_READKEYWORDTAB("PIOINT",NColumns_,HFIDMC_GEOMHEADERINFO_NCOL,const_cast<char *>(tStr.c_str()));
	HFIDMC_READKEYWORDTAB("PIOINT",data_.Header_.NSide_,HFIDMC_GEOMHEADERINFO_NSIDE,const_cast<char *>(tStr.c_str()));
	HFIDMC_READKEYWORDTAB("PIOINT",data_.Header_.PtchSz_,HFIDMC_GEOMHEADERINFO_PTCHSZ,const_cast<char *>(tStr.c_str()));
	HFIDMC_READKEYWORDTAB("PIOINT",data_.Header_.PtchBorder_,HFIDMC_GEOMHEADERINFO_PTCHBORDER,const_cast<char *>(tStr.c_str()));
	HFIDMC_READKEYWORDTAB("PIOINT",data_.Header_.NPatchesPerMainPix_,HFIDMC_GEOMHEADERINFO_NPCTMPIX,const_cast<char *>(tStr.c_str()));
	HFIDMC_READKEYWORDTAB("PIOINT",data_.Header_.NTotalPatches_,HFIDMC_GEOMHEADERINFO_NTOTPTCH,const_cast<char *>(tStr.c_str()));
	HFIDMC_READKEYWORDTAB("PIOINT",data_.Header_.CollListSz_,HFIDMC_GEOMHEADERINFO_CLLISTSZ,const_cast<char *>(tStr.c_str()));
	HFIDMC_READKEYWORDTAB("PIODOUBLE",dummyF,HFIDMC_GEOMHEADERINFO_GALCUT,const_cast<char *>(tStr.c_str()));
	data_.Header_.GalacticCut_			= static_cast<double>(dummyF);
	return true;
}
//
void		WritePatchesGeomInfoHFIDMC::CreatHeader(const PatchGeomType & data)
{
	PIOErr			MyErr;
	std::string		tStr(Wstr2Str(FName_));
	PIOINT			NRows(data.Storage_.size());

	HFIDMC_WRITEKEYWORDTAB("PIOINT",NRows,HFIDMC_GEOMHEADERINFO_NROW,"Table No of rows",const_cast<char *>(tStr.c_str()));
	HFIDMC_WRITEKEYWORDTAB("PIOINT",NColumns_,HFIDMC_GEOMHEADERINFO_NCOL,"Table No of columns",const_cast<char *>(tStr.c_str()));
	HFIDMC_WRITEKEYWORDTAB("PIOINT",data.Header_.NSide_,HFIDMC_GEOMHEADERINFO_NSIDE,HFIDMC_GEOHDCOMMENT_NSIDE,const_cast<char *>(tStr.c_str()));
	HFIDMC_WRITEKEYWORDTAB("PIOINT",data.Header_.PtchSz_,HFIDMC_GEOMHEADERINFO_PTCHSZ,HFIDMC_GEOHDCOMMENT_PTCHSZ,const_cast<char *>(tStr.c_str()));
	HFIDMC_WRITEKEYWORDTAB("PIOINT",data.Header_.PtchBorder_,HFIDMC_GEOMHEADERINFO_PTCHBORDER,HFIDMC_GEOHDCOMMENT_PTCHBORDER,const_cast<char *>(tStr.c_str()));
	HFIDMC_WRITEKEYWORDTAB("PIOINT",data.Header_.NPatchesPerMainPix_,HFIDMC_GEOMHEADERINFO_NPCTMPIX,HFIDMC_GEOHDCOMMENT_NPCTMPIX,const_cast<char *>(tStr.c_str()));
	HFIDMC_WRITEKEYWORDTAB("PIOINT",data.Header_.NTotalPatches_,HFIDMC_GEOMHEADERINFO_NTOTPTCH,HFIDMC_GEOHDCOMMENT_NTOTPTCH,const_cast<char *>(tStr.c_str()));
	HFIDMC_WRITEKEYWORDTAB("PIOINT",data.Header_.CollListSz_,HFIDMC_GEOMHEADERINFO_CLLISTSZ,HFIDMC_GEOHDCOMMENT_CLLISTSZ,const_cast<char *>(tStr.c_str()));
	HFIDMC_WRITEKEYWORDTAB("PIODOUBLE",data.Header_.GalacticCut_,HFIDMC_GEOMHEADERINFO_GALCUT,HFIDMC_GEOHDCOMMENT_GALCUT,const_cast<char *>(tStr.c_str()));
}
//
void		WriteObjResultsHFIDMC::CreatHeader(void)
{
	PIOErr			MyErr;
	std::string		tStr(Wstr2Str(FName_));
	PIOINT			NRows(data_.size());
	PIOSTRING		dummyStr;
	PIOINT			DummyInt(0);

/* might be causing an error on Darwin 
	if(PIOReadKeywordObject(((void *) &(DummyInt) ),dummyStr,INTFILE_HEADERINFO_NROW,"PIOINT",const_cast<char *>(tStr.c_str()),TAB2Dgrp_)==0)
	{NRows += DummyInt;}
*/
	HFIDMC_WRITEKEYWORDTAB("PIOINT",NRows,INTFILE_HEADERINFO_NROW,"Table No of rows",const_cast<char *>(tStr.c_str()));

	if(PIOReadKeywordObject(((void *) &(DummyInt) ),dummyStr,INTFILE_HEADERINFO_NCOL,"PIOINT",const_cast<char *>(tStr.c_str()),TAB2Dgrp_))
	{HFIDMC_WRITEKEYWORDTAB("PIOINT",NColumns_,INTFILE_HEADERINFO_NCOL,"Table No of columns",const_cast<char *>(tStr.c_str()));}
}
//
int			ReadObjResultsHFIDMC::ReadHeader(void)
{
	PIOSTRING		dummyStr;
	std::string		tStr(Wstr2Str(FName_));
	PIOINT			DummyInt(0);

	if(PIOReadKeywordObject(((void *) &(DummyInt) ),dummyStr,INTFILE_HEADERINFO_NROW,"PIOINT",const_cast<char *>(tStr.c_str()),TAB2Dgrp_)==0)
	{NRows_ = DummyInt;}

	if(PIOReadKeywordObject(((void *) &(DummyInt) ),dummyStr,INTFILE_HEADERINFO_NCOL,"PIOINT",const_cast<char *>(tStr.c_str()),TAB2Dgrp_)==0)
	{NColumns_ = DummyInt;}
	else
	{NColumns_ = INTERCATNCOLUMNS;}
	
	return true;
}
//
PeakCollReadbleType::HeaderType *ReadObjResultsHFIDMC::do_Initialize(void)
{
	PIOErr			MyErr;
	std::string		tStr(Wstr2Str(FName_));

	TAB2Dgrp_ = PIOOpenTAB2DGrp(const_cast<char *>(Wstr2Str(GrpName_).c_str()),"r");
	if (!TAB2Dgrp_)
	{
		HFIDMC_errIO(ERROR_COD_HFIDMCERROPENFL,ERROR_MSG_HFIDMCERROPENFL,GrpName_);
	}

	ReadHeader();

	PIOLONG FirstLine;
	PIOLONG LastLine;
		
	if(MyErr = PIOGetTAB2DLines(&FirstLine,&LastLine,const_cast<char *>(tStr.c_str()),TAB2Dgrp_))
	{goto error_do_Initialize;}

	if((LastLine < 0) || !NRows_ || !NColumns_)
	{
		NRows_ = 0;
		return &(data_.Header_);
	}

	PIOLONG	TColumns;

	if((TColumns = PIOGetTAB2DColumnGrp(TAB2Dgrp_))<0)
		goto error_do_Initialize;

	if((TColumns !=  NColumns_) || (NRows_ != (LastLine-FirstLine+1)))
	{
		wchar_t buf[BUFFERMAXCHAR];
		
		PRINTINTOBUFFERFUNCT(buf,BUFFERMAXCHAR,L"Warning! Number of columns or lines in table %ls is different !",tStr.c_str());

		(Zeus::ConManager::Instance())->PrintStr2Console(buf);
		NColumns_	=  TColumns;
		NRows_		=  LastLine-FirstLine + 1;
	}	
	return &(data_.Header_);

error_do_Initialize:
HFIDMC_errIO_TAB2D(MyErr,ERROR_MSG_HFIDMCREADIMGFORMAT,const_cast<char *>(tStr.c_str()));
	return 0;
}
//
PatchGeomType::HeaderType *ReadPatchesGeomInfoHFIDMC::do_Initialize(void)
{
	PIOErr			MyErr;
	std::string		tStr(Wstr2Str(FName_));

	TAB2Dgrp_ = PIOOpenTAB2DGrp(const_cast<char *>(Wstr2Str(GrpName_).c_str()),"r");
	if (!TAB2Dgrp_)
	{
		HFIDMC_errIO(ERROR_COD_HFIDMCERROPENFL,ERROR_MSG_HFIDMCERROPENFL,GrpName_);
	}

	ReadHeader();

	PIOLONG FirstLine;
	PIOLONG LastLine;
		
	if(MyErr = PIOGetTAB2DLines(&FirstLine,&LastLine,const_cast<char *>(tStr.c_str()),TAB2Dgrp_))
	{goto error_do_Initialize;}

	if(LastLine <= 0)
		return &(data_.Header_);

	if(FirstLine > 0)
	{goto error_do_Initialize;}

	PIOLONG	TColumns;

	if((TColumns = PIOGetTAB2DColumnGrp(TAB2Dgrp_))<0)
		goto error_do_Initialize;

	if((TColumns !=  NColumns_) || (NRows_ != (LastLine + 1)))
	{
		wchar_t buf[BUFFERMAXCHAR];
		
		PRINTINTOBUFFERFUNCT(buf,BUFFERMAXCHAR,L"Warning! Number of columns or lines in table %ls is different !",tStr.c_str());

		(Zeus::ConManager::Instance())->PrintStr2Console(buf);
		NColumns_	=  TColumns;
		// NRows_		= (LastLine + 1);
	}	
	return &(data_.Header_);

error_do_Initialize:
HFIDMC_errIO_TAB2D(MyErr,ERROR_MSG_HFIDMCREADIMGFORMAT,const_cast<char *>(tStr.c_str()));
	return 0;
}

//
void		WritePatchesGeomInfoHFIDMC::WriteGeomInfoData(const PatchGeomType & data)
{
	CreatHeader(data);

	if(data.Storage_.empty())
		return;

	PIOSTRING				command;
	PIOLONG					MyErr;
	LArr2D<double>			dataOUT(data.Storage_.size() * NColumns_,NColumns_,SZCATALOGUEDEFAULTVALUE);

	double	*pivOUT(dataOUT.begin());
	double	*pivAuxOUT;

	pivOUT	= dataOUT.begin();

	PatchGeomType::StorageType::const_iterator	pivIN(data.Storage_.begin());
	PatchGeomType::StorageType::const_iterator	const endIN(data.Storage_.end());

	for(;pivIN != endIN;++pivIN,pivOUT += NColumns_)
	{
		pivAuxOUT		= pivOUT;

		*pivAuxOUT++	= static_cast<double>(pivIN->PatchNumber_);
		*pivAuxOUT++	= static_cast<double>(pivIN->BPixel_);
		*pivAuxOUT++	= static_cast<double>(pivIN->SrcXCoord_);
		*pivAuxOUT++	= static_cast<double>(pivIN->SrcYCoord_);
		*pivAuxOUT++	= pivIN->Spin_;
		*pivAuxOUT++	= pivIN->X0Y0Ptg_.theta;
		*pivAuxOUT++	= pivIN->X0Y0Ptg_.phi;
		*pivAuxOUT++	= pivIN->X0Y0_.theta;
		*pivAuxOUT++	= pivIN->X0Y0_.phi;
		*pivAuxOUT++	= pivIN->XLY0_.theta;
		*pivAuxOUT++	= pivIN->XLY0_.phi;
		*pivAuxOUT++	= pivIN->X0YL_.theta;
		*pivAuxOUT++	= pivIN->X0YL_.phi;
		*pivAuxOUT++	= pivIN->XLYL_.theta;
		*pivAuxOUT++	= pivIN->XLYL_.phi;
		*pivAuxOUT++	= pivIN->InitRejectRatio_;
		*pivAuxOUT++	= pivIN->FinalRejectRatio_;
		*pivAuxOUT++	= static_cast<double>(pivIN->PatchValid_);
		*pivAuxOUT++	= pivIN->PredFlux_;
		*pivAuxOUT++	= pivIN->PredRadius_;
		*pivAuxOUT++	= pivIN->SNR_;
		*pivAuxOUT++	= pivIN->ErrRadius_;
		*pivAuxOUT++	= pivIN->ErrFlux_;
		*pivAuxOUT++	= pivIN->ErrPos_;
		*pivAuxOUT++	= static_cast<double>(pivIN->SrcIndex_);
		*pivAuxOUT++	= pivIN->ErrRadiusHigh_;
		*pivAuxOUT++	= pivIN->ErrRadiusLow_;
		*pivAuxOUT++	= pivIN->ErrFluxHigh_;
		*pivAuxOUT++	= pivIN->ErrFluxLow_;
		*pivAuxOUT++	= pivIN->FluxBay_;
		*pivAuxOUT++	= pivIN->RadiusBay_;
		*pivAuxOUT++	= pivIN->SourcePtg_.theta;
		*pivAuxOUT++	= pivIN->SourcePtg_.phi;
		*pivAuxOUT++	= pivIN->QAIN_SrcPtg_.theta;
		*pivAuxOUT++	= pivIN->QAIN_SrcPtg_.phi;
		*pivAuxOUT++	= pivIN->QAIN_CyR500;
		*pivAuxOUT++	= pivIN->QAIN_T500;
		*pivAuxOUT++	= pivIN->QAIN_detectable;
	}

	sprintf(command,"tab=0:%" PRId64 ",0:%" PRId64,(PIOLONG) (NColumns_ - 1),(PIOLONG)(data.Storage_.size() - 1));

	if((MyErr = PIOWriteTAB2DObject(dataOUT.begin(),const_cast<char *>((Wstr2Str(FName_)).c_str()),"PIODOUBLE",command,TAB2Dgrp_)) < 0) 
	{
		HFIDMC_errIO_TAB2D((PIOErr) MyErr,ERROR_MSG_HFIDMCWRITEIMGFORMAT,FName_);
	}

	return;
}
//
int			ReadPatchesGeomInfoHFIDMC::ReadBody(void)
{
	PatchGeomLineType	temp;
	PIOErr				MyErr;
	PIOLONG				NData;
	PIOSTRING			command;
	PIODOUBLE			*tData;
	PIODOUBLE			*tDataPiv;
	PIODOUBLE			*tDataPivAux;
	std::string			tObjName(Wstr2Str(FName_));
	
	data_.Storage_.clear();

	sprintf(command,"tab=0:%" PRId64 ",0:%" PRId64,(PIOLONG)(NColumns_ - 1),(PIOLONG)(NRows_ - 1));

	if((NData = PIOReadTAB2DObject((void **)&tData,const_cast<char *>(tObjName.c_str()),"PIODOUBLE",command,TAB2Dgrp_))< 0)
	{HFIDMC_errIO_TAB2D((PIOErr) NData,ERROR_MSG_HFIDMCREADIMGFORMAT,const_cast<char *>(tObjName.c_str()));}


	if(NData != (NRows_ * NColumns_))
	{
		if((MyErr = PIODeleteTAB2DTable(tData,TAB2Dgrp_))!=0)
			HFIDMC_ReportMemoryLeak(MyErr,ERROR_MSG_HFIDMCDELETEIMGFORMAT,tObjName.c_str());
		HFIDMC_errIO(ERROR_COD_HFIDMCERRIMGWRONGSZ,ERROR_MSG_HFIDMCERRIMGWRONGSZ,FName_);
	}

	tDataPivAux = tData;
	for(int i=0;i<NRows_;++i,tDataPivAux += NColumns_)
	{
		tDataPiv = tDataPivAux;

		temp.PatchNumber_		= toInt(*tDataPiv++);
		temp.BPixel_			= toInt(*tDataPiv++);
		temp.SrcXCoord_			= toInt(*tDataPiv++);
		temp.SrcYCoord_			= toInt(*tDataPiv++);
		temp.Spin_				= *tDataPiv++;
		temp.X0Y0Ptg_.theta		= *tDataPiv++;
		temp.X0Y0Ptg_.phi		= *tDataPiv++;
		temp.X0Y0_.theta		= *tDataPiv++;
		temp.X0Y0_.phi			= *tDataPiv++;
		temp.XLY0_.theta		= *tDataPiv++;
		temp.XLY0_.phi			= *tDataPiv++;
		temp.X0YL_.theta		= *tDataPiv++;
		temp.X0YL_.phi			= *tDataPiv++;
		temp.XLYL_.theta		= *tDataPiv++;
		temp.XLYL_.phi			= *tDataPiv++;
		temp.InitRejectRatio_	= *tDataPiv++;
		temp.FinalRejectRatio_	= *tDataPiv++;
		temp.PatchValid_		=  toInt(*tDataPiv++);
		temp.PredFlux_			= *tDataPiv++;
		temp.PredRadius_		= *tDataPiv++;
		temp.SNR_				= *tDataPiv++;
		temp.ErrRadius_			= *tDataPiv++;
		temp.ErrFlux_			= *tDataPiv++;
		temp.ErrPos_			= *tDataPiv++;
		temp.SrcIndex_			=  toInt(*tDataPiv++);
		temp.ErrRadiusHigh_		= *tDataPiv++;
		temp.ErrRadiusLow_		= *tDataPiv++;
		temp.ErrFluxHigh_		= *tDataPiv++;
		temp.ErrFluxLow_		= *tDataPiv++;
		temp.FluxBay_			= *tDataPiv++;
		temp.RadiusBay_			= *tDataPiv++;
		temp.SourcePtg_.theta	= *tDataPiv++;
		temp.SourcePtg_.phi		= *tDataPiv++;
		temp.QAIN_SrcPtg_.theta = *tDataPiv++;
		temp.QAIN_SrcPtg_.phi	= *tDataPiv++;
		temp.QAIN_CyR500		= *tDataPiv++;
		temp.QAIN_T500			= *tDataPiv++;
		temp.QAIN_detectable	= *tDataPiv++;

		data_.Storage_.push_back(temp);
	}

	if((MyErr = PIODeleteTAB2DTable(tData,TAB2Dgrp_))!=0)
		HFIDMC_ReportMemoryLeak(MyErr,ERROR_MSG_HFIDMCDELETEIMGFORMAT,tObjName.c_str());
	return true;
}
//
int			HFIDMC_ParameterFileReader::DMC_GetValue(const std::wstring& in,DBField& fld)
{
	wchar_t				buffer[BUFFERMAXCHAR];
	long				tInt;
	double				tDouble;

	PRINTINTOBUFFERFUNCT
		(buffer,BUFFERMAXCHAR,L" ; size -> %d ",in.size());
		
	if(in.empty())
	{
		fld = std::wstring(L"");
//		(Zeus::ConManager::Instance())->PrintStr2Console(std::wstring(L"Empty str-> ") + in + std::wstring(buffer));
		return -100;
	}
	if(!Zeus::get_number(in,tInt,CONVDEC))
	{
//		(Zeus::ConManager::Instance())->PrintStr2Console(std::wstring(L"Int -> ") + in  + std::wstring(buffer));
		fld = static_cast<int>(tInt);
		return 1;
	}
	if(!Zeus::get_number(in,tDouble,CONVDEC))
	{
//		(Zeus::ConManager::Instance())->PrintStr2Console(std::wstring(L"Double -> ") + in  + std::wstring(buffer));
		fld = tDouble;
		return 2;
	}
//	(Zeus::ConManager::Instance())->PrintStr2Console(std::wstring(L"Str -> ") + in  + std::wstring(buffer));
	fld = in;
	return 0;
}

//
void		HFIDMC_ParameterFileReader::ReadValuesFromDMC_Pipeline(void)
{
	DBField		tDBField;

	if(DMC_GetValue(Zeus::Achar2Wstr((data_->Header_->GetParContent())->dir_data),tDBField) > -100)
		Insert1element(std::wstring(L"dir_data"),tDBField);

	if(DMC_GetValue(Zeus::Achar2Wstr((data_->Header_->GetParContent())->dir_in_maps),tDBField) > -100)
		Insert1element(std::wstring(L"dir_in_maps"),tDBField);

	if(DMC_GetValue(Zeus::Achar2Wstr((data_->Header_->GetParContent())->dir_pointings),tDBField) > -100)
		Insert1element(std::wstring(L"dir_pointings"),tDBField);

	if(DMC_GetValue(Zeus::Achar2Wstr((data_->Header_->GetParContent())->dir_masks),tDBField) > -100)
		Insert1element(std::wstring(L"dir_masks"),tDBField);

	if(DMC_GetValue(Zeus::Achar2Wstr((data_->Header_->GetParContent())->dir_masks_in),tDBField) > -100)
		Insert1element(std::wstring(L"dir_masks_in"),tDBField);

	if(DMC_GetValue(Zeus::Achar2Wstr((data_->Header_->GetParContent())->dir_out_data),tDBField) > -100)
		Insert1element(std::wstring(L"dir_out_data"),tDBField);

	if(DMC_GetValue(Zeus::Achar2Wstr((data_->Header_->GetParContent())->data_buffer),tDBField) > -100)
		Insert1element(std::wstring(L"data_buffer"),tDBField);

	HFIDMC_INSERTVALUE(GLOBALID_MASKFILENAME);
	HFIDMC_INSERTVALUE(GLOBALID_DATASETNAME);
	HFIDMC_INSERTVALUENOOPT(GLOBALID_NFREQS);
	HFIDMC_INSERTVALUE(GLOBALID_NPRIORPLANES);
	HFIDMC_INSERTVALUE(GLOBALID_NOTALIGNEDOBJS);
	HFIDMC_INSERTVALUE(GLOBALID_PWSSEARCHPOS);
	HFIDMC_INSERTVALUE(GLOBALID_USE2DFORMULA);
	HFIDMC_INSERTVALUE(GLOBALID_NSCALEBINS);
	HFIDMC_INSERTVALUE(GLOBALID_USEBOUNDS);
	HFIDMC_INSERTVALUE(GLOBALID_ASSESS_TYPE);
	HFIDMC_INSERTVALUE(GLOBALID_BOUNDTOL);
	HFIDMC_INSERTVALUE(GLOBALID_SIGMA_THRESHOLD);
	HFIDMC_INSERTVALUE(GLOBALID_SSUB_THRESHOLD);
	HFIDMC_INSERTVALUENOOPT(GLOBALID_SZDETECTION);
	HFIDMC_INSERTVALUE(GLOBALID_SZVIRIALRATIO);
	HFIDMC_INSERTVALUE(GLOBALID_CATALOGUESIGMA);
	HFIDMC_INSERTVALUE(GLOBALID_APODIZEMAPS);
	HFIDMC_INSERTVALUE(GLOBALID_OUTPUTCOORDS);
	HFIDMC_INSERTVALUE(GLOBALID_CATMERGE_AVGT);
	HFIDMC_INSERTVALUE(GLOBALID_GALACTICSIGMA);
	HFIDMC_INSERTVALUE(GLOBALID_OUTUNITS);
	HFIDMC_INSERTVALUE(GLOBALID_FILTERINGSCALES);
	HFIDMC_INSERTVALUE(GLOBALID_SZPROFILE);
	HFIDMC_INSERTVALUE(GLOBALID_CACHESZ);
	HFIDMC_INSERTVALUE(GLOBALID_JFESTIMATTYPE);
	HFIDMC_INSERTVALUE(GLOBALID_OUTPUTRADIUSCAL);
	HFIDMC_INSERTVALUE(GLOBALID_FLUXCALIBCTE);
	HFIDMC_INSERTVALUE(GLOBALID_SRCMAXSCALE);
	HFIDMC_INSERTVALUE(GLOBALID_TWOSTEPSDETECT);
	HFIDMC_INSERTVALUE(GLOBALID_NONBLINDDETECTION);
	HFIDMC_INSERTVALUE(GLOBALID_MELTMAXDIST);
	HFIDMC_INSERTVALUE(GLOBALID_MELTMAXDISTSZLIMIT);
	HFIDMC_INSERTVALUE(GLOBALID_FLUXTHRESHOLD);
	HFIDMC_INSERTVALUE(GLOBALID_OUTPUTLAT);
	HFIDMC_INSERTVALUE(GLOBALID_OUTPUTDEGREES);
	HFIDMC_INSERTVALUE(GLOBALID_GALACTICCUT);
	HFIDMC_INSERTVALUE(GLOBALID_PRIORFLUXMIN);
	HFIDMC_INSERTVALUE(GLOBALID_PRIORFLUXMAX);
	HFIDMC_INSERTVALUE(GLOBALID_PRIORRADDISTMIN);
	HFIDMC_INSERTVALUE(GLOBALID_PRIORRADDISTMAX);
	HFIDMC_INSERTVALUE(GLOBALID_PRIORMASSMIN);
	HFIDMC_INSERTVALUE(GLOBALID_PRIORMASSMAX);
	HFIDMC_INSERTVALUE(GLOBALID_PRIORMAXZ);
	HFIDMC_INSERTVALUE(GLOBALID_PRIORTEMPMIN);
	HFIDMC_INSERTVALUE(GLOBALID_PRIORTEMPMAX);
	HFIDMC_INSERTVALUE(GLOBALID_PRIORGASMASSRATIO);
	HFIDMC_INSERTVALUE(GLOBALID_PRIORSIGMA8);
	HFIDMC_INSERTVALUE(GLOBALID_USEPRESSSCHT);
	HFIDMC_INSERTVALUE(GLOBALID_MN_NLIVEPOINTS);
	HFIDMC_INSERTVALUE(GLOBALID_MN_NINDIVSAMPLES);
	HFIDMC_INSERTVALUE(GLOBALID_MN_FRACTOLEV);
	HFIDMC_INSERTVALUE(GLOBALID_MN_XTRAENLFACT);
	HFIDMC_INSERTVALUE(GLOBALID_PRIORFLUXEXP);
	HFIDMC_INSERTVALUE(GLOBALID_PRIORSRCSCALEMIN);
	HFIDMC_INSERTVALUE(GLOBALID_PRIORSRCSCALEMAX);
	HFIDMC_INSERTVALUE(GLOBALID_PRIORSRCSCALEEXP);
	HFIDMC_INSERTVALUE(GLOBALID_PRIORAVSRCPATCH);
	HFIDMC_INSERTVALUE(GLOBALID_OUTPURITY);
	HFIDMC_INSERTVALUE(GLOBALID_OBSPLANEFREQ0000);
	HFIDMC_INSERTVALUE(GLOBALID_OBSPLANEFREQ0001);
	HFIDMC_INSERTVALUE(GLOBALID_OBSPLANEFREQ0002);
	HFIDMC_INSERTVALUE(GLOBALID_OBSPLANEFREQ0003);
	HFIDMC_INSERTVALUE(GLOBALID_OBSPLANEFREQ0004);
	HFIDMC_INSERTVALUE(GLOBALID_OBSPLANEFREQ0005);
	HFIDMC_INSERTVALUE(GLOBALID_OBSPLANEFREQ0006);
	HFIDMC_INSERTVALUE(GLOBALID_OBSPLANEFREQ0007);
	HFIDMC_INSERTVALUE(GLOBALID_OBSPLANEFREQ0008);
	HFIDMC_INSERTVALUE(GLOBALID_FWHM0000);
	HFIDMC_INSERTVALUE(GLOBALID_FWHM0001);
	HFIDMC_INSERTVALUE(GLOBALID_FWHM0002);
	HFIDMC_INSERTVALUE(GLOBALID_FWHM0003);
	HFIDMC_INSERTVALUE(GLOBALID_FWHM0004);
	HFIDMC_INSERTVALUE(GLOBALID_FWHM0005);
	HFIDMC_INSERTVALUE(GLOBALID_FWHM0006);
	HFIDMC_INSERTVALUE(GLOBALID_FWHM0007);
	HFIDMC_INSERTVALUE(GLOBALID_FWHM0008);

	// MapCut variables
	// Channels' information

	HFIDMC_INSERTVALUENOOPT(MC_FILEMAPNAME0000);
	HFIDMC_INSERTVALUE(MC_FILEFREQ0000);
	HFIDMC_INSERTVALUE(MC_FILEUNITS0000);
	HFIDMC_INSERTVALUE(MC_FILECOORDSYS0000);
	HFIDMC_INSERTVALUE(MC_FILEMAPTYPE0000);
	HFIDMC_INSERTVALUE(MC_FILEUNITCONVFACTOR0000);
	HFIDMC_INSERTVALUE(MC_FILEPROCESSED0000);
	HFIDMC_INSERTVALUENOOPT(MC_FILEFREQID0000);
	HFIDMC_INSERTVALUE(MC_FILEAFWHM0000);
	HFIDMC_INSERTVALUE(MC_FILEPSCAT0000);
	HFIDMC_INSERTVALUE(MC_FILETHSUB0000);
	HFIDMC_INSERTVALUE(MC_FILETHMASK0000);

	HFIDMC_INSERTVALUE(MC_FILEMAPNAME0001);
	HFIDMC_INSERTVALUE(MC_FILEFREQ0001);
	HFIDMC_INSERTVALUE(MC_FILEUNITS0001);
	HFIDMC_INSERTVALUE(MC_FILECOORDSYS0001);
	HFIDMC_INSERTVALUE(MC_FILEMAPTYPE0001);
	HFIDMC_INSERTVALUE(MC_FILEUNITCONVFACTOR0001);
	HFIDMC_INSERTVALUE(MC_FILEPROCESSED0001);
	HFIDMC_INSERTVALUE(MC_FILEFREQID0001);
	HFIDMC_INSERTVALUE(MC_FILEAFWHM0001);
	HFIDMC_INSERTVALUE(MC_FILEPSCAT0001);
	HFIDMC_INSERTVALUE(MC_FILETHSUB0001);
	HFIDMC_INSERTVALUE(MC_FILETHMASK0001);
	
	HFIDMC_INSERTVALUE(MC_FILEMAPNAME0002);
	HFIDMC_INSERTVALUE(MC_FILEFREQ0002);
	HFIDMC_INSERTVALUE(MC_FILEUNITS0002);
	HFIDMC_INSERTVALUE(MC_FILECOORDSYS0002);
	HFIDMC_INSERTVALUE(MC_FILEMAPTYPE0002);
	HFIDMC_INSERTVALUE(MC_FILEUNITCONVFACTOR0002);
	HFIDMC_INSERTVALUE(MC_FILEPROCESSED0002);
	HFIDMC_INSERTVALUE(MC_FILEFREQID0002);
	HFIDMC_INSERTVALUE(MC_FILEAFWHM0002);
	HFIDMC_INSERTVALUE(MC_FILEPSCAT0002);
	HFIDMC_INSERTVALUE(MC_FILETHSUB0002);
	HFIDMC_INSERTVALUE(MC_FILETHMASK0002);
	
	HFIDMC_INSERTVALUE(MC_FILEMAPNAME0003);
	HFIDMC_INSERTVALUE(MC_FILEFREQ0003);
	HFIDMC_INSERTVALUE(MC_FILEUNITS0003);
	HFIDMC_INSERTVALUE(MC_FILECOORDSYS0003);
	HFIDMC_INSERTVALUE(MC_FILEMAPTYPE0003);
	HFIDMC_INSERTVALUE(MC_FILEUNITCONVFACTOR0003);
	HFIDMC_INSERTVALUE(MC_FILEPROCESSED0003);
	HFIDMC_INSERTVALUE(MC_FILEFREQID0003);
	HFIDMC_INSERTVALUE(MC_FILEAFWHM0003);
	HFIDMC_INSERTVALUE(MC_FILEPSCAT0003);
	HFIDMC_INSERTVALUE(MC_FILETHSUB0003);
	HFIDMC_INSERTVALUE(MC_FILETHMASK0003);
	
	HFIDMC_INSERTVALUE(MC_FILEMAPNAME0004);
	HFIDMC_INSERTVALUE(MC_FILEFREQ0004);
	HFIDMC_INSERTVALUE(MC_FILEUNITS0004);
	HFIDMC_INSERTVALUE(MC_FILECOORDSYS0004);
	HFIDMC_INSERTVALUE(MC_FILEMAPTYPE0004);
	HFIDMC_INSERTVALUE(MC_FILEUNITCONVFACTOR0004);
	HFIDMC_INSERTVALUE(MC_FILEPROCESSED0004);
	HFIDMC_INSERTVALUE(MC_FILEFREQID0004);
	HFIDMC_INSERTVALUE(MC_FILEAFWHM0004);
	HFIDMC_INSERTVALUE(MC_FILEPSCAT0004);
	HFIDMC_INSERTVALUE(MC_FILETHSUB0004);
	HFIDMC_INSERTVALUE(MC_FILETHMASK0004);
	
	HFIDMC_INSERTVALUE(MC_FILEMAPNAME0005);
	HFIDMC_INSERTVALUE(MC_FILEFREQ0005);
	HFIDMC_INSERTVALUE(MC_FILEUNITS0005);
	HFIDMC_INSERTVALUE(MC_FILECOORDSYS0005);
	HFIDMC_INSERTVALUE(MC_FILEMAPTYPE0005);
	HFIDMC_INSERTVALUE(MC_FILEUNITCONVFACTOR0005);
	HFIDMC_INSERTVALUE(MC_FILEPROCESSED0005);
	HFIDMC_INSERTVALUE(MC_FILEFREQID0005);
	HFIDMC_INSERTVALUE(MC_FILEAFWHM0005);
	HFIDMC_INSERTVALUE(MC_FILEPSCAT0005);
	HFIDMC_INSERTVALUE(MC_FILETHSUB0005);
	HFIDMC_INSERTVALUE(MC_FILETHMASK0005);
	
	HFIDMC_INSERTVALUE(MC_FILEMAPNAME0006);
	HFIDMC_INSERTVALUE(MC_FILEFREQ0006);
	HFIDMC_INSERTVALUE(MC_FILEUNITS0006);
	HFIDMC_INSERTVALUE(MC_FILECOORDSYS0006);
	HFIDMC_INSERTVALUE(MC_FILEMAPTYPE0006);
	HFIDMC_INSERTVALUE(MC_FILEUNITCONVFACTOR0006);
	HFIDMC_INSERTVALUE(MC_FILEPROCESSED0006);
	HFIDMC_INSERTVALUE(MC_FILEFREQID0006);
	HFIDMC_INSERTVALUE(MC_FILEAFWHM0006);
	HFIDMC_INSERTVALUE(MC_FILEPSCAT0006);
	HFIDMC_INSERTVALUE(MC_FILETHSUB0006);
	HFIDMC_INSERTVALUE(MC_FILETHMASK0006);
	
	HFIDMC_INSERTVALUE(MC_FILEMAPNAME0007);
	HFIDMC_INSERTVALUE(MC_FILEFREQ0007);
	HFIDMC_INSERTVALUE(MC_FILEUNITS0007);
	HFIDMC_INSERTVALUE(MC_FILECOORDSYS0007);
	HFIDMC_INSERTVALUE(MC_FILEMAPTYPE0007);
	HFIDMC_INSERTVALUE(MC_FILEUNITCONVFACTOR0007);
	HFIDMC_INSERTVALUE(MC_FILEPROCESSED0007);
	HFIDMC_INSERTVALUE(MC_FILEFREQID0007);
	HFIDMC_INSERTVALUE(MC_FILEAFWHM0007);
	HFIDMC_INSERTVALUE(MC_FILEPSCAT0007);
	HFIDMC_INSERTVALUE(MC_FILETHSUB0007);
	HFIDMC_INSERTVALUE(MC_FILETHMASK0007);
	
	HFIDMC_INSERTVALUE(MC_FILEMAPNAME0008);
	HFIDMC_INSERTVALUE(MC_FILEFREQ0008);
	HFIDMC_INSERTVALUE(MC_FILEUNITS0008);
	HFIDMC_INSERTVALUE(MC_FILECOORDSYS0008);
	HFIDMC_INSERTVALUE(MC_FILEMAPTYPE0008);
	HFIDMC_INSERTVALUE(MC_FILEUNITCONVFACTOR0008);
	HFIDMC_INSERTVALUE(MC_FILEPROCESSED0008);
	HFIDMC_INSERTVALUE(MC_FILEFREQID0008);
	HFIDMC_INSERTVALUE(MC_FILEAFWHM0008);
	HFIDMC_INSERTVALUE(MC_FILEPSCAT0008);
	HFIDMC_INSERTVALUE(MC_FILETHSUB0008);
	HFIDMC_INSERTVALUE(MC_FILETHMASK0008);

	//	Extra channels

	HFIDMC_INSERTVALUE(MC_FILEMAPNAME0009);
	HFIDMC_INSERTVALUE(MC_FILEFREQ0009);
	HFIDMC_INSERTVALUE(MC_FILEUNITS0009);
	HFIDMC_INSERTVALUE(MC_FILECOORDSYS0009);
	HFIDMC_INSERTVALUE(MC_FILEMAPTYPE0009);
	HFIDMC_INSERTVALUE(MC_FILEUNITCONVFACTOR0009);
	HFIDMC_INSERTVALUE(MC_FILEPROCESSED0009);
	HFIDMC_INSERTVALUE(MC_FILEFREQID0009);
	HFIDMC_INSERTVALUE(MC_FILEAFWHM0009);
	HFIDMC_INSERTVALUE(MC_FILEPSCAT0009);
	HFIDMC_INSERTVALUE(MC_FILETHSUB0009);
	HFIDMC_INSERTVALUE(MC_FILETHMASK0009);

	HFIDMC_INSERTVALUE(MC_FILEMAPNAME0010);
	HFIDMC_INSERTVALUE(MC_FILEFREQ0010);
	HFIDMC_INSERTVALUE(MC_FILEUNITS0010);
	HFIDMC_INSERTVALUE(MC_FILECOORDSYS0010);
	HFIDMC_INSERTVALUE(MC_FILEMAPTYPE0010);
	HFIDMC_INSERTVALUE(MC_FILEUNITCONVFACTOR0010);
	HFIDMC_INSERTVALUE(MC_FILEPROCESSED0010);
	HFIDMC_INSERTVALUE(MC_FILEFREQID0010);
	HFIDMC_INSERTVALUE(MC_FILEAFWHM0010);
	HFIDMC_INSERTVALUE(MC_FILEPSCAT0010);
	HFIDMC_INSERTVALUE(MC_FILETHSUB0010);
	HFIDMC_INSERTVALUE(MC_FILETHMASK0010);
	
	HFIDMC_INSERTVALUE(MC_FILEMAPNAME0011);
	HFIDMC_INSERTVALUE(MC_FILEFREQ0011);
	HFIDMC_INSERTVALUE(MC_FILEUNITS0011);
	HFIDMC_INSERTVALUE(MC_FILECOORDSYS0011);
	HFIDMC_INSERTVALUE(MC_FILEMAPTYPE0011);
	HFIDMC_INSERTVALUE(MC_FILEUNITCONVFACTOR0011);
	HFIDMC_INSERTVALUE(MC_FILEPROCESSED0011);
	HFIDMC_INSERTVALUE(MC_FILEFREQID0011);
	HFIDMC_INSERTVALUE(MC_FILEAFWHM0011);
	HFIDMC_INSERTVALUE(MC_FILEPSCAT0011);
	HFIDMC_INSERTVALUE(MC_FILETHSUB0011);
	HFIDMC_INSERTVALUE(MC_FILETHMASK0011);
	
	HFIDMC_INSERTVALUE(MC_FILEMAPNAME0012);
	HFIDMC_INSERTVALUE(MC_FILEFREQ0012);
	HFIDMC_INSERTVALUE(MC_FILEUNITS0012);
	HFIDMC_INSERTVALUE(MC_FILECOORDSYS0012);
	HFIDMC_INSERTVALUE(MC_FILEMAPTYPE0012);
	HFIDMC_INSERTVALUE(MC_FILEUNITCONVFACTOR0012);
	HFIDMC_INSERTVALUE(MC_FILEPROCESSED0012);
	HFIDMC_INSERTVALUE(MC_FILEFREQID0012);
	HFIDMC_INSERTVALUE(MC_FILEAFWHM0012);
	HFIDMC_INSERTVALUE(MC_FILEPSCAT0012);
	HFIDMC_INSERTVALUE(MC_FILETHSUB0012);
	HFIDMC_INSERTVALUE(MC_FILETHMASK0012);
	
	HFIDMC_INSERTVALUE(MC_FILEMAPNAME0013);
	HFIDMC_INSERTVALUE(MC_FILEFREQ0013);
	HFIDMC_INSERTVALUE(MC_FILEUNITS0013);
	HFIDMC_INSERTVALUE(MC_FILECOORDSYS0013);
	HFIDMC_INSERTVALUE(MC_FILEMAPTYPE0013);
	HFIDMC_INSERTVALUE(MC_FILEUNITCONVFACTOR0013);
	HFIDMC_INSERTVALUE(MC_FILEPROCESSED0013);
	HFIDMC_INSERTVALUE(MC_FILEFREQID0013);
	HFIDMC_INSERTVALUE(MC_FILEAFWHM0013);
	HFIDMC_INSERTVALUE(MC_FILEPSCAT0013);
	HFIDMC_INSERTVALUE(MC_FILETHSUB0013);
	HFIDMC_INSERTVALUE(MC_FILETHMASK0013);
	
	HFIDMC_INSERTVALUE(MC_FILEMAPNAME0014);
	HFIDMC_INSERTVALUE(MC_FILEFREQ0014);
	HFIDMC_INSERTVALUE(MC_FILEUNITS0014);
	HFIDMC_INSERTVALUE(MC_FILECOORDSYS0014);
	HFIDMC_INSERTVALUE(MC_FILEMAPTYPE0014);
	HFIDMC_INSERTVALUE(MC_FILEUNITCONVFACTOR0014);
	HFIDMC_INSERTVALUE(MC_FILEPROCESSED0014);
	HFIDMC_INSERTVALUE(MC_FILEFREQID0014);
	HFIDMC_INSERTVALUE(MC_FILEAFWHM0014);
	HFIDMC_INSERTVALUE(MC_FILEPSCAT0014);
	HFIDMC_INSERTVALUE(MC_FILETHSUB0014);
	HFIDMC_INSERTVALUE(MC_FILETHMASK0014);
	
	HFIDMC_INSERTVALUE(MC_FILEMAPNAME0015);
	HFIDMC_INSERTVALUE(MC_FILEFREQ0015);
	HFIDMC_INSERTVALUE(MC_FILEUNITS0015);
	HFIDMC_INSERTVALUE(MC_FILECOORDSYS0015);
	HFIDMC_INSERTVALUE(MC_FILEMAPTYPE0015);
	HFIDMC_INSERTVALUE(MC_FILEUNITCONVFACTOR0015);
	HFIDMC_INSERTVALUE(MC_FILEPROCESSED0015);
	HFIDMC_INSERTVALUE(MC_FILEFREQID0015);
	HFIDMC_INSERTVALUE(MC_FILEAFWHM0015);
	HFIDMC_INSERTVALUE(MC_FILEPSCAT0015);
	HFIDMC_INSERTVALUE(MC_FILETHSUB0015);
	HFIDMC_INSERTVALUE(MC_FILETHMASK0015);
	
	HFIDMC_INSERTVALUE(MC_FILEMAPNAME0016);
	HFIDMC_INSERTVALUE(MC_FILEFREQ0016);
	HFIDMC_INSERTVALUE(MC_FILEUNITS0016);
	HFIDMC_INSERTVALUE(MC_FILECOORDSYS0016);
	HFIDMC_INSERTVALUE(MC_FILEMAPTYPE0016);
	HFIDMC_INSERTVALUE(MC_FILEUNITCONVFACTOR0016);
	HFIDMC_INSERTVALUE(MC_FILEPROCESSED0016);
	HFIDMC_INSERTVALUE(MC_FILEFREQID0016);
	HFIDMC_INSERTVALUE(MC_FILEAFWHM0016);
	HFIDMC_INSERTVALUE(MC_FILEPSCAT0016);
	HFIDMC_INSERTVALUE(MC_FILETHSUB0016);
	HFIDMC_INSERTVALUE(MC_FILETHMASK0016);
	
	HFIDMC_INSERTVALUE(MC_FILEMAPNAME0017);
	HFIDMC_INSERTVALUE(MC_FILEFREQ0017);
	HFIDMC_INSERTVALUE(MC_FILEUNITS0017);
	HFIDMC_INSERTVALUE(MC_FILECOORDSYS0017);
	HFIDMC_INSERTVALUE(MC_FILEMAPTYPE0017);
	HFIDMC_INSERTVALUE(MC_FILEUNITCONVFACTOR0017);
	HFIDMC_INSERTVALUE(MC_FILEPROCESSED0017);
	HFIDMC_INSERTVALUE(MC_FILEFREQID0017);
	HFIDMC_INSERTVALUE(MC_FILEAFWHM0017);
	HFIDMC_INSERTVALUE(MC_FILEPSCAT0017);
	HFIDMC_INSERTVALUE(MC_FILETHSUB0017);
	HFIDMC_INSERTVALUE(MC_FILETHMASK0017);

	// End of channels information

	HFIDMC_INSERTVALUE(GLOBALID_INTCATNAME);

	HFIDMC_INSERTVALUE(GLOBALID_SZPROFALPHA);
	HFIDMC_INSERTVALUE(GLOBALID_SZPROFBETA);
	HFIDMC_INSERTVALUE(GLOBALID_SZPROFGAMMA);
	HFIDMC_INSERTVALUE(GLOBALID_SZPROF_C500);
	HFIDMC_INSERTVALUE(GLOBALID_SZPROFCY500CYR500);

	HFIDMC_INSERTVALUE(GLOBALID_SZPROFVARALPHAMIN);
	HFIDMC_INSERTVALUE(GLOBALID_SZPROFVARALPHAMAX);
	HFIDMC_INSERTVALUE(GLOBALID_SZPROFVARALPHABIN);
	HFIDMC_INSERTVALUE(GLOBALID_SZPROFVARBETAMIN);
	HFIDMC_INSERTVALUE(GLOBALID_SZPROFVARBETAMAX);
	HFIDMC_INSERTVALUE(GLOBALID_SZPROFVARBETABIN);
	HFIDMC_INSERTVALUE(GLOBALID_SZPROFVARC500MIN);
	HFIDMC_INSERTVALUE(GLOBALID_SZPROFVARC500MAX);
	HFIDMC_INSERTVALUE(GLOBALID_SZPROFVARC500BIN);

	HFIDMC_INSERTVALUENOOPT(MCID_NOFFILES);
	HFIDMC_INSERTVALUENOOPT(MCID_GALACTICCUT);
	HFIDMC_INSERTVALUENOOPT(MCID_MAPCUTTERID);
	HFIDMC_INSERTVALUENOOPT(MCID_COMMANDFLAG);
	HFIDMC_INSERTVALUENOOPT(MCID_PATCHSZ);
	HFIDMC_INSERTVALUE(MCID_PBORDER);
	HFIDMC_INSERTVALUE(MCID_NLINEPTS);
	HFIDMC_INSERTVALUE(MCID_NSIDE);
	HFIDMC_INSERTVALUE(MCID_EPOCH);
	HFIDMC_INSERTVALUE(MCID_MASLENLAGTH);
	HFIDMC_INSERTVALUE(MCID_RMSREJECTLEVEL);
	HFIDMC_INSERTVALUE(MCID_MASKFWHM);
	HFIDMC_INSERTVALUE(MCID_PERCENTREJECT);
	HFIDMC_INSERTVALUE(MCID_PTGSCOORDSYS);
	HFIDMC_INSERTVALUE(MCID_MASKREJECTNAME);
	HFIDMC_INSERTVALUE(MCID_MASKREMOVENAME);
	HFIDMC_INSERTVALUE(MCID_NONBLINDPTGSFILE);

	HFIDMC_INSERTVALUE(GLOBALID_FIRSTPATCH);
	HFIDMC_INSERTVALUE(GLOBALID_LASTPATCH);
	HFIDMC_INSERTVALUE(MCID_SYNCID);

	InsertFinalCatStr();
	InsertOutputCatStr();
}
//
void		HFIDMC_ParameterFileReader::InsertFinalCatStr(void)
{
	std::wstring	Accumulator;

	if((data_->Header_->GetParContent())->flag_final_cat_name)
	{
		Accumulator = Zeus::Achar2Wstr((data_->Header_->GetParContent())->final_cat_name);
	}

	if((data_->Header_->GetParContent())->flag_patch_catalogues)
	{
		for(PIOLONG i=0;i<(data_->Header_->GetParContent())->n_patch_catalogues;++i)
		{
			Accumulator += (Accumulator.empty()?std::wstring():std::wstring(L", ")) + Zeus::Achar2Wstr((data_->Header_->GetParContent())->patch_catalogues[i]);
		}
	}

	DBField		tDBField(Accumulator);
	Insert1element(L"final_cat_name",tDBField);
}
//
void		HFIDMC_ParameterFileReader::InsertOutputCatStr(void)
{
	std::wstring	Accumulator;

	if((data_->Header_->GetParContent())->flag_final_catalogues)
	{
		for(PIOLONG i=0;i<(data_->Header_->GetParContent())->n_final_catalogues;++i)
		{
			Accumulator += (Accumulator.empty()?std::wstring():std::wstring(L", ")) + Zeus::Achar2Wstr((data_->Header_->GetParContent())->final_catalogues[i]);
		}

		DBField		tDBField(Accumulator);
		Insert1element(L"non_blind_ptgs_file",tDBField);
	}
}
//
int			WriteObjResultsHFIDMC::Flush2DB(void)
{

	CreatHeader();

	if(data_.empty())	return 0;

	int					ObIdInitial(ObjectID_);
	PIOSTRING			command;
	PIOLONG				MyErr;
	LArr2D<double>		dataOUT(data_.size() * NColumns_,NColumns_,SZCATALOGUEDEFAULTVALUE);
	std::string			tStr(Wstr2Str(FName_));


	double	*pivOUT(dataOUT.begin());
	double	*pivAuxOUT;

	pivOUT	= dataOUT.begin();

	PeakCollType::const_iterator		pivIN(data_.begin());
	PeakCollType::const_iterator		const endIN(data_.end());
	ScaleLikeNoiseColl::const_iterator	pivS;
	ScaleLikeNoiseColl::const_iterator	endS;


	for(;pivIN != endIN;++pivIN,++ObjectID_,pivOUT += NColumns_)
	{
		pivAuxOUT		= pivOUT;

		*pivAuxOUT++	= static_cast<double>(ObjectID_);
		*pivAuxOUT++	= static_cast<double>(pivIN->PatchNumber_);
		*pivAuxOUT++	= static_cast<double>(pivIN->PK_BayesDetectStat_);
		*pivAuxOUT++	= pivIN->GaussianIndex_;
		*pivAuxOUT++	= pivIN->DetectionSigma_;
		*pivAuxOUT++	= pivIN->SrcAmplNormalised_;
		*pivAuxOUT++	= pivIN->UsedSigmaThreshold_;
		*pivAuxOUT++	= pivIN->JF_UsedSigmaThreshold_;
		*pivAuxOUT++	= pivIN->SrcFlux_;
		*pivAuxOUT++	= pivIN->JF_SrcFlux_.Mode_;
		*pivAuxOUT++	= pivIN->JF_lnRhoTh_;
		*pivAuxOUT++	= pivIN->JF_lnRho_;
		*pivAuxOUT++	= pivIN->JF_lnEvidence_;
		*pivAuxOUT++	= pivIN->JF_lnEvidenceErrBar_;
		*pivAuxOUT++	= pivIN->JF_lnFormFactor_;
		*pivAuxOUT++	= pivIN->JF_lnModelRatio_;
		*pivAuxOUT++	= pivIN->Pos_.YCoord_;
		*pivAuxOUT++	= pivIN->Pos_.XCoord_;
		*pivAuxOUT++	= pivIN->Pos_.JF_YCoord_.Mode_;
		*pivAuxOUT++	= pivIN->Pos_.JF_XCoord_.Mode_;
		*pivAuxOUT++	= pivIN->Pos_.JF_YCoord_.Mean_;
		*pivAuxOUT++	= pivIN->Pos_.JF_XCoord_.Mean_;
		*pivAuxOUT++	= pivIN->GalPt_Colat_;
		*pivAuxOUT++	= pivIN->GalPt_Long_;
		*pivAuxOUT++	= pivIN->SrcFlux_mJys_;
		*pivAuxOUT++	= pivIN->SrcCompt_arcmin2_;
		*pivAuxOUT++	= pivIN->RealParams_.RealScale_;
		*pivAuxOUT++	= pivIN->JF_Radius_.Mode_;
		*pivAuxOUT++	= pivIN->Odds_;
		*pivAuxOUT++	= pivIN->ISNR2_;
		*pivAuxOUT++	= pivIN->CollListIndex_;
		*pivAuxOUT++	= pivIN->ErrorBars_.FluxErrorBar_;
		*pivAuxOUT++	= pivIN->ErrorBars_.TotalPosErrorBar_;
		*pivAuxOUT++	= pivIN->ErrorBars_.RadiusErrorBar_;
		*pivAuxOUT++	= pivIN->y0_;
		*pivAuxOUT++	= pivIN->JF_SrcFlux_.Mean_;
		*pivAuxOUT++	= pivIN->JF_Radius_.Mean_;
		*pivAuxOUT++	= pivIN->CoordPt_Lat_;
		*pivAuxOUT++	= pivIN->CoordPt_Long_;
		*pivAuxOUT++	= pivIN->EquPt_Lat_;
		*pivAuxOUT++	= pivIN->EquPt_Long_;
		*pivAuxOUT++	= pivIN->ErrorBars_.LowFluxErrorBar_;
		*pivAuxOUT++	= pivIN->ErrorBars_.HighFluxErrorBar_;
		*pivAuxOUT++	= pivIN->ErrorBars_.LowRadiusErrorBar_;
		*pivAuxOUT++	= pivIN->ErrorBars_.HighRadiusErrorBar_;
		*pivAuxOUT++	= pivIN->ErrorBars_.DegenCorr_;
		*pivAuxOUT++	= pivIN->ErrorBars_.DegenOrd_;
		*pivAuxOUT++	= pivIN->ErrorBars_.DegenOrdYErr_;
		*pivAuxOUT++	= pivIN->ErrorBars_.DegenSlope_;
		*pivAuxOUT++	= pivIN->ErrorBars_.DegenSlopeYErr_;
//
		if(pivIN->ScaleLikeNoise_.empty())
		{
			*pivAuxOUT++	= -1.0;
			continue;
		}
//		Number of scales is the first value

		*pivAuxOUT++	= static_cast<double>(pivIN->ScaleLikeNoise_.size());

		pivS = pivIN->ScaleLikeNoise_.begin();
		endS = pivIN->ScaleLikeNoise_.end();
		for(;pivS != endS;++pivS)
		{
			*pivAuxOUT++	= pivS->Scale_;
			*pivAuxOUT++	= pivS->Like_;
			*pivAuxOUT++	= pivS->Noise_;		
		}

		if(pivIN->QAResult_[0] < 0.0)
		{
			*pivAuxOUT++	= -1.0;
			continue;
		}

		*pivAuxOUT++	= static_cast<double>(QARESULTSSZ);

		for (int i=0;i<QARESULTSSZ;++i)
		{
			*pivAuxOUT++ = pivIN->QAResult_[i];
		}
	}

	sprintf(command,"tab=0:%" PRId64 ",%" PRId64 ":%" PRId64,(PIOLONG)(NColumns_ - 1),(PIOLONG)(0),(PIOLONG)(data_.size() - 1));

	if((MyErr = PIOWriteTAB2DObject(dataOUT.begin(),const_cast<char *>(tStr.c_str()),"PIODOUBLE",command,TAB2Dgrp_)) < 0) 
	{goto error_WritePeaksColl;}

	return (ObjectID_ - ObIdInitial);

error_WritePeaksColl:
HFIDMC_errIO_TAB2D((PIOErr) MyErr,ERROR_MSG_HFIDMCWRITEIMGFORMAT,FName_);
return (ObjectID_ - ObIdInitial);
}
//
int			ReadObjResultsHFIDMC::ReadBody(void)
{

	PeakCollReadbleType::StorageType::value_type	temp;
	ScaleLikeNoiseColl::value_type					tS;
	PIOErr											MyErr;
	PIOLONG											NData;
	PIODOUBLE										*tData;
	PIODOUBLE										*tDataPiv;
	PIODOUBLE										*tDataPivAux;
	std::string										tObjName(Wstr2Str(FName_));
	int												NScales;

	data_.Storage_.clear();

	if((NData = PIOReadTAB2DObject((void **)&tData,const_cast<char *>(tObjName.c_str()),"PIODOUBLE","tab=*",TAB2Dgrp_))< 0)	
	{HFIDMC_errIO_TAB2D((PIOErr) NData,ERROR_MSG_HFIDMCREADIMGFORMAT,const_cast<char *>(tObjName.c_str()));}

	if(NData != (NRows_ * NColumns_))
	{
		if((MyErr = PIODeleteTAB2DTable(tData,TAB2Dgrp_))!=0)
			HFIDMC_ReportMemoryLeak(MyErr,ERROR_MSG_HFIDMCDELETEIMGFORMAT,tObjName.c_str());
		HFIDMC_errIO(ERROR_COD_HFIDMCERRIMGWRONGSZ,ERROR_MSG_HFIDMCERRIMGWRONGSZ,FName_);
	}

	tDataPivAux = tData;
	for(int i=0;i<NRows_;++i,tDataPivAux += NColumns_)
	{

		tDataPiv								= tDataPivAux;

		temp.DetectID_							= toInt(*tDataPiv++);	
		temp.PatchNumber_						= toInt(*tDataPiv++);	
		temp.PK_BayesDetectStat_				= toInt(*tDataPiv++);	
		temp.GaussianIndex_						= *tDataPiv++;	
		temp.DetectionSigma_					= *tDataPiv++;	
		temp.SrcAmplNormalised_					= *tDataPiv++;	
		temp.UsedSigmaThreshold_				= *tDataPiv++;	
		temp.JF_UsedSigmaThreshold_				= *tDataPiv++;	
		temp.SrcFlux_							= *tDataPiv++;	
		temp.JF_SrcFlux_.Mode_					= *tDataPiv++;	
		temp.JF_lnRhoTh_						= *tDataPiv++;	
		temp.JF_lnRho_							= *tDataPiv++;	
		temp.JF_lnEvidence_						= *tDataPiv++;	
		temp.JF_lnEvidenceErrBar_				= *tDataPiv++;	
		temp.JF_lnFormFactor_					= *tDataPiv++;	
		temp.JF_lnModelRatio_					= *tDataPiv++;	
		temp.Pos_.YCoord_						= *tDataPiv++;	
		temp.Pos_.XCoord_						= *tDataPiv++;	
		temp.Pos_.JF_YCoord_.Mode_				= *tDataPiv++;	
		temp.Pos_.JF_XCoord_.Mode_				= *tDataPiv++;	
		temp.Pos_.JF_YCoord_.Mean_				= *tDataPiv++;	
		temp.Pos_.JF_XCoord_.Mean_				= *tDataPiv++;	
		temp.GalPt_Colat_						= *tDataPiv++;	
		temp.GalPt_Long_						= *tDataPiv++;	
		temp.SrcFlux_mJys_						= *tDataPiv++;	
		temp.SrcCompt_arcmin2_					= *tDataPiv++;	
		temp.RealParams_.RealScale_				= *tDataPiv++;	
		temp.JF_Radius_.Mode_					= *tDataPiv++;	
		temp.Odds_								= *tDataPiv++;	
		temp.ISNR2_								= *tDataPiv++;	
		temp.CollListIndex_						= *tDataPiv++;	
		temp.ErrorBars_.FluxErrorBar_			= *tDataPiv++;	
		temp.ErrorBars_.TotalPosErrorBar_		= *tDataPiv++;	
		temp.ErrorBars_.RadiusErrorBar_			= *tDataPiv++;	
		temp.y0_								= *tDataPiv++;	
		temp.JF_SrcFlux_.Mean_					= *tDataPiv++;	
		temp.JF_Radius_.Mean_					= *tDataPiv++;	
		temp.CoordPt_Lat_						= *tDataPiv++;	
		temp.CoordPt_Long_						= *tDataPiv++;	
		temp.EquPt_Lat_							= *tDataPiv++;	
		temp.EquPt_Long_						= *tDataPiv++;
		temp.ErrorBars_.LowFluxErrorBar_		= *tDataPiv++;	
		temp.ErrorBars_.HighFluxErrorBar_		= *tDataPiv++;	
		temp.ErrorBars_.LowRadiusErrorBar_		= *tDataPiv++;	
		temp.ErrorBars_.HighRadiusErrorBar_		= *tDataPiv++;	
		temp.ErrorBars_.DegenCorr_				= *tDataPiv++;	
		temp.ErrorBars_.DegenOrd_				= *tDataPiv++;	
		temp.ErrorBars_.DegenOrdYErr_			= *tDataPiv++;	
		temp.ErrorBars_.DegenSlope_				= *tDataPiv++;	
		temp.ErrorBars_.DegenSlopeYErr_			= *tDataPiv++;	
//
		temp.ScaleLikeNoise_.clear();

		if((*tDataPiv) > 0.0)
		{
			NScales = Zeus::toInt(*tDataPiv++);
			for(int j=0;j<NScales;++j)
			{
				tS.Scale_ = *tDataPiv++;
				tS.Like_  = *tDataPiv++;
				tS.Noise_ = *tDataPiv++;
				temp.ScaleLikeNoise_.push_back(tS);
			}
		}
//
		for(int i=0;i<QARESULTSSZ;++i)
		{temp.QAResult_[i] = -1.0;}

		if((*tDataPiv ) > 0.0)
		{
			int NQAVals = Zeus::toInt(*tDataPiv++);
			for(int i=0;i<NQAVals;++i)
			{temp.QAResult_[i] = *tDataPiv++;}
		}
//
		data_.Storage_.push_back(temp);
	}

	if((MyErr = PIODeleteTAB2DTable(tData,TAB2Dgrp_))!=0)
	{HFIDMC_ReportMemoryLeak(MyErr,ERROR_MSG_HFIDMCDELETEIMGFORMAT,tObjName.c_str());}

	return true;
}

//
int			ReadNonBlindCatalogueHFIDMC::ReadHeader(void)
{

	CatalogueFormatType::HeaderType	tData;
	PIOErr				MyErr;
	PIOSTRING			dummyStr;
	PIOINT				DummyINT;
	PIODOUBLE			DummyDOUBLE;
	std::string			tFName(Wstr2Str(FName_));

	if((MyErr = PIOReadKeywordObject(((void *) &(DummyINT) ),dummyStr,HFIDMC_NONBLINDHEADER_DEGREES,"PIOINT",const_cast<char *>(tFName.c_str()),TAB2Dgrp_)) == 0)
	{tData.CoordsType_ = DummyINT;}

	if((MyErr = PIOReadKeywordObject(((void *) &(DummyINT) ),dummyStr,HFIDMC_NONBLINDHEADER_COORDSYS,"PIOINT",const_cast<char *>(tFName.c_str()),TAB2Dgrp_)) == 0)
	{tData.CoordSystem_	= static_cast<coordsys>(DummyINT);}

	if((MyErr = PIOReadKeywordObject(((void *) &(DummyDOUBLE) ),dummyStr,HFIDMC_NONBLINDHEADER_EPOCH,"PIODOUBLE",const_cast<char *>(tFName.c_str()),TAB2Dgrp_)) == 0)
	{tData.Epoch_	= DummyDOUBLE;}

	if((MyErr = PIOReadKeywordObject(((void *) &(DummyINT) ),dummyStr,HFIDMC_NONBLINDHEADER_CATNCOL,"PIOINT",const_cast<char *>(tFName.c_str()),TAB2Dgrp_)) == 0)
	{tData.NColumns_	= DummyINT;}

	if((MyErr = PIOReadKeywordObject(((void *) &(DummyINT) ),dummyStr,HFIDMC_NONBLINDHEADER_CATNROW,"PIOINT",const_cast<char *>(tFName.c_str()),TAB2Dgrp_)) == 0)
	{tData.NRows_ = DummyINT;}

	if((MyErr = PIOReadKeywordObject(((void *) &(DummyINT) ),dummyStr,HFIDMC_NONBLINDHEADER_DETCTTYPE,"PIOINT",const_cast<char *>(tFName.c_str()),TAB2Dgrp_)) == 0)
	{tData.DetectionType_	= DummyINT;}

	if((MyErr = PIOReadKeywordObject(((void *) &(DummyINT) ),dummyStr,HFIDMC_NONBLINDHEADER_ESTIMATOR,"PIOINT",const_cast<char *>(tFName.c_str()),TAB2Dgrp_)) == 0)
	{tData.Estimator_ = DummyINT;}

	if((MyErr = PIOReadKeywordObject(((void *) &(DummyINT) ),dummyStr,HFIDMC_NONBLINDHEADER_PRIORTYPE,"PIOINT",const_cast<char *>(tFName.c_str()),TAB2Dgrp_)) == 0)
	{tData.PriorsType_	= DummyINT;}

	if((MyErr = PIOReadKeywordObject(((void *) &(DummyINT) ),dummyStr,HFIDMC_NONBLINDHEADER_SZPROF,"PIOINT",const_cast<char *>(tFName.c_str()),TAB2Dgrp_)) == 0)
	{tData.SZ_params_.SZ_Profile_ = DummyINT;}

	if((MyErr = PIOReadKeywordObject(((void *) &(DummyDOUBLE) ),dummyStr,HFIDMC_NONBLINDHEADER_SZVRATIO,"PIODOUBLE",const_cast<char *>(tFName.c_str()),TAB2Dgrp_)) == 0)
	{tData.SZ_params_.VirialRatio_	= DummyDOUBLE;}

	if((MyErr = PIOReadKeywordObject(((void *) &(DummyDOUBLE) ),dummyStr,HFIDMC_NONBLINDHEADER_ALPHA,"PIODOUBLE",const_cast<char *>(tFName.c_str()),TAB2Dgrp_)) == 0)
	{tData.SZ_params_.MNFW_alpha_	= DummyDOUBLE;}

	if((MyErr = PIOReadKeywordObject(((void *) &(DummyDOUBLE) ),dummyStr,HFIDMC_NONBLINDHEADER_BETA,"PIODOUBLE",const_cast<char *>(tFName.c_str()),TAB2Dgrp_)) == 0)
	{tData.SZ_params_.MNFW_beta_	= DummyDOUBLE;}

	if((MyErr = PIOReadKeywordObject(((void *) &(DummyDOUBLE) ),dummyStr,HFIDMC_NONBLINDHEADER_GAMMA,"PIODOUBLE",const_cast<char *>(tFName.c_str()),TAB2Dgrp_)) == 0)
	{tData.SZ_params_.MNFW_gamma_	= DummyDOUBLE;}

	if((MyErr = PIOReadKeywordObject(((void *) &(DummyDOUBLE) ),dummyStr,HFIDMC_NONBLINDHEADER_C500,"PIODOUBLE",const_cast<char *>(tFName.c_str()),TAB2Dgrp_)) == 0)
	{tData.SZ_params_.MNFW_C500_	= DummyDOUBLE;}

	if((MyErr = PIOReadKeywordObject(((void *) &(DummyDOUBLE) ),dummyStr,HFIDMC_NONBLINDHEADER_CY500CYR500,"PIODOUBLE",const_cast<char *>(tFName.c_str()),TAB2Dgrp_)) == 0)
	{tData.SZ_params_.MNFW_Ratio_CY500CYR500_	= DummyDOUBLE;}

	if((MyErr = PIOReadKeywordObject(((void *) &(DummyDOUBLE) ),dummyStr,HFIDMC_NONBLINDHEADER_FLUXCAL,"PIODOUBLE",const_cast<char *>(tFName.c_str()),TAB2Dgrp_)) == 0)
	{tData.SZ_params_.FluxCalibCte_	= DummyDOUBLE;}

	if((MyErr = PIOReadKeywordObject(((void *) &(DummyDOUBLE) ),dummyStr,HFIDMC_NONBLINDHEADER_RADCAL,"PIODOUBLE",const_cast<char *>(tFName.c_str()),TAB2Dgrp_)) == 0)
	{tData.SZ_params_.RadiusCalCte_	= DummyDOUBLE;}

	if((MyErr = PIOReadKeywordObject(((void *) &(DummyINT) ),dummyStr,HFIDMC_NONBLINDHEADER_COLLSTSZ,"PIOINT",const_cast<char *>(tFName.c_str()),TAB2Dgrp_)) == 0)
	{tData.CollLstSz_ = DummyINT;}

	data_.Header_ = tData;
	return true;
}
//
int			ReadNonBlindCatalogueHFIDMC::ReadBody(void)
{
	PIOErr				MyErr;
	PIOLONG				FRow;
	PIOLONG				NData;
	PIOLONG				ReadData;
	PIODOUBLE			*tData;
	std::string			tObjName(Wstr2Str(FName_));

	if(	(PIOGetTAB2DLines(&FRow,&NData,const_cast<char *>(tObjName.c_str()),TAB2Dgrp_) != 0) ||
		((ReadData = PIOReadTAB2DObject((void **)&tData,const_cast<char *>(tObjName.c_str()),"PIODOUBLE","tab=*",TAB2Dgrp_))< 0)
		)
	{HFIDMC_errIO_TAB2D((PIOErr) NData,ERROR_MSG_HFIDMCREADIMGFORMAT,FName_);}

	if(ReadData == 0)
		return true;

	if((data_.Header_.NRows_ != (NData + 1)))
	{
		wchar_t buf[BUFFERMAXCHAR];
		
		PRINTINTOBUFFERFUNCT(buf,BUFFERMAXCHAR,L"Warning! Number of lines in table %ls is different !",tObjName.c_str());

		(Zeus::ConManager::Instance())->PrintStr2Console(buf);
	}	

	int	res(Convert2CatFormat(data_.Header_.NRows_,tData));

	if((MyErr = PIODeleteIMG2DTable(tData,TAB2Dgrp_))!=0)
		HFIDMC_ReportMemoryLeak(MyErr,ERROR_MSG_HFIDMCDELETEIMGFORMAT,tObjName.c_str());

	if(res < 0)
	{HFIDMC_errIO(ERROR_COD_HFIDMCERRIMGWRONGSZ,ERROR_MSG_HFIDMCERRIMGWRONGSZ,FName_);}		

	return true;
}
//
int			ReadQA_ProfilesLstHFIDMC::ReadBody(void)
{
	PIOErr				MyErr;
	PIOLONG				ReadData;
	PIODOUBLE			*tData;
	std::string			tObjName(Wstr2Str(FName_));

	if((ReadData = PIOReadTAB2DObject((void **)&tData,const_cast<char *>(tObjName.c_str()),"PIODOUBLE","tab=*",TAB2Dgrp_))< 0)
	{HFIDMC_errIO_TAB2D((PIOErr) ReadData,ERROR_MSG_HFIDMCREADIMGFORMAT,FName_);}

	if(ReadData == 0)
		return true;

// 5 is the number of columns

	if((ReadData % 5) || (Convert2CatFormat((ReadData / 5),tData)< 0))
	{
		if((MyErr = PIODeleteIMG2DTable(tData,TAB2Dgrp_))!=0)
			HFIDMC_ReportMemoryLeak(MyErr,ERROR_MSG_HFIDMCDELETEIMGFORMAT,tObjName.c_str());
		HFIDMC_errIO(ERROR_COD_HFIDMCERRIMGWRONGSZ,ERROR_MSG_HFIDMCERRIMGWRONGSZ,FName_);
	}
	else
	{
		if((MyErr = PIODeleteIMG2DTable(tData,TAB2Dgrp_))!=0)
			HFIDMC_ReportMemoryLeak(MyErr,ERROR_MSG_HFIDMCDELETEIMGFORMAT,tObjName.c_str());	
	}
	
	return true;
}
//
int			ReadQA_ColLstHFIDMC::ReadBody(void)
{
	PIOErr				MyErr;
	PIOLONG				ReadData;
	PIODOUBLE			*tData;
	std::string			tObjName(Wstr2Str(FName_));

	if((ReadData = PIOReadTAB2DObject((void **)&tData,const_cast<char *>(tObjName.c_str()),"PIODOUBLE","tab=*",TAB2Dgrp_))< 0)
	{HFIDMC_errIO_TAB2D((PIOErr) ReadData,ERROR_MSG_HFIDMCREADIMGFORMAT,FName_);}

	if(ReadData == 0)
		return true;

// 43 is the number of columns

	if((ReadData % 43) || (Convert2CatFormat((ReadData / 43),tData)< 0))
	{
		if((MyErr = PIODeleteIMG2DTable(tData,TAB2Dgrp_))!=0)
			HFIDMC_ReportMemoryLeak(MyErr,ERROR_MSG_HFIDMCDELETEIMGFORMAT,tObjName.c_str());
		HFIDMC_errIO(ERROR_COD_HFIDMCERRIMGWRONGSZ,ERROR_MSG_HFIDMCERRIMGWRONGSZ,FName_);
	}
	else
	{
		if((MyErr = PIODeleteIMG2DTable(tData,TAB2Dgrp_))!=0)
			HFIDMC_ReportMemoryLeak(MyErr,ERROR_MSG_HFIDMCDELETEIMGFORMAT,tObjName.c_str());	
	}
	
	return true;
}
//
int			ReadQA_ProfilesLstHFIDMC::Convert2CatFormat(PIOLONG NData,PIODOUBLE* tData)
{
// 5 is the number of columns

	QA_ProfilesLstType::StorageType::value_type	temp;
	PIODOUBLE*		piv;

	for(int i=0;i<NData;++i,tData += 5)
	{
		piv										= tData;

		temp.Alpha_ =			*piv++;
		temp.Beta_  =			*piv++;
		temp.Gamma_ =			*piv++;
		temp.C500_  =			*piv++;
		temp.ConvYtot2Y500_ =	*piv++;
		
		data_.Storage_.push_back(temp);
	}

	return true;
}
//
int			ReadQA_ColLstHFIDMC::Convert2CatFormat(PIOLONG NData,PIODOUBLE* tData)
{
// 43 is the number of columns

	QA_CltLstCatalogueType::StorageType::value_type	temp;
	PIODOUBLE*		piv;

	for(int i=0;i<NData;++i,tData += 43)
	{
		piv										= tData;

		temp.IN_theta=			*piv++;
		temp.IN_phi=			*piv++;
		temp.IN_CyR500=			*piv++;
		temp.IN_T500=			*piv++;
		temp.IN_M500=			*piv++;
		temp.IN_z=				*piv++;
		temp.IN_Vrec=			*piv++;
		temp.IN_detectable=		*piv++;
		temp.Cat_GLAT=			*piv++;
		temp.Cat_GLON=			*piv++;
		temp.Cat_yc=			*piv++;
		temp.Cat_yc_error=		*piv++;
		temp.Cat_CY5R500=		*piv++;
		temp.Cat_LOWER_ERR_CY5R500F= *piv++;
		temp.Cat_UPPER_ERR_CY5R500F= *piv++;
		temp.Cat_T500=			*piv++;
		temp.Cat_LOWER_ERR_T500F= *piv++;
		temp.Cat_UPPER_ERR_T500F= *piv++;
		temp.Cat_SNR=			*piv++;
		temp.Cat_Patch_No=		*piv++;
		temp.Cat_DEGEN_PARAM=	*piv++;
		temp.Cat_SNR_30GHZ=		*piv++;
		temp.Cat_SNR_44GHZ=		*piv++;
		temp.Cat_SNR_70GHZ=		*piv++;
		temp.Cat_SNR_100GHZ=	*piv++;
		temp.Cat_SNR_143GHZ=	*piv++;
		temp.Cat_SNR_217GHZ=	*piv++;
		temp.Cat_SNR_353GHZ=	*piv++;
		temp.Cat_SNR_545GHZ=	*piv++;
		temp.Cat_SNR_857GHZ=	*piv++;
		temp.Cat_FLUX5R500_30GHZ= *piv++;
		temp.Cat_FLUX5R500_44GHZ= *piv++;
		temp.Cat_FLUX5R500_70GHZ= *piv++;
		temp.Cat_FLUX5R500_100GHZ= *piv++;
		temp.Cat_FLUX5R500_143GHZ= *piv++;
		temp.Cat_FLUX5R500_217GHZ= *piv++;
		temp.Cat_FLUX5R500_353GHZ= *piv++;
		temp.Cat_FLUX5R500_545GHZ= *piv++;
		temp.Cat_FLUX5R500_857GHZ= *piv++;
		temp.Cat_CHI2=				*piv++;
		temp.Cat_ERR_RAD=			*piv++;
		temp.Cat_CY5R500F=		*piv++;
		temp.Cat_T500F=			*piv++;
		
		data_.Storage_.push_back(temp);
	}

	return true;
}
//
int			ReadNonBlindCatalogueHFIDMC::Convert2CatFormat(PIOLONG NData,PIODOUBLE* tData)
{

	CatalogueFormatType::StorageType::value_type	temp;
	const int NColumns(data_.Header_.NColumns_);
	ScaleLikeNoiseColl::value_type	tS;
	PIODOUBLE*		piv;


	for(int i=0;i<NData;++i,tData += NColumns)
	{
		piv										= tData;

		temp.ID_								= Zeus::toInt(*piv++);
		temp.Patch_								= Zeus::toInt(*piv++);
		temp.lnRho_								= *piv++;
		temp.NormalAmpl_						= *piv++;
		temp.DetectSigma_						= *piv++;
		temp.GalLatDegs_						= *piv++;
		temp.GalLongDegs_						= *piv++;
		temp.PatchGalLatDegs_					= *piv++;
		temp.PatchGalLongDegs_					= *piv++;
		temp.PatchSpin_							= *piv++;
		temp.FluxCompt_							= *piv++;
		temp.FluxComptGLRT_						= *piv++;
		temp.Radius_							= *piv++;
		temp.RadiusGLRT_						= *piv++;
		temp.ErrorBars_.FluxErrorBar_			= *piv++;
		temp.ErrorBars_.RadiusErrorBar_			= *piv++;
		temp.ErrorBars_.TotalPosErrorBar_		= *piv++;
		temp.Gaussianity_						= *piv++;
		temp.lnEvidence_						= *piv++;
		temp.lnPenaltySrc_						= *piv++;
		temp.PatchMF_sigmaSqr_					= *piv++;
		temp.SZ_ConversionCte_					= *piv++;
		temp.ErrorBars_.LowFluxErrorBar_		= *piv++;
		temp.ErrorBars_.HighFluxErrorBar_		= *piv++;
		temp.ErrorBars_.LowRadiusErrorBar_		= *piv++;
		temp.ErrorBars_.HighRadiusErrorBar_		= *piv++;
		temp.ErrorBars_.DegenCorr_				= *piv++;
		temp.ErrorBars_.DegenOrd_				= *piv++;
		temp.ErrorBars_.DegenOrdYErr_			= *piv++;
		temp.ErrorBars_.DegenSlope_				= *piv++;
		temp.ErrorBars_.DegenSlopeYErr_			= *piv++;
		temp.CollLstIndex_						= Zeus::toInt(*piv++);
		
		temp.ScaleLikeNoise_.clear();

		if((*piv ) > 0.0)
		{
			int NScales = Zeus::toInt(*piv++);
			for(int i=0;i<NScales;++i)
			{
				tS.Scale_ = *piv++;
				tS.Like_  = *piv++;
				tS.Noise_ = *piv++;
				temp.ScaleLikeNoise_.push_back(tS);
			}
		}

		for(int i=0;i<QARESULTSSZ;++i)
		{temp.QAResult_[i] = -1.0;}

		if((*piv ) > 0.0)
		{
			int NQAVals = Zeus::toInt(*piv++);
			for(int i=0;i<NQAVals;++i)
			{temp.QAResult_[i] = *piv++;}
		}

		data_.Storage_.push_back(temp);
	}
	return true;
}
//
void		CatalogueWriterHFIDMC::Move2Array2D(const CatalogueFormatType & DataIN,LArr2D<double> & DataOUT)
{
	double		* const endOUT(DataOUT.end());
	double		*pivOUT(DataOUT.begin());
	const int	NFields(NColumns_);

	CatalogueFormatType::StorageType::const_iterator		pivIN(DataIN.Storage_.begin());
	CatalogueFormatType::StorageType::const_iterator		const endIN(DataIN.Storage_.end());

	double	*pivAuxOUT;

	for(;pivIN != endIN;++pivIN,pivOUT += NFields)
	{
		pivAuxOUT			= pivOUT;

		*pivAuxOUT++		= static_cast<double>(pivIN->ID_);
		*pivAuxOUT++		= static_cast<double>(pivIN->Patch_);

		if((pivIN->NormalAmpl_ >= 0) && (pivIN->GalLongDegs_ > -1.0e10))
		{
			*pivAuxOUT++		= pivIN->lnRho_;
			*pivAuxOUT++		= pivIN->NormalAmpl_;
			*pivAuxOUT++		= pivIN->DetectSigma_;
			*pivAuxOUT++		= pivIN->GalLatDegs_;
			*pivAuxOUT++		= pivIN->GalLongDegs_;
			*pivAuxOUT++		= pivIN->PatchGalLatDegs_;
			*pivAuxOUT++		= pivIN->PatchGalLongDegs_;
			*pivAuxOUT++		= pivIN->PatchSpin_;
			*pivAuxOUT++		= pivIN->FluxCompt_;
			*pivAuxOUT++		= pivIN->FluxComptGLRT_;
			*pivAuxOUT++		= pivIN->Radius_;
			*pivAuxOUT++		= pivIN->RadiusGLRT_;
			*pivAuxOUT++		= pivIN->ErrorBars_.FluxErrorBar_;
			*pivAuxOUT++		= pivIN->ErrorBars_.RadiusErrorBar_;
			*pivAuxOUT++		= pivIN->ErrorBars_.TotalPosErrorBar_;
			*pivAuxOUT++		= pivIN->Gaussianity_;
			*pivAuxOUT++		= pivIN->lnEvidence_;
			*pivAuxOUT++		= pivIN->lnPenaltySrc_;
			*pivAuxOUT++		= pivIN->PatchMF_sigmaSqr_;
			*pivAuxOUT++		= pivIN->SZ_ConversionCte_;
			*pivAuxOUT++		= pivIN->ErrorBars_.LowFluxErrorBar_;
			*pivAuxOUT++		= pivIN->ErrorBars_.HighFluxErrorBar_;
			*pivAuxOUT++		= pivIN->ErrorBars_.LowRadiusErrorBar_;
			*pivAuxOUT++		= pivIN->ErrorBars_.HighRadiusErrorBar_;
			*pivAuxOUT++		= pivIN->ErrorBars_.DegenCorr_;
			*pivAuxOUT++		= pivIN->ErrorBars_.DegenOrd_;
			*pivAuxOUT++		= pivIN->ErrorBars_.DegenOrdYErr_;
			*pivAuxOUT++		= pivIN->ErrorBars_.DegenSlope_;
			*pivAuxOUT++		= pivIN->ErrorBars_.DegenSlopeYErr_;
			*pivAuxOUT++		= static_cast<double>(pivIN->CollLstIndex_);

			if(pivIN->ScaleLikeNoise_.empty())
			{
				*pivAuxOUT++ = -1.0;
				continue;
			}

			*pivAuxOUT++		=  static_cast<double>(pivIN->ScaleLikeNoise_.size());

			ScaleLikeNoiseColl::const_iterator	pivS(pivIN->ScaleLikeNoise_.begin());
			ScaleLikeNoiseColl::const_iterator	const endS(pivIN->ScaleLikeNoise_.end());
			for(;pivS != endS;++pivS)
			{
				*pivAuxOUT++ = pivS->Scale_;
				*pivAuxOUT++ = pivS->Like_;
				*pivAuxOUT++ = pivS->Noise_;
			}

			if(pivIN->QAResult_[0] < 0.0)
			{
				*pivAuxOUT++ = -1.0;
				continue;
			}

			*pivAuxOUT++		=  static_cast<double>(QARESULTSSZ);

			for (int i=0;i<QARESULTSSZ;++i)
			{
				*pivAuxOUT++ = pivIN->QAResult_[i];
			}
		}
	}
}
//
void		CatWriterOutputHFIDMC_Orange10::CreateScalesGroup(void)
{

	PIOErr			MyErr;
	std::string		tGrpNameStr(Wstr2Str(ScalesGrpName_));
	int retries(MAXRETRIES);

_retry:
	if(PIOCheckGroup(const_cast<char *>(tGrpNameStr.c_str())))
	{
		if(MyErr = PIOCreateTAB2DGrp(const_cast<char *>(tGrpNameStr.c_str()),3 /*NColumns_*/)) 
		{
			if(--retries >= 0)
			{
				MySleep(1);
				goto _retry;
			}
			HFIDMC_errIO(ERROR_COD_HFIDMCERRCREATFL,ERROR_MSG_HFIDMCERRCREATFL,ScalesGrpName_);
		}
	}

	if (!(ScalesTAB2Dgrp_ = PIOOpenTAB2DGrp(const_cast<char *>(tGrpNameStr.c_str()),"w")))
	{HFIDMC_errIO(ERROR_COD_HFIDMCERROPENFL,ERROR_MSG_HFIDMCERROPENFL,ScalesGrpName_);}
}
//
void		CatWriterOutputHFIDMC_Orange10::Move2Array2D(const OutputFormatType & DataIN,LArr2D<double> & DataOUT)
{
	double			*pivOUT(DataOUT.begin());
	double			* const endOUT(DataOUT.end());
	const			int NFields(GetNColumns());

	double	*pivAuxOUT;

	OutputFormatType::StorageType::const_iterator		pivIN(DataIN.Storage_.begin());
	OutputFormatType::StorageType::const_iterator		const endIN(DataIN.Storage_.end());

	for(int SrcIndex=0;pivIN != endIN;++pivIN,pivOUT += NFields,++SrcIndex)
	{

		if((pivIN->Cat_.NormalAmpl_ < 0) || (pivIN->Cat_.GalLongDegs_ < -1.0e10))
			continue;

		pivAuxOUT			= pivOUT;
		
		*pivAuxOUT++		= PIOVER2 - (pivIN->Cat_.GalLatDegs_/RAD2DEGREE);
		// 1
		*pivAuxOUT++		= pivIN->Cat_.GalLongDegs_ / RAD2DEGREE;
		// 2
		*pivAuxOUT++		= Get_y0(DataIN,pivIN);
		// 3
		*pivAuxOUT++		= Get_y0Err(DataIN,pivIN);
		// 4
		*pivAuxOUT++		= pivIN->Cat_.FluxComptGLRT_;
		// 5
		*pivAuxOUT++		= pivIN->Cat_.ErrorBars_.LowFluxErrorBar_;
		// 6
		*pivAuxOUT++		= pivIN->Cat_.ErrorBars_.HighFluxErrorBar_;
		// 7
		*pivAuxOUT++		= pivIN->Cat_.RadiusGLRT_  *  DataIN.Header_.PCCHeader_.SZ_params_.MNFW_C500_;
		// 8
		*pivAuxOUT++		= pivIN->Cat_.ErrorBars_.LowRadiusErrorBar_ * DataIN.Header_.PCCHeader_.SZ_params_.MNFW_C500_;
		// 9
		*pivAuxOUT++		= pivIN->Cat_.ErrorBars_.HighRadiusErrorBar_ * DataIN.Header_.PCCHeader_.SZ_params_.MNFW_C500_;
		// 10
//		*pivAuxOUT++		= ((pivIN->Cat_.NormalAmpl_ < pivIN->Cat_.DetectSigma_)?pivIN->Cat_.NormalAmpl_:pivIN->Cat_.DetectSigma_);
		*pivAuxOUT++		= pivIN->Cat_.NormalAmpl_;
		// 11
//		*pivAuxOUT++		= static_cast<double>(pivIN->Cat_.Patch_);
		*pivAuxOUT++		= static_cast<double>(SrcIndex);
		// 12
		//*pivAuxOUT++		= pivIN->Cat_.ErrorBars_.DegenSlopeYErr_;
		*pivAuxOUT++		= pivIN->Ext_.SNRF_;
		// 13
		*pivAuxOUT++		= pivIN->Ext_.SNR_030_;
		// 14
		*pivAuxOUT++		= pivIN->Ext_.SNR_044_;
		// 15
		*pivAuxOUT++		= pivIN->Ext_.SNR_070_;
		// 16
		*pivAuxOUT++		= pivIN->Ext_.SNR_100_;
		// 17
		*pivAuxOUT++		= pivIN->Ext_.SNR_143_;
		// 18 
		*pivAuxOUT++		= pivIN->Ext_.SNR_217_;
		// 19
		*pivAuxOUT++		= pivIN->Ext_.SNR_353_;
		// 20
		*pivAuxOUT++		= pivIN->Ext_.SNR_545_;
		// 21
		*pivAuxOUT++		= pivIN->Ext_.SNR_857_;
		// 22
		*pivAuxOUT++		= pivIN->Ext_.Flux5R500_030_;
		// 23
		*pivAuxOUT++		= pivIN->Ext_.Flux5R500_044_;
		// 24
		*pivAuxOUT++		= pivIN->Ext_.Flux5R500_070_;
		// 25
		*pivAuxOUT++		= pivIN->Ext_.Flux5R500_100_;
		// 26
		*pivAuxOUT++		= pivIN->Ext_.Flux5R500_143_;
		// 27
		*pivAuxOUT++		= pivIN->Ext_.Flux5R500_217_;
		// 28
		*pivAuxOUT++		= pivIN->Ext_.Flux5R500_353_;
		// 29
		*pivAuxOUT++		= pivIN->Ext_.Flux5R500_545_;
		// 30
		*pivAuxOUT++		= pivIN->Ext_.Flux5R500_857_;
		// 31
		*pivAuxOUT++		= pivIN->Ext_.CHI2_;
		// 32
		*pivAuxOUT++		= pivIN->Cat_.ErrorBars_.TotalPosErrorBar_;
		// 33
		*pivAuxOUT++		= pivIN->Cat_.FluxCompt_;
		// 34
		*pivAuxOUT++		= pivIN->Cat_.Radius_  * DataIN.Header_.PCCHeader_.SZ_params_.MNFW_C500_;
		// 35

		if(!(pivIN->Cat_.ScaleLikeNoise_.empty()))
		{
#ifdef SCALES1DINDB
			CreateScalesTAB2D(pivIN);
#endif
		}
	}
}
//
void		CatWriterOutputHFIDMC_QA_Contours::Move2Array2D(const OutputFormatType & DataIN,LArr2D<double> & DataOUT)
{
	double			*pivOUT(DataOUT.begin());
	double			* const endOUT(DataOUT.end());
	const			int NFields(GetNColumns());

	double	*pivAuxOUT;

	OutputFormatType::StorageType::const_iterator		pivIN(DataIN.Storage_.begin());
	OutputFormatType::StorageType::const_iterator		const endIN(DataIN.Storage_.end());

	for(;pivIN != endIN;++pivIN,pivOUT += NFields)
	{

		pivAuxOUT			= pivOUT;

		for(int i=0;i<NFields;++i)
		{
			*pivAuxOUT++ = pivIN->Cat_.QAResult_[i];
		}
	}
}
//
void		CatWriterOutputHFIDMC_Orange10::CreateScalesTAB2D(OutputFormatType::StorageType::const_iterator ptr)
{
	PIOLONG			MyErr;
	char			buffer[BUFFERMAXCHAR];

	if(ptr->Cat_.ScaleLikeNoise_.empty())
		return;

	sprintf(buffer,"SZcolat%dlong%dpatchNo%d",Zeus::toInt((90.0-ptr->Cat_.GalLatDegs_) * 60.0),
		Zeus::toInt(ptr->Cat_.GalLongDegs_ * 60.0),ptr->Cat_.Patch_);

// create the object
	if(PIOCheckObject(const_cast<char *>(buffer),ScalesTAB2Dgrp_))
	{
		if(MyErr = ((PIOErr)PIOCreateTAB2DObject(const_cast<char *>(buffer),"PIODOUBLE",ScalesTAB2Dgrp_))) 
		{HFIDMC_errIO_TAB2D(MyErr,ERROR_MSG_HFIDMCCREATEIMGFORMAT,Achar2Wstr(buffer));}
	}
// write data into it

	LArr2D<double>		dataOUT(ptr->Cat_.ScaleLikeNoise_.size() * 3, 3,SZCATALOGUEDEFAULTVALUE);

	double		*pivOUT(dataOUT.begin());

	ScaleLikeNoiseColl::const_iterator	pivS(ptr->Cat_.ScaleLikeNoise_.begin());
	ScaleLikeNoiseColl::const_iterator	const endS(ptr->Cat_.ScaleLikeNoise_.end());

	for(;pivS != endS;++pivS)
	{
		*pivOUT++ = pivS->Scale_;
		*pivOUT++ = pivS->Like_;
		*pivOUT++ = pivS->Noise_;
	}

// commit to DB
	PIOSTRING	command;
	
	sprintf(command,"tab=0:%" PRId64 ",0:%" PRId64,(PIOLONG)(2),(PIOLONG)(ptr->Cat_.ScaleLikeNoise_.size() - 1));

	if((MyErr = PIOWriteTAB2DObject(dataOUT.begin(),const_cast<char *>(buffer),"PIODOUBLE",command,ScalesTAB2Dgrp_)) < 0) 
	{
		HFIDMC_errIO_TAB2D((PIOErr) MyErr,ERROR_MSG_HFIDMCWRITEIMGFORMAT,Achar2Wstr(buffer));
	}
}
//
void		CatWriterOutputHFIDMC_Orange18::Move2Array2D(const OutputFormatType & DataIN,LArr2D<double> & DataOUT)
{
	double			*pivOUT(DataOUT.begin());
	double			* const endOUT(DataOUT.end());
	const			int NFields(GetNColumns());

	double	*pivAuxOUT;

	OutputFormatType::StorageType::const_iterator		pivIN(DataIN.Storage_.begin());
	OutputFormatType::StorageType::const_iterator		const endIN(DataIN.Storage_.end());

	for(;pivIN != endIN;++pivIN,pivOUT += NFields)
	{

		if((pivIN->Cat_.NormalAmpl_ < 0) || (pivIN->Cat_.GalLongDegs_ < -1.0e10))
			continue;

		pivAuxOUT			= pivOUT;

		*pivAuxOUT++		= PIOVER2 - (pivIN->Cat_.GalLatDegs_/RAD2DEGREE);
		// 1
		*pivAuxOUT++		= pivIN->Cat_.GalLongDegs_ / RAD2DEGREE;
		// 2
		*pivAuxOUT++;
		// 3 MHW only
		*pivAuxOUT++;
		// 4 MHW only
		*pivAuxOUT++		= pivIN->Cat_.FluxCompt_ / 1000.0; // Flux in Jys
		// 5
		*pivAuxOUT++		= pivIN->Cat_.ErrorBars_.FluxErrorBar_  / 1000.0;  // Flux in Jys
		// 6
		*pivAuxOUT++		= ((pivIN->Cat_.NormalAmpl_ < pivIN->Cat_.DetectSigma_) ? pivIN->Cat_.NormalAmpl_:pivIN->Cat_.DetectSigma_);
//		*pivAuxOUT++		= pivIN->Cat_.NormalAmpl_;
		// 7
		*pivAuxOUT++		= pivIN->Cat_.ErrorBars_.TotalPosErrorBar_;
		// 8
		*pivAuxOUT++;
		// chi_sq
		// 9 MHW
		*pivAuxOUT++;
		// flux_ratio
		// 10 MHW
		*pivAuxOUT++		= pivIN->Cat_.Radius_;
		// 11
		*pivAuxOUT++		= pivIN->Cat_.ErrorBars_.RadiusErrorBar_;
		// 12
		*pivAuxOUT++		= pivIN->Cat_.Gaussianity_;
		// 13
		*pivAuxOUT++		= -dBCONV * Log1plusLog(pivIN->Cat_.lnRho_);
		// 14
 	}
}
//
int			HFIDMC_RemoveObject(const std::wstring& Dir,const std::wstring& file)
{
	PIOGroup	*Grp;
	PIOErr		MyErr;
	std::string	tDir(Wstr2Str(Dir));
	std::string	tFile(Wstr2Str(file));

	if(!(Grp = PIOOpenVoidGrp(const_cast<char *>(tDir.c_str()),"w")))
		return -1;
	
	MyErr = PIODeleteObject(const_cast<char *>(tFile.c_str()),Grp);

	PIOCloseVoidGrp(&Grp);

	return MyErr;
}
//
#endif //HFIDMC
//
#ifdef	LFIDPC
//
void	LFIDPC_ParameterFileReader::ReadValuesFromLFI_Pipeline(void)
{
	int	NFreq,NPrior,NFiles;

	LFIDPC_INSERTVALUESTR(GLOBALID_DIRIN);
	LFIDPC_INSERTVALUESTR(GLOBALID_DIRINMAPS);
	LFIDPC_INSERTVALUESTR(GLOBALID_POINTINGS);
	LFIDPC_INSERTVALUESTR(GLOBALID_MASKFILENAME);
	LFIDPC_INSERTVALUESTR(GLOBALID_DATASETNAME);
	LFIDPC_INSERTVALUE(GLOBALID_NFREQS,int);
	NFreq	= Params_->find<int>(std::string(MAKECHARARR(GLOBALID_NFREQS)));
	LFIDPC_INSERTVALUE(GLOBALID_NPRIORPLANES,int);
	NPrior	= Params_->find<int>(std::string(MAKECHARARR(GLOBALID_NPRIORPLANES)));
	LFIDPC_INSERTVALUE(GLOBALID_NOTALIGNEDOBJS,int);
	LFIDPC_INSERTVALUE(GLOBALID_PWSSEARCHPOS,int);
	LFIDPC_INSERTVALUE(GLOBALID_USE2DFORMULA,int);
	LFIDPC_INSERTVALUE(GLOBALID_NSCALEBINS,int);
	LFIDPC_INSERTVALUE(GLOBALID_USEBOUNDS,int);
	LFIDPC_INSERTVALUE(GLOBALID_ASSESS_TYPE,int);
	LFIDPC_INSERTVALUE(GLOBALID_BOUNDTOL,double);
	LFIDPC_INSERTVALUE(GLOBALID_SIGMA_THRESHOLD,double);
	LFIDPC_INSERTVALUE(GLOBALID_J_THRESHOLD,double);
	LFIDPC_INSERTVALUE(GLOBALID_SZDETECTION,int);
	LFIDPC_INSERTVALUE(GLOBALID_SZVIRIALRATIO,double);
	LFIDPC_INSERTVALUE(GLOBALID_CATALOGUESIGMA,double);
	LFIDPC_INSERTVALUE(GLOBALID_APODIZEMAPS,int);
	LFIDPC_INSERTVALUE(GLOBALID_OUTPUTCOORDS,int);
	LFIDPC_INSERTVALUE(GLOBALID_CATMERGE_AVGT,int);
	LFIDPC_INSERTVALUE(GLOBALID_GALACTICSIGMA,double);
	LFIDPC_INSERTVALUE(GLOBALID_OUTUNITS,int);
	LFIDPC_INSERTVALUESTR(GLOBALID_FILTERINGSCALES);
	LFIDPC_INSERTVALUE(GLOBALID_SZPROFILE,int);
	LFIDPC_INSERTVALUE(GLOBALID_CACHESZ,double);
	LFIDPC_INSERTVALUE(GLOBALID_JFESTIMATTYPE,int);
	LFIDPC_INSERTVALUE(GLOBALID_OUTPUTRADIUSCAL,double);
	LFIDPC_INSERTVALUE(GLOBALID_FLUXCALIBCTE,double);
	LFIDPC_INSERTVALUE(GLOBALID_SRCMAXSCALE,double);
	LFIDPC_INSERTVALUE(GLOBALID_TWOSTEPSDETECT,int);
	LFIDPC_INSERTVALUE(GLOBALID_NONBLINDDETECTION,int);
	LFIDPC_INSERTVALUE(GLOBALID_MELTMAXDIST,double);
	LFIDPC_INSERTVALUE(GLOBALID_MELTMAXDISTSZLIMIT,double);
	LFIDPC_INSERTVALUE(GLOBALID_FLUXTHRESHOLD,double);
	LFIDPC_INSERTVALUE(GLOBALID_OUTPUTLAT,int);
	LFIDPC_INSERTVALUE(GLOBALID_OUTPUTDEGREES,int);
	LFIDPC_INSERTVALUE(GLOBALID_GALACTICCUT,double);
	LFIDPC_INSERTVALUE(GLOBALID_PRIORFLUXMIN,double);
	LFIDPC_INSERTVALUE(GLOBALID_PRIORFLUXMAX,double);
	LFIDPC_INSERTVALUE(GLOBALID_PRIORRADDISTMIN,double);
	LFIDPC_INSERTVALUE(GLOBALID_PRIORRADDISTMAX,double);
	LFIDPC_INSERTVALUE(GLOBALID_PRIORMASSMIN,double);
	LFIDPC_INSERTVALUE(GLOBALID_PRIORMASSMAX,double);
	LFIDPC_INSERTVALUE(GLOBALID_PRIORMAXZ,double);
	LFIDPC_INSERTVALUE(GLOBALID_PRIORTEMPMIN,double);
	LFIDPC_INSERTVALUE(GLOBALID_PRIORTEMPMAX,double);
	LFIDPC_INSERTVALUE(GLOBALID_PRIORGASMASSRATIO,double);
	LFIDPC_INSERTVALUE(GLOBALID_PRIORSIGMA8,double);
	LFIDPC_INSERTVALUE(GLOBALID_USEPRESSSCHT,int);
	LFIDPC_INSERTVALUE(GLOBALID_MN_NLIVEPOINTS,int);
	LFIDPC_INSERTVALUE(GLOBALID_MN_NINDIVSAMPLES,int);
	LFIDPC_INSERTVALUE(GLOBALID_MN_FRACTOLEV,double);
	LFIDPC_INSERTVALUE(GLOBALID_MN_XTRAENLFACT,double);
	LFIDPC_INSERTVALUE(GLOBALID_PRIORFLUXEXP,double);
	LFIDPC_INSERTVALUE(GLOBALID_PRIORSRCSCALEMIN,double);
	LFIDPC_INSERTVALUE(GLOBALID_PRIORSRCSCALEMAX,double);
	LFIDPC_INSERTVALUE(GLOBALID_PRIORSRCSCALEEXP,double);
	LFIDPC_INSERTVALUE(GLOBALID_PRIORAVSRCPATCH,double);
	LFIDPC_INSERTVALUE(GLOBALID_OUTPURITY,double);
	LFIDPC_INSERTVALUE(GLOBALID_FIRSTPATCH,int);
	LFIDPC_INSERTVALUE(GLOBALID_LASTPATCH,int);
	LFIDPC_INSERTVALUESTR(GLOBALID_INTCATNAME);
	LFIDPC_INSERTVALUESTR(GLOBALID_FINALCATNAME);
	LFIDPC_INSERTVALUE(GLOBALID_SZPROFALPHA,double);
	LFIDPC_INSERTVALUE(GLOBALID_SZPROFBETA,double);
	LFIDPC_INSERTVALUE(GLOBALID_SZPROFGAMMA,double);
	LFIDPC_INSERTVALUE(GLOBALID_SZPROF_C500,double);
	LFIDPC_INSERTVALUE(GLOBALID_SZPROFCY500CYR500,double);

	LFIDPC_INSERTVALUE(MCID_NOFFILES,int);
	NFiles	= Params_->find<int>(std::string(MAKECHARARR(MCID_NOFFILES)));
	LFIDPC_INSERTVALUE(MCID_GALACTICCUT,double);
	LFIDPC_INSERTVALUE(MCID_MAPCUTTERID,int);
	LFIDPC_INSERTVALUESTR(MCID_COMMANDFLAG);
	LFIDPC_INSERTVALUE(MCID_PATCHSZ,int);
	LFIDPC_INSERTVALUE(MCID_PBORDER,int);
	LFIDPC_INSERTVALUE(MCID_NLINEPTS,int);
	LFIDPC_INSERTVALUE(MCID_NSIDE,int);
	LFIDPC_INSERTVALUE(MCID_EPOCH,double);
	LFIDPC_INSERTVALUE(MCID_MASLENLAGTH,double);
	LFIDPC_INSERTVALUE(MCID_RMSREJECTLEVEL,double);
	LFIDPC_INSERTVALUE(MCID_MASKFWHM,double);
	LFIDPC_INSERTVALUE(MCID_PERCENTREJECT,double);
	LFIDPC_INSERTVALUE(MCID_PTGSCOORDSYS,int);
	LFIDPC_INSERTVALUESTR(MCID_MASKREJECTNAME);
	LFIDPC_INSERTVALUESTR(MCID_MASKREMOVENAME);
	LFIDPC_INSERTVALUESTR(MCID_NONBLINDPTGSFILE);

	std::wstring	StrTemp;
	wchar_t			buffer[DOUBLETXTMAXSZ];
	std::wstring	ID;
	std::wstring	ID1;
	std::wstring	ID2;
	std::wstring	ID3;
	std::wstring	IDFWHM;
	std::wstring	mc_freq;
	std::wstring	mc_units;
	std::wstring	mc_coordsys;
	std::wstring	mc_maptype;
	std::wstring	mc_unitconvfact;
	std::wstring	mc_processed;
	std::wstring	mc_file;
	std::wstring	NumberStr;

	for(int i=0;i<NFiles;++i)
	{
		mc_file = mc_freq = mc_units = mc_coordsys = mc_maptype = mc_unitconvfact = mc_processed = std::wstring(MC_FILEIDSTR);
		Zeus::PutNumber2Txt(buffer,i);
		NumberStr = std::wstring(buffer);
		mc_freq			+= NumberStr;
		mc_units		+= NumberStr;
		mc_coordsys		+= NumberStr;
		mc_maptype		+= NumberStr;
		mc_unitconvfact += NumberStr;
		mc_processed	+= NumberStr;
		mc_file			+= NumberStr;

		mc_freq			+= std::wstring(MC_FILEFREQSTR);
		mc_units		+= std::wstring(MC_FILEUNITSSTR);
		mc_coordsys		+= std::wstring(MC_FILECOORDSYSSTR);
		mc_maptype		+= std::wstring(MC_FILEMAPTYPESTR);
		mc_unitconvfact += std::wstring(MC_FILEUNITCONVFACTORSTR);
		mc_processed	+= std::wstring(MC_FILEPROCESSEDSTR);
		mc_file			+= std::wstring(MC_FILEMAPFILENAMESTR);

		Insert1element(mc_freq,DBField(Params_->find<int>(Wstr2Str(mc_freq))));	
		Insert1element(mc_units,DBField(Params_->find<int>(Wstr2Str(mc_units))));	
		Insert1element(mc_coordsys,DBField(Params_->find<int>(Wstr2Str(mc_coordsys))));	
		Insert1element(mc_unitconvfact,DBField(Params_->find<double>(Wstr2Str(mc_unitconvfact))));	
		Insert1element(mc_processed,DBField(Params_->find<int>(Wstr2Str(mc_processed))));	
		Insert1element(mc_maptype,DBField(Params_->find<int>(Wstr2Str(mc_maptype))));
		std::string	mc_filestr(Wstr2Str(mc_file));
		Insert1element(mc_file,DBField(Achar2Wstr((Params_->find<std::string>(mc_filestr)).c_str())));
	}

	for(int i=0;i<NFreq;++i)
	{
		ID			= std::wstring(GLOBAL_OBSPLANEFREQ);
		IDFWHM		= std::wstring(GLOBAL_OBSPLANEFWHM);
		Zeus::PutNumber2Txt(buffer,i);
		ID			+= std::wstring(buffer);
		IDFWHM		+= std::wstring(buffer);
		Insert1element(ID,DBField(Params_->find<int>(Wstr2Str(ID))));	
		Insert1element(IDFWHM,DBField(Params_->find<int>(Wstr2Str(IDFWHM))));	
	}

	for(int i=0;i!=NPrior;++i)
	{
		ID	= std::wstring(GLOBAL_PRIORMAPTYPE);
		ID1	= std::wstring(GLOBAL_PRIORMAPFREQ);
		ID2 = std::wstring(GLOBAL_PRIORMAXL);
		ID3 = std::wstring(GLOBAL_PRIORNPNOISERMS);

		Zeus::PutNumber2Txt(buffer,i);
		StrTemp = std::wstring(buffer);
		ID  += StrTemp;
		ID1 += StrTemp;
		ID2 += StrTemp;
		ID3 += StrTemp;
		Insert1element(ID,DBField(Params_->find<int>(Wstr2Str(ID))));	
		Insert1element(ID1,DBField(Params_->find<int>(Wstr2Str(ID1))));	
		Insert1element(ID2,DBField(Params_->find<double>(Wstr2Str(ID2))));	
		Insert1element(ID3,DBField(Params_->find<double>(Wstr2Str(ID3))));	
	}
}
//
void	ReadPatchesGeomInfoLFIDPC::ReadHeader(void)
{
	inp_->getKey(std::string(HFIDMC_GEOMHEADERINFO_NSIDE),data_.Header_.NSide_);
	inp_->getKey(std::string(HFIDMC_GEOMHEADERINFO_PTCHSZ),data_.Header_.PtchSz_);
	inp_->getKey(std::string(HFIDMC_GEOMHEADERINFO_PTCHBORDER),data_.Header_.PtchBorder_);
	inp_->getKey(std::string(HFIDMC_GEOMHEADERINFO_NPCTMPIX),data_.Header_.NPatchesPerMainPix_);
	inp_->getKey(std::string(HFIDMC_GEOMHEADERINFO_NTOTPTCH),data_.Header_.NTotalPatches_);
	inp_->getKey(std::string(HFIDMC_GEOMHEADERINFO_GALCUT),data_.Header_.GalacticCut_);
}
//
void	ReadPatchesGeomInfoLFIDPC::ReadBody(void)
{
	arr<int>	tempArrINT[GEOMPROPNINTS];
	arr<double> tempArrFLOAT[GEOMPROPNFLOATS];

	inp_->readEntireColumn(std::string(HFIDMC_GEOMLINESINFO_PATCHN),tempArrINT[0]);
	inp_->readEntireColumn(std::string(HFIDMC_GEOMLINESINFO_BASEPIX),tempArrINT[1]);
	inp_->readEntireColumn(std::string(HFIDMC_GEOMLINESINFO_SPIN),tempArrFLOAT[0]);
	inp_->readEntireColumn(std::string(HFIDMC_GEOMLINESINFO_CCCOL),tempArrFLOAT[1]);
	inp_->readEntireColumn(std::string(HFIDMC_GEOMLINESINFO_CCLON),tempArrFLOAT[2]);
	inp_->readEntireColumn(std::string(HFIDMC_GEOMLINESINFO_LDCOL),tempArrFLOAT[3]);
	inp_->readEntireColumn(std::string(HFIDMC_GEOMLINESINFO_LDLON),tempArrFLOAT[4]);
	inp_->readEntireColumn(std::string(HFIDMC_GEOMLINESINFO_RDCOL),tempArrFLOAT[5]);
	inp_->readEntireColumn(std::string(HFIDMC_GEOMLINESINFO_RDLON),tempArrFLOAT[6]);
	inp_->readEntireColumn(std::string(HFIDMC_GEOMLINESINFO_LUCOL),tempArrFLOAT[7]);
	inp_->readEntireColumn(std::string(HFIDMC_GEOMLINESINFO_LULON),tempArrFLOAT[8]);
	inp_->readEntireColumn(std::string(HFIDMC_GEOMLINESINFO_RUCOL),tempArrFLOAT[9]);
	inp_->readEntireColumn(std::string(HFIDMC_GEOMLINESINFO_RULON),tempArrFLOAT[10]);
	inp_->readEntireColumn(std::string(HFIDMC_GEOMLINESINFO_RRINIT),tempArrFLOAT[11]);
	inp_->readEntireColumn(std::string(HFIDMC_GEOMLINESINFO_ARRFIN),tempArrFLOAT[12]);
	inp_->readEntireColumn(std::string(HFIDMC_GEOMLINESINFO_ROTDIT),tempArrINT[2]);
	inp_->readEntireColumn(std::string(HFIDMC_GEOMLINESINFO_PREDFLUX),tempArrFLOAT[13]);
	inp_->readEntireColumn(std::string(HFIDMC_GEOMLINESINFO_PREDRAD),tempArrFLOAT[14]);
	inp_->readEntireColumn(std::string(HFIDMC_GEOMLINESINFO_SNR),tempArrFLOAT[15]);
	inp_->readEntireColumn(std::string(HFIDMC_GEOMLINESINFO_SRCID),tempArrINT[3]);

	PatchGeomLineType	tempCatalogAtom;

	const int	NElements(static_cast<int>(tempArrINT[0].end() - tempArrINT[0].begin()));

	for(int i=0;i<NElements;++i)
	{		
		tempCatalogAtom.PatchNumber_		= (tempArrINT[0])[i];
		tempCatalogAtom.BPixel_				= (tempArrINT[1])[i];
		tempCatalogAtom.Spin_				= static_cast<double>((tempArrFLOAT[0])[i]);
		tempCatalogAtom.X0Y0Ptg_.theta		= static_cast<double>((tempArrFLOAT[1])[i]);
		tempCatalogAtom.X0Y0Ptg_.phi		= static_cast<double>((tempArrFLOAT[2])[i]);
		tempCatalogAtom.X0Y0_.theta			= static_cast<double>((tempArrFLOAT[3])[i]);
		tempCatalogAtom.X0Y0_.phi			= static_cast<double>((tempArrFLOAT[4])[i]);
		tempCatalogAtom.XLY0_.theta			= static_cast<double>((tempArrFLOAT[5])[i]);
		tempCatalogAtom.XLY0_.phi			= static_cast<double>((tempArrFLOAT[6])[i]);
		tempCatalogAtom.X0YL_.theta			= static_cast<double>((tempArrFLOAT[7])[i]);
		tempCatalogAtom.X0YL_.phi			= static_cast<double>((tempArrFLOAT[8])[i]);
		tempCatalogAtom.XLYL_.theta			= static_cast<double>((tempArrFLOAT[9])[i]);
		tempCatalogAtom.XLYL_.phi			= static_cast<double>((tempArrFLOAT[10])[i]);
		tempCatalogAtom.InitRejectRatio_	= static_cast<double>((tempArrFLOAT[11])[i]);
		tempCatalogAtom.FinalRejectRatio_	= static_cast<double>((tempArrFLOAT[12])[i]);
		tempCatalogAtom.RotationDir_		= (tempArrINT[2])[i];
		tempCatalogAtom.PredFlux_			= static_cast<double>((tempArrFLOAT[13])[i]);
		tempCatalogAtom.PredRadius_			= static_cast<double>((tempArrFLOAT[14])[i]);
		tempCatalogAtom.SNR_				= static_cast<double>((tempArrFLOAT[15])[i]);
		tempCatalogAtom.SrcIndex_			= (tempArrINT[3])[i];

		data_.Storage_.push_back(tempCatalogAtom);
	}
}
//
void	ReadNonBlindCatalogueLFIDPC::ReadHeader(void)
{
	int	DummyINT;
	CollHandle_->getKey(std::string(HFIDMC_NONBLINDHEADER_DEGREES),data_.Header_.CoordsType_);
	CollHandle_->getKey(std::string(HFIDMC_NONBLINDHEADER_COORDSYS),DummyINT);
	CollHandle_->getKey(std::string(HFIDMC_NONBLINDHEADER_EPOCH),data_.Header_.Epoch_);
	data_.Header_.CoordSystem_	= static_cast<coordsys>(DummyINT);
}
//
int		CatalogueWriterLFIDPCSZ::GetNextColumn(const CatLineCollType & DataIN,
		const SZPS_ProfParamType& ProfParam,arr<double>& data,std::string& ColumnName)
{
	if(CurrColumn_ == 31)
		return 0;

	Trafo	TransCoordsEqu(2000.0,2000.0,Galactic,Equatorial);
	// J2000 coordinate system, EPOCH = 2000
	pointing		tempEqu;

	double								*temp(data.begin());
	CatLineCollType::const_iterator		pivIN(DataIN.begin());
	CatLineCollType::const_iterator		const endIN(DataIN.end());

	for(;pivIN != endIN;++pivIN,++temp)
	{
		switch(CurrColumn_)
		{
		case 0:
			*temp = pivIN->GalLongRads_ * RAD2DEGREE;
			break;
		case 1:
			*temp = (PIOVER2 - pivIN->GalCoLatRads_) * RAD2DEGREE;;
			break;
		case 2:
			tempEqu = TransCoordsEqu(pointing(pivIN->GalCoLatRads_,pivIN->GalLongRads_));
			*temp = tempEqu.phi * RAD2DEGREE;
			break;
		case 3:
			tempEqu = TransCoordsEqu(pointing(pivIN->GalCoLatRads_,pivIN->GalLongRads_));
			*temp = (PIOVER2 - tempEqu.theta) * RAD2DEGREE;
			break;
		case 4:
			*temp = pivIN->ErrorBars_.TotalPosErrorBar_;
			break;
		case 5:
			*temp = pivIN->FluxComptGLRT_;
			break;
		case 6:
			*temp = pivIN->RadiusGLRT_;
			break;
		case 7:
			*temp = ((pivIN->NormalAmpl_ < pivIN->DetectSigma_)?pivIN->NormalAmpl_:pivIN->DetectSigma_);
			break;
		case 8:
			*temp = pivIN->FluxComptGLRT_ / ProfParam.MNFW_Ratio_CY500CYR500_;
			break;
		case 9:
			*temp = pivIN->RadiusGLRT_ * ProfParam.MNFW_C500_;
			break;
		case 10:
			*temp = pivIN->FluxCompt_;
			break;
		case 11:
			*temp = pivIN->Radius_;
			break;
		case 12:
			*temp = pivIN->FluxCompt_ / ProfParam.MNFW_Ratio_CY500CYR500_;
			break;
		case 13:
			*temp = pivIN->Radius_ * ProfParam.MNFW_C500_;
			break;
		case 14:
			*temp = pivIN->ErrorBars_.FluxErrorBar_;
			break;
		case 15:
			*temp = pivIN->ErrorBars_.RadiusErrorBar_;
			break;
		case 16:
			*temp = pivIN->ErrorBars_.FluxErrorBar_ / ProfParam.MNFW_Ratio_CY500CYR500_;
			break;
		case 17:
			*temp = pivIN->ErrorBars_.RadiusErrorBar_ * ProfParam.MNFW_C500_;
			break;
		default:
			*temp = SZCATALOGUEDEFAULTVALUE;
			break;
		}
	}
	ColumnName = Zeus::Wstr2Str(std::wstring(L"FIELDNO_") + PutNumber2Txt(CurrColumn_));
	return ++CurrColumn_;
}

int		CatalogueWriterLFIDPCPS::GetNextColumn(const CatLineCollType & DataIN,
		const SZPS_ProfParamType& ProfParam,arr<double>& data,std::string& ColumnName)
{
	if(CurrColumn_ == 31)
		return 0;

	Trafo	TransCoordsEqu(2000.0,2000.0,Galactic,Equatorial);
	// J2000 coordinate system, EPOCH = 2000
	pointing		tempEqu;

	double								*temp(data.begin());
	CatLineCollType::const_iterator		pivIN(DataIN.begin());
	CatLineCollType::const_iterator		const endIN(DataIN.end());

	for(;pivIN != endIN;++pivIN,++temp)
	{
		switch(CurrColumn_)
		{
		case 0:
			*temp = pivIN->GalLongRads_ * RAD2DEGREE;
			break;
		case 1:
			*temp = (PIOVER2 - pivIN->GalCoLatRads_) * RAD2DEGREE;
			break;
		case 2:
			tempEqu = TransCoordsEqu(pointing(pivIN->GalCoLatRads_,pivIN->GalLongRads_));
			*temp = tempEqu.phi * RAD2DEGREE;
			break;
		case 3:
			tempEqu = TransCoordsEqu(pointing(pivIN->GalCoLatRads_,pivIN->GalLongRads_));
			*temp = (PIOVER2 - tempEqu.theta) * RAD2DEGREE;
			break;
		case 4:
			*temp = pivIN->ErrorBars_.TotalPosErrorBar_;
			break;
		case 5:
			*temp = pivIN->FluxComptGLRT_;
			break;
		case 6:
			*temp = pivIN->RadiusGLRT_;
			break;
		case 7:
			*temp = ((pivIN->NormalAmpl_ < pivIN->DetectSigma_)?pivIN->NormalAmpl_:pivIN->DetectSigma_);
			break;
		case 8:
			*temp = SZCATALOGUEDEFAULTVALUE;
			break;
		case 9:
			*temp = SZCATALOGUEDEFAULTVALUE;
			break;
		case 10:
			*temp = pivIN->FluxCompt_;
			break;
		case 11:
			*temp = pivIN->Radius_;
			break;
		case 12:
			*temp = SZCATALOGUEDEFAULTVALUE;
			break;
		case 13:
			*temp = SZCATALOGUEDEFAULTVALUE;
			break;
		case 14:
			*temp = pivIN->ErrorBars_.FluxErrorBar_;
			break;
		case 15:
			*temp = pivIN->ErrorBars_.RadiusErrorBar_;
			break;
		case 16:
			*temp = SZCATALOGUEDEFAULTVALUE;
			break;
		case 17:
			*temp = SZCATALOGUEDEFAULTVALUE;
			break;
		default:
			*temp = SZCATALOGUEDEFAULTVALUE;
			break;
		}
	}
	ColumnName = Zeus::Wstr2Str(std::wstring(LFIDPC_FIELDNAMEPREFW) + PutNumber2Txt(CurrColumn_));
	return ++CurrColumn_;
}

//
void	ReadNonBlindCatalogueLFIDPC::ReadBody(void)
{
	arr<int>				tArrINT;
	arr<double>				tArrFLOAT;
	int						*tArrINTptr;
	double					*tArrFLOATptr;

	data_.Storage_.clear();
	data_.Storage_.resize(CollHandle_->columnLength(CollHandle_->columnNumber(std::string(LFIDPC_FIELDNAME("0000")))));

	NonBlingCatType::StorageType::iterator			const pivOrg(data_.Storage_.begin());
	NonBlingCatType::StorageType::iterator			piv;
	NonBlingCatType::StorageType::const_iterator	const End(data_.Storage_.end());

	INTFILEREADVECTFLOATLFI(LFIDPC_FIELDNAME("0000"),ptg_.phi);
	INTFILEREADVECTFLOATLFI(LFIDPC_FIELDNAME("0001"),ptg_.theta);
	INTFILEREADVECTFLOATLFI(LFIDPC_FIELDNAME("0004"),ErrPosition_);
	INTFILEREADVECTFLOATLFI(LFIDPC_FIELDNAME("0005"),PredFluxGLRT_);
	INTFILEREADVECTFLOATLFI(LFIDPC_FIELDNAME("0006"),PredRadiusGLRT_);
	INTFILEREADVECTFLOATLFI(LFIDPC_FIELDNAME("0007"),SNR_);
	INTFILEREADVECTFLOATLFI(LFIDPC_FIELDNAME("0010"),PredFluxBay_);
	INTFILEREADVECTFLOATLFI(LFIDPC_FIELDNAME("0011"),PredRadiusBay_);
	INTFILEREADVECTFLOATLFI(LFIDPC_FIELDNAME("0014"),ErrFlux_);
	INTFILEREADVECTFLOATLFI(LFIDPC_FIELDNAME("0015"),ErrRadius_);

	int i(0);
	for(piv = pivOrg;piv != End;++piv,++i)
	{
		piv->Index_			= i;
		piv->ptg_.theta		= (90.0 - piv->ptg_.theta) / RAD2DEGREE;
		piv->ptg_.phi		/= RAD2DEGREE;
	}
}
//
void	WritePatchesGeomInfoLFIDPC::WriteGeomInfoData(const PatchGeomType & data)
{
	CreatHeader(data);

	if(data.Storage_.empty())	return;

	arr<double>						VectDOUBLE(data.Storage_.size());
	arr<int>						VectINT(data.Storage_.size());
	double							*VectDOUBLEptr;
	int								*VectINTptr;

	PatchGeomType::StorageType::const_iterator	piv;
	PatchGeomType::StorageType::const_iterator	const end;

	INTFILEWRITEVECTINTLFI(HFIDMC_GEOMLINESINFO_PATCHN,PatchNumber_);
	INTFILEWRITEVECTINTLFI(HFIDMC_GEOMLINESINFO_BASEPIX,BPixel_);
	INTFILEWRITEVECTDOUBLELFI(HFIDMC_GEOMLINESINFO_SPIN,Spin_);
	INTFILEWRITEVECTDOUBLELFI(HFIDMC_GEOMLINESINFO_CCCOL,X0Y0Ptg_.theta);
	INTFILEWRITEVECTDOUBLELFI(HFIDMC_GEOMLINESINFO_CCLON,X0Y0Ptg_.phi);
	INTFILEWRITEVECTDOUBLELFI(HFIDMC_GEOMLINESINFO_LDCOL,X0Y0_.theta);
	INTFILEWRITEVECTDOUBLELFI(HFIDMC_GEOMLINESINFO_LDLON,X0Y0_.phi);
	INTFILEWRITEVECTDOUBLELFI(HFIDMC_GEOMLINESINFO_RDCOL,XLY0_.theta);
	INTFILEWRITEVECTDOUBLELFI(HFIDMC_GEOMLINESINFO_RDLON,XLY0_.phi);
	INTFILEWRITEVECTDOUBLELFI(HFIDMC_GEOMLINESINFO_LUCOL,X0YL_.theta);
	INTFILEWRITEVECTDOUBLELFI(HFIDMC_GEOMLINESINFO_LULON,X0YL_.phi);
	INTFILEWRITEVECTDOUBLELFI(HFIDMC_GEOMLINESINFO_RUCOL,XLYL_.theta);
	INTFILEWRITEVECTDOUBLELFI(HFIDMC_GEOMLINESINFO_RULON,XLYL_.phi);
	INTFILEWRITEVECTDOUBLELFI(HFIDMC_GEOMLINESINFO_RRINIT,InitRejectRatio_);
	INTFILEWRITEVECTDOUBLELFI(HFIDMC_GEOMLINESINFO_ARRFIN,FinalRejectRatio_);
	INTFILEWRITEVECTINTLFI(HFIDMC_GEOMLINESINFO_ROTDIT,RotationDir_);
	INTFILEWRITEVECTDOUBLELFI(HFIDMC_GEOMLINESINFO_PREDFLUX,PredFlux_);
	INTFILEWRITEVECTDOUBLELFI(HFIDMC_GEOMLINESINFO_PREDRAD,PredRadius_);
	INTFILEWRITEVECTDOUBLELFI(HFIDMC_GEOMLINESINFO_SNR,SNR_);
	INTFILEWRITEVECTINTLFI(HFIDMC_GEOMLINESINFO_SRCID,SrcIndex_);
}
//
void	WritePatchesGeomInfoLFIDPC::CreatHeader(const PatchGeomType & data)
{
	CollHandle_->setKey(std::string(HFIDMC_GEOMHEADERINFO_NSIDE),data.Header_.NSide_);
	CollHandle_->setKey(std::string(HFIDMC_GEOMHEADERINFO_PTCHSZ),data.Header_.PtchSz_);
	CollHandle_->setKey(std::string(HFIDMC_GEOMHEADERINFO_PTCHBORDER),data.Header_.PtchBorder_);
	CollHandle_->setKey(std::string(HFIDMC_GEOMHEADERINFO_NPCTMPIX),data.Header_.NPatchesPerMainPix_);
	CollHandle_->setKey(std::string(HFIDMC_GEOMHEADERINFO_NTOTPTCH),data.Header_.NTotalPatches_);
	CollHandle_->setKey(std::string(HFIDMC_GEOMHEADERINFO_GALCUT),data.Header_.GalacticCut_);
}
//
int		WriteObjResultsLFIDPC::WritePeaksColl(const PeakCollType& data)
{

	if(data.empty())	return ObjectID_;

	arr<double>						VectDOUBLE(data.size());
	arr<int>						VectINT(data.size());
	double							*VectDOUBLEptr;
	int								*VectINTptr;
	PeakCollType::const_iterator	piv(data.begin());
	PeakCollType::const_iterator	const end(data.end());


	for(VectINTptr = VectINT.begin();piv!=end;++piv,++VectINTptr,++ObjectID_)
	{*VectINTptr = ObjectID_;}
	CollHandle_->appendColumn(std::string(INTFILE_FLDNAME_OBJID),VectINT);

	INTFILEAPPENDVECTINTLFI(INTFILE_FLDNAME_PATCHN,PatchNumber_);
	INTFILEAPPENDVECTINTLFI(INTFILE_FLDNAME_BAYSTAT,PK_BayesDetectStat_);
	INTFILEAPPENDVECTDOUBLELFI(INTFILE_FLDNAME_GAUSLEV,GaussianIndex_);
	INTFILEAPPENDVECTDOUBLELFI(INTFILE_FLDNAME_DETSIGMA,DetectionSigma_);
	INTFILEAPPENDVECTDOUBLELFI(INTFILE_FLDNAME_AMPLNORM,SrcAmplNormalised_);
	INTFILEAPPENDVECTDOUBLELFI(INTFILE_FLDNAME_SIGMTHRES,UsedSigmaThreshold_);
	INTFILEAPPENDVECTDOUBLELFI(INTFILE_FLDNAME_BAYSIGMATH,JF_UsedSigmaThreshold_);
	INTFILEAPPENDVECTDOUBLELFI(INTFILE_FLDNAME_FLUX,SrcFlux_);
	INTFILEAPPENDVECTDOUBLELFI(INTFILE_FLDNAME_BAYFLUXMAX,JF_SrcFlux_.Mode_);
	INTFILEAPPENDVECTDOUBLELFI(INTFILE_FLDNAME_BAYLNRHOTH,JF_lnRhoTh_);
	INTFILEAPPENDVECTDOUBLELFI(INTFILE_FLDNAME_BAYLNRHO,JF_lnRho_);
	INTFILEAPPENDVECTDOUBLELFI(INTFILE_FLDNAME_BAYLNEVID,JF_lnEvidence_);
	INTFILEAPPENDVECTDOUBLELFI(INTFILE_FLDNAME_BAYLNEVIDERR,JF_lnEvidenceErrBar_);
	INTFILEAPPENDVECTDOUBLELFI(INTFILE_FLDNAME_BAYLNFORMFACT,JF_lnFormFactor_);
	INTFILEAPPENDVECTDOUBLELFI(INTFILE_FLDNAME_BAYLNMODRATIO,JF_lnModelRatio_);
	INTFILEAPPENDVECTDOUBLELFI(INTFILE_FLDNAME_POSYCOORD,Pos_.YCoord_);
	INTFILEAPPENDVECTDOUBLELFI(INTFILE_FLDNAME_POSXCOORD,Pos_.XCoord_);
	INTFILEAPPENDVECTDOUBLELFI(INTFILE_FLDNAME_POSBAYYCMAX,Pos_.JF_YCoord_.Mode_);
	INTFILEAPPENDVECTDOUBLELFI(INTFILE_FLDNAME_POSBAYXCMAX,Pos_.JF_XCoord_.Mode_);
	INTFILEAPPENDVECTDOUBLELFI(INTFILE_FLDNAME_POSBAYYCMEAN,Pos_.JF_YCoord_.Mean_);
	INTFILEAPPENDVECTDOUBLELFI(INTFILE_FLDNAME_POSBAYXCMEAN,Pos_.JF_XCoord_.Mean_);
	INTFILEAPPENDVECTDOUBLELFI(INTFILE_FLDNAME_POSGALCOLAT,GalPt_Colat_);
	INTFILEAPPENDVECTDOUBLELFI(INTFILE_FLDNAME_POSGALLONG,GalPt_Long_);
	INTFILEAPPENDVECTDOUBLELFI(INTFILE_FLDNAME_FLUXMJYS,SrcFlux_mJys_);
	INTFILEAPPENDVECTDOUBLELFI(INTFILE_FLDNAME_COMPARC2,SrcCompt_arcmin2_);
	INTFILEAPPENDVECTDOUBLELFI(INTFILE_FLDNAME_RADIUS,RealParams_.RealScale_);
	INTFILEAPPENDVECTDOUBLELFI(INTFILE_FLDNAME_BAYRADMAX,JF_Radius_.Mode_);
	INTFILEAPPENDVECTDOUBLELFI(INTFILE_FLDNAME_ODDS,Odds_);
	INTFILEAPPENDVECTDOUBLELFI(INTFILE_FLDNAME_ISNR2,ISNR2_);
	INTFILEAPPENDVECTDOUBLELFI(INTFILE_FLDNAME_INITODDS,CollListIndex_);
	INTFILEAPPENDVECTDOUBLELFI(INTFILE_FLDNAME_FLUXERR,ErrorBars_.FluxErrorBar_);
	INTFILEAPPENDVECTDOUBLELFI(INTFILE_FLDNAME_POSERR,ErrorBars_.TotalPosErrorBar_);
	INTFILEAPPENDVECTDOUBLELFI(INTFILE_FLDNAME_RADERR,ErrorBars_.RadiusErrorBar_);
	INTFILEAPPENDVECTDOUBLELFI(INTFILE_FLDNAME_FLUXCONV,y0_);
	INTFILEAPPENDVECTDOUBLELFI(INTFILE_FLDNAME_BAYFLUXMEAN,JF_SrcFlux_.Mean_);
	INTFILEAPPENDVECTDOUBLELFI(INTFILE_FLDNAME_BAYRADMEAN,JF_Radius_.Mean_);
	INTFILEAPPENDVECTDOUBLELFI(INTFILE_FLDNAME_POSECLLAT,CoordPt_Lat_);
	INTFILEAPPENDVECTDOUBLELFI(INTFILE_FLDNAME_POSECLLON,CoordPt_Long_);
	INTFILEAPPENDVECTDOUBLELFI(INTFILE_FLDNAME_POSEQULAT,EquPt_Lat_);
	INTFILEAPPENDVECTDOUBLELFI(INTFILE_FLDNAME_POSEQULON,EquPt_Long_);

	return ObjectID_;
}
//
void	ReadObjResultsLFIDPC::ReadBody(void)
{

	arr<int>				tArrINT;
	arr<double>				tArrFLOAT;
	int						*tArrINTptr;
	double					*tArrFLOATptr;

	data_.Storage_.clear();
	data_.Storage_.resize(CollHandle_->columnLength(CollHandle_->columnNumber(std::string(INTFILE_FLDNAME_OBJID))));

	PeakCollReadbleType::StorageType::iterator			const pivOrg(data_.Storage_.begin());
	PeakCollReadbleType::StorageType::iterator			piv;
	PeakCollReadbleType::StorageType::const_iterator	const End(data_.Storage_.end());


	INTFILEREADVECTINTLFI(INTFILE_FLDNAME_OBJID,DetectID_);
	INTFILEREADVECTINTLFI(INTFILE_FLDNAME_PATCHN,PatchNumber_);
	INTFILEREADVECTINTLFI(INTFILE_FLDNAME_BAYSTAT,PK_BayesDetectStat_);
	INTFILEREADVECTFLOATLFI(INTFILE_FLDNAME_GAUSLEV,GaussianIndex_);
	INTFILEREADVECTFLOATLFI(INTFILE_FLDNAME_DETSIGMA,DetectionSigma_);
	INTFILEREADVECTFLOATLFI(INTFILE_FLDNAME_AMPLNORM,SrcAmplNormalised_);
	INTFILEREADVECTFLOATLFI(INTFILE_FLDNAME_SIGMTHRES,UsedSigmaThreshold_);
	INTFILEREADVECTFLOATLFI(INTFILE_FLDNAME_BAYSIGMATH,JF_UsedSigmaThreshold_);
	INTFILEREADVECTFLOATLFI(INTFILE_FLDNAME_FLUX,SrcFlux_);
	INTFILEREADVECTFLOATLFI(INTFILE_FLDNAME_BAYFLUXMAX,JF_SrcFlux_.Mode_);
	INTFILEREADVECTFLOATLFI(INTFILE_FLDNAME_BAYLNRHOTH,JF_lnRhoTh_);
	INTFILEREADVECTFLOATLFI(INTFILE_FLDNAME_BAYLNRHO,JF_lnRho_);
	INTFILEREADVECTFLOATLFI(INTFILE_FLDNAME_BAYLNEVID,JF_lnEvidence_);
	INTFILEREADVECTFLOATLFI(INTFILE_FLDNAME_BAYLNEVIDERR,JF_lnEvidenceErrBar_);
	INTFILEREADVECTFLOATLFI(INTFILE_FLDNAME_BAYLNFORMFACT,JF_lnFormFactor_);
	INTFILEREADVECTFLOATLFI(INTFILE_FLDNAME_BAYLNMODRATIO,JF_lnModelRatio_);
	INTFILEREADVECTFLOATLFI(INTFILE_FLDNAME_POSYCOORD,Pos_.YCoord_);
	INTFILEREADVECTFLOATLFI(INTFILE_FLDNAME_POSXCOORD,Pos_.XCoord_);
	INTFILEREADVECTFLOATLFI(INTFILE_FLDNAME_POSBAYYCMAX,Pos_.JF_YCoord_.Mode_);
	INTFILEREADVECTFLOATLFI(INTFILE_FLDNAME_POSBAYXCMAX,Pos_.JF_XCoord_.Mode_);
	INTFILEREADVECTFLOATLFI(INTFILE_FLDNAME_POSBAYYCMEAN,Pos_.JF_YCoord_.Mean_);
	INTFILEREADVECTFLOATLFI(INTFILE_FLDNAME_POSBAYXCMEAN,Pos_.JF_XCoord_.Mean_);
	INTFILEREADVECTFLOATLFI(INTFILE_FLDNAME_POSGALCOLAT,GalPt_Colat_);
	INTFILEREADVECTFLOATLFI(INTFILE_FLDNAME_POSGALLONG,GalPt_Long_);
	INTFILEREADVECTFLOATLFI(INTFILE_FLDNAME_FLUXMJYS,SrcFlux_mJys_);
	INTFILEREADVECTFLOATLFI(INTFILE_FLDNAME_COMPARC2,SrcCompt_arcmin2_);
	INTFILEREADVECTFLOATLFI(INTFILE_FLDNAME_RADIUS,RealParams_.RealScale_);
	INTFILEREADVECTFLOATLFI(INTFILE_FLDNAME_BAYRADMAX,JF_Radius_.Mode_);
	INTFILEREADVECTFLOATLFI(INTFILE_FLDNAME_ODDS,Odds_);
	INTFILEREADVECTFLOATLFI(INTFILE_FLDNAME_ISNR2,ISNR2_);
	INTFILEREADVECTFLOATLFI(INTFILE_FLDNAME_INITODDS,CollListIndex_);
	INTFILEREADVECTFLOATLFI(INTFILE_FLDNAME_FLUXERR,ErrorBars_.FluxErrorBar_);
	INTFILEREADVECTFLOATLFI(INTFILE_FLDNAME_POSERR,ErrorBars_.TotalPosErrorBar_);
	INTFILEREADVECTFLOATLFI(INTFILE_FLDNAME_RADERR,ErrorBars_.RadiusErrorBar_);
	INTFILEREADVECTFLOATLFI(INTFILE_FLDNAME_FLUXCONV,y0_);
	INTFILEREADVECTFLOATLFI(INTFILE_FLDNAME_BAYFLUXMEAN,JF_SrcFlux_.Mean_);
	INTFILEREADVECTFLOATLFI(INTFILE_FLDNAME_BAYRADMEAN,JF_Radius_.Mean_);
	INTFILEREADVECTFLOATLFI(INTFILE_FLDNAME_POSECLLAT,CoordPt_Lat_);
	INTFILEREADVECTFLOATLFI(INTFILE_FLDNAME_POSECLLON,CoordPt_Long_);
	INTFILEREADVECTFLOATLFI(INTFILE_FLDNAME_POSEQULAT,EquPt_Lat_);
	INTFILEREADVECTFLOATLFI(INTFILE_FLDNAME_POSEQULON,EquPt_Long_);
}
//
#endif //LFIDPC

} // namespace Zeus


