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

#include "ZEUS_InOutBinFile.h"

//------------------------------------

namespace Zeus
{
//
void		CatWriterOutputBinFile_QA_Contours::AddHeaderPCCFields(const OutputFormatType& data)
{
	char			buffer[BUFFERMAXCHAR];
	std::string		StrBuffer;

	std::vector<fitscolumn> cols;

	cols.push_back (fitscolumn("GLAT_inj","RAD",1,planckType<double>()));
	cols.push_back (fitscolumn("GLON_inj","RAD",1,planckType<double>()));
	cols.push_back (fitscolumn("THETA_S","ARCMIN",1,planckType<double>()));
	cols.push_back (fitscolumn("Y_inj","ARCMIN^2",1,planckType<double>()));
	cols.push_back (fitscolumn("HPD_2D","UNITLESS",1,planckType<double>()));
	cols.push_back (fitscolumn("HPD_Y","UNITLESS",1,planckType<double>()));
	cols.push_back (fitscolumn("HPD_THETA","UNITLESS",1,planckType<double>()));
	cols.push_back (fitscolumn("HPD_YTHETA","UNITLESS",1,planckType<double>()));
	cols.push_back (fitscolumn("Ypeak","ARCMIN^2",1,planckType<double>()));
	cols.push_back (fitscolumn("Y_68_low","ARCMIN^2",1,planckType<double>()));
	cols.push_back (fitscolumn("Y_68_up","ARCMIN^2",1,planckType<double>()));
	cols.push_back (fitscolumn("Y_95_low","ARCMIN^2",1,planckType<double>()));
	cols.push_back (fitscolumn("Y_95_up","ARCMIN^2",1,planckType<double>()));
	cols.push_back (fitscolumn("ThetaS_peak","ARCMIN",1,planckType<double>()));
	cols.push_back (fitscolumn("ThetaS_68_low","ARCMIN",1,planckType<double>()));
	cols.push_back (fitscolumn("ThetaS_68_up","ARCMIN",1,planckType<double>()));
	cols.push_back (fitscolumn("ThetaS_95_low","ARCMIN",1,planckType<double>()));
	cols.push_back (fitscolumn("ThetaS_95_up","ARCMIN",1,planckType<double>()));
	cols.push_back (fitscolumn("YTheta_peak","ARCMIN^2",1,planckType<double>()));
	cols.push_back (fitscolumn("YTheta_68_low","ARCMIN^2",1,planckType<double>()));
	cols.push_back (fitscolumn("YTheta_68_up","ARCMIN^2",1,planckType<double>()));
	cols.push_back (fitscolumn("YTheta_95_low","ARCMIN^2",1,planckType<double>()));
	cols.push_back (fitscolumn("YTheta_95_up","ARCMIN^2",1,planckType<double>()));
	cols.push_back (fitscolumn("Ylow_axis","UNITLESS",1,planckType<double>()));
	cols.push_back (fitscolumn("Yhigh_axis","UNITLESS",1,planckType<double>()));
	cols.push_back (fitscolumn("ThetaS_low_axis","UNITLESS",1,planckType<double>()));
	cols.push_back (fitscolumn("ThetaS_high_axis","UNITLESS",1,planckType<double>()));

	FitsHandle_.insert_bintab(cols);
}
//
void		CatWriterOutputBinFile_PCCFITS::AddHeaderPCCFields(const OutputFormatType& data)
{
	char			buffer[BUFFERMAXCHAR];
	std::string		StrBuffer;

	std::vector<fitscolumn> cols;
#ifdef WG5ALLDOUBLES
	cols.push_back (fitscolumn("GLON","DEG",1,planckType<double>()));
	cols.push_back (fitscolumn("GLAT","DEG",1,planckType<double>()));
	cols.push_back (fitscolumn("RAJ2000","DEG",1,planckType<double>()));
	cols.push_back (fitscolumn("DEJ2000","DEG",1,planckType<double>()));
	cols.push_back (fitscolumn("ERR_RAD","ARCMIN",1,planckType<double>()));
	cols.push_back (fitscolumn("CY5R500","ARCMIN^2",1,planckType<double>()));
	cols.push_back (fitscolumn("TS","ARCMIN",1,planckType<double>()));
	cols.push_back (fitscolumn("SNR","UNITLESS",1,planckType<double>()));
	cols.push_back (fitscolumn("CY500","ARCMIN^2",1,planckType<double>()));
	cols.push_back (fitscolumn("T500","ARCMIN",1,planckType<double>()));
	cols.push_back (fitscolumn("CY5R500F","ARCMIN^2",1,planckType<double>()));
	cols.push_back (fitscolumn("TSF","ARCMIN",1,planckType<double>()));
	cols.push_back (fitscolumn("CY500F","ARCMIN^2",1,planckType<double>()));
	cols.push_back (fitscolumn("T500F","ARCMIN",1,planckType<double>()));
	cols.push_back (fitscolumn("ERR_CY5R500F","ARCMIN^2",1,planckType<double>()));
	cols.push_back (fitscolumn("ERR_TSF","ARCMIN",1,planckType<double>()));
	cols.push_back (fitscolumn("ERR_CY500F","ARCMIN^2",1,planckType<double>()));
	cols.push_back (fitscolumn("ERR_T500F","ARCMIN",1,planckType<double>()));
	cols.push_back (fitscolumn("FLUX5R500_100GHZ","DTCMB/TCMB ARCMIN^2",1,planckType<double>()));
	cols.push_back (fitscolumn("SNR_100GHZ","UNITLESS",1,planckType<double>()));
	cols.push_back (fitscolumn("FLUX5R500_143GHZ","DTCMB/TCMB ARCMIN^2",1,planckType<double>()));
	cols.push_back (fitscolumn("SNR_143GHZ","UNITLESS",1,planckType<double>()));
	cols.push_back (fitscolumn("FLUX5R500_217GHZ","DTCMB/TCMB ARCMIN^2",1,planckType<double>()));
	cols.push_back (fitscolumn("SNR_217GHZ","UNITLESS",1,planckType<double>()));
	cols.push_back (fitscolumn("FLUX5R500_353GHZ","DTCMB/TCMB ARCMIN^2",1,planckType<double>()));
	cols.push_back (fitscolumn("SNR_353GHZ","UNITLESS",1,planckType<double>()));
	cols.push_back (fitscolumn("FLUX5R500_545GHZ","DTCMB/TCMB ARCMIN^2",1,planckType<double>()));
	cols.push_back (fitscolumn("SNR_545GHZ","UNITLESS",1,planckType<double>()));
	cols.push_back (fitscolumn("FLUX5R500_857GHZ","DTCMB/TCMB ARCMIN^2",1,planckType<double>()));
	cols.push_back (fitscolumn("SNR_857GHZ","UNITLESS",1,planckType<double>()));
	cols.push_back (fitscolumn("CHI2","UNITLESS",1,planckType<double>()));
	cols.push_back (fitscolumn("SNRF","UNITLESS",1,planckType<double>()));
#else
	cols.push_back (fitscolumn("GLON","DEG",1,planckType<float>()));
	cols.push_back (fitscolumn("GLAT","DEG",1,planckType<float>()));
	cols.push_back (fitscolumn("RAJ2000","DEG",1,planckType<float>()));
	cols.push_back (fitscolumn("DEJ2000","DEG",1,planckType<float>()));
	cols.push_back (fitscolumn("ERR_RAD","ARCMIN",1,planckType<float>()));
	cols.push_back (fitscolumn("CY5R500","ARCMIN^2",1,planckType<double>()));
	cols.push_back (fitscolumn("TS","ARCMIN",1,planckType<float>()));
	cols.push_back (fitscolumn("SNR","UNITLESS",1,planckType<float>()));
	cols.push_back (fitscolumn("CY500","ARCMIN^2",1,planckType<double>()));
	cols.push_back (fitscolumn("T500","ARCMIN",1,planckType<float>()));
	cols.push_back (fitscolumn("CY5R500F","ARCMIN^2",1,planckType<double>()));
	cols.push_back (fitscolumn("TSF","ARCMIN",1,planckType<float>()));
	cols.push_back (fitscolumn("CY500F","ARCMIN^2",1,planckType<double>()));
	cols.push_back (fitscolumn("T500F","ARCMIN",1,planckType<float>()));
	cols.push_back (fitscolumn("ERR_CY5R500F","ARCMIN^2",1,planckType<double>()));
	cols.push_back (fitscolumn("ERR_TSF","ARCMIN",1,planckType<float>()));
	cols.push_back (fitscolumn("ERR_CY500F","ARCMIN^2",1,planckType<double>()));
	cols.push_back (fitscolumn("ERR_T500F","ARCMIN",1,planckType<float>()));
	cols.push_back (fitscolumn("FLUX5R500_100GHZ","DTCMB/TCMB ARCMIN^2",1,planckType<float>()));
	cols.push_back (fitscolumn("SNR_100GHZ","UNITLESS",1,planckType<float>()));
	cols.push_back (fitscolumn("FLUX5R500_143GHZ","DTCMB/TCMB ARCMIN^2",1,planckType<float>()));
	cols.push_back (fitscolumn("SNR_143GHZ","UNITLESS",1,planckType<float>()));
	cols.push_back (fitscolumn("FLUX5R500_217GHZ","DTCMB/TCMB ARCMIN^2",1,planckType<float>()));
	cols.push_back (fitscolumn("SNR_217GHZ","UNITLESS",1,planckType<float>()));
	cols.push_back (fitscolumn("FLUX5R500_353GHZ","DTCMB/TCMB ARCMIN^2",1,planckType<float>()));
	cols.push_back (fitscolumn("SNR_353GHZ","UNITLESS",1,planckType<float>()));
	cols.push_back (fitscolumn("FLUX5R500_545GHZ","DTCMB/TCMB ARCMIN^2",1,planckType<float>()));
	cols.push_back (fitscolumn("SNR_545GHZ","UNITLESS",1,planckType<float>()));
	cols.push_back (fitscolumn("FLUX5R500_857GHZ","DTCMB/TCMB ARCMIN^2",1,planckType<float>()));
	cols.push_back (fitscolumn("SNR_857GHZ","UNITLESS",1,planckType<float>()));
	cols.push_back (fitscolumn("CHI2","UNITLESS",1,planckType<float>()));
	cols.push_back (fitscolumn("SNRF","UNITLESS",1,planckType<float>()));
#endif
	if(!((data.Storage_[0]).Cat_.ScaleLikeNoise_.empty()))
	{
		const int nCols((data.Storage_[0]).Cat_.ScaleLikeNoise_.size());

		for(int i=0;i<nCols;++i)
		{
			sprintf(buffer,"%02d",i);
			StrBuffer = std::string(buffer);
			cols.push_back(fitscolumn(std::string("R500_") + StrBuffer,"ARCMIN",1,planckType<float>()));
			cols.push_back(fitscolumn(std::string("Y500_") + StrBuffer,"ARCMIN^2",1,planckType<float>()));
			cols.push_back(fitscolumn(std::string("ERRY500_") + StrBuffer,"ARCMIN^2",1,planckType<float>()));		
		}

		NColumns_ = cols.size();
	}

	FitsHandle_.insert_bintab(cols);
}
//
void		CatWriterOutputBinFile_PSFITS::AddHeaderPCCFields(void)
{
	std::vector<fitscolumn> cols;

	cols.push_back (fitscolumn("SNR","UNITLESS",1,planckType<float>()));
	cols.push_back (fitscolumn("GAUSSIANITY","UNITLESS",1,planckType<float>()));
	cols.push_back (fitscolumn("LN_RHO","UNITLESS",1,planckType<float>()));
	cols.push_back (fitscolumn("GLON","DEG",1,planckType<float>()));
	cols.push_back (fitscolumn("GLAT","DEG",1,planckType<float>()));
	cols.push_back (fitscolumn("RAJ2000","DEG",1,planckType<float>()));
	cols.push_back (fitscolumn("DEJ2000","DEG",1,planckType<float>()));
	cols.push_back (fitscolumn("ERR_RAD","ARCMIN",1,planckType<float>()));
	cols.push_back (fitscolumn("FLUX_ML","mJys",1,planckType<float>()));
	cols.push_back (fitscolumn("RADIUS_ML","ARCMIN",1,planckType<float>()));
	cols.push_back (fitscolumn("FLUX_B","mJys",1,planckType<float>()));
	cols.push_back (fitscolumn("RADIUS_B","ARCMIN",1,planckType<float>()));
	cols.push_back (fitscolumn("ERR_FLUX","mJys",1,planckType<float>()));
	cols.push_back (fitscolumn("ERR_RADIUS","ARCMIN",1,planckType<float>()));

	FitsHandle_.insert_bintab(cols);
}
//
void		CatWriterOutputFitsFile::CreateHeaderHelper(const OutputFormatType& data)
{

	const int	NCol(GetNColumns());
	const int	NRows(static_cast<int>(data.Storage_.size()));
	std::string	tStr(Wstr2Str(data.Header_.ExtHeader_.DataSetName_));

	FitsHandle_.set_key("COORDSTYPE",static_cast<int>(data.Header_.PCCHeader_.CoordsType_),"Coordinates type");
	FitsHandle_.set_key("COORDSSYS",static_cast<int>(data.Header_.PCCHeader_.CoordSystem_),"Coordinates system");
	FitsHandle_.set_key("EPOCH",static_cast<double>(data.Header_.PCCHeader_.Epoch_),"Coordinates epoch");
	FitsHandle_.set_key("NCOLUMNS",NCol,"Number of columns");
	FitsHandle_.set_key("NROWS",NRows,"Number of rows");
	FitsHandle_.set_key("DETECTTYPE",data.Header_.PCCHeader_.DetectionType_,"Detection type");
	FitsHandle_.set_key("ESTIMATOR",data.Header_.PCCHeader_.Estimator_,"Estimator");
	FitsHandle_.set_key("PRIORTYPE",data.Header_.PCCHeader_.PriorsType_,"Prior type");
	FitsHandle_.set_key("SZ_PROFILE",(data.Header_.PCCHeader_.SZ_params_.SZ_Profile_==2)?std::string("GNFW"):((data.Header_.PCCHeader_.SZ_params_.SZ_Profile_==1)?std::string("BETA"):std::string("UNKNOWN")),"Cluster presure profile");
	FitsHandle_.set_key("VIRIALRATIO",static_cast<double>(data.Header_.PCCHeader_.SZ_params_.VirialRatio_),"Cluster virial ratio");
	FitsHandle_.set_key("GNFW_ALPHA",static_cast<double>(data.Header_.PCCHeader_.SZ_params_.MNFW_alpha_),"Alpha parameter");
	FitsHandle_.set_key("GNFW_BETA",static_cast<double>(data.Header_.PCCHeader_.SZ_params_.MNFW_beta_),"Beta parameter");
	FitsHandle_.set_key("GNFW_GAMMA",static_cast<double>(data.Header_.PCCHeader_.SZ_params_.MNFW_gamma_),"Gamma parameter");
	FitsHandle_.set_key("GNFW_C500",static_cast<double>(data.Header_.PCCHeader_.SZ_params_.MNFW_C500_),"GNFW_C500");
	FitsHandle_.set_key("GNFW_RATIO_CY500CYR500",static_cast<double>(data.Header_.PCCHeader_.SZ_params_.MNFW_Ratio_CY500CYR500_),"GNFW_RATIO_CY500CYR500");
	FitsHandle_.set_key("FLUXCALIB_CTE",static_cast<double>(data.Header_.PCCHeader_.SZ_params_.FluxCalibCte_),"Flux calibration constante");
	FitsHandle_.set_key("RADIUSCALIB_CTE",static_cast<double>(data.Header_.PCCHeader_.SZ_params_.RadiusCalCte_),"Radius calibration constante");
	FitsHandle_.set_key("DATASET",tStr,"Data set name");
}
//
void		CatWriterOutputBinFile_PCCFITS::WriteData2Columns(const OutputFormatType& data,ColumnsCollType& Columns)
{
	const int	tNCols(GetNColumns());
	const OutputFormatType::StorageType& DataSto(data.Storage_);

	Columns.clear();
	Columns.resize(tNCols);

	const	int		DataSz(static_cast<int>(DataSto.size()));
	int				index(0);

	for(int i=0;i<tNCols;++i)
	{
		Columns[i].alloc(DataSz);
		Columns[i].fill(SZCATALOGUEDEFAULTVALUE);
	}

	pointing		tempEqu;
	Trafo			TransCoordsEqu(data.Header_.PCCHeader_.Epoch_,2000.0,data.Header_.PCCHeader_.CoordSystem_,Equatorial);

	int Index(0);
	for(int i=0;i<DataSz;++i,Index=0)
	{
		if((DataSto[i].Cat_.NormalAmpl_<0.0) || (DataSto[i].Cat_.GalLongDegs_ < -1.0e30))
			continue;
		
		(Columns[Index++])[i] = DataSto[i].Cat_.GalLongDegs_;
		(Columns[Index++])[i] = DataSto[i].Cat_.GalLatDegs_;
		tempEqu = TransCoordsEqu(pointing(PIOVER2 - (DataSto[i].Cat_.GalLatDegs_ / RAD2DEGREE),DataSto[i].Cat_.GalLongDegs_ / RAD2DEGREE));
		(Columns[Index++])[i] = tempEqu.phi * RAD2DEGREE;
		(Columns[Index++])[i] = (90.0 - tempEqu.theta * RAD2DEGREE);
		(Columns[Index++])[i] = ((DataSto[i].Cat_.ErrorBars_.TotalPosErrorBar_ < 0.0) ? SZCATALOGUEDEFAULTVALUE : DataSto[i].Cat_.ErrorBars_.TotalPosErrorBar_);
		(Columns[Index++])[i] = DataSto[i].Cat_.FluxComptGLRT_;
		(Columns[Index++])[i] = DataSto[i].Cat_.RadiusGLRT_;
//		(Columns[Index++])[i] = ((DataSto[i].Cat_.NormalAmpl_ < DataSto[i].Cat_.DetectSigma_)?DataSto[i].Cat_.NormalAmpl_:DataSto[i].Cat_.DetectSigma_);
		(Columns[Index++])[i] = DataSto[i].Cat_.NormalAmpl_;
		(Columns[Index++])[i] = DataSto[i].Cat_.FluxComptGLRT_ / data.Header_.PCCHeader_.SZ_params_.MNFW_Ratio_CY500CYR500_;
		(Columns[Index++])[i] = DataSto[i].Cat_.RadiusGLRT_ * data.Header_.PCCHeader_.SZ_params_.MNFW_C500_;
		(Columns[Index++])[i] = DataSto[i].Cat_.FluxCompt_;
		(Columns[Index++])[i] = DataSto[i].Cat_.Radius_;
		(Columns[Index++])[i] = DataSto[i].Cat_.FluxCompt_ / data.Header_.PCCHeader_.SZ_params_.MNFW_Ratio_CY500CYR500_;
		(Columns[Index++])[i] = DataSto[i].Cat_.Radius_ * data.Header_.PCCHeader_.SZ_params_.MNFW_C500_;
		(Columns[Index++])[i] = ((DataSto[i].Cat_.ErrorBars_.FluxErrorBar_ < 0.0) ? SZCATALOGUEDEFAULTVALUE : DataSto[i].Cat_.ErrorBars_.FluxErrorBar_);
		(Columns[Index++])[i] = ((DataSto[i].Cat_.ErrorBars_.RadiusErrorBar_ < 0.0) ? SZCATALOGUEDEFAULTVALUE : DataSto[i].Cat_.ErrorBars_.RadiusErrorBar_);
		(Columns[Index++])[i] = ((DataSto[i].Cat_.ErrorBars_.FluxErrorBar_ < 0.0) ? SZCATALOGUEDEFAULTVALUE : DataSto[i].Cat_.ErrorBars_.FluxErrorBar_ / data.Header_.PCCHeader_.SZ_params_.MNFW_Ratio_CY500CYR500_);
		(Columns[Index++])[i] = ((DataSto[i].Cat_.ErrorBars_.RadiusErrorBar_ < 0.0) ? SZCATALOGUEDEFAULTVALUE : DataSto[i].Cat_.ErrorBars_.RadiusErrorBar_ * data.Header_.PCCHeader_.SZ_params_.MNFW_C500_);
		(Columns[Index++])[i] = DataSto[i].Ext_.Flux5R500_100_;
		(Columns[Index++])[i] = DataSto[i].Ext_.SNR_100_;
		(Columns[Index++])[i] = DataSto[i].Ext_.Flux5R500_143_;
		(Columns[Index++])[i] = DataSto[i].Ext_.SNR_143_;
		(Columns[Index++])[i] = DataSto[i].Ext_.Flux5R500_217_;
		(Columns[Index++])[i] = DataSto[i].Ext_.SNR_217_;
		(Columns[Index++])[i] = DataSto[i].Ext_.Flux5R500_353_;
		(Columns[Index++])[i] = DataSto[i].Ext_.SNR_353_;
		(Columns[Index++])[i] = DataSto[i].Ext_.Flux5R500_545_;
		(Columns[Index++])[i] = DataSto[i].Ext_.SNR_545_;
		(Columns[Index++])[i] = DataSto[i].Ext_.Flux5R500_857_;
		(Columns[Index++])[i] = DataSto[i].Ext_.SNR_857_;
		(Columns[Index++])[i] = DataSto[i].Ext_.CHI2_;
		(Columns[Index++])[i] = DataSto[i].Ext_.SNRF_;

		if(DataSto[i].Cat_.ScaleLikeNoise_.empty())
			continue;
		
		const int NScales(DataSto[i].Cat_.ScaleLikeNoise_.size());

		for(int j=0;j<NScales;++j)
		{
			(Columns[Index++])[i] = (DataSto[i].Cat_.ScaleLikeNoise_[j]).Scale_;
			(Columns[Index++])[i] = (DataSto[i].Cat_.ScaleLikeNoise_[j]).Like_;
			(Columns[Index++])[i] = (DataSto[i].Cat_.ScaleLikeNoise_[j]).Noise_;
		}
	}
}
//
void		CatWriterOutputBinFile_QA_Contours::WriteData2Columns(const OutputFormatType& data,ColumnsCollType& Columns)
{
	const int	tNCols(GetNColumns());
	const OutputFormatType::StorageType& DataSto(data.Storage_);

	Columns.clear();
	Columns.resize(tNCols);

	const	int		DataSz(static_cast<int>(DataSto.size()));
	int				index(0);

	for(int i=0;i<tNCols;++i)
	{
		Columns[i].alloc(DataSz);
		Columns[i].fill(SZCATALOGUEDEFAULTVALUE);
	}

	int Index(0);
	for(int i=0;i<DataSz;++i,Index=0)
	{
		for(int j=0;j<tNCols;++j)
		{
			(Columns[Index++])[i] = DataSto[i].Cat_.QAResult_[j];
		}
	}
}
//
void		CatWriterOutputBinFile_PSFITS::WriteData2Columns(const OutputFormatType& data,ColumnsCollType& Columns)
{
	const int	tNCols(GetNColumns());

	const OutputFormatType::StorageType& DataSto(data.Storage_);

	Columns.clear();
	Columns.resize(tNCols);

	const	int		DataSz(static_cast<int>(DataSto.size()));
	int				index(0);

	for(int i=0;i<tNCols;++i)
	{
		Columns[i].alloc(DataSz);
		Columns[i].fill(SZCATALOGUEDEFAULTVALUE);
	}

	pointing		tempEqu;
	Trafo			TransCoordsEqu(data.Header_.PCCHeader_.Epoch_,2000.0,data.Header_.PCCHeader_.CoordSystem_,Equatorial);

	int Index(0);
	for(int i=0;i<DataSz;++i,Index=0)
	{
		if((DataSto[i].Cat_.NormalAmpl_<0.0) || (DataSto[i].Cat_.GalLongDegs_ < -1.0e30))
			continue;

		(Columns[Index++])[i] = ((DataSto[i].Cat_.NormalAmpl_ < DataSto[i].Cat_.DetectSigma_)?DataSto[i].Cat_.NormalAmpl_:DataSto[i].Cat_.DetectSigma_);
//		(Columns[Index++])[i] = DataSto[i].Cat_.NormalAmpl_;
		(Columns[Index++])[i] = DataSto[i].Cat_.Gaussianity_;
		(Columns[Index++])[i] = DataSto[i].Cat_.lnRho_;
		(Columns[Index++])[i] = DataSto[i].Cat_.GalLongDegs_;
		(Columns[Index++])[i] = DataSto[i].Cat_.GalLatDegs_;
		tempEqu = TransCoordsEqu(pointing(PIOVER2 - (DataSto[i].Cat_.GalLatDegs_ / RAD2DEGREE),DataSto[i].Cat_.GalLongDegs_ / RAD2DEGREE));
		(Columns[Index++])[i] = tempEqu.phi * RAD2DEGREE;
		(Columns[Index++])[i] = (90.0 - tempEqu.theta * RAD2DEGREE);
		(Columns[Index++])[i] = ((DataSto[i].Cat_.ErrorBars_.TotalPosErrorBar_ < 0.0) ? SZCATALOGUEDEFAULTVALUE : DataSto[i].Cat_.ErrorBars_.TotalPosErrorBar_);
		(Columns[Index++])[i] = DataSto[i].Cat_.FluxComptGLRT_;
		(Columns[Index++])[i] = DataSto[i].Cat_.RadiusGLRT_;
		(Columns[Index++])[i] = DataSto[i].Cat_.FluxCompt_;
		(Columns[Index++])[i] = DataSto[i].Cat_.Radius_;
		(Columns[Index++])[i] = ((DataSto[i].Cat_.ErrorBars_.FluxErrorBar_ < 0.0) ? SZCATALOGUEDEFAULTVALUE : DataSto[i].Cat_.ErrorBars_.FluxErrorBar_);
		(Columns[Index++])[i] = ((DataSto[i].Cat_.ErrorBars_.RadiusErrorBar_ < 0.0) ? SZCATALOGUEDEFAULTVALUE : DataSto[i].Cat_.ErrorBars_.RadiusErrorBar_);
	}
}
//
}

