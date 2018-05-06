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

#include "ZEUS_InOutTxtFile.h"

//------------------------------------

namespace Zeus
{
//
void		TxtCollBaseWriter::Initialize(void)
{
	stream_.open(Zeus::Wstr2Str(fname_).c_str(),mode_);
	if(stream_.fail())
		throw Zeus::libException(ERROR_COD_MYINOUT_FILEWRITE,ERROR_MSG_MYINOUT_FILEWRITE,fname_);
	stream_.clear();
	stream_.exceptions(std::ios::failbit | std::ios::badbit);
}

//
void		TxtFileReaderRaw::Open(void)
{
	if(istrm_)
	{
		if(!(istrm_->fail())) goto Open_end;
		errTxtFileReader(ERROR_COD_ZEUSERROPENFL,ERROR_MSG_ZEUSERROPENFL,fname_);
	}
	istrm_ = new std::wifstream(Wstr2Str(fname_).c_str());
	if(!istrm_ || istrm_->fail())
	{
		if(istrm_)
		{delete istrm_;istrm_ = 0;}
		errTxtFileReader(ERROR_COD_ZEUSERROPENFL,ERROR_MSG_ZEUSERROPENFL,fname_);
	}
Open_end:
	Reset();
}

//
TxtFileReaderRaw::TxtFileReaderRaw(const std::wstring& fname,bool NoOpen,const std::locale* loc)
	:istrm_(0),fname_(fname),readLns_(0),parsedLns_(0),loc_(loc),LnOnBuffer_(0)
{
	if(!NoOpen) Open();
}

//
TxtFileReaderRaw::TxtFileReaderRaw(std::wistream* istream,const std::locale* loc)
	:istrm_(0),readLns_(0),parsedLns_(0),loc_(loc),LnOnBuffer_(0)
{
	fname_.clear();
	if(!istream || !(*istream)) errTxtFileReader(ERROR_COD_ZEUSGENERRFIL,ERROR_MSG_ZEUSGENERRFIL);
	istrm_ = istream;
	istrm_->clear();
	oldIoState_ = istrm_->exceptions();
	istrm_->exceptions(std::ios::failbit | std::ios::badbit);
	istrm_->seekg(0);
}
//
bool		TxtFileReaderRaw::ReadNextLn(void)
{
	if(!istrm_) errTxtFileReader(ERROR_COD_ZEUSERROPENFL,ERROR_MSG_ZEUSERROPENFL,fname_);
	wchar_t		StrBuffer[FILESTRBUFFER];

read_again:
	if(!LnOnBuffer_)
	{
		CurrStr_.clear();
		try{
			if(!istrm_->getline(StrBuffer,FILESTRBUFFER))
			{return false;}
			CurrStr_ = std::wstring(StrBuffer);
			++readLns_;
		}
		catch(...){
			if(!istrm_->eof()){errTxtFileReader(ERROR_COD_ZEUSGENERRFIL,ERROR_MSG_ZEUSGENERRFIL);}
			else{return false;}
		}
		if(CurrStr_.empty()) goto read_again;
	}
	else
	{LnOnBuffer_ = 0;}
	return true;
}
//
int			ParameterFileReader::SplitStr(const std::wstring& in,std::wstring& fld,std::wstring& value)
{

	std::wstring temp(FullTrim(in,std::ctype_base::space));
	std::wstring::size_type slen;

	if(temp.empty()) return 0;
	std::wstring::size_type pos(temp.find(COMMENTCHAR));
	if(!pos) return 0;
	if(pos == std::wstring::npos) slen=temp.size();
	else slen=pos;

	std::wstring temp1(temp.substr(0,slen));
	pos = temp1.find(ASSIGNCHAR);
	if(!pos || (pos==std::wstring::npos))
		return -1;
	fld		= ToCase(RightTrim(temp1.substr(0,pos),std::ctype_base::space),true);
	value	= FullTrim(temp1.substr(++pos),std::ctype_base::space);
	return 1;
}

//
int			ParameterFileReader::GetValue(const std::wstring& in,DBField& fld)
{
	long				tInt;
	double				tDouble;

	if(in.empty())
	{fld = std::wstring(L"");return 0;}
	if(!get_number(in,tInt,CONVDEC))
	{fld = static_cast<int>(tInt);return 1;}
	if(!get_number(in,tDouble,CONVDEC))
	{fld = tDouble;return 2;}
	if(in.at(0) == L'@')
	{
		if(in.size() <= 1) return -1;
		fld = ps_VariableType(in.substr(1));
		return 3;
	}
	fld = in;return 0;

}

//
int			ParameterFileReader::do_parseLn(const std::wstring& str)
{
	std::wstring	fld;
	std::wstring	valueStr;
	DBField			value;
	int				stat;

	if((stat = SplitStr(str,fld,valueStr)) < 0)
		errTxtFileReader(ERROR_COD_ZEUSGLBPARAMFL,ERROR_MSG_ZEUSGLBPARAMFL,GetFileName());
	if(!stat) return 1;
	if((stat = GetValue(valueStr,value)) < 0)
		errTxtFileReader(ERROR_COD_ZEUSGLBPARAMFL,ERROR_MSG_ZEUSGLBPARAMFL,GetFileName());
	Insert1element(fld,value);
	return 1;
}
//
int			ReadPatchesGeomInfoTxtFile::do_parseLn(const std::wstring& str)
{
	int							NItems;
	double						ptg_T,ptg_P,X0Y0_T,X0Y0_P,XLY0_T,XLY0_P,X0YL_T,X0YL_P,XLYL_T,XLYL_P,SrcPTG_t,SrcPTG_p;
	T::StorageType::value_type	temp;

	if((NItems = swscanf(str.c_str(),GEOM_LINES_SCANF,
		&(temp.PatchNumber_),
		&(temp.BPixel_),
		&(temp.SrcXCoord_),
		&(temp.SrcYCoord_),
		&(temp.Spin_),
		&ptg_T,
		&ptg_P,
		&X0Y0_T,
		&X0Y0_P,
		&XLY0_T,
		&XLY0_P,
		&X0YL_T,
		&X0YL_P,
		&XLYL_T,
		&XLYL_P,
		&(temp.InitRejectRatio_),
		&(temp.FinalRejectRatio_),
		&(temp.PatchValid_),
		&(temp.PredFlux_),
		&(temp.PredRadius_),
		&(temp.SNR_),
		&(temp.ErrRadius_),
		&(temp.ErrFlux_),
		&(temp.ErrPos_),
		&(temp.SrcIndex_),
		&(temp.ErrRadiusHigh_),
		&(temp.ErrRadiusLow_),
		&(temp.ErrFluxHigh_),
		&(temp.ErrFluxLow_),
		&(temp.FluxBay_),
		&(temp.RadiusBay_),
		&SrcPTG_t,
		&SrcPTG_p,
		&(temp.QAIN_SrcPtg_.theta),
		&(temp.QAIN_SrcPtg_.phi),
		&(temp.QAIN_CyR500),
		&(temp.QAIN_T500),
		&(temp.QAIN_detectable)
		)) != GEOM_LINES_NFIELDS)
		return -1;

	temp.X0Y0Ptg_.theta		= ptg_T;
	temp.X0Y0Ptg_.phi		= ptg_P;
	temp.X0Y0_.theta		= X0Y0_T;
	temp.XLY0_.theta		= XLY0_T;
	temp.X0YL_.theta		= X0YL_T;
	temp.XLYL_.theta		= XLYL_T;
	temp.X0Y0_.phi			= X0Y0_P;
	temp.XLY0_.phi			= XLY0_P;
	temp.X0YL_.phi			= X0YL_P;
	temp.XLYL_.phi			= XLYL_P;
	temp.SourcePtg_.theta	= SrcPTG_t;
	temp.SourcePtg_.phi		= SrcPTG_p;

	data_.Storage_.push_back(temp);
	return NItems;
}
//
void		ReadPatchesGeomInfoTxtFile::ReadHeader(void)
{
	if((!ReadNextLn()) || (swscanf(GetCurrStrRef().c_str(),GEOM_HEADER_SCANF,
		&(data_.Header_.NSide_),
		&(data_.Header_.PtchSz_),
		&(data_.Header_.PtchBorder_),
		&(data_.Header_.NPatchesPerMainPix_),
		&(data_.Header_.NTotalPatches_),
		&(data_.Header_.CollListSz_),
		&(data_.Header_.GalacticCut_)
		)!=GEOM_HEADER_NFIELDS))
		 goto ReadHeader_Error;

	return;
ReadHeader_Error:
errTxtFileReader(ERROR_COD_MYINOUT_NOHEADER,ERROR_MSG_MYINOUT_NOHEADER,GetFileName());
}
//
void		WritePatchesGeomInfoTxtFile::WriteAllLns(const T::StorageType  & coll)
{
	T::StorageType::const_iterator	piv(coll.begin());
	T::StorageType::const_iterator	const end(coll.end());

	for(;piv != end;++piv)
	{Write1Line(*piv);}
}
//
void		WritePatchesGeomInfoTxtFile::WriteHeader(const T::HeaderType& head)
{
	wchar_t	buffer[BUFFERMAXCHAR];
	
	PRINTINTOBUFFERFUNCT (buffer,BUFFERMAXCHAR,GEOM_HEADER_PRINTF,head.NSide_,head.PtchSz_,head.PtchBorder_ ,
		head.NPatchesPerMainPix_,head.NTotalPatches_,head.CollListSz_,head.GalacticCut_);
	
	stream_.WriteLn(std::wstring(buffer));
}

//
void		WritePatchesGeomInfoTxtFile::Write1Line(const T::StorageType::value_type & pt)
{
	wchar_t	buffer[BUFFERMAXCHAR];

	PRINTINTOBUFFERFUNCT (buffer,BUFFERMAXCHAR,GEOM_LINES_PRINTF,
		pt.PatchNumber_,
		pt.BPixel_,
		pt.SrcXCoord_,
		pt.SrcYCoord_,
		pt.Spin_,
		pt.X0Y0Ptg_.theta,
		pt.X0Y0Ptg_.phi,
		pt.X0Y0_.theta,
		pt.X0Y0_.phi,
		pt.XLY0_.theta,
		pt.XLY0_.phi,
		pt.X0YL_.theta,
		pt.X0YL_.phi,
		pt.XLYL_.theta,
		pt.XLYL_.phi,
		pt.InitRejectRatio_,
		pt.FinalRejectRatio_,
		pt.PatchValid_,
		pt.PredFlux_,
		pt.PredRadius_,
		pt.SNR_,
		pt.ErrRadius_,
		pt.ErrFlux_,
		pt.ErrPos_,
		pt.SrcIndex_,
		pt.ErrRadiusHigh_,
		pt.ErrRadiusLow_,
		pt.ErrFluxHigh_,
		pt.ErrFluxLow_,
		pt.FluxBay_,
		pt.RadiusBay_,
		pt.SourcePtg_.theta,
		pt.SourcePtg_.phi,
		pt.QAIN_SrcPtg_.theta,
		pt.QAIN_SrcPtg_.phi,
		pt.QAIN_CyR500,
		pt.QAIN_T500,
		pt.QAIN_detectable
		);
	stream_.WriteLn(std::wstring(buffer));
}
//
int			ReadNonBlindPtgsTxtFile::ReadHeader(void)
{
	int		temp;
	int		NFields;
	if(!(ReadNextLn())) return false;

	CatalogueFormatType::HeaderType	hTemp;

	if((NFields = swscanf(GetCurrStrRef().c_str(),
		CatalogueFormatType::HeaderType::CatHeaderReadStr_,
		&(hTemp.CoordsType_),
		&(temp),
		&(hTemp.Epoch_),
		&(hTemp.NColumns_),
		&(hTemp.NRows_),
		&(hTemp.DetectionType_),
		&(hTemp.Estimator_),
		&(hTemp.PriorsType_),
		&(hTemp.SZ_params_.SZ_Profile_),
		&(hTemp.SZ_params_.VirialRatio_),
		&(hTemp.SZ_params_.MNFW_alpha_),
		&(hTemp.SZ_params_.MNFW_beta_),
		&(hTemp.SZ_params_.MNFW_gamma_),
		&(hTemp.SZ_params_.MNFW_C500_),
		&(hTemp.SZ_params_.MNFW_Ratio_CY500CYR500_),
		&(hTemp.SZ_params_.FluxCalibCte_),
		&(hTemp.SZ_params_.RadiusCalCte_),
		&(hTemp.CollLstSz_)
		)) != CatalogueFormatType::HeaderType::CatHeaderNColumns_)
	{
		data_.Header_				= CatalogueFormatType::HeaderType();
		UnreadLn();
	}
	else
	{
		data_.Header_				= hTemp;
		data_.Header_.CoordSystem_	= static_cast<coordsys>(temp);
	}

	return true;
}
//
int			ReadNonBlindPtgsTxtFile::do_parseLn(const std::wstring& str)
{
	int	NItems;
	wchar_t					buffer[BUFFERMAXCHAR];
	wchar_t					buffer1[BUFFERMAXCHAR];
	std::vector<double>		ScalesColl;

	buffer[0] = L'\0';
	buffer1[0] = L'\0';

	CatalogueFormatType::StorageType::value_type	tLn;
	if((NItems = swscanf(str.c_str(),
		CatalogueFormatType::HeaderType::CatReadStr_,
		&(tLn.ID_),
		&(tLn.Patch_),
		&(tLn.lnRho_),
		&(tLn.NormalAmpl_),
		&(tLn.DetectSigma_),
		&(tLn.GalLatDegs_),
		&(tLn.GalLongDegs_),
		&(tLn.PatchGalLatDegs_),
		&(tLn.PatchGalLongDegs_),
		&(tLn.PatchSpin_),
		&(tLn.FluxCompt_),
		&(tLn.FluxComptGLRT_),
		&(tLn.Radius_),
		&(tLn.RadiusGLRT_),
		&(tLn.ErrorBars_.FluxErrorBar_),
		&(tLn.ErrorBars_.RadiusErrorBar_),
		&(tLn.ErrorBars_.TotalPosErrorBar_),
		&(tLn.Gaussianity_),
		&(tLn.lnEvidence_),
		&(tLn.lnPenaltySrc_),
		&(tLn.PatchMF_sigmaSqr_),
		&(tLn.SZ_ConversionCte_),
		&(tLn.ErrorBars_.LowFluxErrorBar_),
		&(tLn.ErrorBars_.HighFluxErrorBar_),
		&(tLn.ErrorBars_.LowRadiusErrorBar_),
		&(tLn.ErrorBars_.HighRadiusErrorBar_),
		&(tLn.ErrorBars_.DegenCorr_),
		&(tLn.ErrorBars_.DegenOrd_),
		&(tLn.ErrorBars_.DegenOrdYErr_),
		&(tLn.ErrorBars_.DegenSlope_),
		&(tLn.ErrorBars_.DegenSlopeYErr_),
		&(tLn.CollLstIndex_),
		buffer,
		buffer1
		)) < 32)
		return -1;
	
	tLn.ScaleLikeNoise_.clear();

	if(
		Zeus::GetNumbersHomogenColl(std::wstring(buffer),std::wstring(L":"),ScalesColl) &&
		(!(ScalesColl.size() % 3))
	)
	{
		ScaleLikeNoiseColl::value_type			t;
		const int								N_Ele(ScalesColl.size());
		for(int i=0;i<N_Ele;i+=3)
		{
			t.Scale_	= ScalesColl[i];
			t.Like_		= ScalesColl[i+1];
			t.Noise_	= ScalesColl[i+2];

			tLn.ScaleLikeNoise_.push_back(t);
		}
	}

	ScalesColl.clear();
	if(
		Zeus::GetNumbersHomogenColl(std::wstring(buffer1),std::wstring(L":"),ScalesColl) &&
		(ScalesColl.size() == QARESULTSSZ)
	)
	{
		for(int i=0;i<QARESULTSSZ;++i)
		{
			tLn.QAResult_[i] = ScalesColl[i];
		}
	}

	data_.Storage_.push_back(tLn);

	return NItems;
}
//
int			ReadQA_ColLstText::do_parseLn(const std::wstring& str)
{
	int	NItems;

	QA_CltLstCatalogueType::StorageType::value_type	tLn;
	if((NItems = swscanf(str.c_str(),
		L"%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le",
		&(tLn.IN_theta),
		&(tLn.IN_phi),
		&(tLn.IN_CyR500),
		&(tLn.IN_T500),
		&(tLn.IN_M500),
		&(tLn.IN_z),
		&(tLn.IN_Vrec),
		&(tLn.IN_detectable),
		&(tLn.Cat_GLAT),
		&(tLn.Cat_GLON),
		&(tLn.Cat_yc),
		&(tLn.Cat_yc_error),
		&(tLn.Cat_CY5R500),
		&(tLn.Cat_LOWER_ERR_CY5R500F),
		&(tLn.Cat_UPPER_ERR_CY5R500F),
		&(tLn.Cat_T500),
		&(tLn.Cat_LOWER_ERR_T500F),
		&(tLn.Cat_UPPER_ERR_T500F),
		&(tLn.Cat_SNR),
		&(tLn.Cat_Patch_No),
		&(tLn.Cat_DEGEN_PARAM),
		&(tLn.Cat_SNR_30GHZ),
		&(tLn.Cat_SNR_44GHZ),
		&(tLn.Cat_SNR_70GHZ),
		&(tLn.Cat_SNR_100GHZ),
		&(tLn.Cat_SNR_143GHZ),
		&(tLn.Cat_SNR_217GHZ),
		&(tLn.Cat_SNR_353GHZ),
		&(tLn.Cat_SNR_545GHZ),
		&(tLn.Cat_SNR_857GHZ),
		&(tLn.Cat_FLUX5R500_30GHZ),
		&(tLn.Cat_FLUX5R500_44GHZ),
		&(tLn.Cat_FLUX5R500_70GHZ),
		&(tLn.Cat_FLUX5R500_100GHZ),
		&(tLn.Cat_FLUX5R500_143GHZ),
		&(tLn.Cat_FLUX5R500_217GHZ),
		&(tLn.Cat_FLUX5R500_353GHZ),
		&(tLn.Cat_FLUX5R500_545GHZ),
		&(tLn.Cat_FLUX5R500_857GHZ),
		&(tLn.Cat_CHI2),
		&(tLn.Cat_ERR_RAD),
		&(tLn.Cat_CY5R500F),
		&(tLn.Cat_T500F)
		)) != 43)
		return -1;
	
	data_.Storage_.push_back(tLn);

	return NItems;
}
// 
int	ReadQA_ProfilesLstText::do_parseLn(const std::wstring& str)
{
	int	NItems;

	QA_ProfilesLstType::StorageType::value_type	tLn;
	if((NItems = swscanf(str.c_str(),
		L"%le,%le,%le,%le,%le",
		&(tLn.Alpha_),
		&(tLn.Beta_),
		&(tLn.Gamma_),
		&(tLn.C500_),
		&(tLn.ConvYtot2Y500_)
		)) != 5)
		return -1;
	
	data_.Storage_.push_back(tLn);

	return NItems;
}
//
void		CatWriterTxtFile::CreatHeader(const CatalogueFormatType & data)
{
	wchar_t	buffer[BUFFERMAXCHAR];
	
	PRINTINTOBUFFERFUNCT (buffer,BUFFERMAXCHAR,
		CatalogueFormatType::HeaderType::CatHeaderPrintStr_,
		data.Header_.CoordsType_,
		data.Header_.CoordSystem_,
		data.Header_.Epoch_,
		data.Header_.NColumns_,
		static_cast<int>(data.Storage_.size()),
		data.Header_.DetectionType_,
		data.Header_.Estimator_,
		data.Header_.PriorsType_,
		data.Header_.SZ_params_.SZ_Profile_,
		data.Header_.SZ_params_.VirialRatio_,
		data.Header_.SZ_params_.MNFW_alpha_,
		data.Header_.SZ_params_.MNFW_beta_,
		data.Header_.SZ_params_.MNFW_gamma_,
		data.Header_.SZ_params_.MNFW_C500_,
		data.Header_.SZ_params_.MNFW_Ratio_CY500CYR500_,
		data.Header_.SZ_params_.FluxCalibCte_,
		data.Header_.SZ_params_.RadiusCalCte_,
		data.Header_.CollLstSz_
		);

	stream_.WriteLn(std::wstring(buffer));

}
//
std::wstring	CatWriterOutputTxtFile::CreatHeaderHelper(const OutputFormatType& data)
{
	wchar_t	buffer[BUFFERMAXCHAR];

	PRINTINTOBUFFERFUNCT (buffer,BUFFERMAXCHAR,
		L"%4d,%4d,%6.1g,%4d,%4d,%4d,%4d,%4d,%4d,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%ls",
		data.Header_.PCCHeader_.CoordsType_,
		data.Header_.PCCHeader_.CoordSystem_,
		data.Header_.PCCHeader_.Epoch_,
		data.Header_.PCCHeader_.NColumns_,
		static_cast<int>(data.Storage_.size()),
		data.Header_.PCCHeader_.DetectionType_,
		data.Header_.PCCHeader_.Estimator_,
		data.Header_.PCCHeader_.PriorsType_,
		data.Header_.PCCHeader_.SZ_params_.SZ_Profile_,
		data.Header_.PCCHeader_.SZ_params_.VirialRatio_,
		data.Header_.PCCHeader_.SZ_params_.MNFW_alpha_,
		data.Header_.PCCHeader_.SZ_params_.MNFW_beta_,
		data.Header_.PCCHeader_.SZ_params_.MNFW_gamma_,
		data.Header_.PCCHeader_.SZ_params_.MNFW_C500_,
		data.Header_.PCCHeader_.SZ_params_.MNFW_Ratio_CY500CYR500_,
		data.Header_.PCCHeader_.SZ_params_.FluxCalibCte_,
		data.Header_.PCCHeader_.SZ_params_.RadiusCalCte_,
		data.Header_.ExtHeader_.DataSetName_.c_str()
		);

	return std::wstring(buffer);
}

//
int			CatWriterTxtFile::WriteCatLinesColl(const CatalogueFormatType & data)
{
	wchar_t												buffer[BUFFERMAXCHAR];
	std::wstring										tString;
	std::wstring										tString1;
	CatalogueFormatType::StorageType::const_iterator	piv(data.Storage_.begin());
	CatalogueFormatType::StorageType::const_iterator	const end(data.Storage_.end());

	for(;piv != end;++piv)
	{
		if (piv->NormalAmpl_ >= 0.0)
		{
//			if (!(piv->ScaleLikeNoise_.empty()))
			if (0)
			{
				ScaleLikeNoiseColl::const_iterator	pivS(piv->ScaleLikeNoise_.begin());
				ScaleLikeNoiseColl::const_iterator	const endS(piv->ScaleLikeNoise_.end());
				tString.clear();
				for(;pivS != endS;++pivS)
				{
					PRINTINTOBUFFERFUNCT(buffer,BUFFERMAXCHAR,L":%4.6g:%4.6g:%4.6g",pivS->Scale_,pivS->Like_,pivS->Noise_);
					tString += std::wstring(buffer);
				}
				tString = tString.substr(1);
			}
			else
			{
				tString = std::wstring(L"####");
			}

			if(piv->QAResult_[0] >= 0.0)
			{
				tString1.clear();
				for(int i=0;i<QARESULTSSZ;++i)
				{
					PRINTINTOBUFFERFUNCT(buffer,BUFFERMAXCHAR,L":%4.6g",piv->QAResult_[i]);
					tString1 += std::wstring(buffer);
				}
				tString1 = tString1.substr(1);
			}
			else
			{tString1 = std::wstring(L"####");}

			PRINTINTOBUFFERFUNCT
				(buffer,BUFFERMAXCHAR,
				CatalogueFormatType::HeaderType::CatPrintStr_,
				piv->ID_,
				piv->Patch_,
				piv->lnRho_,
				piv->NormalAmpl_,
				piv->DetectSigma_,
				piv->GalLatDegs_,
				piv->GalLongDegs_,
				piv->PatchGalLatDegs_,
				piv->PatchGalLongDegs_,
				piv->PatchSpin_,				
				piv->FluxCompt_,
				piv->FluxComptGLRT_,
				piv->Radius_,
				piv->RadiusGLRT_,
				piv->ErrorBars_.FluxErrorBar_,
				piv->ErrorBars_.RadiusErrorBar_,
				piv->ErrorBars_.TotalPosErrorBar_,
				piv->Gaussianity_,
				piv->lnEvidence_,
				piv->lnPenaltySrc_,
				piv->PatchMF_sigmaSqr_,
				piv->SZ_ConversionCte_,
				piv->ErrorBars_.LowFluxErrorBar_,
				piv->ErrorBars_.HighFluxErrorBar_,
				piv->ErrorBars_.LowRadiusErrorBar_,
				piv->ErrorBars_.HighRadiusErrorBar_,
				piv->ErrorBars_.DegenCorr_,
				piv->ErrorBars_.DegenOrd_,
				piv->ErrorBars_.DegenOrdYErr_,
				piv->ErrorBars_.DegenSlope_,
				piv->ErrorBars_.DegenSlopeYErr_,
				piv->CollLstIndex_,
				tString.c_str(),
				tString1.c_str()
				);
		}
		else
		{
			PRINTINTOBUFFERFUNCT
				(buffer,BUFFERMAXCHAR,
				CatalogueFormatType::HeaderType::CatPrintStr_,
				piv->ID_,
				piv->Patch_,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				-1,
				L"####",
				L"####"
				);				
		}
		WriteLn(buffer);
	}
	return true;
}
//
int			CatWriterOutputTxtFile_PS::WriteCatLinesColl(const OutputFormatType& data)
{
	wchar_t	buffer[BUFFERMAXCHAR];
	double	J2000_LongDegs_;
	double	J2000_LatDegs_;


	OutputFormatType::StorageType::const_iterator	piv(data.Storage_.begin());
	OutputFormatType::StorageType::const_iterator	const end(data.Storage_.end());

	Trafo			TransCoordsEqu(data.Header_.PCCHeader_.Epoch_,2000.0,data.Header_.PCCHeader_.CoordSystem_,Equatorial);
	pointing		tempEqu;

	for(;piv != end;++piv)
	{
		if((piv->Cat_.NormalAmpl_ > 0.0) && (piv->Cat_.GalLongDegs_ > -1.0e10))
		{
			tempEqu = TransCoordsEqu(pointing(PIOVER2 - (piv->Cat_.GalLatDegs_ / RAD2DEGREE),piv->Cat_.GalLongDegs_ / RAD2DEGREE));

			J2000_LongDegs_ = tempEqu.phi * RAD2DEGREE;
			J2000_LatDegs_	= 90.0 - (tempEqu.theta * RAD2DEGREE);

			PRINTINTOBUFFERFUNCT
				(buffer,BUFFERMAXCHAR,
				L"%05d,%05d,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%d,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g",
				piv->Cat_.ID_,
				piv->Cat_.Patch_,
				piv->Cat_.NormalAmpl_,
				piv->Cat_.DetectSigma_,
				piv->Cat_.Gaussianity_,
				piv->Ext_.CHI2_,
				piv->Cat_.lnRho_,
				piv->Cat_.lnEvidence_,
				piv->Cat_.lnPenaltySrc_,
				piv->Cat_.GalLongDegs_,
				piv->Cat_.GalLatDegs_,
				J2000_LongDegs_,
				J2000_LatDegs_,
				piv->Cat_.ErrorBars_.TotalPosErrorBar_,
				piv->Cat_.FluxComptGLRT_,
				piv->Cat_.FluxCompt_,
				piv->Cat_.ErrorBars_.FluxErrorBar_,
				piv->Cat_.RadiusGLRT_,
				piv->Cat_.Radius_,
				piv->Cat_.ErrorBars_.RadiusErrorBar_,
				piv->Ext_.Flux5R500_100_,
				piv->Ext_.SNR_100_,
				piv->Ext_.Flux5R500_143_,
				piv->Ext_.SNR_143_,
				piv->Ext_.Flux5R500_217_,
				piv->Ext_.SNR_217_,
				piv->Ext_.Flux5R500_353_,
				piv->Ext_.SNR_353_,
				piv->Ext_.Flux5R500_545_,
				piv->Ext_.SNR_545_,
				piv->Ext_.Flux5R500_857_,
				piv->Ext_.SNR_857_,
				piv->Ext_.Flux5R500_030_,
				piv->Ext_.SNR_030_,
				piv->Ext_.Flux5R500_044_,
				piv->Ext_.SNR_044_,
				piv->Ext_.Flux5R500_070_,
				piv->Ext_.SNR_070_,
				piv->Ext_.NData_,
				piv->Cat_.ErrorBars_.LowFluxErrorBar_,
				piv->Cat_.ErrorBars_.HighFluxErrorBar_,
				piv->Cat_.ErrorBars_.LowRadiusErrorBar_,
				piv->Cat_.ErrorBars_.HighRadiusErrorBar_,
				piv->Cat_.ErrorBars_.DegenCorr_,
				piv->Cat_.ErrorBars_.DegenOrd_,
				piv->Cat_.ErrorBars_.DegenOrdYErr_,
				piv->Cat_.ErrorBars_.DegenSlope_,
				piv->Cat_.ErrorBars_.DegenSlopeYErr_
				);
		}
		else 
		{
			PRINTINTOBUFFERFUNCT
				(buffer,BUFFERMAXCHAR,
				L"%05d,%05d,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%d,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g",
				piv->Cat_.ID_,
				piv->Cat_.Patch_,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				-1,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE
				);				
		}
		WriteLn(buffer);
	}
	return true;
}
//
int			CatWriterOutputTxtFile_SZ::WriteCatLinesColl(const OutputFormatType& data)
{
	wchar_t	buffer[BUFFERMAXCHAR];
	double	J2000_LongDegs_;
	double	J2000_LatDegs_;

	OutputFormatType::StorageType::const_iterator	piv(data.Storage_.begin());
	OutputFormatType::StorageType::const_iterator	const end(data.Storage_.end());

	Trafo			TransCoordsEqu(data.Header_.PCCHeader_.Epoch_,2000.0,data.Header_.PCCHeader_.CoordSystem_,Equatorial);
	pointing		tempEqu;

	for(;piv != end;++piv)
	{
		if((piv->Cat_.NormalAmpl_ > 0.0) && (piv->Cat_.GalLongDegs_ > -1.0e10))
		{
			tempEqu = TransCoordsEqu(pointing(PIOVER2 - (piv->Cat_.GalLatDegs_ / RAD2DEGREE),piv->Cat_.GalLongDegs_ / RAD2DEGREE));

			J2000_LongDegs_ = tempEqu.phi * RAD2DEGREE;
			J2000_LatDegs_	= 90.0 - (tempEqu.theta * RAD2DEGREE);

			PRINTINTOBUFFERFUNCT
				(buffer,BUFFERMAXCHAR,
				L"%05d,%05d,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%d,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g",
				piv->Cat_.ID_,
				piv->Cat_.Patch_,
				piv->Cat_.NormalAmpl_,
				piv->Cat_.DetectSigma_,
				piv->Cat_.Gaussianity_,
				piv->Ext_.CHI2_,
				piv->Cat_.lnRho_,
				piv->Cat_.lnEvidence_,
				piv->Cat_.lnPenaltySrc_,
				piv->Cat_.GalLongDegs_,
				piv->Cat_.GalLatDegs_,
				J2000_LongDegs_,
				J2000_LatDegs_,
				piv->Cat_.ErrorBars_.TotalPosErrorBar_,
				piv->Cat_.FluxComptGLRT_,
				piv->Cat_.FluxCompt_,
				piv->Cat_.ErrorBars_.FluxErrorBar_,
				piv->Cat_.RadiusGLRT_,
				piv->Cat_.Radius_,
				piv->Cat_.ErrorBars_.RadiusErrorBar_,
				piv->Ext_.Flux5R500_100_,
				piv->Ext_.SNR_100_,
				piv->Ext_.Flux5R500_143_,
				piv->Ext_.SNR_143_,
				piv->Ext_.Flux5R500_217_,
				piv->Ext_.SNR_217_,
				piv->Ext_.Flux5R500_353_,
				piv->Ext_.SNR_353_,
				piv->Ext_.Flux5R500_545_,
				piv->Ext_.SNR_545_,
				piv->Ext_.Flux5R500_857_,
				piv->Ext_.SNR_857_,
				piv->Ext_.Flux5R500_030_,
				piv->Ext_.SNR_030_,
				piv->Ext_.Flux5R500_044_,
				piv->Ext_.SNR_044_,
				piv->Ext_.Flux5R500_070_,
				piv->Ext_.SNR_070_,
				piv->Ext_.NData_,
				piv->Cat_.ErrorBars_.LowFluxErrorBar_,
				piv->Cat_.ErrorBars_.HighFluxErrorBar_,
				piv->Cat_.ErrorBars_.LowRadiusErrorBar_,
				piv->Cat_.ErrorBars_.HighRadiusErrorBar_,
				piv->Cat_.ErrorBars_.DegenCorr_,
				piv->Cat_.ErrorBars_.DegenOrd_,
				piv->Cat_.ErrorBars_.DegenOrdYErr_,
				piv->Cat_.ErrorBars_.DegenSlope_,
				piv->Cat_.ErrorBars_.DegenSlopeYErr_
				);
		}
		else 
		{
			PRINTINTOBUFFERFUNCT
				(buffer,BUFFERMAXCHAR,
				L"%05d,%05d,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%d,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g",
				piv->Cat_.ID_,
				piv->Cat_.Patch_,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				-1,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE,
				SZCATALOGUEDEFAULTVALUE
				);				
		}
		WriteLn(buffer);
	}
	return true;
}
//
int			CatWriterOutputTxtFile_QA_Contours::WriteCatLinesColl(const OutputFormatType& data)
{
	wchar_t	buffer[BUFFERMAXCHAR];
	std::wstring	tString;
	OutputFormatType::StorageType::const_iterator	piv(data.Storage_.begin());
	OutputFormatType::StorageType::const_iterator	const end(data.Storage_.end());

	for(;piv != end;++piv)
	{
		tString.clear();
		for(int i=0;i<QARESULTSSZ;++i)
		{
			PRINTINTOBUFFERFUNCT
				(buffer,BUFFERMAXCHAR,L",%9.6g",piv->Cat_.QAResult_[i]);
			tString += std::wstring(buffer);
		}

		WriteLn((tString.substr(1)).c_str());
	}

	return true;
}
//
int			WriteObjResultsTxtFile::WritePeaksColl(const PeakCollType& data)
{
	wchar_t			buffer[BUFFERMAXCHAR];
	std::wstring	tString;
	std::wstring	tString1;

	PeakCollType::const_iterator	piv(data.begin());
	PeakCollType::const_iterator	const end(data.end());

	for(;piv != end;++piv,++ObjectID_)
	{
//		if (!(piv->ScaleLikeNoise_.empty()))
		if(0)
		{
			ScaleLikeNoiseColl::const_iterator	pivS(piv->ScaleLikeNoise_.begin());
			ScaleLikeNoiseColl::const_iterator	const endS(piv->ScaleLikeNoise_.end());
			tString.clear();
			for(;pivS != endS;++pivS)
			{
				PRINTINTOBUFFERFUNCT(buffer,BUFFERMAXCHAR,L":%4.6g:%4.6g:%4.6g",pivS->Scale_,pivS->Like_,pivS->Noise_);
				tString += std::wstring(buffer);
			}
			tString = tString.substr(1);
		}
		else
		{tString = std::wstring(L"####");}
		if(piv->QAResult_[0] >= 0.0)
		{
			tString1.clear();
			for(int i=0;i<QARESULTSSZ;++i)
			{
				PRINTINTOBUFFERFUNCT(buffer,BUFFERMAXCHAR,L":%4.6g",piv->QAResult_[i]);
				tString1 += std::wstring(buffer);
			}
			tString1 = tString1.substr(1);
		}
		else
		{tString1 = std::wstring(L"####");}

		PRINTINTOBUFFERFUNCT
			(buffer,BUFFERMAXCHAR,PEAKSOUTFILELNFORMAT,
			//1
			ObjectID_,
			//2
			piv->PatchNumber_,
			//3
			static_cast<int>(piv->PK_BayesDetectStat_),
			//4
			piv->GaussianIndex_,
			//5
			piv->DetectionSigma_,
			//6
			piv->SrcAmplNormalised_,
			//7
			piv->UsedSigmaThreshold_,
			//8
			piv->JF_UsedSigmaThreshold_,
			//9
			piv->SrcFlux_,
			//10
			piv->JF_SrcFlux_.Mode_,
			//11
			piv->JF_lnRhoTh_,
			//12
			piv->JF_lnRho_,
			//13
			piv->JF_lnEvidence_,
			//14
			piv->JF_lnEvidenceErrBar_,
			//15
			piv->JF_lnFormFactor_,
			//16
			piv->JF_lnModelRatio_,
			//17
			piv->Pos_.YCoord_,
			//18
			piv->Pos_.XCoord_,
			//19
			piv->Pos_.JF_YCoord_.Mode_,
			//20
			piv->Pos_.JF_XCoord_.Mode_,
			//21
			piv->Pos_.JF_YCoord_.Mean_,
			//22
			piv->Pos_.JF_XCoord_.Mean_,
			//23
			(PIOVER2 - piv->GalPt_Colat_) * RAD2DEGREE,
			//24
			piv->GalPt_Long_  * RAD2DEGREE,
			//25
			piv->SrcFlux_mJys_,
			//26
			piv->SrcCompt_arcmin2_,
			//27
			piv->RealParams_.RealScale_,
			//28
			piv->JF_Radius_.Mode_,
			//29
			piv->Odds_,
			//30
			piv->ISNR2_,
			//31
			piv->CollListIndex_,
			//32
			piv->ErrorBars_.FluxErrorBar_,
			//33
			piv->ErrorBars_.TotalPosErrorBar_,
			//34
			piv->ErrorBars_.RadiusErrorBar_,
			//35
			piv->y0_,
			//36
			piv->JF_SrcFlux_.Mean_,
			//37
			piv->JF_Radius_.Mean_,
			//38
			piv->CoordPt_Lat_,
			//39
			piv->CoordPt_Long_,
			//40
			piv->EquPt_Lat_,
			//41
			piv->EquPt_Long_,
			//42
			piv->ErrorBars_.LowFluxErrorBar_,
			//43
			piv->ErrorBars_.HighFluxErrorBar_,
			//44
			piv->ErrorBars_.LowRadiusErrorBar_,
			//45
			piv->ErrorBars_.HighRadiusErrorBar_,
			//46
			piv->ErrorBars_.DegenCorr_,
			//47
			piv->ErrorBars_.DegenOrd_,
			//48
			piv->ErrorBars_.DegenOrdYErr_,
			//49
			piv->ErrorBars_.DegenSlope_,
			//50
			piv->ErrorBars_.DegenSlopeYErr_,
			//51
			tString.c_str(),
			//52
			tString1.c_str()
			);
		stream_.WriteLn(std::wstring(buffer));
	}
	return static_cast<int>(end-data.begin());
}
//
int			ReadObjResultsTxtFile::do_parseLn(const std::wstring& str)
{
	int						NItems;
	wchar_t					buffer[BUFFERMAXCHAR];
	wchar_t					buffer1[BUFFERMAXCHAR];
	std::vector<double>		ScalesColl;

	PeakCollReadbleType::StorageType::value_type	temp;
	
	if((NItems = swscanf(str.c_str(),PEAKSINFILELNFORMAT,
			&(temp.DetectID_),
			&(temp.PatchNumber_),
			&(temp.PK_BayesDetectStat_),
			&(temp.GaussianIndex_),
			&(temp.DetectionSigma_),
			&(temp.SrcAmplNormalised_),
			&(temp.UsedSigmaThreshold_),
			&(temp.JF_UsedSigmaThreshold_),
			&(temp.SrcFlux_),
			&(temp.JF_SrcFlux_.Mode_),
			&(temp.JF_lnRhoTh_),
			&(temp.JF_lnRho_),
			&(temp.JF_lnEvidence_),
			&(temp.JF_lnEvidenceErrBar_),
			&(temp.JF_lnFormFactor_),
			&(temp.JF_lnModelRatio_),
			&(temp.Pos_.YCoord_),
			&(temp.Pos_.XCoord_),
			&(temp.Pos_.JF_YCoord_.Mode_),
			&(temp.Pos_.JF_XCoord_.Mode_),
			&(temp.Pos_.JF_YCoord_.Mean_),
			&(temp.Pos_.JF_XCoord_.Mean_),
			&(temp.GalPt_Colat_),
			&(temp.GalPt_Long_),
			&(temp.SrcFlux_mJys_),
			&(temp.SrcCompt_arcmin2_),
			&(temp.RealParams_.RealScale_),
			&(temp.JF_Radius_.Mode_),
			&(temp.Odds_),
			&(temp.ISNR2_),
			&(temp.CollListIndex_),
			&(temp.ErrorBars_.FluxErrorBar_),
			&(temp.ErrorBars_.TotalPosErrorBar_),
			&(temp.ErrorBars_.RadiusErrorBar_),
			&(temp.y0_),
			&(temp.JF_SrcFlux_.Mean_),
			&(temp.JF_Radius_.Mean_),
			&(temp.CoordPt_Lat_),
			&(temp.CoordPt_Long_),
			&(temp.EquPt_Lat_),
			&(temp.EquPt_Long_),
			&(temp.ErrorBars_.LowFluxErrorBar_),
			&(temp.ErrorBars_.HighFluxErrorBar_),
			&(temp.ErrorBars_.LowRadiusErrorBar_),
			&(temp.ErrorBars_.HighRadiusErrorBar_),
			&(temp.ErrorBars_.DegenCorr_),
			&(temp.ErrorBars_.DegenOrd_),
			&(temp.ErrorBars_.DegenOrdYErr_),
			&(temp.ErrorBars_.DegenSlope_),
			&(temp.ErrorBars_.DegenSlopeYErr_),
			buffer,
			buffer1
		)) != PEAKS_NFIELDS)
		return -1;

	temp.ScaleLikeNoise_.clear();

	if(
		Zeus::GetNumbersHomogenColl(std::wstring(buffer),std::wstring(L":"),ScalesColl) &&
		(!(ScalesColl.size() % 3))
	)
	{
		ScaleLikeNoiseColl::value_type			t;
		const int								N_Ele(ScalesColl.size());
		for(int i=0;i<N_Ele;i+=3)
		{
			t.Scale_	= ScalesColl[i];
			t.Like_		= ScalesColl[i+1];
			t.Noise_	= ScalesColl[i+2];

			temp.ScaleLikeNoise_.push_back(t);
		}
	
	}
	ScalesColl.clear();
	if(
		Zeus::GetNumbersHomogenColl(std::wstring(buffer1),std::wstring(L":"),ScalesColl) &&
		(ScalesColl.size() == QARESULTSSZ)
	)
	{
		for(int i=0;i<QARESULTSSZ;++i)
		{
			temp.QAResult_[i] = ScalesColl[i];
		}
	}

	temp.GalPt_Colat_	= (90.0 - temp.GalPt_Colat_) / RAD2DEGREE;
	temp.GalPt_Long_	/= RAD2DEGREE;

	data_.Storage_.push_back(temp);
	return NItems;
}
//

}//namespace Zeus

