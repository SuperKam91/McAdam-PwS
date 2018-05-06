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



//-------------------------------------------

#include <memory>
#include "PWS_CatalogueComp.h"
#include "ZEUS_InOut.h"

//-------------------------------------------
//
int		CompressCatalogue::Do_CompressCatalogue(void)
{
	PurgeCatalogue();
	int Nelem(do_Compress(GlbVars_.OutputMeltMaxDist_,GlbVars_.OutputMeltSZMaxLimit_));
	if(GlbVars_.AssessmentKind_ && !(GlbVars_.NonBlindDetection_))
	{
		Catalogue_.erase(RemoveLowProbLns(),Catalogue_.end());
	}
	if(GlbVars_.NonBlindDetection_)
	{
		std::sort(Catalogue_.begin(),Catalogue_.end(),SortCatByPatch);
	}
	else
	{
		std::sort(Catalogue_.begin(),Catalogue_.end(),GlbVars_.AssessmentKind_?SortCatBylnRho:SortCatBySNR);		
	}
	return Nelem;
}
//
int		CompressCatalogue::WriteOutCatalogue(void)
{
	wchar_t	buffer[BUFFERMAXCHAR];

	int	CatSz(static_cast<int>(Catalogue_.size()));
	
	Zeus::CatalogueFormatType tCat;

	tCat.Header_.SZ_params_			= GlbVars_.ProfParam_;

	if(
		((GlbVars_.NonBlindDetection_ % 4) == 3) &&
		(GlbVars_.SZ_ == 1) &&
		(GlbVars_.N_ObsPlanes_ != 1) &&
		(GlbVars_.Jf_Estimator_ == 0) &&
		(GlbVars_.AssessmentKind_ == 0)
		)
	{tCat.Header_.DetectionType_		=  1000;} // QA contours 
	else
	{tCat.Header_.DetectionType_		= GlbVars_.SZ_;}

	tCat.Header_.Estimator_			= GlbVars_.AssessmentKind_; //estimator
	tCat.Header_.PriorsType_		= 0; // Informative cosmological

	tCat.Header_.CollLstSz_			= GeomInfo_.Header_.CollListSz_;

	tCat.Storage_.swap(Catalogue_);

	if(CatWriter_)
	{
		CatWriter_->Write(tCat);
		CatWriter_->Flush();
		PRINTINTOBUFFERFUNCT
			(buffer,BUFFERMAXCHAR,L"Wrote %d lines into -> %ls\n",CatSz,(CatWriter_->GetCollID()).c_str());
		
		(Zeus::ConManager::Instance())->PrintStr2Console(buffer);

		delete CatWriter_;CatWriter_ = 0;
	}


	if(CatWriterExt_)
	{
		CatWriterExt_->Write(tCat);
		CatWriterExt_->Flush();
		PRINTINTOBUFFERFUNCT
			(buffer,BUFFERMAXCHAR,L"Wrote %d lines into -> %ls\n",CatSz,(CatWriterExt_->GetCollID()).c_str());
		
		(Zeus::ConManager::Instance())->PrintStr2Console(buffer);

		delete CatWriterExt_;CatWriterExt_ = 0;
	}



	return CatSz;
}

//
int		CompressCatalogue::Initialize(void)
{
	std::vector<std::wstring>::const_iterator	piv(PeakFileName_.begin());
	std::vector<std::wstring>::const_iterator	const end(PeakFileName_.end());
	PeakColl_.clear();
	Zeus::PeakCollReadbleType	Dummy;

	for(;piv != end; ++piv)
	{
		{
			std::wstring	ExcptStr;
			std::auto_ptr<Zeus::GenCollReader<Zeus::PeakCollReadbleType> > PeakReader(Zeus::GetObjReaderHandler(ContextID_,GlbVars_.DirIn_ , *piv));

			try
			{
				PeakReader->Initialize();
				PeakReader->Read();
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
				(Zeus::ConManager::Instance())->PrintStr2Console(std::wstring(L"Cannot find/read file/object -> ") + PeakReader->GetCollID() + std::wstring(L"\n\n"));
				continue;
			}

			PeakReader->Release(Dummy);
			PeakColl_.insert(PeakColl_.end(),Dummy.Storage_.begin(),Dummy.Storage_.end());

			wchar_t	buffer[BUFFERMAXCHAR];
			swprintf(buffer,BUFFERMAXCHAR,INTERMEDCATNAMEFILE,(int)Dummy.Storage_.size(),(PeakReader->GetCollID()).c_str());
			(Zeus::ConManager::Instance())->PrintStr2Console(buffer);
			Dummy.Storage_.clear();
		}
	}

	CatWriter_ = Zeus::GetCatWriterHandler(ContextID_,CatName_,GlbVars_.DirOut_,Zeus::CatalogueFormatType::HeaderType::CatNColumns_);

	if(CatWriter_->Remove())
	{
		(Zeus::ConManager::Instance())->PrintStr2Console(std::wstring(L"Could not find/delete this object -> ") + CatWriter_->GetCollID());
	}
	else
	{
		(Zeus::ConManager::Instance())->PrintStr2Console(std::wstring(L"Object successfully deleted -> ") + CatWriter_->GetCollID());
	}


#if (defined(HFIDMC) || defined(LFIDPC)) && defined(HFIDMC_EXTOBJECTS)
	if(GlbVars_.Data2Buffer_)
	{
#ifdef WIN32
		std::wstring	CatsExt(L"Cats\\");
#else
		std::wstring	CatsExt(L"Cats/");
#endif
		std::wstring	FullPath(GlbVars_.DirBuffer_ + CatsExt);
		Zeus::CreateDir(FullPath);
		std::wstring	DirDummy;
		std::wstring	tFName(Zeus::ExtractFileName(CatName_,DirDummy));

		CatWriterExt_ = Zeus::GetCatWriterHandler(1000,tFName,FullPath,Zeus::CatalogueFormatType::HeaderType::CatNColumns_);

		if(CatWriterExt_->Remove())
		{
			(Zeus::ConManager::Instance())->PrintStr2Console(std::wstring(L"Could not find/delete this object -> ") + CatWriterExt_->GetCollID());
		}
		else
		{
			(Zeus::ConManager::Instance())->PrintStr2Console(std::wstring(L"Object successfully deleted -> ") + CatWriterExt_->GetCollID());
		}

		CatWriterExt_->Initialize();
	}
#endif

	return CatWriter_->Initialize();
}
//
double	CompressCatalogue::EvalmaxDistSZ_WG5(double MinAggrDist,double MaxDistSZLimit,const Zeus::PeakType& peak)
{
	double ClusterVRadius;
	
	const double Rad(GlbVars_.AssessmentKind_?(GlbVars_.Jf_Estimator_?peak.JF_Radius_.Mean_:peak.JF_Radius_.Mode_):peak.RealParams_.RealScale_);

	if(peak.ErrorBars_.TotalPosErrorBar_ > 0.0)
	{
		ClusterVRadius = peak.ErrorBars_.TotalPosErrorBar_ * 3.0;
	}
	else
	{ClusterVRadius = Rad * GlbVars_.ProfParam_.VirialRatio_;}

	if(ClusterVRadius < MinAggrDist) ClusterVRadius = MinAggrDist;

	if(ClusterVRadius >= MaxDistSZLimit) return MaxDistSZLimit;
	return ClusterVRadius;
}

//
double	CompressCatalogue::EvalmaxDistPS(double MaxDist)
{
	int freq(GlbVars_.FreqsColl_[0].freq_);
	PlanckStaticType::const_iterator	piv(StaticInfoColl_.begin());
	PlanckStaticType::const_iterator	const end(StaticInfoColl_.end());
	for(;(piv != end) && (piv->Freq_ != freq);++piv)
		;
	if(piv == end)
		throw Zeus::libException(ERRCOD_PWS_FREQNOTFOUND,ERRMSG_PWS_FREQNOTFOUND,L"CompressCatalogue::EvalmaxDistPS");

	double  PixSzArcMin(GlbVars_.PixSz_ * RAD2ARCMIN);
	return (PixSzArcMin * static_cast<double>(Zeus::toInt(((piv->AntFWHM_ * MaxDist) / (FWHM2SIGMA * PixSzArcMin)) + 1.0)));
}

//
int		CompressCatalogue::do_Compress(double MaxDist,double MaxDistSZLimit)
{

	const double		FluxConv((GlbVars_.SZ_?(0.001 * GlbVars_.PixSz_ * GlbVars_.PixSz_  * SR2ARCMIN2):
		(GlbVars_.PixSz_ * GlbVars_.PixSz_  * 1.0e9)) *  GlbVars_.ProfParam_.FluxCalibCte_);

	const bool	IsChiSqr((GlbVars_.SZ_ == 1) && (GlbVars_.N_ObsPlanes_ == 1) && ((GlbVars_.NonBlindDetection_%4) == 3) && (GlbVars_.AssessmentKind_ == 0));

	std::sort(PeakColl_.begin(),PeakColl_.end(),SortPeaksByColat);

	Zeus::CatLineType								tCatLine;
	std::vector<Zeus::PeakCollType::iterator>		tCatSortFlux_(PeakColl_.size());
	Zeus::PeakCollType::iterator					PeakCollPiv(PeakColl_.begin());
	std::vector<Zeus::PeakCollType::iterator>::iterator	tCatSortFluxPiv(tCatSortFlux_.begin());
	std::vector<Zeus::PeakCollType::iterator>::iterator	const tCatSortFluxEnd(tCatSortFlux_.end());	

	for(;tCatSortFluxPiv != tCatSortFluxEnd;++PeakCollPiv,++tCatSortFluxPiv)
	{
		*tCatSortFluxPiv = PeakCollPiv;
	}

	if(!(GlbVars_.NonBlindDetection_))
	{
		std::sort(tCatSortFlux_.begin(),tCatSortFlux_.end(),GlbVars_.AssessmentKind_?SortPeaksBylnRho:SortPeaksByNormAmpl);
	}

	Catalogue_.clear();

	std::vector<Zeus::PeakCollType::iterator>::const_iterator	piv(tCatSortFlux_.begin());
	std::vector<Zeus::PeakCollType::iterator>::const_iterator	const end(tCatSortFlux_.end());

	double MaxDistPS(-1.0);
	if(!(GlbVars_.SZ_))
	{
		MaxDistPS = EvalmaxDistPS(MaxDist);
	}

	int	include;

	for(;piv != end;++piv)
	{
		if((*piv)->Hited_)
			continue;

		if(!(GlbVars_.NonBlindDetection_))
		{
			include = Merge1SrsDifPatch(PeakColl_,*(*piv),(GlbVars_.SZ_?EvalmaxDistSZ_WG5(MaxDist,MaxDistSZLimit,*(*piv)):MaxDistPS),GlbVars_.AssessmentKind_);
		}
		else{include = 1;}

		if(include)
		{
			tCatLine.FluxCompt_						= (GlbVars_.Jf_Estimator_?(*piv)->JF_SrcFlux_.Mean_:(*piv)->JF_SrcFlux_.Mode_) * FluxConv;
			tCatLine.FluxComptGLRT_					= (*piv)->SrcFlux_ * FluxConv;
			tCatLine.GalLatDegs_					= (*piv)->GalPt_Colat_;
			tCatLine.GalLongDegs_					= (*piv)->GalPt_Long_;
			tCatLine.NormalAmpl_					= (*piv)->SrcAmplNormalised_;
			tCatLine.DetectSigma_					= (*piv)->DetectionSigma_ ;
			tCatLine.Radius_						= (GlbVars_.Jf_Estimator_?(*piv)->JF_Radius_.Mean_:(*piv)->JF_Radius_.Mode_) * GlbVars_.ProfParam_.RadiusCalCte_;
			tCatLine.RadiusGLRT_					= (*piv)->RealParams_.RealScale_ * GlbVars_.ProfParam_.RadiusCalCte_; 
			tCatLine.lnRho_							= (*piv)->JF_lnRho_;
			tCatLine.lnEvidence_					= (*piv)->JF_lnEvidence_;
			tCatLine.Patch_							= (*piv)->PatchNumber_;
			{
			const Zeus::PatchGeomLineType&	GeomProps(((PlanckInfo::Instance())->GetGeoProps()).Storage_.at(tCatLine.Patch_));

			tCatLine.PatchGalLatDegs_				= (PIOVER2 - GeomProps.X0Y0Ptg_.theta ) * RAD2DEGREE;
			tCatLine.PatchGalLongDegs_				= GeomProps.X0Y0Ptg_.phi * RAD2DEGREE;
//			tCatLine.PatchSpin_						= GetPatchSpin(GeomProps);
			tCatLine.PatchSpin_						= GeomProps.Spin_;

			}
			tCatLine.Gaussianity_					= (*piv)->GaussianIndex_;
			tCatLine.PatchMF_sigmaSqr_				= (IsChiSqr?1.0 / ((tCatLine.Gaussianity_ * tCatLine.Gaussianity_) * (*piv)->ISNR2_):-1.0);
			tCatLine.SZ_ConversionCte_				= (*piv)->y0_;
			tCatLine.lnPenaltySrc_					= ((*piv)->y0_ * tCatLine.FluxComptGLRT_);
			tCatLine.ErrorBars_.FluxErrorBar_		= (*piv)->ErrorBars_.FluxErrorBar_;
			tCatLine.ErrorBars_.TotalPosErrorBar_	= ((GlbVars_.AssessmentKind_)?(((*piv)->PK_BayesDetectStat_== Zeus::PeakType::PK_DET_BAYNOCONV)?-((*piv)->ErrorBars_.TotalPosErrorBar_):((*piv)->ErrorBars_.TotalPosErrorBar_)*2.0):-1.0);
			// position error bar multiplied by 2 to give 2 sigma error bars
			tCatLine.ErrorBars_.RadiusErrorBar_		= ((GlbVars_.AssessmentKind_)?(*piv)->ErrorBars_.RadiusErrorBar_:-1.0);
			tCatLine.ErrorBars_.LowFluxErrorBar_	= ((GlbVars_.AssessmentKind_)?(*piv)->ErrorBars_.LowFluxErrorBar_:-1.0);
			tCatLine.ErrorBars_.HighFluxErrorBar_	= ((GlbVars_.AssessmentKind_)?(*piv)->ErrorBars_.HighFluxErrorBar_:-1.0);
			tCatLine.ErrorBars_.LowRadiusErrorBar_	= ((GlbVars_.AssessmentKind_)?(*piv)->ErrorBars_.LowRadiusErrorBar_:-1.0);
			tCatLine.ErrorBars_.HighRadiusErrorBar_	= ((GlbVars_.AssessmentKind_)?(*piv)->ErrorBars_.HighRadiusErrorBar_:-1.0);
			tCatLine.ErrorBars_.DegenCorr_			= ((GlbVars_.AssessmentKind_)?(*piv)->ErrorBars_.DegenCorr_:-1.0);
			tCatLine.ErrorBars_.DegenOrd_			= ((GlbVars_.AssessmentKind_)?(*piv)->ErrorBars_.DegenOrd_:-1.0);
			tCatLine.ErrorBars_.DegenOrdYErr_		= ((GlbVars_.AssessmentKind_)?(*piv)->ErrorBars_.DegenOrdYErr_:-1.0);
			tCatLine.ErrorBars_.DegenSlope_			= ((GlbVars_.AssessmentKind_)?(*piv)->ErrorBars_.DegenSlope_:-1.0);
			tCatLine.ErrorBars_.DegenSlopeYErr_		= ((GlbVars_.AssessmentKind_)?(*piv)->ErrorBars_.DegenSlopeYErr_:-1.0);
			tCatLine.ScaleLikeNoise_				= (*piv)->ScaleLikeNoise_;
			tCatLine.CollLstIndex_					= Zeus::toInt(((*piv)->CollListIndex_) + 0.5);
			memcpy(tCatLine.QAResult_,(*piv)->QAResult_,QARESULTSSZ*sizeof(double));
			Catalogue_.push_back(tCatLine);
		}
	}

	return static_cast<int>(Catalogue_.size());

}
//
double		CompressCatalogue::GetPatchSpin(const Zeus::PatchGeomLineType&	patch)
{
	/*
	// Uncomment for general spin

	rotmatrix rot;
	pointing ptg(patch.X0Y0Ptg_);
	ptg.theta	-= PIOVER2;
	rot.Make_CPAC_Euler_Matrix(ptg.phi,ptg.theta,0.0);
	rot.Transpose();

	vec3		Vec0(patch.X0Y0_);
	vec3		Vec1(patch.XLY0_);

	Vec1 = Vec1 - Vec0;
	Vec1.x = 0.0;
	Vec1.Normalize();
	double ang((Vec1.z >= 0.0)?std::acos(Vec1.y):-(std::acos(Vec1.y));
	// if and less 1' then return 0.0;
	return (std::abs(ang) < 3.0e-4)?0.0:ang;
*/
	// WARNING
	// always assume iso-latitude patch cutter
	return 0.0;
}
//
int		CompressCatalogue::Merge1SrsDifPatch(Zeus::PeakCollType& Coll,const Zeus::PeakType& element,double MaxDist,int AssessmentType)
{
	Zeus::PeakType			temp;
	double					RankCriterion(AssessmentType?GetRanklnRho(element):GetRankNormAmpl(element));
	int						Include(1);
	
	MaxDist				/= RAD2ARCMIN;

	temp.GalPt_Long_	= element.GalPt_Long_;
	

	temp.GalPt_Colat_	= element.GalPt_Colat_ + MaxDist;
	if(temp.GalPt_Colat_ > PI)  temp.GalPt_Colat_ = PI;
	
	Zeus::PeakCollType::iterator	uBound(std::upper_bound(Coll.begin(),Coll.end(),temp,SortPeaksByColat));

	temp.GalPt_Colat_ = element.GalPt_Colat_ - MaxDist;
	if(temp.GalPt_Colat_ < 0)	temp.GalPt_Colat_ = 0;

	Zeus::PeakCollType::iterator	lBound(std::lower_bound(Coll.begin(),Coll.end(),temp,SortPeaksByColat));

	for(;lBound != uBound;++lBound)
	{
		if(MaxDist < Zeus::Dist2ptgs(element.GalPt_Colat_,element.GalPt_Long_,lBound->GalPt_Colat_,lBound->GalPt_Long_))
			continue;
		lBound->Hited_ = 1;
		if(Include && (RankCriterion < (AssessmentType?GetRanklnRho(*lBound):GetRankNormAmpl(*lBound))))
		{
			Include = 0;
		}
	}
	return Include;
}


//
void	CompressCatalogue::SetOutputFields(double frequency,double SZ_spectralCte)
{
	Zeus::CatLineCollType::iterator			piv(Catalogue_.begin());
	Zeus::CatLineCollType::const_iterator	const end(Catalogue_.end());

	Zeus::PlanckUnitsValuesTranform	UnitsTransfChi(Zeus::PlanckUnitsValuesTranform::BRIGHTNESS,Zeus::PlanckUnitsValuesTranform::THERMO_T,frequency,1.0/CMBTEMP);

	if(piv == end)
		return; //no data

	int			Index(0);
	const double	FluxConv((GlbVars_.SZ_?(0.001 * GlbVars_.PixSz_ * GlbVars_.PixSz_  * SR2ARCMIN2 *  GlbVars_.ProfParam_.FluxCalibCte_) / GlbVars_.ProfParam_.MNFW_Ratio_CY500CYR500_:
		(GlbVars_.PixSz_ * GlbVars_.PixSz_  * 1.0e9)));

	const bool	IsChiSqr((GlbVars_.SZ_ == 1) && (GlbVars_.N_ObsPlanes_ == 1) && ((GlbVars_.NonBlindDetection_%4) == 3) && (GlbVars_.AssessmentKind_ == 0));
	const		double	SZFlux2JB(IsChiSqr?UnitsTransfChi((SZ_spectralCte<0?-1000.0:1000.0)):-1.0e32);
	const		double	Compt2JB(UnitsTransfChi(1000.0*SZ_spectralCte));
	double		SZSigmaConv(UnitsTransfChi((GlbVars_.PixSz_ * GlbVars_.PixSz_  * SR2ARCMIN2)));
	SZSigmaConv	*= SZSigmaConv;


	for(;piv != end;++piv)
	{
		piv->ID_ = Index++;

		if(IsChiSqr)
		{
			piv->FluxCompt_				*= SZFlux2JB;
			piv->FluxComptGLRT_			= piv->FluxCompt_;
			piv->SZ_ConversionCte_		= Compt2JB;
			if(piv->PatchMF_sigmaSqr_ >= 0.0)
			{
				piv->PatchMF_sigmaSqr_ *= SZSigmaConv;
			}
		}
		else
		{
			if(!(piv->ScaleLikeNoise_.empty()))
			{
				Zeus::ScaleLikeNoiseColl::iterator			pivS(piv->ScaleLikeNoise_.begin());
				Zeus::ScaleLikeNoiseColl::const_iterator	const endS(piv->ScaleLikeNoise_.end());

				for(;pivS != endS; ++pivS)
				{
					pivS->Like_		*= FluxConv;
					pivS->Noise_	*= FluxConv;
					pivS->Scale_	*= GlbVars_.ProfParam_.MNFW_C500_;
				}
			}		
		}

		piv->GalLongDegs_	*=	RAD2DEGREE;
		piv->GalLatDegs_	=	90.0 - (piv->GalLatDegs_ * RAD2DEGREE);
	}
}


//
int		CompressCatalogue::PurgeCatalogue(void)
{
	if(!(GlbVars_.NonBlindDetection_) && (MaskMap_.Map().size()!=0))
	{
		RejectUnderMask();
	}

	Zeus::PeakCollType::iterator	newEnd;
	
	if(!(GlbVars_.NonBlindDetection_))
	{
		newEnd = std::remove_if(PeakColl_.begin(),PeakColl_.end(),CatalogueFilterFunctor());
	}
	else
	{
		newEnd = PeakColl_.end();
	}

	if(newEnd != PeakColl_.end())
	{PeakColl_.erase(newEnd,PeakColl_.end());}

	return static_cast<int>(PeakColl_.size());
}

//
void	CompressCatalogue::RejectUnderMask()
{
	Zeus::PeakCollType::iterator				piv(PeakColl_.begin());
	Zeus::PeakCollType::const_iterator		const end(PeakColl_.end());

	for(;piv != end;++piv)
	{
		if(MaskMap_[MaskMap_.ang2pix(pointing(piv->GalPt_Colat_, piv->GalPt_Long_))] < 0.5)
		{
			piv->UnderGalMask_ = 1;
		}
	}
}
//
Zeus::CatLineCollType::iterator	CompressCatalogue::RemoveLowProbLns(void)
{

	std::sort(Catalogue_.begin(),Catalogue_.end(),SortCatBylnRho);

	Zeus::CatLineCollType::iterator		piv(Catalogue_.begin());
	Zeus::CatLineCollType::iterator		pend(Catalogue_.end());
	double							prob(0.0);

	for(int i=1;piv != pend;++piv,++i)
	{
		prob += std::exp(-Zeus::AddLog(0.0,piv->lnRho_));
		if((static_cast<double>(i)*GlbVars_.OutputPurity_ / 100.0) < prob)
			break;
	}
	if(piv != pend) ++piv;
	return piv;
}
