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

/*			
			||
			((peak.GalPt_Colat_ > GalacticCutLow_) && (peak.GalPt_Colat_ < GalacticCutHigh_) && (peak.SrcAmplNormalised_ < MinimumGalacticSigma_))			||

			((peak.GalPt_Colat_  >= MGCLLARGECOLATLOW) && (peak.GalPt_Colat_  < MGCLLARGECOLATHIGH) && (peak.GalPt_Long_   >= MGCLLARGELONGLOW)
			&& (peak.GalPt_Long_   < MGCLLARGELONGHIGH) && (peak.SrcAmplNormalised_ < MinimumGalacticSigma_)) ||

			((peak.GalPt_Colat_  >= MGCLSMALLCOLATLOW) && (peak.GalPt_Colat_  < MGCLSMALLCOLATHIGH) && (peak.GalPt_Long_   >= MGCLSMALLLONGLOW)
			&& (peak.GalPt_Long_   < MGCLSMALLLONGHIGH) && (peak.SrcAmplNormalised_ < MinimumGalacticSigma_))
*/


//---------------------------------------
#ifndef COMPRESSCATALOGH
#define  COMPRESSCATALOGH

#include "PWS_Globals.h"
#include "PWS_GlobalInfoStore.h"

//---------------------------------------

#define MGCLLARGECOLATLOW	2.059488517
#define MGCLLARGECOLATHIGH	2.234021443
#define MGCLLARGELONGLOW	4.78220215
#define MGCLLARGELONGHIGH	4.956735076

#define MGCLSMALLCOLATLOW	2.303834613
#define MGCLSMALLCOLATHIGH	2.391101075
#define MGCLSMALLLONGLOW	5.235987756
#define MGCLSMALLLONGHIGH	5.323254219
#define SIGMAMAXSCALE		7.0
#define RADIUSMINSCALE		0.10


#ifdef WIN32
#define INTERMEDCATNAMEFILE		L"Read %d lines from catalogue %ls\n"
#else
#define INTERMEDCATNAMEFILE		L"Read %d lines from catalogue %S\n"
#endif



inline bool SortCatBySNR(const Zeus::CatLineType& first,const Zeus::CatLineType& sec)
{
	return first.NormalAmpl_ > sec.NormalAmpl_;
}

inline bool SortCatBylnRho(const Zeus::CatLineType& first,const Zeus::CatLineType& sec)
{
	return first.lnRho_ > sec.lnRho_;
}


inline bool SortCatByPatch(const Zeus::CatLineType& first,const Zeus::CatLineType& sec)
{
	return first.Patch_ < sec.Patch_;
}
//
inline double GetGI(double d)
{
	return ((d > 1.0) ? (1.0 / d) : (d*d));
}
//
inline double GetRanklnRho(const Zeus::PeakType& d)
{
	double t(d.JF_lnRho_);

	if(d.PK_BayesDetectStat_ == Zeus::PeakType::PK_DET_BAYNOCONV) t /= 2.0;

	return t * ((d.JF_lnRho_ >= 0.0 )?GetGI(d.GaussianIndex_):1.0/GetGI(d.GaussianIndex_));
}
//
inline double GetRankNormAmpl(const Zeus::PeakType& d)
{
	return d.SrcAmplNormalised_ * GetGI(d.GaussianIndex_);
}
//
inline bool SortPeaksBylnRho(const Zeus::PeakCollType::iterator& first,const Zeus::PeakCollType::iterator& sec)
{
	double JF_lnRho1(first->JF_lnRho_);
	double JF_lnRho2(sec->JF_lnRho_);

	double fGI((first->GaussianIndex_ > 1.0) ? (1.0 /first->GaussianIndex_) : (first->GaussianIndex_*first->GaussianIndex_));
	if(first->PK_BayesDetectStat_ == Zeus::PeakType::PK_DET_BAYNOCONV) JF_lnRho1 /= 2.0;
	if(first->JF_lnRho_ < 0.0) fGI = 1.0 / fGI;
	fGI *= JF_lnRho1;

	double sGI((sec->GaussianIndex_ > 1.0) ? (1.0 /sec->GaussianIndex_) : (sec->GaussianIndex_*sec->GaussianIndex_));
	if(sec->PK_BayesDetectStat_ == Zeus::PeakType::PK_DET_BAYNOCONV) JF_lnRho2 /= 2.0;
	if(sec->JF_lnRho_ < 0.0) sGI = 1.0 / sGI;
	sGI *= JF_lnRho2;

	return fGI  >  sGI;

}
//
inline bool SortPeaksByNormAmpl(const Zeus::PeakCollType::iterator& first,const Zeus::PeakCollType::iterator& sec)
{
	const double fGI((first->GaussianIndex_ > 1.0) ? (1.0 /first->GaussianIndex_) : first->GaussianIndex_*first->GaussianIndex_);
	const double sGI((sec->GaussianIndex_ > 1.0) ? (1.0 /sec->GaussianIndex_) : sec->GaussianIndex_*sec->GaussianIndex_);

	return (first->SrcAmplNormalised_ * fGI)  >  (sec->SrcAmplNormalised_ * sGI);
}
//
inline bool SortPeaksByColat(const Zeus::PeakType& first,const Zeus::PeakType& sec)
{
	return first.GalPt_Colat_ < sec.GalPt_Colat_;
}
//
struct	CatalogueFilterFunctorNonBlind
{
	inline bool operator()(const Zeus::PeakType& peak) const
	{
/*
		if(AssessmentType_ && (peak.PK_BayesDetectStat_ <= Zeus::PeakType::PK_DET_ERROR))
		{
			return true;
		}
*/
		return false;
	}

	CatalogueFilterFunctorNonBlind(int AssessmentType)
		:AssessmentType_(AssessmentType)
	{}

	const int		AssessmentType_;
};


struct CatalogueFilterFunctor
{
	inline bool operator()(const Zeus::PeakType& peak) const
	{

		const double Rad(AssessmentType_?(JF_Estimator_?peak.JF_Radius_.Mean_:peak.JF_Radius_.Mode_):peak.RealParams_.RealScale_);

		if(	(peak.UnderGalMask_) ||
			((peak.GalPt_Colat_ > GalacticCutLow_) && (peak.GalPt_Colat_ < GalacticCutHigh_))
			)
			return true;

// hard constrains
		if(SZ_ && (GlbVars_.OutputLat_ >= 0))
		{
			// no 0 scale sources 
			if((GlbVars_.OutputLat_ & 0x01) && (peak.RealParams_.RealScale_ < GlbVars_.PriorSrcMinScale_))
				return true;
			
			const double	y0(peak.y0_ * peak.SrcFlux_ * FluxConvCteSZ_);
			// y0 not in use. peak.y0_ contains the conversion Cte

			if(AssessmentType_)
			{
				const int Selector(GlbVars_.OutputLat_ >> 1);

				if(Selector == 0)
				{
					if(
						(peak.JF_Radius_.Mean_ <= 0.0)
					||  (std::abs((GlbVars_.SrcMaxScale_ - peak.RealParams_.RealScale_)/ GlbVars_.SrcMaxScale_) < 0.02) 
					||	((peak.JF_Radius_.Mean_ >= 12.0) && (((peak.JF_Radius_.Mean_ - peak.RealParams_.RealScale_)/peak.JF_Radius_.Mean_)  >= 0.2))
					)
						return true;

					goto hardconstrains_end;
				}

				if(
					(peak.GaussianIndex_ < 0.975)
					)
					return true;

				if(Selector <= 1) goto hardconstrains_end;

				if(
					(GlbVars_.AssessmentKind_ && (peak.ErrorBars_.TotalPosErrorBar_ > 2.5))
				||	(peak.GaussianIndex_ < 0.98)
				||	(peak.SrcAmplNormalised_ > 70.0)
				)
					return true;
hardconstrains_end:
				;
			}
		}

		if(AssessmentType_)
		{
			if(
				(peak.PK_BayesDetectStat_ <= Zeus::PeakType::PK_DET_ERROR)
//				|| (peak.JF_lnRho_ < 0.0)
				)
				return true;
		}
		else
		{
			const double flux(peak.SrcFlux_ * (SZ_?FluxConvCteSZ_:FluxCteJys_));
			if(
				(peak.SrcAmplNormalised_ < MinimumSigma_)	||
				(flux < MinimumFluxLevel_)											
				)
				return true;
		}
		return false;
	}

	CatalogueFilterFunctor(void)
		:GlbVars_((PlanckInfo::Instance())->GetGlobalVars()),AssessmentType_(GlbVars_.AssessmentKind_),JF_Estimator_(GlbVars_.Jf_Estimator_),
		SZ_(GlbVars_.SZ_),MinimumFluxLevel_(GlbVars_.OutputFluxThreshold_),MinimumSigma_(GlbVars_.OutputSigma_),
		GalacticCutLow_(PIOVER2 - (GlbVars_.OutputGalCut_ / RAD2DEGREE)),GalacticCutHigh_(PIOVER2 + (GlbVars_.OutputGalCut_ / RAD2DEGREE)),
		MinimumGalacticSigma_(1.0e20),PixSz_(GlbVars_.PixSz_),MaxScale_(GlbVars_.SrcMaxScale_),
		FluxConvCteSZ_(GlbVars_.PixSz_ * GlbVars_.PixSz_  * 0.001 * SR2ARCMIN2 * GlbVars_.ProfParam_.FluxCalibCte_),
		FluxCteJys_(GlbVars_.PixSz_ * GlbVars_.PixSz_  * 1.0e9 * GlbVars_.ProfParam_.FluxCalibCte_)
	{
		if(SZ_)
		{
			MinimumGalacticSigma_ = 1.0e20;
		}
	}

	const GlobalScalarVarsType&		GlbVars_;
	const int		AssessmentType_;
	const int		JF_Estimator_;
	const int		SZ_;
	const double	MinimumFluxLevel_;
	const double	MinimumSigma_;
	const double	GalacticCutLow_;
	const double	GalacticCutHigh_;
	double			MinimumGalacticSigma_;
	const double	PixSz_;
	const double	MaxScale_;
	const double	FluxConvCteSZ_;
	const double	FluxCteJys_	;
};

class CompressCatalogue
{

public:
	CompressCatalogue(int ContextID,const std::vector<std::wstring>& PeakFileName,const std::wstring& CatName,Healpix_Map<HEALPIX_ATOM_PREC>& MaskMap)
		:ContextID_(ContextID),PeakFileName_(PeakFileName),CatName_(CatName),GlbVars_((PlanckInfo::Instance())->GetGlobalVars()),
		StaticInfoColl_((PlanckInfo::Instance())->GetPlanckStaticInfo()),GeomInfo_((PlanckInfo::Instance())->GetGeoProps()),MaskMap_(MaskMap),CatWriter_(0),CatWriterExt_(0)
	{}
	int						Initialize(void);
//
	void					SetOutputFields(double frequency,double SZ_spectralCte);
//
	int						Do_CompressCatalogue(void);
//
	int						WriteOutCatalogue(void);
//
	~CompressCatalogue(void)
	{
		delete CatWriter_;
		delete CatWriterExt_;
	}
//
private:
	double										GetPatchSpin(const Zeus::PatchGeomLineType&	patch);
	int											do_Compress(double MaxDist,double MaxDistSZLimit);
	double										EvalmaxDistSZ_WG5(double MinAggrDist,double MaxDistSZLimit,const Zeus::PeakType& peak);
	double										EvalmaxDistPS(double MaxDist);
	int											PurgeCatalogue(void);

	int		Merge1SrsDifPatch(Zeus::PeakCollType& Coll,const Zeus::PeakType& element,double MaxDist,int AssessmentType);
	Zeus::CatLineCollType::iterator	RemoveLowProbLns(void);
	void						RejectUnderMask(void);

	const GlobalScalarVarsType&					GlbVars_;
	const	PlanckStaticType&					StaticInfoColl_;
	const Zeus::PatchGeomType&					GeomInfo_;
	std::wstring								CatName_;
	const std::vector<std::wstring>&			PeakFileName_;
	int											ContextID_;
	Zeus::GenCollWriter<Zeus::CatalogueFormatType>	*CatWriter_;
	Zeus::GenCollWriter<Zeus::CatalogueFormatType>	*CatWriterExt_;
	Zeus::PeakCollType							PeakColl_;
	Zeus::CatLineCollType						Catalogue_;
	Healpix_Map<HEALPIX_ATOM_PREC>&				MaskMap_;

};


#endif //COMPRESSCATALOGH

