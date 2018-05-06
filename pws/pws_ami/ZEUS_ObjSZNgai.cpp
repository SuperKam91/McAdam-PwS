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

#include "ZEUS_ObjSZNgai.h"

namespace Zeus
{
SrcObjSZNgai::SrcObjSZNgai(const std::wstring& name,double PeriodArc,int NSamples,const SZPS_ProfParamType&	SZPS_ProfParams,
	int ContextID,const std::wstring& DirInMasks,int sync_id,
	const SZ_ProfParamRangeType* SZPS_ProfParamsRange)
	:SrcObject(name,PeriodArc,NSamples),PeriodArc_(PeriodArc),YShift_(0.0),XShift_(0.0),
	SZPS_ProfParams_(SZPS_ProfParams), sync_id_(sync_id), SZPS_ProfParamsRange_(SZPS_ProfParamsRange), profileIndex_(-1)
{
#ifdef PWSMPI
	int		rank;
	int		rc;

	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	if(rank)
	{
		unsigned  int profSz(0);
		unsigned int dSz(NGAIPROFDATASIZE);
		if(SZPS_ProfParamsRange && SZPS_ProfParamsRange->IsInit())
		{
			const int		NBinsAlpha(SZPS_ProfParamsRange->alpha_.nbins_);
			const int		NBinsBeta(SZPS_ProfParamsRange->beta_.nbins_);
			profSz =		(NBinsAlpha*NBinsBeta);
			dSz	*=			(profSz+1);
		}

		Zeus::LArr1D<double>		myArr(dSz,0.0);

		if((rc = MPI_Bcast(myArr.begin(), myArr.getSz(), MPI_DOUBLE,0,MPI_COMM_WORLD))!=MPI_SUCCESS)
			throw Zeus::libMPIException(rc,L"MPI_Bcast Error: receive error",L"SrcObjSZNgai");
		// blocks until it receives OK from 0

		Zeus::LArr1D<double>::iterator	pivArr(myArr.begin());

		pivArr	+= 3;

		for(int i= 0;i<PIXTABLESZ;++i,++pivArr)
		{pixelValues_[i] = *pivArr;}
		for(int i= 0;i<TABLE0PIXSMALL;++i,++pivArr)
		{pixel0Small_[i] = *pivArr;}
		for(int i= 0;i<TABLE0PIXLARGE;++i,++pivArr)
		{pixel0Large_[i] = *pivArr;}

		if(profSz)
		{
			Profiles_.resize(profSz);
			for(int i=0;i<profSz;++i)
			{
				Profiles_[i]._alpha		= *pivArr++;
				Profiles_[i]._beta		= *pivArr++;
				Profiles_[i]._gamma		= *pivArr++;

				for(int j= 0;j<PIXTABLESZ;++j,++pivArr)
				{Profiles_[i].pixelValues_[j] = *pivArr;}
				for(int j= 0;j<TABLE0PIXSMALL;++j,++pivArr)
				{Profiles_[i].pixel0Small_[j] = *pivArr;}
				for(int j= 0;j<TABLE0PIXLARGE;++j,++pivArr)
				{Profiles_[i].pixel0Large_[j] = *pivArr;}
			}
		}
	}
	else
	{
		FillGNFWProfileArrays(SZPS_ProfParamsRange);
		StoreNewlyEvalCoefs(ContextID,DirInMasks);
	}
#else
	if(!ReadGNFWProfileCoefs(ContextID,DirInMasks))
	{
		FillGNFWProfileArrays(SZPS_ProfParamsRange);
		StoreNewlyEvalCoefs(ContextID,DirInMasks);
	}
#endif
	profileIndex_ = -1;
}
//
int		SrcObjSZNgai::ReadGNFWProfileCoefs(int ContextID,const std::wstring& DirInMasks)
{
	wchar_t			buffer[DOUBLETXTMAXSZ];
	std::wstring	ExcptStr;
	std::wstring	tName(GLOBAL_GNFWPROFILECOEFS_DEF);


	unsigned  int	NofProfiles(0);
	double			AlphaStep;
	double			BetaStep;
	int				NBinsAlpha(1);
	int				NBinsBeta(1);

	Profiles_.clear();
	if (SZPS_ProfParamsRange_ && SZPS_ProfParamsRange_->IsInit())
	{
		NBinsAlpha = SZPS_ProfParamsRange_->alpha_.nbins_;
		NBinsBeta = SZPS_ProfParamsRange_->beta_.nbins_;
		AlphaStep = (SZPS_ProfParamsRange_->alpha_.max_ - SZPS_ProfParamsRange_->alpha_.min_) / static_cast<double>((NBinsAlpha <= 1) ? 1 : NBinsAlpha - 1);
		BetaStep = (SZPS_ProfParamsRange_->beta_.max_ - SZPS_ProfParamsRange_->beta_.min_) / static_cast<double>((NBinsBeta <= 1) ? 1 : NBinsBeta - 1);
		NofProfiles =  (NBinsAlpha*NBinsBeta);
	}

	swprintf(buffer,DOUBLETXTMAXSZ,L"%d",sync_id_);

	tName += std::wstring(buffer);
	Zeus::LArr1D<double>	data;

	std::auto_ptr<Zeus::GenCollReader<Zeus::LArr1D<double> > >
		FReader(Zeus::GetVectFileHandlerReader(Loki::Type2Type<Zeus::LArr1D<double> >(),ContextID,
		tName,-1 /* if sz = -1 then get the full binary file*/,
		DirInMasks));

	try
	{
		FReader->Initialize();
		FReader->Read();
		FReader->Release(data);
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

	Zeus::LArr1D<double>::iterator	pivArr(data.begin());
	if(
		(*pivArr		!= SZPS_ProfParams_.MNFW_alpha_) ||
		(*(pivArr+1)	!= SZPS_ProfParams_.MNFW_beta_)	||
		(*(pivArr+2)	!= SZPS_ProfParams_.MNFW_gamma_) ||
		(*(pivArr+3)    != NofProfiles)
	)
		return 0;

	pivArr	+= 4;

	for(int i= 0;i<PIXTABLESZ;++i,++pivArr)
	{pixelValues_[i] = *pivArr;}

	for(int i= 0;i<TABLE0PIXSMALL;++i,++pivArr)
	{pixel0Small_[i] = *pivArr;}

	for(int i= 0;i<TABLE0PIXLARGE;++i,++pivArr)
	{pixel0Large_[i] = *pivArr;}

	if(NofProfiles)
	{
		Profiles_.resize(NBinsAlpha*NBinsBeta);
		for (int i = 0; i<NofProfiles; ++i)
		{
			Profiles_[i]._alpha = SZPS_ProfParamsRange_->alpha_.min_ + (AlphaStep * static_cast<double>(i / NBinsBeta));
			Profiles_[i]._beta = SZPS_ProfParamsRange_->beta_.min_ + (BetaStep * static_cast<double>(i % NBinsBeta));
			Profiles_[i]._gamma = SZPS_ProfParams_.MNFW_gamma_;

			if (
			(Profiles_[i]._alpha != *pivArr++) ||
			(Profiles_[i]._beta != *pivArr++) ||
			(Profiles_[i]._gamma != *pivArr++)
			)
				return 0;
			for (int j = 0; j<PIXTABLESZ; ++j, ++pivArr)
			{Profiles_[i].pixelValues_[j] = *pivArr;}
			for (int j = 0; j<TABLE0PIXSMALL; ++j, ++pivArr)
			{Profiles_[i].pixel0Small_[j] = *pivArr;}
			for (int j = 0; j<TABLE0PIXLARGE; ++j, ++pivArr)
			{Profiles_[i].pixel0Large_[j] = *pivArr;}
		}
	}

	return 1;
}
//
void	SrcObjSZNgai::StoreNewlyEvalCoefs(int ContextID,const std::wstring& DirInMasks)
{
	const unsigned int dSz(Profiles_.empty() ? NGAIPROFDATASIZE : NGAIPROFDATASIZE + ((Profiles_.size() * NGAIPROFVARDATASIZE)));

	Zeus::LArr1D<double>					myArr(dSz,0.0);
	Zeus::LArr1D<double>::iterator			pivArr(myArr.begin());	

	*pivArr++ = SZPS_ProfParams_.MNFW_alpha_;
	*pivArr++ = SZPS_ProfParams_.MNFW_beta_;
	*pivArr++ = SZPS_ProfParams_.MNFW_gamma_;
#ifndef PWSMPI
	*pivArr++ = static_cast<double>(Profiles_.size());
#endif
	for(int i= 0;i<PIXTABLESZ;++i,++pivArr)
	{*pivArr = pixelValues_[i];}
	for(int i= 0;i<TABLE0PIXSMALL;++i,++pivArr)
	{*pivArr = pixel0Small_[i];}
	for(int i= 0;i<TABLE0PIXLARGE;++i,++pivArr)
	{*pivArr = pixel0Large_[i];}
//
	for(int i=0;i<Profiles_.size();++i)
	{
		*pivArr++ = Profiles_[i]._alpha;
		*pivArr++ = Profiles_[i]._beta;
		*pivArr++ = Profiles_[i]._gamma;
		for(int j= 0;j<PIXTABLESZ;++j,++pivArr)
		{*pivArr = Profiles_[i].pixelValues_[j];}
		for(int j= 0;j<TABLE0PIXSMALL;++j,++pivArr)
		{*pivArr = Profiles_[i].pixel0Small_[j];}
		for(int j= 0;j<TABLE0PIXLARGE;++j,++pivArr)
		{*pivArr = Profiles_[i].pixel0Large_[j];}
	}
#ifdef PWSMPI
	// MPI Coefs ready to go.
	int rc;
	if((rc = MPI_Bcast(myArr.begin(), myArr.getSz(), MPI_DOUBLE,0,MPI_COMM_WORLD))!=MPI_SUCCESS)
		throw Zeus::libMPIException(rc,L"MPI_Bcast Error: send error",L"StoreNewlyEvalCoefs");
#else
	wchar_t			buffer[DOUBLETXTMAXSZ];
	std::wstring	tName(GLOBAL_GNFWPROFILECOEFS_DEF);
	swprintf(buffer,DOUBLETXTMAXSZ,L"%d",sync_id_);
	tName += std::wstring(buffer);

	std::auto_ptr<Zeus::GenCollWriter<Zeus::LArr1D<double> > > 
		FWriter(Zeus::GetVectFileHandlerWriter(Loki::Type2Type<Zeus::LArr1D<double> >(),ContextID,tName,DirInMasks));

	FWriter->Initialize();
	FWriter->Write(myArr);
	FWriter->Flush();
#endif //PWSMPI
}
//
void	SrcObjSZNgai::FillGNFWProfileArrays(const SZ_ProfParamRangeType* SZPS_ProfParamsRange)
{
	{
		GNFW_Profile	myProfile(SZPS_ProfParams_.MNFW_alpha_,SZPS_ProfParams_.MNFW_beta_,SZPS_ProfParams_.MNFW_gamma_);
		double			*pixelValuesPtr(pixelValues_);
		double			*pixel0SmallValuesPtr(pixel0Small_);
		double			*pixel0LargeValuesPtr(pixel0Large_);

		(Zeus::ConManager::Instance())->PrintStr2Console(L"Computing GNFW profile coeficients. Wait please.");

		myProfile.Initialise();

		for(int i=0;i<PIXTABLESZ;++i,++pixelValuesPtr)
		{*pixelValuesPtr = myProfile.GetPixValueAtRc(MINPIXDIST + (static_cast<double>(i) * SCALEPIX));}
	
		for(int i=0;i<TABLE0PIXSMALL;++i,++pixel0SmallValuesPtr)
		{*pixel0SmallValuesPtr = myProfile.Get0PixValueFast(MIN0PIXSMALL + (static_cast<double>(i) * SCALE0PIXSMALL));}

		for(int i=0;i<TABLE0PIXLARGE;++i,++pixel0LargeValuesPtr)
		{*pixel0LargeValuesPtr = myProfile.Get0PixValueFast(MIN0PIXLARGE + (static_cast<double>(i) * SCALE0PIXLARGE));}
	}

	if(SZPS_ProfParamsRange && SZPS_ProfParamsRange->IsInit())
	{
		const int		NBinsAlpha(SZPS_ProfParamsRange->alpha_.nbins_);
		const int		NBinsBeta(SZPS_ProfParamsRange->beta_.nbins_);
		const double	AlphaStep((SZPS_ProfParamsRange->alpha_.max_-SZPS_ProfParamsRange->alpha_.min_) / static_cast<double>((NBinsAlpha<=1)?1:NBinsAlpha-1));
		const double	BetaStep((SZPS_ProfParamsRange->beta_.max_-SZPS_ProfParamsRange->beta_.min_) / static_cast<double>((NBinsBeta<=1)?1:NBinsBeta-1));
		Profiles_.clear();
		Profiles_.resize(NBinsAlpha*NBinsBeta);
		for(int i=0;i<NBinsAlpha;++i)
		{
			for(int j=0;j<NBinsBeta;++j)
			{
				{
					GNFW_Profile	myProfile(SZPS_ProfParamsRange->alpha_.min_ + (AlphaStep * static_cast<double>(i)),
										SZPS_ProfParamsRange->beta_.min_ + (BetaStep * static_cast<double>(j)),SZPS_ProfParams_.MNFW_gamma_);

					myProfile.Initialise();

					Profiles_[NBinsBeta*i+j]._alpha = SZPS_ProfParamsRange->alpha_.min_ + (AlphaStep * static_cast<double>(i));
					Profiles_[NBinsBeta*i+j]._beta = SZPS_ProfParamsRange->beta_.min_ + (BetaStep * static_cast<double>(j));
					Profiles_[NBinsBeta*i+j]._gamma = SZPS_ProfParams_.MNFW_gamma_;

					double			*pixelValuesPtr(Profiles_[NBinsBeta*i+j].pixelValues_);
					double			*pixel0SmallValuesPtr(Profiles_[NBinsBeta*i+j].pixel0Small_);
			 		double			*pixel0LargeValuesPtr(Profiles_[NBinsBeta*i+j].pixel0Large_);

					for(int k=0;k<PIXTABLESZ;++k,++pixelValuesPtr)
					{*pixelValuesPtr = myProfile.GetPixValueAtRc(MINPIXDIST + (static_cast<double>(k) * SCALEPIX));}
	
					for(int k=0;k<TABLE0PIXSMALL;++k,++pixel0SmallValuesPtr)
					{*pixel0SmallValuesPtr = myProfile.Get0PixValueFast(MIN0PIXSMALL + (static_cast<double>(k) * SCALE0PIXSMALL));}

					for(int k=0;k<TABLE0PIXLARGE;++k,++pixel0LargeValuesPtr)
					{*pixel0LargeValuesPtr = myProfile.Get0PixValueFast(MIN0PIXLARGE + (static_cast<double>(k) * SCALE0PIXLARGE));}

					wprintf(L".");

				}
			}
		}
	}

	wprintf(L"\n");

	(Zeus::ConManager::Instance())->PrintStr2Console(L"Evaluation of the GNFW coeficients finished.");

	wprintf(L"\n");
}
//
double	SrcObjSZNgai::do_GetObjAtCoord(double YCoord,double XCoord) const
{
	YCoord -= YShift_ ; XCoord -= XShift_;
	if(Rs_ < OBJEPS)
	{
		if((std::abs(YCoord) < PeriodArc_) && (std::abs(XCoord) < PeriodArc_))
			return 1.0;
		else return 0.0;
	}
//
	const int	  tProfVar(IsVarProfile() && (profileIndex_ >=0));
	const double  VRatio(tProfVar ? (SZPS_ProfParams_.R500_Ratio_ *  filtP_.a3_) : SZPS_ProfParams_.VirialRatio_);
	const double dist(std::sqrt(YCoord*YCoord+XCoord*XCoord));
	const double distRc(dist / Rs_);
	const double MaxdistRc((VRatio <= MAXTRSDIST) ? VRatio : MAXTRSDIST);
	const int	 MaxLargePix0(MaxdistRc>MAX0PIXLARGE?MAX0PIXLARGE:MaxdistRc);
//
	if(distRc >= MaxdistRc)
		return 0.0;
//
	if(((std::abs(YCoord) < PeriodArc_) && (std::abs(XCoord) < PeriodArc_)) || (distRc < MINPIXDIST))
	{
		// this is the 0 pixel case
		double tdistRc(PeriodArc_  / (2.0 * Rs_));
		if(tdistRc >= MaxLargePix0)
			return 1;

		if(tdistRc > 2.0)
		{
			double t((tdistRc - MIN0PIXLARGE) / SCALE0PIXLARGE);
			int distIndex(Zeus::toInt(t));

			const double lowIndexV(tProfVar?Profiles_[profileIndex_].pixel0Large_[distIndex]:pixel0Large_[distIndex]);
			return ((((tProfVar?Profiles_[profileIndex_].pixel0Large_[distIndex + 1]:pixel0Large_[distIndex + 1]) - lowIndexV)) * (t - static_cast<double>(distIndex))) + lowIndexV;
		}
		else
		{
			if(tdistRc < MIN0PIXSMALL)
			{return (tProfVar?Profiles_[profileIndex_].pixel0Small_[0]:pixel0Small_[0]);}
			double t((tdistRc - MIN0PIXSMALL) / SCALE0PIXSMALL);
			int distIndex(Zeus::toInt(t));

			const double lowIndexV(tProfVar?Profiles_[profileIndex_].pixel0Small_[distIndex]:pixel0Small_[distIndex]);
			return ((((tProfVar?Profiles_[profileIndex_].pixel0Small_[distIndex + 1]:pixel0Small_[distIndex + 1]) - lowIndexV)) * (t - static_cast<double>(distIndex))) + lowIndexV;
		}
	}
	else
	{
		// All other pixels
		double t((distRc - MINPIXDIST) / SCALEPIX);
		int distIndex(Zeus::toInt(t));

		const double lowIndexV(tProfVar?Profiles_[profileIndex_].pixelValues_[distIndex]:pixelValues_[distIndex]);
		return ((((tProfVar?Profiles_[profileIndex_].pixelValues_[distIndex + 1]:pixelValues_[distIndex + 1]) - lowIndexV)) * (t - static_cast<double>(distIndex))) + lowIndexV;
	}
}
//
void	SrcObjSZNgai::do_GetObjBounds(double& YBound,double& XBound) const
{
	const int	  tProfVar(IsVarProfile() && (profileIndex_ >=0));
	const double  VRatio(tProfVar ? (SZPS_ProfParams_.R500_Ratio_ *  filtP_.a3_) : SZPS_ProfParams_.VirialRatio_);
	const double  MaxLenght(Rs_ * ((VRatio <= MAXTABLEDIST) ? VRatio : MAXTABLEDIST));
	
	if(PeriodArc_ >= MaxLenght)
	{
		YBound = XBound = PeriodArc_;
		YBound += std::abs(YShift_);
		XBound += std::abs(XShift_);
		return;
	}

	double zeroPixVal(do_GetObjAtCoord(0.0,0.0));

	if(zeroPixVal < OBJEPS)
	{
		YBound = XBound = PeriodArc_;
		YBound += std::abs(YShift_);
		XBound += std::abs(XShift_);
		return;
	}
	//
//	zeroPixVal *= boundTol_;
	zeroPixVal *= 1.0e-6;

	const int vectSize(sizeof(pixelValues_)/ sizeof(pixelValues_[0]));
	double currDistAbs(PeriodArc_);
	double currDistArray((currDistAbs - (MINPIXDIST * Rs_)) / (SCALEPIX * Rs_));
	const double DistArrayIncrem(PeriodArc_ / (SCALEPIX * Rs_));
	int	j;

	for(;(currDistAbs < MaxLenght);currDistAbs += PeriodArc_,currDistArray += DistArrayIncrem)
	{
		if((j = Zeus::toInt(currDistArray)) <= 0)
			continue;
		if((j >= vectSize) || ((tProfVar?Profiles_[profileIndex_].pixelValues_[j]:pixelValues_[j]) < zeroPixVal))
			break;
	}
			
	YBound = XBound = currDistAbs;
	YBound += std::abs(YShift_);
	XBound += std::abs(XShift_);
	return;
}
//
}



