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
 *  Copyright (C) 2005, 2013, Pedro Carvalho
 *  Author: Pedro Carvalho
 */
#include "ZEUS_Priors.h"
#include "ZEUS_Priors.inc"

namespace Zeus
{

const double	GaussianPrior::Min_ = -GAUSSIANLIMIT;
const double	GaussianPrior::Max_	= GAUSSIANLIMIT;

//
Priors::~Priors(void)
{}
//
double	GaussianPrior::Approx(double r)
{	
	std::vector<double>::const_iterator	uBound(std::lower_bound(SampleValues_->begin(),SampleValues_->end(),r));
	if(uBound == SampleValues_->end())		return Max_;
	if(uBound == SampleValues_->begin())	return Min_;

	const	double	v0(*(uBound-1));
	const	double	d(*uBound - v0);
	const	double	offset(static_cast<double>(static_cast<int>(uBound - SampleValues_->begin())));
	const	double	Value((offset * Delta_) + Min_);

	return	LinearInterpolation(Value - Delta_,Value,(r - v0)/d);
}
//
void	RadiusNinfPriorLinInterp::FillSamples(std::vector<dPoint>& Samples)
{
	MakeX(Samples);
	FillCDF(Samples);
}
//
void	RadiusNinfPriorLinInterp::MakeX(std::vector<dPoint>& Samples)
{
	dPoint				t;

	Samples.clear();
	t.x_ = min_;
	Samples.push_back(t);
	t.x_ = max_;
	Samples.push_back(t);

	// k / SQRT(2)
	t.x_ = param_/SQRT2;
	if((t.x_ > min_) && (t.x_ < max_))
		Samples.push_back(t);

	// 2·SQRT(3)/(9·k^2)
	const double Max_y(0.3849001794595493 / (param_ * param_));
	const double delta(Max_y/static_cast<double>(NSamples_));

	for(double cur_y=delta;cur_y < Max_y;cur_y += delta)
	{
		t.x_ = LowerBranch(cur_y);
		if((t.x_ > min_) && (t.x_ < max_))
			Samples.push_back(t);

		t.x_ = HigherBranch(cur_y);
		if((t.x_ > min_) && (t.x_ < max_))
			Samples.push_back(t);		
	}
}
//
double	RadiusNinfPriorLinInterp::HigherBranch(double y)
{
	// 0.7598356856509647
	// 3^(-1/4)
	//2.598076211352566
	// 3*SQRT(3) / 2
	//1.732050807567972
	// SQRT(3)
	const double t(param_*param_*y);
	return 0.7598356856509647 * std::sqrt((2.0*std::cos(std::acos(-2.598076211352566*t)/3.0))-1.732050807567972*t)/std::sqrt(y);
}
//
double	RadiusNinfPriorLinInterp::LowerBranch(double y)
{
	// 0.7598356856509647
	// 3^(-1/4)
	//2.598076211352566
	// 3*SQRT(3) / 2
	//1.732050807567972
	// SQRT(3)
	const double t(param_*param_*y);
	return 0.7598356856509647 * std::sqrt((2.0*std::sin(std::asin(2.598076211352566*t)/3.0))-1.732050807567972*t)/std::sqrt(y);
}
//
void	RadiusNinfPriorLinInterp::FillCDF(std::vector<dPoint>& Samples)
{
	std::vector<dPoint>::iterator			piv(Samples.begin());
	std::vector<dPoint>::const_iterator		const end(Samples.end());
	const double	Rmax(std::sqrt(max_*max_ + param_*param_));
	const double	Rmin(std::sqrt(min_*min_ + param_*param_));

	for(;piv != end;++piv)
	{
		const double tx(std::sqrt(piv->x_ * piv->x_ + param_*param_));
		piv->cdf_ = ((Rmax * (tx - Rmin)) / (tx * (Rmax - Rmin)));
	}
}
//
void	GaussianPrior::FillInSampleValues(void)
{
	
	GammaFuncts		tGammaf;
	double	t(0);
	for(int i=0;i<=Samples_;++i,t +=Delta_)
	{
		SampleValues_->at(i) = (1.0 + tGammaf.erf((Min_ + t)/SQRT2))/2.0;
	}

	std::sort(SampleValues_->begin(),SampleValues_->end());
}
//
void	MarginalPrior::do_SetParams(const PriorParamsType& params,double& Min,double& Max)
{
	Pdf_.clear();
	if((params.size()<4) || (params.size() % 2))
		errWrongArgs((params.size()<4)?4:(params.size() >> 1)<<1,params.size());
	const int dSZ(params.size() >> 1);
	if(!(dSZ % 2))
		errWrongArgs((dSZ >> 1)<<1,dSZ);
	for(int i=0;i<(params.size()-2);i+=2)
	{Pdf_.push_back(dPoint(params[i],params[i+3]));}
	std::sort(Pdf_.begin(),Pdf_.end(),SortPDFbyValue);
	Min = Min_ = params[1];
	Max = Pdf_[Pdf_.size()-1].x_;
}
//
void	MarginalPrior::FillSamples(std::vector<dPoint>& Samples)
{
	double Acc(0.0),prev(Min_);
	double t;
	Samples.clear();
	std::vector<dPoint>::const_iterator	piv(Pdf_.begin());
	std::vector<dPoint>::const_iterator	const end(Pdf_.end());
	for(;piv != end;++piv)
	{
		t=(piv->cdf_*(piv->x_-prev));Acc += t;
		Samples.push_back(dPoint(Acc,piv->x_));
		prev=piv->x_;
	}
	if(Acc != 1.0)
	{
		Acc = 1.0 / Acc;
		std::vector<dPoint>::iterator				pivCDF(Samples.begin());
		std::vector<dPoint>::const_iterator	const	endCDF(Samples.end());
		for(;pivCDF != endCDF;++pivCDF)
		{pivCDF->cdf_ *= Acc;}
	}
	Pdf_.clear();
}
//
void	ConditionalPrior::mySetParams(const PriorParamsType& Cdf,const PriorParamsType& Pdf)
{
	ConditionalPrior::SetParams(Cdf);
// Now the PDF
	if(Pdf.size() != PriorsCondPDFSIZE)
		errWrongArgs(PriorsCondPDFSIZE,Pdf.size());
	logX_ = Pdf[0];
	logY_ = Pdf[1];
	PdfA_ = Pdf[2];
	const double Ang(Pdf[3] /RAD2DEGREE);
	a_ = (std::cos(Ang)*std::cos(Ang))/(2.0*Pdf[4]*Pdf[4])+(std::sin(Ang)*std::sin(Ang))/(2.0*Pdf[5]*Pdf[5]);
	b_ = -std::sin(2.0*Ang)/(4.0*Pdf[4]*Pdf[4])+std::sin(2.0*Ang)/(4.0*Pdf[5]*Pdf[5]);
	c_ = (std::sin(Ang)*std::sin(Ang))/(2.0*Pdf[4]*Pdf[4])+(std::cos(Ang)*std::cos(Ang))/(2.0*Pdf[5]*Pdf[5]);
}
//
void	ConditionalPrior::SetParams(const PriorParamsType& Cdf)
{
	const unsigned int tSZ(Distribs_.size());
	const unsigned int minSZ(tSZ<<2);
	const unsigned int condSZ(Cdf.size() / tSZ);

	if((Cdf.size()<minSZ) || (Cdf.size() % 2))
		errWrongArgs((Cdf.size()<minSZ)?minSZ:(Cdf.size() >> 1)<<1,Cdf.size());
	for(unsigned int i=0;i<tSZ;++i)
	{delete Distribs_[i];}
	currConditional_ = 0;
	currPdf_ = Zeus::pwsNaN;
	for(unsigned int i=0;i<tSZ;++i)
	{Distribs_[i] = new MarginalPrior();}

	for(unsigned int i=0;i<tSZ;++i)
	{
		{
			PriorParamsType dummy(Cdf.begin()+(i*condSZ),(Cdf.begin()+condSZ)+(i*condSZ));
			(Distribs_[i])->SetParams(dummy);
		}
	}
}
//
} // namespace Zeus

