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

//----------------------------------

#ifndef ZEUS_PRIORSH
#define ZEUS_PRIORSH
//
//
#include <algorithm>
#include "ZEUS_InOut.h"
#include "ZEUS_Exceptions.h"
#include "ZEUS_PhysicsMath.h"
#include "ZEUS_WorkSpace.h"

#define PriorsCondPDFSIZE			6		
#define PriorsThetaMarginalSIZE		82
#define PriorsYcondThetaSIZE		3280

extern const double PriorsYcondTheta[PriorsYcondThetaSIZE];
extern const double PriorsThetaMarginal[PriorsThetaMarginalSIZE];
extern const double PriorsCondPDF[PriorsCondPDFSIZE];

namespace Zeus
{
//
struct dPoint
{
	double	cdf_;
	double	x_;

	dPoint(void)
		:cdf_(-1.0)
	{}

	dPoint(double cdf_v)
		:cdf_(cdf_v)
	{}

	dPoint(double cdf_v,double x)
		:cdf_(cdf_v),x_(x)
	{}

	inline bool operator<(const dPoint& rhs) const
	{return cdf_ < rhs.cdf_;}
};
//
inline bool SortPDFbyValue(const dPoint& lhs,const dPoint& rhs)
{return lhs.x_ < rhs.x_;}
//
class	Priors
{
	Priors(const Priors& rhs);
	Priors& operator=(const Priors&);
protected:
	inline void errWrongArgs(int Expected, int Delivered)
	{
		wchar_t	buffer[BUFFERMAXCHAR];
		PRINTINTOBUFFERFUNCT
			(buffer,BUFFERMAXCHAR,L"Prior. Wrong number of arguments. Expected -> %d, delivered -> %d\n",Expected,Delivered);

		throw libException(ERROR_COD_ZEUSNUMWRONGARGS,buffer,this);
	}

	virtual double	do_GetSample(double r) = 0;

public:
	typedef	std::vector<double>							PriorParamsType;
	typedef	enum{Uniform,Power,Exponential,Gaussian,RadiusNinf,Marginal,Conditional}	PriorsType;	

	Priors(void)
	{}
	inline	double		GetSample(double r)
	{
		if(r<0.0){r = 0.0;}
		else{if(r>=1.0) r = 1.0 - std::numeric_limits<double>::epsilon();}
		return do_GetSample(r);
	}
	virtual	void		SetParams(const PriorParamsType& params) = 0;
	virtual	void		GetLimits(double& min,double& max) = 0;
	virtual	PriorsType	GetPriorType(void) const = 0;
	virtual	double		GetPDF(double v) const = 0;
	inline	double		operator()(double r)
	{return GetSample(r);}
	virtual	~Priors(void) = 0;
};
//
typedef	std::vector<Zeus::Priors::PriorsType>		PriorsTypeCollType;
//
//
class	PowerlawPrior : public Priors
{
	double	bExponent_,Exponent_;
	double	g_;
	double	a_;
	double	min_,max_,pdfCal_;

	inline void	SetInternalValues(double min,double max,double exponent)
	{
		min_ = min;max_ = max;Exponent_= exponent;
		bExponent_		= 1.0 - exponent;
		g_				= 1.0 / (std::pow(max,bExponent_) - std::pow(min,bExponent_));
		a_				= g_ * std::pow(min,bExponent_);
		pdfCal_			= bExponent_/(std::pow(max,bExponent_) -  std::pow(min,bExponent_));
		bExponent_		= 1.0 / bExponent_;
	}
protected:
	virtual double	do_GetSample(double r)
	{return	(std::pow((a_ + r)/g_,bExponent_));}
public:
	PowerlawPrior(double min=1.0,double max=5.0,double exponent=2.0)
	{SetInternalValues(min,max,exponent);}
	inline virtual	void	SetParams(const PriorParamsType& params)
	{
		if(params.size() < 3) errWrongArgs(3,params.size());
		SetInternalValues(params[0],params[1],params[2]);
	}
	inline virtual	void	GetLimits(double& min,double& max)
	{min = min_;max = max_;}

	inline virtual	PriorsType	GetPriorType(void)  const
	{return Priors::Power;}

	inline virtual	double		GetPDF(double v) const
	{
		return ((v>=min_) && (v<max_))?(std::pow(v,-Exponent_)*pdfCal_):0.0;
	}

	virtual			~PowerlawPrior(void)
	{}
};
//
//
class	ExponencialPrior : public Priors
{
	double	a_;
	double	b_;
	double	min_;
	double	max_;
	double	exp_;
	double	pdfCal_;

	inline void	SetInternalValues(double min,double max,double exponent)
	{
		min_	= min;
		max_	= max;
		exp_	= exponent;
		a_		= exponent *	min;
		b_		= exponent *	max;
		pdfCal_	= (exponent * std::exp(b_ + a_)) /(std::exp(b_)-std::exp(a_));
	}
protected:
	virtual double	do_GetSample(double r)
	{
		const double lnX((r==0.0)?Zeus::logZERO:std::log(r));
		const double t(1.0 - ((r==0.0)?0.0:r));

		return	(a_ + b_ - Zeus::AddLog(lnX+a_,std::log(t)+b_)) / exp_;
	}
public:
	ExponencialPrior(double min=0.0,double max=5.0,double exponent=1.0)
	{SetInternalValues(min,max,exponent);}

	inline virtual	void	SetParams(const PriorParamsType& params)
	{
		if(params.size() < 3) errWrongArgs(3,params.size());
		SetInternalValues(params[0],params[1],params[2]);
	}

	inline virtual	void	GetLimits(double& min,double& max)
	{min = a_/exp_;max = b_/exp_;}

	inline virtual	PriorsType	GetPriorType(void)  const
	{return Priors::Exponential;}

	inline virtual	double		GetPDF(double v) const
	{return ((v>=min_) && (v<max_))?(std::exp(exp_ * v) * pdfCal_):0.0;}

	virtual			~ExponencialPrior(void)
	{}
};
//
//
class	UniformPrior : public Priors
{
	double	min_;
	double	max_;
	double	pdfCal_;

	inline void	SetInternalValues(double min,double max)
	{min_	= min;max_	= max;pdfCal_ = 1.0/(max_ - min_);}
protected:
	virtual double	do_GetSample(double r)
	{return	min_ + (max_ - min_) * r;}
public:
	UniformPrior(double min=0.0,double max=1.0)
	{SetInternalValues(min,max);}

	inline virtual	void	SetParams(const PriorParamsType& params)
	{
		if(params.size() < 2) errWrongArgs(2,params.size());
		SetInternalValues(params[0],params[1]);
	}
	inline virtual	void	GetLimits(double& min,double& max)
	{min = min_;max = max_;}

	inline virtual	PriorsType	GetPriorType(void)  const
	{return Priors::Uniform;}

	inline virtual	double		GetPDF(double v) const
	{return ((v>=min_) && (v<max_))?pdfCal_:0.0;}

	virtual			~UniformPrior(void)
	{}
};
//
//
class	GaussianPrior : public Priors
{
	static const double		Min_,Max_;
	std::vector<double>		*SampleValues_;
	double					mean_;
	double					sigma_;
	double					pdfCal_;
	const int				Samples_;
	const double			Delta_;
	
	inline	void	SetInternalValues(double mean,double sigma)
	{mean_	= mean;sigma_	= sigma;pdfCal_= 1.0/(SQRT2PI * sigma_);}
	inline	double	Rigorous(double r)
	{
		if(r <= 1.0e-16)	return Min_;
		const	double t(1.0 - r);
		if(t <= 1.0e-16)	return Max_;
       	return (SQRT2 * Dierfc(2.0 * t));
	}


	double	Approx(double r);

	void	FillInSampleValues(void);

protected:
	virtual double	do_GetSample(double r)
	{
		double	value((Samples_ <= 0)?Rigorous(r):Approx(r));
		return mean_ + sigma_ * value;
	}
public:
	GaussianPrior(double mean=0.0,double sigma=1.0,int Samples=1000)
		:Samples_(Samples),SampleValues_(0),Delta_((Max_-Min_)/((Samples > 0)?static_cast<double>(Samples):1.0))
	{
		SetInternalValues(mean,sigma);
		if(Samples > 0)
		{
			SampleValues_	= new std::vector<double>(Samples + 1);
			FillInSampleValues();
		}
	}

	inline virtual	void	SetParams(const PriorParamsType& params)
	{
		if(params.size() < 2) errWrongArgs(2,params.size());
		SetInternalValues(params[0],params[1]);
	}
	inline virtual	void	GetLimits(double& min,double& max)
	{max = Max_;min = Min_;}

	inline virtual	PriorsType	GetPriorType(void)  const
	{return Priors::Gaussian;}

	inline virtual	double		GetPDF(double v) const
	{
		const double chi((v - mean_)/sigma_);
		return ((chi>=Min_) && (chi<Max_))?(pdfCal_*std::exp(-0.5*chi*chi)):0.0;
	}

	virtual			~GaussianPrior(void)
	{delete SampleValues_;}
};
//
class	CDF_LinInterp : public Priors
{
protected:
//
	virtual void	FillSamples(std::vector<dPoint>& Samples) = 0;
//
	virtual void	do_SetParams(const PriorParamsType& params,double& Min,double& Max) = 0;
//
	inline virtual double	do_GetSample(double r)
	{
		std::vector<dPoint>::const_iterator	const UBound(std::lower_bound(probSamples_.begin(),probSamples_.end(),dPoint(r)));

		if(UBound == probSamples_.end())	return xMax_;

		const double MinCDF((UBound == probSamples_.begin())?0.0:(UBound-1)->cdf_);
		const double MinX((UBound == probSamples_.begin())?xMin_:(UBound-1)->x_);
		return	LinearInterpolation(MinX,UBound->x_,(r - MinCDF)/(UBound->cdf_ - MinCDF));
	}
//
	inline double	GetSampleIndex(double r,unsigned int& Index)
	{
		std::vector<dPoint>::const_iterator	const UBound(std::lower_bound(probSamples_.begin(),probSamples_.end(),dPoint(r)));

		if(UBound == probSamples_.end())
		{
			Index = (UBound-1)-probSamples_.begin();
			return xMax_;
		}
		Index = UBound-probSamples_.begin();

		const double MinCDF((UBound == probSamples_.begin())?0.0:(UBound-1)->cdf_);
		const double MinX((UBound == probSamples_.begin())?xMin_:(UBound-1)->x_);
		return	LinearInterpolation(MinX,UBound->x_,(r - MinCDF)/(UBound->cdf_ - MinCDF));
	}
//
public:
//
	CDF_LinInterp(void)
		:Priors()
	{}
//
	virtual	void	SetParams(const PriorParamsType& params)
	{
		do_SetParams(params,xMin_,xMax_);
		FillSamples(probSamples_);
		std::sort(probSamples_.begin(),probSamples_.end());
	}
//
	inline virtual	void			GetLimits(double& min,double& max)
	{max = xMax_;min = xMin_;}
//
	virtual		~CDF_LinInterp(void)
	{}

private:	
	std::vector<dPoint>		probSamples_;
	double					xMax_,xMin_;
};
//
class ConditionalPrior : public Priors
{
	friend class				MarginalPrior;
	CDF_LinInterp*				currConditional_;
	double						currPdf_;
	std::vector<CDF_LinInterp*> Distribs_;

	double						logX_;
	double						logY_;
	double						PdfA_;
	double						a_,b_,c_;
//
	inline double						My2D_Pdf(double y,double x) const
	{
		const double tx(x-logX_);
		const double ty(y-logY_);
		return  PdfA_*std::exp(-0.5*(a_*tx*tx+2.0*b_*tx*ty+c_*ty*ty));
	}
//
	virtual	void	mySetParams(const PriorParamsType& Cdf,const PriorParamsType& Pdf);
//
	inline void setCDFConstrains(unsigned int paramsIndex)
	{currConditional_ = ((paramsIndex<Distribs_.size())?Distribs_[paramsIndex]:0);}
//
	inline void setPDFVal(double x)
	{currPdf_ = std::log10(x);}
//
public:
	ConditionalPrior(const PriorParamsType& Cdf,const PriorParamsType& Pdf,int NBins)
		:Priors(),Distribs_(NBins,NULL),currConditional_(NULL),currPdf_(Zeus::pwsNaN)
	{mySetParams(Cdf,Pdf);}
//
	inline virtual double			do_GetSample(double r)
	{return currConditional_?currConditional_->GetSample(r):((r<=0.5)?Distribs_[0]->GetSample(r):Distribs_[Distribs_.size()-1]->GetSample(r));}
//
	virtual	void	SetParams(const PriorParamsType& Cdf);
//
	inline virtual	void	GetLimits(double& min,double& max)
	{
		min = max = Zeus::pwsNaN;
		if(currConditional_)
			currConditional_->GetLimits(min,max);
	}
//
	inline virtual	PriorsType	GetPriorType(void)  const
	{return Priors::Conditional;}
//
	inline virtual	double		GetPDF(double v) const
	{
		const double t(std::log10(v));
		return (Zeus::pwsIsNaN(currPdf_)  || Zeus::pwsIsNaN(t))?0.0:My2D_Pdf(t,currPdf_);
	}
//
	virtual		~ConditionalPrior(void)
	{
		for(unsigned int i=0;i<Distribs_.size();++i)
		{delete Distribs_[i];}
	}
};
//
class MarginalPrior : public CDF_LinInterp
{
	std::vector<dPoint>				Pdf_;
	double							Min_;
	ConditionalPrior				*slaveObj_;
//
protected:
//
	virtual void	FillSamples(std::vector<dPoint>& Samples);
//
	virtual void	do_SetParams(const PriorParamsType& params,double& Min,double& Max);
//
public:
	MarginalPrior(void)
		:CDF_LinInterp(),slaveObj_(NULL)
	{}
//
	MarginalPrior(const PriorParamsType& d,ConditionalPrior *slaveObj=NULL)
		:CDF_LinInterp(),slaveObj_(slaveObj)
	{CDF_LinInterp::SetParams(d);}
//
	inline ConditionalPrior* SetSlaveObject(ConditionalPrior* newObj)
	{
		std::swap(newObj,slaveObj_);
		return newObj;
	}
//
	inline virtual double			do_GetSample(double r)
	{
		unsigned int Index;
		const double v(CDF_LinInterp::GetSampleIndex(r,Index));
		if(slaveObj_)
			slaveObj_->setCDFConstrains(Index);
		return v;
	}
//
	inline virtual	PriorsType		GetPriorType(void)  const
	{return Priors::Marginal;}
//
	inline virtual	double			GetPDF(double v) const
	{
		if(slaveObj_)
			slaveObj_->setPDFVal(v);
		return 1.0;
	}
//
	virtual							~MarginalPrior(void)
	{}
};
//
class	RadiusNinfPriorLinInterp : public CDF_LinInterp
{
//
	double			min_;
	double			max_;
	double			param_;
	double			pdfCal_;
	int				NSamples_;
//
	void			MakeX(std::vector<dPoint>& Samples);
//
	void			FillCDF(std::vector<dPoint>& Samples);
//
	double			LowerBranch(double y);
//
	double			HigherBranch(double y);
//
protected:
//
	virtual void	FillSamples(std::vector<dPoint>& Samples);
//
	inline virtual void	do_SetParams(const PriorParamsType& params,double& Min,double& Max)
	{
		if(params.size() < 3) errWrongArgs(3,params.size());

		Min = min_	= params[0];
		Max = max_	= params[1];
		param_	= params[2];
		const double	a_(std::sqrt(min_*min_+param_*param_));
		const double	b_(std::sqrt(max_*max_+param_*param_));
		pdfCal_	= (a_ * b_) / (b_ - a_);
	}
//
public:
	RadiusNinfPriorLinInterp(void)
		:CDF_LinInterp()
	{}
//
	RadiusNinfPriorLinInterp(double param,double min,double max,int NSamples=500)
		:CDF_LinInterp(),NSamples_(NSamples)
	{
		PriorParamsType params;

		params.push_back(min);
		params.push_back(max);
		params.push_back(param);

		CDF_LinInterp::SetParams(params);
	}
//
	inline virtual	PriorsType	GetPriorType(void)  const
	{return Priors::RadiusNinf;}
//
	inline virtual	double		GetPDF(double v) const
	{
		if((v<min_) || (v>=max_)) return 0.0;

		const double t(v*v+param_*param_);
		return (v*pdfCal_)/(t * std::sqrt(t));
	}
//
	virtual	~RadiusNinfPriorLinInterp(void)
	{}
};
//
class	RadiusNinfPrior : public Priors
{
	double			min_;
	double			max_;
	double			param_;
	double			m_;
	double			M_;
	double			pdfCal_;

	inline void	SetInternalValues(double param,double min,double max)
	{
		min_=min;
		max_=max;
		param_=param;
		m_ = std::sqrt(min_*min_+param_*param_);
		M_ = std::sqrt(max_*max_+param_*param_);
		pdfCal_	= (m_ * M_) / (M_ - m_);
	}
protected:
	inline virtual double	do_GetSample(double r)
	{
		const double v((M_ * m_) / ((r*(m_ - M_)) + M_));
		return std::sqrt(v*v-param_*param_);
	}
public:
	RadiusNinfPrior(double param,double min,double max)
	{SetInternalValues(param,min,max);}

	inline virtual	void	SetParams(const PriorParamsType& params)
	{
		if(params.size() < 3) errWrongArgs(3,params.size());
		SetInternalValues(params[2], params[0], params[1]);
	}
//
	inline virtual	void	GetLimits(double& min,double& max)
	{max = max_; min = min_;}
//
	inline virtual	PriorsType	GetPriorType(void)  const
	{return Priors::RadiusNinf;}
//
	inline virtual	double		GetPDF(double v) const
	{
		if((v<min_) || (v>=max_)) return 0.0;

		const double t(v*v+param_*param_);
		// TODO
#ifdef CONTOURSONEOVERTHETAPRIOR
		return (v*pdfCal_)/(t);
#else
		return (v*pdfCal_)/(t * std::sqrt(t));// this is the right expression
#endif
	}

	virtual			~RadiusNinfPrior(void)
	{}
};
//
} // ZEUS namespace

#endif //ZEUS_PRIORSH

