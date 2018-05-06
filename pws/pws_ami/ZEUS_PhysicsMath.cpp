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

//-----------
#include <complex>
#include "ZEUS_Exceptions.h"
#include "ZEUS_InOut.h"
#include "ZEUS_PhysicsMath.h"
#include "ZEUS_WorkSpace.h"
#include "ZEUS_Debug.h"


namespace Zeus
{

const double	GammaFuncts::coef_[] =
{	
	((double) 76.18009172947146),
	((double)-86.50532032941677),
	((double)24.01409824083091),
	((double)-1.231739572450155),
	((double)0.1208650973866179e-2),
	((double)-0.5395239384953e-5)
};
//
const double	PlanckUnitsValuesTranform::PlanckChFreq_[] = 
{
	((double)100.0),
	((double)143.0),
	((double)217.0),
	((double)353.0),
	((double)545.0),
	((double)857.0),
// use theoretical values for the LFI
	((double)-1.0E32),
	((double)-1.0E32),
	((double)-1.0E32)
//
//	((double)30.0),
//	((double)44.0),
//	((double)70.0)
};
//
const double	PlanckUnitsValuesTranform::KCMB2KRJ[] = 
{
	((double)0.763644),
	((double)0.594176),
	((double)0.3109599),
	((double)0.0700971),
	((double)0.0059441),
	((double)0.00009673),
	((double)0.9787201766),
	((double)0.9504810765),
	((double)0.8783971570)
};
//
/*
const double	PlanckUnitsValuesTranform::MJYSR2KCMB[] = 
{
	((double)0.0041007),
	((double)0.00269099),
	((double)0.00206998),
	((double)0.00347924),
	((double)0.01723783),
	((double)0.44676215),
	((double)0.043485202),
	((double)0.017596543),
	((double)0.0074607856)
};
*/
//
const double	PlanckUnitsValuesTranform::MJYSR2KCMB[] = 
{
	((double)0.00409674814880241400065318552485),
	((double)0.00269010530309604757354384734449),
	((double)0.00206745108793204404051228295433),
	((double)0.003478844831483712709479511239),
	((double)0.01723081444442633294392841409893),
	((double)0.44089257467243224825849637467263),
	((double)0.043485202),
	((double)0.017596543),
	((double)0.0074607856)
};
//
const double	PlanckUnitsValuesTranform::MJYSR2KRJ[] =
{
	((double)0.0031315),
	((double)0.0015989),
	((double)0.0006437),
	((double)0.0002439),
	((double)0.0001025),
	((double)0.0000432),
	((double)pwsNaN),
	((double)pwsNaN),
	((double)pwsNaN)
};
//
const int	PlanckUnitsValuesTranform::PlanckChannelsN(sizeof(PlanckUnitsValuesTranform::PlanckChFreq_)/sizeof(double));

double Dierfc(double y)
{
	double infinity = 5.0;
	double qa = 9.16461398268964E-01,
      	qb = 2.31729200323405E-01, 
      	qc = 4.88826640273108E-01, 
      	qd = 1.24610454613712E-01, 
      	q0 = 4.99999303439796E-01, 
      	q1 = 1.16065025341614E-01, 
      	q2 = 1.50689047360223E-01, 
      	q3 = 2.69999308670029E-01, 
      	q4 = -7.28846765585675E-02,
      	pa = 3.97886080735226000E+00, 
      	pb = 1.20782237635245222E-01, 
      	p0 = 2.44044510593190935E-01, 
      	p1 = 4.34397492331430115E-01, 
      	p2 = 6.86265948274097816E-01, 
      	p3 = 9.56464974744799006E-01, 
      	p4 = 1.16374581931560831E+00, 
      	p5 = 1.21448730779995237E+00, 
      	p6 = 1.05375024970847138E+00, 
      	p7 = 7.13657635868730364E-01, 
      	p8 = 3.16847638520135944E-01, 
      	p9 = 1.47297938331485121E-02, 
      	p10 = -1.05872177941595488E-01, 
      	p11 = -7.43424357241784861E-02,
      	p12 = 2.20995927012179067E-03, 
      	p13 = 3.46494207789099922E-02, 
      	p14 = 1.42961988697898018E-02, 
      	p15 = -1.18598117047771104E-02, 
      	p16 = -1.12749169332504870E-02, 
      	p17 = 3.39721910367775861E-03, 
      	p18 = 6.85649426074558612E-03, 
      	p19 = -7.71708358954120939E-04, 
      	p20 = -3.51287146129100025E-03, 
      	p21 = 1.05739299623423047E-04, 
      	p22 = 1.12648096188977922E-03;
      
	if( y == 0.0 ) return infinity;
	double z, w, u, s, t, x;
	z = y;
	if( y > 1.0 ) z = 2.0 - y;
	w = qa - log( z );
	u = sqrt( w );
	s = ( qc + log( u ) ) / w;
	t = 1.0 / ( u + qb );
	x = u * ( 1.0 - s * ( 0.5 + s * qd ) ) - ( ( ( ( q4 * t + q3 ) * t + q2 ) * t + q1 ) * t + q0 ) * t;
	t = pa / ( pa + x );
	u = t - 0.5;
	s = ( ( ( ( ( ( ( ( ( p22 * u + p21 ) * u + p20 ) * u + p19 ) * u + p18 ) * u + p17 ) * u + p16 ) * 
		u + p15 ) * u + p14 ) * u + p13 ) * u + p12;
	s = ( ( ( ( ( ( ( ( ( ( ( ( s * u + p11 ) * u + p10 ) * u + p9 ) * u + p8 ) * u + p7 ) * u + p6 ) * 
		u + p5 ) * u + p4 ) * u + p3 ) * u + p2 ) * u + p1 ) * u + p0 ) * t - z * exp( x * x - pb );
      	x = x + s * ( 1.0 + x * s );
      	if( y > 1.0 ) x = - x;
	return x;
}

double  GammaFuncts::gcf(double a,double x) const
{
	double an,b,c,d,del,h,gln;

	gln=gammln(a);
	b=x+((double)1.0)-a;
	c=((double)1.0)/GAMMA_FPMIN;
	d=((double)1.0)/b;
	h=d;
	for(int i=1;i<=itMax_;++i)
	{
		an= -static_cast<double>(i)*(static_cast<double>(i)-a);
		b+= ((double)2.0);
		d=an*d+b;
		if(std::abs(d)<GAMMA_FPMIN) d=GAMMA_FPMIN;
		c=b+an/c;
		if(std::abs(c)<GAMMA_FPMIN) c=GAMMA_FPMIN;
		d=((double)1.0)/d;
		del=d*c;
		h *= del;
		if(std::abs(del-((double)1.0))<eps_)
			return std::exp(-x+a*std::log(x)-gln)*h;
	}
	if(UseExcpt_)
		throw libException(ERROR_COD_ZEUSFDIVERGES,ERROR_MSG_ZEUSFDIVERGES L"GammaFuncts::gcf",this);	
	return std::numeric_limits<double>::quiet_NaN();
}
//
double GammaFuncts::betacf(double a,double b,double x) const
{
	int m,m2;
	double aa,c,d,del,h,qab,qam,qap;

	qab=	a+b;
	qap=	a+1.0;
	qam=	a-1.0;
	c=		1.0;
	d=		1.0-(qab*x)/qap;
	if(std::fabs(d) < GAMMA_FPMIN)
		d=GAMMA_FPMIN;
	d=		1.0/d;
	h=		d;

	for (m=1;m<=itMax_;m++)
	{
		m2=	2*m;
		aa=	m*(b-m)*x/((qam+m2)*(a+m2));
		d=	1.0+aa*d;
		if (std::fabs(d) < GAMMA_FPMIN)
			d=GAMMA_FPMIN;
		c=	1.0+(aa/c);
		if (std::fabs(c) < GAMMA_FPMIN)
			c=GAMMA_FPMIN;
		d=	1.0/d;
		h *= (d*c);
		aa= -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
		d=	1.0+aa*d;
		if (std::fabs(d) < GAMMA_FPMIN)
			d=GAMMA_FPMIN;
		c=	1.0+(aa/c);
		if (std::fabs(c) < GAMMA_FPMIN)
			c=GAMMA_FPMIN;
		d=	1.0/d;
		del= d*c;
		h *= del;
		if (std::fabs(del-1.0) < eps_)
			break;
	}
	if(m > itMax_)
	{
		if(UseExcpt_)
			throw libException(ERROR_COD_ZEUSFDIVERGES,ERROR_MSG_ZEUSFDIVERGES L"GammaFuncts::betacf",this);	
		return std::numeric_limits<double>::quiet_NaN();
	}
	return h;
}
//
double GammaFuncts::gser(double a,double x) const
{
	double	sum,del,ap,gln;

	gln=gammln(a);

	if(x<= ((double)0.0))
	{
		if(UseExcpt_)
			throw libException(ERROR_COD_ZEUSOUTBOUNDS,ERROR_MSG_ZEUSOUTBOUNDS,L"Gamma::gser1");
		else x = GAMMA_FPMIN;
	}
		
	ap=a;
	del = sum = ((double)1.0/a);
	for(int n=1;n<itMax_;++n)
	{
		++ap;del *= x/ap;sum += del;
		if(std::abs(del)<std::abs(sum)*eps_)
		{return sum * std::exp(-x+a*std::log(x)-gln);}
	}
	if(UseExcpt_)
		throw libException(ERROR_COD_ZEUSFDIVERGES,ERROR_MSG_ZEUSFDIVERGES L"GammaFuncts::gser",L"Gamma::gser2");
	return std::numeric_limits<double>::quiet_NaN();
}
//
double GammaFuncts::gammln(double xx) const
{
	double x,y,tmp,ser;

	y = x = xx;
	tmp=x+((double)5.5);
	tmp -= (x+((double)0.5))*std::log(tmp);
	ser = ((double) 1.000000000190015);
	for(int j=0;j<=5;++j) ser += coef_[j]/++y;
	return -tmp + ((double)std::log((((double)2.5066282746310005) * ser)/x));
}
//
double	UnitsValueTranform::Set(UnitsType UnitIn,UnitsType UnitOut,double freqGHz,double conv)
{
	bool swap(false);

	if(UnitIn == UnitOut)
	{return (ConvCte_ = conv);}
	if(UnitIn < UnitOut)
	{
		UnitsType temp(UnitIn);
		UnitIn = UnitOut; UnitOut = temp;
		swap = true;
	}
	if(UnitIn == THERMO_T)
	{ConvCte_ =  Thermo2Antenna(freqGHz);}
	else
	{
		if(!UnitOut){ConvCte_ = BRIGHT2ANTENNA / (freqGHz*freqGHz);}
		else{ConvCte_ = Brightness2Thermo(freqGHz);}
	}
	if(swap) ConvCte_ = 1.0/ConvCte_;
	return (ConvCte_ *= conv);
}
//
double	PlanckUnitsValuesTranform::Set(UnitsType UnitIn,UnitsType UnitOut,double freqGHz,double conv)
{
	bool swap(false);

	if(UnitIn == UnitOut)
	{return (ConvCte_ = conv);}
	if(UnitIn < UnitOut)
	{
		UnitsType temp(UnitIn);
		UnitIn = UnitOut; UnitOut = temp;
		swap = true;
	}
	if(UnitIn == THERMO_T)
	{ConvCte_ =  Thermo2Antenna(freqGHz);}
	else
	{
		if(!UnitOut){ConvCte_ = Brightness2Antenna(freqGHz);}
		else{ConvCte_ = Brightness2Thermo(freqGHz);}
	}
	if(swap) ConvCte_ = 1.0/ConvCte_;
	return (ConvCte_ *= conv);
}


//----------------------------------

void	Cosmology::InitMassFunction(MassFunctionsType MassFunction,double MF_MinMass,double MF_MaxMass,double MaxZ,int MF_MaxSamples)
{
	delete MF_PrM_;
	delete MF_PrZM_;
	delete MF_Integrator_;

	MassFunction_ = MassFunction;

	MF_MaxSamples_		= MF_MaxSamples;
	MF_MinMass_			= std::log10(MF_MinMass);
	MF_MaxMass_			= std::log10(MF_MaxMass);
	MF_MaxZ_			= MaxZ;
	DeltaZ_				= MF_MaxZ_/MF_MaxSamples;
	DeltaM_				= (MF_MaxMass_ - MF_MinMass_)/MF_MaxSamples;

	MF_MakeSpace(MF_MaxSamples);
	FillMassValues();
	FillZConditValues();
	//DebugMassFunction1();
}

void	Cosmology::MF_MakeSpace(int MF_MaxSamples)
{
	MF_PrM_		= new MF_PrM_type(MF_MaxSamples);
	MF_PrZM_	= new MF_PrZM_type(MF_MaxSamples);
	MF_PrZM_type::iterator			piv(MF_PrZM_->begin());
	MF_PrZM_type::const_iterator	const end(MF_PrZM_->end());
	for(;piv != end;++piv)
	{
		(*piv).resize(MF_MaxSamples);
	}
	MF_Integrator_	= new SimpsonIntegrator<double,Cosmology>(MF_INITLEVEL,MF_MAXITER,this,&Cosmology::MF_dndMdzdSr_M0,MF_PREC,true);
}

void	Cosmology::FillMassValues(void)
{
	bool			dummy;
	double			CurrMass(MF_MinMass_ + (DeltaM_/2.0));
	MF_PrM_type::iterator				piv(MF_PrM_->begin());
	MF_PrM_type::const_iterator	const	end(MF_PrM_->end());	

	SimpsonIntegrator<double,Cosmology>
		Integrator(MF_INITLEVEL,MF_MAXITER,this,&Cosmology::MF_dndMdSr,MF_PREC,true);

	for(;piv!=end;++piv,CurrMass += DeltaM_)
	{
		*piv = Integrator.Integrate(MF_MinMass_,CurrMass,dummy);
	}
	std::sort(MF_PrM_->begin(),MF_PrM_->end());
	const double NormalCte(*(end - 1));
	piv = MF_PrM_->begin();
	for(;piv != end;++piv)
	{*piv /= NormalCte;}
}

void	Cosmology::FillZConditValues(void)
{
	double			CurrMass(MF_MinMass_ + (DeltaM_/2.0));
	double			CurrZ;
	bool			dummy;

	MF_PrZM_type::iterator				pivZM(MF_PrZM_->begin());
	MF_PrZM_type::const_iterator const	endZM(MF_PrZM_->end());	
	MF_PrZ_type::iterator				pivZ;
	MF_PrZ_type::const_iterator			endZ;	


	for(;pivZM != endZM;++pivZM,CurrMass += DeltaM_)
	{
		pivZ	= pivZM->begin();
		endZ	= pivZM->end();
		CurrZ	= DeltaZ_/2.0;
		SetPiv(std::exp(CurrMass*LN10));
		for(;pivZ != endZ;++pivZ,CurrZ += DeltaZ_)
		{
			*pivZ = MF_Integrator_->Integrate(0.0,CurrZ,dummy);
		}
		std::sort(pivZM->begin(),pivZM->end());
		pivZ	= pivZM->begin();
		const double NormalCte(*(endZ - 1));
		for(;pivZ != endZ;++pivZ)
		{*pivZ /= NormalCte;}
	}
}

void	Cosmology::DebugMassFunction0(void)
{
	double			CurrMass(MF_MinMass_ + (DeltaM_/2.0));
	double			CurrZ(DeltaZ_/2.0);
	Zeus::LArr2D<double>	tArray(MF_MaxSamples_*MF_MaxSamples_,MF_MaxSamples_);

	for(int i=0;i < MF_MaxSamples_;++i, CurrMass += DeltaM_)
	{
		CurrZ = (DeltaZ_/2.0);
		for(int j=0;j < MF_MaxSamples_;++j, CurrZ += DeltaZ_)
		{
			double sigma(-std::log(Sigma(std::exp(CurrMass*LN10),CurrZ)));
			tArray(i,j) = sigma;
		}
	}
	//DumpInOut_2d("sigmavalues",MF_MaxSamples_,MF_MaxSamples_,MF_MaxSamples_,0,tArray.begin(),1.0);
}

void	Cosmology::DebugMassFunction1(void)
{
	double			tDeltaM(2.0 * DeltaM_ );
	double			CurrMass;
	double			CurrZ((DeltaZ_/2.0) + (MF_MaxZ_/2.0));
	Zeus::LArr2D<double>	tArray(((MF_MaxSamples_)<<2)*MF_MaxSamples_,MF_MaxSamples_);
	int i;

	CurrZ = (DeltaZ_/2.0);
	for(i=0;i < ((MF_MaxSamples_)<<2);i += 4, CurrZ += DeltaZ_)
	{
		MassFunction_ = Jenkins;
		for(int k=0;k<2;++k)
		{
			CurrMass = (MF_MinMass_ + (tDeltaM/2.0));
			for(int j=0;j < MF_MaxSamples_;++j, CurrMass += tDeltaM)
			{
				double Sigma0(Sigma(std::exp(CurrMass*LN10),CurrZ));
				double sigma(-std::log(Sigma0));
				tArray(i+(k<<1),j) =   sigma;
				tArray(i+1+(k<<1),j) = GetMassFunction(Sigma0,CurrZ);
			}
			MassFunction_ = PressSchechter;
		}
	}
	//DumpInOut_2d("sigmavaluesPress",(MF_MaxSamples_)<<2,MF_MaxSamples_,MF_MaxSamples_,0,tArray.begin(),1.0);
}

void	Cosmology::InitCosmography(int MaxSamples)
{
	if((CosmographMaxZ_ / Step_) > (MaxSamples-1))
	{
		Step_ = CosmographMaxZ_ / (MaxSamples-1);
	}
	bool dummy(true);
	SimpsonIntegrator<double,Cosmology>
		Integrator(COMOVINGINITLEVEL,COMOVINGINTMAXITER,this,&Cosmology::ComovValues,COMOVINGPREC,true);
	Coefs_.push_back(0.0);
	const double MaxLimitZ(CosmographMaxZ_ + TOLCOSM);
	for(double currvalue = Step_;currvalue <= MaxLimitZ;currvalue += Step_)
	{
		Coefs_.push_back(Integrator.Integrate(0.0,currvalue,dummy));
	}
	LastIndex_ = static_cast<int>(Coefs_.size() - 1);
}

double	Cosmology::GetComovFactor(double RedShift) const
{

	if(RedShift <= 0.0)
	{
		return 0.0;
	}

	if(RedShift >= CosmographMaxZ_)
	{
		if(RedShift < (CosmographMaxZ_+ TOLCOSM))
			return Coefs_[LastIndex_];
		errRedshiftOutOfRange();
	}

	const double t(RedShift / Step_);
	const int index(toInt(t));

	return ((Coefs_[index + 1] - Coefs_[index]) * (t - static_cast<double>(index))) + Coefs_[index]; 
}


/*
double	Cosmology::PS_PressSchechter(double z,double M) const
{
	const double OmegaMz(PS_GetOmegaM(z));
	const double Sigma8z(Sigma8_ * (PS_GetOmegaM_Growth(0)/(PS_GetOmegaM_Growth(z)*(1.0 + z))));
	const double Ro0Bar(OmegaM_ * DensityCrit_0);
	const double ClusterRadius(std::pow(0.2387324 * (M / Ro0Bar), 0.33333333));
	const double BigGamma(OmegaM_ * H_SMALL * 0.98101540 * std::exp(-OmegaBar_ - (std::sqrt(2.0 * H_SMALL)*(OmegaBar_/OmegaM_))));
	const double LitlleGamma((0.3 * BigGamma + 0.2) * (2.92 + std::log10((ClusterRadius*H_SMALL)/ 8.0)));
	//const double RoZ_Bar(Ro0Bar);
	const double RoZ_Bar(Ro0Bar *(1.0 + z)*(1.0 + z)*(1.0 + z));
	const double SigmaRz(Sigma8z * std::pow(((ClusterRadius*H_SMALL)/8.0),-LitlleGamma));

	const double x((1.0 - OmegaMz)/0.7);
	const double A((1.0 - x)*0.27 + x*0.22);
	const double B((1.0 - x)*0.65 + x*0.73);
	const double eps((1.0 - x)*3.77 + x*3.86);
	const double f(A * std::exp(-std::pow(std::abs(std::log(1.0/SigmaRz) + B),eps)));

	return (LitlleGamma*RoZ_Bar*f) / (3.0 * M * M);
}
*/

double	Cosmology::GetMassFunction(double sigma,double z) const
{
	double A,B,eps,x;

	switch(MassFunction_)
	{
	case PressSchechter:
		x	=	(1.0 - GetOmegaMatRedshift(z)) / 0.7;
		A	=	(1.0 - x)*0.27 + x*0.22;
		B	=	((1.0 - x)*0.65 + x*0.73);
		eps	=	((1.0 - x)*3.77 + x*3.86);
		break;
	default:
		A	=	0.315;
		B	=	0.61;
		eps	=	3.8;
		break;
	}

	double SigDebug(- std::log(sigma));
	if((SigDebug < -1.2) || (SigDebug > 1.05))
	{
		//return 0.0;
		//printf("Ooops -> %7.6g\n",SigDebug);
	}

/*

		x = (1.-Omz)/0.7
      		A = (1.-x)*0.27 + x*0.22
	      	B = (1.-x)*0.65 + x*0.73
      		eps = (1.-x)*3.77 + x*3.86
		f = A*exp(-1.*abs(log(1./sigmaz)+B)**eps)



	double SigDebug(- std::log(sigma));
	if((SigDebug < -1.2) || (SigDebug > 1.05))
	{
		printf("Ooops -> %7.6g\n",SigDebug);
	}
*/
	return A * std::exp(-std::pow(std::abs(B - std::log(sigma)),eps));
}


double	Cosmology::MF_GetZSample(int M_Index,double deltaM,double ZSampleNormal,double MSampleNormal) const
{
	if(ZSampleNormal <= TOLCOSM)		return 0.0;
	if(ZSampleNormal >= (1.0-TOLCOSM))	return MF_MaxZ_;
	
	double	d00,d10,d01,d11;
	int		ZIndex0(-1);

	const MF_PrZ_type&				ZValues(*(MF_PrZM_->begin() + M_Index));
	MF_PrZ_type::const_iterator		const beg(ZValues.begin());
	const int ZIndex(static_cast<int>(std::upper_bound(beg,ZValues.end(),ZSampleNormal) - beg));
	d11 = ZValues[ZIndex];
	if(ZIndex)
	{d10 = ZValues[ZIndex-1];}
	else
	{d10 = 0.0;}
	if(M_Index>0)
	{
		const MF_PrZ_type&				ZValues0(*(MF_PrZM_->begin() + (M_Index-1)));
		MF_PrZ_type::const_iterator		const beg0(ZValues0.begin());
		ZIndex0	= (static_cast<int>(std::upper_bound(beg0,ZValues0.end(),ZSampleNormal) - beg0));
		d01 = ZValues0[ZIndex0];
		if(ZIndex0>0) d00 = ZValues0[ZIndex0-1];
		else d00 = 0.0;
	}

	double dy0(-1.0);
	double dy1((ZSampleNormal - d10)/(d11-d10));
	if(ZIndex0 >= 0)
	{dy0 = (ZSampleNormal - d00)/(d01-d00);}

	double	t0,t1;

	t1 = LinearInterpolation(ZIndex2Z(ZIndex-1),ZIndex2Z(ZIndex),dy1);

	if(ZIndex0 < 0) return t1;

	t0 = LinearInterpolation(ZIndex2Z(ZIndex0-1),ZIndex2Z(ZIndex0),dy0);
	
	return LinearInterpolation(t0,t1,deltaM);

}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  logfactorial
// 
// Purpose:   log( k! ) = lookup table up to 99, then logGamma(k+1)
//
// History:   JS         15 May 1998, 24 Mar 2001, 3 Oct 2002, 18 Aug 2003
//-----------------------------------------------------------------------------
// 
double logfactorial(int k)
{
static const double  logfact[100] = {  0.0000000000000000,
            0.0000000000000000,  0.6931471805599453,  1.7917594692280550,
            3.1780538303479456,  4.7874917427820460,  6.5792512120101010,
            8.5251613610654143, 10.6046029027452502, 12.8018274800814696,
           15.1044125730755153, 17.5023078458738858, 19.9872144956618861,
           22.5521638531234229, 25.1912211827386815, 27.8992713838408916,
           30.6718601060806728, 33.5050734501368889, 36.3954452080330536,
           39.3398841871994940, 42.3356164607534850, 45.3801388984769080,
           48.4711813518352239, 51.6066755677643736, 54.7847293981123192,
           58.0036052229805199, 61.2617017610020020, 64.5575386270063311,
           67.8897431371815350, 71.2570389671680090, 74.6582363488301644,
           78.0922235533153106, 81.5579594561150372, 85.0544670175815174,
           88.5808275421976788, 92.1361756036870925, 95.7196945421432025,
           99.3306124547874269,102.9681986145138126,106.6317602606434591,
          110.3206397147573954,114.0342117814617032,117.7718813997450715,
          121.5330815154386340,125.3172711493568951,129.1239336391272149,
          132.9525750356163099,136.8027226373263685,140.6739236482342594,
          144.5657439463448860,148.4777669517730321,152.4095925844973578,
          156.3608363030787852,160.3311282166309070,164.3201122631951814,
          168.3274454484276523,172.3527971391628016,176.3958484069973517,
          180.4562914175437711,184.5338288614494905,188.6281734236715912,
          192.7390472878449024,196.8661816728899940,201.0093163992815267,
          205.1681994826411985,209.3425867525368356,213.5322414945632612,
          217.7369341139542273,221.9564418191303340,226.1905483237275933,
          230.4390435657769523,234.7017234428182677,238.9783895618343230,
          243.2688490029827142,247.5729140961868839,251.8904022097231944,
          256.2211355500095255,260.5649409718632093,264.9216497985528010,
          269.2910976510198225,273.6731242856937041,278.0675734403661429,
          282.4742926876303960,286.8931332954269940,291.3239500942703076,
          295.7666013507606240,300.2209486470141318,304.6868567656687155,
          309.1641935801469219,313.6528299498790618,318.1526396202093268,
          322.6634991267261769,327.1852877037752172,331.7178871969284731,
          336.2611819791984770,340.8150588707990179,345.3794070622668541,
          349.9541180407702369,354.5390855194408088,359.1342053695753988};
    double s, y;

	if(k < 0) return 0.0;
    if( k < 100 )
        return logfact[k];
    y = k + 1.0;
    s = y + 2.269488974204959960;
    s = y + 1.517473649153287398 / s;
    s = y + 1.011523068126841711 / s;
    s = y + 0.525606469002695417 / s;
    s = y + 0.252380952380952380 / s;
    s = y + 0.033333333333333333 / s;
    s =     0.083333333333333333 / s;
    s = s + 0.91893853320467 - y + (y - 0.5) * log(y);
    return s;
}

void	NormalizeProb(double * const buffer,unsigned long sz,unsigned long metric)
{
	double MaxLog(Zeus::logZERO);
	double MinLog(Zeus::logINF);
	wchar_t buf[BUFFERMAXCHAR];
	double logAccuml(0.0);	

	for(unsigned long i=0;i<sz;++i)
	{
		if(buffer[i] > MaxLog) MaxLog = buffer[i];
		if(buffer[i] < MinLog) MinLog = buffer[i];
		logAccuml += buffer[i]; 
	}

	double temp;
	double Accumul(0.0);

	for(unsigned long i=0;i<sz;++i)
	{
		buffer[i] = std::exp(buffer[i] - MaxLog);
		Accumul += buffer[i];
	}

	for(unsigned long i=0;i<sz;++i)
	{
		buffer[i] /= Accumul;
	}
	
//	PRINTINTOBUFFERFUNCT(buf,BUFFERMAXCHAR,L"\nContours, logMin -> %9.3f, logMax -> %9.3f, logAcum -> %9.3f, Calib -> %9.3f\n",MinLog,MaxLog,logAccuml,Accumul);

//	(Zeus::ConManager::Instance())->PrintStr2Console(buf);


/*
// Swaps the axis *********** BIG MISTAKE *********************
	for(unsigned long l=0;l<Nlines;++l)
	{
		*(buffer + (metric * l) + l) /= Accumul;

		for(unsigned long c=(l+1);c<metric;++c)
		{
			temp = *(buffer + (metric * c) + l) / Accumul;
			*(buffer + (metric * c) + l) = *(buffer + (metric * l) + c) / Accumul;
			*(buffer + (metric * l) + c) = temp;
		}
	}
*/
}


}// end of namespace Zeus
//-----------

