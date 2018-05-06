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



//----------------------------------------------------------------
#ifndef PSNAKES_PHYSICSMATHH
#define PSNAKES_PHYSICSMATHH

#include <limits>
#include <cmath>
#include <algorithm>
#include <functional>
#include "PWS_Strings.h"
#include "ZEUS_Exceptions.h"
#include "ZEUS_Debug.h"

//----------------------------------------------------------------

// 5.0 * sqrt(2)
#define	GAUSSIANLIMIT					7.0710678118654752440084436210485
#define GAMMA_FPMIN						((double)1.0e-30)
#define SINC_TINY						((T)1.0e-30)
#define GAMMA_EPS						((double)3.0e-7)
#define GAMMA_ITMAX						100
#define INTG_MIN_PREC					((T)1e-7)
#define ROUND_EPS						((double) 1.0e-30)
#define TOLCOSM							1e-6
#define PI								3.14159265359
#define PI2								9.86960440109
#define PITIMES2SQ						19.7392088022
#define PIOVER2							1.57079632680
#define PITIMES4						12.5663706144
#define PITIMES2						6.28318530718
#define SQRTPI							1.77245385091
#define SQRT2PI							2.506628274631
#define SQRT2							1.41421356237
#define LN10							2.30258509299
#define dBCONV							4.34294481904
#define LNPI							1.1447298858494
#define LN2								0.693147180559945
#define RAD2DEGREE						57.2957795131
#define RAD2ARCMIN						3437.74677079
#define SR2ARCMIN2						11818102.86
#define BRIGHT2ANTENNA					32.54800036199727
#define BRIGHT2THERMO					1.007342360844192e-2
#define CMBTEMP							2.7255
#define LIGHTSPEED						299792458.0
#define BOLTZMANNCTE					1.380650524e-23					
#define PLANCKCTE						6.626069311e-34

// Tcmb = 2.728 
// Freq / 100.0 Ghz units
#define ADIMFREQCTE						1.759243278 // ~ 100.0/56.8
#define	FWHM2SIGMA						2.35482004503
#define ERRPOS95CTE						2.485975524

#define	CMBTHERMO256					6.2271432814507664e-005 //thermo fluctuations through an antenna 41.2 arcmin
#define	CMBTHERMO512					8.0363808564668673e-005	//thermo fluctuations through an antenna 20.6 arcmin

#define COMOVINGSAMPLPER				0.1						
#define COMOVINGMAXSAMPLES				200
#define COMOVINGINITLEVEL				5
#define COMOVINGINTMAXITER				100
#define COMOVINGPREC					5.0e-5
#define MPARSEC2M						3.08568e22
#define MPARSEC2CM						3.08568e24
#define MF_INITLEVEL					5
#define MF_MAXITER						100
#define MF_PREC							1.0e-2
#define MF_MAXSAMPLES					100

//#define PSM_ABUNDANCE_CORRECTION		6.80
// This rather artificial cte corrects the value of the  formula
// to the values of PSM; WHY ? I simply don't know !
#define PSM_ABUNDANCE_CORRECTION		1.00

#define H_SMALL							0.719

//Farhan
//#define H_SMALL							0.7

// Km s^-1 Mpc^-1
#define HUBLECONST						(100.0 * H_SMALL)
//Mpc^3 M_sun^-1 s^-2
#define G_NEWTON						4.518e-48
// Mass = Solar masses; l = Mpc; T = KeV

#define SZNGAIFLUX_Y_CTE				8.1224e-19
// Pressure ergs cm^-3
// SigmaT = 6.65 10^-25  cm^2
// me0C^2 = 8.1872 10^-7 ergs

#define SZBETAFLUX_Y_CTE				1.424e-19
// Mega parsec
const	double	HubleDistance = (3000.0 / H_SMALL);
// M_sun Mpc^3
// Mpc = 3.0857e19 Km
const	double	DensityCrit_0 = ((3.0 * (HUBLECONST/3.0857e19)*(HUBLECONST/3.0857e19))/(8.0*PI*G_NEWTON));

#define OMEGAMDEF						0.256
#define OMEGALDEF						0.744
#define OMEGAKDEF						0.0
#define OMEGABARDEF						(0.0196 / (H_SMALL * H_SMALL))
#define SIGMA8							0.796

#include <vector>
#include "ZEUS_General.h"

namespace Zeus
{

const double logINF  = std::numeric_limits<double>::max()  * std::numeric_limits<double>::epsilon();    // ln(inf) = infinity
const double INF = logINF;
const double logZERO = -std::numeric_limits<double>::max() * std::numeric_limits<double>::epsilon();   // ln(0) = -infinity
const double pwsNaN = std::numeric_limits<double>::quiet_NaN();
const double pwsEps = std::numeric_limits<double>::epsilon();
//
inline int pwsIsNaN(double x)
{
	return (x != x);
}
//
inline int pwsIsNaN(float x)
{
	return (x != x);
}
//
template<typename T>
inline T AddLog(T x,T y)
{return ( x > y ? x+std::log(1.0+std::exp(y - x)): y+std::log(1.0+std::exp(x - y)));}
//
//
/*
template<typename T>
inline T AddLog(T x,T y)
{return ( x > y ? x+std::log(1.0+(((y - x)<-100.0)?0.0:std::exp(y - x))): y+std::log(1.0+(((x - y)<-100.0)?0.0:std::exp(x - y))));}
*/
//
/*
template<typename T>
inline T Log1plusLog(T y)
{return ( y < 0 ? std::log(1.0+ ((y<-100.0)?0.0:std::exp(y))): y +std::log(1.0+ ((y>100.0)?0.0:std::exp(-y))));}
*/
//
template<typename T>
inline T Log1plusLog(T y)
{return ( y < 0 ? std::log(1.0+std::exp(y)): y +std::log(1.0+std::exp(-y)));}
//
/*
template<typename T>
inline T SubLog(T x,T y)
{
	return ((y >= x)?logZERO:x+std::log(1.0- ((y - x)<-100.0?0.0:std::exp(y - x))));
}
*/
//
template<typename T>
inline T SubLog(T x,T y)
{
	return ((y >= x)?logZERO:x+std::log(1.0-std::exp(y - x)));
}
//
/*
template<typename T>
inline T SubLogFrom1(T y)
{
	return std::log(1.0-std::exp(y));
}
*/
//
double logfactorial(int k);
//
void	NormalizeProb(double * const buffer,unsigned long sz,unsigned long metric);
//
inline int	toInt(double d)
{
	double t(std::floor(d));
	return static_cast<int>(t<0?t-ROUND_EPS:t+ROUND_EPS);
}
//
inline int getIntSQRT(int v)
{
	int i(1);
	for(;(i*i)<=v;++i)
		;
	return --i;
}
//
inline int getLessPow2p1(int v)
{
	if(v<=1) return 1;
	int p(0);
	for(;v>0;++p)
	{v>>=1;}
	return (1<<(p-1))+1;
}
//
template<typename T>
inline T	Sq(const T& fst) {return fst*fst;}
//
inline std::complex<double>	Sq(const std::complex<double>& fst)
{return std::complex<double>(fst.real() * fst.real() + fst.imag() * fst.imag(),0.0);}
//
inline std::complex<float>	Sq(const std::complex<float>& fst)
{return std::complex<float>(fst.real() * fst.real() + fst.imag() * fst.imag(),0.0);}
//
template<typename T>
inline T	Mult(const T& fst,const T& sec) {return fst*fst;}
//
inline std::complex<double>	Mult(const std::complex<double>& fst,const std::complex<double>& sec)
{return std::complex<double>(fst.real() * sec.real() + fst.imag() * sec.imag(),fst.real() * sec.imag() - fst.imag() * sec.real());}
//
inline std::complex<float>	Mult(const std::complex<float>& fst,const std::complex<float>& sec)
{return std::complex<float>(fst.real() * sec.real() + fst.imag() * sec.imag(),fst.real() * sec.imag() - fst.imag() * sec.real());}
//
template<typename T>
inline std::complex<T> Conjugate(const std::complex<T>& v)
{return std::complex<T>(v.real(),-v.imag());}
//
template<typename T>
inline T CxPower(const std::complex<T>& v)
{return static_cast<T>(v.real()*v.real()+v.imag()*v.imag());}
//
template<typename T>
inline T Cx_XPower(const std::complex<T>& f,const std::complex<T>& s)
{return static_cast<T>(f.real()*s.real()+f.imag()*s.imag());}
//
template<typename T>
struct PowerFunctor
{
	inline T	operator()(const std::complex<T>& x)
	{
		return static_cast<T>((x.real() * x.real()) + (x.imag() * x.imag()));
	}
};
//
template<typename T>
struct MultFunctor
{
	const T CteValue_;
	explicit MultFunctor(T CteValue)
		:CteValue_(CteValue)
	{}
	inline void	operator()(T& x, int dummy) const
	{
		x *= CteValue_;
	}
};
//
template<typename T>
struct SetValueFunctor
{
	const T CteValue_;
	SetValueFunctor(T CteValue)
		:CteValue_(CteValue)
	{}
	inline void	operator()(T& x, int dummy) const
	{x = CteValue_;}
};
//
template<typename T>
inline T BiLinearInterpolation(T BaseFunct,T YplusFunct,T XplusFunct,T XYplusFunct,double dY,double dX)
{
	double unit_dX(1 - dX);
	double unit_dY(1 - dY);

	return (BaseFunct * (unit_dX * unit_dY)) + (XplusFunct * (dX * unit_dY)) + (YplusFunct * (dY * unit_dX)) + (XYplusFunct * (dX * dY));
}

template<typename T>
inline T LinearInterpolation(T BaseFunct,T plusFunct,T dX)
{
	return (plusFunct - BaseFunct) * dX + BaseFunct;
}
//
double Dierfc(double y);
//
class	GammaFuncts
{
	double					eps_;
	double					betaa_,betab_;
	int						itMax_;
	bool					UseExcpt_;
	static	const double	coef_[];
public:
	 GammaFuncts(bool UseExcpt = true,double eps=GAMMA_EPS, int itMax=GAMMA_ITMAX)
		:eps_(eps),itMax_(itMax),UseExcpt_(UseExcpt),betaa_(pwsNaN),betab_(pwsNaN)
	{}
	inline void		SetParams(double eps, int itMax,bool UseExcpt)
	{eps_=eps;itMax_=itMax;UseExcpt_=UseExcpt;}
	inline void		GetParams(double& eps, int& itMax,bool& UseExcpt) const
	{eps=eps_;itMax=itMax_;UseExcpt=UseExcpt_;}

	inline double	erf(double x) const
	{
		if(std::abs(x)>=((double)GAUSSIANLIMIT))
		{return x<((double)0.0)?-((double)1.0):((double)1.0);}

		return (x<((double)0.0)?-gammp(((double)0.5),x*x):gammp(((double)0.5),x*x));
	}

	inline double	gammp(double a,double x) const
	{
		if((x<((double)0.0)) || (a<=((double)0.0)))
			throw libException(ERROR_COD_ZEUSOUTBOUNDS,ERROR_MSG_ZEUSOUTBOUNDS,L"Gamma::gammp");
		if(x<(a+((double)1.0))) return gser(a,x);
		else return ((double)1.0) - gcf(a,x);
	}
//
	inline int		setbetaAB(double a,double b)
	{
		if((a <= 0.0) || (b <= 0.0))
			return false;
		betaa_ = a;
		betab_ = b;
		return true;
	}
//
	inline double	beta(double a,double b) const
	{
		return std::exp(gammln(a)+gammln(b)-gammln(a+b));
	}
//
	inline double	beta(void) const
	{
		if(pwsIsNaN(betaa_) || pwsIsNaN(betab_))
			throw libException(ERROR_COD_ZEUSOUTBOUNDS,ERROR_MSG_ZEUSOUTBOUNDS,L"Gamma::beta");
		return std::exp(gammln(betaa_)+gammln(betab_)-gammln(betaa_+betab_));
	}
//
	inline double	betainc(double a,double b,double x) const
	{
		double bt;

		if((a <= 0.0) || (b <= 0.0))
			goto betainc_excpt0;

		if(x < 0.0)
		{
			if(UseExcpt_)
				goto betainc_excpt0;
			else
				x = 0.0;
		}
		if(x > 1.0)
		{
			if(UseExcpt_)
				goto betainc_excpt0;
			else
				x = 1.0;
		}
		if(x == 0.0) return 0;
		if(x == 1.0) return 1.0;

		bt = std::exp(gammln(a+b)-gammln(a)-gammln(b)+a*std::log(x)+b*log(1.0-x));
		
		if(x < (a + 1.0)/(a+b+2.0))
			return bt*betacf(a,b,x)/a;
		else
			return 1.0-(bt*betacf(b,a,1.0-x)/b);

betainc_excpt0:
		throw libException(ERROR_COD_ZEUSOUTBOUNDS,ERROR_MSG_ZEUSOUTBOUNDS,L"Gamma::betainc");
	}
//
	inline double	betainc(double x)
	{
		double bt;

		if(pwsIsNaN(betaa_) || pwsIsNaN(betab_))
			goto betainc_excpt1;

		if(x < 0.0)
		{
			if(UseExcpt_)
				goto betainc_excpt1;
			else
				x = 0.0;
		}
		if(x > 1.0)
		{
			if(UseExcpt_)
				goto betainc_excpt1;
			else
				x = 1.0;
		}
		if(x == 0.0) return 0;
		if(x == 1.0) return 1.0;

		bt = std::exp(gammln(betaa_+betab_)-gammln(betaa_)-gammln(betab_)+betaa_*std::log(x)+betab_*log(1.0-x));
		
		if(x < (betaa_ + 1.0)/(betaa_+betab_+2.0))
			return bt*betacf(betaa_,betab_,x)/betaa_;
		else
			return 1.0-(bt*betacf(betab_,betaa_,1.0-x)/betab_);

betainc_excpt1:
		throw libException(ERROR_COD_ZEUSOUTBOUNDS,ERROR_MSG_ZEUSOUTBOUNDS,L"Gamma::betainc2");
	}
//
	double			gammln(double xx) const;
//
private:
	double			gser(double a,double x) const;
//
	double			gcf(double a,double x) const;
//
	double			betacf(double a,double b,double x) const;	
};


template<typename T,typename FUNCTOR>
class	SimpsonIntegrator
{
	bool			UseExcept_;
	int				MaxIter_;
	int				StartLevel_;
	FUNCTOR			*Obj_;
	T				(FUNCTOR::*funct_)(T arg);
	T				errPrec_;
	T				CurrAprox_;

	T	trapzd(T lLimit,T HLimit, int steps)
	{
		int		it,j,tnm;
		T		del,x,sum;

		if(steps == 1) return CurrAprox_ = ((T)0.5) * (HLimit - lLimit) * ((Obj_->*funct_)(lLimit) + (Obj_->*funct_)(HLimit));
		for(it=1,j=1;j<(steps-1);++j) it <<= 1;
		tnm = it ;
		del = (HLimit-lLimit)/tnm;
		x=lLimit+((T)0.5)*del;
		for(sum=((T)0.0),j=1;j<=it;++j,x+=del) sum += (Obj_->*funct_)(x);
		T temp(((HLimit - lLimit)*sum)/tnm); 
		if(steps == StartLevel_) return CurrAprox_ = temp;
		return CurrAprox_ = ((T)0.5)*(temp + CurrAprox_);
	}
public:
	inline void SetParams(int StartingLevel,int MaxIter,T errPrec,FUNCTOR *Obj,T (FUNCTOR::*funct)(T arg))
	{StartLevel_ = StartingLevel;MaxIter_ = MaxIter;Obj_=Obj;funct_ = funct;errPrec_=errPrec;}

	SimpsonIntegrator(bool UseExcept=false):StartLevel_(0),MaxIter_(0),Obj_(0),funct_(0),errPrec_(0),UseExcept_(UseExcept)
	{}
	SimpsonIntegrator(int StartingLevel,int MaxIter,FUNCTOR *Obj,T (FUNCTOR::*funct)(T arg) ,T errPrec,bool UseExcept=false)
		:StartLevel_(StartingLevel),MaxIter_(MaxIter),Obj_(Obj),funct_(funct),errPrec_(errPrec),UseExcept_(UseExcept)
	{}
	inline T GetLastIter(void) const
	{return CurrAprox_;}
	T	Integrate(T LowLimit,T HighLimit,bool& IntOk)
	{
		CurrAprox_ = 0;

		T s,st,ost,os,t_prec;

		ost = os = ((T)-1.0e30);
		for(int j=StartLevel_;j<=MaxIter_;++j)
		{
			st	= trapzd(LowLimit,HighLimit,j);
			s	= (((T)4.0)*st-ost)/((T)3.0);
			t_prec = errPrec_*std::abs(os);
			if(t_prec<INTG_MIN_PREC) t_prec = INTG_MIN_PREC;
			if(std::abs(s-os) < t_prec)
			{IntOk = true;CurrAprox_ = s;return s;}
			os = s;ost = st;
		}
		if(UseExcept_)
		{CurrAprox_ = s;throw libException(ERROR_COD_ZEUSFDIVERGES,ERROR_MSG_ZEUSFDIVERGES L"SimpsonIntegrator",this);}
		IntOk = false;
		return s;
	}
};


template<typename T>
inline T	sinc(T arg)
{
	if(std::abs(arg) < SINC_TINY) return static_cast<T>(1.0);
	return	static_cast<T>(std::sin(PI * arg) / (PI * arg));
}

//
class	UnitsValueTranform
{
	double ConvCte_;
/*
       2  x      
      x ·e       
-----------------
  2·x      x     
 e   - 2·e  + 1

*/
	inline static double Thermo2Antenna(double freqGHz)
	{
		double AdimFreq((freqGHz*ADIMFREQCTE)/100.0);
		double tempexp(std::exp(AdimFreq));
		return (AdimFreq*AdimFreq*tempexp) / ((tempexp-1)*(tempexp-1));
	}

	inline static double Brightness2Thermo(double freqGHz)
	{
		double AdimFreq((freqGHz*ADIMFREQCTE)/100.0);
		double tempexp(std::exp(AdimFreq));
		return BRIGHT2THERMO * ((tempexp-1)*(tempexp-1)) /(tempexp*AdimFreq*AdimFreq*AdimFreq*AdimFreq);
	}
public:
	//2 - Brightness MJy/sr 
	//0 -  AntennaT in K
	//1 - Thermo in K

	enum UnitsType{ANTENNA_T=0,THERMO_T=1,BRIGHTNESS=2};

	UnitsValueTranform(void): ConvCte_(1.0)
	{}
	UnitsValueTranform(UnitsType UnitIn,UnitsType UnitOut,double freqGHz,double conv)
	{Set(UnitIn,UnitOut,freqGHz,conv);}
	double	Set(UnitsType UnitIn,UnitsType UnitOut,double freqGHz,double conv);
	inline double operator() (double value) const
	{return value * ConvCte_;}
};
//
//
class	PlanckUnitsValuesTranform
{
	double ConvCte_;

	static const double	PlanckChFreq_[];
	static const double	KCMB2KRJ[];
	static const double	MJYSR2KCMB[];
	static const double	MJYSR2KRJ[];

	static const int	PlanckChannelsN;

	inline static int GetPlanckChannel(double freqGHz)
	{
#ifdef PLANCKNOMINALCHANNELS
		return -1;
#else
		for(int i=0;i<PlanckChannelsN;++i)
		{
			if(std::abs((PlanckChFreq_[i]-freqGHz)/freqGHz) < 0.05)
				return i;
		}
		return -1;
#endif
	}

	inline static double Thermo2AntennaTheo(double freqGHz)
	{
		double AdimFreq((freqGHz*ADIMFREQCTE)/100.0);
		double tempexp(std::exp(AdimFreq));
		return (AdimFreq*AdimFreq*tempexp) / ((tempexp-1)*(tempexp-1));
	}

	inline static double Brightness2ThermoTheo(double freqGHz)
	{
		double AdimFreq((freqGHz*ADIMFREQCTE)/100.0);
		double tempexp(std::exp(AdimFreq));
		return BRIGHT2THERMO * ((tempexp-1)*(tempexp-1)) /(tempexp*AdimFreq*AdimFreq*AdimFreq*AdimFreq);
	}

	inline static double Brightness2AntennaTheo(double freqGHz)
	{
		return BRIGHT2ANTENNA / (freqGHz*freqGHz);
	}

	inline static double Brightness2Antenna(double freqGHz)
	{
		const int PlnkId(GetPlanckChannel(freqGHz));
		return (PlnkId<0)?Brightness2AntennaTheo(freqGHz):MJYSR2KRJ[PlnkId];
	}

	inline static double Thermo2Antenna(double freqGHz)
	{
		const int PlnkId(GetPlanckChannel(freqGHz));
		return (PlnkId<0)?Thermo2AntennaTheo(freqGHz):KCMB2KRJ[PlnkId];
	}

	inline static double Brightness2Thermo(double freqGHz)
	{
		const int PlnkId(GetPlanckChannel(freqGHz));
		return (PlnkId<0)?Brightness2ThermoTheo(freqGHz):MJYSR2KCMB[PlnkId];
	}

public:
	//2 - Brightness MJy/sr 
	//0 -  AntennaT in K
	//1 - Thermo in K

	enum UnitsType{ANTENNA_T=0,THERMO_T=1,BRIGHTNESS=2};

	PlanckUnitsValuesTranform(void): ConvCte_(1.0)
	{}
	PlanckUnitsValuesTranform(UnitsType UnitIn,UnitsType UnitOut,double freqGHz,double conv)
	{Set(UnitIn,UnitOut,freqGHz,conv);}
	double	Set(UnitsType UnitIn,UnitsType UnitOut,double freqGHz,double conv);
	inline double operator() (double value) const
	{return value * ConvCte_;}
};
//
template<typename T>
inline T Dist2ptgs(T Colat1,T Long1,T Colat2,T Long2)
{
	const	T Lat1(PIOVER2 - Colat1);
	const	T Lat2(PIOVER2 - Colat2);
	const	T sin1(std::sin((Lat2-Lat1)/2.0));
	const	T sin2(std::sin((Long2-Long1)/2.0));

	return 2.0 * std::asin(std::sqrt(sin1 * sin1 + std::cos(Lat1) * std::cos(Lat2) * sin2 * sin2));
}
//
template<typename T>
inline T		FindMaxQuadFromSamples(T fa,T fb,T fc)
{return static_cast<T>(- 0.5 * ((fa - fc)/(((T)2.0)*fb -(fc + fa))));}
//
//inline	double HealPixGetPixSZ(int NSide)
//{return PITIMES4 / (12.0 * static_cast<double>(NSide) * static_cast<double>(NSide));}
//
class Cosmology
{

public:
	enum	MassFunctionsType {PressSchechter=1,Jenkins=2};

	Cosmology(double MaxZ,double Sigma8=SIGMA8,double OmegaM=OMEGAMDEF,double OmegaL=OMEGALDEF,double OmegaK=OMEGAKDEF,
		double OmegaBar=OMEGABARDEF)
		:OmegaM_(OmegaM),OmegaL_(OmegaL),OmegaK_(OmegaK),OmegaBar_(OmegaBar),Sigma8_(Sigma8),MassFunction_(Jenkins),
		CosmographMaxZ_(MaxZ),Step_(COMOVINGSAMPLPER),
		MF_PrM_(0),MF_PrZM_(0),MF_Integrator_(0)
	{}

	void			InitCosmography(int MaxSamples=COMOVINGMAXSAMPLES);
	void			InitMassFunction(MassFunctionsType MassFunction,double MF_MinMass,double MF_MaxMass,double MaxZ,int MF_MaxSamples=MF_MAXSAMPLES);
	double			GetComovFactor(double RedShift) const;
	inline double	GetOmegaMatRedshift(double RedShift) const
	{
		const double t((1.0 + RedShift)*(1.0 + RedShift)*(1.0 + RedShift));
		return (OmegaM_ * t)/(1.0 - OmegaM_ + (t * OmegaM_));
	}

	inline double	growth_ofOmegaMZ(double OmegaMZ) const
	{
		return (2.5 * OmegaMZ) / (0.0142857143 + 1.492857143 * OmegaMZ - ((OmegaMZ*OmegaMZ) / 140.0) + std::pow(OmegaMZ,0.57142857143));			
	}

	inline double	Sigma8atZ(double z) const
	{
		return (Sigma8_ / (1.0 + z)) * (growth_ofOmegaMZ(GetOmegaMatRedshift(z))/growth_ofOmegaMZ(OmegaM_));
	}


	inline double	Sigma(double Mass,double z) const
	{
		const double radius(GetRadiusFromMass(Mass));
		return Sigma8atZ(z) * std::pow((radius*H_SMALL)/8.0,-littleGamma(radius));
	}

	inline double	GetHubleFactor(double RedShift) const
	{
		const double z(RedShift + 1.0);
		return std::sqrt((OmegaM_  * z * z * z + OmegaK_ * z * z + OmegaL_));
	}

	inline double	GetHubleCte(double RedShift) const
	{ return HUBLECONST * GetHubleFactor(RedShift);}

	inline double	GetLineOfSightComDist(double RedShift) const
	{
		return GetComovFactor(RedShift) * HubleDistance;
	}

	inline double	GetTransverseComovDist(double RedShift) const
	{
		const double ComovFactor(GetComovFactor(RedShift));

		if(OmegaK_ == 0.0)
		{
			return ComovFactor * HubleDistance;
		}
		else if(OmegaK_ > 0.0)
		{
			const double t(std::sqrt(OmegaK_));
			return HubleDistance * (1.0 / t)*std::sinh(t * ComovFactor);
		}
		else
		{
			const double t(std::sqrt(std::abs(OmegaK_)));
			return HubleDistance * (1.0 / t)*std::sin(t * ComovFactor);		
		}
	}

	inline double	GetAngularDiameter(double RedShift) const
	{
		return GetTransverseComovDist(RedShift) / (1 + RedShift);
	}

	inline double	GetComovVolume(double RedShift) const
	{
		const double TComovDist(GetTransverseComovDist(RedShift));
		return (HubleDistance * TComovDist * TComovDist) / GetHubleFactor(RedShift);
	}

	void			DebugMassFunction0(void);
	void			DebugMassFunction1(void);

	~Cosmology(void)
	{
		delete MF_PrM_;
		delete MF_PrZM_;
		delete MF_Integrator_;
	}
	inline	double	GetRoZ_Bar(double z) const
	{
		return (OmegaM_ * DensityCrit_0) * (1.0 + z)*(1.0 + z)*(1.0 + z);
	}

	// Number of clusters M,M+dM and V,V+dV at z
	inline 	double	MF_dndMdV_atZ(double z,double M) const
	{
		const double RoZ_Bar(GetRoZ_Bar(z));
		const double RadiusFromMass(GetRadiusFromMass(M));
		const double tlittleGamma(littleGamma(RadiusFromMass));
		const double tSigma(Sigma(M,z));
		const double MassFunction(GetMassFunction(tSigma,z));
		const double TrueValue((tlittleGamma*RoZ_Bar*MassFunction) / (3.0 * M * M));

		//printf("\nZ             => %g\n",z);
		//printf("M             => %g\n",M);
		//printf("Radius        => %g\n",RadiusFromMass);
		//printf("RoZ_Bar       => %g\n",RoZ_Bar);
		//printf("littleGamma   => %g\n",tlittleGamma);
		//printf("Sigma         => %g\n",tSigma);
		//printf("MassFunction  => %g\n",MassFunction);
		//printf("N of clust M,M+dM and V,V+dV at z %g\n",TrueValue);
		return TrueValue / PSM_ABUNDANCE_CORRECTION;
	}

	// Number of clusters M,M+dM   z,z+dz, Sr,Sr+dSr
	inline double	MF_dndMdzdSr(double z,double M)
	{
		const double t(MF_dndMdV_atZ(z,M));
		const double ComovVolume(GetComovVolume(z));
		// dV / (dz dSr) = GetComovVolume(z)
		//printf("Comoving Volume => %g\n",ComovVolume);

		return t*ComovVolume;
	}

	inline void		GetMassFunctionSample(double MNormal,double ZNormal,double& Mass,double& Z) const
	{
		if(!MF_PrM_) errNotInitialised();
		int		index;
		double	deltaM;

		Mass	= MF_GetMSample(index,deltaM,MNormal);
		Z		= MF_GetZSample(index,deltaM,ZNormal,MNormal);
	}

	inline double	GetMassFunctionExpectN(double ZMin,double ZMax,double MMin,double MMax)
	{
		bool dummy;
		if(!MF_PrM_) errNotInitialised();
		MF_HelpZmin_ = ZMin;
		MF_HelpZmax_ = ZMax;
		SimpsonIntegrator<double,Cosmology>
			Integrator(MF_INITLEVEL,MF_MAXITER,this,&Cosmology::MF_dndMdSr_zRange,MF_PREC,true);
		return Integrator.Integrate(std::log10(MMin),std::log10(MMax),dummy);
	}
	
	inline void		ReleaseMassFunction(void)
	{
		delete MF_PrM_;MF_PrM_ = 0; 
		delete MF_PrZM_;MF_PrZM_ = 0;
		delete MF_Integrator_;MF_Integrator_ = 0;
	}

private:
	Cosmology(const Cosmology&);
	Cosmology& operator=(const Cosmology&);

	typedef		std::vector<double>					MF_PrM_type;
	typedef		std::vector<double>					MF_PrZ_type;
	typedef		std::vector<MF_PrZ_type>			MF_PrZM_type;

	inline double	ComovValues(double RedShift)
	{
		return 1.0/GetHubleFactor(RedShift);
	}

	inline  void	errParameterOutOfRange(const std::wstring& ParamName) const
	{
		std::wstring errstring(ERRMSG_PWS_PARAMOUTOFRANGE);
		errstring += std::wstring(L" -> ");
		errstring += ParamName;
		throw Zeus::libException(ERRCOD_PWS_PARAMOUTOFRANGE,errstring,*this);
	}

	inline  void	errRedshiftOutOfRange(void) const
	{errParameterOutOfRange(PHYMATH_REDSHIFT);}
	inline  void	errNotInitialised(void) const
	{
		throw Zeus::libException(ERRCOD_PWS_NOTINITIALISED,ERRMSG_PWS_NOTINITIALISED,*this);
	}

	inline  void	errCurvMustBe0(void) const
	{errParameterOutOfRange(PHYMATH_0CURVONLY);}

	inline	double	BigGamma(void) const
	{
		return OmegaM_ * H_SMALL * (7.29/(CMBTEMP * CMBTEMP)) * std::exp(-OmegaBar_ - (std::sqrt(H_SMALL/0.5) * (OmegaBar_/OmegaM_)));
	}
	inline	double	littleGamma(double radius) const
	{
		return (0.3 * BigGamma() + 0.2) * (2.92 + std::log10((radius*H_SMALL)/8.0));
	}
	inline	double	GetRadiusFromMass(double Mass) const
	{
		return std::pow(0.23873241 * (Mass / (OmegaM_ * DensityCrit_0)),0.333333333);
	}

	double			GetMassFunction(double sigma,double z) const;

	inline	double	MF_dndMdSr(double M)
	{
		bool	dummy;
		SetPiv(std::exp(M*LN10));
		return MF_Integrator_->Integrate(0.0,MF_MaxZ_,dummy);
	}

	inline	double	MF_dndMdSr_zRange(double M)
	{
		bool	dummy;
		SetPiv(std::exp(M*LN10));
		//printf("\nMass             => %g\n",std::exp(M*LN10));
		return MF_Integrator_->Integrate(MF_HelpZmin_,MF_HelpZmax_,dummy);
	}

	inline void		SetPiv(double newPiv)
	{
		piv_ = newPiv;
	}


	inline double	MF_dndMdzdSr_M0(double z)
	{
		return piv_ * LN10 * MF_dndMdzdSr(z,piv_);
	}

/*
	inline double	MF_dndMdzdSr_M0(double z)
	{
		return piv_ * LN10 * MF_Dummy(z,piv_);
	}

	inline double	MF_Dummy(double z,double M)
	{
		return std::sqrt(z)*std::log(M);
	}
*/
	inline	double	MIndex2Mass(int Index) const
	{
		if(Index < 0) return MF_MinMass_; 
		return MF_MinMass_ + DeltaM_ * (0.5 + static_cast<double>(Index));
	}

	inline	double	ZIndex2Z(int Index) const
	{
		if(Index < 0) return 0.0; 
		return DeltaZ_ * (0.5 + static_cast<double>(Index));	
	}

	inline	double	MF_GetMSample(int& MIndex,double& deltaM,double MSampleNormal) const
	{
		const MF_PrM_type& MValues(*MF_PrM_);

		if(MSampleNormal >= 1.0) MSampleNormal = (1.0-TOLCOSM);

		MF_PrM_type::const_iterator			const beg(MValues.begin());
		MIndex = static_cast<int>(std::upper_bound(beg,MValues.end(),MSampleNormal) - beg);
		const double v0(MIndex>0?MValues[MIndex-1]:0.0);
		const double d(MValues[MIndex] - v0);
		deltaM = (MSampleNormal - v0)/d;
		return std::exp(LN10 * LinearInterpolation(MIndex2Mass(MIndex-1),MIndex2Mass(MIndex),deltaM));
	}

	double			MF_GetZSample(int M_Index,double deltaM,double ZSampleNormal,double MSampleNormal) const;
	void			MF_MakeSpace(int MF_MaxSamples);
	void			FillMassValues(void);
	void			FillZConditValues(void);

	const	double	OmegaM_;
	const	double	OmegaL_;
	const	double	OmegaK_;
	const	double	OmegaBar_;
	const	double	Sigma8_;
	const	double	CosmographMaxZ_;
	MassFunctionsType	MassFunction_;
	double			Step_;
	double			piv_;
	int				LastIndex_;

	int				MF_MaxSamples_;
	double			MF_MinMass_;
	double			MF_Log10MinMass_;
	double			MF_MaxMass_;
	double			MF_MaxZ_;
	double			DeltaZ_;
	double			DeltaM_;

	double			MF_HelpZmin_;
	double			MF_HelpZmax_;

	std::vector<double> Coefs_;

	MF_PrM_type							*MF_PrM_;
	MF_PrZM_type						*MF_PrZM_;
	SimpsonIntegrator<double,Cosmology>	*MF_Integrator_;
};

//
template<typename DATA_EL,typename DATA_ACCESSOR,typename WEIGHT_ACCESSOR>
class HistogramMaker
{
public:

	typedef typename	std::vector<DATA_EL>::const_iterator	vIterType;

	typedef HistogramMaker<DATA_EL,DATA_ACCESSOR,WEIGHT_ACCESSOR> meType;

	HistogramMaker(const typename std::vector<DATA_EL>& Data,DATA_ACCESSOR dataAccessor,WEIGHT_ACCESSOR weightAccessor,
		int NBins=40)
		:dataAccessor_(dataAccessor),weightAccessor_(weightAccessor),TotalWeight_(0.0),NBins_(NBins)
	{
		typename std::vector<DATA_EL>::const_iterator	beg(Data.begin());
		typename std::vector<DATA_EL>::const_iterator	const end(Data.end());

		for(;beg!=end;++beg)
		{
			histAux_.push_back(HistAux(beg));
		}

		std::sort(histAux_.begin(),histAux_.end(),HistAuxFunctor(dataAccessor));

		typename std::vector<HistAux>::const_iterator begAux(histAux_.begin());
		typename std::vector<HistAux>::const_iterator const endAux(histAux_.end());

		for(;begAux!=endAux;++begAux)
		{
			TotalWeight_ += weightAccessor(begAux->vAccess_);
		}

		MakeHistogram();
	}
//
	template<typename FUNCT>
	double GetExpectedFunctOverHPD(double HPD,FUNCT f_)
	{
		const double	WeightLimit(TotalWeight_ * HPD);
		double			AccWgh(0.0);
		double			AccVal(0.0);
		double			WCorrect;

		histIterType	beg(Histogram_.begin());
		histIterType	const end(Histogram_.end());

		for(;(beg!=end) && (AccWgh <= WeightLimit);++beg)
		{
			WCorrect = WeightLimit - AccWgh;
			histAuxIterType	begAux(beg->vBeg_);
			histAuxIterType	const endAux(beg->vEnd_);
			if(beg->Weight_ <= WCorrect)
			{WCorrect = 1.0;}
			else
			{WCorrect = (WCorrect / beg->Weight_);}

			for(;begAux!=endAux;++begAux)
			{
				const double w(WCorrect * weightAccessor_(begAux->vAccess_));
				AccVal	+=	(w * f_(dataAccessor_(begAux->vAccess_)));
				AccWgh	+=	w;
			}
		}
		return AccVal / AccWgh;
	}
//
	void GetHPDLimits(double HPD,double& HPD_min,double& HPD_max)
	{
		const double	WeightLimit(TotalWeight_ * HPD);

		HPD_min = INF;
		HPD_max = -INF;
		double	AccWgh(0.0);

		histIterType	beg(Histogram_.begin());
		histIterType	const end(Histogram_.end());

		for(;(beg!=end) && (AccWgh <= WeightLimit);++beg)
		{
			histAuxIterType	begAux(beg->vBeg_);
			histAuxIterType	const endAux(beg->vEnd_);

			for(;(begAux!=endAux) && (AccWgh <= WeightLimit);++begAux)
			{
				const double v(dataAccessor_(begAux->vAccess_));

				AccWgh	+=	weightAccessor_(begAux->vAccess_);
				if(HPD_min > v)	HPD_min = v;
				if(HPD_max < v)	HPD_max = v;
			}
		}
	}
//
private:
	struct HistAux
	{
		vIterType vAccess_;

		HistAux(vIterType& vAccess)
			:vAccess_(vAccess)
		{}
	};

	typedef typename std::vector<HistAux>::const_iterator	histAuxIterType;

	struct HistAtom
	{
		histAuxIterType	vBeg_;
		histAuxIterType	vEnd_;
		double			Weight_;
		double			Density_;
		double			Average_;
		int				Selected_;

		HistAtom()
			:Weight_(pwsNaN),Density_(pwsNaN),Selected_(0)
		{}
		inline bool operator<(const HistAtom& Rhs) const
		{
			return Density_ > Rhs.Density_;
		}
	};

	static std::string DebugGetAvDens(const HistAtom& h)
	{
		char buffer[1024];
		sprintf(buffer,"%9.6g,%9.6g,%9.6g",h.Average_,h.Density_,h.Weight_);
		return std::string(buffer);
	}

	typedef typename std::vector<HistAtom>::const_iterator	histIterType;

	struct HistAuxFunctor
	{
		DATA_ACCESSOR	dataAccessor_;
	
		inline bool operator()(const HistAux& LHSvAcc,const HistAux& RHSvAcc)
		{
			return dataAccessor_(LHSvAcc.vAccess_) < dataAccessor_(RHSvAcc.vAccess_);
		}

		HistAuxFunctor(DATA_ACCESSOR dataAccessor)
			:dataAccessor_(dataAccessor)
		{}
	};

	void MakeHistogram(void)
	{
		const double	wghBin(TotalWeight_ / static_cast<double>(NBins_));
		double			CurrTotalWeight(TotalWeight_);
		double			CurrNBins(static_cast<double>(NBins_));
		double			wghBinCurr;
		double			FstData,LastData, AverageData(0.0);
		double			CurrWgh(0.0);
	
		histAuxIterType begAux(histAux_.begin());
		histAuxIterType const endAux(histAux_.end());

		HistAtom	h;

		while(begAux!=endAux)
		{
			wghBinCurr	=	CurrTotalWeight / CurrNBins;
			h.vBeg_		=	begAux;
			LastData	=	FstData	=	dataAccessor_(h.vBeg_->vAccess_);
			h.Weight_	=	weightAccessor_(begAux->vAccess_);
			h.Average_	=	FstData * h.Weight_;
			++begAux;
toonarrowloop:
			for(;(begAux!=endAux) && (h.Weight_ < wghBinCurr);++begAux)
			{
				const double		w(weightAccessor_(begAux->vAccess_));
				const double		d(dataAccessor_(begAux->vAccess_));

				h.Weight_			+=	w;
				LastData			= d;
				h.Average_			+= (w * d);
			}
			// if interval singular or datum with large weight; rare
			if((begAux!=endAux) && (LastData == FstData))
			{
				wghBinCurr += wghBin;
				goto toonarrowloop;
			}
			
			CurrTotalWeight -= h.Weight_;
			CurrNBins		-= 1.0;

			if(CurrTotalWeight < 0.0) CurrTotalWeight = 0.0;
			if(CurrNBins < 0.5) CurrNBins = 1.0;

			// if just one single datum is left ; rare
			if(begAux == (endAux-1))
			{
				const double		w(weightAccessor_(begAux->vAccess_));
				const double		d(dataAccessor_(begAux->vAccess_));

				h.Weight_	+=	w;
				LastData	=	d;
				h.Average_	+= (w * d);

				++begAux;
			}
			h.vEnd_		=	begAux;
			h.Average_	/=	h.Weight_;
			// if last bin is singular; very rare
			if(LastData == FstData)
			{
				h.Density_ = (AverageData?(h.Weight_ * static_cast<double>(Histogram_.size())) / AverageData : h.Weight_);
			}
			else
			{
				const double l(LastData - FstData);
				h.Density_  = h.Weight_ / l;
				AverageData += l;
			}
			Histogram_.push_back(h);
		}

		std::sort(Histogram_.begin(),Histogram_.end());

		//DebugDumpHistogram();
	}

	void	DebugDumpHistogram(void)
	{
		static int counter(0);
		Zeus::DumpVector<HistAtom>("Histogram",Histogram_,counter,Histogram_.size()+1,meType::DebugGetAvDens);
		counter++;
	}

	std::vector<HistAtom>	Histogram_;
	std::vector<HistAux>	histAux_;
	DATA_ACCESSOR			dataAccessor_;
	WEIGHT_ACCESSOR			weightAccessor_;
	double					TotalWeight_;
	int						NBins_;
};
//
template<typename DATA_EL,typename DATA_ACCESSOR,typename WEIGHT_ACCESSOR>
class LinearRegression
{
public:
	typedef typename	std::vector<DATA_EL>::const_iterator		vIterType;
	typedef LinearRegression<DATA_EL,DATA_ACCESSOR,WEIGHT_ACCESSOR> meType;

	LinearRegression(const typename std::vector<DATA_EL>& Data,DATA_ACCESSOR dataAccessor,WEIGHT_ACCESSOR weightAccessor,
		double xUncertRatio= -1.0,double cLevel=0.95)
		:Data_(Data),dataAccessor_(dataAccessor),weightAccessor_(weightAccessor),
		Quantile_(SQRT2 * Dierfc(2.0 * (1.0-((cLevel<0.60)?0.60:((cLevel>0.9999)?0.9999:cLevel))))),
		xUncertRatio_(xUncertRatio),AvX_(0.0),AvY_(0.0),AvXX_(0.0),AvYY_(0.0),AvXY_(0.0),TotWght_(0.0),
		Slope_(-INF),Ord_(-INF),CorrCoef_(-INF)
	{
		typename std::vector<DATA_EL>::const_iterator	beg(Data_.begin());
		typename std::vector<DATA_EL>::const_iterator	const end(Data_.end());
		double											t;
		std::pair<double,double>						d;

		for(;beg != end;++beg)
		{
			t		= weightAccessor_(beg);
			d		= dataAccessor_(beg);

			AvX_		+= (t * d.first);
			AvY_		+= (t * d.second);
			AvXX_		+= (t * d.first * d.first);
			AvYY_		+= (t * d.second * d.second);
			AvXY_		+= (t * d.first * d.second);
			TotWght_	+=	t;
		}

		AvX_	/= TotWght_; AvY_	/= TotWght_; AvXX_	/= TotWght_;
		AvYY_	/= TotWght_; AvXY_	/= TotWght_;
		Slope_		= (AvXY_ - (AvX_ * AvY_))/(AvXX_ - (AvX_*AvX_));
		Ord_		= AvY_ - (Slope_ * AvX_);
		CorrCoef_	= (AvXY_ - (AvX_ * AvY_))/(std::sqrt(AvXX_ - (AvX_*AvX_))*std::sqrt(AvYY_ - (AvY_*AvY_)));
	}
//
	inline void GetParameters(double& Ord,double& Slope,double& CorrCoef) const
	{
		Slope = Slope_;Ord = Ord_;CorrCoef	= CorrCoef_;
	}
//
	void GetParamUncertain(double& ErrOrd,double& ErrSlope)
	{
		typename std::vector<DATA_EL>::const_iterator	beg(Data_.begin());
		typename std::vector<DATA_EL>::const_iterator	const end(Data_.end());
		double											t;
		std::pair<double,double>						d;
		double											residual(0.0),beta;
	
		if(xUncertRatio_ <= 0.0)
		{
			for(;beg != end;++beg)
			{
				t			=  weightAccessor_(beg);
				d			=  dataAccessor_(beg);
				residual	+= (t * (d.second- Ord_ - (Slope_ * d.first)) * (d.second- Ord_ - (Slope_ * d.first)));
			}
			residual	/= TotWght_;
			beta = std::sqrt(residual/(static_cast<double>(Data_.size()) * (AvXX_ -(AvX_*AvX_))));

			ErrSlope	= beta * Quantile_;
			ErrOrd		= ErrSlope * std::sqrt(AvXX_);
		}
		else
		{
			const double n((AvYY_ - (AvY_*AvY_)) - xUncertRatio_ * (AvXX_ - (AvX_*AvX_)));
			const double d(AvXY_ - (AvX_ * AvY_));

			ErrSlope	= ((n + std::sqrt((n*n) + (4.0 * xUncertRatio_ * d * d)))/(2.0 * d));
			ErrOrd		= (AvY_ - (ErrSlope * AvX_));	
		}
	}
//
private:
	double					AvX_,AvY_,AvXX_,AvYY_,AvXY_,TotWght_;
	DATA_ACCESSOR			dataAccessor_;
	WEIGHT_ACCESSOR			weightAccessor_;
	const double			xUncertRatio_;
	const double			Quantile_;
	double					Slope_;
	double					Ord_;
	double					CorrCoef_;
	const typename std::vector<DATA_EL>& Data_;
};

}
// End of namespace Zeus

#endif // PSNAKES_PHYSICSMATHH

