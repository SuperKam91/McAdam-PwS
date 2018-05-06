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


//------------------------------------------------
#ifndef SOURCEOBJECTSZNGAIH
#define SOURCEOBJECTSZNGAIH

#ifdef PWSMPI
#include <mpi.h>
#endif //PWSMPI

#include "ZEUS_Object.h"

#define MINPIXDIST						0.02
#define	SCALEPIX						0.04

#define MAXTABLEDIST					50.0
//#define MAXTABLEDIST					8.0

#define MAXTRSDIST						10.0

// PIXTABLESZ = MAXTABLEDIST/SCALEPIX  + 1

#define PIXTABLESZ						1251

#define MAX0PIXLARGE					12.0
#define MIN0PIXLARGE					2.0
#define SCALE0PIXLARGE					0.1				
#define TABLE0PIXLARGE					101

#define MIN0PIXSMALL					0.01
#define SCALE0PIXSMALL					0.02
#define TABLE0PIXSMALL					100

#ifdef PWSMPI
#define PARAMETERSIZE					3
#else
#define PARAMETERSIZE					4
#endif

#define NGAIPROFDATASIZE				(PARAMETERSIZE + PIXTABLESZ + TABLE0PIXLARGE + TABLE0PIXSMALL)
#define NGAIPROFVARDATASIZE				(3 + PIXTABLESZ + TABLE0PIXLARGE + TABLE0PIXSMALL)

//
// Canonical set
/*
#define GLOBAL_SZVIRIALRATIO_DEF					5.0
#define GLOBAL_SZPROFALPHA_DEF						1.0620
#define GLOBAL_SZPROFBETA_DEF						5.4807
#define GLOBAL_SZPROFGAMMA_DEF						0.3292
#define GLOBAL_SZPROF_C500_DEF						1.156
#define GLOBAL_SZPROFCY500CYR500_DEF				1.814 
#define PROFILE_VARCORRECTION						1.05875
*/
// Planck set
/*
#define GLOBAL_SZVIRIALRATIO_DEF					5.0
#define GLOBAL_SZPROFALPHA_DEF						1.33
#define GLOBAL_SZPROFBETA_DEF						4.13
#define GLOBAL_SZPROFGAMMA_DEF						0.31
#define GLOBAL_SZPROF_C500_DEF						1.81
#define GLOBAL_SZPROFCY500CYR500_DEF				2.103 
#define PROFILE_VARCORRECTION						3.56733
// Warning beta is < 5
*/

// Canonical *non-standard* set

#define GLOBAL_SZVIRIALRATIO_DEF					5.0
#define GLOBAL_SZPROFALPHA_DEF						1.0510
#define GLOBAL_SZPROFBETA_DEF						5.4905
#define GLOBAL_SZPROFGAMMA_DEF						0.3081
#define GLOBAL_SZPROF_C500_DEF						1.1773
#define GLOBAL_SZPROFCY500CYR500_DEF				1.79643
#define PROFILE_VARCORRECTION						1.0819

/*
// PIP Profiles paper
#define GLOBAL_SZVIRIALRATIO_DEF					5.0
#define GLOBAL_SZPROFALPHA_DEF						1.051
#define GLOBAL_SZPROFBETA_DEF						5.4905
#define GLOBAL_SZPROFGAMMA_DEF						0.3081
#define GLOBAL_SZPROF_C500_DEF						1.1773
#define GLOBAL_SZPROFCY500CYR500_DEF				1.796 
*/
/*
// Coma
#define GLOBAL_SZVIRIALRATIO_DEF					5.0
#define GLOBAL_SZPROFALPHA_DEF						1.21
#define GLOBAL_SZPROFBETA_DEF						3.66
#define GLOBAL_SZPROFGAMMA_DEF						0.0
#define GLOBAL_SZPROF_C500_DEF						2.35
#define GLOBAL_SZPROFCY500CYR500_DEF				2.08653 *Cylindrical*
#define GLOBAL_SZPROFCY500CYR500_DEF				2.54256 *spherical*
*/

/*
// CC set
#define GLOBAL_SZVIRIALRATIO_DEF					5.0
#define GLOBAL_SZPROFALPHA_DEF						1.2223
#define GLOBAL_SZPROFBETA_DEF						5.4905
#define GLOBAL_SZPROFGAMMA_DEF						0.7736
#define GLOBAL_SZPROF_C500_DEF						1.128
#define GLOBAL_SZPROFCY500CYR500_DEF				1.618
*/

/*
// MD set
#define GLOBAL_SZVIRIALRATIO_DEF					5.0
#define GLOBAL_SZPROFALPHA_DEF						1.4063
#define GLOBAL_SZPROFBETA_DEF						5.4905
#define GLOBAL_SZPROFGAMMA_DEF						0.3798
#define GLOBAL_SZPROF_C500_DEF						1.083
#define GLOBAL_SZPROFCY500CYR500_DEF				1.872
*/
/*
// A1413
#define GLOBAL_SZVIRIALRATIO_DEF					5.0
#define GLOBAL_SZPROFALPHA_DEF						0.69
#define GLOBAL_SZPROFBETA_DEF						5.49
#define GLOBAL_SZPROFGAMMA_DEF						0.191
#define GLOBAL_SZPROF_C500_DEF						0.90
#define GLOBAL_SZPROFCY500CYR500_DEF				2.248
*/
/*
// A1914
#define GLOBAL_SZVIRIALRATIO_DEF					5.0
#define GLOBAL_SZPROFALPHA_DEF						0.95
#define GLOBAL_SZPROFBETA_DEF						5.49
#define GLOBAL_SZPROFGAMMA_DEF						0.000
#define GLOBAL_SZPROF_C500_DEF						1.88
#define GLOBAL_SZPROFCY500CYR500_DEF				1.457
*/
/*
// A2034
#define GLOBAL_SZVIRIALRATIO_DEF					5.0
#define GLOBAL_SZPROFALPHA_DEF						1.72
#define GLOBAL_SZPROFBETA_DEF						5.49
#define GLOBAL_SZPROFGAMMA_DEF						0.000
#define GLOBAL_SZPROF_C500_DEF						1.84
#define GLOBAL_SZPROFCY500CYR500_DEF				1.333
*/
/*
// A2218
#define GLOBAL_SZVIRIALRATIO_DEF					5.0
#define GLOBAL_SZPROFALPHA_DEF						0.74
#define GLOBAL_SZPROFBETA_DEF						5.49
#define GLOBAL_SZPROFGAMMA_DEF						0.000
#define GLOBAL_SZPROF_C500_DEF						1.02
#define GLOBAL_SZPROFCY500CYR500_DEF				2.234
*/
/*
// A773
#define GLOBAL_SZVIRIALRATIO_DEF					5.0
#define GLOBAL_SZPROFALPHA_DEF						0.96
#define GLOBAL_SZPROFBETA_DEF						5.49
#define GLOBAL_SZPROFGAMMA_DEF						0.000
#define GLOBAL_SZPROF_C500_DEF						1.25
#define GLOBAL_SZPROFCY500CYR500_DEF				1.894
*/

namespace Zeus
{

class GNFW_Profile
{
	const	double									OneMinusGamma_;
	const	double									Alpha_;
	const	double									AlphaOver2_;
	const	double									MinusGammaOver2_;
	const	double									GammaMinusBetaOverAlpha_;	
	double											RcRaiseAlpha_;
	double											BetaNormalFast_,pivFast_;
	Zeus::GammaFuncts								GammaFunct_;

//
	Zeus::SimpsonIntegrator<double,GNFW_Profile>	*PixIntegrator_;
	Zeus::SimpsonIntegrator<double,GNFW_Profile>	*Pix0AreaIntegrator_;
	Zeus::SimpsonIntegrator<double,GNFW_Profile>	*Pix0AreaIntegratorFast_;
//
	inline double		Get0PixValueAuxFast(double x)
	{
		return GammaFunct_.betainc(pivFast_/(std::pow(std::cos(x),Alpha_) + pivFast_)) * std::cos(x);
	}	
//
	inline double		Get0PixValueAux(double Rc)
	{
		return Rc * GetPixValueAtRc(Rc);
	}	
//
	inline double		GetProfValue(double x)
	{
		return std::pow((1.0 + x*x),MinusGammaOver2_) * std::pow((1.0 + RcRaiseAlpha_ * std::pow((1.0 + x*x),AlphaOver2_)),GammaMinusBetaOverAlpha_);
	}
//
	inline void		SetRc(double Rc)
	{
		RcRaiseAlpha_	= std::pow(Rc,Alpha_);
	}
//
	inline void		SetThetaFast(double a)
	{
		pivFast_	= std::pow(a,Alpha_);
	}
//
public:
	inline	void		Initialise(void)
	{
		PixIntegrator_			= new Zeus::SimpsonIntegrator<double,GNFW_Profile>(10,100,this,&GNFW_Profile::GetProfValue,1.0e-2,true);
		Pix0AreaIntegrator_		= new Zeus::SimpsonIntegrator<double,GNFW_Profile>(10,100,this,&GNFW_Profile::Get0PixValueAux,1.0e-2,true);
		Pix0AreaIntegratorFast_	= new Zeus::SimpsonIntegrator<double,GNFW_Profile>(10,100,this,&GNFW_Profile::Get0PixValueAuxFast,1.0e-2,true);
	}
//
	inline double		GetPixValueAtRc(double Rc)
	{
		bool	Dummy;

		SetRc(Rc);
		return std::pow(Rc,OneMinusGamma_) *  PixIntegrator_->Integrate(0.0,10.0/Rc,Dummy);
	}
//
	inline	double		Get0PixValue(double a)
	{
		bool	Dummy;
		return  2.0 * (Pix0AreaIntegrator_->Integrate(0.001,a,Dummy) / (a * a));
	}
//
	inline	double		Get0PixValueFast(double a)
	{
		bool	Dummy;
		SetThetaFast(a);
		return  (BetaNormalFast_ * (Pix0AreaIntegratorFast_->Integrate(0.0,PIOVER2,Dummy)) / (a * a));
	}
//
	GNFW_Profile(double	Alpha,double Beta,double Gamma)
		:OneMinusGamma_(1.0 - Gamma),Alpha_(Alpha),AlphaOver2_(Alpha / 2.0),
		MinusGammaOver2_(-Gamma/2.0),GammaMinusBetaOverAlpha_((Gamma-Beta) / Alpha),
		RcRaiseAlpha_(-1.0),PixIntegrator_(0), Pix0AreaIntegrator_(0),Pix0AreaIntegratorFast_(0),GammaFunct_()
	{
		GammaFunct_.setbetaAB((3.0-Gamma)/Alpha,(Beta-3.0)/Alpha); 
		BetaNormalFast_ = (2.0*GammaFunct_.beta())/Alpha;
	}
//
	~GNFW_Profile(void)
	{
		delete PixIntegrator_;
		delete Pix0AreaIntegrator_;
		delete Pix0AreaIntegratorFast_;
	}
};

//------------------------------------------------


class SrcObjSZNgai: public SrcObject
{

struct pixelValuesType
{
	double	_alpha;
	double	_beta;
	double	_gamma;

	double	pixelValues_[PIXTABLESZ];
	double	pixel0Small_[TABLE0PIXSMALL];
	double	pixel0Large_[TABLE0PIXLARGE];
};

protected:
	inline virtual	void	do_Initialise(const SrcObjectPropsType& props,double boundTol)
	{
		boundTol_	= boundTol;
		filtP_		= props.ObjParams_;
		Rs_			= props.ObjParams_.RealScale_;
		profileIndex_ = -1;
	}
	virtual	double			do_GetObjAtCoord(double YCoord,double XCoord) const;
	virtual	void			do_GetObjBounds(double& YBound,double& XBound) const;
	inline virtual	void	do_ChangeObjParams(const ObjFilterRealParams& objP,double YShift,double XShift)
	{
		filtP_	= objP;
		YShift_	= YShift;XShift_ = XShift;
		Rs_		= objP.RealScale_;

		if(IsVarProfile())
		{profileIndex_ = filtP_.i0_;}
		else
		{profileIndex_ = -1;}
	}
public:
	SrcObjSZNgai(const std::wstring& name,double PeriodArc,int NSamples,const SZPS_ProfParamType&	SZPS_ProfParams,
		int ContextID,const std::wstring& DirInMasks,int sync_id,
		const SZ_ProfParamRangeType* SZPS_ProfParamsRange=NULL);
//
	virtual ~SrcObjSZNgai(void)
	{}
private:
	double	pixelValues_[PIXTABLESZ];
	double	pixel0Small_[TABLE0PIXSMALL];
	double	pixel0Large_[TABLE0PIXLARGE];

	std::vector<pixelValuesType>	Profiles_;

	double					boundTol_;
	double					PeriodArc_;
	double					YShift_;
	double					XShift_;
	double					Rs_;
	ObjFilterRealParams		filtP_;
	SZPS_ProfParamType		SZPS_ProfParams_;
	const SZ_ProfParamRangeType* SZPS_ProfParamsRange_;
	int						sync_id_;
	int						profileIndex_;

	void		FillGNFWProfileArrays(const SZ_ProfParamRangeType* SZPS_ProfParamsRange);
	void		StoreNewlyEvalCoefs(int ContextID,const std::wstring& DirInMasks);
	int			ReadGNFWProfileCoefs(int ContextID,const std::wstring& DirInMasks);
//
	inline	int	IsVarProfile(void) const
	{
		return (!(pwsIsNaN(filtP_.a0_))) && (!(pwsIsNaN(filtP_.a1_))) && (!(pwsIsNaN(filtP_.a2_))) && (!(pwsIsNaN(filtP_.a3_))) && (filtP_.i0_ >= 0) && !(Profiles_.empty());
	}
//
	inline static void errInvalidArgs(void)
	{
		throw Zeus::libException(ERRCOD_PWS_PARAMOUTOFRANGE,ERRMSG_PWS_PARAMOUTOFRANGE,L"SrcObjSZNgai");
	}
//
	inline static void errCannotFindCoeffs(void)
	{
		throw Zeus::libException(ERRCOD_PWS_CANTFINDFILE,L"Cannot read GNFW coefficients",L"SrcObjSZNgai");
	}

};

}

#endif //SOURCEOBJECTSZNGAIH

