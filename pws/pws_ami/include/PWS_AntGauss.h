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


//---------------------------------
#ifndef ANTENNAGAUSSIANH
#define ANTENNAGAUSSIANH

#include "ZEUS_PhysicsMath.h"
#include "PWS_Antenna.h"

//---------------------------------

class AntennaGaussian : public Antenna
{

protected:
	inline virtual void			do_Initialise(AntennaPropsType& props)
	{
		props.IsCircSymm_	= true;
		double tx(SigmaX_ * FWHM2SIGMA);
		double ty(SigmaY_ * FWHM2SIGMA);
		props.FWHM_			= (tx + ty) / 2.0 ;
	}
	inline virtual	std::complex<double>	do_GetAntennaAtFreq(double freqY,double freqX) const
	{
		return std::exp(CoefY_*freqY*freqY+CoefX_*freqX*freqX);
		//return 1.0;
	}
	inline virtual	void			do_GetAntennaBounds(double& freqY,double& freqX,double boundTol) const
	{
		freqY = std::sqrt(std::log(boundTol)/CoefY_);
		freqX = std::sqrt(std::log(boundTol)/CoefX_);
	}
public:
	AntennaGaussian(const std::wstring& name,double PeriodArc,int NSamples,double Gain,double freq,double FWHM_Y,double FWHM_X)
		: Antenna(name,PeriodArc,NSamples,Gain,freq),SigmaY_(FWHM_Y /FWHM2SIGMA),SigmaX_(FWHM_X /FWHM2SIGMA)
	{
		CoefY_	= -(PITIMES2SQ * SigmaY_ * SigmaY_);
		CoefX_	= -(PITIMES2SQ * SigmaX_ * SigmaX_);
	}

	virtual ~AntennaGaussian(void)
	{}

private:
	const	double	SigmaY_;
	const	double	SigmaX_;
	double			CoefY_;
	double			CoefX_;
};


#endif //ANTENNAGAUSSIANH

