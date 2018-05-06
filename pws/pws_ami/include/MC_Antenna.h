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



//----------------------------------
#ifndef ANTENNAH
#define ANTENNAH

#include "ZEUS_WorkSpace.h"
//----------------------------------


class MC_antenna
{
public:
	enum	AntennaIDType {ANTENNA_GAUSSIAN=0,ANTENNA_WINDOW=1};

	struct AntennaPropsType{
		bool			IsCircSymm_;
		AntennaIDType	AntId_;
		double			PeriodArc_;
		int				NSamples_;
		int				BoundSampleY_;
		int				BoundSampleX_;
		double			ZeroFourMode_;
		double			Gain_;
		double			Freq_;
		double			IntSq_;
		double			Int_;
	};

protected:
	virtual void			do_Initialise(AntennaPropsType& props) = 0;
	virtual	double			do_GetAntennaAtFreq(double freqY,double freqX) const  = 0;
	virtual	void			do_GetAntennaBounds(double& freqY,double& freqX,double boundTol) const = 0;
public:
//: AmplNormCte_(PeriodArc*PeriodArc)
	MC_antenna(double PeriodArc,int	NSamples,double Gain,double freq)
		: AmplNormCte_(1.0)
	{
		AntProps_.PeriodArc_	= PeriodArc;
		AntProps_.NSamples_		= NSamples;
		AntProps_.Gain_			= Gain;
		AntProps_.Freq_			= freq;
	}

	void			Initialise(void);
	inline double	GetAntenaAtSample(int Sy,int Sx) const
	{
		if(Sx < 0) Sx = -Sx;

		if(((AntProps_.BoundSampleY_ >=0 ) || (AntProps_.BoundSampleX_ >= 0)) &&
			((std::abs(Sy) >= AntProps_.BoundSampleY_) || (Sx >= AntProps_.BoundSampleX_))
			) return 0;

		if(Sy < 0) Sy += (AntProps_.NSamples_ << 1); 
		return (buffer_.GetInnerData())(Sy,Sx);
	}
	inline double	GetAntenaAtFreqRaw(double FreqY,double FreqX) const
	{return do_GetAntennaAtFreq(FreqY, FreqX);}
	inline const	AntennaPropsType&	GetAntenaProp(void) const
	{return AntProps_;}

	inline Zeus::FourierPlane<double>& GetAntennaBufferRef(void)
	{return buffer_;}
	inline Zeus::FourierPlane<double> GetAntennaBuffer(void) const
	{return buffer_;}
	inline Zeus::FourierPlane<double> CopyAntennaBuffer(void) const
	{return buffer_;}

	virtual ~MC_antenna(void)
	{}
private:
	MC_antenna(const MC_antenna&);
	MC_antenna& operator=(const MC_antenna&);

	void 	FillEvalAntennaValues(Zeus::FourierPlane<double>& surf,double& Sum,double& SumSq);
	const double					AmplNormCte_;
	AntennaPropsType				AntProps_;
	Zeus::FourierPlane<double>		buffer_;
	double							SamplingFactor_;
};

#endif //ANTENNAH
