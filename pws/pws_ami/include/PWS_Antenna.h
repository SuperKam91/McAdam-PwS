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
#include "ZEUS_Factory.h"
#include "PWS_Globals.h"
//----------------------------------

class Antenna: public Zeus::FactoryObjType
{
protected:
	virtual void					do_Initialise(AntennaPropsType& props) = 0;
	virtual	std::complex<double>	do_GetAntennaAtFreq(double freqY,double freqX) const  = 0;
	virtual	void					do_GetAntennaBounds(double& freqY,double& freqX,double boundTol) const = 0;
public:
	typedef Zeus::FourierPlane<std::complex<double> >	BufferType;
	typedef BufferType::DataInnerType					BufferDataType;

	Antenna(const std::wstring& name,double PeriodArc,int NSamples,double Gain,double freq)
		: Zeus::FactoryObjType(name)
	{
		AntProps_.PeriodArc_	= PeriodArc;
		AntProps_.NSamples_		= NSamples;
		AntProps_.Gain_			= Gain;
		AntProps_.Freq_			= freq;
	}

	virtual void	Initialise(void);

	inline std::complex<double>	GetAntenaAtFreqRaw(double FreqY,double FreqX) const
	{return do_GetAntennaAtFreq(FreqY, FreqX);}

	inline const	AntennaPropsType&	GetAntenaProp(void) const
	{return AntProps_;}

	inline const Zeus::FourierPlane<std::complex<double> >& GetAntennaBufferRef(void) const
	{return buffer_;}

	inline Zeus::FourierPlane<std::complex<double> > GetAntennaBuffer(void) const
	{return buffer_;}

	inline Zeus::FourierPlane<std::complex<double> > CopyAntennaBuffer(void) const
	{return buffer_.Clone();}

	inline void ReleaseAntennaBuffer(void)
	{buffer_.Release();}

	virtual ~Antenna(void)
	{}
private:
	Antenna(const Antenna&);
	Antenna& operator=(const Antenna&);

	void 	FillEvalAntennaValues(Zeus::FourierPlane<std::complex<double> >& surf,std::complex<double>& Sum,std::complex<double>& SumSq);
	AntennaPropsType							AntProps_;
	BufferType									buffer_;
	double										SamplingFactor_;
};

#endif //ANTENNAH

