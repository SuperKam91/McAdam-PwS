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
#ifndef ANTENNAPIXWINDOWH
#define ANTENNAPIXWINDOWH

#include "ZEUS_PhysicsMath.h"
#include "MC_Antenna.h"

//---------------------------------

// TODO - This clas MUST be revised
class AntennaPixWindow : public MC_antenna
{

protected:
	virtual void			do_Initialise(AntennaPropsType& props)
	{
		props.AntId_		= MC_antenna::ANTENNA_WINDOW;
		props.IsCircSymm_	= false;
	}
	virtual	double			do_GetAntennaAtFreq(double freqY,double freqX) const
	{return Zeus::sinc<double>(freqY / SampleFreq_) * Zeus::sinc<double>(freqX / SampleFreq_);}
	virtual	void			do_GetAntennaBounds(double& freqY,double& freqX,double boundTol) const
	{freqY = SampleFreq_;freqX = SampleFreq_;}
public:
	AntennaPixWindow(double PeriodArc,int NSamples,double width)
		: MC_antenna(PeriodArc,NSamples,-1.0,-1.0), SampleFreq_(1.0/PeriodArc)
	{}
private:
	const	double	SampleFreq_;
};


#endif //ANTENNAPIXWINDOWH

