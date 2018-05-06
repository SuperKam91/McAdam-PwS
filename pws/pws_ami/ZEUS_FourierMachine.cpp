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

#include "ZEUS_FourierMachine.h"

//---------------------------------

namespace Zeus
{

void	FourierMachine::Initialize(void)
{
	F_InDataBuffer_		= FFT_InDataType::create(WsSz_,WsSz_);
	F_OutDataBuffer_	= FFT_OutDataType::create(1,WsSz_);
	F_DirectDev_		= FFT_DirectDeviceType::create(F_InDataBuffer_,F_OutDataBuffer_,FFTW_MEASURE | FFTW_DESTROY_INPUT,false,true);
	F_InverseDev_		= FFT_InverseDeviceType::create(F_OutDataBuffer_,F_InDataBuffer_,FFTW_MEASURE | FFTW_DESTROY_INPUT,false,true);
}

void	FourierMachine::Real2Fourier(const Zeus::LArr2D<double>& data,
									 Zeus::LArr2D<std::complex<double> >& spectrum)
{
	double									*fPiv(F_InDataBuffer_->getRawPtr());
	LArr2D<double>::const_iterator	piv(data.begin());
	LArr2D<double>::const_iterator	const end(piv + data.getSz());

	while(piv != end)
	{*fPiv++ = *piv++;}

	F_DirectDev_->DoFFT();

	std::complex<double>	*fOutPiv(F_OutDataBuffer_->getRawPtr());
	LArr2D<std::complex<double> >::iterator
		pivOut(spectrum.begin());
	LArr2D<std::complex<double> >::const_iterator	const
		endOut(pivOut + spectrum.getSz());

	while(pivOut != endOut)
	{*pivOut++ = *fOutPiv++;}
}

void	FourierMachine::Fourier2Real(const Zeus::LArr2D<std::complex<double> >& spectrum,
									 Zeus::LArr2D<double>& data)
{
	std::complex<double>	*fPiv(F_OutDataBuffer_->getRawPtr());
	LArr2D<std::complex<double> >::const_iterator
		piv(spectrum.begin());
	LArr2D<std::complex<double> >::const_iterator	const
		end(piv + spectrum.getSz());

	while(piv != end)
	{*fPiv++ = *piv++;}

	F_InverseDev_->DoFFT();

	double									*fPivOut(F_InDataBuffer_->getRawPtr());
	Zeus::LArr2D<double>::iterator			pivOut(data.begin());
	Zeus::LArr2D<double>::const_iterator	const endOut(pivOut + data.getSz());

	while(pivOut!=endOut)
	{*pivOut++ = *fPivOut++;}
}

void	FourierMachine::Fourier2Real(const Zeus::LArr2D<double>& spectrum,
									 Zeus::LArr2D<double>& data)
{
	std::complex<double>					*fPiv(F_OutDataBuffer_->getRawPtr());
	LArr2D<double>::const_iterator			piv(spectrum.begin());
	LArr2D<double>::const_iterator	const	end(piv + spectrum.getSz());

	while(piv != end)
	{*fPiv++ = std::complex<double>(*piv++,0.0);}

	F_InverseDev_->DoFFT();

	double									*fPivOut(F_InDataBuffer_->getRawPtr());
	Zeus::LArr2D<double>::iterator			pivOut(data.begin());
	Zeus::LArr2D<double>::const_iterator	const endOut(pivOut + data.getSz());

	while(pivOut!=endOut)
	{*pivOut++ = *fPivOut++;}
}

}


