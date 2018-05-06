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


//-------------------------------------
#ifndef MCFOURIERMACHINEH
#define MCFOURIERMACHINEH


#include "ZEUS_StorageBaseManip.h"
#include "ZEUS_WorkSpace.h"
#include <complex>
#include "FFTW_Traits.h"
#include "FFTW_Storage.h"
#include "FFTW_Device.h"

//-------------------------------------


namespace Zeus
{

class FourierMachine
{
	typedef fftw::dataSto<double,2>									FFT_InDataType;
	typedef fftw::dataSto<std::complex<double>,2>					FFT_OutDataType;
	typedef Zeus::ObjHandle<FFT_InDataType>::Type					FFT_DataInHandleType;
	typedef Zeus::ObjHandle<FFT_OutDataType>::Type					FFT_DataOutHandleType;
	typedef fftw::Device<FFT_InDataType,FFT_OutDataType>			FFT_DirectDeviceType;
	typedef fftw::Device<FFT_OutDataType,FFT_InDataType>			FFT_InverseDeviceType;
	typedef ObjHandle<FFT_DirectDeviceType>::Type				FFT_DirectDeviceHandleType;
	typedef ObjHandle<FFT_InverseDeviceType>::Type			FFT_InverseDeviceHandleType;

public:
	FourierMachine(int WsSz)
		:WsSz_(WsSz)
	{}
	void	Initialize(void);

	void	Real2Fourier(const Zeus::LArr2D<double>& data,Zeus::LArr2D<std::complex<double> >& spectrum);

	void	Fourier2Real(const Zeus::LArr2D<std::complex<double> >& spectrum,Zeus::LArr2D<double>& data);

	void	Fourier2Real(const Zeus::LArr2D<double>& spectrum,Zeus::LArr2D<double>& data);

	inline int	GetXComplexBufferSz(void) const
	{return F_OutDataBuffer_->getDim0();}
	inline int  GetRealBufferSz(void) const
	{return WsSz_;}

private:
	FourierMachine(const FourierMachine&);
	FourierMachine& operator=(const FourierMachine&);

	int									WsSz_;
	FFT_DataInHandleType				F_InDataBuffer_;
	FFT_DataOutHandleType				F_OutDataBuffer_;
	FFT_DirectDeviceHandleType			F_DirectDev_;
	FFT_InverseDeviceHandleType			F_InverseDev_;

};

}
#endif // MCFOURIERMACHINEH

