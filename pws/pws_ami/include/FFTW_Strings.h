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



#ifndef FFTW_STRINGS
#define FFTW_STRINGS
#include "FFTW_Globals.h"
#include "ZEUS_Strings.h"

// Compiling time messages

#define CAN_NOT_USE_THIS_CTOR_WITH_THIS_STORAGE					CAN_NOT_USE_THIS_CTOR_WITH_THIS_STORAGE
#define THE_DATA_TYPES_ARE_DIFFERENT							THE_DATA_TYPES_ARE_DIFFERENT
#define THE_TYPES_ARE_INCOMPATIBLE								THE_TYPES_ARE_INCOMPATIBLE
#define DIMENSIONS_CANT_BE_0									DIMENSIONS_CANT_BE_0
#define TO_MANY_DIMENSIONS										TO_MANY_DIMENSIONS_MAX_##FFTW_MAXDIMS

// Runtime messages

#define ERROR_COD_FFTWBASE_OFFSET								100

#define ERROR_COD_FFTWUNDEFINED									ERROR_COD_FFTWBASE_OFFSET
#define ERROR_MSG_FFTWUNDEFINED									L"Undefined error"

#define ERROR_COD_FFTWWRONGDIMS									ERROR_COD_FFTWBASE_OFFSET+1
#define ERROR_MSG_FFTWWRONGDIMS									L"Dimensions must be positive"

#define ERROR_COD_FFTWOUTBOUNDS									ERROR_COD_FFTWBASE_OFFSET+2
#define ERROR_MSG_FFTWOUTBOUNDS									L"Index out of bounds"

#define ERROR_COD_FFTWBADHANDLE									ERROR_COD_FFTWBASE_OFFSET+3
#define ERROR_MSG_FFTWBADHANDLE									L"Invalid handle"

#define ERROR_COD_FFTWDIFDIMENS									ERROR_COD_FFTWBASE_OFFSET+4
#define ERROR_MSG_FFTWDIFDIMENS									L"The number of dims must be equal"

#define ERROR_COD_FFTWDIFLODIMS									ERROR_COD_FFTWBASE_OFFSET+5
#define ERROR_MSG_FFTWDIFLODIMS									L"Lower dims of data must be equal in size"

#define ERROR_COD_FFTWINDATASHT									ERROR_COD_FFTWBASE_OFFSET+6
#define ERROR_MSG_FFTWINDATASHT									L"In data storage not enough"

#define ERROR_COD_FFTWOUTDATSHT									ERROR_COD_FFTWBASE_OFFSET+7
#define ERROR_MSG_FFTWOUTDATSHT									L"In data storage not enough"

#define ERROR_COD_FFTWINDATOVWR									ERROR_COD_FFTWBASE_OFFSET+8
#define ERROR_MSG_FFTWINDATOVWR									L"Input data will be destroyed. Plan earlier"

#define ERROR_COD_FFTWINVALPLAN									ERROR_COD_FFTWBASE_OFFSET+9
#define ERROR_MSG_FFTWINVALPLAN									L"Invalid plan"

#define ERROR_COD_FFTWINDATAREF									ERROR_COD_FFTWBASE_OFFSET+10
#define ERROR_MSG_FFTWINDATAREF									L"Can't dettach. Input data has references"

#define ERROR_COD_FFTWOUTDATALL									ERROR_COD_FFTWBASE_OFFSET+11
#define ERROR_MSG_FFTWOUTDATALL									L"Can't attach. Output data is allocated"

#define ERROR_COD_FFTWNDIMTOBIG									ERROR_COD_FFTWBASE_OFFSET+12
#define ERROR_MSG_FFTWNDIMTOBIG									L"Number of dimensions to big. Max allowed 16"

#define ERROR_COD_FFTWPLANNOTAV									ERROR_COD_FFTWBASE_OFFSET+13
#define ERROR_MSG_FFTWPLANNOTAV									L"There is no plan available for these parameters"

#define ERROR_COD_FFTWNOTINDATA									ERROR_COD_FFTWBASE_OFFSET+14
#define ERROR_MSG_FFTWNOTINDATA									L"Input data does not exists"

#define ERROR_COD_FFTWNDIMISNUL									ERROR_COD_FFTWBASE_OFFSET+15
#define ERROR_MSG_FFTWNDIMISNUL									L"The number of dimensions can't be null"

#define ERROR_COD_FFTWNDIMOUTBO									ERROR_COD_FFTWBASE_OFFSET+16
#define ERROR_MSG_FFTWNDIMOUTBO									L"The number of dimensions is outside of the allowed range"

#define ERROR_COD_FFTWNOTOUTDAT									ERROR_COD_FFTWBASE_OFFSET+17
#define ERROR_MSG_FFTWNOTOUTDAT									L"Output data does not exists"

#define ERROR_COD_FFTWMETHNOTEX									ERROR_COD_FFTWBASE_OFFSET+18
#define ERROR_MSG_FFTWMETHNOTEX									L"Method does not exists"

#define ERROR_COD_FFTWNOTDYNAMI									ERROR_COD_FFTWBASE_OFFSET+19
#define ERROR_MSG_FFTWNOTDYNAMI									L"This is not a FFTW dynamic class"

#define ERROR_COD_FFTWOPTNOTSUP									ERROR_COD_FFTWBASE_OFFSET+20
#define ERROR_MSG_FFTWOPTNOTSUP									L"Option not supported"

#endif //FFTW_STRINGS
