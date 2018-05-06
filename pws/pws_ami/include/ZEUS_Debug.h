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



//---------------------------------------------------------------------------
#ifndef PSNAKES_DEBUGH
#define PSNAKES_DEBUGH

#include <stdio.h>
#include <string.h>
#include <complex>
#include <vector>


#ifdef WIN32
#define	DEBUGPATH	""
#else
#define	DEBUGPATH	""
#endif

namespace	Zeus
{

struct myDebugData
{
	double x0_;
	double x1_;
	double x2_;
	double x3_;
	double x4_;

	myDebugData(void)
		:x0_(0.0),x1_(0.0),x2_(0.0),x3_(0.0),x4_(0.0)
	{}
};

inline std::string myGetString(const myDebugData& d)
{
	char buffer[1024];

	sprintf(buffer,"%g,%g,%g,%g,%g",d.x0_,d.x1_,d.x2_,d.x3_,d.x4_);

	return std::string(buffer);	
}


template<typename T>
void	DumpVector(const char * fname,const std::vector<T>& data ,int Id,int sz,std::string (*extractFn)(const T& data))
{
    FILE *f;
	char buffer[1024];

	sprintf(buffer,"%s_%d_",fname,Id);

#ifdef WIN32
	sprintf_s(buffer+strlen(buffer),300,"%03d",static_cast<int>(((char*)&data)-((char*)0))+time(0));
#else
	sprintf(buffer,"%s%03d",fname,static_cast<int>(((char*)&data)-((char*)0))+time(0));
#endif
	strcat(buffer,".csv");
	f = fopen(buffer,"wt");
	typename std::vector<T>::const_iterator	piv(data.begin());	
	typename std::vector<T>::const_iterator	end(data.end());
	for(int i=0;(piv!=end) && (i < sz);++piv,++i)
	{
		fprintf(f,"%s\n",(extractFn(*piv)).c_str());
	}
    fprintf(f,"\n");
    fclose(f);
}

template<typename T>
void DumpInOut_2d(const char * fname,int YSz,int XSz,int ptMetric,int Offset,const T *out,T scaling);
template<typename T>
void DumpInOut_2d(const char * fname,int YSz,int XSz,int ptMetric,int Offset,const std::complex<T> *out,T scaling,int type);
template<typename T>
void DumpInOutFits_2d(const char * fname,int YSz,int XSz,int ptMetric,int Offset,const T *out,T scaling);
//
template<typename T>
void DumpInOut_1d(const char * fname,int Sz,int Offset,const T * const out,T scaling,int fnumber);

} // End of namespace Zeus

#endif // PSNAKES_DEBUGH

