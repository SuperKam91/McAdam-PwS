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


//-------------------------------------------------------------
#include <time.h>
#include "fitshandle.h"
#include "ZEUS_Debug.h"
//-------------------------------------------------------------

namespace	Zeus
{

template<typename T>
void DumpInOut_1d(const char * fname,int Sz,int Offset,const T * const out,T scaling,int fnumber)
{
    FILE		*f;
	char		buffer[4096];
	static int	Counter(0);

#ifdef WIN32
	sprintf_s(buffer,200,DEBUGPATH "%s_%05d.ext",fname,Counter);
#else
	sprintf(buffer,DEBUGPATH "%s_%05d.ext",fname,Counter);
#endif
	f = fopen(buffer,"wt");
    fprintf(f,"%f",static_cast<double>(out[Offset] * scaling));
    for(int i=1;i<Sz;i++)
	{
		fprintf(f,",%f",static_cast<double>(out[Offset + i]*scaling));
	}
    fprintf(f,"\n");
    fclose(f);
	Counter++;
}

template void DumpInOut_1d(const char * fname,int Sz,int Offset,const float * const out,float scaling,int fnumber);
template void DumpInOut_1d(const char * fname,int Sz,int Offset,const double * const out,double scaling,int fnumber);
template void DumpInOut_1d(const char * fname,int Sz,int Offset,const long double * const out,long double scaling,int fnumber);


template<typename T>
void DumpInOut_2d(const char * fname,int YSz,int XSz,int ptMetric,int Offset,const T *out,T scaling)
{
    FILE *f;
	char		buffer[4096];
	static int	Counter(0);

#ifdef WIN32
	sprintf_s(buffer,200,DEBUGPATH "%s_%05d.ext",fname,Counter);
#else
	sprintf(buffer,DEBUGPATH "%s_%05d.ext",fname,Counter);
#endif
    f = fopen(buffer,"wt");

	T TmpValue;
	out += Offset;
    for(int j=0;j<YSz;j++){
      TmpValue = *(out + (j * ptMetric))*scaling;
      fprintf(f,"%7.6g",static_cast<double>(TmpValue));
      for(int i=1;i<XSz;i++){
            TmpValue = *(out + (j * ptMetric) + i)*scaling;
            fprintf(f," %7.6g",static_cast<double>(TmpValue));
      }
      fprintf(f,"\n");
    }
    fclose(f);
	Counter++;
}

template void DumpInOut_2d(const char * fname,int YSz,int XSz,int ptMetric,int Offset,const float *out,float scaling);
template void DumpInOut_2d(const char * fname,int YSz,int XSz,int ptMetric,int Offset,const double *out,double scaling);

template<typename T>
void DumpInOutFits_2d(const char * fname,int YSz,int XSz,int ptMetric,int Offset,const T *out,T scaling)
{
	fitshandle	*stream_(new fitshandle());
	char		buffer[4096];
	static int	Counter(0);

#ifdef WIN32
	sprintf_s(buffer,4096,"!" DEBUGPATH "%s_%05d.fits",fname,Counter);
#else
	sprintf(buffer,"!" DEBUGPATH "%s_%05d.fits",fname,Counter);
#endif
	stream_->create(std::string(buffer));

	arr2<T>		fitsArr(YSz, XSz);

	out += Offset;
	for(int j=0;j<YSz;j++)
	{
		for(int i=0;i<XSz;i++)
		{
			fitsArr(j, i) = (*((out+i) + (j * ptMetric))) *scaling;
		}
    }
	stream_->insert_image(planckType<T>(), fitsArr);
	stream_->close();
	delete stream_;
	Counter++;
}

template void DumpInOutFits_2d(const char * fname,int YSz,int XSz,int ptMetric,int Offset,const float *out,float scaling);
template void DumpInOutFits_2d(const char * fname,int YSz,int XSz,int ptMetric,int Offset,const double *out,double scaling);

template<typename T>
void DumpInOut_2d(const char * fname,int YSz,int XSz,int ptMetric,int Offset,const std::complex<T> *out,T scaling,int type)
{
    FILE *f;
	char		buffer[4096];
	static int	Counter(0);

#ifdef WIN32
	sprintf_s(buffer,200,DEBUGPATH "%s_%05d.ext",fname,Counter);
#else
	sprintf(buffer,DEBUGPATH "%s_%05d.ext",fname,Counter);
#endif
    f = fopen(buffer,"wt");
    T TmpValue;
	out += Offset;

    for(int j=0;j<YSz;j++){
		if(type) TmpValue = std::abs(*(out + (j * ptMetric)))*scaling;
		else TmpValue = std::arg(*(out + (j * ptMetric)));
		fprintf(f,"%7.6g",static_cast<double>(TmpValue));
		for(int i=1;i<XSz;i++){
			if(type) TmpValue = std::abs(*(out + (j * ptMetric) + i))*scaling;
			else TmpValue = std::arg(*(out + (j * ptMetric) + i));
            fprintf(f," %7.6g",static_cast<double>(TmpValue));
		}
		fprintf(f,"\n");
    }
    fclose(f);
	Counter++;
}

template void DumpInOut_2d(const char * fname,int YSz,int XSz,int ptMetric,int Offset,const std::complex<float> *out,float scaling,int type);
template void DumpInOut_2d(const char * fname,int YSz,int XSz,int ptMetric,int Offset,const std::complex<double> *out,double scaling,int type);

} // end of namespace Zeus

/*
	{
		std::vector<double> tDummy0(PriorsThetaMarginal,PriorsThetaMarginal+PriorsThetaMarginalSIZE);
		Zeus::MarginalPrior*  tMarginal(new Zeus::MarginalPrior(tDummy0));
		std::vector<double> tDummy(PriorsYcondTheta,PriorsYcondTheta+PriorsYcondThetaSIZE);
		std::vector<double> tDummyPDf(PriorsCondPDF,PriorsCondPDF+PriorsCondPDFSIZE);
		Zeus::ConditionalPrior*    tConditional(new Zeus::ConditionalPrior(tDummy,tDummyPDf,(PriorsYcondThetaSIZE / PriorsThetaMarginalSIZE)));
		tMarginal->SetSlaveObject(tConditional);
		Zeus::LArr2D<double> myData(65536,256);
		const double stepY((0.2-0.0002)/256);
		const double stepT((40.0-0.5)/256);
		double Yt,theta;
		for(int i=0;i<256;++i)
		{
			theta = (tMarginal->GetPDF(0.5 + i * stepT));
			for(int j=0;j<256;++j)
			{
				Yt = (tConditional->GetPDF(0.0002 + j * stepY));
				myData(j,i) = theta * Yt;
			}
		}
		Zeus::DumpInOut_2d("ProbData",256,256,256,0,myData.begin(),1.0);
	}

*/

