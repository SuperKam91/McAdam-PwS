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



//----------------------------------------------------------------------
#include "ZEUS_SingletonSerialize.h"
#include "FFTW_Traits.h"
#include "FFTW_Device.h"

//----------------------------------------------------------------------
namespace fftw{

void init_fftw(void)
{
	GetMemType::Initialize();
	freePlanType::Initialize();
	ReleaseMemType::Initialize();
	createPlanType::Initialize();
}

fftw_plan Do_fftw_plan_dft_1d(int n,std::complex<double> *in,std::complex<double> *out,int sign,unsigned int flags)
{return fftw_plan_dft_1d(n,reinterpret_cast<fftw_complex*>(in),reinterpret_cast<fftw_complex*>(out),sign,flags);}
//
fftw_plan Do_fftw_plan_dft_2d(int nx,int ny,std::complex<double> *in,std::complex<double> *out,int sign,unsigned int flags)
{return fftw_plan_dft_2d(nx,ny,reinterpret_cast<fftw_complex*>(in),reinterpret_cast<fftw_complex*>(out),sign,flags);}
//
fftw_plan Do_fftw_plan_dft_3d(int nx,int ny,int nz,std::complex<double> *in,std::complex<double> *out,int sign,unsigned int flags)
{return fftw_plan_dft_3d(nx,ny,nz,reinterpret_cast<fftw_complex*>(in),reinterpret_cast<fftw_complex*>(out),sign,flags);}
//
fftw_plan Do_fftw_plan_dft(int rank,const int *n,std::complex<double> *in,std::complex<double> *out,int sign,unsigned int flags)
{return fftw_plan_dft(rank,n,reinterpret_cast<fftw_complex*>(in),reinterpret_cast<fftw_complex*>(out),sign,flags);}
//
fftw_plan Do_fftw_plan_dft_r2c_1d(int n,double *in,std::complex<double> *out,unsigned int flags)
{return fftw_plan_dft_r2c_1d(n,in,reinterpret_cast<fftw_complex*>(out),flags);}
//
fftw_plan Do_fftw_plan_dft_r2c_2d(int nx,int ny,double *in ,std::complex<double> *out,unsigned int flags)
{return fftw_plan_dft_r2c_2d(nx,ny,in,reinterpret_cast<fftw_complex*>(out),flags);}
//
fftw_plan Do_fftw_plan_dft_r2c_3d(int nx,int ny,int nz,double *in ,std::complex<double> *out,unsigned int flags)
{return fftw_plan_dft_r2c_3d(nx,ny,nz,in,reinterpret_cast<fftw_complex*>(out),flags);}
//
fftw_plan Do_fftw_plan_dft_r2c(int rank,const int *n,double *in ,std::complex<double> *out,unsigned int flags)
{return fftw_plan_dft_r2c(rank,n,in,reinterpret_cast<fftw_complex*>(out),flags);}
//
fftw_plan Do_fftw_plan_dft_c2r_1d(int n,std::complex<double> *in,double *out,unsigned int flags)
{return fftw_plan_dft_c2r_1d(n,reinterpret_cast<fftw_complex*>(in),out,flags);}
//
fftw_plan Do_fftw_plan_dft_c2r_2d(int nx,int ny,std::complex<double> *in,double *out,unsigned int flags)
{return fftw_plan_dft_c2r_2d(nx,ny,reinterpret_cast<fftw_complex*>(in),out,flags);}
//
fftw_plan Do_fftw_plan_dft_c2r_3d(int nx,int ny,int nz,std::complex<double> *in,double *out,unsigned int flags)
{return fftw_plan_dft_c2r_3d(nx,ny,nz,reinterpret_cast<fftw_complex*>(in),out,flags);}
//
fftw_plan Do_fftw_plan_dft_c2r(int rank,const int *n,std::complex<double> *in,double *out,unsigned int flags)
{return fftw_plan_dft_c2r(rank,n,reinterpret_cast<fftw_complex*>(in),out,flags);}
//
fftw_plan Do_fftw_plan_dft_r2r_1d(int n,double *in,double *out,unsigned kind,unsigned int flags)
{return fftw_plan_r2r_1d(n,in,out,static_cast<fftw_r2r_kind>(kind),flags);}
//
fftw_plan Do_fftw_plan_dft_r2r_2d(int nx,int ny,double *in,double *out,unsigned kindx,unsigned kindy,unsigned int flags)
{return fftw_plan_r2r_2d(nx,ny,in,out,static_cast<fftw_r2r_kind>(kindx),static_cast<fftw_r2r_kind>(kindy),flags);}
//
fftw_plan Do_fftw_plan_dft_r2r_3d(int nx,int ny,int nz,double *in,double *out,unsigned kindx,unsigned kindy,unsigned kindz,unsigned int flags)
{return fftw_plan_r2r_3d(nx,ny,nz,in,out,static_cast<fftw_r2r_kind>(kindx),static_cast<fftw_r2r_kind>(kindy),static_cast<fftw_r2r_kind>(kindz),flags);}
//
fftw_plan Do_fftw_plan_dft_r2r(int rank,const int *n,double *in,double *out,unsigned *kind,unsigned int flags)
{return fftw_plan_r2r(rank,n,in,out,reinterpret_cast<fftw_r2r_kind*>(kind),flags);}
//




#ifdef FULFFTW3
fftwf_plan Do_fftw_plan_dft_1d(int n,std::complex<float> *in,std::complex<float> *out,int sign,unsigned int flags)
{return fftwf_plan_dft_1d(n,reinterpret_cast<fftwf_complex*>(in),reinterpret_cast<fftwf_complex*>(out),sign,flags);}
//
fftwl_plan Do_fftw_plan_dft_1d(int n,std::complex<long double> *in,std::complex<long double> *out,int sign,unsigned int flags)
{return fftwl_plan_dft_1d(n,reinterpret_cast<fftwl_complex*>(in),reinterpret_cast<fftwl_complex*>(out),sign,flags);}
//
fftwf_plan Do_fftw_plan_dft_2d(int nx,int ny,std::complex<float> *in,std::complex<float> *out,int sign,unsigned int flags)
{return fftwf_plan_dft_2d(nx,ny,reinterpret_cast<fftwf_complex*>(in),reinterpret_cast<fftwf_complex*>(out),sign,flags);}
//
fftwl_plan Do_fftw_plan_dft_2d(int nx,int ny,std::complex<long double> *in,std::complex<long double> *out,int sign,unsigned int flags)
{return fftwl_plan_dft_2d(nx,ny,reinterpret_cast<fftwl_complex*>(in),reinterpret_cast<fftwl_complex*>(out),sign,flags);}
//
fftwf_plan Do_fftw_plan_dft_3d(int nx,int ny,int nz,std::complex<float> *in,std::complex<float> *out,int sign,unsigned int flags)
{return fftwf_plan_dft_3d(nx,ny,nz,reinterpret_cast<fftwf_complex*>(in),reinterpret_cast<fftwf_complex*>(out),sign,flags);}
//
fftwl_plan Do_fftw_plan_dft_3d(int nx,int ny,int nz,std::complex<long double> *in,std::complex<long double> *out,int sign,unsigned int flags)
{return fftwl_plan_dft_3d(nx,ny,nz,reinterpret_cast<fftwl_complex*>(in),reinterpret_cast<fftwl_complex*>(out),sign,flags);}
//
fftwf_plan Do_fftw_plan_dft(int rank,const int *n,std::complex<float> *in,std::complex<float> *out,int sign,unsigned int flags)
{return fftwf_plan_dft(rank,n,reinterpret_cast<fftwf_complex*>(in),reinterpret_cast<fftwf_complex*>(out),sign,flags);}
//
fftwl_plan Do_fftw_plan_dft(int rank,const int *n,std::complex<long double> *in,std::complex<long double> *out,int sign,unsigned int flags)
{return fftwl_plan_dft(rank,n,reinterpret_cast<fftwl_complex*>(in),reinterpret_cast<fftwl_complex*>(out),sign,flags);}
//
fftwf_plan Do_fftw_plan_dft_r2c_1d(int n,float *in,std::complex<float> *out,unsigned int flags)
{return fftwf_plan_dft_r2c_1d(n,in,reinterpret_cast<fftwf_complex*>(out),flags);}
//
fftwl_plan Do_fftw_plan_dft_r2c_1d(int n,long double *in,std::complex<long double> *out,unsigned int flags)
{return fftwl_plan_dft_r2c_1d(n,in,reinterpret_cast<fftwl_complex*>(out),flags);}
//
fftwf_plan Do_fftw_plan_dft_r2c_2d(int nx,int ny,float *in ,std::complex<float> *out,unsigned int flags)
{return fftwf_plan_dft_r2c_2d(nx,ny,in,reinterpret_cast<fftwf_complex*>(out),flags);}
//
fftwl_plan Do_fftw_plan_dft_r2c_2d(int nx,int ny,long double *in,std::complex<long double> *out,unsigned int flags)
{return fftwl_plan_dft_r2c_2d(nx,ny,in,reinterpret_cast<fftwl_complex*>(out),flags);}
//
fftwf_plan Do_fftw_plan_dft_r2c_3d(int nx,int ny,int nz,float *in ,std::complex<float> *out,unsigned int flags)
{return fftwf_plan_dft_r2c_3d(nx,ny,nz,in,reinterpret_cast<fftwf_complex*>(out),flags);}
//
fftwl_plan Do_fftw_plan_dft_r2c_3d(int nx,int ny,int nz,long double *in,std::complex<long double> *out,unsigned int flags)
{return fftwl_plan_dft_r2c_3d(nx,ny,nz,in,reinterpret_cast<fftwl_complex*>(out),flags);}
//
fftwf_plan Do_fftw_plan_dft_r2c(int rank,const int *n,float *in ,std::complex<float> *out,unsigned int flags)
{return fftwf_plan_dft_r2c(rank,n,in,reinterpret_cast<fftwf_complex*>(out),flags);}
//
fftwl_plan Do_fftw_plan_dft_r2c(int rank,const int *n,long double *in,std::complex<long double> *out,unsigned int flags)
{return fftwl_plan_dft_r2c(rank,n,in,reinterpret_cast<fftwl_complex*>(out),flags);}
//
fftwf_plan Do_fftw_plan_dft_c2r_1d(int n,std::complex<float> *in,float *out,unsigned int flags)
{return fftwf_plan_dft_c2r_1d(n,reinterpret_cast<fftwf_complex*>(in),out,flags);}
//
fftwl_plan Do_fftw_plan_dft_c2r_1d(int n,std::complex<long double> *in,long double *out,unsigned int flags)
{return fftwl_plan_dft_c2r_1d(n,reinterpret_cast<fftwl_complex*>(in),out,flags);}
//
fftwf_plan Do_fftw_plan_dft_c2r_2d(int nx,int ny,std::complex<float> *in,float *out,unsigned int flags)
{return fftwf_plan_dft_c2r_2d(nx,ny,reinterpret_cast<fftwf_complex*>(in),out,flags);}
//
fftwl_plan Do_fftw_plan_dft_c2r_2d(int nx,int ny,std::complex<long double> *in,long double *out,unsigned int flags)
{return fftwl_plan_dft_c2r_2d(nx,ny,reinterpret_cast<fftwl_complex*>(in),out,flags);}
//
fftwf_plan Do_fftw_plan_dft_c2r_3d(int nx,int ny,int nz,std::complex<float> *in,float *out,unsigned int flags)
{return fftwf_plan_dft_c2r_3d(nx,ny,nz,reinterpret_cast<fftwf_complex*>(in),out,flags);}
//
fftwl_plan Do_fftw_plan_dft_c2r_3d(int nx,int ny,int nz,std::complex<long double> *in,long double *out,unsigned int flags)
{return fftwl_plan_dft_c2r_3d(nx,ny,nz,reinterpret_cast<fftwl_complex*>(in),out,flags);}
//
fftwf_plan Do_fftw_plan_dft_c2r(int rank,const int *n,std::complex<float> *in,float *out,unsigned int flags)
{return fftwf_plan_dft_c2r(rank,n,reinterpret_cast<fftwf_complex*>(in),out,flags);}
//
fftwl_plan Do_fftw_plan_dft_c2r(int rank,const int *n,std::complex<long double> *in,long double *out,unsigned int flags)
{return fftwl_plan_dft_c2r(rank,n,reinterpret_cast<fftwl_complex*>(in),out,flags);}
//
fftwf_plan Do_fftw_plan_dft_r2r_1d(int n,float *in,float *out,unsigned kind,unsigned int flags)
{return fftwf_plan_r2r_1d(n,in,out,static_cast<fftwf_r2r_kind>(kind),flags);}
//
fftwl_plan Do_fftw_plan_dft_r2r_1d(int n,long double *in,long double *out,unsigned kind,unsigned int flags)
{return fftwl_plan_r2r_1d(n,in,out,static_cast<fftwl_r2r_kind>(kind),flags);}
//
fftwf_plan Do_fftw_plan_dft_r2r_2d(int nx,int ny,float *in,float *out,unsigned kindx,unsigned kindy,unsigned int flags)
{return fftwf_plan_r2r_2d(nx,ny,in,out,static_cast<fftwf_r2r_kind>(kindx),static_cast<fftwf_r2r_kind>(kindy),flags);}
//
fftwl_plan Do_fftw_plan_dft_r2r_2d(int nx,int ny,long double *in,long double *out,unsigned kindx,unsigned kindy,unsigned int flags)
{return fftwl_plan_r2r_2d(nx,ny,in,out,static_cast<fftwl_r2r_kind>(kindx),static_cast<fftwl_r2r_kind>(kindy),flags);}
//
fftwf_plan Do_fftw_plan_dft_r2r_3d(int nx,int ny,int nz,float *in,float *out,unsigned kindx,unsigned kindy,unsigned kindz,unsigned int flags)
{return fftwf_plan_r2r_3d(nx,ny,nz,in,out,static_cast<fftwf_r2r_kind>(kindx),static_cast<fftwf_r2r_kind>(kindy),static_cast<fftwf_r2r_kind>(kindz),flags);}
//
fftwl_plan Do_fftw_plan_dft_r2r_3d(int nx,int ny,int nz,long double *in,long double *out,unsigned kindx,unsigned kindy,unsigned kindz,unsigned int flags)
{return fftwl_plan_r2r_3d(nx,ny,nz,in,out,static_cast<fftwl_r2r_kind>(kindx),static_cast<fftwl_r2r_kind>(kindy),static_cast<fftwl_r2r_kind>(kindz),flags);}
//
fftwf_plan Do_fftw_plan_dft_r2r(int rank,const int *n,float *in,float *out,unsigned *kind,unsigned int flags)
{return fftwf_plan_r2r(rank,n,in,out,reinterpret_cast<fftwf_r2r_kind*>(kind),flags);}
//
fftwl_plan Do_fftw_plan_dft_r2r(int rank,const int *n,long double *in,long double *out,unsigned *kind,unsigned int flags)
{return fftwl_plan_r2r(rank,n,in,out,reinterpret_cast<fftwl_r2r_kind*>(kind),flags);}
//
#endif //FULFFTW3
}
