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

#ifndef FFTW_TRAITSH
#define FFTW_TRAITSH

#include "fftw3.h"
#include "ZEUS_StorageBaseManip.h"
#include "ZEUS_SingletonSerialize.h"
#include "FFTW_Globals.h"
#include "FFTW_Strings.h"


namespace fftw{

struct NoCheckBounds{
	void check(...) const {}
};

struct CheckBounds
{
	inline void check(int dimension,const int *v,int *dims) const
	{
		for(int i=0;i<dimension;i++){if((v[i]<0) || (v[i] >= dims[i])) errGeneralCheckBounds(ERROR_COD_FFTWOUTBOUNDS,ERROR_MSG_FFTWOUTBOUNDS);}
	}
	inline void check(int i,int dim1) const
	{
		if((i<0) || (i >= dim1)) errGeneralCheckBounds(ERROR_COD_FFTWOUTBOUNDS,ERROR_MSG_FFTWOUTBOUNDS);
	}
	inline void check(int i,int dim1,int j,int dim2) const
	{
		if((i<0) || (i >=dim1) || (j<0) || (j >=dim2)) errGeneralCheckBounds(ERROR_COD_FFTWOUTBOUNDS,ERROR_MSG_FFTWOUTBOUNDS);
	}
	inline void check(int i,int dim1,int j,int dim2,int k,int dim3) const
	{
		if((i<0) || (i >= dim1) || (j<0) || (j >= dim2) || (k<0) || (k >= dim3)) errGeneralCheckBounds(ERROR_COD_FFTWOUTBOUNDS,ERROR_MSG_FFTWOUTBOUNDS);
	}
	void errGeneralCheckBounds(int errNumber,wchar_t* msg) const
	{
		throw Zeus::libException(errNumber,msg,*this);
	}
};


template<typename T,template <class> class LockPolicy = DEFAULT_STORAGE_THREAD_POLICY>
class dataStoGeneral;

template<typename IN_T,typename OUT_T,template <class> class LockPolicy = DEFAULT_STORAGE_THREAD_POLICY>
class Device;

template
<
typename T,
int dimension,
class checkBoundsPolicy = DEFAULT_STORAGE_CHECK_BOUNDS,
template <class> class LockPolicy = DEFAULT_STORAGE_THREAD_POLICY 
>
class dataSto;

////////////////////////////////////////////////////////////////////////////////
// Functions to get the fast memory for the diferent models of storage
// Each overload maps to a function of fftw library
////////////////////////////////////////////////////////////////////////////////

extern void init_fftw(void);

struct GetFastMemObj
{
	inline void getFastMem(long Sz,double** ptrPtrDta)
	{*ptrPtrDta = static_cast<double*>(fftw_malloc(Sz));}
	inline void getFastMem(long Sz,float** ptrPtrDta)
	{*ptrPtrDta = static_cast<float*>(fftwf_malloc(Sz));}
	inline void getFastMem(long Sz,long double** ptrPtrDta)
	{*ptrPtrDta = static_cast<long double*>(fftwl_malloc(Sz));}
};
////////////////////////////////////////////////////////////////////////////////
// Functions to release the fast memory for the diferent models of storage
// Each overload maps to a function of fftw library
////////////////////////////////////////////////////////////////////////////////

struct releaseFastMemObj
{
	inline void releaseFastMem(double* ptrDta)
	{fftw_free(ptrDta);}
	inline void releaseFastMem(float* ptrDta)
	{fftwf_free(ptrDta);}
	inline void releaseFastMem(long double* ptrDta)
	{fftwl_free(ptrDta);}
};

struct freePlanObj
{
	inline void freePlan( fftw_plan plan) { fftw_free(plan);}
	inline void freePlan(fftwf_plan plan) {fftwf_free(plan);}
	inline void freePlan(fftwl_plan plan) {fftwl_free(plan);}
};

inline void exeFFTW( fftw_plan plan) { fftw_execute(plan);}
inline void exeFFTW(fftwf_plan plan) {fftwf_execute(plan);}
inline void exeFFTW(fftwl_plan plan) {fftwl_execute(plan);}

template<typename szT>
struct planTraits;

template<>
struct planTraits<double>
{
	typedef fftw_plan planType;
};

template<>
struct planTraits<float>
{
	typedef fftwf_plan planType;
};

template<>
struct planTraits<long double>
{
	typedef fftwl_plan planType;
};

template<typename T>
struct dataTraits
{
	typedef  T					Type;
	typedef  T					TypeSz;
};

template<typename T>
struct dataTraits<std::complex<T> >
{
	typedef typename std::complex <T>	Type;
	typedef T					TypeSz;
};


fftw_plan Do_fftw_plan_dft_1d(int n,std::complex<double> *in,std::complex<double> *out,int sign,unsigned int flags);
fftwf_plan Do_fftw_plan_dft_1d(int n,std::complex<float> *in,std::complex<float> *out,int sign,unsigned int flags);
fftwl_plan Do_fftw_plan_dft_1d(int n,std::complex<long double> *in,std::complex<long double> *out,int sign,unsigned int flags);
fftw_plan Do_fftw_plan_dft_2d(int nx,int ny,std::complex<double> *in,std::complex<double> *out,int sign,unsigned int flags);
fftwf_plan Do_fftw_plan_dft_2d(int nx,int ny,std::complex<float> *in,std::complex<float> *out,int sign,unsigned int flags);
fftwl_plan Do_fftw_plan_dft_2d(int nx,int ny,std::complex<long double> *in,std::complex<long double> *out,int sign,unsigned int flags);
fftw_plan Do_fftw_plan_dft_3d(int nx,int ny,int nz,std::complex<double> *in,std::complex<double> *out,int sign,unsigned int flags);
fftwf_plan Do_fftw_plan_dft_3d(int nx,int ny,int nz,std::complex<float> *in,std::complex<float> *out,int sign,unsigned int flags);
fftwl_plan Do_fftw_plan_dft_3d(int nx,int ny,int nz,std::complex<long double> *in,std::complex<long double> *out,int sign,unsigned int flags);
fftw_plan Do_fftw_plan_dft(int rank,const int *n,std::complex<double> *in,std::complex<double> *out,int sign,unsigned int flags);
fftwf_plan Do_fftw_plan_dft(int rank,const int *n,std::complex<float> *in,std::complex<float> *out,int sign,unsigned int flags);
fftwl_plan Do_fftw_plan_dft(int rank,const int *n,std::complex<long double> *in,std::complex<long double> *out,int sign,unsigned int flags);
fftw_plan Do_fftw_plan_dft_r2c_1d(int n,double *in,std::complex<double> *out,unsigned int flags);
fftwf_plan Do_fftw_plan_dft_r2c_1d(int n,float *in,std::complex<float> *out,unsigned int flags);
fftwl_plan Do_fftw_plan_dft_r2c_1d(int n,long double *in,std::complex<long double> *out,unsigned int flags);
fftw_plan Do_fftw_plan_dft_r2c_2d(int nx,int ny,double *in,std::complex<double> *out,unsigned int flags);
fftwf_plan Do_fftw_plan_dft_r2c_2d(int nx,int ny,float *in,std::complex<float> *out,unsigned int flags);
fftwl_plan Do_fftw_plan_dft_r2c_2d(int nx,int ny,long double *in,std::complex<long double> *out,unsigned int flags);
fftw_plan Do_fftw_plan_dft_r2c_3d(int nx,int ny,int nz,double *in,std::complex<double> *out,unsigned int flags);
fftwf_plan Do_fftw_plan_dft_r2c_3d(int nx,int ny,int nz,float *in,std::complex<float> *out,unsigned int flags);
fftwl_plan Do_fftw_plan_dft_r2c_3d(int nx,int ny,int nz,long double *in,std::complex<long double> *out,unsigned int flags);
fftw_plan Do_fftw_plan_dft_r2c(int rank,const int *n,double *in,std::complex<double> *out,unsigned int flags);
fftwf_plan Do_fftw_plan_dft_r2c(int rank,const int *n,float *in,std::complex<float> *out,unsigned int flags);
fftwl_plan Do_fftw_plan_dft_r2c(int rank,const int *n,long double *in,std::complex<long double> *out,unsigned int flags);
fftw_plan Do_fftw_plan_dft_c2r_1d(int n,std::complex<double> *in,double *out,unsigned int flags);
fftwf_plan Do_fftw_plan_dft_c2r_1d(int n,std::complex<float> *in,float *out,unsigned int flags);
fftwl_plan Do_fftw_plan_dft_c2r_1d(int n,std::complex<long double> *in,long double *out,unsigned int flags);
fftw_plan Do_fftw_plan_dft_c2r_2d(int nx,int ny,std::complex<double> *in,double *out,unsigned int flags);
fftwf_plan Do_fftw_plan_dft_c2r_2d(int nx,int ny,std::complex<float> *in,float *out,unsigned int flags);
fftwl_plan Do_fftw_plan_dft_c2r_2d(int nx,int ny,std::complex<long double> *in,long double *out,unsigned int flags);
fftw_plan Do_fftw_plan_dft_c2r_3d(int nx,int ny,int nz,std::complex<double> *in,double *out,unsigned int flags);
fftwf_plan Do_fftw_plan_dft_c2r_3d(int nx,int ny,int nz,std::complex<float> *in,float *out,unsigned int flags);
fftwl_plan Do_fftw_plan_dft_c2r_3d(int nx,int ny,int nz,std::complex<long double> *in,long double *out,unsigned int flags);
fftw_plan Do_fftw_plan_dft_c2r(int rank,const int *n,std::complex<double> *in,double *out,unsigned int flags);
fftwf_plan Do_fftw_plan_dft_c2r(int rank,const int *n,std::complex<float> *in,float *out,unsigned int flags);
fftwl_plan Do_fftw_plan_dft_c2r(int rank,const int *n,std::complex<long double> *in,long double *out,unsigned int flags);
fftw_plan Do_fftw_plan_dft_r2r_1d(int n,double *in,double *out,unsigned kind,unsigned int flags);
fftwf_plan Do_fftw_plan_dft_r2r_1d(int n,float *in,float *out,unsigned kind,unsigned int flags);
fftwl_plan Do_fftw_plan_dft_r2r_1d(int n,long double *in,long double *out,unsigned kind,unsigned int flags);
fftw_plan Do_fftw_plan_dft_r2r_2d(int nx,int ny,double *in,double *out,unsigned kindx,unsigned kindy,unsigned int flags);
fftwf_plan Do_fftw_plan_dft_r2r_2d(int nx,int ny,float *in,float *out,unsigned kindx,unsigned kindy,unsigned int flags);
fftwl_plan Do_fftw_plan_dft_r2r_2d(int nx,int ny,long double *in,long double *out,unsigned kindx,unsigned kindy,unsigned int flags);
fftw_plan Do_fftw_plan_dft_r2r_3d(int nx,int ny,int nz,double *in,double *out,unsigned kindx,unsigned kindy,unsigned kindz,unsigned int flags);
fftwf_plan Do_fftw_plan_dft_r2r_3d(int nx,int ny,int nz,float *in,float *out,unsigned kindx,unsigned kindy,unsigned kindz,unsigned int flags);
fftwl_plan Do_fftw_plan_dft_r2r_3d(int nx,int ny,int nz,long double *in,long double *out,unsigned kindx,unsigned kindy,unsigned kindz,unsigned int flags);
fftw_plan Do_fftw_plan_dft_r2r(int rank,const int *n,double *in,double *out,unsigned *kind,unsigned int flags);
fftwf_plan Do_fftw_plan_dft_r2r(int rank,const int *n,float *in,float *out,unsigned *kind,unsigned int flags);
fftwl_plan Do_fftw_plan_dft_r2r(int rank,const int *n,long double *in,long double *out,unsigned *kind,unsigned int flags);

#ifdef FFTW_MULTITHREADING
typedef Zeus::Multi_SingletonHolder<GetFastMemObj>::Type			GetMemType;
typedef Zeus::Multi_SingletonHolder<freePlanObj>::Type				freePlanType;
typedef Zeus::Multi_SingletonHolder<releaseFastMemObj>::Type		ReleaseMemType;
#else
typedef Zeus::Single_SingletonHolder<GetFastMemObj>::Type			GetMemType;
typedef Zeus::Single_SingletonHolder<freePlanObj>::Type				freePlanType;
typedef Zeus::Single_SingletonHolder<releaseFastMemObj>::Type		ReleaseMemType;
#endif //FFTW_MULTITHREADING




struct createPlan
{

void errGeneralCreatePlan(int errNumber,wchar_t* msg) const {throw Zeus::libException(errNumber,msg,*this);}


template <typename T,typename T1>
void CreatePlan(T& dev,int dim,Loki::Type2Type<std::complex<T1> >,Loki::Type2Type<std::complex<T1> >)
{
	int dimBuffer[FFTW_MAXDIMS];
	int dim0(dev.inData_->getDim0());
	int sign(dev.dir_?1:-1);
	typename planTraits<T1>::planType tempPlan;
	std::complex<T1> *outTemp(dev.inPlace_?dev.inData_->getRawPtr():dev.outData_->getRawPtr());
	if(dim > 1){dev.inData_->getHighDims(dimBuffer);dimBuffer[0] = dim0;}

	switch(dim){
		case 1:
			tempPlan = Do_fftw_plan_dft_1d(dim0,dev.inData_->getRawPtr(),outTemp ,sign,dev.optPlan_);
			dev.NormCte_ = dim0;
		break;
		case 2:
			tempPlan = Do_fftw_plan_dft_2d(dim0,dimBuffer[1],dev.inData_->getRawPtr(),outTemp,sign,dev.optPlan_);
			dev.NormCte_ = dim0 * dimBuffer[1];
		break;
		case 3:
			tempPlan = Do_fftw_plan_dft_3d(dim0,dimBuffer[1],dimBuffer[2],dev.inData_->getRawPtr(),outTemp,sign,dev.optPlan_);
			dev.NormCte_ = dim0 * dimBuffer[1]*dimBuffer[2];
		break;
		default:
			tempPlan = Do_fftw_plan_dft(dim,dimBuffer,dev.inData_->getRawPtr(),outTemp,sign,dev.optPlan_);
			dev.NormCte_ = dim0;
			for(int i=1;i<dim;++i) dev.NormCte_ *= dimBuffer[i]; 
		break;
	};
	if(tempPlan)
	{
		if(dev.plan_){(freePlanType::Instance())->freePlan(dev.plan_);}
		dev.plan_= tempPlan;
	}
	else errGeneralCreatePlan(ERROR_COD_FFTWINVALPLAN,ERROR_MSG_FFTWINVALPLAN);
}


template <typename T,typename T1>
void CreatePlan(T& dev,int dim,Loki::Type2Type<T1>,Loki::Type2Type<std::complex<T1> >)
{
	int dimBuffer[FFTW_MAXDIMS];
	int dim0(dev.inData_->getDim0());
	typename planTraits<T1>::planType tempPlan;
	std::complex<T1> *outTemp(dev.inPlace_?reinterpret_cast<std::complex<T1>* >(dev.inData_->getRawPtr()):dev.outData_->getRawPtr());
	if(dim > 1){dev.inData_->getHighDims(dimBuffer);dimBuffer[0] = dim0;}

	switch(dim){
		case 1:
			tempPlan = Do_fftw_plan_dft_r2c_1d(dim0,dev.inData_->getRawPtr(),outTemp,dev.optPlan_);
			dev.NormCte_ = dim0;
		break;
		case 2:
			tempPlan = Do_fftw_plan_dft_r2c_2d(dim0,dimBuffer[1],dev.inData_->getRawPtr(),outTemp,dev.optPlan_);
			dev.NormCte_ = dim0 * dimBuffer[1];
		break;
		case 3:
			tempPlan = Do_fftw_plan_dft_r2c_3d(dim0,dimBuffer[1],dimBuffer[2],dev.inData_->getRawPtr(),outTemp,dev.optPlan_);
			dev.NormCte_ = dim0 * dimBuffer[1]*dimBuffer[2];
		break;
		default:
			tempPlan = Do_fftw_plan_dft_r2c(dim,dimBuffer,dev.inData_->getRawPtr(),outTemp,dev.optPlan_);
			dev.NormCte_ = dim0;
			for(int i=1;i<dim;++i) dev.NormCte_ *= dimBuffer[i]; 
		break;
	};
	if(tempPlan)
	{
		if(dev.plan_){(freePlanType::Instance())->freePlan(dev.plan_);}
		dev.plan_= tempPlan;
	}
	else errGeneralCreatePlan(ERROR_COD_FFTWINVALPLAN,ERROR_MSG_FFTWINVALPLAN);
}

template <typename T,typename T1>
void CreatePlan(T& dev,int dim,Loki::Type2Type<std::complex<T1> >,Loki::Type2Type<T1> )
{
	int dimBuffer[FFTW_MAXDIMS];
	int dim0((2 * dev.inData_->getDim0()) - 1 - (dev.dir_ == true));
	typename planTraits<T1>::planType tempPlan;
	T1* outTemp(dev.inPlace_?reinterpret_cast<T1*>(dev.inData_->getRawPtr()):dev.outData_->getRawPtr());
	if(dim > 1){dev.inData_->getHighDims(dimBuffer);dimBuffer[0] = dim0;}


	switch(dim){
		case 1:
			tempPlan = Do_fftw_plan_dft_c2r_1d(dim0,dev.inData_->getRawPtr(),outTemp,dev.optPlan_);
			dev.NormCte_ = dim0;
		break;
		case 2:
			tempPlan = Do_fftw_plan_dft_c2r_2d(dim0,dimBuffer[1],dev.inData_->getRawPtr(),outTemp,dev.optPlan_);
			dev.NormCte_ = dim0 * dimBuffer[1];
		break;
		case 3:
			tempPlan = Do_fftw_plan_dft_c2r_3d(dim0,dimBuffer[1],dimBuffer[2],dev.inData_->getRawPtr(),outTemp,dev.optPlan_);
			dev.NormCte_ = dim0 * dimBuffer[1] * dimBuffer[2];
		break;
		default:
			tempPlan = Do_fftw_plan_dft_c2r(dim,dimBuffer,dev.inData_->getRawPtr(),outTemp,dev.optPlan_);
			dev.NormCte_ = dim0;
			for(int i=1;i<dim;++i) dev.NormCte_ *= dimBuffer[i]; 
		break;
	};
	if(tempPlan)
	{
		if(dev.plan_){(freePlanType::Instance())->freePlan(dev.plan_);}
		dev.plan_= tempPlan;
	}
	else errGeneralCreatePlan(ERROR_COD_FFTWINVALPLAN,ERROR_MSG_FFTWINVALPLAN);
}

template <typename T,typename T1>
void CreatePlan(T& dev,int dim,Loki::Type2Type<T1>,Loki::Type2Type<T1> )
{
	int dimBuffer[FFTW_MAXDIMS];
	typename planTraits<T1>::planType tempPlan;
	T1* outTemp(dev.inPlace_?dev.inData_->getRawPtr():dev.outData_->getRawPtr());
	int dim0(dev.inData_->getDim0());
	if(dim > 1){dev.inData_->getHighDims(dimBuffer);dimBuffer[0] = dim0;}


	switch(dim){
		case 1:
			tempPlan		=	Do_fftw_plan_dft_r2r_1d(dim0,dev.inData_->getRawPtr(),outTemp,dev.optSym_[0],dev.optPlan_);
			dev.NormCte_	=	getTransfDimLogicalSz(dev.optSym_[0],dim0);
		break;
		case 2:
			tempPlan		=	Do_fftw_plan_dft_r2r_2d(dim0,dimBuffer[1],dev.inData_->getRawPtr(),outTemp,dev.optSym_[0],dev.optSym_[1],dev.optPlan_);
			dev.NormCte_	=	getTransfDimLogicalSz(dev.optSym_[0],dim0)  * getTransfDimLogicalSz(dev.optSym_[1],dimBuffer[1]);
		break;
		case 3:
			tempPlan		=	Do_fftw_plan_dft_r2r_3d(dim0,dimBuffer[1],dimBuffer[2],dev.inData_->getRawPtr(),outTemp,dev.optSym_[0],dev.optSym_[1],dev.optSym_[2],dev.optPlan_);
			dev.NormCte_	=	getTransfDimLogicalSz(dev.optSym_[0],dim0)  * getTransfDimLogicalSz(dev.optSym_[1],dimBuffer[1]) * getTransfDimLogicalSz(dev.optSym_[2],dimBuffer[2]);
		break;
		default:
			tempPlan		=	Do_fftw_plan_dft_r2r(dim,dimBuffer,dev.inData_->getRawPtr(),outTemp,dev.optSym_,dev.optPlan_);
			dev.NormCte_	=	getTransfDimLogicalSz(dev.optSym_[0],dim0);
			for(int i=1;i<dim;++i) dev.NormCte_	*=	getTransfDimLogicalSz(dev.optSym_[i],dimBuffer[i]); 
		break;
	};
	if(tempPlan)
	{
		if(dev.plan_){(freePlanType::Instance())->freePlan(dev.plan_);}
		dev.plan_= tempPlan;
	}
	else errGeneralCreatePlan(ERROR_COD_FFTWINVALPLAN,ERROR_MSG_FFTWINVALPLAN);
}

unsigned getTransfDimLogicalSz(unsigned opt,unsigned physicalSz)
{
	if((opt < 3) || (opt > 10)) errGeneralCreatePlan(ERROR_COD_FFTWOPTNOTSUP,ERROR_MSG_FFTWOPTNOTSUP);
	switch(opt)
	{
	case FFTW_REDFT00 : 
		return 2 * (physicalSz - 1);
	case FFTW_RODFT00 :
		return 2 * (physicalSz + 1);
	default:
		return 2 * physicalSz;
	};
}
};

#ifdef FFTW_MULTITHREADING
typedef Zeus::Multi_SingletonHolder<createPlan>::Type				createPlanType;
#else
typedef Zeus::Single_SingletonHolder<createPlan>::Type				createPlanType;
#endif //FFTW_MULTITHREADING

} //end of namespace fftw

#endif
