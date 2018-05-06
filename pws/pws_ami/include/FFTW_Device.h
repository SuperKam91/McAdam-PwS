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

#ifndef FFTW_deviceH
#define FFTW_deviceH
//---------------------------------------------------------------------------

#include "FFTW_Storage.h"


namespace fftw
{
template<typename inT,typename outT>
int calcN_Out(int Ninp);

template<typename inT>
int calcN_Out(int Ninp,Loki::Type2Type<inT>,Loki::Type2Type<inT>,bool even)
{return Ninp;}

template<typename T>
int calcN_Out(int Ninp,Loki::Type2Type<std::complex<T> >,Loki::Type2Type<T>,bool even)
{return 2 * Ninp - 1 - (even == true);}

template<typename T>
int calcN_Out(int Ninp,Loki::Type2Type<T>,Loki::Type2Type<std::complex<T> >,bool even)
{return (((int)(Ninp/2)) + 1);}

template
<
template <class> class LockPolicy
>
struct DeviceInterface
{
	typedef typename Zeus::ObjHandle<Zeus::HandleStorageBase<LockPolicy> >::Type HandleType;

	virtual void		IntAttachIn(HandleType inData,bool replanNow)=0;
	virtual void		IntAttachOut(HandleType outData,bool replanNow)=0;
	virtual void		IntDetachIn(void)=0;
	virtual void		IntDetachOut(void)=0;
	virtual void		IntDoFFT(bool normalize = false)=0;
	virtual void		IntDoNormalize(bool twice)=0;
	virtual void		IntDoPlan(void)=0;
	DeviceInterface() {}
	virtual ~DeviceInterface(void) throw() {}
};

#ifdef WIN32
#define _FFTWDEVICETYPENAME		typename
#else
#define _FFTWDEVICETYPENAME		
#endif


template
<
typename IN_T,
typename OUT_T,
template <class> class LockPolicy
>
class Device
	: public Zeus::HandleStorageBase<LockPolicy>
{
public:
	typedef Device<IN_T,OUT_T,LockPolicy>								meType;
	typedef IN_T*														inRawPtrType;
	typedef OUT_T*														outRawPtrType;
	typedef typename Zeus::ObjHandle<IN_T>::Type						inDataType;
	typedef typename Zeus::ObjHandle<OUT_T>::Type						outDataType;
	typedef typename IN_T::szType										AtomSzType;
	typedef typename IN_T::dataType										inAtomType;
	typedef typename OUT_T::dataType									outAtomType;
	typedef typename planTraits<AtomSzType>::planType					planType;
	typedef Zeus::HandleStorageBase<LockPolicy>							BaseType;
	typedef typename  BaseType::LockerType								LockerType;

#ifdef __GNUC__
	typedef typename BaseType::template VolatileType<planType>::Type		planVolType;
	typedef typename BaseType::template VolatileType<long>::Type			longVolType;
	typedef typename BaseType::template VolatileType<unsigned>::Type		unsignedVolType;
	typedef typename BaseType::template VolatileType<bool>::Type			boolVolType;
#else
	typedef _FFTWDEVICETYPENAME BaseType::VolatileType<planType>::Type		planVolType;
	typedef _FFTWDEVICETYPENAME	BaseType::VolatileType<long>::Type			longVolType;
	typedef _FFTWDEVICETYPENAME BaseType::VolatileType<unsigned>::Type		unsignedVolType;
	typedef _FFTWDEVICETYPENAME BaseType::VolatileType<bool>::Type			boolVolType;
#endif

	friend struct createPlan;

	static	meType* create(const inDataType& inData,const outDataType& outData,unsigned int options, bool inPlace,bool dir)
	{return new meType(inData,outData,options,inPlace,dir);}
	static	meType* create(const inDataType& inData,const outDataType& outData,unsigned int options, bool inPlace,unsigned int* optSyms)
	{return new meType(inData,outData,options,inPlace,optSyms);}

	inDataType	attachIn(const inDataType& inData, bool rePlanNow=true)
	{
		if(!Zeus::Check_me(inData)) errGeneral(ERROR_COD_FFTWBADHANDLE,ERROR_MSG_FFTWBADHANDLE);
		inDataType Temp;
		{
			LockerType locker(BaseType::getLocker());
			Temp = inData_;
			inData_ = inData;
			needRePlan_ = true;
		}
		if(rePlanNow) DoPlan(); 
		return Temp;
	}
	outDataType	attachOut(const outDataType& outData, bool rePlanNow=true)
	{
		if(!Zeus::Check_me(outData)) errGeneral(ERROR_COD_FFTWBADHANDLE,ERROR_MSG_FFTWBADHANDLE);
		outDataType Temp;
		{
			LockerType locker(BaseType::getLocker());
			Temp = outData_;
			outData_ = outData;
			needRePlan_ = true;
		}
		if(rePlanNow) DoPlan();
		return Temp;
	}
	inDataType	detachIn(void)
	{
		LockerType locker(BaseType::getLocker());
		inDataType Temp(inData_);
		needRePlan_ = true;
		if(Zeus::Check_me(inData_)) Reset(inData_,0);
		return Temp;
	}
	outDataType	detachOut(void)
	{
		LockerType locker(BaseType::getLocker());
		outDataType Temp(outData_);
		needRePlan_ = true;
		if(Zeus::Check_me(outData_)) Reset(outData_,0);
		return Temp;
	}
	void	DoFFT(void)
	{
		int errInd(0);
		LockerType locker(Zeus::HandleStorageBase<LockPolicy>::getLocker());
		if(!(Zeus::Check_me(inData_) && Zeus::Check_me(outData_)))
			errGeneral(ERROR_COD_FFTWBADHANDLE,ERROR_MSG_FFTWBADHANDLE);
		if(needRePlan_)
		{
			if(!(optPlan_ & FFTW_ESTIMATE)) errGeneral(ERROR_COD_FFTWINDATOVWR,ERROR_MSG_FFTWINDATOVWR);
			else DoPlan();
		}
		if(!plan_) errGeneral(ERROR_COD_FFTWINVALPLAN,ERROR_MSG_FFTWINVALPLAN);
		if(inPlace_)
		{
			if(inData_->getNRefs() != 1) errGeneral(ERROR_COD_FFTWINDATAREF,ERROR_MSG_FFTWINDATAREF);
			if(!(inData_->HasMemory())) errGeneral(ERROR_COD_FFTWNOTINDATA,ERROR_MSG_FFTWNOTINDATA);
			if(outData_->HasMemory()) errGeneral(ERROR_COD_FFTWOUTDATALL,ERROR_MSG_FFTWOUTDATALL);
		}
		exeFFTW(plan_);
		if(inPlace_)
		{
			typename OUT_T::LockerType lockOut(outData_->getLocker());
			if(outData_->HasMemory()) errGeneral(ERROR_COD_FFTWOUTDATALL,ERROR_MSG_FFTWOUTDATALL);
			outData_->reattachPtr(inData_->detach(),inData_->getDim0());
		}
		return;
	}

	void	DoPlan(void)
	{
		if(!needRePlan_) return;
		LockerType locker(Zeus::HandleStorageBase<LockPolicy>::getLocker());
		if(!needRePlan_) return;
		if(!(Zeus::Check_me(inData_) && Zeus::Check_me(outData_)))
			errGeneral(ERROR_COD_FFTWBADHANDLE,ERROR_MSG_FFTWBADHANDLE);

		checkData();
		inData_->allocData();
		if(!inPlace_){outData_->allocData();}

		(createPlanType::Instance())->CreatePlan(*this,inData_->getNDims(),Loki::Type2Type<inAtomType>(),Loki::Type2Type<outAtomType>());
		needRePlan_ = false;
	}
	virtual ~Device(void) throw()
	{
		if(plan_){(freePlanType::Instance())->freePlan(plan_);}
		if(optSym_) delete [] optSym_;
	}
protected:
	Device(const inDataType& inData,const outDataType& outData,unsigned int optPlan, bool inPlace,bool dir=true)
		: inData_(inData),outData_(outData),plan_(0),NormCte_(0),optPlan_(optPlan),optSym_(0),needRePlan_(true),inPlace_(inPlace),dir_(dir)
	{
		LOKI_STATIC_CHECK((Loki::IsSameType<AtomSzType,typename OUT_T::szType>::value),THE_DATA_TYPES_ARE_DIFFERENT);
		LOKI_STATIC_CHECK((!((Loki::IsSameType<inAtomType,outAtomType>::value) && (Loki::IsSameType<inAtomType,AtomSzType>::value))),CAN_NOT_USE_THIS_CTOR_WITH_THIS_STORAGE)
		DoPlan();
	}
	Device(const inDataType& inData,const outDataType& outData,unsigned int optPlan, bool inPlace,unsigned int *optSym)
		:inData_(inData),outData_(outData),plan_(0),NormCte_(0),optPlan_(optPlan),optSym_(0),needRePlan_(true),inPlace_(inPlace)
	{
		int n;
		LOKI_STATIC_CHECK((Loki::IsSameType<AtomSzType,typename OUT_T::szType>::value),THE_DATA_TYPES_ARE_DIFFERENT);
		LOKI_STATIC_CHECK(((Loki::IsSameType<inAtomType,outAtomType>::value) && (Loki::IsSameType<inAtomType,AtomSzType>::value)),CAN_NOT_USE_THIS_CTOR_WITH_THIS_STORAGE);
		if(!(Zeus::Check_me(inData_) && Zeus::Check_me(outData_)))
			errGeneral(ERROR_COD_FFTWBADHANDLE,ERROR_MSG_FFTWBADHANDLE);
		optSym_ = new unsigned int[(n = inData->getNDims())];
		for(int i=0;i<n;++i) {optSym_[i] = optSym[i];}
		try{DoPlan();}
		catch(...){delete [] optSym_;throw;}
	}


private:
	void checkData(void) const
	{
		int errInd(0),devDimen(inData_->getNDims());
		long StoSz,StoNeeded;
		int temp[FFTW_MAXDIMS];
		int dimOut0;

		if(devDimen != outData_->getNDims()) goto checkData_error;
		errInd = 1;
		if(devDimen > 1){
			outData_->getHighDims(temp);
			for(int i=1;i<devDimen;++i)
			{
				if(inData_->dims_[i] != temp[i]) goto checkData_error;
			}
		}
		dimOut0 = calcN_Out(inData_->getDim0(),Loki::Type2Type<inAtomType>(),Loki::Type2Type<outAtomType>(),dir_);
		temp[0] = dimOut0;
		StoNeeded = CalcStoSize_Atoms(devDimen,temp,false,false) * sizeof(outAtomType);
		if(inPlace_){StoSz = inData_->getStorageSize(false);}
		else{StoSz = outData_->getStorageSize(false);}
		if(StoSz < StoNeeded)
		{
			if(inPlace_){
				errInd = 2;
				if(inData_->allocData(-1,1 /*change xtra*/)) goto checkData_error;
			}
			else{
				errInd = 3;
				if(outData_->allocData(dimOut0 /*change dim0*/ ,0)) goto checkData_error;
			}
		}
		return;		
checkData_error:
		switch(errInd){
			case 0: errGeneral(ERROR_COD_FFTWDIFDIMENS,ERROR_MSG_FFTWDIFDIMENS);
			case 1:	errGeneral(ERROR_COD_FFTWDIFLODIMS,ERROR_MSG_FFTWDIFLODIMS);
			case 2:	errGeneral(ERROR_COD_FFTWINDATASHT,ERROR_MSG_FFTWINDATASHT);
			case 3:	errGeneral(ERROR_COD_FFTWOUTDATSHT,ERROR_MSG_FFTWOUTDATSHT);
		}
	}
	void errGeneral(int errNumber,wchar_t* msg) const {throw Zeus::libException(errNumber,msg,*this);}

inDataType							inData_;
outDataType							outData_;
planVolType							plan_;
longVolType							NormCte_;
unsignedVolType						optPlan_;
unsignedVolType						*optSym_;
boolVolType							needRePlan_;
boolVolType							inPlace_;
boolVolType							dir_;

};

template
<
typename IN_T,
typename OUT_T,
template <class> class LockPolicy
>
class DeviceGeneral
	: public DeviceInterface<LockPolicy>
	, public Device<IN_T,OUT_T,LockPolicy>

{
public:	
	typedef DeviceInterface<LockPolicy>* 	DeviceGeneralPtrType;
	typedef Device<IN_T,OUT_T,LockPolicy>	BaseDeviceType;
	typedef DeviceInterface<LockPolicy>		BaseInterfaceType;

	static DeviceGeneralPtrType createGeneral(const typename BaseDeviceType::inDataType& inData,const typename BaseDeviceType::outDataType& outData,unsigned int optPlan, bool inPlace,bool dir=true)
	{return new DeviceGeneral(inData,outData,optPlan,inPlace,dir);}
	static DeviceGeneralPtrType createGeneral(const typename BaseDeviceType::inDataType& inData,const typename BaseDeviceType::outDataType& outData,unsigned int optPlan, bool inPlace,unsigned int *optSym)
	{return new DeviceGeneral(inData,outData,optPlan,inPlace,optSym);}
	
	virtual void		IntAttachOut(typename BaseInterfaceType::HandleType outData,bool replanNow)
	{
		if(!outRawPtr_){
			if(!(outRawPtr_ = dynamic_cast<typename BaseDeviceType::outRawPtrType>(GetImplRef(outData)))) 
				errDeviceGeneral(ERROR_COD_FFTWNOTDYNAMI,ERROR_MSG_FFTWNOTDYNAMI);
		}
		typename BaseDeviceType::outDataType outDataSmartPtr_(BaseDeviceType::outDataType::OP::Clone(outRawPtr_));
		Device<IN_T,OUT_T,LockPolicy>::attachOut(outDataSmartPtr_,replanNow);
	}

	virtual void	IntAttachIn(typename BaseInterfaceType::HandleType inData,bool replanNow)
	{

		if(!inRawPtr_){
			if(!(inRawPtr_ = dynamic_cast<typename BaseDeviceType::inRawPtrType>(GetImplRef(inData))))
				errDeviceGeneral(ERROR_COD_FFTWNOTDYNAMI,ERROR_MSG_FFTWNOTDYNAMI);
		}
		typename BaseDeviceType::inDataType inDataSmartPtr_(BaseDeviceType::inDataType::OP::Clone(inRawPtr_));
		Device<IN_T,OUT_T,LockPolicy>::attachIn(inDataSmartPtr_,replanNow);
	}

	virtual void		IntDetachIn(void)
	{
		Device<IN_T,OUT_T,LockPolicy>::detachIn();
	}
	virtual void		IntDetachOut(void)
	{
		Device<IN_T,OUT_T,LockPolicy>::detachOut();
	}
	virtual void		IntDoFFT(bool normalize = false)
	{
		Device<IN_T,OUT_T,LockPolicy>::DoFFT(normalize);
	}
	virtual void		IntDoNormalize(bool twice)
	{
		Device<IN_T,OUT_T,LockPolicy>::DoNormalize(twice);
	}
	virtual void		IntDoPlan(void)
	{
		Device<IN_T,OUT_T,LockPolicy>::DoPlan();
	}

	virtual ~DeviceGeneral(void) throw()
	{}

private:
	DeviceGeneral(const typename BaseDeviceType::inDataType& inData,const typename BaseDeviceType::outDataType& outData,unsigned int optPlan, bool inPlace,bool dir=true)
		: Device<IN_T,OUT_T,LockPolicy>(inData,outData,optPlan,inPlace,dir),inRawPtr_(0),outRawPtr_(0)
	{}
	DeviceGeneral(const typename BaseDeviceType::inDataType& inData,const typename BaseDeviceType::outDataType& outData,unsigned int optPlan, bool inPlace,unsigned int *optSym)
		: Device<IN_T,OUT_T,LockPolicy>(inData,outData,optPlan,inPlace,optSym),inRawPtr_(0),outRawPtr_(0)
	{}
	void errDeviceGeneral(int errNumber,wchar_t* msg) const {throw Zeus::libException(errNumber,msg,*this);}

typename BaseDeviceType::inRawPtrType		inRawPtr_;
typename BaseDeviceType::outRawPtrType		outRawPtr_;

};

template
<
template <class> class LockPolicy = DEFAULT_STORAGE_THREAD_POLICY
>
struct GeneralDeviceFactory
{
	typedef  DeviceInterface<LockPolicy>	DeviceGeneralType;
	typedef  DeviceGeneralType*				DeviceGeneralPtrType;

template<typename IN_T,typename OUT_T>
static DeviceGeneralPtrType createGeneral(const IN_T& inData,const OUT_T& outData,unsigned int optPlan, bool inPlace,bool dir=true)
{
	return DeviceGeneral<typename IN_T::ObjHandledType,typename OUT_T::ObjHandledType,LockPolicy>::createGeneral(inData,outData,optPlan,inPlace,dir);
}

template<typename IN_T,typename OUT_T>
static DeviceGeneralPtrType createGeneral(const IN_T& inData,const OUT_T& outData,unsigned int optPlan, bool inPlace,unsigned int *optSym)
{
	return DeviceGeneral<typename IN_T::ObjHandledType,typename OUT_T::ObjHandledType,LockPolicy>::createGeneral(inData,outData,optPlan,inPlace,optSym);
}

};

} // end of namespace fftw

#endif
