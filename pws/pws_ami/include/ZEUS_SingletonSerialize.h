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



#ifndef SIGLETON_SERIALIZEH
#define SIGLETON_SERIALIZEH

#include <string>

#ifdef WIN32
#include <tchar.h>
#endif //WIN32

#include <complex>
#include <math.h>
#include <typeinfo>
#include <exception>
#include <locale>
#include "Threads.h"
#include "SmartPtr.h"
#include "TypeManip.h"
#include "static_check.h"
#include "ZEUS_Strings.h"
#include "ZEUS_Exceptions.h"


namespace Zeus
{

template<typename InnerHost>
class LockProxyNoThreading
{
public:
	typedef InnerHost*				Host_ptrInnerTy;
	typedef class Empty{} 			LockerTy;

	LockProxyNoThreading(...)
	{}
	~LockProxyNoThreading()
	{}
	Host_ptrInnerTy operator->() const
	{
		if(!obj_) OnNotInitialized();
		return obj_;
	}
	LockProxyNoThreading(const LockProxyNoThreading&)
	{}
	static void initPtr(Host_ptrInnerTy obj)  {obj_= obj;}
private:
	static Host_ptrInnerTy obj_;
	static void OnNotInitialized(void)
	{
		throw Zeus::libException(ERROR_COD_ZEUSSGNOTINIT,ERROR_MSG_ZEUSSGNOTINIT,L"static void LockProxyNoThreading::OnNotInitialized()");
	}

	LockProxyNoThreading& operator=(const LockProxyNoThreading&);
};

template<typename InnerHost>
class LockProxyMultiThreading
{
public:
	typedef typename Loki::ClassLevelLockable<InnerHost>		LockerClassTy;
	typedef typename LockerClassTy::Lock						LockerTy;
	typedef InnerHost*											Host_ptrInnerTy;

	LockProxyMultiThreading()
	{
		classLockerObj_.DoLock();
	}
	~LockProxyMultiThreading()
	{
		classLockerObj_.DoUnLock();
	}
	inline Host_ptrInnerTy operator->() const
	{
		if(!obj_) OnNotInitialized();
		return obj_;
	}
	inline LockProxyMultiThreading(const LockProxyMultiThreading&)
	{
		classLockerObj_.DoLock();
	}
	static void initPtr(Host_ptrInnerTy obj)  {obj_= obj;}
private:
	LockerClassTy				classLockerObj_;
	static void OnNotInitialized(void)
	{
		throw Zeus::libException(ERROR_COD_ZEUSSGNOTINIT,ERROR_MSG_ZEUSSGNOTINIT,L"static void LockProxyNoThreading::OnNotInitialized()");
	}
	static Host_ptrInnerTy		obj_;
	LockProxyMultiThreading& operator=(const LockProxyMultiThreading&);
};

template
<
	typename T,
	template <class> class LockProxy
>
class singleSubstrata
{
public:
	typedef T										InnerTy;
	typedef T*										ptrInnerTy;
	typedef singleSubstrata<T,LockProxy>			meTy;
	typedef LockProxy<T>							myLockerClassTy;
	typedef typename myLockerClassTy::LockerTy		myLockerTy;

	static  myLockerClassTy Instance(void)
	{
		if(!AccessPoint_)
		{
			myLockerTy locker;
			(void)locker; 
			if(!AccessPoint_)
			{
				if(destroyed_){OnDeadReference();}
				else{Create();}
			}
		}
		return myLockerClassTy();
	}
	static void Initialize(void)
	{
		if(!AccessPoint_)
		{
			myLockerTy locker;
			(void)locker; 
			if(!AccessPoint_)
			{
				if(destroyed_){OnDeadReference();}
				else{Create();}
			}
		}
	}

private:
	singleSubstrata()
	{}
    singleSubstrata(const singleSubstrata&);
    singleSubstrata operator=(const singleSubstrata&);
	static void Create()
	{
		static meTy theObject;

		AccessPoint_ = &theObject;
		myLockerClassTy::initPtr(&(AccessPoint_->InnerObj_));
	}
	static void OnDeadReference()
	{
		throw Zeus::libException(ERROR_COD_ZEUSBADHANDLE,ERROR_MSG_ZEUSBADHANDLE,L"static void singleSubstrata::OnDeadReference()");
	}
	virtual ~singleSubstrata()
	{
		myLockerTy locker;
		(void)locker; 
		AccessPoint_ = 0; destroyed_ = true;
	}
	InnerTy							InnerObj_;
    static  meTy*					AccessPoint_;   // Singleton access point
	static  bool					destroyed_;
};

template
<
	typename T,
	template <typename> class LockProxy
>
class singleSubstrataDyn
{
public:
	typedef T										InnerTy;
	typedef T*										ptrInnerTy;
	typedef singleSubstrataDyn<T,LockProxy>			meTy;
	typedef LockProxy<T>							myLockerClassTy;
	typedef typename myLockerClassTy::LockerTy		myLockerTy;

	static  myLockerClassTy Instance(void)
	{
		if(!AccessPoint_) {OnDeadReference(L"static void singleSubstrataDyn::Instance()");}
		return myLockerClassTy();
	}
	static ptrInnerTy Initialize(ptrInnerTy innerObj)
	{
		myLockerTy locker;
		(void)locker; 
		if(!AccessPoint_)
		{
			if(destroyed_){OnDeadReference(L"static void singleSubstrataDyn::Initialize()");}
			else{Create();}
		}
		ptrInnerTy temp(ptrInnerObj_);
		myLockerClassTy::initPtr(innerObj);
		ptrInnerObj_ = innerObj;
		destroyed_ = false;
		return temp;
	}
	static void Release(void)
	{
		myLockerTy locker;
		(void)locker; 
		if(ptrInnerObj_) delete ptrInnerObj_;
		ptrInnerObj_ = 0;
		myLockerClassTy::initPtr(0);
		destroyed_ = true;
	}
	static InnerTy* GetInnerPtr(void) {return ptrInnerObj_;} 
private:
	singleSubstrataDyn(){}
    singleSubstrataDyn(const singleSubstrataDyn&);
    singleSubstrataDyn operator=(const singleSubstrataDyn&);
	static void Create()
	{
		static meTy theObject;
		AccessPoint_ = &theObject;
	}
	static void OnDeadReference(const wchar_t *msg)
	{
		throw Zeus::libException(ERROR_COD_ZEUSBADHANDLE,ERROR_MSG_ZEUSBADHANDLE,msg);
	}
	virtual ~singleSubstrataDyn()
	{
		myLockerTy locker;
		(void)locker; 
		if(ptrInnerObj_) delete ptrInnerObj_;
		AccessPoint_ = 0;destroyed_ = true;
	}
	static  ptrInnerTy				ptrInnerObj_;
	static  bool					destroyed_;
    static  volatile meTy*			AccessPoint_;   // Singleton access point
};


template
<
typename T,
template <class> class LockProxy
>
volatile singleSubstrataDyn<T,LockProxy>* singleSubstrataDyn<T,LockProxy>::AccessPoint_ = 0;

template
<
typename T,
template <class> class LockProxy
>
bool singleSubstrataDyn<T,LockProxy>::destroyed_ = false;

template
<
typename T,
template <class> class LockProxy
>
T* singleSubstrataDyn<T,LockProxy>::ptrInnerObj_ = 0;

template<typename T>
T* LockProxyNoThreading<T>::obj_= 0;

template<typename T>
T* LockProxyMultiThreading<T>::obj_= 0;

template
<
typename T,
template <class> class LockProxy
>
singleSubstrata<T,LockProxy>* singleSubstrata<T,LockProxy>::AccessPoint_ = 0;

template
<
typename T,
template <class> class LockProxy
>
bool singleSubstrata<T,LockProxy>::destroyed_ = false;

template<typename T>
struct Single_SingletonHolder
{
	typedef singleSubstrata<T,LockProxyNoThreading>	Type;	
};

template<typename T>
struct Multi_SingletonHolder
{
	typedef singleSubstrata<T,LockProxyMultiThreading>	Type;	
};

template<typename T>
struct Single_SingletonDynHolder
{
	typedef singleSubstrataDyn<T,LockProxyNoThreading>	Type;	
};

template<typename T>
struct Multi_SingletonDynHolder
{
	typedef singleSubstrataDyn<T,LockProxyMultiThreading>	Type;	
};



} // end of namespace ZEUS


#endif //SIGLETON_SERIALIZEH
