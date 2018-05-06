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
#ifndef ZEUS_STORAGEBASICMANIPH
#define ZEUS_STORAGEBASICMANIPH

#include <string>

#ifdef WIN32
#include <tchar.h>
#endif //WIN32

#include <complex>
#include <math.h>
#include <typeinfo>
#include <exception>
#include <locale>
#include <map>
#include <functional>
#include <map>
#include "Threads.h"
#include "SmartPtr.h"
#include "TypeManip.h"
#include "static_check.h"
#include "ZEUS_LokiVariant.h"
#include "ZEUS_Strings.h"
#include "ZEUS_Exceptions.h"




//---------------------------------------------------------------------------

namespace Zeus
{

typedef Variant<LOKI_TYPELIST_4(Loki::NullType,double,int,std::wstring)> sDBField;

struct ps_VariableType
{
	sDBField		DefValue_;
	std::wstring	VarId_;

	ps_VariableType(void):VarId_(),DefValue_()
	{}
	ps_VariableType(const std::wstring& varid):VarId_(varid),DefValue_()
	{}
	ps_VariableType(const std::wstring& varid,const sDBField& DefValue):VarId_(varid),DefValue_(DefValue)
	{}

};

typedef Variant<LOKI_TYPELIST_5(Loki::NullType,double,int,std::wstring,ps_VariableType)> DBField;

DBField sDBField2DBField(const sDBField& sDB);

template
<
class T,
class H = Loki::NullType
>
struct CollectionDiskStorageHelper
{
	typedef CollectionDiskStorageHelper<T,H>		meType;
	typedef H										HeaderType;
	typedef T										StorageType;
	HeaderType										Header_;
	StorageType										Storage_;

	inline meType& swap(meType& rhs)
	{
		std::swap(Header_,rhs.Header_);
		Storage_.swap(rhs.Storage_);
		return *this;
	}
};

template <class P>
class ZEUSRefCounted
{
public:
    ZEUSRefCounted()
    {}
    
    template <class U>
    ZEUSRefCounted(const ZEUSRefCounted<U>&)
    {}
    
    static P Clone(const P& val)
    {
		if(!val) return 0;
        val->AddRef();
        return val;
    }
    
    static bool Release(const P& val)
    {
		if(!val) return true;
		if(val->Release())
			return false;
		else return true;
	}
    
    enum { destructiveCopy = false };
    
    void Swap(ZEUSRefCounted&)
    {}
};

template<typename T>
struct ObjHandle
{
typedef typename Loki::SmartPtrDef<T,ZEUSRefCounted,Loki::DisallowConversion,Loki::NoCheck,Loki::DefaultSPStorage>::type Type;
};


template<typename T>
typename T::baseIntVolatile Kill_me(Loki::SmartPtr<T,ZEUSRefCounted,Loki::DisallowConversion,Loki::NoCheck,Loki::DefaultSPStorage>& sPt)
{
	typedef typename T::baseIntVolatile myTempIntTy;
	if(GetImpl(sPt) == 0) return 0;
	myTempIntTy NRs(sPt->getNRefs());
	Reset(sPt,0);
	return --NRs;
}

template<typename T>
inline bool Check_me(const Loki::SmartPtr<T,ZEUSRefCounted,Loki::DisallowConversion,Loki::NoCheck,Loki::DefaultSPStorage>& sPt)
{return(GetImpl(sPt) != 0);}


template<typename T>
inline typename Loki::DefaultSPStorage<T>::PointerType GetInnerPtr(const Loki::SmartPtr<T,ZEUSRefCounted,Loki::DisallowConversion,Loki::NoCheck,Loki::DefaultSPStorage>& sPt)
{return GetImpl(sPt);}

template
<
	template <class> class LockPolicy = Loki::SingleThreaded
>
class HandleStorageBase
{
public:
    typedef HandleStorageBase<LockPolicy>										myType;
	typedef LockPolicy<myType>													LockClassType;
	typedef const LockClassType&												LockClassCteRefType;
	typedef typename LockClassType::Lock										LockerType;
	
	template<typename T>
	struct VolatileType
	{
		typedef typename LockPolicy<T>::VolatileType Type;
	};

	typedef typename VolatileType<typename LockClassType::IntType>::Type		baseIntVolatile;

	template<typename FT1>
	friend class ZEUSRefCounted;
	
	template<template <class> class CheckP>
	friend baseIntVolatile Kill_me(Loki::SmartPtr<HandleStorageBase,ZEUSRefCounted,Loki::DisallowConversion,CheckP,Loki::DefaultSPStorage>& );

	inline LockClassCteRefType getLocker() const
	{return LockerClass_;}

	inline baseIntVolatile getNRefs(void) const
	{return NRefs_;}

	// TODO ; Should I put a virtual desctructor in here ?
protected:
	HandleStorageBase() : LockerClass_(),NRefs_(1)
	{}
	HandleStorageBase(const HandleStorageBase&) : LockerClass_(),NRefs_(1)
	{}
	HandleStorageBase& operator=(const HandleStorageBase&)
	{return *this;}

private:
    inline baseIntVolatile Release(void)
	{return LockerClass_.AtomicDecrement(NRefs_);}
	inline void AddRef(void)
	{LockerClass_.AtomicIncrement(NRefs_);}

	LockClassType		LockerClass_;
	baseIntVolatile		NRefs_;
};


template
<
	template <class> class LockPolicy = Loki::SingleThreaded
>
class StaticLocker
{
public:
    typedef StaticLocker<LockPolicy>											myType;
	typedef LockPolicy<myType>													LockClassType;
	typedef const LockClassType&												LCStatCtRefType;
	typedef typename LockClassType::Lock										LockerType;
	
	template<typename T>
	struct VolatileType
	{typedef typename LockPolicy<T>::VolatileType Type;};

	inline LCStatCtRefType getLocker() const
	{return LockerClass_;}
protected:
	StaticLocker() : LockerClass_()
	{}
	~StaticLocker(void)
	{}; 
private:
	LockClassType	LockerClass_;
};


class DeviceBaseType : public HandleStorageBase<Loki::SingleThreaded>
{
public:
	virtual void	Initialise(void) = 0;

	inline std::wstring DevName(void)
	{return Id_;}

	DeviceBaseType(const std::wstring& Id) : HandleStorageBase<Loki::SingleThreaded>(),Id_(Id)
	{}
#if defined(WIN32)
	virtual ~DeviceBaseType(void) = 0
	{}
#else
	// THIS IS NOT COMPLIANT !!!!
	virtual ~DeviceBaseType(void)
	{}
#endif
private:
	const std::wstring	Id_;
};


} // end of namespace Zeus


#endif //ZEUS_STORAGEBASICMANIPH 
