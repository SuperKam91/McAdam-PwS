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
#ifndef FFTW_STORAGEH
#define FFTW_STORAGEH
//---------------------------------------------------------------------------
#include "ZEUS_StorageBaseManip.h"
#include "FFTW_Traits.h"

namespace fftw
{

////////////////////////////////////////////////////////////////////////////////
// class template FFTWRefCounted
// Implementation of the OwnershipPolicy used by SmartPtr
// Adapts FFTWStorage intrusive reference counting to OwnershipPolicy-specific syntax
////////////////////////////////////////////////////////////////////////////////


long CalcStoSize_Atoms(int dimension,const int *dims,bool xtra,bool isReal);
bool checkDims(int dimension,const int *dims);

template<typename T>
void freeFastMem(T* p)
{if(p) (ReleaseMemType::Instance())->releaseFastMem(reinterpret_cast<typename dataTraits<T>::TypeSz*>(p));}

#if defined(WIN32)
#define _FFTWDATASTOGENTYPENAME		typename
#else
#define _FFTWDATASTOGENTYPENAME 
#endif

template
<
	typename T,
	template <class> class LockPolicy
>
class dataStoGeneral: public Zeus::HandleStorageBase<LockPolicy>
{
public:
    typedef T																dataType;
	typedef typename dataTraits<T>::TypeSz									szType; 	
    typedef T*																ptrType;
    typedef T&																refType;
    typedef const T&														consRefType;
    typedef dataStoGeneral<T,LockPolicy>									myType;
	typedef myType															BaseDynType;
	typedef const myType&													cteRefMyType;
	typedef Zeus::HandleStorageBase<LockPolicy>								BaseType;
	typedef typename BaseType::LockerType									LockerType;
#ifdef __GNUC__
	typedef  typename BaseType::template VolatileType<int*>::Type			ptrVolatDimsType;
	typedef  typename BaseType::template VolatileType<ptrType>::Type		ptrVolatStorageType;
	typedef  typename BaseType::template VolatileType<bool>::Type			VolatBoolType;
#else
	typedef _FFTWDATASTOGENTYPENAME BaseType::VolatileType<int*>::Type		ptrVolatDimsType;
	typedef _FFTWDATASTOGENTYPENAME BaseType::VolatileType<ptrType>::Type	ptrVolatStorageType;
	typedef _FFTWDATASTOGENTYPENAME BaseType::VolatileType<bool>::Type		VolatBoolType;
#endif

	template<
		typename IN_T,
		typename OUT_T,
		template <class> class LockPolicyAux
	>
	friend class Device;


public:
	inline bool				HasMemory(void) const {return (dataSto_ != 0);}
	inline ptrType			getRawPtr(void) const{return dataSto_;}
	inline int				getNDims(void) const {return dimension_;}
	inline int				getDim0(void) const {return *dims_;}

	int						getDims(std::vector<int>& d)
	{
        LockerType locker(Zeus::HandleStorageBase<LockPolicy>::getLocker());
		d.clear();
		for(int i=0;i<dimension_;++i){d.push_back(dims_[i]);}
		return dimension_;
	}
	int						getHighDims(int *dims) const
	{
		if(dimension_ <= 1) return 0;
		for(int i=1;i<dimension_;i++) {dims[i] = dims_[i];}
		return dimension_;
	}

	bool					allocData(int Chg0dim = -1,int ChgXtra=-1)
    {
        if(dataSto_) return true;
        LockerType locker(Zeus::HandleStorageBase<LockPolicy>::getLocker());
        if(dataSto_) return true;
		if(Chg0dim > 0)
		{
			if(Chg0dim < 1)			errGeneral(ERROR_COD_FFTWNDIMISNUL,ERROR_MSG_FFTWNDIMISNUL);
			dims_[0] = Chg0dim;
		}
		if(ChgXtra != -1){xtraSp_ = (ChgXtra!=0);}
        szType  *ptrPtrDta;
		long N_Atoms(CalcStoSize_Atoms(dimension_,dims_,xtraSp_,Loki::IsSameType<szType,T>::value));
		(GetMemType::Instance())->getFastMem(N_Atoms * sizeof(T),&ptrPtrDta);
        new(ptrPtrDta) T[N_Atoms];
        dataSto_ = reinterpret_cast<ptrType>(ptrPtrDta);
        return false;
    }
	long					getStorageSize(bool AtomsSzT,bool incXtra=true) const
	{
		long SzAt(CalcStoSize_Atoms(dimension_,dims_,(incXtra?xtraSp_:false),Loki::IsSameType<szType,T>::value));
		return  SzAt *(AtomsSzT?1:sizeof(T));
	}
    ptrType					detach(void)
	{
		if(Zeus::HandleStorageBase<LockPolicy>::getNRefs() != 1) return 0; 
		LockerType locker(Zeus::HandleStorageBase<LockPolicy>::getLocker());
		if(Zeus::HandleStorageBase<LockPolicy>::getNRefs() != 1) return 0; 
		ptrType ptrTemp(dataSto_);
        dataSto_ = 0;
        return ptrTemp;
    }
    bool					freeMem(void)
	{
		if(Zeus::HandleStorageBase<LockPolicy>::getNRefs() != 1) return false; 
		LockerType locker(Zeus::HandleStorageBase<LockPolicy>::getLocker());
		if(Zeus::HandleStorageBase<LockPolicy>::getNRefs() != 1) return false; 
		ptrType ptrTemp(dataSto_);
        dataSto_ = 0;
		if(ptrTemp){(ReleaseMemType::Instance())->releaseFastMem(reinterpret_cast<szType*>(ptrTemp));}
		return true;
	}

    bool					reattachPtr(void* d,int dim0 = -1)
    {
        if(dataSto_) return false;
		LockerType locker(Zeus::HandleStorageBase<LockPolicy>::getLocker());
        if(dataSto_) return false;
        dataSto_ = reinterpret_cast<ptrType>(d);
		if(dim0 > 0) dims_[0] = dim0;
        return true;
    }
	inline dataType			getLinFast(long i) const
	{return *(dataSto_ + i);}
	inline refType			getAddrLinFast(long i) const
	{return *(dataSto_ + i);}
	inline void				putLinFast(long i,consRefType v) const
	{*(dataSto_ + i) = v;}

	virtual dataType		get(const int *v){errMethodNotExists();return 0;}
	virtual void			put(const int *v,consRefType val){errMethodNotExists();}
	virtual dataType		get(int i){errMethodNotExists();return 0;}
	virtual void			put(int i,consRefType v){errMethodNotExists();}
	virtual dataType		get(int i,int j){errMethodNotExists();return 0;}
	virtual void			put(int i,int j,consRefType v){errMethodNotExists();}
	virtual dataType		get(int i,int j,int k){errMethodNotExists();return 0;}
	virtual void			put(int i,int j,int k,consRefType v){errMethodNotExists();}

    dataStoGeneral(int dimension,bool xtraSp): Zeus::HandleStorageBase<LockPolicy>(),dims_(0),xtraSp_(xtraSp),dimension_(dimension),dataSto_(0)
    {
		if(dimension < 1)			 errGeneral(ERROR_COD_FFTWNDIMISNUL,ERROR_MSG_FFTWNDIMISNUL);
		if(dimension > FFTW_MAXDIMS) errGeneral(ERROR_COD_FFTWNDIMTOBIG,ERROR_MSG_FFTWNDIMTOBIG);
		dims_ = new int[dimension];
	}
    virtual  ~dataStoGeneral(void) throw()
    {
		if(dims_) delete [] dims_;
		if(dataSto_){
			(ReleaseMemType::Instance())->releaseFastMem(reinterpret_cast<szType*>(dataSto_));
		}
    }


protected:
	ptrVolatDimsType dims_;
	void errGeneral(int errNumber,wchar_t* msg) const {throw Zeus::libException(errNumber,msg,*this);}
private:
	inline void errMethodNotExists(void) const {errGeneral(ERROR_COD_FFTWMETHNOTEX,ERROR_MSG_FFTWMETHNOTEX);}
    dataStoGeneral(const dataStoGeneral&);
    dataStoGeneral operator=(const dataStoGeneral&);

	ptrVolatStorageType			dataSto_;
	VolatBoolType				xtraSp_;
	int							dimension_;
};

template
<
typename T,
int dimension,
class checkBoundsPolicy,
template <class> class LockPolicy
>
class dataSto
	: public dataStoGeneral<T,LockPolicy>
	, public checkBoundsPolicy
{
public:
	typedef dataSto<T,dimension,checkBoundsPolicy,LockPolicy>		myType;
	typedef dataStoGeneral<T,LockPolicy>							StorageBaseType;
	
	enum{ stoDim = dimension};
	
	template<typename IN_T,typename OUT_T,template <class> class LockPolicyAux>
	friend class Device;

	static myType* create(int dimDyn,int *dims,bool alloc=false,bool xtra=false){return new myType(dimDyn,dims,alloc,xtra);}

	virtual typename StorageBaseType::dataType	get(const int *v)
	{
		checkBoundsPolicy::check(dataStoGeneral<T,LockPolicy>::getNDims(),v,StorageBaseType::dims_);
		if(!StorageBaseType::HasMemory()) {StorageBaseType::allocData();}		
		return StorageBaseType::getLinFast(GetLinIndex(v));
    }
	virtual void put(const int *v,typename StorageBaseType::consRefType val)
	{
		checkBoundsPolicy::check(dataStoGeneral<T,LockPolicy>::getNDims(),v,StorageBaseType::dims_);
		if(!StorageBaseType::HasMemory()){StorageBaseType::allocData();}		
		StorageBaseType::putLinFast(GetLinIndex(v), val);
    }

	inline typename StorageBaseType::dataType		getFast(const int *v) const
	{
		return StorageBaseType::getLinFast(GetLinIndex(v));
    }
	inline typename StorageBaseType::refType		getAddrFast(const int *v) const
	{
		return getAddrLinFast(GetLinIndex(v));
	}
	inline void	putFast(const int *v,typename StorageBaseType::consRefType val) const
	{
		StorageBaseType::putLinFast(GetLinIndex(v), val);
    }
	virtual				~dataSto(void) throw() {}
private:
    dataSto(int dimDyn,int *dims,bool alloc=false,bool xtra=false) :  dataStoGeneral<T,LockPolicy>(dimDyn,xtra),checkBoundsPolicy()
	{
		LOKI_STATIC_CHECK((dimension == 0),DIMENSIONS_CANT_BE_0);
		for(int i=0;i<dimDyn;i++) StorageBaseType::dims_[i]=dims[i];
		if(!checkDims(dimDyn,StorageBaseType::dims_)) StorageBaseType::errGeneral(ERROR_COD_FFTWWRONGDIMS,ERROR_MSG_FFTWWRONGDIMS);
		if (alloc){StorageBaseType::allocData();}
	}
	long GetLinIndex(const int *v)
	{
		long index(v[0]);
		for(int i=1;i<dataStoGeneral<T,LockPolicy>::getNDims();i++){
			index = index * StorageBaseType::dims_[i] + v[i];
		}
		return index;
	}
#if defined(WIN32)
	Loki::Type2Type<bool> dummy;
#endif //WIN32
	dataSto(const dataSto&);
};

template
<
typename T,
class checkBoundsPolicy,
template <class> class LockPolicy
>
class dataSto<T,1,checkBoundsPolicy,LockPolicy>
	: public dataStoGeneral<T,LockPolicy>
	, public checkBoundsPolicy
{
public:
	typedef dataSto<T,1,checkBoundsPolicy,LockPolicy>			myType;
	typedef dataStoGeneral<T,LockPolicy>						StorageBaseType;

	
	enum{ stoDim = 1};
	
	template<typename IN_T,typename OUT_T,template <class> class LockPolicyAux>
	friend class Device;


	static myType* create(int i,bool alloc=false,bool xtra=false){return new myType(i,alloc,xtra);}

	virtual typename StorageBaseType::dataType	get(int i)
	{
		checkBoundsPolicy::check(i,StorageBaseType::dims_[0]);
		if(!StorageBaseType::HasMemory()){StorageBaseType::allocData();}		
		return StorageBaseType::getLinFast(i);
    }
	virtual void	put(int i, typename StorageBaseType::consRefType v)
	{
		checkBoundsPolicy::check(i,StorageBaseType::dims_[0]);
		if(!StorageBaseType::HasMemory()) {StorageBaseType::allocData();}		
		StorageBaseType::putLinFast(i, v);
    }
	
	inline typename StorageBaseType::dataType	getFast(int i) const
	{
		return StorageBaseType::getLinFast(i);
    }

	inline typename StorageBaseType::refType		getAddrFast(int i) const
	{
		return StorageBaseType::getAddrLinFast(i);
	}
	inline void		putFast(int i,typename StorageBaseType::consRefType v)
	{
		StorageBaseType::putLinFast(i, v);
	}
	virtual				~dataSto(void)  throw() {}

private:
    dataSto(int i,bool alloc=false,bool xtra=false) :  dataStoGeneral<T,LockPolicy>(1,xtra),checkBoundsPolicy()
	{
		StorageBaseType::dims_[0] = i;
		if(!checkDims(1,StorageBaseType::dims_)) StorageBaseType::errGeneral(ERROR_COD_FFTWWRONGDIMS,ERROR_MSG_FFTWWRONGDIMS);
		if (alloc){StorageBaseType::allocData();}
	}

#if defined(WIN32)
	Loki::Type2Type<bool> dummy;
#endif //WIN32
	dataSto(const dataSto&);
};

template
<
typename T,
class checkBoundsPolicy,
template <class> class LockPolicy
>
class dataSto<T,2,checkBoundsPolicy,LockPolicy>
	: public dataStoGeneral<T,LockPolicy>
	, public checkBoundsPolicy
{
public:
	typedef dataSto<T,2,checkBoundsPolicy,LockPolicy>			myType;
	typedef dataStoGeneral<T,LockPolicy>						StorageBaseType;
	
	enum{ stoDim = 2};
	
	template<typename IN_T,typename OUT_T,template <class> class LockPolicyAux>
	friend class Device;


	static myType* create(int i,int j,bool alloc=false,bool xtra=false){return new myType(i,j,alloc,xtra);}

	virtual typename StorageBaseType::dataType		get(int i,int j)
	{
		checkBoundsPolicy::check(i,StorageBaseType::dims_[0],j,StorageBaseType::dims_[1]);
		if(!StorageBaseType::HasMemory()){StorageBaseType::allocData();}		
		return StorageBaseType::getLinFast((i*StorageBaseType::dims_[1])+j);
    }
	virtual void	put(int i,int j,typename StorageBaseType::consRefType v)
	{
		checkBoundsPolicy::check(i,StorageBaseType::dims_[0],j,StorageBaseType::dims_[1]);
		if(!StorageBaseType::HasMemory()){StorageBaseType::allocData();}		
		StorageBaseType::putLinFast((i*StorageBaseType::dims_[1])+j,v);
    }

	inline typename StorageBaseType::dataType	getFast(int i,int j) const
	{
		return StorageBaseType::getLinFast((i*StorageBaseType::dims_[1])+j);
    }
	inline typename StorageBaseType::refType	getAddrFast(int i,int j) const
	{
		return getAddrLinFast((i*StorageBaseType::dims_[1])+j);
	}
	inline void		putFast(int i,int j,typename StorageBaseType::consRefType v) const
	{
		StorageBaseType::putLinFast((i*StorageBaseType::dims_[1])+j,v);
    }
	virtual				~dataSto(void)  throw()  {}


private:
    dataSto(int i,int j,bool alloc=false,bool xtra=false) :  dataStoGeneral<T,LockPolicy>(2,xtra),checkBoundsPolicy()
	{
		StorageBaseType::dims_[0] = i;StorageBaseType::dims_[1] = j;
		if(!checkDims(2,StorageBaseType::dims_)) StorageBaseType::errGeneral(ERROR_COD_FFTWWRONGDIMS,ERROR_MSG_FFTWWRONGDIMS);
		if (alloc){StorageBaseType::allocData();}
	}

#if defined(WIN32)
	Loki::Type2Type<bool> dummy;
#endif //WIN32
	dataSto(const dataSto&);
};

template
<
typename T,
class checkBoundsPolicy,
template <class> class LockPolicy
>
class dataSto<T,3,checkBoundsPolicy,LockPolicy>
	: public dataStoGeneral<T,LockPolicy>
	, public checkBoundsPolicy
{
public:
	typedef dataSto<T,3,checkBoundsPolicy,LockPolicy>			myType;
	typedef dataStoGeneral<T,LockPolicy>						StorageBaseType;
	
	enum{ stoDim = 3};
	template<typename IN_T,typename OUT_T,template <class> class LockPolicyAux>
	friend class Device;


	static myType* create(int i,int j,int k,bool alloc=false,bool xtra=false){return new myType(i,j,k,alloc,xtra);}

	virtual typename StorageBaseType::dataType	get(int i,int j,int k)
	{
		checkBoundsPolicy::check(i,StorageBaseType::dims_[0],j,StorageBaseType::dims_[1],k,StorageBaseType::dims_[2]);
		if(!StorageBaseType::HasMemory()){StorageBaseType::allocData();}		
		return StorageBaseType::getLinFast((i*StorageBaseType::dims_[1]+ j)* StorageBaseType::dims_[2]+k);
    }
	virtual void	put(int i,int j,int k,typename StorageBaseType::consRefType v)
	{
		checkBoundsPolicy::check(i,StorageBaseType::dims_[0],j,StorageBaseType::dims_[1],k,StorageBaseType::dims_[2]);
		if(!StorageBaseType::HasMemory()){StorageBaseType::allocData();}		
		StorageBaseType::putLinFast((i*StorageBaseType::dims_[1]+j)*StorageBaseType::dims_[2]+k,v);
    }
	inline typename StorageBaseType::dataType getFast(int i,int j,int k) const
	{
		return StorageBaseType::getLinFast((i*StorageBaseType::dims_[1]+j)*StorageBaseType::dims_[2]+k);
    }
	inline typename StorageBaseType::refType getAddrFast(int i,int j,int k) const
	{
		return getAddrLinFast((i*StorageBaseType::dims_[1]+j)*StorageBaseType::dims_[2]+k);
	}
	inline void putFast(int i,int j,int k,typename StorageBaseType::consRefType v) const
	{
		StorageBaseType::putLinFast((i*StorageBaseType::dims_[1]+j)*StorageBaseType::dims_[2]+k,v);
    }
	virtual				~dataSto(void)  throw() {}

private:

    dataSto(int i,int j,int k, bool alloc=false,bool xtra=false) :  dataStoGeneral<T,LockPolicy>(3,xtra),checkBoundsPolicy()
	{
		StorageBaseType::dims_[0] = i;StorageBaseType::dims_[1] = j;StorageBaseType::dims_[2] = k;
		if(!checkDims(3,StorageBaseType::dims_)) StorageBaseType::errGeneral(ERROR_COD_FFTWWRONGDIMS,ERROR_MSG_FFTWWRONGDIMS);
		if (alloc){StorageBaseType::allocData();}
	}
#if defined(WIN32)
	Loki::Type2Type<bool> dummy;
#endif //WIN32
	dataSto(const dataSto&);
};


template
<
typename T,
typename checkBoundsPolicy = DEFAULT_STORAGE_CHECK_BOUNDS,
template <class> class LockPolicy = DEFAULT_STORAGE_THREAD_POLICY
>
struct DataStorageFactory
{
	static dataStoGeneral<T,LockPolicy> *createStorage(int Ndims,int *dims,bool alloc=false,bool xtra=false)
	{
		if((Ndims < 1) || (Ndims > FFTW_MAXDIMS))
			throw Zeus::libException(ERROR_COD_FFTWNDIMISNUL,ERROR_MSG_FFTWNDIMISNUL,
				L"dataStoGeneral<T,LockPolicy> *createStorage(int Ndims,int *dims,bool alloc=false,bool xtra=false)");

		switch(Ndims){
			case 1:
				return dataSto<T,1,checkBoundsPolicy,LockPolicy>::create(dims[0],alloc,xtra);
			case 2:
				return dataSto<T,2,checkBoundsPolicy,LockPolicy>::create(dims[0],dims[1],alloc,xtra);
			case 3:
				return dataSto<T,3,checkBoundsPolicy,LockPolicy>::create(dims[0],dims[1],dims[2],alloc,xtra);
			default:
				return dataSto<T,0,checkBoundsPolicy,LockPolicy>::create(Ndims,dims,alloc,xtra);
		};
	}
};



} // end of namespace fftw





#endif
