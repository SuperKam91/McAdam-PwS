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
#ifndef ZEUS_WORKSPACEH
#define ZEUS_WORKSPACEH

#include <vector>
#include <complex>
#include <algorithm>
#include "fftw3.h"
#include "TypeManip.h"
#include "ZEUS_Exceptions.h"
#include "ZEUS_StorageBaseManip.h"
#include "ZEUS_PhysicsMath.h"

#define LARRAYNOCOPYSEMANTICS	1

#ifndef WIN32
#undef LARRAYNOCOPYSEMANTICS
#endif

#define INTERPOL1VERYSMALL		1e-8
#define PWSAVSAMPLES			200.0


namespace Zeus
{

namespace Private
{
	template<typename T>
	struct LA_AtomTraits
	{typedef	T				AtomType;};
	template<typename T>
	struct LA_AtomTraits<std::complex<T> >
	{typedef	T				AtomType;};
} //Private
//
template<typename T>
class LargeArrayBase
{
protected:
	unsigned long		MaxSz_;
	T					*dta_;

	inline static void err(void)
	{throw libException(ERROR_COD_ZEUSWRKSPALLOC, ERROR_MSG_ZEUSWRKSPALLOC, L"Zeus::LargeArrayBase");}


	template<typename U>
	inline static U *AllocFastMemory(unsigned long NBytes,Loki::Type2Type<U>)
	{return new U[NBytes / sizeof(U)];}
	inline static double *AllocFastMemory(unsigned long NBytes,Loki::Type2Type<double>)
	{return reinterpret_cast<double*>(fftw_malloc(NBytes));}
	inline static float *AllocFastMemory(unsigned long NBytes,Loki::Type2Type<float>)
	{return reinterpret_cast<float*>(fftwf_malloc(NBytes));}
	inline static long double *AllocFastMemory(unsigned long NBytes,Loki::Type2Type<long double>)
	{return reinterpret_cast<long double*>(fftwl_malloc(NBytes));}

	template<typename U>
	inline static void ReleaseFastMamory(U* ptrDta)
	{delete [] ptrDta;}

	inline static void ReleaseFastMamory(double* ptrDta)
	{if(ptrDta) fftw_free(ptrDta);}
	inline static void ReleaseFastMamory(float* ptrDta)
	{if(ptrDta) fftwf_free(ptrDta);}
	inline static void ReleaseFastMamory(long double* ptrDta)
	{if(ptrDta) fftwl_free(ptrDta);}

	inline static void FreeMem(T *ptrDta)
	{ReleaseFastMamory(reinterpret_cast<typename Private::LA_AtomTraits<T>::AtomType *>(ptrDta));}

	void CopyFrom(const LargeArrayBase& rhs)
	{
		if (!rhs.dta_ || (MaxSz_ != rhs.MaxSz_))
		{
			if(dta_) {FreeMem(dta_);dta_ = 0;}
		}

		if (!rhs.dta_) {dta_ = 0;}
		else {
			T		*tmp_new, *tmp_newPiv;
			const T	*tmp_me(rhs.dta_), * const tmp_meEnd(tmp_me + rhs.MaxSz_);

			if(!dta_ || (MaxSz_ != rhs.MaxSz_)) tmp_new = reinterpret_cast<T*>(AllocFastMemory(rhs.MaxSz_ * sizeof(T),Loki::Type2Type<typename Private::LA_AtomTraits<T>::AtomType>()));
			else tmp_new = dta_;
			for(tmp_newPiv = tmp_new;tmp_me != tmp_meEnd;++tmp_newPiv,++tmp_me)
				*tmp_newPiv = *tmp_me;
			dta_ = tmp_new;
		}
		MaxSz_		= rhs.MaxSz_;
	}
	inline void do_reset(const T& defVal)
	{
		T	*tmp(dta_),* const tmpEnd(dta_ + MaxSz_);
		for(;tmp!=tmpEnd;++tmp)
			*tmp = defVal;
	}

	inline LargeArrayBase(const LargeArrayBase& rhs)
		:MaxSz_(0),dta_(0)
	{CopyFrom(rhs);}

	inline void swap(LargeArrayBase& rhs)
	{std::swap(MaxSz_,rhs.MaxSz_);std::swap(dta_,rhs.dta_);}

	inline void Make(unsigned long MaxSz,T* dta)
	{
		if((MaxSz != MaxSz_) || dta)
		{if(dta_) FreeMem(dta_);dta_=dta;}
		MaxSz_=MaxSz;
	}
	inline void Make(unsigned long MaxSz,const T& defVal,T* dta)
	{
		Make(MaxSz,dta);
		begin();
		do_reset(defVal);
	}

	inline explicit LargeArrayBase(unsigned long MaxSz):MaxSz_(MaxSz),dta_(0)
	{}
	inline explicit LargeArrayBase(unsigned long MaxSz,const T& defVal):MaxSz_(MaxSz),dta_(0)
	{
		if(MaxSz_)
		{
			dta_ = reinterpret_cast<T*>(AllocFastMemory(MaxSz_ * sizeof(T),Loki::Type2Type<typename Private::LA_AtomTraits<T>::AtomType>()));
			do_reset(defVal);
		}
	}

	inline ~LargeArrayBase(void)
	{if(dta_) {FreeMem(dta_);dta_ = 0;}}

public:
	typedef LargeArrayBase<T>	Type;
	typedef T					AtomType;
	typedef T*					iterator;
	typedef const T*			const_iterator;
	typedef unsigned long		SizeType;

	inline unsigned long getSz(void) const
	{return MaxSz_;}

	inline bool IsInitialized(void) const
	{return dta_!= 0;}

	inline T* begin(void)
	{
		if(MaxSz_)
		{
			if(dta_) return dta_;
			else return (dta_ = reinterpret_cast<T*>(AllocFastMemory(MaxSz_ * sizeof(T),Loki::Type2Type<typename Private::LA_AtomTraits<T>::AtomType>())));
		}
		else err();
		return 0;
	}

	inline const T* begin(void) const
	{
		if(!dta_){err();return 0;}
		return dta_;
	}

	inline T* end(void)
	{return begin() + MaxSz_;}

	inline const T* end(void) const
	{return begin() + MaxSz_;}

	inline void releaseDta(void)
	{if(dta_) {FreeMem(dta_);dta_ = 0;}}
	inline T* reset(const T& val = T())
	{
		if(!MaxSz_) return 0;
		begin();
		do_reset(val);
		return dta_;
	}

};

//
template<typename T>
class LArr1D: public LargeArrayBase<T>
{
#ifdef LARRAYNOCOPYSEMANTICS
protected:
#else
public:
#endif //LARRAYNOCOPYSEMANTICS

	inline LArr1D(const LArr1D& rhs)
		:LargeArrayBase<T>(rhs)
	{}
	inline LArr1D& operator=(const LArr1D& rhs)	
	{
		if(this == &rhs) return *this;
		LargeArrayBase<T>::CopyFrom(rhs);
		return *this;
	}

public:
	typedef struct{}		HeaderType; //To make it readable

	inline explicit LArr1D(unsigned long MaxSz=0)
		:LargeArrayBase<T>(MaxSz)
	{}
	inline explicit LArr1D(unsigned long MaxSz,const T& defVal)
		:LargeArrayBase<T>(MaxSz,defVal)
	{}
	inline explicit LArr1D(const std::vector<T>& rhs)
		:LargeArrayBase<T>((unsigned long)rhs.size())
	{
		if(rhs.empty()) return;
		T	*tDestPiv(LargeArrayBase<T>::begin());
		typename std::vector<T>::const_iterator piv(rhs.begin());
		typename std::vector<T>::const_iterator const end(rhs.end());
		for(;piv != end;++piv,++tDestPiv)
		{*tDestPiv = *piv;}
	}

	inline LArr1D& swap(LArr1D& rhs)
	{
		LargeArrayBase<T>::swap(rhs);
		return *this;
	}
	inline LArr1D& Assign(const LArr1D& rhs)
	{return *this = rhs;}

	inline LArr1D& Make(unsigned long MaxSz,T* dta=0)
	{
		LargeArrayBase<T>::Make(MaxSz,dta);
		return *this;
	}
	inline LArr1D& Make(unsigned long MaxSz,const T& defVal,T* dta=0)
	{
		LargeArrayBase<T>::Make(MaxSz,defVal,dta);
		return *this;
	}

	inline LArr1D Clone(void) const
	{return LArr1D(*this);}

	inline T&	operator[] (unsigned long index)
	{return (LargeArrayBase<T>::begin())[index];}

	inline const T&	operator[] (unsigned long index) const
	{return (LargeArrayBase<T>::begin())[index];}
};

//
template<typename T>
class LArr2D: public LargeArrayBase<T>
{
	unsigned long ptrMetric_;
#ifdef LARRAYNOCOPYSEMANTICS
protected:
#else
public:
#endif //LARRAYNOCOPYSEMANTICS
	inline LArr2D(const LArr2D& rhs):LargeArrayBase<T>(rhs),ptrMetric_(rhs.ptrMetric_)
	{}
	inline LArr2D& operator=(const LArr2D& rhs)	
	{
		if(this == &rhs) return *this;
		LargeArrayBase<T>::CopyFrom(rhs);
		ptrMetric_ = rhs.ptrMetric_;
		return *this;
	}
public:
	typedef unsigned long MetricType;
	typedef  struct{}	HeaderType; //To make it readable

	inline explicit LArr2D(unsigned long MaxSz=0,unsigned long ptrMetric=0):LargeArrayBase<T>(MaxSz),ptrMetric_(ptrMetric)
	{}
	inline explicit LArr2D(unsigned long MaxSz,unsigned long ptrMetric,const T& defVal):LargeArrayBase<T>(MaxSz,defVal),ptrMetric_(ptrMetric)
	{}

	inline LArr2D& swap(LArr2D& rhs)
	{
		LargeArrayBase<T>::swap(rhs);
		std::swap(ptrMetric_,rhs.ptrMetric_);
		return *this;
	}
	inline LArr2D& Assign(const LArr2D& rhs)
	{return *this = rhs;}

	inline LArr2D& Make(unsigned long MaxSz,unsigned long ptrMetric,T* dta=0)
	{
		LargeArrayBase<T>::Make(MaxSz,dta);
		ptrMetric_ = ptrMetric;
		return *this;
	}
	inline LArr2D& Make(unsigned long MaxSz,unsigned long ptrMetric,const T& defVal,T* dta=0)
	{
		LargeArrayBase<T>::Make(MaxSz,defVal,dta);
		ptrMetric_ = ptrMetric;
		return *this;
	}

	inline LArr2D Clone(void) const
	{return LArr2D(*this);}
	inline T&	operator() (unsigned long line,unsigned long column)
	{return (LargeArrayBase<T>::begin())[line * ptrMetric_ + column];}

	inline const T&	operator() (unsigned long line,unsigned long column) const
	{return (LargeArrayBase<T>::begin())[line * ptrMetric_ + column];}

	inline T*	ptr2(unsigned long line,unsigned long column)
	{return (LargeArrayBase<T>::begin()) + (line * ptrMetric_ + column);}

	inline const T*	ptr2(unsigned long line,unsigned long column) const
	{return (LargeArrayBase<T>::begin()) + (line * ptrMetric_ + column);}

	inline unsigned long getPtrMetric(void) const
	{return ptrMetric_;}
protected:
	inline unsigned long SetPtrMetric(unsigned long NewPtrMetric)
	{
		unsigned long temp(ptrMetric_);
		ptrMetric_ = NewPtrMetric;
		return temp;
	}
};

//
template<typename T>
class HD_LArr1D : public HandleStorageBase<Loki::SingleThreaded>, public LArr1D<T>
{
#ifdef LARRAYNOCOPYSEMANTICS
protected:
#else
public:
#endif //LARRAYNOCOPYSEMANTICS

	inline HD_LArr1D(const HD_LArr1D& rhs)
		:HandleStorageBase<Loki::SingleThreaded>(),LArr1D<T>(rhs)
	{}

	inline HD_LArr1D& operator=(const HD_LArr1D& rhs)	
	{
		if(this == &rhs) return *this;
		CopyFrom(rhs);
		return *this;
	}

public:
	typedef LArr1D<T>		InnerDataType;

	inline explicit HD_LArr1D(unsigned long MaxSz=0)
		:LArr1D<T>(MaxSz)
	{}

	inline explicit HD_LArr1D(unsigned long MaxSz,const T& defVal)
		:LArr1D<T>(MaxSz,defVal)
	{}

	inline explicit HD_LArr1D(const InnerDataType& rhs)
		:HandleStorageBase<Loki::SingleThreaded>(),LArr1D<T>(rhs)
	{}

	inline explicit HD_LArr1D(const std::vector<T>& rhs)
		:HandleStorageBase<Loki::SingleThreaded>(),LArr1D<T>(rhs)
	{}

	inline HD_LArr1D& swap(HD_LArr1D& rhs)
	{
		LArr1D<T>::swap(rhs);
		return *this;
	}

	inline HD_LArr1D& Assign(const HD_LArr1D& rhs)
	{return *this = rhs;}
	inline HD_LArr1D& Make(unsigned long MaxSz,T* dta=0)
	{LargeArrayBase<T>::Make(MaxSz,dta);return *this;}
	inline HD_LArr1D& Make(unsigned long MaxSz,const T& defVal,T* dta=0)
	{LargeArrayBase<T>::Make(MaxSz,defVal,dta);return *this;}
	inline HD_LArr1D* Clone(void) const
	{return new HD_LArr1D(*this);}
};

//
template<typename T>
class HD_LArr2D : public HandleStorageBase<Loki::SingleThreaded>, public LArr2D<T>
{
#ifdef LARRAYNOCOPYSEMANTICS
protected:
#else
public:
#endif //LARRAYNOCOPYSEMANTICS

	inline HD_LArr2D(const HD_LArr2D& rhs)
		:HandleStorageBase<Loki::SingleThreaded>(),LArr2D<T>(rhs)
	{}

	inline HD_LArr2D& operator=(const HD_LArr2D& rhs)	
	{
		if(this == &rhs) return *this;
		LargeArrayBase<T>::CopyFrom(rhs);
		SetPtrMetric(rhs.getPtrMetric());
		return *this;
	}

public:
	typedef LArr2D<T>		InnerDataType;

	inline explicit HD_LArr2D(unsigned long MaxSz=0,unsigned long ptrMetric=0)
		:HandleStorageBase<Loki::SingleThreaded>(),LArr2D<T>(MaxSz,ptrMetric)
	{}
	inline explicit HD_LArr2D(unsigned long MaxSz,unsigned long ptrMetric,const T& defVal)
		:HandleStorageBase<Loki::SingleThreaded>(),LArr2D<T>(MaxSz,ptrMetric,defVal)
	{}

	inline explicit HD_LArr2D(const InnerDataType& rhs)
		:HandleStorageBase<Loki::SingleThreaded>(),LArr2D<T>(rhs)
	{}

	inline HD_LArr2D& swap(HD_LArr2D& rhs)
	{LArr2D<T>::swap(rhs);return *this;}

	inline HD_LArr2D& Assign(const HD_LArr2D& rhs)
	{return *this = rhs;}
	inline HD_LArr2D& Make(unsigned long MaxSz,unsigned long ptrMetric,T* dta=0)
	{LArr2D<T>::Make(MaxSz,ptrMetric,dta);return *this;}
	inline HD_LArr2D& Make(unsigned long MaxSz,unsigned long ptrMetric,const T& defVal,T* dta=0)
	{LArr2D<T>::Make(MaxSz,ptrMetric,defVal,dta);return *this;}
	inline HD_LArr2D* Clone(void) const
	{return new HD_LArr2D(*this);}

};

//
enum	UseBoundsType{UB_DEFAULT = -1,UB_USE=1,UB_NOUSE=0};

//
template<typename T>
class DataPlaneSurface;
//
struct PlaneBoundsType
{
	static	UseBoundsType		DefaultUseBounds_;
	static	float				DefaultBoundFactor_;

	enum {INVALID_BOUND= -1};
	float		BoundFactor_;
	int			YBound_;
	int			XBound_;

	inline static int SetDefaultUseBounds(int flag)
	{
		int temp(DefaultUseBounds_);
		DefaultUseBounds_ = static_cast<UseBoundsType>(flag==1);
		return temp;
	}

	inline static float SetDefaultBoundingFactor(float bFactor)
	{
		float temp(DefaultBoundFactor_);
		DefaultBoundFactor_ = bFactor;
		return temp;
	}

	inline static float GetDefaultBoundingFactor(void)
	{return DefaultBoundFactor_;}

	inline static UseBoundsType GetDefaultUseBounds(void)
	{return DefaultUseBounds_;}


	inline PlaneBoundsType(void)
		:YBound_(INVALID_BOUND),XBound_(INVALID_BOUND),BoundFactor_(DefaultBoundFactor_)
	{}

	inline PlaneBoundsType(int YBound,int XBound)
		:YBound_(YBound),XBound_(XBound),BoundFactor_(DefaultBoundFactor_)
	{}
	inline PlaneBoundsType(int YBound,int XBound,float	BoundFactor)
		:YBound_(YBound),XBound_(XBound),BoundFactor_(BoundFactor)
	{}
	inline float GetBoundFactor(void) const
	{return BoundFactor_;}
	inline void Reset(void)
	{YBound_ = XBound_ = INVALID_BOUND;BoundFactor_ = DefaultBoundFactor_;}

	inline void ResetBoundsCoords(void)
	{YBound_ = XBound_ = INVALID_BOUND;}
	inline	void SetUseBounds(UseBoundsType& UseBounds) const
	{
		if(UseBounds==UB_DEFAULT)
		{UseBounds = DefaultUseBounds_;}
		if(UseBounds && (YBound_ > 0) && (XBound_ > 0))
		{UseBounds = UB_USE;}
		else{UseBounds = UB_NOUSE;}
	}
	inline int	HasBounds(void) const
	{
		if(DefaultUseBounds_ && (YBound_ > 0) && (XBound_ > 0)) return UB_USE;
		return UB_NOUSE;
	}

	inline void SymmBounds(void)
	{
		if(YBound_ == XBound_)
			return;
		if((XBound_ == INVALID_BOUND) || (YBound_ == INVALID_BOUND))
		{
			YBound_ = XBound_ = INVALID_BOUND;
			return;
		}
		int largerB((YBound_ > XBound_)?YBound_:XBound_);
		YBound_ = XBound_ = largerB;
	}
};

// ----------------------
enum	SumSelectType{SEL_SUM=0x01,SEL_SUMSQ=0x02,SEL_SUMSUMSQ=0x03};
// ----------------------
template<typename T>
class DataPlaneSurface
{
public:
	typedef T													AtomType;
	typedef HD_LArr2D<T>										DataStorageType;
	typedef typename ObjHandle<DataStorageType>::Type			DataHandleType;
	typedef typename HD_LArr2D<T>::InnerDataType				DataInnerType;
protected:
	PlaneBoundsType		Bounds_;
	int					YSz_;
	int					XSz_;

	inline static void	errDataPlaneSurface(int errCode,wchar_t* msg,wchar_t* method)
	{
		throw Zeus::libException(errCode,msg,std::wstring(L"DataPlaneSurface::") + std::wstring(method));
	}

public:
	inline DataPlaneSurface(void)
		:YSz_(-1),XSz_(-1),Bounds_()
	{}
	inline DataPlaneSurface(int YSz,int XSz,bool reset)
		:YSz_(YSz),XSz_(XSz),Bounds_()
	{
		Data_ = DataHandleType(new DataStorageType(((YSz_ < 0) || (XSz_ < 0))?0:YSz_*XSz_,(XSz_ < 0)?0:XSz_));
		if(reset) Data_->reset();
	}

	inline	AtomType* Reset(const AtomType& val = AtomType())
	{
		Bounds_.Reset();
		if(Check_me(Data_))
		{return Data_->reset(val);}
		return 0;
	}

	inline	int  Release(void)
	{
		Bounds_.Reset();YSz_ = XSz_ = -1;
		int nr(Kill_me(Data_));
		return nr;
	}

	inline  bool IsInit(void) const
	{
		if(!Check_me(Data_)) return false;
		return Data_->IsInitialized();
	}

	inline	PlaneBoundsType GetBounds(void) const
	{return Bounds_;}

	inline	void			GetBounds(int& Yb,int& Xb) const
	{
		Xb = Bounds_.XBound_; Yb = Bounds_.YBound_;
	}
	inline void GetSz(int& YSz,int& XSz) const
	{YSz = YSz_;XSz = XSz_;}

	inline	PlaneBoundsType SetBounds(const PlaneBoundsType& NewBounds)
	{
		PlaneBoundsType OldBounds(Bounds_);
		Bounds_ = NewBounds;
		return OldBounds;
	}

	int	HasBounds(void) const
	{return Bounds_.HasBounds();}
	inline float GetBoundFactor(void) const
	{return Bounds_.BoundFactor_;}
	inline	PlaneBoundsType SetBounds(int YBound,int XBound)
	{
		PlaneBoundsType OldBounds(Bounds_);
		Bounds_.XBound_ = XBound;
		Bounds_.YBound_ = YBound;
		return OldBounds;
	}

	inline	PlaneBoundsType ResetBounds(void)
	{
		PlaneBoundsType OldBounds(Bounds_);
		Bounds_.Reset();
		return OldBounds;
	}

	inline	PlaneBoundsType ResetBoundsCoords(void)
	{
		PlaneBoundsType OldBounds(Bounds_);
		Bounds_.ResetBoundsCoords();
		return OldBounds;
	}

	inline  PlaneBoundsType SymmetriseBounds(void)
	{
		PlaneBoundsType OldBounds(Bounds_);
		Bounds_.SymmBounds();
		return OldBounds;
	}

	inline	void MakeNewSz(int YSz,int XSz,bool reset=false)
	{
		Kill_me(Data_);
		Data_ = DataHandleType(new DataStorageType(YSz*XSz,XSz));
		YSz_ = YSz; XSz_ = XSz;
		Bounds_.Reset();
		if(reset) Data_->reset();;
	}
	inline	void Swap(DataPlaneSurface& NewData)
	{
		std::swap(Bounds_,NewData.Bounds_);
		std::swap(YSz_,NewData.YSz_);
		std::swap(XSz_,NewData.XSz_);
		Data_.Swap(NewData.Data_);
	}
	inline	DataInnerType&	GetInnerData(void)
	{
		if(!Check_me(Data_))
		{
			Data_ = DataHandleType(new DataStorageType(0,0));
		}
		return *Data_;
	}
	inline	const DataInnerType&	GetInnerData(void) const
	{return *Data_;}

	inline  DataPlaneSurface	Clone(void) const
	{
		DataPlaneSurface temp(*this);
		if(!Check_me(Data_))
			return temp;
		temp.Data_ = DataHandleType((temp.Data_->Clone()));
		return temp;
	}

	template<typename T0,typename T1>	
	friend	PlaneBoundsType GetMinMaxBounds(const DataPlaneSurface<T0>& fst,const DataPlaneSurface<T1>& sec ,bool MinMax);
private:
	DataHandleType		Data_;
};

// ----------------------
template<typename T>
class VectTriMatxBase
{
public:
	typedef T													AtomType;
	typedef HD_LArr1D<T>										DataStorageType;
	typedef typename ObjHandle<DataStorageType>::Type			DataHandleType;
	typedef typename HD_LArr1D<T>::InnerDataType				DataInnerType;

	inline explicit VectTriMatxBase(int Sz=0)
	{
		if(Sz > 0)
		{Data_ = DataHandleType(new DataStorageType(Sz));}
	}
	inline explicit VectTriMatxBase(DataStorageType* rhs)
		:Data_(DataHandleType(rhs))
	{}
	inline explicit VectTriMatxBase(const DataInnerType& rhs)
		:Data_(DataHandleType(new DataStorageType(rhs)))
	{}
	inline	AtomType*  Reset(void)
	{
		if(Check_me(Data_))
		{return Data_->reset();}
		return 0;
	}
	inline	AtomType*  Release(void)
	{
		int nr(Kill_me(Data_));
		return nr;
	}
	inline	void MakeNewSz(int Sz,bool reset=false)
	{
		Kill_me(Data_);
		Data_ = DataHandleType(new DataStorageType(Sz));
		if(reset) Data_->reset();;
	}
	inline	void Swap(VectTriMatxBase& NewData)
	{Data_.Swap(NewData.Data_);}
	inline	int	 GetSz(void) const
	{
		if(!Check_me(Data_)) return 0;
		return Data_->getSz();
	}
	inline  bool Empty(void) const
	{
		if(!Check_me(Data_)) return false;
		return Data_->IsInitialized();
	}
	inline	const DataInnerType&	GetInnerData(void) const
	{return *Data_;}
	inline	DataInnerType&	GetInnerData(void)
	{
		if(!Check_me(Data_))
		{Data_ = DataHandleType(new DataStorageType(0));}
		return *Data_;
	}
	inline  VectTriMatxBase	Clone(void) const
	{
		if(!Check_me(Data_)) return VectTriMatxBase();
		return VectTriMatxBase(Data_->Clone());
	}
private:
	DataHandleType		Data_;
};
// ----------------------
template<typename T>
class AlgVect: public VectTriMatxBase<T>
{
private:
	inline AlgVect(int end,const VectTriMatxBase<T>& rhs)
		:VectTriMatxBase<T>(rhs),End_(end)
	{}
public:
	inline AlgVect(void)
		:VectTriMatxBase<T>(0),End_(0)
	{}
	inline explicit AlgVect(int sz,bool reset = false)
		:VectTriMatxBase<T>(sz),End_(sz)
	{if(reset) VectTriMatxBase<T>::Reset();}
	inline  void SetEnd(int end)
	{	int t(End_);
		End_ = end;
		return t;
	}
	inline  int GetEnd(void) const
	{return End_;}
	inline	typename VectTriMatxBase<T>::AtomType * Reset(void)
	{
		End_ = VectTriMatxBase<T>::GetSz();
		return VectTriMatxBase<T>::Reset();
	}
	inline	int  Release(void)
	{
		End_ = 0;
		return VectTriMatxBase<T>::Release();
	}
	inline	void MakeNewSz(int Sz,bool reset=false)
	{
		VectTriMatxBase<T>::MakeNewSz(Sz,reset);
		End_ = Sz;
	}
	inline	void Swap(AlgVect& NewData)
	{
		VectTriMatxBase<T>::Swap(NewData);
		std::swap(End_,NewData.End_);
	}
	inline	AlgVect Clone(void) const
	{return AlgVect(End_,VectTriMatxBase<T>::Clone());}

	inline T	operator*(const AlgVect& rhs) const
	{
		if(!End_ || !rhs.End_) return T(0.0);
		typename VectTriMatxBase<T>::DataInnerType::const_iterator			piv(VectTriMatxBase<T>::GetInnerData().begin());
		typename VectTriMatxBase<T>::DataInnerType::const_iterator			pivEnd(piv + ((End_ >= rhs.End_)?rhs.End_:End_));
		typename VectTriMatxBase<T>::DataInnerType::const_iterator			pivRhs(rhs.GetInnerData().begin());
		T																	result(0.0);
//TODO OpenMp
		for(;piv != pivEnd;++piv,++pivRhs)
		{
			result += Zeus::Mult(*piv,*pivRhs);
		}
		return result;
	}
	inline AlgVect	operator+(const AlgVect& rhs) const
	{
		if(!End_ && !rhs.End_) return AlgVect();
		if(!End_) return rhs.Clone();
		if(!rhs.End_) return Clone();
		AlgVect t((End_ >= rhs.End_)?Clone():rhs.Clone());

		typename VectTriMatxBase<T>::DataInnerType::iterator			piv(t.GetInnerData().begin());
		typename VectTriMatxBase<T>::DataInnerType::const_iterator	const pivEnd(piv + ((End_ >= rhs.End_)?rhs.End_:End_));
		typename VectTriMatxBase<T>::DataInnerType::const_iterator	pivOrg((End_ >= rhs.End_)?rhs.GetInnerData().begin():VectTriMatxBase<T>::GetInnerData().begin());
//TODO OpenMp
		for(;piv != pivEnd;++piv,++pivOrg)
		{*piv += *pivOrg;}
		return t;
	}

	inline AlgVect	operator*(const T & rhs) const
	{
		if(!(rhs.End_)) return AlgVect();
		AlgVect t(Clone());
		typename VectTriMatxBase<T>::DataInnerType::iterator				piv(t.GetInnerData().begin());
		typename VectTriMatxBase<T>::DataInnerType::const_iterator		const pivEnd(piv + End_);

//TODO OpenMp
		for(;piv!=pivEnd;++piv)
		{*piv *= rhs;}
		return t;
	}

private:
	int End_;
};

// ----------------------
template<typename T>
struct PowerFunctorVect
{
	T MultCte_;
	inline T	operator()(const AlgVect<std::complex<T> >& x)
	{
		return ((x * x).real()) * MultCte_;
	}
	inline PowerFunctorVect(T MultCte)
		:MultCte_(MultCte)
	{}
};


//-----------------------
template<typename T>
class AlgTriMatrix: public VectTriMatxBase<T>
{
private:
	inline AlgTriMatrix(int metric,const VectTriMatxBase<T>& rhs)
		:VectTriMatxBase<T>(rhs),metric_(metric)
	{}
public:
	inline AlgTriMatrix(void)
		:VectTriMatxBase<T>(0),metric_(0)
	{}
//
	inline explicit AlgTriMatrix(int metric,bool reset = false)
		:VectTriMatxBase<T>((metric * (metric + 1)) >> 1),metric_(metric)
	{if(reset) VectTriMatxBase<T>::Reset();}
//
	inline	typename VectTriMatxBase<T>::AtomType * Reset(void)
	{
		return VectTriMatxBase<T>::Reset();
	}
//
	inline	int  Release(void)
	{
		return VectTriMatxBase<T>::Release();
	}
//
	inline	void MakeNewSz(int metric,bool reset=false)
	{
		VectTriMatxBase<T>::MakeNewSz((metric * (metric + 1)) >> 1,reset);
		metric_ = metric;
	}
//
	inline	void Swap(AlgTriMatrix& NewData)
	{
		VectTriMatxBase<T>::Swap(NewData);
		std::swap(metric_,NewData.metric_);
	}
//
	inline	AlgTriMatrix Clone(void) const
	{return AlgTriMatrix(metric_,VectTriMatxBase<T>::Clone());}
//
	inline  int  GetMetric(void) const
	{return metric_;}
//
//
	template<typename T1>
	inline AlgVect<T1>	operator*(const AlgVect<T1>& rhs) const
	{
		if(!rhs.GetEnd()) return AlgVect<T1>();
		typename VectTriMatxBase<T>::DataInnerType::const_iterator pivMe(VectTriMatxBase<T>::GetInnerData().begin());

		int jumpSz;
		AlgVect<T1> t(metric_);
		t.Reset();
		typename VectTriMatxBase<T>::DataInnerType::const_iterator	pivMeAux;
		typename AlgVect<T1>::DataInnerType::iterator				piv(t.GetInnerData().begin());
		typename AlgVect<T1>::DataInnerType::const_iterator const	pivEnd(piv + t.GetInnerData().getSz());
		typename AlgVect<T1>::DataInnerType::const_iterator			pivRhs(rhs.GetInnerData().begin());
		typename AlgVect<T1>::DataInnerType::const_iterator			pivRhsAux(pivRhs);
		typename AlgVect<T1>::DataInnerType::const_iterator const	pivEndRhs(pivRhs + rhs.GetEnd());

		for(;piv != pivEnd;++piv,++pivMe,pivRhs = rhs.GetInnerData().begin())
		{
			jumpSz = metric_ - 1;
			if(pivRhsAux != pivEndRhs) ++pivRhsAux;
			for(pivMeAux = pivMe;pivRhs != pivRhsAux;++pivRhs,pivMeAux += (jumpSz--))
			{*piv += (*pivRhs * (*pivMeAux));}
		}
		return t;
	}
//
	template<typename T1>
	inline AlgTriMatrix	operator*(const T1 & rhs) const
	{
		AlgTriMatrix t(Clone());
		typename VectTriMatxBase<T>::DataInnerType::iterator			piv(t.GetInnerData().begin());
		typename VectTriMatxBase<T>::DataInnerType::iterator			const pivEnd(piv + t.GetSz());

		for(;piv != pivEnd;++piv)
		{*piv *= rhs;}
		return t;
	}
//
	template<typename T1>
	inline AlgTriMatrix&	ArrayMultInPlace(const T1 & rhs)
	{
		typename VectTriMatxBase<T>::DataInnerType::iterator				pivDest(VectTriMatxBase<T>::GetInnerData().begin());
		typename VectTriMatxBase<T>::DataInnerType::const_iterator			const endDest(pivDest + VectTriMatxBase<T>::GetSz());

		for(;pivDest != endDest;++pivDest)
		{*pivDest *= rhs;}
		return *this;
	}
//
	inline AlgTriMatrix& operator+=(const AlgTriMatrix& rhs)
	{
		typename VectTriMatxBase<T>::DataInnerType::iterator				pivDest;
		typename VectTriMatxBase<T>::DataInnerType::const_iterator			pivSrc(rhs.GetInnerData().begin());
		typename VectTriMatxBase<T>::DataInnerType::const_iterator			const endSrc(pivSrc + VectTriMatxBase<T>::GetSz());

		pivDest = VectTriMatxBase<T>::GetInnerData().begin();
		for(;pivSrc != endSrc;++pivSrc,++pivDest)
		{*pivDest += *pivSrc;}
		
		return *this;
	}
//
	inline AlgTriMatrix& ArrayMultInPlace(const AlgTriMatrix& rhs)
	{
		typename VectTriMatxBase<T>::DataInnerType::iterator				pivDest;
		typename VectTriMatxBase<T>::DataInnerType::const_iterator			pivSrc(rhs.GetInnerData().begin());
		typename VectTriMatxBase<T>::DataInnerType::const_iterator			const endSrc(pivSrc + VectTriMatxBase<T>::GetSz());

		pivDest = VectTriMatxBase<T>::GetInnerData().begin();
		for(;pivSrc != endSrc;++pivSrc,++pivDest)
		{*pivDest *= *pivSrc;}
		return *this;
	}
//
	inline AlgTriMatrix& ArrayDivInPlace(const AlgTriMatrix& rhs)
	{
		typename VectTriMatxBase<T>::DataInnerType::iterator				pivDest;
		typename VectTriMatxBase<T>::DataInnerType::const_iterator			pivSrc(rhs.GetInnerData().begin());
		typename VectTriMatxBase<T>::DataInnerType::const_iterator			const endSrc(pivSrc + VectTriMatxBase<T>::GetSz());

		pivDest = VectTriMatxBase<T>::GetInnerData().begin();
		for(;pivSrc != endSrc;++pivSrc,++pivDest)
		{*pivDest /= *pivSrc;}
		return *this;
	}
//
private:
	int	 metric_;
};
//-----------------------
template<typename T>
struct	AvHolderType
{
	int						Hits_;
	T						Data_;
	AvHolderType(void)
		:Hits_(0)
	{}
};
// ----------------------
template<typename T>
class FourierPlane : public DataPlaneSurface<T>
{
private:
	inline static T   Get1ValByInterPol(T FirstValue,T SecondVal)
	{
		double absFValue(std::abs(FirstValue));
		if(absFValue < INTERPOL1VERYSMALL)
			return FirstValue;

		return FirstValue * (std::abs(SecondVal)/(2.0 * absFValue));
	}

	inline static int AzimuthIndex(int y,int x)
	{
		return toInt(std::sqrt(static_cast<double>(y*y)+static_cast<double>(x*x)) + 0.5);
	}
/*
	inline static int BinWidth(const int FMode)
	{
		const double	rad(static_cast<double>(FMode));
		const int		r(toInt((std::sqrt(rad*rad + (2.0 * PWSAVSAMPLES / PI)) - rad) + 0.5));
		return (r > 0) ? r : 1;
	}
*/

	inline static int BinWidth(const int FMode, const double AvSamples)
	{
		const double	rad(static_cast<double>(FMode));
		const double	radLimit(std::sqrt(AvSamples / PITIMES2));
		const double	radLimitX2(2.0*radLimit);
		double			r;

		if(rad < radLimit)
		{r = radLimitX2;}
		else
		{r = AvSamples / (PI * rad);}
		return ((r > 1.0) ? toInt(r+0.5) : 1);
	}

	inline FourierPlane(const DataPlaneSurface<T>& rhs)
		:DataPlaneSurface<T>(rhs)
	{}
public:
	typedef	typename Private::LA_AtomTraits<T>::AtomType	FourierPlaneStorageType;	

	inline FourierPlane(void)
		:DataPlaneSurface<T>()
	{}
	inline explicit FourierPlane(int Sz,bool reset=false)
		:DataPlaneSurface<T>((Sz << 1),Sz + 1,reset)
	{		
		if(Sz % 2)
			DataPlaneSurface<T>::errDataPlaneSurface(ERROR_COD_DTSURFWRONGSZ,ERROR_MSG_DTSURFWRONGSZ,L"FourierPlane");
	}

	inline FourierPlane(int Sz,int YBound,int XBound,bool reset=true)
		:DataPlaneSurface<T>((Sz << 1),Sz + 1,reset)
	{
		if(Sz % 2)
			DataPlaneSurface<T>::errDataPlaneSurface(ERROR_COD_DTSURFWRONGSZ,ERROR_MSG_DTSURFWRONGSZ,L"FourierPlane");
		DataPlaneSurface<T>::SetBounds(YBound,XBound);
	}

	inline FourierPlane(int Sz,const PlaneBoundsType& Bounds,bool reset=true)
		:DataPlaneSurface<T>((Sz << 1),Sz + 1,reset)
	{
		if(Sz % 2)
			DataPlaneSurface<T>::errDataPlaneSurface(ERROR_COD_DTSURFWRONGSZ,ERROR_MSG_DTSURFWRONGSZ,L"FourierPlane");
		DataPlaneSurface<T>::SetBounds(Bounds);
	}
	inline void	MakeNewSz(int Sz,bool reset=false)
	{
		DataPlaneSurface<T>::MakeNewSz((Sz << 1),Sz + 1,reset);
	}
	void	IntegrateSurface(SumSelectType select,typename DataPlaneSurface<T>::AtomType& Sum,typename DataPlaneSurface<T>::AtomType& SumSq,UseBoundsType UseBounds=UB_DEFAULT) const
	{
		typename DataPlaneSurface<T>::AtomType						value_2X(0.0),value_1X(0.0);
		typename DataPlaneSurface<T>::AtomType						valueSq_2X(0.0),valueSq_1X(0.0);
		const typename DataPlaneSurface<T>::DataInnerType&			InnerData(DataPlaneSurface<T>::GetInnerData());
		typename DataPlaneSurface<T>::DataInnerType::const_iterator	piv(InnerData.begin());
		typename DataPlaneSurface<T>::DataInnerType::MetricType		const ArrMetric(InnerData.getPtrMetric());
		typename DataPlaneSurface<T>::DataInnerType::const_iterator	pivAux;
		int			YBound,XBound,i,j;

		DataPlaneSurface<T>::Bounds_.SetUseBounds(UseBounds);

		if(UseBounds)
		{YBound = DataPlaneSurface<T>::Bounds_.YBound_;XBound = DataPlaneSurface<T>::Bounds_.XBound_;}
		else
		{YBound = ((DataPlaneSurface<T>::YSz_ >> 1) + 1);XBound = DataPlaneSurface<T>::XSz_;}

		typename DataPlaneSurface<T>::DataStorageType::const_iterator	pivBack(piv + (ArrMetric * DataPlaneSurface<T>::YSz_));
		typename DataPlaneSurface<T>::DataStorageType::MetricType		const ArrMetricAux(ArrMetric - 1);
		int																const HalfPlane(DataPlaneSurface<T>::YSz_ >> 1);

// TODO OpenMP
		for(j=0;j<YBound;++j,piv += ArrMetric,pivBack -= ArrMetric)
		{
// TODO OpenMP
			for(i=0,pivAux = piv;i<XBound;++i,++pivAux)
			{
				if(!(i) || (i==ArrMetricAux))
				{
					if(select & SEL_SUM){value_1X += *pivAux;}
					if(select & SEL_SUMSQ){valueSq_1X += Zeus::Sq(*pivAux);}
				}
				else
				{
					if(select & SEL_SUM){value_2X += *pivAux;}
					if(select & SEL_SUMSQ){valueSq_2X += Zeus::Sq(*pivAux);}
				}
			}
			if((j != 0) && (j != HalfPlane))
			{
// TODO OpenMP
				for(i=0,pivAux = pivBack;i<XBound;++i,++pivAux)
				{
					if(!(i) || (i==ArrMetricAux))
					{
						if(select & SEL_SUM){value_1X += *pivAux;}
						if(select & SEL_SUMSQ){valueSq_1X += Zeus::Sq(*pivAux);}
					}
					else
					{
						if(select & SEL_SUM){value_2X += *pivAux;}
						if(select & SEL_SUMSQ){valueSq_2X += Zeus::Sq(*pivAux);}
					}
				}
			}
		}
		if(select & SEL_SUM)	Sum		= ((value_2X   * 2.0) + value_1X) / static_cast<double>(DataPlaneSurface<T>::YSz_ * DataPlaneSurface<T>::YSz_);
		if(select & SEL_SUMSQ)	SumSq	= ((valueSq_2X * 2.0) + valueSq_1X) / static_cast<double>(DataPlaneSurface<T>::YSz_ * DataPlaneSurface<T>::YSz_);
		return;
	}

	inline FourierPlane Clone(void) const
	{return FourierPlane(DataPlaneSurface<T>::Clone());}
	inline	PlaneBoundsType	GetMaxBounds(void) const
	{
		return PlaneBoundsType((DataPlaneSurface<T>::YSz_>> 1) + 1,DataPlaneSurface<T>::XSz_,DataPlaneSurface<T>::Bounds_.BoundFactor_);
	}
	FourierPlane& AverageMode1(void)
	{
		typename DataPlaneSurface<T>::DataInnerType&						InnerData(DataPlaneSurface<T>::GetInnerData());
		typename DataPlaneSurface<T>::DataStorageType::MetricType			const ArrMetric(InnerData.getPtrMetric());
		typename DataPlaneSurface<T>::DataStorageType::iterator				const Org(InnerData.begin());
		double	 Values[5],Minor,Major;
		int		 MajorIndex(0),MinorIndex(0);


		// IMPORTANT TODO NEED TO FIX THIS

		// Zero mode 0
		*Org  = 0.0;

		// Zero nyquist frequency
		*(Org + ((ArrMetric * ArrMetric)-1)) = 0.0;

		Minor = Major = Values[0] = std::abs(*(Org + 1));
		Values[1] = std::abs(*(Org + ArrMetric));
		Values[2] = std::abs(*(Org + ArrMetric + 1));
		Values[3] = std::abs(*(Org + (DataPlaneSurface<T>::YSz_ - 1) * ArrMetric));
		Values[4] = std::abs(*(Org + ((DataPlaneSurface<T>::YSz_ - 1) * ArrMetric) + 1));
		for(int i=1;i<5;++i)
		{
			if(Values[i] > Major)
			{
				Major = Values[i];MajorIndex = i;
			}
			else
			{
				if(Values[i] < Minor)
				{
					Minor = Values[i];MinorIndex = i;
				}
			}
		}

		double Average(0.0);
/*
		Values[MinorIndex] = Values[MajorIndex] = -1.0;
		for(int i=0;i<5;++i)
		{
			if(Values[i] >= 0.0)
			{
				Average += Values[i];
			}
		}
		
		Average /= ((MinorIndex==MajorIndex)?4:3);
*/
		Average = Minor;

		*(Org + 1)								*= (Average / std::abs(*(Org + 1)));
		*(Org + ArrMetric)						*= (Average / std::abs(*(Org + ArrMetric)));
		*(Org + ArrMetric + 1)					*= (Average / std::abs(*(Org + ArrMetric + 1)));
		*(Org + (DataPlaneSurface<T>::YSz_ - 1) * ArrMetric)			*= (Average / std::abs(*(Org + (DataPlaneSurface<T>::YSz_ - 1) * ArrMetric)));
		*(Org + ((DataPlaneSurface<T>::YSz_ - 1) * ArrMetric) + 1)	*= (Average / std::abs(*(Org + ((DataPlaneSurface<T>::YSz_ - 1) * ArrMetric) + 1)));

		return *this;
	}	

	FourierPlane& AzAv(double AV_Samples,int CalibPowerSp,UseBoundsType UseBounds=UB_DEFAULT)
	{
		typename DataPlaneSurface<T>::DataInnerType&						InnerData(DataPlaneSurface<T>::GetInnerData());
		typename DataPlaneSurface<T>::DataStorageType::MetricType			const ArrMetric(InnerData.getPtrMetric());
		LArr1D<AvHolderType<T> >	AvHolder(toInt((SQRT2 * static_cast<double>(ArrMetric))) + 30);
		unsigned long AvHolderSz(AvHolder.getSz());
		LArr1D<AvHolderType<T> >	AvHolderTemp(AvHolderSz);
		typename DataPlaneSurface<T>::DataStorageType::const_iterator		const Org(InnerData.begin());
		typename DataPlaneSurface<T>::DataStorageType::iterator				const OrgDest(InnerData.begin());
		typename DataPlaneSurface<T>::DataStorageType::const_iterator		piv,pivAux;
		typename DataPlaneSurface<T>::DataStorageType::iterator				pivDest,pivDestAux;
		int		 YBound,YBoundInf,XBound,i,j,k1;

		DataPlaneSurface<T>::Bounds_.SetUseBounds(UseBounds);

		if(UseBounds)
		{
			YBound = DataPlaneSurface<T>::Bounds_.YBound_;
			XBound = DataPlaneSurface<T>::Bounds_.XBound_;}
		else
		{
			YBound = ((DataPlaneSurface<T>::YSz_ >> 1) + 1);
			XBound = DataPlaneSurface<T>::XSz_;
		}
		
		if(YBound == ((DataPlaneSurface<T>::YSz_ >> 1) + 1))
		{YBoundInf = -YBound + 2;}
		else
		{YBoundInf = -YBound + 1;}

// initialize original power
		T	OrgPower = (*Org).Clone();
		OrgPower.Reset();
// computes original total power
		for(int j = YBoundInf;j!=YBound;++j)
		{
			if(j<0) {piv = (Org + (ArrMetric * (DataPlaneSurface<T>::YSz_ + j)));}
			else {piv = (Org + (ArrMetric * j));}

			for(int i = 0;i!=XBound;++i)
			{
				if((j==0) && (i==0))
					continue;
				pivAux = piv + i;
				OrgPower += (((!i) || (i == (XBound-1)))?(*pivAux):(*pivAux) * 2.0);
			}
		}

// computes the average cross-power spectrum
		for(j = YBoundInf;j!=YBound;++j)
		{
			if(j<0)
			{piv = (Org + (ArrMetric * (DataPlaneSurface<T>::YSz_ + j)));}
			else
			{piv = (Org + (ArrMetric * j));}

			for(i = 0;i!=XBound;++i)
			{
				pivAux = piv + i;
				k1 = AzimuthIndex(j,i);
				if(!(AvHolder[k1].Hits_))
				{
					AvHolder[k1].Data_ =  (((!i) || (i == (XBound-1)))?(*pivAux).Clone():(*pivAux) * 2.0);
					AvHolder[k1].Hits_ =  (((!i) || (i == (XBound-1)))?1:2);
				}
				else
				{
					AvHolder[k1].Data_ += (((!i) || (i == (XBound-1)))?(*pivAux):(*pivAux) * 2.0);					
					AvHolder[k1].Hits_ += (((!i) || (i == (XBound-1)))?1:2);
				}
			}
		}

		int deltaBin,BinInf,BinSup,SupExtra;
		AvHolderType<T> temp;
#ifdef REMOVELOWFMODES  
		AvHolder[0].Hits_ = 1;
		AvHolder[0].Data_.Reset();
		AvHolder[1].Hits_ = 1;
		AvHolder[1].Data_.Reset();
#endif
		for(int k=0;k<AvHolderSz;++k)
		{
			AvHolderTemp[k].Hits_ = 1;
			// Force initialisation
			AvHolderTemp[k].Data_ = (*Org).Clone();
			AvHolderTemp[k].Data_.Reset();
		}

		for(i=1;i<AvHolderSz;++i)
		{
			deltaBin = BinWidth(i,AV_Samples);
			BinInf = i - (deltaBin >> 1);
			if(BinInf < 1)
			{
				SupExtra = 1 - BinInf;
				BinInf = 1;
			}
			else
			{SupExtra = 0;}

			BinSup = i + ((deltaBin >> 1) + (deltaBin % 2)) + SupExtra;
			if(BinSup > AvHolderSz) BinSup = static_cast<int>(AvHolderSz);

			for(j=BinInf;j<BinSup;++j)
			{
				if(j==BinInf)
				{
					temp.Hits_	= AvHolder[j].Hits_;
					temp.Data_	= (AvHolder[j].Data_).Clone();
				}
				else
				{
					temp.Hits_	+= AvHolder[j].Hits_;
					temp.Data_	+= AvHolder[j].Data_;
				}
			}

			if(BinInf != BinSup)
			{
				AvHolderTemp[i].Hits_ = temp.Hits_;
				AvHolderTemp[i].Data_ = (temp.Data_).Clone();
			}
		}
// Put averaged power spectrum back in place
		for(j = YBoundInf;j!=YBound;++j)
		{
			if(j<0)
			{pivDest = (OrgDest + (ArrMetric * (DataPlaneSurface<T>::YSz_ + j)));}
			else
			{pivDest = (OrgDest + (ArrMetric * j));}

			for(i = 0;i!=XBound;++i)
			{
				pivDestAux = pivDest + i;
				k1 = AzimuthIndex(j,i);
				*pivDestAux = (AvHolderTemp[k1].Data_ * (1.0 / static_cast<double>(AvHolderTemp[k1].Hits_)));
			}
		}

// initialize final power
		T	FinalPower = (*Org).Clone();
		FinalPower.Reset();
// computes final total power
		for(int j = YBoundInf;j!=YBound;++j)
		{
			if(j<0) {piv = (Org + (ArrMetric * (DataPlaneSurface<T>::YSz_ + j)));}
			else {piv = (Org + (ArrMetric * j));}

			for(int i = 0;i!=XBound;++i)
			{
				if((j==0) && (i==0))
					continue;
				pivAux = piv + i;
				FinalPower += (((!i) || (i == (XBound-1)))?(*pivAux):(*pivAux) * 2.0);
			}
		}

		OrgPower.ArrayDivInPlace(FinalPower);

		double* tPiv(OrgPower.GetInnerData().begin());
		// only uses 100,143,217,353
		double	PowerCorr(*tPiv);
		int		index(OrgPower.GetMetric());
		for(int i=1;i<4;++i,--index)
		{tPiv+=index;PowerCorr+= *(tPiv);}
		PowerCorr /= 4.0;

		if(CalibPowerSp)
		{
			for(int j = YBoundInf;j!=YBound;++j)
			{
				pivDestAux = OrgDest + (ArrMetric * ((j<0)?DataPlaneSurface<T>::YSz_ + j:j));

				for(int i = 0;i!=XBound;++i)
				{
					if((j==0) && (i==0))
						continue;
					(*(pivDestAux + i)).ArrayMultInPlace(PowerCorr);
				}
			}
		}

		return *this;
	}

	FourierPlane&	ShiftInPlace(FourierPlaneStorageType yShift,FourierPlaneStorageType xShift,UseBoundsType UseBounds=UB_DEFAULT)
	{
		typename DataPlaneSurface<T>::DataInnerType&						InnerData(DataPlaneSurface<T>::GetInnerData());
		typename DataPlaneSurface<T>::DataStorageType::iterator				piv(InnerData.begin());
		typename DataPlaneSurface<T>::DataStorageType::MetricType			const ArrMetric(InnerData.getPtrMetric());
		int																	YBound,XBound,i,j;
		typename DataPlaneSurface<T>::DataStorageType::iterator				pivAux;
		FourierPlaneStorageType				phi;
		const FourierPlaneStorageType		phiYFactor((-PITIMES2 * yShift) / static_cast<FourierPlaneStorageType>(DataPlaneSurface<T>::YSz_));
		const FourierPlaneStorageType		phiXFactor((-PITIMES2 * xShift) / static_cast<FourierPlaneStorageType>(DataPlaneSurface<T>::YSz_));

		DataPlaneSurface<T>::Bounds_.SetUseBounds(UseBounds);


		if(UseBounds)
		{YBound = DataPlaneSurface<T>::Bounds_.YBound_;XBound = DataPlaneSurface<T>::Bounds_.XBound_;}
		else
		{YBound = XBound = ((DataPlaneSurface<T>::YSz_ >> 1) + 1);}

		for(j=0;j<YBound;++j,piv += ArrMetric)
		{
			const FourierPlaneStorageType yPhi(phiYFactor * static_cast<FourierPlaneStorageType>(j));	
// TODO OpenMP
			for(i=0,pivAux = piv;i<XBound;++i,++pivAux)
			{
				phi = yPhi + (phiXFactor * static_cast<FourierPlaneStorageType>(i));
				*pivAux *= std::complex<FourierPlaneStorageType>(std::cos(phi),std::sin(phi)) ;
			}
		}

		if(YBound == ((DataPlaneSurface<T>::YSz_ >> 1) + 1))
		{--YBound;}

		piv = (InnerData.begin()) + (ArrMetric * (DataPlaneSurface<T>::YSz_ - 1));
// TODO OpenMP
		for(j=1;j<YBound;++j,piv -= ArrMetric)
		{
			const FourierPlaneStorageType yPhi(phiYFactor * static_cast<FourierPlaneStorageType>(-j));	
// TODO OpenMP
			for(i=0,pivAux = piv;i<XBound;++i,++pivAux)
			{
				phi = yPhi + (phiXFactor * static_cast<FourierPlaneStorageType>(i));
				*pivAux *= std::complex<FourierPlaneStorageType>(std::cos(phi),std::sin(phi)) ;
			}
		}
		return *this;
	}
//
	template<typename T1>
	FourierPlane& Transform(T1 & functor,UseBoundsType UseBounds=UB_DEFAULT)
	{
		typename DataPlaneSurface<T>::DataInnerType&						InnerData(DataPlaneSurface<T>::GetInnerData());
		typename DataPlaneSurface<T>::DataStorageType::iterator				piv(InnerData.begin());
		typename DataPlaneSurface<T>::DataStorageType::MetricType			const ArrMetric(InnerData.getPtrMetric());
		int																	YBound,XBound,i,j;
		typename DataPlaneSurface<T>::DataStorageType::iterator				pivAux;

		DataPlaneSurface<T>::Bounds_.SetUseBounds(UseBounds);
	
		if(UseBounds)
		{
			YBound = DataPlaneSurface<T>::Bounds_.YBound_;
			XBound = DataPlaneSurface<T>::Bounds_.XBound_;
			if((YBound == ((DataPlaneSurface<T>::YSz_ >> 1) + 1)) && (XBound == DataPlaneSurface<T>::XSz_))
			{UseBounds = UB_NOUSE;}
		}

		if(!UseBounds)
		{
			typename DataPlaneSurface<T>::DataStorageType::const_iterator	const pivEnd(piv + (ArrMetric * DataPlaneSurface<T>::YSz_));
			typename DataPlaneSurface<T>::DataStorageType::const_iterator	const pivOrg(piv);
// TODO OpenMP
			for(;piv != pivEnd; ++piv)
			{
				functor(*piv,static_cast<int>(piv - pivOrg));
			}
			return *this;
		}
		
		typename DataPlaneSurface<T>::DataStorageType::const_iterator	const pivOrg(piv);

// TODO OpenMP
		for(j=0;j<YBound;++j,piv += ArrMetric)
		{
// TODO OpenMP
			for(i=0,pivAux = piv;i<XBound;++i,++pivAux)
			{
				functor(*pivAux,static_cast<int>(pivAux - pivOrg));
			}
		}

		if(YBound == ((DataPlaneSurface<T>::YSz_ >> 1) + 1))
		{--YBound;}
		piv = (InnerData.begin()) + (ArrMetric * (DataPlaneSurface<T>::YSz_ - 1));
// TODO OpenMP
		for(j=1;j<YBound;++j,piv -= ArrMetric)
		{
// TODO OpenMP
			for(i=0,pivAux = piv;i<XBound;++i,++pivAux)
			{
				functor(*pivAux,static_cast<int>(pivAux - pivOrg));
			}
		}
		return *this;
	}

	FourierPlane& CleanUp(const PlaneBoundsType& Limits)
	{
		if((DataPlaneSurface<T>::Bounds_.DefaultUseBounds_ == UB_NOUSE) ||
		   ((DataPlaneSurface<T>::Bounds_.YBound_ >= Limits.YBound_) &&
		   (DataPlaneSurface<T>::Bounds_.XBound_ >= Limits.XBound_)))
		{return *this;}

		typename DataPlaneSurface<T>::DataInnerType&						InnerData(DataPlaneSurface<T>::GetInnerData());
		typename DataPlaneSurface<T>::DataInnerType::iterator				const pivOrg(InnerData.begin());
		typename DataPlaneSurface<T>::DataInnerType::MetricType				const ArrMetric(InnerData.getPtrMetric());
		typename DataPlaneSurface<T>::DataInnerType::iterator				piv;
		int		YBound,XBound,YOff(1),XStart;

		YBound = DataPlaneSurface<T>::Bounds_.YBound_;
		XBound = DataPlaneSurface<T>::Bounds_.XBound_;

		if(Limits.YBound_ == ((DataPlaneSurface<T>::YSz_ >> 1) + 1))
		{++YOff;}

// TODO OpenMP
		for(int j = -(Limits.YBound_)+YOff;j != Limits.YBound_;++j)
		{
// TODO OpenMP
			if(j<0)
			{piv = (pivOrg + (ArrMetric * (DataPlaneSurface<T>::YSz_ + j)));}
			else
			{piv = (pivOrg + (ArrMetric * j));}
			if(std::abs(j) >= YBound)
			{
				XStart	= 0;
			}
			else
			{
				XStart	= XBound;
			}
// TODO OpenMP
			for(int i = XStart;i !=Limits.XBound_;++i)
			{
				if(i<0)
				{*(piv + (DataPlaneSurface<T>::XSz_ + i)) = 0.0;}
				else
				{*(piv + i) = 0.0;}
			}
		}		
		return *this;	
	}

	template<typename T3,typename T2,typename T1>
	friend FourierPlane<T3> Multiply(Loki::Type2Type<T3>,const FourierPlane<T1>& fst,const FourierPlane<T2>& sec,UseBoundsType UseBounds);
	template<typename T3,typename T2>
	friend FourierPlane<T3> MultiplyInPlace(FourierPlane<T3>& fst,const FourierPlane<T2>& sec,UseBoundsType UseBounds);
	template<typename T1,typename T2>
	friend	FourierPlane<T1> Add(const FourierPlane<T1>& fst,const FourierPlane<T2>& sec,UseBoundsType UseBounds);
	template<typename T1,typename T2>
	friend	FourierPlane<T1> AddInPlace(FourierPlane<T1>& fst,const FourierPlane<T2>& sec,UseBoundsType UseBounds);
	template<typename T2,typename T0,typename T1>
	friend FourierPlane<T2> Transform(Loki::Type2Type<T2>,const FourierPlane<T0>& sec,T1& functor,UseBoundsType UseBounds,int BoundInfo);
};

//-------------------------------------------------------------------------------------------------------------------
template<typename T>
class RealPlane : public DataPlaneSurface<T>
{
public:
	inline RealPlane(void)
		:DataPlaneSurface<T>()
	{}
	inline RealPlane(int Sz,bool reset=false)
		:DataPlaneSurface<T>(Sz,Sz,reset)
	{		
		if(Sz % 2)
			DataPlaneSurface<T>::errDataPlaneSurface(ERROR_COD_DTSURFWRONGSZ,ERROR_MSG_DTSURFWRONGSZ,L"RealPlane");
	}

	inline RealPlane(int Sz,int YBound,int XBound,bool reset=true)
		:DataPlaneSurface<T>(Sz,Sz,reset)
	{
		if(Sz % 2)
			DataPlaneSurface<T>::errDataPlaneSurface(ERROR_COD_DTSURFWRONGSZ,ERROR_MSG_DTSURFWRONGSZ,L"RealPlane");
		DataPlaneSurface<T>::SetBounds(YBound,XBound);
	}

	inline RealPlane(int Sz,const PlaneBoundsType& Bounds,bool reset=true)
		:DataPlaneSurface<T>(Sz,Sz,reset)
	{
		if(Sz % 2)
			DataPlaneSurface<T>::errDataPlaneSurface(ERROR_COD_DTSURFWRONGSZ,ERROR_MSG_DTSURFWRONGSZ,L"RealPlane");
		DataPlaneSurface<T>::SetBounds(Bounds);
	}
	inline void	MakeNewSz(int Sz,bool reset=false)
	{
		DataPlaneSurface<T>::MakeNewSz(Sz,Sz,reset);
	}

	void	IntegrateSurface(SumSelectType select,typename DataPlaneSurface<T>::AtomType& Sum,typename DataPlaneSurface<T>::AtomType& SumSq,UseBoundsType UseBounds=UB_DEFAULT) const
	{
		typename DataPlaneSurface<T>::AtomType								temp;
		const typename DataPlaneSurface<T>::DataInnerType&					InnerData(DataPlaneSurface<T>::GetInnerData());
		typename DataPlaneSurface<T>::DataInnerType::const_iterator			const pivOrg(InnerData.begin());
		typename DataPlaneSurface<T>::DataInnerType::MetricType				const ArrMetric(InnerData.getPtrMetric());
		typename DataPlaneSurface<T>::DataInnerType::const_iterator			piv;
		int		YBound,XBound,YOff(1),XOff(1);

		DataPlaneSurface<T>::Bounds_.SetUseBounds(UseBounds);
		Sum = SumSq = 0.0;
		if(UseBounds)
		{
			YBound = DataPlaneSurface<T>::Bounds_.YBound_;
			XBound = DataPlaneSurface<T>::Bounds_.XBound_;
			if(YBound == ((DataPlaneSurface<T>::YSz_ >> 1) + 1))
			{++YOff;}
			if(XBound == ((DataPlaneSurface<T>::XSz_ >> 1) + 1))
			{++XOff;}
			if((YOff==2) && (XOff==2))
			{UseBounds = UB_NOUSE;}
		}
		if(!UseBounds)		
		{
			typename DataPlaneSurface<T>::DataStorageType::const_iterator	const pivEnd(pivOrg + (ArrMetric * DataPlaneSurface<T>::YSz_));
// TODO OpenMP
			for(piv = pivOrg;piv != pivEnd; ++piv)
			{
				if(select & SEL_SUM){Sum += *piv;}
				if(select & SEL_SUMSQ){SumSq += (*piv * *piv);}
			}
			return;
		}
// TODO OpenMP
		for(int j = -YBound+YOff;j !=YBound;++j)
		{
// TODO OpenMP
			if(j<0)
			{piv = (pivOrg + (ArrMetric * (DataPlaneSurface<T>::YSz_ + j)));}
			else
			{piv = (pivOrg + (ArrMetric * j));}
// TODO OpenMP
			for(int i = -XBound+XOff;i !=XBound;++i)
			{
				if(i<0)
				{temp = *(piv + (DataPlaneSurface<T>::XSz_ + i));}
				else
				{temp = *(piv + i);}

				if(select & SEL_SUM){Sum += temp;}
				if(select & SEL_SUMSQ){SumSq += (temp * temp);}

			}
		}
	}

	inline	RealPlane Clone(void) const
	{
		RealPlane temp(*this);
		DataPlaneSurface<T> DummyGcc(temp.DataPlaneSurface<T>::Clone());
		temp.DataPlaneSurface<T>::Swap(DummyGcc);
		return temp;
	}
	inline	PlaneBoundsType	GetMaxBounds(void) const
	{
		return PlaneBoundsType((DataPlaneSurface<T>::YSz_>> 1) + 1,(DataPlaneSurface<T>::XSz_ >> 1) + 1,DataPlaneSurface<T>::Bounds_.BoundFactor_);
	}
	template<typename T1>
	RealPlane& Transform(T1& functor,int swap,UseBoundsType UseBounds=UB_DEFAULT)
	{
		typename DataPlaneSurface<T>::DataInnerType&			InnerData(DataPlaneSurface<T>::GetInnerData());
		typename DataPlaneSurface<T>::DataInnerType::iterator	const pivOrg(InnerData.begin());
		typename DataPlaneSurface<T>::DataInnerType::MetricType	const ArrMetric(InnerData.getPtrMetric());
		typename DataPlaneSurface<T>::DataInnerType::iterator	piv;
		int														YBound,XBound,offset,offsetAux,initOffsetY,initOffsetX;

		DataPlaneSurface<T>::Bounds_.SetUseBounds(UseBounds);
		if(UseBounds)
		{
			YBound = DataPlaneSurface<T>::Bounds_.YBound_;
			XBound = DataPlaneSurface<T>::Bounds_.XBound_;
			if(YBound == ((DataPlaneSurface<T>::YSz_ >> 1) + 1))
			{initOffsetY = 2;}
			else
			{initOffsetY = 1;}
			if(XBound == ((DataPlaneSurface<T>::XSz_ >> 1) + 1))
			{initOffsetX = 2;}
			else
			{initOffsetX = 1;}
		}
		else
		{
			YBound = ((DataPlaneSurface<T>::YSz_ >> 1) + 1);
			XBound = ((DataPlaneSurface<T>::XSz_ >> 1) + 1);
			initOffsetY = 2;initOffsetX = 2;
		}
// TODO OpenMP
		for(int j = -YBound+initOffsetY;j !=YBound;++j)
		{
			if(swap)
			{
				offset = (j<0?(DataPlaneSurface<T>::YSz_ + j):j);
			}
			else
			{
				offset = (ArrMetric * (j<0?(DataPlaneSurface<T>::YSz_ + j):j));
			}
// TODO OpenMP
			for(int i = -XBound+initOffsetX;i !=XBound;++i)
			{
				if(swap)
				{
					offsetAux = offset + ArrMetric * (i<0?(DataPlaneSurface<T>::XSz_ + i):i);
				}
				else
				{
					offsetAux = offset + (i<0?(DataPlaneSurface<T>::XSz_ + i):i);
				}

				piv = pivOrg + offsetAux;
				functor(*piv,offsetAux);
			}
		}
		return  *this;
	}

	RealPlane& CleanUp(const PlaneBoundsType& Limits)
	{
		if((DataPlaneSurface<T>::Bounds_.DefaultUseBounds_ == UB_NOUSE) ||
		   ((DataPlaneSurface<T>::Bounds_.YBound_ >= Limits.YBound_) &&
		   (DataPlaneSurface<T>::Bounds_.XBound_ >= Limits.XBound_)))
		{return *this;}

		typename DataPlaneSurface<T>::DataInnerType&						InnerData(DataPlaneSurface<T>::GetInnerData());
		typename DataPlaneSurface<T>::DataInnerType::iterator				const pivOrg(InnerData.begin());
		typename DataPlaneSurface<T>::DataInnerType::MetricType				const ArrMetric(InnerData.getPtrMetric());
		typename DataPlaneSurface<T>::DataInnerType::iterator				piv;
		int		YBound,XBound,YOff(1),XOff(1);

		YBound = DataPlaneSurface<T>::Bounds_.YBound_;
		XBound = DataPlaneSurface<T>::Bounds_.XBound_;

		if(Limits.YBound_ == ((DataPlaneSurface<T>::YSz_ >> 1) + 1))
		{++YOff;}
		if(Limits.XBound_ == ((DataPlaneSurface<T>::XSz_ >> 1) + 1))
		{++XOff;}

// TODO OpenMP
		for(int j = -(Limits.YBound_)+YOff;j !=Limits.YBound_;++j)
		{
// TODO OpenMP
			const int absj(std::abs(j));
			if(j<0)
			{piv = (pivOrg + (ArrMetric * (DataPlaneSurface<T>::YSz_ + j)));}
			else
			{piv = (pivOrg + (ArrMetric * j));}
// TODO OpenMP
			for(int i = -(Limits.XBound_)+XOff;i != Limits.XBound_;++i)
			{
				if((absj < YBound) && (std::abs(i) < XBound))
					continue;
				if(i<0)
				{*(piv + (DataPlaneSurface<T>::XSz_ + i)) = 0.0;}
				else
				{*(piv + i) = 0.0;}
			}
		}		
		return *this;	
	}

	// Is not using bounds
	// Use only after reading the maps
	RealPlane& RemoveBias(T& bias)
	{
		const typename DataPlaneSurface<T>::DataInnerType&					InnerArray(DataPlaneSurface<T>::GetInnerData());
		typename DataPlaneSurface<T>::DataInnerType::const_iterator			pixPiv(InnerArray.begin());
		typename DataPlaneSurface<T>::DataInnerType::const_iterator	const	pixEnd(pixPiv + InnerArray.getSz());

//TODO OpenMP
		bias = ((T)0.0);
		for(;pixPiv != pixEnd;++pixPiv)
		{bias	+= *pixPiv;}
		bias /= static_cast<T>(InnerArray.getSz());

//TODO OpenMP
		typename DataPlaneSurface<T>::DataInnerType::iterator	pixPiv2(DataPlaneSurface<T>::GetInnerData().begin());
		for(;pixPiv2 != pixEnd;++pixPiv2)
		{*pixPiv2 -= bias;}

		return  *this;
	}

	inline T&	operator() (long line,long column)
	{
		const int Sz(DataPlaneSurface<T>::YSz_ >> 1);
		if((line < -Sz) && (line > Sz))		line	= Sz;
		if((column < -Sz) && (column > Sz))	column	= Sz;
		if(line < 0)						line	+=	DataPlaneSurface<T>::YSz_;
		if(column < 0)						column	+=	DataPlaneSurface<T>::YSz_;

		return ((DataPlaneSurface<T>::GetInnerData()).begin())[line * DataPlaneSurface<T>::YSz_ + column]; 
	}

	inline const T&	operator() (long line,long column) const
	{
		const int Sz(DataPlaneSurface<T>::YSz_ >> 1);
		if((line < -Sz) && (line > Sz))		line	= Sz;
		if((column < -Sz) && (column > Sz))	column	= Sz;
		if(line < 0)						line	+=	DataPlaneSurface<T>::YSz_;
		if(column < 0)						column	+=	DataPlaneSurface<T>::YSz_;

		return ((DataPlaneSurface<T>::GetInnerData()).begin())[line * DataPlaneSurface<T>::YSz_ + column]; 
	}

	//------------------------	
	template<typename T3,typename T2>
	friend RealPlane<T3> MultiplyInPlace(RealPlane<T3>& fst,const RealPlane<T2>& sec,UseBoundsType UseBounds);
	template<typename T3,typename T2>
	friend RealPlane<T3> AddInPlace(RealPlane<T3>& fst, RealPlane<T2>& sec,UseBoundsType UseBounds);
	template<typename T1>
	friend T1 DotProduct(const RealPlane<T1>& fst,const RealPlane<T1>& sec,UseBoundsType UseBounds);
	template<typename T1>
	friend T1 Correlation(const RealPlane<T1>& MainSurf,const RealPlane<T1>& SecSurf,int YPix,int XPix,UseBoundsType UseBounds);
	template<typename T2,typename T1>
	friend RealPlane<T2>& Transform(T1& functor,RealPlane<T2>& MSurf,int YPix,int XPix,RealPlane<T2>& patch,bool WhichChange,UseBoundsType UseBounds);
	template<typename T1>
	friend RealPlane<T1>& CopyValuesFrom(RealPlane<T1>& lhs,const RealPlane<T1>& rhs,UseBoundsType UseBounds);
};

//------------------------------------------------------------------------------------------------------------------
//Implementations
//
template<typename T1,typename T2>
FourierPlane<T1> MultiplyInPlace(FourierPlane<T1>& fst,const FourierPlane<T2>& sec,UseBoundsType UseBounds=UB_DEFAULT)
{
	if((fst.YSz_ != sec.YSz_) || (fst.XSz_ != sec.XSz_))
		FourierPlane<T1>::errDataPlaneSurface(ERROR_COD_DTSURFDIFOBJS,ERROR_MSG_DTSURFDIFOBJS,L"Multiply");

	typename FourierPlane<T1>::DataStorageType::MetricType			const ArrMetric(fst.GetInnerData().getPtrMetric());
	typename FourierPlane<T1>::DataStorageType::iterator			pivFst(fst.GetInnerData().begin());
	typename FourierPlane<T2>::DataStorageType::const_iterator		pivSec(sec.GetInnerData().begin());
	typename FourierPlane<T1>::DataStorageType::iterator			pivFstAux;
	typename FourierPlane<T2>::DataStorageType::const_iterator		pivSecAux;
	int														YBound,XBound,i,j;
	PlaneBoundsType											tBounds(GetMinMaxBounds(fst,sec,true));

	tBounds.SetUseBounds(UseBounds);

	if(UseBounds)
	{
		YBound = tBounds.YBound_;
		XBound = tBounds.XBound_;
		if((YBound == ((fst.YSz_ >> 1) + 1)) && (XBound == fst.XSz_))
		{UseBounds = UB_NOUSE;}
	}

	if(!UseBounds)
	{
		typename FourierPlane<T1>::DataStorageType::const_iterator	const pivFstEnd(pivFst + (ArrMetric * fst.YSz_));
// TODO OpenMP
		for(;pivFst != pivFstEnd; ++pivFst,++pivSec)
		{*pivFst *= *pivSec;}
		fst.SetBounds(tBounds);
		return fst;
	}
// TODO OpenMP
	for(j=0;j<YBound;++j,pivFst += ArrMetric,pivSec += ArrMetric)
	{
// TODO OpenMP
		for(i=0,pivFstAux = pivFst,pivSecAux = pivSec;
			i<XBound;
			++i,++pivFstAux,++pivSecAux)
		{*pivFstAux *= *pivSecAux;}
	}
	if(YBound == ((fst.YSz_ >> 1) + 1))
	{--YBound;}
	pivFst = fst.GetInnerData().begin() + (ArrMetric * (fst.YSz_ - 1));
	pivSec = sec.GetInnerData().begin() + (ArrMetric * (sec.YSz_ - 1));
// TODO OpenMP
	for(j=1;j<YBound;++j,pivFst -= ArrMetric,pivSec -= ArrMetric)
	{
// TODO OpenMP
		for(i=0,pivFstAux = pivFst,pivSecAux = pivSec;
			i<XBound;
			++i,++pivFstAux,++pivSecAux)
		{*pivFstAux *= *pivSecAux;}
	}
	fst.SetBounds(tBounds);
	return fst;
}

//
template<typename T3,typename T2,typename T1>
FourierPlane<T3> Multiply(Loki::Type2Type<T3>,const FourierPlane<T1>& fst,const FourierPlane<T2>& sec,UseBoundsType UseBounds=UB_DEFAULT)
{
	if((fst.YSz_ != sec.YSz_) || (fst.XSz_ != sec.XSz_))
		FourierPlane<T1>::errDataPlaneSurface(ERROR_COD_DTSURFDIFOBJS,ERROR_MSG_DTSURFDIFOBJS,L"Multiply");

	typename FourierPlane<T1>::DataStorageType::MetricType			const ArrMetric(fst.GetInnerData().getPtrMetric());
	typename FourierPlane<T1>::DataStorageType::const_iterator		pivFst(fst.GetInnerData().begin());
	typename FourierPlane<T2>::DataStorageType::const_iterator		pivSec(sec.GetInnerData().begin());
	typename FourierPlane<T1>::DataStorageType::const_iterator		pivFstAux;
	typename FourierPlane<T2>::DataStorageType::const_iterator		pivSecAux;
	typename FourierPlane<T3>::DataStorageType::iterator			pivDest,pivDestAux;
	int																YBound,XBound,i,j;
	PlaneBoundsType													tBounds(GetMinMaxBounds(fst,sec,true));
	FourierPlane<T3>												tempSurf(fst.XSz_ - 1,true);

	tBounds.SetUseBounds(UseBounds);

	if(UseBounds)
	{
		YBound = tBounds.YBound_;
		XBound = tBounds.XBound_;
		if((YBound == ((fst.YSz_ >> 1) + 1)) && (XBound == fst.XSz_))
		{UseBounds = UB_NOUSE;}
	}

	if(!UseBounds)
	{
		typename FourierPlane<T1>::DataStorageType::const_iterator	const pivFstEnd(pivFst + (ArrMetric * fst.YSz_));
		pivDestAux = tempSurf.GetInnerData().begin();
// TODO OpenMP
		for(;pivFst != pivFstEnd; ++pivFst,++pivSec,++pivDestAux)
		{*pivDestAux = *pivFst * *pivSec;}
		tempSurf.SetBounds(tBounds);
		return tempSurf;
	}
	pivDest = tempSurf.GetInnerData().begin();
// TODO OpenMP
	for(j=0;j<YBound;++j,pivFst += ArrMetric,pivSec += ArrMetric,pivDest += ArrMetric)
	{
// TODO OpenMP
		for(i=0,pivFstAux = pivFst,pivSecAux = pivSec,pivDestAux = pivDest;
			i<XBound;
			++i,++pivFstAux,++pivSecAux,++pivDestAux)
		{*pivDestAux = *pivFstAux * *pivSecAux;}
	}

	if(YBound == ((fst.YSz_ >> 1) + 1))
	{--YBound;}
	pivFst		= fst.GetInnerData().begin() + (ArrMetric * (fst.YSz_ - 1));
	pivSec		= sec.GetInnerData().begin() + (ArrMetric * (sec.YSz_ - 1));
	pivDest		= tempSurf.GetInnerData().begin() + (ArrMetric * (fst.YSz_ - 1));
// TODO OpenMP
	for(j=1;j<YBound;++j,pivFst -= ArrMetric,pivSec -= ArrMetric,pivDest -= ArrMetric)
	{
// TODO OpenMP
		for(i=0,pivFstAux = pivFst,pivSecAux = pivSec,pivDestAux = pivDest;
			i<XBound;
			++i,++pivFstAux,++pivSecAux,++pivDestAux)
		{*pivDestAux = *pivFstAux * *pivSecAux;}
	}
	
	tempSurf.SetBounds(tBounds);
	return tempSurf;
}

//
template<typename T,typename T2>
FourierPlane<T> AddInPlace(FourierPlane<T>& fst,const FourierPlane<T2>& sec,UseBoundsType UseBounds=UB_DEFAULT)
{
	if((fst.YSz_ != sec.YSz_) || (fst.XSz_ != sec.XSz_))
		FourierPlane<T>::errDataPlaneSurface(ERROR_COD_DTSURFDIFOBJSADD,ERROR_MSG_DTSURFDIFOBJSADD,L"Add");

	typename FourierPlane<T>::DataStorageType::MetricType			const ArrMetric(fst.GetInnerData().getPtrMetric());
	typename FourierPlane<T>::DataStorageType::iterator				pivFst(fst.GetInnerData().begin());
	typename FourierPlane<T2>::DataStorageType::const_iterator		pivSec(sec.GetInnerData().begin());
	typename FourierPlane<T>::DataStorageType::iterator				pivFstAux;
	typename FourierPlane<T2>::DataStorageType::const_iterator		pivSecAux;
	int																YBound,XBound,i,j;
	PlaneBoundsType													tBounds(GetMinMaxBounds(fst,sec,false));

	tBounds.SetUseBounds(UseBounds);

	if(UseBounds)
	{
		YBound = tBounds.YBound_;
		XBound = tBounds.XBound_;
		if((YBound == ((fst.YSz_ >> 1) + 1)) && (XBound == fst.XSz_))
		{UseBounds = UB_NOUSE;}
	}

	if(!UseBounds)
	{
		typename FourierPlane<T>::DataStorageType::const_iterator	const pivFstEnd(pivFst + (ArrMetric * fst.YSz_));
// TODO OpenMP
		for(;pivFst != pivFstEnd; ++pivFst,++pivSec)
		{*pivFst += *pivSec;}
		fst.SetBounds(tBounds);
		return fst;
	}

// TODO OpenMP
	for(j=0;j<YBound;++j,pivFst += ArrMetric,pivSec += ArrMetric)
	{
// TODO OpenMP
		for(i=0,pivFstAux = pivFst,pivSecAux = pivSec;
			i<XBound;
			++i,++pivFstAux,++pivSecAux)
		{*pivFstAux += *pivSecAux;}
	}

	if(YBound == ((fst.YSz_ >> 1) + 1))
	{--YBound;}

	pivFst = fst.GetInnerData().begin() + (ArrMetric * (fst.YSz_ - 1));
	pivSec = sec.GetInnerData().begin() + (ArrMetric * (sec.YSz_ - 1));
// TODO OpenMP
	for(j=1;j<YBound;++j,pivFst -= ArrMetric,pivSec -= ArrMetric)
	{
// TODO OpenMP
		for(i=0,pivFstAux = pivFst,pivSecAux = pivSec;
			i<XBound;
			++i,++pivFstAux,++pivSecAux)
		{*pivFstAux += *pivSecAux;}
	}

	fst.SetBounds(tBounds);

	return fst;
}

//
template<typename T,typename T2>
FourierPlane<T> Add(FourierPlane<T>& fst,const FourierPlane<T2>& sec,UseBoundsType UseBounds=UB_DEFAULT)
{
	if((fst.YSz_ != sec.YSz_) || (fst.XSz_ != sec.XSz_))
		FourierPlane<T>::errDataPlaneSurface(ERROR_COD_DTSURFDIFOBJSADD,ERROR_MSG_DTSURFDIFOBJSADD,L"Add");
	FourierPlane<T>	 tempSurf(fst.Clone());
	return AddInPlace(tempSurf,sec,UseBounds);
}

//
template<typename T2,typename T0,typename T1>
FourierPlane<T2> Transform(Loki::Type2Type<T2>,const FourierPlane<T0>& fst,T1& functor,UseBoundsType UseBounds=UB_DEFAULT,int BoundInfo=0)
{
	typename FourierPlane<T0>::DataStorageType::MetricType			const ArrMetric(fst.GetInnerData().getPtrMetric());
	typename FourierPlane<T0>::DataStorageType::const_iterator		pivFst(fst.GetInnerData().begin());
	typename FourierPlane<T0>::DataStorageType::const_iterator		pivFstAux;
	typename FourierPlane<T2>::DataStorageType::iterator			pivDestAux;
	int																YBound,XBound,i,j;
	FourierPlane<T2>												tempSurf(fst.XSz_ - 1,false);
	typename FourierPlane<T2>::DataStorageType::iterator			pivDest(tempSurf.GetInnerData().begin());

	fst.Bounds_.SetUseBounds(UseBounds);

	if(UseBounds)
	{
		YBound = fst.Bounds_.YBound_;
		XBound = fst.Bounds_.XBound_;
		if((YBound == ((fst.YSz_ >> 1) + 1)) && (XBound == fst.XSz_))
		{UseBounds = UB_NOUSE;}
	}

	if(!UseBounds)
	{
		typename FourierPlane<T0>::DataStorageType::const_iterator	const pivEnd(pivFst + (ArrMetric * fst.YSz_));
		typename FourierPlane<T0>::DataStorageType::const_iterator	const pivOrg(pivFst);

// TODO OpenMP
		for(;pivFst != pivEnd; ++pivFst,++pivDest)
		{*pivDest = functor(*pivFst);}
		return tempSurf;
	}

// TODO OpenMP
	for(j=0;j<YBound;++j,pivFst += ArrMetric,pivDest += ArrMetric)
	{
// TODO OpenMP
		for(i=0,pivFstAux = pivFst,pivDestAux = pivDest;i<XBound;++i,++pivFstAux,++pivDestAux)
		{*pivDestAux = functor(*pivFstAux);}
	}

	if(YBound == ((fst.YSz_ >> 1) + 1))
	{--YBound;}

	pivFst  = (fst.GetInnerData().begin()) + (ArrMetric * (fst.YSz_ - 1));
	pivDest = (tempSurf.GetInnerData().begin()) + (ArrMetric * (fst.YSz_ - 1));
// TODO OpenMP
	for(j=1;j<YBound;++j,pivFst -= ArrMetric,pivDest -= ArrMetric)
	{
// TODO OpenMP
		for(i=0,pivFstAux = pivFst,pivDestAux = pivDest;i<XBound;++i,++pivFstAux,++pivDestAux)
		{*pivDestAux = functor(*pivFstAux);}
	}
	return tempSurf;
}
//
template<typename T1,typename T2>
RealPlane<T1> MultiplyInPlace(RealPlane<T1>& fst,const RealPlane<T2>& sec,UseBoundsType UseBounds=UB_DEFAULT)
{
	if((fst.YSz_ != sec.YSz_) || (fst.XSz_ != sec.XSz_))
		RealPlane<T1>::errDataPlaneSurface(ERROR_COD_DTSURFDIFOBJS,ERROR_MSG_DTSURFDIFOBJS,L"Multiply");

	typename RealPlane<T1>::DataStorageType::MetricType			const ArrMetric(fst.GetInnerData().getPtrMetric());
	typename RealPlane<T1>::DataStorageType::iterator			pivFst(fst.GetInnerData().begin());
	typename RealPlane<T2>::DataStorageType::const_iterator		pivSec(sec.GetInnerData().begin());
	typename RealPlane<T1>::DataStorageType::iterator			pivFstAux;
	typename RealPlane<T2>::DataStorageType::const_iterator		pivSecAux;
	typename RealPlane<T1>::DataStorageType::iterator			pivpivFstAux;
	typename RealPlane<T2>::DataStorageType::const_iterator		pivpivSecAux;
	int															YBound,XBound;
	PlaneBoundsType										tBounds(GetMinMaxBounds(fst,sec,true));

	tBounds.SetUseBounds(UseBounds);
	if(UseBounds)
	{
		YBound = tBounds.YBound_;
		XBound = tBounds.XBound_;
		if((YBound == ((fst.YSz_ >> 1) + 1)) && (XBound == ((fst.XSz_ >> 1) + 1)))
		{UseBounds = UB_NOUSE;}
	}

	if(!UseBounds)
	{
		typename RealPlane<T1>::DataStorageType::const_iterator	const pivFstEnd(pivFst + (ArrMetric * fst.YSz_));
// TODO OpenMP
		for(;pivFst != pivFstEnd; ++pivFst,++pivSec)
		{*pivFst *= *pivSec;}
		fst.SetBounds(tBounds);
		return fst;
	}
// TODO OpenMP
	for(int j = -YBound+1;j !=YBound;++j)
	{
// TODO OpenMP
		if(j<0)
		{
			pivFstAux = (pivFst + (ArrMetric * (fst.YSz_ + j)));pivSecAux = (pivSec + (ArrMetric * (sec.YSz_ + j)));
		}
		else
		{
			pivFstAux = (pivFst + (ArrMetric * j));pivSecAux = (pivSec + (ArrMetric * j));
		}
// TODO OpenMP
		for(int i = -XBound+1;i !=XBound;++i)
		{
			if(i<0)
			{
				pivpivFstAux = pivFstAux + (fst.XSz_ + i);pivpivSecAux = pivSecAux + (sec.XSz_ + i);
			}
			else
			{
				pivpivFstAux = pivFstAux + i;pivpivSecAux = pivSecAux + i;
			}

			*pivpivFstAux *= *pivpivSecAux;
		}
	}

	fst.SetBounds(tBounds);
	return fst;
}


//
template<typename T1,typename T2>
RealPlane<T1> AddInPlace(RealPlane<T1>& fst, RealPlane<T2>& sec,UseBoundsType UseBounds=UB_DEFAULT)
{
	if((fst.YSz_ != sec.YSz_) || (fst.XSz_ != sec.XSz_))
		RealPlane<T1>::errDataPlaneSurface(ERROR_COD_DTSURFDIFOBJS,ERROR_MSG_DTSURFDIFOBJS,L"Multiply");

	typename RealPlane<T1>::DataStorageType::MetricType			const ArrMetric(fst.GetInnerData().getPtrMetric());
	typename RealPlane<T1>::DataStorageType::iterator			pivFst(fst.GetInnerData().begin());
	typename RealPlane<T2>::DataStorageType::const_iterator		pivSec(sec.GetInnerData().begin());
	typename RealPlane<T1>::DataStorageType::iterator			pivFstAux;
	typename RealPlane<T2>::DataStorageType::const_iterator		pivSecAux;
	typename RealPlane<T1>::DataStorageType::iterator			pivpivFstAux;
	typename RealPlane<T2>::DataStorageType::const_iterator		pivpivSecAux;
	int													YBound,XBound;
	PlaneBoundsType										tBounds(GetMinMaxBounds(fst,sec,false));

	tBounds.SetUseBounds(UseBounds);

	if(UseBounds)
	{
		YBound = tBounds.YBound_;
		XBound = tBounds.XBound_;
		if((YBound == ((fst.YSz_ >> 1) + 1)) && (XBound == ((fst.XSz_ >> 1) + 1)))
		{UseBounds = UB_NOUSE;}
	}

	if(!UseBounds)
	{
		typename RealPlane<T1>::DataStorageType::const_iterator	const pivFstEnd(pivFst + (ArrMetric * fst.YSz_));
// TODO OpenMP
		for(;pivFst != pivFstEnd; ++pivFst,++pivSec)
		{*pivFst += *pivSec;}
		fst.SetBounds(tBounds);
		return fst;
	}

// TODO OpenMP
	for(int j = -YBound+1;j !=YBound;++j)
	{
// TODO OpenMP
		if(j<0)
		{
			pivFstAux = (pivFst + (ArrMetric * (fst.YSz_ + j)));
			pivSecAux = (pivSec + (ArrMetric * (sec.YSz_ + j)));
		}
		else
		{
			pivFstAux = (pivFst + (ArrMetric * j));
			pivSecAux = (pivSec + (ArrMetric * j));
		}
// TODO OpenMP
		for(int i = -XBound+1;i !=XBound;++i)
		{
			if(i<0)
			{
				pivpivFstAux = pivFstAux + (fst.XSz_ + i);
				pivpivSecAux = pivSecAux + (sec.XSz_ + i);
			}
			else
			{
				pivpivFstAux = pivFstAux + i;
				pivpivSecAux = pivSecAux + i;
			}

			*pivpivFstAux += *pivpivSecAux;
		}
	}

	fst.SetBounds(tBounds);
	return fst;
}

//
template<typename T>
T DotProduct(const RealPlane<T>& fst,const RealPlane<T>& sec,UseBoundsType UseBounds=UB_DEFAULT)
{
	int const fstSz(fst.YSz_);
	int const secSz(sec.YSz_);
	
	typename RealPlane<T>::DataStorageType::MetricType		const ArrMetricShort((fstSz<secSz)?fst.GetInnerData().getPtrMetric():sec.GetInnerData().getPtrMetric());
	typename RealPlane<T>::DataStorageType::MetricType		const ArrMetricLarge((fstSz>secSz)?fst.GetInnerData().getPtrMetric():sec.GetInnerData().getPtrMetric());
	typename RealPlane<T>::DataStorageType::const_iterator	const OrgShort((fstSz<secSz)?fst.GetInnerData().begin():sec.GetInnerData().begin());
	typename RealPlane<T>::DataStorageType::const_iterator	const OrgLarge((fstSz>=secSz)?fst.GetInnerData().begin():sec.GetInnerData().begin());
	typename RealPlane<T>::DataStorageType::const_iterator	pivShort,pivLarge,pivShortAux,pivLargeAux;
	int												YBound,XBound,j,i,YOffSet,XOffSet;
	int												const SzShort((fstSz<secSz)?fstSz:secSz);
	int												const SzLarge((fstSz>=secSz)?fstSz:secSz);
	PlaneBoundsType									tBounds(GetMinMaxBounds(fst,sec,true));
	T												Sum(0.0);

	tBounds.SetUseBounds(UseBounds);

	if(UseBounds)
	{
		YBound = tBounds.YBound_;
		XBound = tBounds.XBound_;
		YOffSet = XOffSet = 1;
		if(YBound == ((SzShort >> 1) + 1))
		{++YOffSet;}
		if(XBound == ((SzShort >> 1) + 1))
		{++XOffSet;}
	}
	else
	{
		YBound = XBound = ((SzShort >> 1) + 1);
		YOffSet = XOffSet = 2;
	}


// TODO OpenMP
	for(j = -YBound+YOffSet;j!=YBound;++j)
	{
		if(j<0)
		{
			pivShort = (OrgShort + (ArrMetricShort * (SzShort + j)));
			pivLarge = (OrgLarge + (ArrMetricLarge * (SzLarge + j)));
		}
		else
		{
			pivShort = (OrgShort + (ArrMetricShort * j));
			pivLarge = (OrgLarge + (ArrMetricLarge * j));
		}
// TODO OpenMP
		for(i = -XBound+XOffSet;i!=XBound;++i)
		{
			if(i<0)
			{
				pivShortAux = pivShort + (SzShort + i);
				pivLargeAux = pivLarge + (SzLarge + i);
			}
			else
			{
				pivShortAux = pivShort + i;
				pivLargeAux = pivLarge + i;
			}

			Sum += (*pivShortAux * *pivLargeAux);
		}
	}
	return Sum;
}
//
template<typename T1>
RealPlane<T1>& CopyValuesFrom(RealPlane<T1>& lhs,const RealPlane<T1>& rhs,UseBoundsType UseBounds=UB_DEFAULT)
{
	typename RealPlane<T1>::DataStorageType::MetricType		const ArrMetricLhs(lhs.GetInnerData().getPtrMetric());
	typename RealPlane<T1>::DataStorageType::MetricType		const ArrMetricRhs(rhs.GetInnerData().getPtrMetric());
	typename RealPlane<T1>::DataStorageType::iterator		const OrgLhs(lhs.GetInnerData().begin());
	typename RealPlane<T1>::DataStorageType::const_iterator	const OrgRhs(rhs.GetInnerData().begin());
	typename RealPlane<T1>::DataStorageType::iterator		pivLhs,pivLhsAux;
	typename RealPlane<T1>::DataStorageType::const_iterator	pivRhs,pivRhsAux;

	int		const	lhsSz(lhs.YSz_);
	int		const	rhsSz(rhs.YSz_);
	int		const	Sz((((lhsSz<rhsSz)?lhsSz:rhsSz) >> 1) + 1);
	int				YBound,XBound,j,i;

	rhs.Bounds_.SetUseBounds(UseBounds);

	if(UseBounds==UB_NOUSE)
	{
		YBound = XBound = Sz;
		lhs.Bounds_.YBound_ = lhs.Bounds_.XBound_ = PlaneBoundsType::INVALID_BOUND;
	}
	else
	{
		YBound = rhs.Bounds_.YBound_;
		if(YBound > Sz) YBound = Sz;
		lhs.Bounds_.YBound_ = YBound;

		XBound = rhs.Bounds_.XBound_;		
		if(XBound > Sz) XBound = Sz;
		lhs.Bounds_.XBound_ = XBound;
	}

	int	const minBoundOffY((YBound == Sz)?2:1);
	int	const minBoundOffX((XBound == Sz)?2:1);

// TODO OpenMP
	for(j = -YBound + minBoundOffY;j!=YBound;++j)
	{
		if(j<0)
		{
			pivLhs		= (OrgLhs + (ArrMetricLhs * (lhsSz + j)));
			pivRhs		= (OrgRhs + (ArrMetricRhs * (rhsSz + j)));
		}
		else
		{
			pivLhs		= (OrgLhs + (ArrMetricLhs * j));
			pivRhs		= (OrgRhs + (ArrMetricRhs * j));
		}
// TODO OpenMP
		for(i = -XBound+minBoundOffX;i!=XBound;++i)
		{
			if(i<0)
			{
				pivLhsAux = pivLhs + (ArrMetricLhs + i);
				pivRhsAux = pivRhs + (ArrMetricRhs + i);
			}
			else
			{
				pivLhsAux	= pivLhs + i;
				pivRhsAux	= pivRhs + i;
			}
			*pivLhsAux = *pivRhsAux;
		}
	}

	return lhs;
}
//
template<typename T>
T Correlation(const RealPlane<T>& MainSurf,const RealPlane<T>& SecSurf,int YPix,int XPix,UseBoundsType UseBounds=UB_DEFAULT)
{
	if((MainSurf.YSz_  < SecSurf.YSz_) || (MainSurf.XSz_ < SecSurf.XSz_))
		RealPlane<T>::errDataPlaneSurface(ERROR_COD_SURFWRONGSZ,ERROR_MSG_SURFWRONGSZ,L"DotProduct");

	typename RealPlane<T>::DataStorageType::MetricType		const ArrMetricShort(SecSurf.GetInnerData().getPtrMetric());
	typename RealPlane<T>::DataStorageType::MetricType		const ArrMetricLarge(MainSurf.GetInnerData().getPtrMetric());
	typename RealPlane<T>::DataStorageType::const_iterator	const OrgShort(SecSurf.GetInnerData().begin());
	typename RealPlane<T>::DataStorageType::const_iterator	const OrgLarge(MainSurf.GetInnerData().begin());
	typename RealPlane<T>::DataStorageType::const_iterator	pivShort,pivLarge,pivShortAux,pivLargeAux;
	int												YBound,XBound,j,i,YOffSet,XOffSet;
	int												const SzShort(SecSurf.YSz_);
	int												const SzLarge(MainSurf.YSz_);
	T												Correlation(0.0);

	SecSurf.Bounds_.SetUseBounds(UseBounds);

	if(UseBounds)
	{
		YBound = SecSurf.Bounds_.YBound_;
		XBound = SecSurf.Bounds_.XBound_;
		YOffSet = XOffSet = 1;
		if(YBound == ((SzShort >> 1) + 1))
		{++YOffSet;}
		if(XBound == ((SzShort >> 1) + 1))
		{++XOffSet;}
	}
	else
	{
		YBound = XBound = ((SzShort >> 1) + 1);
		YOffSet = XOffSet = 2;
	}

	int YMinBound(-YBound+YOffSet + YPix);
	int YMaxBound(YBound + YPix);
	int XMinBound(-XBound+XOffSet + XPix);
	int XMaxBound(XBound + XPix);

	if(YMinBound < 0) YMinBound = 0;
	if(YMaxBound > MainSurf.YSz_) YMaxBound = SzLarge;
	YMinBound -= YPix;
	YMaxBound -= YPix;

	if(XMinBound < 0) XMinBound = 0;
	if(XMaxBound > MainSurf.XSz_) XMaxBound = SzLarge;
	XMinBound -= XPix;
	XMaxBound -= XPix;
// TODO OpenMP
	for(j = YMinBound;j!=YMaxBound;++j)
	{
		if(j<0)
		{pivShort = (OrgShort + (ArrMetricShort * (SzShort + j)));}
		else
		{pivShort = (OrgShort + (ArrMetricShort * j));}

		pivLarge = (OrgLarge + (ArrMetricLarge * (j + YPix)));
// TODO OpenMP
		for(i = XMinBound;i!=XMaxBound;++i)
		{
			if(i<0)
			{pivShortAux = pivShort + (SzShort + i);}
			else
			{pivShortAux = pivShort + i;}
			
			pivLargeAux = pivLarge + (i + XPix);

			Correlation += (*pivShortAux * *pivLargeAux);
		}
	}
	return Correlation;
}
//
template<typename T0,typename T1>	
PlaneBoundsType GetMinMaxBounds(const DataPlaneSurface<T0>& fst,const DataPlaneSurface<T1>& sec ,bool MinMax)
{
	if(!(fst.HasBounds()) || !(sec.HasBounds()))
		return PlaneBoundsType();

	PlaneBoundsType tBounds;
	if(MinMax)
	{
		tBounds.XBound_ = ((fst.Bounds_.XBound_ < sec.Bounds_.XBound_) ? fst.Bounds_.XBound_ : sec.Bounds_.XBound_);
		tBounds.YBound_ = ((fst.Bounds_.YBound_ < sec.Bounds_.YBound_) ? fst.Bounds_.YBound_ : sec.Bounds_.YBound_);
	}
	else
	{
		tBounds.XBound_ = ((fst.Bounds_.XBound_ > sec.Bounds_.XBound_) ? fst.Bounds_.XBound_ : sec.Bounds_.XBound_);
		tBounds.YBound_ = ((fst.Bounds_.YBound_ > sec.Bounds_.YBound_) ? fst.Bounds_.YBound_ : sec.Bounds_.YBound_);
	}
	return tBounds;
}
//
template<typename T,typename T1>
RealPlane<T>& Transform(T1& functor,RealPlane<T>& MSurf,int YPix,int XPix,RealPlane<T>& patch,bool WhichChange,UseBoundsType UseBounds=UB_DEFAULT)
{

	if((MSurf.YSz_  < patch.YSz_) || (MSurf.XSz_ < patch.XSz_))
		RealPlane<T>::errDataPlaneSurface(ERROR_COD_SURFWRONGSZ,ERROR_MSG_SURFWRONGSZ,L"Subtract");

	typename RealPlane<T>::DataStorageType::MetricType		const ArrMetricShort(patch.GetInnerData().getPtrMetric());
	typename RealPlane<T>::DataStorageType::MetricType		const ArrMetricLarge(MSurf.GetInnerData().getPtrMetric());
	typename RealPlane<T>::DataStorageType::iterator			const OrgShort(patch.GetInnerData().begin());
	typename RealPlane<T>::DataStorageType::iterator			const OrgLarge(MSurf.GetInnerData().begin());
	typename RealPlane<T>::DataStorageType::iterator			pivShort,pivShortAux;
	typename RealPlane<T>::DataStorageType::iterator			pivLargeAux,pivLarge;
	int								YBound,XBound,j,i;
	int								const SzShort(patch.YSz_);
	int								const SzLarge(MSurf.YSz_);
	T								temp;

	patch.Bounds_.SetUseBounds(UseBounds);

	if(UseBounds)
	{
		YBound = patch.Bounds_.YBound_;
		XBound = patch.Bounds_.XBound_;
	}
	else
	{
		YBound = XBound = ((SzShort >> 1) + 1);
	}

	int YMinBound(-YBound+1 + YPix);
	int YMaxBound(YBound + YPix);
	int XMinBound(-XBound+1 + XPix);
	int XMaxBound(XBound + XPix);

	if(YMinBound < 0) YMinBound = 0;
	if(YMaxBound > MSurf.YSz_) YMaxBound = MSurf.YSz_;
	YMinBound -= YPix;
	YMaxBound -= YPix;

	if(XMinBound < 0) XMinBound = 0;
	if(XMaxBound > MSurf.XSz_) XMaxBound = MSurf.XSz_;
	XMinBound -= XPix;
	XMaxBound -= XPix;
// TODO OpenMP
	for(j = YMinBound;j!=YMaxBound;++j)
	{
		if(j<0)
		{pivShort = (OrgShort + (ArrMetricShort * (SzShort + j)));}
		else
		{pivShort = (OrgShort + (ArrMetricShort * j));}

		pivLarge = (OrgLarge + (ArrMetricLarge * (j + YPix)));
// TODO OpenMP
		for(i = XMinBound;i!=XMaxBound;++i)
		{
			if(i<0)
			{pivShortAux = pivShort + (SzShort + i);}
			else
			{pivShortAux = pivShort + i;}
			
			pivLargeAux = pivLarge + (i + XPix);

			temp = functor(*pivShortAux,*pivLargeAux);

			if(WhichChange)
			{*pivShortAux = temp;}
			else
			{*pivLargeAux = temp;}
		}
	}
	if(WhichChange)
		return patch;
	return MSurf;
}
//
template<typename T>
class ExtractSurfFromVectorsFunct
{
public:
	inline ExtractSurfFromVectorsFunct(const FourierPlane<AlgVect<T> >& Vects,int Sz,int Plane)
		:Surf_(Sz,Vects.GetBounds(),false),Sz_(Sz),Org_(Surf_.GetInnerData().begin()),Plane_(Plane)
	{
		if(Plane >= (*(Vects.GetInnerData().begin())).GetEnd())
		{
			throw Zeus::libException(ERROR_COD_VECTOROUTOFBOUNDS,ERROR_MSG_VECTOROUTOFBOUNDS,std::wstring(L"ExtractSurfFromVectorsFunctor"));
		}

	}
	inline void operator()(AlgVect<T>& v,int Offset)
	{
		*(Org_ + Offset) = *(v.GetInnerData().begin() + Plane_);
	}
	inline FourierPlane<T> GetSurface(void) const
	{
		return Surf_;
	}
private:
	int									Sz_;
	FourierPlane<T>						Surf_;
	typename FourierPlane<T>::DataStorageType::iterator	const Org_;
	int									Plane_;
};


} // end of namespace Zeus
#endif //ZEUS_WORKSPACEH
