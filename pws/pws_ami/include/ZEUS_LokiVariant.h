////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2001, 2002 by Andrei Alexandrescu
// Permission to use, copy, modify, distribute and sell this software for any 
//     purpose is hereby granted without fee, provided that the above copyright 
//     notice appear in all copies and that both that copyright notice and this 
//     permission notice appear in supporting documentation.
// The author makes no representations about the suitability of this software 
//     for any purpose. It is provided "as is" 
//     without express or implied warranty.
////////////////////////////////////////////////////////////////////////////////

#ifndef VARIANT_INC_
#define VARIANT_INC_

#include <typeinfo>
#include <stdexcept>
#include "Typelist.h"
#include "static_check.h"
#include "ZEUS_Exceptions.h"

template <class TList> struct MaxSize;

template <> 
struct MaxSize<Loki::NullType>
{
    enum { result = 0 };
};

template <class Head, class Tail>
struct MaxSize< Loki::Typelist<Head, Tail> >
{
private:
    enum { tailResult = int(MaxSize<Tail>::result) };
public:
    enum { result = sizeof(Head) > int(tailResult)
        ? sizeof(Head)
        : int(tailResult) };
};

template <typename U> struct Structify { U dummy_; };

template <class TList> struct MakeConst;

template <> struct MakeConst<Loki::NullType>
{
    typedef Loki::NullType Result; // terminator is not const
};

template <typename Head, class Tail>
struct MakeConst< Loki::Typelist<Head, Tail> >
{
private:
    typedef typename MakeConst<Tail>::Result NewTail;
public:
    typedef Loki::Typelist<const Head, NewTail> Result; 
};
////////////////////////////////////////////////////////////////////////////////
// class template Variant
// Implements a discriminated union in C++
////////////////////////////////////////////////////////////////////////////////

template <class TList, class Align =  Structify<void (std::string::*)(void)> >
class Variant
{
private:
	inline void errVariant(void) const
	{throw Zeus::libException(ERROR_COD_ZEUSVARIANTTY,ERROR_MSG_ZEUSVARIANTTY,*this);}

public:
    typedef TList Types;

    // Default constructor
    // Initializes the Variant with a default-constructed object of the 
    // first type in the typelist
    Variant()
    {
        typedef typename TList::Head T;
        new(&buffer_[0]) T;
        lokiVptr_ = &VTableImpl<T>::vTbl_;
    }

    // Copy constructor
    Variant(const Variant& rhs)
    {
        (rhs.lokiVptr_->clone_)(rhs, *this);
    }

    // Converting constructor; accepts any type in the typelist
    // @@@ Suggested simple improvement: accept any type convertible to one of 
    // the types in the typelist. Use Loki::Conversion as a building block
    template <class T>
    explicit Variant(const T& val)
    {
        LOKI_STATIC_CHECK((Loki::TL::IndexOf<TList, T>::value >= 0), 
            Invalid_Type_Used_As_Initializer);
        
        new(&buffer_[0]) T(val);
        lokiVptr_ = &VTableImpl<T>::vTbl_;
    }
    // Create a new Variant of the same type as the source, without copying the 
    // destination into it. The newly created object is default constructed.
    enum CloneTypeOnly { cloneTypeOnly };
    Variant(const Variant& rhs, CloneTypeOnly)
    {
        (rhs.lokiVptr_->cloneTypeOnly_)(rhs, *this);
    }

    // Canonic assignment operator
    Variant& operator=(const Variant& rhs)
    {
        Variant temp(rhs);
        temp.Swap(*this);
		return *this;
    }
    
    // Assignment operator from one of the allowed types
    // This is necessary because the constructor is explicit
    template <class T> 
    Variant& operator=(const T& rhs)
    {
        Variant temp(rhs);
        temp.Swap(*this);
		return *this;
    }
    
    // Assignment from another Variant instantiation
    template <class TList2, typename Align2>
    Variant& operator=(const Variant<TList2, Align2>& rhs)
    {
        Variant temp(rhs);
        temp.Swap(*this);
		return *this;
    }
    
    // ~
    ~Variant()
    {
        (lokiVptr_->destroy_)(*this);
    }

    struct VTable;

    // VTable concrete implementations
    // VTableImpl<T> contains definitions for all of a VTable's pointer to
    // functions.

    template <class T>
    struct VTableImpl
    {
        static const std::type_info& TypeId()
        {
            return typeid(T);
        }
        
        static void Destroy(const Variant& var)
        {
            const T& data = *reinterpret_cast<const T*>(&var.buffer_[0]);
            data.~T();
        }
        
        static void Swap(void* lhs, void* rhs)
        {
            using namespace std;
            swap(*static_cast<T*>(lhs), *static_cast<T*>(rhs));
        }

        static void Clone(const Variant& src, Variant& dest)
        {
            new(&dest.buffer_[0]) T(
                *reinterpret_cast<const T*>(&src.buffer_[0]));
            dest.lokiVptr_ = src.lokiVptr_;
        }
        
        static void CloneTypeOnly(const Variant& src, Variant& dest)
        {
            new(&dest.buffer_[0]) T;
            dest.lokiVptr_ = src.lokiVptr_;
        }
        static VTable vTbl_;
    };

    // VTable structure
    // The essential component of the fake vtable idiom, VTable contains
    // pointers to functions, pointers that will be filled up with addresses
    // of functions in VTableImpl

    struct VTable
    {
        const std::type_info& (*typeId_)();
        void (*destroy_)(const Variant&);
        void (*clone_)(const Variant&, Variant&);
        void (*cloneTypeOnly_)(const Variant&, Variant&);
        void (*swap_)(void* lhs, void* rhs);
    };

public:   // should be private; some compilers prefer 'public' :o}

    enum { neededSize = MaxSize<TList>::result };

    VTable* lokiVptr_;
    union
    {
        Align dummy_;
        char buffer_[neededSize];
    };
    
public:
    void Swap(Variant& rhs)
    {
        // Create a temporary storage
        Variant temp(rhs, cloneTypeOnly);
        // Move rhs into temp
        (rhs.lokiVptr_->swap_)(rhs.buffer_, temp.buffer_);
        rhs.~Variant();
        // Move *this into rhs
        new(&rhs) Variant(*this, cloneTypeOnly);
        (lokiVptr_->swap_)(buffer_, rhs.buffer_);
        this->~Variant();
        // Move temp into *this
        new(this) Variant(temp, cloneTypeOnly);
        (lokiVptr_->swap_)(buffer_, temp.buffer_);
    }
    
    const std::type_info& TypeId() const
    {
        return (lokiVptr_->typeId_)();
    }
    
    template <typename T> T* GetPtr()
    {
        return TypeId() == typeid(T) 
            ? reinterpret_cast<T*>(&buffer_[0])
            : 0;
    }
    
    template <typename T> const T* GetPtr() const
    {
        return TypeId() == typeid(T) 
            ? reinterpret_cast<const T*>(&buffer_[0])
            : 0;
    }
    
    template <typename T> T& Get()
    {
        T* p = GetPtr<T>();
		if (!p) errVariant();
        return *p;
    }
    
    template <typename T> const T& Get() const
    {
        const T* p = GetPtr<T>();
        if (!p) errVariant();
        return *p;
    }
};

template <class TList, class Align> 
template <class T>
typename Variant<TList, Align>::VTable 
Variant<TList, Align>::VTableImpl<T>::vTbl_ = 
{
    &VTableImpl<T>::TypeId,
    &VTableImpl<T>::Destroy,
    &VTableImpl<T>::Clone,
    &VTableImpl<T>::CloneTypeOnly,
    &VTableImpl<T>::Swap
};

#endif
