#ifndef LOKI_THREADS_H_
#define LOKI_THREADS_H_

////////////////////////////////////////////////////////////////////////////////
// macro LOKI_DEFAULT_THREADING
// Selects the default threading model for certain components of Loki
// If you don't define it, it defaults to single-threaded
// All classes in Loki have configurable threading model; LOKI_DEFAULT_THREADING
// affects only default template arguments
////////////////////////////////////////////////////////////////////////////////

// $Header: /planck/cvs/planck/Level2/Task_pkg/HL2_PowellSnakes/src/include/Threads.h,v 3.51 2010/06/17 09:56:06 pcarvalho Exp $


#if defined(WIN32) && (defined(LOKI_CLASS_LEVEL_THREADING) || defined(LOKI_OBJECT_LEVEL_THREADING))
	// threads only on windows 
	#include <Windows.h> 
    #define LOKI_DEFAULT_THREADING_NO_OBJ_LEVEL ::Loki::ClassLevelLockable
    #if defined(LOKI_CLASS_LEVEL_THREADING) 
        #define LOKI_DEFAULT_THREADING ::Loki::ClassLevelLockable
    #else
        #define LOKI_DEFAULT_THREADING ::Loki::ObjectLevelLockable
    #endif
#else
    #define LOKI_DEFAULT_THREADING ::Loki::SingleThreaded
	#define LOKI_DEFAULT_THREADING_NO_OBJ_LEVEL ::Loki::SingleThreaded
#endif

#include <cassert>

namespace Loki
{
////////////////////////////////////////////////////////////////////////////////
// class template SingleThreaded
// Implementation of the ThreadingModel policy used by various classes
// Implements a single-threaded model; no synchronization
////////////////////////////////////////////////////////////////////////////////

    template <class Host>
    class SingleThreaded
    {
    public:
        struct Lock
        {
            Lock() {}
            explicit Lock(const SingleThreaded&) {}
        };
        
		void DoLock() {}
		void DoUnLock() {}

        typedef Host VolatileType;

        typedef int IntType; 

        static IntType AtomicAdd(volatile IntType& lval, IntType val)
        { return lval += val; }
        
        static IntType AtomicSubtract(volatile IntType& lval, IntType val)
        { return lval -= val; }

        static IntType AtomicMultiply(volatile IntType& lval, IntType val)
        { return lval *= val; }
        
        static IntType AtomicDivide(volatile IntType& lval, IntType val)
        { return lval /= val; }
        
        static IntType AtomicIncrement(volatile IntType& lval)
        { return ++lval; }
        
        static IntType AtomicDecrement(volatile IntType& lval)
        { return --lval; }
        
        static void AtomicAssign(volatile IntType & lval, IntType val)
        { lval = val; }
        
        static void AtomicAssign(IntType & lval, volatile IntType & val)
        { lval = val; }
    };
    
#if defined(_WINDOWS_) || defined(_WINDOWS_H) 


////////////////////////////////////////////////////////////////////////////////
// class template ObjectLevelLockable
// Implementation of the ThreadingModel policy used by various classes
// Implements a object-level locking scheme
////////////////////////////////////////////////////////////////////////////////

    template <class Host>
    class ObjectLevelLockable
    {
       mutable CRITICAL_SECTION mtx_;
    public:
        ObjectLevelLockable()
        {
			::InitializeCriticalSection(&mtx_);
        }
        
        ObjectLevelLockable(const ObjectLevelLockable&)
        {
            ::InitializeCriticalSection(&mtx_);
        }

		void DoLock() {::EnterCriticalSection(&mtx_);}
		void DoUnLock() {::LeaveCriticalSection(&mtx_);}

         ~ObjectLevelLockable()
        {
            ::DeleteCriticalSection(&mtx_);
        }

        class Lock;
        friend class Lock;
        
        class Lock
        {
            ObjectLevelLockable const& host_;
            
            Lock(const Lock&);
            Lock& operator=(const Lock&);
        public:
            explicit Lock(const ObjectLevelLockable& host) : host_(host)
            {
                ::EnterCriticalSection(&host_.mtx_);
            }

            ~Lock()
            {
              ::LeaveCriticalSection(&host_.mtx_);
            }
        };FFTW_EXCEPTIONSH

        typedef volatile Host VolatileType;

        typedef LONG IntType; 

        static IntType AtomicIncrement(volatile IntType& lval)
        { return InterlockedIncrement(&const_cast<IntType&>(lval)); }
        
        static IntType AtomicDecrement(volatile IntType& lval)
        { return InterlockedDecrement(&const_cast<IntType&>(lval)); }

        static void AtomicAssign(volatile IntType& lval, IntType val)
        { InterlockedExchange(&const_cast<IntType&>(lval), val); }
        
        static void AtomicAssign(IntType& lval, volatile IntType& val)
        { InterlockedExchange(&lval, val); }
    };

    template<typename Host>
    class ClassLevelLockable
    {
        struct Initializer
        {   
            CRITICAL_SECTION mtx_;
            bool init_;

            Initializer():init_(false)
            {
                ::InitializeCriticalSection(&mtx_);
                init_=true;
            }
            ~Initializer()
            {
                assert(init_);
                ::DeleteCriticalSection(&mtx_);
            }
        };
        
        static Initializer initializer_;

    public:
		void DoLock()
		{
			assert(initializer_.init_);
			::EnterCriticalSection(&initializer_.mtx_);
		}
		void DoUnLock()
		{
            assert(initializer_.init_);
            ::LeaveCriticalSection(&initializer_.mtx_);
		}

        class Lock;
        friend class Lock;
        
        class Lock
        {
            Lock(const Lock&);
            Lock& operator=(const Lock&);
        public:
            Lock()
            {
                assert(initializer_.init_);
                ::EnterCriticalSection(&initializer_.mtx_);
            }
            explicit Lock(const ClassLevelLockable&)
            {
                assert(initializer_.init_);
                ::EnterCriticalSection(&initializer_.mtx_);
            }
            ~Lock()
            {
                assert(initializer_.init_);
                ::LeaveCriticalSection(&initializer_.mtx_);
            }
        };

        typedef volatile Host VolatileType;

        typedef LONG IntType; 

        static IntType AtomicIncrement(volatile IntType& lval)
        { return InterlockedIncrement(&const_cast<IntType&>(lval)); }
        
        static IntType AtomicDecrement(volatile IntType& lval)
        { return InterlockedDecrement(&const_cast<IntType&>(lval)); }
        
        static void AtomicAssign(volatile IntType& lval, IntType val)
        { InterlockedExchange(&const_cast<IntType&>(lval), val); }
        
        static void AtomicAssign(IntType& lval, volatile IntType& val)
        { InterlockedExchange(&lval, val); }
    };

    template<typename Host>
    typename ClassLevelLockable<Host>::Initializer 
    ClassLevelLockable<Host>::initializer_;
    
#else

////////////////////////////////////////////////////////////////////////////////
// class template ObjectLevelLockable
// Implementation of the ThreadingModel policy used by various classes
// Implements a object-level locking scheme
////////////////////////////////////////////////////////////////////////////////

    template <class Host>
    class ObjectLevelLockable
    {
       //mutable CRITICAL_SECTION mtx_;
    public:
        ObjectLevelLockable()
        {
			//::InitializeCriticalSection(&mtx_);
        }
        
        ObjectLevelLockable(const ObjectLevelLockable&)
        {
            //::InitializeCriticalSection(&mtx_);
        }

		void DoLock() {/*::EnterCriticalSection(&mtx_);*/}
		void DoUnLock() {/*::LeaveCriticalSection(&mtx_);*/}

         ~ObjectLevelLockable()
        {
            //::DeleteCriticalSection(&mtx_);
        }

        class Lock;
        friend class Lock;
        
        class Lock
        {
            ObjectLevelLockable const& host_;
            
            Lock(const Lock&);
            Lock& operator=(const Lock&);
        public:
            explicit Lock(const ObjectLevelLockable& host) : host_(host)
            {
                //::EnterCriticalSection(&host_.mtx_);
            }

            ~Lock()
            {
              //::LeaveCriticalSection(&host_.mtx_);
            }
        };

        typedef volatile Host VolatileType;

        typedef int IntType; 
/*
        static IntType AtomicIncrement(volatile IntType& lval)
        { return InterlockedIncrement(&const_cast<IntType&>(lval)); }
        
        static IntType AtomicDecrement(volatile IntType& lval)
        { return InterlockedDecrement(&const_cast<IntType&>(lval)); }

        static void AtomicAssign(volatile IntType& lval, IntType val)
        { InterlockedExchange(&const_cast<IntType&>(lval), val); }
        
        static void AtomicAssign(IntType& lval, volatile IntType& val)
        { InterlockedExchange(&lval, val); }
*/
        
        static IntType AtomicIncrement(volatile IntType& lval)
        { return ++lval; }
        
        static IntType AtomicDecrement(volatile IntType& lval)
        { return --lval; }
        
        static void AtomicAssign(volatile IntType & lval, IntType val)
        { lval = val; }
        
        static void AtomicAssign(IntType & lval, volatile IntType & val)
        { lval = val; }


    };

    template<typename Host>
    class ClassLevelLockable
    {
        struct Initializer
        {   
            //CRITICAL_SECTION mtx_;
            bool init_;

            Initializer():init_(false)
            {
                //::InitializeCriticalSection(&mtx_);
                init_=true;
            }
            ~Initializer()
            {
                assert(init_);
                //::DeleteCriticalSection(&mtx_);
            }
        };
        
        static Initializer initializer_;

    public:
		void DoLock()
		{
			assert(initializer_.init_);
			//::EnterCriticalSection(&initializer_.mtx_);
		}
		void DoUnLock()
		{
            assert(initializer_.init_);
            //::LeaveCriticalSection(&initializer_.mtx_);
		}

        class Lock;
        friend class Lock;
        
        class Lock
        {
            Lock(const Lock&);
            Lock& operator=(const Lock&);
        public:
            Lock()
            {
                assert(initializer_.init_);
                //::EnterCriticalSection(&initializer_.mtx_);
            }
            explicit Lock(const ClassLevelLockable&)
            {
                assert(initializer_.init_);
                //::EnterCriticalSection(&initializer_.mtx_);
            }
            ~Lock()
            {
                assert(initializer_.init_);
                //::LeaveCriticalSection(&initializer_.mtx_);
            }
        };

        typedef volatile Host VolatileType;

        typedef int IntType; 
/*
        static IntType AtomicIncrement(volatile IntType& lval)
        { return InterlockedIncrement(&const_cast<IntType&>(lval)); }
        
        static IntType AtomicDecrement(volatile IntType& lval)
        { return InterlockedDecrement(&const_cast<IntType&>(lval)); }
        
        static void AtomicAssign(volatile IntType& lval, IntType val)
        { InterlockedExchange(&const_cast<IntType&>(lval), val); }
        
        static void AtomicAssign(IntType& lval, volatile IntType& val)
        { InterlockedExchange(&lval, val); }
*/
        static IntType AtomicIncrement(volatile IntType& lval)
        { return ++lval; }
        
        static IntType AtomicDecrement(volatile IntType& lval)
        { return --lval; }
        
        static void AtomicAssign(volatile IntType & lval, IntType val)
        { lval = val; }
        
        static void AtomicAssign(IntType & lval, volatile IntType & val)
        { lval = val; }

    };

    template<typename Host>
    typename ClassLevelLockable<Host>::Initializer 
    ClassLevelLockable<Host>::initializer_;

#endif
}

////////////////////////////////////////////////////////////////////////////////
// Change log:
// June 20, 2001: ported by Nick Thurn to gcc 2.95.3. Kudos, Nick!!!
// July 26, 2005: some asserts by Peter K�mmel
////////////////////////////////////////////////////////////////////////////////


// $Log: Threads.h,v $
// Revision 3.51  2010/06/17 09:56:06  pcarvalho
// PwS v3.51 Pipelined 1.0
//
// Revision 1.13  2005/09/26 07:33:04  syntheticpp
// move macros into LOKI_ namespace
//
// Revision 1.12  2005/07/31 14:00:48  syntheticpp
// make object level threading possible
//
// Revision 1.11  2005/07/28 14:26:09  syntheticpp
// add cvs Header/Log
//
#endif
