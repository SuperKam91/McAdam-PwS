////////////////////////////////////////////////////////////////////////////////
// The Loki Library
// Copyright (c) 2001 by Andrei Alexandrescu
// Copyright (c) 2005 by Peter Kuemmel
// This code DOES NOT accompany the book:
// Alexandrescu, Andrei. "Modern C++ Design: Generic Programming and Design 
//     Patterns Applied". Copyright (c) 2001. Addison-Wesley.
//
// Code covered by the MIT License
// The authors make no representations about the suitability of this software
// for any purpose. It is provided "as is" without express or implied warranty.
////////////////////////////////////////////////////////////////////////////////

// $Header: /planck/cvs/planck/Level2/Task_pkg/HL2_PowellSnakes/src/include/Factory.h,v 3.51 2010/06/17 09:56:05 pcarvalho Exp $

#ifndef LOKI_FACTORYPARM_INC_
#define LOKI_FACTORYPARM_INC_

#include <map>

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable: 4702)
//unreachable code if OnUnknownType throws an exception
#endif


namespace Loki
{

////////////////////////////////////////////////////////////////////////////////
// class template DefaultFactoryError
// Manages the "Unknown Type" error in an object factory
////////////////////////////////////////////////////////////////////////////////

template <typename IdentifierType, class AbstractProduct>
struct DefaultFactoryError
{
    struct Exception : public std::exception
    {
        const char* what() const throw() { return "Unknown Type"; }
    };
    
    static AbstractProduct* OnUnknownType(IdentifierType)
    {
        throw Exception();
    }
};
    
template
<
    class AbstractProduct, 
    typename IdentifierType,
    typename ProductCreator = AbstractProduct* (*)(const void*),
    template<typename, class>
        class FactoryErrorPolicy = DefaultFactoryError
>
class Factory 
    : public FactoryErrorPolicy<IdentifierType, AbstractProduct>
{
public:
    bool Register(const IdentifierType& id, ProductCreator creator)
    {
        return associations_.insert(
            typename IdToProductMap::value_type(id, creator)).second;
    }
    
    bool Unregister(const IdentifierType& id)
    {
        return associations_.erase(id) == 1;
    }
    
    AbstractProduct* CreateObject(const IdentifierType& id,const void* args=0)
    {
        typename IdToProductMap::iterator i = associations_.find(id);
        if (i != associations_.end())
        {
			AbstractProduct* t((i->second)(id,args));
			t->Initialize();
            return t;
        }
        return this->OnUnknownType(id);
    }
    
private:
	typedef std::map<IdentifierType, ProductCreator> IdToProductMap;
    IdToProductMap associations_;
};

} // namespace Loki

#ifdef _MSC_VER
#pragma warning( pop ) 
#endif

////////////////////////////////////////////////////////////////////////////////
// Change log:
// June 20,    2001: ported by Nick Thurn to gcc 2.95.3. Kudos, Nick!!!
// May 08,     2002: replaced const_iterator with iterator so that self-modifying
//                   ProductCreators are supported. Also, added a throw()
//                   spec to what(). Credit due to Jason Fischl.
// February 2, 2003: fixed dependent names - credit due to Rani Sharoni
// March 4,    2003: fixed dependent names - credit due to Ruslan Zasukhin and CW 8.3 
// July 26,    2005: parameter support by Peter Kümmel 
////////////////////////////////////////////////////////////////////////////////

#endif // FACTORY_INC_

// $Log: Factory.h,v $
// Revision 3.51  2010/06/17 09:56:05  pcarvalho
// PwS v3.51 Pipelined 1.0
//
// Revision 1.7  2005/10/05 09:57:37  syntheticpp
// move unreachable code warnings
//
// Revision 1.6  2005/09/26 07:33:04  syntheticpp
// move macros into LOKI_ namespace
//
// Revision 1.5  2005/07/31 14:23:24  syntheticpp
// invert new factory code macro logic to be ReferenceTest more compatible with noncc code
//
// Revision 1.4  2005/07/28 14:26:09  syntheticpp
// add cvs Header/Log
//

