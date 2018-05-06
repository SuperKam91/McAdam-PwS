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
#ifndef ZEUS_FACTORYH
#define ZEUS_FACTORYH

#include "Factory.h"
#include "ZEUS_Strings.h"
#include "ZEUS_StorageBaseManip.h"
#include "ZEUS_SingletonSerialize.h"
#include "ZEUS_Exceptions.h"
//---------------------------------------------------------------------------

namespace Zeus{

template<typename IdentifierType, typename AbstractProduct>
struct ZEUS_FactoryError
{
    struct tmp_Exception : public libException
    {
		tmp_Exception(const IdentifierType& Id) : libException(ERRORMSGZEUS(UNKNONOBJ),Id)
		{}
    };
    
    static AbstractProduct* OnUnknownType(const IdentifierType& Id)
    {
        throw tmp_Exception(Id);
    }
};

typedef DeviceBaseType											FactoryObjType;
typedef FactoryObjType* (*FactoryCreatorType)(const std::wstring&,const void*);
typedef Loki::Factory<FactoryObjType,std::wstring,FactoryCreatorType,ZEUS_FactoryError> FactoryType;

} // end of namespace Zeus

#endif //ZEUS_FACTORYH
