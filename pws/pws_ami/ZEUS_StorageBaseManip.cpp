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

#include "ZEUS_StorageBaseManip.h"

//---------------------------------------------------------------------------

namespace Zeus
{

DBField sDBField2DBField(const sDBField& sDB)
{
	if(sDB.TypeId() == typeid(Loki::NullType)) return DBField();
	if(sDB.TypeId() == typeid(double))			return DBField(sDB.Get<double>());
	if(sDB.TypeId() == typeid(int))				return DBField(sDB.Get<int>());
	return DBField(*(sDB.GetPtr<std::wstring>()));
} 

} //End of namespace Zeus
