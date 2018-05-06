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

//----------------------------


#include <memory>
#include "ZEUS_Exceptions.h"
#include "ZEUS_InOut.h"
#include "ZEUS_InOutTxtFile.h"

#if defined(HFIDMC) || defined(LFIDPC)
#include "ZEUS_InOutPipeline.h"
#endif

#include "ZEUS_GlobalValuesStore.h"

//----------------------------

namespace Zeus
{

int			GlobalValuesStore::PutVar(const std::vector<std::wstring>& ids,const std::vector<DBFieldType>& values)
{
	int NInsert(0);
	std::vector<std::wstring>::const_iterator		pivId(ids.begin());
	std::vector<std::wstring>::const_iterator		const endId(ids.end());
	std::vector<Zeus::DBField>::const_iterator		pivVal(values.begin());		
	std::vector<Zeus::DBField>::const_iterator		const endVal(values.end());

	for(;(pivId != endId) && (pivVal != endVal);++pivId,++pivVal)
	{
		if(!Insert1element(*pivId,*pivVal)) break;
		++NInsert;
	}
	return NInsert;
}

bool		GlobalValuesStore::GetVar(const std::vector<VariableIdType>& ids,std::vector<DBFieldType>& values) const
{
	std::vector<VariableIdType>::const_iterator		pivId(ids.begin());
	std::vector<VariableIdType>::const_iterator		const endId(ids.end());
	Zeus::DBField tempValue;

	values.clear();
	for(;pivId != endId;++pivId)
	{
		if(!Get1var(*pivId,tempValue))
		{errGlobalValuesStore(ERROR_COD_ZEUSGLBPARAMRD,ERROR_MSG_ZEUSGLBPARAMRD,(*pivId).Id_);}
		values.push_back(tempValue);
	}
	return (values.size() == ids.size());
}

bool		GlobalValuesStore::AddDataFromFile(int ContextID,const std::wstring& fname)
{
	std::auto_ptr<GenCollReader<ParamVarsStoreType> >
		FReader(GetParamVarsStoreReaderHandler(Loki::Type2Type<ParamVarsStoreType>(),ContextID,fname));
	FReader->Initialize();
	FReader->Read();
	FReader->Release(DataStore_);
	return true;
}

bool		GlobalValuesStore::Get1var(VariableIdType id,DBFieldType& value) const
{
	ParamVarsStoreType::StorageType::const_iterator i;
		
	do{
		i	= DataStore_.Storage_.find(id.Id_);
		if (i == DataStore_.Storage_.end())
		{
			if(id.DefValue_.TypeId() == typeid(Loki::NullType))
				return false;
			value = id.DefValue_;
		}
		else
		{value = i->second;}

		if(value.TypeId() == typeid(Zeus::ps_VariableType))
		{
			id.Id_			= (value.GetPtr<Zeus::ps_VariableType>())->VarId_;
			id.DefValue_	= Zeus::sDBField2DBField((value.GetPtr<Zeus::ps_VariableType>())->DefValue_);
		}
		else
		{break;}

	} while(true);
	
	return true;
}

void		GlobalValuesStore::errGlobalValuesStore(int errCode,wchar_t* msg,std::wstring xtraName) const
{
	std::wstring errstring(msg);

	errstring += std::wstring(L" -> ");
	errstring += xtraName;
	throw Zeus::libException(errCode,errstring,*this);
}


GenCollReader<ParamVarsStoreType>
*GetParamVarsStoreReaderHandler(Loki::Type2Type<ParamVarsStoreType>,int ContextID,const std::wstring& fn)
{
	switch(ContextID)
	{
	case 1:
#ifdef HFIDMC
	return new HFIDMC_ParameterFileReader(fn);
#else
	break;
#endif
	case 2:
#ifdef LFIDPC
	return new LFIDPC_ParameterFileReader(fn);
#else
	break;
#endif
	case 3:
		// Parameter files are txt files always
	break;
	default:
	break;	
	}
	return new ParameterFileReader(fn);
}


//--------------

} // End of Namespace Zeus

