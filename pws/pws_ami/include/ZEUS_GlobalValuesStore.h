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
#ifndef ZEUS_GLOBALVALUESSTOREH
#define ZEUS_GLOBALVALUESSTOREH

#include "ZEUS_StorageBaseManip.h"
#include "ZEUS_SingletonSerialize.h"
#include "ZEUS_InOut.h"


namespace Zeus
{

struct VariableIdType
{
	std::wstring	Id_;
	Zeus::DBField	DefValue_;
	VariableIdType(const std::wstring& id): Id_(id),DefValue_()
	{}
	VariableIdType(const std::wstring& id,const Zeus::DBField& defvalue): Id_(id),DefValue_(defvalue)
	{}
};

class	GlobalValuesStore
{
public:

	typedef Zeus::DBField		DBFieldType;

	inline bool		PutVar(const std::wstring& id, const Zeus::DBField& value)
	{return Insert1element(id,value);}

	int				PutVar(const std::vector<std::wstring>& ids, const std::vector<Zeus::DBField>& values);
	inline bool		GetVar(const VariableIdType& id, Zeus::DBField& value) const
	{return Get1var(id,value);}

	bool			GetVar(const std::vector<VariableIdType>& ids, std::vector<Zeus::DBField>& values) const;
	inline bool		RemoveElement(const std::wstring& id)
	{return Remove1Elem(id);}
	inline bool		RemoveElement(const VariableIdType& id)
	{return Remove1Elem(id.Id_);}

	bool			AddDataFromFile(int ContextID,const std::wstring& fname);
	inline const ParamVarsStoreType::HeaderType	*GetHeader(void) const;
	inline void		ReleaseData(void)
	{
		DataStore_.Storage_.clear();
	}

	inline bool		InitializeFromFile(const std::wstring& fname,int ContextID = 0)
	{
		ReleaseData();
		return AddDataFromFile(ContextID,fname);
	}

	GlobalValuesStore(void)
	{}

	~GlobalValuesStore(void)
	{
		ReleaseData();	
	}

private:
	GlobalValuesStore(const GlobalValuesStore&);
	GlobalValuesStore& operator=(const GlobalValuesStore&);

	inline bool		Remove1Elem(const std::wstring& id)
	{return DataStore_.Storage_.erase(id) == 1;}
	bool			Get1var(VariableIdType id,DBFieldType& value) const;
	inline void		errGlobalValuesStore(int errCode,wchar_t* msg) const
	{throw Zeus::libException(errCode,msg,*this);}
	void			errGlobalValuesStore(int errCode,wchar_t* msg,std::wstring xtraName) const;
    inline bool		Insert1element(const std::wstring& id,const DBField& value)
	{
		DataStore_.Storage_.erase(id);
		return DataStore_.Storage_.insert(ParamVarsStoreType::StorageType::value_type(id, value)).second;
	}

	ParamVarsStoreType	DataStore_;
};

typedef Zeus::Multi_SingletonHolder<GlobalValuesStore>::Type		MThGlobalVars;
typedef Zeus::Single_SingletonHolder<GlobalValuesStore>::Type		SThGlobalVars;

GenCollReader<ParamVarsStoreType>
*GetParamVarsStoreReaderHandler(Loki::Type2Type<ParamVarsStoreType>,int ContextID,const std::wstring& fn);


} // end of namespace Zeus
#endif //ZEUS_GLOBALVALUESSTOREH
