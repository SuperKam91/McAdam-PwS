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
#ifndef FFTW_EXCEPTIONSH
#define FFTW_EXCEPTIONSH
//---------------------------------------------------------------------------

#include <typeinfo>
#include <exception>
#include "ZEUS_Strings.h"
#include "ZEUS_General.h"

namespace Zeus
{

class libException : public std::exception
{
private:
	int			  errCode_;
	std::wstring  errMsg_;
	std::wstring  errN2_text(int errCode)
	{
	  wchar_t buffer[7];
	  //if(errCode < 0 || errCode >= 1000000) errCode = 0;
#ifdef WIN32
		_snwprintf_s(buffer,7,L"%d",errCode);
#else
		swprintf(buffer,7,L"%d",errCode);
#endif
	  return buffer;
	}
  void buildErrMsg(int ErrorCode ,const std::wstring& tName,const std::wstring& tyName)
{
	errMsg_ =			(std::wstring(ERROR_INTRO_STRING)		+
						tName					+
						ERROR_NUMBER_CODE		+
						errN2_text(ErrorCode)	+ 
						ERROR_TYPE_INFO			+
						tyName);
}

public:
	template<typename T>
	inline libException(int ErrorCode ,const std::wstring& tName,const T& objThis)
	{errCode_ = ErrorCode;buildErrMsg(ErrorCode,tName,Achar2Wstr(typeid(objThis).name()));}
//
	inline libException(int ErrorCode ,const std::wstring& tName,const std::wstring& functName)
	{errCode_ = ErrorCode;buildErrMsg(ErrorCode,tName,functName);}
//
	inline libException(int ErrorCode ,const wchar_t* tName,const wchar_t* functName)
	{errCode_ = ErrorCode;buildErrMsg(ErrorCode,std::wstring(tName),std::wstring(functName));}
//
	inline int	GetErr(void) const
	{return errCode_;}
//
	virtual const char*  what () const throw()
	{
		return Wstr2Str(errMsg_).c_str();
	}
//
	virtual std::wstring  what_Xmsg () const throw()
	{
		return errMsg_;
	}
//
	virtual ~libException(void) throw() {}
};
//
class libMPIException : public libException
{
public:
	inline libMPIException(int ErrorCode ,const std::wstring& tName,const std::wstring& functName)
		:libException(ErrorCode,tName,functName)
	{}
//
	inline libMPIException(int ErrorCode ,const wchar_t* tName,const wchar_t* functName)
		:libException(ErrorCode,tName,functName)
	{}
//
};

} // end of namesapce Zeus

#endif //FFTW_EXCEPTIONSH
