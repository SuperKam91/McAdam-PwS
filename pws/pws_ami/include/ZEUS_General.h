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



#ifndef ZEUS_GENERAL
#define ZEUS_GENERAL

#ifdef	WIN32
#include <process.h>
#else
#include <sys/types.h>
#include <unistd.h>
#endif

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <complex>
#include <math.h>


#include "ZEUS_Strings.h"

#define OBJEPS					1e-10

#define MINIMUMDUOCIMATION		3
//(2^n)

#define CONVDEC std::ios::dec | std::ios::skipws | std::ios::uppercase

namespace Zeus
{
class libException;

char*			Wstr2Achar(const std::wstring& wStr,char* ptCh,const std::locale* const loc=0,char def='?');
std::string		Wstr2Str(const std::wstring& wStr,const std::locale* const loc=0,char def='?');
wchar_t*		Achar2Awchar(const char *cstr,wchar_t* ptWch,const std::locale* const loc=0);
std::wstring	Achar2Wstr(const char *cstr,const std::locale* const loc=0);
std::wstring	LeftTrim(const std::wstring& str,std::ctype_base::mask myMask=std::ctype_base::space,const std::locale* const loc=0);
std::wstring	RightTrim(const std::wstring& str,std::ctype_base::mask myMask=std::ctype_base::space,const std::locale* const loc=0);
std::wstring	FullTrim(const std::wstring& str,std::ctype_base::mask myMask=std::ctype_base::space,const std::locale* const loc=0);
std::wstring	ToCase(const std::wstring& source,bool LowUpperCase,const std::locale* const loc=0);
wchar_t			ToCase(wchar_t c,bool LowUpperCase,const std::locale* const loc=0);
std::string		ToCase(const std::string& source,bool LowUpperCase,const std::locale* const loc=0);
char			ToCase(char c,bool LowUpperCase,const std::locale* const loc=0);
void			MySleep(unsigned int sec);

#if defined(WIN32) || defined(PS_USELOCALE)
template<typename T>
bool get_number(const std::wstring& txt,T& my_value,std::ios_base::fmtflags ConvBits=CONVDEC,const std::locale* const loc=0)
{

	std::basic_stringstream<wchar_t> temp;
	std::locale myStdLoc = std::locale();
	const std::locale*	tmpLoc = &myStdLoc;

	if(loc){temp.imbue(*loc);tmpLoc = loc;}
	std::ios_base::iostate	st(std::ios::goodbit);
	temp.setf(ConvBits,std::ios::basefield);
	std::use_facet<std::num_get<wchar_t,std::wstring::const_iterator> >((const std::locale&)(*tmpLoc)).get(txt.begin(),txt.end(),temp,st,my_value);

	return (st & std::ios_base::failbit) || (~st & std::ios_base::eofbit);
}

#else
template<typename T>
bool get_number(const std::wstring& txt,T& my_value,std::ios_base::fmtflags ConvBits=CONVDEC,const std::locale* const loc=0)
{
	std::wstring	tStr(Zeus::FullTrim(txt));
	const wchar_t 	* const str(tStr.c_str());
	const wchar_t 	* const endStr(str + tStr.size());

	wchar_t 		*endptr;
	my_value = static_cast<T>(wcstod(str, &endptr));
	return (endptr != endStr);	
}

bool get_number(const std::wstring& txt,long& my_value,std::ios_base::fmtflags ConvBits=CONVDEC,const std::locale* const loc=0);
#endif

template<typename T>
bool GetNumbersHomogenColl(const std::wstring& source,const std::wstring& seps,std::vector<T>& Coll,std::ios_base::fmtflags ConvBits=CONVDEC,const std::locale* const loc=0)
{

	T	tempValue;
	std::wstring::size_type	piv(0),pivAux;
	std::wstring::const_iterator strBeg(source.begin()),strEnd(source.end()),strPivBeg,strPivEnd;

	Coll.clear();

scan_again:	
	pivAux = source.find_first_not_of(seps,piv);
	if(pivAux == std::wstring::npos) return true;
	piv = source.find_first_of(seps,pivAux);
	if(piv == std::wstring::npos) {strPivEnd = strEnd;}
	else {strPivEnd = strBeg + piv;}
	strPivBeg = strBeg + pivAux;
	if(get_number(std::wstring(strPivBeg,strPivEnd),tempValue,ConvBits,loc)) return false;
	Coll.push_back(tempValue);
	if(piv == std::wstring::npos) return true;
	goto scan_again;
}

int		SplitSSString(const std::wstring& source,const std::wstring& seps,std::vector<std::wstring>& Coll);

inline	int	GetPID(void)
{
#ifdef WIN32
	return _getpid();
#else
	return getpid();
#endif
}

} // end of namespace Zeus

#endif //ZEUS_GENERAL
