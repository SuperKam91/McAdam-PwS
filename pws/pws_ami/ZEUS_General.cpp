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

#ifdef	WIN32
#include <Windows.h>
#endif

#include <string.h>
#include <stdio.h>

#include "ZEUS_General.h"

//---------------------------------------------------------------------------

namespace Zeus
{

void			MySleep(unsigned int sec)
{
#ifdef WIN32
	Sleep((DWORD)(1000*sec));
#else
	sleep(sec);
#endif
}

int		SplitSSString(const std::wstring& source,const std::wstring& seps,std::vector<std::wstring>& Coll)
{
	std::wstring::size_type	piv(0),pivAux;
	std::wstring::const_iterator strBeg(source.begin()),strEnd(source.end()),strPivBeg,strPivEnd;

	Coll.clear();

scan_again:	
	pivAux = source.find_first_not_of(seps,piv);
	if(pivAux == std::wstring::npos) return Coll.size();
	piv = source.find_first_of(seps,pivAux);
	if(piv == std::wstring::npos) {strPivEnd = strEnd;}
	else {strPivEnd = strBeg + piv;}
	strPivBeg = strBeg + pivAux;	
	Coll.push_back(std::wstring(strPivBeg,strPivEnd));
	if(piv == std::wstring::npos) return Coll.size();
	goto scan_again;
}


char*	Wstr2Achar(const std::wstring& wStr,char* ptCh,const std::locale* const loc,char def)

{
	std::size_t wSz(wStr.size());

	memset(ptCh,0,(wSz+1) * sizeof(char));
	std::use_facet<std::ctype<wchar_t> >((const std::locale&)(loc ? *loc : std::locale())).narrow(wStr.c_str(),wStr.c_str() + wSz,def,ptCh);
	return ptCh;
}

std::string Wstr2Str(const std::wstring& wStr,const std::locale* const loc,char def)

{
  std::size_t wSz(wStr.size());

  char *ptCh(new char[wSz+1]);
  memset(ptCh,0,(wSz+1) * sizeof(char));
  std::use_facet<std::ctype<wchar_t> >((const std::locale&)(loc ? *loc : std::locale())).narrow(wStr.c_str(),wStr.c_str() + wSz,def,ptCh);
  std::string tmp(ptCh,ptCh + wSz);
  delete [] ptCh;
  return tmp;
}


wchar_t* Achar2Awchar(const char *cstr,wchar_t* ptWch,const std::locale* const loc)

{
  std::size_t cSz(strlen(cstr));
  memset(ptWch,0,(cSz+1) * sizeof(wchar_t));
  std::use_facet<std::ctype<wchar_t> >((const std::locale&)(loc ? *loc : std::locale())).widen(cstr,cstr + cSz,ptWch);
  return ptWch;
}

std::wstring Achar2Wstr(const char *cstr,const std::locale* const loc)

{
  std::size_t cSz(strlen(cstr));

  wchar_t* ptWch(new wchar_t[cSz+1]);
  memset(ptWch,0,(cSz+1) * sizeof(wchar_t));
  std::use_facet<std::ctype<wchar_t> >((const std::locale&)(loc ? *loc : std::locale())).widen(cstr,cstr + cSz,ptWch);
  std::wstring tmp(ptWch,ptWch + cSz);
  delete [] ptWch;
  return tmp;
}

std::wstring LeftTrim(const std::wstring& str,std::ctype_base::mask myMask,const std::locale* const loc)

{
	if(str.empty()) return std::wstring();
	const std::ctype<wchar_t>& ct = std::use_facet<std::ctype<wchar_t> >((const std::locale&)(loc ? *loc : std::locale()));

	const wchar_t* tptrBeg(str.c_str());
	const wchar_t* tptrEnd(tptrBeg + str.size());
	if((tptrBeg = ct.scan_not(myMask,tptrBeg,tptrEnd)) == tptrEnd) return std::wstring();
	return std::wstring(tptrBeg,tptrEnd);
}

std::wstring RightTrim(const std::wstring& str,std::ctype_base::mask myMask,const std::locale* const loc)

{
  if(str.empty()) return std::wstring();
  const std::ctype<wchar_t>& ct = std::use_facet<std::ctype<wchar_t> >((const std::locale&)(loc ? *loc : std::locale()));
  std::wstring::const_iterator tptrRev(str.end());
  while(ct.is(myMask,*--tptrRev)) ;
  return std::wstring(str.begin(),++tptrRev);
}

std::wstring FullTrim(const std::wstring& str,std::ctype_base::mask myMask,const std::locale* const loc)

{
	if(str.empty()) return std::wstring();
	const std::ctype<wchar_t>& ct = std::use_facet<std::ctype<wchar_t> >((const std::locale&)(loc ? *loc : std::locale()));
  
	const wchar_t* tptrBeg(str.c_str());
	const wchar_t* tptrEnd(tptrBeg + str.size());

	if((tptrBeg = ct.scan_not(myMask,tptrBeg,tptrEnd)) == tptrEnd) return std::wstring();
	while(ct.is(myMask,*--tptrEnd)) ;
	return std::wstring(tptrBeg,++tptrEnd);
}

std::wstring ToCase(const std::wstring& source,bool LowUpperCase,const std::locale* const loc)
{
	if(source.empty()) return std::wstring();

	std::wstring srcTemp(source);
	wchar_t* tptrBeg(const_cast<wchar_t*>(srcTemp.c_str()));
	wchar_t* tptrEnd(tptrBeg + srcTemp.size());

	const std::ctype<wchar_t>& conv(std::use_facet<std::ctype<wchar_t> >((const std::locale&)(loc ? *loc : std::locale())));

	if(LowUpperCase) {conv.tolower(tptrBeg,tptrEnd);}
	else {conv.toupper(tptrBeg,tptrEnd);}
	return std::wstring(tptrBeg);
}

std::string ToCase(const std::string& source,bool LowUpperCase,const std::locale* const loc)
{
	if(source.empty()) return std::string();

	std::string srcTemp(source);
	char* tptrBeg(const_cast<char*>(srcTemp.c_str()));
	char* tptrEnd(tptrBeg + srcTemp.size());

	const std::ctype<char>& conv(std::use_facet<std::ctype<char> >((const std::locale&)(loc ? *loc : std::locale())));

	if(LowUpperCase) {conv.tolower(tptrBeg,tptrEnd);}
	else {conv.toupper(tptrBeg,tptrEnd);}
	return std::string(tptrBeg);
}

wchar_t	ToCase(wchar_t c,bool LowUpperCase,const std::locale* const loc)
{
	const std::ctype<wchar_t>& conv(std::use_facet<std::ctype<wchar_t> >((const std::locale&)(loc ? *loc : std::locale())));
	return (LowUpperCase?conv.tolower(c):conv.toupper(c));
}

char	ToCase(char c,bool LowUpperCase,const std::locale* const loc)
{
	const std::ctype<char>& conv(std::use_facet<std::ctype<char> >((const std::locale&)(loc ? *loc : std::locale())));
	return (LowUpperCase?conv.tolower(c):conv.toupper(c));
}

bool get_number(const std::wstring& txt,long& my_value,std::ios_base::fmtflags ConvBits,const std::locale* const loc)
{
	std::wstring	tStr(Zeus::FullTrim(txt));
	const wchar_t 	* const str(tStr.c_str());
	const wchar_t 	* const endStr(str + tStr.size());

	wchar_t 		*endptr;
	my_value = wcstol(str, &endptr,10);
	return (endptr != endStr);	
}
//
} //end of namespace Zeus
