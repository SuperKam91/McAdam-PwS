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


//------------------------------------------
#ifndef ZEUSINOUTHEALPIXFITS
#define ZEUSINOUTHEALPIXFITS

#include "fitshandle.h"
#include "ZEUS_InOut.h"

//------------------------------------------

#define	HDNUMDEFAULT				2
#define HEALPIX_TEMP_COLUMN			1


namespace Zeus
{
template<typename T>
class	ReadHealpixMapFits: public GenHealpixReader<T>
{
protected:
	virtual std::wstring	do_GetCollID(void) const
	{return FName_;}
	virtual bool	do_Initialize(Zeus::HealpixHeaderType * hd)
	{
		std::wcout << FName_;
		inp_.open(Zeus::Wstr2Str(FName_));
		inp_.goto_hdu(HDNUMDEFAULT);
		
		if(!hd) return true;

		hd->NSide_ = Healpix_Base::npix2nside((int)inp_.nelems(HEALPIX_TEMP_COLUMN));

		if(hd->Columns_.empty())
			return true;

		HealpixHeaderType::HealpixHeaderItemCollType::iterator			piv(hd->Columns_.begin());
		HealpixHeaderType::HealpixHeaderItemCollType::const_iterator	const end(hd->Columns_.end());

		for(;piv!=end;++piv)
		{
			if(!inp_.key_present(piv->Keyword_))
			{piv->Reset();continue;}
			switch(piv->ValueType_)
			{
			case Zeus::INT_TYPE:
				inp_.get_key(piv->Keyword_,piv->PODValue_.intType_);
				break;
			case Zeus::FLOAT_TYPE:
				inp_.get_key(piv->Keyword_,piv->PODValue_.floatType_);
				break;
			case Zeus::DOUBLE_TYPE:
				inp_.get_key(piv->Keyword_,piv->PODValue_.doubleType_);
				break;
			case Zeus::BOOL_TYPE:
				inp_.get_key(piv->Keyword_,piv->PODValue_.boolType_);
				break;
			case Zeus::STRING_TYPE:
				inp_.get_key(piv->Keyword_,piv->stringValue_);
				break;
			default:
				piv->Reset();
				break;
			}
		}
		return true;
	}

	virtual void	do_Release(void)
	{
		inp_.close();
	}
	virtual bool	do_Read(void)
	{
		read_Healpix_map_from_fits<T>(inp_,*(GenHealpixReader<T>::HPixMap_),HEALPIX_TEMP_COLUMN);
		return true;
	}

public:
	ReadHealpixMapFits(const std::wstring& fname)
		:GenHealpixReader<T>(),FName_(fname + std::wstring(L".fits"))
	{}
	virtual ~ReadHealpixMapFits(void)
	{}
private:
	std::wstring		FName_;
	fitshandle			inp_;
};
//
template<typename T>
class	WriteHealpixMapFits: public GenHealpixWriter<T>
{
protected:
	virtual  std::wstring	do_GetCollID(void) const
	{return FName_;}
//
	virtual bool	do_Initialize(void)
	{return true;}
//
	virtual int	do_Write(const Healpix_Map<T>& data,HealpixHeaderType * hd)
	{
		fitshandle	out;

		out.create(Wstr2Str(FName_));

		if(hd)
		{CreateHeader(*hd,out);}

		write_Healpix_map_to_fits(out,data,planckType<T>());
		return true;
	}
//
public:
	WriteHealpixMapFits(const std::wstring& fname,Healpix_Ordering_Scheme Scheme,coordsys CoordSys)
		:GenHealpixWriter<T>(),FName_(std::wstring(L"!") + fname + std::wstring(L".fits")),Scheme_(Scheme),CoordSys_(CoordSys)
	{}
	virtual ~WriteHealpixMapFits(void)
	{}
private:
	std::wstring				FName_;
	Healpix_Ordering_Scheme		Scheme_;
	coordsys					CoordSys_;
//
	void	CreateHeader(const HealpixHeaderType& hd,fitshandle& out)
	{
		out.set_key<std::string>(std::string(HFIDMC_HEALPIXMAP_SCHEME),GetStrFromOrdering(Scheme_),std::string(""));
		out.set_key<std::string>(std::string(HFIDMC_HEALPIXMAP_COORSYS),GetStrFromCoordsys(CoordSys_),std::string(""));

		if(!(hd.Columns_.empty()))
		{
			HealpixHeaderType::HealpixHeaderItemCollType::const_iterator	piv(hd.Columns_.begin());
			HealpixHeaderType::HealpixHeaderItemCollType::const_iterator	const end(hd.Columns_.end());
			for(;piv != end;++piv)
			{
				switch(piv->ValueType_)
				{
				case INT_TYPE:
					out.set_key<int>(piv->Keyword_,piv->PODValue_.intType_,std::string(""));
					break;
				case FLOAT_TYPE:
					out.set_key<float>(piv->Keyword_,piv->PODValue_.floatType_,std::string(""));
					break;
				case DOUBLE_TYPE:
					out.set_key<double>(piv->Keyword_,piv->PODValue_.doubleType_,std::string(""));
					break;
				case BOOL_TYPE:
					out.set_key<bool>(piv->Keyword_,piv->PODValue_.boolType_,std::string(""));
					break;
				case STRING_TYPE:
					out.set_key<std::string>(piv->Keyword_,piv->stringValue_,std::string(""));
					break;
				default:
					break;
				}
			}
		}
	}
};

}

#endif //ZEUSINOUTHEALPIXFITS

