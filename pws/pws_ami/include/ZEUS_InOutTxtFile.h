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

//------------------------------------
#ifndef ZEUSINOUTTXTFILEH
#define ZEUSINOUTTXTFILEH

#include "ZEUS_InOut.h"
#include "ZEUS_WorkSpace.h"
//------------------------------------

#define PEAKS_NFIELDS			52
#define PEAKSINFILELNFORMAT		L"%d,%d,%d,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%ls %ls"
#define PEAKSOUTFILELNFORMAT	L"%05d,%05d,%2d,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%3.5f,%3.5f,%3.5f,%3.5f,%3.5f,%3.5f,%3.5f,%3.5f,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%3.5f,%3.5f,%3.5f,%3.5f,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%ls %ls"

#define TXTOUTINVALIDFORMAT		L"%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g"

#define ERCSCCATALOGUEFORMAT	L"%3.4f %3.4f %7.6g %7.6g %7.6g %d"
//
#define GEOM_HEADER_NFIELDS		7
#define GEOM_HEADER_SCANF		L"%d,%d,%d,%d,%d,%d,%le"

#define GEOM_LINES_NFIELDS		38
#define GEOM_LINES_SCANF		L"%d,%d,%d,%d,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%d,%le,%le,%le,%le,%le,%le,%d,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le"

#define GEOM_HEADER_PRINTF		L"%4d,%4d,%4d,%4d,%4d,%4d,%7.6g"
#define GEOM_LINES_PRINTF		L"%5d,%2d,%5d,%5d,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%2d,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%2d,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g,%7.6g"
#define ASSIGNCHAR				L'='
#define COMMENTCHAR				L'#'

#define DETECTIONCAT_NFIELDS		8
#define DETECTIONCATNFORMAT			L"%d %le %le %le %le %le %le %le"

namespace Zeus
{
//
class TxtFileReaderRaw
{
protected:
	virtual int					do_parseLn(const std::wstring& str) = 0;
	inline void					errTxtFileReader(int errCode,wchar_t* msg,std::wstring xtraName) const
	{
		std::wstring errstring(msg);

		errstring += std::wstring(L" -> ");
		errstring += xtraName;
		throw Zeus::libException(errCode,errstring,*this);
	}
	inline void					errTxtFileReader(int errCode,wchar_t* msg) const
	{throw Zeus::libException(errCode,msg,*this);}
	inline const std::wstring&	GetCurrStrRef(void) const
	{return CurrStr_;}
	inline std::wstring&		UnreadLn(void)
	{LnOnBuffer_ = 1; return CurrStr_;}

	inline const std::locale*	GetLoc(void) const
	{return loc_;}
	bool						ReadNextLn(void);
	inline std::wstring			GetFileName(void) const
	{return fname_;}
public:
	TxtFileReaderRaw(const std::wstring& fname,bool NoOpen = false,const std::locale* loc = 0);
	TxtFileReaderRaw(std::wistream* istream,const std::locale* loc = 0);
	void						Open(void);
//
	inline int					ParseLn(void)
	{
		int res(0);
		while(!res){
			if(!ReadNextLn()) return 0;
			res = do_parseLn(CurrStr_);
			if(res > 0) ++parsedLns_;
		}
		return res;
	}
//
	inline void		ReadAllFile(void)
	{
		int stat;
		while(stat = ParseLn())
		{
			if(stat <0)
			{
				errTxtFileReader(ERROR_COD_ZEUSERRPARSLN,ERROR_MSG_ZEUSERRPARSLN,
					std::wstring(L" File Name -> ") + fname_ + std::wstring(L", Line -> ") + GetBufferStr());
			}
		}
	}
//
	inline void					GetCurrLnsInfo(int& ReadLn,int& ParsedLns) const
	{
		ReadLn = readLns_;ParsedLns = parsedLns_;
	}
//
	inline const std::wstring	GetBufferStr(void) const
	{return CurrStr_;}
//	
	inline	void				Reset(void)
	{
		istrm_->clear();
		istrm_->exceptions(std::ios::failbit | std::ios::badbit);
		istrm_->seekg(0);
		return;
	}

	virtual						~TxtFileReaderRaw(void)
	{
		if(!(fname_.empty()) && (istrm_)){
			delete istrm_;
			return;
		}
		if(istrm_){
			istrm_->clear();
			istrm_->exceptions(oldIoState_);
		}
	}

private:
	TxtFileReaderRaw(const TxtFileReaderRaw& rhs);
	TxtFileReaderRaw& operator= (const TxtFileReaderRaw& rhs);

	std::wstring				CurrStr_;
	const std::locale*			loc_;
	std::wstring				fname_;
	std::wistream*				istrm_;
	std::ios::iostate			oldIoState_;
	int							readLns_;
	int							parsedLns_;
	int							LnOnBuffer_;

};

//
template<class T>
class ReadTxtFileHelper: public GenCollReader<T>,protected TxtFileReaderRaw
{
protected:
	virtual std::wstring	do_GetCollID(void) const
	{return GetFileName();}
	virtual bool	do_Read(void)
	{ReadAllFile();return true;}

public:
	ReadTxtFileHelper(const std::wstring& fname)
		:Zeus::GenCollReader<T>(),Zeus::TxtFileReaderRaw(fname)
	{}
	virtual ~ReadTxtFileHelper(void)
	{}
};
//
class	TxtCollBaseWriter
{
public:
	TxtCollBaseWriter(const std::wstring& fname,std::ios_base::openmode mode)
		:fname_(fname),mode_(mode)
	{}
	void	Initialize(void);
//
	inline std::wstring GetFileName(void) const
	{return fname_;}
//
	inline  void WriteLn(const std::wstring& ln)
	{
		if(!(stream_.is_open()))
			return;
		stream_ << ln;std::endl(stream_);
	}
//
	inline	void flush(void)
	{if(stream_.is_open()) stream_.flush();}
//
	inline	void Close(void)
	{if(stream_.is_open()) stream_.close();}

//
private:
	TxtCollBaseWriter(const TxtCollBaseWriter& rhs);
	TxtCollBaseWriter& operator= (const TxtCollBaseWriter& rhs);
//
	std::ios_base::openmode		mode_;
	const std::wstring			fname_;
	std::wofstream				stream_;
};

namespace Private
{
	template<typename T>
	struct CheckComplex
	{
		typedef T	CollType;
		static T Parse(typename std::vector<CollType>::const_iterator& piv)
		{return static_cast<T>(*piv);}

	};
	template<typename T>
	struct CheckComplex<std::complex<T> >
	{
		typedef T	CollType;
		static std::complex<CollType> Parse(typename std::vector<CollType>::const_iterator& piv)
		{
			CollType real(*piv);
			++piv;
			return std::complex<CollType>(real,*piv);
		}
	};
}

template<typename T>
class	ReadWrkSpaceTxtFile:  public ReadTxtFileHelper<Zeus::LArr2D<T> >
{
private:
	typedef Zeus::LArr2D<T>									TType;
	typedef typename Private::CheckComplex<T>::CollType		CollType;

	TType	*data_;
	int		YSz_;
	int		Metric_;
	int		CurrLine_;
protected:
//
	virtual	TType&	do_GetData(void)
	{
		return *data_;
	}
//
	virtual void	do_DisposeData(void)
	{delete data_;data_ = 0;}
//
	virtual int		do_parseLn(const std::wstring& str)
	{
		std::vector<CollType> Coll;
		if(!Zeus::GetNumbersHomogenColl(str,std::wstring(L","),Coll)) return -1;
		typename std::vector<CollType>::const_iterator	const beg(Coll.begin());
		typename std::vector<CollType>::const_iterator	piv(beg);
		typename std::vector<CollType>::const_iterator	const end(Coll.end());
		for(unsigned int i=0;(piv != end) && (i!= Metric_);++piv,++i)
		{data_->operator ()(CurrLine_,i) = Private::CheckComplex<T>::Parse(piv) ;}
		++CurrLine_;
		return static_cast<int>(Coll.size());
	}
//
	virtual typename TType::HeaderType * 	do_Initialize(void)
	{
		do_DisposeData();
		TxtFileReaderRaw::Open();
		data_ = new TType(YSz_ * Metric_ ,Metric_);
		data_->reset();
		return reinterpret_cast<typename TType::HeaderType *>(data_);
	} 
//
public:
	ReadWrkSpaceTxtFile(const std::wstring& fname,int YSz,int XSz)
		: ReadTxtFileHelper<Zeus::LArr2D<T> >(fname   + std::wstring(ICAT_FILEEXT)),YSz_(YSz),Metric_(XSz),data_(0),CurrLine_(0)
	{}
//
	virtual ~ReadWrkSpaceTxtFile(void)
	{delete data_;}
};

template<typename T>
class	WriteWrkSpaceTxtFile: public GenCollWriter<Zeus::LArr2D<T> >
{
	typedef Zeus::LArr2D<T>	AtomT;
protected:
//
	inline virtual std::wstring	do_GetCollID(void) const
	{return stream_.GetFileName();}
//
	inline virtual bool	do_Initialize(void)
	{stream_.Initialize();return true;}
//
	inline virtual int	do_Write(const AtomT & data)
	{WorkspaceDump(data);stream_.flush(); return true;}
//
	inline virtual int			do_Flush(void)
	{stream_.flush();return 0;}
//
	inline virtual void			do_DisposeData(void)
	{stream_.Close();return ;}
//
	inline virtual int		do_Remove(void)
	{return std::remove(Wstr2Str(stream_.GetFileName()).c_str());}

public:
	WriteWrkSpaceTxtFile(const std::wstring& fname)
		:Zeus::GenCollWriter<Zeus::LArr2D<T> >(),stream_(fname + std::wstring(ICAT_FILEEXT),std::ios::out)
	{}
	virtual ~WriteWrkSpaceTxtFile(void)
	{}
private:
	int WorkspaceDump(const AtomT & data)
	{
		int				LnSz(data.getPtrMetric()),i;
		const T			*const piv0(data.begin());
		const T			*piv(piv0);
		const T			*const end(piv0 + data.getSz());
		std::wstring	strBuffer;

		for(i=0;piv!=end;++piv)
		{
			strBuffer += Zeus::PutNumber2Txt(*piv);
			if(++i==LnSz)
			{
				i = 0;
				stream_.WriteLn(strBuffer);
				strBuffer.clear();
			;}
			else {strBuffer+= L",";}
		}
		return (int) (piv-piv0);
	}
//
	Zeus::TxtCollBaseWriter stream_;
};
//
class	ParameterFileReader:  public ReadTxtFileHelper<ParamVarsStoreType>
{
	typedef ParamVarsStoreType			TType;

	TType								*data_;

protected:
	virtual	TType&	do_GetData(void)
	{
		return *data_;
	}
	virtual void	do_DisposeData(void)
	{delete data_;data_ = 0;}
	virtual int		do_parseLn(const std::wstring& str);
	virtual TType::HeaderType *do_Initialize(void)
	{
		delete data_;data_ = 0;
		TxtFileReaderRaw::Open();
		data_ = new ParamVarsStoreType;
		return reinterpret_cast<TType::HeaderType *>(data_);
	} 

public:
	ParameterFileReader(const std::wstring& fname)
		: ReadTxtFileHelper<ParamVarsStoreType>(fname + std::wstring(PARAM_FILEEXT)),data_(0)
	{}
	virtual ~ParameterFileReader(void)
	{delete data_;data_ = 0;}
private:
	int				SplitStr(const std::wstring& in,std::wstring& fld,std::wstring& value);
	int				GetValue(const std::wstring& in,DBField& fld);
    inline bool		Insert1element(const std::wstring& id,const DBField& value)
	{
		data_->Storage_.erase(id);
		return data_->Storage_.insert(ParamVarsStoreType::StorageType::value_type(id, value)).second;
	}

};
//
class	ReadPatchesGeomInfoTxtFile: public ReadTxtFileHelper<PatchGeomType>
{
private:
	typedef PatchGeomType	T;

	T		data_;
	void	ReadHeader(void);
protected:
	virtual	T&		do_GetData(void)
	{return data_;}
	virtual void	do_DisposeData(void)
	{data_.Storage_.clear();}
	virtual int		do_parseLn(const std::wstring& str);
	virtual T::HeaderType *do_Initialize(void)
	{
		data_.Storage_.clear();
		TxtFileReaderRaw::Open();
		ReadHeader();
		return &data_.Header_;
	}

public:
	ReadPatchesGeomInfoTxtFile(const std::wstring& fname)
		: ReadTxtFileHelper<PatchGeomType>(fname + std::wstring(ICAT_FILEEXT))
	{}
	virtual ~ReadPatchesGeomInfoTxtFile(void)
	{}
};
//
class	WritePatchesGeomInfoTxtFile: public GenCollWriter<PatchGeomType>
{
	typedef PatchGeomType	T;
protected:
//
	inline virtual std::wstring	do_GetCollID(void) const
	{return stream_.GetFileName();}
//
	inline virtual bool	do_Initialize(void)
	{stream_.Initialize();return true;}
//
	inline virtual int	do_Write(const T & data)
	{
		WriteHeader(data.Header_);
		WriteAllLns(data.Storage_);
		stream_.flush();
		return 1;
	}
//
	inline virtual int			do_Flush(void)
	{stream_.flush();return 0;}
//
	inline virtual void			do_DisposeData(void)
	{stream_.Close();return ;}
//
	inline virtual int		do_Remove(void)
	{return std::remove(Wstr2Str(stream_.GetFileName()).c_str());}
//
public:
	WritePatchesGeomInfoTxtFile(const std::wstring& fname)
		:GenCollWriter<PatchGeomType>(),stream_(fname + std::wstring(ICAT_FILEEXT),std::ios::out)
	{}
	virtual ~WritePatchesGeomInfoTxtFile(void)
	{}
private:
	void	WriteHeader(const T::HeaderType& head);
	void	Write1Line(const T::StorageType::value_type & pt);
	void	WriteAllLns(const T::StorageType  & coll);

	Zeus::TxtCollBaseWriter stream_;
};

//
class	WriteObjResultsTxtFile: public GenCollWriter<PeakCollType>
{
	typedef PeakCollReadbleType	AtomT;
protected:
//
	inline virtual std::wstring	do_GetCollID(void) const
	{return stream_.GetFileName();}
//
	inline virtual bool	do_Initialize(void)
	{stream_.Initialize();return true;}
//
	inline virtual int	do_Write(const PeakCollType & data)
	{
		int res(WritePeaksColl(data));
		stream_.flush();
		return res;
	}
//
	inline virtual int			do_Flush(void)
	{stream_.flush();return 0;}
//
	inline virtual void			do_DisposeData(void)
	{stream_.Close();return ;}
//
	inline virtual int		do_Remove(void)
	{return std::remove(Wstr2Str(stream_.GetFileName()).c_str());}

public:
	WriteObjResultsTxtFile(const std::wstring& fname,int InitObjID)
		:GenCollWriter<PeakCollType>(),stream_(fname + std::wstring(ICAT_FILEEXT),std::ios::app),
		ObjectID_(InitObjID)
	{}
	virtual ~WriteObjResultsTxtFile(void)
	{}
private:
	int		WritePeaksColl(const PeakCollType& data);
	int		ObjectID_;

	Zeus::TxtCollBaseWriter stream_;
};
//
class	ReadObjResultsTxtFile: public ReadTxtFileHelper<PeakCollReadbleType>
{
protected:

	virtual	PeakCollReadbleType&	do_GetData(void)
	{return	data_;}

	virtual void	do_DisposeData(void)
	{data_.Storage_.clear();}

	virtual int		do_parseLn(const std::wstring& str);

	virtual PeakCollReadbleType::HeaderType *do_Initialize(void)
	{
		data_.Storage_.clear();
		TxtFileReaderRaw::Open();
		return reinterpret_cast<PeakCollReadbleType::HeaderType *>(&data_);
	}
public:
	ReadObjResultsTxtFile(const std::wstring& fname)
		: ReadTxtFileHelper<PeakCollReadbleType>(fname + std::wstring(ICAT_FILEEXT))
	{}
	virtual ~ReadObjResultsTxtFile(void)
	{}
private:
	PeakCollReadbleType	data_;
};
//
class	CatWriterTxtFile: public GenCollWriter<CatalogueFormatType>
{
protected:
//
	virtual std::wstring		do_GetCollID(void) const
	{return stream_.GetFileName();}
//
	inline virtual bool			do_Initialize(void)
	{stream_.Initialize();return true;}
//
	inline virtual int			do_Write(const CatalogueFormatType& data)
	{
		CreatHeader(data);
		return WriteCatLinesColl(data);
	}
//
	inline virtual int			do_Flush(void)
	{stream_.flush();return 0;}
//
	inline virtual void			do_DisposeData(void)
	{stream_.Close();return ;}
//
	inline virtual int			do_Remove(void)
	{return std::remove(Wstr2Str(stream_.GetFileName()).c_str());}
//
	inline void					WriteLn(const wchar_t* msg)
	{stream_.WriteLn(std::wstring(msg));}
//
private:
//
	void	CreatHeader(const CatalogueFormatType & data);
//
	int		WriteCatLinesColl(const CatalogueFormatType & data);
//
	TxtCollBaseWriter				stream_;
//
public:
	CatWriterTxtFile(const std::wstring& fname)
		:GenCollWriter<CatalogueFormatType>(),stream_(fname + std::wstring(FCAT_FILEEXT),std::ios::out)
	{}
	virtual ~CatWriterTxtFile(void)
	{}
};
//
class	CatWriterOutputTxtFile: public GenCollWriter<OutputFormatType>
{
protected:
//
	virtual std::wstring		do_GetCollID(void) const
	{return stream_.GetFileName();}
//
	inline virtual bool			do_Initialize(void)
	{stream_.Initialize();return true;}
//
	inline virtual int			do_Write(const OutputFormatType& data)
	{
		CreatHeader(data);
		return WriteCatLinesColl(data);
	}
//
	inline virtual int			do_Flush(void)
	{stream_.flush();return 0;}
//
	inline virtual void			do_DisposeData(void)
	{stream_.Close();return ;}
//
	inline virtual int			do_Remove(void)
	{return std::remove(Wstr2Str(stream_.GetFileName()).c_str());}
//
	inline void					WriteLn(const wchar_t* msg)
	{stream_.WriteLn(std::wstring(msg));}
//
	std::wstring				CreatHeaderHelper(const OutputFormatType& data);
//
	virtual	void				CreatHeader(const OutputFormatType& data) = 0;
//
	virtual	int					WriteCatLinesColl(const OutputFormatType& data) = 0;
//
private:
//
	TxtCollBaseWriter				stream_;
//
public:
	CatWriterOutputTxtFile(const std::wstring& fname)
		:GenCollWriter<OutputFormatType>(),stream_(fname + std::wstring(FCAT_FILEEXT),std::ios::out)
	{}
	virtual ~CatWriterOutputTxtFile(void)
	{}
};
//
class	CatWriterOutputTxtFile_SZ: public CatWriterOutputTxtFile
{
protected:
//
	inline virtual	void				CreatHeader(const OutputFormatType& data)
	{
		WriteLn(CreatHeaderHelper(data).c_str());
	}

//
	virtual	int					WriteCatLinesColl(const OutputFormatType& data);
//
public:
	CatWriterOutputTxtFile_SZ(const std::wstring& fname)
		:CatWriterOutputTxtFile(fname)
	{}
	virtual ~CatWriterOutputTxtFile_SZ(void)
	{}
};
//
class	CatWriterOutputTxtFile_QA_Contours: public CatWriterOutputTxtFile
{
protected:
//
	inline virtual	void				CreatHeader(const OutputFormatType& data)
	{
		WriteLn(CreatHeaderHelper(data).c_str());
	}

//
	virtual	int					WriteCatLinesColl(const OutputFormatType& data);
//
public:
	CatWriterOutputTxtFile_QA_Contours(const std::wstring& fname)
		:CatWriterOutputTxtFile(fname)
	{}
	virtual ~CatWriterOutputTxtFile_QA_Contours(void)
	{}
};
//
class	CatWriterOutputTxtFile_PS: public CatWriterOutputTxtFile
{
protected:
//
	inline virtual	void				CreatHeader(const OutputFormatType& data)
	{
		WriteLn(CreatHeaderHelper(data).c_str());
	}
//
	virtual	int					WriteCatLinesColl(const OutputFormatType& data);
//
public:
	CatWriterOutputTxtFile_PS(const std::wstring& fname)
		:CatWriterOutputTxtFile(fname)
	{}
	virtual ~CatWriterOutputTxtFile_PS(void)
	{}
};
//
class	ReadNonBlindPtgsTxtFile: public Zeus::ReadTxtFileHelper<CatalogueFormatType>
{
protected:

	virtual	CatalogueFormatType&	do_GetData(void)
	{
		data_.Header_.NRows_ = data_.Storage_.size();
		return	data_;
	}

	virtual void	do_DisposeData(void)
	{data_.Storage_.clear();}

	virtual int		do_parseLn(const std::wstring& str);

	virtual CatalogueFormatType::HeaderType *do_Initialize(void)
	{
		data_.Storage_.clear();
		TxtFileReaderRaw::Open();
		if(!(ReadHeader())) return 0;
		return &(data_.Header_);
	}

public:
	ReadNonBlindPtgsTxtFile(const std::wstring& fname)
		: ReadTxtFileHelper<CatalogueFormatType>(fname + std::wstring(FCAT_FILEEXT))
	{}
	virtual ~ReadNonBlindPtgsTxtFile(void)
	{}
private:
	int							ReadHeader(void);
	CatalogueFormatType			data_;
};
//
class	ReadQA_ColLstText: public Zeus::ReadTxtFileHelper<QA_CltLstCatalogueType>
{
protected:

	virtual	QA_CltLstCatalogueType&	do_GetData(void)
	{return	data_;}

	virtual void	do_DisposeData(void)
	{data_.Storage_.clear();}

	virtual int		do_parseLn(const std::wstring& str);

	virtual QA_CltLstCatalogueType::HeaderType *do_Initialize(void)
	{
		data_.Storage_.clear();
		TxtFileReaderRaw::Open();
		return reinterpret_cast<QA_CltLstCatalogueType::HeaderType *>(&data_);
	}

public:
	ReadQA_ColLstText(const std::wstring& fname)
		: ReadTxtFileHelper<QA_CltLstCatalogueType>(fname + std::wstring(FCAT_FILEEXT))
	{}
	virtual ~ReadQA_ColLstText(void)
	{}
private:
	QA_CltLstCatalogueType		data_;
};
//
class	ReadQA_ProfilesLstText: public Zeus::ReadTxtFileHelper<QA_ProfilesLstType>
{
protected:

	virtual	QA_ProfilesLstType&	do_GetData(void)
	{return	data_;}

	virtual void	do_DisposeData(void)
	{data_.Storage_.clear();}

	virtual int		do_parseLn(const std::wstring& str);

	virtual QA_ProfilesLstType::HeaderType *do_Initialize(void)
	{
		data_.Storage_.clear();
		TxtFileReaderRaw::Open();
		return reinterpret_cast<QA_ProfilesLstType::HeaderType *>(&data_);
	}

public:
	ReadQA_ProfilesLstText(const std::wstring& fname)
		: ReadTxtFileHelper<QA_ProfilesLstType>(fname + std::wstring(FCAT_FILEEXT))
	{}
	virtual ~ReadQA_ProfilesLstText(void)
	{}
private:
	QA_ProfilesLstType		data_;
};
//
} // namespace Zeus

#endif //ZEUSINOUTTXTFILEH

