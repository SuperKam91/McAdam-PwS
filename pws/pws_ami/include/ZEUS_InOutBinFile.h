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
#ifndef ZEUSINOUTBINFILEH
#define ZEUSINOUTBINFILEH

#include "alm_fitsio.h"
#include "fitshandle.h"
#include "ZEUS_InOut.h"
#include "ZEUS_WorkSpace.h"
//------------------------------------
											
//
namespace Zeus
{
//
class	CatWriterOutputFitsFile: public GenCollWriter<OutputFormatType>
{
protected:
#ifdef WG5ALLDOUBLES
	typedef	std::vector<arr<double> >	ColumnsCollType;
#else
	typedef	std::vector<arr<float> >	ColumnsCollType;
#endif
//
//
	inline virtual std::wstring	do_GetCollID(void) const
	{return FName_;}
//
	inline virtual bool			do_Initialize(void)
	{
		FitsHandle_.create(Wstr2Str(FName_));
		return true;
	}
//
	inline virtual void			do_DisposeData(void)
	{FitsHandle_.close();}

//
	inline virtual int			do_Flush(void)
	{return 0;}

//
	inline virtual int			do_Write(const OutputFormatType& data)
	{
		CreateHeader(data);

		CatWriterOutputFitsFile::ColumnsCollType Columns;

		WriteData2Columns(data,Columns);
		return WriteColumns2Fits(Columns);
	}
//
	inline virtual int			do_Remove(void)
	{return std::remove(Wstr2Str(FName_.substr(1)).c_str());}
//
	void						CreateHeaderHelper(const OutputFormatType& data);
//
	virtual void				WriteData2Columns(const OutputFormatType& data,ColumnsCollType& Columns) = 0;
//
	virtual void				CreateHeader(const OutputFormatType& data) = 0;
//
	virtual int					GetNColumns(void) const = 0;
//
	fitshandle					FitsHandle_;
//
public:
	CatWriterOutputFitsFile(const std::wstring& fname)
		:GenCollWriter<OutputFormatType>(),FName_(std::wstring(L"!") + fname + std::wstring(L".fits"))
	{}

	virtual ~CatWriterOutputFitsFile(void)
	{}
private:
//
	inline int					WriteColumns2Fits(const ColumnsCollType& Columns)
	{
		const int NCol(Columns.size());

		for(int i=0;i<NCol;++i)
		{
			FitsHandle_.write_column(i+1,Columns[i]);
		}

		return NCol;
	}
//
	std::wstring				FName_;
};
//
//
class	CatWriterOutputBinFile_PCCFITS: public CatWriterOutputFitsFile
{
protected:
	virtual void	WriteData2Columns(const OutputFormatType& data,CatWriterOutputFitsFile::ColumnsCollType& Columns);
//
	inline virtual void	CreateHeader(const OutputFormatType& data)
	{
		AddHeaderPCCFields(data);
		CreateHeaderHelper(data);
	}
//
	inline virtual int					GetNColumns(void) const
	{return NColumns_;}
//
public:
	CatWriterOutputBinFile_PCCFITS(const std::wstring& fname)
		:CatWriterOutputFitsFile(fname),NColumns_(32)
	{}
	virtual ~CatWriterOutputBinFile_PCCFITS(void)
	{}
private:
	void	AddHeaderPCCFields(const OutputFormatType& data);

	int		NColumns_;
};
//
class	CatWriterOutputBinFile_QA_Contours: public CatWriterOutputFitsFile
{
protected:
	virtual void	WriteData2Columns(const OutputFormatType& data,CatWriterOutputFitsFile::ColumnsCollType& Columns);
//
	inline virtual void	CreateHeader(const OutputFormatType& data)
	{
		AddHeaderPCCFields(data);
		CreateHeaderHelper(data);
	}
//
	inline virtual int					GetNColumns(void) const
	{return NColumns_;}
//
public:
	CatWriterOutputBinFile_QA_Contours(const std::wstring& fname)
		:CatWriterOutputFitsFile(fname),NColumns_(QARESULTSSZ)
	{}
	virtual ~CatWriterOutputBinFile_QA_Contours(void)
	{}
private:
	void	AddHeaderPCCFields(const OutputFormatType& data);

	int		NColumns_;
};
//
class	CatWriterOutputBinFile_PSFITS: public CatWriterOutputFitsFile
{
protected:
	virtual void	WriteData2Columns(const OutputFormatType& data,CatWriterOutputFitsFile::ColumnsCollType& Columns);
//
	inline virtual void	CreateHeader(const OutputFormatType& data)
	{
		AddHeaderPCCFields();
		CreateHeaderHelper(data);
	}
//
	inline virtual int					GetNColumns(void) const
	{return 14;}
//
public:
	CatWriterOutputBinFile_PSFITS(const std::wstring& fname)
		:CatWriterOutputFitsFile(fname)
	{}
	virtual ~CatWriterOutputBinFile_PSFITS(void)
	{}
private:
	void	AddHeaderPCCFields(void);
};
//
//
template<typename T>
class	WriteWrkSpaceBinFile: public GenCollWriter<LArr2D<T> >
{
	typedef LArr2D<T>	AtomT;
protected:
	inline virtual std::wstring	do_GetCollID(void) const
	{return FName_;}
//
	inline virtual bool	do_Initialize(void)
	{
		if(!(stream_ = std::fopen(Wstr2Str(FName_).c_str(),"wb")))
			return false;
		return true;
	}
//
	inline virtual void			do_DisposeData(void)
	{Dispose();return ;}
//
	inline virtual int			do_Flush(void)
	{
		if(stream_) std::fflush(stream_);
		return 0;
	}
//
	inline virtual int		do_Remove(void)
	{return std::remove(Wstr2Str(FName_).c_str());}
//
	virtual int	do_Write(const AtomT & data)
	{
		if(!stream_) return -1;
		int	NItems(data.getSz());
		if(NItems == std::fwrite(data.begin(),sizeof(T),NItems,stream_))
		{
			std::fflush(stream_);
			return NItems;
		}
		return -1;
	}
public:
	WriteWrkSpaceBinFile(const std::wstring& fname)
		: GenCollWriter<LArr2D<T> >(),FName_(fname + std::wstring(BIN_FILEEXT)),stream_(0)
	{}
	virtual ~WriteWrkSpaceBinFile(void)
	{Dispose();}
private:
	inline void	Dispose(void)
	{
		if(stream_)
		{std::fclose(stream_);stream_ = 0;}	
	}

	std::wstring	FName_;
	FILE			*stream_;
};
//
template<typename T>
class	WriteWrkSpaceFitsFile: public GenCollWriter<LArr2D<T> >
{
	typedef LArr2D<T>	AtomT;
protected:
	inline virtual std::wstring	do_GetCollID(void) const
	{return FName_;}
//
	inline virtual bool	do_Initialize(void)
	{
		if(!stream_)
		{
			stream_ = new fitshandle();
		}

		stream_->create(Wstr2Str(FName_));

		return true;
	}
//
	inline virtual void			do_DisposeData(void)
	{Dispose();}
//
	inline virtual int			do_Flush(void)
	{
		if(stream_) return 0;
		return -1;
	}
//
	inline virtual int		do_Remove(void)
	{
		if(stream_)
		{
			stream_->delete_file(Wstr2Str(FName_.substr(1)));
			return 0;
		}
		return -1;
	}
//
	virtual int	do_Write(const AtomT & data)
	{
		if(!stream_) return -1;
		
		const	int	 Metric(static_cast<int>(data.getPtrMetric()));
		const	int	 NItems(static_cast<int>(data.getSz()));

		arr2<T>		fitsArr(NItems / Metric, Metric);
		T*			arr2Ptr(fitsArr[0]);
		const T*	myPtr(data.begin());

		for(int i=0;i<NItems;++i,++arr2Ptr,++myPtr)
		{*arr2Ptr = *myPtr;}

		stream_->insert_image(planckType<T>(), fitsArr);

		stream_->set_key("MAXCY5R500",YMax_,"Y axis MAX");
		stream_->set_key("MINCY5R500",YMin_,"Y axis MIN");
		stream_->set_key("MAXRS",XMax_,"X axis MAX");
		stream_->set_key("MINRS",XMin_,"X axis MIN");
		if((aux_.In_Y_ != SZCATALOGUEDEFAULTVALUE) && (aux_.In_Theta_ != SZCATALOGUEDEFAULTVALUE))
		{
			stream_->set_key("Yin",aux_.In_Y_,"Y IN");
			stream_->set_key("Tin",aux_.In_Theta_,"X IN");
		}
		return NItems;

	}
//
public:
	WriteWrkSpaceFitsFile(const std::wstring& fname, double YMax, double YMin, double XMax, double XMin,ContoursOutAux* aux)
		: GenCollWriter<LArr2D<T> >(),FName_(L"!" + fname + std::wstring(FITS_FILEEXT)),YMax_(YMax),YMin_(YMin),
		XMax_(XMax),XMin_(XMin),stream_(0)
	{
		if(aux != NULL) aux_ = *aux;
		stream_ = new fitshandle();
	}
//
	virtual ~WriteWrkSpaceFitsFile(void)
	{Dispose();}
//
private:
	inline void	Dispose(void)
	{
		if(stream_)
		{
			stream_->close();
			delete stream_;
			stream_ = 0;
		}	
	}
//
	double			YMax_;
	double			YMin_;
	double			XMax_;
	double			XMin_;
	std::wstring	FName_;
	fitshandle		*stream_;
	ContoursOutAux	aux_;
};
//
template<typename T>
class	ReadWrkSpaceBinFile:  public GenCollReader<LArr2D<T> >
{
private:
	std::wstring	FName_;
	LArr2D<T>		*data_;
	FILE			*stream_;
	int				YSz_;
	int				Metric_;
protected:

	virtual std::wstring	do_GetCollID(void) const
	{return FName_;}

	virtual bool	do_Read(void)
	{
		int NItems(data_->getSz());
		if(NItems == std::fread(data_->begin(),sizeof(T),NItems,stream_))
			return true;
		return false;
	}

	virtual	LArr2D<T>&	do_GetData(void)
	{return *data_;}

	virtual void	do_DisposeData(void)
	{
		delete data_;data_ = 0;
		if(stream_)
		{std::fclose(stream_);stream_ = 0;}
	}

	virtual typename LArr2D<T>::HeaderType *do_Initialize(void)
	{
		char buffer[BUFFERMAXCHAR];
		do_DisposeData();
		if(!(stream_ = std::fopen(Wstr2Achar(FName_,buffer),"rb")))
			return 0;
		data_ = new LArr2D<T>(YSz_ * Metric_ ,Metric_);
		data_->reset();
		return reinterpret_cast<typename LArr2D<T>::HeaderType *>(data_);
	} 
public:
	ReadWrkSpaceBinFile(const std::wstring& fname,int YSz,int XSz)
		: GenCollReader<LArr2D<T> >(),FName_(fname  + std::wstring(BIN_FILEEXT)),YSz_(YSz),Metric_(XSz),data_(0),stream_(0)
	{}
	virtual ~ReadWrkSpaceBinFile(void)
	{
		do_DisposeData();
	}
};
//
template<typename T>
class	ReadWrkSpaceFitsFile:  public GenCollReader<LArr2D<T> >
{
private:
	std::wstring	FName_;
	LArr2D<T>		*data_;
	fitshandle		*stream_;
	int				YSz_;
	int				Metric_;
protected:

	virtual std::wstring	do_GetCollID(void) const
	{return FName_;}

	virtual bool	do_Read(void)
	{
		if(!stream_) return false;


		arr2<T>		fitsArr(YSz_, Metric_);
		T*			myPtr(data_->begin());
		
		stream_->open(Wstr2Str(FName_));
		stream_->goto_hdu(2); // Images are stored into HDU 2
		stream_->read_image(fitsArr);

		const T*	arr2Ptr(fitsArr[0]);

		const int NItems(data_->getSz());

		for(int i=0;i<NItems;++i,++arr2Ptr,++myPtr)
		{*myPtr = *arr2Ptr;}

		return true;
	}

	virtual	LArr2D<T>&	do_GetData(void)
	{return *data_;}

	virtual void	do_DisposeData(void)
	{
		delete data_;data_ = 0;
		if(stream_)
		{
			stream_->close();
			delete stream_;
			stream_ = 0;
		}	
	}

	virtual typename LArr2D<T>::HeaderType *do_Initialize(void)
	{
		stream_ = new fitshandle();

		data_ = new LArr2D<T>(YSz_ * Metric_ ,Metric_);
		data_->reset();

		return reinterpret_cast<typename LArr2D<T>::HeaderType *>(data_);
	} 
public:
	ReadWrkSpaceFitsFile(const std::wstring& fname,int YSz,int XSz)
		: GenCollReader<LArr2D<T> >(),FName_(fname  + std::wstring(FITS_FILEEXT)),YSz_(YSz),Metric_(XSz),data_(0),stream_(0)
	{}
	virtual ~ReadWrkSpaceFitsFile(void)
	{
		do_DisposeData();
	}
};
//
template<typename T>
class	ReadVectBinFile:  public GenCollReader<LArr1D<T> >
{
private:
	std::wstring	FName_;
	LArr1D<T>		*data_;
	FILE			*stream_;
	int				Sz_;
protected:

	virtual std::wstring	do_GetCollID(void) const
	{return FName_;}

	virtual bool	do_Read(void)
	{
		int NItems(data_->getSz());
		if(NItems == std::fread(data_->begin(),sizeof(T),NItems,stream_))
			return true;
		return false;
	}

	virtual	LArr1D<T>&	do_GetData(void)
	{return *data_;}

	virtual void	do_DisposeData(void)
	{
		delete data_;data_ = 0;
		if(stream_)
		{std::fclose(stream_);stream_ = 0;}
	}

	virtual typename LArr1D<T>::HeaderType *do_Initialize(void)
	{
		int tSZ;
		char buffer[BUFFERMAXCHAR];
		do_DisposeData();
		if (!(stream_ = std::fopen(Wstr2Achar(FName_, buffer), "rb")))
			return 0;
		if (Sz_ < 0)
		{
			std::fseek(stream_, 0, SEEK_END);
			if ((tSZ = std::ftell(stream_)) < 0)
				return 0;
			std::fseek(stream_, 0, SEEK_SET);
		}
		data_ = new LArr1D<T>(Sz_ = (tSZ / sizeof(T)));
		data_->reset();
		return reinterpret_cast<typename LArr1D<T>::HeaderType *>(data_);
	} 
public:
	ReadVectBinFile(const std::wstring& fname,int Sz)
		: GenCollReader<LArr1D<T> >(),FName_(fname  + std::wstring(BIN_FILEEXT)),Sz_(Sz),data_(0),stream_(0)
	{}
	virtual ~ReadVectBinFile(void)
	{
		do_DisposeData();
	}
};
//
template<typename T>
class	WriteVectBinFile: public GenCollWriter<LArr1D<T> >
{
	typedef LArr1D<T>	AtomT;
protected:
	inline virtual std::wstring	do_GetCollID(void) const
	{return FName_;}
//
	virtual bool	do_Initialize(void)
	{
		if(!(stream_ = std::fopen(Wstr2Str(FName_).c_str(),"wb")))
			return false;
		return true;
	}
//
	inline virtual int			do_Flush(void)
	{
		if(stream_) std::fflush(stream_);
		return 0;
	}
//
	inline virtual void			do_DisposeData(void)
	{Dispose();return ;}
//
	inline virtual int		do_Remove(void)
	{return std::remove(Wstr2Str(FName_).c_str());}
//
	inline virtual int	do_Write(const AtomT & data)
	{
		if(!stream_) return -1;
		int	NItems(data.getSz());
		if(NItems == std::fwrite(data.begin(),sizeof(T),NItems,stream_))
		{
			std::fflush(stream_);
			return NItems;
		}
		return -1;
	}
public:
	WriteVectBinFile(const std::wstring& fname)
		: GenCollWriter<LArr1D<T> >(),FName_(fname + std::wstring(BIN_FILEEXT)),stream_(0)
	{}
//
	virtual ~WriteVectBinFile(void)
	{Dispose();}
private:
	inline void	Dispose(void)
	{
		if(stream_)
		{std::fclose(stream_);stream_ = 0;}	
	}

	std::wstring	FName_;
	FILE			*stream_;
};
//



}


#endif //ZEUSINOUTBINFILEH

