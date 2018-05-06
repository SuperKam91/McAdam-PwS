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

#include <memory>

#ifdef LFIDPC
#include "healpix_map_dmcio.h"
#endif

#include "ZEUS_InOut.h"
#include "ZEUS_InOutBinFile.h"
#include "ZEUS_InOutHealpixFits.h"
#include "ZEUS_InOutPipeline.h"
#include "ZEUS_InOutTxtFile.h"

//---------------------------------------------------------------------------

namespace Zeus
{

const wchar_t*	Private::NumTxtFormatType<int>::PrintfFormat	= L"%0*d";
const wchar_t*	Private::NumTxtFormatType<double>::PrintfFormat = L"%.*g";
const wchar_t*	Private::NumTxtFormatType<float>::PrintfFormat = L"%.*g";
const wchar_t*	Private::NumTxtFormatType<std::complex<double> >::PrintfFormat = L"%.*g,%.*g";
const wchar_t*	Private::NumTxtFormatType<std::complex<float> >::PrintfFormat = L"%.*g,%.*g";

const int		PCCFormatHeaderType::CatNColumns_		= 256;
const int		PCCFormatHeaderType::CatHeaderNColumns_ = 18;
const wchar_t*	PCCFormatHeaderType::CatPrintStr_ = L"%05d,%05d,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%05d,%ls %ls";
const wchar_t*	PCCFormatHeaderType::CatHeaderPrintStr_= L"%4d,%4d,%6.1g,%4d,%4d,%4d,%4d,%4d,%4d,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%9.6g,%4d";
const wchar_t*	PCCFormatHeaderType::CatReadStr_ = L"%d,%d,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%le,%d,%ls %ls";
const wchar_t*	PCCFormatHeaderType::CatHeaderReadStr_ = L"%d,%d,%le,%d,%d,%d,%d,%d,%d,%le,%le,%le,%le,%le,%le,%le,%le,%d";

#ifdef HFIDMC
 PIOSTRING		Private::PIOTypeChooserType<int>::PIO_TYPE			= "PIOINT";
 PIOSTRING		Private::PIOTypeChooserType<float>::PIO_TYPE		= "PIOFLOAT";
 PIOSTRING		Private::PIOTypeChooserType<double>::PIO_TYPE		= "PIODOUBLE";
 PIOSTRING		Private::PIOTypeChooserType<bool>::PIO_TYPE			= "PIOINT";
 PIOSTRING		Private::PIOTypeChooserType<std::string>::PIO_TYPE	= "PIOSTRING";
#endif

template<typename T>
void	HealpixHeaderAtomType::Set(const std::string& key,const T& value)
{
	Keyword_		= key;
	ValueType_		= static_cast<HealpixHeaderAtomInnerType>(Private::HHAtomTypeHelper<T>::VALUETYPE);
	Private::HHAtomTypeHelper<T>::Storage(*this) = value;
}

template<typename T>
void	HealpixHeaderAtomType::Set(const T& value)
{
	ValueType_		= static_cast<HealpixHeaderAtomInnerType>(Private::HHAtomTypeHelper<T>::VALUETYPE);
	Private::HHAtomTypeHelper<T>::Storage(*this) = value;
}

template void HealpixHeaderAtomType::Set(const std::string& key,const std::string& value);
template void HealpixHeaderAtomType::Set(const std::string& key,const double& value);
template void HealpixHeaderAtomType::Set(const std::string& key,const int& value);
template void HealpixHeaderAtomType::Set(const std::string& key,const float& value);
template void HealpixHeaderAtomType::Set(const std::string& key,const bool& value);
template void HealpixHeaderAtomType::Set(const std::string& value);
template void HealpixHeaderAtomType::Set(const double& value);
template void HealpixHeaderAtomType::Set(const int& value);
template void HealpixHeaderAtomType::Set(const float& value);
template void HealpixHeaderAtomType::Set(const bool& value);
//
template<typename T>
GenHealpixWriter<T>
*GetGenHealpixWriter(Loki::Type2Type<T>,int ID,const std::wstring& DirName,const std::wstring& MapName,
					 int NSide,Healpix_Ordering_Scheme Scheme,coordsys CoordSys)
{
	switch(ID)
	{
	case ENVID_HFI:
#ifdef HFIDMC
	return new WriteHealpixMapHFIDMC<T>(DirName,MapName,NSide,Scheme,CoordSys);
#else
	break;
#endif
	case ENVID_LFI:
#ifdef LFIDPC
	return new WriteHealpixMapLFIDPC<T>(DirName + MapName,Scheme,CoordSys);
#else
	break;
#endif
	case ENVID_FITS:
		// Healpix files are binary files
		break;
	default:
		break;	
	}

	return new WriteHealpixMapFits<T>(DirName + MapName,Scheme,CoordSys);
}
//
template GenHealpixWriter<float>
*GetGenHealpixWriter(Loki::Type2Type<float>,int ID,const std::wstring& DirName,const std::wstring& MapName,
					 int NSide,Healpix_Ordering_Scheme Scheme,coordsys CoordSys);
//
template GenHealpixWriter<double>
*GetGenHealpixWriter(Loki::Type2Type<double>,int ID,const std::wstring& DirName,const std::wstring& MapName,
					 int NSide,Healpix_Ordering_Scheme Scheme,coordsys CoordSys);
//
template<typename T>
GenHealpixReader<T>
*GetGenHealpixReader(Loki::Type2Type<T>,int ID,const std::wstring& DirName,const std::wstring& MapName)
{
	switch(ID)
	{
	case ENVID_HFI:
#ifdef HFIDMC
		return new ReadHealpixMapHFIDMC<T>(DirName,MapName);
#else
	break;
#endif
	case ENVID_LFI:
#ifdef LFIDPC
		return new ReadHealpixMapLFIDPC<T>(DirName,MapName);
#else
	break;
#endif
	case ENVID_FITS:
		// Healpix files are binary files
		break;
	default:
		break;	
	}

	return new ReadHealpixMapFits<T>(DirName + MapName);
}
//
template GenHealpixReader<float>
*GetGenHealpixReader(Loki::Type2Type<float>,int ID,const std::wstring& DirName,const std::wstring& MapName);
//
template GenHealpixReader<double>
*GetGenHealpixReader(Loki::Type2Type<double>,int ID,const std::wstring& DirName,const std::wstring& MapName);
//
template<typename T>
GenCollReader<LArr2D<T> >
*GetWrkSpFileReaderHandler(Loki::Type2Type<LArr2D<T> >,int ContextID,const std::wstring& fn, int YSz,int Metric,const std::wstring& DirN)
{
	switch(ContextID)
	{
	case 1:
#if defined(HFIDMC)
		return new ReadWrkSpaceHFIDMC<T>(DirN,fn,YSz,Metric);
#else
	break;
#endif
	case 2:
#ifdef LFIDPC
		return new ReadWrkSpaceLFIDPC<T>(DirN,fn,YSz,Metric);
#else
	break;
#endif
	case 3:
		return new ReadWrkSpaceFitsFile<T>(DirN + fn,YSz,Metric);
	default:
	break;	
	}
	return new ReadWrkSpaceFitsFile<T>(DirN + fn,YSz,Metric);
}
//
template GenCollReader<LArr2D<double> >
*GetWrkSpFileReaderHandler(Loki::Type2Type<LArr2D<double> >,int ContextID,const std::wstring& fn, int YSz,int Metric,const std::wstring& DirN);
//
template<typename T>
Zeus::GenCollWriter<LArr2D<T> >
*GetGenCollFileWriterHandler(Loki::Type2Type<LArr2D<T> >,int ID,const std::wstring& fname,const std::wstring& DirN,int AxSz,double PixSz,double PatchCx,double PatchCy,double XMax,double XMin,ContoursOutAux* aux)
{
	switch(ID)
	{
	case ENVID_HFI:
#if defined(HFIDMC)
		return new WriteWrkSpaceHFIDMC<T>(DirN,fname,AxSz,PixSz,PatchCx,PatchCy,XMax,XMin);
#else
	break;
#endif
	case ENVID_LFI:
#ifdef LFIDPC
	return new WriteWrkSpaceLFIDPC<T>(DirN,fname);
#else
	break;
#endif
	case ENVID_FITS:
	default:
	break;	
	}

	Zeus::CreateDir(DirN);

	return new WriteWrkSpaceFitsFile<T>(DirN + fname,PatchCx,PatchCy,XMax,XMin,aux);
}
//
template Zeus::GenCollWriter<LArr2D<float> >
*GetGenCollFileWriterHandler(Loki::Type2Type<LArr2D<float> >,int ID,const std::wstring& fname,const std::wstring& DirN,int AxSz,double PixSz,double PatchCx,double PatchCy,double XMax,double XMin,ContoursOutAux* aux);
//
template Zeus::GenCollWriter<LArr2D<double> >
*GetGenCollFileWriterHandler(Loki::Type2Type<LArr2D<double> >,int ID,const std::wstring& fname,const std::wstring& DirN,int AxSz,double PixSz,double PatchCx,double PatchCy,double XMax,double XMin,ContoursOutAux* aux);
//
template<typename T>
GenCollReader<LArr1D<T> >
*GetVectFileHandlerReader(Loki::Type2Type<LArr1D<T> >,int ContextID,const std::wstring& fn, int Sz,const std::wstring& DirN)
{
	switch(ContextID)
	{
	case 1:
#ifdef HFIDMC
		return new ReadVectHFIDMC<T>(DirN,fn,Sz);
#else
	break;
#endif
	case 2:
#ifdef LFIDPC
		return new ReadVectLFIDPC<T>(DirN,fn,Sz);
#else
	break;
#endif
	case 3:
		return new ReadVectBinFile<T>(DirN + fn,Sz);
	default:
	break;	
	}
	// TODO vector txt file
	return new ReadVectBinFile<T>(DirN + fn,Sz);
}
//
template GenCollReader<LArr1D<double> >
*GetVectFileHandlerReader(Loki::Type2Type<LArr1D<double> >,int ContextID,const std::wstring& fn, int Sz,const std::wstring& DirN);
//
template<typename T>
Zeus::GenCollWriter<LArr1D<T> >
*GetVectFileHandlerWriter(Loki::Type2Type<LArr1D<T> >,int ID,const std::wstring& fname,const std::wstring& DirN)
{
	switch(ID)
	{
	case ENVID_HFI:
#ifdef HFIDMC
	return new WriteVectHFIDMC<T>(DirN,fname);
#else
	break;
#endif
	case ENVID_LFI:
#ifdef LFIDPC
	return new WriteVectLFIDPC<T>(DirN,fname);
#else
	break;
#endif
	case ENVID_FITS:
		return new WriteVectBinFile<T>(DirN + fname);
	default:
	break;	
	}
	// TODO Write txt file
	return new WriteVectBinFile<T>(DirN + fname);
}
//
template Zeus::GenCollWriter<LArr1D<double> >
*GetVectFileHandlerWriter(Loki::Type2Type<LArr1D<double> >,int ID,const std::wstring& fname,const std::wstring& DirN);
//
GenCollWriter<OutputFormatType>
*GetFormatCatWriterHandler(int ContextID,int DetectionType, const std::wstring& fn,const std::wstring& GrpName)
{
	switch(ContextID)
	{
	case 1:
#ifdef HFIDMC
	if(DetectionType==1000)
	{return new CatWriterOutputHFIDMC_QA_Contours(GrpName,fn);}
	else
	{
		if(DetectionType)
		{return new CatWriterOutputHFIDMC_Orange10(GrpName,fn);}
		else
		{return new CatWriterOutputHFIDMC_Orange18(GrpName,fn);}
	}
#else
	break;
#endif
	case 2:
#ifdef LFIDPC
# error LFI interface not usable anymore
#else
	break;
#endif
//
	case 3:
	if(DetectionType==1000)
	{return new CatWriterOutputBinFile_QA_Contours(GrpName + fn);}
	else
	{
		if(DetectionType)
		{return new CatWriterOutputBinFile_PCCFITS(GrpName + fn);}
		else
		{return new CatWriterOutputBinFile_PSFITS(GrpName + fn);}
	}
//
	default:
		break;
	}
//
	if(DetectionType==1000)
	{return new CatWriterOutputTxtFile_QA_Contours(GrpName + fn);}
	else
	{
		if(DetectionType)
		{return new CatWriterOutputTxtFile_SZ(GrpName + fn);}
		else
		{return new CatWriterOutputTxtFile_PS(GrpName + fn);}
	}
}
//
GenCollWriter<CatalogueFormatType>
*GetCatWriterHandler(int ContextID,const std::wstring& fn,const std::wstring& GrpName,int NColumns)
{
	switch(ContextID)
	{
	case 1:
#ifdef HFIDMC
		// catalogues have always 254-256 columns
		return new CatalogueWriterHFIDMC(GrpName,fn,256);
#else
	break;
#endif
	case 2:
#ifdef LFIDPC
	return new CatalogueWriterLFIDPC(GrpName,fn,ProfParam);  // CatWriter for PS Clusters	
#else
	break;
#endif
	case 3:
		return new CatWriterTxtFile(GrpName + fn);
	default:
		break;
	}

	return new CatWriterTxtFile(GrpName + fn); // CatWriter for PS Clusters
}
//
GenCollReader<PeakCollReadbleType>
*GetObjReaderHandler(int ContextID,const std::wstring& dir,const std::wstring& fn)
{
	switch(ContextID)
	{
	case 1:
#ifdef HFIDMC
	return new ReadObjResultsHFIDMC(dir,fn);
#else
	break;
#endif
	case 2:
#ifdef LFIDPC
	return new ReadObjResultsLFIDPC(fn);
#else
	break;
#endif
	case 3:
		break; //TODO use FITS files
	default:
	break;	
	}

	return new ReadObjResultsTxtFile(dir + fn);
}
//
GenCollWriter<PeakCollType>
*GetObjWriterHandler(int ContextID,const std::wstring& Dir,const std::wstring& fn, int InitObjID)
{
	switch(ContextID)
	{
	case 1:
#ifdef HFIDMC
	return new WriteObjResultsHFIDMC(Dir,fn,InitObjID);
#else
	break;
#endif
	case 2:
#ifdef LFIDPC
	return new WriteObjResultsLFIDPC(fn,InitObjID);
#else
	break;
#endif
	case 3:
		// TODO output to FITS file
		break;
	default:
	break;	
	}
	return new WriteObjResultsTxtFile(Dir + fn,InitObjID);
}
//
GenCollReader<PatchGeomType>
*GetPatchGeomFileReaderHandler(Loki::Type2Type<PatchGeomType>,int ContextID,const std::wstring& DirName,const std::wstring& fn)
{
	switch(ContextID)
	{
	case 1:
#ifdef HFIDMC
	return new ReadPatchesGeomInfoHFIDMC(DirName,fn);
#else
	break;
#endif
	case 2:
#ifdef LFIDPC
	return new ReadPatchesGeomInfoLFIDPC(DirName,fn);
#else
	break;
#endif
	case 3:
		break;
	default:
	break;	
	}
	return new ReadPatchesGeomInfoTxtFile(DirName + fn);
}
//
GenCollWriter<PatchGeomType>
*GetGenCollFileWriterHandler(Loki::Type2Type<PatchGeomType>,int ID,const std::wstring& DirName,const std::wstring& fname)
{
	switch(ID)
	{
	case ENVID_HFI:
#ifdef HFIDMC
	return new WritePatchesGeomInfoHFIDMC(DirName,fname);
#else
	break;
#endif
	case ENVID_LFI:
#ifdef LFIDPC
	return new WritePatchesGeomInfoLFIDPC(DirName,fname);
#else
	break;
#endif
	case ENVID_FITS:
		break;
	default:
	break;	
	}
	return new WritePatchesGeomInfoTxtFile(DirName + fname);
}
//
GenCollReader<CatalogueFormatType>
*GetGenCollFileReaderHandler(Loki::Type2Type<CatalogueFormatType>,int ID,const std::wstring& Dname,const std::wstring& fname)
{

	switch(ID)
	{
	case ENVID_HFI:
#ifdef HFIDMC
	return new ReadNonBlindCatalogueHFIDMC(Dname,fname);
#else
	break;
#endif
	case ENVID_LFI:
#ifdef LFIDPC
	return  new ReadNonBlindCatalogueLFIDPC(Dname + fname);
#else
	break;
#endif
	case ENVID_FITS:
	break;
	// No Binary format
	default:
	break;	
	}

	return new ReadNonBlindPtgsTxtFile(Dname + fname);
}
//
GenCollReader<QA_CltLstCatalogueType>
*GetQA_ColLstReaderHandler(int ID,const std::wstring& Dname,const std::wstring& fname)
{
	switch(ID)
	{
	case ENVID_HFI:
#ifdef HFIDMC
	return new ReadQA_ColLstHFIDMC(Dname,fname);
#else
	break;
#endif
	case ENVID_LFI:
#ifdef LFIDPC
#error LFI Not supported anymore
#else
	break;
#endif
	case ENVID_FITS:
	break;
	// No Binary format
	default:
	break;	
	}

	return new ReadQA_ColLstText(Dname + fname);
}
//
GenCollReader<QA_ProfilesLstType>
*GetQA_ProfilesLstReaderHandler(int ID,const std::wstring& Dname,const std::wstring& fname)
{
	switch(ID)
	{
	case ENVID_HFI:
#ifdef HFIDMC
	return new ReadQA_ProfilesLstHFIDMC(Dname,fname);
#else
	break;
#endif
	case ENVID_LFI:
#ifdef LFIDPC
#error LFI Not supported anymore
#else
	break;
#endif
	case ENVID_FITS:
	break;
	// No Binary format
	default:
	break;	
	}

	return new ReadQA_ProfilesLstText(Dname + fname);
}
//
//
void	CreatReadableHealpixFile(int EnvId,const std::wstring& DirName,const std::wstring& FName,
								 coordsys CoordSys,const Healpix_Map<float>& HealpixData)
{
	std::auto_ptr<GenHealpixWriter<float> >	HP_MapWriter(Zeus::GetGenHealpixWriter(Loki::Type2Type<float>(),EnvId,DirName,FName,
		HealpixData.Nside(),HealpixData.Scheme(),CoordSys));
	HP_MapWriter->Initialize();

	HP_MapWriter->Write(HealpixData,0);

}
//
std::wstring&	MakePath(std::wstring& str)
{
	std::wstring::iterator			piv(str.begin());
	std::wstring::const_iterator	const end(str.end());

	for(;piv!=end;++piv){
#ifdef WIN32
		if(*piv == L'/') *piv = L'\\';
#else
		if(*piv == L'\\') *piv = L'/';	
#endif
	}
	return str;
}
//
// 0 if file/object exists

//enum OBJECTTYPE {OBJECTTYPE_MAP=0,OBJECTTYPE_TABLE=1,OBJECTTYPE_VECT=2};

int		RemoveFile(OBJECTTYPE ObjT,const std::wstring& Dir,const std::wstring& file)
{
#if defined(HFIDMC)
	return HFIDMC_RemoveObject(Dir,file);
#elif defined(LFIDPC)
	return LFIDPC_RemoveObject(ObjT,Dir,file);
#else
	std::wstring	ext;
	switch(ObjT)
	{
	case OBJECTTYPE_MAP:
		ext = L".fits";
		break;
	case OBJECTTYPE_TABLE:
		ext = L".csv";
		break;
	case OBJECTTYPE_VECT:
		ext = L".dat";
		break;
	}
	return remove(Wstr2Str(Dir + file + ext).c_str());
#endif
}
//
std::wstring	CorrectDir(const std::wstring& dir,int ForceAppend)
{
	std::wstring	temp(FullTrim(dir));

	if(temp.empty())
	{return temp;}
#ifdef WIN32	
	wchar_t			Suffix(L'\\');
#else
	wchar_t			Suffix(L'/');
#endif
#if defined(HFIDMC) ||  defined(LFIDPC)
	int append(0);
#else
	int append(1);
#endif
	if(ForceAppend) append=1;

	if(append)
	{
		if(*(temp.end()-1) != Suffix)
			temp.push_back(Suffix);
	}
	else
	{
		if(*(temp.end()-1) == Suffix)
			temp = temp.substr(0,temp.size() - 1);
	}

	return temp;
}
//
//
#if	defined(WIN32) && !defined(HFIDMC) && !defined(LFIDPC) && !defined(AMI)
void	ConManagerClass::Initialise(int argc, _TCHAR* argv[],int ExecId)
#else
void	ConManagerClass::Initialise(int argc, char* argv[],int ExecId)
#endif
{
	CloseStreams();
#ifdef PWSMPI
	execID_ = 0;
#else
	execID_ = ExecId;
#endif //PWSMPI
	if(!(execID_ >> 2))
	{
		execID_ = (GetPID() << 2) + (execID_ & 0x03);
	}

#ifdef	HFIDMC

#undef __PIOFUNCTION__
#define __PIOFUNCTION__ "ConManagerClass::Initialise"

	MyErr = PIOInitLogger(&Logger_,argv);
	if (MyErr<0) {
		fprintf(stderr,"%s: Cannot initialize the logger !!!\n",argv[0]);
		exit(MyErr);
	}
/*
	if (argc<2) {
		sprintf(logStr,"Usage : %s <filepar>\n<filepar : Parameter file\n",argv[0]);
		_PIOEXITCAUSE_(Logger_,logStr,-1);
	}
*/
	if(execID_ & ((int) 0x01))
	{
		MsgFile_ = 1;
	}
	if(execID_ & ((int) 0x02))
	{
		ErrFile_ = 1;
	}
#else
//
#ifdef	LFIDPC
	module_startup("r_param", argc, const_cast<const char **>(&(argv[0])), 2, "<init object>");
#endif
	std::string	logFName;

	if(execID_ & ((int) 0x01))
	{
		logFName = GetLogFileName(MSGSLOGFILENAME);
		MsgFile_ = new std::wofstream(logFName.c_str(),std::ios::app);
		if(!MsgFile_ ||  MsgFile_->fail())
			throw Zeus::libException(ERROR_COD_MYINOUT_FILEWRITE,ERROR_MSG_MYINOUT_FILEWRITE,logFName.c_str());
		MsgFile_->clear();
		MsgFile_->exceptions(std::ios::failbit | std::ios::badbit);
	}

	if(execID_ & ((int) 0x02))
	{
		logFName = GetLogFileName(ERRORLOGFILENAME);
		ErrFile_ = new std::wofstream(logFName.c_str(),std::ios::app);
		if(!ErrFile_ ||  ErrFile_->fail())
			throw Zeus::libException(ERROR_COD_MYINOUT_FILEWRITE,ERROR_MSG_MYINOUT_FILEWRITE,logFName.c_str());
		ErrFile_->clear();
		ErrFile_->exceptions(std::ios::failbit | std::ios::badbit);
	}
#endif

}


}

