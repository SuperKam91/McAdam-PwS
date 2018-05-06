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

#ifndef ZEUSINOUTPIPELINEH
#define ZEUSINOUTPIPELINEH

#include "ZEUS_Strings.h"
#include "ZEUS_InOut.h"

//------------------------------------
//

#define	GEOMPROPNINTS							4
#define	GEOMPROPNFLOATS							16
#define NONBLIND_NFIELDS						31
#define	GEOINFONCOLUMNS							256
#define	INTERCATNCOLUMNS						256
#define MAXRETRIES								3

#define HFIDMCCATSGROUPEXTENSION				L"_Cats"
#define HFIDMCNONBCATGROUPEXTENSION				L"_NBCats"
#define HFIDMCINTERMCATGROUPEXTENSION			L"_IntCats"
#define HFIDMCPACTHGEOGROUPEXTENSION			L"_PatGeo"
#define HFIDMCPROFCOEFGROUPEXTENSION			L"_ProfCoefs"

#define	INTFILE_LFIDPC_COLLTYPE					""
#define	IMGCOLL_LFIDPC_COLLTYPE					""
#define	HEALPIXCOLL_LFIDPC_COLLTYPE				""

#define	HFIDMC_INSERTVALUE_AUX(ARGUMENT) if((data_->Header_->GetParContent())->flag_##ARGUMENT && (DMC_GetValue(Zeus::Achar2Wstr((data_->Header_->GetParContent())->ARGUMENT),tDBField) > -100)) \
											Insert1element(std::wstring(MAKEWCHARARR(ARGUMENT)),tDBField)
#define	HFIDMC_INSERTVALUE(ARGUMENT) HFIDMC_INSERTVALUE_AUX(ARGUMENT)

#define	HFIDMC_INSERTVALUENOOPT_AUX(ARGUMENT) if(DMC_GetValue(Zeus::Achar2Wstr((data_->Header_->GetParContent())->ARGUMENT),tDBField) > -100) \
											Insert1element(std::wstring(MAKEWCHARARR(ARGUMENT)),tDBField)
#define	HFIDMC_INSERTVALUENOOPT(ARGUMENT) HFIDMC_INSERTVALUENOOPT_AUX(ARGUMENT)

#define	HFIDMC_INSERTGRPNAME_AUX(PAR_NAME,PWS_NAME) if(DMC_GetValue(Zeus::Achar2Wstr((data_->Header_->GetParContent())->info_##PAR_NAME##.groupName),tDBField) > -100) \
											Insert1element(std::wstring(MAKEWCHARARR(PWS_NAME)),tDBField)
#define	HFIDMC_INSERTGRPNAME(PAR_NAME,PWS_NAME) HFIDMC_INSERTGRPNAME_AUX(PAR_NAME,PWS_NAME)


#define	HFIDMC_READVECTCOLUMN(STORAGE,FNAME,TYPE) MyErr=PIOReadVECTObject((void **)&(STORAGE),FNAME,TYPE,command,VECTgrp_); \
												if (MyErr<0) \
													HFIDMC_errIO_Vect(MyErr,ERROR_MSG_HFIDMCREADVECTFORMAT,FNAME)

#define	HFIDMC_DELETEVECTCOLUMN(STORAGE,FNAME) MyErr=PIODeleteVECTTable(STORAGE,VECTgrp_);	\
												if (MyErr != 0 )	\
													HFIDMC_ReportMemoryLeak((PIOErr) MyErr,ERROR_MSG_HFIDMCDELETEVECTFORMAT,FNAME)

#define	HFIDMC_READKEYWORD(STORAGE,DATANAME,FIELDNAME,OWNER) if((MyErr = PIOReadKeywordObject(((void *) &(DATANAME) ),dummyStr,FIELDNAME,STORAGE,OWNER,VECTgrp_)) != 0) \
																HFIDMC_errIO_Vect(MyErr,ERROR_MSG_HFIDMCREADKEYWORDFORMAT,FIELDNAME)

#define	HFIDMC_READKEYWORDMAP(STORAGE,DATANAME,FIELDNAME,OWNER) if((MyErr = PIOReadKeywordObject(((void *) &(DATANAME) ),dummyStr,FIELDNAME,STORAGE,OWNER,MAPgrp_)) != 0) \
																	HFIDMC_errIO_MAP(MyErr,ERROR_MSG_HFIDMCREADKEYWORDFORMAT,Achar2Wstr(FIELDNAME))

#define	HFIDMC_READKEYWORDTAB(STORAGE,DATANAME,FIELDNAME,OWNER) if((MyErr = PIOReadKeywordObject(((void *) &(DATANAME) ),dummyStr,FIELDNAME,STORAGE,OWNER,TAB2Dgrp_)) != 0) \
																	HFIDMC_errIO_TAB2D(MyErr,ERROR_MSG_HFIDMCREADKEYWORDFORMAT,FIELDNAME)

#define	HFIDMC_WRITEKEYWORD(STORAGE,DATANAME,FIELDNAME,COMMENT,OWNER) if((MyErr = PIOWriteKeywordObject(((void *) &(DATANAME) ),COMMENT,FIELDNAME,STORAGE,OWNER,VECTgrp_)) != 0) \
																		HFIDMC_errIO_Vect(MyErr,ERROR_MSG_HFIDMCWRITEKEYWORDFORMAT,FIELDNAME)

#define	HFIDMC_WRITEKEYWORDMAP(STORAGE,DATANAME,FIELDNAME,COMMENT,OWNER) if((MyErr = PIOWriteKeywordObject(((void *) &(DATANAME) ),COMMENT,FIELDNAME,STORAGE,OWNER,MAPgrp_)) != 0) \
																			HFIDMC_errIO_MAP(MyErr,ERROR_MSG_HFIDMCWRITEKEYWORDFORMAT,Achar2Wstr(FIELDNAME))

#define	HFIDMC_WRITEKEYWORDTAB(STORAGE,DATANAME,FIELDNAME,COMMENT,OWNER) if((MyErr = PIOWriteKeywordObject(((void *) &(DATANAME) ),COMMENT,FIELDNAME,STORAGE,OWNER,TAB2Dgrp_)) != 0) \
																			HFIDMC_errIO_TAB2D(MyErr,ERROR_MSG_HFIDMCWRITEKEYWORDFORMAT,Achar2Wstr(FIELDNAME))

#define	HFIDMC_WRITEKEYWORDIMG2D(STORAGE,DATANAME,FIELDNAME,COMMENT,OWNER) if((MyErr = PIOWriteKeywordObject(((void *) &(DATANAME) ),COMMENT,FIELDNAME,STORAGE,OWNER,IMGgrp_)) != 0) \
																			HFIDMC_errIO_IMG(MyErr,ERROR_MSG_HFIDMCWRITEKEYWORDFORMAT,Achar2Wstr(FIELDNAME))

#define	HFIDMC_READKEYWORDIMG2D(STORAGE,DATANAME,FIELDNAME,OWNER) if((MyErr = PIOReadKeywordObject(((void *) &(DATANAME) ),dummyStr,FIELDNAME,STORAGE,OWNER,IMGgrp_)) != 0) \
																	HFIDMC_errIO_IMG(MyErr,ERROR_MSG_HFIDMCREADKEYWORDFORMAT,FIELDNAME)

#define	LFIDPC_INSERTVALUE_AUX(argument,type) if(Params_->param_present(MAKECHARARR(argument))) \
	{Insert1element(std::wstring(MAKEWCHARARR(argument)),DBField(Params_->find<type>(std::string(MAKECHARARR(argument)))));}

#define	LFIDPC_INSERTVALUE(argument,type) LFIDPC_INSERTVALUE_AUX(argument,type)

#define	LFIDPC_INSERTVALUESTR_AUX(argument) if(Params_->param_present(MAKECHARARR(argument))) \
	{Insert1element(std::wstring(MAKEWCHARARR(argument)),DBField(Zeus::Achar2Wstr((Params_->find<std::string>(std::string(MAKECHARARR(argument)))).c_str())));}

#define	LFIDPC_INSERTVALUESTR(argument) LFIDPC_INSERTVALUESTR_AUX(argument)

#define CREATVECTOBJ(OBJECTNAME,TYPE) if((MyErr = PIOCreateVECTObject(OBJECTNAME,TYPE,VECTgrp_))) \
										HFIDMC_errIO_Vect(MyErr,ERROR_MSG_HFIDMCCREATEVECTFORMAT,OBJECTNAME)


#define GEOFILEWRITEVECTINT(OBJECTNAME,FIELDID) for(piv =data.Storage_.begin(),vectINTpiv=vectINT;piv != end;++piv,++vectINTpiv) \
													{*vectINTpiv = piv->FIELDID;} \
												if((MySZ=PIOWriteVECTObject(vectINT,OBJECTNAME,"PIOINT",command,VECTgrp_))<0) \
													HFIDMC_errIO_Vect((PIOErr) MySZ,ERROR_MSG_HFIDMCWRITEVECTFORMAT,OBJECTNAME)

#define GEOFILEWRITEVECTFLOAT(OBJECTNAME,FIELDID) for(piv =data.Storage_.begin(),vectFLOATpiv=vectFLOAT;piv != end;++piv,++vectFLOATpiv) \
													{*vectFLOATpiv = static_cast<double>(piv->FIELDID);} \
													if((MySZ=PIOWriteVECTObject(vectFLOAT,OBJECTNAME,"PIODOUBLE",command,VECTgrp_))<0) \
													HFIDMC_errIO_Vect((PIOErr) MySZ,ERROR_MSG_HFIDMCWRITEVECTFORMAT,OBJECTNAME) 

#define INTFILEWRITEVECTINT(OBJECTNAME,FIELDID) for(piv =data.begin(),vectINTpiv=vectINT;piv != end;++piv,++vectINTpiv) \
													{*vectINTpiv = piv->FIELDID;} \
												if((MySZ=PIOWriteVECTObject(vectINT,OBJECTNAME,"PIOINT",command,VECTgrp_))<0) \
													HFIDMC_errIO_Vect((PIOErr) MySZ,ERROR_MSG_HFIDMCWRITEVECTFORMAT,OBJECTNAME)

#define INTFILEWRITEVECTFLOAT(OBJECTNAME,FIELDID) for(piv =data.begin(),vectFLOATpiv=vectFLOAT;piv != end;++piv,++vectFLOATpiv) \
													{*vectFLOATpiv = static_cast<double>(piv->FIELDID);} \
													if((MySZ=PIOWriteVECTObject(vectFLOAT,OBJECTNAME,"PIODOUBLE",command,VECTgrp_))<0) \
													HFIDMC_errIO_Vect((PIOErr) MySZ,ERROR_MSG_HFIDMCWRITEVECTFORMAT,OBJECTNAME) 

#define INTFILEREADDELVECTINT(OBJECTNAME,FIELDID) MyErr=PIOReadVECTObject((void **)&(tVectINT),OBJECTNAME,"PIOINT",command,VECTgrp_); \
													if (MyErr<0) HFIDMC_errIO_Vect((PIOErr) MyErr,ERROR_MSG_HFIDMCREADVECTFORMAT,OBJECTNAME); \
													for(piv = pivOrg,tVectINTptr = tVectINT;piv != End;++tVectINTptr,++piv) \
														{piv->FIELDID = *tVectINTptr;} \
													MyErr=PIODeleteVECTTable(tVectINT,VECTgrp_); \
													if (MyErr != 0 ) HFIDMC_ReportMemoryLeak((PIOErr) MyErr,ERROR_MSG_HFIDMCDELETEVECTFORMAT,OBJECTNAME)

#define INTFILEREADDELVECTFLOAT(OBJECTNAME,FIELDID) MyErr=PIOReadVECTObject((void **)&(tVectFLOAT),OBJECTNAME,"PIODOUBLE",command,VECTgrp_); \
													if (MyErr<0) HFIDMC_errIO_Vect((PIOErr) MyErr,ERROR_MSG_HFIDMCREADVECTFORMAT,OBJECTNAME); \
													for(piv = pivOrg,tVectFLOATptr = tVectFLOAT;piv != End;++tVectFLOATptr,++piv) \
														{piv->FIELDID = static_cast<double>(*tVectFLOATptr);} \
													MyErr=PIODeleteVECTTable(tVectFLOAT,VECTgrp_); \
													if (MyErr != 0 ) HFIDMC_ReportMemoryLeak((PIOErr) MyErr,ERROR_MSG_HFIDMCDELETEVECTFORMAT,OBJECTNAME)

#define INTFILEAPPENDVECTDOUBLELFI(OBJECTNAME,FIELDID) for(piv = data.begin(),VectDOUBLEptr = VectDOUBLE.begin();piv!=end;++piv,++VectDOUBLEptr) \
														{*VectDOUBLEptr = piv->FIELDID;} \
														CollHandle_->appendColumn(std::string(OBJECTNAME),VectDOUBLE)

#define INTFILEAPPENDVECTINTLFI(OBJECTNAME,FIELDID) for(piv = data.begin(),VectINTptr = VectINT.begin();piv!=end;++piv,++VectINTptr) \
													 {*VectINTptr = piv->FIELDID;} \
													 CollHandle_->appendColumn(std::string(OBJECTNAME),VectINT)

#define INTFILEREADVECTINTLFI(OBJECTNAME,FIELDID)	CollHandle_->readEntireColumn(std::string(OBJECTNAME),tArrINT); \
													for(tArrINTptr = tArrINT.begin(),piv = pivOrg;piv != End;++piv,++tArrINTptr) \
														piv->FIELDID = *tArrINTptr


#define INTFILEREADVECTFLOATLFI(OBJECTNAME,FIELDID)	CollHandle_->readEntireColumn(std::string(OBJECTNAME),tArrFLOAT); \
													for(tArrFLOATptr = tArrFLOAT.begin(),piv = pivOrg;piv != End;++piv,++tArrFLOATptr) \
														piv->FIELDID = static_cast<double>(*tArrFLOATptr)


#define INTFILEWRITEVECTDOUBLELFI(OBJECTNAME,FIELDID) for(piv = data.Storage_.begin(),VectDOUBLEptr = VectDOUBLE.begin();piv!=end;++piv,++VectDOUBLEptr) \
														{*VectDOUBLEptr = piv->FIELDID;} \
														CollHandle_->writeColumn(std::string(OBJECTNAME),VectDOUBLE)

#define INTFILEWRITEVECTINTLFI(OBJECTNAME,FIELDID) for(piv = data.Storage_.begin(),VectINTptr = VectINT.begin();piv!=end;++piv,++VectINTptr) \
													 {*VectINTptr = piv->FIELDID;} \
													 CollHandle_->writeColumn(std::string(OBJECTNAME),VectINT)


#ifdef WIN32
#define	PRId64 "ld"
#else
#define PRId64 "ld"
#endif
//
namespace Zeus
{
//
#ifdef	HFIDMC
//

int	HFIDMC_RemoveObject(const std::wstring& Dir,const std::wstring& file);

template<typename T>
class	WriteHealpixMapHFIDMC: public GenHealpixWriter<T>
{
protected:
	virtual  std::wstring	do_GetCollID(void) const
	{return DirName_ + std::wstring(L"/") + FName_;}
// 
	virtual bool	do_Initialize(void)
	{
		PIOErr			MyErr;

		Dispose();

		if (!(MAPgrp_ = PIOOpenMAPGrp(const_cast<char *>(Wstr2Str(DirName_).c_str()),"w")))
		{HFIDMC_errIO(ERROR_COD_HFIDMCERROPENFL,ERROR_MSG_HFIDMCERROPENFL,DirName_);}

		if(MyErr = PIOCreateMAPObject(const_cast<char *>(Wstr2Str(FName_).c_str()),Private::PIOTypeChooserType<T>::PIO_TYPE,MAPgrp_))
		{HFIDMC_errIO_MAP(MyErr,ERROR_MSG_HFIDMCIMGMAPCREATEFORMAT,FName_);}

		return true;
	}
//
	virtual int	do_Write(const Healpix_Map<T>& data,HealpixHeaderType * hd)
	{
		if(hd)
		{CreateHeader(*hd);}

		PIOLONG			MyErr;
		PIOSTRING		command;
		std::string		tMapStr(Wstr2Str(FName_));
		
		sprintf(command,"begin=%" PRId64 ";end=%" PRId64,(PIOLONG) 0,static_cast<PIOLONG>(data.Npix() -1));

		if((MyErr = PIOWriteMAPObject((void *)&(data[0]),const_cast<char *>(tMapStr.c_str()),
							Private::PIOTypeChooserType<T>::PIO_TYPE,command,MAPgrp_)) < 0)
		{HFIDMC_errIO_MAP((PIOErr) MyErr,ERROR_MSG_HFIDMCIMGMAPWRITEFORMAT,FName_);}

		return true;
	}
//
public:
	WriteHealpixMapHFIDMC(const std::wstring& DirName,const std::wstring& fname,
		int NSide,Healpix_Ordering_Scheme Scheme,coordsys CoordSys)
		:GenHealpixWriter<T>(),DirName_(DirName),FName_(fname),Scheme_(Scheme),CoordSys_(CoordSys),MAPgrp_(0)
	{
		PIOErr	MyErr;
		std::string	 tDirName(Wstr2Str(DirName_));
		int retries(MAXRETRIES);

_retry:
		if(PIOCheckGroup(const_cast<char *>(tDirName.c_str())))
		{
			if((MyErr = PIOCreateMAPGrp(const_cast<char *>(tDirName.c_str()),
				const_cast<char *>(GetStrFromCoordsys(CoordSys_).c_str()),
				const_cast<char *>(GetStrFromOrdering(Scheme).c_str()),
				NSide)) != 0)
			{
				if(--retries >= 0)
				{
					MySleep(1);
					goto _retry;
				}
				HFIDMC_errIO(ERROR_COD_HFIDMCERRCREATFL,ERROR_MSG_HFIDMCERRCREATFL,DirName_);
			}
		}
	}
	virtual ~WriteHealpixMapHFIDMC(void)
	{Dispose();}
private:
	std::wstring				DirName_;
	std::wstring				FName_;
	Healpix_Ordering_Scheme		Scheme_;
	coordsys					CoordSys_;
	PIOGroup					*MAPgrp_;

	inline void			Dispose(void)
	{
		if(MAPgrp_)
		{PIOCloseMAPGrp(&MAPgrp_);MAPgrp_ = 0;}	
	}

	inline void			HFIDMC_errIO(int errCode,wchar_t* msg,const std::wstring& xtraName) const
	{
		std::wstring errstring(msg);

		errstring += std::wstring(L" -> ");
		errstring += xtraName;
		throw Zeus::libException(errCode,errstring,*this);
	}

	inline void			HFIDMC_errIO_MAP(PIOErr MyErr,char* formatStr,const std::wstring& fieldName) const
	{
		char buff[STRBUFFER];
		sprintf(buff,formatStr,Wstr2Str(fieldName).c_str(),PIOErrMess(MyErr));
		throw Zeus::libException(MyErr,Achar2Wstr(buff),*this);
	}
	
	void	CreateHeader(const HealpixHeaderType& hd)
	{
		PIOErr			MyErr;
		std::string		tFName(Wstr2Str(FName_));
		PIOSTRING		Comment;
		const char*		StrBuff;

		HFIDMC_WRITEKEYWORDMAP("PIOSTRING",Scheme_,HFIDMC_HEALPIXMAP_SCHEME,HFIDMC_HEALPIXMAP_SCHEMEMSG,const_cast<char *>(tFName.c_str()));
		HFIDMC_WRITEKEYWORDMAP("PIOSTRING",CoordSys_,HFIDMC_HEALPIXMAP_COORSYS,HFIDMC_HEALPIXMAP_COORSYSMSG,const_cast<char *>(tFName.c_str()));

		if(!(hd.Columns_.empty()))
		{
			HealpixHeaderType::HealpixHeaderItemCollType::const_iterator	piv(hd.Columns_.begin());
			HealpixHeaderType::HealpixHeaderItemCollType::const_iterator	const end(hd.Columns_.end());
			for(;piv != end;++piv)
			{
				switch(piv->ValueType_)
				{
				case INT_TYPE:
					HFIDMC_WRITEKEYWORDMAP("PIOINT",piv->PODValue_.intType_,const_cast<char *>((piv->Keyword_).c_str()),Comment,const_cast<char *>(tFName.c_str()));
					break;
				case FLOAT_TYPE:
					HFIDMC_WRITEKEYWORDMAP("PIOFLOAT",piv->PODValue_.floatType_,const_cast<char *>((piv->Keyword_).c_str()),Comment,const_cast<char *>(tFName.c_str()));
					break;
				case DOUBLE_TYPE:
					HFIDMC_WRITEKEYWORDMAP("PIODOUBLE",piv->PODValue_.doubleType_,const_cast<char *>((piv->Keyword_).c_str()),Comment,const_cast<char *>(tFName.c_str()));
					break;
				case BOOL_TYPE:
					HFIDMC_WRITEKEYWORDMAP("PIOINT",piv->PODValue_.boolType_,const_cast<char *>((piv->Keyword_).c_str()),Comment,const_cast<char *>(tFName.c_str()));
					break;
				case STRING_TYPE:
					StrBuff = piv->stringValue_.c_str();
					HFIDMC_WRITEKEYWORDMAP("PIOSTRING",StrBuff[0],const_cast<char *>((piv->Keyword_).c_str()),Comment,const_cast<char *>(tFName.c_str()));
					break;
				default:
					break;
				}
			}

		}
	}
};
//
template<typename T>
class	ReadHealpixMapHFIDMC: public GenHealpixReader<T>
{
protected:
//
	virtual std::wstring	do_GetCollID(void) const
	{return  DirName_ + std::wstring(L"/") + MapName_;}
//
	virtual bool	do_Initialize(HealpixHeaderType *hd)
	{
		PIOErr			MyErr;
		PIOSTRING		DMC_helpStr;

		MAPgrp_ = PIOOpenMAPGrp(const_cast<char *>(Wstr2Str(DirName_).c_str()),"r");
		if (!MAPgrp_)
		{HFIDMC_errIO(ERROR_COD_HFIDMCERROPENHEALPX,ERROR_MSG_HFIDMCERROPENHEALPX,DirName_);}
		PIOLONG	tNSide;
		if((tNSide = PIONSideGrp(MAPgrp_)) < 0)
		{HFIDMC_errIO_MAP((PIOErr) tNSide,ERROR_MSG_HFIDMCIMGPROPSFORMAT,MapName_);}
		NSide_		= static_cast<int>(tNSide);
		if ((MyErr = PIOOrderingGrp(DMC_helpStr,MAPgrp_)) < 0)
		{HFIDMC_errIO_MAP(MyErr,ERROR_MSG_HFIDMCIMGPROPSFORMAT,MapName_);}
		Ordering_	= GetOrderingFromStr(std::string(DMC_helpStr));

		if(!hd) return true;

		hd->NSide_ = NSide_;

		if(hd->Columns_.empty())
		{QueryKeys(*hd);}
		else
		{ReadKeys(*hd);}
		return true;
	}

//
	virtual void	do_Release(void)
	{do_DisposeData();}
//
	virtual bool	do_Read(void)
	{
		PIOLONG			MyErr;
		PIOSTRING		command;
		PIOFLOAT		*map_vect;
		std::string		tMapStr(Wstr2Str(MapName_));

		arr<T> myarr(12 * NSide_ * NSide_);
		
		sprintf(command,"begin=%" PRId64 ";end=%" PRId64,(PIOLONG) 0,static_cast<PIOLONG>(myarr.size() -1));

		if((MyErr=PIOReadMAPObject((void **)&map_vect,const_cast<char *>(tMapStr.c_str()),Private::PIOTypeChooserType<T>::PIO_TYPE,command,MAPgrp_)) < 0)
			HFIDMC_errIO_MAP((PIOErr) MyErr,ERROR_MSG_HFIDMCIMGMAPREADFORMAT,MapName_);

		PIOFLOAT		*map_vectPtr(map_vect);
		T				*pivMap(myarr.begin());
		T				* const EndMap(myarr.end());

		for(;pivMap != EndMap;++map_vectPtr,++pivMap)
		{*pivMap = static_cast<T>(*map_vectPtr);}

		if ((MyErr=PIODeleteMAPTable(&map_vect,MAPgrp_)) != 0 )
			HFIDMC_ReportMemoryLeak((PIOErr) MyErr,ERROR_MSG_HFIDMCDELETEIMGFORMAT,tMapStr.c_str());

		GenHealpixReader<T>::HPixMap_->Set(myarr,Ordering_);
		return true;
	}
//
public:
	ReadHealpixMapHFIDMC(const std::wstring& DirName,const std::wstring& MapName)
		:GenHealpixReader<T>(),DirName_(DirName),MapName_(MapName),MAPgrp_(0)
	{}
//
	virtual ~ReadHealpixMapHFIDMC(void)
	{do_DisposeData();}
//
private:
	std::wstring				MapName_;
	std::wstring				DirName_;
	PIOGroup					*MAPgrp_;
	int							NSide_;
	Healpix_Ordering_Scheme		Ordering_;
//
	inline void			HFIDMC_errIO(int errCode,wchar_t* msg,const std::wstring& xtraName) const
	{
		std::wstring errstring(msg);

		errstring += std::wstring(L" -> ");
		errstring += xtraName;
		throw Zeus::libException(errCode,errstring,*this);
	}
//
	inline void			HFIDMC_errIO_MAP(PIOErr MyErr,const char* formatStr,const std::wstring& fieldName) const
	{
		char buff[STRBUFFER];
		sprintf(buff,formatStr,Zeus::Wstr2Str(fieldName).c_str(),PIOErrMess(MyErr));
		throw Zeus::libException(MyErr,Achar2Wstr(buff),*this);
	}
//
	inline void	do_DisposeData(void)
	{
		if(MAPgrp_)
		{PIOCloseMAPGrp(&MAPgrp_);MAPgrp_ = 0;}	
	}
//
	void	QueryKeys(HealpixHeaderType& hd)
	{
		PIOErr						MyErr;
		PIOLONG						NKeys;
		PIOSTRING 					*Keywords(0);
		PIOSTRING 					*Types(0);
		PIOSTRING					dummyStr;
		std::string					tMapName(Wstr2Str(MapName_));
		HFIDMC_ReadKeyBufType		Buffer;

		if((NKeys = PIOKeywordListObject(&Keywords,&Types,const_cast<char *>(tMapName.c_str()),MAPgrp_))<0)
		{HFIDMC_errIO_MAP((PIOErr) NKeys,ERROR_MSG_HFIDMCIMGPROPSFORMAT,MapName_);}
	#undef __PIOFUNCT__ 
	#define __PIOFUNCT__ "PIOReadKeywordObject"

		for(PIOLONG i=0;i<NKeys;++i)
		{
			if((MyErr = PIOReadKeywordObject(((void *)&(Buffer)),dummyStr,Keywords[i],Types[i],const_cast<char *>(tMapName.c_str()),MAPgrp_)) != 0)
			{HFIDMC_errIO_MAP(MyErr,ERROR_MSG_HFIDMCREADKEYWORDFORMAT,Achar2Wstr(Keywords[i]));}
			hd.Columns_.push_back(HealpixHeaderAtomType(std::string(Keywords[i]),GetPIOType(Types[i]),(void *)&Buffer));
		}
	#undef __PIOFUNCT__ 
	#define __PIOFUNCT__ "PIOKeywordListObject"
		if(Keywords)
		{_PIOFREE(Keywords);}
		if(Types)
		{_PIOFREE(Types);}
	}
//
	void	ReadKeys(HealpixHeaderType& hd)
	{
		PIOSTRING					dummyStr;

		std::string					tMapName(Wstr2Str(MapName_));

		HealpixHeaderType::HealpixHeaderItemCollType::iterator	piv(hd.Columns_.begin());
		HealpixHeaderType::HealpixHeaderItemCollType::const_iterator	const end(hd.Columns_.end());

		for(;piv!=end;++piv)
		{
			switch(piv->ValueType_)
			{
			case Zeus::INT_TYPE:
				if(PIOReadKeywordObject((void *) &(piv->PODValue_.intType_),dummyStr,const_cast<char *>(piv->Keyword_.c_str()),
					Private::PIOTypeChooserType<int>::PIO_TYPE,const_cast<char *>(tMapName.c_str()),MAPgrp_) != 0)
				{piv->Reset();}
				break;
			case Zeus::FLOAT_TYPE:
				if(PIOReadKeywordObject((void *) &(piv->PODValue_.floatType_),dummyStr,const_cast<char *>(piv->Keyword_.c_str()),
					Private::PIOTypeChooserType<float>::PIO_TYPE,const_cast<char *>(tMapName.c_str()),MAPgrp_) != 0)
				{piv->Reset();}
				break;
			case Zeus::DOUBLE_TYPE:
				if(PIOReadKeywordObject((void *) &(piv->PODValue_.doubleType_),dummyStr,const_cast<char *>(piv->Keyword_.c_str()),
					Private::PIOTypeChooserType<double>::PIO_TYPE,const_cast<char *>(tMapName.c_str()),MAPgrp_) != 0)
				{piv->Reset();}
				break;
			case Zeus::BOOL_TYPE:
				if(PIOReadKeywordObject((void *) &(piv->PODValue_.boolType_),dummyStr,const_cast<char *>(piv->Keyword_.c_str()),
					Private::PIOTypeChooserType<bool>::PIO_TYPE,const_cast<char *>(tMapName.c_str()),MAPgrp_) != 0)
				{piv->Reset();}
				break;
			case Zeus::STRING_TYPE:
				PIOSTRING	buffStr;

				if(PIOReadKeywordObject((void *) buffStr,dummyStr,const_cast<char *>(piv->Keyword_.c_str()),
					Private::PIOTypeChooserType<bool>::PIO_TYPE,const_cast<char *>(tMapName.c_str()),MAPgrp_) != 0)
				{piv->Reset();}
				else{piv->Set(std::string(buffStr));}
				break;
			default:
				piv->Reset();
				break;
			}
		}
	}
//
};

//
class	HFIDMC_ParameterFileReader: public GenCollReader<ParamVarsStoreType>
{
private:
	std::wstring			FName_;
	ParamVarsStoreType		*data_;
//
	inline bool				Insert1element(const std::wstring& id,const DBField& value)
	{
//		(Zeus::ConManager::Instance())->PrintStr2Console(std::wstring(L"DMC Key -> ") + id);

		data_->Storage_.erase(id);
		return data_->Storage_.insert(ParamVarsStoreType::StorageType::value_type(id, value)).second;
	}
//
	int						DMC_GetValue(const std::wstring& in,DBField& fld);
	void					ReadValuesFromDMC_Pipeline(void);
	void					InsertFinalCatStr(void);
	void					InsertOutputCatStr(void);

protected:
//
	virtual  std::wstring	do_GetCollID(void) const
	{return FName_;}
//
	virtual bool	do_Read(void)
	{
		ReadValuesFromDMC_Pipeline();
		return true;		
	}
//
	virtual	ParamVarsStoreType&	do_GetData(void)
	{return *data_;}
//
	virtual void	do_DisposeData(void)
	{delete data_;data_ = 0;}
//
	virtual ParamVarsStoreType::HeaderType *do_Initialize(void)
	{
#undef __PIOFUNCTION__
#define __PIOFUNCTION__ "HFIDMC_ParameterFileReader:Initialise"

		PIODBObjectLog		*logger((ConManager::Instance())->GetLogger());
		parContent			*t(init_parContent(const_cast<char *>(Wstr2Str(FName_).c_str()),logger));
		if (t == NULL)
		{
			_PIOEXITCAUSE_(logger,"PIE problem, bad parameters returned\n",-1);
			return 0;
		}
		data_	= new ParamVarsStoreType();
		data_->Header_ = HFIDMC_parContent_sptr(new HFIDMC_parContent_helper());
		data_->Header_->SetParContent(t);
		return &(data_->Header_);
	} 
public:
	HFIDMC_ParameterFileReader(const std::wstring& fname)
		: GenCollReader<ParamVarsStoreType>(),FName_(fname),data_(0)
	{}
//
	virtual ~HFIDMC_ParameterFileReader(void)
	{
		do_DisposeData();
	}
};
//
template<typename T>
class		WriteWrkSpaceHFIDMC: public GenCollWriter<LArr2D<T> >
{
	std::wstring	FName_;
	std::wstring	GrpName_;
	PIOGroup		*IMGgrp_;
	PIODOUBLE		PatchCx_;  //YMaxValue_
	PIODOUBLE		PatchCy_;  //YMinValue_	
	PIODOUBLE		XMaxValue_;
	PIODOUBLE		XMinValue_;
	PIOLONG			AxSz_;
	PIOFLOAT		PixSz_;

	inline void	DisposeData(void)
	{
		if(IMGgrp_)
		{PIOCloseIMG2DGrp(&IMGgrp_);IMGgrp_ = 0;}	
	}

	inline void			HFIDMC_errIO(int errCode,wchar_t* msg,const std::wstring& xtraName) const
	{
		std::wstring errstring(msg);

		errstring += std::wstring(L" -> ");
		errstring += xtraName;
		throw Zeus::libException(errCode,errstring,*this);
	}
//
	inline void			HFIDMC_errIO_IMG(PIOErr MyErr,char* formatStr,const std::wstring& fieldName) const
	{
		char buff[STRBUFFER];
		sprintf(buff,formatStr,Wstr2Str(fieldName).c_str(),PIOErrMess(MyErr));
		throw Zeus::libException(MyErr,Achar2Wstr(buff),*this);
	}
//
protected:
	inline virtual std::wstring	do_GetCollID(void) const
	{return GrpName_  + std::wstring(L"/") + FName_;}

	inline virtual int			do_Flush(void)
	{return 0;}
//
	inline virtual void			do_DisposeData(void)
	{DisposeData();return ;}
//
	inline virtual int		do_Remove(void)
	{return HFIDMC_RemoveObject(GrpName_,FName_);}

	inline virtual bool	do_Initialize(void)
	{
		PIOErr			MyErr;
		std::string		tFname(Wstr2Str(FName_));
		std::string		tGrpName(Wstr2Str(GrpName_));

		// Group *MUST* exist

		if (!(IMGgrp_ = PIOOpenIMG2DGrp(const_cast<char *>(tGrpName.c_str()),"w")))
		{HFIDMC_errIO(ERROR_COD_HFIDMCERROPENFL,ERROR_MSG_HFIDMCERROPENFL,GrpName_);}

		if(MyErr = PIOCreateIMG2DObject(const_cast<char *>(tFname.c_str()),Private::PIOTypeChooserType<T>::PIO_TYPE,PatchCx_,PatchCy_,IMGgrp_))
		{HFIDMC_errIO_IMG(MyErr,ERROR_MSG_HFIDMCCREATEIMGFORMAT,FName_);}

		HFIDMC_WRITEKEYWORDIMG2D("PIODOUBLE",PatchCx_,"MAXCY5R500","Max Y5R500",const_cast<char *>(tFname.c_str()));
		HFIDMC_WRITEKEYWORDIMG2D("PIODOUBLE",PatchCy_,"MINCY5R500","Min Y5R500",const_cast<char *>(tFname.c_str()));
		HFIDMC_WRITEKEYWORDIMG2D("PIODOUBLE",XMaxValue_,"MAXRS","Max RS",const_cast<char *>(tFname.c_str()));
		HFIDMC_WRITEKEYWORDIMG2D("PIODOUBLE",XMinValue_,"MINRS","Min RS",const_cast<char *>(tFname.c_str()));

		return true;
	}

	virtual int	do_Write(const LArr2D<T> & data)
	{
		PIOLONG		MySz;
		PIOSTRING	command;
		PIOLONG		ImgLastIndex(data.getPtrMetric() - 1);
		
		sprintf(command,"tab=0:%" PRId64 ",0:%" PRId64,(PIOLONG)ImgLastIndex,(PIOLONG)ImgLastIndex );
		
		if((MySz=PIOWriteIMG2DObject((void *)data.begin(),const_cast<char *>(Wstr2Str(FName_).c_str()),Private::PIOTypeChooserType<T>::PIO_TYPE,command,IMGgrp_))<0)
		{HFIDMC_errIO_IMG((PIOErr) MySz,ERROR_MSG_HFIDMCWRITEIMGFORMAT,FName_);}

		return ImgLastIndex + 1;
	}

public:
	WriteWrkSpaceHFIDMC(const std::wstring& GrpName,const std::wstring& fname,PIOLONG AxSz,PIOFLOAT PixSz,PIODOUBLE PatchCx,PIODOUBLE PatchCy,PIODOUBLE XMaxValue,PIODOUBLE XMinValue)
		:GenCollWriter<LArr2D<T> >(),GrpName_(GrpName),FName_(fname),IMGgrp_(0),AxSz_(AxSz),PixSz_(PixSz),PatchCx_(PatchCx),PatchCy_(PatchCy),XMaxValue_(XMaxValue),XMinValue_(XMinValue)
	{}
//
	virtual ~WriteWrkSpaceHFIDMC(void)
	{DisposeData();}

};

template<typename T>
class	ReadWrkSpaceHFIDMC:  public GenCollReader<LArr2D<T> >
{
private:
	std::wstring	FName_;
	std::wstring	GrpName_;
	LArr2D<T>		*data_;
	int				YSz_;
	int				Metric_;
	PIOGroup		*IMGgrp_;
//
	inline void			HFIDMC_errIO(int errCode,wchar_t* msg,const std::wstring& xtraName) const
	{
		std::wstring errstring(msg);
		errstring += std::wstring(L" -> ");
		errstring += xtraName;
		throw Zeus::libException(errCode,errstring,*this);
	}
//
	inline void			HFIDMC_errIO_IMG(PIOErr MyErr,char* formatStr,const char* fieldName) const
	{
		char buff[STRBUFFER];
		sprintf(buff,formatStr,fieldName,PIOErrMess(MyErr));
		throw Zeus::libException(MyErr,Achar2Wstr(buff),*this);
	}

protected:
	virtual  std::wstring	do_GetCollID(void) const
	{return GrpName_  + std::wstring(L"/") + FName_;}
//
	virtual bool	do_Read(void)
	{
		PIOErr				MyErr;
		PIOLONG				NData;
		T					*tData;
		std::string			tObjName(Wstr2Str(FName_));
		PIOSTRING			command;
		PIOLONG				ImgLastIndex(data_->getPtrMetric() - 1);
		
		sprintf(command,"tab=0:%" PRId64 ",0:%" PRId64,(PIOLONG)ImgLastIndex,(PIOLONG)ImgLastIndex );

		if((NData = PIOReadIMG2DObject((void **)&tData,const_cast<char *>(tObjName.c_str()),Private::PIOTypeChooserType<T>::PIO_TYPE,command,IMGgrp_))< 0)
		{HFIDMC_errIO_IMG((PIOErr) NData,ERROR_MSG_HFIDMCREADIMGFORMAT,tObjName.c_str());}

		if(NData != data_->getSz())
		{
			if((MyErr = PIODeleteIMG2DTable(tData,IMGgrp_))!=0)
				HFIDMC_ReportMemoryLeak(MyErr,ERROR_MSG_HFIDMCDELETEIMGFORMAT,tObjName.c_str());
			HFIDMC_errIO(ERROR_COD_HFIDMCERRIMGWRONGSZ,ERROR_MSG_HFIDMCERRIMGWRONGSZ,FName_);
		}

		T			*tDestPtr(data_->begin());
		T			*tOrgPtr(tData);

		for(PIOLONG i=0;i<NData;++i,++tOrgPtr,++tDestPtr)
		{*tDestPtr = static_cast<T>(*tOrgPtr);}
		

		if((MyErr = PIODeleteIMG2DTable(tData,IMGgrp_))!=0)
			HFIDMC_ReportMemoryLeak(MyErr,ERROR_MSG_HFIDMCDELETEIMGFORMAT,tObjName.c_str());
		return true;
	}
//
	virtual	LArr2D<T>&	do_GetData(void)
	{return *data_;}
//
	virtual void	do_DisposeData(void)
	{
		delete data_;data_ = 0;
		
		if(IMGgrp_)
		{PIOCloseIMG2DGrp(&IMGgrp_);IMGgrp_ = 0;}
	}
//
	virtual typename LArr2D<T>::HeaderType *do_Initialize(void)
	{
		do_DisposeData();
		IMGgrp_ = PIOOpenIMG2DGrp(const_cast<char *>(Zeus::Wstr2Str(GrpName_).c_str()),"r");
		if (!IMGgrp_)
		{HFIDMC_errIO(ERROR_COD_HFIDMCERROPENIMGFL,ERROR_MSG_HFIDMCERROPENIMGFL,GrpName_);}
		data_ = new LArr2D<T>(YSz_ * Metric_ ,Metric_);
		data_->reset();
		return reinterpret_cast<typename LArr2D<T>::HeaderType *>(data_->begin());
	} 
public:
	ReadWrkSpaceHFIDMC(const std::wstring& GrpName,const std::wstring& fname,int YSz,int XSz)
		: GenCollReader<LArr2D<T> >(),GrpName_(GrpName),FName_(fname),YSz_(YSz),Metric_(XSz),data_(0),IMGgrp_(0)
	{}
//
	virtual ~ReadWrkSpaceHFIDMC(void)
	{
		do_DisposeData();
	}
};

//
class	ReadPatchesGeomInfoHFIDMC: public GenCollReader<PatchGeomType>
{
private:
	PIOLONG				NColumns_;
	PIOLONG				NRows_;
	PatchGeomType		data_;
	std::wstring		FName_;
	std::wstring		GrpName_;
	PIOGroup			*TAB2Dgrp_;

	int					ReadHeader(void);
	int					ReadBody(void);
//
	inline void			HFIDMC_errIO(int errCode,wchar_t* msg,const std::wstring& xtraName) const
	{
		std::wstring errstring(msg);

		errstring += std::wstring(L" -> ");
		errstring += xtraName;
		throw Zeus::libException(errCode,errstring,*this);
	}
//
	inline void			HFIDMC_errIO_TAB2D(PIOLONG MyErr,char* formatStr,char* fieldName) const
	{
		char buff[STRBUFFER];
		sprintf(buff,formatStr,fieldName,PIOErrMess((PIOErr) MyErr));
		throw Zeus::libException((int) MyErr,Achar2Wstr(buff),*this);
	}

protected:
	virtual  std::wstring	do_GetCollID(void) const
	{return GrpName_  + std::wstring(L"/") + FName_;}
//
	virtual bool	do_Read(void)
	{
		ReadBody();
		return true;		
	}
//
	virtual	PatchGeomType&	do_GetData(void)
	{return data_;}
//
	virtual void	do_DisposeData(void)
	{
		data_.Storage_.clear();
		if(TAB2Dgrp_)
		{PIOCloseTAB2DGrp(&TAB2Dgrp_);TAB2Dgrp_ = 0;}
	}
//
	virtual PatchGeomType::HeaderType *do_Initialize(void);

public:
	ReadPatchesGeomInfoHFIDMC(const std::wstring& GrpName,const std::wstring& fname)
		: GenCollReader<PatchGeomType>(),FName_(fname),GrpName_(GrpName + std::wstring(HFIDMCPACTHGEOGROUPEXTENSION)),
			NColumns_(0),NRows_(0),TAB2Dgrp_(0)
	{}
//
	virtual ~ReadPatchesGeomInfoHFIDMC(void)
	{do_DisposeData();}
};
//
class	WritePatchesGeomInfoHFIDMC: public GenCollWriter<PatchGeomType>
{
protected:
	virtual  std::wstring	do_GetCollID(void) const
	{return GrpName_  + std::wstring(L"/") + FName_;}

	inline virtual int			do_Flush(void)
	{return 0;}
//
	inline virtual void			do_DisposeData(void)
	{Dispose();return ;}
//
	inline virtual int		do_Remove(void)
	{return HFIDMC_RemoveObject(GrpName_,FName_);}

	inline virtual bool		do_Initialize(void)
	{
		PIOErr			MyErr;
		std::string		tGrpNameStr(Wstr2Str(GrpName_));
		int				retries(MAXRETRIES);

_retry:
		if(PIOCheckGroup(const_cast<char *>(tGrpNameStr.c_str())))
		{
			if(MyErr = PIOCreateTAB2DGrp(const_cast<char *>(tGrpNameStr.c_str()), NColumns_)) 
			{
				if(--retries >= 0)
				{
					MySleep(1);
					goto _retry;
				}
				HFIDMC_errIO(ERROR_COD_HFIDMCERRCREATFL,ERROR_MSG_HFIDMCERRCREATFL,GrpName_);
			}
		}

		if (!(TAB2Dgrp_ = PIOOpenTAB2DGrp(const_cast<char *>(tGrpNameStr.c_str()),"w")))
		{HFIDMC_errIO(ERROR_COD_HFIDMCERROPENFL,ERROR_MSG_HFIDMCERROPENFL,GrpName_);}
		
		if(MyErr = PIOCreateTAB2DObject(const_cast<char *>((Wstr2Str(FName_)).c_str()),"PIODOUBLE",TAB2Dgrp_)) 
		{HFIDMC_errIO_TAB2D(MyErr,ERROR_MSG_HFIDMCCREATEIMGFORMAT,FName_);}

		return true;
	}
//
	virtual int	do_Write(const PatchGeomType & data)
	{WriteGeomInfoData(data); return true;}
//
public:
	WritePatchesGeomInfoHFIDMC(const std::wstring& GrpName,const std::wstring& fname)
		:Zeus::GenCollWriter<PatchGeomType>(),GrpName_(GrpName + std::wstring(HFIDMCPACTHGEOGROUPEXTENSION)),FName_(fname),
		NColumns_(GEOINFONCOLUMNS),TAB2Dgrp_(0)
	{}
//
	virtual ~WritePatchesGeomInfoHFIDMC(void)
	{Dispose();}

private:
//
	inline void			HFIDMC_errIO(int errCode,wchar_t* msg,const std::wstring& xtraName) const
	{
		std::wstring errstring(msg);

		errstring += std::wstring(L" -> ");
		errstring += xtraName;
		throw Zeus::libException(errCode,errstring,*this);
	}
//
	inline void			HFIDMC_errIO_TAB2D(PIOErr MyErr,char* formatStr,const std::wstring& TableName) const
	{
		char buff[STRBUFFER];
		sprintf(buff,formatStr,const_cast<char *>((Wstr2Str(TableName)).c_str()),PIOErrMess(MyErr));
		throw Zeus::libException(MyErr,Achar2Wstr(buff),*this);
	}
//
	inline void			Dispose(void)
	{
		if(TAB2Dgrp_)
		{PIOCloseTAB2DGrp(&TAB2Dgrp_);}
		TAB2Dgrp_	= 0;
	}
//
	void	WriteGeomInfoData(const PatchGeomType & data);
//
	void	CreatStuff(void);
//
	void	CreatHeader(const PatchGeomType & data);

	std::wstring				FName_;
	std::wstring				GrpName_;
	PIOGroup					*TAB2Dgrp_;
	int							NColumns_;
};
//
class	WriteObjResultsHFIDMC: public GenCollWriter<PeakCollType>
{
protected:
	virtual  std::wstring	do_GetCollID(void) const
	{return GrpName_  + std::wstring(L"/") + FName_;}
//
	inline virtual int			do_Flush(void)
	{return Flush2DB();}
//
	inline virtual void			do_DisposeData(void)
	{Dispose();return ;}
//
	inline virtual int		do_Remove(void)
	{return HFIDMC_RemoveObject(GrpName_,FName_);}
//
	inline virtual bool	do_Initialize(void)
	{
		CreatStuff();
		return true;
	}
//
	inline virtual int	do_Write(const PeakCollType & data)
	{
		if(data.empty())	return 0;

		data_.insert(data_.end(),data.begin(),data.end());

		return (int)(data.end() - data.begin());
	}
public:
	WriteObjResultsHFIDMC(const std::wstring& GrpName,const std::wstring& fname,int InitObjID)
		:Zeus::GenCollWriter<PeakCollType>(),GrpName_(GrpName + std::wstring(HFIDMCINTERMCATGROUPEXTENSION)),FName_(fname),
		ObjectID_(InitObjID),NColumns_(INTERCATNCOLUMNS),TAB2Dgrp_(0)
	{}
//

	inline   std::wstring GetObjectID(void) const
	{return do_GetCollID();}
//
	virtual ~WriteObjResultsHFIDMC(void)
	{Dispose();}
private:
//
	inline void			HFIDMC_errIO(int errCode,wchar_t* msg,const std::wstring& xtraName) const
	{
		std::wstring errstring(msg);

		errstring += std::wstring(L" -> ");
		errstring += xtraName;
		throw Zeus::libException(errCode,errstring,*this);
	}
//
	inline void			HFIDMC_errIO_TAB2D(PIOErr MyErr,char* formatStr,const std::wstring& TableName) const
	{
		char buff[STRBUFFER];
		sprintf(buff,formatStr,const_cast<char *>((Wstr2Str(TableName)).c_str()),PIOErrMess(MyErr));
		throw Zeus::libException(MyErr,Achar2Wstr(buff),*this);
	}
//
	inline void			Dispose(void)
	{
		data_.clear();
		if(TAB2Dgrp_)
		{PIOCloseTAB2DGrp(&TAB2Dgrp_);}
		TAB2Dgrp_	= 0;
	}
//
	void				CreatStuff(void);
//
	void				CreatHeader(void);
//
	int					Flush2DB(void);
//
	int					ObjectID_;
	std::wstring		FName_;
	std::wstring		GrpName_;
	PIOGroup			*TAB2Dgrp_;
	int					NRows_;
	int					NColumns_;
	PeakCollType		data_;
};
//
class	CatalogueWriterHFIDMC: public GenCollWriter<CatalogueFormatType>
{
protected:
//
	inline virtual  std::wstring	do_GetCollID(void) const
	{return GrpName_  + std::wstring(L"/") + FName_;}
//
	inline virtual bool			do_Initialize(void)
	{
		CreatStuff();
		return true;
	}
//
	virtual int					do_Write(const CatalogueFormatType & dataIN);
//
	inline virtual int			do_Flush(void)
	{return 0;}
//
	inline virtual void			do_DisposeData(void)
	{Dispose();return ;}
//
	inline virtual int			do_Remove(void)
	{return HFIDMC_RemoveObject(GrpName_,FName_);}

public:
	CatalogueWriterHFIDMC(const std::wstring& GrpName,const std::wstring& fname,int NColumns)
		:GenCollWriter<CatalogueFormatType>(),GrpName_(GrpName + std::wstring(HFIDMCCATSGROUPEXTENSION))
		,FName_(fname),NColumns_(NColumns),TAB2Dgrp_(0)
	{}
//
	virtual ~CatalogueWriterHFIDMC(void)
	{Dispose();}
private:
//
	inline void			HFIDMC_errIO(int errCode,wchar_t* msg,const std::wstring& xtraName) const
	{
		std::wstring errstring(msg);

		errstring += std::wstring(L" -> ");
		errstring += xtraName;
		throw Zeus::libException(errCode,errstring,*this);
	}
//
	inline void			HFIDMC_errIO_TAB2D(PIOErr MyErr,char* formatStr,const std::wstring& TableName) const
	{
		char buff[STRBUFFER];
		sprintf(buff,formatStr,const_cast<char *>((Wstr2Str(TableName)).c_str()),PIOErrMess(MyErr));
		throw Zeus::libException(MyErr,Achar2Wstr(buff),*this);
	}

//
	inline void			Dispose(void)
	{
		if(TAB2Dgrp_)
		{PIOCloseTAB2DGrp(&TAB2Dgrp_);TAB2Dgrp_	= 0;}
	}
//
	void				CreatStuff(void);
//
	void				CreateHeader(const CatalogueFormatType & data);
//
	void				Move2Array2D(const CatalogueFormatType & DataIN,LArr2D<double> & DataOUT);
//
	PIOLONG						NColumns_;
	std::wstring				FName_;
	std::wstring				GrpName_;
	PIOGroup					*TAB2Dgrp_;
};
//
class	CatWriterOutputHFIDMC: public GenCollWriter<OutputFormatType>
{
protected:
//
	inline virtual  std::wstring	do_GetCollID(void) const
	{return GrpName_  + std::wstring(L"/") + FName_;}
//
	inline virtual bool			do_Initialize(void)
	{
		CreatStuff();
		return true;
	}
//
	virtual int					do_Write(const OutputFormatType & dataIN);
//
	inline virtual int			do_Flush(void)
	{return 0;}
//
	inline virtual void			do_DisposeData(void)
	{Dispose();return ;}
//
	inline virtual int			do_Remove(void)
	{return HFIDMC_RemoveObject(GrpName_,FName_);}
//
	void						CreateHeaderHelper(const OutputFormatType& data);
//
	virtual	void				CreateHeader(const OutputFormatType & data) = 0;
//
	virtual	void				Move2Array2D(const OutputFormatType & DataIN,LArr2D<double> & DataOUT) = 0;
//
	virtual	int					GetNColumns(void) const = 0;
//
	inline void			HFIDMC_errIO(int errCode,wchar_t* msg,const std::wstring& xtraName) const
	{
		std::wstring errstring(msg);

		errstring += std::wstring(L" -> ");
		errstring += xtraName;
		throw Zeus::libException(errCode,errstring,*this);
	}
//
	inline void			HFIDMC_errIO_TAB2D(PIOErr MyErr,char* formatStr,const std::wstring& TableName) const
	{
		char buff[STRBUFFER];
		sprintf(buff,formatStr,const_cast<char *>((Wstr2Str(TableName)).c_str()),PIOErrMess(MyErr));
		throw Zeus::libException(MyErr,Achar2Wstr(buff),*this);
	}
//
public:
	CatWriterOutputHFIDMC(const std::wstring& GrpName,const std::wstring& fname)
		:GenCollWriter<OutputFormatType>(),GrpName_(GrpName)
		,FName_(fname),TAB2Dgrp_(0)
	{}
//
	virtual ~CatWriterOutputHFIDMC(void)
	{Dispose();}
private:
//
	inline void			Dispose(void)
	{
		if(TAB2Dgrp_)
		{PIOCloseTAB2DGrp(&TAB2Dgrp_);TAB2Dgrp_	= 0;}
	}
//
	void				CreatStuff(void);
//
	std::wstring		FName_;
	std::wstring		GrpName_;
	PIOGroup			*TAB2Dgrp_;
};
//
// -------------------------------------------------------
class	CatWriterOutputHFIDMC_Orange10: public CatWriterOutputHFIDMC
{

	std::wstring	ScalesGrpName_;
	PIOGroup		*ScalesTAB2Dgrp_;

	inline double			Get_y0(const OutputFormatType & AuxInfo,OutputFormatType::StorageType::const_iterator elem)
	{
		return  (((elem->Cat_.FluxComptGLRT_ > -1.0e+10) && (elem->Cat_.FluxComptGLRT_ < 1.0e+10) && (elem->Cat_.SZ_ConversionCte_ > -1.0e+10) && (elem->Cat_.SZ_ConversionCte_ < 1.0e+10))?elem->Cat_.SZ_ConversionCte_ * elem->Cat_.FluxComptGLRT_:SZCATALOGUEDEFAULTVALUE);
	}
//
	inline double			Get_y0Err(const OutputFormatType & AuxInfo,OutputFormatType::StorageType::const_iterator elem)
	{
		return  (((elem->Cat_.ErrorBars_.FluxErrorBar_ > -1.0e+10) && (elem->Cat_.ErrorBars_.FluxErrorBar_ < 1.0e+10) && (elem->Cat_.SZ_ConversionCte_ > -1.0e+10) && (elem->Cat_.SZ_ConversionCte_ < 1.0e+10))?elem->Cat_.SZ_ConversionCte_ * elem->Cat_.ErrorBars_.FluxErrorBar_:SZCATALOGUEDEFAULTVALUE);
	}
//
	void	CreateScalesGroup(void);
//
	void	CreateScalesTAB2D(OutputFormatType::StorageType::const_iterator ptr);
//
	inline void			DisposeScalesGroup(void)
	{
		if(ScalesTAB2Dgrp_)
		{PIOCloseTAB2DGrp(&ScalesTAB2Dgrp_);ScalesTAB2Dgrp_	= 0;}
	}

protected:
//
	virtual void			Move2Array2D(const OutputFormatType & DataIN,LArr2D<double> & DataOUT);
//
	inline virtual int		GetNColumns(void) const
	{return 35;}

	inline virtual	void	CreateHeader(const OutputFormatType& data)
	{
		CreateHeaderHelper(data);
		if(!((data.Storage_[0]).Cat_.ScaleLikeNoise_.empty()))
		{
#ifdef SCALES1DINDB
			CreateScalesGroup();
#endif
		}
	}
//
public:
	CatWriterOutputHFIDMC_Orange10(const std::wstring& GrpName,const std::wstring& fname)
		:CatWriterOutputHFIDMC(GrpName,fname),ScalesGrpName_(GrpName + L"_SZscales"),ScalesTAB2Dgrp_(0)
	{}
//
	virtual ~CatWriterOutputHFIDMC_Orange10(void)
	{DisposeScalesGroup();}
};
// -------------------------------------------------------
class	CatWriterOutputHFIDMC_QA_Contours: public CatWriterOutputHFIDMC
{

protected:
//
	virtual void			Move2Array2D(const OutputFormatType & DataIN,LArr2D<double> & DataOUT);
//
	inline virtual int		GetNColumns(void) const
	{return QARESULTSSZ;}

	inline virtual	void	CreateHeader(const OutputFormatType& data)
	{
		CreateHeaderHelper(data);
	}
//
public:
	CatWriterOutputHFIDMC_QA_Contours(const std::wstring& GrpName,const std::wstring& fname)
		:CatWriterOutputHFIDMC(GrpName,fname)
	{}
//
	virtual ~CatWriterOutputHFIDMC_QA_Contours(void)
	{}
};
//
class	CatWriterOutputHFIDMC_Orange18: public CatWriterOutputHFIDMC
{
protected:
	virtual void			Move2Array2D(const OutputFormatType & DataIN,LArr2D<double> & DataOUT);
//
	inline virtual int		GetNColumns(void) const
	{return 14;}
//
	inline virtual	void	CreateHeader(const OutputFormatType& data)
	{
		CreateHeaderHelper(data);
	}
//
public:
	CatWriterOutputHFIDMC_Orange18(const std::wstring& GrpName,const std::wstring& fname)
		:CatWriterOutputHFIDMC(GrpName,fname)
	{}
//
	virtual ~CatWriterOutputHFIDMC_Orange18(void)
	{}
};
//
class	ReadObjResultsHFIDMC: public Zeus::GenCollReader<PeakCollReadbleType>
{
private:
	PeakCollReadbleType		data_;
	std::wstring			FName_;
	std::wstring			GrpName_;
	PIOGroup				*TAB2Dgrp_;
	int						NColumns_;
	int						NRows_;

	int						ReadBody(void);
	int						ReadHeader(void);

//
	inline void			HFIDMC_errIO(int errCode,wchar_t* msg,const std::wstring& xtraName) const
	{
		std::wstring errstring(msg);

		errstring += std::wstring(L" -> ");
		errstring += xtraName;
		throw Zeus::libException(errCode,errstring,*this);
	}
//
	inline void			HFIDMC_errIO_TAB2D(PIOLONG MyErr,char* formatStr,char* fieldName) const
	{
		char buff[STRBUFFER];
		sprintf(buff,formatStr,fieldName,PIOErrMess((PIOErr) MyErr));
		throw Zeus::libException((int) MyErr,Achar2Wstr(buff),*this);
	}
//
protected:
//
	virtual std::wstring	do_GetCollID(void) const
	{return  GrpName_  + std::wstring(L"/") + FName_;}
//
	virtual bool	do_Read(void)
	{
		if(!NRows_ || !NColumns_)
			return true;

		ReadBody();
		return true;		
	}
//
	virtual	PeakCollReadbleType&	do_GetData(void)
	{return data_;}
//
	virtual void	do_DisposeData(void)
	{
		data_.Storage_.clear();
		if(TAB2Dgrp_)
		{PIOCloseTAB2DGrp(&TAB2Dgrp_);TAB2Dgrp_ = 0;}
	}
//
	virtual PeakCollReadbleType::HeaderType *do_Initialize(void);
//
public:
	ReadObjResultsHFIDMC(const std::wstring& GrpName,const std::wstring& fname)
		: GenCollReader<PeakCollReadbleType>(),FName_(fname),GrpName_(GrpName + std::wstring(HFIDMCINTERMCATGROUPEXTENSION)),
			NColumns_(0),NRows_(0),TAB2Dgrp_(0)
	{}
//
	virtual ~ReadObjResultsHFIDMC(void)
	{
		do_DisposeData();
	}
};
//
class	ReadNonBlindCatalogueHFIDMC: public GenCollReader<CatalogueFormatType>
{
private:
	CatalogueFormatType							data_;
	std::wstring								FName_;
	std::wstring								DName_;
	PIOGroup									*TAB2Dgrp_;

	int					Convert2CatFormat(PIOLONG NData,PIODOUBLE* tData);
//
	int					ReadHeader(void);
//
	int					ReadBody(void);
//
	inline void			HFIDMC_errIO(int errCode,wchar_t* msg,const std::wstring& xtraName) const
	{
		std::wstring errstring(msg);

		errstring += std::wstring(L" -> ");
		errstring += xtraName;
		throw Zeus::libException(errCode,errstring,*this);
	}
//
	inline void			HFIDMC_errIO_TAB2D(PIOErr MyErr,char* formatStr,const std::wstring& TableName) const
	{
		char buff[STRBUFFER];
		sprintf(buff,formatStr,const_cast<char *>((Wstr2Str(TableName)).c_str()),PIOErrMess(MyErr));
		throw Zeus::libException(MyErr,Achar2Wstr(buff),*this);
	}
//
protected:

	virtual  std::wstring	do_GetCollID(void) const
	{return DName_  + std::wstring(L"/") + FName_;}
//
	virtual bool	do_Read(void)
	{
		ReadBody();
		return true;		
	}
//
	virtual	CatalogueFormatType&	do_GetData(void)
	{
		data_.Header_.NRows_ = data_.Storage_.size();
		return data_;
	}
//
	virtual void	do_DisposeData(void)
	{
		data_.Storage_.clear();
		if(TAB2Dgrp_)
		{PIOCloseTAB2DGrp(&TAB2Dgrp_);}
		TAB2Dgrp_	= 0;
	}
//
	virtual CatalogueFormatType::HeaderType *do_Initialize(void)
	{
		TAB2Dgrp_ = PIOOpenTAB2DGrp(const_cast<char *>(Wstr2Str(DName_).c_str()),"r");
		if (!TAB2Dgrp_)
		{
			HFIDMC_errIO(ERROR_COD_HFIDMCERROPENFL,ERROR_MSG_HFIDMCERROPENFL,DName_);
		}
		ReadHeader();
		return &(data_.Header_);
	}
//
public:
	ReadNonBlindCatalogueHFIDMC(const std::wstring& Dname,const std::wstring& fname)
		: GenCollReader<CatalogueFormatType>(),DName_(Dname),FName_(fname),TAB2Dgrp_(0)
	{}
//
	virtual ~ReadNonBlindCatalogueHFIDMC(void)
	{do_DisposeData();}
};
//
class	ReadQA_ProfilesLstHFIDMC: public GenCollReader<QA_ProfilesLstType>
{
private:
	QA_ProfilesLstType							data_;
	std::wstring								FName_;
	std::wstring								DName_;
	PIOGroup									*TAB2Dgrp_;

	int					Convert2CatFormat(PIOLONG NData,PIODOUBLE* tData);
//
	int					ReadBody(void);
//
	inline void			HFIDMC_errIO(int errCode,wchar_t* msg,const std::wstring& xtraName) const
	{
		std::wstring errstring(msg);

		errstring += std::wstring(L" -> ");
		errstring += xtraName;
		throw Zeus::libException(errCode,errstring,*this);
	}
//
	inline void			HFIDMC_errIO_TAB2D(PIOErr MyErr,char* formatStr,const std::wstring& TableName) const
	{
		char buff[STRBUFFER];
		sprintf(buff,formatStr,const_cast<char *>((Wstr2Str(TableName)).c_str()),PIOErrMess(MyErr));
		throw Zeus::libException(MyErr,Achar2Wstr(buff),*this);
	}
//
protected:

	virtual  std::wstring	do_GetCollID(void) const
	{return DName_  + std::wstring(L"/") + FName_;}
//
	virtual bool	do_Read(void)
	{
		ReadBody();
		return true;		
	}
//
	virtual	QA_ProfilesLstType&	do_GetData(void)
	{return data_;}
//
	virtual void	do_DisposeData(void)
	{
		data_.Storage_.clear();
		if(TAB2Dgrp_)
		{PIOCloseTAB2DGrp(&TAB2Dgrp_);}
		TAB2Dgrp_	= 0;
	}
//
	virtual QA_ProfilesLstType::HeaderType *do_Initialize(void)
	{
		TAB2Dgrp_ = PIOOpenTAB2DGrp(const_cast<char *>(Wstr2Str(DName_).c_str()),"r");
		if (!TAB2Dgrp_)
		{HFIDMC_errIO(ERROR_COD_HFIDMCERROPENFL,ERROR_MSG_HFIDMCERROPENFL,DName_);}
		return &(data_.Header_);
	}

public:
	ReadQA_ProfilesLstHFIDMC(const std::wstring& Dname,const std::wstring& fname)
		: GenCollReader<QA_ProfilesLstType>(),DName_(Dname),FName_(fname),TAB2Dgrp_(0)
	{}
//
	virtual ~ReadQA_ProfilesLstHFIDMC(void)
	{do_DisposeData();}
};
//
class	ReadQA_ColLstHFIDMC: public GenCollReader<QA_CltLstCatalogueType>
{
private:
	QA_CltLstCatalogueType						data_;
	std::wstring								FName_;
	std::wstring								DName_;
	PIOGroup									*TAB2Dgrp_;

	int					Convert2CatFormat(PIOLONG NData,PIODOUBLE* tData);
//
	int					ReadBody(void);
//
	inline void			HFIDMC_errIO(int errCode,wchar_t* msg,const std::wstring& xtraName) const
	{
		std::wstring errstring(msg);

		errstring += std::wstring(L" -> ");
		errstring += xtraName;
		throw Zeus::libException(errCode,errstring,*this);
	}
//
	inline void			HFIDMC_errIO_TAB2D(PIOErr MyErr,char* formatStr,const std::wstring& TableName) const
	{
		char buff[STRBUFFER];
		sprintf(buff,formatStr,const_cast<char *>((Wstr2Str(TableName)).c_str()),PIOErrMess(MyErr));
		throw Zeus::libException(MyErr,Achar2Wstr(buff),*this);
	}
//
protected:

	virtual  std::wstring	do_GetCollID(void) const
	{return DName_  + std::wstring(L"/") + FName_;}
//
	virtual bool	do_Read(void)
	{
		ReadBody();
		return true;		
	}
//
	virtual	QA_CltLstCatalogueType&	do_GetData(void)
	{return data_;}
//
	virtual void	do_DisposeData(void)
	{
		data_.Storage_.clear();
		if(TAB2Dgrp_)
		{PIOCloseTAB2DGrp(&TAB2Dgrp_);}
		TAB2Dgrp_	= 0;
	}
//
	virtual QA_CltLstCatalogueType::HeaderType *do_Initialize(void)
	{
		TAB2Dgrp_ = PIOOpenTAB2DGrp(const_cast<char *>(Wstr2Str(DName_).c_str()),"r");
		if (!TAB2Dgrp_)
		{
			HFIDMC_errIO(ERROR_COD_HFIDMCERROPENFL,ERROR_MSG_HFIDMCERROPENFL,DName_);
		}
		return &(data_.Header_);
	}

public:
	ReadQA_ColLstHFIDMC(const std::wstring& Dname,const std::wstring& fname)
		: GenCollReader<QA_CltLstCatalogueType>(),DName_(Dname),FName_(fname),TAB2Dgrp_(0)
	{}
//
	virtual ~ReadQA_ColLstHFIDMC(void)
	{do_DisposeData();}
};
//
template<typename T>
class	ReadVectHFIDMC:  public GenCollReader<LArr1D<T> >
{
private:
	PIOLONG						Sz_;
	std::wstring				FName_;
	std::wstring				GrpName_;
	PIOGroup					*VECTgrp_;
	LArr1D<T>					*data_;

//
	inline void			HFIDMC_errIO(int errCode,wchar_t* msg,const std::wstring& xtraName) const
	{
		std::wstring errstring(msg);

		errstring += std::wstring(L" -> ");
		errstring += xtraName;
		throw Zeus::libException(errCode,errstring,*this);
	}
//
	inline void			HFIDMC_errIO_Vect(PIOErr MyErr,char* formatStr,char* fieldName) const
	{
		char buff[STRBUFFER];
		sprintf(buff,formatStr,fieldName,PIOErrMess(MyErr));
		throw Zeus::libException(MyErr,Achar2Wstr(buff),*this);
	}
//
protected:

	virtual std::wstring	do_GetCollID(void) const
	{return GrpName_  + std::wstring(L"/") + FName_;}

	virtual bool	do_Read(void)
	{
		PIOLONG			MySz;
		T				*ptrData;
		T				*ptrDataAux;
		PIOSTRING		command;
		std::string		tObjName(Wstr2Str(FName_));

		sprintf(command,"begin=%" PRId64 ";end=%" PRId64,(PIOLONG) 0,(PIOLONG)(Sz_ - 1));

		if((MySz=PIOReadVECTObject((void **)&(ptrData),const_cast<char *>(tObjName.c_str()),Private::PIOTypeChooserType<T>::PIO_TYPE,command,VECTgrp_)) < 0)
		{HFIDMC_errIO_Vect((PIOErr) MySz,ERROR_MSG_HFIDMCREADVECTFORMAT,const_cast<char *>(tObjName.c_str()));}
		
		if(MySz != data_->getSz())
		{
			if((MySz = PIODeleteVECTTable(ptrData,VECTgrp_)) != 0)
				HFIDMC_ReportMemoryLeak((PIOErr) MySz,ERROR_MSG_HFIDMCDELETEVECTFORMAT,const_cast<char *>(tObjName.c_str()));
			HFIDMC_errIO(ERROR_COD_HFIDMCERRIMGWRONGSZ,ERROR_MSG_HFIDMCERRIMGWRONGSZ,FName_);		
		}

		ptrDataAux = ptrData;
		T		*dataDestPtr(data_->begin());

		for(int i=0;i<Sz_;++dataDestPtr,++ptrDataAux,++i)
		{
			*dataDestPtr = *ptrDataAux;
		}

		if((MySz = PIODeleteVECTTable(ptrData,VECTgrp_)) != 0)
			HFIDMC_ReportMemoryLeak((PIOErr) MySz,ERROR_MSG_HFIDMCDELETEVECTFORMAT,const_cast<char *>(tObjName.c_str()));

		return true;
	}

	virtual	LArr1D<T>&	do_GetData(void)
	{return *data_;}

	virtual void	do_DisposeData(void)
	{
		delete data_;data_ = 0;

		if(VECTgrp_)
		{PIOCloseVECTGrp(&VECTgrp_);}
		VECTgrp_	= 0;
	}

	virtual typename LArr1D<T>::HeaderType *do_Initialize(void)
	{
		VECTgrp_ = PIOOpenVECTGrp(const_cast<char *>(Wstr2Str(GrpName_).c_str()),"r");
		if (!VECTgrp_)
		{
			HFIDMC_errIO(ERROR_COD_HFIDMCERROPENFL,ERROR_MSG_HFIDMCERROPENFL,FName_);
		}
		data_ = new LArr1D<T>(Sz_);
		data_->reset();
		return reinterpret_cast<typename LArr1D<T>::HeaderType *>(data_->begin());
		//return reinterpret_cast<typename LArr1D<T>::HeaderType *>(0);
	} 
public:
	ReadVectHFIDMC(const std::wstring& Dname,const std::wstring& fname,int Sz)
		: GenCollReader<LArr1D<T> >(),GrpName_(Dname + std::wstring(HFIDMCPROFCOEFGROUPEXTENSION)),FName_(fname),Sz_(Sz),
		data_(0),VECTgrp_(0)
	{}
	virtual ~ReadVectHFIDMC(void)
	{
		do_DisposeData();
	}
};
//
template<typename T>
class	WriteVectHFIDMC: public GenCollWriter<LArr1D<T> >
{
	typedef LArr1D<T>	AtomT;
protected:
//
	inline virtual std::wstring	do_GetCollID(void) const
	{return GrpName_  + std::wstring(L"/") + FName_;}
//
	inline virtual int			do_Flush(void)
	{return 0;}
//
	inline virtual void			do_DisposeData(void)
	{Dispose();return ;}
//
	inline virtual int		do_Remove(void)
	{return HFIDMC_RemoveObject(GrpName_,FName_);}

	virtual bool	do_Initialize(void)
	{CreatStuff();return true;}

	virtual int	do_Write(const AtomT & data)
	{
		const PIOLONG	NItems(data.getSz());
		PIOSTRING		command;
		PIOLONG			MySZ;

		if(!NItems)
		{return 0;}

		sprintf(command,"begin=%" PRId64 ";end=%" PRId64,(PIOLONG) 0,(PIOLONG)(NItems - 1));

		if(((MySZ=PIOWriteVECTObject((void *) data.begin(),const_cast<char *>(Wstr2Str(FName_).c_str()),"PIODOUBLE",command,VECTgrp_))<0) ||
			(NItems != MySZ)
			)
		{
			HFIDMC_errIO_Vect((PIOErr) MySZ,ERROR_MSG_HFIDMCWRITEVECTFORMAT,std::wstring(L"GNFWProfile coeficients"));
		}

		return NItems;
	}
public:
	WriteVectHFIDMC(const std::wstring& GrpName,const std::wstring& fname)
		: GenCollWriter<LArr1D<T> >(),FName_(fname),GrpName_(GrpName + std::wstring(HFIDMCPROFCOEFGROUPEXTENSION)),VECTgrp_(0)
	{}
	virtual ~WriteVectHFIDMC(void)
	{Dispose();}
private:
//
	inline void			HFIDMC_errIO(int errCode,wchar_t* msg,const std::wstring& xtraName) const
	{
		std::wstring errstring(msg);

		errstring += std::wstring(L" -> ");
		errstring += xtraName;
		throw Zeus::libException(errCode,errstring,*this);
	}
//
	inline void			HFIDMC_errIO_Vect(PIOErr MyErr,char* formatStr,const std::wstring& TableName) const
	{
		char buff[STRBUFFER];
		sprintf(buff,formatStr,const_cast<char *>((Wstr2Str(TableName)).c_str()),PIOErrMess(MyErr));
		throw Zeus::libException(MyErr,Achar2Wstr(buff),*this);
	}
//
	void CreatStuff(void)
	{
		PIOErr			MyErr;
		std::string		tGrpNameStr(Wstr2Str(GrpName_));
		int				retries(MAXRETRIES);

_retry:
		if(PIOCheckGroup(const_cast<char *>(tGrpNameStr.c_str())))
		{
			if(MyErr = PIOCreateVECTGrp(const_cast<char *>(tGrpNameStr.c_str()))) 
			{
				if(--retries >= 0)
				{
					MySleep(1);
					goto _retry;
				}
				HFIDMC_errIO(ERROR_COD_HFIDMCERRCREATFL,ERROR_MSG_HFIDMCERRCREATFL,GrpName_);
			}
		}

		if (!(VECTgrp_ = PIOOpenVECTGrp(const_cast<char *>(tGrpNameStr.c_str()),"w")))
		{HFIDMC_errIO(ERROR_COD_HFIDMCERROPENFL,ERROR_MSG_HFIDMCERROPENFL,GrpName_);}
		if(MyErr = PIOCreateVECTObject(const_cast<char *>((Wstr2Str(FName_)).c_str()),"PIODOUBLE",VECTgrp_)) 
		{HFIDMC_errIO_Vect(MyErr,ERROR_MSG_HFIDMCCREATEVECTFORMAT,FName_);}
	}

	inline void	Dispose(void)
	{
		if(VECTgrp_)
		{PIOCloseVECTGrp(&VECTgrp_);}
		VECTgrp_	= 0;
	}

	std::wstring				FName_;
	std::wstring				GrpName_;
	PIOGroup					*VECTgrp_;
};

//
#endif //HFIDMC
//
#ifdef	LFIDPC
//
template<typename T>
class	WriteHealpixMapLFIDPC: public GenHealpixWriter<T>
{
protected:
	virtual  std::wstring	do_GetCollID(void) const
	{return FName_;}
//
	virtual bool	do_Initialize(void)
	{
		Dispose();
		CollHandle_ = HandleManager.makeHandle();
		CollHandle_->create(Wstr2Str(FName_),std::string(HEALPIXCOLL_LFIDPC_COLLTYPE));
		return true;
	}
//
	virtual int		do_Write(const Healpix_Map<T>& data,HealpixHeaderType * hd)
	{
		if(hd)
		{CreateHeader(*hd);}

		write_Healpix_map_to_dmc(Wstr2Str(FName_),data);
		return true;
	}
//
public:
	WriteHealpixMapLFIDPC(const std::wstring& fname,Healpix_Ordering_Scheme Scheme,coordsys CoordSys)
		:GenHealpixWriter<T>(),FName_(fname),Scheme_(Scheme),CoordSys_(CoordSys),CollHandle_(0)
	{}
	virtual ~WriteHealpixMapLFIDPC(void)
	{Dispose();}
private:
	std::wstring				FName_;
	iohandle					*CollHandle_;
	Healpix_Ordering_Scheme		Scheme_;
	coordsys					CoordSys_;
//
	inline void			Dispose(void)
	{
		if(CollHandle_)
		{delete CollHandle_;CollHandle_ = 0;}
	}
//
	void	CreateHeader(const HealpixHeaderType& hd)
	{
		CollHandle_->setKey<std::string>(std::string(HFIDMC_HEALPIXMAP_SCHEME),GetStrFromOrdering(Scheme_));
		CollHandle_->setKey<std::string>(std::string(HFIDMC_HEALPIXMAP_COORSYS),GetStrFromCoordsys(CoordSys_));
		
		if(!(hd.Columns_.empty()))
		{
			HealpixHeaderType::HealpixHeaderItemCollType::const_iterator	piv(hd.Columns_.begin());
			HealpixHeaderType::HealpixHeaderItemCollType::const_iterator	const end(hd.Columns_.end());
			for(;piv != end;++piv)
			{
				switch(piv->ValueType_)
				{
				case INT_TYPE:
					CollHandle_->setKey<int>(piv->Keyword_,piv->PODValue_.intType_);
					break;
				case FLOAT_TYPE:
					CollHandle_->setKey<float>(piv->Keyword_,piv->PODValue_.floatType_);
					break;
				case DOUBLE_TYPE:
					CollHandle_->setKey<double>(piv->Keyword_,piv->PODValue_.doubleType_);
					break;
				case BOOL_TYPE:
					CollHandle_->setKey<bool>(piv->Keyword_,piv->PODValue_.boolType_);
					break;
				case STRING_TYPE:
					CollHandle_->setKey<std::string>(piv->Keyword_,piv->stringValue_);
					break;
				default:
					break;
				}
			}
		}
	}
};
//
class		CatalogueWriterLFIDPC: public GenCollWriter<CatLineCollType>
{
protected:
	virtual  std::wstring	do_GetCollID(void) const
	{return GrpName_ + FName_;}
//
	virtual bool	do_Initialize(void)
	{
		Dispose();
		CollHandle_ = HandleManager.makeHandle();
		CollHandle_->create(Wstr2Str(GrpName_),std::string(INTFILE_LFIDPC_COLLTYPE));
		return true;
	}
//
	virtual int		GetNextColumn(const CatLineCollType & DataIN,
		const SZPS_ProfParamType& ProfParam,arr<double>& data,std::string& ColumnName) = 0;
//
	virtual int		do_Write(const CatLineCollType & dataIN)
	{
	
		if(dataIN.empty())
			return true;

		const int		CatSz(dataIN.size());
		arr<double>		VectDOUBLE(CatSz);
		std::string		ColumnName;

		while(CatSz == GetNextColumn(dataIN,ProfParam_,VectDOUBLE,ColumnName))
		{
			CollHandle_->writeColumn(ColumnName,VectDOUBLE);
		}
		return true;
	}
//
public:
	CatalogueWriterLFIDPC(const std::wstring& GrpName,const std::wstring& fname,const SZPS_ProfParamType& ProfParam)
		:GenCollWriter<CatLineCollType>(),GrpName_(GrpName),FName_(fname),ProfParam_(ProfParam),CollHandle_(0)
	{}
//
	virtual ~CatalogueWriterLFIDPC(void)
	{Dispose();}
private:
//
	inline void			Dispose(void)
	{
		if(CollHandle_)
		{
			delete CollHandle_;
			CollHandle_ = 0;
		}
	}
//
	iohandle					*CollHandle_;
	std::wstring				FName_;
	std::wstring				GrpName_;
	const SZPS_ProfParamType&	ProfParam_;
};
//
class		CatalogueWriterLFIDPCSZ: public CatalogueWriterLFIDPC
{
	int		CurrColumn_;
protected:
	virtual int		GetNextColumn(const CatLineCollType & DataIN,
		const SZPS_ProfParamType& ProfParam,arr<double>& data,std::string& ColumnName);
public:
	CatalogueWriterLFIDPCSZ(const std::wstring& GrpName,const std::wstring& fname,const SZPS_ProfParamType& ProfParam)
		:CatalogueWriterLFIDPC(GrpName,fname,ProfParam),CurrColumn_(0)
	{}
//
	virtual ~CatalogueWriterLFIDPCSZ(void)
	{}
};
//
class		CatalogueWriterLFIDPCPS: public CatalogueWriterLFIDPC
{
	int		CurrColumn_;
protected:
	virtual int		GetNextColumn(const CatLineCollType & DataIN,
		const SZPS_ProfParamType& ProfParam,arr<double>& data,std::string& ColumnName);
public:
	CatalogueWriterLFIDPCPS(const std::wstring& GrpName,const std::wstring& fname,const SZPS_ProfParamType& ProfParam)
		:CatalogueWriterLFIDPC(GrpName,fname,ProfParam),CurrColumn_(0)
	{}
//
	virtual ~CatalogueWriterLFIDPCPS(void)
	{}
};
//
class		WritePatchesGeomInfoLFIDPC: public GenCollWriter<PatchGeomType>
{
protected:
	virtual  std::wstring	do_GetCollID(void) const
	{return FName_;}
//
	virtual bool	do_Initialize(void)
	{
		Dispose();
		int	tOpen(0);
		CollHandle_ = HandleManager.makeHandle();
		try
		{
			CollHandle_->open(Wstr2Str(FName_),std::string(INTFILE_LFIDPC_COLLTYPE));	
			tOpen = 1;
		}
		catch(...)
		{
			CollHandle_->create(Wstr2Str(FName_),std::string(INTFILE_LFIDPC_COLLTYPE));		
		}
		if(!tOpen)
		{CollHandle_->open(Wstr2Str(FName_),std::string(INTFILE_LFIDPC_COLLTYPE));}
		return true;
	}
//
	virtual int		do_Write(const PatchGeomType & data)
	{WriteGeomInfoData(data); Dispose(); return true;}
//
public:
	WritePatchesGeomInfoLFIDPC(const std::wstring& fname)
		:Zeus::GenCollWriter<PatchGeomType>(),FName_(fname),CollHandle_(0)
	{}
//
	virtual ~WritePatchesGeomInfoLFIDPC(void)
	{Dispose();}
//
private:
//
	inline void			Dispose(void)
	{
		if(CollHandle_)
		{
			delete CollHandle_;
			CollHandle_ = 0;
		}
	}
//
	void	WriteGeomInfoData(const PatchGeomType & data);
//
	void	CreatHeader(const PatchGeomType & data);
//
	std::wstring		FName_;
	iohandle			*CollHandle_;
};
//
template<typename T>
class		ReadHealpixMapLFIDPC: public GenHealpixReader<T>
{
protected:
	virtual std::wstring	do_GetCollID(void) const
	{return MapName_;}
//
	virtual bool	do_Initialize(Zeus::HealpixHeaderType *hd)
	{
		do_DisposeData();
		CollHandle_ = HandleManager.makeHandle();
		CollHandle_->open(Wstr2Str(DirName_));
		if(!hd) return true;
		
		hd->NSide_ = Healpix_Base::npix2nside((int)CollHandle_->columnLength(CollHandle_->columnNumber(Wstr2Str(MapName_))));

		if(hd->Columns_.empty())
			return true;

		HealpixHeaderType::HealpixHeaderItemCollType::iterator	piv(hd->Columns_.begin());
		HealpixHeaderType::HealpixHeaderItemCollType::const_iterator	const end(hd->Columns_.end());

		for(;piv!=end;++piv)
		{
			
			if(!(CollHandle_->keyPresent(piv->Keyword_)))
			{piv->Reset();continue;}
			switch(piv->ValueType_)
			{
			case Zeus::INT_TYPE:
				CollHandle_->getKey(piv->Keyword_,piv->PODValue_.intType_);
				break;
			case Zeus::FLOAT_TYPE:
				CollHandle_->getKey(piv->Keyword_,piv->PODValue_.floatType_);
				break;
			case Zeus::DOUBLE_TYPE:
				CollHandle_->getKey(piv->Keyword_,piv->PODValue_.doubleType_);
				break;
			case Zeus::BOOL_TYPE:
				CollHandle_->getKey(piv->Keyword_,piv->PODValue_.boolType_);
				break;
			case Zeus::STRING_TYPE:
				CollHandle_->getKey(piv->Keyword_,piv->stringValue_);
				break;
			default:
				piv->Reset();
				break;
			}
		}
		return true;
	}
//
	virtual void	do_Release(void)
	{do_DisposeData();}
//
	virtual bool	do_Read(void)
	{
		read_Healpix_map_from_dmc(*CollHandle_, *HPixMap_, Wstr2Str(MapName_));
		return true;
	}
public:
	ReadHealpixMapLFIDPC(const std::wstring& DirName,const std::wstring& MapName)
		:GenHealpixReader<T>(),DirName_(DirName),MapName_(MapName),CollHandle_(0)
	{}
//
	virtual ~ReadHealpixMapLFIDPC(void)
	{do_DisposeData();}
private:
	std::wstring				MapName_;
	std::wstring				DirName_;
	iohandle					*CollHandle_;

	inline void	do_DisposeData(void)
	{
		delete CollHandle_;CollHandle_ = 0;
	}
//
};

//
template<typename T>
class		ReadWrkSpaceLFIDPC:  public GenCollReader<LArr2D<T> >
{
private:
	std::wstring	FName_;
	std::wstring	GrpName_;
	iohandle		*CollHandle_;
	LArr2D<T>		*data_;
	int				YSz_;
	int				Metric_;
//
	inline void			LFIDPC_errIO(int errCode,wchar_t* msg,std::wstring xtraName) const
	{
		std::wstring errstring(msg);
		errstring += std::wstring(L" -> ");
		errstring += xtraName;
		throw Zeus::libException(errCode,errstring,*this);
	}

protected:
	virtual std::wstring	do_GetCollID(void) const
	{return FName_;}
//
	virtual bool	do_Read(void)
	{
		arr<T>					tArrFLOAT;

		CollHandle_->readEntireColumn(Wstr2Str(GrpName_),tArrFLOAT);

		if(tArrFLOAT.size() != data_->getSz())
		{LFIDPC_errIO(ERROR_COD_HFIDMCERRIMGWRONGSZ,ERROR_MSG_HFIDMCERRIMGWRONGSZ,FName_);}

		T		*tDestPtr(data_->begin());
		T		*tArrFLOATptr(tArrFLOAT.begin());
		T		* const tArrFLOAT_End(tArrFLOAT.end());

		for(;tArrFLOATptr != tArrFLOAT_End;++tDestPtr,++tArrFLOATptr)
		{*tDestPtr = static_cast<T>(*tArrFLOATptr);}

		return true;		
	}
//
	virtual	LArr2D<T>&	do_GetData(void)
	{return *data_;}
//
	virtual void	do_DisposeData(void)
	{
		delete data_;data_ = 0;
		delete CollHandle_;CollHandle_ = 0;
	}
//
	virtual typename LArr2D<T>::HeaderType *do_Initialize(void)
	{
		do_DisposeData();
		CollHandle_ = HandleManager.makeHandle();
		CollHandle_->open(Wstr2Str(GrpName_));

		data_ = new LArr2D<T>(YSz_ * Metric_ ,Metric_);
		data_->reset();
		return reinterpret_cast<typename LArr2D<T>::HeaderType *>(data_);
	}
public:
	ReadWrkSpaceLFIDPC(const std::wstring& GrpName,const std::wstring& fname,int YSz,int XSz)
		: GenCollReader<LArr2D<T> >(),GrpName_(GrpName),FName_(fname),YSz_(YSz),Metric_(XSz),CollHandle_(0),data_(0)
	{}
//
	virtual ~ReadWrkSpaceLFIDPC(void)
	{do_DisposeData();}

};
//
class		LFIDPC_ParameterFileReader: public GenCollReader<ParamVarsStoreType>
{
private:
	std::wstring			FName_;
	ParamVarsStoreType		*data_;
	paramfile				*Params_;
//
	inline bool				Insert1element(const std::wstring& id,const DBField& value)
	{
		data_->Storage_.erase(id);
		return data_->Storage_.insert(ParamVarsStoreType::StorageType::value_type(id, value)).second;
	}
//
	void					ReadValuesFromLFI_Pipeline(void);

protected:

	virtual std::wstring	do_GetCollID(void) const
	{return FName_;}
//
	virtual bool	do_Read(void)
	{
		ReadValuesFromLFI_Pipeline();
		return true;		
	}
//
	virtual	ParamVarsStoreType&	do_GetData(void)
	{return *data_;}
//
	virtual void	do_DisposeData(void)
	{delete Params_;Params_ = 0;delete data_;data_ = 0;}
//
	virtual ParamVarsStoreType::HeaderType *do_Initialize(void)
	{
		do_DisposeData();
		iohandle_current:: Manager mng(Wstr2Str(FName_));
		Params_	= new paramfile(mng.getParams());
		return &(data_->Header_);
	} 

public:
	LFIDPC_ParameterFileReader(const std::wstring& fname)
		: GenCollReader<ParamVarsStoreType>(),FName_(fname),data_(0),Params_(0)
	{}
//
	virtual ~LFIDPC_ParameterFileReader(void)
	{
		do_DisposeData();
	}
//
};

//
class		ReadPatchesGeomInfoLFIDPC: public GenCollReader<PatchGeomType>
{
private:
	PatchGeomType		data_;
	std::wstring		FName_;
	iohandle			*inp_;

	void				ReadHeader(void);
	void				ReadBody(void);

protected:
	virtual std::wstring	do_GetCollID(void) const
	{return FName_;}
//
	virtual bool	do_Read(void)
	{
		ReadBody();
		return true;		
	}
//
	virtual	PatchGeomType&	do_GetData(void)
	{return data_;}
//
	virtual void	do_DisposeData(void)
	{
		data_.Storage_.clear();
		delete inp_;inp_ = 0;
	}
//
	virtual PatchGeomType::HeaderType *do_Initialize(void)
	{
		do_DisposeData();
		inp_ = HandleManager.makeHandle();
		inp_->open(Wstr2Str(FName_));
		ReadHeader();
		return &(data_.Header_);
	}
//
public:
	ReadPatchesGeomInfoLFIDPC(const std::wstring& fname)
		: GenCollReader<PatchGeomType>(),FName_(fname),inp_(0)
	{}
//
	virtual ~ReadPatchesGeomInfoLFIDPC(void)
	{do_DisposeData();}
//
};
//
class		ReadNonBlindCatalogueLFIDPC: public GenCollReader<NonBlingCatType>
{
private:
	NonBlingCatType		data_;
	std::wstring		FName_;
	iohandle			*CollHandle_;

	void				ReadHeader(void);
	void				ReadBody(void);

protected:
	virtual std::wstring	do_GetCollID(void) const
	{return FName_;}
//
	virtual bool	do_Read(void)
	{
		ReadBody();
		return true;		
	}
//
	virtual	NonBlingCatType&	do_GetData(void)
	{return data_;}
//
	virtual void	do_DisposeData(void)
	{
		data_.Storage_.clear();
		if(CollHandle_)
		{delete CollHandle_;CollHandle_ = 0;}
	}
//
	virtual NonBlingCatType::HeaderType *do_Initialize(void)
	{
		do_DisposeData();
		CollHandle_ = HandleManager.makeHandle();
		CollHandle_->open(Wstr2Str(FName_));
		ReadHeader();
		return &(data_.Header_);
	}
//
public:
	ReadNonBlindCatalogueLFIDPC(const std::wstring& fname)
		: GenCollReader<NonBlingCatType>(),FName_(fname),CollHandle_(0)
	{}
//
	virtual ~ReadNonBlindCatalogueLFIDPC(void)
	{do_DisposeData();}
//
};
//
class		WriteObjResultsLFIDPC: public Zeus::GenCollWriter<PeakCollType>
{
protected:
	virtual std::wstring	do_GetCollID(void) const
	{return FName_;}
//
	virtual bool	do_Initialize(void)
	{
		Dispose();
		int	tOpen(0);
		CollHandle_ = HandleManager.makeHandle();
		try
		{
			CollHandle_->open(Wstr2Str(FName_),std::string(INTFILE_LFIDPC_COLLTYPE));	
			tOpen = 1;
		}
		catch(...)
		{
			CollHandle_->create(Wstr2Str(FName_),std::string(INTFILE_LFIDPC_COLLTYPE));		
		}
		if(!tOpen)
		{CollHandle_->open(Wstr2Str(FName_),std::string(INTFILE_LFIDPC_COLLTYPE));}
	
		return true;
	}
//
	virtual int		do_Write(const PeakCollType & data)
	{
		int Res(WritePeaksColl(data));
		Dispose();
		return Res;
	}
public:
	WriteObjResultsLFIDPC(const std::wstring& fname,int InitObjID)
		:Zeus::GenCollWriter<PeakCollType>(),FName_(fname),
		ObjectID_(InitObjID),CollHandle_(0)
	{}
//
	virtual ~WriteObjResultsLFIDPC(void)
	{Dispose();}
private:
//
	inline void			Dispose(void)
	{
		if(CollHandle_)
		{
			delete CollHandle_;
			CollHandle_ = 0;
		}
	}
//
	int				WritePeaksColl(const PeakCollType& data);
//
	int					ObjectID_;
	std::wstring		FName_;
	iohandle			*CollHandle_;
};
//
class		ReadObjResultsLFIDPC: public Zeus::GenCollReader<PeakCollReadbleType>
{

private:
	PeakCollReadbleType		data_;
	std::wstring			FName_;
	iohandle				*CollHandle_;
//
	void					ReadBody(void);
protected:
	virtual std::wstring	do_GetCollID(void) const
	{return FName_;}
//
	virtual bool	do_Read(void)
	{
		ReadBody();
		return true;		
	}
//
	virtual	PeakCollReadbleType&	do_GetData(void)
	{return data_;}
//
	virtual void	do_DisposeData(void)
	{
		data_.Storage_.clear();
		delete CollHandle_;CollHandle_ = 0;
	}
//
	virtual PeakCollReadbleType::HeaderType *do_Initialize(void)
	{
		do_DisposeData();
		CollHandle_ = HandleManager.makeHandle();
		CollHandle_->open(Wstr2Str(FName_));
		return &(data_.Header_);
	}
public:
	ReadObjResultsLFIDPC(const std::wstring& fname)
		: GenCollReader<PeakCollReadbleType>(),FName_(fname),CollHandle_(0)
	{}
	virtual ~ReadObjResultsLFIDPC(void)
	{do_DisposeData();}
};
//
template<typename T>
class		WriteWrkSpaceLFIDPC: public GenCollWriter<LArr2D<T> >
{
protected:
	virtual std::wstring	do_GetCollID(void) const
	{return GrpName_ + FName_;}

	virtual bool	do_Initialize(void)
	{
		Dispose();
		int	tOpen(0);
		CollHandle_ = HandleManager.makeHandle();
		try
		{
			CollHandle_->open(Wstr2Str(GrpName_),std::string(IMGCOLL_LFIDPC_COLLTYPE));	
			tOpen = 1;
		}
		catch(...)
		{
			CollHandle_->create(Wstr2Str(GrpName_),std::string(IMGCOLL_LFIDPC_COLLTYPE));
		}
		if(!tOpen)
		{CollHandle_->open(Wstr2Str(GrpName_),std::string(INTFILE_LFIDPC_COLLTYPE));}
	
		return true;
	}

	virtual int		do_Write(const LArr2D<T>& data)
	{
		arr<T>						VectDest(data.getSz());
		T							*destPtr(VectDest.begin());
		const T						*VectOrgPtr(data.begin());
		const T						* const VectOrgEnd(data.end());

		while(VectOrgPtr != VectOrgEnd)
		{*destPtr++ = *VectOrgPtr++;}

		CollHandle_->appendColumn(Wstr2Str(FName_),VectDest);
		return true;
	}
public:
	WriteWrkSpaceLFIDPC(const std::wstring& GrpName,const std::wstring& fname)
		: GenCollWriter<LArr2D<T> >(),GrpName_(GrpName),FName_(fname),CollHandle_(0)
	{}
	virtual ~WriteWrkSpaceLFIDPC(void)
	{Dispose();}
private:
	inline void	Dispose(void)
	{delete CollHandle_;CollHandle_ = 0;}

	std::wstring	GrpName_;
	std::wstring	FName_;
	iohandle		*CollHandle_;
};

#endif	//LFIDPC

} // namespace Zeus

#endif //ZEUSINOUTPIPELINEH

