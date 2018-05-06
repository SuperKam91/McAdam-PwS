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



//----------------------------------
#ifndef OBJFACTORYH
#define OBJFACTORYH

#include "ZEUS_Object.h"
#include "ZEUS_ObjGauss.h"
#include "ZEUS_ObjSZ.h"
#include "ZEUS_ObjSZNgai.h"
//----------------------------------

typedef Zeus::RealPlane<double>		RealPlaneSurfType;

struct	ObjInfo
{
	Zeus::ObjFilterRealParams		RealParams_;
	RealPlaneSurfType				surf_;
	RealPlaneSurfType::AtomType		CentralPixAmpl_;
};

class	ObjFactory
{
public:
//
	struct	ObjFactoryArgs
	{
		int		ObjType_;
		double	PixSz_;
		int		PatchSz_;
		int		SynID_;
		virtual	~ObjFactoryArgs(void)
		{}
	};
//
	struct	ObjFactoryArgsSZ : public ObjFactoryArgs
	{
		Zeus::SZPS_ProfParamType	ProfParams_;
		virtual	~ObjFactoryArgsSZ(void)
		{}
	};
//
	struct	ObjFactoryArgsSZNgai : public ObjFactoryArgsSZ
	{
		int				ContextID_;
		std::wstring	DirInMasks_;
		virtual			~ObjFactoryArgsSZNgai(void)
		{}
	};
//
	ObjFactory(ObjFactoryArgs *args);
//
	ObjInfo	MakeObjectShape(double Radius);
//
	~ObjFactory(void)
	{delete Obj_;}
//
private:
//
	inline  void		errParOutOfRange(void) const
	{
		throw Zeus::libException(ERRCOD_PWS_PARAMOUTOFRANGE,ERRMSG_PWS_PARAMOUTOFRANGE,"ObjFactory::MakeObject");
	}
//
	Zeus::SrcObject		*Obj_;
};


#endif //OBJFACTORYH

