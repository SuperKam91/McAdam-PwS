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

//----------------------------
#include "MC_ObjFactory.h"
//----------------------------

ObjFactory::ObjFactory(ObjFactory::ObjFactoryArgs *args)
{
	switch(args->ObjType_)
	{
	case 0:
		Obj_ = new Zeus::SrcObjGauss(L"GaussianObj",args->PixSz_ * RAD2ARCMIN,args->PatchSz_);
		break;
	case 1:
		{
		ObjFactoryArgsSZ*	ptr;
		if(!(ptr = dynamic_cast<ObjFactoryArgsSZ*>(args)))
		{errParOutOfRange();}
		Obj_ = new Zeus::SrcObjSZ(L"BetaSZ_Obj",ptr->PixSz_ * RAD2ARCMIN,ptr->PatchSz_,ptr->ProfParams_.VirialRatio_);
		}
		break;
	case 2:
		{
		ObjFactoryArgsSZNgai*	ptr;
		if(!(ptr = dynamic_cast<ObjFactoryArgsSZNgai*>(args)))
		{errParOutOfRange();}
		Obj_ = new Zeus::SrcObjSZNgai(L"NgaiSZ_Obj",ptr->PixSz_ * RAD2ARCMIN,ptr->PatchSz_,ptr->ProfParams_,ptr->ContextID_,ptr->DirInMasks_,ptr->SynID_);
		}
		break;
	default:
		errParOutOfRange();
	}

	Obj_->Initialise();
}

ObjInfo		ObjFactory::MakeObjectShape(double Radius)
{
	ObjInfo						temp;
	Zeus::ObjFilterRealParams	objParam(Radius);

	Obj_->ChangeObjParams(objParam,0.0,0.0,true);
	temp.surf_.Swap(Obj_->GetObjectBufferNonConstRef());
	temp.CentralPixAmpl_	= *(temp.surf_.GetInnerData().begin());
	temp.RealParams_		= objParam;
	return temp;
}
