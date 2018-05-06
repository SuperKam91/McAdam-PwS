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

//-----------------------------

#include <memory>
#include "ZEUS_PhysicsMath.h"
#include "MC_PatchFinder.h"
#include "MC_HPixMapPatchCutting.h"

//-----------------------------

void	PatchFinder::InitializeFromFile(void)
{
	std::auto_ptr<Zeus::GenCollReader<Zeus::PatchGeomType> >
		FReader(Zeus::GetPatchGeomFileReaderHandler(Loki::Type2Type<Zeus::PatchGeomType>(),InOutEnvID_,Direct_, std::wstring(PATCH_GEOMPROPSFNAME)));
	FReader->Initialize();
	FReader->Read();
	FReader->Release(PatchGeom_);
	NSide_ = PatchGeom_.Header_.NSide_;
}

int		PatchFinder::FindPatches(const pointing& ptg)
{
	PatchsIDs_.clear();
	Zeus::PatchGeomType::StorageType::const_iterator			piv(PatchGeom_.Storage_.begin());
	Zeus::PatchGeomType::StorageType::const_iterator const		end(PatchGeom_.Storage_.end());

	for(;piv != end;++piv)
	{
		if((!(piv->PatchValid_)) && CheckOnePatch(ptg,*piv))
		{PatchsIDs_.push_back(piv);}
	}
	return static_cast<int>(PatchsIDs_.size());
}

bool	PatchFinder::CheckOnePatch(const pointing& ptg,const Zeus::PatchGeomType::StorageType::value_type& patch)
{

	if((dotprod(vec3(ptg), crossprod(vec3(patch.X0Y0_),vec3(patch.X0YL_))) > 0)		||
		(dotprod(vec3(ptg), crossprod(vec3(patch.X0YL_),vec3(patch.XLYL_))) > 0)	||
		(dotprod(vec3(ptg), crossprod(vec3(patch.XLYL_),vec3(patch.XLY0_))) > 0)	||
		(dotprod(vec3(ptg), crossprod(vec3(patch.XLY0_),vec3(patch.X0Y0_))) > 0)
		)
	return false;

	return true;
}

void	PatchFinder::PrintOnePatch(const Zeus::PatchGeomType::StorageType::value_type& ptLine) const
{
	wchar_t	buffer[BUFFERMAXCHAR];

	PRINTINTOBUFFERFUNCT (buffer,BUFFERMAXCHAR,MC_FINDPATCHFORMSTR,
		ptLine.PatchNumber_,
		ptLine.SrcIndex_,
		ptLine.FinalRejectRatio_,
		ptLine.SNR_,
		static_cast<double>(ptLine.X0Y0Ptg_.phi*RAD2DEGREE),
		static_cast<double>((PIOVER2 - ptLine.X0Y0Ptg_.theta)*RAD2DEGREE)
	);

	(Zeus::ConManager::Instance())->PrintStr2Console(buffer);

	/*
		static_cast<double>((PIOVER2 - ptLine.X0Y0_.theta)*RAD2DEGREE)
		static_cast<double>(ptLine.X0Y0_.phi*RAD2DEGREE),
		static_cast<double>((PIOVER2 - ptLine.X0YL_.theta)*RAD2DEGREE),
		static_cast<double>(ptLine.X0YL_.phi*RAD2DEGREE),
		static_cast<double>((PIOVER2 - ptLine.XLYL_.theta)*RAD2DEGREE),
		static_cast<double>(ptLine.XLYL_.phi*RAD2DEGREE),
		static_cast<double>((PIOVER2 - ptLine.XLY0_.theta)*RAD2DEGREE),
		static_cast<double>(ptLine.XLY0_.phi*RAD2DEGREE)

*/
}

void	PatchFinder::PrintPatchProps(const pointing& ptg)
{
	wchar_t	buffer[BUFFERMAXCHAR];

	PRINTINTOBUFFERFUNCT (buffer,BUFFERMAXCHAR,MC_PTGFORMSTR,
		static_cast<double>(ptg.theta),
		static_cast<double>(ptg.phi),
		static_cast<double>((PIOVER2 - ptg.theta)*RAD2DEGREE),
		static_cast<double>(ptg.phi*RAD2DEGREE)
	);

	(Zeus::ConManager::Instance())->PrintStr2Console(buffer);

	std::sort(PatchsIDs_.begin(),PatchsIDs_.end(),SortPatchIDs);

	PatchIDsCollType::const_iterator	piv(PatchsIDs_.begin());
	PatchIDsCollType::const_iterator	const end(PatchsIDs_.end());

	for(;piv != end;++piv)
	{PrintOnePatch(*(*piv));}
}


void 	PatchFinder::FindPacthesFromPtg(const pointing& ptg)
{
	InitializeFromFile();
	FindPatches(ptg);
	PrintPatchProps(ptg);
}
