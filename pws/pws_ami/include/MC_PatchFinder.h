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



//---------------------------------
#ifndef MCPATCHFINDER
#define MCPATCHFINDER

#include "ZEUS_InOut.h"
#include "MC_HPixMapPatchCutting.h"

//---------------------------------
//#define MC_FINDPATCHFORMSTR	L"Index-> %04d, Patch-> %02d, Reject-> %6.2f, SNR-> %6.2f\nCCTh-> %6.2f, CCPhi-> %6.2f\nC0Th-> %6.2f, C0Phi-> %6.2f, C1Th-> %6.2f, C1Phi-> %6.2f\nC2Th-> %6.2f, C2Phi-> %6.2f, C3Th-> %6.2f, C3Phi-> %6.2f\n\n"

#define MC_FINDPATCHFORMSTR	L"Index-> %04d, Patch-> %02d, Reject-> %6.2f, SNR-> %6.2f LON-> %6.2f, LAT-> %6.2f\n"
#define MC_PTGFORMSTR	L"CoLat(rad)-> %6.4f, Lng(rad)-> %6.4f, Lat(deg)-> %6.2f, Lng(deg)-> %6.2f\n\n"

inline bool SortPatchIDs(Zeus::PatchGeomType::StorageType::const_iterator lhs,Zeus::PatchGeomType::StorageType::const_iterator rhs)
{return lhs->SrcIndex_ < rhs->SrcIndex_;}

class	PatchFinder
{
	typedef Zeus::PatchGeomType::StorageType::const_iterator	PatchsInfoItType;		
	typedef	std::vector<PatchsInfoItType>						PatchIDsCollType;

public:
	PatchFinder(EnvIDsType InOutEnvID,const std::wstring& Direct)
		:InOutEnvID_(InOutEnvID),Direct_(Direct)
	{}
	void FindPacthesFromPtg(const pointing& ptg);
private:
	PatchFinder(const PatchFinder&);
	PatchFinder& operator=(const PatchFinder&);

	inline  void			errPatchFinder(int errCode,const wchar_t* msg) const
	{
		std::wstring errstring(msg);

		errstring += std::wstring(L" -> ");
		errstring += Direct_;
		throw Zeus::libException(errCode,errstring,*this);
	}

	void			InitializeFromFile(void);
	int				FindPatches(const pointing& ptg);
	static bool		CheckOnePatch(const pointing& ptg,const Zeus::PatchGeomType::StorageType::value_type& patch);
	void			PrintPatchProps(const pointing& ptg);
	void			PrintOnePatch(const Zeus::PatchGeomType::StorageType::value_type& ptLine) const;

	const EnvIDsType	InOutEnvID_;
	std::wstring		Direct_;
	Zeus::PatchGeomType		PatchGeom_;
	PatchIDsCollType	PatchsIDs_;
	int					NSide_;
};

#endif //MCPATCHFINDER

