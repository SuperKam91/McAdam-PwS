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
#ifndef POWELLSNAKESH
#define POWELLSNAKESH

#include "PWS_OddsEval.h"
#include "PWS_BrentMinimizer.h"

//------------------------------------------


class PowellSnake
{
public:


	PowellSnake(BrentLineMinimizer *BrentMin,OddsEval *OddsFunct,int NParams,int PowMaxIter,MinArrayAtomType PowTol)
		:Brent_(BrentMin),Odds_(OddsFunct),NParams_(NParams),MaxIter_(PowMaxIter),PTol_(PowTol),SearchDirs_(NParams*NParams,NParams,0.0)
	{
		AllowedSubSpaces_.resize(NParams_);
		InitialDirectionScales_.resize(NParams_);
		SetAllowedSubSpaces(std::vector<int>());
		SetDirectionScales(MinArrayType());
		ResetDirections();
	}

	int					PowellMinimunSearch(MinArrayType& InitGuessFinal,MinArrayAtomType& fMinimum);

	inline void			SetDirectionScales(const MinArrayType& InitDirScales)
	{

		MinArrayType::const_iterator	pivOrg(InitDirScales.begin());
		MinArrayType::iterator			pivDest(InitialDirectionScales_.begin());

		for(int i=0;i<NParams_;++i,++pivDest)
		{
			if(InitDirScales.empty())
			{*pivDest = ((MinArrayAtomType)1.0);}
			else
			{
				*pivDest = *pivOrg;
				++pivOrg;
			}
		}
	}

	inline void 		SetAllowedSubSpaces(const std::vector<int>& Allowed)
	{
		std::vector<int>::const_iterator	pivAllow(Allowed.begin());
		std::vector<int>::iterator			pivAllowSub(AllowedSubSpaces_.begin());
		std::vector<int>::const_iterator	const endAllowSub(AllowedSubSpaces_.end());

		for(;pivAllowSub != endAllowSub;++pivAllowSub)
		{
			if(Allowed.empty())
			{*pivAllowSub = 1;}
			else
			{*pivAllowSub = *pivAllow;++pivAllow;}
		}
	}


private:
	inline void			ResetDirections(void)
	{
		Zeus::LArr2D<double>::iterator	pivSDir(SearchDirs_.begin());
		Zeus::LArr2D<double>::iterator	const endSDir(pivSDir + SearchDirs_.getSz());

		for(;pivSDir != endSDir;++pivSDir)
		{*pivSDir = ((MinArrayAtomType)0.0);}

		pivSDir = SearchDirs_.begin();

		MinArrayType::const_iterator		pivInitDir(InitialDirectionScales_.begin());
		std::vector<int>::const_iterator	pivAllDir(AllowedSubSpaces_.begin());
		std::vector<int>::const_iterator	const endAllDir(AllowedSubSpaces_.end());

		for(;pivAllDir != endAllDir;++pivAllDir,++pivInitDir,pivSDir += (NParams_ + 1))
		{
			if(*pivAllDir) *pivSDir = *pivInitDir;
		}

	}

	inline bool 		CheckDir(Zeus::LArr2D<double>::const_iterator pivDir) const
	{
		Zeus::LArr2D<double>::const_iterator const	endDir(pivDir + NParams_);

		for(;pivDir != endDir;++pivDir)
		{
			if(*pivDir != ((MinArrayAtomType)0.0)) return true;
		}
		return false;
	}

	inline bool 		CheckDir(const MinArrayType& Dir) const
	{
		MinArrayType::const_iterator pivDir(Dir.begin());
		MinArrayType::const_iterator const	endDir(Dir.end());

		for(;pivDir != endDir;++pivDir)
		{
			if(*pivDir != ((MinArrayAtomType)0.0)) return true;
		}
		return false;
	}

	inline int 			FirstNonEmptyDir(void) const
	{
		Zeus::LArr2D<double>::const_iterator	pivDirs(SearchDirs_.begin() + (NParams_ * (NParams_ - 1)));

		for(int i=NParams_-1;i>=0;--i,pivDirs -= NParams_)
		{if(CheckDir(pivDirs)) return i;}
		return -1;
	}

//---------------------------------------------	
	int									NParams_;
	int									MaxIter_;
	MinArrayAtomType					PTol_;
	BrentLineMinimizer					*Brent_;
	OddsEval							*Odds_;
	Zeus::LArr2D<double>				SearchDirs_;
	std::vector<int>					AllowedSubSpaces_;
	MinArrayType						InitialDirectionScales_;
};


#endif //POWELLSNAKESH

