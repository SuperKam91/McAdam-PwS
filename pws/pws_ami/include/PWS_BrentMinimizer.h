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



//----------------------------------------
#ifndef PWSBRENTMINH
#define PWSBRENTMINH

#include "PWS_OddsEval.h"

//----------------------------------------


#define GOLD		((double)1.618034)
#define TINY		((double)1.0e-20)
#define GLIMIT		((double)100.0)
#define CGOLD		((double)0.3819660)
#define ZEPS		((double)1.0e-5)

#define SHFT(a,b,c,d)	{(a) = (b);(b) = (c);(c) = (d);}


class BrentLineMinimizer
{

public:
	BrentLineMinimizer(OddsEval* OddFunct,int NParams,int MaxIter, double BrentTol)
		:OddFunct_(OddFunct),NParams_(NParams),MaxIter_(MaxIter),BrentTol_(BrentTol)
	{
		CurrSearchDir_.resize(NParams);
		CurrArgs_.resize(NParams);
		ZeroArgs_.resize(NParams);
	}

	int							FindLineMinimum(const MinArrayType& searchDir,MinArrayType& Init_Final,MinArrayAtomType& FunctValue);

private:

	int							NParams_;
	double						BrentTol_;
	int							MaxIter_;
	OddsEval*					OddFunct_;
	MinArrayType				CurrSearchDir_;
	MinArrayType				CurrArgs_;
	MinArrayType				ZeroArgs_;

	inline void					SetZeroArgs(const MinArrayType& zeroArgs)
	{
		MinArrayType::const_iterator	pivOrg(zeroArgs.begin());
		MinArrayType::iterator			pivDest(ZeroArgs_.begin());

		for(int i=0;i<NParams_;++i,++pivOrg,++pivDest)
		{*pivDest = *pivOrg;}
	}

	inline void					SetSearchDirs(const MinArrayType& sDirs)
	{
		MinArrayType::const_iterator	pivOrg(sDirs.begin());
		MinArrayType::iterator			pivDest(CurrSearchDir_.begin());

		for(int i=0;i<NParams_;++i,++pivOrg,++pivDest)
		{*pivDest = *pivOrg;}
	}

	inline void					GetFinalResult(MinArrayAtomType finalPos,MinArrayType& finalResult) const
	{
		MinArrayType::const_iterator	pivZero(ZeroArgs_.begin());
		MinArrayType::const_iterator	pivDir(CurrSearchDir_.begin());
		MinArrayType::iterator			pivDest(finalResult.begin());

		for(int i=0;i<NParams_;++i,++pivZero,++pivDir,++pivDest)
		{*pivDest = *pivZero + (finalPos *  *pivDir);}
	}

	inline MinArrayAtomType		Funct_1_dimen(MinArrayAtomType ax)
	{
		MinArrayType::const_iterator	pivZero(ZeroArgs_.begin());
		MinArrayType::const_iterator	pivDir(CurrSearchDir_.begin());
		MinArrayType::iterator			pivArgs(CurrArgs_.begin());

		for(int i=0;i<NParams_;++i,++pivZero,++pivDir,++pivArgs)
		{*pivArgs = *pivZero + (ax * *pivDir);}

		return OddFunct_->OddsValue(NParams_,CurrArgs_);
	}

	inline void					Shft(MinArrayAtomType& a,MinArrayAtomType& b,MinArrayAtomType& c,MinArrayAtomType& d) const 
	{a = b;b = c;c = d;}

	void						Mnbrak(MinArrayAtomType& ax,MinArrayAtomType&  bx,MinArrayAtomType&  cx,MinArrayAtomType&  fa,MinArrayAtomType&  fb,MinArrayAtomType&  fc);
	int							Brent(MinArrayAtomType brAx,MinArrayAtomType brBx,MinArrayAtomType brCx,MinArrayAtomType tol,MinArrayAtomType& xResult,MinArrayAtomType& fxResult);
};



#endif //PWSBRENTMINH

