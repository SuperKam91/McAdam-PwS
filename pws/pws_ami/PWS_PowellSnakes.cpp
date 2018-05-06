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

#include "PWS_PowellSnakes.h"

//------------------------------------------


int		PowellSnake::PowellMinimunSearch(MinArrayType& InitGuessFinal,MinArrayAtomType& fMinimum)
{

	MinArrayAtomType					fValueTemp,fValueTemp1,fValueTemp2,del;
	int									iter,jbig,j,i,dirNempty;
	int									stat;

	MinArrayType						ptBufferWS(NParams_);
	MinArrayType::iterator				const ptBuffer(ptBufferWS.begin());
	MinArrayType						ptBuffer1WS(NParams_);
	MinArrayType::iterator				const ptBuffer1(ptBuffer1WS.begin());
	MinArrayType						dirBufferWS(NParams_);
	MinArrayType::iterator				const dirBuffer(dirBufferWS.begin());
	std::vector<int>::iterator			const AllowedSubSpaces(AllowedSubSpaces_.begin());
	Zeus::LArr2D<double>::iterator		const SearchDirs(SearchDirs_.begin());
	int									NDir;

	ResetDirections();
	for(i=0;i<NParams_;++i)
	{ptBuffer[i]=InitGuessFinal[i];}

	fValueTemp = Odds_->OddsValue(NParams_,InitGuessFinal);

	for(iter=0;;++iter){
		NDir = 0;
		fValueTemp1 = fValueTemp;
		jbig = -1;del = ((MinArrayAtomType)-1.0);
		for(j=0;j<NParams_;++j)
		{
			for(i=0;i<NParams_;++i)
			{
				if(AllowedSubSpaces[i]){dirBuffer[i] = *(SearchDirs + (j * NParams_) + i);}
				else{dirBuffer[i]=((MinArrayAtomType)0.0);}
			}

			if(!CheckDir(dirBufferWS)) continue;

			++NDir;
			fValueTemp2 = fValueTemp;

			if(!(stat = Brent_->FindLineMinimum(dirBufferWS,InitGuessFinal,fValueTemp)))
				return false;

			if(std::abs(fValueTemp2 - fValueTemp) > del)
			{
				del = std::abs(fValueTemp2 - fValueTemp);
				jbig = j;
			}
		}
		if((NDir <=1) ||((((MinArrayAtomType)2.0)*std::abs(fValueTemp1-fValueTemp)) <= (PTol_ * (std::abs(fValueTemp1) + std::abs(fValueTemp)))))
		{
			fMinimum = fValueTemp;
			return true;
		}

		if(MaxIter_ == iter) return false;

		for(j=0;j<NParams_;++j)
		{
			if(AllowedSubSpaces[j])
			{
				ptBuffer1[j] = ((MinArrayAtomType)2.0) * InitGuessFinal[j] - ptBuffer[j];
				dirBuffer[j] = InitGuessFinal[j] - ptBuffer[j];
			}
			else
			{
				ptBuffer1[j] = InitGuessFinal[j];
				dirBuffer[j]=((MinArrayAtomType)0.0);
			}
			ptBuffer[j] = InitGuessFinal[j];
		}

		fValueTemp2 = Odds_->OddsValue(NParams_,ptBuffer1WS);;

		if(fValueTemp2 < fValueTemp1)
		{
			MinArrayAtomType t(((MinArrayAtomType)2.0) * (fValueTemp1 - ((MinArrayAtomType)2.0) * fValueTemp + fValueTemp2) * ((fValueTemp1 - fValueTemp - del) * (fValueTemp1 - fValueTemp - del)) -
				del * ((fValueTemp1 - fValueTemp2)*(fValueTemp1 - fValueTemp2)));
			if(t < ((MinArrayAtomType)0.0))
			{
				if(!(stat = Brent_->FindLineMinimum(dirBufferWS,InitGuessFinal,fValueTemp)))
					return false;

				if((dirNempty = FirstNonEmptyDir()) == -1)
					return false;

				for(i=0;i<NParams_;++i)
				{
					if(jbig != dirNempty) 
					{*(SearchDirs + (jbig * NParams_) + i) =  *(SearchDirs + (dirNempty * NParams_) + i);}
					*(SearchDirs + (dirNempty * NParams_) + i) = dirBuffer[i];
				}
			}
		}
	}

}
