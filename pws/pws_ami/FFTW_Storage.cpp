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


#include "FFTW_Storage.h"

//---------------------------------------------------------------------------
namespace fftw{

long CalcStoSize_Atoms(int dimension,const int *dims,bool xtra,bool isReal)
{
	int xtraSp(dims[0]);
	long blocksTemp(1);

	if(dimension > 0){
		blocksTemp=dims[1];
		for(int i=2;i<dimension;i++)
		{blocksTemp *= dims[i];}
	}
	if(xtra && isReal){xtraSp = 2 * (((int)(xtraSp/2)) + 1);}
	return blocksTemp * xtraSp;
}

bool checkDims(int dimension,const int *dims)
{
	for(int i=0;i<dimension;++i)
	{if(dims[i]<=0) return false;}
	return true;
}



}
