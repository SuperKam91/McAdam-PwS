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
#ifndef ZEUS_BINSH
#define ZEUS_BINSH

#include "ZEUS_General.h"
#include "ZEUS_Exceptions.h"

//---------------------------------------------------------------------------

namespace Zeus
{

template<typename T>
class ValuesBins
{
	typedef	std::vector<T>	StorageType;

	int						NBins_;
	StorageType				Values_;

	inline T	do_CorrectVal(T binT) const
	{
		const double NBindsD(static_cast<double>(NBins_));
		if(binT < 0) binT = -binT;
		if(binT < NBindsD) return binT;
		const double NBindsD2X(2.0 * NBindsD);
		binT = std::fmod(binT,NBindsD2X);
		if(binT < NBindsD) return binT;
		return NBindsD2X - binT;
	}

	inline void errValues_binning(void) const {throw Zeus::libException(ERROR_COD_BINNING,ERROR_MSG_BINNING,*this);}
public:
	inline int				getNBins(void) const
	{return NBins_;}

	inline int	CorrectBin(int Bin) const
	{
		if(Bin < 0) Bin = -Bin;
		if(Bin < NBins_) return Bin;
		Bin %= (NBins_ <<1);
		if(Bin < NBins_) return Bin;
		return (NBins_<<1) - Bin;
	}

	inline T				getValue(int bin) const
	{
		if(NBins_ < 1)	errValues_binning();
		return Values_[CorrectBin(bin)];
	}

	inline T				getValue(T bin) const
	{
		if(NBins_ < 1)	errValues_binning();		
		bin	= do_CorrectVal(bin);
		const int	i(Zeus::toInt(bin));
		if(i>=NBins_) return Values_[NBins_];
		return (Values_[i + 1] - Values_[i]) * (bin - std::floor(bin)) + Values_[i];
	}

	inline T		GetBinInterpol(T bin) const
	{
		typename StorageType::const_iterator	const pivOrg(Values_.begin());
		typename StorageType::const_iterator	piv(std::upper_bound(pivOrg,Values_.end(),bin));
		if(piv == Values_.end())	return  static_cast<double>(NBins_);
		if(piv == pivOrg)			return	0.0;
		
		const	double	v0(*(piv-1));
		const	double	d(*piv - v0);
		const	double	offset(static_cast<double>(static_cast<int>(piv - pivOrg)));
		return ((bin - v0)/d) + (offset - 1);
	}

	inline int		GetBin(T bin) const
	{
		typename StorageType::const_iterator	const pivOrg(Values_.begin());
		typename StorageType::const_iterator	piv(std::upper_bound(pivOrg,Values_.end(),bin));
		if(piv == Values_.end())	return  NBins_;
		if(piv == pivOrg)			return	0;
		return static_cast<int>(((((*piv) - bin) < (bin - (*(piv-1)))) ? piv : piv - 1) - pivOrg);
	}

	inline T		GetNearestScale(T scale) const
	{
		typename StorageType::const_iterator	piv(std::upper_bound(Values_.begin(),Values_.end(),scale));

		if(piv == Values_.begin()) return *piv;
		if(piv == Values_.end()) return *(piv-1);
		
		return ((*piv - scale) < (scale - *(piv-1)))? *piv : *(piv-1);
	}

	inline T		GetNearestScaleButFirst(T scale) const
	{
		typename StorageType::const_iterator	tbeg(Values_.begin()+1);
		typename StorageType::const_iterator	piv(std::upper_bound(tbeg,Values_.end(),scale));

		if(piv == tbeg) return *piv;
		if(piv == Values_.end()) return *(piv-1);
		
		return ((*piv - scale) < (scale - *(piv-1)))? *piv : *(piv-1);
	}


	template<typename FILLER>
	inline void Init(FILLER& filler)
	{
		Values_.clear();
		NBins_ = filler.GetNBins(); Values_.resize(NBins_ + 1);
		for(int i=0;i!=(NBins_ + 1);++i)
		{Values_[i] = filler(i);}
		std::sort(Values_.begin(),Values_.end());
	}

	ValuesBins(void)
		:NBins_(-1)
	{}
};

}

#endif //ZEUS_BINSH
