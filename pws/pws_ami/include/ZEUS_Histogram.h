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
#ifndef ZEUS_HISTOGRAMH
#define ZEUS_HISTOGRAMH

#include <vector>
#include "ZEUS_GaussianRandomGen.h"

/*
This pretends to be a general class for making histograms
*/


namespace	Zeus
{

template<typename HIST_EL>
struct HistResultElem
{
	HIST_EL		Bin_;
	HIST_EL		BinMin_;
	HIST_EL		BinMax_;
	int			Count_;

	HistResultElem(const HIST_EL& Bin)
		:Bin_(Bin),Count_(0)
	{}
};

template<typename HIST_EL>
struct	HistEmptyCellRemoveFunctor
{
	int		MinCount_;
	bool	operator()(const HistResultElem<HIST_EL>& el)
	{
		return (el.Count_ < MinCount_);
	}
	HistEmptyCellRemoveFunctor(int MinCount)
		:MinCount_(MinCount)
	{}
};

template<typename HIST_EL>
struct HistResultElemSortFunctor
{
	bool	(*ptFnComp_)(const HIST_EL& first,const HIST_EL& sec);

	bool	operator()(const HistResultElem<HIST_EL>& fst,const HistResultElem<HIST_EL>& sec)
	{
		return ptFnComp_(fst.Bin_,sec.Bin_);
	}

	HistResultElemSortFunctor(bool	(*ptFnComp)(const HIST_EL& first,const HIST_EL& sec))
		:ptFnComp_(ptFnComp)
	{}
};


template<typename HIST_EL,typename DATA_EL>
class Histogram
{
public:
	typedef HistResultElem<HIST_EL>		HistResultElemType;
	typedef std::vector<HistResultElemType>	HistResultType;

private:
	typedef std::vector<DATA_EL>					DATATYPECOLL;
	typedef typename DATATYPECOLL::const_iterator	DATATYPECOLLIT;
	typedef std::vector<HIST_EL>					DBINSTYPECOLL;

	const DATATYPECOLL&				Data_;
	DBINSTYPECOLL&					Bins_;
	HIST_EL							FstBin_;
	HIST_EL							(*ptFnExtract_)(DATATYPECOLLIT ptr);
	bool							(*ptFnComp_)(const HIST_EL& first,const HIST_EL& sec);
	bool							(*ptFnCompInBin_)(const HIST_EL& first,const HIST_EL& sec);
	int								UnderFlow_;
	int								OverFlow_;

public:
	HIST_EL	DoHistogram(HistResultType& result)
	{
		SetUpResultBins(result);
		DATATYPECOLLIT	piv(Data_.begin());
		DATATYPECOLLIT	const end(Data_.end());
		HistResultElemSortFunctor<HIST_EL>	SortFunctor(ptFnComp_);

		for(;piv != end;++piv)
		{
			HIST_EL		t(ptFnExtract_(piv));
			typename HistResultType::iterator	HigherBinEdge(std::lower_bound(result.begin(),result.end(),HistResultElemType(t),SortFunctor));
			if(HigherBinEdge == result.end())
			{
				++OverFlow_;
				continue;
			}
			if((HigherBinEdge == result.begin()) && (ptFnComp_(HigherBinEdge->Bin_,FstBin_)))
			{
				++UnderFlow_;
				continue;
			}
			if(ptFnCompInBin_)
			{
				if(!(HigherBinEdge->Count_))
				{
					HigherBinEdge->BinMax_ = t;HigherBinEdge->BinMin_ = t;
				}
				else
				{
					if(ptFnCompInBin_(t,HigherBinEdge->BinMin_))
					{HigherBinEdge->BinMin_ = t;}
					if(ptFnCompInBin_(HigherBinEdge->BinMax_,t))
					{HigherBinEdge->BinMax_ = t;}
				}
			}
			(HigherBinEdge->Count_)++;
		}
		return FstBin_;
	}

	HIST_EL	DoHistogram(HistResultType& result,int MinCount)
	{
		HIST_EL	FstBin(DoHistogram(result));
		DoCompressHistogram(result,MinCount);
		return FstBin;
	}

	void	GetExtremes(int& UnderFlow,int& OverFlow) const
	{
		UnderFlow	= UnderFlow_;
		OverFlow	= OverFlow_;
	}

	Histogram(const DATATYPECOLL& Data,DBINSTYPECOLL& Bins,HIST_EL (*ptFnExtract)(DATATYPECOLLIT ptr),
		bool (*ptFnComp)(const HIST_EL& first,const HIST_EL& sec),
		bool (*ptFnCompInBin)(const HIST_EL& first,const HIST_EL& sec))
		:Data_(Data),Bins_(Bins),ptFnExtract_(ptFnExtract),ptFnComp_(ptFnComp),ptFnCompInBin_(ptFnCompInBin),
		UnderFlow_(0),OverFlow_(0)
	{}
private:

	void	DoCompressHistogram(HistResultType& result,int MinCount)
	{
		typename HistResultType::const_iterator	const pivOrg(result.begin());
		typename HistResultType::iterator		piv(result.begin());
		typename HistResultType::iterator		pivAux;
		typename HistResultType::const_iterator	const end(result.end());
		HIST_EL							tBinMin;
		HIST_EL							tBinMax;
		int								CurrCount;

		for(;piv != end;++piv)
		{
			if(piv->Count_ >= MinCount)
			{
				continue;
			}
			CurrCount = 0;
			for(pivAux = piv;pivAux != end;++pivAux)
			{
				CurrCount	+= pivAux->Count_;
				if(ptFnCompInBin_)
				{
					if(pivAux == piv)
					{
						tBinMin = piv->BinMin_;
						tBinMax = piv->BinMax_;
					}
					if(ptFnCompInBin_(pivAux->BinMin_,tBinMin))
					{tBinMin = pivAux->BinMin_;}
					if(ptFnCompInBin_(tBinMax,pivAux->BinMax_))
					{tBinMax = pivAux->BinMax_;}
				}

				if(CurrCount < MinCount)
				{
					continue;
				}

				if(ptFnCompInBin_)
				{
					pivAux->BinMin_		= tBinMin;
					pivAux->BinMax_		= tBinMax;
				}
				pivAux->Count_	= CurrCount;
				piv = pivAux;
				break;
			}
		}

		if(!(result.empty()) && ((end-1)->Count_ < MinCount))
		{
			AdjustLast(result,MinCount);
		}

		PurgeElements(result,MinCount);
		
	}

	void	AdjustLast(HistResultType& result,int MinCount)
	{
		typename HistResultType::const_reverse_iterator		piv(result.rbegin());
		typename HistResultType::const_reverse_iterator		const pivOrg(piv);
		typename HistResultType::const_reverse_iterator		const end(result.rend());
		for(;(piv != end) && (piv->Count_ < MinCount); ++piv)
			;
		if(piv == end)
			return;
		typename HistResultType::iterator		const endF(result.end());
		typename HistResultType::iterator		pivF(endF - ((piv-pivOrg) + 1));
		HIST_EL	tBinMin(pivF->BinMin_);
		HIST_EL	tBinMax(pivF->BinMax_);
		int		CurrCount(pivF->Count_);
		pivF->Count_ = 0;
		// to be removed
		++pivF;
		for(;pivF != endF;++pivF)
		{
			CurrCount	+= pivF->Count_;
			if(ptFnCompInBin_)
			{
				if(ptFnCompInBin_(pivF->BinMin_,tBinMin))
				{tBinMin = pivF->BinMin_;}
				if(ptFnCompInBin_(tBinMax,pivF->BinMax_))
				{tBinMax = pivF->BinMax_;}
			}
		}

		pivF = result.end() - 1;
		if(ptFnCompInBin_)
		{
			pivF->BinMin_		= tBinMin;
			pivF->BinMax_		= tBinMax;
		}
		pivF->Count_	= CurrCount;
	}

	void	SetUpResultBins(HistResultType& result)
	{
		result.clear();
		std::sort(Bins_.begin(),Bins_.end(),ptFnComp_);	
		FstBin_	= Bins_[0];
		typename DBINSTYPECOLL::const_iterator	piv(Bins_.begin() + 1);
		typename DBINSTYPECOLL::const_iterator	const end(Bins_.end());
		for(;piv != end;++piv)
		{
			result.push_back(HistResultElemType(*piv));
		}
	}
	void	PurgeElements(HistResultType& result,int MinCount)
	{
		typename HistResultType::iterator	newEnd(std::remove_if(result.begin(),result.end(), HistEmptyCellRemoveFunctor<HIST_EL>(MinCount)));
		if(newEnd != result.end())
		{
			result.erase(newEnd,result.end());
		}
	}
};


/*
This pretends to be a general class for drawing from an unidimensional distribution
*/

template<typename HIST_EL>
class DrawFrom1dHist
{
private:

	struct	_DrawTemp
	{
		double		Bin_;
		double		Prob_;

		_DrawTemp(double Bin,double Prob)
			:Bin_(Bin),Prob_(Prob)
		{}

		bool	operator<(const _DrawTemp& rhs) const
		{
			return	(Prob_ < rhs.Prob_);
		}
	};

	typedef	std::vector<_DrawTemp>	DrawSupportStruct;	

	DrawSupportStruct	ProbBins_;
	PwSCoreRandGen		RandGen_;

	double				_DoGetSample(double r,int& Index,double& delta)
	{
		typename DrawSupportStruct::const_iterator	const beg(ProbBins_.begin());
		typename DrawSupportStruct::const_iterator	const end(ProbBins_.end());

		typename DrawSupportStruct::const_iterator	vIndex(std::lower_bound(beg,end,_DrawTemp(0.0,r)));

		if(vIndex == end) --vIndex;	

		Index = static_cast<int>(vIndex - beg);
		
		const double v0(ProbBins_[Index-1].Prob_);
		const double d(ProbBins_[Index].Prob_ - v0);

		delta = (r - v0)/d;
		
		return LinearInterpolation(ProbBins_[Index-1].Bin_,ProbBins_[Index].Bin_,delta);
	}


public:
	typedef	std::vector<HistResultElem<HIST_EL> >	HistIN_Type;

	DrawFrom1dHist(HIST_EL fstBin,const HistIN_Type& HistIn,
		double	(*ptFnExtract)(const HIST_EL& d),
		int RanSeed= -1)
	{
		typename HistIN_Type::const_iterator		piv(HistIn.begin());
		typename HistIN_Type::const_iterator		const end(HistIn.end());

		ProbBins_.push_back(_DrawTemp(ptFnExtract(fstBin),0.0));		

		double Acumul(0.0);

		for(;piv != end;++piv)
		{
			Acumul += static_cast<double>(piv->Count_);
			ProbBins_.push_back(_DrawTemp(ptFnExtract(piv->Bin_),Acumul));		
		}

		typename DrawSupportStruct::iterator				pivSec(ProbBins_.begin());
		typename DrawSupportStruct::const_iterator		const endSec(ProbBins_.end());

		const double	div((endSec-1)->Prob_);				
		for(;pivSec != endSec; ++pivSec)
		{
			pivSec->Prob_ /= div;
		}

		ProbBins_[0].Prob_ = -1.0e-8;
		RandGen_.RandInit(RanSeed);
	}

	inline	void		InitInternalRandGen(int Seed)
	{
		RandGen_.RandInit(Seed);
	}
	inline	double		GetNextUnifDraw(void)
	{
		return RandGen_.RandDouble();
	}
	inline	double		GetSample(void)
	{
		int		DummyI;
		double	DummyD;	
		double	r(RandGen_.RandDouble());
		return _DoGetSample(r,DummyI,DummyD);
	}

	inline	double		GetSample(int& Index,double& delta)
	{
		double	r(RandGen_.RandDouble());
		return _DoGetSample(r,Index,delta);	
	}

	inline	double		GetSample(double NormRand)
	{
		int		DummyI;
		double	DummyD;
		return _DoGetSample(NormRand,DummyI,DummyD);
	}

	inline	double		GetSample(double NormRand,int& Index,double& delta)
	{
		return _DoGetSample(NormRand,Index,delta);
	}
};


} //End of namespace Zeus


#endif //ZEUS_HISTOGRAMH

