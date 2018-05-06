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
 *
 *	This class was built based on the code *Nested Sampler* from John Skilling, only with very few minor changes. 
 */


//---------------------------------------------------------------------------
#ifndef BSNAK_GAUSSIANGENH
#define BSNAK_GAUSSIANGENH
//---------------------------------------------------------------------------
#include <time.h>
#include <string>

#ifdef WIN32
#include <tchar.h>
#endif //WIN32

#include <complex>
#include <math.h>
#include <typeinfo>
#include <exception>
#include <locale>
#include "Threads.h"
#include "SmartPtr.h"
#include "TypeManip.h"
#include "static_check.h"
#include "randomc.h"

#include "ZEUS_Strings.h"

#define	NOTSUPPORTEDCOMPILER	NOT_SUPPORTED_COMPILER

namespace Zeus
{

namespace Private
{

template<int v>
struct R32BitsHelper
{
	TYPES_STATIC_CHECK(false, NOT_SUPPORTED_COMPILER);
};

template <> struct R32BitsHelper<2>
{
	typedef		unsigned long			Type;
};

template <> struct R32BitsHelper<4>
{
	typedef		unsigned int			Type;
};

template <> struct R32BitsHelper<8>
{
	typedef		unsigned short int		Type;
};

}

class	JSkillingRandGen
{
	typedef Private::R32BitsHelper<sizeof(int)>::Type		R32BitsType;
	TYPES_STATIC_CHECK(sizeof(R32BitsType) == 4,NOT_SUPPORTED_COMPILER);
	typedef R32BitsType Rand_t[4];

	static const double SHIFT32;    // 2^-32
	static const double HalfPi;     // pi/2
	static const double TwoxPi;     // 2*pi
	static const double SqrPi2;     // sqrt(pi/2)
	static const double Sqr2Pi;     // sqrt(2*pi)
	static const double LogSqr2Pi;	// log(sqrt(2*pi))
    static const double RR;			// 2^62

	Rand_t			Rand;			// Generator state
public:

	inline JSkillingRandGen(int seed)
	{RandomInit(seed);}
	void			RandomInit(int	seed);		// Input seed: +ve = value, -ve = time seed; return seed  
	double			Random(void);		// Output random sample, inside (0.0, 1.0)
private:
	int				RanInt(void);			// Output random sample, in [-2^31, 2^31)
	float			RanFloat(void);
};

template<class BASE_RNG>
class	RandGenBase
{
	BASE_RNG		BaseGen_;
	int				GaussAvailable_;
	double			Gauss2Value_;

	inline RandGenBase(const RandGenBase& );
	inline RandGenBase& operator=(const RandGenBase& rhs);

public:
	RandGenBase(int seed = -1)
		:BaseGen_((seed==-1)?static_cast<int>(time(NULL)):seed),GaussAvailable_(0)
	{}

	inline	void	RandInit(int seed = -1)
	{
		if(seed == -1) seed = static_cast<int>(time(NULL));
		BaseGen_.RandomInit(seed);
	}
	inline	double	RandDouble(void)
	{return BaseGen_.Random();}

	double	RandGauss(void)
	{
		if(GaussAvailable_)
		{
			GaussAvailable_ = 0;
			return Gauss2Value_;
		}
		double	v1,v2,sqr12;
		do
		{
			v1 = 2.0 * BaseGen_.Random() - 1.0;
			v2 = 2.0 * BaseGen_.Random() - 1.0;
			sqr12 = v1*v1+v2*v2;
		}
		while((sqr12==0.0) || (sqr12 >= 1.0));
		const double fac(std::sqrt(-2.0 * std::log(sqr12)/sqr12));
		Gauss2Value_	= v1 * fac;
		GaussAvailable_	= 1;
		return v2 * fac;
	}
};


// This template class combines two different random number generators
// for improved randomness. R1 and R2 are any two different random number
// generator classes.
template <class RG1, class RG2>
class TRandomCombined : private RG1, private RG2 {
public:
   TRandomCombined(int seed) : RG1(seed), RG2(seed+1) {};

   void RandomInit(int seed) {        // re-seed
      RG1::RandomInit(seed);
      RG2::RandomInit(seed+1);
   }

   double Random() {
      double r = RG1::Random() + RG2::Random();
      if (r >= 1.) r -= 1.;
      return r;
   }

   int IRandom(int min, int max){       // output random integer
      // get integer random number in desired interval
      int iinterval = max - min + 1;
      if (iinterval <= 0) return 0x80000000; // error
      int r = int(iinterval * Random()); // truncate
      if (r >= iinterval) r = iinterval-1;
      return min + r;
   }
};


typedef TRandomCombined<CRandomMersenne,CRandomMother>	CurrGenerator;
//typedef JSkillingRandGen								CurrGenerator;

typedef RandGenBase<CurrGenerator>						MultNestRandGen;
typedef RandGenBase<CurrGenerator>						PwSCoreRandGen;

}

#endif //BSNAK_GAUSSIANGENH

