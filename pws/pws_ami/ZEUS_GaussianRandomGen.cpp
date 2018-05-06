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

#include "ZEUS_GaussianRandomGen.h"

namespace Zeus
{

const double JSkillingRandGen::RR			= 4611686018427387904.0;	 // 2^62
const double JSkillingRandGen::SHIFT32		= 1.0 / 4294967296.0;        // 2^-32
const double JSkillingRandGen::HalfPi		= 1.57079632679489661922;    // pi/2
const double JSkillingRandGen::TwoxPi		= 6.28318530717958647688;    // 2*pi
const double JSkillingRandGen::SqrPi2		= 1.25331413731550025120;    // sqrt(pi/2)
const double JSkillingRandGen::Sqr2Pi		= 2.50662827463100050240;    // sqrt(2*pi)
const double JSkillingRandGen::LogSqr2Pi	= 0.91893853320467274177;    // log(sqrt(2*pi))


void		JSkillingRandGen::RandomInit(int seed)
{

	R32BitsType		seed32(static_cast<R32BitsType>(seed));

	Rand[0] = Rand[1] = 0;                     // 64-bit counter
    Rand[2] = 1013904223 + 1664525 * seed32;   // sticky offset
    Rand[3] = 1013904223 + 1664525 * Rand[2];  // extra random integer
}

int 	JSkillingRandGen::RanInt(void)
{
    R32BitsType  m, n;   // 64-bit register
    R32BitsType  u, v, w;
    int       i;
// 64-bit counter, for hashing
    if( !(Rand[0] ++) )
    {
        Rand[1] ++;
// occasional update offset in Rand[2]
        i = Rand[2] >> 1;                    // [0...2147483647]
        i -= i / 24683721;                   // [0...2147483561]
        i++;                                 // [1...2147483562]
        i = 40014 * i - 2147483563 * (i / 53668);
        if( i < 0 )                          // (40014 * i) % 2147483563
            i += 2147483563;                 // [1...2147483562]
        i--;                                 // [0...2147483561]
        i += i / 24683720;                   // [0...2147483647] (with holes)
        Rand[2] = i << 1;                    // even
        if( 1013904223 + 1664525 * i < 0 )   // chance
            Rand[2]++;                       // even or odd
    }
// Two double steps of 64-bit hash
    n = Rand[0] + Rand[2];
    m = Rand[1] + Rand[2];
    w = n ^ 0xbaa96887;
    v = w >> 16;
    w &= 0xffff;
    u = (v - w) * (v + w);
    m ^= (((u >> 16) | (u << 16)) ^ 0xb4f0c4a7) + w * v;
    w = m ^ 0x1e17d32c;
    v = w >> 16;
    w &= 0xffff;
    u = (v - w) * (v + w);
    n ^= (((u >> 16) | (u << 16)) ^ 0x178b0f3c) + w * v;
    w = n ^ 0x03bcdc3c;
    v = w >> 16;
    w &= 0xffff;
    u = (v - w) * (v + w);
    m ^= (((u >> 16) | (u << 16)) ^ 0x96aa3a59) + w * v;
    w = m ^ 0x0f33d1b2;
    v = w >> 16;
    w &= 0xffff;
    u = (v - w) * (v + w);
    n ^= (((u >> 16) | (u << 16)) ^ 0xaa5835b9) + w * v;
    Rand[3] = m;
    return  n;
}

float 	JSkillingRandGen::RanFloat(void)
{
    R32BitsType u((static_cast<R32BitsType>(RanInt()) & 0xfffffe00) ^ 0x00000100);  // switch lowest (2^-24) bit ON
	return static_cast<float>(u) * static_cast<float>(SHIFT32);
}

double  JSkillingRandGen::Random(void)
{
    R32BitsType  hi,lo;
    hi = static_cast<R32BitsType>(RanInt());
    lo = (Rand[3] & 0xfffff000) ^ 0x00000800;

    return  (static_cast<double>(hi) + static_cast<double>(lo) * SHIFT32) * SHIFT32;
}

}

