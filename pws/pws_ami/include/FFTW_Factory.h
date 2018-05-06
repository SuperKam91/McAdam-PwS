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

#ifndef FFTW_factoryH
#define FFTW_factoryH
//---------------------------------------------------------------------------
#include "FFTW_Storage.h"
#define FFTW_FACTORY    (FFTW_Factory::AccessPoint)

// class FFTW_Device;
class FFTW_Factory
{
public:
    // Singleton logic
    static FFTW_Factory& AccessPoint(void)
    {
      if(!AccessPoint_){
        AccessPoint_ = new FFTW_Factory;
      }
      return *AccessPoint_;
    }
    // end of singleton logic

private:
    FFTW_Factory() {;}

    // Singleton protection
    static FFTW_Factory*    AccessPoint_;   // Singleton access point
    FFTW_Factory(const FFTW_Factory&);
    FFTW_Factory operator=(const FFTW_Factory&);
    ~FFTW_Factory() {}

};
#endif //FFTW_factoryH
