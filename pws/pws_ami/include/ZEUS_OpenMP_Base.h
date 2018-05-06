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


#ifndef ZEUS_OPENMPSUPPORTH
#define ZEUS_OPENMPSUPPORTH

#ifdef _OPENMP
#include <omp.h>
#endif

inline bool openmp_enabled()
  {
#ifdef _OPENMP
  return true;
#else
  return false;
#endif
  }

inline int openmp_max_threads()
  {
#ifdef _OPENMP
  return omp_get_max_threads();
#else
  return 1;
#endif
  }

inline int openmp_thread_num()
  {
#ifdef _OPENMP
  return omp_get_thread_num();
#else
  return 0;
#endif
  }
#endif //ZEUS_OPENMPSUPPORTH