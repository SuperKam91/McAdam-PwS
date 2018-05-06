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
 *  Copyright (C) 2005, 2013, Pedro Carvalho
 *  Author: Pedro Carvalho
 */

//----------------------------------
#ifndef QA_ASSESS_CONTOURS
#define QA_ASSESS_CONTOURS

#define QACONTOURSINVALIDVALUE	-1.63750e+30
#define QACONTOURSMARG68		0.68
#define QACONTOURSMARG95		0.95

struct QAContMargType
{
	double probVal_;
	double Val_;

	QAContMargType(void):
		probVal_(QACONTOURSINVALIDVALUE),Val_(QACONTOURSINVALIDVALUE)
	{}
	QAContMargType(double probVal):
		probVal_(probVal),Val_(QACONTOURSINVALIDVALUE)
	{}
};

#ifdef WIN32
// Gnu Fortran
#define QA_CONTOURS_FUNCTION	__qual_assess_contours_MOD_quality_assess_contours	
#else
// Intel Fortran
#define QA_CONTOURS_FUNCTION	qual_assess_contours_mp_quality_assess_contours_
#endif
//
inline bool SortByProb(const QAContMargType& first,const QAContMargType& sec)
{return first.probVal_ > sec.probVal_ ;}
//
void	myQA_ContAssess(int nscales,double input_theta_s,double input_cy5r500, double min_rs,
		  double max_rs, double min_ytot,double  max_ytot,const double* image_contours,double*  qa_results,int peakEstimator=0);

#ifdef __cplusplus
extern "C" {
#endif

      void QA_CONTOURS_FUNCTION(const int* nscales,
		  const double* input_theta_s,const double* input_cy5r500, const double* min_rs,
		  const double* max_rs, const double* min_ytot,const double*  max_ytot,const double* image_contours,double*  qa_results);

#ifdef __cplusplus
}
#endif
#endif //QA_ASSESS_CONTOURS
//----------------------------------

