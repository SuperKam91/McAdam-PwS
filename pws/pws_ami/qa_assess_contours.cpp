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

#include <vector>
#include "ZEUS_PhysicsMath.h"
#include "qa_assess_contours.h"

/*
!! 1 :: HPD_2D     :: dimensionless :: HPD required to enclose the true cluster parameters
!! 2 :: HPD_Y      :: dimensionless :: HPD required to enclose the true cluster Y value (marginalise over theta)
!! 3 :: HPD_THETA  :: dimensionless :: HPD required to enclose the true cluster ThetaS value (marginalise over Y)
!! 4 :: HPD_YTHETA :: dimensionless :: HPD required to enclose the true cluster Y value (on closest grid to true theta)
!! 5 :: Ypeak     :: arcmin^2 :: Peak in Y (peak from maginalised distribution)
!! 6 :: Y_68_low  :: arcmin^2 :: lower limit of the width of HPD which encloses 68 % of the probability when theta is maginalised over
!! 7 :: Y_68_up   :: arcmin^2 :: upper limit of the width of HPD which encloses 68 % of the probability when theta is maginalised over
!! 8 :: Y_95_low  :: arcmin^2 lower limit of the width of HPD which encloses 95 % of the probability when theta is maginalised over
!! 9 :: Y_95_up   :: arcmin^2 upper limit of the width of HPD which encloses 95 % of the probability when theta is maginalised over
!!10 :: ThetaS_peak   :: arcmin :: Peak in ThetaS (peak from maginalised distribution)
!!11 :: ThetaS_68_low :: arcmin :: lower limit of the width of HPD which encloses 68 % of the probability when Y is maginalised over
!!12 :: ThetaS_68_up  :: arcmin  :: upper limit of the width of HPD which encloses 68 % of the probability when Y is maginalised over
!!13 :: ThetaS_95_low :: arcmin :: lower limit of the width of HPD which encloses 95 % of the probability when Y is maginalised over
!!14 :: ThetaS_95_up  :: arcmin :: upper limit of the width of HPD which encloses 95 % of the probability when Y is maginalised over
!!15 :: YTheta_peak   :: arcmin^2 :: Peak in Y at injected value of ThetaS (closest value on grid)
!!16 :: YTheta_68_low :: arcmin^2 :: lower limit of the width of HPD which encloses 68 % of the probability
!!17 :: YTheta_68_up  :: arcmin^2 :: upper limit of the width of HPD which encloses 68 % of the probability
!!18 :: YTheta_95_low :: arcmin^2 :: lower limit of the width of HPD which encloses 95 % of the probability
!!19 :: YTheta_95_up  :: arcmin^2 ::upper limit of the width of HPD which encloses 95 % of the probability
!!
!!20 :: Y_low_axis_prob  :: sum of all points with probabilities >= maximum probability point on the lower limit Y-axis 
!!21 :: Y_high_axis_prob :: sum of all points with probabilities >= maximum probability point on the upper limit Y-axis 
!!22 :: ThetaS_low_axis_prob  :: sum of all points with probabilities >= maximum probability point on the lower limit ThetaS-axis
!!23 :: ThetaS_high_axis_prob :: sum of all points with probabilities >= maximum probability point on the upper limit ThetaS-axis
!!
!! NB HPD set to +10 if input value is off-grid
*/

#define HPD_2D			0
#define HPD_Y			1
#define HPD_THETA		2
#define HPD_YTHETA		3
#define Ypeak			4
#define Y_68_low		5
#define Y_68_up			6
#define Y_95_low		7
#define Y_95_up			8
#define ThetaS_peak		9
#define ThetaS_68_low	10
#define ThetaS_68_up	11
#define ThetaS_95_low	12
#define ThetaS_95_up	13
#define YTheta_peak		14
#define YTheta_68_low	15
#define YTheta_68_up	16
#define YTheta_95_low	17
#define YTheta_95_up	18
// New values
#define Y_low_axis_prob			19
#define Y_high_axis_prob		20
#define ThetaS_low_axis_prob	21
#define ThetaS_high_axis_prob	22


void	myQA_ContAssess(int nscales,double input_theta_s,double input_cy5r500, double min_rs,
		  double max_rs, double min_ytot,double  max_ytot,const double* image_contours,double*  qa_results,int peakEstimator)
{
// initialisation
	std::vector<QAContMargType> Y_Marg(nscales,QAContMargType(0.0));
	std::vector<QAContMargType> ThMarg(nscales,QAContMargType(0.0));
	std::vector<QAContMargType> YThTrue(nscales,QAContMargType(0.0));

	for(int i =Ypeak;i<(YTheta_95_up + 1);i++)
	{*(qa_results + i)= QACONTOURSINVALIDVALUE;}
	qa_results[HPD_2D] = qa_results[HPD_Y] = qa_results[HPD_THETA] = qa_results[HPD_YTHETA] = 10.0; 

	const double	deltaY((max_ytot - min_ytot)/static_cast<double>(nscales));
	const double	deltaTh((max_rs - min_rs)/static_cast<double>(nscales));
	int				ThTrueIndex(-1);
	int				YTrueIndex(-1);

//	Out of boundaries probability
	double Y_low_axis_MaxProb(-1.0);
	double Y_high_axis_MaxProb(-1.0);
	double ThetaS_low_axis_MaxProb(-1.0);
	double ThetaS_high_axis_MaxProb(-1.0);

	for(int i=0;i<nscales;++i)
	{
		if(Y_low_axis_MaxProb < *(image_contours + i)) Y_low_axis_MaxProb = *(image_contours + i);
		if(Y_high_axis_MaxProb < *(image_contours + ((nscales-1) * nscales) + i)) Y_high_axis_MaxProb = *(image_contours + ((nscales-1) * nscales) + i);
		if(ThetaS_low_axis_MaxProb < *(image_contours + (i * nscales))) ThetaS_low_axis_MaxProb = *(image_contours + (i * nscales));
		if(ThetaS_high_axis_MaxProb < *(image_contours + (i * nscales) + (nscales-1))) ThetaS_high_axis_MaxProb = *(image_contours + (i * nscales) + (nscales-1));	
	}
	qa_results[Y_low_axis_prob] = qa_results[Y_high_axis_prob] = qa_results[ThetaS_low_axis_prob] = qa_results[ThetaS_high_axis_prob] = 0.0;
	double	tCurrProb;
	for(int i=0;i<nscales;++i)
	{
		int rowOff(nscales * i);

		for(int j=0;j<nscales;++j)
		{
			tCurrProb = *(image_contours + (rowOff + j));
			if(tCurrProb < Y_low_axis_MaxProb) qa_results[Y_low_axis_prob] += tCurrProb;
			if(tCurrProb < Y_high_axis_MaxProb) qa_results[Y_high_axis_prob] += tCurrProb;
			if(tCurrProb < ThetaS_low_axis_MaxProb) qa_results[ThetaS_low_axis_prob] += tCurrProb;
			if(tCurrProb < ThetaS_high_axis_MaxProb) qa_results[ThetaS_high_axis_prob] += tCurrProb;
		}
	}
// 2D HPD
	if(!((input_theta_s < min_rs) || (input_theta_s >= (max_rs-(deltaTh/2.0)))))
	{
		ThTrueIndex		= Zeus::toInt(((input_theta_s - min_rs)/deltaTh) + 0.5);
//		ThTrueIndex = Zeus::toInt(((input_theta_s - min_rs) / deltaTh));

		if(ThTrueIndex >= nscales) --ThTrueIndex;
		if(!((input_cy5r500 < min_ytot) || (input_cy5r500 >= (max_ytot-(deltaY/2.0)))))
		{
			YTrueIndex = Zeus::toInt(((input_cy5r500 - min_ytot)/deltaY) + 0.5);
//			YTrueIndex = Zeus::toInt(((input_cy5r500 - min_ytot) / deltaY));
			if(YTrueIndex >= nscales) --YTrueIndex;
			const double t2DProbValue(*(image_contours + (YTrueIndex * nscales) + ThTrueIndex));
			double	HDP2DAcc(0.0);
			double	temp;
			for(int i=0;i<nscales;++i)
			{
				int rowOff(nscales * i);
				for(int j=0;j<nscales;++j)
				{
					if((temp=*(image_contours + rowOff + j)) > t2DProbValue)
						HDP2DAcc += temp;
				}		
			}
			qa_results[HPD_2D] = HDP2DAcc;
		}
	}
//
// making the 1d distributions
	for(int i=0;i<nscales;++i)
	{
		int rowOff(nscales * i);
		Y_Marg[i].Val_ = min_ytot + (i * deltaY);
		for(int j=0;j<nscales;++j)
		{
			Y_Marg[i].probVal_	+= *(image_contours + rowOff + j);
			ThMarg[j].probVal_	+= *(image_contours + rowOff + j);
		}		
	}
	for(int i=0;i<nscales;++i)
	{
		ThMarg[i].Val_ = min_rs + (i * deltaTh);
	}

	double	ThMargInTrueVal(1.0e32);
	double	YMargInTrueVal(1.0e32);

	if(ThTrueIndex>=0)
	{
		ThMargInTrueVal			= ThMarg[ThTrueIndex].probVal_;
		qa_results[HPD_THETA]	= 0.0;
	}
	if(YTrueIndex>=0)
	{
		YMargInTrueVal			= Y_Marg[YTrueIndex].probVal_;
		qa_results[HPD_Y]		= 0.0;
	}

	std::sort(Y_Marg.begin(),Y_Marg.end(),SortByProb);
	std::sort(ThMarg.begin(),ThMarg.end(),SortByProb);
//
	qa_results[Y_68_low] = max_ytot + deltaY;
	qa_results[Y_68_up] = min_ytot - deltaY;
	qa_results[Y_95_low] = max_ytot + deltaY;
	qa_results[Y_95_up] = min_ytot - deltaY;
	double	Y_Prob(0.0);
	qa_results[ThetaS_68_low]= max_rs + deltaTh;
	qa_results[ThetaS_68_up]= min_rs - deltaTh;
	qa_results[ThetaS_95_low]= max_rs + deltaTh;
	qa_results[ThetaS_95_up]= min_rs - deltaTh;
	if(peakEstimator)
	{
		qa_results[Ypeak]= 0.0;
		qa_results[ThetaS_peak]= 0.0;
	}
	else
	{
		qa_results[Ypeak]= Y_Marg[0].Val_;
		qa_results[ThetaS_peak]= ThMarg[0].Val_;
	}
	double	ThProb(0.0);

	for(int i=0;i<nscales;++i)
	{
		if(Y_Marg[i].probVal_ >= YMargInTrueVal)
		{qa_results[HPD_Y] += Y_Marg[i].probVal_;}
		if(ThMarg[i].probVal_ >= ThMargInTrueVal)
		{qa_results[HPD_THETA] += ThMarg[i].probVal_;}
		if(peakEstimator)
		{
			qa_results[Ypeak] +=		(Y_Marg[i].probVal_ * Y_Marg[i].Val_);
			qa_results[ThetaS_peak] +=	(ThMarg[i].probVal_ * ThMarg[i].Val_);
		}
		if(Y_Prob < QACONTOURSMARG95)
		{
			if(Y_Marg[i].Val_ < qa_results[Y_95_low])
			{qa_results[Y_95_low] = Y_Marg[i].Val_;}
			if(Y_Marg[i].Val_ > qa_results[Y_95_up])
			{qa_results[Y_95_up] = Y_Marg[i].Val_;}
			if(Y_Prob < QACONTOURSMARG68)
			{
				if(Y_Marg[i].Val_ < qa_results[Y_68_low])
				{qa_results[Y_68_low] = Y_Marg[i].Val_;}
				if(Y_Marg[i].Val_ > qa_results[Y_68_up])
				{qa_results[Y_68_up] = Y_Marg[i].Val_;}
			}
		}
		if(ThProb < QACONTOURSMARG95)
		{
			if(ThMarg[i].Val_ < qa_results[ThetaS_95_low])
			{qa_results[ThetaS_95_low] = ThMarg[i].Val_;}
			if(ThMarg[i].Val_ > qa_results[ThetaS_95_up])
			{qa_results[ThetaS_95_up] = ThMarg[i].Val_;}
			if(ThProb < QACONTOURSMARG68)
			{
				if(ThMarg[i].Val_ < qa_results[ThetaS_68_low])
				{qa_results[ThetaS_68_low] = ThMarg[i].Val_;}
				if(ThMarg[i].Val_ > qa_results[ThetaS_68_up])
				{qa_results[ThetaS_68_up] = ThMarg[i].Val_;}
			}
		}
		Y_Prob += Y_Marg[i].probVal_;
		ThProb += ThMarg[i].probVal_;
	}	
//
	if(ThTrueIndex >= 0)
	{
		double t(0.0);
		for(int i=0;i<nscales;++i)
		{
			YThTrue[i].Val_			= min_ytot + i * deltaY;
			YThTrue[i].probVal_		= *(image_contours + (nscales * i) + ThTrueIndex);
			t						+= YThTrue[i].probVal_;
		}
		for(int i=0;i<nscales;++i)
		{YThTrue[i].probVal_ /= t;}

		double	YThInTrueVal(1.0e32);

		if(YTrueIndex>=0)
		{
			YThInTrueVal			= YThTrue[YTrueIndex].probVal_;
			qa_results[HPD_YTHETA]	= 0.0;
		}
		std::sort(YThTrue.begin(),YThTrue.end(),SortByProb);

		qa_results[YTheta_68_low]	= max_ytot + deltaY;
		qa_results[YTheta_68_up]	= min_ytot - deltaY;
		qa_results[YTheta_95_low]	= max_ytot + deltaY;
		qa_results[YTheta_95_up]	= min_ytot - deltaY;
		qa_results[YTheta_peak]		= YThTrue[0].Val_;
		double	YTh_Prob(0.0);
		for(int i=0;i<nscales;++i)
		{
			if(YThTrue[i].probVal_ >= YThInTrueVal)
			{qa_results[HPD_YTHETA] += YThTrue[i].probVal_;}
			if(YTh_Prob < QACONTOURSMARG95)
			{
				if(YThTrue[i].Val_ < qa_results[YTheta_95_low])
				{qa_results[YTheta_95_low] = YThTrue[i].Val_;}
				if(YThTrue[i].Val_ > qa_results[YTheta_95_up])
				{qa_results[YTheta_95_up] = YThTrue[i].Val_;}
				if(YTh_Prob < QACONTOURSMARG68)
				{
					if(YThTrue[i].Val_ < qa_results[YTheta_68_low])
					{qa_results[YTheta_68_low] = YThTrue[i].Val_;}
					if(YThTrue[i].Val_ > qa_results[YTheta_68_up])
					{qa_results[YTheta_68_up] = YThTrue[i].Val_;}
				}
			}
			YTh_Prob += YThTrue[i].probVal_;
		}
	}
//
	return;
}
