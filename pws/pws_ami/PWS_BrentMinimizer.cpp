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



//--------------------------------------

#include "PWS_BrentMinimizer.h"

//--------------------------------------

void	BrentLineMinimizer::Mnbrak(MinArrayAtomType& ax,MinArrayAtomType&  bx,MinArrayAtomType&  cx,MinArrayAtomType&  fa,MinArrayAtomType&  fb,MinArrayAtomType&  fc)
{
	MinArrayAtomType ulim,u,q,r,fu,dum,tmp;

	fa = Funct_1_dimen(ax);
	fb = Funct_1_dimen(bx);
	if(fb > fa){SHFT(dum, ax, bx, dum);SHFT(dum, fa, fb, dum);}
	cx = bx + GOLD * (bx - ax);
	fc = Funct_1_dimen(cx);
	while(fb > fc)
	{
		r = (bx - ax) * (fb - fc);
		q = (bx - cx) * (fb - fa);
		tmp = q - r;
		if(std::abs(tmp) < TINY) {tmp = ((tmp < 0) ? -TINY:TINY);}
		tmp = ((MinArrayAtomType)2.0) * tmp;
		u = bx - ((((bx - cx) * q - (bx - ax) * r)) / tmp);
		ulim = bx + GLIMIT * (cx - bx);
		if(((u > cx) && (bx > u)) || ((u > bx) && (cx > u))){
			fu = Funct_1_dimen(u);
			if(fu < fc){
				ax=bx;bx=u;fa=fb;fb=fu;
				return;
			}
			else if(fu > fb){
				cx=u;fc=fu;
				return;
			}
			u = cx + GOLD * (cx - bx);
			fu = Funct_1_dimen(u);
		}
		else if(((u > ulim) && (cx > u)) || ((u > cx) && (ulim > u))) {
				fu = Funct_1_dimen(u);
				if(fu<fc){
					SHFT(bx,cx,u,cx+ GOLD * (cx - bx));
					SHFT(fb,fc,fu,Funct_1_dimen(u));
				}
			}
			else if(((u >= ulim) && (ulim >= cx)) || ((cx >= ulim) && (ulim >= u))){
				u = ulim;
				fu =  Funct_1_dimen(u);
			}
				else{
					u = cx + GOLD * (cx -bx);
					fu =  Funct_1_dimen(u);
				}
		SHFT(ax,bx,cx,u);
		SHFT(fa,fb,fc,fu);
	}
}

int		BrentLineMinimizer::Brent(MinArrayAtomType brAx,MinArrayAtomType brBx,MinArrayAtomType brCx,MinArrayAtomType tol,MinArrayAtomType& xResult,MinArrayAtomType& fxResult)
{
	int	iter;

	MinArrayAtomType e(0.0),a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm; 

	a = ((brAx < brCx)?brAx:brCx);
	b = ((brAx > brCx)?brAx:brCx);
	x = w = v = brBx;
	fw = fv = fx = Funct_1_dimen(x);
	for(iter=1;iter<=MaxIter_;++iter){
		xm = (a + b)/((MinArrayAtomType)2.0);
		tol2 = ((MinArrayAtomType)2.0) * (tol1=(tol*std::abs(x)+((MinArrayAtomType)ZEPS)));
		if(std::abs(x - xm) <= ((tol2>=tol?tol2:tol) - ((b - a)/((MinArrayAtomType)2.0)))){
			xResult = x;fxResult = fx;
			return true;
		}
		if(std::abs(e) > tol1){
			r = (x - w) * (fx - fv);q = (x - v) * (fx - fw);
			p = (x - v) * q - (x - w) * r;q = ((MinArrayAtomType)2.0) * (q - r);
			if(q > 0.0) {p = -p;}
			q = std::abs(q);etemp = e;e = d;
			if((std::abs(p) >= std::abs((q*etemp)/((MinArrayAtomType)2.0))) || (p<=q*(a-x)) || (p>=q*(b-x)))
				{d = CGOLD * (e = ((x>=xm)?a-x:b-x));}
			else{
				d = p / q;u = x + d;
				if(((u-a) < tol2) ||((b - u) < tol2)) {d = ((xm>=x)?tol1:-tol1);}
			}
		}
		else{d = CGOLD * (e = ((x>=xm)?a-x:b-x));}
		u = ((std::abs(d) >= tol1)?x+d:x + (d>=0?tol1:-tol1));

		fu = Funct_1_dimen(u);

		if(fu <= fx){
			if(u>=x)
			{a = x;}
			else{b = x;}
			Shft(v,w,x,u);Shft(fv,fw,fx,fu);
		}
		else{
			if(u < x)
			{a = u;}
			else{b = u;}
			if((fu <= fw) || (w == x)){v = w;w = u;fv = fw;fw = fu;}
			else if((fu <= fv) || (v == x) || (v == w)){v = u;fv = fu;}
		}
	}
	xResult = x;fxResult = fx;
	return false;
}

int		BrentLineMinimizer::FindLineMinimum(const MinArrayType& searchDir,MinArrayType& Init_Final,MinArrayAtomType& FunctValue)
{
	MinArrayAtomType	ax,xx,bx,fa,fb,fc,FResTmp,xRes;

	SetSearchDirs(searchDir);	
	SetZeroArgs(Init_Final);
	ax = 0.0;xx = 1.0;
	Mnbrak(ax,xx,bx,fa,fb,fc);

	if(!(Brent(ax,xx,bx,BrentTol_,xRes,FResTmp)))
		return false;

	FunctValue = FResTmp;
	GetFinalResult(xRes,Init_Final);
	return true;
}
