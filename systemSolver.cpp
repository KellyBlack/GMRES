
/* *********************************************************************************
 *
 * Copyright (c) 2014, Kelly Black
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * ********************************************************************************* */



#include "poisson.h"
#include "solution.h"
#include "preconditioner.h"
#include "matrix.h"
//#include "gmres.h"
#include "GMRES.h"

#include <iostream>
#include <cmath>

#define DERIVPOWER 50.0

int main(int argc,char **argv)
{

    Poisson *elliptical = new Poisson;
	Solution *x = new Solution(NUMBER);
	Solution *b = new Solution(NUMBER);
	Preconditioner pre(NUMBER);
	Matrix upper(41,40);

	int lupe;

	int restart = 10;
	int maxIt = 41;
	double tol = 1.0E-8;

	for(lupe=0;lupe<=NUMBER;++lupe)
		{
			double xgrid = elliptical->getX(lupe);
			//b(lupe) = -81.0*M_PI*M_PI*sin(9.0*M_PI*xgrid);
			(*b)(lupe) = 90.0*pow(xgrid,8.0)-2.0;
			(*x)(lupe) = 0.0; //pow(xgrid,10.0)/90.0 - xgrid*xgrid/91;
		}
	(*b)(0) = 0.0;
	(*b)(NUMBER) = -0.0;

	int result= GMRES(elliptical,x,b,maxIt,restart,NUMBER,tol);
	//int result = GMRES(elliptical,x,b,pre,upper,restart,maxIt,tol);
	//Solution y = elliptical*x;

	//std::cout << "max iter: " << restart << " max it: " << maxIt << " result: " << result << std::endl;
#define SOLUTION
#ifdef SOLUTION
//std::cout << "x,approx,true," << result << std::endl;
	for(lupe=0;lupe<=NUMBER;++lupe)
		{
			//if(lupe%5 == 0)
			//	std::cout << std::endl;
			double xgrid = elliptical->getX(lupe);
			std::cout << xgrid << "," 
					  << (*x)(lupe) << "," 
				//<< 56.0*pow(xgrid,6.0)-2.0
								<< (pow(xgrid,10.0)-xgrid*xgrid)
				//<< sin(9.0*M_PI*xgrid) 
								<< std::endl;
		}
#endif

	/*
	for(lupe=0;lupe<=NUMBER;++lupe)
		{
			double xgrid = elliptical.getX(lupe);
			b(lupe) = pow(xgrid,60.0);
		}
	y = elliptical*b;
	for(lupe=0;lupe<=NUMBER;++lupe)
		{
			if(lupe%5 == 0)
				std::cout << std::endl;
			double xgrid = elliptical.getX(lupe);
			std::cout << fabs(y(lupe)-60.0*59.0*pow(xgrid,58.0)) << "  ";
		}
	std::cout << std::endl;
	*/

	return(1);
}

