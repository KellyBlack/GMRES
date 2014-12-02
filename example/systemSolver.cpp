
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
#include "../GMRES.h"

#include <iostream>
#include <cmath>

#define DERIVPOWER 50.0

int main(int argc,char **argv)
{

    Poisson *elliptical = new Poisson;   // The operator to invert.
	Solution *x = new Solution(NUMBER);  // The approximation to calculate.
	Solution *b = new Solution(NUMBER);  // The forcing function for the r.h.s.
	Preconditioner *pre = 
		new Preconditioner(NUMBER);      // The preconditioner for the system.

	int restart = 10;                    // Number of restarts to allow
	int maxIt = 41;                      // Dimension of the Krylov subspace
	double tol = 1.0E-8;                 // How close to make the approximation.

	int lupe;
	for(lupe=0;lupe<=NUMBER;++lupe)
		{
			// initialize the r.h.s to be something we know the
			// solution for. Also, set the initial approximation if
			// wanted.
			double xgrid = elliptical->getX(lupe);
			//b(lupe) = -81.0*M_PI*M_PI*sin(9.0*M_PI*xgrid);
			(*b)(lupe) = 90.0*pow(xgrid,8.0)-2.0;
			(*x)(lupe) = 0.0; //pow(xgrid,10.0)/90.0 - xgrid*xgrid/91;
		}

	// Set the boundary conditions separately.
	(*b)(0) = 0.0;
	(*b)(NUMBER) = -0.0;

	// Find an approximation to the system!
	int result= GMRES(elliptical,x,b,pre,maxIt,restart,tol);

	std::cout << "Iterations: " << result << " residual: " << tol << std::endl;
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

