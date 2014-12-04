
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
	int maxIt   = 800;                   // Dimension of the Krylov subspace
	double tol  = 1.0E-8;                // How close to make the approximation.

	int row;
	int col;
	for(row=0;row<=NUMBER;++row)
		{
			for(col=1;col<NUMBER;++col)
				{
					// initialize the r.h.s to be something we know the
					// solution for. Also, set the initial approximation if
					// wanted.
					double xgrid = elliptical->getX(row);
					double ygrid = elliptical->getX(col);
					(*b)(row,col) = -2.0*(1.0-xgrid*xgrid)-2.0*(1.0-ygrid*ygrid);
					//(*x)(row,col) = 0.0;
				}

			// Set the top and bottom boundary conditions. (This is
			// redundant but may need to be changed for different
			// forcing functions above.
			(*b)(row,0)      = 0.0;
			(*b)(row,NUMBER) = 0.0;
		}

	// Set the left and right boundary conditions separately.
	for(col=0;col<=NUMBER;++col)
		{
			(*b)(0,col)      = 0.0;
			(*b)(NUMBER,col) = 0.0;
		}

	// Find an approximation to the system!
	int result= GMRES(elliptical,x,b,pre,maxIt,restart,tol);

	std::cerr << "Iterations: " << result << " residual: " << tol << std::endl;
#define SOLUTION
#ifdef SOLUTION
//std::cout << "x,approx,true," << result << std::endl;
	for(row=0;row<=NUMBER;++row)
		for(col=0;col<=NUMBER;++col)
		{
			double xgrid = elliptical->getX(row);
			double ygrid = elliptical->getX(col);
			std::cout << xgrid << "," << ygrid << "," 
					  << (*x)(row,col) << "," 
					  << (1.0-xgrid*xgrid)*(1.0-ygrid*ygrid) << "," 
					  << (*b)(row,col) 
					  << std::endl;
		}
#endif


	return(1);
}

