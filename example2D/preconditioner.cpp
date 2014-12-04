
/** *********************************************************************************
 *
 * @file preconditioner.cpp
 * @class Solution
 * @author Kelly Black <kjblack@gmail.com>
 * @version 0.2
 * @copyright BSD 2-Clause License
 *
 * @section LICENSE
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
 * @section DESCRIPTION
 *
 * Class to keep track of the preconditioner for the linearized
 * operator associated with a PDE
 *
 * This is the source file for the Preconditioner class. It includes
 * the basic operations to keep track of the preconditioner and for
 * solving the systems associated with the preconditioner.
 *
 *
 * @brief Basic operations associated with system associated with the
 * preconditioner.
 *
 * ********************************************************************************* */



#include "preconditioner.h"
#include "solution.h"
#include "../util.h"

#include <cmath>


/** ************************************************************************
 * Base constructor  for the Preconditioner class. 
 *
 *
 * @param size The number of grid points used in the approximation.
 * ************************************************************************ */
Preconditioner::Preconditioner(int number)
{
	setN(number);
	// allocate the vector with the lower diagonal matrix entries for
	// the Cholesky decomposition of the second order finite
	// difference operator for the same equation.
	cholesky = ArrayUtils<double>::twotensor(number+1,2);

	// Allocate the vector required to keep the intermediate results
	// of the solver when doing the backwards and forward solve from
	// the Cholesky decomposition.
	intermediate = ArrayUtils<double>::onetensor(number+1);

	// Define the values for the Cholesky decomposition of the finite
	// difference operator. This is the Cholesky decomposition of the
	// second order finite difference approximation of a Helmholtz
	// operator, d^2/dx^2 u + m*u.
	int lupe;
	double m = 0.0;
	double r = 2.0 + m;
	double xnum = (double)number;
	double tmp;
	for(lupe=0;lupe<=number;++lupe)
		{
			//cholesky[lupe][0] = sqrt(r);
			//r = (r*(2+m)-1)/r;
			tmp = sin(M_PI*((double)lupe)/xnum);
			tmp *= tmp;
			cholesky[lupe][0] = -(3.0*tmp*tmp)/((xnum*xnum-1.0)*tmp+3.0);

		}

	for(lupe=1;lupe<=number;++lupe)
		{
			cholesky[lupe][1] = -1.0/cholesky[lupe-1][0];
		}

}


/** ************************************************************************
 *	Copy constructor  for the Preconditioner class. 
 *
 *	@param oldCopy The Preconditioner class member to make a copy of.
 * ************************************************************************ */
Preconditioner::Preconditioner(const Preconditioner& oldCopy)
{
	setN(oldCopy.getN());
	// Allocate the vector for the preconditioner, and copy it over.
	cholesky = ArrayUtils<double>::twotensor(getN()+1,2);
	intermediate = ArrayUtils<double>::onetensor(getN()+1);

	int lupe;
	for(lupe=getN();lupe>=0;--lupe)
		{
			cholesky[lupe][0] = oldCopy.getValue(lupe,0);
			cholesky[lupe][1] = oldCopy.getValue(lupe,1);
		}
}

/** ************************************************************************
 *	Destructor for the Preconditioner class. 
 *  ************************************************************************ */
Preconditioner::~Preconditioner()
{
	ArrayUtils<double>::deltwotensor(cholesky);
	ArrayUtils<double>::delonetensor(intermediate);
}

/** ************************************************************************
 * The method to solve the system of equations associated with the
 * preconditioner.
 * 
 * Returns the value Solution class that is the solution to the
 * preconditioned system.  For the moment we just assume that the
 * preconditioner is the identity matrix. (This needs to be changed!)
 *
 * @param vector The Solution or right hand side of the system.
 * @return A Solution class member that is the solution to the
 *         preconditioned system.
 * ************************************************************************ */
Solution Preconditioner::solve(const Solution &current)
{
	Solution multiplied(current);
	int row;
	int col;
	int lupe;
	int N = current.getN();

	// Apply the Dirichlet boundary conditions on the top and bottom rows.

	for(row=1;row<N;++row)
		for(col=1;col<N;++col)
			multiplied(row,col) = multiplied(row,col)*cholesky[row][0];

	for(col=0;col<=N;++col)
		{
			multiplied(0,col) = current.getEntry(0,col);
			multiplied(N,col) = current.getEntry(N,col);
		}

	for(row=0;row<=N;++row)
		{
			multiplied(row,0) = current.getEntry(row,0);
			multiplied(row,N) = current.getEntry(row,N);
		}


	/*
	// Go through every y column and apply the preconditioner to that
	// column.
	for(row = N-1;row>0;--row)
		{

			// Perform the forward solve to invert the first part of the
			// Cholesky decomposition.
			intermediate[0] = current.getEntry(row,0)/cholesky[0][0];
			for(lupe=1;lupe<=N;++lupe)
				intermediate[lupe] = 
					(current.getEntry(row,lupe)-cholesky[lupe][1]*intermediate[lupe-1])
					/cholesky[lupe][0];

			// Perform the backwards solve for the Cholesky decomposition.
			multiplied(row,N) = intermediate[N]/cholesky[N][0];
			for(lupe=N-1;lupe>=0;--lupe)
				multiplied(row,lupe) = (intermediate[lupe]-multiplied(row,lupe+1)*cholesky[lupe+1][1])
					/cholesky[lupe][0];

			// The previous solves wiped out the boundary conditions. Restore
			// the left and right boundaru condition before sending the result
			// back.
			multiplied(row,0) = current.getEntry(row,0);
			multiplied(row,N) = current.getEntry(row,N);

		}

	// Apply the Dirichlet boundary conditions on the top and bottom rows.
	for(col=0;col<=N;++col)
		{
			multiplied(0,col) = current.getEntry(0,col);
			multiplied(N,col) = current.getEntry(N,col);
		}


	for(col = N-1;col>0;--col)
		{

			// Perform the forward solve to invert the first part of the
			// Cholesky decomposition.
			intermediate[0] = current.getEntry(0,col)/cholesky[0][0];
			for(lupe=1;lupe<=N;++lupe)
				intermediate[lupe] = 
					(current.getEntry(lupe,col)-cholesky[lupe][1]*intermediate[lupe-1])
					/cholesky[lupe][0];

			// Perform the backwards solve for the Cholesky decomposition.
			multiplied(N,col) = intermediate[N]/cholesky[N][0];
			for(lupe=N-1;lupe>=0;--lupe)
				multiplied(lupe,col) = (intermediate[lupe]-multiplied(lupe+1,col)*cholesky[lupe+1][1])
					/cholesky[lupe][0];

			// The previous solves wiped out the boundary conditions. Restore
			// the left and right boundaru condition before sending the result
			// back.
			multiplied(0,col) = current.getEntry(0,col);
			multiplied(N,col) = current.getEntry(N,col);

		}

	// Apply the Dirichlet boundary conditions on the top and bottom rows.
	for(row=0;row<=N;++row)
		{
			multiplied(row,0) = current.getEntry(row,0);
			multiplied(row,N) = current.getEntry(row,N);
		}

	*/



	return(multiplied);
}


