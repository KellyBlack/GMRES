
/** *********************************************************************************
 *
 * @file preconditioner.cpp
 * @class Solution
 * @author Kelly Black <kjblack@gmail.com>
 * @version 0.1
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
#include "util.h"

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
	// allocate the vector with the reciprocol of the diagonal entries
	// of the operator.
	vector = ArrayUtils<double>::onetensor(number+1);

	// Define the values for the vector....
	int lupe;
	vector[0] = 1.0;
	double xnum = (double)number;
  double dxnum = 1.0/((double)number);
	for(lupe=1;lupe<number;++lupe)
		{
			double tmp = sin(M_PI*((double)lupe)*dxnum);
			tmp *= tmp;
			vector[lupe] = (3.0*tmp*tmp)/(-((xnum*xnum-1.0)*tmp+3.0));
		}
	vector[number] = 1.0;
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
	vector = ArrayUtils<double>::onetensor(getN()+1);
	int lupe;
	for(lupe=getN();lupe>=0;--lupe)
		vector[lupe] = oldCopy.getValue(lupe);
}

/** ************************************************************************
 *	Destructor for the Preconditioner class. 
 *  ************************************************************************ */
Preconditioner::~Preconditioner()
{
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
	/*
	int lupe;
	for(lupe=getN();lupe>=0;--lupe)
		multiplied(lupe) *= vector[lupe];
	*/
	return(multiplied);
}


