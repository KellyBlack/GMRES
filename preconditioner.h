#ifndef PRECONDITIONERCLASS
#define PRECONDITIONERCLASS


/** *********************************************************************************
 * @file preconditioner.h
 * @class Preconditioner
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
 * This is the definition (header) file for the Preconditioner class. It
 * includes the definitions for the methods and the data used to keep
 * track of the preconditioner defined for the linearized system.
 *
 *
 * @brief header file for the basic operations associated with the
 * preconditioner for the linearized PDE.
 *
 * ********************************************************************************* */

#define NUMBER 64

class Solution;

class Preconditioner
{

public:
	Preconditioner(int number=NUMBER);              //< Default constructor for the class
	Preconditioner(const Preconditioner& oldCopy);  //< Constructor for making a copy/duplicate
	~Preconditioner();                              //< Destructor for the class

	Solution solve(const Solution &vector);    //< Method to solve the
																						 //< system associated with
																						 //< the preconditioner.

	
	/**
		 Method to set the number of elements to use for the length of the approximation.

		 @param number The number of elements.
		 @return N/A
	 */
	void setN(int number) 
	{
		N = number;
	}

	/**
		 Method to get the number of elements that are used for the approximation.

		 @return The number of grid points used in the approximation.
	 */
	int getN() const
	{
		return(N);
	}

	/**
		 Method to get the value of the preconditioner's vector for a
		 given row.

		 @param row Row in the vector you want to access.
		 @return The value of the vector for the given row.
	 */
	double getValue(int row) const
	{
		return(vector[row]);
	}


protected:


private:

	int N;          //< The number of grid points associated with the approximation.
	double *vector; //< The vector that has the reciprocol of the diagonal entries of the operator.

};




#endif
