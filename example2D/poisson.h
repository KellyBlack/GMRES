#ifndef POISSONCLASS
#define POISSONCLASS


/** *********************************************************************************
 * @file poisson.h
 * @class Poisson
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
 * Class to keep track of the operator and its  linearization  for a PDE.
 *
 * This is the definition (header) file for the Poisson class. It
 * includes the definitions for the methods and the data used to keep
 * track of the operator that defines a given PDE. It is assumed that
 * the operator is in the form Lu=0, and this class keeps track of L
 * and its linearization that is used for a GMRES procedure.
 *
 *
 * @brief header file for the basic operations associated with the
 * operator associated with a PDE.
 *
 * ********************************************************************************* */

#define NUMBER 64

class Solution;

class Poisson
{

public:
	Poisson(int number=NUMBER);        //< Default constructor for the Poisson Class.
	Poisson(const Poisson& oldCopy);   //< Copy constructor for the Poisson Class.
	~Poisson();                        //< Destructor for the Poisson Class.

	// Basic algebraic operators associated with the linearization of the operator.
	Solution operator*(class Solution vector);  //< The linearized operator acting on a given Solution.

	
	/**
		 Method to get the value of the x coordinate for a given row number.

		 @param row The row number you want to access.
		 @return The value of x at the given row number.
	 */
	double getX(int row) const
	{
		return(x[row]);
	}

	/**
		 Method to get the number of elements that are in the approximation.

		 @return The number of elements in the grid.
	 */
	int  getN() const
	{
		return(N);
	}

	/**
		 Method to get one of the elements from the first derivative matrix.

		 @param row The row number in the matrix
		 @param col The column number in the matrix.
		 @return The value within the matrix at the given row and column.
	 */
	inline double getD1(int row,int col) const
	{
		return(d1[row][col]);
	}


	/**
		 Method to get one of the elements from the second derivative matrix.

		 @param row The row number in the matrix
		 @param col The column number in the matrix.
		 @return The value within the matrix at the given row and column.
	 */
	inline double getD2(int row,int col) const
	{
		return(d2[row][col]);
	}

protected:

	// Define the routines that initialize the first and second
	// derivative matrices.
	void cheby1(double **deriv,double *x,int num);   //< Method to define the first derivative matrix
	void cheby2(double **deriv,double *x,int num);   //< Method to define the second derivative matrix


private:

	// Define the resolution of the approximation. Also define the
	// first derivative matrix (d) and the second derivative matrix
	// (d2). The grid points are given by x.
	int N;        //< The number of grid points in the approximation.
	double **d1;  //< Pointer to the first derivative matrix.
	double **d2;  //< Pointer to the second derivative matrix.
	double *x;    //< Pointer to the set of x grid points


};


#endif
