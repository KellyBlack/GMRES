#ifndef SOLUTIONCLASS
#define SOLUTIONCLASS


/** *********************************************************************************
 * @file solution.h
 * @class Solution
 * @author Kelly Black <kjblack@gmail.com>
 * @version 0.1
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
 *
 * @section DESCRIPTION
 *
 * Class to keep track of the approximation to a PDE.
 *
 * This is the definition (header) file for the Solution class. It
 * includes the definitions for the methods and the data used to keep
 * track of the approximation and to implement the basic algebraic
 * operations to perform on the approximation.
 *
 *
 * @brief header file for the basic operations associated with the
 * approximation to the PDE.
 *
 * ********************************************************************************* */

#include "poisson.h"
#include "util.h"

class Solution
{

public:
	Solution(int size=NUMBER);               //< Default constructor for the class
	Solution(const Solution& oldCopy);       //< Constructor for making a copy/duplicate
	~Solution();                             //< Destructor for the class

	// Now define the operators associated with the class.
	double& operator()(int row);                   //< The parenthesis operator for access to data elements
	Solution operator=(const Solution& vector);    //< Assignment operator for copying another Solution
	Solution operator=(const double& value);       //< Assignment operator for assigning a single value to all elements.
	Solution operator+(const Solution& vector);    //< Operator for adding two Solution objects
	Solution operator-(const Solution& vector);    //< Operator for subtracting two Solution objects.
	Solution operator*(const double& value);       //< Operator for scalar multiplication.
	double   operator*(const Solution& vector);    //< Operator for the dot product
	Solution operator*=(const double& value);      //< Operator for scalar multiplication in place.
	Solution operator-=(const Solution& vector);   //< Operator for subtracting another Solution object.
	Solution operator+=(const Solution& vector);   //< Operator for adding another Solution object.


	/** ************************************************************************
	 * The method to set the value of the entry in a row of the solution.
	 * 
	 * Sets the value of the indicated row to the value specified.
	 *
	 * @param value The scalar value (double) set the given row to.
	 * @param row The entry in the vector to change.
	 * @return N/A
	 * ************************************************************************ */
	void setEntry(double value,int row)
	{
		solution[row] = value;
	}


	/** Definition of the dot product of two approximation vectors. */
	static double dot (const Solution& v1,const Solution& v2);
	static double dot (Solution* v1,Solution* v2);

	/** Definition of the l2 norm of an approximation vector. */
	static double norm(const Solution& v1);
	double norm();

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
	   Method to get the number of elements used for the approximation.

	   @return The number of elements in the approximation.
	*/
	inline int getN() const
	{
		return(N);
	}

	/**
	   Method to get the value of the approximation at a certain grid point.

	   @param row The grid point where you want the height of the function.
	   @return The approximation at the given grid point.
	*/
	inline double getEntry(int row) const
	{
		return(solution[row]);
	}

protected:



private:

	// Define the size of the vector and the vector that will contain
	// the information.
	int N;                      //< The number of grid points.
	double *solution = NULL;    //< The vector that contains the approximation.

};


/**
	 Method to define scalar multiplcation for a solution vector.

	 The method to define how to multiply an approximation on the left using a scalar product.

	 @param value  The scalar to multiply the approximation by
	 @param vector The approximation that is being multiplied by the scalar.
	 @return The result of the scalar multiplication on the approximation.
 */
Solution operator*(const double& value,class Solution vector);

#endif
