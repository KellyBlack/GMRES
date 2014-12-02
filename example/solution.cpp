

/** *********************************************************************************
 * @file solution.cpp
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
 *
 * @section DESCRIPTION
 *
 * Class to keep track of the approximation to a PDE.
 *
 * This is the source file for the Solution class. It includes the
 * basic operations to keep track of the approximation and to
 * implement the basic algebraic operations to perform on the
 * approximation.
 *
 *
 * @brief Basic operations associated with the approximation to the PDE.
 *
 * ********************************************************************************* */


//#include <iostream>
#include <cmath>



#include "solution.h"
#include "poisson.h"
#include "../util.h"

/** ************************************************************************
 * Base constructor  for the Solution class. 
 *
 *
 * @param size The length of the vector used in the approximation (optional)
 * ************************************************************************ */
Solution::Solution(int size)
{
	// Set the size of the vector, allocate the space, and zero out the
	// approximation.
	setN(size);
	solution = ArrayUtils<double>::onetensor(size+1);  // allocate the space. 
	// Note that the onetensor routine sets everything to zero so it
	// does not have to be initialized.
}

/** ************************************************************************
 *		Copy constructor  for the Solution class. 
 *
 * @overload
 * @param oldCopy The Solution class member to make a copy of.
 * ************************************************************************ */
Solution::Solution(const Solution& oldCopy)
{
	// Make a copy of the Solution that is passed to me.
	// Set the size of the vector, allocate the space, and then
	// copy the values over.
	int size = oldCopy.getN();
	setN(size);
	solution = ArrayUtils<double>::onetensor(size+1);
	for(;size>=0;--size)
		setEntry(oldCopy.getEntry(size),size);
}

/** ************************************************************************
 *	Destructor for the Solution class. 
 *  ************************************************************************ */
Solution::~Solution()
{
	// delete the approximation.
	if(solution)
		ArrayUtils<double>::delonetensor(solution);
	solution = NULL;
}

/** ************************************************************************
 * The parenthesis operator for the Solution class.
 * 
 * Returns the value of the approximation for the indicated row.
 *
 * @param row The row number to use.
 * @return a double precision value, solution[row]
 * ************************************************************************ */
double& Solution::operator()(int row)
{
	return(solution[row]);
}

/** ************************************************************************
 * The equals operator for the Solution class.
 * 
 * Returns a new Solution object that is identical to the Solution
 * object passed to it.
 *
 * @param vector The Solution argument to copy
 * @return An object from the Solution class.
 * ************************************************************************ */
Solution Solution::operator=(const Solution& vector)
{
	int lupe;

	if(this != &vector)
		{
			for(lupe=getN();lupe>=0;--lupe)
				{
					this->setEntry(vector.getEntry(lupe),lupe);
				}
		}
	return(*this);
}

/** ************************************************************************
 * The equals operator for the Solution class.
 * 
 * Returns a new Solution object in which every entry in the vector is
 * the same as the double precision number passed to it.
 *
 * @overload
 * @param value The value to copy into the vector.
 * @return An object from the Solution class.
 * ************************************************************************ */
Solution Solution::operator=(const double& value)
{
	int lupe;

	for(lupe=getN();lupe>=0;--lupe)
		{
			this->setEntry(value,lupe);
		}

	return(*this);
}

/** ************************************************************************
 * The summation operator for the Solution class.
 * 
 * Returns a new Solution object which is the sum of two Solution
 * objects.
 *
 * @param vector The Solution object to add to the current Solution object.
 * @return An object from the Solution class.
 * ************************************************************************ */
Solution Solution::operator+(const Solution& vector)
{
	Solution result(this->getN());
	int lupe;
	for(lupe=getN();lupe>=0;--lupe)
		{
			result.setEntry(this->getEntry(lupe)+vector.getEntry(lupe),lupe);
		}

	return(result);

}

/** ************************************************************************
 * The subtraction operator for the Solution class.
 * 
 * Returns a new Solution object which is the difference of two Solution
 * objects.
 *
 * @param vector The Solution object to subtract from the current Solution object.
 * @return An object from the Solution class.
 * ************************************************************************ */
Solution Solution::operator-(const Solution& vector)
{
	Solution result(this->getN());
	int lupe;
	for(lupe=getN();lupe>=0;--lupe)
		{
			result.setEntry(this->getEntry(lupe)-vector.getEntry(lupe),lupe);
		}

	return(result);

}

/** ************************************************************************
 * The scalar multiplication operator for the Solution class.
 * 
 * Returns a new Solution object which is the scalar product of an
 * object from the Solution class and a single double precision
 * number.
 *
 * @param value Scalar value to multiply every entry in the current vector.
 * @return An object from the Solution class.
 * ************************************************************************ */
Solution Solution::operator*(const double& value)
{
	Solution result(this->getN());
	int lupe;
	for(lupe=getN();lupe>=0;--lupe)
		{
			result.setEntry(this->getEntry(lupe)*value,lupe);
		}

	return(result);
}

/** ************************************************************************
 * The dot product operation for the Solution class.
 * 
 * Returns the dot product of this object and another object in the
 * Solution class.
 *
 * @overload
 * @param vector The other solution object to use for the dot product.
 * @return the dot product.
 * ************************************************************************ */
double Solution::operator*(const Solution& vector)
{
	double dot = getEntry(0)*vector.getEntry(0);
	int lupe;
	for(lupe=getN();lupe>0;--lupe)
		{
			dot += getEntry(lupe)*vector.getEntry(lupe);
		}
	return(dot);
}

/** ************************************************************************
 * The scalar multiplication operator for "*=" for the Solution class.
 * 
 * Returns a new Solution object which is the scalar product of an
 * object from the Solution class and a single double precision
 * number.
 *
 * @param value Scalar value to multiply every entry in the current vector.
 * @return An object from the Solution class.
 * ************************************************************************ */
Solution Solution::operator*=(const double& value)
{
	int lupe;
	for(lupe=getN();lupe>=0;--lupe)
		{
			this->setEntry(this->getEntry(lupe)*value,lupe);
		}

	return(*this);
}

/** ************************************************************************
 * The subtraction operator for "-=" for the Solution class.
 * 
 * Returns a new Solution object which is found by subtracting each
 * element from the passed Solution object from the current object.
 *
 * @param vector The Solution object to subtract from this object.
 * @return An object from the Solution class.
 * ************************************************************************ */
Solution Solution::operator-=(const Solution& vector)
{
	int lupe;
	for(lupe=getN();lupe>=0;--lupe)
		{
			this->setEntry(this->getEntry(lupe)-vector.getEntry(lupe),lupe);
		}

	return(*this);
}

/** ************************************************************************
 * The addition operator for "+=" for the Solution class.
 * 
 * Returns a new Solution object which is found by adding each
 * element from the passed Solution object to the current object.
 *
 * @param vector The Solution object to add to this object.
 * @return An object from the Solution class.
 * ************************************************************************ */
Solution Solution::operator+=(const Solution& vector)
{
	int lupe;
	for(lupe=getN();lupe>=0;--lupe)
		{
			this->setEntry(this->getEntry(lupe)+vector.getEntry(lupe),lupe);
		}

	return(*this);
}

/** ************************************************************************
 * The scalar multiplication operator for the Solution class.
 * 
 * Returns a new Solution object which is found by multiplying each
 * element from the passed Solution object by the scalar value.
 *
 * @relates Solution
 * @param value The scalar value (Double) to multiply.
 * @param vector The Solution object to multiply the scalar by.
 * @return An object from the Solution class.
 * ************************************************************************ */
Solution operator*(const double &value,class Solution vector)
{
	return(vector*value);
}


/** ************************************************************************
 * The method to find the dot product between two solutions.
 * 
 * Takes two solution vectors and computes their dot product. 
 *
 * @param value The first Solution object to use.
 * @param row The second Solution object to use.
 * @return The scalar dot product.
 * ************************************************************************ */
double Solution::dot(const Solution& v1,const Solution& v2)
{

	double dotProduct = v1.getEntry(0)*v2.getEntry(0);
	int lupe;
	//std::cout << "dot product" << std::endl;
	for(lupe=v1.getN();lupe>0;--lupe)
		{
			dotProduct += v1.getEntry(lupe)*v2.getEntry(lupe);
			//std::cout << dotProduct << "," << v1.getEntry(lupe) << "," << v2.getEntry(lupe) <<std::endl;
		}
	return(dotProduct);
}

/** ************************************************************************
 * The method to find the dot product between two solutions.
 * 
 * Takes two solution vectors and computes their dot product. 
 *
 * @overload
 * @param value The first Solution object to use.
 * @param row The second Solution object to use.
 * @return The scalar dot product.
 * ************************************************************************ */
double Solution::dot(Solution* v1,Solution* v2)
{

	double dotProduct = v1->getEntry(0)*v2->getEntry(0);
	int lupe;
	for(lupe=v1->getN();lupe>0;--lupe)
		{
			dotProduct += v1->getEntry(lupe)*v2->getEntry(lupe);
		}
	return(dotProduct);
}


/** ************************************************************************
 * The method to find the norm of a solution.
 * 
 * Takes a Solution object and computes the l2 norm of the solution.
 *
 * @param v1 The  Solution object to use.
 * @return The norm of the Solution object.
 * ************************************************************************ */
double Solution::norm(const Solution& v1)
{
	double norm = v1.getEntry(0)*v1.getEntry(0);
	int lupe;
	for(lupe=v1.getN();lupe>0;--lupe)
		{
			norm += v1.getEntry(lupe)*v1.getEntry(lupe);
		}
	return(sqrt(norm));
}

/** ************************************************************************
 * The method to find the norm of a solution.
 * 
 * Determines the l2 norm of the associated object.
 *
 * @overload
 * @return The norm of the Solution object.
 * ************************************************************************ */
double Solution::norm()
{
	double norm = getEntry(0)*getEntry(0);
	int lupe;
	for(lupe=getN();lupe>0;--lupe)
		{
			norm += getEntry(lupe)*getEntry(lupe);
		}
	return(sqrt(norm));
}


/** ************************************************************************
 * The method to perform an axpy procedure.
 * 
 * Subtracts multipler*vector from the current approximation.
 *
 * @param The vector to add
 * @param The scalar multiple to add
 * @return N/A
 * ************************************************************************ */
void Solution::axpy(Solution* vector,
					double multiplier)
{
	int lupe;
	for(lupe=getN();lupe>0;--lupe)
		{
			setEntry(getEntry(lupe)+multiplier*vector->getEntry(lupe),lupe);
		}
}
