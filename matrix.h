#ifndef MATRIXCLASS
#define MATRIXCLASS


/** *********************************************************************************
 * @file matrix.h
 * @class Matrix
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
 * @section DESCRIPTION
 *
 * Class to keep track of a matrix. In particular it is a vector of
 * column vectors. It is used by the GMRES routines to build the upper
 * Hessenburg matrix and decompose it using a QR algorithm.
 *
 * This is the definition (header) file for the Matrix class. It
 * includes the definitions for the methods, the accessors, and the
 * data used to keep track of a matrix.
 *
 *
 * @brief Header file for the basic operations associated with a matrix.
 *
 * ********************************************************************************* */

#define NUMBER 64

class Solution;

class Matrix
{

public:
	Matrix(int numberRows=NUMBER,int numberCols=NUMBER); //< Default constructor for the class.
	Matrix(const Matrix& oldCopy);                       //< Constructor used to make a copy of the class.
	~Matrix();                                           //< Destructor for the class.

	// Algebraic operators used for operations on members of this class.
	double& operator()(int row,int column);      //< Access to data elements using parentheses.
	Solution operator*(class Solution vector);   //< Definition of matrix * vector operations.
	Matrix operator=(const Matrix& original);    //< Definition of a copy via assignment.

	// Accessor methods for the class.

	/**
		 Method to set the number of rows to use for the matrix.

		 @param number The number of rows to allocate
		 @return N/A
	 */
	void setNumberRows(const int number)
	{
		numberRows = number;
	}

	/**
		 Method to get the number of rows to use for the matrix.

		 @return The number of rows used in the matrix.
	 */
	int getNumberRows() const
	{
		return(numberRows);
	}

	/**
		 Method to set the number of columns to use for the matrix.

		 @param number The number of columns to allocate
		 @return N/A
	 */
	void setNumberColumns(const int number)
	{
		numberColumns = number;
	}

	/**
		 Method to get the number of columns to use for the matrix.

		 @return The number of columns.
	 */
	int getNumberColumns() const
	{
		return(numberColumns);
	}

	/**
		 Method to set one element within the matrix.

		 @param value The number to set the element to.
		 @param row   The row number to set.
		 @param col   The column number to set.
		 @return N/A
	 */
	void setEntry(const double value,int row,int col)
	{
		matrix[row][col] = value;
	}

	/**
		 Method to get the element within the matrix.

		 @param row The row number to get the element from.
		 @param col The column number to get the element from.
		 @return The value of the entry in the indicated row and column.
	 */
	double getEntry(int row,int col) const
	{
		return(matrix[row][col]);
	}

protected:



private:

	// Define the resolution of the approximation. Also define the
	// first derivative matrix (d) and the second derivative matrix
	// (d2). The grid points are given by x.
	int numberRows;       //< The number of rows in the matrix.
	int numberColumns;    //< The number of columns in the matrix.
	double **matrix;      //< Pointer to the matrix itself.


};


#endif
