
#include <cmath>

#include "matrix.h"
#include "solution.h"
#include "util.h"

/** *********************************************************************************
 * @file matrix.cpp
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
 * This is the file with the code for various methods within the
 * Matrix class. 
 *
 *
 * @brief Code file for the basic operations associated with a matrix.
 *
 * ********************************************************************************* */


/** ************************************************************************
 * Base constructor  for the Matrix class. 
 *
 *
 *	@param numberRows The number of rows to use in the matrix.
 *  @param numberCols The number of columns to use in the matrix.
 * ************************************************************************ */
Matrix::Matrix(int numberRows,int numberCols) 
{
	setNumberRows(numberRows);
	setNumberColumns(numberCols);

	// allocate the space for the matrix.
	// Note that the twotensor routine initializes all the entries
	// to zero.
	matrix = ArrayUtils<double>::twotensor(numberRows+1,numberCols+1); 


}


/** ************************************************************************
 *		Copy constructor  for the Matrix class. 
 *
 *	@param oldCopy The Matrix class member to make a copy of.
 * ************************************************************************ */
Matrix::Matrix(const Matrix& oldCopy)
{
	int lupe;
	int innerLupe;

	setNumberRows(oldCopy.getNumberRows());
	setNumberColumns(oldCopy.getNumberColumns());
	matrix = ArrayUtils<double>::twotensor(getNumberRows()+1,getNumberColumns()+1);
	for(lupe=0;lupe<=getNumberRows();++lupe)
		{
			for(innerLupe=0;innerLupe<=getNumberColumns();++innerLupe)
				{
					setEntry(oldCopy.getEntry(lupe,innerLupe),lupe,innerLupe);
				}
		}	
}

/** ************************************************************************
 *	Destructor for the Solution class. 
 *  ************************************************************************ */
Matrix::~Matrix()
{
	ArrayUtils<double>::deltwotensor(matrix);
}



/** ************************************************************************
 * The parenthesis operator for the Matrix class.
 * 
 * Returns the value of the matrix for the indicated row and column.
 *
 * @param row    The row number to use.
 * @param column The column number to use.
 * @return a double precision value, matrix[row][column].
 * ************************************************************************ */
double& Matrix::operator()(int row,int column)
{
	return(matrix[row][column]);
}


/** ************************************************************************
 * The matrix/vector  multiplication operator for the Matrix class.
 * 
 * Returns a new Solution object which is the matrix/vector product of an
 * object from the matrix class and an object from the Solution class.
 *
 * @param vector  The Solution object to multiply by this matrix.
 * @return The result of the operation, an object from the Solution class.
 * ************************************************************************ */
Solution Matrix::operator*(class Solution vector)
{
	double tmp;
	int lupe;
	int innerLupe;
	Solution result;

	for(lupe=0;lupe<=getNumberRows();++lupe)
		{
			tmp = (*this)(lupe,0)*vector(0);
			for(innerLupe=1;innerLupe<=getNumberColumns();++innerLupe)
				tmp += (*this)(lupe,innerLupe)*vector(innerLupe);
			result.setEntry(tmp,lupe);
		}
	return(result);
}



/** ************************************************************************
 * The equals operator for the Matrix class.
 * 
 * Returns a new Matrix object that is identical to the Matrix
 * object passed to it.
 *
 * @param original The Matrix argument to copy
 * @return An object from the Matrix class.
 * ************************************************************************ */
Matrix Matrix::operator=(const Matrix& original)
{
	int lupe;
	int innerLupe;
	for(lupe=0;lupe<=getNumberRows();++lupe)
		{
			for(innerLupe=0;innerLupe<=getNumberColumns();++innerLupe)
				{
					this->setEntry(original.getEntry(lupe,innerLupe),lupe,innerLupe);
				}
		}	
}

