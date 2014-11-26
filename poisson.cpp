
/* *********************************************************************************
 * @file poisson.cpp
 * @class Poisson
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
 * Class to keep track of the operator associated with a PDE and its linearization.
 *
 * This is the code file for the Poisson class. It includes the code
 * for the methods to keep track of the operator that defines a given
 * PDE. It is assumed that the operator is in the form Lu=0, and this
 * class keeps track of L and its linearization that is used for a
 * GMRES procedure.
 *
 *
 * @brief header file for the basic operations associated with the
 * operator associated with a PDE.
 *
 * ********************************************************************************* */



#include "poisson.h"
#include "solution.h"
#include "util.h"

#include <cmath>


/** ************************************************************************
 * Base constructor  for the Poisson class. 
 *
 *
 *	@param number of grid points to use in the discretization.
 * ************************************************************************ */
Poisson::Poisson(int number)
{
  N  = number;
	d1 = ArrayUtils<double>::twotensor(number+1,number+1);
	d2 = ArrayUtils<double>::twotensor(number+1,number+1);
	x  = ArrayUtils<double>::onetensor(number+1);
	cheby1(d1,x,number);
	cheby2(d2,x,number);
}


/** ************************************************************************
 *		Copy constructor  for the Poisson class. 
 *
 *	@param oldCopy The Poisson class member to make a copy of.
 * ************************************************************************ */
Poisson::Poisson(const Poisson& oldCopy)
{
  N  = oldCopy.getN();
	d1 = ArrayUtils<double>::twotensor(N+1,N+1);
	d2 = ArrayUtils<double>::twotensor(N+1,N+1);
	x  = ArrayUtils<double>::onetensor(N+1);
	int lupe;
	int innerLupe;

  for (lupe=0;lupe<=N;++lupe) // go through every row.
		{
			x[lupe] = oldCopy.getX(lupe);
			for (innerLupe=0;innerLupe<=N;++innerLupe) // fill in the values for this row.
				{
					d1[lupe][innerLupe] = oldCopy.getD1(lupe,innerLupe);
					d2[lupe][innerLupe] = oldCopy.getD2(lupe,innerLupe);
				}
		}
}

/** ************************************************************************
 *	Destructor for the Poisson class. 
 *  ************************************************************************ */
Poisson::~Poisson()
{
	ArrayUtils<double>::deltwotensor(d1);
	ArrayUtils<double>::deltwotensor(d2);
	ArrayUtils<double>::delonetensor(x);
}


/** ************************************************************************
 * The parenthesis operator for the Poisson class.
 * 
 * Returns the value of the coefficient for the linearized operator for the indicated row and column.
 *
 * @param row The row number to use.
 * @param column The column number to use.
 * @return a double precision value, L[row][column]
 * ************************************************************************ */
double& Poisson::operator()(int row,int column)
{
	return(d2[row][column]);
}

/** ************************************************************************
 * The matrix/vector  multiplication operator for the Poisson class.
 * 
 * Returns a new Solution object which is the matrix/vector product of
 * an object from the Poisson class and an object from the Solution
 * class. This is the linearized version of the operator.
 *
 * @param vector  The Solution object to multiply by this matrix.
 * @return The result of the operation, an object from the Solution class.
 * ************************************************************************ */
Solution Poisson::operator*(class Solution vector)
{
	double tmp;
	int lupe;
	int innerLupe;
	Solution result;

	// the first and last row just return the same values, so 
	// there is no need to define the results from that row.
	result.setEntry(vector(0),0);
	for(lupe=1;lupe<getN();++lupe)
		{
			tmp = (*this)(lupe,0)*vector(0);
			for(innerLupe=1;innerLupe<=getN();++innerLupe)
				tmp += (*this)(lupe,innerLupe)*vector(innerLupe);
			result.setEntry(tmp,lupe);
		}
	result.setEntry(vector(getN()),getN());
	return(result);
}


/** ************************************************************************
 * The method to initialize the Chebychev collocation first derivative matrix.
 * 
 * It assumes that the memory for the matrix has been allocated. The
 * values of the matrix are initialized. Note that this does not make
 * use of the current state of the values for the object. This is done
 * to make it as generic as possible, and this method can be called by
 * friends to initialize values of their own matrices.
 *
 * @param deriv  A pointer to the derivative matrix to initialize.
 * @param xVal   The x grid points to initialize.
 * @param num    The number of grid points to use.
 * @return N/A
 * ************************************************************************ */
void  Poisson::cheby1(double **deriv,double *xVal,int num)

/*
      **********************************************
       Subroutine to initialize the Chebychev first
       derivitive matrix.  It also initializes the
       XVAL vector.  Note that up to and including
       (num) subscripts are used.
      **********************************************
*/

{
  int i,j;
  double xnum,dxnum;
  double tmp;

  // helper variables/factors used in various formulas.
  xnum = ((double) num);
  dxnum = 1.0/xnum;

  // define the values of x.
  for (i=1;i<num;++i)
     xVal[i] = cos(M_PI* ((double) i) * dxnum);
  xVal[0] = 1.0;
  xVal[num] = -1.0;

  for (i=0;i<num;++i) // go through every row.
      {
          for (j=i;j<=num;++j) // fill in the values for this row.
              {
                  if ((i==0) && (j==0))
                      {
                          // this is the upper left entry in the matrix
                          deriv[0][0] = (2.0*xnum*xnum + 1.0)/6.0;
                          deriv[num][num] = -deriv[0][0];
                      }

                  else if (i==j)
                      {
                          // this is a diagonal entry in the matrix.
                          tmp = 1.0/sin(M_PI*((double)i)*dxnum);
                          deriv[i][i] = -xVal[i]*tmp*tmp*0.5;
                      }

                  else
                      {
                          // This is an off diagonal entry.
                          deriv[i][j] =
                              0.5/(sin(M_PI*((double)(i+j))*dxnum*0.5)*sin(M_PI*((double)(-i+j))*dxnum*0.5));
                          // Add a mult. factor for the top and bottom
                          // rows as well as the left and right
                          // columns.
                          if (i%num == 0)
                              deriv[i][j] *= 2.0;
                          if (j%num == 0)
                              deriv[i][j] *= 0.5;
                          if ((i+j)%2 == 1)
                              deriv[i][j] *= -1.0;

                          // The matrix is anti-symmetric so fill in
                          // the opposite side of the matrix.
                          deriv[num-i][num-j] = - deriv[i][j];
                      }
              }
      }

}



/** ************************************************************************
 * The method to initialize the Chebychev collocation second derivative matrix.
 * 
 * It assumes that the memory for the matrix has been allocated. The
 * values of the matrix are initialized. Note that this does not make
 * use of the current state of the values for the object. This is done
 * to make it as generic as possible, and this method can be called by
 * friends to initialize values of their own matrices.
 *
 * @param deriv  A pointer to the derivative matrix to initialize.
 * @param xVal   The x grid points to initialize.
 * @param num    The number of grid points to use.
 * @return N/A
 * ************************************************************************ */
void Poisson::cheby2(double **deriv,double *xVal,int num)


/*
     *****************************************************
       Subroutine to initialize the second derivitive
       Chebychev matrix.  It also initializes the X
       vector.  Note that up to and including the num
       subscript is accessed.
    ******************************************************

*/

{
  double xnum,dxnum;
  int i,j;
  double tmp;

  // helper variables/factors used in various formulas.
  xnum = (double) num;
  dxnum = 1.0/xnum;

  // Define the values of x
  for (i=1;i<num;++i)
     xVal[i] = cos(M_PI*((double) i)*dxnum);
  xVal[0] = 1.0;
  xVal[num] = -1.0;

  for (i=0;i<=num/2;++i) // go through each row.
      {
          for (j=0;j<=num;++j) // go through each column.
              {

                  if (((i==0) && (j==0)) ||
                      ((i==num) && (j==num)))
                      // fill in the top left and bottom right entries
                      // in the matrix.
                      deriv[i][j] = (xnum*xnum*xnum*xnum-1.0)/15.0;

                  else if (i==0)
                      {
                          // This is the top row.
                          tmp = sin(M_PI*((double)j)*dxnum*0.5);
                          tmp *= tmp;
                          deriv[i][j] = ((2.0*xnum*xnum+1.0)
                                         *tmp-3.0)/(3.0*tmp*tmp);
                          // The sign of the values are alternating
                          if (j%2 == 1)
                              deriv[i][j] *= -1.0;

                          if ((j==0) || (j==num))
                              // include the multipler for the left and right columns.
                              deriv[i][j] *= 0.5;
                      }

                  else if (i==num)
                      {
                          // This is the bottom row.
                          tmp = cos(M_PI*((double)j)*dxnum*0.5);
                          tmp *= tmp;
                          deriv[i][j] = ((2.0*xnum*xnum+1.0)
                                         *tmp-3.0)/(3*tmp*tmp);

                          // The sign of the values are alternating.
                          if ((num+j)%2 == 1)
                               deriv[i][j] *= -1.0;

                           // Include the multipler for the left and right columns.
                           if ((j==0) || (j==num))
                               deriv[i][j] *= 0.5;

                       }

                   else if (i==j)
                       {
                           // This is the diagonal entry.
                           tmp = sin(M_PI*((double)i)*dxnum);
                           tmp *= tmp;
                           deriv[i][j] = -((xnum*xnum-1.0)*tmp+3.0)
                               /(3.0*tmp*tmp);
                       }

                   else
                       {
                           // This is an off diagonal entry.
                           tmp = sin(M_PI*((double)i)*dxnum);
                           tmp *= sin(M_PI*((double)(i+j))*dxnum*0.5);
                           tmp *= sin(M_PI*((double)(-i+j))*dxnum*0.5);
                           tmp *= tmp;
                           deriv[i][j] = (cos(M_PI*((double)i)*dxnum)*
                                          cos(M_PI*((double)(i+j))*dxnum*0.5) *
                                          cos(M_PI*((double)(i-j))*dxnum*0.5) - 1.0)*0.5/tmp;

                           // The signs of the values are alternating
                           if ((i+j)%2 == 1)
                               deriv[i][j] *= -1.0;

                           // Include the multiplier for the left and right columns.
                           if ((j==0) || (j==num))
                               deriv[i][j] *= 0.5;
                       }

               }
        }


   // The matrix is symmetric so fill in the rest of the matrix.
   for (i=num/2+1;i<=num;++i)
       for (j=0;j<=num;++j)
           deriv[i][j] = deriv[num-i][num-j];

}
