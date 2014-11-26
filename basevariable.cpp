
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


#include "basevariable.h"

#include <math.h>


// Reserve the static variables.
double BaseVariable::d[NUMBER+1][NUMBER+1],BaseVariable::d2[NUMBER+1][NUMBER+1];
double BaseVariable::x[NUMBER+1];



// Define the constructor
BaseVariable::BaseVariable()
{
}

// Initialize the static variables.
void BaseVariable::initializeStaticVariables()
{

    cheby1(d,x,N);
    cheby2(d2,x,N);

}



void  BaseVariable::cheby1(double deriv[N+1][N+1],double xVal[],int num)

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



void BaseVariable::cheby2(double deriv[N+1][N+1],double xVal[],int num)


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


int BaseVariable::lu_decomp(double **a,int *place,int num)
{

/*
     *******************************************************
         Subroutine to perform an L/U decomposition
         in place on the matrix A[][]
     *******************************************************
*/


  double work[N+1];
  register int i,j,k;
  int tmp,pivot;



  // Initialize the place vector. Go through row by row to see which
  // element has the largest entry in the row.
  for (i=0;i<num;++i)
      {
    place[i] = i;
    work[i] = fabs(a[i][0]);
    for (j=1;j<num;++j)
      if (work[i] < fabs(a[i][j]))
        work[i] = fabs(a[i][j]);
      }

  for (i=0;i<num;++i) // go through each row.
      {
          // First decide which column to pivot on for this row.
          // Use scaled partial pivoting.
          pivot = i;
          for (j=i+1;j<num;++j)
              if (fabs(a[place[j]][i])/work[place[j]] >
                  fabs(a[place[pivot]][i])/work[place[pivot]])
                  pivot = j;

          // Swap the pointers in the place vector.
          tmp = place[i];
          place[i] = place[pivot];
          place[pivot] = tmp;

          // bomb out if the diagonal entry is zero.
          // (I know this is bad.... need to fix this to check for "small" numbers.....)
          if (a[place[i]][i] == 0.0)
              {
                  return(0);
              }

          // Go through all the rows that follow and scale the row.
          // Also save the scaling used for the row/column.
          for (j=i+1;j<num;++j)
              if (fabs(a[place[j]][i]) > 1.0E-9)
                  {
                      a[place[j]][i] = a[place[j]][i]/a[place[i]][i];
                      for (k=i+1;k<num;++k)
                          a[place[j]][k] -= a[place[j]][i]*a[place[i]][k];
                  }

      }

  return(1);
}


void BaseVariable::solve_lu(double **a,double *result,double *b,int *place,int num)
/*
         Subroutine to solve the linear equation A*result=b for result.
         A is a matrix that has been decomposed into its L/U
         decomposition.  result and b are vectors where b is given
         and result is returned.  place is the index vector returned
         by the L/U decomposition.  size is the size of the
         matrix as it was declared and num is the number of
         values used in this application.
*/

{
 double c[N+1];
 int i,j;


 // Solve for the intermediate c vector.
 // First get the value for the top most values.
 c[place[0]] = b[place[0]];
 for (i=1;i<num;++i)
     {
         // Now go through and solve for all the values below using
         // the back solve.
         c[place[i]] = b[place[i]];
         for (j=0;j<i;++j)
             c[place[i]] -= a[place[i]][j]*c[place[j]];
     }

 // Solve for the final result value.
 // Solve for the lowest values first.
 result[num-1] = c[place[num-1]]/a[place[num-1]][num-1];
 for (i=num-2;i>=0;--i)
    {
        // go through and perform the forward solve for the previous
        // entries.
        result[i] = c[place[i]];
        for (j=i+1;j<num;++j)
            result[i] -= result[j]*a[place[i]][j];
        result[i] = result[i]/a[place[i]][i];
    }


}

