#ifndef BASEVARIABLE_H
#define BASEVARIABLE_H


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


using namespace std;


#define NUMBER 64



class BaseVariable
{

public:

    // Define the constructor and the routine used to initialize the
    // static variables.
    BaseVariable();
    void initializeStaticVariables();

    // Define the accessors for the approximation.
    // Get one value within the approximation vector
    inline double GetApproximation(int which)
    {
        return(BaseVariable::approx[which]);
    }


    // Set one value within the approximation vector
    inline void SetApproximation(double val,int which)
    {
        BaseVariable::approx[which] = val;
    }

    // Define the accessors for the value of x.
    inline double GetX(int which)
    {
        return(BaseVariable::x[which]);
    }


    // Define the accessors for the value of N
    inline int GetN()
    {
        return(BaseVariable::N);
    }

    // Accessors for the derivative matrices
    static inline double getDerivMat(int row,int col)
    {
        return(d[row][col]);
    }

    static inline double getSecondDerivMat(int row,int col)
    {
        return(d2[row][col]);
    }


    // Define the resolution of the approximation. Also define the
    // first derivative matrix (d) and the second derivative matrix
    // (d2). The grid points are given by x.
    static const int N = NUMBER;
    static double d[NUMBER+1][NUMBER+1],d2[NUMBER+1][NUMBER+1];
    static double x[NUMBER+1];

    // The matrix lhs is used to define the left hand side of an
    // implicitic scheme. (It is also used for the Jacobian used on
    // the left hand side to solve non-linear equations.)
    static double lhs[NUMBER+1][NUMBER+1];

    // Define the routines that initialize the first and second
    // derivative matrices.
    static void cheby1(double deriv[NUMBER+1][NUMBER+1],double x[],int num);
    static void cheby2(double deriv[NUMBER+1][NUMBER+1],double x[],int num);

    // Define the routines to perform an LU decomposition and to solve
    // the resulting system.
    static int  lu_decomp(double **a,int *place,int num);
    static void solve_lu(double **a,double *x,double *b,int *place,int num);

protected:

    // Define the approximation that is being tracked.
    double approx[NUMBER+1];

private:

};


#endif // BASEVARIABLE_H
