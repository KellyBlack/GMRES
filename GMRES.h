
/** *********************************************************************************
 * @file GMRES.h
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
 * This file includes the template functions necessary to implement a
 * restarted GMRES algorithm. This is based on the pseudo code given in the book: 
 * 
 * @BOOK{templates,
 *    AUTHOR = {R. Barrett and M. Berry and T. F. Chan and J. Demmel and J. Donato and 
 *              J. Dongarra and V. Eijkhout and R. Pozo and C. Romine and H. Van der Vorst },
 *    TITLE = {Templates for the Solution of Linear Systems: 
 *                Building Blocks for Iterative Methods, 2nd Edition},
 *    PUBLISHER = {SIAM},
 *    YEAR = {1994},
 *    ADDRESS = {Philadelphia, PA} }
 * 
 * Also, some changes were implemented based on the matlab code by
 * John Burkardt
 * (http://people.sc.fsu.edu/~jburkardt/m_src/toms866/solvers/gmres_r.m
 * accessed 25 November 2015, which is based on a code by Time
 * Kelley.)
 *
 * The idea for using a template came from the IML++ code available at
 * http://math.nist.gov/iml++/ (accessed 26 November 2014). Also, the
 * method for calculating the entries for the Givens rotation matrices
 * came from the IML++ code as well.
 *
 * Additional sources that informed this work:
 * @BOOK{iterativeMethods,
 *     AUTHOR = {Tim Kelley},
 *     TITLE = {Iterative Methods for Linear and Nonlinear Equations},
 *     PUBLISHER = {SIAM}
 *     YEAR={2004},
 *     ADDRESS = {Philadelphia, PA},
 *     ISBN={0898713528} }
 *
 * @BOOK{sparseIterative,
 *       AUTHOR = {Yousef Saad},
 *       TITLE = {Iterative Methods for Sparse Linear Systems},
 *       EDITION = {Second Edition},
 *       PUBLISHER = {SIAM}
 *       YEAR={2003},
 *       ADDRESS = {Philadelphia, PA},
 *       ISBN={0898715342} } 
 *
 * @brief Template files for implementing a GMRES algorithm to solve a linear sytem.
 *
 * ********************************************************************************* */

#include "util.h"
#include <cmath>


/** ************************************************************************
 * Update the current approximation to the solution to the linear
 * system. This assumes that the update is created using a GMRES
 * routine, and the upper Hessenberg matrix has been transformed to an
 * upper diagonal matrix already. Note that it changes the values of
 * the values in the coefficients vector, s, which means that the s
 * vector cannot be reused after this without being re-initialized.
 *
 * @return N/A
 ************************************************************************ */
template <class Approximation, class Double >
void 
Update(Double **H,         //<! The upper diagonal matrix constructed in the GMRES routine.
	   Approximation *x,   //<! The current approximation to the linear system.
	   Double *s,          //<! The vector e_1 that has been multiplied by the Givens rotations.
	   Approximation v[],  //<! The orthogonal basis vectors for the Krylov subspace.
	   int dimension)      //<! The number of vectors in the basis for the Krylov subspace.
{

  // Solve for the coefficients, i.e. solve for c in
  // H*c=s, but we do it in place.
  int lupe;
  for (lupe = dimension; lupe >= 0; --lupe) 
	  {
		  s[lupe] = s[lupe]/H[lupe][lupe];
		  for (int innerLupe = lupe - 1; innerLupe >= 0; --innerLupe)
			  {
				  // Subtract off the parts from the upper diagonal of the
				  // matrix.
				  s[innerLupe] -=  s[lupe]*H[innerLupe][lupe];
			  }
	  }

  // Finally update the approximation.
  for (lupe = 0; lupe <= dimension; lupe++)
	  *x += v[lupe] * s[lupe];
}



/** ************************************************************************
 * Implementation of the restarted GMRES algorithm. Follows the
 * algorithm given in the book Templates for the Solution of Linear
 * Systems: Building Blocks for Iterative Methods, 2nd Edition.
 *
 * @return The number of iterations required. Returns zero if it did 
 *         not converge.
 ************************************************************************ */
template<class Operation,class Approximation,class Double>
int GMRES
(Operation* linearization, //!< Performs the linearization of the PDE on the approximation.
 Approximation* solution,  //!< The approximation to the linear system. (and initial estimate!)
 Approximation* rhs,       //!< the right hand side of the equation to solve.
 int krylovDimension,      //!< The number of vectors to generate in the Krylov subspace.
 int numberRestarts,       //!< Number of times to repeat the GMRES iterations.
 int approxLength,         //!< Total number of points being tracked in an approximation.
 Double  tolerance         //!< How small the residual should be to terminate the GMRES iterations.
 )
{

	// Allocate the space for the givens rotations, and the upper
	// Hessenburg matrix.
	Double **H      = ArrayUtils<Double>::twotensor(krylovDimension+1,krylovDimension);

	// The Givens rotations include the sine and cosine term. The
	// cosine term is in column zero, and the sine term is in column
	// one.
	Double **givens = ArrayUtils<Double>::twotensor(krylovDimension+1,2);

	// The vector s the right hand side for the system that the matrix
	// H satisfies in order to minimize the residual over the Krylov
	// subspace.
	Double *s = ArrayUtils<Double>::onetensor(krylovDimension+1);

	// Determine the residual and allocate the space for the Krylov
	// subspace.
	Approximation residual = (*rhs)-(*linearization)*(*solution);
	Approximation *V       = new Approximation[krylovDimension+1]; 
	Double rho             = residual.norm();
	Double normRHS         = rhs->norm();

	if(normRHS < 1.0E-5)
		normRHS = 1.0;

	// Go through the requisite number of restarts.
	int iteration = 1;
	while( (--numberRestarts >= 0) && (rho > tolerance*normRHS))
		{

			// The first vector in the Krylov subspace is the normalized
			// residual.
			V[0] = residual * (1.0/rho);

			// Need to zero out the s vector in case of restarts
			// initialize the s vector used to estimate the residual.
			for(int lupe=0;lupe<=krylovDimension;++lupe)
				s[lupe] = 0.0;
			s[0] = rho;

			// Go through and generate the pre-determined number of vectors
			// for the Krylov subspace.
			for( iteration=0;iteration<krylovDimension;++iteration)
				{
					// Get the next entry in the vectors that form the basis for
					// the Krylov subspace.
					V[iteration+1] = (*linearization)*V[iteration];

					// Perform the modified Gram-Schmidt to orthogonalize the
					// new vector.
					int row;
					for(row=0;row<=iteration;++row)
						{
							H[row][iteration] = Approximation::dot(V[iteration+1], V[row]);
							//subtract H[row][iteration]*V[row] from the current vector
							V[iteration+1].axpy(&V[row],-H[row][iteration]);
						}

					H[iteration+1][iteration] = V[iteration+1].norm();
					V[iteration+1] *= (1.0/H[iteration+1][iteration]);

					// Apply the Givens Rotations to insure that H is
					// an upper diagonal matrix. First apply previous
					// rotations to the current matrix.
					double tmp;
					for (row = 0; row < iteration; row++)
						{
							tmp = givens[row][0]*H[row][iteration] +
								givens[row][1]*H[row+1][iteration];
							H[row+1][iteration] = -givens[row][1]*H[row][iteration] 
								+ givens[row][0]*H[row+1][iteration];
							H[row][iteration]  = tmp;
						}

					// Figure out the next Givens rotation.
					if(H[iteration+1][iteration] == 0.0)
						{
							// It is already lower diagonal. Just leave it be....
							givens[iteration][0] = 1.0;
							givens[iteration][1] = 0.0;
						}
					else if (fabs(H[iteration+1][iteration]) > fabs(H[iteration][iteration]))
						{
							// The off diagonal entry has a larger
							// magnitude. Use the ratio of the
							// diagonal entry over the off diagonal.
							tmp = H[iteration][iteration]/H[iteration+1][iteration];
							givens[iteration][1] = 1.0/sqrt(1.0+tmp*tmp);
							givens[iteration][0] = tmp*givens[iteration][1];
						}
					else
						{
							// The off diagonal entry has a smaller
							// magnitude. Use the ratio of the off
							// diagonal entry to the diagonal entry.
							tmp = H[iteration+1][iteration]/H[iteration][iteration];
							givens[iteration][0] = 1.0/sqrt(1.0+tmp*tmp);
							givens[iteration][1] = tmp*givens[iteration][0];
						}

					// Apply the new Givens rotation on the
					// new entry in the uppper Hessenberg matrix.
					tmp = givens[iteration][0]*H[iteration][iteration] + 
						givens[iteration][1]*H[iteration+1][iteration];
					H[iteration+1][iteration] = -givens[iteration][1]*H[iteration][iteration] + 
						givens[iteration][0]*H[iteration+1][iteration];
					H[iteration][iteration] = tmp;

					// Finally apply the new Givens rotation on the s
					// vector
					tmp = givens[iteration][0]*s[iteration] + givens[iteration][1]*s[iteration+1];
					s[iteration+1] = -givens[iteration][1]*s[iteration] + givens[iteration][1]*s[iteration+1];
					s[iteration] = tmp;

					rho = fabs(s[iteration+1]);
					if(rho < tolerance*normRHS)
						{
							// We are close enough! Update the approximation.
							Update(H,solution,s,V,iteration);
							ArrayUtils<double>::deltwotensor(givens);
							ArrayUtils<double>::deltwotensor(H);
							ArrayUtils<double>::delonetensor(s);
							delete [] V;
							//tolerance = rho/normRHS;
							return(iteration);
						}

				}

			// We have exceeded the number of iterations. Update the
			// approximation and start over.
			Update(H,solution,s,V,iteration-1);
			residual = (*linearization)*(*solution) - (*rhs);
			rho = residual.norm();
		}


	ArrayUtils<double>::deltwotensor(givens);
	ArrayUtils<double>::deltwotensor(H);
	ArrayUtils<double>::delonetensor(s);
	delete [] V;
	//tolerance = rho/normRHS;

	if(rho < tolerance*normRHS)
		return(iteration);

	return(0);
}

