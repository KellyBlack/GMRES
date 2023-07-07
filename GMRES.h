
/** *********************************************************************************
 * @file GMRES.h
 * @author Kelly Black <kjblack@gmail.com>
 * @version 0.2
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
 * This file includes the template functions necessary to implement a
 * restarted GMRES algorithm. This is based on the pseudo code given
 * in the book Templates for the Solution of Linear Systems: Building
 * Blocks for Iterative Methods, 2nd Edition @cite templates .
 * 
 * Also, some changes were implemented based on the matlab code by
 * John Burkardt @cite BurkardtCode
 * http://people.sc.fsu.edu/~jburkardt/m_src/toms866/solvers/gmres_r.m
 *
 * The idea for using a template came from the IML++ code @cite imlCode
 * http://math.nist.gov/iml++/
 * Also, the method for calculating the entries for the Givens
 * rotation matrices came from the IML++ code as well.
 *
 * Additional sources that informed this work are Tim Kelley's book
 * Iterative Methods for Linear and Nonlinear Equations @cite iterativeMethods
 * Another book is is Yousef Saad's book Iterative Methods for 
 * Sparse Linear Systems  @cite sparseIterative
 *
 * @brief Template files for implementing a GMRES algorithm to solve a linear sytem.
 *
 * @mainpage Overview
 *
 * @latexonly
 * \input{overview}
 * @endlatexonly
 *
 * ********************************************************************************* */

#include "util.h"
#include <cmath>
#include <vector>

/** ************************************************************************
 * Update the current approximation to the solution to the linear
 * system. This assumes that the update is created using a GMRES
 * routine, and the upper Hessenberg matrix has been transformed to an
 * upper diagonal matrix already. Note that it changes the values of
 * the values in the coefficients vector, s, which means that the s
 * vector cannot be reused after this without being re-initialized.
 *
 ************************************************************************ */
template <class Approximation, class Double >
void Update
(Double **H,         //<! The upper diagonal matrix constructed in the GMRES routine.
 Approximation *x,   //<! The current approximation to the linear system.
 Double *s,          //<! The vector e_1 that has been multiplied by the Givens rotations.
 std::vector<Approximation> *v,  //<! The orthogonal basis vectors for the Krylov subspace.
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
	typename std::vector<Approximation>::iterator ptr = v->begin();
  for (lupe = 0; lupe <= dimension; lupe++)
	  *x += (*ptr++) * s[lupe];
}



/** ************************************************************************
 * Implementation of the restarted GMRES algorithm. Follows the
 * algorithm given in the book Templates for the Solution of Linear
 * Systems: Building Blocks for Iterative Methods, 2nd Edition.
 *
 * @return The number of iterations required. Returns zero if it did not converge.
 ************************************************************************ */
template<class Operation,class Approximation,class Preconditioner,class Double>
int GMRES
(Operation* linearization, //!< Performs the linearization of the PDE on the approximation.
 Approximation* solution,  //!< The approximation to the linear system. (and initial estimate!)
 Approximation* rhs,       //!< the right hand side of the equation to solve.
 Preconditioner* precond,  //!< The preconditioner used for the linear system.
 int krylovDimension,      //!< The number of vectors to generate in the Krylov subspace.
 int numberRestarts,       //!< Number of times to repeat the GMRES iterations.
 Double tolerance          //!< How small the residual should be to terminate the GMRES iterations.
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
	std::vector<Approximation> V(krylovDimension+1,
								 Approximation(solution->getN()));
	Approximation residual = precond->solve((*rhs)-(*linearization)*(*solution));
	Double rho             = residual.norm();
	Double normRHS         = rhs->norm();

	// variable for keeping track of how many restarts had to be used.
	int totalRestarts = 0;

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
					V[iteration+1] = precond->solve((*linearization)*V[iteration]);

					// Perform the modified Gram-Schmidt method to orthogonalize
					// the new vector.
					int row;
					typename std::vector<Approximation>::iterator ptr = V.begin();
					for(row=0;row<=iteration;++row)
						{
							H[row][iteration] = Approximation::dot(V[iteration+1], *ptr);
							//subtract H[row][iteration]*V[row] from the current vector
							V[iteration+1].axpy(&(*ptr++),-H[row][iteration]);
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
							Update(H,solution,s,&V,iteration);
							ArrayUtils<double>::deltwotensor(givens);
							ArrayUtils<double>::deltwotensor(H);
							ArrayUtils<double>::delonetensor(s);
							//delete [] V;
							//tolerance = rho/normRHS;
							return(iteration+totalRestarts*krylovDimension);
						}

				} // for(iteration)

			// We have exceeded the number of iterations. Update the
			// approximation and start over.
			totalRestarts += 1;
			Update(H,solution,s,&V,iteration-1);
			residual = precond->solve((*linearization)*(*solution) - (*rhs));
			rho = residual.norm();

		} // while(numberRestarts,rho)


	ArrayUtils<double>::deltwotensor(givens);
	ArrayUtils<double>::deltwotensor(H);
	ArrayUtils<double>::delonetensor(s);
	//delete [] V;
	//tolerance = rho/normRHS;

	if(rho < tolerance*normRHS)
		return(iteration+totalRestarts*krylovDimension);

	return(0);
}

