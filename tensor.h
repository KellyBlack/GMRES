#ifndef TENSOR_HPP
#define TENSOR_HPP


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

#include <cstdlib>

template <class T>
class Tensor
{

public:

    T **twotensor(int n1,int n2);
    void deltwotensor(T **u);

    T *onetensor(int n1);
    void delonetensor(T *u);

protected:

private:

};


template <class T>
T **Tensor<T>::twotensor(int n1,int n2)
{
  T **u;
  register int s,i;

  u = new T*[n1];


  if(u==NULL)
      {
          //printf("Error - twotensor. Could not allocate memory.\n");
          exit(2);
      }

  u[0] = new T[n1*n2];
  if(u[0]==NULL)
      {
          //printf("Error - twotensor. Could not allocate memory for vector.\n");
          exit(2);
      }

  for(s=0;s<n1;++s)
      {
          u[s] = u[0] + s*n2;

          for(i=0;i<n2;++i)
              u[s][i] = 0.0;

      }

  return(u);

}

template <class T>
void Tensor<T>::deltwotensor(T **u)
{

    if(u==NULL)
        return;

    delete u[0];
    delete [] u;

}


template <class T>
T *Tensor<T>::onetensor(int n1)
{
  T *u;
  register int i;

  u = new T[n1];

  if (u==NULL)
      {
          //printf("Error - onetensor. Could not allocate memory.\n");
          exit(2);
      }

  for(i=0;i<n1;++i)
    u[i] = 0.0;

  return(u);

}

template <class T>
void Tensor<T>::delonetensor(T *u) {

    if(u==NULL)
        return;

  delete [] u;

}




#endif // TENSOR_HPP
