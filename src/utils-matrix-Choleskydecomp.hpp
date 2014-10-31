


#ifndef UTILS_MATRIX_CHOLESKYDECOMP_HPP
#define UTILS_MATRIX_CHOLESKYDECOMP_HPP



#include <iostream>
#include <string>
#include <sstream>

#include <cstdlib>

#include "utils.hpp"
#include "utils-math.hpp"
#include "utils-vector.hpp"
#include "utils-matrix.hpp"



namespace utils {


  /* #########################################################################
     Class to store the Cholesky decomposition of a given matrix.
     #########################################################################
  */


  template <typename T>
  class Choleskydecomp {

  private:
    int n;
    Matrix<T> el;

  public:

    Matrix<T> L(void){
      return el;
    }



    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Construct and store L in L * L^T = A.
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    Choleskydecomp(Matrix<T> & A, int & status)
      : n(A.nrows()), el(A) 
    {
      int i,j,k;
      T sum;
      
      status=0;
      if (el.ncols() != el.nrows()){
	std::cout << "Cannot construct Cholesky decomposition of non-square matrix. Exiting." << std::endl;
	exit(EXIT_FAILURE);
      }
      
      for (i=0; i!=n; ++i)
	for (j=0; j!=n; ++j)
	  el.elem(i,j)=0.0;
      

      
      for (i=0; i!=n; ++i){

	for (j=i; j!=n; ++j){
	  //sum=el.elem(i,j);
	  sum=A.elem(i,j);
	  for (k=i-1; k>=0; --k)
	    sum -= el.elem(i,k) * el.elem(j,k);

	  if (i==j){
	    // *******************************************************
	    if (sum <= 0){
	      std::cout << "ERROR: Cholesky decomposition failed (attempted square root of negative value)." << std::endl;
	      status=1;
	      return;
	    }
	    el.elem(i,i) = sqrt(sum);
	    // *******************************************************
	  }
	  else
	    el.elem(j,i) = sum/el.elem(i,i);
	}
      }
      
      for (i=0; i!=n; ++i) for (j=0; j!=i; ++j) el.elem(j,i)=0;


    }






    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Solve A * x = b.
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    void solve(Vector<T> & b, Vector<T> & x){
      int i,k;
      T sum;

      x.resize(b.size());

      if (b.size() != n || x.size() != n){
	std::cout << "Bad sizes of vectors in Cholesky solver. Exiting." << std::endl;
	exit(EXIT_FAILURE);
      }

      // Solve L * y = b, storing y in x.
      for (i=0; i!=n; ++i){
	sum = b[i];
	for (k=i-1; k>=0; --k) sum -= el.elem(i,k) * x[k];
	x[i] = sum/el.elem(i,i);
      }

      // Solve L^T * x = y.
      for (i=n-1; i>=0; --i){
	sum=x[i];
	for (k=i+1; k<n; ++k) sum -= el.elem(k,i)*x[k];
	x[i] = sum/el.elem(i,i);
      }
    }


    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Calculate L * y = b.
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    void elmult(Vector<T> & y, Vector<T> & b){
      int i,j;
      if (b.size() != n || y.size() != n){
	std::cout << "Bad sizes of vectors in Cholesky solver. Exiting." << std::endl;
	exit(EXIT_FAILURE);
      }
      for (i=0; i!=n; ++i){
	b[i]=0;
	for (j=0; j!=i; ++j) b[i] += el.elem(i,j) * y[j];
      }
    }


    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Solve L * y = b.
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    void elsolve(Vector<T> & b, Vector<T> & y){
      int i,j;
      T sum=0;
      
      if (b.size() != n || y.size() != n){
	std::cout << "Bad sizes of vectors in Cholesky solver. Exiting." << std::endl;
	exit(EXIT_FAILURE);
      }
      for (i=0; i!=n; ++i){
	sum=b[i];
	for (j=0; j!=i; ++j) sum -= el.elem(i,j)*y[j];
	  y[i] = sum/el.elem(i,i);
      }
    }
    

    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Calculate the inverse.
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
    void inverse(Matrix<T> & A_inv){
      int i,j,k;
      T sum;
      
      A_inv.resize(n,n);
      for (i=0; i!=n; ++i)
	for (j=0; j!=i; ++j) {
	  sum = (i==j ? 1 : 0);
	  for (k=i-1; k!=j; --k)
	    sum -= el.elem(i,k)*A_inv.elem(j,k);
	  A_inv.elem(j,i) = sum / el.elem(i,i);
	}
      for (i=n-1; i>=0; --i)
	for (j=0; j!=i; ++j){
	  sum = (i<j ? 0 : A_inv.elem(j,i) );
	  for (k=i+1; k!=n; ++k)
	    sum -= el.elem(k,i) * A_inv.elem(j,k);
	  A_inv.elem(i,j) = A_inv.elem(j,i) = sum / el.elem(i,i);
	}
    }


    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
    // Calculate the log. of the determinant.
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    T logdet(){
      T sum=0;
      
      for (int i=0; i!=n; ++i) sum += log(el.elem(i,i));
      return 2*sum;
    }
    
  } ;

}





#endif

