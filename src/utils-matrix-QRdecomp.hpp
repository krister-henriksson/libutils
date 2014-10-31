

//
//  WARNING!!! BUGGY CLASS AND METHODS!!!
//

/*

  DOES NOT WORK IN TEST CASES. IMPLEMENTATION FROM NUMERICAL RECIPES
  IS FAULTY ? TRANSCRIPTION FROM BOOK CHECKED AND SHOULD BE OK.

 */



#ifndef UTILS_MATRIX_QRDECOMP_HPP
#define UTILS_MATRIX_QRDECOMP_HPP



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
     Class to store the QR decomposition of a given matrix.
     #########################################################################
  */


  template <typename T>
  class QRdecomp {

  private:
    int n;
    Matrix<T> qt, r;
    bool sing;

  public:
    QRdecomp(Matrix<T> & A);

    int solve(Vector<T> & b, Vector<T> & x);
    void qtmult(Vector<T> & b, Vector<T> & x);
    int rsolve(Vector<T> & b, Vector<T> & x);

  } ;

}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Construct the QR decomposition.
// A = Q * R
// R is upper triangular
// Q^T * Q = 1
// 
// To solve A*x=b form Q^T*b and solve R*x = Q^T*b
// 
// Restricted to square matrices and no pivoting used in this implementation.
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template <typename T>
utils::QRdecomp<T>::QRdecomp(utils::Matrix<T> & A)
  : n(A.nrows()), qt(n,n,0), r(A), sing(false)
{
  int i,j,k;
  utils::Vector<T> c(n,0), d(n,0);
  T scale,sigma,sum,tau;

  for (k=0; k!=n-1; ++k){
    scale=T(0);
    for (i=k; i!=n; ++i) scale=utils::max(scale, utils::absval(r.elem(i,k)));
    if (utils::fp_are_equal(scale,T(0))){
      sing=true;
      c[k]=d[k]=T(0);
    }
    else {
      for (i=k; i!=n; ++i) r.elem(i,k) /= scale;
      sum=T(0);
      for (i=k; i!=n; ++i) sum += utils::square(r.elem(i,k));
      sigma = utils::sign_nr(sqrt(sum),r.elem(k,k));
      r.elem(k,k) += sigma;
      c[k]=sigma*r.elem(k,k);
      d[k]=-scale*sigma;
      for (j=k+1; j!=n; ++j){
	sum=0;
	for (i=k; i!=n; ++i) sum += r.elem(i,k)*r.elem(i,j);
	tau=sum/c[k];
	for (i=k; i!=n; ++i) r.elem(i,j) -= tau*r.elem(i,k);
      }
    }
  }
  d[n-1]=r.elem(n-1,n-1);
  if (utils::fp_are_equal(d[n-1],T(0))) sing=true;
  for (i=0; i!=n; ++i){
    for (j=0; j!=n; ++j) qt.elem(i,j)=T(0);
    qt.elem(i,i)=1;
  }
  for (k=0; k!=n-1; ++k){
    if (! utils::fp_are_equal(c[k],T(0))){
      for (j=0; j!=n; ++j){
	sum=T(0);
	for (i=k; i!=n; ++i)
	  sum += r.elem(i,k)*qt.elem(i,j);
	sum /= c[k];
	for (i=k; i!=n; ++i)
	  qt.elem(i,j) -= sum*r.elem(i,k);
      }
    }
  }
  for (i=0; i!=n; ++i){
    r.elem(i,i)=d[i];
    for (j=0; j!=i; ++j) r.elem(i,j)=0;
  }
}





// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Solve A*x=b.
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template <typename T>
int utils::QRdecomp<T>::solve(utils::Vector<T> & b,
			      utils::Vector<T> & x){
  if (x.size() != n) x.resize(n);
  qtmult(b,x);
  return rsolve(x,x);
}




// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Perform Q^T * b and store the result in x. Since Q is orthogonal,
// this is equivalent to solving Q*x=b for x.
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template <typename T>
void utils::QRdecomp<T>::qtmult(utils::Vector<T> & b,
				utils::Vector<T> & x){
  int i,j;
  T sum;

  if (x.size() != n) x.resize(n);
  for (i=0; i!=n; ++i){
    sum=T(0);
    for (j=0; j!=n; ++j) sum += qt.elem(i,j)*b[j];
    x[i]=sum;
  }
}





// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Solve R*x=b for x.
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template <typename T>
int utils::QRdecomp<T>::rsolve(utils::Vector<T> & b,
			       utils::Vector<T> & x){
  int i,j;
  T sum;

  if (x.size() != n) x.resize(n);  
  if (sing){
    std::cout << "QR decomposition is singular." << std::endl;
    return 1;
  }
  for (i=n-1; i>=0; --i){
    sum=b[i];
    for (j=i+1; j!=n; ++j) sum -= r.elem(i,j)*x[j];
    x[i]=sum/r.elem(i,i);
  }
  return 0;
}







#endif

