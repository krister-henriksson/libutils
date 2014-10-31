




#ifndef UTILS_MATRIX_LUDECOMP_HPP
#define UTILS_MATRIX_LUDECOMP_HPP



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
     Class to store the LU decomposition of a given matrix.
     #########################################################################
  */


  template <typename T>
  class LUdecomp {

  private:
    int n;
    Matrix<T> LU;
    Vector<int> indx;
    T d;
    Matrix<T> & A_ref; // Reference which is bound to orginal matrix.
    

  public:
    LUdecomp();
    ~LUdecomp();

    LUdecomp(Matrix<T> & A);

    int solve(Vector<T> & b, Vector<T> & x) const; // Solve for single RHS.
    int solve(Matrix<T> & B, Matrix<T> & X) const; // Solve for multiple RHSs.
    int inverse(Matrix<T> & A_inv) const;
    T det() const;
    void improve(Vector<T> & b, Vector<T> & x) const; // Improve solution for single RHS.
    void improve(Matrix<T> & B, Matrix<T> & X) const; // Improve solution for multiple RHSs.
  } ;

}




// Constructor:
template <typename T>
utils::LUdecomp<T>::LUdecomp()
{ }


// Destructor:
template <typename T>
utils::LUdecomp<T>::~LUdecomp()
{ }



// Constructor:
template <typename T>
utils::LUdecomp<T>::LUdecomp(utils::Matrix<T> & A)
  : n(A.nrows()), LU(A), indx(n,0), A_ref(A)
{
  const T TINY = utils::eps_d();
  T big, temp;
  utils::Vector<T> vv(n, 0);
  int i, imax, j, k;

  if (A.nrows() != A.ncols()){
    std::cout << "Cannot construct LU decomposition of non-square matrix. Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }

  //  indx.resize(n);
  d=1.0;

  for (i=0; i!=n; ++i){
    big=0.0;
    for (j=0; j!=n; ++j)
      if ((temp=utils::absval(LU.elem(i,j))) > big) big=temp;
    if (utils::absval(big) < TINY){
      std::cout << "Singular matrix in LU decomposition routine. Exiting." << std::endl;
      exit(EXIT_FAILURE);
    }
    vv[i] = 1.0/big;
  }

  for (k=0; k!=n; ++k){
    big=0.0;
    imax=k;
    for (i=k; i!=n; ++i){
      temp = vv[i] * utils::absval(LU.elem(i,k));
      if (temp>big){
	big=temp;
	imax=i;
      }
    }
    if (k != imax){
      for (j=0; j!=n; ++j){
	temp = LU.elem(imax,j);
	LU.elem(imax,j) = LU.elem(k,j);
	LU.elem(k,j) = temp;
      }
      d = -d;
      vv[imax]=vv[k];
    }
    indx[k]=imax;
    if (utils::absval(LU.elem(k,k))<TINY) LU.elem(k,k)=TINY;
    for (i=k+1; i!=n; ++i){
      LU.elem(i,k) /= LU.elem(k,k);
      temp=LU.elem(i,k);
      for (j=k+1; j!=n; ++j)
	LU.elem(i,j) -= temp * LU.elem(k,j);
    }
  }
}






// Solve for single RHS.
// Solve A x = b.
//
template <typename T>
int utils::LUdecomp<T>::solve(utils::Vector<T> & b,
			      utils::Vector<T> & x) const {
  int i,ii=0,ip,j;
  T sum;

  if (b.size() != n || x.size() != n){
    std::cout << "Vector sizes do not fit the size of the LU decomposition. Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
  for (i=0; i!=n; ++i) x[i]=b[i];
  for (i=0; i!=n; ++i){
    ip = indx[i];
    sum=x[ip];
    x[ip]=x[i];
    if (ii != 0)
      for (j=ii-1; j!=i; ++j) sum -= LU.elem(i,j)*x[j];
    else if (! utils::fp_are_equal(sum,0.0))
      ii=i+1;
    x[i]=sum;
  }
  for (i=n-1; i>=0; --i){
    sum=x[i];
    for (j=i+1; j!=n; ++j) sum -= LU.elem(i,j)*x[j];
    x[i]=sum/LU.elem(i,i);
  }
  return 0;
}



// Solve for multiple RHSs.
// Solve A X = B, X returns A^(-1) B.
//
template <typename T>
int utils::LUdecomp<T>::solve(utils::Matrix<T> & B,
			      utils::Matrix<T> & X) const {
  int i,j,m=B.ncols();
  if (B.nrows() != n || X.nrows() != n || B.ncols() != X.ncols()){
    std::cout << "Matrix sizes do not fit the size of the LU decomposition. Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
  utils::Vector<T> xx(n,0);
  for (j=0; j!=m; ++j){
    for (i=0; i!=n; ++i) xx[i] = B.elem(i,j);
    this->solve(xx,xx);
    for (i=0; i!=n; ++i) X.elem(i,j) = xx[i];
  }
  return 0;
}




template <typename T>
int utils::LUdecomp<T>::inverse(utils::Matrix<T> & A_inv) const {
  int i,j;
  A_inv.resize(n,n);
  for (i=0; i!=n; ++i){
    for (j=0; j!=n; ++j) A_inv.elem(i,j)=0;
    A_inv.elem(i,i)=1;
  }
  return this->solve(A_inv, A_inv);
}




template <typename T>
T utils::LUdecomp<T>::det() const {
  T dd=d;

  for (int i=0; i!=n; ++i) dd *= LU.elem(i,i);
  return dd;
}



// Improve solution for single RHS.
/* Let the exact solution be A * x = b.
   The obtained solution is x+dx, which gives

   A * (d+dx) = b+db (*)

   Subtraction of equations gives A*dx=db (**).
   Solution of (*) for db is

   db = A * (x+dx) - b

   which gives after insertion into (**) that

   A * dx = A * (x+dx) - b = C

   where C is a known matrix when the solution x+dx
   has been obtained.
 */
template <typename T>
void utils::LUdecomp<T>::improve(utils::Vector<T> & b,
				 utils::Vector<T> & x) const {
  utils::Vector<T> bp(n,0);

  bp = A_ref * x - b;
  int status = solve(bp,bp);
  if (status==0) x = x - bp;
  return;
}

template <typename T>
void utils::LUdecomp<T>::improve(utils::Matrix<T> & b,
				 utils::Matrix<T> & x) const {
  utils::Matrix<T> bp(n,n,0);
  utils::Matrix<T> xp(n,n,0);

  bp = A_ref * x - b;
  int status = solve(bp,xp);
  if (status==0) x = x - xp;
  return;
}



#endif

