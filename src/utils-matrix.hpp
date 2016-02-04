



#ifndef UTILS_MATRIX_HPP
#define UTILS_MATRIX_HPP



#include <iostream>
#include <string>
#include <sstream>

#include <cstdlib>

#include "utils.hpp"
#include "utils-math.hpp"
#include "utils-vector.hpp"

using utils::abs;
using utils::Vector;
using utils::swap;


namespace utils {

  /* #########################################################################
     Declaration of mat<T> template class:
     #########################################################################
  */

  template <typename T>
  class Matrix {

  private:
    T *mmat;
    int mnrows;
    int mncols;
    
  public:
    Matrix();
    ~Matrix();

    Matrix(const int Nr, const int Nc);
    Matrix(const int Nr, const int Nc, const T & p);

    Matrix(const Matrix<T> & sv);

    // Operators:
    Matrix<T> & operator=(const Matrix<T> & sv);

    // Storage related functions:
    void resize(const int Nr, const int Nc);
    void reshape(const int Nr, const int Nc);
    // Linear size:
    int nrows() const;
    int ncols() const;
    // Total number of elements:
    int size() const;

    // Matrix element operations:
    inline       T & elem(const int & i, const int & j);
    inline const T & elem(const int & i, const int & j) const;

    // Matrix row and column operations:
    // Set a row:
    void row(const int & i, const Vector<T> & v);
    // Get a row:
    Vector<T> row(const int & i) const;
    // Set a column:
    void col(const int & i, const Vector<T> & v);
    // Get a column:
    Vector<T> col(const int & i) const;

    // Other matrix operations:

    // Set object to unity matrix:
    void unity();
    // Get transpose:
    Matrix<T> transpose();

    // Get rank (number of independent rows):
    int rank() const;

    // Invert matrix:
    /* Matrix inversion can be done by Gauss-Jordan elimination and
       back-substitution. Another method is to perform LU decomposition
       first and then calculate the inverse using the decomposition
       matrices. Both methods can be used to solve a system of linear
       equations, so that the inverse matrix is produced as a side
       effect.
    */
    int solve(Vector<T> & b,
	      Vector<T> & x,
	      Matrix<T> & A_inv) const;
    int solve(Matrix<T> & B,
	      Matrix<T> & X,
	      Matrix<T> & A_inv) const;
    int inverse(Matrix<T> & A_inv) const;
    T det() const;

   
  } ;


  // Handle a + b:
  template <typename T>
  Matrix<T> operator+(const Matrix<T> & a, const Matrix<T> & b);
  
  // Handle a - b:
  template <typename T>
  Matrix<T> operator-(const Matrix<T> & a, const Matrix<T> & b);

  // Handle a * b:
  template <typename T>
  Matrix<T> operator*(const Matrix<T> & a, const Matrix<T> & b);

  // Handle a * vec(b):
  template <typename T>
  Vector<T> operator*(const Matrix<T> & a, const Vector<T> & b);

  // Handle vec(a) * b:
  template <typename T>
  Vector<T> operator*(const Vector<T> & a, const Matrix<T> & b);

  // Handle a * b, b being S type:
  template <typename T, typename S>
  Matrix<T> operator*(const Matrix<T> & a, const S & b);

  // Handle a * b, a being S type:
  template <typename T, typename S>
  Matrix<T> operator*(const S & a, const Matrix<T> & b);

  // Handle a / b, b being S type:
  template <typename T, typename S>
  Matrix<T> operator/(const Matrix<T> & a, const S & b);


  // Handle printing of Matrix<T> objects:
  template <typename T, typename U>
  U & operator<<(U & os, const Matrix<T> & sv);

  
}




/* #########################################################################
   Definitions:
   #########################################################################
*/

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Default ctor:
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template <typename T>
utils::Matrix<T>::Matrix(){
  mmat=0;
  mnrows=0;
  mncols=0;
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Default dtor:
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template <typename T>
utils::Matrix<T>::~Matrix(){
  if (mmat!=0) delete [] mmat;
  mmat=0;
  mnrows=0;
  mncols=0;
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Alternative constructors:
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename T>
utils::Matrix<T>::Matrix(const int Nr, const int Nc){
  mnrows = Nr;
  mncols = Nc;
  int nelem = Nr * Nc;
  mmat = new T [nelem] ();
}
  

template <typename T>
utils::Matrix<T>::Matrix(const int Nr, const int Nc, const T & p){
  mnrows = Nr;
  mncols = Nc;
  int nelem = Nr * Nc;
  mmat = new T [nelem] ();
  for (int i=0; i<nelem; ++i) mmat[i]=p;
}
  
  

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Copy constructor:
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template <typename T>
utils::Matrix<T>::Matrix(const utils::Matrix<T> & sv){
  int Nr = sv.mnrows;
  int Nc = sv.mncols;
  if (Nr>0 && Nc>0){
    int nelem = Nr * Nc;
    mmat = new T [nelem] ();
    mnrows = Nr;
    mncols = Nc;
    for (int i=0; i<nelem; ++i)
      mmat[i] = sv.mmat[i];
  }
  else {
    mmat=0;
    mnrows=0;
    mncols=0;
  }
}





// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Operators:
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename T>
utils::Matrix<T> & utils::Matrix<T>::operator=(const utils::Matrix<T> & sv){
  int Nr = sv.nrows();
  int Nc = sv.ncols();

  if (this==&sv) return *this;

  if (Nr==0 || Nc==0){
    // Source not allocated. Return empty object.
    if (mnrows>0 || mncols>0) delete [] mmat;
    mmat=0;
    mnrows=0;
    mncols=0;
  }
  else {
    // Source is allocated.
    
    if ( (mnrows!=0 && mnrows!=Nr) || (mncols!=0 && mncols!=Nc) ){
      // Destination is allocated and of wrong shape. Delete it.
      delete [] mmat;
      mmat=0;
      mnrows=0;
      mncols=0;
    }


    int nelem = Nr * Nc;

    if (mnrows==0 || mncols==0){
      // Destination not allocated. Allocate it.
      mmat = new T [nelem] ();
      mnrows = Nr;
      mncols = Nc;
    }

    // Fill destination.
    for (int i=0; i<nelem; ++i)
      mmat[i] = sv.mmat[i]; // no problem to access private member mat in object sv here ... ?
    //mmat[Nc*i+j] = sv.elem(i,j); //mat[nrows*i+j]; // no problem to access private member mat in object sv here ... ?
  }
  return *this;
}
  


  


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Storage related functions:
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename T>
void utils::Matrix<T>::resize(const int Nr, const int Nc){
  if (mnrows==0 || mncols==0){
    // Allocate directly:
    mmat = new T [Nr*Nc] ();
    mnrows = Nr;
    mncols = Nc;
  }
  else {
    // Copy elements in original matrix into the same positions in the new one,
    // if they exist in the new matrix. The number of elements do not have to be
    // the same in both matrices.
    int Nrt = (mnrows<Nr) ? mnrows : Nr;
    int Nct = (mncols<Nc) ? mncols : Nc;

    T *mmat2 = new T [Nr*Nc] ();
    for (int i=0; i<Nrt; ++i){
      for (int j=0; j<Nct; ++j){
	mmat2[Nct*i+j] = mmat[Nct*i+j];
      }
    }
    delete [] mmat;

    mmat = new T [Nr*Nc] ();
    for (int i=0; i<Nrt; ++i){
      for (int j=0; j<Nct; ++j){
	mmat[Nct*i+j] = mmat2[Nct*i+j];
      }
    }
    delete [] mmat2;

    mnrows = Nr;
    mncols = Nc;
  }
}





template <typename T>
void utils::Matrix<T>::reshape(const int Nr, const int Nc){
  if (mnrows==0 || mncols==0){
    // Allocate directly:
    mmat = new T [Nr*Nc] ();
    for (int i=0; i<Nr*Nc; ++i) mmat[i]=T();
    mnrows = Nr;
    mncols = Nc;
  }
  else {
    // Copy elements in original matrix into the new one, in the serial order they
    // occur (row, column) in the oiginal matrix.
    if (mnrows*mncols != Nr*Nc){
      std::cout << "Different number of elements in original matrix and requested new matrix: "
		<< mnrows << " " << mncols << " versus "
		<< Nr << " " << Nc << ". Exiting." << std::endl;
      exit(EXIT_FAILURE);
    }

    T *mmat2 = new T [Nr*Nc] ();
    int k=0;
    for (int i=0; i<Nr; ++i){
      for (int j=0; j<Nc; ++j){
	mmat2[Nc*i+j] = mmat[k++];
      }
    }
    delete [] mmat;

    mmat = new T [Nr*Nc] ();
    for (int i=0; i<Nr; ++i){
      for (int j=0; j<Nc; ++j){
	mmat[Nc*i+j] = mmat2[Nc*i+j];
      }
    }
    delete [] mmat2;

    mnrows = Nr;
    mncols = Nc;
  }
}




// Number of rows:
template <typename T>
int utils::Matrix<T>::nrows() const {
  if (mmat!=0) return mnrows;
  else return 0;
}

// Number of columns:
template <typename T>
int utils::Matrix<T>::ncols() const {
  if (mmat!=0) return mncols;
  else return 0;
}


// Total size:
template <typename T>
int utils::Matrix<T>::size() const {
  if (mmat!=0) return mnrows*mncols;
  else return 0;
}




// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Matrix element operations:
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename T>
T & utils::Matrix<T>::elem(const int & i, const int & j) {

  if (i<0 || i>=mnrows){
    std::cout << "First index " << i << " is out of range " << mnrows-1 << ". Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
  if (j<0 || j>=mncols){
    std::cout << "Second index " << j << " is out of range " << mncols-1 << ". Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }

  return mmat[mncols*i+j];
}

template <typename T>
const T & utils::Matrix<T>::elem(const int & i, const int & j) const {
  if (i<0 || i>=mnrows){
    std::cout << "First index " << i << " is out of range " << mnrows-1 << ". Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
  if (j<0 || j>=mncols){
    std::cout << "Second index " << j << " is out of range " << mncols-1 << ". Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
  return mmat[mncols*i+j];
}


/*
template <typename T>
const T & utils::Matrix<T>::elem(const int & i, const int & j) const {
  if (i<0 || i>=mnrows){
    std::cout << "First index " << i << " is out of range. Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
  if (j<0 || j>=mncols){
    std::cout << "Second index " << j << " is out of range. Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
  return mmat[mncols*i+j];
}
*/



// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Matrix row and column operations:
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// Set a row:
template <typename T>
void utils::Matrix<T>::row(const int & irow, const Vector<T> & v){
  int Nc = v.size();

  if (mnrows==0 || mncols==0){
    // Object not allocated. Allocate to same size as input vector.
    mnrows = irow;
    mncols = Nc;
    mmat = new T [mnrows*mncols] ();
  }
  else {
    // Object is allocated.
    if (irow>=mnrows){
      std::cout << "Matrix does not contain the requested row " << irow << ". Exiting." << std::endl;
      exit(EXIT_FAILURE);
    }
    if (Nc!=mncols){
      std::cout << "Cannot set a row of length " << mncols
		<< " since length of the supplied vector is " << Nc << ". Exiting." << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  // Copy input data into object.
  for (int j=0; j<Nc; ++j) mmat[mncols*irow+j] = v[j];
}

// ***********************************************************************************

// Set a column:
template <typename T>
void utils::Matrix<T>::col(const int & icol, const Vector<T> & v){
  int Nr = v.size();

  if (mnrows==0 || mncols==0){
    // Object not allocated. Allocate to same size as input vector.
    mnrows = Nr;
    mncols = icol;
    mmat = new T [mnrows*mncols] ();
  }
  else {
    // Object is allocated.
    if (icol>=mncols){
      std::cout << "Matrix does not contain the requested column " << icol << ". Exiting." << std::endl;
      exit(EXIT_FAILURE);
    }
    if (Nr!=mnrows){
      std::cout << "Cannot set a column of length " << mnrows
		<< " since length of the supplied vector is " << Nr << ". Exiting." << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  // Copy input data into object.
  for (int j=0; j<Nr; ++j) mmat[mncols*j+icol] = v[j];
}

// ***********************************************************************************

// Get a row:
template <typename T>
Vector<T> utils::Matrix<T>::row(const int & irow) const {
  if (mnrows==0 || mncols==0){
    std::cout << "Cannot get a row from an unallocated matrix. Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
  if (irow >= mnrows){
    std::cout << "Specified row " << irow << " is outside matrix, whose last row is " << mnrows-1
	      << ". Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }

  Vector<T> v(mncols);
  for (int j=0; j<mncols; ++j)
    v[j] = mmat[mncols*irow+j];
  return v;
}

// ***********************************************************************************

// Get a column:
template <typename T>
Vector<T> utils::Matrix<T>::col(const int & icol) const {
  if (mnrows==0 || mncols==0){
    std::cout << "Cannot get a column from an unallocated matrix. Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
  if (icol >= mncols){
    std::cout << "Specified column " << icol << " is outside matrix, whose last column is " << mncols-1
	      << ". Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
  Vector<T> v(mnrows);
  for (int j=0; j<mnrows; ++j)
    v[j] = mmat[mncols*j+icol];
  return v;
}







// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Other matrix operations:
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// Set the matrix to unity matrix:
template <typename T>
void utils::Matrix<T>::unity(){
  if (mnrows==0 || mncols==0){
    std::cout << "Cannot set unallocated matrix of unknown size to unity matrix. Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
  if (mnrows!=mncols){
    std::cout << "Cannot set non-square matrix to unity matrix. Exiting." << std::endl;
    exit(EXIT_FAILURE);   
  }

  for (int i=0; i<mnrows; ++i){
    for (int j=0; j<mncols; ++j)
      mmat[mncols*i+j]=0;
    mmat[mncols*i+i]=1;
  }
}



// Get transpose:
template <typename T>
utils::Matrix<T> utils::Matrix<T>::transpose(){
  if (mnrows==0 || mncols==0){
    std::cout << "Cannot transpose unallocated matrix. Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
  utils::Matrix<T> tmpm(mncols, mnrows);
  for (int i=0; i<mnrows; ++i)
    for (int j=0; j<mncols; ++j)
      tmpm.elem(j,i) = elem(i,j);
  return tmpm;
}








template <typename T>
int utils::Matrix<T>::rank() const {
  if      (mnrows==1 || mncols==1) return 1;
  else {
    T elem;
    utils::Matrix<T> M;
    int N, Nr, Nc, i, j, k, m, row, nrowswaps=0, ncolswaps=0;
    T len=-1, maxlen=-1;

    Nr = mnrows;
    Nc = mncols;
    N = (Nr < Nc) ? Nr : Nc;

    M = *this;
    

    //std::cout << "Original matrix:" << std::endl;
    //std::cout << M << std::endl;

    for (k=0; k<N; ++k){

      // Row-pivoting:
      /* Find later row with largest absolute value in column 'k'. */
      row = k;
      for (m=k; m<N; ++m){
	
	len = M.elem(m, k);
	if (len<0) len *= (-1);
	if (m==k || (m>k && len>=maxlen)){
	  maxlen = len;
	  row    = m;
	}
      }

      /* Swap rows 'k' and 'row'. */
      if (k != row){
	T tmp;
	
	for (int i=0; i<Nc; ++i){
	  tmp           = M.elem(k,  i);
	  M.elem(k,  i) = M.elem(row,i);
	  M.elem(row,i) = tmp;
	}
	nrowswaps++;
	//std::cout << "After row swapping:" << std::endl;
	//std::cout << M << std::endl;
      }
      ncolswaps += 0;


      // Build a diagonal matrix:
      for (i=0; i<N; ++i){
	if (i==k) continue;
	
	if (double(utils::absval(M.elem(k,k))) < utils::eps_d()){
	  return k;
	}

	elem = M.elem(i,k)/M.elem(k,k);
	for (j=0; j<Nc; ++j){
	  M.elem(i,j) -= M.elem(k,j) * elem;
	}
      }
      //std::cout << "After row subtraction:" << std::endl;
      //std::cout << M << std::endl;
      

    }


    m=0;
    for (k=0; k<N; ++k)
      if (double(utils::absval(M.elem(k,k))) > utils::eps_d()) m++;
    return m;
  }
}




// Solve for single RHS.
template <typename T>
int utils::Matrix<T>::solve(Vector<T> & b,
			    Vector<T> & x,
			    utils::Matrix<T> & A_inv) const {
  utils::Matrix<T> B(mnrows, 1, 0);
  utils::Matrix<T> X(mnrows, 1, 0);
  int status;

  x.resize(mncols);
  if (b.size() != mnrows || mnrows != mncols){
    std::cout << "Cannot solve a system of linear equations with mismatched dimensions. Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }

  B.col(0, b);
  X.col(0, x);
  status = solve(B, X, A_inv);
  x = X.col(0);
  return status;
}




// Solve for multiple RHSs.
// Solve A X = B.
//
template <typename T>
int utils::Matrix<T>::solve(utils::Matrix<T> & B,
			    utils::Matrix<T> & X,
			    utils::Matrix<T> & A_inv) const {

  int i,icol=-1,irow=-1,j,k,l,ll,n=mnrows,m=B.ncols(), counter;
  T big,dum,pivinv;
  Vector<int> indxc(n,0), indxr(n,0), ipiv(n,0);

  A_inv = *this;
  X     = B;


  for (j=0; j<n; ++j) ipiv[j]=0;

  // ****************************************************************
  // Main loop over columns to be reduced:
  // ****************************************************************
  for (i=0; i<n; ++i){


    // Search for pivot element:
    for (j=0; j<n; ++j){ // Loop over rows.
      if (ipiv[j] == 1) continue;

      big     = 0;
      counter = 0;
      for (k=0; k<n; ++k){ // Loop over columns.
	if (ipiv[k] == 0){
	  if ( counter==0 || (counter>0 && abs(A_inv.elem(j,k)) > abs(big)) ){
	    counter++;
	    big = abs(A_inv.elem(j,k));
	    irow = j;
	    icol = k;
	  }
	}
      }

    }
    ++(ipiv[icol]);
    /* Pivot element is now known. Interchange rows if needed to put pivot
       on the diagonal. Columns are not interchanged, only relabeled.

       indxc[i], column of (i+t):th pivot element, is the (i+1):th column
       that is reduced.
       indxr[i] is the row where the pivot element was originally located.
       If indxr[i] != indxc[i] then there is an implied column interchange.
       With this form the solution in X will end up in the correct order, but
       the inverse matrix will be scrambled by columns.
    */
    if (irow != icol){
      for (l=0; l<n; ++l) swap(A_inv.elem(irow,l), A_inv.elem(icol,l));
      for (l=0; l<m; ++l) swap(X.elem(irow,l),     X.elem(icol,l));
    }
    indxr[i] = irow;
    indxc[i] = icol;



#if 0
    bool retry = true;
    int incr  = 0;
    while (retry){
      incr++;

      if (A_inv.elem(icol,icol) == 0){
	// Try adding a row from below:
	if (icol<n-1){
	  for (l=0; l<n; ++l) A_inv.elem(icol,incr) += A_inv.elem(icol+1,incr);
	  for (l=0; l<m; ++l) X.elem(icol,incr)     += X.elem(icol+1,incr);
	}
	else {
	  std::cout << "Singular matrix." << std::endl;
	  return 1;
	}
      }
      else
	retry = false;
    }
#endif

#if 1
    if (A_inv.elem(icol,icol) == 0){
      std::cout << "Singular matrix." << std::endl;
      return 1;
    }
#endif


    pivinv = 1.0/A_inv.elem(icol,icol);
    A_inv.elem(icol,icol) = 1;
    for (l=0; l<n; ++l) A_inv.elem(icol,l) *= pivinv;
    for (l=0; l<m; ++l) X.elem(icol,l)     *= pivinv;

    // Reduce the rows:
    for (ll=0; ll!=n; ++ll){
      if (ll != icol){
	dum = A_inv.elem(ll,icol);
	A_inv.elem(ll,icol) = 0;
	
	for (l=0; l<n; ++l) A_inv.elem(ll,l) -= A_inv.elem(icol,l) * dum;
	for (l=0; l<m; ++l) X.elem(ll,l)     -= X.elem(icol,l) * dum;
      }
    }

#if 0
    for (j=0; j<n; ++j){ // Loop over rows.
      for (k=0; k<n; ++k){ // Loop over columns.
	std::cout << " " << A_inv.elem(j,k);
      }
      std::cout << std::endl;
    }
#endif


  }
  // ****************************************************************
  // End of main loop
  // ****************************************************************

 
  // Unscramble the solution in terms of column interchanges, in the reverse
  // order that the permutation was built up:
  for (l=n-1; l>=0; --l){
    if (indxr[l] != indxc[l])
      for (k=0; k<n; ++k)
	swap( A_inv.elem(k,indxr[l]), A_inv.elem(k,indxc[l]));
  }
  return 0;

}



template <typename T>
int utils::Matrix<T>::inverse(utils::Matrix<T> & A_inv) const {
  utils::Matrix<T> B(nrows(), 0, 0);
  utils::Matrix<T> X(nrows(), 0, 0);

  return solve(B, X, A_inv);
}




template <typename T>
T utils::Matrix<T>::det() const {

  if (mnrows != mncols){
    std::cout << "Cannot calculate determinant for non-square matrices. Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }

  if      (mnrows==1) return elem(0,0);
  else if (mnrows==2) return elem(0,0)*elem(1,1) - elem(1,0)*elem(0,1);
  else {
    T dett, elem;
    utils::Matrix<T> M;
    int N, i, j, k, m, row, nrowswaps=0, ncolswaps=0;
    T len=-1, maxlen=-1;


    M = *this;
    N = M.nrows();

    //std::cout << "Original matrix:" << std::endl;
    //std::cout << M << std::endl;

    for (k=0; k<N; ++k){

      // Row-pivoting:
      /* Find later row with largest absolute value in column 'k'. */
      row = k;
      for (m=k; m<N; ++m){
	
	len = M.elem(m, k);
	if (len<0) len *= (-1);
	if (m==k || (m>k && len>=maxlen)){
	  maxlen = len;
	  row    = m;
	}
      }

      /* Swap rows 'k' and 'row'. */
      if (k != row){
	T tmp;
	
	for (int i=0; i<N; ++i){
	  tmp           = M.elem(k,  i);
	  M.elem(k,  i) = M.elem(row,i);
	  M.elem(row,i) = tmp;
	}
	nrowswaps++;
	//std::cout << "After row swapping:" << std::endl;
	//std::cout << M << std::endl;
      }
      ncolswaps += 0;


      // Build an upper triangular matrix:
      for (i=0; i<N; ++i){
	if (i==k) continue;
	
	if (double(utils::absval(M.elem(k,k))) < utils::eps_d()){
	  std::cout << "Pivot element is very small. Unable to continue matrix inversion. Exiting." << std::endl;
	  return 1;
	}

	elem = M.elem(i,k)/M.elem(k,k);
	for (j=0; j<N; ++j){
	  M.elem(i,j) -= M.elem(k,j) * elem;
	}
      }
      //std::cout << "After row subtraction:" << std::endl;
      //std::cout << M << std::endl;
      

    }


    dett = pow(-1.0, nrowswaps);
    for (k=0; k<N; ++k)
      dett *= M.elem(k,k);

    return dett;
  }

}










// Handle a + b:
template <typename T>
utils::Matrix<T> utils::operator+(const utils::Matrix<T> & a, const utils::Matrix<T> & b){
  if (a.nrows()!=b.nrows() || a.ncols()!=b.ncols()){
    std::cout << "Cannot add matrices of different shapes: "
	      << a.nrows() << " " << a.ncols() << " versus "
	      << b.nrows() << " " << b.ncols() << ". Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
  int Nr = a.nrows();
  int Nc = a.ncols();
  utils::Matrix<T> r(Nr, Nc);
  for (int i=0; i!=Nr; ++i)
    for (int j=0; j!=Nc; ++j)
      r.elem(i,j) = a.elem(i,j) + b.elem(i,j);
  return r;
}

// Handle a - b:
template <typename T>
utils::Matrix<T> utils::operator-(const utils::Matrix<T> & a, const utils::Matrix<T> & b){
  if (a.nrows()!=b.nrows() || a.ncols()!=b.ncols()){
    std::cout << "Cannot subtract matrices of different shapes: "
	      << a.nrows() << " " << a.ncols() << " versus "
	      << b.nrows() << " " << b.ncols() << ". Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
  int Nr = a.nrows();
  int Nc = a.ncols();
  utils::Matrix<T> r(Nr, Nc);
  for (int i=0; i!=Nr; ++i)
    for (int j=0; j!=Nc; ++j)
      r.elem(i,j) = a.elem(i,j) - b.elem(i,j);
  return r;
}

// Handle a * b:
template <typename T>
utils::Matrix<T> utils::operator*(const utils::Matrix<T> & a, const utils::Matrix<T> & b){
  if (a.ncols()!=b.nrows()){
    std::cout << "Cannot multiply matrices: First matrix is "
	      << a.nrows() << "x" << a.ncols() << " and second matrix is "
	      << b.nrows() << "x" << b.ncols() << ". Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
  int Nr = a.nrows();
  int Nc = b.ncols();
  int N  = a.ncols();
  utils::Matrix<T> r(Nr, Nc, 0);


  for (int i=0; i!=Nr; ++i){
    for (int j=0; j!=Nc; ++j){
      for (int k=0; k!=N; ++k){
	r.elem(i,j) += a.elem(i,k) * b.elem(k,j);
      }
    }
  }
  return r;
}



// Handle a * vec(b):
template <typename T>
Vector<T> utils::operator*(const utils::Matrix<T> & a, const Vector<T> & b){
  if (a.ncols()!=b.size()){
    std::cout << "Cannot multiply matrix and vector: Matrix is "
	      << a.nrows() << "x" << a.ncols() << " and vector is of length "
	      << b.size() << ". Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
  int Nr = a.nrows();
  int Nc = a.ncols();
  Vector<T> r(Nr, 0);

  for (int i=0; i<Nr; ++i){
    for (int j=0; j<Nc; ++j){
      r[i] += a.elem(i,j) * b[j];
    }
  }
  return r;
}



// Handle vec(a) * b:
template <typename T>
Vector<T> utils::operator*(const Vector<T> & a, const utils::Matrix<T> & b){
  if (a.size()!=b.nrows()){
    std::cout << "Cannot multiply matrix and vector: Vector is of length "
	      << a.size() << " and matrix is "
	      << a.nrows() << "x" << a.ncols() << ". Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
  int Nr = a.size();
  int Nc = b.ncols();
  Vector<T> r(Nr,0);

  for (int i=0; i!=Nr; ++i){
    for (int j=0; j!=Nc; ++j){
      r[i] += a[i] * b.elem(i,j);
    }
  }
  return r;
}






// Handle a * b, b being S type:
template <typename T, typename S>
utils::Matrix<T> utils::operator*(const utils::Matrix<T> & a, const S & b){
  int Nr = a.nrows();
  int Nc = a.ncols();
  T bp = T(b);
  utils::Matrix<T> r(Nr, Nc, 0);

  for (int i=0; i!=Nr; ++i)
    for (int j=0; j!=Nc; ++j)
      r.elem(i,j) = a.elem(i,j) * bp;
  return r;
}

// Handle a * b, a being S type:
template <typename T, typename S>
utils::Matrix<T> utils::operator*(const S & a, const utils::Matrix<T> & b){
  int Nr = b.nrows();
  int Nc = b.ncols();
  T ap = T(a);
  utils::Matrix<T> r(Nr, Nc, 0);
  for (int i=0; i!=Nr; ++i)
    for (int j=0; j!=Nc; ++j)
      r.elem(i,j) = ap * b.elem(i,j);
  return r;
}


// Handle a / b, b being S type:
template <typename T, typename S>
utils::Matrix<T> utils::operator/(const utils::Matrix<T> & a, const S & b){
  int Nr = a.nrows();
  int Nc = a.ncols();
  T bp = T(b);
  utils::Matrix<T> r(Nr, Nc, 0);

  for (int i=0; i!=Nr; ++i)
    for (int j=0; j!=Nc; ++j)
      r.elem(i,j) = a.elem(i,j) / bp;
  return r;
}



// Handle printing of Matrix<T> objects:
template <typename T, typename U>
U & utils::operator<<(U & os, const utils::Matrix<T> & sv){
  int Nr = sv.nrows();
  int Nc = sv.ncols();

  os << std::endl;
  for (int i=0; i!=Nr; ++i){
    for (int j=0; j!=Nc; ++j){
      os << "  " << sv.elem(i,j);
      os.clear();
    }
    os << std::endl;
  }
  //os << " ";
  return os;
}


#endif







#if 0

  // Perform partial pivoting on a matrix:
  template <typename T>
  void pivoting(int k,
		Matrix<T> & A,
		Matrix<T> & B,
		Vector<int> & row_id,
		Vector<int> & col_id,
		int & nrowswaps,
		int & ncolswaps
		);



template <typename T>
void utils::pivoting(int k,
		     utils::Matrix<T> & A,
		     utils::Matrix<T> & B,
		     Vector<int> & row_id,
		     Vector<int> & col_id,
		     int & nrowswaps,
		     int & ncolswaps
		     ){
  
  int N, Nr, Nc, m, row, swap1, swap2;
  T len, maxlen;
    
  maxlen = -1;
  row = k;

  Nr = A.nrows();
  Nc = A.ncols();
  N = (Nr < Nc) ? Nr : Nc;



  /* Perform partial pivoting. */
  /* Find later row with largest absolute value in column 'k'. */
  for (m=k; m!=N; ++m){
    len = A.elem(m, k);
    if (len<0) len *= (-1);
    if (m==k || (m>k && len>=maxlen)){
      maxlen = len;
      row    = m;
    }
  }

  /* Swap rows 'k' and 'row'. */
  swap1 = k;
  swap2 = row;

  /* Swapping of rows: */
  if (swap1 != swap2){
    T tmp1;
    int tmpl;

    for (int i=0; i!=Nc; ++i){
      tmp1            = A.elem(swap1,i);
      A.elem(swap1,i) = A.elem(swap2,i);
      A.elem(swap2,i) = tmp1;
	
      tmp1            = B.elem(swap1,i);
      B.elem(swap1,i) = B.elem(swap2,i);
      B.elem(swap2,i) = tmp1;

      col_id[i] = i;
      ncolswaps += 0;
    }
    // Original row row_id[swap1] is moved to row swap2, and vice versa:
    tmpl          = row_id[swap1];
    row_id[swap1] = row_id[swap2];
    row_id[swap2] = tmpl;

    nrowswaps++;
  }



template <typename T>
int utils::Matrix<T>::invertmatrix(utils::Matrix<T> & A_inv) const {
  utils::Matrix<T> M;
  T elem;
  int N, i, j, k, nrowswaps, ncolswaps;
  Vector<int> row_id, col_id;

  if (mnrows != mncols){
    std::cout << "Cannot invert non-square matrices. Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }


  /*  
      M = A (original)
      Row operations are carried out on the augmented matrix I | M, I = unit matrix.
      In the end this matrix has become A_inv | I.
  */

  nrowswaps=0;
  ncolswaps=0;

  N = mnrows;
  M = *this;
  A_inv.resize(N, N);
  row_id.resize(N);
  col_id.resize(N);
  A_inv.unity();
  for (int i=0; i!=N; ++i)
    row_id[i] = i;


  for (k=0; k!=N; ++k){
    /* Pivoting: */
    utils::pivoting(k, M, A_inv, row_id, col_id, nrowswaps, ncolswaps);
    // std::cout << row_id << std::endl;
    
    /* Normalize element M[k][k] to 1.0. */
    /* Avoid dividing with 0.0:          */
    if (double(utils::absval(M.elem(k,k))) < utils::eps_d()){
      std::cout << "Pivot element is very small. Unable to continue matrix inversion. Exiting." << std::endl;
      return 1;
    }
    
    elem = 1.0/M.elem(k, k);
    for (j=0; j!=N; ++j){
      if (j==k) M.elem(k,j) = 1.0;
      else      M.elem(k,j) *= elem;
      
      A_inv.elem(k,j) *= elem;
    }


    /*
      To each row, add some other row k in the matrix
      multiplied with with some element A_bak[i][k]:
    */


    for (i=0; i!=N; ++i){
      if (i==k) continue;

      elem = M.elem(i,k);
      for (j=0; j!=N; ++j){
	M.elem(i,j)     -= M.elem(k,j) * elem;
	A_inv.elem(i,j) -= A_inv.elem(k,j) * elem;
      }
    }
    
    //std::cout << "After pivoting and modification for k = " << k << ":" << std::endl;
    //std::cout << M << std::endl;
  }

 
  
  return 0;
}






template <typename T>
T utils::Matrix<T>::determinant_LM() const {
  if      (mnrows==1) return elem(0,0);
  else if (mnrows==2) return elem(0,0)*elem(1,1) - elem(1,0)*elem(0,1);
  else {
    T det=0;

    //std::cout << "Processing matrix of nrows " << mnrows << std::endl;
 
    for (int i=0; i<mnrows; ++i){
      for (int j=0; j<mnrows; ++j){

	// ----------------------------------------------------------------
	// Make the sub-matrix:
	// ----------------------------------------------------------------
	utils::Matrix<T> tmp(mnrows-1, 0);
	// Row and column which are struck out.
	int iso=i, jso=j;
	// Temporary counters.
	int k, l, it=0, jt=0;


	iso=i;
	jso=j;
	//std::cout<<"Struck out element indices: " << iso << " " << jso <<std::endl;

	it=0;
	for (k=0; k<mnrows; ++k){
	  if (k==iso) continue;
	  //std::cout<<"Made it here 01, k="<<k<<std::endl;

	  jt=0;
	  for (l=0; l<mnrows; ++l){
	    if (l==jso) continue;
	    //std::cout<<"Made it here 01, k,l="<<k<<" "<<l<<std::endl;
	    //std::cout<<"  it, jt = " << it << " " << jt << std::endl;
	    tmp.elem(it, jt) = elem(k, l);
	    jt++;
	  }
	  it++;
	}
	//std::cout << "mnrows=" << mnrows << std::endl;

	det = det + pow(-1.0, i+j) * elem(i, j) * tmp.determinant();
      }
    }
    return det;
  }
}


#endif

