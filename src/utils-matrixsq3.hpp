



#ifndef UTILS_MATRIXSQ3_HPP
#define UTILS_MATRIXSQ3_HPP



#include <iostream>
#include <string>
#include <sstream>

#include <cstdlib>

#include "utils.hpp"
#include "utils-math.hpp"
#include "utils-vector3.hpp"


using utils::Vector3;
using utils::abs;
using utils::swap;


namespace utils {

  template <typename T>
  class MatrixSq3 {

  private:
    T mmat[9];
    
  public:
    MatrixSq3();
    ~MatrixSq3();

    MatrixSq3(const T & p);

    MatrixSq3(const MatrixSq3<T> & sv);

    // Operators:
    MatrixSq3<T> & operator=(const MatrixSq3<T> & sv);

    // MatrixSq3 element operations:
    inline       T & elem(const int & i, const int & j);
    inline const T & elem(const int & i, const int & j) const;

    // MatrixSq3 row and column operations:
    // Set a row:
    void row(const int & i, const Vector3<T> & v);
    // Get a row:
    Vector3<T> row(const int & i) const;
    // Set a column:
    void col(const int & i, const Vector3<T> & v);
    // Get a column:
    Vector3<T> col(const int & i) const;

    // Other matrix operations:

    // Set object to unity matrix:
    void unity();
    // Get transpose:
    MatrixSq3<T> transpose();

    // Get rank (number of independent rows):
    int rank() const;

    // Invert matrix:
    /* MatrixSq3 inversion can be done by Gauss-Jordan elimination and
       back-substitution. Another method is to perform LU decomposition
       first and then calculate the inverse using the decomposition
       matrices. Both methods can be used to solve a system of linear
       equations, so that the inverse matrix is produced as a side
       effect.
    */
    int solve(Vector3<T> & b,
	      Vector3<T> & x,
	      MatrixSq3<T> & A_inv) const;
    int solve(MatrixSq3<T> & B,
	      MatrixSq3<T> & X,
	      MatrixSq3<T> & A_inv) const;
    int inverse(MatrixSq3<T> & A_inv) const;
    T det() const;

   
  } ;


  // Handle a + b:
  template <typename T>
  MatrixSq3<T> operator+(const MatrixSq3<T> & a, const MatrixSq3<T> & b);
  
  // Handle a - b:
  template <typename T>
  MatrixSq3<T> operator-(const MatrixSq3<T> & a, const MatrixSq3<T> & b);

  // Handle a * b:
  template <typename T>
  MatrixSq3<T> operator*(const MatrixSq3<T> & a, const MatrixSq3<T> & b);

  // Handle a * vec(b):
  template <typename T>
  Vector3<T> operator*(const MatrixSq3<T> & a, const Vector3<T> & b);

  // Handle vec(a) * b:
  template <typename T>
  Vector3<T> operator*(const Vector3<T> & a, const MatrixSq3<T> & b);

  // Handle a * b, b being S type:
  template <typename T, typename S>
  MatrixSq3<T> operator*(const MatrixSq3<T> & a, const S & b);

  // Handle a * b, a being S type:
  template <typename T, typename S>
  MatrixSq3<T> operator*(const S & a, const MatrixSq3<T> & b);

  // Handle a / b, b being S type:
  template <typename T, typename S>
  MatrixSq3<T> operator/(const MatrixSq3<T> & a, const S & b);


  // Handle printing of MatrixSq3<T> objects:
  template <typename T, typename U>
  U & operator<<(U & os, const MatrixSq3<T> & sv);

  
}




/* #########################################################################
   Definitions:
   #########################################################################
*/

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Default ctor:
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template <typename T>
utils::MatrixSq3<T>::MatrixSq3(){
  for (int i=0; i<9; ++i) mmat[i]=0;
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Default dtor:
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template <typename T>
utils::MatrixSq3<T>::~MatrixSq3(){
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Alternative constructors:
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template <typename T>
utils::MatrixSq3<T>::MatrixSq3(const T & p){
  for (int i=0; i<9; ++i) mmat[i]=p;
}

  
  

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Copy constructor:
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template <typename T>
utils::MatrixSq3<T>::MatrixSq3(const utils::MatrixSq3<T> & sv){
  for (int i=0; i<9; ++i) mmat[i] = sv.mmat[i];
}





// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Operators:
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename T>
utils::MatrixSq3<T> & utils::MatrixSq3<T>::operator=(const utils::MatrixSq3<T> & sv){
  if (this==&sv) return *this;

  for (int i=0; i<9; ++i) mmat[i] = sv.mmat[i];

  return *this;
}
  


  





// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// MatrixSq3 element operations:
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename T>
T & utils::MatrixSq3<T>::elem(const int & i, const int & j) {
  if (i<0 || i>=3){
    std::cout << "Square matrix 3x3: First index " << i << " is out of range. Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
  if (j<0 || j>=3){
    std::cout << "Square matrix 3x3: Second index " << j << " is out of range. Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
  return mmat[3*i+j];
}

template <typename T>
const T & utils::MatrixSq3<T>::elem(const int & i, const int & j) const {
  if (i<0 || i>=3){
    std::cout << "Square matrix 3x3: First index " << i << " is out of range. Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
  if (j<0 || j>=3){
    std::cout << "Square matrix 3x3: Second index " << j << " is out of range. Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
  return mmat[3*i+j];
}




// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// MatrixSq3 row and column operations:
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// Set a row:
template <typename T>
void utils::MatrixSq3<T>::row(const int & irow, const Vector3<T> & v){
  if (irow<0 || irow>=3){
    std::cout << "Square matrix 3x3: Row index " << irow << " is out of range. Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
  // Copy input data into object.
  for (int j=0; j<3; ++j) mmat[3*irow+j] = v[j];
}

// ***********************************************************************************

// Set a column:
template <typename T>
void utils::MatrixSq3<T>::col(const int & icol, const Vector3<T> & v){
  if (icol<0 || icol>=3){
    std::cout << "Square matrix 3x3: Column index " << icol << " is out of range. Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
  // Copy input data into object.
  for (int j=0; j<3; ++j) mmat[3*j+icol] = v[j];
}

// ***********************************************************************************

// Get a row:
template <typename T>
Vector3<T> utils::MatrixSq3<T>::row(const int & irow) const {
  if (irow<0 || irow>=3){
    std::cout << "Square matrix 3x3: Row index " << irow << " is out of range. Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
  Vector3<T> v(0);
  for (int j=0; j<3; ++j)
    v[j] = mmat[3*irow+j];
  return v;
}

// ***********************************************************************************

// Get a column:
template <typename T>
Vector3<T> utils::MatrixSq3<T>::col(const int & icol) const {
  if (icol<0 || icol>=3){
    std::cout << "Square matrix 3x3: Column index " << icol << " is out of range. Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
  Vector3<T> v(0);
  for (int j=0; j<3; ++j)
    v[j] = mmat[3*j+icol];
  return v;
}







// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Other matrix operations:
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// Set the matrix to unity matrix:
template <typename T>
void utils::MatrixSq3<T>::unity(){
  for (int i=0; i<3; ++i){
    for (int j=0; j<3; ++j)
      mmat[3*i+j]=0;
    mmat[3*i+i]=1;
  }
}



// Get transpose:
template <typename T>
utils::MatrixSq3<T> utils::MatrixSq3<T>::transpose(){
  utils::MatrixSq3<T> tmpm;
  for (int i=0; i<3; ++i)
    for (int j=0; j<3; ++j)
      tmpm.elem(j,i) = elem(i,j);
  return tmpm;
}








template <typename T>
int utils::MatrixSq3<T>::rank() const {
  T elem;
  utils::MatrixSq3<T> M(0);
  int N, Nr, Nc, i, j, k, m, row, nrowswaps=0, ncolswaps=0;
  T len=-1, maxlen=-1;

  Nr = 3;
  Nc = 3;
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





// Solve for single RHS.
template <typename T>
int utils::MatrixSq3<T>::solve(Vector3<T> & b,
			       Vector3<T> & x,
			       utils::MatrixSq3<T> & A_inv) const {
  utils::MatrixSq3<T> B(0);
  utils::MatrixSq3<T> X(0);
  int status;

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
int utils::MatrixSq3<T>::solve(utils::MatrixSq3<T> & B,
			       utils::MatrixSq3<T> & X,
			       utils::MatrixSq3<T> & A_inv) const {

  int i,icol=-1,irow=-1,j,k,l,ll,n=3,m=3, counter;
  T big,dum,pivinv;
  Vector3<int> indxc(0), indxr(0), ipiv(0);

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
int utils::MatrixSq3<T>::inverse(utils::MatrixSq3<T> & A_inv) const {
  utils::MatrixSq3<T> B(0);
  utils::MatrixSq3<T> X(0);

  return solve(B, X, A_inv);
}




template <typename T>
T utils::MatrixSq3<T>::det() const {
  T dett, elem;
  utils::MatrixSq3<T> M(0);
  int N, i, j, k, m, row, nrowswaps=0, ncolswaps=0;
  T len=-1, maxlen=-1;


  M = *this;
  N = 3;

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





// Handle a + b:
template <typename T>
utils::MatrixSq3<T> utils::operator+(const utils::MatrixSq3<T> & a, const utils::MatrixSq3<T> & b){
  utils::MatrixSq3<T> r(0);
  for (int i=0; i<3; ++i)
    for (int j=0; j<3; ++j)
      r.elem(i,j) = a.elem(i,j) + b.elem(i,j);
  return r;
}

// Handle a - b:
template <typename T>
utils::MatrixSq3<T> utils::operator-(const utils::MatrixSq3<T> & a, const utils::MatrixSq3<T> & b){
  utils::MatrixSq3<T> r(0);
  for (int i=0; i<3; ++i)
    for (int j=0; j<3; ++j)
      r.elem(i,j) = a.elem(i,j) - b.elem(i,j);
  return r;
}

// Handle a * b:
template <typename T>
utils::MatrixSq3<T> utils::operator*(const utils::MatrixSq3<T> & a, const utils::MatrixSq3<T> & b){
  utils::MatrixSq3<T> r(0);
  for (int i=0; i<3; ++i){
    for (int j=0; j<3; ++j){
      for (int k=0; k<3; ++k){
	r.elem(i,j) += a.elem(i,k) * b.elem(k,j);
      }
    }
  }
  return r;
}



// Handle a * vec(b):
template <typename T>
Vector3<T> utils::operator*(const utils::MatrixSq3<T> & a, const Vector3<T> & b){
  Vector3<T> r(0);
  for (int i=0; i<3; ++i){
    for (int j=0; j<3; ++j){
      r[i] += a.elem(i,j) * b[j];
    }
  }
  return r;
}



// Handle vec(a) * b:
template <typename T>
Vector3<T> utils::operator*(const Vector3<T> & a, const utils::MatrixSq3<T> & b){
  Vector3<T> r(0);
  for (int i=0; i<3; ++i){
    for (int j=0; j<3; ++j){
      r[i] += a[i] * b.elem(i,j);
    }
  }
  return r;
}






// Handle a * b, b being S type:
template <typename T, typename S>
utils::MatrixSq3<T> utils::operator*(const utils::MatrixSq3<T> & a, const S & b){
  T bp = T(b);
  utils::MatrixSq3<T> r(0);
  for (int i=0; i<3; ++i)
    for (int j=0; j<3; ++j)
      r.elem(i,j) = a.elem(i,j) * bp;
  return r;
}

// Handle a * b, a being S type:
template <typename T, typename S>
utils::MatrixSq3<T> utils::operator*(const S & a, const utils::MatrixSq3<T> & b){
  T ap = T(a);
  utils::MatrixSq3<T> r(0);
  for (int i=0; i<3; ++i)
    for (int j=0; j<3; ++j)
      r.elem(i,j) = ap * b.elem(i,j);
  return r;
}


// Handle a / b, b being S type:
template <typename T, typename S>
utils::MatrixSq3<T> utils::operator/(const utils::MatrixSq3<T> & a, const S & b){
  T bp = T(b);
  utils::MatrixSq3<T> r(0);
  for (int i=0; i<3; ++i)
    for (int j=0; j<3; ++j)
      r.elem(i,j) = a.elem(i,j) / bp;
  return r;
}



// Handle printing of MatrixSq3<T> objects:
template <typename T, typename U>
U & utils::operator<<(U & os, const utils::MatrixSq3<T> & sv){
  os << std::endl;
  for (int i=0; i<3; ++i){
    for (int j=0; j<3; ++j){
      os << "  " << sv.elem(i,j);
      os.clear();
    }
    os << std::endl;
  }
  //os << " ";
  return os;
}


#endif





