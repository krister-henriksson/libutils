



#ifndef UTILS_MAT3_HPP
#define UTILS_MAT3_HPP



#include <iostream>
#include <string>
#include <sstream>

#include <cstdlib>

#include "utils.hpp"




namespace utils {

  /* #########################################################################
     Declaration of Matrix3<T> template class:
     #########################################################################
  */

  template <typename T>
  class Matrix3 {

  private:
    T *mmat3;
    int mn1;
    int mn2;
    int mn3;

  public:
    Matrix3();
    ~Matrix3();


    Matrix3(const int N1, const int N2, const int N3);
    Matrix3(const int N1, const int N2, const int N3, const T & p);

    Matrix3(const Matrix3<T> & sv);

    // Operators:
    Matrix3<T> & operator=(const Matrix3<T> & sv);

    // Storage related functions:
    void resize(const int N1, const int N2, const int N3);
    void reshape(const int N1, const int N2, const int N3);
    // Linear size:
    int n1() const;
    int n2() const;
    int n3() const;
    // Total size:
    int size() const;

    // Matrix element operations:
    inline       T & elem(const int & i, const int & j, const int & k);
    inline const T & elem(const int & i, const int & j, const int & k) const;

    // Other matrix operations:

    // Set object to unity Matrix3rix:
    void unity();

   
  } ;


  // Handle a + b:
  template <typename T>
  Matrix3<T> operator+(const Matrix3<T> & a, const Matrix3<T> & b);
  
  // Handle a - b:
  template <typename T>
  Matrix3<T> operator-(const Matrix3<T> & a, const Matrix3<T> & b);

  // Handle a * b, b being S type:
  template <typename T, typename S>
  Matrix3<T> operator*(const Matrix3<T> & a, const S & b);

  // Handle a * b, a being S type scalar:
  template <typename T, typename S>
  Matrix3<T> operator*(const S & a, const Matrix3<T> & b);

  // Handle a / b, b being S type scalar:
  template <typename T, typename S>
  Matrix3<T> operator/(const Matrix3<T> & a, const S & b);

  // Handle printing of Matrix3<T> objects:
  template <typename T, typename U>
  U & operator<<(U & os, const Matrix3<T> & sv);

  
}




/* #########################################################################
   Definitions:
   #########################################################################
*/

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Default ctor:
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template <typename T>
utils::Matrix3<T>::Matrix3(){
  mmat3=0;
  mn1=0;
  mn2=0;
  mn3=0;
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Default dtor:
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template <typename T>
utils::Matrix3<T>::~Matrix3(){
  if (mmat3!=0) delete [] mmat3;
  mmat3=0;
  mn1=0;
  mn2=0;
  mn3=0;
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Alternative constructors:
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename T>
utils::Matrix3<T>::Matrix3(const int N1, const int N2, const int N3){
  mn1 = N1;
  mn2 = N2;
  mn3 = N3;
  int nelem = N1 * N2 * N3;
  mmat3 = new T [nelem] ();
  for (int i=0; i!=nelem; ++i) mmat3[i]=T(0);
}
  
template <typename T>
utils::Matrix3<T>::Matrix3(const int N1, const int N2, const int N3, const T & p){
  mn1 = N1;
  mn2 = N2;
  mn3 = N3;
  int nelem = N1 * N2 * N3;
  mmat3 = new T [nelem] ();
  for (int i=0; i!=nelem; ++i) mmat3[i]=p;
}
  
  

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Copy constructor:
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template <typename T>
utils::Matrix3<T>::Matrix3(const utils::Matrix3<T> & sv){
  int N1 = sv.mn1;
  int N2 = sv.mn2;
  int N3 = sv.mn3;
  if (N1*N2*N3>0){
    mmat3 = new T [N1*N2*N3] ();
    mn1 = N1;
    mn2 = N2;
    mn3 = N3;
    for (int i=0; i!=N1; ++i)
      for (int j=0; j!=N2; ++j)
	for (int k=0; k!=N3; ++k){
	  mmat3[i*N2*N3 + j*N3 + k] = sv.mmat3[i*N2*N3 + j*N3 + k];
	}
  }
  else {
    mmat3=0;
    mn1=0;
    mn2=0;
    mn3=0;
  }
}





// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Operators:
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename T>
utils::Matrix3<T> & utils::Matrix3<T>::operator=(const utils::Matrix3<T> & sv){
  int N1 = sv.n1();
  int N2 = sv.n2();
  int N3 = sv.n3();

  if (this==&sv) return *this;

  if (N1==0){
    // Source not allocated. Return empty object.
    if (mn1*mn2*mn3>0) delete [] mmat3;
    mmat3=0;
    mn1=0;
    mn2=0;
    mn3=0;
  }
  else {
    // Source is allocated.
    
    if ( (mn1!=0 && mn1!=N1) ||
	 (mn2!=0 && mn2!=N2) ||
	 (mn3!=0 && mn3!=N3) ){
      // Destination is allocated and of wrong shape. Delete it.
      delete [] mmat3;
      mmat3=0;
      mn1=0;
      mn2=0;
      mn3=0;
    }

    if (mn1==0 || mn2==0 || mn3==0){
      // Destination not allocated. Allocate it.
      mmat3 = new T [N1*N2*N3] ();
      mn1 = N1;
      mn2 = N2;
      mn3 = N3;
    }

    // Fill destination.
    for (int i=0; i!=N1; ++i)
      for (int j=0; j!=N2; ++j)
	for (int k=0; k!=N3; ++k){
	  mmat3[i*N2*N3 + j*N3 + k] = sv.mmat3[i*N2*N3 + j*N3 + k];
	}
    //mat3[nrows*i+j]; // no problem to access private member mat3 in object sv here ... ?
    
  }
  return *this;
}
  


  


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Storage related functions:
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename T>
void utils::Matrix3<T>::resize(const int N1, const int N2, const int N3){
  if (mn1==0 || mn2==0 || mn3==0){
    // Allocate directly:
    mmat3 = new T [N1*N2*N3] ();
    mn1 = N1;
    mn2 = N2;
    mn3 = N3;
  }
  else {
    // Copy elements in original matrix into the same positions in the new one,
    // if they exist in the new matrix. The number of elements do not have to be
    // the same in both mat3rices.
    int N1t = (mn1<N1) ? mn1 : N1;
    int N2t = (mn2<N2) ? mn2 : N2;
    int N3t = (mn3<N3) ? mn3 : N3;

    T *mmat32 = new T [N1*N2*N3] ();
    for (int i=0; i!=N1t; ++i)
      for (int j=0; j!=N2t; ++j)
	for (int k=0; k!=N3t; ++k){
	  mmat32[i*N2*N3 + j*N3 + k] = mmat3[i*mn2*mn3 + j*mn3 + k];
	}
    delete [] mmat3;

    mmat3 = new T [N1*N2*N3] ();
    for (int i=0; i!=N1t; ++i){
      for (int j=0; j!=N2t; ++j){
	for (int k=0; k!=N3t; ++k){
	  mmat3[i*N2*N3 + j*N3 + k] = mmat32[i*N2*N3 + j*N3 + k];
	}
      }
    }
    delete [] mmat32;
    
    mn1 = N1;
    mn2 = N2;
    mn3 = N3;
  }
}





template <typename T>
void utils::Matrix3<T>::reshape(const int N1, const int N2, const int N3){
  if (mn1==0 || mn2==0 || mn3==0){
    // Allocate directly:
    mmat3 = new T [N1*N2*N3] ();
    for (int i=0; i!=N1*N2*N3; ++i) mmat3[i]=T();
    mn1 = N1;
    mn2 = N2;
    mn3 = N3;
  }
  else {
    // Copy elements in original matrix into the new one, in the serial order they
    // occur (row, column) in the oiginal matrix.
    if (mn1*mn2*mn3 != N1*N2*N3){
      std::cout << "Matrix reshaping: Different number of elements in original matrix and requested new matrix. Exiting." << std::endl;
      exit(EXIT_FAILURE);
    }

    T *mmat32 = new T [N1*N2*N3] ();
    int counter=0;
    for (int i=0; i!=N1; ++i){
      for (int j=0; j!=N2; ++j){
	for (int k=0; k!=N3; ++k){
	  mmat32[i*N2*N3 + j*N3 + k] = mmat3[counter++];
	}
      }
    }
    delete [] mmat3;

    mmat3 = new T [N1*N2*N3] ();
    for (int i=0; i!=N1; ++i){
      for (int j=0; j!=N2; ++j){
	for (int k=0; k!=N3; ++k){
	  mmat3[i*N2*N3 + j*N3 + k] = mmat32[i*N2*N3 + j*N3 + k];
	}
      }
    }
    delete [] mmat32;

    mn1 = N1;
    mn2 = N2;
    mn3 = N3;
  }
}




// Number of elements in dimension i:
template <typename T>
int utils::Matrix3<T>::n1() const {
  if (mmat3!=0) return mn1;
  else return 0;
}
template <typename T>
int utils::Matrix3<T>::n2() const {
  if (mmat3!=0) return mn2;
  else return 0;
}
template <typename T>
int utils::Matrix3<T>::n3() const {
  if (mmat3!=0) return mn3;
  else return 0;
}
// Total size:
template <typename T>
int utils::Matrix3<T>::size() const {
  if (mmat3!=0) return mn1*mn2*mn3;
  else return 0;
}






// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Matrix element operations:
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename T>
T & utils::Matrix3<T>::elem(const int & i, const int & j, const int & k) {
  if (i<0 || i>=mn1){
    std::cout << "First index " << i << " is out of range. Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
  if (j<0 || j>=mn2){
    std::cout << "Second index " << j << " is out of range. Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
  if (k<0 || k>=mn3){
    std::cout << "Third index " << k << " is out of range. Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }

  return mmat3[i*mn2*mn3 + j*mn3 + k];
}

template <typename T>
const T & utils::Matrix3<T>::elem(const int & i, const int & j, const int & k) const {
  if (i<0 || i>=mn1){
    std::cout << "First index " << i << " is out of range. Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
  if (j<0 || j>=mn2){
    std::cout << "Second index " << j << " is out of range. Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
  if (k<0 || k>=mn3){
    std::cout << "Third index " << k << " is out of range. Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }

  return mmat3[i*mn2*mn3 + j*mn3 + k];
}









// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Other matrix operations:
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// Set the matrix to unity matrix:
template <typename T>
void utils::Matrix3<T>::unity(){
  if (mn1==0 || mn2==0 || mn3==0){
    std::cout << "Cannot set unallocated matrix of unknown size to unity matrix. Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
  if (mn1!=mn2 || mn2!=mn3 || mn1!=mn3){
    std::cout << "Cannot set non-cubic mat3 object to 'unity cube'. Exiting." << std::endl;
    exit(EXIT_FAILURE);   
  }

  for (int i=0; i!=mn1; ++i){
    for (int j=0; j!=mn2; ++j){
      for (int k=0; k!=mn3; ++k){
	mmat3[i*mn2*mn3 + j*mn3 + k] = 0;
      }
    }
  }
  for (int i=0; i!=mn1; ++i)
    mmat3[i*mn2*mn3 + i*mn3 + i] = 1;
}
















// Handle a + b:
template <typename T>
utils::Matrix3<T> utils::operator+(const utils::Matrix3<T> & a, const utils::Matrix3<T> & b){
  if (! ( a.n1()==b.n1() &&
	  a.n2()==b.n2() &&
	  a.n3()==b.n3() ) ){
    std::cout << "Cannot add matrices of different shape. Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
  int N1 = a.n1();
  int N2 = a.n2();
  int N3 = a.n3();
  utils::Matrix3<T> r(N1, N2, N3, 0);
  for (int i=0; i!=N1; ++i)
    for (int j=0; j!=N2; ++j)
      for (int k=0; k!=N3; ++k)
	r.elem(i,j,k) = a.elem(i,j,k) + b.elem(i,j,k);
  return r;
}



// Handle a - b:
template <typename T>
utils::Matrix3<T> utils::operator-(const utils::Matrix3<T> & a, const utils::Matrix3<T> & b){
  if (! ( a.n1()==b.n1() &&
	  a.n2()==b.n2() &&
	  a.n3()==b.n3() ) ){
    std::cout << "Cannot subtract matrices of different shape. Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
  int N1 = a.n1();
  int N2 = a.n2();
  int N3 = a.n3();
  utils::Matrix3<T> r(N1, N2, N3, 0);
  for (int i=0; i!=N1; ++i)
    for (int j=0; j!=N2; ++j)
      for (int k=0; k!=N3; ++k)
	r.elem(i,j,k) = a.elem(i,j,k) - b.elem(i,j,k);
  return r;
}



// Handle a * b, b being S type:
template <typename T, typename S>
utils::Matrix3<T> utils::operator*(const utils::Matrix3<T> & a, const S & b){
  int N1 = a.n1();
  int N2 = a.n2();
  int N3 = a.n3();
  T bp = T(b);
  utils::Matrix3<T> r(N1, N2, N3, 0);
  for (int i=0; i!=N1; ++i)
    for (int j=0; j!=N2; ++j)
      for (int k=0; k!=N3; ++k)
	r.elem(i,j,k) = a.elem(i,j,k) * bp;
  return r;
}

// Handle a * b, a being T type:
template <typename T, typename S>
utils::Matrix3<T> utils::operator*(const S & a, const utils::Matrix3<T> & b){
  int N1 = b.n1();
  int N2 = b.n2();
  int N3 = b.n3();
  T ap = T(a);
  utils::Matrix3<T> r(N1, N2, N3, 0);
  for (int i=0; i!=N1; ++i)
    for (int j=0; j!=N2; ++j)
      for (int k=0; k!=N3; ++k)
	r.elem(i,j,k) = ap * b.elem(i,j,k);
  return r;
}


// Handle a / b, b being T type:
template <typename T, typename S>
utils::Matrix3<T> utils::operator/(const utils::Matrix3<T> & a, const S & b){
  int N1 = a.n1();
  int N2 = a.n2();
  int N3 = a.n3();
  T bp = T(b);
  utils::Matrix3<T> r(N1, N2, N3, 0);
  for (int i=0; i!=N1; ++i)
    for (int j=0; j!=N2; ++j)
      for (int k=0; k!=N3; ++k)
	r.elem(i,j,k) = a.elem(i,j,k) / bp;
  return r;
}




// Handle printing of Matrix3<T> objects:
template <typename T, typename U>
U & utils::operator<<(U & os, const utils::Matrix3<T> & sv){
  int N1 = sv.n1();
  int N2 = sv.n2();
  int N3 = sv.n3();

  for (int i=0; i!=N1; ++i){
    for (int j=0; j!=N2; ++j){
      for (int k=0; k!=N3; ++k){
	os << " " << sv.elem(i,j,k);
	os.clear();
      }
      os << std::endl;
    }
    os << std::endl;
  }
  //os << " ";
  return os;
}


#endif





