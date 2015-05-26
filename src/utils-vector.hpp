


#ifndef UTILS_VECTOR_HPP
#define UTILS_VECTOR_HPP



#include <string>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include <cstdlib>
#include <cmath>


#include "utils.hpp"
#include "utils-string.hpp"
#include "utils-errors.hpp"

#include "utils-vector3.hpp"




using std::cout;
using std::endl;

namespace utils {





  /* #########################################################################
     Declaration of vec<T> template class:
     #########################################################################
  */


  template <typename T>
  class Vector {

  private:
    // Memory pool:
    T *mvec;
    // Total number of elements visible to user:
    int mN;

    // Total capacity:
    int mNcap;
    // [amount of space to allocate when allocation needed]
    // = (1 + mresfrac) * [amount of space used]
    // Default: mresfrac = 1
    // NB: Negative values disallowed
    double mresfrac;


  public:
    Vector();
    ~Vector();

    Vector(const int N);
    Vector(const int N, const T & p);
    Vector(const Vector3<T> & sv); // conversion of Vector3 to Vector
    Vector(const T * p, int N); // Use array e.g. x, as input.

    Vector(const Vector & sv);

    Vector & operator=(const Vector & sv);
    Vector & operator=(const Vector3<T> & sv);

    // Operators:
    inline       T & operator[](const int & i);
    inline const T & operator[](const int & i) const;
    void to_array(T * p);        // explicit conversion of Vector to array
    Vector3<T> to_Vector3(void); // explicit conversion of Vector to Vector3


    inline int size() const;      // Get the number of explicitly stored elements.
    void resize(int N);           // Change the number of explicitly stored elements.
    void push_back(const T & p);   // Add an element to the end.

    // Vector operations:
    T normalize();
    T    magn() const;


    //typedef int size_t;
    //typedef int size_type;
    typedef T * iterator;
    typedef const T * const_iterator;

    T * end() const;
    T * begin() const;

    void   reserve_fraction(double rf);
    double reserve_fraction() const;

    int cap() const;
    void cap(int N);

    void trim();

  } ;


  // -----------------------------------
  // Nonmembers:
  // -----------------------------------

  // Handle a + b:
  template <typename T>
  Vector<T> operator+(const Vector<T> & a, const Vector<T> & b);

  // Handle a - b:
  template <typename T>
  Vector<T> operator-(const Vector<T> & a, const Vector<T> & b);

  // Handle a * b:
  template <typename T>
  T operator*(const Vector<T> & a, const Vector<T> & b);


  // Handle a * B, B being nontype:
  template <typename S, typename T>
  Vector<T> operator*(const Vector<T> & a, const S & B);

  // Handle A * b, A being nontype:
  template <typename S, typename T>
  Vector<T> operator*(const S & A, const Vector<T> & b);

  // Handle a / B, B being nontype:
  template <typename S, typename T>
  Vector<T> operator/(const Vector<T> & a, const S & B);


  // Handle printing of Vector<T> objects:
  template <typename T, typename U>
  U & operator << (U & os, const Vector<T> & sv);


  // Handle a * b, scalar product:
  template <typename T>
  T scalarproduct(const Vector<T> & a, const Vector<T> & b);

  // Handle a x b, vector product:
  template <typename T>
  void vectorproduct(const Vector<T> & a, const Vector<T> & b, Vector<T> & c);

  
}




/* #########################################################################
   Definitions:
   #########################################################################
*/




// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Default ctor:
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template <typename T>
utils::Vector<T>::Vector(){
  // No allocation:
  mvec = 0;
  mN   = 0;
  mresfrac = 1.0;
  mNcap = mN + mresfrac * mN;
}



// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Default dtor:
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template <typename T>
utils::Vector<T>::~Vector(){
  if (mvec != 0) delete [] mvec;
  mvec = 0;
  mN   = 0;
  mresfrac = 1.0;
  mNcap = mN + mresfrac * mN;
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Alternative ctors:
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename T>
utils::Vector<T>::Vector(const int N){
  mN   = ((N < 0) ? 0 : N);
  mvec = 0;
  mresfrac = 1.0;
  mNcap = mN + mresfrac * mN;
  if (mNcap>0)
    mvec = new T [mNcap] ();
}



template <typename T>
utils::Vector<T>::Vector(const int N, const T & p){
  mN   = ((N < 0) ? 0 : N);
  mvec = 0;
  mresfrac = 1.0;
  mNcap = mN + mresfrac * mN;
  if (mNcap>0){
    mvec = new T [mNcap] ();
    for (int i=0; i<mN; ++i) mvec[i] = p;
  }
}


template <typename T>
utils::Vector<T>::Vector(const Vector3<T> & sv){
  mN = 3;
  mvec = 0;
  mresfrac = 1.0;
  mNcap = mN + mresfrac * mN;
  mvec = new T [mNcap] ();
  mvec[0] = sv[0];
  mvec[1] = sv[1];
  mvec[2] = sv[2];
}



template <typename T>
utils::Vector<T>::Vector(const T * p, int N){
  mN = N < 0 ? -N : N;
  mvec = 0;
  mresfrac = 1.0;
  mNcap = mN + mresfrac * mN;
  mvec = new T [mNcap] ();
  try {
    for (int i=0; i<mN; ++i) mvec[i] = p[i];
  }
  catch (std::runtime_error & re){
    cout << "Error: Cannot construct Vector from length-deficient array." << endl;
    exit(EXIT_FAILURE);
  }
}




template <typename T>
void utils::Vector<T>::reserve_fraction(double rf){
  mresfrac = rf;
}

template <typename T>
double utils::Vector<T>::reserve_fraction() const {
  return mresfrac;
}

template <typename T>
int utils::Vector<T>::cap() const {
  return mNcap;
}


// Force a certain capacity, disregard 'reserve fraction' setting:
template <typename T>
void utils::Vector<T>::cap(int N){
  if (N<0) return;

  if (N==0){
    //cout << "requsted cap is 0" << endl;
    if (mvec != 0) delete [] mvec;
    mvec = 0;
    mN = mNcap = 0;
    mresfrac = 1;
    return;
  }
  
  if (N == mNcap){
    //cout << "[equality] cap is " << mNcap << " requested cap is " << N << endl;
    return;
  }

  if (N > mNcap){
    //cout << "[req. more] cap is " << mNcap << " requested cap is " << N << endl;

    if (mNcap == 0){
      //cout << "[req. more] allocing new vector with size " << N << endl;
      mNcap = N;
      mvec = new T [mNcap] ();
    }
    else {


      T *bak = new T [N] ();
      for (int i=0; i<mN; ++i) bak[i] = mvec[i];
      delete [] mvec;
      
      mvec = new T [N] ();
      for (int i=0; i<mN; ++i) mvec[i] = bak[i];
      delete [] bak;
      
      mNcap = N;
    }

    return;
  }
  

  if (N < mNcap){
    // shrink vector
      
    if (mN >= N) mN = N;
      
    T *bak = new T [N] ();
    for (int i=0; i<mN; ++i) bak[i] = mvec[i];
    delete [] mvec;

    mvec = new T [N] ();
    for (int i=0; i<mN; ++i) mvec[i] = bak[i];
    delete [] bak;

    mNcap = N;
      
    return;
  }
}







// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Copy constructor:
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template <typename T>
utils::Vector<T>::Vector(const utils::Vector<T> & sv){

  if (sv.mvec != 0){
    // Source 'sv' is allocated. Allocate destination, with same attributes.
    mN       = sv.mN;
    mresfrac = sv.mresfrac;
    mNcap = mN + mresfrac * mN;
    mvec = new T [mNcap] ();
    for (int i=0; i<((sv.mNcap < mNcap) ? sv.mNcap : mNcap); ++i)
      mvec[i] = sv.mvec[i];
  }
  else {
    mvec = 0;
    mN   = 0;
    mresfrac = 1.0;
    mNcap = mN + mresfrac * mN;
  }
}




// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Operators:
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename T>
utils::Vector<T> & utils::Vector<T>::operator=(const utils::Vector<T> & sv){

  if (this==&sv) return *this;


  if (sv.mvec == 0){
    // Source argument is empty
    // Delete destination, if it is allocated
    if (mvec != 0) delete[] mvec;
    // Make destination empty
    mvec = 0;
    mN   = 0;
    mresfrac = 1.0;
    mNcap = mN + mresfrac * mN;
    return *this;
  }


  // Source is NOT empty

  // Delete destination, if it is allocated
  if (mvec != 0) delete[] mvec;

  mvec = 0;
  mN   = 0;
  mresfrac = 1.0;
  mNcap = mN + mresfrac * mN;
  

  // Allocate destination and copy
  mN       = sv.mN;
  mresfrac = sv.mresfrac;
  mNcap    = sv.mNcap;
  mvec = new T [mNcap] ();

  for (int i=0; i<mN; ++i){
    //std::cout << "i=" << i << std::endl; std::cout.flush();
    mvec[i] = sv.mvec[i]; // no problem to access private member mvec in other object sv here ... ?
  }

  return *this;
}





template <typename T>
utils::Vector<T> & utils::Vector<T>::operator=(const utils::Vector3<T> & sv){

  //  if (this==&sv) return *this;

  // Delete destination, if it is allocated
  if (mvec != 0) delete[] mvec;
  mvec = 0;

  // Allocate destination and copy
  mN       = 3;
  mresfrac = 1.0;
  mNcap = mN + mresfrac * mN;
  mvec = new T [mNcap] ();

  mvec[0] = sv[0];
  mvec[1] = sv[1];
  mvec[2] = sv[2];

  return *this;
}






template <typename T>
T & utils::Vector<T>::operator[](const int & i) {

  /*
  if (i<0 || i>=mN){
    std::cout << "Vector index " << i << " is out of range. Allowed: " << 0 << " to " << mN-1 << ". Exiting." << std::endl;
    //std::cout << "Error in " << __FILE__ << " at line number " << __LINE__ << std::endl;
    exit(EXIT_FAILURE);
  }
  if (mvec==0){
    std::cout << "Vector is not allocated. Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
  */
  return mvec[i];
}

template <typename T>
const T & utils::Vector<T>::operator[](const int & i) const {

  /*
  if (i<0 || i>=mN){
    std::cout << "Vector index " << i << " is out of range. Allowed: " << 0 << " to " << mN-1 << ". Exiting." << std::endl;
    //std::cout << "Error in " << __FILE__ << " at line number " << __LINE__ << std::endl;
    exit(EXIT_FAILURE);
  }
  if (mvec==0){
    std::cout << "Vector is not allocated. Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
  */
  return mvec[i];
}



template <typename T>
void utils::Vector<T>::to_array(T * p){
  try {
    for (int i=0; i<mN; ++i)
      p[i] = mvec[i];
  }
  catch (std::runtime_error & re){
    cout << "Error: Cannot copy Vector(" << mN << ") to length-deficient array." << endl;
    exit(EXIT_FAILURE);
  }
}


template <typename T>
utils::Vector3<T> utils::Vector<T>::to_Vector3(void){
  if (size()<3){
    cout << "Error: Cannot construct Vector3 from length-deficient array." << endl;
    exit(EXIT_FAILURE);
  }
  Vector3<T> r;
  r[0] = mvec[0];
  r[1] = mvec[1];
  r[2] = mvec[2];
  return r;
}






// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// For storage changes/updates/queries:
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++









// Get the number of explicitly stored elements.
template <typename T>
int utils::Vector<T>::size() const {
  return mN;
}



// Change the number of visible elements:
template <typename T>
void utils::Vector<T>::resize(const int N){
  if (N==0){
    cap(0);
    mN = 0;
  }
  else {
    cap(N + mresfrac * N);
    mN = N;
  }
}


template <typename T>
void utils::Vector<T>::trim(){
  cap(mN);
}




// Add an element to the end.
template <typename T>
void utils::Vector<T>::push_back(const T & p){
  if (mN+1 > mNcap){
    //cout << "push_back: allocating more space to vector, now " << mNcap << " will be " << mN+1 + mresfrac * (mN+1) << endl;
    cap(mN+1 + mresfrac * (mN+1));
  }
  //cout << "push_back: mN-1 is " << mN-1 << endl;
  mvec[mN] = p;
  mN++;
}





// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Vector operations:
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename T>
T utils::Vector<T>::normalize() {
  T tmp = 0, len = 0;
  for (int i=0; i<mN; ++i) tmp = tmp + mvec[i]*mvec[i];
  len = sqrt(tmp);
  tmp = 1.0/len;
  for (int i=0; i<mN; ++i) mvec[i] = mvec[i] * tmp;
  return len;
}


template <typename T>
T utils::Vector<T>::magn() const {
  T tmp1 = 0;
  for (int i=0; i<mN; ++i){
    tmp1 = tmp1 + mvec[i] * mvec[i];
  }
  return sqrt(tmp1);
}





template <typename T>
T * utils::Vector<T>::begin() const {
  if (size()>0)
    return &(mvec[0]);
  else
    return 0;
}

template <typename T>
T * utils::Vector<T>::end() const {
  if (size()>0)
    return &(mvec[size()-1]);
  else
    return 0;
}





#if 0
// Handle Vector3 a = Vector b:

template <typename T>
utils::Vector3<T> & utils::operator=(utils::Vector3<T> & a, const utils::Vector<T> & b){
  if (b.size()>=3){
    a[0]=b[0];
    a[1]=b[1];
    a[2]=b[2];
  }
  else {
    cout << "Error: In assignment Vector3 a = Vector b the size of b is less than 3." << endl;
    exit(EXIT_FAILURE);
  }
  return a;
}
#endif




// ------------------------------------------------------------
// Non-members
// ------------------------------------------------------------

// Handle a + b:
template <typename T>
utils::Vector<T> utils::operator+(const utils::Vector<T> & a, const utils::Vector<T> & b){
  if (a.size()!=b.size()){
    std::cout << "Cannot add vectors of different sizes. Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
  Vector<T> r(a.size(), 0);
  for (int i=0; i<a.size(); ++i) r[i] = a[i] + b[i];
  return r;
}

// Handle a - b:
template <typename T>
utils::Vector<T> utils::operator-(const utils::Vector<T> & a, const utils::Vector<T> & b){
  if (a.size()!=b.size()){
    std::cout << "Cannot subtract vectors of different sizes. Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
  Vector<T> r(a.size(), 0);
  for (int i=0; i<a.size(); ++i) r[i] = a[i] - b[i];
  return r;
}


// Handle a * b:
template <typename T>
T utils::operator*(const utils::Vector<T> & a, const utils::Vector<T> & b){
  if (a.size()!=b.size()){
    std::cout << "Cannot multiply vectors of different sizes. Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
  T tmp=0;
  for (int i=0; i<a.size(); ++i) tmp += a[i] * b[i];
  return tmp;
}




// Handle a * B, B being nontype:
template <typename S, typename T>
utils::Vector<T> utils::operator*(const utils::Vector<T> & a, const S & b){
  Vector<T> r(a.size(), 0);
  for (int i=0; i<a.size(); ++i) r[i] = a[i] * T(b);
  return r;
}

// Handle A * b, A being nontype:
template <typename S, typename T>
utils::Vector<T> utils::operator*(const S & a, const utils::Vector<T> & b){
  Vector<T> r(b.size(), 0);
  for (int i=0; i<b.size(); ++i) r[i] = T(a) * b[i];
  return r;
}

// Handle a / B, B being nontype:
template <typename S, typename T>
utils::Vector<T> utils::operator/(const utils::Vector<T> & a, const S & b){
  Vector<T> r(a.size(), 0);
  for (int i=0; i<a.size(); ++i) r[i] = a[i] / T(b);
  return r;
}




// Handle printing of a:
template <typename T, typename U>
U & utils::operator << (U & os, const utils::Vector<T> & sv){
  for (int i=0; i<sv.size(); ++i){
    os << " " << sv[i];
    os.clear();
  }
  os << " ";
  return os;
}



// Handle a * b, scalar product:
template <typename T>
T utils::scalarproduct(const utils::Vector<T> & a, const utils::Vector<T> & b){
  return utils::operator*(a, b);
}




// Handle a x b, vector product:
template <typename T>
void utils::vectorproduct(const utils::Vector<T> & a,
			  const utils::Vector<T> & b,
			  utils::Vector<T> & c){
  if (a.size() != b.size()){
    utils::aborterror("Cannot vector multiply vectors of different dimensions: "
		      + utils::tostring(a.size()) + " and "
		      + utils::tostring(b.size()) );
  }
  if (c.size() != a.size()) c.resize(a.size());

  c[2] = a[0] * b[1] - a[1] * b[0];
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
}






/*
// Handle a * b, vector product:
template <typename T>
utils::vec<T> utils::vectorproduct(const utils::vec<T> & a, const utils::vec<T> & b){
  if (a.size()!=b.size()){
    std::cout << "Cannot multiply vectors of different sizes. Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
  utils::vec<T> res(a.size(), 0);

  for (int i=0; i!=a.size(); ++i){
    for (int j=0; j!=a.size(); ++j){



    }
  }


  return res;
}
*/




#endif

