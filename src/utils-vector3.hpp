


#ifndef UTILS_VECTOR3_HPP
#define UTILS_VECTOR3_HPP



#include <string>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include <cstdlib>
#include <cmath>



#include "utils.hpp"
#include "utils-string.hpp"
#include "utils-errors.hpp"



using std::cout;
using std::endl;




namespace utils {

  template <typename T>
  class Vector3 {

  private:
    // Memory pool:
    T mvec[3];


  public:
    Vector3();
    ~Vector3();
    
    Vector3(const T & p);
    Vector3(const T * p, int N); // Use array e.g. x as input.
    Vector3(T & p1, T & p2, T & p3);
    Vector3(const T & p1, const T & p2, const T & p3);
    
    Vector3(const Vector3 & sv);

    //Vector3 & operator=(const Vector<T> & sv);
    Vector3 & operator=(const Vector3   & sv);

    // Operators:
    inline       T & operator[](const int & i);
    inline const T & operator[](const int & i) const;
    void to_array(T * p);        // explicit conversion of Vector3 to array


    Vector3<T> * This(){ return this; }

    // Vector3 operations:
    T normalize();
    T    magn() const;

    typedef T * iterator;
    typedef const T * const_iterator;

    T * end() const;
    T * begin() const;

  } ;



  // Nonmembers:

  // Handle a * b, scalar product:
  template <typename T>
  T scalarproduct(const Vector3<T> & a, const Vector3<T> & b);

  // Handle a x b, vector product:
  template <typename T>
  void vectorproduct(const Vector3<T> & a, const Vector3<T> & b, Vector3<T> & c);


  // UNARY OPERATORS

  // Handle &a:
  template <typename T>
  Vector3<T> * operator&(Vector3<T> & a);
  // Handle +a:
  template <typename T>
  Vector3<T> & operator+(const Vector3<T> & a);
  // Handle -a:
  template <typename T>
  Vector3<T> operator-(const Vector3<T> & a);

  // Handle ++a:
  template <typename T>
  Vector3<T> operator++(Vector3<T> & a);
  // Handle --a:
  template <typename T>
  Vector3<T> operator--(Vector3<T> & a);
  // Handle a++:
  template <typename T>
  Vector3<T> operator++(Vector3<T> & a, int );
  // Handle a--:
  template <typename T>
  Vector3<T> operator--(Vector3<T> & a, int );


  // BINARY OPERATORS

  // Handle +=:
  template <typename T>
  Vector3<T> & operator+=(Vector3<T> & a, const Vector3<T> & b);
  // Handle -=:
  template <typename T>
  Vector3<T> & operator-=(Vector3<T> & a, const Vector3<T> & b);

  // Handle a *= B, B other type:
  template <typename S, typename T>
  Vector3<T> & operator*=(Vector3<T> & a, const S & b);
  // Handle a /= B, B other type:
  template <typename S, typename T>
  Vector3<T> & operator/=(Vector3<T> & a, const S & b);

  // Handle a < b:
  template <typename T>
  bool operator<(const Vector3<T> & a, const Vector3<T> & b);
  // Handle a > b:
  template <typename T>
  bool operator>(const Vector3<T> & a, const Vector3<T> & b);
  // Handle a == b:
  template <typename T>
  bool operator==(const Vector3<T> & a, const Vector3<T> & b);

  // Handle a + b:
  template <typename T>
  Vector3<T> operator+(const Vector3<T> & a, const Vector3<T> & b);

  // Handle a - b:
  template <typename T>
  Vector3<T> operator-(const Vector3<T> & a, const Vector3<T> & b);

  // Handle a * b:
  template <typename T>
  T operator*(const Vector3<T> & a, const Vector3<T> & b);

  // Handle a * b, b being nontype:
  template <typename S, typename T>
  Vector3<T> operator*(const Vector3<T> & a, const S & b);

  // Handle a * b, a being nontype:
  template <typename S, typename T>
  Vector3<T> operator*(const S & a, const Vector3<T> & b);

  // Handle a / b, b being nontype:
  template <typename S, typename T>
  Vector3<T> operator/(const Vector3<T> & a, const S & b);

  // Handle printing of Vector3<T> objects:
  template <typename T, typename U>
  U & operator << (U & os, const Vector3<T> & sv);

  
}




/* #########################################################################
   Definitions:
   #########################################################################
*/




// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Default ctor:
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template <typename T>
utils::Vector3<T>::Vector3(){
  mvec[0] = mvec[1] = mvec[2] = T(0);
}



// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Default dtor:
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template <typename T>
utils::Vector3<T>::~Vector3(){
  mvec[0] = mvec[1] = mvec[2] = T(0);
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Alternative ctors:
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename T>
utils::Vector3<T>::Vector3(const T & p){
  mvec[0] = mvec[1] = mvec[2] = p;
}


template <typename T>
utils::Vector3<T>::Vector3(const T * p, int N){
  try {
    mvec[0] = p[0];
    mvec[1] = p[1];
    mvec[2] = p[2];
  }
  catch (std::runtime_error & re){
    cout << "Error: Cannot construct Vector3 from length-deficient array." << endl;
    exit(EXIT_FAILURE);
  }
}



template <typename T>
utils::Vector3<T>::Vector3(T & p1, T & p2, T & p3){
  mvec[0] = p1;
  mvec[1] = p2;
  mvec[2] = p3;
}

template <typename T>
utils::Vector3<T>::Vector3(const T & p1, const T & p2, const T & p3){
  mvec[0] = p1;
  mvec[1] = p2;
  mvec[2] = p3;
}



// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Copy constructor:
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template <typename T>
utils::Vector3<T>::Vector3(const utils::Vector3<T> & sv){
  mvec[0] = sv.mvec[0];
  mvec[1] = sv.mvec[1];
  mvec[2] = sv.mvec[2];
}





// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Operators:
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename T>
utils::Vector3<T> & utils::Vector3<T>::operator=(const utils::Vector3<T> & sv){
  if (this==&sv) return *this;
  
  mvec[0] = sv.mvec[0];
  mvec[1] = sv.mvec[1];
  mvec[2] = sv.mvec[2];

  return *this;
}


/*
template <typename T>
utils::Vector3<T> & utils::Vector3<T>::operator=(const utils::Vector<T> & sv){
  if (sv.size()<3){
    cout << "Error: Cannot assign Vector to Vector3 when Vector has size less than 3." << endl;
    exit(EXIT_FAILURE);
  }
  mvec[0] = sv[0];
  mvec[1] = sv[1];
  mvec[2] = sv[2];

  return *this;
}
*/





template <typename T>
T & utils::Vector3<T>::operator[](const int & i) {
  return mvec[i];
}

template <typename T>
const T & utils::Vector3<T>::operator[](const int & i) const {
  return mvec[i];
}



template <typename T>
void utils::Vector3<T>::to_array(T * p){
  try {
    p[0] = mvec[0];
    p[1] = mvec[1];
    p[2] = mvec[2];
  }
  catch (std::runtime_error & re){
    cout << "Error: Cannot copy Vector3 to length-deficient array." << endl;
    exit(EXIT_FAILURE);
  }
}





// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Vector3 operations:
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename T>
T utils::Vector3<T>::normalize() {
  T tmp = mvec[0]*mvec[0] + mvec[1]*mvec[1] + mvec[2]*mvec[2];
  T len = sqrt(tmp);
  tmp = 1.0/len;
  mvec[0] *= tmp;
  mvec[1] *= tmp;
  mvec[2] *= tmp;
  return len;
}


template <typename T>
T utils::Vector3<T>::magn() const {
  T tmp = mvec[0]*mvec[0] + mvec[1]*mvec[1] + mvec[2]*mvec[2];
  return sqrt(tmp);
}


template <typename T>
T * utils::Vector3<T>::begin() const {
  return &(mvec[0]);
}

template <typename T>
T * utils::Vector3<T>::end() const {
  return &(mvec[2]);
}




// Non-members


// Handle a * b, scalar product:
template <typename T>
T utils::scalarproduct(const utils::Vector3<T> & a, const utils::Vector3<T> & b){
  return utils::operator*(a, b);
}


// Handle a x b, vector product:
template <typename T>
void utils::vectorproduct(const utils::Vector3<T> & a,
			  const utils::Vector3<T> & b,
			  utils::Vector3<T> & c){
  c[2] = a[0] * b[1] - a[1] * b[0];
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
}



// UNARY OPERATORS

// Handle &a:
template <typename T>
utils::Vector3<T> * utils::operator&(utils::Vector3<T> & a){
  return a.This();
}
// Handle +a:
template <typename T>
utils::Vector3<T> & utils::operator+(const utils::Vector3<T> & a){
  return a;
}
// Handle -a:
template <typename T>
utils::Vector3<T> utils::operator-(const utils::Vector3<T> & a){
  Vector3<T> r;
  r[0] = - a[0];
  r[1] = - a[1];
  r[2] = - a[2];
  return r;
}
// Handle ++a:
template <typename T>
utils::Vector3<T> utils::operator++(utils::Vector3<T> & a){
  Vector3<T> r;
  r[0] = ++a[0];
  r[1] = ++a[1];
  r[2] = ++a[2];
  return r;
}
// Handle --a:
template <typename T>
utils::Vector3<T> utils::operator--(utils::Vector3<T> & a){
  Vector3<T> r;
  r[0] = --a[0];
  r[1] = --a[1];
  r[2] = --a[2];
  return r;
}
// Handle a++:
template <typename T>
utils::Vector3<T> utils::operator++(utils::Vector3<T> & a, int ){
  Vector3<T> r;
  r[0] = a[0]++;
  r[1] = a[1]++;
  r[2] = a[2]++;
  return r;
}
// Handle a--:
template <typename T>
utils::Vector3<T> utils::operator--(utils::Vector3<T> & a, int ){
  Vector3<T> r;
  r[0] = a[0]--;
  r[1] = a[1]--;
  r[2] = a[2]--;
  return r;
}


// BINARY OPERATORS

// Handle a += b:
template <typename T>
utils::Vector3<T> & utils::operator+=(utils::Vector3<T> & a, const utils::Vector3<T> & b){
  a[0] += b[0];
  a[1] += b[1];
  a[2] += b[2];
  return a;
}
// Handle a -= b:
template <typename T>
utils::Vector3<T> & utils::operator-=(utils::Vector3<T> & a, const utils::Vector3<T> & b){
  a[0] -= b[0];
  a[1] -= b[1];
  a[2] -= b[2];
  return a;
}
// Handle a *= B, B other type:
template <typename S, typename T>
utils::Vector3<T> & utils::operator*=(utils::Vector3<T> & a, const S & b){
  T tmp = T(b);
  a[0] *= tmp;
  a[1] *= tmp;
  a[2] *= tmp;
  return a;
}
// Handle a /= B, B other type:
template <typename S, typename T>
utils::Vector3<T> & utils::operator/=(utils::Vector3<T> & a, const S & b){
  T tmp = T(1)/T(b);
  a[0] *= tmp;
  a[1] *= tmp;
  a[2] *= tmp;
  return a;
}


// Handle a < b:
template <typename T>
bool utils::operator<(const utils::Vector3<T> & a, const utils::Vector3<T> & b){
  int n=0;
  if (a[0] < b[0]) n++;
  if (a[1] < b[1]) n++;
  if (a[2] < b[2]) n++;
  if (n==3) return true;
  else return false;
}
// Handle a > b:
template <typename T>
bool utils::operator>(const utils::Vector3<T> & a, const utils::Vector3<T> & b){
  return operator<(b,a);
}
// Handle a == b:
template <typename T>
bool utils::operator==(const utils::Vector3<T> & a, const utils::Vector3<T> & b){
  int n=0;
  if (a[0] == b[0]) n++;
  if (a[1] == b[1]) n++;
  if (a[2] == b[2]) n++;
  if (n==3) return true;
  else return false;
}


// Handle a + b:
template <typename T>
utils::Vector3<T> utils::operator+(const utils::Vector3<T> & a, const utils::Vector3<T> & b){
  Vector3<T> r;
  r[0] = a[0] + b[0];
  r[1] = a[1] + b[1];
  r[2] = a[2] + b[2];
  return r;
}

// Handle a - b:
template <typename T>
utils::Vector3<T> utils::operator-(const utils::Vector3<T> & a, const utils::Vector3<T> & b){
  Vector3<T> r;
  r[0] = a[0] - b[0];
  r[1] = a[1] - b[1];
  r[2] = a[2] - b[2];
  return r;
}


// Handle a * b:
template <typename T>
T utils::operator*(const utils::Vector3<T> & a, const utils::Vector3<T> & b){
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}


// Handle a * b, b being nontype:
template <typename S, typename T>
utils::Vector3<T> utils::operator*(const utils::Vector3<T> & a, const S & b){
  Vector3<T> r;
  r[0] = a[0] * T(b);
  r[1] = a[1] * T(b);
  r[2] = a[2] * T(b);
  return r;
}

// Handle a * b, a being nontype:
template <typename S, typename T>
utils::Vector3<T> utils::operator*(const S & a, const utils::Vector3<T> & b){
  Vector3<T> r;
  r[0] = T(a) * b[0];
  r[1] = T(a) * b[1];
  r[2] = T(a) * b[2];
  return r;
}

// Handle a / b, b being nontype:
template <typename S, typename T>
utils::Vector3<T> utils::operator/(const utils::Vector3<T> & a, const S & b){
  Vector3<T> r;
  r[0] = a[0] / T(b);
  r[1] = a[1] / T(b);
  r[2] = a[2] / T(b);
  return r;
}


template <typename T, typename U>
U & utils::operator << (U & os, const utils::Vector3<T> & sv){
  os << " " << sv[0]
     << " " << sv[1]
     << " " << sv[2]
     << " ";
  os.clear();
  return os;
}







#endif

