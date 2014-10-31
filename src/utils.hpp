


#ifndef UTILS_HPP
#define UTILS_HPP



namespace utils {



  template <typename T>
  inline
  T max(const T & x, const T & y){
    return x > y ? x : y;
  }

  template <typename T>
  inline
  T min(const T & x, const T & y){
    return x < y ? x : y;
  }

  template <typename T>
  inline
  void swap(T & x, T & y){
    T tmp = x;
    x = y;
    y = tmp;
  }


  template <typename T>
  T min3(const T & a, const T & b, const T & c){
    return ((a < b ? a : b) < c ? (a < b ? a : b) : c);
  }

  template <typename T>
  T max3(const T & a, const T & b, const T & c){
    return ((a > b ? a : b) > c ? (a > b ? a : b) : c);
  }


  template <typename T>
  inline
  T absval(const T & x){
    return (x > T(0)  ? x : -x);
  }

  template <typename T>
  inline
  T abs(const T & x){
    return (x > T(0)  ? x : -x);
  }




  template <typename T>
  inline
  int sign(const T & x){
    return (x > T(0) ? 1 : -1);
  }

  template <typename T>
  inline
  int sign_nr(const T & x, const T & y){
    return (y > T(0) ? absval(x) : -absval(x));
  }


  template <typename T>
  inline
  T square(const T & x){
    return x*x;
  }




  template <typename T>
  void copy_left_2args(T & a, T & b, bool pbc=false){
    if (pbc==false)  a = b;
    else { T tmp(a); a = b; b = tmp; }
  }

  template <typename T>
  void copy_left_3args(T & a, T & b, T & c, bool pbc=false){
    if (pbc==false){ a = b; b = c; }
    else { T tmp(a); a = b; b = c; c = tmp;
    }
  }

  template <typename T>
  void copy_left_4args(T & a, T & b, T & c, T & d, bool pbc=false){
    if (pbc==false){ a = b; b = c; c = d; }
    else { T tmp(a); a = b; b = c; c = d; d = tmp; }
  }




  template <typename T>
  void copy_right_2args(T & a, T & b, bool pbc=false){
    if (pbc==false)  b = a;
    else { T tmp(b); b = a; a = tmp; }
  }

  template <typename T>
  void copy_right_3args(T & a, T & b, T & c, bool pbc=false){
    if (pbc==false){ c = b; b = a; }
    else { T tmp(c); c = b; b = a; a = tmp;
    }
  }

  template <typename T>
  void copy_right_4args(T & a, T & b, T & c, T & d, bool pbc=false){
    if (pbc==false){ d = c; c = b; b = a; }
    else { T tmp(d); d = c; c = b; b = a; a = tmp; }
  }





}

#endif


