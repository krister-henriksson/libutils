


#ifndef UTILS_MATH_HPP
#define UTILS_MATH_HPP



#include <limits>

#include <cfloat>
#include <cmath>




namespace utils {



  // #####################################################################
  // Machine accuracy values
  // #####################################################################
  inline float       eps_f(void){  return FLT_EPSILON; }
  inline double      eps_d(void){  return DBL_EPSILON; }
  inline long double eps_ld(void){ return LDBL_EPSILON; }


  // #####################################################################
  // Calculate n = x/y, where n is an integer.
  // #####################################################################
  inline int fp_divide_to_integer(double x, double y){
    int m = x/y;
    double r1, r2, r3, td;

    td = x - m*y;     if (td<0) r1=-td; else r1=td;
    td = x - (m-1)*y; if (td<0) r2=-td; else r2=td;
    td = x - (m+1)*y; if (td<0) r3=-td; else r3=td;

    if (r1<=r2 && r1<=r3) return m;
    if (r2<=r1 && r2<=r3) return m-1;
    if (r3<=r1 && r3<=r2) return m+1;
  }

  inline double fp_divide(double x, double y){
    int s = ((y<0) ? -1 : 1);
    
    if (s*y < DBL_EPSILON) return x;
    else                   return x/y;
  }




  // #####################################################################
  // Check if two floating point values are equal (paranoid version)
  // #####################################################################
  template <typename T>
  bool fp_are_equal(const T & x, const T & y, const T & tol=-1){
    T eps = std::numeric_limits<T>::epsilon();
    T absx = (x > 0 ? x : -x);
    T absy = (y > 0 ? y : -y);

    eps = (tol > 0 ? tol : eps);

    if (absx <= eps && absy <= eps) return true;

    if (x>=0 && y>=0 &&
	( ( x <= y*(1+eps) && x >= y*(1-eps) ) ||
	  ( y <= x*(1+eps) && y >= x*(1-eps) ) ) ) return true;

    if (x<=0 && y<=0 &&
	( ( absx <= absy*(1+eps) && absx >= absy*(1-eps) ) ||
	  ( absy <= absx*(1+eps) && absy >= absx*(1-eps) ) ) ) return true;

    return false;

    /*
    if ( ((absx <= eps)                      && (absy <= eps))                      ||
	 ((y - absy * eps <= x)              && (x              <= y + absy * eps)) ||
	 ((y - absy * eps <= x - absx * eps) && (x - absx * eps <= y + absy * eps)) ||
	 ((y - absy * eps <= x + absx * eps) && (x + absx * eps <= y + absy * eps)) ||
	 ((x - absx * eps <= y)              && (y              <= x + absx * eps)) ||
	 ((x - absx * eps <= y - absy * eps) && (y - absy * eps <= x + absx * eps)) ||
	 ((x - absx * eps <= y + absy * eps) && (y + absy * eps <= x + absx * eps)) )
      return true;
    else
      return false;
    */
  }


  /*
  template <typename T>
  bool fp_are_equal_within_tol(const T & x, const T & y, const T & tol){
    T absx = (x > 0 ? x : -x);
    T absy = (y > 0 ? y : -y);

    // one possible version: |x-y|/(0.5*(|x|+|y|) < tol

    if ( (2*x - absx * tol <= 2*y + absy * tol) &&
	 (2*x + absx * tol >= 2*y - absy * tol) )
      return true;
    else
      return false;
  }
  */




  // #####################################################################
  // Check for small floating point value
  // #####################################################################
  template <typename T>
  bool fp_is_small(const T & x){
    T eps = std::numeric_limits<T>::epsilon();
    T absx = (x > 0 ? x : -x);

    //eps = sqrt(eps);

    if (absx < eps) return true;
    else            return false;
  }



  template <typename T>
  bool fp_is_small_tol(const T & x, const T eps){
    T absx = (x > 0 ? x : -x);

    if (absx < eps) return true;
    else            return false;
  }




}

#endif



