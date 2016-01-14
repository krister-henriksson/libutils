


#ifndef NR_GOLDEN_HPP
#define NR_GOLDEN_HPP


#include <iostream>
#include <limits>

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "utils.hpp"
#include "utils-math.hpp"



///////////////////////////////////////////////////////////
// Usage:
///////////////////////////////////////////////////////////
//
//   F1dim< ChiSqFunc<Param, double, double> > f1d(fparams.X(), gradient, cs);
//
//   Golden gl(f1d);
//
//   xmin = gld.minimize(a, b, f1d);
//   fmin = gld.fmin;
//
///////////////////////////////////////////////////////////



namespace nr {


  template <typename T>
  class Golden
  {



    //#####################################################################
    // Member variables
    //#####################################################################
  public:
    T & func;
    double ax, bx, cx, fa, fb, fc;
    double xmin, fmin;


    //#####################################################################
    // Constructor
    //#####################################################################
    // NB:  Use constructor to set functor.
    Golden(T & funcc)
      : func(funcc)
    {
      /*
      std::cout << "Size of p  vector: " << func.p.size() << std::endl;
      std::cout << "Size of xi vector: " << func.xi.size() << std::endl;
      */
    }
    


    //#####################################################################
    // Main method 1
    //#####################################################################
    double minimize(double a,
		    double b,
		    double smin,
		    double smax,
		    double tol=-1,
		    bool debug=false,
		    bool debug_bracket=false
		    ){
      // Input: function/functor func and a bracketing triplet ax, bx, cx
      // such that bx is between ax and cx
      // Purpose: Locate minimum of func
      // Output: x value that minimizes func, members xmin and fmin
      // are also set.


      bracket(a, b, smin, smax, tol, debug_bracket);


      // Gets ax, bx, cx, fa, fb, fc
      
      const double R=0.61803399, C=1.0-R;
      double x1, x2;
      double x0=ax;
      double x3=cx;
      using utils::abs;
      int niter=0;


      tol = (tol > 0 ? tol : sqrt(utils::eps_d()));


      //std::cout.precision(std::numeric_limits<double>::digits10);
      //std::cout.precision(15);

      if (abs(cx-bx) > abs(bx-ax)){
	x1=bx;
	x2=bx+C*(cx-bx);
      }
      else {
	x2=bx;
	x1=bx-C*(bx-ax);
      }
      double f1=func(x1);
      double f2=func(x2);

      if (debug)
	std::cout << "golden: Iteration " << niter << ": "
		  << "x1 x2 " << x1 << " " << x2 << " "
		  << "f1 f2 " << f1 << " " << f2 << std::endl;

      // ####################################################################
      while (utils::abs(x3-x0) > tol*(utils::abs(x1)+utils::abs(x2))){
      //      while (utils::fp_are_equal(x1, x2) == false){
	niter++;

	if (debug)
	  std::cout << "golden: Iteration " << niter << ": "
		    << "x1 x2 " << x1 << " " << x2 << " "
		    << "f1 f2 " << f1 << " " << f2 << std::endl;
	

	if (utils::fp_are_equal(x1, x2))
	  break;
	

	if (f2 < f1){
	  shift3(x0, x1, x2, R*x2+C*x3);
	  shift2(f1, f2, func(x2));
	}
	else {
	  shift3(x3, x2, x1, R*x1+C*x0);
	  shift2(f2, f1, func(x1));
	}
      }
      // ####################################################################


      if (f1 < f2){
	xmin=x1;
	fmin=f1;
      }
      else {
	xmin=x2;
	fmin=f2;
      }
      return xmin;
    }





    //#####################################################################
    // Main method 2
    //#####################################################################
    void bracket(const double a,
		 const double b,
		 double smin,
		 double smax,
		 double tol=-1,
		 bool debug=false){
      // Input: distinct initial points a, b and function/functor func
      // Purpose: Search for minimum of func in the downhill direction
      // (determined from function values at points a, b).
      // Output: points ax, bx, cx that bracket the minimum of func
      // and the values of func at these points
      const double GOLD=1.618034, GLIMIT=100.0,
	TINY = std::numeric_limits<double>::epsilon();//utils::eps_d();
      int niter=0;
      double r, q, u, ulim;


      tol = (tol > 0 ? tol : sqrt(utils::eps_d()));


      ax = a;
      bx = b;
      double fu;
      //std::cout << "golden: bracket: made it here 1" << std::endl;
      fa = func(ax);
      //std::cout << "golden: bracket: made it here 2" << std::endl;
      fb = func(bx);
      if (fb > fa){
	utils::swap(fa, fb);
	utils::swap(ax, bx);
      }
      cx = bx + GOLD*(bx-ax);
      /*
      if (ax<bx){
	while (cx<smin || cx>smax) cx = bx + 0.5 * GOLD*(bx-ax);
      }
      else {
	while (cx<smin || cx>smax) cx = bx - 0.5 * GOLD*(bx-ax);
      }
      */
      fc = func(cx);

      std::cout.precision(15);
      
      if (debug)
	std::cout << "golden: bracket: Iteration " << niter << ": "
		  << "ax bx cx "  << ax << " " << bx << " " << cx << " "
		  << "fa fb fc " << fa << " " << fb << " " << fc
		  << std::endl;


      // ####################################################################
      while (fb > fc){

	niter++;
	if (debug)
	  std::cout << "golden: bracket: Iteration " << niter << ": "
 		    << "ax bx cx "  << ax << " " << bx << " " << cx << " "
		    << "fa fb fc " << fa << " " << fb << " " << fc
		    << std::endl;

	r = (bx-ax)*(fb-fc);
	q = (bx-cx)*(fb-fa);
	u = bx - ((bx-cx)*q-(bx-ax)*r)/
	  (2.0*utils::sign_nr(utils::max(utils::abs(q-r),TINY), q-r));
	ulim = bx + GLIMIT*(cx-bx);

	if ((bx-u)*(u-cx) > 0.0){
	  fu = func(u);
	  if (fu < fc){
	    ax=bx; bx=u;
	    fa=fb; fb=fu;
	    return;
	  }
	  else if (fu > fb){
	    cx=u;
	    fc=fu;
	    return;
	  }
	  u = cx + GOLD*(cx-bx);


	  fu = func(u);
	}
	else if ((cx-u)*(u-ulim) > 0.0){
	  fu = func(u);
	  if (fu < fc){
	    shift3(bx, cx, u, u+GOLD*(u-cx));
	    shift3(fb, fc, fu, func(u));
	  }
	}
	else if ((u-ulim)*(ulim-cx) >= 0.0){
	  u=ulim;
	  fu = func(u);
	}
	else {
	  u = cx + GOLD*(cx-bx);


	  fu = func(u);
	}
	shift3(ax, bx, cx, u);
	shift3(fa, fb, fc, fu);
      }
      // ####################################################################

    }


    inline void shift2(double & a, double & b, const double c){
      a=b; b=c;
    }
    inline void shift3(double & a, double & b, double & c, const double d){
      a=b; b=c; c=d;
    }
    inline void move3(double & a, double & b, double & c,
		      const double d, const double e, const double f){
      a=d; b=e; c=f;
    }





  };

}



#endif

