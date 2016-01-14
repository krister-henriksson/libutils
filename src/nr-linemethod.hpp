


#ifndef NR_LINEMETHOD_HPP
#define NR_LINEMETHOD_HPP


#include <iostream>
#include <cstdlib>

#include <limits>

#include "utils.hpp"
#include "utils-math.hpp"
#include "utils-errors.hpp"

#include "nr-golden.hpp"
#include "nr-f1dim.hpp"

#include "param.hpp"

#include "funcfit-basics.hpp"
#include "funcfit-errors.hpp"


using std::cout;
using std::endl;
using std::numeric_limits;
using utils::Vector;
using utils::fp_is_small;
using utils::fp_are_equal;



namespace nr {

  template <typename T>
  class Linemethod
  {
    // - used primarily by classes which need to minimize a function
    // in some direction


    //#####################################################################
    // Member variables
    //#####################################################################
  public:
    Vector<double> p;
    Vector<double> xi;
    T & func;
    funcfit::Minimization_Status status;



    //#####################################################################
    // Constructor
    //#####################################################################
    Linemethod(T & funcc)
      :
      func(funcc)
    {}


    //#####################################################################
    // Main method
    //#####################################################################
    double linmin(Vector<double> & point_io,
		  Vector<double>        & xmin,
		  Vector<double>        & xmax,
		  Vector<parametertype> & xtype,
		  Vector<double> & dir_io,
		  double tol=-1,
		  bool debug=false,
		  bool debug_f1dim=false,
		  bool debug_golden=false,
		  bool debug_bracket=false
		  ){
      // Input: n-dimensional point_io p and direction dir_io.
      // Purpose: Find minimum of supplied function/functor aint direction.
      // Sets p to the minimum point, and replaces dir_io with the vector
      // displacement from point_io to minimum point.
      // Uses bracket() and minimize() from Golden.
      // Output: Returns value of func at minimum point.
      int n;
      double fmin;

      status = funcfit::Minimization_Status();

      // cout.precision(15);

      p  = point_io;
      xi = dir_io;
      n = p.size();
    

      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      xi.normalize();
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




      F1dim<T> f1dim(p, xi, func);
      Golden< F1dim<T> > golden(f1dim);

      f1dim.debug = debug_f1dim;

      Vector<double> sminvec(n,0), smaxvec(n,0);
      double smin, smax;
      double dmax = sqrt( numeric_limits<double>::max() );
      double dmin = -dmax;

      //cout << "dmin dmax : " << dmin << " " << dmax << endl;

      /* f1dim takes scalar argument s, and calculates point q:
	   vec(q) = vec(p) + s * vec(xi)
	 At this point original func is evaluated and returned
	 as the value of f1dim(s).
      */
      for (int i=0; i<n; ++i){
	//cout << "xtype[" << i << "] " << xtype[i] << endl;

	if (xtype[i]==PARAM_FREE){
	  sminvec[i]=dmin;
	  smaxvec[i]=dmax;
	}
	else {
	  sminvec[i] = (xmin[i] - p[i])/xi[i];
	  smaxvec[i] = (xmax[i] - p[i])/xi[i];
	  if (sminvec[i]>smaxvec[i]){
	    double tmp = sminvec[i];
	    sminvec[i] = smaxvec[i];
	    smaxvec[i] = tmp;
	  }
	}
	if (i==0 || (i>0 && sminvec[i]>smin)) smin=sminvec[i];
	if (i==0 || (i>0 && smaxvec[i]<smax)) smax=smaxvec[i];
      }
      //cout << "smin smax : " << smin << " " << smax << endl;

      if (smin > smax || fp_are_equal(smin, smax)){
	// Cannot go anywhere!
	cout << "WARNING: Line minimization: Stuck at starting value due to parameter constraints." << endl;
	point_io = p;
	for (int i=0; i<n; ++i) dir_io[i] = 0.0;
	fmin = func(p);
	status.funcmin = fmin;
	return fmin;
      }


      // cout << "xi now: " << xi << endl;






      double ax=0.0, xx=1.0;


      //xx = 1.0;

      /*
      if (ax<smin) ax=smin;
      if (xx>smax) xx=smax;
      */
      double xx_old=xx, sxmin;
      // *********************************************************************
      while (true){
	try {
	  sxmin = golden.minimize(ax, xx, smin, smax, tol, debug_golden, debug_bracket);
	  fmin = golden.fmin;
	  status.funcmin = golden.fmin;
	}
	catch (funcfit::bad_point & e1){
	  func.reset();

	  // Went too far in the minimization. Restart with smaller limits.
	  if (debug)
	    cout << "Line minimization: Starting value xx " << xx << " too large. Retrying." << endl;
	  xx_old *= 0.5;
	  xx = xx_old;
	  continue;
	}
	break;
      }
      // *********************************************************************

      if (debug)
	cout << "Line minimization: xmin = " << sxmin << endl;


      for (int j=0; j!=n; ++j){
	xi[j] *= sxmin; // actual vectorial step moved
	p[j] += xi[j];
      }
      point_io = p;
      dir_io = xi;


      double hxi = xi.magn();
      if (fp_is_small(hxi)){
	status.small_step = true;
	fmin = golden.fmin;
	status.funcmin = golden.fmin;
	return fmin;
      }
      // *********************************************************************

      return fmin;


    }

  } ;



}

#endif

