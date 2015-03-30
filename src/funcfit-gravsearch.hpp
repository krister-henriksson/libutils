


#ifndef FUNCFIT_GRAVSEARCH_HPP
#define FUNCFIT_GRAVSEARCH_HPP


#include <iostream>
#include <string>
#include <fstream>
#include <limits>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>

#include <boost/format.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

#include "utils.hpp"
#include "utils-math.hpp"
#include "utils-vector.hpp"
#include "utils-string.hpp"
#include "utils-errors.hpp"

#include "param.hpp"

#include "mtwister.hpp"

#include "funcfit-basics.hpp"
#include "funcfit-errors.hpp"



/*
  [1]: Rashedi et al, Information Sciences 179 (2009) 2232-2248
  (modified)


 */

///////////////////////////////////////////////////////////
// Usage:
///////////////////////////////////////////////////////////
//
//   Funcd fd;
//
//   FuncFitDE<Funcd> fgn(fd);
//
//   xmin = fgn.minimize(...);
//   fmin = fgn.status.funcmin;
//
///////////////////////////////////////////////////////////
//
//  class Funcd {
//    // * Return function value:
//    double operator()(void);
//    double operator()(utils::Vector<double> & x);
//    // * Return free parameters:
//    utils::Vector<double> free_parameters(void);
//    utils::Vector<double> free_parameters(const utils::Vector<double> & x);
//    // * Return all parameters:
//    utils::Vector<double> all_parameters(void);
//    utils::Vector<double> all_parameters(const utils::Vector<double> & x);
//
//  }
//
///////////////////////////////////////////////////////////


using namespace utils;
using namespace funcfit;
using std::cout;
using std::endl;
using std::ofstream;
using std::numeric_limits;
using boost::format;
using boost::math::isfinite;


namespace funcfit {


  class Cond_GravSearch {
  public:
    double Gstart;
    double Gend;
    double Grate;
    int NP_fac;

    Cond_GravSearch(){
      Gstart = 100;
      Gend = 1;
      Grate = 0.01;
      NP_fac = 10;
     }

  } ;




  
  template <typename T>
  class GravSearch
  {


    //#####################################################################
    // Member variables
    //#####################################################################
  public:
    double fmin;
    T & func;
    Minimization_Status status;



    //#####################################################################
    // Constructor
    //#####################################################################
    // NB:  Use constructor to set functor. Point and direction are set by
    // NB:  minimize() method.
    GravSearch(T & funcd)
      :
      func(funcd)
    {}



    //#####################################################################
    // Main method
    //#####################################################################
    Vector<double> minimize(Vector<double> & point_in,
			    Vector<double>        & par_min,
			    Vector<double>        & par_max,
			    Vector<parametertype> & par_type,
			    int seed,
			    Cond_Conv  cond_conv  = Cond_Conv(),
			    Cond_Debug cond_debug = Cond_Debug(),
			    Cond_Print cond_print = Cond_Print(),
			    Cond_GravSearch cond_gravsearch = Cond_GravSearch()
			    ){
      // Input: starting point point and constructor-set functor (containing
      // gradient method)
      // Purpose: minimize function, input point will not be changed
      double eps = numeric_limits<double>::epsilon();
      double eps2r = sqrt( numeric_limits<double>::epsilon() );
      double small = sqrt(eps);

      string methodstring("grav-search");
      int niter=0, i,j,k, ixglobmin, ixglobmax;
      double vmax;

      int seed2 = seed < 0 ? time(0) : seed;
      rand_mtwister mtwister( seed2 );

      Vector<double> point = func.free_parameters(point_in);
      Vector<double>        xmin  = func.map_vector_as_free_parameters(par_min);
      Vector<double>        xmax  = func.map_vector_as_free_parameters(par_max);
      Vector<parametertype> xtype = func.map_vector_as_free_parameters(par_type);


      double Gstart = cond_gravsearch.Gstart;
      double Gend   = cond_gravsearch.Gend;
      double Grate  = cond_gravsearch.Grate;
      int NP_fac  = cond_gravsearch.NP_fac;

      int D  = point.size();
      int NP = NP_fac * D;

      status = Minimization_Status();

      bool debug         = cond_debug.debug_fit_level0;
      string prefix_report_debug = cond_debug.prefix_debug_fit_level0;

      bool report_iter    = cond_print.report_iter;
      bool report_warn    = cond_print.report_warn;
      bool report_error   = cond_print.report_error;

      Counters_niter counters_niter = Counters_niter();

      Vector< Vector<double> > xindmin(NP, Vector<double>(D,0));
      Vector< Vector<double> > xindmax(NP, Vector<double>(D,0));
      Vector<double> fxindmin(NP, 0);
      Vector<double> fxindmax(NP, 0);
      
      Vector<double> xglobmin(D, 0);
      Vector<double> xglobmax(D, 0);
      double fxglobmin, fxglobmax;
      double fxglobmin_old, fxglobmax_old;
      
      Vector< Vector<double> > x(NP, Vector<double>(D,0));
      Vector< Vector<double> > v(NP, Vector<double>(D,0));
      Vector< Vector<double> > a(NP, Vector<double>(D,0));
      Vector< Vector<double> > f(NP, Vector<double>(D,0));
      Vector<double> mass(NP);
      Vector<double> fx(NP);
      double masstot, dist, td1,td2,td3, rtol;
      double G = Gstart;

      string dumpfile;
      ofstream fout;

      double dmax = sqrt( numeric_limits<double>::max() );
      double dmin = -dmax;
      for (i=0; i<D; ++i){
	if (xtype[i]==PARAM_FREE){ // no lower or upper limit
	  xmin[i] = dmin;
	  xmax[i] = dmax;
	}
      }


      niter = 0;


      // ###################################################################
      // Create agents
      // ###################################################################
      if (debug) 
	cout << prefix_report_debug
	     << methodstring << ": "
	     << "Creating agents ... " << endl;
      for (i=0; i<NP; ++i){
	// cout << methodstring << ": " << "Agent " << i+1 << " of " << NP << endl;

	// **************************************************************
	while (true){
	  try {
	    for (j=0; j<D; ++j){
	      x[i][j] = xmin[j] + mtwister.unif() * (xmax[j] - xmin[j]);
	      if (x[i][j] - eps2r < xmin[j]) x[i][j] = xmin[j] + eps2r;
	      if (x[i][j] + eps2r > xmax[j]) x[i][j] = xmax[j] - eps2r;

	    }
	    fx[i] = func(x[i]);
	    if ( ! isfinite( fx[i] ) ){
	      i--;
	      continue;
	    }
	  }
	  catch (funcfit::bad_point & err_bad_point){
	    func.reset();
	    cout << "Warning: Bad point " << x[i] << ". Recreating point ..." << endl;
	    cout << "Recommendation: Restate parameter limits and start over." << endl;
	    continue;
	  }
	  break;
	}
	// **************************************************************

	for (j=0; j<D; ++j){
	  xindmin[i][j] = x[i][j];
	  xindmax[i][j] = x[i][j];
	}
	fxindmin[i] = fx[i];
	fxindmax[i] = fx[i];

	if (i==0 || (i>0 && fxglobmin>fx[i])){ ixglobmin = i; fxglobmin = fx[i]; }
	if (i==0 || (i>0 && fxglobmax<fx[i])){ ixglobmax = i; fxglobmax = fx[i]; }
      }
      xglobmin = x[ixglobmin];
      xglobmax = x[ixglobmax];
      fmin = fxglobmin;

      cout << "fxglobmax " << fxglobmax << endl;
      cout << "fxglobmin " << fxglobmin << endl;





      

      if (debug)
	cout << prefix_report_debug
	     << methodstring << ": "
	     << "Starting iterations ... " << endl;



      


      // ###################################################################
      // ###################################################################
      // ###################################################################
      // ###################################################################
      while (true){ // Start of loop over iterations

	if (debug) cout << "Iteration " << niter << ":" << endl;
	if (debug)
	  for (i=0; i<x.size(); ++i)
	    cout << "x[" << i << "]: " << x[i] << endl;




	// *******************************************************************
	// Report
	// *******************************************************************
	if (report_iter){
	  // #################################################################
	  // Last point evaluated by the function is one of the vertices/individuals
	  // (ex. it could be the worst point or the best point). Now need
	  // to report about the best point:
	  Vector<double> pmin(D, 0.0);
	  for (i=0; i<D; ++i) pmin[i] = x[ixglobmin][i];
	  func(pmin); // triggers recalculation of properties
	  // #################################################################

	  printf("%s%s: Iter %4d  Func min %15.8e max %15.8e\n", 
		 cond_print.prefix_report_iter.c_str(),
		 methodstring.c_str(), niter, fxglobmin, fxglobmax);
	  func.report_on_parameters_and_data();
	}

	// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// Check for convergence
	// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	if ( check_conv(x[ixglobmin], true,fxglobmin_old, true,fxglobmin, // use old and current function values
			false,Vector<double>(x[ixglobmin].size(),0),// use gradient
			false,Vector<double>(x[ixglobmin].size(),0),// use step taken
			niter,
			counters_niter,
			cond_conv, cond_debug, cond_print,
			methodstring, status) ){
	  return func.all_parameters(x[ixglobmin]);
	}


	fxglobmin_old = fxglobmin;
	fxglobmax_old = fxglobmax;

	

	// ###################################################################
	// Masses
	// ###################################################################
	if (debug)
	  cout << prefix_report_debug
	       << methodstring << ": "
	       << "Getting masses ... " << endl;
	masstot = 0.0;
	for (i=0; i<NP; ++i){
	  /*
	  //mass[i] = (fx[i] - fxglobmax) / (fxglobmin - fxglobmax);
	  if (i==ixglobmax)
	    mass[i] = small / fxglobmin;
	  else
	    mass[i] = (fx[i] - fxglobmax + small) / fxglobmin;
	  */
	  td1 = fx[i]; if (td1<0) td1 *= -1.0;
	  mass[i] = 1.0/(1.0 + td1);
	  masstot += mass[i];
	}
	masstot = 1.0/masstot;
	//for (i=0; i<NP; ++i) mass[i] *= masstot;

	if (debug)
	  for (i=0; i<NP; ++i){
	    cout << "mass of agents " << i << " is " << mass[i] << endl;
	  }
	


	// ###################################################################
	// Update of gravtitational constant G
	// ###################################################################
	G = Gstart - niter * Grate;
	if (G<Gend) G = Gend;



	// ###################################################################
	// Forces
	// ###################################################################
	if (debug)
	  cout << prefix_report_debug
	       << methodstring << ": "
	       << "Getting forces ... " << endl;
	for (i=0; i<NP; ++i){ // loop over agents
	  for (k=0; k<D; ++k) f[i][k] = 0;
	  dist = 0.0;
	  for (j=0; j<NP; ++j){ // loop over neighboring agents
	    for (k=0; k<D; ++k){
	      // Gravitational force:
	      f[i][k] += mtwister.unif() * G * mass[i] * mass[j] * (x[i][k] - x[j][k]);
	      td1 = x[i][k] - x[j][k];
	      dist += td1*td1;
	    }
	    dist = 1.0/(sqrt(dist) + small);
	    for (k=0; k<D; ++k)
	      f[i][k] *= dist;
	  }
	}


	// ###################################################################
	// Update of velocities and positions
	// Function values
	// ###################################################################
	if (debug)
	  cout << prefix_report_debug
	       << methodstring << ": "
	       << "Updating positions and getting merit function values for agents ... " << endl;
	for (i=0; i<NP; ++i){
	  // cout << methodstring << ": " << "Agent " << i+1 << " of " << NP << endl;

	  // **************************************************************
	  while (true){
	    try {

	      td1 = 0.0;
	      for (k=0; k<D; ++k){
		v[i][k] = mtwister.unif() * v[i][k] + f[i][k]/mass[i];
		vmax = 0.5*(xmax[k] - xmin[k]);
		if (v[i][k] - eps2r < -vmax) v[i][k] = -vmax + eps2r;
		if (v[i][k] + eps2r >  vmax) v[i][k] =  vmax - eps2r;

		td2 = x[i][k];
		x[i][k] = x[i][k] + v[i][k];

		if (x[i][k] - eps2r < xmin[k]) x[i][k] = xmin[k] + eps2r;
		if (x[i][k] + eps2r > xmax[k]) x[i][k] = xmax[k] - eps2r;
		td1 += (x[i][k] - td2)*(x[i][k] - td2);
	      }
	      td1 = sqrt(td1);
	      if (debug) cout << "step taken for agent " << i << " is " << td1 << endl;
	
	      fx[i] = func(x[i]);
	    }
	    catch (funcfit::bad_point & err_bad_point){
	      func.reset();
	      cout << "Warning: Bad point " << x[i] << ". Recreating point ..." << endl;
	      cout << "Recommendation: Restate parameter limits and start over." << endl;
	      continue;
	    }
	    break;
	  }
	  // **************************************************************

	  for (j=0; j<D; ++j){
	    xindmin[i][j] = x[i][j];
	    xindmax[i][j] = x[i][j];
	  }
	  fxindmin[i] = fx[i];
	  fxindmax[i] = fx[i];

	  if (i==0 || (i>0 && fxglobmin>fx[i])){ ixglobmin = i; fxglobmin = fx[i]; }
	  if (i==0 || (i>0 && fxglobmax<fx[i])){ ixglobmax = i; fxglobmax = fx[i]; }
	}
	xglobmin = x[ixglobmin];
	xglobmax = x[ixglobmax];
	fmin = fxglobmin;







	niter++;
	// Go to next iteration.
      }




      fmin = fxglobmin;
      return func.all_parameters(x[ixglobmin]);


    } // End of method

  } ; // End of class


} // End of namespace

#endif



