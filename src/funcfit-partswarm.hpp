


#ifndef FUNCFIT_PARTSWARM_HPP
#define FUNCFIT_PARTSWARM_HPP


#include <iostream>
#include <string>
#include <fstream>
#include <limits>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>

#include <boost/format.hpp>

#include "utils.hpp"
#include "utils-math.hpp"
#include "utils-vector.hpp"
#include "utils-string.hpp"
#include "utils-errors.hpp"

#include "param.hpp"

#include "mtwister.hpp"

#include "funcfit-basics.hpp"
#include "funcfit-errors.hpp"





using namespace utils;
using std::cout;
using std::endl;
using std::ofstream;
using std::numeric_limits;
using boost::format;


namespace funcfit {


  class Cond_PartSwarm {
  public:
    double c1;
    double c2;
    double wrate;
    double wmax;
    double wmin;
    int NP_fac;
    

    Cond_PartSwarm(){
      c1 = 2;
      c2 = 2;
      wrate = 0.01;//0.1; // large value allows for large contribution of old speed
      wmax  = 0.9;
      wmin  = 0.1;
      NP_fac = 10;
    }

  } ;




  
  template <typename T>
  class PartSwarm
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
    PartSwarm(T & funcd)
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
			    Cond_PartSwarm cond_partswarm = Cond_PartSwarm()
			    ){
      // Input: starting point point and constructor-set functor (containing
      // gradient method)
      // Purpose: minimize function, input point will not be changed
      double eps = numeric_limits<double>::epsilon();
      double eps2r = sqrt( numeric_limits<double>::epsilon() );
      string methodstring("particle-swarm");
      int niter=0, i,j,k, ixglobmin, ixglobmax;
      double vmax, w;

      int seed2 = seed < 0 ? time(0) : seed;
      rand_mtwister mtwister( seed2 );

      Vector<double> point = func.free_parameters(point_in);
      Vector<double>        xmin  = func.map_vector_as_free_parameters(par_min);
      Vector<double>        xmax  = func.map_vector_as_free_parameters(par_max);
      Vector<parametertype> xtype = func.map_vector_as_free_parameters(par_type);


      double c1    = cond_partswarm.c1;
      double c2    = cond_partswarm.c2;
      double wrate = cond_partswarm.wrate;
      double wmax  = cond_partswarm.wmax;
      double wmin  = cond_partswarm.wmin;
      int NP_fac  = cond_partswarm.NP_fac;

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
      Vector<double> fx(NP);

      double td1,td2,td3, rtol;
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
      
      // --------------------------------------------------------------
      // --------------------------------------------------------------
      // Create initial population vectors and values
      // --------------------------------------------------------------
      // --------------------------------------------------------------
      if (debug) 
	cout << prefix_report_debug
	     << methodstring << ": "
	     << "Creating initial population vectors and values ... " << endl;
      for (i=0; i<NP; ++i){
	// cout << methodstring << ": " << "Agent " << i+1 << " of " << NP << endl;

	// **************************************************************
	while (true){
	  try {
	    // configuration of individual
	    for (j=0; j<D; ++j){
	      x[i][j] = xmin[j] + mtwister.unif() * (xmax[j] - xmin[j]);
	      if (x[i][j] - eps2r < xmin[j]) x[i][j] = xmin[j] + eps2r;
	      if (x[i][j] + eps2r > xmax[j]) x[i][j] = xmax[j] - eps2r;

	      v[i][j] = 0.5*(xmax[j]-xmin[j]) * (2 * mtwister.unif() - 1);
	    }
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



	// ***************************************************************
	// ***************************************************************
	// Report
	// ***************************************************************
	// ***************************************************************
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
	  //cout << "Chi^2 components: " << func.f() << endl;
	  //cout << "DataY: " << func.DataY() << endl;
	  //cout << "ModelDataY: " << func.ModelDataY() << endl;
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
	// Rescaling of w
	// ###################################################################
	w = wmax - (wmax - wmin)*wrate * niter;
	if (w<wmin) w=wmin;




	// ###################################################################
	// Update
	// ###################################################################
	if (debug)
	  cout << prefix_report_debug
	       << methodstring << ": "
	       << "Getting merit function values for initial population ... " << endl;

	for (i=0; i<NP; ++i){
	  // cout << methodstring << ": " << "Agent " << i+1 << " of " << NP << endl;

	  // **************************************************************
	  while (true){
	    try {
	      for (j=0; j<D; ++j){
		v[i][j] = w * v[i][j]
		  + c1 * mtwister.unif() * (xindmin[i][j] - x[i][j])
		  + c2 * mtwister.unif() * (xglobmin[j]   - x[i][j]);
	    
		vmax = 0.5*(xmax[j]-xmin[j]);
		if (v[i][j] + eps2r >  vmax) v[i][j] =  vmax - eps2r;
		if (v[i][j] - eps2r < -vmax) v[i][j] = -vmax + eps2r;
	      }

	      for (j=0; j<D; ++j){
		x[i][j] = x[i][j] + v[i][j];

		if (x[i][j] - eps2r < xmin[j]) x[i][j] = xmin[j] + eps2r;
		if (x[i][j] + eps2r > xmax[j]) x[i][j] = xmax[j] - eps2r;

	      }
	  
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

	  if (fxindmin[i]>fx[i]){
	    for (j=0; j<D; ++j) xindmin[i][j] = x[i][j];
	    fxindmin[i] = fx[i];
	  }
	  if (fxindmax[i]<fx[i]){
	    for (j=0; j<D; ++j) xindmax[i][j] = x[i][j];
	    fxindmax[i] = fx[i];
	  }

	  if (i==0 || (i>0 && fxglobmin>fx[i])){ ixglobmin = i; fxglobmin = fx[i]; }
	  if (i==0 || (i>0 && fxglobmax<fx[i])){ ixglobmax = i; fxglobmax = fx[i]; }
	}
	xglobmin = x[ixglobmin];
	xglobmax = x[ixglobmax];
	fmin = fxglobmin;


	dumpfile = "points-PS-iter" + tostring(niter) + ".out";
	fout.open(dumpfile.c_str());
	for (i=0; i<NP; ++i){
	  fout << format("%10d  ") % i << " : ";
	  for (j=0; j<D; ++j) fout << format(" %15.8e") % x[i][j];
	  fout << format(" %15.8e") % fx[i] << endl;
	}
	fout.close();







	// Go to next iteration.
	niter++;
      } // End of loop over iterations
      // ###################################################################
      // ###################################################################
      // ###################################################################
      // ###################################################################



      
      fmin = fxglobmin;
      return func.all_parameters(x[ixglobmin]);


    } // End of method.



  };



}





#endif



