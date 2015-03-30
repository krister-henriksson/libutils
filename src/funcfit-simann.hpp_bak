


#ifndef FUNCFIT_SIMANN_HPP
#define FUNCFIT_SIMANN_HPP


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
  [1]: Corana et al, ACM Transactions on Mathematical Software 13 (1987) 262-280

 */

///////////////////////////////////////////////////////////
// Usage:
///////////////////////////////////////////////////////////
//
//   Funcd fd;
//
//   SimAnn<Funcd> fgn(fd);
//
//   xmin = fgn.minimize(...);
//   fmin = fgn.status.fmin;
//
///////////////////////////////////////////////////////////
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


  class Cond_SimAnn {
  public:
    int NTe;
    double rTe;
    int anneal_method;
    double quench_rate;
    int niter_reset;

    Cond_SimAnn(){
      NTe = 4;
      rTe = 0.85;
      anneal_method = 2;
      quench_rate = 0.01;
      niter_reset = 20;
     }

  } ;




  
  template <typename T>
  class SimAnn
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
    SimAnn(T & funcd)
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
			    Vector<double> par_charstepsize,
			    int seed,
			    Cond_Conv  cond_conv  = Cond_Conv(),
			    Cond_Debug cond_debug = Cond_Debug(),
			    Cond_Print cond_print = Cond_Print(),
			    Cond_SimAnn cond_simann = Cond_SimAnn()
			    ){
      // Input: starting point point and constructor-set functor (containing
      // gradient method)
      // Purpose: minimize function, input point will not be changed
      double eps = numeric_limits<double>::epsilon();
      double small = sqrt(eps);

      string methodstring("simulated-annealing");
      int niter=0, i,j,k;

      int seed2 = seed < 0 ? time(0) : seed;
      rand_mtwister mtwister( seed2 );

      Vector<double> point   = func.free_parameters(point_in);
      Vector<double>        xmin  = func.map_vector_as_free_parameters(par_min);
      Vector<double>        xmax  = func.map_vector_as_free_parameters(par_max);
      Vector<parametertype> xtype = func.map_vector_as_free_parameters(par_type);
      Vector<double> xch_ini = func.map_vector_as_free_parameters(par_charstepsize);
      Vector<double> xch = xch_ini;

      int D  = point.size();

      int NTe     = cond_simann.NTe;
      double rTe   = cond_simann.rTe;
      int anneal_method  = cond_simann.anneal_method;
      double quench_rate = cond_simann.quench_rate;
      int niter_reset   = cond_simann.niter_reset;

      double Te, Te0, Te_star;
      //int NT = max((int)100, (int)5*D);
      int NT = max((int)100, (int)5*D);

      status = Minimization_Status();

      bool debug         = cond_debug.debug_fit_level0;
      string prefix_report_debug = cond_debug.prefix_debug_fit_level0;

      bool report_iter    = cond_print.report_iter;
      bool report_warn    = cond_print.report_warn;
      bool report_error   = cond_print.report_error;

      Counters_niter counters_niter = Counters_niter();
      
      Vector<double> x(D), x_try(D), x_old(D), h(D), x_opt(D), v(D,0);
      double fp, fp_try, dfp, fp_opt, fp_old, fp_star, hmagn, td;
      int ci, ck, cm, ch;

      string dumpfile;
      ofstream fout;
      string steptype;
      bool trial_failed, stoch_red;
      int Nfail=0, Nred=0, Nstoch=0;

      niter = 0;
      ch=0;
      ci=0;
      ck=0;
      cm=0;



      double dmax = sqrt( numeric_limits<double>::max() );
      double dmin = -dmax;
      for (i=0; i<D; ++i){
	if (xtype[i]==PARAM_FREE){ // no lower or upper limit
	  xmin[i] = dmin;
	  xmax[i] = dmax;
	}
      }




      // ###################################################################
      // Create population
      // ###################################################################
      if (debug) 
	cout << prefix_report_debug
	     << methodstring << ": "
	     << "Initial point ... " << endl;
      for (i=0; i<D; ++i){
	x[i] = point[i];
	x_opt[i] = x[i];
      }
      fp_opt = fp = func(x);
      // Initial temperature:
      Te0 = abs(fp) / D; //* 1e-3;
      Te = Te0;
      Te_star = Te0;



      if (debug)
	cout << prefix_report_debug
	     << methodstring << ": "
	     << "Starting iterations ... " << endl;



      


      // ###################################################################
      // ###################################################################
      // ###################################################################
      // ###################################################################
      while (true){ // Start of loop over iterations


	func.set_point(x);


	if (debug){
	  for (i=0; i<D; ++i)
	    cout << "x[" << i << "]= " << x[i] << "  xch[" << i << "]= " << xch[i] << endl;

	}




	// *******************************************************************
	// Report
	// *******************************************************************
	if (report_iter){
	  double vb = func.value_barrier();
	  
	  if (niter==0){
	    printf("%s%s: Iter %4d  Func opt %15.8e now %15.8e               "
		   "Func_barrier %15.8e "
		   "Te %15.8e                             "
		   "Nred %4d  Nfail %4d\n",
		   cond_print.prefix_report_iter.c_str(),
		   methodstring.c_str(), niter, fp_opt, fp,
		   vb,
		   Te,
		   Nred, Nfail);
	  }
	  else {
	    printf("%s%s: Iter %4d  Func opt %15.8e now %15.8e change %15.8e "
		   "Func_barrier %15.8e "
		   "Te %15.8e stepsize %15.8e step: %12s  "
		   "Nred %4d  Nfail %4d\n",
		   cond_print.prefix_report_iter.c_str(),
		   methodstring.c_str(), niter, fp_opt, fp, fp-fp_old,
		   vb,
		   Te, hmagn, steptype.c_str(),
		   Nred, Nfail);
	  }
	}


	// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// Check for convergence
	// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	if ( check_conv(x, true,fp_old, true,fp, // use old and current function values
			false,Vector<double>(h.size(),0), // use gradient
			true,h,                           // use step taken
			niter,
			counters_niter,
			cond_conv, cond_debug, cond_print,
			methodstring, status) ){
	  return func.all_parameters(x);
	}








	fp_old = fp;
	x_old  = x;


	// ###################################################################
	// Get trial point:
	// ###################################################################
	stoch_red = false;
	for (ch=0; ch<D; ++ch){

	  // **************************************************************
	  while (true){
	    try {

	      // Start with previous point ...
	      for (i=0; i<D; ++i) x_try[i] = x[i];

	      trial_failed = true;
	      while (trial_failed){
		trial_failed = false;

		// ... and add a trial step in direction 'ch':
		td = (2.0*mtwister.unif()-1.0) * xch[ch];
		if (x_try[ch]+td < xmin[ch]) trial_failed=true;
		if (x_try[ch]+td > xmax[ch]) trial_failed=true;

	      }
	      x_try[ch] += td;
	      fp_try = func(x_try);
	    }
	    catch (funcfit::bad_point & err_bad_point){
	      func.reset();
	      cout << "Warning: Bad point " << x_try << ". Recreating point ..." << endl;
	      cout << "Recommendation: Restate parameter limits and start over." << endl;
	      continue;
	    }
	    break;
	  }
	  // **************************************************************
	  
	  dfp = fp_try - fp;

	  // Check trial function value:
	  if (dfp < 0.0 || dfp < eps){
	    // Reduction successful
	    fp = fp_try;
	    x = x_try;
	    ci++;
	    //Te = 0.90 * Te;
	    xch[ch] *= 1.2;
	  }
	  else {
	    
	    if (Te >= eps && mtwister.unif() < exp( - dfp/Te ) ){
	      // Stochastic try succeeded
	      stoch_red = true;
	      fp = fp_try;
	      x = x_try;
	      ci++;
	      xch[ch] *= 0.8;
	    }
	    else {
	      //Te *= 1.05;
	      xch[ch] *= 0.5;
	    }

	  }


	  // Check if new function value is better than current best value:
	  if (fp <= fp_opt){
	    fp_opt = fp;
	    for (i=0; i<D; ++i) x_opt[i] = x[i];
	  }


	} // End of loop over all directions.


	h = x - x_old;
	hmagn = h.magn();
	if (fp < fp_old || stoch_red){
	  steptype = "reduction";
	  Nred++;
	}
	else {
	  steptype = "NO_reduction";
	  Nfail++;
	}





	niter++;

	//if (niter>1) Te = Te0 / log(niter);




	// ###################################################################
	// Update temperature:
	// ###################################################################
	if (anneal_method==1){
	  if (cm < NTe) continue;

	  Te_star = Te;
	  fp_star = fp;
	  
	  Te = rTe * Te;

	  ck++;
	  cm=0;
	}
	else if (anneal_method==2){
	  Te = Te0 * (1.0 - quench_rate * niter);
	  if (Te<=eps) Te = 0.0;
	}





	if (niter % niter_reset == 0){

	  // ###################################################################
	  // Set x to optimal point:
	  // ###################################################################
	  ci++;
	  if (debug)
	    cout << "Setting state to optimum (optimum point and function value) ..." << endl;
	  for (i=0; i<D; ++i) x[i] = x_opt[i];
	  for (i=0; i<D; ++i) xch[i] = xch_ini[i];
	  fp = fp_opt;

	}


	

	// Go to next iteration.
      }
      // ###################################################################
      // ###################################################################
      // ###################################################################
      // ###################################################################



      fmin = fp_opt;
      return func.all_parameters(x_opt);
      


    } // End of method

  } ; // End of class


} // End of namespace

#endif



