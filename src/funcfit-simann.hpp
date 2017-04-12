


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
    int Nr;
    int NS;
    int NTe;
    double c;
    double rTe;
    double qTe;

    Cond_SimAnn(){
      Nr = -15;  // can edit this
      NS =  20;  // can edit this
      NTe = 100; // can edit this: rate of T and ChiSq descent need to be balanced!
      c = 2;
      rTe = 0.85;
      qTe = 0.67;
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
			    double simann_delta_rel,
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

      Vector<double>        x     = func.free_parameters(point_in);
      Vector<double>        xmin  = func.map_vector_as_free_parameters(par_min);
      Vector<double>        xmax  = func.map_vector_as_free_parameters(par_max);
      Vector<parametertype> xtype = func.map_vector_as_free_parameters(par_type);
      int nx = x.size();

      Vector<double> xch_ini(nx,0);
      Vector<double> xch(nx,0);

      for (int i=0; i<nx; ++i){
	xch_ini[i] = x[i] * simann_delta_rel;
	if (xch_ini[i] < 0) xch_ini[i] *= -1.0;
	if (xch_ini[i]<eps) xch_ini[i] = sqrt(eps);

	if (xch_ini[i] > xmax[i] - xmin[i]) xch_ini[i] = xmax[i] - xmin[i];
      }
      xch = xch_ini;

      status = Minimization_Status();
      bool debug         = cond_debug.debug_fit_level0;
      string prefix_report_debug = cond_debug.prefix_debug_fit_level0;
      bool report_iter    = cond_print.report_iter;
      bool report_warn    = cond_print.report_warn;
      bool report_error   = cond_print.report_error;
      Counters_niter counters_niter = Counters_niter();
      int niterrestart = cond_conv.niterrestart;
      int Nr    = cond_simann.Nr;
      int NS    = cond_simann.NS;
      int NTe    = cond_simann.NTe;
      double cc0 = cond_simann.c;
      double rTe = cond_simann.rTe;
      double qTe = cond_simann.qTe;
      double cool_rate = (1.0 - qTe*rTe)/NTe;

      Vector<double> cc(nx, 0);

      int ci, cj, cm, ch;
      Vector<int> cn(nx, 0);

      Vector<double> x_try(nx), x_old(nx), h(nx), x_opt(nx), v(nx,0);
      double fp, fp_try, dfp, fp_opt, fp_old;
      double Te0, Te1_star, Te2_star, Te;

      double hmagn, td;

      string dumpfile;
      ofstream fout;
      string steptype;
      bool trial_failed, stoch_red;
      int Nfail=0, Nred=0, Nstoch=0;

      niter = 0;
      ci=0;
      cj=0;
      cm=0;
      ch=0;
      for (int i=0; i<nx; ++i) cn[i]=0;


      double dmax = sqrt( numeric_limits<double>::max() );
      double dmin = -dmax;
      for (i=0; i<nx; ++i){
	if (xtype[i]==PARAM_FREE){ // no lower or upper limit
	  xmin[i] = dmin;
	  xmax[i] = dmax;
	}
      }




      if (debug) 
	cout << prefix_report_debug
	     << methodstring << ": "
	     << "Initial point ... " << endl;
      for (i=0; i<nx; ++i){
	x_opt[i] = x[i];
      }
      fp_opt = fp = func(x);

      // Initial temperature:
      Te0   = abs(fp) / nx; //* 1e-3;
      //Te0 = sqrt(abs(fp));
      //Te0   = abs(fp) * nx; //* 1e-3;

      Te1_star = Te0;
      Te2_star = qTe * Te1_star;
      Te = Te1_star * (1.0 - cool_rate * cm);



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
	  for (i=0; i<nx; ++i)
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
	  func.report_on_parameters_and_data();
	}



	// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// Check for convergence
	// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	if ( check_conv(x, true,fp_old, true,fp, // use old and current function values
			false,Vector<double>(nx,0), // use gradient
			true,h,                           // use step taken
			niter,
			counters_niter,
			cond_conv, cond_debug, cond_print,
			methodstring, status) ){
	  return func.all_parameters(x);
	}


	// Should we restart?
        if (niterrestart > 0 && niter % niterrestart == 0){

	  for (int i=0; i<nx; ++i){
	    xch_ini[i] = x[i] * simann_delta_rel;
	    if (utils::abs(xch_ini[i])<eps) xch_ini[i] = sqrt(eps);
	    if (xch_ini[i] < 0) xch_ini[i] *= -1.0;

	    if (xch_ini[i] > xmax[i] - xmin[i]) xch_ini[i] = xmax[i] - xmin[i];
	  }
	  xch = xch_ini;

	}



	x_old  = x;
	fp_old = fp;      
	



	// ###################################################################
	// Get trial point:
	// ###################################################################
	for (ch=0; ch<nx; ++ch){

	  stoch_red = false;
	  
	  // **************************************************************
	  while (true){
	    try {

	      // Start with previous point ...
	      for (i=0; i<nx; ++i) x_try[i] = x[i];

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


	  // ******************************************************
	  // Check if new point is to be accepted or not
	  // ******************************************************
	  if (dfp < 0.0){
	    // --------------------------------------------
	    // Reduction successful
	    // --------------------------------------------
	    fp = fp_try;
	    x = x_try;
	    ci++;
	    cn[ch]++;

	    //Te = 0.90 * Te;
	    //xch[ch] *= 1.2;

	    // Check if new function value is better than current best value:
	    if (fp < fp_opt){
	      fp_opt = fp;
	      for (i=0; i<nx; ++i) x_opt[i] = x[i];
	    }
	  }
	  else {
	    // --------------------------------------------
	    // Stochastic trial acceptance of increased merit function value
	    // --------------------------------------------
	    if (Te >= eps && mtwister.unif() < exp( - dfp/Te  ) ){
	      // Stochastic try succeeded
	      stoch_red = true;
	      fp = fp_try;
	      x = x_try;
	      ci++;
	      cn[ch]++;

	      //xch[ch] *= 0.8;
	    }
	    else {
	      stoch_red = false;
	      //Te *= 1.05;
	      //xch[ch] *= 0.5;
	    }
	  }
	

	  
	} // End of loop over all directions.
	ch=0;
	cj++;
	// ###################################################################
	// Trial point obtained
	// ###################################################################

     


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





	//if (niter>1) Te = Te0 / log(niter);






	if (NS > 0 && cj >= NS){
	  // ###################################################################
	  // Update step vector:
	  // ###################################################################
	  Vector<double> xch_new(xch);
	  for (int u=0; u<nx; ++u){
	    if      (cn[u] - 0.6*NS > 0.0)
	      xch_new[u] = xch[u] * (1.0 + cc[u]*(cn[u]*1.0/NS - 0.6)/0.4);
	    else if (cn[u] - 0.4*NS < 0.0)
	      xch_new[u] = xch[u] / (1.0 + cc[u]*(0.4 - cn[u]*1.0/NS)/0.4);
	    else
	      xch_new[u] = xch[u];
	  }
	  xch = xch_new;
	  cj=0;
	  for (int i=0; i<nx; ++i) cn[i] = 0;


	  /*
	  if (debug)
	    cout << "Setting state to optimum (optimum point and function value) ..." << endl;
	  for (i=0; i<nx; ++i) x[i] = x_opt[i];
	  for (i=0; i<nx; ++i) xch[i] = xch_ini[i];
	  fp = fp_opt;
	  */
	}




	// ###################################################################
	// Update temperature:
	// ###################################################################
	cm++;
	if (cm >= NTe){
	  Te1_star = Te2_star;
	  if (Te1_star < eps){
	    std::cout << methodstring << ": Temperature has reached 0.0, aborting." << std::endl;
	    return func.all_parameters(x_opt);
	  }
	  Te2_star = qTe * Te1_star;
	  cm=0;
	  if (Nr < 0){
	    x  = x_opt;
	    fp = fp_opt;
	  }
	}

	Te = Te1_star * (1.0 - cool_rate * cm);


	if (Nr > 0 && niter % Nr == 0){
	  x  = x_opt;
	  fp = fp_opt;
	}




	niter++;



	

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



