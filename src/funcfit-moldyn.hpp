


#ifndef FUNCFIT_MOLDYN_HPP
#define FUNCFIT_MOLDYN_HPP


#include <iostream>

#include <cstdio>
#include <cstdlib>
#include <string>
#include <cmath>


#include "utils.hpp"
#include "utils-math.hpp"
#include "utils-vector.hpp"

#include "param.hpp"

#include "funcfit-basics.hpp"
#include "funcfit-errors.hpp"

using std::cout;
using std::endl;
using utils::Vector;


namespace funcfit {

  
  template <typename T>
  class MolDynFit
  {


    //#####################################################################
    // Member variables
    //#####################################################################
  public:
    T & func;
    Minimization_Status status;



    //#####################################################################
    // Constructor
    //#####################################################################
    // NB:  Use constructor to set functor. Point and direction are set by
    // NB:  minimize() method.
    MolDynFit(T & funcd)
      :
      func(funcd)
    {}



    //#####################################################################
    // Main method
    //#####################################################################
    Vector<double> minimize(Vector<double>        & point_in,
			    Vector<double>        & par_min,
			    Vector<double>        & par_max,
			    Vector<parametertype> & par_type,
			    Cond_Conv  cond_conv  = Cond_Conv(),
			    Cond_Debug cond_debug = Cond_Debug(),
			    Cond_Print cond_print = Cond_Print()
			    ){
      // Input: starting point point and constructor-set functor (containing
      // gradient method)
      // Purpose: minimize function, input point will not be changed

      double fx, fx_old, gmagn;
      string methodstring;
      int i, j, nx, niter=0;

      status = Minimization_Status();

      bool debug         = cond_debug.debug_fit_level0;
      string prefix_report_debug = cond_debug.prefix_debug_fit_level0;

      bool report_iter    = cond_print.report_iter;
      bool report_warn    = cond_print.report_warn;
      bool report_error   = cond_print.report_error;

      Counters_niter counters_niter = Counters_niter();


      Vector<double>        x     = func.free_parameters(point_in);
      Vector<double>        xmin  = func.map_vector_as_free_parameters(par_min);
      Vector<double>        xmax  = func.map_vector_as_free_parameters(par_max);
      Vector<parametertype> xtype = func.map_vector_as_free_parameters(par_type);

      nx = x.size();

      Vector<double> x_trial(nx, 0.0);
      Vector<double> h(nx, 0.0);
      Vector<double> x_old(nx, 0.0);
      Vector<double> xvel(nx, 0.0), xvel_trial(nx, 0.0);
      Vector<double> xmass(nx, 1.0);
      Vector<double> xacc(nx,  0.0);
      Vector<double> g(nx,  0.0);
      double Ep=0, Ek=0, Et=0;
      double t=0, dt=1.0;


      Vector<double> tv;



      methodstring = "moldyn";




      if (debug)
	cout << prefix_report_debug
		  << methodstring << ": "
		  << "Getting merit function gradient ... " << endl;
      g = func.gradient(x);
      for (i=0; i<nx; ++i) xacc[i] = -1.0 * g[i] / xmass[i];
      gmagn = g.magn();
      if (debug) 
	cout << prefix_report_debug
		  << methodstring << ": "
		  << "Getting merit function value ... " << endl;
      fx = func.value();
      Ep = fx;
      Ek = 0.0;
      Et = Ep + Ek;
      





      // ###################################################################
      // ###################################################################
      // ###################################################################
      // ###################################################################
      while(true){

	// Report
	if (report_iter){
	  if (niter==0){
	    printf("%s%s: Iter %4d  Func %15.8e                 Grad %15.8e\n",
		   cond_print.prefix_report_iter.c_str(),
		   methodstring.c_str(), niter, fx, gmagn);
	    printf("Par.index  Param.  Gradient:\n");
	    for (i=0; i<nx; ++i)
	      printf("%5d  %20.10f  %20.10e\n", i, x[i], g[i]);
	  }
	  else {
	    printf("%s%s: Iter %4d  Func %15.8e  Change %15.8e  Grad %15.8e\n",
		   cond_print.prefix_report_iter.c_str(),
		   methodstring.c_str(), niter, fx, fx-fx_old, gmagn);
	    printf("Par.index  Param.  Gradient  Step:\n");
	    for (int i=0; i<nx; ++i)
	      printf("%5d  %20.10f  %20.10e  %20.10e\n", i, x[i], g[i], h[i]);
	  }
	  cout << "Reporting on parameters and data ..." << endl;
	  func.report_on_parameters_and_data();
	  cout << "End of reporting on parameters and data." << endl;
	}


	cout << "Checking for convergence" << endl;


	// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// Check for convergence
	// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	if ( check_conv(x,
			true,fx_old, true,fx, // use old and current function values
			true,g,               // use gradient
			true,h,              // use step taken
			niter,
			counters_niter,
			cond_conv, cond_debug, cond_print,
			methodstring, status) ){
	  return func.all_parameters(x);
	}

	x_old = x;
	fx_old = fx;

	cout << "Done: Checking for convergence" << endl;	


	// ******************************************************
	// Current iteration:
	// ******************************************************


	cout << "Trying to take a trial step ..." << endl;

	// *************************************************************************
	// *************************************************************************
	// Trial step

	bool step_is_ok=false;
	while ( ! step_is_ok ){


	  try {

	    // Predictor
	    for (i=0; i<nx; ++i){
	      x_trial[i] = x[i] + xvel[i] * dt + 0.5 * xacc[i] * dt*dt;
	      //if (x_trial[i] < xmin[i]) x_trial[i] = xmin[i];
	      //if (x_trial[i] > xmax[i]) x_trial[i] = xmax[i];
	      
	      xvel_trial[i] = xvel[i] + 0.5 * xacc[i] * dt;
	    }
	    for (i=0; i<nx; ++i) h[i] = x_trial[i] - x[i];

	    
	    // Check if trial point is good:
	    if (! func.point_is_good( x_trial )){ // may throw exception bad_point
	      if (report_warn)
		cout << cond_print.prefix_report_warn
		     << methodstring << ": "
		     << "Too large step " << h.magn() << ". Decreasing timestep and retrying." << endl;
	      dt = 0.33333 * dt;
	      continue;
	    }


	    // Check length of trial step:
	    if (fp_is_small( h.magn() )){
	      status.small_step = true;
	      if (report_error)
		cout << cond_print.prefix_report_error
		     << methodstring << ": "
		     << "Too small step " << h.magn() << ". Increasing timestep and retrying." << endl;
	      dt = 2.0 * dt;
	      continue;
	    }


	    step_is_ok = true;

	  }
	  catch (funcfit::bad_point & e1){
	    func.reset();
	    // Went too far. Retry with smaller step.
	    if (report_warn)
	      cout << cond_print.prefix_report_warn
		   << methodstring << ": "
		   << "Bad function value. Retrying with smaller step." << endl;
	    dt = 0.33333 * dt;
	    continue;
	  }

	  //break;
	}


	x = x_trial;
	xvel = xvel_trial;

	// Get force
	cout << "Before getting function value:" << endl;
	tv = func.free_parameters(); cout << "Free parameters" << endl;
	for (i=0; i<nx; ++i) printf("%5d  %20.10f  %20.10e\n", i, tv[i]);
	fx = func(x);
	cout << "Function value: " << fx << endl;

	cout << "Before getting gradient" << endl;
	tv = func.free_parameters(); cout << "Free parameters" << endl;
	for (i=0; i<nx; ++i) printf("%5d  %20.10f  %20.10e\n", i, tv[i]);

	g = func.gradient(x);
	cout << "After getting gradient" << endl;
	tv = func.free_parameters(); cout << "Free parameters" << endl;
	for (i=0; i<nx; ++i) printf("%5d  %20.10f  %20.10e\n", i, tv[i]);



	for (i=0; i<nx; ++i) xacc[i] = -1.0 * g[i] / xmass[i];
	gmagn = g.magn();

	// Corrector
	for (i=0; i<nx; ++i){
	  xvel[i] = xvel[i] + 0.5 * xacc[i] * dt;
	}

	Ep = fx;
	Ek = 0.0; for (i=0; i<nx; ++i) Ek += 0.5 * xmass[i] * xvel[i]*xvel[i];
	Et = Ep + Ek;
	cout << "Total energy : " << Et << endl;
	// This is a different 'h':
	for (i=0; i<nx; ++i) h[i] = x[i] - x_old[i];




	// Go to next iteration.
	niter++;
      }
      // ###################################################################
      // ###################################################################
      // ###################################################################
      // ###################################################################



      return func.all_parameters(x);
    }
    
  } ;

}





#endif

