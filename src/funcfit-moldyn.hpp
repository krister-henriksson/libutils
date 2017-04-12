


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


#define DT_DECR_FAC 0.7
#define DT_INCR_FAC 1.3


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
			    double min_dx,
			    double max_dx,
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

      int niterrestart = cond_conv.niterrestart;

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
      Vector<double> gn(nx,  0.0);
      double Ep=0, Ek=0, Et=0;
      double t=0, dt=1.0;


      Vector<double> tv;



      methodstring = "moldyn";




      if (debug)
	cout << prefix_report_debug
		  << methodstring << ": "
		  << "Getting merit function gradient ... " << endl;
      g = func.gradient(x);
      gn = g; gn.normalize();

      for (i=0; i<nx; ++i) xacc[i] = -1.0*gn[i] / xmass[i];
      gmagn = g.magn();
      if (debug) 
	cout << prefix_report_debug
		  << methodstring << ": "
		  << "Getting merit function value ... " << endl;
      fx = func.value();
      Ep = fx;
      Ek = 0.0;
      Et = Ep + Ek;
      
      double Ek_prev[5];
      Ek_prev[0] = 0.0;
      Ek_prev[1] = 0.0;
      Ek_prev[2] = 0.0;
      Ek_prev[3] = 0.0;
      Ek_prev[4] = Ek;



      // ###################################################################
      // ###################################################################
      // ###################################################################
      // ###################################################################
      while(true){

	if (report_iter){
	  double vb = func.value_barrier();
	  Vector<double> gb = func.gradient_barrier();
	  double gbmagn = gb.magn();

	  if (niter==0){
	    printf("%s%s: Iter %4d   Func %15.8e Func_barrier %15.8e   "
		   "Grad %15.8e Grad_barrier %15.8e\n",
		   cond_print.prefix_report_iter.c_str(),
		   methodstring.c_str(), niter, fx, vb, gmagn, gbmagn);
	    printf("Par.index  Param.  Gradient  Gradient_barrier:\n");
	    for (int i=0; i<nx; ++i)
	      printf("%5d  %20.10f  %20.10e  %20.10e\n", i, x[i], g[i], gb[i]);
	  }
	  else {
	    printf("%s%s: Iter %4d   Func %15.8e Change %15.8e   "
		   "Func_barrier %15.8e   Grad %15.8e Grad_barrier %15.8e\n",
		   cond_print.prefix_report_iter.c_str(),
		   methodstring.c_str(), niter, fx, fx-fx_old,
		   vb, gmagn, gbmagn);
	    printf("Par.index  Param.  Gradient  Gradient_barrier  Step:\n");
	    for (int i=0; i<nx; ++i)
	      printf("%5d  %20.10f  %20.10e  %20.10e  %20.10e\n", i, x[i], g[i], gb[i], h[i]);
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

	  // Adjust time step if needed:
	  while (true){
	    double step_max, step_min, td;

	    for (i=0; i<nx; ++i){
	      //xvel[i] = 0.0;

	      td = xvel[i] * dt + 0.5 * xacc[i] * dt*dt;
	      if (td < 0) td *= -1.0;
	      if (i==0 || (i>0 && td>step_max)) step_max = td;
	      if (i==0 || (i>0 && td<step_min)) step_min = td;
	    }
	    std::cout << "max_dx = " << max_dx
		      << " step_max = " << step_max
		      << " step_min = " << step_min
		      << " dt = " << dt << std::endl;

	    if      (step_max > max_dx) dt *= DT_DECR_FAC;
	    else if (step_min < min_dx) dt *= DT_INCR_FAC;
	    else break;
	  }


	  try {

	    // Predictor
	    for (i=0; i<nx; ++i){
	      //xvel[i]=0;
	      //xvel[i] = 0.0;

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
	      dt = DT_DECR_FAC * dt;
	      continue;
	    }


	    // Check length of trial step:
	    if (fp_is_small( h.magn() )){
	      status.small_step = true;
	      if (report_error)
		cout << cond_print.prefix_report_error
		     << methodstring << ": "
		     << "Too small step " << h.magn() << ". Increasing timestep and retrying." << endl;
	      dt = 1.5 * dt;
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
	    dt = DT_DECR_FAC * dt;
	    continue;
	  }

	  //break;
	}


	x = x_trial;
	xvel = xvel_trial;

	// Get force
	//cout << "Before getting function value:" << endl;
	//tv = func.free_parameters(); cout << "Free parameters" << endl;
	//for (i=0; i<nx; ++i) printf("%5d  %20.10f  %20.10e\n", i, tv[i]);
	fx = func(x);
	//cout << "Function value: " << fx << endl;

	//cout << "Before getting gradient" << endl;
	//tv = func.free_parameters(); cout << "Free parameters" << endl;
	//for (i=0; i<nx; ++i) printf("%5d  %20.10f  %20.10e\n", i, tv[i]);

	g = func.gradient(x);
	gn = g; gn.normalize();

	//cout << "After getting gradient" << endl;
	//tv = func.free_parameters(); cout << "Free parameters" << endl;
	//for (i=0; i<nx; ++i) printf("%5d  %20.10f  %20.10e\n", i, tv[i]);


	for (i=0; i<nx; ++i) xacc[i] = -1.0*gn[i] / xmass[i];
	gmagn = g.magn();

	// Corrector
	for (i=0; i<nx; ++i){
	  //xvel[i]=0;
	  //xvel[i] = 0.0;

	  xvel[i] = xvel[i] + 0.5 * xacc[i] * dt;
	}


	// If forces * momenta < 0 then we just passed a (local) minimum.
	// Remove all velocities.
	double sp = scalarproduct(xacc, xvel);
	if (sp < 0.0){
	  for (i=0; i<nx; ++i) xvel[i] = 0.0;
	}

	Ep = fx;
	Ek = 0.0; for (i=0; i<nx; ++i) Ek += 0.5 * xmass[i] * xvel[i]*xvel[i];

	Ek_prev[0] = Ek_prev[1];
	Ek_prev[1] = Ek_prev[2];
	Ek_prev[2] = Ek_prev[3];
	Ek_prev[3] = Ek_prev[4];
	Ek_prev[4] = Ek;

	if (niter >= 5){
	  double ch = 1.0;
	  if (Ek_prev[0] > 0) ch = Ek_prev[4]/Ek_prev[0];
	  if (ch > 1.0){
	    double adj = sqrt(1.0/ch);
	    for (i=0; i<nx; ++i) xvel[i] *= adj;
	    Ek *= adj;
	    Ek_prev[4] = Ek;
	    std::cout << "Adjusted kinetic energy by factor " << adj << std::endl;
	  }
	}

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

