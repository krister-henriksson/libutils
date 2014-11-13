


#ifndef FUNCFIT_LEVE_MARQ_HPP
#define FUNCFIT_LEVE_MARQ_HPP


#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <cmath>


#include "utils.hpp"
#include "utils-math.hpp"
#include "utils-vector.hpp"

#include "param.hpp"

#include "funcfit-exceptions.hpp"
#include "funcfit-basics.hpp"

using std::cout;
using std::endl;
using utils::Vector;
using utils::Matrix;
using utils::max;
using utils::abs;


namespace funcfit {

  
  template <typename T>
  class LeveMarq
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
    LeveMarq(T & funcd)
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
			    Cond_Conv  cond_conv  = Cond_Conv(),
			    Cond_Debug cond_debug = Cond_Debug(),
			    Cond_Print cond_print = Cond_Print()
			    ){
      // Input: starting point point and constructor-set functor (containing
      // gradient method)
      // Purpose: minimize function, input point will not be changed

      Vector<double> p, f, h, ag;
      Matrix<double> J, JTJ, A, A_inv;
      double fp, fp_old, gmagn, hmagn;
      double tau=1.0, gain, mu, nu=2.0;
      Vector<double> p_trial, f_trial;
      Matrix<double> J_trial;
      double fp_trial;
      string methodstring, choice;
      int n;
      int niter=0;
      double maxe=1.0, tmp;

      status = Minimization_Status();

      bool debug         = cond_debug.debug_fit_level0;
      string prefix_report_debug = cond_debug.prefix_debug_fit_level0;

      bool report_iter    = cond_print.report_iter;
      bool report_warn    = cond_print.report_warn;
      bool report_error   = cond_print.report_error;

      Counters_niter counters_niter = Counters_niter();




      p = func.free_parameters(point_in);
      Vector<double>        xmin  = func.map_vector_as_free_parameters(par_min);
      Vector<double>        xmax  = func.map_vector_as_free_parameters(par_max);
      Vector<parametertype> xtype = func.map_vector_as_free_parameters(par_type);
      n = p.size();

      methodstring = "ls-leve-marq";

      if (debug)
	cout << prefix_report_debug
		  << methodstring << ": "
		  << "Getting auxiliary data for merit function gradient ... " << endl;
      f = func.f(p);
      J = func.J(p);
      if (debug) 
	cout << prefix_report_debug
		  << methodstring << ": "
		  << "Getting merit function gradient ... " << endl;
      ag = -1.0 * func.gradient();// -1.0 * J.transpose() * f;
      gmagn = ag.magn();
      if (debug) 
	cout << prefix_report_debug
		  << methodstring << ": "
		  << "Getting merit function value ... " << endl;
      fp = func.value();//0.5 * f * f;
      JTJ = J.transpose() * J;

      for (int i=0; i!=n; ++i){
	tmp = abs(JTJ.elem(i,i));
	if (i==0 || (i>0 && tmp>maxe))
	  maxe = tmp;
      }
      mu = tau * maxe;






      // ###################################################################
      // ###################################################################
      // ###################################################################
      // ###################################################################
      while (true){


	// Report
	if (report_iter){
	  if (niter==0){
	    printf("%s%s: Iter %4d  Func %15.8e                 Grad %15.8e\n",
		   cond_print.prefix_report_iter.c_str(),
		   methodstring.c_str(), niter, fp, gmagn);
	  }
	  else {
	    printf("%s%s: Iter %4d  Func %15.8e  Change %15.8e  Grad %15.8e "
		   "Step %15.8e  %s\n", 
		   cond_print.prefix_report_iter.c_str(),
		   methodstring.c_str(), niter, fp, fp-fp_old, gmagn,
		   hmagn, choice.c_str());
	  }
	  func.report_on_parameters_and_data();
	}


	// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// Check for convergence
	// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	if (niter>0){
	  if ( check_conv(p, true,fp_old, true,fp, // use old and current function values
			  true,ag,                 // use gradient
			  true,h,                  // use step taken
			  niter,
			  counters_niter,
			  cond_conv, cond_debug, cond_print,
			  methodstring, status) ){
	    return func.all_parameters(p);
	  }
	}






	fp_old = fp;

	// *************************************************************************
	// Solve ((J^T J) + mu I) h = - J^T f for direction h.
	// *************************************************************************
	A = JTJ;
	for (int i=0; i!=n; ++i)
	  A.elem(i,i) = JTJ.elem(i,i) + mu;
	
	if ( A.solve( ag, h, A_inv) == 1 ){
	  status.singular_matrix = true;
	  status.fit_OK = false;
	  return func.all_parameters(p);
	}


	hmagn = h.magn();



	// *************************************************************************
	// *************************************************************************
	// Trial step

	while (true){

	  try {
	    // ------------------------------------------------------------------
	    while (true){
	      p_trial = p + h;
	      if (! func.point_is_good( p_trial )){
		if (report_warn)
		  cout << cond_print.prefix_report_warn
			    << methodstring << ": "
			    << "Too large step " << h.magn() << ". Halving and retrying." << endl;
		h = 0.5 * h;
		
		if (niter > cond_conv.nitermin && fp_is_small(h.magn())){
		  status.small_step = true;
		  if (report_error)
		    cout << cond_print.prefix_report_error
			      << methodstring << ": "
			      << "Too small step. Quitting." << endl;
		  return func.all_parameters(p);
		}
		continue;
	      }
	      break;
	    }


	    if (debug)
	      cout << prefix_report_debug
			<< methodstring << ": "
			<< "Getting trial auxiliary data for merit function gradient ... " << endl;
	    f_trial = func.f(p_trial);
	    // ------------------------------------------------------------------
	  }
	  catch (funcfit::bad_value & e1){
	    // Went too far. Retry with smaller step.
	    if (report_warn)
	      cout << cond_print.prefix_report_warn
			<< methodstring << ": "
			<< "Bad function value. Retrying with smaller step." << endl;
	    h = 0.5 * h;
	    continue;
	  }
	  break;
	}

	if (debug)
	  cout << prefix_report_debug
		    << methodstring << ": "
		    << "Getting trial auxiliary data for merit function gradient ... " << endl;
	J_trial = func.J(p_trial);
	if (debug)
	  cout << prefix_report_debug
		    << methodstring << ": "
		    << "Getting trial merit function value ... " << endl;
	fp_trial = func.value();//0.5 * f_trial * f_trial;

	// *************************************************************************
	// *************************************************************************



	gain = (fp - fp_trial) / (0.5 * h * (mu * h + ag));

	if (gain > 0){
	  // Go to trial point
	  p = p_trial;
	  f = f_trial;
	  J = J_trial;
	  if (debug) 
	    cout << prefix_report_debug
		      << methodstring << ": "
		      << "Getting merit function gradient ... " << endl;
	  ag = -1.0 * func.gradient();//-1.0 * J.transpose() * f;
	  gmagn = ag.magn();
	  if (debug) 
	    cout << prefix_report_debug
		      << methodstring << ": "
		      << "Getting merit function value ... " << endl;
	  fp = fp_trial;
	  //fp = func.value();//0.5 * f * f;
	  JTJ = J.transpose() * J;
	  
	  mu = mu * max( 1.0/3.0, 1.0 - (2*gain-1)*(2*gain-1)*(2*gain-1));
	  nu = 2;
	  choice = "step accepted";
	}
	else {
	  mu = mu * nu;
	  nu = 2 * nu;
	  choice = "step failed";
	}




	// Go to next iteration.
	niter++;
      }
      // ###################################################################
      // ###################################################################
      // ###################################################################
      // ###################################################################


      return func.all_parameters(p);
    }
    
  } ;

}





#endif

