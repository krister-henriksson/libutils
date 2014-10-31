


#ifndef FUNCFIT_GAUSS_NEWTON_HPP
#define FUNCFIT_GAUSS_NEWTON_HPP


#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>


#include "utils.hpp"
#include "utils-math.hpp"
#include "utils-vector.hpp"
#include "nr-linemethod.hpp"

#include "param.hpp"

#include "funcfit-basics.hpp"



namespace funcfit {

  
  template <typename T>
  class GaussNewton
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
    GaussNewton(T & funcd)
      :
      func(funcd)
    {}



    //#####################################################################
    // Main method
    //#####################################################################
    utils::Vector<double> minimize(utils::Vector<double> & point_in,
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



      utils::Vector<double> p, f, h, ag;
      utils::Matrix<double> J, JTJ, JTJ_inv;
      nr::Linemethod<T> linemeth(func);
      double fp, fp_old, gmagn, hmagn;
      std::string methodstring;
      int niter=0;

      status = Minimization_Status();

      bool debug         = cond_debug.debug_fit_level0;
      bool debug_linemeth= cond_debug.debug_fit_level1;
      bool debug_f1dim   = cond_debug.debug_fit_level2;
      bool debug_golden  = cond_debug.debug_fit_level3;
      bool debug_bracket = cond_debug.debug_fit_level4;
      string prefix_report_debug = cond_debug.prefix_debug_fit_level0;

      bool report_iter    = cond_print.report_iter;
      bool report_warn    = cond_print.report_warn;
      bool report_error   = cond_print.report_error;

      Counters_niter counters_niter = Counters_niter();


      p = func.free_parameters(point_in);
      Vector<double>        xmin  = func.map_vector_as_free_parameters(par_min);
      Vector<double>        xmax  = func.map_vector_as_free_parameters(par_max);
      Vector<parametertype> xtype = func.map_vector_as_free_parameters(par_type);


      methodstring = "gauss-newton";

      if (debug)
	std::cout << prefix_report_debug
		  << methodstring << ": "
		  << "Getting merit function value ... " << std::endl;
      fp = func(p);
      if (debug) std::cout << "Getting auxiliary data for merit function gradient ... " << std::endl;
      f = func.f(p);
      J = func.J(p);
      if (debug) 
	std::cout << prefix_report_debug
		  << methodstring << ": "
		  << "Getting merit function gradient ... " << std::endl;
      ag = -1.0 * func.gradient();// -1.0 * J.transpose() * f;
      gmagn = ag.magn();




      // ###################################################################
      // ###################################################################
      // ###################################################################
      // ###################################################################
      while(true){

	if (niter==0){
	  if (report_iter){
	    printf("%s%s: Iter %4d  Func %15.8e                 Grad %15.8e\n",
		   cond_print.prefix_report_iter.c_str(),
		   methodstring.c_str(), niter, fp, gmagn);
	    func.report_on_parameters_and_data();
	  }
	}
	else {
	  if (report_iter){
	    printf("%s%s: Iter %4d  Func %15.8e  Change %15.8e  Grad %15.8e "
		   "Step %15.8e\n", 
		   cond_print.prefix_report_iter.c_str(),
		   methodstring.c_str(), niter, fp, fp-fp_old, gmagn,
		   hmagn);
	    func.report_on_parameters_and_data();
	  }
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
	// Solve (J^T J) h = - J^T f for direction h.
	// *************************************************************************

	JTJ = J.transpose() * J;
	if ( JTJ.solve( ag, h, JTJ_inv) == 1 ){
	  status.singular_matrix = true;
	  status.fit_OK = false;
	  return func.all_parameters(p);
	}




	// *************************************************************************
	// Minimize function aint direction h, starting at point p
	// Updates p, changes h.
	// *************************************************************************
	if (debug) 
	  std::cout << prefix_report_debug
		    << methodstring << ": "
		    << "Performing line minimization ... " << std::endl;
	fp = linemeth.linmin(p, xmin, xmax, xtype, h, cond_conv.steptolrel,
			     debug_linemeth, debug_f1dim, debug_golden, debug_bracket);

	if (niter > cond_conv.nitermin && linemeth.status.small_step){
	  status.small_step = true;
	  if (report_error)
	    std::cout << cond_print.prefix_report_error
		      << methodstring << ": "
		      << "Too small step " << h.magn() << ". Quitting." << std::endl;
	  return func.all_parameters(p);
	}

	if (debug) 
	  std::cout << prefix_report_debug
		    << methodstring << ": "
		    << "Getting auxiliary data for merit function gradient ... " << std::endl;
	f = func.f(p);
	J = func.J(p);
	if (debug)
	  std::cout << prefix_report_debug
		    << methodstring << ": "
		    << "Getting merit function gradient ... " << std::endl;
	ag = -1.0 * func.gradient();// -1.0 * J.transpose() * f;
	gmagn = ag.magn();

	/*
	status.funcmin = fp;
	status.gradmin = gmagn;
	status.niter   = niter;
	*/

	// Step length:
	hmagn = h.magn();




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

