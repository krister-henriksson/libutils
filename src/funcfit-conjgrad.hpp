


#ifndef FUNCFIT_CONJGRAD_HPP
#define FUNCFIT_CONJGRAD_HPP


#include <iostream>
#include <cstdlib>
#include <cmath>


#include "utils.hpp"
#include "utils-math.hpp"
#include "utils-vector.hpp"

#include "param.hpp"

#include "nr-linemethod.hpp"

#include "funcfit-basics.hpp"

using std::cout;
using std::endl;
using std::string;
using utils::abs;
using utils::Vector;



///////////////////////////////////////////////////////////
// Usage:
///////////////////////////////////////////////////////////
//
//   Funcd fd;
//
//   FuncFitCG<Funcd> cgf(fd);
//
//   xmin = cgf.minimize(...);
//   fmin = cgf.status.funcmin;
//
///////////////////////////////////////////////////////////
//
//  class Funcd {
//    // * Return function value:
//    double operator()(void);
//    double operator()(utils::Vector<double> & x);
//    // * Return gradient:
//    utils::Vector<double> gradient(const utils::Vector<double> & x);
//    // * Return free parameters:
//    utils::Vector<double> free_parameters(void);
//    utils::Vector<double> free_parameters(const utils::Vector<double> & x);
//    // * Return all parameters:
//    utils::Vector<double> all_parameters(void);
//    utils::Vector<double> all_parameters(const utils::Vector<double> & x);
//    // * Check parameters:
//    bool bad_point(void);
//    bool bad_point(const utils::Vector<double> & x);
//
//  }
//
///////////////////////////////////////////////////////////


namespace funcfit {

  
  template <typename T>
  class ConjGrad
  {


    //#####################################################################
    // Member variables
    //#####################################################################
  public:
    Vector<double> p;  // needed by Linemethod
    Vector<double> xi; // needed by Linemethod
    T & func;
    Minimization_Status status;


    //#####################################################################
    // Constructor
    //#####################################################################
    // NB:  Use constructor to set functor. Point and direction are set by
    // NB:  minimize() method.
    ConjGrad(T & funcd)
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
			    Cond_Print cond_print = Cond_Print(),
			    bool ver=2
			    ){
      // Input: starting point point and constructor-set functor (containing
      // gradient method)
      // Purpose: minimize function, input point will not be changed



      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      p = func.free_parameters(point_in);
      Vector<double>        xmin  = func.map_vector_as_free_parameters(par_min);
      Vector<double>        xmax  = func.map_vector_as_free_parameters(par_max);
      Vector<parametertype> xtype = func.map_vector_as_free_parameters(par_type);
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      string methodstring = "conjgrad";
      double gmagn, hmagn, fp_old;
      int niter=0;
      nr::Linemethod<T> linemeth(func);
      double gg, dgg, gam, fp;
      int n = p.size();
      Vector<double> g(n,0), h(n,0), dx(n,0);

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

      Vector<double> p_orig;
      double fp_orig;


      /*
      cout << "func.value(p) = " << func.value(p) << endl;
      cout << "func.value( ) = " << func.value( ) << endl;
      cout << "func.gradient(p) = " << func.gradient(p) << endl;
      cout << "func.gradient() = " << func.gradient() << endl;
      */



      xi.resize(n);
      if (debug)
	cout << prefix_report_debug
		  << methodstring << ": "
		  << "Getting merit function value ... " << endl;
      fp = func(p);
      if (debug)
	cout << prefix_report_debug
		  << methodstring << ": "
		  << "Getting merit function gradient ... " << endl;
      xi = func.gradient();
      gmagn = xi.magn();

      for (int j=0; j!=n; ++j){
	g[j] = -xi[j];
	xi[j] = h[j] = g[j];
	dx[j] = 0.0;
      }
      // Step length:
      hmagn = dx.magn();



      // ###################################################################
      // ###################################################################
      // ###################################################################
      // ###################################################################
      while (true){



	if (report_iter){
	  if (niter==0){
	    printf("%s%s: Iter %4d  Func %15.8e                 Grad %15.8e\n",
		   cond_print.prefix_report_iter.c_str(),
		   methodstring.c_str(), niter, fp, gmagn);
	  }
	  else {
	    printf("%s%s: Iter %4d  Func %15.8e  Change %15.8e  Grad %15.8e "
		   "Step %15.8e\n",
		   cond_print.prefix_report_iter.c_str(),
		   methodstring.c_str(), niter, fp, fp-fp_old, gmagn,
		   hmagn);
	  }
	  func.report_on_parameters_and_data();
	}
	/*
	status.funcmin = fp;
	status.gradmin = gmagn;
	status.niter   = niter;
	*/

	// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// Check for convergence
	// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	if (niter>0){
	  if ( check_conv(p, true,fp_old, true,fp, // use old and current function values
			  true,g,                  // use gradient
			  true,dx,                 // use step taken
			  niter,
			  counters_niter,
			  cond_conv, cond_debug, cond_print, methodstring, status) ){
	    return func.all_parameters(p);
	  }
	}
	

	fp_old = fp;

	p_orig  = p;
	fp_orig = fp;



	
	// *************************************************************************
	// Minimize function aint direction xi, starting at point p
	// - uses p, xi, funcd
	// *************************************************************************
	if (debug)
	  cout << prefix_report_debug
		    << methodstring << ": "
		    << "Performing line minimization ... " << endl;

	//func.backup();
	fp = linemeth.linmin(p, xmin, xmax, xtype, xi, cond_conv.steptolrel,
			     debug_linemeth, debug_f1dim, debug_golden, debug_bracket);
	dx = xi;
	// Step length:
	hmagn = dx.magn();


	if (niter > cond_conv.nitermin && linemeth.status.small_step){
	  status.small_step = true;
	  if (report_warn)
	    cout << cond_print.prefix_report_error
		      << methodstring << ": "
		      << "Too small step " << xi.magn() << ". Quitting." << endl;
	  //func.restore();
	  p  = p_orig;
	  fp = fp_orig;
	  return func.all_parameters(p);
	}




	if (debug) 
	  cout << prefix_report_debug
		    << methodstring << ": "
		    << "Getting merit function gradient ... " << endl;
	xi = func.gradient(p);
	gmagn = xi.magn();





	dgg = gg = 0.0;
	if (ver==1){
	  for (int j=0; j!=n; ++j){
	    gg += g[j]*g[j];
	    dgg += xi[j] * xi[j]; // Fletcher-Reeves
	  }
	}
	else {
	  for (int j=0; j!=n; ++j){
	    gg += g[j]*g[j];
	    dgg += (xi[j] + g[j]) * xi[j]; // Polak-Ribiere
	  }
	}

	gam = dgg/gg;
	for (int j=0; j!=n; ++j){
	  g[j] = -xi[j];
	  xi[j] = h[j] = g[j] + gam*h[j];
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

