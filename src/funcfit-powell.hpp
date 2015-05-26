


#ifndef FUNCFIT_POWELL_HPP
#define FUNCFIT_POWELL_HPP


#include <iostream>
#include <limits>

#include <cstdlib>
#include <cmath>


#include "utils.hpp"
#include "utils-math.hpp"
#include "utils-vector.hpp"

#include "param.hpp"

#include "nr-linemethod.hpp"

#include "funcfit-basics.hpp"

using namespace utils;
using std::numeric_limits;
using std::cout;
using std::endl;



namespace funcfit {

  
  template <typename T>
  class Powell
  {


    //#####################################################################
    // Member variables
    //#####################################################################
  public:
    Vector<double> p;
    Vector<double> xi;
    T & func;
    Minimization_Status status;
    double fret;


    //#####################################################################
    // Constructor
    //#####################################################################
    // NB:  Use constructor to set functor. Point and direction are set by
    // NB:  minimize() method.
    Powell(T & funcd)
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

      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      p = func.free_parameters(point_in);
      Vector<double>        xmin  = func.map_vector_as_free_parameters(par_min);
      Vector<double>        xmax  = func.map_vector_as_free_parameters(par_max);
      Vector<parametertype> xtype = func.map_vector_as_free_parameters(par_type);
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      //double eps = numeric_limits<double>::epsilon();
      std::string methodstring = "powell";
      double fp, fp_old, fptt, td1,td2,td3;
      int n = p.size(), niter=0, i,j;
      nr::Linemethod<T> linemeth(func);
      Vector<double> pt(n, 0), ptt(n, 0);
      int ibig;
      double del;
      double hmagn;

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

      int niterrestart = cond_conv.niterrestart;

      /*
      std::cout << "func.value(p) = " << func.value(p) << std::endl;
      std::cout << "func.value( ) = " << func.value( ) << std::endl;
      std::cout << "func.gradient(p) = " << func.gradient(p) << std::endl;
      std::cout << "func.gradient() = " << func.gradient() << std::endl;
      */


      Vector< Vector<double> > ximat(n, Vector<double>(n,0) );
      for (i=0; i<n; ++i) ximat[i][i]=1;

      xi.resize(n);
      if (debug)
	cout << prefix_report_debug
	     << methodstring << ": "
	     << "Getting merit function value ... " << endl;
      fret = fp = func(p);

      for (i=0; i<n; ++i) pt[i] = p[i];



      ibig = 0;
      del = 0; // largest function decrease



      // ###################################################################
      // ###################################################################
      // ###################################################################
      // ###################################################################
      while (true){




	// *******************************************************************
	// Report
	// *******************************************************************
	if (report_iter){
	  if (niter==0){
	    printf("%s%s: Iter %4d  Func %15.8e\n",
		   cond_print.prefix_report_iter.c_str(),
		   methodstring.c_str(), niter, fp);
	  }
	  else {
	    printf("%s%s: Iter %4d  Func %15.8e  Change %15.8e  Step %15.8e\n",
		   cond_print.prefix_report_iter.c_str(),
		   methodstring.c_str(), niter, fp, fp_old-fp, hmagn);
	  }
	  func.report_on_parameters_and_data();
	}

	// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// Check for convergence
	// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	if ( check_conv(p, true,fp_old, true,fp, // use old and current function values
			false,Vector<double>(xi.size(),0), // use gradient
			true,xi,                  // use step taken
			niter,
			counters_niter,
			cond_conv, cond_debug, cond_print,
			methodstring, status) ){
	  return func.all_parameters(p);
	}


	if (niter > cond_conv.nitermin && linemeth.status.small_step){
	  status.small_step = true;
	  if (report_error)
	    cout << cond_print.prefix_report_error
		 << methodstring << ": "
		 << "Too small step " << xi.magn() << ". Quitting." << endl;
	  return func.all_parameters(p);
	}




	fp_old = fp;


	for (i=0; i<n; ++i){
	  // Line minimization in direction i:
	  for (j=0; j<n; ++j) xi[j] = ximat[j][i];
	  if (debug)
	    cout << prefix_report_debug
		 << methodstring << ": "
		 << "Performing line minimization ... " << endl;
	  fptt = fret;
	  fret = linemeth.linmin(p, xmin, xmax, xtype, xi, cond_conv.steptolrel,
				 debug_linemeth, debug_f1dim, debug_golden, debug_bracket);

	  // Step length:
	  hmagn = xi.magn();

	  if (fptt - fret > del){
	    del = fptt - fret;
	    ibig = i+1;
	  }
	}







	for (j=0; j<n; ++j){
	  // Extrapolated point:
	  ptt[j] = 2*p[j] - pt[j];
	  // "Average" direction moved:
	  xi[j] = p[j] - pt[j];
	  // Old starting point:
	  pt[j] = p[j];
	}
	// Function value at extrapolated point:
	if (debug)
	  cout << prefix_report_debug
	       << methodstring << ": "
	       << "Getting merit function value ... " << endl;
	fptt = func(ptt);

	if (fptt < fp){
	  td1 = fp-fret-del;
	  td2 = fp-fptt;
	  td3 = 2.0*(fp - 2*fret + fptt)*td1*td1 - del*td2*td2;
	  if (td3 < 0.0){
	    if (debug)
	      cout << prefix_report_debug
		   << methodstring << ": "
		   << "Performing line minimization ... " << endl;
	    fret = linemeth.linmin(p, xmin, xmax, xtype, xi, cond_conv.steptolrel, debug_linemeth);
	    // Step length:
	    hmagn = xi.magn();

	    for (j=0; j<n; ++j){
	      ximat[j][ibig-1] = ximat[j][n-1];
	      ximat[j][n-1] = xi[j];
	    }
	  }
	}


	fp = fret;

	niter++;
      } // end of loop over iterations
      // ###################################################################
      // ###################################################################
      // ###################################################################
      // ###################################################################



      return func.all_parameters(p);
      
      

    } // end of method

  } ; // end of class


}

#endif








