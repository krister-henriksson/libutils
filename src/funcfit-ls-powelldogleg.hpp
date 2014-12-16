


#ifndef FUNCFIT_POWELLDOGLEG_HPP
#define FUNCFIT_POWELLDOGLEG_HPP


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

using std::cout;
using std::endl;
using utils::Vector;


namespace funcfit {

  
  template <typename T>
  class PowellDogLeg
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
    PowellDogLeg(T & funcd)
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
			    double rt,
			    double rt_min,
			    Cond_Conv  cond_conv  = Cond_Conv(),
			    Cond_Debug cond_debug = Cond_Debug(),
			    Cond_Print cond_print = Cond_Print()
			    ){
      // Input: starting point point and constructor-set functor (containing
      // gradient method)
      // Purpose: minimize function, input point will not be changed

      Vector<double> p, f, af, ag;
      Vector<double> h_gn, h_sd, h, a, b, tmpv1;
      Vector<double> p_trial, f_trial;
      Matrix<double> J, JTJ, JTJ_inv;
      double c, a_sq, ba_sq, tmp1, tmp2;
      double fp, fp_old, gmagn, hmagn;
      double alpha, beta, gain;
      double fp_trial;
      string methodstring, stepmeth, choice;
      int niter=0;

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

      /*
      //cout.precision(15);
      cout.precision(10);
      cout.setf (ios_base::scientific , ios_base::floatfield);
      cout.width(20);
      cout << "cout settings changed" << endl;
      */

      methodstring = "ls-powelldogleg";

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
      ag = -1.0 * func.gradient();//-1.0 * J.transpose() * f;
      gmagn = ag.magn();
      if (debug) 
	cout << prefix_report_debug
		  << methodstring << ": "
		  << "Getting merit function value ... " << endl;
      fp = func.value();//0.5 * f * f;








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
		   methodstring.c_str(), niter, fp, gmagn);
	    printf("Par.index  Param.  Gradient:\n");
	    for (int i=0; i<ag.size(); ++i)
	      printf("%5d  %20.10f  %20.10e\n", i, p[i], -1.0*ag[i]);
	  }
	  else {
	    printf("%s%s: Iter %4d  Func %15.8e  Change %15.8e  Grad %15.8e "
		   "Step %15.8e  Radius %15.8e   %s %s\n", 
		   cond_print.prefix_report_iter.c_str(),
		   methodstring.c_str(), niter, fp, fp-fp_old, gmagn,
		   hmagn, rt, stepmeth.c_str(), choice.c_str());
	    //cout << "Chi^2 components: " << func.f() << endl;
	    //cout << "DataY: " << func.DataY() << endl;
	    //cout << "ModelDataY: " << func.ModelDataY() << endl;
	    printf("Par.index  Param.  Gradient  Step:\n");
	    for (int i=0; i<ag.size(); ++i)
	      printf("%5d  %20.10f  %20.10e  %20.10e\n", i, p[i], -1.0*ag[i], h[i]);
	  }
	  func.report_on_parameters_and_data();
	}

	// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// Check for convergence
	// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	if ( check_conv(p,
			true,fp_old, true,fp, // use old and current function values
			true,ag,              // use gradient
			true,h,               // use step taken
			niter,
			counters_niter,
			cond_conv, cond_debug, cond_print,
			methodstring, status) ){
	  return func.all_parameters(p);
	}

	if (rt < rt_min){
	  if (report_warn)
	    cout << cond_print.prefix_report_warn
		      << methodstring << ": "
		      << "Trust region radius " << rt << " "
		      << "has shrunk below smallest allowed value " << rt_min
		      << endl;
	  status.trustradius_too_small = true;
	  return func.all_parameters(p);
	}






	fp_old = fp;

	
	tmp1  = ag * ag;
	tmpv1 = J * ag;
	tmp2 = tmpv1 * tmpv1;
	alpha = tmp1 / tmp2;


	// *************************************************************************
	// Gauss-Newton direction
	// *************************************************************************

	JTJ = J.transpose() * J;
	if ( JTJ.solve( ag, h_gn, JTJ_inv) == 1 ){
	  status.singular_matrix = true;
	  status.fit_OK = false;
	  return func.all_parameters(p);
	}


	// af = -1.0 * f;
	// J.solve( af, h_gn, JTJ_inv);


 	// *************************************************************************
	// Steepest descent direction
	// *************************************************************************
	h_sd = alpha * ag;


	// *************************************************************************
	// Dog-leg step calculation
	// *************************************************************************
	if (h_gn.magn() < rt){
	  h = h_gn;
	  stepmeth = "GN";
	}
	else if (h_sd.magn() >= rt){
	  h = rt * (h_sd / h_sd.magn());
	  stepmeth = "SD";
	}
	else {
	  a = h_sd;
	  b = h_gn;
	  c = a * (b-a);
	  ba_sq = (b-a)*(b-a);
	  a_sq = a*a;
	  tmp1 = rt*rt - a_sq;
	  tmp2 = sqrt(c*c + ba_sq * tmp1);
	  if (c <= 0){
	    beta = (-c + tmp2)/ba_sq;
	  }
	  else {
	    beta = tmp1 / (c + tmp2);
	  }
	  h = h_sd + beta * (h_gn - h_sd);
	  stepmeth = "DL";
	}




	// *************************************************************************
	// *************************************************************************
	// Trial step

	while (true){
	  int nerr=0;


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
		
		if (fp_is_small(h.magn())){
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

	    f_trial = func.f(p_trial);
	    // ------------------------------------------------------------------
	  }
	  catch (funcfit::bad_point & e1){
	    func.reset();
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
		    << "Getting merit function value ... " << endl;
	fp_trial = func.value();//0.5 * f_trial * f_trial;

	// *************************************************************************
	// *************************************************************************










	tmp1 = fp - fp_trial;
	if (stepmeth=="GN")
	  tmp2 = fp;
	else if (stepmeth=="SD")
	  tmp2 = rt * (2.0*h_sd*h_sd - rt)/(2.0*alpha);
	else
	  tmp2 = 0.5 * alpha * (1-beta)*(1-beta)*ag*ag
	    + beta * (2-beta) * fp;
	
	gain = tmp1 / tmp2;


	if (gain > 0){
	  // Accept step
	  p = p_trial;
	  f = f_trial;
	  fp = fp_trial;

	  if (debug) 
	    cout << prefix_report_debug
		      << methodstring << ": "
		      << "Getting auxiliary data for merit function gradient ... " << endl;	  
	  J = func.J(p);
	  if (debug)
	    cout << prefix_report_debug
		      << methodstring << ": "
		      << "Getting function gradient ... " << endl;
	  ag = -1.0 * func.gradient();//-1.0 * J.transpose() * f;
	  gmagn = ag.magn();
	  
	  choice = "step accepted";
	}
	else {
	  choice = "step failed";
	}	
	
	if (gain > 0.75)
	  rt = max(rt, 3.0 * h.magn());
	else if (gain < 0.25)
	  rt = 0.5 * rt;

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

