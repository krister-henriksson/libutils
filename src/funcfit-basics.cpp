




#include <limits>
#include <string>

#include <cfloat>
#include <cmath>

#include "funcfit-basics.hpp"


using std::string;

using namespace funcfit;
using utils::Vector;



bool funcfit::check_conv(const Vector<double> x,
			 const bool use_fxold, const double fxold,
			 const bool use_fx,    const double fx,
			 const bool use_grad,  const Vector<double> & grad,
			 const bool use_dx,    const Vector<double> & dx,
			 const int niter,
			 Counters_niter & counters_niter,
			 const Cond_Conv  & cond_conv,
			 const Cond_Debug & cond_debug,
			 const Cond_Print & cond_print,
			 string methodstring,
			 Minimization_Status & status
			 ){

  double FUNCTOLABS = cond_conv.functolabs;
  double STEPTOLABS = cond_conv.steptolabs;
  double GRADTOLABS = cond_conv.gradtolabs;
  double FUNCTOLREL = cond_conv.functolrel;
  double STEPTOLREL = cond_conv.steptolrel;

  int NITERMIN_FUNCTOLABS = cond_conv.nitermin_functolabs;
  int NITERMIN_STEPTOLABS = cond_conv.nitermin_steptolabs;
  int NITERMIN_GRADTOLABS = cond_conv.nitermin_gradtolabs;
  int NITERMIN_FUNCTOLREL = cond_conv.nitermin_functolrel;
  int NITERMIN_STEPTOLREL = cond_conv.nitermin_steptolrel;

  int NITERMIN       = cond_conv.nitermin;
  int NITERMAX       = cond_conv.nitermax;
  int NITERMAX_SAMEFUNCVAL = cond_conv.nitermax_samefuncval;

  if (niter <= NITERMIN) return false;


  /*
  bool debug          = cond_print.debug_fit_level0;
  bool debug_linemeth = cond_print.debug_fit_level1;
  bool report_iter    = cond_print.report_iter;
  */
  bool report_warn    = cond_print.report_warn;
  bool report_error   = cond_print.report_error;

  bool report_conv    = cond_conv.report_conv;
    
  double eps = std::numeric_limits<double>::epsilon();

  status.niter = niter;




  // ###############################################################
  // Calculate indicators
  // ###############################################################

  double absfxold, absfx, dfx, dfxtol;
  if (use_fx && use_fxold){
    absfxold = fxold;      if (fxold<0.0) absfxold *= -1.0;
    absfx    = fx;         if (fx   <0.0) absfx    *= -1.0;
    dfx      = fxold - fx; if (dfx  <0.0) dfx *= -1.0;
    dfxtol   = dfx / (0.5*(absfxold + absfx + eps));
  }

  int i,j;
  double absxi, absdxi, ri;

  // ***************************************************************
  bool steptolabs_ok=false;
  if (STEPTOLABS > 0.0 && use_dx){
    j=0;

    for (i=0; i<x.size(); ++i){
      absxi  = dx[i];  if (absxi <0.0) absxi  *= -1.0;
      if (absxi < STEPTOLABS) j++;
    }
    if (j==x.size()) steptolabs_ok=true;
  }
  cout << "steptolabs , j = " << j << endl;

  // ***************************************************************
  bool steptolrel_ok=false;
  if (STEPTOLREL > 0.0 && use_dx){
    j=0;

    for (i=0; i<x.size(); ++i){
      absxi  = x[i];  if (absxi <0.0) absxi  *= -1.0;
      absdxi = dx[i]; if (absdxi<0.0) absdxi *= -1.0;
      if (absxi>eps) ri = absdxi / absxi;
      else           ri = absdxi;
      if (ri < STEPTOLREL) j++;
    }
    if (j==x.size()) steptolrel_ok=true;
  }

  // ***************************************************************
  double absgi;
  bool gradtolabs_ok=false;
  if (GRADTOLABS > 0.0 && use_grad){
    j=0;
    for (i=0; i<grad.size(); ++i){
      absgi = grad[i]; if (absgi<0.0) absgi *= -1.0;
      if (absgi < GRADTOLABS) j++;
    }
    if (j==grad.size()) gradtolabs_ok=true;
  }


  // ###############################################################
  // Counters
  // ###############################################################

  if (use_fx){
    status.funcmin = fx;
    if (FUNCTOLABS > 0.0 &&
	absfx < FUNCTOLABS) counters_niter.niter_functolabs++;
    else counters_niter.niter_functolabs=0;
  }


  if (use_fx && use_fxold){
    if (FUNCTOLREL > 0.0 &&
	dfxtol < FUNCTOLREL) counters_niter.niter_functolrel++;
    else counters_niter.niter_functolrel=0;
  }


  if (use_grad){
    status.gradmin = grad.magn();
    if (GRADTOLABS > 0.0 &&
	gradtolabs_ok) counters_niter.niter_gradtolabs++;
    else counters_niter.niter_gradtolabs=0;
  }


  if (use_dx){
    if (STEPTOLABS > 0.0 &&
	steptolabs_ok) counters_niter.niter_steptolabs++;
    else counters_niter.niter_steptolabs=0;
  }


  if (use_dx){
    if (STEPTOLREL > 0.0 &&
	steptolrel_ok) counters_niter.niter_steptolrel++;
    else counters_niter.niter_steptolrel=0;
  }


  if (use_fx && use_fxold &&
      ( ( fxold*(1+eps)>fx && fxold*(1-eps)<fx )
	||
	( fx*(1+eps)>fxold && fx*(1-eps)<fxold ) ) )
    counters_niter.niter_samefuncval++;
  else
    counters_niter.niter_samefuncval = 0;



  // ###############################################################
  // Check indicators
  // ###############################################################

  if (FUNCTOLABS > 0.0 && use_fx &&
      counters_niter.niter_functolabs >= NITERMIN_FUNCTOLABS &&
      absfx < FUNCTOLABS){
    if (report_conv)
      cout << cond_conv.prefix_report_conv
	   << methodstring << ": "
	   << "Absolute function value " << absfx << " "
	   << "has reached convergence within limit " << FUNCTOLABS
	   << endl;
    status.functolabs_reached = true;
    status.fit_OK = true;
    return true;
  }


  if (FUNCTOLREL > 0.0 && use_fx && use_fxold &&
      counters_niter.niter_functolrel >= NITERMIN_FUNCTOLREL &&
      dfxtol < FUNCTOLREL){
    if (report_conv)
      cout << cond_conv.prefix_report_conv
	   << methodstring << ": "
	   << "Relative change " << dfxtol << " in function value "
	   << "has reached convergence within limit " << FUNCTOLREL
	   << endl;
    status.functolrel_reached = true;
    status.fit_OK = true;
    return true;
  }



  cout << "STEPTOLABS " << STEPTOLABS << endl;
  cout << "counters_niter.niter_steptolabs " << counters_niter.niter_steptolabs << endl;
  cout << "steptolabs_ok " << steptolabs_ok << endl;
  cout << "use_dx " << use_dx << endl;

  if (STEPTOLABS > 0.0 && use_dx &&
      counters_niter.niter_steptolabs >= NITERMIN_STEPTOLABS &&
      steptolabs_ok){
    if (report_conv)
      cout << cond_conv.prefix_report_conv
	   << methodstring << ": "
	   << "All absolute lengths of step coordinates (total " << dx.magn() << ") "
	   << "have reached convergence within limit " << STEPTOLABS
	   << endl;
    status.steptolabs_reached = true;
    status.fit_OK = true;
    return true;
  }


  if (STEPTOLREL > 0.0 && use_dx &&
      counters_niter.niter_steptolrel >= NITERMIN_STEPTOLREL &&
      steptolrel_ok){
    if (report_conv)
      cout << cond_conv.prefix_report_conv
	   << methodstring << ": "
	   << "All relative lengths of step coordinates (total " << dx.magn() << ") "
	   << "have reached convergence within limit " << STEPTOLREL
	   << endl;
    status.steptolrel_reached = true;
    status.fit_OK = true;
    return true;
  }


  if (GRADTOLABS > 0.0 && use_grad &&
      counters_niter.niter_gradtolabs >= NITERMIN_GRADTOLABS &&
      gradtolabs_ok){
    if (report_conv)
      cout << cond_conv.prefix_report_conv
	   << methodstring << ": "
	   << "Absolute gradient value " << grad.magn() << " "
	   << "has reached convergence within limit " << GRADTOLABS
	   << endl;
    status.gradtolabs_reached = true;
    status.fit_OK = true;
    return true;
  }



  cout << "NITERMAX_SAMEFUNCVAL" << NITERMAX_SAMEFUNCVAL << endl;
  cout << "counter_niter.niter_samefuncval " << counters_niter.niter_samefuncval << endl;
  if (NITERMAX_SAMEFUNCVAL > 0 &&
      counters_niter.niter_samefuncval >= NITERMAX_SAMEFUNCVAL){
    if (report_warn)
      cout << cond_print.prefix_report_warn
	   << methodstring << ": "
	   << "Too many iterations with same min/max function values. Exiting." << endl;
    status.nitermax_reached = true;
    return true;
  }


  if (NITERMAX > 0 && 
      niter >= NITERMAX){
    if (report_warn)
      cout << cond_print.prefix_report_warn
	   << methodstring << ": "
	   << "Maximum number " << NITERMAX << " of iterations has been reached. Exiting." << endl;
    status.nitermax_reached = true;
    return true;
  }




  return false;
}


