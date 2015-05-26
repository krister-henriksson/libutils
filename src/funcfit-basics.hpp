


#ifndef FUNCFIT_BASICS_HPP
#define FUNCFIT_BASICS_HPP


#include <limits>
#include <string>

#include <cfloat>
#include <cmath>


#include "utils-vector.hpp"



using std::string;
using utils::Vector;



/* ###########################################################################

Functions to be fitted are given as type

  T func

to the fitting object constructor. The 'func' functor is required
to have at least the following public methods:

()
(p)
gradient()
gradient(p)
free_parameters()
free_parameters(q)
all_parameters()
all_parameters(q)
bad_point()
bad_point(q)
map_vector_as_free_parameters(s)

Here p is vector containing values of free parameters only, and q
can contain free and constrained parameters.
The vector s can be any vector of the same type as the parameters (double).

########################################################################### */


namespace funcfit {


  // ####################################################################
  // ####################################################################
  // ####################################################################
  // CLASSES
  // ####################################################################
  // ####################################################################
  // ####################################################################


  enum Conv_Status { not_Conv,
		     Conv_functolabs, Conv_steptolabs, Conv_gradtolabs,
		     Conv_functolrel, Conv_steptolrel };


  // ####################################################################
  // ####################################################################
  // ####################################################################


  class Cond_Conv {
  public:
    double functolabs;
    double steptolabs;
    double gradtolabs;
    double functolrel;
    double steptolrel;

    int nitermin_functolabs;
    int nitermin_steptolabs;
    int nitermin_gradtolabs;
    int nitermin_functolrel;
    int nitermin_steptolrel;

    int nitermin;
    int nitermax;
    int niterrestart;
    int nitermax_samefuncval;

    bool   report_conv; // report convergence type?
    string prefix_report_conv;


    Cond_Conv(){
      double small = sqrt( std::numeric_limits<double>::epsilon() );
      functolabs =  small;
      functolrel =  100*small;
      gradtolabs =  100*small;
      steptolabs =  small;
      steptolrel =  100*small;

      nitermin_functolabs = 5;
      nitermin_steptolabs = 5;
      nitermin_gradtolabs = 5;
      nitermin_functolrel = 5;
      nitermin_steptolrel = 5;

      nitermin = 3;
      nitermax = 200;
      niterrestart = -50;
      nitermax_samefuncval = 10;

      report_conv = false;
      prefix_report_conv  = "";
    }
  } ;

  // ####################################################################
  // ####################################################################
  // ####################################################################


  class Cond_Debug {
  public:
    bool debug_fit_level0;
    bool debug_fit_level1;
    bool debug_fit_level2;
    bool debug_fit_level3;
    bool debug_fit_level4;

    string prefix_debug_fit_level0;
    string prefix_debug_fit_level1;
    string prefix_debug_fit_level2;
    string prefix_debug_fit_level3;
    string prefix_debug_fit_level4;


    Cond_Debug(){
      debug_fit_level0 = false;
      debug_fit_level1 = false;
      debug_fit_level2 = false;
      debug_fit_level3 = false;
      debug_fit_level4 = false;

      prefix_debug_fit_level0 = "debug_fit_level0: ";
      prefix_debug_fit_level1 = "debug_fit_level1: ";
      prefix_debug_fit_level2 = "debug_fit_level2: ";
      prefix_debug_fit_level3 = "debug_fit_level3: ";
      prefix_debug_fit_level4 = "debug_fit_level4: ";
    }
  } ;



  // ####################################################################
  // ####################################################################
  // ####################################################################


  class Cond_Print {
  public:
    bool report_iter; // report current pogress of minimization process?
    bool report_warn; // report warnings?
    bool report_error; // report errors?

    string prefix_report_iter;
    string prefix_report_warn;
    string prefix_report_error;

    bool report_iter_params; // report current parameter values?
    int  report_iter_params_nmax; // report if # of params is LE (<=) this value


    Cond_Print(){
      report_iter_params = false;
      report_warn = false;
      report_error = true;

      prefix_report_iter  = "";
      prefix_report_warn  = "";
      prefix_report_error = "";

      report_iter_params_nmax = 3;
      report_iter = false;


    }
  } ;

  // ####################################################################
  // ####################################################################
  // ####################################################################


  class Cond_Exec {
  public:
    // Q: Perform auxiliary calculations (derivatives, Hessian, ...) that do not
    // contribute to the value of the merit function?
    // A: Always do this for normal operation.
    bool exec_aux;


    Cond_Exec(){
      exec_aux = true;
    }

  } ;


  // ####################################################################
  // ####################################################################
  // ####################################################################


  class Counters_niter {
  public:
    int niter_functolabs;
    int niter_steptolabs;
    int niter_gradtolabs;
    int niter_functolrel;
    int niter_steptolrel;
    int niter_samefuncval;

    Counters_niter(){
      niter_functolabs = 0;
      niter_steptolabs = 0;
      niter_gradtolabs = 0;
      niter_functolrel = 0;
      niter_steptolrel = 0;
      niter_samefuncval = 0;
    }
  } ;



  // ####################################################################
  // ####################################################################
  // ####################################################################


  // Member of minimization class
  // Diagnose a minimization attempt
  class Minimization_Status {
  public:
    double funcmin;
    double gradmin;

    int   niter;
    bool   fit_OK;

    bool functolabs_reached;
    bool steptolabs_reached;
    bool gradtolabs_reached;
    bool functolrel_reached;
    bool steptolrel_reached;

    bool nitermax_reached;
    bool trustradius_too_small;
    bool small_step;
    bool singular_matrix;
    bool bad_point; // outside parameter limits
    bool bad_value; // e.g. attempting to calculate sqrt(-1.0)

    
    // constructor
    Minimization_Status()
      :
      funcmin( std::numeric_limits<double>::max() ),
      gradmin( std::numeric_limits<double>::max() ),
      niter(0),
      fit_OK(false),
      functolabs_reached(false),
      steptolabs_reached(false),
      gradtolabs_reached(false),
      functolrel_reached(false),
      steptolrel_reached(false),
      nitermax_reached(false),
      trustradius_too_small(false),
      small_step(false),
      singular_matrix(false),
      bad_point(false),
      bad_value(false)
    { }

  } ;


  // ********************************************************************
  // ********************************************************************
  // ********************************************************************
  // FUNCTIONS
  // ********************************************************************
  // ********************************************************************
  // ********************************************************************


  bool check_conv(const Vector<double> x,
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
		  );

  

}




#endif

