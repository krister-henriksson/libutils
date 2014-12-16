

#ifndef FUNCFIT_SIMPLEX_HPP
#define FUNCFIT_SIMPLEX_HPP


#include <iostream>
#include <string>
#include <limits>

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "utils.hpp"
#include "utils-math.hpp"
#include "utils-vector.hpp"
#include "utils-errors.hpp"

#include "param.hpp"

#include "mtwister.hpp"

#include "funcfit-basics.hpp"
#include "funcfit-errors.hpp"



using std::cout;
using std::endl;
using std::numeric_limits;
using namespace funcfit;
using utils::Vector;
using utils::Matrix;
using utils::abs;
using utils::swap;


namespace funcfit {

  
  template <typename T>
  class SimplexFit
  {


    //#####################################################################
    // Member variables
    //#####################################################################
  public:
    int nfunc;
    int mpts;
    int ndim;
    double fmin;
    Vector<double> point;
    Vector<double> dels;
    Vector<double> y; // Function values at the vertices of the simplex.
    Matrix<double> p; // Current simplex.
    T & func;
    Minimization_Status status;



    //#####################################################################
    // Constructor
    //#####################################################################
    // NB:  Use constructor to set functor. Point and direction are set by
    // NB:  minimize() method.
    SimplexFit(T & funcd)
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
			    int seed,
			    Cond_Conv  cond_conv  = Cond_Conv(),
			    Cond_Debug cond_debug = Cond_Debug(),
			    Cond_Print cond_print = Cond_Print()
			    ){
      // Input: starting point point and constructor-set functor (containing
      // gradient method)
      // Purpose: minimize function, input point will not be changed
      double eps = numeric_limits<double>::epsilon(), fp, fp_old, ytry, ysave;
      string methodstring("simplex");
      int niter=0, i,j;

      int seed2 = seed < 0 ? time(0) : seed;
      rand_mtwister mtwister( seed2 );

      status = Minimization_Status();


      double FUNCTOLABS   = cond_conv.functolabs;
      double FUNCTOLREL   = cond_conv.functolrel;
      //double GRADTOLABS   = cond_conv.gradtolabs;
      //double STEPTOLABS   = cond_conv.steptolabs;
      //double STEPTOLREL   = cond_conv.steptolrel;

      int NITERMIN       = cond_conv.nitermin;
      int NITERMAX       = cond_conv.nitermax;

      bool debug         = cond_debug.debug_fit_level0;
      string prefix_report_debug = cond_debug.prefix_debug_fit_level0;

      bool report_iter    = cond_print.report_iter;
      bool report_warn    = cond_print.report_warn;
      bool report_error   = cond_print.report_error;

      bool report_conv    = cond_conv.report_conv;

      //Counters_niter counters_niter = Counters_niter();

      point = func.free_parameters(point_in);
      Vector<double>        xmin  = func.map_vector_as_free_parameters(par_min);
      Vector<double>        xmax  = func.map_vector_as_free_parameters(par_max);
      Vector<parametertype> xtype = func.map_vector_as_free_parameters(par_type);


      ndim = point.size();


      double dmax = sqrt( numeric_limits<double>::max() );
      double dmin = -dmax;
      for (i=0; i<ndim; ++i){
	if (xtype[i]==PARAM_FREE){ // no lower or upper limit
	  xmin[i] = dmin;
	  xmax[i] = dmax;
	}
      }

      int ihi, ilo, inhi;



      p.resize(ndim+1, ndim);
      mpts = ndim+1;
      Vector<double> psum(ndim, 0), pmin(ndim, 0), x(ndim, 0);
      y.resize(mpts);



      if (debug) 
	cout << prefix_report_debug
	     << methodstring << ": "
	     << "Getting vertices and merit function values ... " << endl;

      for (i=0; i<ndim+1; ++i){ // i is for vertex
	cout << methodstring << ": "
	     << "Vertex " << i+1 << " of " << mpts << endl;

	// **************************************************************
	while (true){
	  try {

	    for (j=0; j<ndim; ++j)  // j is for parameter
	      p.elem(i,j) = point[j];

	    if (i > 0){
	      j = i-1;
	      // p.elem(i, i-1) += xch[i-1];

	      p.elem(i,j) = xmin[j] + mtwister.unif() * (xmax[j]-xmin[j]);
	      if (p.elem(i,j) < xmin[j]) p.elem(i,j) = xmin[j];
	      if (p.elem(i,j) > xmax[j]) p.elem(i,j) = xmax[j];
	    }

	    for (j=0; j<ndim; ++j)
	      x[j] = p.elem(i,j);

	    y[i] = func(x);
	  }
	  catch (funcfit::bad_point & err_bad_point){
	    func.reset();
	    cout << "Warning: Bad point " << x << ". Recreating point ..." << endl;
	    cout << "Recommendation: Restate parameter limits and start over." << endl;
	    continue;
	  }
	  break;
	}

      }
      fp = y[0];


      nfunc = 0;
      //get_psum(p, psum);
      for (j=0; j<ndim; ++j){
	double sum=0.0;
	for (i=0; i<mpts; ++i)
	  sum += p.elem(i,j);
	psum[j]=sum;
      }





      
      // ###################################################################
      // ###################################################################
      // ###################################################################
      // ###################################################################
      bool swapped;

      while (true){


	// ********************************************************
	// Determine which point is highest, next-highest and lowest.
	ilo = 0;
	ihi = y[0]>y[1] ? (inhi=1, 0) : (inhi=0, 1);
	for (i=0; i<mpts; ++i){
	  if (y[i] <= y[ilo]) ilo = i;
	  if (y[i] >  y[ihi]){
	    inhi=ihi;
	    ihi=i;
	  }
	  else if (y[i] > y[inhi] && i!=ihi)
	    inhi=i;
	}
	fp = y[ilo];
	// ********************************************************
	swapped=false;


	if (report_iter){
	  // #################################################################
	  // Last point evaluated by the function is one of the vertices
	  // (ex. it could be the worst point or the best point). Now need
	  // to report about the best point:
	  for (i=0; i<ndim; ++i)
	    pmin[i] = p.elem(ilo,i);
	  func(pmin); // triggers recalculation of properties
	  // #################################################################

	  if (niter==0){
	    printf("%s%s: Iter %4d  Func %15.8e\n",
		   cond_print.prefix_report_iter.c_str(),
		   methodstring.c_str(), niter, fp);
	  }
	  else {
	    printf("%s%s: Iter %4d  Func %15.8e  Change %15.8e\n", 
		   cond_print.prefix_report_iter.c_str(),
		   methodstring.c_str(), niter, fp, fp-fp_old);
	  }
	  func.report_on_parameters_and_data();
	}

	status.funcmin = fp;
	status.niter   = niter;





	fp_old = fp;





	if (niter>=NITERMAX) break;



	// ********************************************************
	// Calculate the fractional range from highest to lowest.
	// Return if good enough.
	double td1 = y[ihi]-y[ilo]; if (td1<0) td1 *= -1;
	double td2 = y[ihi]; if (td2<0) td2 *= -1;
	double td3 = y[ilo]; if (td3<0) td3 *= -1;
	double rtol = 2.0 * td1 / (td2 + td3 + eps);
	if (niter > NITERMIN && FUNCTOLREL > 0.0 && rtol < FUNCTOLREL){
	  if (report_conv)
	    cout << cond_conv.prefix_report_conv
		      << methodstring << ": "
		      << "Simplex size has reached convergence within limit " << FUNCTOLREL
		      << endl;
	  status.functolrel_reached = true;
	  if (!swapped){
	    swap(y[0], y[ilo]);
	    for (i=0; i<ndim; ++i){
	      swap(p.elem(0,i), p.elem(ilo, i));
	      pmin[i] = p.elem(0,i);
	    }
	    swapped = true;
	  }
	  fmin = y[0];
	  fp = fmin;
	  status.fit_OK = true;
	  return func.all_parameters(pmin);
	}
	// ********************************************************

	if (niter > NITERMIN && FUNCTOLABS > 0.0 && abs(fp) < FUNCTOLABS){
	  if (report_conv)
	    cout << cond_conv.prefix_report_conv
		      << methodstring << ": "
		      << "Absolute function value " << abs(fp) << " "
		      << "has reached convergence within limit " << FUNCTOLABS
		      << endl;
	  status.functolabs_reached = true;
	  if (!swapped){
	    swap(y[0], y[ilo]);
	    for (i=0; i<ndim; ++i){
	      swap(p.elem(0,i), p.elem(ilo, i));
	      pmin[i] = p.elem(0,i);
	    }
	    swapped = true;
	  }
	  fmin = y[0];
	  fp = fmin;
	  status.fit_OK = true;
	  return func.all_parameters(pmin);
	}






	nfunc += 2;
	
	// ********************************************************
	// Start new iteration.
	// ********************************************************

	// ########################################################
	// Reflect the simplex from he high point.
	ytry = amotry(p, y, psum, ihi, -1.0, xmin, xmax, cond_debug, cond_print);

	if (ytry <= y[ilo]){
	  // Result better than best point. Try some more.
	  // ########################################################
	  if (debug) 
	    cout << prefix_report_debug
		      << methodstring << ": "
		      << "Getting trial merit function value ... " << endl;
	  ytry = amotry(p, y, psum, ihi, 2.0, xmin, xmax, cond_debug, cond_print);
	}
	else if (ytry >= y[inhi]){
	  // Reflected point worst than next highest. Look for intermediate
	  // lower point, i.e. do a contraction.
	  ysave = y[ihi];
	  // ########################################################
	  if (debug)
	    cout << prefix_report_debug
		      << methodstring << ": "
		      << "Getting trial merit function value ... " << endl;
	  ytry  = amotry(p, y, psum, ihi, 0.5, xmin, xmax, cond_debug, cond_print);
	  if (ytry >= ysave){
	    // Can't get rid of high point.
	    // Try contracting around lowest point.
	    for (i=0; i<mpts; ++i){
	      if (i != ilo){
		for (j=0; j<ndim; ++j){
		  psum[j] = 0.5 * (p.elem(i,j) + p.elem(ilo, j));

		  if (psum[j] < xmin[j]) psum[j] = xmin[j];
		  if (psum[j] > xmax[j]) psum[j] = xmax[j];

		  p.elem(i,j) = psum[j];
		}
		y[i] = func(psum);
	      }
	    }
	    nfunc += ndim;

	    //get_psum(p, psum);
	    for (j=0; j<ndim; ++j){
	      double sum=0.0;
	      for (i=0; i<mpts; ++i)
		sum += p.elem(i,j);
	      psum[j]=sum;
	    }

	  }
	}
	else
	  --nfunc;


	


	// Go to next iteration.
	niter++;
      } // End of while loop.
      // ###################################################################
      // ###################################################################
      // ###################################################################
      // ###################################################################



      if (report_warn)
	cout << cond_print.prefix_report_error
	     << methodstring << ": "
	     << "Too many iterations (max = " << NITERMAX << "). Exiting." << endl;
      status.nitermax_reached = true;
      if (!swapped){
	swap(y[0], y[ilo]);
	for (i=0; i<ndim; ++i){
	  swap( p.elem(0,i), p.elem(ilo, i) );
	  pmin[i] = p.elem(0,i);
	}
	swapped = true;
      }
      fmin = y[0];
      fp = fmin;
      return func.all_parameters(pmin);

    } // End of method.




    /*
    inline void get_psum(Matrix<double> & p,
			 Vector<double> & psum){
      for (int j=0; j<ndim; ++j){
	double sum=0.0;
	for (int i=0; i<mpts; ++i)
	  sum += p.elem(i,j);
	psum[j]=sum;
      }
    }
    */


    double amotry(Matrix<double> & p,
		  Vector<double> & y,
		  Vector<double> & psum,
		  const int ihi,
		  const double fac,
		  Vector<double> & xmin,
		  Vector<double> & xmax,
		  Cond_Debug cond_debug,
		  Cond_Print cond_print
		  ){
      Vector<double> ptry(ndim, 0);
      double fac1 = (1.0-fac)/ndim;
      double fac2 = fac1 - fac;
      double ytry;
      funcfit::stationary_point ess;
      int j;

      string prefix_report_debug = cond_debug.prefix_debug_fit_level1;


      for (j=0; j<ndim; ++j){
	ptry[j] = psum[j] * fac1 - p.elem(ihi, j) * fac2;

	if (ptry[j] < xmin[j]) ptry[j] = xmin[j];
	if (ptry[j] > xmax[j]) ptry[j] = xmax[j];
      }
      if (cond_debug.debug_fit_level1)
	cout << prefix_report_debug
	     << "simplexfit: "
	     << "Trial point: " << ptry << endl;

      
      // ********************************************************
      // ********************************************************
      while (true){
	try {
	  if (cond_debug.debug_fit_level1)
	    cout << prefix_report_debug
		      << "simplexfit: "
		      << "Getting trial merit function value ... " << endl;
	  ytry = func(ptry);
	}
	catch (funcfit::bad_point & e1){
	  func.reset();

	  // Went to far. Try again with smaller step.
	  if (cond_print.report_warn)
	    cout << cond_print.prefix_report_warn
		      << "simplexfit: "
		      << "Too large trial step. Halving and retrying." << endl;
	  for (j=0; j<ndim; ++j){
	    ptry[j] *= 0.90;

	    if (ptry[j] < xmin[j]) ptry[j] = xmin[j];
	    if (ptry[j] > xmax[j]) ptry[j] = xmax[j];
	  }
	  
	  if (fp_is_small(ptry.magn())){
	    status.small_step = true;
	    if (cond_print.report_error)
	      cout << cond_print.prefix_report_error
			<< "simplexfit: "
			<< "Too small step. Quitting." << endl;
	    status.small_step = true;
	    throw ess;
	  }
	  continue;
	}
	break;
      }
      // ********************************************************
      // ********************************************************


      if (ytry < y[ihi]){
	y[ihi] = ytry;
	for (j=0; j<ndim; ++j){
	  psum[j] += ptry[j] - p.elem(ihi, j);
	  p.elem(ihi, j) = ptry[j];
	}
      }
      return ytry;
    }

  };



}





#endif

