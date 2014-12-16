


#ifndef FUNCFIT_BEECOLONY_HPP
#define FUNCFIT_BEECOLONY_HPP


#include <iostream>
#include <string>
#include <fstream>
#include <limits>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>

#include <boost/format.hpp>

#include "utils.hpp"
#include "utils-math.hpp"
#include "utils-vector.hpp"
#include "utils-string.hpp"
#include "utils-errors.hpp"

#include "param.hpp"

#include "mtwister.hpp"

#include "funcfit-basics.hpp"
#include "funcfit-errors.hpp"




///////////////////////////////////////////////////////////
// Usage:
///////////////////////////////////////////////////////////
//
//   Funcd fd;
//
//   FuncFitDE<Funcd> fgn(fd);
//
//   xmin = fgn.minimize(...);
//   fmin = fgn.status.funcmin;
//
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////


using namespace utils;
using namespace funcfit;
using std::cout;
using std::endl;
using std::ofstream;
using std::numeric_limits;
using boost::format;


namespace funcfit {


  /*
    [1]: Akay and Karboga, Information sciences 192 (2012) 120-142
    [2]: http://mf.erciyes.edu.tr/abc/
   */


  class Cond_BeeColony {
  public:
    int NP_fac;

    Cond_BeeColony(){
      NP_fac = 10;
    }

  } ;




  
  template <typename T>
  class BeeColony
  {


    //#####################################################################
    // Member variables
    //#####################################################################
  public:
    double fmin;
    T & func;
    Minimization_Status status;



    //#####################################################################
    // Constructor
    //#####################################################################
    // NB:  Use constructor to set functor. Point and direction are set by
    // NB:  minimize() method.
    BeeColony(T & funcd)
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
			    int seed,
			    Cond_Conv  cond_conv  = Cond_Conv(),
			    Cond_Debug cond_debug = Cond_Debug(),
			    Cond_Print cond_print = Cond_Print(),
			    Cond_BeeColony cond_beecolony = Cond_BeeColony()
			    ){
      // Input: starting point point and constructor-set functor (containing
      // gradient method)
      // Purpose: minimize function, input point will not be changed
      double eps = numeric_limits<double>::epsilon();
      string methodstring("bee-colony");
      int niter=0, i,j,k, i1,i2,i3, ixglobmin, ixglobmax;
      double vmax, w;

      int seed2 = seed < 0 ? time(0) : seed;
      rand_mtwister mtwister( seed2 );

      Vector<double> point = func.free_parameters(point_in);
      Vector<double>        xmin  = func.map_vector_as_free_parameters(par_min);
      Vector<double>        xmax  = func.map_vector_as_free_parameters(par_max);
      Vector<parametertype> xtype = func.map_vector_as_free_parameters(par_type);

      int NP_fac  = cond_beecolony.NP_fac;

      int D  = point.size();
      int NP = NP_fac * D;

      status = Minimization_Status();

      bool debug         = cond_debug.debug_fit_level0;
      string prefix_report_debug = cond_debug.prefix_debug_fit_level0;

      bool report_iter    = cond_print.report_iter;
      bool report_warn    = cond_print.report_warn;
      bool report_error   = cond_print.report_error;

      Counters_niter counters_niter = Counters_niter();

      double fxglobmin, fxglobmax;
      double fxglobmin_old, fxglobmax_old;

      Vector< Vector<double> > x(NP, Vector<double>(D,0));
      Vector<double> fx(NP);
      double limit = 0.5*NP*D;
      double fv, fi;
      Vector<double> v(D, 0);
      Vector<double> trials(NP, 0);
      Vector<double> fitness(NP);
      Vector<double> prob(NP);
      Vector<double> xglobmin(D);
      Vector<double> xglobmax(D);

      double td1,td2,td3, rtol;

      string dumpfile;
      ofstream fout;

      double dmax = sqrt( numeric_limits<double>::max() );
      double dmin = -dmax;
      for (i=0; i<D; ++i){
	if (xtype[i]==PARAM_FREE){ // no lower or upper limit
	  xmin[i] = dmin;
	  xmax[i] = dmax;
	}
      }


      niter = 0;
      




      // --------------------------------------------------------------
      // --------------------------------------------------------------
      // Create initial population vectors and values
      // --------------------------------------------------------------
      // --------------------------------------------------------------
      if (debug) 
	cout << prefix_report_debug
	     << methodstring << ": "
	     << "Creating initial population vectors and values ... " << endl;
      for (i=0; i<NP; ++i){
	cout << methodstring << ": "
	     << "Agent " << i+1 << " of " << NP << endl;

	// **************************************************************
	while (true){
	  try {
	    for (j=0; j<D; ++j){
	      x[i][j] = xmin[j] + mtwister.unif() * (xmax[j] - xmin[j]);
	      if (x[i][j] < xmin[j]) x[i][j] = xmin[j];
	      if (x[i][j] > xmax[j]) x[i][j] = xmax[j];
	    }
	    fx[i] = func(x[i]);
	    fitness[i] = 1.0/(1.0 + fx[i]);
	  }
	  catch (funcfit::bad_point & err_bad_point){
	    func.reset();
	    cout << "Warning: Bad point " << x[i] << ". Recreating point ..." << endl;
	    cout << "Recommendation: Restate parameter limits and start over." << endl;
	    continue;
	  }
	  break;
	}
	// **************************************************************

	if (i==0 || (i>0 && fxglobmin>fx[i])){ ixglobmin = i; fxglobmin = fx[i]; }
	if (i==0 || (i>0 && fxglobmax<fx[i])){ ixglobmax = i; fxglobmax = fx[i]; }
      }
      xglobmin = x[ixglobmin];
      xglobmax = x[ixglobmax];
      fmin = fxglobmin;



      dumpfile = "points-BC-iter" + tostring(niter) + ".out";
      fout.open(dumpfile.c_str());
      for (i=0; i<NP; ++i){
	fout << format("%10d  ") % i << " : ";
	for (j=0; j<D; ++j) fout << format(" %15.8e") % x[i][j];
	fout << format(" %15.8e") % fx[i] << endl;
      }
      fout.close();







      if (debug)
	cout << prefix_report_debug
	     << methodstring << ": "
	     << "Starting iterations ... " << endl;


      // ###################################################################
      // ###################################################################
      // ###################################################################
      // ###################################################################
      while (true){ // Start of loop over iterations

	if (debug) cout << "Iteration " << niter << ":" << endl;
	if (debug)
	  for (i=0; i<x.size(); ++i)
	    cout << "x[" << i << "]: " << x[i] << endl;



	// ***************************************************************
	// ***************************************************************
	// Report
	// ***************************************************************
	// ***************************************************************
	if (report_iter){
	  // #################################################################
	  // Last point evaluated by the function is one of the vertices/individuals
	  // (ex. it could be the worst point or the best point). Now need
	  // to report about the best point:
	  Vector<double> pmin(D, 0.0);
	  for (i=0; i<D; ++i) pmin[i] = x[ixglobmin][i];
	  func(pmin); // triggers recalculation of properties
	  // #################################################################

	  printf("%s%s: Iter %4d  Func min %15.8e max %15.8e\n", 
		 cond_print.prefix_report_iter.c_str(),
		 methodstring.c_str(), niter, fxglobmin, fxglobmax);
	  //cout << "Chi^2 components: " << func.f() << endl;
	  //cout << "DataY: " << func.DataY() << endl;
	  //cout << "ModelDataY: " << func.ModelDataY() << endl;
	  func.report_on_parameters_and_data();
	}

	// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// Check for convergence
	// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	if ( check_conv(x[ixglobmin], true,fxglobmin_old, true,fxglobmin, // use old and current function values
			false,Vector<double>(x[ixglobmin].size(),0),// use gradient
			false,Vector<double>(x[ixglobmin].size(),0),// use step taken
			niter,
			counters_niter,
			cond_conv, cond_debug, cond_print,
			methodstring, status) ){
	  return func.all_parameters(x[ixglobmin]);
	}



	fxglobmin_old = fxglobmin;
	fxglobmax_old = fxglobmax;




	// ###################################################################
	// ###################################################################
	// No convergence yet.
	// ###################################################################
	// ###################################################################




	// ###################################################################
	// Employed phase
	// ###################################################################
	if (debug)
	  cout << prefix_report_debug
	       << methodstring << ": "
	       << "Employed phase ... " << endl;

	double probsum=0.0;
	for (i=0; i<NP; ++i){
	  i1 = i2 = -1;
	  while (i1<0 || i1>=D)           i1 = floor(mtwister.unif() * D);
	  while (i2==i || i2<0 || i2>=NP) i2 = floor(mtwister.unif() * NP);
	  v = x[i];
	  v[i1] = x[i][i1] + (2*mtwister.unif()-1) * (x[i][i1] - x[i2][i1]);
	  if (v[i1] < xmin[i1]) v[i1] = xmin[i1];
	  if (v[i1] > xmax[i1]) v[i1] = xmax[i1];
	    
	  fv = func(v);
	  fi = 1.0/(1+fv);
	  if (fi < fitness[i]){
	    trials[i]++;
	  }
	  else {
	    x[i] = v;
	    fx[i] = fv;
	    trials[i]=0;
	    fitness[i] = fi;
	  }
	  probsum += fitness[i];
	}
	// Normalization:
	for (i=0; i<NP; ++i){
	  prob[i] = fitness[i] / probsum;
	}


	// ###################################################################
	// Onlooker phase
	// ###################################################################
	if (debug)
	  cout << prefix_report_debug
	       << methodstring << ": "
	       << "Onlooker phase ... " << endl;

	for (i=0; i<NP; ++i){

	  double q = mtwister.unif();
	  k=-1;
	  for (j=1; j<NP; ++j){	  
	    if (q>=prob[j-1] && q<prob[j]){
	      k = j; break;
	    }
	  }
	  if (k==-1) k=NP-1;

	  i1 = i2 = -1;
	  while (i1<0 || i1>=D)           i1 = floor(mtwister.unif() * D);
	  while (i2==i || i2<0 || i2>=NP) i2 = floor(mtwister.unif() * NP);
	  v = x[k];
	  v[i1] = x[k][i1] + (2*mtwister.unif()-1) * (x[k][i1] - x[i2][i1]);
	  if (v[i1] < xmin[i1]) v[i1] = xmin[i1];
	  if (v[i1] > xmax[i1]) v[i1] = xmax[i1];

	  fv = func(v);
	  fi = 1.0/(1+fv);

	  if (fi < fitness[k]){
	    trials[k]++;
	  }
	  else {
	    x[k] = v;
	    fx[k] = fv;
	    trials[k]=0;
	    fitness[k] = 1.0/(1.0 + fx[k]);
	  }

	}






	// ###################################################################
	// Scout phase
	// ###################################################################
	if (debug)
	  cout << prefix_report_debug
	       << methodstring << ": "
	       << "Scout phase ... " << endl;
	
	for (i=0; i<NP; ++i){
	  cout << methodstring << ": "
	       << "Agent " << i+1 << " of " << NP << endl;

	  if (trials[i]>=limit){
	    // Replace

	    // **************************************************************
	    while (true){
	      try {
		for (j=0; j<D; ++j){
		  x[i][j] = xmin[j] + mtwister.unif() * (xmax[j] - xmin[j]);
		  if (x[i][j] < xmin[j]) x[i][j] = xmin[j];
		  if (x[i][j] > xmax[j]) x[i][j] = xmax[j];
		}
		fx[i] = func(x[i]);
	      }
	      catch (funcfit::bad_point & err_bad_point){
		func.reset();
		cout << "Warning: Bad point " << x[i] << ". Recreating point ..." << endl;
		cout << "Recommendation: Restate parameter limits and start over." << endl;
		continue;
	      }
	      break;
	    }
	    fitness[i] = 1.0/(1.0 + fx[i]);
	    trials[i]=0;
	  }
	}





	dumpfile = "points-BC-iter" + tostring(niter) + ".out";
	fout.open(dumpfile.c_str());
	for (i=0; i<NP; ++i){
	  fout << format("%10d  ") % i << " : ";
	  for (j=0; j<D; ++j) fout << format(" %15.8e") % x[i][j];
	  fout << format(" %15.8e") % fx[i] << endl;
	}
	fout.close();


	// ###################################################################
	// Determine best
	// ###################################################################
	if (debug)
	  cout << prefix_report_debug
	       << methodstring << ": "
	       << "Determining best ... " << endl;
	
	for (i=0; i<NP; ++i){
	  if (i==0 || (i>0 && fxglobmin>fx[i])){ ixglobmin = i; fxglobmin = fx[i]; }
	  if (i==0 || (i>0 && fxglobmax<fx[i])){ ixglobmax = i; fxglobmax = fx[i]; }
	}
	xglobmin = x[ixglobmin];
	xglobmax = x[ixglobmax];
	fmin = fxglobmin;


	niter++;


      } // End of loop over iterations
      // ###################################################################
      // ###################################################################
      // ###################################################################
      // ###################################################################


      fmin = fxglobmin;
      return func.all_parameters(x[ixglobmin]);


    } // End of method.



  };



}





#endif




#if 0

      if (ndiffvec<=0) ndiffvec = 1;
      if ((2*ndiffvec+2)>NP) NP = 2*ndiffvec + 2;

      Vector<int> idiffvec(2*ndiffvec+1, 0);


	  // Code to handle optional number of positive terms in the difference vector (NPTDV):
	  ndt = 2*ndiffvec + 1;
	  if (mutate_best) ndt = 2*ndiffvec;

	  for (ia=0; ia<ndt; ++ia){
	    not_ready=true;
	    while (not_ready){
	      not_ready=false;
	      // get value
	      idiffvec[ia] = floor(mtwister.unif() * NP);
	      if (idiffvec[ia]==i){// || idiffvec[ia]<0 || idiffvec[ia]>=NP-1){
		not_ready=true; continue;
	      }
	      // compare with previous values
	      for (ib=0; ib<ia; ++ib){
		if (idiffvec[ia]==idiffvec[ib]){
		  not_ready=true; break;
		}
	      }
	      // if value was found among the previous ones we are not ready
	      // and will continue inside loop
	    }
	  }

	  // Build difference vector (DV):
	  // cout << "Indices for mutation vector: " << i;
	  for (ia=0; ia<ndiffvec; ++ia){
	    diff = diff + x[idiffvec[ia]];
	    //cout << " " << idiffvec[ia];
	  }
	  for (ia=ndiffvec; ia<2*ndiffvec; ++ia){
	    diff = diff - x[idiffvec[ia]];
	    //cout << " " << idiffvec[ia];
	  }

#endif


