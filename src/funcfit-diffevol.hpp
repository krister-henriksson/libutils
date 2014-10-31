


#ifndef FUNCFIT_DIFFEVOL_HPP
#define FUNCFIT_DIFFEVOL_HPP


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

#include "param.hpp"

#include "mtwister.hpp"

#include "funcfit-basics.hpp"
#include "funcfit-exceptions.hpp"




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
    [1]: Storn, Price, Journal of Global Optimization 11 (1997) 341-359
    [2]: Ali, Törn, Computers & Operations Research 31 (2004) 1703-1725
    [3]: Brest et al, IEEE Transactions on Evolutionary Computation 10 (2006) 646-657
    [4]: Qin et al, IEEE Transactions on Evolutionary Computation 13 (2009) 398-417
    [5]: Fan and Lampinen, Journal of Global Optimization 27 (2003) 105–129
   */


  class Cond_DiffEvol {
  public:
    bool ver_F_rescale;
    bool ver_CR_rescale;
    int mutation_ver;
    int NP_fac;
    double F;
    double CR;
    double Fmin;
    double Flow;
    double Fhigh;
    double tau1;
    double tau2;
    double Mt;

    Cond_DiffEvol(){
      ver_F_rescale  = 0; // 0: no scaling, 1: [2], 2: [3]
      ver_CR_rescale = 0; // 0-1: no scaling,       2: [3]
      mutation_ver = 1;
      NP_fac = 10;
      F    = 0.8;
      CR   = 0.8;
      Fmin  = 0.4;  // used if ver_F_rescale = 1
      Flow  = 0.1; // used if ver_F_rescale = 2
      Fhigh = 0.9; // used if ver_F_rescale = 2
      tau1 = 0.1; // used if ver_CR_rescale = 2
      tau2 = 0.1; // used if ver_CR_rescale = 2
      Mt   = 0.05; // TDE [5]
    }

  } ;
  
  // DE/rand/1
  // v[i] = x[i1] + F * (x[i2] - x[i3]);
  // DE/best/1
  // v[i] = x[ixmin] + F * (x[i1] - x[i2]);
  // DE/rand/2
  // v[i] = x[i1] + F * (x[i2] + x[i3] - x[i4] - x[i5]);
  // ...




  
  template <typename T>
  class DiffEvol
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
    DiffEvol(T & funcd)
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
			    Cond_DiffEvol cond_diffevol = Cond_DiffEvol()
			    ){
      // Input: starting point point and constructor-set functor (containing
      // gradient method)
      // Purpose: minimize function, input point will not be changed
      double eps = numeric_limits<double>::epsilon();
      string methodstring("differential-evolution");
      int niter=0, i,j,k, i1,i2,i3,i4,i5, rj,rni;
      int ixmin, ixmax;
      int ixamin, ixamax;
      double fxmin, fxmax, absfxmin, absfxmax, oldfxmin, oldfxmax;
      double fxamin, fxamax, absfxamin, absfxamax, oldfxamin, oldfxamax;
      double td, ratio1, ratio2;

      int seed2 = seed < 0 ? time(0) : seed;
      rand_mtwister mtwister( seed2 );

      Vector<double> point = func.free_parameters(point_in);
      Vector<double>        xmin  = func.map_vector_as_free_parameters(par_min);
      Vector<double>        xmax  = func.map_vector_as_free_parameters(par_max);
      Vector<parametertype> xtype = func.map_vector_as_free_parameters(par_type);

      bool ver_F_rescale  = cond_diffevol.ver_F_rescale;
      bool ver_CR_rescale = cond_diffevol.ver_CR_rescale;
      int  mutation_ver   = cond_diffevol.mutation_ver;
      int NP_fac  = cond_diffevol.NP_fac;
      double F     = cond_diffevol.F;
      double CR    = cond_diffevol.CR;
      double Fmin  = cond_diffevol.Fmin;
      double Flow  = cond_diffevol.Flow;
      double Fhigh = cond_diffevol.Fhigh;
      double tau1  = cond_diffevol.tau1;
      double tau2  = cond_diffevol.tau2;
      double Mt    = cond_diffevol.Mt;

      int D  = point.size();
      int NP = NP_fac * D;

      status = Minimization_Status();

      bool debug         = cond_debug.debug_fit_level0;
      string prefix_report_debug = cond_debug.prefix_debug_fit_level0;

      bool report_iter    = cond_print.report_iter;
      bool report_warn    = cond_print.report_warn;
      bool report_error   = cond_print.report_error;

      Counters_niter counters_niter = Counters_niter();

      Vector<double> diff(D, 0);
      Vector< Vector<double> > x(NP, Vector<double>(D,0));
      Vector< Vector<double> > x_new(NP, Vector<double>(D,0));
      Vector< Vector<double> > v(NP, Vector<double>(D,0));
      Vector< Vector<double> > u(NP, Vector<double>(D,0));
      Vector<double> fx(NP), fx_new(NP), fu(NP);
      Vector<int> xmin_set_idx(0);

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



      
      // --------------------------------------------------------------
      // --------------------------------------------------------------
      // Create initial population vectors
      // --------------------------------------------------------------
      // --------------------------------------------------------------
      if (debug) 
	cout << prefix_report_debug
	     << methodstring << ": "
	     << "Creating initial population ... " << endl;
      for (i=0; i<NP; ++i){
	// configuration of individual
	for (j=0; j<D; ++j){

	  x[i][j] = xmin[j] + mtwister.unif() * (xmax[j] - xmin[j]);
	  if (x[i][j] < xmin[j]) x[i][j] = xmin[j];
	  if (x[i][j] > xmax[j]) x[i][j] = xmax[j];
	  
	}
      }



      // --------------------------------------------------------------
      // --------------------------------------------------------------
      // Create initial population values
      // --------------------------------------------------------------
      // --------------------------------------------------------------
      if (debug)
	cout << prefix_report_debug
	     << methodstring << ": "
	     << "Getting merit function values for initial population ... " << endl;
      //dumpfile = "DE-points-iter" + tostring(niter) + ".out";
      for (i=0; i<NP; ++i){
	//fout << format("%10ld  ") % i << " : ";
	for (j=0; j<D; ++j) fout << format(" %15.8e") % x[i][j];
	fx[i] = func(x[i]);
	//fout << format(" %15.8e") % fx[i] << endl;
	if (i==0 || (i>0 && fxmin>fx[i])){ ixmin = i; fxmin = fx[i]; }
	if (i==0 || (i>0 && fxmax<fx[i])){ ixmax = i; fxmax = fx[i]; }
      }
      //fout.close();
      xmin_set_idx.resize(0);
      for (i=0; i<NP; ++i)
	if (fx[i]<=fxmin) xmin_set_idx.push_back(i);
      fmin = fx[ixmin];



      

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
	  for (i=0; i<D; ++i) pmin[i] = x[ixmin][i];
	  func(pmin); // triggers recalculation of properties
	  // #################################################################

	  printf("%s%s: Iter %4d  Func min %15.8e max %15.8e\n", 
		 cond_print.prefix_report_iter.c_str(),
		 methodstring.c_str(), niter, fxmin, fxmax);
	  //cout << "Chi^2 components: " << func.f() << endl;
	  //cout << "DataY: " << func.DataY() << endl;
	  //cout << "ModelDataY: " << func.ModelDataY() << endl;
	  func.report_on_parameters_and_data();
	}

	// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// Check for convergence
	// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	if (niter>0){
	  if ( check_conv(x[ixmin], true,oldfxmin, true,fxmin, // use old and current function values
			  false,Vector<double>(x[ixmin].size(),0),// use gradient
			  false,Vector<double>(x[ixmin].size(),0),// use step taken
			  niter,
			  counters_niter,
			  cond_conv, cond_debug, cond_print,
			  methodstring, status) ){
	    return func.all_parameters(x[ixmin]);
	  }
	}



	oldfxmin = fxmin;
	oldfxmax = fxmax;


	// ###################################################################
	// Rescaling of F
	// ###################################################################
       	if (ver_F_rescale==1){
	  absfxmin = fxmin * ((fxmin < 0) ? -1: 1);
	  absfxmax = fxmax * ((fxmax < 0) ? -1: 1);
	  ratio1 = absfxmax / absfxmin;
	  F = Fmin;
	  if (ratio1 < 1){
	    td = 1.0 - ratio1;
	    if (td>F) F = td;
	  }
	  else {
	    td = 1.0 - absfxmin / absfxmax;
	    if (td>F) F = td;
	  }
	}
       	else if (ver_F_rescale==2){	  
	  if (mtwister.unif()<tau1) F  = Flow + mtwister.unif()*Fhigh;
	}

	// ###################################################################
	// Rescaling of CR
	// ###################################################################
       	if (ver_CR_rescale==2){
	  if (mtwister.unif()<tau2) CR = mtwister.unif();
	}







	// ###################################################################
	// Mutation
	// ###################################################################
	if (debug)
	  cout << prefix_report_debug
	       << methodstring << ": "
	       << "Getting mutations ... " << endl;
	for (i=0; i<NP; ++i){

	  i1 = i2 = i3 = i;
	  while (i1==i || i1<0 || i1>=NP)                     i1 = floor(mtwister.unif() * NP);
	  while (i2==i || i2==i1 || i2<0 || i2>=NP)           i2 = floor(mtwister.unif() * NP);
	  while (i3==i || i3==i2 || i3==i1 || i3<0 || i3>=NP) i3 = floor(mtwister.unif() * NP);
	  
	  diff = x[i2] - x[i3];

	  // default to version 1: DE/rand/1/bin
	  v[i] = x[i1] + F * diff;

	  if (mutation_ver==2){
	    // DE/best/1/bin
	    k=0; int p = xmin_set_idx.size();
	    if (p>1){
	      k = floor( mtwister.unif() * p ); while (k< 0) k++; while (k>=p) k--;
	    }
	    v[i] = x[ xmin_set_idx[k] ] + F * diff;
	  }
	  else if (mutation_ver==3){
	    // DE/rand/2/bin
	    i4 = i5 = i;
	    while (i4==i || i4==i3 || i4==i2 || i4==i1 || i4<0 || i4>=NP)
	      i4 = floor(mtwister.unif() * NP);
	    while (i5==i || i5==i4 || i5==i3 || i5==i2 || i5==i1 || i5<0 || i5>=NP)
	      i5 = floor(mtwister.unif() * NP);
	    
	    diff = diff + x[i4] - x[i5];
	    v[i] = x[i1] + F * diff;
	  }
	  else if (mutation_ver==4 && mtwister.unif()<Mt){
	    // [5]
	    double pp = abs(fx[i1]) + abs(fx[i2]) + abs(fx[i3]);
	    double p1 = abs(fx[i1])/pp;
	    double p2 = abs(fx[i2])/pp;
	    double p3 = abs(fx[i3])/pp;
	    v[i] = (x[i1] + x[i2] + x[i3])/3.0
	      + (p2-p1)*(x[i1]-x[i2])
	      + (p3-p2)*(x[i2]-x[i3])
	      + (p1-p3)*(x[i3]-x[i1]);
	  }
	  

	  // Account for parameter limits:
	  for (j=0; j<D; ++j){
	    if (v[i][j] <  xmin[j]) v[i][j] = xmin[j];
	    if (v[i][j] >  xmax[j]) v[i][j] = xmax[j];
	  }

	}




	// ###################################################################
	// Crossover
	// ###################################################################
	if (debug) 
	  cout << prefix_report_debug
	       << methodstring << ": "
	       << "Getting crossovers ... " << endl;
	for (i=0; i<NP; ++i){
	  rni = -1;
	  while (rni<0 || rni>=D) rni = floor(mtwister.unif() * D);

	  for (j=0; j<D; ++j){
	    rj = mtwister.unif();
	    u[i][j] = x[i][j];
	    if (rj <= CR || j == rni) u[i][j] = v[i][j];
	  }
	  fu[i] = func(u[i]);
	}



	// ###################################################################
	// Selection
	// ###################################################################
	if (debug)
	  cout << prefix_report_debug
	       << methodstring << ": "
	       << "Getting selections ... " << endl;
	for (i=0; i<NP; ++i){
	  if (fu[i] < fx[i]){
	    x_new[i] = u[i];
	    fx_new[i] = fu[i];
	  }
	  else {
	    x_new[i] = x[i];
	    fx_new[i] = fx[i];
	  }
	}


	//dumpfile = "DE-points-iter" + tostring(niter) + ".out";
	//fout.open(dumpfile.c_str());
	for (i=0; i<NP; ++i){
	  fx[i] = fx_new[i];
	  x[i]  = x_new[i];

	  if (i==0 || (i>0 && fxmin>fx[i])){ ixmin = i; fxmin = fx[i]; }
	  if (i==0 || (i>0 && fxmax<fx[i])){ ixmax = i; fxmax = fx[i]; }

	  //fout << format("%10ld  ") % i << " : ";
	  //for (j=0; j<D; ++j) fout << format(" %15.8e") % x[i][j];
	  //fout << format(" %15.8e") % fx[i] << endl;
	}
	//fout.close();
	xmin_set_idx.resize(0);
	for (i=0; i<NP; ++i)
	  if (fx[i]<=fxmin) xmin_set_idx.push_back(i);
	fmin = fx[ixmin];


	dumpfile = "points-DE-iter" + tostring(niter) + ".out";
	fout.open(dumpfile.c_str());
	for (i=0; i<NP; ++i){
	  fout << format("%10d  ") % i << " : ";
	  for (j=0; j<D; ++j) fout << format(" %15.8e") % x[i][j];
	  fout << format(" %15.8e") % fx[i] << endl;
	}
	fout.close();








	// Go to next iteration.
	niter++;
      } // End of loop over iterations
      // ###################################################################
      // ###################################################################
      // ###################################################################
      // ###################################################################



      fmin = fx[ixmin];
      return func.all_parameters(x[ixmin]);

    } // End of method.



  };



}





#endif



