

#include <iostream>
#include <boost/format.hpp>

#include "utils.hpp"
#include "utils-math.hpp"
#include "utils-vector.hpp"
#include "utils-errors.hpp"

#include "param.hpp"

#include "funcfit-errors.hpp"


using namespace std;
using boost::format;


/* ##############################################################################

   --------------------------------------------------------
              |                   Xlimuse
   --------------------------------------------------------
              | true                  |    false
   --------------------------------------------------------
              | Xmin < Xmax => free   |
   X:         |                       |   always free
              | Xmin = Xmax => fixed* |
   --------------------------------------------------------
	      | X is constrained      | X is unconstrained
   --------------------------------------------------------
   *) Special case: If Xmin = 0 = Xmax => X is free (unconstrained).

   Note: Barrier/penalty functions are only used with constrained parameters.

   ############################################################################## */



void set_param_type(double xmin, double xmax, parametertype & xtype){
  xtype = PARAM_FIXED;

  if (xmin < xmax) xtype = PARAM_FREE_WITH_LIMITS;
  
  if (utils::fp_is_small(xmin) && utils::fp_is_small(xmax))
    xtype = PARAM_FREE;
}




Param::Param(){ }


Param::Param(const utils::Vector<double> & xi){
  mX     = xi;
  mXtype = utils::Vector<parametertype>(xi.size(), PARAM_FREE);
  mXmin  = utils::Vector<double>(xi.size(), 0);
  mXmax  = utils::Vector<double>(xi.size(), 0);
}

// Copy constructor:
Param::Param(const Param & sv){
  mX     = sv.mX;
  mXtype = sv.mXtype;
  mXmin  = sv.mXmin;
  mXmax  = sv.mXmax;
}

// Operators:
Param & Param::operator=(const Param & sv){
  if (this==&sv) return *this;
  mX     = sv.mX;
  mXtype = sv.mXtype;
  mXmin  = sv.mXmin;
  mXmax  = sv.mXmax;
  return *this;
}


  // **************************************************************************
  // **************************************************************************
  // **
  // **  Methods that may be redefined in derived classes:
  // **
  // **************************************************************************
  // **************************************************************************
  // obj.X(0) returns mX[0] for reading as well as writing to, so this construction
  // enables obj.X(0) = val


  /* IO parameters are all input/output parameters visible to user. Some of these
     parameters are constrained, others are free. Only the free parameters must
     be used in any optimization routine. 
  */



void Param::init(const int N){
  mX     = utils::Vector<double>(N, 0);
  mXtype = utils::Vector<parametertype>(N, PARAM_FREE);
  mXmin  = utils::Vector<double>(N, 0);
  mXmax  = utils::Vector<double>(N, 0);
}

void Param::init(const utils::Vector<double> & xi){
  mX     = xi;
  mXtype = utils::Vector<parametertype>(xi.size(), PARAM_FREE);
  mXmin  = utils::Vector<double>(xi.size(), 0);
  mXmax  = utils::Vector<double>(xi.size(), 0);
}

void Param::init(const utils::Vector<double> & xi,
		 const utils::Vector<double> & ximin,
		 const utils::Vector<double> & ximax){
  mX    = xi;
  mXmin = ximin;
  mXmax = ximax;

  // Assume all paramteres are free:
  mXtype = utils::Vector<parametertype>(xi.size(), PARAM_FREE);
  
  // Check:
  for (int i=0; i<xi.size(); ++i){
    if (mXmin[i] < mXmax[i])
      mXtype[i] = PARAM_FREE_WITH_LIMITS;
    else
      mXtype[i] = PARAM_FIXED;

    if (utils::fp_is_small(mXmin[i]) && utils::fp_is_small(mXmax[i]))
      mXtype[i] = PARAM_FREE;
  }
}


void Param::init(const utils::Vector<double> & xi,
		 const utils::Vector<parametertype> & xitype,
		 const utils::Vector<double> & ximin,
		 const utils::Vector<double> & ximax){
  mX     = xi;
  mXtype = xitype;
  mXmin  = ximin;
  mXmax  = ximax;
}


// ****************************************************************************
// ****************************************************************************
// All parameters have lower and upper limits.
// Lack of any suffix such as 'all' or 'free' for a vector means that it
// contains all parameters.

// Set/get any parameter
double & Param::X(int i){ return mX[i]; }
// Set/get all parameters
utils::Vector<double> & Param::X(void) { return mX; }
 
// Set/get min parameter limit
double & Param::Xmin(int i){ return mXmin[i]; }
// Set/get all min parameter limits
utils::Vector<double> & Param::Xmin(void) { return mXmin; }

// Set/get max parameter limit
double & Param::Xmax(int i){ return mXmax[i]; }
// Set/get all max parameter limits
utils::Vector<double> & Param::Xmax(void) { return mXmax; }

// Set/get parameter type
parametertype & Param::Xtype(int i){ return mXtype[i]; }
// Set/get all parameter types
utils::Vector<parametertype> & Param::Xtype(void){ return mXtype; }




// ****************************************************************************
// ****************************************************************************


bool Param::X_is_free(const int i){
  if (mXtype[i] == PARAM_FREE ||
      mXtype[i] == PARAM_FREE_WITH_LIMITS) return true;
  else return false;
}
bool Param::X_is_fixed(const int i){
  if (mXtype[i] == PARAM_FIXED) return true;
  else return false;
}
bool Param::X_is_unconstrained(const int i){
  if (mXtype[i] == PARAM_FREE) return true;
  else return false;
}
  
int Param::NXfree(void){
  int n=0;
  for (int i=0; i<mX.size(); ++i){
    if (mXtype[i] == PARAM_FREE || mXtype[i] == PARAM_FREE_WITH_LIMITS)
      ++n;
  }
  return n;
}



// ##################################################################################
// ##################################################################################
//
// Update parametrization
// Given any vector (of correct type of course), update all the parameters
// (the full parameter vector)
//
// ##################################################################################
// ##################################################################################
utils::Vector<double> Param::Xupdate(const utils::Vector<double> & xi){
  int NX = mX.size(), NXin = xi.size();

  if (NXin==0) return mX;

  if (NX==NXin){
    // Input vector has full size => contains all parameters. Make a plain copy.
    mX = xi;
  }
  else {
    // Input vector has less than full size => contains free parameters.
    int n=0;
    for (int i=0; i<mX.size(); ++i){
      if (n >= NXin) break;

      if (mXtype[i] == PARAM_FREE || mXtype[i] == PARAM_FREE_WITH_LIMITS){
	mX[i] = xi[n]; // free parameter
	++n;
      }
    }
  }

  funcfit::bad_point ebp;
  if (! Xfree_is_good()){

    for (int i=0; i<mX.size(); ++i){
      if ( mXtype[i] == PARAM_FREE_WITH_LIMITS &&
	   (mX[i]<mXmin[i] || mX[i]>mXmax[i]) ){
	cout << "Parameter x[" << i << "/" << NXfree() << "] is "
	     << format("%20.10e") % mX[i] << ". "
	     << "Limits: "
	     << format("%20.10e") % mXmin[i] << ", "
	     << format("%20.10e") % mXmax[i] << endl;
      }
    }
    
    throw ebp;
  }

  return mX;
}
// ##################################################################################
// ##################################################################################


// Return a vector with all parameters:
utils::Vector<double> Param::X(const utils::Vector<double> & xi){
  Xupdate(xi);
  return mX;
}



// Return a vector with the free parameters:
utils::Vector<double> Param::Xfree(void){
  utils::Vector<double> xout(mX.size(), 0);
  int n=0;
  for (int i=0; i<mX.size(); ++i){
    if (mXtype[i] == PARAM_FREE || mXtype[i] == PARAM_FREE_WITH_LIMITS){
      // free parameter
      xout[n] = mX[i];
      ++n;
    }
  }
  xout.resize(n);
  return xout;
}

// Return a vector with the free parameters:
utils::Vector<double> Param::Xfree(const utils::Vector<double> & xi){
  Xupdate(xi);
  return Xfree();
}


// ##################################################################################

// Inquire if Xfree[] is inside the limits.
bool Param::Xfree_is_good(void){
  for (int i=0; i<mX.size(); ++i){
    if ( mXtype[i] == PARAM_FREE_WITH_LIMITS &&
	 (mX[i]<mXmin[i] || mX[i]>mXmax[i]) ) return false;
  }
  return true;
}

// Inquire if Xfree[] is inside the limits.
bool Param::Xfree_is_good(const utils::Vector<double> & xi){
  Xupdate(xi);
  return Xfree_is_good();
}



// ##################################################################################


// There may be vectors A that have the same size as X (the full vector of all
// parameters). At some places the elements of A corresponding to the free
// parameters may be needed. To simplify operations, a method that converts A
// to a (smaller) vector where elements A[i] corresponds to free parameter Xfree[i]
// is needed.

utils::Vector<double> Param::map_any_as_Xfree(const utils::Vector<double> & any){
  utils::Vector<double> vec(any.size(), 0);
  int n=0;
  for (int i=0; i<mX.size(); ++i){
    if (n >= any.size()) break;

    if (mXtype[i] == PARAM_FREE || mXtype[i] == PARAM_FREE_WITH_LIMITS){
      // free parameter
      vec[n] = any[i];
      ++n;
    }
  }
  vec.resize(n);
  return vec;
}


utils::Vector<parametertype> Param::map_any_as_Xfree(const utils::Vector<parametertype> & any){
  utils::Vector<parametertype> vec(any.size());
  int n=0;
  for (int i=0; i<mX.size(); ++i){
    if (n >= any.size()) break;

    if (mXtype[i] == PARAM_FREE || mXtype[i] == PARAM_FREE_WITH_LIMITS){
      // free parameter
      vec[n] = any[i];
      ++n;
    }
  }
  vec.resize(n);
  return vec;
}




