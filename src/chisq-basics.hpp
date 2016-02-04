


#ifndef CHISQ_BASICS_HPP
#define CHISQ_BASICS_HPP



#include <string>
#include <iostream>
#include <sstream>
#include <limits>

#include <cstdlib>
#include <cmath>


#include "utils.hpp"
#include "utils-math.hpp"
#include "utils-vector.hpp"
#include "utils-matrix.hpp"
#include "utils-errors.hpp"

#include "param.hpp"

#include "funcfit-errors.hpp"


/* ###################################################################
   ###################################################################

  Least-squares merit/cost function object

  Chi^2 = 0.5 * sum_k ( w_k * ( Y_k - MDY_k ) )^2

  w_k   = weight of property k
  u_k   = 1/w_k if uncertainties are used
  Weights are normalized by the finalize_setup() method.
  Y_k   = observed value of property k
  MDY_k = modelled value

  At least these sub-objects are defined:
  f_k  = w_k * ( Y_k - MDY_k )
  J_ki = df_k/dx_i = - w_k * d(MDY_k)/dx_i

  The current point X in parameter space is also a sub-object,
  through the Param type sub-object.

  Most methods have two alternatives: m() and m(p)
  - m() uses currently saved point X
  - m(p) basically sets X=p and then calls m()

  Complication: some methods do not need to update X permanently,
  e.g. only for calculating some measure at a trial point p_t.
  These methods are:
    - d2f_dx1dx2(x)
    - Hessian_exact(x)
    - Hessian(x)
  These methods first calls backup(), then updates X and performs
  calculations, updating MDY,f,J as needed. At the end, restore()
  is called, copying back previously backed up values for X,MDY,f,J.

  set_point(xp) => updates X, MDY,f,J expires
  f() => updates MDY if needed, and updates f
  J() => updates MDY if needed, and updates J
  gradient()



   Example 1:
   (0) Prepare trial step: backup()
   (1) func(p) => X updated, existing MDY,f,J expires, new MDY,f are calculated
   (2) gradient(p) => X not updated, since check notes it's same point, J is recalc.
   or
   (2') gradient() => J is recalculated
   (3) Trial failed: if new trial, goto (1) else goto (5)
   (4) Trial accepted, goto (6)
   (5) restore()
   (6) ...



   ###################################################################
   ################################################################### */



using utils::aborterror;
using utils::Vector;
using utils::Matrix;
using utils::fp_are_equal;
using std::numeric_limits;
using funcfit::bad_point;

enum object_status { expired, updated } ;



template <typename S, typename T, typename U>
class ChiSqFunc {


public:
  typedef Vector<U> (*F_model)(S &, Vector<T> &);
  // use as: typename utils::ChiSqFunc<S,T,U>::F_model
  typedef void (*F_report)(S &, Vector<T> &, Vector<U> &, Vector<U> &);



private:
  S         mParam;
  Vector<T> mDataX;
  Vector<U> mDataY;
  Vector<U> mDataScaleY;
  Vector<U> mDataUncertaintyY;
  Vector<U> mDataWeightY;
  // mModelFuncPointer evaluation: returns model Y data based on
  // parameter object and X data
  typename ChiSqFunc<S,T,U>::F_model  mModelFuncPointer;
  typename ChiSqFunc<S,T,U>::F_report mReportFuncPointer;
  Vector<U> mModelDataY;
  Vector<U> mModelDataY_bak;

  Vector<double> X_bak;
  Vector<double> mf;
  Vector<double> mf_bak;
  Matrix<double> mJ;
  Matrix<double> mJ_bak;

  double mbarrier_scale;
  double mbarrier_scale_bak;
  Vector<double> mV;
  Vector<double> mV_bak;
  Vector<double> mdV_dx;
  Vector<double> mdV_dx_bak;

  bool muse_barrier_rescaling;
  bool muse_scales;

  object_status mstatus_X;
  object_status mstatus_MDY;
  object_status mstatus_f;
  object_status mstatus_J;
  object_status mstatus_V;
  object_status mstatus_dV_dx;

  object_status mstatus_X_bak;
  object_status mstatus_MDY_bak;
  object_status mstatus_f_bak;
  object_status mstatus_J_bak;
  object_status mstatus_V_bak;
  object_status mstatus_dV_dx_bak;


public:
  ChiSqFunc();
  ~ChiSqFunc();

  ChiSqFunc(const S & P,
	    const Vector<T> & DX,
	    const Vector<U> & DY);



  // Copy constructor:
  ChiSqFunc(const ChiSqFunc & sv);

  // Operators:
  ChiSqFunc & operator=(const ChiSqFunc & sv);

  // VERY IMPORTANT OPERATOR ():
  double operator()(void);
  double operator()(const Vector<double> & xi);




  // Setting a parameter point explicitly:
  void set_point(const Vector<double> & xi);

  void reset(const Vector<double> & xi);
  void reset(void);

  void clear_data(void); // Removes DataX, DataY, DataUncertainty, DataWeightY, DataScaleY, ModelDataY

  void finalize_setup();

  void backup();
  void restore();

  void trial_do(const Vector<double> & xi);
  void trial_accept();
  void trial_reject();

  void debug();





  bool           point_is_good(void);
  bool           point_is_good(const Vector<double> & xi);

  bool           point_is_new(const Vector<double> & any);

  Vector<double> free_parameters(void);
  Vector<double> free_parameters(const Vector<double> & any);
  Vector<double> free_parameters_part(const Vector<double> & any);
  bool           free_parameters_part_is_new(const Vector<double> & xifree);

  Vector<double> all_parameters(void);
  Vector<double> all_parameters(const Vector<double> & any);
  Vector<double> nonfree_parameters_part(const Vector<double> & any);
  bool           nonfree_parameters_part_is_new(const Vector<double> & xinonfree);

  Vector<double>        map_vector_as_free_parameters(const Vector<double> & any);
  Vector<parametertype> map_vector_as_free_parameters(const Vector<parametertype> & any);


  double & barrier_scale(void);
  bool   & use_barrier_rescaling(void);

  bool & use_scales(void);

  S & Param(void);
  Vector<T> & DataX(void);
  Vector<U> & DataY(void);
  Vector<U> & DataUncertaintyY(void);
  Vector<U> & DataWeightY(void);

  // Set/get function pointers
  typename ChiSqFunc<S,T,U>::F_model  & ModelFuncPointer(void);
  typename ChiSqFunc<S,T,U>::F_report & ReportFuncPointer(void);

  // Set/get ModelDataY
  Vector<U> & ModelDataY(void);
  Vector<U>   ModelDataY(const Vector<double> & xi);


  void report_on_parameters_and_data();


  Vector<double>   f(const Vector<double> & xi);
  Vector<double> & f(void);
  void             calc_f(void);

  Matrix<double>   J(const Vector<double> & xi);
  Matrix<double> & J(void);

  Vector<double>   V(const Vector<double> & xi);
  Vector<double> & V(void);
  void             calc_V(void);

  Vector<double>   dV_dx(const Vector<double> & xi);
  Vector<double> & dV_dx(void);
  void             calc_dV_dx(void);

  Matrix<double>   df_dx(const Vector<double> & xi);
  Matrix<double> & df_dx(void);
  void             calc_df_dx(void);

  Vector < Matrix<double> > d2f_dx1dx2(const Vector<double> & xi);
  Vector < Matrix<double> > d2f_dx1dx2(void);

  Matrix<double> Hessian_exact(const Vector<double> & xi);
  Matrix<double> Hessian_exact(void);

  Matrix<double> Hessian(const Vector<double> & xi);
  Matrix<double> Hessian(void);




  double value(const Vector<double> & xi);
  double value(void);
  double value_data(const Vector<double> & xi);
  double value_data(void);
  double value_barrier_raw(const Vector<double> & xi);
  double value_barrier_raw(void);
  double adjust_barrier(const Vector<double> & xi);
  double adjust_barrier(void);
  double value_barrier(const Vector<double> & xi);
  double value_barrier(void);


  Vector<double> gradient(const Vector<double> & xi);
  Vector<double> gradient(void);
  Vector<double> gradient_data(const Vector<double> & xi);
  Vector<double> gradient_data(void);
  Vector<double> gradient_barrier(const Vector<double> & xi);
  Vector<double> gradient_barrier(void);










} ;



// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Default ctor:
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template <typename S, typename T, typename U>
ChiSqFunc<S,T,U> :: ChiSqFunc(){
  mModelFuncPointer = 0;
  mReportFuncPointer = 0;
  mbarrier_scale = 1.0;
  muse_barrier_rescaling = true;
  muse_scales = true;
  mstatus_X   = expired;
  mstatus_MDY = expired;
  mstatus_f   = expired;
  mstatus_J   = expired;
  mstatus_V     = expired;
  mstatus_dV_dx = expired;
}



// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Default dtor:
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template <typename S, typename T, typename U>
ChiSqFunc<S,T,U> :: ~ChiSqFunc(){
  mModelFuncPointer = 0;
  mReportFuncPointer = 0;
}



// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Alternative ctors:
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template <typename S, typename T, typename U>
ChiSqFunc<S,T,U> :: ChiSqFunc(const S & P,
			      const Vector<T> & DX,
			      const Vector<U> & DY){
  mParam = P;
  mDataX = DX;
  mDataY = DY;
  mModelFuncPointer = 0;
  mReportFuncPointer = 0;
  mbarrier_scale = 1.0;
  muse_barrier_rescaling = true;
  muse_scales = true;
  mstatus_X   = expired;
  mstatus_MDY = expired;
  mstatus_f   = expired;
  mstatus_J   = expired;
  mstatus_V     = expired;
  mstatus_dV_dx = expired;
}




// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Copy constructor:
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template <typename S, typename T, typename U>
ChiSqFunc<S,T,U> :: ChiSqFunc(const ChiSqFunc<S,T,U> & sv){
  mParam = sv.mParam;
  mDataX = sv.mDataX;
  mDataY = sv.mDataY;
  mDataScaleY = sv.mDataScaleY;
  mDataUncertaintyY = sv.mDataUncertaintyY;
  mDataWeightY = sv.mDataWeightY;
  mModelFuncPointer  = sv.mModelFuncPointer;
  mReportFuncPointer  = sv.mReportFuncPointer;
  mModelDataY     = sv.mModelDataY;
  mModelDataY_bak = sv.mModelDataY_bak;
  mf     = sv.mf;
  mf_bak = sv.mf_bak;
  X_bak  = sv.X_bak;
  mJ     = sv.mJ;
  mJ_bak = sv.mJ_bak;
  mbarrier_scale = sv.mbarrier_scale;
  muse_barrier_rescaling = sv.muse_barrier_rescaling;
  muse_scales = sv.muse_scales;
  mstatus_X   = sv.mstatus_X;
  mstatus_MDY = sv.mstatus_MDY;
  mstatus_f   = sv.mstatus_f;
  mstatus_J   = sv.mstatus_J;
  mstatus_X_bak   = sv.mstatus_X_bak;
  mstatus_MDY_bak = sv.mstatus_MDY_bak;
  mstatus_f_bak   = sv.mstatus_f_bak;
  mstatus_J_bak   = sv.mstatus_J_bak;
  mstatus_V_bak     = sv.mstatus_V_bak;
  mstatus_dV_dx_bak = sv.mstatus_dV_dx_bak;
}





// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Operators:
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template <typename S, typename T, typename U>
ChiSqFunc<S,T,U> & ChiSqFunc<S,T,U> ::operator=(const ChiSqFunc<S,T,U> & sv){

  if (this==&sv) return *this;

  mParam = sv.mParam;
  mDataX = sv.mDataX;
  mDataY = sv.mDataY;
  mDataScaleY = sv.mDataScaleY;
  mDataUncertaintyY = sv.mDataUncertaintyY;
  mDataWeightY = sv.mDataWeightY;
  mModelFuncPointer  = sv.mModelFuncPointer;
  mReportFuncPointer  = sv.mReportFuncPointer;
  mModelDataY     = sv.mModelDataY;
  mModelDataY_bak = sv.mModelDataY_bak;
  mf     = sv.mf;
  mf_bak = sv.mf_bak;
  X_bak  = sv.X_bak;
  mJ     = sv.mJ;
  mJ_bak = sv.mJ_bak;
  mbarrier_scale = sv.mbarrier_scale;
  muse_barrier_rescaling = sv.muse_barrier_rescaling;
  muse_scales = sv.muse_scales;
  mstatus_X   = sv.mstatus_X;
  mstatus_MDY = sv.mstatus_MDY;
  mstatus_f   = sv.mstatus_f;
  mstatus_J   = sv.mstatus_J;
  mstatus_X_bak   = sv.mstatus_X_bak;
  mstatus_MDY_bak = sv.mstatus_MDY_bak;
  mstatus_f_bak   = sv.mstatus_f_bak;
  mstatus_J_bak   = sv.mstatus_J_bak;
  mstatus_V_bak     = sv.mstatus_V_bak;
  mstatus_dV_dx_bak = sv.mstatus_dV_dx_bak;

  return *this;
}



template <typename S, typename T, typename U>
double ChiSqFunc<S,T,U> ::operator()(void){
  return value();
}

template <typename S, typename T, typename U>
double ChiSqFunc<S,T,U> ::operator()(const Vector<double> & xi){

  if (point_is_new(xi)){
    mParam.Xupdate(xi);
    mstatus_X   = updated;
    mstatus_MDY = expired;
    mstatus_f   = expired;
    mstatus_J   = expired;
    mstatus_V     = expired;
    mstatus_dV_dx = expired;
  }
  return value();
}



// ##############################################################################
// ##############################################################################

template <typename S, typename T, typename U>
bool ChiSqFunc<S,T,U> ::point_is_new(const Vector<double> & any){

  if (any.size()==0) return false;

  if (mstatus_X == expired) return true;

  return ( nonfree_parameters_part_is_new( nonfree_parameters_part(any) )
	   || free_parameters_part_is_new( free_parameters_part(any) ) );
}

// ##############################################################################
// ##############################################################################

template <typename S, typename T, typename U>
Vector<double> ChiSqFunc<S,T,U> ::free_parameters(void){
  return mParam.Xfree();
}

// Return the free parameters in input vector, update X if free parameters are new:
template <typename S, typename T, typename U>
Vector<double> ChiSqFunc<S,T,U> ::free_parameters(const Vector<double> & any){

  if (any.size()==0) return mParam.Xfree();

  if (free_parameters_part_is_new( free_parameters_part(any) )){
    // Input free parameters are new. Update.
    mParam.Xupdate(any);

    mstatus_X   = updated;
    mstatus_MDY = expired;
    mstatus_f   = expired;
    mstatus_J   = expired;
    mstatus_V     = expired;
    mstatus_dV_dx = expired;
  }

  return mParam.Xfree(any);
}

// Return the free parameters in input vector:
template <typename S, typename T, typename U>
Vector<double> ChiSqFunc<S,T,U> ::free_parameters_part(const Vector<double> & any){
  int NXin   = any.size();
  int NXfree = mParam.Xfree().size();
  int N = ( (NXin == 0) ? 0 : NXfree );
  Vector<double> Xfree_any(N, 0);

  if (N==0) return Xfree_any;

  // ---------------------------------------------------------------------
  // Extract the free parameters from input vector
  // ---------------------------------------------------------------------
  if (NXin == NXfree){
    // Input vector contains only free parameters, not the other ones.
    for (int i=0; i<NXfree; ++i)
      Xfree_any[i] = any[i];
  }
  else {
    // Input vector contains all parameters. Extract the free ones.
    int NXall = mParam.X().size();
    int n=0;
    for (int i=0; i<NXin; ++i){
      if (n >= NXfree) break;

      if (mParam.Xtype(i) == PARAM_FREE || mParam.Xtype(i) == PARAM_FREE_WITH_LIMITS){
	Xfree_any[n] = any[i]; // free parameter
	++n;
      }
    }
  }
  return Xfree_any;
}

// Check if free parameters in vector are new or not:
template <typename S, typename T, typename U>
bool ChiSqFunc<S,T,U> ::free_parameters_part_is_new(const Vector<double> & xifree){
  int NXfree = mParam.Xfree().size();
  Vector<double> Xfree = mParam.Xfree();

  if (xifree.size() == 0)      return false;
  if (xifree.size() != NXfree) return false;

  // ---------------------------------------------------------------------
  // Compare old and new free parameters
  // ---------------------------------------------------------------------
  double eps = numeric_limits<double>::epsilon(), td, tdsqsum=0;
  int nc=0;
  for (int ii=0; ii<NXfree; ++ii){
    if (fp_are_equal(Xfree[ii], xifree[ii])) nc++;
    /*
    td = Xfree[ii] - xifree[ii];
    td = (td < 0) ? -td : td;
    tdsqsum += td*td;
    */
  }
  /*
  if (tdsqsum > eps) return true;
  else return false;
  */
  if (nc==NXfree) return false;
  else            return true;
}





// ##############################################################################
// ##############################################################################


template <typename S, typename T, typename U>
Vector<double> ChiSqFunc<S,T,U> ::all_parameters(void){
  return mParam.X();
}

// Return all parameters, update X if free parameters are new:
template <typename S, typename T, typename U>
Vector<double> ChiSqFunc<S,T,U> ::all_parameters(const Vector<double> & any){

  if (any.size()==0) return mParam.X();

  // Extract non-free and free parameters.
  if ( (nonfree_parameters_part_is_new( nonfree_parameters_part(any) ))
       ||
       (free_parameters_part_is_new( free_parameters_part(any) )) ){
    mParam.Xupdate(any);

    mstatus_X   = updated;
    mstatus_MDY = expired;
    mstatus_f   = expired;
    mstatus_J   = expired;
    mstatus_V     = expired;
    mstatus_dV_dx = expired;
  }
  return mParam.X(any);
}

// Return the non-free parameters in input vector:
template <typename S, typename T, typename U>
Vector<double> ChiSqFunc<S,T,U> ::nonfree_parameters_part(const Vector<double> & any){
  int NXin   = any.size();
  int NXfree = mParam.Xfree().size();
  int NXall  = mParam.X().size();
  int NXnonfree = NXall - NXfree;

  int N = ( (NXin == 0) ? 0 : NXnonfree );
  Vector<double> Xnonfree_any(N, 0);

  if (N==0) return Xnonfree_any;

  // ---------------------------------------------------------------------
  // Extract the non-free parameters from input vector
  // ---------------------------------------------------------------------
  if (NXin == NXfree){
    // Input vector contains only free parameters. Return null vector.
    return Xnonfree_any;
  }
  else {
    // Input vector contains all parameters. Extract the non-free ones.
    int n=0;
    for (int i=0; i<NXin; ++i){
      if (n >= NXnonfree) break;

      if ( ! (mParam.Xtype(i) == PARAM_FREE || mParam.Xtype(i) == PARAM_FREE_WITH_LIMITS) ){
	Xnonfree_any[n] = any[i]; // non-free parameter
	++n;
      }
    }
  }
  return Xnonfree_any;
}

// Check if non-free parameters in vector are new or not:
template <typename S, typename T, typename U>
bool ChiSqFunc<S,T,U> ::nonfree_parameters_part_is_new(const Vector<double> & xinonfree){
  int NXall  = mParam.X().size();
  int NXfree = mParam.Xfree().size();
  int NXnonfree = NXall - NXfree;
  Vector<double> Xall  = mParam.X();
  Vector<double> Xfree = mParam.Xfree();
  Vector<double> Xnonfree(NXnonfree, 0);

  if (xinonfree.size() == 0)         return false;
  if (xinonfree.size() != NXnonfree) return false;

  // Collect all non-free parameters
  int n=0;
  for (int i=0; i<NXall; ++i){
    if (n >= NXnonfree) break;

    if ( ! (mParam.Xtype(i) == PARAM_FREE || mParam.Xtype(i) == PARAM_FREE_WITH_LIMITS) ){
      Xnonfree[n] = Xall[i]; // non-free parameter
      ++n;
    }
  }

  // ---------------------------------------------------------------------
  // Compare old and new non-free parameters
  // ---------------------------------------------------------------------
  double eps = numeric_limits<double>::epsilon(), td, tdsqsum=0;
  int nc=0;
  for (int ii=0; ii<NXnonfree; ++ii){
    if (fp_are_equal(Xnonfree[ii], xinonfree[ii])) nc++;
    /*
    td = Xnonfree[ii] - xinonfree[ii];
    td = (td < 0) ? -td : td;
    tdsqsum += td*td;
    */
  }
  /*
  if (tdsqsum > eps) return true;
  else return false;
  */
  if (nc==NXnonfree) return false;
  else               return true;
}


// ##############################################################################
// ##############################################################################




template <typename S, typename T, typename U>
Vector<double> ChiSqFunc<S,T,U> ::map_vector_as_free_parameters(const Vector<double> & any){
  return mParam.map_any_as_Xfree(any);
}

template <typename S, typename T, typename U>
Vector<parametertype> ChiSqFunc<S,T,U> ::map_vector_as_free_parameters(const Vector<parametertype> & any){
  return mParam.map_any_as_Xfree(any);
}




template <typename S, typename T, typename U>
bool ChiSqFunc<S,T,U> ::point_is_good(void){
  return mParam.Xfree_is_good();
}

template <typename S, typename T, typename U>
bool ChiSqFunc<S,T,U> ::point_is_good(const Vector<double> & xi){
  if (point_is_new(xi)){
    mParam.Xupdate(xi);
    mstatus_X   = updated;
    mstatus_MDY = expired;
    mstatus_f   = expired;
    mstatus_J   = expired;
    mstatus_V     = expired;
    mstatus_dV_dx = expired;
  }
  return mParam.Xfree_is_good();
}






template <typename S, typename T, typename U>
double &  ChiSqFunc<S,T,U> ::barrier_scale(void){
  return mbarrier_scale;
}

template <typename S, typename T, typename U>
bool &  ChiSqFunc<S,T,U> ::use_barrier_rescaling(void){
  return muse_barrier_rescaling;
}





template <typename S, typename T, typename U>
bool & ChiSqFunc<S,T,U> ::use_scales(void){
  return muse_scales;
}



// ############################################################################
// ############################################################################



// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Set/get parameter object
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename S, typename T, typename U>
S & ChiSqFunc<S,T,U> ::Param(void){
  return mParam;
}



// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Set/get x-data
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename S, typename T, typename U>
Vector<T> & ChiSqFunc<S,T,U> ::DataX(void){
  return mDataX;
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Set/get y-data
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename S, typename T, typename U>
Vector<U> & ChiSqFunc<S,T,U> ::DataY(void){
  return mDataY;
}







// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Set/get errors-y-data
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename S, typename T, typename U>
Vector<U> & ChiSqFunc<S,T,U> ::DataUncertaintyY(void){
  return mDataUncertaintyY;
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Set/get weights-y-data
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename S, typename T, typename U>
Vector<U> & ChiSqFunc<S,T,U> ::DataWeightY(void){
  return mDataWeightY;
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Set/get function pointers
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename S, typename T, typename U>
typename ChiSqFunc<S,T,U>::F_model  & ChiSqFunc<S,T,U> ::ModelFuncPointer(void){
  return mModelFuncPointer;
}


template <typename S, typename T, typename U>
typename ChiSqFunc<S,T,U>::F_report & ChiSqFunc<S,T,U> ::ReportFuncPointer(void){
  return mReportFuncPointer;
}





// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Set/get ModelDataY
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename S, typename T, typename U>
Vector<U> ChiSqFunc<S,T,U> ::ModelDataY(const Vector<double> & xi){
  if (point_is_new(xi)){
    mParam.Xupdate(xi);
    mstatus_X   = updated;
    mstatus_MDY = expired;
    mstatus_f   = expired;
    mstatus_J   = expired;
    mstatus_V     = expired;
    mstatus_dV_dx = expired;
  }
  return ModelDataY();
}

template <typename S, typename T, typename U>
Vector<U> & ChiSqFunc<S,T,U> ::ModelDataY(void){
  if (mstatus_MDY==expired){
    mModelDataY = (mModelFuncPointer)(mParam, mDataX);
    mstatus_MDY=updated;
  }
  return mModelDataY;
}




// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename S, typename T, typename U>
void ChiSqFunc<S,T,U> ::report_on_parameters_and_data(){
  if (mModelDataY.size()==0 || mstatus_MDY==expired){
    mModelDataY = (mModelFuncPointer)(mParam, mDataX);
    mstatus_MDY = updated;
  }

  if (mReportFuncPointer != 0)
    (mReportFuncPointer)(mParam, mDataX, mDataY, mModelDataY);
}











// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Set/get f vector
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename S, typename T, typename U>
Vector<double> ChiSqFunc<S,T,U> ::f(const Vector<double> & xi){
  if (point_is_new(xi)){
    mParam.Xupdate(xi);
    mstatus_X   = updated;
    mstatus_MDY = expired;
    mstatus_f   = expired;
    mstatus_J   = expired;
    mstatus_V     = expired;
    mstatus_dV_dx = expired;
  }
  return f();
}

template <typename S, typename T, typename U>
Vector<double> & ChiSqFunc<S,T,U> ::f(void){
  if (mstatus_f==expired){
    calc_f();
    mstatus_f=updated;
  }
  return mf;
}

template <typename S, typename T, typename U>
void ChiSqFunc<S,T,U> ::calc_f(void){

  int N = mDataY.size();
  if (N==0){
    std::cout << "f calculation error: DataY are not set. Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
  if (mDataScaleY.size()!=N){
    std::cout << "f calculation error: DataScaleY are not set. Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
  if (mDataUncertaintyY.size()!=N && mDataWeightY.size()!=N){
    std::cout << "f calculation error: DataUncertaintyY or DataWeightY are not set. Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }

  if (mf.size() != N) mf.resize(N);

  mModelDataY = ModelDataY();


  for (int i=0; i!=N; ++i){
    double td = double(mDataWeightY[i]);
    if (mDataWeightY[i]<0) td = 1.0/double(mDataUncertaintyY[i]);
    
    mf[i] = td * (double(mDataY[i]) - double(mModelDataY[i])) / double(mDataScaleY[i]);
    
    /*
      cout << "data point " << i << "  uncert " << double(mDataUncertaintyY[i])
      << "  inputdata " << double(mDataY[i])
      << "  preddata " << double(mModelDataY[i])
      << "  datascale " << double(mDataScaleY[i])
      << endl;
    */
    
  }
  
  mstatus_f = updated;
  return;
}











template <typename S, typename T, typename U>
Vector<double>   ChiSqFunc<S,T,U> ::V(const Vector<double> & xi){
  if (point_is_new(xi)){
    mParam.Xupdate(xi);
    mstatus_X   = updated;
    mstatus_MDY = expired;
    mstatus_f   = expired;
    mstatus_J   = expired;
    mstatus_V     = expired;
    mstatus_dV_dx = expired;
  }
  return V();
}

template <typename S, typename T, typename U>
Vector<double> & ChiSqFunc<S,T,U> ::V(void){
  if (mstatus_V==expired){
    calc_V();
    mstatus_V=updated;
  }
  return mV;
}

template <typename S, typename T, typename U>
void             ChiSqFunc<S,T,U> ::calc_V(void){
  // ************************************************************
  // ************************************************************
  // Barrier part: Check all parameters!
  // ************************************************************
  // ************************************************************
  int N=0;
  for (int i=0; i< mParam.X().size(); ++i){
    if (mParam.X_is_unconstrained(i)) continue;
    if (mParam.X_is_fixed(i)) continue;
    N++;
  }
  mV.resize(N);
  if (N==0) return;

  double xi, ai, bi, ci;
  double eps = std::numeric_limits<double>::epsilon();
  double tbs = mbarrier_scale;
  //std::cout << "calc_V(): barrier_scale is " << tbs << std::endl;
  //std::cout << "calc_V(): eps is " << eps << std::endl;
  if (tbs < 0.0) tbs *= -1.0; 
  if (tbs < eps) return;
  funcfit::bad_point ebp;

  int n=0;
  for (int i=0; i< mParam.X().size(); ++i){
    if (mParam.X_is_unconstrained(i)) continue;
    if (mParam.X_is_fixed(i)) continue;

    xi = mParam.X(i);
    ai = mParam.Xmin(i);
    bi = mParam.Xmax(i);

    if (xi<=ai || xi>=bi){
      std::cout << "chisq-basics: calc_V: Warning: Parameter value " << xi
		<< " is outside limits (" << ai << ", " << bi << ")."
		<< " Throwing exception." << std::endl;
      throw ebp;
    }
    /* Penalty/Barrier function is logarithm based:
       U = - mu * sum_i p(xi,ai,bi)
       where
       p(xi,ai,bi) = log(bi-xi) + log(xi-ai) - 2*log(0.5*(bi-ai))
       dp/dxi = -1/(bi-xi) + 1/(xi-ai) = (-xi+ai + bi-xi)/((bi-xi)*(xi-ai))
       = (bi+ai - 2*xi)/((bi-xi)*(xi-ai))
    */
    ci = log(bi-xi) + log(xi-ai) - 2.0*log(0.5*(bi-ai));
    //std::cout << "n = " << n << " constraint (without -mu factor) is " << ci << std::endl;
    mV[n] = ci;
    n++;
    //barrier += -mbarrier_scale * ci;
  }
  // ************************************************************
  // ************************************************************
}

template <typename S, typename T, typename U>
Vector<double>   ChiSqFunc<S,T,U> ::dV_dx(const Vector<double> & xi){
  if (point_is_new(xi)){
    mParam.Xupdate(xi);
    mstatus_X   = updated;
    mstatus_MDY = expired;
    mstatus_f   = expired;
    mstatus_J   = expired;
    mstatus_V     = expired;
    mstatus_dV_dx = expired;
  }
  return dV_dx();
}

template <typename S, typename T, typename U>
Vector<double> & ChiSqFunc<S,T,U> ::dV_dx(void){
  if (mstatus_dV_dx==expired){
    calc_dV_dx();
    mstatus_dV_dx=updated;
  }
  return mdV_dx;
}

template <typename S, typename T, typename U>
void             ChiSqFunc<S,T,U> ::calc_dV_dx(void){
  // ************************************************************
  // ************************************************************
  // Barrier part: Check all parameters!
  // ************************************************************
  // ************************************************************
  int N=0;
  for (int i=0; i< mParam.X().size(); ++i){
    if (mParam.X_is_unconstrained(i)) continue;
    if (mParam.X_is_fixed(i)) continue;
    N++;
  }
  mdV_dx.resize(N);
  if (N==0) return;

  double xi, ci, ai, bi, si, xti, dci_dxi;
  double eps = std::numeric_limits<double>::epsilon();
  double tbs = mbarrier_scale;
  if (tbs < 0.0) tbs *= -1.0; 
  if (tbs < eps) return;
  funcfit::bad_point ebp;

  int n=0;
  for (int i=0; i< mParam.X().size(); ++i){
    if (mParam.X_is_unconstrained(i)) continue;
    if (mParam.X_is_fixed(i)) continue;

    xi = mParam.X(i);
    ai = mParam.Xmin(i);
    bi = mParam.Xmax(i);

    if (xi<=ai || xi>=bi){
      std::cout << "chisq-basics: calc_dV_dx: Warning: Parameter value " << xi
		<< " is outside limits (" << ai << ", " << bi << ")."
		<< " Throwing exception." << std::endl;
      throw ebp;
    }
    /* Penalty/Barrier function is logarithm based:
       U = - mu * sum_i p(xi,ai,bi)
       where
       p(xi,ai,bi) = log(bi-xi) + log(xi-ai) - 2*log(0.5*(bi-ai))
       dp/dxi = -1/(bi-xi) + 1/(xi-ai) = (-xi+ai + bi-xi)/((bi-xi)*(xi-ai))
       = (bi+ai - 2*xi)/((bi-xi)*(xi-ai))
    */
    dci_dxi = -1.0/(bi-xi) + 1.0/(xi-ai);
    mdV_dx[n] = dci_dxi;
    //Xt[n] = - mbarrier_scale * dci_dxi;
    n++;
  }
  // ************************************************************
  // ************************************************************
}











// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Set/get J = df_dx matrix
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename S, typename T, typename U>
Matrix<double> ChiSqFunc<S,T,U> ::J(const Vector<double> & xi){
  return df_dx(xi);
}

template <typename S, typename T, typename U>
Matrix<double> & ChiSqFunc<S,T,U> ::J(void){
  return df_dx();
}



// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Set/get df_dx = J matrix
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename S, typename T, typename U>
Matrix<double> ChiSqFunc<S,T,U> ::df_dx(const Vector<double> & xi){

  if (point_is_new(xi)){
    mParam.Xupdate(xi);
    mstatus_X   = updated;
    mstatus_MDY = expired;
    mstatus_f   = expired;
    mstatus_J   = expired;
    mstatus_V     = expired;
    mstatus_dV_dx = expired;
  }
  return df_dx();
}

template <typename S, typename T, typename U>
Matrix<double> & ChiSqFunc<S,T,U> ::df_dx(void){
  if (mstatus_J==expired){
    calc_df_dx();
    mstatus_J=updated;
  }
  return mJ;
}

template <typename S, typename T, typename U>
void ChiSqFunc<S,T,U> ::calc_df_dx(void){
  int N = mDataY.size();
  Vector<double> Xf = mParam.Xfree();
  int NXf = Xf.size();

  if (mJ.nrows() != N || mJ.ncols() != NXf) mJ.resize(N, NXf);

  if (mDataUncertaintyY.size()!=N && mDataWeightY.size()!=N){
    std::cout << "f calculation error: DataUncertaintyY or DataWeightY are not set. Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
  if (mDataScaleY.size()!=N){
    std::cout << "J calculation error: DataScaleY are not set. Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }

  double eps = utils::eps_d(), tmpx, bakxi, dx;
  Vector<U> MDataY1, MDataY2;

  for (int ix=0; ix<NXf; ++ix){
    bakxi = Xf[ix];
    tmpx  = utils::absval(bakxi);

    dx = pow(eps, 1.0/3.0);
    if (tmpx>eps) dx = tmpx * dx;


    bool do_inc1=true, do_inc2=true;

    // Positive perturbation in parameter:
    try {
      Xf[ix] = bakxi + dx; mParam.Xupdate(Xf); // *** FAILS IF Xf[ix] CLOSE TO UPPER LIMIT !!! ***
    } catch (bad_point & e1){
      do_inc1 = false; Xf[ix] = bakxi; mParam.Xupdate(Xf);
    }
    if (do_inc1) MDataY1 = (mModelFuncPointer)(mParam, mDataX);
    else MDataY1 = ModelDataY();
    
    // Negative perturbation in parameter:
    try {
      Xf[ix] = bakxi - dx; mParam.Xupdate(Xf); // *** FAILS IF Xf[ix] CLOSE TO LOWER LIMIT !!! ***
    } catch (bad_point & e1){
      do_inc2 = false; Xf[ix] = bakxi; mParam.Xupdate(Xf);
    }
    if (do_inc2) MDataY2 = (mModelFuncPointer)(mParam, mDataX);
    else MDataY2 = ModelDataY();

    // Reset:
    Xf[ix] = bakxi;
    mParam.Xupdate(Xf);


    for (int i=0; i<N; ++i){
      double td = mDataWeightY[i];
      if (mDataWeightY[i]<0.0) td = 1.0/mDataUncertaintyY[i];

      double den, num;
      den = double(MDataY1[i]) - double(MDataY2[i]);

      if   (do_inc1 && do_inc2) num = 2*dx;
      else num = dx;

      mJ.elem(i, ix) = - td * (den / num) / double(mDataScaleY[i]);
    }

  }

  mstatus_J = updated;

  return;
}






// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Get d2f_dx1dx2 matrix (f is a vector)
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename S, typename T, typename U>
Vector < Matrix<double> > ChiSqFunc<S,T,U> ::d2f_dx1dx2(const Vector<double> & xi){
  if (point_is_new(xi)){
    mParam.Xupdate(xi);
    mstatus_X   = updated;
    mstatus_MDY = expired;
    mstatus_f   = expired;
    mstatus_J   = expired;
    mstatus_V     = expired;
    mstatus_dV_dx = expired;
  }
  return d2f_dx1dx2();
}


template <typename S, typename T, typename U>
Vector < Matrix<double> > ChiSqFunc<S,T,U> ::d2f_dx1dx2(void){
  int N = mDataY.size();
  Vector<double> Xf = mParam.Xfree();
  int NXf = Xf.size();

  if (mDataUncertaintyY.size()!=N && mDataWeightY.size()!=N){
    std::cout << "f calculation error: DataUncertaintyY or DataWeightY are not set. Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
  if (mDataScaleY.size()!=N){
    std::cout << "d2f_dx1dx2 calculation error: DataScaleY are not set. Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }

  Vector < Matrix<double> > fpp_x1x2;
  double tmpx1, tmpx2, bakx1, bakx2, dx1, dx2;
  Vector<U> MDataY1, MDataY2, MDataY3, MDataY4;
  double eps = utils::eps_d();

  bool do_pt1, do_pt2;
  bool do_inc1, do_inc2, do_inc3, do_inc4;
  bool ok_pos_step1=true, ok_pos_step2=true;
  bool ok_neg_step1=true, ok_neg_step2=true;


  fpp_x1x2.resize(N);
  for (int i=0; i<N; ++i)
    fpp_x1x2[i].resize(NXf, NXf);


  for (int i=0; i<NXf; ++i){
    // Use symmetry:
    for (int j=i+1; j<NXf; ++j){
      if (j>=NXf) continue;

      // First order partial derivatives:

      bakx1 = Xf[i];
      bakx2 = Xf[j];
      tmpx1 = utils::absval(bakx1);
      tmpx2 = utils::absval(bakx2);

      dx1 = pow(eps, 1.0/3.0);
      if (tmpx1>eps) dx1 = tmpx1 * dx1;
      dx2 = pow(eps, 1.0/3.0);
      if (tmpx2>eps) dx2 = tmpx2 * dx2;
      

      ok_pos_step1=true; ok_pos_step2=true;
      ok_neg_step1=true; ok_neg_step2=true;

      // -------------------------------------------------------------------------------
      // Perturbation 1 in parameters:
      do_pt1 = do_pt2 = true;
      try {
	Xf[i] = bakx1 + dx1; mParam.Xupdate(Xf);
      } catch (bad_point & e1){
	do_pt1 = false; Xf[i] = bakx1; mParam.Xupdate(Xf); ok_pos_step1=false;
      }
      try {
	Xf[j] = bakx2 + dx2; mParam.Xupdate(Xf);
      } catch (bad_point & e1){
	do_pt2 = false; Xf[j] = bakx2; mParam.Xupdate(Xf); ok_pos_step2=false;
      }
      if (do_pt1 && do_pt2) do_inc1=true; else do_inc1=false;
      if (!do_inc1) MDataY1 = ModelDataY();
      else MDataY1 = (mModelFuncPointer)(mParam, mDataX);


      // -------------------------------------------------------------------------------
      // Perturbation 2 in parameters:
      do_pt1 = do_pt2 = true;
      try {
	Xf[i] = bakx1 + dx1; mParam.Xupdate(Xf);
      } catch (bad_point & e1){
	do_pt1 = false; Xf[i] = bakx1; mParam.Xupdate(Xf); ok_pos_step1=false;
      }
      try {
	Xf[j] = bakx2 - dx2; mParam.Xupdate(Xf);
      } catch (bad_point & e1){
	do_pt2 = false; Xf[j] = bakx2; mParam.Xupdate(Xf); ok_neg_step2=false;
      }
      if (do_pt1 && do_pt2) do_inc2=true; else do_inc2=false;
      if (!do_inc2) MDataY2 = ModelDataY();
      else MDataY2 = (mModelFuncPointer)(mParam, mDataX);


      // -------------------------------------------------------------------------------
      // Perturbation 3 in parameters:
      do_pt1 = do_pt2 = true;
      try {
	Xf[i] = bakx1 - dx1; mParam.Xupdate(Xf);
      } catch (bad_point & e1){
	do_pt1 = false; Xf[i] = bakx1; mParam.Xupdate(Xf); ok_neg_step1=false;
      }
      try {
	Xf[j] = bakx2 + dx2; mParam.Xupdate(Xf);
      } catch (bad_point & e1){
	do_pt2 = false; Xf[j] = bakx2; mParam.Xupdate(Xf); ok_pos_step2=false;
      }
      if (do_pt1 && do_pt2) do_inc3=true; else do_inc3=false;
      if (!do_inc3) MDataY3 = ModelDataY();
      else MDataY3 = (mModelFuncPointer)(mParam, mDataX);


      // -------------------------------------------------------------------------------
      // Perturbation 4 in parameters:
      do_pt1 = do_pt2 = true;
      try {
	Xf[i] = bakx1 - dx1; mParam.Xupdate(Xf);
      } catch (bad_point & e1){
	do_pt1 = false; Xf[i] = bakx1; mParam.Xupdate(Xf); ok_neg_step1=false;
      }
      try {
	Xf[j] = bakx2 - dx2; mParam.Xupdate(Xf);
      } catch (bad_point & e1){
	do_pt2 = false; Xf[j] = bakx2; mParam.Xupdate(Xf); ok_neg_step2=false;
      }
      if (do_pt1 && do_pt2) do_inc4=true; else do_inc4=false;
      if (!do_inc4) MDataY4 = ModelDataY();
      else MDataY4 = (mModelFuncPointer)(mParam, mDataX);


      // -------------------------------------------------------------------------------

      // Reset:
      Xf[i] = bakx1;
      Xf[j] = bakx2;
      mParam.Xupdate(Xf);

      for (int k=0; k!=N; ++k){
	double td = double(mDataWeightY[k]);
	if (mDataWeightY[k]<0) td = 1.0/double(mDataUncertaintyY[k]);

	double num1=1, num2=1;
	double den, num;

	if (ok_pos_step1 && ok_neg_step1) num1 = 2.0*dx1;
	else num1 = dx1;
	if (ok_pos_step2 && ok_neg_step2) num2 = 2.0*dx2;
	else num2 = dx2;
	num = num1*num2;

	den = double(MDataY1[k]) - double(MDataY2[k]) - double(MDataY3[k]) + double(MDataY4[k]);

	fpp_x1x2[k].elem(i,j) = - td * ( den / num ) / double(mDataScaleY[k]);
	  
	// Symmetry:
	fpp_x1x2[k].elem(j,i) = fpp_x1x2[k].elem(i,j);
      }

    }
  }

  mModelDataY = ModelDataY();


  // Second order derivative:
  for (int i=0; i!=NXf; ++i){
    bakx1 = Xf[i];
    tmpx1 = utils::absval(bakx1);

    dx1 = pow(eps, 1.0/3.0);
    if (tmpx1>eps) dx1 = tmpx1 * dx1;


    bool do_pdx1, do_pdx2, do_ndx1, do_ndx2;

    do_pdx1=true;
    try {
      Xf[i] = bakx1 + dx1; mParam.Xupdate(Xf);
    } catch (bad_point & e1){
      do_pdx1 = false; Xf[i] = bakx1; mParam.Xupdate(Xf);
    }
    do_pdx2=true;
    try {
      Xf[i] = bakx1 + 2*dx1; mParam.Xupdate(Xf);
    } catch (bad_point & e1){
      do_pdx2 = false; Xf[i] = bakx1; mParam.Xupdate(Xf);
    }
    do_ndx1=true;
    try {
      Xf[i] = bakx1 - dx1; mParam.Xupdate(Xf);
    } catch (bad_point & e1){
      do_ndx1 = false; Xf[i] = bakx1; mParam.Xupdate(Xf);
    }
    do_ndx2=true;
    try {
      Xf[i] = bakx1 - 2*dx1; mParam.Xupdate(Xf);
    } catch (bad_point & e1){
      do_ndx2 = false; Xf[i] = bakx1; mParam.Xupdate(Xf);
    }

    
    if (do_pdx1 && do_pdx2 && do_ndx1 && do_ndx2){
      Xf[i] = bakx1 + 2.0 * dx1; mParam.Xupdate(Xf);
      MDataY1 = (mModelFuncPointer)(mParam, mDataX);
      MDataY2 = ModelDataY();
      MDataY3 = ModelDataY();
      Xf[i] = bakx1 - 2.0 * dx1; mParam.Xupdate(Xf);
      MDataY4 = (mModelFuncPointer)(mParam, mDataX);
    }
    else if (do_pdx1 && do_pdx2 && do_ndx1 && !do_ndx2){
      Xf[i] = bakx1 + 2.0 * dx1; mParam.Xupdate(Xf);
      MDataY1 = (mModelFuncPointer)(mParam, mDataX);
      MDataY2 = ModelDataY();
      Xf[i] = bakx1 + 1.0 * dx1; mParam.Xupdate(Xf);
      MDataY3 = (mModelFuncPointer)(mParam, mDataX);
      Xf[i] = bakx1 - 2.0 * dx1; mParam.Xupdate(Xf);
      MDataY4 = (mModelFuncPointer)(mParam, mDataX);
    }
    else if (do_pdx1 && !do_pdx2 && do_ndx1 && do_ndx2){
      Xf[i] = bakx1 + 1.0 * dx1; mParam.Xupdate(Xf);
      MDataY1 = (mModelFuncPointer)(mParam, mDataX);
      Xf[i] = bakx1 - 1.0 * dx1; mParam.Xupdate(Xf);
      MDataY2 = (mModelFuncPointer)(mParam, mDataX);
      MDataY3 = ModelDataY();
      Xf[i] = bakx1 - 2.0 * dx1; mParam.Xupdate(Xf);
      MDataY4 = (mModelFuncPointer)(mParam, mDataX);
    }
    else if (do_pdx1 && do_pdx2 && !do_ndx1 && !do_ndx2){
      Xf[i] = bakx1 + 2.0 * dx1; mParam.Xupdate(Xf);
      MDataY1 = (mModelFuncPointer)(mParam, mDataX);
      Xf[i] = bakx1 + 1.0 * dx1; mParam.Xupdate(Xf);
      MDataY2 = (mModelFuncPointer)(mParam, mDataX);
      MDataY3 = MDataY2;
      MDataY4 = ModelDataY();
    }
    else if (!do_pdx1 && !do_pdx2 && do_ndx1 && do_ndx2){
      MDataY1 = ModelDataY();
      Xf[i] = bakx1 - 1.0 * dx1; mParam.Xupdate(Xf);
      MDataY2 = (mModelFuncPointer)(mParam, mDataX);
      MDataY3 = MDataY2;
      Xf[i] = bakx1 - 2.0 * dx1; mParam.Xupdate(Xf);
      MDataY4 = (mModelFuncPointer)(mParam, mDataX);
    }
    else {
      MDataY1 = ModelDataY();
      MDataY2 = ModelDataY();
      MDataY3 = ModelDataY();
      MDataY4 = ModelDataY();
    }

    // Reset:
    Xf[i] = bakx1;
    mParam.Xupdate(Xf);

    for (int k=0; k!=N; ++k){
      double td = double(mDataWeightY[k]);
      if (mDataWeightY[k]<0) td = 1.0/double(mDataUncertaintyY[k]);

      double den, num;
      den = double(MDataY1[k]) - double(MDataY2[k]) - double(MDataY3[k]) + double(MDataY4[k]);
      num = 4.0 * dx1 * dx1;
      if      ( do_pdx1 &&  do_pdx2 &&  do_ndx1 && !do_ndx2) num = 2.0*dx1*dx1;
      else if ( do_pdx1 && !do_pdx2 &&  do_ndx1 &&  do_ndx2) num = 2.0*dx1*dx1;
      else if ( do_pdx1 &&  do_pdx2 && !do_ndx1 && !do_ndx2) num = 2.0*dx1*dx1;
      else if (!do_pdx1 && !do_pdx2 &&  do_ndx1 &&  do_ndx2) num = 2.0*dx1*dx1;

      fpp_x1x2[k].elem(i,i) = - td * ( den / num ) / double(mDataScaleY[k]);
    }
  }

  return fpp_x1x2;
}







// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Calculate exact Hessian matrix
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename S, typename T, typename U>
Matrix<double> ChiSqFunc<S,T,U> ::Hessian_exact(const Vector<double> & xi){
  if (point_is_new(xi)){
    mParam.Xupdate(xi);
    mstatus_X   = updated;
    mstatus_MDY = expired;
    mstatus_f   = expired;
    mstatus_J   = expired;
    mstatus_V     = expired;
    mstatus_dV_dx = expired;
  }
  return Hessian_exact();
}


template <typename S, typename T, typename U>
Matrix<double> ChiSqFunc<S,T,U> ::Hessian_exact(void){
  int N = mDataY.size();
  int NXf = mParam.NXfree();
  Matrix<double> hessian(NXf, NXf, 0);

  hessian = Hessian();

  Vector < Matrix<double> > fpp_x1x2;
  mf = f();
  fpp_x1x2 = d2f_dx1dx2();

  for (int i=0; i<NXf; ++i){
    for (int j=0; j<NXf; ++j){
      for (int k=0; k<N; ++k){
	hessian.elem(i,j) += mf[k] * fpp_x1x2[k].elem(i,j);
      }
    }
  }

  return hessian;
}





// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Calculate approximate Hessian matrix
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename S, typename T, typename U>
Matrix<double> ChiSqFunc<S,T,U> ::Hessian(const Vector<double> & xi){
  if (point_is_new(xi)){
    mParam.Xupdate(xi);
    mstatus_X   = updated;
    mstatus_MDY = expired;
    mstatus_f   = expired;
    mstatus_J   = expired;
    mstatus_V     = expired;
    mstatus_dV_dx = expired;
  }
  return Hessian();
}

template <typename S, typename T, typename U>
Matrix<double> ChiSqFunc<S,T,U> ::Hessian(void){
  int NXf = mParam.NXfree();
  Matrix<double> hessian(NXf, NXf, 0);

  mJ = J();
  hessian = mJ.transpose() * mJ;

  double eps = std::numeric_limits<double>::epsilon();
  double tbs = mbarrier_scale;
  if (tbs < 0.0) tbs *= -1.0; 
  if (tbs < eps)
    return hessian;


  // Penalty/Barrier function contribution:
  funcfit::bad_point ebp;
  double xi, ai, bi, td1, td2, d2U_dx2;
  int n=0;
  for (int i=0; i< mParam.X().size(); ++i){
    if (mParam.X_is_unconstrained(i)) continue;
    if (mParam.X_is_fixed(i)) continue;

    xi = mParam.X(i);
    ai = mParam.Xmin(i);
    bi = mParam.Xmax(i);
    if (xi<=ai || xi>=bi){
      std::cout << "chisq-basics: Hessian_exact: Warning: Parameter value " << xi
		<< " is outside limits (" << ai << ", " << bi << ")."
		<< " Throwing exception." << std::endl;
      throw ebp;
    }
    td1 = bi-xi;
    td1 *= td1;
    td2 = xi-ai;
    td2 *= td2;
    d2U_dx2 = mbarrier_scale * (1.0/td1 + 1.0/td2);
    hessian.elem(n,n) += d2U_dx2;
    n++;
  }

  return hessian;
}













// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Get ChiSq value
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename S, typename T, typename U>
double ChiSqFunc<S,T,U> ::value(const Vector<double> & xi){
  if (point_is_new(xi)){
    mParam.Xupdate(xi);
    mstatus_X   = updated;
    mstatus_MDY = expired;
    mstatus_f   = expired;
    mstatus_J   = expired;
    mstatus_V     = expired;
    mstatus_dV_dx = expired;
  }
  return value();
}

template <typename S, typename T, typename U>
double ChiSqFunc<S,T,U> ::value(void){
  return value_data() + value_barrier();
}


template <typename S, typename T, typename U>
double ChiSqFunc<S,T,U> ::value_data(const Vector<double> & xi){
  if (point_is_new(xi)){
    mParam.Xupdate(xi);
    mstatus_X   = updated;
    mstatus_MDY = expired;
    mstatus_f   = expired;
    mstatus_J   = expired;
    mstatus_V     = expired;
    mstatus_dV_dx = expired;
  }
  return value_data();
}

template <typename S, typename T, typename U>
double ChiSqFunc<S,T,U> ::value_data(void){
  double chisq_M=0;

  mf = f();
  for (int i=0; i<mf.size(); ++i) chisq_M += mf[i] * mf[i];
  chisq_M *= 0.5;

  return chisq_M;
}

template <typename S, typename T, typename U>
double ChiSqFunc<S,T,U> ::value_barrier_raw(const Vector<double> & xi){
  if (point_is_new(xi)){
    mParam.Xupdate(xi);
    mstatus_X   = updated;
    mstatus_MDY = expired;
    mstatus_f   = expired;
    mstatus_J   = expired;
    mstatus_V     = expired;
    mstatus_dV_dx = expired;
  }
  return value_barrier_raw();
}

template <typename S, typename T, typename U>
double ChiSqFunc<S,T,U> ::value_barrier_raw(void){

  double eps = std::numeric_limits<double>::epsilon();
  double tbs = barrier_scale();
  //std::cout << "barrier_scale is " << tbs << std::endl;
  if (tbs < 0.0) tbs *= -1.0;
  if (tbs < eps) return 0;
  
  double chisq_U=0;
  mV = V();
  //std::cout << "size of V is " << mV.size() << std::endl;
  for (int i=0; i<mV.size(); ++i) chisq_U -= tbs * mV[i];
  return chisq_U;
}

template <typename S, typename T, typename U>
double ChiSqFunc<S,T,U> ::adjust_barrier(const Vector<double> & xi){
  if (point_is_new(xi)){
    mParam.Xupdate(xi);
    mstatus_X   = updated;
    mstatus_MDY = expired;
    mstatus_f   = expired;
    mstatus_J   = expired;
    mstatus_V     = expired;
    mstatus_dV_dx = expired;
  }
  return adjust_barrier();
}

template <typename S, typename T, typename U>
double ChiSqFunc<S,T,U> ::adjust_barrier(void){
  double chisq_M = value_data();
  double chisq_U = value_barrier_raw();

  double eps = std::numeric_limits<double>::epsilon();
  double tbs = barrier_scale();
  if (tbs < 0.0) tbs *= -1.0; 
  if (tbs < eps) tbs = 0;

  double chisq_U_old = chisq_U;
  double f=1.0;
  if (use_barrier_rescaling()){
    if (chisq_U_old > eps) f = chisq_M / chisq_U_old;
    barrier_scale() *= f;
    chisq_U_old *= f;
  }
  return chisq_U_old;
}

template <typename S, typename T, typename U>
double ChiSqFunc<S,T,U> ::value_barrier(const Vector<double> & xi){
  if (point_is_new(xi)){
    mParam.Xupdate(xi);
    mstatus_X   = updated;
    mstatus_MDY = expired;
    mstatus_f   = expired;
    mstatus_J   = expired;
    mstatus_V     = expired;
    mstatus_dV_dx = expired;
  }
  return value_barrier();
}

template <typename S, typename T, typename U>
double ChiSqFunc<S,T,U> ::value_barrier(void){
  return adjust_barrier();
}






// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Get gradient of Chisq wrt parameters x
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename S, typename T, typename U>
Vector<double> ChiSqFunc<S,T,U> ::gradient(const Vector<double> & xi){
  if (point_is_new(xi)){
    mParam.Xupdate(xi);
    mstatus_X   = updated;
    mstatus_MDY = expired;
    mstatus_f   = expired;
    mstatus_J   = expired;
    mstatus_V     = expired;
    mstatus_dV_dx = expired;
  }
  return gradient();
}



template <typename S, typename T, typename U>
Vector<double> ChiSqFunc<S,T,U> ::gradient(void){
  mf = f();
  mJ = J();

  double eps = std::numeric_limits<double>::epsilon();
  double tbs = barrier_scale();
  if (tbs < 0.0) tbs *= -1.0; 
  if (tbs < eps)
    return mJ.transpose() * mf;
  else
    return mJ.transpose() * mf + gradient_barrier();
}


template <typename S, typename T, typename U>
Vector<double> ChiSqFunc<S,T,U> ::gradient_data(const Vector<double> & xi){
  if (point_is_new(xi)){
    mParam.Xupdate(xi);
    mstatus_X   = updated;
    mstatus_MDY = expired;
    mstatus_f   = expired;
    mstatus_J   = expired;
    mstatus_V     = expired;
    mstatus_dV_dx = expired;
  }
  return gradient_data();
}

template <typename S, typename T, typename U>
Vector<double> ChiSqFunc<S,T,U> ::gradient_data(void){
  mf = f();
  mJ = J();
  return mJ.transpose() * mf;
}


template <typename S, typename T, typename U>
Vector<double> ChiSqFunc<S,T,U> ::gradient_barrier(const Vector<double> & xi){
  if (point_is_new(xi)){
    mParam.Xupdate(xi);
    mstatus_X   = updated;
    mstatus_MDY = expired;
    mstatus_f   = expired;
    mstatus_J   = expired;
    mstatus_V     = expired;
    mstatus_dV_dx = expired;
  }
  return gradient_barrier();
}

template <typename S, typename T, typename U>
Vector<double> ChiSqFunc<S,T,U> ::gradient_barrier(void){
  mdV_dx = dV_dx();

  double eps = std::numeric_limits<double>::epsilon();
  double tbs = barrier_scale();
  if (tbs < 0.0) tbs *= -1.0; 

  Vector<double> Xt(mParam.NXfree(), 0);
  for (int i=0; i<mdV_dx.size(); ++i) Xt[i] = - tbs * mdV_dx[i];
  return Xt;
}







// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



// Removes DataX, DataY, DataUncertainty, DataWeightY, DataScaleY, ModelDataY
template <typename S, typename T, typename U>
void ChiSqFunc<S,T,U> ::clear_data(void){
  // mbarrier_scale = 1.0;
  // muse_scales = false;
  mDataX.resize(0);
  mDataY.resize(0);
  mDataUncertaintyY.resize(0);
  mModelDataY.resize(0);

  mstatus_X   = updated;
  mstatus_MDY = expired;
  mstatus_f   = expired;
  mstatus_J   = expired;
  mstatus_V     = expired;
  mstatus_dV_dx = expired;
}






// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename S, typename T, typename U>
void ChiSqFunc<S,T,U> ::finalize_setup(){

  // Check if any prefactor vectors are set:
  if (mDataUncertaintyY.size() != mDataY.size() &&
      mDataWeightY.size()      != mDataY.size()){
    aborterror("Neither uncertainties nor weights are in use! Please fill any or both vectors.");
  }

  // 1. Suppose uncertainty is set. Check for negative values. If such values found,
  // then there must be a weight vector of the same size, and with certain properties.
  if (mDataUncertaintyY.size() == mDataY.size()){
    int neg=0;
    for (int i=0; i<mDataUncertaintyY.size(); ++i){
      if (mDataUncertaintyY[i]<0) neg++;
    }
    // No negative values:
    if (neg==0) mDataWeightY = Vector<U>(mDataY.size(), -1);

    // Negative values means we need the weight vector:
    if (neg>0 && mDataWeightY.size() != mDataY.size()){
      // aborterror("Uncertainty vector OK, but weight vector needed, and not supplied (or of wrong dimension).");
      mDataWeightY = Vector<U>(mDataY.size(), 1);
      for (int i=0; i<mDataUncertaintyY.size(); ++i){
	if (mDataUncertaintyY[i]<0) mDataWeightY[i]=U(1);
      }
    }

    // Check if negative uncertainty values are consistent with positive weight values:
    if (neg>0 && mDataWeightY.size() == mDataY.size()){
      int notOK1=0, notOK2=0;
      for (int i=0; i<mDataUncertaintyY.size(); ++i){
	if (mDataUncertaintyY[i]<0 && mDataWeightY[i]<0) notOK1++;
	if (mDataUncertaintyY[i]>0 && mDataWeightY[i]>0) notOK2++;
      }
      if (notOK1>0 || notOK2>0)
	aborterror("Uncertainty and weight vectors OK, but weight vector elements not consistent with uncertainty elements.");
    }
  }


  // 2. Suppose weight is set. Check for negative values. If such values found,
  // then there must be an uncertainty vector of the same size, and with certain properties.
  if (mDataWeightY.size() == mDataY.size()){
    int neg=0;
    for (int i=0; i<mDataWeightY.size(); ++i){
      if (mDataWeightY[i]<0) neg++;
    }
    // No negative values:
    if (neg==0) mDataUncertaintyY = Vector<U>(mDataY.size(), -1);

    // Negative values means we need the uncertainty vector:
    if (neg>0 && mDataUncertaintyY.size() != mDataY.size()){
      // aborterror("Weight vector OK, but uncertainty vector needed, and not supplied (or of wrong dimension).");
      mDataUncertaintyY = Vector<U>(mDataY.size(), 1);
      for (int i=0; i<mDataWeightY.size(); ++i){
	if (mDataWeightY[i]<0) mDataUncertaintyY[i]=U(1);
      }
    }      

    // Check if negative weight values are consistent with positive uncertainty values:
    if (neg>0 && mDataUncertaintyY.size() == mDataY.size()){
      int notOK1=0, notOK2=0;
      for (int i=0; i<mDataWeightY.size(); ++i){
	if (mDataWeightY[i]<0 && mDataUncertaintyY[i]<0) notOK1++;
	if (mDataWeightY[i]>0 && mDataUncertaintyY[i]>0) notOK2++;
      }
      if (notOK1>0 || notOK2>0)
	aborterror("Weight and uncertainty vectors OK, but uncertainty vector elements not consistent with weight elements.");
    }
  }


  // 3. Normalize weights:
  double tmp=0.0;
  for (int i=0; i<mDataWeightY.size(); ++i)
    if (mDataWeightY[i]>0) tmp += double(mDataWeightY[i]);
  tmp = 1.0/tmp;
  for (int i=0; i<mDataWeightY.size(); ++i)
    if (mDataWeightY[i]>0) mDataWeightY[i] = U(tmp * double(mDataWeightY[i]));








  double eps = std::numeric_limits<double>::epsilon();

  mDataScaleY.resize( mDataY.size() );
  for (int i=0; i<mDataY.size(); ++i)
    mDataScaleY[i] = U(1);

  if (muse_scales){
    // Set scales of DataY
    for (int i=0; i<mDataY.size(); ++i){
      mDataScaleY[i] = ((mDataY[i] < 0) ? -1 : 1) * mDataY[i];
      if (double(mDataScaleY[i]) < eps)
	mDataScaleY[i] = U(eps);
    }
  }


}




template <typename S, typename T, typename U>
void ChiSqFunc<S,T,U> ::set_point(const Vector<double> & xi){
  if (point_is_new(xi)){
    mParam.Xupdate(xi);
    mstatus_X   = updated;
    mstatus_MDY = expired;
    mstatus_f   = expired;
    mstatus_J   = expired;
    mstatus_V     = expired;
    mstatus_dV_dx = expired;
  }
}

// Use this after catching an exception. This resets the point to being a new
// point, so things are recalculated to make sure the internal state is good.
template <typename S, typename T, typename U>
void ChiSqFunc<S,T,U> ::reset(const Vector<double> & xi){
  mstatus_X   = expired;
  mstatus_MDY = expired;
  mstatus_f   = expired;
  mstatus_J   = expired;
  mstatus_V     = expired;
  mstatus_dV_dx = expired;

  mParam.Xupdate(xi);
  mstatus_X   = updated;
  mstatus_MDY = expired;
  mstatus_f   = expired;
  mstatus_J   = expired;
  mstatus_V     = expired;
  mstatus_dV_dx = expired;
}


template <typename S, typename T, typename U>
void ChiSqFunc<S,T,U> ::reset(void){
  mstatus_X   = expired;
  mstatus_MDY = expired;
  mstatus_f   = expired;
  mstatus_J   = expired;
  mstatus_V     = expired;
  mstatus_dV_dx = expired;
}


template <typename S, typename T, typename U>
void ChiSqFunc<S,T,U> ::backup(){
  X_bak = mParam.X();
  mModelDataY_bak = mModelDataY;
  mf_bak = mf;
  mJ_bak = mJ;
  mbarrier_scale_bak = mbarrier_scale;
  mstatus_X_bak   = mstatus_X;
  mstatus_MDY_bak = mstatus_MDY;
  mstatus_f_bak   = mstatus_f;
  mstatus_J_bak   = mstatus_J;
  mstatus_V_bak     = mstatus_V;
  mstatus_dV_dx_bak = mstatus_dV_dx;
}

template <typename S, typename T, typename U>
void ChiSqFunc<S,T,U> ::restore(){
  mParam.Xupdate(X_bak);
  mModelDataY = mModelDataY_bak;
  mf = mf_bak;
  mJ = mJ_bak;
  mbarrier_scale = mbarrier_scale_bak;
  mstatus_X   = mstatus_X_bak;
  mstatus_MDY = mstatus_MDY_bak;
  mstatus_f   = mstatus_f_bak;
  mstatus_J   = mstatus_J_bak;
  mstatus_V     = mstatus_V_bak;
  mstatus_dV_dx = mstatus_dV_dx_bak;
}


template <typename S, typename T, typename U>
void ChiSqFunc<S,T,U> ::trial_do(const Vector<double> & xi){
  backup();
  set_point(xi);
}

template <typename S, typename T, typename U>
void ChiSqFunc<S,T,U> ::trial_accept(){
  /* May need to adjust barrier_scale! This is normally only done when
     value() or operator() is used.
   */
  adjust_barrier();
}

template <typename S, typename T, typename U>
void ChiSqFunc<S,T,U> ::trial_reject(){
  restore();
}




template <typename S, typename T, typename U>
void ChiSqFunc<S,T,U> ::debug(){

  cout << "All  parameters: " << mParam.X() << endl;
  //cout << "Free parameters: " << mParam.Xfree() << endl;
  cout << "Number of DataX, DataY points: " << mDataX.size() << "  " << mDataY.size() << endl;
  cout << "DataY            vector: " << mDataY << endl;
  cout << "DataScaleY       vector: " << mDataScaleY << endl;
  cout << "DataWeightY      vector: " << mDataWeightY << endl;
  cout << "DataUncertaintyY vector: " << mDataUncertaintyY << endl;
  cout << "barrier scale          : " << mbarrier_scale << endl;
  cout << "use_barrier rescaling  : " << muse_barrier_rescaling << endl;
  //cout << "V vector               : " << mV << endl;
  //cout << "dV_dx vector           : " << mdV_dx << endl;
  //cout << "ModelFuncPointer : " << mModelFuncPointer << endl;
  //cout << "ReportFuncPointer: " << mReportFuncPointer << endl;

}




#endif


