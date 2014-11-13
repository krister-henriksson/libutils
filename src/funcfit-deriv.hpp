

#ifndef FUNCFIT_DERIV_HPP
#define FUNCFIT_DERIV_HPP



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
#include "exiterrors.hpp"

#include "param.hpp"



using utils::Vector;
using std::numeric_limits;





namespace funcfit {

  // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  // Get parameter derivative of a function that returns a scalar:
  // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  template <typename T>
  void func_deriv(T              & func,
		  Vector<double> & x,
		  Vector<double> & deriv);
  

}


// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
// Get parameter derivative of a function that returns a scalar:
// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

template <typename T>
void funcfit ::func_deriv(T              & func,
			  Vector<double> & Xin,
			  Vector<double> & deriv){

  Vector<double> Xf = Xin;//func.Param().Xfree();
  int NXf = Xf.size();
  double eps = std::numeric_limits<double>::epsilon();
  double bakxi, tmpx, dx, DY1, DY2;
  int ix;

  deriv.resize(NXf);

  for (int ix=0; ix<NXf; ++ix){
    bakxi = Xf[ix];
    tmpx  = bakxi < 0 ? -bakxi : bakxi;
    dx = pow(eps, 1.0/3.0);
    if (tmpx>eps) dx = tmpx * dx;

    // Positive perturbation in parameter:
    Xf[ix] = bakxi + dx;
    DY1 = func(Xf);

    // Negative perturbation in parameter:
    Xf[ix] = bakxi - dx;
    DY2 = func(Xf);

    // Reset:
    Xf[ix] = bakxi;

    deriv[ix] = (DY1 - DY2)/(2*dx);
  }

  return;
}




















#endif
