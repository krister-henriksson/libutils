


#ifndef NR_F1DIM_HPP
#define NR_F1DIM_HPP


#include <cstdlib>

#include "utils.hpp"
#include "utils-math.hpp"
#include "nr-golden.hpp"




///////////////////////////////////////////////////////////
// Usage:
///////////////////////////////////////////////////////////
//
//   // func is functor or function
//
///////////////////////////////////////////////////////////



namespace nr {

  template <typename T>
  class F1dim
  {

    //#####################################################################
    // Member variables
    //#####################################################################
  public:
    const utils::Vector<double> & p;  // reference, not a local variable
    const utils::Vector<double> & xi; // reference, not a local variable
    int n;
    T & func;
    utils::Vector<double> xt;
    bool debug;


    //#####################################################################
    // Constructor
    //#####################################################################
    F1dim(utils::Vector<double> & pp,
	  utils::Vector<double> & xii,
	  T                     & funcc)
      : p(pp), xi(xii), n(p.size()), func(funcc), xt(n, 0), debug(false)
    {
      /*
      std::cout << "F1dim: size of p  vector: " << p.size() << std::endl;
      std::cout << "F1dim: size of xi vector: " << xi.size() << std::endl;
      std::cout << "F1dim: size of xt vector: " << xt.size() << std::endl;
      */
    }


    //#####################################################################
    // Main method: operator definition
    //#####################################################################
    double operator()(const double x){
      /*
      std::cout << "f1dim: operator(): size of p  vector:" << p.size() << std::endl;
      std::cout << "f1dim: operator(): size of xi vector:" << xi.size() << std::endl;
      std::cout << "f1dim: operator(): size of xt vector:" << xt.size() << std::endl;
      */
      for (int j=0; j!=n; ++j)
	xt[j] = p[j] + x*xi[j];

      //std::cout << "F1dim point now: " << xt << std::endl;
      return func(xt);
    }

  } ;


}

#endif

