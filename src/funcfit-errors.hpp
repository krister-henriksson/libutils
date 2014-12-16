


#ifndef FUNCFIT_ERRORS_HPP
#define FUNCFIT_ERRORS_HPP


#include <exception>

#include <cstdlib>




namespace funcfit {


  class bad_point : public std::exception
  {
    // Point is outside allowed limits or point produces a bad
    // merit function value. E.g. infinites and NaNs occur.
  };


  class stationary_point : public std::exception
  {
  };



}




#endif

