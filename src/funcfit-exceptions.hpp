


#ifndef FUNCFIT_EXCEPTIONS_HPP
#define FUNCFIT_EXCEPTIONS_HPP


#include <exception>

#include <cstdlib>




namespace funcfit {


  class bad_point
    :
    public std::exception
  {
  };

  class bad_value
    :
    public std::exception
  {
  };

  class stationary_point
    :
    public std::exception
  {
  };



}




#endif

