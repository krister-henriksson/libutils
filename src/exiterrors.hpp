


#ifndef EXIT_ERRORS_HPP
#define EXIT_ERRORS_HPP


#include <exception>
#include <iostream>
#include <string>

#include <cstdlib>


namespace exiterrors {


  class bad_input
    :
    public std::exception
  {
  public:
    bad_input(std::string s); // constructor
  };
  

  void aborterror(const std::string s);



} // namespace



#endif

