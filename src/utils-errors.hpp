


#ifndef UTILS_ERRORS_HPP
#define UTILS_ERRORS_HPP


#include <exception>
#include <iostream>
#include <string>

#include <cstdlib>


namespace utils {


  class bad_input : public std::exception
  {
  public:
    bad_input(std::string s); // constructor
  };
  







  void aborterror(const std::string s);



}



#endif

