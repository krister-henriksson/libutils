

#ifndef UTILS_STRING_HPP
#define UTILS_STRING_HPP



#include <string>
#include <vector>
#include <sstream>




namespace utils {


  /* #########################################################################
     Extract substrings from a string:
     #########################################################################
  */

  int get_substrings( const std::string & in,
		       std::vector<std::string> & substr,
		       const std::string & delims );

  

  template <typename T>
  inline std::string tostring(const T & input){
    std::ostringstream strbuf;

    strbuf.clear();
    strbuf << input;
    return strbuf.str();
  }
    
}



#endif


