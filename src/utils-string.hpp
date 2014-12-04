

#ifndef UTILS_STRING_HPP
#define UTILS_STRING_HPP



#include <string>
#include <vector>
#include <sstream>

#include <boost/format.hpp>


using boost::format;
using std::string;
using std::ostringstream;


namespace utils {


  /* #########################################################################
     Extract substrings from a string:
     #########################################################################
  */

  int get_substrings( const std::string & in,
		       std::vector<std::string> & substr,
		       const std::string & delims );

  

  template <typename T>
  string tostring(const T & input){
    ostringstream strbuf;

    strbuf.clear();
    strbuf << input;
    return strbuf.str();
  }



  template <typename T>
  string tostring_fmt(const string & fmt, const T & input){
    ostringstream strbuf;

    strbuf.clear();
    strbuf << format(fmt) % input;
    return strbuf.str();
  }


    
}



#endif


