

#ifndef UTILS_STREAMIO_HPP
#define UTILS_STREAMIO_HPP

#include <iostream>
#include <string>
#include <vector>
#include <sstream>



namespace utils {


  class StreamIOerror {
  public:
    bool eof;
    bool fail;
    bool any;

    StreamIOerror(){ eof=fail=any=false; }
    operator bool(){ return (eof || fail || any); }

  } ;


  /* #########################################################################
     Read line from file:
     #########################################################################
  */
  
  int get_line( std::istream & fin,
		 std::string  & line
		 );




}


#endif



