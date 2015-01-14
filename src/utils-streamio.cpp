


#include <iostream>
#include <string>
#include <vector>
#include <sstream>


#include "utils-streamio.hpp"



/* #########################################################################
   Read line from file:
   #########################################################################
*/

int utils::get_line( std::istream & fin,
		     std::string  & line
		     ){
  int ch;
  int nc=0;

  line.resize(0);

  while ( fin ){
    ch = fin.get();

    // EOF does not indicate failure: When EOF reached 'fin' is true,
    // and the character can be processed in this block. However,
    // next read operation will fail (because EOF reached in the
    // previous operation) to get a character.
      
    // Don't store newline in string.
    if (ch=='\n') break;

    line.push_back(ch);
    nc++;
  }

  return nc;
}




#if 0
  // static int icall=0;
  // std::cout << "get_line call number " << icall << std::endl; icall++;

    if (fin.eof()) ioerr.eof=true; // EOF detected
    if (fin.fail() && !ioerr.eof) ioerr.fail=true; // E.g. corrupted stream
    if (!fin) ioerr.any=true; // Any sort of error
    if (ioerr.eof || ioerr.fail || ioerr.any) break;

  // std::cout << "Read line:" << line << std::endl;

  // Return union of error conditions. User still need to check the size
  // of the returned string for anything that was read before the error
  // was encountered.
#endif

