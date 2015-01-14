


#include <string>
#include <vector>
#include <sstream>

#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>

#include "utils-string.hpp"




/* #########################################################################
   Extract substrings from a string:
   #########################################################################
*/

int utils::get_substrings( const std::string & in,
			    std::vector<std::string> & substr,
			    const std::string & delims ){

  std::string buf(in);

  // Remove leading and trailing whitespace, i.e. use boost::trim:
  boost::trim_if(buf, boost::is_any_of("\t "));

  // Don't bother with boost::split if trimmed input is empty:
  if (buf.size()==0) return 0;
  
  substr.clear();
  substr.resize(0);
  boost::split(substr, buf, boost::is_any_of(delims), boost::token_compress_on);
  return substr.size();
}



bool utils::bool_in_string(const string arg){
  char c = arg[0];

  if (c=='y' || c=='Y' || c=='t' || c=='T' || c=='1') return true;
  if (c!='n' && c!='N' && c!='f' && c!='F' && c!='0') return true;

  return false;
}


