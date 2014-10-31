


#include <exception>
#include <iostream>
#include <string>


#include "exiterrors.hpp"


exiterrors::bad_input::bad_input(std::string s){
  std::cout << "Unexpected input: Did not expect '" << s << "'"
	    << "Aborting." << std::endl;
}


void exiterrors::aborterror(const std::string s){
  std::cout << s << std::endl;
  exit(EXIT_FAILURE);
}


