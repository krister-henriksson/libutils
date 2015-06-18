


#ifndef BOND_HPP
#define BOND_HPP

#include <string>

#include "utils-vector.hpp"

using utils::Vector;



class BondData {
public:
  Vector<double> dist;
  Vector<int>    ndist;

  BondData();
} ;




// ijk: AAA, AAB, ABA, ABB, BAA, BAB, BBA, BBB
class BondAngleData {
public:
  Vector< Vector<double> > costheta_ijk;
  Vector< Vector<int> >    ncostheta_ijk;

  std::string typei;
  std::string typej;
  std::string typek;

  BondAngleData();
} ;






#endif

