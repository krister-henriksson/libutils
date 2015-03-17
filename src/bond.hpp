


#ifndef BOND_HPP
#define BOND_HPP

#include <string>



class BondData {
public:
  double dist;
  int nbonds;

  BondData();
} ;




// ijk: AAA, AAB, ABA, ABB, BAA, BAB, BBA, BBB
class BondAngleData {
public:
  double costheta_ijk[8];
  int ntheta_ijk[8];
  std::string typei[8];
  std::string typej[8];
  std::string typek[8];

  BondAngleData();
} ;




#endif

