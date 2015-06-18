


#include "bond.hpp"



BondData::BondData(){
  dist.resize(0);
  ndist.resize(0);
}



BondAngleData::BondAngleData(){
  costheta_ijk.resize(8);
  for (int i=0; i<8; ++i)  costheta_ijk[i].resize(0);
  ncostheta_ijk.resize(8);
  for (int i=0; i<8; ++i)  ncostheta_ijk[i].resize(0);
}


