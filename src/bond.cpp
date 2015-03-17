


#include "bond.hpp"



BondData::BondData(){
  dist=0.0;
  nbonds=0;
}



BondAngleData::BondAngleData(){
  for (int i=0; i<8; ++i){
    costheta_ijk[i] = 0.0;
    ntheta_ijk[i] = 0;
    typei[i]="";
    typej[i]="";
    typek[i]="";
  }
}







