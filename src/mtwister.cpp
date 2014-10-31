

#include <ctime>
#include <cmath>

#include "mtwister.hpp"



// Constructors:
rand_mtwister::rand_mtwister(){
  mseed = static_cast<unsigned int>(std::time(0));
  mgenerator.seed(mseed);

  mp_uni_dist  = 0;
  mp_rand_unif = 0;
  mp_uni_dist  = new DIST(0,1);
  mp_rand_unif = new GEN(mgenerator, *mp_uni_dist);

  mgauss_ave = 0; mgauss_dev = 1;
  mgauss_r1 = mgauss_r2 = 0;
  mgauss_empty = true;
}

rand_mtwister::rand_mtwister(int s){
  mseed = static_cast<unsigned int>(s);
  mgenerator.seed(mseed);

  mp_uni_dist  = 0;
  mp_rand_unif = 0;
  mp_uni_dist  = new DIST(0,1);
  mp_rand_unif = new GEN(mgenerator, *mp_uni_dist);

  mgauss_ave = 0; mgauss_dev = 1;
  mgauss_r1 = mgauss_r2 = 0;
  mgauss_empty = true;
}


// Destructor
rand_mtwister::~rand_mtwister(){
  if (mp_uni_dist)  delete mp_uni_dist;
  if (mp_rand_unif) delete mp_rand_unif;
  mp_uni_dist = 0;
  mp_rand_unif = 0;

  mgauss_ave = 0; mgauss_dev = 1;
  mgauss_r1 = mgauss_r2 = 0;
  mgauss_empty = true;
}




int rand_mtwister::seed(){ return mseed; }

void rand_mtwister::seed(int s){
  mseed = static_cast<unsigned int>(s);
  mgenerator.seed(mseed);

  if (mp_rand_unif) delete mp_rand_unif;
  mp_rand_unif = new GEN(mgenerator, *mp_uni_dist);
}


double rand_mtwister::unif(){ return (*mp_rand_unif)(); }


double rand_mtwister::gauss(double ave, double dev){
  if (ave != mgauss_ave || dev != mgauss_dev){
    mgauss_r1 = mgauss_r2 = 0; mgauss_empty = true;
    mgauss_ave = ave;
    mgauss_dev = dev;
  }

  if (mgauss_empty == false)
    return mgauss_r2;
  else {
    double u1,u2,r,w;
    while (true){
      u1 = 2.0*(*mp_rand_unif)() - 1.0;
      u2 = 2.0*(*mp_rand_unif)() - 1.0;
      w = u1*u1 + u2*u2;
      if (w <= 1.0) break;
    }
    r = sqrt(-2.0 * log(w));
    mgauss_r1 = (r * u1 / sqrt(w)) * dev + ave;
    mgauss_r2 = (r * u2 / sqrt(w)) * dev + ave;
    mgauss_empty = false;
    return mgauss_r1;
  }
}


