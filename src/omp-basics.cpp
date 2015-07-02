

#include "omp-basics.hpp"



OMP_Info::OMP_Info(){
  mnt_use = 1;
  mnt_max = omp_get_max_threads();
}

// ************************************************************************
// Works inside parallel region, otherwise returns 1.
int  OMP_Info::nt_use(void){ return (mnt_use = omp_get_num_threads()); }
// ************************************************************************

void OMP_Info::nt_use(int n){
  if (n > mnt_max) n=mnt_max;
  omp_set_num_threads(n);
}

int  OMP_Info::nt_max(void){
  return mnt_max;
  //return (mnt_max = omp_get_max_threads());
}

// ************************************************************************
// Works inside parallel region, otherwise returns 0 (master thread).
int  OMP_Info::tid(void){ return (mtid = omp_get_thread_num()); }
// ************************************************************************






int myomp_get_chunksize(size_t sizeoftype){
#ifdef LEVEL1_DCACHE_LINESIZE
  return (double)LEVEL1_DCACHE_LINESIZE/sizeoftype;
#else
  return 64.0/sizeoftype; // reasonable (?) default
#endif
}


