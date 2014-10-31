
#ifndef OMP_BASICS_HPP
#define OMP_BASICS_HPP


#include <omp.h>


class OMP_Info {

private:
  int mnt_max;
  int mnt_use;
  int mtid;

public:
  OMP_Info();

  int  nt_max(void);
  int  nt_use(void);
  void nt_use(int n);
  int  tid(void);

} ;




#endif

