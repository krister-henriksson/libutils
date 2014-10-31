


#include <iostream>
#include "utils-matrix-Choleskydecomp.hpp"

using namespace std;
using namespace utils;


int main(){
  int N=3;
  Matrix<double> M, Lin(N,N,0), L;
  int i,j,k;
  int status;

  for (i=0; i<N; ++i)
    for (j=0; j<=i; ++j)
      Lin.elem(i,j)=1+i+10*j;

  M = Lin * Lin.transpose();

  Choleskydecomp<double> Cdecomp_M(M, status);
  L = Cdecomp_M.L();

  cout << "M = Lin * Lin^T is" << endl;
  cout << M << endl;
  cout << "where Lin is" << endl;
  cout << Lin << endl;
  cout << "Cholesky decomp L is" << endl;
  cout << L << endl;
  cout << "Now L L^T is" << endl;
  cout << L * L.transpose() << endl;



  return EXIT_SUCCESS;
}


