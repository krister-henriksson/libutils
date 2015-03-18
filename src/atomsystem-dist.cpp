

#include <iostream>
#include <fstream>
#include <sstream>
#include <new>
#include <string>
#include <vector>
#include <limits>
#include <algorithm>

#include <boost/format.hpp>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

#include "utils.hpp"
#include "utils-vector.hpp"
#include "utils-vector3.hpp"
#include "utils-matrix.hpp"
#include "utils-matrixsq3.hpp"
#include "utils-matrix-LUdecomp.hpp"
#include "utils-string.hpp"
#include "utils-streamio.hpp"

#include "atomsystem.hpp"
#include "constants.hpp"
#include "utils-errors.hpp"


using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::ofstream;
using std::ifstream;
using std::istringstream;
using std::ios;
using std::numeric_limits;

using utils::Vector;
using utils::Vector3;
using utils::Matrix;
using utils::MatrixSq3;
using utils::LUdecomp;
using utils::get_line;
using utils::get_substrings;
using utils::tostring;
using utils::tostring_fmt;
using utils::aborterror;






// Calculate Cartesian distance vector. Obey periodics, with the cell
// extending between -0.5 * box and 0.5 * box by default. If lowlim>0
// then a distance of 1.0 * box is added to the calculated values to make
// them be between 0 and box.
void AtomSystem::get_atom_distance_vec(const Vector3<double> & r1,
				       const Vector3<double> & r2,
				       Vector3<double> & v,
				       const double lowlim) const {
  double d0,d1,d2;
 
  v[0] = r1[0] - r2[0];
  v[1] = r1[1] - r2[1];
  v[2] = r1[2] - r2[2];

  // Get distance in skew coordinate system, where periodics can be checked:
  d0 = v[0];
  d1 = v[1];
  d2 = v[2];
  if (! isCart){
    d0 = Bravaismatrix_inv.elem(0,0) * v[0]
      + Bravaismatrix_inv.elem(0,1) * v[1]
      + Bravaismatrix_inv.elem(0,2) * v[2];
    d1 = Bravaismatrix_inv.elem(1,0) * v[0]
      + Bravaismatrix_inv.elem(1,1) * v[1]
      + Bravaismatrix_inv.elem(1,2) * v[2];
    d2 = Bravaismatrix_inv.elem(2,0) * v[0]
      + Bravaismatrix_inv.elem(2,1) * v[1]
      + Bravaismatrix_inv.elem(2,2) * v[2];
  }

  // Periodics check:
  double pf1=-0.5, pf2=0.5;
  if (lowlim>0){ pf1 = 0.0; pf2 = 1.0; }

  while (pbc[0] && d0 <  pf1 * boxlen[0]) d0 += boxlen[0];
  while (pbc[0] && d0 >= pf2 * boxlen[0]) d0 -= boxlen[0];
  while (pbc[1] && d1 <  pf1 * boxlen[1]) d1 += boxlen[1];
  while (pbc[1] && d1 >= pf2 * boxlen[1]) d1 -= boxlen[1];
  while (pbc[2] && d2 <  pf1 * boxlen[2]) d2 += boxlen[2];
  while (pbc[2] && d2 >= pf2 * boxlen[2]) d2 -= boxlen[2];



    
  // Get distance in Cartesian coordinate system:
  v[0] = d0;
  v[1] = d1;
  v[2] = d2;
  if (! isCart){
    v[0] = boxdir.elem(0,0) * d0 + boxdir.elem(0,1) * d1 + boxdir.elem(0,2) * d2;
    v[1] = boxdir.elem(1,0) * d0 + boxdir.elem(1,1) * d1 + boxdir.elem(1,2) * d2;
    v[2] = boxdir.elem(2,0) * d0 + boxdir.elem(2,1) * d1 + boxdir.elem(2,2) * d2;
  }
}


// ##############################################################################################



// Get the skew coordinates when Cartesian coordinates are known.
void AtomSystem::get_coords_cart2skew(const Vector3<double> & drc,
				      Vector3<double> & v,
				      const double lowlim) const {

  // Get distance in skew coordinate system, where periodics can be checked:
  v[0] = drc[0];
  v[1] = drc[1];
  v[2] = drc[2];
  if (! isCart){
    v[0] = Bravaismatrix_inv.elem(0,0) * drc[0]
      + Bravaismatrix_inv.elem(0,1) * drc[1]
      + Bravaismatrix_inv.elem(0,2) * drc[2];
    v[1] = Bravaismatrix_inv.elem(1,0) * drc[0]
      + Bravaismatrix_inv.elem(1,1) * drc[1]
      + Bravaismatrix_inv.elem(1,2) * drc[2];
    v[2] = Bravaismatrix_inv.elem(2,0) * drc[0]
      + Bravaismatrix_inv.elem(2,1) * drc[1]
      + Bravaismatrix_inv.elem(2,2) * drc[2];
  }

  // Periodics check:
  double pf1=-0.5, pf2=0.5;
  if (lowlim>0){ pf1 = 0.0; pf2 = 1.0; }
  while (pbc[0] && v[0] <  pf1 * boxlen[0]) v[0] += boxlen[0];
  while (pbc[0] && v[0] >= pf2 * boxlen[0]) v[0] -= boxlen[0];
  while (pbc[1] && v[1] <  pf1 * boxlen[1]) v[1] += boxlen[1];
  while (pbc[1] && v[1] >= pf2 * boxlen[1]) v[1] -= boxlen[1];
  while (pbc[2] && v[2] <  pf1 * boxlen[2]) v[2] += boxlen[2];
  while (pbc[2] && v[2] >= pf2 * boxlen[2]) v[2] -= boxlen[2];
}


// ##############################################################################################


// Get the Cartesian coordinates when skew coordinates are known.
void AtomSystem::get_coords_skew2cart(Vector3<double> & drs,
				      Vector3<double> & v,
				      const double lowlim) const {
  // Periodics check:
  double pf1=-0.5, pf2=0.5;
  if (lowlim>0){ pf1 = 0.0; pf2 = 1.0; }

  while (pbc[0] && drs[0] <  pf1 * boxlen[0]) drs[0] += boxlen[0];
  while (pbc[0] && drs[0] >= pf2 * boxlen[0]) drs[0] -= boxlen[0];
  while (pbc[1] && drs[1] <  pf1 * boxlen[1]) drs[1] += boxlen[1];
  while (pbc[1] && drs[1] >= pf2 * boxlen[1]) drs[1] -= boxlen[1];
  while (pbc[2] && drs[2] <  pf1 * boxlen[2]) drs[2] += boxlen[2];
  while (pbc[2] && drs[2] >= pf2 * boxlen[2]) drs[2] -= boxlen[2];

  // Get distance in Cartesian coordinate system:
  v[0] = drs[0];
  v[1] = drs[1];
  v[2] = drs[2];
  if (! isCart){
    v[0] = boxdir.elem(0,0) * drs[0] + boxdir.elem(0,1) * drs[1] + boxdir.elem(0,2) * drs[2];
    v[1] = boxdir.elem(1,0) * drs[0] + boxdir.elem(1,1) * drs[1] + boxdir.elem(1,2) * drs[2];
    v[2] = boxdir.elem(2,0) * drs[0] + boxdir.elem(2,1) * drs[1] + boxdir.elem(2,2) * drs[2];
  }
}



// ##############################################################################################


/*
  rxyz = (x,y,z) becomes transported by amount
  - f1 in direction a = boxdir.col(0)
  - f2 in direction b = boxdir.col(0)
  - f3 in direction c = boxdir.col(0)
*/
void AtomSystem::translate_cartpos_in_skewspace(Vector3<double> & pos,
						double f1,
						double f2,
						double f3){
  Vector3<double> d;

  // Get coordinates in skew coordinate system:
  d = pos;
  if (! isCart)
    d = Bravaismatrix_inv * pos;

  d[0] += f1;
  d[1] += f2;
  d[2] += f3;

  // Get coordinates in Cartesian coordinate system:
  pos = d;
  if (! isCart)
    pos = boxdir * d;
}



void AtomSystem::handle_pbc_of_positions(const double lowlim){
  int i, nat = natoms();
  double pf1=-0.5, pf2=0.5;
  if (lowlim>0){
    pf1 = 0.0;
    pf2 = 1.0;
  }


#pragma omp parallel for schedule(static)
  for (i=0; i<nat; ++i){

    double drc0, drc1, drc2;
    double drs0, drs1, drs2;

    drc0 = pos[i][0];
    drc1 = pos[i][1];
    drc2 = pos[i][2];

    // Get distance in skew coordinate system, where periodics can be checked:
    drs0 = drc0;
    drs1 = drc1;
    drs2 = drc2;
    if (! isCart){
      drs0 = 0.0
	+ Bravaismatrix_inv.elem(0,0) * drc0
	+ Bravaismatrix_inv.elem(0,1) * drc1
	+ Bravaismatrix_inv.elem(0,2) * drc2;
      drs1 = 0.0
	+ Bravaismatrix_inv.elem(1,0) * drc0
	+ Bravaismatrix_inv.elem(1,1) * drc1
	+ Bravaismatrix_inv.elem(1,2) * drc2;
      drs2 = 0.0
	+ Bravaismatrix_inv.elem(2,0) * drc0
	+ Bravaismatrix_inv.elem(2,1) * drc1
	+ Bravaismatrix_inv.elem(2,2) * drc2;
    }
      
    // Periodics check:
    while (pbc[0] && drs0 <  pf1 * boxlen[0]) drs0 += boxlen[0];
    while (pbc[0] && drs0 >= pf2 * boxlen[0]) drs0 -= boxlen[0];
    while (pbc[1] && drs1 <  pf1 * boxlen[1]) drs1 += boxlen[1];
    while (pbc[1] && drs1 >= pf2 * boxlen[1]) drs1 -= boxlen[1];
    while (pbc[2] && drs2 <  pf1 * boxlen[2]) drs2 += boxlen[2];
    while (pbc[2] && drs2 >= pf2 * boxlen[2]) drs2 -= boxlen[2];

    // Get distance in Cartesian coordinate system:
    drc0 = drs0;
    drc1 = drs1;
    drc2 = drs2;
    if (! isCart){
      drc0 = boxdir.elem(0,0) * drs0 + boxdir.elem(0,1) * drs1 + boxdir.elem(0,2) * drs2;
      drc1 = boxdir.elem(1,0) * drs0 + boxdir.elem(1,1) * drs1 + boxdir.elem(1,2) * drs2;
      drc2 = boxdir.elem(2,0) * drs0 + boxdir.elem(2,1) * drs1 + boxdir.elem(2,2) * drs2;
    }
    pos[i][0] = drc0;
    pos[i][1] = drc1;
    pos[i][2] = drc2;
      
  }

  return;
}


