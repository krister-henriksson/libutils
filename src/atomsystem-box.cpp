

#include <iostream>
#include <fstream>
#include <sstream>
#include <new>
#include <string>
#include <vector>
#include <limits>
#include <algorithm>

//#include <boost/format.hpp>
//#include <boost/range/algorithm/sort.hpp>


#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

#include "utils.hpp"
#include "utils-math.hpp"
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
//#include "bond.hpp"


using std::cout;
using std::endl;
using std::ofstream;
using std::ios;
using std::ifstream;
using std::istringstream;
using std::string;
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
using utils::fp_are_equal;

using namespace constants;




// Set box direction vectors and length of box in the specified direction. The vectors
// are normalized inside the method.
// NOTE: The box direction vectors are stored column-wise in the mboxdir matrix.
void AtomSystem::set_boxdir(const int idir,
			    const double & px,
			    const double & py,
			    const double & pz
			    ){
  //  if (boxdir.nrows()<3 || boxdir.ncols()<3) boxdir.resize(3,3);
  if (idir<0 || idir>=3){
    cout << "Index " << idir << " is out of range. Allowed: 0 to 2." << endl;
    exit(EXIT_FAILURE);
  }

  Vector3<double> tmpv(0.0);
  tmpv[0]=px;
  tmpv[1]=py;
  tmpv[2]=pz;
  tmpv.normalize();
  boxdir.col(idir, tmpv);
}


void AtomSystem::set_boxdir(const int idir, Vector3<double> & p){

  //  if (boxdir.nrows()<3 || boxdir.ncols()<3) boxdir.resize(3,3);
  if (idir<0 || idir>=3){
    cout << "Index " << idir << " is out of range. Allowed: 0 to 2." << endl;
    exit(EXIT_FAILURE);
  }

  p.normalize();
  boxdir.col(idir, p);
}



// Get box direction vectors:
void AtomSystem::get_boxdir(const int idir,
			    Vector3<double> & v) const {
  /*
  if (boxdir.nrows()<3 || boxdir.ncols()<3){
    cout << "Box has not been set yet. Exiting." << endl;
    exit(EXIT_FAILURE);
  }
  */
  if (idir<0 || idir>=3){
    cout << "Index " << idir << " is out of range. Allowed: 0 to 2." << endl;
    exit(EXIT_FAILURE);
  }
  v = boxdir.col(idir);
}



// Normalize the box direction vectors and build the Bravais matrix and its inverse:
void AtomSystem::update_box_geometry(){
  Vector3<double> tmpv(0);

  // Normalize:
  tmpv = boxdir.col(0);  tmpv.normalize();  boxdir.col(0, tmpv);
  tmpv = boxdir.col(1);  tmpv.normalize();  boxdir.col(1, tmpv);
  tmpv = boxdir.col(2);  tmpv.normalize();  boxdir.col(2, tmpv);


  // Calculate the Bravaismatrix:
  /* The mboxdir matrix --- with the box direction vectors stored as columns ---
     is actually the Bravais matrix. We need to keep it up to date and also
     its inverse.
  */

  
  /*
  rvec = x xhat + y yhat + z zhat
       = A ahat + B bhat + C chat

       x = rvec * xhat = A ahat * xhat + B bhat * xhat + C chat * xhat
       y = rvec * yhat = A ahat * yhat + B bhat * yhat + C chat * yhat
       z = rvec * zhat = A ahat * zhat + B bhat * zhat + C chat * zhat

       <=>

      | x |   | ahat * xhat   bhat * xhat   chat * xhat |   | A |
      | y | = | ahat * yhat   bhat * yhat   chat * yhat | * | B |
      | z |   | ahat * zhat   bhat * zhat   chat * zhat |   | C |

       <=>

       rvec_xyz = BM      * rvec_abc
       rvec_abc = BM^(-1) * rvec_xyz

       To get coordinates in xyz space multiply coordinates in abc space by BM.
       To get coordinates in abc space multiply coordinates in xyz space by BM^(-1).
  */

  isCart = false;
  double td0,td1,td2, eps = numeric_limits<double>::epsilon();;
  int k=0;

  // boxdir(0)
  td0 = boxdir.col(0)[0]; if (td0<0) td0 *= -1;
  td1 = boxdir.col(0)[1]; if (td1<0) td1 *= -1;
  td2 = boxdir.col(0)[2]; if (td2<0) td2 *= -1;
  if (td1<eps && td2<eps) ++k;
  // boxdir(1)
  td0 = boxdir.col(1)[0]; if (td0<0) td0 *= -1;
  td1 = boxdir.col(1)[1]; if (td1<0) td1 *= -1;
  td2 = boxdir.col(1)[2]; if (td2<0) td2 *= -1;
  if (td0<eps && td2<eps) ++k;
  // boxdir(2)
  td0 = boxdir.col(2)[0]; if (td0<0) td0 *= -1;
  td1 = boxdir.col(2)[1]; if (td1<0) td1 *= -1;
  td2 = boxdir.col(2)[2]; if (td2<0) td2 *= -1;
  if (td0<eps && td1<eps) ++k;

  if (k == 3) isCart = true;


  /*
  LUdecomp<double> lud ( boxdir );
  lud.inverse(  Bravaismatrix_inv );
  */
  if (! isCart)
    boxdir.inverse( Bravaismatrix_inv );
  else {
    Bravaismatrix_inv.elem(0,0) = 1.0;
    Bravaismatrix_inv.elem(1,0) = 0.0;
    Bravaismatrix_inv.elem(2,0) = 0.0;

    Bravaismatrix_inv.elem(0,1) = 0.0;
    Bravaismatrix_inv.elem(1,1) = 1.0;
    Bravaismatrix_inv.elem(2,1) = 0.0;

    Bravaismatrix_inv.elem(0,2) = 0.0;
    Bravaismatrix_inv.elem(1,2) = 0.0;
    Bravaismatrix_inv.elem(2,2) = 1.0;
  }


  /*
  cout << "Bravaismatrix:" << endl;
  cout << mboxdir << endl;

  cout << "Bravaismatrix_inverse:" << endl;
  cout << mBravaismatrix_inv << endl;
  */
}



void AtomSystem::calc_volume(void){
  Vector3<double> u = boxdir.col(0);
  Vector3<double> v = boxdir.col(1);
  Vector3<double> w = boxdir.col(2);
  Vector3<double> tv(0);
    
  tv[0] = u[1]*v[2] - u[2]*v[1];
  tv[1] = u[2]*v[0] - u[0]*v[2];
  tv[2] = u[0]*v[1] - u[1]*v[0];
    
  double td = tv[0]*w[0] + tv[1]*w[1] + tv[2]*w[2];
  if (td<0) td *= -1;
    
  vol = td * boxlen[0]*boxlen[1]*boxlen[2];
  vol_atom = vol/natoms();
  return;
}



void AtomSystem::calc_closepacked_volume(){
  double drsqmin=100,drsq=0;
  Vector3<double> posi(0), posj(0), dr(0);
  int i,j,ij,counter=0, nat = natoms();


  vol = 0.0;
  for (i=0; i<nat; ++i){
    posi[0] = pos[i][0];
    posi[1] = pos[i][1];
    posi[2] = pos[i][2];
    
    for (ij=0; ij<neighborcollection[i].size(); ++ij){
      j = neighborcollection[i][ij];
      posj[0] = pos[j][0];
      posj[1] = pos[j][1];
      posj[2] = pos[j][2];



      if (isCart){
	drsq = 0.0;
	dr[0] = posi[0] - posj[0];
	while (pbc[0] && dr[0]<-0.5*boxlen[0]) dr[0] += boxlen[0];
	while (pbc[0] && dr[0]>=0.5*boxlen[0]) dr[0] -= boxlen[0];
	drsq += dr[0] * dr[0];
	dr[1] = posi[1] - posj[1];
	while (pbc[1] && dr[1]<-0.5*boxlen[1]) dr[1] += boxlen[1];
	while (pbc[1] && dr[1]>=0.5*boxlen[1]) dr[1] -= boxlen[1];
	drsq += dr[1] * dr[1];
	dr[2] = posi[2] - posj[2];
	while (pbc[2] && dr[2]<-0.5*boxlen[2]) dr[2] += boxlen[2];
	while (pbc[2] && dr[2]>=0.5*boxlen[2]) dr[2] -= boxlen[2];
	drsq += dr[2] * dr[2];
      }
      else {
	get_atom_distance_vec(posi, posj, dr);
	drsq = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
      }
      if (counter==0 || (counter>0 && drsqmin<drsq)){
	drsqmin = drsq;
	counter++;
      }
    }

    vol += 4*PI*drsqmin * sqrt(drsqmin)/3.0;
  }
  vol_atom = vol/nat;
  dr_min   = sqrt(drsqmin);

  return;
}





void AtomSystem::get_bond_list(Vector<BondData> & bond_list,
			       string & name1,
			       string & name2,
			       int & nat_with_bonds,
			       double rc12
			       ){

  // Consider A-B and B-A bonds, since they both contribute to the
  // A-B pair (bond) energy.

  double drsq=0,r,rcutsq;
  Vector3<double> posi(0), posj(0), dr(0);
  int i,j,ij,nat = natoms(),k;
  BondData bond;

  bond_list.resize(0);

  rcutsq = rc12 * rc12;

  nat_with_bonds = 0;


  for (i=0; i<nat; ++i){
    if (matter[i] != name1 || matter[i] != name2) continue;

    nat_with_bonds++;

    posi[0] = pos[i][0];
    posi[1] = pos[i][1];
    posi[2] = pos[i][2];
   
    for (ij=0; ij<neighborcollection[i].size(); ++ij){
      j = neighborcollection[i][ij];
      if ( (matter[i] == name1 && matter[j] != name2) ||
	   (matter[i] == name2 && matter[j] != name1) ) continue;
      /*
      bool c1 = ((matter[i] == name1) && (matter[j] == name2));
      bool c2 = ((matter[j] == name1) && (matter[i] == name2));
      // If c1 or c2 is true this pair should be investigated.
      if ( ! (c1 || c2) ) continue;
      */
      //if (j<i) continue;
      posj[0] = pos[j][0];
      posj[1] = pos[j][1];
      posj[2] = pos[j][2];

      if (isCart){
	drsq = 0.0;
	dr[0] = posi[0] - posj[0];
	while (pbc[0] && dr[0]<-0.5*boxlen[0]) dr[0] += boxlen[0];
	while (pbc[0] && dr[0]>=0.5*boxlen[0]) dr[0] -= boxlen[0];
	drsq += dr[0] * dr[0];
	dr[1] = posi[1] - posj[1];
	while (pbc[1] && dr[1]<-0.5*boxlen[1]) dr[1] += boxlen[1];
	while (pbc[1] && dr[1]>=0.5*boxlen[1]) dr[1] -= boxlen[1];
	drsq += dr[1] * dr[1];
	dr[2] = posi[2] - posj[2];
	while (pbc[2] && dr[2]<-0.5*boxlen[2]) dr[2] += boxlen[2];
	while (pbc[2] && dr[2]>=0.5*boxlen[2]) dr[2] -= boxlen[2];
	drsq += dr[2] * dr[2];
      }
      else {
	get_atom_distance_vec(posi, posj, dr);
	drsq = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
      }
      if (drsq>=rcutsq) continue;

      r = sqrt(drsq);

      bool is_present=false;
      for (k=0; k<bond_list.size(); ++k){
	if (fp_are_equal(bond_list[k].dist, r, 1e-3)){
	  ++bond_list[k].nbonds;
	  is_present=true;
	  break;
	}
      }
      if (!is_present){
	bond.dist   = r;
	bond.nbonds = 1;
	bond_list.push_back(bond);
      }
    }
  }

  bool swapped;

  while (true){
    swapped = false;
    for (k=0; k<bond_list.size()-1; ++k){
      if (bond_list[k].dist > bond_list[k+1].dist){
	BondData tbd = bond_list[k];
	bond_list[k]   = bond_list[k+1];
	bond_list[k+1] = tbd;
	swapped = true; break;
      }
    }
    if (!swapped) break;
  }


  /*
  cout << "List of bond distances 12: " << endl;
  for (k=0; k<bond12_list.size(); ++k) cout << bond12_list[k].dist << endl;
  cout << "List of bond distances 21: " << endl;
  for (k=0; k<bond21_list.size(); ++k) cout << bond21_list[k].dist << endl;
  cout << "List of number of bonds 12: " << endl;
  for (k=0; k<bond12_list.size(); ++k) cout << bond12_list[k].nbonds << endl;
  cout << "List of number of bonds 21: " << endl;
  for (k=0; k<bond21_list.size(); ++k) cout << bond21_list[k].nbonds << endl;
  */

  /* Error: call to swap() from somewhere inside std::sort()
     is ambigous. Ambiguity is introduced via
     using std::ios;
     using std::iostream;
     using std::ostream;
     and possibly other using directives which dumps local swap()
     names into the global namespace !!!
  */

  /*

  std::vector<BondData> vecsort(bond12_list.size());
  for (i=0; i<bond12_list.size(); ++i)
    vecsort[i] = bond12_list[i];
  std::sort( vecsort.begin(), vecsort.end(), compare_bonddata );

  */




}





void AtomSystem::get_bond_angle_list(Vector<BondAngleData> & bondangle_list,
				     string & name1,
				     string & name2,
				     double rc11,
				     double rc22,
				     double rc12
				     ){
  double drsq=0,rcsq;
  Vector3<double> posi(0), posj(0), posk(0), dr(0), drij(0), drik(0);
  int i,j,k,ij,ik,p,nat = natoms(), it;
  BondAngleData bondangle;
  double costheta_ijk;

  bondangle_list.resize(0);

  for (i=0; i<nat; ++i){
    //if (matter[i] != name1 || matter[i] != name2) continue;

    posi[0] = pos[i][0];
    posi[1] = pos[i][1];
    posi[2] = pos[i][2];
   
    for (ij=0; ij<neighborcollection[i].size(); ++ij){
      j = neighborcollection[i][ij];
      if (i==j) continue;

      rcsq = 0.0;
      if      (matter[i] == name1 && matter[j] == name1) rcsq = rc11*rc11;
      else if (matter[i] == name1 && matter[j] == name2) rcsq = rc12*rc12;
      else if (matter[i] == name2 && matter[j] == name1) rcsq = rc12*rc12;
      else if (matter[i] == name2 && matter[j] == name2) rcsq = rc22*rc22;

      /*
      if ( (matter[i] == name1 && matter[j] != name2) ||
	   (matter[i] == name2 && matter[j] != name1) ) continue;
      */

      posj[0] = pos[j][0];
      posj[1] = pos[j][1];
      posj[2] = pos[j][2];
      if (isCart){
	drsq = 0.0;
	dr[0] = posi[0] - posj[0];
	while (pbc[0] && dr[0]<-0.5*boxlen[0]) dr[0] += boxlen[0];
	while (pbc[0] && dr[0]>=0.5*boxlen[0]) dr[0] -= boxlen[0];
	drsq += dr[0] * dr[0];
	dr[1] = posi[1] - posj[1];
	while (pbc[1] && dr[1]<-0.5*boxlen[1]) dr[1] += boxlen[1];
	while (pbc[1] && dr[1]>=0.5*boxlen[1]) dr[1] -= boxlen[1];
	drsq += dr[1] * dr[1];
	dr[2] = posi[2] - posj[2];
	while (pbc[2] && dr[2]<-0.5*boxlen[2]) dr[2] += boxlen[2];
	while (pbc[2] && dr[2]>=0.5*boxlen[2]) dr[2] -= boxlen[2];
	drsq += dr[2] * dr[2];
      }
      else {
	get_atom_distance_vec(posi, posj, dr);
	drsq = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
      }
      if (drsq>=rcsq) continue;
      drij = dr;



      for (ik=0; ik<neighborcollection[i].size(); ++ik){
	k = neighborcollection[i][ik];
	if (i==k) continue;
	if (j==k) continue;

	rcsq = 0.0;
	if      (matter[i] == name1 && matter[k] == name1) rcsq = rc11*rc11;
	else if (matter[i] == name1 && matter[k] == name2) rcsq = rc12*rc12;
	else if (matter[i] == name2 && matter[k] == name1) rcsq = rc12*rc12;
	else if (matter[i] == name2 && matter[k] == name2) rcsq = rc22*rc22;

	//if (matter[k] != name1 || matter[k] != name2) continue;

	posk[0] = pos[k][0];
	posk[1] = pos[k][1];
	posk[2] = pos[k][2];
	if (isCart){
	  drsq = 0.0;
	  dr[0] = posi[0] - posk[0];
	  while (pbc[0] && dr[0]<-0.5*boxlen[0]) dr[0] += boxlen[0];
	  while (pbc[0] && dr[0]>=0.5*boxlen[0]) dr[0] -= boxlen[0];
	  drsq += dr[0] * dr[0];
	  dr[1] = posi[1] - posk[1];
	  while (pbc[1] && dr[1]<-0.5*boxlen[1]) dr[1] += boxlen[1];
	  while (pbc[1] && dr[1]>=0.5*boxlen[1]) dr[1] -= boxlen[1];
	  drsq += dr[1] * dr[1];
	  dr[2] = posi[2] - posk[2];
	  while (pbc[2] && dr[2]<-0.5*boxlen[2]) dr[2] += boxlen[2];
	  while (pbc[2] && dr[2]>=0.5*boxlen[2]) dr[2] -= boxlen[2];
	  drsq += dr[2] * dr[2];
	}
	else {
	  get_atom_distance_vec(posi, posk, dr);
	  drsq = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
	}
	if (drsq>=rcsq) continue;
	drik = dr;



	costheta_ijk = drij * drik / (drij.magn() * drik.magn());

	it = 0;
	if      (matter[i]==name1 && matter[j]==name1 && matter[k]==name1) it=0;
	else if (matter[i]==name1 && matter[j]==name1 && matter[k]==name2) it=1;
	else if (matter[i]==name1 && matter[j]==name2 && matter[k]==name1) it=2;
	else if (matter[i]==name1 && matter[j]==name2 && matter[k]==name2) it=3;
	else if (matter[i]==name2 && matter[j]==name1 && matter[k]==name1) it=4;
	else if (matter[i]==name2 && matter[j]==name1 && matter[k]==name2) it=5;
	else if (matter[i]==name2 && matter[j]==name2 && matter[k]==name1) it=6;
	else if (matter[i]==name2 && matter[j]==name2 && matter[k]==name2) it=7;

	bool is_present=false;
	for (p=0; p<bondangle_list.size(); ++p){
	  if (fp_are_equal(bondangle_list[p].costheta_ijk[it], costheta_ijk, 1e-3)){
	    ++bondangle_list[p].ntheta_ijk[it];
	    is_present=true;
	    break;
	  }
	}
	if (!is_present){
	  bondangle.costheta_ijk[it] = costheta_ijk;
	  bondangle.ntheta_ijk[it]   = 1;
	  bondangle.typei[it] = matter[i];
	  bondangle.typej[it] = matter[j];
	  bondangle.typek[it] = matter[k];
	  bondangle_list.push_back(bondangle);
	}

      }
    }
  }


  for (it=0; it<8; ++it){
    if (name1 == name2 && it>=1) break;

    bool swapped;
    while (true){
      swapped = false;
      for (k=0; k<bondangle_list.size()-1; ++k){
	if (bondangle_list[k].ntheta_ijk[it] < bondangle_list[k+1].ntheta_ijk[it]){
	  BondAngleData tbd = bondangle_list[k];
	  bondangle_list[k]   = bondangle_list[k+1];
	  bondangle_list[k+1] = tbd;
	  swapped = true; break;
	}
      }
      if (!swapped) break;
    }
  }




}









