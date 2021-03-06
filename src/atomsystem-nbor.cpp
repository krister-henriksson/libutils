


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




// Build all neighbor collections:
void AtomSystem:: get_all_neighborcollections(bool debug){

  int ic[3], M[3], icell, j, ix,ix0, iy,iy0, iz,iz0;
  int *head, *list;
  int Ncells, i, nn_tot, nat=natoms();
  double rm, rmsq, drsq, drsq_min=1e10;
  Vector3<double> rnew(0), posi(0), posj(0), dr(0);
  int counter;

  rm = rcut + skint;
  rmsq = rm * rm;
  //cout << "Using rm = " << rm << " for neighborcollection calculation." << endl;

  neighborcollection.resize( nat );

  /* Subcells for neighbor list creation: */
  M[0] = floor( boxlen[0] / rm );
  M[1] = floor( boxlen[1] / rm );
  M[2] = floor( boxlen[2] / rm );

  if (M[0]==0) M[0]++;
  if (M[1]==0) M[1]++;
  if (M[2]==0) M[2]++;
  
  Ncells = M[0] * M[1] * M[2];

  //cout << "Boxsizes " << boxlen[0] << " " << boxlen[1] << " " << boxlen[2] << endl;
  /*
  cout << "The number of cells established in principal box directions are "
       << M[0] << " " << M[1] << " " << M[2] << endl;
  cout << "The total number of cells is " << Ncells << endl;
  */



  /* ###############################################################################
     ###############################################################################

     (1) SMALL CELL CASE

     ###############################################################################     
     ############################################################################### */

  if ( M[0]<=3 || M[1]<=3 || M[2]<=3 ){
    // Box too small. Use the O(N^2) method.
    //cout << "Using O(N^2) method for neighbor list creation." << endl;


    nn_tot = 0;
    counter=0;
    for (i=0; i<nat; ++i){
      neighborcollection[i].resize(0);

      // Cartesian position:
      posi[0] = pos[i][0];
      posi[1] = pos[i][1];
      posi[2] = pos[i][2];

      for (j=0; j<nat; ++j){
	if (i==j) continue;

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



	//cout << posi << posj << drsq << endl;
	if (drsq < rmsq){
	  if (counter==0 || (counter>0 && drsq<drsq_min)) drsq_min = drsq;

	  /*
	  cout << "Atom index " << i << ": Added neighbor atom index " 
		    << j << " at position " << nn << endl;
	  */
	  neighborcollection[i].push_back(j);
	  nn_tot++;
	}
      }
      /*
      cout << "Size of neighborcollection of atom " << i << " is: "
	   << neighborcollection[i].size() << endl;
      */
      ++counter;
    }

    dr_min = sqrt(drsq_min);

    if (debug)
      cout << "AtomSystem: linked-list unused: <N_neighb/atom> " << nn_tot / (1.0*nat)
	   << " dr_min " << dr_min << endl;



    return;

  }

  /* ###############################################################################
     ###############################################################################
     (1) END
     ###############################################################################     
     ############################################################################### */




  /* ###############################################################################
     ###############################################################################

     (2) NORMAL CASE

     ###############################################################################     
     ############################################################################### */

    
  /* ###############################################################################
     Prepare neighbor list.
     ############################################################################### */
    

  
    
  /* Arrays needed for linked-list neighbor list creation: */
  head = new int [Ncells];
  list = new int [nat];
  for (i=0; i<Ncells; i++) head[i] = -1;
  for (i=0; i<nat; i++) list[i] = -1;

  
  /* ###############################################################################
     Build linked list.
     ############################################################################### */
  //cout << "Building linked-list of neighbors ..." << endl;
  
  for (i=0; i<nat; i++){
    // Cartesian position:
    posi[0] = pos[i][0];
    posi[1] = pos[i][1];
    posi[2] = pos[i][2];

    // Get skew coordinates:
    if (isCart){
      rnew[0] = posi[0];
      while (pbc[0] && posi[0] <  0.0      ) posi[0] += boxlen[0];
      while (pbc[0] && posi[0] >= boxlen[0]) posi[0] -= boxlen[0];
      rnew[1] = posi[1];
      while (pbc[1] && posi[1] <  0.0      ) posi[1] += boxlen[1];
      while (pbc[1] && posi[1] >= boxlen[1]) posi[1] -= boxlen[1];
      rnew[2] = posi[2];
      while (pbc[2] && posi[2] <  0.0      ) posi[2] += boxlen[2];
      while (pbc[2] && posi[2] >= boxlen[2]) posi[2] -= boxlen[2];
    }
    else
      get_coords_cart2skew(posi, rnew, 1.0);



    // Get index of the skew coordinate cell where the atom is:
    ic[0] = floor( rnew[0]/boxlen[0] * M[0] );
    ic[1] = floor( rnew[1]/boxlen[1] * M[1] );
    ic[2] = floor( rnew[2]/boxlen[2] * M[2] );

    while (pbc[0] == false && ic[0] <  0)     ic[0] = 0;
    while (pbc[0] == false && ic[0] >= M[0])  ic[0] = M[0]-1;
    while (pbc[0] == true  && ic[0] <  0)     ic[0] += M[0];
    while (pbc[0] == true  && ic[0] >= M[0])  ic[0] -= M[0];

    while (pbc[1] == false && ic[1] <  0)     ic[1] = 0;
    while (pbc[1] == false && ic[1] >= M[1])  ic[1] = M[1]-1;
    while (pbc[1] == true  && ic[1] <  0)     ic[1] += M[1];
    while (pbc[1] == true  && ic[1] >= M[1])  ic[1] -= M[1];

    while (pbc[2] == false && ic[2] <  0)     ic[2] = 0;
    while (pbc[2] == false && ic[2] >= M[2])  ic[2] = M[2]-1;
    while (pbc[2] == true  && ic[2] <  0)     ic[2] += M[2];
    while (pbc[2] == true  && ic[2] >= M[2])  ic[2] -= M[2];

    icell = ic[0] + M[0]*( ic[1] + M[1] * (ic[2]) );
    // Index obtained.
	  
    list[i]     = head[icell];
    head[icell] = i;
  }
  //cout << "Linked-list completed." << endl;

  /*
  nocc_tot = 0;
  for (i=0; i<Ncells; ++i){
    cout << "occupancy of cell " << i << " is " << nocc[i] << endl;
    nocc_tot += nocc[i];
  }
  cout << "Number of atoms is " << natoms() << " and number of atoms in cells is " << nocc_tot << endl;
  */    



      
  /* ###############################################################################
     Build Verlet list.
     ############################################################################### */
  //cout << "Building Verlet-list of neighbors ..." << cout;

  nn_tot = 0;
  counter=0;
  for (i=0; i<nat; ++i){
    neighborcollection[i].resize(0);

    // Cartesian position:
    posi[0] = pos[i][0];
    posi[1] = pos[i][1];
    posi[2] = pos[i][2];

    // Get skew coordinates:
    if (isCart){
      rnew[0] = posi[0];
      while (pbc[0] && posi[0] <  0.0      ) posi[0] += boxlen[0];
      while (pbc[0] && posi[0] >= boxlen[0]) posi[0] -= boxlen[0];
      rnew[1] = posi[1];
      while (pbc[1] && posi[1] <  0.0      ) posi[1] += boxlen[1];
      while (pbc[1] && posi[1] >= boxlen[1]) posi[1] -= boxlen[1];
      rnew[2] = posi[2];
      while (pbc[2] && posi[2] <  0.0      ) posi[2] += boxlen[2];
      while (pbc[2] && posi[2] >= boxlen[2]) posi[2] -= boxlen[2];
    }
    else
      get_coords_cart2skew(posi, rnew, 1.0);

    // cout << "Cartesian and skew positions:" << endl;
    // cout << posi << rnew << endl;

    // Get index of the skew coordinate cell where the atom is:
    ic[0] = floor( rnew[0]/boxlen[0] * M[0] );
    ic[1] = floor( rnew[1]/boxlen[1] * M[1] );
    ic[2] = floor( rnew[2]/boxlen[2] * M[2] );

    while (pbc[0] == false && ic[0] <  0)     ic[0] = 0;
    while (pbc[0] == false && ic[0] >= M[0])  ic[0] = M[0]-1;
    while (pbc[0] == true  && ic[0] <  0)     ic[0] += M[0];
    while (pbc[0] == true  && ic[0] >= M[0])  ic[0] -= M[0];

    while (pbc[1] == false && ic[1] <  0)     ic[1] = 0;
    while (pbc[1] == false && ic[1] >= M[1])  ic[1] = M[1]-1;
    while (pbc[1] == true  && ic[1] <  0)     ic[1] += M[1];
    while (pbc[1] == true  && ic[1] >= M[1])  ic[1] -= M[1];

    while (pbc[2] == false && ic[2] <  0)     ic[2] = 0;
    while (pbc[2] == false && ic[2] >= M[2])  ic[2] = M[2]-1;
    while (pbc[2] == true  && ic[2] <  0)     ic[2] += M[2];
    while (pbc[2] == true  && ic[2] >= M[2])  ic[2] -= M[2];

    /*
    for (k=0; k!=3; ++k){
      ic[k] = floor( rnew[k]/boxlen[k] * M[k] );
      if (pbc[k] == false){
	if      (ic[k] <  0)     ic[k] = 0;
	else if (ic[k] >= M[k])  ic[k] = M[k]-1;
      }
      else {
	while (ic[k] <  0)     ic[k] += M[k];
	while (ic[k] >= M[k])  ic[k] -= M[k];
      }
    }
    */
    

    // printf("Cell indices of atom %ld are %ld %ld %ld\n", i, ic[0], ic[1], ic[2]);
	
    /* Go through this cell and all surrounding cells looking for neighbors. */
    for (ix0=ic[0]-1; ix0<= ic[0]+1; ix0++){
      ix = ix0;

      if    (pbc[0] == false && ix < 0    ) continue;
      if    (pbc[0] == false && ix >= M[0]) continue;
      while (pbc[0] == true  && ix < 0    ) ix += M[0];
      while (pbc[0] == true  && ix >= M[0]) ix -= M[0];


      /*
	if (mpbc[1]==true){
	while (ix < 1)    ix += M[1];
	while (ix > M[1]) ix -= M[1];
	}
	else if (mpbc[1]==false && (ix < 1 || ix > M[1]))
	continue;
      */



      for (iy0=ic[1]-1; iy0<= ic[1]+1; iy0++){
	iy = iy0;
		

	if    (pbc[1] == false && iy < 0    ) continue;
	if    (pbc[1] == false && iy >= M[1]) continue;
	while (pbc[1] == true  && iy < 0    ) iy += M[1];
	while (pbc[1] == true  && iy >= M[1]) iy -= M[1];

	/*
	  if (mpbc[2]==true){
	  while (iy < 1)    iy += M[2];
	  while (iy > M[2]) iy -= M[2];
	  }
	  else if (mpbc[2]==false && (iy < 1 || iy > M[2]))
	  continue;
	*/

		
	for (iz0=ic[2]-1; iz0<= ic[2]+1; iz0++){
	  iz = iz0;

	  if    (pbc[2] == false && iz < 0    ) continue;
	  if    (pbc[2] == false && iz >= M[2]) continue;
	  while (pbc[2] == true  && iz < 0    ) iz += M[2];
	  while (pbc[2] == true  && iz >= M[2]) iz -= M[2];

	  /*
	    if (mpbc[3]==true){
	    while (iz < 1)    iz += M[3];
	    while (iz > M[3]) iz -= M[3];
	    }
	    else if (mpbc[3]==false && (iz < 1 || iz > M[3]))
	    continue;
	  */



		    
	  //		    printf("ic1  ic2  ic3  = %ld  %ld  %ld       ix iy iz = %ld  %ld  %ld\n",
	  //			   ic[1],ic[2],ic[3], ix,iy,iz);
		    


		    
	  /* State cell index. */
	  icell = ix + M[0]*( iy + M[1] * (iz) );
	  //		    printf("icell = %ld\n", icell);

	  // *******************************************************************		      
	  j = head[icell];
	  while (j>=0){

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
	    //cout << posi << posj << drsq << endl;

	    if (i!=j && drsq < rmsq){
	      if (counter==0 || (counter>0 && drsq<drsq_min)) drsq_min = drsq;

	      // Requires neighborcollection array to exist, otherwise fatal error:
	      /*
		if (nn >= mneighborcollection[i].size()){
		mneighborcollection[i].resize(nn+10);
		}
		mneighborcollection[i][nn] = j;
	      */
	      
	      /*
		cout << "Atom index " << i << ": Added neighbor atom index " 
		<< j << " at position " << nn << endl;
	      */
	      
	      neighborcollection[i].push_back(j);
	      nn_tot++;
	    }

	    j = list[j];
	  }
	  // *******************************************************************

	}
      }
    }
    
    // Sort: Only do this if a potential needs it! Saves time ...
    //std::sort( neighborcollection[i].begin(), neighborcollection[i].end() );
    
    //cout << "Size of neighborcollection of atom " << i << " is " << neighborcollection[i].size() << endl;

    /* Proceed to next atom i. */
    ++counter;
  }
  //cout << "Verlet-list completed." << endl;

  dr_min = sqrt(drsq_min);

  if (debug)
    cout << "AtomSystem: linked-list used  : <N_neighb/atom> " << nn_tot / (1.0*nat)
	 << " dr_min " << dr_min << endl;


  //finalize_frame_neighborcollection();

  delete [] head;
  delete [] list;

  return;

}      

