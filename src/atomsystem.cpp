

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




AtomSystem::AtomSystem()
  :
  nat(0),
  time(0.0),
  rcut(0.0),
  drcut(0.0),
  skint(0.0),
  vol(-1.0),
  dr_min(-1.0),
  vol_atom(-1.0),
  iframe(0),
  isCart(true),
  boxlen(0.0),
  boxdir(0.0),
  Bravaismatrix_inv(0.0),
  pbc(true)
{
  nat = 0;

  // Default to Cartesian box:
  boxdir.elem(0,0) = 1.0;
  boxdir.elem(1,1) = 1.0;
  boxdir.elem(2,2) = 1.0;
  Bravaismatrix_inv.elem(0,0) = 1.0;
  Bravaismatrix_inv.elem(1,1) = 1.0;
  Bravaismatrix_inv.elem(2,2) = 1.0;

  sitetype.cap(0);
  type.cap(0);
  itype.cap(0);
  idx.cap(0);
  matter.cap(0);
  field.cap(0);
  pos.cap(0);
  atom_is_fixed.cap(0);
  atom_freedir.cap(0);
  atom_freeplane.cap(0);
  neighborcollection.cap(0);

  use_def_xyz_fmt = false;
  def_xyz_fmt = "xyz";
}


/*
AtomSystem::AtomSystem(const AtomSystem & sys){
  // Default copy constructors for standard vectors don't work properly
  // if the source/destination is not initialized:

  nat = sys.nat;
  time = sys.time;
  rcut = sys.rcut;
  drcut = sys.drcut;
  skint = sys.skint;
  iframe = sys.iframe;
  isCart = sys.isCart;

  boxlen = sys.boxlen;
  boxdir = sys.boxdir;
  Bravaismatrix_inv = sys.Bravaismatrix_inv;
  pbc = sys.pbc;

  sitetype = sys.sitetype;
  type = sys.type;
  itype = sys.itype;
  idx = sys.idx;
  pos = sys.pos;
  matter = sys.matter;
  field = sys.field;

  atom_is_fixed  = sys.atom_is_fixed;
  atom_freedir   = sys.atom_freedir;
  atom_freeplane = sys.atom_freeplane;

  neighborcollection = sys.neighborcollection;

  omp_info = sys.omp_info;

  use_def_xyz_fmt = sys.use_def_xyz_fmt;
  def_xyz_fmt = sys.def_xyz_fmt;
}
*/



AtomSystem::~AtomSystem(){ }



/*
AtomSystem & AtomSystem::operator=(const AtomSystem & sys){
    
  // watch out for self-assignment!
  if (this == &sys) return *this;

  nat = sys.nat;
  time = sys.time;
  rcut = sys.rcut;
  drcut = sys.drcut;
  skint = sys.skint;
  iframe = sys.iframe;
  isCart = sys.isCart;

  boxlen = sys.boxlen;
  boxdir = sys.boxdir;
  Bravaismatrix_inv = sys.Bravaismatrix_inv;
  pbc = sys.pbc;

  sitetype = sys.sitetype;
  type = sys.type;
  itype = sys.itype;
  idx = sys.idx;
  pos = sys.pos;
  matter = sys.matter;
  field = sys.field;

  atom_is_fixed  = sys.atom_is_fixed;
  atom_freedir   = sys.atom_freedir;
  atom_freeplane = sys.atom_freeplane;

  neighborcollection = sys.neighborcollection;

  omp_info = sys.omp_info;

  use_def_xyz_fmt = sys.use_def_xyz_fmt;
  def_xyz_fmt = sys.def_xyz_fmt;

  return *this;
}
*/
  


// ##################################################################################
// ##################################################################################

// Set the number of atoms in the system to 0.
void AtomSystem::clear_all_atoms(){
  nat=0;
  sitetype.resize(0);
  type.resize(0);
  itype.resize(0);
  idx.resize(0);
  pos.resize(0);
  matter.resize(0);
  field.resize(0);
  atom_is_fixed.resize(0);
  atom_freedir.resize(0);
  atom_freeplane.resize(0);
  neighborcollection.resize(0);
}



void AtomSystem::init_atoms(int n){
  sitetype.reserve_fraction(0.0);
  type.reserve_fraction(0.0);
  itype.reserve_fraction(0.0);
  idx.reserve_fraction(0.0);
  pos.reserve_fraction(0.0);
  matter.reserve_fraction(0.0);
  field.reserve_fraction(0.0);
  atom_is_fixed.reserve_fraction(0.0);
  atom_freedir.reserve_fraction(0.0);
  atom_freeplane.reserve_fraction(0.0);
  neighborcollection.reserve_fraction(0.0);

  sitetype.cap(n);
  type.cap(n);
  itype.cap(n);
  idx.cap(n);
  pos.cap(n);
  matter.cap(n);
  field.cap(n);
  atom_is_fixed.cap(n);
  atom_freedir.cap(n);
  atom_freeplane.cap(n);
  neighborcollection.cap(n);
}

void AtomSystem::finalize_atoms(){
  sitetype.trim();
  type.trim();
  itype.trim();
  idx.trim();
  pos.trim();
  matter.trim();
  field.trim();
  atom_is_fixed.trim();
  atom_freedir.trim();
  atom_freeplane.trim();
  neighborcollection.trim();
}


// Add an atom to the system. It is placed at the end of the matoms vector,
// and the mnatoms counter is updated by 1.
int AtomSystem::add_atom(){
  nat++;

  sitetype.resize(nat);
  type.resize(nat);
  itype.resize(nat);
  idx.resize(nat);
  pos.resize(nat);
  matter.resize(nat);
  field.resize(nat);
  atom_is_fixed.resize(nat);
  atom_freedir.resize(nat);
  atom_freeplane.resize(nat);
  neighborcollection.resize(nat);

  return nat-1;
}




