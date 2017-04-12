


#ifndef ATOMSYS_HPP
#define ATOMSYS_HPP


#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "utils.hpp"
#include "utils-vector.hpp"
#include "utils-vector3.hpp"
#include "utils-matrix.hpp"
#include "utils-matrixsq3.hpp"

#include "bond.hpp"

#include "omp-basics.hpp"

using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::ofstream;
using std::ifstream;

using utils::Vector;
using utils::Vector3;
using utils::Matrix;
using utils::MatrixSq3;






class AtomSystem
//  : public CoordSys
{
private:
  int nat;

public:

  // System-related properties:
  double time;
  double rcut;
  double drcut;
  double skint;
  double vol;
  double dr_min;
  double vol_atom;
  int iframe;
  bool isCart;


  Vector3<double> boxlen;
  MatrixSq3<double> boxdir;
  MatrixSq3<double> Bravaismatrix_inv;
  Vector3<bool> pbc;

  // Atom-related properties:
  Vector<int> sitetype;
  Vector<int> type;  // external (arbitrary) type, e.g. specified by user
  Vector<int> itype; // internal type, e.g. index of element in array of element names
  Vector<int> idx;
  Vector<string> matter;
  Vector<string> field;
  Vector< Vector3<double> > pos;
  // constraints
  Vector< bool >           atom_is_fixed;
  Vector< Vector<double> > atom_freedir;
  Vector< Vector<double> > atom_freeplane;
  
  Vector< Vector<int> > neighborcollection;

 
  OMP_Info omp_info;

  bool   use_def_xyz_fmt;
  string def_xyz_fmt;



public:
  AtomSystem();
  virtual ~AtomSystem();

  //AtomSystem(const AtomSystem & sys);
  //AtomSystem & operator=(const AtomSystem & sys);

  // ##################################################################################
  // ##################################################################################

  // Enables sys.clear_all_atoms() before adding any atoms.
  void clear_all_atoms();

  void init_atoms(int n);
  void finalize_atoms();

  // Insert atom info:
  int add_atom(); // returns vector index of atom, use to update position etc.

  inline int natoms() const {
    return pos.size();
  }

  // ##################################################################################
  // ##################################################################################

  // Set box direction vectors and length of box in the specified direction. The vectors
  // are normalized inside the method.
  void set_boxdir(const int idir,
		  const double & px,
		  const double & py,
		  const double & pz);
  void set_boxdir(const int idir, Vector3<double> & p);
  // Get box direction vectors:
  void get_boxdir(const int idir, Vector3<double> & v) const;

  // Normalize the box direction vectors and build the Bravais matrix and its inverse:
  void update_box_geometry();

  void calc_volume(void);

  void calc_closepacked_volume(void);



  // ##################################################################################
  // ##################################################################################

  void get_atom_distance_vec(const Vector3<double> & r1,
			     const Vector3<double> & r2,
			     Vector3<double> & v,
			     const double lowlim=-1) const;
  
  inline double get_atom_distance(const Vector3<double> & r1,
				  const Vector3<double> & r2,
				  const double lowlim=-1) const {
    Vector3<double> drc(0.0);
    get_atom_distance_vec(r1, r2, drc, lowlim);
    return sqrt( drc[0]*drc[0] + drc[1]*drc[1] + drc[2]*drc[2] );
  }

  void get_coords_cart2skew(const Vector3<double> & drc,
			    Vector3<double> & v,
			    const double lowlim=-1) const;

  void get_coords_skew2cart(Vector3<double> & drs,
			    Vector3<double> & v,
			    const double lowlim=-1) const;
  
  void translate_cartpos_in_skewspace(Vector3<double> & pos,
				      double f1,
				      double f2,
				      double f3);

  void handle_pbc_of_positions(const double lowlim=-1);


  // ##################################################################################
  // ##################################################################################

  // Build all neighbor collections:
  void get_all_neighborcollections(bool debug=false);

  // ##################################################################################
  // ##################################################################################


  void dumpframe(ofstream & fout, string fmt="xyz");
  void getframe(ifstream & fin,   string fmt="xyz");

  // ##################################################################################
  // ##################################################################################


  void get_bond_list(BondData & bond_list,
		     BondData & bond_list12,
		     BondData & bond_list21,
		     string & name1,
		     string & name2,
		     int & nat1,
		     int & nat2,
		     int & nbonds,
		     int & nbonds12,
		     int & nbonds21,
		     double rc12
		     );

  void get_bond_angle_list(BondAngleData & bondangle_list,
			   string & name1,
			   string & name2,
			   double rc11,
			   double rc22,
			   double rc1221
			   );



} ;









#endif


