


#ifndef ATOMSYS_HPP
#define ATOMSYS_HPP


#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "utils.hpp"
#include "utils-vector.hpp"
#include "utils-matrix.hpp"

#include "bond.hpp"

#include "omp-basics.hpp"

using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::ofstream;
using std::ifstream;

using utils::Vector;
using utils::Matrix;






class AtomSystem
//  : public CoordSys
{
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


  Vector<double> boxlen;
  Matrix<double> boxdir;
  Matrix<double> Bravaismatrix_inv;
  Vector<bool> pbc;

  // Atom-related properties:
  Vector<int> sitetype;
  Vector<int> type;
  Vector<int> idx;
  Vector<string> matter;
  Vector<string> field;
  Vector< Vector<double> > pos;
  
  Vector< Vector<int> > neighborcollection;
 
  OMP_Info omp_info;

  bool   use_def_xyz_fmt;
  string def_xyz_fmt;



public:
  AtomSystem();
  virtual ~AtomSystem();
  AtomSystem(const AtomSystem & sys);
  
  AtomSystem & operator=(const AtomSystem & sys);

  // ##################################################################################
  // ##################################################################################

  // Enables sys.clear_all_atoms() before adding any atoms.
  void clear_all_atoms();
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
  void set_boxdir(const int idir, Vector<double> & p);
  // Get box direction vectors:
  void get_boxdir(const int idir, Vector<double> & v) const;

  // Normalize the box direction vectors and build the Bravais matrix and its inverse:
  void update_box_geometry();

  void calc_volume(void);

  void calc_closepacked_volume(void);


  // ##################################################################################
  // ##################################################################################

  void get_atom_distance_vec(const Vector<double> & r1,
			     const Vector<double> & r2,
			     Vector<double> & v,
			     const double lowlim=-1) const;
  
  inline double get_atom_distance(const Vector<double> & r1,
				  const Vector<double> & r2,
				  const double lowlim=-1) const {
    Vector<double> drc(3, 0.0);
    get_atom_distance_vec(r1, r2, drc, lowlim);
    return sqrt( drc[0]*drc[0] + drc[1]*drc[1] + drc[2]*drc[2] );
  }

  void get_coords_cart2skew(const Vector<double> & drc,
			    Vector<double> & v,
			    const double lowlim=-1) const;

  void get_coords_skew2cart(Vector<double> & drs,
			    Vector<double> & v,
			    const double lowlim=-1) const;
  
  void translate_cartpos_in_skewspace(Vector<double> & pos,
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


  void get_bond_list(Vector<BondData> & bond_list,
		     string & name1,
		     string & name2,
		     int & nat_with_bonds,
		     double rc12
		     );

  void get_bond_angle_list(Vector<BondAngleData> & bondangle_list,
			   string & name1,
			   string & name2,
			   double rc11,
			   double rc22,
			   double rc12
			   );



} ;









#endif


