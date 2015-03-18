

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






void AtomSystem::dumpframe(ofstream & fout,
			   string fmt
			   ){
  int i;
  string fmt_now;
  const double llim=1.0e-4;
  const double ulim=1.0e+4;
  string formatf  = "%20.10f";
  string formatfs = "%.10f";
  string formate  = "%20.10e";
  double td;


  if (use_def_xyz_fmt) fmt_now = def_xyz_fmt;
  else                 fmt_now = fmt;

  /*
  elem.reserve(5);
  buf.reserve(100);
  */

  //natoms() = AtomSystem::get_natoms()();
  //AtomSystem::get_box(box);


  fout << natoms() << endl;

  if (fmt_now=="xyz"){
    fout << "Frame number " << iframe << " ";

    td = time;
    if (abs(td)<llim || abs(td)>ulim) fout << format(formate) % td;
    else                              fout << format(formatf) % td;
    fout << " fs boxsize ";

    td = boxlen[0];
    if (abs(td)<llim || abs(td)>ulim) fout << format(formate) % td;
    else                              fout << format(formatf) % td;
    fout << " ";

    td = boxlen[1];
    if (abs(td)<llim || abs(td)>ulim) fout << format(formate) % td;
    else                              fout << format(formatf) % td;
    fout << " ";

    td = boxlen[2];
    if (abs(td)<llim || abs(td)>ulim) fout << format(formate) % td;
    else                              fout << format(formatf) % td;
    fout << endl;

  }
  else if (fmt_now=="extxyz"){

    fout << "Lattice=\"";

    string boxstr = "";
    for (int i=0; i<3; ++i){
      Vector3<double> tv(0);
      get_boxdir(i, tv);

      td = boxlen[i] * tv[0];
      if (abs(td)<llim || abs(td)>ulim) boxstr += tostring_fmt(formate,  td) + " ";
      else                              boxstr += tostring_fmt(formatfs, td) + " ";

      td = boxlen[i] * tv[1];
      if (abs(td)<llim || abs(td)>ulim) boxstr += tostring_fmt(formate,  td) + " ";
      else                              boxstr += tostring_fmt(formatfs, td) + " ";

      td = boxlen[i] * tv[2];
      if (abs(td)<llim || abs(td)>ulim) boxstr += tostring_fmt(formate,  td);
      else                              boxstr += tostring_fmt(formatfs, td);

      if (i==2) boxstr += "\"";
      else      boxstr += " ";
    }
    fout << boxstr;
    fout << " Properties=species:S:1:pos:R:3:type:I:1:index:I:1 Time=0.0" << endl;
  }



  for (i = 0; i < natoms(); ++i){
    fout << matter[i] << " ";

    td = pos[i][0];
    if (abs(td)<llim || abs(td)>ulim) fout << format(formate) % td << " ";
    else                              fout << format(formatf) % td << " ";

    td = pos[i][1];
    if (abs(td)<llim || abs(td)>ulim) fout << format(formate) % td << " ";
    else                              fout << format(formatf) % td << " ";

    td = pos[i][2];
    if (abs(td)<llim || abs(td)>ulim) fout << format(formate) % td << " ";
    else                              fout << format(formatf) % td << " ";

    fout << type[i] << " "
	 << idx[i]  << endl;
    
  }
  

  return;
}





void AtomSystem::getframe(ifstream & fin,
			  string fmt
			  ){
  int ttype;
  int ns, atom_counter, tnatoms, iat;
  int tidx;
  double ttime;
  bool last_frame_line = false;
  vector<string> args;
  string line, tmatter, tfield;
  istringstream sstream;
  Vector3<double> tbox(0), tpos(0);

  atom_counter = 0;
  this->clear_all_atoms();



  while (true){

    // Get line from file:
    get_line( fin, line );
    // Error condition means one of two things:
    //   (1) File ended, but read some characters which should be used.
    //   (2) File has already ended, did not read any characters.
    // The following function for extracting substrings will return 0
    // for case (2), and nonzero otherwise.

    // Get words on line:
    ns = get_substrings(line, args, "\t ");

    // for (int i=0; i!=ns; ++i){
    // cout << "Word " << i << " is: " << args[i] << endl;
    // }

    if (ns==1){
      sstream.str(args[0]);
      sstream >> tnatoms;
      sstream.clear();

      // Next call will help avoid unncessary reallocations, don't need to guess
      // appropriate chunk size to optimize memory calls.
      // matoms.cap(tnatoms);
    }
    else if (ns>1 && args[0]=="Frame"){
      // cout << "args[3]: '" << args[3] << "'" << endl;
      // cout << "args[3].size() = " << args[3].size() << endl;

      sstream.str(args[3]);
      // cout << "String stored in stringstream is: '" << sstream.str() << "'" << endl;

      sstream >> ttime;
      // cout << "sstream.fail(): " << sstream.fail() << endl;
      // cout << "sstream.bad() : " << sstream.bad() << endl;
      // cout << "sstream.eof() : " << sstream.eof() << endl;
      // cout << "sstream.good(): " << sstream.good() << endl;
      sstream.clear();
      // cout << "Read time: " << ttime << endl;

      sstream.str(args[6]); sstream >> tbox[0]; sstream.clear();
      sstream.str(args[7]); sstream >> tbox[1]; sstream.clear();
      sstream.str(args[8]); sstream >> tbox[2]; sstream.clear();

      time = ttime;
      this->boxlen[0] = tbox[0];
      this->boxlen[1] = tbox[1];
      this->boxlen[2] = tbox[2];



      if (ns>=10 && args[9]=="boxdir"){
	if (ns<19){
	  cout << "Not enough coordinates for the box directions! Exiting." << endl;
	  exit(EXIT_FAILURE);
	}

	double b1,b2,b3;
	for (int i=0; i<3; ++i){
	  sstream.str(args[10+3*i+0]); sstream >> b1; sstream.clear();
	  sstream.str(args[10+3*i+1]); sstream >> b2; sstream.clear();
	  sstream.str(args[10+3*i+2]); sstream >> b3; sstream.clear();
	  Vector3<double> p(0);
	  p[0]=b1; p[1]=b2; p[2]=b3;
	  this->set_boxdir(i, p);
	  //this->boxdir(i, b1,b2,b3);
	}

      }
    }
    else if (ns>1){
      sstream.str(args[0]); sstream >> tmatter; sstream.clear();
      sstream.str(args[1]); sstream >> tpos[0]; sstream.clear();
      sstream.str(args[2]); sstream >> tpos[1]; sstream.clear();
      sstream.str(args[3]); sstream >> tpos[2]; sstream.clear();
      sstream.str(args[4]); sstream >> ttype; sstream.clear();
      sstream.str(args[5]); sstream >> tidx; sstream.clear();

      tfield.resize(0);
      for (int i=6; i!=ns; ++i){
	tfield += " " + string(args[i]);
      }

      iat = add_atom();
      //cout << "made it here gf 00, tmatter = " << tmatter << endl;
      matter[iat] = tmatter;
      //cout << "made it here gf 00b" << endl;
      pos[iat] = tpos;

      sitetype[iat] = iat;
      type[iat] = ttype;
      idx[iat] = tidx;
      //cout << "made it here gf 01" << endl;
      if (tfield.size()>0) field[iat] = tfield;
      //cout << "made it here gf 02" << endl;
      atom_counter++;

      /*
      cout << "1- ";
      Vector<double> tv;
      tv = at.pos();
      cout << " -1" << endl;

      cout << "Position of atom " << atom_counter-1 << " is:" << endl;
      cout << at.idx() << " " << tidx << tpos << matoms[atom_counter-1].pos() << endl;
      */


      if (atom_counter == tnatoms){
	// Read last atom in this frame.
	//this->natoms(tnatoms);
	last_frame_line = true;
	//matoms.chunksize(0);
	//matoms.trim();
      }

    }



    // ###########################################################
    // NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE 
    // ###########################################################


    
    if (! fin){ // tried to read beyond end of file
      fin.setstate(ios::failbit);
      break;
    }
    
    if (last_frame_line) break;
  }


  // Use the box information to update the geometry.
  this->update_box_geometry();

  /*
  cout << "Box lenghs:" << endl;
  cout << boxlen(0) << endl;
  cout << boxlen(1) << endl;
  cout << boxlen(2) << endl;

  cout << "Box geometry:" << endl;
  cout << boxdir(0) << endl;
  cout << boxdir(1) << endl;
  cout << boxdir(2) << endl;
  */

  return;
}



