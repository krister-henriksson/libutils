

#include <cmath>



#include "constants.hpp"
#include "utils-vecalg.hpp"
#include "utils-errors.hpp"


using namespace constants;
using utils::aborterror;



// #####################################################################
// Get rotation matrices
// #####################################################################

/* ###############################################################
   Rotation around one of the Cartesian axes.
   Rotation axis passes through the origin (0,0,0).
   ###############################################################
 */
void utils::get_rotation_matrix(Matrix<double> & R,
				const double & th,
				const int axis
				){
  if (R.ncols()!=3 && R.nrows()!=3) R.resize(3,3);
  for (int i=0; i<3; ++i) for (int j=0; j<3; ++j) R.elem(i,j)=0.0;
    
  double c = cos(th);
  double s = sin(th);

  if (axis==0){
    // Rotation axis = x axis, goes through (0,0,0).
    R.elem(0,0) = 1;
    R.elem(0,1) = 0;
    R.elem(0,2) = 0;

    R.elem(1,0) = 0;
    R.elem(1,1) = c;
    R.elem(1,2) = -s;

    R.elem(2,0) = 0;
    R.elem(2,1) = s;
    R.elem(2,2) = c;
  }
  else if (axis==1){
    // Rotation axis = y axis, goes through (0,0,0).
    R.elem(0,0) = c;
    R.elem(0,1) = 0;
    R.elem(0,2) = s;
  
    R.elem(1,0) = 0;
    R.elem(1,1) = 1;
    R.elem(1,2) = 0;
  
    R.elem(2,0) = -s;
    R.elem(2,1) = 0;
    R.elem(2,2) = c;
  }
  else {
    // Rotation axis = z axis, goes through (0,0,0).
    R.elem(0,0) = c;
    R.elem(0,1) = -s;
    R.elem(0,2) = 0;
  
    R.elem(1,0) = s;
    R.elem(1,1) = c;
    R.elem(1,2) = 0;
  
    R.elem(2,0) = 0;
    R.elem(2,1) = 0;
    R.elem(2,2) = 1;
  }

  return;
}

void utils::get_improper_rotation_matrix(Matrix<double> & R,
					 const double & th,
					 const int axis
					 ){
  Matrix<double> I(3,3,0), IR(3,3,0);
  I.elem(0,0) = I.elem(1,1) = I.elem(2,2) = -1.0;

  utils::get_rotation_matrix(R, th, axis);
  IR = I * R;
  R = IR;
}



/* ###############################################################
   Simple reflection in a plane represented by its normal vector.
   The normal vector coincides with a Cartesian axis.
   The plane passes through the origin (0,0,0).
   ###############################################################
 */
void utils::get_reflection_matrix(Matrix<double> & R,
				  const int axis
				  ){
  if (R.ncols()!=3 && R.nrows()!=3) R.resize(3,3);
  for (int i=0; i<3; ++i) for (int j=0; j<3; ++j) R.elem(i,j)=0.0;

  if (axis==0){
    R.elem(0,0) = -1;
    R.elem(0,1) = 0;
    R.elem(0,2) = 0;

    R.elem(1,0) = 0;
    R.elem(1,1) = 1;
    R.elem(1,2) = 0;

    R.elem(2,0) = 0;
    R.elem(2,1) = 0;
    R.elem(2,2) = 1;
  }
  else if (axis==1){
    R.elem(0,0) = 1;
    R.elem(0,1) = 0;
    R.elem(0,2) = 0;

    R.elem(1,0) = 0;
    R.elem(1,1) = -1;
    R.elem(1,2) = 0;

    R.elem(2,0) = 0;
    R.elem(2,1) = 0;
    R.elem(2,2) = 1;
  }
  else {
    R.elem(0,0) = 1;
    R.elem(0,1) = 0;
    R.elem(0,2) = 0;

    R.elem(1,0) = 0;
    R.elem(1,1) = 1;
    R.elem(1,2) = 0;

    R.elem(2,0) = 0;
    R.elem(2,1) = 0;
    R.elem(2,2) = -1;
  }

  return;
}





/* ###############################################################
   Rotation around an axis specified by a direction vector.
   Rotation axis passes through the origin (0,0,0).
   ###############################################################
 */
void utils::get_rotation_matrix_u(Matrix<double> & R,
				  const Vector<double> & u,
				  const double & th
				  ){
  if (R.ncols()!=3 && R.nrows()!=3) R.resize(3,3);
  for (int i=0; i<3; ++i) for (int j=0; j<3; ++j) R.elem(i,j)=0.0;

  if (u.size()!=3) aborterror("ERROR: get_rotation_matrix vector u does not contain 3 elements!");

  Vector<double> n = u;
  n.normalize();

  double c = cos(th);
  double s = sin(th);
  double nx = n[0];
  double ny = n[1];
  double nz = n[2];

  Matrix<double> K(3,3,0), K2(3,3,0), I(3,3,0);
  K.elem(0,0) =  0.0;
  K.elem(0,1) = -nz;
  K.elem(0,2) =  ny;
  K.elem(1,0) =  nz;
  K.elem(1,1) =  0.0;
  K.elem(1,2) = -nx;
  K.elem(2,0) = -ny;
  K.elem(2,1) =  nx;
  K.elem(2,2) =  0.0;
  K2 = K * K;
  I.elem(0,0) = I.elem(1,1) = I.elem(2,2) = 1.0;


  R = I + s * K + (1.0 - c) * K2;
  Vector<double> tmp(3,0); tmp[0]=tmp[1]=tmp[2]=0.5;
  //cout << "Rotation matrix times (0.5, 0.5, 0.5) = " << R * tmp << endl;

  return;
}




void utils::get_improper_rotation_matrix_u(Matrix<double> & R,
					   const Vector<double> & u,
					   const double & th
					   ){

  Matrix<double> I(3,3,0), IR(3,3,0);
  I.elem(0,0) = I.elem(1,1) = I.elem(2,2) = -1.0;

  utils::get_rotation_matrix_u(R, u, th);
  IR = I * R;
  R = IR;

  return;
}





/* ###############################################################
   Rotation around an axis specified by a direction vector.
   Rotation axis passes through the point A.
   Vector to be rotated starts at the point B.
   Rotated vector starts at the point 'OB_rot'.
   ###############################################################
 */

// 1. Call: get_rotation_matrix_u(R, axis, th) to get R.
// 2. Call this function using this R:
void utils::rotate_around_axis_u(const Matrix<double> & R,
				 const Vector<double> & axis,
				 const Vector<double> & OA, // axis_point
				 const Vector<double> & vec,
				 const Vector<double> & OB, // vec_point
				 Vector<double> & vec_rot,
				 Vector<double> & OB_rot
				 ){

  Vector<double> AB(3,0), AD(3,0), DB(3,0), DB_rot(3,0);
  Vector<double> kvec(axis);

  // Rotate vector 'vec':
  vec_rot  = R * vec;

  rotate_vec_start_point(R, kvec, OA, OB, OB_rot);
}






/* ###############################################################
   Rotation so that specified vector 'vec' coincides with direction 'uvec'
   after rotation.
   The vector 'vec' will be rotated by an angle around the axis
   'nvec' = 'vec' x 'uvec' (vector product).
   The angle is 'th' = acos( 'vec' * 'uvec' / (|'vec'| * |'uvec'|) ).
   The axis 'nvec' passes through point E.
   Vector to be rotated starts at the point B.
   Rotated vector starts at the point 'OB_rot'.
   ###############################################################
 */

void utils::get_matrix_for_rotation_to_coincide_with_axis(Matrix<double> & R,
							  const Vector<double> & dir,
							  const Vector<double> & vec_to_align
							  ){
  Vector<double> nvec(3,0);
  double th;

  if (R.nrows()!=3 || R.ncols()!=3) R.resize(3,3);

  vectorproduct(vec_to_align, dir, nvec);
  nvec.normalize();

  th = acos( scalarproduct(vec_to_align, dir) / ( vec_to_align.magn() * dir.magn() ));
  if (th<0.0) th += 2*PI;
  // cout << "th = " << th << endl;
  get_rotation_matrix_u(R, nvec, th);

}


// 1. Call: get_matrix_for_rotation_to_coincide_with_axis(R, dir, vec_to_align) to get R.
// 2. Call this using this R:
void utils::rotate_to_coincide_with_axis(Matrix<double> & R,
					 const Vector<double> & dir,
					 const Vector<double> & vec_to_align,
					 // vec_to_align starts at this point:
					 const Vector<double> & OB,
					 // the normal to the plane formed by dir and
					 // vec_to_align passes through this point:
					 const Vector<double> & OA,
					 // rotated vec_to_align will start at this point:
					 Vector<double> & OB_rot
					 ){

  Vector<double> nvec(3,0);

  vectorproduct(vec_to_align, dir, nvec);

  rotate_vec_start_point(R, nvec, OA, OB, OB_rot);
}



void utils::rotate_vec_start_point(const Matrix<double> & R,
				   const Vector<double> & axis,
				   const Vector<double> & axis_point, // OA,
				   const Vector<double> & vec_start,  // OB
				   Vector<double> & vec_start_rot     // OB_rot
				   ){

  Vector<double> nvec(3,0), OA(3,0), OB(3,0), AB(3,0), AD(3,0), DB(3,0), DB_rot(3,0), OB_rot(3,0);

  nvec = axis;
  nvec.normalize();

  OA = axis_point;
  OB = vec_start;

  // Vector from point on axis:
  AB = OB - OA;
  // Projection onto axis:
  AD = scalarproduct(AB, nvec) * nvec;
  // The rejection part:
  DB = AB - AD;
  // Rotate vector:
  DB_rot = R * DB;
  OB_rot = OA + AD + DB_rot;

  vec_start_rot = OB_rot;
}


