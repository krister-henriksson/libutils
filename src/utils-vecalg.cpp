

#include <cmath>



#include "constants.hpp"
#include "utils-errors.hpp"

#include "utils-vecalg.hpp"

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
void utils::get_rotation_matrix(MatrixSq3<double> & R,
				const double & th,
				const int axis
				){
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

void utils::get_improper_rotation_matrix(MatrixSq3<double> & R,
					 const double & th,
					 const int axis
					 ){
  MatrixSq3<double> I(0), IR(0);
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
void utils::get_reflection_matrix(MatrixSq3<double> & R,
				  const int axis
				  ){
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
void utils::get_rotation_matrix_u(MatrixSq3<double> & R,
				  const Vector3<double> & u,
				  const double & th
				  ){
  for (int i=0; i<3; ++i) for (int j=0; j<3; ++j) R.elem(i,j)=0.0;

  Vector3<double> n = u;
  n.normalize();

  double c = cos(th);
  double s = sin(th);
  double nx = n[0];
  double ny = n[1];
  double nz = n[2];

  MatrixSq3<double> K(0), K2(0), I(0);
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
  Vector3<double> tmp(0.5);// tmp[0]=tmp[1]=tmp[2]=0.5;
  //cout << "Rotation matrix times (0.5, 0.5, 0.5) = " << R * tmp << endl;

  return;
}




void utils::get_improper_rotation_matrix_u(MatrixSq3<double> & R,
					   const Vector3<double> & u,
					   const double & th
					   ){

  MatrixSq3<double> I(0), IR(0);
  I.elem(0,0) = I.elem(1,1) = I.elem(2,2) = -1.0;

  utils::get_rotation_matrix_u(R, u, th);
  IR = I * R;
  R = IR;

  return;
}





/* ###############################################################
   Rotation around an axis specified by a direction vector.
   Rotation axis passes through the point A.
   Vector3 to be rotated starts at the point B.
   Rotated vector starts at the point 'OB_rot'.
   ###############################################################
 */

// 1. Call: get_rotation_matrix_u(R, axis, th) to get R.
// 2. Call this function using this R:
void utils::rotate_around_axis_u(const MatrixSq3<double> & R,
				 const Vector3<double> & axis,
				 const Vector3<double> & OA, // axis_point
				 const Vector3<double> & vec,
				 const Vector3<double> & OB, // vec_point
				 Vector3<double> & vec_rot,
				 Vector3<double> & OB_rot
				 ){

  Vector3<double> AB(0), AD(0), DB(0), DB_rot(0);
  Vector3<double> kvec(axis);

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
   Vector3 to be rotated starts at the point B.
   Rotated vector starts at the point 'OB_rot'.
   ###############################################################
 */

void utils::get_matrix_for_rotation_to_coincide_with_axis(MatrixSq3<double> & R,
							  const Vector3<double> & dir,
							  const Vector3<double> & vec_to_align
							  ){
  Vector3<double> nvec(0);
  double th;

  vectorproduct(vec_to_align, dir, nvec);
  nvec.normalize();

  th = acos( scalarproduct(vec_to_align, dir) / ( vec_to_align.magn() * dir.magn() ));
  if (th<0.0) th += 2*PI;
  // cout << "th = " << th << endl;
  get_rotation_matrix_u(R, nvec, th);

}


// 1. Call: get_matrix_for_rotation_to_coincide_with_axis(R, dir, vec_to_align) to get R.
// 2. Call this using this R:
void utils::rotate_to_coincide_with_axis(MatrixSq3<double> & R,
					 const Vector3<double> & dir,
					 const Vector3<double> & vec_to_align,
					 // vec_to_align starts at this point:
					 const Vector3<double> & OB,
					 // the normal to the plane formed by dir and
					 // vec_to_align passes through this point:
					 const Vector3<double> & OA,
					 // rotated vec_to_align will start at this point:
					 Vector3<double> & OB_rot
					 ){

  Vector3<double> nvec(0);

  vectorproduct(vec_to_align, dir, nvec);

  rotate_vec_start_point(R, nvec, OA, OB, OB_rot);
}



void utils::rotate_vec_start_point(const MatrixSq3<double> & R,
				   const Vector3<double> & axis,
				   const Vector3<double> & axis_point, // OA,
				   const Vector3<double> & vec_start,  // OB
				   Vector3<double> & vec_start_rot     // OB_rot
				   ){

  Vector3<double> nvec(0), OA(0), OB(0), AB(0), AD(0), DB(0), DB_rot(0), OB_rot(0);

  nvec = axis;
  nvec.normalize();

  OA = axis_point;
  OB = vec_start;

  // Vector3 from point on axis:
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


