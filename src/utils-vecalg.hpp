


#ifndef UTILS_VECALG_HPP
#define UTILS_VECALG_HPP



#include "utils-matrixsq3.hpp"
#include "utils-vector3.hpp"

using utils::MatrixSq3;
using utils::Vector3;



namespace utils {




  // #####################################################################
  // Get rotation matrices
  // #####################################################################

  void get_rotation_matrix(MatrixSq3<double> & R,
			   const double & th,
			   const int axis);

  void get_improper_rotation_matrix(MatrixSq3<double> & R,
				    const double & th,
				    const int axis);

  void get_reflection_matrix(MatrixSq3<double> & R,
			     const int axis);

  void get_rotation_matrix_u(MatrixSq3<double> & R,
			     const Vector3<double> & u,
			     const double & th
			     );

  void get_improper_rotation_matrix_u(MatrixSq3<double> & R,
				      const Vector3<double> & u,
				      const double & th
				      );

  void rotate_around_axis_u(const MatrixSq3<double> & R,
			    const Vector3<double> & axis,
			    const Vector3<double> & OA, // axis_point
			    const Vector3<double> & vec,
			    const Vector3<double> & OB, // vec_point
			    Vector3<double> & vec_rot,
			    Vector3<double> & OB_rot
			    );

  void get_matrix_for_rotation_to_coincide_with_axis(MatrixSq3<double> & R,
						     const Vector3<double> & dir,
						     const Vector3<double> & vec_to_align
						     );


  void rotate_to_coincide_with_axis(MatrixSq3<double> & R,
				    const Vector3<double> & dir,
				    const Vector3<double> & vec_to_align,
				    // vec_to_align starts at this point:
				    const Vector3<double> & OB,
				    // the normal to the plane formed by dir and
				    // vec_to_align passes through this point:
				    const Vector3<double> & OE,
				    // rotated vec_to_align will start at this point:
				    Vector3<double> & OB_rot
				    );

  void rotate_vec_start_point(const MatrixSq3<double> & R,
			      const Vector3<double> & axis,
			      const Vector3<double> & axis_point, // OA,
			      const Vector3<double> & vec_start,  // OB
			      Vector3<double> & vec_start_rot     // OB_rot
			      );


}

#endif



