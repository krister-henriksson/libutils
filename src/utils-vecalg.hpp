


#ifndef UTILS_VECALG_HPP
#define UTILS_VECALG_HPP



#include "utils-matrix.hpp"
#include "utils-vector.hpp"

using utils::Matrix;
using utils::Vector;



namespace utils {




  // #####################################################################
  // Get rotation matrices
  // #####################################################################

  void get_rotation_matrix(Matrix<double> & R,
			   const double & th,
			   const int axis);

  void get_improper_rotation_matrix(Matrix<double> & R,
				    const double & th,
				    const int axis);

  void get_reflection_matrix(Matrix<double> & R,
			     const int axis);

  void get_rotation_matrix_u(Matrix<double> & R,
			     const Vector<double> & u,
			     const double & th
			     );

  void get_improper_rotation_matrix_u(Matrix<double> & R,
				      const Vector<double> & u,
				      const double & th
				      );

  void rotate_around_axis_u(const Matrix<double> & R,
			    const Vector<double> & axis,
			    const Vector<double> & OA, // axis_point
			    const Vector<double> & vec,
			    const Vector<double> & OB, // vec_point
			    Vector<double> & vec_rot,
			    Vector<double> & OB_rot
			    );

  void get_matrix_for_rotation_to_coincide_with_axis(Matrix<double> & R,
						     const Vector<double> & dir,
						     const Vector<double> & vec_to_align
						     );


  void rotate_to_coincide_with_axis(Matrix<double> & R,
				    const Vector<double> & dir,
				    const Vector<double> & vec_to_align,
				    // vec_to_align starts at this point:
				    const Vector<double> & OB,
				    // the normal to the plane formed by dir and
				    // vec_to_align passes through this point:
				    const Vector<double> & OE,
				    // rotated vec_to_align will start at this point:
				    Vector<double> & OB_rot
				    );

  void rotate_vec_start_point(const Matrix<double> & R,
			      const Vector<double> & axis,
			      const Vector<double> & axis_point, // OA,
			      const Vector<double> & vec_start,  // OB
			      Vector<double> & vec_start_rot     // OB_rot
			      );


}

#endif



