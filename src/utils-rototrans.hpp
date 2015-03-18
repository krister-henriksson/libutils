

#ifndef UTILS_ROTOTRANS_HPP
#define UTILS_ROTOTRANS_HPP



#include "utils-vector3.hpp"
#include "utils-matrixsq3.hpp"
#include "utils-matrix.hpp"

using utils::Vector3;
using utils::MatrixSq3;
using utils::Matrix;


namespace utils {


  class RotoTransOperand {
  private:
    Matrix<double> mmat;

  public:
    Vector3<double> v;
    Vector3<double> OB;
    
    RotoTransOperand();
    RotoTransOperand(const Vector3<double> & iv,
		     const Vector3<double> & iOB = Vector3<double>(0));
    Vector3<double> & vector(void);
    Vector3<double> & vectorpoint(void);

  } ;






  class RotoTransOperator {
  private:
    bool axis_known;
    bool angle_known;
    Matrix<double> mmat;
    void get_rotation_matrix(void);

  public:
    MatrixSq3<double> R;
    Vector3<double> T;
    Vector3<double> n;
    Vector3<double> OA;
    double angle;
    
    RotoTransOperator();
    RotoTransOperator(const Vector3<double> & in,
		      const double th,
		      const Vector3<double> & iOA = Vector3<double>(0));
    RotoTransOperator(const Vector3<double> & iT);


    MatrixSq3<double> & rotation_matrix(void);
    Vector3<double>   & rotation_axis(void);
    double            & rotation_angle(void);
    Vector3<double>   & rotation_axispoint(void);
    Vector3<double>   & translation(void);

    void rotate_vector(Vector3<double> & v,
		       Vector3<double> & vrot);
    void rotate_vector_startpoint(Vector3<double> & P,
				  Vector3<double> & Prot);
    void rotate(RotoTransOperand & opv,
		RotoTransOperand & opv_rot);

  } ;


}


#endif

