

#ifndef UTILS_ROTOTRANS_HPP
#define UTILS_ROTOTRANS_HPP


#include "utils-vector.hpp"
#include "utils-matrix.hpp"

using utils::Vector;
using utils::Matrix;


namespace utils {


  class RotoTransOperand {
  private:
    Matrix<double> mmat;

  public:
    Vector<double> v;
    Vector<double> OB;
    
    RotoTransOperand();
    RotoTransOperand(const Vector<double> & iv,
		     const Vector<double> & iOB = Vector<double>(3,0));
    Vector<double> & vector(void);
    Vector<double> & vectorpoint(void);

  } ;






  class RotoTransOperator {
  private:
    bool axis_known;
    bool angle_known;
    Matrix<double> mmat;
    void get_rotation_matrix(void);

  public:
    Matrix<double> R;
    Vector<double> T;
    Vector<double> n;
    Vector<double> OA;
    double angle;
    
    RotoTransOperator();
    RotoTransOperator(const Vector<double> & in,
		      const double th,
		      const Vector<double> & iOA = Vector<double>(3,0));
    RotoTransOperator(const Vector<double> & iT);


    Matrix<double> & rotation_matrix(void);
    Vector<double> & rotation_axis(void);
    double         & rotation_angle(void);
    Vector<double> & rotation_axispoint(void);
    Vector<double> & translation(void);

    void rotate_vector(Vector<double> & v,
		       Vector<double> & vrot);
    void rotate_vector_startpoint(Vector<double> & P,
				  Vector<double> & Prot);
    void rotate(RotoTransOperand & opv,
		RotoTransOperand & opv_rot);

  } ;


}


#endif

