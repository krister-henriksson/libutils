


#include "utils-rototrans.hpp"


using utils::RotoTransOperand;
using utils::RotoTransOperator;



utils::RotoTransOperand::RotoTransOperand(){
  v  = Vector3<double>(0.0);
  OB = Vector3<double>(0.0);

  mmat.resize(4, 2);
  for (int i=0; i<4; ++i){
    for (int j=0; j<2; ++j){
      mmat.elem(i,j) = 0.0;
    }
  }
}


utils::RotoTransOperand::RotoTransOperand(const Vector3<double> & iv,
					  const Vector3<double> & iOB){
  v = iv;
  OB = iOB;
}


Vector3<double> & utils::RotoTransOperand::vector(void){ return v; }

Vector3<double> & utils::RotoTransOperand::vectorpoint(void){ return OB; }



// #######################################################################
// #######################################################################
// #######################################################################
// #######################################################################



void utils::RotoTransOperator::get_rotation_matrix(void){
  // Need axis and angle before matrix can be specified.
  n.normalize();

  double c = cos(angle);
  double s = sin(angle);
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
}



utils::RotoTransOperator::RotoTransOperator(){
  R  = MatrixSq3<double>(0.0);
  R.elem(0,0) = R.elem(1,1) = R.elem(2,2) = 1.0;
  T  = Vector3<double>(0.0);
  n  = Vector3<double>(0.0);
  OA = Vector3<double>(0.0);
  angle = 0.0;

  mmat.resize(4, 6);
  for (int i=0; i<4; ++i){
    for (int j=0; j<6; ++j){
      mmat.elem(i,j) = 0.0;
    }
  }
  for (int i=0; i<3; ++i){
    for (int j=0; j<3; ++j){
      mmat.elem(i,j) = R.elem(i,j);
    }
  }
}




utils::RotoTransOperator::RotoTransOperator(const Vector3<double> & in,
					    const double theta,
					    const Vector3<double> & iOA
					    ){
  n  = in;
  OA = iOA;
  angle = theta;

  get_rotation_matrix();
}


utils::RotoTransOperator::RotoTransOperator(const Vector3<double> & iT){
  T = iT;
}




MatrixSq3<double> & utils::RotoTransOperator::rotation_matrix(void){ return R; }
Vector3<double>   & utils::RotoTransOperator::rotation_axis(void){ return n; }
double            & utils::RotoTransOperator::rotation_angle(void){ return angle; }
Vector3<double>   & utils::RotoTransOperator::rotation_axispoint(void){ return OA; }
Vector3<double>   & utils::RotoTransOperator::translation(void){ return T; }





void utils::RotoTransOperator::rotate_vector(Vector3<double> & v,
					     Vector3<double> & vrot){
  vrot = R * v;
}


void utils::RotoTransOperator::rotate_vector_startpoint(Vector3<double> & P,
							Vector3<double> & Prot){

  Vector3<double> OB(0), AB(0), AD(0), DB(0), DB_rot(0), OB_rot(0);

  OB = P;

  // Vector3 from point on axis:
  AB = OB - OA;
  // Projection onto axis:
  AD = scalarproduct(AB, n) * n;
  // The rejection part:
  DB = AB - AD;
  // Rotate vector:
  DB_rot = R * DB;
  OB_rot = OA + AD + DB_rot;
  Prot = OB_rot;
}


void utils::RotoTransOperator::rotate(RotoTransOperand & opv,
				      RotoTransOperand & opv_rot){
  Vector3<double> vec  = opv.vector();
  Vector3<double> vecP = opv.vectorpoint();
  
  Vector3<double> vec_rot(0);
  Vector3<double> vecP_rot(0);

  this->rotate_vector(vec, vec_rot);

  this->rotate_vector_startpoint(vecP, vecP_rot);

  opv_rot.vector()      = vec_rot;
  opv_rot.vectorpoint() = vecP_rot;
}




