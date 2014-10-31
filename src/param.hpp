


#ifndef PARAM_HPP
#define PARAM_HPP


#include "utils.hpp"
#include "utils-math.hpp"
#include "utils-vector.hpp"


/* ##############################################################################

   --------------------------------------------------------
              |                   Xlimuse
   --------------------------------------------------------
              | true                  |    false
   --------------------------------------------------------
              | Xmin < Xmax => free   |
   X:         |                       |   always free
              | Xmin = Xmax => fixed* |
   --------------------------------------------------------
	      | X is constrained      | X is unconstrained
   --------------------------------------------------------
   *) Special case: If Xmin = 0 = Xmax => X is free (unconstrained).

   Note: Barrier/penalty functions are only used with constrained parameters.

   ############################################################################## */

enum parametertype { PARAM_FREE_WITH_LIMITS, PARAM_FIXED, PARAM_FREE };


void set_param_type(double xmin, double xmax, parametertype & xtype);



  
class Param {

private:
  utils::Vector<double> mX;
  utils::Vector<parametertype> mXtype;
  utils::Vector<double> mXmin;
  utils::Vector<double> mXmax;

public:
  // Default constructor
  Param();
  // Default destructor
  virtual ~Param(){}

  // Constructor
  Param(const utils::Vector<double> & xi);

  // Copy constructor
  Param(const Param & sv);

  // Operators
  Param & operator=(const Param & sv);



  // **************************************************************************
  // **************************************************************************
  // **
  // **  Methods that may be redefined in derived classes:
  // **
  // **************************************************************************
  // **************************************************************************
  // obj.X(0) returns mX[0] for reading as well as writing to, so this construction
  // enables obj.X(0) = val


  /* IO parameters are all input/output parameters visible to user. Some of these
     parameters are constrained, others are free. Only the free parameters must
     be used in any optimization routine. 
  */

  virtual void init(const int N);

  virtual void init(const utils::Vector<double> & xi);

  virtual void init(const utils::Vector<double> & xi,
		    const utils::Vector<double> & ximin,
		    const utils::Vector<double> & ximax);

  virtual void init(const utils::Vector<double> & xi,
		    const utils::Vector<parametertype> & xitype,
		    const utils::Vector<double> & ximin,
		    const utils::Vector<double> & ximax);


  // ****************************************************************************
  // ****************************************************************************
  // All parameters have lower and upper limits.
  // Lack of any suffix such as 'all' or 'free' for a vector means that it
  // contains all parameters.

  // Set/get any parameter
  virtual double & X(int i);
  // Set/get all parameters
  virtual utils::Vector<double> & X(void);
 
  // Set/get min parameter limit
  virtual double & Xmin(int i);
  // Set/get all min parameter limits
  virtual utils::Vector<double> & Xmin(void);

  // Set/get max parameter limit
  virtual double & Xmax(int i);
  // Set/get all max parameter limits
  virtual utils::Vector<double> & Xmax(void);

  // Set/get parameter type
  virtual parametertype & Xtype(int i);
  // Set/get all parameter types
  virtual utils::Vector<parametertype> & Xtype(void);


  // ****************************************************************************
  // ****************************************************************************
  // Utilities

  virtual bool X_is_free(const int i);
  virtual bool X_is_fixed(const int i);
  virtual bool X_is_unconstrained(const int i);
  
  virtual int NXfree(void);





  // ****************************************************************************
  // ****************************************************************************
  // Given any vector (of correct type of course), update all the parameters
  // (the full parameter vector).

  virtual utils::Vector<double> Xupdate(const utils::Vector<double> & xi);

  // ##################################################################################


  // Return a vector with all parameters:
  virtual utils::Vector<double> X(const utils::Vector<double> & xi);

  // Return a vector with the free parameters:
  virtual utils::Vector<double> Xfree(void);


  // Return a vector with the free parameters:
  virtual utils::Vector<double> Xfree(const utils::Vector<double> & xi);


  // ##################################################################################

  // Inquire if Xfree[] is inside the limits.
  virtual bool Xfree_is_good(void);


  // Inquire if Xfree[] is inside the limits.
  virtual bool Xfree_is_good(const utils::Vector<double> & xi);


  // ##################################################################################


  // There may be vectors A that have the same size as X (the full vector of all
  // parameters). At some places the elements of A corresponding to the free
  // parameters may be needed. To simplify operations, a method that converts A
  // to a (smaller) vector where elements A[i] corresponds to free parameter Xfree[i]
  // is needed.

  virtual utils::Vector<double> map_any_as_Xfree(const utils::Vector<double> & any);
  virtual utils::Vector<parametertype> map_any_as_Xfree(const utils::Vector<parametertype> & any);


  // ##################################################################################



} ;





#endif

