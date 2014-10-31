
#ifndef MTWISTER_HPP
#define MTWISTER_HPP


// http://stackoverflow.com/questions/3617317/encapsulating-boostrandom-for-ease-of-usage-to-replace-rand
// http://tudat.tudelft.nl/projects/tudat/wiki/Boost_random_number_generator
// and standard c++ sources (boost homepage not too helpful!)

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
//#include <boost/generator_iterator.hpp>



class rand_mtwister {
  typedef boost::mt19937                      ENG;
  typedef boost::uniform_real<double>         DIST;
  typedef boost::variate_generator<ENG&,DIST> GEN;

private:
  int mseed;
  ENG    mgenerator;   // Generator type
  DIST * mp_uni_dist;  // Distribution
  GEN  * mp_rand_unif; // Function which provides a new random number in a
                       // sequence (initialized by the seed) each time it's called
  double mgauss_ave;
  double mgauss_dev;
  double mgauss_r1;
  double mgauss_r2;
  bool   mgauss_empty;

  // *******************************************************************
  // Inhibit copy construction and assignment by declaring them private.
  // Missing implementation (definition) will cause compilation error
  // (by design!) if some function tries to call them.
  // *******************************************************************
  rand_mtwister(const rand_mtwister & sv);
  rand_mtwister & operator=(const rand_mtwister & sv);


public:
  rand_mtwister();
  rand_mtwister(int s);
  ~rand_mtwister();

  int seed();
  void seed(int);

  double unif();
  double gauss(double ave=0, double dev=1);
} ;




#endif

