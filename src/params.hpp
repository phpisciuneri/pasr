#ifndef PARAMS_HPP_
#define PARAMS_HPP_

struct Params {
  
  double dt        = 0.0001; // delta time (outer step)
  double tmix      = 0.001;  // mixing time scale
  double tres      = 0.01;   // residency time scale
  double tend      = 0.1;    // time to end simulation
  
  double rtol      = 1.e-7;  // relative error tolerance in direct integration
  double atol      = 1.e-7;  // absolute error tolerance in direct integration
  std::size_t n    = 1;    // number of particles
  int dumpfreq     = 1;      // time-steps to dump result
  
  std::string chem = "lulaw19.cti";
  
};

#endif // PARAMS_HPP_
