#ifndef REACTOR_HPP_
#define REACTOR_HPP_

#include <fstream>
#include <vector>
#include <cantera/zerodim.h>
#include <cantera/IdealGasMix.h>

#include "defs.hpp"

class Reactor
{

public:

  Reactor( const Cantera::IdealGasMix& chem );

  void initialize( size_t n );
  void react();
  void mix();
  size_t n_global_exchange();
  void bc( int n_local_exchange );
  void update_mean();
  void output( real_t t );

private:

  Cantera::IdealGasMix m_chem;         // thermochemistry object
  std::vector<Cantera::IdealGasConstPressureReactor> m_reactor; // cpreactor object
  std::vector< std::unique_ptr<Cantera::ReactorNet> > m_reactorNetwork;
  Cantera::IdealGasMix m_phi_mean;     // mean particle composition
  Cantera::IdealGasMix m_phi_in;       // inlet particle composition
  std::ofstream m_output;              // output file handle
  doublereal m_p0;                     // ambient pressure
  real_t m_decay;                      // mixing time-scale
  real_t m_out_rate;                   // rate of incoming/out-going partciles
  real_t m_out_residual;               // fractional particles
  size_t m_Nin;                        // # particles to insert / time step
  size_t m_nlocal;                     // # of particles local to this rank
  size_t m_nspecies;

};

#endif // REACTOR_HPP_
