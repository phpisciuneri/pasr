#include <boost/progress.hpp>
#include <vector>
#include <ostream>

#include "defs.hpp"
#include "params.hpp"
#include "reactor.hpp"

// +------+
// | PASR |
// +------+
int main()
{

  int size = 1;
  int rank = 0;
  Params params;

  // create reactor
  Cantera::IdealGasMix chem("lulaw19.cti");
  Reactor reactor( chem );

  // initialize reactor
  reactor.initialize( params.n );

  // initialize mean composition
  reactor.update_mean();

  // init progress meter
  std::ofstream null_stream;
  std::ostream* os = ( rank == 0 ) ? &std::cout : &null_stream;
  boost::progress_display prog( unsigned( params.tend / params.dt ), *os );

  // TIME LOOP 
  boost::progress_timer tmr;
  size_t iter = 0;
  std::vector< int >  rank_nexchange( size );
  static rnd_gen_t rand( rnd_engine_t( rnd_engine_t::result_type(0) ), 
    rnd_dist_t( 0, size-1 ) );

  // begin timing
  for (real_t t=0; t<params.tend; t+=params.dt, ++prog, iter++)
  {

    // mixing substep
    reactor.mix();

    // rxn substep
    reactor.react();

    // particle exchange
    int Ninout = reactor.n_global_exchange();
    
    // boundary conditions
    reactor.bc( Ninout );

    // update the particle mean
    // chemical reaction and input/output have changed it
    reactor.update_mean();

    // output stuff
    if ( rank == 0 && iter % params.dumpfreq == 0 )
      reactor.output( t );

  } // time loop

  return 0;

}

/*
// +------------------+
// | SANITY DASHBOARD |
// +------------------+
std::cout << std::endl;
std::cout.width(7);  std::cout << "";
std::cout.width(23); std::cout << "Inlet Composition"; 
std::cout.width(25); std::cout << "Initial Composition";
std::cout.width(25); std::cout << "Mean Composition (t=0)" << std::endl;

// species
for (int i=0; i<chem.nspec(); i++) {
std::cout.width(7);  std::cout << chem.species()[i];
std::cout.width(23); std::cout << phi_in.Y()[i]; 
std::cout.width(25); std::cout << some_phi.Y()[i];
std::cout.width(25); std::cout << particle_mean[i] << std::endl;
}
// temperature
std::cout.width(7);  std::cout << "T (K)";
std::cout.width(23); std::cout << phi_in.T(); 
std::cout.width(25); std::cout << some_phi.T();
std::cout.width(25); std::cout << particle_mean.T() << std::endl;
std::cout << std::endl;
std::cout << "Operating pressure: " << p0 << std::endl;
std::cout << "Total number of particles: " << prm.n << std::endl;
std::cout << "Particles input per time step: " << out_rate << std::endl;
std::cout << "*** All species values are mass fractions. ***" << std::endl;
std::cout << std::endl;
*/
