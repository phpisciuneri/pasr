#include "reactor.hpp"

#include <cantera/IdealGasMix.h>
#include <boost/foreach.hpp>

#include "params.hpp"

Reactor::Reactor( const Cantera::IdealGasMix& chem )
: m_chem(chem), m_phi_mean(chem), m_phi_in(chem), m_p0(Cantera::OneAtm), m_nlocal(0)
{
  
  Params params;

  m_decay = std::exp( -1/params.tmix * params.dt );
  
  m_nspecies = m_chem.nSpecies();

  // Compute outflow
  // # of particles to insert each time step
  m_out_rate = params.dt / params.tres * params.n;
  m_out_residual = 0;
  
  m_reactor.resize( params.n );
  m_reactorNetwork.resize( params.n );

  // Set up solver and arrays
  //m_cpr.solver( m_param.solver );  // ODE solver, set in pasr.in
  //m_cpr.reltol( m_param.rtol );    // ODE relative tolerance
  //m_cpr.abstol( m_param.atol );    // ODE absolute tolerance

}

void Reactor::initialize( std::size_t n)
{

  // set number of local particles
  m_nlocal = n;
  
  // INLET BOUNDARY CONDITION
  m_phi_in.setState_TPX(300.0, Cantera::OneAtm, "CH4:0.5, O2:1.0, N2:3.76");

  // INITIAL CONDITION
  /*
   Here we take the reaction defined in the inlet boundary condition:
     .5CH4 + O2 + 3.76 --> Products
   and we allow it to proceed until it reaches equilibrium
   i.e. until all the fuel or oxidizer is completely consumed
   and the reaction cannot proceed any further
   */
  Cantera::IdealGasMix phi_temp(m_chem);
  phi_temp.setState_TPX(300.0, Cantera::OneAtm, "CH4:0.5, O2:1.0, N2:3.76");
  phi_temp.equilibrate("HP");
  // at this stage:
  //   - the composition, X0, is modified
  //   - the pressure, p0, is unchanged (constant pressure calculation)
  //   - the temperature, T0,  has increased (exothermic reaction)
  
  // CONSTRUCT PARTICLE ARRAY
  // - composition of each particle is the mass fractions and temperature
  //   of the initial condition
  for (std::size_t i=0; i<m_nlocal; i++) {
    m_reactor[i].insert( phi_temp );
    m_reactorNetwork[i].reset( new Cantera::ReactorNet() );
    m_reactorNetwork[i]->addReactor( m_reactor[i] );
  }
  
  std::cout << phi_temp.report() << std::endl;
  std::cout << m_reactor[0].contents().report() << std::endl;

}

void Reactor::react()
{
  
  Params params;

  // rxn substep
  for (std::size_t i=0; i<m_reactorNetwork.size(); i++) {
    //std::cout << m_reactorNetwork[i]->reactor(0).contents().report() << std::endl;
    m_reactorNetwork[i]->setInitialTime(0);
    m_reactorNetwork[i]->advance(params.dt);
    //std::cout << m_reactorNetwork[i]->reactor(0).contents().report() << std::endl;
  }
 
}

void Reactor::mix()
{
  
  std::vector<doublereal> newMassFractions(m_nspecies);
  doublereal newTemperature;

  // mixing substep
  for (std::size_t i=0; i<m_reactor.size(); i++) {
    
    for (std::size_t k=0; k<m_nspecies; k++)
      newMassFractions[k] = m_decay*( m_reactor[i].massFraction(k) - m_phi_mean.massFraction(k) ) + m_phi_mean.massFraction(k);
    
    newTemperature = m_decay*( m_reactor[i].temperature() - m_phi_mean.temperature() ) + m_phi_mean.temperature();
    
    Cantera::ThermoPhase c = m_reactor[i].contents();
    std::cout << c.report() << std::endl;
    c.setState_TPY( newTemperature, m_p0, newMassFractions.data() );
    m_reactor[i].syncState();
  }
}

size_t Reactor::n_global_exchange()
{

  // inlet/outlet
  m_out_residual += m_out_rate;
  m_Nin = static_cast<size_t>( std::floor( m_out_residual ) );
  m_out_residual -= m_Nin; // save fractional particles for next time
  return m_Nin;

}

void Reactor::bc( int n_local_exchange )
{

  // initialize generator for choosing particles at random
  static rnd_gen_t rand( rnd_engine_t( rnd_engine_t::result_type(0) ), 
    rnd_dist_t( 0, m_nlocal-1 ) );

  // randomly replace particles with inlet boundary particles
  for(int i=0; i!=n_local_exchange; ++i)
    m_reactor[ rand() ].insert( m_phi_in );

}

void Reactor::update_mean()
{

  std::vector<doublereal> meanMassFractions(m_nspecies,0);
  doublereal meanTemperature = 0;
  std::size_t n = m_reactor.size();
 
  for (std::size_t i=0; i<n; i++ ) {
    for(std::size_t k=0; k<m_nspecies; k++)
      meanMassFractions[k] += m_reactor[i].massFraction(k) / n;
    meanTemperature += m_reactor[i].temperature() / n;
  }
  
  m_phi_mean.setState_TPY( meanTemperature, Cantera::OneAtm, meanMassFractions.data());
  
  std::cout << m_phi_mean.report() << std::endl;

}

void Reactor::output( real_t t )
{

  static bool first_call = true;

  // create output file & write header
  if ( first_call ) 
  {
    m_output.open("pasr.out", std::ofstream::out );
    m_output << "t "; 
    for (std::size_t i=0; i!=m_nspecies; ++i)
      m_output << m_chem.species(i) << " ";
    //m_output << " T Nin ";
	m_output << " T Nin" << std::endl;
    /*
	for(int i=0; i!=m_chem.nspec(); ++i) 
      m_output << "var_" << m_chem.species(i) << " ";
    m_output << " varT\n";
	*/

    first_call = false;
  }

  m_output << t << " ";
  for (std::size_t i=0; i<m_nspecies; i++)
    m_output << m_phi_mean.massFraction(i) << " ";
  //m_output << m_Nin << " ";
  m_output << m_Nin << std::endl;
  
  // compute variance
  /*
  cpr_t::phi var(m_chem); 
  BOOST_FOREACH( cpr_t::phi& p, m_particles )
    for(int i=0; i!=p.size(); ++i) 
      var[i] += ( p[i] - m_phi_mean[i] )*( p[i] - m_phi_mean[i] );

  for(int i=0; i!=m_phi_mean.size(); ++i) 
    var[i] = std::sqrt( var[i] ) / m_particles.size();

  for(int i=0; i!=m_phi_mean.size(); ++i)
    m_output << var[i] << " ";
  m_output << std::endl;
  */

}
