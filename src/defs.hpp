#ifndef DEFS_HPP_
#define DEFS_HPP_

#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>

// real precision
typedef double real_t;

// random number generation
typedef boost::rand48 rnd_engine_t;
typedef boost::uniform_int< int > rnd_dist_t;
typedef boost::variate_generator< rnd_engine_t, rnd_dist_t > rnd_gen_t;

#endif // DEFS_HPP_
