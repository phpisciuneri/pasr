project( PASR CXX C )
cmake_minimum_required( VERSION 3.3 )

set( CMAKE_CXX_STANDARD 11 )

find_package( LAPACK )

include_directories(
  /usr/local/Cellar/cantera/2.3.0/include
  /usr/local/Cellar/boost/1.64.0_1/include
  /usr/local/Cellar/fmt/3.0.1/include
)

link_directories(
  /usr/local/Cellar/cantera/2.3.0/lib
  /usr/local/Cellar/fmt/3.0.1/lib
)

add_executable( pasr
                src/defs.hpp
		src/params.hpp
		src/pasr.cpp
		src/reactor.cpp
		src/reactor.hpp )

target_link_libraries( pasr
                       cantera
		       fmt
		       ${BLAS_Accelerate_LIBRARY}
		       ${LAPACK_Accelerate_LIBRARY} )
