#ifndef MPI_H
#define MPI_H

#ifdef USING_MPI

#include <mpi.h>

#else
#define MPI_MAX_PROCESSOR_NAME 1
#endif

#include <string>
using std::string;

/******************************************************************************
 * wrap MPI in class for RAII
 * this class should be instantiated ONCE within main()
 * so that it will terminate automatically in destructor
 *
 * MPI environment values can be accessed via static members
 *
 * wrapper functions for MPI_SEND and MPI_RECV
 * help elimniate otherwise cluncky preprocessor checking of USING_MPI
 *****************************************************************************/
class Mpi
{
public:

  // static member variables that hold MPI environment values
  static int world_size;
  static int rank;
  static char processor_name[MPI_MAX_PROCESSOR_NAME];

  static void bcast( string& data );
  static void aggSum( long double* localsum );

  Mpi(void);
  ~Mpi(void);
};

#endif
