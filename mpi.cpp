#include "mpi.h"

/******************************************************************************
 * initiate static member variables
 * which will be filled up appropriately in Mpi constructor
 *****************************************************************************/
int Mpi::world_size = 1;
int Mpi::rank = 0;
char Mpi::processor_name[] = "";

#ifdef USING_MPI

void Mpi::bcast( string& data )
{
  // first find the file length
  unsigned long length = data.length();
  MPI_Bcast( &length, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD );

  // then put file content into buffer, as null terminated string
  char buffer[length+1];
  if( Mpi::rank == 0 ) {
    for( unsigned long i=0; i<length; i++ )
      buffer[i] = data[i];
    buffer[length] = '\0';
  }

  // broadcast and return
  MPI_Bcast( &buffer[0], length+1, MPI_CHAR, 0, MPI_COMM_WORLD );
  data = buffer;
}

void Mpi::aggSum( double* sum )
{
  double agg;
  MPI_Allreduce( sum, &agg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
  *sum = agg;
}

/******************************************************************************
 * MPI RAII segment
 *****************************************************************************/

// constructor initiates MPI environment and stores environment values
Mpi::Mpi(void)
{
  MPI_Init(NULL,NULL);
  MPI_Comm_size( MPI_COMM_WORLD, &world_size );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  int name_len;
  MPI_Get_processor_name( processor_name, &name_len );
}

// automatically handles finalizing
Mpi::~Mpi(void)
{
  MPI_Finalize();
}

/******************************************************************************
 * if not compiling for MPI, then simply skip all MPI issues
 *****************************************************************************/
#else

string Mpi::bcast( char* data, int length ) { return ""; }
void Mpi::aggSum( double* localsum, double* sum ) {}

Mpi::Mpi(void) {}
Mpi::~Mpi(void) {}

#endif
