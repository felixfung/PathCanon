#include "fintegrator.h"
#include "mpi.h"

#include <cmath>

void fIntegrator::input( double newest_history )
{
  history.push_back( newest_history );
}

/******************************************************************************
 * performs fractional integration/derivative
 * this function is the main "meat" of fractional integration
 * where computation takes most time
 * time complexity is O(t) per run of this function
 * where t is the simulation time so far
 * so that the total simulation time complexity over an entire simulation
 * is O(t^2) in this function
 *****************************************************************************/
void fIntegrator::fIntegrate(void)
{
  if( alpha == 0 )
    val = history.back();
  else if( alpha > 0 ) {
/******************************************************************************
 * fractional integration
 *****************************************************************************/
    val = 0;
    int n = history.size();
    int segment_size = (n-1) / Mpi::world_size;
    int segment_begin = Mpi::rank *segment_size +1;
    int segment_end = segment_begin +segment_size;
    if( Mpi::rank == Mpi::world_size -1 )
      segment_end = n;
    #pragma omp parallel for reduction(+:val)
    for( int tau=segment_begin; tau<segment_end; tau++ )
      val += history[tau] *pow( (n-tau)*deltat, alpha-1 );
    Mpi::aggSum( &val );
    val *= deltat/tgamma(alpha);
  }
  else /* alpha < 0 */ {
/******************************************************************************
 * fractional derivation
 *****************************************************************************/
    double a = ceil(alpha) -alpha +1;
    int n = history.size();
    val = 0;
    int N = n-1;
    int segment_size = (N-1) / Mpi::world_size;
    int segment_begin = Mpi::rank *segment_size +1;
    int segment_end = segment_begin +segment_size;
    if( Mpi::rank == Mpi::world_size -1 )
      segment_end = N;
    #pragma omp parallel for reduction(+:val)
    for( int i=segment_begin; i<segment_end; i++ ) {
      double y = history[i+1] -2*history[i] +history[i-1];
      double exp = pow(n-i-1,a+1) +pow(n-i+1,a+1) -2*pow(n-i,a+1);
      val += y *exp;
    }
    Mpi::aggSum( &val );
    val += history[0] *( pow(n,a) *(a+1-n) +pow(n-1,a+1) ) +history[n-1];
    val *= pow(deltat,a) /tgamma(a+2+1) /deltat/deltat;
  }
}

/******************************************************************************
 * constructor and destructor
 *****************************************************************************/
fIntegrator::fIntegrator( double alpha, double deltat )
    : alpha(alpha), deltat(deltat), val(0)
{
  input(0);
}

fIntegrator::~fIntegrator(void)
{
}
