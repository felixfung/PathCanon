#include "ckin.h"
#include <cmath>

void CKin::input( long double input )
{
  // ensure non-negative input
  if( input<0 ) input = 0;

  // following is an efficient way to calculate exponential series
  // of fractional derivatives
  // the terms are alternating in sign
  // where higher order terms tend to be bigger in magnitude
  // than lower order terms
  // so that the series convergence is very tricky:
  // too few terms will result in non-convergence
  // too many terms will result in nan
  // because the term magnitude exceeds the double capacity
  // this is especially true when the length gets long
  // the solution is to use the Lie semigroup property
  // to chain multiple CKin's together
  // at some maximal length for each CKin
  expo_sum = input;
  fd[0]->input( input );
  fd[0]->fIntegrate();
  expo_sum += pow(-length,1) /tgamma(2) *fd[0]->fIntegral();
  for( unsigned int i=1; i<fd.size(); i++ ) {
    fd[i]->input( fd[i-1]->fIntegral() );
    fd[i]->fIntegrate();
    expo_sum += pow(-length,i+1) /tgamma(i+2) *fd[i]->fIntegral();
  }
  if( expo_sum<0 ) expo_sum = 0;
}

CKin::CKin( double alpha, double deltat, double length )
    : alpha(alpha), deltat(deltat), expo_sum(0), length(length)
{
  fd.push_back( new fIntegrator(-alpha,deltat) );
  int expo_terms = 1.5* length /deltat;
  if( expo_terms < 10 )
    expo_terms = 10;
  for( int k=1; k<expo_terms; k++ )
    fd.push_back( new fIntegrator(-alpha,deltat) );
}

CKin::~CKin(void)
{
  for( unsigned int i=0; i<fd.size(); i++ )
    delete fd[i];
}
