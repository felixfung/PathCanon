#include "ckin.h"
#include <cmath>

void CKin::input( double input )
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
  long double expo_term = expo_sum = input;
  for( unsigned int k=1; k<=fd.size(); k++ ) {
    fd[k-1]->input( expo_term );
    fd[k-1]->fIntegrate();
    expo_term = fd[k-1]->fIntegral();
    expo_sum += pow(-length,k) /tgamma(k+1) *expo_term *2;
  }
  /*long double expo_term2 = input;
  for( unsigned int k=1; k<=fd2.size(); k++ ) {
    fd2[k-1]->input( expo_term2 );
    // double integral
    fd2[k-1]->fIntegrate();
    fd2[k-1]->input( fd2[k-1]->fIntegral() );
    fd2[k-1]->fIntegrate();
    expo_term2 = fd2[k-1]->fIntegral();
    expo_sum += pow(-length,k) /tgamma(k+1) *expo_term2;
  }
  double expo_term3 = input;
  for( unsigned int k=1; k<=fd3.size(); k++ ) {
    fd3[k-1]->input( expo_term3 );
    fd3[k-1]->fIntegrate(); fd3[k-1]->fIntegrate(); fd3[k-1]->fIntegrate();
    expo_term3 = fd3[k-1]->fIntegral();
    expo_sum += pow(-length,k) /tgamma(k+1) *expo_term3;
  }*/
  if( expo_sum<0 ) expo_sum = 0;
}

CKin::CKin( double alpha, double deltat, double length )
    : alpha(alpha), deltat(deltat), expo_sum(0), length(length)
{
  fd.push_back( new fIntegrator(-alpha,deltat) );
/*  fd2.push_back( new fIntegrator(-alpha,deltat) );
  fd3.push_back( new fIntegrator(-alpha,deltat) );*/
  // the number of terms in the exponential series
  // empirically found that this linear relation holds
  // up to around length/deltat <= 50
  int expo_terms = 1.5* length /deltat;
  if( expo_terms < 10 )
    expo_terms = 10;
  for( int k=1; k<expo_terms; k++ )
    fd.push_back( new fIntegrator(-alpha,deltat) );
  /*for( int k=1; k<expo_terms*0.2; k++ )
    fd2.push_back( new fIntegrator(-alpha,deltat) );
  for( int k=1; k<expo_terms; k++ )
    fd3.push_back( new fIntegrator(-alpha,deltat) );*/
}

CKin::~CKin(void)
{
  for( unsigned int i=0; i<fd.size(); i++ )
    delete fd[i];
  /*for( unsigned int i=0; i<fd2.size(); i++ )
    delete fd2[i];
  for( unsigned int i=0; i<fd3.size(); i++ )
    delete fd3[i];*/
}
