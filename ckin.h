#ifndef CKIN_H
#define CKIN_H

#include <vector>
using std::vector;

#include "fintegrator.h"

/******************************************************************************
 * Canonical Kinase pathway class
 * given length s, exponent a, and input(t) phi(t),
 * integrate and return p_s(t) = T_s phi(t) = exp(-s(d^a_t)) phi(t)
 *
 * usage: to output the 100 timesteps of a pathway according to stim[t], do
 *   CKin ckin(alpha,length,deltat);
 *   for( int t=0; t<=100; t++ )
 *     ckin.input( stim[t] );
 *     cout<< ckin.kinetics() <<endl;
 *****************************************************************************/
class CKin
{
  CKin();
  CKin(CKin&);
protected:
  double alpha;             // fractional derivative order
  double deltat;            // integration time step

  vector<fIntegrator*> fd, fd2, fd3;  // array of fractional derivatives,
                            // each element corresponding to one term
                            // in exponential series
  long double expo_sum;          // read only value of the final pathway value
  double length;            // pathway length

public:
  double kinetics(void) const { return expo_sum; }
  void input( long double input ); // takes a stimulus during one timestep
  CKin( double alpha, double deltat, double length );
  virtual ~CKin(void);
};

#endif
