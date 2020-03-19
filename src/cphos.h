#ifndef CPHOS_H
#define CPHOS_H

#include "ckin.h"

/******************************************************************************
 * Canonical Phosphatase pathway class
 * given length s, exponent a, and input(t) phi(t),
 * integrate and return p_s(t) = T_s phi(t) = exp(-s(d^a_t)) phi(t)
 *
 * usage: to output the 100 timesteps of a pathway according to stim[t], do
 *   CPhos cphos(alpha,length,deltat);
 *   for( int t=0; t<=100; t++ )
 *     cphos.input( stim[t] );
 *     cout<< cphos.kinetics() <<endl;
 *****************************************************************************/
class CPhos : public CKin
{
  CPhos();
  CPhos(CPhos&);

  fIntegrator fi;
public:
  double kinetics(void) const { return fi.fIntegral(); }
  void input( long double input );
  CPhos( double alpha, double deltat, double length )
      : CKin(alpha,deltat,length), fi(alpha,deltat) {}
  virtual ~CPhos(void) {}
};

#endif
