#ifndef FINTEGRATOR_H
#define FINTEGRATOR_H

#include <vector>
using std::vector;

/******************************************************************************
 * Caputo fractional integration
 * solves for (d^-a/dt) phi(t), given phi(t)
 * when a is -ve, calculate fractional derivative
 * when a is 0, return phi(t)
 *
 * usage: to output the 100 timesteps of a pathway according to stim[t], do
 *   fIntegrator fi(alpha,deltat);
 *   for( int t=0; t<=100; t++ )
 *     fi.input( stim[t] );
 *     fi.fIntegrate();
 *     cout<< fd.fIntegral() <<endl;
 *****************************************************************************/
class fIntegrator
{
  fIntegrator(void);
  fIntegrator(fIntegrator&);

  double alpha;            // fractional derivative order
  double deltat;           // integration time step
  double val;              // read only integral value
  vector<double> history;  // fractional dervative is history dependent
public:
  // provide phi in eqn, x = (d^a/dt) phi
  void input( double newest_history );
  // integrate input
  void fIntegrate(void);

  fIntegrator( double alpha, double deltat );
  virtual ~fIntegrator(void);
  double fIntegral() const { return val; }
};

#endif
