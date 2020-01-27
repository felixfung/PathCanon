#ifndef MOLECULE_H
#define MOLECULE_H

#include <vector>
using std::vector;

#include "ckin.h"
#include "cphos.h"
#include "stimulus.h"

/******************************************************************************
 * Molecule object
 * encapsulates the idea of a molecule with a non-negative kinetics value
 * by taking in contributions from kinase, phosphatase, stimuli
 *
 * usage:
 *
 * Molecule molecule;
 * molecule.addStimulus( some_multiplier, new_stimulus );
 * molecule.addCKin( some_multiplier, new_ckin );
 * molecule.addCPhos( some_multiplier, new_cphos );
 * for( double t=0; t<=ttotal; t+=deltat )
 *   molecule.step(t);
 *****************************************************************************/
class Molecule
{
  vector<double> stimulus_multiplier;
  vector<Stimulus*> stimulus;
  vector<double> ckin_multiplier;
  vector<CKin*> ckin;
  vector<Molecule*> ckin_source;
  vector<double> cphos_multiplier;
  vector<CPhos*> cphos;
  vector<Molecule*> cphos_source;
  double sum;
public:
  double kinetics(void) const { return sum; }
  void step( double t );
  void addStimulus( double multiplier, Stimulus* new_stimulus );
  void addCKin( double multiplier, CKin* new_ckin, Molecule* source );
  void addCPhos( double multiplier, CPhos* new_cphos, Molecule* source );
  Molecule(void): sum(0) {}
  ~Molecule(void) {}
};

#endif
