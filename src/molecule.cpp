#include "molecule.h"

void Molecule::step( double t ) {
  for( unsigned int i=0; i<ckin.size(); i++ )
    ckin[i]->input( ckin_source[i]->kinetics() );
  for( unsigned int i=0; i<cphos.size(); i++ )
    cphos[i]->input( cphos_source[i]->kinetics() );

  sum = 0;
  for( unsigned int i=0; i<stimulus.size(); i++ )
    sum += stimulus_multiplier[i] *stimulus[i]->kinetics(t);
  for( unsigned int i=0; i<ckin.size(); i++ )
    sum += ckin_multiplier[i] *ckin[i]->kinetics();
  for( unsigned int i=0; i<cphos.size(); i++ )
    sum -= cphos_multiplier[i] *cphos[i]->kinetics();
  if( sum < 0 ) sum = 0;
}

void Molecule::addStimulus( double multiplier, Stimulus* new_stimulus )
{
  stimulus_multiplier.push_back( multiplier );
  stimulus.push_back( new_stimulus );
}

void Molecule::addCKin( double multiplier, CKin* new_ckin, Molecule* source )
{
  ckin_multiplier.push_back( multiplier );
  ckin.push_back( new_ckin );
  ckin_source.push_back( source );
}

void Molecule::addCPhos( double multiplier, CPhos* new_cphos, Molecule* source )
{
  cphos_multiplier.push_back( multiplier );
  cphos.push_back( new_cphos );
  cphos_source.push_back( source );
}
