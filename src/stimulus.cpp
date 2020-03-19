#include "stimulus.h"

#include <string>
using std::string;
#include <vector>
using std::vector;

Stimulus* Stimulus::getStimulus
  ( const string& name, const vector<double>& params )
{
  if( name == "Rect" ) {
    if( params.size() != 4 )
      throw string("Rectangular stimulus requires 4 parameters");
    return new Rect( params[0], params[1], params[2], params[3] );
  }
  return NULL;
}
