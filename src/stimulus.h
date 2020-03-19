#ifndef STIMULUS_H
#define STIMULUS_H

#include <string>
using std::string;
#include <vector>
using std::vector;

#include <cmath>

/******************************************************************************
 * Stimulus object
 * abstract class that, given a real valued time
 * return a real valued stimulus
 *
 * also provides a static factory method that,
 * given a stimulus name and parameters, return a Stimulus
 *****************************************************************************/
class Stimulus
{
public:
  virtual double kinetics( double t ) = 0;
  Stimulus(void) {}
  virtual ~Stimulus(void) {}
  static Stimulus* getStimulus
    ( const string& name, const vector<double>& params );
};

class Rect: public Stimulus
{
  double width, period, onset, duration;
public:
  Rect( double width, double period, double onset, double pulses ):
    width(width), period(period), onset(onset), duration( pulses*period ) {}
  virtual double kinetics( double t ) {
    if( ( onset<=t && t<onset+duration )
      && ( 0 <= remainder( t-onset, period )
      && remainder( t-onset,period) < width ) )
        return 1;
    else
      return 0;
  }
};

#endif
