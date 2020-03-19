#ifndef PATHCANON_H
#define PATHCANON_H

#include <vector>
using std::vector;

#include <utility>
using std::pair;

#include <string>
using std::string;

#include <fstream>
using std::ofstream;

#include "molecule.h"

class PathCanon
{
  PathCanon(void);
  PathCanon(PathCanon&);

  double alpha, ttotal, deltat;

  vector<Molecule*> mole;

  void read( const string& content );
  struct Content {
    string content;
    int ptr;
    string readNext(void);
    string expect( const string& delimiters );
    vector< pair< string, vector<double> > > parseObj(void);
    vector<string> parseStrings(void);
    bool isEOF(void);
    bool isAlpha( char c );
    bool isNumber( char c );
    bool isSymbol( char c );
    bool isWhiteSpace( char c );
    Content( const string& filecontent ):
      content(filecontent), ptr(0) {}
  };

  // computing step: step through all molecules
  void step(void);

  // print object, stores molecules to print and procedures to print
  struct Print {
    vector<string> mole_name;
    vector<Molecule*> mole;
    ofstream file;
    void init( const string& filename );
    void step( double t );
  };
  Print print;
public:
  // launch simulation, compute, print, terminate
  void run(void);
  // read and initialize simulation
  PathCanon( const string& content, const string& outputfile );
  ~PathCanon(void);
};

#endif
