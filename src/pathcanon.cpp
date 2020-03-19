#include "pathcanon.h"
#include "mpi.h"

#include <string>
using std::string;

#include <unordered_map>
using std::unordered_map;

#include <fstream>
using std::ofstream;
#include <iomanip>
using std::scientific;
#include <cstdlib>
using std::setw;
#include <iostream>
using std::endl;

#include <stdexcept>
using std::invalid_argument;

/******************************************************************************
 * init, run and step
 *****************************************************************************/
void PathCanon::Print::init( const string& filename )
{
  if( Mpi::rank > 0 ) return;

  file.open( filename.c_str() );
  if( !file )
    throw string("Cannot open file "+string(filename));

  file.precision(14); file<<scientific;

  file << setw(23) <<"Time";
  for( const string& name: mole_name )
    file << setw(23) << name;
  file << endl;
}

void PathCanon::Print::step( double t )
{
  if( Mpi::rank > 0 ) return;
  file << "   " << t;
  for( const Molecule* molecule: mole ){
    file << "   " << molecule->kinetics();
  }
  file << endl;
}

void PathCanon::run(void)
{
  for( double t=0; t<=ttotal; t+=deltat ) {
    for( Molecule* molecule: mole )
      molecule->step(t);
    print.step(t);
  }
}

/******************************************************************************
 * constructor and destructor
 *****************************************************************************/
PathCanon::PathCanon( const string& content, const string& outputfile )
{
  deltat = ttotal = alpha = -1;
  read(content);
  print.init(outputfile);
}

PathCanon::~PathCanon(void)
{
  for( Molecule* del: mole )
    delete del;
}

/******************************************************************************
 * all functions below relate to reading config file, quite tedious
 * the strategy is to first construct all molecules
 * then create pathways and stimuli
 * helper functions exists to read each word in config file
 * in a whitespace independent way, ignoring comments
 * and detects and throws errors
 *****************************************************************************/
void PathCanon::read( const string& filecontent )
{
  Content content(filecontent);

  unordered_map< string, Molecule* > id2mole;
  unordered_map< Molecule*, vector<pair<string,vector<double>>> > mole_params;

  try {
    // read in simulation and molecular descriptions
    while( !content.isEOF() ) {
      string id = content.readNext();
      content.expect(":");
      if( id == "DELTAT" ) {
        if( deltat != -1 )
          throw string("DELTAT definition non-unique");
        double val = stod( content.readNext() );
        content.expect(";") ;
        if( val <= 0 )
          throw string("DELTAT must be positive");
        deltat = val;
      }
      else if( id == "DURATION" ) {
        if( ttotal != -1 )
          throw string("DURATION definition non-unique");
        double val = stod( content.readNext() );
        content.expect(";") ;
        if( val <= 0 )
          throw string("DURATION must be positive");
        ttotal = val;
      }
      else if( id == "ALPHA" ) {
        if( alpha != -1 )
          throw string("ALPHA definition non-unique");
        double val = stod( content.readNext() );
        content.expect(";") ;
        if( val <= 0 )
          throw string("ALPHA must be positive");
        alpha = val;
      }
      else if( id == "PRINT" ) {
        vector<string> mole_name = content.parseStrings();
        for( string name: mole_name ) {
          if( id2mole[name] == NULL )
            throw string("Molecule not found: ")+name;
          print.mole_name.push_back( name );
          print.mole.push_back( id2mole[name] );
        }
      }
      else { // add molecule
        if( id2mole[id] != NULL )
          throw string("Molecule ID already exists: ")+id;
        mole.push_back( new Molecule );
        id2mole[id] = mole.back();
        auto obj_params = content.parseObj();
        mole_params[ mole.back() ] = obj_params;
      }
    }
    // add in pathways
    for( Molecule* molecule: mole ) {
      for( auto source: mole_params[molecule] ) {
        string source_id = source.first;
        double multiplier = source.second[0];
        vector<double> source_params(
          &source.second[1], &source.second.back()+1 );
        Stimulus* stimulus = Stimulus::getStimulus( source_id, source_params );
        if( stimulus != NULL ) {
          molecule->addStimulus( multiplier, stimulus );
        }
        else if( multiplier < 0 ) {
          if( id2mole[ source_id ] == NULL )
            throw string("Molecule "+source_id+" not found");
          double length = source.second[1];
          if( length < 0 )
            throw string("Pathway length must be non-negative: ");
          molecule->addCPhos( -multiplier, new CPhos(alpha,deltat,
            length ), id2mole[ source_id ] );
        }
        else {
          if( id2mole[ source_id ] == NULL )
            throw string("Molecule "+source_id+" not found");
          double length = source.second[1];
          if( length < 0 )
            throw string("Pathway length must be non-negative: ");
          molecule->addCKin( multiplier, new CKin(alpha,deltat,
            length ), id2mole[ source_id ] );
        }
      }
    }
  }
  catch( const string& exception ) {
    if( content.isEOF() )
      throw exception;
    int begin = content.ptr-30; if( begin<0 ) begin = 0;
    int end = content.ptr+30;
    if( end>=(int)filecontent.length() ) end = filecontent.length()-1;
    throw string( exception +string("\n")
      +"Near here:\n"
      +filecontent.substr( begin, content.ptr -begin +1 )
      +string("(*)")
      +filecontent.substr( content.ptr+1, end -content.ptr )
    );
  }
  catch( const invalid_argument& exception ) {
    if( content.isEOF() )
      throw string("Expected number");
    int begin = content.ptr-30; if( begin<0 ) begin = 0;
    int end = content.ptr+30;
    if( end>=(int)filecontent.length() ) end = filecontent.length()-1;
    throw string("Expect number here:\n")
      +filecontent.substr( begin, content.ptr -begin +1 )
      +string("(*)")
      +filecontent.substr( content.ptr+1, end -content.ptr );
  }

  if( deltat == -1 )
    throw string("Time step magnitude must be set via DELTAT");

  if( ttotal == -1 )
    throw string("Simulation duration must be set via DURATION");

  if( alpha == -1 )
    throw string("Pathway integration exponent must be set via ALPHA");

  if( print.mole.empty() )
    throw string("Nothing to print");
}

/******************************************************************************
 * Below are Content functions, used only by PathCanon::read()
 *****************************************************************************/

// read the next word, up to white space or symbol
// if the next character (skipping over white space and comments)
// is not delimiter, throw error
string PathCanon::Content::readNext(void)
{
  bool comment = false;
  while( !isEOF() && (
    isWhiteSpace( content[ptr] )
    || comment || content[ptr]=='|'
  )) {
    if( content[ptr] == '|' )
      comment = !comment;
    ptr++;
  }
  int begin = ptr;
  if( isSymbol( content[ptr] ) )
    ptr++;
  else
    while( !isEOF() && ( isAlpha(content[ptr]) || isNumber(content[ptr]) ) ) {
        //!isSymbol( content[ptr] ) && !isWhiteSpace( content[ptr] ) )
        //if(content[ptr]==' ')
  //std::cout<<"SDLFJKSLDKJFLSKDJ"<<std::endl;
      ptr++;
    }
  string res= content.substr( begin, ptr -begin );
  //std::cout<<begin<<" "<<ptr<<std::endl;
  //std::cout<<res<<std::endl;
  return res;
}

string PathCanon::Content::expect(
  const string& delimiters )
{
  if( isEOF() )
    throw string("Symbols ") +delimiters +string(" expected");

  string next = readNext();
  if( next.length()!=1 )
    throw string("Symbols ") +delimiters +string(" expected");

  bool match = false;
  for( int i=0; i<(int)delimiters.length() && !match; i++ )
    match = next == string(1,delimiters[i]);

  if( !match )
    throw string("Symbols ") +delimiters +string(" expected");

  return next;
}

// read an object specification in the format:
// $*Id( $, $, ..., $ );
// where each $ symbol is a number
// return a pair where the first entry is Id
// and the second is an array of numbers
vector< pair< string, vector<double> > > PathCanon::Content::parseObj(void)
{
  vector< pair< string, vector<double> > > res;
  while(true) {
    string multiplier_str = readNext();
    if( multiplier_str == ";" )
      break;
    double multiplier = stod( multiplier_str );
    expect("*");
    string id = readNext();
    expect("(");
    vector<double> params;
    params.push_back( multiplier );
    string next_param_str;
    while( next_param_str != ")" ) {
      next_param_str = readNext();
      params.push_back( stod( next_param_str ) );
      next_param_str = expect(",)");
    }

    res.push_back( make_pair(id,params) );
  }

  return res;
}

// read in an array of comma separated strings, terminated with semi-colon
vector<string> PathCanon::Content::parseStrings(void)
{
  vector<string> res;
  string delimiter;
  while( delimiter != ";" ) {
    string id = readNext();
    res.push_back(id);
    delimiter = expect(",;");
  }
  return res;
}

bool PathCanon::Content::isEOF(void)
{
  bool comment = false;
  while( ptr < (int)content.length() && (
    isWhiteSpace( content[ptr] )
    || comment || content[ptr]=='|'
  )) {
    if( content[ptr] == '|' )
      comment = !comment;
    ptr++;
  }
return ptr == (int)content.length();
}

bool PathCanon::Content::isAlpha( char c )
{
  return ( 'a'<=c && c<='z' ) || ( 'A'<=c && c<='Z' );
}

bool PathCanon::Content::isNumber( char c )
{
  return ( '0'<=c && c<='9' ) || c=='.' || c=='-' || c=='+';
}

bool PathCanon::Content::isSymbol( char c )
{
  return !isAlpha(c) && !isNumber(c);
}

bool PathCanon::Content::isWhiteSpace( char c )
{
  return c==' ' || c=='\n' || c=='\r' || c=='\t';
}
