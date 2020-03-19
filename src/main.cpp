#include <string>
using std::string;
#include <fstream>
using std::ifstream;
using std::ios;
#include <cerrno>

#include <utility>
using std::pair;
using std::make_pair;

#include <iostream>
using std::cerr;
using std::endl;

#include "mpi.h"
#include "pathcanon.h"

pair<string,string> choosePathFile( int argc, char* argv[] );
string readWholeFile( const string& pathfile );

/******************************************************************************
 * main function
 *****************************************************************************/
int main( int argc, char* argv[] )
{

  // initiate MPI environment (shall be automatically finalized)
  Mpi mpi; 
  // master process reads in config file
  string content, outputfile;
  if( Mpi::rank == 0 ) {
    try {
      auto files = choosePathFile( argc, argv );
      string pathfile = files.first;
      outputfile = files.second;
      content = readWholeFile( pathfile );
    }
    catch( const string& exception ) {
      cerr << exception << endl;
      exit(-1);
    }
  }

  try {
    // every process has the entire solution object
    Mpi::bcast( content );
    PathCanon pathcanon( content, outputfile );
    // launch simulation
    // processes and threads synchronize to solve fIntegrator
    pathcanon.run();
  }
  catch( const string& exception ) {
    cerr << exception << endl;
    exit(-1);
  }

  return 0;
}

/******************************************************************************
 * given the argument list,
 * return a pair of file paths
 * the first being the path config file (default canon.path)
 * the second being the output file (default canon.output)
 *****************************************************************************/
pair<string,string> choosePathFile( int argc, char* argv[] )
{
  if( argc == 1 )
    return make_pair("canon.path","canon.output");
  else if( argc == 2 )
    return make_pair( argv[1], "canon.output" );
  else if( argc == 3 )
    return make_pair( argv[1], argv[2] );
  else
    throw string("PathCanon requires at most two arguments: ")
          +string("path file and output file names");
  return make_pair("","");
}

string readWholeFile( const string& pathfile )
{
  ifstream in( pathfile, ios::in | ios::binary );
  if( in ) {
    string content;
    in.seekg( 0, ios::end );
    content.resize( in.tellg() );
    in.seekg( 0, ios::beg );
    in.read( &content[0], content.size() );
    in.close();
    return content;
  }
  throw(errno);
}
