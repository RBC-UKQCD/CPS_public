#ifndef COMMANDLINE_H
#define COMMANDLINE_H
#include <config.h>
#include <util/data_types.h>

CPS_START_NAMESPACE

/*!
  simple class to provide access to the command-line arguments
  with type conversion (and checking).
*/
class CommandLine
{
private:

  char* cname;
  int    argc;
  char** argv;
  int   count;

private:

  static CommandLine inst;

  CommandLine():
    cname("CommandLine"),
    argc (0),
    argv (0x0),
    count(0)
  {;}

public:

  /*! register the command-line arguments with the class (always do this)*/
  static void is( int _argc, char** _argv ) 
  { 
    inst.argc = _argc;
    inst.argv = _argv;
  }

  static void reset()      { inst.count=0; } 
  static int    num() { return inst.count; } 

  static char* arg         ( int num );
  static int   arg_as_int  ( int num );
  static int   arg_as_hex  ( int num );
  static Float arg_as_Float( int num ); 


  static char* arg         () ;
  static int   arg_as_int  ();
  static int   arg_as_hex  ();
  static Float arg_as_Float();
};

CPS_END_NAMESPACE
#endif
