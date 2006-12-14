#include <config.h>
#include <stdio.h>
#include <util/error.h>
#include <util/command_line.h>

CPS_START_NAMESPACE

CommandLine CommandLine::inst;

int CommandLine::arg_as_int()
{ 
  const int ret(arg_as_int(inst.count));
  inst.count++;
  return ret;
}

int CommandLine::arg_as_hex()
{ 
  const int ret(arg_as_hex(inst.count));
  inst.count++;
  return ret;
}

Float CommandLine::arg_as_Float() 
{
  const Float ret(arg_as_Float(inst.count));
  inst.count++;
  return ret;
}

char* CommandLine::arg() 
{
  char* ret(arg(inst.count));
  inst.count++;
  return ret;
}

char* CommandLine::arg( int num )
{
  if ( num >= (inst.argc-1) || num < 0 )
    {
      ERR.General(inst.cname,"arg()","out of range\n");
    }
  return inst.argv[num+1];
}


Float CommandLine::arg_as_Float( int num )
{
  const char* fname("arg_as_Float");
  if ( num >= (inst.argc-1) || num < 0 )
    {
      ERR.General(inst.cname,fname,"out of range\n");
    }
  double tmp_flt; 
  if ( sscanf(inst.argv[num+1],"%lg",&tmp_flt) != 1 )
    {
      ERR.General(inst.cname,fname, "bad conversion; command-line arg %i\n", num) ;
    }
  return tmp_flt;
}


int CommandLine::arg_as_int( int num )
{
  const char* fname("arg_as_int");
  if ( num >= inst.argc-1 || num < 0 )
    {
      ERR.General(inst.cname,fname,"out of range\n");
    }
  int tmp_int; 
  if ( sscanf(inst.argv[num+1],"%i",&tmp_int) != 1 )
    {
      ERR.General(inst.cname,fname, "bad conversion; command-line arg %i\n", num) ;
    }
  return tmp_int;
}

int CommandLine::arg_as_hex( int num )
{
  const char* fname("arg_as_hex");
  if ( num >= inst.argc-1 || num < 0 )
    {
      ERR.General(inst.cname,fname,"out of range\n");
    }
  int tmp_int; 
  if ( sscanf(inst.argv[num+1],"%x",&tmp_int) != 1 )
    {
      ERR.General(inst.cname,fname, "bad conversion; command-line arg %i\n", num) ;
    }
  return tmp_int;
}

CPS_END_NAMESPACE
