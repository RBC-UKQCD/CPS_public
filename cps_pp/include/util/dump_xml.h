#ifndef DUMP_XML_INC
#define DUMP_XML_INC 1

#include<config.h>
#include <util/data_types.h>

//
//  A simple calls to write xml to a file
//
//

CPS_START_NAMESPACE

class dump_xml 
{
public:

  dump_xml(char filename[], char top_tag[]) ; 

  // methods to write out the 
  void write(int x, char tagname[]) ; 
  void write(Float x, char tagname[]) ;
  void write(Float x[], int dim, char tagname[]) ;

  void close() ; 

  ~dump_xml(); 

private:

  FILE * fp ; 
  char *main_tag ;

} ;


CPS_END_NAMESPACE
#endif
