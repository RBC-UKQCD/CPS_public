//
//
// Simple XML writer class for the testing framework
// This writes out XML to a file.
// There is no check that the XML is 'well formed'.
// (Use the xmllint tool on a linux box to do this)
// 
// The array data structures are specific to the
// xmldiff testing tool.
//

#include <stdio.h>
#include <stdlib.h>	// exit()
#include <math.h>
#include <string.h>

#include <config.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/verbose.h>
#include <util/error.h>
#include <util/error.h>
#include <comms/glb.h>
#include <alg/w_ginfo.h>

#include <util/dump_xml.h>

CPS_START_NAMESPACE


// open at the contructor
// write header information
//

dump_xml::dump_xml(char filename[], char top_tag[])
{
  // open file 

  fp = fopen(filename, "w");
  if( fp == NULL ) {
    ERR.General("","dump_xml","Could not open file %s!\n", filename);
  }



  // write xml header 
  fprintf(fp,"<?xml version=\"1.0\"?>\n") ; 
  fprintf(fp,"<%s>\n",top_tag) ; 

  // store the tag
  int len = 1 + strlen(top_tag) ; 
  main_tag = (char *)  smalloc(len * sizeof(char) )  ;
  strcpy(main_tag,top_tag ) ; 
    
}

// write int 
void dump_xml::write(int x, char tagname[])
{

  fprintf(fp,"<%s>%d</%s>\n",
	 tagname,x,tagname) ; 

}


// write Float 
void dump_xml::write(Float x, char tagname[])
{

  fprintf(fp,"<%s>%e</%s>\n",
	 tagname,x,tagname) ; 

}


//
// write Float_array
//

void dump_xml::write(Float x[], int dim, char tagname[])
{

  fprintf(fp,"<%s>\n",tagname) ; 

  for(int i = 0 ; i < dim ; ++i)
    fprintf(fp,"%e ",x[i]);

  fprintf(fp,"</%s>\n",tagname) ; 

}




// close method 

void dump_xml::close()
{
 fprintf(fp,"</%s>\n",main_tag) ; 
  fclose(fp) ; 
}



dump_xml::~dump_xml()
{

  sfree(main_tag) ; 
}




CPS_END_NAMESPACE
