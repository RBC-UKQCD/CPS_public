#include<config.h>
CPS_START_NAMESPACE
/*
 * class WspectField
 * 
 * calculate and store B and E fields
 */

#ifndef INCLUDED_WSPECT_FIELD_H
#define INCLUDED_WSPECT_FIELD_H

CPS_END_NAMESPACE
#include <util/data_types.h>        // Float, Complex
#include <alg/w_ginfo.h>
CPS_START_NAMESPACE

class Lattice; 

class WspectField : public WspectGinfo{
 
 public:

  WspectField(const Lattice &lat, const WspectFuzzing* fuzz_p);
  ~WspectField();
  void run();
 
  //ACCESSOR
  Matrix getField(int site[],FieldTensorId ftId);
  

 private:

  static char *  d_class_name;
  
  //Buffer used to store Field Tensors
  int ft_size; //size of buffer in IFloats
  Matrix *ft_p; //index [site][FieldTensorId]
 
  //constant references
  const Lattice & d_lat; //for getLink
  const WspectFuzzing *fuzz_p;

  //calculate B and E fields from 
  //original gauge links : fuzz_p=0
  //or fuzzed links      : fuzz_p!=0
  void calcField();
  Matrix calcLclField(int site[], int mu, int nu);
  
 
  //calculate plaq from links(original or fuzzed)
  Matrix calcPlaq(int site[], int musign, int mu, int nusign, int nu); 
  
  //
  Matrix getLink(int x[],int musign,int mu);
};

#endif 

CPS_END_NAMESPACE
