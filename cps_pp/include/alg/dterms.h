//------------------------------------------------------------------
//  
//   dterms.h
//  
//   Header file for the  DTerms class. A class holding the derivative
//   terms acting only on quarks. This is used to build all other derivative
//   terms by appropriate communication
//
//   Based on Chulwoo's optimization proposal.
//  
//   April 2002
//  
//   Kostas Orginos
//  
//------------------------------------------------------------------*/
#ifndef INCLUDED_DTERMS_H
#define INCLUDED_DTERMS_H
#include <stdio.h>
#include <comms/scu.h>
#include <util/data_types.h>
#include <util/site.h>

CPS_START_NAMESPACE


/*!
  This class has a member array holding \f$2^N\f$ Complex numbers
  per lattice Site. N is the number of derivatives the operator has.
  All derivatives are taken to act on the quark i.e.

 */
class DTerms
{
  Complex *term ;
  int  num     ; // Number of terms
  char *cname ;

  void Allocate() ; 
  void Delete() ;
  
public:
  DTerms() ; 
  DTerms(int N) ; 

  Complex& operator()(int i, int site)
  {
    if((i>=num)||(i<0)||(site<0)||(site>=GJP.VolNodeSites()))
      {
	ERR.General(cname, "()", "out of range\n");
	return *term ;
      }
    return *(term + i + num*site) ;
  }
  
 Complex& GetTerm(const int d, const int *vec, Complex& tmp) const;

  ~DTerms(){Delete();}
  
};

CPS_END_NAMESPACE

#endif // !INCLUDED_DTERMS_H  


