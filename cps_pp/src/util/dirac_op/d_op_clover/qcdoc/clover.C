#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2013-04-05 17:51:13 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_clover/qcdoc/clover.C,v 1.8 2013-04-05 17:51:13 chulwoo Exp $
//  $Id: clover.C,v 1.8 2013-04-05 17:51:13 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: clover.C,v $
//  $Revision: 1.8 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_clover/qcdoc/clover.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------

CPS_END_NAMESPACE
#include<util/clover.h>
CPS_START_NAMESPACE

//---------------------------------------------------------------------------
// clover.C:
//         implement struct Clover.
//          
// WARNING:                                                                 
//                                                                          
// This set of routines will work only if the node sublattices have         
// even number of sites in each direction.                                  
//---------------------------------------------------------------------------

CPS_END_NAMESPACE
#include<util/gjp.h>
#include<util/smalloc.h>
#include<util/error.h>
#include<util/verbose.h>
CPS_START_NAMESPACE


//----------------------------------------------------------------------
// external variable definition for clover_mat_mlt_asm.asm
//----------------------------------------------------------------------
CPS_END_NAMESPACE
//#include<comms/nga_reg.h>
CPS_START_NAMESPACE
//extern const unsigned int clover_cram_scratch_addr = CRAM_SCRATCH_ADDR;


//---------------------------------------------------------------------------
// clover_init()  
//---------------------------------------------------------------------------
// Purpose:
//  performs all initializations needed before clover funcs     
//  are called. It sets the addressing related arrays and reserves memory 
//  for the needed temporary buffers. It only needs to be called           
//  once at the begining of the program (or after a clover_end call)       
//  before any number of calls to clover funcs are made.
//---------------------------------------------------------------------------
void clover_init(Clover *clover_p)  
{
  char *cname = "";
  char *fname = "clover_init";  
  VRB.Func(cname,fname);

  // set the local lattice size
  //------------------------------------------------------------------------- 
  clover_p->nsites[0] = GJP.XnodeSites();
  clover_p->nsites[1] = GJP.YnodeSites();
  clover_p->nsites[2] = GJP.ZnodeSites();
  clover_p->nsites[3] = GJP.TnodeSites();
 

  // Do initializations before the wilson library can be used      
  // Initialization involve memory allocation.                       
  //-------------------------------------------------------------------------
  static Wilson wilson_struct;
  clover_p->wilson_p = &wilson_struct;
  wilson_init(clover_p->wilson_p);

  // allocate temporary buffers to be used by DiracOpClover::MatPcDagMatPc
  // and etc.
  //-------------------------------------------------------------------------
  int f_size = SPINOR_SIZE * GJP.VolNodeSites() * sizeof(IFloat);  
  clover_p->frm_buf0 = (IFloat *)smalloc(f_size);
  clover_p->frm_buf1 = clover_p->frm_buf0 + f_size / (2 * sizeof(IFloat));
  bzero((char *)clover_p->frm_buf0,f_size);

  if (clover_p->frm_buf0 == 0 || clover_p->frm_buf1 == 0)
    ERR.Pointer(cname,fname, "frm_buf");
  VRB.Smalloc(cname,fname, "frm_buf", clover_p->frm_buf0, f_size);

  // Set the clover coefficient                                      
  //-------------------------------------------------------------------------
  clover_p->clover_coef = GJP.CloverCoeff();
}


//---------------------------------------------------------------------------
// void clover_end()
//---------------------------------------------------------------------------
// Purpose:
//   free all the memory allocated in clover_init().
//---------------------------------------------------------------------------
void clover_end(Clover *clover_p)
{
  char *cname = "";  
  char *fname = "clover_end";
  VRB.Func(cname,fname);

  // free temporary buffers to be used by DiracOpClover::MatPcDagMatPc
  // and etc.
  //-------------------------------------------------------------------------
  VRB.Sfree(cname,fname, "frm_buf", clover_p->frm_buf0);
  sfree(clover_p->frm_buf0);

  // free memory used by the wilson library.
  //----------------------------------------------------------------
  wilson_end(clover_p->wilson_p);
}


//--------------------------------------------------------------------------
// calculate Y = A X   where  Y, X:   spinors (checkboarded)
//                               A:   compressed Clover matrices (chbd'ed)
//                               n:   number of local sites (chbd'ed)
//--------------------------------------------------------------------------
// Algorithm:
//   A = A0      A1, -A2     A4, -A5    A9, -A10     A16,-A17    A25, -A26
//       1,2     3           6,-7       11,-12       18,-19      27,-28
//       4,5     6,7         8          13,-14       20,-21      29,-30
//       9,10    11,12       13,14      15           22,-23      31,-32
//       16,17   18,19       20,21      22,23        24          33,-34
//       25,26   27,28       29,30      31,32        33,34       35
//--------------------------------------------------------------------------
extern "C" { 
  void clover_mat_mlt_asm(IFloat *out,const IFloat *mat,const IFloat *in,int n);
}
CPS_END_NAMESPACE
#include <stdio.h>
#include <util/dirac_op.h>
#include <time.h>
#include <sys/time.h>
CPS_START_NAMESPACE

void clover_mat_mlt_C(IFloat *out, const IFloat *mat, const IFloat *in, int n) 
{
  clover_mat_mlt_asm(out,mat, in, n) ;
  DiracOp::CGflops += 504*n ;
}


//--------------------------------------------------------------------------
// circular buffer used in clover_mat_mlt_asm.asm
//   const unsigned CBUF_MODE2 = 0xcca52148;  72 words for clover matrix
//   const unsigned CBUF_MODE4 = 0xce318118;  24 words for spinor
//--------------------------------------------------------------------------
//    MODE1 = 0xcb911548 = 11001 01110 01000 10001 0101 0100 1000
//    MODE2 = 0xcca52112 = 11001 10010 10010 10010 0001 0001 0010
//    MODE3 = 0xc98c6106 = 11001 00110 00110 00110 0001 0000 0110
//    MODE4 = 0xcca52112 = 11001 10010 10010 10010 0001 0001 0010
//    MODEa = 0xcca52148 = 11001 10010 10010 10010 0001 0100 1000
//    MODEb = 0xce318118 = 11001 11000 11000 11000 0001 0001 1000 
//    MODEc = 0xcca52124 = 11001 10010 10010 10010 0001 0010 0100
//    MODEd = 0xcb18c10c = 11001 01100 01100 01100 0001 0000 1100
// CBUF_MODEs:       2     1    3   a   b    c   d
// WCNTSET:bit0-7:   18   72    6   72  24  36  12
// gnteeset:8-11:     1    5    1    1   1   1   1
// lookcnt:12-16:    18   17    6   18  24  18  12
// HIGHRQSET:17-21:  18    8    6   18  24
// HISTERHIGH:22-26: 18   14    6   18  24
// WCNT_ENABLE:27:    1
// KICK_ENABLE:28:    0
// CLRBIT:29:         0
// WAITGNTEE:30:      1
// ZFLAG:31:          1  
// purpose:           u   udag chi
//--------------------------------------------------------------------------

CPS_END_NAMESPACE
