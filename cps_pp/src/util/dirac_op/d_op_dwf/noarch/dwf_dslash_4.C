#include<config.h>
#if 0
#include "../sse/sse-dwf_dslash_4.C"
#else
//------------------------------------------------------------------
//
// dwf_dslash_4.C
//
// dwf_dslash_4 is the derivative part of the space-time part of
// the fermion matrix. 
// The in, out fields are defined on the checkerboard lattice
// cb = 0/1 <--> even/odd checkerboard of in field.
// dag = 0/1 <--> Dslash/Dslash^dagger is calculated.
//
//------------------------------------------------------------------

#include<config.h>
#include<util/dwf.h>
#include<util/wilson.h>
#include<util/gjp.h>
#include<util/vector.h>
#include<util/verbose.h>
#include<util/error.h>
#ifdef PARALLEL
#include <comms/sysfunc_cps.h>
#endif

//CK
#include <util/fpconv.h>
#include <util/checksum.h>
#include <util/qioarg.h>

CPS_START_NAMESPACE

#ifdef USE_QMP
//declares external function that is implemented in  util/dirac_op/d_op_wilson/qmp/wilson_dlash_vec.C
void wilson_dslash_vec(IFloat *chi_p_f,
                        IFloat *u_p_f,
                        IFloat *psi_p_f,
                        int cb,
                        int dag,
                        Wilson *wilson_p,
                        int vec_len,
                        unsigned long vec_offset);
#endif

//CK for debugging
static void print_checksum(Float *field, const int &s){
  int nwilson = GJP.VolNodeSites()/2; //how many blocks of 24 floats
  if(GJP.Gparity()) nwilson*=2;
    
  FPConv fp;
  enum FP_FORMAT format = FP_IEEE64LITTLE;
  uint32_t csum(0);
    
  for(int x=0; x<nwilson; x++){
    uint32_t csum_contrib = fp.checksum((char *)(field),24,format);
    csum += csum_contrib;
    field+=24;
  }
    
  QioControl qc;
  csum = qc.globalSumUint(csum);
    
  if(!UniqueID()) printf("s=%d src checksum %u\n",s,csum);
}

void dwf_dslash_4(Vector *out, 
		  Matrix *gauge_field, 
		  Vector *in, 
		  int cb, 
		  int dag, 
		  Dwf *dwf_lib_arg)
{
  int i;
  int ls;
  IFloat *frm_in;
  IFloat *frm_out;
  IFloat *g_field;
  Wilson *wilson_p;
  int size_cb[2];
  int parity=cb;

  //----------------------------------------------------------------
  // Initializations
  //----------------------------------------------------------------
  ls = dwf_lib_arg->ls;
  frm_in = (IFloat *) in;
  frm_out = (IFloat *) out;
  g_field = (IFloat *) gauge_field;
  wilson_p = dwf_lib_arg->wilson_p;
  size_cb[0] = 24*wilson_p->vol[0]; //wilson_p->vol[0] is half of the local lattice volume
  size_cb[1] = 24*wilson_p->vol[1];
  
  if(GJP.Gparity()){ //CK: 2 4d fields on each ls slice
    size_cb[0]*=2;

    size_cb[1]*=2;
  }
#ifndef USE_WILSON_DSLASH_VEC
//#ifndef USE_QMP
  //if USE_QMP declared then we use the vectorised wilson dslash
  //----------------------------------------------------------------
  // Apply 4-dimensional Dslash
  //----------------------------------------------------------------
  int vec_len=1;
  for(i=0; i<ls; i+= vec_len){

    // parity of 4-D checkerboard
    //------------------------------------------------------------
    parity = (i + cb) % 2;

    // Apply on 4-dim "parity" checkerboard part
    //------------------------------------------------------------
    if(vec_len==1){

#if 0
      print_checksum((Float*)frm_in,i); //DEBUG
#endif

      wilson_dslash(frm_out, g_field, frm_in, parity, dag, wilson_p);

#if 0
      print_checksum((Float*)frm_out,i); //DEBUG
#endif


    }else{
#if TARGET == NOARCH
      ERR.NotImplemented("","dwf_dslash_4(..)","wilson_dslash_two() doesn't exist for NOARCH target\n");
#else
      wilson_dslash_two(frm_out, frm_out+size_cb[parity], g_field, frm_in, frm_in+size_cb[parity], parity, 1-parity,dag, wilson_p);
#endif
    }
    frm_in = frm_in + vec_len*size_cb[parity];
    frm_out = frm_out + vec_len*size_cb[parity];
  }
#else
  //----------------------------------------------------------------
  // Apply vectorized 4-dimensional Dslash
  //----------------------------------------------------------------
  
  parity = cb; //CK parity was not intialised previously
  wilson_dslash_vec(frm_out, g_field, frm_in, cb, dag, wilson_p,ls/2,2*size_cb[parity]);

  frm_in = frm_in + size_cb[parity];
  frm_out = frm_out + size_cb[parity];

  parity = 1-cb; //CK parity was not intialised previously
  wilson_dslash_vec(frm_out, g_field, frm_in, 1-cb, dag, wilson_p,ls/2,2*size_cb[parity]);
#endif


}

CPS_END_NAMESPACE
#endif
