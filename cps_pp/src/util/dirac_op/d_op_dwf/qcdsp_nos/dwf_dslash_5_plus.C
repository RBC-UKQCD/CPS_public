#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:46 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_dwf/qcdsp_nos/dwf_dslash_5_plus.C,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
//  $Id: dwf_dslash_5_plus.C,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.6  2001/08/16 10:50:18  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "IFloat".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and IFloat).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.4  2001/07/03 17:01:01  anj
//
//  Multiple minor alterations to change some #include's from referring to
//  files relative to the top-level source directory to referring to files
//  relative to the source-file positions.  This alteration makes the code
//  backwards compatable with the make structure of QCDSP, although this
//  may have to be changed to a more usual form in the future. Anj.
//
//  Revision 1.3  2001/06/21 15:40:17  anj
//  Updated the _TARTAN ifdefs, using PARALLEL instead (where appropriate).Anj
//
//  Revision 1.2  2001/06/19 18:12:40  anj
//  Serious ANSIfication.  Plus, degenerate double64.h files removed.
//  Next version will contain the new nga/include/double64.h.  Also,
//  Makefile.gnutests has been modified to work properly, propagating the
//  choice of C++ compiler and flags all the way down the directory tree.
//  The mpi_scu code has been added under phys/nga, and partially
//  plumbed in.
//
//  Everything has newer dates, due to the way in which this first alteration was handled.
//
//  Anj.
//
//  Revision 1.2  2001/05/25 06:16:05  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: dwf_dslash_5_plus.C,v $
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_dwf/qcdsp_nos/dwf_dslash_5_plus.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// dwf_dslash_5_plus.C
//
// dwf_dslash_5_plus is the derivative part of the 5th direction
// part of the fermion matrix. This routine accumulates the result
// on the out field 
// The in, out fields are defined on the checkerboard lattice.
// The action of this operator is the same for even/odd
// checkerboard fields because there is no gauge field along
// the 5th direction.
// dag = 0/1 <--> Dslash/Dslash^dagger is calculated.
//
//
// Storage order for DWF fermions
//------------------------------------------------------------------
//  
//  |     |
//  | |r| |
//  | |i| |
//  |     |
//  | |r| | = |spin comp|
//  | |i| |
//  |     |
//  | |r| |
//  | |i| |
//  |     |
//  
//  
//  |             |
//  | |spin comp| |
//  |             |
//  |             |
//  | |spin comp| | 
//  |             | = |spinor|
//  |             |
//  | |spin comp| |
//  |             |
//  |             |
//  | |spin comp| |
//  |             |
//  
//  
//  |            |
//  |  |spinor|  |
//  |  |spinor|  |
//  |  |spinor|  |
//  |  |spinor|  |
//  |  |spinor|  |
//  |     .      | = |s-block|   The spinors are arranged in Wilson
//  |     .      |               order with odd - even 4d-checkerboard
//  |     .      |               storage.
//  |evn/odd vol |
//  |     .      |
//  |  |spinor|  |
//  |            |
//  
//  
//  |                |
//  | |s-block even| |  For even chckerboard
//  | |s-block odd|  |
//  | |s-block even| |
//  | |s-block odd|  |
//  |       .        |
//  |       .        |
//  |       .        |
//  |                |
//
//
//  |                |
//  | |s-block odd|  |  For odd chckerboard
//  | |s-block even| |
//  | |s-block odd|  |
//  | |s-block even| |
//  |       .        |
//  |       .        |
//  |       .        |
//  |                |
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include<config.h>
CPS_START_NAMESPACE
#ifdef PARALLEL
CPS_END_NAMESPACE
#include <sysfunc.h>
CPS_START_NAMESPACE
#endif
CPS_END_NAMESPACE
#include<util/dwf.h>
#include<util/gjp.h>
#include<util/vector.h>
#include<util/verbose.h>
#include<util/error.h>
#include<util/smalloc.h>
#include<comms/nga_reg.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------
// Declarations local to this file
//------------------------------------------------------------------
void getPlusFiveData(IFloat *rcv_buf, IFloat *send_buf);
void getMinusFiveData(IFloat* rcv_buf, IFloat* send_buf);
void setScu(void);

//------------------------------------------------------------------
// void dwf_dslash_5_plus(Vector *, Vector *, Float int, Dwf*):  
//------------------------------------------------------------------
void dwf_dslash_5_plus(Vector *out, 
		       Vector *in, 
		       Float mass,
		       int dag, 
		       Dwf *dwf_lib_arg)
{
  int x;
  int s;

// Initializations
//------------------------------------------------------------------
  int local_ls = GJP.SnodeSites(); 
  int s_nodes = GJP.Snodes();
  int s_node_coor = GJP.SnodeCoor();
  int vol_4d_cb = dwf_lib_arg->vol_4d / 2;
  int ls_stride = 24 * vol_4d_cb;
  IFloat *f_in;
  IFloat *f_out;
  IFloat *f_temp;
  IFloat *comm_buf = dwf_lib_arg->comm_buf;
  IFloat two_over_a5 = 2.0 * GJP.DwfA5Inv();
  IFloat neg_mass_two_over_a5 = -2.0 * mass * GJP.DwfA5Inv();
 
  setScu();

// [1 + gamma_5] term (if dag=1 [1 - gamma_5] term)
//
// out[s] = [1 + gamma_5] in[s-1]
//------------------------------------------------------------------
  f_in  = (IFloat *) in;
  f_out = (IFloat *) out;
  if(dag == 1){
    f_in  =  f_in + 12;
    f_out = f_out + 12;
  }
  f_out = f_out + ls_stride; 
  for(s=1; s < local_ls; s++){

    for(x=0; x<vol_4d_cb; x++){

      fTimesV1PlusV2(f_out, two_over_a5, f_in, f_out, 12);

      f_in  =  f_in + 24;
      f_out = f_out + 24;
    }
  }


// [1 + gamma_5] for lower boundary term (if dag=1 [1 - gamma_5] term)
// If there's only one node along fifth direction, no communication
// is necessary; Otherwise data from adjacent node in minus direction
// will be needed.
// If the lower boundary is the s=0 term
// out[0] = - m_f * [1 + gamma_5] in[ls-1]
// else, out[s] = [1 + gamma_5] in[s-1]
//
//------------------------------------------------------------------

  f_in  = (IFloat *) in;  
  f_in = f_in + (local_ls-1)*ls_stride; 
  f_out = (IFloat *) out;
  
  if(dag == 1){
    f_in  =  f_in + 12;
    f_out = f_out + 12;
  }
  
  for(x=0; x<vol_4d_cb; x++){
    
    f_temp = f_in;
    
    if (s_nodes != 1 ) {
      f_temp = comm_buf;
      getMinusFiveData(f_temp, f_in);
    }
    
    if(s_node_coor == 0) { 
      fTimesV1PlusV2(f_out, neg_mass_two_over_a5, f_temp, f_out, 12);
    }
    else {
      fTimesV1PlusV2(f_out, two_over_a5, f_temp, f_out, 12);
    }
    
    f_in  =  f_in + 24;
    f_out = f_out + 24;
  }


// [1 - gamma_5] term (if dag=1 [1 + gamma_5] term)
// 
// out[s] = [1 - gamma_5] in[s+1]
//------------------------------------------------------------------
  f_in  = (IFloat *) in;
  f_out = (IFloat *) out;
  if(dag == 0){
    f_in  =  f_in + 12;
    f_out = f_out + 12;
  }
  f_in = f_in + ls_stride;
  for(s=0; s < local_ls - 1; s++){

    for(x=0; x<vol_4d_cb; x++){

      fTimesV1PlusV2(f_out, two_over_a5, f_in, f_out, 12);

      f_in  =  f_in + 24;
      f_out = f_out + 24;
    }
  }


// [1 - gamma_5] for upper boundary term (if dag=1 [1 + gamma_5] term)
// If there's only one node along fifth direction, no communication
// is necessary; Otherwise data from adjacent node in minus direction
// will be needed.
// If the upper boundary is the s=ls term
// out[ls-1] = - m_f * [1 - gamma_5] in[0]
// else out[s] = [1 - gamma_5] in[s+1]
//
//------------------------------------------------------------------

  f_in  = (IFloat *) in;
  f_out = (IFloat *) out;

  if(dag == 0){
    f_in  =  f_in + 12;
    f_out = f_out + 12;
  }

  f_out = f_out + (local_ls-1)*ls_stride;
  for(x=0; x<vol_4d_cb; x++){

    f_temp = f_in;
   
    if (s_nodes != 1 ) {
      f_temp = comm_buf;
      getPlusFiveData(f_temp, f_in);
    }

    if(s_node_coor == s_nodes - 1) { 
      fTimesV1PlusV2(f_out, neg_mass_two_over_a5, f_temp, f_out, 12);
    }
    else {
      fTimesV1PlusV2(f_out, two_over_a5, f_temp, f_out, 12);
    }
    
    f_in  =  f_in + 24;
    f_out = f_out + 24;
  }


}


#ifdef _TARTAN

static volatile unsigned* dsp_scu_base0x50 = 
        (volatile unsigned* )(DSP_SCU_BASE + 0x50);
static volatile unsigned* dsp_scu_base0x10 = 
        (volatile unsigned* )(DSP_SCU_BASE + 0x10);
static volatile unsigned* dsp_scu_base0x8 = 
        (volatile unsigned* )(DSP_SCU_BASE + 0x8);
static volatile unsigned* dsp_scu_base0x0 = 
        (volatile unsigned* )(DSP_SCU_BASE);

//------------------------------------------------------------------
// getPlusFiveData(IFloat *rcv_buf, IFloat *send_buf)
// Communication routine along the plus fifth direction
// to send/receive 12 words.
// Used instead of the standard getPlusData for increased 
// performance (it bypasses the qos routines).
//------------------------------------------------------------------
void getPlusFiveData(IFloat *rcv_buf, IFloat *send_buf)
{
  *( (IFloat**)(dsp_scu_base0x8+dwfso_wire_map[1]) ) = send_buf;
  *( (IFloat**)(dsp_scu_base0x0+dwfso_wire_map[0]) ) = rcv_buf;
    while( *(dsp_scu_base0x10+dwfso_wire_map[0]) ) ;
    while( *(dsp_scu_base0x10+dwfso_wire_map[1]) ) ;
}


//------------------------------------------------------------------
// getMinusFiveData(IFloat* rcv_buf, IFloat *send_buf,)
// Communication routine along the minus fifth direction
// to send/receive 12 words.
// Used instead of the standard getMinusData for increased 
// performance (it bypasses the qos routines).
//------------------------------------------------------------------
void getMinusFiveData(IFloat* rcv_buf, IFloat *send_buf)
{
  *( (IFloat**)(dsp_scu_base0x8+dwfso_wire_map[0]) ) = send_buf;
  *( (IFloat**)(dsp_scu_base0x0+dwfso_wire_map[1]) ) = rcv_buf;
    while( *(dsp_scu_base0x10+dwfso_wire_map[1]) ) ;
    while( *(dsp_scu_base0x10+dwfso_wire_map[0]) ) ;
}


//------------------------------------------------------------------
// setScu()
//------------------------------------------------------------------
void setScu(void)
{
  *((int*)(dsp_scu_base0x50+dwfso_wire_map[0]))=(1<<22)|(12<<12)|(0&0x0fff);
  *((int*)(dsp_scu_base0x50+dwfso_wire_map[1]))=(1<<22)|(12<<12)|(0&0x0fff);
}


#else

//------------------------------------------------------------------
// Fake communication routine
//------------------------------------------------------------------
void getPlusFiveData(IFloat *rcv_buf, IFloat *send_buf)
{
  for(int i = 0; i < 12; ++i) {
    *rcv_buf++ = *send_buf++;
  }
}

//------------------------------------------------------------------
// Fake communication routine
//------------------------------------------------------------------
void getMinusFiveData(IFloat* rcv_buf, IFloat* send_buf)
{
  for(int i = 0; i < 12; ++i) {
    *rcv_buf++ = *send_buf++;
  }
}

//------------------------------------------------------------------
// Fake communication routine
//------------------------------------------------------------------
void setScu(void)
{
}

#endif


CPS_END_NAMESPACE
