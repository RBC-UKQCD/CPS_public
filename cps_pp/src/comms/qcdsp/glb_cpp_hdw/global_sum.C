#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-01-13 20:39:13 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdsp/glb_cpp_hdw/global_sum.C,v 1.2 2004-01-13 20:39:13 chulwoo Exp $
//  $Id: global_sum.C,v 1.2 2004-01-13 20:39:13 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.1.2.1  2003/11/05 18:12:13  mike
//  Changing directory structure.
//
//  Revision 1.1.1.1  2003/06/22 13:34:47  mcneile
//  This is the cleaned up version of the Columbia Physics System.
//  The directory structure has been changed.
//  The include paths have been updated.
//
//
//  Revision 1.5  2001/08/16 12:54:07  anj
//  Some fixes follosin the float-> IFloat change, mostly of the (variable
//  anme) IFloat_p -> float_p type.  A few fixes to ensure the test
//  programs use the same level of verbosity throughout, and an update of
//  the regression.pl script to make it more useful. Anj
//
//  Revision 1.4  2001/08/16 10:49:54  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "IFloat".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and IFloat).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.2  2001/06/19 18:12:02  anj
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
//  Revision 1.2  2001/05/25 06:16:02  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: global_sum.C,v $
//  $Revision: 1.2 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdsp/glb_cpp_hdw/global_sum.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------



CPS_END_NAMESPACE
#include<global_sum.h>
CPS_START_NAMESPACE

CPS_END_NAMESPACE
#include<global_sum_supp.h>
#include<global_sum_handler.h>
#include<global_sum_info.h>
CPS_START_NAMESPACE


CPS_END_NAMESPACE
#include <string.h>                           // memmove to handle overlapping
CPS_START_NAMESPACE

//--------------------------------------------------------------------------
// inline void GLOBAL_SUM_VERBOSE(const char *format, int i=0) 
//--------------------------------------------------------------------------
// For the purpose of debugging only.
// Cannot use system verbose because of synchronization problem.
//--------------------------------------------------------------------------
#ifdef GLOBAL_SUM_VERBOSE_ON
CPS_END_NAMESPACE
#include <stdio.h>                                   // printf 
CPS_START_NAMESPACE
#endif

inline void GLOBAL_SUM_VERBOSE(const char *format, int i=0) 
{
#ifdef GLOBAL_SUM_VERBOSE_ON
  printf(format, i);
  global_sum_synch();  
#endif
}
  
//--------------------------------------------------------------------------
// inline function
//--------------------------------------------------------------------------
inline void set_int_zero(int *buf, int size) 
{
  for (int i = 0; i < size; ++i)
    buf[i] = 0;
}



//--------------------------------------------------------------------------
// void broadcast(const void *src, void *dst, int n)
//--------------------------------------------------------------------------
void broadcast(const void *src, void *dst, int n)
{
  int debug_buf[9];

  // the master node copies the memory from src to dst iff dst sits between
  // [src, src+n). 
  //------------------------------------------------------------------------
  if ((int)src < (int)dst && (int)dst < (int)src+n && 
      GLOBAL_SUM_INFO_BUF->is_master) {
    memmove(dst, src, n);   
    src = dst;
    GLOBAL_SUM_VERBOSE("brdcst: debug_buf at 0x%x\n", (unsigned int)debug_buf);
  }
  

  // broadcast
  //------------------------------------------------------------------------
  for (int i = 0; i < GLOBAL_SUM_INFO_BUF->max_num_try; ++i) {
    set_int_zero(debug_buf, sizeof(debug_buf)/sizeof(int));

    global_sum_scu_reset(0); 
    GLOBAL_SUM_VERBOSE("brdcst: reset'ed\n");

    global_sum_synch();                    
    GLOBAL_SUM_VERBOSE("brdcst: sync'ed\n");    

    global_sum_brdcst(src, dst, n);             
    GLOBAL_SUM_VERBOSE("brdcst: brdcst'ed\n"); 

    global_sum_scu_poll();                 
    GLOBAL_SUM_VERBOSE("brdcst: poll'ed\n");    

    int error = global_sum_scu_reset(debug_buf); 
    GLOBAL_SUM_VERBOSE("brdcst: reset'ed\n");

    GLOBAL_SUM_VERBOSE("brdcst: lcl_err 0x%x\n", error);
    error = global_sum_glb_or(error);  
    GLOBAL_SUM_VERBOSE("brdcst: glb_or 0x%x\n", error);

    if (!error)
      return;
  }

  // broadcast failed
  //------------------------------------------------------------------------
  global_sum_handler(SCU_ERR, (unsigned int)debug_buf);  
}





//--------------------------------------------------------------------------
// int global_sum(IFloat *float_p)
//--------------------------------------------------------------------------
int global_sum(IFloat *float_p)
{
  int intArray[3];
  int debug_buf[9];  

  GLOBAL_SUM_VERBOSE("global_sum: debug_buf at 0x%x\n", 
		     (unsigned int)debug_buf);
  
  for (int i = 0; i < GLOBAL_SUM_INFO_BUF->max_num_try; ++i) {
    set_int_zero(debug_buf, sizeof(debug_buf)/sizeof(int));

    global_sum_scu_reset(0); 
    GLOBAL_SUM_VERBOSE("global_sum: reset'ed\n");

    global_sum_get_exp(float_p, intArray);    

    global_sum_synch();                        
    GLOBAL_SUM_VERBOSE("global_sum: sync'ed\n");  

    global_sum_max_or_sum(intArray, 0);        
    GLOBAL_SUM_VERBOSE("global_sum: max'ed\n");  

    global_sum_scu_poll();                     
    GLOBAL_SUM_VERBOSE("global_sum: poll'ed\n");  

    global_sum_scu_reset(debug_buf);           
    GLOBAL_SUM_VERBOSE("global_sum: reset'ed\n");  

    global_sum_brdcst(intArray, intArray, 1);  
    GLOBAL_SUM_VERBOSE("global_sum: brdcst'ed\n");  

    global_sum_scu_poll();                     
    GLOBAL_SUM_VERBOSE("global_sum: poll'ed\n");  

    global_sum_scu_reset(debug_buf);          
    GLOBAL_SUM_VERBOSE("global_sum: reset'ed\n");  

    global_sum_ftoi(float_p, intArray);

    global_sum_synch();                        
    GLOBAL_SUM_VERBOSE("global_sum: sync'ed\n"); 

    global_sum_max_or_sum(intArray, 2);        
    GLOBAL_SUM_VERBOSE("global_sum: sum'ed\n");  

    global_sum_scu_poll();                     
    GLOBAL_SUM_VERBOSE("global_sum: poll'ed\n");  

    global_sum_scu_reset(debug_buf);           
    GLOBAL_SUM_VERBOSE("global_sum: reset'ed\n");  

    global_sum_brdcst(intArray, intArray, 2);  
    GLOBAL_SUM_VERBOSE("global_sum: brdcst'ed\n");  

    global_sum_scu_poll();                     
    GLOBAL_SUM_VERBOSE("global_sum: poll'ed\n");  

    int error = global_sum_scu_reset(debug_buf); 
    GLOBAL_SUM_VERBOSE("global_sum: reset'ed\n");

    error = global_sum_glb_or(error);         
    GLOBAL_SUM_VERBOSE("global_sum: err=\n",error);
    
    if (!error) {
      int true_exp = global_sum_itof(float_p, intArray);
      if (true_exp < 128)
	return true_exp;
      global_sum_handler(FLT_OVERFLOW, true_exp);
    }
  }

  global_sum_handler(SCU_ERR, (unsigned int)debug_buf);  
  return -1;  
}


CPS_END_NAMESPACE
