#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-08-18 11:57:46 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdsp/glb_cpp_hdw_nos/global_sum.C,v 1.5 2004-08-18 11:57:46 zs Exp $
//  $Id: global_sum.C,v 1.5 2004-08-18 11:57:46 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.5 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdsp/glb_cpp_hdw_nos/global_sum.C,v $
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
#include <util/qcdio.h>                                   // printf 
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
