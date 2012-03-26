#include<config.h>
#ifdef USE_QMP
//#include<qalloc.h>
CPS_START_NAMESPACE
//-------------------------------------------------------------------
/*!\file
  \brief  Definition of glb_sum_dir routine.

*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qmp/glb_cpp/glb_sum_dir.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//--------------------------------------------------------------
// glb_sum_dir
//
// Sum over all nodes along a direction
// (0,1,2,3,4) <-> (x,y,z,t,s)
//--------------------------------------------------------------
CPS_END_NAMESPACE
#include<comms/glb.h>
#include<comms/scu.h>
#include<util/gjp.h>
#include<comms/double64.h>
#include <comms/sysfunc_cps.h>
CPS_START_NAMESPACE
#ifndef USE_QMP
#define USE_QMP
#endif


static Double64 *transmit_buf = NULL;
static Double64 *receive_buf = NULL;
static Double64 *gsum_buf = NULL;

//----------------------------------------------------------------------
/*!
  \param float_p The number to be summed.
  \param dir The direction in which to sum; one of {0, 1, 2, 3, 4},
  corresponding to {x, y, z, t, s}.
  \post The number pointed to by \a float_p is summed over all nodes along the
  \a dir direction, \e i.e. over each strip of nodes where the grid
  coordinates in all other directions are constant.
  That sum is written back to \a float_p, which is identical on all nodes in
  this strip.

  \ingroup comms
*/
//---------------------------------------------------------------------- 

void glb_sum_dir(Float * float_p, int dir)
{
  int NP[5] = {GJP.Xnodes(), 
	       GJP.Ynodes(), 
	       GJP.Znodes(), 
	       GJP.Tnodes(), 
	       GJP.Snodes()};

  int COOR[5] = {GJP.XnodeCoor(), 
		 GJP.YnodeCoor(), 
		 GJP.ZnodeCoor(), 
		 GJP.TnodeCoor(), 
		 GJP.SnodeCoor()}; 

#ifndef UNIFORM_SEED_NO_COMMS
  if (transmit_buf == NULL){
    QMP_mem_t *mem_t= QMP_allocate_memory (sizeof(Double64));
      transmit_buf = (Double64 *)QMP_get_memory_pointer(mem_t);
  }
  if (receive_buf == NULL){
    QMP_mem_t *mem_t= QMP_allocate_memory (sizeof(Double64));
      receive_buf = (Double64 *)QMP_get_memory_pointer(mem_t);
  }
  if (gsum_buf == NULL){
    QMP_mem_t *mem_t= QMP_allocate_memory (sizeof(Double64));
      gsum_buf = (Double64 *)QMP_get_memory_pointer(mem_t);
  }

  // Sum along dir
  //--------------------------------------------------------------
  *gsum_buf = *float_p;

  *transmit_buf = *gsum_buf;

  int itmp;
  for (itmp = 1; itmp < NP[dir]; itmp++) {

    QMP_msgmem_t msgmem[2];
    QMP_msghandle_t msghandle[2];
    QMP_msghandle_t sndrcv;
    msgmem[0] = QMP_declare_msgmem((void *)transmit_buf, sizeof(Double64));
    msghandle[0] = QMP_declare_send_relative(msgmem[0], dir, 1, 0);
    msgmem[1] = QMP_declare_msgmem((void *)receive_buf, sizeof(Double64));
    msghandle[1] = QMP_declare_receive_relative(msgmem[1], dir, -1, 0);
    sndrcv = QMP_declare_multiple(msghandle, 2);
    
    QMP_start(sndrcv);
    QMP_status_t status = QMP_wait(sndrcv);
    if (status != QMP_SUCCESS)
      QMP_error("Communication error in glb_sum_gimp:%s\n", QMP_error_string(status));
    QMP_free_msghandle(sndrcv);
    QMP_free_msgmem(msgmem[0]);
    QMP_free_msgmem(msgmem[1]);


    *gsum_buf += *receive_buf;
    *transmit_buf = *receive_buf;
  }


  // Broadcast the result of node with dir coordinate == 0
  //--------------------------------------------------------------

  if(COOR[dir] != 0) {
    *gsum_buf = 0;
  }
    
  *transmit_buf = *gsum_buf;
  
  for (itmp = 1; itmp < NP[dir]; itmp++) {
    QMP_msgmem_t msgmem[2];
    QMP_msghandle_t msghandle[2];
    QMP_msghandle_t sndrcv;
    msgmem[0] = QMP_declare_msgmem((void *)transmit_buf, sizeof(Double64));
    msghandle[0] = QMP_declare_send_relative(msgmem[0], dir, 1, 0);
    msgmem[1] = QMP_declare_msgmem((void *)receive_buf, sizeof(Double64));
    msghandle[1] = QMP_declare_receive_relative(msgmem[1], dir, -1, 0);
    sndrcv = QMP_declare_multiple(msghandle, 2);
    
    QMP_start(sndrcv);
    QMP_status_t status = QMP_wait(sndrcv);
    if (status != QMP_SUCCESS)
      QMP_error("Communication error in glb_sum_gimp:%s\n", QMP_error_string(status));
    QMP_free_msghandle(sndrcv);
    QMP_free_msgmem(msgmem[0]);
    QMP_free_msgmem(msgmem[1]);
    
    *gsum_buf += *receive_buf;
    *transmit_buf = *receive_buf;
  }


  *float_p = *gsum_buf;
#endif
//  fprintf(stdout,"glb_sum_dir():NP[%d]=%d COOR[%d]=%d sum=%0.16e \n",dir,NP[dir],dir,COOR[dir],float_p);
}


CPS_END_NAMESPACE
#endif
