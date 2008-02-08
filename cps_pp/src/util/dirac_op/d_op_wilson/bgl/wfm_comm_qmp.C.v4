#include <config.h>
#ifdef USE_QMP
#include <comms/sysfunc_cps.h>
#include <util/wilson.h>
#include <util/verbose.h>
#include <util/error.h>
#include <util/gjp.h>
CPS_START_NAMESPACE

/*--------------------------------------------------------------------------*/
/* Externals                                                                */
/*--------------------------------------------------------------------------*/
extern int wfm_max_numchunk;
extern int wfm_numchunk[];
extern IFloat **wfm_send_ad;
extern IFloat **wfm_recv_ad;
extern IFloat *wfm_s_start[8];
extern IFloat *wfm_r_start[8];
extern unsigned long wfm_blklen[8];
extern unsigned long wfm_numblk[8];
extern unsigned long wfm_stride[8];

static int initted=0;
extern int wilson_initted;

void wfm_comm(){

  char *fname="wfm_comm()";
  const int group = 16;
  void *addr[group];
  size_t blksize[group];
  int  numblk[group];
  ptrdiff_t stride[group];
  int index;

  const int MAX_MSGHANDLE=20;
  if (wfm_max_numchunk/group+1 >MAX_MSGHANDLE)
  ERR.General("",fname,"wfm_max_numchunk(%d)/group+1 >MAX_MSGHANDLE",wfm_max_numchunk);

 static   QMP_msgmem_t send_mem[8][MAX_MSGHANDLE];
 static   QMP_msgmem_t recv_mem[8][MAX_MSGHANDLE];
 static   QMP_msghandle_t send_h[8][MAX_MSGHANDLE];
 static   QMP_msghandle_t recv_h[8][MAX_MSGHANDLE];
 static   int pir=0;

  static int wfm_blocks[8];
if (wilson_initted || !initted){
    VRB.Flow("",fname,"wilson_initted=%d initted=%d\n",wilson_initted,initted);
  for(int dir=0;dir<8;dir++) wfm_blocks[dir]=1;
  for(int ig=0; ig<1; ig++){
//    VRB.Flow("",fname,"ig=%d",ig);
    for(int dir=0;dir<8;dir++){
      int sign=1;
      if(dir>3) sign = -1;
//      VRB.Flow("",fname,"dir=%d",dir);
       int n_site=0;
           addr[n_site] = wfm_s_start[dir];
           blksize[n_site] = wfm_blklen[dir];
           numblk[n_site] = wfm_numblk[dir];
           stride[n_site] = wfm_stride[dir];
           n_site++;
       if (initted){
         QMP_free_msghandle(send_h[dir][ig]);
         QMP_free_msghandle(recv_h[dir][ig]);
         QMP_free_msgmem(send_mem[dir][ig]);
         QMP_free_msgmem(recv_mem[dir][ig]);
       }
       if(n_site>0){
         send_mem[dir][ig]  = QMP_declare_strided_array_msgmem(addr,blksize,numblk,stride,n_site);
         send_h[dir][ig] = QMP_declare_send_relative(send_mem[dir][ig],dir%4,sign,0);
         wfm_blocks[dir]=ig+1;
       }

       int r_site=0;
           addr[r_site] = wfm_r_start[dir];
           blksize[r_site] = wfm_blklen[dir];
           numblk[r_site] = wfm_numblk[dir];
           stride[r_site] = wfm_stride[dir];
           r_site++;
       if (n_site!=r_site)
         ERR.General("",fname,"n_site(%d)!=r_site(%d)\n",n_site,r_site);
//       VRB.Flow("",fname,"n_site=%d r_site=%d",n_site,r_site);
       if(r_site>0){
         recv_mem[dir][ig]  = QMP_declare_strided_array_msgmem(addr,blksize,numblk,stride,r_site);
         recv_h[dir][ig] = QMP_declare_receive_relative(recv_mem[dir][ig],dir%4,-sign,0);
       }


    } // dir
  } // ig
  for(int dir=0;dir<8;dir++) 
     VRB.Flow("",fname,"wfm_blocks[%d]=%d",dir,wfm_blocks[dir]);
  pir = CoorT()%2;
//  pir = 0;
  initted=1;
  wilson_initted=0;
} 

#if 0
   for(int dir=0;dir<8;dir++)
     for ( index=0; index <wfm_numchunk[dir];index++){
           Float *tmp_p = wfm_send_ad[dir+8*index];
           if ( (*tmp_p)*(*tmp_p) >0.0001)
           printf("Node %d: wfm_send_ad[%d][%d]=%e\n",UniqueID(),dir,index,*tmp_p); 
   }
#endif

  int dir_g=4;
  for(int ig=0; ig<1; ig++){
    for(index=0;index<8;index++){
           QMP_start(send_h[index][ig]); QMP_start(recv_h[index][ig]);
    }
    for(index=0;index<8;index++){
         QMP_wait(send_h[index][ig]); QMP_wait(recv_h[index][ig]);
    }
  }

#if 0
   for(int dir=0;dir<8;dir++)
     for ( index=0; index <wfm_numchunk[dir];index++){
           Float *tmp_p = wfm_recv_ad[dir+8*index];
           if ( (*tmp_p)*(*tmp_p) >1e-10)
           printf("Node %d: wfm_recv_ad[%d][%d]=%e\n",UniqueID(),dir,index,*tmp_p); 
   }
#endif
}

CPS_END_NAMESPACE
#endif
