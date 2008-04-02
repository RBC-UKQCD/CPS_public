/****************************************************************************/
/*                                                                          */
/* PAB: 10/1/2001 Scheme for QCDOC                                          */
/* We lay out the 2-spinors in a very different way, but keep the 4spinor   */
/* and gauge layout the same                                                */
/*                                                                          */
/* 4spinor [Npar][Npsite][N4spin][Ncol][Ncomplex]                           */
/* Gauge   [Npar][Npsite][Ndim ][Nrow][Ncol][Ncomplex]                      */
/* 2spinor:                                                                 */
/*         Body[Npsite][Ndim ][Ncol+1][N2spin][Ncomplex]                    */
/*         Sendx[NboundX][Ncol+1][N2spin]                                   */
/*         Sendy[NboundY][Ncol+1][N2spin]                                   */
/*         Sendz[NboundZ][Ncol+1][N2spin]                                   */
/*         Sendt[NboundT][Ncol+1][N2spin]                                   */
/*         Recvx[NboundX][Ncol+1][N2spin]                                   */
/*         Recvy[NboundY][Ncol+1][N2spin]                                   */
/*         Recvz[NboundZ][Ncol+1][N2spin]                                   */
/*         Recvt[NboundT][Ncol+1][N2spin]                                   */
/* Rationale:                                                               */
/*         The PEC has little/no performance hit for scattering operations  */
/*         when a substantial part of the PEC register is written in a lump */
/*         The performance hit on gather operations is large due to prefetch*/
/*         misses.                                                          */
/*         We can scatter during both the project and su3-project stages    */
/*         to perform all site shifts with full bandwidth.                  */
/*         This fills the send buffers                                      */
/*         Sends have linear accesses (actually strided but PEC misses will */
/*         not occur).                                                      */
/*         I believe the X and T can be Block-Strided directly into the face*/
/*         but that at least the Y and Z need to be received into a buffer  */
/*         and then placed in manually since the pattern isn't block/strided*/ 
/*         This will have linear accesses and be fast.                      */
/*         Now, the pay off for this alteration is the following:           */
/*         The both reconstruct-su3 and reconstruct passes have 2 linear    */
/*         access streams: chi and U and psi and chi respectively. Thus they*/
/*         will take essentially no PEC register misses and run at full     */
/*         bandwidth                                                        */
/*                                                                          */
/****************************************************************************/

#include <stdlib.h>
#include <stdio.h>

#include "wfm.h"
bool wfm_bgl_def = true;

void wfm::init(WilsonArg *wilson_p)  /* pointer to Wilson type structure    */
{
  int spinor_words;             /* size of the spinor field on the         */
				/* sublattice checkerboard                 */

  int half_spinor_words;        /* size of the spin-projected "half_spinors*/
                                /* on the sublattice checkerboard including*/
                                /* the communications padding              */

  int slx;                          /* x-direction size of node sublattice */
  int sly;                          /* y-direction size of node sublattice */
  int slz;                          /* z-direction size of node sublattice */
  int slt;                          /* t-direction size of node sublattice */
  int i;
  int mu;

  SloppyPrecision = wilson_p->SloppyPrecision;
  WFM_BGL             = wilson_p->WFM_BGL;


  if ( isBoss() ) printf("wfm::init setting up BG/L MMU state\n");
  mmu_optimise();
  mmu_print();

//  CoreCount( wilson_p->CoreCount );
  CoreCount( 1 );

  if ( WFM_BGL ) PAD_HALF_SPINOR_SIZE = 12;
  else  PAD_HALF_SPINOR_SIZE = 16;

  if ( WFM_BGL && (nthread > 1) && SloppyPrecision ) { 
    if ( isBoss() ) printf("Bagel does not maintain L1 coherence in dual core + single precision mode on BlueGene\n");
    if ( isBoss() ) printf("Get on to IBM to give me access to SWOA MMU options, or even better a non-cache image of DRAM\n");
    if ( isBoss() ) printf("If they give me the tools, I'm happy to do the heroics of mainting sfw coherence\n");
    if ( isBoss() ) printf("Bagel insanity check exiting\n");
    exit(-1);
  }

  IR = wilson_p->instruction_reg_num;
/*--------------------------------------------------------------------------*/
/* Set sublattice direction sizes                                           */
/*--------------------------------------------------------------------------*/
  local_latt[0] =  wilson_p->local_latt[0];
  local_latt[1] =  wilson_p->local_latt[1];
  local_latt[2] =  wilson_p->local_latt[2];
  local_latt[3] =  wilson_p->local_latt[3];
  slx = local_latt[0];
  sly = local_latt[1];
  slz = local_latt[2];
  slt = local_latt[3];

#ifdef USE_COMMS_QMP
  QMP_bool_t qmp_inited=QMP_is_initialized();
  if( !qmp_inited ) { 
	if ( isBoss() ) printf("QMP_not_initialized\n");
        exit(-1);
  }
  const int *ncoor = QMP_get_logical_coordinates();
  base_parity =(ncoor[0]*local_latt[0] 
              + ncoor[1]*local_latt[1]
              + ncoor[2]*local_latt[2]
              + ncoor[3]*local_latt[3])&0x1;

#else
  base_parity = 0;
#endif


/*--------------------------------------------------------------------------*/
/* Set periodic wrap back or not                                            */
/*--------------------------------------------------------------------------*/
  local_comm[0] = wilson_p->local_comm[0];
  local_comm[1] = wilson_p->local_comm[1];
  local_comm[2] = wilson_p->local_comm[2];
  local_comm[3] = wilson_p->local_comm[3];


/*-----------------------------------------------------------------------*/
/* compute the subgrd volume of each chkbd ... at least two local dims   */
/* must be even for this code to be correct.                             */
/*-----------------------------------------------------------------------*/
  vol = (slx * sly * slz * slt)/2;
  
  nbound[0] = (sly * slz * slt)/2; 
  nbound[1] = (slx * slz * slt)/2;
  nbound[2] = (slx * sly * slt)/2;
  nbound[3] = (slx * sly * slz)/2;

  allbound  = nbound[0]
            + nbound[1]
            + nbound[2]
            + nbound[3];

  if ( nbound[0] * slx * 2 != (slx*sly*slz*slt) ) {
    if ( isBoss() ) printf("wfm::init Even x logic bomb\n");
    exit(-1);
  }
  if ( nbound[1] * sly * 2 != (slx*sly*slz*slt) ) {
    if ( isBoss() ) printf("wfm::init Even y logic bomb\n");
    exit(-1);
  }
  if ( nbound[2] * slz * 2 != (slx*sly*slz*slt) ) {
    if ( isBoss() ) printf("wfm::init Even z logic bomb\n");
    exit(-1);
  }
  if ( nbound[3] * slt * 2 != (slx*sly*slz*slt) ) {
    if ( isBoss() ) printf("wfm::init Even t logic bomb\n");
    exit(-1);
  }

  /*------------------------------------------------------------------------*/
  /* Check shape                                                            */
  /*------------------------------------------------------------------------*/
  if ( (slx&1)  ) {
    if ( isBoss() ) printf("Bagel is refusing to run as x-sub latt is odd\n");
    exit(-1);
  }
  if ( (sly&1) &&(slz&1)&&(slt&1)  ) {
    if ( isBoss() ) printf("Bagel is refusing to run as y,z,t sub latts are all odd\n");
    exit(-1);
  }


/*--------------------------------------------------------------------------*/
/* Reserve memory for 1  temporary spinor (needed by mdagm)                 */
/*--------------------------------------------------------------------------*/
  spinor_words = SPINOR_SIZE * vol;

  spinor_tmp = (Float *)ALLOC(spinor_words*sizeof(Float)*2);
//  VRB.Flow(cname,fname,"spinor_tmp=%p\n",spinor_tmp);
#ifdef USE_QALLOC
  // If we used QALLOC, and the ALLOC macro failed we can try 
  // qalloc but without the QFAST flag. Even tho the spinor_tmp is
  // not communicated we leave the QCOMMS bit on in case it puts 
  // spinor tmp into a better place in the memory map
  if(spinor_tmp == 0) {
     if ( isBoss() ) printf("BAGEL: Warning spinor_tmp has spilled out of Edram\n");
     spinor_tmp = (Float *) qalloc(QCOMMS,spinor_words*sizeof(Float)*2); 
  }
#endif  // USE QALLOC

  if(spinor_tmp == 0){
    if ( isBoss() ) printf("wfm::spinor_tmp allocate\n");
    exit(-1);
  }
/*--------------------------------------------------------------------------*/
/* Reserve memory for the 4 forward and 4 backward spin projected half      */ 
/* spinors.                                                                 */
/*--------------------------------------------------------------------------*/


  /*PAB 10/1/2001 */
  half_spinor_words = NMinusPlus * ND * PAD_HALF_SPINOR_SIZE * vol;

#ifndef USE_COMMS_QMP  
  two_spinor = (Float *)ALLOC(half_spinor_words*sizeof(Float));

#ifdef USE_QALLOC
  // If we are using QALLOC and the ALLOC macro failed we can still 
  // try to get slow memory. Leave on the QCOMMS bit for good memory map
  // placement
  if(two_spinor == 0) {
     if ( isBoss() ) printf("BAGEL : warning two spinors have spilled out of Edram\n");
     two_spinor = (Float *)qalloc(QCOMMS,half_spinor_words*sizeof(Float));
  }
#endif // USE_QALLOC

  if(two_spinor == 0){
    if ( isBoss() ) printf("wfm::two_spinor allocate\n");
    exit(-1);
  }

#else

  // Since two spinor is now communicated because of the Tface 
  // receive I have to allocate it in the style of QMP
  two_spinor_mem_t = QMP_allocate_aligned_memory(
                                             half_spinor_words*sizeof(Float),
	                                     WFM_ALIGN_ARG,
					     (QMP_MEM_COMMS|QMP_MEM_FAST));

  if( two_spinor_mem_t == 0x0 ) { 
    // Try slow allocation
    two_spinor_mem_t = QMP_allocate_aligned_memory(
	                                      half_spinor_words*sizeof(Float),
					      WFM_ALIGN_ARG,
                                              QMP_MEM_COMMS);

    if( two_spinor_mem_t == 0x0 ) { 
      if ( isBoss() ) printf("wfm_init::could not allocate two spinor_mem_t\n");
      exit(-1);
    }
  }
  two_spinor = (Float *)QMP_get_memory_pointer(two_spinor_mem_t);
  if (two_spinor == 0x0) { 
    if ( isBoss() ) printf("wfm::init QMP_get_memory_pointer returned NULL pointer from non NULL QMP_mem_t\n");
    exit(-1);
  } 
#endif

 /*--------------------------------------------------------------------------*/
 /* Reserve memory for the 4 forward and 4 backward spin projected half      */
 /* spinors.                                                                 */
 /*--------------------------------------------------------------------------*/
#ifdef USE_COMMS_SPI
  unsigned long long total_length=0;
  unsigned long long face_length[4];
  for ( int pm = 0;pm<2;pm++ ) {
    for ( mu = 0 ; mu < 4 ; mu++) {
      face_length[mu] = PAD_HALF_SPINOR_SIZE * nbound[mu]*sizeof(Float);
      if (face_length[mu]%128)
        face_length[mu] += 128-(face_length[mu]%128);
      total_length += face_length[mu];
    }
    Float *tmp_buf = (Float *)ALLOC(total_length);
    for ( mu = 0 ; mu < 4 ; mu++) {
      recv_bufs[pm][mu] = tmp_buf;
      tmp_buf += (face_length[mu]/sizeof(Float));
//      if (isBoss()) printf("recv_bufs[%d][%d]=%p %d\n", pm,mu, recv_bufs[pm][mu], recv_bufs[pm][mu]-recv_bufs[0][0]);
    }
    tmp_buf = (Float *)SEND_ALLOC(total_length);
    for ( mu = 0 ; mu < 4 ; mu++) {
      send_bufs[pm][mu] = tmp_buf;
      tmp_buf += (face_length[mu]/sizeof(Float));
//      if (isBoss()) printf("send_bufs[%d][%d]=%p %d\n", pm,mu, send_bufs[pm][mu], send_bufs[pm][mu]-send_bufs[0][0]);
    }
  }
#else
  for ( int pm = 0;pm<2;pm++ ) {
    for ( mu = 0 ; mu < 4 ; mu++) {

      half_spinor_words = PAD_HALF_SPINOR_SIZE * nbound[mu];

      // These things are (potentially) communicated so need QMP Style 
      // allocation if using QMP
      //
      // Note: I am allocating the buffers in all directions regardless
      // of whether we are communicating in that dir or not (Copying CPS)
#ifndef USE_COMMS_QMP

      // Not using QMP
      recv_bufs[pm][myyu] = (Float *)ALLOC(half_spinor_words*sizeof(Float));
//      if (isBoss())  printf("recv_bufs[%d][%d]=%p %x\n", pm,mu, recv_bufs[pm][mu], recv_bufs[pm][mu]-recv_bufs[0][0]);
#ifdef  USE_QALLOC

      // If ALLOC fails try slow memory but with QCOMMS bit still set
      if( recv_bufs[pm][mu] == 0x0 ) 
	  recv_bufs[pm][mu] = (Float *)qalloc(QCOMMS, half_spinor_words*sizeof(Float));
#endif
      if(recv_bufs[pm][mu] == 0){
	  if ( isBoss() ) printf("wfm::recv_bufs allocate\n");
	  exit(-1);
      }

      send_bufs[pm][mu]=(Float *)SEND_ALLOC(half_spinor_words*sizeof(Float));
//      if (isBoss()) printf("send_bufs[%d][%d]=%p %x\n", pm,mu, send_bufs[pm][mu], send_bufs[pm][mu]-send_bufs[0][0]);
#ifdef USE_QALLOC

      // If SEND ALLOC macro fails try slow memory but with QNONCACHE bit
      // still set
      if( send_bufs[pm][mu] == 0 ) 
        send_bufs[pm][mu]=(Float *)qalloc(QNONCACHE, half_spinor_words*sizeof(Float));
#endif

      if(send_bufs[pm][mu] == 0){
        if ( isBoss() ) printf("wfm::send_bufs allocate\n");
        exit(-1);
      }
#else
      /* QMP memory allocation: A little involved */
      /* Must allocate "opaque" QMP_mem_t first and then get 
         aligned pointer out of it. It's either what is below or a 
         very complicated send alloc */

      /* Peter in the CPS allocs recv_bufs with ALLOC = QCOMMS|FAST */
      recv_bufs_mem_t[pm][mu] = QMP_allocate_aligned_memory(half_spinor_words*sizeof(Float),
	                                                      WFM_ALIGN_ARG,
							      (QMP_MEM_COMMS|QMP_MEM_FAST));
      if( recv_bufs_mem_t[pm][mu] == 0x0 ) {
        // If QMP_allocate memory fails with FAST, try SLOW but keep COMMS
        recv_bufs_mem_t[pm][mu] = QMP_allocate_aligned_memory(half_spinor_words*sizeof(Float),
								 WFM_ALIGN_ARG,
								 QMP_MEM_COMMS);
        if( recv_bufs_mem_t[pm][mu] == 0x0 ) { 
	  if ( isBoss() ) printf("wfm::init recv_bufs_mem_t[%d][%d]: QMP_allocate_aligned_memory returned NULL\n", pm, mu);
	  exit(-1);
        }
      }
	
      /* Now get the aligned pointer */
      recv_bufs[pm][mu] =(Float *)QMP_get_memory_pointer(recv_bufs_mem_t[pm][mu]);
	
      if( recv_bufs[pm][mu] == 0x0 ) { 
        if ( isBoss() ) printf("wfm::init recv_bufs[%d][%d]: NULL aligned pointer in non NULL QMP_mem_t struct \n", pm, mu);
	exit(-1);
      }

      /* Now do the same for the send bufs */
      /* In CPS Peter allocates as SEND_ALLOC = QNONCACHE | QFAST */
     send_bufs_mem_t[pm][mu] = QMP_allocate_aligned_memory(half_spinor_words*sizeof(Float),
                                                             WFM_ALIGN_ARG,
                                                             (QMP_MEM_NONCACHE|QMP_MEM_FAST));
     if( send_bufs_mem_t[pm][mu] == 0x0 ) {
       // if allocator fails, try slow but still NONCACHE
       send_bufs_mem_t[pm][mu] = QMP_allocate_aligned_memory(half_spinor_words*sizeof(Float),
                                                               WFM_ALIGN_ARG,
                                                               QMP_MEM_NONCACHE);
       if( send_bufs_mem_t[pm][mu] == 0x0 ) {
         if ( isBoss() ) printf("wfm::init: send_bufs_mem_t[%d][%d]: QMP_allocate_aligned_memory returned NULL\n", pm, mu);
         exit(-1);
       }
     }
            
     /* Now get the aligned pointer */
     send_bufs[pm][mu] =(Float *)QMP_get_memory_pointer(send_bufs_mem_t[pm][mu]);
     if( send_bufs[pm][mu] == 0x0 ) {
       if ( isBoss() ) printf("wfm::init send_bufs[%d][%d]: NULL aligned pointer in non NULL QMP_mem_t struct \n", pm, mu);
	exit(-1);
      }

#endif	

    }
  }
#endif //USE_COMMS_SPI



/*----------------------------------------------------------------------*/
/* Build the pointer table                                              */
/*----------------------------------------------------------------------*/
  pointers_init();
  
/*----------------------------------------------------------------------*/
/* Initialise the comms                                                 */
/*----------------------------------------------------------------------*/

  comm_init();

}

void wfm::end (void)
{

  int pm, mu;

  comm_end();
  pointers_end();

  for(pm =0; pm< 2;pm++) {
#ifdef USE_COMMS_SPI
	FREE(send_bufs[pm][0]) ;
	FREE(recv_bufs[pm][0]) ;
#else
    for(mu =0; mu<4 ;mu++) {

	// This is a legit thing to do since we allocate
	// things even if we don't use them. 
#ifndef USE_COMMS_QMP
	FREE(send_bufs[pm][mu]) ;
	FREE(recv_bufs[pm][mu]) ;
#else
	QMP_free_memory(send_bufs_mem_t[pm][mu]);
	QMP_free_memory(recv_bufs_mem_t[pm][mu]);
#endif
    }
#endif
  }

#ifndef USE_COMMS_QMP  
  FREE(two_spinor);
#else
  QMP_free_memory(two_spinor_mem_t);
#endif
  FREE(spinor_tmp);
}

