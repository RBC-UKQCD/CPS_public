/****************************************************************************/
/* 10/6/96                                                                  */
/*                                                                          */
/* wilson_init:                                                             */
/*                                                                          */
/* This routine performs all initializations needed before wilson func      */
/* are called. It sets the addressing related arrays and reserves memory    */
/* for the needed temporary buffers. It only needs to be called             */
/* once at the begining of the program (or after a wilson_end call)         */
/* before any number of calls to wilson funcs are made.                     */
/*                                                                          */
/* PAB: 10/1/2001 Updated for new scheme for QCDOC                          */
/* We lay out the 2-spinors in a very different way, but keep the 4spinor   */
/* and gauge layout the same                                                */
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
/*         In hindsight I'm a little surprised this isn't used on QCDSP.    */
/*                                                                          */
/****************************************************************************/

/*--------------------------------------------------------------------------*/
/* Include header files                                                     */
/*--------------------------------------------------------------------------*/
#include <stdlib.h>
#include <stdio.h>
#include <emalloc.h>
#include <qcdoc_align.h>
#include <util/verbose.h>
#include <util/gjp.h>

#define QCDOC_FORCE_PEC_ALIGN
#include <util/wilson.h>
#include "prec.h"

#define MAX_VOLUME (4*2*2*2)
/*
 * The send buffers - since we use storegathering, we don't need to cache
 */
Float send_comms_buf[8*PAD_HALF_SPINOR_SIZE*MAX_VOLUME] LOCATE("Edramnoncache");

/*
 * Receive buffers - since these are read we need to cache for performance
 * however, by placing in transient we can guarantee they are evicted before
 * reuse even on 2^4, and avoid cache flushing entirely.
 */
Float recv_comms_buf[8*PAD_HALF_SPINOR_SIZE*MAX_VOLUME] LOCATE("Edramtransient");

/*--------------------------------------------------------------------------*/
/* Function declarations                                                    */
/*--------------------------------------------------------------------------*/
void wfm_sublatt_pointers(int slx, 
			  int sly, 
			  int slz, 
			  int slt, 
			  int slatt_eo, 
			  Wilson *wilson_p);

/*=========================================================================*/
/* wilson_init:                                                            */
/*=========================================================================*/

void wilson_init(Wilson *wilson_p)  /* pointer to Wilson type structure    */
{
  char *cname = " ";
  char *fname = "wilson_init(Wilson*)";
  Float * comm_p ;
  VRB.Func(cname,fname);

  int spinor_words;             /* size of the spinor field on the         */
				/* sublattice checkerboard                 */

  int half_spinor_words;        /* size of the spin-projected "half_spinors*/
                                /* on the sublattice checkerboard including*/
                                /* the communications padding              */

  int slx;                          /* x-direction size of node sublattice */
  int sly;                          /* y-direction size of node sublattice */
  int slz;                          /* z-direction size of node sublattice */
  int slt;                          /* t-direction size of node sublattice */
  int slatt_eo;                     /* =0/1 if the sublattice is even/odd. */
  int i;
  int mu;
  int size;

/*--------------------------------------------------------------------------*/
/* Set sublattice direction sizes                                           */
/*--------------------------------------------------------------------------*/
  slx = wilson_p->local_latt[0] = GJP.XnodeSites();
  sly = wilson_p->local_latt[1] = GJP.YnodeSites();
  slz = wilson_p->local_latt[2] = GJP.ZnodeSites();
  slt = wilson_p->local_latt[3] = GJP.TnodeSites();

/*--------------------------------------------------------------------------*/
/* Determine if the sublattice is even or odd from the node coordinates     */
/* If (px,py,pz,pt) are the coordinates of the node and the node            */
/* mesh is of size (nx,ny,nz,nt) then a node is even/odd if its             */
/* lexicographical number =  px + nx * ( py + ny * ( pz + nz * ( pt )))     */
/* is even/odd.                                                             */
/*--------------------------------------------------------------------------*/
/* A runtime system function is needed here to determine (px,py,pz,pt) and  */
/* (nx,ny,nz,nt). For now we set slat_eo = 0 which is a safe choice if      */
/* slx,sly,slz,slt are all even.                                            */
/* ??? */
  slatt_eo = 0;


/*-----------------------------------------------------------------------*/
/* compute the subgrd volume of each chkbd ... at least two local dims   */
/* must be even for this code to be correct.                             */
/*-----------------------------------------------------------------------*/
  wilson_p->vol[0] = (slx * sly * slz * slt)/2;
  wilson_p->vol[1] = (slx * sly * slz * slt)/2;
  
  wilson_p->nbound[0] = (sly * slz * slt)/2; 
  wilson_p->nbound[1] = (slx * slz * slt)/2;
  wilson_p->nbound[2] = (slx * sly * slt)/2;
  wilson_p->nbound[3] = (slx * sly * slz)/2;
  wilson_p->allbound  = wilson_p->nbound[0]
                      + wilson_p->nbound[1]
                      + wilson_p->nbound[2]
                      + wilson_p->nbound[3];

  if ( wilson_p->nbound[0] * slx * 2 != (slx*sly*slz*slt) ) {
    printf("Even x logic bomb\n");
    exit(-1);
  }
  if ( wilson_p->nbound[1] * sly * 2 != (slx*sly*slz*slt) ) {
    printf("Even y logic bomb\n");
    exit(-1);
  }
  if ( wilson_p->nbound[2] * slz * 2 != (slx*sly*slz*slt) ) {
    printf("Even z logic bomb\n");
    exit(-1);
  }
  if ( wilson_p->nbound[3] * slt * 2 != (slx*sly*slz*slt) ) {
    printf("Even t logic bomb\n");
    exit(-1);
  }

/*--------------------------------------------------------------------------*/
/* Reserve memory for 1  temporary spinor (needed by mdagm)                 */
/*--------------------------------------------------------------------------*/
  spinor_words = SPINOR_SIZE * wilson_p->vol[0];

  wilson_p->spinor_tmp = (Float *) emalloc(spinor_words*sizeof(Float)*2);
  if(wilson_p->spinor_tmp == 0){
    printf("wilson_p->spinor_tmp allocate\n");
    exit(-1);
  }
    

/*--------------------------------------------------------------------------*/
/* Reserve memory for the 4 forward and 4 backward spin projected half      */ 
/* spinors.                                                                 */
/*--------------------------------------------------------------------------*/


  /*PAB 10/1/2001 */
  half_spinor_words = 4 * PAD_HALF_SPINOR_SIZE * wilson_p->vol[0];

  wilson_p->af = (Float *) emalloc(half_spinor_words*sizeof(Float));
  wilson_p->ab = (Float *) emalloc(half_spinor_words*sizeof(Float));
  if(wilson_p->af == 0){
    printf("wilson_p->af allocate\n");
    exit(-1);
  }

  if(wilson_p->ab == 0){
    printf("wilson_p->ab allocate\n");
    exit(-1);
  }

  comm_p = recv_comms_buf;
  for ( mu = 0 ; mu < 4 ; mu ++) {

    half_spinor_words = PAD_HALF_SPINOR_SIZE * wilson_p->nbound[mu];

    wilson_p->recv_f[mu] = comm_p;
    comm_p += half_spinor_words;
  }
  for ( mu = 0 ; mu < 4 ; mu ++) {

    half_spinor_words = PAD_HALF_SPINOR_SIZE * wilson_p->nbound[mu];

    wilson_p->recv_b[mu] = comm_p;
    comm_p += half_spinor_words;

  }

  comm_p = send_comms_buf;
  for ( mu = 0 ; mu < 4 ; mu ++) {

    half_spinor_words = PAD_HALF_SPINOR_SIZE * wilson_p->nbound[mu];

    wilson_p->send_f[mu] = comm_p;
    comm_p += half_spinor_words;

  }
  /*---------------------------------------------------------------
   *Put all receive forward recv buffers contiguously, and
   *then all backwards recv buffers contiguously. 
   *This eliminates PEC register misses during face scattering
   *---------------------------------------------------------------
   */
  for ( mu = 0 ; mu < 4 ; mu ++) {

    wilson_p->send_b[mu] = comm_p;
    comm_p += half_spinor_words;


  }

/*----------------------------------------------------------------------*/
/* Build the pointer table                                              */
/*----------------------------------------------------------------------*/

  wfm_sublatt_pointers(slx,sly,slz,slt, slatt_eo, wilson_p);
  
/*----------------------------------------------------------------------*/
/* Initialise the comms                                                 */
/*----------------------------------------------------------------------*/

  wfm_comm_init(wilson_p);

}
























