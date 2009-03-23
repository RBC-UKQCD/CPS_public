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
#include <string.h>
#include <stdio.h>

#include <util/wfm.h>
CPS_START_NAMESPACE

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
  int mu;

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
    printf("wfm::init Even x logic bomb\n");
    exit(-1);
  }
  if ( nbound[1] * sly * 2 != (slx*sly*slz*slt) ) {
    printf("wfm::init Even y logic bomb\n");
    exit(-1);
  }
  if ( nbound[2] * slz * 2 != (slx*sly*slz*slt) ) {
    printf("wfm::init Even z logic bomb\n");
    exit(-1);
  }
  if ( nbound[3] * slt * 2 != (slx*sly*slz*slt) ) {
    printf("wfm::init Even t logic bomb\n");
    exit(-1);
  }

/*--------------------------------------------------------------------------*/
/* Reserve memory for 1  temporary spinor (needed by mdagm)                 */
/*--------------------------------------------------------------------------*/
  spinor_words = SPINOR_SIZE * vol;

  spinor_tmp = (Float *) ALLOC(spinor_words*sizeof(Float)*2);
  if(spinor_tmp == 0)
  	spinor_tmp = (Float *) qalloc(QCOMMS,spinor_words*sizeof(Float)*2);
  if(spinor_tmp == 0){
    printf("wfm::spinor_tmp allocate\n");
    exit(-1);
  }  
  bzero((char *)spinor_tmp,spinor_words*sizeof(Float)*2);
//  printf("spinor_tmp is %x\n",spinor_tmp);


//~~
//~~ twisted mass fermions:  sets WilsonArg.spinor_tmp tp 
//~~ address of temporary spinor in wfm class
//~~    
  wilson_p->spinor_tmp = spinor_tmp;
//~~

/*--------------------------------------------------------------------------*/
/* Reserve memory for the 4 forward and 4 backward spin projected half      */ 
/* spinors.                                                                 */
/*--------------------------------------------------------------------------*/


  /*PAB 10/1/2001 */
  half_spinor_words = NMinusPlus * ND * PAD_HALF_SPINOR_SIZE * vol;
  two_spinor = (Float *)ALLOC(half_spinor_words*sizeof(Float));
  if(two_spinor == 0)
  	two_spinor = (Float *)qalloc(QCOMMS,half_spinor_words*sizeof(Float));
  if(two_spinor == 0){
    printf("wfm::two_spinor allocate\n");
    exit(-1);
  } 
    
  for ( int pm = 0;pm<2;pm++ ) {
    for ( mu = 0 ; mu < 4 ; mu ++) {
      half_spinor_words = PAD_HALF_SPINOR_SIZE * nbound[mu];
      recv_bufs[pm][mu] = (Float *)ALLOC(half_spinor_words*sizeof(Float));
      if(recv_bufs[pm][mu] == 0){
	printf("wfm::recv_bufs allocate\n");
	exit(-1);
      }
      send_bufs[pm][mu]=(Float *)SEND_ALLOC(half_spinor_words*sizeof(Float));
      if(send_bufs[pm][mu] == 0){
	printf("wfm::send_bufs allocate\n");
	exit(-1);
      }
    }
  }


/*----------------------------------------------------------------------*/
/* Build the pointer table                                              */
/*----------------------------------------------------------------------*/
  pointers_init();
//  printf("wfm::Pointers initialised\n");
  
/*----------------------------------------------------------------------*/
/* Initialise the comms                                                 */
/*----------------------------------------------------------------------*/

  comm_init();
//  printf("wfm::comms initialised\n");

}
void wfm::end (void)
{

  comm_end();
  pointers_end();
      
  for(int pm =0;pm<2;pm++) {
    for(int mu =0;mu<4;mu++) {
      FREE(send_bufs[pm][mu]) ;
      FREE(recv_bufs[pm][mu]) ;
    }
  }
  FREE(two_spinor);
  FREE(spinor_tmp);
}

CPS_END_NAMESPACE
