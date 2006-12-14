#include<config.h>
CPS_START_NAMESPACE

/****************************************************************************/
/* 10/7/96                                                                  */
/*                                                                          */
/* wfm_sublatt_pointers:                                                    */
/*                                                                          */
/* This routine sets the arrays used in addressing by the routines in       */
/* wfm_mdagm.                                                               */
/*                                                                          */
/****************************************************************************/

CPS_END_NAMESPACE
#include<util/wilson.h>
#include<util/verbose.h>
CPS_START_NAMESPACE

#define EVEN(x) (~x & 1)

static int ceiling(int x)
{
   /* returns the ceiling of x/2 */
   static int result;
   if( EVEN(x) == 0) result = (x+1)/2;
   if( EVEN(x) == 1) result = x/2;
   return result;
}

static int floor(int x)
{
   /* returns the floor of x/2 */
   static int result;
   if( EVEN(x) == 0) result = (x-1)/2;
   if( EVEN(x) == 1) result = x/2; return result;
}

void wfm_sublatt_pointers(int lx, int ly, int lz, int lt, int slatt,
		       Wilson *wilson_p)
{
  char *cname = " ";
  char *fname = "wfm_sublatt_pointers(...)";
  VRB.Func(cname,fname);

   static int evenlx, evenp;
   static int sx;
   static int chkbd, bpad, shift;
   static int schkbd, SFT, ixx;
   static int acc0, acc1, acc2, acc3, xmax;
   static int it, iz, iy, ix;
   static int ip;
   static int *ptr;

   /* Local copy of the pointer ptr in the Wilson structure  */
   ptr = wilson_p->ptr;

   /* find if lx, ly, lz, lt are even (EVEN=1)*/
   evenlx = EVEN(lx);
   
   sx = (lx + (1 - evenlx))/2;

   /* In the structure PTRTAB initialize the block size, stride, 
      and number of blocks (used for communications). Also in the
      same structure initialize the word that results from packaging 
      them in one word ready to be stored in the SCU DMA registers. 
      This is reduntant information but it is available in convinient
      forms. 
      blklen below indicates the block length as the number of 
      half spinor words = BLOCK = 12 */
   /* X-direction */
   wilson_p->comm_offset[0] = sx * BLOCK;
   wilson_p->comm_stride[0] = sx * BLOCK + 1;
   wilson_p->comm_numblk[0] = ly * lz * lt;
   wilson_p->comm_blklen[0] = 1;
   /* package in one 32 bit word: blklen[10bits]-numblk[10bits]-stride[12bits] */
   wilson_p->comm[0].stride = sx * BLOCK + 1;
   wilson_p->comm[0].numblk = ly * lz * lt;
   wilson_p->comm[0].blklen = 1;

   /* Y-direction */
   wilson_p->comm_offset[1] = sx * ly * BLOCK;
   wilson_p->comm_stride[1] = sx * ly * BLOCK + 1;
   wilson_p->comm_numblk[1] = lz * lt;
   wilson_p->comm_blklen[1] = sx;
   /* package in one 32 bit word: blklen[10bits]-numblk[10bits]-stride[12bits] */
   wilson_p->comm[1].stride = sx * ly * BLOCK + 1;
   wilson_p->comm[1].numblk = lz * lt;
   wilson_p->comm[1].blklen = sx;

   /* Z-direction */
   wilson_p->comm_offset[2] = sx * ly * lz * BLOCK;
   wilson_p->comm_stride[2] = sx * ly * lz * BLOCK + 1;
   wilson_p->comm_numblk[2] = lt;
   wilson_p->comm_blklen[2] = sx * ly;
   /* package in one 32 bit word: blklen[10bits]-numblk[10bits]-stride[12bits] */
   wilson_p->comm[2].stride =  sx * ly * lz * BLOCK + 1;
   wilson_p->comm[2].numblk = lt;
   wilson_p->comm[2].blklen = sx * ly;

   /* T-direction */
   wilson_p->comm_offset[3] = sx * ly * lz * lt * BLOCK;
   wilson_p->comm_stride[3] = sx * ly * lz * lt * BLOCK + 1;
   wilson_p->comm_numblk[3] = 1;
   wilson_p->comm_blklen[3] = sx * ly * lz;
   /* package in one 32 bit word: blklen[10bits]-numblk[10bits]-stride[12bits] */
   wilson_p->comm[3].stride = sx * ly * lz * lt * BLOCK + 1;
   wilson_p->comm[3].numblk = 1;
   wilson_p->comm[3].blklen = sx * ly * lz;


   /************************* INCORRECT HACK ********************************/
   /* compute the subgrd volume of each chkbd */
   wilson_p->vol[0] = sx * ly * lz * lt;
   wilson_p->vol[1] = sx * ly * lz * lt;
   /************************* INCORRECT HACK ********************************/

   /* compute the subgrid volume of each padded temporary in dslash */
   wilson_p->padded_subgrid_vol[0] = (sx+1) * ly * lz * lt;
   wilson_p->padded_subgrid_vol[1] = sx * (ly+1) * lz * lt;
   wilson_p->padded_subgrid_vol[2] = sx * ly * (lz+1) * lt;
   wilson_p->padded_subgrid_vol[3] = sx * ly * lz * (lt+1);

   /* calculate the external variables YZTMAX, OFFSET needed
      for the sublattice loop */
   wilson_p->yztmax = ly * lz * lt - 1;
   wilson_p->offset = ly * lz * lt * 5;

   /* calculate the pointers for the various cases: 
      chkbd=0, 1, bpad=0, 1, shift=0, 1 */ 

   ip = 0;
   for( chkbd=0 ; chkbd <=1; ++chkbd) 
      for(bpad=0 ; bpad <=1; ++bpad) 	 
	 for( shift=0 ; shift <=1; ++shift){
/*
	    VRB.Flow(cname,fname,"\n %d %d %d \n", chkbd, bpad, shift );
	    VRB.Flow(cname,fname,"---------------------------------- \n");
*/

	    /* If we use the coordinate system of the whole lattice then
	       let the coordinates of a point be (X, Y, Z, T). If we use the
	       coordinate system of the given sublattice then let the
	       coordinates of a point on that sublattice be (x, y, z, t). The 
	       formulas use (x, y, z, t) and the parity is defined as 
	       p=(-1)^(x+y+z+t). The true parity is P=(-1)^(X+Y+Z+T).
	       If for (x, y, z, t)=(0, 0, 0, 0) P is not equal to 1 then the 
	       meaning of (even/odd) checkerboard is reversed. 
	       If checkerboard is even chkbd=0.
	       If checkerboard is odd  chkbd=1.
	       If sublattice is even slatt=0 (P=1).
	       If sublattice is odd  slatt=1 (P=-1).
	       Then the "sublattice checkerboard" is:
	       schkbd = slatt .XOR. chkbd */

	    schkbd = (slatt | chkbd) & (1 - (slatt & chkbd));
	    	    
	    /* loop over y, z, t */
	    for( it=0 ; it < lt; ++it){
 	       for( iz=0 ; iz < lz; ++iz){
		  for( iy=0 ; iy < ly; ++iy){

		     /* set the relevant coordinate shifts */
		     if( (bpad == 0) && (shift == 0) ) SFT = 0;
		     if( (bpad == 0) && (shift == 1) ) SFT = 1;
		     if( (bpad == 1) && (shift == 0) ) SFT = 0;
		     if( (bpad == 1) && (shift == 1) ) SFT = -1;

		     /* find the parity of the x-coordinate of the first
			point of the x-line (after any shifting)*/
		     evenp = EVEN(SFT + 0 + iy + iz + it);

		     /* if this point is a point of the given "sublattice
			checkerboard", then ix=0, else ix=1 */
		     if( schkbd == (1 - evenp)) ix = 0, xmax = sx-1;
		     if( schkbd == evenp) ix = 1, xmax = sx-1-(1 - evenlx);

		     /* for shifts along the x (0) direction, 
			the value of the x-coordinate of the first
			point of the x-line of the given "sublattice
			checkerboard" has x-coordinate ix = ix + SFT */
		     if( evenlx == 0 && bpad == 0 ) ixx = ceiling(ix+SFT);
		     if( evenlx == 0 && bpad == 1 ) ixx = floor(ix+SFT);
		     if( evenlx == 1 ) ixx = floor(ix+SFT);
		     acc0 = bpad + ixx + (sx+1) * (iy + ly * ( iz + lz * it));

		     /* for shifts along the y, z, t (1, 2, 3) directions this
			point has x-coordinate ix */
		     if( evenlx == 0 && bpad == 0 ) ixx = ceiling(ix);
		     if( evenlx == 0 && bpad == 1 ) ixx = floor(ix);
		     if( evenlx == 1 ) ixx = floor(ix);
		     acc1=ixx+(sx)*((iy+bpad+SFT) + (ly+1) * ( iz + lz * it));
		     acc2=ixx+(sx)*(iy + ly * ( (iz+bpad+SFT) + (lz+1) * it));
		     acc3=ixx+(sx)*(iy + ly * ( iz + lz * (it+bpad+SFT)));

		     /* set the pointers */
		     ptr[ip++]=BLOCK * acc0;
		     ptr[ip++]=BLOCK * acc1;
		     ptr[ip++]=BLOCK * acc2;
		     ptr[ip++]=BLOCK * acc3;
		     ptr[ip++]=xmax;
/*
		     VRB.Flow(cname,fname, "%4d %4d %4d %4d %4d \n", 
			      ptr[ip-5], ptr[ip-4], ptr[ip-3], 
			      ptr[ip-2], ptr[ip-1] );
*/
		  }
	       }
	    }
	 }

   return;
}
















CPS_END_NAMESPACE
