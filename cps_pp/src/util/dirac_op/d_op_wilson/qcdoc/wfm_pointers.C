/****************************************************************************/
/* PAB 10/1/2001                                                            */
/* QCDOC variant                                                            */
/*                                                                          */
/* wfm_pointers:                                                            */
/*                                                                          */
/* This routine sets the arrays used in addressing by the routines in       */
/* wilson matrix multiply                                                   */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* TwoSpinor *shift_table[chkbd][npsite][NMinusPlus][Nmu]                   */
/*                                                                          */
/*     Maps a site "psite" in chkbd for a 4spinor to an pointer to          */
/*     the site in our private 2spinor-mu array corresponding to the        */
/*     site in NOT(chkbd) that is:                                          */
/*        shifted by -1 in dimension mu when minusplus = 0                  */
/*        shifted by +1 in dimension mu when minusplus = 1                  */
/*                                                                          */
/* shift_table points directly to a location in the -+mu send comms         */
/* buffer if site is on the boundary and we are non-local in dimension mu   */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* TwoSpinor *face_table[chkbd][minusplus][mu][bufsite]                     */
/*                                                                          */
/* Maps a linear comms buffer (indexed bufsite) to the face                 */
/* in a shift in direction +/- mu for  mu-2spinor                           */
/* See wfm_comm.C for usage                                                 */
/*                                                                          */
/****************************************************************************/

#include <stdlib.h>
#include <stdio.h>

#include <util/wfm.h>
CPS_START_NAMESPACE

void wfm::pointers_end(void)
{
  for ( int cb = 0;cb<Ncb;cb++ ) 
    FREE(shift_table[cb]);
  for(int cb=0;cb<2;cb++) {
    for(int pm=0;pm<2;pm++) {
      for(int mu=0; mu<4;mu++) {
	FREE(face_table[cb][pm][mu]);
      }
    }
  }
}

void wfm::pointers_init(void)
{
  int mu;
  int cb, pm;
  int shift;
  int x,y,z,t;
  int lx,ly,lz,lt;
 
  int bound_index[Ncb][NMinusPlus][ND];
  int local_p,shift_p;
  int local_addr[ND];
  int shift_addr[ND];

  int psite,shift_psite;
  int offset,tab,table_size;

  lx = local_latt[0];
  ly = local_latt[1];
  lz = local_latt[2];
  lt = local_latt[3];

  for ( cb = 0 ; cb<2 ; cb++) {
    //ND * plus_minus * hvol * parities
    table_size = NMinusPlus*vol*Nmu;
    shift_table[cb]=(unsigned long *)ALLOC( table_size*sizeof(unsigned long));
    if (!shift_table[cb]){ 
      printf("shift_table=%p\n",shift_table[cb]);exit(-1);}
  }

  for(cb=0;cb<2;cb++){
    for(pm=0;pm<2;pm++){
      for(mu=0; mu<4;mu++){
	face_table[cb][pm][mu]= (unsigned long *)
	  ALLOC(allbound*sizeof(unsigned long));
    if (!face_table[cb][pm][mu]){
      printf("face_table=%p\n",face_table[cb][pm][mu]);exit(-1);}
      }
    }
  }
  
  for( mu = 0 ; mu < Nmu ; mu++ ) {
    bound_index[0][0][mu] = 0;
    bound_index[0][1][mu] = 0;
    bound_index[1][0][mu] = 0;
    bound_index[1][1][mu] = 0;
  }

  /*
   * For now point send buffers beyond end of body.
   * wfm::scale_ptr will remap this to the real send buffer, once
   * it gets assigned with a pointer value.
   */
  send_offset[Minus][0]  = 8* vol; /*We interleave the 8 2spinors*/
  send_offset[Minus][1]  = send_offset[Minus][0] + nbound[0];
  send_offset[Minus][2]  = send_offset[Minus][1] + nbound[1];
  send_offset[Minus][3]  = send_offset[Minus][2] + nbound[2];

  send_offset[Plus][ 0 ] = send_offset[Minus][0] + allbound;
  send_offset[Plus][ 1 ] = send_offset[Minus][1] + allbound;
  send_offset[Plus][ 2 ] = send_offset[Minus][2] + allbound;
  send_offset[Plus][ 3 ] = send_offset[Minus][3] + allbound;

  for ( pm = 0; pm<2; pm++ ) {

    if ( pm == Plus ) { /*forwards or backwards*/
      shift = 1;
    } else { 
      shift = -1;
    }

      /*Loop naively over local lattice*/
    for ( t=0 ; t<lt ; t++ ) {
     local_addr[3] = t;
     for ( z=0 ; z<lz ; z++ ) {
      local_addr[2] = z;
      for ( y=0 ; y<ly ; y++ ) {
       local_addr[1] = y;
       for ( x=0 ; x<lx ; x++ ) {
        local_addr[0] = x;

        local_p = wfm::local_parity(x,y,z,t);
        psite   = wfm::local_psite(local_addr,local_latt);

        for ( mu = 0 ; mu < 4; mu ++) {

          shift_addr[0] = x;
          shift_addr[1] = y;
          shift_addr[2] = z;
          shift_addr[3] = t;

          shift_addr[mu] += shift;
          shift_p = 1-local_p ;

          /*
	   * Offset to the source site in the 2 spinor array is common
	   */
	  tab = wfm::interleave_site(pm,mu,psite);

          /*
	   * Get the offset to the destination site in the 2 spinor array
	   * If destination is in interior, trivial.
	   * Also implement periodic wrap if local_comm[mu].
	   */


	  if (((shift_addr[mu] >= 0 ) && (shift_addr[mu]<local_latt[mu]) )
             || local_comm[mu] ) {

	    /*Local periodicity, does nothing if interior*/
	    shift_addr[mu]= (shift_addr[mu] + local_latt[mu])%local_latt[mu];
	    shift_psite   = local_psite(shift_addr,local_latt);
	    offset = interleave_site(pm,mu,shift_psite) ;
	    shift_table[local_p][tab] = ((unsigned long)two_spinor +
                   offset * PAD_HALF_SPINOR_SIZE * sizeof(Float));

	  } else { /*non-local and we're on the boundary*/


	    /*Local periodicity, does nothing if interior*/
	    shift_addr[mu]= (shift_addr[mu] + local_latt[mu])%local_latt[mu];
	    shift_psite   = local_psite(shift_addr,local_latt);
	    offset = interleave_site(pm,mu,shift_psite) ;
	    /*
	     * The minus face receives something written by the 
	     * plus face. Thus the face table should contain the pointer
	     * to the site we would have sent to in the above case.
	     */
	    face_table[local_p][pm][mu][bound_index[local_p][pm][mu]]
	      = offset;

	    /*
	     * And we use face_table[local_p][pm][mu] in conjunction
	     * with the send buffer[pm][mu] and source parity local_p
	     */
	    offset = bound_index[local_p][pm][mu];

	    shift_table[local_p][tab] = (unsigned long) send_bufs[pm][mu]
	      + offset * PAD_HALF_SPINOR_SIZE * sizeof(Float);

	    bound_index[local_p][pm][mu] ++;
	  }


         
        } /*mu*/
       }/*x*/
      }/*y*/
     }/*z*/
    }/*t*/
  }/*pm*/

  //  int bound_index[Ncb][NMinusPlus][ND];
  for( cb = 0 ; cb < 2 ; cb++ ) {
  for( pm = 0 ; pm < 2 ; pm++ ) {
  for( mu = 0 ; mu < Nmu ; mu++ ) {
    if ( !local_comm[mu] ) { 
      if ( bound_index[cb][pm][mu] != nbound[mu] ) {
	printf("Boundary size mismatch[%d][%d][mu=%d] : %d != %d \b",
	       cb,pm,mu,bound_index[cb][local_p][mu],nbound[mu]);
	exit(-1);
      }
    } else { 
      if ( bound_index[cb][pm][mu] != 0 ) {
	printf("Boundary size mismatch[%d][%d] mu=%d\b",
	       cb,pm,mu);
	exit(-1);
      }
    }
  }
  }
  }
//  scale_ptr();
  return;
}

int wfm::interleave_site(int mp,int mu, int site)
{
  return(mu+mp*ND+site*NMinusPlus*ND);
}

int wfm::local_parity(int x,int y,int z,int t)
{
  return((x+y+z+t) % 2);
}

int wfm::local_psite(int addr[4],int latt[4])
{
  int site;
  site =  addr[0] + latt[0]*(addr[1]
                             + latt[1]*(addr[2]
                                        + latt[2]*addr[3]));
  return(site/2);
}



CPS_END_NAMESPACE
