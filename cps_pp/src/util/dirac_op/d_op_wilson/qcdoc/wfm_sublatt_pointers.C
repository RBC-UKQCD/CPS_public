/****************************************************************************/
/* 10/7/96                                                                  */
/*                                                                          */
/* wfm_sublatt_pointers:                                                    */
/*                                                                          */
/* This routine sets the arrays used in addressing by the routines in       */
/* wfm_mdagm.                                                               */
/*                                                                          */
/* PAB 10/1/2001                                                            */
/* QCDOC variant                                                            */
/*                                                                          */
/* shift_table[chkbd][plusminus][npsite][4]                                 */
/* Maps a site in chkbd to an offset to a site in NOT(chkbd) or to an offset*/
/* to the comms buffer.                                                     */
/*                                                                          */
/* face_table[chkbd][plusminus][mu][bufsite]                                */
/* Maps a linear comms buffer (indexed bufsite) to the (low/received) face  */
/* in a shift in direction plusminus                                        */
/*                                                                          */
/* plusminus=0 => positive mu shift                                         */
/* plusminus=1 => negative mu shift                                         */
/*                                                                          */
/****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <util/wilson.h>
#include "config.h"
#include "prec.h"

int  wfm_local_parity(int x,int y,int z,int t);
int  wfm_local_psite(int addr[4],int latt[4]);
int  wfm_interleave_site(int mu, int site);


void wfm_sublatt_pointers(int lx, int ly, int lz, int lt, int slatt,
		       Wilson *wp)
{
  char *cname = " ";
  char *fname = "wfm_sublatt_pointers(...)";

  int sx, mu;
  int cb, pm;
  int shift;
  int x,y,z,t;
 
  int bound_index_plus[2][4];
  int bound_index_minus[2][4];

  int local_p,shift_p;
  int local_addr[4];
  int shift_addr[4];
  
  int psite,shift_psite;
  int table_size,offset,tab;


  table_size = 2*2 * 4*wp->vol[0];
  table_size = table_size + 2*2*wp->allbound;

  wp->ptr = (unsigned long *)malloc( table_size * sizeof(unsigned long));

  if( wp->ptr == 0){
    printf("wp->ptr\n");
    exit(-1);
  }

  for ( cb = 0 ; cb<2 ; cb++) {
    for ( pm = 0 ; pm<2 ; pm++) {
      wp->shift_table[cb][pm] = &wp->ptr[(2*cb + pm)*wp->vol[0]*4];
    }
  }

  wp->face_table[0][0][0] = & wp->ptr[4*wp->vol[0]*4];
  wp->face_table[0][0][1] = & wp->face_table[0][0][0][wp->nbound[0]];
  wp->face_table[0][0][2] = & wp->face_table[0][0][1][wp->nbound[1]];
  wp->face_table[0][0][3] = & wp->face_table[0][0][2][wp->nbound[2]];

  for(mu=0; mu<4;mu++){
    wp->face_table[0][1][mu] = & wp->face_table[0][0][mu][wp->allbound];
    wp->face_table[1][0][mu] = & wp->face_table[0][1][mu][wp->allbound];
    wp->face_table[1][1][mu] = & wp->face_table[1][0][mu][wp->allbound];
  }

  wp->local_latt[0] = lx;
  wp->local_latt[1] = ly;
  wp->local_latt[2] = lz;
  wp->local_latt[3] = lt;
  
  for( mu = 0 ; mu < 4 ; mu++ ) {
    bound_index_plus[0][mu] = 0;
    bound_index_minus[0][mu] = 0;
    bound_index_plus[1][mu] = 0;
    bound_index_minus[1][mu] = 0;
    wp->local_comm[mu] = 0;
  }

  for ( pm = 0; pm<2; pm++ ) {

    if ( pm == 0 ) { /*forwards or backwards*/
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

        local_p = wfm_local_parity(x,y,z,t);
	psite   = wfm_local_psite(local_addr,wp->local_latt);

	for ( mu = 0 ; mu < 4; mu ++) {
	  shift_addr[0] = x;
	  shift_addr[1] = y;
	  shift_addr[2] = z;
	  shift_addr[3] = t;

	  shift_addr[mu] += shift;
	  shift_p = 1-local_p ; /*If global_latt[mu] == 1 this is wrong*/
	  
	  /*Minus boundary*/
	  if (  shift_addr[mu] < 0 ) {
	    if ( wp->local_comm[mu] ) {

	      shift_addr[mu]= wp->local_latt[mu] - 1; /*Local periodicity*/
	      shift_psite   = wfm_local_psite(shift_addr,wp->local_latt);
	      offset = wfm_interleave_site(mu,shift_psite) ;
	      tab = wfm_interleave_site(mu,psite);

	      wp->shift_table[local_p][pm][tab] = 
		((unsigned long)wp->ab) 
		+ PAD_HALF_SPINOR_SIZE * sizeof(Float) * offset ; 

	    } else {

	      offset = bound_index_minus[local_p][mu];
	      tab = wfm_interleave_site(mu,psite);
	      wp->shift_table[local_p][pm][tab] = 
		((unsigned long)wp->send_b[mu]) 
		+ PAD_HALF_SPINOR_SIZE * sizeof(Float) * offset ; 
	      
	      /*
	       *  We have a point on the boundary, so can set
	       *  an entry in the face_table for the opposite direction
	       */
	      wp->face_table[local_p][1-pm][mu][bound_index_minus[local_p][mu]]
		= wfm_interleave_site(mu,psite);

	      bound_index_minus[local_p][mu] ++;
	    }

	  /*Plus boundary*/
	  } else if ( shift_addr[mu] >= wp->local_latt[mu] ) {

	    if ( wp->local_comm[mu] ) { 

	      shift_addr[mu]= 0; /*Local periodicity*/
	      shift_psite = wfm_local_psite(shift_addr,wp->local_latt);
	      offset = wfm_interleave_site(mu,shift_psite);
	      tab = wfm_interleave_site(mu,psite);

	      wp->shift_table[local_p][pm][tab] = 
		((unsigned long)wp->af)
		+ PAD_HALF_SPINOR_SIZE*sizeof(Float)*offset ; 

	    } else {

	      offset = bound_index_plus[local_p][mu];
	      tab = wfm_interleave_site(mu,psite);
	      wp->shift_table[local_p][pm][tab] = 
		  ((unsigned long)wp->send_f[mu])
		+ PAD_HALF_SPINOR_SIZE*sizeof(Float)*offset ; 

	      /*
	       *  We have a point on the boundary, so can set
	       *  an entry in the face_table
	       */
	      wp->face_table[local_p][1-pm][mu][bound_index_plus[local_p][mu]]
		                     = wfm_interleave_site(mu,psite);

	      bound_index_plus[local_p][mu] ++;
	    }
	  /*Interior*/
	  } else { 

	      shift_psite = wfm_local_psite(shift_addr,wp->local_latt);
	      offset = wfm_interleave_site(mu,shift_psite);
	      tab = wfm_interleave_site(mu,psite);
	      if ( pm == 0 ) {
		wp->shift_table[local_p][pm][tab] = 
		  ((unsigned long)wp->af)+ PAD_HALF_SPINOR_SIZE*sizeof(Float)*offset ; 
	      } else {
		wp->shift_table[local_p][pm][tab] = 
		  ((unsigned long)wp->ab)+ PAD_HALF_SPINOR_SIZE*sizeof(Float)*offset ; 
	      }

	  }/*bound/interior*/
	 
	} /*mu*/
       }/*x*/
      }/*y*/
     }/*z*/
    }/*t*/
  }/*pm*/


  return;
}


int wfm_interleave_site(int mu, int site)
{
  return(mu+site*4);
}

/*
 * int global_parity(int x,int y,int z,int t);
 * Need access to global lattice coordinates for this... how???
 */

int wfm_local_parity(int x,int y,int z,int t)
{
  return((x+y+z+t) % 2);
}

int wfm_local_psite(int addr[4],int latt[4])
{
  int site;
  site =  addr[0] + latt[0]*(addr[1]
			     + latt[1]*(addr[2]
					+ latt[2]*addr[3]));
  return(site/2);
}




