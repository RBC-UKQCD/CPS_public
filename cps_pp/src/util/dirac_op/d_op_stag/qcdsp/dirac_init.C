include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-04-30 12:18:00 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_stag/qcdsp/dirac_init.C,v 1.3 2004-04-30 12:18:00 zs Exp $
//  $Id: dirac_init.C,v 1.3 2004-04-30 12:18:00 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: dirac_init.C,v $
//  $Revision: 1.3 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_stag/qcdsp/dirac_init.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
/***********************************************************/
/*	dirac_init.C	Version 1.1
*************************************************************
	Dong Chen
	Version 1.1	2/4/1994
	Modified	3/9/1994

	Modified by Chulwoo Jung 12/19/96

	void stag_dirac_init()
	This version includes the serial communication.

    local U storage order:	x,y,z,t
    Staggered Fermion phase:	eta_x =   1 
				eta_y = (-1)^x
				eta_z = (-1)^(x+y)
				eta_t = (-1)^(x+y+z)
    Serial Unit Control Register:	t,x,y,z
				(0,t+) (1,t-) (2,x+) (3,x-)
				(4,y+) (5,y-) (6,z+) (7,z-)
    Local X storage order:	t,x,y,z
      X_even and X_odd are stored seperately.

    Notation:
	Even or Odd address tables only relates to the
    result of dirac operator.  For example, if one
    is calculating y_e = D x_o,  then one should use
    even addr tables.  The tables contains the offsets
    of x_o neighbering y_e, relative to the base of x_o.
*/
/***********************************************************/
CPS_END_NAMESPACE
#include<util/dirac_op.h>
//#include<alg/do_arg.h>
#include<util/smalloc.h>
#include<util/gjp.h>
#include<comms/nga_reg.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------
//  RDM  11/10/97
//  Initialization for SCU sytem calls
//------------------------------------------------------------------

extern "C" void StagSCUSetup();

//------------------------------------------------------------------
//  definition of structs
//------------------------------------------------------------------

#define SURF_UX_TAB 3	/* 0 -- surface flag
			   1 -- X offset
			   2 -- U offset			*/
#define INTN_UX_TAB 9	/* 0 -- number of internal U . X's for a site
			   1,8 -- (X offset,U offset) for internal
				x,y,z,t with a maximum 4 pairs 	*/
#define UDAGX_TAB 6	/* 0 -- surface flag
			   1 -- U offset
			   2,5 -- X offset in x,y,z,t order	*/

typedef struct {
  unsigned int	addr : 19;
  unsigned int  rec  :  1;
  unsigned int  send :  1;
} SCU_START_ADDR;

typedef struct {
  unsigned int  stride : 12;
  unsigned int	blk_len : 10;
  unsigned int	num_blk : 10;
} SCU_BLK_DES;

typedef struct {
  int  n_transfer[4];
  struct {
    SCU_START_ADDR  start_addr[8];
    SCU_BLK_DES   blk_des[8];
  } trans_des[2];
}  TRANSFER;

typedef struct {
  int *surf_ux;
  int *intn_ux;
  int *udagx;
  TRANSFER m_trans;
  TRANSFER p_trans;
  int p_send_len[4];
  int p_rec_len[4];
  int m_send_len[4];
  int m_rec_len[4];
} ADDRTAB;

enum{VECT_LEN=6, MATRIX_SIZE=18, SITE_LEN=72};


//------------------------------------------------------------------



#define EVEN_ADDR1(x,y,z,t) (((z)/2)*ntxyb+((z)%2)*ntxyb_e2 \
	+((y)/2)*ntxb+((y)%2)*ntxb_e2+((x)/2)*ntb+((x)%2)*ntb_e2 \
	+(t)/2)*VECT_LEN
#define ODD_ADDR1(x,y,z,t) (((z)/2)*ntxyb+((z)%2)*ntxyb_o2 \
	+((y)/2)*ntxb+((y)%2)*ntxb_o2+((x)/2)*ntb+((x)%2)*ntb_o2 \
	+(t)/2)*VECT_LEN



static void comm_init(int );


//------------------------------------------------------------------
// global variables used by "dirac_serial.asm"
//------------------------------------------------------------------
//DO_ARG do_arg;
ADDRTAB e_tab;
ADDRTAB o_tab;

IFloat *p_rec_buf;
IFloat *p_send_buf;
IFloat *m_rec_buf;
int n_even_sites;
int n_odd_sites;
int nsurf_sites[4];
int srbuf_offset[4];
int scratch_addr = CRAM_SCRATCH_ADDR;




void stag_destroy_dirac_buf()
{
    sfree(p_rec_buf);
    sfree(p_send_buf);
    sfree(m_rec_buf);
    sfree(e_tab.surf_ux);
    sfree(o_tab.surf_ux);
    sfree(e_tab.intn_ux);
    sfree(o_tab.intn_ux);
    sfree(e_tab.udagx);
    sfree(o_tab.udagx);
}
extern "C"
void stag_dirac_init(const void * gauge_field_addr)
{
  int x,y,z,t;
  int base_even;
  int base_even_save;
  int base_odd;
  int base_odd_save;
  int nxb,nyb,nzb,ntb;
  int ntxb,ntxyb;
  int ntb_e2,ntxb_e2,ntxyb_e2;
  int ntb_o2,ntxb_o2,ntxyb_o2;
  int ntb_e2_save,ntxb_e2_save,ntxyb_e2_save;
  int ntb_o2_save,ntxb_o2_save,ntxyb_o2_save;
  int xp1,yp1,zp1,tp1;
  int xm1,ym1,zm1,tm1;
  int buf_size,buf_size1;
  int *e_ux_ptr;
  int *o_ux_ptr;
  int *int_ptr;
  int *addr_ptr;
  unsigned int surf_flag;
  int i,j,k;
  int nb[4];
  int xb[4];
  int send_offset,send_blk_len,send_stride,send_num_blk;
  int rec_addr[8];
  int num,num1,num2;

  //--------------------------------------------------------------
  // initialize the do_arg object
  //--------------------------------------------------------------
//   do_arg.nxb = GJP.XnodeSites();
//   do_arg.nyb = GJP.YnodeSites();
//   do_arg.nzb = GJP.ZnodeSites();
//   do_arg.ntb = GJP.TnodeSites();
  nxb = GJP.XnodeSites();
  nyb = GJP.YnodeSites();
  nzb = GJP.ZnodeSites();
  ntb = GJP.TnodeSites();
  
//   do_arg.base_odd = 0;
  base_odd = 0;
//  do_arg.u_base_addr = (IFloat *)gauge_field_addr;
  int u_base_addr = (IFloat *)gauge_field_addr;

//  int n_sites2 = (do_arg.nxb*do_arg.nyb*do_arg.nzb*do_arg.ntb)/2;
    int n_sites2 = (nxb*nyb*nzb*ntb)/2;
  n_even_sites = n_sites2;
  n_odd_sites = n_sites2;



  //--------------------------------------------------------------
  // Dong's stuff
  //--------------------------------------------------------------
//   nxb = do_arg.nxb;
//   nyb = do_arg.nyb;
//   nzb = do_arg.nzb;
//   ntb = do_arg.ntb;
//  if (do_arg.base_odd)
  if (base_odd)  
    {
    base_even = 0;
    base_odd = 1;
    }
  else
    {
    base_even = 1;
    base_odd = 0;
    }
  ntxb = ntb*nxb;
  ntxyb = ntb*nxb*nyb;
  ntb_e2 = (ntb+base_even)/2;
  ntxb_e2 = (ntxb+base_even)/2;
  ntxyb_e2 = (ntxyb+base_even)/2;
  ntb_o2 = (ntb+base_odd)/2;
  ntxb_o2 = (ntxb+base_odd)/2;
  ntxyb_o2 = (ntxyb+base_odd)/2;

  nsurf_sites[0] = nxb*nyb*nzb;		/* t direction */
  nsurf_sites[1] = nyb*nzb*ntb;		/* x direction */
  nsurf_sites[2] = nxb*nzb*ntb;		/* y direction */
  nsurf_sites[3] = nxb*nyb*ntb;		/* z direction */
  srbuf_offset[1] = 0;
  srbuf_offset[2] = ((nsurf_sites[1]+1)/2) * VECT_LEN;
  srbuf_offset[3] = srbuf_offset[2] + ((nsurf_sites[2]+1)/2) * VECT_LEN;
  srbuf_offset[0] = srbuf_offset[3] + ((nsurf_sites[3]+1)/2) * VECT_LEN;

  buf_size = (nsurf_sites[1]+1)/2 + (nsurf_sites[2]+1)/2
	   + (nsurf_sites[3]+1)/2;
  buf_size1 = (buf_size + (nsurf_sites[0]+1)/2) * VECT_LEN;
  if ( nxb%2 == 0 )
    buf_size += nsurf_sites[0];
  else
    buf_size += (nsurf_sites[0]+1)/2;
  buf_size *= VECT_LEN;

  p_rec_buf = (IFloat *)smalloc( buf_size*sizeof(IFloat) );
  p_send_buf = (IFloat *)smalloc( buf_size1*sizeof(IFloat) );
  m_rec_buf = (IFloat *)smalloc(buf_size1*sizeof(IFloat));

  num1= (((nxb*nyb*nzb*ntb+base_odd)/2 
		- ((nxb-1)*(nyb-1)*(nzb-1)*(ntb-1)+base_odd)/2 )
		  *SURF_UX_TAB+1)*sizeof(int);
  e_tab.surf_ux = (int *)smalloc(num1);

  num2= (((nxb*nyb*nzb*ntb+base_even)/2 
		- ((nxb-1)*(nyb-1)*(nzb-1)*(ntb-1)+base_even)/2)
		  *SURF_UX_TAB+1)*sizeof(int);
  o_tab.surf_ux = (int *)smalloc(num2);

  num1 = n_even_sites*INTN_UX_TAB* sizeof(int);
  e_tab.intn_ux = (int *)smalloc(num1);
  num2 = n_odd_sites*INTN_UX_TAB* sizeof(int);
  o_tab.intn_ux = (int *)smalloc(num2);

  num = n_even_sites*UDAGX_TAB*sizeof(int);
  e_tab.udagx = (int *)smalloc(num);

  num = n_odd_sites*UDAGX_TAB*sizeof(int);
  o_tab.udagx = (int *)smalloc(num);

/*
 *  Address table for X's on x+ surface whose U . X has to be
 *  calculated on node and sent to off-node for accumulation.
 */
  e_ux_ptr = e_tab.surf_ux;
  o_ux_ptr = o_tab.surf_ux;
  for ( z=nzb-1; z>=0; z--)
    for ( y=nyb-1; y>=0; y--)
      for ( x=nxb-1; x>=0; x--)
	for ( t=ntb-1; t>=0; t--)
	  {
	  surf_flag = 0;
	  if (x == nxb-1)	surf_flag |= 0x1;
	  if (y == nyb-1)	surf_flag |= 0x2;
	  if (z == nzb-1)	surf_flag |= 0x4;
	  if (t == ntb-1)	surf_flag |= 0x8;

	  if (surf_flag)
	    {
	    if ((x+y+z+t+base_odd)%2 == 0) 
	      {
	      *o_ux_ptr++ = surf_flag;
	      *o_ux_ptr++ = EVEN_ADDR1(x,y,z,t)+BANK3_BASE;
	      *o_ux_ptr++ = (x+(y+(z+t*nzb)*nyb)*nxb)*SITE_LEN+
	      					BANK2_BASE;
	      }
	    else
	      {
	      *e_ux_ptr++ = surf_flag;
	      *e_ux_ptr++ = ODD_ADDR1(x,y,z,t)+BANK3_BASE;
	      *e_ux_ptr++ = (x+(y+(z+t*nzb)*nyb)*nxb)*SITE_LEN+
	      					BANK2_BASE;
	      }
	    }
	  }
  *e_ux_ptr = 0x10;
  *o_ux_ptr = 0x10;

/*
 *  Address table for internal U . X calculation to be accumulated
 *  on board.
 */
  e_ux_ptr = e_tab.intn_ux;
  o_ux_ptr = o_tab.intn_ux;
  for ( z=0; z<nzb; z++)
    {
    if (z>0)     zm1 = z-1;  else  zm1 = nzb-1;
    for ( y=0; y< nyb; y++)
      {
      if (y>0)     ym1 = y-1;  else  ym1 = nyb-1;
      for ( x=0; x< nxb; x++)
	{
	if (x>0)     xm1 = x-1;  else  xm1 = nxb-1;
	for ( t=0; t<ntb; t++)
	  {
	  if (t>0)     tm1 = t-1;  else  tm1 = ntb-1;

	  if ((x+y+z+t+base_odd)%2 == 0) 
	    {
  	    int_ptr = e_ux_ptr + 1;
	    *e_ux_ptr = 0;
	    if (x != 0)
	      {
	      *e_ux_ptr += 1;
	      *int_ptr++ = ODD_ADDR1(xm1,y,z,t);
	      *int_ptr++ = (xm1+(y+(z+t*nzb)*nyb)*nxb)*SITE_LEN;
	      }
	    if (y != 0)
	      {
	      *e_ux_ptr += 1;
	      *int_ptr++ = ODD_ADDR1(x,ym1,z,t);
	      *int_ptr++ = (x+(ym1+(z+t*nzb)*nyb)*nxb)*SITE_LEN
			      +MATRIX_SIZE;
	      }
	    if (z != 0)
	      {
	      *e_ux_ptr += 1;
	      *int_ptr++ = ODD_ADDR1(x,y,zm1,t);
	      *int_ptr++ = (x+(y+(zm1+t*nzb)*nyb)*nxb)*SITE_LEN
			      +2*MATRIX_SIZE;
	      }
	    if (t != 0)
	      {
	      *e_ux_ptr += 1;
	      *int_ptr++ = ODD_ADDR1(x,y,z,tm1);
	      *int_ptr++ = (x+(y+(z+tm1*nzb)*nyb)*nxb)*SITE_LEN
			      +3*MATRIX_SIZE;
	      }
	    e_ux_ptr += INTN_UX_TAB;
	    }
	  else
	    {
  	    int_ptr = o_ux_ptr + 1;
	    *o_ux_ptr = 0;
	    if (x != 0)
	      {
	      *o_ux_ptr += 1;
	      *int_ptr++ = EVEN_ADDR1(xm1,y,z,t);
	      *int_ptr++ = (xm1+(y+(z+t*nzb)*nyb)*nxb)*SITE_LEN;
	      }
	    if (y != 0)
	      {
	      *o_ux_ptr += 1;
	      *int_ptr++ = EVEN_ADDR1(x,ym1,z,t);
	      *int_ptr++ = (x+(ym1+(z+t*nzb)*nyb)*nxb)*SITE_LEN
			      +MATRIX_SIZE;
	      }
	    if (z != 0)
	      {
	      *o_ux_ptr += 1;
	      *int_ptr++ = EVEN_ADDR1(x,y,zm1,t);
	      *int_ptr++ = (x+(y+(zm1+t*nzb)*nyb)*nxb)*SITE_LEN
			      +2*MATRIX_SIZE;
	      }
	    if (t != 0)
	      {
	      *o_ux_ptr += 1;
	      *int_ptr++ = EVEN_ADDR1(x,y,z,tm1);
	      *int_ptr++ = (x+(y+(z+tm1*nzb)*nyb)*nxb)*SITE_LEN
			      +3*MATRIX_SIZE;
	      }
	    o_ux_ptr += INTN_UX_TAB;
	    }
	  }
	}
      }
    }

/*
 *  even intn_ux
 */
  addr_ptr = e_tab.intn_ux;
  for (i=0; i<n_even_sites; i++)
    {
	addr_ptr++;
	(*addr_ptr++) += BANK3_BASE;
	(*addr_ptr++) += BANK2_BASE;
	(*addr_ptr++) += BANK3_BASE + BANK_SIZE;
	(*addr_ptr++) += BANK2_BASE + BANK_SIZE;
	(*addr_ptr++) += BANK3_BASE;
	(*addr_ptr++) += BANK2_BASE;
	(*addr_ptr++) += BANK3_BASE + BANK_SIZE;
	(*addr_ptr++) += BANK2_BASE + BANK_SIZE;
    }
 
 /*
  *  odd intn_ux
  */
  addr_ptr = o_tab.intn_ux;
  for (i=0; i<n_even_sites; i++)
    {
	addr_ptr++;
	(*addr_ptr++) += BANK3_BASE;
	(*addr_ptr++) += BANK2_BASE;
	(*addr_ptr++) += BANK3_BASE + BANK_SIZE;
	(*addr_ptr++) += BANK2_BASE + BANK_SIZE;
	(*addr_ptr++) += BANK3_BASE;
	(*addr_ptr++) += BANK2_BASE;
	(*addr_ptr++) += BANK3_BASE + BANK_SIZE;
	(*addr_ptr++) += BANK2_BASE + BANK_SIZE;
    }

/*
 *  Address table for Udag . X calculation where all X's are on node.
 */
  e_ux_ptr = e_tab.udagx;
  o_ux_ptr = o_tab.udagx;
  for ( z=0; z<nzb; z++)
    {
    if (z<nzb-1) zp1 = z+1;  else  zp1 = 0;
    for ( y=0; y<nyb; y++)
      {
      if (y<nyb-1) yp1 = y+1;  else  yp1 = 0;
      for ( x=0; x<nxb; x++)
	{
	if (x<nxb-1) xp1 = x+1;  else  xp1 = 0;
	for ( t=0; t<ntb; t++)
	  {
	  if (t<ntb-1) tp1 = t+1;  else  tp1 = 0;

	  if ((x+y+z+t+base_odd)%2 == 0) 
	    {
	    surf_flag = 0;
  	    int_ptr = e_ux_ptr + 1;
	    *int_ptr++ = (x+(y+(z+t*nzb)*nyb)*nxb)*SITE_LEN;
	    if (x == nxb-1)
	      {
	      surf_flag |= 0x1;
	      int_ptr++;
	      }
	    else
	      *int_ptr++ = ODD_ADDR1(xp1,y,z,t);
	    if (y == nyb-1)
	      {
	      surf_flag |= 0x2;
	      int_ptr++;
	      }
	    else
	      *int_ptr++ = ODD_ADDR1(x,yp1,z,t);
	    if (z == nzb-1)
	      {
	      surf_flag |= 0x4;
	      int_ptr++;
	      }
	    else
	      *int_ptr++ = ODD_ADDR1(x,y,zp1,t);
	    if (t == ntb-1)
	      {
	      surf_flag |= 0x8;
	      int_ptr++;
	      }
	    else
	      *int_ptr++ = ODD_ADDR1(x,y,z,tp1);
	    if (x == 0)	surf_flag |= 0x10;
	    if (y == 0)	surf_flag |= 0x20;
	    if (z == 0)	surf_flag |= 0x40;
	    if (t == 0)	surf_flag |= 0x80;

  	    *e_ux_ptr = surf_flag;
  	    e_ux_ptr += UDAGX_TAB;
	    }
	  else
	    {
	    surf_flag = 0;
  	    int_ptr = o_ux_ptr + 1;
	    *(int_ptr++) = (x+(y+(z+t*nzb)*nyb)*nxb)*SITE_LEN;
	    if (x == nxb-1)
	      {
	      surf_flag |= 0x1;
	      int_ptr++;
	      }
	    else
	      *int_ptr++ = EVEN_ADDR1(xp1,y,z,t);
	    if (y == nyb-1)
	      {
	      surf_flag |= 0x2;
	      int_ptr++;
	      }
	    else
	      *int_ptr++ = EVEN_ADDR1(x,yp1,z,t);
	    if (z == nzb-1)
	      {
	      surf_flag |= 0x4;
	      int_ptr++;
	      }
	    else
	      *int_ptr++ = EVEN_ADDR1(x,y,zp1,t);
	    if (t == ntb-1)
	      {
	      surf_flag |= 0x8;
	      int_ptr++;
	      }
	    else
	      *int_ptr++ = EVEN_ADDR1(x,y,z,tp1);
	    if (x == 0)	surf_flag |= 0x10;
	    if (y == 0)	surf_flag |= 0x20;
	    if (z == 0)	surf_flag |= 0x40;
	    if (t == 0)	surf_flag |= 0x80;

  	    *o_ux_ptr = surf_flag;
  	    o_ux_ptr += UDAGX_TAB;
	    }
	  }
	}
      }
    }
 
 /*
  *  Address table for Udag . X calculation where all X's are from
  *  x+ off-node. The 4 comm_init() calls are used to find out how
  *  the X's are sent in block strided moves by the x+ neighboring
  *  boards.
  */
  e_ux_ptr = e_tab.udagx;
  o_ux_ptr = o_tab.udagx;
  nb[0] = ntb;
  nb[1] = nxb;
  nb[2] = nyb;
  nb[3] = nzb;

  comm_init(base_odd);
  for (i=0; i<4; i++)
    {
    rec_addr[i] = e_tab.m_trans.trans_des[0].start_addr[2*i].addr;
    rec_addr[i+4] = o_tab.m_trans.trans_des[0].start_addr[2*i].addr;
    }

  base_even_save = base_even;
  base_odd_save = base_odd;
  ntb_e2_save = ntb_e2;
  ntxb_e2_save = ntxb_e2;
  ntxyb_e2_save = ntxyb_e2;
  ntb_o2_save = ntb_o2;
  ntxb_o2_save = ntxb_o2;
  ntxyb_o2_save = ntxyb_o2;

  for (i=0; i<4; i++)
    {
    base_even = (base_even_save+nb[i])%2;
    base_odd = (base_odd_save+nb[i])%2;
    comm_init( base_odd );
    if (base_odd_save == base_odd)
      {
      ntb_e2 = ntb_e2_save;
      ntxb_e2 = ntxb_e2_save;
      ntxyb_e2 = ntxyb_e2_save;
      ntb_o2 = ntb_o2_save;
      ntxb_o2 = ntxb_o2_save;
      ntxyb_o2 = ntxyb_o2_save;
      }
    else
      {
      ntb_e2 = ntb_o2_save;
      ntxb_e2 = ntxb_o2_save;
      ntxyb_e2 = ntxyb_o2_save;
      ntb_o2 = ntb_e2_save;
      ntxb_o2 = ntxb_e2_save;
      ntxyb_o2 = ntxyb_e2_save;
      }

    for (j=0; j<e_tab.m_trans.n_transfer[i]; j++)
      {
      send_offset = e_tab.m_trans.trans_des[j].start_addr[2*i+1].addr;
      send_blk_len = e_tab.m_trans.trans_des[j].blk_des[2*i+1].blk_len;
      send_stride = e_tab.m_trans.trans_des[j].blk_des[2*i+1].stride;
      send_num_blk = e_tab.m_trans.trans_des[j].blk_des[2*i+1].num_blk;
      for (k=0;k<send_num_blk*(send_blk_len/VECT_LEN);k++)
	{
	for (xb[0]=0; xb[0]<nb[0]; xb[0]++) 
	  for (xb[1]=0; xb[1]<nb[1]; xb[1]++) 
	    for (xb[2]=0; xb[2]<nb[2]; xb[2]++) 
	      for (xb[3]=0; xb[3]<nb[3]; xb[3]++) 
	      {
	      if ((xb[0]+xb[1]+xb[2]+xb[3]+base_odd)%2 != 0)
		if (send_offset == ODD_ADDR1(xb[1],xb[2],xb[3],xb[0]))
		  {
		  if (xb[i] == 0)
		    {
		    xb[i] = nb[i]-1;
		    *(e_ux_ptr + (EVEN_ADDR1(xb[1],xb[2],xb[3],xb[0])
		      /VECT_LEN)*UDAGX_TAB+(i+3)%4+2) = rec_addr[i];
		    }
		  goto next_0;    /* break out of the for loops */
		  }
	      }
	next_0:
	rec_addr[i] += VECT_LEN*sizeof(IFloat);
	if ( (k+1) % (send_blk_len/VECT_LEN) != 0 )
	  send_offset += VECT_LEN;
	else
	  send_offset += VECT_LEN + send_stride - 1;
	}
      }

    for (j=0; j<o_tab.m_trans.n_transfer[i]; j++)
      {
      send_offset = o_tab.m_trans.trans_des[j].start_addr[2*i+1].addr;
      send_blk_len = o_tab.m_trans.trans_des[j].blk_des[2*i+1].blk_len;
      send_stride = o_tab.m_trans.trans_des[j].blk_des[2*i+1].stride;
      send_num_blk = o_tab.m_trans.trans_des[j].blk_des[2*i+1].num_blk;
      for (k=0;k<send_num_blk*(send_blk_len/VECT_LEN);k++)
	{
	for (xb[0]=0; xb[0]<nb[0]; xb[0]++) 
	  for (xb[1]=0; xb[1]<nb[1]; xb[1]++) 
	    for (xb[2]=0; xb[2]<nb[2]; xb[2]++) 
	      for (xb[3]=0; xb[3]<nb[3]; xb[3]++) 
	      {
	      if ((xb[0]+xb[1]+xb[2]+xb[3]+base_odd)%2 == 0)
		if (send_offset == EVEN_ADDR1(xb[1],xb[2],xb[3],xb[0]))
		  {
		  if (xb[i] == 0)
		    {
		    xb[i] = nb[i]-1;
		    *(o_ux_ptr + (ODD_ADDR1(xb[1],xb[2],xb[3],xb[0])
		      /VECT_LEN)*UDAGX_TAB+(i+3)%4+2) = rec_addr[i+4];
		    }
		  goto next_1;    /* break out of the for loops */
		  }
	      }
	next_1:
	rec_addr[i+4] += VECT_LEN*sizeof(IFloat);
	if ( (k+1) % (send_blk_len/VECT_LEN) != 0 )
	  send_offset += VECT_LEN;
	else
	  send_offset += VECT_LEN + send_stride - 1;
	}
      }
    }
  comm_init(base_odd_save);

/*
 *  even udagx
 */
  addr_ptr = e_tab.udagx;
  for (i=0; i<n_even_sites; i++)
    {
	addr_ptr++;
	(*addr_ptr++)+=BANK1_BASE;
	(*addr_ptr++)+=BANK3_BASE;
	(*addr_ptr++)+=BANK3_BASE+BANK_SIZE;
	(*addr_ptr++)+=BANK3_BASE;
	(*addr_ptr++)+=BANK3_BASE+BANK_SIZE;
    }

/*
 *  odd udagx
 */
  addr_ptr = o_tab.udagx;
  for (i=0; i<n_even_sites; i++)
    {
	addr_ptr++;
	(*addr_ptr++)+=BANK1_BASE;
	(*addr_ptr++)+=BANK3_BASE;
	(*addr_ptr++)+=BANK3_BASE+BANK_SIZE;
	(*addr_ptr++)+=BANK3_BASE;
	(*addr_ptr++)+=BANK3_BASE+BANK_SIZE;
    }


//  Added by RDM to set up SCU system calls for staggered cg
  StagSCUSetup();
}



/*************************************************************/
/***************************************************************

    storage order for X's : t x y z.  This is the most
    common case, particularly for zero temperature calculations
    where ntb is large.  When both nxb and nyb are even,
    the amount of data transfered in the t direction are
    doubled than what is necessary so that only one stride is
    needed.  For both nxb, nyb odd, one stride of serial
    communication is sufficient for all directions.

    Notation:
	Even and Odd are refered to the results of dirac operator.
    For example,  when one calculates X_even = D X_odd, then only
    even address table,  e_tab,  is needed in the calculation.

    All length are in terms of IFloats.  All offsets are to be
    added to pointers.  Except for rec.addr and p_send.addr for
    serial units which are absolute addresses.  The m_send.addr
    is an offset to be added to a IFloat pointer.
*/
/*************************************************************/

void comm_init(int base_odd_flag)
{
  int i,length;
  int nxb,nyb,nzb,ntb;
  int base_even,base_odd;
  int top_even,top_odd;

/*
 * This following part is general, independent of the size of the
 * lattice
 */
//   nxb = do_arg.nxb;	nyb = do_arg.nyb;
//   nzb = do_arg.nzb;	ntb = do_arg.ntb;
  nxb = GJP.XnodeSites();
  nyb = GJP.YnodeSites();
  nzb = GJP.ZnodeSites();
  ntb = GJP.TnodeSites();
  if ( base_odd_flag ) 
    {	base_even = 0;	base_odd = 1;	}
  else
    {	base_even = 1;  base_odd = 0;	}
  if ( (base_odd_flag+nxb+nyb+nzb+ntb)%2 ) 
    {	top_even = 0;	top_odd = 1;	}
  else
    {	top_even = 1;  top_odd = 0;	}

/*
  if ( ntb%2 != 0 )
    {
    // printf("comm_init: ntb = %d is not an even number\n",ntb);
    // exit();
    }
  if ( nxb%2 != nyb%2 )
    {
    // printf("comm_init: nxb = %d, nyb= %d are not both even or odd\n",
		nxb, nyb);
    // exit();
    }
*/

  for (i=0; i<4; i++)
    {
    e_tab.p_send_len[i] = ((nsurf_sites[i] + top_odd)/2) * VECT_LEN;
    e_tab.p_rec_len[i] = ((nsurf_sites[i] + top_even)/2) * VECT_LEN;
    e_tab.m_send_len[i] = ((nsurf_sites[i] + base_odd)/2) * VECT_LEN;
    e_tab.m_rec_len[i] = ((nsurf_sites[i] + base_even)/2) * VECT_LEN;
    o_tab.p_send_len[i] = ((nsurf_sites[i] + top_even)/2) * VECT_LEN;
    o_tab.p_rec_len[i] = ((nsurf_sites[i] + top_odd)/2) * VECT_LEN;
    o_tab.m_send_len[i] = ((nsurf_sites[i] + base_even)/2) * VECT_LEN;
    o_tab.m_rec_len[i] = ((nsurf_sites[i] + base_odd)/2) * VECT_LEN;
    }

/*
 *  From here  ntb even, nxb and nyb both even or odd are assumed
 */
  if (nxb%2 == 0)
    {
    e_tab.p_rec_len[0] = nsurf_sites[0]*VECT_LEN;
    o_tab.p_rec_len[0] = nsurf_sites[0]*VECT_LEN;
    }

/*
 *  Transfers in the negative direction	z direction
 */
  length = ((ntb*nxb*nyb)/2)*VECT_LEN;

  e_tab.m_trans.n_transfer[3] = 1;
  e_tab.m_trans.trans_des[0].start_addr[7].send = 1;
  e_tab.m_trans.trans_des[0].start_addr[7].addr = 0;
  e_tab.m_trans.trans_des[0].blk_des[7].num_blk = 1;
  e_tab.m_trans.trans_des[0].blk_des[7].blk_len = length;
  e_tab.m_trans.trans_des[0].blk_des[7].stride = 1;
  e_tab.m_trans.trans_des[0].start_addr[6].rec = 1;
  e_tab.m_trans.trans_des[0].start_addr[6].addr = (int) (p_rec_buf
						+ srbuf_offset[3]);
  e_tab.m_trans.trans_des[0].blk_des[6].num_blk = 1;
  e_tab.m_trans.trans_des[0].blk_des[6].blk_len = length;
  e_tab.m_trans.trans_des[0].blk_des[6].stride = 1;

  o_tab.m_trans.n_transfer[3] = 1;
  o_tab.m_trans.trans_des[0].start_addr[7].send = 1;
  o_tab.m_trans.trans_des[0].start_addr[7].addr = 0;
  o_tab.m_trans.trans_des[0].blk_des[7].num_blk = 1;
  o_tab.m_trans.trans_des[0].blk_des[7].blk_len = length;
  o_tab.m_trans.trans_des[0].blk_des[7].stride = 1;
  o_tab.m_trans.trans_des[0].start_addr[6].rec = 1;
  o_tab.m_trans.trans_des[0].start_addr[6].addr = (int) (p_rec_buf
						+ srbuf_offset[3]);
  o_tab.m_trans.trans_des[0].blk_des[6].num_blk = 1;
  o_tab.m_trans.trans_des[0].blk_des[6].blk_len = length;
  o_tab.m_trans.trans_des[0].blk_des[6].stride = 1;

/*
 *  y direction
 */
  length = (ntb*nxb)/2 * VECT_LEN;

  e_tab.m_trans.n_transfer[2] = 1;
  e_tab.m_trans.trans_des[0].start_addr[5].send = 1;
  e_tab.m_trans.trans_des[0].start_addr[5].addr = 0;
  e_tab.m_trans.trans_des[0].blk_des[5].num_blk = nzb;
  e_tab.m_trans.trans_des[0].blk_des[5].blk_len = length;
  e_tab.m_trans.trans_des[0].blk_des[5].stride = length*(nyb-1) + 1;
  e_tab.m_trans.trans_des[0].start_addr[4].rec = 1;
  e_tab.m_trans.trans_des[0].start_addr[4].addr = (int) (p_rec_buf
						+ srbuf_offset[2]);
  e_tab.m_trans.trans_des[0].blk_des[4].num_blk = 1;
  e_tab.m_trans.trans_des[0].blk_des[4].blk_len = length*nzb;
  e_tab.m_trans.trans_des[0].blk_des[4].stride = 1;

  o_tab.m_trans.n_transfer[2] = 1;
  o_tab.m_trans.trans_des[0].start_addr[5].send = 1;
  o_tab.m_trans.trans_des[0].start_addr[5].addr = 0;
  o_tab.m_trans.trans_des[0].blk_des[5].num_blk = nzb;
  o_tab.m_trans.trans_des[0].blk_des[5].blk_len = length;
  o_tab.m_trans.trans_des[0].blk_des[5].stride = length*(nyb-1) + 1;
  o_tab.m_trans.trans_des[0].start_addr[4].rec = 1;
  o_tab.m_trans.trans_des[0].start_addr[4].addr = (int) (p_rec_buf
						+ srbuf_offset[2]);
  o_tab.m_trans.trans_des[0].blk_des[4].num_blk = 1;
  o_tab.m_trans.trans_des[0].blk_des[4].blk_len = length*nzb;
  o_tab.m_trans.trans_des[0].blk_des[4].stride = 1;

/*
 *  x direction
 */
  length = ntb/2 * VECT_LEN;

  e_tab.m_trans.n_transfer[1] = 1;
  e_tab.m_trans.trans_des[0].start_addr[3].send = 1;
  e_tab.m_trans.trans_des[0].start_addr[3].addr = 0;
  e_tab.m_trans.trans_des[0].blk_des[3].num_blk = nyb*nzb; 
  e_tab.m_trans.trans_des[0].blk_des[3].blk_len = length;
  e_tab.m_trans.trans_des[0].blk_des[3].stride = length*(nxb-1) + 1;
  e_tab.m_trans.trans_des[0].start_addr[2].rec = 1;
  e_tab.m_trans.trans_des[0].start_addr[2].addr = (int) (p_rec_buf
						+ srbuf_offset[1]);
  e_tab.m_trans.trans_des[0].blk_des[2].num_blk = 1;
  e_tab.m_trans.trans_des[0].blk_des[2].blk_len = length*nyb*nzb;
  e_tab.m_trans.trans_des[0].blk_des[2].stride = 1;

  o_tab.m_trans.n_transfer[1] = 1;
  o_tab.m_trans.trans_des[0].start_addr[3].send = 1;
  o_tab.m_trans.trans_des[0].start_addr[3].addr = 0;
  o_tab.m_trans.trans_des[0].blk_des[3].num_blk = nyb*nzb; 
  o_tab.m_trans.trans_des[0].blk_des[3].blk_len = length;
  o_tab.m_trans.trans_des[0].blk_des[3].stride = length*(nxb-1) + 1;
  o_tab.m_trans.trans_des[0].start_addr[2].rec = 1;
  o_tab.m_trans.trans_des[0].start_addr[2].addr = (int) (p_rec_buf
						+ srbuf_offset[1]);
  o_tab.m_trans.trans_des[0].blk_des[2].num_blk = 1;
  o_tab.m_trans.trans_des[0].blk_des[2].blk_len = length*nyb*nzb;
  o_tab.m_trans.trans_des[0].blk_des[2].stride = 1;

/*
 *  Remember:  e_tab contains odd sites to be transfered because
 *  the calculation is organized according to X_e = D X_o.  All
 *  information for compute X_e is contained in e_tab.
 */

/*
 *  t direction
 */
  length = VECT_LEN;

  e_tab.m_trans.n_transfer[0] = 1;
  e_tab.m_trans.trans_des[0].start_addr[1].send = 1;
  e_tab.m_trans.trans_des[0].start_addr[0].rec = 1;
  e_tab.m_trans.trans_des[0].start_addr[0].addr
				= (int) (p_rec_buf + srbuf_offset[0]);

  o_tab.m_trans.n_transfer[0] = 1;
  o_tab.m_trans.trans_des[0].start_addr[1].send = 1;
  o_tab.m_trans.trans_des[0].start_addr[0].rec = 1;
  o_tab.m_trans.trans_des[0].start_addr[0].addr
				= (int) (p_rec_buf + srbuf_offset[0]);

  if ( nxb%2 != 0 )
    {
    e_tab.m_trans.trans_des[0].start_addr[1].addr
				= (ntb/2)*base_even*VECT_LEN;
    e_tab.m_trans.trans_des[0].blk_des[1].num_blk
				= (nxb*nyb*nzb+base_odd)/2; 
    e_tab.m_trans.trans_des[0].blk_des[1].blk_len = length;
    e_tab.m_trans.trans_des[0].blk_des[1].stride = length*(ntb-1) + 1;
    e_tab.m_trans.trans_des[0].blk_des[0].num_blk = 1;
    e_tab.m_trans.trans_des[0].blk_des[0].blk_len 
				= length*(nxb*nyb*nzb+top_even)/2; 
    e_tab.m_trans.trans_des[0].blk_des[0].stride = 1;

    o_tab.m_trans.trans_des[0].start_addr[1].addr
				= (ntb/2)*base_odd*VECT_LEN;
    o_tab.m_trans.trans_des[0].blk_des[1].num_blk
				= (nxb*nyb*nzb+base_even)/2; 
    o_tab.m_trans.trans_des[0].blk_des[1].blk_len = length;
    o_tab.m_trans.trans_des[0].blk_des[1].stride = length*(ntb-1) + 1;
    o_tab.m_trans.trans_des[0].blk_des[0].num_blk = 1;
    o_tab.m_trans.trans_des[0].blk_des[0].blk_len
				= length*(nxb*nyb*nzb+top_odd)/2; 
    o_tab.m_trans.trans_des[0].blk_des[0].stride = 1;
    }
  else
    {
    e_tab.m_trans.trans_des[0].start_addr[1].addr = 0;
    e_tab.m_trans.trans_des[0].blk_des[1].num_blk = nxb*nyb*nzb; 
    e_tab.m_trans.trans_des[0].blk_des[1].blk_len = length;
    e_tab.m_trans.trans_des[0].blk_des[1].stride = length*((ntb-1)/2)+1;
    e_tab.m_trans.trans_des[0].blk_des[0].num_blk = 1;
    e_tab.m_trans.trans_des[0].blk_des[0].blk_len = length*nxb*nyb*nzb; 
    e_tab.m_trans.trans_des[0].blk_des[0].stride = 1;

    o_tab.m_trans.trans_des[0].start_addr[1].addr = 0;
    o_tab.m_trans.trans_des[0].blk_des[1].num_blk = nxb*nyb*nzb; 
    o_tab.m_trans.trans_des[0].blk_des[1].blk_len = length;
    o_tab.m_trans.trans_des[0].blk_des[1].stride = length*((ntb-1)/2)+1;
    o_tab.m_trans.trans_des[0].blk_des[0].num_blk = 1;
    o_tab.m_trans.trans_des[0].blk_des[0].blk_len = length*nxb*nyb*nzb; 
    o_tab.m_trans.trans_des[0].blk_des[0].stride = 1;
    }

/*
 *  Transfers to the positive direction
 */
  for (i=0; i<4; i++)
    {
    e_tab.p_trans.n_transfer[i] = 1;
    e_tab.p_trans.trans_des[0].start_addr[2*i].send = 1;
    e_tab.p_trans.trans_des[0].start_addr[2*i].addr
				= (int) (p_send_buf + srbuf_offset[i]);
    e_tab.p_trans.trans_des[0].blk_des[2*i].num_blk = 1; 
    e_tab.p_trans.trans_des[0].blk_des[2*i].blk_len
				= e_tab.p_send_len[i];
    e_tab.p_trans.trans_des[0].blk_des[2*i].stride = 1;
    e_tab.p_trans.trans_des[0].start_addr[2*i+1].rec = 1;
    e_tab.p_trans.trans_des[0].start_addr[2*i+1].addr
				= (int) (m_rec_buf + srbuf_offset[i]);
    e_tab.p_trans.trans_des[0].blk_des[2*i+1].num_blk = 1;
    e_tab.p_trans.trans_des[0].blk_des[2*i+1].blk_len
				= e_tab.m_rec_len[i];
    e_tab.p_trans.trans_des[0].blk_des[2*i+1].stride = 1;

    o_tab.p_trans.n_transfer[i] = 1;
    o_tab.p_trans.trans_des[0].start_addr[2*i].send = 1;
    o_tab.p_trans.trans_des[0].start_addr[2*i].addr
				= (int) (p_send_buf + srbuf_offset[i]);
    o_tab.p_trans.trans_des[0].blk_des[2*i].num_blk = 1; 
    o_tab.p_trans.trans_des[0].blk_des[2*i].blk_len
				= e_tab.p_send_len[i];
    o_tab.p_trans.trans_des[0].blk_des[2*i].stride = 1;
    o_tab.p_trans.trans_des[0].start_addr[2*i+1].rec = 1;
    o_tab.p_trans.trans_des[0].start_addr[2*i+1].addr
				= (int) (m_rec_buf + srbuf_offset[i]);
    o_tab.p_trans.trans_des[0].blk_des[2*i+1].num_blk = 1;
    o_tab.p_trans.trans_des[0].blk_des[2*i+1].blk_len
				= e_tab.m_rec_len[i];
    o_tab.p_trans.trans_des[0].blk_des[2*i+1].stride = 1;
    }

}


CPS_END_NAMESPACE
