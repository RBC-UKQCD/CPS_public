#include <config.h>
#include <stdio.h>
#include <math.h>
#include <util/dirac_op.h>
#include <util/omp_wrapper.h>
#include <qmp.h>

//#ifdef USE_BFM
#if 0
#include "/bgsys/drivers/ppcfloor/hwi/include/bqc/nd_rese_dcr.h"
#endif

CPS_START_NAMESPACE
/*! \file
  \brief  Routine used internally in the DiracOpWilson class.

  $Id: wilson_dslash_vec.C,v 1.5 2013-05-14 16:56:34 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2013-05-14 16:56:34 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qmp/wilson_dslash_vec.C,v 1.5 2013-05-14 16:56:34 chulwoo Exp $
//  $Id: wilson_dslash_vec.C,v 1.5 2013-05-14 16:56:34 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.5 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qmp/wilson_dslash_vec.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
/***************************************************************************/
/*                                                                         */
/* wilson_dslash: It calculates chi = Dslash * psi, or                     */
/*                                chi = DslashDag * psi, where Dslassh is  */
/* the Wilson fermion Dslash matrix. Dslash is a function of the gauge     */
/* fields u.                                                               */
/* cb = 0/1 denotes input on even/odd checkerboards i.e.                   */
/* cb = 1 --> Dslash_EO, acts on an odd column vector and produces an      */
/* even column vector                                                      */
/* cb = 0 --> Dslash_OE, acts on an even column vector and produces an     */
/* odd column vector,                                                      */
/*                                                                         */
/* dag = 0/1  results into calculating Dslash or DslashDagger.             */
/* lx,ly,lz,lt is the lattice size.                                        */
/*                                                                         */
/* This routine is to be used with scalar machines.                        */
/*                                                                         */
/* WARNING:                                                                 */
/*                                                                          */
/* This set of routines will work only if the node sublattices have         */
/* even number of sites in each direction.                                  */
/*                                                                          */
/***************************************************************************/

CPS_END_NAMESPACE
#include <util/data_types.h>
#include <util/vector.h>
#include <util/gjp.h>
#include <util/wilson.h>
#include <util/error.h>
#include <util/verbose.h>
#include <util/smalloc.h>
#include <util/time_cps.h>
#include <comms/sysfunc_cps.h>
CPS_START_NAMESPACE

#include "wilson_op.h"

//! Access to the elements of the \e SU(3) matrix
/*!
  Gets the element of the \e SU(3) matrix \e u with row \e row,
  column \e col and complex component \e d
*/
#define U(r,row,col,d)  *(u+(r+2*(row+3*(col+3*d))))
//! Access to the elements of a spinor vector.
/*!
  Gets the element of the spinor \e psi with spin \e s,
  colour \e c and complex component \e r
*/
#define PSI(r,c,s)      *(psi +(r+2*(c+3*s)))

//! As above, but the vector is called chi
#define CHI(r,c,s)      *(chi +(r+2*(c+3*s)))
#define TMP(r,c,s)      *(tmp +(r+2*(c+3*s))) 
#define TMP1(r,c,s)     *(tmp1+(r+2*(c+3*s))) 
#define TMP2(r,c,s)     *(tmp2+(r+2*(c+3*s))) 
#define TMP3(r,c,s)     *(tmp3+(r+2*(c+3*s))) 
#define TMP4(r,c,s)     *(tmp4+(r+2*(c+3*s))) 
#define TMP5(r,c,s)     *(tmp5+(r+2*(c+3*s))) 
#define TMP6(r,c,s)     *(tmp6+(r+2*(c+3*s))) 
#define TMP7(r,c,s)     *(tmp7+(r+2*(c+3*s))) 
#define TMP8(r,c,s)     *(tmp8+(r+2*(c+3*s))) 
#define FBUF(r,c,s)     *(fbuf+(r+2*(c+3*s))) 

static int  Printf(char *format,...){}
//#define Printf if ( QMP_is_primary_node() ) printf
//#define Printf printf


static unsigned long called=0;
//static int initted=0;
static int init_len=0;
static double setup=0;
static double local=0;
static double nonlocal=0;
static double qmp=0;

static inline void MOVE_VEC(IFloat *buf, IFloat *psi, int vec_len, unsigned long vec_offset){
        for(int i=0;i<vec_len;i++){
                moveMem(buf, psi, SPINOR_SIZE*sizeof(IFloat) / sizeof(char));
                buf += SPINOR_SIZE;
                psi += vec_offset;
        }
}

static inline void MOVE_VEC2(IFloat *buf, IFloat *psi, int vec_len, unsigned long vec_offset){
        for(int i=0;i<vec_len;i++){
                moveMem(buf, psi, SPINOR_SIZE*sizeof(IFloat) / sizeof(char));
                buf += vec_offset;
                psi += vec_offset;
        }
}



void wilson_dslash_vec(IFloat *chi_p_f, 
			IFloat *u_p_f, 
			IFloat *psi_p_f, 
			int cb,
			int dag,
			Wilson *wilson_p,
			int vec_len,
			unsigned long vec_offset)
{
	char *cname = "";
	char *fname = "wilson_dslash";
	int lx, ly, lz, lt;
	int mu;
//	int r, c, s;
	int vol;
	Float fbuf[SPINOR_SIZE];

	Float time = -dclock();
	for(int i=0;i<SPINOR_SIZE;i++) fbuf[i]=0.;
	
/*--------------------------------------------------------------------------*/
/* Initializations                                                          */
/*--------------------------------------------------------------------------*/
	int sdag;                 /* = +/-1 if dag = 0/1 */
	int cbn;
	Float *chi_p = (Float *) chi_p_f;
	Float *u_p = (Float *) u_p_f;
	Float *psi_p = (Float *) psi_p_f;


	lx = wilson_p->ptr[0];
	ly = wilson_p->ptr[1];
	lz = wilson_p->ptr[2];
	lt = wilson_p->ptr[3];
	vol = wilson_p->vol[0];
	static int local_comm[4];

	if(dag == 0)
		sdag = 1;
	else if(dag == 1)
		sdag = -1;
	else{
		ERR.General(" ",fname,"dag must be 0 or 1");
	}

	if(cb == 0)
		cbn = 1;
	else if(cb == 1)
		cbn = 0;
	else{
		ERR.General(" ",fname,"cb must be 0 or 1");
	}

	Float *ind_buf[8];
	unsigned long ind_nl[8];

	static unsigned long num_nl[8];
	static unsigned long *u_ind[8]; 
	static unsigned long *f_ind[8]; 
	static unsigned long *t_ind[8]; 
	static Float *Send_buf[8];
	static Float *Recv_buf[8];
	static QMP_msgmem_t Send_mem[8];
	static QMP_msgmem_t Recv_mem[8];
	static QMP_msghandle_t Send[8];
	static QMP_msghandle_t Recv[8];


	if(vec_len!=init_len){
		VRB.Result(cname,fname,"init_len(%d)!=vec_len(%d)\n",init_len,vec_len);
		if (init_len>0)
		for (int i =0;i<8;i++){
			if(u_ind[i])sfree(u_ind[i]);
			if(f_ind[i])sfree(f_ind[i]);
			if(Send_buf[i])sfree(Send_buf[i]);
			if(Recv_buf[i])sfree(Recv_buf[i]);
			if(Send_mem[i])QMP_free_msgmem(Send_mem[i]);
			if(Recv_mem[i])QMP_free_msgmem(Recv_mem[i]);
			if(Send[i])QMP_free_msghandle(Send[i]);
			if(Recv[i])QMP_free_msghandle(Recv[i]);
		}
		num_nl[0]=num_nl[4]= vol/lx;
		num_nl[1]=num_nl[5]= vol/ly;
		num_nl[2]=num_nl[6]= vol/lz;
		num_nl[3]=num_nl[7]= vol/lt;
		for (int i=0;i<4;i++){
			if ( GJP.Nodes(i) > 1 )  local_comm[i] = 0;
			else local_comm[i] = 1;
		}
		for (int i =0;i<8;i++){
			mu = i%4;
			int sign = 1 - 2*(i/4);
			if (local_comm[mu]){
				num_nl[i]=0;
				u_ind[i]=NULL;
				f_ind[i]=NULL;
				Send_buf[i]=NULL;
				Recv_buf[i]=NULL;
				Send_mem[i]=NULL;
				Recv_mem[i]=NULL;
				Send[i]=NULL;
				Recv[i]=NULL;
			} else{
				u_ind[i]=(unsigned long *)smalloc(cname,fname,"u_ind[i]",num_nl[i]*sizeof(unsigned long) );
				f_ind[i]=(unsigned long *)smalloc(cname,fname,"f_ind[i]",num_nl[i]*sizeof(unsigned long) );
				t_ind[i]=(unsigned long *)smalloc(cname,fname,"f_ind[i]",num_nl[i]*sizeof(unsigned long) );
				Send_buf[i]=(Float *)smalloc(cname,fname,"Send_buf[i]",num_nl[i]*SPINOR_SIZE*sizeof(Float)*vec_len );
				Recv_buf[i]=(Float *)smalloc(cname,fname,"Recv_buf[i]",num_nl[i]*SPINOR_SIZE*sizeof(Float)*vec_len );
				Send_mem[i] = QMP_declare_msgmem(Send_buf[i], num_nl[i]*SPINOR_SIZE*sizeof(IFloat)*vec_len);
				Recv_mem[i] = QMP_declare_msgmem(Recv_buf[i], num_nl[i]*SPINOR_SIZE*sizeof(IFloat)*vec_len);
				Send[i] = QMP_declare_send_relative(Send_mem[i], mu,-sign, 0);
				Recv[i] = QMP_declare_receive_relative(Recv_mem[i], mu, sign, 0);
			}
		}
		Printf("initted\n");
		init_len=vec_len;
	}
	for (int i =0;i<8;i++){
		ind_nl[i]=0;
		ind_buf[i]=Send_buf[i];
	}

	if (called%100==0)
	VRB.Result(cname,fname,"local_comm=%d %d %d %d\n",local_comm[0],local_comm[1],local_comm[2],local_comm[3]);

//
//  non-local send
//
//
//omp_set_num_threads(8);
//	omp_set_num_threads(64);
#pragma omp parallel for default(shared) private(mu)
	for (int dir=0;dir<8;dir++){
		int x, y, z, t;
		int r, c, s;
		int xp, yp, zp, tp;
		int xm, ym, zm, tm;
		int xyzt;
		int xpyzt, xypzt, xyzpt, xyztp;
		int xmyzt, xymzt, xyzmt, xyztm;
		int parity;
		Float *chi;
		Float *u;
		Float *psi;
		Float tmp[SPINOR_SIZE];
		Float tmp1[SPINOR_SIZE];
		Float tmp2[SPINOR_SIZE];
		Float tmp3[SPINOR_SIZE];
		Float tmp4[SPINOR_SIZE];
		Float tmp5[SPINOR_SIZE];
		Float tmp6[SPINOR_SIZE];
		Float tmp7[SPINOR_SIZE];
		Float tmp8[SPINOR_SIZE];
		for(x=0; x<lx; x++)
		for(y=0; y<ly; y++)
		for(z=0; z<lz; z++)
		for(t=0; t<lt; t++){
		parity = x+y+z+t;
		parity = parity % 2;
		if(parity == cbn){
			Printf("non-local send: %d %d %d %d\n",x,y,z,t);
	
			/* x,y,z,t addressing of cbn checkerboard */
			xyzt = (x/2)+(lx/2)*(y+ly*(z+lz*t));
	
			/* x,y,z,t addressing of cb checkerboard */
			xp = (x+1) % lx;
			yp = (y+1) % ly;
			zp = (z+1) % lz;
			tp = (t+1) % lt;
			xm = x-1 + ( (lx-x)/lx ) * lx;
			ym = y-1 + ( (ly-y)/ly ) * ly;
			zm = z-1 + ( (lz-z)/lz ) * lz;
			tm = t-1 + ( (lt-t)/lt ) * lt;
			xpyzt = (xp/2)+(lx/2)*(y+ly*(z+lz*t));
			xmyzt = (xm/2)+(lx/2)*(y+ly*(z+lz*t));
			xypzt = (x/2)+(lx/2)*(yp+ly*(z+lz*t));
			xymzt = (x/2)+(lx/2)*(ym+ly*(z+lz*t));
			xyzpt = (x/2)+(lx/2)*(y+ly*(zp+lz*t));
			xyzmt = (x/2)+(lx/2)*(y+ly*(zm+lz*t));
			xyztp = (x/2)+(lx/2)*(y+ly*(z+lz*tp));
			xyztm = (x/2)+(lx/2)*(y+ly*(z+lz*tm));
	
			
			chi = chi_p + SPINOR_SIZE * xyzt;
	
			/* 1-gamma_0 */
			/*-----------*/
			mu = 0;
			if((dir==0)&&(x == lx-1) && (!local_comm[mu])){
				u   = u_p + GAUGE_SIZE * ( xyzt  + vol * cbn);
				psi = psi_p + SPINOR_SIZE * xpyzt;
				*(u_ind[0]+ind_nl[0])=xyzt+vol*cbn;
				*(f_ind[0]+ind_nl[0])=xpyzt;
				*(t_ind[0]+ind_nl[0])=xyzt;
				MOVE_VEC(ind_buf[0],psi, vec_len,vec_offset);
				ind_buf[0] += SPINOR_SIZE*vec_len;
				ind_nl[0] ++;
			} 
	
	
			/* 1-gamma_1 */
			/*-----------*/
			u   = u_p + GAUGE_SIZE * ( xyzt  + vol * cbn);
			psi = psi_p + SPINOR_SIZE * xypzt;
			mu = 1;
			if((dir==1)&&(y == ly-1) && (!local_comm[mu])){
				*(u_ind[1]+ind_nl[1])=xyzt+vol*cbn;
				*(f_ind[1]+ind_nl[1])=xypzt;
				*(t_ind[1]+ind_nl[1])=xyzt;
				MOVE_VEC(ind_buf[1],psi, vec_len,vec_offset);
				ind_buf[1] += SPINOR_SIZE*vec_len;
				ind_nl[1] ++;
			}
	
	
			/* 1-gamma_2 */
			/*-----------*/
			u   = u_p + GAUGE_SIZE * ( xyzt  + vol * cbn);
			psi = psi_p + SPINOR_SIZE * xyzpt;
			mu = 2;
			if((dir==2)&&(z == lz-1) && (!local_comm[mu])){
				*(u_ind[2]+ind_nl[2])=xyzt+vol*cbn;
				*(f_ind[2]+ind_nl[2])=xyzpt;
				*(t_ind[2]+ind_nl[2])=xyzt;
				MOVE_VEC(ind_buf[2],psi, vec_len,vec_offset);
				ind_buf[2] += SPINOR_SIZE*vec_len;
				ind_nl[2] ++;
			}
	
	
			/* 1-gamma_3 */
			/*-----------*/
			u   = u_p + GAUGE_SIZE * ( xyzt  + vol * cbn);
			psi = psi_p + SPINOR_SIZE * xyztp;
			mu = 3;
			if((dir==3)&&(t == lt-1) && (!local_comm[mu])){
				*(u_ind[3]+ind_nl[3])=xyzt+vol*cbn;
				*(f_ind[3]+ind_nl[3])=xyztp;
				*(t_ind[3]+ind_nl[3])=xyzt;
				MOVE_VEC(ind_buf[3],psi, vec_len,vec_offset);
				ind_buf[3] += SPINOR_SIZE*vec_len;
				ind_nl[3] ++;
			}
	
	
			/* 1+gamma_0 */
			/*-----------*/
			u   = u_p + GAUGE_SIZE * ( xmyzt  + vol * cb);
			psi = psi_p + SPINOR_SIZE * xmyzt;
			mu = 0;
			if((dir==4)&&x == 0 && (!local_comm[mu])){
				Printf("getMinusData((IFloat *)fbuf, (IFloat *)tmp5, SPINOR_SIZE, 0);\n");
				*(u_ind[4]+ind_nl[4])=xmyzt+vol*cb;
//				VRB.Result(cname,fname,"u_ind[4][ind_nl[4](%d)](%d)=u_cb(%d)\n",ind_nl[4],*(u_ind[4]+ind_nl[4]),xmyzt+vol*cb);
				*(f_ind[4]+ind_nl[4])=xmyzt;
				*(t_ind[4]+ind_nl[4])=xyzt;
				for(int vec_ind=0;vec_ind<vec_len;vec_ind++){
				for(c=0;c<3;c++){
					TMP(0,c,0) = PSI(0,c,0) + sdag * ( -PSI(1,c,3) ); 
					TMP(1,c,0) = PSI(1,c,0) + sdag * (  PSI(0,c,3) ); 
		
					TMP(0,c,1) = PSI(0,c,1) + sdag * ( -PSI(1,c,2) ); 
					TMP(1,c,1) = PSI(1,c,1) + sdag * (  PSI(0,c,2) ); 
		
					TMP(0,c,2) = PSI(0,c,2) + sdag * (  PSI(1,c,1) ); 
					TMP(1,c,2) = PSI(1,c,2) + sdag * ( -PSI(0,c,1) ); 
		
					TMP(0,c,3) = PSI(0,c,3) + sdag * (  PSI(1,c,0) ); 
					TMP(1,c,3) = PSI(1,c,3) + sdag * ( -PSI(0,c,0) ); 
				}
				/* multiply by U_mu */
				mu = 0;
				for(s=0;s<4;s++){
					for(c=0;c<3;c++){
						TMP5(0,c,s) = (  U(0,0,c,mu) * TMP(0,0,s)
								+ U(0,1,c,mu) * TMP(0,1,s)
								+ U(0,2,c,mu) * TMP(0,2,s) 
								+ U(1,0,c,mu) * TMP(1,0,s)
								+ U(1,1,c,mu) * TMP(1,1,s)
								+ U(1,2,c,mu) * TMP(1,2,s) );
						TMP5(1,c,s) = (  U(0,0,c,mu) * TMP(1,0,s)
								+ U(0,1,c,mu) * TMP(1,1,s)
								+ U(0,2,c,mu) * TMP(1,2,s) 
								- U(1,0,c,mu) * TMP(0,0,s)
								- U(1,1,c,mu) * TMP(0,1,s)
								- U(1,2,c,mu) * TMP(0,2,s) );
					}
				}
				moveMem((IFloat *)ind_buf[4], (IFloat *)tmp5, 
				SPINOR_SIZE*sizeof(Float) / sizeof(char));
				ind_buf[4] += SPINOR_SIZE;
				psi = psi + vec_offset;
				}
				ind_nl[4] ++;
	
			}
	
	
	
			/* 1+gamma_1 */
			/*-----------*/
			u   = u_p + GAUGE_SIZE * ( xymzt  + vol * cb);
			psi = psi_p + SPINOR_SIZE * xymzt;
			mu = 1;
			if((dir==5)&&y == 0 && (!local_comm[mu])){
Printf("getMinusData((IFloat *)fbuf, (IFloat *)tmp6, SPINOR_SIZE, 1);\n");
				*(u_ind[5]+ind_nl[5])=xymzt+vol*cb;
//				VRB.Result(cname,fname,"u_ind[5][ind_nl[5](%d)](%d)=u_cb(%d)\n",ind_nl[5],*(u_ind[5]+ind_nl[5]),xymzt  + vol * cb);
				*(f_ind[5]+ind_nl[5])=xymzt;
				*(t_ind[5]+ind_nl[5])=xyzt;
				for(int vec_ind=0;vec_ind<vec_len;vec_ind++){
				for(c=0;c<3;c++){
					TMP(0,c,0) = PSI(0,c,0) + sdag * ( -PSI(0,c,3) ); 
					TMP(1,c,0) = PSI(1,c,0) + sdag * ( -PSI(1,c,3) ); 
		
					TMP(0,c,1) = PSI(0,c,1) + sdag * (  PSI(0,c,2) ); 
					TMP(1,c,1) = PSI(1,c,1) + sdag * (  PSI(1,c,2) ); 
		
					TMP(0,c,2) = PSI(0,c,2) + sdag * (  PSI(0,c,1) ); 
					TMP(1,c,2) = PSI(1,c,2) + sdag * (  PSI(1,c,1) ); 
		
					TMP(0,c,3) = PSI(0,c,3) + sdag * ( -PSI(0,c,0) ); 
					TMP(1,c,3) = PSI(1,c,3) + sdag * ( -PSI(1,c,0) ); 
				}
				/* multiply by U_mu */
				mu = 1;
				for(s=0;s<4;s++){
					for(c=0;c<3;c++){
						TMP6(0,c,s) = (  U(0,0,c,mu) * TMP(0,0,s)
								+ U(0,1,c,mu) * TMP(0,1,s)
								+ U(0,2,c,mu) * TMP(0,2,s) 
								+ U(1,0,c,mu) * TMP(1,0,s)
								+ U(1,1,c,mu) * TMP(1,1,s)
								+ U(1,2,c,mu) * TMP(1,2,s) );
						TMP6(1,c,s) = (  U(0,0,c,mu) * TMP(1,0,s)
								+ U(0,1,c,mu) * TMP(1,1,s)
								+ U(0,2,c,mu) * TMP(1,2,s) 
								- U(1,0,c,mu) * TMP(0,0,s)
								- U(1,1,c,mu) * TMP(0,1,s)
								- U(1,2,c,mu) * TMP(0,2,s) );
					}
				}
				moveMem((IFloat *)ind_buf[5], (IFloat *)tmp6, 
				SPINOR_SIZE*sizeof(Float) / sizeof(char));
Printf("getMinusData((IFloat *)fbuf, (IFloat *)tmp6, SPINOR_SIZE, 1);\n");
				ind_buf[5] += SPINOR_SIZE;
				psi = psi + vec_offset;
				}
				ind_nl[5] ++;
			}
	
	
			/* 1+gamma_2 */
			/*-----------*/
			u   = u_p + GAUGE_SIZE * ( xyzmt  + vol * cb);
			psi = psi_p + SPINOR_SIZE * xyzmt;
			mu = 2;
			if((dir==6)&&z == 0 && (!local_comm[mu])){
				*(u_ind[6]+ind_nl[6])=xyzmt+vol*cb;
				*(f_ind[6]+ind_nl[6])=xyzmt;
				*(t_ind[6]+ind_nl[6])=xyzt;
				for(int vec_ind=0;vec_ind<vec_len;vec_ind++){
				for(c=0;c<3;c++){
					TMP(0,c,0) = PSI(0,c,0) + sdag * ( -PSI(1,c,2) ); 
					TMP(1,c,0) = PSI(1,c,0) + sdag * (  PSI(0,c,2) ); 
		
					TMP(0,c,1) = PSI(0,c,1) + sdag * (  PSI(1,c,3) ); 
					TMP(1,c,1) = PSI(1,c,1) + sdag * ( -PSI(0,c,3) ); 
		
					TMP(0,c,2) = PSI(0,c,2) + sdag * (  PSI(1,c,0) ); 
					TMP(1,c,2) = PSI(1,c,2) + sdag * ( -PSI(0,c,0) ); 
		
					TMP(0,c,3) = PSI(0,c,3) + sdag * ( -PSI(1,c,1) ); 
					TMP(1,c,3) = PSI(1,c,3) + sdag * (  PSI(0,c,1) ); 
				}
				/* multiply by U_mu */
				for(s=0;s<4;s++){
					for(c=0;c<3;c++){
						TMP7(0,c,s) = (  U(0,0,c,mu) * TMP(0,0,s)
								+ U(0,1,c,mu) * TMP(0,1,s)
								+ U(0,2,c,mu) * TMP(0,2,s) 
								+ U(1,0,c,mu) * TMP(1,0,s)
								+ U(1,1,c,mu) * TMP(1,1,s)
								+ U(1,2,c,mu) * TMP(1,2,s) );
						TMP7(1,c,s) = (  U(0,0,c,mu) * TMP(1,0,s)
								+ U(0,1,c,mu) * TMP(1,1,s)
								+ U(0,2,c,mu) * TMP(1,2,s) 
								- U(1,0,c,mu) * TMP(0,0,s)
								- U(1,1,c,mu) * TMP(0,1,s)
								- U(1,2,c,mu) * TMP(0,2,s) );
					}
				}
				Printf("getMinusData((IFloat *)fbuf, (IFloat *)tmp7, SPINOR_SIZE, 2);\n");
				moveMem((IFloat *)ind_buf[6], (IFloat *)tmp7, 
				SPINOR_SIZE*sizeof(Float) / sizeof(char));
				ind_buf[6] += SPINOR_SIZE;
				psi = psi + vec_offset;
				}
				ind_nl[6] ++;
			}
	
	
			/* 1+gamma_3 */
			/*-----------*/
			u   = u_p + GAUGE_SIZE * ( xyztm  + vol * cb);
			psi = psi_p + SPINOR_SIZE * xyztm;
			mu = 3;
			if((dir==7)&&t == 0 && (!local_comm[mu])){
				*(u_ind[7]+ind_nl[7])=xyztm+vol*cb;
				*(f_ind[7]+ind_nl[7])=xyztm;
				*(t_ind[7]+ind_nl[7])=xyzt;
				for(int vec_ind=0;vec_ind<vec_len;vec_ind++){
				for(c=0;c<3;c++){
					TMP(0,c,0) = PSI(0,c,0) + sdag * (  PSI(0,c,2) ); 
					TMP(1,c,0) = PSI(1,c,0) + sdag * (  PSI(1,c,2) ); 
		
					TMP(0,c,1) = PSI(0,c,1) + sdag * (  PSI(0,c,3) ); 
					TMP(1,c,1) = PSI(1,c,1) + sdag * (  PSI(1,c,3) ); 
		
					TMP(0,c,2) = PSI(0,c,2) + sdag * (  PSI(0,c,0) ); 
					TMP(1,c,2) = PSI(1,c,2) + sdag * (  PSI(1,c,0) ); 
		
					TMP(0,c,3) = PSI(0,c,3) + sdag * (  PSI(0,c,1) ); 
					TMP(1,c,3) = PSI(1,c,3) + sdag * (  PSI(1,c,1) ); 
				}
				/* multiply by U_mu */
				for(s=0;s<4;s++){
					for(c=0;c<3;c++){
						TMP8(0,c,s) = (  U(0,0,c,mu) * TMP(0,0,s)
								+ U(0,1,c,mu) * TMP(0,1,s)
								+ U(0,2,c,mu) * TMP(0,2,s) 
								+ U(1,0,c,mu) * TMP(1,0,s)
								+ U(1,1,c,mu) * TMP(1,1,s)
								+ U(1,2,c,mu) * TMP(1,2,s) );
						TMP8(1,c,s) = (  U(0,0,c,mu) * TMP(1,0,s)
								+ U(0,1,c,mu) * TMP(1,1,s)
								+ U(0,2,c,mu) * TMP(1,2,s) 
								- U(1,0,c,mu) * TMP(0,0,s)
								- U(1,1,c,mu) * TMP(0,1,s)
								- U(1,2,c,mu) * TMP(0,2,s) );
					}
				}
				Printf("getMinusData((IFloat *)fbuf, (IFloat *)tmp8, SPINOR_SIZE, 3);\n");
				moveMem((IFloat *)ind_buf[7], (IFloat *)tmp8, 
				SPINOR_SIZE*sizeof(Float) / sizeof(char));
				ind_buf[7] += SPINOR_SIZE;
				psi = psi + vec_offset;
				}
				ind_nl[7] ++;
			}
		}
		}
	}

	for(int i=0;i<8;i++){
		if (ind_nl[i]!=num_nl[i])
			VRB.Result(cname,fname,"ind_nl[%d](%d)!=num_nl[%d](%d)\n",i,ind_nl[i],i,num_nl[i]);
		if (!local_comm[i%4]){
			Printf("QMP_start(Recv[%d])(%p)\n",i,Recv[i]);
			QMP_start(Recv[i]);
			Printf("QMP_start(Send[%d])(%p)\n",i,Send[i]);
			QMP_start(Send[i]);
		}
	}

	time += dclock();
	setup += time;

	time = -dclock();
//	omp_set_num_threads(64);
/*--------------------------------------------------------------------------*/
/* Loop over sites                                                          */
/*--------------------------------------------------------------------------*/
	for(int i=0;i<SPINOR_SIZE;i++) fbuf[i]=0.;
	int index=0;
#pragma omp parallel for default(shared) private(mu)
	for(index = 0; index<vol*2;index++){
//	Printf("wilson_dslash: %d %d %d %d\n",x,y,z,t);
//	if ((called%10000==0) &&(!UniqueID())){
//		printf("wilson_dslash: index=%d thread %d of %d\n",index,omp_get_thread_num(),omp_get_num_threads());
//	}
	int r, c, s;
	int x, y, z, t;
	int xp, yp, zp, tp;
	int xm, ym, zm, tm;
	int xyzt;
	int xpyzt, xypzt, xyzpt, xyztp;
	int xmyzt, xymzt, xyzmt, xyztm;
	int parity;
	Float *chi;
	Float *u;
	Float *psi;
	Float tmp[SPINOR_SIZE];
	Float tmp1[SPINOR_SIZE];
	Float tmp2[SPINOR_SIZE];
	Float tmp3[SPINOR_SIZE];
	Float tmp4[SPINOR_SIZE];
	Float tmp5[SPINOR_SIZE];
	Float tmp6[SPINOR_SIZE];
	Float tmp7[SPINOR_SIZE];
	Float tmp8[SPINOR_SIZE];

	int temp = index;
	x = temp % lx; temp = temp/lx;
	y = temp % ly; temp = temp/ly;
	z = temp % lz; temp = temp/lz;
	t = temp % lt; temp = temp/lt;

if (0) 
if ((called%1000000==0) &&(!UniqueID())){
printf("wilson_dslash: %d %d %d %d %d: thread %d of %d tmp=%p \n",index,x,y,z,t,omp_get_thread_num(),omp_get_num_threads(),tmp);
}


	parity = x+y+z+t;
	parity = parity % 2;
	if(parity == cbn){

		/* x,y,z,t addressing of cbn checkerboard */
		xyzt = (x/2)+(lx/2)*(y+ly*(z+lz*t));
//	VRB.Result(fname,"local", "%d %d %d %d\n",x,y,z,t);

		/* x,y,z,t addressing of cb checkerboard */
		xp = (x+1) % lx;
		yp = (y+1) % ly;
		zp = (z+1) % lz;
		tp = (t+1) % lt;
		xm = x-1 + ( (lx-x)/lx ) * lx;
		ym = y-1 + ( (ly-y)/ly ) * ly;
		zm = z-1 + ( (lz-z)/lz ) * lz;
		tm = t-1 + ( (lt-t)/lt ) * lt;
		xpyzt = (xp/2)+(lx/2)*(y+ly*(z+lz*t));
		xmyzt = (xm/2)+(lx/2)*(y+ly*(z+lz*t));
		xypzt = (x/2)+(lx/2)*(yp+ly*(z+lz*t));
		xymzt = (x/2)+(lx/2)*(ym+ly*(z+lz*t));
		xyzpt = (x/2)+(lx/2)*(y+ly*(zp+lz*t));
		xyzmt = (x/2)+(lx/2)*(y+ly*(zm+lz*t));
		xyztp = (x/2)+(lx/2)*(y+ly*(z+lz*tp));
		xyztm = (x/2)+(lx/2)*(y+ly*(z+lz*tm));

		chi = chi_p + SPINOR_SIZE * xyzt;
        for(int vec_ind=0;vec_ind<vec_len;vec_ind++){
        for(s=0;s<4;s++)
            for(c=0;c<3;c++)
                for(r=0;r<2;r++)
            CHI(r,c,s) = 0.;
            chi += vec_offset;
        }

		/* 1-gamma_0 */
		/*-----------*/
		chi = chi_p + SPINOR_SIZE * xyzt;
		u   = u_p + GAUGE_SIZE * ( xyzt  + vol * cbn);
		psi = psi_p + SPINOR_SIZE * xpyzt;
		mu = 0;
		if((x == lx-1) && (!local_comm[mu])){
//			getPlusData((IFloat *)fbuf, (IFloat *) psi, SPINOR_SIZE, 0);
			psi = fbuf;
		} else
		for(int vec_ind=0;vec_ind<vec_len;vec_ind++){
		for(c=0;c<3;c++){
			TMP(0,c,0) = PSI(0,c,0) - sdag * ( -PSI(1,c,3) ); 
			TMP(1,c,0) = PSI(1,c,0) - sdag * (  PSI(0,c,3) ); 

			TMP(0,c,1) = PSI(0,c,1) - sdag * ( -PSI(1,c,2) ); 
			TMP(1,c,1) = PSI(1,c,1) - sdag * (  PSI(0,c,2) ); 

			TMP(0,c,2) = PSI(0,c,2) - sdag * (  PSI(1,c,1) ); 
			TMP(1,c,2) = PSI(1,c,2) - sdag * ( -PSI(0,c,1) ); 

			TMP(0,c,3) = PSI(0,c,3) - sdag * (  PSI(1,c,0) ); 
			TMP(1,c,3) = PSI(1,c,3) - sdag * ( -PSI(0,c,0) ); 
		}
		/* multiply by U_mu */
		for(s=0;s<4;s++){
			for(c=0;c<3;c++){
				TMP1(0,c,s) = (  U(0,c,0,mu) * TMP(0,0,s)
						+ U(0,c,1,mu) * TMP(0,1,s)
						+ U(0,c,2,mu) * TMP(0,2,s) 
						- U(1,c,0,mu) * TMP(1,0,s)
						- U(1,c,1,mu) * TMP(1,1,s)
						- U(1,c,2,mu) * TMP(1,2,s) );
				TMP1(1,c,s) = (  U(0,c,0,mu) * TMP(1,0,s)
						+ U(0,c,1,mu) * TMP(1,1,s)
						+ U(0,c,2,mu) * TMP(1,2,s) 
						+ U(1,c,0,mu) * TMP(0,0,s)
						+ U(1,c,1,mu) * TMP(0,1,s)
						+ U(1,c,2,mu) * TMP(0,2,s) );
			}
		}
        	for(s=0;s<4;s++)
            for(c=0;c<3;c++)
            for(r=0;r<2;r++)
        	CHI(r,c,s) = TMP1(r,c,s);
		    psi += vec_offset;
            chi += vec_offset;
        }


		/* 1-gamma_1 */
		/*-----------*/
		chi = chi_p + SPINOR_SIZE * xyzt;
		u   = u_p + GAUGE_SIZE * ( xyzt  + vol * cbn);
		psi = psi_p + SPINOR_SIZE * xypzt;
		mu = 1;
		if((y == ly-1) && (!local_comm[mu])){
//			getPlusData((IFloat *)fbuf, (IFloat *) psi, SPINOR_SIZE, 1);
			psi = fbuf;
		} else 
		for(int vec_ind=0;vec_ind<vec_len;vec_ind++){
		for(c=0;c<3;c++){
			TMP(0,c,0) = PSI(0,c,0) - sdag * ( -PSI(0,c,3) ); 
			TMP(1,c,0) = PSI(1,c,0) - sdag * ( -PSI(1,c,3) ); 

			TMP(0,c,1) = PSI(0,c,1) - sdag * (  PSI(0,c,2) ); 
			TMP(1,c,1) = PSI(1,c,1) - sdag * (  PSI(1,c,2) ); 

			TMP(0,c,2) = PSI(0,c,2) - sdag * (  PSI(0,c,1) ); 
			TMP(1,c,2) = PSI(1,c,2) - sdag * (  PSI(1,c,1) ); 

			TMP(0,c,3) = PSI(0,c,3) - sdag * ( -PSI(0,c,0) ); 
			TMP(1,c,3) = PSI(1,c,3) - sdag * ( -PSI(1,c,0) ); 
		}
		/* multiply by U_mu */
		for(s=0;s<4;s++){
			for(c=0;c<3;c++){
				TMP2(0,c,s) = (  U(0,c,0,mu) * TMP(0,0,s)
						+ U(0,c,1,mu) * TMP(0,1,s)
						+ U(0,c,2,mu) * TMP(0,2,s) 
						- U(1,c,0,mu) * TMP(1,0,s)
						- U(1,c,1,mu) * TMP(1,1,s)
						- U(1,c,2,mu) * TMP(1,2,s) );
				TMP2(1,c,s) = (  U(0,c,0,mu) * TMP(1,0,s)
						+ U(0,c,1,mu) * TMP(1,1,s)
						+ U(0,c,2,mu) * TMP(1,2,s) 
						+ U(1,c,0,mu) * TMP(0,0,s)
						+ U(1,c,1,mu) * TMP(0,1,s)
						+ U(1,c,2,mu) * TMP(0,2,s) );
			}
		}
        	for(s=0;s<4;s++)
            for(c=0;c<3;c++)
            for(r=0;r<2;r++)
        	CHI(r,c,s) += TMP2(r,c,s);
		    psi += vec_offset;
            chi += vec_offset;
        }


		/* 1-gamma_2 */
		/*-----------*/
		chi = chi_p + SPINOR_SIZE * xyzt;
		u   = u_p + GAUGE_SIZE * ( xyzt  + vol * cbn);
		psi = psi_p + SPINOR_SIZE * xyzpt;
		mu = 2;
		if((z == lz-1) && (!local_comm[mu])){
			psi = fbuf;
		} else 
		for(int vec_ind=0;vec_ind<vec_len;vec_ind++){
		for(c=0;c<3;c++){
			TMP(0,c,0) = PSI(0,c,0) - sdag * ( -PSI(1,c,2) ); 
			TMP(1,c,0) = PSI(1,c,0) - sdag * (  PSI(0,c,2) ); 

			TMP(0,c,1) = PSI(0,c,1) - sdag * (  PSI(1,c,3) ); 
			TMP(1,c,1) = PSI(1,c,1) - sdag * ( -PSI(0,c,3) ); 

			TMP(0,c,2) = PSI(0,c,2) - sdag * (  PSI(1,c,0) ); 
			TMP(1,c,2) = PSI(1,c,2) - sdag * ( -PSI(0,c,0) ); 

			TMP(0,c,3) = PSI(0,c,3) - sdag * ( -PSI(1,c,1) ); 
			TMP(1,c,3) = PSI(1,c,3) - sdag * (  PSI(0,c,1) ); 
		}
		/* multiply by U_mu */
		for(s=0;s<4;s++){
			for(c=0;c<3;c++){
				TMP3(0,c,s) = (  U(0,c,0,mu) * TMP(0,0,s)
						+ U(0,c,1,mu) * TMP(0,1,s)
						+ U(0,c,2,mu) * TMP(0,2,s) 
						- U(1,c,0,mu) * TMP(1,0,s)
						- U(1,c,1,mu) * TMP(1,1,s)
						- U(1,c,2,mu) * TMP(1,2,s) );
				TMP3(1,c,s) = (  U(0,c,0,mu) * TMP(1,0,s)
						+ U(0,c,1,mu) * TMP(1,1,s)
						+ U(0,c,2,mu) * TMP(1,2,s) 
						+ U(1,c,0,mu) * TMP(0,0,s)
						+ U(1,c,1,mu) * TMP(0,1,s)
						+ U(1,c,2,mu) * TMP(0,2,s) );
			}
		}
        	for(s=0;s<4;s++)
            for(c=0;c<3;c++)
            for(r=0;r<2;r++)
        	CHI(r,c,s) += TMP3(r,c,s);
		    psi += vec_offset;
            chi += vec_offset;
        }


		/* 1-gamma_3 */
		/*-----------*/
		chi = chi_p + SPINOR_SIZE * xyzt;
		u   = u_p + GAUGE_SIZE * ( xyzt  + vol * cbn);
		psi = psi_p + SPINOR_SIZE * xyztp;
		mu = 3;
		if((t == lt-1) && (!local_comm[mu])){
//			getPlusData((IFloat *)fbuf, (IFloat *) psi, SPINOR_SIZE, 3);
			psi = fbuf;
		} else
		for(int vec_ind=0;vec_ind<vec_len;vec_ind++){
		for(c=0;c<3;c++){
			TMP(0,c,0) = PSI(0,c,0) - sdag * (  PSI(0,c,2) ); 
			TMP(1,c,0) = PSI(1,c,0) - sdag * (  PSI(1,c,2) ); 

			TMP(0,c,1) = PSI(0,c,1) - sdag * (  PSI(0,c,3) ); 
			TMP(1,c,1) = PSI(1,c,1) - sdag * (  PSI(1,c,3) ); 

			TMP(0,c,2) = PSI(0,c,2) - sdag * (  PSI(0,c,0) ); 
			TMP(1,c,2) = PSI(1,c,2) - sdag * (  PSI(1,c,0) ); 

			TMP(0,c,3) = PSI(0,c,3) - sdag * (  PSI(0,c,1) ); 
			TMP(1,c,3) = PSI(1,c,3) - sdag * (  PSI(1,c,1) ); 
		}
		/* multiply by U_mu */
		mu = 3;
		for(s=0;s<4;s++){
			for(c=0;c<3;c++){
				TMP4(0,c,s) = (  U(0,c,0,mu) * TMP(0,0,s)
						+ U(0,c,1,mu) * TMP(0,1,s)
						+ U(0,c,2,mu) * TMP(0,2,s) 
						- U(1,c,0,mu) * TMP(1,0,s)
						- U(1,c,1,mu) * TMP(1,1,s)
						- U(1,c,2,mu) * TMP(1,2,s) );
				TMP4(1,c,s) = (  U(0,c,0,mu) * TMP(1,0,s)
						+ U(0,c,1,mu) * TMP(1,1,s)
						+ U(0,c,2,mu) * TMP(1,2,s) 
						+ U(1,c,0,mu) * TMP(0,0,s)
						+ U(1,c,1,mu) * TMP(0,1,s)
						+ U(1,c,2,mu) * TMP(0,2,s) );
			}
		}
        	for(s=0;s<4;s++)
            for(c=0;c<3;c++)
            for(r=0;r<2;r++)
        	CHI(r,c,s) += TMP4(r,c,s);
		    psi += vec_offset;
            chi += vec_offset;
        }


		/* 1+gamma_0 */
		/*-----------*/
		chi = chi_p + SPINOR_SIZE * xyzt;
		u   = u_p + GAUGE_SIZE * ( xmyzt  + vol * cb);
		psi = psi_p + SPINOR_SIZE * xmyzt;
		/* multiply by U_mu */
		mu = 0;
		if(x == 0 && (!local_comm[mu])){
			moveMem((IFloat *)tmp5, (IFloat *)fbuf, 
			SPINOR_SIZE*sizeof(Float) / sizeof(char));
		} else
		for(int vec_ind=0;vec_ind<vec_len;vec_ind++){
		for(c=0;c<3;c++){
			TMP(0,c,0) = PSI(0,c,0) + sdag * ( -PSI(1,c,3) ); 
			TMP(1,c,0) = PSI(1,c,0) + sdag * (  PSI(0,c,3) ); 

			TMP(0,c,1) = PSI(0,c,1) + sdag * ( -PSI(1,c,2) ); 
			TMP(1,c,1) = PSI(1,c,1) + sdag * (  PSI(0,c,2) ); 

			TMP(0,c,2) = PSI(0,c,2) + sdag * (  PSI(1,c,1) ); 
			TMP(1,c,2) = PSI(1,c,2) + sdag * ( -PSI(0,c,1) ); 

			TMP(0,c,3) = PSI(0,c,3) + sdag * (  PSI(1,c,0) ); 
			TMP(1,c,3) = PSI(1,c,3) + sdag * ( -PSI(0,c,0) ); 
		}
		for(s=0;s<4;s++){
			for(c=0;c<3;c++){
				TMP5(0,c,s) = (  U(0,0,c,mu) * TMP(0,0,s)
						+ U(0,1,c,mu) * TMP(0,1,s)
						+ U(0,2,c,mu) * TMP(0,2,s) 
						+ U(1,0,c,mu) * TMP(1,0,s)
						+ U(1,1,c,mu) * TMP(1,1,s)
						+ U(1,2,c,mu) * TMP(1,2,s) );
				TMP5(1,c,s) = (  U(0,0,c,mu) * TMP(1,0,s)
						+ U(0,1,c,mu) * TMP(1,1,s)
						+ U(0,2,c,mu) * TMP(1,2,s) 
						- U(1,0,c,mu) * TMP(0,0,s)
						- U(1,1,c,mu) * TMP(0,1,s)
						- U(1,2,c,mu) * TMP(0,2,s) );
			}
		}
        	for(s=0;s<4;s++)
            for(c=0;c<3;c++)
            for(r=0;r<2;r++)
        	CHI(r,c,s) += TMP5(r,c,s);
		    psi += vec_offset;
            chi += vec_offset;
        }



		/* 1+gamma_1 */
		/*-----------*/
		chi = chi_p + SPINOR_SIZE * xyzt;
		u   = u_p + GAUGE_SIZE * ( xymzt  + vol * cb);
		psi = psi_p + SPINOR_SIZE * xymzt;
		/* multiply by U_mu */
		mu = 1;
		if(y == 0 && (!local_comm[mu])){
//			getMinusData((IFloat *)fbuf, (IFloat *)tmp6, SPINOR_SIZE, 1);
			moveMem((IFloat *)tmp6, (IFloat *)fbuf, 
		SPINOR_SIZE*sizeof(Float) / sizeof(char));
		} else
		for(int vec_ind=0;vec_ind<vec_len;vec_ind++){
		for(c=0;c<3;c++){
			TMP(0,c,0) = PSI(0,c,0) + sdag * ( -PSI(0,c,3) ); 
			TMP(1,c,0) = PSI(1,c,0) + sdag * ( -PSI(1,c,3) ); 

			TMP(0,c,1) = PSI(0,c,1) + sdag * (  PSI(0,c,2) ); 
			TMP(1,c,1) = PSI(1,c,1) + sdag * (  PSI(1,c,2) ); 

			TMP(0,c,2) = PSI(0,c,2) + sdag * (  PSI(0,c,1) ); 
			TMP(1,c,2) = PSI(1,c,2) + sdag * (  PSI(1,c,1) ); 

			TMP(0,c,3) = PSI(0,c,3) + sdag * ( -PSI(0,c,0) ); 
			TMP(1,c,3) = PSI(1,c,3) + sdag * ( -PSI(1,c,0) ); 
		}
		for(s=0;s<4;s++){
			for(c=0;c<3;c++){
				TMP6(0,c,s) = (  U(0,0,c,mu) * TMP(0,0,s)
						+ U(0,1,c,mu) * TMP(0,1,s)
						+ U(0,2,c,mu) * TMP(0,2,s) 
						+ U(1,0,c,mu) * TMP(1,0,s)
						+ U(1,1,c,mu) * TMP(1,1,s)
						+ U(1,2,c,mu) * TMP(1,2,s) );
				TMP6(1,c,s) = (  U(0,0,c,mu) * TMP(1,0,s)
						+ U(0,1,c,mu) * TMP(1,1,s)
						+ U(0,2,c,mu) * TMP(1,2,s) 
						- U(1,0,c,mu) * TMP(0,0,s)
						- U(1,1,c,mu) * TMP(0,1,s)
						- U(1,2,c,mu) * TMP(0,2,s) );
			}
		}
        	for(s=0;s<4;s++)
            for(c=0;c<3;c++)
            for(r=0;r<2;r++)
        	CHI(r,c,s) += TMP6(r,c,s);
		    psi += vec_offset;
            chi += vec_offset;
        }


		/* 1+gamma_2 */
		/*-----------*/
		chi = chi_p + SPINOR_SIZE * xyzt;
		u   = u_p + GAUGE_SIZE * ( xyzmt  + vol * cb);
		psi = psi_p + SPINOR_SIZE * xyzmt;
		/* multiply by U_mu */
		mu = 2;
		if(z == 0 && (!local_comm[mu])){
//			getMinusData((IFloat *)fbuf, (IFloat *)tmp7, SPINOR_SIZE, 2);
			moveMem((IFloat *)tmp7, (IFloat *)fbuf, 
		SPINOR_SIZE*sizeof(Float) / sizeof(char));
		} else
		for(int vec_ind=0;vec_ind<vec_len;vec_ind++){
		for(c=0;c<3;c++){
			TMP(0,c,0) = PSI(0,c,0) + sdag * ( -PSI(1,c,2) ); 
			TMP(1,c,0) = PSI(1,c,0) + sdag * (  PSI(0,c,2) ); 

			TMP(0,c,1) = PSI(0,c,1) + sdag * (  PSI(1,c,3) ); 
			TMP(1,c,1) = PSI(1,c,1) + sdag * ( -PSI(0,c,3) ); 

			TMP(0,c,2) = PSI(0,c,2) + sdag * (  PSI(1,c,0) ); 
			TMP(1,c,2) = PSI(1,c,2) + sdag * ( -PSI(0,c,0) ); 

			TMP(0,c,3) = PSI(0,c,3) + sdag * ( -PSI(1,c,1) ); 
			TMP(1,c,3) = PSI(1,c,3) + sdag * (  PSI(0,c,1) ); 
		}
		for(s=0;s<4;s++){
			for(c=0;c<3;c++){
				TMP7(0,c,s) = (  U(0,0,c,mu) * TMP(0,0,s)
						+ U(0,1,c,mu) * TMP(0,1,s)
						+ U(0,2,c,mu) * TMP(0,2,s) 
						+ U(1,0,c,mu) * TMP(1,0,s)
						+ U(1,1,c,mu) * TMP(1,1,s)
						+ U(1,2,c,mu) * TMP(1,2,s) );
				TMP7(1,c,s) = (  U(0,0,c,mu) * TMP(1,0,s)
						+ U(0,1,c,mu) * TMP(1,1,s)
						+ U(0,2,c,mu) * TMP(1,2,s) 
						- U(1,0,c,mu) * TMP(0,0,s)
						- U(1,1,c,mu) * TMP(0,1,s)
						- U(1,2,c,mu) * TMP(0,2,s) );
			}
		}
        	for(s=0;s<4;s++)
            for(c=0;c<3;c++)
            for(r=0;r<2;r++)
        	CHI(r,c,s) += TMP7(r,c,s);
		    psi += vec_offset;
            chi += vec_offset;
        }


		/* 1+gamma_3 */
		/*-----------*/
		chi = chi_p + SPINOR_SIZE * xyzt;
		u   = u_p + GAUGE_SIZE * ( xyztm  + vol * cb);
		psi = psi_p + SPINOR_SIZE * xyztm;
		/* multiply by U_mu */
		mu = 3;
		if(t == 0 && (!local_comm[mu])){
//			getMinusData((IFloat *)fbuf, (IFloat *)tmp8, SPINOR_SIZE, 3);
			moveMem((IFloat *)tmp8, (IFloat *)fbuf, 
		SPINOR_SIZE*sizeof(Float) / sizeof(char));
		} else
		for(int vec_ind=0;vec_ind<vec_len;vec_ind++){
		for(c=0;c<3;c++){
			TMP(0,c,0) = PSI(0,c,0) + sdag * (  PSI(0,c,2) ); 
			TMP(1,c,0) = PSI(1,c,0) + sdag * (  PSI(1,c,2) ); 

			TMP(0,c,1) = PSI(0,c,1) + sdag * (  PSI(0,c,3) ); 
			TMP(1,c,1) = PSI(1,c,1) + sdag * (  PSI(1,c,3) ); 

			TMP(0,c,2) = PSI(0,c,2) + sdag * (  PSI(0,c,0) ); 
			TMP(1,c,2) = PSI(1,c,2) + sdag * (  PSI(1,c,0) ); 

			TMP(0,c,3) = PSI(0,c,3) + sdag * (  PSI(0,c,1) ); 
			TMP(1,c,3) = PSI(1,c,3) + sdag * (  PSI(1,c,1) ); 
		}
		for(s=0;s<4;s++){
			for(c=0;c<3;c++){
				TMP8(0,c,s) = (  U(0,0,c,mu) * TMP(0,0,s)
						+ U(0,1,c,mu) * TMP(0,1,s)
						+ U(0,2,c,mu) * TMP(0,2,s) 
						+ U(1,0,c,mu) * TMP(1,0,s)
						+ U(1,1,c,mu) * TMP(1,1,s)
						+ U(1,2,c,mu) * TMP(1,2,s) );
				TMP8(1,c,s) = (  U(0,0,c,mu) * TMP(1,0,s)
						+ U(0,1,c,mu) * TMP(1,1,s)
						+ U(0,2,c,mu) * TMP(1,2,s) 
						- U(1,0,c,mu) * TMP(0,0,s)
						- U(1,1,c,mu) * TMP(0,1,s)
						- U(1,2,c,mu) * TMP(0,2,s) );
			}
		}
        	for(s=0;s<4;s++)
            for(c=0;c<3;c++)
            for(r=0;r<2;r++)
        	CHI(r,c,s) += TMP8(r,c,s);
		    psi += vec_offset;
            chi += vec_offset;
        }

	}
	}
	time += dclock();
	local += time;

	time = -dclock();

	for(int i=0;i<8;i++){
		if (!local_comm[i%4]){
			QMP_status_t send_status = QMP_wait(Send[i]);
			if (send_status != QMP_SUCCESS)
				QMP_error("Send failed in wilson_dslash: %s\n", QMP_error_string(send_status));
			QMP_status_t rcv_status = QMP_wait(Recv[i]);
			if (rcv_status != QMP_SUCCESS)
				QMP_error("Receive failed in wilson_dslash: %s\n", QMP_error_string(rcv_status));
		}
		ind_nl[i]=0;
		ind_buf[i]=Recv_buf[i];
	}
	time += dclock();
	qmp += time;

	time = -dclock();

//
// non-local
//

#define USE_TEST2
#ifdef USE_TEST2
	index=0;
	int nl_offset[8];
	for(int i=0;i<8;i++){
		index += num_nl[i];
		nl_offset[i] = index;
	}

#undef NL_OMP
{
#pragma omp parallel for default(shared)
	for( int i=0;i<num_nl[0];i++){ 
//		int i_nl;
//		i_nl=i;
//		i_mu = 0;
	Float tmp[SPINOR_SIZE];
	Float tmp1[SPINOR_SIZE];

		/* 1-gamma_0 */
		/*-----------*/

		Float *chi = chi_p + SPINOR_SIZE * ( *(t_ind[0]+i) );
		Float *u   = u_p + GAUGE_SIZE * ( *(u_ind[0]+i) );
		Float *psi = ind_buf[0] + i* SPINOR_SIZE*vec_len;
		int r, c, s;
		for(int vec_ind=0;vec_ind<vec_len;vec_ind++){
			MINUSX(u,tmp,tmp1,sdag,psi);
			for(s=0;s<4;s++)
			for(c=0;c<3;c++)
			for(r=0;r<2;r++)
				CHI(r,c,s) += TMP1(r,c,s);
			psi += SPINOR_SIZE;
			chi +=vec_offset;
		}
		chi = chi_p + SPINOR_SIZE * ( *(t_ind[4]+i) );
		u   = u_p + GAUGE_SIZE * ( *(u_ind[4]+i) );
		psi = ind_buf[4] + i* SPINOR_SIZE*vec_len;
		for(int vec_ind=0;vec_ind<vec_len;vec_ind++){
			for(s=0;s<4;s++)
			for(c=0;c<3;c++)
			for(r=0;r<2;r++)
				CHI(r,c,s) += PSI(r,c,s);
			psi += SPINOR_SIZE;
			chi +=vec_offset;
		}
	}

#pragma omp parallel for default(shared)
	for( int i=0;i<num_nl[1];i++){ 
//		i_mu = 1;
	Float tmp[SPINOR_SIZE];
	Float tmp1[SPINOR_SIZE];
	Float tmp2[SPINOR_SIZE];

		Float *chi = chi_p + SPINOR_SIZE * ( *(t_ind[1]+i) );
		Float *u   = u_p + GAUGE_SIZE * ( *(u_ind[1]+i) );
		Float *psi = ind_buf[1] + i* SPINOR_SIZE*vec_len;
		int r, c, s;
		for(int vec_ind=0;vec_ind<vec_len;vec_ind++){
			MINUSY(u,tmp,tmp2,sdag,psi);
			for(s=0;s<4;s++)
			for(c=0;c<3;c++)
			for(r=0;r<2;r++)
				CHI(r,c,s) += TMP2(r,c,s);
			psi += SPINOR_SIZE;
			chi +=vec_offset;
		}
		chi = chi_p + SPINOR_SIZE * ( *(t_ind[5]+i) );
		u   = u_p + GAUGE_SIZE * ( *(u_ind[5]+i) );
		psi = ind_buf[5] + i* SPINOR_SIZE*vec_len;
		for(int vec_ind=0;vec_ind<vec_len;vec_ind++){
			for(s=0;s<4;s++)
			for(c=0;c<3;c++)
			for(r=0;r<2;r++)
				CHI(r,c,s) += PSI(r,c,s);
			psi += SPINOR_SIZE;
			chi +=vec_offset;
		}
	}

#pragma omp parallel for default(shared)
	for( int i=0;i<num_nl[2];i++){ 
//		i_mu = 0;
	Float tmp[SPINOR_SIZE];
	Float tmp3[SPINOR_SIZE];

		/* 1-gamma_0 */
		/*-----------*/

		Float *chi = chi_p + SPINOR_SIZE * ( *(t_ind[2]+i) );
		Float *u   = u_p + GAUGE_SIZE * ( *(u_ind[2]+i) );
		Float *psi = ind_buf[2] + i* SPINOR_SIZE*vec_len;
		int r, c, s;
		for(int vec_ind=0;vec_ind<vec_len;vec_ind++){
			MINUSZ(u,tmp,tmp3,sdag,psi);
			for(s=0;s<4;s++)
			for(c=0;c<3;c++)
			for(r=0;r<2;r++)
				CHI(r,c,s) += TMP3(r,c,s);
			psi += SPINOR_SIZE;
			chi +=vec_offset;
		}
		chi = chi_p + SPINOR_SIZE * ( *(t_ind[6]+i) );
		u   = u_p + GAUGE_SIZE * ( *(u_ind[6]+i) );
		psi = ind_buf[6] + i* SPINOR_SIZE*vec_len;
		for(int vec_ind=0;vec_ind<vec_len;vec_ind++){
			for(s=0;s<4;s++)
			for(c=0;c<3;c++)
			for(r=0;r<2;r++)
				CHI(r,c,s) += PSI(r,c,s);
			psi += SPINOR_SIZE;
			chi +=vec_offset;
		}
	}
#pragma omp parallel for default(shared)
	for( int i=0;i<num_nl[3];i++){ 
//		i_mu = 0;
	Float tmp[SPINOR_SIZE];
	Float tmp4[SPINOR_SIZE];

		/* 1-gamma_0 */
		/*-----------*/

		Float *chi = chi_p + SPINOR_SIZE * ( *(t_ind[3]+i) );
		Float *u   = u_p + GAUGE_SIZE * ( *(u_ind[3]+i) );
		Float *psi = ind_buf[3] + i* SPINOR_SIZE*vec_len;
		int r, c, s;
		for(int vec_ind=0;vec_ind<vec_len;vec_ind++){
			MINUST(u,tmp,tmp4,sdag,psi);
			for(s=0;s<4;s++)
			for(c=0;c<3;c++)
			for(r=0;r<2;r++)
				CHI(r,c,s) += TMP4(r,c,s);
			psi += SPINOR_SIZE;
			chi +=vec_offset;
		}
		chi = chi_p + SPINOR_SIZE * ( *(t_ind[7]+i) );
		u   = u_p + GAUGE_SIZE * ( *(u_ind[7]+i) );
		psi = ind_buf[7] + i* SPINOR_SIZE*vec_len;
		for(int vec_ind=0;vec_ind<vec_len;vec_ind++){
			for(s=0;s<4;s++)
			for(c=0;c<3;c++)
			for(r=0;r<2;r++)
				CHI(r,c,s) += PSI(r,c,s);
			psi += SPINOR_SIZE;
			chi +=vec_offset;
		}
	}
}
#endif
//exit(-44);


	
	
	time += dclock();
	nonlocal += time;

	called++;

	if (called%100==0){
		print_flops("wilson_dslash_vec()","local*100",0,local);
		print_flops("wilson_dslash_vec()","nonlocal*100",0,nonlocal);
		print_flops("wilson_dslash_vec()","qmp*100",0,qmp);
		print_flops("wilson_dslash_vec()","setup*100",0,setup);
		local=nonlocal=qmp=setup=0.;
//#ifdef USE_BFM
#if 0
{
	 char link_name[ND_RESE_DCR_num][10] = { "A-", "A+", "B-", "B+", "C-", "C+", "D-", "D+", "E-", "E+", "IO" };
    uint32_t i;
    for (i = 0; i < ND_RESE_DCR_num; i++) {
        uint64_t val_re = DCRReadUser(ND_RESE_DCR(i, RE_LINK_ERR_CNT));
        uint64_t val_retran = DCRReadUser(ND_RESE_DCR(i, SE_RETRANS_CNT));
        if (val_re || val_retran)
          printf("LINK Errors on RESE %s Recv Count = %ld Retran = %ld\n",
                 link_name[i],val_re,val_retran);
//        if (val_re > 3000){
        if (0){
// crash the node
		double *tmp = NULL;
		*tmp = 1.0;
		}
	}
}
#endif
	}
//	VRB.Result(cname,fname,"done");
DiracOp::CGflops += 1320*vol*vec_len;
}

#if 0
void wilson_dslash(IFloat *chi_p_f, 
			IFloat *u_p_f, 
			IFloat *psi_p_f, 
			int cb,
			int dag,
			Wilson *wilson_p)
{
	wilson_dslash_vec(chi_p_f,u_p_f,psi_p_f,cb,dag,wilson_p,1,0);
}
#endif


CPS_END_NAMESPACE

