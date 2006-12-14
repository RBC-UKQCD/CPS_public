#include <config.h>
CPS_START_NAMESPACE
/*--------------------------------------------------------------------*/
/*! \file
  \brief  Routine used internally in the DiracOpWilson class.
*/
/*--------------------------------------------------------------------*/

/***************************************************************************/
/*                                                                         */
/* wfm_mat_trick:                                                          */
/*                                                                         */
/* Do the color multiplication for each forward half spinor af with the    */
/* gauge field. Then expand the half spinors to full spinors (trick)       */
/* add them and store their sum into the full spinor chi.                  */
/*                                                                         */
/* chi      --> 4 component spinors (not padded)                           */
/* u        --> gauge fields (not padded)                                  */
/* af       --> 2 component spproj'ed padded spinors for the backward      */
/*              direction                                                  */
/* sign     --> sign in spin project                                       */
/* wilson_p --> pointer to Wilson structure                                */
/* cb       --> chekerboard                                                */
/*                                                                         */
/* WARNING:                                                                */
/*                                                                         */
/* This set of routines will work only if the node sublattices have        */
/* even number of sites in each direction.                                 */
/*                                                                         */
/***************************************************************************/

CPS_END_NAMESPACE
#include <util/data_types.h>
#include <util/vector.h>
#include <util/wilson.h>
#include <util/error.h>
#include <comms/scu.h>
CPS_START_NAMESPACE



#define CHI(r,c,s)      *(chi_p +(r+2*(c+3*s)))
#define AF0(r,c,s)      *(af0 + af0_ptr_b +(r+2*(c+3*s)))
#define AF1(r,c,s)      *(af1 + af1_ptr_b +(r+2*(c+3*s)))
#define AF2(r,c,s)      *(af2 + af2_ptr_b +(r+2*(c+3*s)))
#define AF3(r,c,s)      *(af3 + af3_ptr_b +(r+2*(c+3*s)))

#define TMP0(r,c,s)     *(tmp0_trick +(r+2*(c+3*s)))
#define TMP1(r,c,s)     *(tmp1_trick +(r+2*(c+3*s)))
#define TMP2(r,c,s)     *(tmp2_trick +(r+2*(c+3*s)))
#define TMP3(r,c,s)     *(tmp3_trick +(r+2*(c+3*s)))

#define U(r,row,col,d)  *(u_p+(r+2*(row+3*(col+3*d))))


IFloat tmp0_trick[HALF_SPINOR_SIZE];
IFloat tmp1_trick[HALF_SPINOR_SIZE];
IFloat tmp2_trick[HALF_SPINOR_SIZE];
IFloat tmp3_trick[HALF_SPINOR_SIZE];


extern void wilson_dslash_csmat(Float *out, Float *u, Float *in, int mu);
extern void wilson_dslash_trick(double *outf, 
				double *in0, double *in1, double *in2, double *in3, 
				double sign, int accum);



void wfm_mat_trick(IFloat *chi, 
		   IFloat *u,
		   IFloat *af0, IFloat *af1, IFloat *af2, IFloat *af3,
		   IFloat sign, 
		   Wilson *wilson_p, 
		   int cb)
{
  char *fname = "wfm_mat_trick";
  int iyzt;
  int ix;
  int ip;
  int c;
  int s;
  int ixb;
  int sx;
  int bpad;
  int shift;
  int *ptr;
  int af0_ptr_a, af1_ptr_a, af2_ptr_a, af3_ptr_a;
  int af0_ptr_b, af1_ptr_b, af2_ptr_b, af3_ptr_b;
  IFloat *chi_p;
  IFloat *u_p;

  /*************************************************************************/
  /* Initializations                                                       */
  /*************************************************************************/

  /* Construct the base pointer for the array of pointers that holds the   */
  /* addresses of the af spinors.                                          */
  bpad  = 0;
  shift = 1;
  ptr = wilson_p->ptr + (shift + 2 * (bpad + 2 * cb)) * wilson_p->offset;

  /*************************************************************************/
  /* Loop over sites                                                       */
  /*************************************************************************/
  chi_p = chi;
  u_p = u;
  ip = 0;
  for(iyzt=0; iyzt <= wilson_p->yztmax; iyzt++){

    af0_ptr_a = ptr[ip++];
    af1_ptr_a = ptr[ip++];
    af2_ptr_a = ptr[ip++];
    af3_ptr_a = ptr[ip++];
    sx        = ptr[ip++];

    for(ix=0; ix <= sx; ix++){

      ixb = ix * BLOCK;
      af0_ptr_b = af0_ptr_a + ixb;
      af1_ptr_b = af1_ptr_a + ixb;
      af2_ptr_b = af2_ptr_a + ixb;
      af3_ptr_b = af3_ptr_a + ixb;


      /*********************************************************************/
      /* Color multiply the af half spinors with U_mu                      */ 
      /*********************************************************************/
      /* mu = 0 */
      wilson_dslash_csmat(tmp0_trick, u_p, (af0+af0_ptr_b), 0);

      /* mu = 1 */
      wilson_dslash_csmat(tmp1_trick, u_p, (af1+af1_ptr_b), 1);

      /* mu = 2 */
      wilson_dslash_csmat(tmp2_trick, u_p, (af2+af2_ptr_b), 2);

      /* mu = 3 */
      wilson_dslash_csmat(tmp3_trick, u_p, (af3+af3_ptr_b), 3);

      u_p = u_p + 72;


      /*********************************************************************/
      /* Expand the half spinors in to a full spinor (trick), where        */
      /* the spin projection was done with [1+sign*gamma].                 */
      /* Store the result into the full spinor (accum=0).                  */
      /*********************************************************************/
      wilson_dslash_trick(chi_p,
			  tmp0_trick,
			  tmp1_trick,
			  tmp2_trick,
			  tmp3_trick,
			  -sign, 0);

      chi_p = chi_p + 24;

    }
  }

}

CPS_END_NAMESPACE
