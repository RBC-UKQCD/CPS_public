#include <config.h>
CPS_START_NAMESPACE
/*--------------------------------------------------------------------*/
/*! \file
  \brief  Routine used internally in the DiracOpWilson class.
*/
/*--------------------------------------------------------------------*/

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

/***************************************************************************/
/* Function declarations                                                   */
/***************************************************************************/
void wfm_comm_forward(IFloat *ab0, IFloat *ab1, IFloat *ab2, IFloat *ab3,
		      Wilson *wilson_p);

void wfm_comm_backward(IFloat *af0, IFloat *af1, IFloat *af2, IFloat *af3,
		       Wilson *wilson_p);

void wfm_spproj(IFloat *af0, IFloat *af1, IFloat *af2, IFloat *af3,
		IFloat *psi, 
		IFloat sign, 
		Wilson *wilson_p, 
		int cb);

void wfm_cmat_spproj(IFloat *ab0, IFloat *ab1, IFloat *ab2, IFloat *ab3,
		     IFloat *u, 
		     IFloat *psi, 
		     IFloat sign, 
		     Wilson *wilson_p, 
		     int cb);

void wfm_mat_trick(IFloat *chi, 
		   IFloat *u,
		   IFloat *af0, IFloat *af1, IFloat *af2, IFloat *af3,
		   IFloat sign, 
		   Wilson *wilson_p, 
		   int cb);

void wfm_trick(IFloat *chi, 
	       IFloat *ab0, IFloat *ab1, IFloat *ab2, IFloat *ab3,
	       IFloat sign, 
	       Wilson *wilson_p, 
	       int cb);






void wilson_dslash(IFloat *chi, 
		   IFloat *u, 
		   IFloat *psi, 
		   int cb,
		   int dag,
		   Wilson *wilson_p)
{
  char *fname = "wilson_dslash";
  IFloat sign;
  IFloat *u_eo[2];
  IFloat *af0;
  IFloat *af1;
  IFloat *af2;
  IFloat *af3;
  IFloat *ab0;
  IFloat *ab1;
  IFloat *ab2;
  IFloat *ab3;
  int cbn;
  

/*--------------------------------------------------------------------------*/
/* Initializations                                                          */
/*--------------------------------------------------------------------------*/
  /* For dag 0/1 set sign to +/- 1 */
  if(dag == 0)
    sign = 1.0;
  else if(dag == 1)
    sign = -1.0;
  else{
    ERR.General(" ",fname,"dag must be 0 or 1");
  }

  /* set cbn = !cb */
  if(cb == 0)
    cbn = 1;
  else if(cb == 1)
    cbn = 0;
  else{
    ERR.General(" ",fname,"cb must be 0 or 1");
  }


  /* Set the pointers for the gauge field on the even and odd checkerboards */
  u_eo[0] = u;
  u_eo[1] = u + GAUGE_SIZE * wilson_p->vol[0];

  /* Set the pointer to the pointer array pointing to the 4 forward half    */
  /* spinors af.                                                            */
  af0 = wilson_p->af[0];
  af1 = wilson_p->af[1];
  af2 = wilson_p->af[2];
  af3 = wilson_p->af[3];

  /* Set the pointer to the pointer array pointing to the 4 backward half   */
  /* spinors ab.                                                            */
  ab0 = wilson_p->ab[0];
  ab1 = wilson_p->ab[1];
  ab2 = wilson_p->ab[2];
  ab3 = wilson_p->ab[3];


  /*------------------------------------------------------------------------*/
  /* Spin project psi with (1-sign*gamma) into the 4 forward half           */
  /* spinors af[4]                                                          */
  /*------------------------------------------------------------------------*/
  wfm_spproj(af0, af1, af2, af3, psi, sign, wilson_p, cb);

  /*------------------------------------------------------------------------*/
  /* Communicate backwards the 4 forward half spinors af[4]                 */
  /*------------------------------------------------------------------------*/
  wfm_comm_backward(af0, af1, af2, af3, wilson_p);

  /*------------------------------------------------------------------------*/
  /* Spin project psi with (1+sign*gamma) and do the color multiplication   */
  /* of the resulting half spinors with the complex conjugate gauge field.  */
  /* Put the result into the 4 backward half spinors ab[4].                 */
  /*------------------------------------------------------------------------*/
  wfm_cmat_spproj(ab0, ab1, ab2, ab3, u_eo[cb], psi, sign, wilson_p, cb);

  /*------------------------------------------------------------------------*/
  /* Communicate forward the 4 backward half spinors ab[4]                  */
  /*------------------------------------------------------------------------*/
  wfm_comm_forward(ab0, ab1, ab2, ab3, wilson_p);

  /*------------------------------------------------------------------------*/
  /* Do the color multiplication for each forward half spinor af with the   */
  /* gauge field. Then expand the half spinors to full spinors (trick)      */
  /* add them and store their sum into the full spinor chi.                 */
  /*------------------------------------------------------------------------*/
  wfm_mat_trick(chi, u_eo[cbn], af0, af1, af2, af3, sign, wilson_p, cb);

  /*------------------------------------------------------------------------*/
  /* Expand the backward half spinors ab to full spinors (trick) add them   */
  /*and then add their sum into the full spinor chi.                        */
  /*------------------------------------------------------------------------*/
  wfm_trick(chi, ab0, ab1, ab2, ab3, sign, wilson_p, cb);

}

CPS_END_NAMESPACE
