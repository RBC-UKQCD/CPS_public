#include<config.h>
CPS_START_NAMESPACE
/****************************************************************************/
/* 10/6/96                                                                  */
/*                                                                          */
/* wilson_int:                                                              */
/*                                                                          */
/* This routine performs all initializations needed before wilson func      */
/* are called. It sets the addressing related arrays and reserves memory    */
/* for the needed temporary buffers. It only needs to be called             */
/* once at the begining of the program (or after a wilson_end call)         */
/* before any number of calls to wilson funcs are made.                     */
/*                                                                          */
/****************************************************************************/

/*--------------------------------------------------------------------------*/
/* Include header files                                                     */
/*--------------------------------------------------------------------------*/
CPS_END_NAMESPACE
#include<util/wilson.h>
#include<util/gjp.h>
#include<util/smalloc.h>
#include<util/verbose.h>
#include<util/error.h>
CPS_START_NAMESPACE




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
  int size;



/*--------------------------------------------------------------------------*/
/* Set sublattice direction sizes                                           */
/*--------------------------------------------------------------------------*/
  slx = GJP.XnodeSites();
  sly = GJP.YnodeSites();
  slz = GJP.ZnodeSites();
  slt = GJP.TnodeSites();

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

/*--------------------------------------------------------------------------*/
/* Reserve memory for the node sublattice pointers                          */
/*--------------------------------------------------------------------------*/
  size = 40*sly*slz*slt*sizeof(int);
  wilson_p->ptr = (int *) smalloc(size);
  if( wilson_p->ptr == 0)
    ERR.Pointer(cname,fname, "ptr");
  VRB.Smalloc(cname,fname,
	      "ptr", wilson_p->ptr, size);

/*--------------------------------------------------------------------------*/
/* Set the node sublattice pointers                                         */
/*--------------------------------------------------------------------------*/
  wfm_sublatt_pointers(slx, sly, slz, slt, slatt_eo, wilson_p);


/*--------------------------------------------------------------------------*/
/* Reserve memory for 2  temporary spinors (nedded by m, mdag and mdagm)    */
/* Memory is allocated into two consecutive chunks. The above routines      */
/* will have to set the pointer of the first temporary to the base          */
/* and of the second in the middle.                                         */
/* WARNING: valid for "even" sublattices only                               */
/*--------------------------------------------------------------------------*/
  spinor_words = 2 * SPINOR_SIZE * wilson_p->vol[0];

  wilson_p->spinor_tmp = (IFloat *) smalloc(spinor_words*sizeof(IFloat));
  if(wilson_p->spinor_tmp == 0)
    ERR.Pointer(cname,fname, "spinor_tmp");
  VRB.Smalloc(cname,fname,
	      "spinor_tmp", wilson_p->spinor_tmp, spinor_words*sizeof(IFloat));
    

/*--------------------------------------------------------------------------*/
/* Reserve memory for the 4 forward and 4 backward spin projected half      */ 
/* spinors.                                                                 */
/*--------------------------------------------------------------------------*/
  for(i=0; i<4; i++){
    half_spinor_words = HALF_SPINOR_SIZE * wilson_p->padded_subgrid_vol[i];

    wilson_p->af[i] = (IFloat *) smalloc(half_spinor_words*sizeof(IFloat));
    if(wilson_p->af[i] == 0)
      ERR.Pointer(cname,fname, "af[i]");
    VRB.Smalloc(cname,fname,
		"af[i]", wilson_p->af[i], half_spinor_words*sizeof(IFloat));

    wilson_p->ab[i] = (IFloat *) smalloc(half_spinor_words*sizeof(IFloat));
    if(wilson_p->ab[i] == 0)
      ERR.Pointer(cname,fname, "ab[i]");
    VRB.Smalloc(cname,fname,
		"ab[i]", wilson_p->ab[i], half_spinor_words*sizeof(IFloat));
  }

}
























CPS_END_NAMESPACE
