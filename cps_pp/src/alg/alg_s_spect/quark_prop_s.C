#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2013-04-05 17:46:30 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_s_spect/quark_prop_s.C,v 1.14 2013-04-05 17:46:30 chulwoo Exp $
//  $Id: quark_prop_s.C,v 1.14 2013-04-05 17:46:30 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: quark_prop_s.C,v $
//  $Revision: 1.14 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_s_spect/quark_prop_s.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
// quark_prop_s.C

CPS_END_NAMESPACE
#include <alg/quark_prop_s.h>
#include <util/smalloc.h>
#include <util/gjp.h>
#include <util/vector.h>	// uDotXEqual()
#include <util/verbose.h>
#include <util/error.h>
#include <string.h>			// memcpy()
#ifdef PARALLEL
#include <comms/sysfunc_cps.h>
#endif
#include <alg/myenum.h>
#include <util/qcdio.h>
CPS_START_NAMESPACE

//------------------------------------------------------------
//  class static members initializations
//------------------------------------------------------------

char QuarkPropS::cname[] = "QuarkPropS";
Float* QuarkPropS::qSrc = 0;	
int QuarkPropS::prop_count = 0;

char QuarkPropSMng::cname[] = "QuarkPropSMng";
int QuarkPropSMng::isInit = 0;
int* QuarkPropSMng::slice = 0;
Float ***QuarkPropSMng::qtab = 0;

//------------------------------------------------------------
//  static variables used inside this file
//------------------------------------------------------------

static int nx[4];

static int vsize;
static int node_origin[4];
static int node_end[4];
        // used when setting up quark source. 
        // and node_origin[] is the global (x,y,z,t)
	// coordinates of the origin of this node 

//------------------------------------------------------------
// functions local to this file
//------------------------------------------------------------

inline int max(int a, int b) { return a > b ? a : b; }
inline int min(int a, int b) { return a < b ? a : b; }

void getNodeOriginEnd()
{
   node_origin[0] =  GJP.XnodeCoor() * nx[0];
   node_origin[1] =  GJP.YnodeCoor() * nx[1];
   node_origin[2] =  GJP.ZnodeCoor() * nx[2];
   node_origin[3] =  GJP.TnodeCoor() * nx[3];
   for (int i = 0 ; i < 4; ++i) {
     node_end[i] = node_origin[i] + nx[i] - 1; 
   }
}

//------------------------------------------------------------
// private function: index the staggered fermion vector correctly
//------------------------------------------------------------

int QuarkPropS::X_OFFSET(const int *x)
{
  return lat.FsiteOffsetChkb(x) * VECT_LEN
 	 + ((x[0]+x[1]+x[2]+x[3]) & 1 ? vsize : 0);
}

//------------------------------------------------------------
//  qSrc is allocated here 
//------------------------------------------------------------

void QuarkPropS::setupQuarkPropS() 
{
    char *fname = "setupQuarkPropS()";
    //--------------------------------------------
    // need to do class initializations?
    //--------------------------------------------

    if(prop_count == 0)
    {
	nx[0] = GJP.XnodeSites();
        nx[1] = GJP.YnodeSites();
        nx[2] = GJP.ZnodeSites();
        nx[3] = GJP.TnodeSites();
 
	vsize = GJP.VolNodeSites() * VECT_LEN / 2;
	qSrc = (Float *)smalloc(2*vsize*sizeof(Float)); 
        if(qSrc == 0)
          ERR.Pointer(cname,fname, "qSrc");
        VRB.Smalloc(cname,fname, "qSrc", qSrc, 2*vsize*sizeof(Float));

	getNodeOriginEnd();
    }


    //--------------------------------------------
    // register this quark propagator
    //--------------------------------------------
    QuarkPropSMng::qadd(this);

    //--------------------------------------------
    // increase the propagator counter
    //--------------------------------------------
    ++prop_count;

}


//------------------------------------------------------------
// CTOR 
//------------------------------------------------------------

QuarkPropS::QuarkPropS(Lattice& lattice, StagQuarkArg& arg) 
: qid(arg.qid), qarg(arg), lat(lattice)
{ VRB.Func(cname,"QuarkPropS(Lattice& , StagQuarkArg& )"); }

//------------------------------------------------------------
// DTOR free the source vector if out of scope
//------------------------------------------------------------

QuarkPropS::~QuarkPropS() 
{
  char *fname = "~QuarkPropS()";
  VRB.Func(cname,"~QuarkPropS()"); 
  prop_count--;
  if(prop_count == 0) {
    VRB.Sfree(cname,fname, "qSrc", qSrc);
    sfree(qSrc);
  }
}

//------------------------------------------------------------
// free the quark propagator memory 
//------------------------------------------------------------
void QuarkPropS::destroyQuarkPropS(int id)
{
  QuarkPropSMng::destroyQuarkPropS(id);
}
	
//------------------------------------------------------------
//  set up Point source
//------------------------------------------------------------

void QuarkPropS::setPntSrc(const int *src_p, int color)
{ 
    VRB.Func(cname,"setPntSrc(const int *, int)"); 

    //-----------------------------------------
    // zero source vector first
    //-----------------------------------------
    Float *src_tmp = qSrc;

    for (int i = 0; i < 2*vsize; ++i) {
      *src_tmp++ = 0.0;  
    }

    //-----------------------------------------
    // determine if the POINT source is on this node
    //-----------------------------------------
    int site[4] = { src_p[0], src_p[1],src_p[2],src_p[3] };
    int isOnNode = 1;

    for (int j = 0; j < 4; j++) {
      isOnNode *= (site[j] >= node_origin[j]);
      isOnNode *= (site[j] <= node_end[j]);
    }
    //-----------------------------------------
    // calculate the offset and set the source
    //-----------------------------------------
    if (isOnNode) {
      for (int k = 0; k < 4; k++) {
        site[k] -= node_origin[k];
      }
      //int soff = X_OFFSET(site);
      int soff = site[0]+GJP.XnodeSites()*(site[1]+GJP.YnodeSites()*(site[2]+GJP.ZnodeSites()*(site[3])));
      *(qSrc+VECT_LEN*soff+2*color) = 1.0;
    }

}

//------------------------------------------------------------
// if there is overlap with src, return 1;
// else return 0;
//------------------------------------------------------------

static int hasOverlapSrc(const int *src_origin, const int *src_end) 
{
    //---------------------------------------
    // No overlap with source on this node
    //---------------------------------------
    for (int i = 0; i < 4; ++i) {
      if ( src_end[i] < node_origin[i] || 
	   src_origin[i] > node_end[i])
	return 0;   
    }

    //---------------------------------------
    // there IS overlap with the source
    //---------------------------------------
    return 1;
}


   
//------------------------------------------------------------
// set up a wall source: Z-wall or 2Z-wall
//------------------------------------------------------------

void QuarkPropS::setWallSrc(Matrix **gm, StagQuarkSrc& qs, int color )
{
    char *fname = "setWallSrc(const Float **, StagQuarkSrc&, int)";

    VRB.Func(cname,"setWallSrc(const Float *, QuarkSrc&, int)"); 
    int i;

    //-----------------------------------------
    // check if the source is valid
    //-----------------------------------------
    int wall = 0;
    int count = 0; 
    for (i = 0; i < 4; ++i) {
      if ( qs.origin[i] == qs.end[i] ) {
        wall = i;
	count++;
      }
    }
    if ( count != 1 ) {
	ERR.General(cname, fname, "Not a 3D-wall source\n");
    }
    if ( qs.dir != wall) {
	ERR.General(cname, fname, "Wall source direction not correct\n");
    }

    //-----------------------------------------
    // zero source vector first
    //-----------------------------------------
    Float *src_tmp = qSrc;

    for (i = 0; i < 2*vsize; ++i) {  
      *src_tmp++ = 0.0;
    }	


    //-----------------------------------------------
    // determine if source has overlap with this node
    // if yes, local coordinates of the source are 
    // calculated
    //-----------------------------------------------
    if (hasOverlapSrc(qs.origin, qs.end)) {   

      int src_origin[4], src_end[4];
      for (i = 0; i < 4; ++i) {
        src_origin[i] = max(qs.origin[i], node_origin[i]) - node_origin[i];
        src_end[i]    = min(qs.end[i], node_end[i]) - node_origin[i];
      } // ASSERT: src_origin[wall] = src_end[wall]
    
      //-------------------------------------------------------
      // walk through the 3D world and set the sources 
      //-------------------------------------------------------
      int dir = qs.dir;
      int dim[4] = { nx[0], nx[1], nx[2], nx[3] };
      dim[dir] = 1;

      int i = (dir + 1)%4;
      int j = (dir + 2)%4;
      int k = (dir + 3)%4;

      int s[4];
      s[dir] = src_origin[dir];

      const Float *g = (Float *)gm[s[dir]] + color*VECT_LEN;
      // Note:
      // base offset to the source slice and color offset
      // CAUTION:
      // g_dagger * delta is
      //
      // | a11+ib11 a12+ib12 a13+ib13|^dagger 	(1)	(0)	(0)
      // | a21+ib21 a22+ib22 a23+ib23| 		(0) or 	(1) or 	(0)
      // | a31+ib31 a32+ib32 a33+ib33| 		(0)	(0)	(1)

      for(s[i] = src_origin[i]; s[i] <= src_end[i]; s[i]++)
        for(s[j] = src_origin[j]; s[j] <= src_end[j]; s[j]++)
          for(s[k] = src_origin[k]; s[k] <= src_end[k]; s[k]++)
          {
	    //--------------------------------------------------
	    // if this is WALLZ source
	    // or WALL2Z source with all 3d coordinates even
	    // need to set the source on this site
	    //--------------------------------------------------
	    if ( qs.type == WALLZ || ( qs.type == WALL2Z)&&
	       ( (s[i]%2 == qs.origin[i]%2) 
		 && (s[j]%2 == qs.origin[j]%2) 
		 && (s[k]%2 == qs.origin[k]%2))) {

              int s3d[4] = { s[0], s[1], s[2], s[3] };
              s3d[dir] = 0;

	      //------------------------------------------------
	      // set the source: g(s)^dagger * delta(s)
	      //------------------------------------------------
              const Float *gp = g + MATRIX_SIZE * g_offset(s3d, dim);
	      Float *src = qSrc + X_OFFSET(s);

	      for(int color = 0; color < 3; color++) {
	        *src++ = *gp++;
	        *src++ = -*gp++;
	      }
	    }
     	  } 
    } 
}

//------------------------------------------------------------
// Calculate all three columns(color) of quark propagator, and
// stored in 3 vectors prop[3]. Three cases are handled here:
//
// 1.  POINT source is used 
//
// 2.  Wall sources used 
//        (D + 2m) G'(x) = delta(x)_Coulomb
//
//     we get G'(x) = g(x)^dagger * G(x)^Coulomb.
// 
//     For local hadron operator using G'(x) is fine
//     because 
//	  
//	  G'(x)G'(x)^dag == G(x)G(x)^dag (meson)
//        e_ijk * g_i(x)*g_j(x)*g_k(x) = det(g) = 1 (nucleon)
// 
// 3.  Wall source used.
//     Solution G(x)^Coulomb is needed because for non-local 
//     hadron operators the identities in case 2 are not true.
//------------------------------------------------------------

void QuarkPropS::getQuarkPropS(char *results)
{
#if TARGET==cpsMPI
    using MPISCU::fprintf;
#endif
    char *fname = "getQuarkPropS(const char *)";
    VRB.Func(cname, fname);

    int cg_iter;
    Float true_res;
    FILE *fp;		// monitoring info of CG

    for (int color = 0; color < 3; ++color)  {
      if ( qarg.src.type == S_QUARK_POINT ) {
	setPntSrc(qarg.src.origin, color);
      }
      else {
	Matrix **gm = lat.FixGaugePtr();
	setWallSrc(gm, qarg.src, color);
      }

      IFloat *src = (IFloat *)qSrc;
      //lat.Convert(CANONICAL,(Vector*)src);
      IFloat *sln = (IFloat *)prop[color];
      cg_iter = lat.FmatInv((Vector *)sln, (Vector *)src, 
			    &(qarg.cg), &true_res, CNV_FRM_YES);
			    //&(qarg.cg), &true_res, CNV_FRM_NO);

      // Added for anisotropic lattices
      vecTimesEquFloat(sln, GJP.XiBare()/GJP.XiV(), 2*vsize); 
      // End modification 

      // Print out monitor info of the inversion
      //---------------------------------------------------------------
      if(results != 0){
        if( NULL == (fp = Fopen(results, "a")) ) {
          ERR.FileA(cname,fname, results);
      	}
    	Fprintf(fp,"%e %d %e\n", IFloat(qarg.cg.mass), cg_iter, 
				 IFloat(true_res));
    	Fclose(fp);
      }
 
      //------------------------------------------------------
      // if it is non-local, realize case 3
      //------------------------------------------------------
      if (qarg.sln == NONLOCAL) {
	Matrix **agm = lat.FixGaugePtr();
        coulomb(sln, agm, qarg.src.dir); 
      }
    }
}

//-------------------------------------------------------------
//  x = g(n1, n2, n3, n4 ) * x(n1, n2, n3, n4)
//-------------------------------------------------------------

void QuarkPropS::coulomb(IFloat *x, Matrix **g, int dir)
{
  char *fname = "coulomb(IFloat *, Matrix **, int)";
  VRB.Func(cname, fname);

  Vector *vr_tmp = (Vector *)smalloc(VECT_LEN * sizeof(IFloat));
  if(vr_tmp == 0)
    ERR.Pointer(cname,fname, "vr_tmp");
  VRB.Smalloc(cname,fname, "vr_tmp", vr_tmp, VECT_LEN * sizeof(IFloat));

  int dim[4] = { nx[0], nx[1], nx[2], nx[3] };
  dim[dir] = 1;

  int i = (dir + 1)%4;
  int j = (dir + 2)%4;
  int k = (dir + 3)%4;

  int s[4];
  for(s[dir] = 0; s[dir] < nx[dir]; s[dir]++) {
    Matrix *g_base = g[s[dir]];

    for(s[i] = 0; s[i] < nx[i]; s[i]++) 
      for(s[j] = 0; s[j] < nx[j]; s[j]++) 
        for(s[k] = 0; s[k] < nx[k]; s[k]++) 
        {
          int s3d[4] = { s[0], s[1], s[2], s[3] };
          s3d[dir] = 0;
          Matrix *gp = g_base + g_offset(s3d, dim);
          IFloat *xp = x + X_OFFSET(s);
          memcpy((IFloat *)vr_tmp, xp, VECT_LEN * sizeof(IFloat));
	  uDotXEqual(xp, (IFloat *)gp, (IFloat *)vr_tmp);
        }
  }
  VRB.Sfree(cname,fname, "vr_tmp",vr_tmp);
  sfree(vr_tmp);
}

//============================================================
// Functions of QuarkPropSMng class
//============================================================
QuarkPropSMng::QuarkPropSMng()
{
  char *fname = "QuarkPropSMng()";
  VRB.Func(cname, fname); 
  //--------------------------------------------
  // class initialization
  //--------------------------------------------
  if (!isInit) {
    qtab = (Float ***)smalloc(MAXNUMQP*sizeof(Float **));
    if(qtab == 0)
      ERR.Pointer(cname,fname, "qtab");
    VRB.Smalloc(cname,fname, "qtab", qtab, MAXNUMQP*sizeof(Float **));

    slice = (int *)smalloc(MAXNUMQP*sizeof(int));
    if(slice == 0)
      ERR.Pointer(cname,fname, "slice");
    VRB.Smalloc(cname,fname, "slice", slice, MAXNUMQP*sizeof(int));
    for (int j = 0; j < MAXNUMQP; j++) {
      qtab[j] = 0;
      slice[j] = 0;
    }
    isInit = 1;
  }
}

//------------------------------------------------------------
//  DTOR
//------------------------------------------------------------
QuarkPropSMng::~QuarkPropSMng()
{
  char *fname = "~QuarkPropSMng()";
  VRB.Func(cname,fname);
 
  QuarkPropSMng::destroyQuarkPropS();
  VRB.Sfree(cname,fname, "qtab",qtab);
  sfree(qtab);
}

//------------------------------------------------------------
// allocate memory and register the quark propagator
//------------------------------------------------------------
void QuarkPropSMng::qadd(QuarkPropS *qp)
{
  char *fname = "qadd(QuarkPropS *)";
  VRB.Func(cname, fname); 
  //--------------------------------------------
  // Allocate memory for this propagator
  //--------------------------------------------
  qp->prop = (Float **)smalloc(3*sizeof(Float *));
  if(qp->prop == 0)
    ERR.Pointer(cname,fname, "qp->prop");
  VRB.Smalloc(cname,fname, "qp->prop", qp->prop, 3*sizeof(Float *));

  int color;
  for (color = 0; color < 3; ++color) {
    qp->prop[color] = (Float *)smalloc(2*vsize*sizeof(Float)); 
    if(qp->prop[color] == 0)
      ERR.Pointer(cname,fname, "qp->prop[color]");
    VRB.Smalloc(cname,fname, "qp->prop[color]", qp->prop[color], 2*vsize*sizeof(Float));
  }

  //--------------------------------------------
  // always zero the propagator when registering
  //--------------------------------------------
  for (color = 0; color < 3; ++color) {
    Float *tmp = qp->prop[color];
    for(int i=0; i<2*vsize; ++i) {
      *tmp++ = 0;
    }
  }

  //--------------------------------------------
  // register prop and keep the prop ptr
  //--------------------------------------------
  int index = qp->qid;
  qtab[index] = qp->prop;
  slice[index] = qp->qarg.src.origin[qp->qarg.src.dir];

  VRB.Flow(cname,fname, "QuarkPropS of id = %d is registered\n", index);
}

//------------------------------------------------------------
// free the quark propagator of qid = id;
// By default (id = -1) all quark propagators are freed.
//------------------------------------------------------------
void QuarkPropSMng::destroyQuarkPropS(int id)
{
  if (id == -1) {
    for (int i = 0; i < MAXNUMQP; i++) {
      release(i);
    }
  }
  else 
    release(id);
}

//------------------------------------------------------------
// release the quark propagator of qid = id
//------------------------------------------------------------
void QuarkPropSMng::release(int id) 
{
  char *fname = "release(int id)";
  VRB.Func(cname, fname); 
  Float **tmp_prop = qtab[id];
  if (tmp_prop) {	
    for (int color = 0; color < 3; ++color) {
      if (tmp_prop[color]) {
        VRB.Sfree(cname,fname, "tmp_prop[color]",tmp_prop[color]);
	sfree(tmp_prop[color]);
      }
    }	
    VRB.Sfree(cname,fname, "tmp_prop", tmp_prop);
    sfree(tmp_prop); 
    qtab[id] = 0;	// turn off the entry light
    VRB.Flow(cname,fname,"QuarkPropS of id = %d is destroyed\n", id);
  }
}


CPS_END_NAMESPACE
