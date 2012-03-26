#include<config.h>
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2012-03-26 13:50:11 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_s_spect/hadron_prop_s.C,v 1.13 2012-03-26 13:50:11 chulwoo Exp $
//  $Id: hadron_prop_s.C,v 1.13 2012-03-26 13:50:11 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: hadron_prop_s.C,v $
//  $Revision: 1.13 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_s_spect/hadron_prop_s.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
// hadron_prop_s.C

#ifndef lint
static char vcid[] = "$Id: hadron_prop_s.C,v 1.13 2012-03-26 13:50:11 chulwoo Exp $";
#endif /* lint */

#include <util/qcdio.h>
#include <alg/hadron_prop_s.h>
#include <util/rcomplex.h>
#include <comms/glb.h>
#include <comms/sysfunc_cps.h>
#include <alg/myenum.h>
CPS_START_NAMESPACE

const int MAX_LEN = 1023;     // max. number of IFloats that can be transfered

char HadronPropS::cname[] = "HadronPropS";

static int nx[4];
static int Nx[4];
static int vsize;
static BndCndType bc[4]; 

//------------------------------------------------------------------
// protected function used in all derived classes
//------------------------------------------------------------------
int HadronPropS::X_OFFSET(const int *x)
{ return lat.FsiteOffsetChkb(x) * VECT_LEN
           + ((x[0]+x[1]+x[2]+x[3]) & 1 ? vsize : 0); }


//------------------------------------------------------------------
// non-pure virtual function, override in NLocalProp class
//------------------------------------------------------------------
void HadronPropS::getNeighbors(int t) {}

//------------------------------------------------------------------
//  CTOR & DTOR
//------------------------------------------------------------------
HadronPropS::HadronPropS( Lattice &lattice, int num, 
		        int direction, int srcslice, int incr )
: lat(lattice),
  dir(direction),
  slice(srcslice),
  n_props(num), 
  stride(incr)
{
    char *fname = "HadronPropS(Lattice &, int, int)";
    nx[0] = GJP.XnodeSites();
    nx[1] = GJP.YnodeSites();
    nx[2] = GJP.ZnodeSites();
    nx[3] = GJP.TnodeSites();

    Nx[0] = GJP.Xnodes();
    Nx[1] = GJP.Ynodes();
    Nx[2] = GJP.Znodes();
    Nx[3] = GJP.Tnodes();

    vsize = GJP.VolNodeSites() * VECT_LEN / 2;

    bc[0] = GJP.Xbc();
    bc[1] = GJP.Ybc(); 
    bc[2] = GJP.Zbc(); 
    bc[3] = GJP.Tbc(); 

    prop = (Complex *)smalloc(nx[dir] * num * sizeof(Complex));
    if(prop == 0)
      ERR.Pointer(cname,fname, "prop");
    VRB.Smalloc(cname,fname, "prop", prop, nx[dir]*num*sizeof(Complex));
}

HadronPropS::~HadronPropS() 
{ 
  char *fname = "~HadronPropS()";
  VRB.Sfree(cname,fname, "prop",prop);
  sfree(prop); 
}

BndCndType HadronPropS::bcd()
{
  return bc[dir];
}

int HadronPropS::propLen() const 
{ return CR * n_props * nx[dir]; }

int HadronPropS::propLenTotal() const 
{ return CR * n_props * nx[dir] * Nx[dir]; }

int HadronPropS::propLenFold() const 
{ return n_props * (nx[dir] * Nx[dir]/2 + 1); }
//------------------------------------------------------------------
//  Compute the contribution from each site and sum over a
//  hyperplane. For nonlocal operators scu transfers are needed
//  for each hyperplane.
//------------------------------------------------------------------

void HadronPropS::getHadronPropS()
{
  int s[4];
  for (int len = 0; len < n_props * nx[dir]; len++)
    prop[len] = 0;

  Complex *currt_p = prop; 

  //---------------------------------------------------------------
  // sum on slice orthogonal to dir (not necessary time)
  // i, j, k are the rest of directions
  //---------------------------------------------------------------

  int i = (dir+1)%4;
  int j = (dir+2)%4;
  int k = (dir+3)%4;

  //---------------------------------------------------------------
  //  For each slice s[dir], sum over the hyperplane on this node 
  //---------------------------------------------------------------

  for (s[dir] = 0; s[dir] < nx[dir]; s[dir]++) { 
    getNeighbors(s[dir]);
    for (s[i] = 0; s[i] < nx[i]; s[i] += stride) {    
      for (s[j] = 0; s[j] < nx[j]; s[j] += stride) { 
        for (s[k] = 0; s[k] < nx[k]; s[k] += stride) {

	  localVal(currt_p, s);   
        }
      }
    }
    currt_p += n_props;	
  }

  //---------------------------------------------------------------
  // sum over all sites on all slices 
  //---------------------------------------------------------------

      slice_sum((Float *)prop, propLen(), dir);
      
      return;
}

void HadronPropS::download_prop(HadronType type, Float *buf)
{
  const Float *tmp = propPtr();
  
  
  collect_prop(type, buf, tmp, n_props * CR, dir, slice);
  return;
}


#if (defined PARALLEL) && (TARGET == QCDOC)
//---------------------------------------------------------------------
// download data from all nodes in the propagator direction to
// node 0 (origin of the lattice)
//---------------------------------------------------------------------

const SCUDir pos_dir[] = { SCU_XP, SCU_YP, SCU_ZP, SCU_TP };
const SCUDir neg_dir[] = { SCU_XM, SCU_YM, SCU_ZM, SCU_TM };

static int isNodeOrigin();
static void storeData(IFloat *buf, const IFloat *term, int len);

void HadronPropS::collect_prop(HadronType type, Float *sum_buf, 
			      const Float *IFloat_p, int unit, 
			      int direction, int t)
{
  char *fname ="collect_prop(HadronType, Float, Float, *IFloat_p, int, int,int)";

  int blcklength = unit*nx[direction]*sizeof(Float);

  IFloat *receive_buf = (IFloat *)smalloc(blcklength*Nx[direction]);

  // added by manke for check

  if(receive_buf == 0)
    ERR.Pointer(cname,fname, "receive_buf");
  VRB.Smalloc(cname,fname, "receive_buf", receive_buf, blcklength*Nx[direction]);

  // end added

  IFloat *receive_buf_p = receive_buf;
  IFloat *transmit_buf_p;
  
  int itmp;		  
	// loop index with 1<= itmp < Nx[i] 

  //--------------------------------------------------------------
  // store data on node 0 first 
  //--------------------------------------------------------------
  if (isNodeOrigin()) {
    memcpy(receive_buf_p, IFloat_p, blcklength);
    receive_buf_p += blcklength/sizeof(Float);
  }
  //--------------------------------------------------------------
  // address of buffer of to be sent first (data on this node) 
  //--------------------------------------------------------------
  transmit_buf_p = (IFloat *)IFloat_p; 

#if USE_QMP
  ERR.General(cname,fname,"Not implemented with QMP\n");
#else
  SCUDirArg send(transmit_buf_p, neg_dir[direction], SCU_SEND, blcklength);
  SCUDirArg recv(receive_buf_p, pos_dir[direction], SCU_REC, blcklength);

  //--------------------------------------------------------------
  // tranmit & receive Nx[direction] - 1 times in snd_dir[direction] direction
  //--------------------------------------------------------------
  for ( itmp = 1; itmp < Nx[direction]; itmp++) {


	 //-----------------------------------------------------------
         // do SCU transfers
	 //-----------------------------------------------------------
         send.StartTrans();
         recv.StartTrans();
         send.TransComplete();
         recv.TransComplete();

	 //-----------------------------------------------------------
         // the received data will be sent out     
	 // the free buffer will be used to receive      
	 //-----------------------------------------------------------
	 send.Addr(transmit_buf_p = receive_buf_p);
	 recv.Addr(receive_buf_p += blcklength/sizeof(Float) );
  }
#endif





  // everything is in receive_buf now

  if(isNodeOrigin()) {
    int i;
    for(i = unit * t; i < unit * nx[direction] * Nx[direction]; i++)
	*sum_buf++ += receive_buf[i];
    if( (type != SMESON) && (type != SMOMMESON) && (bc[dir] == BND_CND_APRD) ) {
      for(i = 0; i < unit * t; i++)
	*sum_buf++ -= receive_buf[i];
    }
    else {
      for(i = 0; i < unit * t; i++)
	*sum_buf++ += receive_buf[i];
    }
  }

  sfree(receive_buf);
}



int isNodeOrigin() 
{
  return (GJP.XnodeCoor() == 0) && (GJP.YnodeCoor() == 0) && 
	 (GJP.ZnodeCoor() == 0) && (GJP.TnodeCoor() == 0);
}

void storeData(IFloat *buf, const IFloat *term, int len)
{
  for(int i = 0; i < len; i++) {
    *buf++ = *term++;
  }
}

#else

void HadronPropS::collect_prop(HadronType type, Float *sum_buf, 
	      const Float *IFloat_p, int unit, int direction, int t)
{
    int i;
    for(i = unit * t; i < unit * nx[direction] * Nx[direction]; i++)
        *sum_buf++ += IFloat_p[i];
    if( (type != SMESON) && (type != SMOMMESON) && (bc[dir] == BND_CND_APRD) ) {
      for(i = 0; i < unit * t; i++)
        *sum_buf++ -= IFloat_p[i];
    }
    else {
      for(i = 0; i < unit * t; i++)
        *sum_buf++ += IFloat_p[i];
    }
}

#endif


CPS_END_NAMESPACE
