#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief Definition of FstagTypes methods.

  $Id: f_stag_t.C,v 1.5 2004-02-16 13:21:42 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-02-16 13:21:42 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/f_stag_types/f_stag_t.C,v 1.5 2004-02-16 13:21:42 zs Exp $
//  $Id: f_stag_t.C,v 1.5 2004-02-16 13:21:42 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: f_stag_t.C,v $
//  $Revision: 1.5 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/f_stag_types/f_stag_t.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// f_stag_types.C
//
// FstagTypes is derived from Lattice and is relevant to
// all fermion classes with Staggered type fermions 
// These classes are derived from FstagTypes
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/enum.h>
#include <util/verbose.h>
#include <util/gjp.h>
#include <util/error.h>
#include <comms/cbuf.h>
CPS_START_NAMESPACE

//------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------
FstagTypes::FstagTypes()
{
  cname = "FstagTypes";
  char *fname = "FstagTypes()";
  VRB.Func(cname,fname);

    e_vsize = VECT_LEN/2;
    for(int i = 0; i < 4; ++i) e_vsize *= node_sites[i];

    f_tmp = (Vector *)smalloc(e_vsize*sizeof(Float));

    xv[0] = node_sites[3]/2;
    xv[1] = (node_sites[3]*node_sites[0])/2;
    xv[2] = (node_sites[3]*node_sites[0]*node_sites[1])/2;


    // Initialize boundary condition on this node
    
    for(int i=0; i<4; i++) bc[i] = 0;

    if(GJP.Xbc() == BND_CND_APRD) bc[0] =
			      GJP.XnodeCoor() == (GJP.Xnodes()-1) ? 1 : 0 ;

    if(GJP.Ybc() == BND_CND_APRD) bc[1] =
			      GJP.YnodeCoor() == (GJP.Ynodes()-1) ? 1 : 0;

    if(GJP.Zbc() == BND_CND_APRD) bc[2] =
			      GJP.ZnodeCoor() == (GJP.Znodes()-1) ? 1 : 0;

    if(GJP.Tbc() == BND_CND_APRD) bc[3] =
			      GJP.TnodeCoor() == (GJP.Tnodes()-1) ? 1 : 0;


    // set up CBUF

    setCbufCntrlReg(1, CBUF_MODE1);
    setCbufCntrlReg(2, CBUF_MODE2);
    setCbufCntrlReg(3, CBUF_MODE3);
    setCbufCntrlReg(4, CBUF_MODE4);
  
}

int FstagTypes::FsiteOffsetChkb(const int *x) const
        { return (x[3]>>1)+xv[0]*x[0]+xv[1]*x[1]+xv[2]*x[2] ; }

//------------------------------------------------------------------
// Returns the number of exact flavors of the matrix that
// is inverted during a molecular dynamics evolution.
//------------------------------------------------------------------
int FstagTypes::ExactFlavors() const
{
  return 4;
}

//------------------------------------------------------------------
// Returns the number of spin components.
//------------------------------------------------------------------
int FstagTypes::SpinComponents() const
{
  return 1;
}

//------------------------------------------------------------------
// Returns the number of fermion field components 
// (including real/imaginary) on a site of the 4-D lattice.
//------------------------------------------------------------------
int FstagTypes::FsiteSize() const
{
  return 2 * Colors() * SpinComponents();  
  // re/im * colors * spin_components
}

//------------------------------------------------------------------
// returns 1 => The fermion fields in the evolution
//      or the CG that inverts the evolution matrix
//      are defined on a single checkerboard (half the 
//      lattice).
//------------------------------------------------------------------
int FstagTypes::FchkbEvl() const
{
  return 1;
}

//------------------------------------------------------------------
// The fermion Hamiltonian of the node sublattice.
// chi must be the solution of Cg with source phi.	       
//------------------------------------------------------------------
Float FstagTypes::FhamiltonNode( Vector *phi,  Vector *chi) {
  char *fname = "FhamiltonNode(V*,V*)";
  VRB.Func(cname,fname);

  return dotProduct((IFloat *)phi, (IFloat *)chi, e_vsize);
}

//------------------------------------------------------------------
// Sets the offsets for the fermion fields on a 
// checkerboard. The fermion field storage order
// is the canonical one. X[I] is the
// ith coordinate where i = {0,1,2,3} = {x,y,z,t}.
//------------------------------------------------------------------
int FstagTypes::FsiteOffset(const int *x) const{

    return x[0]+node_sites[0]*(x[1]+node_sites[1]*(x[2]+node_sites[3]*x[3]));
}


//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
FstagTypes::~FstagTypes()
{
  char *fname = "~FstagTypes()";
  VRB.Func(cname,fname);
  sfree(f_tmp);
}

CPS_END_NAMESPACE
