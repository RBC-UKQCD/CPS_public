#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief Definition of FstagTypes methods.

  $Id: f_stag_t.C,v 1.3 2003-10-23 13:38:59 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2003-10-23 13:38:59 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/f_stag_types/f_stag_t.C,v 1.3 2003-10-23 13:38:59 zs Exp $
//  $Id: f_stag_t.C,v 1.3 2003-10-23 13:38:59 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.2.6.3  2003/09/23 21:42:33  zs
//  Build of the day
//
//  Revision 1.2.6.2  2003/08/29 15:56:31  zs
//  First draft of the asqtad fermion force stuff - it compiles at least!
//
//  Revision 1.2.6.1  2003/08/20 14:39:43  zs
//  Moved some methods from Fstag to FtagTypes
//
//  Revision 1.2  2003/07/24 16:53:54  zs
//  Addition of documentation via doxygen:
//  doxygen-parsable comment blocks added to many source files;
//  New target in makefile and consequent alterations to configure.in;
//  New directories and files under the doc directory.
//
//  Revision 1.2  2001/06/19 18:13:23  anj
//  Serious ANSIfication.  Plus, degenerate double64.h files removed.
//  Next version will contain the new nga/include/double64.h.  Also,
//  Makefile.gnutests has been modified to work properly, propagating the
//  choice of C++ compiler and flags all the way down the directory tree.
//  The mpi_scu code has been added under phys/nga, and partially
//  plumbed in.
//
//  Everything has newer dates, due to the way in which this first alteration was handled.
//
//  Anj.
//
//  Revision 1.2  2001/05/25 06:16:09  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: f_stag_t.C,v $
//  $Revision: 1.3 $
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
#include <util/stag.h>
CPS_START_NAMESPACE

//------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------
FstagTypes::FstagTypes(){

    cname = "FstagTypes";
    char *fname = "FstagTypes()";
    VRB.Func(cname,fname);

    e_vsize = VECT_LEN/2;
    for(int i = 0; i < 4; ++i) e_vsize *= node_sites[i];

    f_tmp = (Vector *)smalloc(e_vsize*sizeof(Float));

    xv[0] = node_sites[3]/2;
    xv[1] = (node_sites[3]*node_sites[0])/2;
    xv[2] = (node_sites[3]*node_sites[0]*node_sites[1])/2;


    //----------------------------------------------------------------
    // Initialize boundary condition on this node
    //----------------------------------------------------------------
    for(int i=0; i<4; i++) bc[i] = 0;

    if(GJP.Xbc() == BND_CND_APRD) bc[0] = GJP.XnodeCoor()
				      == (GJP.Xnodes()-1) ? 1 : 0 ;

    if(GJP.Ybc() == BND_CND_APRD) bc[1] = GJP.YnodeCoor()
				      == (GJP.Ynodes()-1) ? 1 : 0;

    if(GJP.Zbc() == BND_CND_APRD) bc[2] = GJP.ZnodeCoor()
				      == (GJP.Znodes()-1) ? 1 : 0;

    if(GJP.Tbc() == BND_CND_APRD) bc[3] = GJP.TnodeCoor()
				      == (GJP.Tnodes()-1) ? 1 : 0;

    dirac_init(GaugeField());

    // setup CBUF

    CBUF_MODE1 = 0xcb911548;
    CBUF_MODE2 = 0xcca52112;
    CBUF_MODE3 = 0xc98c6106;
    CBUF_MODE4 = 0xcca52112;
    setCbufCntrlReg(1, CBUF_MODE1);
    setCbufCntrlReg(2, CBUF_MODE2);
    setCbufCntrlReg(3, CBUF_MODE3);
    setCbufCntrlReg(4, CBUF_MODE4);
  
}


int FstagTypes::FsiteOffsetChkb(const int *x) const
        { return (x[3]>>1)+xv[0]*x[0]+xv[1]*x[1]+xv[2]*x[2] ; }


//------------------------------------------------------------------
// int ExactFlavors() : 
// Returns the number of exact flavors of the matrix that
// is inverted during a molecular dynamics evolution.
//------------------------------------------------------------------
int FstagTypes::ExactFlavors()
{
  return 4;
}

//------------------------------------------------------------------
// int SpinComponents() : 
// Returns the number of spin components.
//------------------------------------------------------------------
int FstagTypes::SpinComponents()
{
  return 1;
}

//------------------------------------------------------------------
// int FsiteSize() : 
// Returns the number of fermion field components 
// (including real/imaginary) on a site of the 4-D lattice.
//------------------------------------------------------------------
int FstagTypes::FsiteSize()
{
  return 2 * Colors() * SpinComponents();  
  // re/im * colors * spin_components
}

//------------------------------------------------------------------
// int FchkbEvl() :
// returns 1 => The fermion fields in the evolution
//      or the CG that inverts the evolution matrix
//      are defined on a single checkerboard (half the 
//      lattice).
//------------------------------------------------------------------
int FstagTypes::FchkbEvl()
{
  return 1;
}

//------------------------------------------------------------------
// Float FhamiltonNode(Vector *phi, Vector *chi):
// The fermion Hamiltonian of the node sublattice.
// chi must be the solution of Cg with source phi.	       
//------------------------------------------------------------------
Float FstagTypes::FhamiltonNode(Vector *phi, Vector *chi){
  char *fname = "FhamiltonNode(V*,V*)";
  VRB.Func(cname,fname);

  return dotProduct((IFloat *)phi, (IFloat *)chi, e_vsize);
}

//------------------------------------------------------------------
// int FsiteOffset(const int *x):
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
FstagTypes::~FstagTypes(){
  char *fname = "~FstagTypes()";
  VRB.Func(cname,fname);

  destroy_dirac_buf();
  sfree(f_tmp);
}

CPS_END_NAMESPACE
