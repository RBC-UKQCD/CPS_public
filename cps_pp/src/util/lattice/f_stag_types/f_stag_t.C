#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief Definition of FstagTypes methods.

*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /space/cvs/cps/cps++/src/util/lattice/f_stag_types/f_stag_t.C,v $
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

const unsigned FstagTypes::CBUF_MODE1;
const unsigned FstagTypes::CBUF_MODE2;
const unsigned FstagTypes::CBUF_MODE3;
const unsigned FstagTypes::CBUF_MODE4;

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

    f_tmp = (Vector *)smalloc(e_vsize*sizeof(Float), cname, fname, "f_tmp");

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

    /*!
      When a field is stored in an odd-even (checkerboard) STAG order,
      this method converts a site's
      cartesian coordinates into its lattice site index. Note that this
      works for a field defined on a single parity only.
      \param x The cartesian lattice site coordinates.
      \return The single parity lattice site index.
    */
int FstagTypes::FsiteOffsetChkb(const int *x) const{
    return (x[3]>>1)+xv[0]*x[0]+xv[1]*x[1]+xv[2]*x[2] ;
}

    /*!
      When a field is stored in an odd-even (checkerboard) STAG order,
      this method converts a site's
      cartesian coordinates into its lattice site index. Note that this
      works for a field defined over the entire lattice.
      \param x The cartesian lattice site coordinates.
      \return The lattice site index.
    */
int FstagTypes::FsiteOffsetChkb_all(const int *x) const{

    return (x[3]+node_sites[3]*(x[0]+node_sites[0]*(x[1]+node_sites[1]*x[2]))
	+GJP.VolNodeSites()*((x[0]+x[1]+x[2]+x[3])%2))/2;

}


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

    return x[0]+node_sites[0]*(x[1]+node_sites[1]*(x[2]+node_sites[2]*x[3]));
}


/*!< Renormalise the smallest shift into the mass parameters - reduces
  linear algebra in the multi-mass solver.  This optimisation only
  works for staggered type fermions.*/
void FstagTypes::massRenormalise(Float *mass, Float *trueMass, 
				 int degree, Float *shift, 
				 MassRenormaliseDir direction) 
{

  if (direction == RENORM_FORWARDS) {
    *trueMass = *mass;
    *mass = sqrt((*trueMass)*(*trueMass) + shift[0]/4.0);
    Float zeroPole = shift[0];
    for (int j=0; j<degree; j++) shift[j] -= zeroPole;
  } else if (direction == RENORM_BACKWARDS) {
    Float zeroPole = 4.0*((*mass)*(*mass) - (*trueMass)*(*trueMass));
    for (int j=0; j<degree; j++) shift[j] += zeroPole;
    *mass = *trueMass;
  }

}

//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
FstagTypes::~FstagTypes()
{
  char *fname = "~FstagTypes()";
  VRB.Func(cname,fname);
  sfree(f_tmp, cname, fname, "f_tmp");
}

CPS_END_NAMESPACE
