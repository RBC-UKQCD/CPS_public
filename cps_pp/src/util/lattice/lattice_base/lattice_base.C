#include<config.h>
#include<util/qcdio.h>
#if TARGET == QCDOC
#include<qalloc.h>
#endif
/*!\file
  \brief  Lattice class methods.
  
  $Id: lattice_base.C,v 1.73 2013-06-25 12:51:12 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2013-06-25 12:51:12 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/lattice_base/lattice_base.C,v 1.73 2013-06-25 12:51:12 chulwoo Exp $
//  $Id: lattice_base.C,v 1.73 2013-06-25 12:51:12 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: lattice_base.C,v $
//  $Revision: 1.73 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/lattice_base/lattice_base.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// lattice_base.C
//
// Lattice is the base abstract class
//
//------------------------------------------------------------------

#include <util/lattice.h>
#include <util/vector.h>
#include <util/gjp.h>
#include <util/pmalloc.h>
#include <util/verbose.h>
#include <util/error.h>
#include <util/random.h>
#include <util/ReadLatticePar.h>
#include <util/ReadU1LatticePar.h>
#include <util/checksum.h>
#include <util/data_shift.h>
//#include <comms/nga_reg.h>
#include <comms/glb.h>
#include <comms/scu.h>
#include <comms/cbuf.h>
#include<util/time_cps.h>

#ifdef _TARTAN
#include <math64.h>
#else
#include <math.h>
#endif

#ifdef PARALLEL
#include <comms/sysfunc_cps.h>
#endif

#if TARGET == BGL
#include <sys/bgl/bgl_sys_all.h>
#endif

//#include <util/omp_wrapper.h>

CPS_START_NAMESPACE

//------------------------------------------------------------------
//! A macro  defining the opposite direction.
//------------------------------------------------------------------
#define OPP_DIR(dir)  (((dir)+4)&7)
 
//------------------------------------------------------------------
// Strings
//------------------------------------------------------------------

//------------------------------------------------------------------
//! The number of floating point numbers in a 3x3 complex matrix.
enum { MATRIX_SIZE = 18 };
//------------------------------------------------------------------

//------------------------------------------------------------------
// Initialize static variables
//------------------------------------------------------------------
Matrix* Lattice::gauge_field = 0;
Float* Lattice::u1_gauge_field = 0;
int Lattice::is_allocated = 0;
int Lattice::is_initialized = 0;
int Lattice::u1_is_initialized = 0;
StrOrdType Lattice::str_ord = CANONICAL;
Matrix** Lattice::fix_gauge_ptr = 0;
FixGaugeType Lattice::fix_gauge_kind = FIX_GAUGE_NONE;
Float Lattice::fix_gauge_stpCond = 0.;

int Lattice::node_sites[5];
int Lattice::g_dir_offset[4];

int* Lattice::g_upd_cnt = 0 ;
Float Lattice::md_time = (Float) 0.0 ;


static Matrix m_tmp1 CPS_FLOAT_ALIGN;
static Matrix m_tmp2 CPS_FLOAT_ALIGN;
// DRAM temp buffer, used for scu transfer

//------------------------------------------------------------------
// static variables used only inside this file
//------------------------------------------------------------------
//  CRAM temp buffer
#ifdef _TARTAN
static Matrix *mp0 = (Matrix *)CRAM_SCRATCH_ADDR;	// ihdot
static Matrix *mp1 = mp0 + 1;
static Matrix *mp2 = mp1 + 1;
static Matrix *mp3 = mp2 + 1;
static Matrix *mp4 = mp3 + 1;
#else
static Matrix mt0;
static Matrix mt1;
static Matrix mt2;
static Matrix mt3;
static Matrix mt4;
static Matrix *mp0 = &mt0;		// ihdot
static Matrix *mp1 = &mt1;
static Matrix *mp2 = &mt2;
static Matrix *mp3 = &mt3;
static Matrix *mp4 = &mt4;
#endif 


uint64_t Lattice::ForceFlops=0;
int Lattice::scope_lock=0;
//------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------
/*!
  If needed, allocates memory for the gauge field.
  Initialises the random number generator and the gauge configuration.
 */
Lattice::Lattice()
{
  cname = "Lattice";
  const char *fname = "Lattice()";
  int array_size;  // On-node size of the gauge field array.

  VRB.Func(cname,fname);

  StartConfType start_conf_kind = GJP.StartConfKind();

  link_buffer = 0;

  //----------------------------------------------------------------
  // Check and set the scope_lock
  //----------------------------------------------------------------
  if(scope_lock != 0){
    ERR.General(cname,fname,
		"Only one lattice object is allowed to be on scope\n");
  }
  scope_lock = 1;


  if(!(is_allocated)){
    // Allocate memory for the gauge field.
    //--------------------------------------------------------------
    array_size = GsiteSize() * GJP.VolNodeSites() * sizeof(Float);  

    if(start_conf_kind != START_CONF_LOAD ){
//       start_conf_kind!=START_CONF_FILE){
      gauge_field = (Matrix *) pmalloc(array_size);
//     printf("gauge_field=%p\n",gauge_field);
      if( gauge_field == 0) ERR.Pointer(cname,fname, "gauge_field");
      VRB.Pmalloc(cname, fname, "gauge_field", gauge_field, array_size);
      GJP.StartConfLoadAddr(gauge_field);
    }

    //--------------------------------------------------------------
    //   Pmalloc space for the gauge update counter and initialize 
    //   the counter to zero.
    //--------------------------------------------------------------
    char *g_upd_cnt_str = "g_upd_cnt" ;
    g_upd_cnt = (int *) pmalloc(sizeof(int));
    if (g_upd_cnt == 0)
      ERR.Pointer(cname, fname, g_upd_cnt_str);
    VRB.Pmalloc(cname, fname, g_upd_cnt_str, g_upd_cnt, sizeof(int));

    *g_upd_cnt = 0 ;	// will this respect checkpoints with pmalloc???


    // initialize node_sites[] and g_dir_offset[]
    //--------------------------------------------------------------
    node_sites[0] = GJP.XnodeSites();
    node_sites[1] = GJP.YnodeSites();
    node_sites[2] = GJP.ZnodeSites();
    node_sites[3] = GJP.TnodeSites();
    node_sites[4] = GJP.SnodeSites();
    g_dir_offset[0] = 4;
    g_dir_offset[1] = g_dir_offset[0]*GJP.XnodeSites();
    g_dir_offset[2] = g_dir_offset[1]*GJP.YnodeSites();
    g_dir_offset[3] = g_dir_offset[2]*GJP.ZnodeSites();

    // Set is_allocated flag to 1
    //--------------------------------------------------------------
    is_allocated = 1;
  }

  // Initialize the random number generator
  //--------------------------------------------------------------
  LRG.Initialize();

  // Initialize Molecular Dynamics time counter
  //----------------------------------------------------------------
  md_time = 0.0;


  //Load gauge field configuration

  // QED (angles)
  //----------------------------------------------------------------
if ( GJP.ExtInitialized()){
  StartConfType start_u1_conf_kind = GJP.StartU1ConfKind();
    if(start_u1_conf_kind != START_CONF_LOAD ){
      u1_gauge_field = (Float *) pmalloc(array_size/18);
      if( u1_gauge_field == 0) ERR.Pointer(cname,fname, "u1_gauge_field");
      VRB.Pmalloc(cname, fname, "u1_gauge_field", u1_gauge_field, array_size/9);
      GJP.StartU1ConfLoadAddr(u1_gauge_field);
    }
  if(start_u1_conf_kind == START_CONF_ORD){
    SetU1GfieldOrd();
    u1_is_initialized = 1;
    str_ord = CANONICAL;
    VRB.Flow(cname,fname, "Ordered starting U1 configuration\n");
    GJP.StartU1ConfKind(START_CONF_MEM);
  }
  else if(start_u1_conf_kind == START_CONF_FILE){
#if TARGET == NOARCH || TARGET == BGL || TARGET == BGP
    VRB.Flow(cname,fname, "Load starting U1 configuration addr = %x\n",
             gauge_field);
    ReadU1LatticeParallel rd_lat(*this,GJP.StartU1ConfFilename());
    str_ord = CANONICAL;
    u1_is_initialized = 1;
    GJP.StartU1ConfKind(START_CONF_MEM);
#else
    //???
    VRB.Flow(cname,fname, "File starting U1 configuration\n");
    ERR.NotImplemented(cname,fname,
    "Starting U1 config. type  START_CONF_FILE not implemented\n");
#endif
  }
}
  // QCD links
  //----------------------------------------------------------------
  if(start_conf_kind == START_CONF_ORD){
    SetGfieldOrd();
    is_initialized = 1;
    str_ord = CANONICAL;
    VRB.Flow(cname,fname, "Ordered starting configuration\n");
    GJP.StartConfKind(START_CONF_MEM);
  }
  else if(start_conf_kind == START_CONF_DISORD){
    SetGfieldDisOrd();
    is_initialized = 1;
    str_ord = CANONICAL;
    VRB.Flow(cname,fname, "Disordered starting configuration\n");
    GJP.StartConfKind(START_CONF_MEM);
  }
  else if(start_conf_kind == START_CONF_FILE){
//#if TARGET == QCDOC || TARGET == NOARCH || TARGET == BGL || TARGET == BGP
#if 1
//    gauge_field = GJP.StartConfLoadAddr();
    VRB.Flow(cname,fname, "Load starting configuration addr = %x\n",
	     gauge_field);
    ReadLatticeParallel rd_lat(*this,GJP.StartConfFilename());
    str_ord = CANONICAL;
    is_initialized = 1;
    GJP.StartConfKind(START_CONF_MEM);
#else
    //???
    VRB.Flow(cname,fname, "File starting configuration\n");
    ERR.NotImplemented(cname,fname, 
    "Starting config. type  START_CONF_FILE not implemented\n");
#endif
  }
  else if(start_conf_kind == START_CONF_LOAD){
    gauge_field = GJP.StartConfLoadAddr();
    is_initialized = 1;
    str_ord = CANONICAL;
    VRB.Flow(cname,fname, "Load starting configuration addr = %x\n",
	     gauge_field);
    GJP.StartConfKind(START_CONF_MEM);
  }
  else if(start_conf_kind == START_CONF_MEM){
    VRB.Flow(cname,fname, "Memory starting configuration addr = %x\n",
	     gauge_field);
    if(is_initialized == 0){
      if(gauge_field != GJP.StartConfLoadAddr())
	ERR.General(cname,fname, 
	"Starting config. load address is %x but the pmalloc address %x is not the same\n",
		    GJP.StartConfLoadAddr(), gauge_field);
      else
	is_initialized = 1;
    }
    else{
      Convert(CANONICAL);
    }
  }
  else {
    ERR.General(cname,fname, 
		"Unknown starting config. type %d\n", 
		int(start_conf_kind));
  }
#if 1

  // if GJP.Snodes() != 1 check that the gauge field is identical
  // on all s-slices.
  //----------------------------------------------------------------
  if(GJP.Snodes() != 1) {
    VRB.Flow(cname,fname, "Checking gauge field across s-slices\n");
     GsoCheck();
  }

  //----------------------------------------------------------------
  // Added here by Ping for anisotropic lattices 
  // Purpose : 
  //      rescale links in the special anisotropic direction
  //      Ut *= xi.
  // Why Reunitarize() instead of MltFloat():
  //   Since the configuration could either be unitarized
  //   (START_CONF_ORD or START_CONF_DISORD) or have absorbed the
  //   factor xi already (START_CONF_MEM), the simpliest universal way 
  //   is to reunitarize the lattice which also rescales the links
  //   in the special anisotropic direction correctly for all
  //   normalizations.
  //----------------------------------------------------------------
  if (GJP.XiBare() != 1.0)
    Reunitarize();  
#endif

  smeared = 0;
  }


//------------------------------------------------------------------
// Destructor
/*!
  Note that the destructor does not free any memory allocated by the
  constructor for the gauge field. This is a feature, not a bug.
*/
//------------------------------------------------------------------
Lattice::~Lattice()
{
#if 1
  const char *fname = "~Lattice()";
  VRB.Func(cname,fname);
  
  // if GJP.Snodes() != 1 check that the gauge field is identical
  // on all s-slices.
  //----------------------------------------------------------------
  if(GJP.Snodes() != 1) {
    VRB.Flow(cname,fname, "Checking gauge field across s-slices\n");
    GsoCheck();
  }

  //----------------------------------------------------------------
  // Release the scope_lock
  //----------------------------------------------------------------
  scope_lock = 0;
  

  // Undo the scaling of links on anisotropic lattices  -- Ping
  //----------------------------------------------------------------
  // No data manipulations if the scaling factor is 1.0
  // [built into lat.MltFloat].
  MltFloat(1.0 / GJP.XiBare(), GJP.XiDir());
//  sync();
#endif

}

//------------------------------------------------------------------
// const Matrix *GaugeField(void) const:
// Returns the pointer to the gauge field configuration.
//------------------------------------------------------------------
//Matrix *Lattice::GaugeField(void) const
//{
//  return gauge_field;
//}
//------------------------------------------------------------------
// Float *U1GaugeField(void) const:
// Returns the pointer to the u1 gauge field configuration.
//------------------------------------------------------------------
Float *Lattice::U1GaugeField(void) const
{
  return u1_gauge_field;
}

// Mult su3 links by u1 phase factor, exp(i theta_mu/L_mu)
// (i.e., twisted b.c.)
void Lattice::twist_bc(int sign)
{
  // multiply su3 links by u1 links
  // assumes both are in canonical order!
  // link =  exp(-i A_mu Q) where Q is the
  // quark electric charge or twist angle
  int x[4];
  const double PI = 3.1415926535897932384626433832795028842;
  const Float eps = 1e-16;
  Float theta[4];
  for(int mu=0; mu<4; ++mu)
    theta[mu] = sign*2.0*PI*GJP.TwistBc(mu)/GJP.Sites(mu);

  Matrix temp;

  for(x[3]=0; x[3]<GJP.TnodeSites();x[3]++){
    for(x[2]=0; x[2]<GJP.ZnodeSites();x[2]++){
      for(x[1]=0; x[1]<GJP.YnodeSites();x[1]++){
        for(x[0]=0; x[0]<GJP.XnodeSites();x[0]++){
          for(int mu=0; mu<4; ++mu) {
	    if(fabs(theta[mu]) < eps) continue;
            int offset_x = GsiteOffset(x);
            // the link matrix at current site, in direction mu:
            Matrix *su3_link = GaugeField() + offset_x + mu;
	    Rcomplex cc(cos(theta[mu]),sin(theta[mu]));
            vecTimesComplex((IFloat*)&temp, (IFloat)cc.real(), (IFloat)cc.imag(), (const IFloat*)su3_link, 18);
            moveMem((IFloat*)su3_link, (IFloat*)&temp,sizeof(IFloat)*18);
          }//loop(mu)
        }
      }
    }
  }
  // must re-smear the lattice after implementing new bc
  // so reset
  smeared =0;
}


// Mult su3 links by u1 phase factor, exp(i Q theta)
void Lattice::mult_su3_links_by_u1(const Float Q)
{
  // multiply su3 links by u1 links
  // assumes both are in canonical order!
  // link =  exp(-i A_mu Q) where Q is the
  // quark electric charge or twist angle
  int x[4];

  Matrix temp;
  for(x[3]=0; x[3]<GJP.TnodeSites();x[3]++){
    for(x[2]=0; x[2]<GJP.ZnodeSites();x[2]++){
      for(x[1]=0; x[1]<GJP.YnodeSites();x[1]++){
        for(x[0]=0; x[0]<GJP.XnodeSites();x[0]++){
          for(int mu=0; mu<4; ++mu) {
            int offset_x = GsiteOffset(x);
            // the link matrix at current site, in direction mu:
            Matrix *su3_link = GaugeField() + offset_x + mu;
            // the u1 link at current site, in direction mu:
            Float theta = *(U1GaugeField() + offset_x + mu);
            theta *= Q;
            Rcomplex cc(cos(theta),sin(theta));
            vecTimesComplex((IFloat*)&temp, (IFloat)cc.real(), (IFloat)cc.imag(), (const IFloat*)su3_link, 18);
            moveMem((IFloat*)su3_link, (IFloat*)&temp,sizeof(IFloat)*18);
          }//loop(mu)
        }
      }
    }
  }
}

// multiply by exp(i potential) of external B field
// it is linear in coordinates. Includes boundary links
// to ensure flux is same on every plaquette, assunming
// q B = 2 pi n / L_i L_j
void Lattice::mult_su3_links_by_mag_field(Float* B0)
{
  // multiply su3 links by u1 links
  // assumes both are in canonical order!
  // link =  exp(-i A_mu Q) where Q is the
  // quark electric charge or twist angle
  // \vec B = \vec const, so \vec A = 1/3 \vec B x \vec x
  // Oops! Really \vec A = 1/2 \vec B x \vec x
  int x[4];

  int xshift = GJP.XnodeCoor()*GJP.XnodeSites();
  int yshift = GJP.YnodeCoor()*GJP.YnodeSites();
  int zshift = GJP.ZnodeCoor()*GJP.ZnodeSites();
  int Lx = GJP.Xnodes()*GJP.XnodeSites();
  int Ly = GJP.Ynodes()*GJP.YnodeSites();
  int Lz = GJP.Znodes()*GJP.ZnodeSites();
  Matrix temp;

  printf("MAGNETIC FIELD IN GAUGE FIELD: %g %g %g\n",B0[0],B0[1],B0[2]);

  for(x[3]=0; x[3]<GJP.TnodeSites();x[3]++){
    for(x[2]=0; x[2]<GJP.ZnodeSites();x[2]++){
      for(x[1]=0; x[1]<GJP.YnodeSites();x[1]++){
        for(x[0]=0; x[0]<GJP.XnodeSites();x[0]++){
          // only for spatial links, A_t = 0
          for(int mu=0; mu<3; ++mu) {
            int offset_x = GsiteOffset(x);
            // the link matrix at current site, in direction mu:
            Matrix *su3_link = GaugeField() + offset_x + mu;
            // the linear vector potential, compactified
            int gx = x[0]+xshift;
            int gy = x[1]+yshift;
            int gz = x[2]+zshift;
            Float theta;
            switch (mu) {
            case 0:
              theta = (Float)gz*B0[1]-(Float)gy*B0[2];
              // fix up the boundary
              if(gx==Lx-1){
                theta -= (Float)gy*B0[2]*Lx;
                theta += (Float)gz*B0[1]*Lx;
              }
              break;
            case 1:
              theta = (Float)gx*B0[2]-(Float)gz*B0[0];
              // fix up the boundary
              if(gy==Ly-1){
                theta += (Float)gx*B0[2]*Ly;
                theta -= (Float)gz*B0[0]*Ly;
              }
              break;
            case 2:
              theta = (Float)gy*B0[0]-(Float)gx*B0[1];
              // fix up the boundary
              if(gz==Lz-1){
                theta -= (Float)gx*B0[1]*Lz;
                theta += (Float)gy*B0[0]*Lz;
              }
              break;
            default:
              ERR.General(cname, "mult_su3_links_by_const_mag_field", "Someone messed up!\n") ;
            }
            //
            theta *= 0.5;
            Rcomplex cc(cos(theta),sin(theta));
            vecTimesComplex((IFloat*)&temp, (IFloat)cc.real(), (IFloat)cc.imag(), (const IFloat*)su3_link, 18);
            moveMem((IFloat*)su3_link, (IFloat*)&temp,sizeof(IFloat)*18);
          }//loop(mu)
        }
      }
    }
  }
}


//------------------------------------------------------------------
/*! Copies the array pointed to by u into the gauge configuration.
   \param u The array to be copied from.
   \post The gauge configuration ia a copy of the array \a u
*/

// void GaugeField(Matrix *u) const:
// Copies the array pointed to by gauge_field to the
// array pointed to by u.
//------------------------------------------------------------------
void Lattice::GaugeField(Matrix *u)
{
  const char *fname = "GaugeField(M*)";
  VRB.Func(cname,fname);
  int size;

  size = GsiteSize() * GJP.VolNodeSites() * sizeof(IFloat);  

  // Copy from u to gauge_field
  //----------------------------
  moveMem((IFloat *) gauge_field, (IFloat *) u, size);
  smeared = 0;
}

void Lattice::U1GaugeField(Float *u)
{
  char *fname = "U1GaugeField(F*)";
  VRB.Func(cname,fname);
  int size;

  size = U1GsiteSize() * GJP.VolNodeSites() * sizeof(IFloat);

  // Copy from u to gauge_field
  //----------------------------
  moveMem((IFloat *) u1_gauge_field, (IFloat *) u, size);
  smeared = 0;
}

//------------------------------------------------------------------
/*!
  \param u An array
  \post The array \a u contains a copy of the gauge configuration.
 */

// void CopyGaugeField(Matrix* u) const:
// Copies the array pointed to by gauge_field to the
// array pointed to by u.
//------------------------------------------------------------------
void Lattice::CopyGaugeField(Matrix* u)
{
  const char *fname = "CopyGaugeField(M*)";
  VRB.Func(cname,fname);
  int size;

  size = GsiteSize() * GJP.VolNodeSites() * sizeof(IFloat);  

  // Copy from gauge_field to u
  //----------------------------
  moveMem((IFloat *) u, (IFloat *) gauge_field, size);
}

void Lattice::CopyU1GaugeField(Float* u)
{
  char *fname = "CopyU1GaugeField(f*)";
  VRB.Func(cname,fname);
  int size;

  size = U1GsiteSize() * GJP.VolNodeSites() * sizeof(IFloat);

  // Copy from gauge_field to u
  //----------------------------
  moveMem((IFloat *) u, (IFloat *) u1_gauge_field, size);
}

//------------------------------------------------------------------
// Compares the gauge configuration to the lattice instance;
// Returns 1 if identical, 0 if not
//------------------------------------------------------------------
int Lattice::CompareGaugeField(Matrix* u)
{
  const char *fname = "CompareGaugeField(M*)";
  VRB.Func(cname,fname);

  int m_size = GsiteSize() * GJP.VolNodeSites() ;
  Float *g1 = (Float*)u;
  Float *g2 = (Float*)gauge_field;
  int val=1;
  for (int i=0; i<m_size; i++) {
    if (*(g1+i) != *(g2+i) ){
	double diff = *(g1+i) - *(g2+i);
        fprintf(stderr,"Node %d: u[%d](%0.16e)!=gauge_field[%d](%0.16e), rel. diff=%e \n", UniqueID(),i,*(g1+i),i,*(g2+i),diff/(*(g1+1)) );
       val=0;
       exit(-3);
    }
  }

  return val;

}

//------------------------------------------------------------------
// StrOrd(): returns the value of the current storage order.
//------------------------------------------------------------------
StrOrdType Lattice::StrOrd(void) const
{
  return str_ord;
}


//------------------------------------------------------------------
// int Colors();
// Returns the number of colors.	  
//------------------------------------------------------------------
int Lattice::Colors(void) const
{
  return GJP.Colors();
}


//------------------------------------------------------------------
/*!
  \return The number of components (real numbers) of the gauge field
  at each lattice site, \a i.e. taking into account the number of complex,
  colour, and lattice direction components.
 */
//------------------------------------------------------------------
int Lattice::GsiteSize(void)
{
  return  2 * Colors() * Colors() * 4; //re/im * colors*colors * dim
}
int Lattice::U1GsiteSize(void)
{
  return  4; // dim
}


//--------------------------------------------------------------------------
// const Matrix *
// Lattice::GetLink(const int *site, int dir) const
//--------------------------------------------------------------------------
/*!
  Get a link at specified coordinates and direction, whether on node
  or off node. The coordinates are defined relative to the local lattice
  origin.
  
  \param site The lattice coordinates[x,y,z,t]
  which could be out of range, \a i.e.located off-node.
  \param dir The direction 0, 1, 2 or 3 for U_x, U_y, U_z and U_t respectively.
  \return a pointer to the link. If off-node, it points to a block of
  static memory. Be careful to use it!
*/
// GRF Notes:
//  - modified from PC's DiracOpClover::GetLink()
//  - now properly handles the case where link is two or more nodes in
//    negative direction
//  - currently, the code only works if the lattice is in the canonical
//    order.  This could be fixed by making GsiteOffset work correctly
//    for any supported storage order (not hard!).
//  - reminder that this code assumes that the number of sites per node
//    is the same on each node in a given direction.  Otherwise, this
//    won't work properly.
//--------------------------------------------------------------------------
const Matrix *
Lattice::GetLink(const int *site, int dir) const
{
const char *fname = "GetLink()";
//VRB.Func(cname,fname);

  // offset out-of-range coordinates site[] into on_node_site[]
  // in order to locate the link
  //------------------------------------------------------------------------
  int on_node_site[4];
  int on_node = 1;
  const Matrix *on_node_link;
  {
    for (int i = 0; i < 4; ++i) {
      on_node_site[i] = site[i] ;
      while (on_node_site[i] < 0) {
        on_node_site[i] += node_sites[i] ;
      }
      on_node_site[i] %= node_sites[i];
      if (on_node_site[i] != site[i]) {
        on_node = 0;
      }
    }
    on_node_link = gauge_field + GsiteOffset(on_node_site) + dir ;
  }

#ifndef PARALLEL
  return on_node_link;
#endif

  // send to the destination node if the site is off-node
  //------------------------------------------------------------------------
  if (on_node) {
    return on_node_link;
  } else {
    Matrix send = *on_node_link;
    Matrix &recv = m_tmp1 ;
    for (int i = 0; i < 4; ++i) {
      while (site[i] != on_node_site[i]) {
        if (site[i] < 0) {
          getMinusData((IFloat *)&recv, (IFloat *)&send, MATRIX_SIZE, i);
          on_node_site[i] -= node_sites[i];
        } else {
          getPlusData((IFloat *)&recv, (IFloat *)&send, MATRIX_SIZE, i);
          on_node_site[i] += node_sites[i];
        }
        send = recv;
      }
    }
    return &recv ;
  }
}


// get U_mu(x+v)

const Matrix* Lattice::
GetLinkOld(Matrix *g_offset, const int *x, int v, int mu) const
{
    if(x[v] == node_sites[v]-1) {  // off node
        getPlusData((IFloat *)&m_tmp1,
	    (IFloat *)(g_offset-x[v]*g_dir_offset[v]+mu),
	    MATRIX_SIZE, v);
	return &m_tmp1;
    } else {
        return g_offset+g_dir_offset[v]+mu;
    }
}

// U1 field
const Float* Lattice::
GetLinkOld(Float *g_offset, const int *x, int v, int mu) const
{
    if(x[v] == node_sites[v]-1) {  // off node
        getPlusData((IFloat *)&m_tmp1,
            (IFloat *)(g_offset-x[v]*g_dir_offset[v]+mu),
            1, v);
        return (Float*)&m_tmp1;
    } else {
        return g_offset+g_dir_offset[v]+mu;
    }
}


const unsigned CBUF_MODE2 = 0xcca52112;
const unsigned CBUF_MODE4 = 0xcca52112;


//------------------------------------------------------------------
/*!
  The staple sum around the link \f$ U_mu(x) \f$ is
\f[
   \sum_{\nu \neq \mu}[
              U_\nu(x+\mu) U^\dagger_\mu(x+\nu) U^\dagger_\nu(x)
           +  U^\dagger_\nu(x+\mu-\nu) U^\dagger_\mu(x-\nu) U_\nu(x-\nu)]     
\f]
  \param stap The computed staple sum.
  \param x The coordinates of the lattice site 
  \param mu The link direction 
*/
//------------------------------------------------------------------
void Lattice::Staple(Matrix& stap, int *x, int mu)
{
//const char *fname = "Staple(M&,i*,i)";
//VRB.Func(cname,fname);

  // set cbuf
  setCbufCntrlReg(2, CBUF_MODE2);
  setCbufCntrlReg(4, CBUF_MODE4);

  const Matrix *p1;
  int offset_x = GsiteOffset(x);
  Matrix *g_offset = GaugeField()+offset_x;

  for(int nu = 0; nu < 4; ++nu) {
    if(nu != mu) {

      //----------------------------------------------------------
      // mp3 = U_u(x+v)~
      //----------------------------------------------------------
      p1 = GetLinkOld(g_offset, x, nu, mu);
      mp3->Dagger((IFloat *)p1);


      //----------------------------------------------------------
      // p1 = &U_v(x+u)
      //----------------------------------------------------------
      p1 = GetLinkOld(g_offset, x, mu, nu);


      //----------------------------------------------------------
      // mp2 = U_v(x+u) U_u(x+v)~
      //----------------------------------------------------------
      mDotMEqual((IFloat *)mp2, (const IFloat *)p1,
		 (const IFloat *)mp3);
      

      //----------------------------------------------------------
      //  mp3 = U_v(x)~
      //----------------------------------------------------------
      mp3->Dagger((IFloat *)(g_offset+nu));
      

      //----------------------------------------------------------
      // calculate U_v(x+u)*U_u(x+v)~*U_v(x)~ = mp2 * mp3
      //----------------------------------------------------------
      if( nu == 0  ||  (mu==0 && nu==1) )
	mDotMEqual((IFloat *)&stap, (const IFloat *)mp2,
		   (const IFloat *)mp3);
      else
	mDotMPlus((IFloat *)&stap, (const IFloat *)mp2,
		  (const IFloat *)mp3);
      


      
      //----------------------------------------------------------
      //  calculate U_v(x+u-v)~ U_u(x-v)~ U_v(x-v)
      //----------------------------------------------------------
      int off_pv = (x[nu] == 0) ?
              (node_sites[nu]-1)*g_dir_offset[nu]
	    : -g_dir_offset[nu];

      Matrix *g_offpv = g_offset+off_pv;


      //----------------------------------------------------------
      // p1 = U_v(x+u-v)
      // mp3 = U_v(x+u-v)
      //----------------------------------------------------------
      p1 = GetLinkOld(g_offpv, x, mu, nu);
      moveMem((IFloat *)mp3, (IFloat *)p1, 
	      MATRIX_SIZE * sizeof(IFloat));


      //----------------------------------------------------------
      // mp2 = U_u(x-v) U_v(x+u-v)
      //----------------------------------------------------------
      mDotMEqual((IFloat *)mp2, (const IFloat *)(g_offpv+mu),
		 (const IFloat *)mp3);
      

      //----------------------------------------------------------
      // mp3 = U_v(x+u-v)~ U_u(x-v)~ = mp2~
      //----------------------------------------------------------
      mp3->Dagger((IFloat *)mp2);


      //----------------------------------------------------------
      // mp2 = U_v(x-v)
      //----------------------------------------------------------
      moveMem((IFloat *)mp2, (const IFloat *)(g_offpv+nu),
	  MATRIX_SIZE * sizeof(IFloat));
      

      //----------------------------------------------------------
      // stap += mp3 * mp2
      //----------------------------------------------------------
      if(x[nu] == 0) {	// x-v off node
	mDotMEqual((IFloat *)&m_tmp1, (const IFloat *)mp3,
		   (const IFloat *)mp2);
	
	// m_tmp2 = U_v(x+u-v)*U_u(x-v)*U_v(x-v)^dag
	getMinusData((IFloat *)&m_tmp2, (IFloat *)&m_tmp1,
		     MATRIX_SIZE, nu);
	
	//stap += m_tmp2;
	vecAddEquVec((IFloat *)&stap, (IFloat *)&m_tmp2,
		MATRIX_SIZE);
	
      } else {
	mDotMPlus((IFloat *)&stap, (const IFloat *)mp3,
		  (const IFloat *)mp2);
        // dummy read
	*((IFloat *)mp2) = *((IFloat *)g_offpv);
      }
    }
  }
}

// U1 plaquette
Float Lattice::ReU1Plaq(int *x, int mu, int nu) const
{
//char *fname = "Staple(C&,i*,i)";
//VRB.Func(cname,fname);

  // set cbuf
  setCbufCntrlReg(2, CBUF_MODE2);
  setCbufCntrlReg(4, CBUF_MODE4);

  int offset_x = GsiteOffset(x);
  Float *g_offset = U1GaugeField()+offset_x;

  //----------------------------------------------------------
  // U_u(x)
  //----------------------------------------------------------
  Float p1 = *((IFloat *)(g_offset+mu));

  //----------------------------------------------------------
  // U_v(x+u)
  //----------------------------------------------------------
  p1 += *(GetLinkOld(g_offset, x, mu, nu));

  //----------------------------------------------------------
  // U_u(x+v)^dagger
  //----------------------------------------------------------
  p1 -= *(GetLinkOld(g_offset, x, nu, mu));

  //----------------------------------------------------------
  // U_v(x)^dagger
  //----------------------------------------------------------
  p1 -= *((IFloat *)(g_offset+nu));

  return cos(p1);

}

//------------------------------------------------------------------
// RectStaple(Matrix& stap, int *x, int mu):
// It calculates the rectangle staple field at x, mu.
/*! The 5-link rectangle staple sum around the link \f$ U_\mu(x) \f$ is:
\f[
 \sum_{\nu \neq \mu}\left[\right.     
 U_\mu(x+\mu) U_\nu(x+2\mu) U^\dagger_\mu(x+\mu+\nu)
 U^\dagger_\mu(x+\nu)  U^\dagger_\nu(x) \f]\f[
 + U_\mu(x+\mu)    U^\dagger_\nu(x+2\mu-\nu) U^\dagger_\mu(x+\mu-\nu) 
 U^\dagger_\mu(x-\nu)  U_\nu(x-\nu) \f]\f[
 + U_\nu(x+\mu)    U^\dagger_\mu(x+\nu)  U^\dagger_\mu(x-\mu+\nu)
 U^\dagger_\nu(x-\mu) U_\mu(x-\mu)  \f]\f[
 + U^\dagger_\nu(x+\mu-\nu) U^\dagger_\mu(x-\nu)  U^\dagger_\mu(x-\mu-\nu)
 U_\nu(x-\mu-\nu) U_\mu(x-\mu) \f]\f[
 + U_\nu(x+\mu)    U_\nu(x+\mu+\nu)   U^\dagger_\mu(x+2\nu)
 U^\dagger_\nu(x+\nu)  U^\dagger_\nu(x) \f]\f[
 + U^\dagger_\nu(x+\mu-\nu) U^\dagger_\nu(x+\mu-2\nu) U^\dagger_\mu(x-2\nu)
 U_\nu(x-2\nu)  U_\nu(x-\nu)       
\left.\right]
\f]

  \param x The coordinates of the lattice site 
  \param mu The link direction
  \param rect The computed staple sum.
*/
//------------------------------------------------------------------
void Lattice::RectStaple(Matrix& rect, int *x, int mu)
{
//const char *fname = "RectStaple(M&,i*,i)" ;
//VRB.Func(cname, fname) ;
//VRB.Debug(cname,fname, "rect %3i %3i %3i %3i ; %i\n",
//          x[0], x[1], x[2], x[3], mu) ;

  int link_site[4] ;

  // set CBUF
  setCbufCntrlReg(4, CBUF_MODE4) ;

  //----------------------------------------------------------------------------
  // do a dummy read from the DRAM image controlled by CBUF mode ctrl reg 0
  // to guarantee CBUF will start a new proces  on the next read from
  // a DRAM image controlled by CBUF mode ctrl reg 4.
  //----------------------------------------------------------------------------
#ifdef _TARTAN
  *((unsigned *)mp2) = *((unsigned *)0x2000) ;
#endif

  const Matrix *p1;

  rect.ZeroMatrix() ;

  for (int i=0; i<4; ++i) link_site[i] = x[i] ;

  for(int nu = 0; nu < 4; ++nu)
  if(nu != mu) {

    //----------------------------------------------------------
    // mp4 = U_v(x)~
    //----------------------------------------------------------
    mp4->Dagger((IFloat *)GetLink(link_site, nu)) ;

    //----------------------------------------------------------
    // mp3 = U_mu(x+v)~
    //----------------------------------------------------------
    ++(link_site[nu]) ;
    mp3->Dagger((IFloat *)GetLink(link_site, mu)) ;

    //----------------------------------------------------------
    // mp2 = U_u(x+v)~ U_v(x)~
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp2, (const IFloat *)mp3, (const IFloat *)mp4) ;

    //----------------------------------------------------------
    // mp4 = U_u(x+u+v)~
    //----------------------------------------------------------
    ++(link_site[mu]) ;
    mp4->Dagger((IFloat *)GetLink(link_site, mu));

    //----------------------------------------------------------
    // mp3 = U_u(x+u+v)~ U_u(x+v)~ U_v(x)~
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp3, (const IFloat *)mp4, (const IFloat *)mp2);

    //----------------------------------------------------------
    // p1 = &U_v(x+2u)
    //----------------------------------------------------------
    ++(link_site[mu]) ;
    --(link_site[nu]) ;
    p1 = GetLink(link_site, nu) ;

    //----------------------------------------------------------
    // mp4 = U_v(x+2u) U_u(x+u+v)~ U_u(x+v)~ U_v(x)~
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp4, (const IFloat *)p1,
               (const IFloat *)mp3);

    //----------------------------------------------------------
    // mp2 = U_u(x+u)
    //----------------------------------------------------------
    --(link_site[mu]) ;
    moveMem((IFloat *)mp2, (const IFloat *)GetLink(link_site, mu),
            MATRIX_SIZE * sizeof(IFloat)) ;

    //----------------------------------------------------------
    // rect += U_u(x+u) U_v(x+2u) U_u(x+u+v)~ U_u(x+v)~ U_v(x)~
    //----------------------------------------------------------
    mDotMPlus((IFloat *)&rect, (const IFloat *)mp2, (const IFloat *)mp4);



    //----------------------------------------------------------
    // mp4 = U_v(x+2u-v)~
    //----------------------------------------------------------
    ++(link_site[mu]) ;
    --(link_site[nu]) ;
    mp4->Dagger((IFloat *)GetLink(link_site, nu)) ;

    //----------------------------------------------------------
    // mp3 = U_u(x+u) U_v(x+2u-v)~
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp3, (const IFloat *)mp2, (const IFloat *)mp4) ;

    //----------------------------------------------------------
    // mp4 = U_u(x+u-v)~
    //----------------------------------------------------------
    --(link_site[mu]) ;
    mp4->Dagger((IFloat *)GetLink(link_site, mu)) ;

    //----------------------------------------------------------
    // mp2 = U_u(x+u) U_v(x+2u-v)~ U_u(x+u-v)~
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp2, (const IFloat *)mp3, (const IFloat *)mp4) ;

    //----------------------------------------------------------
    // mp4 = U_u(x-v)~
    //----------------------------------------------------------
    --(link_site[mu]) ;
    mp4->Dagger((IFloat *)GetLink(link_site, mu)) ;

    //----------------------------------------------------------
    // mp3 = U_u(x+u) U_v(x+2u-v)~ U_u(x+u-v)~ U_u(x-v)~
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp3, (const IFloat *)mp2, (const IFloat *)mp4) ;

    //----------------------------------------------------------
    // mp2 = U_v(x-v)
    //----------------------------------------------------------
    moveMem((IFloat *)mp2, (const IFloat *)GetLink(link_site, nu),
            MATRIX_SIZE * sizeof(IFloat)) ;

    //----------------------------------------------------------
    // rect += U_u(x+u) U_v(x+2u-v)~ U_u(x+u-v)~ U_u(x-v)~ U_v(x-v)
    //----------------------------------------------------------
    mDotMPlus((IFloat *)&rect, (const IFloat *)mp3, (const IFloat *)mp2) ;



    //----------------------------------------------------------
    // p1 = &U_v(x-2v)
    //----------------------------------------------------------
    --(link_site[nu]) ;
    p1 = GetLink(link_site, nu) ;

    //----------------------------------------------------------
    // mp3 = U_v(x-2v) U_v(x-v)
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp3, (const IFloat *)p1,
               (const IFloat *)mp2) ;

    //----------------------------------------------------------
    // mp4 = U_u(x-2v)~
    //----------------------------------------------------------
    mp4->Dagger((IFloat *)GetLink(link_site, mu)) ;

    //----------------------------------------------------------
    // mp2 = U_u(x-2v)~ U_v(x-2v) U_v(x-v)
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp2, (const IFloat *)mp4, (const IFloat *)mp3) ;

    //----------------------------------------------------------
    // mp4 = U_v(x+u-2v)~
    //----------------------------------------------------------
    ++(link_site[mu]) ;
    mp4->Dagger((IFloat *)GetLink(link_site, nu)) ;

    //----------------------------------------------------------
    // mp3 = U_v(x+u-2v)~ U_u(x-2v)~ U_v(x-2v) U_v(x-v)
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp3, (const IFloat *)mp4, (const IFloat *)mp2) ;

    //----------------------------------------------------------
    // mp2 = U_v(x+u-v)~
    //----------------------------------------------------------
    ++(link_site[nu]) ;
    mp2->Dagger((IFloat *)GetLink(link_site, nu)) ;

    //----------------------------------------------------------
    // rect += U_v(x+u-v)~ U_v(x+u-2v)~ U_u(x-2v)~ U_v(x-2v) U_v(x-v)
    //----------------------------------------------------------
    mDotMPlus((IFloat *)&rect, (const IFloat *)mp2, (const IFloat *)mp3) ;



    //----------------------------------------------------------
    // mp4 = U_u(x-v)~
    //----------------------------------------------------------
    --(link_site[mu]) ;
    mp4->Dagger((IFloat *)GetLink(link_site, mu)) ;

    //----------------------------------------------------------
    // mp3 =  U_v(x+u-v)~ U_u(x-v)~
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp3, (const IFloat *)mp2, (const IFloat *)mp4) ;

    //----------------------------------------------------------
    // mp4 = U_u(x-u-v)~
    //----------------------------------------------------------
    --(link_site[mu]) ;
    mp4->Dagger((IFloat *)GetLink(link_site, mu)) ;

    //----------------------------------------------------------
    // mp2 = U_v(x+u-v)~ U_u(x-v)~ U_u(x-u-v)~
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp2, (const IFloat *)mp3, (const IFloat *)mp4) ;

//  GRF: this code may cause a hang
//
//  //----------------------------------------------------------
//  // p1 = &U_v(x-u-v)
//  //----------------------------------------------------------
//  p1 = GetLink(link_site, nu) ;
//
//  //----------------------------------------------------------
//  // mp3 = U_v(x+u-v)~ U_u(x-v)~ U_u(x-u-v)~ U_v(x-u-v)
//  //----------------------------------------------------------
//  mDotMEqual((IFloat *)mp3, (const IFloat *)mp2,
//             (const IFloat *)p1) ;
//
//  GRF: here is a work-around

    //----------------------------------------------------------
    // mp4 = U_v(x-u-v)
    //----------------------------------------------------------
    moveMem((IFloat *)mp4,
            (const IFloat *)GetLink(link_site, nu),
            MATRIX_SIZE*sizeof(IFloat)) ;

    //----------------------------------------------------------
    // mp3 = U_v(x+u-v)~ U_u(x-v)~ U_u(x-u-v)~ U_v(x-u-v)
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp3, (const IFloat *)mp2, (const IFloat *)mp4) ;

//  GRF: end work-around

    //----------------------------------------------------------
    // mp2 = U_u(x-u)
    //----------------------------------------------------------
    ++(link_site[nu]) ;
    moveMem((IFloat *)mp2, (const IFloat *)GetLink(link_site, mu),
            MATRIX_SIZE * sizeof(IFloat)) ;

    //----------------------------------------------------------
    // rect += U_v(x+u-v)~ U_u(x-v)~ U_u(x-u-v)~ U_v(x-u-v) U_u(x-u)
    //----------------------------------------------------------
    mDotMPlus((IFloat *)&rect, (const IFloat *)mp3, (const IFloat *)mp2) ;



    //----------------------------------------------------------
    // mp4 = U_v(x-u)~
    //----------------------------------------------------------
    mp4->Dagger((IFloat *)GetLink(link_site, nu)) ;

    //----------------------------------------------------------
    // mp3 = U_v(x-u)~ U_u(x-u)
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp3, (const IFloat *)mp4, (const IFloat *)mp2) ;

    //----------------------------------------------------------
    // mp4 = U_u(x-u+v)~
    //----------------------------------------------------------
    ++(link_site[nu]) ;
    mp4->Dagger((IFloat *)GetLink(link_site, mu)) ;

    //----------------------------------------------------------
    // mp2 = U_u(x-u+v)~ U_v(x-u)~ U_u(x-u)
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp2, (const IFloat *)mp4, (const IFloat *)mp3) ;

    //----------------------------------------------------------
    // mp4 = U_u(x+v)~
    //----------------------------------------------------------
    ++(link_site[mu]) ;
    mp4->Dagger((IFloat *)GetLink(link_site, mu)) ;

    //----------------------------------------------------------
    // mp3 = U_u(x+v)~ U_u(x-u+v)~ U_v(x-u)~ U_u(x-u)
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp3, (const IFloat *)mp4, (const IFloat *)mp2) ;

    //----------------------------------------------------------
    // mp2 = U_v(x+u)
    //----------------------------------------------------------
    ++(link_site[mu]) ;
    --(link_site[nu]) ;
    moveMem((IFloat *)mp2, (const IFloat *)GetLink(link_site, nu),
            MATRIX_SIZE * sizeof(IFloat)) ;

    //----------------------------------------------------------
    // rect += U_v(x+u) U_u(x+v)~ U_u(x-u+v)~ U_v(x-u)~ U_u(x-u)
    //----------------------------------------------------------
    mDotMPlus((IFloat *)&rect, (const IFloat *)mp2, (const IFloat *)mp3) ;



//  GRF: this code may cause a hang
//
//  //----------------------------------------------------------
//  // p1 = &U_v(x+u+v)
//  //----------------------------------------------------------
//  ++(link_site[nu]) ;
//  p1 = GetLink(link_site, nu) ;
//
//  //----------------------------------------------------------
//  // mp3 = U_v(x+u) U_v(x+u+v)
//  //----------------------------------------------------------
//  mDotMEqual((IFloat *)mp3, (const IFloat *)mp2,
//             (const IFloat *)p1) ;
//
//  GRF: here is a work-around

    //----------------------------------------------------------
    // mp4 = U_v(x+u+v)
    //----------------------------------------------------------
    ++(link_site[nu]) ;
    moveMem((IFloat *)mp4,
            (const IFloat *)GetLink(link_site, nu),
            MATRIX_SIZE*sizeof(IFloat)) ;

    //----------------------------------------------------------
    // mp3 = U_v(x+u) U_v(x+u+v)
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp3, (const IFloat *)mp2, (const IFloat *)mp4) ;

//  GRF: end work-around

    //----------------------------------------------------------
    // mp4 = U_u(x+2v)~
    //----------------------------------------------------------
    --(link_site[mu]) ;
    ++(link_site[nu]) ;
    mp4->Dagger((IFloat *)GetLink(link_site, mu)) ;

    //----------------------------------------------------------
    // mp2 = U_v(x+u) U_v(x+u+v) U_u(x+2v)~
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp2, (const IFloat *)mp3, (const IFloat *)mp4) ;

    //----------------------------------------------------------
    // mp4 = U_v(x+v)~
    //----------------------------------------------------------
    --(link_site[nu]) ;
    mp4->Dagger((IFloat *)GetLink(link_site, nu)) ;

    //----------------------------------------------------------
    // mp3 = U_v(x+u) U_v(x+u+v) U_u(x+2v)~ U_v(x+v)~
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp3, (const IFloat *)mp2, (const IFloat *)mp4) ;

    //----------------------------------------------------------
    // mp2 = U_v(x)~
    //----------------------------------------------------------
    --(link_site[nu]) ;
    mp2->Dagger((IFloat *)GetLink(link_site, nu)) ;

    //----------------------------------------------------------
    // rect += U_v(x+u) U_v(x+u+v) U_u(x+2v)~ U_v(x+v)~  U_v(x)~
    //----------------------------------------------------------
    mDotMPlus((IFloat *)&rect, (const IFloat *)mp3, (const IFloat *)mp2) ;

    //----------------------------------------------------------
    // dummy read to switch CBUF banks for looping
    //----------------------------------------------------------
    *((IFloat *)mp4) = *((IFloat *)p1) ;
  }
}


//------------------------------------------------------------------
// void Lattice::Plaq(Matrix &plaq, int *x, int mu, int nu) const 
//------------------------------------------------------------------
// Purpose: 
//    calculates the plaquette U_u(x) U_v(x+u) U_u(x+v)~ U_v(x)~
// Added by Ping to help code debugging, may be more useful later on.
//------------------------------------------------------------------
/*!
  The plaquette is
\f[
  U_\mu(x) U_\nu(x+\nu) U^\dagger_\mu(x+\nu) U^\dagger_\nu(x)
\f]

  \param plaq The computed plaquette.
  \param x the coordinates of the lattice site at the start of the plaquette
  \param mu The first plaquette direction.
  \param nu The second plaquette direction; should be different from \a mu.
*/
void Lattice::Plaq(Matrix &plaq, int *x, int mu, int nu) const 
{
  // set cbuf
  setCbufCntrlReg(2, CBUF_MODE2);
  setCbufCntrlReg(4, CBUF_MODE4);
 
  const Matrix *p1;
 
  //----------------------------------------
  //  "g_offset" points to the links
  //  at site "x"
  //----------------------------------------
  Matrix *g_offset = GaugeField()+GsiteOffset(x);
 
  //----------------------------------------
  //  mp3 = U_u(x) U_v(x+u)
  //    p1 = &U_v(x+u) --> mp2
  //----------------------------------------
  p1 = GetLinkOld(g_offset, x, mu, nu);
  moveMem((IFloat *)mp2, (const IFloat *)p1,
          MATRIX_SIZE * sizeof(IFloat));
  mDotMEqual((IFloat *)mp3, (const IFloat *)(g_offset+mu),
             (const IFloat *)mp2);
  
  //----------------------------------------
  //  mp1 = (U_v(x) U_u(x+v))~
  //    p1 = &U_u(x+v) --> mp1
  //    mp2 = U_v(x) U_u(x+v)
  //    mp1 = mp2~
  //----------------------------------------
  p1 = GetLinkOld(g_offset, x, nu, mu);
  moveMem((IFloat *)mp1, (const IFloat *)p1,
          MATRIX_SIZE * sizeof(IFloat));
  mDotMEqual((IFloat *)mp2, (const IFloat *)(g_offset+nu),
             (const IFloat *)mp1);
  mp1->Dagger((IFloat *)mp2);
 
  mDotMEqual((IFloat *)(&plaq), (const IFloat *)mp3, (const IFloat *)mp1);
}  
 

//------------------------------------------------------------------
/*!
  The plaquette is
\f[
  U_\mu(x) U_\nu(x+\nu) U^\dagger_\mu(x+\nu) U^\dagger_\nu(x)
\f]

  \param x the coordinates of the lattice site at the start of the plaquette
  \param mu The first plaquette direction
  \param nu The second plaquette direction; should be different from \a mu.
  \return  The real part of the trace of the plaquette.
*/
Float Lattice::ReTrPlaq(int *x, int mu, int nu) const
{
//  const char *fname = "ReTrPlaq(i*,i,i) const";
//  VRB.Func(cname,fname);

  // set cbuf
  setCbufCntrlReg(2, CBUF_MODE2);
  setCbufCntrlReg(4, CBUF_MODE4);

  const Matrix *p1;

  //----------------------------------------
  //  "g_offset" points to the links
  //  at site "x"
  //----------------------------------------
  Matrix *g_offset = GaugeField()+GsiteOffset(x);


  //----------------------------------------
  //  mp3 = U_u(x) U_v(x+u)
  //	p1 = &U_v(x+u) --> mp2
  //----------------------------------------
  p1 = GetLinkOld(g_offset, x, mu, nu);
  moveMem((IFloat *)mp2, (const IFloat *)p1,
	  MATRIX_SIZE * sizeof(IFloat));
  mDotMEqual((IFloat *)mp3, (const IFloat *)(g_offset+mu),
	     (const IFloat *)mp2);


  //----------------------------------------
  //  mp1 = (U_v(x) U_u(x+v))~
  //	p1 = &U_u(x+v) --> mp1
  //	mp2 = U_v(x) U_u(x+v)
  //	mp1 = mp2~
  //----------------------------------------
  p1 = GetLinkOld(g_offset, x, nu, mu);
  moveMem((IFloat *)mp1, (const IFloat *)p1,
	  MATRIX_SIZE * sizeof(IFloat));
  mDotMEqual((IFloat *)mp2, (const IFloat *)(g_offset+nu),
	     (const IFloat *)mp1);
  mp1->Dagger((IFloat *)mp2);


  
  mDotMEqual((IFloat *)mp2, (const IFloat *)mp3, (const IFloat *)mp1);
  return mp2->ReTr();
}


//------------------------------------------------------------------
/*!
  At a site \a x and in the \f$ \mu-\nu \f$plane , the plaquette is

\f[
  U_\mu(x) U_\nu(x+\nu) U^\dagger_\mu(x+\nu) U^\dagger_\nu(x)
\f]
  
  The sum is over all local lattice sites and all six \f$ \mu-\nu \f$ planes.

  \return The summed real trace of the plaquette.
*/
//------------------------------------------------------------------
Float Lattice::SumReTrPlaqNode(void) const
{
  const char *fname = "SumReTrPlaqNode() const";
  sync();
  VRB.Func(cname,fname);
  
  Float sum = 0.0;
  int x[4];
  
  for(x[0] = 0; x[0] < node_sites[0]; ++x[0]) {
    for(x[1] = 0; x[1] < node_sites[1]; ++x[1]) {
      for(x[2] = 0; x[2] < node_sites[2]; ++x[2]) {
	for(x[3] = 0; x[3] < node_sites[3]; ++x[3]) {
	  
	  for (int mu = 0; mu < 3; ++mu) {
	    for(int nu = mu+1; nu < 4; ++nu) {
//		VRB.Flow(cname,fname,"%d %d %d %d %d %d\n",x[0],x[1],x[2],x[3],mu,nu);
	      sum += ReTrPlaq(x,mu,nu);
	    }
	  }
	}
      }
    }
  }
//  sync();
  VRB.FuncEnd(cname,fname);
  return sum;
}

Float Lattice::SumReU1PlaqNode() const
{
  char *fname = "SumReU1PlaqNode() const";
  sync();
  VRB.Func(cname,fname);

  Float sum = 0.0;
  int x[4];

  for(x[0] = 0; x[0] < node_sites[0]; ++x[0]) {
    for(x[1] = 0; x[1] < node_sites[1]; ++x[1]) {
      for(x[2] = 0; x[2] < node_sites[2]; ++x[2]) {
        for(x[3] = 0; x[3] < node_sites[3]; ++x[3]) {

          for (int mu = 0; mu < 3; ++mu) {
            for(int nu = mu+1; nu < 4; ++nu) {
//              VRB.Flow(cname,fname,"%d %d %d %d %d %d\n",x[0],x[1],x[2],x[3],mu,nu);
              sum += ReU1Plaq(x,mu,nu);
            }
          }
        }
      }
    }
  }
  sync();
  VRB.FuncEnd(cname,fname);
  return sum;
}

//------------------------------------------------------------------
/*!
  At a site \a x and in the \f$ \mu-\nu \f$ plane, the plaquette is
\f[  
  U_\mu(x) U_\nu(x+\nu) U^\dagger_\mu(x+\nu) U^\dagger_\nu(x)
\f]
  
  The sum is over all lattice sites and all six \f$ \mu-\nu \f$planes.

  \return The globally summed real trace of the plaquette.
*/
//------------------------------------------------------------------
Float Lattice::SumReTrPlaq(void) const
{
  const char *fname = "SumReTrPlaq() const";
  VRB.Func(cname,fname);

  Float sum = SumReTrPlaqNode();
  glb_sum(&sum);
//  printf("sum= %0.18e\n",sum);
  VRB.FuncEnd(cname,fname);
  return sum;
}

Float Lattice::SumReU1Plaq(void) const
{
  char *fname = "SumReTrPlaq() const";
  VRB.Func(cname,fname);

  Float sum = SumReU1PlaqNode();
  glb_sum(&sum);
  return sum;
}

//-----------------------------------------------------------------------------
/*!
  The rectangle at site \a x in the \f$ \mu-\nu \f$ plane with the long axis
  of the rectangle in the \f$ \mu \f$ direction is:
\f[
  U_\mu(x) U_\mu(x+\mu) U_\nu(x+2\mu) U^\dagger_\mu(x+\mu+\nu) U^\dagger_\mu(x+\nu) U^\dagger_\nu(x)
\f]

  \param x the coordinates of the lattice site at the start of the rectangle
  \param mu The first rectangle direction.
  \param nu The second rectangle direction; should be different from \a mu.
  \return The computed rectangle
*/

//
//-----------------------------------------------------------------------------
Float Lattice::ReTrRect(int *x, int mu, int nu) const
{
    const char *fname = "ReTrRect(i*,i,i) const";
//  VRB.Func(cname,fname);

  int link_site[4] ;

  // set CBUF
  setCbufCntrlReg(4, CBUF_MODE4);

  //---------------------------------------------------------------------------
  // do a dummy read from the DRAM image controlled by CBUF mode ctrl reg 0
  // to guarantee CBUF will start a new proces  on the next read from
  // a DRAM image controlled by CBUF mode ctrl reg 4.
  //---------------------------------------------------------------------------
#ifdef _TARTAN
  *((unsigned *)mp2) = *((unsigned *)0x2000) ;
#endif

  const Matrix *p1;

  for (int i=0; i<4; ++i) link_site[i] = x[i] ;

  if (nu == mu) {
    ERR.General(cname, fname, "(mu == nu) not allowed.\n") ;
  } else {

    //----------------------------------------------------------
    // mp4 = U_v(x)~
    //----------------------------------------------------------
    mp4->Dagger((IFloat *)GetLink(link_site, nu)) ;

    //----------------------------------------------------------
    // mp3 = U_u(x+v)~
    //----------------------------------------------------------
    ++(link_site[nu]) ;
    mp3->Dagger((IFloat *)GetLink(link_site, mu)) ;

    //----------------------------------------------------------
    // mp2 = U_u(x+v)~ U_v(x)~
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp2, (const IFloat *)mp3, (const IFloat *)mp4) ;

    //----------------------------------------------------------
    // mp4 = U_u(x+u+v)~
    //----------------------------------------------------------
    ++(link_site[mu]) ;
    mp4->Dagger((IFloat *)GetLink(link_site, mu));

    //----------------------------------------------------------
    // mp3 = U_u(x+u+v)~ U_u(x+v)~ U_v(x)~
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp3, (const IFloat *)mp4, (const IFloat *)mp2);

    //----------------------------------------------------------
    // p1 = &U_v(x+2u)
    //----------------------------------------------------------
    ++(link_site[mu]) ;
    --(link_site[nu]) ;
    p1 = GetLink(link_site, nu) ;

    //----------------------------------------------------------
    // mp2 = U_v(x+2u) U_u(x+u+v)~ U_u(x+v)~ U_v(x)~
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp2, (const IFloat *)p1,
               (const IFloat *)mp3);

    //----------------------------------------------------------
    // p1 = &U_u(x+u)
    //----------------------------------------------------------
    --(link_site[mu]) ;
    p1 = GetLink(link_site, mu) ;

    //----------------------------------------------------------
    // mp3 = U_u(x+u) U_v(x+2u) U_u(x+u+v)~ U_u(x+v)~ U_v(x)~
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp3, (const IFloat *)p1,
               (const IFloat *)mp2);

    //----------------------------------------------------------
    // p1 = &U_u(x)
    //----------------------------------------------------------
    p1 = GetLink(x, mu) ;

    //----------------------------------------------------------
    // mp2 = U_u(x) U_u(x+u) U_v(x+2u) U_u(x+u+v)~ U_u(x+v)~ U_v(x)~
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp2, (const IFloat *)p1,
               (const IFloat *)mp3);

  }
  return mp2->ReTr();
}


//-----------------------------------------------------------------------------
/*!
  The rectangle at site \a x in the \f$ \mu-\nu \f$ plane with the long axis
  of the rectangle in the \a \mu direction is:
\f[
U_\mu(x) U_\mu(x+\mu) U_\nu(x+2\mu) U^\dagger_\mu(x+\mu+\nu)
U^\dagger_\mu(x+\nu) U^\dagger_\nu(x)
\f]
  
The sum is over all local lattice sites and all sixteen \f$ \mu-\nu \f$
combinations.

  \return The summed real trace of the rectangle.
*/
//-----------------------------------------------------------------------------
Float Lattice::SumReTrRectNode(void) const
{
//const char *fname = "SumReTrRectNode() const";
//VRB.Func(cname,fname);

  Float sum = 0.0 ;
  int x[4] ;

  for(x[0] = 0; x[0] < node_sites[0]; ++x[0])
  for(x[1] = 0; x[1] < node_sites[1]; ++x[1])
  for(x[2] = 0; x[2] < node_sites[2]; ++x[2])
  for(x[3] = 0; x[3] < node_sites[3]; ++x[3])
  for(int mu = 0; mu < 4; ++mu)
  for(int nu = 0; nu < 4; ++nu) {
    if (mu != nu) {
      sum += ReTrRect(x,mu,nu);
    }
  }
  return sum;
}


//-----------------------------------------------------------------------------
/*!
  The rectangle at site \a x in the \f$\mu-\nu \f$ plane with the long axis
  of the rectangle in the \f$ \mu \f$ direction is:
\f[
    U_\mu(x) U_\mu(x+\mu) U_\nu(x+2\mu)
    U^\dagger_\mu(x+\mu+\nu) U^\dagger_\mu(x+\nu) U^\dagger_\nu(x)
\f]
  
  The sum is over all lattice sites and all sixteen \f$ \mu-\nu \f$ combinations.

  \return The globally summed real trace of the rectangle.
*/
//-----------------------------------------------------------------------------
Float Lattice::SumReTrRect(void) const
{
  const char *fname = "SumReTrRect() const" ;
  VRB.Func(cname, fname) ;

  Float sum = SumReTrRectNode() ;
  glb_sum(&sum) ;
  return sum ;
}

//-------------------------------------------------------------------
/*!
  Given the starting site x, the directions of each step on the path
  and the number of steps. calculate the path ordered product of 
  all the links along the path and take the real part of the trace.
  Each direction is one of 0, 1, 2, 3, 4, 5, 6 or 7} corresponding to
  the directions X, Y, Z, T, -X, -Y, -Z and -T respectively.

  \param x The coordinates of the starting point of the path
  \param dir The list of directions.
  \param length The number of links in the path.
  \return The real part of the trace of the product along the path.
  
  \a N.B. The user is responsible for defining directions that close the loop.
*/
//--------------------------------------------------------------------
Float Lattice::ReTrLoop(const int *x, const int *dir,  int length)
{
  const char *fname = "ReTrLoop(i*,i,i)";
  VRB.Func(cname, fname) ;

  const unsigned CBUF_MODE4 = 0xcca52112;
  const unsigned CBUF_MODE2 = 0xcca52112;

  setCbufCntrlReg(2, CBUF_MODE2);
  setCbufCntrlReg(4, CBUF_MODE4);

  mp3->ZeroMatrix();

  PathOrdProdPlus(*mp3, x, dir, length);

  return mp3->ReTr() ;
}



//-----------------------------------------------------------------------------
/*!
  The cube loop is
\f[
      U_\mu(x) U_\nu(x+\mu) U_\rho(x+\mu+\nu) U^\dagger_\mu(x+\mu+\nu+\rho)
       U^\dagger_\nu(x+\nu+\rho) U^\dagger_\rho(x+\rho)
\f]

The sum runs over all positive values of \f$ \mu\f$, \f$ \nu>\mu \f$ and
\f$ \rho>\nu \f$
     The real part of the trace of this loop is summed over all local
     lattice sites \a x. 

     \return The locally summed real trace of the cube.
     \todo Check this code.
*/
//-----------------------------------------------------------------------------
Float Lattice::SumReTrCubeNode(void) 
{
  const char *fname = "SumReTrCubeNode()";
  VRB.Func(cname, fname) ;

  Float sum = 0.0 ;
  int x[4] ;
  int dir[6] ;

  for(x[0] = 0; x[0] < node_sites[0]; ++x[0])
    for(x[1] = 0; x[1] < node_sites[1]; ++x[1])
      for(x[2] = 0; x[2] < node_sites[2]; ++x[2])
	for(x[3] = 0; x[3] < node_sites[3]; ++x[3])
	  for(int mu = 0; mu < 4; mu++)
	    {
	      dir[0] = mu ;
	      dir[3] = OPP_DIR(mu) ;    
	      for(int nu = 0; nu < 4; nu++) 
                if(nu!=mu)
		{
		   dir[1] = nu ;
		   dir[4] = OPP_DIR(nu) ;
		   for(int rho = mu+1; rho < 4; rho++)
		     if(rho!=mu && rho!=nu)
		     {
		       dir[2] = rho ;
		       dir[5] = OPP_DIR(rho) ;     
		       sum += ReTrLoop(x,dir,6);
		     }
#if 1
		   dir[1] = OPP_DIR(nu) ;
		   dir[4] = nu ;
		   for(int rho = mu+1; rho < 4; rho++)
		     if(rho!=mu && rho!=nu)
		     {
		       dir[2] = rho ;
		       dir[5] = OPP_DIR(rho) ;     
		       sum += 0.33333333333333333333333*ReTrLoop(x,dir,6);
		     }
#endif
		}
	    }
//  ClearAllBufferedLink();
  return sum;
}


//-----------------------------------------------------------------------------
/*!
  The cube loop is

\f[
U_\mu(x) U_\nu(x+\mu) U_\rho(x+\mu+\nu) U^\dagger_\mu(x+\mu+\nu+\rho)
U^\dagger_\nu(x+\nu+\rho) U^\dagger_\rho(x+\rho)
\f]     
     The sum runs over all positive values of
     \f$ \mu, \nu>\mu \f$ and \f$ \rho>\nu \f$
     The real part of the trace of this loop is summed over all 
     lattice sites \a x. 

  \return The globally summed real trace of the cube.
*/
//-----------------------------------------------------------------------------
Float Lattice::SumReTrCube(void)
{
  const char *fname = "SumReTrCube()" ;
  VRB.Func(cname, fname) ;

  Float sum = SumReTrCubeNode() ;
  glb_sum(&sum) ;
  return sum ;
}


//------------------------------------------------------------------
// Added by Ping for anisotropic lattices
//------------------------------------------------------------------
enum {NUM_SPACE_PLAQ = 3, //!< Number of planes in a 3-dimensional lattice slice.
      NUM_TIME_PLAQ = 3, //!< Number of planes in a 3-dimensional lattice slice.
      NUM_COLORS = 3,    //!< Number of colours (again).
      NUM_DIM = 4        //!< Number of lattice dimensions.
};

//------------------------------------------------------------------
// Float AveReTrPlaqNodeNoXi(void) const
//------------------------------------------------------------------
// Normalization:  1 for ordered links
// Added by Ping for anisotropic lattices
/*!
  At a site \a x and in the \f$ \mu-\nu plane\f$, the plaquette is

\f[
  U_\mu(x) U_\nu(x+\nu) U^\dagger_\mu(x+\nu) U^\dagger_\nu(x)
\f]
  
  The sum is over all local lattice sites and all three \f$x \mu-\nu \f$ planes
  where neither \f$ \mu \f$nor \f$ \nu \f$ is the anisotropic direction.

  \return The real trace of the plaquette averaged over local sites,
  planes and colours.
*/

//------------------------------------------------------------------
Float Lattice::AveReTrPlaqNodeNoXi() const
{
  const char *fname = "AveReTrPlaqNodeNoXi()";
  VRB.Func(cname,fname);
  
  Float sum = 0.0;
  int x[4];
  
  for(x[0] = 0; x[0] < node_sites[0]; ++x[0]) {
    for(x[1] = 0; x[1] < node_sites[1]; ++x[1]) {
      for(x[2] = 0; x[2] < node_sites[2]; ++x[2]) {
        for(x[3] = 0; x[3] < node_sites[3]; ++x[3]) {
          for (int mu = 0; mu < NUM_DIM; ++mu) {
            for(int nu = mu+1; nu < NUM_DIM; ++nu) {
	      if (mu != GJP.XiDir() && nu != GJP.XiDir())
		sum += ReTrPlaq(x,mu,nu);
            }
          }
        }
      }
    }
  }
  
  return (sum / (GJP.VolNodeSites() * NUM_SPACE_PLAQ * NUM_COLORS));
}

//------------------------------------------------------------------
// Float AveReTrPlaqNodeXi(void) const
//------------------------------------------------------------------
// Normalization:  1 for ordered links
// Added by Ping for anisotropic lattices
/*!
  At a site \a x and in the \f$ \mu-\nu \f$ plane, the plaquette is

\f[
  U_\mu(x) U_\nu(x+\nu) U^\dagger_\mu(x+\nu) U^\dagger_\nu(x)
\f]
  
  The sum is over all local lattice sites and all three \f$ \mu-\nu \f$ planes
  where one of \f$ \mu \f$ or \f$\nu \f$ is the anisotropic direction. The bare
  anisotropy is taken into account here.

  \return The real trace of the plaquette averaged over local sites, planes and colours.
*/
//------------------------------------------------------------------
Float Lattice::AveReTrPlaqNodeXi() const
{
  const char *fname = "AveReTrPlaqNodeXi()";
  VRB.Func(cname,fname);
  
  Float sum = 0.0;
  int x[4];
  
  for(x[0] = 0; x[0] < node_sites[0]; ++x[0]) {
    for(x[1] = 0; x[1] < node_sites[1]; ++x[1]) {
      for(x[2] = 0; x[2] < node_sites[2]; ++x[2]) {
        for(x[3] = 0; x[3] < node_sites[3]; ++x[3]) {
          for (int mu = 0; mu < NUM_DIM; ++mu) {
            for(int nu = mu+1; nu < NUM_DIM; ++nu) {
	      if (mu == GJP.XiDir() || nu == GJP.XiDir())
		sum += ReTrPlaq(x,mu,nu);
            }
          }
        }
      }
    }
  }
  
  return (sum / (GJP.VolNodeSites() * NUM_TIME_PLAQ * NUM_COLORS *
		 GJP.XiBare() * GJP.XiBare()));
}

//------------------------------------------------------------------
// Float AveReTrPlaqNoXi(void) const
//------------------------------------------------------------------
// Normalization:  1 for ordered links
// Added by Ping for anisotropic lattices
/*!
  At a site \a x and in the \f$ \mu-\nu \f$ plane, the plaquette is
\f[
  U_\mu(x) U_\nu(x+\nu) U^\dagger_\mu(x+\nu) U^\dagger_\nu(x)
\f]
  
  The sum is over all lattice sites and all three \f$ \mu-\nu \f$ planes
  where neither \f$ \mu \f$ nor \f$ \nu \f$ is the anisotropic direction.

  \return The real trace of the plaquette averaged over sites, planes and colours.
*/
//------------------------------------------------------------------
Float Lattice::AveReTrPlaqNoXi() const
{
  const char *fname = "AveReTrPlaqNoXi()";
  VRB.Func(cname,fname);

  Float sum = AveReTrPlaqNodeNoXi();
  glb_sum(&sum);
  return (sum * GJP.VolNodeSites() / GJP.VolSites());
}

//------------------------------------------------------------------
// Float AveReTrPlaqXi(void) const
//------------------------------------------------------------------
// Normalization:  1 for ordered links
// Added by Ping for anisotropic lattices
/*!
  At a site \a x and in the \f$ \mu-\nu \f$ plane, the plaquette is
\f[
  U_\mu(x) U_\nu(x+\nu) U^\dagger_\mu(x+\nu) U^\dagger_\nu(x)
\f]
  
  The sum is over all local lattice sites and all three \a \mu-\nu planes
  where one of \f$ \mu \f$ or \f$ \nu \f$ is the anisotropic direction. The bare
  anisotropy is taken into account here.

  \return The real trace of the plaquette averaged over local sites, planes and colours.
*/
//------------------------------------------------------------------
Float Lattice::AveReTrPlaqXi() const
{
  const char *fname = "AveReTrPlaqXi()";
  VRB.Func(cname,fname);
  
  Float sum = AveReTrPlaqNodeXi();
  glb_sum(&sum);
  return (sum * GJP.VolNodeSites() / GJP.VolSites());
}

//------------------------------------------------------------------
// Added by Schmidt for anisotropic lattices / Thermo
//------------------------------------------------------------------
enum {NUM_SPACE_RECT = 6, //!< Number of planes in a 3-dimensional lattice slice times two.
      NUM_TIME_RECT = 3, //!< Number of planes in a 3-dimensional lattice slice.
};

//------------------------------------------------------------------
// Float AveReTrRectNodeNoXi(void) const
//------------------------------------------------------------------
// Normalization:  1 for ordered links
// Added by Schmidt for anisotropic lattices
/*!
  At a site \a x and in the \f$ \mu-\nu plane\f$, the rectangle with the long axis of the rectangle in the \f$ \mu \f$ direction is:

\f[
  U_\mu(x) U_\mu(x+\mu) U_\nu(x+2\mu) U^\dagger_\mu(x+\mu+\nu) U^\dagger_\mu(x+\nu) U^\dagger_\nu(x)
\f]
  
  The sum is over all local lattice sites and all three \f$x \mu-\nu \f$ planes
  where neither \f$ \mu \f$nor \f$ \nu \f$ is the anisotropic direction and all 
  two rectangles in this plane.

  \return The real trace of the rectangle averaged over local sites,
  planes and colours.

*/
//------------------------------------------------------------------
Float Lattice::AveReTrRectNodeNoXi() const
{
  const char *fname = "AveReTrRectNodeNoXi()";
  VRB.Func(cname,fname);
  
  Float sum = 0.0;
  int x[4];
  
  for(x[0] = 0; x[0] < node_sites[0]; ++x[0]) {
    for(x[1] = 0; x[1] < node_sites[1]; ++x[1]) {
      for(x[2] = 0; x[2] < node_sites[2]; ++x[2]) {
        for(x[3] = 0; x[3] < node_sites[3]; ++x[3]) {
          for (int mu = 0; mu < NUM_DIM; ++mu) {
            for(int nu = mu+1; nu < NUM_DIM; ++nu) {
	      if (mu != GJP.XiDir() && nu != GJP.XiDir()) {
		sum += ReTrRect(x,mu,nu);
		sum += ReTrRect(x,nu,mu);
	      }
            }
          }
        }
      }
    }
  }
  
  return (sum / (GJP.VolNodeSites() * NUM_SPACE_RECT * NUM_COLORS));
}

//------------------------------------------------------------------
// Float AveReTrRectNodeXi1(void) const
//------------------------------------------------------------------
// Normalization:  1 for ordered links
// Added by Schmidt for anisotropic lattices / Thermo
/*!
  At a site \a x and in the \f$ \mu-\nu plane\f$, the rectangle with the long axis of the rectangle in the \f$ \mu \f$ direction is:

\f[
  U_\mu(x) U_\mu(x+\mu) U_\nu(x+2\mu) U^\dagger_\mu(x+\mu+\nu) U^\dagger_\mu(x+\nu) U^\dagger_\nu(x)
\f]
  
  The sum is over all local lattice sites and all three \f$x \mu-\nu \f$ planes
  where the \f$ \nu \f$ (short axis of the rectangle) is the anisotropic direction. The bare anisotropy is taken into account here. 

  \return The real trace of the rectangle averaged over local sites,
  planes and colours.

*/
//------------------------------------------------------------------
Float Lattice::AveReTrRectNodeXi1() const
{
  const char *fname = "AveReTrRectNodeXi1()";
  VRB.Func(cname,fname);
  
  Float sum = 0.0;
  int x[4];
  
  for(x[0] = 0; x[0] < node_sites[0]; ++x[0]) {
    for(x[1] = 0; x[1] < node_sites[1]; ++x[1]) {
      for(x[2] = 0; x[2] < node_sites[2]; ++x[2]) {
        for(x[3] = 0; x[3] < node_sites[3]; ++x[3]) {
          for (int mu = 0; mu < NUM_DIM; ++mu) {
            for(int nu = mu+1; nu < NUM_DIM; ++nu) {
	      if (nu == GJP.XiDir()) 
		sum += ReTrRect(x,mu,nu);
            }
          }
        }
      }
    }
  }
  
  return (sum / (GJP.VolNodeSites() * NUM_TIME_RECT * NUM_COLORS *
		 GJP.XiBare() * GJP.XiBare()));
}

//------------------------------------------------------------------
// Float AveReTrRectNodeXi2(void) const
//------------------------------------------------------------------
// Normalization:  1 for ordered links
// Added by Schmidt for anisotropic lattices / Thermo
/*!
  At a site \a x and in the \f$ \mu-\nu plane\f$, the rectangle with the long axis of the rectangle in the \f$ \mu \f$ direction is:

\f[
  U_\mu(x) U_\mu(x+\mu) U_\nu(x+2\mu) U^\dagger_\mu(x+\mu+\nu) U^\dagger_\mu(x+\nu) U^\dagger_\nu(x)
\f]
  
  The sum is over all local lattice sites and all three \f$x \mu-\nu \f$ planes
  where the \f$ \mu \f$ (long axis of the rectangle) is the anisotropic direction. The bare anisotropy is taken into account here.

  \return The real trace of the rectangle averaged over local sites,
  planes and colours.

*/
//------------------------------------------------------------------
Float Lattice::AveReTrRectNodeXi2() const
{
  const char *fname = "AveReTrRectNodeXi2()";
  VRB.Func(cname,fname);
  
  Float sum = 0.0;
  int x[4];
  
  for(x[0] = 0; x[0] < node_sites[0]; ++x[0]) {
    for(x[1] = 0; x[1] < node_sites[1]; ++x[1]) {
      for(x[2] = 0; x[2] < node_sites[2]; ++x[2]) {
        for(x[3] = 0; x[3] < node_sites[3]; ++x[3]) {
          for (int mu = 0; mu < NUM_DIM; ++mu) {
            for(int nu = mu+1; nu < NUM_DIM; ++nu) {
	      if (mu == GJP.XiDir()) 
		sum += ReTrRect(x,mu,nu);
            }
          }
        }
      }
    }
  }
  
  return (sum / (GJP.VolNodeSites() * NUM_TIME_RECT * NUM_COLORS *
		 GJP.XiBare() * GJP.XiBare() * GJP.XiBare() * GJP.XiBare()));
}

//------------------------------------------------------------------
// Float AveReTrRectNoXi(void) const
//------------------------------------------------------------------
// Normalization:  1 for ordered links
// Added by Schmidt for anisotropic lattices / Thermo
/*!
  At a site \a x and in the \f$ \mu-\nu plane\f$, the rectangle with the long axis of the rectangle in the \f$ \mu \f$ direction is:

\f[
  U_\mu(x) U_\mu(x+\mu) U_\nu(x+2\mu) U^\dagger_\mu(x+\mu+\nu) U^\dagger_\mu(x+\nu) U^\dagger_\nu(x)
\f]
  
  The sum is over all lattice sites and all three \f$x \mu-\nu \f$ planes
  where neither \f$ \mu \f$nor \f$ \nu \f$ is the anisotropic direction and all 
  two rectangles in this plane.

  \return The real trace of the rectangle averaged over local sites,
  planes and colours.

*/
//------------------------------------------------------------------
Float Lattice::AveReTrRectNoXi() const
{
  const char *fname = "AveReTrRectNoXi()";
  VRB.Func(cname,fname);

  Float sum = AveReTrRectNodeNoXi();
  glb_sum(&sum);
  return (sum * GJP.VolNodeSites() / GJP.VolSites());
}

//------------------------------------------------------------------
// Float AveReTrRectXi1(void) const
//------------------------------------------------------------------
// Normalization:  1 for ordered links
// Added by Schmidt for anisotropic lattices / Thermo
/*!
  At a site \a x and in the \f$ \mu-\nu plane\f$, the rectangle with the long axis of the rectangle in the \f$ \mu \f$ direction is:

\f[
  U_\mu(x) U_\mu(x+\mu) U_\nu(x+2\mu) U^\dagger_\mu(x+\mu+\nu) U^\dagger_\mu(x+\nu) U^\dagger_\nu(x)
\f]
  
  The sum is over all lattice sites and all three \f$x \mu-\nu \f$ planes
  where the \f$ \nu \f$ (short axis of the rectangle) is the anisotropic direction. The bare anisotropy is taken into account here. 

  \return The real trace of the rectangle averaged over local sites,
  planes and colours.

*/
//------------------------------------------------------------------
Float Lattice::AveReTrRectXi1() const
{
  const char *fname = "AveReTrRectXi1()";
  VRB.Func(cname,fname);
  
  Float sum = AveReTrRectNodeXi1();
  glb_sum(&sum);
  return (sum * GJP.VolNodeSites() / GJP.VolSites());
}

//------------------------------------------------------------------
// Float AveReTrRectXi2(void) const
//------------------------------------------------------------------
// Normalization:  1 for ordered links
// Added by Schmidt for anisotropic lattices / Thermo
/*!
  At a site \a x and in the \f$ \mu-\nu plane\f$, the rectangle with the long axis of the rectangle in the \f$ \mu \f$ direction is:

\f[
  U_\mu(x) U_\mu(x+\mu) U_\nu(x+2\mu) U^\dagger_\mu(x+\mu+\nu) U^\dagger_\mu(x+\nu) U^\dagger_\nu(x)
\f]
  
  The sum is over all lattice sites and all three \f$x \mu-\nu \f$ planes
  where the \f$ \mu \f$ (long axis of the rectangle) is the anisotropic direction. The bare anisotropy is taken into account here. 

  \return The real trace of the rectangle averaged over local sites,
  planes and colours.

*/
//------------------------------------------------------------------
Float Lattice::AveReTrRectXi2() const
{
  const char *fname = "AveReTrRectXi2()";
  VRB.Func(cname,fname);
  
  Float sum = AveReTrRectNodeXi2();
  glb_sum(&sum);
  return (sum * GJP.VolNodeSites() / GJP.VolSites());
}


//------------------------------------------------------------------
// Lattice::MltFloat(Float factor, int dir = 3)
// Added by Ping for anisotropic lattices
//------------------------------------------------------------------
// Canonical order:    U[t][z][y][x][x,y,z,t]
// Wilson order:       U[t][z][y][x][x,y,z,t] EVEN sites first!
//        Wilson spinors: checkerboarded and  ODD  sites first!
// Staggered order:    U[t][z][y][x][x,y,z,t]
//     same as canonical order up to some phase factor
// G_WILSON_HB:        U[t][z][y][x][t,x,y,z]
  /*!
    \param factor The real scale factor.
    \param dir The direction index of the links to be scaled; 
    \a dir = 0, 1, 2 or 3 for direction X, Y, Z and T respectively (all kinds
    of storage order are handled correctly).
    \post The  gauge field  links in direction \a dir are scaled.
  */
//------------------------------------------------------------------
void Lattice::MltFloatImpl(Float factor, int dir)
{
  if (str_ord == G_WILSON_HB) 
    ++dir;
  dir %= 4;
 
  int nodeSites = GJP.VolNodeSites();
  IFloat *u = (IFloat *)GaugeField() + dir * MATRIX_SIZE;
 
  for (int i = 0; i < nodeSites; ++i) {
    vecTimesEquFloat(u, factor, MATRIX_SIZE);    
    u += 4*MATRIX_SIZE;    
  }
}


//------------------------------------------------------------------
// EvolveGfield(Matrix *mom, Float step_size):
// It evolves the gauge field by step_size using
// the canonical momentum mom
/*!
  Updates each gauge link U according to the canonical equation of motion

  <em>
  U(t+dt) = exp(i dt H) U(t)
  </em>
  
  A ninth order Horner expansion is used to compute the exponential.
  \param mom The multiple \a iH of the conjugate momentum field \a H.
  \param step_size The molecular dynamics time-step \a dt size used in the
  numerical integration of the equations of motion.  
*/  
//------------------------------------------------------------------
void Lattice::EvolveGfield(Matrix *mom, Float step_size)
{
    const char *fname = "EvolveGfield(M*,F)";
    VRB.Func(cname,fname);

    setCbufCntrlReg(4, CBUF_MODE4);

    int n_links = 4 * GJP.VolNodeSites();

#ifdef UNIFORM_SEED_TESTING
    VRB.Result(cname,fname,"gauge checksum(before) = %p\n",
               global_checksum((Float *)GaugeField(),n_links*MATRIX_SIZE));
#endif

    // checksuming local gauge matrices before update
    //-------------------------------------------------
    unsigned long loc_sum = local_checksum((Float *)GaugeField(),n_links*MATRIX_SIZE);

    // checksuming local momentum matrices
    //----------------------------------------------
    loc_sum = local_checksum((Float *)mom,n_links*MATRIX_SIZE);
    CSM.SaveCsum(CSUM_EVL_MOM,loc_sum);

    Matrix *curU_p = GaugeField();

    // Hantao: no problem with this since there are no thread reductions.
#pragma omp parallel for
    for(int i = 0; i < n_links; ++i) {
        Matrix t1, t2, t3;

        t1 = mom[i];
        t1 *= step_size;
        t2 = t1;

        for(int j = 9; j > 1; --j) {
            // t3 = 1 + (1/j) * t2
            oneMinusfTimesMatrix((Float *)&t3, -1./j, (const Float *)&t2, 18);
            // t2 = t1 * t3
            t2.DotMEqual(t1, t3);
        }

        // t3 = 1 + t2
        oneMinusfTimesMatrix((Float *)&t3, -1., (const Float *)&t2, 18);
        t2 = curU_p[i];
        // U' = t3 * U
        curU_p[i].DotMEqual(t3, t2);
    }

    // checksuming local gauge matrices after update
    //------------------------------------------------
    loc_sum = local_checksum((Float *)GaugeField(),n_links*MATRIX_SIZE);
    CSM.SaveCsum(CSUM_EVL_LAT,loc_sum);

#ifdef UNIFORM_SEED_TESTING
    VRB.Result(cname,fname,"gauge checksum(after) = %p\n",
               global_checksum((Float *)GaugeField(),n_links*MATRIX_SIZE));
#endif

    smeared = 0;
}


//------------------------------------------------------------------
// Float MomHamiltonNode(Matrix *momentum):
/*!
  The local kinetic energy is 
  \f[
  \sum_{x, \mu} \frac{1}{2} Tr H_\mu(x)^2
  \f]
  where the \a H is the conjugate momentum field and the
  sum is over all directions \f$ \mu \f$ and all local lattice sites \a x.
  \param momentum  The antihermitian momentum field \a iH.
  \return the local kinetic energy.
 */
//------------------------------------------------------------------
Float Lattice::MomHamiltonNode(Matrix *momentum){
  const char *fname = "MomHamiltonNode(M*)";
  VRB.Func(cname,fname);

  Float ham = 0.0;
  int n_links = 4 * GJP.VolNodeSites();

  // We remove the openmp directive since it may give different
  // answers even on an identical set of data.

  //#pragma omp parallel for reduction(+:ham)
  for(int i = 0; i < n_links; ++i) {
    ham += momentum[i].NegHalfTrSquare();
  }
  
  return ham;
}


//------------------------------------------------------------------
// void Reunitarize():
// Re-unitarize the gauge field configuration.
//------------------------------------------------------------------
// Modefied by Ping on January 1999 for anisotropic lattices.
//------------------------------------------------------------------
/*!
  \post The gauge field is reunitarised.
 */
void Lattice::Reunitarize(void)
{
  const char *fname = "Reunitarize()";
  VRB.Func(cname,fname);
  int i;
  int links;
  Matrix *u;

  links = 4 * GJP.VolNodeSites();
  u = GaugeField();

  for(i=0; i<links; i++){
    u[i].Unitarize();
  }

  // Modified here by Ping for anisotropic lattices
  //----------------------------------------------------------------
  MltFloat(GJP.XiBare(), GJP.XiDir());  
  smeared = 0;
}


//------------------------------------------------------------------
// void Lattice::Reunitarize(Float &dev, Float &max_diff):
/*! 
 \param dev  
 \param max_diff  
 \post The gauge field is reunitarised.
 \post The (averaged) L-2 norm of the resulting change in the gauge field,
\f[
            \surd\{ \sum_i [ (U(i) - V(i))^2 ] / (Vol\times 4\times 18) \}
\f] 
 where U and V are the gauge field before and after 
 reunitarization,  the index \e i runs over all components
 (link direction, local lattice site and colour indices) of the gauge field.
 and \e Vol is the local lattice volume,
 is assigned to \a dev.
\post The L-infinity norm of the resulting change in the gauge field,
\f[
             \max_i( |U(i) - V(i)| )
\f]
 is assigned to \a max_diff.
*/
//------------------------------------------------------------------
// Modefied by Ping on January 1999 for anisotropic lattices.
//------------------------------------------------------------------
void Lattice::Reunitarize(Float &dev, Float &max_diff)
{
  const char *fname = "Reunitarize(F&,F&)";
  VRB.Func(cname,fname);

  // Modified here by Ping for anisotropic lattices
  //----------------------------------------------------------------
  MltFloat(1.0 / GJP.XiBare(), GJP.XiDir());    


  int i,j;
  int links;
  Matrix *u;
  Matrix tmp;
  Float *tmp_p = (Float *)&tmp;

  links = 4 * GJP.VolNodeSites();
  u = GaugeField();

  dev = 0.0;
  max_diff = 0.0;

  for(i=0; i<links; i++){
    tmp = u[i];
    u[i].Unitarize();
    tmp -= u[i];
    for(j=0; j<18; j++){
      dev = dev + tmp_p[j] * tmp_p[j];
      if(tmp_p[j]*tmp_p[j] > max_diff*max_diff){
	max_diff = tmp_p[j];
      }
    }
  }
  dev = dev / Float(18 * links);
  dev = double(sqrt(double(dev)));

  //----------------------------------------------------------------
  // Modified here by Ping for anisotropic lattices
  //----------------------------------------------------------------
  MltFloat(GJP.XiBare(), GJP.XiDir());    
  smeared = 0;
}


//------------------------------------------------------------------
/*!
  \param delta_h The energy difference.
  \param accept The acceptance probability
  \return True (1) if accepted, false (0) otherwise.

  If \a delta_h is greater than or equal to 20 the routine always rejects.
  \post The acceptance probability is written to \a accept.
 */
//------------------------------------------------------------------
int Lattice::MetropolisAccept(Float delta_h, Float *accept) 
{
  const char *fname = "MetropolisAccept(F)";
  VRB.Func(cname,fname);
  int node_id;
  Float flag = 0.0;

  node_id  = GJP.XnodeCoor();
  node_id += GJP.YnodeCoor();
  node_id += GJP.ZnodeCoor();
  node_id += GJP.TnodeCoor();

  // check that delta_h is the same across s-slices
  // (relevant if GJP.Snodes() != 1)
  if(GJP.Snodes() != 1) {
    VRB.Flow(cname,fname, "Checking Delta H across s-slices\n");
    SoCheck(delta_h);
  }

  if (node_id == 0) {
    if(delta_h <= 0.0) {
      flag = 1.0;
      *accept = 1.0;
//    } else if (delta_h < 20.0) {
    } else {
      LRG.AssignGenerator(0, 0, 0, 0);
      LRG.SetInterval(1, 0);
      IFloat exp_mdh;
      IFloat rand_num;

#ifdef _TARTAN
      exp_mdh = double ( exp(- double(delta_h) ) );
#else
      exp_mdh = exp(-delta_h);
#endif
      *accept = (Float)exp_mdh;
      rand_num = LRG.Urand();
      if( rand_num <= exp_mdh ) flag = 1.0;
    }
  }

  // broadcast flag through all nodes
  glb_sum(&flag);

  // check that flag is the same across s-slices
  // (relevant if GJP.Snodes() != 1)
  if(GJP.Snodes() != 1) {
    VRB.Flow(cname,fname, "Checking metropolis flag across s-slices\n");
    SoCheck(flag);
  }

  if(flag < 0.2) {
    return 0;       // Reject
  }
  else {
    return 1;       // Accept
  }

}


//------------------------------------------------------------------
// int MetropolisAccept(Float delta_h, Float *accept):
// 0 reject, 1 accept. If delta_h < 0 it accepts unconditionally.
/*!
  \param delta_h The energy difference.
  \return True (1) if accepted, false (0) otherwise.

  If \a delta_h is greater than or equal to 20 the routine always rejects.
  (Included for backwards compatability)
 */
//------------------------------------------------------------------
int Lattice::MetropolisAccept(Float delta_h)
{
  Float accept;
  return MetropolisAccept(delta_h, &accept);
}



//------------------------------------------------------------------
// void RandGaussAntiHermMatrix(Matrix *mat, Float sigma):
/*!
  Creates a field of antihermitian 3x3 complex matrices with each complex
  element drawn at random from a gaussian distribution with zero mean.
  Hence the matrices are distributed according to

  <em>
  exp[- Tr(mat^2)/(2 sigma2)]
  </em>
  \param mat The field.
  \param sigma2 The variance of the gaussian distribution.
*/
//------------------------------------------------------------------
void Lattice::RandGaussAntiHermMatrix(Matrix *mat, Float sigma2)
{
  const char *fname = "RandGaussAntiHermMatrix(M*,F)";
  VRB.Func(cname,fname);

#if TARGET == QCDSP
    IFloat *a = (IFloat *)CRAM_SCRATCH_ADDR;
#else
    IFloat a[8];
#endif

    LRG.SetSigma(0.5*sigma2);

    Matrix *p = mat;
    for(int n = 0; n < GJP.VolNodeSites(); n++) {
	LRG.AssignGenerator(n);
	for(int j = 0; j < 4; j++) {
	    for(int i = 0; i < 8; ++i) {
		*(a+i) = LRG.Grand(FOUR_D);
	    }
	    p->AntiHermMatrix(a);
	    p++;
	}
    }
}


//--------------------------------------------------------------------------
/*!
  The field is defined on all lattice sites.
  \param frm  A vector.
  \param sigma2  The variance of the gaussian distribution from which the
  vector elements will be drawn.
  \param frm_dim This should be set to ::FOUR_D if the lattice is  
  4-dimensional. The default is ::FIVE_D, \e i.e. a 5-dimensional lattice for
  domain-wall fermions.
  \post The real and imaginary parts of each element of this vector are drawn
  at random from a gaussian distribution with mean 0 and variance \a sigma2. 
 */  
//--------------------------------------------------------------------------
void Lattice::RandGaussVector(Vector *frm, Float sigma2,
                          FermionFieldDimension frm_dim)
{
  RandGaussVector(frm, sigma2, 2, frm_dim);
}
//--------------------------------------------------------------------------
/*!
  The field is defined on all sites of a 5-dimensional lattice for
  domain-wall fermions.
  \param frm  A vector.
  \param sigma2  The variance of the gaussian distribution from which the
  vector elements will be drawn.
  \post The real and imaginary parts of each element of this vector are drawn
  at random from a gaussian distribution with mean 0 and variance \a sigma2. 
 */  
//--------------------------------------------------------------------------
void Lattice::RandGaussVector(Vector *frm, Float sigma2)
{
  RandGaussVector(frm, sigma2, 2, CANONICAL, FIVE_D);
}

void Lattice::RandGaussVector(Vector *frm, Float sigma2, int
chkbds,FermionFieldDimension frm_field_dim  )
{
  if (FstagType())
	RandGaussVector(frm, sigma2, chkbds, STAG, frm_field_dim);
  else
	RandGaussVector(frm, sigma2, chkbds, CANONICAL, frm_field_dim);
}
 
//--------------------------------------------------------------------------
/*!
  \param frm  A vector.
  \param sigma2  The variance of the gaussian distribution from which the
  vector elements will be drawn.
  \param num_chkbds This should be set to 2 if the field is defined on all
  lattice sites in canonical order or 1 if the field is defined on lattice
  sites of a single parity.
  \param frm_dim This should be set to ::FOUR_D if the lattice is  
  4-dimensional. The default is ::FIVE_D, \e i.e. a 5-dimensional lattice for
  domain-wall fermions.
  \post The real and imaginary parts of each element of this vector are drawn
  at random from a gaussian distribution with mean 0 and variance \a sigma2. 
 */  
//--------------------------------------------------------------------------
void Lattice::RandGaussVector(Vector * frm, Float sigma2, int num_chkbds,
             StrOrdType str, FermionFieldDimension frm_dim /* = FIVE_D */ )
{
  const char *fname = "RandGaussVector(Vector *, Float, int, FermionFieldDimension)";
  VRB.Func(cname, fname);

  int vec_size = 2 * Colors() * SpinComponents();

  int s_node_sites = GJP.SnodeSites();
  if(frm_dim == FOUR_D
     || s_node_sites == 0
     // FIXME: checking Fclass() is a bad idea, replace it with something more reasonable.
     || (Fclass() != F_CLASS_DWF && Fclass() != F_CLASS_BFM)) {
    s_node_sites = 1; frm_dim = FOUR_D;
  }
  LRG.SetSigma(sigma2);

  IFloat * ptr = (IFloat *) frm;
  int checker, i, k, s, x[4];
  IFloat sum=0.0,square=0.0;


  if(num_chkbds == 2) {
    for(checker = 0; checker < 2; checker++)
    for(s = 0; s < s_node_sites; s++) {
      if( (s % 2) == checker) {
        for(x[3] = 0; x[3] < GJP.TnodeSites(); x[3]++)
        for(x[2] = 0; x[2] < GJP.ZnodeSites(); x[2]++)
        for(x[1] = 0; x[1] < GJP.YnodeSites(); x[1]++)
        for(x[0] = 0; x[0] < GJP.XnodeSites(); x[0]++) {
//	printf("%d %d %d %d %d \n",x[0],x[1],x[2],x[3],s);
          LRG.AssignGenerator(x[0],x[1],x[2],x[3],s);
//	printf("%d %d %d %d %d \n",x[0],x[1],x[2],x[3],s);
          for(k = 0; k < vec_size; k++) {
            *(ptr++) = LRG.Grand(frm_dim);
          }
        }
      }
    }
  }
  else if(num_chkbds == 1) {
    if (str == STAG){           
      for(x[2] = 0; x[2] < GJP.ZnodeSites(); x[2]++)     // z
      for(x[1] = 0; x[1] < GJP.YnodeSites(); x[1]++)     // y
      for(x[0] = 0; x[0] < GJP.XnodeSites(); x[0]++)     // x
      for(x[3] = 0; x[3] < GJP.TnodeSites(); x[3] += 2) {   // t
        LRG.AssignGenerator(x);
        for(k = 0; k < FsiteSize(); k++) {
          *(ptr) = LRG.Grand(frm_dim);
          sum += *ptr;
          square += (*ptr)*(*ptr);
          ptr++;
        }
      }
    } else {
      for(i = 0; i < GJP.VolNodeSites()*s_node_sites; i+=2) {
        LRG.AssignGenerator(i);
        for(k = 0; k < vec_size; k++) {
          *(ptr) = LRG.Grand(frm_dim);
          sum += *ptr;
          square += (*ptr)*(*ptr);
 	  ptr++;
        }
      }
    }
  }
#if 0
  glb_sum_five(&sum);
  glb_sum_five(&square);
  printf("sum=%0.18e square=%0.18e\n",sum,square);
#endif
}


//------------------------------------------------------------------
/*
  \post Each gauge field link is set to the identity matrix.
*/
//------------------------------------------------------------------
void Lattice::SetGfieldOrd(void){
  const char *fname = "SetGfieldOrd()";
  VRB.Func(cname,fname);
  int i;
  int links;
  Matrix *u;

  links = 4 * GJP.VolNodeSites();
  u = GaugeField();

  for(i=0; i<links; i++){
    u[i].UnitMatrix();
  }
}

//------------------------------------------------------------------
/*
  \post Each u1 field(angle) is set to zero.
*/
//------------------------------------------------------------------
void Lattice::SetU1GfieldOrd(void){
  char *fname = "SetGfieldOrd()";
  VRB.Func(cname,fname);
  int i;
  int links;
  Float *u;

  links = 4 * GJP.VolNodeSites();
  u = U1GaugeField();

  for(i=0; i<links; i++){
    u[i]=0.0;
  }
}


//------------------------------------------------------------------
/*
  \post Each gauge field link is set to a random SU(3) matrix.
*/
//------------------------------------------------------------------
void Lattice::SetGfieldDisOrd(void){
#if 1
  const char *fname = "SetGfieldDisOrd()";
  VRB.Func(cname,fname);

  int site_size = GsiteSize();
  LRG.SetInterval(1, -1);

  IFloat *pmat=(IFloat *) gauge_field;

  for(int i=0; i<GJP.VolNodeSites(); i++) {
    LRG.AssignGenerator(i);
    for(int k = 0; k < site_size; k++) {
      *(pmat++) = LRG.Urand(FOUR_D);
//      printf("i=%d *pmat=%e\n",i,*(pmat-1));
    }
  }
  Reunitarize();
#endif
  smeared=0;

}


//------------------------------------------------------------------
// Counter related functions (g_upd_cnt, md_time)
//------------------------------------------------------------------

/*!
  \return The number of gauge field updates that have been performed.
*/
int Lattice::GupdCnt(void)
{
  return *g_upd_cnt ;
}

/*!
  Sets the the initial value of the counter of the number of gauge field
  updates that have been performed.
  \param set_val The inital value for the counter.
  \return The inital value for the counter.
*/
int Lattice::GupdCnt(int set_val)
{
  *g_upd_cnt = set_val ;
  return *g_upd_cnt ;
}

/*!
  \param inc_val The number by which to inrease the gauge field update counter.
  \return The new value of the gauge field update counter.
*/
int Lattice::GupdCntInc(int inc_val)
{
  *g_upd_cnt += inc_val ;
  return *g_upd_cnt ;
}

//! The molecular dynamics time counter.
/*!
  \return The number of timesteps completed so far in a molecular dynamics
  trajectory.
*/
Float Lattice::MdTime(void)
{
  return md_time ;
}

//! Sets the value of the molecular dynamics time counter.
/*!
  Sets the number of timesteps completed so far in a molecular dynamics
  trajectory.
*/
Float Lattice::MdTime(Float set_val)
{
  md_time = set_val ;
  return md_time ;
}

//! Increments the value of the molecular dynamics time counter.
/*!
  Increments the number of timesteps completed so far in a molecular dynamics
  trajectory.
*/
Float Lattice::MdTimeInc(Float inc_val)
{
  md_time += inc_val ;
  return md_time ;
}


//------------------------------------------------------------------
//! Returns an array of pointers to the gauge fixed hyperplanes.
/*!
  The array are allocated with FixGaugeAllocate.
*/
//------------------------------------------------------------------
Matrix **Lattice::FixGaugePtr(void){
  return fix_gauge_ptr;
}


//------------------------------------------------------------------
// Returns fix_gauge_kind (the kind of gauge).
//------------------------------------------------------------------
FixGaugeType Lattice::FixGaugeKind(void){
  return fix_gauge_kind;
}

Float Lattice::FixGaugeStopCond(void){
  return fix_gauge_stpCond;
}


//==================================================================
// Gauge action related virtual functions implemented 
// as dummy functions that exit with an error.
// These functions are "truly" implemented  
// only at some derived class or classes but not
// at all of the derived classes and this is why they are not
// implemented as pure virtual.
//==================================================================

//------------------------------------------------------------------
// v_out = gamma_5 v_in
//------------------------------------------------------------------
void Lattice::Gamma5(Vector *v_out, Vector *v_in, int num_sites){
  const char *fname = "Gamma5";
  ERR.NotImplemented(cname,fname);
}

//==================================================================
// Fermion action related virtual functions implemented
// as dummy functions that exit with an error.
// These functions are "truly" implemented  
// only at some derived class or classes but not
// at all of the derived classes and this is why they are not
// implemented as pure virtual.
//==================================================================

//------------------------------------------------------------------
// This is "truly" implemented only in the Fdwf derived class
//------------------------------------------------------------------
void Lattice::Ffour2five(Vector *five, Vector *four, 
			 int s_r, int s_l, int Ncb)
{
  const char *fname = "Ffour2five";
  ERR.NotImplemented(cname,fname);
}

//------------------------------------------------------------------
// This is "truly" implemented only in the Fdwf derived class
//------------------------------------------------------------------
void Lattice::Ffive2four(Vector *four, Vector *five, 
			 int s_r, int s_l, int Ncb)
{
  const char *fname = "Ffive2four";
  ERR.NotImplemented(cname,fname);
}

//----------------------------------------------------------------------
// dir_flag is flag which takes value 0 when all direction contribute to D
// 1 - when only the special anisotropic direction contributes to D,
// 2 - when all  except the special anisotropic direction.
// Currently this function is implemented only in the Fstag class
//------------------------------------------------------------------
void Lattice::Fdslash(Vector *f_out, Vector *f_in, CgArg *cg_arg, 
		    CnvFrmType cnv_frm, int dir_flag)
{
  const char *fname = "Fdslash";
  ERR.NotImplemented(cname,fname);
}

//----------------------------------------------------------------------
// order is the order of the derivative  with respect to the chemical
// potential. 
// Currently this function is implemented only in the Fstag class
//------------------------------------------------------------------

void Lattice::FdMdmu(Vector *f_out, Vector *f_in, CgArg *cg_arg, 
		    CnvFrmType cnv_frm, int order)
{
  const char *fname = "FdMdmu";
  ERR.NotImplemented(cname,fname);
}

//~~------------------------------------------------------------------
//~~ Added by DHR for twisted-mass fermions
//~~ For most functions, twisted mass parameter passed in cg_arg
//~~ For functions not using cg_arg, added explicit twisted-mass
//~~ parameter to function call parameters
//~~ These functions implemented in lattice_base as an error, and
//~~ "truly" implemented only in the F_wilsonTm derived class
//~~------------------------------------------------------------------

Float Lattice::SetPhi(Vector *phi, Vector *frm1, Vector *frm2,
    		Float mass, Float epsilon, DagType dag) { 
  const char *fname = "SetPhi(V*,V*,V*,F,F,DagType)";
  ERR.NotImplemented(cname,fname);
  return ((Float) 0.0);
};

ForceArg Lattice::EvolveMomFforce(Matrix *mom, Vector *frm, 
				 Float mass, Float epsilon, Float step_size) {
  const char *fname = "EvolveMomFforce(M*,V*,F,F,F)";
  ERR.NotImplemented(cname,fname);
  return ForceArg(0.0,0.0,0.0);
};

ForceArg Lattice::EvolveMomFforce(Matrix *mom, Vector *phi, Vector *eta,
		 Float mass, Float epsilon, Float step_size) {
  const char *fname = "EvolveMomFforce(M*,V*,V*,F,F,F)";
  ERR.NotImplemented(cname,fname);
  return ForceArg(0.0,0.0,0.0);
};

Float Lattice::BhamiltonNode(Vector *boson, Float mass, Float epsilon)  {
  const char *fname = "BhamiltonNode(V*,F,F)";
  ERR.NotImplemented(cname,fname);
  return ((Float) 0.0);
};


//------------------------------------------------------------------
//void *Lattice::Aux0Ptr(void):
// Returns the general purpose auxiliary pointer 0;
//------------------------------------------------------------------
void *Lattice::Aux0Ptr(void){
  return aux0_ptr;
}

//------------------------------------------------------------------
//void *Lattice::Aux1Ptr(void):
// Returns the general purpose auxiliary pointer 1;
//------------------------------------------------------------------
void *Lattice::Aux1Ptr(void){
  return aux1_ptr;
}


//------------------------------------------------------------------
/*!
  Checks that the gauge field is
  identical (using a checksum and the average plaquette)
  on each slice of the lattice perpendicular to the 5th direction
  and local in the 5th direction.
  Obviously this is always the case when the entire 5th direction is local.
  If any of the node slices fail to match the program exits with an error.
*/
//------------------------------------------------------------------
void Lattice::GsoCheck(void){
  const char *fname = "GsoCheck()";
  VRB.Func(cname,fname);
  int s;
  IFloat rcv_buf[2];
  IFloat snd_buf[2];
  IFloat failed_value;
  Float failed_flag;
  int failed_slice;
  int s_nodes = GJP.Snodes();

  // If GJP.Snodes() != 1 do the check
  //----------------------------------------------------------------
  if(s_nodes != 1) {
#if 1
    ERR.General(cname,fname,"s_nodes(%d)!=1!\n",s_nodes);
#endif

    // Calculate plaquette
    //--------------------------------------------------------------
    IFloat plaq = SumReTrPlaqNode();
    
    // Compare plaquettes along fifth direction and set failed_flag
    //--------------------------------------------------------------
    failed_flag = 0.0;
    failed_slice = 0;
    snd_buf[0] = plaq;
    for(s=1; s < s_nodes; s++){
      getMinusData(rcv_buf, snd_buf, 2, 4);
      if(rcv_buf[0] != plaq) {
	printf("plaq=%e rcv_buf=%e\n",plaq,rcv_buf[0]);
	failed_flag = 1.0;
	ERR.General(cname,fname,"GsoCheck: PLAQUETTE TEST FAILED");
	if(failed_slice == 0) {
	  failed_slice = s;
	  failed_value = rcv_buf[0];
	}
      }
      snd_buf[0] = rcv_buf[0];
    }
    glb_sum_five(&failed_flag);
    if(failed_flag > 0.1) {
      ERR.General(cname,fname, 
      "Node (%d,%d,%d,%d,%d): plaquette test failed:\n\tlocal plaq = %e, but at -%d slices plaq = %e\n",
		  GJP.XnodeCoor(),
		  GJP.YnodeCoor(),
		  GJP.ZnodeCoor(),
		  GJP.TnodeCoor(),
		  GJP.SnodeCoor(),
		  plaq,
		  failed_slice,
		  failed_value);
    }
    
    // Calculate gauge field checksum
    //--------------------------------------------------------------
    Float *u = (Float *) gauge_field;
    int size = GsiteSize() * GJP.VolNodeSites();  
    IFloat checksum = 0.0;
    for(int ic=0; ic < size; ic++){
      checksum = checksum + u[ic];
    }

    // Compare checksum along fifth direction and set failed_flag
    //--------------------------------------------------------------
    failed_flag = 0.0;
    failed_slice = 0;
    snd_buf[0] = checksum;
    for(s=1; s < s_nodes; s++){
      getMinusData(rcv_buf, snd_buf, 2, 4);
      if(rcv_buf[0] != checksum) {
	failed_flag = 1.0;
#ifdef _TARTAN
	InterruptExit(-1, "GsoCheck: CHECKSUM TEST FAILED");
#endif
	if(failed_slice == 0) {
	  failed_slice = s;
	  failed_value = rcv_buf[0];
	}
      }
      snd_buf[0] = rcv_buf[0];
    }
    glb_sum_five(&failed_flag);
    if(failed_flag > 0.1) {
      ERR.General(cname,fname, 
      "Node (%d,%d,%d,%d,%d): checksum test failed:\n\tlocal checksum = %e, but at -%d slices checksum = %e\n",
		  GJP.XnodeCoor(),
		  GJP.YnodeCoor(),
		  GJP.ZnodeCoor(),
		  GJP.TnodeCoor(),
		  GJP.SnodeCoor(),
		  checksum,
		  failed_slice,
		  failed_value);
    }

  }

  VRB.Flow(cname,fname,"Plaquette, checksum test successful\n");

}



//------------------------------------------------------------------
/*!
  Checks that a floating point number is
  identical (using a checksum and the average plaquette)
  on each slice of the lattice perpendicular to the 5th direction
  and local in the 5th direction. 
  Obviously this is always the case when the entire 5th direction is local.
  If any of the node slices fail to match the program exits with an error.
  \param num The number to check.
*///------------------------------------------------------------------
void Lattice::SoCheck(Float num){
  const char *fname = "SoCheck()";
  VRB.Func(cname,fname);
  int s;
  IFloat rcv_buf;
  IFloat snd_buf;
  IFloat failed_value;
  Float failed_flag;
  int failed_slice;
  int s_nodes = GJP.Snodes();
  IFloat number = num;

  // If GJP.Snodes() != 1 do the check
  //----------------------------------------------------------------
  if(s_nodes != 1) {
    
    // Compare numbers along fifth direction and set failed_flag
    //--------------------------------------------------------------
    failed_flag = 0.0;
    failed_slice = 0;
    snd_buf = number;
    for(s=1; s < s_nodes; s++){
      getMinusData(&rcv_buf, &snd_buf, 1, 4);
      if(rcv_buf != number) {
	failed_flag = 1.0;
#ifdef _TARTAN
	InterruptExit(-1, "GsoCheck: CHECKSUM TEST FAILED");
#endif
	if(failed_slice == 0) {
	  failed_slice = s;
	  failed_value = rcv_buf;
	}
      }
      snd_buf = rcv_buf;
    }
    glb_sum_five(&failed_flag);
    if(failed_flag > 0.1) {
      ERR.General(cname,fname,
      "Node (%d,%d,%d,%d,%d): SoCheck test failed:\n\tlocal number = %e, but at -%d slices number = %e\n",
		  GJP.XnodeCoor(),
		  GJP.YnodeCoor(),
		  GJP.ZnodeCoor(),
		  GJP.TnodeCoor(),
		  GJP.SnodeCoor(),
		  number,
		  failed_slice,
		  failed_value);
    }

  }

  VRB.Flow(cname,fname,"SoCheck test successful\n");

}

void Lattice::Shift()
{
   GDS.Shift(gauge_field, GJP.VolNodeSites()*4*sizeof(Matrix));
}

#if 0
void Lattice::ForceMagnitude(Matrix *mom, Matrix *mom_old, 
			     Float mass, Float dt, char *type) {

  const char *fname = "ForceMagnitude(Matrix*, Matrix*, Float, Float, char*)";
  FILE *fp;


  int g_size = GJP.VolNodeSites() * GsiteSize();
  
  vecMinusEquVec( (IFloat*)mom_old, (const IFloat*)mom, g_size);
  Float force;
  if (Fclass() == F_CLASS_DWF)
    force = sqrt(((Vector*)mom_old)->NormSqGlbSum(g_size));
  else
    force = sqrt(((Vector*)mom_old)->NormSqGlbSum4D(g_size));

  // Print out monitor info
  //---------------------------------------------------------------
  if( (fp = Fopen("force.dat", "a")) == NULL ) {
    ERR.FileA(cname,fname, "force.dat");
  }
  Fprintf(fp,"%s = %e mass = %f dt = %f\n", type, IFloat(force), 
	  (IFloat)mass, (IFloat)dt);
  Fclose(fp);

}
#endif

int Lattice::FmatEvlMInv(Vector **f_out, Vector *f_in, Float *shift,
                            int Nshift, int isz, CgArg *cg_arg,
                            CnvFrmType cnv_frm, MultiShiftSolveType type,
                            Float *alpha, Vector **f_out_d){
    CgArg **cg_arg_p = new CgArg *[Nshift];
    for(int i =0;i<Nshift;i++) cg_arg_p[i] = cg_arg;
 
    int iter = FmatEvlMInv(f_out,f_in,shift,Nshift,isz,cg_arg_p,cnv_frm,
                               type,alpha, f_out_d);
    delete[] cg_arg_p;
    return iter;
}

//by Qi Liu
int Lattice::eig_FmatInv(Vector **V, const int vec_len, Float *M, const int nev, const int m, float **U, Rcomplex *invH, const int def_len, const Float *restart,const int restart_len, Vector *f_out, Vector *f_in, CgArg *cg_arg, Float *true_res, CnvFrmType cnv_frm , PreserveType prs_f_in)
{
	 char *cname = "Lattice";
	 const char *fname = "FmatInvProj()";
	 ERR.General(cname,fname,"Only have code for dwf class not others so this is not a pure virtual function\n");
}

unsigned long Lattice::GsiteOffset(const int *x, const int dir) const{
  const char *fname="GsiteOffset(*i,i)";
  int parity = (x[0]+x[1]+x[2]+x[3])%2;
  int vol = GJP.VolNodeSites();
  unsigned long index;
  switch(StrOrd()){
  case WILSON:
// XYZT ordering, checkerboarded, even first
    index = x[3];
    for(int i=2;i>=0;i--) index = x[i]+GJP.NodeSites(i)*index;
    index = (index+vol*parity)/2;  
// dir also XYZT
    index = index * 4 + dir; 
    break;
  default:
    ERR.NotImplemented(cname,fname);
  }
  return index;
}


inline void compute_coord(int x[4], const int hl[4], const int low[4], int i)
{
    x[0] = i % hl[0] + low[0]; i /= hl[0];
    x[1] = i % hl[1] + low[1]; i /= hl[1];
    x[2] = i % hl[2] + low[2]; i /= hl[2];
    x[3] = i % hl[3] + low[3];
}

// ----------------------------------------------------------------
// BondCond: toggle boundary condition on/off for gauge field. Based
// on code from src/util/dirac_op/d_op_base/comsrc/dirac_op_base.C
//
// The gauge field must be in CANONICAL order.
// ----------------------------------------------------------------
void Lattice::BondCond()
{
    Matrix *u_base = this->GaugeField();

    for(int mu = 0; mu < 4; ++mu) {
        if(GJP.NodeBc(mu) != BND_CND_APRD) continue;

        int low[4] = { 0, 0, 0, 0 };
        int high[4] = { GJP.XnodeSites(), GJP.YnodeSites(),
                        GJP.ZnodeSites(), GJP.TnodeSites() };
        low[mu] = high[mu] - 1;

        int hl[4] = { high[0] - low[0], high[1] - low[1],
                      high[2] - low[2], high[3] - low[3] };

        const int hl_sites = hl[0] * hl[1] * hl[2] * hl[3];

#pragma omp parallel for
        for(int i = 0; i < hl_sites; ++i) {
            int x[4];
            compute_coord(x, hl, low, i);

            int off = mu + 4 * (x[0] + high[0] *
                                (x[1] + high[1] *
                                 (x[2] + high[2] * x[3])));
            u_base[off] *= -1.;
        }
    }
}

CPS_END_NAMESPACE
