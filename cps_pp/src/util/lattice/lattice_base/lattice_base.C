#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:47 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/lattice_base/lattice_base.C,v 1.1.1.1 2003-06-22 13:34:47 mcneile Exp $
//  $Id: lattice_base.C,v 1.1.1.1 2003-06-22 13:34:47 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.9  2002/12/04 17:16:27  zs
//  Merged the new 2^4 RNG into the code.
//  This new RNG is implemented in the LatRanGen class.
//  The following algorithm and utility classes are affected:
//
//  AlgEig                  Fdwf
//  AlgGheatBath            Fstag
//  AlgHmd                  GlobalJobParameter
//  AlgNoise                Lattice
//  AlgPbp                  Matrix
//  AlgThreept              RandomGenerator
//                          Vector
//
//  Revision 1.8  2001/09/06 11:51:47  anj
//  Minor modifications to the test suite, e.g. standardizing the
//  verbosity and such.  Collected the output from the original and the
//  latest versions using the new test suite, and checked them. Anj
//
//  Revision 1.7  2001/08/16 10:50:37  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "IFloat".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and IFloat).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.4  2001/07/03 17:01:03  anj
//
//  Multiple minor alterations to change some #include's from referring to
//  files relative to the top-level source directory to referring to files
//  relative to the source-file positions.  This alteration makes the code
//  backwards compatable with the make structure of QCDSP, although this
//  may have to be changed to a more usual form in the future. Anj.
//
//  Revision 1.3  2001/06/21 15:40:18  anj
//  Updated the _TARTAN ifdefs, using PARALLEL instead (where appropriate).Anj
//
//  Revision 1.2  2001/06/19 18:13:29  anj
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
//  Revision 1.2  2001/05/25 06:16:10  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: lattice_base.C,v $
//  $Revision: 1.1.1.1 $
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

CPS_END_NAMESPACE
#include<config.h>
#include<util/lattice.h>
#include<util/vector.h>
#include<util/gjp.h>
#include<util/pmalloc.h>
#include<util/verbose.h>
#include<util/error.h>
#include<util/random.h>
#include<comms/nga_reg.h>
#include<comms/glb.h>
#include<comms/scu.h>
#include<comms/cbuf.h>
CPS_START_NAMESPACE

#ifdef _TARTAN
CPS_END_NAMESPACE
#include <math64.h>
CPS_START_NAMESPACE
#else
CPS_END_NAMESPACE
#include <math.h>
CPS_START_NAMESPACE
#endif

#ifdef PARALLEL
CPS_END_NAMESPACE
#include <sysfunc.h>
CPS_START_NAMESPACE
#endif

//------------------------------------------------------------------
// Define the opposite dir macro for the  Pathorder routines
//------------------------------------------------------------------
#define OPP_DIR(dir)  (((dir)+4)&7)
 
//------------------------------------------------------------------
// Strings
//------------------------------------------------------------------

//------------------------------------------------------------------
enum { MATRIX_SIZE = 18 };
//------------------------------------------------------------------

//------------------------------------------------------------------
// Initialize static variables
//------------------------------------------------------------------
Matrix* Lattice::gauge_field = 0;
int Lattice::is_allocated = 0;
int Lattice::is_initialized = 0;
StrOrdType Lattice::str_ord = CANONICAL;
Matrix** Lattice::fix_gauge_ptr = 0;
FixGaugeType Lattice::fix_gauge_kind = FIX_GAUGE_NONE;

int Lattice::node_sites[5];
int Lattice::g_dir_offset[4];

int* Lattice::g_upd_cnt = 0 ;
Float Lattice::md_time = (Float) 0.0 ;


static Matrix m_tmp1, m_tmp2;	// DRAM temp buffer
				// used for scu transfer

//------------------------------------------------------------------
// static variables used only inside this file
//------------------------------------------------------------------
static IFloat invs[10] = { -100.0,	// not used (1/0)
	-1.0,	-0.5,	-1./3.,	-0.25,	-0.2,
	-1./6.,	-1./7.,	-0.125,	-1./9. };


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



//------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------
Lattice::Lattice()
{
  cname = "Lattice";
  char *fname = "Lattice()";
  int array_size;  // On-node size of the gauge field array.
  RandomGenerator ran_gen; // Random generator

  VRB.Func(cname,fname);

  StartConfType start_conf_kind = GJP.StartConfKind();
  link_buffer = 0;


  if(!(is_allocated)){
    // Allocate memory for the gauge field.
    //--------------------------------------------------------------
    array_size = GsiteSize() * GJP.VolNodeSites() * sizeof(Float);  
    if(start_conf_kind != START_CONF_LOAD){
      gauge_field = (Matrix *) pmalloc(array_size);
      if( gauge_field == 0)
	ERR.Pointer(cname,fname, "gauge_field");
      VRB.Pmalloc(cname, fname,
		  "gauge_field",
		  gauge_field,
		  array_size);
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
    //???
    VRB.Flow(cname,fname, "File starting configuration\n");
    ERR.NotImplemented(cname,fname, 
    "Starting config. type  START_CONF_FILE not implemented\n");
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
}


//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
Lattice::~Lattice()
{
  char *fname = "~Lattice()";
  VRB.Func(cname,fname);
  
  // if GJP.Snodes() != 1 check that the gauge field is identical
  // on all s-slices.
  //----------------------------------------------------------------
  if(GJP.Snodes() != 1) {
    VRB.Flow(cname,fname, "Checking gauge field across s-slices\n");
    GsoCheck();
  }

  // Undo the scaling of links on anisotropic lattices  -- Ping
  //----------------------------------------------------------------
  // No data manipulations if the scaling factor is 1.0
  // [built into lat.MltFloat].
  MltFloat(1.0 / GJP.XiBare(), GJP.XiDir());
}


//------------------------------------------------------------------
// const Matrix *GaugeField(void) const:
// Returns the pointer to the gauge field configuration.
//------------------------------------------------------------------
Matrix *Lattice::GaugeField(void) const
{
  return gauge_field;
}


//------------------------------------------------------------------
// void GaugeField(Matrix *u) const:
// Copies the array pointed to by u to the array pointed to 
// by gauge_field.
//------------------------------------------------------------------
void Lattice::GaugeField(Matrix *u)
{
  char *fname = "GaugeField(M*)";
  VRB.Func(cname,fname);
  int size;

  size = GsiteSize() * GJP.VolNodeSites() * sizeof(IFloat);  

  // Copy from u to gauge_field
  //----------------------------
  moveMem((IFloat *) gauge_field, (IFloat *) u, size);
}


//------------------------------------------------------------------
// void CopyGaugeField(Matrix* u) const:
// Copies the array pointed to by gauge_field to the
// array pointed to by u.
//------------------------------------------------------------------
void Lattice::CopyGaugeField(Matrix* u)
{
  char *fname = "CopyGaugeField(M*)";
  VRB.Func(cname,fname);
  int size;

  size = GsiteSize() * GJP.VolNodeSites() * sizeof(IFloat);  

  // Copy from gauge_field to u
  //----------------------------
  moveMem((IFloat *) u, (IFloat *) gauge_field, size);
}

//------------------------------------------------------------------
// StrOrd(): returns the value of the current storage order.
//------------------------------------------------------------------
StrOrdType Lattice::StrOrd(void)
{
  return str_ord;
}


//------------------------------------------------------------------
// int Colors();
// Returns the number of colors.	  
//------------------------------------------------------------------
int Lattice::Colors(void)
{
  return GJP.Colors();
}


//------------------------------------------------------------------
// int GsiteSize() : 
// Returns the number of gauge field components 
// (including real/imaginary) on a site of the 4-D lattice.
//------------------------------------------------------------------
int Lattice::GsiteSize(void)
{
  return  2 * Colors() * Colors() * 4; //re/im * colors*colors * dim
}


//--------------------------------------------------------------------------
// const Matrix *
// Lattice::GetLink(const int *site, int dir) const
//--------------------------------------------------------------------------
// Purpose:
//   get a link at specified coordinates and direction, whether on node
//   or off node.
// Arguments:
//   site:    lattice coordinates[x,y,z,t]
//            which could be out-of-range, i.e., is located off-node.
//   dir:     0,1,2,3 for Ux, Uy, Uz, Ut
//   return:  a pointer to the link. If off-node, it points to
//            a block of static memory. Be careful to use it!
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
//char *fname = "GetLink()";
//VRB.Func(cname,fname);

//VRB.Debug(cname, fname, "link %3i %3i %3i %3i ; %i\n",
//  site[0], site[1], site[2], site[3], dir) ;

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
//VRB.FuncEnd(cname, fname) ;
  return on_node_link;
#endif

  // send to the destination node if the site is off-node
  //------------------------------------------------------------------------
  if (on_node) {
//  VRB.FuncEnd(cname, fname) ;
    return on_node_link;
  } else {
    Matrix send = *on_node_link;
    Matrix &recv = m_tmp1 ;
    for (int i = 0; i < 4; ++i) {
      while (site[i] != on_node_site[i]) {
        if (site[i] < 0) {
          getMinusData((IFloat *)&recv, (IFloat *)&send, sizeof(recv), i);
          on_node_site[i] -= node_sites[i];
        } else {
          getPlusData((IFloat *)&recv, (IFloat *)&send, sizeof(recv), i);
          on_node_site[i] += node_sites[i];
        }
        send = recv;
      }
    }
//  VRB.FuncEnd(cname, fname) ;
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



const unsigned CBUF_MODE2 = 0xcca52112;
const unsigned CBUF_MODE4 = 0xcca52112;


//------------------------------------------------------------------
// Staple(Matrix& stap, int *x, int mu):
// It calculates the staple field at x, mu.
// The staple field is:
//
//	V_u(x) = \sum_v(!=u) {
//		U_v(x+u) U_u(x+v)~ U_v(x)~
//	     +  U_v(x+u-v)~ U_u(x-v)~ U_v(x-v)  }
//
//------------------------------------------------------------------
void Lattice::Staple(Matrix& stap, int *x, int mu)
{
//char *fname = "Staple(M&,i*,i)";
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
      mp3->Dagger((IFloat *)p1+BANK4_BASE+BANK_SIZE);


      //----------------------------------------------------------
      // p1 = &U_v(x+u)
      //----------------------------------------------------------
      p1 = GetLinkOld(g_offset, x, mu, nu);


      //----------------------------------------------------------
      // mp2 = U_v(x+u) U_u(x+v)~
      //----------------------------------------------------------
      mDotMEqual((IFloat *)mp2, (const IFloat *)p1+BANK2_BASE,
		 (const IFloat *)mp3);
      

      //----------------------------------------------------------
      //  mp3 = U_v(x)~
      //----------------------------------------------------------
      mp3->Dagger((IFloat *)(g_offset+nu)+BANK4_BASE);
      

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
      moveMem((IFloat *)mp3, (IFloat *)p1+BANK2_BASE, 
	      MATRIX_SIZE * sizeof(IFloat));


      //----------------------------------------------------------
      // mp2 = U_u(x-v) U_v(x+u-v)
      //----------------------------------------------------------
      mDotMEqual((IFloat *)mp2, (const IFloat *)(g_offpv+mu)+BANK4_BASE,
		 (const IFloat *)mp3);
      

      //----------------------------------------------------------
      // mp3 = U_v(x+u-v)~ U_u(x-v)~ = mp2~
      //----------------------------------------------------------
      mp3->Dagger((IFloat *)mp2);


      //----------------------------------------------------------
      // mp2 = U_v(x-v)
      //----------------------------------------------------------
      moveMem((IFloat *)mp2, (const IFloat *)(g_offpv+nu)+BANK2_BASE,
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
	vecAddEquVec((IFloat *)&stap, (IFloat *)&m_tmp2+BANK4_BASE,
		MATRIX_SIZE);
	
      } else {
	mDotMPlus((IFloat *)&stap, (const IFloat *)mp3,
		  (const IFloat *)mp2);
        // dummy read
	*((IFloat *)mp2) = *((IFloat *)g_offpv+BANK4_BASE);
      }
    }
  }
}


//------------------------------------------------------------------
// RectStaple(Matrix& stap, int *x, int mu):
// It calculates the rectangle staple field at x, mu.
// The rectangle staple field is:
//
// \sum_{v != u} {
//     U_u(x+u)    U_v(x+2u)    U_u(x+u+v)~ U_u(x+v)~  U_v(x)~
//   + U_u(x+u)    U_v(x+2u-v)~ U_u(x+u-v)~ U_u(x-v)~  U_v(x-v)
//   + U_v(x+u-v)~ U_v(x+u-2v)~ U_u(x-2v)~  U_v(x-2v)  U_v(x-v)
//   + U_v(x+u-v)~ U_u(x-v)~    U_u(x-u-v)~ U_v(x-u-v) U_u(x-u)
//   + U_v(x+u)    U_u(x+v)~    U_u(x-u+v)~ U_v(x-u)~  U_u(x-u)
//   + U_v(x+u)    U_v(x+u+v)   U_u(x+2v)~  U_v(x+v)~  U_v(x)~
// }
//
//------------------------------------------------------------------
void Lattice::RectStaple(Matrix& rect, int *x, int mu)
{
//char *fname = "RectStaple(M&,i*,i)" ;
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
    mp4->Dagger((IFloat *)GetLink(link_site, nu)+BANK4_BASE) ;

    //----------------------------------------------------------
    // mp3 = U_u(x+v)~
    //----------------------------------------------------------
    ++(link_site[nu]) ;
    mp3->Dagger((IFloat *)GetLink(link_site, mu)+BANK4_BASE+BANK_SIZE) ;

    //----------------------------------------------------------
    // mp2 = U_u(x+v)~ U_v(x)~
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp2, (const IFloat *)mp3, (const IFloat *)mp4) ;

    //----------------------------------------------------------
    // mp4 = U_u(x+u+v)~
    //----------------------------------------------------------
    ++(link_site[mu]) ;
    mp4->Dagger((IFloat *)GetLink(link_site, mu)+BANK4_BASE);

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
    mDotMEqual((IFloat *)mp4, (const IFloat *)p1+BANK4_BASE+BANK_SIZE,
               (const IFloat *)mp3);

    //----------------------------------------------------------
    // mp2 = U_u(x+u)
    //----------------------------------------------------------
    --(link_site[mu]) ;
    moveMem((IFloat *)mp2, (const IFloat *)GetLink(link_site, mu)+BANK4_BASE,
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
    mp4->Dagger((IFloat *)GetLink(link_site, nu)+BANK4_BASE+BANK_SIZE) ;

    //----------------------------------------------------------
    // mp3 = U_u(x+u) U_v(x+2u-v)~
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp3, (const IFloat *)mp2, (const IFloat *)mp4) ;

    //----------------------------------------------------------
    // mp4 = U_u(x+u-v)~
    //----------------------------------------------------------
    --(link_site[mu]) ;
    mp4->Dagger((IFloat *)GetLink(link_site, mu)+BANK4_BASE) ;

    //----------------------------------------------------------
    // mp2 = U_u(x+u) U_v(x+2u-v)~ U_u(x+u-v)~
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp2, (const IFloat *)mp3, (const IFloat *)mp4) ;

    //----------------------------------------------------------
    // mp4 = U_u(x-v)~
    //----------------------------------------------------------
    --(link_site[mu]) ;
    mp4->Dagger((IFloat *)GetLink(link_site, mu)+BANK4_BASE+BANK_SIZE) ;

    //----------------------------------------------------------
    // mp3 = U_u(x+u) U_v(x+2u-v)~ U_u(x+u-v)~ U_u(x-v)~
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp3, (const IFloat *)mp2, (const IFloat *)mp4) ;

    //----------------------------------------------------------
    // mp2 = U_v(x-v)
    //----------------------------------------------------------
    moveMem((IFloat *)mp2, (const IFloat *)GetLink(link_site, nu)+BANK4_BASE,
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
    mDotMEqual((IFloat *)mp3, (const IFloat *)p1+BANK4_BASE+BANK_SIZE,
               (const IFloat *)mp2) ;

    //----------------------------------------------------------
    // mp4 = U_u(x-2v)~
    //----------------------------------------------------------
    mp4->Dagger((IFloat *)GetLink(link_site, mu)+BANK4_BASE) ;

    //----------------------------------------------------------
    // mp2 = U_u(x-2v)~ U_v(x-2v) U_v(x-v)
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp2, (const IFloat *)mp4, (const IFloat *)mp3) ;

    //----------------------------------------------------------
    // mp4 = U_v(x+u-2v)~
    //----------------------------------------------------------
    ++(link_site[mu]) ;
    mp4->Dagger((IFloat *)GetLink(link_site, nu)+BANK4_BASE+BANK_SIZE) ;

    //----------------------------------------------------------
    // mp3 = U_v(x+u-2v)~ U_u(x-2v)~ U_v(x-2v) U_v(x-v)
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp3, (const IFloat *)mp4, (const IFloat *)mp2) ;

    //----------------------------------------------------------
    // mp2 = U_v(x+u-v)~
    //----------------------------------------------------------
    ++(link_site[nu]) ;
    mp2->Dagger((IFloat *)GetLink(link_site, nu)+BANK4_BASE) ;

    //----------------------------------------------------------
    // rect += U_v(x+u-v)~ U_v(x+u-2v)~ U_u(x-2v)~ U_v(x-2v) U_v(x-v)
    //----------------------------------------------------------
    mDotMPlus((IFloat *)&rect, (const IFloat *)mp2, (const IFloat *)mp3) ;



    //----------------------------------------------------------
    // mp4 = U_u(x-v)~
    //----------------------------------------------------------
    --(link_site[mu]) ;
    mp4->Dagger((IFloat *)GetLink(link_site, mu)+BANK4_BASE+BANK_SIZE) ;

    //----------------------------------------------------------
    // mp3 =  U_v(x+u-v)~ U_u(x-v)~
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp3, (const IFloat *)mp2, (const IFloat *)mp4) ;

    //----------------------------------------------------------
    // mp4 = U_u(x-u-v)~
    //----------------------------------------------------------
    --(link_site[mu]) ;
    mp4->Dagger((IFloat *)GetLink(link_site, mu)+BANK4_BASE) ;

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
//             (const IFloat *)p1+BANK4_BASE+BANK_SIZE) ;
//
//  GRF: here is a work-around

    //----------------------------------------------------------
    // mp4 = U_v(x-u-v)
    //----------------------------------------------------------
    moveMem((IFloat *)mp4,
            (const IFloat *)GetLink(link_site, nu)+BANK4_BASE+BANK_SIZE,
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
    moveMem((IFloat *)mp2, (const IFloat *)GetLink(link_site, mu)+BANK4_BASE,
            MATRIX_SIZE * sizeof(IFloat)) ;

    //----------------------------------------------------------
    // rect += U_v(x+u-v)~ U_u(x-v)~ U_u(x-u-v)~ U_v(x-u-v) U_u(x-u)
    //----------------------------------------------------------
    mDotMPlus((IFloat *)&rect, (const IFloat *)mp3, (const IFloat *)mp2) ;



    //----------------------------------------------------------
    // mp4 = U_v(x-u)~
    //----------------------------------------------------------
    mp4->Dagger((IFloat *)GetLink(link_site, nu)+BANK4_BASE+BANK_SIZE) ;

    //----------------------------------------------------------
    // mp3 = U_v(x-u)~ U_u(x-u)
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp3, (const IFloat *)mp4, (const IFloat *)mp2) ;

    //----------------------------------------------------------
    // mp4 = U_u(x-u+v)~
    //----------------------------------------------------------
    ++(link_site[nu]) ;
    mp4->Dagger((IFloat *)GetLink(link_site, mu)+BANK4_BASE) ;

    //----------------------------------------------------------
    // mp2 = U_u(x-u+v)~ U_v(x-u)~ U_u(x-u)
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp2, (const IFloat *)mp4, (const IFloat *)mp3) ;

    //----------------------------------------------------------
    // mp4 = U_u(x+v)~
    //----------------------------------------------------------
    ++(link_site[mu]) ;
    mp4->Dagger((IFloat *)GetLink(link_site, mu)+BANK4_BASE+BANK_SIZE) ;

    //----------------------------------------------------------
    // mp3 = U_u(x+v)~ U_u(x-u+v)~ U_v(x-u)~ U_u(x-u)
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp3, (const IFloat *)mp4, (const IFloat *)mp2) ;

    //----------------------------------------------------------
    // mp2 = U_v(x+u)
    //----------------------------------------------------------
    ++(link_site[mu]) ;
    --(link_site[nu]) ;
    moveMem((IFloat *)mp2, (const IFloat *)GetLink(link_site, nu)+BANK4_BASE,
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
//             (const IFloat *)p1+BANK4_BASE+BANK_SIZE) ;
//
//  GRF: here is a work-around

    //----------------------------------------------------------
    // mp4 = U_v(x+u+v)
    //----------------------------------------------------------
    ++(link_site[nu]) ;
    moveMem((IFloat *)mp4,
            (const IFloat *)GetLink(link_site, nu)+BANK4_BASE+BANK_SIZE,
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
    mp4->Dagger((IFloat *)GetLink(link_site, mu)+BANK4_BASE) ;

    //----------------------------------------------------------
    // mp2 = U_v(x+u) U_v(x+u+v) U_u(x+2v)~
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp2, (const IFloat *)mp3, (const IFloat *)mp4) ;

    //----------------------------------------------------------
    // mp4 = U_v(x+v)~
    //----------------------------------------------------------
    --(link_site[nu]) ;
    mp4->Dagger((IFloat *)GetLink(link_site, nu)+BANK4_BASE+BANK_SIZE) ;

    //----------------------------------------------------------
    // mp3 = U_v(x+u) U_v(x+u+v) U_u(x+2v)~ U_v(x+v)~
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp3, (const IFloat *)mp2, (const IFloat *)mp4) ;

    //----------------------------------------------------------
    // mp2 = U_v(x)~
    //----------------------------------------------------------
    --(link_site[nu]) ;
    mp2->Dagger((IFloat *)GetLink(link_site, nu)+BANK4_BASE) ;

    //----------------------------------------------------------
    // rect += U_v(x+u) U_v(x+u+v) U_u(x+2v)~ U_v(x+v)~  U_v(x)~
    //----------------------------------------------------------
    mDotMPlus((IFloat *)&rect, (const IFloat *)mp3, (const IFloat *)mp2) ;

    //----------------------------------------------------------
    // dummy read to switch CBUF banks for looping
    //----------------------------------------------------------
    *((IFloat *)mp4) = *((IFloat *)p1+BANK4_BASE+BANK_SIZE) ;
  }
//VRB.FuncEnd(cname, fname) ;
}


//------------------------------------------------------------------
// void Lattice::Plaq(Matrix &plaq, int *x, int mu, int nu) const 
//------------------------------------------------------------------
// Purpose: 
//    calculates the plaquette U_u(x) U_v(x+u) U_u(x+v)~ U_v(x)~
// Added by Ping to help code debugging, may be more useful later on.
//------------------------------------------------------------------
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
  moveMem((IFloat *)mp2, (const IFloat *)p1+BANK4_BASE+BANK_SIZE,
          MATRIX_SIZE * sizeof(IFloat));
  mDotMEqual((IFloat *)mp3, (const IFloat *)(g_offset+mu)+BANK4_BASE,
             (const IFloat *)mp2);
  
  //----------------------------------------
  //  mp1 = (U_v(x) U_u(x+v))~
  //    p1 = &U_u(x+v) --> mp1
  //    mp2 = U_v(x) U_u(x+v)
  //    mp1 = mp2~
  //----------------------------------------
  p1 = GetLinkOld(g_offset, x, nu, mu);
  moveMem((IFloat *)mp1, (const IFloat *)p1+BANK4_BASE+BANK_SIZE,
          MATRIX_SIZE * sizeof(IFloat));
  mDotMEqual((IFloat *)mp2, (const IFloat *)(g_offset+nu)+BANK4_BASE,
             (const IFloat *)mp1);
  mp1->Dagger((IFloat *)mp2);
 
  mDotMEqual((IFloat *)(&plaq), (const IFloat *)mp3, (const IFloat *)mp1);
}  
 






//------------------------------------------------------------------
// ReTrPlaq(int *x, int mu, int nu):
// It calculates the real part of the trace of the plaquette 
// field at site x, mu, nu with mu < nu.
// The plaquette field is:
//
//   U_u(x) U_v(x+u) U_u(x+v)~ U_v(x)~
//
//------------------------------------------------------------------
Float Lattice::ReTrPlaq(int *x, int mu, int nu) const
{
//  char *fname = "ReTrPlaq(i*,i,i) const";
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
  moveMem((IFloat *)mp2, (const IFloat *)p1+BANK4_BASE+BANK_SIZE,
	  MATRIX_SIZE * sizeof(IFloat));
  mDotMEqual((IFloat *)mp3, (const IFloat *)(g_offset+mu)+BANK4_BASE,
	     (const IFloat *)mp2);


  //----------------------------------------
  //  mp1 = (U_v(x) U_u(x+v))~
  //	p1 = &U_u(x+v) --> mp1
  //	mp2 = U_v(x) U_u(x+v)
  //	mp1 = mp2~
  //----------------------------------------
  p1 = GetLinkOld(g_offset, x, nu, mu);
  moveMem((IFloat *)mp1, (const IFloat *)p1+BANK4_BASE+BANK_SIZE,
	  MATRIX_SIZE * sizeof(IFloat));
  mDotMEqual((IFloat *)mp2, (const IFloat *)(g_offset+nu)+BANK4_BASE,
	     (const IFloat *)mp1);
  mp1->Dagger((IFloat *)mp2);


  
  mDotMEqual((IFloat *)mp2, (const IFloat *)mp3, (const IFloat *)mp1);
  return mp2->ReTr();
}


//------------------------------------------------------------------
// SumReTrPlaqNode():
// It calculates the sum of the real part of the trace of the 
// plaquette field at each site of the node sublattice.
//------------------------------------------------------------------
Float Lattice::SumReTrPlaqNode(void) const
{
  char *fname = "SumReTrPlaqNode() const";
  VRB.Func(cname,fname);
  
  Float sum = 0.0;
  int x[4];
  
  for(x[0] = 0; x[0] < node_sites[0]; ++x[0]) {
    for(x[1] = 0; x[1] < node_sites[1]; ++x[1]) {
      for(x[2] = 0; x[2] < node_sites[2]; ++x[2]) {
	for(x[3] = 0; x[3] < node_sites[3]; ++x[3]) {
	  
	  for (int mu = 0; mu < 3; ++mu) {
	    for(int nu = mu+1; nu < 4; ++nu) {
	      sum += ReTrPlaq(x,mu,nu);
	    }
	  }
	}
      }
    }
  }
  
  return sum;
}


//------------------------------------------------------------------
// SumReTrPlaq():
// It calculates the sum of the real part of the trace of the 
// plaquette field at each site of the whole lattice
//------------------------------------------------------------------
Float Lattice::SumReTrPlaq(void) const
{
  char *fname = "SumReTrPlaq() const";
  VRB.Func(cname,fname);

  Float sum = SumReTrPlaqNode();
  glb_sum(&sum);
  return sum;
}


//-----------------------------------------------------------------------------
// ReTrRect(int *x, int mu, int nu):
//   It calculates the real part of the trace of the rectangle
//   field at site x, in the (mu, nu) plane with the long axis
//   of the rectangle in the mu direction.
//   The rectangle field is:
//
//   U_u(x) U_u(x+u) U_v(x+2u) U_u(x+u+v)~ U_u(x+v)~ U_v(x)~
//
//-----------------------------------------------------------------------------
Float Lattice::ReTrRect(int *x, int mu, int nu) const
{
    char *fname = "ReTrRect(i*,i,i) const";
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
    mp4->Dagger((IFloat *)GetLink(link_site, nu)+BANK4_BASE) ;

    //----------------------------------------------------------
    // mp3 = U_u(x+v)~
    //----------------------------------------------------------
    ++(link_site[nu]) ;
    mp3->Dagger((IFloat *)GetLink(link_site, mu)+BANK4_BASE+BANK_SIZE) ;

    //----------------------------------------------------------
    // mp2 = U_u(x+v)~ U_v(x)~
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp2, (const IFloat *)mp3, (const IFloat *)mp4) ;

    //----------------------------------------------------------
    // mp4 = U_u(x+u+v)~
    //----------------------------------------------------------
    ++(link_site[mu]) ;
    mp4->Dagger((IFloat *)GetLink(link_site, mu)+BANK4_BASE);

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
    mDotMEqual((IFloat *)mp2, (const IFloat *)p1+BANK4_BASE+BANK_SIZE,
               (const IFloat *)mp3);

    //----------------------------------------------------------
    // p1 = &U_u(x+u)
    //----------------------------------------------------------
    --(link_site[mu]) ;
    p1 = GetLink(link_site, mu) ;

    //----------------------------------------------------------
    // mp3 = U_u(x+u) U_v(x+2u) U_u(x+u+v)~ U_u(x+v)~ U_v(x)~
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp3, (const IFloat *)p1+BANK4_BASE,
               (const IFloat *)mp2);

    //----------------------------------------------------------
    // p1 = &U_u(x)
    //----------------------------------------------------------
    p1 = GetLink(x, mu) ;

    //----------------------------------------------------------
    // mp2 = U_u(x) U_u(x+u) U_v(x+2u) U_u(x+u+v)~ U_u(x+v)~ U_v(x)~
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp2, (const IFloat *)p1+BANK4_BASE+BANK_SIZE,
               (const IFloat *)mp3);

  }
  return mp2->ReTr();
}


//-----------------------------------------------------------------------------
// SumReTrRectNode()
//   It calculates the sum of the real part of the trace of the
//   rectangle field at each site of the node sublattice.
//-----------------------------------------------------------------------------
Float Lattice::SumReTrRectNode(void) const
{
//char *fname = "SumReTrRectNode() const";
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
// SumReTrRect()
//   It calculates the sum of the real part of the trace of the
//   rectangle field at each site of the whole lattice
//-----------------------------------------------------------------------------
Float Lattice::SumReTrRect(void) const
{
  char *fname = "SumReTrRect() const" ;
  VRB.Func(cname, fname) ;

  Float sum = SumReTrRectNode() ;
  glb_sum(&sum) ;
  return sum ;
}

//-------------------------------------------------------------------
// ReTrLoop(int *x, int *dir,int length):
//   It calculates the real part of the trace of the loop at site x
//   specified by the list of directions in dir
//   length is the length of the loop 
//
// Warning!!:
//   The user is responcible for handing in directions that close a loop!
//--------------------------------------------------------------------
Float Lattice::ReTrLoop(const int *x, const int *dir,  int length) 
{
  char *fname = "ReTrLoop(i*,i,i)";
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
// SumReTrCubeNode()
//   It calculates the sum of the real part of the trace of the
//   Cube field at each site of the node sublattice.
//
//   The Cube loop is: mu nu rho -mu -nu -rho
//-----------------------------------------------------------------------------
Float Lattice::SumReTrCubeNode(void) 
{
  char *fname = "SumReTrCubeNode()";
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
	      dir[5] = OPP_DIR(mu) ;
	      for(int nu = mu+1; nu < 4; nu++) 
		{
		   dir[1] = nu ;
		   dir[4] = OPP_DIR(nu) ;
		   for(int rho = nu+1; rho < 4; rho++)
		     {
		       dir[2] = rho ;
		       dir[3] = OPP_DIR(rho) ;
		       sum += ReTrLoop(x,dir,6);
		     }
		}
	    }
  return sum;
}


//-----------------------------------------------------------------------------
// SumReTrCube()
//   It calculates the sum of the real part of the trace of the
//   Cube field at each site of the whole lattice
//-----------------------------------------------------------------------------
Float Lattice::SumReTrCube(void)
{
  char *fname = "SumReTrCube()" ;
  VRB.Func(cname, fname) ;

  Float sum = SumReTrCubeNode() ;
  glb_sum(&sum) ;
  return sum ;
}


//------------------------------------------------------------------
// Added by Ping for anisotropic lattices
//------------------------------------------------------------------
enum {NUM_SPACE_PLAQ = 3, NUM_TIME_PLAQ = 3,
      NUM_COLORS = 3, 
      NUM_DIM = 4
};

//------------------------------------------------------------------
// Float AveReTrPlaqNodeNoXi(void) const
//------------------------------------------------------------------
// Normalization:  1 for ordered links
// Average over plaq's perpendicular to the special anisotropic dir.
// Added by Ping for anisotropic lattices
//------------------------------------------------------------------
Float Lattice::AveReTrPlaqNodeNoXi() const
{
  char *fname = "AveReTrPlaqNodeNoXi()";
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
// Average over plaq's parallel to the special anisotropic dir.
// Added by Ping for anisotropic lattices
//------------------------------------------------------------------
Float Lattice::AveReTrPlaqNodeXi() const
{
  char *fname = "AveReTrPlaqNodeXi()";
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
// Average over plaq's perpendicular to the special anisotropic dir.
// Added by Ping for anisotropic lattices
//------------------------------------------------------------------
Float Lattice::AveReTrPlaqNoXi() const
{
  char *fname = "AveReTrPlaqNoXi()";
  VRB.Func(cname,fname);

  Float sum = AveReTrPlaqNodeNoXi();
  glb_sum(&sum);
  return (sum * GJP.VolNodeSites() / GJP.VolSites());
}

//------------------------------------------------------------------
// Float AveReTrPlaqXi(void) const
//------------------------------------------------------------------
// Normalization:  1 for ordered links
// Average over plaq's perpendicular to the special anisotropic dir.
// Added by Ping for anisotropic lattices
//------------------------------------------------------------------
Float Lattice::AveReTrPlaqXi() const
{
  char *fname = "AveReTrPlaqXi()";
  VRB.Func(cname,fname);
  
  Float sum = AveReTrPlaqNodeXi();
  glb_sum(&sum);
  return (sum * GJP.VolNodeSites() / GJP.VolSites());
}


//------------------------------------------------------------------
// Lattice::MltFloat(Float factor, int dir = 3)
//------------------------------------------------------------------
// Purpose:
//   U_dir(x) => factor * U_dir(x)
// Arguments:
//   dir:    {0,1,2,3} = {x,y,z,t}.
// Added by Ping for anisotropic lattices
//------------------------------------------------------------------
// Canonical order:    U[t][z][y][x][x,y,z,t]
// Wilson order:       U[t][z][y][x][x,y,z,t] EVEN sites first!
//        Wilson spinors: checkerboarded and  ODD  sites first!
// Staggered order:    U[t][z][y][x][x,y,z,t]
//     same as canonical order up to some phase factor
// G_WILSON_HB:        U[t][z][y][x][t,x,y,z]
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
//------------------------------------------------------------------
void Lattice::EvolveGfield(Matrix *mom, Float step_size){
  char *fname = "EvolveGfield(M*,F)";
  VRB.Func(cname,fname);

  setCbufCntrlReg(4, CBUF_MODE4);

  int n_links = 4 * GJP.VolNodeSites();

  Matrix *curU_p = GaugeField();

  for(int i = 0; i < n_links; ++i) {

    // *mp1 = *(mom+i)
    moveMem((IFloat *)mp1, (IFloat *)(mom+i)+BANK4_BASE+BANK_SIZE,
    	MATRIX_SIZE * sizeof(IFloat));

    *mp1 *= step_size;	// t*iH
    *mp2 = *mp1;

    for(int j = 9; j > 1; --j) {

      // mp3 = 1 + 1/j*mp2
      oneMinusfTimesMatrix((IFloat *)mp3, invs[j],
			   (const IFloat *)mp2, 18);

      // mp2 = mp1 * mp3
      mDotMEqual((IFloat *)mp2, (const IFloat *)mp1,
		 (const IFloat *)mp3);

    }


    // mp3 = 1 + mp2
    oneMinusfTimesMatrix((IFloat *)mp3, -1, (const IFloat *)mp2, 18);
    

    // need to copy
    moveMem((IFloat *)mp2, (IFloat *)(curU_p+i)+BANK4_BASE,
	    MATRIX_SIZE * sizeof(IFloat));

    // U' = mp3 U
    mDotMEqual((IFloat *)(curU_p+i),(const IFloat *)mp3,
	       (const IFloat *)mp2);
  }
}


//------------------------------------------------------------------
// Float MomHamiltonNode(Matrix *momentum):
// The conjugate momentum Hamiltonian of the node sublattice.
//------------------------------------------------------------------
Float Lattice::MomHamiltonNode(Matrix *momentum){
  char *fname = "MomHamiltonNode(M*)";
  VRB.Func(cname,fname);

  Float ham = 0.0;
  
  int n_links = 4 * GJP.XnodeSites() * GJP.YnodeSites()
    * GJP.ZnodeSites() * GJP.TnodeSites();
  
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
void Lattice::Reunitarize(void)
{
  char *fname = "Reunitarize()";
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
}


//------------------------------------------------------------------
// void Lattice::Reunitarize(Float &dev, Float &max_diff):
// Re-unitarize the gauge field configuration
// and return:
// dev = sqrt( Sum_i [ (U(i) - V(i))^2 ] / (Vol*4*18) ),
// max_diff = Max_i[ |U(i) - V(i)| ]
// where U(i), V(i) is the gauge field before and after 
// reunitarization. The index i runs over all components of
// the gauge field.
//------------------------------------------------------------------
// Modefied by Ping on January 1999 for anisotropic lattices.
//------------------------------------------------------------------
void Lattice::Reunitarize(Float &dev, Float &max_diff)
{
  char *fname = "Reunitarize(F&,F&)";
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
}


//------------------------------------------------------------------
// int MetropolisAccept(Float delta_h):
// 0 reject, 1 accept. If delta_h < 0 it accepts unconditionally.
//------------------------------------------------------------------
int Lattice::MetropolisAccept(Float delta_h)
{
  char *fname = "MetropolisAccept(F)";
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
    } else if(delta_h < 20.0) {
      LRG.AssignGenerator(0, 0, 0, 0);
      LRG.SetInterval(1, 0);
      IFloat exp_mdh;
      IFloat rand_num;

#ifdef _TARTAN
      exp_mdh = double ( exp(- double(delta_h) ) );
#else
      exp_mdh = exp(-delta_h);
#endif
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
// void RandGaussAntiHermMatrix(Matrix *mat, Float sigma):
// It produces an anti-Hermitian matrix for each site of 
// the lattice with random entries weighted according to 
// exp(- Tr(mat^2) / (2 * sigma2)
//------------------------------------------------------------------
void Lattice::RandGaussAntiHermMatrix(Matrix *mat, Float sigma2)
{
  char *fname = "RandGaussAntiHermMatrix(M*,F)";
  VRB.Func(cname,fname);

#ifdef _TARTAN
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
		*(a+i) = LRG.Grand();
	    }
	    p->AntiHermMatrix(a);
	    p++;
	}
    }

}


//--------------------------------------------------------------------------
// void Lattice::RandGaussVector(Vector * frm, Float sigma,
//                                 FermionFieldDimension frm_dim = FIVE_D)
// Puts a Complex Three-Vector at each point of the lattice
// with random numbers weighted according to
// exp(-x^2 / (2 * sigma2)).
//--------------------------------------------------------------------------
void Lattice::RandGaussVector(Vector *frm, Float sigma2,
                          FermionFieldDimension frm_dim)
{
  RandGaussVector(frm, sigma2, 2, frm_dim);
}
//--------------------------------------------------------------------------
// void Lattice::RandGaussVector(Vector * frm, Float sigma,
//                                 FermionFieldDimension frm_dim = FIVE_D)
// Puts a Complex Three-Vector at each point of the lattice
// with random numbers weighted according to
// exp(-x^2 / (2 * sigma2)).
//--------------------------------------------------------------------------
void Lattice::RandGaussVector(Vector *frm, Float sigma2)
{
  RandGaussVector(frm, sigma2, 2, FIVE_D);
}

//--------------------------------------------------------------------------
// void Lattice::RandGaussVector(Vector * frm, int number_of_checkerboards,
//              Float sigma, FermionFieldDimension frm_dim /* = FIVE_D */)
// Puts a Complex Three-Vector at each point of the lattice
// with random numbers weighted according to
// exp(-x^2 / (2 * sigma2)).
//--------------------------------------------------------------------------
void Lattice::RandGaussVector(Vector * frm, Float sigma2, int num_chkbds,
                        FermionFieldDimension frm_dim /* = FIVE_D */)
{
  char * fname = "RandGaussVector(Vector *, Float, int, FermionFieldDimension)";
  VRB.Func(cname, fname);

  int vec_size = 2 * Colors() * SpinComponents();

  int s_node_sites = GJP.SnodeSites();
  if(frm_dim == FOUR_D || s_node_sites == 0  || Fclass() != F_CLASS_DWF) {
    s_node_sites = 1;
  }
  LRG.SetSigma(sigma2);

  IFloat * ptr = (IFloat *) frm;
  int checker, i, k, s, x[4];


  if(num_chkbds == 2) {
    for(checker = 0; checker < 2; checker++)
    for(s = 0; s < s_node_sites; s++) {
      if( (s % 2) == checker) {
        for(x[3] = 0; x[3] < GJP.TnodeSites(); x[3]++)
        for(x[2] = 0; x[2] < GJP.ZnodeSites(); x[2]++)
        for(x[1] = 0; x[1] < GJP.YnodeSites(); x[1]++)
        for(x[0] = 0; x[0] < GJP.XnodeSites(); x[0]++) {
          LRG.AssignGenerator(x);
          for(k = 0; k < vec_size; k++) {
            *(ptr++) = LRG.Grand();
          }
        }
      }
    }
  }
  else if(num_chkbds == 1) {
    for(i = 0; i < GJP.VolNodeSites(); i+=2) {
      LRG.AssignGenerator(i);
      for(k = 0; k < FsiteSize(); k++) {
        *(ptr++) = LRG.Grand();
      }
    }
  }
}


//------------------------------------------------------------------
// Sets the gauge field to the identity
//------------------------------------------------------------------
void Lattice::SetGfieldOrd(void){
  char *fname = "SetGfieldOrd()";
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
// Sets the gauge field to disordered (random) values
//------------------------------------------------------------------
void Lattice::SetGfieldDisOrd(void){
  char *fname = "SetGfieldDisOrd()";
  VRB.Func(cname,fname);

  int site_size = GsiteSize();
  LRG.SetInterval(1, -1);

  IFloat *pmat=(IFloat *) gauge_field;

  for(int i=0; i<GJP.VolNodeSites(); i++) {
    LRG.AssignGenerator(i);
    for(int k = 0; k < site_size; k++) {
      *(pmat++) = LRG.Urand();
    }
  }
  Reunitarize();

}


//------------------------------------------------------------------
// Counter related functions (g_upd_cnt, md_time)
//------------------------------------------------------------------
int Lattice::GupdCnt(void)
{
  return *g_upd_cnt ;
}

int Lattice::GupdCnt(int set_val)
{
  *g_upd_cnt = set_val ;
  return *g_upd_cnt ;
}

int Lattice::GupdCntInc(int inc_val)
{
  *g_upd_cnt += inc_val ;
  return *g_upd_cnt ;
}

Float Lattice::MdTime(void)
{
  return md_time ;
}

Float Lattice::MdTime(Float set_val)
{
  md_time = set_val ;
  return md_time ;
}

Float Lattice::MdTimeInc(Float inc_val)
{
  md_time += inc_val ;
  return md_time ;
}


//------------------------------------------------------------------
// Returns fix_gauge_ptr (pointer to an array of pointers 
// that point to the various gauge fixed hyperplanes.
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
  char *fname = "Gamma5";
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
			 int s_r, int s_l)
{
  char *fname = "Ffour2five";
  ERR.NotImplemented(cname,fname);
}

//------------------------------------------------------------------
// This is "truly" implemented only in the Fdwf derived class
//------------------------------------------------------------------
void Lattice::Ffive2four(Vector *four, Vector *five, 
			 int s_r, int s_l)
{
  char *fname = "Ffive2four";
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
  char *fname = "Fdslash";
  ERR.NotImplemented(cname,fname);
}

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
// void GsoCheck(void):
// If GJP.Snodes() == 1 it just returns.
// If GJP.Snodes() != 1 it checks that the "spread-out"
// gauge field is identical along all s-slices by comparing 
// the checksum and plaquette value. This situation arises for
// DWF with the s direction spread out across many processors.
//------------------------------------------------------------------
void Lattice::GsoCheck(void){
  char *fname = "GsoCheck()";
  VRB.Func(cname,fname);
  int s;
  IFloat rcv_buf;
  IFloat snd_buf;
  IFloat failed_value;
  Float failed_flag;
  int failed_slice;
  int s_nodes = GJP.Snodes();

  // If GJP.Snodes() != 1 do the check
  //----------------------------------------------------------------
  if(s_nodes != 1) {

    // Calculate plaquette
    //--------------------------------------------------------------
    IFloat plaq = SumReTrPlaqNode();
    
    // Compare plaquettes along fifth direction and set failed_flag
    //--------------------------------------------------------------
    failed_flag = 0.0;
    failed_slice = 0;
    snd_buf = plaq;
    for(s=1; s < s_nodes; s++){
      getMinusData(&rcv_buf, &snd_buf, 1, 4);
      if(rcv_buf != plaq) {
	failed_flag = 1.0;
#ifdef _TARTAN
	InterruptExit(-1, "GsoCheck: PLAQUETTE TEST FAILED");
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
    snd_buf = checksum;
    for(s=1; s < s_nodes; s++){
      getMinusData(&rcv_buf, &snd_buf, 1, 4);
      if(rcv_buf != checksum) {
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
// void SoCheck(Float num):
// If GJP.Snodes() == 1 it just returns.
// If GJP.Snodes() != 1 it checks that the value of num
// is identical along all s-slices. This situation arises for
// DWF with the s direction spread out across many processors.
// If the num along all s-slices is not identical it exits
// with an error.
//------------------------------------------------------------------
void Lattice::SoCheck(Float num){
  char *fname = "SoCheck()";
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





CPS_END_NAMESPACE
