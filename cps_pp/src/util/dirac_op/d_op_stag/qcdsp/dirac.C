#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-01-13 20:39:37 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_stag/qcdsp/dirac.C,v 1.3 2004-01-13 20:39:37 chulwoo Exp $
//  $Id: dirac.C,v 1.3 2004-01-13 20:39:37 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.2.4.1  2003/11/06 20:22:21  cwj
//  *** empty log message ***
//
//  Revision 1.1.1.1  2003/11/04 05:05:05  chulwoo
//
//  starting again
//
//
//  Revision 1.1.1.1  2003/06/22 13:34:46  mcneile
//  This is the cleaned up version of the Columbia Physics System.
//  The directory structure has been changed.
//  The include paths have been updated.
//
//
//  Revision 1.5  2001/08/16 10:50:21  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "IFloat".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and IFloat).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.3  2001/06/28 14:34:22  anj
//
//  The core ANSIfication should now be complete.  There are a few
//  remaining issues, but this version should compile anywhere and be
//  backward compatable with QCDSP (although this requires the top source
//  directory (.../phys/ to be added to the include path).
//
//  The serial GCC version has also been tested, and all test programs
//  appear to behave as they should (not to imply that they all work, but
//  I believe those that should work are ok).  There are minor differences
//  in the results due to rounding, (see example pbp_gccsun.dat files),
//  but that is all.
//
//  Anj.
//
//  Revision 1.2  2001/06/19 18:12:46  anj
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
//  Revision 1.2  2001/05/25 06:16:06  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: dirac.C,v $
//  $Revision: 1.3 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_stag/qcdsp/dirac.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//-------------------------------------------------------------------
//  2/6/99 RDM
//  A hopefully faster implementation of the staggered dirac operator.
//  Currently this uses twice the memory of earlier versions, since
//  the lattice is stored twice.
//
//  4/3/99 SUI 
//  Port into physics system
//-------------------------------------------------------------------

CPS_END_NAMESPACE
#include<util/gjp.h>
#include<comms/scu.h>
#include<comms/glb.h>
#include<util/lattice.h>
#include<util/dirac_op.h>
#include<util/vector.h>
#include <sysfunc.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
CPS_START_NAMESPACE

//---------------------------------------------------------------------
//  Bit 24 of chi pointers is set if the address is not absolute
//---------------------------------------------------------------------

static const int not_abs_addr = 0x1000000;

enum{VECT_LEN=6, MATRIX_SIZE=18, SITE_LEN=72};

extern unsigned int cmv_vector;

//---------------------------------------------------------------------
//  Arrays for send direction, receive direction, lattice size per
//  node, coordinate (for site where cluster of 8 matrices
//  is being assembled) and the coordinates of a nearest neighbor site.
//
//  0 == T, 1 == X, 2 == Y, 3 == Z
//--------------------------------------------------------------------

static SCUDir scudir[] =
{
  SCU_TP, SCU_XP, SCU_YP, SCU_ZP, SCU_TM, SCU_XM, SCU_YM, SCU_ZM
};

static int size[4];
static int coord[4];
static int coord_nn[4];
static int vol;

//---------------------------------------------------------------------
//  uc[0] points to a cluster of 8 matrices per even site, arranged so
//  that all parallel transport of spinors is accomplished by the
//  same matrix times vector function.  uc[1] is the same for odd sites.
//
//  At each site, the matrices are ordered as tp, xp, yp, zp, tm,
//  xm, ym, zm (p = plus, m = minus).
//---------------------------------------------------------------------

static IFloat * uc[2];


//---------------------------------------------------------------------
//  pointers to storage area for color vectors from tp, xp, yp, zp, tm,
//  xm, ym, zm (p = plus, m = minus).  Indexed as 0-7
//---------------------------------------------------------------------

static IFloat * chi_off_node[8];


//------------------------------------------------------------------
//  pointer to array of pointers telling where color vectors are
//  located for cluster arrangement of lattice (even and odd).
//
//  chi[0] points to a list of addresses for chi's needed
//  to get an even site result from application of D.  chi[0] has
//  bit 24 set if the address is not absolute.  chi[1] is similar
//  for odd sites.
//
//------------------------------------------------------------------

static IFloat ** chi[2];


//------------------------------------------------------------------
//  Values for send and receive transfers. 0-7 correspond to
//  tp, xp, yp, zp, tm, xm, ym, zm (p = plus, m = minus).
//
//  Rarg (for SCU receives) never changes, since it receives into
//  the buffers for the specified direction.
//------------------------------------------------------------------

SCUDirArg * Xarg[8];
SCUDirArg * Rarg[8];

//------------------------------------------------------------------
//  The offset into the even or odd checkboard chi's used for a send
//  from a given face.
//------------------------------------------------------------------

static int Xoffset[8];

//------------------------------------------------------------------
//  The plus and minus T faces send one word at a time.  Both faces
//  send even and odd chi's.  ToffsetP[0] is plus face, even chi's.
//  The counts are the number of transfers for these cases.
//------------------------------------------------------------------

static int * ToffsetP[2];
static int * ToffsetM[2];
static int countP[2];
static int countM[2];

//-------------------------------------------------------------------
//  Given a lexical value for gauge fields, set the coordinates.
//  Return 1 if odd, 0 if even
//-------------------------------------------------------------------

static int SetCoord( int sg );

//-------------------------------------------------------------------
//  Given an index into an eo buffer, calculate the coordinates,
//  given whether the buffer is even or odd.
//-------------------------------------------------------------------


//---------------------------------------------------------------------
//  Find nearest neighbor coordinate for coordinates given.  Nearest
//  neighbor coordinates are placed in coord_nn, which are always
//  on-node.  Function returns 0 if nearest neighbor is on-node,
//  1 if off-node.
//
//  nn = 0 to 7.  tp, xp, yp, zp, tm, xm, ym, zm respectively.
//---------------------------------------------------------------------

static int CoordNN( int nn );


//---------------------------------------------------------------------
//  Return lexical value for links from coordinates c
//---------------------------------------------------------------------

static int LexGauge( int * c );


//---------------------------------------------------------------------
//  Return lexical value for vectors from coordinates c
//---------------------------------------------------------------------

static int LexVector( int * c );


//-------------------------------------------------------------------
//  Given a coordinate and a surface ( 0 = t, 1 = x, 2 = y, 3 = x )
//  calculate the offset into the receive buffers (chi_off_node);
//-------------------------------------------------------------------

static int LexSurface( int * cc, int surface );

static IFloat * gauge_field_addr;

//-------------------------------------------------------------------
//  Called by the lattice constructor
//  Fermion initializations: pointer tables, 
//-------------------------------------------------------------------
void stag_dirac_init(const void * gauge_u )
{
  gauge_field_addr = ( IFloat * ) gauge_u;

  int c,i,j,m,n;

  int blklen[4];
  int numblk[4];
  int stride[4];

  //-------------------------------------------------------------------
  //  sg is a lexical index for (t,z,y,x) where x runs fastest. This is
  //  the gauge field order produced by convert for staggered fermions.
  //
  //    sg = x + L_x * ( y + L_y * ( z + L_z * t ) )
  //
  //  sc is a lexical index for (t,x,y,z) where t runs fastest.  The
  //  even and odd staggered color vectors are stored with indices
  //  running in this order, except that even sites come before odd
  //  sites.
  //  
  //    sc = t + L_t * ( x + L_x * ( y + L_y * z ) )
  //
  //  When the gauge fields are copied into the 8-link clusters,
  //  sg/2 indexes the clusters (both even and odd)
  //
  //  Similarly the color vectors are indexed by sc/2 for both even
  //  and odd blocks.  Even and odd blocks have a different base
  //  address.
  //-------------------------------------------------------------------

  int sg, sc;

  //-----------------------------------------------------------
  //  If t + x + y + z is odd, odd = 1.  Otherwise it is 0.
  //-----------------------------------------------------------

  int odd;

  //-----------------------------------------------------------
  //  The physics system storage order has vector indices as
  //  0-3, x,y,z,t.  Our vector indices run 0-3 as t,x,y,z.
  //  nn is used to hold physics system values for our index,
  //  given by n.
  //-----------------------------------------------------------

  size[0] = GJP.TnodeSites();
  size[1] = GJP.XnodeSites();
  size[2] = GJP.YnodeSites();
  size[3] = GJP.ZnodeSites();

  vol = size[0] * size[1] * size[2] * size[3];

  //-----------------------------------------------------------------
  //  Allocate 8 receive buffers for off-node vectors
  //-----------------------------------------------------------------

  for ( i = 0; i < 8; i++ )
    chi_off_node[i] = ( IFloat * ) smalloc(
	   VECT_LEN * vol * sizeof( IFloat ) / ( 2 * size[ i % 4 ] ) );

  //-----------------------------------------------------------------
  //  Space for storage of pointers to chi's.  9 pointers per site,
  //  but split into even and odd groups.
  //-----------------------------------------------------------------

  for ( i = 0; i < 2; i++ )
    chi[i] = ( IFloat ** ) smalloc ( 9 * vol / 2 );

  //-----------------------------------------------------------------
  //  Loop over all sites.  set up pointers to vector field
  //-----------------------------------------------------------------

  for ( sg = 0; sg < vol; sg++ ) {

    //-----------------------------------------------------------
    //  Find coordinate of site, setting odd flag.
    //-----------------------------------------------------------

    odd = SetCoord( sg );

    //-----------------------------------------------------------------
    //  Now set up chi pointers.  First do vector at site of cluster,
    //  then the 8 nearest neighbors.  0xffffff indicates that an
    //  offset must be added.
    //-----------------------------------------------------------------

    m = 9 * ( sg / 2 ) + 8;

    *( chi[ odd ]  +  m ) =
      ( IFloat * ) ( VECT_LEN * ( LexVector( coord ) / 2 ) ) +
      not_abs_addr;

    for ( n = 0; n < 8; n++ ) {

      if ( n % 2 ) c = 0xa00000;
      else c = 0xb00000;

      m = 9 * ( sg / 2 ) + n;

      if ( CoordNN( n ) ) {		//  off-node

        *( chi[ odd ]  +  m )
          = chi_off_node[n]
	  + VECT_LEN * ( LexSurface( coord_nn, n%4 ) / 2 ) + c;

      }
      else {

        *( chi[ odd ]  +  m )
          = ( IFloat * ) ( VECT_LEN * ( LexVector( coord_nn ) / 2 ) )+c
	  + not_abs_addr;
      }
    }
  }

  //-------------------------------------------------------------------
  //  Set up SCU buffer parameters.  T direction is special, since
  //  the block-strided move will not work here.
  //-------------------------------------------------------------------

  blklen[0] = VECT_LEN;
  blklen[1] = VECT_LEN * size[0] / 2;
  blklen[2] = VECT_LEN * size[0] * size[1] / 2;
  blklen[3] = VECT_LEN * size[0] * size[1] * size[2] / 2;

  numblk[0] = 1;
  numblk[1] = size[2] * size[3];
  numblk[2] = size[3];
  numblk[3] = 1;
  
  stride[0] = 1;
  stride[1] = VECT_LEN * size[0] * ( size[1] - 1 ) / 2 + 1;
  stride[2] = VECT_LEN * size[0] * size[1] * ( size[2] - 1 ) / 2 + 1;
  stride[3] = 1;

  //-------------------------------------------------------------------
  //  Index i says data has been received from TP, XP, YP, ZP, TM, XM,
  //  YM, ZM
  //-------------------------------------------------------------------

  for ( i = 0; i < 8; i++ ) {

    Rarg[i] = ( SCUDirArg * ) smalloc ( sizeof( SCUDirArg ) );
    Xarg[i] = ( SCUDirArg * ) smalloc ( sizeof( SCUDirArg ) );

    j = i % 4;

    Rarg[i]->Init( chi_off_node[i], scudir[i], SCU_REC,
      VECT_LEN * vol / ( 2 * size[j] ) );

    Xarg[i]->Init( ( void * ) 0, scudir[i], SCU_SEND, blklen[j],
      numblk[j], stride[j] );

  }

  //-------------------------------------------------------------------
  //  Need send offsets for various transfers.  The index for
  //  sends is TM, XM, YM, ZM, TP, XP, YP, ZP, since the
  //  transfers are indexed by the node data is received from.
  //-------------------------------------------------------------------

  Xoffset[0] = 0;
  Xoffset[1] = VECT_LEN * size[0] * (size[1] - 1) / 2;
  Xoffset[2] = VECT_LEN * size[0] * size[1] * (size[2] - 1) / 2;
  Xoffset[3] = VECT_LEN * size[0] * size[1] * size[2] * (size[3]-1) / 2;
  Xoffset[4] = 0;
  Xoffset[5] = 0;
  Xoffset[6] = 0;
  Xoffset[7] = 0;

  //-------------------------------------------------------------------
  //  Calculate offsets for T transfers done one word at a time.
  //  We have plus (P) transfers for both the even and odd
  //  checkerboards.  Same for minus (M) transfers.
  //-------------------------------------------------------------------


  for ( i = 0; i < 2; i++ ) {

    ToffsetP[i] = ( int * ) smalloc ( size[1] * size[2] * size[3] *
      sizeof( int ) / 2 );

    ToffsetM[i] = ( int * ) smalloc ( size[1] * size[2] * size[3] *
      sizeof( int ) / 2 );

    countP[i] = 0;
    countM[i] = 0;
  }

  for ( sg = 0; sg < vol; sg++ ) { 

    odd = SetCoord( sg );
    sc = LexVector( coord );

    if ( coord[0] == 0 ) {
      *( ToffsetM[ odd ] + countM[ odd ] ) = VECT_LEN * ( sc / 2 );
      countM[ odd ]++;
    }

    if ( coord[0] == size[0] - 1 ) {
      *( ToffsetP[ odd ] + countP[ odd ] ) = VECT_LEN * ( sc / 2 );
      countP[ odd ]++;
    }
  }
}

//---------------------------------------------------------------------
//---------------------------------------------------------------------


void stag_destroy_dirac_buf()
{
  int i;

  for ( i = 0; i < 2; i++ ) {
    sfree(ToffsetP[i]);
    sfree(ToffsetM[i]);
    sfree(chi[i]);
  }
    
  for ( i = 0; i < 8; i++ ) {
    sfree(Rarg[i]);
    sfree(Xarg[i]);
    sfree(chi_off_node[i]);
  }  
}
//-------------------------------------------------------------------
//  Given a lexical value for gauge fields, set the coordinates.
//  Return 1 if odd, 0 if even
//-------------------------------------------------------------------

static int SetCoord( int sg )
{
  coord[1] =   sg % size[1];
  coord[2] = ( sg / size[1] ) % size[2];
  coord[3] = ( sg / ( size[1] * size[2] ) ) % size[3];
  coord[0] = ( sg / ( size[1] * size[2] * size[3] ) ) % size[0];

  return ( coord[0] + coord[1] + coord[2] + coord[3] ) % 2;
}
//-------------------------------------------------------------------
//  Prepare a copy of gauge fields and set up 
//  related pointer for Dirac
//-------------------------------------------------------------------

void stag_dirac_init_g()
{

  IFloat * u = gauge_field_addr;
  int c,i,j,m,n,r;
  int off_node;

  //-------------------------------------------------------------------
  //  sg is a lexical index for (t,z,y,x) where x runs fastest. This is
  //  the gauge field order produced by convert for staggered fermions.
  //
  //    sg = x + L_x * ( y + L_y * ( z + L_z * t ) )
  //
  //  sc is a lexical index for (t,x,y,z) where t runs fastest.  The
  //  even and odd staggered color vectors are stored with indices
  //  running in this order, except that even sites come before odd
  //  sites.
  //  
  //    sc = t + L_t * ( x + L_x * ( y + L_y * z ) )
  //
  //  When the gauge fields are copied into the 8-link clusters,
  //  sg/2 indexes the clusters (both even and odd)
  //
  //  Similarly the color vectors are indexed by sc/2 for both even
  //  and odd blocks.  Even and odd blocks have a different base
  //  address.
  //-------------------------------------------------------------------

  int sg;

  //-----------------------------------------------------------
  //  If t + x + y + z is odd, odd = 1.  Otherwise it is 0.
  //-----------------------------------------------------------

  int odd;

  //-----------------------------------------------------------
  //  The physics system storage order has vector indices as
  //  0-3, x,y,z,t.  Our vector indices run 0-3 as t,x,y,z.
  //  nn is used to hold physics system values for our index,
  //  given by n.
  //-----------------------------------------------------------

  int nn;

  //-----------------------------------------------------------
  //  Once all the index arithmetic is finished, v points to
  //  the initial gauge field matrix.  w points to where it should
  //  be stored.
  //-----------------------------------------------------------

  IFloat * v;
  IFloat * w;

  //-----------------------------------------------------------
  //  SCU transfer structure to get links from off node and a
  //  location where one link matrix can be stored.
  //-----------------------------------------------------------

  SCUDirArg X((void *) 0, SCU_TP, SCU_SEND, 1);
  SCUDirArg R((void *) 0, SCU_TP, SCU_REC, 1);
  IFloat mtmp[18];

  size[0] = GJP.TnodeSites();
  size[1] = GJP.XnodeSites();
  size[2] = GJP.YnodeSites();
  size[3] = GJP.ZnodeSites();

  vol = size[0] * size[1] * size[2] * size[3];

  //-----------------------------------------------------------------
  //  Allocate space for two copies of the gauge fields on this node
  //-----------------------------------------------------------------

  for ( i = 0; i < 2; i++ )
    uc[i] = ( IFloat * ) smalloc( SITE_LEN * vol * sizeof(IFloat) ); 

  //-----------------------------------------------------------------
  //  Loop over all sites.  First rearrange gauge field for this
  //  site and then set up pointers to vector field
  //-----------------------------------------------------------------

  for ( sg = 0; sg < vol; sg++ ) {

    //-----------------------------------------------------------
    //  Find coordinate of site, setting odd flag.
    //-----------------------------------------------------------

    odd = SetCoord( sg );

    //-----------------------------------------------------------
    //  Copy links in positive directions.  These are all on-node
    //  Take Hermitian conjugate while doing this
    //-----------------------------------------------------------

    for ( n = 0; n < 4; n++ ) {

      nn = ( n - 1 + 4 ) % 4;		//  nn is physics system order

      v = u + SITE_LEN * sg + MATRIX_SIZE * nn;

      w = uc[odd] + 2 * SITE_LEN * ( sg / 2 ) + 6 * n;

      for ( r = 0; r < 3; r++ ) {
        for ( c = 0; c < 3; c++ ) {

          i = 48*r + 2*c;

          j = 6*c + 2*r;

          w[i] = v[j];
          w[i+1] = -v[j+1];
        }
      }

    }


    //-----------------------------------------------------------------
    //  Copy links in negative directions.  Some are off-node
    //  Add an overall minus sign
    //-----------------------------------------------------------------

    for ( n = 0; n < 4; n++ ) {

      nn = ( n - 1 + 4 ) % 4;		//  nn is physics system order

      w = uc[odd] + 2*SITE_LEN * (sg/2) + 24 + 6 * n;

      //---------------------------------------------------------------
      //  Find index for link in -n direction
      //---------------------------------------------------------------

      off_node = CoordNN( n + 4 );

      v = u + SITE_LEN * LexGauge( coord_nn ) + MATRIX_SIZE * nn;

      //-----------------------------------------------------------
      //  If at edge of lattice, get link from neighbor
      //-----------------------------------------------------------

      if ( off_node ) {

        X.Init( v, scudir[n], SCU_SEND, 18 );
        R.Init( mtmp, scudir[n+4], SCU_REC, 18 );
        SCUTrans( &X );
        SCUTrans( &R );
        SCUTransComplete();

        v = mtmp;
      }

      for ( i = 0; i < 18; i++ ) {

	m = i % 6 + 48 * ( i / 6 );

        w[m] = -v[i];
      }

    }
  }
}

void stag_destroy_dirac_buf_g(void)
{
  sfree(uc[0]);
  sfree(uc[1]);
}


//---------------------------------------------------------------------
//  Find nearest neighbor coordinate for coordinates given.  Nearest
//  neighbor coordinates are placed in coord_nn, which are always
//  on-node.  Function returns 0 if nearest neighbor is on-node,
//  1 if off-node.
//
//  nn = 0 to 7.  tp, xp, yp, zp, tm, xm, ym, zm respectively.
//---------------------------------------------------------------------

static int CoordNN( int nn )
{
  int i;
  int off_node = 0;
  int m;

  //-----------------------------------------------------------
  //  Find index for link in -n direction
  //-----------------------------------------------------------

  for ( i = 0; i < 4; i++ )
    coord_nn[i] = coord[i];

  m = nn % 4;

  if ( ( nn / 4 == 0 && coord_nn[ m ] == size[ m ] - 1 ) ||
       ( nn / 4 == 1 && coord_nn[ m ] == 0             ) ) off_node = 1;

  coord_nn [ m ] = ( coord_nn [ m ] + 1 - 2 * ( nn / 4 ) + size[m] )
		 % size[m]; 

  return off_node;
}

//---------------------------------------------------------------------
//  Return lexical value for links from coordinates c
//---------------------------------------------------------------------

static int LexGauge( int * c )
{
  return c[1] + size[1] * ( c[2] + size[2] * ( c[3] + size[3] * c[0] ));
}

//---------------------------------------------------------------------
//  Return lexical value for vectors from coordinates c
//---------------------------------------------------------------------

static int LexVector( int * c )
{
  return c[0] + size[0] * ( c[1] + size[1] * ( c[2] + size[2] * c[3] ));
}

//-------------------------------------------------------------------
//  Given a coordinate and a surface ( 0 = t, 1 = x, 2 = y, 3 = x )
//  calculate the offset into the receive buffers (s?) .
//-------------------------------------------------------------------

static int LexSurface( int * cc, int surface )
{
 
  int i;
  int s[4];
  int c[4];

  for ( i = 0; i < 4; i++ ) {
    s[i] = size[i];
    c[i] = cc[i];
  }

  s[ surface ] = 1;
  c[ surface ] = 0;

  return c[0] + s[0] * ( c[1] + s[1] * ( c[2] + s[2] * c[3] ));
}

//-------------------------------------------------------------------
//  add_flag = 0 :       b = D a
//  add_flag = 1 :       b += D a
//
//  a_odd    = 1 :	 b even;  a odd
//  a_odd    = 0 :	 b odd ;  b even
//
//   D a = \sum_u( U^dag_u(x) a(x+u) - U(x-u) a(x-u) )
//-------------------------------------------------------------------

// These two external functions are defined in dirac_serial.asm
//
extern void dirac_SCU( SCUDirArg ** Xarg, SCUDirArg ** Rarg, IFloat * a,
  int a_odd, int * Xoffset, int ** ToffsetP, int ** ToffsetM,
  int * countP, int * countM );

extern void dirac_cmv( int sites, IFloat ** chi,
  IFloat * u, int a, int b, int add_flag, IFloat * cmv_vector );

extern "C"	// to comply with the dirac interface by Dong Chen
void stag_dirac(IFloat* b, IFloat* a, int a_odd, int add_flag)
{
  int i;
  int odd = 1 - a_odd;

  //-----------------------------------------------------------------
  //  Transfer chi's on faces.  First Send plus and receive minus
  //-----------------------------------------------------------------

  dirac_SCU( Xarg, Rarg, a, a_odd, Xoffset, ToffsetP, ToffsetM, countP,
    countM );

  unsigned int * cbuf = ( unsigned int * ) 0x815801;

  for ( i = 0; i < 4; i++ )
    cbuf[i] = 0xc2088100;

  dirac_cmv( 9 * (vol/2), chi[odd], uc[odd] + 0x900000, (int) a,
     (int) b, add_flag, ( IFloat * ) &cmv_vector );
}

CPS_END_NAMESPACE
