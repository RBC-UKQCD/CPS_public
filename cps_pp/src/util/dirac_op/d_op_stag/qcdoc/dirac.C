/*!\file
  Staggered Dirac operator for QCDOC

  $Id: dirac.C,v 1.17 2010/08/16 20:56:18 chulwoo Exp $
*/
//-------------------------------------------------------------------
//   12/27/01 Calin Cristian
//   Staggered Dirac operator for QCDOC. Communications and computations
//   are overlapped.
//   Uses many functions already implemented by RDM for QCDSP 
//
//-------------------------------------------------------------------

#include <util/gjp.h>
#include <util/qcdio.h>
#include <comms/scu.h>
#include <comms/glb.h>
#include <util/lattice.h>
#include <util/dirac_op.h>
#include <util/vector.h>
#include <stdio.h>
#include <comms/sysfunc_cps.h>
#include <qalloc.h>

#undef CPP

CPS_START_NAMESPACE

const char *chi_l_filename = CWDPREFIX("chi_l.h");
const char *chi_nl_filename = CWDPREFIX("chi_nl.h");

void dirac_cmv_l( int sites, long chi, long u, long a, long tmpfrm);
void dirac_cmv_nl( int sites, long chi, long u, long a, long tmpfrm);

extern "C" void save_reg(long intbuf, long dbuf );
extern "C" void restore_reg(long intbuf, long dbuf );
extern "C" void copy_buffer(int n, long src, long dest, long ptable);
extern "C" void dirac_cmv( int sites, long chi, long u, 
                       long a, long tmpfrm);
extern "C" void dirac_sum_acc( int sites, long chi, long tmpfrm, 
                       long b);
extern "C" void dirac_sum( int sites, long chi, long tmpfrm, 
                       long b);
//flushes nflush * 192 bytes + 2 cache lines
extern "C" void flush_cache_spinor(int nflush, long flush_buffer);
//flushes nflush * 128 bytes + 2 cache lines
extern "C" void flush_cache(int nflush, long flush_buffer);

enum{VECT_LEN=6, MATRIX_SIZE=18, SITE_LEN=72, NUM_DIR=8};
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
static int non_local_chi;
static int local_chi;
//static int buffer_flush[8];
static int nflush;

//---------------------------------------------------------------------
//  uc_l[0] points to a cluster of matrices per even site for local 
//  computations , arranged so that all parallel transport of spinors 
//  is accomplished by the same matrix times vector function.  
//  uc_l[1] is the same for odd sites.
//  uc_nl[0] is the same for even non-local sites.
//  uc_nl[1] is the same for odd non-local sites.
//---------------------------------------------------------------------

static IFloat * uc_l[2];
static IFloat * uc_nl[2];

//------------------------------------------------------------------
// Allocate these arrays dynamicaly once the cache, noncached eDRAM
// malloc are available (should be changed according to volume)
//------------------------------------------------------------------
static IFloat *tmpfrm;
const int MAX_TBUF_LEN = 1024;
#if 1
static IFloat *Tbuffer[2];
#else
static IFloat Tbuffer[2][MAX_TBUF_LEN] PEC_ALIGN LOCATE("Edramnoncache");
#endif
//-------------------------------------------------------------------
// end of stack based arrays which should be heap based 
//-------------------------------------------------------------------
//static int intreg[18];
//static IFloat dreg[18];
static int * ToffsetP[2];
static int * ToffsetM[2];
static int countP[2];
static int countM[2];

//---------------------------------------------------------------------
//  pointers to storage area for color vectors from tp, xp, yp, zp, tm,
//  xm, ym, zm (p = plus, m = minus).  Indexed as 0-7
//---------------------------------------------------------------------
#if 1
static IFloat * chi_off_node[8];
#else
static IFloat  chi_off_node[8][MAX_TBUF_LEN] PEC_ALIGN LOCATE("Edramnoncache");
#endif


//------------------------------------------------------------------
//  pointer to array of pointers telling where color vectors are
//  located for cluster arrangement of lattice (even and odd).
//  chi_l[0] points to a list of local adresses for chi's needed 
//  to get an even site result from application of D and 
//  to a temporary storage area where U_mu * chi's for each direction 
//  are stored. 
//  chi_l[1] same for odd site. 
//  chi_nl[0] same for computations involving non-local spinors
//  chi[0] has 9 pointers per even site, the first is the address where
//  the result of the  addition of a cluster of 8 color vectors 
//  (specified by the following 8 pointers) is stored
//  chi[1] is the same for odd sites 
//------------------------------------------------------------------

static IFloat ** chi_l[2];
static IFloat ** chi_nl[2];
static IFloat ** chi[2];
//------------------------------------------------------------------
//  Values for send and receive transfers. 0-7 correspond to
//  tp, xp, yp, zp, tm, xm, ym, zm (p = plus, m = minus).
//
//  Rarg (for SCU receives) never changes, since it receives into
//  the buffers for the specified direction.
//------------------------------------------------------------------

static SCUDirArgIR * SCUarg[16];
static SCUDirArgMulti * SCUmulti;

//------------------------------------------------------------------
//  The offset into the even or odd checkboard chi's used for a send
//  from a given face.
//------------------------------------------------------------------

static int Xoffset[NUM_DIR];

//-------------------------------------------------------------------
//  Given a lexical value for gauge fields, set the coordinates.
//  Return 1 if odd, 0 if even
//-------------------------------------------------------------------

static int SetCoord( int sg );

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
//  Fermion initializations: pointer tables 
//-------------------------------------------------------------------
static int initted = 0;
extern "C" void stag_dirac_init(const void * gauge_u )
{
  gauge_field_addr = ( IFloat * ) gauge_u;
  int i,j,m,n;
  int blklen[NUM_DIR/2];
  int numblk[NUM_DIR/2];
  int stride[NUM_DIR/2];
  int local_count[2];
  int non_local_count[2];

  int x[NUM_DIR/2];
  char *cname = "";
  char *fname = "stag_dirac_init(const void *gauge)";
  if (initted !=0) {
    Fprintf(stderr,"stag_dirac_init already initted\n");
    return;
  }
  VRB.Func(cname,fname);
  initted = 1;

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
  VRB.Result(cname,fname,"vol=%d\n",vol);
  non_local_chi = 2*(size[0]*size[1]*size[2] + size[1]*size[2]*size[3]+
    size[2]*size[3]*size[0] + size[3]*size[0]*size[1]);
  local_chi = NUM_DIR*vol - non_local_chi;

  //-------------------------------------------------------------
  // flush_cache_spinor() function will flush 192 bytes * nflush 
  //-------------------------------------------------------------
  nflush = vol/8;

#if 0
  if (vol>16000)
  tmpfrm = (IFloat *) smalloc ( 8 * vol/2 * VECT_LEN * sizeof(IFloat),
				cname,fname, "tmpfrm");
  else
  tmpfrm = (IFloat *) fmalloc ( 8 * vol/2 * VECT_LEN * sizeof(IFloat),
				cname,fname, "tmpfrm");
#endif

  


  //-----------------------------------------------------------------
  //  Allocate 8 receive buffers for off-node vectors
  //-----------------------------------------------------------------
  
  for ( i = 0; i < NUM_DIR; i++ ){
#if 1
    chi_off_node[i] = ( IFloat * ) fmalloc(cname,fname,"chi_off_node[i]",
 VECT_LEN * vol * sizeof( IFloat ) / ( 2 * size[ i % 4 ] ) );    
    if(chi_off_node[i] == 0)
      ERR.Pointer(cname,fname, "chi_off_node[i]");
#else
    if( (vol/size[i%4])*VECT_LEN/2 >MAX_TBUF_LEN ){
	ERR.General(cname,fname,"chi_off_node size overflow\n");
    }
#endif
  }

  //-----------------------------------------------------------------
  //  Space for storage of pointers to chi's.  2 pointers per site,
  //  but split into even and odd groups for the first part of the 
  //  computation (parallel transport of spinors). 9 pointers per site
  //  to obtain the result of the application of the dirac operator
  //-----------------------------------------------------------------
  

  for ( i = 0; i < 2; i++ ){
      VRB.Result(cname,fname,"local_chi=%d sizeof(IFloat)=%d\n",local_chi,
		 sizeof(IFloat));
      chi[i] = (IFloat **) fmalloc(9 * vol/2 * sizeof(IFloat *),
				   cname,fname, "chi[i]");
      chi_l[i] = ( IFloat ** ) fmalloc(2*(local_chi/2)*sizeof(IFloat *),
				       cname,fname, "chi_l[i]");
      chi_nl[i] = (IFloat ** ) fmalloc(2*(non_local_chi/2)*sizeof(IFloat *),
				       cname,fname, "chi_nl[i]");
  }
  
  for ( i = 0; i < 2; i++){
    local_count[i] = 0;
    non_local_count[i] = 0;
  }
  //-----------------------------------------------------------------
  // Assembly written for double precision only, check sizeof(IFloat)
  //-----------------------------------------------------------------
  if ( sizeof(IFloat) != sizeof(double)){
     ERR.General(cname, fname, 
		 "Assembly functions implemented only for double precision!");
  }

  //-----------------------------------------------------------------
  //  Loop over all directions
  //-----------------------------------------------------------------
  for ( n = 0; n < NUM_DIR; n++ ) {
    //-----------------------------------------------------------------
    //  Loop over all sites
    //-----------------------------------------------------------------
    
    for (x[3] = 0; x[3] < size[3]; x[3]++){
      for (x[2] = 0; x[2] < size[2]; x[2]++){
	for (x[1] = 0; x[1] < size[1]; x[1]++){
	  for (x[0] = 0; x[0] < size[0]; x[0]++){  

	    for (i = 0; i < 4 ; i++) coord[i] = x[i];
	    odd = ( coord[0] + coord[1] + coord[2] + coord[3] ) % 2;
	    sg = coord[1] + size[1] * ( coord[2] + size[2] * ( coord[3] +
			       			       size[3] * coord[0] ));
	    m = (NUM_DIR + 1) * (sg/2);
	    if ( CoordNN( n ) ) {		//  off-node
	      //----------------------------------------------------------
	      // Assembly written for double precision only, multiplication
	      // by sizeof(double) done to avoid a bitshift inside the
	      // high performance code
	      //----------------------------------------------------------
	      //pointer to source field (offset in the receive buffer)
	      *( chi_nl[ odd ]  +  2 * non_local_count[ odd ] )
	    = chi_off_node[n] + VECT_LEN * ( LexSurface( coord_nn, n%4 ) / 2 );
	      // pointer to temporary field where U*chi is stored
		*( chi_nl[ odd ] + 2 * non_local_count[ odd ] + 1) = 
		( IFloat *) ( VECT_LEN * (NUM_DIR * int(sg/2) + n )
		              * sizeof(IFloat));
	      // pointer to the above temporary field 
	      *( chi[ odd ] + m + n + 1) = 
		( IFloat *) ( VECT_LEN * (NUM_DIR * int(sg/2) + n) 
		              * sizeof(IFloat));
	      // Pointer to solution field
	      *( chi[ odd ]  +  m ) = 
		( IFloat * ) ( VECT_LEN * (LexVector( coord ) / 2 ) 
			       * sizeof(IFloat));  
	      non_local_count[odd]++; 
	    }
	    else{//on node
	      //pointer to source field
	      *( chi_l[ odd ]  +  2 * local_count[ odd ] )
		= ( IFloat * ) ( VECT_LEN * ( LexVector( coord_nn ) / 2 )
		                 * sizeof(IFloat));
	      // pointer to temporary field where U*chi is stored
	      *( chi_l[ odd ]  +  2 * local_count[ odd ] + 1) = 
		( IFloat * ) ( VECT_LEN * (NUM_DIR * int(sg/2) + n)
		               * sizeof(IFloat));
	      // pointer to the above temporary field
	      *( chi[ odd ] + m + n + 1) = 
		( IFloat *) ( VECT_LEN * (NUM_DIR * int(sg/2) + n)
		              * sizeof(IFloat));
	      // pointer to solution field
	      *( chi[ odd ]  +  m ) = 
		( IFloat * ) ( VECT_LEN * (LexVector( coord ) / 2 ) 
		               * sizeof(IFloat));
	      local_count[odd]++; 
	    }
	  }
	}
      }
    }
  }

#if 0
  char buf[200];

  sprintf(buf,"chi.h");
  int fd = open(buf,O_CREAT|O_TRUNC|O_RDWR,00644);

  for(j=0;j<2;j++){
    sprintf(buf,"IFloat * chi%d[] LOCATE(\"edramtransient\") = {\n",j); 
    write(fd,buf,strlen(buf));
    sprintf(buf," (IFloat *) %d",*(chi[j])); 
    write(fd,buf,strlen(buf));
    for(i=1;i< 9*vol/2;i++){
      sprintf(buf,",\n (IFloat *) %d",*(chi[j]+i)); 
      write(fd,buf,strlen(buf));
    }
    sprintf(buf,"\n};\n"); 
    write(fd,buf,strlen(buf));
  }
  close(fd);
#endif

#if 0
  char filename[200];
  sprintf(filename,"%s_%d%d%d%d%d%d",
  chi_l_filename,CoorX(), CoorY(), CoorZ(), CoorT(), CoorS(), CoorW());
  FILE *fp = Fopen(filename,"w");
  for(j=0;j<2;j++){
    Fprintf(fp,"IFloat * chi_l%d[] LOCATE(\"edramtransient\") = {\n",j); 
    Fprintf(fp," (IFloat *) %d",*(chi_l[j])); 
    for(i=1;i< 2*(local_chi/2);i++){
      Fprintf(fp,",\n (IFloat *) %d",*(chi_l[j]+i)); 
    }
    Fprintf(fp,"\n};\n"); 
  }
  Fclose(fp);
#endif

#if 0
  char filename[200];
  sprintf(filename,"%s_%d%d%d%d%d%d",
  chi_nl_filename, CoorX(), CoorY(), CoorZ(), CoorT(), CoorS(), CoorW());
  FILE *fp = Fopen(filename,"w");
  for(j=0;j<2;j++){
    Fprintf(fp,"IFloat * chi_nl%d[] LOCATE(\"edramtransient\") = {\n",j); 
    Fprintf(fp," (IFloat *) 0x%x",*(chi_nl[j])); 
    for(i=1;i< 2*(non_local_chi/2);i++){
      Fprintf(fp,",\n (IFloat *) 0x%x",*(chi_nl[j]+i)); 
    }
    Fprintf(fp,"\n};\n"); 
  }

  Fclose(fp);

#endif 

  //-------------------------------------------------------------------
  //  Set up SCU buffer parameters.  T direction is special, since
  //  the block-strided move will not work here.
  //-------------------------------------------------------------------

  blklen[0] = VECT_LEN * sizeof(IFloat) * size[1] * size[2] * size[3] / 2;
  blklen[1] = VECT_LEN * sizeof(IFloat) * size[0] / 2;
  blklen[2] = VECT_LEN * sizeof(IFloat) * size[0] * size[1] / 2;
  blklen[3] = VECT_LEN * sizeof(IFloat) * size[0] * size[1] * size[2] / 2;
  

  numblk[0] = 1; 
  numblk[1] = size[2] * size[3];
  numblk[2] = size[3];
  numblk[3] = 1;
  
  stride[0] = 0;
  stride[1] = (VECT_LEN * size[0] * ( size[1] - 1 ) / 2)*sizeof(IFloat);
  stride[2] = (VECT_LEN * size[0] * size[1] * ( size[2] - 1 ) / 2)*
    sizeof(IFloat) ;
  stride[3] = 0;
  //-------------------------------------------------------------------
  //  Calculate offsets for T transfers done one word at a time.
  //  We have plus (P) transfers for both the even and odd
  //  checkerboards.  Same for minus (M) transfers.
  //-------------------------------------------------------------------


  for ( i = 0; i < 2; i++ ) {
#if 1
    Tbuffer[i] = (IFloat *) qalloc (QFAST|QNONCACHE, size[1] * size[2] * size[3] *
				    VECT_LEN * sizeof( IFloat ) / 2);
    if(!Tbuffer) ERR.Pointer(cname, fname, "Tbuffer");
#else
    if( size[1]*size[2]*size[3]*VECT_LEN/2 >MAX_TBUF_LEN ){
	ERR.General(cname,fname,"Tbuffer size overflow\n");
    }
#endif
    ToffsetP[i] = ( int * ) fmalloc ( size[1] * size[2] * size[3] *
      sizeof( int ) / 2 );

    ToffsetM[i] = ( int * ) fmalloc ( size[1] * size[2] * size[3] *
      sizeof( int ) / 2 );

    countP[i] = 0;
    countM[i] = 0;
  }

//  printf("dirac_init: Set up SCU parameters\n");
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
//  printf("dirac_init: Set up SCU parameters\n");

  //-------------------------------------------------------------------
  //  Index i says data has been received from TP, XP, YP, ZP, TM, XM,
  //  YM, ZM
  //-------------------------------------------------------------------

//  for(i=0;i<4;i++)
//   printf("blklen numblk stride [%d]= %d %d %d\n",i, blklen[i],numblk[i],stride[i]);

  for ( i = 0; i < NUM_DIR; i++ ) {
    j = i % (NUM_DIR/2);
//      SCUarg[i + 8] = new SCUDirArgIR;
//      printf("%d: %p %d\n",i+8,chi_off_node[i],blklen[j]*numblk[j]);
      SCUarg[i + 8]  = new SCUDirArgIR(chi_off_node[i], scudir[i], SCU_REC, 
		    blklen[j]*numblk[j], 1, 0, IR_5);
//		    VECT_LEN * sizeof(IFloat) * vol / ( 2 * size[j] ), 
//			       1, 0, IR_5);
//      buffer_flush[i] = VECT_LEN * sizeof(IFloat) * vol/ (384 * size[j]);
//send arguments
//   SCUarg[i+8]->Print();
    if ((i == 0) || ( i == 4)){
      SCUarg[i] = new SCUDirArgIR(Tbuffer[(4 - i)/4], scudir[i], SCU_SEND, 
		       blklen[j], numblk[j], stride[j], IR_5 );
    }
    else{ 
      SCUarg[i] = new SCUDirArgIR(Tbuffer[0], scudir[i], SCU_SEND, 
		       blklen[j], numblk[j], stride[j], IR_5 );
    }
//   SCUarg[i]->Print();
//    printf("SCUarg[%d] done\n",i);
  }
//  for(i = 0;i<2*NUM_DIR;i++) SCUarg[i]->Print();
  SCUmulti = new SCUDirArgMulti();
  SCUmulti->Init(SCUarg, 16);
//  for(i = 0;i<2*NUM_DIR;i++) SCUarg[i]->Print();
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

//  print("dirac_init: Done\n");

}


extern "C" void stag_destroy_dirac_buf()
{
  int i;
//  ffree(tmpfrm);
  delete SCUmulti;
  
  for ( i = 0; i < 2; i++ ) {
    qfree(Tbuffer[i]);
    ffree(ToffsetP[i]);
    ffree(ToffsetM[i]);
    ffree(chi_nl[i]);
    ffree(chi[i]);
    ffree(chi_l[i]);
  }
    
  for ( i = 0; i < NUM_DIR; i++ ) {
    delete SCUarg[i];
    delete SCUarg[i+8];
    ffree(chi_off_node[i]);
  }  
  initted=0;
}
//-------------------------------------------------------------------
//  Given a lexical value for gauge fields, set the coordinates.
//  Return 1 if odd, 0 if even
//-------------------------------------------------------------------

static int SetCoord( int sg )
{
  coord[1] =   sg % size[1];
  coord[2] = ((int)( sg / size[1] )) % size[2];
  coord[3] = ((int)( sg / ( size[1] * size[2] ) )) % size[3];
  coord[0] = ((int)( sg / ( size[1] * size[2] * size[3] ) )) % size[0];

  return ( coord[0] + coord[1] + coord[2] + coord[3] ) % 2;
}

//-------------------------------------------------------------------
//  Prepare a copy of gauge fields and set up 
//  related pointer for Dirac
//-------------------------------------------------------------------

extern "C" void stag_dirac_init_g()
{

  IFloat * u = gauge_field_addr;
  int c,i,j,n,r;
  int off_node;
  int local_count[2];
  int nflush_g = 1;
  int non_local_count[2];  
  int x[NUM_DIR/2];
  char *cname = "DiracOpStag";
  char *fname = "dirac_init_g()";

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

#if 0
  SCUDirArg X;
  SCUDirArg R;
#endif

#if 1
  if (vol>16000)
  tmpfrm = (IFloat *) smalloc ( 8 * vol/2 * VECT_LEN * sizeof(IFloat),
				cname,fname, "tmpfrm");
  else
  tmpfrm = (IFloat *) fmalloc ( 8 * vol/2 * VECT_LEN * sizeof(IFloat),
				cname,fname, "tmpfrm");
#endif
  IFloat *mtmp;
  mtmp = (IFloat *)smalloc(18*sizeof(IFloat));
//  printf("mtmp=%p\n",mtmp);

  VRB.Func(cname,fname);
//  print("dirac_init_g start\n");
  size[0] = GJP.TnodeSites();
  size[1] = GJP.XnodeSites();
  size[2] = GJP.YnodeSites();
  size[3] = GJP.ZnodeSites();

  vol = size[0] * size[1] * size[2] * size[3];
  non_local_chi = 2*(size[0]*size[1]*size[2] + size[1]*size[2]*size[3]+
    size[2]*size[3]*size[0] + size[3]*size[0]*size[1]);
  local_chi = NUM_DIR*vol - non_local_chi;

  //-----------------------------------------------------------------
  //  Allocate space for two copies of the gauge fields on this node
  //-----------------------------------------------------------------
  for ( i = 0; i < 2; i++ ){
      uc_l[i] = (IFloat*) fmalloc(MATRIX_SIZE* local_chi/2 * sizeof(IFloat),
				  cname,fname, "uc_l[i]");
      uc_nl[i] = (IFloat*) fmalloc( cname,fname, "uc_nl[i]", 
       MATRIX_SIZE* non_local_chi/2 * sizeof(IFloat));
  }

  for ( i = 0; i < 2; i++){
    local_count[i] = 0;
    non_local_count[i] = 0;
  }
  
  //-----------------------------------------------------------
  //  Copy links in positive directions.  These are all on-node
  //  Take Hermitian conjugate while doing this
  //-----------------------------------------------------------
  
  for ( n = 0; n < NUM_DIR/2; n++ ) {
 
  //-----------------------------------------------------------------
  //  Loop over all sites.  First rearrange gauge field for this
  //  site and then set up pointers to vector field
  //-----------------------------------------------------------------

    for (x[3] = 0; x[3] < size[3]; x[3]++){
      for (x[2] = 0; x[2] < size[2]; x[2]++){
	for (x[1] = 0; x[1] < size[1]; x[1]++){
	  for (x[0] = 0; x[0] < size[0]; x[0]++){
	    for (i = 0; i < 4 ; i++) coord[i] = x[i];
	    odd = ( coord[0] + coord[1] + coord[2] + coord[3] ) % 2;
	    sg = coord[1] + size[1] * ( coord[2] + size[2] * ( coord[3] +
						   size[3] * coord[0] ));
	    nn = ( n - 1 + 4 ) % 4;		//  nn is physics system order
	    v = u + SITE_LEN * sg + MATRIX_SIZE * nn;

	    if ( CoordNN( n ) ) {		//  off-node

	      w = uc_nl[odd] + MATRIX_SIZE * non_local_count[odd];
	      for ( r = 0; r < 3; r++ ) {
		for ( c = 0; c < 3; c++ ) {

		  i = 6*r + 2*c;

		  j = 6*c + 2*r;

		  w[i] = v[j];
		  w[i+1] = -v[j+1];
		}
	      }

	      non_local_count[odd]++;
	    }
	    else{ //on-node
	      w = uc_l[odd] + MATRIX_SIZE * local_count[odd];
	      for ( r = 0; r < 3; r++ ) {
		for ( c = 0; c < 3; c++ ) {
		  
		  i = 6*r + 2*c;
		  
		  j = 6*c + 2*r;
		  
		  w[i] = v[j];
		  w[i+1] = -v[j+1];
		}
	      }
	      local_count[odd]++;
	    }
	  }
	}
      }
    }
  }

    //-----------------------------------------------------------------
    //  Copy links in negative directions.  Some are off-node
    //  Add an overall minus sign, setting dummy address to be changed
    //  via Addr() later.
    //-----------------------------------------------------------------
    for ( n = 0; n < NUM_DIR/2; n++ ) {
  SCUDirArgIR X( uc_l[0], scudir[n], SCU_SEND, 18 * sizeof(IFloat));
  SCUDirArgIR R( uc_l[0], scudir[n+4], SCU_REC, 18 * sizeof(IFloat));

      for (x[3] = 0; x[3] < size[3]; x[3]++){
	for (x[2] = 0; x[2] < size[2]; x[2]++){
	  for (x[1] = 0; x[1] < size[1]; x[1]++){
	    for (x[0] = 0; x[0] < size[0]; x[0]++){

	      for (i = 0; i < 4 ; i++) coord[i] = x[i];
	      odd = ( coord[0] + coord[1] + coord[2] + coord[3] ) % 2;
	      sg = coord[1] + size[1] * ( coord[2] + size[2] * ( coord[3] +
						   size[3] * coord[0] ));
	      nn = ( n - 1 + 4 ) % 4;		//  nn is physics system order

	      //---------------------------------------------------------------
	      //  Find index for link in -n direction
	      //---------------------------------------------------------------

	      off_node = CoordNN( n + 4 );
	      
	      v = u + SITE_LEN * LexGauge( coord_nn ) + MATRIX_SIZE * nn;


	      //-----------------------------------------------------------
	      //  If at edge of lattice, get link from neighbor
	      //-----------------------------------------------------------

	      if ( off_node ) {
                //make sure link is in main memory
                flush_cache(nflush_g, (long)v);
                //invalidate receive buffer
                flush_cache(nflush_g, (long)mtmp);
                //initialize transfer
  X.Addr(v);
  R.Addr(mtmp);
		X.SlowStartTrans();
		R.SlowStartTrans();
		X.TransComplete();
		R.TransComplete();
		v = mtmp;
	      }

	      if ( off_node ) {		//  off-node
		w = uc_nl[odd] + MATRIX_SIZE * non_local_count[odd];
		
		for ( i = 0; i < 18; i++ ) {
		  
		  w[i] = -v[i];
		}
		non_local_count[odd]++;
	      }
	      else{
		w = uc_l[odd] + MATRIX_SIZE * local_count[odd];

		for ( i = 0; i < 18; i++ ) {
		  
		  w[i] = -v[i];

		}
		local_count[odd]++;
	      }
	    }
	  }
	}
      }
    }
    ffree(mtmp);
//    printf("Copying links done\n");

#if 0
  char buf[200];

  sprintf(buf,"uc_l.h");
  int fd = open(buf,O_CREAT|O_TRUNC|O_RDWR,00644);
  for(j=0;j<2;j++){
    sprintf(buf,"IFloat uc_l%d[] LOCATE(\"edramtransient\") = {\n",j); 
    write(fd,buf,strlen(buf));
    sprintf(buf,"  %1.12e",*(uc_l[j])); 
    write(fd,buf,strlen(buf));
    for(i=1;i< MATRIX_SIZE * (local_chi/2);i++){
      sprintf(buf,",\n  %1.12e",*(uc_l[j]+i)); 
      write(fd,buf,strlen(buf));
    }
    sprintf(buf,"\n};\n"); 
    write(fd,buf,strlen(buf));
  }
  close(fd);

  sprintf(buf,"uc_nl.h");
  fd = open(buf,O_CREAT|O_TRUNC|O_RDWR,00644);
  for(j=0;j<2;j++){
    sprintf(buf,"IFloat uc_nl%d[] LOCATE(\"edramtransient\") = {\n",j); 
    write(fd,buf,strlen(buf));
    sprintf(buf,"  %1.12e",*(uc_nl[j])); 
    write(fd,buf,strlen(buf));
    for(i=1;i< MATRIX_SIZE * (local_chi/2);i++){
      sprintf(buf,",\n  %1.12e",*(uc_nl[j]+i)); 
      write(fd,buf,strlen(buf));
    }
    sprintf(buf,"\n};\n"); 
    write(fd,buf,strlen(buf));
  }
  close(fd);

#endif 
  
}

extern "C" void stag_destroy_dirac_buf_g(void)
{
  int i;
  ffree(tmpfrm);
  for ( i = 0; i < 2; i++){
  ffree(uc_l[i]);
  ffree(uc_nl[i]);
  }
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


extern "C"
void stag_dirac(Vector* b_v, Vector* a_v, int a_odd, int add_flag, int dir_flag)
{
  IFloat* b = (IFloat *)b_v;
  IFloat* a = (IFloat *)a_v;

  int c =0;
  int odd = 1 - a_odd;
  //-----------------------------------------------------------------
  //  Transfer chi's on faces.  
  //-----------------------------------------------------------------

  copy_buffer(countM[a_odd], (long)a, (long)Tbuffer[0], (long)ToffsetM[a_odd]);
  copy_buffer(countP[a_odd], (long)a, (long)Tbuffer[1], (long)ToffsetP[a_odd]);

  //make sure spinor field is in main memory before starting transfers
//  flush_cache_spinor(nflush, (long)a);

  SCUarg[1]->Addr( a + Xoffset[1]);
  SCUarg[5]->Addr( a + Xoffset[5]);
  SCUarg[2]->Addr( a + Xoffset[2]);
  SCUarg[6]->Addr( a + Xoffset[6]);
  SCUarg[3]->Addr( a + Xoffset[3]);
  SCUarg[7]->Addr( a + Xoffset[7]);
//    printf("Setting offsets done\n");
  
  sys_cacheflush(0);
//  for(int i = 0;i<2*NUM_DIR;i++) SCUarg[i]->Print();
//  exit(-34);
  for(int i = 0;i<2*NUM_DIR;i++) SCUarg[i]->Assert();
  SCUmulti->SlowStartTrans();
//  printf("SCUmulti started\n");
//  save_reg((long)intreg, (long)dreg);
  //-----------------------------------------------------------------
  //do first local computations
  //-----------------------------------------------------------------


#ifdef CPP
  dirac_cmv_l( local_chi/2, (long)chi_l[odd], (long)uc_l[odd], 
	       (long)a, (long)tmpfrm);
#else
  dirac_cmv( local_chi/2, (long)chi_l[odd], (long)uc_l[odd], 
	       (long)a, (long)tmpfrm);
#endif


  //-----------------------------------------------------------------
  // check to see if transfers are done
  //-----------------------------------------------------------------
//  for(int i = 0;i<2*NUM_DIR;i++) SCUarg[i]->Assert();
  SCUmulti->TransComplete();
//  printf("SCUmulti ended\n");

  //-----------------------------------------------------------------
  //do the computations involving non-local spinors
  //-----------------------------------------------------------------

#if 1
  dirac_cmv( non_local_chi/2, (long)chi_nl[odd], (long)uc_nl[odd], 
		(long)c, (long)tmpfrm);
#else
  dirac_cmv_nl( non_local_chi/2, (long)chi_nl[odd], (long)uc_nl[odd], 
		(long)c, (long)tmpfrm);
#endif


  //-----------------------------------------------------------------
  //do the sum of 8 temporary vectors at each lattice site
  //-----------------------------------------------------------------

  if ( add_flag == 0){
    dirac_sum( vol/2, (long)chi[odd], (long)tmpfrm, (long)b);
  }
  else{
    dirac_sum_acc( vol/2, (long)chi[odd], (long)tmpfrm, (long)b);
  }

  DiracOp::CGflops += 285*vol;


}

CPS_END_NAMESPACE

