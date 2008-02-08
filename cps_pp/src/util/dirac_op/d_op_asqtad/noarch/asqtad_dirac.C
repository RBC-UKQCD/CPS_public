//-------------------------------------------------------------------
//  $Id: asqtad_dirac.C,v 1.20 2008-02-08 18:35:07 chulwoo Exp $
//
//    12/21/02 HueyWen Lin, Chulwoo Jung
//
//   Asqtad Dirac operator for QCDOC. Communications and computations
//   are overlapped.
//   Uses many functions implemented by CJ for QCDOC
//
//
//-------------------------------------------------------------------
#include <config.h>
#include <util/gjp.h>
#include <comms/scu.h>
#include <comms/glb.h>
#include <util/lattice.h>
#include <util/dirac_op.h>
#include <util/vector.h>
#include <stdio.h>
#include <stdlib.h>
//#include <unistd.h>
#include <fcntl.h>
//#include <time_cps.h>
#include <math.h>
CPS_START_NAMESPACE

#undef num //temperary test for uc_l and uc_nl
//~~~~~~~~~~~~~

#define fat_link
#define naik
#define staple5
#define staple7


/*****************************************************************
 SIMUL switched on/off the hack CJ put in to help speed up the
simulation. if SIMUL is undefined, program writes the temporary arraies
for dirac operator. If SIMUL is defined, it will include (pre-calcuated)
temporary arraies and skip the generation of arraies.
******************************************************************/
/****************************************************************
CPP is a switch for using C++ routine for dirac_cmv.
*****************************************************************/

typedef struct gauge_agg{
  int src;
  int dest;
  IFloat mat[18];
} gauge_agg ;

void dirac_sum2_acc_cpp(int s, long chi, long tmpfrm, long b);
void dirac_cmv_cpp( int sites, long chi, long u,long a, long tmpfrm);
void dirac_cmv_jcw_agg_cpp( int sites, long chi, long u,long a, long tmpfrm);
void dirac_sum2_cpp(int s, long chi, long tmpfrm, long b);
void copy_buffer_cpp(int n, long src, long dest, long ptable);

extern "C" void copy_buffer(int n, long src, long dest, long ptable);

enum{VECT_LEN=6, VECT_LEN2=6, MATRIX_SIZE=18, SITE_LEN=72, NUM_DIR=8};
//---------------------------------------------------------------------
//  Arrays for send direction, receive direction, lattice size per
//  node, coordinate (for site where cluster of 8 matrices
//  is being assembled) and the coordinates of a nearest neighbor site.
//
//  0 == T, 1 == X, 2 == Y, 3 == Z
//--------------------------------------------------------------------


static int size[4];
static int split; 
static int coord[4];
static int coord_nn[4];
static int coord_knn[4];
static int vol;
static int non_local_chi_3;
static int non_local_chi;
static int local_chi;
static int local_chi_3;
static int buffer_order[3][4]; // never used ??
static int buffer_flush[3][8]; // never used ?? asign space in dirac_init but never be used after that.
static int nflush;

//---------------------------------------------------------------------
//  uc_l[0] points to a cluster of matrices per even site for local 
//  computations , arranged so that all parallel transport of spinors 
//  is accomplished by the same matrix times vector function.
//  The volume is devided by 2 areas: nn and 3rd nn part.
//  uc_l[1] is the same for odd sites.
//  uc_nl[0] is the same for even non-local sites.
//  uc_nl[1] is the same for odd non-local sites.
//---------------------------------------------------------------------
static gauge_agg * uc_l_agg[2];
static gauge_agg * uc_nl_agg[2];
static IFloat * uc_l[2];
static IFloat * uc_nl[2];

//------------------------------------------------------------------
// Allocate these arrays dynamicaly once the cache, noncached eDRAM
// malloc are available (should be changed according to volume)
//------------------------------------------------------------------

static IFloat *tmpfrm;

static IFloat *Tbuffer[3][2];

//-------------------------------------------------------------------
// end of stack based arrays which should be heap based
//-------------------------------------------------------------------
static int intreg[18] ;
static IFloat dreg[18] ;
static int * ToffsetP[3][2];
static int * ToffsetM[3][2];
static int countP[3][2];
static int countM[3][2];
//---------------------------------------------------------------------
//  pointers to storage area for color vectors from tp, xp, yp, zp, tm,
//  xm, ym, zm (p = plus, m = minus).  Indexed as 0-7
//---------------------------------------------------------------------
static IFloat * chi_off_node_total;
static IFloat * chi_off_node[3][8];
static IFloat * chi_off_node_p[3][8];
//------------------------------------------------------------------
//  pointer to array of pointers telling where color vectors are
//  located for cluster arrangement of lattice (even and odd).
//  chi_l[0] points to a list of local adresses for chi's needed 
//  to get an even site result from application of D and
//  to a temporary storage area where U_mu * chi's for each direction 
//  are stored.
//  the first half of storage size local_chi is for the nn,
//   while the second half of storage size local_chi_3 is for the 3rd nn.
//  chi_l[1] same for odd site.
//  chi_nl[0] same for computations involving non-local spinors
//------------------------------------------------------------------
static IFloat ** chi_l[2];
static IFloat ** chi_nl[2];

//------------------------------------------------------------------
//  Values for send and receive transfers. 0-7 correspond to
//  tp, xp, yp, zp, tm, xm, ym, zm (p = plus, m = minus).
//
//  Rarg (for SCU receives) never changes, since it receives into
//  the buffers for the specified direction.
//  SCUarg[3] is set up for different Xoffset and Toffset.
//------------------------------------------------------------------

#if 0
SCUDirArgIR * SCUarg[2*NUM_DIR];
SCUDirArgIR SCUargIR[2*NUM_DIR];
SCUDirArgMulti * SCUmulti;
SCUDirArgMulti SCUmultiIR;

SCUDirArgIR * SCUarg_1[2*NUM_DIR];
SCUDirArgIR SCUarg_1IR[2*NUM_DIR];
SCUDirArgMulti * SCUmulti_1;
SCUDirArgMulti SCUmulti_1IR;

SCUDMAInst *SCUDMAarg_p[NUM_DIR*4];
SCUDMAInst SCUDMAarg[NUM_DIR*4];

SCUDirArgIR * SCUarg_2[2*NUM_DIR];
SCUDirArgIR SCUarg_2IR[2*NUM_DIR];
SCUDirArgMulti * SCUmulti_2;
SCUDirArgMulti SCUmulti_2IR;
#endif

//------------------------------------------------------------------
//  The offset into the even or odd checkboard chi's and chi3's used for a send
//  from a given face.
//------------------------------------------------------------------

static int Xoffset[3][NUM_DIR];

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

static int CoordkNN( int nn, int k );

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

static  int blklen[NUM_DIR/2];
static  int numblk[NUM_DIR/2];
static  int stride[NUM_DIR/2];

int k;
//-------------------------------------------------------------------
//  Called by the lattice constructor
//  Fermion initializations: pointer tables 
//-------------------------------------------------------------------
extern "C"
void asqtad_dirac_init(Fasqtad * lat )
{
  gauge_field_addr = ( IFloat * ) lat->GaugeField();
  int r,c,i,j,m,n,nn;
  int local_count[2];
  int non_local_count[2];
  int local_count_3[2];
  int non_local_count_3[3][2];
  int off_node;
  int x[NUM_DIR/2];
  char buf[200];
  int fd;
  char *cname = "DiracOpAsqtad";
  char *fname = "dirac_init(const void *gauge)";

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
  split =  (  (size[0]>2 &&size[1]>2 && size[2]>2 && size[3] >2 ) ? 0 : 1 ); 
//  split =1;

  vol = size[0] * size[1] * size[2] * size[3];

  non_local_chi = 2*(size[0]*size[1]*size[2] + size[1]*size[2]*size[3]+
    size[2]*size[3]*size[0] + size[3]*size[0]*size[1]);

  non_local_chi_3 = 0;
  for(i=0;i<4;i++){
    if(size[i]<4)
      non_local_chi_3 += 4*vol/size[i];
    else
      non_local_chi_3 += 6*vol/size[i];
  }


  local_chi = NUM_DIR*vol - non_local_chi;
  local_chi_3 = NUM_DIR*vol - non_local_chi_3;
  nflush = vol/8;
  tmpfrm = (IFloat*)smalloc(NUM_DIR*2*vol/2 * VECT_LEN2 * sizeof(IFloat),
			    "tmpfrm",fname,cname);


  //-----------------------------------------------------------------
  //  Allocate 8 receive buffers for off-node vectors
  //-----------------------------------------------------------------


  chi_off_node_total = (IFloat*) 
    smalloc(3*non_local_chi*VECT_LEN*sizeof(IFloat)/2,
	    "chi_off_node_total",fname,cname);
  for ( j= 0; j < 3; j++ ){
    chi_off_node[j][0] = &(chi_off_node_total[ j*non_local_chi* VECT_LEN/2 ] ); 
    chi_off_node_p[j][0] = (IFloat *)(sizeof (IFloat)*j*non_local_chi* VECT_LEN/2 );
    for ( i = 1; i < NUM_DIR; i++ ){
      chi_off_node[j][i] = chi_off_node[j][i-1]+vol/(2*size[i%4])*VECT_LEN;
      chi_off_node_p[j][i] = chi_off_node_p[j][i-1]+vol/(2*size[i%4])*VECT_LEN;
    }
  }
  //-----------------------------------------------------------------
  //  Space for storage of pointers to chi's.  2 pointers per site,
  //  but split into even and odd groups for the first part of the
  //  computation (parallel transport of spinors). 17 pointers per site
  //  to obtain the result of the application of the dirac operator
  //
  //  The size of chi_l and chi_nl has been doubled for additional U*Chi and UUU*chi
  //-----------------------------------------------------------------

#ifndef SIMUL
  for ( i = 0; i < 2; i++ ){
    chi_l[i] = (IFloat**)smalloc((2*(local_chi+local_chi_3)/2)*sizeof(IFloat*),
				 "chi_l[i]", fname, cname);
    chi_nl[i] = (IFloat**)
      smalloc((2*(non_local_chi+non_local_chi_3)/2)*sizeof(IFloat*),
	      "chi_l[i]", fname, cname);
  }
#endif
  
  for ( i = 0; i < 2; i++){
    local_count[i] = 0;
    non_local_count[i] = 0;
    local_count_3[i] = 0;
    for( n = 0; n < 3; n++ )  non_local_count_3[n][i] = 0;

  }


#ifndef SIMUL
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

	    sg = coord[1] + size[1] * ( coord[2] + size[2] * ( coord[3] + size[3] * coord[0] ));
	    m = (2* NUM_DIR + 1) * (sg/2);

	    if ( CoordNN( n ) ) {		//  off-node
	      //----------------------------------------------------------
	      // Assembly written for double precision only, multiplication
	      // by sizeof(double) done to avoid a bitshift inside the
	      // high performance code
	      //----------------------------------------------------------
	      //pointer to source field (offset in the receive buffer)

	      *( chi_nl[ odd ]  +  2 * non_local_count[ odd ] )
		= ( IFloat * ) ( VECT_LEN * ( LexVector( coord_nn ) / 2 ) * sizeof(IFloat));

	      // pointer to temporary field where U*chi is stored
	      *( chi_nl[ odd ] + 2*non_local_count[ odd ] +1 ) =
				( IFloat * ) ( VECT_LEN2 * (LexVector( coord ) / 2 * 2*NUM_DIR+n )  * sizeof(IFloat));
	      non_local_count[odd]++;
	    }
	}
       }
      }
     }
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

	    sg = coord[1] + size[1] * ( coord[2] + size[2] * ( coord[3] + size[3] * coord[0] ));
	    m = (2* NUM_DIR + 1) * (sg/2);

	    if ( !CoordNN( n ) ) {		//  off-node
#ifndef SIMUL
	      //pointer to source field
	      *( chi_l[ odd ]  +  2 * local_count[ odd ] )
		= ( IFloat * ) ( VECT_LEN * ( LexVector( coord_nn ) / 2 ) * sizeof(IFloat));
	      // pointer to temporary field where U*chi is stored
	      *( chi_l[ odd ]  +  2 * local_count[ odd ] + 1) =
					( IFloat * ) ( VECT_LEN2 * (LexVector( coord ) / 2 *2*NUM_DIR+n) * sizeof(IFloat));
#endif
	      local_count[odd]++;
	    }  //else-on_node case
	}
       }
      }
     }
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

	    sg = coord[1] + size[1] * ( coord[2] + size[2] * ( coord[3] + size[3] * coord[0] ));
	    m = (2* NUM_DIR + 1) * (sg/2);



   //******chi3*********

  if ( CoordkNN( n, 3 ) ) {		//  chi3_off-node
	
    for (j=0; j<3; j++){
       if((coord_knn[n%4]==j&&n<4)||(coord_knn[n%4]==(size[n%4]-1-j)&& n>3)) {
	 *( chi_nl[ odd ] + 2 * non_local_count_3[j][odd]+(j+1)*non_local_chi )
	= ( IFloat * ) ( VECT_LEN * ( LexVector( coord_knn ) / 2 )
		                 * sizeof(IFloat));
        *( chi_nl[odd] + 2 * non_local_count_3[j][odd] + (j+1)*non_local_chi+1 ) = ( IFloat * ) ( VECT_LEN2 * (LexVector( coord ) / 2*2*NUM_DIR+n+8 )  * sizeof(IFloat));
	non_local_count_3[j][odd]++;
       }
    }

	      		       
          }  // end of chi3_off-node
	}
       }
      }
     }
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

	    sg = coord[1] + size[1] * ( coord[2] + size[2] * ( coord[3] + size[3] * coord[0] ));
	    m = (2* NUM_DIR + 1) * (sg/2);

  if ( !CoordkNN( n, 3 ) ) {		//  chi3_off-node
#ifndef SIMUL
	      //pointer to source field
	      *( chi_l[ odd ]  +  2 * local_count_3[ odd ] + local_chi )
		= ( IFloat * ) ( VECT_LEN * ( LexVector( coord_knn ) / 2 )
		                 * sizeof(IFloat));
	      // pointer to temporary field where U*chi is stored
	      *( chi_l[ odd ]  +  2 * local_count_3[ odd ] + local_chi + 1) =
	  	( IFloat * ) ( VECT_LEN2 * (LexVector( coord ) / 2 *2*NUM_DIR+n+8) * sizeof(IFloat));
#endif
	      local_count_3[odd]++;
	    }  //else-on_node case


	  } // for x[0] loop
	}
      }
    }// for x[3] loop
  }//for n loop
#endif /* SIMUL */
  FILE *fp;

#ifndef SIMUL
#if 0
  fp=Fopen("chi_l.h","w");
  for(j=0;j<2;j++){
    Fprintf(fp,"IFloat * chi_l%d[] LOCATE(\"edramtransient\") = {\n",j);
    Fprintf(fp," (IFloat *) %d",*(chi_l[j]));
    for(i=1;i< 2*((local_chi+local_chi_3)/2);i++){
      Fprintf(fp,",\n (IFloat *) %d",*(chi_l[j]+i));
    }
    Fprintf(fp,"\n};\n");
  }
  Fclose(fp);

  fp=Fopen("chi_nl.h","w");
    Fprintf(fp,"testing\n");
  for(j=0;j<2;j++){
    Fprintf(fp,"IFloat * chi_nl%d[] LOCATE(\"edramtransient\") = {\n",j);
    Fprintf(fp," (IFloat *) %d",*(chi_nl[j]));
    for(i=1;i< 2*((non_local_chi+ non_local_chi_3)/2);i++){
      Fprintf(fp,",\n (IFloat *) %d",*(chi_nl[j]+i));
    }
    Fprintf(fp,"\n};\n");
  }

  Fclose(fp);
#endif

#endif /* SIMUL */






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
  stride[2] = (VECT_LEN * size[0] * size[1] * ( size[2] - 1 ) / 2)* sizeof(IFloat) ;
  stride[3] = 0;

  //-------------------------------------------------------------------
  //  Calculate offsets for T transfers done one word at a time.
  //  We have plus (P) transfers for both the even and odd
  //  checkerboards.  Same for minus (M) transfers.
  //-------------------------------------------------------------------

  for ( k = 0; k < 3; k++ ) {
  for ( i = 0; i < 2; i++ ) {
#ifndef SIMUL_TBUF
    Tbuffer[k][i] = (IFloat*) 
      smalloc(size[1]*size[2]*size[3]*VECT_LEN*sizeof(IFloat)/2,
	      "Tbuffer[k][i]", fname, cname);
    ToffsetP[k][i] = (int*) smalloc(size[1]*size[2]*size[3]*sizeof(int)/2,
				    "ToffsetP[k][i]", fname, cname);
    ToffsetM[k][i] = (int*) smalloc(size[1]*size[2]*size[3]*sizeof(int)/2,
				    "ToffsetM[k][i]", fname, cname);
#endif
    countP[k][i] = 0;
    countM[k][i] = 0;
  }
 }
  for ( sg = 0; sg < vol; sg++ ) {

    odd = SetCoord( sg );
    sc = LexVector( coord );

    for ( int j = 0; j < 3; j++ ) {
      if ( coord[0] == j ) {
        *( ToffsetM[j][ odd ] + countM[j][ odd ] ) = VECT_LEN * ( sc / 2 );
        countM[j][ odd ]++;
      }
      if ( coord[0] == size[0] - 1 -j ) {
        *( ToffsetP[j][ odd ] + countP[j][ odd ] ) = VECT_LEN * ( sc / 2 );
        countP[j][ odd ]++;
      }
    }//end of j loop
  } //end of sg loop
  
  int vol3 = (size[1] * size[2] * size[3])/2;

  //-------------------------------------------------------------------
  //  Index i says data has been received from TP, XP, YP, ZP, TM, XM,
  //  YM, ZM
  //-------------------------------------------------------------------

  for ( i = 0; i < NUM_DIR; i++ ) {
    j = i % (NUM_DIR/2);
#if 0
      SCUarg[i + NUM_DIR] = &(SCUargIR[i+NUM_DIR]);
      SCUarg[i + NUM_DIR] ->Init(chi_off_node[2][i], scudir[i], SCU_REC,
		    VECT_LEN * sizeof(IFloat) * vol / ( 2 * size[j] ), 1, 0, IR_7);
      SCUarg[i + NUM_DIR] ->Assert();

    SCUDMAarg_p[(i+NUM_DIR)*2]  = &(SCUDMAarg[(i+NUM_DIR)*2]);
    SCUDMAarg_p[(i+NUM_DIR)*2] ->Init(chi_off_node[0][i],
      VECT_LEN * sizeof(IFloat) * vol / ( 2 * size[j] ), 1, 0);
    SCUDMAarg_p[(i+NUM_DIR)*2+1]  = &(SCUDMAarg[(i+NUM_DIR)*2+1]);
    SCUDMAarg_p[(i+NUM_DIR)*2+1] ->Init(chi_off_node[1][i],
      VECT_LEN * sizeof(IFloat) * vol / ( 2 * size[j] ), 1, 0);

    if( split ){
      SCUarg_1[i + NUM_DIR] = &(SCUarg_1IR[i+NUM_DIR]);
      SCUarg_1[i + NUM_DIR] ->Init(scudir[i],SCU_REC, &SCUDMAarg_p[(i+NUM_DIR)*2],1, IR_5);
      SCUarg_1[i + NUM_DIR] ->Assert();
      SCUarg_2[i + NUM_DIR] = &(SCUarg_2IR[i+NUM_DIR]);
      SCUarg_2[i + NUM_DIR] ->Init(scudir[i],SCU_REC, &SCUDMAarg_p[(i+NUM_DIR)*2+1],1, IR_6);
      SCUarg_2[i + NUM_DIR] ->Assert();
    } else {
      SCUarg_1[i + NUM_DIR] = &(SCUarg_1IR[i+NUM_DIR]);
      SCUarg_1[i + NUM_DIR] ->Init(scudir[i],SCU_REC, &SCUDMAarg_p[(i+NUM_DIR)*2],2, IR_5);
      SCUarg_1[i + NUM_DIR] ->Assert();
    }

      // change the size of buffer_flush or add 2 more??
      // never used ??  // 12288 or 384?
      buffer_flush[0][i] = VECT_LEN * sizeof(IFloat) * vol/ (12288 * size[j]);
      buffer_flush[1][i] = VECT_LEN * sizeof(IFloat) * vol/ (12288 * size[j]);
      buffer_flush[2][i] = VECT_LEN * sizeof(IFloat) * vol/ (12288 * size[j]);

    //send arguments
    if ((i == 0) || ( i == 4)){
      SCUarg[i] = &(SCUargIR[i]);
      if (size[j] >4)
        SCUarg[i] ->Init (Tbuffer[2][(4 - i)/4], scudir[i], SCU_SEND,
		       blklen[j], numblk[j], stride[j], IR_7 );
      else if (size[j] >2) 
        SCUarg[i] ->Init (Tbuffer[1][i/4], scudir[i], SCU_SEND,
		       blklen[j], numblk[j], stride[j], IR_7 );
      else
        SCUarg[i] ->Init (chi_off_node[0][4-i], scudir[i], SCU_SEND,
		       blklen[j], numblk[j], stride[j], IR_7 );
      SCUarg[i] ->Assert();

      SCUDMAarg_p[i*2] = &(SCUDMAarg[i*2]);
      SCUDMAarg_p[i*2] ->Init(Tbuffer[0][(4 - i)/4], 
		       blklen[j], numblk[j], stride[j]);
      SCUDMAarg_p[i*2+1] = &(SCUDMAarg[i*2+1]);

      if(size[j]>2)
        SCUDMAarg_p[i*2+1] ->Init(Tbuffer[1][(4 - i)/4], 
		       blklen[j], numblk[j], stride[j]);
      else
        SCUDMAarg_p[i*2+1] ->Init(Tbuffer[0][i/4], 
		       blklen[j], numblk[j], stride[j]);
    if( split ){
      SCUarg_1[i] = &(SCUarg_1IR[i]);
      SCUarg_1[i] ->Init(scudir[i],SCU_SEND,&SCUDMAarg_p[i*2],1,IR_5);
      SCUarg_1[i] ->Assert();
      SCUarg_2[i] = &(SCUarg_2IR[i]);
      SCUarg_2[i] ->Init(scudir[i],SCU_SEND,&SCUDMAarg_p[i*2+1],1,IR_6);
      SCUarg_2[i] ->Assert();
    } else {
      SCUarg_1[i] = &(SCUarg_1IR[i]);
      SCUarg_1[i] ->Init(scudir[i],SCU_SEND,&SCUDMAarg_p[i*2],2,IR_5);
      SCUarg_1[i] ->Assert();
    }

    }
    else{
      SCUarg[i] = &(SCUargIR[i]);
      if(size[j] >2)
        SCUarg[i] ->Init(( void * ) 0, scudir[i], SCU_SEND,
		       blklen[j], numblk[j], stride[j], IR_7 );
      else
        SCUarg[i] ->Init(chi_off_node[0][(i+4)%NUM_DIR], scudir[i], SCU_SEND,
           VECT_LEN * sizeof(IFloat) * vol / ( 2 * size[j]),1,0, IR_7 );
      SCUarg[i] ->Assert();

      SCUDMAarg_p[i*2] = &(SCUDMAarg[i*2]);
      SCUDMAarg_p[i*2] ->Init((void *)0, 
		       blklen[j], numblk[j], stride[j]);
      SCUDMAarg_p[i*2+1] = &(SCUDMAarg[i*2+1]);
      SCUDMAarg_p[i*2+1] ->Init((void *)0, 
		       blklen[j], numblk[j], stride[j]);
      if( split ){
        SCUarg_1[i] = &(SCUarg_1IR[i]);
        SCUarg_1[i] -> Init(scudir[i],SCU_SEND,&SCUDMAarg_p[i*2],1,IR_5);
        SCUarg_1[i] ->Assert();
        SCUarg_2[i] = &(SCUarg_2IR[i]);
        SCUarg_2[i] -> Init(scudir[i],SCU_SEND,&SCUDMAarg_p[i*2+1],1,IR_6);
        SCUarg_2[i] ->Assert();
      } else {
        SCUarg_1[i] = &(SCUarg_1IR[i]);
        SCUarg_1[i] -> Init(scudir[i],SCU_SEND,&SCUDMAarg_p[i*2],2,IR_5);
        SCUarg_1[i] ->Assert();
      }
    }
#endif
  }// end of NUM_DIR loop


#if 0
  SCUmulti = &SCUmultiIR;
  SCUmulti->Init(SCUarg, 2*NUM_DIR);
  SCUmulti_1 = &SCUmulti_1IR;
  SCUmulti_1->Init(SCUarg_1, 2*NUM_DIR);
  if( split ){
    SCUmulti_2 = &SCUmulti_2IR;
    SCUmulti_2->Init(SCUarg_2, 2*NUM_DIR);
  }
#endif

  //-------------------------------------------------------------------
  //  Need send offsets for various transfers.  The index for
  //  sends is TM, XM, YM, ZM, TP, XP, YP, ZP, since the
  //  transfers are indexed by the node data is received from.
  //-------------------------------------------------------------------
for ( k = 0; k < 3; k++ ) {
  Xoffset[k][0] = 0;
  Xoffset[k][1] = VECT_LEN * size[0] * (size[1] - 1-k) / 2;
  Xoffset[k][2] = VECT_LEN * size[0] * size[1] * (size[2] - 1-k) / 2;
  Xoffset[k][3] = VECT_LEN * size[0] * size[1] * size[2] * (size[3]-1-k) / 2;
  Xoffset[k][4] = 0;
  Xoffset[k][5] = VECT_LEN * size[0] *k /2;
  Xoffset[k][6] = VECT_LEN * size[0] * size[1]*k /2;
  Xoffset[k][7] = VECT_LEN * size[0] * size[1] * size[2]*k /2;
} // end of k loop
}

extern "C"
void asqtad_destroy_dirac_buf()
{
  char *cname = "";
  char *fname = "asqtad_destroy_dirac_buf()";
  int i,k;

  for ( i = 0; i < 2; i++ ) {
    for ( k = 0; k < 3; k++ ) {
#ifndef SIMUL_TBUF
      sfree(Tbuffer[k][i],"Tbuffer[k][i]",fname,cname);
      sfree(ToffsetP[k][i],"ToffsetP[k][i]",fname,cname);
      sfree(ToffsetM[k][i],"ToffsetM[k][i]",fname,cname);
#endif
    }
#ifndef SIMUL
    sfree(chi_nl[i],"chi_nl[i]",fname,cname);
    sfree(chi_l[i],"chi_l[i]",fname,cname);
#endif
  }
    
   sfree(chi_off_node_total,"chi_off_node_total",fname,cname);
  for ( i = 0; i < NUM_DIR; i++ ) {
#if 0
    delete SCUarg[i];
    delete SCUarg[i+8];
    delete SCUarg_1[i];
    delete SCUarg_1[i+8];
    for(k=0;k<4;k++)
      delete SCUDMAarg_p[i*4+k];
#endif
  }
  sfree(tmpfrm,"tmpfrm",fname,cname);
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

void TransfP(int off_node, int nflush_g,IFloat * v,IFloat * mtmp

, int n ){
	getPlusData(mtmp,v,18,(n+1)%4);
 }

void TransfM(int off_node, int nflush_g,IFloat * v,IFloat * mtmp, int n ){
	getMinusData(mtmp,v,18,(n+1)%4);
 }
void DaggerM (IFloat * w_t1,IFloat * v){
   	int i, j, c, r;
     	for ( r = 0; r < 3; r++ ) {
			for ( c = 0; c < 3; c++ ) {
        	  i = 6*r + 2*c;
			  j = 6*c + 2*r;
			  w_t1[i] = v[j];
			  w_t1[i+1] = -v[j+1];
           }
	       }
 }

void Parallel( IFloat * InMatrix,IFloat *OutMatrix  , int *coord, int dir_1, int dir_2, IFloat *v, IFloat * u,int multi_flag ,int tranfs_flag )
     {
 int nn, nflush_g=1, off_node; 	
 IFloat stp3_1[18], mtmp[18] ;
 if ( dir_1 <4){      //dir_1 positive
    //U_mu(x+mu)
      nn = ( dir_1 - 1 + 4 ) % 4;
      v = u + SITE_LEN *(LexGauge( coord ))+ MATRIX_SIZE * nn;
      DaggerM(stp3_1, v); 
 } else{       //dir_1 negative
    //U_mu(x+mu)~
      nn = ( dir_1 - 1 + 4 ) % 4;
      v = u + SITE_LEN *(LexGauge( coord ))+ MATRIX_SIZE * nn;
      for (int i=0; i<18; i++){
	stp3_1[i] = v[i];
      }
 }//end of dir_1
 
if ( multi_flag )
     mDotMEqual( OutMatrix, stp3_1, InMatrix ) ;
 else {
 for (int i=0; i<18; i++){
 OutMatrix[i] = stp3_1[i];
 }
 }


if ( tranfs_flag){
  if(dir_2<4){ 
     
      off_node = CoordNN( dir_2 );
      coord[ dir_2]= coord_nn[ dir_2];// x+nu-> x

      TransfM( off_node,  nflush_g, OutMatrix, mtmp,dir_2%4  );//transfer OutMatrix to x     
  }else {
     dir_2= dir_2%4; 
      off_node = CoordNN( dir_2+4 );
      coord[ dir_2]= coord_nn[ dir_2];// x+nu-> x

      TransfP( off_node,  nflush_g, OutMatrix, mtmp,dir_2%4  );//transfer OutMatrix to x     
  
  }
}

     } //Parallel


void Staple3_PP( IFloat * InMatrix,IFloat *OutMatrix  , int *coord, int n,
		 int nu, IFloat *v, IFloat * u, int nflush_g,int sum_flag )
    { 
 int  nn, i, off_node;
IFloat  stp3_0[18], stp3_2[18],  stp3_4[18], mtmp[18];

  CoordNN( n ); // 
  coord[n]= coord_nn[n];// x+mu

  //get U_nu(x+mu)~ and transfer it to  x+nu+mu
   Parallel( InMatrix,stp3_0, coord, nu+4, nu, v,  u, 1, 1 );

   // transfer U_nu(x+mu)~ to  x+nu
     off_node = CoordNN( n+4 );// x+mu +nu-> x+nu
     coord[n]= coord_nn[n];// x+nu
     TransfP( off_node,  nflush_g, stp3_0, mtmp,  n);//transfer stp3_0 to  x+nu
 
 //get U_mu(x+nu)& transfer U_mu(x+nu)*U_nu(x+mu)~ from x+nu to x site  
    Parallel( stp3_0,stp3_2, coord, n, nu+4, v,  u, 1, 1 );
//get U_nu(x+mu) & mulyiply U_nu(x+mu)*U_mu(x+nu)*U_nu(x+mu)~
    Parallel( stp3_2,stp3_4, coord, nu, 0, v,  u, 1, 0 );

    for ( i = 0; i < 18; i++ ) {
               OutMatrix[i] = stp3_4[i];
           }

    } //Staple3_PP

void Staple3_PN( IFloat * InMatrix,IFloat *OutMatrix  , int *coord, int n,
		 int nu, IFloat *v, IFloat * u, int nflush_g,int sum_flag ) 
{
int nn, i, off_node;
IFloat  stp3_0[18], stp3_2[18],  stp3_4[18], mtmp[18];

    
     CoordNN( n );
     coord[n]= coord_nn[n];// x+mu
     CoordNN( nu+4 );
     coord[nu]= coord_nn[nu];// x+mu-nu

 //U_nu(x+mu) at x+mu-nu & transfer it to  x-nu
      Parallel( InMatrix,stp3_0, coord, nu, n+4, v,  u, 1, 1 );

//U_mu(x-nu) *U_nu(x+mu)  at x-nu
       Parallel(stp3_0 ,stp3_2, coord, n, 0, v,  u, 1, 0 );

//U_nu(x-nu)~ * U_mu(x-nu) *U_nu(x+mu) & transfer it to  x
        Parallel(stp3_2 ,stp3_4, coord, nu+4, nu, v,  u, 1, 1 );

    for ( i = 0; i < 18; i++ ) {
               OutMatrix[i] = stp3_4[i];
           }

 
     
}//Staple3_PN


void Staple3_NP( IFloat * InMatrix,IFloat *OutMatrix  , int *coord, int n,int nu, IFloat *v, IFloat * u, int nflush_g,int sum_flag )
    { 
 int  nn, i, off_node;
IFloat  stp3_0[18], stp3_2[18],  stp3_4[18], mtmp[18];
  

  CoordNN( n+4 ); // 
  coord[n]= coord_nn[n];// x-mu

  //get U_nu(x-mu)~ and transfer it to  x+nu-mu
   Parallel( InMatrix,stp3_0, coord, nu+4, nu, v,  u, 1, 1 );

 //get U_mu(x-mu+nu)~& transfer U_mu(x-mu+nu~)*U_nu(x-mu)~ from x+nu-mu to x+nu site  
    Parallel( stp3_0,stp3_2, coord, n+4, n, v,  u, 1, 1 );

    //transfer U_mu(x-mu+nu~)*U_nu(x-mu)~ from x+nu to x site  
  CoordNN( nu+4 ); // 
  coord[nu]= coord_nn[nu];// x
 TransfP( off_node,  nflush_g, stp3_2 , mtmp, nu  );
  // TransfM( off_node,  nflush_g, stp3_2 , mtmp, nu  );
//get U_nu(x) &  U_nu(x)*U_mu(x+nu-mu)*U_nu(x-mu)~
    Parallel( stp3_2,stp3_4, coord, nu, 0, v,  u, 1, 0 );


    for ( i = 0; i < 18; i++ ) {
   
               OutMatrix[i] = stp3_4[i];
           }



}

void Staple3_NN( IFloat * InMatrix,IFloat *OutMatrix  , int *coord, int n,int nu, IFloat *v, IFloat * u, int nflush_g, int sum_flag )
    { 

 int nn, i, off_node;
IFloat  stp3_0[18], stp3_2[18],  stp3_4[18], mtmp[18];
 
     
     CoordNN( n+4 );
     coord[n]= coord_nn[n];// x-mu
     CoordNN( nu+4 );
     coord[nu]= coord_nn[nu];// x-mu-nu

 //*U_nu(x-mu-nu ) at x-mu-nu 
      Parallel( InMatrix,stp3_0, coord, nu, 0, v,  u, 1, 0 );
      // Parallel( InMatrix,stp3_0, coord, nu, 0, v,  u, 0, 0 );

//U_mu(x-mu-nu)~ *U_nu(x-mu-nu )at x-mu-nu & transfer it to  x-nu
       Parallel(stp3_0 ,stp3_2, coord, n+4, n, v,  u, 1, 1 );
 
       // get U_nu(x-nu)~ & transfer U_nu(x-nu)~ * U_mu(x-mu-nu)~ *U_nu(x-mu-nu ) to  x
        Parallel(stp3_2 ,stp3_4, coord, nu+4, nu, v,  u, 1, 1 );

    for ( i = 0; i < 18; i++ ) {

               OutMatrix[i] = stp3_4[i];
           
           }

   
}//Staple3_NN

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//first P refers to n dir positive, the second P refers to  ro dir positive 
//while the nu direction could be taken from 0 to 7 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void Staple5_PP( IFloat * InMatrix,IFloat *OutMatrix  , int *coord, int n,
		 int nu,int ro, IFloat *v, IFloat * u, int nflush_g )
 { 
    int  nn, i, off_node;
 IFloat stp5_0[18], stp5_1[18], stp5_2[18], 
   stp5_3[18] , stp5_4[18], stp5_5[18], mtmp[18];
     //U_ro(x+mu)~ and transfer it to  x+ro+mu
     CoordNN( n ); //
     coord[n]= coord_nn[n];// x+mu
      nn = ( ro - 1 + 4 ) % 4;
      v = u + SITE_LEN *(LexGauge( coord ))+ MATRIX_SIZE * nn;
    mDotMEqual( stp5_0, v, InMatrix  ) ;    

     off_node = CoordNN( ro );// x+mu -> x+mu+ro
     coord[ro]= coord_nn[ro];// x+mu+ro
     TransfM( off_node,  nflush_g, stp5_0, mtmp,  ro); //transfer stp5_0 to  x+mu+ro

     CoordNN( n+4 ); // 
     coord[n]= coord_nn[n];//  x+mu+ro -> x+ro
 
 
 ///U_nu(x+ro)  *U_mu(x+ro+nu)  * U_nu(x+mu+ro)~ * U_ro(x+mu)~  at x+ro 
 if( nu<4)
    Staple3_PP(stp5_0,stp5_4, coord, n,nu, v, u, nflush_g,1 );
else    
    Staple3_PN(stp5_0,stp5_4, coord, n,nu%4, v, u, nflush_g,1 );

 

   //transfer it to  x 
     off_node = CoordNN( ro +4 );//  x+ro ->x
     coord[ro]= coord_nn[ro];//
     TransfP( off_node,  nflush_g, stp5_4, mtmp, ro);

  ///U_ro(x) * U_nu(x+ro)  *U_mu(x+ro+nu)  * U_nu(x+mu+ro)~ * U_ro(x+mu)~     
      nn = ( ro - 1 + 4 ) % 4;
      v = u + SITE_LEN *(LexGauge( coord ))+ MATRIX_SIZE * nn;
      DaggerM(stp5_1,v);
       mDotMEqual( OutMatrix, stp5_1,  stp5_4) ;


 }

//n, ro 
void Staple5_PN( IFloat * InMatrix,IFloat *OutMatrix  , int *coord, int n,
		 int nu,int ro, IFloat *v, IFloat * u, int nflush_g )
{

    int  nn, i, off_node;
 IFloat stp5_0[18], stp5_1[18], stp5_2[18], 
   stp5_3[18] , stp5_4[18], stp5_5[18], mtmp[18];

  

    //U_ro(x+mu-ro) and transfer it to  x+mu-ro
     CoordNN( n ); //
     coord[n]= coord_nn[n];// x+mu
     CoordNN( ro+4 ); //
     coord[ro]= coord_nn[ro];// x+mu-ro

      nn = ( ro - 1 + 4 ) % 4;
      v = u + SITE_LEN *(LexGauge( coord ))+ MATRIX_SIZE * nn;
      DaggerM(stp5_1, v);
      mDotMEqual( stp5_0, stp5_1 , InMatrix  ) ; 
 
      /*     //transfer stp5_0 to  x+mu-ro+nu
     off_node = CoordNN( nu  );//  x+mu-ro ->x+mu-ro+nu
     coord[nu]= coord_nn[nu];//
     TransfP( off_node,  nflush_g, stp5_0, mtmp, nu); 

     CoordNN( nu+4 ); // 
     coord[nu]= coord_nn[nu];//  x+mu-ro+nu -> x+mu-ro*/
     CoordNN( n+4 ); // 
     coord[n]= coord_nn[n];//  x+mu-ro -> x-ro



 if( nu<4)
    Staple3_PP(stp5_0,stp5_4, coord, n,nu, v, u, nflush_g,1 );
else    
    Staple3_PN(stp5_0,stp5_4, coord, n,nu%4, v, u, nflush_g,1 );


      
      nn = ( ro - 1 + 4 ) % 4;
      v = u + SITE_LEN *(LexGauge( coord ))+ MATRIX_SIZE * nn;     
      mDotMEqual( OutMatrix, v,  stp5_4) ;

   //transfer it to  x 
     off_node = CoordNN( ro  );//  x-ro ->x
     coord[ro]= coord_nn[ro];//
     TransfM( off_node,  nflush_g, OutMatrix, mtmp, ro);

}

void Staple5_NP( IFloat * InMatrix,IFloat *OutMatrix  , int *coord, int n,
		 int nu,int ro, IFloat *v, IFloat * u, int nflush_g )
 { 
    int  nn, i, off_node;
 IFloat stp5_0[18], stp5_1[18], stp5_2[18], 
   stp5_3[18] , stp5_4[18], stp5_5[18], mtmp[18];
 

     //U_ro(x-mu)~ and transfer it to  x+ro+mu
     CoordNN( n+4 ); //
     coord[n]= coord_nn[n];// x-mu
      nn = ( ro - 1 + 4 ) % 4;
      v = u + SITE_LEN *(LexGauge( coord ))+ MATRIX_SIZE * nn;
      mDotMEqual( stp5_0, v, InMatrix  ) ; 

     off_node = CoordNN( ro );// x-mu -> x+mu+ro
     coord[ro]= coord_nn[ro];// x-mu+ro
     TransfM( off_node,  nflush_g, stp5_0, mtmp,  ro); //transfer stp5_0 to  x-mu+ro

     CoordNN( n ); // 
     coord[n]= coord_nn[n];//  x-mu+ro -> x+ro
 
 
 ///U_nu(x+ro)  *U_mu(x+ro+nu)  * U_nu(x+mu+ro)~ * U_ro(x+mu)~  at x+ro 
 if( nu<4)
    Staple3_NP(stp5_0,stp5_4, coord, n,nu, v, u, nflush_g,1 );
else    
    Staple3_NN(stp5_0,stp5_4, coord, n,nu%4, v, u, nflush_g,1 );

 

   //transfer it to  x 
     off_node = CoordNN( ro +4 );//  x+ro ->x
     coord[ro]= coord_nn[ro];//
     TransfP( off_node,  nflush_g, stp5_4, mtmp, ro);

  ///U_ro(x) * U_nu(x+ro)  *U_mu(x+ro+nu)  * U_nu(x+mu+ro)~ * U_ro(x+mu)~     
      nn = ( ro - 1 + 4 ) % 4;
      v = u + SITE_LEN *(LexGauge( coord ))+ MATRIX_SIZE * nn;
      DaggerM(stp5_1,v);
       mDotMEqual(OutMatrix, stp5_1,  stp5_4) ;



 }
void Staple5_NN( IFloat * InMatrix,IFloat *OutMatrix  , int *coord, int n,
		 int nu,int ro, IFloat *v, IFloat * u, int nflush_g )
 { 
    int  nn, i, off_node;
 IFloat stp5_0[18], stp5_1[18], stp5_2[18], 
   stp5_3[18] , stp5_4[18], stp5_5[18], mtmp[18];


   
    //U_ro(x-mu-ro) and transfer it to  x+mu-ro
     CoordNN( n+4 ); //
     coord[n]= coord_nn[n];// x-mu
     CoordNN( ro+4 ); //
     coord[ro]= coord_nn[ro];// x-mu-ro

      nn = ( ro - 1 + 4 ) % 4;
      v = u + SITE_LEN *(LexGauge( coord ))+ MATRIX_SIZE * nn;
      DaggerM(stp5_1, v);
      mDotMEqual( stp5_0,stp5_1 , InMatrix  ) ; 

     CoordNN( n ); // 
     coord[n]= coord_nn[n];//  x-mu-ro -> x-ro



 if( nu<4)
    Staple3_NP(stp5_0,stp5_4, coord, n,nu, v, u, nflush_g,1 );
else    
    Staple3_NN(stp5_0,stp5_4, coord, n,nu%4, v, u, nflush_g,1 );


      
      nn = ( ro - 1 + 4 ) % 4;
      v = u + SITE_LEN *(LexGauge( coord ))+ MATRIX_SIZE * nn;     
      mDotMEqual( OutMatrix,v,  stp5_4) ;

   //transfer it to  x 
     off_node = CoordNN( ro  );//  x-ro ->x
     coord[ro]= coord_nn[ro];//
     TransfM( off_node,  nflush_g, OutMatrix, mtmp, ro);


 }

//#if 0
#ifdef staple7

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//first P refers to n dir positive, the second P refers to  de dir positive 
//while the nu & ro direction could be taken from 0 to 7 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void Staple7_PP( IFloat * InMatrix,IFloat *OutMatrix  , int *coord, int n,
		 int nu,int ro,int de, IFloat *v, IFloat * u, int nflush_g )
 { 
    int  nn, i, off_node;
 IFloat stp5_0[18], stp5_1[18], stp5_2[18], 
   stp5_3[18] , stp5_4[18], stp5_5[18], mtmp[18];
     //U_ro(x+mu)~ and transfer it to  x+ro+mu
     CoordNN( n ); //
     coord[n]= coord_nn[n];// x+mu
      nn = ( de - 1 + 4 ) % 4;
      v = u + SITE_LEN *(LexGauge( coord ))+ MATRIX_SIZE * nn;
    mDotMEqual( stp5_0, v, InMatrix  ) ;    

     off_node = CoordNN( de );// x+mu -> x+mu+ro
     coord[de]= coord_nn[de];// x+mu+ro
     TransfM( off_node,  nflush_g, stp5_0, mtmp,  de); //transfer stp5_0 to  x+mu+ro

     CoordNN( n+4 ); // 
     coord[n]= coord_nn[n];//  x+mu+de -> x+de
 
 
 ///U_nu(x+ro)  *U_mu(x+ro+nu)  * U_nu(x+mu+ro)~ * U_ro(x+mu)~  at x+ro 
 if( ro<4)
    Staple5_PP(stp5_0,stp5_4, coord, n,nu,ro, v, u, nflush_g );
else    
    Staple5_PN(stp5_0,stp5_4, coord, n,nu,ro%4, v, u, nflush_g );

 

   //transfer it to  x 
     off_node = CoordNN( de +4 );//  x+de ->x
     coord[de]= coord_nn[de];//
     TransfP( off_node,  nflush_g, stp5_4, mtmp, de);

   
      nn = ( de - 1 + 4 ) % 4;
      v = u + SITE_LEN *(LexGauge( coord ))+ MATRIX_SIZE * nn;
      DaggerM(stp5_1,v);
       mDotMEqual( OutMatrix, stp5_1,  stp5_4) ;


 }


void Staple7_PN( IFloat * InMatrix,IFloat *OutMatrix  , int *coord, int n,
		 int nu,int ro,int de, IFloat *v, IFloat * u, int nflush_g )
{

    int  nn, i, off_node;
 IFloat stp5_0[18], stp5_1[18], stp5_2[18], 
   stp5_3[18] , stp5_4[18], stp5_5[18], mtmp[18];

  

    //U_de(x+mu-de) and transfer it to  x+mu-de
     CoordNN( n ); //
     coord[n]= coord_nn[n];// x+mu
     CoordNN( de+4 ); //
     coord[de]= coord_nn[de];// x+mu-de

      nn = ( de - 1 + 4 ) % 4;
      v = u + SITE_LEN *(LexGauge( coord ))+ MATRIX_SIZE * nn;
      DaggerM(stp5_1, v);
      mDotMEqual( stp5_0, stp5_1 , InMatrix  ) ; 
 
     CoordNN( n+4 ); // 
     coord[n]= coord_nn[n];//  x+mu-de -> x-de



 if( ro<4)
    Staple5_PP(stp5_0,stp5_4, coord, n,nu,ro, v, u, nflush_g);
 
else    
    Staple5_PN(stp5_0,stp5_4, coord, n,nu,ro%4, v, u, nflush_g );


      
      nn = ( de - 1 + 4 ) % 4;
      v = u + SITE_LEN *(LexGauge( coord ))+ MATRIX_SIZE * nn;     
      mDotMEqual( OutMatrix, v,  stp5_4) ;

   //transfer it to  x 
     off_node = CoordNN( de  );//  x-de ->x
     coord[de]= coord_nn[de];//
     TransfM( off_node,  nflush_g, OutMatrix, mtmp, de);

}


void Staple7_NP( IFloat * InMatrix,IFloat *OutMatrix  , int *coord, int n,
		 int nu,int ro,int de, IFloat *v, IFloat * u, int nflush_g )
{

 
    int  nn, i, off_node;
 IFloat stp5_0[18], stp5_1[18], stp5_2[18], 
   stp5_3[18] , stp5_4[18], stp5_5[18], mtmp[18];


 /*    //U_de(x-mu)~ and transfer it to  x+de+mu
     CoordNN( n+4 ); //
     coord[n]= coord_nn[n];// x-mu
      nn = ( de - 1 + 4 ) % 4;
      v = u + SITE_LEN *(LexGauge( coord ))+ MATRIX_SIZE * nn;
      mDotMEqual( stp5_0, v, InMatrix  ) ; 

     off_node = CoordNN( de );// x-mu -> x+mu+de
     coord[ro]= coord_nn[de];// x-mu+de
     TransfM( off_node,  nflush_g, stp5_0, mtmp,  de); //transfer stp5_0 to  x-mu+ro

     CoordNN( n ); // 
     coord[n]= coord_nn[n];//  x-mu+de -> x+de
 
 

 if( ro<4)
    Staple5_NP(stp5_0,stp5_4, coord, n,nu, ro, v, u, nflush_g );
else    
    Staple5_NN(stp5_0,stp5_4, coord, n,nu, ro%4, v, u, nflush_g );

 

   //transfer it to  x 
     off_node = CoordNN( de +4 );//  x+ro ->x
     coord[de]= coord_nn[de];//
     TransfP( off_node,  nflush_g, stp5_4, mtmp, de);

  ///U_ro(x) * U_nu(x+ro)  *U_mu(x+ro+nu)  * U_nu(x+mu+ro)~ * U_ro(x+mu)~     
      nn = ( de - 1 + 4 ) % 4;
      v = u + SITE_LEN *(LexGauge( coord ))+ MATRIX_SIZE * nn;
      DaggerM(stp5_1,v);
       mDotMEqual(OutMatrix, stp5_1,  stp5_4) ;

 */
 

     //U_de(x-mu)~ and transfer it to  x+de+mu
     CoordNN( n+4 ); //
     coord[n]= coord_nn[n];// x-mu
      nn = ( de - 1 + 4 ) % 4;
      v = u + SITE_LEN *(LexGauge( coord ))+ MATRIX_SIZE * nn;
      mDotMEqual( stp5_0, v, InMatrix  ) ; 

     off_node = CoordNN( de );// x-mu -> x+mu+de
     coord[de]= coord_nn[de];// x-mu+de
     TransfM( off_node,  nflush_g, stp5_0, mtmp,  de); //transfer stp5_0 to  x-mu+de

     CoordNN( n ); // 
     coord[n]= coord_nn[n];//  x-mu+de -> x+de
 
 
 
 if( ro<4)
    Staple5_NP(stp5_0,stp5_4, coord, n,nu,ro, v, u, nflush_g);
else    
    Staple5_NN(stp5_0,stp5_4, coord, n,nu, ro%4, v, u, nflush_g );



   //transfer it to  x 
     off_node = CoordNN( de +4 );//  x+de ->x
     coord[de]= coord_nn[de];//
     TransfP( off_node,  nflush_g, stp5_4, mtmp, de);

   
      nn = ( de - 1 + 4 ) % 4;
      v = u + SITE_LEN *(LexGauge( coord ))+ MATRIX_SIZE * nn;
      DaggerM(stp5_1,v);
       mDotMEqual(OutMatrix, stp5_1,  stp5_4) ;

 

 }


void Staple7_NN( IFloat * InMatrix,IFloat *OutMatrix  , int *coord, int n,
		 int nu,int ro,int de, IFloat *v, IFloat * u, int nflush_g )
  { 
    int  nn, i, off_node;
 IFloat stp5_0[18], stp5_1[18], stp5_2[18], 
   stp5_3[18] , stp5_4[18], stp5_5[18], mtmp[18];


   
    //U_de(x-mu-de) and transfer it to  x+mu-de
     CoordNN( n+4 ); //
     coord[n]= coord_nn[n];// x-mu
     CoordNN( de+4 ); //
     coord[de]= coord_nn[de];// x-mu-de

      nn = ( de - 1 + 4 ) % 4;
      v = u + SITE_LEN *(LexGauge( coord ))+ MATRIX_SIZE * nn;
      DaggerM(stp5_1, v);
      mDotMEqual( stp5_0,stp5_1 , InMatrix  ) ; 

     CoordNN( n ); // 
     coord[n]= coord_nn[n];//  x-mu-de -> x-de



 if( ro<4)
    Staple5_NP(stp5_0,stp5_4, coord, n,nu,ro, v, u, nflush_g );
else    
    Staple5_NN(stp5_0,stp5_4, coord, n,nu,ro%4, v, u, nflush_g );


      
      nn = ( de - 1 + 4 ) % 4;
      v = u + SITE_LEN *(LexGauge( coord ))+ MATRIX_SIZE * nn;     
      mDotMEqual( OutMatrix,v,  stp5_4) ;

   //transfer it to  x 
     off_node = CoordNN( de  );//  x-de ->x
     coord[de]= coord_nn[de];//
     TransfM( off_node,  nflush_g, OutMatrix, mtmp, de);


  
}

#endif //staple7
extern "C"
void asqtad_dirac_init_g()
{

  IFloat * u = gauge_field_addr;
  int c,i,j,m,n,r;
  int off_node;
  int local_count[2];
  int nflush_g = 2;
  int non_local_count[2];
  int local_count_3[2];
int non_local_count_3[3][2];
  int x[NUM_DIR/2];
  char *cname = "";
  char *fname = "asqtad_dirac_init_g()";

  //--------------------------------------------------------------------
  // c1 -> one link; c2 -> 3-link; c3 -> 3-link staple; c5 -> 5-link staple;
  // c7 -> 7-link staple; c6 -> 5-link "straight" staple
  //--------------------------------------------------------------------
 
  IFloat c1 = GJP.KS_coeff();
  IFloat c2 = GJP.Naik_coeff();
  IFloat c3 = -GJP.staple3_coeff();
  IFloat c5 = GJP.staple5_coeff();
  IFloat c7 = -GJP.staple7_coeff();
  IFloat c6 = GJP.Lepage_coeff(); 

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
  //  be stored.  same for w3 (for UUU matrix).
  //-----------------------------------------------------------

  IFloat * v;
  IFloat * w, wp1[18];
  IFloat * w3;

  IFloat  w_t1[18], w_t2[18], w_t3[18];
  IFloat mtmp[18] ;

  //-----------------------------------------------------------
  //  SCU transfer structure to get links from off node and a
  //  location where one link matrix can be stored.
  //-----------------------------------------------------------
  VRB.Func(cname,fname);
#if 0
  size[0] = GJP.TnodeSites();
  size[1] = GJP.XnodeSites();
  size[2] = GJP.YnodeSites();
  size[3] = GJP.ZnodeSites();

  vol = size[0] * size[1] * size[2] * size[3];

  non_local_chi = 2*(size[0]*size[1]*size[2] + size[1]*size[2]*size[3]+
    size[2]*size[3]*size[0] + size[3]*size[0]*size[1]);

  non_local_chi_3 = 6*(size[0]*size[1]*size[2] + size[1]*size[2]*size[3]+
    size[2]*size[3]*size[0] + size[3]*size[0]*size[1]);

  local_chi = NUM_DIR*vol - non_local_chi;
  local_chi_3 = NUM_DIR*vol - non_local_chi_3;
#endif

#ifdef SIMUL_U
  uc_l[0] = uc_l0;
  uc_l[1] = uc_l1;
  uc_nl[0] = uc_nl0;
  uc_nl[1] = uc_nl1;
#endif

  //-----------------------------------------------------------------
  //  Allocate space for two copies of the gauge fields on this node
  //-----------------------------------------------------------------
  for ( i = 0; i < 2; i++ ){

#ifndef SIMUL_U
    uc_l[i]  = (IFloat*)
      smalloc(MATRIX_SIZE*((local_chi+local_chi_3)/2)*sizeof(IFloat),
	      "uc_l[i]", fname, cname);
    for(int j=0;j<MATRIX_SIZE*(local_chi+local_chi_3)/2;j++)
	uc_l[i][j]=0.;
#endif
    uc_l_agg[i]  = (gauge_agg*)
      smalloc(((local_chi+local_chi_3)/2)*sizeof(gauge_agg),
	      "uc_l_agg[i]", fname, cname);

 
#ifndef SIMUL_U
    uc_nl[i] = (IFloat*)
      smalloc(MATRIX_SIZE*((non_local_chi+non_local_chi_3)/2) * sizeof(IFloat),
	      "uc_nl[i]", fname, cname);
    for(int j=0;j<MATRIX_SIZE*(non_local_chi+non_local_chi_3)/2;j++)
	uc_nl[i][j]=0.;
#endif
    uc_nl_agg[i] = (gauge_agg*)
      smalloc(((non_local_chi+non_local_chi_3)/2)*sizeof(gauge_agg),
	      "uc_nl_agg[i]", fname, cname);

  }

#ifndef SIMUL_U



  for ( i = 0; i < 2; i++){
    local_count[i] = 0;
    non_local_count[i] = 0;
    local_count_3[i]=0;
   for (j=0; j<3; j++) non_local_count_3[j][i] = 0;
       

  }

  //-----------------------------------------------------------
  //  Copy links in POSTIVE directions.  These are all on-node
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
	    sg = coord[1] + size[1] * ( coord[2] + size[2] * ( coord[3] +  size[3] * coord[0] ));

    //---------------------------------------------------------------
     //U_mu(x) ~
     //---------------------------------------------------------------
       nn = ( n - 1 + 4 ) % 4;		//  nn is physics system order
       v = u + SITE_LEN *(LexGauge( coord ))+ MATRIX_SIZE * nn;
	    //v = u + SITE_LEN * sg + MATRIX_SIZE * nn;
       //storage the uc fields in the same order of chi's
      if ( CoordNN( n ) ) {		// chi(x+mu) off-node
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
	else {
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
      } // end of else

 for (i = 0; i < 4 ; i++) coord[i] = x[i];
 CoordNN( n );
 coord[n]=coord_nn[n];  // x+mu
 CoordNN( n );
 coord[n]=coord_nn[n]; // x+2m

 //U_mu(x+mu) ~
 nn =  ( n - 1 + 4 ) % 4;          
 v = u + SITE_LEN * (LexGauge( coord )) + MATRIX_SIZE * nn;
 DaggerM( w_t2, v);  

  off_node= CoordNN( n+4 ); // x+2mu -> x+mu
  TransfP( off_node,  nflush_g, w_t2, mtmp,  n);

 //U_mu(x+mu) ~
  coord[n]=coord_nn[n];
  nn =  ( n - 1 + 4 ) % 4;          
 v = u + SITE_LEN * (LexGauge( coord )) + MATRIX_SIZE * nn;
 DaggerM( w_t1, v); 

  mDotMEqual( w_t3, w_t1, w_t2  ) ; //on x+mu node

  off_node= CoordNN( n+4 ); // x+mu -> x
  coord[n]=coord_nn[n];
  TransfP( off_node,  nflush_g, w_t3, mtmp,  n);
     //---------------------------------------------------------------
    //U_mu(x) ~U_mu(x+mu) ~U_mu(x+2mu) ~   (1 transformation  max)
     //---------------------------------------------------------------


   for (i = 0; i < 4 ; i++) coord[i] = x[i];
   if ( CoordkNN (n, 3)) {         // 3rd nn off-node
  
   for (j=0; j<3; j++){

    if ((coord_knn[n%4]==j)) {
	w3 = uc_nl[odd]+MATRIX_SIZE*(non_local_count_3[j][odd]+(j+1)* non_local_chi/2);
	non_local_count_3[j][odd]++;
       }
   }

     }
   else {                      // 3rd nn on-node
    w3 = uc_l[odd] + MATRIX_SIZE * (local_chi/2 + local_count_3[odd]);
    local_count_3[odd]++;
   }


  

   mDotMEqual(w3 ,  w ,  w_t3  ) ;  

  for ( i = 0; i < 18; i++ ) {
    w3[i] = c2 * w3[i]; 

  }
   // mDotMEqual(w3 ,  w_t3 ,  w  ) ;  //wrong one




 int nu;
 IFloat  stp3_0[18],  stp3_2[18],  stp3_4[18], stp3_5[18];
 IFloat stp5_0[18], stp5_1[18], stp5_2[18], stp5_3[18] , stp5_4[18], stp5_5[18];
 IFloat UnitMatrix[18], stp5_f[18], stp3_f[18] ; 
 IFloat stp7_0[18], stp7_f[18], stp7_2[18];
 IFloat stp6_f[18];

   for ( i = 0; i < 18; i++ ) {
               stp3_5[i] = 0.;
	       stp5_f[i] = 0.;
	       stp6_f[i] = 0.;
	       stp3_f[i] = 0.;
	       stp7_f[i] = 0.;
	       UnitMatrix[i] = 0.;
           }


    UnitMatrix[0]=1.;
    UnitMatrix[8]=1.;
    UnitMatrix[16]=1.;

for (nu=0; nu < 8 ; nu++){
 if  (nu%4 != n){
    for (i = 0; i < 4 ; i++) coord[i] = x[i]; 
if (nu<4){
      Staple3_PP (UnitMatrix,stp3_4, coord, n,nu, v,u,nflush_g,1 );
      Staple5_PP (UnitMatrix,stp5_2, coord, n,nu, nu, v,u,nflush_g );
}
else{ 
     Staple3_PN (UnitMatrix,stp3_4, coord, n,nu%4, v,u,nflush_g,1 );
     Staple5_PN (UnitMatrix,stp5_2, coord, n,nu,nu%4, v,u,nflush_g );
}

         for ( i = 0; i < 18; i++ ) {
               stp3_5[i] += stp3_4[i];
               stp6_f[i] += stp5_2[i];
           }
 }}





for (int ro=0; ro < 4 ; ro++){  //positive ro direction
if  (ro%4 != n)
   {  
 for (nu=0; nu < 8 ; nu++){ // 8 direction
if  ( ( (nu%4)!= n )&&( (nu%4)!=(ro%4)) )
   {  
     for (i = 0; i < 4 ; i++) coord[i] = x[i]; 

  Staple5_PP (UnitMatrix ,stp5_2 , coord, n, nu, ro, v, u, nflush_g  );

       for ( i = 0; i < 18; i++ ) {
	 stp5_f[i] += stp5_2[i];}
}
}}}	

 for (int ro=0; ro < 4 ; ro++){ //negative ro direction
if  (ro != n)
   {  
 for (nu=0; nu < 8 ; nu++){ //8 direction
if  ( ( (nu%4)!= n )&&( (nu%4)!=(ro%4)) )
   {  
     for (i = 0; i < 4 ; i++) coord[i] = x[i]; 

     Staple5_PN(UnitMatrix ,stp5_2 , coord, n, nu, ro, v, u, nflush_g  );

       for ( i = 0; i < 18; i++ ) {
	 stp5_f[i] += stp5_2[i];}
   }
 }
   }
 }



for (int de=0; de < 4 ;de++){  //positive de direction
  if  (de%4 != n ){

for (int ro=0; ro < 8 ; ro++){  //all ro direction
if  (ro%4 != n && ro%4 != de%4)
   {  
 for (nu=0; nu < 8 ; nu++){ // 8 direction
if  (nu%4 != n && nu%4 != ro%4 && nu%4 != de%4)
   {  
     for (i = 0; i < 4 ; i++) coord[i] = x[i]; 

     Staple7_PP (UnitMatrix ,stp7_2 , coord, n, nu, ro,de, v, u, nflush_g  );
 
       for ( i = 0; i < 18; i++ ) {
	 stp7_f[i] += stp7_2[i];}
}
 }}} }}


for (int de=0; de < 4 ;de++){  //negative de direction
  if  (de%4 != n){

for (int ro=0; ro < 8 ; ro++){  //all ro direction
if  (ro%4 != n && ro%4 != de%4)
   {  
 for (nu=0; nu < 8 ; nu++){ // 8 direction
if  (nu%4 != n && nu%4 != ro%4 && nu%4 != de%4)
   {  
     for (i = 0; i < 4 ; i++) coord[i] = x[i]; 

     Staple7_PN (UnitMatrix ,stp7_2 , coord, n, nu, ro,de, v, u, nflush_g  );
 
       for ( i = 0; i < 18; i++ ) {
	 stp7_f[i] += stp7_2[i];}
}
 }}} }}











 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 // add the fat link term to one link
 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
    for ( i = 0; i < 18; i++ ) {

      w[i]*= c1;
      w[i] = w[i] + c3 * stp3_5[i]+ c5 * stp5_f[i] + c6
*stp6_f[i] + c7 * stp7_f[i] ;
      /*      w[i]+=stp3_5[i];
        w[i]+=stp5_f[i];
	 w[i]+=stp7_f[i];
*/
           }



        }//end of for x[0] loop
       }
      }
    } //end of for x[3] loop
  }//end of the n-loop (sum_postive mu _

 

   //-----------------------------------------------------------------
    //  Copy links in NEGATIVE directions.  Some are off-node
    //  Add an overall minus sign
    //-----------------------------------------------------------------
    for ( n = 0; n < NUM_DIR/2; n++ ) {
 
  for (x[3] = 0; x[3] < size[3]; x[3]++){
	for (x[2] = 0; x[2] < size[2]; x[2]++){
	  for (x[1] = 0; x[1] < size[1]; x[1]++){
	    for (x[0] = 0; x[0] < size[0]; x[0]++){

	      for (i = 0; i < 4 ; i++) coord[i] = x[i];
	      odd = ( coord[0] + coord[1] + coord[2] + coord[3] ) % 2;

	      nn = ( n - 1 + 4 ) % 4;		//  nn is physics system order
      	//---------------------------------------------------------------
	      //U_mu (x-mu)
	      //---------------------------------------------------------------
        off_node = CoordNN( n + 4 );

	      v = u + SITE_LEN * LexGauge( coord_nn ) + MATRIX_SIZE * nn;

        if ( off_node ) {		//  off-node
	  w = uc_nl[odd] + MATRIX_SIZE * non_local_count[odd];
           non_local_count[odd]++;
				TransfM( off_node,  nflush_g, v, mtmp,  n);
				v = mtmp;
	      }
	else{
		  w = uc_l[odd] + MATRIX_SIZE * local_count[odd];
		  local_count[odd]++;
	      }  //end of "else"

	//U_mu (x-mu) -> x  
 for ( i = 0; i < 18; i++ ) {
         	w[i] = -v[i];
           }

 for (i = 0; i < 4 ; i++) coord[i] = x[i];
 CoordNN( n + 4);
 coord[n]=coord_nn[n];  // x-mu
 CoordNN( n + 4);
 coord[n]=coord_nn[n];  // x-2mu
 off_node=CoordNN( n + 4); //coord_nn= x-3mu

 //U_mu (x-3mu) to x2mu
     v = u + SITE_LEN * LexGauge( coord_nn ) + MATRIX_SIZE * nn;
     for ( i = 0; i < 18; i++ ) {
         	w_t2[i] = v[i];
           }

     TransfM( off_node,  nflush_g, w_t2, mtmp,  n);
 //U_mu (x-2mu)
     v = u + SITE_LEN * LexGauge( coord ) + MATRIX_SIZE * nn;
     for ( i = 0; i < 18; i++ ) {
         	w_t1[i] = v[i];
           }
  // /U_mu (x-2mu)*U_mu (x-3mu) on x-2mu to x
    mDotMEqual(w_t3, w_t1, w_t2  ) ; //on x-2mu
    off_node= CoordkNN( n,2 ); // x-2mu -> x //U_mu (x-mu) transfered to x too 
    TransfM( off_node,  nflush_g, w_t3, mtmp,  n);

   //---------------------------------------------------------------
   //U_mu(x-mu)*U_mu(x-2mu)*U_mu(x-3mu)
   //---------------------------------------------------------------
   for (i = 0; i < 4 ; i++) coord[i] = x[i];
   if ( CoordkNN (n + 4, 3)) {      // 3rd nn off-node

   for (j=0; j<3; j++){

    if ((coord_knn[n%4]==(size[n%4]-1-j))) {
	w3 = uc_nl[odd]+MATRIX_SIZE*(non_local_count_3[j][odd]+(j+1)* non_local_chi/2);
	non_local_count_3[j][odd]++;
       }
   }
   
   }
   else {                     // 3rd nn om-node
    w3 = uc_l[odd] + MATRIX_SIZE * (local_count_3[odd] + local_chi/2);
    local_count_3[odd]++;
   }

   // mDotMEqual( mtmp,  w_t1,  w_t2) ;
   mDotMEqual( w3 ,  w ,  w_t3 ) ;

  for ( i = 0; i < 18; i++ ) {
    w3[i] = c2 * w3[i]; 

  }


 int nu;
 IFloat  stp3_0[18],stp3_1[18],  stp3_2[18], stp3_3[18], stp3_4[18], stp3_5[18];
IFloat stp5_0[18], stp5_1[18], stp5_2[18], stp5_3[18] , stp5_4[18], stp5_5[18];
IFloat UnitMatrix[18], stp5_f[18], stp3_f[18], stp7_2[18],stp7_f[18],stp6_f[18] ; 




   for ( i = 0; i < 18; i++ ) {
               stp3_5[i] = 0;
	       stp5_f[i] = 0;
	       stp3_f[i] = 0;
	       stp6_f[i] = 0;
	       stp7_f[i] = 0;
	       UnitMatrix[i] = 0; 
   
           }

    UnitMatrix[0]=1;
    UnitMatrix[8]=1;
    UnitMatrix[16]=1;


for (nu=0; nu < 8 ; nu++){
 if  (nu%4 != n){
    for (i = 0; i < 4 ; i++) coord[i] = x[i]; 
if (nu<4){
      Staple3_NP (UnitMatrix,stp3_4, coord, n,nu, v,u,nflush_g,1 );
      Staple5_NP (UnitMatrix,stp5_2, coord, n,nu, nu, v,u,nflush_g);
}
else {
     Staple3_NN (UnitMatrix,stp3_4, coord, n,nu%4, v,u,nflush_g,1 );
     Staple5_NN (UnitMatrix,stp5_2, coord, n,nu,nu%4, v,u,nflush_g);
}

         for ( i = 0; i < 18; i++ ) {
               stp3_5[i] += stp3_4[i];
               stp6_f[i] += stp5_2[i];
           }
 }}


 for (int ro=0; ro < 4 ; ro++){  //positive ro direction
if  (ro%4 != n)
   {  
 for (nu=0; nu < 8 ; nu++){ // 8 direction
if  (nu%4 != n && nu%4 != ro%4)
   {  
     for (i = 0; i < 4 ; i++) coord[i] = x[i]; 
 Staple5_NP (UnitMatrix ,stp5_2 , coord, n, nu, ro, v, u, nflush_g  );

 //summation
       for ( i = 0; i < 18; i++ ) {
          stp5_f[i] += stp5_2[i]; }
  

   }
 }
   
   }}	
		


 for (int ro=0; ro < 4 ; ro++){ //nrgative ro direction
if  (ro != n)
   {  
 for (nu=0; nu < 8 ; nu++){ //8 direction
if  (nu%4 != n && nu%4 != ro%4)
   {  
     for (i = 0; i < 4 ; i++) coord[i] = x[i]; 

 Staple5_NN(UnitMatrix ,stp5_2 , coord, n, nu, ro, v, u, nflush_g  );

       for ( i = 0; i < 18; i++ ) {
	 stp5_f[i] += stp5_2[i];}


   }
 }
   }
 }


	    

for (int de=0; de < 4 ;de++){  //positive de direction
  if  (de%4 != n){

for (int ro=0; ro < 8 ; ro++){  //all ro direction
if  (ro%4 != n && ro%4 != de%4)
   {  
      for (nu=0; nu < 8 ; nu++){ // 8 direction
	//for (nu=0; nu < 4 ; nu++){ // positive direction
if  (nu%4 != n && nu%4 != ro%4 && nu%4 != de%4)
   {  
     for (i = 0; i < 4 ; i++) coord[i] = x[i]; 

        Staple7_NP (UnitMatrix ,stp7_2 , coord, n, nu, ro,de, v, u, nflush_g  );

     /*   //U_de(x-mu)~ and transfer it to  x+de+mu
     CoordNN( n+4 ); //
     coord[n]= coord_nn[n];// x-mu
      nn = ( de - 1 + 4 ) % 4;
      v = u + SITE_LEN *(LexGauge( coord ))+ MATRIX_SIZE * nn;
      mDotMEqual( stp5_0, v, UnitMatrix  ) ; 

     off_node = CoordNN( de );// x-mu -> x+mu+de
     coord[de]= coord_nn[de];// x-mu+de
     TransfM( off_node,  nflush_g, stp5_0, mtmp,  de); //transfer stp5_0 to  x-mu+de

     CoordNN( n ); // 
     coord[n]= coord_nn[n];//  x-mu+de -> x+de
 
 
 
 if( ro<4)
    Staple5_NP(stp5_0,stp5_4, coord, n,nu,ro, v, u, nflush_g);
else    
    Staple5_NN(stp5_0,stp5_4, coord, n,nu, ro%4, v, u, nflush_g );



   //transfer it to  x 
     off_node = CoordNN( de +4 );//  x+de ->x
     coord[de]= coord_nn[de];//
     TransfP( off_node,  nflush_g, stp5_4, mtmp, de);

   
      nn = ( de - 1 + 4 ) % 4;
      v = u + SITE_LEN *(LexGauge( coord ))+ MATRIX_SIZE * nn;
      DaggerM(stp5_1,v);
       mDotMEqual(stp7_2, stp5_1,  stp5_4) ;

     */
 
       for ( i = 0; i < 18; i++ ) {
	 stp7_f[i] += stp7_2[i];}
}
 }}} }}


 

for (int de=0; de < 4 ;de++){  //negative de direction
  if  (de%4 != n){

for (int ro=0; ro < 8 ; ro++){  //all ro direction
if  (ro%4 != n && ro%4 != de%4)
   {
     //for (nu=0; nu < 4 ; nu++){ // positive direction  
 for (nu=0; nu < 8 ; nu++){ // 8 direction
if  (nu%4 != n && nu%4 != ro%4 && nu%4 != de%4)
   {  
     for (i = 0; i < 4 ; i++) coord[i] = x[i]; 

     Staple7_NN (UnitMatrix ,stp7_2 , coord, n, nu, ro,de, v, u, nflush_g  );
 
       for ( i = 0; i < 18; i++ ) {
	 stp7_f[i] += stp7_2[i];}
}
 }}} }}






 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 // add the fat link term to one link
 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
    for ( i = 0; i < 18; i++ ) {

      w[i]*= c1;
      //w[i] = w[i] + c3 * stp3_5[i]+ c5 * stp5_f[i] + c7 * stp7_f[i] ;
      w[i] = w[i] - c3 * stp3_5[i]- c5 * stp5_f[i] - c6*stp6_f[i] - c7 * stp7_f[i] ; //test
      /*
      w[i]+=stp3_5[i];
        w[i]+=stp5_f[i];
	 w[i]+=stp7_f[i];
*/
           }











 

 













     }}} }//end of all the x-for loop
   }//end of the for-n-loop //all negative direction

#endif // SIMUL 

  int fd;
  char buf[200];
  FILE *fp;
#ifndef SIMUL_U
#if 0
  fp=Fopen("uc_l.h","w");
  for(j=0;j<2;j++){
    Fprintf(fp,"IFloat uc_l%d[] LOCATE(\"edramnormal\") = {\n",j);
    for(i=0;i< MATRIX_SIZE * ((local_chi+ local_chi_3)/2);i++){
      Fprintf(fp,"%0.4e, ",*(uc_l[j]+i));
      if( i%6 == 5){
        Fprintf(fp,"\n");
      }
    }
    Fprintf(fp,"\n};\n"); 
  }
  Fclose(fp);
#endif
#endif

#ifndef SIMUL_AGG

  gauge_agg *temp = (gauge_agg *)
    smalloc(12*vol*sizeof(gauge_agg),"temp",fname,cname);
  int *num_ind = (int *)smalloc(6*vol*sizeof(int),"num_ind",fname,cname);
  int src;
  for(j=0;j<2;j++){
    for(i=0;i<vol;i++) num_ind[i]=0;
    for(i=0;i< ((local_chi+ local_chi_3)/2);i++){
      src = (unsigned long)chi_l[j][2*i];
      if (src%(VECT_LEN*sizeof(IFloat))!=0){
        ERR.General(cname,fname,"src = %d\n",src);
        exit(1);
      }
      src = src/(VECT_LEN*sizeof(IFloat));
      if(src > vol/2) {
        ERR.General(cname,fname,"src[%d](%d) > vol/2\n",i,src);
        exit(1);
      }
      temp[src*NUM_DIR*2+num_ind[src]].src = (unsigned long)chi_l[j][2*i];
      temp[src*NUM_DIR*2+num_ind[src]].dest= (unsigned long)chi_l[j][2*i+1];
      for(k=0;k<MATRIX_SIZE;k++){
        temp[src*NUM_DIR*2+num_ind[src]].mat[k] = uc_l[j][i*MATRIX_SIZE+k];
      }
      num_ind[src]++;
    }
    int index = 0;
    int div = 2;
    for( n = 0; n*div<NUM_DIR*2;n++)
    for(i=0;i<vol/2;i++){
      for(m=0;m<div;m++){
        if(num_ind[i] > n*div+m) {
           uc_l_agg[j][index].src = temp[i*NUM_DIR*2+n*div+m].src;
           uc_l_agg[j][index].dest = temp[i*NUM_DIR*2+n*div+m].dest;
           for(k=0;k<18;k++)
              uc_l_agg[j][index].mat[k] = temp[i*NUM_DIR*2+n*div+m].mat[k];
           index++;
        }
      }
    }
  }


#if 0
  struct gauge_agg *agg_p;
  sprintf(buf,"uc_l_agg.h");
  fd = open(buf,O_CREAT|O_TRUNC|O_RDWR,00644);
  for(j=0;j<2;j++){
    sprintf(buf,"struct gauge_agg uc_l_agg%d[] LOCATE(\"edramtransient\") = {\n",j);
    write(fd,buf,strlen(buf));
    for(i=0;i< ((local_chi+ local_chi_3)/2);i++){
      agg_p = &(uc_l_agg[j][i]);
      sprintf(buf,"{%d,%d,{\n",agg_p->src,agg_p->dest);
      write(fd,buf,strlen(buf));
      for(k=0;k<18;k++){
        sprintf(buf,"%0.8e",agg_p->mat[k]);
        write(fd,buf,strlen(buf));
	if(k!=17){
        sprintf(buf,", ");
        write(fd,buf,strlen(buf));
	}
        if( k%6 == 5){
          sprintf(buf,"\n");
          write(fd,buf,strlen(buf));
        }
      }
      sprintf(buf,"}},\n");
      write(fd,buf,strlen(buf));
    }
    sprintf(buf,"\n};\n"); 
    write(fd,buf,strlen(buf));
  }
  close(fd);
#endif
#else

#if 0
  sprintf(buf,"uc_l_agg.h");
  fd = open(buf,O_CREAT|O_TRUNC|O_RDWR,00644);
  for(j=0;j<2;j++){
    sprintf(buf,"struct gauge_agg uc_l_agg%d[] LOCATE(\"edramtransient\") = {\n",j);
    write(fd,buf,strlen(buf));
    for(i=0;i< ((local_chi+ local_chi_3)/2);i++){
      sprintf(buf,"{%d,%d,{\n",*(chi_l[j]+2*i),*(chi_l[j]+2*i+1));
      write(fd,buf,strlen(buf));
      uc_l_agg[j][i].src = (unsigned long)chi_l[j][2*i];
      uc_l_agg[j][i].dest = (unsigned long)chi_l[j][2*i+1];
      for(k=0;k<18;k++){
        sprintf(buf,"%0.8e",*(uc_l[j]+i*18+k));
        write(fd,buf,strlen(buf));
        uc_l_agg[j][i].mat[k] = uc_l[j][i*18+k];
	if(k!=17){
        sprintf(buf,", ");
        write(fd,buf,strlen(buf));
	}
        if( k%6 == 5){
          sprintf(buf,"\n");
          write(fd,buf,strlen(buf));
        }
      }
      sprintf(buf,"}},\n");
      write(fd,buf,strlen(buf));
    }
    sprintf(buf,"\n};\n"); 
    write(fd,buf,strlen(buf));
  }
  close(fd);
#endif

#endif //SIMUL_AGG

#ifndef SIMUL_U
#if 0
  fp = Fopen("uc_nl.h","w");
  for(j=0;j<2;j++){
    Fprintf(fp,"IFloat uc_nl%d[] LOCATE(\"edramnormal\") = {\n",j); 
    for(i=0;i< MATRIX_SIZE * ((non_local_chi+non_local_chi_3 )/2);i++){
      Fprintf(fp,"%0.4e, ",*(uc_nl[j]+i));
      if( i%6 == 5){
        Fprintf(fp,"\n");
      }
    }
    Fprintf(fp,"\n};\n"); 
  }
  Fclose(fp);
#endif
#endif

#ifndef SIMUL_AGG

#if 1
  int max = 12;
  for(j=0;j<2;j++){
    for(i=0;i<non_local_chi*3/2;i++) num_ind[i]=0;
    for(i=0;i< ((non_local_chi+ non_local_chi_3)/2);i++){
      src = (unsigned long)chi_nl[j][2*i];
      if (src%(VECT_LEN*sizeof(IFloat))!=0){
        ERR.General(cname,fname,"src = %d\n",src);
        exit(1);
      }
      src = src/(VECT_LEN*sizeof(IFloat));
      if(src > non_local_chi*3/2) {
        ERR.General(cname,fname,"src(%d) > non_local_chi*3\n",src);
        exit(1);
      }
      temp[src*max+num_ind[src]].src = (unsigned long)chi_nl[j][2*i];
      temp[src*max+num_ind[src]].dest= (unsigned long)chi_nl[j][2*i+1];
      for(k=0;k<18;k++)
        temp[src*max+num_ind[src]].mat[k] = uc_nl[j][i*18+k];
      num_ind[src]++;
      if(num_ind[src]>max){
        ERR.General(cname,fname,"num_ind[%d](%d) > %d \n",src, num_ind[src],max);
        exit(1);
      }
    }

    int index = 0;
    int div = 12;
    for( n=0; n*div<max; n++)
    for(i=0;i<non_local_chi*3/2;i++){
      for(m=0;m<div;m++){
        if(num_ind[i] > n*div+m) {
           uc_nl_agg[j][index] = temp[i*max+n*div+m];
           index++;
        }
      }
    }

  }
  sfree(temp,"temp",fname,cname);
  sfree(num_ind,"num_ind",fname,cname);

#if 0
  sprintf(buf,"uc_nl_agg.h");
  fd = open(buf,O_CREAT|O_TRUNC|O_RDWR,00644);
  for(j=0;j<2;j++){
    sprintf(buf,"struct gauge_agg uc_nl_agg%d[] LOCATE(\"edramtransient\") = {\n",j);
    write(fd,buf,strlen(buf));
    for(i=0;i< ((non_local_chi+ non_local_chi_3)/2);i++){
      agg_p = &(uc_nl_agg[j][i]);
      sprintf(buf,"{%d,%d,{\n",agg_p->src,agg_p->dest);
      write(fd,buf,strlen(buf));
      for(k=0;k<18;k++){
        sprintf(buf,"%0.8e",agg_p->mat[k]);
        write(fd,buf,strlen(buf));
	if(k!=17){
        sprintf(buf,", ");
        write(fd,buf,strlen(buf));
	}
        if( k%6 == 5){
          sprintf(buf,"\n");
          write(fd,buf,strlen(buf));
        }
      }
      sprintf(buf,"}},\n");
      write(fd,buf,strlen(buf));
    }
    sprintf(buf,"\n};\n"); 
    write(fd,buf,strlen(buf));
  }
  close(fd);
#endif
#else
#if 0
  sprintf(buf,"uc_nl_agg.h");
  fd = open(buf,O_CREAT|O_TRUNC|O_RDWR,00644);
  for(j=0;j<2;j++){
    sprintf(buf,"struct gauge_agg uc_nl_agg%d[] LOCATE(\"edramtransient\") = {\n",j);
    write(fd,buf,strlen(buf));
    for(i=0;i< ((local_chi+ local_chi_3)/2);i++){
      sprintf(buf,"{%d,%d,{\n",*(chi_nl[j]+2*i),*(chi_nl[j]+2*i+1));
      uc_nl_agg[j][i].src = (unsigned long)chi_nl[j][2*i];
      uc_nl_agg[j][i].dest = (unsigned long)chi_nl[j][2*i+1];
      write(fd,buf,strlen(buf));
      for(k=0;k<18;k++){
        sprintf(buf,"%0.8e",*(uc_nl[j]+i*18+k));
        write(fd,buf,strlen(buf));
        uc_nl_agg[j][i].mat[k] = uc_nl[j][i*18+k];
	if(k!=17){
        sprintf(buf,", ");
        write(fd,buf,strlen(buf));
	}
        if( k%6 == 5){
          sprintf(buf,"\n");
          write(fd,buf,strlen(buf));
        }
      }
      sprintf(buf,"}},\n");
      write(fd,buf,strlen(buf));
    }
    sprintf(buf,"\n};\n"); 
    write(fd,buf,strlen(buf));
  }
  close(fd);
#endif
#endif

#endif // SIMUL_AGG

  VRB.FuncEnd(cname,fname);
}

extern "C"
void asqtad_destroy_dirac_buf_g(void)
{
  int i;
  char *cname = "";
  char *fname = "asqtad_destroy_dirac_buf_g()";
  VRB.Func(cname,fname);
  for ( i = 0; i < 2; i++){
#ifndef SIMUL_U
    sfree(uc_l[i],"uc_l[i]",fname,cname);
    sfree(uc_nl[i],"uc_nl[i]",fname,cname);
#endif
    sfree(uc_l_agg[i],"uc_l_agg[i]",fname,cname);
    sfree(uc_nl_agg[i],"uc_nl_agg[i]",fname,cname);
  }
  VRB.FuncEnd("",fname);
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

  for ( i = 0; i < 4; i++ )
    coord_nn[i] = coord[i];

  m = nn % 4;
  for ( i = 0; i < 4; i++ ){
        coord_nn [ i ] = coord[i];
   }//end of i-loop

     if (nn<4 ) coord_nn[m]+=1;
     else coord_nn[m]-=1;
     while (coord_nn[m] < 0) {
        coord_nn[m] += size[m] ;
         }
       coord_nn[m] %= size[m] ;
    if (nn<4){ if (coord_nn[m] != (coord[m]+1)) off_node = 1;}
       else {if (coord_nn[m] != (coord[m]-1)) off_node = 1;}
  return off_node;
}

static int CoordkNN( int nn, int k )
{
  int i;
  int off_node = 0;
  int m;

  for ( i = 0; i < 4; i++ )
    coord_knn[i] = coord[i];

  m = nn % 4;
  for ( i = 0; i < 4; i++ ){
        coord_knn [ i ] = coord[i];
   }//end of i-loop

     if (nn<4 ) coord_knn[m]+=k;
     else coord_knn[m]-=k;
     while (coord_knn[m] < 0) {
        coord_knn[m] += size[m] ;
         }
       coord_knn[m] %= size[m] ;
 
  if (nn<4 ){ if (coord_knn[m] != (coord[m]+k)) off_node = 1;}
       else {if (coord_knn[m] != (coord[m]-k)) off_node = 1;}
  return off_node;
}

/*static int Coord3NN( int nn )
{
  int i;
  int off_node = 0;
  int m;

  m = nn % 4;
  for ( i = 0; i < 4; i++ ){
        coord_3nn [ i ] = coord[i];
   }//end of i-loop

   if (nn<4 ) coord_3nn[m]+=3;
     else coord_3nn[m]-=3;

     while (coord_3nn[m] < 0) {
         coord_3nn[m] += size[m] ;
         }
       coord_3nn[m] %= size[m] ;

    if (coord_3nn[m] != (coord[m]+3) ||coord_3nn[m] != (coord[m]-3) ) off_node = 1;
 return off_node;
}
*/
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

extern "C" void dirac_comm_assert()
{
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
void asqtad_dirac(IFloat* b, IFloat* a, int a_odd, int add_flag)
{
//  int i,j,k,c = (unsigned long) chi_off_node_total;
  int i,j,k;
  long c = (long) a;
  int odd = 1 - a_odd;
  //-----------------------------------------------------------------
  //  Transfer chi's on faces.  
  //-----------------------------------------------------------------

#if 0
  copy_buffer_cpp(countM[0][a_odd], (long)a, (long)Tbuffer[0][0], (long)ToffsetM[0][a_odd]);
  copy_buffer_cpp(countP[0][a_odd], (long)a, (long)Tbuffer[0][1], (long)ToffsetP[0][a_odd]);
  if(size[0]>2) {
  copy_buffer_cpp(countM[1][a_odd], (long)a, (long)Tbuffer[1][0], (long)ToffsetM[1][a_odd]);
  copy_buffer_cpp(countP[1][a_odd], (long)a, (long)Tbuffer[1][1], (long)ToffsetP[1][a_odd]);
  if(size[0]>4) {
  copy_buffer_cpp(countM[2][a_odd], (long)a, (long)Tbuffer[2][0], (long)ToffsetM[2][a_odd]);
  copy_buffer_cpp(countP[2][a_odd], (long)a, (long)Tbuffer[2][1], (long)ToffsetP[2][a_odd]);
  }
  }
#endif

  //make sure spinor field is in main memory before starting transfers

#if 0
  for(i=1;i<4;i++)
  if(size[i] >2 ){
    SCUarg[i]->Addr( a + Xoffset[2][i]);
    SCUarg[i+4]->Addr( a + Xoffset[2][i+4]);
  }

  save_reg((long)intreg, (long)dreg);
 //*************

  void *addr[2];
  if(split) {
    for(i=1;i<4;i++){
      addr[0] = (void *)(a+Xoffset[0][i]);
      addr[1] = (void *)(a+Xoffset[1][i]);
      SCUarg_1[i]->Addr( addr,1);
      SCUarg_2[i]->Addr( addr+1,1);
      addr[0] = (void *)(a+Xoffset[0][i+4]);
      addr[1] = (void *)(a+Xoffset[1][i+4]);
      SCUarg_1[i+4]->Addr( addr,1);
      SCUarg_2[i+4]->Addr( addr+1,1);
    }
  } else {
    for(i=1;i<4;i++){
      addr[0] = (void *)(a+Xoffset[0][i]);
      addr[1] = (void *)(a+Xoffset[1][i]);
      SCUarg_1[i]->Addr( addr,2);
      addr[0] = (void *)(a+Xoffset[0][i+4]);
      addr[1] = (void *)(a+Xoffset[1][i+4]);
      SCUarg_1[i+4]->Addr( addr,2);
    }
  }

  SCUmulti_1->StartTrans();

  save_reg((long)intreg, (long)dreg);
#endif

  //-----------------------------------------------------------------
  //do first local computations
  //-----------------------------------------------------------------
  IFloat *fp0,*fp1, *uu;
  int *ch;


  dirac_cmv_jcw_agg_cpp( (local_chi + local_chi_3)/2, (long)0, (long)uc_l_agg[odd],
	       (long)a, (long)tmpfrm);


  //-----------------------------------------------------------------
  // check to see if transfers are done and start another transfer
  //-----------------------------------------------------------------



  //-----------------------------------------------------------------
  //do the computations involving "chi" non-local spinors
  //-----------------------------------------------------------------

  dirac_cmv_jcw_agg_cpp( non_local_chi_3/2, (long)0, (long)&(uc_nl_agg[odd][0]), (long)c, (long)tmpfrm);

 

  //-----------------------------------------------------------------
  //do the computations involving chi3 non-local spinors
  //----------------------------------------------------------------

  dirac_cmv_jcw_agg_cpp( non_local_chi/2, (long)0, (long)&(uc_nl_agg[odd][non_local_chi_3/2]), (long)c, (long)tmpfrm);

  //printf ("the computations involving chi3 non-local spinors done \n");
  // check to see if transfers are done
  //-----------------------------------------------------------------
  //do the sum of 16 temporary vectors at each lattice site
  //              ^^^ change must be made in  dirac_sum**
  //-----------------------------------------------------------------

  if ( add_flag == 0){
    dirac_sum2_cpp( vol/2, (long)0, (long)tmpfrm, (long)b);
  }
  else{

    dirac_sum2_acc_cpp( vol/2, (long)0, (long)tmpfrm, (long)b);

  }

}

CPS_END_NAMESPACE




