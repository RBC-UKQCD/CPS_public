
//    12/02 HueyWen Lin, Chulwoo Jung
//
//   Asqtad Dirac operator for QCDOC. Communications and computations
//   are overlapped.
//   Uses functions implemented by CJ for QCDOC
//
//
//-------------------------------------------------------------------
#include <config.h>
#include <stdio.h>
CPS_END_NAMESPACE
#include <util/gjp.h>
#include <util/time.h>
#include <comms/scu.h>
#include <comms/glb.h>
#include <util/lattice.h>
#include <util/dirac_op.h>
#include <util/vector.h>
#include <util/pt.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <time.h>
#include <math.h>
#include <qcdocos.h>
#include <qalloc.h>
#include <ppc_lib.h>
CPS_START_NAMESPACE

const char *uc_l_filename = CWDPREFIX("uc_l.h");
const char *uc_nl_filename = CWDPREFIX("uc_nl.h");
const char *Toffset_filename = CWDPREFIX("Toffset.h");
const char *uc_l_agg_filename = CWDPREFIX("uc_l_agg.h");
const char *uc_nl_agg_filename = CWDPREFIX("uc_nl_agg.h");
const char *chi_l_filename = CWDPREFIX("chi_l.h");
const char *chi_nl_filename = CWDPREFIX("chi_nl.h");
#define IND_AGG //aggregate index for uc and chi

/*****************************************************************
 SIMUL switched on/off the hack CJ put in to help speed up the
simulation. if SIMUL is undefined, program writes the temporary arraies
for dirac operator. If SIMUL is defined, it will include (pre-calcuated)
temporary arraies and skip the generation of arraies.
******************************************************************/

#define USE_NEW
#undef CPP
/****************************************************************
CPP is a switch for using C++ routine for dirac_cmv.
*****************************************************************/

void dirac_sum_acc_cpp(int s, long chi, long tmpfrm, long b);
void dirac_cmv_jcw_agg_cpp( int sites, long chi, long u,long a, long tmpfrm);
void dirac_sum2_64_cpp( int sites, long chi, long tmpfrm,long b);

extern "C" void save_reg(long intbuf, long dbuf );
extern "C" void restore_reg(long intbuf, long dbuf );
extern "C" void copy_buffer(int n, long src, long dest, long ptable);
extern "C" void dirac_cmv_jcw_agg( int sites, long chi, long u,long a, long tmpfrm);
extern "C" void dirac_cmv_agg_4( int sites, long chi, long u,long a, long tmpfrm);
extern "C" void dirac_sum2_64( int sites, long chi, long tmpfrm,long b);
extern "C" void flush_cache_spinor(int nflush, long flush_buffer);
extern "C" void flush_cache(int nflush, long flush_buffer);

enum{VECT_LEN=6, VECT_LEN2=8, MATRIX_SIZE=18, SITE_LEN=72, NUM_DIR=8};
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
static int split; 
static int isplit; 
static int coord[4];
static int coord_nn[4];
static int coord_knn[4];
static int vol;
static int non_local_chi_3[4];
static int non_local_chi;
static int local_chi;
static int local_chi_3;
static int nflush;
static int odd_num= 0;

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

SCUDirArgIR SCUarg[2][2*NUM_DIR];
SCUDirArgMulti SCUmulti[2];

SCUDirArgIR SCUarg_1[2][2*NUM_DIR];
SCUDirArgMulti SCUmulti_1[2];

SCUDMAInst *SCUDMAarg_p[2][NUM_DIR*4];
SCUDMAInst SCUDMAarg[2][NUM_DIR*4];

SCUDirArgIR SCUarg_2[2][2*NUM_DIR];
SCUDirArgMulti SCUmulti_2[2];

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

static int NotParallel( int mu, int nu){
	int dif = mu-nu;
	if (dif==0 ||dif==-4||dif==4)return 0;
	else return 1;
}

static IFloat * gauge_field_addr;

int k;

static Fasqtad *lat_pt;
void set_pt (Fasqtad *lat)
{
  lat_pt = lat;
}
//-------------------------------------------------------------------
//  Called by the lattice constructor
//  Fermion initializations: pointer tables 
//-------------------------------------------------------------------
extern "C"
void asqtad_dirac_init(const void * gauge_u )
{
  gauge_field_addr = ( IFloat * ) gauge_u;
  int r,c,i,j,m,n,nn;
  int blklen[NUM_DIR/2];
  int numblk[NUM_DIR/2];
  int stride[NUM_DIR/2];
  int local_count[2];
  int non_local_count[2];
  int local_count_3[2];
 int non_local_count_3[3][2];
 int off_node;
  int x[NUM_DIR/2];
  char buf[200];
  int fd;
  int scu_irs[2][3];
  char *cname = "DiracOpAsqtad";
  char *fname = "asqtad_dirac_init(const void *gauge)";

  scu_irs[0][0]=IR_5;
  scu_irs[0][1]=IR_6;
  scu_irs[0][2]=IR_7;
  scu_irs[1][0]=IR_8;
  scu_irs[1][1]=IR_9;
  scu_irs[1][2]=IR_10;


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
  split =  (  vol>64 ? 0 : 1 ); 
//  split =0;

  non_local_chi = 2*(size[0]*size[1]*size[2] + size[1]*size[2]*size[3]+
    size[2]*size[3]*size[0] + size[3]*size[0]*size[1]);

  non_local_chi_3[0] = 0;
  non_local_chi_3[1] = 0;
  for(i=0;i<4;i++){
    if(size[i]>2)
      non_local_chi_3[1] += 2*vol/size[i];
  }
  non_local_chi_3[2]  = non_local_chi_3[1]+non_local_chi;
  non_local_chi_3[3]  = non_local_chi_3[2]+non_local_chi;


  local_chi = NUM_DIR*vol - non_local_chi;
  local_chi_3 = NUM_DIR*vol - (non_local_chi_3[3]);
  isplit = (non_local_chi_3[1]+non_local_chi)/2;

//  printf("local_chi=%d local_chi_3=%d\n",local_chi,local_chi_3);
//  printf("non_local_chi=%d non_local_chi_3[1]=%d\n",non_local_chi,non_local_chi_3[1]);
//  printf("non_local_chi_3[2]=%d non_local_chi_3[3]=%d\n",non_local_chi_3[2],non_local_chi_3[3]);
  //-------------------------------------------------------------
  // flush_cache_spinor() function will flush 192 bytes * nflush 
  //-------------------------------------------------------------
  nflush = vol/8;

  tmpfrm = (IFloat *) qalloc (QFAST|QCOMMS,NUM_DIR*2 * vol/2 * VECT_LEN2 * sizeof(IFloat));

  if(tmpfrm == NULL){ 
    tmpfrm = (IFloat *) qalloc (QCOMMS,NUM_DIR*2 * vol/2 * VECT_LEN2 * sizeof(IFloat));
    printf("tmpfrm is allocated at (%p),length 0x%x \n",tmpfrm,NUM_DIR*vol*VECT_LEN2*sizeof(IFloat));
  }
  if(tmpfrm == 0) 
    ERR.Pointer(cname,fname, "tmpfrm");


  //-----------------------------------------------------------------
  //  Allocate 8 receive buffers for off-node vectors
  //-----------------------------------------------------------------
  
  if(vol>1024) chi_off_node_total=NULL;
    else
    chi_off_node_total = ( IFloat * ) qalloc(QFAST|QCOMMS, 3*non_local_chi*
      VECT_LEN * sizeof( IFloat ) / 2 );
  if(chi_off_node_total == NULL){ 
    chi_off_node_total = ( IFloat * ) qalloc(QCOMMS, 3*non_local_chi*
      VECT_LEN * sizeof( IFloat ) / 2 );
    printf("chi_off_node_total is allocated at DDR (%p)\n",chi_off_node_total);
  }
    if(chi_off_node_total == 0)
      ERR.Pointer(cname,fname, "chi_off_node_total");
 for ( j= 0; j < 3; j++ ){
    chi_off_node[j][0] = &(chi_off_node_total[ non_local_chi*j* VECT_LEN/2 ] ); 
//    fprintf(stderr,"chi_off_node[%d][0] =%d\n",j,(chi_off_node[j][0]-chi_off_node_total)/VECT_LEN);
    chi_off_node_p[j][0] = (IFloat *)(sizeof (IFloat)*non_local_chi*j* VECT_LEN/2 );
  for ( i = 1; i < NUM_DIR; i++ ){
    chi_off_node[j][i] = chi_off_node[j][i-1]+vol/(2*size[(i-1)%4])*VECT_LEN;
//    fprintf(stderr,"chi_off_node[%d][%d] =%d\n",j,i,(chi_off_node[j][i]-chi_off_node_total)/VECT_LEN);
    chi_off_node_p[j][i] = chi_off_node_p[j][i-1]+vol/(2*size[(i-1)%4])*VECT_LEN;
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

  for ( i = 0; i < 2; i++ ){
    chi_l[i] = ( IFloat ** ) smalloc((2*(local_chi+local_chi_3)/2)*sizeof(IFloat *));
   if(chi_l[i] == NULL)
	ERR.Pointer(cname,fname, "chi_l[i]");

    chi_nl[i] = (IFloat ** ) smalloc((2* (non_local_chi + non_local_chi_3[3])/2)*sizeof(IFloat *));

      if(chi_nl[i] == NULL)
	ERR.Pointer(cname,fname, "chi_nl[i]");
  }
  

  for ( i = 0; i < 2; i++){
    local_count[i] = 0;
    non_local_count[i] = 0;
    local_count_3[i] = 0;
    for( n = 0; n < 3; n++ )  non_local_count_3[n][i] = 0;

  }

  for ( k = 0; k < 3; k++ ) {
  for ( i = 0; i < 2; i++ ) {
    Tbuffer[k][i] = (IFloat *) qalloc (QFAST|QNONCACHE, size[1] * size[2] * size[3] * VECT_LEN * sizeof( IFloat ) / 2);
   if(Tbuffer[k][i] == NULL)
	ERR.Pointer(cname,fname, "Tbuffer[i][j]");
    ToffsetP[k][i] = ( int * ) qalloc (0,  size[1] * size[2] * size[3] *  sizeof( int ) / 2 );
   if(ToffsetP[k][i] == NULL)
	ERR.Pointer(cname,fname, "TOffsetP[i][j]");
    ToffsetM[k][i] = ( int * ) qalloc (0,  size[1] * size[2] * size[3] *  sizeof( int ) / 2 );
   if(ToffsetM[k][i] == NULL)
	ERR.Pointer(cname,fname, "TOffsetM[i][j]");
    countP[k][i] = 0;
    countM[k][i] = 0;
  }
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
	     		= chi_off_node_p[0][n] + VECT_LEN * ( LexSurface( coord_nn, n%4 ) / 2 );
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
//  printf("chi_l=0x%x\nchi_l[0]=%x\n",chi_l,chi_l[0]);
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
	      //pointer to source field
	      *( chi_l[ odd ]  +  2 * local_count[ odd ] )
		= ( IFloat * ) ( VECT_LEN * ( LexVector( coord_nn ) / 2 ) * sizeof(IFloat));
	      // pointer to temporary field where U*chi is stored
	      *( chi_l[ odd ]  +  2 * local_count[ odd ] + 1) =
					( IFloat * ) ( VECT_LEN2 * (LexVector( coord ) / 2 *2*NUM_DIR+n) * sizeof(IFloat));
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
	 *( chi_nl[ odd ] + 2 * non_local_count_3[j][odd]+non_local_chi+non_local_chi_3[j] )
            =chi_off_node_p[j][n]+ VECT_LEN * ( LexSurface( coord_knn,n%4)/2 ); 
        *( chi_nl[odd] + 2 * non_local_count_3[j][odd] + non_local_chi+non_local_chi_3[j]+1 )
 = ( IFloat * ) ( VECT_LEN2 * (LexVector( coord ) / 2*2*NUM_DIR+n+8 )  * sizeof(IFloat));
//	printf("%d %d %d %d %d chi_nl[%d][%d]=%d j=%d\n",x[0],x[1],x[2],x[3],n,odd,non_local_count_3[j][odd]+(j+1)*non_local_chi/2,LexSurface(coord_knn,n%4)/2,j);
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
	      //pointer to source field
	      *( chi_l[ odd ]  +  2 * local_count_3[ odd ] + local_chi )
		= ( IFloat * ) ( VECT_LEN * ( LexVector( coord_knn ) / 2 )
		                 * sizeof(IFloat));
	      // pointer to temporary field where U*chi is stored
	      *( chi_l[ odd ]  +  2 * local_count_3[ odd ] + local_chi + 1) =
	  	( IFloat * ) ( VECT_LEN2 * (LexVector( coord ) / 2 *2*NUM_DIR+n+8) * sizeof(IFloat));
	      local_count_3[odd]++;
	    }  //else-on_node case


	  } // for x[0] loop
	}
      }
    }// for x[3] loop
  }//for n loop
  FILE *fp;

#if 0
  fp=fopen(chi_l_filename,"w");
  for(j=0;j<2;j++){
    fprintf(fp,"IFloat * chi_l%d[] LOCATE(\"edramtransient\") = {\n",j);
    fprintf(fp," (IFloat *) %d",*(chi_l[j]));
    for(i=1;i< 2*((local_chi+local_chi_3)/2);i++){
      fprintf(fp,",\n (IFloat *) %d",*(chi_l[j]+i));
    }
    fprintf(fp,"\n};\n");
  }
  fclose(fp);

  fp=fopen(chi_nl_filename,"w");
  for(j=0;j<2;j++){
    fprintf(fp,"IFloat * chi_nl%d[] LOCATE(\"edramtransient\") = {\n",j);
    fprintf(fp," (IFloat *) %d",*(chi_nl[j]));
    for(i=1;i< 2*((non_local_chi+ non_local_chi_3[3])/2);i++){
      fprintf(fp,",\n (IFloat *) %d",*(chi_nl[j]+i));
    }
    fprintf(fp,"\n};\n");
  }

  fclose(fp);
#endif

#if 0
  for(i=0;i<2;i++){
	printf("local_count[%d]=%d non_local_count[%d]=%d\n",i,local_count[i],i,non_local_count[i]);
	printf("local_count_3[%d]=%d non_local_count_3[0][%d]=%d\n",i,local_count_3[i],i,non_local_count_3[0][i]);
  	for(j=1;j<3;j++)
	printf("non_local_count_3[%d][%d]=%d\n",j,i,non_local_count_3[j][i]);
  }
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
  stride[2] = (VECT_LEN * size[0] * size[1] * ( size[2] - 1 ) / 2)* sizeof(IFloat) ;
  stride[3] = 0;

  //-------------------------------------------------------------------
  //  Calculate offsets for T transfers done one word at a time.
  //  We have plus (P) transfers for both the even and odd
  //  checkerboards.  Same for minus (M) transfers.
  //-------------------------------------------------------------------
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

#if 0
  fp= fopen(Toffset_filename,"w");
  for ( k = 0; k < 3; k++ ) 
  for ( i = 0; i < 2; i++ ) {
    fprintf(fp,"countP[%d][%d]=%d\n",k,i,countP[k][i]);
    fprintf(fp,"int ToffsetP%d%d[] LOCATE(\"edramnormal\") = {\n",k,i);
    for( j = 0;j<vol3;j++){
      fprintf(fp, "%d,\n",ToffsetP[k][i][j]);
    }
    fprintf(fp,"};\n");
    fprintf(fp,"countM[%d][%d]=%d\n",k,i,countM[k][i]);
    fprintf(fp,"int ToffsetM%d%d[] LOCATE(\"edramnormal\") = {\n",k,i);
    for( j = 0;j<vol3;j++){
      fprintf(fp, "%d,\n",ToffsetM[k][i][j]);
    }
    fprintf(fp,"};\n");
  }
  fclose(fp);
#endif


  //-------------------------------------------------------------------
  //  Index i says data has been received from TP, XP, YP, ZP, TM, XM,
  //  YM, ZM
  //-------------------------------------------------------------------

#if 0
	for(int kk = 0;kk<3;kk++){
	printf("Tbuffer[%d][0]=%p\n",kk,Tbuffer[kk][0]);
	printf("Tbuffer[%d][1]=%p\n",kk,Tbuffer[kk][1]);
	}
#endif

  for( odd=0;odd<2;odd++){
    for ( i = 0; i < NUM_DIR; i++ ) {
      j = i % (NUM_DIR/2);
  
        SCUarg[odd][i + NUM_DIR].Init(chi_off_node[2][i], scudir[i], SCU_REC,
  		    VECT_LEN * sizeof(IFloat) * vol / ( 2 * size[j] ), 1, 0, scu_irs[odd][2]);
        SCUarg[odd][i + NUM_DIR].Assert();
  
      SCUDMAarg_p[odd][(i+NUM_DIR)*2]  = new SCUDMAInst;
  //      printf("SCUDMAarg_p[%d]=%p\n",(i+NUM_DIR)*2,SCUDMAarg_p[(i+NUM_DIR)*2]);
      SCUDMAarg_p[odd][(i+NUM_DIR)*2] ->Init(chi_off_node[0][i],
        VECT_LEN * sizeof(IFloat) * vol / ( 2 * size[j] ), 1, 0);
      SCUDMAarg_p[odd][(i+NUM_DIR)*2+1]  = new SCUDMAInst;
      SCUDMAarg_p[odd][(i+NUM_DIR)*2+1] ->Init(chi_off_node[1][i],
        VECT_LEN * sizeof(IFloat) * vol / ( 2 * size[j] ), 1, 0);
  
      if( split ){
        SCUarg_1[odd][i + NUM_DIR].Init(scudir[i],SCU_REC, &SCUDMAarg_p[odd][(i+NUM_DIR)*2],1, scu_irs[odd][0]);
        SCUarg_1[odd][i + NUM_DIR].Assert();
        SCUarg_2[odd][i + NUM_DIR].Init(scudir[i],SCU_REC, &SCUDMAarg_p[odd][(i+NUM_DIR)*2+1],1, scu_irs[odd][1]);
        SCUarg_2[odd][i + NUM_DIR].Assert();
      } else {
        SCUarg_1[odd][i + NUM_DIR].Init(scudir[i],SCU_REC, &SCUDMAarg_p[odd][(i+NUM_DIR)*2],2, scu_irs[odd][0]);
        SCUarg_1[odd][i + NUM_DIR].Assert();
      }
  
      //send arguments
      if ((i == 0) || ( i == 4)){
        if (size[j] >4)
          SCUarg[odd][i].Init (Tbuffer[2][(4 - i)/4], scudir[i], SCU_SEND,
  		       blklen[j], numblk[j], stride[j], scu_irs[odd][2] );
        else if (size[j] >2) 
          SCUarg[odd][i].Init (Tbuffer[1][i/4], scudir[i], SCU_SEND,
  		       blklen[j], numblk[j], stride[j], scu_irs[odd][2] );
        else
          SCUarg[odd][i].Init (chi_off_node[0][4-i], scudir[i], SCU_SEND,
  		       blklen[j], numblk[j], stride[j], scu_irs[odd][2] );
        SCUarg[odd][i].Assert();
  
        SCUDMAarg_p[odd][i*2] = new SCUDMAInst;
        SCUDMAarg_p[odd][i*2] ->Init(Tbuffer[0][(4 - i)/4], 
  		       blklen[j], numblk[j], stride[j]);

        SCUDMAarg_p[odd][i*2+1] = new SCUDMAInst;
        if(size[j]>2)
          SCUDMAarg_p[odd][i*2+1] ->Init(Tbuffer[1][(4 - i)/4], 
  		       blklen[j], numblk[j], stride[j]);
        else
          SCUDMAarg_p[odd][i*2+1] ->Init(Tbuffer[0][i/4], 
  		       blklen[j], numblk[j], stride[j]);
  
        if( split ){
          SCUarg_1[odd][i].Init(scudir[i],SCU_SEND,&SCUDMAarg_p[odd][i*2],1,scu_irs[odd][0]);
          SCUarg_1[odd][i].Assert();
          SCUarg_2[odd][i].Init(scudir[i],SCU_SEND,&SCUDMAarg_p[odd][i*2+1],1,scu_irs[odd][1]);
          SCUarg_2[odd][i].Assert();
        } else {
          SCUarg_1[odd][i].Init(scudir[i],SCU_SEND,&SCUDMAarg_p[odd][i*2],2,scu_irs[odd][0]);
          SCUarg_1[odd][i].Assert();
        }
      }
      else{
        if(size[j] >2)
  //
  // put in chi_off_node_total to pass communicable test. CJ
  //
          SCUarg[odd][i].Init(chi_off_node_total, scudir[i], SCU_SEND,
  		       blklen[j], numblk[j], stride[j], scu_irs[odd][2] );
        else
          SCUarg[odd][i].Init(chi_off_node[0][(i+4)%NUM_DIR], scudir[i], SCU_SEND,
             VECT_LEN * sizeof(IFloat) * vol / ( 2 * size[j]),1,0, scu_irs[odd][2] );

        SCUarg[odd][i].Assert();
  
        SCUDMAarg_p[odd][i*2] = new SCUDMAInst;
        SCUDMAarg_p[odd][i*2] ->Init(chi_off_node_total, 
  		       blklen[j], numblk[j], stride[j]);
        SCUDMAarg_p[odd][i*2+1] = new SCUDMAInst;
        SCUDMAarg_p[odd][i*2+1] ->Init(chi_off_node_total, 
  		       blklen[j], numblk[j], stride[j]);
        if( split ){
          SCUarg_1[odd][i].Init(scudir[i],SCU_SEND,&SCUDMAarg_p[odd][i*2],1,scu_irs[odd][0]);
          SCUarg_1[odd][i].Assert();
          SCUarg_2[odd][i]. Init(scudir[i],SCU_SEND,&SCUDMAarg_p[odd][i*2+1],1,scu_irs[odd][1]);
          SCUarg_2[odd][i].Assert();
        } else {
          SCUarg_1[odd][i].Init(scudir[i],SCU_SEND,&SCUDMAarg_p[odd][i*2],2,scu_irs[odd][0]);
          SCUarg_1[odd][i].Assert();
        }
      }
    }// end of NUM_DIR loop
  
    SCUDirArgIR *SCUarg_p[2*NUM_DIR];
  
//    SCUmulti[odd] = new SCUDirArgMulti;
    for(i = 0;i<2*NUM_DIR;i++) SCUarg_p[i] = &(SCUarg[odd][i]);
    SCUmulti[odd].Init(SCUarg_p, 2*NUM_DIR);
//    SCUmulti_1[odd] = new SCUDirArgMulti;
    for(i = 0;i<2*NUM_DIR;i++) SCUarg_p[i] = &(SCUarg_1[odd][i]);
    SCUmulti_1[odd].Init(SCUarg_p, 2*NUM_DIR);
    if( split ){
//      SCUmulti_2[odd] = new SCUDirArgMulti;
      for(i = 0;i<2*NUM_DIR;i++) SCUarg_p[i] = &(SCUarg_2[odd][i]);
      SCUmulti_2[odd].Init(SCUarg_p, 2*NUM_DIR);
    }
  } // end of odd loop

  //-------------------------------------------------------------------
  //  Need send offsets for various transfers.  The index for
  //  sends is TM, XM, YM, ZM, TP, XP, YP, ZP, since the
  //  transfers are indexed by the node data is received from.
  //-------------------------------------------------------------------
// printf("Xoffset=%p\n",Xoffset);
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
printf("done\n");fflush(stdout);
}

extern "C"
void asqtad_destroy_dirac_buf()
{
  int i,k;

  qfree (tmpfrm);
  qfree(chi_off_node_total);

  for ( i = 0; i < 2; i++ ) {
    sfree(chi_nl[i]);
    sfree(chi_l[i]);
    for ( k = 0; k < 3; k++ ) {
      qfree(Tbuffer[k][i]);
      qfree(ToffsetP[k][i]);
      qfree(ToffsetM[k][i]);
    }
  }
    

  for (int odd = 0; odd<2;odd++){
//    delete SCUmulti[odd];
//    delete SCUmulti_1[odd];
 //   if(split)
//    delete SCUmulti_2[odd];
    for ( i = 0; i < NUM_DIR; i++ ) {
      for(k=0;k<4;k++)
        delete SCUDMAarg_p[odd][i*4+k];
    }
  }
}
//-------------------------------------------------------------------
//  Given a lexical value for gauge fields, set the coordinates.
//  Return 1 if odd, 0 if even
//-------------------------------------------------------------------

static int SetCoord( int sg )
{
  coord[0] =   sg % size[0];
  coord[1] = ((int)( sg / size[0] )) % size[1];
  coord[2] = ((int)( sg / ( size[0] * size[1] ) )) % size[2];
  coord[3] = ((int)( sg / ( size[0] * size[1] * size[2] ) ));

  return ( coord[0] + coord[1] + coord[2] + coord[3] ) % 2;
}

//-------------------------------------------------------------------
//  Prepare a copy of gauge fields and set up 
//  related pointer for Dirac
//-------------------------------------------------------------------

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
  char *cname = "DiracOpAsqtad";
  char *fname = "asqtad_dirac_init_g()";
  VRB.Func(cname,fname);

  //--------------------------------------------------------------------
  // c1 -> one link; c2 -> 3-link; c3 -> 3-link staple; c5 -> 5-link staple;
  // c7 -> 7-link staple; c6 -> 5-link "straight" staple
  //--------------------------------------------------------------------
 
  IFloat c1 = GJP.KS_coeff();
  IFloat c2 = GJP.Naik_coeff();
  IFloat c3 = GJP.staple3_coeff();
  IFloat c5 = GJP.staple5_coeff();
  IFloat c7 = GJP.staple7_coeff();
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

  //-----------------------------------------------------------------
  //  Allocate space for two copies of the gauge fields on this node
  //-----------------------------------------------------------------
  for ( i = 0; i < 2; i++ ){

    uc_l[i]  = (IFloat*)smalloc(MATRIX_SIZE*((local_chi+local_chi_3)/2)*sizeof(IFloat));
     if(uc_l[i] == 0){
       ERR.Pointer(cname,fname, "uc_l[i]"); exit(3);
	}
    for(j=0;j<MATRIX_SIZE*(local_chi+local_chi_3)/2;j++) uc_l[i][j]=0.;
    uc_l_agg[i]  = (gauge_agg *)qalloc(QFAST,((local_chi+local_chi_3)/2)*sizeof(gauge_agg));
   if(uc_l_agg[i] == 0){
      uc_l_agg[i]  = (gauge_agg *)qalloc(QCOMMS,((local_chi+local_chi_3)/2)*sizeof(gauge_agg));
      printf("uc_l_agg[%d] is allocated at DDR (%p)\n",i,uc_l_agg[i]);
    }

     if(uc_l[i] == 0)
       ERR.Pointer(cname,fname, "uc_l[i]");
     if(uc_l_agg[i] == 0)
       ERR.Pointer(cname,fname, "uc_l_agg[i]");
 
    uc_nl[i] = (IFloat*)smalloc( MATRIX_SIZE * ((non_local_chi +non_local_chi_3[3])/2) * sizeof(IFloat) );

    uc_nl_agg[i]  = (gauge_agg*)qalloc(QFAST,((non_local_chi+non_local_chi_3[3])/2)*sizeof(gauge_agg));
    if(uc_nl_agg[i] == 0){
      uc_nl_agg[i]  = (gauge_agg*)qalloc(QCOMMS,((non_local_chi+non_local_chi_3[3])/2)*sizeof(gauge_agg));
      printf("uc_nl_agg[%d] is allocated at DDR (%p)\n",i,uc_nl_agg[i]);
    }

    if(uc_nl[i] == 0){
       ERR.Pointer(cname,fname, "uc_nl[i]");
     }
    if(uc_nl_agg[i] == 0){
       ERR.Pointer(cname,fname, "uc_nl_agg[i]");
     }
  }

FILE *fp;
{
  ParTransAsqtad pt(*lat_pt);
  Float *gauge_p = (Float *)gauge_field_addr;
  Matrix *result[NUM_DIR];
#ifdef USE_NEW
  Matrix *Unit = new Matrix[vol];
  Matrix *P3 = new Matrix[vol];
  Matrix *P3mu = new Matrix[vol];
  Matrix *P5 = new Matrix[vol];
  Matrix *P5mu = new Matrix[vol];
  Matrix *P7 = new Matrix[vol];
  Matrix *P7mu = new Matrix[vol];
  for(j = 0;j<NUM_DIR;j++)
     result[j] = new Matrix[vol];
#else
  Matrix *Unit = (Matrix *)qalloc(QFAST|QCOMMS,sizeof(Matrix)*vol);
  Matrix *P3 = (Matrix *)qalloc(QFAST|QCOMMS,sizeof(Matrix)*vol);
  Matrix *P3mu = (Matrix *)qalloc(QFAST|QCOMMS,sizeof(Matrix)*vol);
  Matrix *P5 = (Matrix *)qalloc(QFAST|QCOMMS,sizeof(Matrix)*vol);
  Matrix *P5mu = (Matrix *)qalloc(QFAST|QCOMMS,sizeof(Matrix)*vol);
  Matrix *P7 = (Matrix *)qalloc(QFAST|QCOMMS,sizeof(Matrix)*vol);
  Matrix *P7mu = (Matrix *)qalloc(QFAST|QCOMMS,sizeof(Matrix)*vol);
  for(j = 0;j<NUM_DIR;j++)
     result[j] = (Matrix *)qalloc(QFAST|QCOMMS,sizeof(Matrix)*vol);
#endif
  Matrix *P6 = P5;
  Matrix *P6mu = P5mu;
  Matrix *Pmumu = P7;
  Matrix *Pmumumu = P7mu;
  for(j = 0;j<vol;j++) Unit[j].UnitMatrix();
     result[j] = new Matrix[vol];
  for(j = 0;j<NUM_DIR;j++)
     for(int k = 0;k<vol;k++) result[j][k].ZeroMatrix();
  IFloat *temp;
  int dirs[] = {6,0,2,4,7,1,3,5}; //mapping between ParTrans and DiracOpAsqtad
  Matrix *min[NUM_DIR],*mout[NUM_DIR];
  for(int mu = 0;mu<NUM_DIR;mu++){
    pt.run(1,&(result[mu]),&Unit,&dirs[mu]);
    for(int nu = 0;nu<NUM_DIR;nu++)
    if(NotParallel(mu,nu)){
      pt.run(1,&P3,&Unit,&dirs[nu]);
      pt.run(1,&P3mu,&P3,&dirs[mu]);
      for(int rho = 0;rho<NUM_DIR;rho++)
      if(NotParallel(mu,rho) && NotParallel(nu,rho)){
        pt.run(1,&P5,&P3,&dirs[rho]);
        pt.run(1,&P5mu,&P5,&dirs[mu]);
        for(int sigma = 0;sigma<NUM_DIR;sigma++)
        if(NotParallel(mu,sigma) && NotParallel(nu,sigma)&&NotParallel(rho,sigma)){
          pt.run(1,&P7,&P5,&dirs[sigma]);
          pt.run(1,&P7mu,&P7,&dirs[mu]);
          int sig_n = (sigma+4)%8;
          pt.run(1,&P7,&P7mu,&dirs[sig_n]);
          fTimesV1PlusV2((IFloat*)P5mu,c7/c5,(IFloat*)P7,(IFloat*)P5mu,vol*18);
        }
        int rho_n = (rho+4)%8;
        pt.run(1,&P5,&P5mu,&dirs[rho_n]);
        fTimesV1PlusV2((IFloat*)P3mu,c5/c3,(IFloat*)P5,(IFloat*)P3mu,vol*18);
      }
      pt.run(1,&P6,&P3,&dirs[nu]);
      pt.run(1,&P6mu,&P6,&dirs[mu]);
      int nu_n = (nu+4)%8;
      pt.run(1,&P6,&P6mu,&dirs[nu_n]);
      fTimesV1PlusV2((IFloat*)P3mu,c6/c3,(IFloat*)P6,(IFloat*)P3mu,vol*18);

      pt.run(1,&P3,&P3mu,&dirs[nu_n]);
      fTimesV1PlusV2((IFloat*)result[mu],c3/c1,(IFloat*)P3,(IFloat*)result[mu],vol*18);
    }
  }

  for ( i = 0; i < 2; i++){
    local_count[i] = 0;
    non_local_count[i] = 0;
    local_count_3[i]=0;
   for (j=0; j<3; j++) non_local_count_3[j][i] = 0;
  }
  //-----------------------------------------------------------------
  //  Loop over all sites.  First rearrange gauge field for this
  //  site and then set up pointers to vector field
  //-----------------------------------------------------------------

  Matrix *tmp,*tmp2;

  for ( n = 0; n < NUM_DIR/2; n++ ) {
    for (x[3] = 0; x[3] < size[3]; x[3]++)
    for (x[2] = 0; x[2] < size[2]; x[2]++)
    for (x[1] = 0; x[1] < size[1]; x[1]++)
    for (x[0] = 0; x[0] < size[0]; x[0]++){
   	 for (i = 0; i < 4 ; i++) coord[i] = x[i];
	    odd = ( coord[0] + coord[1] + coord[2] + coord[3] ) % 2;
	    sg = coord[1] + size[1] * ( coord[2] + size[2] * ( coord[3] +  size[3] * coord[0] ));
      if ( CoordNN( n ) ) {		// chi(x+mu) off-node
        tmp = (Matrix *)(uc_nl[odd] + MATRIX_SIZE * non_local_count[odd]);
        *tmp = result[n][LexGauge(coord)]; 
	*tmp *= (Float)c1;
        non_local_count[odd]++;
      }
      else {
        tmp = (Matrix *)(uc_l[odd] + MATRIX_SIZE * local_count[odd]);
        *tmp = result[n][LexGauge(coord)]; 
	*tmp *= (Float)c1;
        local_count[odd]++;
      }
    }
  }
  for ( n = 0; n < NUM_DIR/2; n++ ) {
    for (x[3] = 0; x[3] < size[3]; x[3]++)
    for (x[2] = 0; x[2] < size[2]; x[2]++)
    for (x[1] = 0; x[1] < size[1]; x[1]++)
    for (x[0] = 0; x[0] < size[0]; x[0]++){
   	 for (i = 0; i < 4 ; i++) coord[i] = x[i];
	    odd = ( coord[0] + coord[1] + coord[2] + coord[3] ) % 2;
	    sg = coord[1] + size[1] * ( coord[2] + size[2] * ( coord[3] +  size[3] * coord[0] ));
      if ( CoordNN( n+4 ) ) {		// chi(x+mu) off-node
        tmp = (Matrix *)(uc_nl[odd] + MATRIX_SIZE * non_local_count[odd]);
//printf("tmp=%p result[%d][%d]=%p\n",tmp,n+4,LexGauge(coord),&(result[n+4][LexGauge(coord)]));fflush(stdout);
        *tmp = result[n+4][LexGauge(coord)]; 
	*tmp *= (Float)-c1;
        non_local_count[odd]++;
      }
      else {
        tmp = (Matrix *)(uc_l[odd] + MATRIX_SIZE * local_count[odd]);
//printf("tmp=%p result[%d][%d]=%p\n",tmp,n+4,LexGauge(coord),&(result[n+4][LexGauge(coord)]));fflush(stdout);
        *tmp = result[n+4][LexGauge(coord)]; 
	*tmp *= (Float)-c1;
        local_count[odd]++;
      }
    }
  }
  for(int mu = 0;mu<NUM_DIR;mu++){
    pt.run(1,&(result[mu]),&Unit,&dirs[mu]);
    pt.run(1,&Pmumu,&result[mu],&dirs[mu]);
    pt.run(1,&(result[mu]),&Pmumu,&dirs[mu]);
  }

  for ( n = 0; n < NUM_DIR/2; n++ ) {
    for (x[3] = 0; x[3] < size[3]; x[3]++)
    for (x[2] = 0; x[2] < size[2]; x[2]++)
    for (x[1] = 0; x[1] < size[1]; x[1]++)
    for (x[0] = 0; x[0] < size[0]; x[0]++){
   	 for (i = 0; i < 4 ; i++) coord[i] = x[i];
	    odd = ( coord[0] + coord[1] + coord[2] + coord[3] ) % 2;
      if ( CoordkNN( n,3 ) ){ 	// chi(x+mu) off-node
        for(int j=0;j<3;j++)
          if(coord_knn[n%4]==j){
            tmp = (Matrix *)(uc_nl[odd] + MATRIX_SIZE * (non_local_count_3[j][odd]+(non_local_chi_3[j]+non_local_chi)/2));
            *tmp = result[n][LexGauge(coord)]; 
	    int index = non_local_count_3[j][odd]+(j+1)*non_local_chi/2;
//        printf("x = %d %d %d %d n= %d uc_nl[%d][%d]=%e ",coord[0],coord[1],coord[2],coord[3],n,odd,index,*((Float *)tmp));
            *tmp *= (Float)c2;
//	printf("tmp[0] =%e(%p)\n",*((Float *)tmp),tmp);
            non_local_count_3[j][odd]++;
          }
      } else {
        tmp = (Matrix *)(uc_l[odd] + MATRIX_SIZE * (local_chi/2+local_count_3[odd]));
        *tmp = result[n][LexGauge(coord)]; 
	    int index = local_count_3[odd]+local_chi/2;
//        printf("x = %d %d %d %d n= %d uc_l[%d][%d]=%e ",coord[0],coord[1],coord[2],coord[3],n,odd,index,*((Float *)tmp));
	*tmp *= (Float)c2;
//	printf("tmp[0] =%e(%p)\n",*((Float *)tmp),tmp);
        local_count_3[odd]++;
      }
    }
  }
  for ( n = 0; n < NUM_DIR/2; n++ ) {
    for (x[3] = 0; x[3] < size[3]; x[3]++)
    for (x[2] = 0; x[2] < size[2]; x[2]++)
    for (x[1] = 0; x[1] < size[1]; x[1]++)
    for (x[0] = 0; x[0] < size[0]; x[0]++){
   	 for (i = 0; i < 4 ; i++) coord[i] = x[i];
	    odd = ( coord[0] + coord[1] + coord[2] + coord[3] ) % 2;
      if ( CoordkNN( n+4,3 ) ){ 	// chi(x+mu) off-node
        for(int j=0;j<3;j++)
          if(coord_knn[n%4]==(size[n%4]-1-j)){
            tmp = (Matrix *)(uc_nl[odd] + MATRIX_SIZE * (non_local_count_3[j][odd]+(non_local_chi_3[j]+non_local_chi)/2));
	    int index = non_local_count_3[j][odd]+(j+1)*non_local_chi/2;
            *tmp = result[n+4][LexGauge(coord)]; 
//        printf("x = %d %d %d %d n= %d uc_nl[%d][%d]=%e ",coord[0],coord[1],coord[2],coord[3],n+4,odd,index,*((Float *)tmp));
            *tmp *= (Float)-c2;
//	printf("tmp[0] =%e(%p)\n",*((Float *)tmp),tmp);
            non_local_count_3[j][odd]++;
          }
      } else {
        tmp = (Matrix *)(uc_l[odd] + MATRIX_SIZE * (local_chi/2+local_count_3[odd]));
        *tmp = result[n+4][LexGauge(coord)]; 
	    int index = local_count_3[odd]+local_chi/2;
//        printf("x = %d %d %d %d n= %d uc_l[%d][%d]=%e ",coord[0],coord[1],coord[2],coord[3],n+4,odd,index,*((Float *)tmp));
	*tmp *= (Float)-c2;
//	printf("tmp[0] =%e(%p)\n",*((Float *)tmp),tmp);
        local_count_3[odd]++;
      }
    }
  }

#ifdef USE_NEW
  delete[] Unit;
  delete[] P3;
  delete[] P3mu;
  delete[] P5;
  delete[] P5mu;
  delete[] P7;
  delete[] P7mu;
  for(int j = 0;j<NUM_DIR;j++)
  delete[] result[j];
#else
  qfree( Unit);
  qfree( P3);
  qfree( P3mu);
  qfree( P5);
  qfree( P5mu);
  qfree( P7);
  qfree( P7mu);
  for(int j = 0;j<NUM_DIR;j++)
  qfree( result[j]);
#endif
}

  int fd;
  char buf[200];
#if 0
  fp=fopen(uc_l_filename,"w");
  for(j=0;j<2;j++){
    fprintf(fp,"IFloat uc_l%d[] LOCATE(\"edramnormal\") = {\n",j);
    for(i=0;i< MATRIX_SIZE * ((local_chi+ local_chi_3)/2);i++){
      fprintf(fp,"%0.4e, ",*(uc_l[j]+i));
      if( i%6 == 5){
        fprintf(fp,"\n");
      }
    }
    fprintf(fp,"\n};\n"); 
  }
  fclose(fp);
#endif

  gauge_agg *temp = new gauge_agg[12*vol];
  int num_ind[vol*6];
  int src;
//  printf("chi_l=0x%x\nchi_l[0][0]=%d\n",chi_l,chi_l[0][0]);
  for(j=0;j<2;j++){
    for(i=0;i<vol;i++) num_ind[i]=0;
    for(i=0;i< ((local_chi+ local_chi_3)/2);i++){
      src = (int)chi_l[j][2*i];
//      printf("chi_l[%d][%d]=%d %d\n",j,i*2,chi_l[j][2*i],src);
      if (src%(VECT_LEN*sizeof(IFloat))!=0){
        ERR.General(cname,fname,"src = %d\n",src);
        exit(1);
      }
      src = src/(VECT_LEN*sizeof(IFloat));
//      printf("src[%d]=%d\n",i,src);
      if(src > vol/2) {
        ERR.General(cname,fname,"src[%d](%d) > vol/2\n",i,src);
        exit(1);
      }
      temp[src*NUM_DIR*2+num_ind[src]].src = (int)chi_l[j][2*i];
      temp[src*NUM_DIR*2+num_ind[src]].dest= (int)chi_l[j][2*i+1];
      for(k=0;k<MATRIX_SIZE;k++){
        temp[src*NUM_DIR*2+num_ind[src]].mat[k] = uc_l[j][i*MATRIX_SIZE+k];
      }
      num_ind[src]++;
    }
    int index = 0;
    int div = 4;
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
      if (n==0) odd_num +=num_ind[i]%div;
    }
  }

  gauge_agg *agg_p;
#if 0
  fp=fopen(uc_l_agg_filename,"w");
  for(j=0;j<2;j++){
    fprintf(fp,"struct gauge_agg uc_l_agg%d[] LOCATE(\"edramtransient\") = {\n",j);
    for(i=0;i< ((local_chi+ local_chi_3)/2);i++){
      agg_p = &(uc_l_agg[j][i]);
      fprintf(fp,"{%d,%d,{\n",agg_p->src,agg_p->dest);
      for(k=0;k<18;k++){
        fprintf(fp,"%0.8e",agg_p->mat[k]);
	if(k!=17){
        fprintf(fp,", ");
	}
        if( k%6 == 5){
          fprintf(fp,"\n");
        }
      }
      fprintf(fp,"}},\n");
    }
    fprintf(fp,"\n};\n"); 
  }
  fclose(fp);
#endif



#if 0
  fp = fopen(uc_nl_filename,"w");
  for(j=0;j<2;j++){
    fprintf(fp,"IFloat uc_nl%d[] LOCATE(\"edramnormal\") = {\n",j); 
    for(i=0;i< MATRIX_SIZE * ((non_local_chi+non_local_chi_3[3] )/2);i++){
      fprintf(fp,"%0.4e, ",*(uc_nl[j]+i));
      if( i%6 == 5){
        fprintf(fp,"\n");
      }
    }
    fprintf(fp,"\n};\n"); 
  }
  fclose(fp);
#endif

  for(j=0;j<2;j++){
    for(i=0;i<non_local_chi*3/2;i++) num_ind[i]=0;
    for(i=0;i< ((non_local_chi+ non_local_chi_3[3])/2);i++){
      src = (int)chi_nl[j][2*i];
      if (src%(VECT_LEN*sizeof(IFloat))!=0){
        ERR.General(cname,fname,"src = %d\n",src);
        exit(1);
      }
      src = src/(VECT_LEN*sizeof(IFloat));
      if(src > non_local_chi*3/2) {
        ERR.General(cname,fname,"src(%d) > non_local_chi*3\n",src);
        exit(1);
      }
      temp[src*2+num_ind[src]].src = (int)chi_nl[j][2*i];
      temp[src*2+num_ind[src]].dest= (int)chi_nl[j][2*i+1];
      for(k=0;k<18;k++)
        temp[src*2+num_ind[src]].mat[k] = uc_nl[j][i*18+k];
      num_ind[src]++;
      if(num_ind[src]>2){
        ERR.General(cname,fname,"num_ind[%d](%d) > 2 \n",src, num_ind[src]);
        exit(1);
      }
    }

    int index = 0;
    int div = 2;
    for( n=0; n*div<2; n++)
    for(i=0;i<non_local_chi*3/2;i++){
      for(m=0;m<div;m++){
        if(num_ind[i] > n*div+m) {
           uc_nl_agg[j][index] = temp[i*2+n*div+m];
           index++;
        }
      }
    }
  }
  delete[] temp;

#if 0
  fp = fopen(uc_nl_agg_filename,"w");
  for(j=0;j<2;j++){
    fprintf(fp,"struct gauge_agg uc_nl_agg%d[] LOCATE(\"edramtransient\") = {\n",j);
    for(i=0;i< ((non_local_chi+ non_local_chi_3[3])/2);i++){
      agg_p = &(uc_nl_agg[j][i]);
      fprintf(fp,"{%d,%d,{\n",agg_p->src,agg_p->dest);
      for(k=0;k<18;k++){
        fprintf(fp,"%0.8e",agg_p->mat[k]);
	if(k!=17){
        fprintf(fp,", ");
	}
        if( k%6 == 5){
          fprintf(fp,"\n");
        }
      }
      fprintf(fp,"}},\n");
    }
    fprintf(fp,"\n};\n"); 
  }
  fclose(fp);
#endif



}

extern "C"
void asqtad_destroy_dirac_buf_g(void)
{
  int i;
  for ( i = 0; i < 2; i++){
  sfree(uc_l[i]);
  sfree(uc_nl[i]);
  qfree(uc_l_agg[i]);
  qfree(uc_nl_agg[i]);
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
#if 0
     while (coord_knn[m] < 0) {
        coord_knn[m] += size[m] ;
         }
       coord_knn[m] %= size[m] ;
#endif
 
  if (nn<4 ){ if (coord_knn[m] >(size[m]-1)) {coord_knn[m] -=size[m];off_node = 1;}}
       else {if (coord_knn[m] <0 ) {coord_knn[m]  +=size[m];off_node = 1;}}
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
  int i,odd;
  for(odd=0;odd<2;odd++)
  for(i=0;i<2*NUM_DIR;i++){
     SCUarg[odd][i].Assert();
     SCUarg_1[odd][i].Assert();
  }
}

//-------------------------------------------------------------------
//  add_flag = 0 :       b = D a
//  add_flag = 1 :       b += D a
//
//  a_odd    = 1 :	 b even;  a odd
//  a_odd    = 0 :	 b odd ;  a even
//
//   D a = \sum_u( U^dag_u(x) a(x+u) - U(x-u) a(x-u) )
//-------------------------------------------------------------------


static int called=0;
static unsigned long address[]={0x0,0x0};
extern "C"
void asqtad_dirac(IFloat* b, IFloat* a, int a_odd, int add_flag)
{
  int i,j,k;
  int odd = 1 - a_odd;
  long c = (long) chi_off_node_total;
  long uc_l = (long)uc_l_agg[odd];
  long uc_nl = (long)uc_nl_agg[odd];
  long uc_nl2 = (long)&(uc_nl_agg[odd][isplit]);
  long cache_p0=uc_l;
  long cache_p1=cache_p0+32;
  long cache_p2=cache_p1+32;
  long cache_p3=cache_p2+32;
  long cache_p4=cache_p3+32;
  long cache_p5=cache_p4+32;
  long cache_p6=cache_p5+32;
  long cache_p7=cache_p6+32;
  int num_flops;
  Float dtime;
  struct timeval start,end;
  

#undef PROFILE
#ifdef PROFILE
  num_flops = 0;
  gettimeofday(&start,NULL);
#endif
  if( (unsigned long)a != address[odd]){
//    printf("a[%d] = %p -> %p\n",odd,address[odd],a);
    address[odd] = (unsigned long)a;
  
    for(i=1;i<4;i++)
    if(size[i] >2 ){
      SCUarg[odd][i].Addr( a + Xoffset[2][i]);
      SCUarg[odd][i+4].Addr( a + Xoffset[2][i+4]);
    }
  
    void *addr[2];
    if(split) {
      for(i=1;i<4;i++){
        addr[0] = (void *)(a+Xoffset[0][i]);
        addr[1] = (void *)(a+Xoffset[1][i]);
        SCUarg_1[odd][i].Addr( addr,1);
        SCUarg_2[odd][i].Addr( addr+1,1);
        addr[0] = (void *)(a+Xoffset[0][i+4]);
        addr[1] = (void *)(a+Xoffset[1][i+4]);
        SCUarg_1[odd][i+4].Addr( addr,1);
        SCUarg_2[odd][i+4].Addr( addr+1,1);
      }
    } else {
      for(i=1;i<4;i++){
        addr[0] = (void *)(a+Xoffset[0][i]);
        addr[1] = (void *)(a+Xoffset[1][i]);
        SCUarg_1[odd][i].Addr( addr,2);
        addr[0] = (void *)(a+Xoffset[0][i+4]);
        addr[1] = (void *)(a+Xoffset[1][i+4]);
        SCUarg_1[odd][i+4].Addr( addr,2);
      }
    }
  } // if address[odd] != a


  //-----------------------------------------------------------------
  //  Transfer chi's on faces.  
  //-----------------------------------------------------------------

  copy_buffer(countM[0][a_odd], (long)a, (long)Tbuffer[0][0], (long)ToffsetM[0][a_odd]);
  copy_buffer(countP[0][a_odd], (long)a, (long)Tbuffer[0][1], (long)ToffsetP[0][a_odd]);
  if(size[0]>2) {
  copy_buffer(countM[1][a_odd], (long)a, (long)Tbuffer[1][0], (long)ToffsetM[1][a_odd]);
  copy_buffer(countP[1][a_odd], (long)a, (long)Tbuffer[1][1], (long)ToffsetP[1][a_odd]);
  }
  if(size[0]>4) {
  copy_buffer(countM[2][a_odd], (long)a, (long)Tbuffer[2][0], (long)ToffsetM[2][a_odd]);
  copy_buffer(countP[2][a_odd], (long)a, (long)Tbuffer[2][1], (long)ToffsetP[2][a_odd]);
  }
#if 0
	if(called<2)
    for (int i = 0; i<(vol*VECT_LEN/2);i++)
	if (fabs(a[i]) >1e-10)
		fprintf(stderr,"a[%d]=%e\n",i,a[i]);

	if(called<2)
    for (int i = 0; i<1;i++)
    for (int j = 0; j<2;j++)
	for (int k = 0; k<(size[1]*size[2]*size[3]*VECT_LEN)/2;k++){
	if (fabs(*(Tbuffer[i][j]+k)) >1e-10)
		fprintf(stderr,"Tbuffer[%d][%d][%d]=%e\n",i,j,k,*(Tbuffer[i][j]+k));
    }
#endif
#ifdef PROFILE
  gettimeofday(&end,NULL);
  print_flops(0,&start,&end);
#endif

  sys_cacheflush(0);

  //make sure spinor field is in main memory before starting transfers

  SCUmulti_1[odd].StartTrans();

#undef PROFILE
#ifdef PROFILE
  num_flops = 33*(local_chi + local_chi_3);
  gettimeofday(&start,NULL);
#endif

  //-----------------------------------------------------------------
  //do first local computations
  //-----------------------------------------------------------------
  dcbt(cache_p0);cache_p0=uc_nl;
  dcbt(cache_p1);cache_p1=cache_p0+32;
  dcbt(cache_p2);cache_p2=cache_p1+32;
  dcbt(cache_p3);cache_p3=cache_p2+32;
  dcbt(cache_p4);cache_p4=cache_p3+32;
  dcbt(cache_p5);cache_p5=cache_p4+32;
  dcbt(cache_p6);cache_p6=cache_p5+32;
#ifdef CPP
  dirac_cmv_jcw_agg_cpp( (local_chi + local_chi_3)/2, (long)0, uc_l, (long)a, (long)tmpfrm);
#else
  if (odd_num)
  dirac_cmv_jcw_agg( (local_chi + local_chi_3)/2, (long)0, (long)uc_l_agg[odd],
	       (long)a, (long)tmpfrm);
  else
  dirac_cmv_agg_4( (local_chi + local_chi_3)/2, (long)0, uc_l, (long)a, (long)tmpfrm);
#endif


  //-----------------------------------------------------------------
  // check to see if transfers are done and start another transfer
  //-----------------------------------------------------------------

#ifdef PROFILE
  gettimeofday(&end,NULL);
  print_flops(num_flops,&start,&end);
#endif

  SCUmulti_1[odd].TransComplete();

if(split) {

  SCUmulti_2[odd].StartTrans();

#undef PROFILE
#ifdef PROFILE
  num_flops = 66*isplit;
  gettimeofday(&start,NULL);
#endif

  dcbt(cache_p0);cache_p0=uc_nl2;
  dcbt(cache_p1);cache_p1=cache_p0+32;
  dcbt(cache_p2);cache_p2=cache_p1+32;
  dcbt(cache_p3);cache_p3=cache_p2+32;
  dcbt(cache_p4);cache_p4=cache_p3+32;
  dcbt(cache_p5);cache_p5=cache_p4+32;
  dcbt(cache_p6);cache_p6=cache_p5+32;
#ifdef CPP
  dirac_cmv_jcw_agg_cpp( isplit, (long)0, (long)&(uc_nl_agg[odd][0]), (long)c, (long)tmpfrm);
#else
  dirac_cmv_jcw_agg( isplit, (long)0, uc_nl, c, (long)tmpfrm);
#endif

#ifdef PROFILE
  gettimeofday(&end,NULL);
  print_flops(num_flops,&start,&end);
#endif

  SCUmulti_2[odd].TransComplete();

  SCUmulti[odd].StartTrans();

#ifdef PROFILE
  num_flops = 33*(non_local_chi);
  gettimeofday(&start,NULL);
#endif

#ifdef CPP
  dirac_cmv_jcw_agg_cpp( non_local_chi/2, (long)0, (long)&(uc_nl_agg[odd][isplit]), (long)c, (long)tmpfrm);
#else
  dirac_cmv_jcw_agg( non_local_chi/2, (long)0, uc_nl2, (long)c, (long)tmpfrm);
#endif //#ifdef CPP


} else {

  SCUmulti[odd].StartTrans();

  //-----------------------------------------------------------------
  //do the computations involving "chi" non-local spinors
  //-----------------------------------------------------------------


#ifdef PROFILE
  num_flops = 33*(non_local_chi + non_local_chi_3[2]);
  gettimeofday(&start,NULL);
#endif

#ifdef CPP
  dirac_cmv_jcw_agg_cpp( (non_local_chi+non_local_chi_3[2])/2, (long)0, (long)&(uc_nl_agg[odd][0]), (long)c, (long)tmpfrm);
#else
  dirac_cmv_jcw_agg( (non_local_chi+non_local_chi_3[2])/2, (long)0, (long)&(uc_nl_agg[odd][0]), (long)c, (long)tmpfrm);
#endif

}


#ifdef PROFILE
  gettimeofday(&end,NULL);
  print_flops(num_flops,&start,&end);
#endif
  SCUmulti[odd].TransComplete();

 

  //-----------------------------------------------------------------
  //do the computations involving chi3 non-local spinors
  //----------------------------------------------------------------


#undef PROFILE
#ifdef PROFILE
//  num_flops = 1146*vol/2;
  num_flops = 33*non_local_chi;
  dtime = -dclock();
#endif

#ifdef CPP
  dirac_cmv_jcw_agg_cpp( non_local_chi/2, (long)0, (long)&(uc_nl_agg[odd][(non_local_chi+non_local_chi_3[2])/2]) , (long)c, (long)tmpfrm);
#else
  dirac_cmv_jcw_agg( non_local_chi/2, (long)0, (long)&(uc_nl_agg[odd][(non_local_chi+non_local_chi_3[2])/2]) , (long)c, (long)tmpfrm);
#endif

#ifdef PROFILE
  dtime +=dclock();
  printf("dirac_cmv_jcw_agg::%ld flops/%e seconds = %e MFlops\n",num_flops,dtime,(Float)num_flops/(dtime*1e6));
#endif

#if 0
	if(called<2)
    for (int i = 0; i<3*non_local_chi* VECT_LEN / 2;i++ ){
	if (fabs(chi_off_node_total[i]) >1e-10)
		fprintf(stderr,"chi_off_node_total[%d]=%e\n",i,chi_off_node_total[i]);
    }
	called++;
#endif


  //printf ("the computations involving chi3 non-local spinors done \n");
  // check to see if transfers are done
  //-----------------------------------------------------------------
  //do the sum of 16 temporary vectors at each lattice site
  //              ^^^ change must be made in  dirac_sum**
  //-----------------------------------------------------------------


#undef PROFILE
#ifdef PROFILE
  num_flops = 45*vol;
  dtime = -dclock();
#endif
  if ( add_flag == 0){
    dirac_sum2_64( vol/2, (long)0, (long)tmpfrm, (long)b);
  }
  else{
    dirac_sum_acc_cpp( vol/2, (long)0, (long)tmpfrm, (long)b);
  }
#ifdef PROFILE
  dtime +=dclock();
  printf("dirac_cmv_jcw_agg::%ld flops/%e seconds = %e MFlops\n",num_flops,dtime,(Float)num_flops/(dtime*1e6));
#endif

}

CPS_END_NAMESPACE
