//  $Source: /space/cvs/cps/cps++/src/util/dirac_op/d_op_asqtad/qcdoc/asqtad_dirac.C,v $
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2006/06/11 05:35:06 $
//  $Header: /space/cvs/cps/cps++/src/util/dirac_op/d_op_asqtad/qcdoc/asqtad_dirac.C,v 1.30 2006/06/11 05:35:06 chulwoo Exp $
//  $Id: asqtad_dirac.C,v 1.30 2006/06/11 05:35:06 chulwoo Exp $
//  $Name: v5_0_16_hantao_io_test_v7 $
//  $Locker:  $
//  $RCSfile: asqtad_dirac.C,v $
//  $Revision: 1.30 $
//  $Source: /space/cvs/cps/cps++/src/util/dirac_op/d_op_asqtad/qcdoc/asqtad_dirac.C,v $
//  $State: Exp $
//
//    12/02 HueyWen Lin, Chulwoo Jung
//
//   Asqtad Dirac operator for QCDOC. Communications and computations
//   are overlapped.
//   Uses functions implemented by CJ for QCDOC
//
//
//-------------------------------------------------------------------

#include <string.h>
#include <sys/time.h>
#include "asq_data_types.h"
#include "asqtad_int.h"
#include <qcdocos/scu_mmap.h>
#include <ppc_lib.h>
#include <qcdoc.h>
#include <math.h>

Float asq_print_flops(char *cname, char *fname, unsigned long long nflops, struct
timeval *start, struct timeval *end);

void matrix::Dagger(const Float* a)
{
    u[0]  = a[0];   u[1]  = -a[1];
    u[6]  = a[2];   u[7]  = -a[3];
    u[12] = a[4];   u[13] = -a[5];
    u[2]  = a[6];   u[3]  = -a[7];
    u[8]  = a[8];   u[9]  = -a[9];
    u[14] = a[10];  u[15] = -a[11];
    u[4]  = a[12];  u[5]  = -a[13];
    u[10] = a[14];  u[11] = -a[15];
    u[16] = a[16];  u[17] = -a[17];
}


#if 0
const char *uc_l_filename = "uc_l.h";
const char *uc_nl_filename = "uc_nl.h";
const char *Toffset_filename = "Toffset.h";
const char *uc_l_agg_filename = "uc_l_agg.h";
const char *uc_nl_agg_filename = "uc_nl_agg.h";
const char *chi_l_filename = "chi_l.h";
const char *chi_nl_filename = "chi_nl.h";
#endif

/*****************************************************************
 SIMUL switched on/off the hack CJ put in to help speed up the
simulation. if SIMUL is undefined, program writes the temporary arraies
for dirac operator. If SIMUL is defined, it will include (pre-calcuated)
temporary arraies and skip the generation of arraies.
******************************************************************/

/****************************************************************
CPP is a switch for using C++ routine for dirac_cmv.
*****************************************************************/

void asqd_sum_acc_cpp(int s, long chi, long tmpfrm, long b);

#define NEW_SPLIT

#undef CPP
#ifdef CPP
void dirac_sum2_64_cpp( int sites, long chi, long tmpfrm,long b);
void dirac_cmv_jcw_agg_cpp( int sites, long chi, long u,long a, long tmpfrm);
#define asq_cmv(A,B,C,D) dirac_cmv_jcw_agg_cpp(A,0,B,C,D)
#define asq_dsum(A,B,C,D) dirac_sum2_64_cpp(A,B,C,D)
#else
#ifdef ASQD_SINGLE
extern "C" void asq_cmv_s( int sites, long u,long a, unsigned long tmpfrm);
extern "C" void asq_cmv_4_s( int sites, long chi,long u,long a, unsigned long tmpfrm);
#define asq_cmv_4(A,B,C,D) asq_cmv_4_s(A,0,B,C,D)
#define asq_cmv(A,B,C,D) asq_cmv_s(A,B,C,D)
extern "C" void asq_dsum_s( int sites, unsigned long tmpfrm,long b);
#define asq_dsum(A,B,C) asq_dsum_s(A,B,C)
extern "C" void copy_buffer_s(int n, long src, long dest, long ptable);
#define copy_buffer(A,B,C,D) copy_buffer_s(A,B,C,D)
#else
extern "C" void  dirac_cmv_jcw_agg( int sites, long chi, long u,long a, long tmpfrm);
extern "C" void  dirac_cmv_agg_4( int sites, long chi, long u,long a, unsigned long tmpfrm);
extern "C" void asq_cmv( int sites, long u,long a, unsigned long tmpfrm);
#define asq_cmv_4(A,B,C,D) dirac_cmv_agg_4(A,0,B,C,D)
extern "C" void asq_dsum( int sites, unsigned long tmpfrm,long b);
extern "C" void copy_buffer(int n, long src, long dest, long ptable);
#endif
#endif

#if 0
extern "C" void copy_buffer_cpp(int n, long src, long dest, long ptable);
#define copy_buffer(A,B,C,D) copy_buffer_cpp(A,B,C,D)
#endif

enum{VECT_LEN=6, VECT_LEN2=8, MATRIX_SIZE=18, SITE_LEN=72, NUM_DIR=8, N=4};

static inline int PAD (int a){
  if (a%2) return a+1;
  else return a;
}

static inline int PAD4 (int a){
  if (a%4) return a+4-(a%4);
  else return a;
}

//---------------------------------------------------------------------
//  Arrays for send direction, receive direction, lattice size per
//  node, coordinate (for site where cluster of 8 matrices
//  is being assembled) and the coordinates of a nearest neighbor site.
//
//  0 == T, 1 == X, 2 == Y, 3 == Z
//--------------------------------------------------------------------

SCUDir AsqD::scudir[] =
{
  SCU_TP, SCU_XP, SCU_YP, SCU_ZP, SCU_TM, SCU_XM, SCU_YM, SCU_ZM
};


static int Rotate (int mu, int i){
	int mu_p = (mu+i)%4;
	if( mu >= 4)
		mu_p += 4;
	return mu_p;
}
static int NotParallel( int mu, int nu){
	int dif = mu-nu;
	if (dif==0 ||dif==-4||dif==4)return 0;
	else return 1;
}

static Float * gauge_field_addr;

int k;

//-------------------------------------------------------------------
//  Called by the lattice constructor
//  Fermion initializations: pointer tables 
//-------------------------------------------------------------------
void AsqD::init(AsqDArg *arg)
{

//  printf("AsqD::init()\n");
  gauge_field_addr = ( Float * ) arg->gauge_u;
  int i,j,k,m,n; 
  int blklen[NUM_DIR/2];
  int numblk[NUM_DIR/2];
  int stride[NUM_DIR/2];
  int local_count[2];
  int non_local_count[2];
  int local_count_3[2];
  int non_local_count_3[3][2];
  int x[NUM_DIR/2];
  int scu_irs[2][3];
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

  for(int i = 0;i<4;i++) size[i] = arg->size[i];
  for(int i = 0;i<4;i++) coor[i] = arg->coor[i];
  non_local_dirs=0;
  for(int i = 0;i<4;i++){
     NP[i] = arg->NP[i];
     if (NP[i]>1) non_local_dirs++;
//     printf("size[%d] coor[%d] NP[%d]= %d %d %d \n",i,i,i,size[i],coor[i],NP[i]);
  }
  c1 = arg->c1;
  c2 = arg->c2;
  c3 = -arg->c3;
  c5 = arg->c5;
  c7 = -arg->c7;
  c6 = arg->c6;
  for(int i = 0;i<4;i++){
    fat[i] = (matrix *)arg->Fat[i];
    naik[i] = (matrix *)arg->Naik[i];
    naik_m[i] = (matrix *)arg->NaikM[i];
  }

  node_odd=0;
  for(int i = 0;i<4;i++)
    node_odd += size[i]*coor[i];
  node_odd = node_odd%2;
//  fprintf(stdout,"node_odd=%d\n",node_odd);

  vol = size[0] * size[1] * size[2] * size[3];
  f_size_cb = vol*3;
//  split =  (  vol>64 ? 0 : 1 ); 
  split = 1;
  if (vol >1000 || (size[0]%2) || (size[1]%2) || (size[2]%2) || (size[3]%2) )
  split = 0;
//  printf("vol=%d split=%d\n",vol,split);


  int area[2][4];
  for(i=0;i<4;i++){
     int total= vol/size[i];
     area[0][i] = PAD(total)/2;
     area[1][i] = total-area[0][i];
//     printf("area[0][%d]=%d area[1][%d]=%d\n",i,area[0][i],i,area[1][i]);
  }

//  non_local_chi = 2*(size[0]*size[1]*size[2] + size[1]*size[2]*size[3]+
//    size[2]*size[3]*size[0] + size[3]*size[0]*size[1]);
  for(int i =0;i<2;i++)
  non_local_chi[i] = 0;
  for(int i =0;i<4;i++)
  if(NP[i]>1){
     if (size[i]%2){
       non_local_chi[0] += 2*area[0][i];
       non_local_chi[1] += 2*area[1][i];
     } else {
       non_local_chi[0] += area[0][i]+area[1][i];
       non_local_chi[1] += area[0][i]+area[1][i];
     }
  }
#if 0
  for(int i =0;i<2;i++)
    printf("non_local_chi[%d]=%d  ",i,non_local_chi[i]);
  printf("\n");
#endif
 
  for(int i =0;i<2;i++){
  non_local_chi_3[i][0] = 0;
  non_local_chi_3[i][1] = 0;
  }

  for(i=0;i<4;i++)
  if(NP[i]>1 && (size[i]>2) ){
     if (size[i]%2){
       non_local_chi_3[0][1] += 2*area[0][i];
       non_local_chi_3[1][1] += 2*area[1][i];
     } else {
       non_local_chi_3[0][1] += area[0][i]+area[1][i];
       non_local_chi_3[1][1] += area[0][i]+area[1][i];
     }
  }

  for(int i =0;i<2;i++){
  non_local_chi_3[i][2]  = non_local_chi_3[i][1]+non_local_chi[1-i];
  non_local_chi_3[i][3]  = non_local_chi_3[i][2]+non_local_chi[i];
  }
#if 0
  for(int j =0;j<4;j++){
  for(int i =0;i<2;i++)
    printf("non_local_chi_3[%d][%d]=%d  ",i,j,non_local_chi_3[i][j]);
    printf("\n");
  }
#endif

  int chk_vol[2];
  
  chk_vol[0] = PAD(vol)/2;
  chk_vol[1] = vol-chk_vol[0];

  for(int i =0;i<2;i++){
  local_chi[i] = NUM_DIR*chk_vol[i] - non_local_chi[i];
  local_chi_3[i] = NUM_DIR*chk_vol[i] - (non_local_chi_3[i][3]);
  }
  
  isplit  = (non_local_chi_3[0][1]+non_local_chi[0]);
  if(isplit %2){
    printf("isplit=%d\n",isplit);exit(-4);
  }

#if 0
  for(int i =0;i<2;i++)
    printf("local_chi[%d]=%d  ",i,local_chi[i]);
  printf("\n");
  for(int i =0;i<2;i++)
    printf("local_chi_3[%d]=%d  ",i,local_chi_3[i]);
  printf("\n");
#endif

  for(int i =0;i<2;i++){
    comp_l[i] = local_chi[i]+local_chi_3[i];
    comp_nl[i] = non_local_chi[i]+non_local_chi_3[i][2];
    comp_nl_2[i] = non_local_chi[i];
  }


  //-------------------------------------------------------------
  // flush_cache_spinor() function will flush 192 bytes * nflush 
  //-------------------------------------------------------------
  nflush = vol/8;

  //-----------------------------------------------------------------
  //  Allocate 8 receive buffers for off-node vectors
  //-----------------------------------------------------------------
  if(vol>4096) chi_off_node_total=NULL;
    else
    chi_off_node_total = ( Float * ) qalloc(QFAST|QCOMMS,
      3*PAD(non_local_chi[0])* VECT_LEN * sizeof( Float )  );
  if(chi_off_node_total == NULL){ 
    chi_off_node_total = ( Float * ) qalloc(QCOMMS, 
      3*PAD(non_local_chi[0])*VECT_LEN * sizeof( Float )  );
    printf("chi_off_node_total is allocated at DDR (%p)\n",chi_off_node_total);
  }
    if(chi_off_node_total == 0)
      PointerErr(cname,fname, "chi_off_node_total");
//  printf("chi_off_node_total=%p\n",chi_off_node_total);

 for (int k= 0; k < 2; k++ )
 for ( j= 0; j < 3; j++ ){
    chi_off_node[k][j][0] = &(chi_off_node_total[ PAD(non_local_chi[0])*j* VECT_LEN ] ); 
    chi_off_node_p[k][j][0] = (Float *)(sizeof (Float)*PAD(non_local_chi[0])*j* VECT_LEN );
  for ( i = 1; i < NUM_DIR; i++ ){
    int offset = 0;
    int dir = (i-1)%4;
    if (NP[dir]>1){
      offset = area[k][dir];
      if((i-1)==0) offset = area[(k+j+1)%2][0]; //SCU_TP
      if((i-1)==4) offset = area[(k+j)%2][0]; //SCU_TM
    }
    offset *= VECT_LEN;
    chi_off_node[k][j][i] = chi_off_node[k][j][i-1]+offset;
    chi_off_node_p[k][j][i] = chi_off_node_p[k][j][i-1]+offset;
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
    chi_l[i] = ( Float ** ) Alloc((2*(local_chi[i]+local_chi_3[i]))*sizeof(Float *));
   if(chi_l[i] == NULL)
	PointerErr(cname,fname, "chi_l[i]");

    chi_nl[i] = (Float ** ) Alloc((2* (non_local_chi[i] + non_local_chi_3[i][3]))*sizeof(Float *));

      if(chi_nl[i] == NULL)
	PointerErr(cname,fname, "chi_nl[i]");
  }
  

  for ( i = 0; i < 2; i++){
    local_count[i] = 0;
    non_local_count[i] = 0;
    local_count_3[i] = 0;
    for( n = 0; n < 3; n++ )  non_local_count_3[n][i] = 0;

  }

  for ( int k = 0; k < 3; k++ ) {
  for ( i = 0; i < 2; i++ ) {
//  if(vol>1024) Tbuffer[k][i]=NULL;
//    else
    Tbuffer[k][i] = (Float *) qalloc (QFAST|QNONCACHE, PAD4(area[0][0])* VECT_LEN * sizeof( Float ) );

  if( Tbuffer[k][i] == NULL){
    Tbuffer[k][i] = (Float *) qalloc (QCOMMS, PAD4(area[0][0])* VECT_LEN * sizeof( Float ) );
    printf("Tbuffer[%d][%d] is allocated at DDR (%p)\n",k,i,Tbuffer[k][i]);
  }

   if(Tbuffer[k][i] == NULL)
	PointerErr(cname,fname, "Tbuffer[i][j]");
    ToffsetP[k][i] = ( int * ) qalloc (0,  PAD4(area[0][0])* sizeof( int ) );
   if(ToffsetP[k][i] == NULL)
	PointerErr(cname,fname, "TOffsetP[i][j]");
    ToffsetM[k][i] = ( int * ) qalloc (0,  PAD4(area[0][0])* sizeof( int ) );
   if(ToffsetM[k][i] == NULL)
	PointerErr(cname,fname, "TOffsetM[i][j]");
    countP[k][i] = 0;
    countM[k][i] = 0;
  }
 }
// printf("Tbuffer done\n");


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

	    if ( CoordNN( n )  && (NP[n%4]>1) ) {	//  off-node
	      //----------------------------------------------------------
	      // Assembly written for double precision only, multiplication
	      // by sizeof(double) done to avoid a bitshift inside the
	      // high performance code
	      //----------------------------------------------------------
	      //pointer to source field (offset in the receive buffer)

	      *( chi_nl[ odd ]  +  2 * non_local_count[ odd ] )
	     		= chi_off_node_p[odd][0][n] + VECT_LEN * ( LexSurface( coord_nn, n%4 ) / 2 );
	      // pointer to temporary field where U*chi is stored
	      *( chi_nl[ odd ] + 2*non_local_count[ odd ] +1 ) =
				( Float * ) ( VECT_LEN2 * (LexVector( coord ) / 2 * 2*NUM_DIR+n )  * sizeof(Float));
	      non_local_count[odd]++;
	    }
    }
    }
    }
    }
  }

  for(int i = 0;i<2;i++)
  if (non_local_count[i] != non_local_chi[i]){
    printf(" non_local_count[%d]=%d\n",i,non_local_count[i]);
    exit(-4);
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

	    if ( !CoordNN( n )  || NP[n%4] <2 ) {		//  off-node
	      //pointer to source field
	      *( chi_l[ odd ]  +  2 * local_count[ odd ] )
		= ( Float * ) ( VECT_LEN * ( LexVector( coord_nn ) / 2 ) * sizeof(Float));
	      // pointer to temporary field where U*chi is stored
	      *( chi_l[ odd ]  +  2 * local_count[ odd ] + 1) =
					( Float * ) ( VECT_LEN2 * (LexVector( coord ) / 2 *2*NUM_DIR+n) * sizeof(Float));
	      local_count[odd]++;
	    }  //else-on_node case
    }
    }
    }
    }
  }

  for(int i = 0;i<2;i++)
  if (local_count[i] != local_chi[i]){
    printf(" local_count[%d]=%d\n",i,local_count[i]);
    exit(-4);
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

  if ( NP[n%4] >1 && CoordkNN( n, 3 ) ) {		//  chi3_off-node
    for (j=0; j<3; j++){
       int offset = non_local_chi[odd]+non_local_chi_3[odd][j];
       if((coord_knn[n%4]==j&&n<4)||(coord_knn[n%4]==(size[n%4]-1-j)&& n>3)) {
	 *( chi_nl[ odd ] + 2 * non_local_count_3[j][odd]+2*offset )
            =chi_off_node_p[odd][j][n]+ VECT_LEN * ( LexSurface( coord_knn,n%4)/2 ); 
        *( chi_nl[odd] + 2 * non_local_count_3[j][odd] + 2*offset +1 )
 = ( Float * ) ( VECT_LEN2 * (LexVector( coord ) / 2*2*NUM_DIR+n+8 )  * sizeof(Float));
	non_local_count_3[j][odd]++;
       }
    }

	      		       
          }  // end of chi3_off-node
    }
    }
    }
    }
  }

  for(int i = 0;i<2;i++)
  for(int j = 0;j<3;j++)
    if (non_local_count_3[j][i] != (non_local_chi_3[i][j+1]-non_local_chi_3[i][j]) ){
    printf("non_local_count_3[%d][%d]=%d\n",i,j,non_local_count_3[j][i]);
    printf("non_local_chi_3[%d][%d]=%d\n",i,j+1,non_local_chi_3[i][j+1]);
    printf("non_local_chi_3[%d][%d]=%d\n",i,j,non_local_chi_3[i][j]);
    exit(-4);
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

  if ( !CoordkNN( n, 3 )  || NP[n%4] <2 ) {		//  chi3_off-node
	      //pointer to source field
	      *( chi_l[ odd ]  +  2 * (local_count_3[ odd ] + local_chi[odd]) )
		= ( Float * ) ( VECT_LEN * ( LexVector( coord_knn ) / 2 )
		                 * sizeof(Float));
	      // pointer to temporary field where U*chi is stored
	      *( chi_l[ odd ]  +  2 * (local_count_3[ odd ] + local_chi[odd]) + 1) =
	  	( Float * ) ( VECT_LEN2 * (LexVector( coord ) / 2 *2*NUM_DIR+n+8) * sizeof(Float));
	      local_count_3[odd]++;
	    }  //else-on_node case


	  } // for x[0] loop
	}
      }
    }// for x[3] loop
  }//for n loop

  for(int i = 0;i<2;i++)
  if (local_count_3[i] != local_chi_3[i]){
    printf(" local_count_3[%d]=%d\n",i,local_count_3[i]);
    exit(-4);
  }

  //-------------------------------------------------------------------
  //  Set up SCU buffer parameters.  T direction is special, since
  //  the block-strided move will not work here.
  //-------------------------------------------------------------------

//  blklen[0] = VECT_LEN * sizeof(Float) * size[1] * size[2] * size[3] / 2;
  blklen[1] = VECT_LEN * sizeof(Float) * size[0] / 2;
  blklen[2] = VECT_LEN * sizeof(Float) * size[0] * size[1] / 2;
  blklen[3] = VECT_LEN * sizeof(Float) * size[0] * size[1] * size[2] / 2;

  numblk[0] = 1;
  numblk[1] = size[2] * size[3];
  numblk[2] = size[3];
  numblk[3] = 1;

  stride[0] = 0;
  stride[1] = (VECT_LEN * size[0] * ( size[1] - 1 ) / 2)*sizeof(Float);
  stride[2] = (VECT_LEN * size[0] * size[1] * ( size[2] - 1 ) / 2)* sizeof(Float) ;
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
      if ( NP[0] >1 && coord[0] == j ) {
        *( ToffsetM[j][ odd ] + countM[j][ odd ] ) = VECT_LEN * ( sc / 2 );
        countM[j][ odd ]++;
      }
      if ( NP[0] >1 && coord[0] == size[0] - 1 -j ) {
        *( ToffsetP[j][ odd ] + countP[j][ odd ] ) = VECT_LEN * ( sc / 2 );
        countP[j][ odd ]++;
      }
    }//end of j loop
  } //end of sg loop
  for(int i = 0;i<2;i++)
  for(int j = 0;j<3;j++){
//   printf("countM[%d][%d] = %d ",j,i,countM[j][i]);
//   printf("countP[%d][%d] = %d\n",j,i,countP[j][i]);
    if (NP[0]>1 && (countM[j][i] != area[(i+j)%2][0])&&(size[0]>j*2)) {printf("countM\n"); exit(-4);}
    if (NP[0]>1 && (countP[j][i] != area[(i+1+j)%2][0])&&(size[0]>j*2)){printf("countP\n"); exit(-4);}
    int pad_len = 4-(countM[j][i]%4);
    if (pad_len==4) pad_len=0;
    int *last = ToffsetM[j][i]+ countM[j][i]-1;
    for(int pad =0;pad<pad_len;pad++) *(last+pad+1) = *last;
    countM[j][i] +=pad_len;

    pad_len = 4-(countP[j][i]%4);
    if (pad_len==4) pad_len=0;
    last = ToffsetP[j][i]+ countP[j][i]-1;
    for(int pad =0;pad<pad_len;pad++) *(last+pad+1) = *last;
    countP[j][i] +=pad_len;
//    printf("countM[%d][%d] = %d ",j,i,countM[j][i]);
//    printf("countP[%d][%d] = %d\n",j,i,countP[j][i]);
  }
  

  //-------------------------------------------------------------------
  //  Index i says data has been received from TP, XP, YP, ZP, TM, XM,
  //  YM, ZM
  //-------------------------------------------------------------------

  int frm_len = VECT_LEN * sizeof(Float);

  for( odd=0;odd<2;odd++){
    int comms=0;
    for ( i = 0; i < NUM_DIR; i++ ) {
    j = i % (NUM_DIR/2);
    if(NP[j]>1){
      int buf_len = frm_len*vol/(2*size[j]) ;
      if (i==0) buf_len = frm_len*area[(1+odd)%2][0];
      if (i==4) buf_len = frm_len*area[(odd)%2][0];
      SCUarg[odd][2*comms].Init(chi_off_node[odd][2][i], scudir[i], SCU_REC,
  	    buf_len, 1, 0, scu_irs[odd][2]);

      if (i==0) buf_len = frm_len*area[(1+odd)%2][0];
      if (i==4) buf_len = frm_len*area[(odd)%2][0];
      SCUDMAarg[odd][(i+NUM_DIR)*2].Init(chi_off_node[odd][0][i],
        buf_len, 1, 0);

      if (i==0) buf_len = frm_len*area[(odd)%2][0];
      if (i==4) buf_len = frm_len*area[(1+odd)%2][0];
      SCUDMAarg[odd][(i+NUM_DIR)*2+1].Init(chi_off_node[odd][1][i],
        buf_len, 1, 0);
      SCUDMAInst *temp[2];
      temp[0] = &SCUDMAarg[odd][(i+NUM_DIR)*2];
      temp[1] = &SCUDMAarg[odd][(i+NUM_DIR)*2+1];
      if( split ){
        SCUarg_1[odd][2*comms].Init(scudir[i],SCU_REC, 
          temp,1, scu_irs[odd][0]);
        SCUarg_2[odd][2*comms].Init(scudir[i],SCU_REC, 
          temp+1,1,scu_irs[odd][1]);
      } else {
        SCUarg_1[odd][2*comms].Init(scudir[i],SCU_REC, 
          temp,2, scu_irs[odd][0]);
      }
  
      //send arguments
      if ((i == 0) || ( i == 4)){
        int buf_len;
        if (i==0) buf_len = frm_len*area[odd%2][0];
        else buf_len = frm_len*area[(1+odd)%2][0];
        if (size[j] >4)
          SCUarg[odd][2*comms+1].Init (Tbuffer[2][(4 - i)/4], scudir[i], SCU_SEND,
  		       buf_len, 1, 0,scu_irs[odd][2] );
        else if (size[j] >2) 
          SCUarg[odd][2*comms+1].Init (Tbuffer[1][i/4], scudir[i], SCU_SEND,
  		       buf_len, 1, 0,scu_irs[odd][2] );
        else
          SCUarg[odd][2*comms+1].Init (chi_off_node[odd][0][4-i], scudir[i], SCU_SEND,
  		       buf_len, 1, 0,scu_irs[odd][2] );
  
        if (i==0) buf_len = frm_len*area[odd%2][0];
        else buf_len = frm_len*area[(1+odd)%2][0];
        SCUDMAarg[odd][i*2].Init(Tbuffer[0][(4 - i)/4], 
  		       buf_len, 1, 0);

        if (i==0) buf_len = frm_len*area[(1+odd)%2][0];
        else buf_len = frm_len*area[(odd)%2][0];
        if(size[j]>2)
          SCUDMAarg[odd][i*2+1].Init(Tbuffer[1][(4 - i)/4], 
  		       buf_len, 1, 0);
        else
          SCUDMAarg[odd][i*2+1].Init(Tbuffer[0][i/4], 
  		       buf_len, 1, 0);
  
        SCUDMAInst *temp[2];
        temp[0]=&(SCUDMAarg[odd][i*2]);
        temp[1]=&(SCUDMAarg[odd][i*2+1]);
        if( split ){
          SCUarg_1[odd][2*comms+1].Init(scudir[i],SCU_SEND,temp,1,scu_irs[odd][0]);
          SCUarg_2[odd][2*comms+1].Init(scudir[i],SCU_SEND,temp+1,1,scu_irs[odd][1]);
        } else {
          SCUarg_1[odd][2*comms+1].Init(scudir[i],SCU_SEND,temp,2,scu_irs[odd][0]);
        }
      }
      else{
        if(size[j] >2)
//
// put in chi_off_node_total to pass communicable test. CJ
//
          SCUarg[odd][2*comms+1].Init(chi_off_node_total, scudir[i], SCU_SEND,
  		       blklen[j], numblk[j], stride[j], scu_irs[odd][2] );
        else
          SCUarg[odd][2*comms+1].Init(chi_off_node[odd][0][(i+4)%NUM_DIR], scudir[i], SCU_SEND,
             VECT_LEN * sizeof(Float) * vol / ( 2 * size[j]),1,0, scu_irs[odd][2] );

        SCUDMAarg[odd][i*2].Init(chi_off_node_total, 
  		       blklen[j], numblk[j], stride[j]);
        SCUDMAarg[odd][i*2+1].Init(chi_off_node_total, 
  		       blklen[j], numblk[j], stride[j]);
        SCUDMAInst *temp[2];
        temp[0]=&SCUDMAarg[odd][i*2];
        temp[1]=&SCUDMAarg[odd][i*2+1];
        if( split ){
          SCUarg_1[odd][2*comms+1].Init(scudir[i],SCU_SEND,temp,1,scu_irs[odd][0]);
          SCUarg_2[odd][2*comms+1]. Init(scudir[i],SCU_SEND,temp+1,1,scu_irs[odd][1]);
        } else {
          SCUarg_1[odd][2*comms+1].Init(scudir[i],SCU_SEND,temp,2,scu_irs[odd][0]);
        }
      }
      comms++;
    }//end of if NP[j]>1
    }// end of NUM_DIR loop
    if (comms!=(2*non_local_dirs)) exit(-55);
  
    SCUDirArgIR *SCUarg_p[2*NUM_DIR];
  
    if(comms){
      for(i = 0;i<2*comms;i++) SCUarg_p[i] = &(SCUarg[odd][i]);
      SCUmulti[odd].Init(SCUarg_p, 2*comms);
      for(i = 0;i<2*comms;i++) SCUarg_p[i] = &(SCUarg_1[odd][i]);
      SCUmulti_1[odd].Init(SCUarg_p, 2*comms);
      if( split ){
        for(i = 0;i<2*comms;i++) SCUarg_p[i] = &(SCUarg_2[odd][i]);
        SCUmulti_2[odd].Init(SCUarg_p, 2*comms);
      }
    }
  } // end of odd loop

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

void AsqD::destroy_buf()
{
  int i,k;

  qfree(chi_off_node_total);

  for ( i = 0; i < 2; i++ ) {
    Free(chi_nl[i]);
    Free(chi_l[i]);
    for ( k = 0; k < 3; k++ ) {
      qfree(Tbuffer[k][i]);
      qfree(ToffsetP[k][i]);
      qfree(ToffsetM[k][i]);
    }
  }
    
}
//-------------------------------------------------------------------
//  Given a lexical value for gauge fields, set the coordinates.
//  Return 1 if odd, 0 if even
//-------------------------------------------------------------------

int AsqD::SetCoord( int sg )
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
void AsqD::init_g(Float *frm_p,Float **fat_p,Float **naik_p, Float **naikm_p)
{

  int i,j,m,n;
  int local_count[2];
  int non_local_count[2];
  int local_count_3[2];
  int non_local_count_3[3][2];
  int x[NUM_DIR/2];
  char *fname = "asqtad_dirac_init_g()";
//  printf("%s::%s\n",cname,fname);

  address[0] = NULL;
  address[1] = NULL;

  //--------------------------------------------------------------------
  // c1 -> one link; c2 -> 3-link; c3 -> 3-link staple; c5 -> 5-link staple;
  // c7 -> 7-link staple; c6 -> 5-link "straight" staple
  //--------------------------------------------------------------------
 
//  printf("fat_p=%p naik_p=%p naikm_p=%p\n",fat_p,naik_p,naikm_p);
  for(int i =0;i<4;i++){
    if(fat_p) fat[i] = (matrix *)fat_p[i];
    if(naik_p) naik[i] = (matrix *)naik_p[i];
    if (naikm_p) naik_m[i] = (matrix *)naikm_p[i];
//    printf("fat=%p naik=%p naik_m=%p\n",fat[i],naik[i],naik_m[i]);
  }
  frm_tmp = frm_p;
//  printf("frm_tmp=%p\n",frm_tmp);
 
  //-----------------------------------------------------------
  //  If t + x + y + z is odd, odd = 1.  Otherwise it is 0.
  //-----------------------------------------------------------

  int odd;

  //-----------------------------------------------------------
  //  SCU transfer structure to get links from off node and a
  //  location where one link matrix can be stored.
  //-----------------------------------------------------------

  //-----------------------------------------------------------------
  //  Allocate space for two copies of the gauge fields on this node
  //-----------------------------------------------------------------
  for ( i = 0; i < 2; i++ ){
    long long buf_len = local_chi[i] + local_chi_3[i];
    buf_len *= MATRIX_SIZE;

    uc_l[i]  = (Float*)Alloc(buf_len*sizeof(Float));
    if(uc_l[i] == 0){
       PointerErr(cname,fname, "uc_l[i]"); exit(3);
    }

    buf_len = PAD(non_local_chi[i]+non_local_chi_3[i][3]);
    buf_len *= MATRIX_SIZE;

    uc_nl[i] = (Float*)Alloc(buf_len*sizeof(Float));
    if(uc_nl[i] == 0){
       PointerErr(cname,fname, "uc_nl[i]");
    }
  }


{
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

  matrix *tmp;
  matrix * Fat;
  matrix *nl[2]; 
  matrix *l[2]; 
  for(i = 0;i<2;i++){
    nl[i] = (matrix *)(uc_nl[i]);
    l[i] = (matrix *)(uc_l[i]);
  }
  for ( n = 0; n < NUM_DIR/2; n++ ) {
    Fat = (matrix *)fat[n];
    int coor=0;
    for (coord[3] = 0; coord[3] < size[3]; coord[3]++)
    for (coord[2] = 0; coord[2] < size[2]; coord[2]++)
    for (coord[1] = 0; coord[1] < size[1]; coord[1]++)
    for (coord[0] = 0; coord[0] < size[0]; coord[0]++){
      odd = ( coord[0] + coord[1] + coord[2] + coord[3] ) % 2;
      int off_node = CoordNN(n);
      if ( NP[n]>1 && off_node ) {		// chi(x+mu) off-node
        *(nl[odd]) = Fat[LexGauge(coord)]; nl[odd]++;
        non_local_count[odd]++;
      } else {
        *(l[odd]) = Fat[LexGauge(coord)]; l[odd]++;
        local_count[odd]++;
      }
      coor++;
    }
  }
  Float *rcv_mat = (Float *)qalloc(QFAST|QNONCACHE,sizeof(matrix));
  if(!rcv_mat) PointerErr(cname,fname, "rcv_mat");
  SCUDir snd_dirs[]={SCU_TP,SCU_XP,SCU_YP,SCU_ZP};
  SCUDir rcv_dirs[]={SCU_TM,SCU_XM,SCU_YM,SCU_ZM};
  sys_cacheflush(0);
  for ( n = 0; n < NUM_DIR/2; n++ ) {
    Fat = (matrix *)fat[n];
//    printf("Fat=%p\n",Fat);
    SCUDirArgIR snd;
    SCUDirArgIR rcv;
    if(NP[n]>1){
      snd.Init(Fat,snd_dirs[n],SCU_SEND,sizeof(matrix));
      rcv.Init(rcv_mat,rcv_dirs[n],SCU_REC,sizeof(matrix));
    }
    for (x[3] = 0; x[3] < size[3]; x[3]++)
    for (x[2] = 0; x[2] < size[2]; x[2]++)
    for (x[1] = 0; x[1] < size[1]; x[1]++)
    for (x[0] = 0; x[0] < size[0]; x[0]++){
      for (i = 0; i < 4 ; i++) coord[i] = x[i];
      odd = ( coord[0] + coord[1] + coord[2] + coord[3] ) % 2;
      int off_node = CoordNN(n+4);
      int *y = coord;
//      printf("coord = %d %d %d %d\n",y[0],y[1],y[2],y[3]);
      y = coord_nn;
//      printf("coord_nn = %d %d %d %d\n",y[0],y[1],y[2],y[3]);
      if ( NP[n]>1 && off_node) {		// chi(chi-mu) off-node
        snd.Addr(Fat+LexGauge(coord_nn));
		snd.StartTrans();rcv.StartTrans();
		snd.TransComplete();rcv.TransComplete();
        nl[odd] ->Dagger((const Float *)rcv_mat ); 
	(nl[odd])->Negate();
        nl[odd]++;
        non_local_count[odd]++;
      } else {
        l[odd] ->Dagger((const Float *)(Fat+LexGauge(coord_nn)) ); 
	    (l[odd])->Negate(); l[odd]++;
        local_count[odd]++;
      }
    }
  } // n
  for(i = 0;i<2;i++)
  if (local_count[i] != local_chi[i]){
    printf(" local_count[%d]=%d\n",i,local_count[i]);
    exit(-4);
  }
  for(i = 0;i<2;i++)
  if (non_local_count[i] != non_local_chi[i]){
    printf(" non_local_count[%d]=%d\n",i,non_local_count[i]);
    exit(-4);
  }

  matrix * Naik;
  for ( n = 0; n < NUM_DIR/2; n++ ) {
    Naik = (matrix *)naik[n];
    for (x[3] = 0; x[3] < size[3]; x[3]++)
    for (x[2] = 0; x[2] < size[2]; x[2]++)
    for (x[1] = 0; x[1] < size[1]; x[1]++)
    for (x[0] = 0; x[0] < size[0]; x[0]++){
   	  for (i = 0; i < 4 ; i++) coord[i] = x[i];
	  odd = ( coord[0] + coord[1] + coord[2] + coord[3] ) % 2;
          int off_node = CoordkNN( n,3 );
      if (NP[n]>1&& off_node ){ 	// chi(x+mu) off-node
        for(int j=0;j<3;j++)
          if(coord_knn[n%4]==j){
            tmp = (matrix *)(uc_nl[odd] + MATRIX_SIZE * 
              (non_local_count_3[j][odd]+non_local_chi_3[odd][j]+
              non_local_chi[odd]));
            *tmp = Naik[LexGauge(coord)]; 
            non_local_count_3[j][odd]++;
          }
      } else {
        tmp = (matrix *)(uc_l[odd] + MATRIX_SIZE * (local_chi[odd]+local_count_3[odd]));
        *tmp = Naik[LexGauge(coord)]; 
        local_count_3[odd]++;
      }
    }
  }

  if (naik_m[0]){
  for ( n = 0; n < NUM_DIR/2; n++ ) {
    Naik = (matrix *)naik_m[n]; 
    for (x[3] = 0; x[3] < size[3]; x[3]++)
    for (x[2] = 0; x[2] < size[2]; x[2]++)
    for (x[1] = 0; x[1] < size[1]; x[1]++)
    for (x[0] = 0; x[0] < size[0]; x[0]++){
   	  for (i = 0; i < 4 ; i++) coord[i] = x[i];
	  odd = ( coord[0] + coord[1] + coord[2] + coord[3] ) % 2;
          int off_node = CoordkNN( n+4,3 );
      if (NP[n]>1 && off_node ){ 	// chi(x+mu) off-node
        for(int j=0;j<3;j++)
          if(coord_knn[n%4]==(size[n%4]-1-j)){
            tmp = (matrix *)(uc_nl[odd] + MATRIX_SIZE * 
              (non_local_count_3[j][odd]+non_local_chi_3[odd][j]+
              non_local_chi[odd]));
            *tmp = Naik[LexGauge(coord)]; 
            non_local_count_3[j][odd]++;
          }
      } else {
        tmp = (matrix *)(uc_l[odd] + MATRIX_SIZE * 
          (local_chi[odd]+ local_count_3[odd]));
        *tmp = Naik[LexGauge(coord)]; 
        local_count_3[odd]++;
      }
    }
  }
  } else {
//    printf("USING NAIK INSTEAD OF NAIK_M\n");
    for(i = 0;i<4;i++)
    if (size[i]<3) {
      printf("Asqd::size[%d](%d) <3\n",i,size[i]);
      exit(13);
    }
    sys_cacheflush(0);
    for ( n = 0; n < NUM_DIR/2; n++ ) {
      Naik = (matrix *)naik[n];
      SCUDirArgIR snd_naik;
      SCUDirArgIR rcv_naik;
      if(NP[n]>1){
      snd_naik.Init(Naik,snd_dirs[n],SCU_SEND,sizeof(matrix));
      rcv_naik.Init(rcv_mat,rcv_dirs[n],SCU_REC,sizeof(matrix));
      }
      for (x[3] = 0; x[3] < size[3]; x[3]++)
      for (x[2] = 0; x[2] < size[2]; x[2]++)
      for (x[1] = 0; x[1] < size[1]; x[1]++)
      for (x[0] = 0; x[0] < size[0]; x[0]++){
//        printf("dir = %d x = %d %d %d %d\n", n,x[0],x[1],x[2],x[3]);
     	for (i = 0; i < 4 ; i++) coord[i] = x[i];
  	    odd = ( coord[0] + coord[1] + coord[2] + coord[3] ) % 2;
            int off_node = CoordkNN( n+4,3);
        if ( NP[n]>1 && off_node ) {		// chi(chi-mu) off-node
          for(int j=0;j<3;j++)
          if(coord_knn[n%4]==(size[n%4]-1-j)){
            snd_naik.Addr(Naik+LexGauge(coord_knn));
            snd_naik.Assert();rcv_naik.Assert();
            snd_naik.StartTrans();rcv_naik.StartTrans();
  	  	    snd_naik.TransComplete();rcv_naik.TransComplete(100);
            tmp = (matrix *)(uc_nl[odd] + MATRIX_SIZE * 
              (non_local_count_3[j][odd]+(non_local_chi_3[odd][j]+
              non_local_chi[odd])));
#if 0
            printf("non_local_count_3[%d][%d](%d) uc_nl=%p tmp=%p\n",j,odd,non_local_count_3[j][odd], uc_nl[odd],tmp);
            printf("non_local_chi_3[%d][%d]=%d non_local_chi[%d]=%d\n",
            odd,j,non_local_chi_3[odd][j],odd,non_local_chi[odd]);
#endif
            tmp->Dagger((const Float *)rcv_mat ); 
  	        tmp->Negate();
            non_local_count_3[j][odd]++;
          }
        } else {
          tmp = (matrix *)(uc_l[odd] + MATRIX_SIZE * (local_chi[odd]+local_count_3[odd]));
          tmp->Dagger((const Float *)(Naik+LexGauge(coord_knn)) ); 
  	      tmp->Negate();
          local_count_3[odd]++;
        }
      }
    }
  }
  for(int i = 0;i<2;i++)
  if (local_count_3[i] != local_chi_3[i]){
    printf("local_count_3[%d]=%d\n",i,local_count_3[i]);
     exit(-4);
  }
  for(int j = 0;j<3;j++)
  for(int i = 0;i<2;i++)
  if (non_local_count_3[j][i] != (non_local_chi_3[i][j+1]-non_local_chi_3[i][j]) ){
    printf("non_local_count_3[%d][%d]=%d\n",i,j,non_local_count_3[j][i]);
    printf("non_local_chi_3[%d][%d]=%d\n",i,j+1,non_local_chi_3[i][j+1]);
    printf("non_local_chi_3[%d][%d]=%d\n",i,j,non_local_chi_3[i][j]);
     exit(-4);
  }
  qfree (rcv_mat);

}

  tmpfrm = NULL;
//  if (vol<1296)
  tmpfrm = (Float *) qalloc (QFAST|QCOMMS,NUM_DIR*2 * vol/2 * VECT_LEN2 * sizeof(Float));
  if(tmpfrm == NULL){ 
    tmpfrm = (Float *) qalloc (QCOMMS,NUM_DIR*2 * vol/2 * VECT_LEN2 * sizeof(Float));
    printf("tmpfrm is allocated at DDR(%p),length 0x%x \n",tmpfrm,NUM_DIR*vol*VECT_LEN2*sizeof(Float));
  }
  if(tmpfrm == 0) 
    PointerErr(cname,fname, "tmpfrm");
  if (qalloc_is_fast(tmpfrm) && (sizeof(Float)==4))  //single precision
    tmpfrm2 = (unsigned long)tmpfrm - 0xb0000000 + 0x9c000000;
  else
    tmpfrm2 = (unsigned long)tmpfrm;
//  printf("tmpfrm=%p tmpfrm2=%p\n",tmpfrm,tmpfrm2);

  bzero( (char *)tmpfrm, NUM_DIR*2 * vol/2 * VECT_LEN2 * sizeof(Float));
  for(i=0;i<2;i++){
    int buf_len = PAD(local_chi[i]+local_chi_3[i]);
    buf_len *= sizeof(gauge_agg);
    if(vol> 4096) uc_l_agg[i]=NULL;
    else
      uc_l_agg[i]  = (gauge_agg *)qalloc(QFAST,buf_len);
    if(uc_l_agg[i] == NULL){
      uc_l_agg[i]  = (gauge_agg *)qalloc(QCOMMS,buf_len);
      printf("uc_l_agg[%d] is allocated at DDR (%p)\n",i,uc_l_agg[i]);
    }
    buf_len = PAD(non_local_chi[i]+non_local_chi_3[i][3]);
    buf_len *= sizeof(gauge_agg);
    if(uc_l_agg[i] == 0)
      PointerErr(cname,fname, "uc_l_agg[i]");
    uc_nl_agg[i]  = (gauge_agg*)qalloc(QFAST,buf_len);
    if(uc_nl_agg[i] == 0){
      uc_nl_agg[i]  = (gauge_agg*)qalloc(QCOMMS,buf_len);
      printf("uc_nl_agg[%d] is allocated at DDR (%p)\n",i,uc_nl_agg[i]);
    }
    if(uc_nl_agg[i] == 0){
      PointerErr(cname,fname, "uc_nl_agg[i]");
    }
  }

  int temp_len = NUM_DIR*vol;
  if ( (PAD(non_local_chi[0])*6)>temp_len) temp_len = PAD(non_local_chi[0])*6;

  gauge_agg *temp = (gauge_agg *)Alloc(temp_len*sizeof(gauge_agg));
//  printf("temp=%p\n",temp);
  int num_ind[NUM_DIR*vol];
  int src;
  for(j=0;j<2;j++){
    for(i=0;i<vol;i++) num_ind[i]=0;
    for(i=0;i< (local_chi[j]+ local_chi_3[j]);i++){
      src = (int)chi_l[j][2*i];
      if (src%(VECT_LEN*sizeof(Float))!=0){
        printf("%s::%s: src = %d\n",cname,fname,src);
        exit(1);
      }
      src= src/(VECT_LEN*sizeof(Float));
      if(src > vol/2) {
        printf("%s::%s: src[%d](%d) > vol/2\n",cname,fname,i,src);
        exit(1);
      }
      int dest = (int)chi_l[j][2*i+1]/(16*VECT_LEN2*sizeof(Float));
      if(dest > vol/2) {
        printf("%s::%s: dest_l[%d](%d) > vol/2\n",cname,fname,i,dest);
        exit(1);
      }
      if ((src*NUM_DIR*2+num_ind[src])>temp_len){
        printf("%s::%s: index[%d](%d) > NUM_DIR*vol\n",cname,fname,i,
        src*NUM_DIR*2+num_ind[src]);
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
#ifndef NEW_SPLIT
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
//      if (n==0) odd_num +=num_ind[i]%div;
    }
#else
    for( n = 0; n*div<NUM_DIR*2;n++)
    for(i=0;i<vol/2;i++){
      if( num_ind[i] >= (n+1)*div )
      for(m=0;m<div;m++){
           uc_l_agg[j][index].src = temp[i*NUM_DIR*2+n*div+m].src;
           uc_l_agg[j][index].dest = temp[i*NUM_DIR*2+n*div+m].dest;
           for(k=0;k<18;k++)
              uc_l_agg[j][index].mat[k] = temp[i*NUM_DIR*2+n*div+m].mat[k];
           index++;
      }
    }
    even_l[j] = index;
//    printf("even_l[%d]=%d\n",j,even_l[j]);
    for( n = 0; n*div<NUM_DIR*2;n++)
    for(i=0;i<vol/2;i++){
      if( (num_ind[i] > n*div) && (num_ind[i] <(n+1)*div) )
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
    odd_l[j] = index-even_l[j];
//    printf("odd_l[%d]=%d\n",j,odd_l[j]);
#endif
    if (index != comp_l[j]){
      printf("index(%d) = comp_l[%d](%d)\n",
       index, j,comp_l[j]);
      exit(1);
    }
  }


  for(j=0;j<2;j++){
    for(i=0;i<PAD(non_local_chi[j])*3;i++) num_ind[i]=0;
    for(i=0;i< ((non_local_chi[j]+ non_local_chi_3[j][3]));i++){
      src = (int)chi_nl[j][2*i];
      if (src%(VECT_LEN*sizeof(Float))!=0){
        printf("%s::%s: src = %d\n",cname,fname,src);
        exit(1);
      }
      src = src/(VECT_LEN*sizeof(Float));
      if(src > PAD(non_local_chi[0])*3) {
        printf("%s::%s: src(%d) > non_local_chi*3\n",cname,fname,src);
        exit(1);
      }
      int dest = (int)chi_nl[j][2*i+1]/(16*VECT_LEN2*sizeof(Float));
      if(dest > vol/2) {
        printf("%s::%s: dest_nl[%d](%d) > vol/2\n",cname,fname,i,dest);
        exit(1);
      }
      if ((src*2+num_ind[src])>temp_len){
        printf("%s::%s: index_nl[%d](%d) > NUM_DIR*vol\n",cname,fname,i,
        src*2+num_ind[src]);
        exit(1);
      }
      temp[src*2+num_ind[src]].src = (int)chi_nl[j][2*i];
      temp[src*2+num_ind[src]].dest= (int)chi_nl[j][2*i+1];
      for(k=0;k<18;k++)
        temp[src*2+num_ind[src]].mat[k] = uc_nl[j][i*18+k];
      num_ind[src]++;
      if(num_ind[src]>2){
        printf("%s::%s: num_ind[%d](%d) > 2 \n",cname,fname,src, num_ind[src]);
        exit(1);
      }
    }

    int index = 0;
    int div = 2;
    for( n=0; n*div<2; n++)
    for(i=0;i<PAD(non_local_chi[0])*3;i++){
      for(m=0;m<div;m++){
        if(num_ind[i] > n*div+m) {
           uc_nl_agg[j][index] = temp[i*2+n*div+m];
           index++;
        }
      }
    }
    if (index != comp_nl[j] + comp_nl_2[j]){
      printf("index(%d) = comp_nl[%d](%d) + comp_nl_2[%d](%d)\n",
       index, j,comp_nl[j],j,comp_nl_2[j]);
      exit(1);
    }
  }

  for ( i = 0; i < 2; i++){
    if(comp_l[i]%2 ==1){
      uc_l_agg[comp_l[i]] = uc_l_agg[comp_l[i]-1];
      comp_l[i]++;
    }
    int total = comp_nl[i] + comp_nl_2[i];
    if(total%2 == 1){
      uc_nl_agg[total] = uc_nl_agg[total-1];
      comp_nl_2[i]++;
    }
    if(comp_nl[i]%2==1){
      comp_nl[i]--;
      comp_nl_2[i]++;
    }
  }
//  delete[] temp;
  Free(temp);

  for ( i = 0; i < 2; i++){
  Free(uc_l[i]);
  Free(uc_nl[i]);
  }
//  printf("init_g() done\n");
  
}

extern "C"
void AsqD::destroy_buf_g(void)
{
  int i;
  for ( i = 0; i < 2; i++){
  qfree(uc_l_agg[i]);
  qfree(uc_nl_agg[i]);
  }
  qfree (tmpfrm);
}
//---------------------------------------------------------------------
//  Find nearest neighbor coordinate for coordinates given.  Nearest
//  neighbor coordinates are placed in coord_nn, which are always
//  on-node.  Function returns 0 if nearest neighbor is on-node,
//  1 if off-node.
//
//  nn = 0 to 7.  tp, xp, yp, zp, tm, xm, ym, zm respectively.
//---------------------------------------------------------------------

int AsqD::CoordNN( int nn )
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

int AsqD::CoordkNN( int nn, int k )
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

int AsqD::LexGauge( int * c )
{
  return c[1] + size[1] * ( c[2] + size[2] * ( c[3] + size[3] * c[0] ));
}

//---------------------------------------------------------------------
//  Return lexical value for vectors from coordinates c
//---------------------------------------------------------------------

int AsqD::LexVector( int * c )
{
  return c[0] + size[0] * ( c[1] + size[1] * ( c[2] + size[2] * c[3] ));
}

//-------------------------------------------------------------------
//  Given a coordinate and a surface ( 0 = t, 1 = x, 2 = y, 3 = x )
//  calculate the offset into the receive buffers (s?) .
//-------------------------------------------------------------------

int AsqD::LexSurface( int * cc, int surface )
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

  int temp =  c[0] + s[0] * ( c[1] + s[1] * ( c[2] + s[2] * c[3] ));
  return temp;
}

void AsqD::comm_assert()
{
  int i,odd;
  for(odd=0;odd<2;odd++)
  for(i=0;i<2*non_local_dirs;i++){
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



#undef PROFILE
void AsqD::dirac(Float* b, Float* a, int a_odd, int add_flag)
{
#ifdef PROFILE
  int num_flops;
  Float dtime;  
  struct timeval start,end;
  
  num_flops = 0;
  gettimeofday(&start,NULL);
#endif

  int i;
  int odd = 1 - a_odd;
  long c = (long) chi_off_node_total;
  long uc_l = (long)uc_l_agg[odd];
  long uc_nl = (long)uc_nl_agg[odd];
  long uc_nl2 = (long)&(uc_nl_agg[odd][isplit]);

  if( (unsigned long)a != address[odd]){
    address[odd] = (unsigned long)a;
  
    int comms=0;
    if (NP[0]>1) comms++;
    for(i=1;i<4;i++){
      if(size[i] >2 ){
        SCUarg[odd][2*comms+1].Addr( a + Xoffset[2][i]);
        SCUarg[odd][2*(comms+non_local_dirs)+1].Addr( a + Xoffset[2][i+4]);
      }
      if (NP[i]>1) comms++;
    }
  
    void *addr[2];
    if(split) {
      comms=0;
      if (NP[0]>1) comms++;
      for(i=1;i<4;i++){
        addr[0] = (void *)(a+Xoffset[0][i]);
        addr[1] = (void *)(a+Xoffset[1][i]);
        SCUarg_1[odd][2*comms+1].Addr( addr,1);
        SCUarg_2[odd][2*comms+1].Addr( addr+1,1);
        addr[0] = (void *)(a+Xoffset[0][i+4]);
        addr[1] = (void *)(a+Xoffset[1][i+4]);
        SCUarg_1[odd][2*(comms+non_local_dirs)+1].Addr( addr,1);
        SCUarg_2[odd][2*(comms+non_local_dirs)+1].Addr( addr+1,1);
        comms++;
      }
    } else {
      comms=0;
      if (NP[0]>1) comms++;
      for(i=1;i<4;i++){
        addr[0] = (void *)(a+Xoffset[0][i]);
        addr[1] = (void *)(a+Xoffset[1][i]);
        SCUarg_1[odd][2*comms+1].Addr( addr,2);
        addr[0] = (void *)(a+Xoffset[0][i+4]);
        addr[1] = (void *)(a+Xoffset[1][i+4]);
        SCUarg_1[odd][2*(comms+non_local_dirs)+1].Addr( addr,2);
        comms++;
      }
    }
  } // if address[odd] != a


  //-----------------------------------------------------------------
  //  Transfer chi's on faces.  
  //-----------------------------------------------------------------

#if 0
  for(int i=0;i<2;i++)
  for(int j=0;j<3;j++)
  for(int k=0;k<(countM[j][i]*VECT_LEN);k++)
    if ( fabs(*(Tbuffer[j][i]+k))>1e-10)
    printf("Tbuffer[%d][%d][%d]=%e\n",j,i,k,*(Tbuffer[j][i]+k));
  printf("copy_buffer\n");
#endif
  if (NP[0]>1){
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
  }
  sys_cacheflush(0);
#if 0
  for(int i=0;i<2;i++)
  for(int j=0;j<3;j++)
  for(int k=0;k<(countM[j][i]*VECT_LEN);k++)
    if ( fabs(*(Tbuffer[j][i]+k))>1e-10)
    printf("Tbuffer[%d][%d][%d]=%e\n",j,i,k,*(Tbuffer[j][i]+k));
#endif
#ifdef PROFILE
  gettimeofday(&end,NULL);
  asq_print_flops(cname,"copy_buffer",0,&start,&end);
#endif


  //make sure spinor field is in main memory before starting transfers

#ifdef PROFILE
  gettimeofday(&start,NULL);
#endif
  SCUmulti_1[odd].StartTrans();
#ifdef PROFILE
  gettimeofday(&end,NULL);
  asq_print_flops(cname,"StartTrans()",0,&start,&end);
#endif

//  printf("uc=%p a=%p tmpfrm=%p\n",uc_l_agg[odd],a,tmpfrm);
#ifdef PROFILE
  num_flops = 66*comp_l[odd];
  gettimeofday(&start,NULL);
#endif

//  printf("comp_l=%d comp_nl=%d comp_nl_2=%d\n",comp_l[odd],comp_nl[odd],comp_nl_2[odd]);
  //-----------------------------------------------------------------
  //do first local computations
  //-----------------------------------------------------------------
#ifdef NEW_SPLIT
  asq_cmv_4( even_l[odd], (long)uc_l_agg[odd], (long)a, tmpfrm2);
  if (odd_l[odd])
    asq_cmv( odd_l[odd], (long)(uc_l_agg[odd]+even_l[odd]), (long)a, (long)tmpfrm2);
#else
  asq_cmv( comp_l[odd], (long)uc_l_agg[odd], (long)a, (long)tmpfrm2);
#endif


  //-----------------------------------------------------------------
  // check to see if transfers are done and start another transfer
  //-----------------------------------------------------------------

#ifdef PROFILE
  gettimeofday(&end,NULL);
  asq_print_flops(cname,"asq_cmv",num_flops,&start,&end);
#endif

#ifdef PROFILE
  gettimeofday(&start,NULL);
#endif
  SCUmulti_1[odd].TransComplete();
#ifdef PROFILE
  gettimeofday(&end,NULL);
  asq_print_flops(cname,"TrnasComplete()",0,&start,&end);
#endif

if(split) {

#ifdef PROFILE
  gettimeofday(&start,NULL);
#endif
  SCUmulti_2[odd].StartTrans();
#ifdef PROFILE
  gettimeofday(&end,NULL);
  asq_print_flops(cname,"StartTrans()",0,&start,&end);
#endif


#ifdef PROFILE
  num_flops = 66*isplit;
  gettimeofday(&start,NULL);
#endif

  asq_cmv( isplit, uc_nl, c, (long)tmpfrm2);

#ifdef PROFILE
  gettimeofday(&end,NULL);
  asq_print_flops(cname,"asq_cmv",num_flops,&start,&end);
#endif
  SCUmulti_2[odd].TransComplete();

  SCUmulti[odd].StartTrans();

#ifdef PROFILE
  num_flops = 66*(non_local_chi[odd]);
  gettimeofday(&start,NULL);
#endif

  asq_cmv( non_local_chi[odd], uc_nl2, (long)c, (long)tmpfrm2);

} else {

#ifdef PROFILE
  gettimeofday(&start,NULL);
#endif
  SCUmulti[odd].StartTrans();
#ifdef PROFILE
  gettimeofday(&end,NULL);
  asq_print_flops(cname,"StartTrans()",0,&start,&end);
#endif

  //-----------------------------------------------------------------
  //do the computations involving "chi" non-local spinors
  //-----------------------------------------------------------------


#ifdef PROFILE
  num_flops = 66*comp_nl[odd];
  gettimeofday(&start,NULL);
#endif

  asq_cmv( comp_nl[odd], (long)&(uc_nl_agg[odd][0]), (long)c, (long)tmpfrm2);
}


#ifdef PROFILE
  gettimeofday(&end,NULL);
  asq_print_flops(cname,"asq_cmv",num_flops,&start,&end);
#endif

#ifdef PROFILE
  gettimeofday(&start,NULL);
#endif
  SCUmulti[odd].TransComplete();

#ifdef PROFILE
  gettimeofday(&end,NULL);
  asq_print_flops(cname,"TrnasComplete()",0,&start,&end);
#endif

 

  //-----------------------------------------------------------------
  //do the computations involving chi3 non-local spinors
  //----------------------------------------------------------------


#ifdef PROFILE
  num_flops = 66*comp_nl_2[odd];
  gettimeofday(&start,NULL);
#endif

  asq_cmv( comp_nl_2[odd], (long)&(uc_nl_agg[odd][comp_nl[odd]]) , (long)c, (long)tmpfrm2);

#ifdef PROFILE
  gettimeofday(&end,NULL);
  asq_print_flops(cname,"asq_cmv",num_flops,&start,&end);
#endif


  // check to see if transfers are done

  //-----------------------------------------------------------------
  //do the sum of 16 temporary vectors at each lattice site
  //              ^^^ change must be made in  dirac_sum**
  //-----------------------------------------------------------------


#ifdef PROFILE
  num_flops = 45*vol;
  gettimeofday(&start,NULL);
#endif
  if ( add_flag == 0){
    asq_dsum( vol/2, (long)tmpfrm, (long)b);
  }
  else{
    asqd_sum_acc_cpp( vol/2, (long)0, (long)tmpfrm, (long)b);
  }
#ifdef PROFILE
  gettimeofday(&end,NULL);
  asq_print_flops(cname,"asq_dsum",num_flops,&start,&end);
#endif

}

