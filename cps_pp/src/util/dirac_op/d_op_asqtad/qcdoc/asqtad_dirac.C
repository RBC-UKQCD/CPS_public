//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2005-03-07 00:22:22 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_asqtad/qcdoc/asqtad_dirac.C,v 1.18 2005-03-07 00:22:22 chulwoo Exp $
//  $Id: asqtad_dirac.C,v 1.18 2005-03-07 00:22:22 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: asqtad_dirac.C,v $
//  $Revision: 1.18 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_asqtad/qcdoc/asqtad_dirac.C,v $
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

#include <util/asqtad_int.h>

#if 0
void AsqD::PointerErr(char *cname, char *fname, char *vname){
  printf("%s::%s: %s not allocated\n",cname,fname,vname);
  exit(-1);
}
#endif

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


const char *uc_l_filename = "uc_l.h";
const char *uc_nl_filename = "uc_nl.h";
const char *Toffset_filename = "Toffset.h";
const char *uc_l_agg_filename = "uc_l_agg.h";
const char *uc_nl_agg_filename = "uc_nl_agg.h";
const char *chi_l_filename = "chi_l.h";
const char *chi_nl_filename = "chi_nl.h";

/*****************************************************************
 SIMUL switched on/off the hack CJ put in to help speed up the
simulation. if SIMUL is undefined, program writes the temporary arraies
for dirac operator. If SIMUL is defined, it will include (pre-calcuated)
temporary arraies and skip the generation of arraies.
******************************************************************/

#undef CPP
/****************************************************************
CPP is a switch for using C++ routine for dirac_cmv.
*****************************************************************/

void asqd_sum_acc_cpp(int s, long chi, long tmpfrm, long b);
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

enum{VECT_LEN=6, VECT_LEN2=8, MATRIX_SIZE=18, SITE_LEN=72, NUM_DIR=8, N=4};

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

#if 0
static Fasqtad *lat_pt;
void set_pt (Fasqtad *lat)
{
  lat_pt = lat;
}
#endif
//-------------------------------------------------------------------
//  Called by the lattice constructor
//  Fermion initializations: pointer tables 
//-------------------------------------------------------------------
void AsqD::init(AsqDArg *arg)
{

//  arg = asq_arg;
//  *this = *asq_arg;
  gauge_field_addr = ( Float * ) arg->gauge_u;
  int i,j,m,n; 
  int blklen[NUM_DIR/2];
  int numblk[NUM_DIR/2];
  int stride[NUM_DIR/2];
  int local_count[2];
  int non_local_count[2];
  int local_count_3[2];
  int non_local_count_3[3][2];
  int x[NUM_DIR/2];
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

  for(int i = 0;i<4;i++) size[i] = arg->size[i];
  for(int i = 0;i<4;i++) NP[i] = arg->NP[i];
  for(int i = 0;i<4;i++) coor[i] = arg->coor[i];
  c1 = arg->c1;
  c2 = arg->c2;
  c3 = arg->c3;
  c5 = arg->c5;
  c7 = arg->c7;
  c6 = arg->c6;
  fat = (matrix *)arg->Fat;
  naik = (matrix *)arg->Naik;
  naik_m = (matrix *)arg->NaikM;

  vol = size[0] * size[1] * size[2] * size[3];
  f_size_cb = vol*3;
//VRB.Flow(cname,fname,"vol=%d\n",vol);
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

  //-----------------------------------------------------------------
  //  Allocate 8 receive buffers for off-node vectors
  //-----------------------------------------------------------------
  if(vol>1024) chi_off_node_total=NULL;
    else
    chi_off_node_total = ( Float * ) qalloc(QFAST|QCOMMS, 3*non_local_chi*
      VECT_LEN * sizeof( Float ) / 2 );
  if(chi_off_node_total == NULL){ 
    chi_off_node_total = ( Float * ) qalloc(QCOMMS, 3*non_local_chi*
      VECT_LEN * sizeof( Float ) / 2 );
    printf("chi_off_node_total is allocated at DDR (%p)\n",chi_off_node_total);
  }
    if(chi_off_node_total == 0)
      PointerErr(cname,fname, "chi_off_node_total");

 for ( j= 0; j < 3; j++ ){
    chi_off_node[j][0] = &(chi_off_node_total[ non_local_chi*j* VECT_LEN/2 ] ); 
//    Fprintf(stderr,"chi_off_node[%d][0] =%d\n",j,(chi_off_node[j][0]-chi_off_node_total)/VECT_LEN);
    chi_off_node_p[j][0] = (Float *)(sizeof (Float)*non_local_chi*j* VECT_LEN/2 );
  for ( i = 1; i < NUM_DIR; i++ ){
    chi_off_node[j][i] = chi_off_node[j][i-1]+vol/(2*size[(i-1)%4])*VECT_LEN;
//    Fprintf(stderr,"chi_off_node[%d][%d] =%d\n",j,i,(chi_off_node[j][i]-chi_off_node_total)/VECT_LEN);
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
    chi_l[i] = ( Float ** ) Alloc((2*(local_chi+local_chi_3)/2)*sizeof(Float *));
   if(chi_l[i] == NULL)
	PointerErr(cname,fname, "chi_l[i]");

    chi_nl[i] = (Float ** ) Alloc((2* (non_local_chi + non_local_chi_3[3])/2)*sizeof(Float *));

      if(chi_nl[i] == NULL)
	PointerErr(cname,fname, "chi_nl[i]");
  }
  

  for ( i = 0; i < 2; i++){
    local_count[i] = 0;
    non_local_count[i] = 0;
    local_count_3[i] = 0;
    for( n = 0; n < 3; n++ )  non_local_count_3[n][i] = 0;

  }

  for ( k = 0; k < 3; k++ ) {
  for ( i = 0; i < 2; i++ ) {
  if(vol>1024) Tbuffer[k][i]=NULL;
    else
    Tbuffer[k][i] = (Float *) qalloc (QFAST|QNONCACHE, size[1] * size[2] * size[3] * VECT_LEN * sizeof( Float ) / 2);

  if( Tbuffer[k][i] == NULL)
    Tbuffer[k][i] = (Float *) qalloc (QCOMMS, size[1] * size[2] * size[3] * VECT_LEN * sizeof( Float ) / 2);

   if(Tbuffer[k][i] == NULL)
	PointerErr(cname,fname, "Tbuffer[i][j]");
    ToffsetP[k][i] = ( int * ) qalloc (0,  size[1] * size[2] * size[3] *  sizeof( int ) / 2 );
   if(ToffsetP[k][i] == NULL)
	PointerErr(cname,fname, "TOffsetP[i][j]");
    ToffsetM[k][i] = ( int * ) qalloc (0,  size[1] * size[2] * size[3] *  sizeof( int ) / 2 );
   if(ToffsetM[k][i] == NULL)
	PointerErr(cname,fname, "TOffsetM[i][j]");
    countP[k][i] = 0;
    countM[k][i] = 0;
  }
 }

#if 0
  //-----------------------------------------------------------------
  // Assembly written for double precision only, check sizeof(Float)
  //-----------------------------------------------------------------
  if ( sizeof(Float) != sizeof(double)){
     ERR.General(cname, fname, 
		 "Assembly functions implemented only for double precision!");
  }
#endif


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
				( Float * ) ( VECT_LEN2 * (LexVector( coord ) / 2 * 2*NUM_DIR+n )  * sizeof(Float));
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
 = ( Float * ) ( VECT_LEN2 * (LexVector( coord ) / 2*2*NUM_DIR+n+8 )  * sizeof(Float));
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
		= ( Float * ) ( VECT_LEN * ( LexVector( coord_knn ) / 2 )
		                 * sizeof(Float));
	      // pointer to temporary field where U*chi is stored
	      *( chi_l[ odd ]  +  2 * local_count_3[ odd ] + local_chi + 1) =
	  	( Float * ) ( VECT_LEN2 * (LexVector( coord ) / 2 *2*NUM_DIR+n+8) * sizeof(Float));
	      local_count_3[odd]++;
	    }  //else-on_node case


	  } // for x[0] loop
	}
      }
    }// for x[3] loop
  }//for n loop

#if 0
  FILE *fp;
  fp=Fopen(chi_l_filename,"w");
  for(j=0;j<2;j++){
    Fprintf(fp,"Float * chi_l%d[] LOCATE(\"edramtransient\") = {\n",j);
    Fprintf(fp," (Float *) %d",*(chi_l[j]));
    for(i=1;i< 2*((local_chi+local_chi_3)/2);i++){
      Fprintf(fp,",\n (Float *) %d",*(chi_l[j]+i));
    }
    Fprintf(fp,"\n};\n");
  }
  Fclose(fp);

  fp=Fopen(chi_nl_filename,"w");
  for(j=0;j<2;j++){
    Fprintf(fp,"Float * chi_nl%d[] LOCATE(\"edramtransient\") = {\n",j);
    Fprintf(fp," (Float *) %d",*(chi_nl[j]));
    for(i=1;i< 2*((non_local_chi+ non_local_chi_3[3])/2);i++){
      Fprintf(fp,",\n (Float *) %d",*(chi_nl[j]+i));
    }
    Fprintf(fp,"\n};\n");
  }

  Fclose(fp);
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

  blklen[0] = VECT_LEN * sizeof(Float) * size[1] * size[2] * size[3] / 2;
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
  

#if 0
  int vol3 = (size[1] * size[2] * size[3])/2;
  fp= Fopen(Toffset_filename,"w");
  for ( k = 0; k < 3; k++ ) 
  for ( i = 0; i < 2; i++ ) {
    Fprintf(fp,"countP[%d][%d]=%d\n",k,i,countP[k][i]);
    Fprintf(fp,"int ToffsetP%d%d[] LOCATE(\"edramnormal\") = {\n",k,i);
    for( j = 0;j<vol3;j++){
      Fprintf(fp, "%d,\n",ToffsetP[k][i][j]);
    }
    Fprintf(fp,"};\n");
    Fprintf(fp,"countM[%d][%d]=%d\n",k,i,countM[k][i]);
    Fprintf(fp,"int ToffsetM%d%d[] LOCATE(\"edramnormal\") = {\n",k,i);
    for( j = 0;j<vol3;j++){
      Fprintf(fp, "%d,\n",ToffsetM[k][i][j]);
    }
    Fprintf(fp,"};\n");
  }
  Fclose(fp);
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
  		    VECT_LEN * sizeof(Float) * vol / ( 2 * size[j] ), 1, 0, scu_irs[odd][2]);
//        SCUarg[odd][i + NUM_DIR].Assert();
      SCUDMAarg_p[odd][(i+NUM_DIR)*2]  = new SCUDMAInst;
  //      printf("SCUDMAarg_p[%d]=%p\n",(i+NUM_DIR)*2,SCUDMAarg_p[(i+NUM_DIR)*2]);
      SCUDMAarg_p[odd][(i+NUM_DIR)*2] ->Init(chi_off_node[0][i],
        VECT_LEN * sizeof(Float) * vol / ( 2 * size[j] ), 1, 0);
      SCUDMAarg_p[odd][(i+NUM_DIR)*2+1]  = new SCUDMAInst;
      SCUDMAarg_p[odd][(i+NUM_DIR)*2+1] ->Init(chi_off_node[1][i],
        VECT_LEN * sizeof(Float) * vol / ( 2 * size[j] ), 1, 0);
      if( split ){
        SCUarg_1[odd][i + NUM_DIR].Init(scudir[i],SCU_REC, &SCUDMAarg_p[odd][(i+NUM_DIR)*2],1, scu_irs[odd][0]);
//        SCUarg_1[odd][i + NUM_DIR].Assert();
        SCUarg_2[odd][i + NUM_DIR].Init(scudir[i],SCU_REC, &SCUDMAarg_p[odd][(i+NUM_DIR)*2+1],1, scu_irs[odd][1]);
//        SCUarg_2[odd][i + NUM_DIR].Assert();
      } else {
        SCUarg_1[odd][i + NUM_DIR].Init(scudir[i],SCU_REC, &SCUDMAarg_p[odd][(i+NUM_DIR)*2],2, scu_irs[odd][0]);
//        SCUarg_1[odd][i + NUM_DIR].Assert();
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
//        SCUarg[odd][i].Assert();
  
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
//          SCUarg_1[odd][i].Assert();
          SCUarg_2[odd][i].Init(scudir[i],SCU_SEND,&SCUDMAarg_p[odd][i*2+1],1,scu_irs[odd][1]);
//          SCUarg_2[odd][i].Assert();
        } else {
          SCUarg_1[odd][i].Init(scudir[i],SCU_SEND,&SCUDMAarg_p[odd][i*2],2,scu_irs[odd][0]);
//          SCUarg_1[odd][i].Assert();
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
             VECT_LEN * sizeof(Float) * vol / ( 2 * size[j]),1,0, scu_irs[odd][2] );

//        SCUarg[odd][i].Assert();
  
        SCUDMAarg_p[odd][i*2] = new SCUDMAInst;
        SCUDMAarg_p[odd][i*2] ->Init(chi_off_node_total, 
  		       blklen[j], numblk[j], stride[j]);
        SCUDMAarg_p[odd][i*2+1] = new SCUDMAInst;
        SCUDMAarg_p[odd][i*2+1] ->Init(chi_off_node_total, 
  		       blklen[j], numblk[j], stride[j]);
        if( split ){
          SCUarg_1[odd][i].Init(scudir[i],SCU_SEND,&SCUDMAarg_p[odd][i*2],1,scu_irs[odd][0]);
//          SCUarg_1[odd][i].Assert();
          SCUarg_2[odd][i]. Init(scudir[i],SCU_SEND,&SCUDMAarg_p[odd][i*2+1],1,scu_irs[odd][1]);
//          SCUarg_2[odd][i].Assert();
        } else {
          SCUarg_1[odd][i].Init(scudir[i],SCU_SEND,&SCUDMAarg_p[odd][i*2],2,scu_irs[odd][0]);
//          SCUarg_1[odd][i].Assert();
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
void AsqD::init_g(Float *frm_p,Float *fat_p,Float *naik_p, Float *naikm_p)
{

  int i,j,m,n;
  int local_count[2];
  int non_local_count[2];
  int local_count_3[2];
  int non_local_count_3[3][2];
  int x[NUM_DIR/2];
  char *cname = "DiracOpAsqtad";
  char *fname = "asqtad_dirac_init_g()";
  printf("%s::%s\n",cname,fname);
//  VRB.Func(cname,fname);

  //--------------------------------------------------------------------
  // c1 -> one link; c2 -> 3-link; c3 -> 3-link staple; c5 -> 5-link staple;
  // c7 -> 7-link staple; c6 -> 5-link "straight" staple
  //--------------------------------------------------------------------
 
#if 0
  Float c1 = GJP.KS_coeff();
  Float c2 = GJP.Naik_coeff();
  Float c3 = GJP.staple3_coeff();
  Float c5 = GJP.staple5_coeff();
  Float c7 = GJP.staple7_coeff();
  Float c6 = GJP.Lepage_coeff();
#else
#endif
  printf("fat=%p naik=%p naik_m=%p\n",fat,naik,naik_m);
  if (fat_p) fat = (matrix *)fat_p;
  if (naik_p) naik = (matrix *)naik_p;
  if (naikm_p) naik_m = (matrix *)naikm_p;
  printf("fat=%p naik=%p naik_m=%p\n",fat,naik,naik_m);
  frm_tmp = frm_p;
 

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

//  int sg;

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


  //-----------------------------------------------------------
  //  Once all the index arithmetic is finished, v points to
  //  the initial gauge field matrix.  w points to where it should
  //  be stored.  same for w3 (for UUU matrix).
  //-----------------------------------------------------------

//  Float * v;
//  Float * w, wp1[18];
//  Float * w3;

//  Float  w_t1[18], w_t2[18], w_t3[18];
//  Float mtmp[18] ;

  //-----------------------------------------------------------
  //  SCU transfer structure to get links from off node and a
  //  location where one link matrix can be stored.
  //-----------------------------------------------------------

  //-----------------------------------------------------------------
  //  Allocate space for two copies of the gauge fields on this node
  //-----------------------------------------------------------------
  for ( i = 0; i < 2; i++ ){

    uc_l[i]  = (Float*)Alloc(MATRIX_SIZE*((local_chi+local_chi_3)/2)*sizeof(Float));
     if(uc_l[i] == 0){
       PointerErr(cname,fname, "uc_l[i]"); exit(3);
	}
    for(j=0;j<MATRIX_SIZE*(local_chi+local_chi_3)/2;j++) uc_l[i][j]=0.;
     if(uc_l[i] == 0)
       PointerErr(cname,fname, "uc_l[i]");
    uc_nl[i] = (Float*)Alloc( MATRIX_SIZE * ((non_local_chi +non_local_chi_3[3])/2) * sizeof(Float) );
    if(uc_nl[i] == 0){
       PointerErr(cname,fname, "uc_nl[i]");
     }

  }


#define USE_SMEAR
//  lat_pt->Smear();

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
//Matrix * Fat = lat_pt->Fields(0);
matrix * Fat = fat;
matrix *nl[2]; 
matrix *l[2]; 
for(i = 0;i<2;i++){
nl[i] = (matrix *)(uc_nl[i]);
l[i] = (matrix *)(uc_l[i]);
}
  for ( n = 0; n < NUM_DIR/2; n++ ) {
    int coor=0;
    for (coord[3] = 0; coord[3] < size[3]; coord[3]++)
    for (coord[2] = 0; coord[2] < size[2]; coord[2]++)
    for (coord[1] = 0; coord[1] < size[1]; coord[1]++)
    for (coord[0] = 0; coord[0] < size[0]; coord[0]++){
//   	 for (i = 0; i < 4 ; i++) coord[i] = x[i];
	    odd = ( coord[0] + coord[1] + coord[2] + coord[3] ) % 2;
      if ( CoordNN( n ) ) {		// chi(x+mu) off-node
        *(nl[odd]) = Fat[LexGauge(coord)]; nl[odd]++;
        non_local_count[odd]++;
      }
      else {
        *(l[odd]) = Fat[LexGauge(coord)]; l[odd]++;
        local_count[odd]++;
      }
      coor++;
    }
    Fat += vol;
  }
//  Fat = lat_pt->Fields(0);
    Fat = fat;
//  Float rcv_mat[18];
//  printf("sizeof(Matrix)=%d\n",sizeof(Matrix));
  Float *rcv_mat = (Float *)qalloc(QFAST|QNONCACHE,sizeof(matrix));
  SCUDir snd_dirs[]={SCU_TP,SCU_XP,SCU_YP,SCU_ZP};
  SCUDir rcv_dirs[]={SCU_TM,SCU_XM,SCU_YM,SCU_ZM};
  sys_cacheflush(0);
  for ( n = 0; n < NUM_DIR/2; n++ ) {
    SCUDirArgIR snd(Fat,snd_dirs[n],SCU_SEND,sizeof(matrix));
    SCUDirArgIR rcv(rcv_mat,rcv_dirs[n],SCU_REC,sizeof(matrix));
    for (x[3] = 0; x[3] < size[3]; x[3]++)
    for (x[2] = 0; x[2] < size[2]; x[2]++)
    for (x[1] = 0; x[1] < size[1]; x[1]++)
    for (x[0] = 0; x[0] < size[0]; x[0]++){
   	 for (i = 0; i < 4 ; i++) coord[i] = x[i];
	    odd = ( coord[0] + coord[1] + coord[2] + coord[3] ) % 2;
      if ( CoordNN( n+4 ) ) {		// chi(chi-mu) off-node
        snd.Addr(Fat+LexGauge(coord_nn));
		snd.StartTrans();rcv.StartTrans();
		snd.TransComplete();rcv.TransComplete();
        nl[odd] ->Dagger((const Float *)rcv_mat ); 
//	    *(nl[odd]) *= (Float)-1.;
	    (nl[odd])->Negate();
        nl[odd]++;
        non_local_count[odd]++;
      }
      else {
#if 1
        l[odd] ->Dagger((const Float *)(Fat+LexGauge(coord_nn)) ); 
//	    *l[odd] *= (Float)-1.; l[odd]++;
	    (l[odd])->Negate(); l[odd]++;
#else
        tmp = (matrix *)(uc_l[odd] + MATRIX_SIZE * local_count[odd]);
        tmp ->Dagger((const Float *)(Fat+LexGauge(coord_nn)) ); 
	    *tmp *= (Float)-1.;
#endif
        local_count[odd]++;
      }
    }
    Fat += vol;
  }

matrix * Naik = naik;
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
            tmp = (matrix *)(uc_nl[odd] + MATRIX_SIZE * (non_local_count_3[j][odd]+(non_local_chi_3[j]+non_local_chi)/2));
        *tmp = Naik[LexGauge(coord)]; 
            non_local_count_3[j][odd]++;
          }
      } else {
        tmp = (matrix *)(uc_l[odd] + MATRIX_SIZE * (local_chi/2+local_count_3[odd]));
        *tmp = Naik[LexGauge(coord)]; 
        local_count_3[odd]++;
      }
    }
    Naik += vol;
  }

  if (naik_m){
//  if (0){
  Naik = naik_m;
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
            tmp = (matrix *)(uc_nl[odd] + MATRIX_SIZE * (non_local_count_3[j][odd]+(non_local_chi_3[j]+non_local_chi)/2));
        *tmp = Naik[LexGauge(coord)]; 
            non_local_count_3[j][odd]++;
          }
      } else {
        tmp = (matrix *)(uc_l[odd] + MATRIX_SIZE * (local_chi/2+local_count_3[odd]));
        *tmp = Naik[LexGauge(coord)]; 
        local_count_3[odd]++;
      }
    }
    Naik += vol;
  }
  } else {
    printf("USING NAIK INSTEAD OF NAIK_M\n");
    for(i = 0;i<4;i++)
    if (size[i]<3) {
      printf("Asqd::size[%d](%d) <3\n",i,size[i]);
      exit(13);
    }
    Naik = naik;
    sys_cacheflush(0);
    for ( n = 0; n < NUM_DIR/2; n++ ) {
      SCUDirArgIR snd(Naik,snd_dirs[n],SCU_SEND,sizeof(matrix));
      SCUDirArgIR rcv(rcv_mat,rcv_dirs[n],SCU_REC,sizeof(matrix));
      for (x[3] = 0; x[3] < size[3]; x[3]++)
      for (x[2] = 0; x[2] < size[2]; x[2]++)
      for (x[1] = 0; x[1] < size[1]; x[1]++)
      for (x[0] = 0; x[0] < size[0]; x[0]++){
     	for (i = 0; i < 4 ; i++) coord[i] = x[i];
  	    odd = ( coord[0] + coord[1] + coord[2] + coord[3] ) % 2;
        if ( CoordkNN( n+4,3 ) ) {		// chi(chi-mu) off-node
          for(int j=0;j<3;j++)
          if(coord_knn[n%4]==(size[n%4]-1-j)){
            snd.Addr(Naik+LexGauge(coord_knn));
            snd.StartTrans();rcv.StartTrans();
  	  	  snd.TransComplete();rcv.TransComplete();
            tmp = (matrix *)(uc_nl[odd] + MATRIX_SIZE * (non_local_count_3[j][odd]+(non_local_chi_3[j]+non_local_chi)/2));
            tmp->Dagger((const Float *)rcv_mat ); 
  	        tmp->Negate();
            non_local_count_3[j][odd]++;
          }
        }
        else {
          tmp = (matrix *)(uc_l[odd] + MATRIX_SIZE * (local_chi/2+local_count_3[odd]));
          tmp->Dagger((const Float *)(Naik+LexGauge(coord_knn)) ); 
  	      tmp->Negate();
          local_count_3[odd]++;
        }
      }
      Naik += vol;
    }
  }
  qfree (rcv_mat);

}

#if 1
  tmpfrm = NULL;
  if (vol<=4096)
  tmpfrm = (Float *) qalloc (QFAST|QCOMMS,NUM_DIR*2 * vol/2 * VECT_LEN2 * sizeof(Float));
  if(tmpfrm == NULL){ 
    tmpfrm = (Float *) qalloc (QCOMMS,NUM_DIR*2 * vol/2 * VECT_LEN2 * sizeof(Float));
    printf("tmpfrm is allocated at (%p),length 0x%x \n",tmpfrm,NUM_DIR*vol*VECT_LEN2*sizeof(Float));
  }
  if(tmpfrm == 0) 
    PointerErr(cname,fname, "tmpfrm");
  //printf("tmpfrm=%p\n",tmpfrm);
#endif

  for(i=0;i<2;i++){
    if(vol> 4096) uc_l_agg[i]=NULL;
    else
      uc_l_agg[i]  = (gauge_agg *)qalloc(QFAST,((local_chi+local_chi_3)/2)*sizeof(gauge_agg));
    if(uc_l_agg[i] == NULL){
      uc_l_agg[i]  = (gauge_agg *)qalloc(QCOMMS,((local_chi+local_chi_3)/2)*sizeof(gauge_agg));
      printf("uc_l_agg[%d] is allocated at DDR (%p)\n",i,uc_l_agg[i]);
    }
    if(uc_l_agg[i] == 0)
      PointerErr(cname,fname, "uc_l_agg[i]");
    uc_nl_agg[i]  = (gauge_agg*)qalloc(QFAST,((non_local_chi+non_local_chi_3[3])/2)*sizeof(gauge_agg));
    if(uc_nl_agg[i] == 0){
      uc_nl_agg[i]  = (gauge_agg*)qalloc(QCOMMS,((non_local_chi+non_local_chi_3[3])/2)*sizeof(gauge_agg));
      printf("uc_nl_agg[%d] is allocated at DDR (%p)\n",i,uc_nl_agg[i]);
    }
    if(uc_nl_agg[i] == 0){
      PointerErr(cname,fname, "uc_nl_agg[i]");
    }
  }

  gauge_agg *temp = new gauge_agg[12*vol];
  int num_ind[vol*6];
  int src;
//  printf("chi_l=0x%x\nchi_l[0][0]=%d\n",chi_l,chi_l[0][0]);
  for(j=0;j<2;j++){
    for(i=0;i<vol;i++) num_ind[i]=0;
    for(i=0;i< ((local_chi+ local_chi_3)/2);i++){
      src = (int)chi_l[j][2*i];
//      printf("chi_l[%d][%d]=%d %d\n",j,i*2,chi_l[j][2*i],src);
      if (src%(VECT_LEN*sizeof(Float))!=0){
        printf("%s::%s: src = %d\n",cname,fname,src);
        exit(1);
      }
      src = src/(VECT_LEN*sizeof(Float));
//      printf("src[%d]=%d\n",i,src);
      if(src > vol/2) {
        printf("%s::%s: src[%d](%d) > vol/2\n",cname,fname,i,src);
//        ERR.General(cname,fname,"src[%d](%d) > vol/2\n",i,src);
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


  for(j=0;j<2;j++){
    for(i=0;i<non_local_chi*3/2;i++) num_ind[i]=0;
    for(i=0;i< ((non_local_chi+ non_local_chi_3[3])/2);i++){
      src = (int)chi_nl[j][2*i];
      if (src%(VECT_LEN*sizeof(Float))!=0){
        printf("%s::%s: src = %d\n",cname,fname,src);
//        ERR.General(cname,fname,"src = %d\n",src);
        exit(1);
      }
      src = src/(VECT_LEN*sizeof(Float));
      if(src > non_local_chi*3/2) {
        printf("%s::%s: src(%d) > non_local_chi*3\n",cname,fname,src);
//        ERR.General(cname,fname,"src(%d) > non_local_chi*3\n",src);
        exit(1);
      }
      temp[src*2+num_ind[src]].src = (int)chi_nl[j][2*i];
      temp[src*2+num_ind[src]].dest= (int)chi_nl[j][2*i+1];
      for(k=0;k<18;k++)
        temp[src*2+num_ind[src]].mat[k] = uc_nl[j][i*18+k];
      num_ind[src]++;
      if(num_ind[src]>2){
        printf("%s::%s: num_ind[%d](%d) > 2 \n",cname,fname,src, num_ind[src]);
//        ERR.General(cname,fname,"num_ind[%d](%d) > 2 \n",src, num_ind[src]);
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

  for ( i = 0; i < 2; i++){
  Free(uc_l[i]);
  Free(uc_nl[i]);
  }
  printf("init_g done\n");
  
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

  return c[0] + s[0] * ( c[1] + s[1] * ( c[2] + s[2] * c[3] ));
}

void AsqD::comm_assert()
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



//static unsigned long address[]={0x0,0x0};
void AsqD::dirac(Float* b, Float* a, int a_odd, int add_flag)
{
  int i;
  int odd = 1 - a_odd;
  long c = (long) chi_off_node_total;
  long uc_l = (long)uc_l_agg[odd];
  long uc_nl = (long)uc_nl_agg[odd];
  long uc_nl2 = (long)&(uc_nl_agg[odd][isplit]);
  

#undef PROFILE
#ifdef PROFILE
  int num_flops;
  Float dtime;  
  struct timeval start,end;
  
  num_flops = 0;
  gettimeofday(&start,NULL);
#endif
  if( (unsigned long)a != address[odd]){
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
#ifdef PROFILE
  gettimeofday(&end,NULL);
  print_flops(0,&start,&end);
#endif

  sys_cacheflush(0);

  //make sure spinor field is in main memory before starting transfers

  SCUmulti_1[odd].StartTrans();
  SCUmulti_1[odd].TransComplete();

#undef PROFILE
#ifdef PROFILE
  num_flops = 33*(local_chi + local_chi_3);
  gettimeofday(&start,NULL);
#endif

  //-----------------------------------------------------------------
  //do first local computations
  //-----------------------------------------------------------------
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


if(split) {

  SCUmulti_2[odd].StartTrans();
  SCUmulti_2[odd].TransComplete();


#undef PROFILE
#ifdef PROFILE
  num_flops = 66*isplit;
  gettimeofday(&start,NULL);
#endif

#ifdef CPP
  dirac_cmv_jcw_agg_cpp( isplit, (long)0, (long)&(uc_nl_agg[odd][0]), (long)c, (long)tmpfrm);
#else
  dirac_cmv_jcw_agg( isplit, (long)0, uc_nl, c, (long)tmpfrm);
#endif

#ifdef PROFILE
  gettimeofday(&end,NULL);
  print_flops(num_flops,&start,&end);
#endif

  SCUmulti[odd].StartTrans();
  SCUmulti[odd].TransComplete();

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
  SCUmulti[odd].TransComplete();

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
		Fprintf(stderr,"chi_off_node_total[%d]=%e\n",i,chi_off_node_total[i]);
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
    asqd_sum_acc_cpp( vol/2, (long)0, (long)tmpfrm, (long)b);
  }
#ifdef PROFILE
  dtime +=dclock();
  printf("dirac_cmv_jcw_agg::%ld flops/%e seconds = %e MFlops\n",num_flops,dtime,(Float)num_flops/(dtime*1e6));
#endif
//  DiracOp::CGflops +=573*vol;

//  VRB.FuncEnd("","asqtad_dirac");
}

//CPS_END_NAMESPACE
