/*! \file
  \brief  Definition of parallel transport definitions for QCDOC.
  
  $Id: pt.C,v 1.15 2005-01-13 07:46:20 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2005-01-13 07:46:20 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/parallel_transport/pt_base/qcdoc/pt.C,v 1.15 2005-01-13 07:46:20 chulwoo Exp $
//  $Id: pt.C,v 1.15 2005-01-13 07:46:20 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: pt.C,v $
//  $Revision: 1.15 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/parallel_transport/pt_base/qcdoc/pt.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
#include <config.h>
#include <util/pt_int.h>
#include <qalloc.h>
#if 0
#include <util/gjp.h>
#include <util/pt.h>
#include <util/time.h>
#include <sysfunc.h>
#include <comms/scu.h>
#include <stdio.h>
#endif

#undef CPP
//CPS_START_NAMESPACE
//void dirac_cmv_jcw_agg_cpp( int sites, long chi, long u,long in, long out);

//External function definitions
extern "C"{
  void cmm_agg(gauge_agg *chi, matrix *phi,matrix *result, int counter);
  void cmm_agg_cpp( int sites, long chi, long u,long in, long out);
  void cmv_agg_cpp( int sites, long u,long in, long out);
  void pt_asqtad_agg( int sites, long chi, long u,long in, long out);
  void copy_buffer(int n, long src, long dest, long ptable);
  // Assembler copying routines
  void copy_matrix(IFloat *res, IFloat *src, int *length, 
		   unsigned long *res_ptr, unsigned long *src_ptr);
  void copy_gauge(IFloat *res, struct gauge_agg *src, int *length,
		  unsigned long *res_ptr);
  // This is perhaps overkill but gives a couple of extra flops
  // cross_look - all input fields are lookup and sum to result
  void cross_look(IFloat *result, Float *fac, const IFloat *chi, const IFloat *phi,  
	     int counter, unsigned long *src, unsigned long *dest);
  // cross_lin - one input field is linear and sum to result
  void cross_lin(IFloat *result, Float *fac, const IFloat *chi, const IFloat *phi,  
	     int counter, unsigned long *dest, unsigned long *dest);
  // cross_over_look - all input fields are lookup and overwrite result
  void cross_over_look(IFloat *result, Float *fac, const IFloat *chi, const IFloat *phi,  
	     int counter, unsigned long *src, unsigned long *dest);
  // cross_over_lin - one input field is linear and overwrite result
  void cross_over_lin(IFloat *result, Float *fac, const IFloat *chi, const IFloat *phi,  
	     int counter, unsigned long *dest, unsigned long *dest);

  //---------------------------------------------------------------------------
  //matrix multiply for checkerboarded fields
  void cmm_agg_cpp_cb(int sites, long u, long in, long out, IFloat *gauge_field, int pad=0);

  //matrix vector multiply for checkerboarded fields
  void cmv_agg_cpp_cb(int sites, long u, long in, long out, IFloat * gauge_field, int pad=0);
  //---------------------------------------------------------------------------
}

#if 0
//Number of dimensions
static const int NDIM=4;
//Maximum length of parallel transport.  3 links for the Naik term
static const int MAX_HOP=3;
//Local volume size in each of four directions
static int size[NDIM];
//SU(3) matrix has 18 floating point numbers
//SU(3) vector has 6 floating point numbers
enum {GAUGE_LEN=18,VECT_LEN=6, VECT_LEN2=6};
int vol = 1;

//gauge_agg holds source and destination indices, as well as the SU(3)
//link matrix.  One for local parallel transport, another for non-local
static gauge_agg *uc_l[2*SCUMachDim];
static gauge_agg *uc_nl[2*SCUMachDim];

//---------------------------------------------------------------------------
//Holds source,destination indexes for the matrix multiplication.
//Also holds index for gauge field, and whether gauge field needs to be
//conjugated
//First array index = 0 for even parity block, =1 for odd parity block
//
//uc_nl_cb_pre holds information for pre-multiplication of the fields
//that are transported in the positive direction.

static gauge_agg_cb *uc_l_cb[2][2*SCUMachDim];
static gauge_agg_cb *uc_nl_cb[2][2*SCUMachDim];
static gauge_agg_cb *uc_nl_cb_pre[2][SCUMachDim];

//--------------------------------------------------------------------------

static hop_pointer *hp_l[MAX_HOP][2*SCUMachDim];
static hop_pointer *hp_nl[MAX_HOP][2*SCUMachDim];

static unsigned long *src_l[MAX_HOP][2*SCUMachDim];
static unsigned long *dest_l[MAX_HOP][2*SCUMachDim];
static unsigned long *src_nl[MAX_HOP][2*SCUMachDim];
static unsigned long *dest_nl[MAX_HOP][2*SCUMachDim];

//Length of block of data for SCU communication
static int blklen[2*SCUMachDim];

//Number of blocks of data
static int numblk[2*SCUMachDim];

//Stride between blocks
static int stride[2*SCUMachDim];

//number of parallel transports that can be done locally
static int local_chi[2*SCUMachDim];

//Parallel transports that require non-local communication
static int non_local_chi[2*SCUMachDim];

//Initial offset for the data when using SCU communication
static int offset[2*SCUMachDim];

//--------------------------------------------------------------
//Checkerboarded data
static int blklen_cb[2*SCUMachDim];
static int numblk_cb[2*SCUMachDim];
static int stride_cb[2*SCUMachDim];
static int local_chi_cb[2*SCUMachDim];
static int non_local_chi_cb[2*SCUMachDim];
static int offset_cb[2*SCUMachDim];

//This determines whether the gauge links that are stored
//are the normal, canonical gauge links U_mu(x)
//or the conjugated versions, U_mu(x).dagger()
//For the staggered storage, the links are conjugated.

DagType conjugated;
//--------------------------------------------------------------

//Buffer for receiving data via SCU
static IFloat *rcv_buf[2*6];
static IFloat *rcv_buf2[2*6];

//--------------------------------------------------------------
//Send buffer transfers in the positive direction
static IFloat *snd_buf_cb[6];

//Buffer for transfer in the negative T direction.  This is needed because
//the block-stride communication does not work for communication
//in the T direction on a checkerboard lattice
static IFloat *snd_buf_t_cb;

//List of indexes for the vectors that are transferred when
//communicating in the negative T direction
static int *Toffset[2];
//--------------------------------------------------------------

//Pointer to the gauge field
static IFloat *gauge_field_addr;

//SCU communication parameters
static SCUDirArgIR *SCUarg[MAX_HOP][4*SCUMachDim];
static SCUDirArgIR *SCUarg_mat[MAX_HOP][4*SCUMachDim];
static SCUDirArgIR *SCUarg2[MAX_HOP][4*SCUMachDim];

//--------------------------------------------------------------
//Checkerboarded SCU
static SCUDirArgIR *SCUarg_cb[4*SCUMachDim];
static SCUDirArgIR *SCUarg_mat_cb[4*SCUMachDim];
//--------------------------------------------------------------

//Function primitives
static void (*Copy) (IFloat *dest, IFloat *src);
static void (*DagCopy) (IFloat *dest, IFloat *src);
static int (*LexVector)(int *x);
static int (*LexGauge) (int *x,int mu);

//-------------------------------------------------------------
//Added for checkerboarded parallel transport
static int (*LexVector_cb)(int *x);
//------------------------------------------------------------
#endif

int PT::size[NDIM];
int PT::vol;
//dest=src
void PT::cpy (IFloat *dest, IFloat *src){
  for(int i=0;i<18;i++)
    dest[i]=src[i];
}

//dest=src.Dagger()
void PT::dag_cpy (IFloat *dest, IFloat *src){
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++){
      dest[2*(3*i+j)]=src[2*(3*j+i)];
      dest[2*(3*i+j)+1]=-src[2*(3*j+i)+1];
    }
}

//Returns lexical index associated with coordinate x[4]
//where the 0th coordinate runs fastest, 3rd coordinate runs slowest
int PT::lex_xyzt(int *x){
//  printf("lex_xyzt(%d %d %d %d)\n",x[0],x[1],x[2],x[3]);
  int result = x[0] + size[0]*(x[1]+size[1]*(x[2]+size[2]*x[3]));
  return result;
}

//Returns checkerboard index associated with coordinate x[4]
int PT::lex_xyzt_cb_o(int *x){
//  printf("lex_xyzt_cb_o(%d %d %d %d)\n",x[0],x[1],x[2],x[3]);
  int result = x[0] + size[0]*(x[1]+size[1]*(x[2]+size[2]*x[3]));
  if ( (x[0]+x[1]+x[2]+x[3])%2 == 0) result = result/2+vol/2;
  else result = result/2;
  return result;
}

//---------------------------------------------------------------------------
//Returns index for fields in the STAG storage order on lattice
//sites of a given parity
int PT:: lex_txyz_cb(int *x)
{
//  printf("lex_txyz_cb(%d %d %d %d)\n",x[0],x[1],x[2],x[3]);
  int result = x[3]+size[3]*(x[0]+size[0]*(x[1]+size[1]*x[2]));
  return result/2;
}
//---------------------------------------------------------------------------

//Returns index associated with x[4] for txyz ordering
int PT:: lex_txyz(int *x){
  return  (x[3] + size[3]*(x[0]+size[0]*(x[1]+size[1]*x[2])))/2 ;
}

//Returns first index associated with the surface x[3] = 0, on
//a checkerboarded lattice
int PT:: LexSurface(int *x){
  return  (x[0]+size[0]*(x[1]+size[1]*x[2]))/2 ;
}

//Returns index associated with gauge link in the mu direction 
//and coordinate x
int PT::lex_g_xyzt(int *x, int mu){
  int temp =  lex_xyzt(x);
  return (temp*NDIM + mu);
}

//Returns index associated with gauge link in the mu direction and 
//coordinate x for checkerboarded storage
int PT:: lex_g_xyzt_cb_o(int *x, int mu){
  int temp =  lex_xyzt_cb_o(x);
  return (temp*NDIM + mu);
}

// Calculate the required offset given the direction and hop
int PT::set_offset(int dir, int hop) {

  // if positive direction then start at 0
  if (dir%2 == 0) return 0;

  int temp=1;
  int offset=0;
  for(int i=0;i<dir/2+1;i++){
    offset = temp*(size[i]-hop);
    temp *= size[i];
  }
  return offset;

}

void PT::set_hop_pointer() {

  char *fname = "set_hop_pointer()";

//  VRB.Func("PT",fname);
  //Actual memory usage of vectors
  int vlen = VECT_LEN*sizeof(IFloat);
  int vlen2 =VECT_LEN2*sizeof(IFloat);

  int x[NDIM], nei[NDIM];
  
  //Counts how many parallel transports of given length and direction are local
  //and non-local, respectively
  int hp_local_count[MAX_HOP][2*NDIM];
  int hp_non_local_count[MAX_HOP][2*NDIM];
  int hop, i;

#if 0
  //Local volume size in four directions
  size[0] = GJP.XnodeSites();
  size[1] = GJP.YnodeSites();
  size[2] = GJP.ZnodeSites();
  size[3] = GJP.TnodeSites();
#endif

  //Initialize local and non-local hop counters.
      printf("start\n");
  for (hop=0; hop<MAX_HOP; hop++) {
    for (i=0; i<2*NDIM; i++) {
      hp_non_local_count[hop][i] = 0;
      hp_local_count[hop][i] = 0;
    }
  }
  
  //For a given length of the parallel transport
  for (hop = 1; hop <= MAX_HOP; hop++) {
    hop_pointer **h_l = hp_l[hop-1];
    hop_pointer **h_nl = hp_nl[hop-1];

    //Local and non-local counts for given length of the hop
    int *local_count = hp_local_count[hop-1];
    int *non_local_count = hp_non_local_count[hop-1];

    //Loop over all directions
    for (i=0; i<NDIM; i++) {

      //Total number of sites that require non-local communication
      int non_local_check = hop*non_local_chi[i*2];
      //Total number of sites where parallel transport can be done locally
      int local_check = vol - non_local_check;

      //Loop through all the sites on the lattice
      //nei represents the coordinates of the neighboring site.
      for(x[3]=0,nei[3]=0;x[3]<size[3];x[3]++,nei[3]++)
	for(x[2]=0,nei[2]=0;x[2]<size[2];x[2]++,nei[2]++)
	  for(x[1]=0,nei[1]=0;x[1]<size[1];x[1]++,nei[1]++)
	    for(x[0]=0,nei[0]=0;x[0]<size[0];x[0]++,nei[0]++){

	      //This is the parallel transport of the field in the 
	      //negative direction to another node
	      //"Positive hop" because the link variable points in the 
	      //positive direction, even though the resulting field is 
	      //"transported" in the negative direction
	      // positive direction

	      if(x[i] < hop){
		//This calculates the neighbor coordinate
		nei[i] = size[i]-hop+x[i];  

		//Sets the index for source and destination
		(h_nl[2*i]+non_local_count[2*i])->src = non_local_count[2*i]*vlen;
		(h_nl[2*i]+non_local_count[2*i])->dest = LexVector(nei)*vlen2;

		//Increments the non-local count
		non_local_count[i*2]++;

		//Make sure we haven't gone over the non non-local check
		if (non_local_count[i*2]>non_local_check)
		  fprintf(stderr,
			"%s:non_local_count[%d](%d)>non_local_check[%d](%d)\n",
			 fname,2*i,non_local_count[2*i],2*i,non_local_check);
		//The rest of the parallel transports in the local volume can 
		//be handled locally
	      } else {
		//Calculate the new coordinate
		nei[i] = x[i]-hop;

		//Calculate the index for the source and the destination
		(h_l[2*i]+local_count[2*i])->src = LexVector(x)*vlen;
		(h_l[2*i]+local_count[2*i])->dest = LexVector(nei)*vlen2;
		
		//Increment the local count
		local_count[i*2]++;
		//Make sure we haven't exceeded the number of local sites
		if (local_count[i*2]>local_check)
		  fprintf(stderr,"%s:local_count[%d](%d)>local_check[%d](%d)\n",
			      fname,2*i,local_count[2*i],2*i,local_check);
	      }
	      
	      //Consider hopping in the negative direction, which is parallel 
	      //transport in the positive direction
	      // negative direction
	      if(x[i] >= (size[i]-hop)){
		//Calculate the non-local coordinate for this hop
		nei[i] = x[i]+hop-size[i];
		//Calculate source and destination indices
		(h_nl[2*i+1]+non_local_count[2*i+1])->src = non_local_count[2*i+1]*vlen;
		(h_nl[2*i+1]+non_local_count[2*i+1])->dest = LexVector(nei)*vlen2;

		//Increment the non-local count, check that bounds have not 
		//been exceeded
		non_local_count[i*2+1]++;
		if (non_local_count[i*2]>non_local_check)
		  fprintf(stderr,"%s:non_local_count[%d](%d)>non_local_check[%d](%d)\n",
			      fname,2*i,non_local_count[2*i],2*i,non_local_check);
	      } else {
		//Calculate the local coordinate for this hop
		nei[i] = x[i]+hop;
		//Calculate source and destination indices
		(h_l[2*i+1]+local_count[2*i+1])->src = LexVector(x)*vlen;
		(h_l[2*i+1]+local_count[2*i+1])->dest = LexVector(nei)*vlen2;
		//Increment local count, check that bounds not exceeded
		local_count[i*2+1]++;
		if (local_count[i*2]>local_check)
		  fprintf(stderr,"%s:local_count[%d](%d)>local_check[%d](%d)\n",
			      fname,2*i,local_count[2*i],2*i,local_check);
	      }
	      // Need to reset the neighbour pointer
	      nei[i] = x[i];
	    }
    }
  }
//  VRB.Func("PT",fname);
//  exit(44);
}


//Initialization of Parallel Transport class
void PT::init(PTArg *pt_arg)
{
  char *cname = "";
  char *fname = "pt_init()";
//  VRB.Func("",fname);
  int i, j, x[NDIM],nei[NDIM];
  int local_count[2*NDIM];
  int non_local_count[2*NDIM];
  int vlen = VECT_LEN*sizeof(IFloat); //size of incoming vector
  int vlen2 =VECT_LEN2*sizeof(IFloat); //size of outgoing vector (maybe different to optimize for QCDOC PEC)

  //---------------------------------------------------------------------------
  int local_count_cb[2][2*NDIM];
  int non_local_count_cb[2][2*NDIM];
  //---------------------------------------------------------------------------

  size[0] = pt_arg->size[0];
  size[1] = pt_arg->size[1];
  size[2] = pt_arg->size[2];
  size[3] = pt_arg->size[3];
  gauge_field_addr = pt_arg->gauge_field_addr;
  g_str_ord = pt_arg->g_str_ord;
  g_conj = pt_arg->g_conj;
  v_str_ord = pt_arg->v_str_ord;
  v_str_ord_cb = pt_arg->v_str_ord_cb;
  evenodd = pt_arg->evenodd;
  prec = pt_arg->prec;

  switch(g_str_ord){
    case PT_XYZT:
      LexGauge = lex_g_xyzt;
      break;
    case PT_XYZT_CB_O:
      LexGauge = lex_g_xyzt_cb_o;
      break;
    default:
      fprintf(stderr,"PT::init got invalid g_str_ord\n");
      break;
  }

  switch(v_str_ord){
    case PT_XYZT:
      LexVector = lex_xyzt;
      break;
    case PT_XYZT_CB_O:
      LexVector = lex_xyzt_cb_o;
      break;
    default:
      fprintf(stderr,"PT::init got invalid v_str_ord\n");
      break;
  }

  switch(v_str_ord_cb){
    case PT_TXYZ:
      LexVector_cb = lex_txyz_cb;
      break;
    default:
      fprintf(stderr,"PT::init got invalid v_str_ord_cb\n");
      break;
  }
  
  if (g_conj) { Copy = &dag_cpy; DagCopy = &cpy; conjugated = PT_DAG_YES;}
  else        { Copy = &cpy; DagCopy = &dag_cpy; conjugated = PT_DAG_NO;}

#if 0
  //Size of local volume in all four directions
  size[0] = GJP.XnodeSites();
  size[1] = GJP.YnodeSites();
  size[2] = GJP.ZnodeSites();
  size[3] = GJP.TnodeSites();

  //Location of the gauge field
  gauge_field_addr = (IFloat *) lat.GaugeField();

  //Retrive the storage order
  StrOrdType str_ord = lat.StrOrd();

  //Define the appropriate function to calculate the indices for the fermion
  //fields and the guage links.

  //Canonical ordering xyzt without checkerboarding
  if (str_ord == PT_CANONICAL){
    Copy = cpy; DagCopy = dag_cpy;
    LexVector = lex_xyzt;
    LexGauge = lex_g_xyzt;

    //------------------------------------------------------------------------
    LexVector_cb = lex_txyz_cb;
    conjugated = DAG_NO;
    //------------------------------------------------------------------------

    //Wilson ordering is checkerboarded
  } else if (str_ord == WILSON){
    Copy = cpy; DagCopy = dag_cpy;
    LexVector = lex_xyzt_cb_o;
    LexGauge = lex_g_xyzt_cb_o;

    //-----------------------------------------------------------------------
    LexVector_cb = lex_txyz_cb;
    conjugated = DAG_NO;
    //-----------------------------------------------------------------------

    //Staggered ordering
  } else if (str_ord == STAG){
    Copy = dag_cpy; DagCopy = cpy;
    LexVector = lex_xyzt;
    LexGauge = lex_g_xyzt;

    //-----------------------------------------------------------------------
    LexVector_cb = lex_txyz_cb;
    conjugated = DAG_YES;
    //----------------------------------------------------------------------

  } else
    fprintf(stderr,"storage ordering not implemented");
#endif

  //For the fastest changing index, data must be sent in many short messages
  //For the slowest changing index, the boundaries of the hypersurface are 
  //stored together in large blocks, so a few long messages can be sent.

  blklen[0] = blklen[1]= vlen;
  for(i=1;i<NDIM;i++) {blklen[2*i+1] = blklen[2*i] = blklen[2*i-1]*size[i-1]; }

  numblk[2*NDIM-1]=numblk[2*NDIM-2]=1;
  for(i=NDIM-2;i>=0;i--) {numblk[i*2+1] = numblk[2*i] = numblk[2*i+2]*size[i+1]; }

  //The stride length is longer when the blocks are large
  for(i=0;i<NDIM*2;i++)  stride[i] = blklen[i]* (size[i/2]-1);

  //Calculate the local volume
  vol = 1;
  for(i=0; i<NDIM;i++) vol *= size[i];

  //Calculate the number of local and non-local parallel transports
  //are needed in each direction
  for(i=0; i<NDIM;i++){
    non_local_chi[2*i+1] = non_local_chi[2*i] = vol/size[i];
    local_chi[2*i+1] = local_chi[2*i] = vol - non_local_chi[2*i];
  }

  //---------------------------------------------------------------------------
  //Calculation of block length, number of blocks, and stride for 
  //checkerboarded storage

  //Block length for checkerboarded scheme is similar to canonical scheme
  //However, there must be a few modifications.
  //
  //For the T (3rd) direction, a block strided communciation will not work
  //As a result, pointers to the appropriate fields will be aggregated
  //into a send buffer before being sent as a single block.
  //
  //For transfers in all other directions, a block-strided move is allowed.
  //In these cases, the fastest changing index (X) requires many short messages
  //while the slowest changing index (Z) can be transfered in one block
  
  blklen_cb[6] = blklen_cb[7] = vol*vlen/(2*size[3]);
  numblk_cb[6] = numblk_cb[7] = 1;

  blklen_cb[4] = blklen_cb[5] = vol*vlen/(2*size[2]);
  numblk_cb[4] = numblk_cb[5] = 1;

  blklen_cb[2] = blklen_cb[3] = vol*vlen/(2*size[1]*size[2]);
  numblk_cb[2] = blklen_cb[3] = size[2];

  blklen_cb[0] = blklen_cb[1] = vol*vlen/(2*size[0]*size[1]*size[2]);
  numblk_cb[0] = numblk_cb[1] = size[1]*size[2];

  //The stride is also similar
  stride_cb[0] = stride_cb[1] = 0;
  for(i = 0; i < NDIM*2; i++)
      stride_cb[i] = blklen_cb[i] * (size[i/2]-1);
  stride_cb[6] = stride_cb[7] = 0;

  //Calculate the number of local and non-local parallel transports
  for(i = 0; i < NDIM; i++)
    {
      non_local_chi_cb[2*i+1] = non_local_chi_cb[2*i] = vol/(2*size[i]);
      local_chi_cb[2*i+1] = local_chi_cb[2*i] = vol/2 - non_local_chi_cb[2*i];
    }

  //---------------------------------------------------------------------------

  for(i=0; i<2*NDIM;i++){
    local_count[i]=non_local_count[i]=0;
    //-------------------------------------------------------------------------
    for(int parity = 0; parity < 2; parity++)
      local_count_cb[parity][i] = non_local_count_cb[parity][i] = 0;
    //------------------------------------------------------------------------

    //For small volumes, we can allocate the memory for the gauge aggregate
    //in the faster part of memory
    if(vol> 4096){
      uc_l[i]=uc_nl[i]=NULL;
    } else {	
      uc_l[i] = (gauge_agg *)qalloc(QFAST,sizeof(gauge_agg)*(1+local_chi[i]));
      uc_nl[i] = (gauge_agg *)qalloc(QFAST,sizeof(gauge_agg)*(1+non_local_chi[i]));
    }

    //Otherwise, we will have to allocate using smalloc
    if(uc_l[i]==NULL) {
      uc_l[i] = (gauge_agg *)Alloc(cname,fname,"uc_l[i]",sizeof(gauge_agg)*(1+local_chi[i]));
    }
    if(uc_nl[i]==NULL) {
      uc_nl[i] = (gauge_agg *)Alloc(cname,fname,"uc_nl[i]",sizeof(gauge_agg)*(1+non_local_chi[i]));
    }

    //-------------------------------------------------------------------------
    //Allocate memory for gauge_agg_cb
    for(int parity = 0; parity < 2; parity++)
    {
      uc_l_cb[parity][i] = (gauge_agg_cb *)Alloc(cname,fname,"uc_l_cb[parity][i]",sizeof(gauge_agg_cb)*(1+local_chi_cb[i]));
      uc_nl_cb[parity][i] = (gauge_agg_cb *)Alloc(cname,fname,"uc_nl_cb[parity][i]",sizeof(gauge_agg_cb)*(1+non_local_chi_cb[i]));
    }
    //-------------------------------------------------------------------------

    // This buffer is actually overkill, but ensures will work if
    // shift_field is called with hop>1
    rcv_buf[i] = (IFloat *)qalloc(QCOMMS,3*MAX_HOP*non_local_chi[i]*vlen);
    if(rcv_buf[i]==NULL)PointerErr("",fname,"rcv_buf[i]");

    //Used buffer used in vvpd
    rcv_buf2[i] = (IFloat *)qalloc(QCOMMS,MAX_HOP*non_local_chi[i]*vlen);
    if(rcv_buf2[i]==NULL)PointerErr("",fname,"rcv_buf2[i]");
  }

  //---------------------------------------------------------------------------
  //Allocate memory for send buffer
  for(i=0; i<NDIM;i++)
    {
      snd_buf_cb[i] = (IFloat *)qalloc(QCOMMS,3*non_local_chi_cb[2*i+1]*vlen);
      if(snd_buf_cb[i]==NULL)PointerErr("",fname,"snd_buf_cb[i]");
    }
  snd_buf_t_cb = (IFloat *)qalloc(QCOMMS,3*non_local_chi_cb[6]*vlen);
  if(snd_buf_t_cb==NULL) PointerErr("",fname,"snd_buf_t_cb");

  for(i = 0; i < 2;i++)
    Toffset[i] = (int *)Alloc(cname,fname,"Toffset[parity]",non_local_chi_cb[6]*sizeof(int));

  //Allocate memory for the gauge_agg_cb used for matrix pre-multiplication
  for(i = 0; i< NDIM;i++)
    for(int parity = 0; parity<2;parity++)
      uc_nl_cb_pre[parity][i] = (gauge_agg_cb *)Alloc(cname,fname,"uc_nl_cb_pre[parity][i]",sizeof(gauge_agg_cb)*(1+non_local_chi_cb[2*i+1]));

  int parity = 0;
  //---------------------------------------------------------------------------
 
  //Calculate source and destination indices for gauge aggregates (see set_hop_pointer())
  //Only for one hop
   for(x[3]=0,nei[3]=0;x[3]<size[3];x[3]++,nei[3]++)
    for(x[2]=0,nei[2]=0;x[2]<size[2];x[2]++,nei[2]++)
      for(x[1]=0,nei[1]=0;x[1]<size[1];x[1]++,nei[1]++)
	for(x[0]=0,nei[0]=0;x[0]<size[0];x[0]++,nei[0]++){
	  for(i=0;i<NDIM;i++){
	    
	    //printf("%d %d %d %d %d\n",x[0],x[1],x[2],x[3],i);
	    // positive direction
	    //This is for transport of a vector in the negative direction
	    //An even index for uc_nl, uc_l, uc_nl_cb, uc_l_cb corresponds
	    //to parallel transport in the negative direction
	    
	    if(x[i] == 0){
	      nei[i] = size[i]-1;
	      (uc_nl[2*i]+non_local_count[2*i])->src = non_local_count[2*i]*vlen;
	      (uc_nl[2*i]+non_local_count[2*i])->dest = LexVector(nei)*vlen2;
	      non_local_count[i*2]++;
	      if (non_local_count[i*2]>non_local_chi[i*2])
		fprintf(stderr,"%s:non_local_count[%d](%d)>non_local_chi[%d](%d)\n",
			fname,2*i,non_local_count[2*i],2*i,non_local_chi[2*i]);
	    } 
	    else {
	      nei[i] = x[i]-1;
	      if(local_count[2*i]<0) fprintf(stderr,"%s:local_count[%d]=%d]n",
			fname,2*i,local_count[2*i]);
	      (uc_l[2*i]+local_count[2*i])->src = LexVector(x)*vlen;
	      (uc_l[2*i]+local_count[2*i])->dest = LexVector(nei)*vlen2;
	      local_count[i*2]++;
	      if (local_count[i*2]>local_chi[i*2])
		fprintf(stderr,"%s:local_count[%d](%d)>local_chi[%d](%d)\n",
			fname,2*i,local_count[2*i],2*i,local_chi[2*i]);
	    }
	    // negative direction
	    //This is parallel transport in the positive direction
	    //An odd index for uc_l, uc_nl, uc_l_cb,uc_nl_cb corresponds to
	    //transport in the positive direction
	    if(x[i] == (size[i] -1)){
	      nei[i] = 0;
	      (uc_nl[2*i+1]+non_local_count[2*i+1])->src = non_local_count[2*i+1]*vlen;
	      (uc_nl[2*i+1]+non_local_count[2*i+1])->dest = LexVector(nei)*vlen2;
	      non_local_count[i*2+1]++;
	      if (non_local_count[i*2+1]>non_local_chi[i*2+1])
		fprintf(stderr,"%s:non_local_count[%d](%d)>non_local_chi[%d](%d)\n",
			fname,2*i+1,non_local_count[2*i+1],2*i+1,non_local_chi[2*i+1]);
	    } 
	    else {
	      nei[i] = x[i]+1;
	      if(local_count[2*i+1]<0) fprintf(stderr,"%s:local_count[%d]=%d]n",
			fname,2*i+local_count[2*i+1]);
	      (uc_l[2*i+1]+local_count[2*i+1])->src = LexVector(x)*vlen;
	      (uc_l[2*i+1]+local_count[2*i+1])->dest = LexVector(nei)*vlen2;
	      local_count[i*2+1]++;
	      if (local_count[i*2+1]>local_chi[i*2+1])
		fprintf(stderr,"%s:local_count[%d](%d)>local_chi[%d](%d)\n",
			fname,2*i+1,local_count[2*i+1],2*i+1,local_chi[2*i+1]);
	    }
	    nei[i] = x[i];
	  }
	}

   //--------------------------------------------------------------------------
  //Calculate source and destination indices for gauge aggregates
  //Only for one hop
   for(x[2]=0,nei[2]=0;x[2]<size[2];x[2]++,nei[2]++)
     for(x[1]=0,nei[1]=0;x[1]<size[1];x[1]++,nei[1]++)
       for(x[0]=0,nei[0]=0;x[0]<size[0];x[0]++,nei[0]++)
	 for(x[3]=0,nei[3]=0;x[3]<size[3];x[3]++,nei[3]++)
	   {

	    parity = (x[0]+x[1]+x[2]+x[3])%2;
	    //Calculate offsets for transfers in the negative T direction
	    if(x[3] == 0)
	      *(Toffset[parity] +non_local_count_cb[parity][6]) = LexVector_cb(x);

	    for(i=0;i<NDIM;i++){

	    //printf("%d %d %d %d %d\n",x[0],x[1],x[2],x[3],i);
	    // positive direction
	    //This is for transport of a vector in the negative direction
	    //An even index for uc_nl, uc_l, uc_nl_cb, uc_l_cb corresponds
	    //to parallel transport in the negative direction

	    if(x[i] == 0)
	      {
	      nei[i] = size[i]-1;

	      //The src and dest indexes index the Vector, and do not include 
	      //information for the size of the vector, nor
	      //the size of the IFloat object.  This is to allow runt-time 
	      //adjustment of these parameters
	      //
	      //src - Source index in the receive buffer coming from a 
	      //      positive adjacent node
	      //dest - Destination index for the transported vector field
	      //dest2 - Destination index for the padded vector field
	      //gauge_index - Index of the SU(3) gauge link needed to 
	      //              transport the field
	      //dagger - determines if the gauge link needs to be conjugated

	      (uc_nl_cb[parity][2*i]+non_local_count_cb[parity][2*i])->src = non_local_count_cb[parity][2*i];
	      (uc_nl_cb[parity][2*i]+non_local_count_cb[parity][2*i])->dest = LexVector_cb(nei);
	      (uc_nl_cb[parity][2*i]+non_local_count_cb[parity][2*i])->dest2 = LexVector_cb(nei)*8+2*i;
	      (uc_nl_cb[parity][2*i]+non_local_count_cb[parity][2*i])->gauge_index = LexGauge(nei,i)*GAUGE_LEN;
	      (uc_nl_cb[parity][2*i]+non_local_count_cb[parity][2*i])->dagger = (int) conjugated;

	      non_local_count_cb[parity][i*2]++;
	      if(non_local_count_cb[parity][i*2]>non_local_chi_cb[i*2])
		fprintf(stderr,
			"%s:non_local_count_cb[%d][%d](%d)>non_local_chi_cb[%d](%d)\n",
			fname,parity,2*i,non_local_count_cb[parity][2*i],2*i,non_local_chi[2*i]);
	    } 
	    else 
	      {
	      nei[i] = x[i]-1;

	      (uc_l_cb[parity][2*i]+local_count_cb[parity][2*i])->src = LexVector_cb(x);
	      (uc_l_cb[parity][2*i]+local_count_cb[parity][2*i])->dest = LexVector_cb(nei);
	      (uc_l_cb[parity][2*i]+local_count_cb[parity][2*i])->dest2 = LexVector_cb(nei)*8+2*i;
	      (uc_l_cb[parity][2*i]+local_count_cb[parity][2*i])->gauge_index = LexGauge(nei,i)*GAUGE_LEN;
	      (uc_l_cb[parity][2*i]+local_count_cb[parity][2*i])->dagger = (int) conjugated;
	      local_count_cb[parity][i*2]++;
	      if(local_count_cb[parity][i*2]>local_chi_cb[i*2])
		fprintf(stderr,"%s:local_count_cb[%d][%d](%d)>local_chi_cb[%d](%d)\n",parity,2*i,local_count_cb[parity][2*i],2*i,local_chi[2*i]);
	      }

	    // negative direction
	    //This is parallel transport in the positive direction
	    //An odd index for uc_l, uc_nl, uc_l_cb,uc_nl_cb corresponds to
	    //transport in the positive direction

	    if(x[i] == (size[i] -1))
	      {
	      nei[i] = 0;

	      //src - Source index in the receive buffer
	      //dest - Destination index for the transported vector field
	      //In only this case, the field is transported pre-multiplied by
	      //the gauge link.  As a result, gauge_index and dagger are not 
	      //strictly necessary.
	      //
	      //However, we do need to specify another gauge aggregate that 
	      //will contain the information necessary
	      //for the pre-multiplication of the SU(3) link matrix
	      
	      (uc_nl_cb[parity][2*i+1]+non_local_count_cb[parity][2*i+1])->src = non_local_count_cb[parity][2*i+1];
	      (uc_nl_cb[parity][2*i+1]+non_local_count_cb[parity][2*i+1])->dest = LexVector_cb(nei);
	      (uc_nl_cb[parity][2*i+1]+non_local_count_cb[parity][2*i+1])->dest2 = LexVector_cb(nei)*8+2*i+1;
	      (uc_nl_cb[parity][2*i+1]+non_local_count_cb[parity][2*i+1])->gauge_index = LexGauge(x,i)*GAUGE_LEN;
	      (uc_nl_cb[parity][2*i+1]+non_local_count_cb[parity][2*i+1])->dagger = 1-(int)conjugated;

	      (uc_nl_cb_pre[parity][i]+non_local_count_cb[parity][2*i+1])->src = LexVector_cb(x);
	      (uc_nl_cb_pre[parity][i]+non_local_count_cb[parity][2*i+1])->dest = non_local_count_cb[parity][2*i+1];
	      (uc_nl_cb_pre[parity][i]+non_local_count_cb[parity][2*i+1])->dest2 = LexVector_cb(nei)*8+2*i+1;
	      (uc_nl_cb_pre[parity][i]+non_local_count_cb[parity][2*i+1])->gauge_index = LexGauge(x,i)*GAUGE_LEN;
	      (uc_nl_cb_pre[parity][i]+non_local_count_cb[parity][2*i+1])->dagger = 1-(int)conjugated;

	      non_local_count_cb[parity][i*2+1]++;
	      if(non_local_count_cb[parity][i*2+1]>non_local_chi_cb[i*2+1])
		fprintf(stderr,"%s:non_local_count_cb[%d][%d](%d)>non_local_chi_cb[%d](%d)\n",parity,2*i+1,non_local_count_cb[parity][2*i+1],2*i+1,non_local_chi[2*i+1]);
	      } 
	    else 
	      {
	      nei[i] = x[i]+1;

	      (uc_l_cb[parity][2*i+1]+local_count_cb[parity][2*i+1])->src = LexVector_cb(x);
	      (uc_l_cb[parity][2*i+1]+local_count_cb[parity][2*i+1])->dest = LexVector_cb(nei);
	      (uc_l_cb[parity][2*i+1]+local_count_cb[parity][2*i+1])->dest2 = LexVector_cb(nei)*8+2*i+1;
	      (uc_l_cb[parity][2*i+1]+local_count_cb[parity][2*i+1])->gauge_index = LexGauge(x,i)*GAUGE_LEN;
	      (uc_l_cb[parity][2*i+1]+local_count_cb[parity][2*i+1])->dagger = 1-(int)conjugated;

	      local_count_cb[parity][i*2+1]++;
	      if(local_count_cb[parity][i*2+1]>local_chi_cb[i*2+1])
		fprintf(stderr,"%s:local_count_cb[%d][%d](%d)>local_chi_cb[%d](%d)\n",parity,2*i+1,local_count_cb[parity][2*i+1],2*i+1,local_chi[2*i+1]);
	    }
	    nei[i] = x[i];
	  }
	}
   //--------------------------------------------------------------------------

  //Sets bits in uc_l and uc_nl to zero
  for(i=0;i<NDIM;i++){
    memset(uc_l[2*i]+local_count[2*i],0,sizeof(gauge_agg));
    memset(uc_l[2*i+1]+local_count[2*i+1],0,sizeof(gauge_agg));
    memset(uc_nl[2*i]+non_local_count[2*i],0,sizeof(gauge_agg));
    memset(uc_nl[2*i+1]+non_local_count[2*i+1],0,sizeof(gauge_agg));
  }

  //Calculate offsets?
  //For even array index (transfer in the negative  direction) the
  //offset is 0
  //For odd indexes, the offsets are:
  //offset[1] = size[0]-1
  //offset[3] = size[0]*(size[1]-1)
  //offset[5] = size[0]*size[1]*(size[2]-1)
  //offset[7] = size[0]*size[1]*size[2]*(size[3]-1)
  //These offsets correspond to the starting index for data transfer
  //in the positive direction
  int temp=1;
  for(i=0;i<NDIM;i++){
    offset[2*i]  = 0;
    offset[2*i+1] = temp*(size[i]-1);
    temp *= size[i];
  }

  //-------------------------------------------------------------------
  temp = 1;
  for(i = 0;i<NDIM;i++)
    {
      offset_cb[2*i] = 0;
      offset_cb[2*i+1] = stride_cb[2*i+1]/sizeof(IFloat);
    }
  //-------------------------------------------------------------------

  // Allocate memory for hop pointers
  for (j=0; j<MAX_HOP; j++) {
  
    for(i=0; i<2*NDIM; i++){
      
      //Calculate the number of local and non-local sites needed 
      //j+1 is the length of the hop
      //i indicates the communication direction
      int nl_size = (j+1)*non_local_chi[i];
      int l_size = vol - nl_size;
      if (l_size>0){
        hp_l[j][i] = (hop_pointer*) Alloc(cname,fname,"hp_l[j][i]",l_size*sizeof(hop_pointer));
        src_l[j][i] = (unsigned long*)qalloc(0,l_size*sizeof(unsigned long));
        dest_l[j][i] = (unsigned long*)qalloc(0,l_size*sizeof(unsigned long));
      }
=======
      
      hp_nl[j][i] = (hop_pointer*) Alloc(cname,fname,"hp_nl[j][i]",nl_size*sizeof(hop_pointer));
      dest_nl[j][i] = (unsigned long*)qalloc(0,nl_size*sizeof(unsigned long));
      src_nl[j][i] = (unsigned long*)qalloc(0,nl_size*sizeof(unsigned long));
    }
  }
  
  set_hop_pointer();

  //Calculate the indices for the source and destination
  for (j=0; j<MAX_HOP; j++) {
    for(i=0; i<2*NDIM; i++){
      int nl_size = (j+1)*non_local_chi[i]+1;
      int l_size = vol - nl_size+2;
      if (l_size>2)
      for (int s=0; s<l_size; s++) {
	src_l[j][i][s] = hp_l[j][i][s].src/(VECT_LEN*sizeof(IFloat));
	dest_l[j][i][s] = hp_l[j][i][s].dest/(VECT_LEN2*sizeof(IFloat));
      }
      for (int s=0; s<nl_size; s++) {
	src_nl[j][i][s] = hp_nl[j][i][s].src/(VECT_LEN*sizeof(IFloat));
	dest_nl[j][i][s] = hp_nl[j][i][s].dest/(VECT_LEN2*sizeof(IFloat));
      }
    }
  }

	
}

//Free memory associated with the parallel transport of the fermions
void PT::delete_buf(){
  char *fname = "pt_delete()";
//  VRB.Func("",fname);
	
  for(int i = 0; i < 2*NDIM; i++){
    qfree(uc_l[i]);
    qfree(uc_nl[i]);
    //--------------------------------------------------------------------
    for(int parity = 0; parity < 2; parity++)
      {
	Free(uc_l_cb[parity][i]);
	Free(uc_nl_cb[parity][i]);
      }
    //-------------------------------------------------------------------
    qfree(rcv_buf[i]);
    qfree(rcv_buf2[i]);
  }

  //-----------------------------------------------------------------------
  for(int i = 0; i < NDIM; i++)
    {
      qfree(snd_buf_cb[i]);
      for(int parity = 0; parity < 2; parity++)
	Free(uc_nl_cb_pre[parity][i]);
    }
  Free(snd_buf_t_cb);
  for(int i = 0; i < 2; i++)
    Free(Toffset[i]);    
  //-----------------------------------------------------------------------

  for (int hop=0; hop<MAX_HOP; hop++) {
    for(int i = 0; i < 2*NDIM; i++){
      int nl_size = (hop+1)*non_local_chi[i];
      int l_size = vol - nl_size;
      if (l_size>0){
        Free(hp_l[hop][i]);
        qfree(src_l[hop][i]);
        qfree(dest_l[hop][i]);
      }
      Free(hp_nl[hop][i]);
      qfree(src_nl[hop][i]);
      qfree(dest_nl[hop][i]);
    }
  }
}

//Free memory associated with gauge parallel transport
void PT::delete_g_buf(){
  char *fname = "pt_delete_g()";
//  VRB.Func("",fname);
  for(int hop = 0; hop < MAX_HOP; hop++)
    for(int i = 0; i < 4*NDIM; i++)
      {
	delete SCUarg[hop][i];
	delete SCUarg2[hop][i];
	delete SCUarg_mat[hop][i];
      }

  //----------------------------------------------------------------------
  //Checkerboarding
  for(int i = 0; i < 4*NDIM; i++)
    {
      delete SCUarg_cb[i];
      delete SCUarg_mat_cb[i];
    }
  //---------------------------------------------------------------------
}

void PT::init_g(void){
  int x[NDIM], nei[NDIM];
  int local_count[2*NDIM];
  int non_local_count[2*NDIM];
  int i;

  char *fname = "init_g()";
//  VRB.Func("",fname);
  for(i=0; i<2*NDIM;i++){
    local_count[i]=non_local_count[i]=0;
  }
  //Location of gauge field
  IFloat *u = gauge_field_addr;

  //Send and receive directions
  SCUDir rcv_dir[]={SCU_XP, SCU_XM, SCU_YP, SCU_YM, SCU_ZP, SCU_ZM,SCU_TP,SCU_TM};
  SCUDir snd_dir[]={SCU_XM, SCU_XP, SCU_YM, SCU_YP, SCU_ZM, SCU_ZP,SCU_TM,SCU_TP};

  //Temporary buffer (allocated on cache) that receives an SU(3) matrix
  IFloat *rcv_mat = (IFloat *)qalloc(QFAST|QNONCACHE,18*sizeof(IFloat));
  sys_cacheflush(0);


  for(i=0;i<NDIM;i++){
    //Initialize SCUDirArg for sending and receiving the one-hop term
    SCUDirArgIR snd(u,snd_dir[i*2+1],SCU_SEND,sizeof(matrix));
    SCUDirArgIR rcv(rcv_mat,rcv_dir[i*2+1],SCU_REC,sizeof(matrix));

    for(x[3]=0,nei[3]=0;x[3]<size[3];x[3]++,nei[3]++)
      for(x[2]=0,nei[2]=0;x[2]<size[2];x[2]++,nei[2]++)
	for(x[1]=0,nei[1]=0;x[1]<size[1];x[1]++,nei[1]++)
	  for(x[0]=0,nei[0]=0;x[0]<size[0];x[0]++,nei[0]++){
	    // positive direction
	    //this is a hop in the positive direction, meaning data must be sent
	    //in the negative direction.
	    if(x[i] == 0){
	      //Calculate the appropriate coordinate on the adjacent node
	      nei[i] = size[i]-1;  
	      //Copy the appropriate matrix from u to uc_nl
	      Copy((uc_nl[2*i]+non_local_count[2*i])->mat, u+LexGauge(nei,i)*GAUGE_LEN);
	      non_local_count[i*2]++;
	    } else {
	      //Calculate the appropriate neighbor coordinate on the local node
	      nei[i] = x[i]-1;
	      //Copy from u to uc_l
	      Copy((uc_l[2*i]+local_count[2*i])->mat, u+LexGauge(nei,i)*GAUGE_LEN);
	      local_count[i*2]++;
	    }
	    // negative direction
	    if(x[i] == (size[i] -1)){
	      nei[i] = 0;
#if 0
	      getMinusData( rcv_mat, u+LexGauge(x,i)*GAUGE_LEN, GAUGE_LEN, i);
#else
	      //Send the appropriate matrix
	      snd.Addr(u+LexGauge(x,i)*GAUGE_LEN);
	      //Send the transmission, prepare to receive
	      snd.StartTrans();rcv.StartTrans();
	      //Complete the send and receive
	      snd.TransComplete();rcv.TransComplete();
#endif
	      //Copy to uc_nl from the received matrix
	      DagCopy((uc_nl[2*i+1]+non_local_count[2*i+1])->mat, rcv_mat);
	      non_local_count[i*2+1]++;
	    } else {
	      //Calculate the appropriate neighbor coordinate on the local volume
	      nei[i] = x[i]+1;
	      //Copy from u to uc_l
	      DagCopy((uc_l[2*i+1]+local_count[2*i+1])->mat, u+LexGauge(x,i)*GAUGE_LEN);
	      local_count[i*2+1]++;
	    }
	    nei[i] = x[i];
	  } // x[]
  } // i

#if 0
  for(i=0;i<2*NDIM;i++) {
    SCUarg[i*2] = new SCUDirArgIR;
    SCUarg[i*2]->Init((void *)rcv_buf[i],rcv_dir[i],SCU_REC,non_local_chi[i]*VECT_LEN*sizeof(IFloat),1,0,IR_9);
    SCUarg[i*2+1] = new SCUDirArgIR;
    //
    //  inputs a dummy but valid address to pass syscall test and changed later, CJ
    //
    SCUarg[i*2+1]->Init((void *)rcv_buf[i],snd_dir[i],SCU_SEND,blklen[i],numblk[i],stride[i],IR_9);
    SCUarg_mat[i*2] = new SCUDirArgIR;
    SCUarg_mat[i*2]->Init((void *)rcv_buf[i],rcv_dir[i],SCU_REC,non_local_chi[i]*VECT_LEN*sizeof(IFloat)*3,1,0,IR_9);
    SCUarg_mat[i*2+1] = new SCUDirArgIR;
    SCUarg_mat[i*2+1]->Init((void *)rcv_buf[i],snd_dir[i],SCU_SEND,blklen[i]*3,numblk[i],stride[i]*3,IR_9);
  }
#else
  //Loop over all possible communication directions
  for(i=0;i<2*NDIM;i++) {

      for (int hop=1; hop<=MAX_HOP; hop++) {
      //Initialize SCUArg to receive fermion fields
      SCUarg[hop-1][i*2] = new SCUDirArgIR;
      SCUarg[hop-1][i*2]->Init((void *)rcv_buf[i],rcv_dir[i],SCU_REC,
			       hop*non_local_chi[i]*VECT_LEN*sizeof(IFloat),1,0,IR_9);
      //Initialize SCUArg to send fermion field
      SCUarg[hop-1][i*2+1] = new SCUDirArgIR;
      SCUarg[hop-1][i*2+1]->Init((void *)rcv_buf[i],snd_dir[i],SCU_SEND,hop*blklen[i],
				 numblk[i],stride[i]+(1-hop)*blklen[i],IR_9);
    
      //Receive for Matrices
      SCUarg_mat[hop-1][i*2] = new SCUDirArgIR;
      SCUarg_mat[hop-1][i*2]->Init((void *)rcv_buf[i],rcv_dir[i],SCU_REC,
				   3*hop*non_local_chi[i]*VECT_LEN*sizeof(IFloat),1,0,IR_9);
      //send for matrices
      SCUarg_mat[hop-1][i*2+1] = new SCUDirArgIR;
      SCUarg_mat[hop-1][i*2+1]->Init((void *)rcv_buf[i],snd_dir[i],SCU_SEND,
				     3*hop*blklen[i],numblk[i],
				     3*(stride[i]+(1-hop)*blklen[i]),IR_9);
      //Receive
      SCUarg2[hop-1][i*2] = new SCUDirArgIR;
      SCUarg2[hop-1][i*2]->Init((void *)rcv_buf2[i],rcv_dir[i],SCU_REC,
				hop*non_local_chi[i]*VECT_LEN*sizeof(IFloat),1,0,IR_9);
      //Send
      SCUarg2[hop-1][i*2+1] = new SCUDirArgIR;
      SCUarg2[hop-1][i*2+1]->Init((void *)rcv_buf2[i],snd_dir[i],SCU_SEND,hop*blklen[i],
				  numblk[i],stride[i]+(1-hop)*blklen[i],IR_9);
      }

      //-----------------------------------------------------------------------
      //Initialize SCUArg to receive fermion fields
      SCUarg_cb[i*2] = new SCUDirArgIR;
      SCUarg_cb[i*2]->Init((void *)rcv_buf[i],rcv_dir[i],SCU_REC,
			       non_local_chi_cb[i]*VECT_LEN*sizeof(IFloat),1,0,IR_9);
      //Initialize SCUArg to send fermion field
      SCUarg_cb[i*2+1] = new SCUDirArgIR;

      if(i%2)
	SCUarg_cb[i*2+1]->Init((void *)rcv_buf[i],snd_dir[i],SCU_SEND,
					non_local_chi_cb[i]*VECT_LEN*sizeof(IFloat),1,0,IR_9);
      else 
	  SCUarg_cb[i*2+1]->Init((void *)rcv_buf[i],snd_dir[i],SCU_SEND,blklen_cb[i],
				 numblk_cb[i],stride_cb[i],IR_9);
    
      //Receive for Matrices
      SCUarg_mat_cb[i*2] = new SCUDirArgIR;
      SCUarg_mat_cb[i*2]->Init((void *)rcv_buf[i],rcv_dir[i],SCU_REC,
				   3*non_local_chi_cb[i]*VECT_LEN*sizeof(IFloat),1,0,IR_9);
      //send for matrices
      SCUarg_mat_cb[i*2+1] = new SCUDirArgIR;
      if(i%2)
	SCUarg_mat_cb[i*2+1]->Init((void *)rcv_buf[i],snd_dir[i],SCU_SEND,
					3*non_local_chi_cb[i]*VECT_LEN*sizeof(IFloat),1,0,IR_9);
      else
	SCUarg_mat_cb[i*2+1]->Init((void *)rcv_buf[i],snd_dir[i],SCU_SEND,
				     3*blklen_cb[i],numblk_cb[i],
				     3*stride_cb[i],IR_9);
      //-----------------------------------------------------------------------
  }
#endif
  qfree(rcv_mat);
//  VRB.FuncEnd("",fname);
}


//End Initialization
//----------------------------------------------------------------------------

//Parallel transport of a matrix defined on one half of the
//checkerboaded lattice
//
//Parameters
//
//n - The number of direction in which to perform the parallel transport
//mout - Result of the parallel transport, on sites with opposite parity of min
//min - Initial field, defined on sites with only one parity
//dir - a list of the n directions in which the field will be transported
//cb - Checkerboard parity of the vector min

#undef PROFILE
void PT::mat_cb(int n, IFloat **mout, IFloat **min, const int *dir, int
parity)
{
  //List of the different directions
  int wire[n];
  int i;
  //SCUDirArgs for sending and receiving in the n directions
  SCUDirArgIR *SCUarg_p[2*n];
  SCUDirArgMulti SCUmulti;
  static int call_num = 0;
//  int parity = (int) cb;
  int vlen = VECT_LEN;
  int vlen2 = VECT_LEN2;

  call_num++;
  
  //Name our function
  char *fname="pt_mat_cb()";
//  VRB.Func("",fname);
  
  //Set the transfer directions
  //If wire[i] is even, then we have communication in the negative direction
  //If wire[i] is odd, then we have communication in the positive direction
  for(i=0;i<n;i++)
    wire[i]=dir[i];

#ifdef PROFILE
  Float dtime  = - dclock();
#endif

  //If wire[i] is odd, then we have parallel transport in the
  //positive direction.  In this case, multiplication by the link matrix is
  //done before the field is transferred over to the adjacent node
  //
  //If we have transfer in the negative T direction (wire[i] = 6), then
  //we have to copy the appropriate fields to a send buffer
  for(i=0;i<n;i++)
    {
      if(wire[i]%2)
	  cmm_agg_cpp_cb(non_local_chi_cb[wire[i]],(long)uc_nl_cb_pre[parity][wire[i]/2],(long)min[i],(long)snd_buf_cb[wire[i]/2],gauge_field_addr);
      else if(wire[i] == 6)
	{
	  for(int j = 0; j < non_local_chi_cb[6];j++)
	    for(int d = 0; d < GAUGE_LEN; d++)
	      *(snd_buf_t_cb + j*GAUGE_LEN + d) = *(min[i] + 3 * *(Toffset[parity]+j)*vlen+d);
	}
    }


  for(i=0;i<n;i++)
    {
      //Calculate the starting address for the data to be sent
      IFloat *addr = min[i] + GAUGE_LEN * offset_cb[wire[i]];
      //This points to the appropriate SCUDirArg for receiving
      SCUarg_p[2*i] = SCUarg_mat_cb[2*wire[i]];
      //This points to the appropriate SCUDirArg for sending
      SCUarg_p[2*i+1] = SCUarg_mat_cb[2*wire[i]+1];

      //Set the send address
      if(wire[i]%2)
	SCUarg_p[2*i+1]->Addr((void *)snd_buf_cb[wire[i]/2]);
      else if(wire[i] == 6)
	SCUarg_p[2*i+1]->Addr((void *)snd_buf_t_cb);
      else
	SCUarg_p[2*i+1]->Addr((void *)addr);
    }

  SCUmulti.Init(SCUarg_p,2*n);

  //Begin transmission
  SCUmulti.SlowStartTrans();

  //Do local calculations
  for(i=0;i<n;i++)
    {
      cmm_agg_cpp_cb(local_chi_cb[wire[i]],(long)uc_l_cb[parity][wire[i]],(long)min[i],(long)mout[i],gauge_field_addr);
    }

  //End transmission
  SCUmulti.TransComplete();

  //If wire[i] is even, then we have transport in the negative direction
  //In this case, the vector field is multiplied by the SU(3) link matrix
  //after all communication is complete
  for(i=0;i<n;i++)
    {
      if(!(wire[i]%2))
	{
	  cmm_agg_cpp_cb(non_local_chi_cb[wire[i]],(long)uc_nl_cb[parity][wire[i]],(long)rcv_buf[wire[i]],(long)mout[i],gauge_field_addr);
	}
      //Otherwise we have parallel transport in the positive direction.
      //In this case, the received data has already been pre-multiplied
      //All we need to do is to put the transported field in the correct place
      else
	{
	  IFloat *fp0,*fp1;
	  int destination, source;
	  //Place the data in the receive buffer into the result vector
	  for(int s=0;s<non_local_chi_cb[wire[i]];s++)
	    {
	      source = uc_nl_cb[parity][wire[i]][s].src*vlen;
	      fp0 = (IFloat *)(rcv_buf[wire[i]])+3*source;
	      destination = uc_nl_cb[parity][wire[i]][s].dest*vlen2;
	      fp1 = (IFloat *)(mout[i])+3*destination;
	      for(int d = 0;d<GAUGE_LEN;d++)
		*(fp1+d) = *(fp0+d);
	    }
	}
    }
#ifdef PROFILE
  dtime +=dclock();
  print_flops("",fname,99*vol*n,dtime);
#endif
//  ParTrans::PTflops +=99*n*vol;
}

//Parallel transport of a vector defined on single parity sites
//
//Parameters
//
//n - number of directions in which to parallel transport
//vout - Transported vector
//vin - Initial vector
//dir - list of directions in which to transport the vectors
//cb - Parity of the sites where the vectors are defined
//gauge - Pointer to block of gauge fields in STAG order

//Normal parallel transport with normal gauge fields
#undef PROFILE
void PT::vec_cb(int n, IFloat **vout, IFloat **vin, const int *dir, int
 parity)
{
  PT::vec_cb_norm(n,vout,vin,dir,parity,gauge_field_addr);
}

//Normal parallel transport, but with user-specified gauge fields
#undef PROFILE
void PT::vec_cb(int n, IFloat **vout, IFloat **vin, const int *dir, int
parity, IFloat * new_gauge_field)
{
  vec_cb_norm(n,vout,vin,dir,parity,new_gauge_field);
}

//Padded parallel transport with normal gauge fields
#undef PROFILE
void PT::vec_cb(int n, IFloat *vout, IFloat **vin, const int *dir, int
parity, int pad)
{
  vec_cb_pad(n,vout,vin,dir,parity,pad,gauge_field_addr);
}

//Padded parallel transport, but with user-specified gauge fields
#undef PROFILE
void PT::vec_cb(int n, IFloat *vout, IFloat **vin, const int *dir, int
parity, int pad, IFloat * new_gauge_field)
{
  vec_cb_pad(n,vout,vin,dir,parity,pad,new_gauge_field);
}

#undef PROFILE
void PT::vec_cb_norm(int n, IFloat **vout, IFloat **vin, const int *dir,int parity, IFloat * gauge)
{
  //List of the different directions
  int wire[n];
  int i;
  //SCUDirArgs for sending and receiving in the n directions
  SCUDirArgIR *SCUarg_p[2*n];
  SCUDirArgMulti SCUmulti;
  static int call_num = 0;
  int vlen = VECT_LEN;
  int vlen2 = VECT_LEN;

  #if 0
  for(int s = 0; s < GJP.VolNodeSites(); s++)
      for(int t = 0; t < 4; t++)
  	{
  	  printf("site = %d, direction = %d\n",s,t);
  	  for(int u = 0; u < 9; u++)
  	    printf("%e %e\n",*(gauge+4*GAUGE_LEN*s + GAUGE_LEN*t + 2*u),*(gauge+4*GAUGE_LEN*s + GAUGE_LEN*t + 2*u+1));
  	}
  #endif

  call_num++;
  
  //Name our function
  char *fname="pt_1vec_cb_norm()";
//  VRB.Func("",fname);
  
  //Set the transfer directions
  //If wire[i] is even, then we have communication in the negative direction
  //If wire[i] is odd, then we have communication in the positive direction
  for(i=0;i<n;i++)
    wire[i]=dir[i];

#ifdef PROFILE
  Float dtime  = - dclock();
#endif

  //If wire[i] is odd, then we have parallel transport in the
  //positive direction.  In this case, the matrix multiplication is
  //done before the field is transferred over to the adjacent node
  //
  //If we have transfer in the negative T direction (wire[i] = 6) then
  //we have to copy the appropriate fields into the send buffer
  for(i=0;i<n;i++)
    {
      if(wire[i]%2)
	{
	  //printf("dir = %d, pre-mulitply\n", wire[i]);
	  cmv_agg_cpp_cb(non_local_chi_cb[wire[i]],(long)uc_nl_cb_pre[parity][wire[i]/2],(long)vin[i],(long)snd_buf_cb[wire[i]/2],gauge);
	}
      else if(wire[i] == 6)
	{
	  for(int j = 0; j < non_local_chi_cb[6];j++)
	    for(int d = 0; d < VECT_LEN; d++)
	      *(snd_buf_t_cb + j*VECT_LEN + d) = *(vin[i] + *(Toffset[parity]+j)*vlen+d);
	}
    }


  for(i=0;i<n;i++)
    {
      //Calculate the starting address for the data to be sent
      IFloat *addr = vin[i] + VECT_LEN * offset_cb[wire[i]];
      //This points to the appropriate SCUDirArg for receiving
      SCUarg_p[2*i] = SCUarg_cb[2*wire[i]];
      //This points to the appropriate SCUDirArg for sending
      SCUarg_p[2*i+1] = SCUarg_cb[2*wire[i]+1];

      //Set the send address
      if(wire[i]%2)
	SCUarg_p[2*i+1]->Addr((void *)snd_buf_cb[wire[i]/2]);
      else if(wire[i] == 6)
	SCUarg_p[2*i+1]->Addr((void *)snd_buf_t_cb);
      else
	SCUarg_p[2*i+1]->Addr((void *)addr);
    }

  SCUmulti.Init(SCUarg_p,2*n);

  //Begin transmission
  SCUmulti.SlowStartTrans();

  //for(int i = 0; i < 4; i++)
  //  for(int j = 0; j < VECT_LEN*non_local_chi_cb[2*i+1]; j++)
  //    {
  //	printf("snd_buf_cb[%d][%d] = %e\n",i,j,*(snd_buf_cb[i]+j));
  //	printf("rcv_buf[%d][%d] = %e\n",(2*i+1),j,*(rcv_buf[2*i+1]+j));
  //   }

  //Do local calculations
  for(i=0;i<n;i++)
    {
	cmv_agg_cpp_cb(local_chi_cb[wire[i]],(long)uc_l_cb[parity][wire[i]],(long)vin[i],(long)vout[i],gauge);
    }

  //End transmission
  SCUmulti.TransComplete();

  //If wire[i] is even, then we have transport in the negative direction.
  //In this case, the vector field is multiplied by the SU(3) link matrix
  //after all communication is complete
  for(i=0;i<n;i++)
    {
      if(!(wire[i]%2))
	{

	  #if 0
	  int destination, source;
	  printf("wire[%d] = %d\n",i,wire[i]);
	  for(int s = 0; s< non_local_chi_cb[wire[i]];s++)
	    {
	      source = ((int)(s*VECT_LEN*sizeof(IFloat)/blklen_cb[wire[i]]))*(stride_cb[wire[i]]+blklen_cb[wire[i]])/(sizeof(IFloat))+(s%(blklen_cb[wire[i]]/(VECT_LEN*sizeof(IFloat))))*VECT_LEN;
	      destination = s*VECT_LEN;
	      	      for(int d = 0;d<VECT_LEN;d++)
	      {
	        printf("*(vin[%d]+%d+%d) = %e\n",i,source,d,*((IFloat *)vin[i] + source+d));
	        printf("*(rcv_buf[%d]+%d+%d) = %e\n",wire[i],destination,d,*((IFloat *)rcv_buf[wire[i]]+destination+d));
	       }
	   }
	  #endif
	    cmv_agg_cpp_cb(non_local_chi_cb[wire[i]],(long)uc_nl_cb[parity][wire[i]],(long)rcv_buf[wire[i]],(long)vout[i],gauge);
	}
      //Otherwise we have parallel transport in the positive direction.
      //In this case, the received data has already been pre-multiplied
      //All we need to do is to put the transported field in the correct place
      else
	{
	  IFloat *fp0,*fp1;
	  int destination, source;
	  //Place the data in the receive buffer into the result vector
	  for(int s=0;s<non_local_chi_cb[wire[i]];s++)
	    {
	      source = uc_nl_cb[parity][wire[i]][s].src*vlen;
	      destination = uc_nl_cb[parity][wire[i]][s].dest*vlen2;
	      fp0 = (IFloat *)(rcv_buf[wire[i]])+source;
	      fp1 = (IFloat *)(vout[i])+destination;
	      for(int d = 0;d<VECT_LEN;d++)
		*(fp1+d) = *(fp0+d);

	    }
	}
    }
#ifdef PROFILE
  dtime +=dclock();
  print_flops("",fname,33*vol*n,dtime);
#endif
//  ParTrans::PTflops +=33*n*vol;
}

#undef PROFILE
void PT::vec_cb_pad(int n, IFloat *vout, IFloat **vin, const int *dir,int
parity,int pad, IFloat * gauge)
{
  //List of the different directions
  int wire[n];
  int i;
  //SCUDirArgs for sending and receiving in the n directions
  SCUDirArgIR *SCUarg_p[2*n];
  SCUDirArgMulti SCUmulti;
  static int call_num = 0;
  int vlen = VECT_LEN;
  int vlen2 = VECT_LEN;
  if(pad)
    vlen2 = 8;

  #if 0
  for(int s = 0; s < GJP.VolNodeSites(); s++)
      for(int t = 0; t < 4; t++)
  	{
  	  printf("site = %d, direction = %d\n",s,t);
  	  for(int u = 0; u < 9; u++)
  	    printf("%e %e\n",*(gauge+4*GAUGE_LEN*s + GAUGE_LEN*t + 2*u),*(gauge+4*GAUGE_LEN*s + GAUGE_LEN*t + 2*u+1));
  	}
  #endif

  call_num++;
  
  //Name our function
  char *fname="pt_1vec_cb_pad()";
//  VRB.Func("",fname);
  
  //Set the transfer directions
  //If wire[i] is even, then we have communication in the negative direction
  //If wire[i] is odd, then we have communication in the positive direction
  for(i=0;i<n;i++)
    wire[i]=dir[i];

#ifdef PROFILE
  Float dtime  = - dclock();
#endif

  //If wire[i] is odd, then we have parallel transport in the
  //positive direction.  In this case, the matrix multiplication is
  //done before the field is transferred over to the adjacent node
  //
  //If we have transfer in the negative T direction (wire[i] = 6) then
  //we have to copy the appropriate fields into the send buffer
  for(i=0;i<n;i++)
    {
      if(wire[i]%2)
	{
	  cmv_agg_cpp_cb(non_local_chi_cb[wire[i]],(long)uc_nl_cb_pre[parity][wire[i]/2],(long)vin[i],(long)snd_buf_cb[wire[i]/2],gauge,0);
	}
      else if(wire[i] == 6)
	{
	  for(int j = 0; j < non_local_chi_cb[6];j++)
	    for(int d = 0; d < VECT_LEN; d++)
	      *(snd_buf_t_cb + j*VECT_LEN + d) = *(vin[i] + *(Toffset[parity]+j)*vlen+d);
	}
    }
  for(i=0;i<n;i++)
    {
      //Calculate the starting address for the data to be sent
      IFloat *addr = vin[i] + VECT_LEN * offset_cb[wire[i]];
      //This points to the appropriate SCUDirArg for receiving
      SCUarg_p[2*i] = SCUarg_cb[2*wire[i]];
      //This points to the appropriate SCUDirArg for sending
      SCUarg_p[2*i+1] = SCUarg_cb[2*wire[i]+1];

      //Set the send address
      if(wire[i]%2)
	SCUarg_p[2*i+1]->Addr((void *)snd_buf_cb[wire[i]/2]);
      else if(wire[i] == 6)
	SCUarg_p[2*i+1]->Addr((void *)snd_buf_t_cb);
      else
	SCUarg_p[2*i+1]->Addr((void *)addr);
    }

  SCUmulti.Init(SCUarg_p,2*n);

  //Begin transmission
  SCUmulti.SlowStartTrans();

  //  for(int i = 0; i < 4; i++)
  //  for(int j = 0; j < VECT_LEN*non_local_chi_cb[2*i+1]; j++)
  //    {
  //	printf("snd_buf_cb[%d][%d] = %e\n",i,j,*(snd_buf_cb[i]+j));
  //	printf("rcv_buf[%d][%d] = %e\n",(2*i+1),j,*(rcv_buf[2*i+1]+j));
  //   }

  //Do local calculations
  for(i=0;i<n;i++)
    {
      //printf("local, dir = %d\n",i);
      cmv_agg_cpp_cb(local_chi_cb[wire[i]],(long)uc_l_cb[parity][wire[i]],(long) vin[i],(long)vout,gauge,pad);
    }

  //End transmission
  SCUmulti.TransComplete();

  //If wire[i] is even, then we have transport in the negative direction.
  //In this case, the vector field is multiplied by the SU(3) link matrix
  //after all communication is complete
  for(i=0;i<n;i++)
    {
      if(!(wire[i]%2))
	{

	  #if 0
	  int destination, source;
	  printf("wire[%d] = %d\n",i,wire[i]);
	  for(int s = 0; s< non_local_chi_cb[wire[i]];s++)
	    {
	      source = ((int)(s*VECT_LEN*sizeof(IFloat)/blklen_cb[wire[i]]))*(stride_cb[wire[i]]+blklen_cb[wire[i]])/(sizeof(IFloat))+(s%(blklen_cb[wire[i]]/(VECT_LEN*sizeof(IFloat))))*VECT_LEN;
	      destination = s*VECT_LEN;
	      	      for(int d = 0;d<VECT_LEN;d++)
	      {
	        printf("*(vin[%d]+%d+%d) = %e\n",i,source,d,*((IFloat *)vin[i] + source+d));
	        printf("*(rcv_buf[%d]+%d+%d) = %e\n",wire[i],destination,d,*((IFloat *)rcv_buf[wire[i]]+destination+d));
	       }
	   }
	  #endif
	  //printf("post-multiply, dir = %d\n",i);
	    cmv_agg_cpp_cb(non_local_chi_cb[wire[i]],(long)uc_nl_cb[parity][wire[i]],(long)rcv_buf[wire[i]],(long)vout,gauge,pad);
	}
      //Otherwise we have parallel transport in the positive direction.
      //In this case, the received data has already been pre-multiplied
      //All we need to do is to put the transported field in the correct place
      else
	{
	  IFloat *fp0,*fp1;
	  int destination, source;
	  //Place the data in the receive buffer into the result vector
	  for(int s=0;s<non_local_chi_cb[wire[i]];s++)
	    {
	      source = uc_nl_cb[parity][wire[i]][s].src*vlen;
	      fp0 = (IFloat *)(rcv_buf[wire[i]])+source;
		  destination = uc_nl_cb[parity][wire[i]][s].dest2*vlen2;
		  fp1 = vout+destination;
	      for(int d = 0;d<VECT_LEN;d++)
		*(fp1+d) = *(fp0+d);

	    }
	}
    }
#ifdef PROFILE
  dtime +=dclock();
  print_flops("",fname,33*vol*n,dtime);
#endif
//  ParTrans::PTflops +=33*n*vol;
}


//-----------------------------------------------------------------------------

//Parallel transport of a matrix. through one hop.
//The matrix min is parallel transported and the result is placed in mout
#undef PROFILE
void PT::mat(int n, matrix **mout, matrix **min, const int *dir){
    
  int wire[n];
  int i;
  SCUDirArgIR *SCUarg_p[2*n];
  SCUDirArgMulti SCUmulti;
  static int call_num = 0;

  call_num++;
  char *fname="pt_mat()";
//  VRB.Func("",fname);
	
  for(i=0;i<n;i++) wire[i] = dir[i]; 
#ifdef PROFILE
  Float dtime  = - dclock();
#endif
  for(i=0;i<n;i++) {
    //Calculate the address for transfer in a particular direction
    Float * addr = ((Float *)min[i]+GAUGE_LEN*offset[wire[i]]);
    //This should point to the appropriate SCUDirArg for receiving
    SCUarg_p[2*i] = SCUarg_mat[0][2*wire[i]];
    //This points to the appropriate SCUDirArg for sending
    SCUarg_p[2*i+1] = SCUarg_mat[0][2*wire[i]+1];
    //Reset the send address
    SCUarg_p[2*i+1]->Addr((void *)addr);
  }
  SCUmulti.Init(SCUarg_p,n*2);
  //Start transmission
  SCUmulti.SlowStartTrans();

  //Initerleaving of local computation of matrix multiplication
  for(i=0;i<n;i++){
#ifdef CPP
    cmm_agg_cpp(local_chi[wire[i]],0, (long)uc_l[wire[i]], (long)min[i],(long)mout[i]);
#else
    cmm_agg(uc_l[wire[i]],min[i],mout[i],local_chi[wire[i]]/2);
#endif
  }
  //End transmission
  SCUmulti.TransComplete();
  //Do non-local computations
  for(i=0;i<n;i++){ 
#ifdef CPP
    cmm_agg_cpp(non_local_chi[wire[i]],0, (long)uc_nl[wire[i]], (long)rcv_buf[wire[i]],(long)mout[i]);
#else
    cmm_agg(uc_nl[wire[i]],(matrix *)rcv_buf[wire[i]],mout[i],non_local_chi[wire[i]]/2);
#endif
  }
#ifdef PROFILE
  dtime +=dclock();
  print_flops("",fname,198*vol*n,dtime);
#endif
//  ParTrans::PTflops +=198*n*vol;
}


#undef PROFILE
//Parallel transport of a vector through one hop
void PT::vec(int n, IFloat **vout, IFloat **vin, const int *dir){
  int i;
  static int call_num=0;
  SCUDirArgIR *SCUarg_p[2*n];
  call_num++;
  //for(int s = 0; s < GJP.VolNodeSites(); s++)
  //  {
  //    for(int t = 0; t < 4; t++)
  //	{
  //	  printf("site = %d, direction = %d\n",s,t);
  //	  for(int u = 0; u < 9; u++)
  //	    printf("%e %e\n",*(gauge_field_addr+4*GAUGE_LEN*s + GAUGE_LEN*t + 2*u),*(gauge_field_addr+4*GAUGE_LEN*s + GAUGE_LEN*t + 2*u+1));
  //	}
  //  }

#ifdef PROFILE
  Float dtime  = - dclock();
#endif
  int wire[n];
  SCUDirArgMulti SCUmulti;

  char *fname="pt_1vec";
//  VRB.Func("",fname);
	
  for(i=0;i<n;i++) wire[i] = dir[i]; // from (x,y,z,t) to (t,x,y,z)
  for(i=0;i<n;i++) {
    IFloat * addr = (vin[i]+VECT_LEN*offset[wire[i]]);
    SCUarg_p[2*i] = SCUarg[0][2*wire[i]];
    SCUarg_p[2*i+1] = SCUarg[0][2*wire[i]+1];
    SCUarg_p[2*i+1]->Addr((void *)addr);
  }
  SCUmulti.Init(SCUarg_p,n*2);
  SCUmulti.SlowStartTrans();
	
#ifndef CPP
  for(i=0;i<n;i++) pt_asqtad_agg(local_chi[wire[i]],0, (long)uc_l[wire[i]], (long)vin[i],(long)vout[i]);
#else
  for(i=0;i<n;i++) 
    cmv_agg_cpp(local_chi[wire[i]],(long)uc_l[wire[i]], (long)vin[i],(long)vout[i]);
#endif
  SCUmulti.TransComplete();
	
#ifndef CPP
  for(i=0;i<n;i++) pt_asqtad_agg(non_local_chi[wire[i]],0, (long)uc_nl[wire[i]], (long)rcv_buf[wire[i]],(long)vout[i]);
#else
  for(i=0;i<n;i++) 
    cmv_agg_cpp(non_local_chi[wire[i]],(long)uc_nl[wire[i]], (long)rcv_buf[wire[i]],(long)vout[i]);
#endif

#ifdef PROFILE
  dtime +=dclock();
  print_flops("",fname,66*n*vol,dtime);
#endif
//  ParTrans::PTflops +=66*n*vol;
}

/*! 
  Computes sum[x] = vect[x] vect[x + hop dir]^dagger
  where the sum is over n_vect vectors and the hop is in a forward direction.
*/
void PT::vvpd(IFloat **vect, int n_vect, const int *dir, 
	     int n_dir, int hop, IFloat **sum){
  char *fname = "pt_vvpd()";
//  VRB.Func("",fname);
  int i, s, v;
  Float f = 2.0;
  int wire[n_dir];
  for(i=0;i<n_dir;i++) wire[i] = dir[i]; // from (x,y,z,t) to (t,x,y,z)

  SCUDirArgIR *SCUarg_p[2*n_dir];  
  SCUDirArgIR *SCUarg_p2[2*n_dir];  

  // Only do communication in forward direction
  for(i=0;i<n_dir;i++) {
    if ( size[wire[i]] <hop)
      fprintf(stderr, 
		"%s:size(%d) in direction %d is smaller than the hop(%d)\n",
		fname,size[wire[i]],wire[i],hop);
    SCUarg_p[2*i] = SCUarg[hop-1][4*wire[i]];
    SCUarg_p[2*i+1] = SCUarg[hop-1][4*wire[i]+1];
    SCUarg_p2[2*i] = SCUarg2[hop-1][4*wire[i]];
    SCUarg_p2[2*i+1] = SCUarg2[hop-1][4*wire[i]+1];
  }

  for(v=0; v<n_vect; v++){
    SCUDirArgMulti SCUmulti;

    if (v%2==0) {
      for(i=0;i<n_dir;i++)
	SCUarg_p[2*i+1]->Addr((void *)(vect[v]+set_offset(2*wire[i], hop)));

      // Start communication
      SCUmulti.Init(SCUarg_p,2*n_dir);
    } else {
      for(i=0;i<n_dir;i++)
	SCUarg_p2[2*i+1]->Addr((void *)(vect[v]+set_offset(2*wire[i], hop)));

      // Start communication
      SCUmulti.Init(SCUarg_p2,2*n_dir);
    }
    SCUmulti.SlowStartTrans();

    // Perform non-local calculation for previous v
    if (v>0)
      if (v==1) {
	for(i=0; i<n_dir; i++) 
	  cross_over_lin(sum[i], &f, vect[v-1],rcv_buf[2*wire[i]], hop*non_local_chi[2*wire[i]],
		src_nl[hop-1][2*wire[i]], dest_nl[hop-1][2*wire[i]]);
      } else if (v%2==1) {
	for(i=0; i<n_dir; i++) 
	  cross_lin(sum[i], &f, vect[v-1],rcv_buf[2*wire[i]], hop*non_local_chi[2*wire[i]],
		src_nl[hop-1][2*wire[i]], dest_nl[hop-1][2*wire[i]]);
      } else {
	for(i=0; i<n_dir; i++) 
	  cross_lin(sum[i], &f,vect[v-1],rcv_buf2[2*wire[i]], hop*non_local_chi[2*wire[i]],
		src_nl[hop-1][2*wire[i]], dest_nl[hop-1][2*wire[i]]);
      }
    
    // Perform local calculation for current v
    if (v==0)
      for(i=0; i<n_dir; i++)
	cross_over_look(sum[i], &f, vect[v], vect[v], vol-hop*non_local_chi[2*wire[i]], 
			src_l[hop-1][2*wire[i]], dest_l[hop-1][2*wire[i]]);
    else
      for(i=0; i<n_dir; i++)
	cross_look(sum[i], &f, vect[v], vect[v], vol-hop*non_local_chi[2*wire[i]], 
	      src_l[hop-1][2*wire[i]], dest_l[hop-1][2*wire[i]]);

    // Finalise communication
    SCUmulti.TransComplete();

  }

  if (v==1) {
    for(i=0; i<n_dir; i++) 
      cross_over_lin(sum[i], &f, vect[v-1],rcv_buf[2*wire[i]], hop*non_local_chi[2*wire[i]],
	    src_nl[hop-1][2*wire[i]], dest_nl[hop-1][2*wire[i]]);
  } else if (v%2==1) {
    for(i=0; i<n_dir; i++) 
      cross_lin(sum[i], &f, vect[v-1],rcv_buf[2*wire[i]], hop*non_local_chi[2*wire[i]],
	    src_nl[hop-1][2*wire[i]], dest_nl[hop-1][2*wire[i]]);
  } else {
    for(i=0; i<n_dir; i++) 
      cross_lin(sum[i], &f,vect[v-1],rcv_buf2[2*wire[i]], hop*non_local_chi[2*wire[i]],
	    src_nl[hop-1][2*wire[i]], dest_nl[hop-1][2*wire[i]]);
  }
  
//  ParTrans::PTflops += 90*n_vect*n_dir*vol;
  
}

//! u[x] = v[x+dir] for n_dir forward or backward directions dir.
void PT::shift_field(IFloat **v, const int *dir, int n_dir,
		    int hop, IFloat **u){

  int i, length;
  int wire[n_dir];
  for (i=0; i<n_dir;i++) wire[i] = dir[i];
  SCUDirArgMulti SCUmulti;
  SCUDirArgIR *SCUarg_p[2*n_dir];

  for (i=0; i<n_dir; i++) {
    SCUarg_p[2*i] = SCUarg_mat[hop-1][2*wire[i]];
    SCUarg_p[2*i+1] = SCUarg_mat[hop-1][2*wire[i]+1];
    SCUarg_p[2*i+1]->Addr((void *)(v[i]+set_offset(wire[i], hop)));
  }

  SCUmulti.Init(SCUarg_p,2*n_dir);
  SCUmulti.SlowStartTrans();

  for (i=0; i<n_dir; i++) {
    length = vol-hop*non_local_chi[wire[i]];
    copy_matrix(u[i],v[i],&length,dest_l[hop-1][wire[i]],
		src_l[hop-1][wire[i]]);
  }

  SCUmulti.TransComplete();
  
  for (i=0; i<n_dir; i++) {
    length = hop*non_local_chi[wire[i]];
    copy_matrix(u[i],(IFloat*)rcv_buf[wire[i]],&length,
		dest_nl[hop-1][wire[i]],src_nl[hop-1][wire[i]]);
  }

}

//! u[-/+nu](x) = U_[-/+nu](x) 
void PT::shift_link(IFloat **u, const int *dir, int n_dir){

  char *fname = "pt_shift_link()";
//  VRB.Func("",fname);
  int length;
  for (int i=0; i<n_dir; i++) {
    
    length = local_chi[dir[i]];
    copy_gauge(u[i],uc_l[dir[i]],&length,dest_l[0][dir[i]]);
    length = non_local_chi[dir[i]];
    copy_gauge(u[i],uc_nl[dir[i]],&length,dest_nl[0][dir[i]]);
    
  }

}

CPS_END_NAMESPACE
