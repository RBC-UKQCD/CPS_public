#ifndef _PT_INT_H_
#define _PT_INT_H_
#include <config.h>
#include <util/gauge_agg.h>
//#include <util/lattice.h>
#include <util/amalloc.h>
#include <util/asqtad_int.h>
/*!\file
  \brief Declaration of functions used by the parallel transport classes.

  $Id: pt_int.h,v 1.10 2005-04-05 06:44:43 chulwoo Exp $
  Why are (at least some of) these not class methods?
*/
//#include <util/lattice.h>
#if 0
CPS_START_NAMESPACE
void pt_init(Lattice &lat);  //!< Initialization for parallel transporters
void pt_init_g();
void pt_delete();
void pt_delete_g();
void pt_mat(int n, Float **mout, Float **min, const int *dir);
void pt_1vec(int n, Float **vout, Float **vin, int const *dir);
void pt_2vec(int n, Float **vout, Float **vin, const int *dir);
//void pt_set_hop_pointer();
int pt_offset(int dir, int hop);

void pt_vvpd(Float **vect, int n_vect, const int *dir,
	     int n_dir, int hop, Float **sum);
void pt_shift_field(Float **v, const int *dir, int n_dir,
		    int hop, Float **u);
void pt_shift_link(Float **u, const int *dir, int n_dir);

//---------------------------------------------------------------
//Checkerboarding methods
void pt_mat_cb(int n, Float **mout, Float **min, const int *dir, ChkbType cb);  //!<Parallel transport for checkerboarded Matrix fields
void pt_mat_cb(int n, Float **mout, Float **min, const int *dir, ChkbType cb,new_gauge_field);  //!<Parallel transport for checkerboarded Matrix fields

void pt_1vec_cb(int n, Float **vout, Float **vin, const int *dir, ChkbType cb); //!<Parallel transport for checkerboarded Vector fields
void pt_1vec_cb(int n, Float **vout, Float **vin, const int *dir, ChkbType cb, Float * new_gauge_field); //!<Parallel transport for checkerboarded Vector fields
void pt_1vec_cb(int n, Float *vout, Float **vin, const int *dir, ChkbType cb, int pad); //!<Parallel transport for padded checkerboarded Vector fields
void pt_1vec_cb(int n, Float *vout, Float **vin, const int *dir, ChkbType cb, int pad, Float * new_gauge_field); //!<Parallel transport for padded checkerboarded Vector fields
void pt_1vec_cb_norm(int n, Float **vout, Float **vin, const int *dir,ChkbType cb, Float * gauge);
void pt_1vec_cb_pad(int n, Float *vout, Float **vin, const int *dir,ChkbType cb,int pad, Float * gauge);
//---------------------------------------------------------------
#endif


#if 0
struct gauge_agg{
  int src;
  int dest;
  Float mat[18];
};

//-----------------------------------------------------------------------
//Checkerboarded parallel transport

struct gauge_agg_cb{
  //Index for initial position of field
  int src;
  //Index for "transported field"
  int dest;
  //Index for the padded field.  This assumes that the field will be
  //transported in all 8 possible directions.  Instead of storing
  //each direction in a single block, fields transported to each
  //lattice site from different directions are stored together
  //This allows for faster linear combination in the p4 action.
  int dest2;
  //Index for the gauge link
  int gauge_index;
  //Determines if the gauge link needs to be complex conjugated
  int dagger;
};

//-----------------------------------------------------------------------

struct hop_pointer {
  int src;
  int dest;
};
#endif

#include <qcdocos/scu_dir_arg.h>
struct PTArg {
  int size[4];
  int local[4];
  Float * gauge_field_addr;
  int g_str_ord;
  int g_conj;
  int v_str_ord;
  int v_str_ord_cb;
  int evenodd;
  int prec;
};

enum {
PT_DAG_YES=1,PT_DAG_NO=0,
PT_CHKB_EVEN=0,PT_CHKB_ODD=1,
PT_XYZT = 0, PT_TXYZ, PT_XYZT_CB_O,
PT_EVEN, PT_ODD
};



class PT  {
  private:
  char *cname;
  static int size[4];
  int local[4];
  int g_str_ord;
  int g_conj;
  int v_str_ord;
  int v_str_ord_cb;
  int evenodd;
  int prec;
  int non_local_dirs;
  IFloat * gauge_txyz;
	enum {NDIM = 4,
		  MAX_HOP = 3,
    	  GAUGE_LEN=18,
		  VECT_LEN=6, 
		  VECT_LEN2=6};
//gauge_agg holds source and destination indices, as well as the SU(3)
//link matrix.  One for local parallel transport, another for non-local
    gauge_agg *uc_l[2*NDIM];
    gauge_agg *uc_nl[2*NDIM];

//---------------------------------------------------------------------------
//Holds source,destination indexes for the matrix multiplication.
//Also holds index for gauge field, and whether gauge field needs to be
//conjugated
//First array index = 0 for even parity block, =1 for odd parity block
//
//uc_nl_cb_pre holds information for pre-multiplication of the fields
//that are transported in the positive direction.

    gauge_agg_cb *uc_l_cb[2][2*NDIM];
    gauge_agg_cb *uc_nl_cb[2][2*NDIM];
    gauge_agg_cb *uc_nl_cb_pre[2][NDIM];

    gauge_agg_cb *uc_l_pad_cb[2][2*NDIM];
    gauge_agg_cb *uc_nl_pad_cb[2][2*NDIM];
    gauge_agg_cb *uc_nl_pad_cb_pre[2][NDIM];

//--------------------------------------------------------------------------

    hop_pointer *hp_l[MAX_HOP][2*NDIM];
    hop_pointer *hp_nl[MAX_HOP][2*NDIM];

    unsigned long *src_l[MAX_HOP][2*NDIM];
    unsigned long *dest_l[MAX_HOP][2*NDIM];
    unsigned long *src_nl[MAX_HOP][2*NDIM];
    unsigned long *dest_nl[MAX_HOP][2*NDIM];

//Length of block of data for SCU communication
    int blklen[2*NDIM];

//Number of blocks of data
    int numblk[2*NDIM];

//Stride between blocks
    int stride[2*NDIM];

//number of parallel transports that can be done locally
    int local_chi[2*NDIM];

//Parallel transports that require non-local communication
    int non_local_chi[2*NDIM];

//Initial offset for the data when using SCU communication
    int offset[2*NDIM];

//--------------------------------------------------------------
//Checkerboarded data
    int blklen_cb[2*NDIM];
    int numblk_cb[2*NDIM];
    int stride_cb[2*NDIM];
    int local_chi_cb[2*NDIM];
    int non_local_chi_cb[2*NDIM];
    int offset_cb[2*NDIM];

//This determines whether the gauge links that are stored
//are the normal, canonical gauge links U_mu(x)
//or the conjugated versions, U_mu(x).dagger()
//For the staggered storage, the links are conjugated.

int conjugated;
//--------------------------------------------------------------

//Buffer for receiving data via SCU
    Float *rcv_buf[2*6];
    Float *rcv_buf2[2*6];

//--------------------------------------------------------------
//Send buffer transfers in the positive direction
    Float *snd_buf_cb[6];

//Buffer for transfer in the negative T direction.  This is needed because
//the block-stride communication does not work for communication
//in the T direction on a checkerboard lattice
    Float *snd_buf_t_cb;

//List of indexes for the vectors that are transferred when
//communicating in the negative T direction
    int *Toffset[2];
//--------------------------------------------------------------

//Pointer to the gauge field
    Float *gauge_field_addr;

//SCU communication parameters
    SCUDirArgIR *SCUarg[MAX_HOP][4*NDIM];
    SCUDirArgIR *SCUarg_mat[MAX_HOP][4*NDIM];
    SCUDirArgIR *SCUarg2[MAX_HOP][4*NDIM];

//--------------------------------------------------------------
//Checkerboarded SCU
    SCUDirArgIR *SCUarg_cb[4*NDIM];
    SCUDirArgIR *SCUarg_mat_cb[4*NDIM];
//--------------------------------------------------------------

//Function primitives
    void (*Copy) (Float *dest, Float *src);
    void (*DagCopy) (Float *dest, Float *src);
    int (*LexVector)(int *x);
    int (*LexVector_cb)(int *x);
    int (*LexGauge) (int *x,int mu);
    int (*LexGauge2) (int *x, int mu);
    static void cpy (Float *dest, Float *src);
    static void dag_cpy (Float *dest, Float *src);
    static int lex_xyzt(int *x);
    static int lex_xyzt_cb_o(int *x);
    static int lex_txyz_cb(int *x);
    static int lex_txyz(int *x);
    int LexSurface(int *x);
    static int lex_g_xyzt(int *x, int mu);
    static int lex_g_xyzt_cb_o(int *x, int mu);
    static int lex_g_txyz(int *x, int mu);
    static int set_offset(int dir, int hop);
    void set_hop_pointer();
  public:
    PT() {};
    ~PT() {};
    static  int vol;
    void init (PTArg *pt_arg);
	void init_g();
	void delete_buf();
	void delete_g_buf();

	void mat_cb_norm(int n, Float **mout, Float **min, const int *dir, 
			 int parity, Float * gauge);
	void mat_cb(int n, Float **mout, Float **min, const int *dir, 
		    int parity, Float * new_gauge_field);
	void mat_cb(int n, Float **mout, Float **min, const int *dir, 
		    int parity);

	void vec_cb(int n, Float *vout, Float **vin, const int *dir,
		    int parity, int pad);
	void vec_cb(int n, Float **vout, Float **vin, const int *dir,
		    int parity, Float * new_gauge_field);
	void vec_cb(int n, Float **vout, Float **vin, const int *dir,
		    int parity);
	void vec_cb(int n, Float *vout, Float **vin, const int *dir,
		    int parity, int pad, Float * new_gauge_field);
	void vec_cb_norm(int n, Float **vout, Float **vin, const int
			 *dir,int parity, Float * gauge);
	void vec_cb_pad(int n, Float *vout, Float **vin, const int
			*dir,int parity, Float * gauge);

	void mat(int n, matrix **mout, matrix **min, const int *dir);
	void vec(int n, Float **vout, Float **vin, const int *dir);
	void vvpd(Float **vect, int n_vect, const int *dir,
		  int n_dir, int hop, Float **sum);
	void shift_field(Float **v, const int *dir, int n_dir,
			 int hop, Float **u);
	void shift_link(Float **u, const int *dir, int n_dir);
	void asqtad_force(AsqDArg *asq_arg, matrix *mom, Float *X, Float mass, Float dt);
	void PointerErr(char *cname, char *fname, char *vname){
	  printf("%s::%s: %s not allocated\n",cname,fname,vname);
  	  exit(-1);
    }
  void *Alloc(int request){
    return qalloc(QCOMMS,request);
  }
  void *Alloc(char *cname, char *fname, char *vname, int request){
    if (request<0){ 
      printf("Alloc(): %s::%s: %s %d bytes\n",cname,fname,vname,request);
      exit(-42);
    }
    if (request==0) return NULL;
    return qalloc(QCOMMS,request);
  }
  void Free(void *p){
    if (p) qfree(p);
  }
  void force_product_sum(vector *v, vector *w, Float coeff, matrix *f);
  void update_momenta(matrix **force, Float dt, matrix *mom);
  void asqtad_fat(AsqDArg *asq_arg, matrix *fatlink);
  void asqtad_long(AsqDArg *asq_arg, matrix *longlink, matrix *longlink_m = NULL);
};
#endif
