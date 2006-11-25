#ifndef ASQTAD_INT_H
#define ASQTAD_INT_H
#include <stdlib.h>
#ifdef SCIDAC
#include <asq_data_types.h>
#include <qcdocos/scu_dir_arg.h>
#include <gauge_agg.h>
#define USE_SCU
#else //SCIDAC
#include <config.h>
#include <util/asq_data_types.h>
#include <util/gauge_agg.h>
#ifndef USE_QMP
#define USE_SCU
#define USE_QALLOC
#endif //USE_QMP
#if TARGET == QCDOC
#define USE_QALLOC
#endif
#endif //SCIDAC

#ifdef USE_QALLOC
#include <qalloc.h>
#else
#include <malloc.h>
#endif
#ifdef USE_QMP
#include <qmp.h>
#endif

class matrix{
  private:
  Float u[18];
  public:
    matrix(){};
    ~matrix(){};
    void Dagger(const Float *a);
    void Negate(){
      for(int i =0;i<18;i++) u[i] = -u[i];
    }
	void fTimesV1Plus(Float f,matrix &m){
      for(int i =0;i<18;i++) u[i] += f*m.u[i];
	}
	matrix & operator= (Float x){
      for(int i =0;i<18;i++) u[i] =0.;
	  u[0]=x;u[8]=x;u[16]=x;
	  return *this;
    }
    void TrLessAntiHermMatrix()
{

    Float *p = (Float *)u;
    *p = *(p+8) = *(p+16)=0.;
    Float tmp = 0.5*(p[2] - p[6]);
    p[2]=tmp; p[6] = -tmp;
    tmp = 0.5*(p[3] + p[7]);
    p[3]=tmp; p[7] = tmp;
    tmp = 0.5*(p[4] - p[12]);
    p[4]=tmp; p[12] = -tmp;
    tmp = 0.5*(p[5] + p[13]);
    p[5]=tmp; p[13] = tmp;
    tmp = 0.5*(p[10] - p[14]);
    p[10]=tmp; p[14] = -tmp;
    tmp = 0.5*(p[11] + p[15]);
    p[11]=tmp; p[15] = tmp;

    static Float inv3 = 0.3333333333333333333333333;
    IFloat c = inv3 * (*(p+1) + *(p+9) + *(p+17));
    p[1] -= c;
    p[9] -= c;
    p[17] -= c;
}
};

typedef Float vector;

inline Float NormSqNode(vector *src,int size){
  Float sum = 0.;
  for(int i = 0;i<size;i++) sum += src[i]*src[i];
  return sum;
}

inline void CopyVec(vector *dest,vector *src, int size){
  for(int i = 0;i<size;i++) dest[i] = src[i];
}

inline void VecMinusEquVec(vector *dest,vector *src, int size){
  for(int i = 0;i<size;i++) dest[i] -= src[i];
}

enum AsqInvParity {EVEN=0,ODD=1,EVENODD=2};

typedef struct InvArg{
  Float mass;
  Float stop_rsd;
  int niter;
  int evenodd;
  int restart;
  double final_rsq;
  int final_iter;
  double final_sec;
  double final_flop;
} InvArg;

typedef struct AsqDArg{
  int size[4];  
  int NP[4];  
  int coor[4];  
  Float *gauge_u;
  Float c1; // one-link
  Float c2; // Naik
  Float c3; // 3staple
  Float c5; // 5staple
  Float c7; // 7staple
  Float c6; // LePage
  Float *Fat[4];
  Float *Naik[4];
  Float *NaikM[4];
} AsqDArg;

class AsqD : public AsqDArg{
  private:
    char *cname;
    Float *frm_tmp;
    matrix *fat[4];
    matrix *naik[4];
    matrix *naik_m[4];

#ifdef USE_SCU
    static SCUDir scudir[8];
#endif
    static int cg_called;
    int node_odd;
    int non_local_dirs;
    int split;
    int isplit;
    int coord[4];
    int coord_nn[4];
    int coord_knn[4];
    int vol;
    int f_size_cb;
    int non_local_chi_3[2][4];
    int non_local_chi[2];
    int local_chi[2];
    int local_chi_3[2];
    int nflush;
    int odd_num;
enum{VECT_LEN=6, VECT_LEN2=8, MATRIX_SIZE=18, SITE_LEN=72, NUM_DIR=8,
N=4};
//---------------------------------------------------------------------
//  uc_l[0] points to a cluster of matrices per even site for local 
//  computations , arranged so that all parallel transport of spinors 
//  is accomplished by the same matrix times vector function.
//  The volume is devided by 2 areas: nn and 3rd nn part.
//  uc_l[1] is the same for odd sites.
//  uc_nl[0] is the same for even non-local sites.
//  uc_nl[1] is the same for odd non-local sites.
//---------------------------------------------------------------------
    gauge_agg * uc_l_agg[2];
    gauge_agg * uc_nl_agg[2];
    Float * uc_l[2];
    Float * uc_nl[2];

//------------------------------------------------------------------
// Allocate these arrays dynamicaly once the cache, noncached eDRAM
// malloc are available (should be changed according to volume)
//------------------------------------------------------------------

    Float *tmpfrm;
    unsigned long tmpfrm2;

    Float *Tbuffer[3][2];

//-------------------------------------------------------------------
// end of stack based arrays which should be heap based
//-------------------------------------------------------------------
    int * ToffsetP[3][2];
    int * ToffsetM[3][2];
    int countP[3][2];
    int countM[3][2];
//---------------------------------------------------------------------
//  pointers to storage area for color vectors from tp, xp, yp, zp, tm,
//  xm, ym, zm (p = plus, m = minus).  Indexed as 0-7
//---------------------------------------------------------------------
    Float * chi_off_node_total;
    Float * chi_off_node[2][3][8];
    Float * chi_off_node_p[2][3][8];
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
    Float ** chi_l[2];
    Float ** chi_nl[2];

//------------------------------------------------------------------
//  Values for send and receive transfers. 0-7 correspond to
//  tp, xp, yp, zp, tm, xm, ym, zm (p = plus, m = minus).
//
//  Rarg (for SCU receives) never changes, since it receives into
//  the buffers for the specified direction.
//  SCUarg[3] is set up for different Xoffset and Toffset.
//------------------------------------------------------------------

#ifdef USE_SCU
SCUDirArgIR SCUarg[2][2*NUM_DIR];
SCUDirArgMulti SCUmulti[2];

SCUDirArgIR SCUarg_1[2][2*NUM_DIR];
SCUDirArgMulti SCUmulti_1[2];

//SCUDMAInst *SCUDMAarg_p[2][NUM_DIR*4];
SCUDMAInst SCUDMAarg[2][NUM_DIR*4];

SCUDirArgIR SCUarg_2[2][2*NUM_DIR];
SCUDirArgMulti SCUmulti_2[2];
#endif

//------------------------------------------------------------------
//  The offset into the even or odd checkboard chi's and chi3's used for a send
//  from a given face.
//------------------------------------------------------------------

    int Xoffset[3][NUM_DIR];

    unsigned long address[2];

    int comp_l[2];
    int even_l[2];
    int odd_l[2];
    int comp_nl[2];
    int comp_nl_2[2];

//-------------------------------------------------------------------
//  Given a lexical value for gauge fields, set the coordinates.
//  Return 1 if odd, 0 if even
//-------------------------------------------------------------------

    int SetCoord( int sg );

//---------------------------------------------------------------------
//  Find nearest neighbor coordinate for coordinates given.  Nearest
//  neighbor coordinates are placed in coord_nn, which are always
//  on-node.  Function returns 0 if nearest neighbor is on-node,
//  1 if off-node.
//
//  nn = 0 to 7.  tp, xp, yp, zp, tm, xm, ym, zm respectively.
//---------------------------------------------------------------------

    int CoordNN( int nn );

    int CoordkNN( int nn, int k );

//---------------------------------------------------------------------
//  Return lexical value for links from coordinates c
//---------------------------------------------------------------------

    int LexGauge( int * c );

//---------------------------------------------------------------------
//  Return lexical value for vectors from coordinates c
//---------------------------------------------------------------------

    int LexVector( int * c );

//-------------------------------------------------------------------
//  Given a coordinate and a surface ( 0 = t, 1 = x, 2 = y, 3 = x )
//  calculate the offset into the receive buffers (chi_off_node);
//-------------------------------------------------------------------

//extern "C" void vaxpy3(Matrix *a,Float *b,Matrix *c,Matrix *d,int nvec);

    int LexSurface( int * cc, int surface );
  public:
  AsqD() {cname = "AsqD";};
  ~AsqD() {};
  int Size(int dir){return size[dir];}
  int Vol() {return vol;}
  void init(AsqDArg  *arg);
  void init_g(Float *frm_p,Float **fat_p=NULL, Float **naik_p=NULL, 
    Float **naikm_p=NULL);
  void destroy_buf();
  void destroy_buf_g();
  void comm_assert();
  void dirac(Float* b, Float* a, int a_odd, int add_flag);

#if TARGET == QCDOC
  void *Alloc(int request){
    void *ptr = qalloc(QCOMMS,request);
    if (!ptr) ptr = qalloc(0,request);
    if (!ptr){ printf("AsqD::Alloc failed\n"); exit(-4);}
    return ptr;
  }
  void Free(void *p){
    qfree(p);
  }
#else
  void *Alloc(int request){
    void *ptr = malloc(request);
    if (!ptr){ printf("AsqD::Alloc failed\n"); exit(-4);}
    return ptr;
  }
  void Free(void *p){
    free(p);
  }
#endif

  void PointerErr(char *cname, char *fname, char *vname){
    printf("%s::%s: %s not allocated\n",cname,fname,vname);
    exit(-1);
  }

  int NodeOdd(){return node_odd;}
  int InvCg(InvArg *inv_arg, Float *out, Float *in, Float *true_res, int odd=0);
  void MdagM(Float *mass_sq, Float *out, Float *in, int odd, Float *dot_prd = 0);
  void Dslash(Float *out, Float *in);
  void Sum( Float *sum);
};

extern AsqD *asqd_p;
#endif
