#include <config.h>
#ifdef PROFILE
#include <stdio.h>
#endif
#include <util/gjp.h>
#include <util/time.h>
#include <comms/scu.h>
#include <comms/glb.h>
#include <util/lattice.h>
#include <util/dirac_op.h>
#include <util/vector.h>
#include <util/pt.h>
#include <qcdocos.h>
#include <qalloc.h>
#include <ppc_lib.h>

CPS_START_NAMESPACE

enum{VECT_LEN=6, VECT_LEN2=8, MATRIX_SIZE=18, SITE_LEN=72, NUM_DIR=8, N=4};

static int size[4];
static int coord[4];
static int vol;
static Vector * fermion[8];
static Vector * tmp_frm[8];
static Vector * tmp_frm2[8];
static IFloat * knight_onetwo;
static IFloat * knight_twoone;
static IFloat * smeared_onelink;
static IFloat * smeared_gauge;
static IFloat c_onelink;
static IFloat c_knight;

static int LexVector(int * x);
static int LexGauge(int * x, int dir);
static void cpy(IFloat *dest, IFloat *src,int len);

static Fp4 *lat_pt;
void set_pt(Fp4 *lat)
{
  lat_pt = lat;
}

extern "C"
void p4_dirac_init(const void * gauge_u)
{
  //Using the convention consistent for the staggered parallel transporter
  //that {0,1,2,3} = {X,Y,Z,T}
  size[0] = GJP.XnodeSites();
  size[1] = GJP.YnodeSites();
  size[2] = GJP.ZnodeSites();
  size[3] = GJP.TnodeSites();
  vol = size[0]*size[1]*size[2]*size[3];
  c_onelink = GJP.p4_KS_coeff();
  c_knight = GJP.p4_knight_coeff();
}

extern "C"
void p4_dirac_init_g()
{
  char *cname = "DiracOpP4";
  char *fname = "p4_dirac_init_g()";
  if (vol < 4096)
    {
      knight_onetwo = (IFloat *) qalloc (QFAST|QCOMMS,NUM_DIR * vol/2 * VECT_LEN2 * sizeof(IFloat));
      knight_twoone = (IFloat *) qalloc (QFAST|QCOMMS,NUM_DIR * vol/2 * VECT_LEN2 * sizeof(IFloat));
      smeared_onelink = (IFloat *) qalloc (QFAST|QCOMMS,NUM_DIR * vol/2 * VECT_LEN2 * sizeof(IFloat));
      for(int i = 0; i < 8; i++)
	{
	  tmp_frm[i] = (Vector *) qalloc (QFAST|QCOMMS,vol/2*VECT_LEN*sizeof(IFloat));
	  tmp_frm2[i] = (Vector *) qalloc (QFAST|QCOMMS,vol/2*VECT_LEN*sizeof(IFloat));
	}
    }
  else if(1)
    {;
      knight_onetwo = (IFloat *) qalloc (QCOMMS,NUM_DIR * vol/2 * VECT_LEN2 * sizeof(IFloat));
      knight_twoone = (IFloat *) qalloc (QCOMMS,NUM_DIR * vol/2 * VECT_LEN2 * sizeof(IFloat));
      smeared_onelink = (IFloat *) qalloc (QCOMMS,NUM_DIR * vol/2 * VECT_LEN2 * sizeof(IFloat));
      for(int i = 0; i < 8; i++)
	{
	  tmp_frm[i] = (Vector *) qalloc (QCOMMS,vol/2*VECT_LEN*sizeof(IFloat));
	  tmp_frm2[i] = (Vector *) qalloc (QCOMMS,vol/2*VECT_LEN*sizeof(IFloat));
	}
    }
  else
    {
      knight_onetwo = (IFloat *) smalloc (NUM_DIR * vol/2 * VECT_LEN2 * sizeof(IFloat));
      knight_twoone = (IFloat *) smalloc (NUM_DIR * vol/2 * VECT_LEN2 * sizeof(IFloat));
      smeared_onelink = (IFloat *) smalloc (NUM_DIR * vol/2 * VECT_LEN2 * sizeof(IFloat));
      for(int i = 0; i < 8; i++)
	{
	  tmp_frm[i] = (Vector *) smalloc (vol/2*VECT_LEN*sizeof(IFloat));
	  tmp_frm2[i] = (Vector *) smalloc (vol/2*VECT_LEN*sizeof(IFloat));
	}
    }
  smeared_gauge = (IFloat *) smalloc(vol*SITE_LEN*sizeof(IFloat));

  if(smeared_gauge == 0)
    ERR.Pointer(cname,fname, "smeared_gauge");
  if(smeared_onelink == 0)
    ERR.Pointer(cname,fname, "smeared_one_link");
  if(knight_onetwo == 0) 
    ERR.Pointer(cname,fname, "knight_onetwo");
  if(knight_twoone == 0) 
    ERR.Pointer(cname,fname, "knight_twoone");

  //Copy the smeared links into smeared_gauge
  lat_pt->Smear();
  Matrix *Fat = lat_pt->Fields(0);
  IFloat *fp0;
  IFloat *fp1;
  for(coord[3] = 0; coord[3] < size[3]; coord[3]++)
    for(coord[2] = 0; coord[2] < size[2]; coord[2]++)
      for(coord[1] = 0; coord[1] < size[1]; coord[1]++)
	for(coord[0] = 0; coord[0] < size[0]; coord[0]++)
	  for(int i = 0; i < N; i++)
	    {
	      fp0 = smeared_gauge + LexGauge(coord, i)*MATRIX_SIZE;
	      fp1 = (IFloat *)(Fat + i*vol + LexVector(coord));
	      cpy(fp0,fp1,MATRIX_SIZE);
	    }
	    
}

extern "C"
void p4_destroy_dirac_buf()
{
}

extern "C"
void p4_destroy_dirac_buf_g()
{
  if(1)
    {
      qfree(knight_onetwo);
      qfree(knight_twoone);
      qfree(smeared_onelink);
      sfree(smeared_gauge);
      for(int i = 0; i < 8; i++)
	{
	  qfree(tmp_frm[i]);
	  qfree(tmp_frm2[i]);
	}
    }
  else
    {
      sfree(knight_onetwo);
      sfree(knight_twoone);
      sfree(smeared_onelink);
      sfree(smeared_gauge);
      for(int i = 0; i < 8; i++)
	{
	  sfree(tmp_frm[i]);
	  sfree(tmp_frm2[i]);
	}
    }
}

static void cpy(IFloat *dest, IFloat *src, int len)
{
  for(int i = 0; i < len; i++)
    *(dest + i) = *(src + i);
}

//Return index for CANONICAL ordering of fields
static int LexVector(int * x)
{
  return x[0]+size[0]*(x[1]+size[1]*(x[2]+size[2]*x[3]));
}

//Return index for CANONICAL or STAG ordering of gauge links
static int LexGauge(int * x, int dir)
{
  return N*LexVector(x) + dir;
}

extern "C"
void p4_dirac(Vector *f_out, Vector *f_in, int cb, int dag)
{
  ParTransStaggered_cb pt(*(lat_pt));
  int dir[8] = {0,1,2,3,4,5,6,7};
  IFloat * fp0,*fp1,*fp2,*fp3,*fp4;
  int i,j,k,n = 0;
  for(int i = 0; i < 8; i++)
    fermion[i] = f_in;

  int nflops=0;
  Float dtime = -dclock();
  ParTrans::PTflops = 0;

  //First step, gather one link and two link parallel transports from
  //all eight directions
  pt.run(8,knight_onetwo,fermion,dir,(ChkbType) cb,1);

  pt.run(8,tmp_frm,fermion,dir,(ChkbType) cb);
  pt.run(8,knight_twoone,tmp_frm,dir,(ChkbType) (1-cb),1);

  //Do appropriate linear combinations for the initial one link term
  for(i = 0; i < 8; i++)
    for(n = 0; n < vol/2; n++)
      {
	fp1 = knight_onetwo + n*NUM_DIR*VECT_LEN2;
	fp0 = (IFloat *)(tmp_frm[i]) + n*VECT_LEN;

	for(j = 0; j < VECT_LEN; j++)
	  *(fp0+j) = 0;

	for(k = 0; k < VECT_LEN; k++)
	  for(j = 0; j < 8; j++)
	    if(i/2 != j/2)
	      {
		if(j%2==0)
		  *(fp0+k) += *(fp1+j*VECT_LEN2 + k);
		else
		  *(fp0+k) -= *(fp1+j*VECT_LEN2 + k);
	      }
	nflops += 8*(vol/2)*6*6;
      }

  //Transport these linear combinations by two links
  pt.run(8,tmp_frm2,tmp_frm,dir,(ChkbType) (1-cb));
  pt.run(8,knight_onetwo,tmp_frm2,dir,(ChkbType) cb, 1);

  //Do linear combination for the initial two-link terms
  for(i = 0; i < 8; i++)
    for(n = 0; n < vol/2; n++)
      {
	fp1 = knight_twoone + n*NUM_DIR*VECT_LEN2;
	fp0 = (IFloat *)(tmp_frm[i]) + n*VECT_LEN;

	for(j = 0; j < VECT_LEN; j++)
	  *(fp0+j) = 0;

	for(k = 0; k < VECT_LEN; k++)
	  for(j = 0; j < 8; j++)
	    if(i/2 != j/2)
	      *(fp0+k) += *(fp1+j*VECT_LEN2 + k);
	nflops += 8*(vol/2)*6*6;
      }

  //Transport these by one link
  pt.run(8,knight_twoone,tmp_frm,dir,(ChkbType) cb, 1);

  //Calculate the smeared one link term
  pt.run(8,smeared_onelink,fermion,dir,(ChkbType) cb,1,smeared_gauge);

  //Sum up contributions from the knight's move term
  //and the smeared one link term and place them in fout

   for(n = 0; n < vol/2; n++)
    {
      fp3 = knight_onetwo + n*NUM_DIR*VECT_LEN2;
      fp2 = smeared_onelink + n*NUM_DIR*VECT_LEN2;
      fp1 = knight_twoone + n*NUM_DIR*VECT_LEN2;
      fp0 = (IFloat *)(f_out) + n*VECT_LEN;

      for(j = 0; j < VECT_LEN; j++)
	*(fp0+j) = 0;

      for(k = 0; k < VECT_LEN; k++)
	{
	  fp4 = fp0 + k;
	  for(j = 0; j < 8; j++)
	    {
	      *(fp4) += *(fp3+j*VECT_LEN2 + k)*c_knight/2.0;
	      if(j%2==0)
		{
		  *(fp4) += *(fp1+j*VECT_LEN2 + k)*c_knight/2.0;
		  *(fp4) += *(fp2+j*VECT_LEN2 + k);
		}
	      else
		{
		  *(fp4) -= *(fp1+j*VECT_LEN2 + k)*c_knight/2.0;
		  *(fp4) -= *(fp2+j*VECT_LEN2 + k);
		}
	    }
	}
      nflops += (vol/2)*6*8*3;
    }
   dtime += dclock();
   nflops += ParTrans::PTflops;

   //Parallel Transport flops:
   //7 parallel transports * (264*vol) flops/transport = 1848*vol
   //Linear combination flops:
   //144*vol + 144*vol + 72*vol = 360*vol
   //Total flop count = 2208*vol

   DiracOp::CGflops += 7*264*vol+360*vol;
}

CPS_END_NAMESPACE
