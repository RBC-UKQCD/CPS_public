#include<config.h>
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2008-02-08 18:35:08 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/fix_gauge/main.C,v 1.9 2008-02-08 18:35:08 chulwoo Exp $
//  $Id: main.C,v 1.9 2008-02-08 18:35:08 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: main.C,v $
//  $Revision: 1.9 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/fix_gauge/main.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------

/*-----------------------------------------------------------------
 * File main.C. Version 14.1. Last modified on 97/12/18 at 19:28:16.
 *        Yuriy Zhestkov, Columbia University.
 * E-mail: zhestkov@phys.columbia.edu
 *-----------------------------------------------------------------
 */


#include<config.h>
#include<util/lattice.h>
#include<util/gjp.h>
#include<util/verbose.h>
#include<util/error.h>
#include<alg/alg_fix_gauge.h>
#include<alg/do_arg.h>
#include<alg/hmd_arg.h>
#include<alg/common_arg.h>
#include<util/random.h>
#include <math.h>

#ifdef PARALLEL
#include <comms/sysfunc_cps.h>
#endif

#ifndef M_PI            // for TARTAN compiler
#define M_PI 3.14159265358979323846
#endif

#define X_LINK 0
#define Y_LINK 1
#define Z_LINK 2
#define T_LINK 3

CPS_START_NAMESPACE
Complex I = Complex(0,1);

class XXX
{
public:

  static const int CheckFreq;
  static const int Dimension;
  static const int SiteSize;
  static const Float SmallFloat;

public:

  static int  Coor4d(int dir);

  static int  Node_Size_in_Dir(int dir);
  static int  Num_Nodes_in_Dir(int dir);
};

inline Matrix operator * (const Matrix& m1, const Matrix& m2)
{ Matrix r; r.DotMEqual(m1,m2); return r; }

//inline Matrix operator * (const Matrix& m, const Complex& c)
//{ Matrix r = m; for(int i=0; i<9; i++) ((Complex*)&r)[i] *= c; return r; }

//inline Matrix operator * (const Complex& c, const Matrix& m)
//{ return m*c; }

//inline Matrix operator + (const Matrix& m1, const Matrix& m2)
//{ Matrix r = m1; r += m2; return r; }

//inline Matrix operator - (const Matrix& m1, const Matrix& m2)
//{ Matrix r = m1; r -= m2; return r; }


#define NX (XXX::Node_Size_in_Dir(0))
#define NY (XXX::Node_Size_in_Dir(1))
#define NZ (XXX::Node_Size_in_Dir(2))
#define NT (XXX::Node_Size_in_Dir(3))

#define IND(x,y,z,t,l) ((((((t+NT)%NT)*NZ+((z+NZ)%NZ))*NY+   \
			  ((y+NY)%NY))*NX+((x+NX)%NX))*4+l)


static long seed = 314159227l;

Float ran345(long *idum);

Float rnd(Float rng)
{
//  const unsigned long msk = (1ul<<15)-1;
  Float res; // = (rand()&msk)*rng/(Float(msk)+1);
  res = rng * ran345(&seed);
  return res;
}

void srnd(long s)
{
  seed = s;
}

void p(Matrix x)
{
  for(int i=0; i<3; i++)
    {
      for(int j=0; j<3; j++)
	{
	  Complex xx = x(i,j);
	  if(fabs(real(xx))<1e-5)
	    xx=Complex(0,imag(xx));
	  if(fabs(imag(xx))<1e-5)
	    xx=Complex(real(xx),0);
	  VRB.Debug( "(%e,%e)", real(xx), imag(xx));
	  if(j<2)
	    VRB.Debug( "       ");
	}
      VRB.Debug("\n");
    }
  VRB.Debug("\n");
}

Matrix u1(Float a)
{
  Matrix U;
  U(0,0)=cos(a/2);   U(0,1)=I*sin(a/2); U(0,2)=0;
  U(1,0)=I*sin(a/2); U(1,1)=cos(a/2);   U(1,2)=0;
  U(2,0)=0;          U(2,1)=0;          U(2,2)=1;
  return U;
}

Matrix u2(Float a)
{
  Matrix U;
  U(0,0)=cos(a/2); U(0,1)=-sin(a/2); U(0,2)=0;
  U(1,0)=sin(a/2); U(1,1)=cos(a/2);  U(1,2)=0;
  U(2,0)=0;        U(2,1)=0;         U(2,2)=1;
  return U;
}

Matrix u3(Float a)
{
  Matrix U;
  U(0,0)=cos(a/2)+I*sin(a/2); U(0,1)=0;                   U(0,2)=0;
  U(1,0)=0;                   U(1,1)=cos(a/2)-I*sin(a/2); U(1,2)=0;
  U(2,0)=0;                   U(2,1)=0;                   U(2,2)=1;
  return U;
}

Matrix u4(Float a)
{
  Matrix U;
  U(0,0)=cos(a/2);   U(0,1)=0; U(0,2)=I*sin(a/2);
  U(1,0)=0;          U(1,1)=1; U(1,2)=0;
  U(2,0)=I*sin(a/2); U(2,1)=0; U(2,2)=cos(a/2);
  return U;
}

Matrix u5(Float a)
{
  Matrix U;
  U(0,0)=cos(a/2);  U(0,1)=0;  U(0,2)=sin(a/2);
  U(1,0)=0;         U(1,1)=1;  U(1,2)=0;
  U(2,0)=-sin(a/2); U(2,1)=0;  U(2,2)=cos(a/2);
  return U;
}

Matrix u6(Float a)
{
  Matrix U;
  U(0,0)=1; U(0,1)=0;          U(0,2)=0;
  U(1,0)=0; U(1,1)=cos(a/2);   U(1,2)=I*sin(a/2);
  U(2,0)=0; U(2,1)=I*sin(a/2); U(2,2)=cos(a/2);
  return U;
}

Matrix u7(Float a)
{
  Matrix U;
  U(0,0)=1; U(0,1)=0;         U(0,2)=0;
  U(1,0)=0; U(1,1)=cos(a/2);  U(1,2)=sin(a/2);
  U(2,0)=0; U(2,1)=-sin(a/2); U(2,2)=cos(a/2);
  return U;
}

Matrix u8(Float a)
{
  Matrix U;
  U(0,0)=cos(a/sqrt(3))+I*sin(a/sqrt(3));   U(0,1)=0;         U(0,2)=0;
  U(1,0)=0;        U(1,1)=cos(a/sqrt(3))+I*sin(a/sqrt(3)); U(1,2)=0;
  U(2,0)=0;        U(2,1)=0; U(2,2)=cos(2*a/sqrt(3))-I*sin(2*a/sqrt(3));
  return U;
}

Matrix (*u[8])(Float) = {u1,u2,u3,u4,u5,u6,u7,u8};

Matrix urnd(Float rng=2*M_PI)
{
  Matrix uu; uu = Complex(1);
  for(int i=0; i<8; i++)
    {
      Float r = rnd(rng);
      uu = uu*u[i](r);
    }

  Matrix D;
  D.Dagger(uu);

  return D;
}


Matrix *L, *M;
Matrix **G, **Gp;

CPS_END_NAMESPACE

USING_NAMESPACE_CPS

int main()
{
  char *cname = " ";
  char *fname = "main()";

  //----------------------------------------------------------------
  // Initializes all Global Job Parameters
  //----------------------------------------------------------------
  DoArg do_arg;

  do_arg.x_node_sites = 2;
  do_arg.y_node_sites = 5;
  do_arg.z_node_sites = 3;
  do_arg.t_node_sites = 4;
  do_arg.s_node_sites = 2;

#ifdef PARALLEL
  do_arg.x_nodes = 2; //SizeX();
  do_arg.y_nodes = 2; //SizeY();
  do_arg.z_nodes = 2; //SizeZ();
  do_arg.t_nodes = 2; //SizeT();
  do_arg.s_nodes = 1;
#else
  do_arg.x_nodes = 2;
  do_arg.y_nodes = 2;
  do_arg.z_nodes = 2;
  do_arg.t_nodes = 2;
  do_arg.s_nodes = 1;
#endif 

  do_arg.x_bc = BND_CND_PRD;
  do_arg.y_bc = BND_CND_PRD;
  do_arg.z_bc = BND_CND_PRD;
  do_arg.t_bc = BND_CND_APRD;
  do_arg.start_conf_kind = START_CONF_DISORD;
  do_arg.start_seed_kind = START_SEED_FIXED;
  do_arg.beta = 4.8;
  do_arg.dwf_height = 0.9;

  GJP.Initialize(do_arg);


  //----------------------------------------------------------------
  // After having set GJP, can call verbose routines
  //----------------------------------------------------------------

  VRB.Func(cname, fname);

  //----------------------------------------------------------------
  // Initialize argument structures
  //----------------------------------------------------------------
  CommonArg common_arg;
  HmdArg hmd_arg;
  
  hmd_arg.n_frm_masses = 1;
  hmd_arg.frm_mass[0] = 0.1;
  hmd_arg.n_bsn_masses = 0;
  hmd_arg.max_num_iter[0] = 5;
  hmd_arg.stop_rsd[0] = 1.0E-12;
  hmd_arg.step_size = 0.02;
  hmd_arg.steps_per_traj = 20;
  hmd_arg.metropolis = METROPOLIS_YES;

  //----------------------------------------------------------------
  // Run HMC Phi
  //----------------------------------------------------------------
  {



    //----------------------------------------------------------------
    // TESTING !!!
    //----------------------------------------------------------------



    int plns[] = {0};
    int npln = sizeof(plns)/sizeof(plns[0]);
    
    int t, x, y, z, ii;
    
    FixGaugeType normdir = FIX_GAUGE_COULOMB_Y;
    
    VRB.Debug( "Generating lattice . . .\n");
    
    L = (Matrix*) smalloc(4*NX*NY*NZ*NT*sizeof(Matrix));
    M = (Matrix*) smalloc(NX*NY*NZ*NT*sizeof(Matrix));
    
    int gsz = (normdir!=FIX_GAUGE_LANDAU) ? XXX::Node_Size_in_Dir(normdir) : 1;
    G = (Matrix**) smalloc(npln*sizeof(Matrix*));
    Gp = (Matrix**) smalloc(npln*sizeof(Matrix*));
    
    if((M==NULL) || (L==NULL) || (G==NULL) || (Gp==NULL))
      ERR.Pointer("","", "MLGGp");
    
    for(ii=0; ii<npln; ii++)
      {
	G[ii] =  (Matrix*) smalloc(NX*NY*NZ*NT/gsz * sizeof(Matrix));
	Gp[ii] = (Matrix*) smalloc(NX*NY*NZ*NT/gsz * sizeof(Matrix));
	if((G[ii] == NULL) || (Gp[ii]==NULL))
	  ERR.Pointer("","", "GGp[]");
      }
    
    for(t=0; t<NT; t++)
      for(x=0; x<NX; x++)
	for(y=0; y<NY; y++)
	  for(z=0; z<NZ; z++)
	    {
	      L[IND(x,y,z,t,T_LINK)] = urnd();
	      L[IND(x,y,z,t,X_LINK)] = urnd();
	      L[IND(x,y,z,t,Y_LINK)] = urnd();
	      L[IND(x,y,z,t,Z_LINK)] = urnd();
	    }
    
    {
      VRB.Debug("Fixing . . .\n");

      GwilsonFnone lat;

      lat.GaugeField(L);
      lat.FixGaugeAllocate(normdir,npln,plns);
      int itnum = lat.FixGauge(1e-11,10000);
      VRB.Debug("Iternum = %d\n", itnum);
      
      VRB.Debug("Resulting Gauge Fixing Matrices:\n\n");
    
      if(itnum > 0)
	for(int cnt=0; cnt<NX*NY*NZ*NT/gsz; cnt++)
	  p(Gp[0][cnt]=lat.FixGaugePtr()[0][cnt]);
    
      lat.FixGaugeFree();
    }

    // Norman's test

    VRB.Debug("---------- Norman's Test\n");
    
    y=5;
    for(t=0; t<NT; t++)
      for(x=0; x<NX; x++)
	for(z=0; z<NZ; z++)
	  {
	    M[IND(x,0,z,t,0)/4] = urnd();
	    Matrix D; D.Dagger(M[IND(x,0,z,t,0)/4]);
	    VRB.Debug("(%d,%d,%d,%d)\n",x,y,z,t);
	    L[IND(x,y,z,t,T_LINK)] 
	      = D * L[IND(x,y,z,t,T_LINK)];
	    L[IND(x,y,z,t,X_LINK)] 
	      = D * L[IND(x,y,z,t,X_LINK)];
	    L[IND(x,y,z,t,Y_LINK)] 
	      = D * L[IND(x,y,z,t,Y_LINK)];
	    L[IND(x,y,z,t,Z_LINK)] 
	      = D * L[IND(x,y,z,t,Z_LINK)];
	    
	    L[IND(x-1,y,z,t,X_LINK)]
	      = L[IND(x-1,y,z,t,X_LINK)] * M[IND(x,0,z,t,0)/4];
	    L[IND(x,y-1,z,t,Y_LINK)] 
	      = L[IND(x,y-1,z,t,Y_LINK)] * M[IND(x,0,z,t,0)/4];
	    L[IND(x,y,z-1,t,Z_LINK)] 
	      = L[IND(x,y,z-1,t,Z_LINK)] * M[IND(x,0,z,t,0)/4];
	    L[IND(x,y,z,t-1,T_LINK)] 
	      = L[IND(x,y,z,t-1,T_LINK)] * M[IND(x,0,z,t,0)/4];
	  }
    
    { 
      VRB.Debug("Fixing Gauge .........\n");
      GwilsonFnone lat;
      
      lat.GaugeField(L);
      lat.FixGaugeAllocate(normdir,npln,plns);
      int itnum = lat.FixGauge(1e-11,10000);
      VRB.Debug("Iternum = %d\n", itnum);
    
      if(itnum > 0)
	for(int cnt=0; cnt<NX*NY*NZ*NT/gsz; cnt++)
	  G[0][cnt]=lat.FixGaugePtr()[0][cnt];
      
      lat.FixGaugeFree();
    }
    
    VRB.Debug("Now Checking G'(n)M(n)G(n)+ = V(n) ............ \n");
    
    y=5;
    for(t=0; t<NT; t++)
      for(x=0; x<NX; x++)
	for(z=0; z<NZ; z++)
	  {
	    Matrix D; D.Dagger(Gp[0][(t*NZ+z)*NX+x]);
	    Matrix DM; DM.Dagger(M[IND(x,0,z,t,0)/4]);
	    p(G[0][(t*NZ+z)*NX+x]*DM*D);
	  }
    
    VRB.Debug("------- HOW WAS IT?????\n");
    
    // end of Norman test

    return 0;
  }

}


/****************
    VRB.Debug("Transforming the lattice to the Coulomb gauge\n");

    y=5;
    for(t=0; t<NT; t++)
      for(z=0; z<NZ; z++)
	for(x=0; x<NX; x++)
	  {
	    Matrix g = G[0][(t*NZ+z)*NX+x];
	    Matrix D; D.Dagger(g);
	    
	    Matrix gg = G[0][(t*NZ+z)*NX+size_t(x+1)%NX];
	    L[IND(x,y,z,t,X_LINK)] = gg*L[IND(x,y,z,t,X_LINK)]*D;
	    
	    gg = G[0][(t*NZ+size_t(z+1)%NZ)*NX+x];
	    L[IND(x,y,z,t,Z_LINK)] = gg*L[IND(x,y,z,t,Z_LINK)]*D;
	    
	    gg = G[0][(size_t(t+1)%NT*NZ+z)*NX+x];
	    L[IND(x,y,z,t,T_LINK)] = gg*L[IND(x,y,z,t,T_LINK)]*D;
	  }
    
    
    VRB.Debug("Fixing Gauge .........\n");
    lat.GaugeField(L);
    int itnum = lat.FixGauge(normdir,G,npln,plns);
    VRB.Debug("Iternum = %d\n", itnum);
    
    VRB.Debug("Resulting Gauge Fixing Matrices (BETTER BE 1's):\n\n");
    
    if(itnum > 0)
      for(int cnt=0; cnt<NX*NY*NZ*NT/gsz; cnt++)
	p(G[0][cnt]);
    
    VRB.Debug("Checking that the gauge is Coulomb ........\n");
    
    
    y=5;
    for(t=0; t<NT; t++)
      for(z=0; z<NZ; z++)
	for(x=0; x<NX; x++)
	  {
	    Matrix B, D;
	    D.Dagger(L[IND(x,y,z,t-1,T_LINK)]);
	    B  = L[IND(x,y,z,t,T_LINK)] + D;
	    D.Dagger(L[IND(x-1,y,z,t,X_LINK)]);
	    B += L[IND(x,y,z,t,X_LINK)] + D;
	    D.Dagger(L[IND(x,y,z-1,t,Z_LINK)]);
	    B += L[IND(x,y,z,t,Z_LINK)] + D;
	    
	    p(B);
	    D.Dagger(B);
	    B -= D;
	    B -= B.Tr()/3;
	    
	    Float nrm = 0;
	    for(int j=0; j<3; j++)
	      for(int k=0; k<3; k++)
		nrm += norm(B(j,k));
	    VRB.Debug("is at (%d,%d,%d,%d) has norm= %e\n\n",x,y,z,t,nrm);
	  }
    return 0;
  }
}
q************/


#define MBIG 1000000000L
#define MSEED 1618033L
#define MZ 0L
#define FAC 1.0e-9			/* 	1.0/MBIG 	*/

static long ma[55] ;
static long mk,mj  ;
static int i,k,ii ;
static int inext,inextp ;
static int iff = 0 ;

Float ran345(long *idum)
{
  Float temp;
	if ( *idum < 0 || iff == 0 )
	{
		iff = 1 ;
		if ( *idum < 0 ) *idum = -(*idum) ;
		mj = MSEED - *idum ;
		mj = mj%MBIG ;
		ma[54] = mj ;
		mk = 1L ;
		for( i = 1; i < 55; i++)
		{
			ii = (21*i)%55 - 1 ;
			ma[ii] = mk ;
			mk = mj - mk ;
			if ( mk < MZ )
				mk += MBIG ;
			mj = ma[ii] ;
		}
		for( k = 0; k < 4; k++)
		{
			for( i = 0; i < 55; i++)
			{
				ma[i] -= ma[(i+31)%55] ;
				if ( ma[i] < MZ )
					ma[i] = ma[i] + MBIG ;
			}
		}
		inext = -1 ;
		inextp = 30 ;
		*idum = 1L ;
	}
	inext = inext + 1 ;
	if ( inext == 55 ) inext = 0 ;
	inextp = inextp + 1 ;
	if ( inextp == 55 ) inextp = 0 ;
	mj = ma[inext] - ma[inextp] ;
	if ( mj < MZ ) mj += MBIG ;
	ma[inext] = mj ;
	temp = (Float )mj * FAC ;
	return temp;
}



