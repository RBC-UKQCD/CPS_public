#include <stdio.h>
#include <math.h>
#include <alg/common_arg.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/smalloc.h>
#include <util/vector.h>
#include <util/verbose.h>
#include <util/error.h>
#include <comms/scu.h>
#include <alg/alg_rnd_gauge.h>
#include <alg/no_arg.h>

CPS_START_NAMESPACE

#define X_LINK 0
#define Y_LINK 1
#define Z_LINK 2
#define T_LINK 3

//====================================
// change them to produce a paramter
//====================================

#define NX (GJP.XnodeSites())
#define NY (GJP.YnodeSites())
#define NZ (GJP.ZnodeSites())
#define NT (GJP.TnodeSites())
#define XSHIFT (GJP.XnodeCoor()*GJP.XnodeSites())
#define YSHIFT (GJP.YnodeCoor()*GJP.YnodeSites())
#define ZSHIFT (GJP.ZnodeCoor()*GJP.ZnodeSites())
#define TSHIFT (GJP.TnodeCoor()*GJP.TnodeSites())

#define IND(x,y,z,t,l) ((((((t+NT)%NT)*NZ+((z+NZ)%NZ))*NY+   \
			  ((y+NY)%NY))*NX+((x+NX)%NX))*4+l)

inline int index( int x, int y, int z, int t )
{
  return (
	  (
	   ( 
	    (t+NT)%NT)*NZ +((z+NZ)%NZ) 
	   )*NY+((y+NY)%NY)
	  )*NX
	    +((x+NX)%NX);
}



AlgRandomGauge::AlgRandomGauge( Lattice& latt , CommonArg *c_arg )
  :Alg(latt,c_arg),
   theta(1)
{
  cname       = "AlgRandomGauge";
  char *fname = "AlgRandomGauge(Lattice&)";
  VRB.Func(cname,fname);
}




void AlgRandomGauge::run()
{
  //-------------------------
  // Set the Lattice pointer
  //-------------------------

  Lattice& lat( AlgLattice() );
    
  //----------------------------------------------
  //
  // FixGaugePtr return Lattice::fix_gauge_ptr
  // which is of type Matrix **;
  // The memory for this is usually allocated by 
  //
  // Lattice::FixGaugeAllocate(....)
  //
  // The top level array is of size GaugePtrLen
  // 
  //   - this depends on the gauge being fixed
  //
  // while the next level is of size;
  //
  //    MemBlockLen * sizeof(Matrix)      
  //
  // where MemBlockLen is a function of the 
  // number of hyperplanes that need to be
  // fixed ( it's *just* a volume factor)
  //
  // and deallocated by
  //
  //           Lattice::FixGaugeFree();
  //
  //-------------------------------------------------  
  //
  // For our purposes one gauge transformation matrix 
  // for every lattice site is needed. This is the same as is
  // used in Landau gauge fixing. 
  //
  //-------------------------------------------------

  lat.FixGaugeAllocate( FIX_GAUGE_LANDAU, 0, 0 );

  //-------------------------------------------------
  // The random gauge transformation matrix
  //-------------------------------------------------
  
  const int nt(NT);
  const int nx(NX);
  const int ny(NY);
  const int nz(NZ);


  LRG.SetInterval(0,3);

  for( int t(0); t<nt ; t++ )
    {
      for( int z(0); z<nz ; z++ )
	{
	  for( int y(0); y<ny ; y++ )
	    {
	      for( int x(0); x<nx ; x++ ) 
		{
		  //---------------------------------------------
		  // Generate random gauge transformation
		  //---------------------------------------------
		  
		  Matrix& gmat(lat.FixGaugePtr()[0][ index(x,y,z,t) ]);
		  
		  const IFloat dum( LRG.Urand()  ); // between 0 and 3
                  
                  if      ( dum >= 2 )
                    UpperRandMatrix( theta, gmat );
                  else if ( dum >=1  )
                    LowerRandMatrix( theta, gmat );
                  else
                    MixedRandMatrix( theta, gmat );


		  //-----------------------------------------
		  // just to be paranoid: make sure that the 
                  // gauge transformation is unitary
		  //-----------------------------------------
		  gmat.Unitarize();

		} 
	    }
	}
    } //end loop over sites
}




//------------------------------------------------------------------
// Fixes the gauge of the gauge fields using preset gauge fixing
// matrices
//------------------------------------------------------------------

void AlgRotateGauge::run()
{
  char *fname = "run()";
  VRB.Func(cname,fname);
  
  //-------------------------
  // Set the Lattice pointer
  //-------------------------

  Lattice& lat( AlgLattice() );
  

  //---------------------------------------
  // check gauge fixing matrice are
  // allocated with one fixing matrix per
  // lattice site landau gauge 
  //---------------------------------------

  if ( lat.FixGaugeKind() != FIX_GAUGE_LANDAU )
    {
      ERR.General(cname,fname,"Gauge fixing matrices not allocated correctly\n");
    }

  //----------------------------
  // pointer to the gauge field
  //----------------------------

  Matrix *L( lat.GaugeField() );

  //---------------------
  // get node parameters
  //---------------------

  int node_sites[4]={ NX , NY , NZ ,NT };

  //------------------------------
  // holds the offset to "flatten"  
  // an array indexed by x,y,z
  // and t
  //=============================
  // m_dir_offset[0]=1;
  // m_dir_offset[1]=NX;
  // m_dir_offset[2]=NX*NY;
  // m_dir_offset[3]=NX*NY*NZ;
  //=============================

  int  m_dir_offset[4];
  m_dir_offset[0]   = 1  ;

  for(int mu=1; mu<4; mu++)
    {
      m_dir_offset[mu] = m_dir_offset[mu-1] * node_sites[mu-1];
    }
  
  //-----------------------------------------
  // temp matrices for rotated gauge fields
  //-----------------------------------------
  
  Matrix *M;
  M = (Matrix *) smalloc(4*NX*NY*NZ*NT*sizeof(Matrix));
    
  Matrix   gl; // "left-hand side"  gauge rotation 
  Matrix   gr; // "right-hand side" gauge rotation
  Matrix mtmp; // temp matrix
  
  int st[4]       ;
  int offset_st   ;
  Matrix *m_offset;

  //---------------------------------------------
  // Fix Lattice
  //---------------------------------------------
  
  for( int t=0; t<NT; t++)
    {
      for(int z=0; z<NZ; z++)
	{
	  for(int y=0; y<NY; y++)
	    {
	      for(int x=0; x<NX; x++) 
		{
		  
		  st[0]=x; 
		  st[1]=y; 
		  st[2]=z; 
		  st[3]=t;
		  
		  offset_st = st[0]*m_dir_offset[0] 
		    + st[1]*m_dir_offset[1]
		      + st[2]*m_dir_offset[2] 
			+ st[3]*m_dir_offset[3];
		  
		  //===============================================
		  // offset_st = x + y*NX + z*NX*NY + t*NX*NY*NZ
		  //===============================================
		  // this is the flattened local space co-ordinate
		  // (the same as index(x,y,z,t) )
		  //===============================================


		  //=======================================================
		  // this is the local position of the gauge tranformation
		  // matrix
		  //=======================================================
		  
		  m_offset = lat.FixGaugePtr()[0] + offset_st;
		  
		  for(int nu=0; nu<4; nu ++)
		    {
		      
		      gl = *m_offset;
		      
		      //----------------------------------
		      // get the next gauge tranformation
		      // in the direction nu
		      //----------------------------------

		      mtmp = GetMat(  m_offset, 
                                      st, 
                                      nu, 
                                      node_sites, 
                                      m_dir_offset );
		      
		      gr.Dagger( mtmp );
		      
		      //============================
		      // gauge transform the link
		      //============================

		      mtmp.DotMEqual( gl, L[IND(x,y,z,t,nu)] );
		      
		      M[IND(x,y,z,t,nu)].DotMEqual( mtmp, gr );

		    } // end loop over directions
		} // x
	    }// y 
	}// z
    } // t

  lat.GaugeField(M);
  lat.ClearAllBufferedLink();


  sfree(M);
}

//=============================================================
// gets the gauge transformation matrix from the neighbouring
// link. Knows about node boundaries etc..
//=============================================================
Matrix AlgRotateGauge::GetMat( 
			   Matrix *m_offset, 
			   const int *x, 
			   int nu,
			   int *node_sites, 
			   int *m_dir_offset ) 
{
  static Matrix m_tmp1;

  if( x[nu] == node_sites[nu]-1) 
    {
      getPlusData((IFloat *)&m_tmp1,
		  (IFloat *)(m_offset-x[nu]*m_dir_offset[nu]), 
		  18, 
		  nu );

      return m_tmp1;
    } 
  else 
    {
      return *(m_offset+m_dir_offset[nu]);
    }
  
}


//----------------------------------------------------------------
// Construct an SU(3) matrix by sticking an SU(2)
// matrix in the upper/lower/mixed corner 
//----------------------------------------------------------------
//
// agrument is a Float "theta" and the matrix generated 
// has the su2 matrix
//
//  cos(theta) -i nz sin(theta)  , -(inx+ny)sin(theta)
//  (-inx+ny)sin(theta)          , cos(theta) + i nz sin (theta) ..,
//
// where n* are the components of a 3d unit vector,
//
//----------------------------------------------------------------


const int su2_index[][3]= { {0,1,2},
                            {0,2,1},
                            {1,2,0} };


void RandMatrix(  Float theta , Matrix& U , int sub )
{
  //=========================
  // calculate random
  // three dimensional
  // unit vector
  //=========================
  
  Float nx( LRG.Urand() );
  Float ny( LRG.Urand() );
  Float nz( LRG.Urand() );
  
  const Float norm(1.0/sqrt(nx*nx+ny*ny+nz*nz));
  
  nx*=norm;
  ny*=norm;
  nz*=norm;

  //==========================
  // set up the elements
  //==========================

  const Float Cos( cos(theta) );
  const Float Sin( sin(theta) );
  const Complex el00( Cos     , -nz*Sin );
  const Complex el01( -ny*Sin , -nx*Sin );
  const Complex el10(  ny*Sin , -nx*Sin );
  const Complex el11( Cos     ,  nz*Sin );
  
  //=================
  // Fill the Matrix
  //=================

  U.ZeroMatrix();
  U(su2_index[sub][0],su2_index[sub][0]) = el00;
  U(su2_index[sub][0],su2_index[sub][1]) = el01;
  U(su2_index[sub][1],su2_index[sub][0]) = el10;
  U(su2_index[sub][1],su2_index[sub][1]) = el11;
  U(su2_index[sub][2],su2_index[sub][2]) = Complex(1,0);
 
}


//==================================
// print out a matrix for debugging
//==================================

void printMatrix( const Matrix& x )
{
  for (int i=0;i<3;++i)
    {
      for (int j=0;j<3;++j)
	{
	  printf("%d %d : %e %e \n",i,j,real(x(i,j)),imag(x(i,j)));
	}
    }
}


CPS_END_NAMESPACE
