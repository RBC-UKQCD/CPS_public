/*!\file
  \brief Implementation of AlgRandomGauge and AlgRotateGauge methods.

  $Id: alg_rnd_gauge.C,v 1.3 2004-09-02 16:52:48 zs Exp $
*/

#include <math.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/smalloc.h>
#include <util/vector.h>
#include <util/verbose.h>
#include <util/error.h>
#include <comms/scu.h>
#include <alg/alg_rnd_gauge.h>


CPS_START_NAMESPACE

//====================================
// change them to produce a paramter
//====================================



int AlgRandomGauge::index( int x, int y, int z, int t )
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
    const char *fname = "AlgRandomGauge";
    VRB.Func(cname,fname);

    NX = GJP.XnodeSites();
    NY = GJP.YnodeSites();
    NZ = GJP.ZnodeSites();
    NT = GJP.TnodeSites();
  
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
  

  LRG.SetInterval(0,3); 

  for( int t(0); t<NT ; t++ )
    {
      for( int z(0); z<NZ ; z++ )
	{
	  for( int y(0); y<NY ; y++ )
	    {
	      for( int x(0); x<NX ; x++ ) 
		{
		  //---------------------------------------------
		  // Generate random gauge transformation
		  //---------------------------------------------
		  
		  Matrix& gmat(lat.FixGaugePtr()[0][ index(x,y,z,t) ]);

		  LRG.AssignGenerator(x, y, z, t);
		  
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

void AlgRandomGauge::RandMatrix(  Float theta , Matrix& U , int sub )
{

    const int su2_index[][3]= { {0,1,2},
				{0,2,1},
				{1,2,0} };
    
  //=========================
  // calculate random three dimensional unit vector
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

void AlgRandomGauge::printMatrix( const Matrix& x )
{
    for (int i=0;i<3;++i)
	for (int j=0;j<3;++j)
	    printf("%d %d : %e %e \n",i,j,real(x(i,j)),imag(x(i,j)));


}




int AlgRotateGauge::index(int x, int y, int z, int t, int l){
    return 
    (((((t+NT)%NT)*NZ+((z+NZ)%NZ))*NY+   
       ((y+NY)%NY))*NX+((x+NX)%NX))*4+l;
}
	


//------------------------------------------------------------------
// Fixes the gauge of the gauge fields using preset gauge fixing
// matrices
//------------------------------------------------------------------
/*!
  \pre Gauge transformation matrices must be present inthe Lattice object.
*/
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

    NX = GJP.XnodeSites();
    NY = GJP.YnodeSites();
    NZ = GJP.ZnodeSites();
    NT = GJP.TnodeSites();

    int node_sites[4]={ NX , NY , NZ ,NT };

  //------------------------------
  // holds the offset to "flatten"  
  // an array indexed by x,y,z and t
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

		      mtmp.DotMEqual( gl, L[index(x,y,z,t,nu)] );
		      
		      M[index(x,y,z,t,nu)].DotMEqual( mtmp, gr );

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
			   const Matrix *m_offset, 
			   const int *x, 
			   int nu,
			   const int *node_sites, 
			   const int *m_dir_offset ) 
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


/*!
  Fill one of three su(2) subgroups with

  cos(theta) -i nz sin(theta)  , -(inx+ny)sin(theta)
  (-inx+ny)sin(theta)          , cos(theta) + i nz sin (theta) ..,

  where n* are the components of a 3d unit vector, and theta
  is a specified argument.
  
  other co-ord 1 on diagonal, zero otherwise.
*/



CPS_END_NAMESPACE
