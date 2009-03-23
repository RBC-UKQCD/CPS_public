#include <config.h>
#include <stdio.h>
#include <util/gjp.h>
#include <util/site.h>
#include <util/qcdio.h>
#include <alg/alg_tcharge.h>
#include <comms/glb.h>
CPS_START_NAMESPACE

//----------------------------------------------------------
//
// alg_tcharge.C
// 
// measures a simple defintion of the
// topological charge
//  
//  Using Clover Leaf   c.f   hep-lat/010610  Eq (6) -- (10)
//-----------------------------------------------------------


/*!
  takes the imaginary part of a matrix
*/
void ZeroReal(Matrix& m)
{
  Matrix temp;
  temp.Dagger( m );
  m-=temp;
  m*=0.5;
}

/*!
  constructs 

  \f[
  \frac{1}{32 \pi^2} \epsilon_{\mu \nu \eta \nu } {\rm Tr} \left\{ F_{\mu \nu} F_{\eta \nu} \right\}
  \f]

  passed two arrays of matrices holding the approximations for F:

  plaqs[0]  = F_01
  plaqs[1]  = F_02
  plaqs[2]  = F_03
  plaqs[3]  = F_12
  plaqs[4]  = F_13
  plaqs[5]  = F_23
  
*/

Complex MkTop( Matrix plaqs1[], Matrix plaqs2[] )
{
  const Float nfactor(-1.0/( 4 * 3.141592654 * 3.141592654 ));

  Matrix Top;

  // negative weight

  Top.DotMEqual( plaqs1[1] , plaqs2[4] );
  
  Top *= -1.0 ;
  // positive weight 
  
  Top.DotMPlus ( plaqs1[2] , plaqs2[3] );
  Top.DotMPlus ( plaqs1[5] , plaqs2[0] );
  
  return nfactor*Top.Tr();
}



/*!
  Calculate Clover leaf (1x1 size) SU(3) Matrix 
  Sheikholeslami, Wohlert, NPB259, 572 (1985)  Eq. (2.9)
*/
void CloverLeaf(Lattice& lattice, Matrix& pl,  int* pos, int mu, int nu)
{
   Matrix P0,P1,P2,P3;

   P0.ZeroMatrix();
   P1.ZeroMatrix();
   P2.ZeroMatrix();
   P3.ZeroMatrix();

   // each direction could be {0,1,2,3,4,5,6,7} coresponding to
   // the directions {n_x, n_y, n_z, n_t, -n_x, -n_y, -n_z, -n_t}


   int dirs0[4]={mu,nu, mu+4, nu+4};
   lattice.PathOrdProdPlus(P0, pos, dirs0, 4);


   int dirs1[4]={nu+4,mu+4, nu, mu};
   lattice.PathOrdProdPlus(P1, pos, dirs1, 4);

   int dirs2[4]={nu,mu+4, nu+4, mu};
   lattice.PathOrdProdPlus(P2, pos, dirs2, 4);

   int dirs3[4]={mu,nu+4, mu+4, nu};
   lattice.PathOrdProdPlus(P3, pos, dirs3, 4);


   
   P0 -=  P1;
   P0 +=  P2;
   P0 -=  P3;
   P0 *= 0.25;
   
   moveMem((Float*) &pl,(Float*) &P0, 18 * sizeof(Float) );

}


// Calculate Clover leaf (2x1, 1x2 size) SU(3) Matrix 
// Sheikholeslami, Wohlert, NPB259, 572 (1985)  Eq. (2.9)
// hep-lat/010610  Eq (8)
void CloverLeafRect(Lattice& lattice, Matrix& pl,  int* pos, int mu, int nu)
{
   Matrix P0,P1,P2,P3;


   // 1x2 size
   P0.ZeroMatrix();
   P1.ZeroMatrix();
   P2.ZeroMatrix();
   P3.ZeroMatrix();

   // each direction could be {0,1,2,3,4,5,6,7} coresponding to
   // the directions {n_x, n_y, n_z, n_t, -n_x, -n_y, -n_z, -n_t}


   {
     int dirs0[6]={mu,mu, nu, mu+4,mu+4, nu+4};
     lattice.PathOrdProdPlus(P0, pos, dirs0, 6);
     

     int dirs1[6]={nu+4,mu+4,mu+4, nu, mu,mu};
     lattice.PathOrdProdPlus(P1, pos, dirs1, 6);

     int dirs2[6]={nu,mu+4,mu+4, nu+4, mu,mu};
     lattice.PathOrdProdPlus(P2, pos, dirs2, 6);

     int dirs3[6]={mu,mu, nu+4, mu+4,mu+4, nu};
     lattice.PathOrdProdPlus(P3, pos, dirs3, 6);
   }

   P0 -=  P1;
   P0 +=  P2;
   P0 -=  P3;
   

   moveMem((Float*) &pl,(Float*) &P0, 18 * sizeof(Float) );
   
   // 2x1  size
   P0.ZeroMatrix();
   P1.ZeroMatrix();
   P2.ZeroMatrix();
   P3.ZeroMatrix();

   // each direction could be {0,1,2,3,4,5,6,7} coresponding to
   // the directions {n_x, n_y, n_z, n_t, -n_x, -n_y, -n_z, -n_t}


   {
     int dirs0[6]={mu,nu,nu, mu+4, nu+4,nu+4};
     lattice.PathOrdProdPlus(P0, pos, dirs0, 6);
     

     int dirs1[6]={nu+4,nu+4,mu+4, nu,nu, mu};
     lattice.PathOrdProdPlus(P1, pos, dirs1, 6);

     int dirs2[6]={nu,nu,mu+4, nu+4,nu+4, mu};
     lattice.PathOrdProdPlus(P2, pos, dirs2, 6);
     
     int dirs3[6]={mu,nu+4,nu+4, mu+4, nu, nu};
     lattice.PathOrdProdPlus(P3, pos, dirs3, 6);
   }

   P0 -=  P1;
   P0 +=  P2;
   P0 -=  P3;

   pl += P0;
   pl *= 1.0/16.0;
}

void CloverLeaf1x3(Lattice& lattice, Matrix& pl,  int* pos, int mu, int nu)
{
   Matrix P0,P1,P2,P3;


   // 1x3 size
   P0.ZeroMatrix();
   P1.ZeroMatrix();
   P2.ZeroMatrix();
   P3.ZeroMatrix();

   // each direction could be {0,1,2,3,4,5,6,7} coresponding to
   // the directions {n_x, n_y, n_z, n_t, -n_x, -n_y, -n_z, -n_t}


   {
     int dirs0[8]={mu,mu, mu,
		   nu, 
		   mu+4,mu+4,mu+4,
		   nu+4};
     lattice.PathOrdProdPlus(P0, pos, dirs0, 8);
     

     int dirs1[8]={nu+4,
		   mu+4,mu+4, mu+4,
		   nu, 
		   mu,mu,mu };
     lattice.PathOrdProdPlus(P1, pos, dirs1, 8);

     int dirs2[8]={nu,
		   mu+4,mu+4,mu+4,
		   nu+4, 
		   mu,mu,mu};
     lattice.PathOrdProdPlus(P2, pos, dirs2, 8);

     int dirs3[8]={mu,mu,mu,
		   nu+4, 
		   mu+4,mu+4, mu+4,
		   nu};
     lattice.PathOrdProdPlus(P3, pos, dirs3, 8);
   }

   P0 -=  P1;
   P0 +=  P2;
   P0 -=  P3;
   

   moveMem((Float*) &pl,(Float*) &P0, 18 * sizeof(Float) );
   
   // 3x1  size
   P0.ZeroMatrix();
   P1.ZeroMatrix();
   P2.ZeroMatrix();
   P3.ZeroMatrix();

   // each direction could be {0,1,2,3,4,5,6,7} coresponding to
   // the directions {n_x, n_y, n_z, n_t, -n_x, -n_y, -n_z, -n_t}
   {
     int dirs0[8]={mu,
		   nu,nu,nu,
		   mu+4, 
		   nu+4,nu+4,nu+4};
     lattice.PathOrdProdPlus(P0, pos, dirs0, 8);
     
     int dirs1[8]={nu+4,nu+4,nu+4,
		   mu+4, 
		   nu,nu,nu,
		   mu};
     lattice.PathOrdProdPlus(P1, pos, dirs1, 8);

     int dirs2[8]={nu,nu,nu,
		   mu+4, 
		   nu+4,nu+4,nu+4,
		   mu};
     lattice.PathOrdProdPlus(P2, pos, dirs2, 8);
     
     int dirs3[8]={mu,
		   nu+4,nu+4,nu+4,
		   mu+4, 
		   nu, nu, nu};
     lattice.PathOrdProdPlus(P3, pos, dirs3, 8);
   }

   P0 -=  P1;
   P0 +=  P2;
   P0 -=  P3;

   pl += P0;
   pl *= 1.0/24.0;
}

void CloverLeaf2x2(Lattice& lattice, Matrix& pl,  int* pos, int mu, int nu)
{
  Matrix P0,P1,P2,P3;
  // 1x2 size
  P0.ZeroMatrix();
  P1.ZeroMatrix();
  P2.ZeroMatrix();
  P3.ZeroMatrix();
  
  // each direction could be {0,1,2,3,4,5,6,7} coresponding to
  // the directions {n_x, n_y, n_z, n_t, -n_x, -n_y, -n_z, -n_t}
  int dirs0[8]={mu,mu, nu, nu, mu+4,mu+4, nu+4, nu+4};
  lattice.PathOrdProdPlus(P0, pos, dirs0, 8);
    
  int dirs1[8]={nu+4, nu+4, mu+4, mu+4, nu, nu, mu,mu };
  lattice.PathOrdProdPlus(P1, pos, dirs1, 8);
  
  int dirs2[8]={nu, nu, mu+4, mu+4, nu+4, nu+4, mu, mu };
  lattice.PathOrdProdPlus(P2, pos, dirs2, 8);
  
  int dirs3[8]={mu,mu, nu+4, nu+4, mu+4, mu+4, nu, nu };
  lattice.PathOrdProdPlus(P3, pos, dirs3, 8);
  
  P0 -=  P1;
  P0 +=  P2;
  P0 -=  P3;
  P0 *= 1.0/16;
   
  moveMem((Float*) &pl,(Float*) &P0, 18 * sizeof(Float) );
  
}

void CloverLeaf3x3(Lattice& lattice, Matrix& pl,  int* pos, int mu, int nu)
{
  Matrix P0,P1,P2,P3;
  // 1x2 size
  P0.ZeroMatrix();
  P1.ZeroMatrix();
  P2.ZeroMatrix();
  P3.ZeroMatrix();
  
  // each direction could be {0,1,2,3,4,5,6,7} coresponding to
  // the directions {n_x, n_y, n_z, n_t, -n_x, -n_y, -n_z, -n_t}
  int dirs0[12]={ mu, mu, mu,
		  nu, nu, nu,
		  mu+4, mu+4, mu+4,
		  nu+4, nu+4, nu+4 };
  lattice.PathOrdProdPlus(P0, pos, dirs0, 12);
    
  int dirs1[12]={nu+4, nu+4, nu+4,
		 mu+4, mu+4, mu+4,
		 nu,   nu,   nu, 
		 mu,   mu,   mu   };
  lattice.PathOrdProdPlus(P1, pos, dirs1, 12);
  
  int dirs2[12]={nu, nu, nu,
		 mu+4, mu+4, mu+4,
		 nu+4, nu+4, nu+4,
		 mu,   mu  , mu   };
  lattice.PathOrdProdPlus(P2, pos, dirs2, 12);
  
  int dirs3[12]={mu  , mu,   mu,
		 nu+4, nu+4, nu+4, 
		 mu+4, mu+4, mu+4, 
		 nu  , nu  , nu };
  lattice.PathOrdProdPlus(P3, pos, dirs3, 12);
  
  P0 -=  P1;
  P0 +=  P2;
  P0 -=  P3;
  P0 *= 1./(9*4);
   
  moveMem((Float*) &pl,(Float*) &P0, 18 * sizeof(Float) );
}


// number of clover-leaf functions
const int nfunc(5);

typedef void (*leaf_function)(Lattice&, Matrix&,  int*, int, int);

// map to the functions used
leaf_function leaf_map[5] = { &CloverLeaf,
			      &CloverLeafRect,
			      &CloverLeaf2x2,
			      &CloverLeaf3x3,
			      &CloverLeaf1x3 };
// map to the names
const char* names[5] = { "1x1",
			 "1x2",
			 "2x2",
			 "3x3",
			 "1x3" };


void AlgTcharge::run()
{
  Lattice& lattice( AlgLattice() );  

  Float tmat[nfunc][nfunc];
  for (int f1(0);f1<nfunc;f1++)
    for (int f2(0);f2<nfunc;f2++)
      tmat[f1][f2] = 0;

  // sum over lattice
  Site nloop;
  
  while ( nloop.LoopsOverNode() )
    {
      // Array of imaginary parts of the plaquettes
      // at a given site
      
      // plaqs[0]  = F_01
      // plaqs[1]  = F_02
      // plaqs[2]  = F_03
      // plaqs[3]  = F_12
      // plaqs[4]  = F_13
      // plaqs[5]  = F_23

      Matrix plaqs[nfunc][6];
      //
      // fill plaqs with the full plaquettes
      // - then zero the real parts
      //

      int mu;
      int nu;
      int index(0);

      for (mu=0;mu<3;++mu)
        {
          for (nu=mu+1;nu<4;nu++)
            { 
	      for (int f(0);f<nfunc;f++)
		{
		  (*(leaf_map[f]))( lattice, plaqs[f][index], nloop.pos(), mu, nu );
		  ZeroReal(plaqs[f][index]);
		}
	      index++;
            }
        }
      for (int f1(0);f1<nfunc;f1++)
	{
	  for (int f2(f1);f2<nfunc;f2++)
	    {
	      tmat[f1][f2] += MkTop(plaqs[f1],plaqs[f2]).real();
	    }
	}
    }
  
  // global sum the approximations
  for (int f1(0);f1<nfunc;f1++)
    {
      for (int f2(f1);f2<nfunc;f2++)
	{
	  glb_sum( &tmat[f1][f2] );
	}
    }
  
  // Print out results
  //----------------------------------------------------------------

  if(common_arg->filename != 0)
    {
      char *fname = "alg_tcharge()";
      FILE *fp;
      if( (fp = Fopen(common_arg->filename, "a")) == NULL ) {
        ERR.FileA(cname,fname,common_arg->filename);
      }
      Fprintf(fp,"AlgTcharge:\n");
      Fprintf(fp,"nleaf : %i\n",nfunc);
      for (int f(0);f<nfunc;f++)
	Fprintf(fp,"   %i : %s\n",f,names[f]);
      for (int f1(0);f1<nfunc;f1++)
	{
	  for (int f2(f1);f2<nfunc;f2++)
	    {
	      Fprintf(fp,"%i %i : %15e\n",f1,f2,tmat[f1][f2]);
	    }
	}
      Fclose(fp);
    }

}
CPS_END_NAMESPACE

