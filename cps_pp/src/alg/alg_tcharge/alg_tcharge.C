#include <stdio.h>
#include <util/gjp.h>
#include <util/site.h>
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


void AlgTcharge::run()
{
  Lattice& lattice( AlgLattice() );  

  Float tcharge_node     (0);
  Float tchargeCross_node(0);
  Float tchargeRect_node (0);
  
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

      Matrix plaqs1[6];
      Matrix plaqs2[6];
      
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
	      CloverLeaf    ( lattice, plaqs1[index], nloop.pos(), mu, nu );
              ZeroReal(plaqs1[index]);
              
              CloverLeafRect( lattice, plaqs2[index], nloop.pos(), mu, nu );
              ZeroReal(plaqs2[index]);
                            
              index++;
            }
        }

      /*
        construct the matrix that must be traced to
        give the topological charge
      */
      const Float clover( MkTop(plaqs1,plaqs1).real() );
      const Float rect  ( MkTop(plaqs2,plaqs2).real() );
      const Float cross ( MkTop(plaqs1,plaqs2).real() );
      
        
      /*
        keep running sum for node
      */
      tcharge_node     += clover;
      tchargeCross_node+= cross ;
      tchargeRect_node += rect  ;

    }
  
  glb_sum( &tcharge_node      ) ;
  glb_sum( &tchargeCross_node ) ;
  glb_sum( &tchargeRect_node  ) ;

  charge_clov = tcharge_node;
  charge_rect = tchargeRect_node;
  
  // construct O(a^2) improved definition of Q by cancelling the
  // o(a^2) errors of Q^R and Q^C
  const Float c0(+5.0/3.0);
  const Float c1(-2.0/3.0);
  charge = c0*charge_clov + c1*charge_rect;
  
  // construct O(a^2) improved of Q from O(a^2) definition
  // of F.
  charge_2 = 
      (25.0/9.0) * tcharge_node 
    + (4.0 /9.0) * tchargeRect_node 
    - (20.0/9.0) * tchargeCross_node; 

  // Print out results
  //----------------------------------------------------------------

  if(common_arg->results != 0)
    {
      char *fname = "alg_tcharge()";
      char *cname = "alg_tcharge()";
      FILE *fp;
      if( (fp = fopen((char *)common_arg->results, "a")) == NULL ) {
        ERR.FileA(cname,fname, (char *)common_arg->results);
      }
      fprintf(fp,"AlgTcharge   :\n Clov      : %e\n", tcharge_node);
      fprintf(fp," Rect      : %e\n", tchargeRect_node );
      fprintf(fp," Cross     : %e\n", tchargeCross_node);
      fprintf(fp," Tot(a^2) Q: %e \n",charge    );
      fprintf(fp," Tot(a^2) F: %e \n",charge_2 );

      fclose(fp);
    }

}
CPS_END_NAMESPACE

