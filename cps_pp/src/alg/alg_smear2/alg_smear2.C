//--------------------------------------------------------------------
/*!\file
  \brief Implementation of smearing class methods.

  AlgSmear, AlgApeSmear, AlgKineticSmear and AlgHypSmear classes.
  
  $Id: alg_smear2.C,v 1.2 2006-05-30 20:32:27 chulwoo Exp $
*/
//--------------------------------------------------------------------
#include <config.h>
#include <math.h>
#include <util/qcdio.h>
#include <alg/alg_smear2.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/error.h>
#include <util/site.h>
#include <util/link_buffer.h>
#include <util/smalloc.h>

CPS_START_NAMESPACE

/*!
  get su2 submatrix of x and return the su3 matrix y that
  has the inverse of this matrix in the relevant row and column
*/
void sub( Matrix& x, Matrix& y , int ind )
{

  const int su2_index[][3]= { {0,1,2},
                              {0,2,1},
                              {1,2,0} };
  y=x;
  const int zero_rc(su2_index[ind][2]);
  const int i1     (su2_index[ind][0]);
  const int i2     (su2_index[ind][1]);
  // zero out the row and column not used
  int i;
  for (i=0;i<3;i++) { y(zero_rc,i)=Complex(0,0); }
  for (i=0;i<3;i++) { y(i,zero_rc)=Complex(0,0); }
  y(zero_rc,zero_rc)=Complex(1,0);
  // project onto SU(2)
  Float p0(x(i1,i1).real() + x(i2,i2).real());
  Float p1(x(i1,i2).imag() + x(i2,i1).imag());
  Float p2(x(i1,i2).real() - x(i2,i1).real());
  Float p3(x(i1,i1).imag() - x(i2,i2).imag());
  const Float psqr( sqrt(p0*p0 + p1*p1 + p2*p2 + p3*p3) ); 
  Float ipsqr;
  if ( psqr == 0. ) { ipsqr=1;     }
  else              { ipsqr=1/psqr;} 
  p0*=ipsqr; p1*=ipsqr; p2*=ipsqr; p3*=ipsqr;
  // fill with inverse
  y(i1,i1) = Complex( p0,-p3);
  y(i2,i2) = Complex( p0, p3);
  y(i1,i2) = Complex(-p2,-p1);
  y(i2,i1) = Complex( p2,-p1);
}

int su3_proj( Matrix& x , Float tolerance )
{
  // usually takes ~5 hits, so just exit
  // if hits the max, as something is 
  // probably very wrong.
  const int max_iter(10000);
  Matrix tmp  ;
  Matrix inv  ;
  Matrix y    ;
  Matrix ycopy;
  Matrix xdag ; xdag.Dagger(x);
  y.UnitMatrix();
  tmp = xdag;
  Float old_tr( xdag.ReTr() );
  int i,j;
  for (i=0;i<max_iter;i++)
    {
      // loop over su2 subgroups
      Float diff(-1);
      for (j=0;j<3;j++)
        {
          sub(tmp,inv,j);
          ycopy = y ; 
          y  .DotMEqual( inv, ycopy );
          tmp.DotMEqual( y, xdag );
          const Float tr (tmp.ReTr());
          const Float dtr(tr-old_tr );
          if ( dtr > diff ) { diff = dtr ; }
          old_tr = tr;
        } 
      // for single precision the difference seems
      // to never get below 1e-7 (not too suprising)
      if (diff < tolerance) { break; }
    }
  if ( i == max_iter )
    {
      // hit the max iterations
      ERR.General("","su3_proj","Max iterations");
    }
  x=y; 
  return i;
} 
 
/*!
  \param latt The Lattice object containg the gauge field with which smearing
  is done.
  \param c_arg Container for generic parameters. .
  \param su3_proj Whether or not to project the smeared link on to the SU(3)
  manifold.
*/
AlgSmear2::AlgSmear2( Lattice&   lat,
		      CommonArg* ca ,
		      int su3_proj ):
  Alg(lat,ca),
  bool_su3_proj(su3_proj),
  tolerance    (1e-6),
  orthog       (-1)
{
  cname = "AlgSmear2";
  lat_back = new Matrix[GJP.VolNodeSites()*4];
  if ( lat_back == 0x0 ) { ERR.Pointer(cname, cname,"lat_back"); }
}

AlgSmear2::~AlgSmear2()
{
  delete[] lat_back;
}


void AlgSmear2::run()
{
  Lattice& lattice(AlgLattice());
  
  //-----------------------------------------------------
  // Make a copy of the lattice. This will be the one
  // that is smeared. It will be swapped with the one in
  // the lattice class at the end of run()
  
  lattice.CopyGaugeField(lat_back);
  
  Site nloop;
  
  while ( nloop.LoopsOverNode() )
    {
      int mu;
      //      for (mu=0;mu<4;++mu)
      for (mu=0;mu<3;++mu)
        {
          // check that this isn't along a direction we are not
          // smearing in

          if ( mu != get_orthog() )
            {

              // link to be smeared
              
              Matrix& link(*(lat_back + 4*nloop.Index() + mu ));
              
              /*
                the matrix to be added to this link, take
                this from the gauge field the lattice
                class knows about. This allows the link
                buffer to be used and also means that smearing
                a link will not change the staple of the
                neighbouring links
              */
              smear_link2(link,nloop.pos(),mu);
              
              // project the matrix back down to SU(3) (if needed)
              
              if ( bool_su3_proj )  { su3_proj( link , tolerance ); }
            } // smear in this direction ?
        } // direction
    } // spatial position
  
  // copy smeared configuration to the lattice

  lattice.GaugeField(lat_back);
  lattice.ClearAllBufferedLink();

}

/*!
  \param latt The Lattice containing the gauge field.
  \param link The accumulated three staple.
  \param pos The coordinates of the link.
  \param u  The direction of the link.
  \param orth A direction in which no smearing is done. 
*/
void three_staple2( Lattice& latt,  Matrix& link , 
                   int *pos, int u, int orth )
{
  Matrix acumulate_mp; acumulate_mp.ZeroMatrix();
  int dir[3];
  //loop over all directions (+ve and -ve)

  for ( int v=0; v<8; v++ )
    {      
      // except the ones the link is aligned with 
      //      if((v&3)==u || (v&3) == orth)  continue; 
      if((v&3) == u)  continue; 
      if((v&3) == orth)  continue; 
      if((v&3) == 3)  continue; 
      
      const int v1((v+4)&7); // direction opposite to v
      
      dir[0] = v;
      dir[1] = u;
      dir[2] = v1;
      
      latt.PathOrdProdPlus(acumulate_mp, pos, dir, 3); 
    }
  // 18 is the matrix size
  moveMem((Float *) &link, (Float*)&acumulate_mp, 18*sizeof(Float));
}

AlgApeSmear2::AlgApeSmear2(Lattice&     lat,
            CommonArg*   ca ,
            ApeSmearArg* asa ):
  AlgSmear2(lat,ca,1),
  cname("AlgApeSmear2")
{
  c = asa->coef;
  set_tol   (asa->tolerance);
  set_orthog(asa->orthog);
}
  
/*!
  If an output file is specified in the CommonArg argument, then
  the smearing coefficients are written to the file.

*/
void AlgApeSmear2::run()
{
  if(common_arg->filename != 0){
    FILE* f = Fopen(common_arg->filename, "a");
    if(!f) ERR.FileA(cname, "run", common_arg->filename);
    Fprintf(f,"AlgApeSmear2: coef = %e \n",c);
    Fclose(f);
  }
  AlgSmear2::run();
}



void AlgApeSmear2::smear_link2(Matrix& link,
                             int*     pos,
                             int       mu)
{
  Lattice& lattice(AlgLattice());
  Matrix stap;
  three_staple2(lattice,stap,pos,mu,get_orthog());
  link*=(1.0-c);
  stap*=c/6.0;
  link+=stap;
}

CPS_END_NAMESPACE
