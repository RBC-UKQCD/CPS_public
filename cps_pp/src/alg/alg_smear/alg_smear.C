//--------------------------------------------------------------------
/*!\file
  \brief Implementation of smearing class methods.

  AlgSmear, AlgApeSmear, AlgKineticSmear and AlgHypSmear classes.
  
  $Id: alg_smear.C,v 1.7 2008-02-12 18:16:30 chulwoo Exp $
*/
//--------------------------------------------------------------------
#include <config.h>
#include <math.h>
#include <util/qcdio.h>
#include <alg/alg_smear.h>
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
AlgSmear::AlgSmear( Lattice&   lat,
                    CommonArg* ca ,
                    int su3_proj ):
  Alg(lat,ca),
  bool_su3_proj(su3_proj),
  tolerance    (1e-6),
  orthog       (-1)
{
  cname = "AlgSmear";
  lat_back = new Matrix[GJP.VolNodeSites()*4];
  if ( lat_back == 0x0 ) { ERR.Pointer(cname, cname,"lat_back"); }
}

AlgSmear::~AlgSmear()
{
  delete[] lat_back;
}


void AlgSmear::run()
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
      for (mu=0;mu<4;++mu)
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
              smear_link(link,nloop.pos(),mu);
              
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
void three_staple( Lattice& latt,  Matrix& link , 
                   int *pos, int u, int orth )
{
  Matrix acumulate_mp; acumulate_mp.ZeroMatrix();
  int dir[3];
  //loop over all directions (+ve and -ve)

  for ( int v=0; v<8; v++ )
    {      
      // except the ones the link is aligned with 
      if((v&3)==u || (v&3) == orth )  continue; 
      
      const int v1((v+4)&7); // direction opposite to v
      
      dir[0] = v;
      dir[1] = u;
      dir[2] = v1;
      
      latt.PathOrdProdPlus(acumulate_mp, pos, dir, 3); 
    }
  // 18 is the matrix size
  moveMem((Float *) &link, (Float*)&acumulate_mp, 18*sizeof(Float));
}



void five_staple ( Lattice& latt,  Matrix& link ,
                   int *pos, int u , int orth )
{
  Matrix acumulate_mp; acumulate_mp.ZeroMatrix();

  int v;
  int dir[5];
  // loop over all directions (+ve and -ve)
  for(v=0; v<8; v++)
    {
      // except the ones the link is aligned with 
      if((v&3)==u || (v&3) == orth ) { continue; }
      
      const int v1((v+4)&7); // direction opposite to v
      
      // loop over all directions (+ve and -ve)
      for(int w=0; w<8; w++)
        {
          /*
            except the ones aligned with either
            the link or the v direction
          */
          if( (w&3) == u || (w&3) == (v&3) || (w&3) == orth ) { continue; }
        
          const int w1((w+4)&7); // direction opposite to w
          
          // the chair (but not by Guofengs definition).
          dir[0] = v; 
          dir[1] = w;
          dir[2] = u;
          dir[3] = w1;
          dir[4] = v1;
          
          latt.PathOrdProdPlus(acumulate_mp, pos, dir, 5); 
        }
    }// all directions 

  // 18 is the matrix size
  moveMem((Float *) &link, (Float*)&acumulate_mp, 18*sizeof(Float));
}


void seven_staple( Lattice& latt,  Matrix& link ,
                   int *pos, int u , int orth )
{
  Matrix acumulate_mp; acumulate_mp.ZeroMatrix();
  int v;
  int w;
  int x;
  int dir[7];
  // loop over all directions (+ve and -ve)
  for(v=0; v<8; v++)
    {
      // except the ones the link is aligned with 
      if((v&3)==u || (v&3)==orth ) { continue; }
      
      const int v1((v+4)&7); // direction opposite to v
      
      // loop over all directions (+ve and -ve)
      for(w=0; w<8; w++)
        {
          /*
            except the ones aligned with either
            the link or the v direction
          */
          if( (w&3) == u || (w&3) == (v&3) || (w&3) == orth ) { continue; }
          
          const int w1((w+4)&7); // direction opposite to w
          
          // loop over all directions (+ve and -ve)
          for (x=0;x<8;x++)
            {
              /*
                except the ones aligned with
                either the link, v or w
              */
              if( (x&3) == u     || 
                  (x&3) == (v&3) ||
                  (x&3) == (w&3) ||
                  (x&3) == orth    ) { continue; }
              
              const int x1((x+4)&7); // direction opposite to x
              
              // the chair (or at least a subset of)
              dir[0] = v; 
              dir[1] = w;
              dir[2] = x;
              dir[3] = u;
              dir[4] = x1;
              dir[5] = w1; 
              dir[6] = v1;
                                          
              latt.PathOrdProdPlus(acumulate_mp, pos, dir, 7); 
            }
        }// all directions 
    } // all directions
  // 18 is the matrix size
  moveMem((Float *) &link, (Float*)&acumulate_mp, 18*sizeof(Float));
}


void lepage_staple( Lattice& latt,  Matrix& link ,
                    int *pos, int u , int orth )
{
  Matrix acumulate_mp; acumulate_mp.ZeroMatrix();
  int dir[5];
  //loop over all directions (+ve and -ve)
  int v;
  for ( v=0; v<8; v++ )
    {      
      // except the ones the link is aligned with 
      if((v&3)==u||(v&3)==orth) { continue; }
      
      const int v1((v+4)&7); // direction opposite to v
      
      dir[0] = v;
      dir[1] = v;
      dir[2] = u;
      dir[3] = v1;
      dir[4] = v1;
      
      latt.PathOrdProdPlus(acumulate_mp, pos, dir, 5); 
    }
  // 18 is the matrix size
  moveMem((Float *) &link, (Float*)&acumulate_mp, 18*sizeof(Float));
}


AlgApeSmear::AlgApeSmear(Lattice&     lat,
            CommonArg*   ca ,
            ApeSmearArg* asa,
	    int 	 in_bool_su3_proj):
  AlgSmear(lat,ca,in_bool_su3_proj),
  cname("AlgApeSmear")
{
  c = asa->coef;
  set_tol   (asa->tolerance);
  set_orthog(asa->orthog);
}
  
/*!
  If an output file is specified in the CommonArg argument, then
  the smearing coefficients are written to the file.

*/
void AlgApeSmear::run()
{
  if(common_arg->filename != 0){
    FILE* f = Fopen(common_arg->filename, "a");
    if(!f) ERR.FileA(cname, "run", common_arg->filename);
    Fprintf(f,"AlgApeSmear: coef = %e ",c);
    // YA changed AlgApeSmear being able for no projection, print in that case
    if( ! ifSu3Proj() )
      Fprintf(f,"with NO SU(3) projection");
    Fprintf(f,"\n");
    Fclose(f);
  }
  AlgSmear::run();
}



void AlgApeSmear::smear_link(Matrix& link,
                             int*     pos,
                             int       mu)
{
  Lattice& lattice(AlgLattice());
  Matrix stap;
  three_staple(lattice,stap,pos,mu,get_orthog());
  link*=(1.0-c);
  stap*=c/6.0;
  link+=stap;
}


AlgKineticSmear::AlgKineticSmear(Lattice&   lat,
                                 CommonArg* ca ,
                                 KineticSmearArg* ksa ):
  AlgSmear(lat,ca,0),
  cname("AlgKineticSmear")
{
  set_orthog(ksa->orthog);
  _coef[0] = ksa->single_link;
  _coef[1] = ksa->three_link;
  _coef[2] = ksa->five_link;
  _coef[3] = ksa->seven_link;
  _coef[4] = ksa->lepage;
}


/*!
  If an output file is specified in the CommonArg argument, then
  the smearing coefficients are written to the file.

*/
void AlgKineticSmear::run()
{
  if(common_arg->filename != 0){
    FILE* f = Fopen(common_arg->filename, "a");
    if(!f) ERR.FileA(cname, "run", common_arg->filename);
    for (int i=0;i<5;++i) Fprintf(f,"coef %2i : %e \n",i,(Float)_coef[i]);
    Fclose(f);
  }
  AlgSmear::run();
}


void AlgKineticSmear::smear_link( Matrix& link,
                                  int*    pos,
                                  int      mu )
{

  typedef void (*staple_func)(Lattice&,Matrix&,int*,int,int);
    
  Lattice& lattice(AlgLattice());
  staple_func funcs[]={ &three_staple,
                        &five_staple ,
                        &seven_staple,
                        &lepage_staple };
  link*=_coef[0];
  Matrix stap; 

  for (int i=1;i<5;++i)
    {
      if ( _coef[i] != 0 ) continue;
      
      (*(funcs[i-1]))(lattice,stap,pos,mu,get_orthog());
      stap*=_coef[i];
      link+=stap;
      
    }
}


AlgHypSmear::AlgHypSmear( Lattice&     lat, 
                          CommonArg*   ca ,
                          HypSmearArg* hsa ): 
  AlgSmear(lat,ca,1),
  cname("AlgHypSmear")
{
  set_tol   (hsa->tolerance);
  set_orthog(hsa->orthog);
  c1 = hsa->c1;
  c2 = hsa->c2;
  c3 = hsa->c3;
}
  
const Matrix AlgHypSmear::GetLink( Lattice& lat, const int* x, int mu )
{
  int link_site[4];
  int i;
  for (i=0;i<4;i++){ link_site[i] = x[i]; }
  const int abs_dir ( mu & 3  );
  const int dir_sign( mu >> 2  );
  link_site[abs_dir] -=dir_sign;

  if( dir_sign){
    Matrix tmp;
    tmp.Dagger( *(lat.GetBufferedLink( link_site, abs_dir )) );
    return tmp;
  } else
    return *(lat.GetBufferedLink( link_site, abs_dir ));
}


void AlgHypSmear::get_vbar( Matrix& link, int *pos, int mu, int nu, int rho )
{
  Lattice& lat(AlgLattice());
  Matrix accum; accum.ZeroMatrix();
  int v;
  int dir[3];
  for ( v=0;v<8;v++)
    {
      const int mv(v&3);
      if (mv==(mu&3)||mv==(nu&3)||mv==(rho&3)) { continue; }
      const int v1((v+4)&7); // direction opposite to v
      dir[0] = v;
      dir[1] = mu;
      dir[2] = v1;
      lat.PathOrdProdPlus(accum, pos, dir, 3);
    }
  link =  GetLink(lat,pos,mu);
  link *= (1-c3);
  accum*= c3/2;
  link+=accum;
  // project down to su3
  su3_proj( link , get_tol());
}

void AlgHypSmear::get_vtilde( Matrix& link , int *pos_in, int mu_in, int nu )
{
  Lattice& latt(AlgLattice());
  link.ZeroMatrix();
  Matrix tmp1,tmp2,tmp3,tmp4;
  Matrix stap;
  
  const int sign_mu_in( mu_in >> 2  );
  const int mu ( mu_in & 3  );

  int pos[4] = { pos_in[0],pos_in[1],pos_in[2],pos_in[3] };
  if(sign_mu_in) pos[mu]--;

  // +1 in mu direction
  int pos_mu[4] = { pos[0],pos[1],pos[2],pos[3] };
  pos_mu[mu]++;

  int dir[2]; dir[0]=-1;
  int i;
  // dir1 and dir2 should be the two directions orthogonal to
  // mu and nu 
  for ( i=0;i<4;i++) 
    { 
      if ( i==(mu&3) || i==(nu&3) ) { continue; }
      if ( dir[0] < 0 ) { dir[0] = i ; }
      else { dir[1] = i; }
    }

  stap.ZeroMatrix();
  for(i=0;i<2;i++)
    {
      // forward 
      get_vbar(tmp1,pos   ,dir[i],nu    ,mu); pos[dir[i]]++;
      get_vbar(tmp2,pos   ,mu    ,dir[i],nu); pos[dir[i]]--;
      get_vbar(tmp3,pos_mu,dir[i],nu    ,mu);
     
      tmp4.DotMEqual( tmp1, tmp2 );
      tmp1.Dagger   ( tmp3 );
      stap.DotMPlus ( tmp4, tmp1 );
    
      // backwards 
      get_vbar(tmp1,pos   ,(dir[i]+4)&7,nu,mu); pos[dir[i]]--;
      get_vbar(tmp2,pos   ,mu,(dir[i]+4)&7,nu); pos[dir[i]]++;
      get_vbar(tmp3,pos_mu,(dir[i]+4)&7,nu,mu);
      
      tmp4.DotMEqual( tmp1, tmp2 );
      tmp1.Dagger   ( tmp3 );
      stap.DotMPlus ( tmp4, tmp1 );
    }
  
  link = GetLink(latt,pos,mu);
  link*=(1-c2);
  stap*=c2/4;
  link+=stap;
  // project down to su3
  su3_proj( link, get_tol() );

  if( sign_mu_in ){
    tmp1.Dagger( link );
    link = tmp1;
  }
}



void AlgHypSmear::smear_link(Matrix& link,
                             int*     pos,
                             int       mu)
{
  Lattice& latt(AlgLattice());
  link.ZeroMatrix();
  Matrix tmp1,tmp2,tmp3,tmp4;
  Matrix stap;
  
  // +1 in mu direction
  int pos_mu[4];
  pos_mu[0] = pos[0];
  pos_mu[1] = pos[1];
  pos_mu[2] = pos[2];
  pos_mu[3] = pos[3];
  pos_mu[mu]++;

  int dir[3]; dir[0]=-1; dir[1]=-1;
  int i;
  // dir1-3 should be orthogonal to mu
  for ( i=0;i<4;i++) 
    { 
      if ( i==(mu&3) ) { continue; }
      if      ( dir[0] < 0 ) { dir[0] = i; }
      else if ( dir[1] < 0 ) { dir[1] = i; } 
      else                   { dir[2] = i; }
    }

  stap.ZeroMatrix();
  for ( i=0;i<3;i++)
    {
      // forward 
      get_vtilde(tmp1,pos,dir[i],mu);

      pos[dir[i]]++;
      get_vtilde(tmp2,pos,mu,dir[i]);

      pos[dir[i]]--;
      get_vtilde(tmp3,pos_mu,dir[i],mu);
      
      tmp4.DotMEqual( tmp1, tmp2 );
      tmp1.Dagger   ( tmp3 );
      stap.DotMPlus ( tmp4, tmp1 );
      
      // backwards 
      get_vtilde(tmp1,pos,(dir[i]+4)&7,mu);

      pos[dir[i]]--;
      get_vtilde(tmp2,pos,mu,(dir[i]+4)&7);

      pos[dir[i]]++;
      get_vtilde(tmp3,pos_mu,(dir[i]+4)&7,mu);
      
      tmp4.DotMEqual( tmp1, tmp2 );
      tmp1.Dagger   ( tmp3 );
      stap.DotMPlus ( tmp4, tmp1 );
    }
  link = GetLink(latt,pos,mu);
  link*=(1-c1);
  stap*=c1/6;
  link+=stap;
  
  // don't su3 project here because we'll do it in run()
}

/*!
  If an output file is specified in the CommonArg argument, then
  the smearing coefficients are written to the file.

  \pre The smearing coefficents should be set with set_c1, set_2 and set_c3.
*/
void AlgHypSmear::run()
{
  //  if ( get_orthog() >=0 || get_orthog() <4 )
  //   ^ YA: this will bring anything to ERR
  if ( get_orthog() >=0 && get_orthog() <4 )
      ERR.General(cname, "run",
		  "Bad value %d for orthogonal direction", get_orthog());
   

  if(common_arg->filename != 0){
      FILE* f = Fopen(common_arg->filename, "a");
      if(!f) ERR.FileA(cname, "run", common_arg->filename);
      Fprintf(f,"AlgHypSmear hit: c1=%e  c2=%e  c3=%e \n",c1,c2,c3);
      Fclose(f);
  }

  AlgSmear::run();
}  
CPS_END_NAMESPACE
