#include<config.h>

/*!\file
  \brief Helper routines for the fermionic force term

  $Id: force_product_sum.C,v 1.2 2004-05-10 15:26:55 zs Exp $
*/
//--------------------------------------------------------------------

#include <util/lattice.h>
#include <util/gjp.h>

USING_NAMESPACE_CPS




// f += coeff * U_dir(x) * v(x)^dagger
// Note that STAG order stores hermitian conjugate of links.

void Fasqtad::force_product_sum(const Matrix *v,  int dir,
				IFloat coeff, Matrix *f){

    Matrix m, md;
    int s[4];
    
    for(s[3]=0; s[3]<node_sites[3]; s[3]++)
	for(s[2]=0; s[2]<node_sites[2]; s[2]++)
	    for(s[1]=0; s[1]<node_sites[1]; s[1]++)
		for(s[0]=0; s[0]<node_sites[0]; s[0]++){
		    int  n = FsiteOffset(s); 
		    m.DotMEqual(v[n], *(GaugeField()+GsiteOffset(s)+dir));
		    md.Dagger(m);
#ifdef MILC_COMPATIBILITY
		    if(parity(s)==CHKB_ODD) md *= -coeff;
#else
		    if(parity(s)==CHKB_EVEN) md *= -coeff; 
#endif
		    else md *= coeff;
		    f[n] += md;
		}

}


// f += u(vw)^dagger

void Fasqtad::force_product_sum(const Matrix *u, const Matrix *v,
				const Matrix *w, IFloat coeff, Matrix *f){

    Matrix m, md;
    int s[4];
    
    for(s[3]=0; s[3]<node_sites[3]; s[3]++)
	for(s[2]=0; s[2]<node_sites[2]; s[2]++)
	    for(s[1]=0; s[1]<node_sites[1]; s[1]++)
		for(s[0]=0; s[0]<node_sites[0]; s[0]++){
		    int n = FsiteOffset(s);
		    m.DotMEqual(v[n], w[n]);
		    md.Dagger(m);
#ifdef MILC_COMPATIBILITY
		    if(parity(s)==CHKB_ODD) md *= -coeff;
#else
		    if(parity(s)==CHKB_EVEN) md *= -coeff;
#endif
		    else md *= coeff;
		    f[n].DotMPlus(u[n], md);
		}    

}


// f += coeff v w^dagger

void Fasqtad::force_product_sum(const Matrix *v, const Matrix *w,
				    IFloat coeff, Matrix *f){

    Matrix m, md;
    int s[4];
    
    for(s[3]=0; s[3]<node_sites[3]; s[3]++)
	for(s[2]=0; s[2]<node_sites[2]; s[2]++)
	    for(s[1]=0; s[1]<node_sites[1]; s[1]++)
		for(s[0]=0; s[0]<node_sites[0]; s[0]++){
		    int n = FsiteOffset(s);
		    md.Dagger(w[n]);
		    m.DotMEqual(v[n], md);
#ifdef MILC_COMPATIBILITY
		    if(parity(s)==CHKB_ODD) m *= -coeff;
#else
		    if(parity(s)==CHKB_EVEN) m *= -coeff;
#endif
		    else m *= coeff;
		    f[n] += m;
		}

}

// f += coeff (v w)^dagger

void Fasqtad::force_product_d_sum(const Matrix *v, const Matrix *w,
				    IFloat coeff, Matrix *f){

    Matrix m, md;
    int n, s[4];
    
    for(s[3]=0; s[3]<node_sites[3]; s[3]++)
	for(s[2]=0; s[2]<node_sites[2]; s[2]++)
	    for(s[1]=0; s[1]<node_sites[1]; s[1]++)
		for(s[0]=0; s[0]<node_sites[0]; s[0]++){
		    n = FsiteOffset(s);
		    m.DotMEqual(v[n], w[n]);
		    md.Dagger(m);
#ifdef MILC_COMPATIBILITY
		    if(parity(s)==CHKB_ODD) md *= -coeff;
#else
		    if(parity(s)==CHKB_EVEN) md *= -coeff;
#endif
		    else md *= coeff;
		    f[n] += md;
		}

}



// The outer product of v and w, multiplied by the parity factor and
// coeffient, and added to the force f.
// N.B. The force term is multiplied by -1 on ODD parity sites. This the
// the MILC convention and is here to test against the MILC code.
// This should eventually be changed to EVEN for the CPS.


void Fasqtad::force_product_sum(const Vector *v, const Vector *w,
				    IFloat coeff, Matrix *f){

    Matrix m;
    int n, s[4];
    
    for(s[3]=0; s[3]<GJP.TnodeSites(); s[3]++)
	for(s[2]=0; s[2]<GJP.ZnodeSites(); s[2]++)
	    for(s[1]=0; s[1]<GJP.YnodeSites(); s[1]++)
		for(s[0]=0; s[0]<GJP.XnodeSites(); s[0]++){
		    n = FsiteOffset(s);
		    m.Cross2(v[n], w[n]);
// Note the extra factor of 2 from Matrix::Cross2; this is corrected by
// the factor of 1/2 in Matrix::TrLessAntiHermMatrix used in
// Fasqtad::update_momenta.
#ifdef MILC_COMPATIBILITY
		    if(parity(s)==CHKB_ODD) m *= -coeff;    
#else
		    if(parity(s)==CHKB_EVEN) m *= -coeff;    
#endif
		    else m *= coeff;
		    f[n] += m;
		}

}
