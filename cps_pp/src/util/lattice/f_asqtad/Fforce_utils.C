#include<config.h>
//----------------------------------------------------------------------
/*!\file
  \brief  Routines used in the Asqtad RHMC fermion force calculation

  $Id: Fforce_utils.C,v 1.4 2004-07-01 17:43:46 chulwoo Exp $
*/
//----------------------------------------------------------------------


#include <util/lattice.h>
#include <util/gjp.h>
#include <comms/scu.h>
#include<comms/glb.h>

USING_NAMESPACE_CPS


// Computes sum[x] = vect[x] vect[x + hop dir]^dagger
// where the sum is over n_vect vectors and the hop is in a forward direction.
// Quick (to write, at least) implementation; I do n_dir directions dir
// at once in this routine in case anyone wants to optimise this.

void Fasqtad::vvpd(Vector **vect, int n_vect,
		   const int *dir, int n_dir, int hop, Matrix **sum){
  
  Vector vp;
  int s[4];
  Matrix m;

  for(s[3] = 0; s[3]<node_sites[3]; s[3]++)
    for(s[2] = 0; s[2]<node_sites[2]; s[2]++)
      for(s[1] = 0; s[1]<node_sites[1]; s[1]++)
	for(s[0] = 0; s[0]<node_sites[0]; s[0]++) {
	  int  n = FsiteOffset(s);
	  
	  for(int d=0; d<n_dir; d++){ 
	    sum[d][n].ZeroMatrix();  

	    if(s[dir[d]]+hop>node_sites[dir[d]]-1){
	      
	      int save_coord = s[dir[d]];
	      s[dir[d]] = (s[dir[d]]+hop)%node_sites[dir[d]];
	      int np = FsiteOffset(s);
	      s[dir[d]] = save_coord;
	      
	      for(int v=0; v<n_vect; v++){
		getPlusData((IFloat*)&vp, (IFloat*)&vect[v][np],
			    FsiteSize(), dir[d]);
		m.Cross2(vect[v][n], vp);
		// Note the factor of 2 from Matrix::Cross2; this is corrected by the factor
		// of 1/2 in Matrix::TrLessAntiHermMatrix used in Fasqtad::update_momenta.
		sum[d][n] += m;
	      }
	      
	    }else{
	      
	      s[dir[d]] += hop;
	      int np = FsiteOffset(s);
	      s[dir[d]] -= hop;
	      
	      for(int v=0; v<n_vect; v++){
		m.Cross2(vect[v][n], vect[v][np]);
		// see comment above.
		sum[d][n] += m;
	      }
	    }
	  }
	}
  
}



// u[x] = v[x+dir] for n_dir forward or backward directions dir.

void Fasqtad::shift_field(Matrix **v, const int *dir, int n_dir,
			  int hop, Matrix **u){

    Matrix m;
    int s[4], spm[4];
    
    for(s[3]=0; s[3]<node_sites[3]; s[3]++)
	for(s[2]=0; s[2]<node_sites[2]; s[2]++)
	    for(s[1]=0; s[1]<node_sites[1]; s[1]++)
		for(s[0]=0; s[0]<node_sites[0]; s[0]++){
		    int  n = FsiteOffset(s);
		    
		    for(int d=0; d<n_dir; d++){
			int mu = dir[d];
			bool backwards = mu%2 ? true : false;
			mu /= 2;

			if(backwards){
			    if(s[mu]-hop<0){

				int save_coord = s[mu];
				s[mu] -= hop;
				while(s[mu]<0) s[mu] += node_sites[mu];
				int np = FsiteOffset(s);
				s[mu] = save_coord;

				getMinusData((IFloat*)&u[d][n],
					     (IFloat*)&v[d][np],
					     MATRIX_SIZE, mu);
			    
			    }else{

				s[mu] -= hop;
				int np = FsiteOffset(s);
				s[mu] += hop;
				
				u[d][n] = v[d][np];

			    }
			}else{
			    if(s[mu]+hop>node_sites[mu]-1){

				int save_coord = s[mu];
				s[mu] = (s[mu]+hop)%node_sites[mu];
				int np = FsiteOffset(s);
				s[mu] = save_coord;

				getPlusData((IFloat*)&u[d][n],
					    (IFloat*)&v[d][np],
					    MATRIX_SIZE, mu);
			    
			    }else{

				s[mu] += hop;
				int np = FsiteOffset(s);
				s[mu] -= hop;

				u[d][n] = v[d][np];

			    }
			}
		    }
		}

}


// u[-/+nu](x) = U_[-/+nu](x) 

void Fasqtad::shift_link(Matrix **u, const int *dir, int n_dir){

    int s[4];
    Matrix m; 
    
    for(s[3]=0; s[3]<node_sites[3]; s[3]++)
	for(s[2]=0; s[2]<node_sites[2]; s[2]++)
	    for(s[1]=0; s[1]<node_sites[1]; s[1]++)
		for(s[0]=0; s[0]<node_sites[0]; s[0]++){

		    for(int d=0; d<n_dir; d++){
			int nu = dir[d];
			bool backwards = nu%2 ? true : false;
			nu /= 2;

			if(backwards){
			    if(s[nu]==0){
				s[nu] = node_sites[nu]-1;
				getMinusData((IFloat*)&m,
					     (IFloat*)(GaugeField()+
						       GsiteOffset(s)+ nu),
					     MATRIX_SIZE, nu);
				s[nu] = 0;
			    }else{
				s[nu]--;
				m = *(GaugeField()+GsiteOffset(s)+nu);
				s[nu]++;				
			    }
			}else{
			    m.Dagger(*(GaugeField()+GsiteOffset(s)+nu));
			}
			u[d][FsiteOffset(s)] = m;
		    }
		}

}


// The update of the momentum with the force

void Fasqtad::update_momenta(Matrix **force, IFloat dt, Matrix *mom, int
mul_parity){

    Matrix mf, mfd;
    double dt_tmp;

    int s[4];
    for(s[3]=0; s[3]<node_sites[3]; s[3]++)
	for(s[2]=0; s[2]<node_sites[2]; s[2]++)
	    for(s[1]=0; s[1]<node_sites[1]; s[1]++)
		for(s[0]=0; s[0]<node_sites[0]; s[0]++){

		    Matrix *ip = mom+GsiteOffset(s);

		    for (int mu=0; mu<4; mu++){			
			mf = force[mu][FsiteOffset(s)];
			mfd.Dagger((IFloat*)&mf);
			mf.TrLessAntiHermMatrix(mfd);
#ifdef MILC_COMPATIBILITY
			mf *= 0.5;	
#endif
			if(mul_parity && parity(s) == CHKB_ODD) dt_tmp =-dt;
			else  dt_tmp = dt;
			fTimesV1PlusV2((IFloat*)(ip+mu), dt_tmp, (IFloat*)&mf,
				       (IFloat*)(ip+mu),
				       MATRIX_SIZE);
		    }
		}
	ForceFlops += GJP.VolNodeSites()*54;	
    
}


// The parity of the lattice site n

ChkbType Fasqtad::parity(const int *n){

    int d, p;
    for(p=0, d=0; d<4; d++) p += n[d];
    return  p%2?CHKB_ODD:CHKB_EVEN;
    
}

