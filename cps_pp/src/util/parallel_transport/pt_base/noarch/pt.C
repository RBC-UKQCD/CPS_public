/*! \file
  \brief  Functions used by the ParTransAsqtad class.
  
  $Id: pt.C,v 1.6 2004-09-02 17:00:05 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-09-02 17:00:05 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/parallel_transport/pt_base/noarch/pt.C,v 1.6 2004-09-02 17:00:05 zs Exp $
//  $Id: pt.C,v 1.6 2004-09-02 17:00:05 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: pt.C,v $
//  $Revision: 1.6 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/parallel_transport/pt_base/noarch/pt.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
#include <config.h>
#include <util/gjp.h>
#include <util/pt.h>

CPS_START_NAMESPACE

enum {MATRIX_SIZE = 18,NUM_DIR=8};
static int vol;
static Lattice  *Lat;
static int node_sites[5];
void pt_mat(int N, IFloat **fout, IFloat **fin, const int *dir){

    IFloat *u, *v;
    int s[4], spm[4];
    Matrix ud;
    Matrix *vin[NUM_DIR],*vout[NUM_DIR];


    for(int d=0; d<N; d++){
	vin[d] = (Matrix *)fin[d];
	vout[d] = (Matrix *)fout[d];

	int mu = dir[d];
	bool backwards = mu%2 ? true : false;
        mu /= 2;

	for(spm[3]=s[3]=0; s[3]<GJP.TnodeSites(); spm[3]++, s[3]++)
	    for(spm[2]=s[2]=0; s[2]<GJP.ZnodeSites(); spm[2]++, s[2]++)
		for(spm[1]=s[1]=0; s[1]<GJP.YnodeSites(); spm[1]++, s[1]++)
		    for(spm[0]=s[0]=0; s[0]<GJP.XnodeSites(); spm[0]++, s[0]++) {
			if(backwards){
			    spm[mu] = (s[mu]-1+GJP.NodeSites(mu))%GJP.NodeSites(mu);
			    u = (IFloat*)(Lat->GaugeField()+Lat->GsiteOffset(spm)+mu);
			    v = (IFloat*)&vin[d][Lat->FsiteOffset(spm)];
			    mDotMEqual((IFloat*)&vout[d][Lat->FsiteOffset(s)],
				       u, v);
			}else{
			    spm[mu] = (s[mu]+1)%GJP.NodeSites(mu);
			    ud.Dagger(*(Lat->GaugeField()+Lat->GsiteOffset(s)+mu));
			    u = (IFloat*)&ud;
			    v = (IFloat*)&vin[d][Lat->FsiteOffset(spm)];
			    // STAG stores hermitian conjugate links
			    mDotMEqual((IFloat*)&vout[d][Lat->FsiteOffset(s)],
				      u, v);			    
			}
			spm[mu] = s[mu];

		    }
			
    }
    ParTrans::PTflops += 198*N*vol;
}

void pt_1vec(int N, IFloat **fout, IFloat **fin, const int *dir){

    IFloat *u, *v;
    int s[4], spm[4];
    Vector *vin[NUM_DIR],*vout[NUM_DIR];

    for(int d=0; d<N; d++){
        vin[d] = (Vector *)fin[d];
        vout[d] = (Vector *)fout[d];

	int mu = dir[d];
	bool backwards = mu%2 ? true : false;
        mu /= 2;

	for(spm[3]=s[3]=0; s[3]<GJP.TnodeSites(); spm[3]++, s[3]++)
	    for(spm[2]=s[2]=0; s[2]<GJP.ZnodeSites(); spm[2]++, s[2]++)
		for(spm[1]=s[1]=0; s[1]<GJP.YnodeSites(); spm[1]++, s[1]++)
		    for(spm[0]=s[0]=0; s[0]<GJP.XnodeSites(); spm[0]++, s[0]++) {
			if(backwards){
			    spm[mu] = (s[mu]-1+GJP.NodeSites(mu))%GJP.NodeSites(mu);
			    u = (IFloat*)(Lat->GaugeField()+Lat->GsiteOffset(spm)+mu);
			    v = (IFloat*)&vin[d][Lat->FsiteOffset(spm)];
			    uDotXEqual((IFloat*)&vout[d][Lat->FsiteOffset(s)],
				      u, v);			    
			}else{
			    spm[mu] = (s[mu]+1)%GJP.NodeSites(mu);
			    u = (IFloat*)(Lat->GaugeField()+Lat->GsiteOffset(s)+mu);
			    v = (IFloat*)&vin[d][Lat->FsiteOffset(spm)];
			    // STAG stores hermitian conjugate links
			    uDagDotXEqual((IFloat*)&vout[d][Lat->FsiteOffset(s)],
				      u, v);			    
			}
			spm[mu] = s[mu];

		    }
			
    }
    ParTrans::PTflops += 198*N*vol;
}

/*! 
  Computes sum[x] = vect[x] vect[x + hop dir]^dagger
  where the sum is over n_vect vectors and the hop is in a forward direction.
*/
void pt_vvpd(Vector **vect, int n_vect,
	     const int *dir, int n_dir, int hop, Matrix **sum){
  
  int s[4];
  Matrix m;
  
  Vector vp;

  for(s[3] = 0; s[3]<node_sites[3]; s[3]++)
    for(s[2] = 0; s[2]<node_sites[2]; s[2]++)
      for(s[1] = 0; s[1]<node_sites[1]; s[1]++)
	for(s[0] = 0; s[0]<node_sites[0]; s[0]++) {
	  int n = Lat->FsiteOffset(s);
	  
	  for(int d=0; d<n_dir; d++){ 
	    sum[d][n].ZeroMatrix();  
	    if(s[dir[d]]+hop>node_sites[dir[d]]-1){
	      
	      int save_coord = s[dir[d]];
	      s[dir[d]] = (s[dir[d]]+hop)%node_sites[dir[d]];
	      int np = Lat->FsiteOffset(s);
	      s[dir[d]] = save_coord;
	      
	      for(int v=0; v<n_vect; v++){
		getPlusData((IFloat*)&vp, (IFloat*)&vect[v][np],
			    Lat->FsiteSize(), dir[d]);
		m.Cross2(vect[v][n], vp);
		// Note the factor of 2 from Matrix::Cross2; this is corrected by the factor
		// of 1/2 in Matrix::TrLessAntiHermMatrix used in Fasqtad::update_momenta.
		sum[d][n] += m;
	      }
	    }else{
	      s[dir[d]] += hop;
	      int np = Lat->FsiteOffset(s);
	      s[dir[d]] -= hop;
	      
	      for(int v=0; v<n_vect; v++){
		m.Cross2(vect[v][n], vect[v][np]);
		// see comment above.
		sum[d][n] += m;
	      }
	    }
	  }
	}

  ParTrans::PTflops += 90*n_vect*n_dir*vol;
}



//! u[x] = v[x+dir] for n_dir forward or backward directions dir.
void pt_shift_field(Matrix **v, const int *dir, int n_dir,
		    int hop, Matrix **u){
    

    int s[4];
  
    for(int d=0; d<n_dir; d++){
	int mu = dir[d];
	bool backwards = mu%2 ? true : false;
	mu /= 2;
    
	if(backwards){
	    for(s[3]=0; s[3]<node_sites[3]; s[3]++)
		for(s[2]=0; s[2]<node_sites[2]; s[2]++)
		    for(s[1]=0; s[1]<node_sites[1]; s[1]++)
			for(s[0]=0; s[0]<node_sites[0]; s[0]++){
			    int n = Lat->FsiteOffset(s);
	      
			    if(s[mu]-hop<0){
		
				int save_coord = s[mu];
				s[mu] -= hop;
				while(s[mu]<0) s[mu] += node_sites[mu];
				int np = Lat->FsiteOffset(s);
				s[mu] = save_coord;
				getMinusData((IFloat*)&u[d][n], (IFloat*)&v[d][np], MATRIX_SIZE, mu);
		
			    }else{
		
				s[mu] -= hop;
				int np = Lat->FsiteOffset(s);
				s[mu] += hop;
				u[d][n] = v[d][np];
		
			    }
			}
	} else {
	    for(s[3]=0; s[3]<node_sites[3]; s[3]++)
		for(s[2]=0; s[2]<node_sites[2]; s[2]++)
		    for(s[1]=0; s[1]<node_sites[1]; s[1]++)
			for(s[0]=0; s[0]<node_sites[0]; s[0]++){
			    int n = Lat->FsiteOffset(s);
	    
			    if(s[mu]+hop>node_sites[mu]-1){
		
				int save_coord = s[mu];
				s[mu] = (s[mu]+hop)%node_sites[mu];
				int np = Lat->FsiteOffset(s);
				s[mu] = save_coord;
				getPlusData((IFloat*)&u[d][n], (IFloat*)&v[d][np], MATRIX_SIZE, mu);
		
			    }else{
		
				s[mu] += hop;
				int np = Lat->FsiteOffset(s);
				s[mu] -= hop;
				u[d][n] = v[d][np];
		
			    }
			}
	}
    }
  
}

//! u[-/+nu](x) = U_[-/+nu](x) 
void pt_shift_link(Matrix **u, const int *dir, int n_dir){

    int s[4];
    Matrix m; 
  
    for(int d=0; d<n_dir; d++){
	int nu = dir[d];
	bool backwards = nu%2 ? true : false;
	nu /= 2;
    
	if(backwards){
	    for(s[3]=0; s[3]<node_sites[3]; s[3]++)
		for(s[2]=0; s[2]<node_sites[2]; s[2]++)
		    for(s[1]=0; s[1]<node_sites[1]; s[1]++)
			for(s[0]=0; s[0]<node_sites[0]; s[0]++){	      
			    if(s[nu]==0){
				s[nu] = node_sites[nu]-1;
				getMinusData((IFloat*)&m, 
					     (IFloat*)(Lat->GaugeField()+ Lat->GsiteOffset(s)+ nu), MATRIX_SIZE, nu);
				s[nu] = 0;
			    }else{
				s[nu]--;
				m = *(Lat->GaugeField()+Lat->GsiteOffset(s)+nu);
				s[nu]++;				
			    }
			    u[d][Lat->FsiteOffset(s)] = m;
			}
	} else {
      
	    for(s[3]=0; s[3]<node_sites[3]; s[3]++)
		for(s[2]=0; s[2]<node_sites[2]; s[2]++)
		    for(s[1]=0; s[1]<node_sites[1]; s[1]++)
			for(s[0]=0; s[0]<node_sites[0]; s[0]++){
			    m.Dagger(*(Lat->GaugeField()+Lat->GsiteOffset(s)+nu));
			    u[d][Lat->FsiteOffset(s)] = m;
			}
      
	}
    
    }
  
}



void pt_init(Lattice &latt) {
  Lat = &latt;
  vol = GJP.VolNodeSites();

  node_sites[0] = GJP.XnodeSites();
  node_sites[1] = GJP.YnodeSites();
  node_sites[2] = GJP.ZnodeSites();
  node_sites[3] = GJP.TnodeSites();
  node_sites[4] = GJP.SnodeSites();
}

void pt_delete(){
}

void pt_delete_g(){
}

void pt_init_g(){
}

CPS_END_NAMESPACE
	
