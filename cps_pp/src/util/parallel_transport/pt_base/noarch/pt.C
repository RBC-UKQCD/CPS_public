#include <config.h>
/*! \file
  \brief  Functions used by the ParTransAsqtad class.
  
  $Id: pt.C,v 1.12 2012-08-02 21:20:01 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2012-08-02 21:20:01 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/parallel_transport/pt_base/noarch/pt.C,v 1.12 2012-08-02 21:20:01 chulwoo Exp $
//  $Id: pt.C,v 1.12 2012-08-02 21:20:01 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: pt.C,v $
//  $Revision: 1.12 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/parallel_transport/pt_base/noarch/pt.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
#include <util/gjp.h>
#include <util/pt.h>
#include <comms/scu.h>

CPS_START_NAMESPACE

enum {MATRIX_SIZE = 18,NUM_DIR=8, VECT_LEN = 6};
static int vol;
static Lattice  *Lat;
static int node_sites[5];
//--------------------------------------------------------------------

#ifndef USE_QMP
void pt_mat(int N, IFloat **fout, IFloat **fin, const int *dir){

    IFloat *u, *v;
    int s[4], spm[4];
    Matrix ud;
    Matrix *vin[NUM_DIR],*vout[NUM_DIR];
    Matrix m, mp;


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
			  if(s[mu] == 0) {
			    spm[mu] = GJP.NodeSites(mu)-1;
			    //Get gauge field and matrix field
			      getMinusData((IFloat*)&m, (IFloat*)(Lat->GaugeField()+Lat->GsiteOffset(spm)+mu), MATRIX_SIZE, mu);
			      getMinusData((IFloat*)&mp, (IFloat*)(&vin[d][Lat->FsiteOffset(spm)]), MATRIX_SIZE, mu);
			      // STAG stores hermitian conjugate links
			      mDotMEqual((IFloat*)&vout[d][Lat->FsiteOffset(s)],
					    (IFloat *)&m, (IFloat *)&mp);		    			    
			  }
			  else {
			    spm[mu] = s[mu]-1;
			    u = (IFloat*)(Lat->GaugeField()+Lat->GsiteOffset(spm)+mu);
			    v = (IFloat*)&vin[d][Lat->FsiteOffset(spm)];
			    mDotMEqual((IFloat*)&vout[d][Lat->FsiteOffset(s)],
				       u, v);
			  }
			}else{
			  if(s[mu] == GJP.NodeSites(mu)-1) {
			    spm[mu] = 0;
			    //Get only matrix field
			    getPlusData((IFloat*)&mp, (IFloat*)(&vin[d][Lat->FsiteOffset(spm)]), MATRIX_SIZE, mu);
			    ud.Dagger(*(Lat->GaugeField()+Lat->GsiteOffset(s)+mu));
			    u = (IFloat*)&ud;
			    // STAG stores hermitian conjugate links
			    mDotMEqual((IFloat*)&vout[d][Lat->FsiteOffset(s)],
				      u, (IFloat *)&mp);	
			  }
			  else {
			    spm[mu] = spm[mu]+1;
			    ud.Dagger(*(Lat->GaugeField()+Lat->GsiteOffset(s)+mu));
			    u = (IFloat*)&ud;
			    v = (IFloat*)&vin[d][Lat->FsiteOffset(spm)];
			    mDotMEqual((IFloat*)&vout[d][Lat->FsiteOffset(s)],
				       u, v);	
			  }
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

    Matrix m;
    Vector vp;

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
			  if(s[mu] == 0) {
			    spm[mu] = GJP.NodeSites(mu)-1;
			    //Get gauge field and fermion field
			      getMinusData((IFloat*)&m, (IFloat*)(Lat->GaugeField()+Lat->GsiteOffset(spm)+mu), MATRIX_SIZE, mu);
			      getMinusData((IFloat*)&vp, (IFloat*)(&vin[d][Lat->FsiteOffset(spm)]), VECT_LEN, mu);
			      // STAG stores hermitian conjugate links
			      uDotXEqual((IFloat*)&vout[d][Lat->FsiteOffset(s)],
					    (IFloat *)&m, (IFloat *)&vp);		    			    
			  }
			  else {
			    spm[mu] = s[mu]-1;
			    u = (IFloat*)(Lat->GaugeField()+Lat->GsiteOffset(spm)+mu);
			    v = (IFloat*)&vin[d][Lat->FsiteOffset(spm)];
			    uDotXEqual((IFloat*)&vout[d][Lat->FsiteOffset(s)],
				       u, v);
			  }
			}else{
			  if(s[mu] == GJP.NodeSites(mu)-1) {
			    spm[mu] = 0;
			    //Get only fermion field
			    getPlusData((IFloat*)&vp, (IFloat*)(&vin[d][Lat->FsiteOffset(spm)]), VECT_LEN, mu);
			    u = (IFloat*)(Lat->GaugeField()+Lat->GsiteOffset(s)+mu);
			    // STAG stores hermitian conjugate links
			    uDagDotXEqual((IFloat*)&vout[d][Lat->FsiteOffset(s)],
				      u, (IFloat *)&vp);	
			  }
			  else {
			    spm[mu] = spm[mu]+1;
			    u = (IFloat*)(Lat->GaugeField()+Lat->GsiteOffset(s)+mu);
			    v = (IFloat*)&vin[d][Lat->FsiteOffset(spm)];
			    uDagDotXEqual((IFloat*)&vout[d][Lat->FsiteOffset(s)],
				       u, v);	
			  }
			}
			spm[mu] = s[mu];

		    }
			
#if 0
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
#endif
			
    }
    ParTrans::PTflops += 66*N*vol;
}
#endif

/*! 
  Computes sum[x] = vect2[x] vect[x +/- hop dir]^dagger
  where the sum is over n_vect vectors and the hop is in the forward 
  or negative direction.
*/

void pt_vvpd(IFloat **vect2_v, IFloat ***vect_v, int n_vect, const int *dir, int n_dir, int hop, IFloat **sum_m, int overwrite){
  //Do nothing
  char *cname = "PT";
  char *fname = "pt_vvpd(F***,I,I*,I,I,F**,I)";

  Vector **vect2 = (Vector **) vect2_v;
  Vector ***vect = (Vector ***) vect_v;
  Matrix **sum = (Matrix **) sum_m;
  int s[4];
  Matrix m;
  
  Vector vp;

  for(s[3] = 0; s[3]<node_sites[3]; s[3]++)
    for(s[2] = 0; s[2]<node_sites[2]; s[2]++)
      for(s[1] = 0; s[1]<node_sites[1]; s[1]++)
	for(s[0] = 0; s[0]<node_sites[0]; s[0]++) {
	  int n = Lat->FsiteOffset(s);
	  
	  for(int d=0; d<n_dir; d++){ 
	    if (overwrite == 1)
	      sum[d][n].ZeroMatrix();

	    int mu = dir[d];
	    bool backwards = mu%2 ? true : false;
	    mu /= 2;

	    //Transport fields "forward"
	    if (backwards) {
	      if(s[mu]-hop<0){
		
		int save_coord = s[mu];
		s[mu] = (s[mu]-hop+2*node_sites[mu])%node_sites[mu];
		int np = Lat->FsiteOffset(s);
		s[mu] = save_coord;
		
		for(int v=0; v<n_vect; v++){
		  getMinusData((IFloat*)&vp, (IFloat*)&vect[v][d][np],
			      Lat->FsiteSize(), mu);
		  m.Cross2(vect2[v][n], vp);
		  // Note the factor of 2 from Matrix::Cross2; this is corrected by the factor
		  // of 1/2 in Matrix::TrLessAntiHermMatrix used in Fasqtad::update_momenta.
		  sum[d][n] += m;
		}
	      }else{
		s[mu] -= hop;
		int np = Lat->FsiteOffset(s);
		s[mu] += hop;
		
		for(int v=0; v<n_vect; v++){
		  m.Cross2(vect2[v][n], vect[v][d][np]);
		  // see comment above.
		  sum[d][n] += m;
		}
	      }
	    }
	    else {
	      
	      if(s[mu]+hop>node_sites[mu]-1){
		
		int save_coord = s[mu];
		s[mu] = (s[mu]+hop)%node_sites[mu];
		int np = Lat->FsiteOffset(s);
		s[mu] = save_coord;
		
		for(int v=0; v<n_vect; v++){
		  getPlusData((IFloat*)&vp, (IFloat*)&vect[v][d][np],
			      Lat->FsiteSize(), mu);
		  m.Cross2(vect2[v][n], vp);
		  // Note the factor of 2 from Matrix::Cross2; this is corrected by the factor
		  // of 1/2 in Matrix::TrLessAntiHermMatrix used in Fasqtad::update_momenta.
		  sum[d][n] += m;
		}
	      }else{
		s[mu] += hop;
		int np = Lat->FsiteOffset(s);
		s[mu] -= hop;
		
		for(int v=0; v<n_vect; v++){
		  m.Cross2(vect2[v][n], vect[v][d][np]);
		  // see comment above.
		  sum[d][n] += m;
		}
	      }
	    }
	  }
	}
  
  ParTrans::PTflops += 90*n_vect*n_dir*vol;
}

/*! 
  Computes sum[x] = vect[x] vect[x + hop dir]^dagger
  where the sum is over n_vect vectors and the hop is in a forward direction.
*/
void pt_vvpd(IFloat **vect_v, int n_vect,
	     const int *dir, int n_dir, int hop, IFloat **sum_m){
  
  Vector **vect = (Vector **) vect_v;
  Matrix **sum = (Matrix **) sum_m;
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


#ifndef USE_QMP
void pt_shift_field_vec(IFloat **v_v, const int *dir, int n_dir,
                    int hop, IFloat **u_v){
  char *fname = "pt_shift_field_vec(F**,I*,I,I,F**)";
  char *cname = "PT";

  Vector **v = (Vector **)v_v;
  Vector **u = (Vector **)u_v;
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
		getMinusData((IFloat*)&u[d][n], (IFloat*)&v[d][np], VECT_LEN, mu);
		
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
		getPlusData((IFloat*)&u[d][n], (IFloat*)&v[d][np], VECT_LEN, mu);
		
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



//! u[x] = v[x+dir] for n_dir forward or backward directions dir.
void pt_shift_field(IFloat **v_m, const int *dir, int n_dir,
		    int hop, IFloat **u_m){
    

    Matrix **v = (Matrix **)v_m;
    Matrix **u = (Matrix **)u_m;
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
void pt_shift_link(IFloat **u_m, const int *dir, int n_dir){

    Matrix **u = (Matrix **)u_m;
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

//----------------------------------------------------------------------------

//Parallel transport of a matrix defined on one half of the checkerboarded
//lattice
//
//Parameters
//
//n - The number of direction in which to perform the parallel transport
//mout - Result of the parallel transport, on sites with opposite parity of min
//min - Initial field, defined on sites with only one parity
//dir - a list of the n directions in which the field will be transported
//cb - Checkerboard parity of the vector min

void pt_mat_cb(int n, IFloat **mout, IFloat **min, const int *dir, ChkbType cb)
{
  pt_mat_norm(n,mout,min,dir,cb,(IFloat *)Lat->GaugeField());
}

void pt_mat_cb(int n, IFloat **mout, IFloat **min, const int *dir, ChkbType cb, IFloat * new_gauge_field)
{
  pt_mat_norm(n,mout,min,dir,cb,new_gauge_field);
}

void pt_mat_norm(int n, IFloat **mout, IFloat **min, const int *dir, ChkbType cb, IFloat * gauge)
{
  int s[4], spm[4];
  Matrix * vin[NUM_DIR];
  Matrix * vout[NUM_DIR];
  Matrix * udag;
  IFloat * u, *v;
  int parity = (int) cb;

  //Temporaries to store transferred data
  Matrix m, mp;
  
  for(int i = 0; i < n; i++)
    {
        vin[i] = (Matrix *)min[i];
        vout[i] = (Matrix *)mout[i];
	bool backwards = dir[i]%2 ? true : false;
	int mu = dir[i]/2;

	//Need to transport fields in the forward direction
	if (backwards) {

	for(spm[3]=s[3]=0; s[3]<GJP.TnodeSites(); spm[3]++, s[3]++)
	    for(spm[2]=s[2]=0; s[2]<GJP.ZnodeSites(); spm[2]++, s[2]++)
		for(spm[1]=s[1]=0; s[1]<GJP.YnodeSites(); spm[1]++, s[1]++)
		    for(spm[0]=s[0]=0; s[0]<GJP.XnodeSites(); spm[0]++, s[0]++)
		      {
			if(parity == (s[0]+s[1]+s[2]+s[3])%2)
			  {
			    //If we are going off-node for the gauge field
			    if (s[mu] == 0) {
			      spm[mu] = GJP.NodeSites(mu)-1;
			      //Get the gauge field and the matrix field
			      getMinusData((IFloat*)&m, (IFloat*)(gauge+(Lat->GsiteOffset(spm)+mu)*MATRIX_SIZE), MATRIX_SIZE, mu);
			      getMinusData((IFloat*)&mp, (IFloat*)(&vin[i][Lat->FsiteOffsetChkb(spm)]), MATRIX_SIZE, mu);
			      // STAG stores hermitian conjugate links
			      mDotMEqual((IFloat*)&vout[i][Lat->FsiteOffsetChkb(s)],
					    (IFloat *)&m, (IFloat *)&mp);		    
			    }
			    //Else we don't need to transfer anything
			    else {
			      spm[mu] = s[mu]-1;
			      u = gauge+(Lat->GsiteOffset(spm)+mu)*MATRIX_SIZE;
			      v = (IFloat*)&vin[i][Lat->FsiteOffsetChkb(spm)];
			      mDotMEqual((IFloat*)&vout[i][Lat->FsiteOffsetChkb(s)],u,v);
			    }
			  }
			spm[mu] = s[mu];  
		      }
	}
	//Transport fields in the "negative" direction
	else {

	for(spm[3]=s[3]=0; s[3]<GJP.TnodeSites(); spm[3]++, s[3]++)
	    for(spm[2]=s[2]=0; s[2]<GJP.ZnodeSites(); spm[2]++, s[2]++)
		for(spm[1]=s[1]=0; s[1]<GJP.YnodeSites(); spm[1]++, s[1]++)
		    for(spm[0]=s[0]=0; s[0]<GJP.XnodeSites(); spm[0]++, s[0]++)
		      {
			if(parity == (s[0]+s[1]+s[2]+s[3])%2)
			  {
			    if (s[mu] == GJP.NodeSites(mu)-1) {
			      spm[mu] = 0;
			      //Only get the matrix field, the gauge field is local
			      getPlusData((IFloat*)&mp, (IFloat*)(&vin[i][Lat->FsiteOffsetChkb(spm)]), MATRIX_SIZE, mu);
			      u = gauge+(Lat->GsiteOffset(s)+mu)*MATRIX_SIZE;
			      m.Dagger(u);
			      // STAG stores hermitian conjugate links
			      mDotMEqual((IFloat*)&vout[i][Lat->FsiteOffsetChkb(s)],
					    (IFloat *)&m, (IFloat *)&mp);		 
			    }
			    //We don't need to transfer anything
			    else {
			      spm[mu] = s[mu]+1;
			      u = gauge+(Lat->GsiteOffset(s)+mu)*MATRIX_SIZE;
			      m.Dagger(u);
			      v = (IFloat*)&vin[i][Lat->FsiteOffsetChkb(spm)];
			      mDotMEqual((IFloat*)&vout[i][Lat->FsiteOffsetChkb(s)],(IFloat *)&m,v);
			    }
			  }
			spm[mu] = s[mu];  
		      }
	}
    }
}

//Parallel transport of a vector defined on single parity sites
//
//Parameters
//
//n - number of directions in which to parallel transport
//vout - Transported vector
//vin - Initial vector
//dir - list of directions in which to transport the vectors
//cb - Parity of the sites where the vectors are defined
//gauge - Pointer to block of gauge fields in STAG order

//Normal parallel transport with normal gauge fields
#undef PROFILE
void pt_1vec_cb(int n, IFloat **vout, IFloat **vin, const int *dir, ChkbType cb)
{
  pt_1vec_cb_norm(n,vout,vin,dir,cb,(IFloat *)Lat->GaugeField());
}

//Normal parallel transport, but with user-specified gauge fields
#undef PROFILE
void pt_1vec_cb(int n, IFloat **vout, IFloat **vin, const int *dir, ChkbType cb, IFloat * new_gauge_field)
{
  pt_1vec_cb_norm(n,vout,vin,dir,cb,new_gauge_field);
}

//Padded parallel transport with normal gauge fields
#undef PROFILE
void pt_1vec_cb(int n, IFloat *vout, IFloat **vin, const int *dir, ChkbType cb, int pad)
{
  pt_1vec_cb_pad(n,vout,vin,dir,cb,pad,(IFloat *)Lat->GaugeField());
}

//Padded parallel transport, but with user-specified gauge fields
#undef PROFILE
void pt_1vec_cb(int n, IFloat *vout, IFloat **vin, const int *dir, ChkbType cb, int pad, IFloat * new_gauge_field)
{
  pt_1vec_cb_pad(n,vout,vin,dir,cb,pad,new_gauge_field);
}


void pt_1vec_cb_norm(int n, IFloat **fout, IFloat **fin, const int *dir,ChkbType cb, IFloat * gauge)
{
  int s[4], spm[4];
  Vector * vin[NUM_DIR];
  Vector * vout[NUM_DIR];
  IFloat * u, *v;
  int parity = (int) cb;

  //Temporaries to store transferred data
  Matrix m;
  Vector vp;

  for(int i = 0; i < n; i++)
    {
        vin[i] = (Vector *)fin[i];
        vout[i] = (Vector *)fout[i];
	bool backwards = dir[i]%2 ? true : false;
	int mu = dir[i]/2;

	//Need to transport fields in the forward direction
	if (backwards) {

	for(spm[3]=s[3]=0; s[3]<GJP.TnodeSites(); spm[3]++, s[3]++)
	    for(spm[2]=s[2]=0; s[2]<GJP.ZnodeSites(); spm[2]++, s[2]++)
		for(spm[1]=s[1]=0; s[1]<GJP.YnodeSites(); spm[1]++, s[1]++)
		    for(spm[0]=s[0]=0; s[0]<GJP.XnodeSites(); spm[0]++, s[0]++)
		      {
			if(parity != (s[0]+s[1]+s[2]+s[3])%2)
			  {
			    //If we are going off-node for the gauge field
			    if (s[mu] == 0) {
			      spm[mu] = GJP.NodeSites(mu)-1;
			      //Get the gauge field and the fermion field
			      getMinusData((IFloat*)&m, (IFloat*)(gauge+(Lat->GsiteOffset(spm)+mu)*MATRIX_SIZE), MATRIX_SIZE, mu);
			      getMinusData((IFloat*)&vp, (IFloat*)(&vin[i][Lat->FsiteOffsetChkb(spm)]), VECT_LEN, mu);
			      // STAG stores hermitian conjugate links
			      uDotXEqual((IFloat*)&vout[i][Lat->FsiteOffsetChkb(s)],
					 (IFloat *)&m, (IFloat *)&vp);		    
			    }
			    //Else we don't need to transfer anything
			    else {
			      spm[mu] = s[mu]-1;
			      u = gauge+(Lat->GsiteOffset(spm)+mu)*MATRIX_SIZE;
			      v = (IFloat*)&vin[i][Lat->FsiteOffsetChkb(spm)];
			      uDotXEqual((IFloat*)&vout[i][Lat->FsiteOffsetChkb(s)],u,v);
			    }
			  }
			spm[mu] = s[mu];  
		      }
	}
	//Transport fields in the "negative" direction
	else {

	for(spm[3]=s[3]=0; s[3]<GJP.TnodeSites(); spm[3]++, s[3]++)
	    for(spm[2]=s[2]=0; s[2]<GJP.ZnodeSites(); spm[2]++, s[2]++)
		for(spm[1]=s[1]=0; s[1]<GJP.YnodeSites(); spm[1]++, s[1]++)
		    for(spm[0]=s[0]=0; s[0]<GJP.XnodeSites(); spm[0]++, s[0]++)
		      {
			if(parity != (s[0]+s[1]+s[2]+s[3])%2)
			  {
			    if (s[mu] == GJP.NodeSites(mu)-1) {
			      spm[mu] = 0;
			      //Only get the fermion field, the gauge field is local
			      getPlusData((IFloat*)&vp, (IFloat*)(&vin[i][Lat->FsiteOffsetChkb(spm)]), VECT_LEN, mu);
			      u = gauge+(Lat->GsiteOffset(s)+mu)*MATRIX_SIZE;
			      // STAG stores hermitian conjugate links
			      uDagDotXEqual((IFloat*)&vout[i][Lat->FsiteOffsetChkb(s)],
					    u, (IFloat *)&vp);		 
			    }
			    //We don't need to transfer anything
			    else {
			      spm[mu] = s[mu]+1;
			      u = gauge+(Lat->GsiteOffset(s)+mu)*MATRIX_SIZE;
			      v = (IFloat*)&vin[i][Lat->FsiteOffsetChkb(spm)];
			      uDagDotXEqual((IFloat*)&vout[i][Lat->FsiteOffsetChkb(s)],u,v);
			    }
			  }
			spm[mu] = s[mu];  
		      }
	}

    }
}


void pt_1vec_cb_pad(int n, IFloat *fout, IFloat **fin, const int *dir,ChkbType cb,int pad, IFloat * gauge)
{
  int s[4], spm[4];
  Vector * vin[NUM_DIR];
  IFloat * vout = fout;
  IFloat * u, *v;
  int parity = (int) cb;

  //Temporaries to store transferred data
  Matrix m;
  Vector vp;

  for(int i = 0; i < n; i++)
    {
        vin[i] = (Vector *)fin[i];
	bool backwards = dir[i]%2 ? true : false;
	int mu = dir[i]/2;

	//Need to transport fields in the forward direction
	if (backwards) {

	for(spm[3]=s[3]=0; s[3]<GJP.TnodeSites(); spm[3]++, s[3]++)
	    for(spm[2]=s[2]=0; s[2]<GJP.ZnodeSites(); spm[2]++, s[2]++)
		for(spm[1]=s[1]=0; s[1]<GJP.YnodeSites(); spm[1]++, s[1]++)
		    for(spm[0]=s[0]=0; s[0]<GJP.XnodeSites(); spm[0]++, s[0]++)
		      {
			if(parity != (s[0]+s[1]+s[2]+s[3])%2)
			  {
			    //If we are going off-node for the gauge field
			    if (s[mu] == 0) {
			      spm[mu] = GJP.NodeSites(mu)-1;
			      //Get the gauge field and the fermion field
			      getMinusData((IFloat*)&m, (IFloat*)(gauge+(Lat->GsiteOffset(spm)+mu)*MATRIX_SIZE), MATRIX_SIZE, mu);
			      getMinusData((IFloat*)&vp, (IFloat*)(&vin[i][Lat->FsiteOffsetChkb(spm)]), VECT_LEN, mu);
			      // STAG stores hermitian conjugate links
			      uDotXEqual((IFloat*)&vout[8*(Lat->FsiteOffsetChkb(s)*8+i)],
					 (IFloat *)&m, (IFloat *)&vp);		    
			    }
			    //Else we don't need to transfer anything
			    else {
			      spm[mu] = s[mu]-1;
			      u = gauge+(Lat->GsiteOffset(spm)+mu)*MATRIX_SIZE;
			      v = (IFloat*)&vin[i][Lat->FsiteOffsetChkb(spm)];
			      uDotXEqual((IFloat*)&vout[8*(Lat->FsiteOffsetChkb(s)*8+i)],u,v);
			    }
			  }
			spm[mu] = s[mu];  
		      }
	}
	//Transport fields in the "negative" direction
	else {

	for(spm[3]=s[3]=0; s[3]<GJP.TnodeSites(); spm[3]++, s[3]++)
	    for(spm[2]=s[2]=0; s[2]<GJP.ZnodeSites(); spm[2]++, s[2]++)
		for(spm[1]=s[1]=0; s[1]<GJP.YnodeSites(); spm[1]++, s[1]++)
		    for(spm[0]=s[0]=0; s[0]<GJP.XnodeSites(); spm[0]++, s[0]++)
		      {
			if(parity != (s[0]+s[1]+s[2]+s[3])%2)
			  {
			    if (s[mu] == GJP.NodeSites(mu)-1) {
			      spm[mu] = 0;
			      //Only get the fermion field, the gauge field is local
			      getPlusData((IFloat*)&vp, (IFloat*)(&vin[i][Lat->FsiteOffsetChkb(spm)]), VECT_LEN, mu);
			      u = gauge+(Lat->GsiteOffset(s)+mu)*MATRIX_SIZE;
			      // STAG stores hermitian conjugate links
			      uDagDotXEqual((IFloat*)&vout[8*(Lat->FsiteOffsetChkb(s)*8+i)],
					    u, (IFloat *)&vp);	
			    }
			    //We don't need to transfer anything
			    else {
			      spm[mu] = s[mu]+1;
			      u = gauge+(Lat->GsiteOffset(s)+mu)*MATRIX_SIZE;
			      v = (IFloat*)&vin[i][Lat->FsiteOffsetChkb(spm)];
			      uDagDotXEqual((IFloat*)&vout[8*(Lat->FsiteOffsetChkb(s)*8+i)],u,v);
			    }
			  }
			spm[mu] = s[mu];  
		      }
	}

    }
}

//-----------------------------------------------------------------------------
#endif //USE_QMP

CPS_END_NAMESPACE
	
