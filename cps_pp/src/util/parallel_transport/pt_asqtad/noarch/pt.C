/*! \file
  \brief  Definition of ParTransAsqtad class methods.
  
  $Id: pt.C,v 1.4 2004-08-05 19:00:56 mclark Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mclark $
//  $Date: 2004-08-05 19:00:56 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/parallel_transport/pt_asqtad/noarch/pt.C,v 1.4 2004-08-05 19:00:56 mclark Exp $
//  $Id: pt.C,v 1.4 2004-08-05 19:00:56 mclark Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: pt.C,v $
//  $Revision: 1.4 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/parallel_transport/pt_asqtad/noarch/pt.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
#include <config.h>
#include <util/gjp.h>
#include <util/pt.h>

USING_NAMESPACE_CPS

/*!
  \param latt The lattice object containing the gauge field on which the
  parallel tranport operates.
  \post The gauge field storage order is converted to STAG order if it is not
  already.
*/
ParTransAsqtad::ParTransAsqtad(Lattice &latt): ParTransStagTypes(latt) {

  cname = "ParTransAsqtad";
  char *fname = "ParTransAsqtad(Lattice&)";
  VRB.Func(cname,fname);
  converted = CNV_FRM_NO;
  if(lat.StrOrd()==CANONICAL){
    lat.Convert(STAG);
    converted = CNV_FRM_YES;
  }

}

/*!
  \post The gauge field storage order is converted back to CANONICAL order
  if that was how it was when this object was created.
*/
ParTransAsqtad::~ParTransAsqtad() {
    
  char *fname = "~ParTransAsqtad()";
  VRB.Func(cname,fname);

  // Leave the gauge field as we found it
  if(converted==CNV_FRM_YES) lat.Convert(CANONICAL);  

}


void ParTransAsqtad::run(int N, Matrix **vout, Matrix **vin, const int *dir){

  IFloat *u, *v;
  int s[4], spm[4];
  Matrix ud;

  for(int d=0; d<N; d++){

    int mu = dir[d];
    bool backwards = mu%2 ? true : false;
    mu /= 2;

    for(spm[3]=s[3]=0; s[3]<GJP.TnodeSites(); spm[3]++, s[3]++)
      for(spm[2]=s[2]=0; s[2]<GJP.ZnodeSites(); spm[2]++, s[2]++)
	for(spm[1]=s[1]=0; s[1]<GJP.YnodeSites(); spm[1]++, s[1]++)
	  for(spm[0]=s[0]=0; s[0]<GJP.XnodeSites(); spm[0]++, s[0]++) {
	    if(backwards){
	      spm[mu] = (s[mu]-1+GJP.NodeSites(mu))%GJP.NodeSites(mu);
	      u = (IFloat*)(lat.GaugeField()+lat.GsiteOffset(spm)+mu);
	      v = (IFloat*)&vin[d][lat.FsiteOffset(spm)];
	      mDotMEqual((IFloat*)&vout[d][lat.FsiteOffset(s)],
			 u, v);
	    }else{
	      spm[mu] = (s[mu]+1)%GJP.NodeSites(mu);
	      ud.Dagger(*(lat.GaugeField()+lat.GsiteOffset(s)+mu));
	      u = (IFloat*)&ud;
	      v = (IFloat*)&vin[d][lat.FsiteOffset(spm)];
	      // STAG stores hermitian conjugate links
	      mDotMEqual((IFloat*)&vout[d][lat.FsiteOffset(s)],
			 u, v);			    
	    }
	    spm[mu] = s[mu];

	  }
			
  }
    
}

void ParTransAsqtad::run(int N, Vector **vout, Vector **vin, const int *dir){

  IFloat *u, *v;
  int s[4], spm[4];

  for(int d=0; d<N; d++){

    int mu = dir[d];
    bool backwards = mu%2 ? true : false;
    mu /= 2;

    for(spm[3]=s[3]=0; s[3]<GJP.TnodeSites(); spm[3]++, s[3]++)
      for(spm[2]=s[2]=0; s[2]<GJP.ZnodeSites(); spm[2]++, s[2]++)
	for(spm[1]=s[1]=0; s[1]<GJP.YnodeSites(); spm[1]++, s[1]++)
	  for(spm[0]=s[0]=0; s[0]<GJP.XnodeSites(); spm[0]++, s[0]++) {
	    if(backwards){
	      spm[mu] = (s[mu]-1+GJP.NodeSites(mu))%GJP.NodeSites(mu);
	      u = (IFloat*)(lat.GaugeField()+lat.GsiteOffset(spm)+mu);
	      v = (IFloat*)&vin[d][lat.FsiteOffset(spm)];
	      uDotXEqual((IFloat*)&vout[d][lat.FsiteOffset(s)],
			 u, v);			    
	    }else{
	      spm[mu] = (s[mu]+1)%GJP.NodeSites(mu);
	      u = (IFloat*)(lat.GaugeField()+lat.GsiteOffset(s)+mu);
	      v = (IFloat*)&vin[d][lat.FsiteOffset(spm)];
	      // STAG stores hermitian conjugate links
	      uDagDotXEqual((IFloat*)&vout[d][lat.FsiteOffset(s)],
			    u, v);			    
	    }
	    spm[mu] = s[mu];

	  }
			
  }

}

/*! 
  Computes sum[x] = vect[x] vect[x + hop dir]^dagger
  where the sum is over n_vect vectors and the hop is in a forward direction.
*/
void ParTransAsqtad::vvpd(Vector **vect, int n_vect,
			  const int *dir, int n_dir, int hop, Matrix **sum){
  
  int s[4];
  Matrix m;
  
  Vector vp;

  for(s[3] = 0; s[3]<node_sites[3]; s[3]++)
    for(s[2] = 0; s[2]<node_sites[2]; s[2]++)
      for(s[1] = 0; s[1]<node_sites[1]; s[1]++)
	for(s[0] = 0; s[0]<node_sites[0]; s[0]++) {
	  int n = lat.FsiteOffset(s);
	  
	  for(int d=0; d<n_dir; d++){ 
	    sum[d][n].ZeroMatrix();  
	    if(s[dir[d]]+hop>node_sites[dir[d]]-1){
	      
	      int save_coord = s[dir[d]];
	      s[dir[d]] = (s[dir[d]]+hop)%node_sites[dir[d]];
	      int np = lat.FsiteOffset(s);
	      s[dir[d]] = save_coord;
	      
	      for(int v=0; v<n_vect; v++){
		getPlusData((IFloat*)&vp, (IFloat*)&vect[v][np],
			    lat.FsiteSize(), dir[d]);
		m.Cross2(vect[v][n], vp);
		// Note the factor of 2 from Matrix::Cross2; this is corrected by the factor
		// of 1/2 in Matrix::TrLessAntiHermMatrix used in Fasqtad::update_momenta.
		sum[d][n] += m;
	      }
	    }else{
	      s[dir[d]] += hop;
	      int np = lat.FsiteOffset(s);
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



//! u[x] = v[x+dir] for n_dir forward or backward directions dir.
void ParTransAsqtad::shift_field(Matrix **v, const int *dir, int n_dir,
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
	      int n = lat.FsiteOffset(s);
	      
	      if(s[mu]-hop<0){
		
		int save_coord = s[mu];
		s[mu] -= hop;
		while(s[mu]<0) s[mu] += node_sites[mu];
		int np = lat.FsiteOffset(s);
		s[mu] = save_coord;
		getMinusData((IFloat*)&u[d][n], (IFloat*)&v[d][np], MATRIX_SIZE, mu);
		
	      }else{
		
		s[mu] -= hop;
		int np = lat.FsiteOffset(s);
		s[mu] += hop;
		u[d][n] = v[d][np];
		
	      }
	    }
    } else {
      for(s[3]=0; s[3]<node_sites[3]; s[3]++)
	for(s[2]=0; s[2]<node_sites[2]; s[2]++)
	  for(s[1]=0; s[1]<node_sites[1]; s[1]++)
	    for(s[0]=0; s[0]<node_sites[0]; s[0]++){
	      int n = lat.FsiteOffset(s);
	    
	      if(s[mu]+hop>node_sites[mu]-1){
		
		int save_coord = s[mu];
		s[mu] = (s[mu]+hop)%node_sites[mu];
		int np = lat.FsiteOffset(s);
		s[mu] = save_coord;
		getPlusData((IFloat*)&u[d][n], (IFloat*)&v[d][np], MATRIX_SIZE, mu);
		
	      }else{
		
		s[mu] += hop;
		int np = lat.FsiteOffset(s);
		s[mu] -= hop;
		u[d][n] = v[d][np];
		
	      }
	    }
    }
  }
  
}

//! u[-/+nu](x) = U_[-/+nu](x) 
void ParTransAsqtad::shift_link(Matrix **u, const int *dir, int n_dir){

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
			     (IFloat*)(lat.GaugeField()+ lat.GsiteOffset(s)+ nu), MATRIX_SIZE, nu);
		s[nu] = 0;
	      }else{
		s[nu]--;
		m = *(lat.GaugeField()+lat.GsiteOffset(s)+nu);
		s[nu]++;				
	      }
	      u[d][lat.FsiteOffset(s)] = m;
	    }
    } else {
      
      for(s[3]=0; s[3]<node_sites[3]; s[3]++)
	for(s[2]=0; s[2]<node_sites[2]; s[2]++)
	  for(s[1]=0; s[1]<node_sites[1]; s[1]++)
	    for(s[0]=0; s[0]<node_sites[0]; s[0]++){
	      m.Dagger(*(lat.GaugeField()+lat.GsiteOffset(s)+nu));
	      u[d][lat.FsiteOffset(s)] = m;
	    }
      
    }
    
  }
  
}

// Private methods used in the QCDOC implementation.

void ParTransAsqtad::pt_init(const void *gauge_u) {
    ERR.NotImplemented(cname, "pt_init");
}

void ParTransAsqtad::pt_delete(){
    ERR.NotImplemented(cname, "pt_delete");
}

void ParTransAsqtad::pt_delete_g(){
    ERR.NotImplemented(cname, "pt_delete_g");
}

void ParTransAsqtad::pt_init_g(){
    ERR.NotImplemented(cname, "pt_init_g");
}

int ParTransAsqtad::Offset(int dir, int hop){
    ERR.NotImplemented(cname, "Offset");
}

void ParTransAsqtad::setHopPointer(){
    ERR.NotImplemented(cname, "setHopPointer");
}
	
