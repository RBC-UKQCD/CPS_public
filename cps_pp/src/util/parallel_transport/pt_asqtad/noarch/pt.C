/*! \file
  \brief  Definition of ParTransAsqtad class methods.
  
  $Id: pt.C,v 1.3 2004-06-04 21:14:14 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-06-04 21:14:14 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/parallel_transport/pt_asqtad/noarch/pt.C,v 1.3 2004-06-04 21:14:14 chulwoo Exp $
//  $Id: pt.C,v 1.3 2004-06-04 21:14:14 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: pt.C,v $
//  $Revision: 1.3 $
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
	
