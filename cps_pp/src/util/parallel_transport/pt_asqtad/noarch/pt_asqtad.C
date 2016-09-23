/*! \file
  \brief  Definition of ParTransAsqtad class methods.
  
  $Id: pt_asqtad.C,v 1.4 2004/12/21 19:02:41 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004/12/21 19:02:41 $
//  $Header: /space/cvs/cps/cps++/src/util/parallel_transport/pt_asqtad/noarch/pt_asqtad.C,v 1.4 2004/12/21 19:02:41 chulwoo Exp $
//  $Id: pt_asqtad.C,v 1.4 2004/12/21 19:02:41 chulwoo Exp $
//  $Name: v5_0_16_hantao_io_test_v7 $
//  $Locker:  $
//  $RCSfile: pt_asqtad.C,v $
//  $Revision: 1.4 $
//  $Source: /space/cvs/cps/cps++/src/util/parallel_transport/pt_asqtad/noarch/pt_asqtad.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
#include <config.h>
#include <util/gjp.h>
#include <util/pt.h>

CPS_START_NAMESPACE

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

#if 0
static Lattice  *Lat;
void pt_mat(int N, IFloat **fout, IFloat **fin, const int *dir){

    IFloat *u, *v;
    int s[4], spm[4];
    Matrix ud;
    Matrix *vin[N],*vout[N];

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
    
}

void pt_1vec(int N, IFloat **fout, IFloat **fin, const int *dir){

    IFloat *u, *v;
    int s[4], spm[4];
    Vector *vin[N],*vout[N];


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

}

// Private methods used in the QCDOC implementation.

void pt_init(Lattice &latt) {
  Lat = &latt;

}

void pt_delete(){
}

void pt_delete_g(){
}

void pt_init_g(){
}
#endif

CPS_END_NAMESPACE
	
