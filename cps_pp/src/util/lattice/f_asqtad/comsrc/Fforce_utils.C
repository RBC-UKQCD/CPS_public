#include<config.h>
//----------------------------------------------------------------------
/*!\file
  \brief  Routines used in the Asqtad RHMC fermion force calculation

  $Id: Fforce_utils.C,v 1.2 2004-08-18 11:58:02 zs Exp $
*/
//----------------------------------------------------------------------


#include <util/lattice.h>
#include <util/gjp.h>
#include <comms/scu.h>
#include<comms/glb.h>

USING_NAMESPACE_CPS


//! The update of the momentum with the force

void Fasqtad::update_momenta(Matrix **force, IFloat dt, Matrix *mom) {

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
			if(parity(s) == CHKB_ODD) dt_tmp =-dt;
			else dt_tmp = dt;
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

