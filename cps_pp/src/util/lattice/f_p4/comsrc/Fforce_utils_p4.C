#include<config.h>
//----------------------------------------------------------------------
/*!\file
  \brief  Routines used in the Asqtad RHMC fermion force calculation

  $Id: Fforce_utils_p4.C,v 1.4 2006-04-13 18:20:37 chulwoo Exp $
*/
//----------------------------------------------------------------------


#include <util/lattice.h>
#include <util/gjp.h>
#include <comms/scu.h>
#include<comms/glb.h>

USING_NAMESPACE_CPS


//! The update of the momentum with the force

ForceArg Fp4::update_momenta(Matrix **force, IFloat dt, Matrix *mom) {

  Float L1 = 0.0;
  Float L2 = 0.0;
  Float Linf = 0.0;
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
	    Float norm = mf.norm();
	    Float tmp = sqrt(norm);
	    L1 += tmp;
	    L2 += norm;
	    Linf = (tmp>Linf ? tmp : Linf);
	  }
	}
  ForceFlops += GJP.VolNodeSites()*54;	
    
  glb_sum(&L1);
  glb_sum(&L2);
  glb_max(&Linf);

  L1 /= 4.0*GJP.VolSites();
  L2 /= 4.0*GJP.VolSites();

  return ForceArg(dt*L1, dt*sqrt(L2), dt*Linf);
 
}


// The parity of the lattice site n

ChkbType Fp4::parity(const int *n){

  int d, p;
  for(p=0, d=0; d<4; d++) p += n[d];
  return  p%2?CHKB_ODD:CHKB_EVEN;
    
}

