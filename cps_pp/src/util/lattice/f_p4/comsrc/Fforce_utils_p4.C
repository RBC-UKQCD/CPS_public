#include<config.h>
//----------------------------------------------------------------------
/*!\file
  \brief  Routines used in the Asqtad RHMC fermion force calculation

  $Id: Fforce_utils_p4.C,v 1.3 2006-03-29 21:21:18 chulwoo Exp $
*/
//----------------------------------------------------------------------


#include <util/lattice.h>
#include <util/gjp.h>
#include <comms/scu.h>
#include<comms/glb.h>

USING_NAMESPACE_CPS


//! The update of the momentum with the force

Float Fp4::update_momenta(Matrix **force, IFloat dt, Matrix *mom) {

  Float Fdt = 0.0;
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
			   (IFloat*)(ip+mu), MATRIX_SIZE);
	    Fdt += dt_tmp*dt_tmp*dotProduct((Float*)&mf, (Float*)&mf, 18);
	  }
	}
  ForceFlops += GJP.VolNodeSites()*54;
  glb_sum(&Fdt);
  return sqrt(Fdt);
}


// The parity of the lattice site n

ChkbType Fp4::parity(const int *n){

  int d, p;
  for(p=0, d=0; d<4; d++) p += n[d];
  return  p%2?CHKB_ODD:CHKB_EVEN;
    
}

