#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:46 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_stag/noarch/dirac.C,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
//  $Id: dirac.C,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.4  2001/08/16 10:50:19  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "IFloat".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and IFloat).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.2  2001/06/19 18:12:44  anj
//  Serious ANSIfication.  Plus, degenerate double64.h files removed.
//  Next version will contain the new nga/include/double64.h.  Also,
//  Makefile.gnutests has been modified to work properly, propagating the
//  choice of C++ compiler and flags all the way down the directory tree.
//  The mpi_scu code has been added under phys/nga, and partially
//  plumbed in.
//
//  Everything has newer dates, due to the way in which this first alteration was handled.
//
//  Anj.
//
//  Revision 1.2  2001/05/25 06:16:05  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: dirac.C,v $
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_stag/noarch/dirac.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//-------------------------------------------------------------------
//  dirac.C
//
//  A pure C++ code.
//-------------------------------------------------------------------

CPS_END_NAMESPACE
#include<util/gjp.h>
#include<comms/scu.h>
#include<comms/glb.h>
#include<util/lattice.h>
#include<util/dirac_op.h>
#include<util/vector.h>
CPS_START_NAMESPACE

CPS_END_NAMESPACE
#include<util/stag.h>
CPS_START_NAMESPACE

enum{VECT_LEN=6, MATRIX_SIZE=18, SITE_LEN=72};



//-------------------------------------------------------------------
//  static Vectors and Matrices
//-------------------------------------------------------------------
static Vector vtmp0, vtmp1;
static int nx[4];
static int nb[4];
static int xv[3];
static const Matrix* curU_p;



//-------------------------------------------------------------------
//  initialize dirac
//-------------------------------------------------------------------
void dirac_init(const void * gauge_field_addr)
{
    //-----------------------------------------------------------
    //  nx[4]
    //-----------------------------------------------------------
    nx[0] = GJP.XnodeSites();
    nx[1] = GJP.YnodeSites();
    nx[2] = GJP.ZnodeSites();
    nx[3] = GJP.TnodeSites();

    //-----------------------------------------------------------
    //  nb[4]
    //-----------------------------------------------------------
    nb[0] = 4;
    nb[1] = nb[0]*nx[0];
    nb[2] = nb[1]*nx[1];
    nb[3] = nb[2]*nx[2];

    //-----------------------------------------------------------
    //  xv[3]
    //-----------------------------------------------------------
    xv[0] = nx[3]/2;
    xv[1] = (nx[3]*nx[0])/2;
    xv[2] = (nx[3]*nx[0]*nx[1])/2;

    //-----------------------------------------------------------
    //  pointer to the links
    //-----------------------------------------------------------
    curU_p = (Matrix *)gauge_field_addr;
}


inline int u_offset(const int *x)
{ return nb[0]*x[0]+nb[1]*x[1]+nb[2]*x[2]+nb[3]*x[3]; }

inline int x_offset(const int *x)
{ return (x[3]>>1)+xv[0]*x[0]+xv[1]*x[1]+xv[2]*x[2]; }



//-------------------------------------------------------------------
//  add_flag = 0 :       b = D a
//  add_flag = 1 :       b += D a
//
//  a_odd    = 1 :	 b even;  a odd
//  a_odd    = 0 :	 b odd ;  a even
//
//  D a = \sum_u( U^dag_u(x) a(x+u) - U(x-u) a(x-u) )
//  dir_flag is flag which takes value 0 when all direction contribute to D,
//  1 - when only the special anisotropic direction contributes to D,
//  2 - when all  except the special anisotropic direction. 
//------------------------------------------------------------------

extern "C"
void dirac(IFloat* b, IFloat* a, int a_odd, int add_flag, int dir_flag)
{
  int x[4], nu=GJP.XiDir();
  const Matrix *uoff;
  Vector *boff;
  
  Vector *vp0;
  Vector *vp1;
  const Matrix *mp0;
  
  for(x[0] = 0; x[0] < nx[0]; ++x[0]) {
    for(x[1] = 0; x[1] < nx[1]; ++x[1]) {
      for(x[2] = 0; x[2] < nx[2]; ++x[2]) {
	for(x[3] = 0; x[3] < nx[3]; ++x[3]) {
	  
	  int parity = (x[0]+x[1]+x[2]+x[3])%2;
	  if((a_odd && parity==0) || (a_odd==0 && parity)){
	    
	    uoff = curU_p+u_offset(x);
	    boff = (Vector *)b+x_offset(x);
	    
	    int iter = 0;
	    for (int mu = 0; mu < 4; ++mu) {
	      if ( (dir_flag==0) || 
		   (dir_flag==1 && mu==nu) || 
		   (dir_flag==2 && mu!=nu))  {
		
		//-------------------------------------------
		//  calculate U^dag_u(x) a(x+u)
		//-------------------------------------------
		if(x[mu] == nx[mu]-1) { 	// x+mu off node
		  x[mu] = 0;
		  getPlusData((IFloat *)&vtmp0, a+x_offset(x)*VECT_LEN,
			      VECT_LEN, mu);
		  x[mu] = nx[mu]-1;
		  vp0 = &vtmp0;
		  
		} else { 			// x+mu on node
		  x[mu]++;
		  vp0 = (Vector *)a+x_offset(x);
		  x[mu]--;
		}
		
		mp0 = uoff+mu;
		
		if(iter == 0 && add_flag == 0)
		  uDagDotXEqual((IFloat *)boff, (const IFloat *)mp0,
				(const IFloat *)vp0);
		else
		  uDagDotXPlus((IFloat *)boff,(const IFloat *)mp0,
			       (const IFloat *)vp0);
		
		
		
		//-------------------------------------------
		//  calculate U_u(x-u) a(x-u)
		//-------------------------------------------
		if(x[mu] == 0) { 		// x-mu off node
		  x[mu] = nx[mu]-1;
		  mp0 = uoff+x[mu]*nb[mu]+mu;
		  vp1 = (Vector *)a+x_offset(x);
		  x[mu] = 0;
		  
		  uDotXEqual((IFloat *)&vtmp0, (const IFloat *)mp0,
			     (const IFloat *)vp1);
		  
		  getMinusData((IFloat *)&vtmp1, (IFloat *)&vtmp0,
			       VECT_LEN, mu);
		  
		  *boff -= vtmp1;
		  
		} else { 			// x-mu on node
		  
		  x[mu]--;
		  mp0 = uoff-nb[mu]+mu;
		  vp0 = (Vector *)a+x_offset(x);
		  x[mu]++;
		  
		  uDotXMinus((IFloat *)boff,(const IFloat *)mp0,
			     (const IFloat *)vp0);
		}

		iter++;

	      }
	    }
	  }
	}
      }
    }
  }
}




void destroy_dirac_buf()
{
// do nothing
}

CPS_END_NAMESPACE
