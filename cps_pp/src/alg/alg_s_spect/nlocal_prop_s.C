#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2008/02/08 18:35:05 $
//  $Header: /space/cvs/cps/cps++/src/alg/alg_s_spect/nlocal_prop_s.C,v 1.7 2008/02/08 18:35:05 chulwoo Exp $
//  $Id: nlocal_prop_s.C,v 1.7 2008/02/08 18:35:05 chulwoo Exp $
//  $Name: v5_0_16_hantao_io_test_v7 $
//  $Locker:  $
//  $RCSfile: nlocal_prop_s.C,v $
//  $Revision: 1.7 $
//  $Source: /space/cvs/cps/cps++/src/alg/alg_s_spect/nlocal_prop_s.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
// nlocal_prop_s.C
CPS_END_NAMESPACE
#include <alg/nlocal_prop_s.h>
#include <alg/quark_prop_s.h>
#include <util/gjp.h>
#include <util/smalloc.h>
CPS_START_NAMESPACE
		// smalloc()
CPS_END_NAMESPACE
#include <comms/scu.h>	
CPS_START_NAMESPACE
		// getMinusData(), getMinus2Data(), getMinus3Data()
CPS_END_NAMESPACE
#include <alg/myenum.h>
CPS_START_NAMESPACE

#ifdef PARALLEL
CPS_END_NAMESPACE
#include <comms/sysfunc_cps.h>
CPS_START_NAMESPACE
#endif

char NLocalPropS::cname[] = "NLocalPropS";

//-----------------------------------------------------------------------
// local to this file
//-----------------------------------------------------------------------
static int nx[4];

inline int isOffNode(int *s)
{ return (s[0] < 0 || s[1] < 0 || s[2] < 0); }


//-----------------------------------------------------------------------
// get determinant for 3x3 matrix with 3 columns: v0, v1, v2
//-----------------------------------------------------------------------
Complex NLocalPropS::determinant(Complex *v0, Complex *v1, Complex *v2)
{
  return v0[0]*v1[1]*v2[2] + v0[2]*v1[0]*v2[1] + v0[1]*v1[2]*v2[0] -	
	 v0[1]*v1[0]*v2[2] - v0[2]*v1[1]*v2[0] - v0[0]*v1[2]*v2[1]; 	
}


//-----------------------------------------------------------------------
// general term in non-local hadron propagator
//-----------------------------------------------------------------------
Complex NLocalPropS::element(Float **p0, Float **p1, Float **p2) {
   return 
   determinant((Complex *)p0[0], (Complex *)p1[1], (Complex *)p2[2]) +
   determinant((Complex *)p0[2], (Complex *)p1[0], (Complex *)p2[1]) +
   determinant((Complex *)p0[1], (Complex *)p1[2], (Complex *)p2[0]) -
   determinant((Complex *)p0[1], (Complex *)p1[0], (Complex *)p2[2]) -
   determinant((Complex *)p0[2], (Complex *)p1[1], (Complex *)p2[0]) -
   determinant((Complex *)p0[0], (Complex *)p1[2], (Complex *)p2[1]); 
}


//-----------------------------------------------------------------------
//  fetch quark propagators from neighboring nodes in the Minus spatial
//  directions
//-----------------------------------------------------------------------

void NLocalPropS::getNeighbors(int t)
{
  if(isDegQuarks()) {
    transfer(qp0, buffer0, t);
  }
  else {
    transfer(qp0, buffer0, t);
    transfer(qp1, buffer1, t);
    transfer(qp2, buffer2, t);
  }
}

void NLocalPropS::transfer(Float **qp, Float **buffer, int t)
{
  //------------------------------------------------------------
  // * 3 kinds (7 blocks) off-node quark propagators 
  // * to be transfered, the coordinates of origins are 
  // * listed in the array rcv_x[7], and the corresponding 
  // * propagator to send FROM this node have coordinates 
  // * snd_x[7].
  //------------------------------------------------------------

  int rcv_x[7][4] = { 	{ -1, -1, -1, t }, 
		      	{ 0, -1, -1, t  }, 
			{ -1, 0, -1, t  }, 
			{ -1, -1, 0, t  }, 
		      	{ 0, 0, -1, t   }, 
			{ 0, -1, 0, t   }, 
			{ -1, 0, 0, t   }	};
  int snd_x[7][4];

  for (int i = 0; i < 7; ++i)
    for (int j = 0; j < 4; ++j) {
      snd_x[i][j] = (rcv_x[i][j] + nx[j]) % nx[j];
  }

  Float *rcv, *snd;
  int color;

  //-----------------------------------------------------------	
  // * transfer G(-1, -1, -1, t): 3 times SCU
  //-----------------------------------------------------------	
  for (color = 0; color < 3; color++) {
    rcv = buffer[color] + map(rcv_x[0]);
    snd = qp[color] + X_OFFSET(snd_x[0]);

    getMinus3Data((IFloat *)rcv, (IFloat *)snd, VECT_LEN, HDM_T);
  }

  //-----------------------------------------------------------	
  // * transfer from 3 axes: x, y, z
  // * G(x, -1, -1, t ), G(-1, y, -1, t), G(-1, -1, z, t). 
  // * 2 times SCU
  //-----------------------------------------------------------	
  for (int mu = 0; mu < 3; mu++) {
    for (int x = 0; x < nx[mu]; x++) {
      snd_x[mu+1][mu] = rcv_x[mu+1][mu] = x;

      for (color = 0; color < 3; color++) {
        rcv = buffer[color] + map(rcv_x[mu+1]);
        snd = qp[color] + X_OFFSET(snd_x[mu+1]);
        getMinus2Data((IFloat *)rcv, (IFloat *)snd,VECT_LEN,(mu+1)%3,(mu+2)%3);
      }
    }
  } 

  //-----------------------------------------------------------	
  // * transfer from 3 planes, counted as 01(xy), 02(xz) and 12(yz)
  // * 1 time SCU			  2(z),   1(y),      0(x)
  //-----------------------------------------------------------	
  for (int k = 1; k < 3; k++)
    for (int j = 0; j < k; j++) {	// plane j-k

      for (int xk = 0; xk < nx[k]; xk++) { // in k dir
	snd_x[k+j+3][k] = rcv_x[k+j+3][k] = xk;

	for (int xj = 0; xj < nx[j]; xj++) {//in j 
	  snd_x[k+j+3][j] = rcv_x[k+j+3][j] = xj;

          for (color = 0; color < 3; color++) {
            rcv = buffer[color] + map(rcv_x[k+j+3]);
	    snd = qp[color] + X_OFFSET(snd_x[k+j+3]); 
	    
            getMinusData((IFloat *)rcv, (IFloat *)snd, VECT_LEN, 3-k-j );
	  }
        }
      }
  }
}

//-----------------------------------------------------------------------
// map function for the off-node site quark propagator 
// which is stored in the three buffers(by color)
//-----------------------------------------------------------------------
int NLocalPropS::map(int *x)
{
  int offset0 = 1; 
  int offset1 = offset0 + nx[0];
  int offset2 = offset1 + nx[1];
  
  int offset01 = offset2 + nx[2];
  int offset02 = offset01 + nx[0] * nx[1];
  int offset12 = offset02 + nx[0] * nx[2];

  if ((x[0] == -1) && (x[1] == -1) && (x[2] == -1))
     return 0;

  //	 * 1-D 

  else if ((x[1] == -1) && (x[2] == -1)) 
     return (x[0] + offset0)*VECT_LEN;

  else if ((x[0] == -1) && (x[2] == -1)) 
     return (x[1] + offset1)*VECT_LEN;

  else if ((x[0] == -1) && (x[1] == -1)) 
     return (x[2] + offset2)*VECT_LEN;

  //	 * 2-D

  else if (x[2] == -1)
     return (x[0] + x[1]*nx[0]+ offset01)*VECT_LEN;

  else if (x[1] == -1)
     return (x[0] + x[2]*nx[0]+ offset02)*VECT_LEN;

  else if (x[0] == -1)
     return (x[1] + x[2]*nx[1]+ offset12)*VECT_LEN;

  //	 * NOT off-node site

  else 
     return -1;
}

//------------------------------------------------------------------
// local Hypercube sum
//------------------------------------------------------------------

void NLocalPropS::localVal(Complex *currp, int *s)
{
  int color;
  int x[3];
  for (x[0] = -1; x[0] <= 1; x[0] += 2)
    for (x[1] = -1; x[1] <= 1; x[1] += 2)
      for (x[2] = -1; x[2] <= 1; x[2] += 2) {

        int sq[7][4] = {
			{s[0]+x[0], s[1],      s[2],      s[3]},   // x
                       	{s[0],      s[1]+x[1], s[2],      s[3]},   // y
                       	{s[0],      s[1],      s[2]+x[2], s[3]},   // z
                       	{s[0]+x[0], s[1]+x[1], s[2],      s[3]},   // x+y
                       	{s[0]+x[0], s[1],      s[2]+x[2], s[3]},   // x+z
                       	{s[0],      s[1]+x[1], s[2]+x[2], s[3]},   // y+z
                       	{s[0]+x[0], s[1]+x[1], s[2]+x[2], s[3]}};  // x+y+z

        Float *p0[7][3];
        Float *p1[7][3];
        Float *p2[7][3];

        if(isDegQuarks()) {
          for (int i = 0; i < 7; i++) {
            if (isOffNode(sq[i])) {
	      for (color = 0; color < 3; color++) {
                p2[i][color] = p1[i][color] = 
		p0[i][color] = buffer0[color] + map(sq[i]);
    	      }
            }
            else {
              for (color = 0; color < 3; color++) {
                p2[i][color] = p1[i][color] = 
	        p0[i][color] = qp0[color] + X_OFFSET(sq[i]);
              }
            }
          }
	}
        else {
          for (int i = 0; i < 7; i++) {
            if (isOffNode(sq[i])) {
	      for (color = 0; color < 3; color++) {
                p0[i][color] = buffer0[color] + map(sq[i]);
		p1[i][color] = buffer1[color] + map(sq[i]);
		p2[i][color] = buffer2[color] + map(sq[i]);
    	      }
            }
            else {
              for (color = 0; color < 3; color++) {
                p0[i][color] = qp0[color] + X_OFFSET(sq[i]); 
		p1[i][color] = qp1[color] + X_OFFSET(sq[i]);
	        p2[i][color] = qp2[color] + X_OFFSET(sq[i]);
              }
            }
          }
        }

        //-----------------------------------------------------------------
        // calculate the components: 4 is magic number of signs
        //-----------------------------------------------------------------
        Complex c012 =   element(p0[0], p1[1], p2[2]);
        Complex c034 = - element(p0[0], p1[3], p2[4]);
        Complex c135 =   element(p0[1], p1[3], p2[5]);
        Complex c245 = - element(p0[2], p1[4], p2[5]);
 
        //-----------------------------------------------------------------
        // combine the components and construct delta0 , delta1, nucleon16,
        //-----------------------------------------------------------------
        currp[0] += c012;
           // Delta0! factor excluded
        currp[1] += c135 + c034 + c245;
           // Delta1! factor excluded
  
        currp[2] += c034 - c135;
           // Nucleon16 factor excluded
        currp[3] += c034 + c135 - 2*c245;
           // Nucleon16 factor excluded
      }
}

//------------------------------------------------------------------
// CTOR
//------------------------------------------------------------------

NLocalPropS::NLocalPropS(Lattice& lattice, StagNonLocalArg& narg)
: HadronPropS(lattice, 4, narg.dir, QuarkPropSMng::srcSlice(narg.qid0), 2),
  qp0(QuarkPropSMng::prop(narg.qid0)),
  qp1(QuarkPropSMng::prop(narg.qid1)),
  qp2(QuarkPropSMng::prop(narg.qid2))
{
    char *fname = "NLocalPropS(Lattice&, StagNonLocalArg&)";
    nx[0] = GJP.XnodeSites();
    nx[1] = GJP.YnodeSites();
    nx[2] = GJP.ZnodeSites();
    nx[3] = GJP.TnodeSites();

    int i = 0; int j = 1; int k = 2;

    int bufferSize = VECT_LEN * (1 + nx[i] + nx[j] + nx[k] + 
		     nx[i]*nx[j] + nx[i]*nx[k] + nx[j]*nx[k]) ;

    for(int color = 0; color < 3; color++) {
      if(isDegQuarks()) {
        buffer0[color] = (Float *)smalloc(bufferSize*sizeof(Float));
        if(buffer0[color] == 0)
          ERR.Pointer(cname,fname, "buffer0[color]");
        VRB.Smalloc(cname,fname, "buffer0[color]", buffer0[color], 
		  bufferSize*sizeof(Float));
      }
      else {
        buffer0[color] = (Float *)smalloc(bufferSize*sizeof(Float));
        buffer1[color] = (Float *)smalloc(bufferSize*sizeof(Float));
        buffer2[color] = (Float *)smalloc(bufferSize*sizeof(Float));
        if(buffer0[color] == 0 || buffer1[color] == 0 || buffer2[color] == 0)
          ERR.Pointer(cname,fname, "buffer[color]");
        VRB.Smalloc(cname,fname, "buffer0[color]", buffer0[color], 
		  bufferSize*sizeof(Float));
        VRB.Smalloc(cname,fname, "buffer1[color]", buffer1[color], 
		  bufferSize*sizeof(Float));
        VRB.Smalloc(cname,fname, "buffer2[color]", buffer2[color], 
		  bufferSize*sizeof(Float));
      }
    }

}

//------------------------------------------------------------------
// DTOR
//------------------------------------------------------------------
NLocalPropS::~NLocalPropS()
{
  char *fname = "~NLocalPropS()";
  if(isDegQuarks()) {
    for(int color = 0; color < 3; color++) {
      VRB.Sfree(cname,fname, "buffer0[color]",buffer0[color]);
      sfree(buffer0[color]);
    }
  }
  else {
    for(int color = 0; color < 3; color++) {
      VRB.Sfree(cname,fname, "buffer0[color]",buffer0[color]);
      sfree(buffer0[color]);
      VRB.Sfree(cname,fname, "buffer1[color]",buffer1[color]);
      sfree(buffer1[color]);
      VRB.Sfree(cname,fname, "buffer2[color]",buffer2[color]);
      sfree(buffer2[color]);
    }
  }
}


CPS_END_NAMESPACE
