#include <config.h>
#include <util/random.h>
CPS_START_NAMESPACE
//------------------------------------------------------------------
/*!\file
  \brief Wrapper for saving/restoring RNG states

  $Id: rng_state.C,v 1.2 2005-08-02 18:07:43 chulwoo Exp $
*/
//--------------------------------------------------------------------
/*
  $Author: chulwoo $
  $Date: 2005-08-02 18:07:43 $
  $Header: /home/cvs/cps/cps++/src/alg/alg_hmd/alg_hmc_rhmc.C,v 1.18
2005/06/16
07:18:55 chulwoo Exp $
  $Id: rng_state.C,v 1.2 2005-08-02 18:07:43 chulwoo Exp $
  $Name: not supported by cvs2svn $
  $Locker:  $
  $RCSfile: rng_state.C,v $
  $Revision: 1.2 $
  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/random/comsrc/rng_state.C,v $
  $State: Exp $
*/

LRGState::LRGState(){

  cname = "LRGState";
  char *fname = "LRGState()";

  VRB.Func(cname,fname);
  // Allocate memory for the initial rng state
  //----------------------------------------------------------------
  rng4d  = (unsigned int**) smalloc(LRG.NStates(FOUR_D)*sizeof(unsigned int*), 
	cname, fname, "rng4d");
  rng5d  = (unsigned int**) smalloc(LRG.NStates()*sizeof(unsigned int*),
	cname, fname, "rng5d");
  for (int i=0; i<LRG.NStates(FOUR_D); i++)
    rng4d[i] = (unsigned int*) smalloc(LRG.StateSize()*sizeof(unsigned int), 
	cname, fname, "rng4d[i]");
  for (int i=0; i<LRG.NStates(); i++)
    rng5d[i] = (unsigned int*) smalloc(LRG.StateSize()*sizeof(unsigned int), 
	cname, fname, "rng5d[i]");
  
}

LRGState::~LRGState(){

  char *fname = "~LRGState()";
  VRB.Func(cname,fname);
   // Free memory for the initial rng state
    //----------------------------------------------------------------
    for (int i=0; i<LRG.NStates(); i++)
      sfree(rng5d[i], cname, fname, "rng5d[i]");
    for (int i=0; i<LRG.NStates(FOUR_D); i++)
      sfree(rng4d[i], cname, fname, "rng4d[i]");
   
    sfree(rng5d, cname, fname, "rng5d");
    sfree(rng4d, cname, fname, "rng4d");

}

void LRGState::GetStates(){
   char *fname = "GetState()";
   VRB.Func(cname,fname);
   LRG.GetStates(rng5d, FIVE_D);
   LRG.GetStates(rng4d, FOUR_D);
}

void LRGState::SetStates(){
   char *fname = "SetState()";
   VRB.Func(cname,fname);
   LRG.SetStates(rng5d, FIVE_D);
   LRG.SetStates(rng4d, FOUR_D);
}

CPS_END_NAMESPACE
