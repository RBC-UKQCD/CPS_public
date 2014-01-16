#include<config.h>
CPS_START_NAMESPACE
//------------------------------------------------------------------
/*!\file
  \brief Definitions of the AlgLanczos class.

  $Id: alg_lanczos.h,v 1.2 2013-04-05 17:46:30 chulwoo Exp $

*/
//------------------------------------------------------------------


#ifndef INCLUDED_ALG_LANCZOS_H
#define INCLUDED_ALG_LANCZOS_H           //!< Prevent multiple inclusion.

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/smalloc.h>
#include <util/pmalloc.h>
#include <alg/alg_base.h>
#include <alg/common_arg.h>
#include <alg/lanczos_arg.h>
#include <util/eigen_container.h>
CPS_START_NAMESPACE

//! A class implementing a Lanczos eigenvalue solver for the fermion matrix.
/*! \ingroup alg */
class AlgLanczos : public Alg
{
 private:
  char *cname;
  
  LanczosArg *alg_lanczos_arg;
  // The argument structure for the eig algorithm
  
  int Ncb;       
  // Number of checkerboards for fermion field (1 or 2)
  
  Vector **eigenv;
  // The eigenvectors (initial and final)
  
  Float *lambda;
  // The eigenvalues (final)

  EigenCache *ecache;

public:

  AlgLanczos(Lattice & latt, CommonArg *c_arg, LanczosArg *arg, EigenCache* a_cache);

    virtual ~AlgLanczos();

    void run(int init_flag, int ncompress, char* comp_file );
  //    void run(Float **lambda);
  //    void runLanczos();//Float *eval, Vector **evec);

 int NumChkb( RitzMatType ritz_mat)
 {
   // Determine the number of checkerboards
   int Ncb=0;
   char *fname = "NumChkb(RitzType)";
   switch(ritz_mat)
   {
   case MAT_HERM:
   case MATDAG_MAT:
   case NEG_MATDAG_MAT:
   case MATDAG_MAT_NORM:
   case NEG_MATDAG_MAT_NORM:
     Ncb = 2;
     break;
 
   case MATPC_HERM:
   case MATPCDAG_MATPC:
   case MATPCDAG_MATPC_SHIFT:
   case NEG_MATPCDAG_MATPC:
     Ncb = 1;
     break;
 
   default:
     ERR.General("",fname,"RitzMatOper %d not implemented", ritz_mat);
   }
   return Ncb;
 }

};

#endif


CPS_END_NAMESPACE
