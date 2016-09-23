#include<config.h>
CPS_START_NAMESPACE
//------------------------------------------------------------------
/*!\file 
  \brief  PAB... Definitions of the AlgMeas class.

  $Id: alg_meas.h,v 1.2 2005/05/21 09:41:33 chulwoo Exp $
*/
//------------------------------------------------------------------


#ifndef INCLUDED_ALG_MEAS_H
#define INCLUDED_ALG_MEAS_H

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/smalloc.h>
#include <util/pmalloc.h>
#include <alg/alg_base.h>
#include <alg/meas_arg.h>
#include <alg/common_arg.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------
//! A class implementing many other algorithms in a scriptable manner
/*!
  Aim is to provide a VML interface for driving stuff without
  recompilation. This will be driven by a main programme that
  loops over configurations for a smooth and seamless interface to
  the CPS valence measurement system.

  \ingroup alg
*/
//------------------------------------------------------------------
    /*
     * Factory for making lattice of various types
     */
class LatticeFactory {
 public:
  static Lattice & Create (FclassType F, GclassType g);
  static void Destroy(void);
  static Lattice *lat_p;
};

class AlgMeas : public Alg
{
 private:
    char *cname;

    // The argument structure containing the task list
    MeasArg *alg_meas_arg;

 public:
    AlgMeas(CommonArg *c_arg, MeasArg *arg);

    virtual ~AlgMeas();


    void run(void);

    void RunTask(MeasTask * Task);
    void Document(char *dir,MeasTask * Task);
    char * Dirname(char *);
    void TruncateFile(char *);
    void TruncateWspectFiles(void);
};



#endif





CPS_END_NAMESPACE
