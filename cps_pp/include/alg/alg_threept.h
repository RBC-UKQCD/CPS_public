#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Definitions of the AlgThreePt class.

  $Id: alg_threept.h,v 1.3 2004-08-18 11:57:35 zs Exp $
*/
//------------------------------------------------------------------

#ifndef INCLUDED_ALG_3PT_H
#define INCLUDED_ALG_3PT_H

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/smalloc.h>
#include <util/pmalloc.h>
#include <alg/alg_base.h>
#include <alg/qpropw.h>
#include <alg/common_arg.h>
#include <alg/threept_arg.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------
//! Class implementing the calculation of meson three point functions.
/*!
  This uses Wilson, clover or domain wall fermions.

  \ingroup alg
*/
//------------------------------------------------------------------
class AlgThreePt : public Alg
{
 private:
    char* cname;

    ThreePtArg* alg_ThreePt_arg;
        // The argument structure for the
        // three point calculation
    int f_size;
        // Node checkerboard size of the fermion field

 public:
    AlgThreePt(Lattice & latt, CommonArg* c_arg, ThreePtArg* arg);

    virtual ~AlgThreePt();

    void run(void);
 
    void figure8(QPropWWallSrc& prop, QPropWWallSrc& prop2);
    void figure8_mix(QPropWWallSrc& prop, QPropWWallSrc& prop2);
    void eye(QPropWWallSrc& prop, QPropWWallSrc& prop2,
        QPropWRandWallSrc& prop3, int t_Op, int t_nt);
    void eye_mix_c4(QPropWWallSrc& prop, QPropWWallSrc& prop2,
        QPropWRandWallSrc& prop3, int t_Op, int t_nt);
    void eye_mix_c31(QPropWWallSrc& prop, QPropWWallSrc& prop2,
        QPropWRandWallSrc& prop3, int t_Op, int t_nt);
    void k_to_vac(QPropWWallSrc& prop, QPropWWallSrc& prop2,
        QPropWRandWallSrc& prop3, int t_Op);
    void k_to_vac_mix_c3(QPropWWallSrc& prop, QPropWWallSrc& prop2,
        QPropWRandWallSrc& prop3, int t_Op);
    void k_to_vac_mix_c21(QPropWWallSrc& prop, QPropWWallSrc& prop2,
        QPropWRandWallSrc& prop3, int t_Op);
    void spectrum(QPropWWallSrc& prop, QPropWWallSrc& prop2);
    void wall_spectrum(QPropWWallSrc& prop, QPropWWallSrc& prop2);

};


#endif

CPS_END_NAMESPACE
