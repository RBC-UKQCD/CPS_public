#include<config.h>
CPS_START_NAMESPACE
//------------------------------------------------------------------
//
// alg_threept.h
//
// Header file for all alg classes relevant to Wilson-type fermion
// three point functions. The type of glue or fermion is given as
// an argument of type Lattice& to the constructor. If the 
// fermion type is not F_CLASS_WILSON or F_CLASS_CLOVER or 
// F_CLASS_DWF the constructors exit with a general error.
//
//------------------------------------------------------------------


#ifndef INCLUDED_ALG_3PT_H
#define INCLUDED_ALG_3PT_H

CPS_END_NAMESPACE
#include<util/lattice.h>
#include<util/smalloc.h>
#include<util/pmalloc.h>
#include<alg/alg_base.h>
#include<alg/qpropw.h>
#include<alg/common_arg.h>
#include<alg/threept_arg.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------
//
// AlgThreePt is derived from Alg and is relevant to  
// meson three point functions with Wilson type 
// fermions (i.e. Wilson, Clover, Dwf). 
// The type of fermion is determined by the argument to the 
// constructor. If the fermion type is not F_CLASS_WILSON or 
// F_CLASS_CLOVER or F_CLASS_DWF the constructors exit with a 
// general error.
//
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
 
    void AlgThreePt::figure8(QPropWWallSrc& prop, QPropWWallSrc& prop2);
    void AlgThreePt::figure8_mix(QPropWWallSrc& prop, QPropWWallSrc& prop2);
    void AlgThreePt::eye(QPropWWallSrc& prop, QPropWWallSrc& prop2,
        QPropWRandWallSrc& prop3, int t_Op, int t_nt);
    void AlgThreePt::eye_mix_c4(QPropWWallSrc& prop, QPropWWallSrc& prop2,
        QPropWRandWallSrc& prop3, int t_Op, int t_nt);
    void AlgThreePt::eye_mix_c31(QPropWWallSrc& prop, QPropWWallSrc& prop2,
        QPropWRandWallSrc& prop3, int t_Op, int t_nt);
    void AlgThreePt::k_to_vac(QPropWWallSrc& prop, QPropWWallSrc& prop2,
        QPropWRandWallSrc& prop3, int t_Op);
    void AlgThreePt::k_to_vac_mix_c3(QPropWWallSrc& prop, QPropWWallSrc& prop2,
        QPropWRandWallSrc& prop3, int t_Op);
    void AlgThreePt::k_to_vac_mix_c21(QPropWWallSrc& prop, QPropWWallSrc& prop2,
        QPropWRandWallSrc& prop3, int t_Op);
    void AlgThreePt::spectrum(QPropWWallSrc& prop, QPropWWallSrc& prop2);
    void AlgThreePt::wall_spectrum(QPropWWallSrc& prop, QPropWWallSrc& prop2);

};


#endif
CPS_END_NAMESPACE
