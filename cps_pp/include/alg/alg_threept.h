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

#include <alg/common_arg.h>
#include <alg/alg_base.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/smalloc.h>
#include <util/vector.h>
#include <util/verbose.h>
#include <util/error.h>
#include <comms/glb.h>
#include <alg/qpropw.h>
#include <alg/threept_arg.h>


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
 
    void AlgThreePt::figure8(QPropWWallSrc& prop1, QPropWWallSrc& prop2,
			     QPropWWallSrc& prop3, QPropWWallSrc& prop4);
    void AlgThreePt::figure8_mix(QPropWWallSrc& prop1, QPropWWallSrc& prop2,
				 QPropWWallSrc& prop3, QPropWWallSrc& prop4);
    void AlgThreePt::figure8_mix2(QPropWWallSrc& prop1, QPropWWallSrc& prop2,
				  QPropWWallSrc& prop3, QPropWWallSrc& prop4);

    void AlgThreePt::eye(QPropWWallSrc& prop, QPropWWallSrc& prop2,
			 QPropWRandSlabSrc& prop3, WilsonMatrix& spect, 
			 int t_Op, int t_Op_2);
    void AlgThreePt::eye_mix_c4(QPropWWallSrc& prop, QPropWWallSrc& prop2,
				QPropWRandSlabSrc& prop3, WilsonMatrix& spect,
				int t_Op, int t_Op_2);
    void AlgThreePt::eye_mix_c31(QPropWWallSrc& prop, QPropWWallSrc& prop2,
				 QPropWRandSlabSrc& prop3, WilsonMatrix& spect,
				 int t_Op, int t_Op_2);

    void AlgThreePt::k_to_vac(QPropWWallSrc& prop, QPropWWallSrc& prop2,
			      QPropWRandSlabSrc& prop3, 
			      int t_Op, int t_Op_2);
    void AlgThreePt::k_to_vac_mix_c3(QPropWWallSrc& prop, QPropWWallSrc& prop2,
				     QPropWRandSlabSrc& prop3, 
				     int t_Op, int t_Op_2);
    void AlgThreePt::k_to_vac_mix_c21(QPropWWallSrc& prop, 
				      QPropWWallSrc& prop2,
				      QPropWRandSlabSrc& prop3, 
				      int t_Op, int t_Op_2);

    void AlgThreePt::spectrum(QPropWWallSrc& prop, QPropWWallSrc& prop2);
    void AlgThreePt::wall_spectrum(QPropWWallSrc& prop, QPropWWallSrc& prop2);

    // Baryon functions

    Rcomplex AlgThreePt::prop_nucleon(int dd, const WilsonMatrix& p1, 
				      const WilsonMatrix& p2);

    void AlgThreePt::Bspectrum(QPropWWallSrc& prop, QPropWWallSrc& prop2,
			       int t_0);
    void AlgThreePt::Bspectrum_mix(QPropWWallSrc& prop, QPropWWallSrc& prop2, 
				   int t_0);
    void AlgThreePt::wall_Bspectrum(QPropWWallSrc& prop, QPropWWallSrc& prop2,
				    int t_0);
    void AlgThreePt::Axial_connect(QPropWWallSrc& prop, QPropWWallSrc& prop2,
				   int t_sink);
    void AlgThreePt::Vector_connect(QPropWWallSrc& prop, QPropWWallSrc& prop2,
				    int t_sink);
    void AlgThreePt::Scalar_connect(QPropWWallSrc& prop, QPropWWallSrc& prop2, 
				    int t_sink);
    void AlgThreePt::Tensor_connect(QPropWWallSrc& prop, QPropWWallSrc& prop2, 
				    int t_sink);

};


#endif



