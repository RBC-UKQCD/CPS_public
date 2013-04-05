#include<config.h>
CPS_START_NAMESPACE

#ifndef INCLUDED_ALG_VACPOLSTAG_H
#define INCLUDED_ALG_VACPOLSTAG_H           //!< Prevent multiple inclusion.

CPS_END_NAMESPACE

#include <util/lattice.h>
#include <util/smalloc.h>
#include <alg/alg_base.h>
#include <alg/common_arg.h>
#include <alg/quark_prop_s.h>
#include <alg/s_spect_arg.h>
#include <util/fft.h>

CPS_START_NAMESPACE

//! A class implementing a Lanczos eigenvalue solver for the fermion matrix.
/*! \ingroup alg */
class AlgVacPolStag : public Alg
{
 private:
    char *cname;

    StagQuarkArg sq_arg;
        // The argument structure for the staggered propagator
    StagQuarkArg sq_arg_twisted;
        // The argument structure for another staggered propagator
    CommonArg svp_common_arg;

 public:

    AlgVacPolStag(Lattice & latt, CommonArg *c_arg, StagQuarkArg *sq_arg);
    AlgVacPolStag(Fp4 & latt, CommonArg *c_arg, StagQuarkArg *sq_arg);
    AlgVacPolStag(Fp4 & latt, CommonArg *c_arg, StagQuarkArg *sq_arg, StagQuarkArg *sq_arg2);

    void VacPolStagConsCons(Rcomplex* VP); // compute the vacuum polarization, FFT'd on sink end
    void VacPolStagConsCons(); // compute the vacuum polarization, FFT'd on sink end
    void VacPolStagConsConsAMA(int Start, int Inc, int tStart, int tInc); // compute the vacuum polarization, FFT'd on src end
    void VacPolStagConsConsTwistedBC(Rcomplex* VacPol);
    void CheckVacPolStagConsCons();
    void CheckVacPolStagConsConsTwistedBC();

    void  FFT4(Float* fpin, const int nleft, const int nright,
               int flag_dist_back=1, int fft_forward=1  )
    {
      const int data_size = GJP.VolNodeSites()*nleft*nright;
      //< total data length in units of Float
      
      int nr_stride = 1;
      for(int mu=0; mu<4; ++mu) {
        const int n  = GJP.NodeSites(mu);
        const int nr = nright*nr_stride;
        const int nl = data_size/n/nr;
        FFT_one_dir(mu, fpin, nl,n,nr,  fft_forward, flag_dist_back);
        nr_stride *= n;
      }
    }

    virtual ~AlgVacPolStag();
};

#endif


CPS_END_NAMESPACE
