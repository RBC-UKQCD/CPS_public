//------------------------------------------------------------------
//
// alg_muon.h
//
// T. Blum 8/07
//
// Header file for all alg classes relevant to Wilson-type fermion
// prop' spectrum. The type of glue or fermion is given as
// an argument of type Lattice& to the constructor. If the 
// fermion type is not F_CLASS_WILSON or F_CLASS_CLOVER or 
// F_CLASS_DWF the constructors exit with a general error.
//
//------------------------------------------------------------------


#ifndef INCLUDED_ALG_MUON_H
#define INCLUDED_ALG_MUON_H

#include <alg/alg_base.h>
#include <alg/qpropw.h>
#include <util/lattice.h>
#include <util/fft.h>
#include <alg/muon_arg.h>
#include <alg/eigcg_arg.h>

CPS_START_NAMESPACE

//------------------------------------------------------------------
//
// AlgMuon is derived from Alg and is relevant to  
//
//------------------------------------------------------------------
class AlgMuon : public Alg
{
 private:
    char* cname;

  
    MuonArg* alg_muon_arg;
    // The argument structure for the
    // three point calculation
    QPropWArg* qp_arg;
    // The argument structure for the
    // propagator calculation
    EigCGArg* eigcg_arg;
    // The argument structure for the
    // propagator calculation
    
    QPropW* prop ;
    
    void GetThePropagator(int,Float) ;
    
 public:

    AlgMuon(Lattice & latt, CommonArg* c_arg, MuonArg* arg, QPropWArg* qarg);
    AlgMuon(Lattice & latt, CommonArg* c_arg, MuonArg* arg, QPropWArg* qarg, EigCGArg* ecgarg);

    
    virtual ~AlgMuon();

    void  FFT4(Float* fpin, const int nleft, const int nright,
	       int flag_dist_back=1, int fft_forward=1  )
    {
#ifndef USE_FFTW
      ERR.General(cname,"FFTW4()","Needs FFTW");
#else

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
#endif
    }
    void  FFTXYZ(Float* fpin, const int nleft, const int nright,
		 int flag_dist_back=1, int fft_forward=1  )
    {
#ifndef USE_FFTW
      ERR.General(cname,"FFTWXYZ()","Needs FFTW");
#else
      const int data_size = GJP.VolNodeSites()*nleft*nright;
      //< total data length in the unit of Float
      
      int nr_stride = 1;
      for(int mu=0; mu<3; ++mu) {
	const int n  = GJP.NodeSites(mu);
	const int nr = nright*nr_stride;
	const int nl = data_size/n/nr;
	FFT_one_dir(mu, fpin, nl,n,nr,  fft_forward, flag_dist_back);
	nr_stride *= n;
      }
#endif
    }

    void Mres();
    void MuonLine(int *out_mom, QPropWMomSrc &min, Rcomplex* MuonLine);
    void MuonLineCons(int *out_mom, QPropWMomSrc &min, Rcomplex* MuonLine);
    double photon_prop(int *q,int mu,double momX,double momT);
    void print_vacpol_muonline(void);
    void print_vacpol(void);
    void run(void);
    void run_2(void); //calcs vp and line, and calls run_lbl
    void run_3(void); //calcs vp and line, and calls run_jk_sub_step2(...)
    void run_lbl(int *out_mom,
		 int t_op,
		 Rcomplex *Muon,
		 Rcomplex *vacpol);
    void increment_avg_vacpol(int t_op, Rcomplex *vacpol);
    void increment_avg_Muon(int *out_mom, Rcomplex *vacpol);
    void run_block(void);
    void run_jk_sub(void);
    int run_jk_sub_step1(void);//returns number of confs in sum
    int run_jk_sub_step1_mom(void);//returns number of confs in sum
    void run_jk_sub_step2(int nconfs);
    void run_jk_sub_step2(int*,int,Rcomplex*,Rcomplex*,Rcomplex*,Rcomplex*);
    void run_jk_sub_step2_mom(int nconfs);
    void run_lbl(void);
    void run_print_sub(void);
    void run_sub(void);
    void Pimunu_sub(int t_op);
    void VacPol(int t_op, Rcomplex* VacPol);
    void VacPolConsLoc(int t_op, Rcomplex* VacPol);
    void VacPolConsLocTwistedBC(int t_op, Rcomplex* VacPol);
    void VacPolConsLocTest(int t_op, Rcomplex* VacPol, int* mom);
    void VacPolConsLocLoopPtSrc(int t_op);
    void VacPolConsLocPtSrc(int t_op,int *xsrc);
    void VacPolConsLocLowMode(int t_op, int *mom, Rcomplex* VacPol,
			      Vector* eval, Vector** evec, int Neig);
    void VacPolConsLocLowModeDisc(int t_op, int *mom, Rcomplex* VacPol,
				  Vector* eval, Vector** evec, int Neig);
    void VacPolConsLocRandSrc(int t_op, int *mom, Rcomplex* VacPol);
    void VacPolConsLocQIO(int t_op, Rcomplex* VacPol);
    void two_point(void);
    void two_point(int *mom, QPropWMomSrc &prop, int conj=0);
    void two_pointv2(void);

    void AxialConsLoc(int t_op); //for checking
};

CPS_END_NAMESPACE

#endif
