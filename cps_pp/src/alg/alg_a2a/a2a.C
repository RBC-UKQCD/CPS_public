#ifdef USE_MPI
//mpi headers
#warning "WARNING : USING MPI"
#include<mpi.h>
#endif

#include <util/verbose.h>
#include <util/vector.h>
#include <util/lattice.h>
#include <alg/lanc_arg.h>
#include <alg/a2a_arg.h>
#include <vector>
#include <fftw3.h>

#ifdef USE_BFM
#include <alg/eigen/Krylov_5d.h>
#endif

#include <alg/a2a/a2a_params.h>


CPS_START_NAMESPACE

A2Aparams::A2Aparams(const A2AArg &_args): args(_args){
  if(args.src_width <= 0) {
    ERR.General("A2Aparams","A2Aparams", "Invalid number for source width (value = %d).\n", args.src_width);
  }
  if(GJP.TnodeSites() % args.src_width != 0) {
    ERR.General("A2Aparams","A2Aparams", "Local t size(%d) is not a multiple of source width(%d).\n", GJP.TnodeSites(), args.src_width);
  }
  nspincolor = 12;
  nflavors = (GJP.Gparity() ? 2:1);
  ntblocks = GJP.Tnodes()*GJP.TnodeSites()/ args.src_width;

  ndilute =  ntblocks * nspincolor* nflavors;      
  nhits = args.nhits;
  nh = nhits * ndilute;

  nl = args.nl;
    
  nv = nl + nh;
}


CPS_END_NAMESPACE
