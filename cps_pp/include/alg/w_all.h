#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-08-18 11:57:36 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/alg/w_all.h,v 1.4 2004-08-18 11:57:36 zs Exp $
//  $Id: w_all.h,v 1.4 2004-08-18 11:57:36 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: w_all.h,v $
//  $Revision: 1.4 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/alg/w_all.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------

//
// w_spect.h:
//      
//
// classes defined in this header file
//   
//   WspectGinfo     --  implemented in w_ginfo.C
//
//     A support class to provide some global information to the rest
//     of this package. It is designed that all the rest classes be
//     derived from this class without any memory or executation
//     overhead.
//
//   WspectHyperRectangle
//                   --  implemented in w_hyper_rect.C
//
//   WspectMomenta   --  implemented in w_momenta.C 
//
//   WspectQuark     --  implemented in w_quark.C
//
//   WspectMesons    --  implemented in w_meson.C
//
//   WspectBaryon    --  implemented in w_baryon.C
//
//
//   added by Thomas and Xiaodong for extended mesons
//   WspectExtendedMesons  -- implemented in w_ext_mesons.C
//   WspectExtendedMesonsBE
//
//   WspectFuzzing
//
// Iff interested in the implementation details [calculation procedures,
// private member functions and quark propagator storage order]:
//          go to the end of this file for more comments.
// 

#ifndef INCLUDED_W_ALL
#define INCLUDED_W_ALL


CPS_END_NAMESPACE
#include <stdio.h>                                // FILE *
#include <util/data_types.h>        // Float, Complex
#include <alg/w_spect_arg.h>        // SourceKind
CPS_START_NAMESPACE
// added by Xiaodong & Thomas
CPS_END_NAMESPACE
#include <alg/w_ginfo.h>
#include <alg/w_ext_mesons.h>
#include <alg/w_hyper_rect.h>
#include <alg/w_ferm_vec.h>
#include <alg/w_quark.h>
#include <alg/w_mesons.h>
#include <alg/w_baryon.h>
#include <alg/w_fuzzing.h>
#include <alg/w_field.h>
#include <alg/w_ext_mesons.h>
#include <alg/w_ext_mesonBE.h>
#include <alg/w_axialcurr.h>
CPS_START_NAMESPACE


//---------------------------------------------------------------------------
// Forward Declarations  -- classes defined in other translation units
//---------------------------------------------------------------------------
class Lattice;                      // defined in util/include/lattice.h
class CgArg;                        // defined in  alg/include/cg_arg.h
class CommonArg;                    // defined in  alg/include/common_arg.h
class AlgWspect;                    // defined in  alg/include/alg_w_spect.h


//---------------------------------------------------------------------------
// Forward Declarations  -- classes defined here, implemented in *.C
//---------------------------------------------------------------------------
class WspectGinfo;                      
class WspectHyperRectangle;                      
class WspectQuark;                     
class WspectMesons;                     
class WspectBaryon;                     
//Added by Xiaodong  and Thomas
class WspectExtendedMesons;



#endif // ! _INCLUDED_W_ALL










//---------------------------------------------------------------------------
// STORAGE ORDER - quark propagator from source y to sink x
//---------------------------------------------------------------------------
//
//    float[y_Dirac][y_Color][x_T][x_Z][x_Y][x_X][x_Dirac][x_Color][Complex]
//
// or briefly
//
//    float[Dy][Cy][T][Z][Y][X][Dx][Cx][2].
//
//
//---------------------------------------------------------------------------
// ARGUMENTs in the private calculation functions of Mesons and Baryons
//---------------------------------------------------------------------------
//
// D1x:  Dirac index for quark propagator 1 at sink   position x
//
// D1y:  Dirac index for quark propagator 1 at source position y
//
// D2x:  Dirac index for quark propagator 2 at sink   position x
//
// D2y:  Dirac index for quark propagator 2 at source position y
//
// D3x:  Dirac index for quark propagator 3 at sink   position x
//
// D3y:  Dirac index for quark propagator 3 at source position y
//
// a_local_site_offset_for_quark_prop:
//       offset in words associated with a local site lcl[] for
//       the access to the quark propagators.
//       It equals to WspectGinfo::siteOffset(lcl) * SPINORs.
//
// result:    RETURNED WITH ACCUMMULATION (!)
//                             ^^^^^^^^^^^^^
//       result += whatever number the [Dirac or Color] algebra gives
//
//
//---------------------------------------------------------------------------
// CALCULATION PROCEDURE
//---------------------------------------------------------------------------
//
// CALCULATION PROCEDURE -- MESONs
//   step 1: color algebra        [same for 16 mesons]
//   step 2: Momentum Projection  [same for 16 mesons]
//   step 3: Dirac Algebra        [different for each meson]
//   finish: Everything()
//
// CALCULATION PROCEDURE -- BARYONs
//   step 1: color algebra
//   step 2: Dirac Algebra
//   step 3: Momentum Projection
//   step 4: Dirac Projection
//   finish: Everything()
//
//
// ColorAlgebra -- all [inner] color indexes disappear thereafter
//     Mesons:    result +=  Trace (QuarkProp1 QuarkProp2^dag)
//     Baryon:    result +=  two color tensors * ....
//
//
// DiracAlgebra -- all  inner  Dirac indexes disappear thereafter
//
//
// MomProject      -- sum over space with right momentum.
//
//---------------------------------------------------------------------------


CPS_END_NAMESPACE
