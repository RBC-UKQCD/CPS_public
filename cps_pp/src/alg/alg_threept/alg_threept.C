#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief Definitions of the AlgThreePt class methods.
  
  $Id: alg_threept.C,v 1.5 2004-06-04 21:14:00 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-06-04 21:14:00 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_threept/alg_threept.C,v 1.5 2004-06-04 21:14:00 chulwoo Exp $
//  $Id: alg_threept.C,v 1.5 2004-06-04 21:14:00 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: alg_threept.C,v $
//  $Revision: 1.5 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_threept/alg_threept.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// alg_threept.C
//
// AlgThreePt is derived from Alg and is relevant to  
// three point correlation functions with Wilson-type 
// fermions (i.e. Wilson, Clover, Dwf). 
// The type of fermion is determined by the argument to the 
// constructor but it should not be a non Wilson type lattice.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <stdio.h>
#include <alg/alg_threept.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/smalloc.h>
#include <util/vector.h>
#include <util/verbose.h>
#include <util/error.h>
#include <comms/glb.h>
#include <alg/wilson_matrix.h>
#include <alg/spin_matrix.h>
#ifdef PARALLEL
#include <comms/sysfunc.h>
#endif
CPS_START_NAMESPACE

//------------------------------------------------------------------
/*!
  \param latt The lattice on which to compute the three-point function.
  \param c_arg The common argument structure for all algorithms.
  \param arg The parameters specific to this algorithm.
 */
//------------------------------------------------------------------
AlgThreePt::AlgThreePt(Lattice& latt, 
		     CommonArg *c_arg,
		     ThreePtArg *arg) : 
		     Alg(latt, c_arg) 
{
  cname = "AlgThreePt";
  char *fname = "AlgThreePt(L&,CommonArg*,ThreePtArg*)";
  VRB.Func(cname,fname);

  // Initialize the argument pointer
  //----------------------------------------------------------------
  if(arg == 0)
    ERR.Pointer(cname,fname, "arg");
  alg_ThreePt_arg = arg;


  // Print out input parameters
  //----------------------------------------------------------------
  VRB.Input(cname,fname,
	    "stop_rsd = %g\n",IFloat (alg_ThreePt_arg->cg.stop_rsd));
  VRB.Input(cname,fname,
	    "max_num_iter = %d\n", alg_ThreePt_arg->cg.max_num_iter);
  VRB.Input(cname,fname,
	    "mass = %g\n",IFloat (alg_ThreePt_arg->cg.mass));

  //???
}


//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
AlgThreePt::~AlgThreePt() {
  char *fname = "~AlgThreePt()";
  VRB.Func(cname,fname);

  //???
}


//------------------------------------------------------------------
//
//------------------------------------------------------------------
void AlgThreePt::run()
{
#if TARGET==cpsMPI
    using MPISCU::fprintf;
#endif
  char *fname = "run()";
  VRB.Func(cname,fname);

  // Set the Lattice pointer
  //----------------------------------------------------------------
  Lattice& lat = AlgLattice();

  // get the props
  //----------------------------------------------------------------
  CgArg *cg = &(alg_ThreePt_arg->cg);

  int t_0 = alg_ThreePt_arg->t_src;
	//the src time slice for the u and s quarks
  int t_Op = alg_ThreePt_arg->t_Op;
	//the operator time slice
  int t_Op_2 = alg_ThreePt_arg->t_Op_2;
	//the 2nd operator time slice 
  int t_sink = alg_ThreePt_arg->t_sink;
       // sink time for spectator quark
       // also the src time for the d quark

  GJP.Xbc(BND_CND_PRD);
  GJP.Ybc(BND_CND_PRD);
  GJP.Zbc(BND_CND_PRD);
  GJP.Tbc(BND_CND_PRD);

  // Do the first mass
  cg->mass=alg_ThreePt_arg->mass[0];
  if(common_arg->results != 0){
    FILE *fp;
    if( (fp = fopen((char *)common_arg->results, "a")) == NULL ) {
      ERR.FileA(cname,fname, (char *)common_arg->results);
    }
    fprintf(fp, "MASS= %e\n",cg->mass);
    fclose(fp);
  }

  // Take 0.5*(propp+propm)->prop
  GJP.Tbc(BND_CND_PRD);
  QPropWWallSrc* propp= new QPropWWallSrc( lat, cg, t_0, common_arg );
  GJP.Tbc(BND_CND_APRD);
  QPropWWallSrc* propm= new QPropWWallSrc( lat, cg, t_0, common_arg );
  QPropWWallSrc prop( lat, cg, propp, propm);
  delete propp;
  delete propm;
 
  // Take 0.5*(propp+propm)->prop
  GJP.Tbc(BND_CND_PRD);
  propp= new QPropWWallSrc( lat, cg, t_sink, common_arg );
  GJP.Tbc(BND_CND_APRD);
  propm= new QPropWWallSrc( lat, cg, t_sink, common_arg );
  QPropWWallSrc prop2( lat, cg, propp, propm);
  delete propp;
  delete propm;

  GJP.Tbc(BND_CND_PRD);
  QPropWRandWallSrc prop3( lat, cg, t_Op, alg_ThreePt_arg->seed, common_arg ); 

  // meson 2 and 3 point functions
  //----------------------------------------------------------------
  spectrum(prop, prop);
  spectrum(prop2, prop2);
  wall_spectrum(prop, prop);
  wall_spectrum(prop2, prop2);

  figure8(prop, prop2);
  figure8_mix(prop, prop2);
 
  eye(prop, prop2, prop3, t_Op, t_sink);
  eye_mix_c4(prop, prop2, prop3, t_Op, t_sink);
  eye_mix_c31(prop, prop2, prop3, t_Op, t_sink);

  k_to_vac(prop, prop2, prop3, t_Op);
  k_to_vac_mix_c3(prop, prop2, prop3, t_Op);
  k_to_vac_mix_c21(prop, prop2, prop3, t_Op);

  // do another random source
  GJP.Tbc(BND_CND_PRD);
  QPropWRandWallSrc prop7( lat, cg, t_Op_2, alg_ThreePt_arg->seed, common_arg ); 

  eye(prop, prop2, prop7, t_Op_2, t_sink);
  eye_mix_c4(prop, prop2, prop7, t_Op_2, t_sink);
  eye_mix_c31(prop, prop2, prop7, t_Op_2, t_sink);

  k_to_vac(prop, prop2, prop7, t_Op_2);
  k_to_vac_mix_c3(prop, prop2, prop7, t_Op_2);
  k_to_vac_mix_c21(prop, prop2, prop7, t_Op_2);

  //Do the rest
  for(int m=1; m<alg_ThreePt_arg->num_masses; m++){

    cg->mass=alg_ThreePt_arg->mass[m];
    if(common_arg->results != 0){
      FILE *fp;
      if( (fp = fopen((char *)common_arg->results, "a")) == NULL ) {
        ERR.FileA(cname,fname, (char *)common_arg->results);
      }
      fprintf(fp, "MASS= %e\n",cg->mass);
      fclose(fp);
    }

    // Take 0.5*(propp+propm)->prop
    GJP.Tbc(BND_CND_PRD);
    propp= new QPropWWallSrc( lat, cg, t_0, common_arg );
    GJP.Tbc(BND_CND_APRD);
    propm= new QPropWWallSrc( lat, cg, t_0, common_arg );
    QPropWWallSrc prop4( lat, cg, propp, propm);
    delete propp;
    delete propm;
 
    // Take 0.5*(propp+propm)->prop
    GJP.Tbc(BND_CND_PRD);
    propp= new QPropWWallSrc( lat, cg, t_sink, common_arg );
    GJP.Tbc(BND_CND_APRD);
    propm= new QPropWWallSrc( lat, cg, t_sink, common_arg );
    QPropWWallSrc prop5( lat, cg, propp, propm);
    delete propp;
    delete propm;


    GJP.Tbc(BND_CND_PRD);
    QPropWRandWallSrc prop6( lat, cg, t_Op, alg_ThreePt_arg->seed, common_arg ); 

    //----------------------------------------------------------------

    // meson 2 and 3 point functions
    //----------------------------------------------------------------
    spectrum(prop4, prop4);
    spectrum(prop5, prop5);
    wall_spectrum(prop4, prop4);
    wall_spectrum(prop5, prop5);
	// degenerate
    spectrum(prop, prop4);
    spectrum(prop2, prop5);
    wall_spectrum(prop, prop4);
    wall_spectrum(prop2, prop5);
	// non-degenerate

    figure8(prop4, prop5);
    figure8_mix(prop4, prop5);

    eye(prop4, prop5, prop6, t_Op, t_sink);
    eye_mix_c4(prop4, prop5, prop6, t_Op, t_sink);
    eye_mix_c31(prop4, prop5, prop6, t_Op, t_sink);

    eye(prop, prop2, prop6, t_Op, t_sink);
    eye_mix_c4(prop, prop2, prop6, t_Op, t_sink);
    eye_mix_c31(prop, prop2, prop6, t_Op, t_sink);
    	// charm loop

    k_to_vac(prop, prop4, prop3, t_Op);
    k_to_vac_mix_c3(prop, prop4, prop3, t_Op);
    k_to_vac_mix_c21(prop, prop4, prop3, t_Op);
    k_to_vac(prop, prop4, prop7, t_Op_2);
    k_to_vac_mix_c3(prop, prop4, prop7, t_Op_2);
    k_to_vac_mix_c21(prop, prop4, prop7, t_Op_2);
	// u, d loop
    k_to_vac(prop, prop4, prop6, t_Op);
    k_to_vac_mix_c3(prop, prop4, prop6, t_Op);
    k_to_vac_mix_c21(prop, prop4, prop6, t_Op);
	// s loop

    prop6.Delete();
    GJP.Tbc(BND_CND_PRD);
    QPropWRandWallSrc prop8( lat, cg, t_Op_2, alg_ThreePt_arg->seed, common_arg ); 

    eye(prop4, prop5, prop8, t_Op_2, t_sink);
    eye_mix_c4(prop4, prop5, prop8, t_Op_2, t_sink);
    eye_mix_c31(prop4, prop5, prop8, t_Op_2, t_sink);
 
    eye(prop, prop2, prop8, t_Op_2, t_sink);
    eye_mix_c4(prop, prop2, prop8, t_Op_2, t_sink);
    eye_mix_c31(prop, prop2, prop8, t_Op_2, t_sink);
        // charm loop
    k_to_vac(prop, prop4, prop8, t_Op_2);
    k_to_vac_mix_c3(prop, prop4, prop8, t_Op_2);
    k_to_vac_mix_c21(prop, prop4, prop8, t_Op_2);
        // s loop

    //----------------------------------------------------------------

  } // End masses

}

//
// "figure eight" three-point functions
//----------------------------------------------------------------
void AlgThreePt::figure8(QPropWWallSrc& prop, QPropWWallSrc& prop2)
  {

  cname = "AlgThreePt";
  char *fname = "figure8()";

  int time_size=GJP.Tnodes()*GJP.TnodeSites();
  Rcomplex* aa=(Rcomplex*)smalloc(time_size*sizeof(Rcomplex));
  if(aa == 0)
    ERR.Pointer(cname,fname, "aa");
  VRB.Smalloc(cname,fname,
              "aa",aa, time_size*sizeof(Rcomplex));
	// axial-vector^2
  Rcomplex* vv=(Rcomplex*)smalloc(time_size*sizeof(Rcomplex));
  if(vv == 0)
    ERR.Pointer(cname,fname, "vv");
  VRB.Smalloc(cname,fname,
              "vv",vv, time_size*sizeof(Rcomplex));
	// vector^2
  Rcomplex* aa2=(Rcomplex*)smalloc(time_size*sizeof(Rcomplex));
  if(aa2 == 0)
    ERR.Pointer(cname,fname, "aa2");
  VRB.Smalloc(cname,fname,
              "aa2",aa2, time_size*sizeof(Rcomplex));
  if(aa2 == 0)
    ERR.Pointer(cname,fname, "aa2");
  VRB.Smalloc(cname,fname,
              "aa2",aa2, time_size*sizeof(Rcomplex));
	// axial-vector^2
  Rcomplex* vv2=(Rcomplex*)smalloc(time_size*sizeof(Rcomplex));
  if(vv2 == 0)
    ERR.Pointer(cname,fname, "vv2");
  VRB.Smalloc(cname,fname,
              "vv2",vv2, time_size*sizeof(Rcomplex));
	// vector^2
  WilsonMatrix temp, temp2;

  int t;
  for( t=0; t<time_size; t++) 
	vv[t] = aa[t] = vv2[t] = aa2[t] = 0.0;

  int shift_t=GJP.TnodeCoor()*GJP.TnodeSites();
  int vol= (GJP.VolNodeSites()/GJP.TnodeSites());
  for(int dir=0; dir< 4; dir++){
    for(int i=0; i< GJP.VolNodeSites(); i++){
      temp=prop[i];
      temp.gr(-5).gr(dir);	// mult by anti-quark g_5, then operator gamma
      temp.hconj(); 		// anti-quark conjugation
      temp2=prop2[i];
      temp2.gr(-5).gr(dir);	// mult by anti-quark g_5, then operator gamma
      temp2.hconj();		// anti-quark conjugation
      t=i/vol;
      t+=shift_t;
      vv[t]+=Trace(temp2,prop2[i])*Trace(temp,prop[i]);
      vv2[t]+=Trace(temp2*prop2[i], temp*prop[i]);
      temp.gl(-5);		// operator gamma
      temp2.gl(-5);		// mult by anti-quark g_5, then operator gamma
      aa[t]+=Trace(temp2,prop2[i])*Trace(temp,prop[i]);
      aa2[t]+=Trace(temp2*prop2[i], temp*prop[i]);
    }
  }

  // Global sums and Output the correlators
  for(t=0; t<time_size; t++){
    slice_sum((Float*)&aa[t], 2, 99);
    slice_sum((Float*)&vv[t], 2, 99);
    slice_sum((Float*)&aa2[t], 2, 99);
    slice_sum((Float*)&vv2[t], 2, 99);
  }
  // Print out results
  //----------------------------------------------------------------
  if(common_arg->results != 0){
    FILE *fp;
    if( (fp = fopen((char *)common_arg->results, "a")) == NULL ) {
      ERR.FileA(cname,fname, (char *)common_arg->results);
    }
    for(t=0; t<time_size; t++){
      fprintf(fp," F8TRTR_ %d = %e %e\t%e %e\n", t,
	    vv[t].real(), vv[t].imag(),
    	    aa[t].real(), aa[t].imag());
      fprintf(fp," F8TR_ %d = %e %e\t%e %e\n", t,
   	    vv2[t].real(), vv2[t].imag(),
            aa2[t].real(), aa2[t].imag());
    }
    fclose(fp);
  }
  sfree(vv2);
  sfree(aa2);
  sfree(vv);
  sfree(aa);
}

//
// "figure eight" three-point functions with mixed color indices
//----------------------------------------------------------------
void AlgThreePt::figure8_mix(QPropWWallSrc& prop, QPropWWallSrc& prop2)
  {
#if TARGET==cpsMPI
    using MPISCU::fprintf;
#endif

  cname = "AlgThreePt";
  char *fname = "figure8_mix()";

  int time_size=GJP.Tnodes()*GJP.TnodeSites();
  Rcomplex* aa=(Rcomplex*)smalloc(time_size*sizeof(Rcomplex));
  if(aa == 0)
    ERR.Pointer(cname,fname, "aa");
  VRB.Smalloc(cname,fname,
              "aa",aa, time_size*sizeof(Rcomplex));
        // axial-vector^2
  Rcomplex* vv=(Rcomplex*)smalloc(time_size*sizeof(Rcomplex));
  if(vv == 0)
    ERR.Pointer(cname,fname, "vv");
  VRB.Smalloc(cname,fname,
              "vv",vv, time_size*sizeof(Rcomplex));
        // vector^2
  WilsonMatrix temp, temp2;
  SpinMatrix smat, smat2;
 
  int t;
  for(t=0; t<time_size; t++) vv[t] = aa[t] = 0.0;

  int shift_t=GJP.TnodeCoor()*GJP.TnodeSites();
  int vol= (GJP.VolNodeSites()/GJP.TnodeSites());
  for(int dir=0; dir< 4; dir++){
    for(int i=0; i< GJP.VolNodeSites(); i++){
      temp=prop[i];
      temp.gr(-5).gr(dir);	// mult by anti-quark g_5, then operator gamma
      temp.hconj(); 		// anti-quark conjugation
      temp2=prop2[i];
      temp2.gr(-5).gr(dir);	// mult by anti-quark g_5, then operator gamma
      temp2.hconj();		// anti-quark conjugation
      smat=ColorTrace(temp,prop[i]);
      smat2=ColorTrace(temp2,prop2[i]);
      t=i/vol;
      t+=shift_t;
      vv[t]+=Tr(smat,smat2);
      temp.gl(-5);		// operator gamma
      temp2.gl(-5);		// mult by anti-quark g_5, then operator gamma
      smat=ColorTrace(temp,prop[i]);
      smat2=ColorTrace(temp2,prop2[i]);
      aa[t]+=Tr(smat,smat2);
    }
  }
  // Global sums and Output the correlators
  for(t=0; t<time_size; t++){
    slice_sum((Float*)&aa[t], 2, 99);
    slice_sum((Float*)&vv[t], 2, 99);
  }
  // Print out results
  //----------------------------------------------------------------
  if(common_arg->results != 0){
    FILE *fp;
    if( (fp = fopen((char *)common_arg->results, "a")) == NULL ) {
      ERR.FileA(cname,fname, (char *)common_arg->results);
    }
    for(t=0; t<time_size; t++){
      fprintf(fp,"F8_MX_ %d = %e %e\t%e %e\n", t,
            vv[t].real(), vv[t].imag(),
            aa[t].real(), aa[t].imag());
    }
    fclose(fp);
  }

  sfree(vv);
  sfree(aa);

}


//
// "eye" two and three point functions
//----------------------------------------------------------------
void AlgThreePt::eye(QPropWWallSrc& prop, QPropWWallSrc& prop2,
	QPropWRandWallSrc& prop3, int t_Op, int t_sink)
  { 
#if TARGET==cpsMPI
    using MPISCU::fprintf;
#endif

  cname = "AlgThreePt";
  char *fname = "eye()";

  int time_size=GJP.Tnodes()*GJP.TnodeSites();
  Rcomplex* aa=(Rcomplex*)smalloc(time_size*sizeof(Rcomplex));
  if(aa == 0)
    ERR.Pointer(cname,fname, "aa");
  VRB.Smalloc(cname,fname,
              "aa",aa, time_size*sizeof(Rcomplex));
        // axial-vector^2
  Rcomplex* vv=(Rcomplex*)smalloc(time_size*sizeof(Rcomplex));
  if(vv == 0)
    ERR.Pointer(cname,fname, "vv");
  VRB.Smalloc(cname,fname,
              "vv",vv, time_size*sizeof(Rcomplex));
        // vector^2
  Rcomplex* aa2=(Rcomplex*)smalloc(time_size*sizeof(Rcomplex));
  if(aa2 == 0)
    ERR.Pointer(cname,fname, "aa2");
  VRB.Smalloc(cname,fname,
              "aa2",aa2, time_size*sizeof(Rcomplex));
  if(aa2 == 0)
    ERR.Pointer(cname,fname, "aa2");
  VRB.Smalloc(cname,fname,
              "aa2",aa2, time_size*sizeof(Rcomplex));
        // axial-vector^2
  Rcomplex* vv2=(Rcomplex*)smalloc(time_size*sizeof(Rcomplex));
  if(vv2 == 0)
    ERR.Pointer(cname,fname, "vv2");
  VRB.Smalloc(cname,fname,
              "vv2",vv2, time_size*sizeof(Rcomplex));
        // vector^2
  Rcomplex* scl=(Rcomplex*)smalloc(time_size*sizeof(Rcomplex));
  if(scl == 0)
    ERR.Pointer(cname,fname, "scl");
  VRB.Smalloc(cname,fname,
              "scl",scl, time_size*sizeof(Rcomplex));
  Rcomplex* pbp=(Rcomplex*)smalloc(time_size*sizeof(Rcomplex));
  if(pbp == 0)
    ERR.Pointer(cname,fname, "pbp");
  VRB.Smalloc(cname,fname,
              "pbp",pbp, time_size*sizeof(Rcomplex));
 	// scalar, vacuum saturation

  WilsonMatrix temp, temp2, temp3, temp5;
	// props
  WilsonMatrix temp4 = (Float)0.0;
	// spectator quark

  int t;
  for(t=0; t<time_size; t++) 
     vv[t] = aa[t] = vv2[t] = aa2[t] = 0.;
  for(t=0; t<time_size; t++) 
     scl[t] = pbp[t] = 0.;

  for(int i=0; i< GJP.VolNodeSites(); i++){
	t=i/(GJP.VolNodeSites()/GJP.TnodeSites());
#ifdef PARALLEL
	t += CoorT()*GJP.TnodeSites();
#endif
	if(t != t_sink) continue;
	temp4+=prop[i];
		// u spectator quark 
		// summed over sink time slice
  }

#ifdef PARALLEL
  wilson_matrix wmat = temp4.wmat();
  slice_sum((Float*)&wmat, 288, 99);
	// 99 is to trick it into summing 
	// over all 4 slices
  temp4=wmat;
#endif

  temp4.gr(-5);
	// source gamma_5

#ifdef PARALLEL
  int my_node = CoorT();
#endif
  int sink_node = (int)(t_Op/GJP.TnodeSites());
  int sink_t_on_node = t_Op%GJP.TnodeSites();

  int shift_t=GJP.TnodeCoor()*GJP.TnodeSites();
  int vol= (GJP.VolNodeSites()/GJP.TnodeSites());
  for(int dir=0; dir< 4; dir++){
    for(int i=0; i< GJP.VolNodeSites(); i++){
#ifdef PARALLEL
      if(my_node != sink_node)continue;
#endif
      t=i/vol;
      if(t != sink_t_on_node) continue;
      temp=prop[i];
      temp.gr(-5).hconj();		// s quark x anti-quark gamma_5
      temp2=prop2[i];
      t+=shift_t;
      if(dir==0)scl[t] += Trace(temp * temp4, temp2);
      temp2.gr(dir);			// d quark x Op gamma_mu
      temp3=prop3[i];
      int i_wall = i%vol;
      Rcomplex cc = conj(prop3.rand_src(i_wall));
      if(dir==0)pbp[t] += cc * temp3.Trace();
	  // contract with random source to project out pupil
      temp3.gr(dir);			// pupil quark x Op gamma_mu
      vv[t] += cc * Trace(temp * temp4, temp2) * temp3.Trace();
      vv2[t] += cc * Trace(temp * temp4, temp2 * temp3 );
      temp5 = prop[i];
      temp2.gr(-5);			// d quark x Op gamma_5	
      temp3.gr(-5);			// pupil quark x Op gamma_5
      aa[t] += cc * Trace(temp * temp4, temp2) * temp3.Trace();
      aa2[t] += cc * Trace(temp * temp4, temp2 * temp3 );
    }
  }
  // Global sums and Output the correlators
  for(t=t_Op; t<=t_Op; t++){
    slice_sum((Float*)&aa[t], 2, 99);
    slice_sum((Float*)&vv[t], 2, 99);
    slice_sum((Float*)&aa2[t], 2, 99);
    slice_sum((Float*)&vv2[t], 2, 99);
    slice_sum((Float*)&scl[t], 2, 99);
    slice_sum((Float*)&pbp[t], 2, 99);
  }
  // Print out results
  //----------------------------------------------------------------
  if(common_arg->results != 0){
    FILE *fp;
    if( (fp = fopen((char *)common_arg->results, "a")) == NULL ) {
      ERR.FileA(cname,fname, (char *)common_arg->results);
    }
    for(t=t_Op; t<=t_Op; t++){
      fprintf(fp," I3TRTR_ %d = %e %e\t%e %e\n", t,
            vv[t].real(), vv[t].imag(),
            aa[t].real(), aa[t].imag());
      fprintf(fp," I3TR_ %d = %e %e\t%e %e\n", t,
            vv2[t].real(), vv2[t].imag(),
            aa2[t].real(), aa2[t].imag());
      fprintf(fp," VS3_ %d = %e %e\t%e %e\n", t,
	   scl[t].real(), scl[t].imag(),
	   pbp[t].real(), pbp[t].imag());
    }
    fclose(fp);
  }

  sfree(pbp);
  sfree(scl);
  sfree(vv2);
  sfree(aa2);
  sfree(vv);
  sfree(aa);

}
void AlgThreePt::k_to_vac(QPropWWallSrc& prop, QPropWWallSrc& prop2,
	QPropWRandWallSrc& prop3, int t_Op)
  { 
#if TARGET==cpsMPI
    using MPISCU::fprintf;
#endif

  cname = "AlgThreePt";
  char *fname = "k_to_vac()";

  int time_size=GJP.Tnodes()*GJP.TnodeSites();
  Rcomplex* av=(Rcomplex*)smalloc(time_size*sizeof(Rcomplex));
  if(av == 0)
    ERR.Pointer(cname,fname, "av");
  VRB.Smalloc(cname,fname,
              "av",av, time_size*sizeof(Rcomplex));
  Rcomplex* va=(Rcomplex*)smalloc(time_size*sizeof(Rcomplex));
  if(va == 0)
    ERR.Pointer(cname,fname, "va");
  VRB.Smalloc(cname,fname,
              "va",va, time_size*sizeof(Rcomplex));
  Rcomplex* av2=(Rcomplex*)smalloc(time_size*sizeof(Rcomplex));
  if(av2 == 0)
    ERR.Pointer(cname,fname, "av2");
  VRB.Smalloc(cname,fname,
              "av2",av2, time_size*sizeof(Rcomplex));
  Rcomplex* va2=(Rcomplex*)smalloc(time_size*sizeof(Rcomplex));
  if(va2 == 0)
    ERR.Pointer(cname,fname, "va2");
  VRB.Smalloc(cname,fname,
              "va2",va2, time_size*sizeof(Rcomplex));
  	// k to vac correlators

  WilsonMatrix temp, temp3, temp5;

  int t;
  for(t=0; t<time_size; t++) 
     av2[t] = va2[t] = av[t] = va[t] = 0.;

  int my_node = GJP.TnodeCoor();
  int sink_node = (int)(t_Op/GJP.TnodeSites());
  int sink_t_on_node = t_Op%GJP.TnodeSites();

  int shift_t=GJP.TnodeCoor()*GJP.TnodeSites();
  int vol= (GJP.VolNodeSites()/GJP.TnodeSites());
  for(int dir=0; dir< 4; dir++){
    for(int i=0; i< GJP.VolNodeSites(); i++){
#ifdef PARALLEL
      if(my_node != sink_node)continue;
#endif
      t=i/vol;
      if(t != sink_t_on_node) continue;
      t+=shift_t;
      temp=prop2[i];
      temp.gr(-5).hconj();		// s quark x anti-quark gamma_5
      temp3=prop3[i];
      temp3.gr(dir);			// pupil quark x Op gamma_mu
      temp5 = prop[i];
      temp5.gr(dir).gr(-5);		// d quark for k->vac x Op gamma's
      int i_wall = i%vol;
      Rcomplex cc = conj(prop3.rand_src(i_wall));
	  // contract with random source to project out pupil
      av2[t] += cc * Trace(temp5 * temp3, temp);
      av[t] += cc * Trace(temp , temp5) * temp3.Trace();
      temp3.gr(-5);			// pupil quark x Op gamma_5
      temp5 = prop[i];
      temp5.gr(dir);			// d quark for k->vac x Op gamma's
      va2[t] += cc * Trace(temp * temp5, temp3);
      va[t] += cc * Trace(temp , temp5) * temp3.Trace();
    }
  }
  // Global sums and Output the correlators
  for(t=t_Op; t<=t_Op; t++){
    slice_sum((Float*)&va[t], 2, 99);
    slice_sum((Float*)&av[t], 2, 99);
    slice_sum((Float*)&va2[t], 2, 99);
    slice_sum((Float*)&av2[t], 2, 99);
  }
  // Print out results
  //----------------------------------------------------------------
  if(common_arg->results != 0){
    FILE *fp;
    if( (fp = fopen((char *)common_arg->results, "a")) == NULL ) {
      ERR.FileA(cname,fname, (char *)common_arg->results);
    }
    for(t=t_Op; t<=t_Op; t++){
      fprintf(fp," I2TRTR_ %d = %e %e\t%e %e\n", t,
	   -va[t].real(), -va[t].imag(),
	   -av[t].real(), -av[t].imag());
      fprintf(fp," I2TR_ %d = %e %e\t%e %e\n", t,
	   -va2[t].real(), -va2[t].imag(),
	   -av2[t].real(), -av2[t].imag());
    }
    fclose(fp);
  }

  sfree(va2);
  sfree(av2);
  sfree(va);
  sfree(av);

}

//
// mixed "eye" two and three point functions
// Color 4 prop trace
//----------------------------------------------------------------
void AlgThreePt::eye_mix_c4(QPropWWallSrc& prop, QPropWWallSrc& prop2,
	QPropWRandWallSrc& prop3, int t_Op, int t_sink)
  { 

#if TARGET==cpsMPI
    using MPISCU::fprintf;
#endif
  cname = "AlgThreePt";
  char *fname = "eye_mix_c4()";

  int time_size=GJP.Tnodes()*GJP.TnodeSites();
  Rcomplex* aa=(Rcomplex*)smalloc(time_size*sizeof(Rcomplex));
  if(aa == 0)
    ERR.Pointer(cname,fname, "aa");
  VRB.Smalloc(cname,fname,
              "aa",aa, time_size*sizeof(Rcomplex));
        // axial-vector^2
  Rcomplex* vv=(Rcomplex*)smalloc(time_size*sizeof(Rcomplex));
  if(vv == 0)
    ERR.Pointer(cname,fname, "vv");
  VRB.Smalloc(cname,fname,
              "vv",vv, time_size*sizeof(Rcomplex));
        // vector^2

  WilsonMatrix temp, temp2, temp3;
  WilsonMatrix temp4 = (Float)0.0;
	// spectator quark
  Matrix cmat, cmat2;

  int t;
  for(t=0; t<time_size; t++) 
	vv[t] = aa[t] = 0.;

  for(int i=0; i< GJP.VolNodeSites(); i++){
    t=i/(GJP.VolNodeSites()/GJP.TnodeSites());
#ifdef PARALLEL
    t += CoorT()*GJP.TnodeSites();
#endif
    if(t != t_sink) continue;
    temp4+=prop[i];
  }

#ifdef PARALLEL
  wilson_matrix wmat = temp4.wmat();
  slice_sum((Float*)&wmat, 288, 5);
	// 5 is to trick it into summing 
	// over all 4 slices
  temp4=wmat;
#endif

  temp4.gr(-5);
	// source gamma_5

#ifdef PARALLEL
  int my_node = CoorT();
#endif
  int sink_node = (int)(t_Op/GJP.TnodeSites());
  int sink_t_on_node = t_Op%GJP.TnodeSites();

  int shift_t=GJP.TnodeCoor()*GJP.TnodeSites();
  int vol= (GJP.VolNodeSites()/GJP.TnodeSites());
  for(int dir=0; dir< 4; dir++){
    for(int i=0; i< GJP.VolNodeSites(); i++){
#ifdef PARALLEL
      if(my_node != sink_node)continue;
#endif
      t=i/vol;
      if(t != sink_t_on_node) continue;
      temp=prop[i];
      temp.gr(-5).hconj();		// s quark x anti-quark gamma_5
      temp2=prop2[i];
      temp2.gr(dir);			// d quark x Op gamma_mu
      temp3=prop3[i];
      temp3.gr(dir);			// pupil quark x Op gamma_mu
      int i_wall = i%vol;
      Rcomplex cc = conj(prop3.rand_src(i_wall));
	  // contract with random source to project out pupil
      cmat= SpinTrace(temp, temp4, temp2);
      cmat2= SpinTrace(temp3);
      t+=shift_t;
      vv[t] += cc * Tr(cmat,cmat2);
      temp2.gr(-5);			// d quark x Op gamma_5	
      temp3.gr(-5);			// pupil quark x Op gamma_5
      cmat= SpinTrace(temp, temp4, temp2);
      cmat2= SpinTrace(temp3);
      aa[t] += cc * Tr(cmat,cmat2);
    }
  }

  // Global sums and Output the correlators
  for(t=t_Op; t<=t_Op; t++){
    slice_sum((Float*)&vv[t], 2, 99);
    slice_sum((Float*)&aa[t], 2, 99);
  }
  // Print out results
  //----------------------------------------------------------------
  if(common_arg->results != 0){
    FILE *fp;
    if( (fp = fopen((char *)common_arg->results, "a")) == NULL ) {
      ERR.FileA(cname,fname, (char *)common_arg->results);
    }
    for(t=t_Op; t<=t_Op; t++){
      fprintf(fp,"I3TR_MX_ %d = %e %e\t%e %e\n", t,
            vv[t].real(), vv[t].imag(),
            aa[t].real(), aa[t].imag());
    }
    fclose(fp);
  }
  sfree(vv);
  sfree(aa);

}
void AlgThreePt::k_to_vac_mix_c3(QPropWWallSrc& prop, QPropWWallSrc& prop2,
	QPropWRandWallSrc& prop3, int t_Op)
  { 

#if TARGET==cpsMPI
    using MPISCU::fprintf;
#endif
  cname = "AlgThreePt";
  char *fname = "k_to_vac_mix_c3()";

  int time_size=GJP.Tnodes()*GJP.TnodeSites();
  Rcomplex* av=(Rcomplex*)smalloc(time_size*sizeof(Rcomplex));
  if(av == 0)
    ERR.Pointer(cname,fname, "av");
  VRB.Smalloc(cname,fname,
              "av",av, time_size*sizeof(Rcomplex));
  Rcomplex* va=(Rcomplex*)smalloc(time_size*sizeof(Rcomplex));
  if(va == 0)
    ERR.Pointer(cname,fname, "va");
  VRB.Smalloc(cname,fname,
              "va",va, time_size*sizeof(Rcomplex));
        // k to vac correlators

  WilsonMatrix temp, temp3, temp5;
  Matrix cmat, cmat2;

  int t;
  for(t=0; t<time_size; t++) 
	av[t] = va[t]=0.;

  int my_node = GJP.TnodeCoor();
  int sink_node = (int)(t_Op/GJP.TnodeSites());
  int sink_t_on_node = t_Op%GJP.TnodeSites();

  int shift_t=GJP.TnodeCoor()*GJP.TnodeSites();
  int vol= (GJP.VolNodeSites()/GJP.TnodeSites());
  for(int dir=0; dir< 4; dir++){
    for(int i=0; i< GJP.VolNodeSites(); i++){
#ifdef PARALLEL
      if(my_node != sink_node)continue;
#endif
      t=i/vol;
      if(t != sink_t_on_node) continue;
      temp=prop2[i];
      temp.gr(-5).hconj();		// s quark x anti-quark gamma_5
      temp3=prop3[i];
      temp3.gr(dir);			// pupil quark x Op gamma_mu
      int i_wall = i%vol;
      Rcomplex cc = conj(prop3.rand_src(i_wall));
	  // contract with random source to project out pupil
      cmat2= SpinTrace(temp3);
      temp5 = prop[i];
      temp5.gr(dir).gr(-5);		// d quark for k->vac x Op gamma's
      cmat= SpinTrace(temp, temp5);
      t+=shift_t;
      av[t] += cc * Tr(cmat,cmat2);
      temp3.gr(-5); 		        // pupil quark x Op gamma_5
      cmat2= SpinTrace(temp3);
      temp5 = prop[i];
      temp5.gr(dir);			// d quark for k->vac x Op gamma's
      cmat= SpinTrace(temp, temp5);
      va[t] += cc * Tr(cmat,cmat2);
    }
  }

  // Global sums and Output the correlators
  for(t=t_Op; t<=t_Op; t++){
    slice_sum((Float*)&va[t], 2, 99);
    slice_sum((Float*)&av[t], 2, 99);
  }
  // Print out results
  //----------------------------------------------------------------
  if(common_arg->results != 0){
    FILE *fp;
    if( (fp = fopen((char *)common_arg->results, "a")) == NULL ) {
      ERR.FileA(cname,fname, (char *)common_arg->results);
    }
    for(t=t_Op; t<=t_Op; t++){
      fprintf(fp,"I2TR_MX_ %d = %e %e\t%e %e\n", t,
           -va[t].real(), -va[t].imag(),
           -av[t].real(), -av[t].imag());
    }
    fclose(fp);
  }
  sfree(va);
  sfree(av);

}
//
// mixed "eye" two and three point functions
// Color 3 prop trace times 1 prop trace
//----------------------------------------------------------------
void AlgThreePt::eye_mix_c31(QPropWWallSrc& prop, QPropWWallSrc& prop2,
	QPropWRandWallSrc& prop3, int t_Op, int t_sink)
  { 

#if TARGET==cpsMPI
    using MPISCU::fprintf;
#endif
  cname = "AlgThreePt";
  char *fname = "eye_mix_c31()";

  int time_size=GJP.Tnodes()*GJP.TnodeSites();
  Rcomplex* aa=(Rcomplex*)smalloc(time_size*sizeof(Rcomplex));
  if(aa == 0)
    ERR.Pointer(cname,fname, "aa");
  VRB.Smalloc(cname,fname,
              "aa",aa, time_size*sizeof(Rcomplex));
        // axial-vector^2
  Rcomplex* vv=(Rcomplex*)smalloc(time_size*sizeof(Rcomplex));
  if(vv == 0)
    ERR.Pointer(cname,fname, "vv");
  VRB.Smalloc(cname,fname,
              "vv",vv, time_size*sizeof(Rcomplex));
        // vector^2
  WilsonMatrix temp, temp2, temp3;
  WilsonMatrix temp4 = (Float)0.0;
	// spectator quark
  SpinMatrix smat, smat2;

  int t;
  for(t=0; t<time_size; t++) 
	vv[t] = aa[t] = 0.;

  for(int i=0; i< GJP.VolNodeSites(); i++){
    t=i/(GJP.VolNodeSites()/GJP.TnodeSites());
#ifdef PARALLEL
    t += CoorT()*GJP.TnodeSites();
#endif
    if(t != t_sink) continue;
    temp4+=prop[i];
	// u spectator quark summed over time slice
  }

#ifdef PARALLEL
  wilson_matrix wmat = temp4.wmat();
  slice_sum((Float*)&wmat, 288, 99);
	// 99 is to trick it into summing 
	// over all 4 slices
  temp4=wmat;
#endif

  temp4.gr(-5);
	// source gamma_5

#ifdef PARALLEL
  int my_node = CoorT();
#endif
  int sink_node = (int)(t_Op/GJP.TnodeSites());
  int sink_t_on_node = t_Op%GJP.TnodeSites();

  int shift_t=GJP.TnodeCoor()*GJP.TnodeSites();
  int vol= (GJP.VolNodeSites()/GJP.TnodeSites());
  for(int dir=0; dir< 4; dir++){
    for(int i=0; i< GJP.VolNodeSites(); i++){
#ifdef PARALLEL
      if(my_node != sink_node)continue;
#endif
      t=i/vol;
      if(t != sink_t_on_node) continue;
      temp=prop[i];
      temp.gr(-5).hconj();		// s quark x anti-quark gamma_5
      temp2=prop2[i];
      temp2.gr(dir);			// d quark x Op gamma_mu
      temp3=prop3[i];
      temp3.gr(dir);			// pupil quark x Op gamma_mu
      int i_wall = i%vol;
      Rcomplex cc = conj(prop3.rand_src(i_wall));
	  // contract with random source to project out pupil
      smat= ColorTrace(temp, temp4, temp2);
      smat2= ColorTrace(temp3);
      t+=shift_t;
      vv[t] += cc * Tr(smat,smat2);
      temp2.gr(-5);			// d quark x Op gamma_5	
      temp3.gr(-5);			// pupil quark x Op gamma_5
      smat= ColorTrace(temp, temp4, temp2);
      smat2= ColorTrace(temp3);
      aa[t] += cc * Tr(smat,smat2);
    }
  }
  // Global sums and Output the correlators
  for(t=t_Op; t<=t_Op; t++){
    slice_sum((Float*)&vv[t], 2, 99);
    slice_sum((Float*)&aa[t], 2, 99);
  }
  // Print out results
  //----------------------------------------------------------------
  if(common_arg->results != 0){
    FILE *fp;
    if( (fp = fopen((char *)common_arg->results, "a")) == NULL ) {
      ERR.FileA(cname,fname, (char *)common_arg->results);
    }
    for(t=t_Op; t<=t_Op; t++){
      fprintf(fp,"I3TRTR_MX_ %d = %e %e\t%e %e\n", t,
            vv[t].real(), vv[t].imag(),
            aa[t].real(), aa[t].imag());
    }
    fclose(fp);
  }
  sfree(vv);
  sfree(aa);

}
void AlgThreePt::k_to_vac_mix_c21(QPropWWallSrc& prop, QPropWWallSrc& prop2,
	QPropWRandWallSrc& prop3, int t_Op)
  { 

#if TARGET==cpsMPI
    using MPISCU::fprintf;
#endif
  cname = "AlgThreePt";
  char *fname = "k_to_vac_mix_c21()";

  int time_size=GJP.Tnodes()*GJP.TnodeSites();
  Rcomplex* av=(Rcomplex*)smalloc(time_size*sizeof(Rcomplex));
  if(av == 0)
    ERR.Pointer(cname,fname, "av");
  VRB.Smalloc(cname,fname,
              "av",av, time_size*sizeof(Rcomplex));
  Rcomplex* va=(Rcomplex*)smalloc(time_size*sizeof(Rcomplex));
  if(va == 0)
    ERR.Pointer(cname,fname, "va");
  VRB.Smalloc(cname,fname,
              "va",va, time_size*sizeof(Rcomplex));
        // k to vac correlators
  WilsonMatrix temp, temp3, temp5;
	// props
  SpinMatrix smat, smat2;

  int t;
  for(t=0; t<time_size; t++) 
	av[t] = va[t]=0.;

  int my_node = GJP.TnodeCoor();
  int sink_node = (int)(t_Op/GJP.TnodeSites());
  int sink_t_on_node = t_Op%GJP.TnodeSites();

  int shift_t=GJP.TnodeCoor()*GJP.TnodeSites();
  int vol= (GJP.VolNodeSites()/GJP.TnodeSites());
  for(int dir=0; dir< 4; dir++){
    for(int i=0; i< GJP.VolNodeSites(); i++){
#ifdef PARALLEL
      if(my_node != sink_node)continue;
#endif
      t=i/vol;
      if(t != sink_t_on_node) continue;
      temp=prop2[i];
      temp.gr(-5).hconj();		// s quark x anti-quark gamma_5
      temp3=prop3[i];
      temp3.gr(dir);			// pupil quark x Op gamma_mu
      int i_wall = i%vol;
      Rcomplex cc = conj(prop3.rand_src(i_wall));
	  // contract with random source to project out pupil
      smat2= ColorTrace(temp3);
      t+=shift_t;
      temp5 = prop[i];
      temp5.gr(dir).gr(-5);		// d quark for k->vac x Op gamma's
      smat= ColorTrace(temp, temp5);
      av[t] += cc * Tr(smat,smat2);
      temp3.gr(-5);			// pupil quark x Op gamma_5
      smat2= ColorTrace(temp3);
      temp5 = prop[i];
      temp5.gr(dir);			// d quark for k->vac x Op gamma's
      smat= ColorTrace(temp, temp5);
      va[t] += cc * Tr(smat,smat2);
    }
  }
  // Global sums and Output the correlators
  for(t=t_Op; t<=t_Op; t++){
    slice_sum((Float*)&va[t], 2, 99);
    slice_sum((Float*)&av[t], 2, 99);
  }
  // Print out results
  //----------------------------------------------------------------
  if(common_arg->results != 0){
    FILE *fp;
    if( (fp = fopen((char *)common_arg->results, "a")) == NULL ) {
      ERR.FileA(cname,fname, (char *)common_arg->results);
    }
    for(t=t_Op; t<=t_Op; t++){
      fprintf(fp,"I2TRTR_MX_ %d = %e %e\t%e %e\n", t,
           -va[t].real(), -va[t].imag(),
           -av[t].real(), -av[t].imag());
    }
    fclose(fp);
  }
  sfree(va);
  sfree(av);

}

//
// two-point functions
//------------------------------------------------------------------
void AlgThreePt::spectrum(QPropWWallSrc& prop, QPropWWallSrc& prop2)
  { 

#if TARGET==cpsMPI
    using MPISCU::fprintf;
#endif
  cname = "AlgThreePt";
  char *fname = "spectrum()";

  int time_size=GJP.Tnodes()*GJP.TnodeSites();
  Rcomplex* trace=(Rcomplex*)smalloc(time_size*sizeof(Rcomplex));
  if(trace == 0)
    ERR.Pointer(cname,fname, "trace");
  VRB.Smalloc(cname,fname,
              "trace",trace, time_size*sizeof(Rcomplex));
  Rcomplex* trace2=(Rcomplex*)smalloc(time_size*sizeof(Rcomplex));
  if(trace2 == 0)
    ERR.Pointer(cname,fname, "trace2");
  VRB.Smalloc(cname,fname,
              "trace2",trace2, time_size*sizeof(Rcomplex));
  Rcomplex* trace3=(Rcomplex*)smalloc(time_size*sizeof(Rcomplex));
  if(trace3 == 0)
    ERR.Pointer(cname,fname, "trace3");
  VRB.Smalloc(cname,fname,
              "trace3",trace3, time_size*sizeof(Rcomplex));

  int t;
  for(t=0; t<time_size; t++)
	trace[t] = trace2[t] = trace3[t] = 0.0;

  WilsonMatrix temp;
 
  int shift_t=GJP.TnodeCoor()*GJP.TnodeSites();
  int vol= (GJP.VolNodeSites()/GJP.TnodeSites());
  for(int i=0; i< GJP.VolNodeSites(); i++){
    t=i/vol;
    temp=prop2[i];
    temp.hconj();
	// hermitean conjugate
    t+=shift_t;
    trace[t] += Trace(prop[i], temp);
	// pseudo-scalar correlator
    temp=prop2[i];
    temp.gr(3);
	// t dir
    temp.hconj();
	// hermitean conjugate
    trace2[t] += Trace(prop[i], temp);
	// axial correlator
    temp.gl(2).gl(-5).gr(3).gr(2).gr(-5);
    trace3[t] += Trace(prop[i], temp);
        // vector correlator

  }
  // Global sums and Output the correlators
  for(t=0; t<time_size; t++){
    slice_sum((Float*)&trace[t], 2, 99);
    slice_sum((Float*)&trace2[t], 2, 99);
    slice_sum((Float*)&trace3[t], 2, 99);
  }
  // Print out results
  //----------------------------------------------------------------
  if(common_arg->results != 0){
    FILE *fp;
    if( (fp = fopen((char *)common_arg->results, "a")) == NULL ) {
      ERR.FileA(cname,fname, (char *)common_arg->results);
    }
    for(t=0; t<time_size; t++){
      fprintf(fp, " PSEUDO %d   %e %e %e %e\n", t,
	    trace[t].real(), trace[t].imag(),
	    trace2[t].real(), trace2[t].imag());
      fprintf(fp," RHO3232 %d  %e %e\n", t,
	    trace3[t].real(), trace3[t].imag());
    }
    fclose(fp);
  }
  sfree(trace3);
  sfree(trace2);
  sfree(trace);

}
void AlgThreePt::wall_spectrum(QPropWWallSrc& prop, QPropWWallSrc& prop2)
  { 

#if TARGET==cpsMPI
    using MPISCU::fprintf;
#endif
  cname = "AlgThreePt";
  char *fname = "wall_spectrum()";

  int time_size=GJP.Tnodes()*GJP.TnodeSites();
  Rcomplex* trace=(Rcomplex*)smalloc(time_size*sizeof(Rcomplex));
  if(trace == 0)
    ERR.Pointer(cname,fname, "trace");
  VRB.Smalloc(cname,fname,
              "trace",trace, time_size*sizeof(Rcomplex));
  WilsonMatrix* wprop=(WilsonMatrix*)smalloc(time_size*sizeof(WilsonMatrix));
  if(wprop == 0)
    ERR.Pointer(cname,fname, "wprop");
  VRB.Smalloc(cname,fname,
              "wprop",wprop, time_size*sizeof(Rcomplex));
  WilsonMatrix* wprop2=(WilsonMatrix*)smalloc(time_size*sizeof(WilsonMatrix));
  if(wprop2 == 0)
    ERR.Pointer(cname,fname, "wprop2");
  VRB.Smalloc(cname,fname,
              "wprop2",wprop2, time_size*sizeof(Rcomplex));

  int t;
  for(t=0; t<time_size; t++)
	trace[t] = 0.0;
  for(t=0; t<time_size; t++)
	wprop[t] = wprop2[t] = 0.0;

  // Wall sink (sum on node)
  for(int i=0; i< GJP.VolNodeSites(); i++){
    t=i/(GJP.VolNodeSites()/GJP.TnodeSites());
    t += GJP.TnodeCoor()*GJP.TnodeSites();
    wprop[t]+=prop[i];
    wprop2[t]+=prop2[i];
        // quark, antiquark summed over time slice
  }
 
  // Wall sink (sum over nodes)
  for(t=0; t<time_size; t++){
    wilson_matrix wmat = wprop[t].wmat();
    wilson_matrix wmat2 = wprop2[t].wmat();
    slice_sum((Float*)&wmat, 288, 99);
    slice_sum((Float*)&wmat2, 288, 99);
    WilsonMatrix temp=wmat;
    WilsonMatrix temp2=wmat2;
    temp2.hconj();
	// hermitean conjugate
    trace[t] += Trace(temp2, temp);
	// pseudo-scalar correlator
  }
  // Print out results
  //----------------------------------------------------------------
  if(common_arg->results != 0){
    FILE *fp;
    if( (fp = fopen((char *)common_arg->results, "a")) == NULL ) {
      ERR.FileA(cname,fname, (char *)common_arg->results);
    }
    for(t=0; t<time_size; t++){
      fprintf(fp, "WALLPSEUDO %d  %e %e\n", t,
	    trace[t].real(), trace[t].imag());
    }
    fclose(fp);
  }
  sfree(wprop2);
  sfree(wprop);
  sfree(trace);

}

CPS_END_NAMESPACE
