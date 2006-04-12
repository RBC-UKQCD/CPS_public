#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2006-04-12 22:08:21 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_w_spect/alg_w_spect.C,v 1.14 2006-04-12 22:08:21 chulwoo Exp $
//  $Id: alg_w_spect.C,v 1.14 2006-04-12 22:08:21 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: alg_w_spect.C,v $
//  $Revision: 1.14 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_w_spect/alg_w_spect.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------

//--------------------------------------------------------------------------
// alg_w_spect.C
//
// AlgWspect is derived from Alg and is relevant to  
// spectroscopy with Wilson type fermions (i.e. Wilson, Clover, Dwf). 
// The type of fermion is determined by the argument to the 
// constructor but it should not be a non Wilson type lattice.
//--------------------------------------------------------------------------


CPS_END_NAMESPACE
#include <alg/alg_w_spect.h>         // class AlgWspect
#include <alg/common_arg.h>          // class CommonArg
#include <alg/w_spect_arg.h>         // class WspectArg
#include <util/lattice.h>            // class Lattice
#include <util/verbose.h>            // VRB
#include <util/error.h>              // ERR
#include <alg/w_all.h>                                 // class Wspect*
CPS_START_NAMESPACE

#define DEBUG_ALG_W_SPECT


//---------------------------------------------------------------------------
// For the purpose of debugging or timing during code upgrade
//---------------------------------------------------------------------------

//#define TIMING_ALG_W_SPECT

#ifdef  TIMING_ALG_W_SPECT
CPS_END_NAMESPACE
  #include <time.h>                 // clock()
CPS_START_NAMESPACE
#endif


//--------------------------------------------------------------------------
// Static data members
//--------------------------------------------------------------------------
char * AlgWspect::d_class_name = "AlgWspect";
int    AlgWspect::d_counter    = 1;
int    AlgWspect::d_count_step = 1;


//--------------------------------------------------------------------------
// Static member function
//--------------------------------------------------------------------------
void
AlgWspect::SetCounter(int counter, int step)  
{ 
  d_counter = counter;  
  d_count_step = step; 
}
  
//--------------------------------------------------------------------------
// AlgWspect CTOR
//--------------------------------------------------------------------------
AlgWspect::AlgWspect(Lattice& latt, 
		     CommonArg *c_arg,
		     WspectArg *w_arg,
		     CgArg *cg,
		     int n_quark_masses) 
  : Alg(latt, c_arg),
    d_arg_p(w_arg),
    cg_arg_p(cg),
    d_num_args(n_quark_masses)
{
//  cg_arg = cg;

  // Obtain an instance of the support class WspectGinfo
  //------------------------------------------------------------------------
  WspectGinfo g_info;
  VRB.Func(d_class_name, g_info.ctor_str);


  // Check the lattice type
  //------------------------------------------------------------------------
  if ((latt.Fclass() != F_CLASS_WILSON) &&
      (latt.Fclass() != F_CLASS_CLOVER) &&
      (latt.Fclass() != F_CLASS_DWF))        {
    ERR.General(d_class_name,
		g_info.ctor_str,
		g_info.wrong_type_str);
  }

  // Check the WspectArg pointer
  //------------------------------------------------------------------------
  if (d_arg_p == 0) {
    ERR.Pointer(d_class_name,
		g_info.ctor_str, 
		"arg");
  }
  if (cg_arg_p == 0) {
    ERR.Pointer(d_class_name,
                g_info.ctor_str,
                "arg");
  }


  // Check the num of masses, prop_dir and AOTS
  //------------------------------------------------------------------------
  int prop_dir   = d_arg_p->prop_dir;
  int aots_start = d_arg_p->aots_start;
  int aots_step  = d_arg_p->aots_step;
  int aots_num   = d_arg_p->aots_num;
  if (d_num_args < 1 ||
      prop_dir < 0   || prop_dir >= g_info.LORENTZs ||
      aots_start+aots_step*(aots_num-1) >= g_info.glb_sites[prop_dir]) {
    ERR.General(d_class_name, 
		g_info.ctor_str, 
		g_info.out_range_str);
  }

  // Check the prop consistency among different quark masses
  //------------------------------------------------------------------------
  // Ping:  will we use different aots for non-degenerate quarks
  //        and then mix them up to form a hadron???
  for (int i = 0; i < d_num_args; ++i) {
    if (d_arg_p[i].prop_dir    != prop_dir   ||
	d_arg_p[i].aots_start  != aots_start ||
	d_arg_p[i].aots_step   != aots_step  ||
	d_arg_p[i].aots_num    != aots_num) {
      ERR.General(d_class_name, 
		  g_info.ctor_str, 
		  g_info.inconsistent_str);
    }
  }
}


//--------------------------------------------------------------------------
// Destructor
//--------------------------------------------------------------------------
AlgWspect::~AlgWspect() {
}


//--------------------------------------------------------------------------
//
//--------------------------------------------------------------------------
void AlgWspect::run()
{
#ifdef  TIMING_ALG_W_SPECT
  int quark_b, quark_e;
  int meson_b, meson_m, meson_e, meson_bmp, meson_mmp, meson_emp;
  int nucleon_b, nucleon_m, nucleon_e;
  int total_b = clock();
  int total_e;  
#endif
	  CgArg cg = *cg_arg_p;
	  char *fname = "run()";
	  VRB.Func(d_class_name,fname);

	  // printf("in AlgWspect::run \n");

	  WspectOutput * output = (WspectOutput *)common_arg->results;
	  

	  // Set the Lattice pointer
	  //------------------------------------------------------------------------
	  Lattice& lat = AlgLattice();

	  int src_slice = d_arg_p->aots_start;
	  int src_slice_step = d_arg_p->aots_step;
	  int src_slice_end  = src_slice + src_slice_step * d_arg_p->aots_num;

	  VRB.Result(d_class_name,fname,"%d %d %d \n",src_slice,src_slice_step, src_slice_end);

	  for ( ; src_slice < src_slice_end; src_slice += src_slice_step) {


	    // Calculate quark propagator
	    //----------------------------------------------------------------------
	    // Ping:  certainly more work here to be done about the desired
	    //        combinations of non-degenerate quarks.
	    //        Presumably, more arguments will have to be passed in.
	    //        One way: for three flavors, [100] means use only q1
	    //                 to caculate spectrum.
	    //        Also some care needed to get the scope (CTOR and DTOR) 
	    //        of each quark propagator right.
	    //    const WspectQuark & q2 = q;
	    //    const WspectQuark & q3 = q;

	    // Xiaodong & Thomas:
	    // Modified to calculate also extended mesons
	    // q1 is the usual propagator(no source operator), which can be
	    // used to pass propagation direction and src_slice infomation
	    // to spectrum class

	    // there is a problem here --> check !
	    VRB.Result(d_class_name,fname,"prop_dir = %d , src_slice = %d \n",d_arg_p->prop_dir, src_slice);

	    WspectHyperRectangle hyperRect(d_arg_p->prop_dir, src_slice);    

	#ifdef  TIMING_ALG_W_SPECT
	    quark_b = clock();
	#endif

    VRB.Result(d_class_name,fname,"created quark q1 \n");
	    // create local quark propagator
	     WspectQuark q1(lat, output->cg, output->pbp,
			  output->mid_point, output->a0_p, d_arg_p[0], cg,hyperRect);
		  
      VRB.Result(d_class_name,fname,"finished quark q1 \n");
#ifdef  TIMING_ALG_W_SPECT
    quark_e = clock();
#endif

    //Note: for ExtendedMesons, do only zero momentum projection
    WspectMomenta  mom(hyperRect, q1.SourceCenter2(), d_arg_p->num_mom - 1);
    //    mom.dumpData();


    // Calculate LOCAL meson CORRELATOR
    // added control by Thomas and Xiaodong
    //----------------------------------------------------------------------
    {
#ifdef  TIMING_ALG_W_SPECT
      meson_b = clock();
#endif
     if(d_arg_p->normal_mesons_on) {
        WspectMesons mes(q1, q1, hyperRect, mom);
#if 0
	q1.dumpData("qprop.dat");
#endif

#ifdef  TIMING_ALG_W_SPECT
        meson_m = clock();
#endif
        //write data to files
        mes.print(output);
      }
#ifdef  TIMING_ALG_W_SPECT
      meson_e = clock();
#endif
    } //end of normal mesons
   

    // Calculate <\Delta J^5 \bar q1 \gamma^5 q1> with middle point sink
    // changed
    //----------------------------------------------------------------------
    if (lat.Fclass() == F_CLASS_DWF && output->mid_point)
      {
#ifdef  TIMING_ALG_W_SPECT
	meson_bmp = clock();
#endif
        WspectMesons mes(q1.Data_SP1(), q1.Data_SP2(), hyperRect, mom);
	
#ifdef  TIMING_ALG_W_SPECT
	meson_mmp = clock();
#endif
	
	mes.print_mp(output->mid_point);	
	
#ifdef  TIMING_ALG_W_SPECT
	meson_emp = clock();
#endif
      }
    
    // Calculate nucleon and delta's
    //----------------------------------------------------------------------
    if (d_arg_p->baryons_on) {
      
      {
#ifdef  TIMING_ALG_W_SPECT
	nucleon_b = clock();
#endif
	WspectBaryon nuc(q1, q1, q1, hyperRect,
			 WspectBaryon::NUCLEON_CONSTI,
			 WspectBaryon::NUCLEON_DIRAC);
#ifdef  TIMING_ALG_W_SPECT
	nucleon_m = clock();
#endif

	nuc.print(output->nucleon, output->fold);	

#ifdef  TIMING_ALG_W_SPECT
	nucleon_e = clock();
#endif   
      }

      
      {
	WspectBaryon nucPrime(q1, q1, q1, hyperRect,
			      WspectBaryon::NUCLEON_CONSTI,
			      WspectBaryon::UnitUnit);
	nucPrime.print(output->nucleon_prime, output->fold);	
      }
      
      
      {
	WspectBaryon deltaX(q1, q1, q1, hyperRect,
			    WspectBaryon::DELTA_CONSTI,
			    WspectBaryon::DELTAX_DIRAC);
	deltaX.print(output->delta_x, output->fold);	
      }
      
      
      {
	WspectBaryon deltaY(q1, q1, q1, hyperRect,
			    WspectBaryon::DELTA_CONSTI,
			    WspectBaryon::DELTAY_DIRAC);
	deltaY.print(output->delta_y, output->fold);      
      }
      
      
      {
	WspectBaryon deltaZ(q1, q1, q1, hyperRect,
			    WspectBaryon::DELTA_CONSTI,
			    WspectBaryon::DELTAZ_DIRAC);
	deltaZ.print(output->delta_z, output->fold); 
      }
      

      {
	WspectBaryon deltaT(q1, q1, q1, hyperRect,
			    WspectBaryon::DELTA_CONSTI,
			    WspectBaryon::DELTAT_DIRAC);
	deltaT.print(output->delta_t, output->fold); 
      }
    } //end if(baryons_on) 
    
    // Increment the counter
    d_counter += d_count_step;
  } // end of for(sc_slice,..)

#ifdef  TIMING_ALG_W_SPECT
  total_e = clock();
  printf("Total: %d = [%d - %d]\n", total_e - total_b, total_e, total_b);
  printf("Quark: %d = [%d - %d]\n", quark_e - quark_b, quark_e, quark_b);
  printf("Meson: \t%d = [%d - %d] = \n\tcalc %d = [%d - %d] + \n\tprint %d = [%d - %d]\n", 
	 meson_e - meson_b, meson_e, meson_b,
	 meson_m - meson_b, meson_m, meson_b,
	 meson_e - meson_m, meson_e, meson_m);
  printf("Nucleon:\t%d = [%d - %d] = \n\tcalc %d = [%d - %d] + \n\tprint %d = [%d - %d]\n", 
	 nucleon_e - nucleon_b, nucleon_e, nucleon_b,
	 nucleon_m - nucleon_b, nucleon_m, nucleon_b,
	 nucleon_e - nucleon_m, nucleon_e, nucleon_m);  
#endif 
   VRB.FuncEnd(d_class_name,fname);

}

CPS_END_NAMESPACE
