//------------------------------------------------------------------
//
// alg_fourier_prop.C
//
// AltFourierProp is derived from Alg 
// here we fourier transform point source quark propagators
// to momentum space. very very slowly.
//
//------------------------------------------------------------------

#include <config.h>     
#include <stdlib.h>     // exit()
#include <stdio.h>
#include <math.h>
#include <time.h>

#include <alg/common_arg.h>
#include <alg/alg_fourier_prop.h>
#include <alg/fourierprop_arg.h>
#include <util/lattice.h>
#include <util/data_types.h>
#include <util/gjp.h>
#include <util/smalloc.h>
#include <util/vector.h>
#include <util/verbose.h>
#include <util/error.h>
#include <alg/do_arg.h>
#include <comms/glb.h>
#include <comms/sysfunc_cps.h>
#include <alg/qpropw.h>
#include <alg/qpropw_arg.h>
#include <alg/wilson_matrix.h>
#include <alg/spin_matrix.h>
#include <util/site.h>
#include <util/qcdio.h>

CPS_START_NAMESPACE

void AlgFourierProp::print_header( FILE* fp , const char* name )
{
  Fprintf(fp,"*************************************\n");
  Fprintf(fp,"* %30s  \n",name);
  Fprintf(fp,"*************************************\n");
  Fprintf(fp,"* Number of Momenta: %d \n",(alg_FourierProp_arg->plist.size()));
  Fprintf(fp,"* Mass             : %e \n",alg_FourierProp_arg->cg.mass);
  Fprintf(fp,"* Volume (x->t)    : %d %d %d %d \n",
          GJP.XnodeSites()*GJP.Xnodes(),
          GJP.YnodeSites()*GJP.Ynodes(),
          GJP.ZnodeSites()*GJP.Znodes(),
	  GJP.TnodeSites()*GJP.Tnodes() );
  Fprintf(fp,"* Boundary (x->t)  : %d %d %d %d \n",bc[0],bc[1],bc[2],bc[3]);
  Fprintf(fp,"*************************************\n");
}


AlgFourierProp::AlgFourierProp(Lattice& latt,
			       CommonArg *c_arg,
			       FourierPropArg *arg) :
  Alg(latt, c_arg)
{
  cname = "AlgFourierProp";
  char *fname = "AlgFourierProp(L&,CommonArg*,FourierPropArg*)";
  VRB.Func(cname,fname);

  if(arg == 0)
    ERR.Pointer(cname,fname, "arg");
  alg_FourierProp_arg = arg;

  if ( GJP.Xbc() == BND_CND_PRD ) { bc[0] = 0; } else { bc[0]=1; }
  if ( GJP.Ybc() == BND_CND_PRD ) { bc[1] = 0; } else { bc[1]=1; }
  if ( GJP.Zbc() == BND_CND_PRD ) { bc[2] = 0; } else { bc[2]=1; }
  if ( GJP.Tbc() == BND_CND_PRD ) { bc[3] = 0; } else { bc[3]=1; }
}

AlgFourierProp::~AlgFourierProp() {
  char *fname = "~AlgFourierProp()";
  VRB.Func(cname,fname);
}

void AlgFourierProp::calcmom(const FourMom& p)
{
  char *fname="calcmom(const FourMom&)";
  VRB.Func(cname,fname);

  Float   p1( p.x() + 0.5*bc[0] );
  Float   p2( p.y() + 0.5*bc[1] );
  Float   p3( p.z() + 0.5*bc[2] );
  Float   p4( p.t() + 0.5*bc[3] );
 
  const Float PI(3.141592654);
 
  p1 *= 2.0*PI/(GJP.XnodeSites()*GJP.Xnodes());
  p2 *= 2.0*PI/(GJP.YnodeSites()*GJP.Ynodes());
  p3 *= 2.0*PI/(GJP.ZnodeSites()*GJP.Znodes());
  p4 *= 2.0*PI/(GJP.TnodeSites()*GJP.Tnodes());

  Site site;

  while ( site.LoopsOverNode() )
    {
      const Float px( site.physX()*p1 );
      const Float py( site.physY()*p2 );
      const Float pz( site.physZ()*p3 );
      const Float pt( site.physT()*p4 );
      
      const Float pdotx( px + py + pz + pt );


      Rcomplex fact( cos(pdotx), -sin(pdotx) );

      //--------------------------------------------
      //  Using AddMult is logically equivalent to
      //     momprop += fact * (*prop)[i]; ,
      //  but cuts down on the number of loops over
      //  the data. The time taken is reduced by
      //  a factor of two.
      //--------------------------------------------
      
      momprop.AddMult(fact,(*prop)[site.Index()]);
      
    }
}

//------------------------------------------------------------------
// Fourier Transform for all momenta
//------------------------------------------------------------------

void AlgFourierProp::ft(FILE *fp,MomentaList& mlist)
{
  char *fname="ft(FILE*,MomentaList&)";
  VRB.Func(cname,fname);
  
  //----------------------------------------------------------------
  // begin loop over momenta
  //----------------------------------------------------------------
  int nmom;
  for( nmom=0; nmom < mlist.size(); nmom++ ) 
    {
      //------------------------------------------------------------
      // Calculate single momenta
      //------------------------------------------------------------
      
      momprop=Float(0.0);
      
      calcmom(mlist[nmom]);
      
      //------------------------------------------------------------
      // Global sum
      //------------------------------------------------------------
      
      // 288 = 4 x 3(src) x 4 x 3(snk) x 2(re/im)
      // the final argument means which direction to exclude.
      // here 99 means NO exclusion and end up to be a global sum

      slice_sum((Float*)&momprop, 288, 99);
      
      //------------------------------------------------------------
      // Output to file
      //------------------------------------------------------------
      
      int s1,c1,s2,c2;
      Fprintf(fp,"MOMENTUM= %d %d %d %d \n",
              mlist[nmom].x(),
              mlist[nmom].y(),
              mlist[nmom].z(),
              mlist[nmom].t() );
      for( s1=0; s1<4; ++s1){ 
	for( c1=0; c1<3; ++c1){
	  for( s2=0; s2<4; ++s2){ 
	    for( c2=0; c2<3; ++c2){
	      Fprintf(fp, "%-25.15e %-25.15e \n",
		      (Float)momprop.wmat().d[s1].c[c1].d[s2].c[c2].real(),
		      (Float)momprop.wmat().d[s1].c[c1].d[s2].c[c2].imag());
	    }
	  }
	}
      }
    }   // end loop over momenta
}


void AlgFourierProp::run()
{
  char *fname = "run()";
  VRB.Func(cname,fname);

  //----------------------------------------------------------------
  // Set the Lattice pointer
  //----------------------------------------------------------------

  Lattice& lat( AlgLattice() );

  //----------------------------------------------------------------
  // set the variables
  //----------------------------------------------------------------

  QPropWArg qarg;
  qarg.cg = alg_FourierProp_arg->cg;
  qarg.x  = alg_FourierProp_arg->x_src;
  qarg.y  = alg_FourierProp_arg->y_src;
  qarg.z  = alg_FourierProp_arg->z_src;
  qarg.t  = alg_FourierProp_arg->t_src;
  qarg.gauge_fix_snk  = 1;
  
  //----------------------------------------------------------------
  // Fetch propagator
  //----------------------------------------------------------------
  
  QPropWGFPointSrc prop0( lat, &qarg, common_arg );
  
  prop0.Run();

  //----------------------------------------------------------------
  // set current propagator
  //----------------------------------------------------------------

  prop=&prop0;
  
  //----------------------------------------------------------------
  // output header. 
  //----------------------------------------------------------------
  
  FILE* fp(NULL);
  
  if( (fp = Fopen((char *)alg_FourierProp_arg->results, "a")) == NULL ) {
    ERR.FileA(cname,fname, (char *)alg_FourierProp_arg->results);
  }
  print_header(fp,"AlgFourierProp");
  
  ft(fp,alg_FourierProp_arg->plist);
  
  Fclose(fp);
}


//-------------------------------------------------------------------
// 
//  AlgFourierPropDis:
//
//
//-------------------------------------------------------------------
void AlgFourierPropDis::run()
{
  char *fname = "run()";
  VRB.Func(cname,fname);

  //----------------------------------------------------------------
  // Set the Lattice pointer
  //----------------------------------------------------------------

  Lattice& lat = AlgLattice();

  //----------------------------------------------------------------
  // get the props
  //----------------------------------------------------------------

  //----------------------------------------------------------------
  // set the point source
  //----------------------------------------------------------------

  QPropWArg qarg;
  qarg.cg = alg_FourierProp_arg->cg;
  qarg.x  = alg_FourierProp_arg->x_src;
  qarg.y  = alg_FourierProp_arg->y_src;
  qarg.z  = alg_FourierProp_arg->z_src;
  qarg.t  = alg_FourierProp_arg->t_src;
  qarg.gauge_fix_snk  = 1;
  
  FILE *fp,*fp2;
  if( (fp = Fopen((char *)alg_FourierProp_arg->results, "a")) == NULL ) 
    {
      ERR.FileA(cname,fname, (char *)alg_FourierProp_arg->results);
    }
  if( (fp2 = Fopen((char *)alg_FourierProp_arg->results2, "a")) == NULL ) 
    {
      ERR.FileA(cname,fname, (char *)alg_FourierProp_arg->results2);
    }
  
  //------------------------------------------------------------------
  // Calculate momentum space propagator
  //------------------------------------------------------------------
  //
  // Enter new scope so memory for propagator is freed before going
  // on to invert again for the disconnected part.
  //
  { 
    //----------------------------------------------------------------
    // Fetch propagator
    //----------------------------------------------------------------
    
    printf("Invert Dirac Op.\n");

    QPropWGFPointSrc prop0( lat, &qarg, common_arg);
    
    prop0.Run();

    printf("Inversion complete!\n");
    
    //----------------------------------------------------------------
    // set current propagator
    //----------------------------------------------------------------
    
    prop=&prop0;
    
    //----------------------------------------------------------------
    // output headers. 
    //----------------------------------------------------------------
    
    
    print_header(fp,"AlgFourierPropDis");
    
    Fprintf(fp2,"*************************************\n");
    Fprintf(fp2,"*    AlgFourierPropDis  Header      *\n");
    Fprintf(fp2,"*************************************\n");
    Fprintf(fp2,"Mass             : %e \n",alg_FourierProp_arg->cg.mass);
    Fprintf(fp2,"Volume (x->t)    : %d %d %d %d \n"
            ,GJP.XnodeSites()*GJP.Xnodes()
            ,GJP.YnodeSites()*GJP.Ynodes()
            ,GJP.ZnodeSites()*GJP.Znodes()
            ,GJP.TnodeSites()*GJP.Tnodes());
    Fprintf(fp2,"Boundary (x->t)  : %d %d %d %d \n",bc[0],bc[1],bc[2],bc[3]);
    Fprintf(fp2,"Number of Fixed  : %d \n",alg_FourierProp_arg->smom.size());
    Fprintf(fp2,"*************************************\n");
    
    //----------------------------------------------
    // output quark loop contribution to fp2
    // This is just the propagator at the origin 
    // (i=0) on node 0.
    //----------------------------------------------
    
    Fprintf(fp2,"*************************************\n");
    Fprintf(fp2,"* quark loop                        *\n");
    Fprintf(fp2,"*************************************\n");
    
    for(int s1=0;s1<4;++s1){
      for(int c1=0;c1<3;++c1){
        for(int s2=0;s2<4;++s2){
          for(int c2=0;c2<3;++c2){
            Fprintf(fp2, "%-25.15e %-25.15e \n",
                    (Float)prop0[0].wmat().d[s1].c[c1].d[s2].c[c2].real(),
                    (Float)prop0[0].wmat().d[s1].c[c1].d[s2].c[c2].imag());
          }
        }
      }
    }
    
    //----------------------------------------------------------------
    // fourier tranform with NO input momenta (point-point prop)
    // (non-fixed momenta list is DIFFERENT for each FIXED momenta)
    //----------------------------------------------------------------
        
    for (int s=0;s<alg_FourierProp_arg->smom.size();++s)
      {
        Fprintf(fp2,"*************************************\n");
        Fprintf(fp2,"* point-point                       *\n");
        Fprintf(fp2,"*************************************\n");
        Fprintf(fp2,"Fixed Mom             : %d %d %d %d \n"
                ,alg_FourierProp_arg->smom[s].x()
                ,alg_FourierProp_arg->smom[s].y()
                ,alg_FourierProp_arg->smom[s].z()
                ,alg_FourierProp_arg->smom[s].t());
        Fprintf(fp2,"Number of Non-Fixed   : %d \n",
                alg_FourierProp_arg->smom_comp[s].size());
        Fprintf(fp2,"*************************************\n");
        
        
        ft(fp2,alg_FourierProp_arg->smom_comp[s]);
        
      }
    


    //-----------------------------------------------
    // fourier for stand-alone momenta list
    //-----------------------------------------------
    
    ft(fp,alg_FourierProp_arg->plist);

  } //exit trivial scope for point-point correlator, memory should be released now
  
  Fclose(fp);



  
  //------------------------------------------------------------------
  // calculate the 'disconnected' contribution
  // (Prop with input momentum at source)
  //------------------------------------------------------------------
  
  if ( discon )
    {

      Fprintf(fp2,"*************************************\n");
      Fprintf(fp2,"*    AlgFourierPropDis  Header      *\n");
      Fprintf(fp2,"*************************************\n");
      Fprintf(fp2,"Mass             : %e \n",alg_FourierProp_arg->cg.mass);
      Fprintf(fp2,"Number of Fixed  : %d \n",alg_FourierProp_arg->smom.size());
      Fprintf(fp2,"*************************************\n");
      
      
      //----------------------------
      // loop over input momenta
      //----------------------------
      
      
      for (int s=0;s<alg_FourierProp_arg->smom.size();++s)
        {
          //----------------------------------------------------------------
          // output header. 
          //----------------------------------------------------------------
          
          Fprintf(fp2,"*************************************\n");
          Fprintf(fp2,"* Disconnected                      *\n");
          Fprintf(fp2,"*************************************\n");
          Fprintf(fp2,"Fixed Mom             : %d %d %d %d \n"
                  ,alg_FourierProp_arg->smom[s].x()
                  ,alg_FourierProp_arg->smom[s].y()
                  ,alg_FourierProp_arg->smom[s].z()
                  ,alg_FourierProp_arg->smom[s].t());
          Fprintf(fp2,"Number of Non-Fixed   : %d \n",
                  alg_FourierProp_arg->smom_comp[s].size());
          Fprintf(fp2,"*************************************\n");
          
          //----------------------------------------------------------------
          // make propagator with fixed input momenta
          //----------------------------------------------------------------
          
          // this call does not run the inverter
          
          QPropWLandauGaugeVolumeSrc prop0( lat, 
                                            &qarg,
                                            common_arg );
          // set fixed momenta
          prop0.SetMomenta(alg_FourierProp_arg->smom[s].as_array());
          
          // invert propagator
          prop0.Run();
          
          prop=&prop0;
          
          //-------------------------------------------------------------------
          // fourier tranform of the non-fixed momenta with input momenta set
          //-------------------------------------------------------------------
                    
          ft(fp2,alg_FourierProp_arg->smom_comp[s]);

        } // s

    } // if (discon)

  Fclose(fp2);
  
  if(common_arg->results != 0)
    {
      FILE *fp;
      if( (fp = Fopen((char *)common_arg->results, "a")) == NULL ) {
        ERR.FileA(cname,fname, (char *)common_arg->results);
      }
      Fprintf(fp, "fourier prop:finished mass %e \n",alg_FourierProp_arg->cg.mass );
      Fclose(fp);
    }
}

CPS_END_NAMESPACE
