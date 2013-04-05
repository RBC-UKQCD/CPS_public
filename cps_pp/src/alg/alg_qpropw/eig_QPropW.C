#include <config.h> 
#include <stdlib.h>     // exit()
#include <stdio.h>
#include <string.h>
#include <alg/common_arg.h>
#include <comms/glb.h>
#include <comms/scu.h>
// #include <util/data_io.h>
#include <comms/sysfunc_cps.h>

#include <fcntl.h>      // read and write control flags,
#include <unistd.h>     // close(). These are needed for io parts to
                        // compile on PCs
#include <alg/qpropw.h>
#include <util/qcdio.h>

#include <util/qioarg.h>
#include <util/sproj_tr.h>
#include <util/time_cps.h>
#include <util/vector_template_func.h>

//YA
#include <alg/alg_plaq.h>
#include <alg/alg_smear.h>
#include <alg/no_arg.h>
#include <util/dirac_op.h>

#define VOLFMT QIO_VOLFMT

CPS_START_NAMESPACE
//eig_int_CG. by Qi Liu

void invert_H_matrix(Rcomplex *data, int n)
{
	int actualsize=n;
    if (actualsize <= 0) return;  // sanity check
    if (actualsize == 1) {data[0]=1.0/data[0];return;}  // must be of dimension >= 2

    for (int i=1; i < actualsize; i++) data[i] /= data[0]; // normalize row 0

    for (int i=1; i < actualsize; i++)  
	{ 
      for (int j=i; j < actualsize; j++)  { // do a column of L
        Rcomplex sum = 0.0;
        for (int k = 0; k < i; k++)  
            sum += data[j*actualsize+k] * data[k*actualsize+i];
        data[j*actualsize+i] -= sum;
        }
      if (i == actualsize-1) continue;
      for (int j=i+1; j < actualsize; j++)  {  // do a row of U
        Rcomplex sum = 0.0;
        for (int k = 0; k < i; k++)
            sum += data[i*actualsize+k]*data[k*actualsize+j];
        data[i*actualsize+j] = 
           (data[i*actualsize+j]-sum) / data[i*actualsize+i];
        }
      }
    for ( int i = 0; i < actualsize; i++ )  // invert L
      for ( int j = i; j < actualsize; j++ )  {
        Rcomplex x = 1.0;
        if ( i != j ) {
          x = 0.0;
          for ( int k = i; k < j; k++ ) 
              x -= data[j*actualsize+k]*data[k*actualsize+i];
          }
        data[j*actualsize+i] = x / data[j*actualsize+j];
        }
    for ( int i = 0; i < actualsize; i++ )   // invert U
      for ( int j = i; j < actualsize; j++ )  {
        if ( i == j ) continue;
        Rcomplex sum = 0.0;
        for ( int k = i; k < j; k++ )
            sum += data[k*actualsize+j]*( (i==k) ? 1.0 : data[i*actualsize+k] );
        data[i*actualsize+j] = -sum;
        }
    for ( int i = 0; i < actualsize; i++ )   // final inversion
      for ( int j = 0; j < actualsize; j++ )  {
        Rcomplex sum = 0.0;
        for ( int k = ((i>j)?i:j); k < actualsize; k++ )  
            sum += ((j==k)?1.0:data[j*actualsize+k])*data[k*actualsize+i];
        data[j*actualsize+i] = sum;
        }
}

void QPropW::eig_Run(Vector **V, const int vec_len, Float *M, Float max_eig, const int nev, const int m, float **U, Rcomplex *H, const int max_def_len, int &def_len, const Float *restart, const int restart_len, const bool always_restart, const int do_rerun, const Float precision) 
{

   char *fname = "eig_Run()";
   VRB.Func(cname, fname);

   // Set the node size of the full (non-checkerboarded) fermion field
   //----------------------------------------------------------------
   //int f_size = GJP.VolNodeSites() * Lat.FsiteSize()/GJP.SnodeSites();
   int iter;
   Float true_res;
   if(def_len>max_def_len)
   {
	   VRB.Warn(cname,fname,"def_len>max_def_len, reset to equal");
	   def_len=max_def_len;
   }

   int Nspins = 4; // Number of spin components to be done
   
   // Flag set if sequential propagator 
   int seq_src = ((SrcType()==PROT_U_SEQ)||
				  (SrcType()==PROT_D_SEQ)||
				  (SrcType()==MESSEQ)      );

   if (DoHalfFermion()) Nspins = 2;

   // does prop exist? Assume it does not.
   int do_cg = 1;

   /****************************************************************
     The code below is temporarily isolated for purposes of merging
     with CPS main branch,  12/09/04, Oleg Loktik
   -------------------- Quarantine starts --------------------------

   if (pfs_file_exists(Arg.file)){
     // read data into prop
     RestoreQProp(Arg.file,PROP); // Only restores stuff into prop 

     // multiply source by 1/2(1+gamma_t). If the propagator  
     // on disk is half fermion it does nothing. Otherwise it gets
     // converted to the half fermion propagator. 
     if (Arg.DoHalfFermion) 
       for (int s(0);s<GJP.VolNodeSites();s++)
	 prop[s].PParProjectSink() ;
      
     do_cg = 0;
   }
   -------------------- End of quarantine -------------------------*/ 

  //-----------------------------------------------------------------
  // TY Add Start
   // we need to store the source
   Float *save_source=NULL;
  // TY Add End
  //-----------------------------------------------------------------
   WilsonMatrix *read_prop=NULL;
   WilsonMatrix *save_prop=NULL;

   if (do_cg) {

     Allocate(PROP); 
     if (StoreMidprop()) Allocate(MIDPROP);
     
     FermionVectorTp src;
     FermionVectorTp sol;
     FermionVectorTp midsol;
     
  //-----------------------------------------------------------------
  // TY Add Start
  // For conserved axial current
     int glb_walls = GJP.TnodeSites()*GJP.Tnodes();
     int fsize = glb_walls * sizeof(Float);
     conserved = (Float *) smalloc(fsize);
     if (!conserved)
       ERR.Pointer(cname, fname, "d_conserved_p");
     VRB.Smalloc(cname, fname, "d_conserved_p", conserved, fsize);
  
     Float* flt_p = (Float *)conserved;
     for ( int i = 0; i < glb_walls; i++) *flt_p++ = 0.0;

     spnclr_cnt = 0;

     // we need to store the source
     if (qp_arg.save_prop || do_rerun ){ 
       save_source = (Float*)smalloc(GJP.VolNodeSites()*288*sizeof(Float));
       if (save_source == 0) ERR.Pointer(cname, fname, "save_source");
       VRB.Smalloc(cname, fname, "save_source", save_source,
		   GJP.VolNodeSites() * 288*sizeof(Float));
     }
     // TY Add End
     //-----------------------------------------------------------------

     // in case we do a rerun, we also need to store a propagator
     if(do_rerun){
       read_prop = (WilsonMatrix*)smalloc(GJP.VolNodeSites()*sizeof(WilsonMatrix));
       if (read_prop == 0) ERR.Pointer(cname, fname, "pr3op");
       VRB.Smalloc(cname, fname, "read_prop", read_prop,
		   GJP.VolNodeSites() * sizeof(WilsonMatrix));     
       
#ifdef USE_QIO
       qio_readPropagator readPropQio(qp_arg.file, QIO_FULL_SOURCE, read_prop, save_source,
				      GJP.argc(), GJP.argv(), VOLFMT);
#endif //USE_QIO

       if(AlgLattice().Fclass() == F_CLASS_DWF){
	 
	 Float renFac = 5. - GJP.DwfHeight();
	 
	 for(int ii(0); ii <  GJP.VolNodeSites(); ++ii)
	   *(read_prop +ii) *= renFac;
       }
       
     }
     
     //-----------------------------------------------------------------
     // M. Lightman
     // For m_res
     j5q_pion = (Float *) smalloc(fsize);
     if (!j5q_pion)
       ERR.Pointer(cname, fname, "d_j5q_pion_p");
     VRB.Smalloc(cname, fname, "d_j5q_pion_p", j5q_pion, fsize);
     
     flt_p = (Float *) j5q_pion;
     for ( int i = 0; i < glb_walls; i++) *flt_p++ = 0.0;
     // End M. Lightman
     //-----------------------------------------------------------------
     

     int *low=NULL;
     if(nev>0)low=(int *)smalloc(2*nev*sizeof(int));
     Rcomplex *invH=(Rcomplex *)smalloc(max_def_len*max_def_len*sizeof(Rcomplex));
     Vector *AU=(Vector *)fmalloc(cname,fname,"AU",vec_len * sizeof(Float));
     Lattice& Lat = this->AlgLattice() ;
     
     for (int spn=0; spn < Nspins; spn++)
       for (int col=0; col < GJP.Colors(); col++) {
	 
	 // initial guess (Zero)
	 sol.ZeroSource();
	 
	 if(!do_rerun){
	   SetSource(src,spn,col);
	   
	   // store the source
	   if (qp_arg.save_prop) {
	     for(int index(0); index < GJP.VolNodeSites(); ++index)
	       for(int mm(0); mm < 4; ++ mm)
		 for(int cc(0); cc < GJP.Colors(); ++cc){
		   // now same ordering as propagator [volume][spin][color][solution_spin][solution_color][ReIm]
		   *(save_source + 288*index + 72*mm + 24*cc + 6*spn + 2*col)       = src[24*index + 6*mm + 2*cc];
		   *(save_source + 288*index + 72*mm + 24*cc + 6*spn + 2*col + 1)   = src[24*index + 6*mm + 2*cc+1];
		 }
	   }
	 }
	 else{ // rerun
	   for(int index(0); index < GJP.VolNodeSites(); ++index){
	     
	     WilsonMatrix *tmp_mat= (WilsonMatrix *) save_source+index;
	     
	     src.CopyWilsonMatSink(index, spn, col,*tmp_mat);
	   }
	 }
	 
	 if ((DoHalfFermion())&&(!seq_src)) // Rotate to chiral basis
	   src.DiracToChiral();

	 // Get the prop
	 VRB.Debug(cname,fname,"Before CG in QpropW.Run() \n");

	 //calculate invH from H (def_len*def_len matrix)
	 if(def_len>0)
	   {	    
	     //initialize H if necessary
	     if(nev == 0){
	       Vector *temp=(Vector *)fmalloc(cname,fname,"temp",vec_len * sizeof(Float));
	       DiracOpDwf dwf_aux(Lat, NULL, NULL, &(qp_arg.cg),CNV_FRM_NO);
	       for(int i=0;i<max_def_len;i++){
		 for(int j=0;j<vec_len;++j) *((Float*)temp + j) = U[i][j];
		 dwf_aux.MatPcDagMatPc(AU,temp);
		 for(int j=0;j<i;j++)
		   {
		     Float c_r, c_i;
		     compDotProduct<float,Float>(&c_r, &c_i, U[j], (Float *)AU,vec_len);
		     glb_sum_five(&c_r);
		     glb_sum_five(&c_i);
		     H[j*max_def_len+i]=Complex(c_r,c_i);
		   }
		 for(int j=0;j<i;j++)H[i*max_def_len+j]=conj(H[j*max_def_len+i]);
		 H[i*max_def_len+i]=temp->ReDotProductGlbSum(AU,vec_len);
	       }
	       sfree(temp);
	     }
	     // Hinv
	     for(int i=0;i<def_len;i++)
	       for(int j=0;j<def_len;j++)
		 invH[i*def_len+j]=H[i*max_def_len+j];
	     invert_H_matrix(invH, def_len);
	   }
	 if(!always_restart && nev>0 && def_len<max_def_len)
	   {
	     eig_CG(V, vec_len, M, nev, m, U, invH, def_len, restart, 0, src, sol, midsol, iter, true_res);
	   }
	 else{
	   if(def_len==max_def_len)eig_CG(V, vec_len, M, 0, 0, U, invH, def_len, restart,restart_len, src, sol, midsol, iter, true_res);
	   else eig_CG(V, vec_len, M, nev, m, U, invH, def_len, restart,restart_len, src, sol, midsol, iter, true_res);
	 }
	 
	 if(def_len<max_def_len)
	   {
	     Float *fvptr=NULL;
	     DiracOpDwf dwf_aux(Lat, NULL, NULL, &(qp_arg.cg),CNV_FRM_NO);
	     //add low modes from the lowest first
	     int s;
	     for(s=0;s<2*nev;s++)low[s]=s;
	     for(int i=0;i<2*nev-1;i++)
	       {
		 for(int j=2*nev-1;j>i;j--)
		   {
		     if(M[low[j]]<M[low[j-1]]){s=low[j];low[j]=low[j-1];low[j-1]=s;}
		   }
	       }
	     for(int i=0;i<2*nev;i++)
	       {
		 VRB.Result(cname,fname,"eigen value %d is %e \n",i,M[low[i]]);
		 //update deflation space U from V
		 if(M[low[i]]<max_eig && M[low[i]]>1e-30 && def_len<max_def_len) //remember to set M[i>=rank] to zero
		   {
		     //update H 
		     //U[def_len]->CopyVec(V[low[i]],vec_len);
		     fvptr = (Float *)V[low[i]];
		     for(int ii=0;ii<vec_len;ii++){U[def_len][ii]=(float)(fvptr[ii]);fvptr[ii]=(Float)(U[def_len][ii]);}//Make low accuracy
		     dwf_aux.MatPcDagMatPc(AU,V[low[i]]);
		     //for(int j=0;j<def_len;j++)H[j*max_def_len+def_len]=U[j]->CompDotProductGlbSum(AU,vec_len);
		     for(int j=0;j<def_len;j++)
		       {
			 Float c_r, c_i;
			 compDotProduct<float,Float>(&c_r, &c_i, U[j], (Float *)AU,vec_len);
			 glb_sum_five(&c_r);
			 glb_sum_five(&c_i);
			 H[j*max_def_len+def_len]=Complex(c_r,c_i);
		       }
		     for(int j=0;j<def_len;j++)H[def_len*max_def_len+j]=conj(H[j*max_def_len+def_len]);
		     H[def_len*max_def_len+def_len]=V[low[i]]->ReDotProductGlbSum(AU,vec_len);
		     
		     def_len++;	
		   }
	       }
	   }
	 //gauge fix solution
	 FixSol(sol);
	 if (StoreMidprop()) FixSol(midsol);
	 
	 // Collect solutions in propagator.
	 LoadRow(spn,col,sol,midsol);
	 
	 if (DoHalfFermion()) {// copy spin 0 to spin 1 and spin 2 to spin 3
	   int spn2 = spn + 2;
	   if (seq_src) {
	     LoadRow(spn2,col,sol,midsol);
	   } else { // Regular propagator zero the extra components
	     src.ZeroSource();
	     LoadRow(spn2,col,src,src);
	   }
	 }
		 
	 if (common_arg->results != 0) {
	   FILE *fp;
	   if ((fp = Fopen((char *)common_arg->results, "a")) == NULL) {
	     ERR.FileA(cname,fname, (char *)common_arg->results);
	   }
	   Fprintf(fp, "Cg iters = %d true residual = %e\n",
		   iter, (Float)true_res);
	   Fclose(fp);
	 }
	 
       } // End spin-color loop
	sfree(invH);
	sfree(AU);
	if(nev>0)sfree(low);
     
     // Rotate the source indices to Chiral basis if needed
     if ((DoHalfFermion())&&(!seq_src)) {	
       for (int s=0;s<GJP.VolNodeSites();s++)
	 prop[s].SinkChiralToDirac(); // multiply by V^\dagger
       
       if (StoreMidprop())
	 for (int s=0;s<GJP.VolNodeSites();s++)
	   midprop[s].SinkChiralToDirac(); // multiply by V^\dagger
     }
   }

   //-----------------------------------------------------------------
   // TY Add Start
   // Print out conserved axial results
   int time_size = GJP.TnodeSites()*GJP.Tnodes();
   for(int t(0);t<time_size;t++)
     slice_sum((Float*)&conserved[t], 1, 99);
   if(common_arg->results != 0){
     FILE *fp;
     if( (fp = Fopen((char *)common_arg->results, "a")) == NULL ) {
       ERR.FileA(cname,fname, (char *)common_arg->results);
     }
     Fprintf(fp,"Conserved Axial w_spect\n");
     for(int t=0; t<time_size; t++){
       Fprintf(fp,"%d = %.16e\n", t,conserved[t]);
     }
     Fclose(fp);
   }
   sfree(conserved);
   // TY Add End
   //-----------------------------------------------------------------

   //-----------------------------------------------------------------
   // M. Lightman
   // Print out J5q Pion contraction
   for(int t(0);t<time_size;t++)
     slice_sum((Float*)&j5q_pion[t], 1, 99);
   if(common_arg->results != 0){
     FILE *fp1;
     if( (fp1 = Fopen((char *)common_arg->results, "a")) == NULL ) {
       ERR.FileA(cname,fname, (char *)common_arg->results);
     }
     Fprintf(fp1,"J5q Pion Contraction\n");
     for(int t=0; t<time_size; t++){
       Fprintf(fp1,"%d = %.16e\n", t, j5q_pion[t]);
     }
    Fclose(fp1);
   }
   sfree(j5q_pion);
   // End M. Lightman
   //-----------------------------------------------------------------

   // save prop
   if (do_cg && qp_arg.save_prop) {
     
     char propType[256], sourceType[256], propOutfile[256];

     char gfixInfo[256];

     switch ( AlgLattice().FixGaugeKind() ){

     case FIX_GAUGE_NONE:
       sprintf(gfixInfo,"no GF");
       break;

     case FIX_GAUGE_LANDAU:
       sprintf(gfixInfo,"Landau GF, StpCnd=%0.0E", AlgLattice().FixGaugeStopCond());
       break;

     case FIX_GAUGE_COULOMB_T:
       sprintf(gfixInfo,"Coulomb(T) GF, StpCnd=%0.0E", AlgLattice().FixGaugeStopCond());
       break;

     default:
       sprintf(gfixInfo,"UNKNOWN GF");

     }


     char fermionInfo[256];
     
     switch (AlgLattice().Fclass() ){

     case F_CLASS_DWF:
       sprintf(fermionInfo,"DWF, Ls=%i, M5=%0.2f",GJP.Sites(4),GJP.DwfHeight());
       break;
	 
     case F_CLASS_NONE:
       sprintf(fermionInfo,"NO FERMION TYPE");
       break;

     case F_CLASS_STAG:
       sprintf(fermionInfo,"staggered fermion");
       break;

     case F_CLASS_WILSON:
       sprintf(fermionInfo,"Wilson fermion");
       break;

     case F_CLASS_CLOVER:
       sprintf(fermionInfo,"Clover fermion");
       break;
 	
     case F_CLASS_ASQTAD:
       sprintf(fermionInfo,"aSqTad fermion");
       break;

     case F_CLASS_P4: 
       sprintf(fermionInfo,"P4 fermion");
       break;

     default:
       sprintf(fermionInfo,"UNKNOWN FERMION TYPE");
      
     }

       

     sprintf(propType,"4D propagator, mass=%0.3f, StpCond=%0.0E,\nBC=%s%s%s%s,\n%s,\n%s", 
	     qp_arg.cg.mass, qp_arg.cg.stop_rsd,
	     ((GJP.Xbc()==BND_CND_PRD) ? "P" : "A"),
	     ((GJP.Ybc()==BND_CND_PRD) ? "P" : "A"),
	     ((GJP.Zbc()==BND_CND_PRD) ? "P" : "A"),
	     ((GJP.Tbc()==BND_CND_PRD) ? "P" : "A"),   
	     gfixInfo, fermionInfo
	     );
    
     //sprintf(sourceType, "fullSource");
     // be a bit more sophisticated
     sprintf(sourceType,"%s-source at t=%i",SourceType_map[SrcType()].name ,SourceTime());
   
     if(!do_rerun)
       sprintf(propOutfile,qp_arg.file);
     else
       sprintf(propOutfile,"%s.rewrite",qp_arg.file);

     //in case of DWF, renormalize first
     if(AlgLattice().Fclass() == F_CLASS_DWF){
       
       Float renFac = 1./(5. - GJP.DwfHeight());
       
       save_prop = (WilsonMatrix*)smalloc(GJP.VolNodeSites()*sizeof(WilsonMatrix));
       if (save_prop == 0) ERR.Pointer(cname, fname, "pr3op");
       VRB.Smalloc(cname, fname, "save_prop", save_prop,
		   GJP.VolNodeSites() * sizeof(WilsonMatrix));
       
       for(int ii(0); ii <  GJP.VolNodeSites(); ++ii)
	 *(save_prop + ii) = renFac * prop[ii];

       
     }
     else
       save_prop = &prop[0];
     
#ifdef USE_QIO
     Float qio_time = -dclock();
     
     // always writes the full 4D source
     //qio_writePropagator writePropQio(propOutfile, QIO_FULL_SOURCE, save_prop, save_source,
     //			      qp_arg.ensemble_id, qp_arg.ensemble_label, qp_arg.seqNum, propType, sourceType,
     //			      GJP.argc(), GJP.argv(), VOLFMT);
     
     // write a t-slice/slices or hypercube in some cases for the source
     qio_writePropagator writePropQio;

     writePropQio.setHeader(qp_arg.ensemble_id, qp_arg.ensemble_label, qp_arg.seqNum, propType, sourceType);
     
     // just writing one time-slice?
     if( (SrcType() == POINT) || (SrcType() == VOLUME) || (SrcType() == BOX) || (SrcType() == WALL) ){
       VRB.Flow(cname,fname," source-type: %s only write t-slice %i to file\n",SourceType_map[SrcType()].name ,SourceTime());
       writePropQio.setSourceTslice(SourceTime());
     }

     writePropQio.write_12pairs(propOutfile, QIO_FULL_SOURCE, save_prop, save_source, VOLFMT);
     qio_time +=dclock();
     print_time("QPropW::Run","qio_writePropagator",qio_time);


#endif // USE_QIO
     
     if(AlgLattice().Fclass() == F_CLASS_DWF)
       sfree(save_prop);
     
     // the old storage function
     //SaveQProp(qp_arg.file,PROP); 
   }

   
   sfree(save_source);
   
   if(do_rerun){
     
     //now compare prop and read_prop
     
     Float errCnt(0.);
     Float sumerr(0.);
     
     for(int index(0); index < GJP.VolNodeSites(); ++index){
       
       WilsonMatrix mat_read, mat_calc;
       
       mat_calc = prop[index];
       mat_read = *(read_prop + index);
       
       for(int s_src(0); s_src < 4; ++s_src)
	 for(int c_src(0); c_src < 3; ++c_src)
	   for(int s_snk(0); s_snk < 4; ++s_snk)
	     for(int c_snk(0); c_snk < 3; ++c_snk){
	       
	       Complex tmp_calc = mat_calc(s_snk,c_snk,s_src,c_src);
	       Complex tmp_read = mat_read(s_snk,c_snk,s_src,c_src);
	       
	       Float diff;
	       diff = ( fabs(tmp_calc.real()-tmp_read.real()) + fabs(tmp_calc.imag()-tmp_read.imag()) )/ sqrt((tmp_calc.real()*tmp_calc.real() + tmp_calc.imag()*tmp_calc.imag()));
	       
	       if( diff > precision){
		 
		 errCnt += 1.0;
		 sumerr+=diff;
		 VRB.Result(cname,fname,"mismatch propagator: index %i snk %i %i src %i %i\n %f: (%f,%f) <-> (%f,%f)\n",
			index, s_snk,c_snk,s_src,c_src,
			diff, tmp_calc.real(), tmp_calc.imag(), tmp_read.real(), tmp_read.imag() );
	       }
	       
	     }
     }
     
     
     glb_sum_five(&errCnt);
     glb_sum_five(&sumerr);
     Float averr=sumerr/errCnt;
     
     if( fabs(errCnt) > 0.){
       VRB.Result(cname,fname," ReRun prop. with TOTAL NUMBER OF ERRORS: %f\n",errCnt);
       VRB.Result(cname,fname," Average error: %e\n",averr);
       VRB.Result(cname,fname," The precision is set at: %e\n",precision);
     } else {
       VRB.Result(cname,fname," ReRun prop. successfully!\n");
       VRB.Result(cname,fname," The precision is set at: %e\n",precision);
     }
     
   }
   
   if (qp_arg.save_prop || do_rerun )   
     sfree(read_prop);
   
}


// Do conjugate gradient
void QPropW::eig_CG(Vector **V, const int vec_len, Float *M, const int nev, const int m, float **U, Rcomplex *invH, const int def_len, const Float *restart,const int restart_len, FermionVectorTp& source, FermionVectorTp& sol, FermionVectorTp& midsol, int& iter, Float& true_res) 
{

  char *fname = "eig_CG(source&, sol&, midsol&, int&, Float&)";
  VRB.Func(cname, fname);

  /*
  CgArg cg_arg;
  cg_arg.mass = 0.03;
  cg_arg.max_num_iter = 1000;
  cg_arg.stop_rsd = 1.0e-8;
  */

  Lattice& Lat = AlgLattice();

  // Set the node size of the full (non-checkerboarded) fermion field
  //----------------------------------------------------------------
  int ls = GJP.SnodeSites();
  int ls_glb = GJP.SnodeSites()*GJP.Snodes();
  int f_size = GJP.VolNodeSites() * Lat.FsiteSize()/GJP.SnodeSites();
  int f_size_5d = f_size * ls;

  // Do inversion
  //----------------------------------------------------------------
  if (Lat.Fclass() == F_CLASS_DWF) {
    Vector *src_4d    = (Vector*)source.data();
    Vector *sol_4d    = (Vector*)sol.data();
    Vector *midsol_4d = (Vector*)midsol.data();
    Vector *src_5d    = (Vector*)smalloc(f_size_5d * sizeof(IFloat));
    Vector *sol_5d    = (Vector*)smalloc(f_size_5d * sizeof(IFloat));
    if (src_5d == 0)
      ERR.Pointer(cname,fname, "src_5d");
    VRB.Smalloc(cname,fname, "src_5d", src_5d, f_size_5d * sizeof(IFloat));
    if (sol_5d == 0)
      ERR.Pointer(cname,fname, "sol_5d");
    VRB.Smalloc(cname,fname, "sol_5d", sol_5d, f_size_5d * sizeof(IFloat));

    //Get the lattice form the Alg base class
    Lattice& Lat = this->AlgLattice() ;

    Lat.Ffour2five(src_5d, src_4d, 0, ls_glb-1);
    Lat.Ffour2five(sol_5d, sol_4d, ls_glb-1, 0);

	iter = Lat.eig_FmatInv(V, vec_len, M, nev,m,U,invH, def_len, restart, restart_len, sol_5d, src_5d, &(qp_arg.cg), &true_res, CNV_FRM_YES, PRESERVE_NO);
        
   /****************************************************************
     The code below is temporarily isolated for purposes of merging
     with CPS main branch,  12/09/04, Oleg Loktik 
   -------------------- Quarantine starts --------------------------

     iter = Lat.FmatInv4dSrc(sol_5d, src_5d, 0, ls_glb-1, &(Arg.cg), &true_res, 
    	                     CNV_FRM_YES, PRESERVE_NO);
   -------------------- End of quarantine -------------------------*/ 

  //-----------------------------------------------------------------
  // TY Add Start
    if(qp_arg.save_ls_prop) 
       for(int nls(0);nls<GJP.SnodeSites();nls++) SaveQPropLs(sol_5d, qp_arg.file, nls);
    spnclr_cnt++;

    MeasConAxialOld(sol_5d);
  // TY Add End
  //-----------------------------------------------------------------

  //-----------------------------------------------------------------
  // M. Lightman
    MeasJ5qPion(sol_5d);
  // End M. Lightman
  //-----------------------------------------------------------------

    // prop on walls
    Lat.Ffive2four(sol_4d, sol_5d, ls_glb-1, 0);
    // midpoint prop
    if (StoreMidprop())
      Lat.Ffive2four(midsol_4d, sol_5d, ls_glb/2-1, ls_glb/2);

    VRB.Sfree(cname,fname, "sol_5d", sol_5d);
    sfree(sol_5d);
    VRB.Sfree(cname,fname, "src_5d", src_5d);
    sfree(src_5d);

  } else {
    iter = Lat.FmatInv((Vector*)sol.data(),(Vector*)source.data(),
					   &(qp_arg.cg), &true_res, CNV_FRM_YES, PRESERVE_NO);
  }

}
CPS_END_NAMESPACE
