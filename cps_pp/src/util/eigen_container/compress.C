#include<util/verbose.h>
#include<util/eig_io.h>
#define MAX_EVEC_PRINT_NORM 8
#define FP16_WIDTH_24 2.0833333333333
#define FP16_COEF_EXP_SHARE_FLOATS 10
#define FP16_WIDTH_COEF (2.0*(1.0 + 1.0 / (double)FP16_COEF_EXP_SHARE_FLOATS))
// Adopting eig-comp from C. Lehner, based on multi-grid lanczos arxiv:1710.06884

CPS_START_NAMESPACE

// prec = sizeof(float or double)
int EvecWriter::writeCompressedVector 
( const char *dir, OPT *V, struct evec_write &warg, std::vector<Float>  &evals){

//  const char *cname="";
  const char *fname="";
  static Timer timer (cname,fname);
  timer.start();
  int neig = evals.size();
  int ID = UniqueID();
#pragma omp parallel
  {
#pragma omp single
    {
      nthreads = omp_get_num_threads();
    }
  }
  float *raw_in = (float*)V;

  VRB.Result(cname,fname,"%s\n%d threads\n\n",header,nthreads);
//  size_t vol4d, vol5d;
  size_t f_size, f_size_coef_block, nkeep_fp16;
//, f_size_block, 

//    fprintf(stderr,"Arguments: sx sy sz st s5 bx by bz bt b5 nkeep fileindex filesperdir bigendian nkeep_single_prec vrb_nkeep_res vrb_evec_res\n");

  {
    for(int i=0;i<5;i++) args.s[i]=GJP.NodeSites(i);
    for(int i=0;i<5;i++) args.b[i]=warg.b[i];
    args.nkeep=warg.nkeep;
    args.nkeep_single=warg.nkeep_single;
//    args.filesperdir=warg.filesperdir;
    warg.findex=ID;
 
   if(warg.bigendian) bigendian=true;

    VRB.Result(cname,fname,"Parameters:\n");
    for (int i=0;i<5;i++)
      VRB.Result(cname,fname,"s[%d] = %d\n",i,args.s[i]);
    for (int i=0;i<5;i++)
      VRB.Result(cname,fname,"b[%d] = %d\n",i,args.b[i]);
    VRB.Result(cname,fname,"nkeep = %d\n",args.nkeep);
    VRB.Result(cname,fname,"file_index = %10.10d\n",warg.findex);
    VRB.Result(cname,fname,"n_dir = %d\n",warg.n_dir);
    VRB.Result(cname,fname,"big_endian = %d\n",warg.bigendian);
    VRB.Result(cname,fname,"nkeep_single = %d\n",args.nkeep_single);

    vol4d = args.s[0] * args.s[1] * args.s[2] * args.s[3];
    vol5d = vol4d * args.s[4];
    f_size = vol5d / 2 * 24;

    VRB.Result(cname,fname,"f_size = %d\n",f_size);

    // sanity check
    args.blocks = 1;
    for (int i=0;i<5;i++) {
      if (args.s[i] % args.b[i]) {
	fprintf(stderr,"Invalid blocking in dimension %d\n",i);
	return 72;
      }

      args.nb[i] = args.s[i] / args.b[i];
      args.blocks *= args.nb[i];
    }

    f_size_block = f_size / args.blocks;

    VRB.Result(cname,fname,"number of blocks = %d\n",args.blocks);
    VRB.Result(cname,fname,"f_size_block = %d\n",f_size_block);

    VRB.Result(cname,fname,"Internally using sizeof(OPT) = %d\n",sizeof(OPT));

  }

  {
//    off_t size = ftello(f);

//    if ( size % (sizeof(float)*f_size) ) {
//      fprintf(stderr,"Invalid file size\n");
//      return 2;
//    }

//    neig = ( size / sizeof(float) / f_size );
    size_t size = f_size*sizeof(float)*neig;

    f_size_coef_block = neig * 2 * args.nkeep;

    VRB.Result(cname,fname,"We have %d eigenvectors \n",neig);

    printf("Size of operating coefficient data in GB: %g\n", (double)f_size_coef_block * (double)args.blocks / 1024./1024./1024. * sizeof(OPT));

    nkeep_fp16 = args.nkeep - args.nkeep_single;
    if (nkeep_fp16 < 0)
      nkeep_fp16 = 0;
    VRB.Result(cname,fname,"neig nkeep_single nkee_fp16=%d %d %d\n",neig,  args.nkeep_single, nkeep_fp16);

    // estimate of compression

    {
      double size_of_coef_data = (neig * (double)args.blocks * (args.nkeep_single * 4 + nkeep_fp16 * FP16_WIDTH_COEF)*2)  / 1024. / 1024. / 1024.;
      double size_of_evec_data = ((args.nkeep - nkeep_fp16)* f_size * 4 + nkeep_fp16 * f_size * FP16_WIDTH_24)  / 1024. / 1024. / 1024.;
      double size_orig = (double)size  / 1024. / 1024. / 1024.;
      double size_of_comp = size_of_coef_data+size_of_evec_data;
      printf("--------------------------------------------------------------------------------\n");
//      printf("Original size:     %g GB\n",size_orig);
      printf("Compressed size:   %g GB  (%g GB coef, %g GB evec)\n",size_of_comp,size_of_coef_data,size_of_evec_data);
      printf("Compressed to %.4g%% of original\n",size_of_comp / size_orig * 100.);
      printf("--------------------------------------------------------------------------------\n");
    }

    double t1 = dclock();

//    double size_in_gb = f_size * (double)neig * sizeof(float) / 1024. / 1024. / 1024.;

//    printf("Read %.4g GB in %.4g seconds at %.4g GB/s\n", size_in_gb, t1-t0,size_in_gb / (t1-t0) );

//    uint32_t crc_comp = crc32_fast(raw_in,(size_t)f_size * neig * sizeof(float),0);
    uint32_t crc_comp = crc32_loop(0,(unsigned char *)raw_in,(size_t)f_size * neig * sizeof(float));

    double t2 = dclock();

    printf("Node %d: Computed CRC32: %X   (in %.4g seconds)\n",ID,crc_comp,t2-t1);

//    fclose(f);

    VRB.Result(cname,fname,"Fixing endian-ness\n");

    // fix endian if needed
// turning off threading for sum testing. Eventually unnecessary 
//#pragma omp parallel for
    for (int j=0;j<neig;j++) {
      float* segm = &raw_in[ (size_t)f_size * j ];
      fix_float_endian(segm, f_size);

        if (j<MAX_EVEC_PRINT_NORM) {
        Float sum=sp_single(segm,segm,f_size).real();
        glb_sum(&sum);
	printf("Norm %d: %g\n",j, sum);
      }
    }


  }

  //
  // Status: 
  //  loaded uncompressed eigenvector slot and verified it
  //

    double t0 = dclock();
  {

    // create block memory
    block_data.resize(args.blocks);
    for (int i=0;i<args.blocks;i++)
      block_data[i].resize(f_size_block * neig);    

    
    //
#pragma omp parallel 
    {
      for (int nev=0;nev<neig;nev++) {
	float* raw_in_ev = &raw_in[ (size_t)f_size * nev ];
	
//        Float norm=0.;// check as CJ is not smart enough
#pragma omp for 
	for (size_t idx=0;idx<vol4d;idx++) {
	  int pos[5], pos_in_block[5], block_coor[5];
	  index_to_pos(idx,pos,args.s);
	  
	  int parity = (pos[0] + pos[1] + pos[2] + pos[3]) % 2;
	  if (parity == 1) {
	    
	    for (pos[4]=0;pos[4]<args.s[4];pos[4]++) {
	      pos_to_blocked_pos(pos,pos_in_block,block_coor);
	      
	      int bid = pos_to_index(block_coor, args.nb);
	      int ii = pos_to_index(pos_in_block, args.b) / 2;

	      OPT* dst = &block_data[bid][ ii*24 + (size_t)f_size_block * nev ];
	      
	      int co;
	      for (co=0;co<12;co++) {
		float* in=&raw_in_ev[ get_cps_index(pos,co,args.s) ];
		dst[2*co + 0] = in[0]; // may convert precision depending on OPT
		dst[2*co + 1] = in[1];
//		norm += in[0]*in[0]+in[1]*in[1];
	      }
	    }
	  }
	}
//        glb_sum(&norm);
//        VRB.Result(cname,fname,"norm %d=%e\n",nev,norm);
      }
    }

    double t1 = dclock();

    printf("Created block structure in %.4g seconds\n",t1-t0);

    // simple test
    for(int test_ev=0;test_ev<neig;test_ev++)
    {
//      int test_ev = neig - 1;
      float* raw_in_ev = &raw_in[ (size_t)f_size * test_ev ];
      float nrm = sp(raw_in_ev,raw_in_ev,f_size).real();

      int i;
      double nrm_blocks = 0.0;
      for (i=0;i<args.blocks;i++) {
	OPT* in_ev = &block_data[i][ (size_t)f_size_block * test_ev ];
	nrm_blocks += sp(in_ev,in_ev,f_size_block).real();
      }

      VRB.Debug(cname,fname,"Difference of norms after blocking: %g - %g = %g\n",nrm,nrm_blocks,nrm-nrm_blocks);

      if (fabs(nrm - nrm_blocks) > 1e-5) {
	fprintf(stderr,"Unexpected error in creating blocks\n");
//	return 91;
      }
    }
  }

  // Now do Gram-Schmidt
  {
    // create block memory
    block_data_ortho.resize(args.blocks);
    for (int i=0;i<args.blocks;i++)
      block_data_ortho[i].resize(f_size_block * args.nkeep);    
    
    double t0 = dclock();

    int nevmax = args.nkeep;

    double flops = 0.0;
    double bytes = 0.0;
#pragma omp parallel shared(flops,bytes)
    {
      double flopsl = 0.0;
      double bytesl = 0.0;
#define COUNT_FLOPS_BYTES(f,b) flopsl += (f) + 1; bytesl += (b) + 2;
      // #define COUNT_FLOPS_BYTES(f,b)

#pragma omp for
      for (int nb=0;nb<args.blocks;nb++) {
	
	for (int iev=0;iev<nevmax;iev++) {
	  
	  OPT* orig = &block_data[nb][ (int64_t)f_size_block * iev ];
	  OPT* res = &block_data_ortho[nb][ (int64_t)f_size_block * iev ];
	  
	  memcpy(res,orig,sizeof(OPT)*f_size_block); 	  COUNT_FLOPS_BYTES(f_size_block,2*f_size_block*sizeof(OPT));

	  for (int jev=0;jev<iev;jev++) {
	    
	    OPT* ev_j = &block_data_ortho[nb][ (int64_t)f_size_block * jev ];
	    
	    // res = |i> - <j|i> |j>
	    // <j|res>
	    std::complex<OPT> res_j = sp_single(ev_j,res,f_size_block);  
	    COUNT_FLOPS_BYTES(8 / 2 * f_size_block, 2*f_size_block*sizeof(OPT)); // 6 per complex multiply, 2 per complex add -> 8 / 2 = 4

	    caxpy_single(res,- res_j,ev_j,res,f_size_block); 
	    COUNT_FLOPS_BYTES(8 / 2 * f_size_block, 3*f_size_block*sizeof(OPT));
	  }
	  
	  // normalize
	  std::complex<OPT> nrm = sp_single(res,res,f_size_block); 
	  COUNT_FLOPS_BYTES(8 / 2 * f_size_block,2*f_size_block*sizeof(OPT));

	  scale_single(res, (OPT)(1.0 / sqrt(nrm.real())),f_size_block); 
	  COUNT_FLOPS_BYTES(f_size_block,2*f_size_block*sizeof(OPT));
	  
	}
      }

#pragma omp critical
      {
	flops += flopsl + 1;
	bytes += bytesl + 2;
      }
    }
    
    double t1 = dclock();

    printf("Gram-Schmidt took %.4g seconds (%g Gflops/s, %g GB/s)\n",t1-t0,flops / (t1-t0) / 1000./1000./1000.,
	   bytes / (t1-t0) / 1024./1024./1024.);

  }


  // Get coefficients and create graphs
  {
    // create block memory
    block_coef.resize(args.blocks);
    for (int i=0;i<args.blocks;i++)
      block_coef[i].resize(f_size_coef_block);    

//    double t0 = dclock();

    if (!warg.vrb_nkeep_res && !warg.vrb_evec_res) {
      printf("Do not display convergence, use fast codepath for obtaining coefficients\n");

#pragma omp parallel for
      for (int nb=0;nb<args.blocks;nb++) {
	for (int j=0;j<neig;j++) {
	  for (int i=0;i<args.nkeep;i++) {
//      printf("Do not display convergence, use fast codepath for obtaining coefficients %d %d %d\n",nb,j,i);
	    get_coef(nb,i,j);
	  }
	}
      }

    } else {

      printf("Slow codepath to display convergence\n");
      
      for (int j=0;j<neig;j++) {
	
	double norm_j = -1.0;
	
	// only compute norm if needed for verbosity
	if ((j < args.nkeep && warg.vrb_nkeep_res < args.nkeep) ||
	    !(j % warg.vrb_evec_res))
	  norm_j = norm_of_evec(block_data,j);
	
	for (int i=0;i<args.nkeep;i++) {
//      printf("Slow codepath to display convergence %d %d\n",j,i);
	  
	  if (i == j && !(i % warg.vrb_nkeep_res))
	    printf("Node %d: nkeep_residuum %d = %g\n",ID,i,norm_of_evec(block_data,j) / norm_j);
	  
//#pragma omp parallel for
	  for (int nb=0;nb<args.blocks;nb++) {
	    get_coef(nb,i,j);
	  }
	  
	}
	
	if (!(j % warg.vrb_evec_res))
	  printf("Node %d: evec_residuum %d = %g\n",ID,j,norm_of_evec(block_data,j) / norm_j);
      }
    }



  }
  double t1 = dclock();
  printf("Node %d: Computing block-coefficients took %.4g seconds\n",ID,t1-t0);

  // write result
  {
    uint32_t crc = 0x0;
    off_t begin_fp16_evec;
    off_t begin_coef;
    char buf[1024];
    
//    cps::sync();
    // write data
    int concur = warg.concur;
    if(concur<1) concur=1;
    int n_block = GJP.TotalNodes()/concur;
    if(n_block<1) n_block=1;
    int fperdir =  GJP.TotalNodes()/warg.n_dir;
    if (GJP.TotalNodes()%warg.n_dir) fperdir +=1;
    if(fperdir<1 ) fperdir=1;

    Float temp=1.;
    VRB.Result(cname,fname,"concur=%d n_block=%d file_per_dir=%d \n",concur,n_block,fperdir);
    if(!ID) mkdir(dir,0755);
      cps::sync();
    glb_sum(&temp);
//for(int i = 0; i< GJP.TotalNodes();i++){
for(int i = 0; i< n_block;i++){
      cps::sync();
    glb_sum(&temp);
    if( (ID%n_block)==i ) 
    {
      int n_dir = warg.findex/fperdir;
      sprintf(buf,"%s/%2.2d",dir,n_dir);
//      sprintf(buf,"%s/%2.2d",dir,warg.findex / warg.filesperdir,warg.findex);
      mkdir(buf,0755);
      sprintf(buf,"%s/%2.2d/%10.10d.compressed",dir,n_dir,warg.findex);
      printf("writing to %s/%2.2d/%10.10d.compressed\n",dir,n_dir,warg.findex);
      FILE* f = fopen(buf,"w+b");
      if (!f) {
	ERR.General(cname,fname,"Could not open %s for writing!\n",buf);
//	return 1;
      }
      
//      int nb;
      
      size_t _t = (size_t)f_size_block * (args.nkeep - nkeep_fp16);
      for (int nb=0;nb<args.blocks;nb++){
//        printf("Node %d: write_floats %d\n",ID,nb);
	write_floats(f,crc,  &block_data_ortho[nb][0], _t );
      }
      
      begin_fp16_evec = ftello(f);
      
      for (int nb=0;nb<args.blocks;nb++){
//        printf("Node %d: write_floats_fp16 %d\n",ID,nb);
	write_floats_fp16(f,crc,  &block_data_ortho[nb][ _t ], (int64_t)f_size_block * nkeep_fp16, 24 );
      }
      
      begin_coef = ftello(f);

      // write coefficients of args.nkeep_single as floats, higher coefficients as fp16
      
      
      for (int j=0;j<neig;j++)
	for (int nb=0;nb<args.blocks;nb++) {
//        printf("Node %d: write_floats %d %d \n",ID,j,nb);
	  write_floats(f,crc,  &block_coef[nb][2*args.nkeep*j], 2*(args.nkeep - nkeep_fp16) );
	  write_floats_fp16(f,crc,  &block_coef[nb][2*args.nkeep*j + 2*(args.nkeep - nkeep_fp16) ], 2*nkeep_fp16 , FP16_COEF_EXP_SHARE_FLOATS);
	}
      
      fclose(f);
      printf("writing to %s/%2.2d/%10.10d.compressed done \n",dir,n_dir,warg.findex);
    }
}
      cps::sync();
    glb_sum(&temp);

if(!ID){
      sprintf(buf,"%s/eigen-values.txt",dir);
      FILE* f = fopen(buf,"wt");
      fprintf(f,"%d\n",neig);
      for(int i=0;i<neig;i++)
      fprintf(f,"%0.14E\n",evals[i]);
      fclose(f);
}
      cps::sync();
    glb_sum(&temp);
     
      
    // write meta data
    {
//      sprintf(buf,"%s/%2.2d/%10.10d.meta",dir,warg.findex / warg.filesperdir,warg.findex);
      sprintf(buf,"%s/metadata.txt",dir);
for(int node = 0; node< GJP.TotalNodes();node++){
      cps::sync();
    glb_sum(&temp);
if (node==ID){
//printf("Node %d got the slot\n",i);
if (!node){
      FILE* f = fopen(buf,"wt");
      if (!f) {
	fprintf(stderr,"Could not open %s for writing!\n",buf);
	return 1;
      }

      for (int i=0;i<5;i++)
	fprintf(f,"s[%d] = %d\n",i,args.s[i]);
      for (int i=0;i<5;i++)
	fprintf(f,"b[%d] = %d\n",i,args.b[i]);
      for (int i=0;i<5;i++)
	fprintf(f,"nb[%d] = %d\n",i,args.nb[i]);
      fprintf(f,"neig = %d\n",neig);
      fprintf(f,"nkeep = %d\n",args.nkeep);
      fprintf(f,"nkeep_single = %d\n",args.nkeep_single);
      fprintf(f,"blocks = %d\n",args.blocks);
      fprintf(f,"FP16_COEF_EXP_SHARE_FLOATS = %d\n",FP16_COEF_EXP_SHARE_FLOATS);
//      fprintf(f,"crc32[%d] = %X\n",warg.findex,crc);
//      fprintf(f,"index = %d\n",warg.findex);

      fclose(f);
}
      FILE* f = fopen(buf,"a");
      fprintf(f,"crc32[%d] = %X\n",node,crc);
      fclose(f);
      if(VRB.Level(VERBOSE_FLOW_LEVEL))
      printf("Node %d done\n",node);

}
}
    }
    Float t2= dclock();
    printf("Node %d: File I/O took %.4g seconds\n",ID,t2-t1);
  }

  // Cleanup
//  free(raw_in);

  timer.stop(true);
  return 0;
}


CPS_END_NAMESPACE
