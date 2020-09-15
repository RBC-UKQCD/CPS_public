#ifndef LATTICE_DOUBLE_T_H
#define LATTICE_DOUBLE_T_H

#include<util/lat_cont.h>
#include<util/dwf.h>
CPS_START_NAMESPACE

//Take a lattice resident in memory and double the temporal length by duplication
class LatticeTimeDoubler{
  LatticeContainer lat_orig;
  LRGState lrg_orig;
  size_t lat_size_orig; //units of float
  int t_nodesites_orig;
  int t_sites_orig;
  int t_orig_off; //original base timeslice for this node
  int nodesites_3d;
    
  DoArg do_arg_orig;

  static void docopy(Float *to_base, Float const* from_base, const int copy_toff, const int t_nodesites_orig, const int nodesites_3d){
    for(int f=0;f<GJP.Gparity()+1;f++){
      for(int t_orig_buf=0; t_orig_buf<t_nodesites_orig; t_orig_buf++){
#pragma omp parallel for
	for(int x3d=0; x3d<nodesites_3d;x3d++){
	  Float *to = to_base + 18*4*(x3d + nodesites_3d*(t_orig_buf + copy_toff + GJP.TnodeSites()*f));
	  Float const *from = from_base + 18*4*(x3d + nodesites_3d*(t_orig_buf + t_nodesites_orig*f));
	  memcpy((void*)to, (void*)from, 18*4*sizeof(Float));
	}
      }
    }
  }

  static void fixup_bfm(Lattice &lat){
#ifdef USE_BFM
    bool is_fbfm = (lat.Fclass() == F_CLASS_BFM || lat.Fclass() == F_CLASS_BFM_TYPE2);
    
    //End fbfm bfm instances
    if(is_fbfm){
      if(!UniqueID()){ printf("LatticeTimeDoubler : Fixing up Fbfm\n"); fflush(stdout); }
      Fbfm &fb = dynamic_cast<Fbfm&>(lat);
      fb.bd.end();
      if(Fbfm::use_mixed_solver){
	fb.bf.comm_init();
	fb.bf.end();
      }
    }

    //Fixup QDP
    QDP::multi1d<int> nrow(Nd);  
    for(int i = 0;i<Nd;i++) nrow[i] = GJP.Sites(i);
    QDP::Layout::setLattSize(nrow);
    QDP::Layout::create();
    if(!UniqueID()) printf("LatticeTimeDoubler: QDP reinitialized with latt size (%d,%d,%d,%d) and node size (%d,%d,%d,%d)\n",QDP::Layout::lattSize()[0],QDP::Layout::lattSize()[1],QDP::Layout::lattSize()[2],QDP::Layout::lattSize()[3],
			   QDP::Layout::subgridLattSize()[0],QDP::Layout::subgridLattSize()[1],QDP::Layout::subgridLattSize()[2],QDP::Layout::subgridLattSize()[3]);

    //Reinitialize fbfm bfm instances
    if(is_fbfm){
      Fbfm &fb = dynamic_cast<Fbfm&>(lat);

      for(int i=0;i<2;i++){
	Fbfm::bfm_args[i].node_latt[0] = QDP::Layout::subgridLattSize()[0];
	Fbfm::bfm_args[i].node_latt[1] = QDP::Layout::subgridLattSize()[1];
	Fbfm::bfm_args[i].node_latt[2] = QDP::Layout::subgridLattSize()[2];
	Fbfm::bfm_args[i].node_latt[3] = QDP::Layout::subgridLattSize()[3];
      }	

      fb.bd.init(Fbfm::bfm_args[Fbfm::current_arg_idx]);
      if(Fbfm::use_mixed_solver) {
	fb.bd.comm_end();
	fb.bf.init(Fbfm::bfm_args[Fbfm::current_arg_idx]);
	fb.bf.comm_end();
	fb.bd.comm_init();
      }
      fb.ImportGauge(); //import new gauge field (doesn't apply time BC)
    }

    //Reinitialize FDwf
    if(lat.Fclass() == F_CLASS_DWF){
      void* ptr = lat.FdiracOpInitPtr();
      dwf_end((Dwf *)ptr);
      dwf_init((Dwf *)ptr);
    }


#endif
  }

  void doDoublingSingleNode(Lattice &lat){    //one node in t direction
    Float* lat_p = (Float*)lat.GaugeField();
    Float* lat_orig_p = (Float*)lat_orig.GaugeField();
    for(int half=0;half<2;half++){
      docopy(lat_p, lat_orig_p, 0, t_nodesites_orig,nodesites_3d);
      docopy(lat_p, lat_orig_p, t_nodesites_orig, t_nodesites_orig,nodesites_3d);
    }
  }


  void doDoublingMultiNode(Lattice &lat){ 
    //Copy/double the lattice
    //| 0 1 | 2 3 | 4 5 | 6 7 |  Orig lattice time dir
    //| 0 1 2 3 | 4 5 6 7 | 0 1 2 3 | 4 5 6 7 |  Desired output
    Float* lat_orig_p = (Float*)lat_orig.GaugeField();
    int t_sites_doubled = GJP.TnodeSites()*GJP.Tnodes();
    if(t_sites_doubled != 2*t_sites_orig) ERR.General("LatticeTimeDoubler","doubleLattice","Time size appears not to have doubled!\n");

    int lat_size_dbl = 18*4*GJP.VolNodeSites()*(GJP.Gparity()+1);

    //What orig lat time range does this node need?
    int t_off_need_start = (GJP.TnodeCoor()*GJP.TnodeSites()) % t_sites_orig;   // 0 1 2 3 | 4 5 6 7 | 0 1 2 3 | 4 5 6 7 |
    int t_off_need_lessthan = t_off_need_start + GJP.TnodeSites();

    Float* buf1 = (Float*)malloc(lat_size_orig * sizeof(Float));
    Float* buf2 = (Float*)malloc(lat_size_orig * sizeof(Float));
    memcpy((void*)buf1, (void*)lat_orig.GaugeField(), lat_size_orig*sizeof(Float));

    Float *send = buf2;
    Float *recv = buf1;

    //What orig latt time range is presently in the buffer?
    int t_off_buf_start = t_orig_off;
    int t_off_buf_lessthan = t_orig_off + t_nodesites_orig;
    
    //Do the doubling
    int halves_got[2] = {0,0};

    Float* lat_p = (Float*)lat.GaugeField();

    for(int shift = 0; shift < GJP.Tnodes(); shift++){
      // if(!GJP.XnodeCoor() && !GJP.YnodeCoor() && !GJP.ZnodeCoor() && !GJP.SnodeCoor()){ 
      // 	printf("Shift %d, TnodeCoor %d  buf contains range %d .. %d. I want %d .. %d\n",shift,GJP.TnodeCoor(),t_off_buf_start,t_off_buf_lessthan-1, t_off_need_start,t_off_need_lessthan-1);
      // 	fflush(stdout);
      // }

      //There are 3 cases:
      //1) t_off_buf_start < t_off_need_start || t_off_buf_start >= t_off_need_lessthan    in which case we do nothing
      //2) t_off_buf_start == t_off_need_start  in which case we copy to the first time-half of the nodes doubled lattice
      //3) t_off_buf_start == t_off_need_start + t_nodesites_orig  in which case we copy to the second time-half of the nodes doubled lattice
      if(t_off_buf_start == t_off_need_start && !halves_got[0]){
	// if(!GJP.XnodeCoor() && !GJP.YnodeCoor() && !GJP.ZnodeCoor() && !GJP.SnodeCoor()){ 
	//   printf("TnodeCoor %d buf contains first half of desired range %d .. %d\n",GJP.TnodeCoor(),t_off_need_start,t_off_need_lessthan-1);
	//   fflush(stdout);
	// }
	docopy(lat_p,recv,0,t_nodesites_orig,nodesites_3d);
	halves_got[0] = 1;
      }else if(t_off_buf_start == t_off_need_start + t_nodesites_orig && !halves_got[1]){
	// if(!GJP.XnodeCoor() && !GJP.YnodeCoor() && !GJP.ZnodeCoor() && !GJP.SnodeCoor()){ 
	//   printf("TnodeCoor %d buf contains second half of desired range %d .. %d\n",GJP.TnodeCoor(),t_off_need_start,t_off_need_lessthan-1);
	//   fflush(stdout);
	// }
	docopy(lat_p,recv,t_nodesites_orig,t_nodesites_orig,nodesites_3d);
	halves_got[1] = 1;
      }

      Float* tmp = send;
      send = recv;
      recv = tmp;
      getPlusData(recv, send, lat_size_orig, 3); //send left
      t_off_buf_start = (t_off_buf_start + t_nodesites_orig) % t_sites_orig;
      t_off_buf_lessthan = t_off_buf_start + t_nodesites_orig;
    }
    Float hv_count = 2-halves_got[0]-halves_got[1];
    glb_sum_five(&hv_count);
    if(hv_count != 0.0) ERR.General("LatticeTimeDoubler","doubleLattice","Some nodes are missing halves!\n");
    
    free(buf1);
    free(buf2);
  }


  void printSanityCheck(Lattice &lattice){
    int lat_size = 18*4*GJP.VolNodeSites()*(GJP.Gparity()+1);
    int Lt = GJP.TnodeSites()*GJP.Tnodes();

    for(int f=0;f<GJP.Gparity()+1;f++){
      Float zerothelems[Lt/2];
      for(int i=0;i<Lt/2;i++) zerothelems[i] = 0.;

      if(!GJP.XnodeCoor() && !GJP.YnodeCoor() && !GJP.ZnodeCoor() && !GJP.SnodeCoor()){
	for(int t=0;t<GJP.TnodeSites()/2;t++){
	  int off = 18*4*GJP.XnodeSites()*GJP.YnodeSites()*GJP.ZnodeSites()*(t + GJP.TnodeSites()/2*f);
	  zerothelems[GJP.TnodeSites()*GJP.TnodeCoor()/2+t] = ((Float const*)OrigGaugeField())[off];
	}
      }
      if(!UniqueID()) printf("LatticeTimeDoubler : Original elements with (x,y,z)=(0,0,0) along t direction with flavor %d\n",f);
      for(int t=0;t<Lt/2;t++){
	glb_sum_five(&zerothelems[t]);
	if(!UniqueID()){
	  if(t % (GJP.TnodeSites()/2) == 0) printf("------\n");
	  printf("%d %g\n",t,zerothelems[t]);
	}
      }
    }

    for(int f=0;f<GJP.Gparity()+1;f++){
      Float zerothelems[Lt];
      for(int i=0;i<Lt;i++) zerothelems[i] = 0.;

      if(!GJP.XnodeCoor() && !GJP.YnodeCoor() && !GJP.ZnodeCoor() && !GJP.SnodeCoor()){
	for(int t=0;t<GJP.TnodeSites();t++){
	  int off = 18*4*GJP.XnodeSites()*GJP.YnodeSites()*GJP.ZnodeSites()*(t + GJP.TnodeSites()*f);
	  zerothelems[GJP.TnodeSites()*GJP.TnodeCoor()+t] = ((Float*)lattice.GaugeField())[off];
	}
      }
      if(!UniqueID()) printf("LatticeTimeDoubler : Elements with (x,y,z)=(0,0,0) along t direction (first half second half) with flavor %d\n",f);
      for(int t=0;t<Lt/2;t++){
	glb_sum_five(&zerothelems[t]);
	glb_sum_five(&zerothelems[t+Lt/2]);
	if(!UniqueID()){ 
	  if(t % GJP.TnodeSites() == 0) printf("------\n");
	  printf("%d %g | %d %g\n",t,zerothelems[t],t+Lt/2,zerothelems[t+Lt/2]);
	}
      }
    }
  }



 public:

  Matrix const* OrigGaugeField(){
    return lat_orig.GaugeField();
  }

  void doubleLattice(Lattice &lat, DoArg do_arg){
    if(GJP.Bc(3) != BND_CND_PRD) ERR.General("LatticeTimeDoubler","doubleLattice","Not implemented for non-periodic fermion temporal BCs\n");

    //First store a copy of the current state
    lat_orig.Get(lat);
    lrg_orig.GetStates();
    lat_size_orig = 18*4*GJP.VolNodeSites()*(GJP.Gparity()+1);
    do_arg_orig = do_arg;

    t_nodesites_orig = GJP.TnodeSites();
    t_sites_orig = GJP.TnodeSites()*GJP.Tnodes();
    t_orig_off = GJP.TnodeSites()*GJP.TnodeCoor(); //original base timeslice for this node
    nodesites_3d = GJP.VolNodeSites()/GJP.TnodeSites();
    
    //Setup double-t environment
    lat.FreeGauge(); //Free the lattice memory. We don't need to free the LRG as I have a separate function that re-initializes that
    do_arg.t_sites *= 2; //note this does not change the stored do_arg of the main program

    GJP.Initialize(do_arg);
    lat.AllocGauge(); //double-length gauge field, uninitialized
    LRG.Reinitialize();

    if(GJP.Tnodes()>1) doDoublingMultiNode(lat);
    else doDoublingSingleNode(lat);

    fixup_bfm(lat);
    printSanityCheck(lat);
  }

};


CPS_END_NAMESPACE
#endif

