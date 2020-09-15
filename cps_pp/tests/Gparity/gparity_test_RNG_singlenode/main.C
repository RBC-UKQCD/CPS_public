#include<config.h>

#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>

#include<util/lattice.h>
#include<util/random.h>
#include<util/time_cps.h>

#include<alg/alg_hmc.h>
#include<alg/common_arg.h>
#include<alg/hmc_arg.h>
#include<alg/hmd_arg.h>

#include<alg/alg_int.h>
#include<alg/int_arg.h>

#include<alg/no_arg.h>
#include<alg/do_arg.h>
#include<alg/alg_plaq.h>
#include<alg/alg_pbp.h>
#include<alg/pbp_arg.h>
#include<alg/alg_remez.h>
#include<alg/alg_wline.h>

#include<util/gjp.h>
#include<util/verbose.h>
#include<util/error.h>
#include<util/qcdio.h>
#include<util/WriteLatticePar.h>
#include<util/ReadLatticePar.h>
#include<util/qioarg.h>
#include<util/command_line.h>
//#include <sys/bgl/bgl_sys_all.h>

#undef USE_SCU_CHECKSUMS
#ifdef USE_SCU_CHECKSUMS
#include <qcdocos/scu_checksum.h>
#endif

#include<sstream>
#include<algorithm>
#include<omp.h>
#include<cassert>
#include<set>
//--------------------------------------------------------------

USING_NAMESPACE_CPS
using namespace std;

DoArg do_arg;

void setupDoArg(DoArg &do_arg, const int size[], const int ngp){
  do_arg.x_sites = size[0];
  do_arg.y_sites = size[1];
  do_arg.z_sites = size[2];
  do_arg.t_sites = size[3];
  do_arg.s_sites = size[4];
  do_arg.x_node_sites = 0;
  do_arg.y_node_sites = 0;
  do_arg.z_node_sites = 0;
  do_arg.t_node_sites = 0;
  do_arg.s_node_sites = 0;
  do_arg.x_nodes = 0;
  do_arg.y_nodes = 0;
  do_arg.z_nodes = 0;
  do_arg.t_nodes = 0;
  do_arg.s_nodes = 0;
  do_arg.updates = 0;
  do_arg.measurements = 0;
  do_arg.measurefreq = 0;
  do_arg.cg_reprod_freq = 10;
  do_arg.x_bc = BND_CND_PRD;
  do_arg.y_bc = BND_CND_PRD;
  do_arg.z_bc = BND_CND_PRD;
  do_arg.t_bc = BND_CND_APRD;
  do_arg.start_conf_kind = START_CONF_ORD;
  do_arg.start_conf_load_addr = 0x0;
  do_arg.start_seed_kind = START_SEED_FIXED;
  do_arg.start_seed_filename = "../rngs/ckpoint_rng.0";
  do_arg.start_conf_filename = "../configurations/ckpoint_lat.0";
  do_arg.start_conf_alloc_flag = 6;
  do_arg.wfm_alloc_flag = 2;
  do_arg.wfm_send_alloc_flag = 2;
  do_arg.start_seed_value = 83209;
  do_arg.beta =   2.25;
  do_arg.c_1 =   -3.3100000000000002e-01;
  do_arg.u0 =   1.0000000000000000e+00;
  do_arg.dwf_height =   1.8000000000000000e+00;
  do_arg.dwf_a5_inv =   1.0000000000000000e+00;
  do_arg.power_plaq_cutoff =   0.0000000000000000e+00;
  do_arg.power_plaq_exponent = 0;
  do_arg.power_rect_cutoff =   0.0000000000000000e+00;
  do_arg.power_rect_exponent = 0;
  do_arg.verbose_level = -1202;//VERBOSE_DEBUG_LEVEL; //-1202;
  do_arg.checksum_level = 0;
  do_arg.exec_task_list = 0;
  do_arg.xi_bare =   1.0000000000000000e+00;
  do_arg.xi_dir = 3;
  do_arg.xi_v =   1.0000000000000000e+00;
  do_arg.xi_v_xi =   1.0000000000000000e+00;
  do_arg.clover_coeff =   0.0000000000000000e+00;
  do_arg.clover_coeff_xi =   0.0000000000000000e+00;
  do_arg.xi_gfix =   1.0000000000000000e+00;
  do_arg.gfix_chkb = 1;
  do_arg.asqtad_KS =   0.0000000000000000e+00;
  do_arg.asqtad_naik =   0.0000000000000000e+00;
  do_arg.asqtad_3staple =   0.0000000000000000e+00;
  do_arg.asqtad_5staple =   0.0000000000000000e+00;
  do_arg.asqtad_7staple =   0.0000000000000000e+00;
  do_arg.asqtad_lepage =   0.0000000000000000e+00;
  do_arg.p4_KS =   0.0000000000000000e+00;
  do_arg.p4_knight =   0.0000000000000000e+00;
  do_arg.p4_3staple =   0.0000000000000000e+00;
  do_arg.p4_5staple =   0.0000000000000000e+00;
  do_arg.p4_7staple =   0.0000000000000000e+00;
  do_arg.p4_lepage =   0.0000000000000000e+00;

  BndCndType* bcs[3] = { &do_arg.x_bc , &do_arg.y_bc, &do_arg.z_bc };
  if(ngp > 3 || ngp < 0){
    ERR.General("","","Number of GPBC must be between 0 and 3\n");
  }
  for(int i=0;i<ngp;i++)
    *bcs[i] = BND_CND_GPARITY;

}


static const uint32_t kCrc32Table[256] = {
    0x00000000, 0x77073096, 0xee0e612c, 0x990951ba,
    0x076dc419, 0x706af48f, 0xe963a535, 0x9e6495a3,
    0x0edb8832, 0x79dcb8a4, 0xe0d5e91e, 0x97d2d988,
    0x09b64c2b, 0x7eb17cbd, 0xe7b82d07, 0x90bf1d91,
    0x1db71064, 0x6ab020f2, 0xf3b97148, 0x84be41de,
    0x1adad47d, 0x6ddde4eb, 0xf4d4b551, 0x83d385c7,
    0x136c9856, 0x646ba8c0, 0xfd62f97a, 0x8a65c9ec,
    0x14015c4f, 0x63066cd9, 0xfa0f3d63, 0x8d080df5,
    0x3b6e20c8, 0x4c69105e, 0xd56041e4, 0xa2677172,
    0x3c03e4d1, 0x4b04d447, 0xd20d85fd, 0xa50ab56b,
    0x35b5a8fa, 0x42b2986c, 0xdbbbc9d6, 0xacbcf940,
    0x32d86ce3, 0x45df5c75, 0xdcd60dcf, 0xabd13d59,
    0x26d930ac, 0x51de003a, 0xc8d75180, 0xbfd06116,
    0x21b4f4b5, 0x56b3c423, 0xcfba9599, 0xb8bda50f,
    0x2802b89e, 0x5f058808, 0xc60cd9b2, 0xb10be924,
    0x2f6f7c87, 0x58684c11, 0xc1611dab, 0xb6662d3d,
    0x76dc4190, 0x01db7106, 0x98d220bc, 0xefd5102a,
    0x71b18589, 0x06b6b51f, 0x9fbfe4a5, 0xe8b8d433,
    0x7807c9a2, 0x0f00f934, 0x9609a88e, 0xe10e9818,
    0x7f6a0dbb, 0x086d3d2d, 0x91646c97, 0xe6635c01,
    0x6b6b51f4, 0x1c6c6162, 0x856530d8, 0xf262004e,
    0x6c0695ed, 0x1b01a57b, 0x8208f4c1, 0xf50fc457,
    0x65b0d9c6, 0x12b7e950, 0x8bbeb8ea, 0xfcb9887c,
    0x62dd1ddf, 0x15da2d49, 0x8cd37cf3, 0xfbd44c65,
    0x4db26158, 0x3ab551ce, 0xa3bc0074, 0xd4bb30e2,
    0x4adfa541, 0x3dd895d7, 0xa4d1c46d, 0xd3d6f4fb,
    0x4369e96a, 0x346ed9fc, 0xad678846, 0xda60b8d0,
    0x44042d73, 0x33031de5, 0xaa0a4c5f, 0xdd0d7cc9,
    0x5005713c, 0x270241aa, 0xbe0b1010, 0xc90c2086,
    0x5768b525, 0x206f85b3, 0xb966d409, 0xce61e49f,
    0x5edef90e, 0x29d9c998, 0xb0d09822, 0xc7d7a8b4,
    0x59b33d17, 0x2eb40d81, 0xb7bd5c3b, 0xc0ba6cad,
    0xedb88320, 0x9abfb3b6, 0x03b6e20c, 0x74b1d29a,
    0xead54739, 0x9dd277af, 0x04db2615, 0x73dc1683,
    0xe3630b12, 0x94643b84, 0x0d6d6a3e, 0x7a6a5aa8,
    0xe40ecf0b, 0x9309ff9d, 0x0a00ae27, 0x7d079eb1,
    0xf00f9344, 0x8708a3d2, 0x1e01f268, 0x6906c2fe,
    0xf762575d, 0x806567cb, 0x196c3671, 0x6e6b06e7,
    0xfed41b76, 0x89d32be0, 0x10da7a5a, 0x67dd4acc,
    0xf9b9df6f, 0x8ebeeff9, 0x17b7be43, 0x60b08ed5,
    0xd6d6a3e8, 0xa1d1937e, 0x38d8c2c4, 0x4fdff252,
    0xd1bb67f1, 0xa6bc5767, 0x3fb506dd, 0x48b2364b,
    0xd80d2bda, 0xaf0a1b4c, 0x36034af6, 0x41047a60,
    0xdf60efc3, 0xa867df55, 0x316e8eef, 0x4669be79,
    0xcb61b38c, 0xbc66831a, 0x256fd2a0, 0x5268e236,
    0xcc0c7795, 0xbb0b4703, 0x220216b9, 0x5505262f,
    0xc5ba3bbe, 0xb2bd0b28, 0x2bb45a92, 0x5cb36a04,
    0xc2d7ffa7, 0xb5d0cf31, 0x2cd99e8b, 0x5bdeae1d,
    0x9b64c2b0, 0xec63f226, 0x756aa39c, 0x026d930a,
    0x9c0906a9, 0xeb0e363f, 0x72076785, 0x05005713,
    0x95bf4a82, 0xe2b87a14, 0x7bb12bae, 0x0cb61b38,
    0x92d28e9b, 0xe5d5be0d, 0x7cdcefb7, 0x0bdbdf21,
    0x86d3d2d4, 0xf1d4e242, 0x68ddb3f8, 0x1fda836e,
    0x81be16cd, 0xf6b9265b, 0x6fb077e1, 0x18b74777,
    0x88085ae6, 0xff0f6a70, 0x66063bca, 0x11010b5c,
    0x8f659eff, 0xf862ae69, 0x616bffd3, 0x166ccf45,
    0xa00ae278, 0xd70dd2ee, 0x4e048354, 0x3903b3c2,
    0xa7672661, 0xd06016f7, 0x4969474d, 0x3e6e77db,
    0xaed16a4a, 0xd9d65adc, 0x40df0b66, 0x37d83bf0,
    0xa9bcae53, 0xdebb9ec5, 0x47b2cf7f, 0x30b5ffe9,
    0xbdbdf21c, 0xcabac28a, 0x53b39330, 0x24b4a3a6,
    0xbad03605, 0xcdd70693, 0x54de5729, 0x23d967bf,
    0xb3667a2e, 0xc4614ab8, 0x5d681b02, 0x2a6f2b94,
    0xb40bbe37, 0xc30c8ea1, 0x5a05df1b, 0x2d02ef8d,
}; // kCrc32Table

class Crc32
{
public:
    Crc32() { Reset(); }
    ~Crc32() throw() {}
    void Reset() { _crc = (uint32_t)~0; }
    void AddData(const uint8_t* pData, const uint32_t length)
    {
        uint8_t* pCur = (uint8_t*)pData;
        uint32_t remaining = length;
        for (; remaining--; ++pCur)
            _crc = ( _crc >> 8 ) ^ kCrc32Table[(_crc ^ *pCur) & 0xff];
    }
    const uint32_t GetCrc32() { return ~_crc; }

private:
    uint32_t _crc;
};








bool RNG_equal(UGrandomGenerator &a, UGrandomGenerator &b){
  // int aa[57], bb[57];
  // a.store(aa); b.store(bb);
  int const *aa = a.getState();
  int const *bb = b.getState();
  
  for(int i=0;i<55;i++) if(aa[i] != bb[i]) return false;
  return true;
}

int main(int argc, char *argv[])
{ 
  Start(&argc,&argv);
  const char *cname=argv[0];
  const char *fname="main()";

  int ngp;
  {
    std::stringstream ss; ss << argv[1]; 
    ss >> ngp;
  }
  if(!UniqueID()) printf("Number of GPBC: %d\n",ngp);

  int size[] = {4,4,4,4,4};
  bool save_lrg(false);
  char *save_lrg_file;
  bool load_lrg(false);
  char *load_lrg_file;
  
  int nthreads = 1;
  
  int i=2;
  while(i<argc){
    char* cmd = argv[i];  
    if( strncmp(cmd,"-latt",10) == 0){
      if(i>argc-6){
	printf("Did not specify enough arguments for 'latt' (require 5 dimensions)\n"); exit(-1);
      }
      for(int xx=0;xx<5;xx++){
	std::stringstream ss; ss << argv[i+1+xx]; ss >> size[xx];
      }
      i+=6;  
    }else if( strncmp(cmd,"-save_lrg",15) == 0){
      if(i==argc-1){ printf("-save_lrg requires an argument\n"); exit(-1); }
      save_lrg=true;
      save_lrg_file = argv[i+1];
      i+=2;
    }else if( strncmp(cmd,"-load_lrg",15) == 0){
      if(i==argc-1){ printf("-load_lrg requires an argument\n"); exit(-1); }
      load_lrg=true;
      load_lrg_file = argv[i+1];
      i+=2;
    }else if( strncmp(cmd,"-nthread",8) == 0){
      std::stringstream ss; ss << argv[i+1];
      ss >> nthreads;
      if(!UniqueID()){ printf("Setting number of threads to %d\n",nthreads); }
      i+=2;
    }else{
      if(UniqueID()==0) printf("Unrecognised argument: %s\n",cmd);
      exit(-1);
    }
  }

  if(!UniqueID()) printf("Lattice size is %d %d %d %d %d\n",size[0],size[1],size[2],size[3],size[4]);

  DoArg do_arg;
  setupDoArg(do_arg,size,ngp);

  GJP.Initialize(do_arg);
  LRG.Initialize();

  omp_set_num_threads(nthreads);

  int nodes = 1;
  for(int i=0;i<5;i++) nodes *= GJP.Nodes(i);
  
  if(nodes > 1) ERR.General("","","This test code is designed to be run on a local machine!\n");

  if(load_lrg){
    if(UniqueID()==0) printf("Loading RNG state from %s\n",load_lrg_file);
    LRG.Read(load_lrg_file,32);
  }
  if(save_lrg){
    if(UniqueID()==0) printf("Writing RNG state to %s\n",save_lrg_file);
    LRG.Write(save_lrg_file,32);
  }

  int rngs_dir[5];
  for(int i=0;i<5;i++) rngs_dir[i] = GJP.NodeSites(i)/2;

  int vol5d = GJP.VolNodeSites()*GJP.SnodeSites()/32;
  int vol4d = GJP.VolNodeSites()/16;

  int max_threads = omp_get_max_threads();
  std::vector<  std::vector<std::pair<int,int> >   > thread_matches(max_threads);


  //Check 5d generators
  printf("Starting comparison with %d threads\n",max_threads); fflush(stdout);
  Crc32 crc;
  std::map<uint32_t,std::vector<int> > test_matches; //compile by checksum

  assert( sizeof(int) % sizeof(uint32_t) == 0 );

  for(int i=0;i<vol5d*2;i++){
    UGrandomGenerator & rng_i_gen = LRG.UGrandGen(i);
    int const *rng_i = rng_i_gen.getState();
    crc.Reset();
    crc.AddData( (const uint8_t*)rng_i, 55*sizeof(int)/sizeof(uint32_t) );
    test_matches[crc.GetCrc32()].push_back(i);
  }
  std::vector<std::pair<int,int> > all_matches;

  for(std::map<uint32_t,std::vector<int> >::const_iterator it = test_matches.begin(); it != test_matches.end(); it++){
    for(int ii=0;ii<it->second.size();ii++){
      for(int jj=ii+1;jj<it->second.size();jj++){
	int i= it->second[ii];
	int j= it->second[jj];
	UGrandomGenerator & rng_i = LRG.UGrandGen(i);
	UGrandomGenerator & rng_j = LRG.UGrandGen(j);
	if(RNG_equal(rng_i,rng_j)) all_matches.push_back(std::pair<int,int>(i,j));
      }
    }
  }


  //    void AddData(const uint8_t* pData, const uint32_t length)
// #pragma omp parallel for
//   for(int i=0;i<vol5d*2;i++){
//     int me = omp_get_thread_num();
//     if(!me){ printf("%d\n",i); fflush(stdout); }

//     UGrandomGenerator & rng_i = LRG.UGrandGen(i);

//     for(int j=i+1;j<vol5d*2;j++){
//       UGrandomGenerator & rng_j = LRG.UGrandGen(j);

//       if(RNG_equal(rng_i,rng_j)) thread_matches[me].push_back(std::pair<int,int>(i,j));
//     }
//   }

  // printf("Combining\n");


  printf("Sorting\n");
  std::sort(all_matches.begin(), all_matches.end());
  

#if 0
  printf("Printing\n");
  for(int m=0;m<all_matches.size();m++){
    int i = all_matches[m].first;
    int j = all_matches[m].second;

    int xi[5];
    int rem = i;
    for(int ss=0;ss<4;ss++){ xi[ss] = rem % rngs_dir[ss]; rem /= rngs_dir[ss]; }

    //Flavor is stacked inside s loop     rgen_pos = x + hx[0] * (y + hx[1] * (z + hx[2] * (t +hx[3] * (field_idx+ nstacked*s))));
    int fi = rem % 2; rem /= 2;
    xi[4] = rem;

    /////////
    int xj[5];
    rem = j;
    for(int ss=0;ss<4;ss++){ xj[ss] = rem % rngs_dir[ss]; rem /= rngs_dir[ss]; }
      
    int fj = rem % 2; rem /= 2;
    xj[4] = rem;

    /////
    printf("Found match (%d,%d,%d,%d,%d; %d) <-> (%d,%d,%d,%d,%d; %d)\n", xi[0],xi[1],xi[2],xi[3],xi[4],fi,  xj[0],xj[1],xj[2],xj[3],xj[4],fj);
  }

#else

  //Instead of printing all of them (there are a lot!), find compile just the unique deltas
  std::map<std::string,int> unique_deltas;

  for(int m=0;m<all_matches.size();m++){
    int i = all_matches[m].first;
    int j = all_matches[m].second;

    int xi[5];
    int rem = i;
    for(int ss=0;ss<4;ss++){ xi[ss] = rem % rngs_dir[ss]; rem /= rngs_dir[ss]; }

    //Flavor is stacked inside s loop     rgen_pos = x + hx[0] * (y + hx[1] * (z + hx[2] * (t +hx[3] * (field_idx+ nstacked*s))));
    int fi = rem % 2; rem /= 2;
    xi[4] = rem;

    /////////
    int xj[5];
    rem = j;
    for(int ss=0;ss<4;ss++){ xj[ss] = rem % rngs_dir[ss]; rem /= rngs_dir[ss]; }
      
    int fj = rem % 2; rem /= 2;
    xj[4] = rem;


    int xdelta[5];
    for(int d=0;d<5;d++) xdelta[d] = xj[d]-xi[d];
    std::ostringstream os; 
    os << "f1=" << fi << " f2=" << fj << " dx=(";
    for(int d=0;d<5;d++){
      os << xdelta[d]; 
      if(d<4) os << ", ";
    }
    os << ")";

    std::string dstr = os.str();

    if(!unique_deltas.count(dstr)) unique_deltas[dstr] = 1;
    else unique_deltas[dstr]++;
  }

  printf("Unique patterns:\n");
  for(std::map<std::string,int>::const_iterator it = unique_deltas.begin(); it != unique_deltas.end(); it++){
    std::cout << it->first << " count " << it->second << "\n";
  }



#endif





  // int size_5d_lcl = GJP.Gparity() ? 2:1;
  // int size_5d_glb = size_5d_lcl;
  // int sizes_glb[5];
  // int sizes_lcl[5];
  // int node_rng_offset[5];

  // for(int i=0;i<5;i++){
  //   sizes_lcl[i] = GJP.NodeSites(i)/2;
  //   size_5d_lcl *= sizes_lcl[i];

  //   sizes_glb[i] = GJP.NodeSites(i)*GJP.Nodes(i)/2;
  //   size_5d_glb *= sizes_glb[i];

  //   node_rng_offset[i] = GJP.NodeCoor(i)*sizes_lcl[i];
  // }
  // int f_off = size_5d_glb/2; //offset to second flavor




  // if(!UniqueID()){ printf("Generating random field\n"); fflush(stdout); }

  // //Store a global 5D field and just populate this node's segment
  // Float* field = (Float*)malloc(size_5d_glb*sizeof(Float));
  // memset(field,0, size_5d_glb*sizeof(Float)); //zero it!

  // //Generate a 5D random field on each 2^4 hypercube of local volume
  // for(int s=0;s<GJP.SnodeSites();s+=2){
  //   for(int t=0;t<GJP.TnodeSites();t+=2){
  //     for(int z=0;z<GJP.ZnodeSites();z+=2){
  // 	for(int y=0;y<GJP.YnodeSites();y+=2){
  // 	  for(int x=0;x<GJP.XnodeSites();x+=2){  
  // 	    int x_glb_rng[5] = { x/2 + node_rng_offset[0],
  // 				 y/2 + node_rng_offset[1],
  // 				 z/2 + node_rng_offset[2],
  // 				 t/2 + node_rng_offset[3],
  // 				 s/2 + node_rng_offset[4]};
  // 	    int xoff_5d = x_glb_rng[4];
  // 	    for(int i=3;i>=0;i--){
  // 	      xoff_5d *= sizes_glb[i];
  // 	      xoff_5d += x_glb_rng[i];
  // 	    }

  // 	    LRG.AssignGenerator(x,y,z,t,s,0);
  // 	    field[xoff_5d] = LRG.Urand(FIVE_D);
	    
  // 	    if(GJP.Gparity()){
  // 	      LRG.AssignGenerator(x,y,z,t,s,1);
  // 	      field[xoff_5d+f_off] = LRG.Urand(FIVE_D);
  // 	    }
  // 	  }
  // 	}
  //     }
  //   }
  // }

  // if(!UniqueID()){ printf("Running check\n"); fflush(stdout); }
  // std::map<Float,int> rand_map;  //map of random number to RNG global coordinate
  // for(int i=0;i<size_5d_glb;i++){
  //   glb_sum_five(field+i); //now all nodes have the full 5D global field

  //   int rem = i;
  //   int x[5];
  //   for(int xx=0;xx<5;xx++){
  //     x[xx] = rem % sizes_glb[xx]; rem /= sizes_glb[xx];
  //   }
  //   int f = rem;
      
  //   Float val = field[i];
  //   if(!UniqueID()){ printf("%d %g\n",i,val); fflush(stdout); }

  //   std::map<Float,int>::const_iterator it = rand_map.find(val);

  //   if(it != rand_map.end()){
  //     rem = it->second;
  //     int x_other[5];
  //     for(int xx=0;xx<5;xx++){
  //   	x_other[xx] = rem % sizes_glb[xx]; rem /= sizes_glb[xx];
  //     }
  //     int f_other = rem;

  //     if(!UniqueID()){ 
  //   	printf("Found match in random numbers (%g) between global RNG sites (%d,%d,%d,%d,%d; %d) and (%d,%d,%d,%d,%d; %d)\n",
  //   	       field[i],
  //   	       x[0],x[1],x[2],x[3],x[4],f,
  //   	       x_other[0],x_other[1],x_other[2],x_other[3],x_other[4],f_other
  //   	       );
  //     }
  //   }else rand_map[val] = i;

  // }

  // free(field);

  
  End();
  if(!UniqueID()) printf("End of program\n");
  return(0);
}

