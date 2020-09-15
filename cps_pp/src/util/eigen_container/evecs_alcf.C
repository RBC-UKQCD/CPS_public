#include <config.h>

#include<util/gjp.h>
#include<util/eigen_container.h>

//USING_NAMESPACE_CPS
using namespace std;


namespace cps {
//uint32_t crc32_fast(const void* data, size_t length, uint32_t previousCrc32);


int get_alcf_index( int x, int y, int z, int t, int* latt, int s, int ivec, int co ) {
  //
  
  int ls = GJP.SnodeSites();
  int vol_4d_oo = latt[0]*latt[1]*latt[2]*latt[3] / 2;
  int vol_5d = vol_4d_oo * ls;
  
  int NtHalf = latt[3] / 2;
  int simd_coor = t / NtHalf;
  int regu_coor = (x + latt[0] * (y + latt[1] * ( z + latt[2] * (t % NtHalf) ) )) / 2;
  int regu_vol  = vol_4d_oo / 2;
  return ivec * vol_5d * 24       
    + regu_coor * ls * 48
    + s * 48
    + co * 4  
    + simd_coor * 2;
}


int get_index( int x, int y, int z, int t, int* latt, int s, int ivec, int co ) {
  //
  
  int ls = GJP.SnodeSites();
  int vol_4d_oo = latt[0]*latt[1]*latt[2]*latt[3] / 2;
  int vol_5d = vol_4d_oo * ls;
  
  int regu_coor = (x + latt[0] * (y + latt[1] * ( z + latt[2] * t ) )) / 2;
  int regu_vol  = vol_4d_oo / 2;
  return ivec * vol_5d * 24       
    + s * vol_4d_oo * 24
    + regu_coor * 24
    + co * 2;
}


void alcf_evecs_save(char* dest,EigenCache* ec, int nkeep) {

  // format string
  char* dest_sep = strrchr(dest,'/');
  *dest_sep = '\0';
  char* meta_tag = dest_sep + 1;

  FILE* f;
  char fn[2048];
  int i;

  if (!UniqueID()) {
    printf("STATUS: Writing %d eigenvectors to %s in ALCF format\n",
           nkeep,dest);
  }

  // 1) write eigen-values.txt
  if (!UniqueID()) {
    sprintf(fn,"%s/eigen-values.txt",dest);
    f=fopen(fn,"wt");
    if (f) {
      fprintf(f,"%d\n",nkeep);
      for (i=0;i<nkeep;i++)
        fprintf(f,"%.20E\n",ec->eval_address()[i]);
      fclose(f);
    }
  }


  // 2) write eigenvectors in 32 directories
  int id = UniqueID();
  int slots_total = GJP.Nodes(0) * GJP.Nodes(1) * GJP.Nodes(2) * GJP.Nodes(3);
  int dir_id = id / (slots_total / 32);
  int nwrite_at_once = 32;
  int ngroups = slots_total / nwrite_at_once;
  int igroup;
  
  uint32_t crc32_all = 0x0;
  for (igroup=0;igroup<ngroups;igroup++) {
    Float test = 0.0;
    if (id % ngroups == igroup) {
      // write now
      sprintf(fn,"%s/%.2d",dest,dir_id);
      mkdir(fn,0777);
      sprintf(fn,"%s/%.2d/%.10d",dest,dir_id,id);
      f=fopen(fn,"wt");
      if (f) {
        int size_float = GJP.VolNodeSites()*GJP.SnodeSites()/2*12*2;
        float* ev = new float[size_float];
        memset(ev,0,sizeof(float)*size_float);
        for (i=0;i<nkeep;i++) {
//          float* ev_in = &((float*)ec->evec_address())[size_float*i];
          float* ev_in = (float*)ec->vec_ptr(i);
          
          // Keep in mind: without re-ordering ev<>ev_in the norm is OK when reading in, with it is not!
          {
            int x,y,z,t,s,co;
            int latt[] = { GJP.NodeSites(0), GJP.NodeSites(1), GJP.NodeSites(2), GJP.NodeSites(3) };
            for (x=0;x<GJP.NodeSites(0);x++)
              for (y=0;y<GJP.NodeSites(1);y++)
                for (z=0;z<GJP.NodeSites(2);z++)
                  for (t=0;t<GJP.NodeSites(3);t++) 
                    if (( x + y + z + t ) % 2 == 1)
                      for (s=0;s<GJP.SnodeSites();s++)
                        for (co=0;co<12;co++) {
                          int idx = get_index( x, y, z, t, latt, s, 0, co );
                          int idx_alcf = get_alcf_index( x, y, z, t, latt, s, 0, co );
                          
                          ev[idx_alcf] = ev_in[idx];
                          ev[idx_alcf + 1] = ev_in[idx + 1];
                        }
          }

//          crc32 = crc32_fast(ev, sizeof(float)*size_float, crc32);
          crc32_all = crc32(crc32_all, (const Bytef*) ev, sizeof(float)*size_float);
	  {
	    int nr = 0;
	    while (fwrite(ev,sizeof(float)*size_float,1,f) != 1) {
	      printf("WRITE-WARNING: node %d failed to write, retry %d in 5 seconds\n",UniqueID(),nr++);
	      sleep(5);
	    }
	  }
        }
        delete[] ev;

        fclose(f);
      }
      glb_sum(&test);
    }
  }

  // 3) write metadata.txt
  if (!UniqueID()) {
    sprintf(fn,"%s/metadata.txt",dest);
    f=fopen(fn,"wt");
    if (f) {
      fprintf(f,"tag: %s\n",meta_tag);
      fprintf(f,"nodes: %d %d %d %d\n",
              GJP.Nodes(0),GJP.Nodes(1),
              GJP.Nodes(2),GJP.Nodes(3));
      fprintf(f,"lattice: %d %d %d %d\n",
              GJP.Sites(0),GJP.Sites(1),
              GJP.Sites(2),GJP.Sites(3));
      fclose(f);
    }
  }

  // 4) write checksums.txt
  sprintf(fn,"%s/checksums.txt",dest);
  if (!UniqueID()) {
    f=fopen(fn,"wt");
    if (f)
      fprintf(f,"XXXXXXXX\n\n");
  } else {
    f = 0;
  }

  for (i=0;i<slots_total;i++) {
    uint32_t this_crc32 = 0x0;
    if (i == id)
      this_crc32 += crc32_all;
    uint32_t rcv_crc32;
    MPI_Allreduce(&this_crc32, &rcv_crc32, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    if (f)
      fprintf(f,"%X\n",rcv_crc32);
  }
  
  if (f)
    fclose(f);
}

};
