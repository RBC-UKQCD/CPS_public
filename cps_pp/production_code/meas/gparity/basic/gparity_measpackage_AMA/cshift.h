#ifndef _CK_CSHIFT_H
#define _CK_CSHIFT_H

CPS_START_NAMESPACE

//C-shift any 4d field with canonical ordering in the time direction
//(should be trivial to generalize!)
void TbarrelShift4D(Float* data, const int site_size, const int tshift){
  if(tshift == 0) return;

  const int plane_size = GJP.VolNodeSites()/GJP.TnodeSites();
  const int bulk_size = GJP.VolNodeSites()-plane_size;

  const int nf = GJP.Gparity() ? 2:1;
  Float *tplane_send = (Float*)malloc(site_size * plane_size * nf * sizeof(Float));
  Float *tplane_recv = (Float*)malloc(site_size * plane_size * nf * sizeof(Float));
  Float *bulk = (Float*)malloc(site_size * bulk_size * nf * sizeof(Float));
 
  const int tsend = tshift > 0 ? GJP.TnodeSites()-1 : 0;
  const int trecv = tshift > 0 ? 0 : GJP.TnodeSites()-1;

  const int tbulk_start = tshift > 0 ? 0 : 1;
  const int tbulk_lessthan = tshift > 0 ? GJP.TnodeSites()-1 : GJP.TnodeSites();
  
  const int shift_one = tshift > 0 ? 1 : -1;

  const int gp_foff_bulk = site_size*bulk_size;
  const int gp_foff_plane = site_size*plane_size;
  const int gp_foff_orig = site_size*GJP.VolNodeSites();

  const int nshift = tshift > 0 ? tshift : -tshift;

  for(int shift = 0; shift < nshift; shift++){
    //Pull out the time plane we are sending and the rest of the *bulk*
#pragma omp parallel for
    for(int x3d=0;x3d<plane_size;x3d++){
      for(int f=0;f<nf;f++){
	int tb = 0;
	for(int t=tbulk_start; t<tbulk_lessthan; t++){
	  Float *orig_off = data + f*gp_foff_orig + site_size*(x3d + plane_size * t);
	  Float *bulk_off = bulk + f*gp_foff_bulk + site_size*(x3d + plane_size * tb);
	  memcpy(bulk_off,orig_off,site_size*sizeof(Float));
	  ++tb;
	}
	Float *orig_off = data + f*gp_foff_orig + site_size*(x3d + plane_size * tsend);
	Float *plane_off = tplane_send + f*gp_foff_plane + site_size*x3d;
	memcpy(plane_off,orig_off,site_size*sizeof(Float));
      }
    }
    if(tshift > 0)
      getMinusData(tplane_recv,tplane_send,nf*site_size*plane_size,3); //send in +t direction
    else
      getPlusData(tplane_recv,tplane_send,nf*site_size*plane_size,3); //send in +t direction
        
#pragma omp parallel for
    for(int x3d=0;x3d<plane_size;x3d++){
      for(int f=0;f<nf;f++){
	int tb = 0;
	for(int t=tbulk_start + shift_one; t<tbulk_lessthan + shift_one; t++){ //put the bulk back but shifted one in time
	  Float *bulk_off = bulk + f*gp_foff_bulk + site_size*(x3d + plane_size * tb);
	  Float *orig_off = data + f*gp_foff_orig + site_size*(x3d + plane_size * t);
	  memcpy(orig_off,bulk_off,site_size*sizeof(Float));
	  ++tb;
	}
	Float *plane_off = tplane_recv + f*gp_foff_plane + site_size*x3d;
	Float *orig_off = data + f*gp_foff_orig + site_size*(x3d + plane_size * trecv);
	memcpy(orig_off,plane_off,site_size*sizeof(Float));
      }
    }
  }

  free(tplane_send);
  free(tplane_recv);
  free(bulk);
}





void Tshift4D(Float* data, const int site_size, const int tshift){
  TbarrelShift4D(data, site_size, tshift);
}


CPS_END_NAMESPACE

#endif
