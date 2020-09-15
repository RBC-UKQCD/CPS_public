#ifndef _MODE_MAPPING_IMPL_H
#define _MODE_MAPPING_IMPL_H

/////////Work out the unmapping for a single type
template<int Depth, typename DilutionType>
struct computeModeUnmapping{
  typedef typename IndexVector<Depth>::Type VectorType;
  static void doit(VectorType &v,modeIndexSet &coord, const DilutionType &dil){
    int nidx = IndexConvention<Depth>::getNidx();
    v.resize(nidx);
    for(int i=0;i<nidx;i++){
      IndexConvention<Depth>::set(coord,i);
      computeModeUnmapping<Depth-1,DilutionType>::doit(v[i],coord, dil);
    }
  }
};

template<typename DilutionType>
struct computeModeUnmapping<0,DilutionType>{
  typedef typename IndexVector<0>::Type VectorType;

  static void doit(VectorType &v,modeIndexSet &coord, const DilutionType &dil){
    std::vector<int> &std_to_dil = v.first;
    std::vector<bool> &non_zeroes = v.second;
    dil.getIndexMapping(std_to_dil, non_zeroes, coord);
  }
};

//Compute the mode map between types recursively

//First move through packed nested coords   //spincolor=3, flavor=2, time=1
template<int DepthPacked, int DepthUnpacked, typename PackedType, typename UnpackedType>
struct computeModeMap{
  typedef typename IndexVector<DepthPacked>::Type VectorTypePacked;
  typedef typename IndexVector<DepthUnpacked>::Type VectorTypeUnpacked;

  typedef typename IndexTensor<DepthPacked, DepthUnpacked>::Type TensorType;
  
  static void doit(TensorType &v,const VectorTypePacked &packed_unmap, const VectorTypeUnpacked &unpacked_unmap){
    int nidx = IndexConvention<DepthPacked>::getNidx();
    v.resize(nidx);
    for(int i=0;i<nidx;i++)
      computeModeMap<DepthPacked-1,DepthUnpacked, PackedType,UnpackedType>::doit(v[i],packed_unmap[i],unpacked_unmap);
  }
};
template<int DepthUnpacked, typename PackedType, typename UnpackedType>
struct computeModeMap<0,DepthUnpacked,PackedType,UnpackedType>{ //gotten down to the base modes for packed, start on unpacked
  typedef typename IndexVector<0>::Type VectorTypePacked;
  typedef typename IndexVector<DepthUnpacked>::Type VectorTypeUnpacked;

  typedef typename IndexTensor<0, DepthUnpacked>::Type TensorType;
  
  static void doit(TensorType &v,const VectorTypePacked &packed_unmap, const VectorTypeUnpacked &unpacked_unmap){
    int nidx = IndexConvention<DepthUnpacked>::getNidx();
    v.resize(nidx);
    for(int i=0;i<nidx;i++)
      computeModeMap<0,DepthUnpacked-1, PackedType,UnpackedType>::doit(v[i],packed_unmap,unpacked_unmap[i]);
  }
};
template<typename PackedType, typename UnpackedType>
struct computeModeMap<0,0,PackedType,UnpackedType>{ //gotten down to the base modes for packed, start on unpacked
  typedef typename IndexVector<0>::Type VectorTypePacked;
  typedef typename IndexVector<0>::Type VectorTypeUnpacked;

  typedef typename IndexTensor<0, 0>::Type TensorType; //mode * mode
  
  static void doit(TensorType &v,const VectorTypePacked &packed_unmap, const VectorTypeUnpacked &unpacked_unmap){
    //Fully unpack both and find the overlap between the sets of non-zero indices
    const std::vector<bool> &non_zeroes_packed = packed_unmap.second;
    const std::vector<bool> &non_zeroes_unpacked = unpacked_unmap.second;

    const std::vector<int> &std_to_packed = packed_unmap.first;
    const std::vector<int> &std_to_unpacked = unpacked_unmap.first;

    std::vector<bool> overlap;
    compute_overlap(overlap, non_zeroes_packed, non_zeroes_unpacked);

    int n_std = overlap.size();
    assert(std_to_packed.size() == n_std);
    assert(std_to_unpacked.size() == n_std);

    v.resize(0); v.reserve(n_std);

    for(int i=0;i<n_std;i++)
      if(overlap[i]){
	std::pair<int,int> idx_pair(std_to_packed[i], std_to_unpacked[i]);
	v.push_back(idx_pair);
      }
  }
};



template<typename PackedType, typename UnpackedType>
void ModeMapping<PackedType,UnpackedType>::compute(TensorType &idx_map, const A2Aparams &p){
  modeIndexSet tmp;
  const PackedType &packed = static_cast<const PackedType &>(p);
  const UnpackedType &unpacked = static_cast<const UnpackedType &>(p);
    
  //Completely unpack the 'packed' and 'unpacked' (or *less* packed anyway - it doesn't have to be a fully unpacked index) type
  VectorTypePacked packed_unmap;
  computeModeUnmapping<DepthPacked,PackedType>::doit(packed_unmap,tmp,packed);

  VectorTypeUnpacked unpacked_unmap;
  computeModeUnmapping<DepthUnpacked,UnpackedType>::doit(unpacked_unmap,tmp,unpacked);
  
  computeModeMap<DepthPacked,DepthUnpacked, PackedType,UnpackedType>::doit(idx_map,packed_unmap,unpacked_unmap);
}




#endif
