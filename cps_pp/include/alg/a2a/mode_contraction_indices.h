#ifndef _MODE_CONTRACTION_INDICES_H
#define _MODE_CONTRACTION_INDICES_H

#include<alg/a2a/mode_mapping.h>

CPS_START_NAMESPACE

//When we contract over modes between two meson fields with different numbers of row/column modes (i.e. because one is packed in time) we need to figure out the mapping between the two
//This class computes and stores the unpacked indices making the contraction transparent

//Deduce which of two dilutions is the most 'packed' (least number of diluted indices)
template<typename DilA, typename DilB>
struct _getMostDilutedType{
  enum { LorR = (int)DilA::UndilutedIndices < (int)DilB::UndilutedIndices ? 0 : 1 };
  typedef typename _selectLR<LorR,DilA,DilB>::Type Type;
};
template<typename DilA, typename DilB>
struct _getLeastDilutedType{
  enum { LorR = (int)DilA::UndilutedIndices >= (int)DilB::UndilutedIndices ? 0 : 1 };
  typedef typename _selectLR<LorR,DilA,DilB>::Type Type;
};

//Predefine meson field
template<typename mf_Float,template <typename> class A2AfieldL,template <typename> class A2AfieldR>
class A2AmesonField;


//Get an index
template<int PackedDepth, int UnpackedDepth>
struct getIndex{
  typedef typename IndexTensor<PackedDepth,UnpackedDepth>::Type TensorType;
  static const std::pair<int,int> & doit(const int &i, const modeIndexSet &packed_coord, const modeIndexSet &unpacked_coord, const TensorType &mode_map){
    int val_packed = IndexConvention<PackedDepth>::get(packed_coord); 
    return getIndex<PackedDepth-1,UnpackedDepth>::doit(i, packed_coord, unpacked_coord,mode_map[val_packed] );
  }
};
template<int UnpackedDepth>
struct getIndex<0,UnpackedDepth>{
  typedef typename IndexTensor<0,UnpackedDepth>::Type TensorType;
  static const std::pair<int,int> & doit(const int &i, const modeIndexSet &packed_coord, const modeIndexSet &unpacked_coord, const TensorType &mode_map){
    int val_unpacked = IndexConvention<UnpackedDepth>::get(unpacked_coord); 
    return getIndex<0,UnpackedDepth-1>::doit(i, packed_coord, unpacked_coord,mode_map[val_unpacked] );
  }
};
template<>
struct getIndex<0,0>{
  typedef typename IndexTensor<0,0>::Type TensorType;
  static const std::pair<int,int> & doit(const int &i, const modeIndexSet &packed_coord, const modeIndexSet &unpacked_coord, const TensorType &mode_map){
    return mode_map[i];
  }
};

//Get number of overlapping indices
template<int PackedDepth, int UnpackedDepth>
struct _getNindices{
  typedef typename IndexTensor<PackedDepth,UnpackedDepth>::Type TensorType;
  static int doit(const modeIndexSet &packed_coord, const modeIndexSet &unpacked_coord, const TensorType &mode_map){
    int val_packed = IndexConvention<PackedDepth>::get(packed_coord); 
    return _getNindices<PackedDepth-1,UnpackedDepth>::doit(packed_coord, unpacked_coord,mode_map[val_packed] );
  }
};
template<int UnpackedDepth>
struct _getNindices<0,UnpackedDepth>{
  typedef typename IndexTensor<0,UnpackedDepth>::Type TensorType;
  static int doit(const modeIndexSet &packed_coord, const modeIndexSet &unpacked_coord, const TensorType &mode_map){
    int val_unpacked = IndexConvention<UnpackedDepth>::get(unpacked_coord); 
    return _getNindices<0,UnpackedDepth-1>::doit(packed_coord, unpacked_coord,mode_map[val_unpacked] );
  }
};
template<>
struct _getNindices<0,0>{
  typedef typename IndexTensor<0,0>::Type TensorType;
  static int doit(const modeIndexSet &packed_coord, const modeIndexSet &unpacked_coord, const TensorType &mode_map){
    return mode_map.size();
  }
};


template<typename dilutionL, typename dilutionR>
class ModeContractionIndices{
  //Determine the least diluted 'packed' type. This will be the one whose indices we loop over
  typedef typename _getLeastDilutedType<dilutionL,dilutionR>::Type PackedType;
  typedef typename _getMostDilutedType<dilutionL,dilutionR>::Type UnpackedType;

  typedef typename ModeMapping<PackedType,UnpackedType>::TensorType TensorType;

  enum { DepthPacked = PackedType::UndilutedIndices };
  enum { DepthUnpacked = UnpackedType::UndilutedIndices };

  TensorType mode_map;

  int packedlr;
public:
  ModeContractionIndices(): packedlr( (int)_getLeastDilutedType<dilutionL,dilutionR>::LorR ){}

  ModeContractionIndices(const A2Aparams &a2a_params): packedlr( (int)_getLeastDilutedType<dilutionL,dilutionR>::LorR ){ 
    compute(a2a_params);
  }

  void compute(const A2Aparams &a2a_params){
    ModeMapping<PackedType,UnpackedType>::compute(mode_map,a2a_params);
  } 
  
  const int & packedLorR() const{ return packedlr; }

  const int &getLeftIndex(const int &i, const modeIndexSet &left_coord, const modeIndexSet &right_coord) const{
    const modeIndexSet &packed_coord = packedlr == 0 ? left_coord : right_coord;
    const modeIndexSet &unpacked_coord = packedlr == 0 ? right_coord : left_coord;

    const std::pair<int,int> &idx_pair = getIndex<DepthPacked,DepthUnpacked>::doit(i, packed_coord, unpacked_coord, mode_map);
    return packedlr == 0 ? idx_pair.first : idx_pair.second;
  }
  const int &getRightIndex(const int &i, const modeIndexSet &left_coord, const modeIndexSet &right_coord) const{
    const modeIndexSet &packed_coord = packedlr == 0 ? left_coord : right_coord;
    const modeIndexSet &unpacked_coord = packedlr == 0 ? right_coord : left_coord;

    const std::pair<int,int> &idx_pair = getIndex<DepthPacked,DepthUnpacked>::doit(i, packed_coord, unpacked_coord, mode_map);
    return packedlr == 1 ? idx_pair.first : idx_pair.second;
  }

  void getBothIndices(int &il, int &ir, const int &i, const modeIndexSet &left_coord, const modeIndexSet &right_coord) const{
    const modeIndexSet &packed_coord = packedlr == 0 ? left_coord : right_coord;
    const modeIndexSet &unpacked_coord = packedlr == 0 ? right_coord : left_coord;

    const std::pair<int,int> &idx_pair = getIndex<DepthPacked,DepthUnpacked>::doit(i, packed_coord, unpacked_coord, mode_map);
    if(packedlr == 0){ il = idx_pair.first; ir = idx_pair.second; }
    else{ ir = idx_pair.first; il = idx_pair.second; } 
  }

  int getNindices(const modeIndexSet &left_coord, const modeIndexSet &right_coord) const{
    const modeIndexSet &packed_coord = packedlr == 0 ? left_coord : right_coord;
    const modeIndexSet &unpacked_coord = packedlr == 0 ? right_coord : left_coord;

    return _getNindices<DepthPacked,DepthUnpacked>::doit(packed_coord, unpacked_coord, mode_map);
  }

};

CPS_END_NAMESPACE

#endif
