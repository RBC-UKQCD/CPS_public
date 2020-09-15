#ifndef _MODE_MAPPING_H
#define _MODE_MAPPING_H

CPS_START_NAMESPACE

//As we dilute we unpack in the following order: spincolor , flavor, time. We assign these indices  spincolor=3, flavor=2, time=1
template<int TypeIndex>
struct IndexConvention{};

template<>
struct IndexConvention<3>{
  static int getNidx(){ return 12; }
  static std::string getStr(){ return "spin-color"; }
  static void set(modeIndexSet &into, const int &val){ into.spin_color = val; } 
  static const int & get(const modeIndexSet &from){ return from.spin_color; } 
};
template<>
struct IndexConvention<2>{
  static int getNidx(){ return GJP.Gparity() ? 2:1; }
  static std::string getStr(){ return "flavor"; }  
  static void set(modeIndexSet &into, const int &val){ into.flavor = val; } 
  static const int & get(const modeIndexSet &from){ return from.flavor; } 
};
template<>
struct IndexConvention<1>{
  static int getNidx(){ return GJP.Tnodes()*GJP.TnodeSites(); }
  static std::string getStr(){ return "time"; }
  static void set(modeIndexSet &into, const int &val){ into.time = val; } 
  static const int & get(const modeIndexSet &from){ return from.time; } 
};


//To store an index mapping we need a large number of matrices
typedef std::vector<std::pair<int,int> > ModeMapType;

//Unmapping a single index
template<int Depth>
struct IndexVector{
  typedef typename IndexVector<Depth-1>::Type SubType;
  typedef std::vector<SubType> Type;
};
template<>
struct IndexVector<0>{
  typedef std::pair< std::vector<int>, std::vector<bool> > Type;
};

//We want a big set of nested vectors where the first Ldepth vectors are associated with the unmappings for the left index, and the remaining Rdepth with the right index
template<int Ldepth, int Rdepth>
struct IndexTensor{
  typedef typename IndexTensor<Ldepth-1,Rdepth>::Type SubType;
  typedef std::vector<SubType> Type;
};
template<int Rdepth>
struct IndexTensor<0,Rdepth>{
  typedef typename IndexTensor<0,Rdepth-1>::Type SubType;
  typedef std::vector<SubType> Type;
};
template<>
struct IndexTensor<0,0>{
  typedef ModeMapType Type; //mode * mode
};


template<typename PackedType, typename UnpackedType>
class ModeMapping{
public:
  enum { DepthPacked = PackedType::UndilutedIndices };
  enum { DepthUnpacked = UnpackedType::UndilutedIndices };

  typedef typename IndexVector<DepthPacked>::Type VectorTypePacked;
  typedef typename IndexVector<DepthUnpacked>::Type VectorTypeUnpacked;

  typedef typename IndexTensor<DepthPacked, DepthUnpacked>::Type TensorType;

  static void compute(TensorType &idx_map, const A2Aparams &p);
};



#include<alg/a2a/mode_mapping_impl.h>

CPS_END_NAMESPACE


#endif
