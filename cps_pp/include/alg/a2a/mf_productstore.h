#ifndef  _MF_PRODUCT_STORE_H
#define  _MF_PRODUCT_STORE_H

CPS_START_NAMESPACE

//Try to avoid recomputing products of meson fields by re-using wherever possible
template<typename mf_Policies>
class MesonFieldProductStore{
  typedef A2AmesonField<mf_Policies,A2AvectorWfftw,A2AvectorVfftw> MfType;
  typedef A2AmesonField<mf_Policies,A2AvectorWfftw,A2AvectorVfftw> const* MfPtr;
  typedef std::pair<MfPtr,MfPtr> KeyType;
  typedef std::map<KeyType,MfType> MapType;

  MapType products; //compute the product for each time separation independently
  const MfType & compute(MfPtr a, MfPtr b, const bool node_local){
    KeyType key(a,b);
    int Lt = GJP.TnodeSites()*GJP.Tnodes();
    MfType &into = products[key];
    mult(into,*a,*b,node_local);
    return into;
  }
public:
  //Product 'compute' is multi-node so please don't use this inside a node-specific piece of code
  const MfType & getProduct(const MfType &a, const MfType &b, const bool node_local = false){
    KeyType key(&a,&b);
    typename MapType::const_iterator pp = products.find(key);
    if(pp != products.end()) return pp->second;
    else return compute(&a,&b,node_local);
  }
  //Const version fails if product has not been recomputed
  const MfType & getPrecomputedProduct(const MfType &a, const MfType &b) const{
    KeyType key(&a,&b);
    typename MapType::const_iterator pp = products.find(key);
    if(pp != products.end()){
      return pp->second;
    }
    ERR.General("MesonFieldProductStore","getProduct (const version)","Product not pre-computed\n");
  } 
  
  int size() const{ return products.size(); }

  //Return total size in bytes of all stored meson fields
  size_t byte_size() const{
    size_t size = 0;
    for(typename MapType::const_iterator it = products.begin(); it != products.end(); it++){
      size += it->second.byte_size();
    }
    return size;
  }

};

CPS_END_NAMESPACE

#endif
