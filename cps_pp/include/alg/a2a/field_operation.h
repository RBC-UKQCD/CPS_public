#ifndef _FIELD_OPERATION_H
#define _FIELD_OPERATION_H

CPS_START_NAMESPACE

//Perform some operation on a field
template <typename FieldType>
class fieldOperation{
public:
  virtual void operator()(const FieldType &in, FieldType &out) = 0;
};

template <typename FieldType>
class twist: public fieldOperation<FieldType>{
  int p[3];
public:
  twist(const int _p[3]){
    memcpy(p,_p,3*sizeof(int));
  }

  void operator()(const FieldType &in, FieldType &out){
    //Gauge fix and apply phase in parallel (i.e. don't parallelize over modes)
    out = in;
#ifndef MEMTEST_MODE
    out.applyPhase(p,true);
#endif
  }
};


template <typename FieldType>
class gaugeFixAndTwist: public fieldOperation<FieldType>{
  int p[3];
  Lattice *lat;
  
public:
  gaugeFixAndTwist(const int _p[3], Lattice &_lat): lat(&_lat){ memcpy(p,_p,3*sizeof(int)); }

  void operator()(const FieldType &in, FieldType &out){
    //Gauge fix and apply phase in parallel (i.e. don't parallelize over modes)
    out = in;
#ifndef MEMTEST_MODE
    out.gaugeFix(*lat,true);
    out.applyPhase(p,true);
#endif
  }
};


//Apply the inverse of the above 
template <typename FieldType>
class reverseGaugeFixAndTwist: public fieldOperation<FieldType>{
  int mp[3];
  Lattice *lat;
  
public:
  reverseGaugeFixAndTwist(const int p[3], Lattice &_lat): lat(&_lat){ for(int i=0;i<3;i++) mp[i] = -p[i]; }

  void operator()(const FieldType &in, FieldType &out){
    out = in;
#ifndef MEMTEST_MODE
    out.applyPhase(mp,true); //apply - the phase
    out.gaugeFix(*lat,true,true); //last bool is optional dagger
#endif
  }
};




CPS_END_NAMESPACE

#endif
