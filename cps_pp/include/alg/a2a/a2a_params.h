#ifndef _A2A_PARAMS_H
#define _A2A_PARAMS_H

CPS_START_NAMESPACE

class A2Aparams{
protected:
  A2AArg args;
  int nl; //number of low modes
  int nhits; //number of stochastic hits per site for high modes
  int nflavors; //number of flavors
  int nspincolor; //number of spin*color
  int ntblocks; //number of blocked timeslices
  int ndilute; //number of high mode stochastic dilutions (e.g. if spin/color/timeslice diluting then  ndilute = 12*Lt 
  int nh; //total number of high mode fields = ndilute * nhits
  int nv; //total number of v fields =  nl + ndilute * nhits

public:
  const int &getNl() const{ return nl; }
  const int &getNhits() const{ return nhits; }
  const int &getNflavors() const{ return nflavors; }
  const int &getNspinColor() const{ return nspincolor; }
  const int &getNtBlocks() const{ return ntblocks; }
  const int &getNdilute() const{ return ndilute; }
  const int &getNh() const{ return nh; }
  const int &getNv() const{ return nv; }

  A2Aparams(): nl(0),nhits(0),nflavors(0),nspincolor(0),ntblocks(0),ndilute(0),nh(0),nv(0){}
  A2Aparams(const A2AArg &_args);

  bool paramsEqual(const A2Aparams &r) const{
    return (nl == r.nl  &&  nhits == r.nhits  &&  nflavors == r.nflavors  &&  ntblocks == r.ntblocks  &&  ndilute == r.ndilute  &&  nh == r.nh  &&  nv == r.nv);
  }

  const A2AArg &getArgs() const{ return args; }

  std::string print() const{
    std::ostringstream os; os << "nl=" << nl << " nhits=" << nhits << " nflavors=" << nflavors 
			      << " nspincolor=" << nspincolor << " ntblocks=" << ntblocks << " ndilute=" << ndilute << " nh=" << nh << " nv=" << nv;
    return os.str();
  }


};

CPS_END_NAMESPACE

#endif
