#ifndef _WALLSINK_PROP_H
#define _WALLSINK_PROP_H

CPS_START_NAMESPACE


template<typename MatrixType>
class siteMatrixCompute{};

template<>
class siteMatrixCompute<SpinColorFlavorMatrix>{
  PropWrapper prop;
  bool setup;
public:
  siteMatrixCompute():setup(false){}

  bool isSetup() const{ return setup; }

  void setProp(const PropWrapper &_prop){
    prop = _prop; setup = true;
  }
  void siteMatrix(SpinColorFlavorMatrix &into, const int site){
    prop.siteMatrix(into,site);
  }
  void multGFmat(SpinColorFlavorMatrix &what, const int site, Lattice &lat){
    Matrix const* gfmat_f0 = lat.FixGaugeMatrix(site,0);
    Matrix const* gfmat_f1 = lat.FixGaugeMatrix(site,1);
    if(gfmat_f0 == NULL || gfmat_f1 == NULL) ERR.General("siteMatrixCompute<SpinColorFlavorMatrix>","multGFmat","No gauge fixing matrix for site %d",site);
    what(0,0).LeftTimesEqual(*gfmat_f0);
    what(0,1).LeftTimesEqual(*gfmat_f0);
    what(1,0).LeftTimesEqual(*gfmat_f1);
    what(1,1).LeftTimesEqual(*gfmat_f1);
  }

  void latticeSum(SpinColorFlavorMatrix &what){
    for(int i=0;i<2;i++)
      for(int j=0;j<2;j++){
	WilsonMatrix &wm = what(i,j);
	Float* w = (Float*)wm.ptr();
	static const int size = 2*12*12;
	slice_sum(w, size, 99); //99 is a *magic* number (we are abusing slice_sum here)
      }
  }
};

template<>
class siteMatrixCompute<WilsonMatrix>{
  const QPropW *prop;
  bool setup;
public:
  siteMatrixCompute():setup(false){}

  bool isSetup() const{ return setup; }

  void setProp(const QPropW &_prop){
    prop = &_prop; setup = true;
  }
  void setProp(const PropWrapper &_prop){
    prop = _prop.getPtr(0); setup = true;
  }

  void siteMatrix(WilsonMatrix &into, const int site){
    into = prop->SiteMatrix(site);
  }
  void multGFmat(WilsonMatrix &what, const int site, Lattice &lat){
    Matrix const* gfmat = lat.FixGaugeMatrix(site);
    if(gfmat == NULL) ERR.General("siteMatrixCompute<WilsonMatrix>","multGFmat","No gauge fixing matrix for site %d",site);
    what.LeftTimesEqual(*gfmat);
  }

  void latticeSum(WilsonMatrix &what){
    Float* w = (Float*)what.ptr();
    static const int size = 2*12*12;
    slice_sum(w, size, 99); //99 is a *magic* number (we are abusing slice_sum here)
  }  
};



template<typename MatrixType>
class WallSinkProp: public siteMatrixCompute<MatrixType>{
  std::vector<MatrixType> result;

  void global_coord(const int site, int into_vec[4]){    
    int rem = site;
    for(int i=0;i<4;i++){
      into_vec[i] = rem % GJP.NodeSites(i) + GJP.NodeCoor(i)*GJP.NodeSites(i);
      rem /= GJP.NodeSites(i);
    }
  }
  bool gauge_fix_sink;

public:
  WallSinkProp(const bool _gauge_fix_sink = true): gauge_fix_sink(_gauge_fix_sink){ }
  
  //Can provide an optional sink momentum in lattice units, for which the phase exp(-ip.x) is then applied
  void compute(Lattice &lat, const double *p = NULL){
    if(!this->isSetup()) ERR.General("WallSinkProp","compute","Class has not been set up\n");
    const int global_T = GJP.TnodeSites()*GJP.Tnodes();
    const int local_T = GJP.TnodeSites();
    const int local_toff = GJP.TnodeCoor()*local_T;
    
    const int nthread = omp_get_max_threads();
    
    std::vector<std::vector<MatrixType> > thread_mats(global_T); //[t][thread]
    for(int t=0;t<global_T;t++){
      thread_mats[t].resize(nthread);
      for(int thr=0;thr<nthread;thr++)
	thread_mats[t][thr] = 0.0;
    }

#pragma omp parallel for
    for(int x=0;x<GJP.VolNodeSites();x++){
      int x_glb[4]; global_coord(x,x_glb);
      MatrixType mat; this->siteMatrix(mat, x);
      
      //Gauge fix sink
      if(gauge_fix_sink){
	if(!lat.FixGaugeKind()) ERR.General("WallSinkProp","compute","Lattice is not gauge fixed!\n");
	this->multGFmat(mat,x,lat);
      }

      if(p!=NULL){
	Float pdotx = 0.0;
	for(int i=0;i<3;i++) pdotx += p[i]*x_glb[i];
	mat *= Complex(cos(pdotx),-sin(pdotx));
      }

      thread_mats[x_glb[3]][omp_get_thread_num()] += mat;
    }

    //Thread sum
    result.resize(global_T);
#pragma omp parallel for
    for(int t=0;t<global_T;t++){
      result[t] = thread_mats[t][0];
      for(int thr=1;thr<nthread;thr++){
	result[t] += thread_mats[t][thr];
      }
    }
    thread_mats.clear();
    
    //Lattice sum
    for(int t=0;t<global_T;t++)
      this->latticeSum(result[t]);
  }

  void compute(Lattice &lat, const ThreeMomentum &p){
    Float pp[3]; p.latticeUnits(pp);
    compute(lat,pp);
  }


  const MatrixType & operator()(const int t_glb) const{
    return result[t_glb];
  }
};



CPS_END_NAMESPACE

#endif
