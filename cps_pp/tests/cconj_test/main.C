#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <util/qcdio.h>
#ifdef PARALLEL
#include <comms/sysfunc_cps.h>
#endif
#include <comms/scu.h>
#include <comms/glb.h>

#include <util/lattice.h>
#include <util/time_cps.h>
#include <alg/do_arg.h>
#include <alg/no_arg.h>
#include <alg/common_arg.h>
#include <alg/hmd_arg.h>
#include <alg/alg_plaq.h>
#include <alg/alg_rnd_gauge.h>
#include <alg/threept_arg.h>
#include <alg/threept_prop_arg.h>
#include <alg/alg_threept.h>
#include <util/smalloc.h>
#if(0==1)
 #include <ReadLattice.h>
 #include <WriteLattice.h>
#endif

#include <util/ReadLatticePar.h>
#include <util/WriteLatticePar.h>

#include <util/command_line.h>
#include <sstream>

#include<unistd.h>
#include<config.h>

#include <alg/qpropw.h>
#include <alg/qpropw_arg.h>

#include <alg/alg_fix_gauge.h>
#include <alg/fix_gauge_arg.h>

#include <util/data_shift.h>

#include <alg/prop_attribute_arg.h>
#include <alg/gparity_contract_arg.h>
#include <alg/propmanager.h>
#include <alg/alg_gparitycontract.h>

using namespace std;
USING_NAMESPACE_CPS

void unit_matrix(WilsonMatrix &m){
  for(int i=0;i<4;i++){
    for(int j=0;j<4;j++){
      for(int a=0;a<3;a++){
  	for(int b=0;b<3;b++){
	  if(i==j && a == b) m(i,a,j,b) = Complex(1.0,0.0);
  	  else m(i,a,j,b) = Complex(0.0,0.0);
  	}
      }
    }
  }
}
void print_spin(const WilsonMatrix &m){
  printf("\n");
  for(int i=0;i<4;i++){
    for(int j=0;j<4;j++){
      double re = m(i,0,j,0).real();
      double im = m(i,0,j,0).imag();
      if(re == 0.0 && im == 1.0) printf("i ");
      else if(re == 0.0 && im == -1.0) printf("-i ");
      else if(re == 1.0 && im == 0.0) printf("1 ");
      else if(re == -1.0 && im == 0.0) printf("-1 ");
      else if(re == 0.0 && im == 0.0) printf("0 ");
      else { printf("print can only do +-1 or +-i, got (%f,%f)\n",re,im); exit(-1); }
    }
    printf("\n");
  }
  printf("\n");
}

int main(int argc,char *argv[])
{
  Start(&argc,&argv);
  CommandLine::is(argc,argv);
  
  WilsonMatrix m;
  unit_matrix(m);

  m.ccr(1); 
  printf("C = \n");
  print_spin(m);

  WilsonMatrix g2g4;
  unit_matrix(g2g4);
  
  //CPS gamma matrices have opposite sign than they should in gamma^1 and gamma^3.
  g2g4.gr(1).gr(3);
  g2g4*= Complex(-1.0,0.0);

  print_spin(g2g4);

  //check C* = C   C^-1 = -C
  WilsonMatrix Cinv; unit_matrix(Cinv); Cinv.ccr(-1);
  WilsonMatrix Cstar; unit_matrix(Cstar); Cstar.ccr(1); Cstar.cconj();

  printf("C* = \n");
  print_spin(Cstar);
  
  printf("C^-1 = \n");
  print_spin(Cinv);

  //check C^2 = -1
  WilsonMatrix C2; unit_matrix(C2);
  C2.ccr(1).ccr(1);
  printf("C^2 = \n");
  print_spin(C2);

  //check ccl == ccr
  WilsonMatrix Cl; unit_matrix(Cl);
  Cl.ccl(1);
  printf("C(left)I = \n");
  print_spin(Cl);

  //check C^2 = -1 alternate
  WilsonMatrix C2A; unit_matrix(C2A);
  C2A.ccl(1).ccr(1);
  printf("C^2 alternate = \n");
  print_spin(C2A);


  //test C gamma^mu C^-1 = [gamma^mu]^T
  for(int mu=0;mu<4;mu++){
    WilsonMatrix gamma; unit_matrix(gamma); gamma.gr(mu);
    printf("gamma^%d = ",mu+1);
    print_spin(gamma);

    WilsonMatrix CgC = gamma;
    CgC.ccl(-1).ccr(-1); //CPS code has minus sign wrong in ccl
    printf("C gamma^%d C^-1 = ",mu+1);

    print_spin(CgC);

    WilsonMatrix gT = gamma;
    gT.transpose();

    printf("(gamma^%d)^T =",mu+1);
    print_spin(gT);

    WilsonMatrix mCgstarC = gamma*Complex(-1.0,0.0); mCgstarC.cconj();
    mCgstarC.ccl(-1).ccr(1);
    
    printf("-C (gamma^%d )*C = ",mu+1);
    print_spin(mCgstarC);

    WilsonMatrix t1 = gamma; t1.cconj();
    t1.ccl(-1).ccr(-1);
    
    printf("C (gamma^%d )*C^-1 = ",mu+1);
    print_spin(t1);



    printf("\n\n");
  }
  
      
  return 0;
}




  // prop_args.props.props_len = 1;
  // prop_args.props.props_val = new PropagatorArg[1];

  // prop_args.props.props_val[0].attributes.attributes_len = 2;
  // prop_args.props.props_val[0].attributes.attributes_val = new AttributeSet[2];
  // prop_args.props.props_val[0].attributes.attributes_val[0].type = POINT_SOURCE_ATTR;
  // prop_args.props.props_val[0].attributes.attributes_val[1].type = MOMENTUM_ATTR;

  // prop_args.Encode("parg.vml","prop_arg");



  
  // Matrix *gauge = lattice.GaugeField();
  // int volnodesites = GJP.VolNodeSites();
  
  //int pos[4];

  // for(int t=0;t<GJP.NodeSites(3);t++){
  //   pos[3] = t;
  //   for(int z=0;z<GJP.NodeSites(2);z++){
  //     pos[2] = z;
  //     for(int y=0;y<GJP.NodeSites(1);y++){
  // 	pos[1] = y;
  // 	for(int x=0;x<GJP.NodeSites(0);x++){
  // 	  pos[0] = x;

  // 	  for(int mu=0;mu<4;mu++){
  // 	    Matrix *U = const_cast<Matrix *>(lattice.GetLink(pos, mu));
  // 	    Matrix *Uconj = const_cast<Matrix *>(lattice.GetLink(pos, mu,1));

  // 	    for(int i=0;i<3;i++){
  // 	    	for(int j=0;j<3;j++){
  // 	    	  const Complex &c = (*U)(i,j);
  // 		  const Complex &cstar = (*Uconj)(i,j);
  // 		  if(c.real() != cstar.real() || c.imag() != -cstar.imag() ){
  // 		    ERR.General("main","()","Non-conj problem at %d %d %d %d; %d :  %d %d\n",x,y,z,t,mu,i,j);
  // 		  }
  // 	    	}
  // 	    }
  // 	  }
	    

  // 	}
  //     }
  //   }
  // }


  // printf("U:\n");
  // for(int i=0;i<4*volnodesites;i++){
  //   printf("i=%d ",i);
  //   for(int j=0;j<18;j++){
  //     if(gauge[i].elem(j)>=0)
  // 	printf("%d: %.3f ",j,gauge[i].elem(j));
  //     else
  // 	printf("%d:%.3f ",j,gauge[i].elem(j));
  //   }
  //   printf("\n");
  // }
  // if(GJP.Gparity()){
  //   printf("U*:\n");
  //   for(int i=4*volnodesites;i<8*volnodesites;i++){
  //     printf("i=%d ",i-4*volnodesites);
  //     for(int j=0;j<18;j++){
	
  // 	if(j%2==0 && gauge[i].elem(j) != gauge[i-4*volnodesites].elem(j)){
  // 	  printf("Error re %d,%d : %.3f %.3f\n",i,j,gauge[i].elem(j),gauge[i-4*volnodesites].elem(j));
  // 	}else if(j%2==1 && gauge[i].elem(j) != gauge[i-4*volnodesites].elem(j)){
  // 	  printf("Error im %d,%d : %.3f %.3f\n",i,j,gauge[i].elem(j),gauge[i-4*volnodesites].elem(j));
  // 	}

  // 	if(gauge[i].elem(j)>=0)
  // 	  printf("%d: %.3f ",j,gauge[i].elem(j));
  // 	else
  // 	  printf("%d:%.3f ",j,gauge[i].elem(j));
  //     }
  //     printf("\n");
  //   }
  // }
  
  // Matrix orig[8*volnodesites];
  // for(int i=0;i<8*volnodesites;i++) orig[i] = gauge[i];

  // lattice.Convert(WILSON);
  // lattice.Convert(CANONICAL);

  // gauge = lattice.GaugeField();
  
  // for(int i=0;i<8*volnodesites;i++){
  //   bool linepass =true;
  //   for(int j=0;j<18;j++){
  //     if(gauge[i].elem(j) != orig[i].elem(j)){
  // 	printf("Err %d %d; %e %e\n",i,j,gauge[i].elem(j),orig[i].elem(j));
  // 	linepass = false;
  //     }
  //   }
  //   if(linepass) printf("Matrix %d pass\n",i);
  // }
  

  // exit(0);
  // printf("Post convert\n");

  // printf("cb0:\n");
  // int cbsize = 2*volnodesites;

  // for(int i=0;i<cbsize;i++){
  //   printf("U  i=%d ",i);
  //   for(int j=0;j<18;j++){
  //     if(gauge[i].elem(j)>=0)
  // 	printf("%d: %.3f ",j,gauge[i].elem(j));
  //     else
  // 	printf("%d:%.3f ",j,gauge[i].elem(j));
  //   }
  //   printf("\n");

  //   if(GJP.Gparity()){
  //     printf("U* i=%d ",i+cbsize);
  //     for(int j=0;j<18;j++){
  // 	if(gauge[i+cbsize].elem(j)>=0)
  // 	  printf("%d: %.3f ",j,gauge[i+cbsize].elem(j));
  // 	else
  // 	  printf("%d:%.3f ",j,gauge[i+cbsize].elem(j));
  //     }
  //     printf("\n\n");

  //   }


  // }
  // printf("cb1:\n");
  // int offset = cbsize;
  // if(GJP.Gparity()) offset*=2;

  // for(int i=offset;i<offset+cbsize;i++){
  //   printf("U  i=%d ",i);
  //   for(int j=0;j<18;j++){
  //     if(gauge[i].elem(j)>=0)
  // 	printf("%d: %.3f ",j,gauge[i].elem(j));
  //     else
  // 	printf("%d:%.3f ",j,gauge[i].elem(j));
  //   }
  //   printf("\n");

  //   if(GJP.Gparity()){
  //     printf("U* i=%d ",i+cbsize);
  //     for(int j=0;j<18;j++){
  // 	if(gauge[i+cbsize].elem(j)>=0)
  // 	  printf("%d: %.3f ",j,gauge[i+cbsize].elem(j));
  // 	else
  // 	  printf("%d:%.3f ",j,gauge[i+cbsize].elem(j));
  //     }
  //     printf("\n\n");

  //   }


  // }


  // printf("Gauge field:\n\n");

  // int base_loc[4];
  // for(int i=0;i<4;i++) base_loc[i] = GJP.NodeCoor(i)*GJP.NodeSites(i); //absolute pos of start point
    
  // //int pos[4];

  // for(int t=0;t<GJP.NodeSites(3);t++){
  //   pos[3] = t;
  //   for(int z=0;z<GJP.NodeSites(2);z++){
  //     pos[2] = z;
  //     for(int y=0;y<GJP.NodeSites(1);y++){
  // 	pos[1] = y;
  // 	for(int x=0;x<GJP.NodeSites(0);x++){
  // 	  pos[0] = x;

  // 	  for(int mu=0;mu<4;mu++){
  // 	    Matrix *U = const_cast<Matrix *>(lattice.GetLink(pos, mu));
  // 	    printf("%d %d %d %d; %d:\n",x+base_loc[0],y+base_loc[1],z+base_loc[2],t+base_loc[3],mu);
  // 	    for(int i=0;i<3;i++){
  // 	    	for(int j=0;j<3;j++){
  // 	    	  const Complex &c = (*U)(i,j);
  // 	    	  printf("[%.3f %.3f] ",c.real(), c.imag());
  // 	    	}
  // 	    	printf("\n");
  // 	    }
  // 	    printf("\n");
  // 	  }
	    

  // 	}
  //     }
  //   }
  // }

  // if(GJP.Gparity()){
  //   printf("Conjugated Gauge field:\n\n");

  //   for(int t=0;t<GJP.NodeSites(3);t++){
  //     pos[3] = t;
  //     for(int z=0;z<GJP.NodeSites(2);z++){
  // 	pos[2] = z;
  // 	for(int y=0;y<GJP.NodeSites(1);y++){
  // 	  pos[1] = y;
  // 	  for(int x=0;x<GJP.NodeSites(0);x++){
  // 	    pos[0] = x;
	    
  // 	    for(int mu=0;mu<4;mu++){
  // 	      Matrix *U = const_cast<Matrix *>(lattice.GetLink(pos, mu, 1));
  // 	      printf("%d %d %d %d; %d:\n",x+base_loc[0],y+base_loc[1],z+base_loc[2],t+base_loc[3],mu);
  // 	      for(int i=0;i<3;i++){
  // 	    	for(int j=0;j<3;j++){
  // 	    	  const Complex &c = (*U)(i,j);
  // 	    	  printf("[%.3f %.3f] ",c.real(), c.imag());
  // 	    	}
  // 	    	printf("\n");
  // 	      }
  // 	      printf("\n");
  // 	    }
	    

  // 	  }
  // 	}
  //     }
  //   }

  // }



  // // PropagatorContainer prop;
  // //  prop.setup(prop_args.props.props_val[0]);

   

  //  QPropW *src = &  PropManager::getProp(prop_args.props.props_val[0].generics.tag).getProp(lattice);  

  //  //QPropW *src = &prop.getProp(lattice);


  // // exit(0);

  // // CgArg cg;
  // // cg.mass =  0.5;
  // // cg.max_num_iter = 5000;
  // // cg.stop_rsd =   1.0000000000000000e-06;
  // // cg.true_rsd =   1.0000000000000000e-06;
  // // cg.RitzMatOper = NONE;
  // // cg.Inverter = CG;
  // // cg.bicgstab_n = 0;


  // // QPropWArg qpropw_arg;
  // // qpropw_arg.cg = cg;
  // // qpropw_arg.x = 0;
  // // qpropw_arg.y = 0;
  // // qpropw_arg.z = 0;
  // // qpropw_arg.t = 0;
  // // if(gparity){
  // //   qpropw_arg.flavor = 0; //point source on d field
  // // }
  // // qpropw_arg.ensemble_label = "ens";
  // // qpropw_arg.ensemble_id = "ens_id";
  // // qpropw_arg.StartSrcSpin = 0;
  // // qpropw_arg.EndSrcSpin = 4;
  // // qpropw_arg.StartSrcColor = 0;
  // // qpropw_arg.EndSrcColor = 3;

  // // qpropw_arg.save_prop = 1;
  // // qpropw_arg.file = "prop.dat";
  
  // // CommonArg common_arg("label","filename");
  
  // // QPropW* src = new QPropWPointSrc(lattice, &qpropw_arg, &common_arg);
  

  
  

  // if(UniqueID()==0) printf("Inversion finished\n");

  // //int id = UniqueID();
  // //char propfile[50];
  // //sprintf(&propfile,"prop.$d.dat",id);

  // FILE * ftxt = Fopen(ADD_ID,"prop_txt","w");

  // //int pos[4];
  
  // int shift_x = GJP.XnodeCoor()*GJP.XnodeSites();
  // int shift_y = GJP.YnodeCoor()*GJP.YnodeSites();
  // int shift_z = GJP.ZnodeCoor()*GJP.ZnodeSites();
  // int shift_t = GJP.TnodeCoor()*GJP.TnodeSites();
  // //Local lattice dimensions:
  // int size_x = GJP.XnodeSites();
  // int size_y = GJP.YnodeSites();
  // int size_z = GJP.ZnodeSites();
  // int size_t = GJP.TnodeSites();
  // int size_xy = size_x*size_y;
  // int vol = (GJP.VolNodeSites()/GJP.TnodeSites()); // =size_x*size_y_size_z
  // //Global lattice dimensions
  // int Size_X = GJP.Xnodes()*GJP.XnodeSites();
  // int Size_Y = GJP.Ynodes()*GJP.YnodeSites();
  // int Size_Z = GJP.Znodes()*GJP.ZnodeSites();
  // int Size_T = GJP.Tnodes()*GJP.TnodeSites();


  // for (int i=0; i<GJP.VolNodeSites(); i++) {
  //   //Global coordinates
  //   pos[3] = i/vol + shift_t;
  //   pos[2] = (i%vol)/size_xy + shift_z;
  //   pos[1] = (i%size_xy)/size_x + shift_y;
  //   pos[0] = i%size_x + shift_x;
    
  //   //printf("site (%d,%d,%d,%d)\n",pos[0],pos[1],pos[2],pos[3]);  
  //   Fprintf(ADD_ID,ftxt,"site (%d,%d,%d,%d)\n",pos[0],pos[1],pos[2],pos[3]);
    
  //   WilsonMatrix &m = (*src)[i];
  //   for(int s1=0;s1<4;s1++){
  //     for(int c1=0;c1<3;c1++){
  // 	for(int s2=0;s2<4;s2++){
  // 	  for(int c2=0;c2<3;c2++){
  // 	    Complex &val = m(s1,c1,s2,c2);
  // 	    IFloat re = val.real();
  // 	    //if(re!=0.0){
  // 	      //printf("u %d %d %d %d : %f %f\n",s1,c1,s2,c2,val.real(),val.imag());
  // 	      Fprintf(ADD_ID,ftxt,"u %d %d %d %d : %f %f\n",s1,c1,s2,c2,val.real(),val.imag());
  // 	      //}
  // 	  }
  // 	}
  //     }
  //   }
  // }
  // if(GJP.Gparity()){
  //   for (int i=0; i<GJP.VolNodeSites(); i++) {
  //     pos[3] = i/vol + shift_t;
  //     pos[2] = (i%vol)/size_xy + shift_z;
  //     pos[1] = (i%size_xy)/size_x + shift_y;
  //     pos[0] = i%size_x + shift_x;

  //     //printf("site (%d,%d,%d,%d)\n",pos[0],pos[1],pos[2],pos[3]);  
  //     Fprintf(ADD_ID,ftxt,"site (%d,%d,%d,%d)\n",pos[0],pos[1],pos[2],pos[3]);

  //     WilsonMatrix &m = (*src)[i+GJP.VolNodeSites()];
  //     for(int s1=0;s1<4;s1++){
  // 	for(int c1=0;c1<3;c1++){
  // 	  for(int s2=0;s2<4;s2++){
  // 	    for(int c2=0;c2<3;c2++){
  // 	      Complex &val = m(s1,c1,s2,c2);
  // 	      IFloat re = val.real();
  // 	      //if(re!=0.0){
  // 		//printf("d %d %d %d %d : %f %f\n",s1,c1,s2,c2,val.real(),val.imag());
  // 		Fprintf(ADD_ID,ftxt,"d %d %d %d %d : %f %f\n",s1,c1,s2,c2,val.real(),val.imag());
  // 		//}
  // 	    }
  // 	  }
  // 	}
  //     }
  //   }
  // }

  // Fclose(ADD_ID,ftxt);







// class PropagatorContainer{
// protected:
//   QPropW *prop;
//   AttributeContainer* attributes[50]; //positions are mapped by the integer values of the enum  AttrType with no duplicate entries
//   AttributeContainer* findAttr(const AttrType &type) const; //returns NULL if not present
// public:
//   void add(const AttributeContainer &p);
//   void setup(PropagatorArg &arg);
//   template<typename T>
//   bool getAttr(T* &to) const;
  
//   void readProp(); //loads prop if on disk, does nothing otherwise
//   void calcProp(Lattice &latt);
//   void deleteProp(); //deletes the prop from memory, used to save space. Will be automatically recalculated if getProp is called again

//   QPropW & getProp(Lattice &latt); //get the prop, calculate or load if necessary

//   bool tagEquals(const char* what); //check if tag is equal to input

//   int flavor() const;

//   PropagatorContainer();
//   ~PropagatorContainer();
// };

// PropagatorContainer::PropagatorContainer(): prop(NULL){ for(int i=0;i<50;i++) attributes[i]=NULL;}


// void PropagatorContainer::setup(PropagatorArg &arg){
//   for(int i=0;i<50;i++) if(attributes[i]) delete attributes[i];
  
//   //special case, generic attributes are always present and maintained outside of the array
//   AttributeContainer* p = new AttributeContainer;
//   p->type = GENERIC_PROP_ATTR;
//   p->AttributeContainer_u.generic_prop_attr = arg.generics;
//   attributes[ (int)p->type ] = p;

//   printf("Got %d additional attributes\n",arg.attributes.attributes_len);

//   for(int i=0;i<arg.attributes.attributes_len;i++){
//     //no duplicates, later entries of same type overwrite earlier
//     p = &arg.attributes.attributes_val[i];
//     int idx = (int)p->type;
//     if(attributes[idx]) delete attributes[idx];
//     attributes[idx] = new AttributeContainer(*p); //make a copy
//   }
// }

// void PropagatorContainer::add(const AttributeContainer &p){
//   int idx = (int)p.type;
//   if(attributes[idx]!=NULL) delete attributes[idx]; //no duplicates
//   attributes[idx] = new AttributeContainer(p);
// }

// AttributeContainer*  PropagatorContainer::findAttr(const AttrType &type) const{
//   return attributes[(int)type];
// }

// template<typename T>
// bool PropagatorContainer::getAttr(T* &to) const{
//   AttributeContainer *c = findAttr(T::getType());
//   if(c==NULL) return false;
//   to = reinterpret_cast<T*>(&(c->AttributeContainer_u)); //C99 standard: A pointer to a union object, suitably converted, points to each of its members
//   return true;
// }

// PropagatorContainer::~PropagatorContainer(){
//   for(int i=0;i<50;i++) if(attributes[i]) delete attributes[i];
//   if(prop!=NULL) delete prop;
// }


// void PropagatorContainer::readProp(){
//   if(prop) return; //don't load if already inverted
//   PropIOAttrArg *io;
//   if(!getAttr(io)) return;
//   if(!io->prop_on_disk) return;

//   //load prop with info from io
//   prop->ReLoad(io->qio_filename);
// }

// void PropagatorContainer::calcProp(Lattice &latt){
//   //function acts as factory for QPropW objects depending on the attribute objects
//   if(prop) return; //don't calculate twice

//   char *cname = "PropagatorContainer";
//   char *fname = "calcProp()";

//   CommonArg c_arg("label","filename");//find out what this does!

//   GenericPropAttrArg *generics;
//   if(!getAttr(generics)) ERR.General(cname,fname,"Propagator attribute list does not contain a GenericPropAttr\n");

//   printf("Calculating propagator %s\n",generics->tag);

//   CgArg cg;
//   cg.mass =  generics->mass;
//   cg.max_num_iter = 5000;
//   cg.stop_rsd =   1.0000000000000000e-06;
//   cg.true_rsd =   1.0000000000000000e-06;
//   cg.RitzMatOper = NONE;
//   cg.Inverter = CG;
//   cg.bicgstab_n = 0;

//   QPropWArg qpropw_arg;
//   qpropw_arg.cg = cg;
//   qpropw_arg.x = 0;
//   qpropw_arg.y = 0;
//   qpropw_arg.z = 0;
//   qpropw_arg.t = 0;
//   qpropw_arg.flavor = 0; //default on d field  
//   qpropw_arg.ensemble_label = "ens";
//   qpropw_arg.ensemble_id = "ens_id";
//   qpropw_arg.StartSrcSpin = 0;
//   qpropw_arg.EndSrcSpin = 4;
//   qpropw_arg.StartSrcColor = 0;
//   qpropw_arg.EndSrcColor = 3;

//   //fill out qpropw_arg arguments
//   GparityFlavorAttrArg *flav;
//   if(getAttr(flav)) qpropw_arg.flavor = flav->flavor;

//   PropIOAttrArg *io;
//   if(getAttr(io)){ 
//     qpropw_arg.save_prop = io->save_to_disk;
//     qpropw_arg.file = io->qio_filename;
//   }

//   PointSourceAttrArg *pt;
//   WallSourceAttrArg *wl;
//   if(getAttr(pt) && getAttr(wl) )
//     ERR.General(cname,fname,"Propagator %s cannot have both WallSourceAttrArg and PointSourceAttrArg\n",generics->tag);

//   if(getAttr(pt)){
//     qpropw_arg.x = pt->pos[0];
//     qpropw_arg.y = pt->pos[1];
//     qpropw_arg.z = pt->pos[2];
//     qpropw_arg.t = pt->pos[3];
    
//     //assemble the arg objects
//     prop = new QPropWPointSrc(latt,&qpropw_arg,&c_arg);

//     MomentumAttrArg *mom;
//     if(getAttr(mom)){
//       //for a point source with momentum, apply the e^-ipx factor at the source location 
//       //such that it can be treated in the same was as a momentum source
//       const Float Pi_const=3.141592654;
//       Float pdotx=0.0;
//       pdotx +=((Float) mom->p[0]*pt->pos[0])/((Float) GJP.XnodeSites()*GJP.Xnodes() )
// 	    + ((Float) mom->p[1]*pt->pos[1])/((Float) GJP.YnodeSites()*GJP.Ynodes() )
//       	    + ((Float) mom->p[2]*pt->pos[2])/((Float) GJP.ZnodeSites()*GJP.Znodes() );
//       Rcomplex mom_fac(cos(pdotx),sin(pdotx));
//       for(int i=0;i<GJP.VolNodeSites();i++) (*prop)[i] *= mom_fac; 
//     }
//   }else if(getAttr(wl)){
//     qpropw_arg.t = wl->t;
    
//     MomentumAttrArg *mom;
//     if(getAttr(mom)){
//       prop = new QPropWMomSrc(latt,&qpropw_arg,mom->p,&c_arg);
//     }else{
//       prop = new QPropWWallSrc(latt,&qpropw_arg,&c_arg);
//     }
//   }else{
//     ERR.General(cname,fname,"Propagator %s has no source type AttrArg\n",generics->tag);
//   }
  
//   if(getAttr(io)){
//     io->prop_on_disk = true; //QPropW saves the prop
//   }
// }


// int PropagatorContainer::flavor() const{
//   GparityFlavorAttrArg *flav;
//   if(getAttr(flav)) return flav->flavor;
//   return 0;
// }

// QPropW & PropagatorContainer::getProp(Lattice &latt){
//   if(!prop){
//     readProp();
//     calcProp(latt); //will calculate if prop was not read
//   }
//   return *prop;
// }

// bool PropagatorContainer::tagEquals(const char* what){
//   GenericPropAttrArg *generics;
//   if(!getAttr(generics)) ERR.General("PropagatorContainer","tagEquals(const char* what)","Propagator attribute list does not contain a GenericPropAttr\n");
//   if(strcmp(generics->tag,what)==0) return true;
//   return false;
// }

// void PropagatorContainer::deleteProp(){
//   if(prop) delete prop;
// }


// //container for multiple props with destructor that deletes all props
// class PropVector{
// public:
//   const static int MAX_SIZE = 100;
//   PropagatorContainer & operator[](const int &idx);
//   PropagatorContainer & addProp(PropagatorArg &arg);
//   const int &size() const;
//   PropVector();
//   ~PropVector();
// private:
//   PropagatorContainer *props[MAX_SIZE];
//   int sz;
// };
// PropagatorContainer & PropVector::operator[](const int &idx){ return *props[idx]; }
// PropVector::PropVector(): sz(0){ for(int i=0;i<MAX_SIZE;i++) props[i] = NULL; }
// PropVector::~PropVector(){ for(int i=0;i<MAX_SIZE;i++) if(props[i]!=NULL) delete props[i]; }
  
// const int &PropVector::size() const{ return sz; }

// PropagatorContainer & PropVector::addProp(PropagatorArg &arg){
//   if(sz==MAX_SIZE){ ERR.General("PropVector","addProp(PropagatorArg &arg)","Reached maximum number of allowed propagators: %d\n",MAX_SIZE); }
//   PropagatorContainer* p = new PropagatorContainer;
//   p->setup(arg);
//   props[sz++] = p;
//   return *p;
// }

// //static class for managing propagators
// class PropManager{
// private:
//   static PropVector props;
// public:
//   static PropagatorContainer & getProp(const char *tag);
//   static PropagatorContainer & addProp(PropagatorArg &arg); //this does not actually invert the propagator until getProp is called on the PropagatorContainer or using the above function on the vector
//   static void setup(JobPropagatorArgs &prop_args);
// };
// PropVector PropManager::props;

// PropagatorContainer & PropManager::getProp(const char *tag){
//   for(int i=0;i<props.size();i++) if(props[i].tagEquals(tag)) return props[i];
//   ERR.General("PropManager","getProp(const char *tag, Lattice &latt)","Prop '%s' does not exist!\n",tag);
// };

// PropagatorContainer & PropManager::addProp(PropagatorArg &arg){
//   return props.addProp(arg);
// }

// void PropManager::setup(JobPropagatorArgs &prop_args){
//   for(int i=0;i< prop_args.props.props_len; i++){
//     addProp(prop_args.props.props_val[i]);
//   }
// }


// class CorrelationFunction{
//   char *label;
//   int ncontract;
//   int time_size;
//   Rcomplex** wick; //the wick contractions, each a function of t

//   void sumLattice(); //sum each element of wick[contraction][t] over all nodes
// public:
//   void setNcontractions(const int &n);
//   void write(const char *file);
//   void write(FILE *fp);
//   Rcomplex & operator()(const int &contraction_idx, const int &t);
//   CorrelationFunction(const char *_label);
//   ~CorrelationFunction();
// };
// CorrelationFunction::CorrelationFunction(const char *_label): wick(NULL),ncontract(0){
//   time_size=GJP.Tnodes()*GJP.TnodeSites();
//   label = new char[strlen(_label)];
//   strcpy(label,_label);
// }
// void CorrelationFunction::setNcontractions(const int &n){
//   if(ncontract!=0){
//     sfree(wick[0]);
//     sfree(wick);
//   }
//   ncontract = n;
//   wick = (Rcomplex**)smalloc(n*sizeof(Rcomplex*)); //allocate array of pointers (I wish we could just use std::vector here!)

//   //contiguous memory slot
//   Rcomplex *stack = (Rcomplex*)smalloc(n*time_size*sizeof(Rcomplex));

//   for(int i=0;i<n;i++){
//     wick[i] = stack+i*time_size;
//     for(int t=0;t<time_size;t++){ wick[i][t].real(0.0); wick[i][t].imag(0.0); }
//   }
// }
// CorrelationFunction::~CorrelationFunction(){
//   if(ncontract!=0){
//     sfree(wick[0]); //free the contiguous memory block holding the data
//     sfree(wick); //delete the pointer array holding this
//   }
// }
// Rcomplex & CorrelationFunction::operator()(const int &contraction_idx, const int &t){
//   return wick[contraction_idx][t];
// }
// void CorrelationFunction::sumLattice(){
//   for(int i=0;i<ncontract;i++){
//     for(int t=0;t<time_size;t++){
//       slice_sum( (Float*)&wick[i][t], 2, 99); //2 for re/im, 99 is a *magic* number (we are abusing slice_sum here)
//     }
//   }
// }
// void CorrelationFunction::write(const char *file){
//   FILE *fp;
//   if ((fp = Fopen(file, "w")) == NULL) {
//     ERR.FileW("CorrelationFunction","write(const char *file)",file);
//   }
//   write(fp);
//   Fclose(fp);
// }
// void CorrelationFunction::write(FILE *fp){
//   sumLattice(); //sum the correlation function over all nodes

//   Fprintf(fp,"%s\n",label);
//   Fprintf(fp,"%d contractions\n",ncontract);

//   for(int c=0;c<ncontract;c++){
//     Fprintf(fp,"Contraction %d\n",c);

//     for(int t=0; t<time_size; t++)
//       Fprintf(fp, "%d %.16e %.16e\n",t, wick[c][t].real(), wick[c][t].imag());
//   }
// }


  
// class AlgGparityContract : public Alg{
// private:
//   char *cname;
//   GparityContractArg *args;
// public:
//   AlgGparityContract(Lattice & latt, CommonArg& c_arg, GparityContractArg& arg);
//   void run();
//   void spectrum(const GparityMeasurement &measargs);

//   void pi_plus(const char *q_f0_tag, const char *q_f1_tag, CorrelationFunction &corr);
//   void pi_minus(const char *q_f0_tag, const char *q_f1_tag, CorrelationFunction &corr);
// };
// AlgGparityContract::AlgGparityContract(Lattice & latt, CommonArg& c_arg, GparityContractArg& arg): Alg(latt,&c_arg), args(&arg){ cname = "AlgGparityContract"; }

// void AlgGparityContract::run(){
//   for(int i=0;i<args->meas.meas_len;i++){
//     spectrum(args->meas.meas_val[i]);   
//   }
// }

// void AlgGparityContract::spectrum(const GparityMeasurement &measargs){
//   char * prop_f[3]; bool got_f[3] = {false,false,false};

//   PropagatorContainer &p1 = PropManager::getProp(measargs.prop_1);
//   PropagatorContainer &p2 = PropManager::getProp(measargs.prop_2);
  
//   int f1 = p1.flavor(); int f2 = p2.flavor();
//   prop_f[f1] = measargs.prop_1; got_f[f1] = true;
//   prop_f[f2] = measargs.prop_2; got_f[f2] = true;

//   char out[100];

//   if(got_f[0] && got_f[1]){
//     sprintf(out,"%s pi^+",measargs.label_stub);
//     CorrelationFunction c_piplus(out);
//     sprintf(out,"%s pi^-",measargs.label_stub);
//     CorrelationFunction c_piminus(out);

//     if(UniqueID()==0) printf("Doing pi^+ contractions\n");
//     pi_plus(prop_f[0],prop_f[1],c_piplus);
//     if(UniqueID()==0) printf("Doing pi^- contractions\n");
//     pi_minus(prop_f[0],prop_f[1],c_piminus);
    
//     sprintf(out,"%s_pi_plus.dat",measargs.file_stub);
//     c_piplus.write(out);
//     sprintf(out,"%s_pi_minus.dat",measargs.file_stub);
//     c_piminus.write(out);
//   }else{
//     ERR.General(cname, "spectrum(const char *prop_1, const char *prop_2)", "Could not find any valid contractions for propagators %s and %s\n",measargs.prop_1,measargs.prop_2);
//   }

// }

// void AlgGparityContract::pi_plus(const char *q_f0_tag, const char *q_f1_tag, CorrelationFunction &corr){
//   int shift_t = GJP.TnodeCoor()*GJP.TnodeSites();
//   int spatial_vol = GJP.VolNodeSites()/GJP.TnodeSites();

//   corr.setNcontractions(2);

//   WilsonMatrix trace_a, trace_b, tmp;
//   PropagatorContainer &q_f0_pc = PropManager::getProp(q_f0_tag);
//   PropagatorContainer &q_f1_pc = PropManager::getProp(q_f1_tag);

//   QPropW &q_f0_qpw = q_f0_pc.getProp(AlgLattice());
//   QPropW &q_f1_qpw = q_f1_pc.getProp(AlgLattice());
  
//   Rcomplex tmpc;

//   for(int i=0;i<GJP.VolNodeSites();i++){
//     int t = i/spatial_vol + shift_t;

//     // //DEBUG
//     // int shift_x = GJP.XnodeCoor()*GJP.XnodeSites();
//     // int shift_y = GJP.YnodeCoor()*GJP.YnodeSites();
//     // int shift_z = GJP.ZnodeCoor()*GJP.ZnodeSites();
//     // int shift_t = GJP.TnodeCoor()*GJP.TnodeSites();
//     // //Local lattice dimensions:
//     // int size_x = GJP.XnodeSites();
//     // int size_y = GJP.YnodeSites();
//     // int size_z = GJP.ZnodeSites();
//     // int size_t = GJP.TnodeSites();
//     // int size_xy = size_x*size_y;
//     // int vol = (GJP.VolNodeSites()/GJP.TnodeSites()); // =size_x*size_y_size_z
    
//     // int pos[4];
//     // pos[3] = i/vol + shift_t;
//     // pos[2] = (i%vol)/size_xy + shift_z;
//     // pos[1] = (i%size_xy)/size_x + shift_y;
//     // pos[0] = i%size_x + shift_x;

//     // printf("Piplus (%d,%d,%d,%d)\n",pos[0],pos[1],pos[2],pos[3]);
//     // //DEBUG


//     //Doing contraction 1 of 2
//     Rcomplex &wick0 = corr(0,t);
//     tmpc.real(1.0); tmpc.imag(0.0);

//     //Doing trace 1 of 1
//     trace_a = q_f0_qpw.SiteMatrix(i,0);
//     trace_a.hconj();
//     trace_a.gl(3).gl(1).gl(-5);
//     trace_a.gr(-5).gr(1).gr(3);
//     trace_b = q_f1_qpw.SiteMatrix(i,1);
//     trace_b.cconj();
//     tmpc*= Trace(trace_a,trace_b);

//     wick0 += tmpc;
//     //printf("Wick0 -> %e,%e\n",wick0.real(),wick0.imag());

//     //Doing contraction 2 of 2
//     Rcomplex &wick1 = corr(1,t);
//     tmpc.real(1.0); tmpc.imag(0.0);

//     //Doing trace 1 of 1
//     trace_a = q_f1_qpw.SiteMatrix(i,0);
//     trace_a.hconj();
//     trace_a.gl(3).gl(1).gl(-5);
//     trace_a.gr(-5).gr(1).gr(3);
//     trace_b = q_f0_qpw.SiteMatrix(i,1);
//     trace_b.cconj();
//     tmpc*= Trace(trace_a,trace_b);

//     wick1 += tmpc;
//     //printf("Wick1 -> %e,%e\n",wick1.real(),wick1.imag());
//   }
// }

// void AlgGparityContract::pi_minus(const char *q_f0_tag, const char *q_f1_tag, CorrelationFunction &corr){
//   int shift_t = GJP.TnodeCoor()*GJP.TnodeSites();
//   int spatial_vol = GJP.VolNodeSites()/GJP.TnodeSites();

//   corr.setNcontractions(2);

//   WilsonMatrix trace_a, trace_b, tmp;
//   PropagatorContainer &q_f0_pc = PropManager::getProp(q_f0_tag);
//   PropagatorContainer &q_f1_pc = PropManager::getProp(q_f1_tag);

//   QPropW &q_f0_qpw = q_f0_pc.getProp(AlgLattice());
//   QPropW &q_f1_qpw = q_f1_pc.getProp(AlgLattice());

//   Rcomplex tmpc;
  
//   for(int i=0;i<GJP.VolNodeSites();i++){
//     int t = i/spatial_vol + shift_t;

//     //Doing contraction 1 of 2
//     Rcomplex &wick0 = corr(0,t); 
//     tmpc.real(1.0); tmpc.imag(0.0);

//     //Doing trace 1 of 1
//     trace_a = q_f0_qpw.SiteMatrix(i,0);
//     trace_a.gl(3).gl(1).gl(-5);
//     trace_a.gr(-5).gr(1).gr(3);
//     trace_b = q_f1_qpw.SiteMatrix(i,1);
//     trace_b.transpose();
//     tmpc*= Trace(trace_a,trace_b);

//     wick0 += tmpc;

//     //Doing contraction 2 of 2
//     Rcomplex &wick1 = corr(1,t); 
//     tmpc.real(1.0); tmpc.imag(0.0);

//     //Doing trace 1 of 1
//     trace_a = q_f0_qpw.SiteMatrix(i,1);
//     trace_a.gr(-5).gr(1).gr(3);
//     trace_b = q_f1_qpw.SiteMatrix(i,0);
//     trace_b.transpose();
//     trace_b.gr(-5).gr(1).gr(3);
//     tmpc*= Trace(trace_a,trace_b);

//     wick1 += tmpc;
//   }
// }

