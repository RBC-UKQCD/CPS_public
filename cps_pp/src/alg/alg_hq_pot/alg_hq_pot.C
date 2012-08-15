#include <config.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/verbose.h>
#include <util/error.h>
#include <alg/alg_fix_gauge.h>
#include <alg/do_arg.h>
#include <alg/hmd_arg.h>
#include <alg/common_arg.h>
#include <util/random.h>
#include <math.h>
#include <util/qcdio.h>
#include <util/lattice.h>
#include <util/smalloc.h>
#include <util/vector.h>
#include <comms/scu.h>
#include <comms/glb.h>
#include <alg/alg_plaq.h>
#include <alg/alg_hq_pot.h>
#include <util/time_cps.h>
#include <util/pt.h>
#include <string.h>

#ifdef PARALLEL
#include <comms/sysfunc_cps.h>
#endif

CPS_START_NAMESPACE

#define X_LINK 0
#define Y_LINK 1
#define Z_LINK 2
#define T_LINK 3

Complex I = Complex(0,1);

//inline Matrix operator * (const Matrix& m1, const Matrix& m2)
//{ Matrix r; r.DotMEqual(m1,m2); return r; }

void p(Matrix x)
{
  for(int i=0; i<3; i++)
    {
      for(int j=0; j<3; j++)
	{
	  Complex xx = x(i,j);
	  if(fabs(real(xx))<1e-10)
	    xx=Complex(0,imag(xx));
	  if(fabs(imag(xx))<1e-10)
	    xx=Complex(real(xx),0);
	  printf( "(%e,%e)", real(xx), imag(xx));
	  if(j<2)
	    printf( "       ");
	}
      printf("\n");
    }
  printf("\n");
}



#define NX (GJP.XnodeSites())
#define NY (GJP.YnodeSites())
#define NZ (GJP.ZnodeSites())
#define NT (GJP.TnodeSites())

#define INDp(x,y,z,t) (((((t+NT)%NT)*NZ+((z+NZ)%NZ))*NY+   \
			  ((y+NY)%NY))*NX+((x+NX)%NX))

#define IND(x,y,z,t,l) ((((((t+NT)%NT)*NZ+((z+NZ)%NZ))*NY+   \
			  ((y+NY)%NY))*NX+((x+NX)%NX))*4+l)


//------------------------------------------------------------------
/*!  
  \param latt The Lattice object containg the gauge field on which to compute the %Wilson lines.
  \param c_arg Container for generic parameters. .
  \param arg Empty parameter container.
*/
//------------------------------------------------------------------
AlgHQPotential::AlgHQPotential(Lattice& latt, 
	     CommonArg *c_arg,
	     NoArg *arg) : 
	     Alg(latt, c_arg) 
{
  cname = "AlgHQPotential";
  const char *fname = "AlgHQPotential";
  VRB.Func(cname,fname);

  // Initialize the argument pointer
  //----------------------------------------------------------------
  if(arg == 0)
    ERR.Pointer(cname,fname, "arg");
  alg_HQPotential_arg = arg;

}


//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
AlgHQPotential::~AlgHQPotential() {
  const char *fname = "~AlgHQPotential";
  VRB.Func(cname,fname);
}

//-----------------------------------------------------------------
// Calculating the wline correlators, storing the dist-correlators
// and also storing the weight of each distances for later calculation.
// Here dir->direction of coulomb gauge fixing, n -> number of separation
// to be calculated, m -> starting t plane
//-----------------------------------------------------------------
void AlgHQPotential::run(int dir, int n, int m)
{
  const char *fname = "run";
  VRB.Func(cname,fname);

  if ((dir<0)||(dir>3)){
    printf("Error:: direction should be 0,1,2,3\n");
    return;
  }
  
  // Set the Lattice pointer
  //----------------------------------------------------------------
  Lattice& lat = AlgLattice();
  
  // Set up the gauge fixing and other required conditions
  //---------------------------------------------------------------

  int num_nodes[4]
    = { GJP.Xnodes(), GJP.Ynodes(), GJP.Znodes(), GJP.Tnodes() } ;
  
  int node_sites[4]
    = { GJP.XnodeSites(), GJP.YnodeSites(), GJP.ZnodeSites(), GJP.TnodeSites() } ;
  
  int Size[4];
  for (int i=0; i<4; i++){
    Size[i] = num_nodes[i]*node_sites[i];
  }
  
  int *plns;
  plns = (int*) smalloc(Size[dir]*sizeof(int));
  
  for (int i=0; i<Size[dir]; i++){
    plns[i] = i;
  }
  
  int npln = Size[dir];
  
  FixGaugeType normdir;
  
  if (dir==3) normdir = FIX_GAUGE_COULOMB_T;
  else if (dir==0) normdir = FIX_GAUGE_COULOMB_X;
  else if (dir==1) normdir = FIX_GAUGE_COULOMB_Y;
  else normdir = FIX_GAUGE_COULOMB_Z;
      
  //----------------------------------------------------------------------
  //initialize the parameters need to gauge fixing ---------------------
  //----------------------------------------------------------------------
  Matrix *L;
  Matrix **Gp;
  int ii;
  
  int volume = NX*NY*NZ*NT;
  
  L = (Matrix*) smalloc(4*volume*sizeof(Matrix));
  Gp = (Matrix**) smalloc(npln*sizeof(Matrix*));
  
  for(ii=0; ii<npln; ii++)
    Gp[ii] = (Matrix*) smalloc(volume/node_sites[dir] * sizeof(Matrix));
  
  
  //-----------------------------------------------------------------------------
  //GAUGE FIXING
  //-----------------------------------------------------------------------------
  Float bf_gf_time = dclock();
  
  VRB.Debug("Fixing . . .\n");
  
  lat.FixGaugeAllocate(normdir,npln,plns);
  int itnum = lat.FixGauge(1e-6,20000);
  VRB.Debug("Iternum = %d\n", itnum);
  
  VRB.Debug("Resulting Gauge Fixing Matrices:\n\n");
  
  if(itnum > 0)
    for (int slice=0; slice<node_sites[dir]; slice++)
      for(int cnt=0; cnt<volume/node_sites[dir]; cnt++)
	Gp[slice][cnt]=lat.FixGaugePtr()[slice][cnt];
  
  p(Gp[0][0]);

  lat.FixGaugeFree();
  
  //-------------------------------------------------------------------------------
  //-------------------------------------------------------------------------------
  // TRY TO Transform the Lattice to the Coulomb Gauge and store it ---------------
  //-------------------------------------------------------------------------------

  int tt,xx,yy,zz,temp_ind;
  int slice_ind[3];
  
  //slice_ind[3] stores the 3 directions on the 'dir' slice with indices increasing 
  temp_ind = 0;
  for (int i=0; i<4; i++){
    if (i!=dir) {
      slice_ind[temp_ind] = i;
      temp_ind++;
    }
  }
  
  printf("dir == %d \n",dir);
  printf("slice index == %d %d %d\n",slice_ind[0],slice_ind[1],slice_ind[2]);
  
  // the dummy node_sites for each dummy dirction
  int NN[4];
  NN[0] = node_sites[slice_ind[0]];
  NN[1] = node_sites[slice_ind[1]];
  NN[2] = node_sites[slice_ind[2]];
  NN[3] = node_sites[dir];
  
  int s[4];
  for (int i=0; i<3; i++)
    s[i] = slice_ind[i];
  s[3] = dir;

  //---------------------------------------------------------------------------------
  //copy the old lattice config to matrix array L, transform L and then copy back
  //---------------------------------------------------------------------------------
  
  lat.CopyGaugeField(L);  
  int x[4];
  
  // xx yy zz tt are dummy position vector, as tt represents the gfixing direction  
  for (x[3]=0; x[3]<node_sites[3]; x[3]++)
    for (x[2]=0; x[2]<node_sites[2]; x[2]++)
      for (x[1]=0; x[1]<node_sites[1]; x[1]++)
	for (x[0]=0; x[0]<node_sites[0]; x[0]++)
	  {
	    xx = x[slice_ind[0]];
	    yy = x[slice_ind[1]];
	    zz = x[slice_ind[2]];
	    tt = x[dir];
	
	    Matrix g = Gp[tt][(zz*NN[1]+yy)*NN[0]+xx];
	    Matrix D; 
	    
	    Matrix gg ;
	    Matrix transmit;
	    
	    //----------------- T Direction ----------------------------
	    if (tt+1<NN[3]) gg = Gp[tt+1][(zz*NN[1]+yy)*NN[0]+xx];
	    else { 
	      transmit = Gp[0][(zz*NN[1]+yy)*NN[0]+xx];
	      getPlusData((IFloat *)&gg, (IFloat *)&transmit,
			  sizeof(Matrix)/sizeof(IFloat), dir) ;
	    }
	    D.Dagger(gg);
	    L[IND(x[0],x[1],x[2],x[3],dir)] = g*L[IND(x[0],x[1],x[2],x[3],dir)]*D;
	    
	    //----------------- Z Direction ----------------------------
	    if (zz+1<NN[2]) gg = Gp[tt][((zz+1)*NN[1]+yy)*NN[0]+xx]; 
	    else {
	      transmit = Gp[tt][((zz+1)%NN[2]*NN[1]+yy)*NN[0]+xx]; 
	      getPlusData((IFloat *)&gg, (IFloat *)&transmit,
			  sizeof(Matrix)/sizeof(IFloat), slice_ind[2]) ;
	    }
	    D.Dagger(gg);
	    L[IND(x[0],x[1],x[2],x[3],slice_ind[2])] = g*L[IND(x[0],x[1],x[2],x[3],slice_ind[2])]*D;
	    
	    //----------------- Y Direction ----------------------------
	    if (yy+1<NN[1]) gg = Gp[tt][(zz*NN[1]+(yy+1))*NN[0]+xx];
	    else {
	      transmit = Gp[tt][(zz*NN[1]+(yy+1)%NN[1])*NN[0]+xx];
	      getPlusData((IFloat *)&gg, (IFloat *)&transmit,
			  sizeof(Matrix)/sizeof(IFloat), slice_ind[1]) ;
	    }
	    D.Dagger(gg);
	    L[IND(x[0],x[1],x[2],x[3],slice_ind[1])] = g*L[IND(x[0],x[1],x[2],x[3],slice_ind[1])]*D;
	    
	    //----------------- X Direction ----------------------------
	    if (xx+1<NN[0]) gg = Gp[tt][(zz*NN[1]+yy)*NN[0]+(xx+1)];
	    else {
	      transmit = Gp[tt][(zz*NN[1]+yy)*NN[0]+(xx+1)%NN[0]];
	      getPlusData((IFloat *)&gg, (IFloat *)&transmit,
			  sizeof(Matrix)/sizeof(IFloat), slice_ind[0]) ;
	    }
	    D.Dagger(gg);
	    L[IND(x[0],x[1],x[2],x[3],slice_ind[0])] = g*L[IND(x[0],x[1],x[2],x[3],slice_ind[0])]*D;
	  }
  
  
  lat.GaugeField(L);
  
  Float af_gf_time = dclock();
  
  printf("gauge fixing takes %f minutes\n",(double)(af_gf_time-bf_gf_time)/60);

  //--------------------------------Free L and Gp -------------------------------------
  sfree(L);
  
  for (ii=0; ii<npln; ii++)
    sfree(Gp[ii]);
  
  sfree(Gp);
  sfree(plns);
  

  //-----------------------------------------------------------------------------------
  // Compute the Wilson Line for the 'dir' direction with the separation specified
  //-----------------------------------------------------------------------------------
  
  //Matrix *gauge=lat.GaugeField() ;
  
  //----------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------
  int Max = node_sites[0]*num_nodes[0]/2+1;
  Max = (Max>9)?9:Max;
  int Max2 = Max*Max;
  int Max3 = Max*Max2;
  
  Matrix **v;
  Matrix **v1,**v2,**v3,**v4;
  //Matrix *v1,*v2;
  //Matrix **u;
  
  //u  = (Matrix**) smalloc(Max*sizeof(Matrix*));
  
  //v1 = (Matrix*) smalloc(volume*sizeof(Matrix));
  //v2 = (Matrix*) smalloc(volume*sizeof(Matrix));
  
  v  = (Matrix**) smalloc(n*sizeof(Matrix*));
  v1 = (Matrix**) smalloc(n*sizeof(Matrix*));
  v2 = (Matrix**) smalloc(n*sizeof(Matrix*));
  v3 = (Matrix**) smalloc(n*sizeof(Matrix*));
  v4 = (Matrix**) smalloc(n*sizeof(Matrix*));
  
  //for (int i=0; i<Max; i++)
  //u[i]  = (Matrix*) smalloc(volume*sizeof(Matrix));
  
  for (int i=0; i<n; i++){
    v[i]  = (Matrix*) smalloc(volume*sizeof(Matrix));
    v1[i] = (Matrix*) smalloc(volume*sizeof(Matrix));  
    v2[i] = (Matrix*) smalloc(volume*sizeof(Matrix));  
    v3[i] = (Matrix*) smalloc(volume*sizeof(Matrix));  
    v4[i] = (Matrix*) smalloc(volume*sizeof(Matrix));  
  }
  
  for (x[0]=0; x[0]<node_sites[0]; x[0]++)
    for (x[1]=0; x[1]<node_sites[1]; x[1]++)  
      for (x[2]=0; x[2]<node_sites[2]; x[2]++)
	for (x[3]=0; x[3]<node_sites[3]; x[3]++)
	  v1[0][INDp(x[0],x[1],x[2],x[3])] = 1.0;
  
  ParTransGauge pt(lat);
  int DIR;
  printf("ParTransGauge Enterted:::::::");
  

  DIR = 2*dir;

  //--------------wilson lines-------------
  for (int i=0; i<m-1; i++){
    pt.run(1,&(v2[0]),&(v1[0]),(const int*)&DIR);
    // memcpy(v1[0],v2[0],volume*sizeof(Matrix));
    for (x[0]=0; x[0]<node_sites[0]; x[0]++)
      for (x[1]=0; x[1]<node_sites[1]; x[1]++)
	for (x[2]=0; x[2]<node_sites[2]; x[2]++)
	  for (x[3]=0; x[3]<node_sites[3]; x[3]++)
	    v1[0][INDp(x[0],x[1],x[2],x[3])]=v2[0][INDp(x[0],x[1],x[2],x[3])];
  }
  
  for (int i=0; i<n; i++){
    pt.run(1,&(v[i]),&(v1[0]),(const int*)&DIR);
    //  memcpy(v1[0],v[i],volume*sizeof(Matrix));
    for (x[0]=0; x[0]<node_sites[0]; x[0]++)
      for (x[1]=0; x[1]<node_sites[1]; x[1]++)
	for (x[2]=0; x[2]<node_sites[2]; x[2]++)
	  for (x[3]=0; x[3]<node_sites[3]; x[3]++)
	    v1[0][INDp(x[0],x[1],x[2],x[3])]=v[i][INDp(x[0],x[1],x[2],x[3])];
  }
  
  //-------------- correlators ------------
  int tmp[4],local[4],y[4];
  int *dist_count;
  Float **wline_pair;
  
  wline_pair = (Float**) smalloc(n*sizeof(Float*));
  for (int i=0; i<n; i++)
    wline_pair[i] = (Float*) smalloc(Max3*sizeof(Float));
  
  dist_count = (int*) smalloc(Max3*sizeof(int));
  
  for (int i=0; i<Max3; i++){
    dist_count[i] = 0;
    for (int j=0; j<n; j++){
      wline_pair[j][i] = 0.0;
    }
  }
  
  //************************************************************************************************
  // -------Method 1-------------------
  //************************************************************************************************
  
  
  // //   for (int ii=0; ii<n; ii++){
  
  // //     DIR = 2*s[0];
  // //     memcpy(u[0],v[ii],volume*sizeof(Matrix));
  // //     memcpy(v2,v[ii],volume*sizeof(Matrix));
  
  // //     for (int i=1; i<Max; i++){
  // //       pt.shift_field(&v2,(const int*)&DIR,1,1,&(u[i]));
  // //       memcpy(v2,u[i],volume*sizeof(Matrix));
  // //     }
  
  // //     for (int i=0; i<Max3; i++)
  // //       dist_count[i]=0;
  
  //   for (x[s[0]]=0; x[s[0]]<Max; x[s[0]]++)
  //     for (x[s[1]]=0; x[s[1]]<Max; x[s[1]]++)
  //       for (x[s[2]]=0; x[s[2]]<Max; x[s[2]]++){
  
  // 	tmp[0] = x[s[0]];
  // 	tmp[1] = x[s[1]];
  // 	tmp[2] = x[s[2]];
  
  // 	if ((tmp[0]+tmp[1]+tmp[2])!=0){
  // 	  int change;
  // 	  for (int i=0; i<2;i++)
  // 	    for (int j=i+1; j<3; j++)
  // 	      if (tmp[i]<tmp[j]) {
  // 		change = tmp[i];
  // 		tmp[i] = tmp[j];
  // 		tmp[j] = change;
  // 	      }
  
  // 	  dist_count[tmp[0]*Max2+tmp[1]*Max+tmp[2]]++;
  
  // 	  for (ii=0; ii<n; ii++){
  // 	    //memcpy(v2,v[ii],volume*sizeof(Matrix));
  // 	    for (y[0]=0; y[0]<node_sites[0]; y[0]++)
  // 	      for (y[1]=0; y[1]<node_sites[1]; y[1]++)
  // 		for (y[2]=0; y[2]<node_sites[2]; y[2]++)
  // 		  for (y[3]=0; y[3]<node_sites[3]; y[3]++)
  // 		    v2[INDp(y[0],y[1],y[2],y[3])]=v[ii][INDp(y[0],y[1],y[2],y[3])];
  // 	    //v2[INDp(y[0],y[1],y[2],y[3])]=u[x[s[0]]][INDp(y[0],y[1],y[2],y[3])];
  
  // 	    for (int i=0; i<3; i++){
  // 	    //for (int i=1; i<3; i++){
  // 	      DIR = 2*s[i];
  // 	      for (int j=0; j<x[s[i]]; j++){
  // 		pt.shift_field(&v2,(const int*)&DIR,1,1,&v1);
  
  // 		//memcpy(v2,v1,volume*sizeof(Matrix));
  // 		for (y[0]=0; y[0]<node_sites[0]; y[0]++)
  // 		  for (y[1]=0; y[1]<node_sites[1]; y[1]++)
  // 		    for (y[2]=0; y[2]<node_sites[2]; y[2]++)
  // 		      for (y[3]=0; y[3]<node_sites[3]; y[3]++)
  // 			v2[INDp(y[0],y[1],y[2],y[3])]=v1[INDp(y[0],y[1],y[2],y[3])];
  
  
  // 	      }
  // 	    }
  
  // 	    for (y[0]=0; y[0]<node_sites[0]; y[0]++)
  // 	      for (y[1]=0; y[1]<node_sites[1]; y[1]++)  
  // 		for (y[2]=0; y[2]<node_sites[2]; y[2]++)
  // 		  for (y[3]=0; y[3]<node_sites[3]; y[3]++){
  // 		    Matrix AAA;
  // 		    AAA.Dagger((const Float*)&(v1[INDp(y[0],y[1],y[2],y[3])]));
  // 		    Matrix BBB;
  // 		    BBB.DotMEqual(v[ii][INDp(y[0],y[1],y[2],y[3])],AAA);
  // 		    wline_pair[ii][tmp[0]*Max2+tmp[1]*Max+tmp[2]] += (BBB.Char3()).real()/volume;
  // 		  }
  // 	  }
  //       	}
  //       }
  
  //   printf("DONE_______________________DONE______________________________\n");
  
  //**************************************************************************************************
  //   --------Method 2 -------------------------------
  //**************************************************************************************************
  
  int nn[4];
  
  const int TransMax = 4096;
  int trans_N = volume*sizeof(Matrix)/sizeof(Float)/TransMax;
  int trans_res = volume*sizeof(Matrix)/sizeof(Float)%TransMax;

  printf("??????????????????????????????????\n");

  for (int i=0; i<n; i++){
    //memcpy(v1[i],v[i],volume*sizeof(Matrix));
    //memcpy(v4[i],v[i],volume*sizeof(Matrix));
    for (x[0]=0; x[0]<node_sites[0]; x[0]++)
      for (x[1]=0; x[1]<node_sites[1]; x[1]++)
	for (x[2]=0; x[2]<node_sites[2]; x[2]++)
	  for (x[3]=0; x[3]<node_sites[3]; x[3]++){
	    v1[i][INDp(x[0],x[1],x[2],x[3])]=v[i][INDp(x[0],x[1],x[2],x[3])];
	    //v2[i][INDp(x[0],x[1],x[2],x[3])]=v[i][INDp(x[0],x[1],x[2],x[3])];
	    //v3[i][INDp(x[0],x[1],x[2],x[3])]=v[i][INDp(x[0],x[1],x[2],x[3])];
	    v4[i][INDp(x[0],x[1],x[2],x[3])]=v[i][INDp(x[0],x[1],x[2],x[3])];
	  }
  }
  
  for (nn[0]=0; nn[0]<num_nodes[s[0]]/2+1;nn[0]++){ 
    
    for (int i=0; i<n; i++)
      //memcpy(v3[i],v4[i],volume*sizeof(Matrix));
      for (x[0]=0; x[0]<node_sites[0]; x[0]++)
	for (x[1]=0; x[1]<node_sites[1]; x[1]++)
	  for (x[2]=0; x[2]<node_sites[2]; x[2]++)
	    for (x[3]=0; x[3]<node_sites[3]; x[3]++)
	      v3[i][INDp(x[0],x[1],x[2],x[3])]=v4[i][INDp(x[0],x[1],x[2],x[3])];
    
    
    for (nn[1]=0; nn[1]<num_nodes[s[1]]/2+1; nn[1]++){
      
      for (int i=0; i<n; i++)
	//memcpy(v2[i],v3[i],volume*sizeof(Matrix));
	for (x[0]=0; x[0]<node_sites[0]; x[0]++)
	  for (x[1]=0; x[1]<node_sites[1]; x[1]++)
	    for (x[2]=0; x[2]<node_sites[2]; x[2]++)
	      for (x[3]=0; x[3]<node_sites[3]; x[3]++)
		v2[i][INDp(x[0],x[1],x[2],x[3])]=v3[i][INDp(x[0],x[1],x[2],x[3])];
      
      
      for (nn[2]=0; nn[2]<num_nodes[s[2]]/2+1; nn[2]++){
	
	int Dispmnt[3];
	Dispmnt[0] = nn[0]*node_sites[s[0]];
	Dispmnt[1] = nn[1]*node_sites[s[1]];
	Dispmnt[2] = nn[2]*node_sites[s[2]];
	//printf(" X Y Z ==== %d %d %d \n",Dispmnt[0],Dispmnt[1],Dispmnt[2]);

	for(x[s[0]]=0; x[s[0]]<node_sites[s[0]]; x[s[0]]++)
	  for(x[s[1]]=0; x[s[1]]<node_sites[s[1]]; x[s[1]]++)
	    for(x[s[2]]=0; x[s[2]]<node_sites[s[2]]; x[s[2]]++)
	      for(y[s[0]]=0; y[s[0]]<node_sites[s[0]]; y[s[0]]++)
		for(y[s[1]]=0; y[s[1]]<node_sites[s[1]]; y[s[1]]++)
		  for(y[s[2]]=0; y[s[2]]<node_sites[s[2]]; y[s[2]]++){
	
		    tmp[0] = y[s[0]]-x[s[0]]+Dispmnt[0];
		    tmp[1] = y[s[1]]-x[s[1]]+Dispmnt[1];
		    tmp[2] = y[s[2]]-x[s[2]]+Dispmnt[2];
		    int sum = tmp[0]+tmp[1]+tmp[2];
		    
		    if ((tmp[0]>=0)&&(tmp[1]>=0)&&(tmp[2]>=0)&&(sum>0)&&(tmp[0]<Max)&&(tmp[1]<Max)&&(tmp[2]<Max)) {
		      int change;
		      for (int i=0; i<2;i++)
			for (int j=i+1; j<3; j++)
			  if (tmp[i]<tmp[j]) {
			    change = tmp[i];
			    tmp[i] = tmp[j];
			    tmp[j] = change;
			  }
		      
		      dist_count[tmp[0]*Max2+tmp[1]*Max+tmp[2]]++;
		      
		      for (int i=0; i<n; i++){
			for (y[dir]=0; y[dir]<node_sites[dir]; y[dir]++){
			  x[dir]=y[dir];
			  Matrix AAA;
			  AAA.Dagger((const Float*)&(v2[i][INDp(y[0],y[1],y[2],y[3])]));
			  Matrix BBB;
			  BBB.DotMEqual(v[i][INDp(x[0],x[1],x[2],x[3])],AAA);
			  wline_pair[i][tmp[0]*Max2+tmp[1]*Max+tmp[2]] += (BBB.Char3()).real()/node_sites[dir];
			}
		      }
		      
		    }
		  }
	
	for (int i=0; i<n; i++){
	  for (int pp=0; pp<trans_N; pp++)
	    getPlusData((IFloat *)((IFloat*)v1[i]+pp*TransMax), (IFloat *)((IFloat*)v2[i]+pp*TransMax),TransMax, s[2]) ;
	  
	  if (trans_res>0)
	    getPlusData((IFloat *)((IFloat*)v1[i]+trans_N*TransMax), (IFloat *)((IFloat*)v2[i]+trans_N*TransMax),
			volume*sizeof(Matrix)/sizeof(Float)-trans_N*TransMax, s[2]) ;
	}
	
	for (int i=0; i<n; i++)
	  //memcpy(v2[i],v1[i],volume*sizeof(Matrix));
	  for (x[0]=0; x[0]<node_sites[0]; x[0]++)
	    for (x[1]=0; x[1]<node_sites[1]; x[1]++)
	      for (x[2]=0; x[2]<node_sites[2]; x[2]++)
		for (x[3]=0; x[3]<node_sites[3]; x[3]++)
		  v2[i][INDp(x[0],x[1],x[2],x[3])]=v1[i][INDp(x[0],x[1],x[2],x[3])];
      }
      
      for (int i=0; i<n; i++){
	for (int pp=0; pp<trans_N; pp++){
	  getPlusData((IFloat *)((IFloat*)v1[i]+pp*TransMax), (IFloat *)((IFloat*)v3[i]+pp*TransMax),TransMax, s[1]) ;
	}
	if (trans_res>0)
	  getPlusData((IFloat *)((IFloat*)v1[i]+trans_N*TransMax), (IFloat *)((IFloat*)v3[i]+trans_N*TransMax),
		      volume*sizeof(Matrix)/sizeof(Float)-trans_N*TransMax, s[1]) ;
      }
      
      for (int i=0; i<n; i++)
	//memcpy(v3[i],v1[i],volume*sizeof(Matrix));
	for (x[0]=0; x[0]<node_sites[0]; x[0]++)
	  for (x[1]=0; x[1]<node_sites[1]; x[1]++)
	    for (x[2]=0; x[2]<node_sites[2]; x[2]++)
	      for (x[3]=0; x[3]<node_sites[3]; x[3]++)
		v3[i][INDp(x[0],x[1],x[2],x[3])]=v1[i][INDp(x[0],x[1],x[2],x[3])];
    }
    
    for (int i=0; i<n; i++){
      for (int pp=0; pp<trans_N; pp++){
	getPlusData((IFloat *)((IFloat*)v1[i]+pp*TransMax), (IFloat *)((IFloat*)v4[i]+pp*TransMax),TransMax, s[0]) ;
      }
      if (trans_res>0)
	getPlusData((IFloat *)((IFloat*)v1[i]+trans_N*TransMax), (IFloat *)((IFloat*)v4[i]+trans_N*TransMax),
		    volume*sizeof(Matrix)/sizeof(Float)-trans_N*TransMax, s[0]) ;
    }
    
    
    for (int i=0; i<n; i++)
      //memcpy(v4[i],v1[i],volume*sizeof(Matrix));
      for (x[0]=0; x[0]<node_sites[0]; x[0]++)
	for (x[1]=0; x[1]<node_sites[1]; x[1]++)
	  for (x[2]=0; x[2]<node_sites[2]; x[2]++)
	    for (x[3]=0; x[3]<node_sites[3]; x[3]++)
	      v4[i][INDp(x[0],x[1],x[2],x[3])]=v1[i][INDp(x[0],x[1],x[2],x[3])];
  }
  
  printf("??????????????????????????????????\n");
  //-------------global sum and average--------------------------------------------------
  for (int k=0; k<n; k++)
    for (int i=0; i<Max3;i++)
      glb_sum((Float *)(wline_pair[k]+i));
  
  
  for (int k=0; k<n; k++)
    for (int i=0; i<Max3; i++)
      wline_pair[k][i] /= (Float)num_nodes[0]*num_nodes[1]*num_nodes[2]*num_nodes[3];
  
  //---------------print out -------------------------------------------------------------
  int ind = 0;
  int* d2;
  Float** crltor;
  
  
  crltor = (Float**) smalloc(n*sizeof(Float*));
  for (int k=0; k<n; k++)
    crltor[k] = (Float*) smalloc(Max3*sizeof(Float));
  d2 = (int*) smalloc(Max3*sizeof(int));
  
  FILE *fp;
  fp = Fopen(common_arg->filename,"w");
  
  VRB.Debug("begin to print out the results::::::::::::\n");
  printf("print out the results---------\n");
  
  for (x[0]=1;x[0]<Max;x[0]++)
    for (x[1]=0;x[1]<x[0]+1;x[1]++)
      for (x[2]=0;x[2]<x[1]+1;x[2]++)
	{
	  int pos = x[0]*Max2+x[1]*Max+x[2];
	  d2[ind] = x[0]*x[0]+x[1]*x[1]+x[2]*x[2];
	  for (int k=0; k<n; k++)
	    crltor[k][ind] = wline_pair[k][pos]/dist_count[pos];
	  
	  Fprintf(fp, "%0.16e  %d ", sqrt((Float)d2[ind]),
		  dist_count[pos]);
	  
	  for (int k=0; k<n; k++)
	    Fprintf(fp, "   %0.16e   ", crltor[k][ind]);
	  
	  Fprintf(fp, "\n");
	  
	  ind++;
	}

  printf("finished printing to file------------\n");
  
  Fclose(fp);

  Float af_cal_time = dclock();
  
  printf("wilson correlator takes %f minutes\n",(double)(af_cal_time-af_gf_time)/60);
  
    
  //----------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------

  //   for (int i=0; i<Max; i++)
  //     sfree(u[i]);
  
  //   sfree(u);

  for (int i=0; i<n; i++){
    sfree(v1[i]);
    sfree(v2[i]);
    sfree(v3[i]);
    sfree(v4[i]);
    sfree(v[i]);
  }
  
  sfree(v1);
  sfree(v2);
  sfree(v3);
  sfree(v4);
  sfree(v);
  
  sfree(d2);

  for (ii=0; ii<n; ii++)
    sfree(crltor[ii]);
  
  sfree(crltor);

  for (ii=0; ii<n; ii++)
    sfree(wline_pair[ii]);
  
  sfree(wline_pair);
  sfree(dist_count);
  
}

CPS_END_NAMESPACE
