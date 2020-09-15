#include<config.h>

#include<alg/alg_remez.h>
#include<util/command_line.h>
#include<util/qcdio.h>
#include<util/error.h>
#include<assert.h>
#include<vector>

USING_NAMESPACE_CPS

class ChevApprox{

private:
	double *ChevCoeffs=NULL;
	double high, low;
        int N;
public: 
      ChevApprox(int N_, double high_, double low_):
      N(N_),high(high_),low(low_),ChevCoeffs(NULL)
      { BuildApprox();}
      

double Chev (int n, double x){
  double Tp1, Tn = x, Tm1 = 1.;
  if(n==0) return Tm1;
//  else if(n==1) return Tn;
  else {
    for(int i=1;i<n; i++){
      Tp1 = 2.*x*Tn - Tm1;
      Tm1 = Tn;
      Tn = Tp1;
    }
  }
  return Tn;
}

double ChevRoots (int n, int k){
   assert(k<n);
   assert(k>=0);
   return cos(M_PI*((double)k+0.5)/(double)n);
}

double func(double x){
  return 1./x;
}

double unscaled(double theta){return low+(theta+1)*0.5*(high-low);}
double scaled(double x){return (x-low)/(high-low)*2.-1.; }

double BuildApprox (){
   if(ChevCoeffs) delete ChevCoeffs;
   ChevCoeffs = new double[N];
   for(int j=0;j<N;j++){
   ChevCoeffs[j]=0.;
   for(int k=0;k<N;k++){
      double theta = ChevRoots(N,k);
      double x = unscaled(theta);
      ChevCoeffs[j] += func(x)* Chev(j,theta);
//      printf("theta=%e func(%e)=%e Chev=%e\n",theta,x,func(x),Chev(j,theta));
   }
   if (j==0) ChevCoeffs[j] *= 1./double(N);
   else ChevCoeffs[j] *= 2./double(N);
   printf("ChevCoeffs[%d]=%e\n",j,ChevCoeffs[j]);
   }
}

double approx (double x){
   double theta = scaled(x);
   double answer=0;
   for(int k=0;k<N;k++){
	answer += ChevCoeffs[k]*Chev(k,theta);
   }
  return answer;
}

};

class NewtonApprox{

private:
	double *Coeffs=NULL;
//	double *Points=NULL;
        std::vector<double> Points;
	double high, low;
        int N;
public: 
      NewtonApprox(std::vector<double> points_)
      { 
      Points.resize(points_.size()); N = points_.size();
      for(int i =0;i<N;i++) Points[i] = points_[i]; 
      BuildApprox();}
      

double Newton (int n, double x){
  double Tp1, Tn = x, Tm1 = 1.;
  if(n==0) return Tm1;
//  else if(n==1) return Tn;
  else {
    for(int i=0;i<n; i++){
      Tn = Tm1*(x-Points[i]);
      Tm1 = Tn;
//      Tn = Tp1;
    }
  }
  return Tn;
}

double Roots (int n, int k){
   assert(k<n);
   assert(k>=0);
   return Points[k];
}

double func(double x){
  return 1./x;
}

double unscaled(double theta){return theta;}
double scaled(double x){return x; }

double BuildApprox (){
   if(Coeffs) delete Coeffs;
   Coeffs = new double[N];
   for(int j=0;j<N;j++){
      if (j==0) Coeffs[j]=func(Roots(N,0));
      else {
           double x= Roots(N,j);
           Coeffs[j]=(func(x)-approx(j-1,x))/Newton(j,x);
      }
   printf("Coeffs[%d]=%e\n",j,Coeffs[j]);
   }
}

double approx (int n, double x){
   double theta = scaled(x);
   double answer=0;
   for(int k=0;k<=n;k++){
	answer += Coeffs[k]*Newton(k,theta);
   }
  return answer;
}

};

#if 0
double  CgApprox (int N, double x,double scale){
//being really thick
double src = 1., sol=0.;
double res = src - x* sol;
double dir = res;
if (N==100){
#include "../talk/data/defl100/eig100/coeffs.h"
for(int i=0;i<b.size();i++){
  sol = a[i]*dir+sol;
  res = -a[i]*(x*scale)*dir+ res;
  dir = b[i]*dir+res;
}
} else if (N==300){
#include "../talk/data/defl300/eig100/coeffs.h"
for(int i=0;i<b.size();i++){
  sol = a[i]*dir+sol;
  res = -a[i]*(x*scale)*dir+ res;
  dir = b[i]*dir+res;
} 
}else
ERR.General("","CgApprox","Not implemented\n");

return sol*scale;
}
#endif

void use_omegas( const char *name, std::vector<std::complex<double> > &bs, std::vector<std::complex<double> > &omegas) 
{
  const char* fname = "use_omegas";
  VRB.Result("", fname, "%s: \n",name);

  assert(bs.size() == omegas.size());
  
  for(int s = 0; s < bs.size(); s++) {
    std::complex<double> b_s = 0.5 * (1.0 / omegas[s] + 1.0);
    std::complex<double> c_s = 0.5 * (1.0 / omegas[s] - 1.0);

    VRB.Result("", fname, "b[s=%d] = %0.10e + i %0.10e, c[s=%d] = %0.10e + i %0.10e\n", s, b_s.real(), b_s.imag(), s, c_s.real(), c_s.imag());

    bs[s] = b_s;
  }
}

std::complex<double> zmob(double x, std::vector<std::complex<double> > &bs)
{
	int N = bs.size();
	std::complex<double> num,den;
		num = 1. +(2. *bs[0]-1.)*x;
		den = 1. -(2. *bs[0]-1.)*x;
	for(int i =1;i<N;i++){
		num *= 1. +(2. *bs[i]-1.)*x;
		den *= 1. -(2. *bs[i]-1.)*x;
//		printf("zmob %d: bs=%g %g num den=%g %g\n",i,bs[i].real(),bs[i].imag(),std::abs(num),std::abs(den));
	}
//	if (std::abs(temp)<1.){
	if (1){
		return (num - den)/(num+den);
//	} else {
//		temp = 1./temp;
//		return (1. - temp)/(1.+temp);
	}

}

std::vector<std::complex<double> > omegas(14);
std::vector<std::complex<double> > bs14(14);
std::vector<std::complex<double> > omega12(12);
std::vector<std::complex<double> > bs12(12);
std::vector<std::complex<double> > bs10(10);
std::vector<std::complex<double> > bs8(8);
std::vector<std::complex<double> > omega8_12(8);
std::vector<std::complex<double> > bs8_12(8);
std::vector<std::complex<double> > omega10_12(10);
std::vector<std::complex<double> > bs10_12(10);
std::vector<std::complex<double> > mob12(12),mob12_32(12),mob24(24),mob32(32),shamir32(32),mob16(16);

double print_diff_all(double x,FILE *fp,FILE *fp2){

	std::complex<double>  m12=zmob(x,mob12);
	std::complex<double>  m24=zmob(x,mob24);
	std::complex<double>  m32=zmob(x,mob32);
	std::complex<double>  zmob14=zmob(x,bs14);
	std::complex<double>  zmob12=zmob(x,bs12);
	std::complex<double>  zmob10=zmob(x,bs10);
	std::complex<double>  zmob8=zmob(x,bs8);
	std::complex<double>  zmob8_12=zmob(x,bs8_12);
	std::complex<double>  zmob10_12=zmob(x,bs10_12);
	std::complex<double>  s32=zmob(x,shamir32);
	std::complex<double>  m12_32=zmob(x,mob12_32);
	std::complex<double>  m16=zmob(x,mob16);
    printf("x=%g mob32 =%g zmob14= %g mob24=%g zmob12= %g zmob10= %g zmob8=%g mob12=%g zmob8_12=%g zmob10_12=%g shamir32 = %g mob12_32 = %g m16 = %g \n",
	x, m32.real(), zmob14.real(),
        m24.real(), zmob12.real(), zmob10.real(), zmob8.real(),
	m12.real(),zmob8_12.real(),zmob10_12.real(),
	s32.real(),m12_32.real(),m16.real());
    if(fp) Fprintf(fp, "%0.16e\t%0.16e\t%0.16e\t%0.16e\t%0.16e\t%0.16\t%0.16e\t%0.16e\t%0.16e\t%0.16e\t%0.16e\t%0.16e\t%0.16e\t%0.16e\n",
	x,m32.real(),zmob14.real(),
	m24.real(),zmob12.real(),zmob10.real(),zmob8.real(),
	m12.real(),zmob8_12.real(),zmob10_12.real(),
	s32.real(),m12_32.real(),m16.real());
    if (fp2)Fprintf(fp2, "%0.16e\t%0.16e\t%0.16e\t%0.16e\t%0.16e\t%0.16e\t%0.16e\t%0.16e\t%0.16e\t%0.16e\t%0.16e\t%0.16e\t%0.16e\n",
	x,m32.real(),zmob14.real()-m32.real(),
	m24.real(),zmob12.real()-m24.real(),zmob10.real()-m24.real(),zmob8.real()-m24.real(),
	m12.real(),zmob8_12.real()-m12.real(),zmob10_12.real()-m12.real(),
	s32.real(),m12_32.real()-s32.real(),m16.real()-s32.real());
}

double print_diff(double x, 
	std::vector< std::vector<std::complex<double> > > mob,
	FILE *fp,FILE *fp2){

	std::vector< std::complex<double>  > results;
	std::cout <<"x= "<< x;  
	int size = mob.size();
	if(fp) Fprintf(fp,"%0.16e",x);
	if(fp2) Fprintf(fp2,"%0.16e",x);
	for(int i=0;i<size;i++){
		results.push_back(zmob(x,mob[i]));
		std::cout << " " << results[i] ;
	if(fp) Fprintf(fp,"\t%0.16e",results[i]);
	if(fp2) Fprintf(fp2,"\t%0.16e",results[i]-results[0]);
	}
	std::cout << std::endl;
	if(fp) Fprintf(fp,"\n");
	if(fp2) Fprintf(fp2,"\n");
}

#define CHEV
main(int argc, char *argv[])
{
	Start(&argc,&argv);
	CommandLine::is(argc,argv);
//  int N = CommandLine::arg_as_int();

    // concentrated
    omegas[0] = 1.395566319;
    omegas[1] = 0.8734555838;
    omegas[2] = 0.5976868047;
    omegas[3] = 0.3575610385;
    omegas[4] = 0.1765657746;
    omegas[5] = std::complex<double>(0.09298492219, +0.02404746984);
    omegas[6] = std::complex<double>(0.06049178877, +0.06547113295);
    omegas[7] = std::complex<double>(0.06049178877, -0.06547113295);
    omegas[8] = std::complex<double>(0.09298492219, -0.02404746984);
    omegas[9] = 0.119596635;
    omegas[10] = 0.2566858057;
    omegas[11] = 0.4745730595;
    omegas[12] = 0.6986654036;
    omegas[13] = 1.217162021;

	use_omegas("omega(14)", bs14,omegas);


    omega12[0]=1.375459773201338;
    omega12[1]=1.2029223219041347;
    omega12[2]=0.9477017499097049;
    omega12[3]=0.6985264047680767;
    omega12[4]=0.49568739129368533;
    omega12[5]=0.34412317778856183;
    omega12[6]=0.23503212589677097;
    omega12[7]=0.158126475421311;
    omega12[8]=std::complex<double>(0.11767285755366402,-0.027042564411322537);
    omega12[9]=std::complex<double>(0.11767285755366402,0.027042564411322537);
    omega12[10]=std::complex<double>(0.07803725152427325,-0.08118001690070928);
    omega12[11]=std::complex<double>(0.07803725152427325,0.08118001690070928);

	use_omegas("omega12",bs12,omega12);

    // concentrated
    bs10[0] = 8.4292038368159705e-01;
    bs10[1] = 9.2289979238280184e-01;
    bs10[2] = 1.1017200769981794e+00;
    bs10[3] = 1.4219097980542994e+00;
    bs10[4] = 1.9620523417564424e+00;
    bs10[5] = 2.8654191667525488e+00;
    bs10[6] = 4.4659153528626341e+00;
    bs10[7] = 5.5498080139636414e+00;
    bs10[8] = std::complex<double>(4.9320961582039766e+00, -3.5559998543638791e+00);
    bs10[9] = std::complex<double>(4.9320961582039766e+00, 3.5559998543638791e+00);

    // concentrated
    bs8[0] = 8.7405859026591415e-01;
    bs8[1] = 1.0190035946640621e+00;
    bs8[2] = 1.3658219772391325e+00;
    bs8[3] = 2.0514624685283716e+00;
    bs8[4] = 3.3661546566870935e+00;
    bs8[5] = 5.7994094935966318e+00;
    bs8[6] = std::complex<double>(6.7467786395784248e+00, -3.5543607236727506e+00);
    bs8[7] = std::complex<double>(6.7467786395784248e+00,  3.5543607236727506e+00);

    omega8_12[0]=1.4001203912480307;
    omega8_12[1]=1.1773978904661793;
    omega8_12[2]=0.7726975190839666;
    omega8_12[3]=0.4783206657578016;
    omega8_12[4]=0.283406239978277;
    omega8_12[5]=0.18359391336232933;
    omega8_12[6]=std::complex<double>(0.13830171268829916,-0.09133465746279675);
    omega8_12[7]=std::complex<double>(0.13830171268829916,0.09133465746279675);

	use_omegas("omega8_12",bs8_12,omega8_12);


    omega10_12[0]=std::complex<double>(1.1583495090856435,-0.29856482818559954);
    omega10_12[1]=0.9630215561588858;
    omega10_12[2]=0.4915802748744991;
    omega10_12[3]=std::complex<double>(0.23183097646302267,-0.03752095928313601);
    omega10_12[4]=std::complex<double>(0.17904585055176475,-0.1526089168867371);
    omega10_12[5]=std::complex<double>(0.17904585055176475,0.1526089168867371);
    omega10_12[6]=std::complex<double>(0.23183097646302267,0.03752095928313601);
    omega10_12[7]=0.33122433041804167;
    omega10_12[8]=0.7080735727535874;
    omega10_12[9]=std::complex<double>(1.1583495090856435,0.29856482818559954);

	use_omegas("omega10_12",bs10_12,omega10_12);

     
std::vector<std::complex<double> > omega10_24_4(10);
std::vector<std::complex<double> > bs10_24_4(10);
std::vector<std::complex<double> > mob24_4(24);

    omega10_24_4[0]=1.1060500168471108;
    omega10_24_4[1]=0.8628776580450203;
    omega10_24_4[2]=0.5714600222920464;
    omega10_24_4[3]=0.35186106038230286;
    omega10_24_4[4]=0.21025375146992942;
    omega10_24_4[5]=0.12380708297838017;
    omega10_24_4[6]=0.07154614533685795;
    omega10_24_4[7]=0.04627048410634415;
    omega10_24_4[8]=Complex(0.03263720419156931,-0.023047672895747478);
    omega10_24_4[9]=Complex(0.03263720419156931,0.023047672895747478);

	use_omegas("omega10_24_4",bs10_24_4,omega10_24_4);

std::vector<std::complex<double> > omega12_24_4(12);
std::vector<std::complex<double> > bs12_24_4(12);

    omega12_24_4[0]=1.0903256131299373;
    omega12_24_4[1]=0.9570283702230611;
    omega12_24_4[2]=0.7048886040934104;
    omega12_24_4[3]=0.48979921782791747;
    omega12_24_4[4]=0.328608311201356;
    omega12_24_4[5]=0.21664245377015995;
    omega12_24_4[6]=0.14121112711957107;
    omega12_24_4[7]=0.0907785101745156;
    omega12_24_4[8]=Complex(0.05608303440064219,-0.007537158177840385);
    omega12_24_4[9]=Complex(0.05608303440064219,0.007537158177840385);
    omega12_24_4[10]=Complex(0.0365221637144842,-0.03343945161367745);
    omega12_24_4[11]=Complex(0.0365221637144842,0.03343945161367745);

	use_omegas("omega12_24_4",bs12_24_4,omega12_24_4);


//  ChevApprox chev(N,high, low);
  for(int i=0;i<12;i++) mob12[i] = 1.5;
  for(int i=0;i<16;i++) mob16[i] = 1.5;
  for(int i=0;i<12;i++) mob12_32[i] = 1.8333333333333;
  for(int i=0;i<24;i++) mob24[i] = 1.5;
  for(int i=0;i<24;i++) mob24_4[i] = 2.5;
  for(int i=0;i<32;i++) mob32[i] = 1.5;
  for(int i=0;i<32;i++) shamir32[i] = 1.0;

  std::vector<std::vector<std::complex<double> > > list;
  list.push_back(mob24_4);
  list.push_back(bs10_24_4);
  list.push_back(bs12_24_4);

    FILE * fp = Fopen(CommandLine::arg(),"w");
    FILE * fp2 = Fopen(CommandLine::arg(),"w");
//    double scale = CommandLine::arg_as_Float();
#define LOOP
#define DELTA 1e-3
#define RANGE 3
	double x=0;
  while (x>-100){
	scanf("%lf",&x);
	print_diff_all(x,NULL,NULL);
	print_diff(x,list, NULL,NULL);
  }
    x = - RANGE ;
    while (x < RANGE ){
//	print_diff_all(x,NULL,NULL);
	print_diff(x,list,fp,fp2);
    x += DELTA ;
  }
  Fclose(fp);
  Fclose(fp2);
}
