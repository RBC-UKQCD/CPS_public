#include <vector>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
void bisec(std::vector<double> c,std::vector<double> b,int n, int m1,int m2,double eps1,double relfeh,std::vector<double> &x,int *z,double *eps2)
{
	double *wu;
	double h,q,x1,xu,x0,xmin,xmax; 
	int i,a,k;
	assert(c.size()>n);
	assert(b.size()>(n+1));

	b[1]=0.0;
	wu=(double *)malloc(sizeof(double)*(n+2));
	xmin=c[n]-fabs(b[n]);
	xmax=c[n]+fabs(b[n]);
	for(i=1;i<n;i++){
		h=fabs(b[i])+fabs(b[i+1]);
		if(c[i]+h>xmax) xmax= c[i]+h;
		if(c[i]-h<xmin) xmin= c[i]-h;
	}
	xmax *=2.;
//	printf("xmax=%e xmin=%e\n",xmax,xmin);
	*eps2=relfeh*((xmin+xmax)>0.0 ? xmax : -xmin);
	if(eps1<=0.0) eps1=*eps2;
	*eps2=0.5*eps1+7.0*(*eps2);
	x0=xmax;
	for(i=m1;i<=m2;i++){
		x[i]=xmax;
		wu[i]=xmin;
	}
	z=0;
	for(k=m2;k>=m1;k--){
		xu=xmin;
		i=k;
		do{
			if(xu<wu[i]){
				xu=wu[i];
				i=m1-1;
			}
			i--;
		}while(i>=m1);
		if(x0>x[k]) x0=x[k];
		while((x0-xu)>2*relfeh*(fabs(xu)+fabs(x0))+eps1){
			x1=(xu+x0)/2;
			z +=1;
			a=0;
			q=1.0;
			for(i=1;i<=n;i++){
			    q=c[i]-x1-((q!=0.0)? b[i]*b[i]/q:fabs(b[i])/relfeh);
				if(q<0) a++;
			}
//			printf("x1=%e a=%d\n",x1,a);
			if(a<k){
				if(a<m1){
				xu=x1;
				wu[m1]=x1;
				}else {
				xu=x1;
				wu[a+1]=x1;
				if(x[a]>x1) x[a]=x1;
				}
			}else x0=x1;
		}
	x[k]=(x0+xu)/2;
	}
    free(wu);
}
