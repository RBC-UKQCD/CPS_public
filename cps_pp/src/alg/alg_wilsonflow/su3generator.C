#include<config.h>
#include<math.h>
#include<util/vector.h>
//#include<iostream>

//using namespace cps;
CPS_START_NAMESPACE
//using namespace std;

extern const Float SU3_lambda[8][18]={
	{
		0.0,0.0,1.0,0.0,0.0,0.0,
		1.0,0.0,0.0,0.0,0.0,0.0,
		0.0,0.0,0.0,0.0,0.0,0.0
	},
	{
		0.0,0.0,0.0,-1.0,0.0,0.0,
		0.0,1.0,0.0,0.0,0.0,0.0,
		0.0,0.0,0.0,0.0,0.0,0.0
	},
	{
		1.0,0.0,0.0,0.0,0.0,0.0,
		0.0,0.0,-1.0,0.0,0.0,0.0,
		0.0,0.0,0.0,0.0,0.0,0.0
	},
	{
		0.0,0.0,0.0,0.0,1.0,0.0,
		0.0,0.0,0.0,0.0,0.0,0.0,
		1.0,0.0,0.0,0.0,0.0,0.0
	},
	{
		0.0,0.0,0.0,0.0,0.0,-1.0,
		0.0,0.0,0.0,0.0,0.0,0.0,
		0.0,1.0,0.0,0.0,0.0,0.0
	},
	{
		0.0,0.0,0.0,0.0,0.0,0.0,
		0.0,0.0,0.0,0.0,1.0,0.0,
		0.0,0.0,1.0,0.0,0.0,0.0
	},
	{
		0.0,0.0,0.0,0.0,0.0,0.0,
		0.0,0.0,0.0,0.0,0.0,-1.0,
		0.0,0.0,0.0,1.0,0.0,0.0
	},
	{
		1.0/sqrt(3.0),0.0,0.0,0.0,0.0,0.0,
		0.0,0.0,1.0/sqrt(3.0),0.0,0.0,0.0,
		0.0,0.0,0.0,0.0,-2.0/sqrt(3.0),0.0
	}
};

void generateSU3(Matrix & mat, Float Q[8])
{
	//return exp( i lambda_a/2 * Q[a] ), checked with mathematica.
	//follow the algorithm from physical review D 69,054501 (2004) paper
	Float invsqrt3 = 1.0/sqrt(3);
	Float c0 = (invsqrt3 * ((Q[0]*Q[0]*Q[7]+Q[1]*Q[1]*Q[7]+Q[2]*Q[2]*Q[7])*3.0-Q[7]*Q[7]*Q[7]) 
		- 0.5*invsqrt3 * (Q[3]*Q[3]*Q[7]+Q[4]*Q[4]*Q[7]+Q[5]*Q[5]*Q[7]+Q[6]*Q[6]*Q[7]) * 3.0
		+ 0.5 * (Q[0]*Q[3]*Q[5]+Q[0]*Q[4]*Q[6] - Q[1]*Q[3]*Q[6] + Q[1]*Q[4]*Q[5]) * 6.0
		+ 0.5 * (Q[2]*Q[3]*Q[3]+Q[2]*Q[4]*Q[4]-Q[2]*Q[5]*Q[5]-Q[2]*Q[6]*Q[6])*3.0 )/12.0;
	Float c1 = 0.25 * (Q[0]*Q[0]+Q[1]*Q[1]+Q[2]*Q[2]+Q[3]*Q[3]+Q[4]*Q[4]+Q[5]*Q[5]+Q[6]*Q[6]+Q[7]*Q[7]);
	if(c1<1e-20){mat=1.0;return;}
#if 0
	cout<<"c0 = "<<c0<<endl;
	cout<<"c1 = "<<c1<<endl;
#endif
	bool c0negative=false;
	if(c0<0){c0=-c0;c0negative=true;}

	Float c0max = 2*pow(c1/3.0, 1.5);
	Float theta = acos(c0/c0max);
	Float u = sqrt(c1/3.0)*cos(theta/3.0);
	Float w = sqrt(c1)*sin(theta/3.0);
	Float w2 = w*w;
	Float u2 = u*u;
	Float cosw = cos(w);
	Float cosu = cos(u);
	Float sinu = sin(u);
	Float cos2u = cos(2*u);
	Float sin2u = sin(2*u);

	//calculate sin(w)/w accurately. checked to double accuracy
	Float xi0w;
	if(fabs(w)>0.01)xi0w = sin(w)/w;
	else xi0w = 1.0 - w2/6.0*(1-w2/20.0*(1-w2/42.0));

	Rcomplex h0,h1,h2;
	h0 = (u2-w2)*Rcomplex(cos2u,sin2u) + Rcomplex(cosu,-sinu)*Rcomplex(8*u2*cosw,2*u*(3*u2+w2)*xi0w);
	h1 = 2*u*Rcomplex(cos2u,sin2u) - Rcomplex(cosu,-sinu)*Rcomplex(2*u*cosw, -(3*u2-w2)*xi0w);
	h2 = Rcomplex(cos2u,sin2u) - Rcomplex(cosu,-sinu)*Rcomplex(cosw,3*u*xi0w);

	Rcomplex f0,f1,f2;
	Float fdenominator = 9*u2-w2;
	f0 = h0/fdenominator;
	f1 = h1/fdenominator;
	f2 = h2/fdenominator;

	if(c0negative)
	{
		f0 = conj(f0);
		f1 = -conj(f1);
		f2 = conj(f2);
		c0 = -c0;
		c0negative = false;
	}
	
	mat = f0 + 2.0/3.0*c1*f2;
	Rcomplex ua[8];
	for(int i=0;i<8;i++)ua[i]=f1*Q[i];
	ua[0] += 0.5*f2*(invsqrt3*Q[0]*Q[7]+0.5*Q[3]*Q[5]+0.5*Q[4]*Q[6])*2.0;
	ua[1] += 0.5*f2*(invsqrt3*Q[1]*Q[7]-0.5*Q[3]*Q[6]+0.5*Q[4]*Q[5])*2.0;
	ua[2] += 0.5*f2*(invsqrt3*Q[2]*Q[7]*2.0+0.5*Q[3]*Q[3]+0.5*Q[4]*Q[4]-0.5*Q[5]*Q[5]-0.5*Q[6]*Q[6]);
	ua[3] += 0.5*f2*(-0.5*invsqrt3*Q[3]*Q[7]*2.0+0.5*Q[0]*Q[5]*2.0 - 0.5*Q[1]*Q[6]*2.0+0.5*Q[2]*Q[3]*2.0);
	ua[4] += 0.5*f2*(-0.5*invsqrt3*Q[4]*Q[7]*2.0+0.5*Q[0]*Q[6]*2.0 + 0.5*Q[1]*Q[5]*2.0 + 0.5*Q[2]*Q[4]*2.0);
	ua[5] += 0.5*f2*(-0.5*invsqrt3*Q[5]*Q[7]*2.0+0.5*Q[0]*Q[3]*2.0 + 0.5*Q[1]*Q[4]*2.0 - 0.5*Q[2]*Q[5]*2.0);
	ua[6] += 0.5*f2*(-0.5*invsqrt3*Q[6]*Q[7]*2.0+0.5*Q[0]*Q[4]*2.0 - 0.5*Q[1]*Q[3]*2.0 - 0.5*Q[2]*Q[6]*2.0);
	ua[7] += 0.5*f2*(invsqrt3*(Q[0]*Q[0]+Q[1]*Q[1]+Q[2]*Q[2]-Q[7]*Q[7])-0.5*invsqrt3*(Q[3]*Q[3]+Q[4]*Q[4]+Q[5]*Q[5]+Q[6]*Q[6]));

	Rcomplex tmp[9];
	Matrix *matpointer=(Matrix*)tmp;
	for(int i=0;i<8;i++){
		for(int j=0;j<9;j++)tmp[j]=0.5*ua[i]*((Rcomplex*)(SU3_lambda[i]))[j]; 
		mat += *matpointer;
	} 
}

CPS_END_NAMESPACE
