void invcg_xp_update(double *out1, double *out2, double *A, double *B, double *mult, double *add, int size)
{
	//out1(x part)=A*mult+add(x part);
	//out2=B*mult+add(r part);
	double *x=out1;
	double *p=out2;
	double *p_pre=mult;
	double *x_pre=add;
	double *r=add+12;
	int i,j;
	//x=a*p_pre+x_pre;
	//p=b*p_pre+r;
	for(i=0;i<size;i++)
	{
		{
			*(x++) = (*A)*(*p_pre) + (*(x_pre++));
			*(p++) = (*B)*(*p_pre) + (*(r++));
			*p_pre++;
			*(x++) = (*A)*(*p_pre) + (*(x_pre++));
			*(p++) = (*B)*(*p_pre) + (*(r++));
			*p_pre++;
			*(x++) = (*A)*(*p_pre) + (*(x_pre++));
			*(p++) = (*B)*(*p_pre) + (*(r++));
			*p_pre++;
			*(x++) = (*A)*(*p_pre) + (*(x_pre++));
			*(p++) = (*B)*(*p_pre) + (*(r++));
			*p_pre++;
			*(x++) = (*A)*(*p_pre) + (*(x_pre++));
			*(p++) = (*B)*(*p_pre) + (*(r++));
			*p_pre++;
			*(x++) = (*A)*(*p_pre) + (*(x_pre++));
			*(p++) = (*B)*(*p_pre) + (*(r++));
			*p_pre++;
			*(x++) = (*A)*(*p_pre) + (*(x_pre++));
			*(p++) = (*B)*(*p_pre) + (*(r++));
			*p_pre++;
			*(x++) = (*A)*(*p_pre) + (*(x_pre++));
			*(p++) = (*B)*(*p_pre) + (*(r++));
			*p_pre++;
			*(x++) = (*A)*(*p_pre) + (*(x_pre++));
			*(p++) = (*B)*(*p_pre) + (*(r++));
			*p_pre++;
			*(x++) = (*A)*(*p_pre) + (*(x_pre++));
			*(p++) = (*B)*(*p_pre) + (*(r++));
			*p_pre++;
			*(x++) = (*A)*(*p_pre) + (*(x_pre++));
			*(p++) = (*B)*(*p_pre) + (*(r++));
			*p_pre++;
			*(x++) = (*A)*(*p_pre) + (*(x_pre++));
			*(p++) = (*B)*(*p_pre) + (*(r++));
			*p_pre++;
		}
		x+=12;
		x_pre+=12;
		r+=12;
	}

}
