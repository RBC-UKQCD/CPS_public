void invcg_r_norm(double *resa, double *scale, double *mult, double *add, int ncvec, double *norm)
{
	//resa=scale*mult+add;
	//norm=resa*resa;
	*norm=0.0;
	double *py=resa;
	double *px=mult;
	double *pb=add;
	int i,j;
	//y=a*x+b
	for(i=0;i<ncvec;i++)
	{
		//call 12 times then skip 12 for py and pb
		{
			*py=(*scale)*(*(px++)) + (*(pb++));
			*norm+=(*py)*(*py);
			py++;
			*py=(*scale)*(*(px++)) + (*(pb++));
			*norm+=(*py)*(*py);
			py++;
			*py=(*scale)*(*(px++)) + (*(pb++));
			*norm+=(*py)*(*py);
			py++;
			*py=(*scale)*(*(px++)) + (*(pb++));
			*norm+=(*py)*(*py);
			py++;
			*py=(*scale)*(*(px++)) + (*(pb++));
			*norm+=(*py)*(*py);
			py++;
			*py=(*scale)*(*(px++)) + (*(pb++));
			*norm+=(*py)*(*py);
			py++;
			*py=(*scale)*(*(px++)) + (*(pb++));
			*norm+=(*py)*(*py);
			py++;
			*py=(*scale)*(*(px++)) + (*(pb++));
			*norm+=(*py)*(*py);
			py++;
			*py=(*scale)*(*(px++)) + (*(pb++));
			*norm+=(*py)*(*py);
			py++;
			*py=(*scale)*(*(px++)) + (*(pb++));
			*norm+=(*py)*(*py);
			py++;
			*py=(*scale)*(*(px++)) + (*(pb++));
			*norm+=(*py)*(*py);
			py++;
			*py=(*scale)*(*(px++)) + (*(pb++));
			*norm+=(*py)*(*py);
			py++;
		}

		py+=12;
		pb+=12;
	}
}
