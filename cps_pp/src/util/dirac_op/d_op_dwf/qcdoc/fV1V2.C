extern "C" {
void FtV1pV2Skip(double *_a, double b, const double *_c,
	const double *_d, int len )
{
	double a0;
	double a1;
	double a2;
	double a3;
	double a4;
	double a5;

	for(int i=0; i < len ; ++i) 
	{
		a0 = b * _c[0] + _d[0];
		a1 = b * _c[1] + _d[1];
		a2 = b * _c[2] + _d[2];
		a3 = b * _c[3] + _d[3];
		a4 = b * _c[4] + _d[4];
		a5 = b * _c[5] + _d[5];

		_a[0] = a0;
		_a[1] = a1;
		_a[2] = a2;
		_a[3] = a3;
		_a[4] = a4;
		_a[5] = a5;

		a0 = b * _c[6] + _d[6];
		a1 = b * _c[7] + _d[7];
		a2 = b * _c[8] + _d[8];
		a3 = b * _c[9] + _d[9];
		a4 = b * _c[10] + _d[10];
		a5 = b * _c[11] + _d[11];

		_a[6] = a0;
		_a[7] = a1;
		_a[8] = a2;
		_a[9] = a3;
		_a[10] = a4;
		_a[11] = a5;

		_a +=24;		
		_c +=24;		
		_d +=24;		
	}
}
}
