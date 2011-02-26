//89 seconds. 82Mflops
#ifndef BLOCK_SIZE
#define BLOCK_SIZE ((int) 4)
#endif
void basic_dgemm (const int KB, const int M, const int N, const int K, double **A, const double *B, double *C, const int NB)
{
	int i, j, k;
	register const double *Bp;
	register const double **pA=A;
	switch(K)
	{
		case 1:
			for (i = 0; i < M; ++i) 
			{
				Bp=B;
				for (j = 0; j < N; ++j) 
				{
					*(C + i*NB + j) += pA[0][i]*(*(Bp++));
				}
			}
			break;
		case 2:
			for (i = 0; i < M; ++i) 
			{
				Bp=B;
				for (j = 0; j < N; ++j) 
				{
					register double cij = *(C + i*NB + j);
					cij += pA[0][i]*(*(Bp++));
					cij += pA[1][i]*(*(Bp++));
					*(C + i*NB + j) = cij;
				}
			}
			break;
		case 3:
			for (i = 0; i < M; ++i) 
			{
				Bp=B;
				for (j = 0; j < N; ++j) 
				{
					register double cij = *(C + i*NB + j);
					cij += pA[0][i]*(*(Bp++));
					cij += pA[1][i]*(*(Bp++));
					cij += pA[3][i]*(*(Bp++));
					*(C + i*NB + j) = cij;
				}
			}
			break;
		case 4:
			for (i = 0; i < M; ++i) 
			{
				Bp=B;
				for (j = 0; j < N; ++j) 
				{
					register double cij = *(C + i*NB + j);
					cij += pA[0][i]*(*(Bp++));
					cij += pA[1][i]*(*(Bp++));
					cij += pA[2][i]*(*(Bp++));
					cij += pA[3][i]*(*(Bp++));
					*(C + i*NB + j) = cij;
				}
			}
			break;
		default:
			break;
	}
}

void unrolled_dgemm (const int lda, double **A, const double *B, double *C, const int N)
{
	const double * pB=B;
	register const double A0_0 = A[0][0];
	register const double B0_0 = *(pB++);
	register const double A0_1 = A[1][0];
	register const double B1_0 = *(pB++);
	register const double A0_2 = A[2][0];
	register const double B2_0 = *(pB++);
	register const double A0_3 = A[3][0];
	register const double B3_0 = *(pB++);
	register const double A1_0 = A[0][1];
	register const double B0_1 = *(pB++);
	register const double A1_1 = A[1][1];
	register const double B1_1 = *(pB++);
	register const double A1_2 = A[2][1];
	register const double B2_1 = *(pB++);
	register const double A1_3 = A[3][1];
	register const double B3_1 = *(pB++);
	register const double A2_0 = A[0][2];
	register const double B0_2 = *(pB++);
	register const double A2_1 = A[1][2];
	register const double B1_2 = *(pB++);
	register const double A2_2 = A[2][2];
	register const double B2_2 = *(pB++);
	register const double A2_3 = A[3][2];
	register const double B3_2 = *(pB++);
	register const double A3_0 = A[0][3];
	register const double B0_3 = *(pB++);
	register const double A3_1 = A[1][3];
	register const double B1_3 = *(pB++);
	register const double A3_2 = A[2][3];
	register const double B2_3 = *(pB++);
	register const double A3_3 = A[3][3];
	register const double B3_3 = *(pB++);

	register double c0_0 = *(C + 0*N + 0);
	register double c0_1 = *(C + 0*N + 1);
	register double c0_2 = *(C + 0*N + 2);
	register double c0_3 = *(C + 0*N + 3);
	register double c1_0 = *(C + 1*N + 0);
	register double c1_1 = *(C + 1*N + 1);
	register double c1_2 = *(C + 1*N + 2);
	register double c1_3 = *(C + 1*N + 3);
	register double c2_0 = *(C + 2*N + 0);
	register double c2_1 = *(C + 2*N + 1);
	register double c2_2 = *(C + 2*N + 2);
	register double c2_3 = *(C + 2*N + 3);
	register double c3_0 = *(C + 3*N + 0);
	register double c3_1 = *(C + 3*N + 1);
	register double c3_2 = *(C + 3*N + 2);
	register double c3_3 = *(C + 3*N + 3);

	c0_0 += A0_0 * B0_0;
	c0_0 += A0_1 * B1_0;
	c0_0 += A0_2 * B2_0;
	c0_0 += A0_3 * B3_0;
	*(C + 0*N + 0) = c0_0;

	c0_1 += A0_0 * B0_1;
	c0_1 += A0_1 * B1_1;
	c0_1 += A0_2 * B2_1;
	c0_1 += A0_3 * B3_1;
	*(C + 0*N + 1) = c0_1;

	c0_2 += A0_0 * B0_2;
	c0_2 += A0_1 * B1_2;
	c0_2 += A0_2 * B2_2;
	c0_2 += A0_3 * B3_2;
	*(C + 0*N + 2) = c0_2;

	c0_3 += A0_0 * B0_3;
	c0_3 += A0_1 * B1_3;
	c0_3 += A0_2 * B2_3;
	c0_3 += A0_3 * B3_3;
	*(C + 0*N + 3) = c0_3;

	c1_0 += A1_0 * B0_0;
	c1_0 += A1_1 * B1_0;
	c1_0 += A1_2 * B2_0;
	c1_0 += A1_3 * B3_0;
	*(C + 1*N + 0) = c1_0;

	c1_1 += A1_0 * B0_1;
	c1_1 += A1_1 * B1_1;
	c1_1 += A1_2 * B2_1;
	c1_1 += A1_3 * B3_1;
	*(C + 1*N + 1) = c1_1;

	c1_2 += A1_0 * B0_2;
	c1_2 += A1_1 * B1_2;
	c1_2 += A1_2 * B2_2;
	c1_2 += A1_3 * B3_2;
	*(C + 1*N + 2) = c1_2;

	c1_3 += A1_0 * B0_3;
	c1_3 += A1_1 * B1_3;
	c1_3 += A1_2 * B2_3;
	c1_3 += A1_3 * B3_3;
	*(C + 1*N + 3) = c1_3;

	c2_0 += A2_0 * B0_0;
	c2_0 += A2_1 * B1_0;
	c2_0 += A2_2 * B2_0;
	c2_0 += A2_3 * B3_0;
	*(C + 2*N + 0) = c2_0;

	c2_1 += A2_0 * B0_1;
	c2_1 += A2_1 * B1_1;
	c2_1 += A2_2 * B2_1;
	c2_1 += A2_3 * B3_1;
	*(C + 2*N + 1) = c2_1;

	c2_2 += A2_0 * B0_2;
	c2_2 += A2_1 * B1_2;
	c2_2 += A2_2 * B2_2;
	c2_2 += A2_3 * B3_2;
	*(C + 2*N + 2) = c2_2;

	c2_3 += A2_0 * B0_3;
	c2_3 += A2_1 * B1_3;
	c2_3 += A2_2 * B2_3;
	c2_3 += A2_3 * B3_3;
	*(C + 2*N + 3) = c2_3;

	c3_0 += A3_0 * B0_0;
	c3_0 += A3_1 * B1_0;
	c3_0 += A3_2 * B2_0;
	c3_0 += A3_3 * B3_0;
	*(C + 3*N + 0) = c3_0;

	c3_1 += A3_0 * B0_1;
	c3_1 += A3_1 * B1_1;
	c3_1 += A3_2 * B2_1;
	c3_1 += A3_3 * B3_1;
	*(C + 3*N + 1) = c3_1;

	c3_2 += A3_0 * B0_2;
	c3_2 += A3_1 * B1_2;
	c3_2 += A3_2 * B2_2;
	c3_2 += A3_3 * B3_2;
	*(C + 3*N + 2) = c3_2;

	c3_3 += A3_0 * B0_3;
	c3_3 += A3_1 * B1_3;
	c3_3 += A3_2 * B2_3;
	c3_3 += A3_3 * B3_3;
	*(C + 3*N + 3) = c3_3;
}
//A(M*K) B(K*N)= C(M*N)
void matrix_dgemm (const int M,const int N, const int K, double **A, const double *B, double *C)
{
	const int m_cblocks = M / BLOCK_SIZE + (M%BLOCK_SIZE ? 1 : 0);
	const int k_cblocks = K / BLOCK_SIZE + (K%BLOCK_SIZE ? 1 : 0);
	const int n_cblocks = N / BLOCK_SIZE + (N%BLOCK_SIZE ? 1 : 0);
	double *pt;
	int ci, cj, ck;
	for (ci = 0; ci < m_cblocks; ++ci) {
		const int i = ci * BLOCK_SIZE;
		const int M_cblock = (i+BLOCK_SIZE>M? M-i:BLOCK_SIZE);
		for (cj = 0; cj < n_cblocks; ++cj) {
			const int j = cj * BLOCK_SIZE;
			const int N_cblock = (j+BLOCK_SIZE>N? N-j:BLOCK_SIZE);
			for (ck = 0; ck < k_cblocks; ++ck) {
				const int k = ck * BLOCK_SIZE;
				const int K_cblock = (k+BLOCK_SIZE>K? K-k:BLOCK_SIZE);
				if((M_cblock==BLOCK_SIZE)&&(N_cblock==BLOCK_SIZE)&&(K_cblock==BLOCK_SIZE))
				{
					//unrolled_dgemm (K, A + i*K + k, B + k + j*K, C + i*N + j, N);
					unrolled_dgemm (K, A+k, B + k*BLOCK_SIZE+ j*K, C + i*N + j, N);
				}
				else
				{
					//basic_dgemm (K,M_cblock, N_cblock, K_cblock, A + i*K + k, B + k + j*K, C + i*N + j,N);
					basic_dgemm (K,M_cblock, N_cblock, K_cblock, A + k, B + k*N_cblock+ j*K, C + i*N + j,N);
				}
			}
		}
	//replace A with C
	pt=C+i*N;
	for(ck=0;ck<M_cblock;ck++)
	for(cj=0;cj<N;cj++)
			A[cj][ck]=(*(pt++));

	for(cj=0;cj<K;cj++)A[cj]+=M_cblock;
	}
}

