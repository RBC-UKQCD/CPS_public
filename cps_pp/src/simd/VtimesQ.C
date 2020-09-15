//#define xlc
void VtimesQ(double *QZ, int N, double **V, double *row_tmp_in, double *row_tmp_out, int row_start, int row_step, int row_end) 
{
#ifdef xlc
	vector4double v1;
	vector4double v2;
	vector4double v3;
#else
#endif
	for(int r = row_start; r <= row_end; r += row_step) {
		int m = r + row_step <= row_end ? row_step : row_end - r; // number of rows in current iteration
		
		// Copy V into row_tmp_in
		for (int c2 = 0; c2 < N; c2++) {
			int c2m = c2 * row_step;
			for (int i = 0; i < m; i++) {
				row_tmp_in[c2m+i] = V[c2][r+i];
			}
			for (int i = m; i < row_step; i++) {
				row_tmp_in[c2m+i] = 0;
			}
		}

		// Linear transformation: row_tmp_out = row_tmp_in * Q
		for(int c1 = 0; c1 < N; c1++) {
			int c1m = c1 * row_step;
			for (int i = 0; i < row_step; i++) {
				row_tmp_out[c1m+i] = 0;
			}
			for(int c2 = 0; c2 < N; c2++) { 
				int c2m = c2 * row_step;
#ifdef xlc
				v3 = vec_lds(0, QZ+c1*N+c2);  
#else
				double qze = QZ[c1*N+c2];
#endif
				for (int i = 0; i < row_step; i += 16) {
					double *tmp0 = row_tmp_in  + c2m + i;
					double *tmp1 = row_tmp_out + c1m + i;
					{
#ifdef xlc
						v1 = vec_lda(0, tmp1);  
						v2 = vec_lda(0, tmp0);  
						v1 = vec_madd(v3, v2, v1);    // axpy 
						vec_st(v1, 0, tmp1);    

						v1 = vec_lda(0, tmp1+4);  
						v2 = vec_lda(0, tmp0+4);  
						v1 = vec_madd(v3, v2, v1);    // axpy 
						vec_st(v1, 0, tmp1+4);    

						v1 = vec_lda(0, tmp1+8);  
						v2 = vec_lda(0, tmp0+8);  
						v1 = vec_madd(v3, v2, v1);    // axpy 
						vec_st(v1, 0, tmp1+8);    

						v1 = vec_lda(0, tmp1+12);  
						v2 = vec_lda(0, tmp0+12);  
						v1 = vec_madd(v3, v2, v1);    // axpy 
						vec_st(v1, 0, tmp1+12);    
#else
						tmp1[0] += tmp0[0] * qze;
						tmp1[1] += tmp0[1] * qze;
						tmp1[2] += tmp0[2] * qze;
						tmp1[3] += tmp0[3] * qze;
						tmp1[4] += tmp0[4] * qze;
						tmp1[5] += tmp0[5] * qze;
						tmp1[6] += tmp0[6] * qze;
						tmp1[7] += tmp0[7] * qze;
						tmp1[8] += tmp0[8] * qze;
						tmp1[9] += tmp0[9] * qze;
						tmp1[10] += tmp0[10] * qze;
						tmp1[11] += tmp0[11] * qze;
						tmp1[12] += tmp0[12] * qze;
						tmp1[13] += tmp0[13] * qze;
						tmp1[14] += tmp0[14] * qze;
						tmp1[15] += tmp0[15] * qze;
#endif
					}
				}
			}
		}

		// Copy row_tmp_out into V
		for(int c = 0; c < N; c++) {
			int cm = c * row_step;
			for (int i = 0; i < m; i++) {
				V[c][r+i] = row_tmp_out[cm+i];
			}
		}
	}
}

//#define xlc
//#ifdef xlc
//		vector4double v1;
//		vector4double v2;
//		vector4double v3;
//
//		for(long i = 0; i < task; i++) {
//			v1 = vec_lds(i*4, xv);  // Loads b[i], b[i+1], b[i+2], and b[i+3].
//			v2 = vec_lds(i*4, yv);  // Loads c[i], c[i+1], c[i+2], and c[i+3].
//			v3 = vec_lds(i*4, av);  // Loads c[i], c[i+1], c[i+2], and c[i+3].
//			v1 = vec_madd(v1, v2, v3);    // axpy 
//			vec_st(v3, i*4, yv);    // Stores the result to yv[i], yv[i+1], yv[i+2], and yv[i+3].
//		}
//#else
//		for(long i = 0; i < task; i ++) {
//			long off = i * 4;
//			yv[off    ] += av[off    ] * xv[off    ];
//			yv[off + 1] += av[off + 1] * xv[off + 1];
//			yv[off + 2] += av[off + 2] * xv[off + 2];
//			yv[off + 3] += av[off + 3] * xv[off + 3];
//		}
//#endif
//
void VtimesQ(double *QZ, int N, float** V, double *row_tmp_in, double *row_tmp_out, int row_start, int row_step, int row_end)
{
#ifdef xlc
	vector4double v1;
	vector4double v2;
	vector4double v3;
#else
#endif
	for(int r = row_start; r <= row_end; r += row_step) {
		int m = r + row_step <= row_end ? row_step : row_end - r; // number of rows in current iteration

		// Copy V into row_tmp_in
		for (int c2 = 0; c2 < N; c2++) {
			int c2m = c2 * row_step;
			for (int i = 0; i < m; i++) {
				row_tmp_in[c2m+i] = V[c2][r+i];
			}
			for (int i = m; i < row_step; i++) {
				row_tmp_in[c2m+i] = 0;
			}
		}

		// Linear transformation: row_tmp_out = row_tmp_in * Q
		for(int c1 = 0; c1 < N; c1++) {
			int c1m = c1 * row_step;
			for (int i = 0; i < row_step; i++) {
				row_tmp_out[c1m+i] = 0;
			}
			for(int c2 = 0; c2 < N; c2++) {
				int c2m = c2 * row_step;
#ifdef xlc
				v3 = vec_lds(0, QZ+c1*N+c2);
#else
				double qze = QZ[c1*N+c2];
#endif
				for (int i = 0; i < row_step; i += 16) {
					double *tmp0 = row_tmp_in  + c2m + i;
					double *tmp1 = row_tmp_out + c1m + i;
					{
#ifdef xlc
						v1 = vec_lda(0, tmp1);
						v2 = vec_lda(0, tmp0);
						v1 = vec_madd(v3, v2, v1);    // axpy
						vec_st(v1, 0, tmp1);

						v1 = vec_lda(0, tmp1+4);
						v2 = vec_lda(0, tmp0+4);
						v1 = vec_madd(v3, v2, v1);    // axpy
						vec_st(v1, 0, tmp1+4);

						v1 = vec_lda(0, tmp1+8);
						v2 = vec_lda(0, tmp0+8);
						v1 = vec_madd(v3, v2, v1);    // axpy
						vec_st(v1, 0, tmp1+8);

						v1 = vec_lda(0, tmp1+12);
						v2 = vec_lda(0, tmp0+12);
						v1 = vec_madd(v3, v2, v1);    // axpy
						vec_st(v1, 0, tmp1+12);
#else
						tmp1[0] += tmp0[0] * qze;
						tmp1[1] += tmp0[1] * qze;
						tmp1[2] += tmp0[2] * qze;
						tmp1[3] += tmp0[3] * qze;
						tmp1[4] += tmp0[4] * qze;
						tmp1[5] += tmp0[5] * qze;
						tmp1[6] += tmp0[6] * qze;
						tmp1[7] += tmp0[7] * qze;
						tmp1[8] += tmp0[8] * qze;
						tmp1[9] += tmp0[9] * qze;
						tmp1[10] += tmp0[10] * qze;
						tmp1[11] += tmp0[11] * qze;
						tmp1[12] += tmp0[12] * qze;
						tmp1[13] += tmp0[13] * qze;
						tmp1[14] += tmp0[14] * qze;
						tmp1[15] += tmp0[15] * qze;
#endif
					}
				}
			}
		}

		// Copy row_tmp_out into V
		for(int c = 0; c < N; c++) {
			int cm = c * row_step;
			for (int i = 0; i < m; i++) {
				V[c][r+i] = row_tmp_out[cm+i];
			}
		}
	}
}
