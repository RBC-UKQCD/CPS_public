
//! Access to the elements of the \e SU(3) matrix
/*!
  Gets the element of the \e SU(3) matrix \e u with row \e row,
  column \e col and complex component \e d
*/
#define U(r,row,col,d)  *(u+(r+2*(row+3*(col+3*d))))
//! Access to the elements of a spinor vector.
/*!
  Gets the element of the spinor \e psi with spin \e s,
  colour \e c and complex component \e r
*/
#define PSI(r,c,s)      *(psi +(r+2*(c+3*s)))

//! As above, but the vector is called chi
#define CHI(r,c,s)      *(chi +(r+2*(c+3*s)))
#define TMP(r,c,s)      *(tmp +(r+2*(c+3*s))) 
#define TMP1(r,c,s)     *(tmp1+(r+2*(c+3*s))) 
#define TMP2(r,c,s)     *(tmp2+(r+2*(c+3*s))) 
#define TMP3(r,c,s)     *(tmp3+(r+2*(c+3*s))) 
#define TMP4(r,c,s)     *(tmp4+(r+2*(c+3*s))) 
#define TMP5(r,c,s)     *(tmp5+(r+2*(c+3*s))) 
#define TMP6(r,c,s)     *(tmp6+(r+2*(c+3*s))) 
#define TMP7(r,c,s)     *(tmp7+(r+2*(c+3*s))) 
#define TMP8(r,c,s)     *(tmp8+(r+2*(c+3*s))) 
#define FBUF(r,c,s)     *(fbuf+(r+2*(c+3*s))) 

#if 1
static inline void PLUSX( Float *u, Float *tmp5,int sdag, Float *psi){
	int c,s,mu;
	Float tmp[SPINOR_SIZE];
				for(c=0;c<3;c++){
					TMP(0,c,0) = PSI(0,c,0) + sdag * ( -PSI(1,c,3) ); 
					TMP(1,c,0) = PSI(1,c,0) + sdag * (  PSI(0,c,3) ); 
		
					TMP(0,c,1) = PSI(0,c,1) + sdag * ( -PSI(1,c,2) ); 
					TMP(1,c,1) = PSI(1,c,1) + sdag * (  PSI(0,c,2) ); 
		
					TMP(0,c,2) = PSI(0,c,2) + sdag * (  PSI(1,c,1) ); 
					TMP(1,c,2) = PSI(1,c,2) + sdag * ( -PSI(0,c,1) ); 
		
					TMP(0,c,3) = PSI(0,c,3) + sdag * (  PSI(1,c,0) ); 
					TMP(1,c,3) = PSI(1,c,3) + sdag * ( -PSI(0,c,0) ); 
				}
				/* multiply by U_mu */
				mu = 0;
				for(s=0;s<4;s++){
					for(c=0;c<3;c++){
						TMP5(0,c,s) = (  U(0,0,c,mu) * TMP(0,0,s)
								+ U(0,1,c,mu) * TMP(0,1,s)
								+ U(0,2,c,mu) * TMP(0,2,s) 
								+ U(1,0,c,mu) * TMP(1,0,s)
								+ U(1,1,c,mu) * TMP(1,1,s)
								+ U(1,2,c,mu) * TMP(1,2,s) );
						TMP5(1,c,s) = (  U(0,0,c,mu) * TMP(1,0,s)
								+ U(0,1,c,mu) * TMP(1,1,s)
								+ U(0,2,c,mu) * TMP(1,2,s) 
								- U(1,0,c,mu) * TMP(0,0,s)
								- U(1,1,c,mu) * TMP(0,1,s)
								- U(1,2,c,mu) * TMP(0,2,s) );
					}
				}
}
	
static inline void PLUSY( Float *u, Float *tmp6,int sdag, Float *psi){
	int c,s,mu;
	Float tmp[SPINOR_SIZE];
			mu = 1;
				for(c=0;c<3;c++){
					TMP(0,c,0) = PSI(0,c,0) + sdag * ( -PSI(0,c,3) ); 
					TMP(1,c,0) = PSI(1,c,0) + sdag * ( -PSI(1,c,3) ); 
		
					TMP(0,c,1) = PSI(0,c,1) + sdag * (  PSI(0,c,2) ); 
					TMP(1,c,1) = PSI(1,c,1) + sdag * (  PSI(1,c,2) ); 
		
					TMP(0,c,2) = PSI(0,c,2) + sdag * (  PSI(0,c,1) ); 
					TMP(1,c,2) = PSI(1,c,2) + sdag * (  PSI(1,c,1) ); 
		
					TMP(0,c,3) = PSI(0,c,3) + sdag * ( -PSI(0,c,0) ); 
					TMP(1,c,3) = PSI(1,c,3) + sdag * ( -PSI(1,c,0) ); 
				}
				/* multiply by U_mu */
				mu = 1;
				for(s=0;s<4;s++){
					for(c=0;c<3;c++){
						TMP6(0,c,s) = (  U(0,0,c,mu) * TMP(0,0,s)
								+ U(0,1,c,mu) * TMP(0,1,s)
								+ U(0,2,c,mu) * TMP(0,2,s) 
								+ U(1,0,c,mu) * TMP(1,0,s)
								+ U(1,1,c,mu) * TMP(1,1,s)
								+ U(1,2,c,mu) * TMP(1,2,s) );
						TMP6(1,c,s) = (  U(0,0,c,mu) * TMP(1,0,s)
								+ U(0,1,c,mu) * TMP(1,1,s)
								+ U(0,2,c,mu) * TMP(1,2,s) 
								- U(1,0,c,mu) * TMP(0,0,s)
								- U(1,1,c,mu) * TMP(0,1,s)
								- U(1,2,c,mu) * TMP(0,2,s) );
					}
				}
}
		
static inline void PLUSZ( Float *u, Float *tmp7,int sdag, Float *psi){
	int c,s,mu;
	Float tmp[SPINOR_SIZE];
	
			mu = 2;
				for(c=0;c<3;c++){
					TMP(0,c,0) = PSI(0,c,0) + sdag * ( -PSI(1,c,2) ); 
					TMP(1,c,0) = PSI(1,c,0) + sdag * (  PSI(0,c,2) ); 
		
					TMP(0,c,1) = PSI(0,c,1) + sdag * (  PSI(1,c,3) ); 
					TMP(1,c,1) = PSI(1,c,1) + sdag * ( -PSI(0,c,3) ); 
		
					TMP(0,c,2) = PSI(0,c,2) + sdag * (  PSI(1,c,0) ); 
					TMP(1,c,2) = PSI(1,c,2) + sdag * ( -PSI(0,c,0) ); 
		
					TMP(0,c,3) = PSI(0,c,3) + sdag * ( -PSI(1,c,1) ); 
					TMP(1,c,3) = PSI(1,c,3) + sdag * (  PSI(0,c,1) ); 
				}
				/* multiply by U_mu */
				for(s=0;s<4;s++){
					for(c=0;c<3;c++){
						TMP7(0,c,s) = (  U(0,0,c,mu) * TMP(0,0,s)
								+ U(0,1,c,mu) * TMP(0,1,s)
								+ U(0,2,c,mu) * TMP(0,2,s) 
								+ U(1,0,c,mu) * TMP(1,0,s)
								+ U(1,1,c,mu) * TMP(1,1,s)
								+ U(1,2,c,mu) * TMP(1,2,s) );
						TMP7(1,c,s) = (  U(0,0,c,mu) * TMP(1,0,s)
								+ U(0,1,c,mu) * TMP(1,1,s)
								+ U(0,2,c,mu) * TMP(1,2,s) 
								- U(1,0,c,mu) * TMP(0,0,s)
								- U(1,1,c,mu) * TMP(0,1,s)
								- U(1,2,c,mu) * TMP(0,2,s) );
					}
				}
}
	
static inline void PLUST( Float *u, Float *tmp8,int sdag, Float *psi){
	int c,s,mu;
	Float tmp[SPINOR_SIZE];
	
			mu = 3;
				for(c=0;c<3;c++){
					TMP(0,c,0) = PSI(0,c,0) + sdag * (  PSI(0,c,2) ); 
					TMP(1,c,0) = PSI(1,c,0) + sdag * (  PSI(1,c,2) ); 
		
					TMP(0,c,1) = PSI(0,c,1) + sdag * (  PSI(0,c,3) ); 
					TMP(1,c,1) = PSI(1,c,1) + sdag * (  PSI(1,c,3) ); 
		
					TMP(0,c,2) = PSI(0,c,2) + sdag * (  PSI(0,c,0) ); 
					TMP(1,c,2) = PSI(1,c,2) + sdag * (  PSI(1,c,0) ); 
		
					TMP(0,c,3) = PSI(0,c,3) + sdag * (  PSI(0,c,1) ); 
					TMP(1,c,3) = PSI(1,c,3) + sdag * (  PSI(1,c,1) ); 
				}
				/* multiply by U_mu */
				for(s=0;s<4;s++){
					for(c=0;c<3;c++){
						TMP8(0,c,s) = (  U(0,0,c,mu) * TMP(0,0,s)
								+ U(0,1,c,mu) * TMP(0,1,s)
								+ U(0,2,c,mu) * TMP(0,2,s) 
								+ U(1,0,c,mu) * TMP(1,0,s)
								+ U(1,1,c,mu) * TMP(1,1,s)
								+ U(1,2,c,mu) * TMP(1,2,s) );
						TMP8(1,c,s) = (  U(0,0,c,mu) * TMP(1,0,s)
								+ U(0,1,c,mu) * TMP(1,1,s)
								+ U(0,2,c,mu) * TMP(1,2,s) 
								- U(1,0,c,mu) * TMP(0,0,s)
								- U(1,1,c,mu) * TMP(0,1,s)
								- U(1,2,c,mu) * TMP(0,2,s) );
					}
				}
}

static inline void MINUSX( Float *u, Float *tmp, Float *tmp1,int sdag, Float *psi){
	int c,s,mu;
//	Float tmp[SPINOR_SIZE];
		mu = 0;
		for(c=0;c<3;c++){
			TMP(0,c,0) = PSI(0,c,0) - sdag * ( -PSI(1,c,3) ); 
			TMP(1,c,0) = PSI(1,c,0) - sdag * (  PSI(0,c,3) ); 

			TMP(0,c,1) = PSI(0,c,1) - sdag * ( -PSI(1,c,2) ); 
			TMP(1,c,1) = PSI(1,c,1) - sdag * (  PSI(0,c,2) ); 

			TMP(0,c,2) = PSI(0,c,2) - sdag * (  PSI(1,c,1) ); 
			TMP(1,c,2) = PSI(1,c,2) - sdag * ( -PSI(0,c,1) ); 

			TMP(0,c,3) = PSI(0,c,3) - sdag * (  PSI(1,c,0) ); 
			TMP(1,c,3) = PSI(1,c,3) - sdag * ( -PSI(0,c,0) ); 
		}
		/* multiply by U_mu */
		for(s=0;s<4;s++){
			for(c=0;c<3;c++){
				TMP1(0,c,s) = (  U(0,c,0,mu) * TMP(0,0,s)
						+ U(0,c,1,mu) * TMP(0,1,s)
						+ U(0,c,2,mu) * TMP(0,2,s) 
						- U(1,c,0,mu) * TMP(1,0,s)
						- U(1,c,1,mu) * TMP(1,1,s)
						- U(1,c,2,mu) * TMP(1,2,s) );
				TMP1(1,c,s) = (  U(0,c,0,mu) * TMP(1,0,s)
						+ U(0,c,1,mu) * TMP(1,1,s)
						+ U(0,c,2,mu) * TMP(1,2,s) 
						+ U(1,c,0,mu) * TMP(0,0,s)
						+ U(1,c,1,mu) * TMP(0,1,s)
						+ U(1,c,2,mu) * TMP(0,2,s) );
			}
		}
}


static inline void MINUSY( Float *u, Float *tmp, Float *tmp2,int sdag, Float *psi){
	int c,s,mu;
//	Float tmp[SPINOR_SIZE];
		mu = 1;
		for(c=0;c<3;c++){
			TMP(0,c,0) = PSI(0,c,0) - sdag * ( -PSI(0,c,3) ); 
			TMP(1,c,0) = PSI(1,c,0) - sdag * ( -PSI(1,c,3) ); 

			TMP(0,c,1) = PSI(0,c,1) - sdag * (  PSI(0,c,2) ); 
			TMP(1,c,1) = PSI(1,c,1) - sdag * (  PSI(1,c,2) ); 

			TMP(0,c,2) = PSI(0,c,2) - sdag * (  PSI(0,c,1) ); 
			TMP(1,c,2) = PSI(1,c,2) - sdag * (  PSI(1,c,1) ); 

			TMP(0,c,3) = PSI(0,c,3) - sdag * ( -PSI(0,c,0) ); 
			TMP(1,c,3) = PSI(1,c,3) - sdag * ( -PSI(1,c,0) ); 
		}
		/* multiply by U_mu */
		for(s=0;s<4;s++){
			for(c=0;c<3;c++){
				TMP2(0,c,s) = (  U(0,c,0,mu) * TMP(0,0,s)
						+ U(0,c,1,mu) * TMP(0,1,s)
						+ U(0,c,2,mu) * TMP(0,2,s) 
						- U(1,c,0,mu) * TMP(1,0,s)
						- U(1,c,1,mu) * TMP(1,1,s)
						- U(1,c,2,mu) * TMP(1,2,s) );
				TMP2(1,c,s) = (  U(0,c,0,mu) * TMP(1,0,s)
						+ U(0,c,1,mu) * TMP(1,1,s)
						+ U(0,c,2,mu) * TMP(1,2,s) 
						+ U(1,c,0,mu) * TMP(0,0,s)
						+ U(1,c,1,mu) * TMP(0,1,s)
						+ U(1,c,2,mu) * TMP(0,2,s) );
			}
		}
}

static inline void MINUSZ( Float *u, Float *tmp, Float *tmp3,int sdag, Float *psi){
	int c,s,mu;

		mu = 2;
		for(c=0;c<3;c++){
			TMP(0,c,0) = PSI(0,c,0) - sdag * ( -PSI(1,c,2) ); 
			TMP(1,c,0) = PSI(1,c,0) - sdag * (  PSI(0,c,2) ); 

			TMP(0,c,1) = PSI(0,c,1) - sdag * (  PSI(1,c,3) ); 
			TMP(1,c,1) = PSI(1,c,1) - sdag * ( -PSI(0,c,3) ); 

			TMP(0,c,2) = PSI(0,c,2) - sdag * (  PSI(1,c,0) ); 
			TMP(1,c,2) = PSI(1,c,2) - sdag * ( -PSI(0,c,0) ); 

			TMP(0,c,3) = PSI(0,c,3) - sdag * ( -PSI(1,c,1) ); 
			TMP(1,c,3) = PSI(1,c,3) - sdag * (  PSI(0,c,1) ); 
		}
		/* multiply by U_mu */
		for(s=0;s<4;s++){
			for(c=0;c<3;c++){
				TMP3(0,c,s) = (  U(0,c,0,mu) * TMP(0,0,s)
						+ U(0,c,1,mu) * TMP(0,1,s)
						+ U(0,c,2,mu) * TMP(0,2,s) 
						- U(1,c,0,mu) * TMP(1,0,s)
						- U(1,c,1,mu) * TMP(1,1,s)
						- U(1,c,2,mu) * TMP(1,2,s) );
				TMP3(1,c,s) = (  U(0,c,0,mu) * TMP(1,0,s)
						+ U(0,c,1,mu) * TMP(1,1,s)
						+ U(0,c,2,mu) * TMP(1,2,s) 
						+ U(1,c,0,mu) * TMP(0,0,s)
						+ U(1,c,1,mu) * TMP(0,1,s)
						+ U(1,c,2,mu) * TMP(0,2,s) );
			}
		}
}

static inline void MINUST( Float *u, Float *tmp, Float *tmp4,int sdag, Float *psi){
	int c,s,mu;
//	Float tmp[SPINOR_SIZE];

		mu = 3;
		for(c=0;c<3;c++){
			TMP(0,c,0) = PSI(0,c,0) - sdag * (  PSI(0,c,2) ); 
			TMP(1,c,0) = PSI(1,c,0) - sdag * (  PSI(1,c,2) ); 

			TMP(0,c,1) = PSI(0,c,1) - sdag * (  PSI(0,c,3) ); 
			TMP(1,c,1) = PSI(1,c,1) - sdag * (  PSI(1,c,3) ); 

			TMP(0,c,2) = PSI(0,c,2) - sdag * (  PSI(0,c,0) ); 
			TMP(1,c,2) = PSI(1,c,2) - sdag * (  PSI(1,c,0) ); 

			TMP(0,c,3) = PSI(0,c,3) - sdag * (  PSI(0,c,1) ); 
			TMP(1,c,3) = PSI(1,c,3) - sdag * (  PSI(1,c,1) ); 
		}
		/* multiply by U_mu */
		mu = 3;
		for(s=0;s<4;s++){
			for(c=0;c<3;c++){
				TMP4(0,c,s) = (  U(0,c,0,mu) * TMP(0,0,s)
						+ U(0,c,1,mu) * TMP(0,1,s)
						+ U(0,c,2,mu) * TMP(0,2,s) 
						- U(1,c,0,mu) * TMP(1,0,s)
						- U(1,c,1,mu) * TMP(1,1,s)
						- U(1,c,2,mu) * TMP(1,2,s) );
				TMP4(1,c,s) = (  U(0,c,0,mu) * TMP(1,0,s)
						+ U(0,c,1,mu) * TMP(1,1,s)
						+ U(0,c,2,mu) * TMP(1,2,s) 
						+ U(1,c,0,mu) * TMP(0,0,s)
						+ U(1,c,1,mu) * TMP(0,1,s)
						+ U(1,c,2,mu) * TMP(0,2,s) );
			}
		}
}
#endif

