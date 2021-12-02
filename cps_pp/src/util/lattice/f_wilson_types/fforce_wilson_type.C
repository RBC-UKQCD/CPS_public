// -*- mode:c++;c-basic-offset:4 -*-

#include <util/gjp.h>
#include <util/error.h>
#include <util/vector.h>
#include <util/omp_wrapper.h>
#include <util/sproj_tr.h>
#include <comms/scu.h>
#include <comms/glb.h>
#include <alg/force_arg.h>
#include <util/wilson.h>
#include <util/lattice/fforce_wilson_type.h>
#include <util/lattice/bfm_evo_aux.h>


USING_NAMESPACE_CPS;

static void compute_coord(long x[4], const long hl[4], const long low[4], long i)
{
    x[0] = i % hl[0] + low[0]; i /= hl[0];
    x[1] = i % hl[1] + low[1]; i /= hl[1];
    x[2] = i % hl[2] + low[2]; i /= hl[2];
    x[3] = i % hl[3] + low[3];
}

static long idx_4d(const long x[4], const long lx[4])
{
    long ret = 0;
    for(int i = 3; i >= 0; --i) {
        ret = ret * lx[i] + x[i];
    }
    return ret;
}

static long idx_4d_surf(const long x[4], const long lx[4], int mu)
{
    long ret = 0;
    for(int i = 3; i >= 0; --i) {
        if(i == mu) continue;
        ret = ret * lx[i] + x[i];
    }
    return ret;
}

static void thread_work_partial(long nwork, int me, int nthreads,
                                long &mywork, long &myoff)
{
    long basework = nwork / nthreads;
    long backfill = nthreads - (nwork % nthreads);
    mywork = (nwork + me) / nthreads;
    myoff  = basework * me;
    if ( me > backfill ) 
        myoff += (me-backfill);
}

void updateForce(ForceArg *f_arg, const Matrix &m)
{
    Float a2 = m.norm();
    Float a = sqrt(a2);
    
    f_arg->L1 += a;
    f_arg->L2 += a2;
    f_arg->Linf = f_arg->Linf > a ? f_arg->Linf : a;
}

static Float vsum (Float *v){
    Float sum=0.;
#if 1
    for(int i =0;i<SPINOR_SIZE;i++) sum +=*(v+i)**(v+i);
#else
    sum = *v;
#endif
    return sum;
}
// Calculate fermion force on a specific site, also do the
// summation over s direction (if s > 1).
static void do_site_force(Matrix *force, const Matrix &gauge,
                          Float *v1, Float *v1p,
                          Float *v2, Float *v2p,
                          int mu, int ls)
{
    const char *fname="do_site_force()";
    Matrix t1, t2;
    
    Float *t1f = (Float *)&t1;
    Float *t2f = (Float *)&t2;
    VRB.Flow("",fname,"ls=%d\n",ls);
    VRB.Flow("",fname,"v1=%g\n",vsum(v1));
    VRB.Flow("",fname,"v1p=%g\n",vsum(v1p));
    VRB.Flow("",fname,"v2=%g\n",vsum(v2));
    VRB.Flow("",fname,"v2p=%g\n",vsum(v2p));
    

    // sproj_tr[mu](   (Float *)&t1, v1p, v2, ls, 0, 0);
    // sproj_tr[mu+4]( (Float *)&t2, v2p, v1, ls, 0, 0);
    
    switch(mu) {
    case 0:
      sprojTrXm(t1f, v1p, v2, ls, 0, 0);
      sprojTrXp(t2f, v2p, v1, ls, 0, 0);
      break;
    case 1:
      sprojTrYm(t1f, v1p, v2, ls, 0, 0);
      sprojTrYp(t2f, v2p, v1, ls, 0, 0);
      break;
    case 2:
      sprojTrZm(t1f, v1p, v2, ls, 0, 0);
      sprojTrZp(t2f, v2p, v1, ls, 0, 0);
      break;
    default:
      sprojTrTm(t1f, v1p, v2, ls, 0, 0);
      sprojTrTp(t2f, v2p, v1, ls, 0, 0);
    }

    VRB.Flow("",fname,"t1=%g\n",t1.Norm());
    VRB.Flow("",fname,"t2=%g\n",t2.Norm());
    t1 += t2;
    force->DotMEqual(gauge, t1);
    VRB.Debug("",fname,"force=%g\n",force->Norm());
}

//CK: Same as above, only it calculates the contribution of the second G-parity flavour
static void do_site_force_f1(Matrix *force, const Matrix &gauge,
                          Float *v1, Float *v1p,
                          Float *v2, Float *v2p,
                          int mu, int ls)
{
    Matrix t1, t2;
    
    Float *t1f = (Float *)&t1;
    Float *t2f = (Float *)&t2;

    // sproj_tr[mu](   (Float *)&t1, v1p, v2, ls, 0, 0);
    // sproj_tr[mu+4]( (Float *)&t2, v2p, v1, ls, 0, 0);
    
    switch(mu) {
    case 0:
      sprojTrXp(t1f, v1, v2p, ls, 0, 0);
      sprojTrXm(t2f, v2, v1p, ls, 0, 0);
      break;
    case 1:
      sprojTrYp(t1f, v1, v2p, ls, 0, 0);
      sprojTrYm(t2f, v2, v1p, ls, 0, 0);
      break;
    case 2:
      sprojTrZp(t1f, v1, v2p, ls, 0, 0);
      sprojTrZm(t2f, v2, v1p, ls, 0, 0);
      break;
    default:
      sprojTrTp(t1f, v1, v2p, ls, 0, 0);
      sprojTrTm(t2f, v2, v1p, ls, 0, 0);
    }

    t1 += t2;
    bfm_evo_aux::mStarDotMTransEqual( (Float*)force, (Float*)&gauge, t1f);
}



FforceWilsonType::FforceWilsonType(Matrix *momentum, Matrix *gauge_field,
                                   Float *vec1, Float *vec2, int Ls, Float dt)
:cname("FforceWilsonType")
{
    lcl[0] = GJP.XnodeSites();
    lcl[1] = GJP.YnodeSites();
    lcl[2] = GJP.ZnodeSites();
    lcl[3] = GJP.TnodeSites();
    // This has to be passed in as a parameter (say, we are using 4D
    // fermions together with 5D fermions, then GJP.SnodeSites() is
    // not reliable).
    lcl[4] = Ls;
    vol_5d = (size_t)lcl[0] * lcl[1] * lcl[2] * lcl[3] * lcl[4];
    f_size=(size_t)SPINOR_SIZE * vol_5d ;

    mom = momentum;
    gauge = gauge_field;
    v1 = vec1;
    v2 = vec2;
    Vector *v_tmp = (Vector *)v1;
    VRB.Debug(cname,cname,"v1=%g\n",v_tmp->NormSqGlbSum(f_size));
    v_tmp = (Vector *)v2;
    VRB.Debug(cname,cname,"v2=%g\n",v_tmp->NormSqGlbSum(f_size));
    coef = dt;
    if(cps::GJP.Gparity1fX()){
	if(!UniqueID()) printf("FforceWilsonType::FforceWilsonType : Gparity1fX dt *= 2\n");
	coef*=2; //This is for testing between 1f and 2f approaches
    }
    bufsize = 0;
    for(int mu = 0; mu < 4; ++mu) {
        surf_size[mu] = SPINOR_SIZE * (vol_5d / lcl[mu]);
	if(GJP.Gparity()) surf_size[mu] *= 2; //stack second flavour after the first, offset for second flavour == SPINOR_SIZE * (vol_5d / lcl[mu]);
	VRB.Debug(cname,cname,"surf_size[%d]=%d\n",mu,surf_size[mu]);
        v1so[mu] = bufsize;
        v2so[mu] = bufsize + surf_size[mu];
        surf_size[mu] *= 2;
        bufsize += surf_size[mu];
    }

    sndbuf = new Float[bufsize];
    rcvbuf = new Float[bufsize];

    if(sndbuf == NULL || rcvbuf == NULL) {
        ERR.Pointer("FforceWilsonType", "FforceWilsonType", "snd/rcv buf");
    }
}

FforceWilsonType::~FforceWilsonType()
{
    delete[] sndbuf;
    delete[] rcvbuf;
}

// copy 3d surface data from v4d to v3d in mu direction, use this
// function to fill the buffer v4d before communication.
//
// NOTE:
//
// 1. v4d is assumed to be in sxyzt order, i.e., the s index
// changes fastest.
//
// 2. we always send in negative direction.
void FforceWilsonType::collect_surface(int mu)
{
    const char *fname="collect_surface(i)";
    long l[4] = { 0, 0, 0, 0 };
    long h[4] = { lcl[0], lcl[1], lcl[2], lcl[3] };
    h[mu] = 1;

    const long block = SPINOR_SIZE * lcl[4];
    const long sites = h[0] * h[1] * h[2] * h[3];

    Float *v1s = sndbuf + v1so[mu];
    Float *v2s = sndbuf + v2so[mu];

    if(!GJP.Gparity()){
#pragma omp parallel for 
	for(long i = 0; i < sites; ++i) {
	    long x[4];
	    compute_coord(x, h, l, i);
	    long o4d = idx_4d(x, lcl);
	    long o3d = idx_4d_surf(x, lcl, mu);
	    VRB.Debug(cname,fname,"memcpy(%p+%d*%d,%p+%d*%d,sizeof(Float) *%d\n",v1s,o3d,block,v1,o4d,block,block);
	    memcpy(v1s + o3d * block, v1  + o4d * block, sizeof(Float) * block);
	    VRB.Debug(cname,fname,"memcpy(%p+%d*%d,%p+%d*%d,sizeof(Float) *%d\n",v2s,o3d,block,v2,o4d,block,block);
	    memcpy(v2s + o3d * block, v2  + o4d * block, sizeof(Float) * block);
	}
    }else{
	//G-parity flavour twist at global lattice boundary
	const long f_off = lcl[0] * lcl[1] * lcl[2] * lcl[3];
	const long f_off_s = surf_size[mu]/4/block;

	int fsites = sites*2;
#pragma omp parallel for 
	for(long fi = 0; fi < fsites; ++fi) {
	    int flav = fi / sites;
	    long i = fi % sites;

	    long x[4];
	    compute_coord(x, h, l, i);

	    long o4d = idx_4d(x, lcl) + flav * f_off;

	    int dest_flav = flav;
	    if(GJP.Bc(mu) == BND_CND_GPARITY && GJP.NodeCoor(mu) == 0) dest_flav = 1-flav; //swap destination flavour
	    	    
	    long o3d = idx_4d_surf(x, lcl, mu) + dest_flav * f_off_s;
	    
	    memcpy(v1s + o3d * block, v1  + o4d * block,
		   sizeof(Float) * block);
	    memcpy(v2s + o3d * block, v2  + o4d * block,
		   sizeof(Float) * block);
	}	
    }
}

void FforceWilsonType::comm(int mu)
{
    getPlusData(rcvbuf + v1so[mu], sndbuf + v1so[mu],
                surf_size[mu], mu);
}

ForceArg FforceWilsonType::do_internal(int mu, int nthreads)
{
    long low[4] = { 0, 0, 0, 0 };
    long high[4] = { lcl[0], lcl[1], lcl[2], lcl[3] };
    --high[mu];
    const long sites = high[0] * high[1] * high[2] * high[3];
    const long block = SPINOR_SIZE * lcl[4];

    const long f_off = lcl[0] * lcl[1] * lcl[2] * lcl[3]; //G-parity flavour offset

    int me = omp_get_thread_num();
    long mywork, myoff;
    // some threads are used in communication
    thread_work_partial(sites, me, nthreads, mywork, myoff);

    ForceArg ret;
    for(long i = 0; i < mywork; ++i) {
        long x[4];
        compute_coord(x, high, low, i + myoff);
        long off_4d = idx_4d(x, lcl);
        long gid = mu + 4 * off_4d;
        long fid = block * off_4d;

        long y[4] = {x[0], x[1], x[2], x[3]};
        ++y[mu];
        long fidp = block * idx_4d(y, lcl);

        Matrix force;
        do_site_force(&force, gauge[gid],
                      v2 + fid, v2 + fidp,
                      v1 + fid, v1 + fidp, mu, lcl[4]);

        force.TrLessAntiHermMatrix();
        force *= -coef;
        *(mom + gid) += force;
        
	if(GJP.Gparity()){
	    //Calculate and add the force from the second G-parity flavour
	    fid += block * f_off;
	    fidp += block * f_off;
	    int gid_f1 = gid + 4 * f_off; //U* links

	    Matrix force_f1;
	    do_site_force_f1(&force_f1, gauge[gid_f1],
			  v2 + fid, v2 + fidp,
			  v1 + fid, v1 + fidp, mu, lcl[4]);

	    force_f1.TrLessAntiHermMatrix();
	    force_f1 *= -coef;
	    *(mom + gid) += force_f1;

	    //for reporting purposes, add force and force_f1 together such that the L* values are for the total force from this action
	    Matrix total_force = force; total_force += force_f1;
	    updateForce(&ret, total_force);

	    //Update momentum for the second flavour using complex conjugate of force
	    for(int c=1;c<18;c+=2){ ((Float*)&force)[c]*=-1; ((Float*)&force_f1)[c]*=-1; }
	    *(mom + gid_f1) += force;
	    *(mom + gid_f1) += force_f1;
	}else{
	    updateForce(&ret, force);
	}

    }

    return ret;
}

ForceArg FforceWilsonType::do_surface(int mu, int nthreads)
{
    long low[4] = { 0, 0, 0, 0 };
    long high[4] = { lcl[0], lcl[1], lcl[2], lcl[3] };
    low[mu] = lcl[mu] - 1;
    high[mu] = lcl[mu];
    const long hl[4] = {high[0] - low[0],
                        high[1] - low[1],
                        high[2] - low[2],
                        high[3] - low[3] };
    const long sites = hl[0] * hl[1] * hl[2] * hl[3];
    const long block = SPINOR_SIZE * lcl[4];

    const long f_off = lcl[0] * lcl[1] * lcl[2] * lcl[3];
    const long f_off_s = surf_size[mu]/4/block;

    int me = omp_get_thread_num();
    long mywork, myoff;
    thread_work_partial(sites, me, nthreads, mywork, myoff);

    Float *v1_s = rcvbuf + v1so[mu];
    Float *v2_s = rcvbuf + v2so[mu];

    ForceArg ret;
    for(long i = 0; i < mywork; ++i) {
        long x[4];
        compute_coord(x, hl, low, i + myoff);

        long off_4d = idx_4d(x, lcl);
        long gid = mu + 4 * off_4d;
        long fid = block * off_4d;
        long fid_s = block * idx_4d_surf(x, lcl, mu);

        Matrix force;
        do_site_force(&force, gauge[gid],
                      v2 + fid, v2_s + fid_s,
                      v1 + fid, v1_s + fid_s, mu, lcl[4]);

        force.TrLessAntiHermMatrix();
        force *= -coef;
        *(mom + gid) += force;

	if(GJP.Gparity()){
	    //Calculate and add the force from the second G-parity flavour
	    fid += block * f_off;
	    fid_s += block * f_off_s; //surface data layout:  | v1_f0 | v1_f1 | v2_f0 | v2_f1 |

	    int gid_f1 = gid + 4 * f_off; //U* links

	    Matrix force_f1;
	    do_site_force_f1(&force_f1, gauge[gid_f1],
			  v2 + fid, v2_s + fid_s,
			  v1 + fid, v1_s + fid_s, mu, lcl[4]);

	    force_f1.TrLessAntiHermMatrix();
	    force_f1 *= -coef;
	    *(mom + gid) += force_f1;

	    //for reporting purposes, add force and force_f1 together such that the L* values are for the total force from this action
	    Matrix total_force = force; total_force += force_f1;
	    updateForce(&ret, total_force);

	    //Update momentum for the second flavour using complex conjugate of force
	    for(int c=1;c<18;c+=2){ ((Float*)&force)[c]*=-1; ((Float*)&force_f1)[c]*=-1; }
	    *(mom + gid_f1) += force;
	    *(mom + gid_f1) += force_f1;
	}else{
	    updateForce(&ret, force);
	}

    }
    return ret;
}

ForceArg FforceWilsonType::run()
{
    VRB.Func(cname,"run()");
    for(int mu = 0; mu < 4; ++mu) {
        collect_surface(mu);
    }

    // single threaded comm
    for(int mu = 0; mu < 4; ++mu) {
        comm(mu);
    }

    ForceArg ret;

#pragma omp parallel
    {
        int nthreads = omp_get_num_threads();
        ForceArg f_arg; // threaded

        for(int mu = 0; mu < 4; ++mu) {
            f_arg.combine(do_internal(mu, nthreads));
            f_arg.combine(do_surface(mu, nthreads));
        }

#pragma omp critical
        {
            ret.combine(f_arg);
        }
    }

    //GPARITY TESTING: COMPARE 1F AND 2F METHODS (NOT USED IN PRODUCTION CODE)
    if(cps::GJP.Gparity1fX()){ //this is not a threaded region
	printf("Patching up 1f G-parity force\n");

	//want  p_0' = p_0 + delta p_0 + cconj(delta p_1)
	//      p_1' = p_1 + delta p_1 + cconj(delta p_0)
	//we did p_i' = p_i + 2 * delta p_i
	//and we know p_1 = cconj(p_0)
	//so we now do  p_0' = 0.5* p_0' + 0.5* cconj(p_1') = 0.5*p_0 + delta p_0 + cconj(0.5*p_1 + delta p_1) = p0 + delta p_0 + cconj(delta p_1)
	//so we now do  p_1' = 0.5* p_1' + 0.5* cconj(p_0') = .....
	//to fix this
	int momsz = 4*18*cps::GJP.VolNodeSites();
	Float *buf = (Float *)pmalloc(momsz * sizeof(Float) );

	for(int mu=0; mu<4; mu++){
	    for(int ii=0;ii<momsz;ii++) buf[ii] = 0.0;

	    //Communicate \delta p from first half onto second half and vice versa
	    Float *data_buf = (Float*)mom;
	    Float *send_buf = data_buf;
	    Float *recv_buf = buf;

	    if(cps::GJP.Xnodes()>1){
		//pass between nodes
		for(int i=0;i<cps::GJP.Xnodes()/2;i++){
		    cps::getMinusData((Float *)recv_buf, (Float *)send_buf, momsz , 0);
		    data_buf = recv_buf;
		    recv_buf = send_buf;
		    send_buf = data_buf;
		}
	    }else{
		//shift mom[mu] field by xsites/2
		for(size_t i=0;i<cps::GJP.VolNodeSites();i++){
		    //i = (x + Lx*(y+Ly*(z+Lz*t) ) )
		    int x = i % cps::GJP.XnodeSites();
		    int pos_rem = i/cps::GJP.XnodeSites(); //(y+Ly*(z+Lz*t)

		    int x_from = (x + cps::GJP.XnodeSites()/2) % cps::GJP.XnodeSites();

		    int i_from = 18*mu + 18*4*(x_from + cps::GJP.XnodeSites()*pos_rem);
		    int i_to = 18*mu + 18*4*i;

		    for(int j=0;j<18;j++) buf[i_to+j] = ((Float*)mom)[i_from+j];
		}
		data_buf = buf;
	    }
	    for(int i=0;i<cps::GJP.VolNodeSites();i++){ //do fixup step
		int mat_off = 18*mu + 18*4*i;
      
		for(int j=0;j<18;j++){
		    if(j%2==0) ((Float*)mom)[mat_off+j] = ((Float*)mom)[mat_off+j]/2.0 + data_buf[mat_off+j]/2.0;
		    else ((Float*)mom)[mat_off+j] = ((Float*)mom)[mat_off+j]/2.0 - data_buf[mat_off+j]/2.0;
		}
	    }
	}

	pfree(buf);
    }


    ret.glb_reduce();
    VRB.FuncEnd(cname,"run()");
    return ret;
}
