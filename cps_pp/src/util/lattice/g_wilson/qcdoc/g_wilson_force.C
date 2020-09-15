#include <config.h>
#include <math.h>
#include <util/gjp.h>
#include <util/lattice.h>
#include <util/smalloc.h>
#include <util/pt.h>
#include <util/time_cps.h>
#include <cassert>
CPS_START_NAMESPACE
#define PROFILE
//------------------------------------------------------------------
// EvolveMomGforce(Matrix *mom, Float dt):
// It evolves the canonical momentum mom by dt
// using the pure gauge force.
//------------------------------------------------------------------
static const Float invs3 = -1. / 3.;
ForceArg Gwilson::EvolveMomGforce (Matrix * mom, Float dt)
{
  char *fname = "EvolveMomGforce(M*,F)";
  VRB.Func (cname, fname);
  static Matrix mt0;
  static Matrix *mp0 = &mt0;

  Float L1 = 0.0;
  Float L2 = 0.0;
  Float Linf = 0.0;

#ifdef PROFILE
  Float time = -dclock ();
  ForceFlops = 0;
  ParTrans::PTflops = 0;
#endif
  static int vol = GJP.VolNodeSites ();
  int mu, nu;
  int dirs_p[] = { 0, 2, 4, 6, 0, 2, 4 };
  int dirs_m[] = { 1, 3, 5, 7, 1, 3, 5 };

  int if_block = 0;
  if (SigmaBlockSize () > 0)
    if_block = 1;

  const int N = 4;
//  const Float SMALL = 0.001;
  Float tmp = GJP.Beta () * invs3;
  Matrix *Unit = (Matrix *) fmalloc (vol * sizeof (Matrix));
  Matrix *tmp1[N];
  Matrix *tmp2[N];
  Matrix *result[4];
  for (int i = 0; i < 4; i++) {
    result[i] = (Matrix *) fmalloc (vol * sizeof (Matrix));
  }
  for (int i = 0; i < N; i++) {
    tmp1[i] = (Matrix *) fmalloc (vol * sizeof (Matrix));
    tmp2[i] = (Matrix *) fmalloc (vol * sizeof (Matrix));
#ifdef C11
    memset ((char *) tmp2[i], 0, vol * sizeof (Matrix));
#else
    bzero ((char *) tmp2[i], vol * sizeof (Matrix));
#endif
  }
  for (int i = 0; i < vol; i++)
    Unit[i] = 1.;
  Matrix *Units[N];
  for (int i = 0; i < N; i++)
    Units[i] = Unit;
  if (!if_block && (fabs (delta_beta) > 0))
    ERR.General (cname, fname, "Not implemented for non-blocked noisy MC\n");

  LatData Plaqs (1);
  if (if_block) {
    ParTransGauge pt (*this);
    for (mu = 0; mu < 4; mu++)
      for (nu = mu + 1; nu < 4; nu++) {
	pt.run (1, tmp1, Units, dirs_m + mu);
	pt.run (1, result, tmp1, dirs_m + nu);
	pt.run (1, tmp1, result, dirs_p + mu);
	pt.run (1, result, tmp1, dirs_p + nu);
#pragma omp parallel for
	for (int i = 0; i < vol; i++) {
	  Float *tmp_f = Plaqs.Field (i);
	  Float re_tr = (result[0] + i)->ReTr ();
	  assert (re_tr >= 0.);
	  if (mu == 0 && nu == 1)
	    *tmp_f = re_tr;
	  else
	    *tmp_f += re_tr;
	  if (i == 0)
	    VRB.Result (cname, fname, "ReTr(Plaq)[%d][%d][0]=%0.12e\n", mu, nu,
			re_tr);
	}
      }
  }


  if (if_block) {
    int x[4];
    for (x[0] = 0; x[0] < node_sites[0]; x[0] += sigma_blocks[0])
      for (x[1] = 0; x[1] < node_sites[1]; x[1] += sigma_blocks[1])
	for (x[2] = 0; x[2] < node_sites[2]; x[2] += sigma_blocks[2])
	  for (x[3] = 0; x[3] < node_sites[3]; x[3] += sigma_blocks[3]) {
	    int offset[4], x_tmp[4];
	    Float re_tr_plaq = 0.;
	    for (offset[0] = 0; offset[0] < sigma_blocks[0]; offset[0] += 1)
	      for (offset[1] = 0; offset[1] < sigma_blocks[1]; offset[1] += 1)
		for (offset[2] = 0; offset[2] < sigma_blocks[2]; offset[2] += 1)
		  for (offset[3] = 0; offset[3] < sigma_blocks[3];
		       offset[3] += 1) {
		    for (int i = 0; i < 4; i++)
		      x_tmp[i] = x[i] + offset[i];
		    re_tr_plaq += *(Plaqs.Field (GsiteOffset (x_tmp) / 4));
		  }

	    for (offset[0] = 0; offset[0] < sigma_blocks[0]; offset[0] += 1)
	      for (offset[1] = 0; offset[1] < sigma_blocks[1]; offset[1] += 1)
		for (offset[2] = 0; offset[2] < sigma_blocks[2]; offset[2] += 1)
		  for (offset[3] = 0; offset[3] < sigma_blocks[3];
		       offset[3] += 1) {
		    for (int i = 0; i < 4; i++)
		      x_tmp[i] = x[i] + offset[i];
		    *(Plaqs.Field (GsiteOffset (x_tmp) / 4)) = re_tr_plaq;
		  }
	  }

  }



  {
    ParTransGauge pt (*this);
    LatMatrix Scale[N];
    Matrix *tmp3[N];
    for (int i = 0; i < N; i++)
      tmp3[i] = Scale[i].Mat ();
    Float beta = GJP.Beta ();


    for (nu = 1; nu < 4; nu++) {
      for (int i = 0; i < N; i++) {
#pragma omp parallel for
	for (int j = 0; j < vol; j++) {
	  int sigma = *(SigmaField () + SigmaOffset (j, i, (i + nu) % 4));
	  Float re_tr_plaq = *(Plaqs.Field (j));
	  if (sigma == 0)
	    *(tmp3[i] + j) = 1.0 + delta_beta / beta * DeltaSDer (re_tr_plaq);
	  else {
	    Float exponent = DeltaS (re_tr_plaq);
//                      if (exponent< SMALL) exponent = SMALL;
	    *(tmp3[i] + j) =
	      1.0 -
	      delta_beta / (beta * (exp (exponent) - 1.0)) *
	      DeltaSDer (re_tr_plaq);
	  }

	}
      }
      pt.run (N, tmp1, tmp3, dirs_m + nu);
      pt.run (N, result, tmp1, dirs_m);
      pt.run (N, tmp1, result, dirs_p + nu);
      for (int i = 0; i < N; i++) {
	tmp2[i]->FTimesV1PlusV2 (tmp, tmp1[i], tmp2[i], vol);
      }
      pt.run (N, tmp1, Units, dirs_p + nu);
      for (int i = 0; i < N; i++) {
#pragma omp parallel for
	for (int j = 0; j < vol; j++) {
	  int sigma = *(SigmaField () + SigmaOffset (j, i, (i + nu) % 4));
	  Float re_tr_plaq = *(Plaqs.Field (j));
	  if (sigma == 0)
	    *(tmp1[i] + j) *= 1.0 + delta_beta / beta * DeltaSDer (re_tr_plaq);
	  else {
	    Float exponent = DeltaS (re_tr_plaq);
//                      if (exponent< SMALL) exponent = SMALL;
	    *(tmp1[i] + j) *=
	      1.0 -
	      delta_beta / (beta * (exp (exponent) - 1.0)) *
	      DeltaSDer (re_tr_plaq);
	  }

	}
      }
      pt.run (N, result, tmp1, dirs_m);
      pt.run (N, tmp1, result, dirs_m + nu);
      for (int i = 0; i < N; i++) {
	tmp2[i]->FTimesV1PlusV2 (tmp, tmp1[i], tmp2[i], vol);
      }
      ForceFlops += vol * 12 * N;
    }
    pt.run (N, result, tmp2, dirs_p);
  }
#if 1
#pragma omp parallel for default(shared) private(mu) reduction(+:L1,L2)
  for (int index = 0; index < 4 * vol; index++) {
    Matrix mp1;
    int i = index % vol;
    mu = index / vol;
    Matrix *mtmp = (result[mu] + i);
    mp1.Dagger ((IFloat *) mtmp);
    mtmp->TrLessAntiHermMatrix (mp1);
    IFloat *ihp = (IFloat *) (mom + i * 4 + mu);	//The gauge momentum
//    IFloat *dotp = (IFloat *)mp0;
    IFloat *dotp2 = (IFloat *) (result[mu] + (i));
    assert (mtmp->norm () >= 0.);
    if (i < 4)
      VRB.Result (cname, fname, "Gforce[%d][%d]=%0.12e\n", mu, i,
		  mtmp->norm ());
    fTimesV1PlusV2Single (ihp, dt, dotp2, ihp, 18);	//Update the gauge momentum
    Float norm = ((Matrix *) dotp2)->norm ();
    Float tmp = sqrt (norm);
    L1 += tmp;
    L2 += norm;
  }


#else
  Matrix mp1;
  for (mu = 0; mu < 4; mu++) {
    Matrix *mtmp = result[mu];
    for (int i = 0; i < vol; i++) {
      mtmp->TrLessAntiHermMatrix ();
      mtmp++;
    }
  }

  int x[4];

  for (x[0] = 0; x[0] < GJP.XnodeSites (); ++x[0]) {
    for (x[1] = 0; x[1] < GJP.YnodeSites (); ++x[1]) {
      for (x[2] = 0; x[2] < GJP.ZnodeSites (); ++x[2]) {
	for (x[3] = 0; x[3] < GJP.TnodeSites (); ++x[3]) {

	  int uoff = GsiteOffset (x);

	  for (int mu = 0; mu < 4; ++mu) {

	    IFloat *ihp = (IFloat *) (mom + uoff + mu);
	    IFloat *dotp2 = (IFloat *) (result[mu] + (uoff / 4));
	    fTimesV1PlusV2 (ihp, dt, dotp2, ihp, 18);
	    Float norm = ((Matrix *) dotp2)->norm ();
	    Float tmp = sqrt (norm);
	    L1 += tmp;
	    L2 += norm;
	    Linf = (tmp > Linf ? tmp : Linf);
	  }
	}
      }
    }
  }
#endif
  ForceFlops += vol * 60;
  ForceFlops += vol * 144;

#ifdef PROFILE
  time += dclock ();
  print_flops (cname, fname, ForceFlops + ParTrans::PTflops, time);
#endif

  ffree (Unit);
  for (int i = 0; i < N; i++) {
    ffree (tmp1[i]);
    ffree (tmp2[i]);
  }
  for (int i = 0; i < 4; i++)
    ffree (result[i]);

  glb_sum (&L1);
  glb_sum (&L2);
  glb_max (&Linf);

  L1 /= 4.0 * GJP.VolSites ();
  L2 /= 4.0 * GJP.VolSites ();

  return ForceArg (dt * L1, dt * sqrt (L2), dt * Linf);

}

CPS_END_NAMESPACE
