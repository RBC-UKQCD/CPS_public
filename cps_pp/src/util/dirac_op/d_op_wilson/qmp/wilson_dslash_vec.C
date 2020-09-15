#include <config.h>
#include <stdio.h>
#include <math.h>
#include <qmp.h>
#include <util/dirac_op.h>
#include <util/omp_wrapper.h>
#include <util/timer.h>

#warning "Using vectorised wilson dslash"

CPS_START_NAMESPACE
/*! \file
  \brief  Routine used internally in the DiracOpWilson class.

*/
/***************************************************************************/
/*                                                                         */
/* wilson_dslash: It calculates chi = Dslash * psi, or                     */
/*                                chi = DslashDag * psi, where Dslassh is  */
/* the Wilson fermion Dslash matrix. Dslash is a function of the gauge     */
/* fields u.                                                               */
/* cb = 0/1 denotes input on even/odd checkerboards i.e.                   */
/* cb = 1 --> Dslash_EO, acts on an odd column vector and produces an      */
/* even column vector                                                      */
/* cb = 0 --> Dslash_OE, acts on an even column vector and produces an     */
/* odd column vector,                                                      */
/*                                                                         */
/* dag = 0/1  results into calculating Dslash or DslashDagger.             */
/* lx,ly,lz,lt is the lattice size.                                        */
/*                                                                         */
/* This routine is to be used with scalar machines.                        */
/*                                                                         */
/* WARNING:                                                                 */
/*                                                                          */
/* This set of routines will work only if the node sublattices have         */
/* even number of sites in each direction.                                  */
/*                                                                          */
/***************************************************************************/
  CPS_END_NAMESPACE
#include <util/data_types.h>
#include <util/vector.h>
#include <util/gjp.h>
#include <util/wilson.h>
#include <util/error.h>
#include <util/verbose.h>
#include <util/smalloc.h>
#include <util/time_cps.h>
#include <comms/sysfunc_cps.h>
  CPS_START_NAMESPACE
#include "wilson_op.h"
static int Printf (char *format, ...)
{
}

//#define Printf if ( QMP_is_primary_node() ) printf
//#define Printf printf

//! Access to the elements of the \e SU(3) matrix
/*!
  Gets the element of the \e SU(3) matrix \e u with row \e row,
  column \e col and complex component \e d
*/
#define UELEM(u,r,row,col,d) *(u+(r+2*(row+3*(col+3*d))))
//! Access to the elements of a spinor vector.
/*!
  Gets the element of the spinor \e f with spin \e s,
  colour \e c and complex component \e r
*/
#define FERM(f,r,c,s) *(f +(r+2*(c+3*s)))

static unsigned long called = 0;
//static int initted=0;
static int init_len = 0;
static double setup = 0;
static double local = 0;
static double nonlocal = 0;
static double qmp = 0;
//keep track of whether initializations were performed in the G-parity framework. If they were and we switch out of G-parity
//mode, we want to reset the initialization status
static int gparity_init_status = 0;

static inline void MOVE_VEC (IFloat * buf, IFloat * psi, int vec_len,
                             unsigned long vec_offset)
{
  for (int i = 0; i < vec_len; i++) {
    moveMem (buf, psi, SPINOR_SIZE * sizeof (IFloat) / sizeof (char));
    buf += SPINOR_SIZE;
    psi += vec_offset;
  }
}

static inline void MOVE_VEC2 (IFloat * buf, IFloat * psi, int vec_len,
                              unsigned long vec_offset)
{
  for (int i = 0; i < vec_len; i++) {
    moveMem (buf, psi, SPINOR_SIZE * sizeof (IFloat) / sizeof (char));
    buf += vec_offset;
    psi += vec_offset;
  }
}



void wilson_dslash_vec (IFloat * chi_p_f,
                        IFloat * u_p_f,
                        IFloat * psi_p_f,
                        int cb,
                        int dag,
                        Wilson * wilson_p,
                        int vec_len, unsigned long vec_offset)
{
  //if(!UniqueID()){ printf("Running vectorised wilson dslash\n"); fflush(stdout); }


  char *cname = "";
  char *fname = "wilson_dslash_vec";
  static Timer timer_setup (fname,"setup()");
  static Timer timer_local (fname,"local()");
  static Timer timer_nl (fname,"nonlocal()");
  static Timer timer_qmp (fname,"qmp()");

  int lx, ly, lz, lt;
//  int mu;
  //    int r, c, s;
  int vol;

  int temp_size = SPINOR_SIZE;
  if (GJP.Gparity ())
    temp_size *= 2;

  Float fbuf[temp_size];

//  Float dtime = -dclock (true);
  timer_setup.start();
  for (int i = 0; i < temp_size; i++)
    fbuf[i] = 0.;

  /*--------------------------------------------------------------------------*/
  /* Initializations                                                          */
  /*--------------------------------------------------------------------------*/
  int sdag;                     /* = +/-1 if dag = 0/1 */
  int cbn;
  Float *chi_p = (Float *) chi_p_f;
  Float *u_p = (Float *) u_p_f;
  Float *psi_p = (Float *) psi_p_f;


  lx = wilson_p->ptr[0];
  ly = wilson_p->ptr[1];
  lz = wilson_p->ptr[2];
  lt = wilson_p->ptr[3];
  vol = wilson_p->vol[0];

  int lattsz[4] = { lx, ly, lz, lt };

  static int local_comm[4];

  if (dag == 0)
    sdag = 1;
  else if (dag == 1)
    sdag = -1;
  else {
    ERR.General (" ", fname, "dag must be 0 or 1");
  }

  //cb is the fermion parity, cbn is the gauge field parity
  if (cb == 0)
    cbn = 1;
  else if (cb == 1)
    cbn = 0;
  else {
    ERR.General (" ", fname, "cb must be 0 or 1");
  }

  Float *ind_buf[8];
  unsigned long ind_nl[8];

  static unsigned long num_nl[8];       //number non-local sites
  static unsigned long *u_ind[8];
  static unsigned long *f_ind[8];
  static unsigned long *t_ind[8];
  static Float *Send_buf[8];
  static Float *Recv_buf[8];
  QMP_msgmem_t Send_mem[8];
  QMP_msgmem_t Recv_mem[8];
//  QMP_msghandle_t Send[8];
//  QMP_msghandle_t Recv[8];
  QMP_msghandle_t multiple_p[16];
  QMP_msghandle_t multiple;
  int n_dir=0;


  //reset setup variables if G-parity status changes
  if (gparity_init_status == 1 && !GJP.Gparity () || gparity_init_status == 0
      && GJP.Gparity ())
    init_len = -1;

  if (vec_len != init_len) {
    //These are only initialized once unless the vec_len changes
    //For G-parity tests where we change the volume, I have added a reset for when GJP.Gparity changes above
    //coupled with the initialization status below
    if (GJP.Gparity ())
      gparity_init_status = 1;
    else
      gparity_init_status = 0;

    VRB.Result (cname, fname, "init_len(%d)!=vec_len(%d)\n", init_len, vec_len);
    if (init_len != 0)
      for (int i = 0; i < 8; i++) {
        if (u_ind[i]) sfree (u_ind[i]);
        if (f_ind[i]) sfree (f_ind[i]);
        if (Send_buf[i]) sfree (Send_buf[i]);
        if (Recv_buf[i]) sfree (Recv_buf[i]);
      }

    num_nl[0] = num_nl[4] = vol / lx;
    num_nl[1] = num_nl[5] = vol / ly;
    num_nl[2] = num_nl[6] = vol / lz;
    num_nl[3] = num_nl[7] = vol / lt;

    if (!UniqueID ()) {
      printf
        ("Lattice is %d x %d x %d x %d, volume is %d, cb volume %d. Input 'vol' is %d.\n",
         GJP.XnodeSites (), GJP.YnodeSites (), GJP.ZnodeSites (),
         GJP.TnodeSites (), GJP.VolNodeSites (), GJP.VolNodeSites () / 2, vol);
//      fflush(stdout);
    }

    for (int i = 0; i < 4; i++) {
      if (GJP.Nodes (i) > 1) {
        local_comm[i] = 0;
        if (!UniqueID ()) {
          printf ("Non-local in direction %d, sites on wall %d\n", i,
                  num_nl[i]);
          fflush (stdout);
        }
      }

      else
        local_comm[i] = 1;
    }
    for (int i = 0; i < 8; i++) {
      int mu = i % 4;
      int sign = 1 - 2 * (i / 4);       //1,1,1,1,-1,-1,-1,-1
      if (local_comm[mu]) {
        num_nl[i] = 0;
        u_ind[i] = NULL;
        f_ind[i] = NULL;
        Send_buf[i] = NULL;
        Recv_buf[i] = NULL;
      } else {
        u_ind[i] = (unsigned long *) smalloc (cname, fname, "u_ind[i]", num_nl[i] * sizeof (unsigned long));
        f_ind[i] = (unsigned long *) smalloc (cname, fname, "f_ind[i]", num_nl[i] * sizeof (unsigned long));
        t_ind[i] = (unsigned long *) smalloc (cname, fname, "f_ind[i]", num_nl[i] * sizeof (unsigned long));
        size_t buf_size = num_nl[i] * SPINOR_SIZE * sizeof (Float) * vec_len;
        if (GJP.Gparity ()) buf_size *= 2;        //use second half of buffer for f1
        Send_buf[i] = (Float *) smalloc (cname, fname, "Send_buf[i]", buf_size);
        Recv_buf[i] = (Float *) smalloc (cname, fname, "Recv_buf[i]", buf_size);
      }
    }
    VRB.Flow (cname,fname,"initted\n");
    init_len = vec_len;
  }
  n_dir=0;
  for (int i = 0; i < 8; i++) {
      int mu = i % 4;
      int sign = 1 - 2 * (i / 4);       //1,1,1,1,-1,-1,-1,-1
      if (!local_comm[i%4]) {
        size_t buf_size = num_nl[i] * SPINOR_SIZE * sizeof (Float) * vec_len;
        Send_mem[n_dir] = QMP_declare_msgmem (Send_buf[i], buf_size);
        multiple_p[n_dir*2] = QMP_declare_send_relative (Send_mem[n_dir], mu, -sign, 0);
        Recv_mem[n_dir] = QMP_declare_msgmem (Recv_buf[i], buf_size);
        multiple_p[n_dir*2+1] = QMP_declare_receive_relative (Recv_mem[n_dir], mu, sign, 0);
	n_dir ++;
      }
  }
  if(n_dir>0) multiple = QMP_declare_multiple(multiple_p,2*n_dir);
  for (int i = 0; i < 8; i++) {
    ind_nl[i] = 0;
    ind_buf[i] = Send_buf[i];
    //if(!local_comm[i%4]){ printf("Node %d, non-local mu = %d send buf %i at %p\n",UniqueID(),i%4,i,ind_buf[i]); fflush(stdout); }
  }

  //printf("Node %d, local_comm=%d %d %d %d\n",UniqueID(),local_comm[0],local_comm[1],local_comm[2],local_comm[3]); fflush(stdout);

  if (called % 100000 == 0) {
    VRB.Result (cname, fname, "local_comm=%d %d %d %d\n", local_comm[0],
                local_comm[1], local_comm[2], local_comm[3]);
    if (called > 0)
      called = 0;
  }

  int u_cboff = vol;            //vol is the checkerboard volume, i.e. half the 4d volume
  if (GJP.Gparity ())
    u_cboff *= 2;               //[(f0 odd)(f1 odd)(f0 even)(f1 even)]  each bracket is one cbvolume

  //
  //  non-local send
  //
  //
  GJP.SetNthreads ();

#pragma omp parallel for default(shared)
  for (int dir = 0; dir < 8; dir++) {
    int x, y, z, t;
    int r, c, s;
    int xyzt;
    int parity;
    int mu;
    Float *chi;
    Float *u;
    Float *psi;

    //CK:for G-parity
    Float *chi_f1;
    Float *u_f1;
    Float *psi_f1;

    Float tmp[temp_size];
    Float tmp1[temp_size];
    Float tmp2[temp_size];
    Float tmp3[temp_size];
    Float tmp4[temp_size];
    Float tmp5[temp_size];
    Float tmp6[temp_size];
    Float tmp7[temp_size];
    Float tmp8[temp_size];

    Float *temps[4][2] =
      { {tmp1, tmp5}, {tmp2, tmp6}, {tmp3, tmp7}, {tmp4, tmp8} };

    for (x = 0; x < lx; x++) {
      for (y = 0; y < ly; y++) {
        for (z = 0; z < lz; z++) {
          for (t = 0; t < lt; t++) {
            parity = x + y + z + t;
            parity = parity % 2;
            if (parity == cbn) {
              //for Dslash_xy\psi_y, y is opposite parity to that of sites at which \psi is taken from 
              //printf("Node %d, non-local send: %d %d %d %d\n",UniqueID(),x,y,z,t); fflush(stdout); 

              /* x,y,z,t addressing of cbn checkerboard */
              xyzt = (x / 2) + (lx / 2) * (y + ly * (z + lz * t));
              int pos[4] = { x, y, z, t };

              chi = chi_p + SPINOR_SIZE * xyzt;
              if (GJP.Gparity ())
                chi_f1 = chi_p + SPINOR_SIZE * (xyzt + vol);
              int psi_bufoff = num_nl[dir] * SPINOR_SIZE * vec_len;     //offset for G-parity f1 pointers in buffer

              if (dir < 4) {    //backwards send. Fermion site is on the left-most boundary. Send \psi{x+\mu} to x
                mu = dir;
                //printf("Node %d, dir %d, mu %d, 1, local_comm[mu] = %d\n",UniqueID(),dir,mu,local_comm[mu]); fflush(stdout); 
                int posp[4] = { x, y, z, t };
                posp[mu] = (posp[mu] + 1) % lattsz[mu];
                int posp_xyzt =
                  (posp[0] / 2) + (lx / 2) * (posp[1] +
                                              ly * (posp[2] + lz * posp[3]));
                //printf("Node %d, dir %d, mu %d, 2, local_comm[mu] = %d\n",UniqueID(),dir,mu,local_comm[mu]); fflush(stdout); 

                if ((pos[mu] == lattsz[mu] - 1) && !local_comm[mu]) {   //posp[mu]==0
                  //place fermions in buffer to prepare for backwards send
                  if (!GJP.Gparity ()) {
                    psi = psi_p + SPINOR_SIZE * posp_xyzt;

                    // printf("Node %d, non-local dir %d, moving data from %p to %p\n",UniqueID(),dir,psi,ind_buf[dir]); fflush(stdout); 
                    // printf("Node %d, local_comm[%d] = %d\n", UniqueID(),mu,local_comm[mu]);
                    // printf("Node %d, Testing locations\n",UniqueID());
                    // for(int ii=0;ii<SPINOR_SIZE;ii++){
                    //   printf("Node %d, psi[%d] = %f\n",UniqueID(),ii,psi[ii]);
                    //   fflush(stdout);
                    //   printf("Node %d, ind_buf[%d] = %f\n",UniqueID(),ii,ind_buf[dir][ii]);
                    //   fflush(stdout);
                    // }
                    // printf("Node %d, Done testing locations\n",UniqueID()); fflush(stdout);
                    MOVE_VEC (ind_buf[dir], psi, vec_len, vec_offset);

                    //printf("Node %d, dir %d, 4\n",dir,UniqueID()); fflush(stdout); 

                  } else {
                    int psi_f0_bufoff = 0;
                    int psi_f1_bufoff = psi_bufoff;
                    if (GJP.Bc (mu) == BND_CND_GPARITY
                        && GJP.NodeCoor (mu) == 0) {
                      //implement G-parity twist at boundary. Note this data is sent in the *minus* direction
                      psi_f0_bufoff = psi_bufoff;
                      psi_f1_bufoff = 0;
                    }

                    psi = psi_p + SPINOR_SIZE * posp_xyzt;
                    MOVE_VEC (ind_buf[dir] + psi_f0_bufoff, psi, vec_len,
                              vec_offset);

                    psi_f1 = psi_p + SPINOR_SIZE * (posp_xyzt + vol);
                    MOVE_VEC (ind_buf[dir] + psi_f1_bufoff, psi_f1, vec_len,
                              vec_offset);
                  }

                  *(u_ind[dir] + ind_nl[dir]) = xyzt + u_cboff * cbn;
                  *(f_ind[dir] + ind_nl[dir]) = posp_xyzt;
                  *(t_ind[dir] + ind_nl[dir]) = xyzt;
                  ind_buf[dir] += SPINOR_SIZE * vec_len;
                  ind_nl[dir]++;
                }

              } else {          //Forwards send, fermion is on right-most boundary. Send P_- U^\dagger_{x-\mu} \psi_{x-\mu} to x
                /* 1+gamma_mu */
                /*-----------*/
                mu = dir - 4;

                int posm[4] = { x, y, z, t };
                posm[mu] =
                  posm[mu] - 1 +
                  ((lattsz[mu] - posm[mu]) / lattsz[mu]) * lattsz[mu];
                int posm_xyzt =
                  (posm[0] / 2) + (lx / 2) * (posm[1] +
                                              ly * (posm[2] + lz * posm[3]));

                u = u_p + GAUGE_SIZE * (posm_xyzt + u_cboff * cb);
                psi = psi_p + SPINOR_SIZE * posm_xyzt;

                if (GJP.Gparity ()) {
                  u_f1 = u_p + GAUGE_SIZE * (posm_xyzt + u_cboff * cb + vol);
                  psi_f1 = psi_p + SPINOR_SIZE * (posm_xyzt + vol);
                }

                if ((pos[mu] == 0) && !local_comm[mu]) {        //posm on right-most boundary
                  Printf
                    ("getMinusData((IFloat *)fbuf, (IFloat *)tmp%d, SPINOR_SIZE, %d);\n",
                     dir + 1, mu);
                  *(u_ind[dir] + ind_nl[dir]) = posm_xyzt + u_cboff * cb;
                  *(f_ind[dir] + ind_nl[dir]) = posm_xyzt;
                  *(t_ind[dir] + ind_nl[dir]) = xyzt;
                  for (int vec_ind = 0; vec_ind < vec_len; vec_ind++) {
                    PLUSMU (mu, u, tmp, temps[mu][1], sdag, psi);

                    if (!GJP.Gparity ()) {
                      moveMem ((IFloat *) ind_buf[dir], (IFloat *) temps[mu][1],
                               SPINOR_SIZE * sizeof (Float) / sizeof (char));
                      psi = psi + vec_offset;
                    } else {
                      //both psi and u are at the same site. G-parity twist should swap which component of chi gets each contribution
                      PLUSMU (mu, u_f1, tmp, temps[mu][1] + SPINOR_SIZE, sdag,
                              psi_f1);

                      int psi_f0_bufoff = 0;
                      int psi_f1_bufoff = psi_bufoff;
                      if (GJP.Bc (mu) == BND_CND_GPARITY
                          && GJP.NodeCoor (mu) == GJP.Nodes (mu) - 1) {
                        //implement G-parity twist at boundary. Note this data is sent in the *plus* direction
                        psi_f0_bufoff = psi_bufoff;
                        psi_f1_bufoff = 0;
                      }

                      moveMem ((IFloat *) (ind_buf[dir] + psi_f0_bufoff),
                               (IFloat *) temps[mu][1],
                               SPINOR_SIZE * sizeof (Float) / sizeof (char));
                      psi = psi + vec_offset;

                      moveMem ((IFloat *) (ind_buf[dir] + psi_f1_bufoff),
                               (IFloat *) (temps[mu][1] + SPINOR_SIZE),
                               SPINOR_SIZE * sizeof (Float) / sizeof (char));
                      psi_f1 = psi_f1 + vec_offset;
                    }

                    ind_buf[dir] += SPINOR_SIZE;

                  }
                  ind_nl[dir]++;
                }
              }
            }
          }
        }
      }
    }
  }

  for (int i = 0; i < 8; i++) {
    if (ind_nl[i] != num_nl[i])
      VRB.Result (cname, fname, "ind_nl[%d](%d)!=num_nl[%d](%d)\n", i,
                  ind_nl[i], i, num_nl[i]);
#if 0
    if (!local_comm[i % 4]) {
      Printf ("QMP_start(Recv[%d])(%p)\n", i, Recv[i]);
      QMP_start (Recv[i]);
      Printf ("QMP_start(Send[%d])(%p)\n", i, Send[i]);
      QMP_start (Send[i]);
    }
#endif
  }
  if(n_dir>0) QMP_start (multiple);

  timer_setup.stop();
//  dtime += dclock (true);
//  setup += dtime;
  timer_local.start();

//  dtime = -dclock ();
  GJP.SetNthreads ();

  /*--------------------------------------------------------------------------*/
  /* Loop over sites                                                          */
  /*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/* Loop over sites                                                          */
/*--------------------------------------------------------------------------*/
  for (int i = 0; i < SPINOR_SIZE; i++)
    fbuf[i] = 0.;
  int index = 0;
#pragma omp parallel for default(shared)
  for (index = 0; index < vol * 2; index++) {
    //  Printf("wilson_dslash: %d %d %d %d\n",x,y,z,t);
    //  if ((called%10000==0) &&(!UniqueID())){
    //          printf("wilson_dslash: index=%d thread %d of %d\n",index,omp_get_thread_num(),omp_get_num_threads());
    //  }
    int r, c, s;
    int x, y, z, t;
    int xyzt;
    int parity;
    int mu;
    Float *chi;
    Float *u;
    Float *psi;

    //CK: for G-parity
    Float *chi_f1;
    Float *u_f1;
    Float *psi_f1;

    Float tmp[temp_size];
    Float tmp1[temp_size];
    Float tmp2[temp_size];
    Float tmp3[temp_size];
    Float tmp4[temp_size];
    Float tmp5[temp_size];
    Float tmp6[temp_size];
    Float tmp7[temp_size];
    Float tmp8[temp_size];

    Float *temps[4][2] =
      { {tmp1, tmp5}, {tmp2, tmp6}, {tmp3, tmp7}, {tmp4, tmp8} };

    int temp = index;
    x = temp % lx;
    temp = temp / lx;
    y = temp % ly;
    temp = temp / ly;
    z = temp % lz;
    temp = temp / lz;
    t = temp % lt;
    temp = temp / lt;

    if (0)
      if ((called % 1000000 == 0) && (!UniqueID ())) {
        printf ("wilson_dslash: %d %d %d %d %d: thread %d of %d tmp=%p \n",
                index, x, y, z, t, omp_get_thread_num (),
                omp_get_num_threads (), tmp);
      }


    parity = x + y + z + t;
    parity = parity % 2;
    if (parity == cbn) {

      /* x,y,z,t addressing of cbn checkerboard */
      xyzt = (x / 2) + (lx / 2) * (y + ly * (z + lz * t));
      //        VRB.Result(fname,"local", "%d %d %d %d\n",x,y,z,t);

      chi = chi_p + SPINOR_SIZE * xyzt;
      if (GJP.Gparity ())
        chi_f1 = chi_p + SPINOR_SIZE * (xyzt + vol);

      for (int vec_ind = 0; vec_ind < vec_len; vec_ind++) {
        for (s = 0; s < 4; s++)
          for (c = 0; c < 3; c++)
            for (r = 0; r < 2; r++)
              FERM (chi, r, c, s) = 0.;
        chi += vec_offset;
      }

      if (GJP.Gparity ()) {
        for (int vec_ind = 0; vec_ind < vec_len; vec_ind++) {
          for (s = 0; s < 4; s++)
            for (c = 0; c < 3; c++)
              for (r = 0; r < 2; r++)
                FERM (chi_f1, r, c, s) = 0.;
          chi_f1 += vec_offset;
        }
      }

      int pos[4] = { x, y, z, t };

      for (mu = 0; mu < 4; mu++) {
        //1-gamma_mu
        int posp[4] = { x, y, z, t };
        posp[mu] = (posp[mu] + 1) % lattsz[mu];
        int posp_xyzt =
          (posp[0] / 2) + (lx / 2) * (posp[1] + ly * (posp[2] + lz * posp[3]));

        chi = chi_p + SPINOR_SIZE * xyzt;
        u = u_p + GAUGE_SIZE * (xyzt + u_cboff * cbn);
        psi = psi_p + SPINOR_SIZE * posp_xyzt;

        if (GJP.Gparity ()) {
          chi_f1 = chi_p + SPINOR_SIZE * (xyzt + vol);
          u_f1 = u_p + GAUGE_SIZE * (xyzt + u_cboff * cbn + vol);
          psi_f1 = psi_p + SPINOR_SIZE * (posp_xyzt + vol);
        }

        if ((pos[mu] == lattsz[mu] - 1) && !local_comm[mu]) {
          //not sure if this actually *does* anything!
          psi = fbuf;
          if (GJP.Gparity ())
            psi_f1 = fbuf + SPINOR_SIZE;
        } else {
          //internal site or local comms for edge sites

          for (int vec_ind = 0; vec_ind < vec_len; vec_ind++) {
            if (!GJP.Gparity ()) {
              MINUSMU (mu, u, tmp, temps[mu][0], sdag, psi);

              for (s = 0; s < 4; s++)
                for (c = 0; c < 3; c++)
                  for (r = 0; r < 2; r++)
                    if (mu == 0)
                      FERM (chi, r, c, s) = FERM (temps[mu][0], r, c, s);
                    else
                      FERM (chi, r, c, s) += FERM (temps[mu][0], r, c, s);

              psi += vec_offset;
              chi += vec_offset;

            } else {            //G-parity
              Float *psi_f0_use = psi;
              Float *psi_f1_use = psi_f1;

              if (GJP.Bc (mu) == BND_CND_GPARITY && pos[mu] == lattsz[mu] - 1) {
                //psi crosses G-parity boundary, do flavour twist
                psi_f0_use = psi_f1;
                psi_f1_use = psi;
              }
              MINUSMU (mu, u, tmp, temps[mu][0], sdag, psi_f0_use);
              MINUSMU (mu, u_f1, tmp, temps[mu][0] + SPINOR_SIZE, sdag,
                       psi_f1_use);

              for (s = 0; s < 4; s++) {
                for (c = 0; c < 3; c++) {
                  for (r = 0; r < 2; r++) {
                    if (mu == 0) {
                      FERM (chi, r, c, s) = FERM (temps[mu][0], r, c, s);
                      FERM (chi_f1, r, c, s) =
                        FERM (temps[mu][0] + SPINOR_SIZE, r, c, s);
                    } else {
                      FERM (chi, r, c, s) += FERM (temps[mu][0], r, c, s);
                      FERM (chi_f1, r, c, s) +=
                        FERM (temps[mu][0] + SPINOR_SIZE, r, c, s);
                    }
                  }
                }
              }

              psi += vec_offset;
              chi += vec_offset;

              psi_f1 += vec_offset;
              chi_f1 += vec_offset;
            }

          }                     //vec_ind loop
        }                       //local comms loop
      }                         //mu loop

      for (mu = 0; mu < 4; mu++) {
        /* 1+gamma_mu */
        /*-----------*/
        int posm[4] = { x, y, z, t };
        posm[mu] =
          posm[mu] - 1 + ((lattsz[mu] - posm[mu]) / lattsz[mu]) * lattsz[mu];
        int posm_xyzt =
          (posm[0] / 2) + (lx / 2) * (posm[1] + ly * (posm[2] + lz * posm[3]));

        chi = chi_p + SPINOR_SIZE * xyzt;
        u = u_p + GAUGE_SIZE * (posm_xyzt + u_cboff * cb);      //note opposite parity
        psi = psi_p + SPINOR_SIZE * posm_xyzt;

        if (GJP.Gparity ()) {
          chi_f1 = chi_p + SPINOR_SIZE * (xyzt + vol);
          u_f1 = u_p + GAUGE_SIZE * (posm_xyzt + u_cboff * cb + vol);   //note opposite parity
          psi_f1 = psi_p + SPINOR_SIZE * (posm_xyzt + vol);
        }

        if (pos[mu] == 0 && (!local_comm[mu])) {
          moveMem ((IFloat *) temps[mu][1], (IFloat *) fbuf,
                   SPINOR_SIZE * sizeof (Float) / sizeof (char));

          if (GJP.Gparity ()) {
            moveMem ((IFloat *) (temps[mu][1] + SPINOR_SIZE),
                     (IFloat *) (fbuf + SPINOR_SIZE),
                     SPINOR_SIZE * sizeof (Float) / sizeof (char));
          }
        } else {
          for (int vec_ind = 0; vec_ind < vec_len; vec_ind++) {
            if (!GJP.Gparity ()) {
              PLUSMU (mu, u, tmp, temps[mu][1], sdag, psi);
              for (s = 0; s < 4; s++)
                for (c = 0; c < 3; c++)
                  for (r = 0; r < 2; r++)
                    FERM (chi, r, c, s) += FERM (temps[mu][1], r, c, s);

              psi += vec_offset;
              chi += vec_offset;
            } else {
              //here the gauge fields and psi are drawn from the same site but placed into chi at the next site over
              //hence for the G-parity twist we need to swap over the parts contributing to psi_f0 and psi_f1

              PLUSMU (mu, u, tmp, temps[mu][1], sdag, psi);
              PLUSMU (mu, u_f1, tmp, temps[mu][1] + SPINOR_SIZE, sdag, psi_f1);

              Float *chi_f0_use = chi;
              Float *chi_f1_use = chi_f1;

              if (GJP.Bc (mu) == BND_CND_GPARITY && pos[mu] == 0) {
                //contribution to chi crosses G-parity boundary, do flavour twist
                chi_f0_use = chi_f1;
                chi_f1_use = chi;
              }

              for (s = 0; s < 4; s++) {
                for (c = 0; c < 3; c++) {
                  for (r = 0; r < 2; r++) {
                    FERM (chi_f0_use, r, c, s) += FERM (temps[mu][1], r, c, s);
                    FERM (chi_f1_use, r, c, s) +=
                      FERM (temps[mu][1] + SPINOR_SIZE, r, c, s);
                  }
                }
              }

              psi += vec_offset;
              chi += vec_offset;

              psi_f1 += vec_offset;
              chi_f1 += vec_offset;
            }


          }                     //vec_ind loop
        }                       //if local comms
      }                         //mu loop

    }                           //parity==cbn
  }
  timer_local.stop();
//  dtime += dclock (true);
//  local += dtime;

  timer_qmp.start();
//  dtime = -dclock ();

  for (int i = 0; i < 8; i++) {
#if 0
    if (!local_comm[i % 4]) {
      QMP_status_t send_status = QMP_wait (Send[i]);
      if (send_status != QMP_SUCCESS)
        QMP_error ("Send failed in wilson_dslash: %s\n",
                   QMP_error_string (send_status));
      QMP_status_t rcv_status = QMP_wait (Recv[i]);
      if (rcv_status != QMP_SUCCESS)
        QMP_error ("Receive failed in wilson_dslash: %s\n",
                   QMP_error_string (rcv_status));
    }
#endif
    ind_nl[i] = 0;
    ind_buf[i] = Recv_buf[i];
  }
  if(n_dir>0){
  QMP_status_t send_status = QMP_wait (multiple);
      if (send_status != QMP_SUCCESS)
        QMP_error ("QMP_multiple failed in wilson_dslash: %s\n", QMP_error_string (send_status));
  }
  timer_qmp.stop();
//  dtime += dclock ();
//  qmp += dtime;

  timer_nl.start();
//  dtime = -dclock ();

  //
  // non-local
  //

#define USE_TEST2
#ifdef USE_TEST2
  index = 0;
  int nl_offset[8];
  for (int i = 0; i < 8; i++) {
    index += num_nl[i];
    nl_offset[i] = index;
  }

#undef NL_OMP
  {
    int i_mu;

    for (i_mu = 0; i_mu < 4; i_mu++) {
#pragma omp parallel for default(shared)
      for (int i = 0; i < num_nl[i_mu]; i++) {
        Float tmp[temp_size];
        Float tmp1[temp_size];

        Float *chi;
        Float *u;
        Float *psi;

        Float *chi_f1;
        Float *u_f1;
        Float *psi_f1;

        chi = chi_p + SPINOR_SIZE * (*(t_ind[i_mu] + i));
        u = u_p + GAUGE_SIZE * (*(u_ind[i_mu] + i));
        psi = ind_buf[i_mu] + i * SPINOR_SIZE * vec_len;

        if (GJP.Gparity ()) {
          int psi_bufoff = num_nl[i_mu] * SPINOR_SIZE * vec_len;        //offset for G-parity f1 pointers in buffer
          chi_f1 = chi_p + SPINOR_SIZE * (*(t_ind[i_mu] + i) + vol);
          u_f1 = u_p + GAUGE_SIZE * (*(u_ind[i_mu] + i) + vol);
          psi_f1 = ind_buf[i_mu] + psi_bufoff + i * SPINOR_SIZE * vec_len;
        }

        int r, c, s;
        for (int vec_ind = 0; vec_ind < vec_len; vec_ind++) {
          MINUSMU (i_mu, u, tmp, tmp1, sdag, psi);
          for (s = 0; s < 4; s++)
            for (c = 0; c < 3; c++)
              for (r = 0; r < 2; r++)
                FERM (chi, r, c, s) += FERM (tmp1, r, c, s);
          psi += SPINOR_SIZE;
          chi += vec_offset;
        }

        if (GJP.Gparity ()) {
          for (int vec_ind = 0; vec_ind < vec_len; vec_ind++) {
            MINUSMU (i_mu, u_f1, tmp, tmp1 + SPINOR_SIZE, sdag, psi_f1);
            for (s = 0; s < 4; s++)
              for (c = 0; c < 3; c++)
                for (r = 0; r < 2; r++)
                  FERM (chi_f1, r, c, s) += FERM (tmp1 + SPINOR_SIZE, r, c, s);
            psi_f1 += SPINOR_SIZE;
            chi_f1 += vec_offset;
          }
        }


        chi = chi_p + SPINOR_SIZE * (*(t_ind[i_mu + 4] + i));
        u = u_p + GAUGE_SIZE * (*(u_ind[i_mu + 4] + i));
        psi = ind_buf[i_mu + 4] + i * SPINOR_SIZE * vec_len;

        if (GJP.Gparity ()) {
          int psi_bufoff = num_nl[i_mu + 4] * SPINOR_SIZE * vec_len;    //offset for G-parity f1 pointers in buffer
          chi_f1 = chi_p + SPINOR_SIZE * (*(t_ind[i_mu + 4] + i) + vol);
          u_f1 = u_p + GAUGE_SIZE * (*(u_ind[i_mu + 4] + i) + vol);
          psi_f1 = ind_buf[i_mu + 4] + psi_bufoff + i * SPINOR_SIZE * vec_len;
        }

        for (int vec_ind = 0; vec_ind < vec_len; vec_ind++) {
          for (s = 0; s < 4; s++)
            for (c = 0; c < 3; c++)
              for (r = 0; r < 2; r++)
                FERM (chi, r, c, s) += FERM (psi, r, c, s);
          psi += SPINOR_SIZE;
          chi += vec_offset;
        }

        if (GJP.Gparity ()) {
          for (int vec_ind = 0; vec_ind < vec_len; vec_ind++) {
            for (s = 0; s < 4; s++)
              for (c = 0; c < 3; c++)
                for (r = 0; r < 2; r++)
                  FERM (chi_f1, r, c, s) += FERM (psi_f1, r, c, s);
            psi_f1 += SPINOR_SIZE;
            chi_f1 += vec_offset;
          }
        }


      }
    }

  }

#endif

//  dtime += dclock (true);
//  nonlocal += dtime;
  timer_nl.stop();

  called++;
#if 0
  if (called % 100 == 0) {
    print_flops ("wilson_dslash_vec()", "local*1000", 0, local);
    print_flops ("wilson_dslash_vec()", "nonlocal*1000", 0, nonlocal);
    print_flops ("wilson_dslash_vec()", "qmp*1000", 0, qmp);
    print_flops ("wilson_dslash_vec()", "setup*1000", 0, setup);
    local = nonlocal = qmp = setup = 0.;
  }
#endif
  DiracOp::CGflops += 1320 * vol * vec_len;

      for (int i = 0; i < n_dir; i++) {
        QMP_free_msgmem (Send_mem[i]);
        QMP_free_msgmem (Recv_mem[i]);
//        QMP_free_msghandle (multiple_p[i*2]);
//        QMP_free_msghandle (multiple_p[i*2+1]);
      }
      if(n_dir>0) QMP_free_msghandle (multiple);
}

CPS_END_NAMESPACE
