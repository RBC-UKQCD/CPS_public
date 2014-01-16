#line 1555 "dwf.nw"
struct neighbor {
    int              size;            /* size of site table */
    int              inside_size;     /* number of inside sites */
    int              boundary_size;   /* number of boundary sites */
    int              snd_size[2*DIM]; /* size of send buffers in 8 dirs */
    int              rcv_size[2*DIM]; /* size of receive buffers */
    int             *snd[2*DIM];      /* i->x translation for send buffers */
    int             *inside;          /* i->x translation for inside sites */
    struct boundary *boundary;        /* i->x,mask translation for boundary */
    struct site     *site;            /* x->site translation for sites */
    vHalfFermion    *snd_buf[2*DIM];  /* Send buffers */
    vHalfFermion    *rcv_buf[2*DIM];  /* Receive buffers */

    int              qmp_size[4*DIM]; /* sizes of QMP buffers */
    void            *qmp_xbuf[4*DIM]; /* QMP snd/rcv buffer addresses */
    vHalfFermion    *qmp_buf[4*DIM];  /* send and receive buffers for QMP */
    QMP_msgmem_t     qmp_mm[4*DIM];   /* msgmem's for send and receive */
    int              Nx;              /* number of msegs */

    int              qmp_smask;       /* send flags for qmp_sh[] */
    QMP_msghandle_t  qmp_handle;      /* common send & receive handle */
};
#line 1581 "dwf.nw"
struct boundary {
    int   index;   /* x-index of this boundary site */
    int   mask;    /* bitmask of the borders */
};
#line 1590 "dwf.nw"
struct site {
  int Uup;           /* up-links are Uup, Uup+1, Uup+2, Uup+3 */
  int Udown[DIM];    /* four down-links */
  int F[2*DIM];      /* eight neighboring fermions on the other sublattice */
};
#line 4180 "dwf.nw"
extern struct neighbor even_odd;
extern struct neighbor odd_even;
