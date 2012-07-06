/*! BfmArg is used to pass parameters to bfm.

    Note: the corresponding class in bfm is called bfmarg.

    This one is the equivalence in CPS, but the two are not entirely
    the same.

    Check the bfm package to find what each parameter does. */

class BfmArg { 
      BfmSolverType solver; 
      /*! IMPORTANT: BfmSolverType is not the same as the BfmSolver in the bfm package. BfmSolverType is defined in enum.x. Basically it adds a BFM_ prefix to the corresponding names in the BfmSolver enum. */

      int solveMobiusDminus;
      int Ls;
      Float M5;
      Float mass;
      Float twistedmass;
      Float Csw;
      
      /* node_latt, local_comm and ncoor are omitted */

      Float zolo_delta;
      Float zolo_lo;
      Float zolo_hi;
      
      int list_engine;
      int list_length;
      int powerloop;
      int time_report_iter;
      int max_iter;
      Float residual;
      int precon_5d;

      /* watchfile is omitted, don't know what it is used for. */
      int reproduce;
      int reproduce_checksum;
      int reproduce_master_check;
      int threads;
      int verbose;
      int onepluskappanorm;
      Float mobius_scale;
};
