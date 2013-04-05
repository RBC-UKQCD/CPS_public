2012-02-14

Command line arguments are 
  ./cps-lanczos.x  ./ $LMA_SHIFTS $do_cg $do_eigv_read \
  -qmp-geom $geom

*LMA_SHIFTs,  set 0 if you don't want to do LMA or AMA
 For example, on 16^3 x 32 lattice,  and LMA_SHIFTS="2 2 2 4" then  source points are
 (0,0,0,0), (8,0,0,0) ... (0,0,0,8), ... (8,8,8,24), 
  in total 2x2x2x4=32 source points

* do_cg :  if it will do the deflated CG after lanczos to compute several hadron 2pt or not
* do_eigv_read : if it will read the eigenvector files and cache in memeory before computing CG, LMA and AMA.


The executable will do the following forur works :

1.  Lanczos

   lanczos_arg.vml  : specify the parameters for lanczos

      RitzMat_lanczos  : matrix used in the polynomial, which is used in the lanczos steps
      RitzMat_convcheck : matrix used in the convergence check and the eigenvalue computation
      mass 
      nk_lanczos_vectors : number of wanted eigenvectors
      np_lanczos_vectors : number of unwanted eigenvectors
      eigen_shift : shift amount of targetting spectrum. e.g if the eigen_shift = 0.01, then eigenvalues closest to 0.01, rather than 0.0,  will be computed.  If it's 0, it will be slightly faster.
      stop_residual : convergence criteria for norm, not norm squared, of residue vector in eigen equation,   e.g.  1e-10
      maxiters : max number of restarting
      results  results file, to report the number of iteration etc.
      file : directory name to store the eigenvalue and eigenvectors
     matpoly_arg : place holder, leave it to zero.

  cheby_arg.vml   : specify the parameters for polynomial of matrix in lanczos

      Npol : degree of polynomial used in the lanczos steps
        params[0] :  this value should be slightly above the  largest magnitude of wanted eigenvalues
                     e.g. when nk_lanczos_vector = 50,  and expected  lambda(49) == 0.07, then optimal params[0] would be 0.1 or so
        params[1] : must be larger than maximum magnitude of  eigenvalues

   Lanczos part will be skipped  if  nk_lanczos_vectors == 0  in lanczos_arg.vml. 

2. Full CG propagator (if command line argment do_cg is non zero)

   qpropw_arg.vml   (point source only for now, but should be easy to extend)
      mass  : mass of quark propagator (don't confuse with mass in lanczos_arg, which is used for eigenvector computation)
    
     fname_eigen :  directory from which eigenvectors will be read for deflation / LMA / AMA
     neig  : number of eigenvectors to be used  (again don't mixed up with lanczos_arg's  nk_lanczos_vectors) 
     eigen_shift :  this is a place-holder, rather than user's input, so leave it to zero.
     double ama_stop_rsd=-1.0


3.  LMA  (if  LMA_SHIFT[0] *LMA_SHIFT[1] *LMA_SHIFT[2] *LMA_SHIFT[3]  > 0 )

    do the LMA, the source location will be shifted from the source location specified by qpropw_arg's x,y,z,t.
    Amount of shift will be equally spaced in x,y,z,t direction from the original source point.
    number of source point in each of four direction is specified  the command line arguments LMA_SHIFT[0..3],
       For example, on 16^3 x 32 lattice,  and LMA_SHIFTS="2 2 2 4" then  source points are
         (0,0,0,0), (8,0,0,0) ... (0,0,0,8), ... (8,8,8,24), 
        in total 2x2x2x4=32 source points


4.  AMA  (if  ama_stop_rsd > 0  in qpropw_arg.cg )

    It will do the AMA with same source points described in LMA section.




 
