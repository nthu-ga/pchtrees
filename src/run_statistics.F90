module Run_Statistics
      ! Variables holding information on run-time statistics
      integer NI
      parameter (NI=60) ! Maximum number of iterations.
      !
      ! Integers
      integer icall(4),icall_sf,max_inode,nc1,nc2,ncall,nfail,nhist(20),niter_tab(NI),n_nelder_mead
#ifdef DEBUG
      data icall /0,0,0,0/
      data ncall /0/
      data nc1 /0/
      data nc2 /0/
#endif
end module Run_Statistics
