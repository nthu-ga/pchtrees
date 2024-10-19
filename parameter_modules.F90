module Cosmological_Parameters
  !
  ! Floats
  real h0,lambda0,omega0,omegab,CMB_T0
end module Cosmological_Parameters

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

module Halo_Mass_Function
  !
  ! Integers
  integer massfun
  integer, parameter :: Halo_Mass_Function_PS=0,Halo_Mass_Function_SMT=1,Halo_Mass_Function_J2000=2
end module Halo_Mass_Function
