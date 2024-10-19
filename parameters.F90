! Prompt for and read in all model parameters.
! Many are then distributed to the subroutines that need to know them
! via a set of modules.

module Power_Spectrum_Parameters
  ! Variables used to hold properties of the power spectrum
  !
  ! Array dimensions
  integer Transfer_Function_Table_N_Max
  parameter (Transfer_Function_Table_N_Max=1800)
  !
  ! Integers
  integer igwave,ireset,itrans,nktab,NKTABMAX,NSPL,Trans_Func_Table_N_Points
  !
  ! Array dimensions
  parameter(NKTABMAX=1000)
  !
  ! Floats
  real dndlnk,gamma,kref,lnktab(NKTABMAX),lnpktab(NKTABMAX),mwdm,nspec,sigma8,scla,sclm
  real Transfer_Function_Table_lnk(Transfer_Function_Table_N_Max),Transfer_Function_Table_lnTk(Transfer_Function_Table_N_Max)
  !
  ! Logicals
  logical WDMrun
  !
  ! Characters
  character pkinfile*1024,splinefile*220,tffile*1024
  !
  ! Parameters
  parameter (NSPL=200)
end module Power_Spectrum_Parameters

module Time_Parameters
  ! Parameters used in making binary splits in the merger tree
  real eps1,eps2
  integer istep
end module Time_Parameters
