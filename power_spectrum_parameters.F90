module Power_Spectrum_Parameters
  implicit none
  ! Variables used to hold properties of the power spectrum
  
  ! Array dimensions
  integer, parameter :: Transfer_Function_Table_N_Max = 1800
  
  ! Integers
  integer :: igwave,ireset,itrans,nktab,Trans_Func_Table_N_Points
  
  ! Array dimensions
  integer, parameter :: NKTABMAX=1000
  integer, parameter :: NSPL=200

  ! Floats
  real :: gamma
  real :: dndlnk,kref,lnktab(NKTABMAX),lnpktab(NKTABMAX),mwdm,nspec,sigma8,scla,sclm
  real :: Transfer_Function_Table_lnk(Transfer_Function_Table_N_Max),Transfer_Function_Table_lnTk(Transfer_Function_Table_N_Max)
  
  ! Logicals
  logical :: WDMrun
  
  ! Characters
  character(len=1024) :: pkinfile
  character(len=220)  :: splinefile
  character(len=1024) :: tffile
  
end module Power_Spectrum_Parameters


