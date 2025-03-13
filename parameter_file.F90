module Parameter_File

  public :: parse_parameter_file
  public :: pa_runtime, pa_output, pa_tree

  integer, parameter :: UNIT_STDOUT = 6
  private :: UNIT_STDOUT

  type Parameters_Runtime
    character(len=:), allocatable :: data_path

    ! Initial random seed
    integer :: iseed = -8635
    ! Maximum trees per file
    integer :: max_trees_per_file = 1000
  end type Parameters_Runtime
  type(Parameters_Runtime) :: pa_runtime

  type Parameters_Output
    character(len=:), allocatable :: file_base
    ! Number of levels in the tree
    integer :: nlev = 10
    ! Mass resolution
    real    :: mres = 1.0e+08
  end type Parameters_Output
  type(Parameters_Output) :: pa_output
 
  type Parameters_Cosmology
    real :: omega0  = 0.25
    real :: lambda0 = 0.75
    real :: h0      = 0.73
    real :: omegab  = 0.04

    ! For Eisenstein and Hu CDM transfer function one must specify the CMB temperature
    real :: CMB_T0  = 2.73  ! K
    
  end type Parameters_Cosmology
  type(Parameters_Cosmology) :: pa_cosmo

  type Parameters_Powerspectrum   
    ! itrans=-1  !indicates use transfer function tabulated in file pkinfile
    ! itrans=1   !indicates use BBKS CDM transfer function with specified Gamma and Omega0
    ! itrans=2   !indicates use Bond & Efstathiou CDM transfer function with specified Gamma and Omega0
    ! itrans=3   !indicates use Eisenstein and Hu CDM transfer function with specified Omega0, Omegab and h0
    integer :: itrans = 1

    !Tabulated Millennium Simulation linear P(k)
    character(len=:), allocatable :: pkinfile 

    ! Primordial P(k) parameters (ignored if itrans=-1)
    real :: nspec  = 1.0
    real :: dndlnk = 0.0
    real :: kref   = 1.0

    ! Power spectrum amplitude set regardless of other parameters
    real :: sigma8  = 0.9  

    ! Computed at runtime
    ! Omega_m x h ignoring effect of baryons
    real :: gamma   
  end type Parameters_Powerspectrum
  type(Parameters_Powerspectrum) :: pa_powerspec

  type Parameters_Tree
    ! Parameters used to modify the merger rate used in split.F90
    ! to be slightly different to the standard Press-Schechter formula
    !
    ! The modify factor is
    ! G0 [sigma(m1)/sigma(m2)]^gamma_1 [w/sigma(m2)]^gamma_2

    ! Parameters of the Merger Tree Algorithm as defined in Parkinson, Cole
    ! and Helly (2007arXiv0708.138 version 3 and in MNRAS paper).

    ! These values supercede the values given in the original astro-ph posting
    ! due to a small error in the code being identified. In this version of the
    ! code the error has been rectified and the fits redone. Using this code
    ! and these new parameters will produce near identical results to the old
    ! code with the old parameters.
    real :: G0 = 0.57
    real :: gamma_1 = 0.38
    real :: gamma_2 = -0.01
    real :: eps1 = 0.1
    real :: eps2 = 0.1
  end type Parameters_Tree
  type(Parameters_Tree) :: pa_tree

contains

  subroutine parse_parameter_file(file_name, dump_parameters_unit)
    use, intrinsic :: iso_fortran_env, only: stderr => error_unit
    use :: tomlf
    implicit none
    
    character(len=*), intent(in) :: file_name
    integer, intent(in), optional :: dump_parameters_unit
    
    logical :: dump_parameters
    integer :: dpu

    integer                       :: fu, rc
    logical                       :: file_exists
    type(toml_table), allocatable :: table
    type(toml_table), pointer     :: child

    dump_parameters = present(dump_parameters_unit)
    dpu = merge(dump_parameters_unit, 6, dump_parameters)

    inquire (file=file_name, exist=file_exists)
    if (.not. file_exists) then
      write (stderr, '("Error: Parameter file ", a, " not found")') file_name
      stop
    end if

    open(action='read', file=file_name, iostat=rc, newunit=fu)

    if (rc /= 0) then
      write (stderr, '("Error: Reading TOML file ", a, " failed")') file_name
      stop
    end if

    call toml_parse(table, fu)
    close(fu)

    if (.not. allocated(table)) then
      write (stderr, '("Error: Parsing failed (allocating table)")')
      stop
    end if

    ! Get [runtime] section.
    call get_value(table, 'runtime', child, requested=.false.)
    runtime_parameters: if (associated(child)) then
      call get_value(child, 'data_path', pa_runtime%data_path, './')
      call get_value(child, 'iseed',     pa_runtime%iseed, -8365)
      call get_value(child, 'max_trees_per_file',  pa_runtime%max_trees_per_file, 1000)
    end if runtime_parameters
    
    if (dump_parameters) then
      write(*,*)
      write(*,*) '[runtime]'
      call print_kv('data_path', './')
      call print_kv('iseed', pa_runtime%iseed)
      call print_kv('max_trees_per_file', pa_runtime%max_trees_per_file)
    end if

    ! Get [output] section.
    call get_value(table, 'output', child, requested=.false.)
    output_parameters: if (associated(child)) then
      call get_value(child, 'file_base', pa_output%file_base, './output_tree')
      call get_value(child, 'nlev', pa_output%nlev, 10)
      call get_value(child, 'mres', pa_output%mres, 1.0e+8)
    end if output_parameters

    if (dump_parameters) then
      write(*,*)
      write(*,*) '[output]'
      call print_kv('file_base', './output_tree')
      call print_kv('nlev', pa_output%nlev)
      call print_kv('mres', pa_output%mres)
    end if

    ! Get [cosmology] section
    call get_value(table, 'cosmology', child, requested=.false.)
    cosmo_parameters: if (associated(child)) then
      call get_value(child, 'omega0',  pa_cosmo%omega0)
      call get_value(child, 'lambda0', pa_cosmo%lambda0)
      call get_value(child, 'h0',      pa_cosmo%h0)
      call get_value(child, 'omegab',  pa_cosmo%omegab)
      call get_value(child, 'cmb_T0',  pa_cosmo%CMB_T0) 
    end if cosmo_parameters

    if (dump_parameters) then
      write(*,*)
      write(*,*) '[cosmology]'
      call print_kv('omega0',  pa_cosmo%omega0)
      call print_kv('lambda0', pa_cosmo%lambda0)
      call print_kv('h0',      pa_cosmo%h0)
      call print_kv('omegab',  pa_cosmo%omegab)
      call print_kv('CMB_T0',  pa_cosmo%CMB_T0)
    end if

    ! Get [powerspec] section
    call get_value(table, 'powerspec', child, requested=.false.)
    powerspec_parameters: if (associated(child)) then
      call get_value(child, 'pkinfile', pa_powerspec%pkinfile, 'pk_Mill.dat')
    
      call get_value(child, 'itrans', pa_powerspec%itrans)
      call get_value(child, 'nspec',  pa_powerspec%nspec)
      call get_value(child, 'dndlnk', pa_powerspec%dndlnk)
      call get_value(child, 'kref',   pa_powerspec%kref)
      call get_value(child, 'sigma8', pa_powerspec%sigma8) 
    end if powerspec_parameters

    ! Compute powerspec gamma
    pa_powerspec%gamma = pa_cosmo%omega0*pa_cosmo%h0 

    if (dump_parameters) then
      write(*,*)
      write(*,*) '[powerspec]'
      call print_kv('itrans',   pa_powerspec%itrans)
      call print_kv('pkinfile', pa_powerspec%pkinfile)
      call print_kv('nspec',    pa_powerspec%nspec)
      call print_kv('dndlnk',   pa_powerspec%dndlnk)
      call print_kv('kref',     pa_powerspec%kref)
      call print_kv('sigma8',   pa_powerspec%sigma8)
    end if 

    ! Get [tree] section.
    call get_value(table, 'tree', child, requested=.false.)
    tree_parameters: if (associated(child)) then
      call get_value(child, 'G0',      pa_tree%G0,       0.57)
      call get_value(child, 'gamma_1', pa_tree%gamma_1,  0.38)
      call get_value(child, 'gamma_2', pa_tree%gamma_2, -0.01)
      call get_value(child, 'eps1',    pa_tree%eps1,     0.1)
      call get_value(child, 'eps2',    pa_tree%eps2,     0.1)
    end if tree_parameters

    if (dump_parameters) then
      write(*,*)
      write(*,*) '[tree]'
      call print_kv('G0', pa_tree%G0)
      call print_kv('gamma_1', pa_tree%gamma_1)
      call print_kv('gamma_2', pa_tree%gamma_2)
      call print_kv('eps1', pa_tree%eps1)
      call print_kv('eps2', pa_tree%eps2)
    end if

    write(*,*)
  end subroutine parse_parameter_file
  
  ! ############################################################ 
  subroutine print_kv(k, v, unit)
    implicit none

    character(len=*), intent(in) :: k
    class(*), intent(in) :: v
    integer, intent(in), optional :: unit
    
    integer :: output_unit
  
    ! Use provided unit, else default to stdout (unit 6)
    output_unit = merge(unit, 6, present(unit))

    select type(v)
    type is (integer)
      write(output_unit, '(1X, A, A, I0)') trim(k), ' = ', v
    type is (real)
      write(output_unit, '(1X, A, A, G0)') trim(k), ' = ', v
    type is (character(len=*))
      write(output_unit, '(1X, A, A, A, A)') trim(k), ' = "', trim(v), '"'
    class default
      write(output_unit, *) 'Unsupported type'
    end select
  end subroutine print_kv

end module Parameter_File
