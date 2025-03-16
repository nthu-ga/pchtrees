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
    ! use :: tomlf
    use TinyTOML
    implicit none
 
    character(len=*), intent(in) :: file_name
    integer, intent(in), optional :: dump_parameters_unit
    
    logical :: dump_parameters
    integer :: dpu

    integer                       :: fu, rc
    logical                       :: file_exists

    type(toml_object) :: toml_content
    type(toml_object) :: section

    dump_parameters = present(dump_parameters_unit)
    dpu = merge(dump_parameters_unit, 6, dump_parameters)

    inquire (file=file_name, exist=file_exists)
    if (.not. file_exists) then
      write (stderr, '("Error: Parameter file ", a, " not found")') file_name
      stop
    end if

    toml_content = parse_file(file_name)

    ! Get [runtime] section.
    section = toml_content%get("runtime")
    call read_value(section%get("data_path",          error=.false.), pa_runtime%data_path, default='./')
    call read_value(section%get("iseed",              error=.false.), pa_runtime%iseed, default=-8365)
    call read_value(section%get("max_trees_per_file", error=.false.), pa_runtime%max_trees_per_file, default=1000)

    if (dump_parameters) then
      write(*,*)
      write(*,*) '[runtime]'
      call print_kv('data_path', pa_runtime%data_path)
      call print_kv('iseed', pa_runtime%iseed)
      call print_kv('max_trees_per_file', pa_runtime%max_trees_per_file)
    end if

    ! Get [output] section.
    section = toml_content%get("output")
    call read_value(section%get('file_base', error=.false.), pa_output%file_base, default='./output_tree')
    call read_value(section%get('nlev', error=.false.), pa_output%nlev, default=10)
    call read_value(section%get('mres', error=.false.), pa_output%mres, default=1.0e+8)

    if (dump_parameters) then
      write(*,*)
      write(*,*) '[output]'
      call print_kv('file_base', pa_output%file_base)
      call print_kv('nlev', pa_output%nlev)
      call print_kv('mres', pa_output%mres)
    end if

    ! Get [cosmology] section
    section = toml_content%get("cosmology")
    call read_value(section%get('omega0',  error=.false.),  pa_cosmo%omega0)
    call read_value(section%get('lambda0', error=.false.), pa_cosmo%lambda0)
    call read_value(section%get('h0',      error=.false.),      pa_cosmo%h0)
    call read_value(section%get('omegab',  error=.false.),  pa_cosmo%omegab)
    call read_value(section%get('cmb_T0',  error=.false.),  pa_cosmo%CMB_T0) 

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
    section = toml_content%get("powerspec")
    call read_value(section%get('itrans', error=.false.),   pa_powerspec%itrans)
    call read_value(section%get('pkinfile', error=.false.), pa_powerspec%pkinfile, default='pk_Mill.dat')
    call read_value(section%get('nspec', error=.false.),    pa_powerspec%nspec)
    call read_value(section%get('dndlnk', error=.false.),   pa_powerspec%dndlnk)
    call read_value(section%get('kref', error=.false.),     pa_powerspec%kref)
    call read_value(section%get('sigma8', error=.false.),   pa_powerspec%sigma8) 

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
    section = toml_content%get("tree")
    call read_value(section%get('G0', error=.false.),      pa_tree%G0,     default= 0.57)
    call read_value(section%get('gamma_1', error=.false.), pa_tree%gamma_1,default= 0.38)
    call read_value(section%get('gamma_2', error=.false.), pa_tree%gamma_2,default=-0.01)
    call read_value(section%get('eps1', error=.false.),    pa_tree%eps1,   default= 0.1)
    call read_value(section%get('eps2', error=.false.),    pa_tree%eps2,   default= 0.1)

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
