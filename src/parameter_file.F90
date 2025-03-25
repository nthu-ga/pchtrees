module Parameter_File

  public :: parse_parameter_file
  public :: pa_runtime, pa_output, pa_tree

  integer, parameter :: UNIT_STDOUT = 6
  private :: UNIT_STDOUT

  ! Default parameters 

  integer, parameter :: PA_RUNTIME_ISEED_DEF = -8635
  integer, parameter :: PA_OUTPUT_MAX_TREES_PER_FILE_DEF = 1000

  integer, parameter :: PA_OUTPUT_NLEV_DEF = 10
  real,    parameter :: PA_OUTPUT_MRES_DEF = 1.0e+08

  real,    parameter :: PA_COSMOLOGY_OMEGA0_DEF  = 0.25
  real,    parameter :: PA_COSMOLOGY_LAMBDA0_DEF = 0.75
  real,    parameter :: PA_COSMOLOGY_H0_DEF      = 0.73
  real,    parameter :: PA_COSMOLOGY_OMEGAB_DEF  = 0.04
  real,    parameter :: PA_COSMOLOGY_CMB_T0_DEF  = 2.73  ! K

  integer, parameter :: PA_POWERSPEC_ITRANS_DEF = 1
  real,    parameter :: PA_POWERSPEC_NSPEC_DEF  = 1.0
  real,    parameter :: PA_POWERSPEC_DNDLNK_DEF = 0.0
  real,    parameter :: PA_POWERSPEC_KREF_DEF   = 1.0
  real,    parameter :: PA_POWERSPEC_SIGMA8_DEF = 0.9  

  real,    parameter :: PA_TREE_G0_DEF      = 0.57
  real,    parameter :: PA_TREE_GAMMA_1_DEF = 0.38
  real,    parameter :: PA_TREE_GAMMA_2_DEF = -0.01
  real,    parameter :: PA_TREE_EPS1_DEF    = 0.1
  real,    parameter :: PA_TREE_EPS2_DEF    = 0.1

  ! Enum for output format
  integer, parameter :: OUTPUT_HDF5 = 1
  integer, parameter :: OUTPUT_JET  = 2

  ! Classes and global isntances for the parameter file sections
  type Parameters_Runtime
    character(len=:), allocatable :: data_path
    ! Initial random seed
    integer :: iseed = PA_RUNTIME_ISEED_DEF
      end type Parameters_Runtime
  type(Parameters_Runtime) :: pa_runtime

  type Parameters_Output
    ! Type of output
    integer :: output_format
    ! Base of file name (e.g. <file_base>.nnn.hdf5)
    character(len=:), allocatable :: file_base
    ! Extension
    character(len=:), allocatable :: file_ext
    ! Number of levels in the tree
    integer :: nlev = PA_OUTPUT_NLEV_DEF 
    ! Mass resolution
    real    :: mres = PA_OUTPUT_MRES_DEF
    ! Output time list
    character(len=:), allocatable :: aexp_list
    logical :: have_aexp_list
    ! Maximum trees per file
    integer :: max_trees_per_file = PA_OUTPUT_MAX_TREES_PER_FILE_DEF
  end type Parameters_Output
  type(Parameters_Output) :: pa_output
 
  type Parameters_Cosmology
    real :: omega0  = PA_COSMOLOGY_OMEGA0_DEF
    real :: lambda0 = PA_COSMOLOGY_LAMBDA0_DEF
    real :: h0      = PA_COSMOLOGY_H0_DEF
    real :: omegab  = PA_COSMOLOGY_OMEGAB_DEF

    ! For Eisenstein and Hu CDM transfer function one must specify the CMB temperature
    real :: CMB_T0  = PA_COSMOLOGY_CMB_T0_DEF  ! K
  end type Parameters_Cosmology
  type(Parameters_Cosmology) :: pa_cosmo

  type Parameters_Powerspectrum   
    ! itrans=-1  !indicates use transfer function tabulated in file pkinfile
    ! itrans=1   !indicates use BBKS CDM transfer function with specified Gamma and Omega0
    ! itrans=2   !indicates use Bond & Efstathiou CDM transfer function with specified Gamma and Omega0
    ! itrans=3   !indicates use Eisenstein and Hu CDM transfer function with specified Omega0, Omegab and h0
    integer :: itrans = PA_POWERSPEC_ITRANS_DEF

    !Tabulated Millennium Simulation linear P(k)
    character(len=:), allocatable :: pkinfile 

    ! Primordial P(k) parameters (ignored if itrans=-1)
    real :: nspec  = PA_POWERSPEC_NSPEC_DEF
    real :: dndlnk = PA_POWERSPEC_DNDLNK_DEF
    real :: kref   = PA_POWERSPEC_KREF_DEF

    ! Power spectrum amplitude set regardless of other parameters
    real :: sigma8 = PA_POWERSPEC_SIGMA8_DEF

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

  subroutine parse_parameter_file(file_name_in, dump_parameters_unit)
    use, intrinsic :: iso_fortran_env, only: stderr => error_unit
    ! use :: tomlf
    use TinyTOML
    implicit none
 
    character(len=*), intent(in) :: file_name_in
    character(len=100) :: file_name
    integer, intent(in), optional :: dump_parameters_unit
    
    logical :: dump_parameters
    integer :: dpu

    ! integer :: fu, rc
    logical :: file_exists, dump_only=.false.

    type(toml_object) :: toml_content
    type(toml_object) :: section
    type(toml_object) :: temp_keyval

    character(len=100) :: temp_filename
    integer :: temp_unit
    integer :: io_err

    dump_parameters = present(dump_parameters_unit)
    dpu = merge(dump_parameters_unit, 6, dump_parameters)

    inquire (file=file_name_in, exist=file_exists)
    if (.not. file_exists) then
      if (dump_parameters) then
        if (dump_parameters_unit.eq.6) then
          dump_only = .true.
        endif
      else
        write (stderr, '("Error: Parameter file ", a, " not found")') file_name_in
        stop
      endif
    end if

    ! An ugly hack to write a dummy pf
    if (dump_only) then
      ! Generate a unique temporary file name
      temp_filename = 'tmpfile.dat'
      open(newunit=temp_unit, file=temp_filename, status='new', action='readwrite')
      write(temp_unit, '(A)') '[runtime]' 
      write(temp_unit, '(A)') '[output]' 
      write(temp_unit, '(A)') '[cosmology]' 
      write(temp_unit, '(A)') '[powerspec]' 
      write(temp_unit, '(A)') '[tree]' 
      file_name = trim(temp_filename)
      close(temp_unit)

      write(dpu,*)
      write(dpu,*) '# This is an example pchtrees paramter file.'
      write(dpu,*) '# Some of the options will need to be supplied by you,'
      write(dpu,*) '# e.g. data_path, file_base.'
    else
      file_name = trim(file_name_in)
    endif

    ! Parse the TOML file
    toml_content = parse_file(file_name)

    if (dump_only) then
      io_err = unlink(temp_filename)
    endif

    ! Get [runtime] section.
    section = toml_content%get("runtime")
    call read_value(section%get("data_path", error=.false.), &
      & pa_runtime%data_path, default='./data')
    call read_value(section%get("iseed", error=.false.), &
      & pa_runtime%iseed, default=PA_RUNTIME_ISEED_DEF)
    
    if (dump_parameters) then
      write(*,*)
      write(*,*) '[runtime]'
      call print_kv('data_path', pa_runtime%data_path)
      call print_kv('iseed', pa_runtime%iseed)
    end if

    ! Get [output] section.
    section = toml_content%get("output")
    call read_value(section%get('file_base', error=.false.), & 
      & pa_output%file_base, default='./output_tree')
    call read_value(section%get('nlev', error=.false.), & 
      & pa_output%nlev, default=PA_OUTPUT_NLEV_DEF)
    call read_value(section%get('mres', error=.false.), & 
      & pa_output%mres, default=PA_OUTPUT_MRES_DEF)
    call read_value(section%get("max_trees_per_file", error=.false.), &
      & pa_output%max_trees_per_file, default=PA_OUTPUT_MAX_TREES_PER_FILE_DEF)

    ! Output format
    temp_keyval = section%get('output_format', error=.false.)
    select case(temp_keyval%value)
    case ('hdf5')
      pa_output%output_format = OUTPUT_HDF5
      pa_output%file_ext = 'hdf5'
    case ('jet')
      pa_output%output_format = OUTPUT_JET
      pa_output%file_ext = 'jet.txt'
    case default
      if (dump_parameters) then
        ! Set a dummy value; we want to write a default string as output
        ! and then quit, but we need to keep going here to process the
        ! other parameters.
        pa_output%output_format = -1
      else
        write(*,*) 'Warning: invalid output format:', temp_keyval%value
        stop
      end if
    end select

    ! Optional aexp list
    temp_keyval = section%get('aexp_list', error=.false.)
    pa_output%have_aexp_list = .false.
    select case(temp_keyval%error_code)
    case (KEY_NOT_FOUND)
      ! No axep list, ok
      ! Not really needed, but no "pass" in Fortran...
      pa_output%have_aexp_list = .false.
    case (SUCCESS)
      call read_value(temp_keyval, pa_output%aexp_list)
      ! Only sanity check is that an empty value is
      ! counted as no value
      if (len(trim(pa_output%aexp_list)).gt.0) then
        pa_output%have_aexp_list = .true.
      else
        pa_output%have_aexp_list = .false.
      endif
    case default
      write(*,*) 'Bad news!'
      stop
    end select

    if (dump_parameters) then
      write(*,*)
      write(*,*) '[output]'
      if (pa_output%output_format.ge.0) then
        call print_kv('output_format', pa_output%output_format)
      else
        ! A value of -1 indicates we're writing default arguments
        ! so we should dump a default string, not the enum value.
        call print_kv('output_format', 'jet')
      endif
      call print_kv('file_base', pa_output%file_base)
      call print_kv('nlev', pa_output%nlev)
      call print_kv('mres', pa_output%mres)
      if (pa_output%have_aexp_list) then
        call print_kv('aexp_list', pa_output%aexp_list)
      endif
    end if

    ! Get [cosmology] section
    section = toml_content%get("cosmology")
    call read_value(section%get('omega0',  error=.false.), &
      &   pa_cosmo%omega0, default=PA_COSMOLOGY_OMEGA0_DEF)
    call read_value(section%get('lambda0', error=.false.), &
      &  pa_cosmo%lambda0, default=PA_COSMOLOGY_LAMBDA0_DEF)
    call read_value(section%get('h0',      error=.false.), &
      &   pa_cosmo%h0, default=PA_COSMOLOGY_H0_DEF)
    call read_value(section%get('omegab',  error=.false.), &
      &   pa_cosmo%omegab, default=PA_COSMOLOGY_OMEGAB_DEF)
    call read_value(section%get('cmb_T0',  error=.false.), &
      &   pa_cosmo%CMB_T0, default=PA_COSMOLOGY_CMB_T0_DEF) 

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
    call read_value(section%get('itrans', error=.false.),   &
      & pa_powerspec%itrans, default=PA_POWERSPEC_ITRANS_DEF)
    call read_value(section%get('pkinfile', error=.false.), &
      & pa_powerspec%pkinfile, default='pk_Mill.dat')
    call read_value(section%get('nspec', error=.false.),    &
      & pa_powerspec%nspec,  default=PA_POWERSPEC_NSPEC_DEF)
    call read_value(section%get('dndlnk', error=.false.),   &
      & pa_powerspec%dndlnk, default=PA_POWERSPEC_DNDLNK_DEF)
    call read_value(section%get('kref', error=.false.),     &
      & pa_powerspec%kref,   default=PA_POWERSPEC_KREF_DEF)
    call read_value(section%get('sigma8', error=.false.),   &
      & pa_powerspec%sigma8, default=PA_POWERSPEC_SIGMA8_DEF) 

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
    call read_value(section%get('G0', error=.false.),     &
      & pa_tree%G0,     default=PA_TREE_G0_DEF)
    call read_value(section%get('gamma_1', error=.false.),&
      & pa_tree%gamma_1,default=PA_TREE_GAMMA_1_DEF)
    call read_value(section%get('gamma_2', error=.false.),&
      & pa_tree%gamma_2,default=PA_TREE_GAMMA_2_DEF)
    call read_value(section%get('eps1', error=.false.),   &
      & pa_tree%eps1,   default=PA_TREE_EPS1_DEF)
    call read_value(section%get('eps2', error=.false.),   &
      & pa_tree%eps2,   default=PA_TREE_EPS2_DEF)

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
