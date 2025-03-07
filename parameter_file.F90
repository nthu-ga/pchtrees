module Parameter_File

  public :: parse_parameter_file
  public :: pa_runtime, pa_output, pa_tree

  type Parameters_Runtime
    character(len=:), allocatable :: data_path
  end type Parameters_Runtime
  type(Parameters_Runtime) :: pa_runtime

  type Parameters_Output
    character(len=:), allocatable :: file_path
    ! Number of levels in the tree
    integer :: nlev = 10
    ! Mass resolution
    real    :: mres = 1.0e+08
    ! Initial random seed
    integer :: iseed = -8635

  end type Parameters_Output
  type(Parameters_Output) :: pa_output
 
  type Parameters_Cosmology
    real :: omega0  = 0.25
    real :: lambda0 = 0.75
    real :: h0      = 0.73
    real :: omegab  = 0.04
    real :: sigma8  = 0.9   !power spectrum amplitude set regardless of other parameters

    ! Computed at runtime
    ! Omega_m x h ignoring effect of baryons
    real :: Gamma   
  end type Parameters_Cosmology
  type(Parameters_Cosmology) :: pa_cosmo

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

  subroutine parse_parameter_file(file_name)
    use, intrinsic :: iso_fortran_env, only: stderr => error_unit
    use :: tomlf
    implicit none
    ! character(len=*), parameter :: FILE_NAME = 'sample.toml'
    character(len=*), intent(in) :: file_name

    integer                       :: fu, rc
    logical                       :: file_exists
    type(toml_table), allocatable :: table
    type(toml_table), pointer     :: child

    inquire (file=file_name, exist=file_exists)

    if (.not. file_exists) then
      write (stderr, '("Error: TOML file ", a, " not found")') file_name
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
      print '(2a)', 'Data path: ', pa_runtime%data_path
    end if runtime_parameters

    ! Get [output] section.
    call get_value(table, 'runtime', child, requested=.false.)
    output_parameters: if (associated(child)) then
      call get_value(child, 'file_path', pa_output%file_path, './output_tree.hdf5')
      call get_value(child, 'nlev', pa_output%nlev, 10)
      call get_value(child, 'mres', pa_output%mres, 1.0e+8)
      call get_value(child, 'iseed', pa_output%iseed, -8365)

      print '(2a)',      'Output file path: ', pa_output%file_path
      print '(a,i10)',   'Tree levels: ', pa_output%nlev
      print '(a,e10.4)', 'Mass resolution: ', pa_output%mres
      print '(a,i10)',   'Random seeds: ', pa_output%iseed
    end if output_parameters

    ! Get [cosmology] section
    call get_value(table, 'cosmology', child, requested=.false.)
    cosmo_parameters: if (associated(child)) then
      call get_value(child, 'omega0',  pa_cosmo%omega0)
      call get_value(child, 'lambda0', pa_cosmo%lambda0)
      call get_value(child, 'h0',      pa_cosmo%h0)
      call get_value(child, 'omegab',  pa_cosmo%omegab)
      call get_value(child, 'sigma8',  pa_cosmo%sigma8) 
    end if cosmo_parameters
    ! Compute Gamma
    pa_cosmo%Gamma = pa_cosmo%omega0*pa_cosmo%h0 

    ! Get [tree] section.
    call get_value(table, 'tree', child, requested=.false.)
    tree_parameters: if (associated(child)) then
      call get_value(child, 'G0',      pa_tree%G0,       0.57)
      call get_value(child, 'gamma_1', pa_tree%gamma_1,  0.38)
      call get_value(child, 'gamma_2', pa_tree%gamma_2, -0.01)
      call get_value(child, 'eps1',    pa_tree%eps1,     0.1)
      call get_value(child, 'eps2',    pa_tree%eps2,     0.1)
    end if tree_parameters

  end subroutine parse_parameter_file

end module Parameter_File
