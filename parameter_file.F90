! example.f90
module Parameter_File

  public :: parse_parameter_file
  public :: runtime, pa_tree

  type Parameters_Runtime
    character(len=:), allocatable :: data_path
  end type Parameters_Runtime

  type(Parameters_Runtime) :: runtime

  type Parameters_Tree
   real :: G0 = 0.57
   real :: gamma_1 = 0.38
   real :: gamma_2 = -0.01
   real :: eps1 = 0.1
   real :: eps2 = 0.1
  end type Parameters_Tree

  type(Parameters_Tree) :: pa_tree

  contains

subroutine parse_parameter_file()
    use, intrinsic :: iso_fortran_env, only: stderr => error_unit
    use :: tomlf
    implicit none
    character(len=*), parameter :: FILE_NAME = 'sample.toml'

    integer                       :: fu, rc
    logical                       :: file_exists
    type(toml_table), allocatable :: table
    type(toml_table), pointer     :: child

    inquire (file=FILE_NAME, exist=file_exists)

    if (.not. file_exists) then
      write (stderr, '("Error: TOML file ", a, " not found")') FILE_NAME
      stop
    end if

    open (action='read', file=FILE_NAME, iostat=rc, newunit=fu)

    if (rc /= 0) then
      write (stderr, '("Error: Reading TOML file ", a, " failed")') FILE_NAME
      stop
    end if

    call toml_parse(table, fu)
    close (fu)

    if (.not. allocated(table)) then
      write (stderr, '("Error: Parsing failed")')
      stop
    end if

    ! Get [runtime] section.
    call get_value(table, 'runtime', child, requested=.false.)

    if (associated(child)) then
      ! Output server address.
      call get_value(child, 'data_path', runtime%data_path, './')
      print '(2a)', 'Data path: ', runtime%data_path
    end if

    ! Get [tree] section.
    call get_value(table, 'treee', child, requested=.false.)

    if (associated(child)) then
      ! Parameters used to modify the merger rate used in split.F90
      ! to be slightly different to the standard Press-Schechter formula
      !
      ! The modify factor is
      ! G0 [sigma(m1)/sigma(m2)]^gamma_1 [w/sigma(m2)]^gamma_2

      call get_value(child, 'G0',      pa_tree%G0,       0.57)
      call get_value(child, 'gamma_1', pa_tree%gamma_1,  0.38)
      call get_value(child, 'gamma_2', pa_tree%gamma_2, -0.01)
      call get_value(child, 'eps1',    pa_tree%eps1,     0.1)
      call get_value(child, 'eps2',    pa_tree%eps2,     0.1)

    end if
 end subroutine parse_parameter_file

end module Parameter_File
