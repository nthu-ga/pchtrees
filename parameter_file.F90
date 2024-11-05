! example.f90
module Parameter_File

  public :: parse_parameter_file
  public :: runtime

  type Parameters_Runtime
    character(len=:), allocatable :: data_path
  end type Parameters_Runtime

  type(Parameters_Runtime) :: runtime

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
end subroutine parse_parameter_file

end module Parameter_File
