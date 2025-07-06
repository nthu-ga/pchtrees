module file_utils
  implicit none

contains

  subroutine delete_file(filename, deleted, ierr)
    ! Delete a file if it exists
    character(len=*), intent(in)   :: filename
    logical, intent(out)           :: deleted
    integer, intent(out), optional :: ierr

    integer, parameter :: unit = 99
    integer :: ios
    logical :: exists

    ! Check if the file exists
    inquire(file=filename, exist=exists)
    if (.not.exists) then
      ! The file didn't exist
      deleted = .false.
      if (present(ierr)) ierr = 0
      return
    end if

    call execute_command_line("test -f " // filename, exitstat=ierr)
    if (ierr.ne.0) then
      write(0,*) 'Trying to delete a directory: ', filename
      deleted = .false.
      if (present(ierr)) ierr = 0
      return
    end if

    ! Try opening the file
    open(unit=unit, file=filename, status='OLD', iostat=ios)
    if (ios == 0) then
      ! Close unit (deletes the file if open succeeded)
      close(unit, status='DELETE', iostat=ios)
      if (ios == 0) then
        deleted = .true.
      else
        deleted = .false.
      end if
    else
      ! The file could not be opened
      deleted = .false. 
    endif

    if (present(ierr)) ierr = ios
  end subroutine delete_file

end module file_utils
