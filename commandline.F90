module Commandline
  implicit none
    
  character(len=64) :: arg_pf_path
  logical :: found_pf_path = .false.

contains 

  subroutine read_command_line_args()
    implicit none
    
    integer :: i, nargs
    character(len=64) :: argval 

    nargs = command_argument_count()
    write(*,*) 'Have', nargs, 'arguments'
    i = 1
    loop_over_args: do while (i.le.nargs)
      call get_command_argument(i,argval)
      is_keyword_arg: if (argval(1:1).eq.'-') then

        ! Keyword argument 
        select case (trim(argval))
        case ('--params')
          if (i + 1 <= nargs) then
            if (.not.found_pf_path) then
              write(*,*) 'Keyword argument: path'
              call get_command_argument(i + 1, arg_pf_path)
              found_pf_path = .true.
              i = i + 1
            else
              write(*,*) 'Argument given as keyword and positional: path'
            end if
          end if
        case default
          write(*,*) "Unknown command line argument!"
          write(*,*) trim(argval)
        end select

      else

        ! Positional argument
        if (.not.found_pf_path) then
          write(*,*) 'Positional argument: path'
          call get_command_argument(i, arg_pf_path)
          found_pf_path= .true.
          i = i + 1
        else
          write(*,*) "Too many positional arguments!"
          stop
        end if

      end if is_keyword_arg

      i = i+1
    end do loop_over_args

  end subroutine read_command_line_args

end module Commandline
