module Commandline
  implicit none
    
  character(len=64) :: arg_pf_path
  character(len=64) :: arg_ntrees
  character(len=64) :: arg_mphalo
  character(len=64) :: arg_ahalo
  character(len=64) :: arg_zmax
  character(len=64) :: arg_nlev

  logical :: found_pf_path = .false.
  logical :: found_ntrees  = .false.
  logical :: found_mphalo  = .false.
  logical :: found_ahalo   = .false.
  logical :: found_zmax    = .false.
  logical :: found_nlev    = .false.
  logical :: found_switch_defaults = .false.
  logical :: found_switch_verbose  = .false.

contains 

  subroutine read_command_line_args()
    implicit none
    
    integer :: i, nargs
    character(len=64) :: argval 

    nargs = command_argument_count()
    write(*,*) 'Have', nargs, 'arguments'

    if (nargs.eq.0) call usage()
  
    i = 1
    loop_over_args: do while (i.le.nargs)
      call get_command_argument(i,argval)
      is_keyword_arg: if (argval(1:1).eq.'-') then

        ! Keyword argument 
        select case (trim(argval))
        case ('--defaults')
          ! Print default options, then stop
          if (.not.found_switch_defaults) then
            write(*,*) 'Keyword argument: defaults'
            found_switch_defaults = .true.
            i = i + 1
          end if
        case ('--verbose')
          ! Print more output
          if (.not.found_switch_verbose) then
            write(*,*) 'Keyword argument: verbose'
            found_switch_verbose = .true.
            i = i + 1
          end if
        case ('--params')
          ! Parameter file path
          if (i + 1 <= nargs) then
            if (.not.found_pf_path) then
              write(*,*) 'Keyword argument: path'
              call get_command_argument(i + 1, arg_pf_path)
              found_pf_path = .true.
              i = i + 2
            else
              write(*,*) 'Argument given as keyword and positional: path'
              stop
            end if
          end if
        case ('--ntrees')
          ! Number of trees to generate
          if (i + 1 <= nargs) then
            if (.not.found_ntrees) then
              write(*,*) 'Keyword argument: number of trees'
              call get_command_argument(i + 1, arg_ntrees)
              found_ntrees = .true.
              i = i + 2
            else
              write(*,*) 'Argument given as keyword and positional: ntrees'
              stop
            end if
          end if
        case ('--mphalo')
          ! Root mass of halos to growe
          if (i + 1 <= nargs) then
            if (.not.found_mphalo) then
              write(*,*) 'Keyword argument: target root mass (Msol)'
              call get_command_argument(i + 1, arg_mphalo)
              found_mphalo = .true.
              i = i + 2
            else
              write(*,*) 'Argument given as keyword and positional: mphalo'
              stop
            end if
          end if
        case ('--ahalo')
          ! Expansion factor at root of treee
          if (i + 1 <= nargs) then
            if (.not.found_ahalo) then
              write(*,*) 'Keyword argument: expansion factor at root of tree'
              call get_command_argument(i + 1, arg_ahalo)
              found_ahalo = .true.
              i = i + 2
            else
              write(*,*) 'Argument given as keyword and positional: ahalo'
              stop
            end if
          end if
        case ('--zmax')
          ! Highest redshift in tree
          if (i + 1 <= nargs) then
            if (.not.found_zmax) then
              write(*,*) 'Keyword argument: highest redshift in tree'
              call get_command_argument(i + 1, arg_zmax)
              found_zmax = .true.
              i = i + 2
            else
              write(*,*) 'Argument given as keyword and positional:zmax'
              stop
            end if
          end if
        case ('--nlev')
          ! Override nlev from parameter file
          if (i + 1 <= nargs) then
            if (.not.found_nlev) then
              write(*,*) 'Keyword argument: highest redshift in tree'
              call get_command_argument(i + 1, arg_nlev)
              found_nlev = .true.
              i = i + 2
            else
              write(*,*) 'Argument given as keyword and positional:nlev'
              stop
            end if
          end if
        case default
          write(*,*) "Unknown command line argument!"
          write(*,*) trim(argval)
          stop
        end select

      else

        ! Positional argument
        if (.not.found_pf_path) then
          write(*,*) 'Positional argument: path'
          call get_command_argument(i, arg_pf_path)
          found_pf_path = .true.
          i = i + 1
        else if (.not.found_ntrees) then
          write(*,*) 'Positional argument: ntrees'
          call get_command_argument(i, arg_ntrees)
          found_ntrees = .true.
          i = i + 1
        else if (.not.found_mphalo) then
          write(*,*) 'Positional argument: mphalo'
          call get_command_argument(i, arg_mphalo)
          found_mphalo = .true.
          i = i + 1
        else if (.not.found_ahalo) then
          write(*,*) 'Positional argument: ahalo'
          call get_command_argument(i, arg_ahalo)
          found_ahalo = .true.
          i = i + 1
        else if (.not.found_zmax) then
          write(*,*) 'Positional argument: zmax'
          call get_command_argument(i, arg_zmax)
          found_zmax = .true.
          i = i + 1
        else
          write(*,*) "Too many positional arguments!"
          stop
        end if
      end if is_keyword_arg
    end do loop_over_args

    ! Check for valid parameters 
    if (.not.found_ntrees) then
      write(*,*) "Missing number of trees to generate, using default (1)"
      arg_ntrees = "1"
    endif

    if (.not.found_mphalo) then
      write(*,*) "Missing target halo mass, using default (1e12 Msol)"
      arg_mphalo = "1.0e+12" 
    endif

    if (.not.found_ahalo) then
      write(*,*) "Missing expansion factor at root of trees, using default (1.0)"
      arg_ahalo = "1.0" 
    endif

    if (.not.found_zmax) then
      write(*,*) "Missing highest redshift in tree, using default (4.0)"
      arg_zmax = "4.0" 
    endif
  end subroutine read_command_line_args

  subroutine usage()
    implicit none

    write(*,*) 
    write(*,*) 'Usage: ./pchtrees parameter_file_path ntrees mphalo ahalo zmax [options]'
    write(*,*) 
    write(*,*) 'Options (positional or by keyword):' 
    write(*,*) 'path   (--path  ) : path to parameter file in TOML format'
    write(*,*) 'ntrees (--ntrees) : integer number of trees to generate (1)'
    write(*,*) 'mphalo (--mphalo) : target mass of tree root notes (1e12 Msol)'
    write(*,*) 'ahalo  (-ahalo)   : Expansion factor at root of tree (1.0)'
    write(*,*) 'zmax   (--zmax)   : highest redshift in tree (4.0)' 
    write(*,*) 
    write(*,*) 'Options (keyword only):'
    write(*,*) '--nlev : number of levels in tree'
    write(*,*) 

  end subroutine usage
end module Commandline
