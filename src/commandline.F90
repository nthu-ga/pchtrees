module Commandline
  implicit none
    
  character(len=64) :: arg_pf_path
  character(len=64) :: arg_ntrees
  character(len=64) :: arg_mphalo
  character(len=64) :: arg_ahalo
  character(len=64) :: arg_zmax
  character(len=64) :: arg_mmax
  character(len=64) :: arg_nlev

  logical :: found_pf_path = .false.
  logical :: found_ntrees  = .false.
  logical :: found_mphalo  = .false.
  logical :: found_ahalo   = .false.
  logical :: found_zmax    = .false.
  logical :: found_nlev    = .false.
  logical :: found_switch_defaults = .false.
  logical :: found_switch_verbose  = .false.
  logical :: found_task_no_output_trees  = .false.
  logical :: found_task_process_first_order_progenitors = .false.

  logical :: found_mmax = .false.
  logical :: found_switch_loguniform  = .false.

  logical :: is_kw_pf_path = .false.
  logical :: is_kw_ntrees  = .false.
  logical :: is_kw_mphalo  = .false.
  logical :: is_kw_ahalo   = .false.
  logical :: is_kw_zmax    = .false.
  logical :: is_kw_nlev    = .false.

contains 

  subroutine read_command_line_args()
    implicit none
    
    integer :: i, nargs
    character(len=64) :: argval 

    nargs = command_argument_count()

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
            found_switch_defaults = .true.
            i = i + 1
          end if
        case ('--verbose')
          ! Print more output
          if (.not.found_switch_verbose) then
            found_switch_verbose = .true.
            i = i + 1
          end if
        case ('--loguniform')
          ! Random sample uniformly in log10 mass
          if (.not.found_switch_loguniform) then
            found_switch_loguniform = .true.
            i = i + 1
          end if
        case ('--no-output-trees')
          ! Do not write any tree ouput 
          if (.not.found_task_no_output_trees) then
            found_task_no_output_trees = .true.
            i = i + 1
          end if
        case ('--process-first-order-progenitors')
          ! Do not write any tree ouput 
          if (.not.found_task_process_first_order_progenitors) then
            found_task_process_first_order_progenitors = .true.
            i = i + 1
          end if
        case ('--params')
          ! Parameter file path
          if (i + 1 <= nargs) then
            if (.not.found_pf_path) then
              call get_command_argument(i + 1, arg_pf_path)
              is_kw_pf_path = .true.
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
              call get_command_argument(i + 1, arg_ntrees)
              is_kw_ntrees = .true.
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
              call get_command_argument(i + 1, arg_mphalo)
              is_kw_mphalo = .true.
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
              call get_command_argument(i + 1, arg_ahalo)
              is_kw_ahalo = .true.
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
              call get_command_argument(i + 1, arg_zmax)
              is_kw_zmax = .true.
              found_zmax = .true.
              i = i + 2
            else
              write(*,*) 'Argument given as keyword and positional: zmax'
              stop
            end if
          end if
        case ('--mmax')
          ! Maximum halo masse
          if (i + 1 <= nargs) then
            if (.not.found_mmax) then
              call get_command_argument(i + 1, arg_mmax)
              ! This argument can only be a keyword
              found_mmax = .true.
              i = i + 2
            else
              write(*,*) 'Argument given as keyword and positional: mmax'
              stop
            end if
          end if
        case ('--nlev')
          ! Override nlev from parameter file
          if (i + 1 <= nargs) then
            if (.not.found_nlev) then
              call get_command_argument(i + 1, arg_nlev)
              is_kw_nlev = .true.
              found_nlev = .true.
              i = i + 2
            else
              write(*,*) 'Argument given as keyword and positional: nlev'
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
          call get_command_argument(i, arg_pf_path)
          found_pf_path = .true.
          i = i + 1
        else if (.not.found_ntrees) then
          call get_command_argument(i, arg_ntrees)
          found_ntrees = .true.
          i = i + 1
        else if (.not.found_mphalo) then
          call get_command_argument(i, arg_mphalo)
          found_mphalo = .true.
          i = i + 1
        else if (.not.found_ahalo) then
          call get_command_argument(i, arg_ahalo)
          found_ahalo = .true.
          i = i + 1
        else if (.not.found_zmax) then
          call get_command_argument(i, arg_zmax)
          found_zmax = .true.
          i = i + 1
        else
          write(*,*) "Too many positional arguments!"
          stop
        end if
      end if is_keyword_arg
    end do loop_over_args

    if (found_switch_verbose) then
      write(*,*) 'Have', nargs, 'arguments'
      if (found_switch_verbose)  call write_kw_or_pos('verbose',  .true.)
      if (found_switch_defaults) call write_kw_or_pos('defaults', .true.)  
      if (found_pf_path)         call write_kw_or_pos('path',   is_kw_pf_path)
      if (found_ntrees)          call write_kw_or_pos('ntrees', is_kw_ntrees)
      if (found_mphalo)          call write_kw_or_pos('mphalo (target root mass, Msol)',   is_kw_mphalo)
      if (found_ahalo)           call write_kw_or_pos('ahalo (root expansion factor',      is_kw_ahalo)
      if (found_zmax)            call write_kw_or_pos('zmax (highest redshift in tree)',   is_kw_zmax)
      if (found_nlev)            call write_kw_or_pos('nlev (number of tree levels)',      is_kw_zmax)
      
      if (found_mmax)            call write_kw_or_pos('mmax (upper limit of mass sampling range)', .true.)
    end if

    if (.not.found_switch_defaults) then
      ! Check for valid parameters (unless we're printing defaults) 
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
    end if
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
    write(*,*) 'ahalo  (--ahalo)  : Expansion factor at root of tree (1.0)'
    write(*,*) 'zmax   (--zmax)   : highest redshift in tree (4.0)' 
    write(*,*) 
    write(*,*) 'Options (keyword only):'
    write(*,*) '--nlev : number of levels in tree'
    write(*,*) '--mmax : upper limit of mass sampling range'
    write(*,*) '--loguniform : random uniform sampling in log10 mass'
    write(*,*) 

  end subroutine usage

  subroutine write_kw_or_pos(argname, is_kw)
    character(len=*), intent(in) :: argname
    logical, intent(in) :: is_kw

    if (is_kw) then
        write(*,*) 'Keyword argument: ', trim(argname)
      else
        write(*,*) 'Positional argument: ', trim(argname)
    endif
  end subroutine write_kw_or_pos

end module Commandline
