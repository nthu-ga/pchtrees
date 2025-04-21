program tree
  use Defined_Types ! defined_types.F90
  use Cosmological_Parameters ! cosmological_parameters.F90
  use Power_Spectrum_Parameters ! parameters.F90
  use Runtime_Parameters
  use Time_Parameters ! parameters.F90
  use Tree_Memory_Arrays ! memory_modules.F90
  use Tree_Memory_Arrays_Passable ! memory_modules.F90
  use Tree_Routines ! tree_routines.F90
  use Modified_Merger_Tree ! modified_merger_tree.F90
  use Overdensity
  use Parameter_File
  use IO
  use Commandline
#ifdef WITH_HDF5
  use HDF5
#endif
  implicit none

  type (TreeNode), pointer :: This_Node
 
  integer :: j
  integer :: itree, ntrees
  integer :: ilev, nlev
  ! APC seems this is not used?
  ! integer, parameter :: long = selected_real_kind(9,99)
  real,    allocatable :: wlev(:)
  real,    allocatable :: alev(:)
  integer, allocatable :: ifraglev(:)
  real :: mphalo,ahalo,sigmacdm,zmin,zmax
  integer :: ierr,nhalomax,nhalo
  integer :: iter, iseed
  EXTERNAL sigmacdm,split

  integer, allocatable :: nhalolev(:)
  integer, allocatable :: jphalo(:)
  
  ! Tree counting example
  ! integer :: i, ncount
  ! integer :: inode

  ! Trees table
  integer, allocatable :: trees_nhalos(:)

  ! File-reading temps
  integer :: unit_num
  real    :: temp

  ! Random sampling
  real :: mphalo_max, mphalo_min = 0
  real :: rand
  integer, allocatable :: rand_seed(:)
  integer :: seed_size
 
  ! Start a new output file every N trees
  integer :: ifile
  integer :: nfiles
  integer :: first_tree_in_file
  integer, allocatable :: ntrees_per_file(:)
  integer, allocatable :: nhalos_per_file(:)
  character(len=1024)  :: file_path
  character(len=5)     :: file_ext

#ifdef WITH_HDF5
  ! HDF5 output
  integer(hsize_t) :: N_min, N_max
  integer(hsize_t) :: N_min_pfop, N_max_pfop
  integer :: hdferr
#else
  integer :: N_min, N_max
#endif

  logical :: TASK_OUTPUT_TREES = .true.
 
  ! Special task
  logical :: TASK_PROCESS_FIRST_ORDER_PROGENITORS = .false.
  character(len=1024)  :: file_path_pfop
  integer :: nfop
  integer, allocatable :: trees_nfop(:)
  real, allocatable    :: trees_mroot(:)
  real, parameter      :: FOP_MASS_LIMIT = 1e10

  file_path_pfop = './pfop.hdf5'

  ! Parse the command line
  call read_command_line_args()

  ! APC: very cautious test for an empty argument string
  ! (obvious test for '' failed on mac)
  if (arg_pf_path(1:1).eq.char(0)) then
    arg_pf_path = ''
  else
    write(*,*) 'Parameters read from: ', trim(arg_pf_path)
  end if

  ! Process the parameter file
  if (found_switch_defaults) then
    ! Special case, just print default parameters in TOML format
    call parse_parameter_file(trim(arg_pf_path), dump_parameters_unit=6)
    stop
  else
    call parse_parameter_file(trim(arg_pf_path))
  end if

  ! Set variables from command line
  read(arg_ntrees, '(I10)', IOSTAT=ierr) ntrees ! Number of trees
  read(arg_mphalo, *, IOSTAT=ierr) mphalo ! Halo mass at base of tree

  if (found_task_no_output_trees) then
    write(*,*) 'Task: DO NOT OUTPUT TREES'
    TASK_OUTPUT_TREES = .false.
  else
    write(*,*) 'Task: OUTPUT TREES'
  endif

  if (found_task_process_first_order_progenitors) then
    write(*,*) 'Task: PROCESS FIRST ORDER TREES'
    TASK_PROCESS_FIRST_ORDER_PROGENITORS = .true.
  endif

  if (found_mmax) then
    read(arg_mmax, *, IOSTAT=ierr) mphalo_max ! Halo mass at base of tree

    if (mphalo_max.lt.mphalo) then
      write(*,*) 'Invalid upper limit for mass sampling: ', mphalo_max
      stop
    end if

    ! APC: in this case we random-uniformly sample the halo masses over the range 
    ! mphalo < M < mmax rather than using the same mass for each tree.
    mphalo_min = mphalo
    
    if (found_switch_loguniform) then
      mphalo_min = log10(mphalo)
      mphalo_max = log10(mphalo_max)
    endif
  end if

  ! Set invalid defaults
  nlev  = -1
  ahalo = -1
  zmax  = -1

  ! If we have an expansion factor list, we don't care about these
  ! command line values (although they still get set to the defaults)
  if (.not.pa_output%have_aexp_list) then
    ! Note that arg_ahalo and arg_zmax are assigned default values already, so
    ! no need to check if these arguments were actually specified here.
    read(arg_ahalo, *, IOSTAT=ierr) ahalo ! Root expansion factor
    read(arg_zmax, *, IOSTAT=ierr) zmax ! Higest redshift in tree
  end if

  ! Can override parmeterfile nlev on the commandline
  if (found_nlev) then
    read(arg_nlev, *, IOSTAT=ierr) nlev ! Higest redshift in tree
  else
    nlev = pa_output%nlev
  endif 

  ! Set up output file type options
  if (TASK_OUTPUT_TREES) then
    if (pa_output%output_format.eq.OUTPUT_HDF5) then
      file_ext = ".hdf5"
    elseif (pa_output%output_format.eq.OUTPUT_JET) then
      file_ext = '.txt'
    endif  
  endif

  write(*,'(1x,a,i10)')   'Trees to generate                : ', ntrees
  write(*,'(1x,a,g10.3)') 'Target halo mass (Msol)          : ', mphalo

   if (pa_output%have_aexp_list.and.((nlev.gt.0).or.found_ahalo.or.found_zmax)) then 
     write(*,*)
     write(*,*) 'Warning: parameter file specifies an expansion factor list but'
     write(*,*) '         nlev, ahalo and/or zmax were passed on the command'
     write(*,*) '         line or in the parameter file. The expansion factor'
     write(*,*) '         list has priority.'
     write(*,*)
  else
    write(*,'(1x,a,f10.3)') 'Expansion factor at root of tree : ', ahalo
    write(*,'(1x,a,f10.3)') 'Maximum redshift for tree nodes  : ', zmax
    write(*,'(1x,a,i10)')   'Tree levels                      : ', nlev
  end if
  
  ! Set the variables in the cosmological parameters module from the parameter
  ! file values.

  ! TODO these module variables should be better named with a prefix or
  ! something.

  ! FUTURE: consider a neater / more robust way to do this initialization, but
  ! without introducing unnecessary module dependencies.
  
  omega0  = pa_cosmo%omega0
  lambda0 = pa_cosmo%lambda0
  h0      = pa_cosmo%h0
  omegab  = pa_cosmo%omegab
  CMB_T0  = pa_cosmo%CMB_T0

  ! Set variables in the power spectrum parameters module.
  !
  ! FUTURE: see above for cosmological parameters.

  itrans = pa_powerspec%itrans
  nspec  = pa_powerspec%nspec
  dndlnk = pa_powerspec%dndlnk
  kref   = pa_powerspec%kref
  gamma  = pa_powerspec%gamma
  sigma8 = pa_powerspec%sigma8

  ! Set initial random seed
  iseed0 = pa_runtime%iseed

  ! APC: I think this is OK, since the tree generator is using its own RNG.
  call random_seed(size=seed_size)
  allocate(rand_seed(seed_size))
  rand_seed(:) = iseed0
  call random_seed(put=rand_seed)

  ! Set up the array of redshifts at which the tree is to be stored
  if (pa_output%have_aexp_list) then
    write(*,*)
    write(*,*) "Using expansion factor list: ", trim(pa_output%aexp_list)

    ! First pass: Count the number of lines
    nlev = 0
    open(newunit=unit_num, file=pa_output%aexp_list, status="old", action="read", iostat=ierr)
    if (ierr /= 0) then
        write(*,*) "Error opening expansion factor list!"
        stop
    end if

    do
        read(unit_num, *, iostat=ierr) temp
        if (ierr /= 0) exit  ! Exit loop on error (end of file)
        nlev = nlev + 1
    end do
    close(unit_num)
  end if

  ! Allocate timestep arrays
  if (nlev > 1) then
    allocate(wlev(nlev))
    allocate(alev(nlev))
    allocate(ifraglev(nlev))
  else
    write(*,*) 'Invalid nlev = ', nlev
    write(*,*) 'Check the expansion factor list!'
    stop
  endif 

  if (pa_output%have_aexp_list) then
    ! Second pass: read expansion factors
    open(newunit=unit_num, file=pa_output%aexp_list, status="old", action="read", iostat=ierr)
    do ilev=1,nlev
        read(unit_num, *) alev(ilev)
    end do
    close(unit_num)

    ! Ensure descending order
    call insertion_sort_desc(alev,nlev)
  else
    ! Tree levels uniform between z=0 and zmax
    ! alev(1) is the latest time (largest aexp)

    ! APC sure there must be a better way to do this...
    zmin = (1.0/ahalo) -1 
    do ilev=1,nlev    
      alev(ilev) = 1.0/(1.0+zmin+(zmax-zmin)*real(ilev-1)/real(nlev-1))
    end do
  end if
 
  ! Print the output times to console
  if (found_switch_verbose) then
    write(*,*) 
    write(*,*) 'The redshifts at which the tree will be stored:'
    do ilev=1,nlev
        write(*,'(1x,a2,1x,f6.3,1x,a,f6.3)') 'z=',(1/alev(ilev)) -1.0, &
          & 'at which deltcrit=', deltcrit(alev(ilev))
    end do
  end if

  ! Define root level
  !
  ! APC: The PCH code was a little uncelar on this point. In principle one could
  ! have output levels `alev` that did not include the tree root `ahalo`, but
  ! the `make_trees` routine enforces alev(1) = ahalo (or has uncertain
  ! behaviour if this is not the case). 

  ! In fact, the original code doesn't use the a0 (=ahalo) argument to
  ! make_trees at all.

  ! Here we force this to be the case:
  ahalo = alev(1)

#ifdef WITH_HDF5
  ! Set up HDF5 output file
  call h5open_f(hdferr)
#endif

  ! Allocate tree workspace
  allocate(trees_nhalos(ntrees), source=0)
  allocate(nhalolev(nlev),       source=0)
  allocate(jphalo(nlev),         source=0)

  ierr     = 1 ! initial error status us to control make_tree()
  nhalomax = 0 ! initialise
  nhalo    = 0
  
  nfiles = ceiling(real(ntrees) / real(pa_output%max_trees_per_file))
  allocate(ntrees_per_file(nfiles), source=0)
  allocate(nhalos_per_file(nfiles), source=0)
  
  ! Set up first output file
  ifile = 1
  first_tree_in_file = 1
  write(file_path, '(A, A, I3.3, A, A)') trim(pa_output%file_base), '.', ifile, '.', trim(pa_output%file_ext)
  write(*,*) 'Output file:   ', trim(file_path)
  write(*,*) 'Output format: ', pa_output%output_format

  ! APC FIXME estimate these numbers better
  N_min = 100
  N_max = 1000

  if (TASK_OUTPUT_TREES) then
    if (pa_output%output_format.eq.OUTPUT_HDF5) then
#ifdef WITH_HDF5
      call create_hdf5_output(file_path, N_min, N_max) 
#else
      if (found_switch_verbose) then
        write(*,*) 'Not creating HDF5 file, HDF5 output not supported in this build'
      end if
#endif
    else if (pa_output%output_format.eq.OUTPUT_JET) then
      call create_jet_output(file_path)
    else
      write(*,*) 'FATAL'
      stop
    end if
  endif

  if (TASK_PROCESS_FIRST_ORDER_PROGENITORS) then
    n_min_pfop = 100
    n_max_pfop = 1e7
    call create_hdf5_output_process_first_order_progenitors(file_path_pfop,n_min_pfop,n_max_pfop)
    allocate(trees_nfop(ntrees), source=0)
    allocate(trees_mroot(ntrees), source=0.0)
  endif

  ! Start generating trees
  generate_trees: do itree=1,ntrees
    iter = 1   

    ! Choose the mass for this tree

    ! APC: we're using the global RNG here, I think the tree generator uses
    ! its own RNG...
    if (found_mmax) then
      call random_number(rand)  ! Generates rand in [0,1)
      if (found_switch_loguniform) then
        mphalo  = 10**(mphalo_min + (mphalo_max - mphalo_min) * rand)
      else
        mphalo  = mphalo_min + (mphalo_max - mphalo_min) * rand
      end if
    endif

    ! If we run out of allocated memory, which is flagged
    ! by ierr=1 or ierr=2, then we do another iteration 
    ! with more allocated memory.

    build_tree: do while ((ierr.ne.0).or.(iter.eq.1))
      if (iter.eq.1) iseed0 = iseed0 - 19 ! Advance seed for new tree
      iseed = iseed0

      ! Allocate memory
      ! If needed, increase the amount of memory allocated
      call Memory(nhalo,nhalomax,ierr,nlev,mphalo,&
        & pa_output%mres)
      do j=1, nhalomax, 1
        MergerTree_Aux(j)%index = j
      end do
      MergerTree => MergerTree_Aux  ! Maps MergerTree to allocated 

      ! Build the tree
      write(*,*) 'Making tree:', mphalo
      call make_tree(mphalo,ahalo,pa_output%mres,alev,nlev,&
        & iseed,split,sigmacdm, &
        & nhalomax,ierr,nhalo,nhalolev,jphalo,wlev,ifraglev)

      iter=iter+1
    end do build_tree

    ! Write the tree to active output file
    ! Label trees with 0-based index itree-1
    This_Node => MergerTree(1)
  
    ! Special functions to derive properties from trees
    if (TASK_PROCESS_FIRST_ORDER_PROGENITORS) then
      ! This sets nfop
      call process_first_order_progenitors(file_path_pfop, itree-1, This_Node, FOP_MASS_LIMIT, nfop)
      trees_nfop(itree) = nfop
      trees_mroot(itree) = This_Node%mhalo
    endif

    ! Only proceed past here if we're writing tree output
    ! Otherwise, continue with the next iteration.
    if (.not.TASK_OUTPUT_TREES) then
      cycle
    endif
    
    if (pa_output%output_format.eq.OUTPUT_HDF5) then
#ifdef WITH_HDF5
      call write_tree_hdf5(file_path, itree-1, this_node, &
        & nhalo, nlev)
#else
      write(*,*) "Not writing tree -- hdf5 format specified, but no HDF5 support!"
#endif
    elseif (pa_output%output_format.eq.OUTPUT_JET) then
      call write_tree_jet(file_path, itree-1, this_node, &
        & nhalo, nlev, alev, nhalos_per_file(ifile))
    else
      write(*,*) 'FATAL'
      stop
    endif

    write(*,*) 'Wrote tree',itree, nhalo, This_Node%mhalo
    ntrees_per_file(ifile) = ntrees_per_file(ifile) + 1

    ! Record number of halos for trees table
    trees_nhalos(itree) = nhalo
    nhalos_per_file(ifile) = nhalos_per_file(ifile) + nhalo

    ! Close and refresh output file if needed
    write_file: if ((ntrees_per_file(ifile).eq.pa_output%max_trees_per_file).or.(itree.eq.ntrees)) then
      write(*,*) "Writing trees", first_tree_in_file, itree

#ifdef WITH_HDF5
      if (pa_output%output_format.eq.OUTPUT_HDF5) then
        ! Write the aexp list
        call write_output_times(file_path, alev)
       
        ! Write the trees table
        call write_tree_table(file_path, trees_nhalos(first_tree_in_file:itree))

        ! Write parameters
        call write_parameters(file_path)
      end if
#endif

      ! Set up the next file
      need_more_files: if (itree.lt.ntrees) then
        ifile = ifile + 1
        first_tree_in_file = itree + 1

        write(file_path, '(A, A, I3.3, A, A)') trim(pa_output%file_base), '.', ifile, '.', trim(pa_output%file_ext)
        write(*,*) 'Output file:   ', trim(file_path)
        if (pa_output%output_format.eq.OUTPUT_HDF5) then
#ifdef WITH_HDF5
          call create_hdf5_output(file_path, N_min, N_max) 
#else
          if (found_switch_verbose) then
            write(*,*) 'Not creating HDF5 file, HDF5 output not supported in this build'
          end if
#endif
        else if (pa_output%output_format.eq.OUTPUT_JET) then
          call create_jet_output(file_path)
        end if
      end if need_more_files
    end if write_file

    !   You might want to insert your own code here and pass it the
    !   tree.

    !   Write out the information for the first couple of
    !   halos in the tree
    
    ! write(0,*) 'Example information from the tree:'
    ! This_Node => MergerTree(1)
    ! write(0,*) 'Base node:'
    ! write(0,*) '  mass=',This_node%mhalo,' z= ',1.0/alev(This_node%jlevel)-1.0,' number of progenitors ',This_node%nchild
    ! This_Node => This_node%child !move to first progenitor
    ! write(0,*) 'First progenitor:'
    ! write(0,*) '  mass=',This_node%mhalo,' z= ',1.0/alev(This_node%jlevel)-1.0
    ! This_Node => This_node%sibling !move to 2nd progenitor
    ! write(0,*) '  mass=',This_node%mhalo,' z= ',1.0/alev(This_node%jlevel)-1.0

  end do generate_trees

  if (TASK_PROCESS_FIRST_ORDER_PROGENITORS) then
    ! Write the header
    call write_header(file_path_pfop, '/Header', nlev-1, & 
      & sum(trees_nfop), sum(trees_nfop), &
      & ntrees, ntrees, 1)

    ! Write the per-tree data
    call write_tree_table_process_first_order_progenitors(file_path_pfop, trees_nfop, trees_mroot)
  endif

  if (TASK_OUTPUT_TREES) then
    ! The number of files we have written
    if (.not.(ifile.eq.nfiles)) then
      write(*,*) 'FATAL: mismatch in expected number of files written'
      stop
    end if

    ! Write the trees table
    ! all write_tree_table(pa_output%file_path, trees_nhalos)
  
#ifdef WITH_HDF5
    ! Loop over each output file and write the header, which needs information
    ! about all the files.

    write(*,*) 'Writing output file headers...'
    do ifile=1,nfiles
      write(file_path, '(A, A, I3.3, A)') pa_output%file_base, '.', ifile, ".hdf5"
      write(*,*) trim(file_path)

      ! Write the header
      call write_header(file_path, '/Header', nlev-1, & 
        & nhalos_per_file(ifile), sum(trees_nhalos), &
        & ntrees_per_file(ifile), ntrees, nfiles)
    end do
#endif
  endif 

  deallocate(nhalos_per_file)
  deallocate(ntrees_per_file)
  deallocate(trees_nhalos)
  deallocate(wlev,alev,ifraglev)
  deallocate(nhalolev)
  deallocate(jphalo)

#ifdef WITH_HDF5
  ! Tidy up HDF5
  call h5close_f(hdferr)
#endif

contains

  subroutine process_first_order_progenitors(filename, tree_id, root_node, fop_mass_limit, &
      & nfop)
    implicit none
    character(len=*), intent(IN) :: filename
    type (TreeNode), pointer :: root_node, This_Node, child_node
    integer, intent(IN) :: tree_id
    real,    intent(IN) :: fop_mass_limit
    integer, intent(OUT) :: nfop

    integer :: iprog
    integer :: n_merge

    real, allocatable :: mprog(:), zprog(:), mhost(:)

    integer, parameter :: GUESS_NPROG = 1e6

    allocate(mprog(GUESS_NPROG), source=-1.0)
    allocate(mhost(GUESS_NPROG), source=-1.0)
    allocate(zprog(GUESS_NPROG), source=-1.0)

    ! This will also be used to count the progenitor nodes
    ! so start at zero
    iprog = 0

    This_Node => root_node
    main_branch: do while (associated(This_Node))
      n_merge = 0
      has_children: if (associated(This_Node%child)) then
        child_node => This_Node%child
        children: do while (associated(child_node))
          above_mass_limit: if (child_node%mhalo.ge.fop_mass_limit) then 
            iprog = iprog + 1
            if (iprog.eq.GUESS_NPROG) then
              write(*,*) 'FATAL: increase GUESS_NPROG!'
              STOP
            endif

            ! Store properties of first level progenitor
            mprog(iprog) = child_node%mhalo
            mhost(iprog) = This_Node%mhalo
            zprog(iprog) = 1.0/alev(This_node%jlevel) - 1
          endif above_mass_limit 

          ! Move to sibling
          if (associated(child_node%sibling)) then
            child_node => child_node%sibling
          else
            child_node => null()
          endif
        end do children
      endif has_children
      
      ! Move along main branch
      if (associated(This_Node%child)) then
        This_Node => This_Node%child 
      else
        This_Node => null()
      endif
    end do main_branch

    ! Output
    nfop = iprog
    call write_pfop_hdf5(filename, tree_id, mprog(1:iprog), mhost(1:iprog), zprog(1:iprog))

    ! In principle this memory could be re-used...
    deallocate(mprog)
    deallocate(mhost)
    deallocate(zprog)

  end subroutine process_first_order_progenitors

  subroutine insertion_sort_desc(arr, n)
  ! A ChatGPT descending order insertion sort...
      implicit none
      integer, intent(in) :: n
      real, intent(inout) :: arr(n)
      integer :: i, j
      real :: key

      do i = 2, n
          key = arr(i)
          j = i - 1
          do while (j > 0 .and. arr(j) < key)  ! Reverse comparison
              arr(j + 1) = arr(j)
              j = j - 1
          end do
          arr(j + 1) = key
      end do
  end subroutine insertion_sort_desc

end program tree
