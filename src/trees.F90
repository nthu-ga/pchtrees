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
  use HDF5
  use IO
  use Commandline
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
  real :: mphalo,ahalo,sigmacdm,zmax
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

  ! Start a new output file every N trees
  integer :: ifile
  integer :: nfiles
  integer :: first_tree_in_file
  integer, allocatable :: ntrees_per_file(:)
  integer, allocatable :: nhalos_per_file(:)
  character(len=1024)  :: file_path
  
  ! HDF5 output
  integer(hsize_t) :: N_min, N_max
  integer :: hdferr

  ! Set variables from command line
  call read_command_line_args()

  read(arg_ntrees, '(I10)', IOSTAT=ierr) ntrees ! Number of trees
  read(arg_mphalo, *, IOSTAT=ierr) mphalo ! Halo mass at base of tree
  read(arg_ahalo, *, IOSTAT=ierr) ahalo ! Root expansion factor
  read(arg_zmax, *, IOSTAT=ierr) zmax ! Higest redshift in tree

  ! Can override parmeterfile nlev on the commandline
  if (found_nlev) then
    read(arg_nlev, *, IOSTAT=ierr) nlev ! Higest redshift in tree
  else
    nlev = pa_output%nlev
  endif 

  write(*,'(1x,a,i10)')   'Trees to generate                : ', ntrees
  write(*,'(1x,a,g10.3)') 'Target halo mass (Msol)          : ', mphalo
  write(*,'(1x,a,f10.3)') 'Expansion factor at root of tree : ', ahalo
  write(*,'(1x,a,f10.3)') 'Maximum redshift for tree nodes  : ', zmax
  write(*,'(1x,a,i10)')   'Tree levels                      : ', nlev

  ! Process the parameter file
  if (found_switch_defaults) then
    ! Special case, just print default parameters in TOML format
    call parse_parameter_file(trim(arg_pf_path), dump_parameters_unit=6)
    stop
  else
    call parse_parameter_file(trim(arg_pf_path))
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

  ! Set up the array of redshifts at which the tree is to be stored
  write(*,*) 
  write(*,*) 'The redshifts at which the tree will be stored:'
  allocate(wlev(nlev))
  allocate(alev(nlev))
  allocate(ifraglev(nlev))

  ! Specify output/storage times of the merger tree
  do ilev=1,nlev  !tree levels uniform between z=0 and zmax
    alev(ilev) = 1.0/(1.0+zmax*real(ilev-1)/real(nlev-1))
    if (found_switch_verbose) then
      write(0,'(a2,1x,f6.3,1x,a,f6.3)')'z=',(1/alev(ilev)) -1.0, &
        & 'at which deltcrit=', deltcrit(alev(ilev))
    end if
  end do

  ! Set up HDF5 output file
    call h5open_f(hdferr)

  ! Allocate tree workspace
  allocate(trees_nhalos(ntrees), source=0)
  allocate(nhalolev(nlev),       source=0)
  allocate(jphalo(nlev),         source=0)

  ierr     = 1 ! initial error status us to control make_tree()
  nhalomax = 0 ! initialise
  nhalo    = 0
  
  nfiles = ceiling(real(ntrees) / real(pa_runtime%max_trees_per_file))
  allocate(ntrees_per_file(nfiles), source=0)
  allocate(nhalos_per_file(nfiles), source=0)
  
  ! Set up first output file
  ifile = 1
  first_tree_in_file = 1
  write(file_path, '(A, A, I3.3, A)') trim(pa_output%file_base),'.', ifile, ".hdf5"
  write(*,*) trim(file_path)

  ! FIXME estimate these numbers better
  N_min = 100
  N_max = 1000
  call create_hdf5_output(file_path, N_min, N_max) 

  ! Start generating trees
  generate_trees: do itree=1,ntrees
    iter = 1   

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
      call make_tree(mphalo,ahalo,pa_output%mres,alev,nlev,&
        & iseed,split,sigmacdm, &
        & nhalomax,ierr,nhalo,nhalolev,jphalo,wlev,ifraglev)

      iter=iter+1
    end do build_tree

    ! Write the tree to active output file
    ! Label trees with 0-based index itree-1
    This_Node => MergerTree(1)
    call write_tree_hdf5(file_path, itree-1, This_Node, &
      & nhalo, nlev)

    write(0,*) 'Wrote tree',itree, nhalo, This_Node%mhalo
    ntrees_per_file(ifile) = ntrees_per_file(ifile) + 1

    ! Record number of halos for trees table
    trees_nhalos(itree) = nhalo
    nhalos_per_file(ifile) = nhalos_per_file(ifile) + nhalo

    ! Close and refresh output file if needed
    if ((ntrees_per_file(ifile).eq.pa_runtime%max_trees_per_file).or.(itree.eq.ntrees)) then
      write(*,*) "Writing trees", first_tree_in_file, itree

      ! Write the aexp list
      call write_output_times(file_path, alev)
     
      ! Write the trees table
      call write_tree_table(file_path, trees_nhalos(first_tree_in_file:itree))

      ! Write parameters
      call write_parameters(file_path)
      
      ! Set up the next file
      if (itree.lt.ntrees) then
        ifile = ifile + 1
        first_tree_in_file = itree + 1
        write(file_path, '(A, A, I3.3, A)') trim(pa_output%file_base),'.', ifile, ".hdf5"
        write(*,*) trim(file_path)
        call create_hdf5_output(file_path, N_min, N_max) 
      end if
    end if 

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

  ! The number of files we have written
  if (.not.(ifile.eq.nfiles)) then
    write(*,*) 'FATAL: mismatch in expected number of files written'
    stop
  end if

  ! Write the trees table
  ! all write_tree_table(pa_output%file_path, trees_nhalos)
  
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

  deallocate(nhalos_per_file)
  deallocate(ntrees_per_file)
  deallocate(trees_nhalos)
  deallocate(wlev,alev,ifraglev)
  deallocate(nhalolev)
  deallocate(jphalo)

  ! Tidy up HDF5
  call h5close_f(hdferr)

end program tree
