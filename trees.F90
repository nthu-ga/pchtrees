program tree
  use Defined_Types ! defined_types.F90
  use Cosmological_Parameters ! cosmological_parameters.F90
  use Power_Spectrum_Parameters ! parameters.F90
  use Tree_Memory_Arrays ! memory_modules.F90
  use Tree_Memory_Arrays_Passable ! memory_modules.F90
  use Time_Parameters ! parameters.F90
  use Tree_Routines ! tree_routines.F90
  use Modified_Merger_Tree ! modified_merger_tree.F90
  use Overdensity
  use Parameter_File
  use HDF5
  use IO
  use Commandline
  implicit none

  type (TreeNode), pointer :: This_Node
  integer :: i,j,ncount
  integer :: itree, ntree
  ! APC seems this is not used?
  ! integer, parameter :: long = selected_real_kind(9,99)
  real, allocatable  :: wlev(:),alev(:)
  integer, allocatable  :: ifraglev(:)
  real :: mphalo,ahalo,sigmacdm,zmax
  integer :: ierr,nhalomax,nhalo,ilev
  integer, allocatable :: nhalolev(:),jphalo(:)
  integer :: iter, iseed0, iseed
  EXTERNAL sigmacdm,split
  real :: dc

  ! Treewalk counter
  integer :: inode

  ! Trees table
  integer, allocatable :: trees_nhalos(:)

  ! HDF5 output
  integer(hsize_t) :: N_min, N_max
  integer :: hdferr

  ! Setup
  call read_command_line_args()
  call parse_parameter_file(trim(arg_pf_path))
  call h5open_f(hdferr)
  
  ! Set the cosmology parameters from the parameter file values
  ! 
  ! FUTURE: consider a neater / more robust way to do this, but 
  ! without introducing unnecessary module dependencies.
  !
  h0      = pa_cosmo%h0
  omega0  = pa_cosmo%omega0
  lambda0 = pa_cosmo%lambda0
  omegab  = pa_cosmo%omegab
  CMB_T0  = pa_cosmo%CMB_T0
  Gamma   = pa_cosmo%Gamma

  ! Mass of halo for which the tree is to be grown. The mass resolution of the
  ! tree and the number of trees to grow. 
  
  mphalo=1.0e+14  !halo mass at base of tree
  ntree=2         !number of trees

! Cosmological and Power Spectrum parameters
! (passed in module  Cosmological_Parameters and Power_Spectrum_Parameters)

 pkinfile='pk_Mill.dat' !Tabulated Millennium Simulation linear P(k)
! itrans=-1  !indicates use transfer function tabulated in file pkinfile
  itrans=1   !indicates use BBKS CDM transfer function with specified Gamma and Omega0
!  itrans=2   !indicates use Bond & Efstathiou CDM transfer function with specified Gamma and Omega0
! itrans=3   !indicates use Eisenstein and Hu CDM transfer function with specified Omega0, Omegab and h0
! CMB_T0=2.73 !For Eisenstein and Hu CDM transfer function one must specify the CMB temperature

!Set primordial P(k) parameters (ignored if itrans=-1)
 nspec=1.0     !primoridial power spectrum spectral index
 dndlnk=0.0    !allow running spectral index by setting ne.0
 kref=1.0      !pivot point for running index

  ierr = 1     !initial error status us to control make_tree()
  nhalomax = 0 !initialise
  nhalo = 0
  
  ! Set initial seed
  iseed0 = pa_runtime%iseed0

  ! Set up the array of redshifts at which the tree is to be stored
  write(0,*) 'The redshifts at which the tree will be stored:'
  allocate(wlev(pa_output%nlev))
  allocate(alev(pa_output%nlev))
  allocate(ifraglev(pa_output%nlev))

  ! Specify output/storage times of the merger tree
  ahalo=1.0       !expansion factor at base of tree
  zmax=4.0        !maximum redshift of stored tree
  do ilev=1,pa_output%nlev  !tree levels uniform between z=0 and zmax
    alev(ilev)=1.0/(1.0+zmax*real(ilev-1)/real(pa_output%nlev-1))
    dc = deltcrit(alev(ilev))
    write(0,'(a2,1x,f6.3,1x,a,f6.3)')'z=',(1/alev(ilev)) -1.0,'at which deltcrit=',dc
  end do

  ! FIXME estimate these numbers better
  N_min = 100
  N_max = 1000
  call create_hdf5_output(pa_output%file_path, N_min, N_max) 

  allocate(trees_nhalos(ntree))

  ! Allocate tree workspace
  allocate(nhalolev(pa_output%nlev))
  allocate(jphalo(pa_output%nlev))

  ! Start generating trees
  generate_trees: do itree=1,ntree
    iter = 1   

    ! If we run out of allocated memory, which is flagged
    ! by ierr=1 or ierr=2, then we do another iteration 
    ! with more allocated memory.

    build_tree: do while ((ierr.ne.0).or.(iter.eq.1))
      if (iter.eq.1) iseed0 = iseed0 - 19 ! Advance seed for new tree
      iseed = iseed0

      ! Allocate memory
      ! If needed, increase the amount of memory allocated
      call Memory(nhalo,nhalomax,ierr,pa_output%nlev,mphalo,&
        & pa_output%mres)
      do j=1, nhalomax, 1
        MergerTree_Aux(j)%index = j
      end do
      MergerTree => MergerTree_Aux  !Maps MergerTree to allocated 

      ! Build the tree
      call make_tree(mphalo,ahalo,pa_output%mres,alev,pa_output%nlev,&
        & iseed,split,sigmacdm, &
        & nhalomax,ierr,nhalo,nhalolev,jphalo,wlev,ifraglev)

      iter=iter+1
    end do build_tree

    !    You might want to insert your own code here and pass it the
    !    tree.

    ! Write the tree to output
    ! Label trees with 0-based index itree-1
    This_Node => MergerTree(1)
    call write_tree_hdf5(pa_output%file_path, itree-1, This_Node, &
      & nhalo, pa_output%nlev)
    
    write(0,*) 'Made tree',itree, nhalo, This_Node%mhalo
    
    ! Record number of halos for trees table
    trees_nhalos(itree) = nhalo
    
    !   Write out the information for the first couple of
    !   halos in the tree
    !write(0,*) 'Example information from the tree:'
    !This_Node => MergerTree(1)
    !write(0,*) 'Base node:'
    !write(0,*) '  mass=',This_node%mhalo,' z= ',1.0/alev(This_node%jlevel)-1.0,' number of progenitors ',This_node%nchild
    !This_Node => This_node%child !move to first progenitor
    !write(0,*) 'First progenitor:'
    !write(0,*) '  mass=',This_node%mhalo,' z= ',1.0/alev(This_node%jlevel)-1.0
    !This_Node => This_node%sibling !move to 2nd progenitor
    !write(0,*) '  mass=',This_node%mhalo,' z= ',1.0/alev(This_node%jlevel)-1.0

  end do generate_trees

  ! Write the trees table
  call write_tree_table(pa_output%file_path, trees_nhalos)

  ! Write the aexp list
  call write_output_times(pa_output%file_path, alev)

  ! Write the header
  ! Currently only support single file output
  call write_header(pa_output%file_path, '/Header', pa_output%nlev-1, & 
    & sum(trees_nhalos), sum(trees_nhalos), ntree, ntree, 1)
 
  ! Write parameters
  call write_parameters(pa_output%file_path)

  deallocate(trees_nhalos)
  deallocate(wlev,alev,ifraglev)
  deallocate(nhalolev)
  deallocate(jphalo)

  ! Tidy up HDF5
  call h5close_f(hdferr)

end program tree
