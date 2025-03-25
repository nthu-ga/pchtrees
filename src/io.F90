module IO
#ifdef WITH_HDF5
  use hdf5
#endif
  use Tree_Memory_Arrays_Passable ! memory_modules.F90
  use Tree_Routines ! tree_routines.F90
  use Cosmological_Parameters
  use Power_Spectrum_Parameters
  use Runtime_Parameters
  use Parameter_File
  use Overdensity
  implicit none

#ifdef WITH_HDF5
  interface write_1d_array
      module procedure write_1d_array_real, write_1d_array_integer
  end interface write_1d_array

  character(*), parameter :: DSET_TREE_INDEX_PCH = "/TreeHalos/TreeIndexPCH"
  character(*), parameter :: DSET_TREE_INDEX = "/TreeHalos/TreeIndex"
  character(*), parameter :: DSET_TREE_BRANCH = "/TreeHalos/TreeBranch"
  character(*), parameter :: DSET_TREE_ID = "/TreeHalos/TreeID"
  character(*), parameter :: DSET_TREE_FIRST_PROGENITOR = "/TreeHalos/TreeFirstProgenitor"
  character(*), parameter :: DSET_TREE_MAIN_PROGENITOR = "/TreeHalos/TreeMainProgenitor"
  character(*), parameter :: DSET_TREE_NEXT_PROGENITOR = "/TreeHalos/TreeNextProgenitor"
  character(*), parameter :: DSET_TREE_FIRST_DESCENDANT = "/TreeHalos/TreeFirstDescendant"
  character(*), parameter :: DSET_TREE_DESCENDANT = "/TreeHalos/TreeDescendant"
  character(*), parameter :: DSET_TREE_SNAPNUM = "/TreeHalos/SnapNum"
  character(*), parameter :: DSET_TREE_GROUP_M_CRIT200 = "/TreeHalos/Group_M_Crit200"
  character(*), parameter :: DSET_TREE_SUBHALO_MASS = "/TreeHalos/SubhaloMass"
#endif

contains

  subroutine write_tree_jet(filename, tree_id, Tree_Root, nnodes, nlevels, alev, idx_offset)
    !
    ! Writes tree in "JET" text format
    !
    character(len=*), intent(IN) :: filename
    type (TreeNode), pointer, intent(IN) :: Tree_Root
    integer, intent(IN) :: tree_id, nnodes, nlevels
    real, dimension(:), intent(in) :: alev
    integer, intent(IN), optional :: idx_offset

    integer :: i
    integer :: inode, ibranch
    integer :: nodes_written_so_far
    integer :: unit_num, ierr

    ! Treewalk pointers
    type (TreeNode), pointer :: This_Node

    real :: zred, mass
    integer :: idx_node, idx_desc, idx_main_prog

    if (present(idx_offset)) then
      nodes_written_so_far = idx_offset
    else
      nodes_written_so_far = 0
    end if

    ! FIXME might want to pass an open fileid

    ! JET format
    !column 0 = redshift of the halo 
    !column 1 = halo mass 
    !column 2 = halo id (in practice the line number of the file starting from 1)
    !column 3 = id of the first descendent (at lower redshift) of the halo. -1 if no
    !descendent
    !column 4 = id of the most massive progenitor of the halo (at higher redshift),
    !-1 if the halo is a first leaf (i.e. has no progenitor)

    ! FIXME
    ! Following the spec, the tree is traversed in depth first order with 
    ! walk_tree, but *written* in order of the MergerTree array...

    !write(1,*) 1.0/alev(MergerTree(count)%jlevel)-1.0,
    !        MergerTree(count)%mhalo, MergerTree(count)%index+nodes_written_so_far,&
    !                  & temp_parent, temp_flag

    ! open file 
    ! FIXME
    open(newunit=unit_num, file=filename, status="old", action="append", iostat=ierr)

    This_Node => Tree_Root
    inode   = 1
    ibranch = 1
    do while (associated(This_Node))

      zred     = 1.0/alev(MergerTree(inode)%jlevel)
      mass     = MergerTree(inode)%mhalo
      idx_node = MergerTree(inode)%index + nodes_written_so_far

      if (associated(This_Node%parent)) then
        idx_desc = This_node%parent%index + nodes_written_so_far
      else
        idx_desc = -1
      endif

      if (associated(This_Node%child)) then
        idx_main_prog = This_Node%child%index + nodes_written_so_far
      else
        idx_main_prog = -1
      end if

      write(unit_num,*) zred, mass, idx_node, idx_desc, idx_main_prog

      inode = inode + 1
      This_Node => Walk_Tree(This_Node)
    end do

    close(unit_num)

  end subroutine write_tree_jet

#ifdef WITH_HDF5

  ! ############################################################ 
  subroutine create_hdf5_output(filename, N_min, N_max)
    !
    ! Set up the output file
    ! Output file must not exist (no automatic overwrite)
    !
    implicit none
    character(len=*), intent(in) :: filename
    integer(hsize_t), intent(in) :: N_min, N_max  ! Estimated min/max append size

    logical :: file_exists
    
    integer(hid_t) :: file_id
    integer(hsize_t) :: dataset_type
    integer :: hdferr
   
    ! Open file
    inquire(file=trim(filename), exist=file_exists)
    if (file_exists) then
      write(*,*) 'Remove existing output file: ', trim(filename)
      stop
    endif 
    call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, file_id, hdferr)
    call check_hdf5_err(hdferr,"Error creating file",trim(filename))
  
    ! Create integer arrays
    dataset_type = H5T_NATIVE_INTEGER

    ! The index in the array order assigned by PCHTrees
    call create_extensible_dataset(filename, DSET_TREE_INDEX_PCH, dataset_type, &
      & N_min, N_max, hdferr)

    ! A simple 0..N index (will differ from TreeIndexRaw if nodes are written
    ! in depth-first order
    call create_extensible_dataset(filename, DSET_TREE_INDEX, dataset_type, &
      & N_min, N_max, hdferr)
 
    ! A 0..N label for each branch of the tree (0: main branch)
    call create_extensible_dataset(filename, DSET_TREE_BRANCH, dataset_type, &
      & N_min, N_max, hdferr)
 
    ! A 0..N label for each tree
    call create_extensible_dataset(filename, DSET_TREE_ID, dataset_type, &
      & N_min, N_max, hdferr)

    ! A 0..N label for each branch of the tree (0: main branch)
    call create_extensible_dataset(filename, DSET_TREE_FIRST_PROGENITOR, dataset_type, &
      & N_min, N_max, hdferr)

    ! FIXME
    call create_extensible_dataset(filename, DSET_TREE_NEXT_PROGENITOR, dataset_type, &
      & N_min, N_max, hdferr)
 
    ! FIXME
    call create_extensible_dataset(filename, DSET_TREE_MAIN_PROGENITOR, dataset_type, &
      & N_min, N_max, hdferr)
  
    ! FIXME
    call create_extensible_dataset(filename, DSET_TREE_FIRST_DESCENDANT, dataset_type, &
      & N_min, N_max, hdferr)

    ! FIXME
    call create_extensible_dataset(filename, DSET_TREE_DESCENDANT, dataset_type, &
      & N_min, N_max, hdferr)
 
    ! FIXME
    call create_extensible_dataset(filename, DSET_TREE_SNAPNUM, dataset_type, &
      & N_min, N_max, hdferr)
 
    ! Create float arrays
    dataset_type = H5T_NATIVE_REAL

    ! FIXME
    call create_extensible_dataset(filename, DSET_TREE_GROUP_M_CRIT200, dataset_type, &
      & N_min, N_max, hdferr)
 
    ! FIXME
    call create_extensible_dataset(filename, DSET_TREE_SUBHALO_MASS, dataset_type, &
      & N_min, N_max, hdferr)
  
  end subroutine create_hdf5_output
  
  ! ############################################################
  subroutine write_output_times(filename, alev)
    implicit none
    
    character(len=*), intent(in) :: filename
    real, dimension(:), intent(in) :: alev

    real, allocatable :: output_time_property(:)
    integer :: i

    ! Write expansion factors at each output time
    call write_1d_array_real(filename, '/OutputTimes/ExpansionFactor', alev)

    allocate(output_time_property(size(alev)))

    ! Write critical density at each output time
    do i=1,size(alev)
      output_time_property = deltcrit(alev(i))
    end do
    call write_1d_array_real(filename, '/OutputTimes/DeltaCrit', output_time_property)

    ! Write redshifts
    do i=1,size(alev)
      output_time_property = (1.0/alev(i))-1.0
    end do
    call write_1d_array_real(filename, '/OutputTimes/Redshift', output_time_property)

    deallocate(output_time_property)

  end subroutine write_output_times

  ! ############################################################
  subroutine write_tree_table(filename, tree_lengths)
    implicit none
    
    character(len=*), intent(in) :: filename
    integer, dimension(:), intent(in) :: tree_lengths

    integer, allocatable :: tree_property(:)
    integer :: i

    call write_1d_array_integer(filename, '/TreeTable/Length', tree_lengths)

    allocate(tree_property(size(tree_lengths)))

    ! Create offsets
    tree_property(1) = 0
    do i=2,size(tree_lengths)
      tree_property(i) = tree_lengths(i-1)
    end do
    call write_1d_array_integer(filename, '/TreeTable/StartOffset', tree_property)

    ! Create tree ids
    do i=1,size(tree_lengths)
      tree_property(i) = i - 1 ! O-based 
    end do
    call write_1d_array_integer(filename, '/TreeTable/TreeID', tree_property)

    deallocate(tree_property)

  end subroutine write_tree_table

  ! ############################################################ 
  subroutine write_parameters(filename)
    implicit none

    character(len=*), intent(in)  :: filename      ! HDF5 file name

    integer(hid_t) :: file_id, group_id
    integer :: hdferr

    ! Open file
    call open_existing_file(filename, file_id, 'write_parameters')
    
    ! Create or open group
    call h5gcreate_f(file_id, '/Parameters', group_id, hdferr)
    call h5gopen_f(file_id, '/Parameters', group_id, hdferr)

    ! Write attributes
    call write_group_attr(group_id, 'run_iseed', iseed0)
    
    call write_group_attr(group_id, 'tree_G0', pa_tree%G0)
    call write_group_attr(group_id, 'tree_gamma_1', pa_tree%gamma_1)
    call write_group_attr(group_id, 'tree_gamma_2', pa_tree%gamma_2)
    call write_group_attr(group_id, 'tree_eps_1', pa_tree%eps1)
    call write_group_attr(group_id, 'tree_eps_2', pa_tree%eps2)
  
    ! From cosmology_parameters
    call write_group_attr(group_id, 'cosmo_h0', pa_cosmo%h0)
    call write_group_attr(group_id, 'cosmo_omega0', pa_cosmo%omega0)
    call write_group_attr(group_id, 'cosmo_lambda0', pa_cosmo%lambda0)
    call write_group_attr(group_id, 'cosmo_omegab', pa_cosmo%omegab)
    call write_group_attr(group_id, 'cosmo_CMBT0', pa_cosmo%CMB_T0)

    ! From power_spectrum_parameters
    call write_group_attr(group_id, 'pspec_itrans', itrans)
    call write_group_attr(group_id, 'pspec_gamma', gamma)
    call write_group_attr(group_id, 'pspec_kref', kref)
    call write_group_attr(group_id, 'pspec_nspec', nspec)
    call write_group_attr(group_id, 'pspec_sigma8', sigma8)
    call write_group_attr(group_id, 'pspec_infile', trim(pkinfile))
    call write_group_attr(group_id, 'pspec_splinefile', trim(splinefile))
    call write_group_attr(group_id, 'pspec_tffile', trim(tffile))

    ! Close resources
    call h5gclose_f(group_id, hdferr)
    call h5fclose_f(file_id, hdferr)

  end subroutine write_parameters

  ! ############################################################ 
  subroutine write_header(filename, group_name, &
      & last_snapshot, nhalos_thisfile, nhalos_total, & 
      & ntrees_thisfile, ntrees_total, nfiles)
    implicit none

    character(len=*), intent(in)  :: filename      ! HDF5 file name
    character(len=*), intent(in)  :: group_name    ! Name of the group

    integer, intent(in) :: last_snapshot, nfiles
    integer, intent(in) :: nhalos_total, nhalos_thisfile
    integer, intent(in) :: ntrees_total, ntrees_thisfile

    integer(hid_t) :: file_id, group_id
    integer :: hdferr

    ! Open file
    call open_existing_file(filename, file_id, 'write_header')
    
    ! Create or open group
    call h5gcreate_f(file_id, group_name, group_id, hdferr)
    call h5gopen_f(file_id, group_name, group_id, hdferr)

    ! Write attributes
    call write_group_attr(group_id, 'LastSnapShotNr', last_snapshot)
    call write_group_attr(group_id, 'Nhalos_ThisFile', nhalos_thisfile)
    call write_group_attr(group_id, 'Nhalos_Total', nhalos_total)
    call write_group_attr(group_id, 'Ntrees_ThisFile', ntrees_thisfile)
    call write_group_attr(group_id, 'Ntrees_Total', ntrees_total)
    call write_group_attr(group_id, 'NumFiles', nfiles)

    ! Close resources
    call h5gclose_f(group_id, hdferr)
    call h5fclose_f(file_id, hdferr)

  end subroutine write_header

  ! ############################################################
  subroutine write_group_attr(group_id, attr_name, attr_value)
    implicit none
    integer(hid_t), intent(in) :: group_id
    character(*), intent(in) :: attr_name
    class(*), intent(in) :: attr_value

    integer(hid_t) :: attr_id, space_id, dtype_id
    integer(hsize_t) :: attr_dims(1)
    integer :: hdferr

    ! Create dataspace for scalar attributes
    attr_dims = (/ 1 /)
    call h5screate_simple_f(1, attr_dims, space_id, hdferr)

    select type(attr_value)
    type is (integer)
      call h5acreate_f(group_id, attr_name, H5T_NATIVE_INTEGER, space_id, attr_id, hdferr)
      call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, attr_value, attr_dims, hdferr)
    type is (real)
      call h5acreate_f(group_id, attr_name, H5T_NATIVE_REAL, space_id, attr_id, hdferr)
      call h5awrite_f(attr_id, H5T_NATIVE_REAL, attr_value, attr_dims, hdferr)
    type is (character(len=*))
      call h5tcopy_f(H5T_NATIVE_CHARACTER, dtype_id, hdferr)
      call h5tset_size_f(dtype_id, int(len(attr_value), kind=hsize_t), hdferr)
      call h5acreate_f(group_id, attr_name, dtype_id, space_id, attr_id, hdferr)
      call h5awrite_f(attr_id, dtype_id, attr_value, attr_dims, hdferr)
     class default
      write(*,*) "Unknown attribute type!" 
      write(*,*) "Attribute:", trim(attr_name)
      stop
    end select

    call h5aclose_f(attr_id, hdferr)
    call h5sclose_f(space_id, hdferr)
  end subroutine write_group_attr

  ! ############################################################
  subroutine write_tree_hdf5(filename, tree_id, Tree_Root, nnodes, nlevels)
    implicit none

    ! FIXME may want to pass open fileid

    character(len=*), intent(IN) :: filename
    type (TreeNode), pointer, intent(IN) :: Tree_Root
    integer, intent(IN) :: tree_id, nnodes, nlevels

    ! HDF5 file
    integer :: hdferr

    ! Counters
    integer :: i
    integer :: inode, inode_child, ibranch

    ! Treewalk pointers
    type (TreeNode), pointer :: This_Node, Child_Node    

    ! Output arrays (generic)
    integer, allocatable :: pch_to_df_index(:)
    integer, allocatable :: tree_id_array(:)
    integer, allocatable :: tree_index(:)
    integer, allocatable :: tree_index_pch(:)
    integer, allocatable :: tree_branch(:)
    integer, allocatable :: tree_first_progenitor(:)
    integer, allocatable :: tree_next_progenitor(:)
    integer, allocatable :: tree_first_descendant(:)
    integer, allocatable :: tree_descendant(:)
    integer, allocatable :: tree_snapnum(:)
    real,    allocatable :: tree_mass(:)

    ! Nodes are written in depth first order
    ! First output treewalk
    ! - count nodes
    !This_Node => Tree_Root
    !nnodes = 0
    !do while (associated(This_Node))
    !  nnodes = nnodes + 1
    !  This_Node => Walk_Tree(This_Node)
    !end do

    ! Write TreeID
    allocate(tree_id_array(nnodes))
    tree_id_array(:) = tree_id
    call append_to_dataset(filename, DSET_TREE_ID, tree_id_array, hdferr)
    deallocate(tree_id_array)

    ! Write TreeIndex (simple 0..N-1)
    allocate(tree_index(nnodes))
    tree_index(:) = (/ (i, i = 0, nnodes - 19) /)
    call append_to_dataset(filename, DSET_TREE_INDEX, tree_index, hdferr)
    deallocate(tree_index)

    ! Create map of original index to depth first index
    allocate(pch_to_df_index(nnodes))
    inode = 1
    This_Node => Tree_Root
    do while (associated(This_Node))
      pch_to_df_index(This_Node%index) = inode
      inode = inode + 1 
      This_Node => Walk_Tree(This_Node)
    end do

    ! Scratchspace to write properties
    allocate(tree_index_pch(nnodes), source=-1)
    allocate(tree_branch(nnodes), source=-1)
    allocate(tree_first_progenitor(nnodes), source=-1)
    allocate(tree_next_progenitor(nnodes), source=-1)
    allocate(tree_first_descendant(nnodes), source=-1)
    allocate(tree_descendant(nnodes), source=-1)
    allocate(tree_snapnum(nnodes), source=-1)
    allocate(tree_mass(nnodes), source=-1.0)
    
    ! Notes on indexing
    !
    ! Node%index is set when the tree is built, and not used afterwards even
    ! though the nodes are re-ordered in the array at the end of the build step. 
    !
    ! We could probably abuse %index to hold the index in the order we actually
    ! write, but no obvious need to do this.
    !
    ! Since the tree is output in DF order, there is no need to store pointer
    ! indices for DF order.
    !
    ! TODO: might want to suppport different orderings in future.
    ! Could make the mapping index array more generic.
    ! - PCH: logical grouping by siblings
    ! - Branch order
    ! - Depth first
    ! - Same order as Gadget  
    ! - etc.
  
    ! Second output treewalk
    ! Write nodes in depth-first order (but zero-based)
    This_Node => Tree_Root
    inode   = 1
    ibranch = 1
    
    do while (associated(This_Node))
      tree_index_pch(inode) = This_Node%index - 1 ! 0-based
      
      ! Record branch
      tree_branch(inode) = ibranch - 1 ! 0-based
      
      ! Record other simple properties
      tree_snapnum(inode) = nlevels - This_Node%jlevel
      tree_mass(inode) = This_Node%mhalo

      ! Record progenitor
      ! Increment branch index after processing leaf nodes
      if (associated(This_Node%child)) then
        Child_Node => This_Node%child
        
        tree_first_progenitor(inode) = pch_to_df_index(Child_Node%index) - 1 ! 0-based 
        tree_first_descendant(inode) = pch_to_df_index(This_Node%index) - 1 ! 0-based 
      
        ! Create next progenitor links for non-leaf nodes
        do while (associated(Child_Node%sibling)) 
          inode_child = pch_to_df_index(Child_Node%index)
          if (inode_child.gt.nnodes) then
            write(*,*) Child_Node%index, inode_child, nnodes
            write(*,*) "FAIL"
            stop
          endif
          ! The next progenitor is the sibling of the current child
          tree_next_progenitor(inode_child) = pch_to_df_index(Child_Node%sibling%index) - 1 ! 0-based
          ! All children descend to the same progenitor
          tree_descendant(inode_child) = inode - 1 !  0-based

          if (associated(Child_Node%sibling)) Child_Node => Child_Node%sibling
        end do
      else
        ibranch = ibranch + 1
      endif

      inode = inode + 1
      This_Node => Walk_Tree(This_Node)
    end do

    ! The PCHTrees index
    call append_to_dataset(filename, DSET_TREE_INDEX_PCH, tree_index_pch, hdferr)
    deallocate(tree_index_pch)

    ! The branch index
    call append_to_dataset(filename, DSET_TREE_BRANCH, tree_branch, hdferr)
    deallocate(tree_branch)

    ! First progenitor
    call append_to_dataset(filename, DSET_TREE_FIRST_PROGENITOR, tree_first_progenitor, hdferr)
    ! Main progenitor (always equal to first progenitor in PS trees)
    call append_to_dataset(filename, DSET_TREE_MAIN_PROGENITOR, tree_first_progenitor, hdferr)
    deallocate(tree_first_progenitor)

    ! Next progenitor
    call append_to_dataset(filename, DSET_TREE_NEXT_PROGENITOR, tree_next_progenitor, hdferr)
    deallocate(tree_next_progenitor)

    ! First descendant
    call append_to_dataset(filename, DSET_TREE_FIRST_DESCENDANT, tree_first_descendant, hdferr)
    deallocate(tree_first_descendant)

    ! Descendant
    call append_to_dataset(filename, DSET_TREE_DESCENDANT, tree_descendant, hdferr)
    deallocate(tree_descendant)

    ! Snapnum
    call append_to_dataset(filename, DSET_TREE_SNAPNUM, tree_snapnum, hdferr)
    deallocate(tree_snapnum)

    ! Mass (these definitions are the same for EPS trees)
    call append_to_dataset(filename, DSET_TREE_SUBHALO_MASS, tree_mass, hdferr)
    call append_to_dataset(filename, DSET_TREE_GROUP_M_CRIT200, tree_mass, hdferr)
    deallocate(tree_mass)

    ! - We can guess the size for the number of trees we want and then shrink
    ! - We probably want to have the option to sample the mass function
    ! - Fix parameters
    ! - Tidy up junk

    write(*,*) 'Wrote file'
    
    deallocate(pch_to_df_index)

  end subroutine write_tree_hdf5

  ! ############################################################
  subroutine create_extensible_dataset(filename, dataset_name, dataset_type, &
      & N_min, N_max, hdferr)
    !
    ! Creates a 1-d extensible dataset
    ! File must exist
    !
    implicit none

    character(len=*), intent(in) :: filename, dataset_name
    integer(hsize_t), intent(in) :: dataset_type, N_min, N_max  ! Estimated min/max append size
    integer, intent(out) :: hdferr

    integer(hid_t) :: file_id, dset_id, dspace_id, dcpl_id, lcpl_id
    integer(hsize_t) :: dims(1), maxdims(1), chunk_dims(1), N_avg

    ! Open file for reading
    call open_existing_file(filename, file_id, 'create_extensible_dataset')
        
    ! Compute geometric mean for chunk size
    ! This casts double to int
    N_avg = INT(sqrt(real(N_min * N_max, kind=8)))

    ! Round to a power of 2 for better alignment
    chunk_dims = (/ 2**nint(log(real(N_avg, kind=8))/log(2.0)) /)

    ! Ensure chunk is at least min append size
    chunk_dims = max(chunk_dims, N_min)
  
    ! Define initial and max dimensions (1D dataset)
    dims = (/ 0 /)   ! Start empty
    maxdims = (/ H5S_UNLIMITED_F /)  ! Allow unlimited growth

    ! Create dataspace
    call h5screate_simple_f(1, dims, dspace_id, hdferr, maxdims)

    ! Create dataset creation property list and enable chunking
    call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, hdferr)
    call h5pset_chunk_f(dcpl_id, 1, chunk_dims, hdferr)

    ! Create link creation property list
    call h5pcreate_f(H5P_LINK_CREATE_F, lcpl_id, hdferr)
    ! Enable create_intermediate_group property
    call h5pset_create_inter_group_f(lcpl_id, 1, hdferr)

    ! Create dataset
    call h5dcreate_f(file_id, dataset_name, dataset_type, dspace_id, &
      & dset_id, hdferr, dcpl_id, lcpl_id, H5P_DEFAULT_F)
    call check_hdf5_err(hdferr,"Error creating dataset",info=trim(dataset_name))

    ! Close resources
    call h5pclose_f(dcpl_id, hdferr)
    call h5pclose_f(lcpl_id, hdferr)
    call h5sclose_f(dspace_id, hdferr)
    call h5dclose_f(dset_id, hdferr)
    call h5fclose_f(file_id, hdferr)
  end subroutine create_extensible_dataset

  ! ############################################################
  subroutine append_to_dataset(filename, dataset_name, new_data, hdferr)
    !
    ! Append data to an existing extensible dataset
    !
    implicit none
    character(len=*), intent(in) :: filename, dataset_name
    class(*), dimension(:), intent(in) :: new_data
    integer, intent(out) :: hdferr
    
    integer(hsize_t) :: new_size

    integer(hid_t) :: file_id, dset_id, dspace_id, memspace_id
    integer(hsize_t) :: dims(1), new_dims(1), start(1), ncount(1)

    ! Open file for reading
    call open_existing_file(filename, file_id, 'append_to_dataset')
        
    ! Open dataset
    call h5dopen_f(file_id, dataset_name, dset_id, hdferr)

    ! Get current dataset size
    call h5dget_space_f(dset_id, dspace_id, hdferr)
    call h5sget_simple_extent_dims_f(dspace_id, dims, new_dims, hdferr)

    ! Compute new dataset size
    new_size = size(new_data)
    new_dims(1) = dims(1) + new_size
    call h5dset_extent_f(dset_id, new_dims, hdferr)

    ! Select hyperslab (newly appended portion)
    start  = (/ dims(1) /)   ! Start writing at the previous last index
    ncount = (/ new_size /)  ! Number of elements to append
    call h5dget_space_f(dset_id, dspace_id, hdferr)
    call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, start, ncount, hdferr)

    ! Define memory space
    call h5screate_simple_f(1, ncount, memspace_id, hdferr)

    ! Write new data
    select type(new_data)
    type is (real)
      call h5dwrite_f(dset_id, H5T_NATIVE_REAL, new_data, ncount, hdferr, &
        & memspace_id, dspace_id)
    type is (integer)
      call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, new_data, ncount, hdferr, &
        & memspace_id, dspace_id)
    class default
      write(*,*) "Unknown dataset type!" 
      stop
    end select

    ! Close resources
    call h5sclose_f(memspace_id, hdferr)
    call h5sclose_f(dspace_id, hdferr)
    call h5dclose_f(dset_id, hdferr)
    call h5fclose_f(file_id, hdferr)
  end subroutine append_to_dataset

  ! ############################################################
  subroutine write_1d_array_integer(filename, dataset_path, data)
    implicit none
    character(len=*), intent(in) :: filename       ! File name
    character(len=*), intent(in) :: dataset_path   ! Full HDF5 dataset path
    integer, dimension(:), intent(in) :: data      ! 1D array data
   
    integer(hid_t) :: file_id, dset_id, dspace_id, plist_id
    integer(hsize_t), dimension(1) :: dims
    integer :: hdferr
    
    ! Open file for reading
    call open_existing_file(filename, file_id, 'write_1d_array')

    dims(1) = size(data)
    call h5screate_simple_f(1, dims, dspace_id, hdferr)
    call check_hdf5_err(hdferr,"Error creating dataspace")
    
    ! Create link creation property list
    call h5pcreate_f(H5P_LINK_CREATE_F, plist_id, hdferr)
    ! Enable create_intermediate_group property
    call h5pset_create_inter_group_f(plist_id, 1, hdferr)

    call h5dcreate_f(file_id, trim(dataset_path), H5T_NATIVE_INTEGER, dspace_id, dset_id, hdferr, &
      & H5P_DEFAULT_F, plist_id, H5P_DEFAULT_F)
    call check_hdf5_err(hdferr,"Error creating dataset")
    
    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, data, dims, hdferr)
    call check_hdf5_err(hdferr,"Error writing dataset")
    
    call h5dclose_f(dset_id, hdferr)
    call h5sclose_f(dspace_id, hdferr)
    call h5fclose_f(file_id, hdferr)
  
  end subroutine write_1d_array_integer

  ! ############################################################
  subroutine write_1d_array_real(filename, dataset_path, data)
    implicit none
    character(len=*), intent(in) :: filename       ! File name
    character(len=*), intent(in) :: dataset_path   ! Full HDF5 dataset path
    real, dimension(:), intent(in) :: data         ! 1D array
   
    integer(hid_t) :: file_id, dset_id, dspace_id, plist_id
    integer(hsize_t), dimension(1) :: dims
    integer :: hdferr
    
    ! Open file for reading
    call open_existing_file(filename, file_id, 'write_1d_array')

    dims(1) = size(data)
    call h5screate_simple_f(1, dims, dspace_id, hdferr)
    call check_hdf5_err(hdferr,"Error creating dataspace")
    
    ! Create link creation property list
    call h5pcreate_f(H5P_LINK_CREATE_F, plist_id, hdferr)
    ! Enable create_intermediate_group property
    call h5pset_create_inter_group_f(plist_id, 1, hdferr)

    call h5dcreate_f(file_id, trim(dataset_path), H5T_NATIVE_REAL, dspace_id, dset_id, hdferr, &
      & H5P_DEFAULT_F, plist_id, H5P_DEFAULT_F)
    call check_hdf5_err(hdferr,"Error creating dataset")
    
    call h5dwrite_f(dset_id, H5T_NATIVE_REAL, data, dims, hdferr)
    call check_hdf5_err(hdferr,"Error writing dataset")
    
    call h5dclose_f(dset_id, hdferr)
    call h5sclose_f(dspace_id, hdferr)
    call h5fclose_f(file_id, hdferr)
  
  end subroutine write_1d_array_real

  ! ############################################################  
  subroutine open_existing_file(filename, file_id, caller_name)
    implicit none
    character(len=*), intent(in) :: filename
    character(len=*), intent(in), optional :: caller_name
    integer(hid_t) :: file_id
    integer :: hdferr
    logical :: file_exists

    inquire(file=trim(filename), exist=file_exists)
    if (.not.file_exists) then
      if (present(caller_name)) then
        write(*,*) '[',caller_name,'] Output file does not exist:', filename
      else
        write(*,*) 'Output file does not exist:', filename
      endif
      stop
    endif
    call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, hdferr)
    call check_hdf5_err(hdferr,"Error opening file",info=trim(filename))

  end subroutine open_existing_file

  ! ############################################################  
  subroutine check_hdf5_err(hdferr, msg, info, fatal)
    implicit none
    ! Check and report HDF5 errors

    integer, intent(in) :: hdferr
    character(len=*), intent(in) :: msg
    character(len=*), intent(in), optional :: info
    logical, intent(in), optional :: fatal
    
    logical :: is_fatal

    ! Print this message only if there's an actual error
    if (hdferr.ne.0) then
      write(*,*) 'HDF5 error:', hdferr
      write(*,*) ' ', TRIM(msg)

      ! Print an extra info string
      if (present(info)) then
        write(*,*) ' ', info
      endif
    endif
    
    ! By default, errors are fatal
    is_fatal = .true.
    if (present(fatal)) then
      is_fatal = fatal
    else 
      is_fatal = (hdferr.ne.0)
    endif

    if (is_fatal) then
      write(*,*) 'Stopping due to HDF5 error'
      stop
    endif
  end subroutine check_hdf5_err
#endif

end module io
