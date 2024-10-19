program tree

 use Defined_Types ! defined_types.F90
 use Cosmological_Parameters ! parameter_modules.F90
 use Power_Spectrum_Parameters ! parameters.F90
 use Tree_Memory_Arrays ! memory_modules.F90
 use Tree_Memory_Arrays_Passable ! memory_modules.F90
 use Time_Parameters ! parameters.F90
 use Tree_Routines ! tree_routines.F90
 use Modified_Merger_Tree ! modified_merger_tree.F90

implicit none

  type (TreeNode), pointer :: This_Node
  integer :: i,j,count,ntree
  integer, parameter :: long = selected_real_kind(9,99)
  real, allocatable  :: wlev(:),alev(:)
  integer, allocatable  :: ifraglev(:)
  real :: mphalo,mres,ahalo,deltcrit,sigmacdm,zmax
  integer, parameter:: nlev=10 !number of levels of the tree to store
  integer :: ierr,nhalomax,nhalo,nhalolev(nlev),jphalo(nlev),ilev
  integer :: iter,iseed0,iseed
  EXTERNAL deltcrit,sigmacdm,split
  real :: dc
!
!Mass of halo for which the tree is to be grown. The mass resolution of
!the tree and the number of trees to grow. 
 mphalo=1.0e+14  !halo mass at base of tree
 mres = 1.0e+08  !mass resolution
 ntree=2         !number of trees

! Parameters of the Merger Tree Algorithm as defined in 
! Parkinson, Cole and Helly (2007arXiv0708.138 version 3 and in MNRAS paper)
! These values supercede the values given in the
! original astro-ph posting due to a small error in the code being
! identified. In this version of the code the error has been rectified
! and the fits redone. Using this code and these new parameters will
! produce near identical results to the old code with the old parameters.
! (passed in module Modified_Merger_Tree and Time_Parameters)
 G0=0.57
 gamma_1=0.38
 gamma_2=-0.01
 eps1=0.1        
 eps2=0.1        

!
! Cosmological and Power Spectrum parameters
! (passed in module  Cosmological_Parameters and Power_Spectrum_Parameters)

 pkinfile='pk_Mill.dat' !Tabulated Millennium Simulation linear P(k)
! itrans=-1  !indicates use transfer function tabulated in file pkinfile
  itrans=1   !indicates use BBKS CDM transfer function with specified Gamma and Omega0
!  itrans=2   !indicates use Bond & Efstathiou CDM transfer function with specified Gamma and Omega0
! itrans=3   !indicates use Eisenstein and Hu CDM transfer function with specified Omega0, Omegab and h0
! CMB_T0=2.73 !For Eisenstein and Hu CDM transfer function one must specify the CMB temperature
 omega0=0.25 
 lambda0=0.75
 h0=0.73
 omegab=0.04
 Gamma=omega0*h0  ! Omega_m.h  ignoring effect of baryons


!Set primordial P(k) parameters (ignored if itrans=-1)
 nspec=1.0     !primoridial power spectrum spectral index
 dndlnk=0.0    !allow running spectral index by setting ne.0
 kref=1.0      !pivot point for running index


 sigma8=0.9   !power spectrum amplitude set regardless of other parameters



!
!
  ierr=1     !initial error status us to control make_tree()
  nhalomax=0 !initialise
  nhalo=0
  iseed0=-8635 !random number seed
  iseed=iseed0

! Set up the array of redshifts at which the tree is to be stored
  write(0,*) 'The redshifts at which the tree will be stored:'
  allocate(wlev(nlev),alev(nlev),ifraglev(nlev))
! Specify output/storage times of the merger tree
  ahalo=1.0       !expansion factor at base of tree
  zmax=4.0        !maximum redshift of stored tree
  do ilev=1,nlev  !tree levels uniform between z=0 and zmax
     alev(ilev)=1.0/(1.0+zmax*real(ilev-1)/real(nlev-1))
     dc = deltcrit(alev(ilev))
     write(0,'(a2,1x,f6.3,1x,a,f6.3)')'z=',(1/alev(ilev)) -1.0,'at which deltcrit=',dc
  end do


!Start generating trees
do i=1,ntree

 iter=1   !if we run out of allocated memory, which is flagged
     !by ierr=1 or ierr=2 then we do another iteration 
     !with more allocated memory
    do while (ierr.ne.0 .or. iter.eq.1) 
     if (iter.eq.1) iseed0=iseed0-19  !advance seed for new tree
     iseed=iseed0
     !if needed increase the amount of memory allocated
        call Memory(nhalo,nhalomax,ierr,nlev,mphalo,mres)
        do j = 1, nhalomax, 1
           MergerTree_Aux(j)%index = j
        end do
        MergerTree => MergerTree_Aux  !Maps MergerTree to allocated 

     call make_tree(mphalo,ahalo,mres,alev,nlev,iseed,split,sigmacdm,deltcrit,&
             & nhalomax,ierr,nhalo,nhalolev,jphalo,wlev,ifraglev)
     iter=iter+1
     end do

     write(0,*)'made a tree',i
write(0,*) 'Omega_0=',omega0,'Lambda_0=',lambda0
write(0,*) 'sigma_8=',sigma8,'Gamma=',Gamma

     write(0,*)'Counting the number of nodes in the tree.'
!
!    You might want to insert your own code here and pass it the
!    tree.

     This_Node => MergerTree(1)
     count = 0
     do while (associated(This_Node))
        count = count + 1
        This_Node => Walk_Tree(This_Node)
     end do
     write(0,'(a,i3,a,i8)') 'number of nodes in tree',i,' is',count

!   Write out the information for the first couple of
!   halos in the tree
     write(0,*) 'Example information from the tree:'
     This_Node => MergerTree(1)
     write(0,*) 'Base node:'
     write(0,*) '  mass=',This_node%mhalo,' z= ',1.0/alev(This_node%jlevel)-1.0,' number of progenitors ',This_node%nchild
     This_Node => This_node%child !move to first progenitor
     write(0,*) 'First progenitor:'
     write(0,*) '  mass=',This_node%mhalo,' z= ',1.0/alev(This_node%jlevel)-1.0
     This_Node => This_node%sibling !move to 2nd progenitor
     write(0,*) '  mass=',This_node%mhalo,' z= ',1.0/alev(This_node%jlevel)-1.0


end do

deallocate(wlev,alev,ifraglev)

end program tree
