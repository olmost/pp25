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
  integer :: i,j,count,main_count,ntree, temp_parent, temp_flag
  integer, parameter :: long = selected_real_kind(9,99)
  
  integer, allocatable  :: ifraglev(:)
  real :: mphalo,mphalo_exp,mres,ahalo,deltcrit,sigmacdm,zmax
  integer, parameter:: nlev=286 !number of levels of the tree to store
  integer :: lines_in_zfile
  integer :: ierr,nhalomax,nhalo,nhalolev(nlev),jphalo(nlev),ilev
  real, allocatable  :: wlev(:), alev(:)
  integer :: iter,iseed0,iseed
  EXTERNAL deltcrit,sigmacdm,split
  real :: dc
!
!Mass of halo for which the tree is to be grown. The mass resolution of
!the tree and the number of trees to grow. 
 mphalo_exp = 14   !smallest halo mass
 !mphalo=10.0**mphalo_exp  !halo mass at base of tree
 mres = 10.0**(8)  !mass resolution
 ntree=1    !number of trees
 main_count = 0 !total number of nodes up to the last full tree
 
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

! Defining the cosmological parameters from Planck 2015
 omega0=0.25 
 lambda0=0.75
 h0=0.73
 omegab=0.04
 Gamma=omega0*h0  ! Omega_m.h  ignoring effect of baryons

!Set primordial P(k) parameters (ignored if itrans=-1)
 nspec=1     !primordial power spectrum spectral index
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


 !Set up the array of redshifts at which the tree is to be stored
  !write(0,*) 'The redshifts at which the tree will be stored:'
  allocate(wlev(nlev),alev(nlev),ifraglev(nlev))
 
! open(15, file='zlist_22_to_4_step10myr.dat', status='old')
!do i=1,nlev
!    read(15,*) wlev(i)
!enddo

!Specify output/storage times of the merger tree
  ahalo=1.0/(wlev(1)+1.0)       !expansion factor at base of tree
  zmax=wlev(size(wlev))        !maximum redshift of stored tree

  wlev = (/1.002,1.007,1.012,1.017,1.022,1.027,1.032,1.037,&
    1.042,1.047,1.053,1.058,1.063,1.068,1.074,1.079,1.084,1.09,1.095,1.101,1.106,1.112,1.117,1.123,1.129,&
    1.134,1.14,1.146,1.152,1.157,1.163,1.169,1.175,1.181,1.187,1.193,1.199,1.205,1.211,1.217,1.224,1.23,&
    1.236,1.242,1.249,1.255,1.262,1.268,1.275,1.281,1.288,1.295,1.301,1.308,1.315,1.322,1.329,1.335,1.342,&
    1.35,1.357,1.364,1.371,1.378,1.385,1.393,1.4,1.408,1.415,1.423,1.43,1.438,1.446,1.453,1.461,1.469,1.477,&
    1.485,1.493,1.501,1.509,1.518,1.526,1.534,1.543,1.551,1.56,1.568,1.577,1.586,1.595,1.604,1.613,1.622,&
    1.631,1.64,1.649,1.659,1.668,1.677,1.687,1.697,1.706,1.716,1.726,1.736,1.746,1.757,1.767,1.777,1.788,&
    1.798,1.809,1.82,1.83,1.841,1.852,1.863,1.875,1.886,1.897,1.909,1.921,1.932,1.944,1.956,1.969,1.981,&
    1.993,2.006,2.018,2.031,2.044,2.057,2.07,2.083,2.097,2.11,2.124,2.138,2.152,2.166,2.18,2.195,2.209,&
    2.224,2.239,2.254,2.269,2.285,2.3,2.316,2.332,2.348,2.365,2.381,2.398,2.415,2.432,2.45,2.467,2.485,&
    2.503,2.521,2.54,2.558,2.577,2.597,2.616,2.636,2.656,2.676,2.697,2.717,2.739,2.76,2.782,2.804,2.826,&
    2.849,2.872,2.895,2.919,2.943,2.967,2.992,3.017,3.043,3.069,3.095,3.122,3.149,3.177,3.205,3.233,3.262,&
    3.292,3.322,3.353,3.384,3.415,3.448,3.481,3.514,3.548,3.583,3.618,3.654,3.691,3.728,3.767,3.806,3.846,&
    3.886,3.928,3.97,4.013,4.058,4.103,4.149,4.196,4.245,4.294,4.345,4.397,4.45,4.504,4.56,4.617,4.676,4.736,&
    4.798,4.862,4.927,4.994,5.064,5.135,5.208,5.283,5.361,5.442,5.524,5.61,5.698,5.79,5.884,5.982,6.083,6.189,&
    6.298,6.411,6.529,6.652,6.779,6.913,7.052,7.197,7.349,7.508,7.675,7.851,8.036,8.23,8.436,8.653,8.883,9.127,&
    9.387,9.665,9.961,10.28,10.622,10.991,11.392,11.827,12.303,12.826,13.403,14.045,14.763,15.574,16.498,17.562,&
    18.805,20.279,22.062 /)

  do ilev=1,nlev  !tree levels uniform between z=0 and zmax
     alev(ilev)=1.0/(1.0+wlev(ilev))
     dc = deltcrit(alev(ilev))
     write(0, *) ilev
     write(0,'(a2,1x,f6.3,1x,a,f6.3)')'z=',(1/alev(ilev)) -1.0,'at which deltcrit=',dc
  end do

! open output file
open(1, file = 'mt_parkinson_z20_z1_step20myr_1halo_14msun_8mres.dat')

!Start generating trees
do i=1,ntree
 count = 0
 mphalo = 10.0**mphalo_exp
 write(0,*) "tree is = ",i, " halo mass is = ",mphalo
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

    !write(0,*) 'trees.F90, ', 'bella2', mphalo, ahalo, mres, nhalomax, iter
     call make_tree(mphalo,ahalo,mres,alev,nlev,iseed,split,sigmacdm,deltcrit,&
             & nhalomax,ierr,nhalo,nhalolev,jphalo,wlev,ifraglev)
    !write(0,*) 'trees.F90, ', 'bella3'
     iter=iter+1
     end do

    !write(0,*)'made a tree',i
    !write(0,*) 'Omega_0=',omega0,'Lambda_0=',lambda0
    !write(0,*) 'sigma_8=',sigma8,'Gamma=',Gamma

     !write(0,*)'Counting the number of nodes in the tree.'
!
!    You might want to insert your own code here and pass it the
!    tree.

     This_Node => MergerTree(1)
     do while (associated(This_Node))
        count = count + 1
    
    !   Writing output file    
        if (MergerTree(count)%nchild == 0) then
        temp_flag = -1
        else 
        temp_flag = MergerTree(count)%child%index+main_count
        endif

        if (count == 1) then 
        temp_parent = -1
        else
        temp_parent = MergerTree(count)%parent%index+main_count
        endif

        !Writing output file
        write(1,*) 1.0/alev(MergerTree(count)%jlevel)-1.0, MergerTree(count)%mhalo, MergerTree(count)%index+main_count,&
        & temp_parent, temp_flag
        !write(0,*) 1.0/alev(MergerTree(count)%jlevel)-1.0, MergerTree(count)%mhalo, MergerTree(count)%index+main_count,&
        !& temp_parent, temp_flag
        
        This_Node => Walk_Tree(This_Node)
     end do
     write(0,'(a,i3,a,i8,a,i8)') 'number of nodes in tree',i,' is',count,"  total nodes is ",main_count


!   Write out the information for the first couple of
!   halos in the tree
!     write(0,*) 'Example information from the tree:'
!     This_Node => MergerTree(1)
!     write(0,*) 'Base node:'
!     write(0,*) '  mass=',This_node%mhalo,' z= ',1.0/alev(This_node%jlevel)-1.0,' number of progenitors !',This_node%nchild
!     This_Node => This_node%child !move to first progenitor
!     write(0,*) 'First progenitor:'
!     write(0,*) '  mass=',This_node%mhalo,' z= ',1.0/alev(This_node%jlevel)-1.0
!     This_Node => This_node%sibling !move to 2nd progenitor
!     write(0,*) '  mass=',This_node%mhalo,' z= ',1.0/alev(This_node%jlevel)-1.0


mphalo_exp = mphalo_exp + 0.05
main_count = main_count+count
end do
close(1)

deallocate(wlev,alev,ifraglev)

end program tree
