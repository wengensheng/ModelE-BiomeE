!@sum BiomeE model: vegetation dynamics module "BiomeE_VEG"
!@ derived from LM3-PPA (standalone BiomeE), with a simplifed
!@ single cohort version for NASA GISS modelE.
!@author: Ensheng Weng
!---------------------------------------------------------------------
#define BiomeE_DEBUG
#define saturatedN
!---------------------------------------------------------------------
      module BiomeE_VEG
!@sum Routines to calculate growth, reproduction, and mortality rates in an patch.
!     Weng, Updated: 10-23-2018, 11/17/2021
      use ent_const
      use ent_types
      use ent_pfts, only: BMEspdata
      use cohorts, only: cohort_construct,copy_cohort,
     &           Put_c1_into_c2,MergeableCohorts
      implicit none

      public BiomeE_params_init, BiomeE_patch_init,Cohort_birth
      public BiomeE_FullDemography, BiomeE_SingleTree
      public BiomeE_soil_BGC
      !public BiomeE_biogeography
      !public ptenv

!=========================================================================
! Weng Nov 07 2016
      ! for pre-describing initial cohorts
      type initialcohorts
         integer :: pft
         ! given size
         real*8 :: n
         real*8 :: woodC
         real*8 :: woodN
         real*8 :: labC
         real*8 :: labN
         real*8 :: height
         real*8 :: dbh
         ! initialize as seedlings
         real*8 :: C_seed
         real*8 :: N_seed
      end type

!* --------------------------------------------------
! Set PFT trait initial values for all the PFTs
      ! BiomeE PFT types:
      ! Ent No.                          BiomeE No.
      ! 2:  Evergreen broadleaf;            1
      ! 4:  Evergreen needleleaf;           2
      ! 6:  Cold deciduous broadleaf;       3
      ! 7:  Drought deciduous broadleaf;    4
      ! 8:  Cold deciduous needleleaf;      5
      ! 9:  Cold shrub;                     6
      ! 10: Arid srhub;                     7
      ! 11: C3 grass;                       8
      ! 12: C4 grass;                       9
      integer,parameter :: N_BMEPFT = 9
      integer,parameter :: pftB2E(N_BMEPFT)= ! PFT No.
     &      (/2, 4, 6, 7, 8, 9, 10, 11, 12/)

      logical,dimension(N_BMEPFT) :: Woody =
     &  (/.True.,.True.,.True.,.True.,.True.,.True.,
     &    .True.,.False.,.False./)
      !1: evergreen, 2: cold deciduous, 3: drought deciduous, 4: cold & drought
      integer*8,dimension(N_BMEPFT) :: phenotype = ! Added by Weng, 2021-05-28
     &      (/1, 1, 2, 3, 2, 2, 4, 2, 4/)
      !leaftype - 1=broadleaf, 2=needleleaf, 3=monocot (not crops)
      integer*8,dimension(N_BMEPFT) :: leaftype = ! Added by Weng, 2021-05-28
     &      (/1, 2, 1, 1, 2, 1, 1, 3, 3/)
      ! photosynthetic pathway, 1-C3, 2-C4
      integer*8,dimension(N_BMEPFT)  :: pst =
     &      (/1, 1, 1, 1, 1, 1, 1, 1, 2/)

      ! Key traits for define PFTs
      real*8,dimension(N_BMEPFT) :: LMA = ! 0.02 ! kgC m-2 leaf
     &      (/7.0d-2, 14.d-2, 2.5d-2, 3.0d-2,3.0d-2,
     &        2.5d-2, 3.0d-2, 2.5d-2, 2.5d-2/)
      real*8,dimension(N_BMEPFT) :: rhoWood =  !350.  ! wood carbon density, kgC m-3
     &   (/360.,300.,350.,250.,300.,400.,400.,90.,90./)
!     &   (/360.,360.,360.,360.,360.,360.,360.,120.,120./)
      real*8,dimension(N_BMEPFT) :: LNbase = 0.6E-3    !basal leaf Nitrogen per unit area, kg N/m2, (Rubisco)
      real*8,dimension(N_BMEPFT) :: CNstrucL = 60.0 ! leaf structural tissues, 175

      ! Photosynthesis
      real*8,dimension(N_BMEPFT) :: PARabsorb =
     &      (/0.90, 0.93, 0.90, 0.90, 0.93, 0.89, 0.89, 0.86, 0.86 /)
      real*8,dimension(N_BMEPFT)  :: Vcmax =
!     &    (/54., 30.,  40., 35., 43., 25., 17., 43., 42. /) ! original
     &     (/18., 18., 22., 20., 20., 18., 18., 20., 15. /) ! tuned
      real*8,dimension(N_BMEPFT) :: m =
     &      (/6., 6., 8., 8., 8., 8., 8., 8., 5./)
      real*8,dimension(N_BMEPFT) :: b =
     &      (/2.d-3,2.d-3,2.d-3,2.d-3,2.d-3,2.d-3,2.d-3,8.d-3,2.d-3/)

      real*8,dimension(N_BMEPFT) :: Nfixrate0 = 0.0

      real*8,dimension(N_BMEPFT) :: CrownLAImax = ! max crown leaf area index
     &   (/4.8, 4.8, 4.5, 4.5, 4.0, 3.0, 3.0, 2.5, 2.5 /)
!  PFT       1,   2,   3,   4,    5,    6,   7,   8,   9

      real*8,dimension(N_BMEPFT) :: cLAI_light = 5.0
      real*8,dimension(N_BMEPFT) :: H0_cL =
     &   (/6.0, 6.0, 6.0, 6.0, 6.0, 3.0, 3.0, 0.5, 0.5/)

      ! Sensitivity to soil water
      real*8,dimension(N_BMEPFT) :: hwilt =
     &      (/ -153.d0, -153.d0, -500.d0, -500.d0, -100.d0,
     &         -153.d0, -2030.d0, -2030.d0, -2030.d0 /)
      real*8,dimension(N_BMEPFT) :: sstar =
     &      (/0.66, 0.66, 0.66, 0.66, 0.66, 0.66, 0.66, 0.50, 0.50 /)
      real*8,dimension(N_BMEPFT) :: swilt =
     &      (/0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05 /)

!  PFT       1,   2,   3,   4,    5,    6,   7,   8,   9
      real*8,dimension(N_BMEPFT) :: rhoFR = 150.    ! fine root carbon density, kgC m-3
      real*8,dimension(N_BMEPFT) :: rootLS = 0.8d0  ! fine root lifespan, year
      real*8,dimension(N_BMEPFT) :: root_r = 2.9E-4  ! radius of fine roots, m
      real*8,dimension(N_BMEPFT) :: root_perm = 1.0d-5  ! fine root membrane permeability, kg/(m3 s)
      real*8,dimension(N_BMEPFT) :: root_zeta = 0.29  !
      real*8,dimension(N_BMEPFT) :: alphaHT =  !20.  ! scalar from DBH to height of a tree
!     PFT   1,   2,   3,    4,   5,   6,   7 ,   8,   9
     & (/ 35., 30., 30., 35., 30., 20., 20., 10., 10./)
      real*8,dimension(N_BMEPFT) :: alphaCA =  ! 120. ! scalar from DBH to crown area
!     PFT  1,    2,   3,    4,   5,   6,   7 ,   8,   9
     & (/ 120.,120.,120.,120.,120.,150.,150., 60., 60./)
      real*8,dimension(N_BMEPFT) :: alphaBM = 8150. ! &  ! scalar from DBH to biomass via Farrior et al. 2013. ! not used!
!  PFT  1,   2,   3,  4,   5,   6,   7,   8,   9
       ! (/2500.,2500.,2500.,2500.,2500.,2500.,2500.,2500.,2500./)

      real*8,dimension(N_BMEPFT) :: thetaHT = 0.5 ! &  ! exponents of DBH to height
      real*8,dimension(N_BMEPFT) :: thetaCA = 1.5 ! &  ! exponents of DBH to crown area
      real*8,dimension(N_BMEPFT) :: thetaBM = 2.5 ! &  ! exponents of DBH to biomass
      real*8,dimension(N_BMEPFT) :: f_taper = 0.7 ! taper factor, from a cylinder to a tree
      real*8,dimension(N_BMEPFT) :: K_cambium = 0.1d0
      real*8,dimension(N_BMEPFT) :: f_Acam = 1.7
      real*8,dimension(N_BMEPFT) :: Acrown0    = 0.04 ! m2, plant minimum space
      real*8,dimension(N_BMEPFT) :: tauNSC     = ! 2.0, NSC use strategy
!     PFT  1,    2,   3,    4,   5,   6,   7 ,   8,   9
     & (/ 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 3.0, 3.0/)
      real*8,dimension(N_BMEPFT) :: phiRL = 2.0 ! fine root to leaf area ratio
!     PFT  1,   2,   3,  4,   5,   6,   7,   8,   9
!     & (/ 2.69, 2.69, 2.69, 2.69, 2.69, 2.69, 2.69, 2.69, 2.69/)

      real*8,dimension(N_BMEPFT) :: phiCSA = 0.54E-4 ! ratio of sapwood to leaf area
!  PFT  1,   2,   3,  4,   5,   6,   7,   8,   9
!     & (/ 2.54e-4, 2.54e-4, 2.54e-4, 2.54e-4, 2.54e-4, 2.54e-4, 2.54e-4, 2.54e-4, 2.54e-4/)

      real*8,dimension(N_BMEPFT) :: r_seed =  ! ratio of C to seed and rhizome for grasses
!  PFT      1,   2,   3,    4,   5,   6,    7,   8,   9
     &  (/0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 0.2, 0.6, 0.6/)

      real*8,dimension(N_BMEPFT) :: seedlingC0 =  !0.1 ! seedling total carbon, 0.025
!  PFT      1,     2,    3,    4,    5,    6,    7,    8,    9
     &  (/0.1,   0.1,  0.1,  0.1,  0.1, 0.05, 0.05, 0.005, 0.005/)
      real*8,dimension(N_BMEPFT) :: mu0 =  ! 0.02 ! &    ! annual mortality rate
     &     (/.025, .025, .025, .025, .025, .02, .02, .02, .02/)

      real*8,dimension(N_BMEPFT) :: D0mu =  ! 0.8 ! &    ! Mortality curve parameter
     &     (/1.5, 1.5, 1.5, 1.5, 1.5, 0.5, 0.5, 0.0, 0.0/)

      real*8,dimension(N_BMEPFT) :: A_sd =  ! 6.0 ! &    ! Mortality multifactor for seedlings
     &     (/8.d0, 8.d0, 8.d0, 8.d0, 8.d0, 19.d0, 19.d0, 0.d0, 0.d0/)

      real*8,dimension(N_BMEPFT) :: B_sd =  ! -60.d0 ! & ! Mortality sensitivity for seedlings
     &     (/-25., -25., -25., -25., -25., -25., -25., -60., -60./)

         !real*8,dimension(N_BMEPFT)    :: CNleaf0
      real*8,dimension(N_BMEPFT) :: CNroot0 = 30.0
      real*8,dimension(N_BMEPFT) :: CNsw0   = 150.0
      real*8,dimension(N_BMEPFT) :: CNwood0 = 150.0
      real*8,dimension(N_BMEPFT) :: CNseed0 = 15.0

      real*8,dimension(N_BMEPFT) :: Tc_pheno0 =  ! Critical temperature of leaf falling
!  PFT    2,    4,     6,    7,    8,    9,    10,   11,  12
     &  (/15.0, -80.0, 15.0, 15.0, 15.0, 15.0, 15.0, 5.0, 5.0/)
      real*8,dimension(N_BMEPFT) :: gdd_par1 = ! 30.d0   !50.d0   ! -68.d0
!  PFT    2,   4,   6,   7,   8,   9,   10,  11,  12
     &  (/20., 0.0, 50., 20., 50., 30., 30., 30., 30./)
      real*8,dimension(N_BMEPFT) :: gdd_par2 = ! 800. ! 650.d0  ! 638.d0
!  PFT    2,    4,   6,    7,    8,    9,    10,   11,  12
     &  (/200., 0.0, 800., 200., 800., 800., 800., 800., 800./)
      real*8,dimension(N_BMEPFT) :: gdd_par3 = -0.02d0 ! -0.01d0

      ! Plant hydraulics (for drought phenology)
      real*8,dimension(N_BMEPFT) :: betad_ON =  ! soil Betad for leaf ON
!  PFT      2,    4,   6,    7,  8,    9,   10,  11,  12
     &  (/0.0, 0.0, 0.0, 0.2, 0.1, 0.2, 0.2, 0.2, 0.2/)
      real*8,dimension(N_BMEPFT) :: betad_OFF =  ! soil Betad for leaf OFF
!  PFT      2,    4,   6,    7,   8,   9,   10,  11,  12
     &  (/0.0, 0.0, 0.0, 0.2, 0.0, 0.1, 0.1, 0.2, 0.2/)

! Soil biogeochemical cycles
      real*8 :: N_dpst = 0.0 ! 2.4E-3 ! 2.4 gN yr-1 m-2, Nitrogen deposition

      !         fineL, structuralL, microbial, fast, slow
      real*8 :: K0SOM(5)  = (/0.8, 0.25, 6.0, 1.0, 0.5/) ! turnover rate of SOM pools (yr-1)
      real*8 :: CN0SOM(5) = (/50., 150., 10., 15., 40./) ! target CN ratios of SOM
      real*8 :: K_n    = 8.d0  ! mineral Nitrogen turnover rate
      real*8 :: f_M2SOM= 0.8d0 ! the ratio of C and N returned to SOM from microbes
      real*8 :: CUE0   = 0.4  ! default microbial CUE
      real*8 :: CUEf0  = 0.45 ! CUE for high quality litter
      real*8 :: CUEs0  = 0.15 ! CUE for low quality litter
      ! Nitrogen loss
      real*8 :: fDON = 0.0 ! 0.02  ! fraction of DON production in decomposition
      real*8 :: etaN = 0.0 ! 0.25  ! leaching rate

!     -------------for debug------------------
      integer, public :: ijdebug
!      # Bonanza Creek (77, 14): 63.92°, -145.38° -> J: 77; I: 14
!      # Manitoba (73, 33);
!      # Harvard Forest (67,44);
!      # Oak Ridge (63, 39)
!      # Brazil Tapajos (44, 51: LAI too small) (use: 43, 51, );
!      # Konza LTER (65,34)
!      # Sevilleta LTER (Lat 34.36, Long -106.88), Jm=63; Im=30
!      # US-Wkg, Walnut Gulch Kendall Grassland (31.74, -109.94), Jm=61; Im=28
!      #jm = #77 #73  #67 #63 #65 #63  #61 #43
!      #im = #14 #33  #44 #39 #32 #30  #28 #51
      integer,parameter :: N_grids = 8
      integer :: grids(N_grids) =      ! Grids ijdebug
     &      (/14077,33073,44067,39063,32065,30063,28061,51043/)

!     initial soil pools and plant cohorts -- to be read in fromd data file
      real*8 :: initial_slowC    = 0.0 ! 0.5 ! 5.0 ! kgC m-2
      real*8 :: initial_passiveC = 0.0 ! 0.2 ! 2.5 ! kgC m-2
      real*8 :: QN_mineral0      = 16.0d-3 ! 0.002 ! kgN m-2, to define ecosystem N
      ! iniital plant cohorts
      integer :: N_iniCH = 1  ! number of initial cohorts Weng 2017-09-5
      integer :: iniPFT(9) = (/4, 4, 6, 7, 8, 9, 10, 11, 12/)
      real*8,dimension(N_BMEPFT)  :: iniDen = 0.04d0 ! individuals/m2

      ! BiomeE PFT types:
      ! 2:  Evergreen broadleaf; 4:  Evergreen needleleaf;
      ! 6:  Cold deciduous broadleaf; 7:  Drought deciduous broadleaf;
      ! 8:  Cold deciduous needleleaf;
      ! 9:  Cold shrub; 10: Arid srhub;
      ! 11: C3 grass; 12: C4 grass;
      real*8  :: iniC(9)   = (/0.02,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1/) ! kgC m-2
      type(initialcohorts) :: iniCH(9) = ! 4:evergreen tree; 6:deciduous tree
     &   (/initialcohorts(4,0.2,0.5,0.0002,0.2,0.02,2.,0.01,0.01,0.001),
     &     initialcohorts(6,0.2,0.5,0.0002,0.2,0.02,2.,0.01,0.2,0.02),
     &     initialcohorts(4,0.2,0.5,0.0002,0.2,0.02,2.,0.01,0.2,0.02),
     &     initialcohorts(4,0.2,0.5,0.0002,0.2,0.02,2.,0.01,0.2,0.02),
     &     initialcohorts(4,0.2,0.5,0.0002,0.2,0.02,2.,0.01,0.2,0.02),
     &     initialcohorts(4,0.2,0.5,0.0002,0.2,0.02,2.,0.01,0.2,0.02),
     &     initialcohorts(4,0.2,0.5,0.0002,0.2,0.02,2.,0.01,0.2,0.02),
     &     initialcohorts(4,0.2,0.5,0.0002,0.2,0.02,2.,0.01,0.2,0.02),
     &     initialcohorts(4,0.2,0.5,0.0002,0.2,0.02,2.,0.01,0.2,0.02)
     & /)
      type(climate_env) :: ptenv

        !namelist /pft_traits_nml/ &
        !           LMA, CrownLAImax, &
        !           rhoWood, f_taper, K_cambium,  &
        !           rhoFR, root_r, root_zeta,root_perm, &
        !           phiRL,phiCSA, &
        !           alphaHT, thetaHT, &
        !           alphaCA, thetaCA, &
        !           alphaBM, thetaBM, &
        !            r_seed, &
        !           seedlingC0

      contains

!====================================================================
!=========== Subroutines updated 09/21/2021, Ensheng Weng ===========
!====================================================================
      subroutine BiomeE_FullDemography(pp)
!@sum Vegetation dynamics driver, default, full demography
!@ Author: Ensheng Weng, 11/17/2021 updated
        use ent_types, only: max_cohorts
        implicit none
        type(patch),   intent(inout) :: pp

        !------local var -----------
        type(cohort), pointer :: cop
        integer :: nCohorts, i

         ! count cohorts
         cop => pp%tallest
         nCohorts = 0
         do while(ASSOCIATED(cop))
            nCohorts = nCohorts + 1
            cop  => cop%shorter
         enddo
         pp%nCohorts = nCohorts

          ! Get C from labile C pool
          call BiomeE_PlantGrowth(pp)
          call BiomeE_turnover_all(pp)
          call BiomeE_veg_update(pp)
          ! Grass proportional thinning
          call BiomeE_GrassThinning(pp)
          ! update patch LAI
          call get_patch_LAI (pp)
          ! call Ent-PPA diagnostics, including daily and yearly
          call BiomeE_diagnostics(pp)

          ! ----- Call yearly subroutines at the first day after growing season
          ! use phenology flags ? Weng, 4/26/2018
          !if(ecop%gdd .eq. 0.0 .and. pp%tallest%C_seed > 0.000001) then
          if( (pp%cellptr%hemi ==1 .and. pp%doy ==365)
     &    .or.(pp%cellptr%hemi ==2 .and. pp%doy ==181))then
              ! Call yearly-step subroutines to set cohorts for the next year

!              if(pp%nCohorts<max_cohorts)
!     &              call BiomeE_reprod_seedling (pp)
              call BiomeE_reprod_seedling (pp)
              call BiomeE_Plant_mortality (pp)

              ! Cohort resetting
              call Sorting_cohorts (pp)
              ! Grass proportional thinning
              call BiomeE_GrassThinning(pp)
              call Remove_ghost_cohorts (pp)
              call Merge_cohorts (pp)
              ! Layering
              call Sorting_cohorts (pp)
              call Layering_cohorts (pp)
              !call Sorting_cohorts (pp)

              ! Update max crownLAI
              ! call CrownLAImax_update(pp)

          endif ! Weng's annual vegetation dynamics

          ! update doy
          pp%doy = pp%doy + 1
          if(pp%doy.eq.366) pp%doy = 1
          !When the first day of model run is not Jan. 1
          if(pp%cellptr%doy==365 .and. pp%doy<365)pp%doy=1
      end subroutine BiomeE_FullDemography

!==============================================================================
      subroutine BiomeE_SingleTree(pp)
      !@sum Vegetation dynamics for single tree model
      ! Weng, updated: 11/07/2021
        use ent_types, only: max_cohorts
        implicit none
        type(patch),   intent(inout) :: pp

        !------local var -----------
        type(cohort), pointer :: cop
        integer :: nCohorts, i
        logical :: do_PPA_diag
        do_PPA_diag = .TRUE. ! .FALSE. !

        ! Plant grwoth, mortality, and organization
        call BiomeE_PlantGrowth(pp)
        call BiomeE_turnover_all(pp)
        call BiomeE_turnover_Wood(pp) ! continuous mortality
        call BiomeE_veg_update(pp)
        call BiomeE_selfthinning(pp)

        ! call Ent-PPA diagnostics, including daily
        if(do_PPA_diag)call BiomeE_diagnostics(pp)

        ! ----- Call yearly subroutines at the first day after growing season
        ! use phenology flags ? Weng, 4/26/2018
        !if(ecop%gdd .eq. 0.0 .and. pp%tallest%C_seed > 0.000001) then
        if( (pp%cellptr%hemi ==1 .and. pp%doy ==365)
     &  .or.(pp%cellptr%hemi ==2 .and. pp%doy ==181))then
            ! Reproduction, Self-organization, and Cohort resetting
            call Sorting_cohorts (pp)
            call BiomeE_reprod_samesized(pp)
            !zero patch diagnostics
            pp%N_min_yr = 0.d0
        endif ! annual vegetation dynamics
        ! update patch LAI
        call get_patch_LAI (pp)
        ! update doy
        pp%doy = pp%doy + 1
        if(pp%doy ==366) pp%doy = 1
        !!When the first day of model run is not Jan. 1
        if(pp%cellptr%doy==365 .and. pp%doy<365)pp%doy=365
      end subroutine BiomeE_SingleTree

!==============================================================================

      subroutine BiomeE_patch_init(pp)
!@sum initialize plant cohorts and soil pools
! called in "init_simple_entcell" in entcells.f

        implicit none
        type(patch),   intent(inout) :: pp

        !------local var -----------
        type(cohort), pointer :: cop, newc
        integer :: pft
        real*8  :: C_seed,N_seed
        integer :: EntPFT(10),nCohorts, n, i
        ! Initialize pft parameters
        EntPFT = -999
        !call BiomeE_params_init() ! called in ent_init_params
        ! Initialize soil pools
        pp%mineralN = QN_mineral0 ! 0.2 ! KgN m-2
        pp%SOC(3) = 0.0
        pp%SOC(4) = initial_slowC
        pp%SOC(5) = initial_slowC

        pp%SON(3) = pp%SOC(3)/CN0SOM(3)
        pp%SON(4) = pp%SOC(4)/CN0SOM(4)
        pp%SON(5) = pp%SOC(5)/CN0SOM(5)

        pp%Tpool (:,4:12,:) = 0.d0
        pp%Tpool(1,11,1) = initial_slowC ! kgC m-2
        pp%Tpool(1,12,1) = initial_passiveC  ! kgC m-2
        do n = NLIVE+1,NPOOLS
          do i=1,N_CASA_LAYERS
            pp%Tpool(2,n,i)=pp%Tpool(1,n,i)/CN0(n)
          enddo
        enddo
        ! Remove all existing cohorts
        cop => pp%tallest
        nCohorts = 0
        do while(ASSOCIATED(cop))
            nCohorts = nCohorts + 1
            EntPFT(nCohorts)=cop%pft
            call Remove_a_cohort (pp,cop)
            cop => pp%tallest
        enddo
        pp%tallest => null()
        pp%shortest => null()

! Use Ent veg map PFT
#ifdef Do_VegMap_Ent
        if(nCohorts>0)  then
          write(933,*)'Do_VegMap_Ent',nCohorts,(EntPFT(i),i=1,nCohorts)
          do i=1,nCohorts
             if(EntPFT(i)>0) iniPFT(i) = EntPFT(i)
          enddo
          N_iniCH = nCohorts
        endif
#endif
        ! add new cohorts
        do i =1, N_iniCH
           call Cohort_birth(pp,iniPFT(i),
     &                  iniC(i),iniC(i)/15.0)

!           if(ijdebug==44067)then
!              cop => pp%tallest
!              write(932,'(a50,I9,L4,9(E12.4,","))')
!     &        'initial PFT:,Pheno,n,dbh,h,Acrown,c_lab,sw',
!     &        iniPFT(i),cop%PhenoON,cop%n,cop%dbh,cop%h,
!     &        cop%crownarea,cop%C_lab,cop%C_sw
!           endif

        enddo
        pp%nCohorts = N_iniCH
        call Sorting_cohorts (pp)
      end subroutine BiomeE_patch_init

!==============================================================================
! Add this subroutine to parameter and initial state initialization
      subroutine BiomeE_params_init() !initialize_PPABMEspdata()
        use ent_const
        use ent_pfts, only: BMEspdata,pfpar
      ! -------local variables
        integer :: i,j,n
        real*8 :: lnscl
        integer :: io           ! i/o status for the namelist
        integer :: ierr         ! error code, returned by i/o routines
        integer :: nml_unit

!   read in plant traits from a namelist file
!     !nml_unit = get_unit()  ! check how to get a file unit in Ent
!     !open(nml_unit, file='plant_traits.nml', form='formatted', action='read', status='old')
!     !do
!     !   read (nml_unit, nml=pft_traits_nml, iostat=io, end=10)
!     !   if (check_nml_error (io, 'pft_traits_nml')==0) exit ! from loop
!     !enddo
!! 10  close (nml_unit)

        associate (sp => BMEspdata)
        !-------- from pfpar to BMEspdata -------!
        do i=1, N_PFT
          sp(i)%nf      = pfpar(i)%nf !Canopy nitrogen factor
          !CASA and CLM parameters
          sp(i)%r       = pfpar(i)%r
          sp(i)%lrage   = pfpar(i)%lrage
          sp(i)%woodage = pfpar(i)%woodage
          sp(i)%lit_C2N = pfpar(i)%lit_C2N
          sp(i)%lignin  = pfpar(i)%lignin
          sp(i)%croot_ratio   = pfpar(i)%croot_ratio
        enddo

        ! BiomeE plant species assignment
        do i=1, N_BMEPFT
            j = pftB2E(i)
            sp(j)%woody       = woody(i)
            sp(j)%phenotype   = phenotype(i)
            sp(j)%leaftype    = leaftype(i)

            sp(j)%hwilt       = hwilt(i) ! Wilting point matric potential (m)
            sp(j)%sstar       = sstar(i) ! Rel. soil moist at stress onset (Rodriguez-Iturbe)
            sp(j)%swilt       = swilt(i) ! Normalized soil water at wilting point (dim'less)

            sp(j)%pst         = pst(i)
            sp(j)%LMA         = LMA(i)
            sp(j)%PARabsorb   = PARabsorb(i)
            sp(j)%Vcmax       = Vcmax(i)
            sp(j)%m           = m(i)
            sp(j)%b           = b(i)

            sp(j)%LNbase      = 1.5d0 * Vcmax(i)/1.0d5 ! LNbase(i)
            sp(j)%Nfixrate0   = Nfixrate0(i)
            sp(j)%CNstrucL    = CNstrucL(i)
            sp(j)%CrownLAImax = CrownLAImax(i)
            sp(j)%cLAI_light  = cLAI_light(i)
            sp(j)%H0_cL       = H0_cL(i)

            sp(j)%rhoWood      = rhoWood(i)
            sp(j)%f_taper      = f_taper(i)
            sp(j)%K_cambium    = K_cambium(i)
            sp(j)%f_Acam       = f_Acam(i)
            sp(j)%Acrown0      = Acrown0(i)
            sp(j)%tauNSC       = tauNSC(i)

            sp(j)%rhoFR        = rhoFR(i)
            sp(j)%rootLS       = rootLS(i)
            sp(j)%root_r       = root_r(i)
            sp(j)%root_zeta    = root_zeta(i)
            sp(j)%root_perm    = root_perm(i)

            sp(j)%phiRL       = phiRL(i)
            sp(j)%phiCSA      = phiCSA(i)

            sp(j)%alphaHT     = alphaHT(i)
            sp(j)%thetaHT     = thetaHT(i)
            sp(j)%alphaCA     = alphaCA(i)
            sp(j)%thetaCA     = thetaCA(i)

            sp(j)%thetaBM     = thetaBM(i)

            sp(j)%r_seed      = r_seed(i)
            sp(j)%seedlingC0  = seedlingC0(i) ! kgC / seedling
            sp(j)%mu0         = mu0(i)
            sp(j)%D0mu        = D0mu(i)
            sp(j)%A_sd        = A_sd(i)
            sp(j)%B_sd        = B_sd(i)

            sp(j)%CNroot0 = CNroot0(i)
            sp(j)%CNsw0   = CNsw0(i)
            sp(j)%CNwood0 = CNwood0(i)
            sp(j)%CNseed0 = CNseed0(i)
            ! Phenology
            sp(j)%Tc_pheno0 = Tc_pheno0(i)
            sp(j)%gdd_par1  = gdd_par1(i)
            sp(j)%gdd_par2  = gdd_par2(i)
            sp(j)%gdd_par3  = gdd_par3(i)
            sp(j)%betad_ON  = betad_ON(i)
            sp(j)%betad_OFF = betad_OFF(i)

            ! ---------- derived traits ---------------------------
            sp(j)%LeafLS = MAX(c_LLS*sp(j)%LMA, 1.d0) ! Years
            sp(j)%LNA = sp(j)%LNbase + sp(j)%LMA/sp(j)%CNstrucL
            sp(j)%Nleaf = sp(j)%LNA * 1000.d0 ! !g-N/m2[leaf]
            sp(j)%SRA = 2.d0/(sp(j)%root_r * sp(j)%rhoFR)
            sp(j)%CNleaf0 = sp(j)%LMA / sp(j)%LNA
            sp(j)%alphaBM = PI/4. * 1.25 *
     &          sp(j)%rhoWood * sp(j)%f_taper * sp(j)%alphaHT

        enddo

#ifdef BiomeE_DEBUG
      do i = 1,N_BMEPFT
        j = pftB2E(i)
        write(98,'(A13,I3,A13)') '=======  pft=',j,'  ========='
        write(98,*) 'pst      ',sp(j)%pst
        write(98,*) 'woody    ',sp(j)%woody
        write(98,*) 'phenotype',sp(j)%phenotype
        write(98,*) 'leaftype ',sp(j)%leaftype
        write(98,'(A9,F10.2)') 'Vcmax',sp(j)%Vcmax
        write(98,'(A9,F10.2)') 'H0_cL',sp(j)%H0_cL
        write(98,*)'  '
      enddo
#endif

      end associate

        !----- Initialize output files -----------
        call write_outputfile_headers()

      end subroutine BiomeE_params_init

!================================================================================
!======================Subroutines for BiomeE VegDyn ============================
!================================================================================
      subroutine BiomeE_PlantGrowth(pp)
      ! for plant growth, including carbon fetch and allocation
      ! @Ensheng Weng, 2018-10-31
      implicit none
      type(patch), intent(inout) :: pp

      ! ---- local vars
      type(cohort), pointer :: cop  ! current cohort
      integer   :: nCohorts
      !type(PlantSP),pointer :: sp
      real*8 :: f2Wood_min = 0.1d0  ! minimum allocation to wood biomass from carbon_gain
      real*8 :: G_LFR     ! amount of carbon spent on leaf and root growth
      real*8 :: dC_LF, dC_FR, dC_SW, dC_SD  ! tendencies of growth, kgC/individual
      real*8 :: dN_LF, dN_FR, dN_SW, dN_SD
      real*8 :: dC_LRS    ! = dC_LF+dC_FR+dC_SD
      real*8 :: dDBH,dCA,dHeight   ! tendency of growt
      real*8 :: bl_max,br_max ! maximum leaf and fine root biomass at given CA
      real*8 :: cLAImax ! maximum crown LAI at given height: =CrownLAImax*(h+0.25)/(h+sp%H0_cL)
      real*8 :: r_Nsup, N_supply, N_demand, extraN
      real*8 :: C_push, C_demand, N_push
      real*8 :: fNr ! fraction of NSN per day that are available for growth maximumly
      real*8 :: growthC, resp_growth
      !real*8 :: fNSC_push,fLFR_demand ! make these two variables to PFT-specific parameters
      real*8 :: copC, copN
      ! for size update
      real*8 :: structuralC
      real*8 :: Aleafmax,  CSAsw, CSAtot, CSAwd
      real*8 :: DBHwd, BSWmax, dC_WD, dC_BM
      logical :: PhenoON = .True. !

      ! ----- Fetch carbon from C_lab -------------
      !looping through cohorts and Fetch C from labile C pool if it is in the growing seasaon
      nCohorts =0
      cop => pp%tallest
      do while(ASSOCIATED(cop))
         nCohorts = nCohorts + 1
         ! fectch carbon from labile carbon pool
         associate(sp => BMEspdata(cop%pft))
         cop%crownarea = dbh2Crownarea (cop%pft,cop%dbh)

         PhenoON=(cop%phenostatus==2 .OR. cop%phenostatus==3)

         !! Reset grass density in the first day of a growing season
         !if((PhenoON).and.(.not.cop%PhenoON).and.(.not.sp%woody))then !
         !    copC=btotal(cop) * cop%n
         !    copN=Ntotal(cop) * cop%n
         !    call reset_grass(cop, copC, copN)
         !endif

         ! Update pheno status
         cop%PhenoON = PhenoON

         ! Put C from labile C to C_for_growth
         ! C_lab has excluded resp_grwoth, i.e., only NPP is added to C_lab in respauto_physio.f
!         update leaf age
         if(cop%C_fol>0.d-4)then
            cop%leafage = cop%leafage + 1.d0/365.d0
         else
            cop%leafage = 0.d0
         endif
         ! Growth control variables
         cLAImax = sp%CrownLAImax
     &           * (cop%h+0.25)/(cop%h+sp%H0_cL) ! LAImax vs. Height
     &           / Max(1.d0, (cop%layer-1)*2.d0+1.d0) ! Low LAImax in understory
         bl_max = max(cLAImax*sp%LMA * cop%crownarea, 0.2d-2) ! kgC/tree
         br_max = max(sp%phiRL*bl_max/(sp%LMA*sp%SRA),0.2d-2) ! kgC/tree
         cop%NSCstar = 1.d0 * (bl_max + br_max)          ! kgC/tree

         if (PhenoON .and. cop%C_lab>0.2d0*cop%NSCstar)then ! growing season
             C_demand = (Max(bl_max - cop%C_fol,0.d0) +
     &                   Max(br_max - cop%C_froot,0.d0))/15.d0
             C_push   = sp%tauNSC/365.d0 *
     &              max(cop%C_lab-0.5d0*cop%NSCstar, 0.d0)
             growthC  = Max(MIN(0.03*cop%C_lab,C_demand+C_push),0.d0)

             ! ------------- Growth --------------------------
             ! calculate the carbon allocated for the growth of leaves and roots
             ! and distribute it between fine roots and leaves
             G_LFR=MIN((1.d0 - f2Wood_min) * growthC,
     &             Max(bl_max+br_max-cop%C_fol-cop%C_froot,0.d0))
             dC_LF = min(G_LFR, max(0.d0,
     &          (G_LFR*bl_max+bl_max*cop%C_froot-br_max*cop%C_fol)/
     &          (bl_max + br_max)))
             dC_FR = G_LFR - dC_LF
             dC_SD = sp%r_seed*(growthC - G_LFR)
             dC_SW = growthC - G_LFR - dC_SD

!!             ----Nitrogen effects on allocations-----
             ! fraction of NSN per day that are available for growth maximumly
             fNr = max(growthC/cop%C_lab,0.08d0)
!!             Nitrogen demand by leaves, roots, and seeds (Their C/N ratios are fixed.)
             N_demand = dC_LF/sp%CNleaf0 + dC_FR/sp%CNroot0
     &                + dC_SD/sp%CNseed0 + dC_SW/sp%CNsW0
!!              Nitrogen available for all tisues, including wood
             N_supply= MIN(MAX(0.d0, fNr*cop%N_lab),N_demand)
!!              same ratio reduction for leaf, root, and seed if(N_supply < N_demand)
             if(N_supply < N_demand)then
                r_Nsup = N_supply / N_demand
                dC_LRS = dC_LF+dC_FR+dC_SD
                dC_FR  =  r_Nsup * dC_FR
                dC_LF  =  r_Nsup * dC_LF
                dC_SD  =  r_Nsup * dC_SD
                dC_SW  =  r_Nsup * dC_SW ! dC_SW + (1.d0 - r_Nsup) * dC_LRS
                N_demand = N_supply
             endif
!              update plant pools
             dC_BM = dC_LF + dC_FR + dC_SW + dC_SD

             cop%C_lab   = cop%C_lab   - dC_BM
             cop%C_fol   = cop%C_fol   + dC_LF
             cop%C_froot = cop%C_froot + dC_FR
             cop%C_sw    = cop%C_sw    + dC_SW
             cop%C_seed  = cop%C_seed  + dC_SD
             cop%leafarea= cop%C_fol/sp%LMA ! ??
!              update leaf age
             if(dC_LF>0.d0)
     &         cop%leafage = (1.0 - dC_LF/cop%C_fol) * cop%leafage
!!             update nitrogen pools, Nitrogen allocation
             dN_LF = dC_LF/sp%CNleaf0
             dN_FR = dC_FR/sp%CNroot0
             dN_SD = dC_SD/sp%CNseed0
             dN_SW = N_supply - dN_LF - dN_FR - dN_SD
             cop%N_lab   = cop%N_lab   - N_supply
             cop%N_fol   = cop%N_fol   + dN_LF
             cop%N_froot = cop%N_froot + dN_FR
             cop%N_seed  = cop%N_seed  + dN_SD
             cop%N_sw    = cop%N_sw    + dN_SW
!              Return excessiive Nitrogen in SW back to NSN
             if(cop%N_sw > cop%C_sw/sp%CNsw0)then
                extraN    = cop%N_sw  - cop%C_sw/sp%CNsw0
                cop%N_lab = cop%N_lab + extraN
                cop%N_sw  = cop%N_sw  - extraN
             endif
         else ! non-growing season
             dC_BM = 0.d0
         endif
        end associate
        cop%C_for_growth = dC_BM
        resp_growth      = 0.33d0 * dC_BM

         !* Tissue growth respiration is distributed over day
         !* subtracted at physical time step in canopy biophysics
         !* module with R_auto. The unit is kgC m-2 time-1
         cop%C_growth = resp_growth ! kg-C tree-1 step-1
         cop%R_growth = resp_growth/(24.d0*3600.d0)  ! flux, kg-C tree-1 s-1
!!         for debug
!         if(ijdebug==44067) !101061 !44067
!     &       write(932,'(a35,I9,L4,9(E12.4))')
!     &       'DOY,Pheno,n,c_lab,sw,h,growthC,Rgr,dC_BM',
!     &        pp%doy,cop%PhenoON,cop%n,cop%C_lab,cop%C_sw,
!     &        cop%h,growthC,cop%R_growth,dC_BM

        ! next cohort
        cop  => cop%shorter
      enddo

      end subroutine BiomeE_PlantGrowth

!======================================================================================
!, Weng, 2018-10-29, Updated 2021-10-23
! Normal turnover of leaves, grass stems, and fine roots in a growing season
      subroutine BiomeE_turnover_all(pp)
      implicit none
      type(patch), intent(inout) :: pp

! ---- local vars
      type(cohort), pointer :: cop  ! current cohort
      logical :: dormant
      real*8 :: alpha_LF, alpha_WD, alpha_FR ! turnover rate, fraction per day
      real*8 :: dC_LF, dC_SW, dC_HW, dC_FR   ! Turnover of tissues
      real*8 :: dN_LF, dN_HW, dN_SW, dN_FR   ! N Turnover
      real*8 :: bl_max, L_fall   ! leaf fall at the end of growing season, kgC tree-1 day-1
      real*8 :: f_fall = 0.2d0 ! per day, at the end of a growing season
      real*8 :: dt_day2year = 1.d0/365.d0 ! time step

      cop => pp%tallest
      do while(ASSOCIATED(cop))
         associate(sp => BMEspdata(cop%pft))

         ! Calculate normal turnover of leaves, fine roots, and grass stems
         !alpha_LF = dt_day2year / sp%leafLS ! leave turnover rate, fraction per day
         alpha_LF = MIN(0.2,Max(2.0*cop%leafage/sp%leafLS-1.0,0.0)) ! leave turnover rate, fraction per day
         alpha_FR = dt_day2year / sp%rootLS ! root turnover rate, fraction per day

         ! Normal leaf and fine root turnover
         dC_LF  = cop%C_fol * alpha_LF
         dN_LF = dC_LF/sp%CNleaf0
         dC_FR = cop%C_froot * alpha_FR
         dN_FR = cop%N_froot * alpha_FR

         ! Grass stems turnover with leaves
         if((.not. sp%woody).and.(cop%C_sw>0.2d-2))then
            dC_SW = cop%C_sw * alpha_LF
            dC_HW = cop%C_hw * alpha_LF
            dN_SW = cop%N_sw * alpha_LF
            dN_HW = cop%N_hw * alpha_LF
         else
            dC_SW = 0.; dC_HW = 0.
            dN_SW = 0.; dN_HW = 0.
         endif

         ! Leaf and grass stems fall at the end of a growing season
         dormant = cop%phenostatus==1 .OR. cop%phenostatus==4
         L_fall = 0.
         if(dormant .and. cop%C_fol > 0.d0)then
            ! Leaf fall and update C_lab
            bl_max=max(sp%CrownLAImax*(cop%h+0.25)/(cop%h+sp%H0_cL)
     &            * sp%LMA * cop%crownarea, 0.2d-2) ! kgC/tree
            L_fall = MAX(f_fall * cop%C_fol,
     &               MIN(cop%C_fol, 0.1d0*bl_max))
            ! Update leaf turnover
            dC_LF = MIN(cop%C_fol, dC_LF+L_fall)
            dN_LF = dC_LF/sp%CNleaf0

            ! Grass stem fall, Weng, 2021-01-18
            if(.not. sp%woody)then
               ! Stem falls and Keep a minimum stem biomass for grasses
               if(cop%C_sw>0.2d-2)then
                ! Resorption to NSC and NSN pools
                cop%C_lab = cop%C_lab+ f_fall*cop%C_sw*0.4  ! rhizome
                cop%N_lab = cop%N_lab+ f_fall*cop%N_sw*0.4
                ! Update grass stems
                cop%C_sw  = cop%C_sw - f_fall*cop%C_sw*0.4  ! rhizome
                cop%N_sw  = cop%N_sw - f_fall*cop%N_sw*0.4
                ! Litter
                dC_SW = f_fall * cop%C_sw * 0.6
                dN_SW = f_fall * cop%N_sw * 0.6
               endif
            endif
         endif

         ! Update plant carbon pools
         cop%C_fol   = cop%C_fol   - dC_LF
         cop%C_sw    = cop%C_sw    - dC_SW
         cop%C_hw    = cop%C_hw    - dC_HW
         cop%C_froot = cop%C_froot - dC_FR

         cop%N_fol   = cop%N_fol   - dN_LF
         cop%N_sw    = cop%N_sw    - dN_SW
         cop%N_hw    = cop%N_hw    - dN_HW
         cop%N_froot = cop%N_froot - dN_FR

!         update leaf age
         if(cop%C_fol>0.d-4)then
           cop%leafage = (1.0 - dC_LF/cop%C_fol)*cop%leafage
         else
           cop%leafage = 0.d0
         endif

!         ! Put grass roots back to labile pools at the end of a growing season
!         if((.not. sp%woody)
!     &       .and. cop%phenostatus==4
!     &       .and. cop%C_froot > 0.d0)then
!                ! Put root carbon and nitrogen to labile pools
!                cop%C_lab = cop%C_lab + cop%C_froot
!                cop%N_lab = cop%N_lab + cop%N_froot
!                cop%C_froot = 0.d0
!                cop%N_froot = 0.d0
!         endif

!         ! Update soil C and N
!         ! Put the dead tissue's carbon to soil pools, pp%Tpool (1,12,1)
         call turn2SoilBiomeE(pp, cop,
     &       dC_LF, dC_SW + dC_HW, dC_FR,
     &       dN_LF, dN_SW + dN_HW, dN_FR)

         end associate

         cop => cop%shorter
      enddo

      end subroutine BiomeE_turnover_all

!======================================================================================
! Wood turnover, for single cohort modeling as replacement of yearly mortality
! Weng, 10-23-2021
      subroutine BiomeE_turnover_Wood(pp)
      implicit none
      type(patch), intent(inout) :: pp

! ---- local vars
      type(cohort), pointer :: cop  ! current cohort
      real*8 :: alpha_WD ! turnover rate, fraction per day
      real*8 :: dC_SW, dC_HW, dN_HW, dN_SW  ! Turnover of tissues
      real*8 :: dt_year = 1.d0/365.d0    ! time step

      cop => pp%tallest
      do while(ASSOCIATED(cop))
         associate(sp => BMEspdata(cop%pft))
         if(sp%woody) then ! for woody plants
            cop%mu = mortality_rate(cop) ! per year
            alpha_WD = cop%mu * dt_year
            ! Normal Wood turnover
            dC_SW = cop%C_sw * alpha_WD
            dN_SW = cop%N_sw * alpha_WD
            dC_HW = cop%C_hw * alpha_WD
            dN_HW = cop%N_hw * alpha_WD

            ! Update plant carbon pools
            cop%C_sw    = cop%C_sw - dC_SW
            cop%C_hw    = cop%C_hw - dC_HW
            cop%N_sw    = cop%N_sw - dN_SW
            cop%N_hw    = cop%N_hw - dN_HW

            ! Put the dead tissue's carbon to soil pools, pp%Tpool (1,12,1)
            call turn2SoilBiomeE(pp, cop,
     &         0.d0, dC_SW + dC_HW, 0.d0,
     &         0.d0, dN_SW + dN_HW, 0.d0)
         endif
         end associate
         ! next cohort
         cop => cop%shorter
      enddo

      end subroutine BiomeE_turnover_Wood
!==============================================================================
! Update crown area, height, dbh, LAI, and conversion of heartwood
      subroutine BiomeE_veg_update(pp)
      implicit none
      type(patch), intent(inout) :: pp

! ---- local vars
      type(cohort), pointer :: cop  ! current cohort
      integer   :: nCohorts, nLayer
      !real*8 :: dn ! changes in density
      real*8 :: structuralC
      real*8 :: Aleafmax, bl_max, CSAsw, CSAtot, CSAwd
      real*8 :: DBHwd, BSWmax, dC_WD, dN_WD

      ! Update plan structure
      nCohorts =0
      pp%LAIlayer = 0.d0
      cop => pp%tallest
      do while(ASSOCIATED(cop))
         nCohorts = nCohorts + 1
         associate (sp => BMEspdata(cop%pft))

!        update breast height diameter given increase of bsw
         call treeBM2structure(cop)
         call prescr_calc_rootprof(cop%fracroot,cop%pft) ! added by Weng, 09/09/2021
         nLayer = MIN(N_CAImax,cop%layer)
         pp%LAIlayer(nLayer)=pp%LAIlayer(nLayer)+
     &                          cop%leafarea * cop%n

!        conversion of sapwood to heartwood
         if(sp%woody .and. cop%C_sw > 1.0d-2)then
           Aleafmax = sp%CrownLAImax * cop%crownarea
           CSAsw    = Aleafmax * sp%phiCSA ! * cop%h
           CSAtot   = PI * (cop%DBH/2.d0)**2
           CSAwd    = max(0.d0, CSAtot - CSAsw)
           DBHwd    = 2.d0 * sqrt(CSAwd/PI)
           BSWmax   = sp%alphaBM *
     &          (cop%DBH**sp%thetaBM - DBHwd**sp%thetaBM)
           dC_WD   = max(cop%C_sw - BSWmax, 0.d0)
           dN_WD   = cop%N_sw * dC_WD/cop%C_sw
           ! Update wood C and N pools
           cop%C_hw  = cop%C_hw + dC_WD
           cop%C_sw  = cop%C_sw - dC_WD
           cop%N_hw  = cop%N_hw + dN_WD
           cop%N_sw  = cop%N_sw - dN_WD
         endif

         end associate
         ! next cohort
         cop  => cop%shorter
      enddo

      end subroutine BiomeE_veg_update

!======================================================================================
      subroutine BiomeE_selfthinning(pp)
      ! Selfthinning removes the individuals that cannot be hold.
      ! Daily time step.
      implicit none
      type(patch),intent(inout) :: pp

      !-----local--------
      type(cohort),pointer :: cop
      real*8  :: totalCA, CAcrown, dn
      integer :: L = 1

      totalCA  = 0.d0
      cop => pp%tallest
      do while(ASSOCIATED(cop))
         associate(sp => BMEspdata(cop%pft))
         CAcrown = Max((1.d0+fgap)*cop%crownarea,
     &             sp%Acrown0)
         totalCA = totalCA + cop%n * CAcrown
         cop%layer = L
         if(totalCA .gt. L )then
             ! selfthinning
             dn = MIN(cop%n,(totalCA - L) / CAcrown)
             ! Kill some trees in this cohort to let totalCA = L
             cop%n = cop%n - dn
             cop%LAI = cop%leafarea * cop%n
             totalCA = totalCA - dn * CAcrown
             ! Put the dead trees carbon to soil pools, pp%Tpool (1,12,1)
             call BiomeE_Plant2Soil(cop,pp,dn)
         endif !'totalCA .gt. L'
         end associate
         ! go to the next cohort
         cop  => cop%shorter
      enddo
      pp%totalCA  = totalCA

      end subroutine BiomeE_selfthinning

!======================================================================================
      subroutine BiomeE_GrassThinning(pp)
      ! Proportionally reduce each grass cohort density if grass cohorts total
      ! crown area is greater than land area
      ! Added by Weng, 2021-11-16
      implicit none
      type(patch),intent(inout) :: pp
      !-----local--------
      type(cohort),pointer :: cop
      real*8 :: GrassCA, totalCA, dn
      real*8 :: r_rem ! Ratio of grasses should be removed
      real*8 :: GCmax = 1.d0 ! Maximum grass coverage

      ! Grass cohorts and area calculation
      GrassCA = 0.d0
      cop => pp%tallest
      do while(ASSOCIATED(cop))
         associate(sp => BMEspdata(cop%pft))
         if(.not. BMEspdata(cop%pft)%woody)then
            GrassCA = GrassCA + cop%n *
     &         Max((1.d0+fgap)*cop%crownarea,sp%Acrown0)
         endif
         end associate
         ! go to the next cohort
         cop  => cop%shorter
      enddo

      ! Proportionally thinning
      if(GrassCA > GCmax)then
        r_rem = (GrassCA - GCmax) / GrassCA
        totalCA = 0.d0
        cop => pp%tallest
        do while(ASSOCIATED(cop))
          if(.not. BMEspdata(cop%pft)%woody)then
            dn = r_rem * cop%n
            cop%n = cop%n - dn
            cop%LAI = cop%leafarea * cop%n
            call BiomeE_Plant2Soil(cop,pp,dn)
          endif
          totalCA = totalCA + Max((1.0+fgap)*cop%crownarea,
     &             BMEspdata(cop%pft)%Acrown0) * cop%n

          ! go to the next cohort
          cop  => cop%shorter
        enddo
        pp%totalCA  = totalCA
      endif

      end subroutine BiomeE_GrassThinning

!=================================================================
! updated according to Weng's allometry 11-21-2016
      subroutine reset_grass(cop, copC,copN)
      ! Weng: reset grass density at the first day of a growing season
      use ent_pfts, only: COVEROFFSET
      type(cohort),intent(inout) :: cop
      real*8, intent(in) :: copC,copN
      !--local----------------
      real*8 :: leafarea
      real*8 :: structuralC,seedlingN
      integer :: pft

      ! new grass seedlings
      pft = cop%pft
      associate(sp => BMEspdata(cop%pft) )
      if((.not.sp%woody).and.copC>0.0.and.copN>0.0)then
       cop%n = copC/sp%seedlingC0
       cop%C_lab = 0.8 * sp%seedlingC0
       structuralC = 0.15 * sp%seedlingC0

       cop%dbh   = structuralC2dbh (pft,structuralC)
       cop%h     = dbh2Height (pft,cop%dbh)
       cop%crownarea = dbh2Crownarea (pft,cop%dbh)
       cop%crown_dx  = dbh2CrownRadius (pft,cop%dbh)
       cop%crown_dy  = crown_radius_vert(pft,cop%h,cop%crown_dx)

       cop%C_sw  = structuralC
       cop%C_hw  = 0.0
       cop%C_fol = 0.025 * sp%seedlingC0
       cop%C_croot = 0.0 ! 0.25 * (cop%C_sw + cop%C_hw)
       cop%C_froot = 0.025 * sp%seedlingC0
       cop%C_seed  = 0.0 ! Must be zero!!

       ! N pools
       seedlingN = copN * sp%seedlingC0/copC
       cop%N_sw  = MIN(0.5*seedlingN,structuralC/sp%CNsw0)
       cop%N_hw  = 0.0
       cop%N_fol = 0.025 * sp%seedlingC0/sp%CNleaf0
       cop%N_croot = 0.0 ! 0.25 * (cop%C_sw + cop%C_hw)
       cop%N_froot = 0.025 * sp%seedlingC0/sp%CNroot0
       cop%N_seed = 0.0
       cop%N_lab  = seedlingN-(cop%N_sw+cop%N_fol+cop%N_froot)

       cop%leafarea = cop%C_fol/sp%LMA
       cop%leafage = 0.d0
       cop%LAI      = cop%leafarea * cop%n
       cop%root_d = sp%root_r  ! Check if "root_d" is root radius !!!!
       !call rootprofile(cop)
       call prescr_calc_rootprof(cop%fracroot,cop%pft+COVEROFFSET)
       !cop%fracroot(:) = fracroot(:)  !?
       cop%Ci = 0.0 !Ci
       cop%GCANOPY =  0.0
       cop%GPP =  0.0
       cop%NPP =  0.0
       cop%R_auto = 0.0
       cop%R_root = 0.0
       cop%N_uptake =  0.0
       cop%N_up_yr =  0.0
       cop%C_to_Nfix = 0.0
       cop%betad_10d = 0.0
       cop%CB_d = 0.0
       cop%llspan = sp%leafLS
      endif
      end associate
      end subroutine reset_grass

!======================================================================================
      subroutine BiomeE_reprod_seedling(pp)
! Reproductions of cohorts
! The idea is to define an allocatable array, put each PFT's seeds into these array,
!  and then generate new seedling cohorts
! After call BiomeE_reprod_seedling, one must immediately call cohort_merging in case an
! cohort with cop%n=0 is generated.
      type(patch),intent(inout) :: pp
      !-----local--------
      type(cohort),pointer :: cop
      integer :: i,j,PFT,matchflag
      integer :: nCohorts,nFC ! number of PFTs that have fecundity carbon
      integer, dimension(N_PFT) :: ParentPFT
      real*8,dimension(N_PFT,2) :: seedpool ! 1 seed Carbon , 2 seed Nitrogen
      real*8 :: copC, copN, minCseed = 1.d-3 ! KgC, minimum C in seeds for reproduction
      logical :: reproducible

      ! Reset grass when pheno is OFF (Tropical grasses, which might be phenoON,
      ! will generate seedlings in the next step)
      cop => pp%tallest
      do while(ASSOCIATED(cop))
         associate(sp => BMEspdata(cop%pft) )
         if((.not.cop%PhenoON).and.(.not.sp%woody))then !
             copC=btotal(cop) * cop%n
             copN=Ntotal(cop) * cop%n
             call reset_grass(cop, copC, copN)
         endif
         end associate
        cop => cop%shorter
      enddo

!     Check how many PFTs are in these cohorts and
!     put C_seed of each cohort into PFT-specific Seed pool
      nFC = 0 ! number of PFTs that have cohorts with C_seed>minCseed
      seedpool = 0.d0
      cop => pp%tallest
      ! Record the PFT of the first cohort
      if(ASSOCIATED(cop))then
        minCseed = 1.d-2 * BMEspdata(cop%pft)%seedlingC0
        if(cop%C_seed >= minCseed)then
           ParentPFT(1) = cop%PFT
           nFC = nFC + 1
        endif
      endif

      ! Looping through cohorts and put their C_seed to the seed pool
      do while(ASSOCIATED(cop))
        minCseed = 1.d-2 * BMEspdata(cop%pft)%seedlingC0
        reproducible = (cop%C_seed >= minCseed)
!     &     .and.((BMEspdata(cop%pft)%woody.and.cop%age >4.d0)
!     &      .or. (BMEspdata(cop%pft)%woody.and.cop%age >0.d0))
        if(reproducible)then
          matchflag = 0
          do i=1,nFC
            if(cop%PFT == ParentPFT(i))then
               seedpool(i,1)=seedpool(i,1) + cop%C_seed*cop%n
               seedpool(i,2)=seedpool(i,2) + cop%N_seed*cop%n
               cop%C_seed = 0.0
               cop%N_seed = 0.0
               matchflag = 1
               exit
            endif
          enddo
          if(matchflag==0)then ! a new PFT
            nFC = nFC+1
            ParentPFT(nFC)  = cop%PFT
            seedpool(nFC,1) = seedpool(i,1)+cop%C_seed*cop%n ! seed carbon
            seedpool(nFC,2) = seedpool(i,2)+cop%N_seed*cop%n ! seed nitrogen
            cop%C_seed = 0.0
            cop%N_seed = 0.0
          endif
        endif
        pp%shortest => cop ! make sure it points to the last one
        cop => cop%shorter
      enddo

      ! from seeds to seedling cohorts
      do i=1,nFC
         PFT = ParentPFT(i)
         call Cohort_birth(pp,PFT,seedpool(i,1),seedpool(i,2))
      enddo

      ! Cohorts counting
      if(nFC >= 1)then
        nCohorts = 0
        cop => pp%tallest
        do while(ASSOCIATED(cop))
          nCohorts = nCohorts + 1
          pp%shortest => cop
          cop  => cop%shorter
        enddo
        pp%nCohorts = nCohorts
      endif

      end subroutine BiomeE_reprod_seedling

!=================================================================
      subroutine Cohort_birth(pp,pft,C_seed,N_seed)  ! Weng, 11-08-2016
!   create one new seedling cohort of pft on a patch (pp) based
!   on the PFT information and seed carbon
!   and then put it in the appropriate place in the chain of cohorts
      type(patch),intent(inout) :: pp
      integer,intent(in) :: pft
      real*8, intent(in) :: C_seed,N_seed
      !--local----------------
      type(cohort),pointer :: newc
      real*8 :: totCseedling

      ! Construct, Initialize, and attach to the patch
      call cohort_construct(newc, pp, pft)
      call initialize_seedling(pp,newc,pft,C_seed,N_seed)
      call attach_seedling_cohort(pp,newc) ! Put it into the pathch and sort
      nullify(newc) ! newc => null()
      end subroutine Cohort_birth

!=================================================================
      subroutine attach_seedling_cohort(pp,newc)
!@sum Weng, 11-09-2016
! put the new seedling cohort to the last
! The patch must have parent cohorts
      type(patch),intent(inout) :: pp
      type(cohort),pointer,intent(inout):: newc
      !--local----------------
      type(cohort),pointer :: cop

      ! if this is the only one
      if(.not. ASSOCIATED(pp%tallest)) then
         pp%tallest => newc
         pp%shortest => newc
         newc%chID = 1
         newc%taller => null()
         newc%shorter => null()
         newc%csptaller => null()
         newc%cspshorter => null()
      else ! Put the new cohort to the end of the link
         newc%chID = pp%shortest%chID + 1
         pp%shortest%shorter => newc
         newc%taller => pp%shortest
         newc%shorter => null()
         newc%cspshorter => null()
         ! update the pp%shortest
         pp%shortest => newc
         ! update the taller co-species cohort
         newc%csptaller  => null()
         cop => newc%taller
         do while (ASSOCIATED(cop))
          if (cop%pft .eq. newc%pft) then
             newc%csptaller => cop
             cop%cspshorter => newc
             cop => null()
             exit !exit loop
          else
             cop => cop%taller
          endif
         enddo
      endif
      pp%nCohorts = pp%nCohorts + 1

      end subroutine attach_seedling_cohort

!=================================================================
! updated according to Weng's allometry 11-21-2016
      subroutine initialize_seedling (pp,cop,pft,C_seed,N_seed)
      ! Weng: initialize a new seedling cohort with information of PFT, 11-09-2016
      ! must define a PFT_data to contain the information of all PFTs
      ! To do: replace this subroutine with a parameter array of seedling cohorts (Type cohort)
      use ent_pfts, only: COVEROFFSET
      type(patch),intent(inout) :: pp
      type(cohort),pointer :: cop
      integer,intent(in) :: pft
      real*8, intent(in) :: C_seed,N_seed
      !--local----------------
      type(cohort),pointer :: cp
      real*8 :: leafarea
      real*8 :: structuralC,seedlingN

      ! Copy parent values to the new cohort
      cp => pp%tallest
      do while(ASSOCIATED(cp))
        if(cp%pft == cop%pft)then
           call copy_cohort(cp, cop)
           exit
        endif
        cp =>cp%shorter
      enddo
      ! Assign specific values
      associate(sp => BMEspdata(pft) )
      ! new seedling cohort
      cop%pft = pft
      cop%n = C_seed/sp%seedlingC0
      cop%age = 0.d0
      cop%leafage = 0.d0
      cop%C_lab = 0.8 * sp%seedlingC0
      structuralC = 0.15 * sp%seedlingC0

      cop%dbh   = structuralC2dbh (pft,structuralC)
      cop%h     = dbh2Height (pft,cop%dbh)
      cop%crownarea = dbh2Crownarea (pft,cop%dbh)
      cop%crown_dx  = dbh2CrownRadius (pft,cop%dbh) !crown_dx
      cop%crown_dy  = crown_radius_vert(pft, cop%h,cop%crown_dx)

      cop%C_sw  = structuralC
      cop%C_hw  = 0.0
      cop%C_fol = 0.025 * sp%seedlingC0
      cop%C_croot = 0.0 ! 0.25 * (cop%C_sw + cop%C_hw)
      cop%C_froot = 0.025 * sp%seedlingC0
      cop%C_seed  = 0.0

      ! N pools
      seedlingN = N_seed/C_seed * sp%seedlingC0
      cop%N_sw  = MIN(0.5*seedlingN, structuralC/sp%CNsw0)
      cop%N_hw  = 0.0
      cop%N_fol = 0.025 * sp%seedlingC0/sp%CNleaf0
      cop%N_croot = 0.0 ! 0.25 * (cop%C_sw + cop%C_hw)
      cop%N_froot = 0.025 * sp%seedlingC0/sp%CNroot0
      cop%N_seed = 0.0
      cop%N_lab  = seedlingN - (cop%N_sw+cop%N_fol+cop%N_froot)

      cop%leafarea = cop%C_fol/sp%LMA
      cop%LAI      = cop%leafarea * cop%n
      cop%root_d = sp%root_r  ! Check if "root_d" is root radius !!!!
      !call rootprofile(cop)
      call prescr_calc_rootprof(cop%fracroot,cop%pft + COVEROFFSET)
      !cop%clump = clump
      cop%LMA =  sp%LMA
      cop%Ci = 0.0 !Ci
      cop%GCANOPY =  0.0
      cop%GPP =  0.0
      cop%NPP =  0.0
      cop%R_auto = 0.0
      cop%R_root = 0.0
      cop%N_uptake =  0.0
      cop%N_up_yr =  0.0
      cop%C_to_Nfix = 0.0
      cop%llspan = sp%leafLS
      !Soil water stress !???
!      if(pp%nCohorts>=1)
!     &      cop%stressH2Ol(:) = pp%tallest%stressH2Ol(:) ! Wrong here!!

      end associate
      end subroutine initialize_seedling

!==========================================================================
      subroutine BiomeE_reprod_samesized(pp)
      ! It only generates the exact same individuals as parents
      ! This is a lumped module for reproduction . We do NOT simulate
      ! the reproduction of seedlings in this module.
      ! Annual time step.
      implicit none
      type(patch),intent(inout) :: pp

      !-----local--------
      type(cohort),pointer :: cop
      integer :: i,j,L
      integer :: matchflag
      integer :: STAT
      integer :: nCohorts=0
      integer :: PFT = -999
      real*8  :: plantC, plantN, n_new
      real*8  :: N_demand, N_left

      !--------- Reproduction ---------------
      !looping through cohorts and convert their C_seed to plants
      !by increasing cop%n
      cop => pp%tallest
      do while(ASSOCIATED(cop))
         !Carbon content of a current individual of this cohort
         plantC = cop%C_lab+cop%C_fol+cop%C_sw+cop%C_hw+cop%C_froot
         plantN = cop%N_lab+cop%N_fol+cop%N_sw+cop%N_hw+cop%N_froot
         n_new = cop%C_seed*cop%n/plantC
         N_demand = n_new * plantN
         N_left = cop%N_seed*cop%n - N_demand

         ! Update density and seed pools
         if(cop%C_seed>0.0)
     &      cop%n  = cop%n + cop%C_seed*cop%n/plantC
         cop%N_lab = cop%N_lab + N_left/cop%n ! put the left N back to C_lab
         cop%C_seed = 0.0
         cop%N_seed = 0.0

         cop => cop%shorter
      enddo

      end subroutine BiomeE_reprod_samesized

!====================================================================
      subroutine BiomeE_Plant_mortality (pp)
! Mortality reduces cohort plant density, operated at patch and called yearly.
! After call BiomeE_reprod_seedling, one must immediately call cohort_merging
! in case a cohort with cop%n=0 is generated. Weng, 2016-12-01
! Updated: 10/15/2021, Weng
      type(patch),intent(inout) :: pp
      !-----local--------
      type(cohort),pointer :: cop
      real*8  :: dt_year ! time stem normalized to a year
      real*8  :: dn ! the number of dead trees due to mortality

      dt_year = 1.d0 ! If the time step is daily, it is 1.0/365.0
      cop => pp%tallest
      do while(ASSOCIATED(cop))
         associate(sp => BMEspdata(cop%pft))
         ! Cohort-specific mortality rate
         cop%mu = mortality_rate(cop)
         dn = (1. - exp(-cop%mu*dt_year)) * cop%n
         cop%n  = cop%n - dn
         cop%LAI = cop%leafarea * cop%n

!         Put the dead trees' carbon to soil pools, pp%Tpool (1,12,1)
         call BiomeE_Plant2Soil(cop,pp,dn)

         end associate

         cop => cop%shorter
      enddo
      end subroutine BiomeE_Plant_mortality
!============================================================
      real*8 function mortality_rate(cop) result(mu)
      !@sum calculate cohort mortality, Ensheng Weng, 11/07/2021
      type(cohort),intent(inout) :: cop

      !-------local var -------------
      real*8 :: f_L, f_S, f_D, expD
      real*8 :: m_S ! Mortality multifactor for size effects
      real*8 :: A_D ! Sensitivity to dbh

      m_S = 3.d0
      A_D = 8.d0
      associate(sp => BMEspdata(cop%pft))
      !mu = sp%mu0*(1.0 + 3.0*(cop%layer-1))
      if(sp%woody) then ! Trees
         f_L =  SQRT(cop%layer - 1.d0) ! Layer effects (0~ infinite)
         f_S = 1. + sp%A_sd * exp(sp%B_sd*cop%dbh) ! Seedling mortality

         expD = exp(A_D * (cop%dbh - sp%D0mu))
         f_D = 1. + m_S * expD / (1.d0 + expD) ! Size effects (big D)

         mu =  Min(0.5d0, sp%mu0*(1.d0+f_L*f_S)*f_D) ! per year
       else ! Grass
         mu =  sp%mu0
       endif
       end associate
      end function mortality_rate

!===========================================================
      subroutine BiomeE_Plant2Soil(cop,pp,dn)
!@sum put a whole trees carbon and nitrogen to soil after mortality
!@ author: Ensheng Weng, 02/27/2021
      type(patch), intent(inout) :: pp
      type(cohort),intent(in) :: cop
      real*8,intent(in) :: dn ! number of dead trees

        ! Put the dead trees' carbon to litter pools
        pp%SOC(1) = pp%SOC(1)
     &           + (cop%C_fol + cop%C_froot
     &           + cop%C_seed + cop%C_lab) * dn
        pp%SOC(2) = pp%SOC(2) + (cop%C_sw+cop%C_hw)*dn

        ! Put the dead trees' nitrogen to litter pools
        pp%SON(1) = pp%SON(1)
     &            + (cop%N_fol + cop%N_froot
     &            + cop%N_seed + cop%N_lab) * dn
        pp%SON(2) = pp%SON(2) + (cop%N_sw+cop%N_hw)*dn

      end subroutine BiomeE_Plant2Soil

!=================================================================================
      subroutine turn2soilBiomeE(pp,cop,
     &       dC_LF,dC_WD,dC_FR,
     &       dN_LF,dN_WD,dN_FR)
!@sum put plant turnover C and N to soil pools
!@ author: Ensheng Weng, 02/27/2021

      type(patch), intent(inout) :: pp
      type(cohort),intent(in) :: cop  ! current cohort
      real*8,intent(in) :: dC_LF  ! Turnover of leaves
      real*8,intent(in) :: dC_WD  ! turnover of wood, represents mortality
      real*8,intent(in) :: dC_FR  ! Turnover of fine roots
      real*8,intent(in) :: dN_LF  ! N Turnover of leaves
      real*8,intent(in) :: dN_WD  ! N Turnover of stem
      real*8,intent(in) :: dN_FR  ! N Turnover of fine roots

      ! Put the dead tissue's carbon & nitrogen to litter pools
      pp%SOC(1) = pp%SOC(1) + (dC_LF + dC_FR) * cop%n
      pp%SOC(2) = pp%SOC(2) + dC_WD * cop%n

      pp%SON(1) = pp%SON(1) + (dN_LF+dN_FR) * cop%n
      pp%SON(2) = pp%SON(2) + dN_WD * cop%n

      end subroutine turn2SoilBiomeE

!==================================================================
      subroutine BiomeE_soil_BGC(dtsec,pp,update_day)
!@auth  Weng, 02/28/2021, Soil Biogeochemical processes
!@ Adapted from Weng et al. 2017 Global Change Biology.
       real*8, intent(in) :: dtsec
       type(patch), intent(inout) :: pp
       logical, intent(in) :: update_day
       !-------- local var ----------
       real*8 :: dt_yr ! dtsec/seconds_in_a_year
       real*8 :: Soilmoist(N_CASA_LAYERS) !(vol fraction) by CASA layers.
       real*8 :: Soiltemp(N_CASA_LAYERS)  !(C) by CASA layers.
       real*8 :: d_C(5), d_N(5), newM(5),N_min_tot
       real*8 :: runoff ! kg m-2 /step
       real*8 :: minC = 1.d-6 ! Minimum C changes
       real*8 :: N_loss
       real*8 :: DON_fast,DON_slow,d_DON ! Dissolved organic N loss, kg N m-2 step-1
       real*8 :: A  ! decomp rate reduction due to moisture and temperature
       integer :: i

       do i=1, 5
          if(isnan(pp%SOC(i))) pp%SOC(i)=0.d0
          if(isnan(pp%SON(i))) pp%SON(i)=0.d0
          pp%SON(i)=MIN(pp%SON(i),2.d0*pp%SOC(i)/CN0SOM(i))
       enddo

       ! do nothing if no vegetation
       if (.not. ASSOCIATED(pp%tallest)) return

       dt_yr = dtsec/(365.*24.*3600.)

       call Soillayer_convert_Ent(pp%cellptr%Soilmoist(:),
     &      SOILDEPTH_m, Soilmoist)
       Soilmoist(:) = Soilmoist(:)*pp%cellptr%soil_Phi !Convert from rel.sat. to vol fraction.
       call Soillayer_convert_Ent(pp%cellptr%Soiltemp(:),
     &     SOILDEPTH_m, Soiltemp)
       runoff = 0.0 ! pp%runoff  !* dt_yr !kgH2O m-2 yr-1 ->kgH2O m-2/time step

       ! Environmental scalar
       A = Max(0.d0, A_function(Soiltemp(1),soilmoist(1)))
       if(isnan(A)) A = 0.d0

       ! Convert litters to SOM
       do i=1,2
          d_C(i)      = pp%SOC(i) * K0SOM(i) * dt_yr
          d_N(i)      = pp%SON(i) * K0SOM(i) * dt_yr
          if(isnan(d_C(i))) d_C(i)=0.d0
          if(isnan(d_N(i))) d_N(i)=0.d0
          pp%SOC(i)   = pp%SOC(i)   - d_C(i)
          pp%SON(i)   = pp%SON(i)   - d_N(i)
          pp%SOC(i+3) = pp%SOC(i+3) + d_C(i)
          pp%SON(i+3) = pp%SON(i+3) + d_N(i)
       enddo
       pp%Soil_resp = 0.d0

       !! Decomposition of three SOM pools, check NAN
       !d_C(i) = pp%SOC(i) * (1.0 - exp(-A*K0SOM(i) * dt_yr))
       do i=3,5
          d_C(i) = pp%SOC(i) * A * K0SOM(i) * dt_yr
          d_N(i) = pp%SON(i) * A * K0SOM(i) * dt_yr
          if(d_C(i) <= minC)then
             d_C(i) = 0.d0
             d_N(i) = 0.d0
          endif
       enddo

       ! New Microbes from pools 4 and 5
       newM(4) = MIN(d_C(4)*CUEf0, d_N(4)*CN0SOM(3))
       newM(5) = MIN(d_C(5)*CUEs0, d_N(5)*CN0SOM(3))
       newM(3) = (1.d0-f_M2SOM) * (newM(4)+newM(5))

       ! update C and N pools
       pp%SOC(3) = pp%SOC(3) - d_C(3) + newM(3)
       pp%SOC(4) = pp%SOC(4) - d_C(4) + f_M2SOM*newM(4)
       pp%SOC(5) = pp%SOC(5) - d_C(5) + f_M2SOM*newM(5)

       ! Update Nitrogen pools
       pp%SON(3) = pp%SON(3) - d_N(3) + newM(3)/CN0SOM(3)
       pp%SON(4) = pp%SON(4) - d_N(4) + f_M2SOM*newM(4)/CN0SOM(3) ! - DON_fast
       pp%SON(5) = pp%SON(5) - d_N(5) + f_M2SOM*newM(5)/CN0SOM(3) ! - DON_slow

       ! Mineralized N
       N_min_tot =  d_N(3)+d_N(4)+d_N(5)-(newM(4)+newM(5))/CN0SOM(3)

       ! Heterotrophic respiration: decomposition of litters and SOM, kgC m-2 s-1
       pp%Soil_resp = (d_C(3)+d_C(4)+d_C(5)-newM(4)-newM(5)) / dtsec

       !! DON loss, revised by Weng. 2016-03-03  ??
       !! Assume it is proportional to decomposition rates (papers?)
       !DON_fast = fDON * d_C(4)*NCfast * (etaN*runoff)
       !DON_slow = fDON * d_C(5)*NCslow * (etaN*runoff)
       !d_DON = DON_fast + DON_slow


       !N_loss = pp%mineralN * A*K_n*dt_yr ! Turned off
       ! N_loss = MAX(0.,pp%mineralN) * (1. - exp(0.0 - etaN*runoff - A*K_n*dt_yr))
       N_loss = 0.d0 ! pp%mineralN*MIN(0.2,(A * K_n * dt_yr + etaN*runoff))

       ! Update mineral N pool
       pp%mineralN = pp%mineralN + N_dpst*dt_yr - N_loss + N_min_tot

       ! Update respiration and N mineralization
       pp%Rh_daily  = pp%Rh_daily + pp%Soil_resp * dtsec
       pp%N_min_day = pp%N_min_day + N_min_tot
       pp%N_min_yr  = pp%N_min_yr  + N_min_tot

       ! Remove NAN
       do i=1, 5
          if(isnan(pp%SOC(i))) pp%SOC(i)=0.d0
          if(isnan(pp%SON(i))) pp%SON(i)=0.d0
       enddo

       ! Update CASA pools (Tpool)
       pp%Tpool(1,(/4,5,10,11,12/),1) = pp%SOC(:)
       pp%Tpool(2,(/4,5,10,11,12/),1) = pp%SON(:)

#ifdef saturatedN
       ! Fix mineral N for test
       pp%mineralN = 50.0d-3 ! kgN m-2 ! Mineral N saturated, just for testing
#endif

       !! Plant N uptake, changed to hourly, Weng, 2021-05-25
       call BiomeE_N_uptake(pp,dtsec,Soiltemp(1),soilmoist(1))

      end subroutine BiomeE_soil_BGC

!=====================================================
! Weng, 2021-05-25
      subroutine BiomeE_N_uptake (pp, dtsec, tsoil, theta)
      type(patch), intent(inout) :: pp
      real*8, intent(in) :: dtsec
      real*8, intent(in) :: tsoil ! average temperature of soil, C
      real*8, intent(in) :: theta ! average soil wetness, unitless

! ---- local vars
      type(cohort), pointer :: cop  ! current cohort pointer
      real*8  :: seconds_per_hour = 3600.0
      real*8  :: rho_N_up0 = 1.d-2 ! maximum N uptake rate (hourly), fraction of the total mineral N
      real*8  :: N_roots0  = 0.1  ! root biomass at half max N-uptake rate,kg C m-2
      real*8  :: fnsnmax   = 0.2 ! a parameter for NSNmax
      real*8  :: NSNmax
      real*8  :: N_up, totNup, avgNup    ! kgN m-2
      real*8  :: N_roots   ! root density
      integer :: i

      !! Nitrogen uptake parameter
      ! It considers competition here. How much N one can absorp depends on
      ! how many roots it has and how many roots other individuals have.

        pp%N_uptake = 0.d0
        if(pp%mineralN > 0.d0)then
           N_Roots  = 0.d0
           cop => pp%tallest
           do while(ASSOCIATED(cop))
              NSNmax = fnsnmax * cop%C_lab
              cop%N_uptake = 0.d0 ! In case N_roots = 0.d0
#ifdef saturatedN
              cop%N_uptake = Max(0.d0, NSNmax - cop%C_lab)
              cop%N_lab    = NSNmax
              cop%N_up_yr  = cop%N_up_yr + cop%N_uptake
              ! patch nitrogen update
              N_up = cop%N_uptake * cop%n
              pp%N_uptake  = pp%N_uptake  + N_up
              pp%N_up_day  = pp%N_up_day  + N_up
              pp%N_up_yr   = pp%N_up_yr   + N_up
#endif

              if(cop%N_lab < NSNmax)
     &              N_Roots = N_Roots + cop%C_froot*cop%n

              cop => cop%shorter
           enddo

           ! M-M equation for Nitrogen absoption, McMurtrie et al. 2012, Ecology & Evolution
           ! rate at given root biomass and period of time
           if(N_roots>1.d-3)then
              ! Add a temperature response equation herefor rho_N_up0 (Zhu Qing 2016)
              totNup = pp%mineralN
     &               * rho_N_up0 * dtsec/seconds_per_hour
     &               * N_roots/(N_roots0+N_roots)
                    !* (1.-exp(-rho_N_up0 * dtsec/seconds_per_hour* N_roots/(N_roots0+N_roots))
!     &                 * exp(9000.0*(1/298.16-1/(273.16+tsoil))) ! kgN m-2 time step-1
              avgNup = totNup / N_roots ! kgN time step-1 kg roots-1

              ! Nitrogen uptaken by each cohort, N_uptake
              cop => pp%tallest
              do while(ASSOCIATED(cop))
                 NSNmax = fnsnmax * cop%C_lab
                 if(cop%N_lab < NSNmax)then
                     cop%N_uptake = cop%C_froot * avgNup  !
                     cop%N_lab    = cop%N_lab   + cop%N_uptake
                     cop%N_up_yr  = cop%N_up_yr + cop%N_uptake

                     ! patch nitrogen update
                     N_up = cop%N_uptake * cop%n
                     pp%N_uptake  = pp%N_uptake  + N_up
                     pp%N_up_day  = pp%N_up_day  + N_up
                     pp%N_up_yr   = pp%N_up_yr   + N_up
                 endif

                 cop => cop%shorter
              enddo
              !! Update soil mineral N
              pp%mineralN = pp%mineralN - pp%N_uptake
           endif ! N_roots>0
        endif ! pp%mineralN > 0.0

      end subroutine BiomeE_N_uptake

! =============================================================================
!     Added by Weng 08-02-2017
      subroutine CrownLAImax_update(pp)
!@sum  used for updating LAImax according to mineral N in soil
!       There is a problem here: all PFTs' crown LAImax are updated!!
      type(patch),   intent(inout) :: pp

      !------local var -----------
      type(cohort), pointer :: cop
        real*8 :: LAImin, LAIfixedN, LAImineralN
        real*8 :: LAI_Nitrogen
        real*8 :: fixedN
        logical:: fixedN_based
        integer :: i
        ! Calculating LAI max based on mineral N or mineralN + fixed N
        fixedN_based = .True. ! .False. !
        LAImin       = 0.25
        if(pp%previousN<0.0)then
            pp%previousN = pp%annualN ! for the first year
        else
            pp%previousN = 0.8 * pp%previousN + 0.2 * pp%annualN
        endif
        do i=1, N_PFT
            associate (sp => BMEspdata(i))

            LAIfixedN   = 0.5 * sp%Nfixrate0 * sp%CNleaf0 * sp%leafLS
            LAImineralN = 0.5 * pp%previousN*sp%CNleaf0*sp%leafLS/sp%LMA
            LAI_nitrogen = LAIfixedN + LAImineralN

            sp%CrownLAImax = MAX(LAImin,
     &              MIN(LAI_nitrogen,sp%cLAI_light))

            end associate
        enddo
      end subroutine CrownLAImax_update

!==========================================================================
! Merge similar cohorts in a patch. Sorting before call this subroutine
      subroutine Merge_cohorts (pp)
      type(patch), intent(inout) :: pp

! ---- local vars
      type(cohort), pointer :: cop  ! current cohort
      type(cohort), pointer :: cop2, csp2 ! the pointers of current cohort
      real*8  :: w1, w2 ! weight for cohorter merge
      integer :: i,j,k,nCohorts
      integer :: c2ID, nMerged

      w1 = 1.0
      w2 = 1.0
      call update_csp_pointers (pp)
      nMerged = 0
      cop   => pp%tallest
      csp2 => cop%cspshorter
      !looping through cohorts to find out the cohorts to be merged
      do while(ASSOCIATED(cop).and.ASSOCIATED(csp2))

         if(MergeableCohorts(csp2,cop))then
             c2ID = csp2%chID
             call Put_c1_into_c2 (csp2,cop)
             ! update derived variables
             call treeBM2structure (cop)

             call Remove_a_cohort (pp,csp2)
             nMerged = nMerged + 1
             !write(6,*)'a cohort merged into =',c2ID,cop%chID
             cop%chID = c2ID
         endif
         cop => cop%shorter
         if(associated(cop))csp2 => cop%cspshorter
      enddo

      nCohorts = 0
      cop   => pp%tallest
      !looping through cohorts to find out the cohorts to be merged
      do while(ASSOCIATED(cop))
         nCohorts = nCohorts + 1
         cop => cop%shorter
      enddo
!      write(6,*)'pp%n, nCC, nMerged',pp%nCohorts,nCohorts,nMerged
      !call sys_flush(6)
      pp%nCohorts = pp%nCohorts - nMerged
      call update_csp_pointers (pp)

      end subroutine Merge_cohorts

!=================================================================================
! Weng, 2016-12-09
      subroutine Remove_a_cohort (pp,cop)
      ! remove 'cop' from pp cohorts
      type(patch),  intent(inout) :: pp ! The patch that hosts the cohorts
      type(cohort), pointer :: cop  ! current cohort

      ! Re-arrange the pointers pointed to this cohort
      if(ASSOCIATED(cop%taller).and.ASSOCIATED(cop%shorter))then ! Not the tallest and not the shortest
          !link the taller one and the shorter on
          cop%taller%shorter => cop%shorter
          cop%shorter%taller => cop%taller
      elseif(ASSOCIATED(cop%shorter))then   ! the tallest cohort.
          pp%tallest => cop%shorter ! update the patch's tallest pointer
          pp%tallest%taller =>  null()

      elseif(ASSOCIATED(cop%taller))then !  the shortest one
          pp%shortest => cop%taller ! Update the patch's shortest pointer
          pp%shortest%shorter => null()

      else ! the only one
          pp%tallest => null()
          pp%shortest => null()
      endif

      ! for the same species pointers
      if(ASSOCIATED(cop%csptaller))
     &    cop%csptaller%cspshorter => cop%cspshorter
      if(ASSOCIATED(cop%cspshorter))
     &    cop%cspshorter%csptaller => cop%csptaller

      cop%taller => null()
      cop%shorter => null()
      cop%csptaller => null()
      cop%cspshorter => null()

      ! Release the memory  ?????? leaking !!!!
      deallocate(cop)
      cop => null()
      !nullify(cop)

      end subroutine Remove_a_cohort

! ============================================================================
      subroutine Remove_ghost_cohorts (pp)
!@sum Remove low density cohorts and low NSC cohorts, Weng, 2021-11-17
      type(patch), intent(inout) :: pp
! ---- local vars
      type(cohort), pointer :: cop  ! current cohort
      type(cohort), pointer :: cop2 ! pointing to the next cohort
      real*8, parameter :: mindensity = 5.0E-5
      real*8 :: loss_alive,loss_wood
      integer :: i,j,k,nCohorts

      !looping through cohorts to find out low density cohorts
      cop => pp%tallest
      nCohorts =0
      do while(ASSOCIATED(cop))
         nCohorts = nCohorts + 1
         cop2  => cop%shorter
         if(pp%nCohorts>1 .and.
     &     (cop%n<mindensity.OR.cop%C_lab<0.2d0*cop%NSCstar))then ! C_lab negative cohort
!              Put the dead trees' carbon to soil pools, pp%Tpool (1,12,1)
               call BiomeE_Plant2Soil(cop,pp,cop%n)

             ! Re-arrange the pointers pointed to this cohort
             call Remove_a_cohort (pp,cop)
             pp%nCohorts = pp%nCohorts - 1
         endif  !(cop%n < mindensity), killing a cohort

         cop => cop2
      enddo

      ! Revive the ghost trees, Weng, 11/17/2021
       cop => pp%tallest
       if(pp%nCohorts==1 .and.cop%C_lab<0.2d0*cop%NSCstar)then
          cop%C_lab=cop%C_lab + 2.d0*cop%NSCstar
          cop%C_hw =Max(cop%C_hw - 2.d0*cop%NSCstar, 0.d0)
       endif
      end subroutine Remove_ghost_cohorts

!=================================================================================
! rank cohorts in descending order by height.
      subroutine Sorting_cohorts (pp)
      type(patch),intent(inout) :: pp

! ---- local vars ---------------------
      type(cohort), pointer :: cop,tmp  ! current cohort
      type(cohort), pointer :: cop1, cop2 ! taller and shorter cohorts
      logical :: skip_section
      integer :: nCohorts,i, n

      skip_section = .false.

      !looping through cohorts to find out how many cohorts are in this patch
      ! and find out the tallest and shortest cohorts
      pp%tallest%taller => null()
      pp%shortest%shorter => null()
      cop1 => pp%tallest
      cop2 => pp%shortest
      nCohorts =0
      cop => pp%tallest
      do while(ASSOCIATED(cop))
         nCohorts = nCohorts + 1
         ! Nullify csp pointers
         cop%cspshorter => null()
         cop%csptaller  => null()
         ! cop1 => tallest, cop2 => shortest
         if(cop1%h .lt. cop%h) cop1 => cop ! point to the taller
         if(cop2%h .gt. cop%h) cop2 => cop ! point to the shorter
         cop  => cop%shorter
      enddo
      pp%nCohorts = nCohorts

      ! update the tallest cohort
      if(cop1%h > pp%tallest%h)then
         ! link cop1's taller and shorter cohorts
         if(associated(cop1%shorter))
     &        cop1%shorter%taller => cop1%taller
         if(associated(cop1%taller))
     &        cop1%taller%shorter => cop1%shorter
         if(.not. associated(cop1%shorter))! cop1 is in the tail
     &        pp%shortest => cop1%taller

         ! put the tallest (cop1) in the head
         pp%tallest%taller => cop1
         cop1%shorter => pp%tallest
         pp%tallest => cop1
         pp%tallest%taller => null()
      endif

      ! update shortest cohort
      if(cop2%h < pp%shortest%h)then
         ! link cop2's taller and shorter cohorts
         if(associated(cop2%shorter))
     &       cop2%shorter%taller=>cop2%taller
         if(associated(cop2%taller))
     &       cop2%taller%shorter=>cop2%shorter
         ! Put 'cop2' in the tail
         pp%shortest%shorter => cop2
         cop2%taller => pp%shortest
         pp%shortest => cop2
         pp%shortest%shorter => null()
      endif
      cop1 => null()
      cop2 => null()

!     sorting those in between
      if(pp%nCohorts .ge. 4)then
        cop => pp%tallest%shorter
        nCohorts=1
        do while(ASSOCIATED(cop%shorter))
          cop1 => cop
          tmp => cop%shorter
          nCohorts = nCohorts + 1
          n=nCohorts
          ! Test dead loop
          do while(ASSOCIATED(tmp%shorter))
             n = n + 1
             if(tmp%h > cop1%h) cop1 => tmp ! cop1 always points to the taller one
             !call sys_flush(6)
             if(n>pp%nCohorts+5) stop
             tmp => tmp%shorter
          enddo ! tmp

          ! move cop1 before cop if cop%h < cop1%h
          if(cop1%h > cop%h)then
              ! link cop1's taller and shorter cohorts
              cop1%shorter%taller=>cop1%taller
              cop1%taller%shorter=>cop1%shorter

              ! Insert cop1 before cop
              cop1%shorter => cop
              cop1%taller  => cop%taller
              cop%taller%shorter => cop1
              cop%taller => cop1

!              !call sys_flush(6)
          else
              cop => cop%shorter
          endif

        enddo ! cop
      endif ! pp%nCohorts >= 4

      ! find out csp taler and shorter cohorts (same PFT)
      call update_csp_pointers (pp)


      end subroutine Sorting_cohorts

!=================================================================================
! Layering according to height and crown area
      subroutine Layering_cohorts (pp)
!----- Weng 12-19-2016 --------------
      use ent_const, only: fgap,tolerance
      type(patch), intent(inout) :: pp ! input patch
    ! ---- local vars
      type(cohort),pointer :: cop, newc

      real*8  :: totalCA ! total crown area
      real*8  :: dn
      real*8  :: CAcrown  ! the area taken by a crown (= (1+fgap)*crownarea
      integer :: L ! layer index (top-down)
      integer :: nCohorts,i

      newc => null()
      nCohorts = 0
      totalCA  = 0.0
      L        = 1
      cop => pp%tallest
      do while(ASSOCIATED(cop))
         nCohorts = nCohorts + 1
         CAcrown = Max((1.d0+fgap)*cop%crownarea,
     &             BMEspdata(cop%pft)%Acrown0)
         totalCA = totalCA + CAcrown * cop%n
         cop%layer = L
         if(totalCA .gt. L )then
            dn = (totalCA - L) / CAcrown ! plants that are put to the next layer
            L = L + 1 ! update ayer number
            ! Split the last cohort in this layer:Generate a "newc" to stay
            ! in the upper layer and "cop" stays in the lower layer
            ! creat a new cohort that is the same with the current one
            call cohort_construct(newc,pp,cop%pft)
            call copy_cohort(cop, newc)
            !! in case the cohort is the tallest
            if(nCohorts ==1) pp%tallest=>newc

            !! update new cohort
            newc%n     = cop%n - dn
            newc%layer = L - 1
            newc%shorter    => cop
            newc%cspshorter => cop
            if(ASSOCIATED(cop%taller))then
                cop%taller%shorter => newc
                newc%taller        => cop%taller
            endif
            if(ASSOCIATED(cop%csptaller))then
                cop%csptaller%cspshorter => newc
                newc%csptaller    => cop%csptaller
            endif

            !! Update the original cohort, which has been put in the next layer
            cop%n  = cop%n - newc%n
            cop%taller    => newc
            cop%csptaller => newc ! same sp
            cop%layer = L

            ! Nullify newc and Count the new cohort
            newc => null()
            nCohorts = nCohorts + 1
         endif !'totalCA .gt. L'
         ! go to the next cohort
         cop  => cop%shorter
      enddo
      pp%nCohorts = nCohorts
      pp%totalCA  = totalCA
      pp%crownlayers = L

      !call update_csp_pointers (pp)
      end subroutine Layering_cohorts

!=================================================================================
! rank cohorts in descending order by height.
! Layering is included in this subroutine
      subroutine update_csp_pointers (pp)
      type(patch),intent(inout) :: pp

! ---- local vars ---------------------
      type(cohort), pointer :: cop  ! current cohort
      type(cohort), pointer :: cop1 ! taller and shorter cohorts
      integer :: nCohorts,i, n

      ! nullify tallest and shortest
      pp%tallest%taller => null()
      pp%shortest%shorter => null()
      nCohorts =0
      cop => pp%tallest
      do while(ASSOCIATED(cop))
         nCohorts = nCohorts + 1
         ! Nullify csp pointers
         cop%cspshorter => null()
         cop%csptaller  => null()

         cop  => cop%shorter
      enddo
      pp%nCohorts = nCohorts

      ! find out csp taler and shorter cohorts (same PFT)
      cop => pp%tallest
      do while(ASSOCIATED(cop%shorter))
         cop1 => cop%shorter
         do while(ASSOCIATED(cop1))
             if(cop%pft.eq.cop1%pft) then
                cop%cspshorter => cop1
                cop1%csptaller => cop
                exit ! exit "do while(ASSOCIATED(cop1))"
             else
                cop1 =>cop1%shorter
             endif
         enddo
         cop  => cop%shorter
      enddo

      end subroutine update_csp_pointers

!=================================================================================
!------------------------------------
!----- Weng 03-31-2017 --------------
! Calculate Patch LAI
      subroutine get_patch_LAI (pp)
      type(patch), intent(inout) :: pp ! input patch
    ! type(patch),pointer :: pp

      ! ---- local vars
      type(cohort),pointer :: cop
      integer :: nCohorts,i

      !-----------------
      pp%LAI = 0.0
      nCohorts = 0
      cop => pp%tallest
      do while(ASSOCIATED(cop))
         nCohorts = nCohorts + 1
         cop%LAI = cop%leafarea * cop%n
         pp%LAI = pp%LAI + cop%leafarea * cop%n

         ! go to the next cohort
         cop  => cop%shorter
      enddo
      pp%nCohorts = nCohorts

      ! Update max LAI for a year, just for reporting
      pp%LAIyr = Max(pp%LAIyr,pp%LAI)

      end subroutine get_patch_LAI

!==============================================================================
      subroutine BiomeE_diagnostics(pp)
!@sum output daily Ent-PPA diagnostics

        implicit none
        type(patch),   intent(inout) :: pp
        !------local var -----------
        type(cohort), pointer :: cop
        real*8  :: VegC,VegN,runyrs
        integer :: nCohorts, n, i
        integer :: iu896,iu898,iu897,iu899

          ! Summarize patch level fluxes and pools
          pp%C_lab = 0.
          pp%C_fol = 0.
          pp%C_w   = 0.
          pp%C_froot = 0.
          pp%C_seed = 0.
          pp%N_lab =0.
          pp%N_fol = 0.
          pp%N_w = 0.
          pp%N_froot = 0.
          pp%N_seed = 0.
          ! fluxes
          pp%GPP_daily = 0.0
          pp%NPP_daily = 0.0
          pp%Ra_daily = 0.0

          runyrs = pp%age/(365.0*24.0*3600.0)
          ! calculating diagnostics
          cop => pp%tallest
          do while(ASSOCIATED(cop))
            ! Cohort yearly
            cop%GPP_yearly = cop%GPP_yearly + cop%GPP_daily
            cop%NPP_yearly = cop%NPP_yearly + cop%NPP_daily
            cop%Ra_yearly  = cop%Ra_yearly  + cop%Ra_daily
            ! Patch daily
            pp%GPP_daily = pp%GPP_daily + cop%GPP_daily * cop%n
            pp%NPP_daily = pp%NPP_daily + cop%NPP_daily * cop%n
            pp%Ra_daily = pp%Ra_daily   + cop%Ra_daily * cop%n
            ! Patch plant pools
            pp%C_lab =pp%C_lab + cop%C_lab * cop%n
            pp%C_fol =pp%C_fol + cop%C_fol * cop%n
            pp%C_w   =pp%C_w   + (cop%C_sw + cop%C_hw) * cop%n
            pp%C_froot =pp%C_froot + cop%C_froot * cop%n
            pp%C_seed =pp%C_seed + cop%C_seed * cop%n
            pp%N_lab =pp%N_lab   + cop%N_lab * cop%n
            pp%N_fol =pp%N_fol   + cop%N_fol * cop%n
            pp%N_w =pp%N_w + (cop%N_sw + cop%N_hw) * cop%n
            pp%N_froot = pp%N_froot + cop%N_froot * cop%n
            pp%N_seed = pp%N_seed + cop%N_seed * cop%n

            cop  => cop%shorter
          enddo
          pp%GPP_yearly = pp%GPP_yearly + pp%GPP_daily
          pp%NPP_yearly = pp%NPP_yearly + pp%NPP_daily
          pp%Ra_yearly  = pp%Ra_yearly  + pp%Ra_daily
          pp%Rh_yearly  = pp%Rh_yearly  + pp%Rh_daily

          iu896 = 896
          iu898 = 898
#ifndef SITE_TEST
          !Harvard Forest: 44067; tapajos: 51043
          !grids(8)=(/14077,33073,44067,39063,32065,30063,28061,51043/)
          if(ANY(grids==ijdebug).and.runyrs>=150 .and.runyrs<181.)then
             do i=1, N_grids
                if(ijdebug == grids(i))then
                   iu896 = 896*10 + i
                   iu898 = 898*10 + i
                   exit
                endif
             enddo
#endif
            ! output cohorts
            write(iu898,'(4(I5,","))')pp%nCohorts,
     &        pp%cellptr%year1,pp%cellptr%doy,pp%doy

            cop => pp%tallest
            do while(ASSOCIATED(cop))
              write(iu898,'(4(I7,","),L4,",",I8,",",40(E12.4,","))')
     &           ijdebug,pp%doy,cop%pft,cop%layer,
     &           cop%PhenoON,cop%phenostatus,cop%phenofactor,cop%n,
     &           cop%dbh,cop%h,cop%crownarea,cop%leafarea,cop%LAI,
     &           cop%GPP_daily*1000,cop%Ra_daily*1000,
     &           cop%C_lab,cop%C_fol,
     &           cop%C_sw,cop%C_hw,cop%C_froot,cop%C_seed,
     &           cop%N_up_yr*1000, cop%N_lab*1000,cop%N_fol*1000,
     &           cop%N_sw*1000,cop%N_hw*1000,
     &           cop%N_froot*1000,cop%N_seed*1000
              cop  => cop%shorter
            enddo

            write(iu896,'(1(F12.4,","),2(I7,","),96(E12.4,","))')
     &          pp%age/(365.0*24.0*3600.0),ijdebug,pp%doy,
     &          pp%GPP_daily,pp%NPP_daily,pp%Ra_daily,pp%Rh_daily,
     &          pp%N_up_day, pp%N_min_day,
     &          pp%C_lab,pp%C_fol,
     &          pp%C_w,pp%C_froot,pp%C_seed,
     &          (pp%SOC(i),i=1,N_SOM),
     &          (pp%SON(i)*1.0d3,i=1,N_SOM),
     &          pp%mineralN,
     &          pp%LAI

#ifndef SITE_TEST
          endif ! grid number
#endif
          ! ----- Output yearly diagnostics ----------------------
          !if(ecop%gdd .eq. 0.0 .and. pp%tallest%C_seed > 0.000001) then
          if(pp%doy .eq. 365)then ! Print the tallest cohort again, for test

            pp%annualN = pp%N_min_yr + N_dpst ! for next year's LAI

            iu897 = 897
            iu899 = 899
#ifndef SITE_TEST
            !Harvard Forest: 44067; tapajos: 51043
            !grids(8)=(/14077,33073,44067,39063,32065,30063,28061,51043/)
            if(ANY(grids == ijdebug)) then
              do i=1, N_grids
                if(ijdebug == grids(i))then
                   iu897 = 897*10 + i
                   iu899 = 899*10 + i
                   exit
                endif
              enddo
#endif
              ! output the states of cohorts to a file
              write(iu899,'(1(I4,","))')pp%nCohorts
              cop => pp%tallest
              do while(ASSOCIATED(cop))
                 write(iu899,'(5(I7,","),50(F15.6,","))')ijdebug,
     &             pp%doy,cop%pft,cop%layer,cop%phenostatus,
     &             cop%n,cop%dbh,cop%h,cop%crownarea,
     &             cop%leafarea,cop%LAI,
     &             cop%GPP_yearly,cop%Ra_yearly,
     &             cop%C_lab,cop%C_fol,
     &             cop%C_sw,cop%C_hw,cop%C_froot,cop%C_seed,
     &             cop%N_up_yr*1000, cop%N_lab*1000,cop%N_fol*1000,
     &             cop%N_sw*1000,cop%N_hw*1000,
     &             cop%N_froot*1000,cop%N_seed*1000,cop%mu
                 cop  => cop%shorter
              enddo

              VegC=pp%C_lab+pp%C_fol+pp%C_w+pp%C_froot+pp%C_seed
              VegN=pp%N_lab+pp%N_fol+pp%N_w+pp%N_froot+pp%N_seed
              write(iu897,'((I7,","),7(F15.6,","),40(F15.6,","))')
     &          ijdebug,runyrs, pp%GPP_yearly,pp%NPP_yearly,
     &          pp%Ra_yearly,pp%Rh_yearly,pp%N_up_yr,pp%N_min_yr,
     &          pp%C_lab,pp%C_fol,pp%C_w,pp%C_froot,pp%C_seed,VegC,
     &          (pp%SOC(i),i=1,N_SOM),sum(pp%SOC),
     &          pp%N_lab,pp%N_fol,pp%N_w,pp%N_froot,pp%N_seed,VegN,
     &          (pp%SON(i),i=1,N_SOM),sum(pp%SON),
     &          pp%mineralN,pp%LAIyr

#ifndef SITE_TEST
              endif ! grid number
#endif
          endif ! output

          ! zeroing
          if(pp%doy .eq. 365)then
              ! zero yearly diagnostics
              cop => pp%tallest
              do while(ASSOCIATED(cop))
                 cop%GPP_yearly = 0.0
                 cop%NPP_yearly = 0.0
                 cop%Ra_yearly  = 0.0
                 cop%N_up_yr = 0.0
                 cop  => cop%shorter
              enddo
              pp%LAIyr = 0.d0
              pp%GPP_yearly = 0.0
              pp%NPP_yearly = 0.0
              pp%Ra_yearly  = 0.0
              pp%Rh_yearly = 0.0
              pp%N_min_yr = 0.0
              pp%N_up_yr  = 0.0
          endif ! doy=365

          ! zero daily variables, Weng, 2017-02-07
          cop => pp%tallest
          do while(ASSOCIATED(cop))
                cop%GPP_daily = 0.0
                cop%NPP_daily = 0.0
                cop%Ra_daily  = 0.0
                cop  => cop%shorter
          enddo
          ! Zero patch fluxes, daily
          pp%GPP_daily = 0.0
          pp%NPP_daily = 0.0
          pp%Ra_daily  = 0.0
          pp%Rh_daily  = 0.0
          pp%N_min_day = 0.0
          pp%N_up_day  = 0.0

      end subroutine BiomeE_diagnostics

!==============================================================================
      subroutine write_outputfile_headers()

      integer :: iu896,iu898,iu897,iu899,i

      iu896 = 896
      iu898 = 898
      iu897 = 897
      iu899 = 899

      !grids(8)=(/14077,33073,44067,39063,32065,30063,28061,51043/)
#ifndef SITE_TEST
      do i=1, N_grids
          iu896 = 896*10 + i
          iu898 = 898*10 + i
          iu897 = 897*10 + i
          iu899 = 899*10 + i
          write(*,*)'output file name',iu897,iu899
#endif
          write(iu897,'(60(a9,","))') ! Patch, Yearly
     &          'ij','age','GPP','NPP','Rauto','Rh','N_up','N_min', ! fluxes
                ! Plant and soil carbon pools
     &          'C_lab','C_fol','C_w','C_froot','C_seed','VegC',
     &          'fineL','struL','Micr','fast','slow','SOC',
                ! Plant and soil nitrogen pools
     &          'N_lab', 'N_fol','N_w','N_froot','N_seed','VegN',
     &          'Nfine','Nstru','Nmic','Nfast','Nslow','SON', ! Soil Nitrogen
     &          'mineralN', 'LAI','BiomeE'! ,'pheno' ! Phenology

          write(iu899,'(4(a5,","),60(a9,","))') ! Cohort, Yearly
     &          'ij','doy','PFT','layer','pheno',
     &          'density','dbh','height','Acrown','Aleaf','LAI',
     &          'GPP','Rauto',
     &          'C_lab', 'C_fol','C_sw','C_hw','C_froot','C_seed',
     &          'Nup','N_lab', 'N_fol','N_sw','N_hw','N_froot','N_seed',
     &          'mu'
         write(iu896,'(80(a9,","))') ! Patch, Daily
     &        'age','ij','doy',
     &        'GPP','NPP','Rauto','Rh', ! Carbon fluxes
     &        'N_up','N_min',                 ! Nitrogen fluxes
     &        'C_lab','C_fol',
     &        'C_w','C_froot','C_seed',
     &        'fineL','struL','Micr','fast','slow',
     &        'NfineL','NstruL','Nmicr','Nfast','Nslow',
     &        'mineralN','LAI' ! Phenology

         write(iu898,'(4(a6,","),40(a9,","))') ! Cohort, Daily
     &     'ij','doy','PFT','layer','PhenoON','status','factor',
     &     'density','dbh','height','Acrown','Aleaf','LAI',
     &     'GPP','Rauto',
     &     'C_lab', 'C_fol','C_sw','C_hw','C_froot','C_seed',
     &     'Nup','N_lab', 'N_fol','N_sw','N_hw','N_froot','N_seed'

#ifndef SITE_TEST
      enddo
#endif

!     ! Header for the cohort output
      write(798,'(4(a5,","),30(a10,","))') ! cohort (individual), hourly or half-hourly,
     &        'ch','doy','PFT','layer',
     &        'runyrs','density','Acrown',
     &        'GPP','NPP','Rauto',
     &        'LAI','pheno','LSM_PPA'

      write(799,'(1(a5,","),50(a10,","))') ! patch level, hourly
     &        'doy','age',
     &        'GPP','NPP','Rauto','Rh',
     &        'LAI','C_L','C_W','C_FR','Trans_SW',
     &        'Ci','Gcanopy','IPP','Ctot','C_growth',
     &        'Nuptake','betad'



      end subroutine write_outputfile_headers

!============================================================================
! Below are the allometry equations from Weng et al. 2015, Biogeosciences.
! Added by Weng 11-17-2016
! dbh2structuralC ,dbh2Height , dbh2Crownarea ,dbh2CrownRadius ,
! structuralC2dbh , Height2dbh , dCstrctural2dDBH ,dDBH2dHeight ,
! dDBH2dCrownarea
!========================================================================
!===========Allometry functions============
! calculate tree height, DBH, height, and crown area by bwood and denstiy
! The allometry equations are from Forrior et al. 2013
!   HT = alphaHT * DBH ** (gamma-1)   ! DBH --> Height
!   CA = alphaCA * DBH ** gamma       ! DBH --> Crown Area
!   BM = alphaBM * DBH ** (gamma + 1) ! DBH --> tree biomass
!=========================================================================

      subroutine treeBM2structure(cop)
      type(cohort), intent(inout) :: cop

      !-------local-------
      real*8    :: structuralC ! total structural C per individual, kgC/tree

      associate(sp =>BMEspdata(cop%pft))
      structuralC = cop%C_sw + cop%C_hw
      cop%dbh     = (structuralC/sp%alphaBM)**(1.d0/sp%thetaBM)
      cop%h       = sp%alphaHT * cop%dbh ** sp%thetaHT
      cop%crownarea = sp%alphaCA * cop%dbh ** sp%thetaCA
      cop%crown_dx = SQRT(1.d0/PI*sp%alphaCA*cop%dbh**sp%thetaCA)
      cop%crown_dy = crown_radius_vert(cop%pft,cop%h,cop%crown_dx)

      cop%leafarea = cop%C_fol/sp%LMA
      cop%CrownLAI = cop%leafarea/cop%crownarea
      cop%LAI      = cop%leafarea * cop%n ! Weng: I don't agree with this definition !!!
      cop%C_croot  = 0.25 * structuralC

      end associate
      end subroutine treeBM2structure

!**************************************************************************

      real*8 function dbh2structuralC(pft,dbh)result(structuralC)
      ! return structuralC
!@sum (gC/plant) Carbon in structural tissues as a function of dbh
!@+   Structural tissue = Sapwood + Heart wood + coarse roots

      integer,intent(in) :: pft
      real*8, intent(in) :: dbh !(m)

      associate ( sp => BMEspdata(pft) )
      structuralC = sp%alphaBM * dbh ** sp%thetaBM  !kgC/tree
      end associate

      end function  dbh2structuralC

!**************************************************************************
      real*8 function dbh2Height(pft,dbh) result(Height)
!@sum Tree height as a function of dbh

      integer,intent(in) :: pft
      real*8, intent(in) :: dbh !(m)

      associate ( sp => BMEspdata(pft) )
      Height = sp%alphaHT * dbh ** sp%thetaHT  !m
      end associate

      end function  dbh2Height
!**************************************************************************
      real*8 function dbh2Crownarea(pft,dbh)result (Crownarea)
!@sum Crown area as a function of DBH

      integer,intent(in) :: pft
      real*8, intent(in) :: dbh !(m)

      associate ( sp => BMEspdata(pft) )
      Crownarea = sp%alphaCA * dbh ** sp%thetaCA  ! m^2
      end associate

      end function  dbh2Crownarea

!**************************************************************************
      real*8 function dbh2CrownRadius(pft,dbh)result(CrownRadius)

!@sum Crown diameter (m) as a function of dbh, x
      integer,intent(in) :: pft
      real*8, intent(in) :: dbh !(m)

      associate ( sp => BMEspdata(pft) )
      CrownRadius = SQRT(1.d0/PI*sp%alphaCA*dbh**sp%thetaCA)  ! m
      end associate

      end function  dbh2CrownRadius

!**************************************************************************
      real*8 function crown_radius_vert(pft, h,crx)
!@sum Vertical crown radius (m) from allometry.
!@+   Subject to change.  Currently allows tall ellipsoid to spherical but not
!@+   oblate crowns. May want to allow oblate for understory crowns.
      integer :: pft
      real*8 :: h, crx !Tree height, crown horizontal radius

         ! Set PFt-specific parameters
         crown_radius_vert = min(2.7*crx, 0.49d0*h)

      end function crown_radius_vert

!**************************************************************************
      real*8 function structuralC2dbh(pft,structuralC) result(dbh)
!@sum (gC/plant) DBH as a function of Carbon in structural tissues
!@+   Structural tissue = Sapwood + Heart wood + coarse roots

      integer,intent(in) :: pft
      real*8, intent(in) :: structuralC !(kg/tree)

      associate ( sp => BMEspdata(pft) )
      dbh = (structuralC/sp%alphaBM) ** (1.d0/sp%thetaBM) ! m
      end associate

      end function  structuralC2dbh

!**************************************************************************
      real*8 function Height2dbh(pft,height) result(dbh)
!@sum DBH as a function of tree height

      integer,intent(in) :: pft
      real*8, intent(in) :: height !(m)

      associate ( sp => BMEspdata(pft) )
      dbh = (1.d0/sp%alphaHT*height)**(1.d0/sp%thetaHT) ! m
      end associate

      end function  Height2dbh

!**************************************************************************
      real*8 function height2LAmax(pft,height) result (LAmax)
!@sum max leaf area as a function of height, only for top layer trees

      integer,intent(in) :: pft
      real*8, intent(in) :: height !(m)
      !-------local-------
      real*8 :: dbh

      associate ( sp => BMEspdata(pft) )
      dbh = (height/sp%alphaHT)**(1.d0/sp%thetaHT)
      LAmax = sp%CrownLAImax * sp%alphaCA * dbh ** sp%thetaCA  ! m^2
      end associate

      end function  height2LAmax

!**************************************************************************
      real*8 function dCstrctural2dDBH(cop,dCstructural)
     &               result(dDBH)
!@sum DBH growth when the Structural C pool is added carbon 'dCstrctural'

      type(cohort),intent(inout) :: cop
      real*8, intent(in) :: dCstructural !(kg C), new carbon allocated to
                    !structural tissues (active part of SW and CR)

      associate ( sp => BMEspdata(cop%pft) )
      dDBH = dCstructural /
     &      (sp%thetaBM * sp%alphaBM * cop%DBH**(sp%thetaBM-1.d0))
      end associate

      end function dCstrctural2dDBH

!**************************************************************************
      real*8 function dDBH2dHeight(cop,dDBH) result(dHeight)
!@sum Height growth given DBH growth (dDBH)

      type(cohort),intent(inout) :: cop
      real*8, intent(in) :: dDBH

      associate ( sp => BMEspdata(cop%pft) )
      dHeight = sp%thetaHT * sp%alphaHT *
     &      cop%DBH**(sp%thetaHT-1.d0) * dDBH
      end associate

      end function dDBH2dHeight

!**************************************************************************
      real*8 function dDBH2dCrownarea(cop,dDBH) result(dCrownarea)
!@sum Crown area growth given DBH growth (dDBH)

      type(cohort),intent(inout) :: cop
      real*8, intent(in) :: dDBH

      associate ( sp => BMEspdata(cop%pft) )
      dCrownarea  = sp%thetaCA * sp%alphaCA *
     &      cop%DBH**(sp%thetaCA-1.d0) * dDBH
      end associate

      end function dDBH2dCrownarea

!**************************************************************************
!     Added by Weng, 05-30-2017
      subroutine rootprofile(cop)
!@sum Return fractions of roots in soil layer
!@+   for a single plant. (Rosenzweig & Abrampoulos, 1997)
!@+   ! Individual level, pft-specific.
!     Weng et al. 2015

      use ent_pfts, only: COVEROFFSET, aroot, broot
      real*8 :: rootprof(N_DEPTH)
      integer :: ncov !plant functional type + COVEROFFSET
      type(cohort),intent(inout) :: cop

      !-----Local variables------------------
      real*8,parameter :: dz_soil(1:6)=  !N_DEPTH
     &     (/  0.99999964d-01,  0.17254400d+00,
     &     0.29771447d+00,  0.51368874d+00,  0.88633960d+00,
     &     0.15293264d+01 /)
      integer :: n,l
      real*8 :: z, frup,frdn

      call prescr_calc_rootprof(cop%fracroot,cop%pft + COVEROFFSET)

      end subroutine rootprofile

!===================================================================
      subroutine prescr_calc_rootprof(rootprof, ncov)
!@sum Return prescribed array rootprof of fractions of roots in soil layer
!@+   for single cover type. (Rosenzweig & Abrampoulos, 1997)
      use ent_pfts, only: COVEROFFSET, aroot, broot
      real*8 :: rootprof(:)
      integer :: ncov !plant functional type + COVEROFFSET
      !-----Local variables------------------
      real*8,parameter :: dz_soil(1:6)=  !N_DEPTH
     &     (/  0.99999964d-01,  0.17254400d+00,
     &     0.29771447d+00,  0.51368874d+00,  0.88633960d+00,
     &     0.15293264d+01 /)
      integer :: n,l
      real*8 :: z, frup,frdn

c**** calculate root fraction afr averaged over vegetation types
      !Initialize zero
      do l=1,N_DEPTH
        rootprof(l) = 0.0
      end do
      do n=1,N_DEPTH
        if (dz_soil(n) <= 0.0) exit !Get last layer w/roots in it.
      end do
      n=n-1
      z=0.
      frup=0.
      do l=1,n
        z=z+dz_soil(l)
        frdn=aroot(ncov)*z**broot(ncov) !cumulative root distrib.
        !frdn=min(frdn,one)
        frdn=min(frdn,1d0)
        if(l.eq.n)frdn=1.
        rootprof(l) = frdn-frup
        frup=frdn
      end do
      !Return rootprof(:)
      end subroutine prescr_calc_rootprof

!***********************************************************************
      subroutine Soillayer_convert_Ent(s, depthm, Savglayer)
!@sum Calculates layer-weighted average of soil variables to produce
!@+   average values for 0-30 cm and, optional, 30-100 cm.
      implicit none
      real*8 :: s(N_DEPTH)            !Any intensive (non-extensive) variable by LSM soil layer.
      real*8,intent(in) :: depthm(N_DEPTH) !(m) Depths of bottoms of LSM soil layers, increasing.
      real*8, intent(out) :: Savglayer(N_CASA_LAYERS)
      !---Local-----
#ifdef NCASA2
      real*8, parameter :: ENTD(N_CASA_LAYERS) = ( / 0.30d0, 1.0d0 /)
#else
      real*8, parameter :: ENTD(N_CASA_LAYERS) = ( / 0.30d0 /)
#endif
      integer :: k,ek
      real*8 :: sk,ddk,depthup
      real*8 :: savg, dsum
      logical :: entdpassed !Flag for if Ent soil layer depth was passed.

      k=1
      do ek = 1,N_CASA_LAYERS
         entdpassed = .false.
         savg = 0.d0
         dsum = 0.d0
         ddk = depthm(k)
         do while ((.not.entdpassed).and.(k.le.N_DEPTH))
            if (k.eq.1) then !Top layer
               depthup = 0.d0
            else !Lower layers
               depthup = depthm(k-1)
            endif
            if (depthm(k).le.ENTD(ek)) then
               ddk = depthm(k)- depthup
            else
                ddk = ENTD(ek) - depthup
                entdpassed = .true.
            endif
            savg = savg + ddk*s(k)
            dsum = dsum + ddk
            k = k+1
         end do
         Savglayer(ek) = savg/dsum
      end do

      end subroutine Soillayer_convert_Ent

! ============================================================================
! The combined reduction in decomposition rate as a funciton of TEMP and MOIST
! Based on CENTURY Parton et al 1993 GBC 7(4):785-809 and Bolker's copy of
! CENTURY code
      function A_function(soilt, theta) result(A)
        real*8 :: A                 ! return value, resulting reduction in decomposition rate
        real*8, intent(in) :: soilt ! deg C
        real*8, intent(in) :: theta

        real*8 :: soil_temp ! temperature of the soil, deg C
        real*8 :: Td        ! rate multiplier due to temp
        real*8 :: Wd        ! rate reduction due to mositure

        ! coefficeints and terms used in temperaturex term
        real*8 :: Topt,Tmax,t1,t2,tshl,tshr

        soil_temp = soilt !- 273.16

        ! EFFECT OF TEMPERATURE , ! from Bolker's century code
        Tmax=45.0;
        if (soil_temp > Tmax) soil_temp = Tmax;
        Topt=35.0;
        tshr=0.2; tshl=2.63;
        t1=(Tmax-soil_temp)/(Tmax-Topt);
        t2=exp((tshr/tshl)*(1.-t1**tshl));
        Td=t1**tshr*t2;

        if (soil_temp > -10) Td=Td+0.05;
        if (Td > 1.) Td=1.;

        ! EFFECT OF MOISTURE
        ! Linn and Doran, 1984, Soil Sci. Amer. J. 48:1267-1272
        ! This differs from the Century Wd
        ! was modified by slm/ens based on the figures from the above paper
        !     (not the reported function)

        if(theta <= 0.3) then
           Wd = 0.2;
        elseif(theta <= 0.6) then
           Wd = 0.2+0.8*(theta-0.3)/0.3
        else
           Wd = 1.0 ! exp(2.3*(0.6-theta)); ! Weng, 2016-11-26
        endif

        A = (Td*Wd); ! the combined (multiplicative) effect of temp and water
               ! on decomposition rates
      end function A_function
! ============================================================================
      function btotal(cop)
        real*8 :: btotal ! returned value
        type(cohort), intent(in) :: cop

        btotal = cop%C_lab+cop%C_fol+ cop%C_froot
     &         + cop%C_sw + cop%C_hw+ cop%C_seed
      end function
! ============================================================================
      function Ntotal(cop)
        real*8 :: Ntotal ! returned value
        type(cohort), intent(in) :: cop

        Ntotal =  cop%N_lab+cop%N_sw+cop%N_hw
     &          + cop%N_froot+cop%N_seed
      end function

!=================================================================================
! ===== Get PFTs based on climate conditions =================
! ================= Place holder only ====================
      subroutine BiomeE_biogeography(pp,ptenv)
        type(patch),   intent(inout) :: pp
        type(climate_env),intent(inout) :: ptenv
        !------local var -----------
        type(cohort), pointer :: cop
        real*8  :: totalC

        cop => pp%tallest
        do while(ASSOCIATED(cop))
            !totalC = btotal(cop)
            !call update_species(ptenv,totalC)
            !cop%pft = ptenv%spp

            cop => cop%shorter
        enddo

      end subroutine BiomeE_biogeography
!=================================================================================

      end module BiomeE_VEG
