#include "rundeck_opts.h"

      module ent_types
      use ent_const

      implicit none

! the following two parameters are basically needed for integrity
! checks (to avoid infinite loops when processing linked lists)
!@var MAX_PATCHES maximal number of patches per cell
!@var MAX_COHORTS maximal number of cohorts per patch
        integer, parameter :: MAX_PATCHES=32, MAX_COHORTS=64


!****************************************************************************
!*       TYPE DECLARATIONS
!****************************************************************************
      type timestruct
         !Structure for holding time variables
         integer :: year
         integer :: month
         integer :: jday  !day of year
         integer :: day   !day of month, not of year
         integer :: hour  !0-23
         integer :: minute
         integer :: seconds
      end type timestruct
!=========================================================================
! Weng July 30 2019
! for PFT climate envelopes
      type climate_env
        real*8  :: t_ann   ! annual-mean temperature, degK
        real*8  :: t_cold  ! average temperature of the coldest month, degK
        real*8  :: p_ann   ! annual-mean precipitation, mm/yr
        real*8  :: cm      ! number of cold months
        real*8  :: ncd     ! number of chilling days
        integer :: landuse ! pasture, crops, natural vegetation,secondary
        integer :: lifeform  ! 0 for grasses, 1 for trees
        integer :: pstype    ! photosynthetic pathway: C3=1, C4=2
        integer :: phenotype ! 0 for deciduous, 1 for evergreen
        integer :: spp       ! 4 species only
      end type
!=========================================================================
!=====================================================
      type PlantSP
         integer :: pst = 1 ! Photosynth type 1=C3, 2=C4
         logical :: woody = .True. !Woody, FALSE=NO, TRUE=YES.
         integer :: phenotype = 1 ! evergreen (1),cold deciduous (2),
                             ! drought deciduous (3),
                             ! cold/drought deciduous (4),annual (5)

         ! leaftype and plantform are not used in BiomeE
         integer :: leaftype = 1 !1=broadleaf, 2=needleleaf, 3=monocot
         !integer :: plantform = 4 ! GRASS=1; HERB=2; SHRUB=3;TRE=4; BARE=5

         ! traits of leaf, fine root, and wood
         real*8  :: LMA = 0.02    ! leaf mass per area (kgC m-2 leaf)
         real*8  :: leafLS = 1.0 ! leaf lifespan (year)
         real*8  :: LNA    ! Nitrogen per unit leaf area, kg N/m2
         real*8  :: LNbase = 0.6d-3 ! basal leaf Nitrogen per unit area, kg N/m2, (Rubisco)
         real*8  :: Nfixrate0 ! N fixation rate (kgN yr-1 per leaf area ? unit?)
         real*8  :: CNstrucL = 175.d0 ! CN ratio of leaf structural tissues,175
         real*8  :: CrownLAImax = 4.5  ! maximum leaf layers a crown can have
         real*8  :: H0_cL       = 6.d0 ! cLmax = cLAImax* H/(H0_cL + H)
         real*8  :: cLAI_light  = 4.5  ! Light limited crown LAI
         real*8  :: rhoWood = 350.    ! wood density (kgC m-3)
         real*8  :: rhoFR   = 150.    ! fine root carbon density (kgC m-3)
         real*8  :: rootLS  = 0.8     ! fine root lifespan (year)
         real*8  :: root_r  = 2.9d-4  ! radius of the fine roots, m
         real*8  :: root_zeta = 0.29  ! vertical profile parameter: (1/zeta * exp(-z/zeta)
         real*8  :: SRA         ! specific root area (m2 kgC-1 root)
         real*8  :: root_perm = 1.d-5  ! fine root membrane permeability , kg/(m3 s)
         real*8  :: f_taper = 0.7d0 ! Tree structure taper factor
         real*8  :: K_cambium = 0.1d0 ! Cambium respiration rate (kgC m-2 Cambium yr-1)
         real*8  :: f_Acam = 1.7d0 ! Cambium area taper factor
         real*8  :: Acrown0 = 0.04d0 ! m2, minimum plant space
         real*8  :: tauNSC  = 2.5 ! year, NSC use strategy
         real*8  :: phiRL = 0.8      ! ratio of fine root to leaf area
         real*8  :: phiCSA= 0.5d-4    ! ratio of sapwood CSA to target leaf area

         ! Allometry
         real*8  :: alphaHT = 20.d0
         real*8  :: thetaHT = 0.5d0  ! height = alphaHT * DBH ** thetaHT
         real*8  :: alphaCA = 150.d0 !
         real*8  :: thetaCA = 1.5 ! crown area = alphaCA * DBH ** thetaCA
         real*8  :: alphaBM = 8150.0
         real*8  :: thetaBM = 2.5 ! biomass = alphaBM * DBH ** thetaBM

         real*8  :: r_seed    = 0.1   ! the ratio of Cgain to seeds each stem
         real*8  :: seedlingHeight   ! Height of seedling (m), Weng 2016-11-16
         real*8  :: seedlingLA       ! Seedling leaf area (m2/seedling)
         real*8  :: seedlingC0 = 0.5 !  Seedling carbon (kgC/seedling) ! 0.5 for testing resubmission

         real*8  :: prob_g = 1.0     ! germination probability
         real*8  :: prob_e = 1.0     ! establishment probability
         real*8  :: mu0    = 0.02    ! annual mortality rate in canopy
         real*8  :: D0mu   = 0.8     ! m, mortality curve as a function of DBH
         real*8  :: A_SD   = 6.d0    ! Max multifactor of mortality rate for seedlings
         real*8  :: B_SD   = -60.d0  ! Decreasing rate relative to DBH for seedlings
!        for Nitrogen model
         real*8  :: CNleaf0 = 40.0
         real*8  :: CNroot0 = 30.0
         real*8  :: CNsw0   = 150.0
         real*8  :: CNwood0 = 150.0
         real*8  :: CNseed0 = 15.0

         ! Critical temperature for initiating leaf fall
         real*8  :: Tc_pheno0 = 12.d0 ! a function of LMA
         real*8  :: gdd_par1 = 30. ! gdd_threshold = gdd_par1 + gdd_par2*exp(gdd_par3*cop%ncd)
         real*8  :: gdd_par2 = 800.
         real*8  :: gdd_par3 = -0.02
         ! Plant hydraulics (for phenology, added by Weng, 02/20/2021)
         real*8  :: betad_ON  = 0.3 ! 0.6 ! drought-limited leaf onset
         real*8  :: betad_OFF = 0.1 ! 0.3 ! drought-limited leaf offset

         ! Copied here from pftype (Weng, 02/21/2021)

         ! from pspartype, Weng, 02/08/2021
         real*8 :: PARabsorb     !Leaf PAR absorptance (fraction)
         ! Farquhar/Ball-Berry parameters
         real*8 :: Vcmax         !Maximum photosynthetic capacity (umol m-2 s-1)
         real*8 :: m             !Slope of Ball-Berry equation
         real*8 :: b             !Intercept of Ball-Berry equation (mol m-2 s-1)
         real*8 :: Nleaf         !g-N/m2[leaf]

         !Soil water limits
         real*8 :: hwilt  !Wilting point matric potential (m)
         real*8 :: sstar  !Rel. soil moist at stress onset (Rodriguez-Iturbe)
         real*8 :: swilt  !Normalized soil water at wilting point (dim'less)

         !!!!---------- never used in BiomeE ---------!!!!
         real*8 :: nf !Canopy nitrogen factor (dimensionless) (Kull and Kruijt)

         !CASA parameters, CLM parameters
         real*8 :: sla !Specific leaf area (m^2 leaf area/kg C)
         real*8 :: r       !CLM respiration parameter (gC/gN)
         real*8 :: lrage   !CASA Turnover time of leaves and roots (years)
         real*8 :: woodage !CASA Turnover time of stems (years)
         real*8 :: lit_C2N !CASA litcn_casa (C:N ratio) IS THIS FOLIAGE&ROOTS?
         real*8 :: lignin  !CASA lignin (UNITS?  lignin content of ??)
         real*8 :: croot_ratio !Coarse roots:Stem mass woody ratio

         !* Parameters for plant allomteries ! not used in BiomeE
         real*8 :: b1Cf !para 1 for allometric relation btw DBH & foliage C
         real*8 :: b2Cf !para 2 for allometric relation btw DBH & foliage C
         real*8 :: b1Cd !para 1 for allometric relation btw DBH & dead C
         real*8 :: b2Cd !para 2 for allometric relation btw DBH & dead C
         real*8 :: b1Ht !para 1 for allometric relation btw DBH & height
         real*8 :: b2Ht !para 2 for allometric relation btw DBH & height

      end type PlantSP

!****************************************************************************
      type pftype
         !Parameters specific to plant functional types
         integer :: pst ! Photosynth type 1=C3, 2=C4
         logical :: woody !Woody, FALSE=NO, TRUE=YES.
         integer :: leaftype !1=broadleaf, 2=needleleaf, 3=monocot (not crops)
         real*8 :: hwilt  !Wilting point matric potential (m)
         real*8 :: sstar  !Rel. soil moist at stress onset (Rodriguez-Iturbe)
         real*8 :: swilt  !Normalized soil water at wilting point (dim'less)
         real*8 :: nf !Canopy nitrogen factor (dimensionless) (Kull and Kruijt)
         !CASA parameters, CLM parameters
         real*8 :: sla !Specific leaf area (m^2 leaf area/kg C)
         real*8 :: r      !CLM respiration parameter (gC/gN)
         real*8 :: lrage !CASA Turnover time of leaves and roots (years)
         real*8 :: woodage !CASA Turnover time of stems (years)
         real*8 :: lit_C2N !CASA litcn_casa (C:N ratio) IS THIS FOLIAGE&ROOTS?
         real*8 :: lignin  !CASA lignin (UNITS?  lignin content of ??)
         real*8 :: croot_ratio !Coarse roots:Stem mass woody ratio

         !Phenology parameter - KIM
         !* Parameter for phenology
         !phenotype - phenological type
         !            evergreen (1),
         !            cold deciduous (2),
         !            drought deciduous (3),
         !            cold/drought deciduous (4),
         !            annual  (5)
         integer :: phenotype !phenological types
         !* Parameters for plant allomteries
         real*8 :: b1Cf !para 1 for allometric relation btw DBH & foliage C
         real*8 :: b2Cf !para 2 for allometric relation btw DBH & foliage C
         real*8 :: b1Cd !para 1 for allometric relation btw DBH & dead C
         real*8 :: b2Cd !para 2 for allometric relation btw DBH & dead C
         real*8 :: b1Ht !para 1 for allometric relation btw DBH & height
         real*8 :: b2Ht !para 2 for allometric relation btw DBH & height

         ! from pspartype, Weng, 02/08/2021
         real*8 :: PARabsorb     !Leaf PAR absorptance (fraction)
         !Photosynthesis/Conductance - Farquhar/Ball-Berry parameters
         real*8 :: Vcmax         !Maximum photosynthetic capacity (umol m-2 s-1)
         real*8 :: m             !Slope of Ball-Berry equation
         real*8 :: b             !Intercept of Ball-Berry equation (mol m-2 s-1)
         real*8 :: Nleaf         !g-N/m2[leaf]

      end type pftype
!****************************************************************************

      type cohort
         type(entcelltype),pointer :: cellptr => null() !Pointer to ent grid cell
         type(patch),pointer :: pptr    => null() !Pointer to patch
         type(cohort),pointer :: taller => null() !Pointer to next tallest cohort
         type(cohort),pointer :: shorter => null() !Pointer to next shortest cohort
         type(cohort),pointer :: csptaller => null() !Pointer to next taller conspecific
         type(cohort),pointer :: cspshorter => null() !Pointer to next shorter conspecfic
         integer :: chID       ! cohort ID
         integer :: pft        ! PFT number
         real*8 :: n           ! Density of individuals in cohort (#/m^2)
         real*8 :: age         ! cohort age (year), Weng 11-10-2016
         real*8 :: leafage     ! mean leaf age in a crown (year), Weng, 2021-09-20)
         real*8 :: NSCstar     ! Target C_lab, for contoling growth and resp
         integer :: layer=1    ! crown position in canopy, Weng 12-19-2016
         !* PFT PARAMETERS
         ! Only need to index array of pftypes.

         !* NITROGEN status */
         !@var nm   Mean cohort nitrogen (g/m2[leaf]) over whole grid cell
         real*8 :: nm
         !@var Ntot Total cohort nitrogen (g/m[ground]2).
         real*8 :: Ntot
         !@var LAI Total cohort leaf area index (m2[leaf]/m2[ground])
         real*8 :: LAI               !*
         !@var LMA Leaf mass per leaf area (gC/m2)
         real*8 :: LMA
         !real*8 :: LA            ! Leaf area (m2[leaf]/individual)

         !* ALL QUANTITIES BELOW ARE FOR AN INDIVIDUAL *!

         !* GEOMETRY - trees:  GORT ellipsoids, grasses:leaf only
         real*8 :: h              !* Height (m)
         real*8 :: crown_dx       ! Crown horizontal radius (m)
         real*8 :: crown_dy       ! Crown vertical radius (m)
         real*8 :: dbh            ! Stem diameter at breast height (cm)
         real*8 :: root_d         ! Root half spheroid diameter (m)
         real*8 :: clump          ! Leaf clumping parameter (TBA)
         real*8 :: leafarea   ! m2, Weng, 2016-11-07
         real*8 :: crownarea  ! m2, Weng, 2016-11-07
         real*8 :: CrownLAI ! leafarea/crownarea, Weng, 2016-11-29
         real*8 :: fracroot(N_DEPTH) ! Fraction of roots in soil layer

         !* BIOMASS POOLS (KgC/single plant) ! Ensheng will change them all to kgC/plant
         real*8 :: C_fol          ! Foliage carbon
         real*8 :: C_sw           ! Sapwood carbon (=rho_wood*sw_vol) (units?)
         real*8 :: C_hw           ! Dead stem (heartwood) carbon
         real*8 :: C_croot        ! Coarse root carbon
         real*8 :: C_froot        ! Fine root carbon
         real*8 :: C_lab          ! Labile stored carbon
         real*8 :: C_seed=0.      ! Seed carbon,     Weng 11-07-2016

         real*8 :: N_lab          ! Labile stored nitrogen
         real*8 :: N_fol          ! Foliage nitrogen
         real*8 :: N_sw           ! Sapwood nitrogen
         real*8 :: N_hw           ! Dead stem (heartwood) nitrogen
         real*8 :: N_froot        ! Fine root nitrogen
         real*8 :: N_croot        ! Coarse root nitrogen
         real*8 :: N_branch       ! Branch nitrogen, Weng 11-07-2016
         real*8 :: N_seed         ! Seed nitrogen,   Weng 11-07-2016

         !* growth and ecological dynamics rate variable
         real*8 :: mu             ! Mortality rate,  Weng 11-07-2016, fraction of population per year

         !* FLUXES (for whole cohort over area cover)
         real*8 :: Ci             !*Internal foliage CO2 (mol/m3) !!Cohort level
         real*8 :: gcanopy        ! Conductance of water vapor/cohort (m/s)
         real*8 :: GPP            ! GPP flux/cohort/area cover (kg-C/m2/s)
         real*8 :: IPP            ! Isoprene emission flux (kg-C/m2/s)
         real*8 :: NPP            ! NPP flux/cohort/area cover (kg-C/m2/s)
         real*8 :: R_auto         ! Autotrophic respiration/cohort/area (kg-C/m2/s)
                                  ! = growth(Acan) + maint(fol,sapwood,root)
         real*8 :: R_root         ! Root respiration/cohort/area  (kg-C/m2/s) -PK 5/15/07
         real*8 :: leafRd        ! !umol individual-1 s-1
         !real*8 :: NSNmax
         real*8 :: N_up           ! N uptake from soil/cohort/area (kg-N/m2/s)
         real*8 :: N_uptake       ! N uptake from soil (kgN/tree/time step)
         real*8 :: C_litter       ! C in litterfall
         real*8 :: N_litter       ! N in litterfall
         real*8 :: C_to_Nfix      ! Carbon flux to N fixers symbionts

         !* PHENOLOGY - KIM
         real*8 :: phenofactor    = 0.d0 !phenofactor_c * phenofactor_d
         real*8 :: phenofactor_c  = 0.d0 !Cold deciduousness
         real*8 :: phenofactor_d  = 0.d0 !Drought deciduousness
         integer :: phenostatus    = 1
         real*8 :: betad_10d ! 10-day running average of betad
         real*8 :: CB_d !daily carbon balance
         real*8 :: turnover_amp
         real*8 :: llspan
         real*8 :: Sacclim ! state of acclimation/frost hardiness [deg C]
         !* Weng, 01/29/2021
         real*8 :: gdd = 0.d0
         real*8 :: sgdd = 0.d0
         real*8 :: ncd = 0.d0 ! number of cold days in non-growing season
         real*8 :: ngd = 0.d0 ! number of growing days
         real*8 :: alt = 0.d0    ! accumulative low temperature
         real*8 :: phenoF = 0.d0 ! foliage growth & fall status
         logical :: PhenoON = .False.       ! Pheno ON or OFF

         !* PHYSIOLOGICAL STATUS *!  !NYK
         real*8 :: stressH2O !* fraction stress factor, 0=stressed, 1=no stress
         real*8 :: stressH2Ol(N_DEPTH) ! Water stress in layers.
         real*8 :: senescefrac  !Net fraction of foliage that is litterfall.
         !* Additional C accounting
         real*8 :: C_growth  !* Daily tissue growth respiration (kg-C/m2-cohort/day)
                             !*  Save C_growth to restart to distribute flux over the day.
                             ! - this is remainig respiration carbon
         real*8 :: R_growth !respiration flux, kgC s-1 tree-1 
         real*8 :: C_growth_flux
         real*8 :: C_for_growth   ! C from labile C that will be used for building plant tissues
         real*8 :: resp_growth    ! Resp cost of growth = 0.33333 * C_for_growth, Weng 2016-11-18
         real*8 :: Ctot   !plant total C, kgC/tree
!        for reporting, Weng 12/02/2016
         real*8 :: GPP_daily=0.        ! GPP flux/tree (kg-C/tree/day)
         real*8 :: NPP_daily=0.        ! NPP flux/tree (kg-C/tree/day)
         real*8 :: Ra_daily=0.
         real*8 :: GPP_yearly=0.        ! GPP flux/tree (kg-C/tree)
         real*8 :: NPP_yearly=0.        ! NPP flux/tree (kg-C/tree)
         real*8 :: Ra_yearly=0.
         real*8 :: N_up_yr = 0.0        ! Annual N uptake(kg-N/tree)

      end type cohort

!****************************************************************************
      type canradtype
!         !Arrays in height levels in the canopy
         real*8,pointer :: h_lai(:)  => null()  !(m) height levels according to cumulative LAI
         real*8,pointer :: h_coh(:)  => null()  !(m) height levels acoording to cohorts
         real*8,pointer :: heights(:) => null() !(m) height levels at canopy layer boundaries (#layers+1)
         real*8,pointer :: LAI(:) => null()     !LAI within height level
         real*8,pointer :: f_sun(:) => null() !sunlit fraction profile
         real*8,pointer :: f_sha(:) => null() !shaded fraction profile
         real*8,pointer :: T_sun(:) => null() !sunlit PAR transmittance profile
         real*8,pointer :: T_sha(:) => null() !shaded PAR transmittance profile
         real*8,pointer :: I_sun(:) => null() !sunlit PAR absorbance profile
         real*8,pointer :: I_sha(:) => null() !shaded PAR absorbance profile
!         !Whole-canopy foliage clumping factor
         real*8 :: GORTclump
         real*8,pointer :: crad_heights(:) => null() !(m) height levels for incident light levels
      end type canradtype

!****************************************************************************
      type patch
         real*8 :: age                !*Patch age (years)
         real*8 :: area               !*Patch area (fraction of entcell)
         type(entcelltype),pointer:: cellptr => null() !Pointer to grid cell
         type(patch),pointer :: older => null() !Pointer to next older patch
         type(patch),pointer :: younger=> null()  !Pointer to next younger patch
         type(cohort),pointer :: tallest => null() !Pointer to tallest cohort
         type(cohort),pointer :: shortest => null() !Pointer to shortest cohort
         type(climate_env) :: ptenv ! climate envelope for PFTs ! Weng, 07/30/2019
         integer :: nCohorts = 0  ! How many cohorts in this patch, Weng 11-08-2016
         integer :: doy      = 1  ! Weng: for monitoring days, temporary

         !*- - - - - - - Cohorts summary variables - - - - - - - - - - -*!
         !  Intensive properties (e.g. geometry, LMA) are averages weighted by
         ! total number of individuals or LAI.
         !  Extensive properties (e.g. biomass, Ntot) are totals per m2 ground

         !* DIAGNOSTICS - NITROGEN and LEAF status */
         !@var nm   Mean leaf nitrogen (g/m2[leaf])
         real*8 :: nm
         !@var Ntot Total cohort nitrogen (g/m[ground]2).
         real*8 Ntot
         !@var LMA Leaf mass per leaf area (gC/m2)
         real*8 LMA
         !@var LAI Total cohort leaf area index (m2[leaf]/m2[ground])
         real*8 LAI
         real*8 LAIyr ! Yearly maximumu LAI, for diagnostics
         !@var LAIpft LAI by cover type.
         real*8,pointer :: LAIpft(:) => null() !(N_COVERTYPES)

         ! Weng 12-19-2016, patch summary
         real*8   :: totalCA     ! total crown area (crown area index), Weng 12-19-16
         real*8   :: LAIlayer(N_CAImax) ! for the LAI of each crown layer, Weng 05-05-17
         integer  :: crownlayers ! Crown layers (ceiling(totalCA))

         !* GEOMETRY - trees:  GORT ellipsoids, grasses:leaf only
         real*8 :: h              !* mean Height (m)
         real*8 :: crown_dx       ! Crown horizontal radius (m)
         real*8 :: crown_dy       ! Crown vertical radius (m)
!         real*8 :: dbh           ! Stem diameter at breast height (m)
!         real*8 :: root_d         ! Root half spheroid diameter (m)
         real*8 :: clump          ! Leaf clumping parameter (TBA)
         real*8,pointer :: fracroot(:) => null() ! Fraction of roots in soil layer

         !* DIAGNOSTICS - BIOMASS POOLS (kg-C/m^2-ground = sum kg-C/cohorts)
         real*8 :: C_fol          ! Foliage carbon (=LMA*LAI = kgC/m2-gnd)
         real*8 :: N_fol          ! Foliage nitrogen (gN/m2-gnd)
         real*8 :: C_w           ! Sapwood+hardwood carbon (=rho_wood*sw_vol) (units?)
         real*8 :: N_w           ! Sapwood+hardwood nitrogen
         real*8 :: C_lab          ! Labile stored carbon
         real*8 :: N_lab          ! Labile stored nitrogen
         real*8 :: C_froot        ! Fine root carbon
         real*8 :: N_froot        ! Fine root nitrogen
         real*8 :: C_root        ! Fine+coarse root carbon
         real*8 :: N_root        ! Fine+coarse root nitrogen
         real*8 :: C_seed
         real*8 :: N_seed

         !* EXPORT - FLUXES (whole patch)
         real*8 :: Ci             !*Internal foliage CO2 (mol/m3)
         real*8 :: GCANOPY        ! Conductance of water vapor (m/s)
         !* DIAGNOSTICS - FLUXES
         real*8 :: GPP            ! GPP flux (kg-C/m2/s)
         real*8 :: IPP            ! Isoprene flux (kg-C/m2/s)
         real*8 :: NPP            ! NPP flux (kg-C/m2/s)
         real*8 :: R_auto         ! Autotrophic respiration (kg-C/m2/s)
                                  ! = growth(Acan) + maint(fol,sapwood,root)
         real*8 :: R_root         ! Root respiration (kg-C/m2/s) -PK 5/15/07
         real*8 :: N_up           ! N uptake from soil(kg-N/m2/s)
         real*8 :: GPP_daily,NPP_daily,Ra_daily
         real*8 :: GPP_yearly,NPP_yearly,Ra_yearly
         real*8 :: C_litter       ! C in litterfall
         real*8 :: N_litter       ! N in litterfall
!         real*8 :: C_to_Nfix      ! Carbon flux to N fixers symbionts
         !- - - - - - - end of cohort summary variables - - - - - - - - - - - - -

         !- - - - - - - Patch total - - - - - - - - - - - - - - - - - - - - - - -

         !* IMPORT-PRESCRIBED, EXPORT-SIMULATED
         real*8 :: z0              !Roughness length, average over patch
         !* EXPORT
         real*8 :: albedo(N_BANDS) !total given Id,Ii, snow-free except grd
         real*8 :: albedodir(N_BANDS) !direct black sky, snow-free except grd
         real*8 :: albedodif(N_BANDS) !diffuse white sky, snow-free except grd
         real*8 :: betad             !Water stress  # CALC FROM Soilmoist & SSTAR by PFT
         real*8,pointer :: betadl(:) => null() !Water stress in layers.
         real*8 :: TRANS_SW          !Transmittance of shortwave radiation (400-2500 nm) to the ground (fraction)
         real*8 :: CO2flux           !Net CO2 flux up (kg-C/m2-gnd/s)
         !* DIAGNOSTICS - soil
         real*8 :: Soil_resp         !soil resp flux (kg-C/m2/s) -PK 6/14/06,changed umol to kg-NK 07/28/06
         real*8 :: Rh_daily = 0.d0  ! Rh per day, kgC m-2 day-1
         real*8 :: Rh_yearly = 0.d0 ! Rh per year, kgC m-2 yr-1
         ! soil C and N pools: 1 Carbon, 2 Nitrogen;  1 Surfmet, 2 Surfstr, 3 Soilmet, 4 Soilstr, 5 CWD, 6 Surfmic, 7 Soilmic, 8 Slow, 9 Passive
         real*8, dimension(PTRACE,NPOOLS,N_CASA_LAYERS) :: Tpool !(g-C/m^2, CASA Tpools, single cell) !added dim N_CASA_LAYERS -PK

         !* IMPORT - Variables calculated by GCM/EWB - downscaled from grid cell

      !use soil moisture (and temperature) for 2 CASA layers:    -PK
      !0-30 cm, 30-100 cm (second might change to 30-200 cm)
      !**might change this and soiltemp to dynamically allocated arrays** -PK 7/07
!         real*8 :: Soilmoist(N_CASA_LAYERS) !Soil moisture (volumetric fraction)
         real*8 :: Soilmoist(N_DEPTH) !Soil moisture (rel. sat. vol. fraction)

!         real*8 :: N_deposit    !N deposition (kgN/m2)

         !* Variables for biophysics and biogeochemistry
         type(canradtype) :: crad !Data structure for light profile

         !* Reproductive pools *!
         real*8,pointer :: Reproduction(:) => null() !Reproductive/seed pools array of length N_PFT (kgC/m2-patch)
             ! Weng: I don't use this array because a patch does not have enough infomration of PFTs

         !* Disturbance values
         real*8 :: fuel
         real*8 :: ignition_rate
         real*8 :: lambda1(T_SUB) !Site-averaged fire dist. rate during year
         real*8 :: disturbance_rate(N_DIST_TYPES)

         !* Soil data (needed for albedo computation)
         integer soil_type      ! 1 - sand (bright) ; 2 - dirt (dark)
#ifdef ENT_ALBEDO_1
         real*8 :: soil_albedo(N_BANDS) ! albedo of soil only
#endif


#ifdef NEWDIAG
        !* Rates of change - patch total
         real*8 :: dadt              !Rate of change of patch age = 1
         real*8 :: dpdt              !Rate of change of patch area
         real*8 :: dwdt              !Rate of change of available soil water
#endif
         real*8 :: SOC(N_SOM)
         real*8 :: SON(N_SOM) ! finleL, structuralL, micro, fast, slow
         real*8 :: mineralN ! Mineralized N (kgN/m2)

         real*8 :: N_min_day = 0.        ! N mineralized per day
         real*8 :: N_min_yr = 0.         ! N mineralized per year
         real*8 :: N_uptake = 0.         ! total N absorbed by roots at a time step
         real*8 :: N_up_day = 0.         ! Daily N uptake
         real*8 :: N_up_yr = 0.          ! Yearly N uptake
         real*8 :: previousN = -1.0      ! For the first year
         real*8 :: annualN
         ! diags and hacks
         real*8 :: Ctot  !plant total carbon kg-C/m2 ! Weng
         real*8 :: C_growth
      end type patch


!****************************************************************************
      type entcelltype
!        real*8 :: long, lat      !longitude, latitude
!         integer :: longi, latj    !grid cell i,j
         real*8 :: area         !Area km^2
         type(patch), pointer:: youngest => null()
         type(patch), pointer:: oldest => null()
         integer :: nPatches ! number of patches in a cell grid ! Weng, 2020-03-13
         integer :: ijdebug ! Weng, 08/28/2020
         integer :: year0=0  ! current year
         integer :: year1=0  ! receiving the year of ModelEclock
         integer :: doy  =0  ! dayOfYear
         integer :: hemi =1 ! 1: north, 2: south
         integer :: steps_d ! steps after last daily update ! Weng, 12/14/2020
         logical :: BiomeE_Veg_not_initialized = .true.
         !*- - - - - - - Cohorts summary variables - - - - - - - - - - -*!
         !* Per vegetated ground area of entcell ** excludes bare soil area.
         !  Intensive properties (e.g. geometry, LMA) are averages weighted by
         ! total number of individuals or leaf area.
         !  Extensive properties (e.g. biomass, Ntot) are totals per m2 ground

         !* IMPORT-PRESCRIBED, EXPORT-SIMULATED - NITROGEN and LEAF status */
         !@var nm   Mean leaf nitrogen (g/m2[leaf]).
         real*8 :: nm
         !@var Ntot Total cohort nitrogen (g/m[ground]2) for vegetated patches only.
         real*8 :: Ntot
         !@var LMA Leaf mass per leaf area (gC/m2)
         real*8 :: LMA
         !@var LAI Leaf area index (m2[leaf]/m2[ground]) for vegetated patches only.
         real*8 :: LAI
         !@var LAIpft LAI by cover type.
         real*8,pointer :: LAIpft(:) => null()!(N_COVERTYPES)
         !real*8 :: LA            ! Leaf area (m2[leaf]/individual)

         !* IMPORT-PRESCRIBED, EXPORT-SIMULATED - GEOMETRY - trees:  GORT ellipsoids, grasses:leaf only
         real*8 :: h              !* mean Height (m)
!         real*8 :: crown_dx       ! Crown horizontal radius (m)
!         real*8 :: crown_dy       ! Crown vertical radius (m)
!         real*8 :: dbh            ! Stem diameter at breast height (m)
!         real*8 :: root_d         ! Root half spheroid diameter (m)
!         real*8 :: clump          ! Leaf clumping parameter (TBA)
         real*8,pointer :: fracroot(:) => null() ! Fraction of roots in soil layer

         !*  IMPORT-PRESCRIBED, EXPORT-SIMULATED - BIOMASS POOLS (g-C/m^2)
         real*8 :: C_fol          ! Foliage carbon (=LMA*LAI = kgC/m2-gnd)
         real*8 :: N_fol          ! Foliage nitrogen (gN/m2-gnd)
         real*8 :: C_w           ! Sapwood+dead wood carbon (=rho_wood*sw_vol) (units?)
         real*8 :: N_w           ! Sapwood+dead wood nitrogen
         real*8 :: C_lab          ! Labile stored carbon
         real*8 :: N_lab          ! Labile stored nitrogen
         real*8 :: C_froot        ! Fine root carbon
         real*8 :: N_froot        ! Fine root nitrogen
         real*8 :: C_root        ! Fine+coarse root carbon
         real*8 :: N_root        ! Fine+coarse root nitrogen

         !* IMPORT/EXPORT PUBLIC - FLUXES)
         real*8 :: Ci             !*Internal foliage CO2 (mol/m3)
         real*8 :: GCANOPY        ! Conductance of water vapor (m/s)
         !* EXPORT - FLUXES
         real*8 :: GPP            ! GPP flux (kg-C/m2/s)
         real*8 :: IPP            ! Isoprene flux (kg-C/m2/s)
         real*8 :: NPP            ! NPP flux (kg-C/m2/s)
         real*8 :: R_auto         ! Autotrophic respiration (kg-C/m2/s)
                                  ! = growth(Acan) + maint(fol,sapwood,root)
         real*8 :: R_root         ! Root respiration (kg-C/m2/s) -PK 5/15/07
         real*8 :: N_up           ! N uptake from soil(kg-N/m2/s)
!         real*8 :: C_litter       ! C in litterfall
!         real*8 :: N_litter       ! N in litterfall
!         real*8 :: C_to_Nfix      ! Carbon flux to N fixers symbionts
         !- - - - - - - end of cohort summary variables - - - - - - - - - - - - - - - - - - -

         !- - - - - -  Patch-level summary values - PHYSICAL ------------------
         !* EXPORT - from radiative transfer
         real*8 :: z0              !Roughness length, average over patch
         real*8 :: albedo(N_BANDS) !Albedo total dir + dif (black sky + white sky)
         real*8 :: albedodir(N_BANDS) !Albedo direct-beam black sky component.
         real*8 :: albedodif(N_BANDS) !Albedo diffuse white sky component.
         real*8 :: betad             !Water stress  # CALC FROM Soilmoist & SSTAR by PFT
         real*8,pointer :: betadl(:) => null() !Water stress in layers.
         real*8 :: TRANS_SW  !Transmittance of shortwave radiation to the ground (fraction)
         real*8 :: CO2flux           !Net CO2 flux up (kg-C/m2-gnd/s)
         !* DIAGNOSTICS - soil
         real*8 :: Soil_resp         !soil resp flux (kg-C/m2/s) -PK 6/14/06,changed umol to kg-NK 07/28/06
         real*8, dimension(PTRACE,NPOOLS,N_CASA_LAYERS) :: Tpool !(g-C/m^2, CASA Tpools, single cell) !added dim N_CASA_LAYERS -PK

         !* Disturbance values
         real*8 :: fuel
         real*8 :: ignition_rate
         real*8 :: lambda1(T_SUB) !Site-averaged fire dist. rate during year
         real*8 :: disturbance_rate(N_DIST_TYPES)

         !- - - - - - Entcell-level variables - - - - - - - - - - - - - - - -
         !* IMPORT/EXPORT - vegetation
         real*8 :: fv            ! vegetation fraction of entcell
         real*8 :: heat_capacity ! total veg. heat capacity
         real*8 :: fwet_canopy   ! fraction of canopy that is wet
         !* IMPORT - SOIL
         real*8 :: soil_Phi      !Soil porosity (m3/m3)
         real*8 :: soil_dry     !Soil wetness "when evapotranspiration stops"
         real*8 :: soildepth    !Soil depth (m)
         real*8 :: theta_max    !Saturated soil water volume (m/m)
         real*8 :: k_sat        !Saturated hydraulic conductivity
         real*8 :: root_Phi     !Infiltration factor promoted by roots (units?)

         !SOIL - CONSTANTS
         !Soil textures for CASA -PK
         real*8 :: soil_texture(N_SOIL_TEXTURES) ! fractions of soil textures
!         real*8 clayfrac  !fractional clay content (passed from GHY.f)
!     real*8 sandfrac  !fractional sand content (also from GHY.f)
#ifdef ENT_ALBEDO_1
         real*8 :: soil_albedo(N_BANDS) !Always set, but used only with ENT_ALBEDO
         real*8 :: snow_albedo(N_BANDS,2) !set only for ENT_ALBEDO
#endif
         !IMPORT - METEOROLOGICAL STATE VARIABLES
         !Cell-level summary values - CALCULATED BY GCM/EWB OR OFF-LINE FILE
         real*8 :: TairC ! Air temperature (Clesius) !KIM-to drive phenology
         real*8 :: Precp ! mm/s, Precipitation, Weng for Biogeography, 07/29/2019
         real*8 :: TcanopyC     !Canopy temperatue (Celsius)
         real*8 :: Qf           !*Foliage surface vapor mixing ratio (kg/kg)
         real*8 :: P_mbar       !Atmospheric pressure (mb)
         real*8 :: Ca           !@Atmos CO2 conc at surface height (mol/m3).
         !next two now explicitly depth-structured (see above) -PK
!         real*8 :: Soilmoist(N_CASA_LAYERS) !Soil moisture (volumetric fraction)
!         real*8 :: Soiltemp(N_CASA_LAYERS)  !Soil temperature (Celsius)
         real*8 :: Soilmoist(N_DEPTH) !Soil moisture (volumetric fraction)
         real*8 :: Soiltemp(N_DEPTH)  !Soil temperature (Celsius)
         real*8,pointer :: Soilmp(:) => null() !Soil matric potential (m)
         real*8,pointer :: fice(:) => null() !Fraction of soil layer that is ice
#ifdef ENT_ALBEDO_1
!         real*8 :: snow_fracb,snow_fracv !Snow cover fraction, 1-bare, 2-veg
!         real*8 :: snow_depthb,snow_depthv !Snow depth (m),  1-bare, 2-veg
         real*8 :: snow_frac(2) !Snow cover fraction, 1-bare, 2-veg
         real*8 :: snow_depth(2) !Snow depth (m), 1-bare, 2-veg
#endif
         real*8 :: Ch           !Ground to surface heat transfer coefficient
         real*8 :: U            !Surface layer wind speed (m s-1)


         !Radiation - IMPORT STATE VARIABLES
         !may later be broken down into hyperspectral increments.
         ! in an array
!         real*8 :: Ivis          !Incident visible  (W m-2)
!         real*8 :: Idir          !Incident direct visible  (W m-2)
!         real*8 :: IPAR         !Incident PAR 400-700 nm (W m-2)
         real*8 :: IPARdir        !Incident direct PAR (W m-2)
         real*8 :: IPARdif        !Incident diffuse PAR (W m-2)
         real*8 :: CosZen         !cos of solar zenith angle

         !PHENOLOGY - KIM
         real*8 :: soiltemp_10d
         real*8 :: airtemp_10d
         real*8 :: paw_10d
         real*8 :: par_10d
         real*8 :: gdd
         real*8 :: ncd
         real*8 :: sgdd
         real*8 :: daylength(2) !previous & present day
         integer :: fall
         type(climate_env) :: ptenv ! climate envelope for PFTs ! Weng, 07/30/2019

         ! diags and hacks
         real*8 :: Ctot   !plant total C, kgC/tree
         real*8 :: C_growth

!!! hack !!! - just to try master phenology
         real*8 :: ld
         real*8 :: light
      end type entcelltype


!****************************************************************************
!      type entdatatype
!        real longmin, longmax, latmin, latmax
!        integer longi,latj
!        type(timestruct),pointer :: tt      !Greenwich Mean Time
!        type (entcelltype),pointer :: grid(:,:)
!      end type entdatatype

!****************************************************************************
      type ent_config
      ! this type should contain all parameters that describe the run
      ! i.e. flags, array dimensions etc. They assumed to be constant
      ! during the run but may change from run to run
        logical do_soilresp       ! do soil respiration
        logical do_phenology_activegrowth
        logical do_PPA_dynamics_EW ! Weng, 04/16/2017
        logical do_demographic_EW ! Weng, 07/22/2019
        logical do_structuralgrowth
        logical do_frost_hardiness
        logical do_patchdynamics
        logical do_init_geo
!        logical mixed_veg
      end type ent_config

! For plant physiology
      !=====DECLARED TYPES======!
      !------------------------------------------------------------
      type :: canraddrv
         !Canopy radiative transfer
         !@var sigma Leaf scattering coefficient (?unitless).
         real*8 :: sigma        !=0.2D0
         real*8  sqrtexpr  !This just calculated from sigma.
         !@var kdf Canopy extinction coeff. for diffuse radiation (unitless).
         real*8 :: kdf          !=0.71D0
         !@var rhor Canopy reflectivity (?unitless).
         real*8 rhor
         !@var kbl Canopy extinction coeff. for black leaves (unitless).
         real*8 kbl

         !Vegetation geometry, biology
         !UPDATE WITH CANOPY RADIATIVE TRANSFER SCHEME
         integer :: pft
         !real*8 :: leafalbedo
         real*8 :: canalbedo
         real*8 :: LAI  !canopy LAI

         !Radiation
!         real*8 :: solarzen !Solar zenith angle (radians)
         real*8 :: CosZen !cos(solarzen) = sin(solarelev)
         real*8 :: I0df, I0dr   !Direct and diffuse incident PAR at top of canopy (umol m-2 s-1)
      end type canraddrv

      ! Parameters calculated from pftpar
      type photosynthpar       !Calculated values from pft-dependent pspartypes
      integer :: pft           !Plant functional type.  1-C3 grassland
      real*8 :: PARabsorb      !Leaf PAR absorptance (fraction)
      real*8 :: Vcmax           !Maximum photosynthetic capacity (umol m-2 s-1)
      real*8 :: Kc              !Michaelis-Menten parameter for CO2 (Pa)
      real*8 :: Ko              !Michaelis-Menten parameter for O2 (Pa)
      real*8 :: Gammastar       !CO2 compensation point (Pa)
      real*8 :: m               !Slope of Ball-Berry equation
      real*8 :: b               !Intercept of Ball-Berry equation (mol m-2 s-1)
      real*8 :: Nleaf           !g-N/m^2[leaf] - May want to take this from Tpool instead.
      real*8 :: stressH2O       !Water stress factor (fraction, 1=no stress)
      logical :: first_call     !For optimizing run
      real*8 :: Ac              !Save Ac for calculation only once per timestep.
      real*8 :: As              !Save As for calculation only once per timestep.
      logical :: reset_ci_cubic1
      end type photosynthpar

! Copied from "module FarquharBBpspar", Weng, 2021-01-19
      ! photosynthesis drivers
      type psdrvtype
        real*8 :: ca            !Surface CO2 mole fraction (umol mol-1)
        real*8 :: ci            !Leaf internal CO2 mole fraction (umol mol-1)
        real*8 :: Tc            !Canopy (foliage) temperature (Celsius)
        real*8 :: Pa            !Atmospheric pressure (Pa)
        real*8 :: rh            !Relative humidity (fraction)
      end type psdrvtype
      ! Photosythesis parameters
      type pspartype
        integer :: pst          !Photosynth type.  1-C3, 2=C4
        real*8 :: PARabsorb     !Leaf PAR absorptance (fraction)
        !Photosynthesis/Conductance - Farquhar/Ball-Berry parameters
        real*8 :: Vcmax         !Maximum photosynthetic capacity (umol m-2 s-1)
        real*8 :: m             !Slope of Ball-Berry equation
        real*8 :: b             !Intercept of Ball-Berry equation (mol m-2 s-1)
        real*8 :: Nleaf         !g-N/m2[leaf]
      end type pspartype

!=========================================================================
      end module ent_types
