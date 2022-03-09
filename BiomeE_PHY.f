#ifdef ENT_STANDALONE_DIAG
#define  DEBUG  1
#endif
!#define USE_NR_SOLVER_FOR_FBB
      module BiomeE_PHY !canopyspitters
      ! added by Ensheng Weng, 2018-11-03
      ! including physiology, photosynthesis, respiration
      ! from three files of Ent: canopyspitters.f, respauto_physio.f,
      ! FBBphotosynthesis.f, and phenology.f

      use ent_types
      use ent_const
      use ent_pfts
      use ent_debug_mod, only : ent_d
      use BiomeE_VEG, only: ijdebug
      implicit none

      save

      public clim_stats
      public BiomeE_CanopyModel
      public BiomeE_Phenology
      ! Parameters
      real*8,parameter :: IPARMINMOL=LOW_PAR_LIMIT  !umol m-2 s-1
      real*8,parameter :: O2frac=.20900 !fraction of Pa, partial pressure.

      !=====photosynthesis and respiration constants =====!
      real*8,parameter :: ciMIN = 1.d-8  !Small error
!      real*8,parameter :: Kelvin = 273.15d0
      real*8,parameter :: Rgas = 8.314510d0 !gasc, gas constant (8.314510 J/mol K)
      real*8,parameter :: Kc = 30.0       !Michaelis-Menten coeff.constant for CO2 (Pa), Collatz (1991)
      real*8,parameter :: Ko = 3.d4       !Michaelis-Menten coeff.constant for O2 (Pa), Collatz (1991)
      real*8,parameter :: KcQ10 = 2.1d0   !Kc Q10 exponent, Collatz (1991)
      real*8,parameter :: KoQ10 = 1.2d0   !Ko Q10 exponent, Collatz (1991)
      real*8,parameter :: gamma_r = 0.015d0 ! Rdark = 0.015d0 * Vcmax ! S. von Caemmerer (2000) Biochemical Models of Leaf Photosynthesis

      !*****************Phenology model constants **************************
      ! Only update phenology status, Ensheng Weng, 2018-11-04
      ! temperature constraint for cold-deciduous PFTs (Botta et al. 1997)
      ! T0_gdd5: base temperature to calculate growing degree days (gdd)
      ! gdd_par1/2/3: paramters to estimate the threshold for gdd
      ! gdd_threshold = gdd_par1 + gdd_par2*exp(gdd_par3*ncd)
      real*8, parameter :: T0_gdd5 = 2.d0 ! 5.d0
      real*8, parameter :: T0_chill = 5.d0
      real*8, parameter :: gdd_par1 = 50.d0   ! -68.d0
      real*8, parameter :: gdd_par2 = 800. ! 650.d0  !800.d0  ! 638.d0
      real*8, parameter :: gdd_par3 = -0.025d0 ! -0.01d0
      !*sgdd_threshold & length - tuning parameters (tuned for Barrow)
      real*8, parameter :: T0_sgdd = 0.d0
      real*8, parameter :: sgdd_threshold = 100.d0

      private
      !=====GLOBAL VARIABLES (MODULE ONLY)=====!
      type(photosynthpar) :: pspar !Photosynthesis parameters.

      contains

!================================================================
!================ Phenology =====================================
!================================================================
      subroutine BiomeE_Phenology(pp)
!@sum Update statstics for phneology_update at a DAILY time step.
!@sum Weng condensed it into one subroutine
!     updated by Weng, 2021-09-21
      use ent_const
      type(patch) :: pp
      !--Local-----
      type(cohort), pointer :: cop
      real*8 :: temp_10d, betad_10d !10 day running mean of temp & betad (calculated with stressH2O)
      real*8 :: gdd_c ! gdd or sgdd, depending on plant form
      real*8 :: T10d_crit, gdd_thld  ! GDD treshold (White et al.)
      !real*8 :: soil_beta   ! soil moisture index

      !soil_beta = SUM(pp%cellptr%Soilmoist(1:3))/3
      cop => pp%tallest
      do while(ASSOCIATED(cop))
         associate(sp =>BMEspdata(cop%pft))
         ! Climatic variables for phenology
         if(sp%woody)then
           temp_10d = pp%cellptr%airtemp_10d
           gdd_c   = cop%gdd
         else
           temp_10d = pp%cellptr%soiltemp_10d
           gdd_c   = cop%sgdd
         endif
         betad_10d = cop%betad_10d

         ! GDD threshold for leaf green-up
         gdd_thld = sp%gdd_par1 +
     &           sp%gdd_par2 * exp(sp%gdd_par3*cop%ncd)
         ! Critical temperature trigering offset of phenology
         T10d_crit = sp%Tc_pheno0 !Tmax0
     &        - 5.0 * exp(-0.05*max(-15.0,(cop%ngd-90.0)))
         ! Persistence of low temperature
         if(temp_10d > T10d_crit)then
            cop%alt=0.0
         elseif(cop%phenofactor_c>0.99)then
            cop%alt=Max(-999.0, cop%alt + temp_10d - T10d_crit)
         endif

         !-----------Calculate phenofactor_c: from 0 to 1---
         ! Leaf ON
         if (gdd_c.gt.gdd_thld.and.cop%phenofactor_c.lt.1.0)then
            cop%phenofactor_c=min(1.0,cop%phenofactor_c+0.1)
            if (cop%phenofactor_c == 1.d0) cop%ncd = 0.d0
            cop%alt = 0.d0
         endif
         ! Leaf OFF
         if((temp_10d.lt.T10d_crit) .and. (cop%alt<-30.))then
            cop%phenofactor_c=max(0.,cop%phenofactor_c-0.1)
            if (cop%phenofactor_c > 0.d0) then ! .eq. 0.d0
              ! zero gdd for the next cycle of thermal accumulation, 02/22/21
              ! before moving to the next day to prevent triggering leaf onset
              cop%gdd  = 0.d0
              cop%sgdd = 0.d0
            else
              cop%alt = 0.d0
              cop%ngd = 0.d0
            endif
         endif

         ! Calculate phenofactor_d with soil water drought index
         if(cop%phenostatus<3 .and. betad_10d>sp%betad_ON)then
            cop%phenofactor_d = min(1.d0,cop%phenofactor_d+0.1)
         elseif((cop%phenostatus==3 .and.betad_10d<sp%betad_OFF)
     &      .or.(cop%phenostatus==4 .and.cop%phenofactor_d<1.0))then
            cop%phenofactor_d = max(0.d0,cop%phenofactor_d-0.1)
         endif

         !---------Update phenostatus-------------
         ! phenostatus translated from phenofactor:
         !1 - no leaf (0);   2 - growing leaf (0 to 1)
         !3 - full leaf (1); 4 - leaf senescence (1 to 0)
         cop%phenofactor = MIN(cop%phenofactor_c,cop%phenofactor_d)

         if(cop%phenofactor > 0.99d0)then
            cop%phenostatus = 3
         elseif(cop%phenostatus==1 .and. cop%phenofactor>=0.05)then
            cop%phenostatus = 2
         elseif(cop%phenostatus==3 .and. cop%phenofactor < 0.99)then
            cop%phenostatus = 4
         elseif(cop%phenostatus==4 .and. cop%phenofactor < 0.01)then
            cop%phenostatus = 1
         endif

         end associate

!! for debug
!#ifndef SITE_TEST
!         if(ijdebug == 44067)then  ! 101061 44067 !Harvard Forest
!#endif
!           write(935,'(a60, ",", L4, ",", I8, ",",20(f8.2,","))')
!     &      'pheno,stat,fct,fct_c,fct_d,gdd,sgdd,alt,ncd,LAI,H',
!     &       cop%PhenoON,cop%phenostatus,cop%phenofactor,
!     &       cop%phenofactor_c,cop%phenofactor_d,
!     &       cop%gdd,cop%sgdd,cop%alt,cop%ncd,cop%lai,cop%h
!#ifndef SITE_TEST
!         endif
!#endif

         ! to the next cohort
         cop => cop%shorter
      enddo

      end subroutine BiomeE_Phenology

!================================================================
      subroutine clim_stats(dtsec,ecp,config,dailyupdate)
!@sum Calculate climate statistics such as 10 day running mean
!@+   Called by ent.f every physical time step.
      use BiomeE_VEG, only : Soillayer_convert_Ent
      real*8,intent(in)  :: dtsec           !dt in seconds
      type(entcelltype) :: ecp
      type(ent_config)  :: config
      logical, intent(in) :: dailyupdate
      !-----local--------
      type(patch), pointer :: pp
      type(cohort), pointer :: cop
      !*local variables for entcell-level envrionment variable
      real*8 :: airtemp        !air temperature degC
      real*8 :: soiltemp       !soil temperature degC
      real*8 :: par            !photosynthetic active radiation (PAR)
      !*local variables for entcell-level variables,
      !*updated in this subroutine
      real*8 :: Tair_10d    !10 day running mean of air temperature
      real*8 :: Tsoil_10d   !10 day running mean of soil temperature
      real*8 :: par_10d        !10 day running mean of PAR
      real*8 :: gdd            !growing degree days, based on air temperature
      real*8 :: ncd            !number of chilling days, based on air temperature
      real*8 :: sgdd           !growing degree days, based on soil temperature
      !*PAR-limited phenology parameters
      real*8 :: par_crit       !PAR threshold
      logical :: par_limit     !logical whether PAR-limited phenology parameterization is applied or not for certain PFTs
      real*8 :: turnover0      !turnover amplitude, calculated with the phenology parameterization
      real*8 ::  llspan0       !leaf life span, calculated with the phenology parameterization
      !*soil temperature for CASA layers
      real*8 :: Soiltemp2layer(N_CASA_LAYERS)


      airtemp = ecp%TairC
      call Soillayer_convert_Ent(ecp%Soiltemp(:), SOILDEPTH_m,
     &     Soiltemp2layer)
      soiltemp = Soiltemp2layer(1)

      Tsoil_10d = ecp%soiltemp_10d
      Tair_10d = ecp%airtemp_10d
      par_10d = ecp%par_10d
      gdd = ecp%gdd
      ncd = ecp%ncd
      sgdd = ecp%sgdd

      !*10-day running mean of Air Temperature
      Tair_10d = running_mean(dtsec, 10.d0, airtemp, Tair_10d)

      !*10-day running mean of Soil Temperature
      Tsoil_10d = running_mean(dtsec, 10.d0, soiltemp, Tsoil_10d)

      !*10-day running mean of PAR
      par = ecp%IPARdif + ecp%IPARdir !total PAR is the sum of diffused and direct PARs
      par_10d = running_mean(dtsec, 10.d0, par, par_10d)

      !*daylength
      if (ecp%CosZen > 0.d0) then
         ecp%daylength(2) = ecp%daylength(2) + dtsec/60.d0
      end if

      !*GDD & NCD - Update Once a day
      if (dailyupdate .and. ecp%steps_d > 12)then
         !!!! cohort-specific gdd, sgdd, and ncd, Weng, 02/03/2021
         pp=>ecp%oldest
         do while (ASSOCIATED(pp))
           cop => pp%tallest
           do while(ASSOCIATED(cop))
            ! Thermal accumulate regardless of pheno status
            cop%gdd = cop%gdd +MAX(0.d0,(Tair_10d - T0_gdd5))
            cop%sgdd= cop%sgdd+MAX(0.d0,(Tsoil_10d- T0_sgdd))
            cop%gdd = MIN(cop%gdd, 9125.d0 ) ! just to avoid a huge number
            cop%sgdd= MIN(cop%sgdd, 9125.d0) ! 365 * 25
            if (cop%PhenoON) then
             cop%ngd = Min(366., cop%ngd + 1.)
            else ! cop%PhenoON = .false.
              if(Tair_10d < T0_chill)then
                 cop%ncd = cop%ncd + 1.d0
                 ! Keep gdd and sgdd zero in early chilling period
                 if(cop%ncd < 20.d0)then
                   cop%gdd  = 0.d0
                   cop%sgdd = 0.d0
                 endif ! ncd
              endif
            endif ! PhenoON

            cop => cop%shorter
           enddo ! all cothors in a patch
           pp => pp%younger
         enddo ! all patches in a cell

         !**************************************************************
         !*If the season is fall or not (if fall, it's 1; else, it's 0)
         !1) it is to control the phenological status
         !2) it is determined according to the daylength increases or not
         if (NInt(ecp%daylength(2)) .lt. NInt(ecp%daylength(1)))then
            ecp%fall = 1
         elseif (NInt(ecp%daylength(2)).gt.NInt(ecp%daylength(1)))then
            ecp%fall = 0
         end if

      endif  ! daily update

      !============================================================
      !!!! fast step calculation, dtsec
      !============================================================
      pp => ecp%oldest
      do while (ASSOCIATED(pp))
        cop => pp%tallest
        do while(ASSOCIATED(cop))
          !*10-day running mean of stressH2O (betad)
          cop%betad_10d = running_mean(dtsec, 10.d0,
     &                    cop%stressH2O, cop%betad_10d)

          !======== photosynthetic acclimation factor for evergreen veg
          if (config%do_frost_hardiness) then
             if (((BMEspdata(cop%pft)%phenotype.eq.EVERGREEN).and.
     &           (BMEspdata(cop%pft)%leaftype.eq.NEEDLELEAF)).or.
     &           (BMEspdata(cop%pft)%phenotype.eq.COLDDECID).or.
     &           (BMEspdata(cop%pft)%phenotype.eq.ARCTIC)) then
                call photosyn_acclim(dtsec,Tair_10d,cop%Sacclim)
             else
                cop%Sacclim = 25.d0 !Force no cold hardening, mild temperature.
             endif
          else
              cop%Sacclim = 25.d0 !Force no cold hardening, mild temperature.
          endif

          cop => cop%shorter
        enddo
        pp => pp%younger
      enddo

      ecp%soiltemp_10d = Tsoil_10d
      ecp%airtemp_10d  = Tair_10d
      ecp%par_10d = par_10d
      ecp%gdd = ecp%oldest%tallest%gdd
      ecp%ncd = ecp%oldest%tallest%ncd
      ecp%sgdd = ecp%oldest%tallest%sgdd

      end subroutine clim_stats

!========================================================================
      subroutine photosyn_acclim(dtsec,Ta,Sacc)
!@sum Calculate Sacclim accumulator state used to calculate
!@+   acclimation/frost hardiness for boreal coniferous forests.
!@+   Called by clim_stats during update of phenology climate stats.
!@+   Based on Repo et al (1990),
!@+   Hanninen & Kramer (2007),and  Makela et al (2006)
      implicit none
      real*8,intent(in) :: dtsec ! time step size [sec]
      real*8,intent(in) :: Ta ! air temperature [deg C]
      real*8,intent(inout) :: Sacc ! state of acclimation [deg C]

      !----Local-----
      real*8,parameter :: tau_inv = 2.22222e-6
                          ! inverse of time constant of delayed
                          ! response to ambient temperature [sec] = 125 hr
                          ! Makela et al (2004) for Scots pine

!      Use a first-order Euler scheme
       Sacc = Sacc + dtsec*(tau_inv)*(Ta - Sacc)

!     Predictor-corrector method requires temperature from next timestep
!       Sacc_old = Sacc
!       Sacc = Sacc_old + dtsec*(1/tau_acclim)*(Ta - Sacc )
!       Sacc = Sacc_old + ((1/tau_acclim)*(Ta - Sacc_old)+
!     &                     (1/tau_acclim)*(Ta_next - Sacc))*0.5d*dtsec

      end subroutine photosyn_acclim

!============================================================================
!=============== Main photosynthesis ========================================
!============================================================================
      subroutine BiomeE_CanopyModel(pp,dtsec)
!@sum photosynth_cond  Main routine to set up drivers and calculate
!@sum canopy scale fluxes.
!@+   Version that calls Farquhar-Ball-Berry leaf biophysics.
!@+   Calculates photosynthesis, autotrophic respiration, H2O conductance,
!@+   looping through cohorts.
!@+   Inputs:  met drivers, radiation, Ca, Tcanopy
!@+   Outputs:  GPP, NPP, autotrophic respiration components
!@+   Use PPA's assumptions to calculate radiative transfer.
!@+   Update labile C pool (kg C/tree), Weng 2017-4-5, Weng 2018-4-16
      use ent_const
      use ent_types
      use patches, only : patch_print
      use physutil, only:  QSAT

      implicit none
      real*8, intent(in) :: dtsec
      type(patch),pointer :: pp
      !----Local----------------!
      type(cohort),pointer :: cop
      type(psdrvtype) :: psdrvpar !Met biophysics drivers, except for radiation.
      real*8 :: ci_umol !umol/mol, Leaf internal CO2
      real*8 :: ca_umol !umol/mol, Ambient air CO2
      real*8 :: TsurfK, TcanK, TsoilK, Pa, RH !,rh
      real*8 :: CosZen !,betad
      real*8 :: IPAR   !Incident PAR 400-700 nm (W m-2)
      real*8 :: fdir   !Fraction of IPAR that is direct
      real*8 :: Gb     !Boundary layer conductance of water vapor(mol m-2 s-1)
      real*8 :: fdry_pft_eff ! pft-specific effective dry canoopy fraction
      real*8 :: Anet,Atot,Rd    !umol m-2 s-1
      real*8 :: Iemis ! Nadine's isoprene emission umol m-2 s-1
      real*8 :: Ciavg,GCANOPY,TRANS_SW ! Ci,NPP !,R_auto
      real*8 :: molconc_to_umol
      type(canraddrv) :: cradpar
      ! Weng, 05-04-2017
      real*8 :: IPARlayer(10)   ! Incident PAR 400-700 nm (W m-2) for each layer
      real*8 :: TRANS_SW_layer ! PAR for each cohort
      real*8 :: crownLAI,actualLAI,kappa  ! Accumulative LAI above the lower layer
      type(canraddrv) :: layerPar
      integer :: i

      kappa = 0.5d0 ! just for test
      ! Check if bare soil
      if ( .NOT.ASSOCIATED(pp%tallest)) then ! bare soil
        pp%TRANS_SW = 1.d0
        return
      endif

      if (( pp%tallest%pft.eq.0).or.(pp%tallest%pft > N_PFT)) then
        print *,"photosynth_cond: wrong pft = ", pp%tallest%pft
        call patch_print(6,pp,"ERROR ")
        call stop_model("photosynth_cond: wrong pft",255)
      endif

      !* ZERO SOME OUTPUT VARIABLES AT PATCH LEVEL
      pp%TRANS_SW = 1.d0 !Case of zero LAI.
      Iemis = undef !NKDEBUG
      !* Time-stepped outputs:  CNC, Ci, Qf.


      !* Radiation drivers *!
      IPAR = pp%cellptr%IPARdir + pp%cellptr%IPARdif
      if (pp%cellptr%IPARdir.eq.0.d0) then
        fdir = 0.d0
      else
        fdir = pp%cellptr%IPARdir / IPAR
      endif
      CosZen = pp%cellptr%CosZen
      Pa = pp%cellptr%P_mbar * 100.d0

      !* Other photosynthesis drivers *!
      !Set up psdrvpar - pack in met drivers.
      Gb = pp%cellptr%Ch*pp%cellptr%U* Pa/
     &     (gasc*(pp%cellptr%TairC+KELVIN)) !m/s * N/m2 * mol K/J * 1/K = mol/m2/s

      molconc_to_umol = gasc * (pp%cellptr%TcanopyC + KELVIN)/Pa * 1d6
      ca_umol = 370.0d0 ! pp%cellptr%Ca * molconc_to_umol  !Convert to umol/mol or ppm.
      ci_umol = pp%cellptr%Ci * molconc_to_umol  !This is solved for is pscubic in FBBphotosynthesis.f.  Replace with dummy initialization.
      TcanK = pp%cellptr%TcanopyC + KELVIN
      TsurfK = pp%cellptr%TairC + KELVIN
      TsoilK = pp%cellptr%Soiltemp(1) + KELVIN
      RH = min(1.d0,
     &     max(pp%cellptr%Qf,0.d0)/Qsat(TcanK,
     &     2500800.d0 - 2360.d0*(TsurfK-KELVIN),Pa/100.d0))

      ! put met forcing data to psdrvpar
      call biophysdrv_setup(ca_umol,ci_umol,
     &     pp%cellptr%TcanopyC, Pa, RH,
     &     psdrvpar)  !Equation for latent heat of vaporizaiton comes from ..?

      ! IPAR for each layer of crown, Weng
      IPARlayer(1) = IPAR ! top layer
      do i=1, pp%crownlayers
         ! Beer's law, just for test
         !IPARlayer(i+1) = IPARlayer(i) * exp(-kappa*pp%LAIlayer(i))

         ! use radiative transfer scheme of Spitters model
         actualLAI = pp%LAIlayer(i) /(1.0d0-fgap)
         !Set up cradpar and calculate diffuse and direct PAR.
         call canopy_rad_setup(pp%tallest%pft,
     &     CosZen,fdir,IPAR,
     &     actualLAI,pp%albedo(1),layerPar)
         call canopy_transmittance(TRANS_SW,CosZen,fdir,layerPar)
         IPARlayer(i+1) = (TRANS_SW*(1.0d0-fgap) + fgap) * IPARlayer(i) ! assume fgap PAR go to next layer
      enddo

      !* LOOP THROUGH COHORTS *!
      Ciavg = 0.d0
      cop => pp%tallest
      do while (ASSOCIATED(cop))
        !* Assign vegpar
        if(cop%LAI.gt.0.d0) then
          ! fraction of wet leaves
          if (BMEspdata(cop%pft)%leaftype.eq.BROADLEAF) then
            ! stomata on underside of leaves so max stomatal blocking = 0
            fdry_pft_eff = 1.d0
          elseif(BMEspdata(cop%pft)%leaftype.eq.NEEDLELEAF) then
            ! Bosveld & Bouten (2003) max stomatal blocking = 1/3
            fdry_pft_eff = 1.d0 - min(pp%cellptr%fwet_canopy, 0.333d0)
          else
             fdry_pft_eff = 1.d0
          endif

          !KIM - water_stress4 uses Soilmoist as a satured fraction
          cop%stressH2O = water_stress4(cop%pft, N_DEPTH,
     i          pp%cellptr%Soilmoist(:),
     &          pp%tallest%fracroot, pp%cellptr%fice(:),
     &          cop%stressH2Ol(:))

          call calc_Pspar(dtsec,cop%pft
     i         ,psdrvpar%Pa,psdrvpar%Tc
     i         ,O2frac*psdrvpar%Pa
     i         ,cop%stressH2O,cop%Sacclim,cop%llspan)

!         just for test: Atot of the first ch >> the others
          if(cop%crownarea>0.d0)then
               crownLAI = cop%leafarea/cop%crownarea
          else
               crownLAI = 0.0 ! cop%leafarea*cop%n ! cop%LAI !/cop%n/cop%crownarea
          endif

          call canopyfluxes(dtsec, pp%tallest%pft ! cop%pft ! cohortfluxesPPA
     &         ,pp%albedo(1)
     &         ,crownLAI,crownLAI
     &         ,IPARlayer(cop%layer) * 4.05d0  ! IPAR*4.05d0
     &         ,CosZen,fdir
     &         ,Gb
     &         ,psdrvpar
     &         ,GCANOPY,Anet,Atot,Rd !NOTE: Ci should be cohort level
     &         ,Iemis
     &         ,TRANS_SW)       !NOTE:  Should include stressH2O.

         !* Assign outputs to cohort *!
         !* Account for wet leaves with pp%cellptr%fwet_canopy.
          cop%Ci = psdrvpar%ci * !ci is in mole fraction
     &         psdrvpar%Pa/(gasc * (psdrvpar%Tc+KELVIN)) !mol m-3
          !GCANOPY = GCANOPY * cop%stressH2O ! added by Weng for test: no water stress on Atot, but keep it for transpiration
          cop%GCANOPY = GCANOPY*fdry_pft_eff*(gasc*TsurfK)/Pa  !Convert mol-H2O/m2/s to m/s
          cop%GPP = Atot * fdry_pft_eff * 0.012d-6*cop%crownarea ! !umol/m2/s to kg-C/tree/s
          cop%IPP = Iemis * 0.0600d-6*cop%crownarea  !umol m-2 s-1 to kgC-isoprene/tree/s
       else                     !Zero LAI or no light
          GCANOPY = 0.d0
          Anet = 0.d0
          Atot = 0.d0
          TRANS_SW = 1.d0
          cop%GCANOPY=0.d0 !May want minimum conductance for stems.
          cop%Ci = EPS
          cop%GPP = 0.d0
          cop%IPP = 0.d0
          !!! why not just zero ??? -IA
          Rd = 0.d0 ! gamma_r * pspar%Vcmax * cop%LAI ! Weng ?
       endif

       !* Update cohort respiration components, NPP, C_lab

       ! Added by Weng, 2017-02-08
       cop%leafRd = Max(0.d0,
     &      MIN(0.04*cop%C_lab,Rd*cop%crownarea)) ! per tree
       call BiomeE_PlantResp(cop, dtsec,
     &      TcanK,TsoilK,
     &      pp%cellptr%airtemp_10d+KELVIN,
     &      pp%cellptr%soiltemp_10d+KELVIN)

       !! update C_lab here, for test, Weng, 04-03-2017
        cop%C_lab = cop%C_lab + cop%NPP * dtsec ! /cop%n ! Weng, 5-5-2017
       ! The following is to disconnect plant grwoth from physiology, Weng 2017-04-03
       !  cop%C_lab = cop%C_lab + 0.6d-5 * cop%crownarea

       ! for daily GPP and NPP per tree, Weng, 2017-02-07
        cop%GPP_daily =cop%GPP_daily+cop%GPP*dtsec
        cop%NPP_daily =cop%NPP_daily+cop%NPP*dtsec
        cop%Ra_daily  =cop%Ra_daily +cop%R_auto*dtsec

        Ciavg = Ciavg + cop%Ci*cop%LAI
        !set values for debugging
        ent_d%Anet(cop%pft) = Anet * 0.012d-6
        ent_d%Atot(cop%pft) = Atot * 0.012d-6
        ent_d%Rd(cop%pft) = Rd * 0.012d-6
        ent_d%GCANOPY(cop%pft) = GCANOPY
        ent_d%TRANS_SW(cop%pft) = TRANS_SW ! /(3600*24*1000.d0)

        cop => cop%shorter
      end do ! cohort by cohort

      !* Patch-level OUTPUTS *!
      if ( pp%LAI > 0.d0 ) then
        pp%Ci = Ciavg/pp%LAI
      else
        pp%Ci = 0.d0
      endif

      ! for patch level albedo and TRANS_SW in Ent when there is only one cohort, Weng
      call canopy_rad_setup(pp%tallest%pft,CosZen,fdir,IPAR, !!pft not needed but need cradpar
     &     pp%LAI,pp%albedo(1),cradpar)
      call canopy_transmittance(TRANS_SW,CosZen,fdir,cradpar)
      pp%TRANS_SW = TRANS_SW

      end subroutine BiomeE_CanopyModel


!============================================================
! Weng, 2018-11-03:  It's the begining of plant growth for BiomeE model
! Weng, 5/11/2017: all the fluxes are per tree (defined by input)
! kgC/s (per tree)
! Never do unit converting in this subroutine!!!!!
      subroutine BiomeE_PlantResp(cop,dtsec,
     &     TcanopyK,TsoilK,TairK_10d,TsoilK_10d)

!@sum Updates plant autotrophic respiration, NPP, and C_lab
!@sum (All in kg-C/tree/s)
      implicit none
      type(cohort),intent(inout) :: cop
      real*8,intent(in) :: dtsec
      real*8,intent(in) :: TcanopyK
      real*8,intent(in) :: TsoilK
      real*8,intent(in) :: TairK_10d
      real*8,intent(in) :: TsoilK_10d
      !----Local-----
      integer :: pft
      real*8 :: Resp_fol, Resp_sw, Resp_froot
      real*8 :: Resp_maint, Resp_can, Resp_growth
      real*8 :: totLFR,r_NSC
      real*8 :: rootCN
      real*8 :: facclim,K_cambium,Acambium
      real*8 :: umols2kgC

      pft = cop%pft
      associate(sp =>BMEspdata(cop%pft))

      umols2kgC = 0.012D-6 ! * ndens !Convert from umol/s to kgC/s (per tree).
      facclim = frost_hardiness(cop%Sacclim)
      rootCN = 30.d0 ! BMEspdata(pft)%LMA*1000.d0/BMEspdata(pft)%Nleaf

      !* Maintenance respiration - leaf + sapwood + storage

      ! Sapwood respiration as a function of cambium area, Weng 2017-12-14
      K_cambium = sp%K_cambium
      Acambium  = sp%f_Acam * PI*cop%DBH*cop%h ! cambium area
      totLFR    = cop%C_fol + cop%C_froot + 0.005d0
      r_NSC     = Max(MIN((cop%C_lab-0.2*cop%NSCstar)/totLFR,
     &            1.d0),0.d0)

      !!! Rd is already per unit area ! converted per tree, Weng 5-5-2017
      Resp_fol = 1.2d-8 * r_NSC * cop%leafRd  !Rd is already per tree, convert umol to kgC ! Weng, 2020-02-03
      Resp_sw = umols2kgC * r_NSC *
     &    Resp_sapwood_maint(pft,K_cambium,Acambium,
     &                   TcanopyK,TairK_10d,facclim)

      ! fine root respiration
      Resp_froot = umols2kgC * r_NSC *
     &   Resp_cpool_maint(pft,cop%C_froot,rootCN,
     &            TsoilK,TsoilK_10d,facclim)

      end associate

      !Total maintainence respiration
      Resp_maint = Resp_fol + Resp_sw + Resp_froot

      ! Update Rauto and NPP
      Resp_growth = cop%R_growth !kgC s-1 tree-1. Weng 05-05-2017

      !* Growth respiration tied to GPP; compensates for tissue growth respir.
      ! Turned off by Weng 2021-11-17: R_growth = max(0.d0, growth_r*(Acan - Rmaint) - Rtgrowth)
      !Resp_can = Resp_can_growth(pft,cop%GPP,Resp_maint,Resp_growth)
      Resp_can = 0.d0 ! Weng, 11/17/2021
      cop%R_root = Resp_froot
      cop%R_auto = Resp_maint + Resp_growth + Resp_can

      ! Detach Rauto from plant physiology, test only!! ! Weng, 2021-11-17
      cop%R_auto = 0.45d0 * cop%GPP

      ! Tree NPP
      cop%NPP    = cop%GPP - cop%R_auto !kg-C/tree/s

      end subroutine BiomeE_PlantResp

!---------------------------------------------------------------------------
      subroutine canopyfluxes(dt, pft,
     i     canalbedo,LAIcanopy,LAIcohort,
     i     IPAR,CosZen,fdir,Gb,psdrvpar,
     o     Gs,Anet,Atot,Rd,Iemis,TRANS_SW)
!@sum canopyfluxes Calculates cohort photosynthesis and conductance with
!@sum Farqhuar et al. (1980) photosynthesis, Ball-Berry stomatal conductance,
!@um  and Spitters (1986, 1987) canopy radiation (sunlit, shaded leaves).
!@sum Integrates vertically over the canopy with Simpson's Rule.
!@sum Ci is updated at the canopy level using the canopy boundary layer
!@sum conductance as in Friend and Kiang (2005).
!@+

!@sum Spitters (1986) canopy radiative transfer (sunlit/shaded leaves), and
!@+   Simpson's Rule (Price, Numerical Recipes) for canopy layering.
!@+   Call photosynthesis/conductance routines from other module to scale
!@+   up leaf-level fluxes to the canopy.

!@+   If PAR is not directly available, the following conversions may be used:
      !  From total shortwave (W m-2) to PAR (umol m-2 s-1) (Monteith & Unsworth):
      !          PAR(umol m-2 s-1) = 2.3(umol/J)*SW(W m-2)
      !  From PAR (W m-2) to PAR (umol m-2 s-1) (U.Maryland, Dept. of Met., PAR Project),
      !  suggest nominal 485 nm for conversion, which gives:
      !          PAR(umol m-2 s-1) = 4.05(umol/J) * PAR(W m-2)
      !  Dye, D.G. (2004) suggests a slightly different conversion:
      !          PAR(umol m-2 s-1) = 4.56(umol/J) * PAR(W m-2)

      implicit none
      real*8,intent(in) :: dt   !time step (seconds)
      integer,intent(in) :: pft !Plant functional type, 1-C3 grassland
      real*8,intent(in) :: canalbedo !Whole-canopy albedo
      real*8,intent(in) :: LAIcanopy  !Whole-canopy LAI
      real*8,intent(in) :: LAIcohort  !cohort LAI
      real*8,intent(in) :: IPAR !Incident PAR (umol m-2 s-1)

      real*8,intent(in) :: CosZen !Cosine ofSolar zenith angle
      real*8,intent(in) :: fdir !Fraction of PAR that is direct
      real*8,intent(in) :: Gb   !Canopy boundary layer conductance of water vapor (mol m-2 s-1)
      type(psdrvtype) :: psdrvpar !Photosynthesis drivers, except for radiation.
      real*8,intent(inout) :: Gs !Canopy stomatal conductance of water vapor (mol m-2 s-1)
      real*8,intent(out) :: Anet !Leaf net photosynthesis (micromol m-2 s-1)
      real*8,intent(out) :: Atot !Leaf gross photosynthesis (micromol m-2 s-1)
      real*8,intent(out) :: Rd  !Leaf respiration (umol/m2/s)
      real*8,intent(out) :: TRANS_SW !Transmittance of shortwave to ground surface.
      real*8,intent(out) :: Iemis ! Leaf isoprene emission Nadine (micromol m-2 s-1)
!      integer,intent(in) :: if_ci !FLAG 0-don't calculate ci, 1-calculate ci

      !----Local var-------
      type(canraddrv) :: cradpar
      real*8 :: Gsint    !Gs canopy conductance from qsimp integration (mol/m2/s)
      real*8 :: Rdint    !Rd canopy foliage respiration form qsimp integration (umol/m2/s)
      real*8 :: Iint     !Iint Nadine's isoprene emission
      !real*8 :: Ciconc          !Leaf internal CO2 mole fraction calculate (umol mol-1)

      !Set up cradpar and calculate diffuse and direct PAR.
      call canopy_rad_setup(pft,CosZen,fdir,IPAR,
     &     LAIcanopy,canalbedo,cradpar)

      !Calculate net photosynthesis, integrate vertically with Simpson's Rule.
!      call qsimp(cradpar%LAI,cradpar,psdrvpar,Gb,Atot,Gsint,Rdint,Iint)
      call qsimp(LAIcohort,cradpar,psdrvpar,Gb,Atot,Gsint,Rdint,Iint)

      call canopy_transmittance(TRANS_SW,CosZen,fdir,cradpar)

      !Calculate conductance and ci at the canopy level
!      call Gs_bound(dt, LAI,Gsint, Gs) !Limit rate of change of Gs.
!      Ciconc = calc_Ci_canopy(ca,Gb,Gs,Anet,LAI,IPAR)

      !Different options to solve for Gs and Ci:
      !If Ciconc < CiMIN, then Gs was too big, so limit Gs to allow
      !only Ciconc >= CiMIN
!      if (Ciconc < CiMIN) then
!        Ciconc = CiMIN
!        Gs =  Gs_from_Ci(Anet, ca, Gb, Gs, ci,IPAR)
!      else
!        ci = Ciconc
!      endif

      Rd = Rdint
      Gs = Gsint
      Anet = Atot - Rd
      Iemis = Iint

      end subroutine canopyfluxes

!################## PHOTOSYNTHESIS #########################################

      subroutine photosynth_sunshd(
     &      Lcum,crp,psd,Gb,Aleaf,gsleaf,Rdleaf,Ileaf)
!@sum photosynth_sunshd  Calculates sunlit and shaded leaf fluxes.
!@auth N.Y.Kiang with Spitters parameters
      implicit none
      real*8,intent(in) :: Lcum    !Cumulative LAI from top of canopy (m2/m2)
      type(canraddrv) :: crp  !Canopy radiation parameters
      type(psdrvtype) :: psd   !Photosynthesis met drivers
      real*8,intent(in) :: Gb   !Leaf boundary layer conductance (mol/m2/s)
      real*8,intent(out) :: Aleaf !Leaf Net assimilation of CO2 in layer (umol m-2 s-1)
      real*8,intent(out) :: gsleaf !Leaf Conductance of water vapor in layer (mol m-2 s-1)
      real*8,intent(out) :: Rdleaf !Leaf respiration (umol m-2 s-1)
      real*8,intent(out) :: Ileaf ! Leaf isoprene emission (umol m-2 s-1)
      !------Local---------
      real*8 fsl  !Fraction of sunlit foliage in layer at Lc (unitless).
      real*8 Isl !PAR absorbed by sunlit foliage at Lc (umol/m[foliage]2/s).
      real*8 Ish !PAR absorbed by shaded foliage at Lc (umol/m[foliage]2/s).
      real*8 Asl !Anet from sunlit leaves (umol/m2/s)
      real*8 Ash !Anet from shaded leaves (umol/m2/s)
      real*8 gssl,gssh !Leaf conductance, sunlit,shaded (mol-H2O/m2/s)
      real*8 Rdsl, Rdsh
      real*8 Iemisl  ! Nadine's isop emis from sunlit leaves (umol/m2/s)
      real*8 Iemiss  ! Nadine's isop emis from shaded leaves (umol/m2/s)
      integer :: sunlitshaded !1-sunlit, 2-shaded

      call canopy_rad(Lcum,crp,Isl,Ish,fsl)

      !Calculate photosynthesis and stomatal conductance.
      sunlitshaded = 1
      call pscondleaf(crp%pft,Isl,psd,Gb,gssl,Asl,Rdsl,sunlitshaded,
     & Iemisl)
!      write(992,*) 'shaded'
      sunlitshaded = 2
      call pscondleaf(crp%pft,Ish,psd,Gb,gssh,Ash,Rdsh,sunlitshaded,
     & Iemiss)

      Aleaf  = fsl*Asl    + (1.0d0 - fsl)*Ash
      gsleaf = fsl*gssl   + (1.0d0 - fsl)*gssh
      Rdleaf = fsl*Rdsl   + (1.0d0 - fsl)*Rdsh
      Ileaf  = fsl*Iemisl + (1.0d0 - fsl)*Iemiss

      end subroutine photosynth_sunshd

!################# CANOPY CONDUCTANCE ######################################

      subroutine Gs_bound(dt, LAI, Gsnew, Gsinout)
!@sum Gs_bound Limit rate of change of Gs (umol m-2 s-1)
!@+   Call to this was commented out - Igor?
!@auth A.Friend
      ! Required change in canopy conductance to reach equilibrium (m/s).
      implicit none
      real*8,intent(in) :: dt, LAI
      real*8,intent(in) :: Gsnew !Canopy conductance of curr time step (mol m-2 s-1)
      real*8,intent(inout) :: Gsinout !Bounded canopy conductance (mol m-2 s-1)
      !---Local----!
      real*8 :: Gsold           !Canopy conductance at prev time step (mol m-2 s-1)
      real*8, parameter :: rhoH2O = 998.2 !Density of water (1000 kg/m^3)
      real*8, parameter :: MW_H2O = 18.015 !Molecular weight of water (g/mol)
      !real*8, parameter :: ghi = 0.006d0*rhoh2o*1000./MW_H2O !Upper limit of gs leaf (mol m-2 s-1)
      !real*8, parameter :: glo = 0.000001d0*rhoH2O*1000./MW_H2O !Lower limit of gs leaf (mol m-2 s-1), See Ball and Berry paper.
      real*8, parameter :: ghi = 333.0 !Conversion from 6 mm s-1 upper limit.(mol m-2 s-1)
      real*8, parameter :: glo = .015 !Temperate grassland. Korner (1994) (mol m-2 s-1)

      real*8 :: dGs, dGs_max

      dGs=Gsnew-Gsinout
      Gsold = Gsinout
      Gsinout = Gsnew
      !nu Limit Gs change over timestep because of guard cell mechanics (m/s)
      dGs_max=dt*LAI*(ghi-glo)/1800.0D0
      if( dGs.gt.dGs_max) Gsinout = Gsold + dGs_max
      if(-dGs.gt.dGs_max) Gsinout = Gsold - dGs_max
      ! Biological limits of absolute Gs (m/s).
      if(Gsinout.gt.ghi*LAI) Gsinout=ghi*LAI
      if(Gsinout.lt.glo*LAI) Gsinout=glo*LAI

      end subroutine Gs_bound

!---------------------------------------------------------------------------
      function Gs_from_Ci(Anet,Ca,Gb,Gs,Ci,IPAR) Result(gsout)
!@sum Inversion of calc_Ci_canopy
      implicit none
      real*8 :: Anet !Leaf net assimilation of CO2 (umol m-2 s-1)
      real*8 :: ca !Ambient air CO2 mole fraction at surface reference height (umol mol-1)
      real*8 :: gb !Canopy boundary layer conductance of water vapor (mol m-2 s-1)
      real*8 :: gs !Canopy Stomatal conductance of water vapor(mol m-2 s-1)
      real*8 :: ci !Leaf internal CO2 concentration (umol mol-1)
      real*8 :: IPAR !Incident PAR (umol m-2 s-1)
      real*8 :: gsout !Adjusted Gs (mol m-2 s-1)
      !----- Local ------

      if (IPAR.lt.IPARMINMOL) then  !Stomates closed
        gsout = gs
      else
        gsout =  1.65/((Ca - Ci)/(Anet) - 1.37/gb )
      endif

      end function Gs_from_Ci

!################# RADIATIVE TRANSFER ######################################
      subroutine canopy_rad_setup(
     i     pft, CosZen, fdir,IPAR,LAI,canalbedo,
     o     crp)
!@sum Calculate diffuse and direct PAR, canopy albedo, if unknown, from given
!@sum solar zenith angle, and other canopy radiation parameters.

      implicit none

! NOTES FOR OBTAINING PAR
! If only incident shortwave is available, convert
!  total shortwave (W/m2) to PAR (umol/m2/s).
!  2.3 umol/J = shortwave to fraction that is PAR (Monteith & Unsworth).
!      PAR=2.3d0*max(srht,zero)/(1.0D0-0.08D0)
! Current: Replaced back-calculation with actual incident PAR at surface.

! GENERAL PAR CONVERSIONS.
! nyk  For 400-700 nm, energy is 3.3-6 umol photons/J.
!      Top-of-atmosphere flux-weighted average of energy is 4.54 umol/J.
!      Univ. Maryland, Dept. of Met., PAR Project, suggests nominal
!      485 nm for conversion, which gives 4.05 umol/J.

      !Input parameters
      integer,intent(in) :: pft
!      real*8,intent(in) :: solarzen !Solar zenith angle
      real*8,intent(in) :: CosZen !Cosine of Solar zenith angle
      real*8,intent(in) :: fdir   !Fraction of shortwave radiation that is direct
      real*8,intent(in) :: IPAR    !PAR incident at top of canopy (umol/m2/s)
      real*8,intent(in) :: LAI    !Leaf area index of whole canopy
      real*8,intent(in) :: canalbedo !Canopy albedo, if known

      !Output parameters
      type(canraddrv) :: crp
      !real*8,intent(out) :: I0df, I0dr !Diffuse and direct PAR incident at top of canopy (umol m-2 s-1)
      !----Local var-------
      real*8 :: sbeta  !Sine of solar zenith angle
      ! CONSTANTS
      real*8,parameter :: sigma=0.2d0  !Leaf scattering coefficient
      real*8,parameter :: kdf=0.71d0  !Canopy extinction coefficient
      real*8,parameter :: EPS=1.d-8   !Small, to prevent div zero.

      sbeta = CosZen ! cos(solarzen) = sin(solarelev)

!      crp%solarzen = acos(CosZen)
      crp%CosZen = CosZen
! Direct beam PAR incident on canopy (umol/m2/s).
      crp%I0dr=fdir*IPAR
! Diffuse PAR incident on canopy (umol/m2/s).
      crp%I0df=(1.0D0-fdir)*IPAR
! Canopy reflectivity depends on sun angle. Or take prescribed albedos.
      crp%sigma = sigma
      crp%kdf = kdf
      crp%sqrtexpr=sqrt(1.0D0-crp%sigma)
! Canopy extinction coefficient for direct beam radiation depends on
! sun angle (unitless).
      crp%kbl=0.5D0*crp%kdf/(0.8D0*crp%sqrtexpr*sbeta+EPS)

      !Hack canopy reflectivity.  This assumes a closed canopy!
      if (canalbedo.eq.0.0) then !Use if radiation gets initialized late.
        crp%rhor=((1.0D0-crp%sqrtexpr)/(1.0D0+crp%sqrtexpr))*
     &       (2.0D0/(1.0D0+1.6D0*sbeta))
        crp%canalbedo = crp%rhor
      else
       !Prescribed veg albedos until have canopy scheme.
        crp%rhor = canalbedo
        crp%canalbedo = crp%rhor
      end if

      crp%pft = pft
      crp%LAI = LAI !canopy LAI

      end subroutine canopy_rad_setup

!-----------------------------------------------------------------------------
      subroutine canopy_rad(
!@sum Calculate PAR on sunlit and shaded foliage in layer Lc
     i     Lc,                  !Cumulative LAI at layer from top of canopy
     i     crp,                !veg parameters
     o     Isla,                !PAR on shaded foliage in layer Lc
     o     Isha,                 !PAR on sunlit foliage in layer Lic
     o     fsl)                 !Fraction of sunlit foliage at Lc (unitless).

      implicit none
! Get incident radiation on layer, diffuse/direct, sunlit/shaded.
      real*8,intent(in) :: Lc
      type(canraddrv) :: crp
!@var Isla PAR absorbed by sunlit foliage at Lc (umol/m[foliage]2/s).
      real*8,intent(out) :: Isla
!@var Isha PAR absorbed by shaded foliage at Lc (umol/m[foliage]2/s).
      real*8,intent(out) :: Isha
!@var fsl Fraction of sunlit foliage in layer at Lc (unitless).
      real*8,intent(out) :: fsl

      !-------Local vars--------------------------------------------------
      real*8 :: sbeta           !cos(solarzen) = sin(solarelev)
      real*8 :: I0df, I0dr !Diffuse and direct IPAR (umol m-2 s-1)
!@var Idfa Diffuse PAR in canopy at Lc (umol/m[ground]2/s).
      real*8 Idfa
!@var Idra Direct PAR in canopy at Lc (umol/m[ground]2/s).
      real*8 Idra
!@var Idrdra Direct PAR at Lc that remains direct (umol/m[ground]2/s).
      real*8 Idrdra
      real*8,parameter :: zero = 0.d0
      real*8 :: Ish, Isl !Incident rather than absorbed PAR (umol/m2/s).

      sbeta = crp%CosZen !cos(solarzen) = sin(solarelev)
      I0df = crp%I0df
      I0dr = crp%I0dr
! When sun is above the horizon.
      if(sbeta.gt.EPS)then
!? Diffuse PAR in canopy at Lc (umol/m[ground]2/s). ?Old comment.
! Diffuse PAR that is absorbed at Lc (umol/m2-gnd/s). Eq. 10.
        Idfa=(1.0D0-crp%rhor)*I0df*crp%kdf *exp(-crp%kdf*Lc)
!? Direct PAR in canopy at Lc (umol/m[ground]2/s). ?Old comments.
! Direct PAR that is absorbed at Lc (umol/m2-gnd/s). Eq. 11.
        Idra=(1.0D0-crp%rhor)*I0dr* crp%sqrtexpr * crp%kbl *
     $        exp(-crp%sqrtexpr*crp%kbl*Lc)
!? Direct PAR at Lc that remains direct (umol/m[ground]2/s).?Old comment
! Non-scattered part of direct flux that is absorbed (umol/m2-gnd/s). Eq. 12.
        Idrdra=(1.0D0-crp%sigma)*I0dr * crp%kbl *
     $        exp(-crp%sqrtexpr*crp%kbl*Lc)
! PAR absorbed by shaded foliage at Lc (umol/m[foliage]2/s). Diffuse + Direct portion that is scattered (becomes diffuse)
        Isha=Idfa+(Idra-Idrdra)
! PAR aborbed by sunlit foliage at Lc (umol/m[foliage]2/s)._
        Isla=Isha+(1.0D0-crp%sigma)* crp%kbL * I0dr

!######## HACK - above are absorbed radiation, not incident ###########
!######## Replace with Spitters' equations for incident radiation #####
!        Isl = (1.0D0-crp%rhor)*I0dr*exp(-crp%kbL*crp%sqrtexpr*Lc)
!        Ish = (1.0D0-crp%rhor)*I0df*exp(-crp%kdf*Lc)
!        Isla = I0dr + Ish  !## HACK - this is to pass out Isl, using Isla var.
!        Isha = Ish !## HACK - this is to pass out Ish, using Isha var.
! Fraction of sunlit foliage in layer at Lc (unitless).
        fsl=exp(-crp%kbL*Lc)
        if(Isha.le.0.0D0)Isha=0.0D0
        if(Isla.le.0.0D0)Isla=0.0D0
        if(fsl.lt.zero)fsl=zero
      else
        Isha=0.0D0
        Isla=0.0D0
        fsl=zero
      endif

      end subroutine canopy_rad

!-----------------------------------------------------------------------
      subroutine canopy_transmittance(TRANS_SW, sbeta,fdir,crp)
!@sum Calculate the transmittance (fraction) of radiation through the canopy
!@+   to the soil surface. Verbatim from Spitters et al. (1986,1987)
!     Transmission of shortwave radiation through canopy. This has errors.
!----------------------------------------------------------------------!
      implicit none
!----------------------------------------------------------------------!
!Parameters -------------------------------------------------------------
!@var trans_sw Total canopy transmittance of shortwave (fraction)
      real*8, intent(out) :: TRANS_SW
!@var sbeta Cos of solar zenith angle
      real*8, intent(in) :: sbeta
!@var fdir Fraction of surface visible radiation that is direct
      real*8, intent(in) :: fdir
! Canopy radiation parameters
      type(canraddrv) :: crp

!Local variables --------------------------------------------------------
! Get assigned from crp passed in.
      real*8 sigma,sqrtexpr,rhor,kbl,kdf,alai
!@var absdf Absorbance of diffuse SW radiation (fraction)
      real*8 absdf
!@var absdr Absorbance of direct SW radiation (fraction)
      real*8 absdr
!@var absdrdr Absorbance of direct SW that remains direct (fraction)
      real*8 absdrdr
!@var abssh Absorbance of SW through shaded foliage (fraction)
      real*8 abssh
!@var abssl Absorbance of SW through sunlit foliage (fraction)
      real*8 abssl
!@var fracsl Fraction of leaves in canopy that are sunlit (fraction)
      real*8 fracsl
!------------------------------------------------------------------------
      sigma=crp%sigma
      sqrtexpr=crp%sqrtexpr
      rhor=crp%rhor
      kbl=crp%kbl
      kdf=crp%kdf
      alai=crp%LAI

! Diffuse shortwave in canopy (umol/m[ground]2/s).
      if(sbeta.gt.zero)then
        absdf=(1-fdir)*(1.0D0-rhor)*kdf*exp(-kdf*alai)
! Direct shortwave in canopy (umol/m[ground]2/s).
        absdr=fdir*(1.0D0-rhor)*sqrtexpr*kbl*exp(-sqrtexpr*kbl*alai)
! Direct shortwave that remains direct (umol/m[ground]2/s).
        absdrdr=absdr*(1.0D0-sigma)*kbl*exp(-sqrtexpr*kbl*alai)
! Shortwave penetrating shaded foliage (umol/m[foliage]2/s).
        abssh=absdf + (absdr-absdrdr)
! Shortwave penetrating sunlit foliage (umol/m[foliage]2/s).
        abssl=abssh+(1.0D0-sigma)*kbL*fdir
        if(abssh.lt.0.0D0) abssh=0.0D0
        if(abssl.lt.0.0D0) abssl=0.0D0
        if(abssh.gt.1.0D0) abssh=1.0D0
        if(abssl.gt.0.0D0) abssl=1.0D0
        fracsl=exp(-kbl*alai)
      else
        abssh=0.001D0
        abssl=0.0D0
        fracsl=0
      endif

!      TRANS_SW = 1.d0 - ((1.d0-fracsl)*abssh + fracsl*abssl)
      TRANS_SW = ((1.d0-fracsl)*abssh + fracsl*abssl)

!      write(910,*) TRANS_SW,alai,sigma,kdf, kbl,sqrtexpr,rhor,
!     &     abssh, abssl, sbeta


      end subroutine canopy_transmittance

!#############################################################################
!###################### INTEGRATION ##########################################

      subroutine qsimp(Xlim,crp,psp,Gb,S,Sg,Sr,Si)
!----------------------------------------------------------------------!
!@sum qsimp Numerical routine to calculate canopy photosynthesis by
!@+   increasing the number of
!@+   layers in the integration until the result (S) changes by less than
!@+    0.1 umol[CO2]/m2/s.
!@auth A.D.Friend
!----------------------------------------------------------------------!
      implicit none
!----------------------------------------------------------------------!
!Passed parameters
      real*8,intent(in) :: Xlim !cohort LAI
      type(canraddrv) :: crp    !contains canopy LAI
      type(psdrvtype) :: psp
      real*8,intent(in) :: Gb
      !real*8,intent(in) :: ca, ci,T,P,RH
      real*8,intent(out) :: S,Sg,Sr,Si
!----------------------------------------------------------------------!
!Local variables
      integer, parameter :: MAXIT=6
      real*8,parameter :: ERRLIM=0.1d0
      integer IT, layers
      real*8 :: A, B, OST,OS,ST, ERROR
      real*8 :: Ac,Bc
      real*8 :: OSTg, OSg, STg !NK
      real*8 :: OSTr, OSr, STr !NK
      real*8 :: OSTi, OSi, STi !NK + NU
!----------------------------------------------------------------------!
! Calculate canopy radiative transfer and photosynthesis in increasing
!   number of layers.
! NOTE:  ERROR is only calculated on photosynthesis value (ST).
!        Corresponding conductance value is then returned also (STg).
!        and isoprene emission value (STi)
      A=0.0D0
      B=crp%LAI  !Canopy LAI.  Ideally, qsimp should not see inside crp.
      Ac=0.0d0
      Bc=Xlim  !Cohort LAI. Assumes spans canopy height. Later:make different heights.
      OST=-1.D30
      OS= -1.D30
      OSTg=-1.D30 !NK
      OSg= -1.D30 !NK
      OSTr=-1.D30 !NK
      OSr= -1.D30 !NK
      OSTi= -1.D30 !NK + NU
      OSi= -1.D30 !NK + NU
      layers=1
      do 11 IT=1,MAXIT
         CALL TRAPZD(A,B,Ac,Bc,ST,STg,STr,STi,IT,layers,crp,psp,Gb)
         S=(4.D0*ST-OST)/3.D0
         Sg=(4.D0*STg-OSTg)/3.D0 !NK
         Sr=(4.D0*STr-OSTr)/3.D0 !NK
         Si=(4.D0*STi-OSTi)/3.D0 !NK+NU
         ERROR=ABS(S-OS)
         IF (ERROR.lt.ERRLIM) RETURN
         OS=S
         OST=ST
         OSg=Sg !NK
         OSTg=STg !NK
         OSr=Sr !NK
         OSTr=STr !NK
         OSi=Si !NK + NU
         OSTi=STi !NK + NU

         if(IT.gt.1) layers=layers*2
   11 enddo
      return
      end subroutine qsimp

!======================================================================
      subroutine trapzd(L1,L2,L1c,L2c,S,Sg,Sr,Si,N,layers,crp,psp,Gb)
!----------------------------------------------------------------------!
!@sum Integrates canopy photosynthesis over canopy layers using Simpson's
!@+   Rule (Press et al., 19??).
!----------------------------------------------------------------------!

      implicit none
!----------------------------------------------------------------------!
      real*8 L1,L2 !LAI of whole canopy
      real*8 L1c,L2c,S,Sg,Sr,Si !LAI and sums of cohort
      integer N,layers
      type(canraddrv) :: crp
      type(psdrvtype) :: psp
      real*8,intent(in) :: Gb
      !-----Local------------------
      integer L
      real*8 DEL,X,SUM,SUMg,SUMr,SUMi,RCL
      real*8 DELc
      real*8 Aleaf1 !Mean net photosynthesis at layer bottom (umol[CO2]/m2/s).
      real*8 Aleaf2 !Mean net photosynthesis at layer top (umol[CO2]/m2/s).
      real*8 gleaf1, gleaf2 !Mean conductance at L1, L2 (umol[H2O]/m2/s).
      real*8 Rdleaf1,Rdleaf2 !Mean leaf respiration at L1, L2 (umol[CO2]/m2/s)
      real*8 Ileaf1,Ileaf2 ! mean isop emis at L1,L2 (umol[C]/m2/s)

      if(N.eq.1)then
         call photosynth_sunshd(L1,crp,psp,Gb,Aleaf1,gleaf1,Rdleaf1,
     * Ileaf1)
         call photosynth_sunshd(L2,crp,psp,Gb,Aleaf2,gleaf2,Rdleaf2,
     * Ileaf2)
!         S=0.5D0*(L2-L1)*(Aleaf1+Aleaf2)
!         Sg=0.5d0*(L2-L1)*(gleaf1+gleaf2)
!         Sr=0.5d0*(L2-L1)*(Rdleaf1+Rdleaf2)
         S=0.5D0*(L2c-L1c)*(Aleaf1+Aleaf2) !Cohort LAI for leaf processes.
         Sg=0.5d0*(L2c-L1c)*(gleaf1+gleaf2)
         Sr=0.5d0*(L2c-L1c)*(Rdleaf1+Rdleaf2)
         Si=0.5d0*(L2c-L1c)*(Ileaf1+Ileaf2)
      else
         RCL=layers  ! convert to real*8
         DEL=(L2-L1)/RCL !Canopy LAI for radiative transfer.
         DELc=(L2c-L1c)/RCL !Cohort LAI for radiative transfer.
         X=L1+0.5D0*DEL
         SUM=0.D0
         SUMg=0.d0
         SUMr=0.d0
         SUMi=0.d0
         do 11 L=1,layers
           call photosynth_sunshd(X,crp,psp,Gb,Aleaf1,gleaf1,Rdleaf1,
     *  Ileaf1)
           SUM = SUM + Aleaf1
           SUMg = SUMg + gleaf1
           SUMr = SUMr + Rdleaf1
           SUMi = SUMi + Ileaf1
           X=X+DEL
   11    continue
!         S=0.5D0*(S+(L2c-L1c)*SUM/RCL)
!         Sg=0.5D0*(Sg+(L2c-L1c)*SUMg/RCL)
!         Sr=0.5D0*(Sr+(L2c-L1c)*SUMr/RCL)
         S=0.5D0*(S+(L2-L1)*SUM/RCL)
         Sg=0.5D0*(Sg+(L2-L1)*SUMg/RCL)
         Sr=0.5D0*(Sr+(L2-L1)*SUMr/RCL)
         Si=0.5D0*(Si+(L2-L1)*SUMi/RCL)
      endif
      return
      end subroutine trapzd

!==========================================================
!================AUTOTROPHIC RESPIRATION===================
      real*8 function Rdark(Vcmax)
!@sum Rdark  Leaf dark respiration, Rd (umol m-2_leaf s-1)
!@+   From S. von Caemmerer (2000) Biochemical Models of Leaf Photosynthesis,
!@+   CSIRO book.
!@auth N.Y.Kiang

      real*8, intent(in) :: Vcmax

      Rdark = 0.015d0 * Vcmax !von Caemmerer book.

      end function Rdark

!---------------------------------------------------------------------!
      real*8 function Resp_cpool_maint(pft,C,CN,T_k,T_k_10d,facclim)
     &     Result(R_maint)
!@sum Maintenance respiration for a plant carbon pool of size C (umol/plant/s)
!@+   Based on biomass amount (total N in pool). From CLM3.0.
!@auth N.Y.Kiang
      !C3 vs. C4:  Byrd et al. (1992) showed no difference in maintenance
      ! respiration costs between C3 and C4 leaves in a lab growth study.
      ! Also, maintenance (dark) respiration showed no relation to
      ! leaf nitrogen content (assimilation and growth respiration did
      ! respond to leaf N content). In lab conditions, leaf dark respiration
      ! was about 1 umol-CO2 m-2 s-1 for an N range of ~70 to 155 mmol-N m-2.
! Re CLM parameters:
!       CLM calculates this per individual*population/area_fraction
!       to give flux per area of pft cover rather than per ground area.

!Other versions:
!        !*Original CASA *!
!     &       exp(308.56d0*(1/56.02d0 - (1/(T_k-227.13d0)))) *
!     &       ugBiomass_per_gC/ugBiomass_per_umolCO2
!        !*Acclimation vertical shift*! 56.02 = 10+273.15-227.13.  76.02 = 30+273.15-227.13.
!        !*Acclimation horizontal shift.
!     &       exp(308.56d0*
!     &       (1/56.02d0
!     &       - (1/(T_k-min(30.d0,max(10.d0,T_k_10d))+10.d0-227.13d0))))
!     &       * ugBiomass_per_gC/ugBiomass_per_umolCO2

      implicit none
      integer :: pft            !Plant functional type.
      real*8 :: C               !g-C/individual ! Weng, 01/11/2017: now, it's kgC.
                                !Can be leaf, stem, or root pools.
      real*8 :: CN              !C:N ratio of the respective pool
      real*8 :: T_k             !Temperature of canopy (Kelvin)
      real*8 :: T_k_10d         !Temperature of air 10-day average (Kelvin)
                                !  Should be canopy temp, but 10-day avg. okay.
      real*8 :: facclim         !frost-hardiness stress factor

      !---Local-------
      real*8,parameter :: k_CLM = 6.34d-07 !(s-1) rate from CLM.
      real*8,parameter :: ugBiomass_per_gC = 2.d6
      real*8,parameter :: ugBiomass_per_umolCO2 = 28.5 ! ?
      real*8,parameter :: umolCO2_per_kgC = 1000.0/12 * 1.d6 ! added by Weng, 01/10/2017

      if (T_k>228.15d0) then    ! set to cut-off at 45 deg C
        R_maint = facclim * BMEspdata(pft)%r * k_CLM * (C/CN) * ! unit of C is kgC
     &       exp(308.56d0*
     &       (1/min(max(56.02d0,T_k_10d-227.13d0),76.02d0)
     &       - (1/(T_k-227.13d0))))
     &       * umolCO2_per_kgC  ! Weng, 01/11/2017
      else
         R_maint = 0.d0
      endif
      end function Resp_cpool_maint
!---------------------------------------------------------------------!

      real*8 function Resp_sapwood_maint(pft,
     &  K_cambium,Acambium,T_k,T_k_10d,facclim)
     &     Result(R_maint)
!@sum Maintenance respiration for trees' sapwood (umol/plant/s)
!@+   Based on Weng et al. 2015, LM3-PPA
!@auth Ensheng Weng, 12-14-2017
      implicit none
      integer :: pft            ! Plant functional type.
      real*8 :: K_cambium, Acambium        ! m2/individual
      real*8 :: T_k             ! Temperature of canopy (Kelvin)
      real*8 :: T_k_10d         ! Temperature of air 10-day average (Kelvin)
                                ! Should be canopy temp, but 10-day avg. okay.
      real*8 :: facclim         !frost-hardiness stress factor

      !---Local-------
      real*8,parameter :: umolCO2_per_kgC = 1000.0/12 * 1.d6 ! added by Weng, 01/10/2017

      if (T_k>228.15d0) then    ! set to cut-off at 45 deg C
         R_maint  = K_cambium * Acambium * ! unit of C is kgC/yr/tree
     &       facclim * exp(308.56d0*
     &       (1/min(max(56.02d0,T_k_10d-227.13d0),76.02d0)
     &       - (1/(T_k-227.13d0))))
     &       * umolCO2_per_kgC ! umolC/tree/year
     &       / (365.0d0 *24.0d0*3600.0d0) ! year-1 --> s-1

      else
         R_maint = 0.d0
      endif
      end function Resp_sapwood_maint

!---------------------------------------------------------------------!
      real*8 function Resp_can_growth(pft,Acan,Rmaint,Rtgrowth)
     &     Result(R_growth)
!@sum Growth (light) respiration (units as input for Acan and Rmaint).
!@auth N.Y.Kiang
      !Based on photosynthetic activity. See Amthor (2000) review of
      ! Mcree - de Wit - Penning de Vries - Thornley respiration paradigms.
      ! See also Ruimy et al. (1996) analysis of growth_r.
      !Fixed to min 0.d0 like ED2. - NYK
      !Amthor (2000) range 0.39-0.77 for annual TOTAL respiration/GPP.
      integer :: pft
      real*8 :: Acan !Canopy photosynthesis rate (mass/m2/s)(or any units)
      real*8 :: Rmaint !Canopy maintenance respiration rate (mass/m2/s)
      real*8 :: Rtgrowth !Growth respiration from tissue growth (mass/m2/s)
      real*8 :: growth_r !pft-dependent. E.g.CLM3.0-0.25, ED2 conifer-0.53, ED2 hw-0.33

!      if (BMEspdata(pft)%leaftype.eq.NEEDLELEAF) then
      if (BMEspdata(pft)%woody) then
         growth_r = 0.30d0      !old 0.40 Total resp/GPP, low end, Amthor (2000)
                                !Guess 0.30 is portion that is growth_r.
      else
         growth_r = 0.28d0      !0.28 Value from Ruimy et al. (1996),model
      endif

      !If growth_r is actually growth_r component, then:
      !R_growth = max(0.d0, growth_r*(Acan))
      !If growth_r is actually total, so substract other resp components.
      R_growth = max(0.d0, growth_r*(Acan - Rmaint) - Rtgrowth) !Original
      !R_growth = max(0.d0, growth_r*(Acan - Rmaint - Rtgrowth)) !Rev test

      end function Resp_can_growth

!---------------------------------------------------------------------------
      function water_stress4(pft, nlayers, thetarel,
     &     fracroot, fice, betadl) Result(betad)
!@sum Rodriguez-Iturbe et al. (2001) water stress function.
!@+   Version if input is relative soil water content (saturated fraction).
!@+   Version if water access is equal throughout the root zone (not waited by root mass distr.).
!@+      Distribution betadl(:) scaled so that sum betad equals betak stress of least stress layer.
!@auth N.Y.Kiang
      implicit none
      integer,intent(in) :: pft  !Plant functional type number.
      integer,intent(in) :: nlayers !Number of soil layers
!      real*8,intent(in) ::  thetas(:) !Soil vol. water (vol.water/vol.soil)
      real*8,intent(in) ::  thetarel(:) !Relative soil vol. water (vol.water/vol. saturated)
      real*8,intent(in) :: fracroot(:) !Fraction of roots in layer
      real*8,intent(in) :: fice(:)  !Fraction of ice in layer
      real*8,intent(out) :: betadl(:) !Water stress in layers
      real*8 :: betad !Stress value, 0-1, 1=no stress, weighted by layers

      !Local vars
      integer :: k
      real*8 :: s  !Normalized soil moisture, s=thetas/thetasat
      real*8 :: betak  !Stress value for layer k
      real*8 :: rootzone(N_DEPTH)
      real*8 :: betadlsum

      !4. Rodriguez-Iturbe, Laio, & Porporato (2001 set) water stress fn
      !   with uniform depth access, and stress equal to the least stressed soil layer.
      betad = 0.d0
      rootzone(:) = 1.d0 ! added by Weng, 09/09/2021, for avoiding permanent zero fo betad
      do k=4,nlayers ! Weng changed 1 to 4, original is k=1
        if (fracroot(k).gt.0.d0) then
           rootzone(k) = 1.d0
        else
           rootzone(k) = 0.d0
        endif
      end do

      do k = 1,nlayers
        s = thetarel(k)
        ! Updated by Weng, 2021/11/09
        betak = Max(0.d0,Min(1.d0,(s-BMEspdata(pft)%swilt)/
     &         (BMEspdata(pft)%sstar-BMEspdata(pft)%swilt))) !Just linear

        !betadl(k) = (1.d0-fice(k))*rootzone(k)*betak
        betadl(k) = rootzone(k)*betak !Weng added, 2021/09/13
        betad = max( betad, betadl(k))  !Stress is of least stressed layer
      end do
      betadlsum = sum(betadl(:))
      if ( betadlsum < 1.d-12 ) then  ! to avoid 0/0 divisions
        betadl(:) = 0.d0
        betad = 0.d0
      else
        betadl(:) = betadl(:) * betad/betadlsum  !Scale so that betadl(:) sums to betad.
      endif
      if (betad < EPS2) betad=0.d0

      end function water_stress4
!---------------------------------------------------------------------------

      real*8 function frost_hardiness(Sacclim) Result(facclim)
!@sum frost_hardiness.  Calculate factor for adjusting photosynthetic capacity
!@+   due to frost hardiness phenology.
!@+   Based on Repo et al (1990), Hanninen & Kramer (2007),
!@+   and Makela et al (2006)
!@auth M.Puma
      real*8,intent(in) :: Sacclim
!      real*8 :: facclim ! acclimation/frost hardiness factor [0. to 1.]
      !----Local-----
      real*8,parameter :: Tacclim=-5.93d0 ! threshold temperature for photosynthesis [deg C]
      !real*8,parameter :: Tacclim=-3.d0 ! Best tune for Hyytiala
                        ! Site specific thres. temp.: state of photosyn.acclim
                        ! Hyytiala Scots Pine, -5.93 deg C Makela et al (2006)
      !real*8,parameter :: a_const=0.0595 ! factor to convert from Sacclim [degC] to facclim [-]
                        ! Site specific; conversion (1/Sacclim_max)=1/16.8115
                        ! estimated by using the max S from Hyytiala 1998
      real*8, parameter :: a_const = 0.1d0 !Closer tune for Hyytiala

      if (Sacclim > Tacclim) then ! photosynthesis occurs
         facclim = a_const * (Sacclim-Tacclim)
         if (facclim > 1.d0) facclim = 1.d0
!      elseif (Sacclim < -1E10)then !UNDEFINED
      elseif (Sacclim.eq.UNDEF)then !UNDEFINED
         facclim = 1.d0         ! no acclimation for this pft and/or simualtion
      else
         facclim = 0.01d0       ! arbitrary min value so that photosyn /= zero
      endif
      !facclim = 1.d0 ! Weng for test, 08/09/2020
      end function frost_hardiness

!-----------------------------------------------------------------------------
      real*8 function running_mean(dtsec,numd,var,var_mean)
!@sum Function for expiring running mean for phenological climate statistics.
!@+   Daily average, taking into account time step of function call.
      real*8, intent(in) :: dtsec
      real*8, intent(in) :: numd !number of days for running mean
      real*8, intent(in) :: var
      real*8, intent(in) :: var_mean
      real*8 :: zweight

      zweight=exp(-1.d0/(numd*86400.d0/dtsec))
      running_mean=zweight*var_mean+(1.d0-zweight)*var

      end function running_mean

!-----------------------------------------------------------------------------
      real*8 function snowmask_frac(ncover,snow_frac, snow_depth)
!@sum Fraction of patch (veg or bare) masked by snow.
!@+   Affects conductance like canopy wetness fraction
!@+   Follows VTFRAC calculation originally in ALBEDO.f.
!@+   Condition on snow_depth imitates ALBEDO.f.
      use ent_pfts
      implicit none
      integer,intent(in) :: ncover !Number of cover type, can be bare of veg
      real*8,intent(in) :: snow_frac(:), snow_depth(:) !1=bare, 2=veg
      !----Local-------
      integer :: pft

      if (snow_depth(1) + snow_depth(2) > 1.D-04 * 10.d0) then
      ! Condition imitates ALBEDO.f
         pft = ncover - COVEROFFSET
         if (pft.lt.COVEROFFSET.or.pft.gt.N_PFT) then !is bare soil
            snowmask_frac = snow_frac(1)*
     &           (1.d0 - EXP(-snow_depth(1)/VTMASKENT(ncover)) )
         else                   !is vegetated
            snowmask_frac = snow_frac(2)*
     &           (1.d0 - EXP(-snow_depth(2)/VTMASKENT(ncover)) )
         endif
      else
         snowmask_frac = 0.d0
      endif
      end function snowmask_frac

C---------------------------------------------------------------------------
      function calc_Ci_canopy(Ca,Gb, Gs,Anet,LAI,IPAR) Result(ci)
!@sum Foliage internal CO2 conc (mol mol-3) assuming diffusive flux of CO2
!@sum is at steady-state with biochemical uptake by photosynthesis and
!@sum that there is zero leaf boundary layer resistance (infinite gb),
!@sum and that there is no leaf cuticular conductance of CO2.
!@sum Full equation:  ci = ca - Anet*(1.37/gb + 1.65/gs)
!@sum 1.37 = ratio of diffusivities of CO2 and water vapor in laminar flow
!@sum       in the leaf boundary layer
!@sum 1.65 = ratio of diffusivities of CO2 and water vapor in still air at
!@sum       the leaf surface
!@sum (Monteith, 1995;  Kiang, 2003;  Collatz 1991)
!@sum ##### NOTE: Parameter Ball_b should actually be passed in with pspar.
!@sum #####       Have to set up for generic pspar for different photosynthesis
!@sum #####       routines.  (NK)
!@+   FOR DIAGNOSTIC/TESTING, NOT USED FOR REGULAR RUNS.

      implicit none

      real*8,intent(in) :: ca !CO2 mole fraction at surface reference height (umol mol-1)
      real*8,intent(in) :: Gb !Canopy boundary layer conductance of water vapor (mol m-2 s-1)
      real*8,intent(in) :: Gs !Canopy Stomatal conductance of water vapor(mol m-2 s-1)
      real*8,intent(in) :: Anet !Leaf net assimilation of CO2 (umol m-2 s-1)
      real*8,intent(in) :: LAI !Leaf area index
      real*8,intent(in) :: IPAR !Incident PAR (umol m-2 s-1)
      real*8 :: ci              !Leaf internal CO2 concentration (umol mol-1)
      !----Local------
      real*8,parameter :: MINPARMOL=50  !umol m-2 s-1
      real*8,parameter :: Ball_b = 0.01 !mol m-2 s-1

      if (IPAR.lt.MINPARMOL) then  !Stomates closed
        ci = ca - Anet*1.37/Ball_b
      else
        ci = ca - Anet*(1.37/Gb + 1.65/Gs) !LAI cancels in numerator and denominator.
      endif

      end function calc_Ci_canopy

!========================Photosynthesis =============================
      subroutine init_ci(ca, ci)
!@sum init_ci  Initialize leaf internal CO2 concentration.
      implicit none
      real*8,intent(in) :: ca   !Ambient air CO2 concentration (umol mol-1)
      real*8,intent(inout) :: ci !Leaf internal CO2 mole fraction  (umol mol-1)

      !ci should be initialized. For rule of thumb, initialize to typically
      !observed ratio:
      ci = 0.7d0*ca

      end subroutine init_ci

!-----------------------------------------------------------------------------
      subroutine biophysdrv_setup(ca,ci,Tc,Pa,rh,psdrvpar)
!@sum Set up met drivers for photosynthesis at every physical time step.
      implicit none
      real*8,intent(in) :: ca, ci, Tc, Pa, rh
      type(psdrvtype),intent(out) :: psdrvpar

      psdrvpar%ca = ca
      psdrvpar%ci = ci
      psdrvpar%Tc = Tc
      psdrvpar%Pa = Pa
      psdrvpar%rh = rh

      end subroutine biophysdrv_setup

!-----------------------------------------------------------------------------
      subroutine pscondleaf(pft,IPAR,psd,Gb,gsout,Aout,Rdout
     &     ,sunlitshaded,ISPout)
!@sum pscondleaf  Main routine to obtain leaf-level photosynthesis,
!@+   stomatal conductance, and dark respiration.
!@+   See Photosynth_analyticsoln for units.
      implicit none
      integer,intent(in) :: pft
      real*8,intent(in) :: IPAR !umol m-2 s-1. Absorbed PAR. Should APAR.
      type(psdrvtype) :: psd
      real*8,intent(in) :: Gb !mol m-2 s-1
      real*8,intent(out) :: gsout, Aout, Rdout !ci recorded in psd
      real*8,intent(out) :: ISPout
      integer,intent(in) :: sunlitshaded
      !---Local---
      real*8 :: ci, cs
      !real*8,parameter :: LOW_LIGHT_LIMIT = 2.5d0 !umol m-2 s-1.  Nobel 1999, lower light limit for green plants is 0.7 W m-2 ~ 3 umol m-2 s-1.
      real*8,parameter :: LOW_LIGHT_LIMIT = 3.0d0 !umol m-2 s-1.  Nobel 1999, lower light limit for green plants is 0.7 W m-2 ~ 3 umol m-2 s-1.

        call Photosynth_analyticsoln(pft,IPAR,psd%ca,ci,
     &     psd%Tc,psd%Pa,psd%rh,Gb,gsout,Aout,Rdout,sunlitshaded,
     &  ISPout)

        psd%ci = ci     !Ball-Berry:  ci is analytically solved.  F-K: ci saved between time steps.

      end subroutine pscondleaf

!-----------------------------------------------------------------------------
      subroutine Photosynth_analyticsoln(pft,IPAR,ca,ci,Tl,Pa,rh,gb,
     o     gs,Atot,Rd,sunlitshaded,isp)
!@sum Photosynth_analyticsoln  Selects the correct root for the
!@+   solution to the cubic coupled equations of  Farquhar photosynthesis and
!@+   Ball-Berry conductance, including dark respiration.
!@+   ci is solved for analytically for each of the limiting cases.
!@+   Outputs gs, Atot, Rd. May output also other VOC fluxes.
!@auth  N.Y.Kiang, I.Aleinov
      use ent_pfts, only : BMEspdata
      implicit none
      integer,intent(in) :: pft !Plant functional type, 1-C3 grassland
      real*8,intent(in) :: IPAR !Absorbed PAR.  WRONG OLD COMMENT:Incident PAR (umol m-2 s-1)
      real*8,intent(in) :: ca   !Ambient air CO2 mole fraction (umol mol-1)
      real*8,intent(in) :: Tl   !Leaf temperature (Celsius)
      real*8,intent(in) :: rh   !Relative humidity
      real*8,intent(in) :: Pa   !Pressure (Pa)
      real*8,intent(in) :: gb   !Leaf boundary layer conductance of water vapor (mol m-2 s-1)
      integer,intent(in) :: sunlitshaded !For diagnostic outputs only.
      real*8,intent(out) :: ci   !Leaf internal CO2 concentration (umol mol-1)
      real*8,intent(out) :: gs  !Leaf stomatal conductance (mol-H2O m-2 s-1)
      real*8,intent(out) :: Atot !Leaf gross photosynthesis (CO2 uptake, micromol m-2 s-1)
      real*8,intent(out) :: Rd  !Dark = above-ground growth + maintenance respiration (umol m-2 s-1)
      real*8,intent(out) :: isp ! Isoprene emission (umol C m-2 s-1)
        !---Local----
      real*8,parameter :: O2pres=20900.d0 !O2 partial pressure in leaf (Pa) Not exactly .209*101325.
      real*8 :: cie, cic, cis   !Leaf internal CO2 (umol mol-1)
!      real*8 :: Je1, Jc1, Js1   !Assimilation of CO2, 3 limiting cases
      real*8 :: Anet            !Net assimilation of CO2 = Atot - aboveground respir (umol m-2 s-1)
      real*8 :: Aiso            ! Rate of photosynthesis for isoprene emissions (umol m-2 s-1)
      real*8 :: cs   !CO2 mole fraction at the leaf surface (umol mol-1)
      real*8 :: Ae, Ac, As      !* These are Anet!
      real*8, save :: a1c=1.d30, f1c=-1.d30
      real*8 :: a1e, f1e
      real*8, parameter :: alpha=.08d0 !Intrinsic quantum efficiency for CO2 uptake
#ifdef PS_BVOC
      logical, parameter :: need_isoprene = .true.
#else
      logical, parameter :: need_isoprene = .false.
#endif
      integer, save :: counter = 0
      counter = counter + 1

      !write(888,*) "counter=", counter
      Rd = gamma_r * pspar%Vcmax

      if ( IPAR < .000001d0 ) then
        !print *, 'IPAR<.000001d0',Rd,ca,gb  !NK DEBUG
        Atot = 0.d0
        Anet = - Rd
        cs = ca - Anet*1.37d0/gb
        gs = BallBerry(Anet, rh, cs, pspar)
        ci = cs - Anet/(gs/1.65d0)
        isp = 0.d0
        return
      endif


      !* Photosynthetic rate limited by RuBP saturation
      !* Assimilation is of the form a1*(Ci - Gammastar)/(e1*Ci + f)
!      call Ci_Jc(ca,gb,rh,IPAR,Pa,pspar, Rd,O2pres, cic, Jc1)
      ! Jc_RuBP = pspar%Vcmax*(Cip - pspar%Gammastar)/
      !           (Cip + pspar%Kc*(1 + O2/pspar%Ko))

      !Assimilation is of the form a1*(Ci - Gammastar)/(e1*Ci + f)
      if ( pspar%first_call ) then
           !NK DEBUG
         if (BMEspdata(pspar%pft)%pst.eq.C3) then
            a1c = pspar%Vcmax
            f1c = pspar%Kc*(1.d0 + O2pres/pspar%Ko) * 1.d06/Pa !umol/mol
           !NK DEBUG
            !call ci_cubic (ca,rh,gb,Pa,Rd,a1c,f1c,pspar,Axxx)
            call ci_cubic(ca,rh,gb,Pa,Rd,a1c,f1c,pspar,Ac)
            !if ( Ac >= -Rd ) write(578,*) Axxx, Ac, Ac - Axxx
            !write(888,*) "Ac", ca,rh,gb,Pa,Rd,a1,f1,pspar,Ac
         else !C4 photosynthesis
           !NK DEBUG
            Ac = pspar%Vcmax - Rd
         endif
         pspar%Ac = Ac
         !pspar%first_call = .false. !Probably bug-prone, but reset after As
      else
         Ac = pspar%Ac
      endif

!      call Ci_Je(ca,gb,rh,IPAR,Pa, pspar, Rd, cie, Je1)
      ! Photosynthetic rate limited by light electron transport (umol m-2 s-1)
      ! Je_light = (pspar%PARabsorb*IPAR)*alpha*(Cip-pspar%Gammastar)/
      !            (Cip+2*pspar%Gammastar)

      !Assimilation is of the form a1*(ci - Gammastar.umol)/(e1*ci + f1)

      if (BMEspdata(pspar%pft)%pst.eq.C3) then
         !a1 = pspar%PARabsorb*IPAR*alpha
        a1e = IPAR*alpha        !### HACK:  IPAR from canopyspitters.f is APAR.  When we switch to Wenze's canopyrad, then leaf PARabsorb will be used -NK ###
        f1e = 2*pspar%Gammastar * 1.d06/Pa !Convert from Pa to umol/mol

        if ( a1e < a1c .or.
     &       f1e > f1c .or.
     &       need_isoprene ) then
            !call ci_cubic (ca,rh,gb,Pa,Rd,a1e,f1e,pspar,Axxx)
          call ci_cubic(ca,rh,gb,Pa,Rd,a1e,f1e,pspar,Ae)
            !write(888,*) "Ae", ca,rh,gb,Pa,Rd,a1,f1,pspar,Ae
cddd        call ci_cubic1(ca,rh,gb,Pa,Rd,a1,f1,pspar,Axxx)
cddd        write(579,*) Ae, Axxx
cddd        if ( Ae > 0.d0 ) write(578,*) Axxx - Ae
        !if ( Ae >= -Rd ) write(578,*) Axxx, Ae, Ae - Axxx
        else
          Ae = 1.d30
        endif
      else                      !C4 photosynthesis
        Ae = IPAR*alpha - Rd
      endif


      !* Photosynthetic rate limited by utilization of photosynthetic products:
      !* (umol m-2 s-1)  Triosphosphate (TPU limitation for C3,
      !*                 PEP carboxylase limitation for C4.
      if (pspar%first_call) then
        if (BMEspdata(pspar%pft)%pst.eq.C3) then
           !call Ci_Js(ca,gb,rh,IPAR,Pa,pspar,Rd, cis, Js1)
           !Js_sucrose = pspar%Vcmax/2.d0
          As = pspar%Vcmax/2.d0 - Rd !Anet
           !write(888,*) "As", As
        else                    !C4 photosynthesis
            !As = 4000.d0*pspar%Vcmax*ci - Rd
          call Asnet_C4(ca,rh,gb,Rd,pspar,As) !This is Anet
        endif
        pspar%As = As
        pspar%first_call = .false.
      else
        As = pspar%As
      endif

      !Anet = min(Ae, Ac, As)
      !Limit flux for numerical stability, keep cs>0.2*ca
      Anet = min(Ae, Ac, As, 0.8d0*ca*gb/1.37d0)
      Atot = Anet + Rd
      Aiso = Ae + Rd

      if (Atot.lt.0.d0) then
      ! can only happen if ca < Gammastar . Does it make sense? -Yes-NK
#ifdef OFFLINE
         write(997,*) "Error, Atot<0.0:",Atot,Ae,Ac,As,ca,gb,rh,IPAR
     &        ,Pa,pspar,sunlitshaded, pspar%Gammastar * 1.d06/Pa
#endif
         Atot = 0.d0
         Anet = - Rd
!!       ci = pspar%Gammastar * 1.d06/Pa
!!       gs = 0. ! MK: setting to 0 to avoid erratic results

      endif

      cs = ca - Anet*1.37d0/gb
      gs = BallBerry(Anet, rh, cs, pspar)
      ci = cs - Anet/(gs/1.65d0)


#ifdef PS_BVOC
         call Voccalc(pft,pa,ca,ci,Tl,pspar%Gammastar,
     & isp,Aiso)
#else
       isp=0.0d0
#endif

      !write(888,*) "gs,ci,cs", gs,ci,cs

      end subroutine Photosynth_analyticsoln

!-----------------------------------------------------------------------------
      subroutine Voccalc(pft,pa,ca,ci,Tl,Gammastar,isp,Aiso)
!@sum Voccalc  Isoprene emissions coupled to photosynthesis
!@auth Nadine Unger
      use ent_pfts

      implicit none
      integer,intent(in) :: pft !Plant functional type, 1-C3 grassland
      real*8,intent(in) :: ca   !Ambient air CO2 mole fraction (umol mol-1)
      real*8,intent(in) :: Pa   !Pressure (Pa)
      real*8,intent(in) :: Tl   !Leaf temperature (Celsius)
      real*8,intent(in) :: Gammastar
      real*8,intent(in) :: Aiso   !(umol m-2 s-1)
      real*8,intent(in) :: ci
      real*8,intent(out) :: isp ! isoprene emission (umol C m-2 s-1)
      type(photosynthpar) :: pspar
        !---Local----
      integer, parameter :: numpft = 8
      real*8 :: gammamol,fact
      real*8 :: IBASER, Y_alpha, kapco2
      real*8 :: tauiso
C April 2009 values
c      real*8, parameter, dimension(numpft) :: Y_eps = !
c     & (/0.0d0,2.10d-02,8.24d-02,6.48d-02,1.08d-01,
c     & 4.44d-02,1.38d-01,0.0d0/)
C New values June 2009
C      real*8, parameter, dimension(numpft) :: Y_eps = !
C     & (/0.0d0,1.91d-02,7.19d-02,5.13d-02,8.79d-02,
C     & 2.18d-02,8.35d-02,0.0d0/)

      gammamol =  Gammastar * 1.d06/Pa !Convert from Pa to umol/mol

C Y_alpha, Y_eps unitless

      Y_alpha=(ci-gammamol)/(6.0*(4.67*ci+9.33*gammamol))

c      isp = Y_eps(pft)*Aiso*Y_alpha

!!!      isp = BMEspdata(pft)%Y_eps*Aiso*Y_alpha
!!hack
      isp = 0

C Include CO2 effects

       kapco2 = (0.7*370.0)/(0.7*ca)

C Include temperature effects

       tauiso = exp(0.1*(Tl-30.0))


C Include seasonal effects? Add later.
C Note can switch on and off kapco2

       isp = isp*kapco2*tauiso

      end subroutine Voccalc

!-----------------------------------------------------------------------------
      function calc_CO2compp(O2,Kc,Ko,Tl) Result(Gammastar)
!@sum calc_CO2compp CO2 compensation point in absence of dark respiration (Pa)

      implicit none
      real*8,intent(in) :: O2 !O2 partial pressure in leaf (Pa)
      real*8,intent(in) :: Kc   !Michaelis-Menten parameter for CO2 (Pa)
      real*8,intent(in) :: Ko   !Michaelis-Menten parameter for O2 (Pa)
      real*8,intent(in) :: Tl !Leaf temperature (Celsius)
      real*8 :: Gammastar  !CO2 compensation point (Pa)
      !----Local-----
      real*8,parameter :: tau=2600.d0  !CO2/O2-specificity ratio

!      Gammastar = O2/(2.d0*tau*Q10fn(0.57d0,Tl)) !Collatz (A3)
!      Gammastar = O2*Q10fn(1.75,Tl)/(2.d0*tau) !Collatz (A3) Same as above, !KcQ10/KoQ10 = 2.1/1.2 = 1.75 = 1/.57
!      Gammastar = 0.5d0*(Kc/Ko)*0.21*O2 !CLM Tech Note. Gives smaller Gammastar than Collatz.

C Nadine - use previous T-dep version
#ifdef PS_BVOC
        Gammastar = O2/(2.d0*tau*Q10fn(0.57d0,Tl)) !Collatz (A3)
#else
       Gammastar = 0.5d0*(Kc/Ko)*0.21*O2 !CLM Tech Note. Gives smaller Gammastar than Collatz.
#endif

      end function calc_CO2compp
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

      function arrhenius(Tcelsius,c1,c2) Result(arrh)
!@sum arrhenius Arrhenius response to temperature for biological kinetics.
!@+   From David Medvigy's lphys.f90
      implicit none
      real*8 :: Tcelsius
      real*8 :: c1,c2 !Process-specific temperature response parameters.
      real*8 :: arrh

      arrh = c1*exp(c2*(1.d0/288.15d0-1.d0/(Tcelsius+Kelvin)))
      return
      end function arrhenius
!=================================================
      function Q10fn(Q10par,Tcelsius) Result(Q10factor)
!@sum Q10fn  Q10 function, biological response to temperature.
!@+   From Collatz, et al. (1991)
      implicit none
      real*8 :: Q10par, Tcelsius !parameter, temperature
      real*8 :: Q10factor

      Q10factor = Q10par**((Tcelsius-25.d0)/10.d0)

      end function Q10fn
!=================================================

      function Tresponse(c,deltaH,Tcelsius) Result(Tfactor)
!@sum Tresponse  Arrhenius temperature response function that accounts for
!@+   activation energy.
!@+   From Bernacchi, et al. (2001).
      implicit none
      real*8,intent(in) :: c !Scaling factor
      real*8,intent(in) :: deltaH !Activation energy
      real*8,intent(in) :: Tcelsius !Temperature (Celsius)
      real*8 :: Tfactor

      Tfactor = exp(c - deltaH/(Rgas * (Tcelsius + Kelvin)))

      end function Tresponse
!=================================================

!#define ENT_CHECK_C4_SOLUTION
      subroutine  Asnet_C4(ca,rh,gb,Rd,pspar,Asnet)
!@sum Asnet_C4 PEP carboxlase-limited carbon assimilation for C4 photosynthesis
!@+   After Collatz, and CLM's correction of the coefficient
!@+   Returns Asnet = Astot - Rd
!@+   Solving for Asnet via the equations:
!@+   1) Asnet = Astot - Rd
!@+           = (ca - ci)/[(1.37*rb + 1.65*rs)]
!@+           = (ca - cs)/(1.37*rb)
!@+           = (cs - ci)/(1.65*rs)
!@+   2) 1/rs = gs = m*A*rh/cs + b
!@+   3) Astot = 4000.d0*pspar%Vcmax*(ci*1e-06)
!@+        !4000 is CLM, Collatz had 1800. Convert ci from umol/mol to mol/mol
!@+   Do subsitutions to eliminate cs and ci and solve for Asnet.

      implicit none
      real*8,intent(in) :: ca  !Surface air CO2 concentration (umol/mol)
      real*8,intent(in) :: rh   !Relative humidity
      real*8,intent(in) :: gb   !Leaf boundary layer conductance of water vapor (mol m-2 s-1)
      real*8,intent(in) :: Rd   !Leaf mitochondrial respiration (umol m-2 s-1)
      type(photosynthpar) :: pspar
      real*8,intent(out) :: Asnet  !Net assimilation of carbon (umol m-2 s-1)
      !---Local----
      real*8 :: ci !Leaf internal CO2 concentration (umol/mol)
      real*8 :: K1, K2, K3, K4
      real*8 :: X, Y, Z
      real*8 :: b0, a0, c0, sqrtop
      real*8 :: Aspos, Asneg
#ifdef ENT_CHECK_C4_SOLUTION
      real*8 ::  cs, gs, Acheck, Asnet1, Asnet2
#endif

      K1 = pspar%m * rh
      !print *,'K1',K1
      K2 = 4000.d0*pspar%Vcmax * 1.d-06
      !print *, 'K2',K2
      K3 = 1.37d0/gb
      !print *, 'K3',K3
      K4 = 1.65d0*ca

      X = (K1 - pspar%b*K3)*(1/K2 + K3) - 1.65d0*K3

      Y = ca*(pspar%b/K2 - K1 + 1.65d0 + 2.d0*K3*pspar%b)
     &     + Rd/K2*(K1 - pspar%b*K3)

      Z = pspar%b * ca * ( Rd/K2 - ca )
      sqrtop = Y**2.d0 - 4.d0*X*Z
      !ERROR CHECK
      if ((sqrtop.lt.0.d0).or.(X.eq.0.d0)) then
         !Asnet = 0.8d0*ca*gb/1.37d0  !Set to upper limit
         Asnet = -Rd  !Stomatal shutdown caused by m=0, all ice.
         return
      endif
      Asnet = max (-Rd
     &     ,(-Y + sqrt(sqrtop))/(2.d0*X)) !Positive root is max.
      !Asnet = min(Asnet,  0.8d0*ca*gb/1.37d0)

#ifdef ENT_CHECK_C4_SOLUTION
      if ( Asnet < 0.d0 ) return ! willnot check

      !!! checking the solution
      ci = (Asnet + Rd) / (4000.d0*pspar%Vcmax*1.d-06)
      cs = ca - Asnet * 1.37d0/gb
      gs = pspar%m*Asnet*rh/cs + pspar%b
      Acheck = (cs - ci)*gs/1.65d0

      if ( abs(Acheck-Asnet) > 1.d-10 .or. cs < 1.d-5 ) then
      ! &     .or. Asnet >  0.8d0*ca*gb/1.37d0 ) then
        write(889,'(E15.5, 10f15.3)') Acheck-Asnet,
     &       Asnet, 0.8d0*ca*gb/1.37d0,
     &       cs, ci
      endif

#endif

      end subroutine Asnet_C4

!=================================================

#ifndef USE_NR_SOLVER_FOR_FBB
      subroutine ci_cubic(ca,rh,gb,Pa,Rd,a1,f1,pspar,A)
!@sum ci_cubic Analytical solution for cubic equation of coupled
!@+   Ball-Berry/Farquhar stomatal conductance/photosynthesis.
!@+   Version that uses analytical equation solution.
!@+   Solves for Anet.
!@auth I.Aleinov
      !@sum ci (umol/mol)
      !@sum For the case of assimilation being of the form:
      !@sum         A = Anet = a*(Cip - Gammastar)/(e*Cip + f) - Rd
      !@sum Numerator and denominator are converted from (Pa/Pa) to (umol mol-1)/(umol mol-1)
      !@sum         A = Anet = a1*(ci - gammamol) /(e1*ci + fmol) - Rd
      !@sum where gammamol = Gammastar*1d06/Pa, fmol = f1 = f*1d06/Pa

      implicit none
      real*8 :: ca              !Ambient air CO2 concentration (umol mol-1)
      real*8 :: rh              !Relative humidity
      real*8 :: gb              !Leaf boundary layer conductance of water vapor (mol m-2 s-1)
      real*8 :: Pa              !Pressure (Pa)
      real*8 :: Rd              !Leaf mitochondrial respiration (umol m-2 s-1)
      real*8 :: a1              !Coefficient in linear Farquhar equ.
      real*8 :: f1              !Coefficient in linear Farquhar equ.
      real*8, intent(out) :: A
      type(photosynthpar) :: pspar
      !----Local----
      real*8, parameter :: S_ATM=1.37d0  ! diffusivity ratio H2O/CO2 (atmosph.)
      real*8, parameter :: S_STOM=1.65d0 ! diffusivity ratio H2O/CO2 (stomatal)
      real*8 :: Ra, b, K, gamol, A_d_asymp
      real*8 :: X, Y, Z, Y1  ! tmp vars
      real*8 :: c3, c2, c1, c   !Coefficients of the cubic of ci (c3*ci^3 + c2*ci^2 + c1*ci + c)
      real*8 :: cixx(3) ! solutions of cubic
!      real*8 :: cs, Rs ! needed to compute ci
      integer :: nroots, i
!      real*8 ci

      Ra = 1/gb * S_ATM
      b = pspar%b / S_STOM
      K = pspar%m * rh / S_STOM
      gamol = pspar%Gammastar * 1.d06/Pa !Convert Pa to umol/mol
      A_d_asymp = - b*Ca / (K - b*Ra) ! asymptotic val of A from diffusion eq.

      ! first check some special cases
      if ( A_d_asymp >= 0.d0 ) then
        ! this can happen only for very low humidity
        ! probably should never happen in the real world, but if it does,
        ! this case should be considered separately
        !!print *,"!!! A_d_asymp >= 0.d0 !!!", A_d_asymp
        A_d_asymp = -1.d30 !!! hack
        !!print *,"K<b*Ra: m,rh,b,Ra:",pspar%m,rh,b,Ra
        !!call stop_model("ci_cubic: rh too small ?",255)
      endif

      ! dependence on e1 if needed
cddd      Y= f1/e1
cddd      X= -a1/e1 * (gamol+f1/e1)
cddd      Z= a1/e1 -Rd
      Y= f1
      X= -a1 * (gamol+f1)
      Z= a1 -Rd

      if ( Z + X/(Ca+Y) > 0.d0 ) then
        ! Farquhar curve is above zero. May have solution A > 0
        c = -(b*Ca*(X + (Ca + Y)*Z))
        c1 = Ca*Z - K*(X + Ca*Z + Y*Z) +
     &       b*(Ca**2 + Ca*(Y + 2*Ra*Z) + Ra*(X + Y*Z))
        c2 = Ca*(-1 + K - 2*b*Ra) + K*(Y + Ra*Z) - Ra*(b*Y + Z + b*Ra*Z)
        c3 = Ra*(1 - K + b*Ra)

        call cubicroot(c3, c2, c1, c, cixx, nroots)

        !!print *,"roots= ", cixx(1:nroots)

        ! find minimal root above the asymptotic value
        A = 1.d30
        do i=1,nroots
          if ( cixx(i) < A .and. cixx(i) > A_d_asymp ) A = cixx(i)
        enddo
        if ( A == 1.d30 )  then
          print *," m,rh,b,Ra:",pspar%m,rh,b,Ra
          print *,"ca,gb,Pa:",ca,gb,Pa
          print *,"pspar:",pspar
          print *," A_d_asymp,K,gamol,f1,a1,Rd",
     &         A_d_asymp,K,gamol,f1,a1,Rd
          print *,"c3,c2,c1,c", c3,c2,c1,c
          print *,"nroots,cixx",nroots,cixx(1:nroots)
          call stop_model("ci_cubic: no solution",255)
        endif

        if ( A >= 0 ) then
cddd          cs = ca - A*Ra
cddd          Rs = 1.d0 / (K*A/cs + b)
cddd          ci = cs - A*Rs
cddd          ! just in case, check consistency
cddd          if ( ci < 0.d0 ) call stop_model("ci_cubic: ci<0",255)
cddd          if ( cs < 0.d0 ) call stop_model("ci_cubic: cs<0",255)
cddd          !!print *,'QQQQ ',A,ci
          return
        endif

      endif

      ! if we got here then A<0 : have to solve quaratic equation

      Y1 = Y + ca
      c2 = Ra + 1.d0/b
      c1 = - (Y1 + c2*Z)
      c  = X + Y1*Z

      ! just in case,
      if (  c1*c1 - 4.d0*c2*c < 0.d0 )
     &     call stop_model("ci_cubic: no solution to quadratic",255)
      A = ( - c1 - sqrt( c1*c1 - 4.d0*c2*c ) ) / ( 2.d0 * c2 )
cddd      cs = ca - A*Ra
cddd      Rs = 1.d0 / ( b)
cddd      ci = cs - A*Rs
cddd      !!print *,"q ", ci, cs, A, Rs, Ra
cddd      ! just in case, check consistency
cddd      if ( ci < 0.d0 ) call stop_model("ci_cubic: q: ci<0",255)
cddd      if ( cs < 0.d0 ) call stop_model("ci_cubic: q: cs<0",255)
cddd      !!print *,'QQQQ ',A,ci

      end subroutine ci_cubic


!=================================================
      subroutine cubicroot(a,b,c,d,x,n)
!@sum cubicroot  Solve cubic equation: a x^3 + b x^2 + c x + d = 0 *!
!@auth I.Aleinov
      !* Written by Igor Aleinov from solution by Cardano in
      !* Korn, Korn, Mathematical Handbook.
      implicit none
      real*8,intent(in) :: a,b,c,d  ! coefficients of cubic
      real*8, intent(out) :: x(:)   ! results ( 0-3 roots )
      integer, intent(out) :: n     ! number of roots
      real*8 :: x0,x1,x2
      real*8 :: a0,a1,a2,Q1,R1,D1
      real*8, parameter :: EPS0 = 1.d-8 ! 1.d-15
      real*8, parameter :: one3rd = 1.d0/3.d0
      real*8 :: arg, S, T
      complex*16 :: ST

      !print *,"cubicroot:",a,b,c,d

      if (abs(a) < (abs(b)+abs(c)+abs(d))*EPS0 ) then
        if (abs(b) < (abs(c)+abs(d))*EPS0) then
          if (abs(c) < abs(d)*EPS0) then
            write(*,*) "Internal Error in Cardano: no solution."
            stop
          endif
          x0 = -d/c
          x(1) = x0
          !write(*,*) "Cardano: returning",x0
          n = 1
        else
          !write(*,*) "What's this?"
          D1 = c*c - 4.d0*b*d

          if (D1 > 0.d0) then
            Q1 = sqrt(D1)
            x0 = (-c + Q1) / (2.d0 * b)
            x1 = (-c - Q1) / (2.d0 * b)
            !return
            n = 2
          else if (D1.eq.0.) then
            x0 = -c / (2.d0 * b)
            x1 = x0
            n = 1
          else
            x0 = -c /(2.d0 *b)
            x1 = sqrt(-D1) / (2.d0* b)
            n = 0
          end if
        end if
        !print *,"CX1",x0,x1
        !x = max(x0,x1)
        x(1) = x0
        x(2) = x1
      else
        a2 = b/a
        a1 = c/a
        a0 = d/a
        Q1 = (3.d0 * a1 - a2*a2 ) / 9.d0
        R1 = (9.d0 * a2 * a1 - 27.d0 * a0 - 2.d0 * a2*a2*a2) /54.d0
        D1 = Q1*Q1*Q1 + R1*R1
        !write(*,*) "abcda2a1a0Q1R1D1",a,b,c,d,a2,a1,a0,Q1,R1,D1
        if (D1 > 0.d0) then       !* only one real root *!
          !write(*,*) "One real root."
          arg = R1 + sqrt(D1)
          S = sign(1.d0, arg) * (abs(arg)**one3rd)
          arg = R1 - sqrt(D1)
          T = sign(1.d0, arg) * (abs(arg)**one3rd)
          x0 = -a2/3.d0 + S + T
          x1 = -a2/3.d0 - (S+T)*0.5d0
          x2 = sqrt(3.d0) * (S-T)*0.5d0
          !print *,"CX2",x0,x1,x2
          n = 1
        else if (D1.eq.0.) then !* two roots coincide * *!
          !write(*,*) "Two roots coincide."
          S = sign(1.d0, R1) * (abs(R1)**one3rd)
          x0 = -a2/3.d0 + 2.d0*S
          x1 = -a2/3.d0 - S
          x2 = x1
          !print *,"CX3",x0,x1,x2
          n =2
        else                    !* three different real roots *!
          !call CRtCube( R1, sqrt(-D1), S, T)
          !write(*,*) "Three different real roots. a2R1D1ST",a2,R1,D1,S,T
          ST = ( cmplx(R1, sqrt(-D1),kind(1.d0)) )**one3rd
          S = real (ST)
          T = aimag(ST)
          x0 = -a2/3.d0 + 2.d0*S
          x1 = -a2/3.d0 - S + sqrt(3.d0)*T
          x2 = -a2/3.d0 - S - sqrt(3.d0)*T
          !print *,"CX4",x0,x1,x2
          n = 3
        end if
        !x = max(x0,x1,x2)
        !x = x2
        x(1) = x0
        x(2) = x1
        x(3) = x2
      end if
      end subroutine cubicroot
#endif

!=================================================

      subroutine calc_Pspar(dtsec,pft,Pa,Tl,O2pres,stressH2O,
     &                      Sacclim,llspan)
!@sum calc_Pspar Collatz photosynthesis parameters in data structure pspar.
!@    pspar is GLOBAL TO MODULE.
      !Later need to replace these with von Caemmerer book Arrhenius
      !function sensitivities (her Table 2.3)
      implicit none
      integer,intent(in) :: pft   !Plant functional type, 1=C3 grassland
      real*8,intent(in) :: dtsec
      real*8,intent(in) :: Pa     !Atmospheric pressure (Pa)
      real*8,intent(in) :: Tl     !Leaf temperature (Celsius)
      real*8,intent(in) :: O2pres !O2 partial pressure in leaf (Pa)
      real*8,intent(in) :: stressH2O
      real*8,intent(in) :: Sacclim !state of acclimation/frost hardiness
      real*8,intent(in) :: llspan !mean leaf life span
!      type(photosynthpar),intent(inout) :: pspar !Moved to global to module.
      integer :: p

      !----Local-----
      real*8 :: facclim ! acclimation/forst hardiness factor [-]
      real*8 :: tk, fparlimit !light(i.e.,PAR) control
      integer, save :: counter = 0
      counter = counter + 1
      tk = tl + 273.2d0
      facclim = frost_hardiness(Sacclim)

!      fparlimit = par_phenology(pft,llspan)
      fparlimit = 1.d0

!!! this var is not reproducible on restart, please figure out why
!      fparlimit = 1.d0 ! seems to be ok now
      p = pft
      pspar%pft = pft
      pspar%PARabsorb = BMEspdata(p)%PARabsorb !Collatz et al. (1991)
!      pspar%Vcmax = BMEspdata(p)%Vcmax/(1 + exp((-220.e03+703.*(Tl+Kelvin))
!     &     /(Rgas*(Tl+Kelvin))))
! Commented by Weng, 08/17/2020
!      pspar%Vcmax = BMEspdata(p)%Vcmax * Q10fn(2.21d0, Tl)
!     &            * facclim * fparlimit
!!hack
!      pspar%Kc = Kc*Q10fn(KcQ10,Tl) !(Collatz, eq. A12)
!      pspar%Ko = Ko*Q10fn(KoQ10,Tl) !(Collatz, eq. A12)
      !  in BiomeE (LM3-PPA)
      ! ko=0.248    * exp(35948/Rgas*(1.0/298.2-1.0/tl))*p_sea/p_surf ! Weng, 2013-01-10
      ! kc=0.000404 * exp(59356/Rgas*(1.0/298.2-1.0/tl))*p_sea/p_surf ! Weng, 2013-01-10
      ! vm=spdata(pft)%Vmax*exp(24920/Rgas*(1.0/298.2-1.0/tl))
      pspar%Vcmax = BMEspdata(p)%Vcmax
     &      * exp(24920.d0/Rgas*(1.d0/298.2d0-1.d0/tk))
     &      * facclim * fparlimit

      pspar%Kc = 0.0004d0 *exp(59356.d0/Rgas*(1.d0/298.2d0-1.d0/tk))
      pspar%Ko = 0.248d0  *exp(35948.d0/Rgas*(1.d0/298.2d0-1.d0/tk))

      pspar%Gammastar = calc_CO2compp(O2pres,pspar%Kc,pspar%Ko,Tl) !(Pa) (Collatz)
      !Slope of Ball-Berry equation (Collatz)
      pspar%m = BMEspdata(p)%m * stressH2O ! Weng, remove water stresses on m
      !pspar%m = BMEspdata(p)%m     !Slope of Ball-Berry equation (Collatz)
      pspar%b = BMEspdata(p)%b     !Intercept of Ball-Berry equation (mol m-2 s-1) (Collatz)
      pspar%Nleaf = BMEspdata(p)%Nleaf !g-N/m^2[leaf] Needed for foliar respiration.
      !pspar%Nleaf = BMEspdata(p)%Nleaf * phenology factor !Here can adjust Nleaf according
                                !to foliage N pools or phenology factor.
      pspar%stressH2O = stressH2O

      pspar%first_call = .true.
!      pspar%As = 0.d0 !Unnecessary but zero anyway
!      pspar%Ac = 0.d0 !Unnecessary but zero anyway
      pspar%reset_ci_cubic1 = .true.

      end subroutine calc_Pspar

!-----------------------------------------------------------------------------

      function BallBerry(Anet, rh, cs, pspar) Result (gsw)
!@sum Ball-Berry (1987) model of leaf stomatal conductance of
!@    water vapor, gsw (mol m-2 s-1)
!@auth N.Y.Kiang
      implicit none
      real*8,intent(in) :: Anet !Net assimilation of CO2 (umol m-2 s-1)
      real*8,intent(in) :: rh   !Relative humidity (fractional ratio)
      real*8,intent(in) :: cs   !Leaf surface CO2 mole fraction (umol mol-1)
      type(photosynthpar) :: pspar
      real*8 :: gsw !Leaf conductance of water vapor (mol m-2 s-1)
      !----Local-----

      ! just in case check cs (remove after debugging ?)
      if ( cs <= 0.d0 ) call stop_model("BallBerry: cs <= 0", 255)
      gsw = pspar%m*Anet*rh/cs + pspar%b
      if (gsw < pspar%b) gsw = pspar%b

      end function BallBerry


!-----------------------------------------------------------------------------
! Igor's codes
!-----------------------------------------------------------------------------

#ifdef USE_NR_SOLVER_FOR_FBB
      subroutine ci_cubic(ca,rh,gb,Pa,Rd,a1,f1,pspar,A)
!@sum ci_cubic Numerical solution for cubic equation of coupled
!@+   Ball-Berry/Farquhar stomatal conductance/photosynthesis.
!@+   Solves for Atot.
!@+   Version that uses Newton-Raphson solver
!@auth I.Aleinov
      implicit none
      real*8 :: ca              !Ambient air CO2 concentration (umol mol-1)
      real*8 :: rh              !Relative humidity
      real*8 :: gb              !Leaf boundary layer conductance of water vapor (mol m-2 s-1)
      real*8 :: Pa              !Pressure (Pa)
      real*8 :: Rd              !Leaf mitochondrial respiration (umol m-2 s-1)
      real*8 :: a1              !Coefficient in linear Farquhar equ.
      real*8 :: f1              !Coefficient in linear Farquhar equ.
      real*8, intent(out) :: A
      type(photosynthpar) :: pspar
      !----Local----
      real*8, parameter :: S_ATM=1.37d0  ! diffusivity ratio H2O/CO2 (atmosph.)
      real*8, parameter :: S_STOM=1.65d0 ! diffusivity ratio H2O/CO2 (stomatal)
      real*8 :: Ra, b, K, gamol, A_d_asymp
      real*8 :: x1, x2, xacc, x2tmp, x2save
      integer :: numit, counter=0
      save

      if ( pspar%reset_ci_cubic1 ) then
        pspar%reset_ci_cubic1 = .false.
cddd      if ( pspar%reset_ci_cubic1 == .false. ) then
cddd        write(777,*) counter, Ra, b, K, gamol, A_d_asymp,
cddd     &       x1, x2, xacc, x2tmp
cddd      endif
      Ra = 1/gb * S_ATM
      b = pspar%b / S_STOM
      K = pspar%m * rh / S_STOM
      gamol = pspar%Gammastar * 1.d06/Pa !Convert Pa to umol/mol
      A_d_asymp = - b*Ca / (K - b*Ra) ! asymptotic val of A from diffusion eq.

cddd      ! first check some special cases
cddd      if ( A_d_asymp >= 0.d0 ) then
cddd        ! this can happen only for very low humidity
cddd        ! probably should never happen in the real world, but if it does,
cddd        ! this case should be considered separately
cddd        !!print *,"!!! A_d_asymp >= 0.d0 !!!", A_d_asymp
cddd      !!!  A_d_asymp = -1.d30 !!! hack
cddd        !!print *,"K<b*Ra: m,rh,b,Ra:",pspar%m,rh,b,Ra
cddd        !!call stop_model("ci_cubic: rh too small ?",255)
cddd      endif

      !x1 = 0.d0
      x1 = -Rd
      x2save = ca/Ra
      x2tmp =  b*ca / (1.d0 - K + b*Ra)
      if( x2tmp > 0.d0 ) x2save = min( x2save, x2tmp )
      x2tmp = A_d_asymp
      if( x2tmp > 0.d0 ) x2save = min( x2save, x2tmp )
      x2save = x2save - .0000001d0
      x2 = min( x2, a1 - Rd)
      xacc = .0001d0
      !xacc = .01d0
cddd      if ( pspar%reset_ci_cubic1 == .false. ) then
cddd        write(778,*) counter, Ra, b, K, gamol, A_d_asymp,
cddd     &       x1, x2, xacc, x2tmp
cddd      endif
      endif
      x2 = min( x2save, a1 - Rd)
      A = rtsafe(A_eqn, x1,x2,xacc,  Ra, b, K, gamol,  ca, a1, f1, Rd
     &     , numit)
      !write(577,*) numit

      end subroutine ci_cubic

!---------------------------------------------------------------------------
      subroutine A_eqn(A, f, df,  Ra, b, K1, gamol,  ca, a1, f1, Rd )
!@sum Calculates the f coefficient in the coupled equ of photosynth/cond.
!@+   Igor, what is this solving for?? For f and df?
!@+   I.Aleinov
      real*8 A, f, df
      real*8 Ra, b, K1, gamol,  ca, a1, f1, Rd
      !---
      real*8, parameter :: S_ATM=1.37d0  ! diffusivity ratio H2O/CO2 (atmosph.)
      real*8, parameter :: S_STOM=1.65d0 ! diffusivity ratio H2O/CO2 (stomatal)
      real*8 cs, ci, dci
      real*8 byAKbcs, bycif1, K

      if ( A > 0.d0 ) then
        K = K1
      else
        K = 0.d0
      endif
      !write(579,*) "start A_eqn", A
      cs = ca - A*Ra
      byAKbcs = 1.d0/(A*K + b*cs)
      ci = cs * ( 1.d0 - A*byAKbcs )

      bycif1 = 1.d0/(ci+f1)
      f = A - ( a1*(ci-gamol)*bycif1 -Rd)

      dci = -Ra*( 1.d0 - A*byAKbcs )
     &     + cs*(- byAKbcs + A*byAKbcs*byAKbcs*(K-b*Ra) )
      df = 1 - a1*(f1+gamol)*bycif1*bycif1 * dci

      !write(579,*) "stop A_eqn", f, df

      end subroutine A_eqn

!---------------------------------------------------------------------------
      subroutine A_eqn_0(A, f, Ra, b, K1, gamol,  ca, a1, f1, Rd )
!@sum Calculates coefficients in equation for coupled photosynth/cond.
!@auth I.Aleinov
      real*8 A, f
      real*8 Ra, b, K1, gamol,  ca, a1, f1, Rd
      !---
      real*8, parameter :: S_ATM=1.37d0  ! diffusivity ratio H2O/CO2 (atmosph.)
      real*8, parameter :: S_STOM=1.65d0 ! diffusivity ratio H2O/CO2 (stomatal)
      real*8 cs, ci, dci
      real*8 byAKbcs, bycif1, K

      if ( A > 0.d0 ) then
        K = K1
      else
        K = 0.d0
      endif
      !write(579,*) "start A_eqn", A
      cs = ca - A*Ra
      byAKbcs = 1.d0/(A*K + b*cs)
      ci = cs * ( 1.d0 - A*byAKbcs )

      bycif1 = 1.d0/(ci+f1)
      f = A - ( a1*(ci-gamol)*bycif1 -Rd)

      !write(579,*) "stop A_eqn", f, df

      end subroutine A_eqn_0

!-----------------------------------------------------------

      FUNCTION rtsafe(funcd,x1,x2,xacc,  Ra, b, K, gamol,ca
     &     , a1, f1, Rd, numit )
!@sum Newton-Raphson solver (Numerical Recipes)
!@auth   I.Aleinov
      INTEGER MAXIT
      REAL*8 rtsafe,x1,x2,xacc
      real*8 Ra, b, K, gamol,  ca, a1, f1, Rd
      integer numit
      EXTERNAL funcd
      PARAMETER (MAXIT=100)
      INTEGER j
      REAL*8 df,dx,dxold,f,fh,fl,temp,xh,xl

      numit = 0

      call A_eqn_0(x1,fl,  Ra, b, K, gamol,  ca, a1, f1, Rd)
      call A_eqn_0(x2,fh,  Ra, b, K, gamol,  ca, a1, f1, Rd)
      if((fl.gt.0..and.fh.gt.0.).or.(fl.lt.0..and.fh.lt.0.)) then
        rtsafe = -1.d30
        return ! for now return 0
        !call stop_model('root must be bracketed in rtsafe',255)
      endif
      if(fl.eq.0.)then
        rtsafe=x1
        return
      else if(fh.eq.0.)then
        rtsafe=x2
        return
      else if(fl.lt.0.)then
        xl=x1
        xh=x2
      else
        xh=x1
        xl=x2
      endif
      !rtsafe=.5*(x1+x2)
      rtsafe=x1
      dxold=abs(x2-x1)
      dx=dxold
      call funcd(rtsafe,f,df,  Ra, b, K, gamol,  ca, a1, f1, Rd)
      do 11 j=1,MAXIT
        numit = j
        if(((rtsafe-xh)*df-f)*((rtsafe-xl)*df-f).ge.0..or. abs(2.*
     *f).gt.abs(dxold*df) ) then
          dxold=dx
          dx=0.5*(xh-xl)
          rtsafe=xl+dx
          if(xl.eq.rtsafe)return
        else
          dxold=dx
          dx=f/df
          temp=rtsafe
          rtsafe=rtsafe-dx
          if(temp.eq.rtsafe)return
        endif
        if(abs(dx).lt.xacc) return
        call funcd(rtsafe,f,df,  Ra, b, K, gamol,  ca, a1, f1, Rd)
        if(f.lt.0.) then
          xl=rtsafe
        else
          xh=rtsafe
        endif
11    continue
      call stop_model('rtsafe exceeding maximum iterations',255)
      return
      END FUNCTION rtsafe
#endif

!========================================================================
      end module BiomeE_PHY
