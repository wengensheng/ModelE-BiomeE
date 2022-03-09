#include "rundeck_opts.h"
      module ent
!@sum replaced by BiomeE processes, Weng 2018-11-05
! 2019-06-24 for test

      use ent_const
      use ent_types
      use BiomeE_VEG, only: ijdebug
      implicit none
      private
      save

      public BiomeE_integrate
      public update_veg_structure
      public ent_biophysics
      !public BiomeE_update_veg

      contains
      !*********************************************************************
      subroutine BiomeE_integrate(dtsec, ecp, update_day, config) ! BiomeE_integrate
!@sum Main routine to control BiomeE plant physiology and soil biogeochemistry.

      use cohorts
      use patches
      use BiomeE_PHY, only : BiomeE_CanopyModel,clim_stats
      use BiomeE_VEG, only : BiomeE_soil_BGC
      use entcells, only : summarize_entcell, entcell_print

      implicit none

      real*8 :: dtsec  !dt in seconds
      real*8 :: runyrs  ! years of model run, Weng for output
      !type(timestruct),pointer :: tt !Time in year.fraction, Greenwich Mean Time
      type(entcelltype) :: ecp
      logical :: update_day
      type(ent_config) :: config

      !-----local--------
      integer :: patchnum, npatch
      type(patch),pointer :: pp, pp_tmp
      type(cohort), pointer :: cop ! Weng
      integer :: n,i

!     Start here!
      ijdebug = ecp%ijdebug
      ecp%steps_d = ecp%steps_d + 1

      call clim_stats(dtsec,ecp,config,update_day)
      !* Dynamic vegetation, Weng, 2017-02-07, 2018-4-16
      ! 'update_day' means the start of a new day
      if (update_day .and. ecp%steps_d > 12) then
        call update_veg_structure(ecp, config) ! call BiomeE_update_veg(ecp)

        ecp%steps_d = 0
      endif

      ! Plant physiological and Soil biogeochemical processes
      pp => ecp%oldest
      do while (ASSOCIATED(pp))
        call BiomeE_CanopyModel(pp,dtsec) ! Physiological processes
        call BiomeE_soil_BGC(dtsec,pp,update_day)
        !call soil_bgc_CN(dtsec,pp,update_day)  ! Soil biogeochemical processes
        pp%CO2flux = -pp%NPP + pp%Soil_resp
        pp%age = pp%age + dtsec

        pp => pp%younger
      end do ! all patches hourly update

      ! Summarize
      call summarize_entcell(ecp)

!*********** DIAGNOSTICS FOR PLOTTING ********************!
#ifdef ENT_STANDALONE_DIAG
      patchnum = 0
      pp => ecp%oldest
      do while (ASSOCIATED(pp))
        patchnum = patchnum + 1
        runyrs = pp%age/(365.0*24.0*3600.0)
        if(runyrs >= 0.0 .and. runyrs < 3.0)then
          cop => pp%tallest
          write(798,'(4(I3,","))')pp%nCohorts
          do while(ASSOCIATED(cop))
            write(798,'(4(I5,","),9(F10.5,","))')
     &                cop%chID,pp%doy,cop%pft,cop%layer,
     &                runyrs,cop%n,cop%crownarea,
     &                cop%GPP*1000.0*3600.0,
     &                cop%NPP*1000.0*3600.0,
     &                cop%R_auto*1000.0*3600.0,
     &                cop%LAI,cop%phenostatus
                 cop  => cop%shorter
          enddo

          write(799,'(1(I5,","),30(E12.6,","))') !Fluxes are positive up.
     &     pp%doy,runyrs,
     &     pp%GPP*1000.*3600.,   pp%NPP*1000.*3600.,
     &     pp%R_auto*1000.*3600.,pp%Soil_resp*1000.*3600.,
     &     pp%LAI,pp%C_fol,pp%C_w,pp%C_froot,
     &     pp%TRANS_SW,
     &     pp%Ci, pp%GCANOPY,pp%IPP,
     &     pp%c_growth,pp%N_up,pp%betad

        endif
        pp => pp%younger
      end do ! all patches hourly update
#endif
!*********************************************************!
      call debug_diags(ecp) ! Weng, 2020-09-04

      end subroutine BiomeE_integrate

!============================================================
      subroutine update_veg_structure(ecp, config) ! BiomeE_update_veg
!@sum Update vegetation structure at the end of day

      use BiomeE_PHY,only : BiomeE_Phenology
      use entcells,only : summarize_entcell, entcell_print
      use BiomeE_VEG,only: BiomeE_FullDemography,BiomeE_SingleTree

      implicit none
      type(entcelltype) :: ecp
      type(ent_config) :: config
      !-----local--------
      type(patch),pointer :: pp

      !config%do_demographic_EW = .false.

      !* Loop through patches
      pp => ecp%oldest
      do while (ASSOCIATED(pp))

        ! update phenology status and vegetation structure (BiomeE)
        call BiomeE_Phenology(pp)
        ! update vegetation
        if(config%do_demographic_EW)then
          call BiomeE_FullDemography(pp)
        else
          call BiomeE_SingleTree(pp)
        endif
        pp => pp%younger
      enddo

      call summarize_entcell(ecp)

      ecp%daylength(1) = ecp%daylength(2)
      ecp%daylength(2) = 0.d0

      end subroutine update_veg_structure

!*********************************************************************
      subroutine debug_diags(ecp)
      use ent_debug_mod
      type(entcelltype), intent(in) :: ecp
      !---
      type(patch),pointer :: pp
      type(cohort), pointer :: cop

      !------local var ---------
      real*8 :: s_per_yr
      real*8 :: defacc_scale= 1.d0/(3600*24*1000.d0)  !!! hack
      real*8 :: area, n, Npatch,A_PFT(16)
      integer i, k, pft
      logical :: Vegetated

      Vegetated = .False.
      s_per_yr = 3600.d0 * 24 * 365
      A_PFT(:)  = 0.d0
      ent_d%vf(:) =  0.d0 
      ent_d%LAI(:) = 0.d0 
      ent_d%GPP(:) = 0.d0
      ent_d%R_auto(:) = 0.d0

      ent_d%total(:) = 0.d0
      ent_d%C_lab(:) = 0.d0
      ent_d%C_fol(:) = 0.d0
      ent_d%C_wd(:) = 0.d0
      ent_d%P_h(:) = 0.d0
      ent_d%C_froot(:) = 0.d0
      ent_d%P_dbh(:) = 0.d0
      ent_d%phenofactor(:) = 0.d0
      ent_d%betad(:) = 0.d0
      ent_d%C_soil(:) = 0.d0
      ent_d%Resp_soil(:) = 0.d0

      pp => ecp%oldest
      do while (associated(pp))
        area = pp%area
        pft = 0
        cop => pp%tallest
        do while (associated(cop))
          pft = cop%pft
          n = cop%n
          Vegetated = .True.
          Npatch = cop%n * pp%area ! n*1.d-3 * defacc_scale      ! was : n*area*1.d-3
          !set values for debugging
          A_PFT(pft)  = A_PFT(pft)  + pp%area
          ent_d%vf(pft) = ent_d%vf(pft) + pp%area  ! /(3600*24*1000.d0)
          ent_d%LAI(pft) = ent_d%LAI(pft) + cop%leafarea*Npatch ! /(3600*24*1000.d0)
          ent_d%GPP(pft) = ent_d%GPP(pft) + cop%GPP*Npatch
          ent_d%R_auto(pft) = ent_d%R_auto(pft) + cop%R_auto*Npatch

          ent_d%C_lab(pft) = ent_d%C_lab(pft) + cop%C_lab*Npatch
          ent_d%C_fol(pft) = ent_d%C_fol(pft) + cop%C_fol*Npatch
          ent_d%C_wd(pft) = ent_d%C_wd(pft)
     &                    + (cop%C_sw+cop%C_hw)*Npatch
          ent_d%P_h(pft) = Max(ent_d%P_h(pft), cop%h) 
          ent_d%C_froot(pft) = ent_d%C_froot(pft) + cop%C_froot*Npatch
          ent_d%P_dbh(pft) = Max(ent_d%P_dbh(pft), cop%dbh) 
          ent_d%phenofactor(pft) = ent_d%phenofactor(pft)
     &         + cop%phenofactor * pp%area !*defacc_Npatch
          ent_d%betad(pft) = ent_d%betad(pft)
     &         + cop%stressH2O * pp%area ! *defacc_scale

          cop => cop%shorter
        end do

        !!! patch-level variables
        if ( pft > 0) then      ! skip cells with no vegetation
          do i=1,N_CASA_LAYERS
            do k=(NLIVE+1),NPOOLS
              ent_d%C_soil(pft)=ent_d%C_soil(pft)
     &             + pp%Tpool(CARBON,k,i) * area ! * defacc_scale ! *area *1.d-3
            enddo
          enddo
          ent_d%Resp_soil(pft)=ent_d%Resp_soil(pft)+pp%Soil_resp*area
        endif ! PFT > 0

        pp => pp%younger
      end do ! all patches

      ! Area-normalized ! Weng, 1/1/2021
      do pft=1, 16
         if(A_PFT(pft)>0.d0)then
          ent_d%vf(pft) = ent_d%vf(pft)  /A_PFT(pft)
          ent_d%LAI(pft) = ent_d%LAI(pft)/A_PFT(pft)
          ent_d%GPP(pft) = ent_d%GPP(pft)/A_PFT(pft)       * s_per_yr
          ent_d%R_auto(pft) = ent_d%R_auto(pft)/A_PFT(pft) * s_per_yr

          ent_d%C_lab(pft) = ent_d%C_lab(pft)/A_PFT(pft)
          ent_d%C_fol(pft) = ent_d%C_fol(pft)/A_PFT(pft)
          ent_d%C_wd(pft) = ent_d%C_wd(pft)/A_PFT(pft)
          ent_d%C_froot(pft) = ent_d%C_froot(pft)/A_PFT(pft)
          ent_d%phenofactor(pft) = ent_d%phenofactor(pft) /A_PFT(pft)
          ent_d%betad(pft) = ent_d%betad(pft)/A_PFT(pft)
          ent_d%C_soil(pft) = ent_d%C_soil(pft)/A_PFT(pft)
          ent_d%Resp_soil(pft)=ent_d%Resp_soil(pft)/A_PFT(pft)*s_per_yr
         endif
      enddo
      if(Vegetated)then ! ra017
!        ent_d%total(:) = ent_d%C_lab(:) + ent_d%C_fol(:)
!     &     + ent_d%C_wd(:) ! + ent_d%C_hw(:)
!     &     + ent_d%C_froot(:) ! + ent_d%C_croot(:)
!!     &     + ent_d%C_soil(:)
        ent_d%total(1)  = ecp%LAI
        ent_d%total(2)  = ecp%h
        ent_d%total(3)  = ecp%GPP       * s_per_yr
        ent_d%total(4)  = ecp%R_auto    * s_per_yr
        ent_d%total(5)  = ecp%Soil_resp * s_per_yr

        ent_d%total(6)  = ecp%C_lab
        ent_d%total(7)  = ecp%C_fol
        ent_d%total(8)  = ecp%C_froot
        ent_d%total(9)  = ecp%Ctot
        ent_d%total(10) = sum(ecp%Tpool(1,4:12,1))

        ent_d%total(11) = ecp%N_lab
        ent_d%total(12) = ecp%N_fol
        ent_d%total(13) = ecp%N_froot
        ent_d%total(14) = ecp%Ntot
        ent_d%total(15) = sum(ecp%Tpool(2,4:12,1))
        ent_d%total(16) = ecp%area
      endif
!!!!!! hack hack hack
!!!!!! using phenofactor to output snow albedo
#ifdef ENT_ALBEDO_DIAGS_HACK
      ent_d%phenofactor(:) = 0.d0
      ent_d%phenofactor(1:6) = ecp%snow_albedo(1:6,1)*defacc_scale
      ent_d%phenofactor(7:12) = ecp%snow_albedo(1:6,2)*defacc_scale
#endif

      end subroutine debug_diags

      !*********************************************************************
      subroutine ent_biophysics(dtsec, ecp, config)
!@sum  Photosynthesis CO2 uptake and conductance of water vapor.
!@+    This optional routine may be used by a GCM or land surface driver
!@+    to calculate only surface conductance without biogeochemistry.
!@+    Useful for implicit schemes, but not used in the default setup.
!@+    If do_soilresp, then also  soil respiration for net CO2 fluxes.
      implicit none
      real*8 :: dtsec  !dt in seconds
      type(entcelltype) :: ecp
      type(ent_config) :: config
      !---Local--------
      type(patch),pointer :: pp
      integer :: patchnum

      call stop_model("ent_biophysics: obsolete code",255)

      end subroutine ent_biophysics
      !*********************************************************************

      end module ent

