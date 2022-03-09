#include "rundeck_opts.h"

      module ent_prescribed_updates
!@sum Routines for updating prescribed (Matthews) vegetation. 
!@+   These routines work on the entcell level or lower.

      use ent_types

      implicit none
      private

      public entcell_vegupdate

      contains

!************************************************************************

      subroutine entcell_vegupdate(ecp, hemi, jday
     &     ,do_giss_phenology, do_giss_lai, do_giss_albedo !,mixed_veg
     &     ,do_update_crops
     &     ,laidata, hdata, albedodata, cropsdata, init )
!@sum Main routine for prescribed vegetation structure updates
!@+   at the entcell level (and down). 
!@+   DAILY TIME STEP ASSUMED. Coordinate-dependent
!@+   All var parameters except entcell are optional. 
!@+   Var arrays have pointer attribute to provide a way to tell the 
!@+   program that an argument is actually optional and not missing
!@+   (see how it is used in ent_prescribe_vegupdate)
      use ent_pfts, only: COVEROFFSET
      use entcells,only : summarize_entcell
      use ent_prescr_veg, only : Matthews_calc_veg_albedo
      implicit none
      type(entcelltype),pointer :: ecp
      integer,intent(in) :: jday
      integer,intent(in) :: hemi
      logical, intent(in) :: do_giss_phenology
      logical, intent(in) :: do_giss_lai
      logical, intent(in) :: do_giss_albedo
!      logical, intent(in) :: mixed_veg
      logical, intent(in) :: do_update_crops
      real*8,  pointer :: laidata(:)  !Array of length N_PFT
      real*8,  pointer :: hdata(:)  !Array of length N_PFT
      real*8,  pointer :: albedodata(:,:)
      real*8,  pointer :: cropsdata
      logical, intent(in) :: init
      !----Local------
      type(patch),pointer :: pp

      !* 1. Update crops to get right patch/cover distribution.
      !*     NOTE:  CARBON CONSERVATION NEEDS TO BE CALCULATED FOR CHANGING VEG/CROP COVER ##
      !* 2. Update height to get any height growth (with GISS veg, height is static)
      !* 3. Update LAI, and accumulate litter from new LAI and growth/senescence.
      !*       Cohort litter is accumulated to the patch level.
      !*    3a. If external LAI, then litter is calculated based on that external LAI change.
      !*    3b. If GISS prescribed LAI, then new LAI is calculated, and then litter.
      !* 4. Update albedo based on new vegetation structure.

      !write(*,*) 'Got here before calculating giss albedo' !## NK DEBUG ##
      !* ALBEDO *!
      if (do_giss_albedo) then 
         if ( associated(albedodata) ) then
            call entcell_update_albedo(ecp, albedodata)
            !print *, "update albedodata from array"
         else
            pp => ecp%oldest
            do while (ASSOCIATED(pp))
               ! update if have vegetation or not prognostic albedo
               if ( ASSOCIATED(pp%tallest).and.do_giss_albedo )
     &              call Matthews_calc_veg_albedo(hemi,
     &              pp%tallest%pft+COVEROFFSET, jday, pp%albedo)
               pp => pp%younger
            end do
         !print *, "update albedodata hemi"
         endif
      endif !else ACTS already updates albedo: must have RAD_MODEL=GORT
      
      call summarize_entcell(ecp)

!====================================================================
      end subroutine entcell_vegupdate

      subroutine entcell_update_albedo( ecp,
     i    albedodata)
!@sum sets prescribed albedo in vegetated patches of the cell (skips bare soil)
!@+   This subroutine assumes one cohort per patch !!!
      type(entcelltype) :: ecp
      real*8,intent(in) :: albedodata(N_BANDS,N_PFT) !@var albedo for all PFTs 
      !-----Local---------
      type(patch), pointer :: pp  !@var p current patch

      pp => ecp%oldest      
      do while ( associated(pp) )
        ! update albedo of vegetated patches only
        if ( associated(pp%tallest) ) then ! have vegetation
          if ( pp%tallest%pft > N_PFT .or.  pp%tallest%pft < 1 )
     &         call stop_model("entcell_update_albedo: bad pft",255)
          pp%albedo(1:N_BANDS) = albedodata(1:N_BANDS, pp%tallest%pft)
        endif
        pp => pp%younger
      enddo

      end subroutine entcell_update_albedo

!******************************************************************************

      end module ent_prescribed_updates
