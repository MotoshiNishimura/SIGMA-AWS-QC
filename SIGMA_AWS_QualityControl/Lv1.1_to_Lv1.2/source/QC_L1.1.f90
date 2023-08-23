!=======================================================================
! This module summarizes the Quality Control process
! to create Level 1.2 dataset from Level 1.1 dataset of SIGMA-AWS observation data.
! Motoshi Nishimura @ National Institute of Polar Research
!=======================================================================

module QC_L11_12
  implicit none
  
  real(8), parameter :: atemp1A_max = 7.2d0
  real(8), parameter :: atemp2A_max = 7.2d0
  real(8), parameter :: atemp1B_max = 10.7d0
  real(8), parameter :: atemp1A_min = -49.9d0
  real(8), parameter :: atemp2A_min = -49.9d0
  real(8), parameter :: atemp1B_min = -40.5d0
  real(8), parameter :: ratio_nir = 0.5151d0     ! Fraction of the near-infrared spectral domain at the top of atmosphere (Wehrli, 1985)
  real(8), parameter :: UL_S_albedo = 0.95d0     ! upper limit of snow albedo (Aoki et al., 2003, 2011)  
  real(8), parameter :: UL_S_nir_albedo = 0.90d0 ! upper limit of snow near-infrared albedo (Aoki et al., 2003, 2011)
 
  real(8), parameter :: abstemp = 273.15d0         ! absolute temperature of 0 degC
  real(8), parameter :: eps = 0.98d0               ! snow emissivity [no dimension]
  real(8), parameter :: sigma = 5.67d-8            ! Stefan-Bolzmann constant[W/m^2/K^4]
  real(8), parameter :: mask = -9999.0d0           ! constant value repredenting data mask
  real(8), parameter :: manual_mask = -8888.0d0    ! constant value repredenting data manual mask
  real(8), parameter :: mask_threshold = -8887.0d0 ! mask or mannual mask threthold value

  contains

!==========================================
!>>> subroutines below
!==========================================
!-----------------------------------------------------------------------
! subroutines  : QC (Initial Control)
! version      : 1.0
! subject      : Initial control of SIGMA-AWS data
! editer       : Motoshi Nishimura @ National Institute of Polar Research
!-----------------------------------------------------------------------
subroutine QC_wind_speed (ndata, site_name, ws1, ws2)
  implicit none

  character(*), intent(in) :: site_name
  integer(8), intent(in)   :: ndata
  real(8), intent(inout)   :: ws1(ndata), ws2(ndata)
  
  real(8), parameter :: ws1A_max = 23.9d0
  real(8), parameter :: ws2A_max = 25.5d0
  real(8), parameter :: ws1B_max = 21.9d0

  integer(8) :: i

  do i = 1, ndata
    if (site_name == 'SIGMA-A') then
      if (ws1(i) < 0.0d0 .or. ws1(i) > ws1A_max + 15.0d0) ws1(i) = mask
      if (ws2(i) < 0.0d0 .or. ws2(i) > ws2A_max + 15.0d0) ws2(i) = mask
    end if
    if (site_name == 'SIGMA-B') then
      if (ws1(i) < 0.0d0 .or. ws1(i) > ws1B_max + 15.0d0) ws1(i) = mask
    end if
  end do

end subroutine QC_wind_speed

!----------------------------------------------------------------------------------------------------
subroutine QC_wind_direction (ndata, ws1, ws2, wd1, wd2)
  implicit none

  integer(8), intent(in) :: ndata
  real(8), intent(in)    :: ws1(ndata), ws2(ndata)
  real(8), intent(inout) :: wd1(ndata), wd2(ndata)

  integer(8) :: i

  do i = 1, ndata
     if (ws1(i) > 0.0d0) then
       if (wd1(i) <= 0.0d0 .or. wd1(i) > 360.0d0) wd1(i) = mask
     else if (ws1(i) == 0.0d0) then
       if (wd1(i) /= 0.0d0) then
         wd1(i) = mask
       end if
     end if

     if (ws2(i) > 0.0d0) then
       if (wd2(i) <= 0.0d0 .or. wd2(i) > 360.0d0) wd2(i) = mask
     else if (ws2(i) == 0.0d0) then
       if (wd2(i) /= 0.0d0) then
         wd2(i) = mask
       end if
     end if
  end do

end subroutine QC_wind_direction

!----------------------------------------------------------------------------------------------------
subroutine wind_direction_correction (ndata, site_name, wd1, wd2)
  implicit none

  character(*), intent(in) :: site_name
  integer(8), intent(in)   :: ndata
  real(8), intent(inout)   :: wd1(ndata), wd2(ndata)

  integer(8) :: i

  ! wind direction correction for SIGMA-A
  do i = 1, ndata
    if (site_name == 'SIGMA-A') then
      if (i >= 9334 .and. i <= 16940) then ! 22:00 on Jul. 24, 2013 - 20:00 on Jun. 6, 2014 20:00
        if (wd1(i) >= 0.0d0) wd1(i) = wd1(i) + 25.0d0
        if (wd2(i) >= 0.0d0) wd2(i) = wd2(i) + 25.0d0
        if (wd1(i) >= 360.0d0) wd1(i) = wd1(i) - 360.0d0
        if (wd2(i) >= 360.0d0) wd2(i) = wd2(i) - 360.0d0
      end if
      if (i >= 16941) then ! 21:00 on Jun. 6, 2014 - 
        if (wd1(i) >= 0.0d0) wd1(i) = wd1(i) + 10.0d0
        if (wd2(i) >= 0.0d0) wd2(i) = wd2(i) + 10.0d0
        if (wd1(i) >= 360.0d0) wd1(i) = wd1(i) - 360.0d0
        if (wd2(i) >= 360.0d0) wd2(i) = wd2(i) - 360.0d0     
      end if
    end if
  end do

end subroutine wind_direction_correction

!----------------------------------------------------------------------------------------------------
subroutine QC_air_temperature (ndata, site_name, atemp1, atemp2)
  implicit none

  character(*), intent(in) :: site_name
  integer(8), intent(in)   :: ndata
  real(8), intent(inout)   :: atemp1(ndata), atemp2(ndata)

  integer(8) :: i

  real(8), parameter :: atemp_margin = 10.0d0
  real(8), parameter :: specific_mask_value1 = -50.0d0
  real(8), parameter :: specific_mask_value2 = -11.8d0
  real(8), parameter :: specific_mask_value3 = 11.8d0

  do i = 1, ndata
    if (atemp1(i) == specific_mask_value1) atemp1(i) = mask

    if (site_name == 'SIGMA-A') then
      if (atemp1(i) < atemp1A_min - atemp_margin .or. atemp1(i) > atemp1A_max + atemp_margin) atemp1(i) = mask
      if (atemp2(i) < atemp2A_min - atemp_margin .or. atemp2(i) > atemp2A_max + atemp_margin) atemp2(i) = mask

      if (atemp2(i) == specific_mask_value1) atemp2(i) = mask
      if (i >= 45313) then ! 1:00 on Sep. 1, 2017 -
        if (atemp1(i) == specific_mask_value2) atemp1(i) = mask
        if (atemp1(i) == specific_mask_value3) atemp1(i) = mask
      end if
    end if

    if (site_name == 'SIGMA-B') then
      if (atemp1(i) < atemp1B_min - atemp_margin .or. atemp1(i) > atemp1B_max + atemp_margin) atemp1(i) = mask
    end if
  end do

end subroutine QC_air_temperature

!----------------------------------------------------------------------------------------------------
subroutine QC_air_humidity (ndata, arh1, arh2)
  implicit none

  integer(8), intent(in) :: ndata
  real(8), intent(inout) :: arh1(ndata), arh2(ndata)

  integer(8) :: i

  do i = 1, ndata
     if (arh1(i) <= 0.0d0 .or. arh1(i) > 100.0d0) arh1(i) = mask
     if (arh2(i) <= 0.0d0 .or. arh2(i) > 100.0d0) arh2(i) = mask
  end do

end subroutine QC_air_humidity

!----------------------------------------------------------------------------------------------------
subroutine QC_shortwave_radiation (ndata, site_name, SW_TOA, solz, solz_slope, swd, swu)
  implicit none

  character(*), intent(in) :: site_name
  integer(8), intent(in)   :: ndata
  real(8), intent(in)      :: SW_TOA(ndata)
  real(8), intent(in)      :: solz(ndata), solz_slope(ndata)
  real(8), intent(inout)   :: swd(ndata), swu(ndata)

  integer(8) :: i

  do i = 1, ndata
    if (swd(i) > maxval(SW_TOA)) swd(i) = mask
    if (swu(i) > maxval(SW_TOA) * UL_S_albedo) swu(i) = mask
    if (swd(i) < mask_threshold) swd(i) = mask
    if (swu(i) < mask_threshold) swu(i) = mask

    if (site_name == 'SIGMA-A') then
      if (swd(i) > mask_threshold .or. swu(i) > mask_threshold) then
        if (swd(i) < 0.0d0) then
          if (solz(i) < 90.0d0) then
            swd(i) = mask
          else
            swd(i) = 0.0d0
          end if
        end if
        if (swu(i) < 0.0d0) then
          if (solz(i) < 90.0d0) then
            swu(i) = mask
          else
            swu(i) = 0.0d0
          end if
        end if
      end if
    end if

    if (site_name == 'SIGMA-B') then
      if (swd(i) > mask_threshold .or. swu(i) > mask_threshold) then
        if (swd(i) < 0.0d0) then
          if (solz_slope(i) < 90.0d0) then
            swd(i) = mask
          else
            swd(i) = 0.0d0
          end if
        end if
        if (swu(i) < 0.0d0) then
          if (solz_slope(i) < 90.0d0) then
            swu(i) = mask
          else
            swu(i) = 0.0d0
          end if
        end if
      end if
    end if
  end do

end subroutine QC_shortwave_radiation

!----------------------------------------------------------------------------------------------------
subroutine QC_nir_infrared_radiation (ndata, site_name, SW_TOA, solz, solz_slope, nird, niru)
  implicit none

  character(*), intent(in) :: site_name
  integer(8), intent(in)   :: ndata
  real(8), intent(in)      :: SW_TOA(ndata)
  real(8), intent(in)      :: solz(ndata), solz_slope(ndata)
  real(8), intent(inout)   :: nird(ndata), niru(ndata)

  integer(8) :: i

  do i = 1, ndata
    if (nird(i) > maxval(SW_TOA) * ratio_nir) nird(i) = mask
    if (niru(i) > maxval(SW_TOA) * ratio_nir * UL_S_nir_albedo) niru(i) = mask
    if (nird(i) < mask_threshold) nird(i) = mask
    if (niru(i) < mask_threshold) niru(i) = mask

    if (site_name == 'SIGMA-A') then
      if (nird(i) > mask_threshold .or. niru(i) > mask_threshold) then
        if (nird(i) < 0.0d0) then
          if (solz(i) < 90.0d0) then
            nird(i) = mask
          else
            nird(i) = 0.0d0
          end if
        end if
        if (niru(i) < 0.0d0) then
          if (solz(i) < 90.0d0) then
            niru(i) = mask
          else
            niru(i) = 0.0d0
          end if
        end if
      end if
    end if

    if (site_name == 'SIGMA-B') then
      if (nird(i) > mask_threshold .or. niru(i) > mask_threshold) then
        if (nird(i) < 0.0d0) then
          if (solz_slope(i) < 90.0d0) then
            nird(i) = mask
          else
            nird(i) = 0.0d0
          end if
        end if
        if (niru(i) < 0.0d0) then
          if (solz_slope(i) < 90.0d0) then
            niru(i) = mask
          else
            niru(i) = 0.0d0
          end if
        end if
      end if
    end if
  end do

end subroutine QC_nir_infrared_radiation

!----------------------------------------------------------------------------------------------------
subroutine QC_longwave_radiation (ndata, site_name, lwd, lwu)
  implicit none

  character(*), intent(in) :: site_name
  integer(8), intent(in)   :: ndata
  real(8), intent(inout)   :: lwd(ndata), lwu(ndata)

  integer(8) :: i
  real(8), parameter :: sst_max = 10.0d0 ! max snow temperature including observation error etc.
  real(8), parameter :: emiss_cloud = 1.0d0 ! atmospheric emissivity in cloudy condition

  do i = 1, ndata
    if (site_name == 'SIGMA-A') then
      if (lwd(i) < 0.0d0 .or. lwd(i) > (sigma * emiss_cloud * (atemp2A_max + abstemp) ** 4.0d0)) lwd(i) = mask

      if (i >= 45313) then ! 1:00 on Sep. 1, 2017 -
        if (lwd(i) == 140.6d0) lwd(i) = mask
        
        if (lwu(i) == 128.5d0) lwu(i) = mask
        if (lwu(i) == 130.0d0) lwu(i) = mask
        if (lwu(i) == 131.5d0) lwu(i) = mask
        if (lwu(i) == 133.0d0) lwu(i) = mask
        if (lwu(i) == 134.5d0) lwu(i) = mask
        if (lwu(i) == 137.6d0) lwu(i) = mask
        if (lwu(i) == 139.1d0) lwu(i) = mask
        if (lwu(i) == 140.6d0) lwu(i) = mask
        if (lwu(i) == 264.5d0) lwu(i) = mask
      end if

      if (i >= 48141 .and. i <= 51393) then
        if (lwd(i) == 123.9d0) lwd(i) = mask
        if (lwd(i) == 125.4d0) lwd(i) = mask
        if (lwd(i) == 126.9d0) lwd(i) = mask
        if (lwd(i) == 134.3d0) lwd(i) = mask
        if (lwd(i) == 138.1d0) lwd(i) = mask
        if (lwd(i) == 139.3d0) lwd(i) = mask
        if (lwd(i) == 264.5d0) lwd(i) = mask
        
        if (lwu(i) == 123.9d0) lwu(i) = mask
        if (lwu(i) == 125.4d0) lwu(i) = mask
        if (lwu(i) == 126.9d0) lwu(i) = mask
      end if
      if (i >= 55009 .and. i <= 59160) then
        if (lwd(i) == 123.9d0) lwd(i) = mask
        if (lwd(i) == 125.4d0) lwd(i) = mask
        if (lwd(i) == 126.9d0) lwd(i) = mask
        if (lwd(i) == 138.1d0) lwd(i) = mask
        if (lwd(i) == 139.3d0) lwd(i) = mask
        if (lwd(i) == 262.0d0) lwd(i) = mask
        if (lwd(i) == 263.3d0) lwd(i) = mask
        if (lwd(i) == 264.5d0) lwd(i) = mask
        
        if (lwu(i) == 126.9d0) lwu(i) = mask
      end if
      if (i >= 65107 .and. i <= 67944) then 
        if (lwd(i) == 123.9d0) lwd(i) = mask
        if (lwd(i) == 125.4d0) lwd(i) = mask
        if (lwd(i) == 126.9d0) lwd(i) = mask
        if (lwd(i) == 138.1d0) lwd(i) = mask
        if (lwd(i) == 139.3d0) lwd(i) = mask
        if (lwd(i) == 139.4d0) lwd(i) = mask
        if (lwd(i) == 263.3d0) lwd(i) = mask
      end if
    end if
    
    if (site_name == 'SIGMA-B') then
      if (lwd(i) < 0.0d0) lwd(i) = mask
      if (lwd(i) > (sigma * emiss_cloud * (atemp1B_max + abstemp) ** 4.0d0)) lwd(i) = mask
    end if

    if (lwu(i) < 0.0d0) lwu(i) = mask
    if (lwu(i) > (sigma * eps * (sst_max + abstemp) ** 4.0d0)) lwu(i) = mask
     
    if (lwd(i) > -9998.0d0 .and. lwu(i) > -9998.0d0) then
      if (abs(lwd(i) - lwu(i)) < 1.0d0) then
        lwd(i) = -9998.0d0
        lwu(i) = -9998.0d0
      end if
    end if
  end do

end subroutine QC_longwave_radiation

!----------------------------------------------------------------------------------------------------
subroutine QC_snow_depth (ndata, site_name, atemp, sd)
  implicit none

  character(*), intent(in) :: site_name
  integer(8), intent(in)   :: ndata
  real(8), intent(inout)   :: atemp(ndata)
  real(8), intent(inout)   :: sd(ndata)

  integer(8) :: i
  real(8)    :: sd_diff0
  
  real(8), parameter :: specific_mask_value1 = -500.0d0
  real(8), parameter :: specific_mask_value2 = 193.0d0
  real(8), parameter :: specific_mask_value3 = 300.0d0
  real(8), parameter :: pole_length1 = 130.0d0
  real(8), parameter :: pole_length2 = 147.0d0
  real(8), parameter :: pole_length3 = 130.0d0
  real(8), parameter :: set_height_diff1 = 45.1d0
  real(8), parameter :: pole_adjustment_length = 0.8d0
  real(8), parameter :: sd_median1 = 58.5d0
  real(8), parameter :: sd_median2 = 163.9d0
  real(8), parameter :: sd_median3 = 238.0d0
  real(8), parameter :: sd_median4 = 513.0d0
  real(8), parameter :: margin_sd = 100.0d0
  real(8), parameter :: margin_coefficient = 1.5d0

  do i = 1, ndata
    if (site_name == 'SIGMA-A') then
      if (sd(i) == specific_mask_value1) sd(i) = mask
      if (sd(i) == specific_mask_value2) sd(i) = mask
      if (sd(i) == specific_mask_value3) sd(i) = mask

      sd_diff0 = sd(323) - sd(347)
      if (i >= 326 .and. i <= 9330) then
        if (sd(i) > -200.0d0) sd(i) = sd(i) + sd_diff0
      end if  
    
      ! maintenance record
      if (sd(i) > mask_threshold) then
        ! Jul. 24, 2013
        if (i >= 9330 .and. i <= 9346) sd(i) = mask ! maintenance time
        if (i >= 9347 .and. i <= ndata) sd(i) = sd(i) + pole_length1 ! 11:00 on Jul. 25, 2013 -
        ! Jun. 6, 2014
        if (i >= 16931 .and. i <= 16939) sd(i) = mask ! maintenance time
        if (i >= 16940 .and. i <= ndata)  sd(i) = sd(i) + pole_length2 ! 20:00 on Jun. 6, 2014 -
        ! May 26, 2017
        if (i >= 42889 .and. i <= 42977) sd(i) = mask ! maintenance time
        if (i >= 42978 .and. i <= ndata)  sd(i) = sd(i) + pole_length3 ! 18:00 on May 26, 2017 -
      end if

      if (i <= 9346) then ! - 11:00 on Jul. 25, 2013
        if (sd(i) < sd_median1 - margin_sd .or. sd(i) > sd_median1 + margin_sd) sd(i) = mask
      end if
      if (i >= 9347 .and. i <= 16939) then ! 11:00 on Jul. 25, 2013 - 20:00 on Jun. 6, 2014
        if (sd(i) < sd_median2 - margin_sd .or. sd(i) > sd_median2 + margin_sd) sd(i) = mask
      end if
      if (i >= 16940 .and. i <= 42976) then ! 20:00 on Jun. 6, 2014 - 18:00 on May 26, 2017
        if (sd(i) < sd_median3 - margin_sd * 1.5d0 .or. sd(i) > sd_median3 + margin_sd * 1.5d0) sd(i) = mask
      end if
      if (i >= 42977 .and. i <= ndata) then ! 18:00 on May 26, 2017 -
        if (sd(i) < sd_median4 - margin_sd * 1.5d0 .or. sd(i) > sd_median4 + margin_sd * 1.5d0) sd(i) = mask
      end if
    end if

    if (site_name == 'SIGMA-B') then
      if (sd(i) == -112.0d0) sd(i) = mask
      if (sd(i) == -111.0d0) sd(i) = mask
      if (sd(i) == -110.0d0) sd(i) = mask
      if (sd(i) == -109.0d0) sd(i) = mask
      if (sd(i) == 123.0d0) sd(i) = mask
      if (sd(i) == 124.0d0) sd(i) = mask
      if (sd(i) == 125.0d0) sd(i) = mask    
      if (sd(i) == 126.0d0) sd(i) = mask
      if (sd(i) == 256.0d0) sd(i) = mask
      if (sd(i) == 300.0d0) sd(i) = mask

      ! 1:00 on Apr. 26, 2021 - 0:00 on Jul. 8, 2022
      if (i >= 76849 .and. i<= 87360) then
        if (sd(i) >= 123.0d0 .and. sd(i) <= 127.0d0) sd(i) = mask
      end if

      ! maintenance record
      ! Jul. 10, 2022
      if (i >= 87422 .and. i <= 87426) sd(i) = mask ! maintenance time

      if (sd(i) > mask_threshold) then
        ! 14:00 on Jul. 10, 2023 - 14:00 on Jul. 25, 2022
        if (i >= 87422 .and. i < 87783) sd(i) = sd(i) - set_height_diff1

        ! 15:00 on Jul. 25, 2022
        if (i >= 87783) then
          sd(i) = 300.0d0 - (sd(i) + set_height_diff1)! data output change due to an Argos algorithm change
          if (atemp(i) > mask_threshold) sd(i) = sd(i) * ((atemp(i) + abstemp) / abstemp) ** 0.5d0
        end if
        
        ! 22:00 on Aug. 7, 2022 -
        if (i >= 88102) sd(i) = sd(i) - pole_adjustment_length 
      end if
    end if  
  end do

end subroutine QC_snow_depth

!----------------------------------------------------------------------------------------------------
subroutine QC_air_press (ndata, site_name, press)
  implicit none

  character(*), intent(in) :: site_name
  integer(8), intent(in)   :: ndata
  real(8), intent(inout)   :: press(ndata)

  integer(8) :: i
  
  real(8), parameter :: press_meanA = 833.06d0
  real(8), parameter :: press_meanB = 894.19d0
  real(8), parameter :: press_margin = 100.0d0

  do i = 1, ndata
    if (site_name == 'SIGMA-A') then
      if (press(i) < press_meanA - press_margin .or. press(i) > press_meanA + press_margin) press(i) = mask
    end if
    if (site_name == 'SIGMA-B') then
      if (press(i) < press_meanB - press_margin .or. press(i) > press_meanB + press_margin) press(i) = mask
    end if
  end do

end subroutine QC_air_press

!----------------------------------------------------------------------------------------------------
subroutine QC_snow_temperature (ndata, st1, st2, st3, st4, st5, st6)
  implicit none

  integer(8), intent(in) :: ndata
  real(8), intent(inout) :: st1(ndata), st2(ndata), st3(ndata), st4(ndata), st5(ndata), st6(ndata)

  integer(8) :: i

  real(8), parameter :: lowest_atemp = atemp1A_min
  real(8), parameter :: lower_margin = 10.0d0
  real(8), parameter :: upper_thershold = 0.2d0
  real(8), parameter :: specific_mask_value1 = -50.0d0
  real(8), parameter :: specific_mask_value2 = -11.8d0
  
  do i = 1, ndata
    ! 3 days after maintenance: 13:00 on Jul. 26, 2013 - 12:00 on Jul. 29, 2013
    if (i >= 9373 .and. i <= 9444) then
       st3(i) = mask
       st4(i) = mask
    end if
    
    if (st1(i) < lowest_atemp - lower_margin .or. st1(i) > upper_thershold) st1(i) = mask
    if (st2(i) < lowest_atemp - lower_margin .or. st2(i) > upper_thershold) st2(i) = mask
    if (st3(i) < lowest_atemp - lower_margin .or. st3(i) > upper_thershold) st3(i) = mask
    if (st4(i) < lowest_atemp - lower_margin .or. st4(i) > upper_thershold) st4(i) = mask
    if (st5(i) < lowest_atemp - lower_margin .or. st5(i) > upper_thershold) st5(i) = mask
    if (st6(i) < lowest_atemp - lower_margin .or. st6(i) > upper_thershold) st6(i) = mask

    if (st1(i) == specific_mask_value1 .or. st1(i) == specific_mask_value2) st1(i) = mask
    if (st2(i) == specific_mask_value1 .or. st2(i) == specific_mask_value2) st2(i) = mask
    if (st3(i) == specific_mask_value1 .or. st3(i) == specific_mask_value2) st3(i) = mask
    if (st4(i) == specific_mask_value1 .or. st4(i) == specific_mask_value2) st4(i) = mask
    if (st5(i) == specific_mask_value1 .or. st5(i) == specific_mask_value2) st5(i) = mask
    if (st6(i) == specific_mask_value1 .or. st6(i) == specific_mask_value2) st6(i) = mask

  end do

end subroutine QC_snow_temperature
!==========================================
!>>> subroutines end
!==========================================
end module QC_L11_12


!>>> Revision records
! ver.1.0    (date: Sep. 8, 2022) subroutine snow_temperature was added
! ver.2.0    (date: Sep. 9, 2022) Subroutines for the all elements have been completed to be added
