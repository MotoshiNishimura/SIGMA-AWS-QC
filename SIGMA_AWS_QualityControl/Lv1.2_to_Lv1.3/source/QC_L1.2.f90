!=======================================================================
! This module summarizes the Quality Control process
! to create Level 1.3 dataset from Level 1.2 dataset of SIGMA-AWS observation data.
! Motoshi Nishimura @ National Institute of Polar Research
!=======================================================================
module QC_L12_13
  implicit none

  real(8), parameter :: ratio_nir = 0.5151d0 ! Fraction of the near-infrared spectral domain at the top of atmosphere (Wehrli, 1985)
  real(8), parameter :: UL_S_albedo = 0.95d0 ! upper limit of snow albedo (Aoki et al., 2003, 2011)
  real(8), parameter :: LL_DS_albedo = 0.60d0 ! lower limit of dry snow albedo (Aoki et al., 2003, 2011)
  real(8), parameter :: LL_WS_albedo = 0.40d0 ! lower limit of wet snow albedo (Aoki et al., 2003, 2011)
  real(8), parameter :: LL_DI_albedo = 0.10d0 ! lower limit of dark ice albedo (Aoki et al., 2003, 2011)
  real(8), parameter :: UL_S_nir_albedo = 0.90d0 ! upper limit of snow near-infrared albedo (Aoki et al., 2003, 2011)
  real(8), parameter :: LL_DS_nir_albedo = 0.50d0 ! lower limit of dry snow near-infrared albedo (Aoki et al., 2003, 2011)
  real(8), parameter :: LL_WS_nir_albedo = 0.30d0 ! lower limit of wet snow near-infrared albedo (Aoki et al., 2003, 2011)
  real(8), parameter :: co_transA = 0.881d0 ! Maximum atmospheric transimissivity at SIGMA-A site
  real(8), parameter :: co_transB = 0.872d0 ! Maximum atmospheric transimissivity at SIGMA-B site

  real(8), parameter :: mask = -9999.0d0  !constant value repredenting data mask
  real(8), parameter :: radiation_sensor_mask = -9998.0d0 ! radiation error data when the sensor is covered with snow or ice  
  real(8), parameter :: manual_mask = -8888.0d0  !constant value repredenting data manual mask
  real(8), parameter :: mask_threshold = -8887.0d0  !mask or mannual mask threthold value

  contains

!==========================================
!>>> subroutines below
!==========================================
!-----------------------------------------------------------------------
! subroutines  : QC (Secondary Control)
! version      : 1.0
! subject      : Secondary control of SIGMA-AWS data
! editer       : Motoshi Nishimura @ National Institute of Polar Research
!-----------------------------------------------------------------------
subroutine QC_wind_speed (ndata, wd1, wd2, ws1, ws2)
  implicit none

  integer(8), intent(in) :: ndata
  real(8), intent(inout) :: wd1(ndata), wd2(ndata)
  real(8), intent(inout) :: ws1(ndata), ws2(ndata)

  integer(8) :: i, j

  do i = 1, ndata
    ! Detecting of icing and freezing on wind sensor
    if (i >= 6) then
      if (ws1(i-5) == 0.0d0 .and. ws1(i-4) == 0.0d0 .and. ws1(i-3) == 0.0d0 .and. &
&         ws1(i-2) == 0.0d0 .and. ws1(i-1) == 0.0d0 .and. ws1(i) == 0.0d0) then
        do j = 0, 5
          ws1(i-j) = mask 
          wd1(i-j) = mask
        end do
      end if

      if (ws2(i-5) == 0.0d0 .and. ws2(i-4) == 0.0d0 .and. ws2(i-3) == 0.0d0 .and. &
&         ws2(i-2) == 0.0d0 .and. ws2(i-1) == 0.0d0 .and. ws2(i) == 0.0d0) then
        do j = 0, 5
          ws2(i-j) = mask 
          wd2(i-j) = mask
        end do
      end if
    end if
  end do

end subroutine QC_wind_speed

!----------------------------------------------------------------------------------------------------
subroutine QC_air_temperature (ndata, site_name, atemp1, atemp2)
  implicit none

  character(*), intent(in) :: site_name
  integer(8), intent(in)   :: ndata
  real(8), intent(inout)   :: atemp1(ndata), atemp2(ndata)

  integer(8) :: i
  real(8), parameter :: median_delta_atemp = 0.30d0
  real(8), parameter :: SD_atemp = 1.698d0

  do i = 1, ndata
   if (site_name == 'SIGMA-A') then
      if (i == 42976 .or. i == 42977 .or. i == 42978 .or. i == 43124) atemp1(i) = manual_mask
      if (i == 4901) atemp2(i) = manual_mask
      if (i >= 42976 .and. i <= 43146) atemp2(i) = manual_mask
       
      if (atemp1(i) > mask_threshold .and. atemp2(i) > mask_threshold) then
        if (abs(atemp1(i) - atemp2(i)) > (median_delta_atemp + SD_atemp * 3.0d0)) atemp1(i) = mask
      end if
    
      if (i >= 45313) then ! 1:00 on Sep. 1, 2017 -
        if (atemp1(i) > mask_threshold .and. atemp2(i) > mask_threshold) then
          if (abs(atemp1(i) - atemp2(i)) > (median_delta_atemp + SD_atemp * 1.0d0)) atemp1(i) = mask
        end if
      end if
    end if
    
    if (site_name == 'SIGMA-B') then
      if (i == 8391) atemp1(i) = manual_mask
    end if
  end do

end subroutine QC_air_temperature

!----------------------------------------------------------------------------------------------------
subroutine QC_air_humidity (ndata, site_name, arh1, arh2)
  use subroutines_L12_13
  implicit none

  character(*), intent(in) :: site_name
  integer(8), intent(in)   :: ndata
  real(8), intent(inout)   :: arh1(ndata), arh2(ndata)

  integer(8) :: i
  real(8)    :: delta_arh(ndata)
  
  integer(8), parameter :: iflag_refer_former_parameter = 1
  real(8), parameter    :: median_delta_arh = -1.0d0
  real(8), parameter    :: SD_arh = 1.897d0

  do i = 1, ndata    
    if (site_name == 'SIGMA-A') then
      if (arh1(i) > mask_threshold .and. arh2(i) > mask_threshold) then
        if (abs(arh1(i) - arh2(i)) > (median_delta_arh + SD_arh * 3.0d0)) arh1(i) = mask
      end if
    end if

    if (site_name == 'SIGMA-B') then
      if (arh1(i) == 0.0d0) arh1(i) = mask
    end if
  end do
     
  if (site_name == 'SIGMA-B') then
    call refer_formar_parameter(iflag_refer_former_parameter, arh1, ndata, delta_arh)
    do i = 1, ndata
      if (i == 70882 .or. i == 70887 .or. i == 70906 .or. i == 71095) arh1(i) = manual_mask
    end do
  end if
  
end subroutine QC_air_humidity

!----------------------------------------------------------------------------------------------------
subroutine QC_shortwave_radiation (ndata, site_name, SW_TOA, solz, solz_slope, sola, swd, swd_slope, swu)
  use subroutines_L12_13
  implicit none

  character(*), intent(in) :: site_name
  integer(8), intent(in)   :: ndata
  real(8), intent(in)      :: SW_TOA(ndata), solz(ndata), sola(ndata)
  real(8), intent(inout)   :: swd(ndata), swd_slope(ndata), swu(ndata), solz_slope(ndata)

  integer(8) :: i
  real(8)    :: Id_slope(ndata), Is_slope(ndata), diffu_ratio(ndata)

  real(8), parameter :: median_swd_A = 0.0d0
  real(8), parameter :: median_swu_A = 0.40d0
  real(8), parameter :: SD_swd_A = 2.719d0
  real(8), parameter :: SD_swu_A = 2.169d0
  real(8), parameter :: median_swd_B = 0.0d0
  real(8), parameter :: median_swu_B = 0.36d0
  real(8), parameter :: SD_swd_B = 2.521d0
  real(8), parameter :: SD_swu_B = 2.273d0

  do i = 1, ndata
    if (site_name == 'SIGMA-A') then
      if (solz(i) > 90.0d0) then
        if (swd(i) > mask_threshold .and. swd(i) < median_swd_A + SD_swd_A * 3.0d0) swd(i) = 0.0d0
        if (swu(i) > mask_threshold .and. swu(i) < median_swu_A + SD_swu_A * 3.0d0) swu(i) = 0.0d0
      else
        if (SW_TOA(i) /= 0.0d0) then
          if (swd(i) > (co_transA * SW_TOA(i))) swd(i) = mask
        end if
      end if

      if (swd(i) > 0.0d0 .and. swu(i) > 0.0d0) then
        if (swd(i) < swu(i)) swd(i) = radiation_sensor_mask
      end if    

      if (i == 26752 .or. i == 34673 .or. i == 42975 .or. i == 54011 .or. i == 54012 .or. i == 60761 .or. i == 61890) swd(i) = manual_mask
      if (i == 69327 .or. i == 69328) swd(i) = manual_mask
      
      if (i == 25721 .or. i == 26107 .or. i == 26108 .or. i == 26200 .or. i == 27113 .or. i == 34912 .or. i == 34935) swu(i) = manual_mask
      if (i == 35560 .or. i == 41876 .or. i == 41877 .or. i == 43242 .or. i == 51426 .or. i == 51427 .or. i == 51428) swu(i) = manual_mask
      if (i == 52721 .or. i == 60714 .or. i == 60761 .or. i == 60904 .or. i == 60905 .or. i == 60914 .or. i == 60915) swu(i) = manual_mask
      if (i == 61028 .or. i == 61172 .or. i == 62656 .or. i == 68151 .or. i == 68152 .or. i == 68153 .or. i == 68631) swu(i) = manual_mask
      if (i == 69930) swu(i) = manual_mask
      
      if (i >= 70441) then
        swd(i) = mask
        swu(i) = mask
      end if
 
      swd_slope(i) = mask
      solz_slope(i) = mask 
    end if

    if (site_name == 'SIGMA-B') then
      if (solz_slope(i) > 90.0d0) then
        if (swd(i) > mask_threshold .and. swd(i) < median_swd_B + SD_swd_B * 3.0d0) swd(i) = 0.0d0
        if (swu(i) > mask_threshold .and. swu(i) < median_swu_B + SD_swu_B * 3.0d0) swu(i) = 0.0d0
      else
        if (SW_TOA(i) /= 0.0d0) then
          if (swd(i) > (co_transB * SW_TOA(i))) swd(i) = mask
        end if
      end if

      if (i == 43335 .or. i == 43336) swu(i) = manual_mask

      call slope_correction &
        &  (swd(i), solz(i), sola(i), SW_TOA(i), swd_slope(i), Id_slope(i), Is_slope(i), solz_slope(i), diffu_ratio(i))

      if (swd(i) > 0.0d0 .and. swu(i) > 0.0d0) then
        if (swd(i) < swu(i)) swd(i) = radiation_sensor_mask
      end if
      if (swd_slope(i) > 0.0d0 .and. swu(i) > 0.0d0) then
        if (swd_slope(i) < swu(i)) swd_slope(i) = radiation_sensor_mask
      end if
    end if
  end do

end subroutine QC_shortwave_radiation

!----------------------------------------------------------------------------------------------------
subroutine QC_nir_infrared_radiation (ndata, site_name, SW_TOA, solz, solz_slope, nird, niru)
  implicit none

  character(*), intent(in) :: site_name
  integer(8), intent(in)   :: ndata
  real(8), intent(in)      :: SW_TOA(ndata), solz(ndata), solz_slope(ndata)
  real(8), intent(inout)   :: nird(ndata), niru(ndata)

  integer(8) :: i
  real(8)    :: nir_TOA(ndata)

  real(8), parameter :: assumed_max_albedo = 0.90d0
  real(8), parameter :: median_nird_A = 0.0d0
  real(8), parameter :: SD_nird_A = 1.733d0
  real(8), parameter :: median_niru_A = 0.0d0
  real(8), parameter :: SD_niru_A = 1.578d0

  do i = 1, ndata
    nir_TOA(i) = SW_TOA(i) * ratio_nir
    if (solz(i) > 90.0d0) then
      if (nird(i) > mask_threshold .and. nird(i) < median_nird_A + SD_nird_A * 3.0d0) nird(i) = 0.0d0
      if (niru(i) > mask_threshold .and. niru(i) < median_niru_A + SD_niru_A * 3.0d0) niru(i) = 0.0d0
    else
      if (nir_TOA(i) /= 0.0d0) then
        if (nird(i) > (co_transA * nir_TOA(i))) nird(i) = mask
      end if
    end if

    if (nird(i) > 0.0d0 .and. niru(i) > 0.0d0) then
      if (nird(i) < niru(i)) nird(i) = radiation_sensor_mask
    end if

    if (site_name == 'SIGMA-A') then
      if (i == 26200 .or. i == 26752 .or. i == 45901 .or. i == 60761 .or. i == 60809 .or. i == 60904) nird(i) = manual_mask
      if (i == 61809) nird(i) = manual_mask
    
      if (i == 5873 .or. i == 5874 .or. i == 5875 .or. i == 5897 .or. i == 5898 .or. i == 5899) niru(i) = manual_mask
      if (i == 5900 .or. i == 5901 .or. i == 5921 .or. i == 5922 .or. i == 5944 .or. i == 5945) niru(i) = manual_mask
      if (i == 5946 .or. i == 5947 .or. i == 5948 .or. i == 6015 .or. i == 6016 .or. i == 6017) niru(i) = manual_mask
      if (i == 6018 .or. i == 6019 .or. i == 6020 .or. i == 6064 .or. i == 6065 .or. i == 6066) niru(i) = manual_mask
      if (i == 9344 .or. i == 9368 .or. i == 9375 .or. i == 9377 .or. i == 9379 .or. i == 9382) niru(i) = manual_mask
      
      if (i == 10603 .or. i == 11153 .or. i == 11154 .or. i == 11155 .or. i == 11156 .or. i == 18089) niru(i) = manual_mask
     
      if (i == 25721 .or. i == 26107 .or. i == 26108 .or. i == 26200 .or. i == 27113) niru(i) = manual_mask
     
      if (i == 31936 .or. i == 31937 .or. i == 31938 .or. i == 32175 .or. i == 32176 .or. i == 32177) niru(i) = manual_mask
      if (i == 32178 .or. i == 32179 .or. i == 33929 .or. i == 34935) niru(i) = manual_mask
     
      if (i == 40867 .or. i == 41876 .or. i == 41877 .or. i == 42953 .or. i == 43242) niru(i) = manual_mask
     
      if (i == 51426 .or. i == 51427 .or. i == 51428 .or. i == 52721 .or. i == 54973 .or. i == 54974) niru(i) = manual_mask
     
      if (i == 60714 .or. i == 60761 .or. i == 60904 .or. i == 60905 .or. i == 61028 .or. i == 61029) niru(i) = manual_mask
      if (i == 61172 .or. i == 61890 .or. i == 62656 .or. i == 67072 .or. i == 67095 .or. i == 68151) niru(i) = manual_mask
      if (i == 68152 .or. i == 68153 .or. i == 68631) niru(i) = manual_mask

      if (i >= 70441) then
        nird(i) = mask
        niru(i) = mask
      end if
    end if 

    if (site_name == 'SIGMA-B') then
    end if
  end do

end subroutine QC_nir_infrared_radiation

!----------------------------------------------------------------------------------------------------
subroutine QC_albedo (ndata, site_name, month, swd, swd_slope, swu, SW_TOA, solz, solz_slope, lwd, albedo)
  use subroutines_L12_13
  implicit none

  character(*), intent(in) :: site_name
  integer(8), intent(in)   :: ndata
  integer(8), intent(in)   :: month(ndata)
  real(8), intent(in)      :: swd(ndata), swd_slope(ndata), swu(ndata), SW_TOA(ndata)
  real(8), intent(in)      :: solz(ndata), solz_slope(ndata), lwd(ndata)
  real(8), intent(inout)   :: albedo(ndata)

  integer(8) :: i
  real(8), parameter :: median_swd_solz85_A = 28.6d0
  real(8), parameter :: median_swd_solz80_A = 61.7d0
  real(8), parameter :: median_swd_solz85_B = 24.4d0
  real(8), parameter :: median_swd_solz80_B = 55.6d0
  real(8), parameter :: median_tau_A = 0.483d0
  real(8), parameter :: SD_tau_A = 0.184d0
  real(8), parameter :: median_tau_B = 0.344d0
  real(8), parameter :: SD_tau_B = 0.458d0
  
  do i = 1, ndata
    if (site_name == 'SIGMA-A') then
      if (swd(i) <= 0.0d0 .or. swu(i) <= 0.0d0) then
        albedo(i) = mask
      else
        albedo(i) = swu(i) / swd(i)
      end if

      if (swd(i) < median_swd_solz85_A) albedo(i) = mask
      if (solz(i) > 85.0d0) albedo(i) = mask

      ! in warm period (from May to September): surface condition is assumed to be dry or wet snow
      if (month(i) == 5 .or. month(i) == 6 .or. month(i) == 7 .or. month(i) == 8 .or. month(i) == 9) then
        if (albedo(i) < LL_WS_albedo .or. albedo(i) > UL_S_albedo) albedo(i) = mask
        if (solz(i) > 80.0d0) then
          if (swd(i) > median_swd_solz80_A .and. albedo(i) < 0.60d0) albedo(i) = mask
          if (swd(i) / SW_TOA(i) > (median_tau_A + SD_tau_A * 1.0d0) .and. albedo(i) < 0.60d0) albedo(i) = mask
        end if
 
      ! in cold period (from October to April): surface condition is assumed to be dry snow
      else
        if (albedo(i) < 0.60d0 .or. albedo(i) > UL_S_albedo) albedo(i) = mask
         if (solz(i) > 80.0d0) then
           if (swd(i) > median_swd_solz80_A .and. albedo(i) < 0.70d0) albedo(i) = mask
           if (swd(i) / SW_TOA(i) > (median_tau_A + SD_tau_A * 1.0d0) .and. albedo(i) < 0.70d0) albedo(i) = mask
        end if
      end if

      if (lwd(i) == radiation_sensor_mask) albedo(i) = mask
      if (i == 1891 .or. i == 25208 .or. i == 27876) albedo(i) = manual_mask
    end if

    if (site_name == 'SIGMA-B') then
      if (swd_slope(i) <= 0.0d0 .or. swu(i) <= 0.0d0) then
        albedo(i) = mask
      else
        albedo(i) = swu(i) / swd_slope(i)
      end if
				
      if (swd_slope(i) < median_swd_solz85_B) albedo(i) = mask
      if (solz_slope(i) > 85.0d0) albedo(i) = mask

      ! in warm period (from May to September): surface condition is assumed to be dry/wet snow or bare/dark ice
      if (month(i) == 5 .or. month(i) == 6 .or. month(i) == 7 .or. month(i) == 8 .or. month(i) == 9) then  
        if (albedo(i) < LL_DI_albedo .or. albedo(i) > UL_S_albedo) albedo(i) = mask
        if (solz_slope(i) > 80.0d0) then
          if (swd_slope(i) > median_swd_solz80_B .and. albedo(i) < 0.30d0) albedo(i) = mask
          if (swd_slope(i) / SW_TOA(i) > (median_tau_B + SD_tau_B * 1.0d0) .and. albedo(i) < 0.30d0) albedo(i) = mask
        end if
        
      ! in cold period (from October to April): surface condition is assumed to be dry snow        
      else
        if (albedo(i) < LL_WS_albedo .or. albedo(i) > UL_S_albedo) albedo(i) = mask
        if (solz_slope(i) > 80.0d0) then
          if (swd_slope(i) > median_swd_solz80_B .and. albedo(i) < 0.60d0) albedo(i) = mask
          if (swd_slope(i) / SW_TOA(i) > (median_tau_B + SD_tau_B * 1.0d0) .and. albedo(i) < 0.60d0) albedo(i) = mask
        end if
      end if

      if (lwd(i) == radiation_sensor_mask) albedo(i) = mask  
      if (i == 1579 .or. i == 58146 .or. i == 58218) albedo(i) = manual_mask
    end if
  end do

end subroutine QC_albedo

!----------------------------------------------------------------------------------------------------
subroutine QC_nir_infrared_albedo (ndata, site_name, month, nird, niru, solz, lwd, nir_albedo)
  implicit none

  character(*), intent(in) :: site_name
  integer(8), intent(in)   :: ndata, month(ndata)
  real(8), intent(in)      :: nird(ndata), niru(ndata), solz(ndata), lwd(ndata)
  real(8), intent(inout)   :: nir_albedo(ndata)

  integer(8) :: i
  real(8), parameter :: median_nird_solz85 = 11.2d0
  real(8), parameter :: median_nird_solz80 = 26.3d0

  do i = 1, ndata
    if (site_name == 'SIGMA-A') then
      if (nird(i) <= 0.0d0 .or. niru(i) <= 0.0d0) then
        nir_albedo(i) = mask
      else
        nir_albedo(i) = niru(i) / nird(i)
      end if

      if (nird(i) < median_nird_solz85) nir_albedo(i) = mask
      if (solz(i) > 85.0d0) nir_albedo(i) = mask

      ! in warm period (from May to September): surface condition is assumed to be dry or wet snow
      if (month(i) == 5 .or. month(i) == 6 .or. month(i) == 7 .or. month(i) == 8 .or. month(i) == 9) then
        if (nir_albedo(i) < LL_WS_nir_albedo .or. nir_albedo(i) > UL_S_nir_albedo) nir_albedo(i) = mask
        if (solz(i) > 80.0d0) then
          if (nird(i) > median_nird_solz80 .and. nir_albedo(i) < 0.50d0) nir_albedo(i) = mask
        end if

      ! in cold period (from October to April): surface condition is assumed to be dry snow        
      else
        if (nir_albedo(i) < LL_DS_nir_albedo .or. nir_albedo(i) > UL_S_nir_albedo) nir_albedo(i) = mask
        if (solz(i) > 80.0d0) then
          if (nird(i) > median_nird_solz80 .and. nir_albedo(i) < 0.60d0) nir_albedo(i) = mask
        end if
      end if

      if (lwd(i) == radiation_sensor_mask) nir_albedo(i) = mask
    end if

    if (site_name == 'SIGMA-B') then
      nir_albedo(i) = mask
    end if

  end do

end subroutine QC_nir_infrared_albedo

!----------------------------------------------------------------------------------------------------
subroutine QC_albedo_acc (ndata, site_name, month, albedo_acc, nir_albedo_acc)
  implicit none

  character(*), intent(in) :: site_name
  integer(8), intent(in)   :: ndata
  integer(8), intent(in)   :: month(ndata)
  real(8), intent(inout)   :: albedo_acc(ndata)
  real(8), intent(inout)   :: nir_albedo_acc(ndata)

  integer(8) :: i

  do i = 1, ndata
    ! in warm period (from May to September):
    if (month(i) == 5 .or. month(i) == 6 .or. month(i) == 7 .or. month(i) == 8 .or. month(i) == 9) then
      if (site_name == 'SIGMA-A') then
        if (albedo_acc(i) < LL_WS_albedo .or. albedo_acc(i) > UL_S_albedo) albedo_acc(i) = mask
        if (nir_albedo_acc(i) < LL_WS_nir_albedo .or. nir_albedo_acc(i) > UL_S_nir_albedo) nir_albedo_acc(i) = mask
      end if
      if (site_name == 'SIGMA-B') then
        if (albedo_acc(i) < LL_DI_albedo .or. albedo_acc(i) > UL_S_albedo) albedo_acc(i) = mask
      end if
      
    ! in cold period (from October to April):
    else
      if (site_name == 'SIGMA-A') then
        if (albedo_acc(i) < LL_DS_albedo .or. albedo_acc(i) > UL_S_albedo) albedo_acc(i) = mask
        if (nir_albedo_acc(i) < LL_DS_nir_albedo .or. nir_albedo_acc(i) > UL_S_nir_albedo) nir_albedo_acc(i) = mask
      end if
      if (site_name == 'SIGMA-B') then
        if (albedo_acc(i) < LL_WS_albedo .or. albedo_acc(i) > UL_S_albedo) albedo_acc(i) = mask
        nir_albedo_acc(i) = mask
      end if
    end if
  end do

end subroutine QC_albedo_acc

!----------------------------------------------------------------------------------------------------
subroutine nir_infrared_fraction (ndata, site_name, swd, nird, solz, nir_frac)
  implicit none

  character(*), intent(in) :: site_name
  integer(8), intent(in)   :: ndata
  real(8), intent(in)      :: swd(ndata)
  real(8), intent(in)      :: nird(ndata)
  real(8), intent(in)      :: solz(ndata)
  real(8), intent(inout)   :: nir_frac(ndata)

  integer(8) :: i

  do i = 1, ndata
    if (site_name == 'SIGMA-A') then
      if (solz(i) <= 90.0d0) then
        if (nird(i) > 0.0d0 .or. swd(i) > 0.0d0) then
          nir_frac(i) = nird(i) / swd(i)
          if (nir_frac(i) < 0.10d0) nir_frac(i) = mask
          if (nir_frac(i) > ratio_nir) nir_frac(i) = mask
        else
          nir_frac(i) = mask
        end if
      else
        nir_frac(i) = mask
      end if
    end if

    if (site_name == 'SIGMA-B') then
      nir_frac(i) = mask
    end if  
  end do

end subroutine nir_infrared_fraction

!----------------------------------------------------------------------------------------------------
subroutine QC_longwave_radiation (ndata, site_name, atemp1, atemp2, lwd, lwu)
  use subroutines_L12_13
  implicit none

  character(*), intent(in) :: site_name
  integer(8), intent(in)   :: ndata
  real(8), intent(in)      :: atemp1(ndata), atemp2(ndata)
  real(8), intent(inout)   :: lwd(ndata), lwu(ndata)

  integer(8) :: i
  real(8) :: atemp(ndata)
  real(8) :: lw_standard

  real(8), parameter :: arh_threshold = 92.61d0
  real(8), parameter :: median_delta_lwd = 11.76d0
  real(8), parameter :: median_delta_lwu = 47.89d0
  real(8), parameter :: SD_lwd = 41.40d0
  real(8), parameter :: SD_lwu = 21.22d0

  if (site_name == 'SIGMA-A') then
    atemp = atemp2
  end if
  if (site_name == 'SIGMA-B') then
    atemp = atemp1
  end if

  
  do i = 1, ndata
  
    if (site_name == 'SIGMA-A') then
      if (atemp(i) > mask_threshold) then
        ! Brock and Arnold (2000)
        lw_standard = standard_longwave_radiation (atemp(i), 0.5d0)
        if (lwd(i) > mask_threshold) then

          if (abs(lwd(i) - lw_standard) > (median_delta_lwd + SD_lwd * 2)) lwd(i) = mask 
          if (i >= 41803 .and. i <= 43182) then ! 19:00 on Apr. 7, 2017 - 6:00 on Jun. 4, 2017
            if (abs(lwd(i) - lw_standard) > (median_delta_lwd + SD_lwd)) lwd(i) = mask
          end if
          if (i >= 45313) then ! 1:00 on Sep. 1, 2017 -
            if (abs(lwd(i) - lw_standard) > (median_delta_lwd + SD_lwd)) lwd(i) = mask
          end if
        end if
        
        if (lwu(i) > mask_threshold) then
          if (abs(lwu(i) - lw_standard) > (median_delta_lwu + SD_lwu * 2)) lwu(i) = mask
        end if          
      else
        lw_standard = mask
      end if
      
      if (i == 43064 .or. i == 43065 .or. i == 43138 .or. i == 45895 .or. i == 46370 .or. i == 46707 .or. i == 46708) lwd(i) = manual_mask
      
      if (i == 41803 .or. i == 41967 .or. i == 41970 .or. i == 41996 .or. i == 42021 .or. i == 42028 .or. i == 42032) lwu(i) = manual_mask
      if (i == 42033 .or. i == 42078 .or. i == 42094 .or. i == 42097 .or. i == 42102 .or. i == 42119 .or. i == 42160) lwu(i) = manual_mask
      if (i == 42129 .or. i == 42130 .or. i == 42138 .or. i == 42139 .or. i == 42140 .or. i == 42141 .or. i == 42142) lwu(i) = manual_mask
      if (i == 42158 .or. i == 42161 .or. i == 42165 .or. i == 42196 .or. i == 42197 .or. i == 42207 .or. i == 42214) lwu(i) = manual_mask
      if (i == 42216 .or. i == 42220 .or. i == 42225 .or. i == 42226 .or. i == 42229 .or. i == 42237 .or. i == 42238) lwu(i) = manual_mask
      if (i == 42239 .or. i == 42240 .or. i == 42241 .or. i == 42245 .or. i == 42247 .or. i == 42253 .or. i == 42255) lwu(i) = manual_mask
      if (i == 42256 .or. i == 42258 .or. i == 42259 .or. i == 42260 .or. i == 42279 .or. i == 42280 .or. i == 42281) lwu(i) = manual_mask
      if (i == 42282 .or. i == 42285 .or. i == 42286 .or. i == 42287 .or. i == 42289 .or. i == 42290 .or. i == 42291) lwu(i) = manual_mask
      if (i == 42292 .or. i == 42293 .or. i == 42294 .or. i == 42295 .or. i == 42296 .or. i == 42297 .or. i == 42298) lwu(i) = manual_mask
      if (i == 42299 .or. i == 42301 .or. i == 42302 .or. i == 42305 .or. i == 42306 .or. i == 42309 .or. i == 42313) lwu(i) = manual_mask
      if (i == 42314 .or. i == 42315 .or. i == 42316 .or. i == 42319 .or. i == 42321 .or. i == 42324 .or. i == 42326) lwu(i) = manual_mask
      if (i == 42328 .or. i == 42331 .or. i == 42334 .or. i == 42336 .or. i == 42337 .or. i == 42341 .or. i == 42345) lwu(i) = manual_mask
      if (i == 42347 .or. i == 42348 .or. i == 42350 .or. i == 42360 .or. i == 42361 .or. i == 42362 .or. i == 42364) lwu(i) = manual_mask
      if (i == 42365 .or. i == 42366 .or. i == 42367 .or. i == 42368 .or. i == 42369 .or. i == 42371 .or. i == 42376) lwu(i) = manual_mask
      if (i == 42377 .or. i == 42387 .or. i == 42397 .or. i == 42398 .or. i == 42409 .or. i == 42424 .or. i == 42425) lwu(i) = manual_mask
      if (i == 42426 .or. i == 42468 .or. i == 42477 .or. i == 42478 .or. i == 42479 .or. i == 42522 .or. i == 42523) lwu(i) = manual_mask
      if (i == 42524 .or. i == 42532 .or. i == 42544 .or. i == 42550 .or. i == 43009 .or. i == 43013 .or. i == 43064) lwu(i) = manual_mask
      if (i == 43065 .or. i == 43106 .or. i == 43138 .or. i == 47349 .or. i == 47363 .or. i == 47364 .or. i == 47365) lwu(i) = manual_mask
      if (i == 48626 .or. i == 48789 .or. i == 49561 .or. i == 49562 .or. i == 49563 .or. i == 49565 .or. i == 51391) lwu(i) = manual_mask
      if (i == 59034 .or. i == 59038 .or. i == 59039 .or. i == 59290 .or. i == 59291 .or. i == 59292 .or. i == 59463) lwu(i) = manual_mask
      if (i == 59471 .or. i == 63195 .or. i == 63312 .or. i == 63347) lwu(i) = manual_mask
    end if
    
    if (site_name == 'SIGMA-B') then
      if (atemp(i) > mask_threshold) then
        ! Brock and Arnold (2000)
        lw_standard = standard_longwave_radiation (atemp(i), 0.5d0)
      else
        lw_standard = mask
      end if
    end if  
  end do

end subroutine QC_longwave_radiation

!----------------------------------------------------------------------------------------------------
subroutine QC_snow_depth (ndata, site_name, sd)
  use subroutines_L12_13
  implicit none

  character(*), intent(in) :: site_name
  integer(8), intent(in)   :: ndata
  real(8), intent(inout)   :: sd(ndata)

  integer(8) :: i
  integer(8) :: referring_period
  real(8)    :: regression_sd
  real(8)    :: sd_mean(ndata)

  real(8), parameter :: median_delta_sdA1 = 0.15d0
  real(8), parameter :: median_delta_sdA2 = 0.07d0
  real(8), parameter :: SD_sdA1 = 16.02d0
  real(8), parameter :: SD_sdA2 = 2.41d0
  real(8), parameter :: median_delta_sdB1 = -0.01d0
  real(8), parameter :: median_delta_sdB2 = -0.12d0
  real(8), parameter :: SD_sdB1 = 19.55d0
  real(8), parameter :: SD_sdB2 = 5.51d0
  
  real(8), parameter :: slope_sd = 0.006359d0
  real(8), parameter :: intercept_sd = 166.2894d0
  real(8), parameter :: SD_sdA_reg = 26.55d0
  
  if (site_name == 'SIGMA-A') then
    do i = 1, ndata
      sd(19) = manual_mask
      if (i >= 27769 .and. i <= 35808) then ! 1:00 on Sep. 1, 2015 - 0:00 on Aug. 1, 2016 
        if (sd(i) > 280.0d0) sd(i) = manual_mask
      end if
    end do

    call before_average(1, sd, ndata, sd_mean, referring_period)
    do i = referring_period, ndata
      if (sd(i) > mask_threshold.and. sd_mean(i) > mask_threshold) then
        if (abs(sd(i) - sd_mean(i)) > median_delta_sdA1 + SD_sdA1) sd(i) = mask
      end if
    end do

    do i = 1, ndata
      if (i >= 45313) then ! 1:00 on Sep. 1, 2017 -
        regression_sd = i * slope_sd + intercept_sd
        if (abs(sd(i) - regression_sd) > SD_sdA_reg) sd(i) = mask
      end if
    end do

    call before_average(4, sd, ndata, sd_mean, referring_period)
    do i = referring_period, ndata
      if (sd(i) > mask_threshold.and. sd_mean(i) > mask_threshold) then
        if (abs(sd(i) - sd_mean(i)) > median_delta_sdA2 + SD_sdA2 * 3) sd(i) = mask
      end if
    end do

    do i = 1, ndata
      if (i == 14249 .or. i == 14250 .or. i == 14299 .or. i == 14341 .or. i == 14342 .or. i == 14343 .or. i == 14353) sd(i) = manual_mask
      if (i == 14524 .or. i == 14538 .or. i == 14550 .or. i == 14552 .or. i == 14556 .or. i == 14567 .or. i == 14586) sd(i) = manual_mask
      if (i == 14587 .or. i == 14588 .or. i == 14857 .or. i == 14858 .or. i == 14859 .or. i == 15612 .or. i == 18741) sd(i) = manual_mask
      if (i == 18742 .or. i == 19351 .or. i == 19859 .or. i == 19944) sd(i) = manual_mask
      
      if (i == 20093 .or. i == 20094 .or. i == 20095 .or. i == 20096 .or. i == 20097 .or. i == 20649 .or. i == 20650) sd(i) = manual_mask
      if (i == 20651 .or. i == 20653 .or. i == 20654 .or. i == 20730 .or. i == 20758 .or. i == 20761 .or. i == 20813) sd(i) = manual_mask
      if (i == 20814 .or. i == 20815 .or. i == 20819 .or. i == 20821 .or. i == 20898 .or. i == 20982 .or. i == 22825) sd(i) = manual_mask
      if (i == 23255 .or. i == 23548 .or. i == 23849 .or. i == 23889 .or. i == 23913 .or. i == 23916 .or. i == 23919) sd(i) = manual_mask
      if (i == 24063 .or. i == 24069 .or. i == 24080 .or. i == 24102 .or. i == 24103 .or. i == 24434 .or. i == 24435) sd(i) = manual_mask
      if (i == 24438 .or. i == 24443 .or. i == 24656 .or. i == 24657 .or. i == 24658 .or. i == 25156 .or. i == 25157) sd(i) = manual_mask
      if (i == 25252 .or. i == 25255 .or. i == 25264 .or. i == 25299 .or. i == 25300 .or. i == 28336 .or. i == 28340) sd(i) = manual_mask
      if (i == 28355 .or. i == 28591 .or. i == 28625 .or. i == 28762 .or. i == 29114 .or. i == 29238 .or. i == 29507) sd(i) = manual_mask
      if (i == 29654 .or. i == 29863 .or. i == 29865 .or. i == 29869 .or. i == 29870) sd(i) = manual_mask

      if (i == 30452 .or. i == 30676 .or. i == 31059 .or. i == 31062 .or. i == 31063 .or. i == 31064 .or. i == 31065) sd(i) = manual_mask
      if (i == 31067 .or. i == 31091 .or. i == 31301 .or. i == 32117 .or. i == 34439 .or. i == 34441 .or. i == 34444) sd(i) = manual_mask
      if (i == 34445 .or. i == 34487 .or. i == 34488 .or. i == 34493 .or. i == 34494 .or. i == 34496 .or. i == 34497) sd(i) = manual_mask
      if (i == 34518 .or. i == 34536 .or. i == 34543 .or. i == 34619 .or. i == 34910 .or. i == 34914 .or. i == 34921) sd(i) = manual_mask
      if (i == 34927 .or. i == 35033 .or. i == 35043 .or. i == 35069 .or. i == 35114 .or. i == 35141 .or. i == 35143) sd(i) = manual_mask
      if (i == 35151 .or. i == 35159 .or. i == 35160 .or. i == 35164 .or. i == 35168 .or. i == 35174 .or. i == 35183) sd(i) = manual_mask
      if (i == 35185 .or. i == 35236 .or. i == 35329 .or. i == 35391 .or. i == 35394 .or. i == 35512 .or. i == 35543) sd(i) = manual_mask
      if (i == 35547 .or. i == 35551 .or. i == 35553 .or. i == 37825) sd(i) = manual_mask
      
      if (i == 43718 .or. i == 46752 .or. i == 47757 .or. i == 47778 .or. i == 47782 .or. i == 47806 .or. i == 47927) sd(i) = manual_mask
      if (i == 48006 .or. i == 48088 .or. i == 48137 .or. i == 48218 .or. i == 48234 .or. i == 48235 .or. i == 48372) sd(i) = manual_mask
      if (i == 48383 .or. i == 48705 .or. i == 48858 .or. i == 48931 .or. i == 48941 .or. i == 49252 .or. i == 49308) sd(i) = manual_mask
      if (i == 49351 .or. i == 49405 .or. i == 49448 .or. i == 49452 .or. i == 49502) sd(i) = manual_mask
      
      if (i == 50495 .or. i == 51022 .or. i == 51288 .or. i == 51643 .or. i == 51645 .or. i == 51667 .or. i == 51512) sd(i) = manual_mask
      if (i == 51613 .or. i == 51909 .or. i == 52013 .or. i == 52021 .or. i == 52022 .or. i == 52023 .or. i == 52024) sd(i) = manual_mask
      if (i == 52025 .or. i == 52026 .or. i == 52027 .or. i == 52909 .or. i == 52964 .or. i == 52965 .or. i == 53084) sd(i) = manual_mask
      if (i == 53090 .or. i == 53101 .or. i == 53103 .or. i == 53109 .or. i == 53133 .or. i == 53135 .or. i == 53136) sd(i) = manual_mask
      if (i == 53138 .or. i == 53139 .or. i == 53140 .or. i == 53209 .or. i == 53238 .or. i == 53239 .or. i == 53242) sd(i) = manual_mask
      if (i == 53256 .or. i == 53259 .or. i == 53264 .or. i == 53265 .or. i == 53266 .or. i == 53269 .or. i == 53273) sd(i) = manual_mask
      if (i == 53274 .or. i == 53276 .or. i == 53278 .or. i == 53283 .or. i == 53284 .or. i == 53288 .or. i == 53289) sd(i) = manual_mask
      if (i == 53292 .or. i == 53306 .or. i == 53308 .or. i == 53441 .or. i == 53449 .or. i == 53459 .or. i == 53503) sd(i) = manual_mask
      if (i == 53524 .or. i == 53617 .or. i == 53955 .or. i == 54524 .or. i == 55165 .or. i == 55428 .or. i == 55588) sd(i) = manual_mask
      if (i == 55624 .or. i == 55659 .or. i == 55803 .or. i == 55871 .or. i == 55994 .or. i == 56035 .or. i == 56963) sd(i) = manual_mask
      if (i == 57034 .or. i == 57510 .or. i == 57743 .or. i == 58183 .or. i == 58196 .or. i == 58287 .or. i == 58358) sd(i) = manual_mask
      if (i == 58375 .or. i == 58381 .or. i == 58425 .or. i == 58474 .or. i == 58478 .or. i == 58675 .or. i == 58721) sd(i) = manual_mask
      if (i == 58722 .or. i == 58725 .or. i == 58735 .or. i == 58746 .or. i == 58820 .or. i == 59359 .or. i == 59667) sd(i) = manual_mask
      if (i == 59749 .or. i == 59752 .or. i == 59759 .or. i == 59760) sd(i) = manual_mask
      
      if (i == 62897 .or. i == 63496 .or. i == 63515 .or. i == 63573 .or. i == 63583 .or. i == 63691 .or. i == 63746) sd(i) = manual_mask
      if (i == 63752 .or. i == 63961 .or. i == 64001 .or. i == 64194 .or. i == 64277 .or. i == 64606 .or. i == 64638) sd(i) = manual_mask
      if (i == 64645 .or. i == 64654 .or. i == 64674 .or. i == 64706 .or. i == 64781 .or. i == 64800 .or. i == 64852) sd(i) = manual_mask
      if (i == 64853 .or. i == 64858 .or. i == 64873 .or. i == 64894 .or. i == 65005 .or. i == 65106 .or. i == 65162) sd(i) = manual_mask
      if (i == 65177 .or. i == 65288 .or. i == 65369 .or. i == 65511 .or. i == 65518 .or. i == 65525 .or. i == 65548) sd(i) = manual_mask
      if (i == 65560 .or. i == 65569 .or. i == 65699 .or. i == 65721 .or. i == 67752 .or. i == 67898 .or. i == 68257) sd(i) = manual_mask
      if (i == 68498 .or. i == 68642 .or. i == 68820 .or. i == 68898 .or. i == 68487 .or. i == 68798 .or. i == 68848) sd(i) = manual_mask
      if (i == 69221 .or. i == 69282 .or. i == 69363 .or. i == 69375 .or. i == 69376 .or. i == 69567 .or. i == 69588) sd(i) = manual_mask
      if (i == 69601 .or. i == 69677 .or. i == 69920 .or. i == 69935 .or. i == 69999) sd(i) = manual_mask
      
      if (i == 70055 .or. i == 70058 .or. i == 70063 .or. i == 70225 .or. i == 70242 .or. i == 70250 .or. i == 70252) sd(i) = manual_mask
      if (i == 70253 .or. i == 70259 .or. i == 70263 .or. i == 70265 .or. i == 70266 .or. i == 70268 .or. i == 70285) sd(i) = manual_mask
      if (i == 70313 .or. i == 70399 .or. i == 70405 .or. i == 70435) sd(i) = manual_mask
    end do
  end if

  if (site_name == 'SIGMA-B') then
    call before_average(1, sd, ndata, sd_mean, referring_period)
    do i = referring_period, ndata
      if (sd(i) > mask_threshold .and. sd_mean(i) > mask_threshold) then
        if (abs(sd(i) - sd_mean(i)) > median_delta_sdB1 + SD_sdB1 * 2) sd(i) = mask
      end if
    end do

    call before_average(4, sd, ndata, sd_mean, referring_period)
    do i = referring_period, ndata
      if (sd(i) > mask_threshold .and. sd_mean(i) > mask_threshold) then
        if (abs(sd(i) - sd_mean(i)) > median_delta_sdB2 + SD_sdB2 * 3) sd(i) = mask
      end if
    end do

    do i = 1, ndata
      if (i == 1138 .or. i == 1139 .or. i == 1140 .or. i == 1141 .or. i == 1142 .or. i == 1143 .or. i == 1144) sd(i) = manual_mask
      if (i == 1145 .or. i == 1146 .or. i == 1147 .or. i == 1148 .or. i == 1149 .or. i == 1150 .or. i == 1151) sd(i) = manual_mask
      if (i == 1153 .or. i == 1239 .or. i == 1240 .or. i == 2410 .or. i == 2411 .or. i == 2432 .or. i == 2433) sd(i) = manual_mask
      if (i == 2885 .or. i == 2886 .or. i == 4196 .or. i == 4204 .or. i == 4205 .or. i == 4206 .or. i == 4211) sd(i) = manual_mask
      if (i == 4216 .or. i == 4217 .or. i == 4218 .or. i == 4527 .or. i == 4551 .or. i == 4730 .or. i == 4802) sd(i) = manual_mask
      if (i == 4808 .or. i == 4811 .or. i == 4820 .or. i == 4825 .or. i == 4909 .or. i == 5473 .or. i == 5474) sd(i) = manual_mask
      if (i == 5480 .or. i == 5733 .or. i == 6102 .or. i == 6247 .or. i == 6757 .or. i == 6758 .or. i == 6849) sd(i) = manual_mask
      if (i == 6850 .or. i == 6869 .or. i == 7683 .or. i == 7698 .or. i == 7699 .or. i == 8000 .or. i == 8001) sd(i) = manual_mask
      if (i == 8003 .or. i == 9109 .or. i == 9110 .or. i == 9418 .or. i == 9445 .or. i == 9446 .or. i == 9487) sd(i) = manual_mask
      if (i == 9488 .or. i == 9919 .or. i == 9920 .or. i == 9921 .or. i == 9922 .or. i == 9923 .or. i == 9924) sd(i) = manual_mask
      if (i == 9925 .or. i == 9926 .or. i == 9927 .or. i == 9928 .or. i == 9929 .or. i == 9930 .or. i == 9931) sd(i) = manual_mask
      if (i == 9932 .or. i == 9933) sd(i) = manual_mask
      
      if (i == 10034 .or. i == 10710 .or. i == 10711 .or. i == 11864 .or. i == 11865 .or. i == 11866 .or. i == 11867) sd(i) = manual_mask
      if (i == 11868 .or. i == 11869 .or. i == 11956 .or. i == 11961 .or. i == 11962 .or. i == 11963 .or. i == 12142) sd(i) = manual_mask
      if (i == 12143 .or. i == 12144 .or. i == 12943 .or. i == 14868 .or. i == 15114 .or. i == 16173 .or. i == 16267) sd(i) = manual_mask
      if (i == 16938 .or. i == 16939 .or. i == 16945 .or. i == 16965 .or. i == 16966 .or. i == 16967 .or. i == 16968) sd(i) = manual_mask
      if (i == 16969 .or. i == 16970 .or. i == 16971 .or. i == 16972 .or. i == 16973 .or. i == 16974 .or. i == 16975) sd(i) = manual_mask
      if (i == 16976 .or. i == 16977 .or. i == 16978 .or. i == 16979 .or. i == 16980 .or. i == 16981 .or. i == 16982) sd(i) = manual_mask
      if (i == 17005 .or. i == 17013 .or. i == 17014 .or. i == 17038 .or. i == 17076 .or. i == 18741 .or. i == 18867) sd(i) = manual_mask
      if (i == 19106 .or. i == 19107 .or. i == 19108 .or. i == 19109 .or. i == 19110) sd(i) = manual_mask
      
      if (i == 20279 .or. i == 20280 .or. i == 20414 .or. i == 23457 .or. i == 23461 .or. i == 23464 .or. i == 24795) sd(i) = manual_mask
      if (i == 27891 .or. i == 27903 .or. i == 27904 .or. i == 27905 .or. i == 27906 .or. i == 27907 .or. i == 27908) sd(i) = manual_mask
      if (i == 27909 .or. i == 27914 .or. i == 27915 .or. i == 28167 .or. i == 28168 .or. i == 28506 .or. i == 28508) sd(i) = manual_mask
      if (i == 28748 .or. i == 29097 .or. i == 29098 .or. i == 29099 .or. i == 29777 .or. i == 29778) sd(i) = manual_mask
      
      if (i == 30073 .or. i == 30075 .or. i == 30904 .or. i == 36747 .or. i == 36754 .or. i == 36756 .or. i == 36758) sd(i) = manual_mask
      if (i == 36760 .or. i == 36761 .or. i == 36762 .or. i == 36825 .or. i == 36826 .or. i == 36834 .or. i == 37751) sd(i) = manual_mask
      if (i == 37753 .or. i == 37756 .or. i == 37759 .or. i == 38383 .or. i == 38387 .or. i == 38390 .or. i == 38391) sd(i) = manual_mask
      if (i == 38392 .or. i == 38422 .or. i == 38427 .or. i == 38428 .or. i == 38432 .or. i == 38688 .or. i == 38706) sd(i) = manual_mask
      if (i == 38707 .or. i == 38715 .or. i == 38716 .or. i == 38722 .or. i == 38723 .or. i == 39067 .or. i == 39068) sd(i) = manual_mask

      if (i == 40498 .or. i == 40500 .or. i == 40503 .or. i == 42166 .or. i == 42974 .or. i == 42984 .or. i == 42985) sd(i) = manual_mask
      if (i == 42986 .or. i == 42987 .or. i == 42988 .or. i == 42989 .or. i == 42990 .or. i == 42991 .or. i == 42992) sd(i) = manual_mask
      if (i == 42993 .or. i == 42998 .or. i == 43001 .or. i == 43002 .or. i == 43003 .or. i == 43006 .or. i == 43007) sd(i) = manual_mask
      if (i == 43008 .or. i == 43009 .or. i == 43010 .or. i == 43011 .or. i == 43012 .or. i == 43016 .or. i == 43018) sd(i) = manual_mask
      if (i == 43019 .or. i == 43031 .or. i == 43033 .or. i == 43034 .or. i == 43035 .or. i == 43036 .or. i == 43037) sd(i) = manual_mask
      if (i == 43038 .or. i == 43039 .or. i == 43040 .or. i == 43041 .or. i == 43042 .or. i == 43043 .or. i == 43044) sd(i) = manual_mask
      if (i == 43045 .or. i == 43047 .or. i == 43067 .or. i == 43068 .or. i == 43069 .or. i == 43071 .or. i == 43075) sd(i) = manual_mask
      if (i == 43076 .or. i == 48536 .or. i == 48537 .or. i == 48653 .or. i == 49617 .or. i == 49986) sd(i) = manual_mask
      
      if (i == 50028 .or. i == 50030 .or. i == 50032 .or. i == 50554 .or. i == 50555 .or. i == 51422 .or. i == 52181) sd(i) = manual_mask
      if (i == 52182 .or. i == 52183 .or. i == 52186 .or. i == 52204 .or. i == 52265 .or. i == 52266 .or. i == 52269) sd(i) = manual_mask
      if (i == 52271 .or. i == 52326 .or. i == 52329 .or. i == 52330 .or. i == 52331 .or. i == 52332 .or. i == 52333) sd(i) = manual_mask
      if (i == 52334 .or. i == 52335 .or. i == 52338 .or. i == 52339 .or. i == 52820 .or. i == 52821 .or. i == 52993) sd(i) = manual_mask
      if (i == 53010 .or. i == 54102 .or. i == 55501 .or. i == 55575 .or. i == 55582 .or. i == 56355 .or. i == 56356) sd(i) = manual_mask
      if (i == 56434 .or. i == 58402 .or. i == 58636 .or. i == 58639 .or. i == 58641 .or. i == 58642 .or. i == 58643) sd(i) = manual_mask
      if (i == 58644 .or. i == 58645) sd(i) = manual_mask
      
      if (i == 61302 .or. i == 61320 .or. i == 61321 .or. i == 61322 .or. i == 61325 .or. i == 61349 .or. i == 64777) sd(i) = manual_mask
      if (i == 64778 .or. i == 64809 .or. i == 66076 .or. i == 66077 .or. i == 66078 .or. i == 66079 .or. i == 66080) sd(i) = manual_mask
      if (i == 66081 .or. i == 66083 .or. i == 66084 .or. i == 66091 .or. i == 66092 .or. i == 66096 .or. i == 66099) sd(i) = manual_mask
      if (i == 66101 .or. i == 67922 .or. i == 67932 .or. i == 67935 .or. i == 67936 .or. i == 67937 .or. i == 67938) sd(i) = manual_mask
      if (i == 67939 .or. i == 67940 .or. i == 67941 .or. i == 67943 .or. i == 67945 .or. i == 67946 .or. i == 67947) sd(i) = manual_mask
      if (i == 68099 .or. i == 68100 .or. i == 68107 .or. i == 68127 .or. i == 68130 .or. i == 68132 .or. i == 68133) sd(i) = manual_mask
      if (i == 68134 .or. i == 68135 .or. i == 68138 .or. i == 68141 .or. i == 68667 .or. i == 68677 .or. i == 68690) sd(i) = manual_mask
      if (i == 68691 .or. i == 68692 .or. i == 68854 .or. i == 68855 .or. i == 68856 .or. i == 68857 .or. i == 68858) sd(i) = manual_mask
      if (i == 68859 .or. i == 68860 .or. i == 68863 .or. i == 68864 .or. i == 68867 .or. i == 68959 .or. i == 68960) sd(i) = manual_mask
      if (i == 68961 .or. i == 68977 .or. i == 68978 .or. i == 68984 .or. i == 69003 .or. i == 69138 .or. i == 69142) sd(i) = manual_mask
      if (i == 69143 .or. i == 69144 .or. i == 69146 .or. i == 69150 .or. i == 69167 .or. i == 69206 .or. i == 69232) sd(i) = manual_mask
      if (i == 69251 .or. i == 69253 .or. i == 69257) sd(i) = manual_mask
      
      if (i == 71199 .or. i == 71279 .or. i == 71280 .or. i == 71285 .or. i == 71522 .or. i == 71533 .or. i == 73616) sd(i) = manual_mask
      if (i == 74096 .or. i == 74105 .or. i == 74555 .or. i == 75025 .or. i == 75027 .or. i == 75028 .or. i == 75048) sd(i) = manual_mask
      if (i == 75389 .or. i == 75480 .or. i == 75518 .or. i == 75519 .or. i == 75879 .or. i == 75880 .or. i == 75882) sd(i) = manual_mask
      if (i == 75883 .or. i == 75884 .or. i == 75885 .or. i == 75886 .or. i == 76221 .or. i == 76222 .or. i == 76238) sd(i) = manual_mask
      if (i == 76846 .or. i == 76847 .or. i == 76849 .or. i == 76850 .or. i == 76851 .or. i == 76858 .or. i == 76873) sd(i) = manual_mask
      if (i == 76877 .or. i == 76888 .or. i == 76944 .or. i == 76948 .or. i == 76995 .or. i == 77021 .or. i == 77023) sd(i) = manual_mask
      if (i == 77024 .or. i == 77026 .or. i == 77046 .or. i == 77069 .or. i == 77072 .or. i == 77073 .or. i == 77101) sd(i) = manual_mask
      if (i == 77130 .or. i == 77133 .or. i == 77134 .or. i == 77135 .or. i == 77354 .or. i == 77363 .or. i == 77366) sd(i) = manual_mask
      if (i == 77368 .or. i == 77369 .or. i == 77373 .or. i == 77375 .or. i == 77397 .or. i == 77404 .or. i == 77533) sd(i) = manual_mask
      if (i == 77534 .or. i == 77556 .or. i == 77574 .or. i == 77579 .or. i == 77620 .or. i == 77622 .or. i == 77623) sd(i) = manual_mask
      if (i == 77647 .or. i == 77653 .or. i == 77654 .or. i == 77655 .or. i == 77691 .or. i == 77692 .or. i == 77693) sd(i) = manual_mask
      if (i == 77695 .or. i == 77804 .or. i == 77945 .or. i == 77950 .or. i == 77952 .or. i == 77953 .or. i == 77955) sd(i) = manual_mask
      if (i == 77956 .or. i == 77957 .or. i == 78876 .or. i == 78877 .or. i == 78883 .or. i == 78885 .or. i == 78886) sd(i) = manual_mask
      if (i == 78947 .or. i == 78949 .or. i == 78953 .or. i == 78955 .or. i == 78956 .or. i == 79030 .or. i == 79033) sd(i) = manual_mask
      if (i == 79035 .or. i == 79044 .or. i == 79661 .or. i == 79662 .or. i == 79663 .or. i == 79664 .or. i == 79665) sd(i) = manual_mask
      if (i == 79666 .or. i == 79667) sd(i) = manual_mask
      
      if (i == 80598 .or. i == 80599 .or. i == 80600 .or. i == 80601 .or. i == 80603 .or. i == 80604 .or. i == 80708) sd(i) = manual_mask
      if (i == 80709 .or. i == 80710 .or. i == 80711 .or. i == 80712 .or. i == 80713 .or. i == 80714 .or. i == 80715) sd(i) = manual_mask
      if (i == 80716 .or. i == 80717 .or. i == 80718 .or. i == 80719 .or. i == 80720 .or. i == 80721 .or. i == 80722) sd(i) = manual_mask
      if (i == 80723 .or. i == 80735 .or. i == 80736 .or. i == 80737 .or. i == 80738 .or. i == 80739 .or. i == 80740) sd(i) = manual_mask
      if (i == 80741 .or. i == 80742 .or. i == 80743 .or. i == 80744 .or. i == 80745 .or. i == 80746 .or. i == 80747) sd(i) = manual_mask
      if (i == 80748 .or. i == 80749 .or. i == 80750 .or. i == 80751 .or. i == 80752 .or. i == 80753 .or. i == 80754) sd(i) = manual_mask
      if (i == 80755 .or. i == 80756 .or. i == 80757 .or. i == 80758 .or. i == 80759 .or. i == 80760 .or. i == 80761) sd(i) = manual_mask
      if (i == 80762 .or. i == 80763 .or. i == 80764 .or. i == 80765 .or. i == 80766 .or. i == 80943 .or. i == 80944) sd(i) = manual_mask
      if (i == 80945 .or. i == 80946 .or. i == 80947 .or. i == 80948 .or. i == 80949 .or. i == 80958 .or. i == 80959) sd(i) = manual_mask
      if (i == 80960 .or. i == 80961 .or. i == 80962 .or. i == 80963 .or. i == 80964 .or. i == 80965 .or. i == 80966) sd(i) = manual_mask
      if (i == 80967 .or. i == 80968 .or. i == 80969 .or. i == 80970 .or. i == 81180 .or. i == 81181 .or. i == 81182) sd(i) = manual_mask
      if (i == 81183 .or. i == 81183 .or. i == 81184 .or. i == 81185 .or. i == 81186 .or. i == 81187 .or. i == 81188) sd(i) = manual_mask
      if (i == 81189 .or. i == 81190 .or. i == 81191 .or. i == 81368 .or. i == 81369 .or. i == 81370 .or. i == 81378) sd(i) = manual_mask
      if (i == 81380 .or. i == 81382 .or. i == 81421 .or. i == 81627 .or. i == 81628 .or. i == 84041 .or. i == 84049) sd(i) = manual_mask
      if (i == 84051 .or. i == 84057 .or. i == 84058 .or. i == 84063 .or. i == 84065 .or. i == 84068 .or. i == 85116) sd(i) = manual_mask
      if (i == 85174 .or. i == 86073 .or. i == 86076 .or. i == 86299 .or. i == 86301 .or. i == 86302 .or. i == 86376) sd(i) = manual_mask
      if (i == 86397 .or. i == 86398 .or. i == 86425 .or. i == 86426 .or. i == 86428 .or. i == 86437 .or. i == 86449) sd(i) = manual_mask
      if (i == 86689 .or. i == 86719 .or. i == 87284 .or. i == 87329 .or. i == 87334 .or. i == 87338 .or. i == 87339) sd(i) = manual_mask
      if (i == 87340 .or. i == 87344 .or. i == 87348 .or. i == 87349 .or. i == 87350) sd(i) = manual_mask     
    end do
  end if

end subroutine QC_snow_depth

!----------------------------------------------------------------------------------------------------
subroutine QC_snow_temperature (ndata, site_name, st1_depth, st2_depth, st3_depth, st4_depth, st5_depth, st6_depth, st1, st2, st3, st4, st5, st6, st1_mean, st2_mean, st3_mean, st4_mean, st5_mean, st6_mean)
  use subroutines_L12_13
  implicit none

  character(*), intent(in) :: site_name
  integer(8), intent(in) :: ndata
  real(8), intent(in)    :: st1_depth(ndata), st2_depth(ndata), st3_depth(ndata), st4_depth(ndata), st5_depth(ndata), st6_depth(ndata)
  real(8), intent(inout) :: st1(ndata), st2(ndata), st3(ndata), st4(ndata), st5(ndata), st6(ndata)
  real(8), intent(out)   :: st1_mean(ndata), st2_mean(ndata), st3_mean(ndata), st4_mean(ndata), st5_mean(ndata), st6_mean(ndata)

  integer(8) :: i
  integer :: iflag_before_average
  integer(8) :: referring_period

  ! New version: Sep. 13rd, 2022-
  real(8), parameter :: median_delta_st4_st1 = -0.60d0
  real(8), parameter :: median_delta_st4_st2 = -0.60d0
  real(8), parameter :: median_delta_st4_st3 = -0.20d0
  real(8), parameter :: median_delta_st4_st5 = 0.70d0
  real(8), parameter :: median_delta_st4_st6 = 0.30d0
  real(8), parameter :: SD_st4_st1 = 3.495d0
  real(8), parameter :: SD_st4_st2 = 3.390d0
  real(8), parameter :: SD_st4_st3 = 1.709d0
  real(8), parameter :: SD_st4_st5 = 2.815d0
  real(8), parameter :: SD_st4_st6 = 3.252d0

  real(8), parameter :: median_delta_st1_mean = -0.003d0
  real(8), parameter :: median_delta_st2_mean = 0.000d0
  real(8), parameter :: median_delta_st3_mean = -0.004d0
  real(8), parameter :: median_delta_st4_mean = -0.004d0
  real(8), parameter :: median_delta_st5_mean = -0.019d0
  real(8), parameter :: median_delta_st6_mean = -0.008d0
  real(8), parameter :: SD_st1_mean = 0.225d0
  real(8), parameter :: SD_st2_mean = 0.226d0
  real(8), parameter :: SD_st3_mean = 0.231d0
  real(8), parameter :: SD_st4_mean = 0.392d0
  real(8), parameter :: SD_st5_mean = 0.459d0
  real(8), parameter :: SD_st6_mean = 0.681d0

  if (site_name == 'SIGMA-A') then
    do i = 1, ndata
      if (i < 16940) then
        st5(i) = mask
        st6(i) = mask
      end if

      if (st1_depth(i) > mask_threshold .and. st1_depth(i) < -1.0d0) st1(i) = mask
      if (st2_depth(i) > mask_threshold .and. st2_depth(i) < -1.0d0) st2(i) = mask
      if (st3_depth(i) > mask_threshold .and. st3_depth(i) < -1.0d0) st3(i) = mask
      if (st4_depth(i) > mask_threshold .and. st4_depth(i) < -1.0d0) st4(i) = mask
      if (st5_depth(i) > mask_threshold .and. st5_depth(i) < -1.0d0) st5(i) = mask
      if (st6_depth(i) > -9998.0d0 .and. st6_depth(i) < 1.0d0) st6(i) = mask
      if (i >= 9373 .and. i <= 9444) then ! 3 days after maintenance 
        st3(i) = manual_mask
        st4(i) = manual_mask
      end if
      if (i >= 37273) then ! Oct 1, 2016-
        if (st4(i) > mask_threshold .and. st1(i) > mask_threshold) then
          if (abs(st4(i) - st1(i)) > (median_delta_st4_st1 + SD_st4_st1 * 3.0d0)) st1(i) = mask
        end if
        if (st4(i) > mask_threshold .and. st2(i) > mask_threshold) then
          if (abs(st4(i) - st2(i)) > (median_delta_st4_st2 + SD_st4_st2 * 3.0d0)) st2(i) = mask
        end if
        if (st4(i) > mask_threshold .and. st3(i) > mask_threshold) then
          if (abs(st4(i) - st3(i)) > (median_delta_st4_st3 + SD_st4_st3 * 1.0d0)) st3(i) = mask
        end if
        if (st4(i) > mask_threshold .and. st5(i) > mask_threshold) then
          if (abs(st4(i) - st5(i)) > (median_delta_st4_st5 + SD_st4_st5 * 2.0d0)) st5(i) = mask
        end if
        ! if (st4(i) > mask_threshold .and. st6(i) > mask_threshold) then
        !   if (abs(st4(i) - st6(i)) > (median_delta_st4_st6 + SD_st4_st6 * 3.0d0)) st6(i) = mask
        ! end if
      end if
    end do
      
    iflag_before_average = 2
    call before_average(iflag_before_average, st1, ndata, st1_mean, referring_period)
    call before_average(iflag_before_average, st2, ndata, st2_mean, referring_period)
    call before_average(iflag_before_average, st3, ndata, st3_mean, referring_period)
    call before_average(iflag_before_average, st4, ndata, st4_mean, referring_period)
    call before_average(iflag_before_average, st5, ndata, st5_mean, referring_period)
    call before_average(iflag_before_average, st6, ndata, st6_mean, referring_period)

    do i = 45313, ndata ! 1:00 on Sep. 1, 2017 -
      if (st1(i) > mask_threshold .and. st1_mean(i) > mask_threshold) then
        if (abs(st1(i) - st1_mean(i)) > (median_delta_st1_mean + SD_st1_mean * 1.0d0)) st1(i) = mask
      end if
      if (st2(i) > mask_threshold .and. st2_mean(i) > mask_threshold) then
        if (abs(st2(i) - st2_mean(i)) > (median_delta_st2_mean + SD_st2_mean * 1.0d0)) st2(i) = mask
      end if
      if (st3(i) > mask_threshold .and. st3_mean(i) > mask_threshold) then
        if (abs(st3(i) - st3_mean(i)) > (median_delta_st3_mean + SD_st3_mean * 1.0d0)) st3(i) = mask
      end if
      ! if (st4(i) > mask_threshold .and. st4_mean(i) > mask_threshold) then
      !   if (abs(st4(i) - st4_mean(i)) > (median_delta_st4_mean + SD_st5_mean * 1.0d0)) st4(i) = mask
      ! end if
      if (st5(i) > mask_threshold .and. st5_mean(i) > mask_threshold) then
        if (abs(st5(i) - st5_mean(i)) > (median_delta_st5_mean + SD_st5_mean * 1.0d0)) st5(i) = mask
      end if
      ! if (st6(i) > mask_threshold .and. st6_mean(i) > mask_threshold) then
      !   if (abs(st6(i) - st6_mean(i)) > (median_delta_st6_mean + SD_st6_mean * 1.0d0)) st6(i) = mask
      ! end if

      if (i == 49322 .or. i == 49333 .or. i == 49386 .or. i == 49427 .or. i == 49459 .or. i == 49461 .or. i == 49478 .or. &
        & i == 49479 .or. i == 49501 .or. i == 49506 .or. i == 49551 .or. i == 49576 .or. i == 49625 .or. i == 49627 .or. &
        & i == 49634) then
        st1(i) = manual_mask
      end if

      if (i == 66886 .or. i == 66893) then
        st2(i) = manual_mask
      end if

      if (i == 48225 .or. i == 48263 .or. i == 48417 .or. i == 48418 .or. i == 48433 .or. i == 48439 .or. i == 48484) st3(i) = manual_mask
      if (i == 48626 .or. i == 48808 .or. i == 48811 .or. i == 48824 .or. i == 48825 .or. i == 48859 .or. i == 48906) st3(i) = manual_mask
      if (i == 48942 .or. i == 48970 .or. i == 49029 .or. i == 49045 .or. i == 49069 .or. i == 49071 .or. i == 49095) st3(i) = manual_mask
      if (i == 49098 .or. i == 49109 .or. i == 49162 .or. i == 49218 .or. i == 49219 .or. i == 49276 .or. i == 49286) st3(i) = manual_mask
      if (i == 49289 .or. i == 49300 .or. i == 49342 .or. i == 49472 .or. i == 49497 .or. i == 49516 .or. i == 49522) st3(i) = manual_mask
      if (i == 49533 .or. i == 49613 .or. i == 49624 .or. i == 49629 .or. i == 49632 .or. i == 49636 .or. i == 49697) st3(i) = manual_mask
      if (i == 49709 .or. i == 49718 .or. i == 49723 .or. i == 49789 .or. i == 49865 .or. i == 49910 .or. i == 49928) st3(i) = manual_mask
      if (i == 49970 .or. i == 49987 .or. i == 49999) st3(i) = manual_mask
     
      if (i == 50051 .or. i == 50076 .or. i == 50269 .or. i == 50272 .or. i == 50283 .or. i == 57460 .or. i == 57487) st3(i) = manual_mask
      if (i == 57670 .or. i == 58345 .or. i == 58352 .or. i == 58776 .or. i == 58795 .or. i == 58803 .or. i == 58804) st3(i) = manual_mask
      if (i == 58805 .or. i == 58809 .or. i == 58815 .or. i == 58824 .or. i == 58832 .or. i == 58835 .or. i == 58994) st3(i) = manual_mask
      if (i == 59105 .or. i == 59121 .or. i == 59134 .or. i == 59141 .or. i == 59143 .or. i == 59155 .or. i == 59157) st3(i) = manual_mask
      if (i == 59451 .or. i == 59618) st3(i) = manual_mask
      
      if (i == 64898 .or. i == 64989 .or. i == 65185 .or. i == 65347 .or. i == 65370 .or. i == 65459 .or. i == 65465) st3(i) = manual_mask
      if (i == 65467 .or. i == 65505 .or. i == 65625 .or. i == 65652 .or. i == 65679 .or. i == 65685 .or. i == 65697) st3(i) = manual_mask
      if (i == 65765 .or. i == 65773 .or. i == 66083 .or. i == 66169 .or. i == 66639 .or. i == 66467 .or. i == 66510) st3(i) = manual_mask
      if (i == 67190 .or. i == 67252 .or. i == 67411 .or. i == 68183 .or. i == 68715 .or. i == 69927 .or. i == 69929) st3(i) = manual_mask
      if (i == 69937) st3(i) = manual_mask
     
      if (i == 57926 .or. i == 58146 .or. i == 64409 .or. i == 69375 .or. i == 69379 .or. i == 69412 .or. i == 69413) st5(i) = manual_mask
      if (i == 69424) st5(i) = manual_mask
    end do
  end if

  if (site_name == 'SIGMA-B') then
    st1 = mask
    st2 = mask
    st3 = mask
    st4 = mask
    st5 = mask
    st6 = mask
  end if

end subroutine QC_snow_temperature

!----------------------------------------------------------------------------------------------------
subroutine QC_air_press (ndata, site_name, press)
  use subroutines_L12_13
  implicit none

  character(*), intent(in) :: site_name
  integer(8), intent(in)   :: ndata
  real(8), intent(inout)   :: press(ndata)

  integer(8) :: i
  real(8)    :: delta_press(ndata)
  
  integer(8), parameter :: iflag_refer_former_parameter = 2

  call refer_formar_parameter(iflag_refer_former_parameter, press, ndata, delta_press)
  if (site_name == 'SIGMA-A') then
    do i = 1, ndata
      if (i == 49076 .or. i == 65063) press(i) = manual_mask
    end do
  end if

end subroutine QC_air_press
!==========================================
!>>> subroutines end
!==========================================
end module QC_L12_13


!>>> Revision records
! ver.1.0    (date: Sep. 8, 2022) subroutine snow_temperature was added
