!=======================================================================
! This module summarizes the subroutines used in the Quality Control process (Secondary Control) of SIGMA-AWS observation data.
! Motoshi Nishimura @ National Institute of Polar Research
!=======================================================================
module subroutines_L12_13
  implicit none

  real(8), parameter :: abstemp = 273.15d0 ! absolute temperature of 0 degC
  real(8), parameter :: eps = 0.98d0       ! snow emissivity [no dimension]
  real(8), parameter :: sigma = 5.67d-8    ! Stefan-Bolzmann constant[W/m^2/K^4]  

  real(8), parameter :: mask = -9999.0d0           ! constant value repredenting data mask
  real(8), parameter :: manual_mask = -8888.0d0    ! constant value repredenting data manual mask
  real(8), parameter :: mask_threshold = -8887.0d0 ! mask or mannual mask threthold value  
  
  contains

!==========================================
!>>> subroutines below
!==========================================
subroutine slope_correction (swd, solz, sola, SW_TOA, swd_slope, Id_slope, Is_slope, solz_slope, diffu_ratio)
  implicit none
 
  real(8), intent(in)  :: swd
  real(8), intent(in)  :: solz
  real(8), intent(in)  :: sola
  real(8), intent(in)  :: SW_TOA
  real(8), intent(out) :: swd_slope
  real(8), intent(out) :: Id_slope
  real(8), intent(out) :: Is_slope
  real(8), intent(out) :: solz_slope
  real(8), intent(out) :: diffu_ratio

  real(8) :: trans_ratio
  real(8) :: solar_zenith_slope_radian
  
  real(8), parameter :: pi = atan(1.0d0) * 4.0d0 !circular constant [no dimention]  
  real(8), parameter :: slope_angle = 4.019799d0
  real(8), parameter :: slope_azm = 220.33951d0


  solar_zenith_slope_radian = &
&   acos(cos(slope_angle / 180.0d0 * pi) * cos(solz / 180.0d0 * pi) + &
&   sin(slope_angle / 180.0d0 * pi) * sin(solz / 180.0d0 * pi) * cos((sola - slope_azm) / 180.0d0 * pi))
  
  solz_slope = solar_zenith_slope_radian * 180.0 / pi
  
  if (swd > mask_threshold) then
    trans_ratio = swd / SW_TOA

    ! Hock and Holmgren (2005)
    if (trans_ratio <= 0.15d0) then
      diffu_ratio = 1.0d0
    else if (trans_ratio > 0.15d0 .and. trans_ratio < 0.8d0) then
      diffu_ratio = 0.929d0 + 1.134d0 * trans_ratio - 5.111d0 * trans_ratio ** 2.0d0 + 3.106d0 * trans_ratio ** 3.0d0
    else if (trans_ratio >= 0.8d0) then
      diffu_ratio = 0.15d0
    end if

    ! Jonsell et al. (2003)
    if (solz < 90.0d0 .and. solz_slope < 90.0d0) then
      Id_slope = swd * (1.0d0 - diffu_ratio) * cos(solar_zenith_slope_radian) / cos(solz / 180.0d0 * pi)
    else
      Id_slope = 0.0d0
    end if
    Is_slope = swd * diffu_ratio
    swd_slope = Id_slope + Is_slope
  
  else
    Id_slope = mask
    Is_slope = mask
    swd_slope = mask
  end if

end subroutine slope_correction

!------------------------------------------------------------------------
! subroutine   : Shortwave radiation for top of atmosphere
! version      : 1.0
! subject      : calculating SW_TOA
! editer       : Motoshi Nishimura @ National Institute of Polar Research
!------------------------------------------------------------------------
subroutine SW_TOA_calculation (solar_zenith, year, month, day, SW_TOA)
  implicit none

  integer(8), intent(in) :: year, month, day
  real(8), intent(in)    :: solar_zenith
  real(8), intent(out)   :: SW_TOA

  integer(8) :: i, j
  integer(8) :: nday(12)
  integer(8) :: day_months
  integer(8) :: day_of_year
  integer(8) :: day_number
  real(8)    :: std_time
  real(8)    :: phi
  real(8)    :: d
  
  real(8), parameter :: pi = atan(1.0d0) * 4.0d0 !circular constant [no dimention]  
  real(8), parameter :: solar_const = 1361.0d0
  real(8), parameter :: annual_days = 365.25d0
  real(8), parameter :: eccentricity = 0.01637d0
 
  data (nday(i), i = 1, 12) / 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /

  if (mod(year, 4) == 0) nday(2) = 29
  if (mod(year, 100) == 0) nday(2) = 28
  if (mod(year, 400) == 0) nday(2) = 29
 
  day_months = 0
  if (month >= 2) then
    do j = 2, month
      day_months = day_months + nday(j-1)  
    end do
  end if

  day_of_year =  day + day_months    
  day_number = day_of_year - 2    ! the day from the perihelion (assumed as 3, January)

  std_time = day_number / annual_days 
  phi = 2.0d0 * pi * std_time + 2.0d0 * eccentricity * sin(2.0d0 * pi * std_time)
  d = (1.0d0 + eccentricity * cos(phi)) / (1.0d0 - eccentricity ** 2.0d0)
 
  if (solar_zenith > 90.0d0) then
    SW_TOA = 0.0d0
  else 
    SW_TOA = solar_const * d ** 2.0d0 * cos(solar_zenith / 180.0d0 * pi)
  end if

end subroutine SW_TOA_calculation

!------------------------------------------------------------------------
! subroutine   : accumulated albedo
! version      : 1.0
! subject      : calculating daily accumulated albedo
!                Daily accumulated albedo is an albedo calculated from the integrated shortwave radiation over the past 24 hours
!                at any given time.
! reference    : van den Broeke et al. (2004, JGR Atmospheres)
! editer       : Motoshi Nishimura @ National Institute of Polar Research
!------------------------------------------------------------------------
subroutine accumulated_albedo (ndata, swd, swu, acc_albedo)
  implicit none

  integer(8), intent(in) :: ndata
  real(8), intent(in)    :: swd(ndata)
  real(8), intent(in)    :: swu(ndata)
  real(8), intent(out)   :: acc_albedo(ndata)
 
  integer(8) :: i, ii
  integer(8) :: mask_num
  real(8)    :: swd_tmp(ndata)
  real(8)    :: swu_tmp(ndata)
 
  do i = 1, ndata
    swd_tmp(i) = swd(i)
    swu_tmp(i) = swu(i)
  end do

  do i = 1, 23
    if (i == 1) acc_albedo(i) = swu(i) / swd(i)
    if (i >= 2) then
      mask_num = 0
      if (swd(i) <= 0.0d0 .or. swu(i) <= 0.0d0) then
        mask_num = mask_num + 1
        swd_tmp(i) = 0.0d0
        swu_tmp(i) = 0.0d0
      end if

      if (mask_num >= 12) then
        acc_albedo(i) = mask
      else
        if (sum(swd_tmp(1:i)) > 0.0d0) then
          acc_albedo(i) = sum(swu_tmp(1:i)) / sum(swd_tmp(1:i))
        else
          acc_albedo(i) = mask
        end if
      end if
    end if
  end do

  do i = 24, ndata
    mask_num = 0
    do ii = i-23, i
      if (swd(ii) <= 0.0d0 .or. swu(ii) <= 0.0d0) then
        mask_num = mask_num + 1
        swd_tmp(ii) = 0.0d0
        swu_tmp(ii) = 0.0d0
      end if
    end do
 
    if (mask_num >= 12) then
      acc_albedo(i) = mask
    else
      if (sum(swd_tmp(i-23:i)) > 0.0d0) then
        acc_albedo(i) = sum(swu_tmp(i-23:i)) / sum(swd_tmp(i-23:i))
      else
        acc_albedo(i) = mask
      end if
    end if
  end do

end subroutine accumulated_albedo

!------------------------------------------------------------------------
! subroutine   : standard longwave radiation
! version      : 1.0
! subject      : calculating atmospheric emission of longwave radiation using air temperature
! reference    : Brock and Arnold (2000)
! editer       : Motoshi Nishimura @ National Institute of Polar Research
!------------------------------------------------------------------------
function standard_longwave_radiation (atemp, cloud_cover) result(lw_standard)
  implicit none

  real(8), intent(in) :: atemp
  real(8), intent(in) :: cloud_cover
  real(8) :: lw_standard

  real(8) :: emiss_clearsky
  real(8) :: emiss_air
  
  real(8), parameter :: cloud_type_constant = 0.26d0 ! Braithwaite and Olsen (1990)

  emiss_clearsky = 8.733d0 * 1.0d-3 * (atemp + abstemp) ** 0.788d0 ! Ohmura (1981)
  emiss_air = (1.0d0 + cloud_type_constant * cloud_cover) * emiss_clearsky
  lw_standard = emiss_air * sigma * (atemp + abstemp) ** 4.0d0


end function standard_longwave_radiation

!------------------------------------------------------------------------
! subroutine   : snow surface temperature
! version      : 1.0
! subject      : calculating surface temperature using longwave radiation
! editer       : Motoshi Nishimura @ National Institute of Polar Research
!------------------------------------------------------------------------
function snow_surface_temperature (lwd, lwu) result(sst)
  implicit none

  real(8), intent(in) :: lwd, lwu
  real(8) :: sst

  if (lwd > 0.0d0 .and. lwu > 0.0d0) then
    sst = ((lwu - (1.0d0 - eps) * lwd) / (eps * sigma)) ** 0.25d0
    sst = sst - abstemp

    if (sst > 0.0d0) then
      sst = 0.0d0
    end if
  else
    sst = mask
  end if

end function snow_surface_temperature
!------------------------------------------------------------------------
! subroutine   : sensor height
! version      : 1.0
! subject      : calculating set height of the air temperature and humidity and wind speed sensors above surface
! editer       : Motoshi Nishimura @ National Institute of Polar Research
!------------------------------------------------------------------------
subroutine sensor_height (ndata, site_name, snow_depth, sensor_height1, sensor_height2)
  implicit none

  character(*), intent(in) :: site_name
  integer(8), intent(in)   :: ndata
  real(8), intent(in)      :: snow_depth(ndata)
  real(8), intent(out)     :: sensor_height1(ndata), sensor_height2(ndata)

  integer(8) :: i
  real(8), parameter :: pole_length1 = 130.0d0
  real(8), parameter :: pole_length2 = 147.0d0
  real(8), parameter :: pole_length3 = 130.0d0

  do i = 1, ndata
    if (site_name == 'SIGMA-A') then
      sensor_height1(1) = 3.0d0

      if (snow_depth(i) > mask_threshold) then
        sensor_height1(i) = sensor_height1(1) - snow_depth(i) * 1.0d-2
        ! 11:00 on Jul. 25, 2013 -
        if (i >= 9347 .and. i <= ndata) sensor_height1(i) = sensor_height1(i) + pole_length1 * 1.0d-2
        ! 20:00 on Jun. 6, 2014
        if (i >= 16940 .and. i <= ndata) sensor_height1(i) = sensor_height1(i) + pole_length2 * 1.0d-2
        ! 18:00 on May 26, 2017 -
        if (i >= 42978 .and. i <= ndata) sensor_height1(i) = sensor_height1(i) + pole_length3 * 1.0d-2
      else
        sensor_height1(i) = mask
      end if

      if (sensor_height1(i) > mask_threshold) then
        sensor_height2(i) = sensor_height1(i) + 3.0d0
      else
        sensor_height2(i) = mask
      end if
    end if

    if (site_name == 'SIGMA-B') then
      sensor_height1(1) = 3.0d0
      sensor_height2(i) = mask

      if (snow_depth(i) > mask_threshold) then
        sensor_height1(i) = sensor_height1(1) - snow_depth(i) * 1.0d-2
      else
        sensor_height1(i) = mask
      end if
    end if
  end do

end subroutine sensor_height
!------------------------------------------------------------------------
! subroutine   : sensor depth of snow temperature sensors
! version      : 1.0
! subject      : calculating snow temperature senseor depth in the snowpack
! editer       : Motoshi Nishimura @ National Institute of Polar Research
!------------------------------------------------------------------------
subroutine snow_temperature_sensor_depth (ndata, site_name, snow_depth, st1_depth, st2_depth, st3_depth, st4_depth, st5_depth, st6_depth)
  implicit none

  character(*), intent(in) :: site_name
  integer(8), intent(in)   :: ndata
  real(8), intent(in)      :: snow_depth(ndata)
  real(8), intent(out)     :: st1_depth(ndata), st2_depth(ndata), st3_depth(ndata)
  real(8), intent(out)     :: st4_depth(ndata), st5_depth(ndata), st6_depth(ndata)

  integer(8) :: i
  real(8)    :: st1_depth_tmp(ndata), st2_depth_tmp(ndata), st3_depth_tmp(ndata)
  real(8)    :: st4_depth_tmp(ndata), st5_depth_tmp(ndata), st6_depth_tmp(ndata)

  integer(8), parameter :: maintenance_time1 = 9375  ! 26, July 2013 15:00-
  integer(8), parameter :: maintenance_time2 = 16940 !  6, June 2014 20:00-
  real(8), parameter    :: st_depth_threshold = -1.0d0

  do i = 1, ndata ! 29, June 2012 19:00 initial setting
    if (site_name == 'SIGMA-A') then
      if (snow_depth(i) > mask_threshold) then
        st1_depth(i) = snow_depth(i) + 100.0d0 - snow_depth(1) ! initial depth = 100 cm
        st2_depth(i) = snow_depth(i) + 70.0d0 - snow_depth(1) ! initial depth = 70 cm
        st3_depth(i) = snow_depth(i) + 40.0d0 - snow_depth(1) ! initial depth = 40 cm
        st4_depth(i) = snow_depth(i) + 5.0d0 - snow_depth(1) ! initial depth = 5 cm

        if (i >= maintenance_time1) then 
          st3_depth(i) = snow_depth(i) + 46.0d0 - snow_depth(maintenance_time1) ! initial depth = 40 cm
          st4_depth(i) = snow_depth(i) + 16.0d0 - snow_depth(maintenance_time1) ! initial depth = 5 cm
        end if  
        if (i >= maintenance_time2) then 
          st5_depth(i) = snow_depth(i) + 5.0d0 - snow_depth(maintenance_time2) ! initial depth = 5 cm
          st6_depth(i) = snow_depth(i) - 45.0d0 - snow_depth(maintenance_time2) ! initial depth = -45 cm (45 cm above from the surface)
        else
          st5_depth(i) = mask
          st6_depth(i) = mask
        end if
      else
        st1_depth(i) = mask
        st2_depth(i) = mask
        st3_depth(i) = mask
        st4_depth(i) = mask
        st5_depth(i) = mask
        st6_depth(i) = mask
      end if

      st1_depth_tmp(i) = st1_depth(i)
      st2_depth_tmp(i) = st2_depth(i)
      st3_depth_tmp(i) = st3_depth(i)
      st4_depth_tmp(i) = st4_depth(i)
      st5_depth_tmp(i) = st5_depth(i)
      st6_depth_tmp(i) = st6_depth(i)

      if (i >= maintenance_time2) then
        if (st6_depth(i) < st_depth_threshold) st6_depth(i) = -9997.0d0
      end if
    end if
  end do

  if (site_name == 'SIGMA-B') then
    st1_depth = mask
    st2_depth = mask
    st3_depth = mask
    st4_depth = mask
    st5_depth = mask
    st6_depth = mask
  end if

  call flag_snow_sensor_depth (ndata, site_name, maintenance_time2, st6_depth_tmp, st6_depth)

end subroutine snow_temperature_sensor_depth

!----------------------------------------------------------------------------------------------------
subroutine flag_snow_sensor_depth (ndata, site_name, maintenance_time, st_depth_tmp, st_depth)
  implicit none

  character(*), intent(in) :: site_name
  integer(8), intent(in)   :: ndata
  integer(8), intent(in)   :: maintenance_time
  real(8), intent(in)      :: st_depth_tmp(ndata)
  real(8), intent(inout)   :: st_depth(ndata)

  integer(8) :: i
  integer(8) :: mask_count

  real(8), parameter :: st_depth_threshold = -1.0d0
  real(8), parameter :: st_flag = -9997.0d0 ! a flag indicating snow sensor is suspected to be located above, not below the surface

  if (site_name == 'SIGMA-A') then
    i = maintenance_time
    do
      mask_count = 0
      do
        if (st_depth_tmp(i+mask_count) < mask_threshold) then
          mask_count = mask_count + 1
         if (i + mask_count == ndata) exit
        else
          exit
        end if
      end do

      if (mask_count > 0) then
        ! snow temperature sensor is positioned under the surface
        if (st_depth_tmp(i+mask_count) > st_depth_threshold .and. st_depth_tmp(i-1) > st_depth_threshold) then
          st_depth(i:i+mask_count-1) = mask

        ! not be able to distinguish snow temperature sensor is positioned under the surface or not
        else if (st_depth_tmp(i+mask_count) > st_depth_threshold .and. st_depth_tmp(i-1) < st_depth_threshold) then
          st_depth(i:i+mask_count-1) = st_flag
        
        ! snow temperature sensor is positioned above the surface
        else if (st_depth_tmp(i+mask_count) < st_depth_threshold .and. st_depth_tmp(i-1) < st_depth_threshold) then
          st_depth(i:i+mask_count-1) = st_flag
          if (st_depth_tmp(i+mask_count) < mask_threshold) st_depth(i:i+mask_count) = mask
        
        ! not be able to distinguish snow temperature sensor is positioned under the surface or not
        else if (st_depth_tmp(i+mask_count) < st_depth_threshold .and. st_depth_tmp(i-1) > st_depth_threshold) then
          st_depth(i:i+mask_count-1) = st_flag
          if (st_depth_tmp(i+mask_count) < mask_threshold) st_depth(i:i+mask_count) = mask
        end if
        i = i + mask_count
      else
        i = i + 1
      end if
      if (i == ndata) exit
    end do
  end if
end subroutine flag_snow_sensor_depth

!------------------------------------------------------------------------
! subroutine   : before_average
! version      : 1.0 (October 2, 2020)
! version      : 1.1 referring period was changed (October 5, 2020)
! subject      : calculating averages over a period of time going back an arbitrary time step from the data at a given time
! editer       : Motoshi Nishimura @ National Institute of Polar Research
!------------------------------------------------------------------------
subroutine before_average(iflag_before_average, rparam, ndata, mean, referring_period)
  implicit none

  integer, intent(in)     :: iflag_before_average
  integer(8), intent(in)  :: ndata
  real(8), intent(in)     :: rparam(ndata)
  integer(8), intent(out) :: referring_period
  real(8), intent(out)    :: mean(ndata)

  integer(8) :: i, ii
  integer(8) :: masking_threshold
  integer(8) :: mask_num
  real(8)    :: rparam_tmp(ndata)
 
  if (iflag_before_average == 1) then ! for snow_depth
    referring_period = 72
    masking_threshold = 48
  else if (iflag_before_average == 2) then ! for snow_temp
    referring_period = 72
    masking_threshold = 48
  else if (iflag_before_average == 3) then ! for press
    referring_period = 8
    masking_threshold = 6
  else if (iflag_before_average == 4) then ! for 2nd snow depth
    referring_period = 72
    masking_threshold = 60
  else if (iflag_before_average == 5) then ! for atmospheric relative humidity
    referring_period = 3
    masking_threshold = 2
  end if

  do i = 1, ndata
    rparam_tmp(i) = rparam(i)
  end do

  do i = referring_period, ndata
    mask_num = count(rparam(i-(referring_period-1):i) < mask_threshold)

    do ii = i-(referring_period-1), i
      if (rparam_tmp(ii) < mask_threshold) rparam_tmp(ii) = 0.0d0
    end do

    if (mask_num >= masking_threshold) then
      mean(i) = mask
    else
      mean(i) = sum(rparam_tmp(i-(referring_period-1):i)) / (referring_period - mask_num)
    end if
  end do

end subroutine before_average

!------------------------------------------------------------------------
! subroutine   : refering before parameter
! version      : 3.0
! subject      : calculating the diferenece between data(i) - data(i-x)
! editer       : Motoshi Nishimura @ National Institute of Polar Research
!------------------------------------------------------------------------
subroutine refer_formar_parameter(param_flag, rparam, ndata, delta_rparam)
  implicit none

  integer(8), intent(in) :: ndata
  integer(8), intent(in) :: param_flag
  real(8), intent(inout) :: rparam(ndata) 
  real(8), intent(out)   :: delta_rparam(ndata)

  integer(8) :: i, ii, iii, iv
  integer(8) :: mask_count1, mask_count2
  real(8)    :: rparam_tmp(ndata)
  real(8)    :: rparam_mean(ndata)!, delta_rparam(ndata)
  real(8)    :: delta_threshold
  
  integer(8), parameter :: referring_period = 3!, masking_threshold = 2

  if (param_flag == 1) delta_threshold = 30.0d0 ! for air humidity
  if (param_flag == 2) delta_threshold = 20.0d0 ! for atmospheric pressure
 
  i = 2
  do
    rparam_tmp(i) = rparam(i)

    ii = 0
    do
      if (rparam(i+ii) < mask_threshold) then
        ii = ii + 1
        if (i + ii == ndata) exit
      else
        exit
      end if
    end do
    mask_count1 = ii

    iii = 1 
    do
      if (rparam(i+ii-iii) < mask_threshold) then
        iii = iii + 1
      else
        exit
      end if
    end do
  
    if ((i+ii-iii) > 4) then
      mask_count2 = 0
      do iv = 1, referring_period
        if (rparam(i+ii-iii-(iv-1)) < mask_threshold) then
          mask_count2 = mask_count2 + 1
          rparam_tmp(i+ii-iii-(iv-1)) = 0.0d0
        end if
      end do
      rparam_mean(i+ii) = sum(rparam_tmp(i+ii-iii-(referring_period-1):i+ii-iii)) / (referring_period-mask_count2)
    end if

    if ((i+ii-iii) > 4) then
      delta_rparam(i+ii) = rparam(i+ii) - rparam_mean(i+ii)
      if (abs(delta_rparam(i+ii)) > delta_threshold) rparam(i+ii) = mask
    end if

    if (ii > 1) then
      i = i + ii
    else
      i = i + 1
    end if
    if (i == ndata .or. i+mask_count1 == ndata) exit
  end do

end subroutine refer_formar_parameter
!==========================================
!>>> subroutines end
!==========================================
end module subroutines_L12_13
