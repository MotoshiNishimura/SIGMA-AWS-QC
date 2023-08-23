!=============================================================================
! Main program for making hourly SIGMA-AWS data
! This program performs Quality Control (QC) to mask the abnormal data and creates Level 1.2 dataset from Level 1.1 dataset.
! Masashi Niwano (@ Meteorological Research Institute, Japan) created the base of this program
! and Motoshi Nishimura (@ National Institute of Polar Research, Japan) developed it.
! Motoshi Nishimura @ national Institute of Polar Research
!=============================================================================
program makedata
  use data_input_output
  use sol_info, only : solar
  use QC_L11_12
  implicit none

  !==============
  !>>> variables
  !==============
  character(13)  :: input_info
  character(7)   :: site_name
  character(5)   :: cyear
  character(2)   :: cmonth
  character(2)   :: cday
  character(22)  :: time_stamp
  character(100) :: lat
  character(100) :: lon
  character(100) :: time
  character(100) :: ymd
  character(11), allocatable :: ctime(:)
  character(22), allocatable :: cdate_time(:)
  character(50), allocatable :: cell(:,:)
  character(:), allocatable  :: infile
  character(:), allocatable  :: outfile1
  
  integer(8) :: i
  integer(8) :: ndata
  integer(8) :: column_number
  integer(8), allocatable :: year(:)
  integer(8), allocatable :: month(:)
  integer(8), allocatable :: day(:)

  ! input data
  real(8) :: solz
  real(8) :: sola
  real(8), allocatable :: atemp1(:), atemp2(:)
  real(8), allocatable :: arh1(:), arh2(:)
  real(8), allocatable :: ws1(:), ws2(:)
  real(8), allocatable :: wd1(:), wd2(:)
  real(8), allocatable :: swd(:), swu(:)
  real(8), allocatable :: lwd(:), lwu(:)
  real(8), allocatable :: nird(:), niru(:)
  real(8), allocatable :: press(:)
  real(8), allocatable :: snow_depth(:), sr50raw(:)
  real(8), allocatable :: sensor_height(:)
  real(8), allocatable :: st1(:), st2(:), st3(:), st4(:), st5(:), st6(:)
  real(8), allocatable :: cnr4t(:)
  real(8), allocatable :: tiltx(:), tilty(:)
  real(8), allocatable :: initt(:), bat(:)
  real(8), allocatable :: qual(:)
  real(8), allocatable :: solar_zenith(:), solar_azimuth(:)
  
  ! output data
  real(8), allocatable :: SW_TOA(:)
  real(8), allocatable :: swd_slope(:)
  real(8), allocatable :: Id_slope(:)
  real(8), allocatable :: Is_slope(:)
  real(8), allocatable :: solar_zenith_slope(:)
  real(8), allocatable :: diffu_ratio(:)
  
  integer(8), parameter :: param_num = 100
  
  !>>> Initial setting
  ! input_info: data period and site name >>> 12-22-SIGMA-B => 2012-2022, SIGMA-B site
  read (*,*) input_info
  print *, 'input data Information: ', input_info
  cyear = input_info(:5)
  site_name = input_info(7:)
  infile = cyear // '.csv'
  
  !>>> set output file name & site
  if (site_name == 'SIGMA-A') then
    outfile1 = site_name // '_' // cyear // '_met_Lv1.2.csv'
    lat = '78:03:06'
    lon = '-292:-22:-18' ! [E] == [67:37:42 W]
  else if (site_name == 'SIGMA-B') then
    outfile1 = site_name // '_' // cyear // '_met_Lv1.2.csv'
    lat = '77:31:06'
    lon = '-290:-56:-17' ! [E] = [69:03:43 W]
  else
    stop
  end if
  
  print *, '=========================================================================='
  print *, '  Lunched Program Information                                             '
  print *, '=========================================================================='
  print *, '[Observation site and Data Level]'
  if (site_name == 'SIGMA-A') print *, '>>> Observation site: SIGMA-A'
  if (site_name == 'SIGMA-B') print *, '>>> Observation site: SIGMA-B' 
  print *, '>>> Making Level 1.2: Initial Control (easy Quality Control)'
  print *, '----------------------------------------------------------------------' 
  
  call data_number (infile, ndata)
  
  allocate(character(50) :: cell(ndata+1,param_num))
  allocate(cdate_time(ndata), ctime(ndata))
  allocate(year(ndata), month(ndata), day(ndata))
  allocate(atemp1(ndata), arh1(ndata), ws1(ndata), wd1(ndata))
  allocate(atemp2(ndata), arh2(ndata), ws2(ndata), wd2(ndata))
  allocate(swd(ndata), swu(ndata), lwd(ndata), lwu(ndata), nird(ndata), niru(ndata))
  allocate(snow_depth(ndata), press(ndata), sensor_height(ndata))
  allocate(st1(ndata), st2(ndata), st3(ndata), st4(ndata), st5(ndata), st6(ndata))
  allocate(cnr4t(ndata), tiltx(ndata), tilty(ndata), initt(ndata), bat(ndata))
  allocate(sr50raw(ndata), qual(ndata))
  allocate(solar_zenith(ndata), solar_azimuth(ndata))
  allocate(SW_TOA(ndata), swd_slope(ndata), Id_slope(ndata), Is_slope(ndata), solar_zenith_slope(ndata), diffu_ratio(ndata))
  
  !>>> reading data file
  call data_input (infile, cdate_time, &
                 & ws1, wd1, atemp1, arh1, swd, &
                 & swu, lwd, lwu, nird, niru, &
                 & snow_depth, press, st1, st2, st3, &
                 & cnr4t, tiltx, tilty, initt, bat, &
                 & ws2, wd2, atemp2, arh2, st4, &
                 & st5, st6, sr50raw, qual)
  !==========================================================================================
  print *, '--------------------------------------------------------------' 
  print *,  '[Main Program] Main procedures begin.'
  
  !>>> setting time stamps
  do i = 1, ndata
    time_stamp = cdate_time(i)
    cyear = time_stamp(1:4)
    cmonth = time_stamp(6:7)
    cday = time_stamp(9:10)
    ctime(i) = time_stamp(12:16)
    
    read(cyear,'(I4.4)') year(i)
    read(cmonth,'(I2.2)') month(i)
    read(cday,'(I2.2)') day(i)
    
    ymd = trim(cyear) // trim(cmonth) // trim(cday)
    time = trim(ctime(i)) // ':00'

    if (time(2:2) == ':') then
      time = '0' // trim(time)
    end if

    !>>> calculate solar position information     
    call solar(ymd, time, lat, lon, solz, sola)
    solar_zenith(i) = solz
    solar_azimuth(i) = sola
    call SW_TOA_calculation (solar_zenith(i), year(i), month(i), day(i), SW_TOA(i))
    call slope_correction (swd(i), solar_zenith(i), solar_azimuth(i), SW_TOA(i), swd_slope(i), Id_slope(i), Is_slope(i), &
&                          solar_zenith_slope(i), diffu_ratio(i))
     
  end do

  !----------------------------------------------------------------------------
  !================================
  ! Initial QC
  !================================
  call QC_wind_speed (ndata, site_name, ws1, ws2)
  call QC_wind_direction (ndata, ws1, ws2, wd1, wd2)
  call wind_direction_correction (ndata, site_name, wd1, wd2)
  call QC_air_temperature (ndata, site_name, atemp1, atemp2)
  call QC_air_humidity (ndata, arh1, arh2)
  call QC_shortwave_radiation (ndata, site_name, SW_TOA, solar_zenith, solar_zenith_slope, swd, swu)
  call QC_nir_infrared_radiation (ndata, site_name, SW_TOA, solar_zenith, solar_zenith_slope, nird, niru)
  call QC_longwave_radiation (ndata, site_name, lwd, lwu)
  if (site_name == 'SIGMA-A') call QC_snow_depth (ndata, site_name, atemp2, snow_depth)
  if (site_name == 'SIGMA-B') call QC_snow_depth (ndata, site_name, atemp1, snow_depth)
  call QC_air_press (ndata, site_name, press)
  if (site_name == 'SIGMA-A') call QC_snow_temperature (ndata, st1, st2, st3, st4, st5, st6)
  
  print *, '[Main Program:] All QC procedures... Done'

  !----------------------------------------------------------------------------
  !>>> Level 1.2 data set output
  call data_writing (site_name, ndata, cell, cdate_time, &
                   & ws1, wd1, atemp1, arh1, ws2, &
                   & wd2, atemp2, arh2, swd, swu, &
                   & lwd, lwu, nird, niru, snow_depth, &
                   & press, st1, st2, st3, st4, &
                   & st5, st6, solar_zenith, solar_azimuth, column_number) 
              
  call data_output (ndata, column_number, outfile1, cell)  
  !----------------------------------------------------------------------------

  deallocate(cell)
  deallocate(cdate_time, ctime)
  deallocate(year, month, day)
  deallocate(atemp1, arh1, ws1, wd1)
  deallocate(atemp2, arh2, ws2, wd2)
  deallocate(swd, swu, lwd, lwu, nird, niru)
  deallocate(press, snow_depth, sensor_height)
  deallocate(st1, st2, st3, st4, st5, st6)
  deallocate(cnr4t, tiltx, tilty, initt, bat)
  deallocate(sr50raw, qual)
  deallocate(solar_zenith, solar_azimuth)
  deallocate(SW_TOA, swd_slope, Id_slope, Is_slope, solar_zenith_slope, diffu_ratio)

  print *, '[Main Program] >>> Program end ...'
  print *, '==============================================================================='

end program makedata

!>>> subroutines below
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
  
  real(8), parameter :: mask = -9999.0d0           ! constant value repredenting data mask
  real(8), parameter :: mask_threshold = -8887.0d0 ! mask or mannual mask threthold value


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

!----------------------------------------------------------------------------------------------------
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
