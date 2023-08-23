!=============================================================================
! Main program for making hourly SIGMA-AWS data
! This program performs Quality Control (QC) to mask the abnormal data and creates Level 1.3 dataset from Level 1.2 dataset.
! Level 1.3 dataset shall be the dataset from which all outliers have been masked by the QC.
! Masashi Niwano (@ Meteorological Research Institute, Japan) created the base of this program
! and Motoshi Nishimura (@ National Institute of Polar Research, Japan) developed it.
! Motoshi Nishimura @ national Institute of Polar Research
!=============================================================================
program makedata
  use data_input_output
  use subroutines_L12_13
  use QC_L12_13
  implicit none

  !==============
  !>>> variables
  !==============
  character(13) :: input_info
  character(7)  :: site_name  
  character(5)  :: cyear
  character(2)  :: cmonth
  character(2)  :: cday
  character(22) :: time_stamp
  character(5)  :: time  
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
  real(8), allocatable :: atemp1(:), atemp2(:)
  real(8), allocatable :: arh1(:), arh2(:)
  real(8), allocatable :: ws1(:), ws2(:)
  real(8), allocatable :: wd1(:), wd2(:)
  real(8), allocatable :: swd(:), swu(:)
  real(8), allocatable :: lwd(:), lwu(:)
  real(8), allocatable :: nird(:), niru(:)
  real(8), allocatable :: press(:)
  real(8), allocatable :: snow_depth(:)
  real(8), allocatable :: sensor_height1(:), sensor_height2(:)
  real(8), allocatable :: st1(:), st2(:), st3(:), st4(:), st5(:), st6(:)
  real(8), allocatable :: solz(:), sola(:)

  ! output data
  real(8), allocatable :: albedo(:), albedo_acc(:)
  real(8), allocatable :: nir_albedo(:), nir_albedo_acc(:)
  real(8), allocatable :: nir_frac(:)
  real(8), allocatable :: sst(:)
  real(8), allocatable :: SW_TOA(:)
  real(8), allocatable :: swd_slope(:)
  real(8), allocatable :: Id_slope(:)
  real(8), allocatable :: Is_slope(:)
  real(8), allocatable :: solz_slope(:)
  real(8), allocatable :: diffu_ratio(:)  
  real(8), allocatable :: nir_max(:)
  real(8), allocatable :: st1_depth(:), st2_depth(:), st3_depth(:), st4_depth(:), st5_depth(:), st6_depth(:)
  real(8), allocatable :: st1_mean(:), st2_mean(:), st3_mean(:), st4_mean(:), st5_mean(:), st6_mean(:)  
  real(8), allocatable :: sd_mean(:)
  real(8), allocatable :: press_mean(:)
  real(8), allocatable :: delta_press(:)
  real(8), allocatable :: lw_standard(:)
    
  integer(8), parameter :: param_num = 100
  
  !>>> Initial setting
  ! input_info: data period and site name >>> 12-22-SIGMA-B => 2012-2022, SIGMA-B site
  read (*,*) input_info
  print *, 'input data Information: ', input_info
  cyear = input_info(:5)
  site_name = input_info(7:)
  infile = cyear // '.csv'

  !>>> set output file name & site
  print *, '=========================================================================='
  print *, '  Lunched Program Information                                             '
  print *, '=========================================================================='
  print *, '[Observation site and Data Level]'
  if (site_name == 'SIGMA-A') print *, '>>> Observation site: SIGMA-A'
  if (site_name == 'SIGMA-B') print *, '>>> Observation site: SIGMA-B' 
  print *, '>>> Making Level 1.3: Secondary Control (Completed Quality Control)'

  if (site_name == 'SIGMA-A') outfile1 = site_name // '_' // cyear // '_met_Lv1.3.csv'
  if (site_name == 'SIGMA-B') outfile1 = site_name // '_' // cyear // '_met_Lv1.3.csv'

  call data_number (infile, ndata)

  allocate(character(50) :: cell(ndata+1,param_num))
  allocate(cdate_time(ndata), ctime(ndata))
  allocate(year(ndata), month(ndata), day(ndata))
  allocate(atemp1(ndata), arh1(ndata), ws1(ndata), wd1(ndata))
  allocate(atemp2(ndata), arh2(ndata), ws2(ndata), wd2(ndata))
  allocate(swd(ndata), swu(ndata), lwd(ndata), lwu(ndata), lw_standard(ndata), nird(ndata), niru(ndata))
  allocate(press(ndata), snow_depth(ndata), sensor_height1(ndata), sensor_height2(ndata))
  allocate(sst(ndata), st1(ndata), st2(ndata), st3(ndata), st4(ndata), st5(ndata), st6(ndata))
  allocate(solz(ndata), sola(ndata))
  allocate(albedo(ndata), albedo_acc(ndata), nir_albedo(ndata), nir_albedo_acc(ndata), nir_frac(ndata))
  allocate(SW_TOA(ndata), nir_max(ndata))
  allocate(swd_slope(ndata), solz_slope(ndata), Id_slope(ndata), Is_slope(ndata), diffu_ratio(ndata))
  allocate(st1_depth(ndata), st2_depth(ndata), st3_depth(ndata), st4_depth(ndata), st5_depth(ndata), st6_depth(ndata))
  allocate(st1_mean(ndata), st2_mean(ndata), st3_mean(ndata), st4_mean(ndata), st5_mean(ndata), st6_mean(ndata))
  allocate(sd_mean(ndata), press_mean(ndata), delta_press(ndata))

  !>>> reading data file
  if (site_name == 'SIGMA-A') call data_inputA (infile, cdate_time, &
                                  & ws1, wd1, atemp1, arh1, ws2, &
                                  & wd2, atemp2, arh2, swd, swu, &
                                  & lwd, lwu, nird, niru, snow_depth, &
                                  & press, st1, st2, st3, st4, &
                                  & st5, st6, solz, sola)
  
  if (site_name == 'SIGMA-B') call data_inputB (infile, cdate_time, &
                                  & ws1, wd1, atemp1, arh1, swd, &
                                  & swu, lwd, lwu, snow_depth, press, &
                                  & solz, sola)  
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
    
    time = trim(ctime(i)) // ':00'
    if (time(2:2) == ':') then
      time = '0' // trim(time)
    end if
  end do

  !----------------------------------------------------------------------------
  !================================
  ! Secondary QC
  !================================
  do i = 1, ndata
    call SW_TOA_calculation (solz(i), year(i), month(i), day(i), SW_TOA(i))
  end do

  call QC_wind_speed (ndata, wd1, wd2, ws1, ws2)
  call QC_air_temperature (ndata, site_name, atemp1, atemp2)
  call QC_air_humidity (ndata, site_name, arh1, arh2)
  call QC_shortwave_radiation (ndata, site_name, SW_TOA, solz, solz_slope, sola, swd, swd_slope, swu)
  call QC_nir_infrared_radiation (ndata, site_name, SW_TOA, solz, solz_slope, nird, niru)
  call QC_albedo (ndata, site_name, month, swd, swd_slope, swu, SW_TOA, solz, solz_slope, lwd, albedo)
  call QC_nir_infrared_albedo (ndata, site_name, month, nird, niru, solz, lwd, nir_albedo)
  call nir_infrared_fraction (ndata, site_name, swd, nird, solz, nir_frac)
  call QC_longwave_radiation (ndata, site_name, atemp1, atemp2, lwd, lwu)
  call QC_snow_depth (ndata, site_name, snow_depth)
  call sensor_height (ndata, site_name, snow_depth, sensor_height1, sensor_height2)
  call snow_temperature_sensor_depth (ndata, site_name, snow_depth, st1_depth, st2_depth, st3_depth, st4_depth, st5_depth, st6_depth)
  call QC_snow_temperature (ndata, site_name, st1_depth, st2_depth, st3_depth, st4_depth, st5_depth, st6_depth, st1, st2, st3, st4, st5, st6, st1_mean, st2_mean, st3_mean, st4_mean, st5_mean, st6_mean)
  call QC_air_press (ndata, site_name, press)

  do i = 1, ndata
    sst(i) = snow_surface_temperature (lwd(i), lwu(i))
  end do  

  ! accumulated albedo
  if (site_name == 'SIGMA-A') then
    call accumulated_albedo (ndata, swd, swu, albedo_acc)
    call accumulated_albedo (ndata, nird, niru, nir_albedo_acc)
  end if
  if (site_name == 'SIGMA-B') then
    call accumulated_albedo (ndata, swd_slope, swu, albedo_acc)
    call accumulated_albedo (ndata, nird, niru, nir_albedo_acc)
  end if
  call QC_albedo_acc (ndata, site_name, month, albedo_acc, nir_albedo_acc)
  call QC_albedo_acc (ndata, site_name, month, albedo_acc, nir_albedo_acc)

  print *, '[Main Program:] All QC procedures... Done'

  !----------------------------------------------------------------------------
  !>>> Level 1.3 data set output
  call data_writing (site_name, ndata, cell, cdate_time, &
                   & ws1, wd1, atemp1, arh1, ws2, &
                   & wd2, atemp2, arh2, swd, swu, &
                   & lwd, lwu, nird, niru, snow_depth, &
                   & press, st1, st2, st3, st4, &
                   & st5, st6, solz, sola, albedo, &
                   & albedo_acc, nir_albedo, nir_albedo_acc, nir_frac, sensor_height1, &
                   & sensor_height2, sst, st1_depth, st2_depth, st3_depth, &
                   & st4_depth, st5_depth, st6_depth, solz_slope, swd_slope, column_number)

  call data_output (ndata, column_number, outfile1, cell)
  !----------------------------------------------------------------------------

  deallocate(cell)
  deallocate(cdate_time, ctime)
  deallocate(year, month, day)  
  deallocate(atemp1, arh1, ws1, wd1, atemp2, arh2, ws2, wd2)
  deallocate(swd, swu, lwd, lwu, lw_standard, nird, niru)
  deallocate(press, snow_depth, sensor_height1, sensor_height2)
  deallocate(sst, st1, st2, st3, st4, st5, st6)
  deallocate(solz, sola)
  deallocate(albedo, albedo_acc, nir_albedo, nir_albedo_acc, nir_frac)
  deallocate(SW_TOA, nir_max, swd_slope, Id_slope, Is_slope, solz_slope, diffu_ratio)
  deallocate(st1_depth, st2_depth, st3_depth, st4_depth, st5_depth, st6_depth)
  deallocate(st1_mean, st2_mean, st3_mean, st4_mean, st5_mean, st6_mean)
  deallocate(sd_mean, press_mean, delta_press)

  print *, '[Main Program] >>> Program end ...'
  print *, '==============================================================================='

end program makedata
