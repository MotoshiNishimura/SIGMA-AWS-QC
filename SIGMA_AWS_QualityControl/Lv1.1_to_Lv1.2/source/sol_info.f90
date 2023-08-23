module sol_info

  private
  public :: solar
  
contains
  !----------------------------------------------------------------------------
  subroutine solar(ymd, time, lat, lon, solz, sola)
    !===============================================================
    !
    ! read (year-month-day, time, lattitude, longitude), and 
    ! return solar zenith angle and solar azimuthal angle
    ! in character.
    !
    ! *** input like ***
    ! ymd: 20081212
    ! time: 12:12:12
    ! lat: 36:03:10
    ! lon: -140:-07:-35
    ! ******************
    !
    ! this is based on solz.30min.f (made by Dr.Aoki)
    ! 
    ! Masashi Niwano @ Meteorological Research Institute
    ! 
    ! ver0.1 05/23/2010
    !
    !===============================================================
    
    ! parameters
    implicit none
    character(len=100) :: ymd, time, lat, lon ! in
    character(len=4) :: year, month, day, hour, minu, sec
    character(len=4) :: lat1, lat2, lat3, lon1, lon2, lon3 
    integer :: iyear, imonth, iday, i = 0, j, n
    integer :: ihour, iminu, isec
    real(8) :: rlat1, rlat2, rlat3, rlon1, rlon2, rlon3
    real(8) :: rtime, rlat, rlon
    real(8) :: alf, dlt
    real(8) :: a1(19), a2(19), a3(19), pi, rd, dg 
    integer :: nday(12)
    integer :: nyy, nmm, ndd, nhh, nmi, nss, ndayl
    real(8) :: w, f, utj, t, tt, gst, al, p
    real(8) :: ah, solz, wk1, hs, sola, wk2
    integer :: nw, nx, ny, nr, ns, nz
    
    ! for Japan !
    ! integer, parameter :: ndhh = 9, ndmi = 0
    ! for GrIS
    integer, parameter :: ndhh = 0, ndmi = 0 ! input data is in GMT    

    data (nday(j), j = 1, 12) / 31, 28, 31, 30, 31, 30, &
         &                      31, 31, 30, 31, 30, 31 /
    
    data (a1(j), j=1,19) / 1.9159D0, 0.0200D0, 0.0048D0, 0.0020D0, 0.0018D0, &
         &                 0.0018D0, 0.0015D0, 0.0013D0, 0.0008D0, 0.0007D0, &
         &                 0.0007D0, 0.0006D0, 0.0005D0, 0.0005D0, 0.0004D0, &
         &                 0.0004D0, 0.0004D0, 0.0003D0, 0.0003D0 /
    
    data (a2(j), j=1,19) / 356.531D0, 353.06D0, 68.64D0, 285.0D0, 334.2D0, &
         &                   293.7D0,  242.4D0, 211.1D0, 208.0D0,  53.5D0, &
         &                    12.1D0,  239.1D0,  10.1D0,  99.1D0, 264.8D0, &
         &                   233.8D0,   18.1D0, 349.6D0, 241.2D0 /
    
    data (a3(j),j=1,19) / 359.991D0,719.981D0,-19.341D0,329.64D0,-4452.67D0, &
         &                  -0.20D0, 450.37D0, 225.18D0,659.29D0,   90.38D0, &
         &                 -30.35D0, 337.18D0,  -1.50D0,-22.81D0,  315.56D0, &
         &                 299.30D0,  720.02D0, 1079.97D0, -44.43D0 / 
    
    
    pi = acos(-1.0D0)
    rd = pi / 180.0D0
    dg = 1.0D0 / rd
    
    ! write(*,*) 'pi = ', pi
    
    !===============================================================
    ! read ymd    
    !  write(*,*) 'Input a day like 20080101'
    !  read(*,*) ymd
    
    year = ymd(1:4)
    month = ymd(5:6)
    day = ymd(7:8)
    
    read(year,*) iyear
    
    i = 0
    j = len(trim(month))
    do n = 1, j
       i = i * 10 + ichar(month(n:n))-48
       imonth = i
    end do
    
    i = 0
    j = len(trim(day))
    do n = 1, j
       i = i * 10 + ichar(day(n:n))-48
       iday = i
    end do
    
    ! write(*,*) 'ymd = ', ymd
    ! write(*,*) 'year = ', year, iyear
    ! write(*,*) 'month = ', month, imonth
    ! write(*,*) 'day = ', day, iday 
    
    !===============================================================
    ! read time
    !  write(*,*) 'Input a time like 12:12:12'
    !  read(*,*) time
    
    hour = time(1:2)
    minu = time(4:5)
    sec = time(7:8)
    
    i = 0
    j = len(trim(hour))
    do n = 1, j
       i = i * 10 + ichar(hour(n:n))-48
       ihour = i
    end do
    
    i = 0
    j = len(trim(minu))
    do n = 1, j
       i = i * 10 + ichar(minu(n:n))-48
       iminu = i
    end do

    i = 0
    j = len(trim(sec))
    do n = 1, j
       i = i * 10 + ichar(sec(n:n))-48
       isec = i
    end do
    
    rtime = ihour + ( iminu + isec / 60.0_8 ) / 60.0_8
    
    ! write(*,*) 'hour = ', hour, ihour
    ! write(*,*) 'min = ', minu, iminu
    ! write(*,*) 'sec = ', sec, isec
    ! write(*,*) 'rtime =', rtime
  
    !===============================================================
    ! conversion of ymd and time to world time
    if(mod(iyear, 4) .EQ. 0) then
       nday(2) = 29
    end if
    
    if(mod(iyear, 100) .EQ. 0) then
       nday(2) = 28
    end if
    
    if(mod(iyear, 400) .EQ. 0) then
       nday(2) = 29
    end if
    
    nyy = iyear
    nmm = imonth
    ndd = iday
    nhh = ihour - ndhh
    nmi = iminu - ndmi
    nss = isec
    
    if(nmi .LT. 0) then
       nmi = nmi + 60
       nhh = nhh - 1
    end if
    
    if(nhh .LT. 0) then
       nhh = nhh + 24
       ndd = ndd - 1
    end if
    
    if(ndd .EQ. 0) then
       if((nmm - 1) .EQ. 0) then
          ndayl = 31
       else 
          ndayl = nday(nmm - 1)
       end if
       ndd = ndayl
       nmm = nmm - 1
    end if
    
    if(nmm .EQ. 0) then
       nmm = 12
       nyy = nyy - 1
    end if
    
    if((mod(nyy,4) .EQ. 0) .AND. (mod(nyy,100) .NE. 0)) then
       nday(2) = 29
    end if
    
    !===============================================================
    ! sequencial day(NZ) from Jan.1.1975
    w = (nyy - 1900D0) / 4.0D0
    nw = int(w)
    f = w - nw
    ny = int(1461.0D0 * w)
    nx = int((nmm + 7.0D0) / 10.0D0)
    nr = int(1.0D0 - f)
    ns = int(0.44D0 * (nmm + 4.4D0))
    nz = ny + 31 * nmm + ndd + (nx - 1) * nr - nx * ns - 27424
    
    ! write(*,*) 'nz = ', nz
    
    !===============================================================
    ! read lat&lon
    !  write(*,*) 'Input lat like 36:03:10'
    !  read(*,*) lat
    
    lat1 = lat(1:2)
    lat2 = lat(4:5)
    lat3 = lat(7:8)
    
    i = 0
    j = len(trim(lat1))
    do n = 1, j
       i = i * 10 + ichar(lat1(n:n))-48
       rlat1 = i
    end do
    
    i = 0
    j = len(trim(lat2))
    do n = 1, j
       i = i * 10 + ichar(lat2(n:n))-48
       rlat2 = i
    end do
    
    i = 0
    j = len(trim(lat3))
    do n = 1, j
       i = i * 10 + ichar(lat3(n:n))-48
       rlat3 = i
    end do
    
    rlat = rlat1 + ( rlat2 + rlat3 / 60.0_8 ) / 60.0_8
    
    !  write(*,*) 'Input lon like -140:-07:-36'
    !  read(*,*) lon
    
    lon1 = lon(2:4)
    lon2 = lon(7:8)
    lon3 = lon(11:12)
    
    i = 0
    j = len(trim(lon1))
    do n = 1, j
       i = i * 10 + ichar(lon1(n:n))-48
       rlon1 = i * (-1)
    end do
    
    i = 0
    j = len(trim(lon2))
    do n = 1, j
       i = i * 10 + ichar(lon2(n:n))-48
       rlon2 = i * (-1)
    end do
    
    i = 0
    j = len(trim(lon3))
    do n = 1, j
       i = i * 10 + ichar(lon3(n:n))-48
       rlon3 = i * (-1)
    end do
        
    ! for GrIS
    rlon = rlon1 + ( rlon2 + rlon3 / 60.0_8 ) / 60.0_8

    ! write(*,*) 'lat1 = ', lat1, rlat1
    ! write(*,*) 'lat2 = ', lat2, rlat2
    ! write(*,*) 'lat3 = ', lat3, rlat3
    ! write(*,*) 'rlat = ', rlat
    ! write(*,*) 'lon1 = ', lon1, rlon1
    ! write(*,*) 'lon2 = ', lon2, rlon2
    ! write(*,*) 'lon3 = ', lon3, rlon3
    ! write(*,*) 'rlon = ', rlon
    
    !===============================================================
    ! main
    ! ut expressed as the mod of day
    utj = nhh / 24.0D0 + nmi / 1440.0D0 + nss /86400.0D0
    
    ! time parameter
    t = (nz + utj) / 365.25D0
    tt = t + (0.0286D0 * t + 1.407D0) * 1.0D-6 
    
    ! Greenwich day
    gst = 99.0361D0 + 360.00770D0 * t + 360.0D0 * utj &
         & + 0.0044D0 * sin((68.6D0 - 19.3D0 * t) * rd) &
         & + 0.0003D0 * sin((18.0D0 + 720D0 * t) * rd)
    
    ! koukei of sun
    al = 279.0415D0 + 360.00769D0 * tt
    al = al + (a1(1) - 0.00005D0 * tt ) * sin((a2(1) + a3(1) * tt) * rd)
  
    do j = 2, 19
       al = al + a1(j) * sin((a2(j) + a3(j) * tt) * rd)
    end do
    
    al = al - 0.0057D0
    
    ! koukeikeisya
    p = 23.44254D0 - 0.00013D0 * tt &
         & + 0.00256D0 * cos((248.6D0 - 19.3D0 * tt) * rd) &
         & + 0.00016D0 * cos((198.0D0 + 720.0D0 * tt) * rd)
    
    ! sekkei, sekii
    alf = atan(tan(al * rd) * cos(p * rd)) * dg
    
    if (cos(al * rd) .LT. 0.0) then
       alf = alf + 180.0D0
    end if
    
    dlt = asin(sin(al * rd) * sin(p * rd)) * dg
    
    ! jikaku, solar zeith, and azimuth
    ah = gst - alf - rlon
    solz = acos(sin(dlt * rd) * sin(rlat * rd) &
         & + cos(dlt * rd) * cos(rlat * rd) * cos(ah * rd))
    solz = solz * dg
    
    wk1 = tan(dlt * rd) * cos(rlat * rd) -cos(ah * rd) * sin(rlat * rd)
    if (wk1 .EQ. 0.0) then
       hs = ah - int(ah / 360.0D0) * 360D0
       if(hs .GT. 180.0D0) then
          sola = 90.0D0
       else 
          sola = 270.0D0
       end if
    else 
       wk2 = -sin(ah * rd)
       sola = atan(wk2 / wk1)
       sola = sola * dg
       if(wk1 .LT. 0.0D0) then
          sola = sola + 180.0D0
       end if
       if((wk1 .GT. 0.0D0) .AND. (wk2 .LT. 0.0D0)) then
          sola = sola + 360.0D0
       end if
    end if
    
    !===============================================================
    
    ! write(*,*) 'solar azimuthal angle = ', sola
    
    ! write(csolzen, '(f7.2)') solz
    ! write(csolazi, '(f7.2)') sola
    
    ! write(*,*) 'csolzen = ', csolzen
    ! write(*,*) 'csolazi = ', csolazi
    
    return
  end subroutine solar
  
end module sol_info

