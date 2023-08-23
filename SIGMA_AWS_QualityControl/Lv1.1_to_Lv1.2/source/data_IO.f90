!=======================================================================
! This module summarizes the data input and output process for the Quality Control (Initial Control)
! Motoshi Nishimura @ National Institute of Polar Research
!=======================================================================

module data_input_output
  implicit none
  contains

subroutine data_number (infile, ndata)
  implicit none
  
  character(*), intent(in) :: infile
  integer(8), intent(out)  :: ndata

  integer(8) :: ios

  print *, ''
  print *, '--------------------------------------------------------------' 
  print *, '[subroutine: data_number] Scanning Data process begins.' 
  open(10, iostat=ios, file=infile, status='old', form='formatted')
    if (ios /= 0) then
      print *, 'I/O status  = ', ios
      print *, '================================================='
      print *, 'Error1-1!! Failed to open input_file: ', infile
      print *, '================================================='
      stop
    end if

    ndata = 0
    do
      read(10, fmt=*, iostat=ios)
      if (ios > 0) stop 'Error'
      if (ios < 0) exit
      ndata = ndata + 1
    end do
    
  close(10)
  ndata = ndata - 1

  print *, 'A number of data is', ndata 
  print *, '[subroutine: data_number] Scanning Data process... Done' 
  
end subroutine data_number

!---------------------------------------------------------------------------------------------------------------------------  
subroutine data_input (infile, date_time, &
                     & param1, param2, param3, param4, param5, param6, param7, param8, param9, param10, &
                     & param11, param12, param13, param14, param15, param16, param17, param18, param19, param20, &
                     & param21, param22, param23, param24, param25, param26, param27, param28, param29)
  implicit none
  
  character(*), intent(in)  :: infile
  character(*), intent(out) :: date_time(:)
  real(8), intent(out)      :: param1(:), param2(:), param3(:), param4(:), param5(:)
  real(8), intent(out)      :: param6(:), param7(:), param8(:), param9(:), param10(:)
  real(8), intent(out)      :: param11(:), param12(:), param13(:), param14(:), param15(:)
  real(8), intent(out)      :: param16(:), param17(:), param18(:), param19(:), param20(:)
  real(8), intent(out)      :: param21(:), param22(:), param23(:), param24(:), param25(:)
  real(8), intent(out)      :: param26(:), param27(:), param28(:), param29(:)
  
  integer(8) :: ios
  integer(8) :: i

  print *, '--------------------------------------------------------------' 
  print *, '[subroutine: data_input] Data input process begins.' 
   
  open (10, iostat=ios, file=infile, status='old', form='formatted')
    if (ios /= 0) then
      write (*,*) 'I/O status = ', ios
      write (*,*) 'Error2-1!! Failed to open input_file: ', infile
      stop
    end if
    read (10,'()')  !'()' -> omit header line

    do i = 1, size(date_time)
      read(10,*) date_time(i), &
               & param1(i), param2(i), param3(i), param4(i), param5(i), &
               & param6(i), param7(i), param8(i), param9(i), param10(i), &
               & param11(i), param12(i), param13(i), param14(i), param15(i), &
               & param16(i), param17(i), param18(i), param19(i), param20(i), &
               & param21(i), param22(i), param23(i), param24(i), param25(i), &
               & param26(i), param27(i), param28(i), param29(i)
      if (i == 1) then
        print *, i, date_time(i), &
                  & param1(i), param2(i), param3(i), param4(i), param5(i), &
                  & param6(i), param7(i), param8(i), param9(i), param10(i), &
                  & param11(i), param12(i), param13(i), param14(i), param15(i), &
                  & param16(i), param17(i), param18(i), param19(i), param20(i), &
                  & param21(i), param22(i), param23(i), param24(i), param25(i), &
                  & param26(i), param27(i), param28(i), param29(i)
      end if
    end do
    
    print *, '================================================================'
    print *, 'Start date/time: ', date_time(1)
    print *, 'End   date/time: ', date_time(size(date_time))
    print *, '================================================================'
  close(10)
  
  print *, '[subroutine: data_input] Data input process... Done' 

end subroutine data_input

!---------------------------------------------------------------------------------------------------------------------------
subroutine data_writing (site_name, ndata, cell, date_time, &
                       & param1, param2, param3, param4, param5, param6, param7, param8, param9, param10, &
                       & param11, param12, param13, param14, param15, param16, param17, param18, param19, param20, &
                       & param21, param22, param23, param24, column_number)
  implicit none
  
  character(*), intent(inout) :: cell (:,:)
  character(*), intent(in)    :: site_name
  character(*), intent(in)    :: date_time(:)
  integer(8), intent(in)      :: ndata
  real(8), intent(in)         :: param1(:), param2(:), param3(:), param4(:), param5(:)
  real(8), intent(in)         :: param6(:), param7(:), param8(:), param9(:), param10(:)
  real(8), intent(in)         :: param11(:), param12(:), param13(:), param14(:), param15(:)
  real(8), intent(in)         :: param16(:), param17(:), param18(:), param19(:), param20(:)
  real(8), intent(in)         :: param21(:), param22(:), param23(:), param24(:)
  integer(8), intent(out)     :: column_number
    
  integer(8) :: i

  print *, ''
  print *, '--------------------------------------------------------------' 
  print *, '[subroutine: data_writing] Data writing process begins.' 
  
  do i = 1, size(date_time)
    cell(i+1,1) = date_time(i)
  end do

  !>>> appending parameters
  cell(1,1) = 'date'
	! column number: initial column number of the dataset
  column_number = 2
  if (site_name == 'SIGMA-A') then
    call column (cell, column_number, 'U1', ndata, param1)
    call column (cell, column_number, 'WD1', ndata, param2)
    call column (cell, column_number, 'T1', ndata, param3)
    call column (cell, column_number, 'RH1', ndata, param4)
    call column (cell, column_number, 'U2', ndata, param5)
    call column (cell, column_number, 'WD2', ndata, param6)
    call column (cell, column_number, 'T2', ndata, param7)
    call column (cell, column_number, 'RH2', ndata, param8)
    call column (cell, column_number, 'SWd', ndata, param9)
    call column (cell, column_number, 'SWu', ndata, param10)
    call column (cell, column_number, 'LWd', ndata, param11)
    call column (cell, column_number, 'LWu', ndata, param12)
    call column (cell, column_number, 'NIRd', ndata, param13)
    call column (cell, column_number, 'NIRu', ndata, param14)
    call column (cell, column_number, 'sh', ndata, param15)
    call column (cell, column_number, 'Pa', ndata, param16)
    call column (cell, column_number, 'st1', ndata, param17)
    call column (cell, column_number, 'st2', ndata, param18)
    call column (cell, column_number, 'st3', ndata, param19)
    call column (cell, column_number, 'st4', ndata, param20)
    call column (cell, column_number, 'st5', ndata, param21)
    call column (cell, column_number, 'st6', ndata, param22)
    call column (cell, column_number, 'solz', ndata, param23)
    call column (cell, column_number, 'sola', ndata, param24)
  end if
  
  if (site_name == 'SIGMA-B') then
    call column (cell, column_number, 'U1', ndata, param1)
    call column (cell, column_number, 'WD1', ndata, param2)
    call column (cell, column_number, 'T1', ndata, param3)
    call column (cell, column_number, 'RH1', ndata, param4)
    call column (cell, column_number, 'SWd', ndata, param9)
    call column (cell, column_number, 'SWu', ndata, param10)
    call column (cell, column_number, 'LWd', ndata, param11)
    call column (cell, column_number, 'LWu', ndata, param12)
    call column (cell, column_number, 'sh', ndata, param15)
    call column (cell, column_number, 'Pa', ndata, param16)
    call column (cell, column_number, 'solz', ndata, param23)
    call column (cell, column_number, 'sola', ndata, param24)
  end if
  
  print *, '[subroutine: data_writing] Data writing process... Done'

end subroutine data_writing

!---------------------------------------------------------------------------------------------------------------------------
subroutine column (cell, n, param_name, ndata, param)
  implicit none

  character(*), intent(inout) :: cell(:,:)
  integer(8), intent(inout)   :: n
  character(*), intent(in)    :: param_name
  integer(8), intent(in)      :: ndata
  real(8), intent(in)         :: param(:)

  integer(8) :: i

  ! n: column_number, i: line number
  cell(1, n) = param_name

  do i = 1, ndata
    write (cell(i+1,n), '(f12.5)') param(i)
  end do
  n = n + 1

end subroutine column

!---------------------------------------------------------------------------------------------------------------------------
subroutine data_output (ndata, column_number, outfile1, cell)
  implicit none
  
  character(*), intent(in)    :: outfile1
  character(*), intent(inout) :: cell (:,:)
  integer(8), intent(in)      :: ndata
  integer(8), intent(in)      :: column_number
  
  integer(8) :: ios
  integer(8) :: i, j

  print *, '--------------------------------------------------------------'   
  print *, '[subroutine: data_output] Data output process begins.' 
  open (20, iostat=ios, file = outfile1, status = 'replace', form = 'formatted')

  if (ios /= 0) then
    print *, 'I/O status = ', ios
    print *, 'Error20!! Failed to open output_file: ', outfile1
    stop
  end if
  
  print *, 'a number of column: ', column_number
  if (column_number > 100) then
    print *, 'ERROR20!! A number of columns exceeds the number of file export columns.'
    print *, 'Please check the subroutine data output'
    stop
  end if 

  do i = 1, ndata+1
    write (20, 2000) (trim(cell(i,j)), j = 1, 100)
    2000 format (A, 99(',', A))
  end do

  close (20)
  print *, '[subroutine: data_output] Data output process... Done' 
  print *, '--------------------------------------------------------------' 


end subroutine data_output

end module data_input_output
