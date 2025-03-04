program LocalGridExtremes
!
! Compile with:
!
! module load compiler/gcc-7.3 development/netcdf
! gfortran $(nf-config --fflags) $(nf-config --flibs) XXX.f90 -o XXX.Abs
use netcdf
implicit none

integer ncid0 ,ncid1,ncid2, ncid3, ncid4
integer rcode, dimidix, dimidjx, dimidtx, varid, ix, ixs, jx, jxs, tx
character(len=3) :: dimnameix
character(len=3) :: dimnamejx
character(len=4) :: dimnametX
character(len=9) :: dimnameix2
character(len=8) :: dimnamejx2
real dist12,ppoints4, ppoints6, ppoints8, ppoints10
integer i, j, t, time_index, ys, ye, ms, me, n, p, s, tt, oo 
integer nummax, yyyy_int,mm_int,dd_int,hh_int
integer line, linemax, lineheader, year_int, month_int
integer tacc, totaltx
parameter (linemax = 934010, lineheader = 16, nummax=500)
real TPpercentile1, TPpercentile2, TPpercentile3, WDpercentile

character*130 fileout0, fileout1, fileout2, fileout3, fileout4, &
              fileout5, fileout6, &
              fileout1a, fileout2a, fileout3a 
character*130 ifile0, ifile1, ifile2, ifile3, ifile4
character*2 mm, dd, hh, month
character*4 year
character*4 TPpercentile_str1, TPpercentile_str2, TPpercentile_str3, WDpercentile_str

integer iun, iun2, iun3, iun4, l_nj, l_ni, ftyp, ezgdef_fmem, gdll
real scale_factor1, add_offset1, missing_value1
real scale_factor2, add_offset2, missing_value2
real scale_factor3, add_offset3, missing_value3
real dx
logical global, ll_double

! Internal parameters
real, parameter :: earthr = 6371.22 ! in kilometers
real            :: deg2rad, glons0conv


real, parameter :: R_threshold  =  1000.    ! in kilometers (1100km used in Owen et al., 2021)
real            :: distc
real, parameter :: Di_threshold  =  55.     ! in degree (set looser)
real, parameter :: Dj_threshold  =  15.     ! in degree (set looser)
logical :: ETCpresence
! Threshold data
real,       dimension (:,:)  , allocatable :: TPper1, TPper2, TPper3
real,       dimension (:)    , allocatable :: glats0, glons0
real,       dimension (:,:,:), allocatable :: tpo,tp

integer TPpercent_nx1, TPpercent_nx2, TPpercent_nx3, WDpercent_nx
integer,    dimension (:,:)  , allocatable :: TPorder1, TPorder2, TPorder3
integer,    dimension (:)  , allocatable :: num

! Storm data
integer jday0, jday, ic, icmax, icmin, jc, jcmax, jcmin
integer ixwest, jxnorth
real latc0, lonc0
real,    dimension(:,:), allocatable :: latc, lonc 

! Output data
real,       dimension (:,:,:), allocatable :: TP_ex1, TP_ex2, TP_ex3
integer,    dimension (:,:,:), allocatable :: TP_extime1, TP_extime2, TP_extime3, &
                                              ETC_ex1, ETC_ex2, ETC_ex3
!=========
deg2rad = acos(-1d0)/180d0
!=========

! Period of interest 
ys=2001
ye=2020
ms=1
me=12

call date_to_julian(ys,ms,1,jday0)

! Total time steps (hourly data) during the entire period
! and the corresponding sample numbers for the selected threshold
totaltx     =  7305*24 ! from 01/01/2001 to 31/12/2020 (included)/

print*,'totaltx:', totaltx
! -----------------------------------------

! Read in the tracked storm data

allocate (num(totaltx))
allocate (latc(totaltx,nummax))
allocate (lonc(totaltx,nummax))
num = 0
latc = -9999.
lonc = -9999.

open (10, FILE='ETC_identified_n_tracking_output_2000_2020.txt', STATUS='old')

  do line = 1, linemax

    if (line .le. lineheader) then
        read(10,*)     !header
    else
        read(10,100)  yyyy_int, mm_int, dd_int, hh_int, latc0, lonc0

        if (yyyy_int .ge. ys) then
        
            call date_to_julian(yyyy_int,mm_int,dd_int,jday)
            tacc = 1+(jday-jday0)*24+hh_int   !tacc=1 indicates 1979 Jan 01 00UTC

            num(tacc) = num(tacc)+1
            latc(tacc, num(tacc)) = latc0
            lonc(tacc, num(tacc)) = lonc0

        endif
    endif

  enddo

close(10)

100 format (24x, i4, i2, i2, i2, 3x, f5.2, 2x, f6.2)


!===Setup parameters===

  TPpercentile1 = 98
  TPpercentile2 = 99
  TPpercentile3 = 99.9

  TPpercentile_str1 = '98p0'
  TPpercentile_str2 = '99p0'
  TPpercentile_str3 = '99p9'

  jxs = 201
  ixs = 521
  ! 0.25 degree per grid length, we want starting from 240: ixwest=961
  ! 230: ixwest=921
  ixwest  = 841 ! The western boundary of the targeted region in the ERA5 domain
  ! 0.25 degree per grid length, we want starting from 80 N (j start:41) to 25N: gxnorth=41
  ! 75 N: jxnorth= 61
  jxnorth = 61  ! The northern boundary of the targeted region in the ERA5 domain


!  WDpercent_nx  = 1+int((100.-WDpercentile)/100*float(totaltx))
  TPpercent_nx1 = int((100.-TPpercentile1)/100*float(totaltx))
  TPpercent_nx2 = int((100.-TPpercentile2)/100*float(totaltx))
  TPpercent_nx3 = int((100.-TPpercentile3)/100*float(totaltx))

  print*,'TPpercent_num1:', TPpercent_nx1
  print*,'TPpercent_num2:', TPpercent_nx2
  print*,'TPpercent_num3:', TPpercent_nx3
!  print*,'WDpercent_num:', WDpercent_nx

!======================

! input data: Percentile values in TP and WDSP
ifile1 = 'TPp98p0_2001_2020_ERA5.nc'
ifile2 = 'TPp99p0_2001_2020_ERA5.nc'
ifile3 = 'TPp99p9_2001_2020_ERA5.nc'
!
!------------------------------------
! Read in the threshold value at each grid point
  allocate (TPper1(ixs, jxs))
  allocate (TPper2(ixs, jxs))
  allocate (TPper3(ixs, jxs))

! Open input files
  rcode = nf90_open(ifile1, nf90_nowrite, ncid1)
  if (rcode /= nf90_noerr) call handle_err(rcode)
  rcode = nf90_open(ifile2, nf90_nowrite, ncid2)
  if (rcode /= nf90_noerr) call handle_err(rcode)
  rcode = nf90_open(ifile3, nf90_nowrite, ncid3)
  if (rcode /= nf90_noerr) call handle_err(rcode)

  !TP:
  call check(nf90_inq_varid(ncid1,'p98p0',varid))
  call check(nf90_get_var(ncid1,varid,TPper1))

  call check(nf90_inq_varid(ncid2,'p99p0',varid))
  call check(nf90_get_var(ncid2,varid,TPper2))

  call check(nf90_inq_varid(ncid3,'p99p9',varid))
  call check(nf90_get_var(ncid3,varid,TPper3))

  call check(nf90_close(ncid1))
  call check(nf90_close(ncid2))
  call check(nf90_close(ncid3))


!--------------------------------------------------------------------------
! Read in the source data files 

tacc = 0 
do year_int  = ys, ye
do month_int = ms, me

    print*, 'year=',year_int,'month=',month_int

    write(year,  '(i4)') year_int
    write(month, '(i2.2)') month_int

    ifile1 = '/home/archive/REANALYSES/ERA5/1h/tp/ll/nc4/'//year//'/'//month//'/era5_tp_ll_'//year//''//month//'_1h.nc4'

    ! Open input file
    rcode = nf90_open(ifile1, nf90_nowrite, ncid1)
    if (rcode /= nf90_noerr) call handle_err(rcode)

    call check(nf90_inq_dimid(ncid1,'time',dimidtx))
    call check(nf90_inquire_dimension(ncid1,dimidtx,dimnametx,tx))

    ! [Below: only declare once]
    if (year_int.eq.ys .and. month_int.eq.ms) then

      allocate (glons0(ixs)) 
      allocate (glats0(jxs)) 
      call check(nf90_inq_varid(ncid1,'longitude',varid))
      call check(nf90_get_var(ncid1,varid,glons0, start=(/ixwest/), count=(/ixs/)))
      call check(nf90_inq_varid(ncid1,'latitude',varid))
      call check(nf90_get_var(ncid1,varid,glats0, start=(/jxnorth/) , count=(/jxs/))) 
      print*, 'The starting grid point in the geographic domain in the IMERG data'
      print*, 'glats0(1)=',glats0(1),'glats0(jxs)',glats0(jxs)
      print*, 'glons0(1)=',glons0(1),'glons0(ixs)',glons0(ixs)
      !
      ! [Output variables]
      ! Extreme Values
      allocate (TP_ex1(ixs,jxs,TPpercent_nx1))
      allocate (TP_ex2(ixs,jxs,TPpercent_nx2))
      allocate (TP_ex3(ixs,jxs,TPpercent_nx3))

      ! Extreme Timesteps (remeber to indicate the starting date in the output!)
      allocate (TP_extime1(ixs,jxs,TPpercent_nx1))
      allocate (TP_extime2(ixs,jxs,TPpercent_nx2))
      allocate (TP_extime3(ixs,jxs,TPpercent_nx3))
     
      ! Number/Index of extremes  
      allocate (TPorder1(ixs,jxs))
      allocate (TPorder2(ixs,jxs))
      allocate (TPorder3(ixs,jxs))

      ! ETC association (1 indicates the extreme is associated with one or more ETCs in 1000 km; 0 not)
      allocate (ETC_ex1(ixs,jxs,TPpercent_nx1))
      allocate (ETC_ex2(ixs,jxs,TPpercent_nx2))
      allocate (ETC_ex3(ixs,jxs,TPpercent_nx3))

      ETC_ex1=0
      ETC_ex2=0
      ETC_ex3=0
      TP_extime1=0
      TP_extime2=0
      TP_extime3=0
      TP_ex1=missing_value1
      TP_ex2=missing_value1
      TP_ex3=missing_value1

    endif

    allocate (tpo(ixs,jxs,tx))
    allocate (tp(ixs,jxs,tx))
    
    !TP:
    call check(nf90_inq_varid(ncid1,'tp',varid))
    call check(nf90_get_var(ncid1,varid,tpo, start=(/ixwest,jxnorth,1/), count=(/ixs,jxs,tx/)))
    call check(nf90_get_att(ncid1,varid,'scale_factor',scale_factor1))
    call check(nf90_get_att(ncid1,varid,'add_offset',add_offset1))
    call check(nf90_get_att(ncid1,varid,'missing_value',missing_value1))
    tp=(tpo*scale_factor1+add_offset1)*1000.

    deallocate(tpo)

    call check(nf90_close(ncid1))

    !==============

 
    do t = 1, tx
       tacc=tacc+1
       !print*,'tacc=',tacc
    do i = 1, ixs
       if ( glons0(i) .lt. 0 ) then
               glons0conv = 360. + glons0(i)
       else
               glons0conv = glons0(i)
       endif
    do j = 1, jxs
         if (tacc.eq.1) then 
             TPorder1(i,j) = 0
             TPorder2(i,j) = 0 
             TPorder3(i,j) = 0 
         endif
       ! === Check ETCs presence in a distance of R from (i, j) at tacc ===

       ETCpresence = .FALSE.
       do n = 1, num(tacc)
       !print*,'n=',n
          if (abs(latc(tacc,n)-glats0(j)).le.12 .and. &
              abs(lonc(tacc,n)-glons0conv).le.20 ) then  
              call calc_dist(glats0(j),glons0(i),latc(tacc,n),lonc(tacc,n),dist12)
              if ( dist12.le.R_threshold ) then
                   ETCpresence = .TRUE. 
              endif
          endif
       enddo               
       ! === For hourly precipitation extreme ===
         if ( tp(i,j,t).gt.TPper1(i,j) .and. TPorder1(i,j).lt.TPpercent_nx1) then
            TPorder1(i,j)                  = TPorder1(i,j)+1  !in time order
            TP_ex1(i,j,TPorder1(i,j))      = tp(i,j,t)
            TP_extime1(i,j,TPorder1(i,j))  = tacc
            if (ETCpresence) ETC_ex1(i,j,TPorder1(i,j)) = 1 
         endif
         if ( tp(i,j,t).gt.TPper2(i,j) .and. TPorder2(i,j).lt.TPpercent_nx2 ) then
            TPorder2(i,j)                  = TPorder2(i,j)+1  !in time order
            TP_ex2(i,j,TPorder2(i,j))      = tp(i,j,t)
            TP_extime2(i,j,TPorder2(i,j))  = tacc
            if (ETCpresence) ETC_ex2(i,j,TPorder2(i,j)) = 1 
         endif
         if ( tp(i,j,t).gt.TPper3(i,j) .and. TPorder3(i,j).lt.TPpercent_nx3) then
            TPorder3(i,j)                  = TPorder3(i,j)+1  !in time order
            TP_ex3(i,j,TPorder3(i,j))      = tp(i,j,t)
            TP_extime3(i,j,TPorder3(i,j))  = tacc
            if (ETCpresence) ETC_ex3(i,j,TPorder3(i,j)) = 1 
         endif
       ! =======  

    enddo ! j-loop
    enddo ! i-loop
    enddo ! t-loop

    deallocate(tp)
    !print*,'deallocation...'
enddo
enddo

do i = 1, ixs
do j = 1, jxs
print*,'TPorder1(i,j)=',TPorder1(i,j)
print*,'TPorder2(i,j)=',TPorder2(i,j)
print*,'TPorder3(i,j)=',TPorder3(i,j)
enddo
enddo


!===Write out data===========
! output data

fileout1='/PrecExtremes_p'//TPpercentile_str1//'&
        valtimeETC_2001_2020_ERA5.nc'
fileout2='/PrecExtremes_p'//TPpercentile_str2//'&
        valtimeETC_2001_2020_ERA5.nc'
fileout3='/PrecExtremes_p'//TPpercentile_str3//'&
        valtimeETC_2001_2020_ERA5.nc'

call writegrid(fileout1, glons0, glats0, &
               TP_ex1, TP_extime1, ETC_ex1, &
               ixs, jxs, TPpercent_nx1, missing_value1)
call writegrid(fileout2, glons0, glats0, &
               TP_ex2, TP_extime2, ETC_ex2, &
               ixs, jxs, TPpercent_nx2, missing_value1)
call writegrid(fileout3, glons0, glats0, &
               TP_ex3, TP_extime3, ETC_ex3, &
               ixs, jxs, TPpercent_nx3, missing_value1)

!stop
contains
!end

!====================================================================================================
subroutine date_to_julian(yyyy, mm, dd, jday)
  implicit none 
  integer, intent ( in) :: yyyy, mm, dd
  integer, intent (out) :: jday
  !converts date to julian date 
  !formula from: https://books.google.ca/books?id=ZwaehQpHHKAC&pg=PA314&lpg=PA314&dq=intrinsic+function+fortran+julian+date&source=bl&ots=LyScHGhs-s&sig=ACfU3U1hDj1TgEDYbJqvQDWtwKhc6Iv2sA&hl=en&sa=X&ved=2ahUKEwi7v6H54Jv2AhUPVt8KHR1iAx0Q6AF6BAgYEAM#v=onepage&q=intrinsic%20function%20fortran%20julian%20date&f=false

  
  jday=dd-32075+1461*(yyyy+4800+(mm-14)/12)/4+367*(mm-2-((mm-14)/12)*12)/12 &
       -3*((yyyy+4900+(mm-14)/12)/100)/4

  return
end 
!----------------------------------------------------------
subroutine julian_to_date(yyyy, mm, dd, jday)
  implicit none
  integer, intent (out) :: yyyy, mm, dd
  integer, intent ( in) :: jday
  integer l, n
  !converts date to julian date
  !formula from: https://books.google.ca/books?id=ZwaehQpHHKAC&pg=PA314&lpg=PA314&dq=intrinsic+function+fortran+julian+date&source=bl&ots=LyScHGhs-s&sig=ACfU3U1hDj1TgEDYbJqvQDWtwKhc6Iv2sA&hl=en&sa=X&ved=2ahUKEwi7v6H54Jv2AhUPVt8KHR1iAx0Q6AF6BAgYEAM#v=onepage&q=intrinsic%20function%20fortran%20julian%20date&f=false

  l = jday + 68569
  n = 4*l/146097
  l = l - (146097*n+3)/4
  yyyy = 4000*(l+1)/1461001
  l = l - 1461*yyyy/4 + 31
  mm = 80*l/2447
  dd = l - 2447*mm/80
  l = mm/11
  mm = mm + 2 - 12*l
  yyyy = 100*(n-49) + yyyy + 1

  return
end
!----------------------------------------------------------
subroutine writegrid(fileout, glons0, glats0, var1, var2, var3, ix, jx, px, missing) 
      
  use netcdf
  implicit none
  character(len=130) fileout
  real ,allocatable, dimension(:,:,:) ::  var1 
  integer, allocatable, dimension(:,:,:) :: var2, var3

  real ,allocatable, dimension(:) :: glats0, glons0
  integer, dimension(3) :: dimids
  integer, dimension(1) :: dimid_ix, dimid_jx
  integer ix, jx, kx, tx, px, sx
  integer ncid5, varid1, varid2, varid3, varid4, varid5 
  integer i_dimid, j_dimid, t_dimid, s_dimid, p_dimid
  integer i_varid, j_varid, t_varid
  real missing

!
!Creat the netCDF file
!  
!  call check(nf90_create(fileout, NF90_CLOBBER, ncid5))
!  call check(nf90_create(fileout,cmode=or(NF90_CLOBBER,NF90_64BIT_OFFSET),&
  call check(nf90_create(fileout,cmode=or(NF90_CLOBBER,NF90_NETCDF4),&
  & ncid=ncid5))
!
!Define the dimensions
!
  call check(nf90_def_dim(ncid5,'longitude', ix, i_dimid))
  call check(nf90_def_dim(ncid5,'latitude', jx, j_dimid))
  call check(nf90_def_dim(ncid5,'percent_num', px, p_dimid))
!
!Define coordinate variables  
!
  call check(nf90_def_var(ncid5,'longitude', NF90_FLOAT, i_dimid, i_varid))
  call check(nf90_def_var(ncid5,'latitude' , NF90_FLOAT, j_dimid, j_varid))
  dimids = (/ i_dimid, j_dimid, p_dimid /)
!
!Define variable 
!
  call check(nf90_def_var(ncid5,'Extreme_values', NF90_FLOAT, dimids, varid1))
  call check(nf90_def_var(ncid5,'Extreme_time',   NF90_INT, dimids, varid2))
  call check(nf90_def_var(ncid5,'ETC_presence',   NF90_INT, dimids, varid3))
  call check(nf90_put_att(ncid5, varid1, "missing_value", missing))
  call check(nf90_put_att(ncid5, varid1, "units", "m/s"))
  call check(nf90_put_att(ncid5, varid2, "definition", "time 1 indicates 1979 Jan 01 00UTC")) 
  call check(nf90_put_att(ncid5, varid3, "definition", "1 indicates one or more ETCs in a distance of 1000 km")) 
 call check(nf90_enddef(ncid5))
!
!Write Data  
!
  call check(nf90_put_var(ncid5, i_varid,  glons0))
  call check(nf90_put_var(ncid5, j_varid,  glats0))
  call check(nf90_put_var(ncid5,  varid1, var1))
  call check(nf90_put_var(ncid5,  varid2, var2))
  call check(nf90_put_var(ncid5,  varid3, var3))
  call check(nf90_close(ncid5))

end subroutine writegrid

!----------------------------------------------------------
!SUBROUTINES (from Katja)

! Subroutine to calculate the distance in kilometers between
! two points (lat1,lon1) and (lat2,lon2) on the sphere

subroutine calc_dist (lat1, lon1, lat2, lon2, dist12)

  implicit none

  real  lat1, lon1, lat2, lon2, dist12
  real  rlat1, rlon1, rlat2, rlon2
  real  x1(3), x2(3)
  real  d, pi, a

! Define some constants

! PI = 3.14159...
  pi = 2.*ASIN(1.)

! Mean radius of the Earth in kilometers
  a  = 6371.22!E3

! Convert to radians
  rlon1=lon1
  rlon2=lon2
  if ( rlon1 .lt. 0 ) rlon1 = 360. + rlon1
  if ( rlon2 .lt. 0 ) rlon2 = 360. + rlon2

  rlat1=lat1*pi/180.
  rlon1=rlon1*pi/180.
  rlat2=lat2*pi/180.
  rlon2=rlon2*pi/180.

! Locate points in Cartesian space

  x1(1)=COS(rlat1)*COS(rlon1)
  x1(2)=COS(rlat1)*SIN(rlon1)
  x1(3)=SIN(rlat1)

  x2(1)=COS(rlat2)*COS(rlon2)
  x2(2)=COS(rlat2)*SIN(rlon2)
  x2(3)=SIN(rlat2)

! Find shortest distance in Cartesian space (divided by Earth's radius)

  d=SQRT((x1(1)-x2(1))**2+(x1(2)-x2(2))**2+(x1(3)-x2(3))**2)

! Find distance following Earth's surface

  dist12 = 2.*a*ASIN(d/2.)

  return

end

!----------------------------------------------------------
subroutine checkopen(rcode,filein)
implicit none
integer::rcode
character(len=120)::filein

 if( rcode == 0) then
    write(6,*) ''
 else
    write(6,*) ' error opening netcdf file ', filein, 'rcode=',rcode
    stop
 endif

end subroutine checkopen
!----------------------------------------------------------
subroutine handle_err(status)
  integer, intent ( in) :: status

  if(status /= nf90_noerr) then
    print *, trim(nf90_strerror(status))
    stop "Stopped"
  end if
end subroutine handle_err
!----------------------------------------------------------
subroutine check(status)
implicit none
integer, intent ( in) :: status

 if(status /= nf90_noerr) then
    print *, trim(nf90_strerror(status))
    stop "Stopped"
 end if

end subroutine check
!=======================================================================

end program LocalGridExtremes
