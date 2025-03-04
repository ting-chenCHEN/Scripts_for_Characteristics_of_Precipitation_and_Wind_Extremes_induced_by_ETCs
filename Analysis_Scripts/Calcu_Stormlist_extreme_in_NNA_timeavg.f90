program LocalGridExtremes
!
! Compile with:
!
! module load compiler/gcc-7.3 development/netcdf
! gfortran $(nf-config --fflags) $(nf-config --flibs) XXX.f90 -o XXX.Abs
use netcdf
implicit none

integer ncid1,ncid2, ncid3, ncid4
integer rcode, dimidix, dimidjx, dimidtx, varid, ix, ixs, jx, jxs, tx
character(len=3) :: dimnameix
character(len=3) :: dimnamejx
character(len=4) :: dimnametX
character(len=9) :: dimnameix2
character(len=8) :: dimnamejx2
real dist12,ppoints4, ppoints6, ppoints8, ppoints10
integer i, j, p, t, time_index, ys, ys0, ye, ms, me, n, e, f, s, tt, oo 
integer num, point, nummax,pointmax
integer line, linemax, lineheader, year_int, month_int
integer tacc, totaltx, TPpercentile, WDpercentile
parameter (linemax = 934013, lineheader = 16)
parameter (nummax  = 12359, pointmax=525)

character*130 fileout1, fileout2, fileout3
character*130 ifile1, ifile2, ifile3, ifile4
character*2 mm, dd, hh, month
character*4 year
character*4 TPpercentile_str, WDpercentile_str

integer iun, iun2, iun3, iun4, l_nj, l_ni, ftyp, ezgdef_fmem, gdll
real scale_factor1, add_offset1, missing_value1
real scale_factor2, add_offset2, missing_value2
real scale_factor3, add_offset3, missing_value3
real dx, NNA_area
logical global, ll_double
! Internal parameters
real, parameter :: earthr = 6371.22 ! in kilometers
real            :: deg2rad


real, parameter :: R_threshold  =  1000.    ! in kilometers (1100km used in Owen et al., 2021)
real            :: distc
real, parameter :: Di_threshold  =  55.     ! in degree (set looser)
real, parameter :: Dj_threshold  =  15.     ! in degree (set looser)

! Threshold data
real,       dimension (:,:)  , allocatable :: TPper, WDper
real,       dimension (:)    , allocatable :: glats0, glons0
real,       dimension (:,:,:), allocatable :: WDex, TPex
integer,    dimension (:,:,:), allocatable :: WDtime, TPtime

integer TPpercent_nx, WDpercent_nx

! Storm data
integer yyyy_int, dd_int, hh_int, mm_int, dd2
real,          dimension (:,:),  allocatable :: latc, lonc
integer,       dimension (:,:),  allocatable :: tindex 
integer,       dimension (:),    allocatable :: pointnum, nummaxupto, nummin
real,          dimension (:),    allocatable :: Wex_avg, Pex_avg, Wex_area, Pex_area
real,          dimension (:),    allocatable :: pmin_peak,vors_peak
real,          dimension (:,:,:),allocatable :: Wex_mask, Pex_mask
real,          dimension (:,:,:),allocatable :: Wextreme_weighted, Pextreme_weighted
integer jday0, jday, ic, icmax, icmin, jc, jcmax, jcmin, ic1, ic2, jc1, jc2
integer ixwest, jxnorth
real totalgrid
real pmin, vorsavg
character*37 intensity
logical Wcounted, Pcounted, Ccounted, COMPOUND



!=========
deg2rad = acos(-1d0)/180d0
!=========

!===Setup parameters===
  ys=2001
  ye=2020
  ms=1
  me=12

! Julian day for the first day of the interested period
  call date_to_julian(ys,ms,1,jday0) 

! Total time steps (hourly data) during the entire period
! and the corresponding sample numbers for the selected threshold
  totaltx     = 7305*24

  print*,'totaltx:', totaltx

  ! For all extremes files
  TPpercentile_str='99p0' 
  WDpercentile_str='99p0'

  ! Remeber to change the threhsold files BELOW!!!
  ! ---->
! input data: Percentile thresholds in TP and WDSP
  ifile1 = 'WSp'//WDpercentile_str//'_2001_2020_ERA5.nc'
  ifile2 = 'TPp'//TPpercentile_str//'_2001_2020_ERA5.nc'

! Read in the threshold value at each grid point
  
  rcode = nf90_open(ifile1, nf90_nowrite, ncid1)
  if (rcode /= nf90_noerr) call handle_err(rcode)
  rcode = nf90_open(ifile2, nf90_nowrite, ncid2)
  if (rcode /= nf90_noerr) call handle_err(rcode)

  call check(nf90_inq_dimid(ncid1,'longitude',dimidix))
  call check(nf90_inquire_dimension(ncid1,dimidix,dimnameix,ix))
  call check(nf90_inq_dimid(ncid1,'latitude',dimidjx))
  call check(nf90_inquire_dimension(ncid1,dimidjx,dimnamejx,jx))


  allocate (WDper(ix, jx))
  allocate (TPper(ix, jx))

  !WDSP:
  call check(nf90_inq_varid(ncid1,'p'//WDpercentile_str//'',varid))
  call check(nf90_get_var(ncid1,varid,WDper))
  !TP:
  call check(nf90_inq_varid(ncid2,'p'//TPpercentile_str//'',varid))
  call check(nf90_get_var(ncid2,varid,TPper))

  call check(nf90_close(ncid1))
  call check(nf90_close(ncid2))


! input data: All extreme values in TP and WDSP
ifile1 = 'WindExtremes_p'//TPpercentile_str//'valtimeETC_2001_2020_ERA5.nc'
ifile2 = 'PrecExtremes_p'//TPpercentile_str//'valtimeETC_2001_2020_ERA5.nc'
!------------------------------------
! Read in the threshold value at each grid point

! Open input files
  rcode = nf90_open(ifile1, nf90_nowrite, ncid1)
  if (rcode /= nf90_noerr) call handle_err(rcode)
  rcode = nf90_open(ifile2, nf90_nowrite, ncid2)
  if (rcode /= nf90_noerr) call handle_err(rcode)

  call check(nf90_inq_dimid(ncid1,'percent_num',dimidtx))
  call check(nf90_inquire_dimension(ncid1,dimidtx,dimnametx,WDpercent_nx))
  call check(nf90_inq_dimid(ncid2,'percent_num',dimidtx))
  call check(nf90_inquire_dimension(ncid2,dimidtx,dimnametx,TPpercent_nx))
  
  allocate (WDex(ix,jx,WDpercent_nx))
  allocate (WDtime(ix,jx,WDpercent_nx))
  allocate (TPex(ix,jx,TPpercent_nx))
  allocate (TPtime(ix,jx,TPpercent_nx))

  !WDSP:
  call check(nf90_inq_varid(ncid1,'Extreme_values',varid))
  call check(nf90_get_var(ncid1,varid,WDex))
  call check(nf90_inq_varid(ncid1,'Extreme_time',varid))
  call check(nf90_get_var(ncid1,varid,WDtime))

  !TP:
  call check(nf90_inq_varid(ncid2,'Extreme_values',varid))
  call check(nf90_get_var(ncid2,varid,TPex))
  call check(nf90_inq_varid(ncid2,'Extreme_time',varid))
  call check(nf90_get_var(ncid2,varid,TPtime))
  
  allocate (glons0(ix)) 
  allocate (glats0(jx)) 
  call check(nf90_inq_varid(ncid1,'longitude',varid))
  call check(nf90_get_var(ncid1,varid,glons0))
  call check(nf90_inq_varid(ncid1,'latitude',varid))
  call check(nf90_get_var(ncid1,varid,glats0)) 

  call check(nf90_close(ncid1))
  call check(nf90_close(ncid2))
  
  print*,'extreme values reading finished'

!---------------
! Read in the tracked storm data
  allocate (latc(nummax, pointmax))
  allocate (lonc(nummax, pointmax))
  allocate (tindex(nummax, pointmax))
  allocate (pointnum(nummax))
  allocate (pmin_peak(nummax))
  allocate (vors_peak(nummax))
  allocate (nummaxupto(totaltx))
  allocate (nummin(totaltx))
  nummaxupto = 1
  nummin = 1

  pmin_peak = 1100  !(hPa)
  vors_peak = 0        
  ! input data: Storm location
  open (10, FILE='ETC_identified_n_tracking_output_2000_2020.txt', STATUS='old')

  do line = 1, linemax

    if (line .le. lineheader) then
        read(10,*)     !header
    else
        read(10,100)  num, point, yyyy_int, mm_int, dd_int, hh_int, &
                      latc(num,point), lonc(num,point), pmin, vorsavg
        pointnum(num) = point
        !convert the DATE TIME to the same time index starting from ys, ms for later purpose
        call date_to_julian(yyyy_int,mm_int,dd_int,jday) 
        tindex(num,point)= 1+(jday-jday0)*24+hh_int   !tindex=1 indicates 2001 Jan 01 00UTC

        if (tindex(num,point).lt.1) then
        else
          nummaxupto(tindex(num,point))=num
          if (nummin(tindex(num,point)).eq.1) nummin(tindex(num,point))=num
        endif

        if (pmin.lt.pmin_peak(num))    pmin_peak(num)=pmin
        if (vorsavg.gt.vors_peak(num)) vors_peak(num)=vorsavg
      
  endif
  enddo

  close(10)


  100 format (3x, i5, 2x, i3, 11x, i4, i2, i2, i2, 3x, f5.2, 2x, f6.2, 15x, &
             es11.5, 15x, es11.5) 
           

  print*,'storm data reading finished'
! stop

! [Output variables]
! From the perspective of each tracked storm
  allocate (Wex_avg(nummax))
  allocate (Pex_avg(nummax))
  allocate (Wex_area(nummax))
  allocate (Pex_area(nummax))
  allocate (Wex_mask(ix,jx,nummax))
  allocate (Pex_mask(ix,jx,nummax))
  allocate (Wextreme_weighted(ix, jx, nummax))
  allocate (Pextreme_weighted(ix, jx, nummax))
  Wex_mask=0.
  Pex_mask=0.
  Wex_area=0.
  Pex_area=0.
  Wex_avg=0.
  Pex_avg=0.
  Wextreme_weighted=0.
  Pextreme_weighted=0.

  !------------------------------
  ! Targeted NNA region : 
  ! lat : 40 - 60N
  ! lon : 105 - 60W ( 255 - 300) 
  !------------------------------

  do i = 1, ix
  do j = 1, jx

  if (i.eq.1.and.j.eq.1) print*,'start calculating the association with ETCs...'
     ! === For winds ===

     IF (glats0(j).ge.40. .and. glats0(j).le.60. .and. glons0(i).ge.255 .and. glons0(i).le.300.) THEN

     do e = 1, WDpercent_nx

       Wcounted=.FALSE.
       ! Associated with one or more ETC?
         do n = nummin(WDtime(i,j,e)), nummaxupto(WDtime(i,j,e))
         do p = 1, pointnum(n)
              if (tindex(n,p).eq.WDtime(i,j,e)) then !At the same time
              call calc_dist(glats0(j),glons0(i),latc(n,p),lonc(n,p),dist12)
              if ( dist12.le.R_threshold ) then
                 !------------------
                 !--For stroms record
                    if (WDex(i,j,e).gt.0.) then
                    Wex_mask(i,j,n)= Wex_mask(i,j,n) +1.  
                    Wextreme_weighted(i,j,n)=Wextreme_weighted(i,j,n)+&
                                          (WDex(i,j,e)-WDper(i,j))
                    endif
                 !------------------
               endif
               endif
         enddo ! p-loop
         enddo ! n-loop

     enddo ! e-loop 

  if (i.eq.ix.and.j.eq.jx) print*,'finish calculating for wind events...'
     ! === For precipitation (only) ===
     do e = 1, TPpercent_nx
       ! Associated with one or more ETC?
         do n = nummin(TPtime(i,j,e)), nummaxupto(TPtime(i,j,e))
         do p = 1, pointnum(n)
              if (tindex(n,p).eq.TPtime(i,j,e)) then !At the same time
              call calc_dist(glats0(j),glons0(i),latc(n,p),lonc(n,p),dist12)
              if ( dist12.le.R_threshold ) then
                 !------------------
                 !--For stroms record
                    if (TPex(i,j,e).gt.0.) then
                    Pex_mask(i,j,n)=Pex_mask(i,j,n) +1. 
                    Pextreme_weighted(i,j,n)=Pextreme_weighted(i,j,n)+&
                                          (TPex(i,j,e)-TPper(i,j))
                    endif
                 !------------------

               endif
               endif
         enddo ! p-loop
         enddo ! n-loop
     enddo ! e-loop 
 
     ENDIF

  enddo
  enddo
 


! ===========================================
  print*, 'all initial calculation ends...'

  NNA_area =0.

  do i = 1, ix
  do j = 1, jx
        IF (glats0(j).ge.40. .and. glats0(j).le.60. .and. glons0(i).ge.255 .and. glons0(i).le.300.) THEN
           NNA_area=NNA_area+1.
        ENDIF

  enddo
  enddo

  do n = 1, nummax

     Wex_area(n)=0.
     Pex_area(n)=0.

     do i = 1, ix
     do j = 1, jx
        
        if (Wex_mask(i,j,n).gt.0.5) Wex_area(n)=Wex_area(n)+1.
        if (Pex_mask(i,j,n).gt.0.5) Pex_area(n)=Pex_area(n)+1.
     
        if (Wex_mask(i,j,n).gt.0.5) Wex_avg(n)=Wex_avg(n)+Wextreme_weighted(i,j,n)/Wex_mask(i,j,n)/NNA_area
        if (Pex_mask(i,j,n).gt.0.5) Pex_avg(n)=Pex_avg(n)+Pextreme_weighted(i,j,n)/Pex_mask(i,j,n)/NNA_area
        
     enddo
     enddo

  enddo
  
  print*, 'Writing out data... :'

!===Write out data===========
! output data

open (20, FILE='/pampa/chen/local_extremes/&
        Stormlist_inNNA_2001_2020_exceedance_WDp'//WDpercentile_str//'&
        _TPp'//TPpercentile_str//'_timeavg.txt',&
        STATUS='replace')


do n = 1, nummax
        write(20,300) n, pmin_peak(n), vors_peak(n),Wex_avg(n), Pex_avg(n), Wex_area(n), Pex_area(n)
enddo

close(20)


300 format (1x, i5, 2x, f7.2, 2x, es11.4, 2x, f8.5, 2x, f8.5, 2x, f7.0, 2x, f7.0 &
           )



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
