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
integer jday0, jday0_ans, jday0_ane, jday
integer yyyy_int, mm_int, dd_int
character*130 fileout1, fileout2, fileout3
character*130 ifile1, ifile2, ifile3, ifile4
character*2 mm, dd, hh, month
character*4 year
character*4 TPpercentile_str, WDpercentile_str
character*3 TPpercentile_str0, WDpercentile_str0

integer iun, iun2, iun3, iun4, l_nj, l_ni, ftyp, ezgdef_fmem, gdll
real scale_factor1, add_offset1, missing_value1
real scale_factor2, add_offset2, missing_value2
real scale_factor3, add_offset3, missing_value3
real dx
logical global, ll_double
logical :: ws6htp, ETCpresence
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
real,       dimension (:,:,:), allocatable :: WDex, TPex, WDcases_sea, TPcases_sea, COcases_sea
integer,    dimension (:,:,:), allocatable :: WDtime, TPtime, WDetc, TPetc

integer TPpercent_nx, WDpercent_nx

! Output data
real,      dimension (:,:,:), allocatable :: &
                                             Compound_ETC_ws6htp,  &
                                             Compound_noETC_ws6htp

!=========
deg2rad = acos(-1d0)/180d0
!=========

!===Setup parameters===
  ys=2001
  ye=2020
  ms=1
  me=12

! Julian day for the first day of the interested period
  call date_to_julian(2001,1,1,jday0) !where the extreme data time started 
  ! WDtime = 1 indicates 01-01-2001 00 UTC
  ! WDtime-hh equals jday0
  ! WDtime = x (original value) -> x-hh+jday0 (converted to real julian dayy)
  call date_to_julian(ys,ms, 1,jday0_ans) 
  call date_to_julian(ye,me,31,jday0_ane) 

! Total time steps (hourly data) during the entire period
! and the corresponding sample numbers for the selected threshold
  totaltx     = 7305*24

  print*,'totaltx:', totaltx

  ! For all extremes files
  TPpercentile_str0='990' 
  WDpercentile_str0='990'
  TPpercentile_str='99p0' 
  WDpercentile_str='99p0'

  ! Remeber to change the threhsold files BELOW!!!
  ! ---->
! input data: Percentile thresholds in TP and WDSP
  ifile1 = 'WSp'//WDpercentile_str0//'_2001_2020_ERA5.nc'
  ifile2 = 'TPp'//TPpercentile_str0//'_2001_2020_ERA5.nc'

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

  print*,'percentile threshold reading finished'
! input data: All extreme values in TP and WDSP
ifile1 = 'WindExtremes_1h_p'&
        //WDpercentile_str//'valtimeETC_2001_2020_ERA5.nc'
ifile2 = 'PrecExtremes_1h_p'&
        //TPpercentile_str//'valtimeETC_2001_2020_ERA5.nc'
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
  
  if (WDpercent_nx.ne.TPpercent_nx) print*,'Inconsistent number of extremes between WDSP and TP. Double check!'

  allocate (WDex(ix,jx,WDpercent_nx))
  allocate (WDtime(ix,jx,WDpercent_nx))
  allocate (WDetc(ix,jx,WDpercent_nx))
  allocate (TPex(ix,jx,TPpercent_nx))
  allocate (TPtime(ix,jx,TPpercent_nx))
  allocate (TPetc(ix,jx,TPpercent_nx))

  !WDSP:
  call check(nf90_inq_varid(ncid1,'Extreme_values',varid))
  call check(nf90_get_var(ncid1,varid,WDex))
  call check(nf90_inq_varid(ncid1,'Extreme_time',varid))
  call check(nf90_get_var(ncid1,varid,WDtime))
  call check(nf90_inq_varid(ncid1,'ETC_presence',varid))
  call check(nf90_get_var(ncid1,varid,WDetc))

  !TP:
  call check(nf90_inq_varid(ncid2,'Extreme_values',varid))
  call check(nf90_get_var(ncid2,varid,TPex))
  call check(nf90_inq_varid(ncid2,'Extreme_time',varid))
  call check(nf90_get_var(ncid2,varid,TPtime))
  call check(nf90_inq_varid(ncid2,'ETC_presence',varid))
  call check(nf90_get_var(ncid2,varid,TPetc))
  
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

! [Output variables]

! From station
  
  allocate (Compound_ETC_ws6htp(ix,jx,4))    !total number of compound events 
  !(if ETC-associated WS occurs at t and there's at least a TP occurring in the following 6hours)
  allocate (Compound_noETC_ws6htp(ix,jx,4))    !total number of compound events 
  !(if nonETC-associated TP occurs at t and there's at least a WS occurring in the following 6hours)

  Compound_ETC_ws6htp =0.
  Compound_noETC_ws6htp =0.
  
  do i = 1, ix
  do j = 1, jx

!-----------------------------------------
     do e = 1, WDpercent_nx

       ws6htp = .FALSE.
       ETCpresence = .FALSE.

       ! Convert the extremem timestep to julian day:
       jday = int((WDtime(i,j,e)-1)/24)+jday0
       ! Only calculate when the time of extreme is 
       ! after  2001-01-01 00UTC (period of interes)
       ! before 2020-12-31 23UTC (period of interes)
       IF (jday.ge.jday0_ans .and. jday.le.jday0_ane+23) THEN

       ! Identify the season:
         call julian_to_date( yyyy_int, mm_int, dd_int, jday)
              if (mm_int.eq.6  .or. mm_int.eq.7 .or. mm_int.eq.8  ) s=1
              if (mm_int.eq.9  .or. mm_int.eq.10.or. mm_int.eq.11 ) s=2
              if (mm_int.eq.12 .or. mm_int.eq.1 .or. mm_int.eq.2  ) s=3
              if (mm_int.eq.3  .or. mm_int.eq.4 .or. mm_int.eq.5  ) s=4

         !======================
         do f = 1, TPpercent_nx

            IF ( TPtime(i,j,f).ge. (WDtime(i,j,e)-6) .and. &
               & TPtime(i,j,f).lt.(WDtime(i,j,e)+6) ) THEN 
                     ws6htp = .TRUE.
                     if (TPetc(i,j,f).gt.0.5 .or. WDetc(i,j,e).gt.0.5) ETCpresence=.TRUE.
            ENDIF
         enddo ! f-loop 
         !======================

       IF (ws6htp) THEN
           if (ETCpresence) then
             Compound_ETC_ws6htp(i,j,s)= Compound_ETC_ws6htp(i,j,s)+1
           else
             Compound_noETC_ws6htp(i,j,s)= Compound_noETC_ws6htp(i,j,s)+1
           endif
       ENDIF


       ENDIF !Period after 2001 
                
     enddo ! e-loop 
!-----------------------------------------
  enddo ! j-loop
  enddo ! i-loop
 

  print*, 'Writing out data... :'

!===Write out data===========

fileout1='CompoundExtremes_p'//WDpercentile_str//'&
        _ETC_association_relaxto12h_2001_2020_ERA5.nc'

call writegrid(fileout1,glons0,glats0,  &
               Compound_ETC_ws6htp, &
               Compound_noETC_ws6htp, &
               ix, jx, 4)

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
subroutine writegrid(fileout, glons0, glats0, var1, var2,  ix, jx, sx) 
      
  use netcdf
  implicit none
  character(len=130) fileout
  real ,allocatable, dimension(:,:,:) ::  var1, var2
  real ,allocatable, dimension(:) :: glats0, glons0
  integer, dimension(3) :: dimids
  integer, dimension(1) :: dimid_ix, dimid_jx
  integer ix, jx, kx, tx, px, sx
  integer ncid5, varid1, varid2, varid3, varid4, varid5, varid6 
  integer i_dimid, j_dimid, t_dimid, s_dimid, p_dimid
  integer i_varid, j_varid, t_varid
!
!Creat the netCDF file
!  
!  call check(nf90_create(fileout, NF90_CLOBBER, ncid5))
  call check(nf90_create(fileout,cmode=or(NF90_CLOBBER,NF90_64BIT_OFFSET),&
  & ncid=ncid5))
!
!Define the dimensions
!
  call check(nf90_def_dim(ncid5,'longitude', ix, i_dimid))
  call check(nf90_def_dim(ncid5,'latitude', jx, j_dimid))
  call check(nf90_def_dim(ncid5,'season', sx, s_dimid))
!
!Define coordinate variables  
!
  call check(nf90_def_var(ncid5,'longitude', NF90_FLOAT, i_dimid, i_varid))
  call check(nf90_def_var(ncid5,'latitude' , NF90_FLOAT, j_dimid, j_varid))
  dimids = (/ i_dimid, j_dimid, s_dimid /)
!Define variable 
!
  call check(nf90_def_var(ncid5,'ETC_Compound_counts_12h',   NF90_FLOAT, dimids, varid1))
  call check(nf90_def_var(ncid5,'noETC_Compound_counts_12h', NF90_FLOAT, dimids, varid2))
  call check(nf90_enddef(ncid5))
!Write Data  
!
  call check(nf90_put_var(ncid5, i_varid,  glons0))
  call check(nf90_put_var(ncid5, j_varid,  glats0))
  call check(nf90_put_var(ncid5,  varid1, var1))
  call check(nf90_put_var(ncid5,  varid2, var2))
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

! ==============================================================================
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
