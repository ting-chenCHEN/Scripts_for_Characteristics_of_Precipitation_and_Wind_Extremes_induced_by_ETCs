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
integer e, i, j, t, time_index, ys, ye, ms, me, n, p, s, tt, oo 
integer num, point, nummax,pointmax
integer line, linemax, lineheader, year_int, month_int
integer tacc, totaltx
integer yyyy_int, mm_int, dd_int
integer jday0, jday
real TPpercentile

character*130 fileout0, fileout1, fileout2, fileout3, fileout4, &
              fileout5, fileout6 
character*130 ifile, ifile1, ifile2, ifile3, ifile4, &
              ifilein1, filein2, filein3, &
              ifilein1a, filein2a, filein3a 
character*2 mm, dd, hh, month
character*4 year
character*4 TPpercentile_str

integer iun, iun2, iun3, iun4, l_nj, l_ni, ftyp, ezgdef_fmem, gdll
real scale_factor1, add_offset1, missing_value1
real scale_factor2, add_offset2, missing_value2
real scale_factor3, add_offset3, missing_value3
real dx
logical global, ll_double
! Internal parameters
real, parameter :: earthr = 6371.22 ! in kilometers
real            :: deg2rad

! Threshold data
real,       dimension (:,:)  , allocatable :: TPper
! Input data
real,       dimension (:)    , allocatable :: glats0, glons0
real,       dimension (:,:,:), allocatable :: TPex, TPtime

integer TPpercent_nx


! Output data
real,      dimension (:,:), allocatable :: &
        Exceedance, Excount
  
real,      dimension (:,:,:), allocatable :: &
        seaExceedance, seaExcount
!=========
deg2rad = acos(-1d0)/180d0
!=========

!===Setup parameters===

  TPpercentile = 99

  TPpercentile_str = '99p0'
  

  jxs = 201
  ixs = 521

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

  TPpercent_nx = 1+int((100.-TPpercentile)/100*float(totaltx))

  print*,'TPpercent_num:', TPpercent_nx

!======================

! input data: Percentile values in TP

ifile = 'TPp'//TPpercentile_str//'_2001_2020_ERA5.nc'
ifilein1='PrecExtremes_p'//TPpercentile_str//'valtimeETC_2001_2020_ERA5.nc'

!------------------------------------
! Read in the threshold value at each grid point
!  allocate (WDper(ixs, jxs))
  allocate (TPper(ixs, jxs))

! Open input files
  rcode = nf90_open(ifile, nf90_nowrite, ncid1)
  if (rcode /= nf90_noerr) call handle_err(rcode)
  rcode = nf90_open(ifilein1, nf90_nowrite, ncid2)
  if (rcode /= nf90_noerr) call handle_err(rcode)

  allocate (glons0(ixs))
  allocate (glats0(jxs)) 
  call check(nf90_inq_varid(ncid2,'longitude',varid))
  call check(nf90_get_var(ncid2,varid,glons0))
  call check(nf90_inq_varid(ncid2,'latitude',varid))
  call check(nf90_get_var(ncid2,varid,glats0))
  
  allocate (TPex(ixs,jxs,TPpercent_nx))
  allocate (TPtime(ixs,jxs,TPpercent_nx))
 
  !TP:
  call check(nf90_inq_varid(ncid1,'p'//TPpercentile_str//'',varid))
  call check(nf90_get_var(ncid1,varid,TPper))

  call check(nf90_inq_varid(ncid2,'Extreme_values',varid))
  call check(nf90_get_var(ncid2,varid,TPex))

  call check(nf90_inq_varid(ncid2,'Extreme_time',varid))
  call check(nf90_get_var(ncid2,varid,TPtime))

  call check(nf90_close(ncid1))
  call check(nf90_close(ncid2))
  !call check(nf90_close(ncid3))

!---------------
! [Output variables]

  allocate (seaExceedance(ixs,jxs,4))
  allocate (Exceedance(ixs,jxs))
  allocate (seaExcount(ixs,jxs,4))
  allocate (Excount(ixs,jxs))
  seaExceedance=0.
  Exceedance=0.
  seaExcount=0.
  Excount=0.

!--------------------------------------------------------------------------
print*,'TPpercent_nx=',TPpercent_nx
do e = 1, TPpercent_nx

    do i = 1, ixs
    do j = 1, jxs
     
            if (TPex(i,j,e).ge. 0.) then
                ! Identify the season:
              jday = int((TPtime(i,j,e)-1)/24)+jday0
              call julian_to_date( yyyy_int, mm_int, dd_int, jday)
    
              if (mm_int.eq.6  .or. mm_int.eq.7 .or. mm_int.eq.8  ) s=1
              if (mm_int.eq.9  .or. mm_int.eq.10.or. mm_int.eq.11 ) s=2
              if (mm_int.eq.12 .or. mm_int.eq.1 .or. mm_int.eq.2  ) s=3
              if (mm_int.eq.3  .or. mm_int.eq.4 .or. mm_int.eq.5  ) s=4

              seaExcount(i,j,s) = seaExcount(i,j,s)+1.
              seaExceedance(i,j,s) =  seaExceedance(i,j,s) + TPex(i,j,e)-TPper(i,j)
          
              Excount(i,j) = Excount(i,j)+1.
              Exceedance(i,j) =  Exceedance(i,j) + TPex(i,j,e)-TPper(i,j)
            !else
            ! print*,'non positive values for tpex:',TPex(i,j,e),'at e=',e,'i=',i,'j=',j
 
            endif

    enddo ! j-loop
    enddo ! i-loop

enddo ! e-loop

do i = 1, ixs
do j = 1, jxs

do s = 1, 4
   seaExceedance(i,j,s) = seaExceedance(i,j,s) / seaExcount(i,j,s)
enddo
   Exceedance(i,j) = Exceedance(i,j) / Excount(i,j)

enddo
enddo


! ===========================================

!===Write out data===========

fileout1='Averaged_exceedance_TP_p'//TPpercentile_str//'_ERA5.nc'

call writegrid(fileout1,glons0,glats0, &
               seaExcount, seaExceedance, &
               Excount, Exceedance, & 
               ixs, jxs, 4)


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
subroutine writegrid(fileout, glons0, glats0, seaExcount,seaExceedance, Excount, Exceedance, &
                ix, jx, sx) 
      
  use netcdf
  implicit none
  character(len=130) fileout
  real ,allocatable, dimension(:,:,:) ::  seaExcount, seaExceedance
  real, allocatable, dimension(:,:) :: Excount, Exceedance

  real ,allocatable, dimension(:) :: glats0, glons0
  integer, dimension(2) :: dimids
  integer, dimension(3) :: dimids2, dimids3
  integer, dimension(1) :: dimid_ix, dimid_jx
  integer ix, jx, kx, tx, px, sx
  integer ncid5, varid1, varid2, varid3, varid4, varid5 
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
  dimids = (/ i_dimid, j_dimid /)
  dimids2 = (/ i_dimid, j_dimid, s_dimid /)
!Define variable 
  call check(nf90_def_var(ncid5,'seaExcount', NF90_FLOAT, dimids2, varid1))
  call check(nf90_def_var(ncid5,'seaExceedance', NF90_FLOAT, dimids2, varid2))
  call check(nf90_def_var(ncid5,'Excount', NF90_FLOAT, dimids, varid3))
  call check(nf90_def_var(ncid5,'Exceedance', NF90_FLOAT, dimids, varid4))
  call check(nf90_enddef(ncid5))
!Write Data  
!
  call check(nf90_put_var(ncid5, i_varid,  glons0))
  call check(nf90_put_var(ncid5, j_varid,  glats0))
  call check(nf90_put_var(ncid5,  varid1, seaExcount))
  call check(nf90_put_var(ncid5,  varid2, seaExceedance))
  call check(nf90_put_var(ncid5,  varid3, Excount))
  call check(nf90_put_var(ncid5,  varid4, Exceedance))
  call check(nf90_close(ncid5))

end subroutine writegrid
!----------------------------------------------------------
subroutine writegrida(fileout, glons0, glats0, var2,  ix, jx, px) 
      
  use netcdf
  implicit none
  character(len=130) fileout
  real ,allocatable, dimension(:,:,:) ::  var1
  integer ,allocatable, dimension(:,:,:) ::  var2

  real ,allocatable, dimension(:) :: glats0, glons0
  integer, dimension(3) :: dimids2, dimids
  integer, dimension(1) :: dimid_ix, dimid_jx
  integer ix, jx, kx, tx, px, sx
  integer ncid5, varid1, varid2, varid3
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
  call check(nf90_def_dim(ncid5,'percent_num', px, p_dimid))
!
!Define coordinate variables  
!
  call check(nf90_def_var(ncid5,'longitude', NF90_FLOAT, i_dimid, i_varid))
  call check(nf90_def_var(ncid5,'latitude' , NF90_FLOAT, j_dimid, j_varid))
  dimids = (/ i_dimid, j_dimid, s_dimid /)
  dimids2 = (/ i_dimid, j_dimid, p_dimid /)
!Define variable 
!
!  call check(nf90_def_var(ncid5,'longitude', NF90_FLOAT, dimid_ix, varid1))
!  call check(nf90_def_var(ncid5,'latitude ', NF90_FLOAT, dimid_jx, varid2))
  call check(nf90_def_var(ncid5,'Extreme_time', NF90_FLOAT, dimids2, varid2))
  call check(nf90_enddef(ncid5))
!Write Data  
!
  call check(nf90_put_var(ncid5, i_varid,  glons0))
  call check(nf90_put_var(ncid5, j_varid,  glats0))
  call check(nf90_put_var(ncid5,  varid2, var2))
  call check(nf90_close(ncid5))

end subroutine writegrida
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
