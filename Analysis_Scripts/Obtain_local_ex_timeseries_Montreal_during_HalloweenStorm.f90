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
integer i, j, p, t, time_index, ys, ye, ms, me, n, e, f, s, tt, oo 
integer num, point, nummax,pointmax
integer line, year_int, month_int
integer tacc, tacc0,totaltx, dirc
parameter (nummax  = 1, pointmax=127)

character*130 fileout1, fileout2, fileout3
character*130 ifile1, ifile2, ifile3, ifile4
character*2 mm, dd, hh, month
character*4 year

integer iun, iun2, iun3, iun4, l_nj, l_ni, ftyp, ezgdef_fmem, gdll
real scale_factor1, add_offset1, missing_value1
real scale_factor2, add_offset2, missing_value2
real scale_factor3, add_offset3, missing_value3
real dx
! Internal parameters
real, parameter :: earthr = 6371.22 ! in kilometers
real            :: deg2rad

real, parameter :: R_threshold  =  2000.    ! in kilometers (1100km used in Owen et al., 2021)
real            :: distc
real, parameter :: missing = -9999.
integer, parameter :: missingint = -9999

! Threshold data
real,       dimension (:)    , allocatable :: glats0, glons0

! Source data
real,       dimension (:,:), allocatable :: u10o, v10o, tpo, u10, v10, wdsp10, tp

! Storm data
integer yyyy_int, dd_int, hh_int, mm_int, yyyy_intpre, mm_intpre
integer,       dimension (:),    allocatable :: pointnum
integer jday0, jday, ic, icmax, icmin, jc, jcmax, jcmin, ic1, ic2, jc1, jc2
integer ixwest, jxnorth, numoriginal
real totalgrid
real latc, lonc

! Output data
real, dimension(:,:,:,:), allocatable:: distocen,tpseries,wdseries, lifetime
integer, dimension(:,:), allocatable:: STORMn
integer, dimension(:,:,:), allocatable:: duration
integer, dimension(:,:,:,:), allocatable:: dirofetc
logical, dimension(:,:,:), allocatable:: STORMcounted
!=========
deg2rad = acos(-1d0)/180d0
!=========

!ix=161 !total 441: from 230 to 340, with 0.25 interval
!jx=81 !total 201: from 25 to 75, with 0.25 interval

! in this program, we target the new NNA region, defined as
! 40-60 N, 
! 260-300


!---
!Sub sample:
ix=521
jx=201

ixwest  = 1   !1145: lon=286  ;1145+6 lon=287.5 
jxnorth = 1    !179:  lat=45.5 ;179+8  lat=47.5


print*,'test...'

!---------------
! Read in the tracked storm data
allocate (pointnum(nummax))

open (10, FILE='HalloweenStorm_track.txt', STATUS='old')

  do line = 1, 127

        read(10,100)  num, numoriginal, point, yyyy_int, mm_int, dd_int, hh_int, &
                      latc, lonc !, vorm(num,point)
        
              pointnum(num) = point

        IF (numoriginal .eq. 11649) THEN   
           print*,'point=',point,'lonc=',lonc 
           if (point.eq.1) then

              print*,'open new files'
              write(year,  '(i4)') yyyy_int
              write(month, '(i2.2)') mm_int

             !Open files for TP, U, V
              ifile1 = 'ERA5_2019OctNov_TP_10mUV.nc'
              rcode = nf90_open(ifile1, nf90_nowrite, ncid1)
              if (rcode /= nf90_noerr) call handle_err(rcode)
              call check(nf90_inq_dimid(ncid1,'time',dimidtx))
              call check(nf90_inquire_dimension(ncid1,dimidtx,dimnametx,tx))
              allocate (glons0(ix))
              allocate (glats0(jx))

              call check(nf90_inq_varid(ncid1,'longitude',varid))
              call check(nf90_get_var(ncid1,varid,glons0, start=(/ixwest/), count=(/ix/)))
              call check(nf90_inq_varid(ncid1,'latitude',varid))
              call check(nf90_get_var(ncid1,varid,glats0, start=(/jxnorth/) , count=(/jx/)))
              do i = 1, ix
              if (glons0(i).lt.0.) glons0(i)=360+glons0(i)
              enddo

              ! [Output variables]
 
              allocate (STORMn(ix,jx))
              allocate (distocen(ix,jx,nummax, pointmax))
              allocate (tpseries(ix,jx,nummax, pointmax))
              allocate (wdseries(ix,jx,nummax, pointmax))
              allocate (lifetime(ix,jx,nummax, pointmax))
              allocate (duration(ix,jx,nummax))
              allocate (STORMcounted(ix,jx,nummax))
              allocate (dirofetc(ix,jx,nummax, pointmax))

              STORMn = 0
              STORMcounted = .FALSE.
              distocen=missing
              tpseries=missing
              wdseries=missing
              lifetime=missing
              duration=missingint
              dirofetc= missingint

           endif

           !time_index=1 corresponds to 01.10.2019 00UTC
           if (mm_int.eq.10) then
               time_index= ((dd_int-1)*24+hh_int)+1
           else
               time_index= ((dd_int-1)*24+hh_int)+31*24+1
           endif
           allocate (tpo(ix,jx))
           allocate (tp(ix,jx))
           !TP:
           call check(nf90_inq_varid(ncid1,'tp',varid))
           call check(nf90_get_var(ncid1,varid,tpo, start=(/ixwest,jxnorth,time_index/), count=(/ix,jx,1/)))
           call check(nf90_get_att(ncid1,varid,'scale_factor',scale_factor1))
           call check(nf90_get_att(ncid1,varid,'add_offset',add_offset1))
           call check(nf90_get_att(ncid1,varid,'missing_value',missing_value1))
           tp=(tpo*scale_factor1+add_offset1)*1000.
           deallocate(tpo)

           allocate (u10o(ix,jx))
           allocate (u10(ix,jx))

           !U10:
           call check(nf90_inq_varid(ncid1,'u10',varid))
           call check(nf90_get_var(ncid1,varid,u10o, start=(/ixwest,jxnorth,time_index/), count=(/ix,jx,1/)))
           call check(nf90_get_att(ncid1,varid,'scale_factor',scale_factor2))
           call check(nf90_get_att(ncid1,varid,'add_offset',add_offset2))
           call check(nf90_get_att(ncid1,varid,'missing_value',missing_value2))
           u10=(u10o*scale_factor2+add_offset2)

           deallocate(u10o)
           allocate (v10o(ix,jx))
           allocate (v10(ix,jx))

           !V10:
           call check(nf90_inq_varid(ncid1,'v10',varid))
           call check(nf90_get_var(ncid1,varid,v10o, start=(/ixwest,jxnorth,time_index/), count=(/ix,jx,1/)))
           call check(nf90_get_att(ncid1,varid,'scale_factor',scale_factor3))
           call check(nf90_get_att(ncid1,varid,'add_offset',add_offset3))
           call check(nf90_get_att(ncid1,varid,'missing_value',missing_value3))
           v10=(v10o*scale_factor3+add_offset3)

           deallocate(v10o)
           allocate (wdsp10(ix,jx))
           wdsp10 = ((u10**2.)+(v10**2.))**0.5

           deallocate(u10)
           deallocate(v10)
       
           !<---Finished collecting the source data, tp, wdsp10, for this current time frame
    
           num=1
           

           do i = 1, ix
           do j = 1, jx

              call calc_dist(glats0(j),glons0(i),latc,lonc,dist12)
              if ( dist12.le.R_threshold ) then
                   lifetime(i,j,num,point)=point

                   ! Calculate the statistics for the local grid cells
                   ! How many storms (the same storm num# is counted once only)
                   ! affect this grid cells?
                   if (Stormcounted(i,j,num)) then
                      duration(i,j,num)=duration(i,j,num)+1
                   else
                      STORMn(i,j)=STORMn(i,j)+1
                      if (duration(i,j,num).eq.missingint) duration(i,j,num)=0
                      duration(i,j,num)=duration(i,j,num)+1
                      Stormcounted(i,j,num)=.TRUE.
                   endif
                      
                   !write out the specific details for targeted grid cells:
                   !IF (i.ge.igs .and. i.le.jge .and. j.ge.jgs .annd. j.le.jge) THEN
                       distocen(i,j,num,point)=dist12 !distance to the storm center
                       tpseries(i,j,num,point)=tp(i,j)
                       wdseries(i,j,num,point)=wdsp10(i,j)
                       call dirloccen(glats0(j),glons0(i),latc,lonc,dirc)
                       dirofetc(i,j,num,point)=dirc
                   !ENDIF

              endif
           enddo
           enddo

           deallocate(tp)
           deallocate(wdsp10)
           
        ENDIF ! (numoriginal .eq. 11649)   

  enddo

  close(10)

call check(nf90_close(ncid1))

  100 format (1x, i5, 3x,i5, 5x, i3, 4x, i4, i2, i2, i2, 4x, f5.2, 3x, f6.2) !,&


!Write out data

print*,'write out data...'

fileout1='Timeseries_WSnTP_local_record_Montreal_2019HalloweenStorm.nc'

call writegrid(fileout1,glons0,glats0, &
               STORMn, duration, &
               distocen, dirofetc, tpseries, wdseries, lifetime, &
               ix, jx, nummax, pointmax)

  deallocate(STORMn)
  deallocate(Stormcounted)
  deallocate(duration)
  deallocate(distocen)
  deallocate(dirofetc)
  deallocate(tpseries)
  deallocate(wdseries)
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
subroutine writegrid(fileout, glons0, glats0, var1, var2, &
                      var3, var4, var5,var6, var7, ix, jx, nx, px) 
      
  use netcdf
  implicit none
  character(len=130) fileout
  integer, allocatable, dimension(:,:) ::  var1
  integer, allocatable, dimension(:,:,:) ::  var2
  integer, allocatable, dimension(:,:,:,:) ::  var4
  real, allocatable, dimension(:,:,:,:) ::  var3, var5, var6, var7

  real ,allocatable, dimension(:) :: glats0, glons0
  integer, dimension(2) :: dimids2
  integer, dimension(3) :: dimids3
  integer, dimension(4) :: dimids4
  integer, dimension(1) :: dimid_ix, dimid_jx
  integer ix, jx, kx, tx, px, nx
  integer ncid5, varid1, varid2, varid3, varid4, varid5, varid6, varid7
  integer i_dimid, j_dimid, t_dimid, n_dimid, p_dimid
  integer i_varid, j_varid, t_varid
  real, parameter :: missing = -9999.
  integer, parameter :: missingint = -9999
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
  call check(nf90_def_dim(ncid5,'stormnum_newNNA', nx, n_dimid))
  call check(nf90_def_dim(ncid5,'stormpoint', px, p_dimid))
!
!Define coordinate variables  
!
  call check(nf90_def_var(ncid5,'longitude', NF90_FLOAT, i_dimid, i_varid))
  call check(nf90_def_var(ncid5,'latitude' , NF90_FLOAT, j_dimid, j_varid))
  dimids2 = (/ i_dimid, j_dimid /)
  dimids3 = (/ i_dimid, j_dimid, n_dimid /)
  dimids4 = (/ i_dimid, j_dimid, n_dimid, p_dimid /)
!Define variable 
!
  call check(nf90_def_var(ncid5,'Storm_count', NF90_INT, dimids2, varid1))
  call check(nf90_def_var(ncid5,'Duration_of_storm', NF90_INT, dimids3, varid2))
  call check(nf90_def_var(ncid5,'Distance_to_storm', NF90_FLOAT, dimids4, varid3))
  call check(nf90_def_var(ncid5,'Re_direction_of_storm', NF90_INT, dimids4, varid4))
  call check(nf90_def_var(ncid5,'TP_series', NF90_FLOAT, dimids4, varid5))
  call check(nf90_def_var(ncid5,'WS_series', NF90_FLOAT, dimids4, varid6))
  call check(nf90_def_var(ncid5,'lifetime', NF90_FLOAT, dimids4, varid7))
  call check(nf90_put_att(ncid5, varid2, "missing_value", missingint))
  call check(nf90_put_att(ncid5, varid3, "missing_value", missing))
  call check(nf90_put_att(ncid5, varid4, "missing_value", missingint))
  call check(nf90_put_att(ncid5, varid4, "title", "0=N,the storm is located to the N of this grid cell; 6=W"))
  call check(nf90_put_att(ncid5, varid5, "missing_value", missing))
  call check(nf90_put_att(ncid5, varid6, "missing_value", missing))
  call check(nf90_put_att(ncid5, varid7, "missing_value", missing))
  call check(nf90_enddef(ncid5))
!Write Data  
!
  call check(nf90_put_var(ncid5, i_varid,  glons0))
  call check(nf90_put_var(ncid5, j_varid,  glats0))
  call check(nf90_put_var(ncid5,  varid1, var1))
  call check(nf90_put_var(ncid5,  varid2, var2))
  call check(nf90_put_var(ncid5,  varid3, var3))
  call check(nf90_put_var(ncid5,  varid4, var4))
  call check(nf90_put_var(ncid5,  varid5, var5))
  call check(nf90_put_var(ncid5,  varid6, var6))
  call check(nf90_put_var(ncid5,  varid7, var7))
  call check(nf90_close(ncid5))

end subroutine writegrid
!----------------------------------------------------------
subroutine dirloccen (latg, long, latetc, lonetc, dirc)
  !find out the relative direction of the storm with respect to the grid cell: 
  !if "W" (west), it means the storm is located to the west of the grid cell (i,j)
  !if dirc = 0 :N
  !if dirc = 1 :NE
  !if dirc = 2 :E
  !if dirc = 3 :SE
  !if dirc = 4 :S
  !if dirc = 5 :SW
  !if dirc = 6 :W
  !if dirc = 7 :NW
  !if dirc = 8 :C (center overlappting with the grid point)

  implicit none

  real  latg, long, latetc, lonetc
  integer dirc
 
  real degree
  real, parameter :: pi=3.14159

  if ((latg-latetc).eq.0 .and. (long-lonetc).eq.0) then
    dirc = 8
  else
    degree= atan2(latetc-latg,lonetc-long)*180./pi
    if (degree.gt.-22.5 .and. degree.le.22.5) then
        dirc = 2
    elseif (degree.gt.22.5 .and. degree.le.67.5) then
        dirc = 1
    elseif (degree.gt.67.5 .and. degree.le.112.5) then
        dirc = 0
    elseif (degree.gt.112.5 .and. degree.le.157.5) then
        dirc = 7
    elseif (degree.gt.157.5 .or. degree.le.-157.5) then
        dirc = 6
    elseif (degree.gt.-157.5 .and. degree.le.-112.5) then
        dirc = 5
    elseif (degree.gt.-112.5 .and. degree.le.-67.5) then
        dirc = 4
    elseif (degree.gt.-67.5 .and. degree.le.-22.5) then
        dirc = 3
    endif
  endif


end subroutine dirloccen       
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
