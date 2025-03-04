program LocalGridExtremes

! ulimit -t unlimited 
!
! Compile with:
!
! module load compiler/gcc-7.3 development/netcdf
! gfortran $(nf-config --fflags) $(nf-config --flibs) XXX.f90 -o XXX.Abs
use netcdf
implicit none

integer ncid1,ncid1a,ncid1b,ncid2, ncid3, ncid4, ncid5, ncid6
integer rcode, dimidix, dimidjx, dimidtx, varid, ix, ixs, jx, jxs, tx
character(len=3) :: dimnameix
character(len=3) :: dimnamejx
character(len=4) :: dimnametX
character(len=9) :: dimnameix2
character(len=8) :: dimnamejx2
real beta,dist12,ppoints4, ppoints6, ppoints8, ppoints10
integer s_pre,i, j, p, t, time_index, ys, ye, ms, me, n, e, f, s, tt, ee, oo 
integer num, point, nummax,pointmax, num_pre, num_pre0
integer line, linemax, lineheader, year_int, month_int, linenum_for2021
integer TPpercent_nx, tacc, tacc0,totaltx, dirc
parameter (linemax = 934010, lineheader = 16) !<- to start from the storm#630(included) after 2001
parameter (linenum_for2021=45539) !<- first in 20210101
parameter (nummax  = 12359, pointmax=525)

character*130 fileout0, fileout1, fileout2, fileout3
character*130 ifile1, ifile1a, ifile1b, ifile2, ifile3, ifile4, ifile5, ifile6
character*2 mm, dd, hh, month
character*4 year
character*4 TPpercentile_str, WDpercentile_str

integer iun, iun2, iun3, iun4, l_nj, l_ni, ftyp, ezgdef_fmem, gdll
real scale_factor1, add_offset1, missing_value1
real scale_factor2, add_offset2, missing_value2
real scale_factor3, add_offset3, missing_value3
real dx, expoint
! Internal parameters
real, parameter :: earthr = 6371.22 ! in kilometers
real            :: deg2rad

real, parameter :: R_threshold  =  1000.    ! in kilometers (1100km used in Owen et al., 2021)
real            :: distc
real, parameter :: missing = -9999.
real, parameter :: prinv = 0.5 ! interval for precipitation exceedance PDF (0.5*72 total)
integer, parameter :: missingint = -9999
integer, parameter :: tscales= 6 !1h,3h,6h,12h,24h, 48h
integer, parameter :: timeranges= 72 !1h->36h
integer, parameter :: intranges= 72 !1h->36h
integer :: startcal 
! Threshold data
real,       dimension (:,:)  , allocatable :: TPper1
real,       dimension (:)    , allocatable :: glats0, glons0

! Source data
real,       dimension (:,:), allocatable :: tpo, tp1, TPduration, TPduration_avg, addnum_all, TPexceedance_avg
real, dimension(:,:,:), allocatable :: TPduration_season, addnum_season, TPexceedance_season

! Storm data
integer ETC_extreme, num_extreme
integer yyyy_int, dd_int, hh_int, mm_int, yyyy_intpre, mm_intpre
integer jday0, jday, ic, icmax, icmin, jc, jcmax, jcmin, ic1, ic2, jc1, jc2
integer ixwest, jxnorth, numoriginal
real totalgrid
real latc, lonc

! Output data
real, dimension(:,:), allocatable:: cont, TPdupdf, TPdupdf_C,TPdupdf_season, TPdupdf_season_C, &
                                          TPexpdf, TPexpdf_season
real, dimension(:), allocatable:: TPdupdf_avg, TPdupdf_avg_C, TPexpdf_avg
integer, dimension(:), allocatable:: num_extreme_sea
!=========
deg2rad = acos(-1d0)/180d0
!=========

num_pre0=11649
num_pre=num_pre0
s = 3  !<- start from DJF


!---
!Sub sample:
ix=521
jx=201

ixwest  = 841   !lon=230 to 240 (included)
jxnorth =  61   !lat=75 to 25 (included) 

!===Setup parameters===
ys=2019
ye=2020
ms=10
me=12

ETC_extreme=0


! Total time steps (hourly data) during the entire period
! and the corresponding sample numbers for the selected threshold
  totaltx     = 7305*24
  print*,'totaltx:', totaltx

    ! For all extremes files
  TPpercentile_str='99p0'

  ! Remeber to change the threhsold files BELOW!!!
  ! ---->
! input data: Percentile thresholds in TP and WDSP
  ifile1 = 'TPp'//TPpercentile_str//'_2001_2020_ERA5.nc'

! Read in the threshold value at each grid point

  rcode = nf90_open(ifile1, nf90_nowrite, ncid1)
  if (rcode /= nf90_noerr) call handle_err(rcode)

  call check(nf90_inq_dimid(ncid1,'longitude',dimidix))
  call check(nf90_inquire_dimension(ncid1,dimidix,dimnameix,ix))
  call check(nf90_inq_dimid(ncid1,'latitude',dimidjx))
  call check(nf90_inquire_dimension(ncid1,dimidjx,dimnamejx,jx))

!  print*,'ix=',ix,'jx=',jx

  allocate (TPper1(ix, jx))

  !TP:
  call check(nf90_inq_varid(ncid1,'p'//TPpercentile_str//'',varid))
  call check(nf90_get_var(ncid1,varid,TPper1))

  call check(nf90_close(ncid1))

!---------------

  ! input data: Storm location
  open (10, FILE='ETC_identified_n_tracking_output_2000_2020.txt', STATUS='old')
  yyyy_intpre = 1999 !random start
  mm_intpre = 1 !random start

  startcal = 0

  do line = 1, linemax

    if (line .le. lineheader) then
        read(10,*)     !header
    else
        read(10,100)  num, point, yyyy_int, mm_int, dd_int, hh_int, &
                      latc, lonc !, vorm(num,point)
        if (num.lt.num_pre0 .or. num.gt.num_pre0+1) cycle 

        IF (yyyy_int.eq.yyyy_intpre .and. mm_int.eq.mm_intpre) THEN
        ! no need to open a new ifile (same ETC track) 
        ELSE
           startcal=startcal+1
           yyyy_intpre=yyyy_int
           mm_intpre=mm_int

           if (yyyy_int.ne.ys.and.mm_int.ne.ms) then
              call check(nf90_close(ncid1a))
           endif

           write(year,  '(i4)') yyyy_int
           write(month, '(i2.2)') mm_int

           if (num.eq.num_pre0 .or. num.ne.num_pre) then
             s_pre = s
             if (mm_int.eq.6  .or. mm_int.eq.7 .or. mm_int.eq.8  ) s=1
             if (mm_int.eq.9  .or. mm_int.eq.10.or. mm_int.eq.11 ) s=2
             if (mm_int.eq.12 .or. mm_int.eq.1 .or. mm_int.eq.2  ) s=3
             if (mm_int.eq.3  .or. mm_int.eq.4 .or. mm_int.eq.5  ) s=4
           endif

           !Open files for TP, U, V
           print*,'open new files'
           ifile1a = '/home/archive/REANALYSES/ERA5/1h/tp/ll/nc4/'//year//'/'//month//'/era5_tp_ll_'//year//''//month//'_1h.nc4'

           rcode = nf90_open(ifile1a, nf90_nowrite, ncid1a)
           if (rcode /= nf90_noerr) call handle_err(rcode)

           call check(nf90_inq_dimid(ncid1a,'time',dimidtx))
           call check(nf90_inquire_dimension(ncid1a,dimidtx,dimnametx,tx))
           if (startcal.eq.1) then
              print*, 'reading glons0 glats0...'
              allocate (glons0(ix))
              allocate (glats0(jx))
              call check(nf90_inq_varid(ncid1a,'longitude',varid))
              call check(nf90_get_var(ncid1a,varid,glons0, start=(/ixwest/), count=(/ix/)))
              call check(nf90_inq_varid(ncid1a,'latitude',varid))
              call check(nf90_get_var(ncid1a,varid,glats0, start=(/jxnorth/) , count=(/jx/)))

              ! [Output variables]
 
              allocate (TPduration(ix,jx))
              allocate (TPduration_avg(ix,jx))
              allocate (TPduration_season(ix,jx,4))
              allocate (TPexceedance_avg(ix,jx))
              allocate (TPexceedance_season(ix,jx,4))
              allocate (addnum_all(ix,jx))
              allocate (addnum_season(ix,jx,4))
              allocate (cont(ix,jx))
              allocate (TPdupdf(nummax,timeranges))
              allocate (TPexpdf(nummax,timeranges))
              allocate (TPdupdf_C(nummax,timeranges))
              TPduration = missing
              TPduration_avg = missing
              TPduration_season= missing
              TPexceedance_avg = missing
              TPexceedance_season= missing
              TPdupdf = 0.
              TPexpdf = 0.
              TPdupdf_C = 0.
              cont=0.

              allocate (TPdupdf_avg(timeranges))
              allocate (TPexpdf_avg(timeranges))
              allocate (TPdupdf_season(timeranges,4))
              allocate (TPexpdf_season(timeranges,4))
              allocate (TPdupdf_avg_C(timeranges))
              allocate (TPdupdf_season_C(timeranges,4))
              TPdupdf_avg = 0.
              TPdupdf_season = 0.
              TPexpdf_avg = 0.
              TPexpdf_season = 0.
              TPdupdf_avg_C = 0.
              TPdupdf_season_C = 0.

              num_extreme=0
              allocate (num_extreme_sea(4))
              num_extreme_sea=0
           endif

        ENDIF 

           
           time_index= ((dd_int-1)*24+hh_int)+1
           allocate (tpo(ix,jx))
           allocate (tp1(ix,jx))
           !U10:
           call check(nf90_inq_varid(ncid1a,'tp',varid))
           call check(nf90_get_var(ncid1a,varid, tpo, start=(/ixwest,jxnorth,time_index/), count=(/ix,jx,1/)))
           call check(nf90_get_att(ncid1a,varid,'scale_factor',scale_factor1))
           call check(nf90_get_att(ncid1a,varid,'add_offset',add_offset1))
           call check(nf90_get_att(ncid1a,varid,'missing_value',missing_value1))
           tp1= (tpo*scale_factor1+add_offset1)*1000.
           deallocate(tpo)

       
           !<---Finished collecting the source data, tp, wdsp10, for this current time frame

           !------------------------------------------------------------------------
           !After reading the entire lifetime of EACH (previous) ETC:
           if (num.ne.num_pre) then
               
               do i = 1, ix
               do j = 1, jx

                  !Only the cyclone that caused extremes are counted here...
                  !====== FIRST, WE STORE EXTRME DURATION AND THE CORRESPONDING EXCEEDANCE ======
                  !here, the cont(i,j) is already the accumulated exceedance at the same grid point over the lifetime of this ETC
                  !when the extreme is caused. So cont(i,j)/TPduration(i,j) is the temporal averaged intensity of the exceedance 
                  !for a given grid point during this ETC.
                   
                  IF (glats0(j).le.60.) THEN !To be compared with IMERG data... only applied to the pdf calculation
                  
                  do tt = 1, timeranges-1

                      if (TPduration(i,j).gt.float(tt-1) .and. TPduration(i,j).le.float(tt)) then
                          TPdupdf(num_pre,tt)=TPdupdf(num_pre,tt)+1/expoint                                  !PDF of duration (the
                                                                                                     !accumulated grid points of the
                                                                                                     !same duration) 

                          TPdupdf_C(num_pre,tt)=TPdupdf_C(num_pre,tt)+cont(i,j)/expoint              !accumulated exceedance
                      endif

                  enddo

                  if (TPduration(i,j).ge.float(timeranges)) then
                      TPdupdf(num_pre,timeranges)=TPdupdf(num_pre,timeranges)+1/expoint
                      TPdupdf_C(num_pre,timeranges)=TPdupdf_C(num_pre,timeranges)+cont(i,j)/expoint
                  endif
                  
                  !====== SECOND, WE CALULATE PDF OF EXCEEDANCE ======
                  do ee = 1, intranges-1

                      if ((cont(i,j)/TPduration(i,j).gt.float(ee-1)*prinv) .and. &
                          (cont(i,j)/TPduration(i,j).le.float(ee)*prinv)) then
                          TPexpdf(num_pre,ee)=TPexpdf(num_pre,ee)+1/expoint                                  !PDF of exceedance
                      endif

                  enddo
                   
                  if (cont(i,j)/TPduration(i,j).ge.float(intranges)*prinv) then
                          TPexpdf(num_pre,intranges)=TPexpdf(num_pre,intranges)+1/expoint                     !PDF of exceedance
                  endif

                  ENDIF ! for regions to the south of 60N

                  !======= THIRD, WE ACCUMULATE SOME DIAGNOSTICS FOR MORE THAN ONE ETCS ======
                  ! for maps output 
                  if (TPduration(i,j).ge.0.) then
                   if (TPduration_avg(i,j).eq.missing) then
                           TPduration_avg(i,j)=0.
                           TPexceedance_avg(i,j)=0.
                           addnum_all(i,j)=0.      ! the total number of ECTs attributing to the extremes for this local grid
                   endif
                   TPduration_avg(i,j)=TPduration_avg(i,j)+TPduration(i,j)
                   TPexceedance_avg(i,j)=TPexceedance_avg(i,j)+cont(i,j)
                   addnum_all(i,j)=addnum_all(i,j)+1.    

                   if (TPduration_season(i,j,s_pre).eq.missing) then
                           TPduration_season(i,j,s_pre)=0.
                           TPexceedance_season(i,j,s_pre)=0.
                           addnum_season(i,j,s_pre)=0.
                   endif
                   TPduration_season(i,j,s_pre)=TPduration_season(i,j,s_pre)+TPduration(i,j)
                   TPexceedance_season(i,j,s_pre)=TPexceedance_season(i,j,s_pre)+cont(i,j)
                   addnum_season(i,j,s_pre)=addnum_season(i,j,s_pre)+1.
                  endif
                  !-------------------------

                  
                  !reset after finishing ALL timesteps for A GIVEN ETC
                  TPduration(i,j)=missing
                  cont(i,j)=0.
                  

               enddo
               enddo
               
               expoint =0.!reset
              
                     

           endif 
           !------------------------------------------------------------------------


           do i = 1, ix
           do j = 1, jx
 
              call calc_dist(glats0(j),glons0(i),latc,lonc,dist12)
              if ( dist12.le.R_threshold ) then
                   if (tp1(i,j).gt. TPper1(i,j)) then
                       IF (glats0(j).le.60. .and. TPduration(i,j).eq.missing) THEN !To be compared with IMERG data... only applied to the pdf calculation
                               expoint = expoint +1. !total points in the domain covered by extremes given an ETC(only ocunted once during one
                                                 !ETC event)
                       ENDIF
                       if (TPduration(i,j).eq.missing) TPduration(i,j)=0.
                       TPduration(i,j)=TPduration(i,j)+1.
                       cont(i,j)=cont(i,j)+tp1(i,j)-TPper1(i,j) ! accumulated/total exceedance during the passage of one ETC
                       ETC_extreme=1 !1 means true; 0 means false
                   endif
              endif

           enddo
           enddo


           !========ADDING the PDF over many ETCs in different seasons (will take averages at the end)============
           
           if (num.ne.num_pre) then

               do t = 1, timeranges
               TPdupdf_avg(t)=TPdupdf_avg(t)+TPdupdf(num_pre,t)
               TPexpdf_avg(t)=TPexpdf_avg(t)+TPexpdf(num_pre,t)
               TPdupdf_season(t,s_pre)=TPdupdf_season(t,s_pre)+TPdupdf(num_pre,t)
               TPexpdf_season(t,s_pre)=TPexpdf_season(t,s_pre)+TPexpdf(num_pre,t)
               TPdupdf_avg_C(t)=TPdupdf_avg_C(t)+TPdupdf_C(num_pre,t)
               TPdupdf_season_C(t,s_pre)=TPdupdf_season_C(t,s_pre)+TPdupdf_C(num_pre,t)
               enddo

               if (ETC_extreme.eq.1) then
                   num_extreme=num_extreme+1
                   num_extreme_sea(s_pre)=num_extreme_sea(s_pre)+1
                   !reset
                   ETC_extreme=0
               endif
               
           endif
           !-------------------------------------------------------------------

           num_pre=num
           
           deallocate(tp1)

    endif !line 
  enddo

  close(10)
  print*,'close(10)'



  do i = 1, ix
  do j = 1, jx
     if (TPduration_avg(i,j).ne.missing) then
     TPduration_avg(i,j)=TPduration_avg(i,j)/addnum_all(i,j)     !Averaged duration of extreme on each grid over all ETCs
     TPexceedance_avg(i,j)=TPexceedance_avg(i,j)/addnum_all(i,j)   !Averaged exeedance of extreme on each grid over all ETCs
     endif
     do s = 1, 4
     if (TPduration_season(i,j,s).ne.missing) then
     TPduration_season(i,j,s)=TPduration_season(i,j,s)/addnum_season(i,j,s)
     TPexceedance_season(i,j,s)=TPexceedance_season(i,j,s)/addnum_season(i,j,s)
     endif
     enddo
  enddo
  enddo

  do t = 1, timeranges
     TPdupdf_avg(t)=TPdupdf_avg(t)/float(num_extreme)
     TPexpdf_avg(t)=TPexpdf_avg(t)/float(num_extreme)
     TPdupdf_avg_C(t)=TPdupdf_avg_C(t)/float(num_extreme)
     do s = 1, 4
     if (num_extreme_sea(s).ne.0) then
     TPdupdf_season(t,s)=TPdupdf_season(t,s)/float(num_extreme)
     TPexpdf_season(t,s)=TPexpdf_season(t,s)/float(num_extreme)
     TPdupdf_season_C(t,s)=TPdupdf_season_C(t,s)/float(num_extreme)
     endif
     enddo
  enddo

  100 format (3x, i5, 2x, i3, 11x, i4, i2, i2, i2, 3x, f5.2, 2x, f6.2, 2x )

               
!----------------------------------------------------------


!Write out data

print*,'write out data...'
fileout0='TPExtremes_p'//TPpercentile_str//'_duration_exceedance_distribution_Halloween_ERA5.nc'
fileout1='PDF_TPexceedance_duration_p'//TPpercentile_str//'_NNA_HalloweenStorm.txt'

call writegrid(fileout0,glons0,glats0, TPduration_avg, TPduration_season, TPexceedance_avg, TPexceedance_season, ix, jx, 4)

open(11,file=fileout1,status='replace')

do t=1,timeranges

write(11,*) TPdupdf_avg(t), TPdupdf_avg_C(t), TPexpdf_avg(t)

enddo


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
subroutine writegrid(fileout, glons0, glats0, var1, var2, var3, var4,&
                     ix, jx, sx) 
      
  use netcdf
  implicit none
  character(len=130) fileout
  real, allocatable, dimension(:,:) ::  var1, var3
  real, allocatable, dimension(:,:,:) ::  var2, var4

  real ,allocatable, dimension(:) :: glats0, glons0
  integer, dimension(2) :: dimids2
  integer, dimension(3) :: dimids3
  integer, dimension(4) :: dimids4
  integer, dimension(1) :: dimid_ix, dimid_jx
  integer ix, jx, kx, tx, px, nx, sx
  integer ncid5, varid1, varid2, varid3, varid4, varid5, varid6
  integer i_dimid, j_dimid, s_dimid, n_dimid, p_dimid
  integer i_varid, j_varid, s_varid
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
  call check(nf90_def_dim(ncid5,'season', sx, s_dimid))
!
!Define coordinate variables  
!
  call check(nf90_def_var(ncid5,'longitude', NF90_FLOAT, i_dimid, i_varid))
  call check(nf90_def_var(ncid5,'latitude' , NF90_FLOAT, j_dimid, j_varid))
  dimids2 = (/ i_dimid, j_dimid /)
  dimids3 = (/ i_dimid, j_dimid, s_dimid /)
!Define variable 
!
  call check(nf90_def_var(ncid5,'extreme_duration_all', NF90_FLOAT, dimids2, varid1))
  call check(nf90_def_var(ncid5,'extreme_duration_season', NF90_FLOAT, dimids3, varid2))
  call check(nf90_def_var(ncid5,'extreme_exceedance_all', NF90_FLOAT, dimids2, varid3))
  call check(nf90_def_var(ncid5,'extreme_exceedance_season', NF90_FLOAT, dimids3, varid4))
  call check(nf90_put_att(ncid5, varid1, "missing_value", missing))
  call check(nf90_put_att(ncid5, varid2, "missing_value", missing))
  call check(nf90_put_att(ncid5, varid3, "missing_value", missing))
  call check(nf90_put_att(ncid5, varid4, "missing_value", missing))
  call check(nf90_enddef(ncid5))
!Write Data  
!
  call check(nf90_put_var(ncid5, i_varid,  glons0))
  call check(nf90_put_var(ncid5, j_varid,  glats0))
  call check(nf90_put_var(ncid5,  varid1, var1))
  call check(nf90_put_var(ncid5,  varid2, var2))
  call check(nf90_put_var(ncid5,  varid3, var3))
  call check(nf90_put_var(ncid5,  varid4, var4))
 
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
