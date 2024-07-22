!!=============================================================
!!=============================================================
!!module: atm_external_atmosphere
!! NAME
!!  atm_external_atmosphere (RModolo)
!!
!! Read and load an external atmosphere file
!!
!! NOTE
module atm_external_atmosphere
 use defs_basis
 use defs_species
 use defs_atmospheretype
! use defs_parametre
 use m_writeout
#ifdef HAVE_NETCDF 
    use defs_basic_cdf
    use netcdf
#endif


#include "q-p_common.h"

 implicit none
 
 contains
 
 !!=============================================================
 !!subroutine: external_atmosphere/load_atmosphere_mars_LMD
 !!
 !! FUNCTION
 !!  Read and load the Martian atmosphere from the GCM LMD model (Forget et al, 1999)
 !!   
 subroutine load_atmosphere_mars_LMD(Spe,ncm,gstep,s_min_loc,resistivity,density_O,density_CO2)
 use defs_variable,only : viscosity
 use defs_parametre,only :dt,atmospherename
 use defs_mpitype,only     : mpiinfo

 integer, intent(in) :: ncm(3)
 real(dp),intent(in) :: gstep(3),s_min_loc(3)
 real(dp),intent(inout) :: resistivity(:,:,:),density_O(:,:,:),density_CO2(:,:,:)
 type(species_type),intent(in) :: Spe 
 character(len=500) :: msg    
 integer :: nAlti, ncid,StId,nTime,nLong,nLati,nHoru,varfieldId,nvar(3),varid
 integer :: nAltiKp1,i,ilat,ilon,iHoru,ir,ip,it,ialt,idone,ilatSub,ilonSub,idone0
 integer,dimension(nf90_max_var_dims) :: DimfieldId

 real(dp),allocatable :: Time(:), Latitude(:), Longitude(:) ,longsubsol(:)
 real(dp),allocatable :: zls(:), dsm(:), dec(:)
 real(dp),allocatable :: sza(:,:) 
 real(dp),dimension(:,:,:),allocatable :: Altitude,nCO2,nO,nCO2p,nO2p,Tn, &
                                          Temp_elec,Temp_elec0
 integer :: ii,jj,kk
 integer :: ilat_GCM,ilon_GCM,ialt_GCM
 real(dp) :: x_cdr,y_cdr,z_cdr,ss(3),Xpc,Ypc,radius,r_cdr
 real(dp),save :: Obliqui = 0.43 ! Obliquitee of the planet
 real(dp) :: clock, szamin,cs0, cs, cmcs, coscmcs,cosclock
 real(dp) :: diff_Lat,diff_Long,diff_Alt,Lon_GEO,Lat_GEO,Alt_GEO
!real(dp),intent(out),dimension(1:nr_GCM,1:nt_GCM,1:np_GCM) :: nCO2_n,nO_n,nO2_p,nCO2_p,T_n,T_e                                          
 
  __WRT_DEBUG_IN("load_atmosphere_mars_LMD")
 
 call wrt_double(6,"Reading file : "//atmospherename,wrtscreen,wrtdisk)
print *,'Atmosphere file :',atmospherename
!-- Open NetCDF file
stId = nf90_open(trim(atmospherename),nf90_nowrite, ncid)
call test_cdf(stId)
call get_simple_dimens_cdf(ncid,"Time",nTime)
call get_simple_dimens_cdf(ncid,"Longitude",nLong)
call get_simple_dimens_cdf(ncid,"Latitude",nLati)

nHoru = nLong * nLati
StId = nf90_inq_varId(ncid,"Altitude",varfieldId)
StId = nf90_Inquire_Variable(ncid,varfieldId,dimids = DimfieldId)
StId = nf90_Inquire_Dimension(ncid,DimfieldId(1), len=nvar(1))
StId = nf90_Inquire_Dimension(ncid,DimfieldId(2), len=nvar(2))
StId = nf90_Inquire_Dimension(ncid,DimfieldId(3), len=nvar(3))

nAlti = nvar(2)
  write(msg,'(2a,4(a,a17,i12))')&
       & ch10," ___________________ Atmosphere file Information  _________________",&
       & ch10, "   nAlti  = ",nAlti,&
       & ch10, "   nLong  = ",nLong,&
       & ch10, "   nLati  = ",nLati,&
       & ch10, "   nTime  = ",nTime
call wrt_double(qp_out,msg,wrtscreen,wrtdisk)
allocate(Time(nTime))
call get_simple_variable_cdf(ncid,"Time",Time(:))
  write(msg,'(2a,2(a,a17,e12.5))')&
       & ch10," ___________________ Atmosphere Time  _________________",&
       & ch10, "   Start  = ",minval(Time(:)),&
       & ch10, "   Stop  = ",maxval(Time(:))
call wrt_double(qp_out,msg,wrtscreen,wrtdisk)       
allocate(zls(nTime))
call get_simple_variable_cdf(ncid,"zls",zls(:))
write(msg,'(2a,2(a,a20,e12.5))')&
       & ch10," ______________ Atmosphere Solar Longitude (Ls)  _____________",&
       & ch10, " Start (degrees)  = ",minval(zls(:)*180.0_dp/pi),&
       & ch10, " Stop             = ",maxval(zls(:)*180.0_dp/pi)
call wrt_double(qp_out,msg,wrtscreen,wrtdisk)
allocate(dsm(nTime))
call get_simple_variable_cdf(ncid,"dsm",dsm(:))
write(msg,'(2a,2(a,a20,e12.5))')&
       & ch10," ___________________ Atmosphere Sun distance  _________________",&
       & ch10, " Start (AU)  = ",dsm(1),&
       & ch10, " Stop        = ",dsm(nTime)
call wrt_double(qp_out,msg,wrtscreen,wrtdisk)
allocate(dec(nTime))
call get_simple_variable_cdf(ncid,"dec",dec(:))
write(msg,'(2a,2(a,a20,e12.5))')&
       & ch10," ___________________ Atmosphere Declination  _________________",&
       & ch10, " Start (AU)  = ",dec(1),&
       & ch10, " Stop        = ",dec(nTime)
call wrt_double(qp_out,msg,wrtscreen,wrtdisk)
allocate(sza(nHoru,nTime))
call get_simple_variable_cdf(ncid,"mu0",sza(:,:))
sza(1:nHoru,1:nTime) = acos(sza(1:nHoru,1:nTime))
write(msg,'(2a,2(a,a20,e12.5))')&
       & ch10," ___________________ Atmosphere SZA  _________________",&
       & ch10, " Min (AU)  = ",minval(sza(:,:)),&
       & ch10, " Max       = ",maxval(sza(:,:))
call wrt_double(qp_out,msg,wrtscreen,wrtdisk)
allocate(longsubsol(nTime))
call get_simple_variable_cdf(ncid,"longsubsol",longsubsol(:))
write(msg,'(2a,2(a,a20,e12.5))')&
         & ch10," ___________________ Atmosphere Subsolar longitude  _________",&
         & ch10, " Start (AU)  = ",longsubsol(1),&
         & ch10, " Stop        = ",longsubsol(nTime)
call wrt_double(qp_out,msg,wrtscreen,wrtdisk)
allocate(Longitude(nLong))
call get_simple_variable_cdf(ncid,"longitude",Longitude(:))
write(msg,'(2a,2(a,a20,e12.5))')&
       & ch10," ___________________ Atmosphere Longitude  _________________",&
       & ch10, " Min (degrees)  = ",minval(Longitude(:)),&
       & ch10, " Max            = ",maxval(Longitude(:))
call wrt_double(qp_out,msg,wrtscreen,wrtdisk)
allocate(Latitude(nLati))
call get_simple_variable_cdf(ncid,"Latitude",Latitude(:))
write(msg,'(2a,2(a,a20,e12.5))')&
       & ch10," ___________________ Atmosphere Latitude  _________________",&
       & ch10, " Min (degrees)  = ",minval(Latitude(:)),&
       & ch10, " Max            = ",maxval(Latitude(:))
call wrt_double(qp_out,msg,wrtscreen,wrtdisk)
allocate(Altitude(nHoru,nAlti,nTime))
stId = nf90_inq_varid(ncid,'Altitude', varid)
call test_cdf(stId)
call get_simple_variable_cdf(ncid,"Altitude",Altitude)
write(msg,'(2a,2(a,a20,e12.5))')&
       & ch10," ___________________ Atmosphere Altitude  _________________",&
       & ch10, " Min (km)  = ",minval(Altitude(:,:,:)),&
       & ch10, " Max       = ",maxval(Altitude(:,:,:))
call wrt_double(qp_out,msg,wrtscreen,wrtdisk)
allocate(nCO2(nHoru,nAlti,nTime))
call get_simple_variable_cdf(ncid,"co2",nCO2)
write(msg,'(2a,2(a,a20,e12.5))')&
       & ch10," ___________________ Atmosphere nCO2  _________________",&
       & ch10, " Min (cm-3)  = ",minval(nCO2(:,:,:)),&
       & ch10, " Max         = ",maxval(nCO2(:,:,:))
call wrt_double(qp_out,msg,wrtscreen,wrtdisk)
allocate(nO(nHoru,nAlti,nTime))
call get_simple_variable_cdf(ncid,"o",nO)
write(msg,'(2a,2(a,a20,e12.5))')&
       & ch10," ___________________ Atmosphere nO  _________________",&
       & ch10, " Min (cm-3)  = ",minval(nO(:,:,:)),&
       & ch10, " Max         = ",maxval(nO(:,:,:))
call wrt_double(qp_out,msg,wrtscreen,wrtdisk)

Latitude(1:nLati)  = Latitude(1:nLati)*pi/180.0_dp
Longitude(1:nLong) = Longitude(1:nLong)*pi/180.0_dp
stId = nf90_close(ncid)
print *,'End reading atmosphere file :',atmospherename
! Philisophy
! for each grid pont (MSO or simulation coord)
! we check if the altitude is lower than 200 km (exobase)
! if no => nothing to do, the density is determined by the exospheric model
! if yes = > from the MSO position we convert it in GEO coord.
!		we find the closest grid point in the GCM input and affect it to the density

! first we determine the location of subsolar latitude and longitude of the GCM
! we also determoine angle cs
! we use only the last GCM time diagnostic (nTime)


cosclock = cos(Obliqui)/cos(dec(nTime))
if (abs(cosclock) <= 1.) then
  clock   = acos(cos(Obliqui)/cos(dec(nTime))) 
else
  if (cosclock > 1) clock = 0.
  if (cosclock < -1) clock = pi
endif

szamin = 3.0_dp*pi
cs0    = 2._dp*pi*(0.5_dp - Time(nTime))
if (cs0.gt. pi) cs0 = cs0 - 2._dp*pi
if (cs0.lt.-pi) cs0 = cs0 + 2._dp*pi
do ilat=1,nLati
  do ilon=1,nLong
    iHoru   = ilon + (ilat - 1)*nLong
    if (sza(iHoru,nTime).lt.szamin) then
      szamin  = sza(iHoru,nTime)
      ilatSub = ilat
      ilonSub = ilon
!   write(6,'(2(1x,i3),3(1x,e12.5))')ilatSub,ilonSub,szamin,Latitude(ilat),Longitude(ilon)
    endif
!  write(6,'(2(1x,i3),3(1x,e12.5))')ilat,ilon,sza(iHoru,iTime)*180.0/Pi,Latitude(ilat)*180.0/Pi,Longitude(ilon)*180.0/Pi
   if (cos(Latitude(ilat)).gt.1.E-06) then
       coscmcs = (cos(sza(iHoru,nTime)) - sin(Latitude(ilat))*sin(dec(nTime))) &
     & / (cos(Latitude(ilat))*cos(dec(nTime)))
!  write(6,'(5(1x,e12.5))')dcos(sza(iHoru,iTime)),dsin(Latitude(ilat)),dsin(dec(iTime)),dcos(Latitude(ilat)),dcos(dec(iTime))
       cmcs    = acos(coscmcs)
       if (Longitude(ilon)-Longitude(ilonSub).gt.pi) cmcs = 2.0*Pi - cmcs
       cs = Longitude(ilon) - cmcs
       if (cs.gt.pi)  cs = pi
       if (cs.lt.-pi) cs = -pi             
    endif 
  enddo
enddo

!print *,'cs,cs0,longsubsol(nTime)',cs,cs0,longsubsol(nTime)
cs0 = longsubsol(nTime)
print *,'clock,cs0,dec(nTime)',clock,cs0,dec(nTime)

! reset to zero previously charged neutral density
ss(:) = zero
density_O(:,:,:) = zero
density_CO2(:,:,:) = zero


 do kk = 1,ncm(3)-1
   do jj = 1,ncm(2)-1
     do ii = 1,ncm(1)-1
           ss(3) = (real((kk-1),dp))*gstep(3) + s_min_loc(3) 
           ss(2) = (real((jj-1),dp))*gstep(2) + s_min_loc(2)
           ss(1) = (real((ii-1),dp))*gstep(1) + s_min_loc(1)      

           radius = sqrt(dot_product(ss-Spe%P%centr,ss-Spe%P%centr))

!!!!!!!! Change by F. Leblanc 01/2015 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          if ((radius <= Spe%P%radius+220.d0/Spe%ref%c_omegapi).and.(radius >= Spe%P%radius)) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             x_cdr = -(ss(1)-Spe%P%centr(1))
             y_cdr = -(ss(2)-Spe%P%centr(2))
             z_cdr =  (ss(3)-Spe%P%centr(3))
              r_cdr = sqrt(x_cdr**2+y_cdr**2+z_cdr**2)
              x_cdr = x_cdr/r_cdr
              y_cdr = y_cdr/r_cdr
              z_cdr = z_cdr/r_cdr             
             ! x_cdr,y_cdr and z_cdr are the grid point position in MSO frame
             ! we determine the latitude and longitude in the GEO (GCM) frame
             Lat_GEO = asin(sin(dec(nTime))*x_cdr+sin(clock)*cos(dec(nTime))*y_cdr + &
                & cos(clock)*cos(dec(nTime))*z_cdr)
             Xpc = cos(dec(nTime))*cos(cs0)*x_cdr + &
                & (-sin(clock)*sin(dec(nTime))*cos(cs0)-cos(clock)*sin(cs0))*y_cdr + &
                & (-cos(clock)*sin(dec(nTime))*cos(cs0)+sin(clock)*sin(cs0))*z_cdr
             Ypc = cos(dec(nTime))*sin(cs0)*x_cdr + &
                & (-sin(clock)*sin(dec(nTime))*sin(cs0)+cos(clock)*cos(cs0))*y_cdr + &
                & (-cos(clock)*sin(dec(nTime))*sin(cs0)-sin(clock)*cos(cs0))*z_cdr
             Lon_GEO = atan(Ypc/Xpc)  
            if (Xpc.lt.0._dp) then
              Lon_GEO = Lon_GEO + pi
            else
               if (Ypc.lt.0._dp) Lon_GEO = Lon_GEO + 2._dp*pi
            endif
            if (Lon_GEO.gt.pi) Lon_GEO = Lon_GEO - 2._dp*pi
            if (abs(Lon_GEO).gt.pi) then
                write(6,'(3(1x,i3),a,4(1x,e12.5))')ii,jj,kk,' Lon_GEO = ',Lon_GEO,Xpc,Ypc,atan(Ypc/Xpc)
                stop
            endif
            ! we determine the index for the corresponding longitude in the GCM
            ilon = 1
            ilat = 1
            diff_Long = 3.*pi
            diff_Lat  = 3.*pi
            do ilon_GCM = 2,nLong
              if ((Longitude(ilon_GCM-1)<Lon_GEO) .and. (Longitude(ilon_GCM) >= Lon_GEO)) ilon = ilon_GCM
              !if (abs(Longitude(ilon_GCM)-Lon_GEO) < diff_Long) then
              !   diff_Long = abs(Longitude(ilon_GCM)-Lon_GEO)
              !   ilon = ilon_GCM
              ! endif
            enddo
            do ilat_GCM = 2,nLati
              
              if ((Latitude(ilat_GCM-1) > Lat_GEO) .and. (Latitude(ilat_GCM) <= Lat_GEO)) ilat = ilat_GCM
              !print *,'Lattitude GCM, LAT_GEO,ilat',Latitude(ilat_GCM),Lat_GEO,ilat
              !if (abs(Latitude(ilat_GCM)-Lat_GEO) < diff_Lat) then
              !	 diff_Lat = abs(Latitude(ilat_GCM)-Lat_GEO)
              !	 ilat = ilat_GCM
              ! endif
            enddo
            iHoru   = ilon + (ilat - 1)*nLong
            
!!!!!!!! Change by F. Leblanc 01/2015 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!          if ((radius <= Spe%P%radius+maxval(Altitude(iHoru,:,nTime))/Spe%ref%c_omegapi).and.(radius >= Spe%P%radius)) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            Alt_GEO = (radius-Spe%P%radius)*Spe%ref%c_omegapi
            ialt = 1
            diff_alt = 4000.
            do ialt_GCM = 1,nAlti
              if (abs(Altitude(iHoru,ialt_GCM,nTime)-Alt_GEO) < diff_alt) then
                diff_alt = abs(Altitude(iHoru,ialt_GCM,nTime)-Alt_GEO)
                ialt = ialt_GCM
              endif  
            enddo
            density_O(ii,jj,kk) = nO(iHoru,ialt,nTime)
            density_CO2(ii,jj,kk) = nCO2(iHoru,ialt,nTime)
         endif ! end on the altitude check      
     enddo
   enddo
 enddo

 write(6,'(a,2(1x,e12.5))')' Loading Atmosphere Min, Max CO2 = ',minval(density_CO2(:,:,:)),maxval(density_CO2(:,:,:))
 write(6,'(a,2(1x,e12.5))')' Loading Atmosphere Min, Max O   = ',minval(density_O(:,:,:)),maxval(density_O(:,:,:))

deallocate(Time,zls,sza,dsm,dec,Longitude,Latitude,Altitude,nCO2,nO)
!deallocate(nO,nCO2p,nO2p,Tn)
__WRT_DEBUG_OUT("load_atmosphere_mars_LMD")
 end subroutine load_atmosphere_mars_LMD
 
 
 !!=============================================================
 !!subroutine: external_atmosphere/load_atmosphere_mars_Michigan
 !!
 !! FUNCTION
 !!  Read and load the Martian atmosphere from the GCM Michigan model (Bougher et al, 2000)
 !!   
 subroutine load_atmosphere_mars_Michigan(Spe,ncm,gstep,s_min_loc,resistivity,density_O,density_CO2)
 use defs_variable,only : viscosity
 use defs_parametre,only :dt,atmospherename
 use defs_mpitype,only     : mpiinfo

 integer, intent(in) :: ncm(3)
 real(dp),intent(in) :: gstep(3),s_min_loc(3)
 real(dp),intent(inout) :: resistivity(:,:,:),density_O(:,:,:),density_CO2(:,:,:)
 type(species_type),intent(in) :: Spe 
 character(len=500) :: msg 
 
 integer :: ncid,StId
 integer :: iunit,nr,nt,np,n_alt,n_theta,n_phi
 integer :: ii,jj,kk,iii,jjj,kkk
 real(dp) :: radius,altitude_km,r_planet_km,inv_r_exo_km,nrm
 real(dp) :: ss(3),x_cdr,y_cdr,z_cdr,tan_p,tan_t
 real(dp),dimension(:),allocatable :: Latitude, Longitude,Altitude,Distance
 real(dp),dimension(:,:,:),allocatable :: T,n_O,n_CO2
 real(dp) :: diff_alt,diff_theta,diff_phi,g
 real(dp) :: H_scale_height_O,H_scale_height_CO2,z,alt_km
 real(dp) :: dO_0,dO_1,eO_0,eO_1,hO_0,hO_1
 logical :: find_alt,find_phi,find_theta

 
   __WRT_DEBUG_IN("load_atmosphere_mars_Michigan")
 
 ! add hot O corona
 ! for Solar max
! dO_0 =5.e2 
! dO_1 =5.e4 
! eO_0 =-2.48_dp 
! eO_1 =-1.86_dp 
! hO_0 =800._dp 
! hO_1 =310_dp 
 
 ! for Solar min
  dO_0 =5.e2 
  dO_1 =9.e4 
  eO_0 =-3.08_dp 
  eO_1 =-1.62_dp 
  hO_0 =100._dp 
  hO_1 =100_dp 
 
 
 
 call wrt_double(6,"Reading file : "//atmospherename,wrtscreen,wrtdisk)
 
 
 !-- Open NetCDF file
 stId = nf90_open(trim(atmospherename),nf90_nowrite, ncid)
 call test_cdf(stId)
 call get_simple_variable_cdf(ncid,"nr",nr)
 call get_simple_variable_cdf(ncid,"nt",nt)
 call get_simple_variable_cdf(ncid,"np",np)
 
 write(msg,'(2a,3(a,a20,i7))')&
        & ch10," ___________________ Atmosphere Michigan information  _________________",&
        & ch10," nr  = ",nr,&
        & ch10," nt  = ",nt,&
        & ch10," np  = ",np
 call wrt_double(qp_out,msg,wrtscreen,wrtdisk) 
 
 allocate(Altitude(nr));        Altitude(:) = zero
 allocate(Latitude(nt));        Latitude(:) = zero
 allocate(Longitude(np));       Longitude(:) = zero
 
 call get_simple_variable_cdf(ncid,"Altitude",Altitude)
 call get_simple_variable_cdf(ncid,"Latitude",Latitude)
 call get_simple_variable_cdf(ncid,"Longitude",Longitude)
 write(msg,'(2a,2(a,a20,e12.5))')&
        & ch10," ___________________ Atmosphere Altitude  _________________",&
        & ch10, " Min (km)       = ",minval(Altitude(:)),&
        & ch10, " Max            = ",maxval(Altitude(:))
 call wrt_double(qp_out,msg,wrtscreen,wrtdisk)
 write(msg,'(2a,2(a,a20,e12.5))')&
       & ch10," ___________________ Atmosphere Longitude  _________________",&
       & ch10, " Min (degrees)  = ",minval(Longitude(:)),&
       & ch10, " Max            = ",maxval(Longitude(:))
 call wrt_double(qp_out,msg,wrtscreen,wrtdisk)
 write(msg,'(2a,2(a,a20,e12.5))')&
       & ch10," ___________________ Atmosphere Latitude  _________________",&
       & ch10, " Min (degrees)  = ",minval(Latitude(:)),&
       & ch10, " Max            = ",maxval(Latitude(:))
 call wrt_double(qp_out,msg,wrtscreen,wrtdisk)
 
 allocate(T(nr,np,nt)); T(:,:,:) = zero
 allocate(n_O(nr,np,nt));       n_O(:,:,:) = zero
 allocate(n_CO2(nr,np,nt));     n_CO2(:,:,:) = zero
 call get_simple_variable_cdf(ncid,"Temperature",T)
 call get_simple_variable_cdf(ncid,"Density_O",n_O)
 call get_simple_variable_cdf(ncid,"Density_CO2",n_CO2)
 
 write(msg,'(2a,2(a,a20,e12.5))')&
        & ch10," ___________________ Atmosphere CO2  _________________",&
        & ch10, " Min (cm-3)  = ",minval(n_CO2(:,:,:)),&
        & ch10, " Max         = ",maxval(n_CO2(:,:,:))
 call wrt_double(qp_out,msg,wrtscreen,wrtdisk)

 write(msg,'(2a,2(a,a20,e12.5))')&
        & ch10," ___________________ Atmosphere O  _________________",&
        & ch10, " Min (cm-3)  = ",minval(n_O(:,:,:)),&
        & ch10, " Max         = ",maxval(n_O(:,:,:))
 call wrt_double(qp_out,msg,wrtscreen,wrtdisk)

 write(msg,'(2a,2(a,a20,e12.5))')&
        & ch10," ___________________ Atmosphere T  _________________",&
        & ch10, " Min (K)  = ",minval(T(:,:,:)),&
        & ch10, " Max      = ",maxval(T(:,:,:))
 call wrt_double(qp_out,msg,wrtscreen,wrtdisk)

 
 !--Close the file
 stId = nf90_close(ncid); call test_cdf(stId)
 
 allocate(Distance(nr));        Distance(:) =zero
 Distance = Altitude/(Spe%ref%c_omegapi)+Spe%P%radius
 Longitude = Longitude*pi/180._dp
 Latitude = Latitude*pi/180._dp

  
 write(msg,'(2a,2(a,a20,e12.5))')&
        & ch10," ___________________ Atmosphere Distance  _________________",&
        & ch10, " Min (c/wpi)       = ",minval(Distance(:)),&
        & ch10, " Max            = ",maxval(Distance(:))
 call wrt_double(qp_out,msg,wrtscreen,wrtdisk)
 write(msg,'(2a,2(a,a20,e12.5))')&
       & ch10," ___________________ Atmosphere Longitude  _________________",&
       & ch10, " Min (rad)  = ",minval(Longitude(:)),&
       & ch10, " Max            = ",maxval(Longitude(:))
 call wrt_double(qp_out,msg,wrtscreen,wrtdisk)
 write(msg,'(2a,2(a,a20,e12.5))')&
       & ch10," ___________________ Atmosphere Latitude  _________________",&
       & ch10, " Min (rad)  = ",minval(Latitude(:)),&
       & ch10, " Max            = ",maxval(Latitude(:))
 call wrt_double(qp_out,msg,wrtscreen,wrtdisk)
 
  density_O(:,:,:) = 1.e-10
  density_CO2(:,:,:) = 1.e-10
 
 !--Main loop
 do kk = 1,ncm(3)-1
   do jj = 1,ncm(2)-1
     do ii = 1,ncm(1)-1
!        do kkk=-4,4
!        do jjj=-4,4
!        do iii=-4,4
!           ss(3) = (real((kk-1),dp)+real(kkk,dp)*0.1111111)*gstep(3) + s_min_loc(3) 
!           ss(2) = (real((jj-1),dp)+real(jjj,dp)*0.1111111)*gstep(2) + s_min_loc(2)
!           ss(1) = (real((ii-1),dp)+real(iii,dp)*0.1111111)*gstep(1) + s_min_loc(1)      
           ss(3) = (real((kk-1),dp))*gstep(3) + s_min_loc(3) 
           ss(2) = (real((jj-1),dp))*gstep(2) + s_min_loc(2)
           ss(1) = (real((ii-1),dp))*gstep(1) + s_min_loc(1)      

           radius = sqrt(dot_product(ss-Spe%P%centr,ss-Spe%P%centr))
!          if (radius <= Distance(nr)+1000.*Spe%P%radius/(Spe%ref%c_omegapi)) then
          if ((radius <= Distance(nr)+1000./Spe%ref%c_omegapi).and.(radius >= Distance(1))) then
           x_cdr = -(ss(1)-Spe%P%centr(1))
           y_cdr = -(ss(2)-Spe%P%centr(2))
           z_cdr =  (ss(3)-Spe%P%centr(3))
           tan_p = atan(y_cdr/x_cdr)
           tan_t = atan(sqrt(z_cdr**2/(x_cdr**2+y_cdr**2)))
 ! check conditions to obtain the good angles
 !          if (((x_cdr == 0.).and.(y_cdr == 0.)).and.(z_cdr > 0.))  tan_t = pi/2._dp
 !          if (((x_cdr == 0.).and.(y_cdr == 0.)).and.(z_cdr < 0.))  tan_t = -pi/2._dp
 !          if (z_cdr < 0.) tan_t = -tan_t
 !          
 !          !if ((x_cdr > 0.).and. (y_cdr == 0.)) tan_p = 0._dp
 !          
 !          if ((x_cdr < 0.).and. (y_cdr < 0.)) tan_p = pi/2._dp+tan_p
 !          if ((x_cdr < 0.).and. (y_cdr == 0.)) tan_p = pi
 !          if ((x_cdr < 0.).and. (y_cdr >0.)) tan_p = pi+tan_p
 !          
 !          if ((x_cdr == 0.).and.(y_cdr > 0.)) tan_p = pi/2._dp
 !          if ((x_cdr == 0.).and.(y_cdr < 0.)) tan_p = -pi/2._dp
 !          if ((x_cdr == 0.).and.(y_cdr == 0.)) tan_p = 0._dp
 ! check conditions to obtain the good angles
            if (((x_cdr == 0.).and.(y_cdr == 0.)).and.(z_cdr > 0.))  tan_t = pi/2._dp
            if (((x_cdr == 0.).and.(y_cdr == 0.)).and.(z_cdr < 0.))  tan_t = -pi/2._dp
            if (z_cdr < 0.) tan_t = -tan_t
            
            !if ((x_cdr > 0.).and. (y_cdr == 0.)) tan_p = 0._dp
            
            if ((x_cdr < 0.).and. (y_cdr < 0.)) tan_p = -pi+tan_p
            if ((x_cdr < 0.).and. (y_cdr == 0.)) tan_p = pi
            if ((x_cdr < 0.).and. (y_cdr >0.)) tan_p = pi+tan_p
            
            if ((x_cdr == 0.).and.(y_cdr > 0.)) tan_p = pi/2._dp
            if ((x_cdr == 0.).and.(y_cdr < 0.)) tan_p = -pi/2._dp
           if ((x_cdr == 0.).and.(y_cdr == 0.)) tan_p = 0._dp
 !          print *,'Coming IN'
 !!! check up to there
           
           n_alt = 2
           n_phi = 2
           n_theta = 2
           diff_alt = radius-Distance(1)
           diff_phi = abs(tan_p-Longitude(1))
           diff_theta = abs(tan_t - Latitude(1))
           find_alt = .false.
           find_theta = .false.
           find_phi = .false.

           ! Find the good phi
           do while ((find_phi .eqv. .false.).and.(n_phi <= np))
             if (abs(tan_p-Longitude(n_phi)) < abs(diff_phi)) then
               diff_phi = abs(tan_p-Longitude(n_phi))
               n_phi = n_phi+1
             else
               find_phi =.true.
             endif
           enddo
           ! Find the good theta
           do while ((find_theta .eqv. .false.).and.(n_theta <= nt))
             if (abs(tan_t-Latitude(n_theta)) < abs(diff_theta)) then
               diff_theta = abs(tan_t-Latitude(n_theta))
               n_theta = n_theta+1
             else
               find_theta =.true.
             endif
           enddo
          ! if (((find_alt==.true.).and.(find_phi==.true.)).and.(find_theta == .true.)) then
          if (find_phi .eqv. .false.) n_phi = 1
          if (find_theta .eqv. .false.) n_theta = 1
          
          
          if (radius <=Distance(nr)) then
          ! Find the good altitude
             do while ((find_alt .eqv. .false.).and.(n_alt <= nr))
               if (abs(radius-Distance(n_alt)) < abs(diff_alt)) then
                   diff_alt = radius-Distance(n_alt)
                   n_alt = n_alt+1
               else
                   find_alt =.true.
               endif
             enddo
          
          
          if (find_alt .eqv. .true.)   then
             density_O(ii,jj,kk) = density_O(ii,jj,kk)+ n_O(n_alt,n_phi,n_theta)
             density_CO2(ii,jj,kk) = density_CO2(ii,jj,kk)+ n_CO2(n_alt,n_phi,n_theta)         
            ! print *,'Load atm',density_O(ii,jj,kk), density_CO2(ii,jj,kk)
          endif   
          endif 
          if (radius > Distance(nr)) then 
            g = M_Mars*G_grav/((Distance(nr)*Spe%ref%c_omegapi*1.e3)**2)
            H_scale_height_O = (kb_JK*T(nr,n_phi,n_theta))/(g*16.*amu_pmass)
            H_scale_height_CO2 = (kb_JK*T(nr,n_phi,n_theta))/(g*44.*amu_pmass)
            
            z = (radius-Distance(nr))*Spe%ref%c_omegapi*1.e3
             
            density_O(ii,jj,kk) = n_O(nr,n_phi,n_theta)*exp(-z/H_scale_height_O)
            density_CO2(ii,jj,kk) = n_CO2(nr,n_phi,n_theta)*exp(-z/H_scale_height_CO2)
            

            
          endif

           
 !          if ((radius < altitude(1))) then
 !            g = M_Mars*G_grav/((altitude(2)*Spe%ref%c_omegapi*1.e3)**2)
 !            
 !            H_scale_height_O = (kb_JK*tempO(2,n_theta,n_phi))/(g*16.*amu_pmass)
 !            H_scale_height_CO2 = (kb_JK*tempCO2(2,n_theta,n_phi))/(g*44.*amu_pmass)
 !            
 !            n0_CO2 = densCO2(2,n_theta,n_phi)*exp((altitude(2)*Spe%ref%c_omegapi-rad_mars)*1.e3/H_scale_height_CO2)
 !            n0_O = densO(2,n_theta,n_phi)*exp((altitude(2)*Spe%ref%c_omegapi-rad_mars)*1.e3/H_scale_height_O)
 !            
 !            density_O(ii,jj,kk) = n0_O*exp(-(radius*Spe%ref%c_omegapi-rad_mars)*1.e3/H_scale_height_O)
 !            density_CO2(ii,jj,kk) = n0_CO2*exp(-(radius*Spe%ref%c_omegapi-rad_mars)*1.e3/H_scale_height_CO2)
 !
 !                   
 !          endif  
           if (radius < Spe%P%radius) density_O(ii,jj,kk) = 1.e-10
           if (radius < Spe%P%radius) density_CO2(ii,jj,kk) = 1.e-10
             
                   
           !if ((radius >= altitude(1)).and.(radius <= altitude(nr))) then
            ! n_alt = minloc(abs(altitude - radius),dim=1)
            ! n_phi = minloc(abs(tan_p-tan(phi)),dim=1)
            ! n_theta = minloc(abs(tan_t-tan(theta)),dim=1)
             !print *,'Altitude :',n_alt,n_theta,n_phi,dens(n_alt,n_theta,n_phi)
            ! density_O(ii,jj,kk) = dens(n_alt,n_theta,n_phi)
           !endif
           endif   !on the altitude
            if (radius > Distance(nr)+200./Spe%ref%c_omegapi) then
                        alt_km = (radius-Spe%P%radius)*Spe%ref%c_omegapi
                        density_O(ii,jj,kk) = density_O(ii,jj,kk) + &
                              & dO_0*(alt_km/hO_0)**eO_0 + dO_1*(alt_km/hO_1)**eO_1
            endif           
!           enddo
!      	  enddo
!      	  enddo
          enddo
        enddo
     enddo
  
  write(msg,'(2a,2(a,a20,e12.5))')&
         & ch10," ___________________ Atmosphere CO2 loaded  _________________",&
         & ch10, " Min (cm-3)  = ",minval(density_CO2(:,:,:)),&
         & ch10, " Max         = ",maxval(density_CO2(:,:,:))
  call wrt_double(qp_out,msg,wrtscreen,wrtdisk)
 
  write(msg,'(2a,2(a,a20,e12.5))')&
         & ch10," ___________________ Atmosphere O loaded _________________",&
         & ch10, " Min (cm-3)  = ",minval(density_O(:,:,:)),&
         & ch10, " Max         = ",maxval(density_O(:,:,:))
 call wrt_double(qp_out,msg,wrtscreen,wrtdisk)
 
 
 deallocate(Altitude,Latitude,Longitude)
 deallocate(T,n_O,n_CO2)
   __WRT_DEBUG_OUT("load_atmosphere_mars_Michigan")
  
 end subroutine load_atmosphere_mars_Michigan

end module atm_external_atmosphere
