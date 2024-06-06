!!=============================================================
!!=============================================================
!!module: atm_magnetic_fields
!! NAME
!!  atm_magnetic_fields (SHess)
!!
!! Gathers all the planet magnetic field subroutines
!!
module atm_magnetic_fields

 use defs_basis
 use defs_species
 use defs_parametre
 use defs_arr3Dtype
 use m_writeout

#include "q-p_common.h"
 implicit none
 private

 public ::                   &
        add_dipole_generic,  &
        add_dipquad,        &
        preserve_dip_generic, &
        preserve_dipquad,    &
        add_multipole_generic
 private ::&
        convert,&
        calc_legendre

contains

 !!=============================================================
 !!routine: atm_magnetic_fields/add_dipole_generic
 !! NAME
 !!  add_dipole_generic (RModolo,GM Chanteur, E. Richer, S. Hess, MMancini,RAllioux)
 !!
 !! FUNCTION
 !!  Contains Dipole moment calculation generic
 !!
 subroutine add_dipole_generic(Bfield,ncm,Spe,gstep,s_min_loc,bme,inclinaison,phase)
 use defs_arr3Dtype

  integer, intent(in) :: ncm(3)
  type(arr3Dtype),intent(inout) :: Bfield
  real(dp),intent(in) :: gstep(3),s_min_loc(3),bme,inclinaison,phase
  type(species_type),intent(in) :: Spe

  !local
  integer :: ii,jj,kk
  real(dp) :: radius_dip,radius_pl,rinclinaison,rphase
  real(dp),dimension(3) :: ss,moment_dip,b_dip,moment_dip_u

  __WRT_DEBUG_IN("add_dipole_generic")  
  rinclinaison = inclinaison*deg_to_rad
  rphase = phase*deg_to_rad

  !--Initialisation
  !--Dipolar Moment in planetary units (Spe%ref%mag) mu0/4pi*M
   moment_dip_u(1) = sin(rinclinaison)*cos(rphase)
   moment_dip_u(2) = sin(rinclinaison)*sin(rphase)
   moment_dip_u(3) = cos(rinclinaison)
   moment_dip = moment_dip_u*bme/Spe%ref%mag*(Spe%P%radius)**three
  !--Main loop
  !--Relative distance from the dipole position in m
  do kk = 1,ncm(3)-1
   ss(3) = (real((kk-1),dp)*gstep(3) + s_min_loc(3)-Spe%P%dipole(3))          !Spe%P%centr(3))
   do jj = 1,ncm(2)-1
    ss(2) = (real((jj-1),dp)*gstep(2) + s_min_loc(2)-Spe%P%dipole(2))         !Spe%P%centr(2))
    do ii = 1,ncm(1)-1
     ss(1) = (real((ii-1),dp)*gstep(1) + s_min_loc(1)-Spe%P%dipole(1))        !Spe%P%centr(1))
     !--Distance from the center in radius of planet
     radius_dip = sqrt(dot_product(ss,ss))
         radius_pl = sqrt(dot_product(ss+Spe%P%dipole-Spe%P%centr,ss+Spe%P%dipole-Spe%P%centr))
     if (radius_pl >= .75*Spe%P%radius) then
      b_dip = (three*ss*dot_product(ss,moment_dip)/(radius_dip*radius_dip)-moment_dip)/radius_dip**three
      Bfield%x(ii,jj,kk) = Bfield%x(ii,jj,kk) + b_dip(1)
      Bfield%y(ii,jj,kk) = Bfield%y(ii,jj,kk) + b_dip(2)
      Bfield%z(ii,jj,kk) = Bfield%z(ii,jj,kk) + b_dip(3)
     else
      Bfield%x(ii,jj,kk) = zero
      Bfield%y(ii,jj,kk) = zero
      Bfield%z(ii,jj,kk) = zero
     endif
    enddo
   enddo
  enddo

  __WRT_DEBUG_OUT("add_dipole_generic")
 end subroutine add_dipole_generic
 
  !!=============================================================
 !!routine: atm_magnetic_fields/add_dipquad
 !! NAME
 !!  add_dipquad_generic (RModolo,GM Chanteur, E. Richer, S. Hess, MMancini,RAllioux)
 !!
 !! FUNCTION
 !!  Contains Dipole and Quadrupole moment calculation generic
 !!
 subroutine add_dipquad(Bfield,ncm,Spe,gstep,s_min_loc,bme,inclinaison,phase,cd,cq)
 use defs_arr3Dtype

  integer, intent(in) :: ncm(3)
  type(arr3Dtype),intent(inout) :: Bfield
  real(dp),intent(in) :: gstep(3),s_min_loc(3),bme,inclinaison,phase,cd,cq
  type(species_type),intent(in) :: Spe

  !local
  integer :: ii,jj,kk
  real(dp) :: radius,rinclinaison,rphase
  real(dp),dimension(3) :: ss,b_quad,b_2,moment_dip,b_dip,moment_dip_u

  __WRT_DEBUG_IN("add_dipquad")  
  rinclinaison = inclinaison*deg_to_rad
  rphase = phase*deg_to_rad

  !--Initialisation
  !--Dipolar Moment in planetary units (Spe%ref%mag) mu0/4pi*M
   moment_dip_u(1) = sin(rinclinaison)*cos(rphase)
   moment_dip_u(2) = sin(rinclinaison)*sin(rphase)
   moment_dip_u(3) = cos(rinclinaison)
   moment_dip = moment_dip_u*bme/Spe%ref%mag*(Spe%P%radius)**three

  !--Main loop
  !--Relative distance from the planet centre in m
  do kk = 1,ncm(3)-1
   ss(3) = (real((kk-1),dp)*gstep(3) + s_min_loc(3)-Spe%P%centr(3))          !Spe%P%centr(3))
   do jj = 1,ncm(2)-1
    ss(2) = (real((jj-1),dp)*gstep(2) + s_min_loc(2)-Spe%P%centr(2))         !Spe%P%centr(2))
    do ii = 1,ncm(1)-1
     ss(1) = (real((ii-1),dp)*gstep(1) + s_min_loc(1)-Spe%P%centr(1))        !Spe%P%centr(1))
     !--Distance from the center in radius of planet
     radius = sqrt(dot_product(ss,ss))
     if (radius >= .75*Spe%P%radius) then
          b_dip = (three*ss*dot_product(ss,moment_dip)/(radius*radius)-moment_dip)/radius**three
          
      b_quad(1) = -1*(-5.*radius**(-7)*(ss(1)**2+ss(2)**2-2.*ss(3)**2)*ss(1)+2.*radius**(-5)*ss(1))
          b_quad(2) = -1*(-5.*radius**(-7)*(ss(1)**2+ss(2)**2-2.*ss(3)**2)*ss(2)+2.*radius**(-5)*ss(2))
          b_quad(3) = -1*(-5.*radius**(-7)*(ss(1)**2+ss(2)**2-2.*ss(3)**2)*ss(3)-4.*radius**(-5)*ss(3))
          b_quad = b_quad*bme/Spe%ref%mag*(Spe%P%radius)**4
          
          b_2 = b_quad*cq+b_dip*cd
          
      Bfield%x(ii,jj,kk) = Bfield%x(ii,jj,kk) + b_2(1)
      Bfield%y(ii,jj,kk) = Bfield%y(ii,jj,kk) + b_2(2)
      Bfield%z(ii,jj,kk) = Bfield%z(ii,jj,kk) + b_2(3)
     else
      Bfield%x(ii,jj,kk) = zero
      Bfield%y(ii,jj,kk) = zero
      Bfield%z(ii,jj,kk) = zero
     endif
    enddo
   enddo
  enddo

  __WRT_DEBUG_OUT("add_dipquad")
 end subroutine add_dipquad
 
 !!=============================================================
 !!routine: atm_magnetic_fields/preserve_dip_generic
 !! NAME
 !!  preserve_dip_generic (E. Richer)
 !!
 !! FUNCTION
 !!  Contains Dipole moment calculation generic 
 !!
 subroutine preserve_dip_generic(Bfield,ncm,Spe,gstep,s_min_loc,bme,inclinaison,phase,percentr)
 use defs_arr3Dtype

  integer, intent(in) :: ncm(3)
  type(arr3Dtype),intent(inout) :: Bfield
  real(dp),intent(in) :: gstep(3),s_min_loc(3),bme,inclinaison,phase
  real(dp),intent(in) :: percentr
  type(species_type),intent(in) :: Spe

  !local
  integer :: ii,jj,kk
  real(dp) :: radius_dip,radius_pl,rinclinaison,rphase,rss,fracr
  real(dp),dimension(3) :: ss,moment_dip,b_dip,moment_dip_u

  __WRT_DEBUG_IN("preserve_dip_generic")  
  rinclinaison = inclinaison*deg_to_rad
  rphase = phase*deg_to_rad
  fracr = percentr/100.
  !--Initialisation
  !--Dipolar Moment in planetary units (Spe%ref%mag) mu0/4pi*M
   moment_dip_u(1) = sin(rinclinaison)*cos(rphase)
   moment_dip_u(2) = sin(rinclinaison)*sin(rphase)
   moment_dip_u(3) = cos(rinclinaison)
   moment_dip = moment_dip_u*bme/Spe%ref%mag*(Spe%P%radius)**three
  !--Main loop
  !--Relative distance from the dipole position in m
  do kk = 1,ncm(3)-1
   ss(3) = (real((kk-1),dp)*gstep(3) + s_min_loc(3)- Spe%P%dipole(3))             !Spe%P%centr(3))
   if (ss(3)>=Spe%P%radius) CYCLE
   do jj = 1,ncm(2)-1
    ss(2) = (real((jj-1),dp)*gstep(2) + s_min_loc(2)-Spe%P%dipole(2))             !Spe%P%centr(2))
        !rss = sqrt(ss(3)**2+ss(2)**2)
        if (ss(2)>=Spe%P%radius) CYCLE
    do ii = 1,ncm(1)-1
     ss(1) = (real((ii-1),dp)*gstep(1) + s_min_loc(1)-Spe%P%dipole(1))            !Spe%P%centr(1))
     !--Distance from the center in radius of planet
     radius_dip = sqrt(dot_product(ss,ss))
         radius_pl = sqrt(dot_product(ss+Spe%P%dipole-Spe%P%centr,ss+Spe%P%dipole-Spe%P%centr))
     if ((radius_pl >= fracr*Spe%P%radius).and.(radius_pl<Spe%P%radius)) then
      b_dip = (three*ss*dot_product(ss,moment_dip)/(radius_dip*radius_dip)-moment_dip)/radius_dip**three
      Bfield%x(ii,jj,kk) = b_dip(1)
      Bfield%y(ii,jj,kk) = b_dip(2)
      Bfield%z(ii,jj,kk) = b_dip(3)
     endif
         if (radius_pl < fracr*Spe%P%radius) then
      Bfield%x(ii,jj,kk) = zero
      Bfield%y(ii,jj,kk) = zero
      Bfield%z(ii,jj,kk) = zero
     endif
    enddo
   enddo
  enddo

  __WRT_DEBUG_OUT("preserve_dip_generic")
 end subroutine preserve_dip_generic
 
   !!=============================================================
 !!routine: atm_magnetic_fields/preserve_dipquad
 !! NAME
 !!  add_dipquad_generic (RModolo,GM Chanteur, E. Richer, S. Hess, MMancini,RAllioux)
 !!
 !! FUNCTION
 !!  Contains Dipole and Quadrupole moment calculation generic
 !!
 subroutine preserve_dipquad(Bfield,ncm,Spe,gstep,s_min_loc,bme,inclinaison,phase,percentr,cd,cq)
 use defs_arr3Dtype

  integer, intent(in) :: ncm(3)
  type(arr3Dtype),intent(inout) :: Bfield
  real(dp),intent(in) :: gstep(3),s_min_loc(3),bme,inclinaison,phase,cd,cq
  type(species_type),intent(in) :: Spe

  !local
  integer :: ii,jj,kk
  real(dp) :: radius,rinclinaison,rphase,fracr
  real(dp),intent(in) :: percentr
  real(dp),dimension(3) :: ss,b_quad,b_2,moment_dip,b_dip,moment_dip_u

  __WRT_DEBUG_IN("preserve_dipquad")  
  rinclinaison = inclinaison*deg_to_rad
  rphase = phase*deg_to_rad
  
  fracr = percentr/100.

  !--Initialisation
  !--Dipolar Moment in planetary units (Spe%ref%mag) mu0/4pi*M
   moment_dip_u(1) = sin(rinclinaison)*cos(rphase)
   moment_dip_u(2) = sin(rinclinaison)*sin(rphase)
   moment_dip_u(3) = cos(rinclinaison)
   moment_dip = moment_dip_u*bme/Spe%ref%mag*(Spe%P%radius)**three

  !--Main loop
  !--Relative distance from the planet centre in m
  do kk = 1,ncm(3)-1
   ss(3) = (real((kk-1),dp)*gstep(3) + s_min_loc(3)-Spe%P%centr(3))          !Spe%P%centr(3))
   do jj = 1,ncm(2)-1
    ss(2) = (real((jj-1),dp)*gstep(2) + s_min_loc(2)-Spe%P%centr(2))         !Spe%P%centr(2))
    do ii = 1,ncm(1)-1
     ss(1) = (real((ii-1),dp)*gstep(1) + s_min_loc(1)-Spe%P%centr(1))        !Spe%P%centr(1))
     !--Distance from the center in radius of planet
     radius = sqrt(dot_product(ss,ss))
     if ((radius >= fracr*Spe%P%radius).and.(radius<Spe%P%radius)) then
          b_dip = (three*ss*dot_product(ss,moment_dip)/(radius*radius)-moment_dip)/radius**three
          
      b_quad(1) = -1*(-5.*radius**(-7)*(ss(1)**2+ss(2)**2-2.*ss(3)**2)*ss(1)+2.*radius**(-5)*ss(1))
          b_quad(2) = -1*(-5.*radius**(-7)*(ss(1)**2+ss(2)**2-2.*ss(3)**2)*ss(2)+2.*radius**(-5)*ss(2))
          b_quad(3) = -1*(-5.*radius**(-7)*(ss(1)**2+ss(2)**2-2.*ss(3)**2)*ss(3)-4.*radius**(-5)*ss(3))
          b_quad = b_quad*bme/Spe%ref%mag*(Spe%P%radius)**4
          
          b_2 = b_dip*cd+cq*b_quad
          
      Bfield%x(ii,jj,kk) = b_2(1)
      Bfield%y(ii,jj,kk) = b_2(2)
      Bfield%z(ii,jj,kk) = b_2(3)
         endif
     if (radius < fracr*Spe%P%radius) then
      Bfield%x(ii,jj,kk) = zero
      Bfield%y(ii,jj,kk) = zero
      Bfield%z(ii,jj,kk) = zero
     endif
    enddo
   enddo
  enddo

  __WRT_DEBUG_OUT("preserve_dipquad")
 end subroutine preserve_dipquad


 !!=============================================================
 !!routine: atm_magnetic_fields/calc_legendre
 !! NAME
 !!  calc_legendre (S. Hess)
 !!
 !! tableau de polynome de legendre Yml et des derivees en theta
 subroutine calc_legendre(n_sph,l,p,d_p,theta)
 integer,intent(in) :: l,n_sph
 real(dp),intent(inout),dimension(0:n_sph,0:n_sph) ::p,d_p
 real(dp),intent(in) ::theta
 real(dp) ::somx2,somxpdx2,xpdx,x,Y2P_coeff
 integer ::i,j,k,n
  x=theta
  xpdx=x+0.001_dp
  x=cos(x)
  xpdx=cos(xpdx)

  somx2=sqrt((1._dp-x)*(1._dp+x))
  somxpdx2=sqrt((1._dp-xpdx)*(1._dp+xpdx))

  p(0:l,0)=sqrt(2.); d_p(0:l,0)=sqrt(2.)
  do k=1,l
        do i=k,l
        p(i,k)   = p(i,k-1)*(-(1._dp+2._dp*real(k-1,dp))*somx2)*(-1._dp)/sqrt(real(i-k+1,dp)*real(i+k,dp))
        d_p(i,k) = d_p(i,k-1)*(-(1._dp+2._dp*real(k-1,dp))*somxpdx2)*(-1._dp)/sqrt(real(i-k+1,dp)*real(i+k,dp))
        enddo
  enddo
  p(0:l,0)=1.; d_p(0:l,0)=1.
  do  k=0,l-1 
        p(k+1,k)  = x*(2._dp*real(k,dp)+1._dp)*p(k+1,k)
        d_p(k+1,k) = xpdx*(2._dp*real(k,dp)+1._dp)*d_p(k+1,k)
  enddo
  if (l.ge.2) then 
  do n=2,l 
        do k=0,l-n 
                p(k+n,k)=(x*(2._dp*real(k+n,dp)-1._dp)*p(k+n-1,k)/&
        & sqrt(2._dp*real(k,dp)/real(n,dp)+1._dp)-&
        & (real(2*k+n,dp)-1._dp)*p(k+n-2,k)/sqrt(2._dp*real(k,dp)/&
        & real(n,dp)+1._dp)/sqrt(2._dp*real(k,dp)/real(n-1,dp)+1._dp))/&
        & real(n,dp)
                d_p(k+n,k)=(xpdx*(2._dp*real(k+n,dp)-1._dp)*d_p(k+n-1,k)/&
        & sqrt(2._dp*real(k,dp)/real(n,dp)+1._dp)-&
        & (real(2*k+n,dp)-1._dp)*d_p(k+n-2,k)/sqrt(2._dp*real(k,dp)/&
        & real(n,dp)+1._dp)/sqrt(2._dp*real(k,dp)/real(n-1,dp)+1._dp))/&
        & real(n,dp)
        enddo
  enddo
  endif
  d_p=(d_p-p)*1E3
 end subroutine calc_legendre

 !!=============================================================
 !!routine: atm_magnetic_fields/add_multipole_generic
 !! NAME
 !!  add_multipole_generic (S. Hess)
 !!
 !! FUNCTION
 !!  Contains Multipole moment calculation generic
 !!
subroutine add_multipole_generic(Bfield,ncm,Spe,gstep,s_min_loc,n_sph,g,h)


  integer, intent(in) :: ncm(3),n_sph
  type(arr3Dtype),intent(inout) :: Bfield
  real(dp),intent(in) :: gstep(3),s_min_loc(3)
  real(dp),intent(in),dimension(1:n_sph,0:n_sph) ::g,h
  type(species_type),intent(in) :: Spe
   
  !local
  real(dp),dimension(0:n_sph,0:n_sph) ::p,d_p
  real(dp),dimension(0:n_sph) :: Co,Si
  real(dp) :: Obliqui,cs0,dec,clock,x_cdr,y_cdr,z_cdr,&
        & r_cdr,Lon_GEO,Lat_GEO,cosclock,Xpc,Ypc
  integer :: ii,jj,kk,i,j,lmax
  real(dp) :: radius,ss(3),b_sph(3),phi,theta,rp,b_xyz(3),satur,thetai,phii

  __WRT_DEBUG_IN("add_multipole_generic")
  Bfield%x(:,:,:)=0.;   Bfield%y(:,:,:)=0.;   Bfield%z(:,:,:)=0.;  Co(0)=1._dp ; Si(0)=0._dp;
   ! specificity for Mars
   Obliqui = 0.43
   cs0 = Spe%P%ssl*deg_to_rad
   dec = Spe%P%sslat*deg_to_rad
   cosclock = cos(Obliqui)/cos(dec)
   if (abs(cosclock)<=1) then
      clock = acos(cos(Obliqui)/cos(dec))
   else
     if (cosclock > 1) clock=0.
     if (cosclock<-1) clock = pi
   endif
  print *,'B CF: clock,dec,cs0',clock,dec,cs0

  !--Main loop
  do kk = 1,ncm(3)-1
   ss(3) = (real((kk-1),dp)*gstep(3) + s_min_loc(3)-Spe%P%centr(3))
   do jj = 1,ncm(2)-1
    ss(2) = (real((jj-1),dp)*gstep(2) + s_min_loc(2)-Spe%P%centr(2))
    do ii = 1,ncm(1)-1
     ss(1) = (real((ii-1),dp)*gstep(1) + s_min_loc(1)-Spe%P%centr(1))
     !--Distance from the center in radius of planet
     b_sph(:)=0.
     radius = sqrt(dot_product(ss,ss))
     lmax=min(int(100.*(Spe%P%radius/radius)**3),&
        & int(100.*(radius/Spe%P%radius)**20)) ! limits the precision as a function of distance
     lmax=max(lmax,1);  lmax=min(lmax,n_sph)
     if ((radius >= .85*Spe%P%radius).and.(lmax.gt.0)) then
            phii=atan2(ss(2),ss(1))
            thetai=atan2(sqrt(ss(1)**2+ss(2)**2),ss(3))
!  This parameter should be specifed only for Mars !
           if (Spe%P%ssl >-1) then
! conversion to MSO (normalized
              r_cdr = sqrt(ss(1)**2+ss(2)**2+ss(3)**2)
              x_cdr = -ss(1)/r_cdr
              y_cdr = -ss(2)/r_cdr
              z_cdr = ss(3)/r_cdr             
             ! x_cdr,y_cdr and z_cdr are the grid point position in MSO frame
             ! we determine the latitude and longitude in the GEO (GCM) frame
             Lat_GEO = asin(sin(dec)*x_cdr+sin(clock)*cos(dec)*y_cdr + &
                & cos(clock)*cos(dec)*z_cdr)
             Xpc = cos(dec)*cos(cs0)*x_cdr + &
                & (-sin(clock)*sin(dec)*cos(cs0)-cos(clock)*sin(cs0))*y_cdr + &
                & (-cos(clock)*sin(dec)*cos(cs0)+sin(clock)*sin(cs0))*z_cdr
             Ypc = cos(dec)*sin(cs0)*x_cdr + &
                & (-sin(clock)*sin(dec)*sin(cs0)+cos(clock)*cos(cs0))*y_cdr + &
                & (-cos(clock)*sin(dec)*sin(cs0)-sin(clock)*cos(cs0))*z_cdr
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
            phi = Lon_GEO; theta = pi/2-Lat_GEO;
           else
              phi = phii; theta = thetai;
           endif
            call calc_legendre(n_sph,lmax,p,d_p,theta); p=p/Spe%ref%mag; d_p=d_p/Spe%ref%mag
        do i=1,lmax 
            Co(i)=cos(real(i,dp)*phi);   Si(i)=sin(real(i,dp)*phi);   RP=radius/Spe%P%radius; if (RP.lt.1._dp) RP=1._dp 
                !limits the magnetic field inside the planet
            RP=(1._dp/RP)**(i+2)
                do j=0,i
                        b_sph(1)=b_sph(1)+real(i+1,dp)*(P(I,J)*(G(i,j)*Co(j)+H(i,j)*Si(j)))*RP
                        b_sph(2)=b_sph(2)-(D_P(I,J)*(G(i,j)*Co(j)+H(i,j)*Si(j)))*RP
                        b_sph(3)=b_sph(3)+(P(I,J)*real(j,dp)*(G(i,j)*Si(j)-H(i,j)*Co(j)))*RP
                enddo
        enddo
        if (theta.ne.0) then 
        b_sph(3)=b_sph(3)/sin(theta) 
        else 
        b_sph(3)=0.
        endif
        b_xyz=convert(b_sph,thetai,phii)
        Bfield%x(ii,jj,kk)=Bfield%x(ii,jj,kk)+b_xyz(1)
        Bfield%y(ii,jj,kk)=Bfield%y(ii,jj,kk)+b_xyz(2)
        Bfield%z(ii,jj,kk)=Bfield%z(ii,jj,kk)+b_xyz(3)
     endif
    enddo
   enddo
  enddo

  __WRT_DEBUG_OUT("add_multipole_generic")
 end subroutine add_multipole_generic

function convert(brtp,theta,phi)
real(dp),dimension(1:3),intent(in):: brtp
real(dp),intent(in)::phi,theta
real(dp),dimension(1:3)::convert
convert(3)=cos(theta)*brtp(1)-sin(theta)*brtp(2)
convert(1)=sin(theta)*cos(phi)*brtp(1)+cos(theta)*cos(phi)*brtp(2)-sin(phi)*brtp(3)
convert(2)=sin(theta)*sin(phi)*brtp(1)+cos(theta)*sin(phi)*brtp(2)+cos(phi)*brtp(3)
end function

end module atm_magnetic_fields

