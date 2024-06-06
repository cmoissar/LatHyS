!!=============================================================
!!=============================================================
module part_moment

 use defs_basis
 use defs_arr3Dtype
 use defs_particletype
 use m_writeout,only       : wrt_debug
 use m_timing,only         : time_get
 use field_cond_limit,only       : cond_limit_func  
#include "q-p_common.h"

 implicit none
 private

 private::          &
      momad3r,      &
      Cmtsp3       

 public::           &
      momtin,       &
      momt3d,       &
      momad,        &
      Amtsp3,       &
      Bmtsp3,       &
      Dmtsp3,       &
      Emtsp3,       &
      Fmtsp3  

contains
 !!#####################################################################

 !********************************* MOMTIN ************************************
 subroutine momtin(particule,dn,dna,vel,vela,gstep,s_min_loc,iwr,nptot)
  !--Calculation of moments in 3-D array.
  ! Charge and currents in DNA and UA.
  ! MASS density DN ; Coherent velocity U

  use defs_grid
  use defs_mpitype,only     : mpiinfo
  use field_lissage,only          : smth_func

  type(particletype),intent(in) :: particule(:)
  real(dp),intent(inout),dimension(:,:,:) :: dn,dna
  type(arr3Dtype),intent(inout) :: vel,vela
  real(dp),intent(in) :: gstep(3),s_min_loc(3)
  integer,intent(in) :: iwr(:)
  integer,intent(in) :: nptot

  __WRT_DEBUG_IN("momtin")

  !--Initialise moment arrays to 0 now done in Cmtsp3

  !--Collect species moments in usa
  call Cmtsp3 (particule,dna,vela,dn,vel,1,nptot,gstep,s_min_loc,&
       &       mpiinfo,nc1)

  !--Smoothing
  if(iwr(6)==1)then
   call smth_func(dna,   nc,mpiinfo)
   call smth_func(vela%x,nc,mpiinfo)
   call smth_func(vela%y,nc,mpiinfo)
   call smth_func(vela%z,nc,mpiinfo)
   call smth_func(dn   ,nc,mpiinfo)
   call smth_func(vel%x,nc,mpiinfo)
   call smth_func(vel%y,nc,mpiinfo)
   call smth_func(vel%z,nc,mpiinfo)
  end if

  !--Normalise to get coherent velocity
  !--NB(GC) dmin n'est pas utilise!
  where (dn < zero) 
   vel%x = zero;  vel%y = zero;  vel%z = zero
  elsewhere
   vel%x = vel%x/dn;   vel%y = vel%y/dn;   vel%z = vel%z/dn;
  end where

  __WRT_DEBUG_OUT("momtin")
 end subroutine momtin
 !*************************** END MOMTIN *****************************************

 !************************************* MOMT3D ************************************
 subroutine momt3d(imode)
  !--Calculation of moments in 2-D array.
  ! Charge and currents in DNA and UA.
  ! IMODE + 1 : Moments at t=0 (DNA,UA)
  ! IMODE = 2 : Arrays for current advance DN and U, DNF and UF
  ! IMODE = 3 : MASS density DN ; Coherent velocity U

  use defs_parametre,only   : gstep
  use defs_variable,only    : dn,dna,dn_e_pl,dn_e_incdt,dnf,vel,vela,velf,particule,nptot,&
       &                      s_min_loc,iwr
  use defs_grid
  use defs_mpitype,only     : mpiinfo
  use field_lissage,only          : smth_func


  integer,intent(in) :: imode

  __WRT_DEBUG_IN("momt3d")

  !--Initialise moment arrays to 0
  select case(imode)

   !----------------------------------------------------------------
   ! IMODE=1 
   !----------------------------------------------------------------
  case(1)
   !--4-current in dna,vela(1),...
   !call Amtsp3(particule,dna,vela,1,nptot,gstep,s_min_loc,mpiinfo,nc1)
   call Fmtsp3(particule,dna,dn_e_pl,dn_e_incdt,vela,1,nptot,gstep,s_min_loc,mpiinfo,nc1)

   if(iwr(6)==1)then
    call smth_func(dna,   nc,mpiinfo)
    call smth_func(vela%x,nc,mpiinfo)
    call smth_func(vela%y,nc,mpiinfo)
    call smth_func(vela%z,nc,mpiinfo)
   endif

   !----------------------------------------------------------------
   ! IMODE=2 
   !----------------------------------------------------------------
  case(2)
   !--4-current in dnf,velf(1,:,:,:),... and Arrays for current advance in dn,vel(1),...
   call Bmtsp3(particule,dnf,velf,dn,vel,1,nptot,gstep,s_min_loc,mpiinfo,nc1)

   if(iwr(6)==1)then 
    call smth_func(dnf,   nc,mpiinfo)
    call smth_func(velf%x,nc,mpiinfo)
    call smth_func(velf%y,nc,mpiinfo)
    call smth_func(velf%z,nc,mpiinfo)
   end if

   !----------------------------------------------------------------
   ! IMODE=3 (data output for analysis)
   !----------------------------------------------------------------
  case(3)
   !--4-current and momentum in dn,vel(1,:,:,:),...
   call Cmtsp3(particule,dna,vela,dn,vel,1,nptot,gstep,s_min_loc&
        &,mpiinfo,nc1)

   if(iwr(6)==1)then
    call smth_func(dna,   nc,mpiinfo)
    call smth_func(vela%x,nc,mpiinfo)
    call smth_func(vela%y,nc,mpiinfo)
    call smth_func(vela%z,nc,mpiinfo)
   end if

   !--Normalise to get coherent velocity
   !--NB(GC) dmin n'est pas utilise!
   where (dn <= zero) 
    vel%x = zero;  vel%y = zero;  vel%z = zero
   elsewhere
    vel%x = vel%x/dn;   vel%y = vel%y/dn;   vel%z = vel%z/dn;
   end where

  case default 
   print *, ' MOMT3D argument invalide'
  end select

  __WRT_DEBUG_OUT("momt3d")
 end subroutine momt3d
 !*********************************** MOMT3D ****************************************

 !************************************ MOMAD ****************************************
 subroutine momad()
  !--Half-step advance of currents
  ! The "special" moments are in (dn,vel(1),..) on the cell corner
  ! (uf) is advanced to (u) using (dn,u : "advanced currents")
  ! Interpolation trilineaire au sommet des cellules

  use defs_parametre,only        : dt,idisp,ipe,resis,dmin,gstep
  use defs_mpitype,only     : mpiinfo
  use defs_grid               
  use defs_variable,only    : dn,dnf,vel,velf,vela,&
       &                      te,pe,Efield,Bfield, &
       &                      resistivity,e_conv,e1_conv,&
       &                      rmu0,iwr
  use field_pe
  use field_lissage,only    : smth_func
  use field_e

  real(dp) :: dth
  
  __WRT_DEBUG_IN("momad")
  __GETTIME(42,1)!--Timer start

  call pecalc(ipe,te,dnf,pe)

  call ecalc3(Efield,Bfield,vela,&
       &      dnf,pe,resistivity,&
       &      e_conv,e1_conv,gstep,&
       &      dmin,rmu0,resis,&
       &      idisp,nc1, &
       &      mpiinfo)

  !--interpolate E to B-type grid (use USA as storage array)
  !--only when no NGP calculation
  call MEfield(Efield,nc1)

  dth = half*dt

  call momad3r(dn,vel,velf,Efield,Bfield,nc1,dth)
  
  !--Smoothing (optional)
  if(iwr(6) == 1)then
   call smth_func(vel%x,nc,mpiinfo)
   call smth_func(vel%y,nc,mpiinfo)
   call smth_func(vel%z,nc,mpiinfo)
  end if
  
  __GETTIME(42,2)!--Timer stop
  __WRT_DEBUG_OUT("momad")
 end subroutine momad
 !********************************* END MOMAD ***************************************

 !******************************** MOMAD3R ******************************************
 subroutine momad3r(dn,vel,velf,Efield,Bfield,nc1,dt)
  !--Cell-corner E and dn,u,uf,B identical grid

  real(dp),intent(in) :: dn(:,:,:)
  type(arr3Dtype),intent(inout) :: vel
  type(arr3Dtype),intent(in) :: velf
  type(arr3Dtype),intent(in) :: Efield,Bfield
  integer,intent(in) :: nc1(3)
  real(dp),intent(in) :: dt
  real(dp) :: veln(3),b_i(3)

  integer :: ii,jj,kk

  __WRT_DEBUG_IN("momad3r")
  __GETTIME(67,1)!--Timer start

! !$OMP PARALLEL PRIVATE(veln,b_i)
! !$OMP DO
  do kk = 1,nc1(3)
   do jj = 1,nc1(2)
    do ii = 1,nc1(1)
     veln(1) = vel%x(ii,jj,kk)
     veln(2) = vel%y(ii,jj,kk)
     veln(3) = vel%z(ii,jj,kk)

     b_i(1) = Bfield%x(ii,jj,kk)
     b_i(2) = Bfield%y(ii,jj,kk)
     b_i(3) = Bfield%z(ii,jj,kk)

     vel%x(ii,jj,kk) = velf%x(ii,jj,kk)+& 
          &         dt*( dn(ii,jj,kk)*Efield%x(ii,jj,kk)+&
          &              veln(2)*b_i(3)-veln(3)*b_i(2))

     vel%y(ii,jj,kk) = velf%y(ii,jj,kk)+&
          &         dt*( dn(ii,jj,kk)*Efield%y(ii,jj,kk)+&
          &              veln(3)*b_i(1)-veln(1)*b_i(3))

     vel%z(ii,jj,kk) = velf%z(ii,jj,kk)+&
          &         dt*( dn(ii,jj,kk)*Efield%z(ii,jj,kk)+&
          &              veln(1)*b_i(2)-veln(2)*b_i(1))

    enddo
   enddo
  enddo
! !$OMP END DO
! !$OMP END PARALLEL

  __GETTIME(67,2)!--Timer stop
  __WRT_DEBUG_OUT("momad3r")
 end subroutine momad3r
 !****************************** END MOMAD3R **************************************

 !********************************** AMTSP3 *****************************************
 subroutine Amtsp3(particule,dna,vela,n1,n2,gstep,s_min,&
      &            infompi,nc1)
  !--Calcul de la densite de charge (rho) au temps 0    -> dna
  !--et des composantes du courant (Ji) au temps 0      -> vela%x,vela%y,vela%z

  use defs_mpitype,only     : mpitype

  integer,intent(in) :: n1,n2
  integer,intent(in) :: nc1(3)
  real(dp),intent(in) :: gstep(3),s_min(3)
  type(mpitype),intent(in) :: infompi
  type(particletype),intent(in) :: particule(:)
  real(dp),intent(inout) :: dna(:,:,:)
  type(arr3Dtype),intent(inout) ::  vela

  integer :: nn,ijk(3)
  real(dp) :: sqp
  real(dp) :: w(nb_voisins)
  real(dp),dimension(3) :: s_f,s_a,v_p,s_m,gstep_inv

  __WRT_DEBUG_IN("Amtsp3")
  __GETTIME(68,1)!--Timer start

  !--Species partial moment arrays = 0
  dna = zero
  vela%x = zero;  vela%y = zero;  vela%z = zero

  gstep_inv = one/gstep

  do nn = n1,n2
   !--On collecte la composante X,Y,Z de la vitesse
   v_p = particule(nn)%vel
   sqp = particule(nn)%char !--On collecte la charge de la particule
   
   !--Relative position in cell s_f=(xf,yf,zf) center of the particule
   s_m = one + (particule(nn)%pos-s_min)*gstep_inv
   
   ijk = int(s_m)

#ifdef HAVE_DEBUG   
   if(any(ijk>nc1)) then
    print *,"ERROR Amtsp3",infompi%me,nn
    stop
   endif
#endif
   
   s_f = s_m-real(ijk,dp)
   
   !--Sequence of indices of B at cell corners
   !--Trilinear weight
   s_a = one-s_f

   w(1) = s_a(1)*s_a(2)*s_a(3)*sqp
   w(2) = s_f(1)*s_a(2)*s_a(3)*sqp
   w(3) = s_a(1)*s_f(2)*s_a(3)*sqp
   w(4) = s_f(1)*s_f(2)*s_a(3)*sqp
   w(5) = s_a(1)*s_a(2)*s_f(3)*sqp
   w(6) = s_f(1)*s_a(2)*s_f(3)*sqp
   w(7) = s_a(1)*s_f(2)*s_f(3)*sqp
   w(8) = s_f(1)*s_f(2)*s_f(3)*sqp

   !--Number density (in DN)
   dna(ijk(1)  ,ijk(2)  ,ijk(3)  ) = dna(ijk(1)  ,ijk(2)  ,ijk(3)  ) + w(1)
   dna(ijk(1)+1,ijk(2)  ,ijk(3)  ) = dna(ijk(1)+1,ijk(2)  ,ijk(3)  ) + w(2)
   dna(ijk(1)  ,ijk(2)+1,ijk(3)  ) = dna(ijk(1)  ,ijk(2)+1,ijk(3)  ) + w(3)
   dna(ijk(1)+1,ijk(2)+1,ijk(3)  ) = dna(ijk(1)+1,ijk(2)+1,ijk(3)  ) + w(4)
   dna(ijk(1)  ,ijk(2)  ,ijk(3)+1) = dna(ijk(1)  ,ijk(2)  ,ijk(3)+1) + w(5)
   dna(ijk(1)+1,ijk(2)  ,ijk(3)+1) = dna(ijk(1)+1,ijk(2)  ,ijk(3)+1) + w(6)
   dna(ijk(1)  ,ijk(2)+1,ijk(3)+1) = dna(ijk(1)  ,ijk(2)+1,ijk(3)+1) + w(7)
   dna(ijk(1)+1,ijk(2)+1,ijk(3)+1) = dna(ijk(1)+1,ijk(2)+1,ijk(3)+1) + w(8)

   !--vx-collection  (in UX)
   vela%x(ijk(1)  ,ijk(2)  ,ijk(3)  ) = vela%x(ijk(1)  ,ijk(2)  ,ijk(3)  ) + ((w(1)*v_p(1)))
   vela%x(ijk(1)+1,ijk(2)  ,ijk(3)  ) = vela%x(ijk(1)+1,ijk(2)  ,ijk(3)  ) + ((w(2)*v_p(1)))
   vela%x(ijk(1)  ,ijk(2)+1,ijk(3)  ) = vela%x(ijk(1)  ,ijk(2)+1,ijk(3)  ) + ((w(3)*v_p(1)))
   vela%x(ijk(1)+1,ijk(2)+1,ijk(3)  ) = vela%x(ijk(1)+1,ijk(2)+1,ijk(3)  ) + ((w(4)*v_p(1)))
   vela%x(ijk(1)  ,ijk(2)  ,ijk(3)+1) = vela%x(ijk(1)  ,ijk(2)  ,ijk(3)+1) + ((w(5)*v_p(1)))
   vela%x(ijk(1)+1,ijk(2)  ,ijk(3)+1) = vela%x(ijk(1)+1,ijk(2)  ,ijk(3)+1) + ((w(6)*v_p(1)))
   vela%x(ijk(1)  ,ijk(2)+1,ijk(3)+1) = vela%x(ijk(1)  ,ijk(2)+1,ijk(3)+1) + ((w(7)*v_p(1)))
   vela%x(ijk(1)+1,ijk(2)+1,ijk(3)+1) = vela%x(ijk(1)+1,ijk(2)+1,ijk(3)+1) + ((w(8)*v_p(1)))

   vela%y(ijk(1)  ,ijk(2)  ,ijk(3)  ) = vela%y(ijk(1)  ,ijk(2)  ,ijk(3)  ) + ((w(1)*v_p(2)))
   vela%y(ijk(1)+1,ijk(2)  ,ijk(3)  ) = vela%y(ijk(1)+1,ijk(2)  ,ijk(3)  ) + ((w(2)*v_p(2)))
   vela%y(ijk(1)  ,ijk(2)+1,ijk(3)  ) = vela%y(ijk(1)  ,ijk(2)+1,ijk(3)  ) + ((w(3)*v_p(2)))
   vela%y(ijk(1)+1,ijk(2)+1,ijk(3)  ) = vela%y(ijk(1)+1,ijk(2)+1,ijk(3)  ) + ((w(4)*v_p(2)))
   vela%y(ijk(1)  ,ijk(2)  ,ijk(3)+1) = vela%y(ijk(1)  ,ijk(2)  ,ijk(3)+1) + ((w(5)*v_p(2)))
   vela%y(ijk(1)+1,ijk(2)  ,ijk(3)+1) = vela%y(ijk(1)+1,ijk(2)  ,ijk(3)+1) + ((w(6)*v_p(2)))
   vela%y(ijk(1)  ,ijk(2)+1,ijk(3)+1) = vela%y(ijk(1)  ,ijk(2)+1,ijk(3)+1) + ((w(7)*v_p(2)))
   vela%y(ijk(1)+1,ijk(2)+1,ijk(3)+1) = vela%y(ijk(1)+1,ijk(2)+1,ijk(3)+1) + ((w(8)*v_p(2)))

   vela%z(ijk(1)  ,ijk(2)  ,ijk(3)  ) = vela%z(ijk(1)  ,ijk(2)  ,ijk(3)  ) + ((w(1)*v_p(3)))
   vela%z(ijk(1)+1,ijk(2)  ,ijk(3)  ) = vela%z(ijk(1)+1,ijk(2)  ,ijk(3)  ) + ((w(2)*v_p(3)))
   vela%z(ijk(1)  ,ijk(2)+1,ijk(3)  ) = vela%z(ijk(1)  ,ijk(2)+1,ijk(3)  ) + ((w(3)*v_p(3)))
   vela%z(ijk(1)+1,ijk(2)+1,ijk(3)  ) = vela%z(ijk(1)+1,ijk(2)+1,ijk(3)  ) + ((w(4)*v_p(3)))
   vela%z(ijk(1)  ,ijk(2)  ,ijk(3)+1) = vela%z(ijk(1)  ,ijk(2)  ,ijk(3)+1) + ((w(5)*v_p(3)))
   vela%z(ijk(1)+1,ijk(2)  ,ijk(3)+1) = vela%z(ijk(1)+1,ijk(2)  ,ijk(3)+1) + ((w(6)*v_p(3)))
   vela%z(ijk(1)  ,ijk(2)+1,ijk(3)+1) = vela%z(ijk(1)  ,ijk(2)+1,ijk(3)+1) + ((w(7)*v_p(3)))
   vela%z(ijk(1)+1,ijk(2)+1,ijk(3)+1) = vela%z(ijk(1)+1,ijk(2)+1,ijk(3)+1) + ((w(8)*v_p(3)))

  enddo
  __GETTIME(68,2)!--Timer stop
  
  !--Boundary conditions
  call cond_limit_func(dna,vela,infompi,nc1)
  
  __WRT_DEBUG_OUT("Amtsp3")
 end subroutine Amtsp3
 !**************************** AMTSP3 ***********************************************

 !**************************** BMTSP3 ***********************************************
 subroutine Bmtsp3(particule,dnf,velf,dn,vel,n1,n2,gstep,s_min,&
      &           infompi,nc1)
  !--calcul de la densite de charge (rho) au temps 1/2 -> dnf
  !  des composantes du courant au temps 1/2     -> velf%x,velf%y,velf%z
  !    de la quantite lambda(x1/2,v0)            -> dn
  !     des composantes de la quantite gamma(x1/2,v0) -> vel%x,vel%y,vel%z

  use defs_mpitype,only     : mpitype

  integer,intent(in) :: n1,n2
  integer,intent(in) :: nc1(3)
  real(dp),intent(in) :: gstep(3),s_min(3)
  type(mpitype),intent(in) :: infompi
  type(particletype),intent(in) :: particule(:)
  real(dp),intent(inout) :: dn(:,:,:)
  real(dp),intent(inout) :: dnf(:,:,:)
  type(arr3Dtype),intent(inout) ::  vel,velf

  integer :: nn,ijk(3)
  real(dp) :: w(nb_voisins),w0(nb_voisins)
  real(dp) :: sqp,smp,q2m
  real(dp),dimension(3) :: gstep_inv,s_f,s_a,v_p,s_m

  __WRT_DEBUG_IN("Bmtsp3")
  __GETTIME(69,1)!--Timer start

  !--Species partial moment arrays = 0
  dn = zero; vel%x = zero; vel%y = zero; vel%z = zero
  dnf = zero; velf%x = zero; velf%y = zero; velf%z = zero

  gstep_inv = one/gstep

  do nn = n1,n2
   
   !--On collecte la composante X,Y,Z de la vitesse
   v_p = particule(nn)%vel
   smp = particule(nn)%mass !--On collecte la masse de la particule
   sqp = particule(nn)%char !--On collecte la charge de la particule
   q2m = (sqp*sqp)/smp

   !--Relative position in cell s_f=(xf,yf,zf) center of the particule
   s_m = one +(particule(nn)%pos-s_min)*gstep_inv
 
   ijk = int(s_m)

#ifdef HAVE_DEBUG   
   if(any(ijk>nc1)) then
    print *,"ERROR Bmtsp3",infompi%me,nn
    print *,"Particule%vel ",particule(nn)%vel
    print *,"Particule%pos ",particule(nn)%pos
    print *,"Particule%mass ",particule(nn)%mass
    print *,"Particule%char ",particule(nn)%char
    print *,"Particule%dir ",particule(nn)%dir
    print *,"ijk ",ijk
    print *,"nc1 ",nc1
    stop
   endif
#endif

   s_f = s_m-real(ijk,dp)

   !--Sequence of indices of B at cell corners
   !--Trilinear weights
   s_a = one-s_f

   w0(1) = s_a(1)*s_a(2)*s_a(3)
   w0(2) = s_f(1)*s_a(2)*s_a(3)
   w0(3) = s_a(1)*s_f(2)*s_a(3)
   w0(4) = s_f(1)*s_f(2)*s_a(3)
   w0(5) = s_a(1)*s_a(2)*s_f(3)
   w0(6) = s_f(1)*s_a(2)*s_f(3)
   w0(7) = s_a(1)*s_f(2)*s_f(3)
   w0(8) = s_f(1)*s_f(2)*s_f(3)

   w = w0*sqp

   !--Number density (in DN)
   dnf(ijk(1)  ,ijk(2)  ,ijk(3)  ) = dnf(ijk(1)  ,ijk(2)  ,ijk(3)  ) + w(1)
   dnf(ijk(1)+1,ijk(2)  ,ijk(3)  ) = dnf(ijk(1)+1,ijk(2)  ,ijk(3)  ) + w(2)
   dnf(ijk(1)  ,ijk(2)+1,ijk(3)  ) = dnf(ijk(1)  ,ijk(2)+1,ijk(3)  ) + w(3)
   dnf(ijk(1)+1,ijk(2)+1,ijk(3)  ) = dnf(ijk(1)+1,ijk(2)+1,ijk(3)  ) + w(4)
   dnf(ijk(1)  ,ijk(2)  ,ijk(3)+1) = dnf(ijk(1)  ,ijk(2)  ,ijk(3)+1) + w(5)
   dnf(ijk(1)+1,ijk(2)  ,ijk(3)+1) = dnf(ijk(1)+1,ijk(2)  ,ijk(3)+1) + w(6)
   dnf(ijk(1)  ,ijk(2)+1,ijk(3)+1) = dnf(ijk(1)  ,ijk(2)+1,ijk(3)+1) + w(7)
   dnf(ijk(1)+1,ijk(2)+1,ijk(3)+1) = dnf(ijk(1)+1,ijk(2)+1,ijk(3)+1) + w(8)

   !--vx-collection  (in UX)
   w = w0*sqp*v_p(1)
   velf%x(ijk(1)  ,ijk(2)  ,ijk(3)  ) = velf%x(ijk(1)  ,ijk(2)  ,ijk(3)  ) + w(1)
   velf%x(ijk(1)+1,ijk(2)  ,ijk(3)  ) = velf%x(ijk(1)+1,ijk(2)  ,ijk(3)  ) + w(2)
   velf%x(ijk(1)  ,ijk(2)+1,ijk(3)  ) = velf%x(ijk(1)  ,ijk(2)+1,ijk(3)  ) + w(3)
   velf%x(ijk(1)+1,ijk(2)+1,ijk(3)  ) = velf%x(ijk(1)+1,ijk(2)+1,ijk(3)  ) + w(4)
   velf%x(ijk(1)  ,ijk(2)  ,ijk(3)+1) = velf%x(ijk(1)  ,ijk(2)  ,ijk(3)+1) + w(5)
   velf%x(ijk(1)+1,ijk(2)  ,ijk(3)+1) = velf%x(ijk(1)+1,ijk(2)  ,ijk(3)+1) + w(6)
   velf%x(ijk(1)  ,ijk(2)+1,ijk(3)+1) = velf%x(ijk(1)  ,ijk(2)+1,ijk(3)+1) + w(7)
   velf%x(ijk(1)+1,ijk(2)+1,ijk(3)+1) = velf%x(ijk(1)+1,ijk(2)+1,ijk(3)+1) + w(8)

   w = w0*sqp*v_p(2)
   velf%y(ijk(1)  ,ijk(2)  ,ijk(3)  ) = velf%y(ijk(1)  ,ijk(2)  ,ijk(3)  ) + w(1)
   velf%y(ijk(1)+1,ijk(2)  ,ijk(3)  ) = velf%y(ijk(1)+1,ijk(2)  ,ijk(3)  ) + w(2)
   velf%y(ijk(1)  ,ijk(2)+1,ijk(3)  ) = velf%y(ijk(1)  ,ijk(2)+1,ijk(3)  ) + w(3)
   velf%y(ijk(1)+1,ijk(2)+1,ijk(3)  ) = velf%y(ijk(1)+1,ijk(2)+1,ijk(3)  ) + w(4)
   velf%y(ijk(1)  ,ijk(2)  ,ijk(3)+1) = velf%y(ijk(1)  ,ijk(2)  ,ijk(3)+1) + w(5)
   velf%y(ijk(1)+1,ijk(2)  ,ijk(3)+1) = velf%y(ijk(1)+1,ijk(2)  ,ijk(3)+1) + w(6)
   velf%y(ijk(1)  ,ijk(2)+1,ijk(3)+1) = velf%y(ijk(1)  ,ijk(2)+1,ijk(3)+1) + w(7)
   velf%y(ijk(1)+1,ijk(2)+1,ijk(3)+1) = velf%y(ijk(1)+1,ijk(2)+1,ijk(3)+1) + w(8)

   w = w0*sqp*v_p(3)
   velf%z(ijk(1)  ,ijk(2)  ,ijk(3)  ) = velf%z(ijk(1)  ,ijk(2)  ,ijk(3)  ) + w(1)
   velf%z(ijk(1)+1,ijk(2)  ,ijk(3)  ) = velf%z(ijk(1)+1,ijk(2)  ,ijk(3)  ) + w(2)
   velf%z(ijk(1)  ,ijk(2)+1,ijk(3)  ) = velf%z(ijk(1)  ,ijk(2)+1,ijk(3)  ) + w(3)
   velf%z(ijk(1)+1,ijk(2)+1,ijk(3)  ) = velf%z(ijk(1)+1,ijk(2)+1,ijk(3)  ) + w(4)
   velf%z(ijk(1)  ,ijk(2)  ,ijk(3)+1) = velf%z(ijk(1)  ,ijk(2)  ,ijk(3)+1) + w(5)
   velf%z(ijk(1)+1,ijk(2)  ,ijk(3)+1) = velf%z(ijk(1)+1,ijk(2)  ,ijk(3)+1) + w(6)
   velf%z(ijk(1)  ,ijk(2)+1,ijk(3)+1) = velf%z(ijk(1)  ,ijk(2)+1,ijk(3)+1) + w(7)
   velf%z(ijk(1)+1,ijk(2)+1,ijk(3)+1) = velf%z(ijk(1)+1,ijk(2)+1,ijk(3)+1) + w(8)

   !--Number density (in DN)
   w = w0*q2m
   dn(ijk(1)  ,ijk(2)  ,ijk(3)  ) = dn(ijk(1)  ,ijk(2)  ,ijk(3)  ) + w(1)
   dn(ijk(1)+1,ijk(2)  ,ijk(3)  ) = dn(ijk(1)+1,ijk(2)  ,ijk(3)  ) + w(2)
   dn(ijk(1)  ,ijk(2)+1,ijk(3)  ) = dn(ijk(1)  ,ijk(2)+1,ijk(3)  ) + w(3)
   dn(ijk(1)+1,ijk(2)+1,ijk(3)  ) = dn(ijk(1)+1,ijk(2)+1,ijk(3)  ) + w(4)
   dn(ijk(1)  ,ijk(2)  ,ijk(3)+1) = dn(ijk(1)  ,ijk(2)  ,ijk(3)+1) + w(5)
   dn(ijk(1)+1,ijk(2)  ,ijk(3)+1) = dn(ijk(1)+1,ijk(2)  ,ijk(3)+1) + w(6)
   dn(ijk(1)  ,ijk(2)+1,ijk(3)+1) = dn(ijk(1)  ,ijk(2)+1,ijk(3)+1) + w(7)
   dn(ijk(1)+1,ijk(2)+1,ijk(3)+1) = dn(ijk(1)+1,ijk(2)+1,ijk(3)+1) + w(8)

   !--vx-collection  (in UX)
   w = w0*q2m*v_p(1)
   vel%x(ijk(1)  ,ijk(2)  ,ijk(3)  ) = vel%x(ijk(1)  ,ijk(2)  ,ijk(3)  ) + w(1)
   vel%x(ijk(1)+1,ijk(2)  ,ijk(3)  ) = vel%x(ijk(1)+1,ijk(2)  ,ijk(3)  ) + w(2)
   vel%x(ijk(1)  ,ijk(2)+1,ijk(3)  ) = vel%x(ijk(1)  ,ijk(2)+1,ijk(3)  ) + w(3)
   vel%x(ijk(1)+1,ijk(2)+1,ijk(3)  ) = vel%x(ijk(1)+1,ijk(2)+1,ijk(3)  ) + w(4)
   vel%x(ijk(1)  ,ijk(2)  ,ijk(3)+1) = vel%x(ijk(1)  ,ijk(2)  ,ijk(3)+1) + w(5)
   vel%x(ijk(1)+1,ijk(2)  ,ijk(3)+1) = vel%x(ijk(1)+1,ijk(2)  ,ijk(3)+1) + w(6)
   vel%x(ijk(1)  ,ijk(2)+1,ijk(3)+1) = vel%x(ijk(1)  ,ijk(2)+1,ijk(3)+1) + w(7)
   vel%x(ijk(1)+1,ijk(2)+1,ijk(3)+1) = vel%x(ijk(1)+1,ijk(2)+1,ijk(3)+1) + w(8)

   w = w0*q2m*v_p(2)
   vel%y(ijk(1)  ,ijk(2)  ,ijk(3)  ) = vel%y(ijk(1)  ,ijk(2)  ,ijk(3)  ) + w(1)
   vel%y(ijk(1)+1,ijk(2)  ,ijk(3)  ) = vel%y(ijk(1)+1,ijk(2)  ,ijk(3)  ) + w(2)
   vel%y(ijk(1)  ,ijk(2)+1,ijk(3)  ) = vel%y(ijk(1)  ,ijk(2)+1,ijk(3)  ) + w(3)
   vel%y(ijk(1)+1,ijk(2)+1,ijk(3)  ) = vel%y(ijk(1)+1,ijk(2)+1,ijk(3)  ) + w(4)
   vel%y(ijk(1)  ,ijk(2)  ,ijk(3)+1) = vel%y(ijk(1)  ,ijk(2)  ,ijk(3)+1) + w(5)
   vel%y(ijk(1)+1,ijk(2)  ,ijk(3)+1) = vel%y(ijk(1)+1,ijk(2)  ,ijk(3)+1) + w(6)
   vel%y(ijk(1)  ,ijk(2)+1,ijk(3)+1) = vel%y(ijk(1)  ,ijk(2)+1,ijk(3)+1) + w(7)
   vel%y(ijk(1)+1,ijk(2)+1,ijk(3)+1) = vel%y(ijk(1)+1,ijk(2)+1,ijk(3)+1) + w(8)

   w = w0*q2m*v_p(3)
   vel%z(ijk(1)  ,ijk(2)  ,ijk(3)  ) = vel%z(ijk(1)  ,ijk(2)  ,ijk(3)  ) + w(1)
   vel%z(ijk(1)+1,ijk(2)  ,ijk(3)  ) = vel%z(ijk(1)+1,ijk(2)  ,ijk(3)  ) + w(2)
   vel%z(ijk(1)  ,ijk(2)+1,ijk(3)  ) = vel%z(ijk(1)  ,ijk(2)+1,ijk(3)  ) + w(3)
   vel%z(ijk(1)+1,ijk(2)+1,ijk(3)  ) = vel%z(ijk(1)+1,ijk(2)+1,ijk(3)  ) + w(4)
   vel%z(ijk(1)  ,ijk(2)  ,ijk(3)+1) = vel%z(ijk(1)  ,ijk(2)  ,ijk(3)+1) + w(5)
   vel%z(ijk(1)+1,ijk(2)  ,ijk(3)+1) = vel%z(ijk(1)+1,ijk(2)  ,ijk(3)+1) + w(6)
   vel%z(ijk(1)  ,ijk(2)+1,ijk(3)+1) = vel%z(ijk(1)  ,ijk(2)+1,ijk(3)+1) + w(7)
   vel%z(ijk(1)+1,ijk(2)+1,ijk(3)+1) = vel%z(ijk(1)+1,ijk(2)+1,ijk(3)+1) + w(8)

  enddo

  __GETTIME(69,2)!--Timer stop
  !--Boundary conditions
  call cond_limit_func(dnf,velf,infompi,nc1)
  call cond_limit_func(dn,vel,infompi,nc1)

  __WRT_DEBUG_OUT("Bmtsp3")
 end subroutine Bmtsp3
 !******************************* END BMTSP3 ****************************************


 !******************************* CMTSP3 ********************************************
 subroutine Cmtsp3(particule,dna,vela,dn,vel,n1,n2,gstep,s_min,&   
      &            infompi,&
      &            nc1)

  !--Calcul de la densite de masse (rhom) au temps...   -> dn
  !--et des composantes de quantit de mouvement (Pi) au temps...-> vel%x,vel%y,vel%z

  use defs_mpitype,only     : mpitype

  integer,intent(in) :: n1,n2
  integer,intent(in) :: nc1(3)
  real(dp),intent(in) :: gstep(3),s_min(3)
  type(mpitype),intent(in) :: infompi
  type(particletype),intent(in) :: particule(:)
  real(dp),intent(inout) :: dn(:,:,:)
  real(dp),intent(inout) :: dna(:,:,:)
  type(arr3Dtype),intent(inout) :: vela,vel

  integer :: nn,ijk(3)
  real(dp) :: w(nb_voisins)
  real(dp) :: sqp,smp
  real(dp),dimension(3) :: gstep_inv,s_m,s_f,s_a,v_p,v_tmp

  __WRT_DEBUG_IN("Cmtsp3")
  __GETTIME(70,1)!--Timer start

  !--Species partial moment arrays = 0
  dn = zero; vel%x = zero; vel%y = zero; vel%z = zero
  dna = zero; vela%x = zero; vela%y = zero; vela%z = zero

  gstep_inv = one/gstep

  do nn = n1,n2
   !--On collecte la composante X,Y,Z de la vitesse
   v_p = particule(nn)%vel
   smp = particule(nn)%mass !--On collecte la masse de la particule
   sqp = particule(nn)%char !--On collecte la charge de la particule

   !--Relative position in cell s_f=(xf,yf,zf) center of the particule
   s_m = one +(particule(nn)%pos-s_min)*gstep_inv

   ijk = int(s_m)

#ifdef HAVE_DEBUG   
   if(any(ijk>nc1)) then
    print *,"ERROR Cmtsp3",infompi%me,nn
    stop
   endif
#endif

   s_f = s_m-real(ijk,dp)

   !--Sequence of indices of B at cell corners
   !--Trilinear weights
   s_a = one-s_f

   w(1) = s_a(1)*s_a(2)*s_a(3)
   w(2) = s_f(1)*s_a(2)*s_a(3)
   w(3) = s_a(1)*s_f(2)*s_a(3)
   w(4) = s_f(1)*s_f(2)*s_a(3)
   w(5) = s_a(1)*s_a(2)*s_f(3)
   w(6) = s_f(1)*s_a(2)*s_f(3)
   w(7) = s_a(1)*s_f(2)*s_f(3)
   w(8) = s_f(1)*s_f(2)*s_f(3)

   !--Number density (in DN)
   dna(ijk(1)  ,ijk(2)  ,ijk(3)  ) = dna(ijk(1)  ,ijk(2)  ,ijk(3)  ) + w(1)*sqp
   dna(ijk(1)+1,ijk(2)  ,ijk(3)  ) = dna(ijk(1)+1,ijk(2)  ,ijk(3)  ) + w(2)*sqp
   dna(ijk(1)  ,ijk(2)+1,ijk(3)  ) = dna(ijk(1)  ,ijk(2)+1,ijk(3)  ) + w(3)*sqp
   dna(ijk(1)+1,ijk(2)+1,ijk(3)  ) = dna(ijk(1)+1,ijk(2)+1,ijk(3)  ) + w(4)*sqp
   dna(ijk(1)  ,ijk(2)  ,ijk(3)+1) = dna(ijk(1)  ,ijk(2)  ,ijk(3)+1) + w(5)*sqp
   dna(ijk(1)+1,ijk(2)  ,ijk(3)+1) = dna(ijk(1)+1,ijk(2)  ,ijk(3)+1) + w(6)*sqp
   dna(ijk(1)  ,ijk(2)+1,ijk(3)+1) = dna(ijk(1)  ,ijk(2)+1,ijk(3)+1) + w(7)*sqp
   dna(ijk(1)+1,ijk(2)+1,ijk(3)+1) = dna(ijk(1)+1,ijk(2)+1,ijk(3)+1) + w(8)*sqp

   !--vx-collection  (in UX)
   v_tmp = sqp*v_p
   vela%x(ijk(1)  ,ijk(2)  ,ijk(3)  ) = vela%x(ijk(1)  ,ijk(2)  ,ijk(3)  ) + ((w(1)*v_tmp(1)))
   vela%x(ijk(1)+1,ijk(2)  ,ijk(3)  ) = vela%x(ijk(1)+1,ijk(2)  ,ijk(3)  ) + ((w(2)*v_tmp(1)))
   vela%x(ijk(1)  ,ijk(2)+1,ijk(3)  ) = vela%x(ijk(1)  ,ijk(2)+1,ijk(3)  ) + ((w(3)*v_tmp(1)))
   vela%x(ijk(1)+1,ijk(2)+1,ijk(3)  ) = vela%x(ijk(1)+1,ijk(2)+1,ijk(3)  ) + ((w(4)*v_tmp(1)))
   vela%x(ijk(1)  ,ijk(2)  ,ijk(3)+1) = vela%x(ijk(1)  ,ijk(2)  ,ijk(3)+1) + ((w(5)*v_tmp(1)))
   vela%x(ijk(1)+1,ijk(2)  ,ijk(3)+1) = vela%x(ijk(1)+1,ijk(2)  ,ijk(3)+1) + ((w(6)*v_tmp(1)))
   vela%x(ijk(1)  ,ijk(2)+1,ijk(3)+1) = vela%x(ijk(1)  ,ijk(2)+1,ijk(3)+1) + ((w(7)*v_tmp(1)))
   vela%x(ijk(1)+1,ijk(2)+1,ijk(3)+1) = vela%x(ijk(1)+1,ijk(2)+1,ijk(3)+1) + ((w(8)*v_tmp(1)))

   vela%y(ijk(1)  ,ijk(2)  ,ijk(3)  ) = vela%y(ijk(1)  ,ijk(2)  ,ijk(3)  ) + ((w(1)*v_tmp(2)))
   vela%y(ijk(1)+1,ijk(2)  ,ijk(3)  ) = vela%y(ijk(1)+1,ijk(2)  ,ijk(3)  ) + ((w(2)*v_tmp(2)))
   vela%y(ijk(1)  ,ijk(2)+1,ijk(3)  ) = vela%y(ijk(1)  ,ijk(2)+1,ijk(3)  ) + ((w(3)*v_tmp(2)))
   vela%y(ijk(1)+1,ijk(2)+1,ijk(3)  ) = vela%y(ijk(1)+1,ijk(2)+1,ijk(3)  ) + ((w(4)*v_tmp(2)))
   vela%y(ijk(1)  ,ijk(2)  ,ijk(3)+1) = vela%y(ijk(1)  ,ijk(2)  ,ijk(3)+1) + ((w(5)*v_tmp(2)))
   vela%y(ijk(1)+1,ijk(2)  ,ijk(3)+1) = vela%y(ijk(1)+1,ijk(2)  ,ijk(3)+1) + ((w(6)*v_tmp(2)))
   vela%y(ijk(1)  ,ijk(2)+1,ijk(3)+1) = vela%y(ijk(1)  ,ijk(2)+1,ijk(3)+1) + ((w(7)*v_tmp(2)))
   vela%y(ijk(1)+1,ijk(2)+1,ijk(3)+1) = vela%y(ijk(1)+1,ijk(2)+1,ijk(3)+1) + ((w(8)*v_tmp(2)))

   vela%z(ijk(1)  ,ijk(2)  ,ijk(3)  ) = vela%z(ijk(1)  ,ijk(2)  ,ijk(3)  ) + ((w(1)*v_tmp(3)))
   vela%z(ijk(1)+1,ijk(2)  ,ijk(3)  ) = vela%z(ijk(1)+1,ijk(2)  ,ijk(3)  ) + ((w(2)*v_tmp(3)))
   vela%z(ijk(1)  ,ijk(2)+1,ijk(3)  ) = vela%z(ijk(1)  ,ijk(2)+1,ijk(3)  ) + ((w(3)*v_tmp(3)))
   vela%z(ijk(1)+1,ijk(2)+1,ijk(3)  ) = vela%z(ijk(1)+1,ijk(2)+1,ijk(3)  ) + ((w(4)*v_tmp(3)))
   vela%z(ijk(1)  ,ijk(2)  ,ijk(3)+1) = vela%z(ijk(1)  ,ijk(2)  ,ijk(3)+1) + ((w(5)*v_tmp(3)))
   vela%z(ijk(1)+1,ijk(2)  ,ijk(3)+1) = vela%z(ijk(1)+1,ijk(2)  ,ijk(3)+1) + ((w(6)*v_tmp(3)))
   vela%z(ijk(1)  ,ijk(2)+1,ijk(3)+1) = vela%z(ijk(1)  ,ijk(2)+1,ijk(3)+1) + ((w(7)*v_tmp(3)))
   vela%z(ijk(1)+1,ijk(2)+1,ijk(3)+1) = vela%z(ijk(1)+1,ijk(2)+1,ijk(3)+1) + ((w(8)*v_tmp(3)))

   !--Number density (in DN)
   dn(ijk(1)  ,ijk(2)  ,ijk(3)  ) = dn(ijk(1)  ,ijk(2)  ,ijk(3)  ) + w(1)*smp
   dn(ijk(1)+1,ijk(2)  ,ijk(3)  ) = dn(ijk(1)+1,ijk(2)  ,ijk(3)  ) + w(2)*smp
   dn(ijk(1)  ,ijk(2)+1,ijk(3)  ) = dn(ijk(1)  ,ijk(2)+1,ijk(3)  ) + w(3)*smp
   dn(ijk(1)+1,ijk(2)+1,ijk(3)  ) = dn(ijk(1)+1,ijk(2)+1,ijk(3)  ) + w(4)*smp
   dn(ijk(1)  ,ijk(2)  ,ijk(3)+1) = dn(ijk(1)  ,ijk(2)  ,ijk(3)+1) + w(5)*smp
   dn(ijk(1)+1,ijk(2)  ,ijk(3)+1) = dn(ijk(1)+1,ijk(2)  ,ijk(3)+1) + w(6)*smp
   dn(ijk(1)  ,ijk(2)+1,ijk(3)+1) = dn(ijk(1)  ,ijk(2)+1,ijk(3)+1) + w(7)*smp
   dn(ijk(1)+1,ijk(2)+1,ijk(3)+1) = dn(ijk(1)+1,ijk(2)+1,ijk(3)+1) + w(8)*smp

   !--vx-collection  (in UX)
   v_tmp = smp*v_p
   vel%x(ijk(1)  ,ijk(2)  ,ijk(3)  ) = vel%x(ijk(1)  ,ijk(2)  ,ijk(3)  ) + ((w(1)*v_tmp(1)))
   vel%x(ijk(1)+1,ijk(2)  ,ijk(3)  ) = vel%x(ijk(1)+1,ijk(2)  ,ijk(3)  ) + ((w(2)*v_tmp(1)))
   vel%x(ijk(1)  ,ijk(2)+1,ijk(3)  ) = vel%x(ijk(1)  ,ijk(2)+1,ijk(3)  ) + ((w(3)*v_tmp(1)))
   vel%x(ijk(1)+1,ijk(2)+1,ijk(3)  ) = vel%x(ijk(1)+1,ijk(2)+1,ijk(3)  ) + ((w(4)*v_tmp(1)))
   vel%x(ijk(1)  ,ijk(2)  ,ijk(3)+1) = vel%x(ijk(1)  ,ijk(2)  ,ijk(3)+1) + ((w(5)*v_tmp(1)))
   vel%x(ijk(1)+1,ijk(2)  ,ijk(3)+1) = vel%x(ijk(1)+1,ijk(2)  ,ijk(3)+1) + ((w(6)*v_tmp(1)))
   vel%x(ijk(1)  ,ijk(2)+1,ijk(3)+1) = vel%x(ijk(1)  ,ijk(2)+1,ijk(3)+1) + ((w(7)*v_tmp(1)))
   vel%x(ijk(1)+1,ijk(2)+1,ijk(3)+1) = vel%x(ijk(1)+1,ijk(2)+1,ijk(3)+1) + ((w(8)*v_tmp(1)))

   vel%y(ijk(1)  ,ijk(2)  ,ijk(3)  ) = vel%y(ijk(1)  ,ijk(2)  ,ijk(3)  ) + ((w(1)*v_tmp(2)))
   vel%y(ijk(1)+1,ijk(2)  ,ijk(3)  ) = vel%y(ijk(1)+1,ijk(2)  ,ijk(3)  ) + ((w(2)*v_tmp(2)))
   vel%y(ijk(1)  ,ijk(2)+1,ijk(3)  ) = vel%y(ijk(1)  ,ijk(2)+1,ijk(3)  ) + ((w(3)*v_tmp(2)))
   vel%y(ijk(1)+1,ijk(2)+1,ijk(3)  ) = vel%y(ijk(1)+1,ijk(2)+1,ijk(3)  ) + ((w(4)*v_tmp(2)))
   vel%y(ijk(1)  ,ijk(2)  ,ijk(3)+1) = vel%y(ijk(1)  ,ijk(2)  ,ijk(3)+1) + ((w(5)*v_tmp(2)))
   vel%y(ijk(1)+1,ijk(2)  ,ijk(3)+1) = vel%y(ijk(1)+1,ijk(2)  ,ijk(3)+1) + ((w(6)*v_tmp(2)))
   vel%y(ijk(1)  ,ijk(2)+1,ijk(3)+1) = vel%y(ijk(1)  ,ijk(2)+1,ijk(3)+1) + ((w(7)*v_tmp(2)))
   vel%y(ijk(1)+1,ijk(2)+1,ijk(3)+1) = vel%y(ijk(1)+1,ijk(2)+1,ijk(3)+1) + ((w(8)*v_tmp(2)))

   vel%z(ijk(1)  ,ijk(2)  ,ijk(3)  ) = vel%z(ijk(1)  ,ijk(2)  ,ijk(3)  ) + ((w(1)*v_tmp(3)))
   vel%z(ijk(1)+1,ijk(2)  ,ijk(3)  ) = vel%z(ijk(1)+1,ijk(2)  ,ijk(3)  ) + ((w(2)*v_tmp(3)))
   vel%z(ijk(1)  ,ijk(2)+1,ijk(3)  ) = vel%z(ijk(1)  ,ijk(2)+1,ijk(3)  ) + ((w(3)*v_tmp(3)))
   vel%z(ijk(1)+1,ijk(2)+1,ijk(3)  ) = vel%z(ijk(1)+1,ijk(2)+1,ijk(3)  ) + ((w(4)*v_tmp(3)))
   vel%z(ijk(1)  ,ijk(2)  ,ijk(3)+1) = vel%z(ijk(1)  ,ijk(2)  ,ijk(3)+1) + ((w(5)*v_tmp(3)))
   vel%z(ijk(1)+1,ijk(2)  ,ijk(3)+1) = vel%z(ijk(1)+1,ijk(2)  ,ijk(3)+1) + ((w(6)*v_tmp(3)))
   vel%z(ijk(1)  ,ijk(2)+1,ijk(3)+1) = vel%z(ijk(1)  ,ijk(2)+1,ijk(3)+1) + ((w(7)*v_tmp(3)))
   vel%z(ijk(1)+1,ijk(2)+1,ijk(3)+1) = vel%z(ijk(1)+1,ijk(2)+1,ijk(3)+1) + ((w(8)*v_tmp(3)))
  enddo

  __GETTIME(70,2)!--Timer stop
  !--Boundary conditions
  call cond_limit_func(dna,vela,infompi,nc1)
  call cond_limit_func(dn,vel,infompi,nc1)

  __WRT_DEBUG_OUT("Cmtsp3")
 end subroutine Cmtsp3
 !************************************ END CMTSP3 **********************************

 !******************************* DMTSP3 ********************************************
 subroutine Dmtsp3(particule,dn,vel,n1,n2, &
      &            gstep,s_min, &
      &            infompi,&
      &            nc1)


  use defs_mpitype,only     : mpitype

  integer,intent(in) :: n1,n2
  integer,intent(in) :: nc1(3)
  real(dp),intent(in) :: gstep(3),s_min(3)
  type(mpitype),intent(in) :: infompi
  type(particletype),intent(in) :: particule(:)
  type(arr3Dtype),intent(inout) :: vel
  real(dp),intent(inout) :: dn(:,:,:)

  integer :: nn,ijk(3)
  real(dp) :: smp
  real(dp) :: w(nb_voisins)
  real(dp),dimension(3) :: gstep_inv,s_m,s_f,s_a,v_p

  __WRT_DEBUG_IN("Dmtsp3")
  __GETTIME(7,1)!--Timer start

  !--Species partial moment arrays = 0
  dn = zero
  vel%x = zero;  vel%y = zero;  vel%z = zero

  gstep_inv = one/gstep

  do nn = n1,n2
   !--On collecte la composante X,Y,Z de la vitesse
   v_p = particule(nn)%vel
   smp = particule(nn)%mass !--On collecte la masse de la particule
   
   !--Relative position in cell s_f=(xf,yf,zf) center of the particule
   s_m = one +(particule(nn)%pos-s_min)*gstep_inv

   ijk = int(s_m)

#ifdef HAVE_DEBUG   
   if(any(ijk>nc1)) then
    print *,"ERROR Dmtsp3",infompi%me,nn
    stop
   endif
#endif

   s_f = s_m-real(ijk,dp)

   !--Sequence of indices of B at cell corners
   !--Trilinear weights
   s_a = one-s_f

   w(1) = s_a(1)*s_a(2)*s_a(3)*smp
   w(2) = s_f(1)*s_a(2)*s_a(3)*smp
   w(3) = s_a(1)*s_f(2)*s_a(3)*smp
   w(4) = s_f(1)*s_f(2)*s_a(3)*smp
   w(5) = s_a(1)*s_a(2)*s_f(3)*smp
   w(6) = s_f(1)*s_a(2)*s_f(3)*smp
   w(7) = s_a(1)*s_f(2)*s_f(3)*smp
   w(8) = s_f(1)*s_f(2)*s_f(3)*smp

   !--Number density (in DN)
   dn(ijk(1)  ,ijk(2)  ,ijk(3)  ) = dn(ijk(1)  ,ijk(2)  ,ijk(3)  ) + w(1)
   dn(ijk(1)+1,ijk(2)  ,ijk(3)  ) = dn(ijk(1)+1,ijk(2)  ,ijk(3)  ) + w(2)
   dn(ijk(1)  ,ijk(2)+1,ijk(3)  ) = dn(ijk(1)  ,ijk(2)+1,ijk(3)  ) + w(3)
   dn(ijk(1)+1,ijk(2)+1,ijk(3)  ) = dn(ijk(1)+1,ijk(2)+1,ijk(3)  ) + w(4)
   dn(ijk(1)  ,ijk(2)  ,ijk(3)+1) = dn(ijk(1)  ,ijk(2)  ,ijk(3)+1) + w(5)
   dn(ijk(1)+1,ijk(2)  ,ijk(3)+1) = dn(ijk(1)+1,ijk(2)  ,ijk(3)+1) + w(6)
   dn(ijk(1)  ,ijk(2)+1,ijk(3)+1) = dn(ijk(1)  ,ijk(2)+1,ijk(3)+1) + w(7)
   dn(ijk(1)+1,ijk(2)+1,ijk(3)+1) = dn(ijk(1)+1,ijk(2)+1,ijk(3)+1) + w(8)

   !--vx-collection  (in UX)
   vel%x(ijk(1)  ,ijk(2)  ,ijk(3)  ) = vel%x(ijk(1)  ,ijk(2)  ,ijk(3)  ) + ((w(1)*v_p(1)))
   vel%x(ijk(1)+1,ijk(2)  ,ijk(3)  ) = vel%x(ijk(1)+1,ijk(2)  ,ijk(3)  ) + ((w(2)*v_p(1)))
   vel%x(ijk(1)  ,ijk(2)+1,ijk(3)  ) = vel%x(ijk(1)  ,ijk(2)+1,ijk(3)  ) + ((w(3)*v_p(1)))
   vel%x(ijk(1)+1,ijk(2)+1,ijk(3)  ) = vel%x(ijk(1)+1,ijk(2)+1,ijk(3)  ) + ((w(4)*v_p(1)))
   vel%x(ijk(1)  ,ijk(2)  ,ijk(3)+1) = vel%x(ijk(1)  ,ijk(2)  ,ijk(3)+1) + ((w(5)*v_p(1)))
   vel%x(ijk(1)+1,ijk(2)  ,ijk(3)+1) = vel%x(ijk(1)+1,ijk(2)  ,ijk(3)+1) + ((w(6)*v_p(1)))
   vel%x(ijk(1)  ,ijk(2)+1,ijk(3)+1) = vel%x(ijk(1)  ,ijk(2)+1,ijk(3)+1) + ((w(7)*v_p(1)))
   vel%x(ijk(1)+1,ijk(2)+1,ijk(3)+1) = vel%x(ijk(1)+1,ijk(2)+1,ijk(3)+1) + ((w(8)*v_p(1)))

   vel%y(ijk(1)  ,ijk(2)  ,ijk(3)  ) = vel%y(ijk(1)  ,ijk(2)  ,ijk(3)  ) + ((w(1)*v_p(2)))
   vel%y(ijk(1)+1,ijk(2)  ,ijk(3)  ) = vel%y(ijk(1)+1,ijk(2)  ,ijk(3)  ) + ((w(2)*v_p(2)))
   vel%y(ijk(1)  ,ijk(2)+1,ijk(3)  ) = vel%y(ijk(1)  ,ijk(2)+1,ijk(3)  ) + ((w(3)*v_p(2)))
   vel%y(ijk(1)+1,ijk(2)+1,ijk(3)  ) = vel%y(ijk(1)+1,ijk(2)+1,ijk(3)  ) + ((w(4)*v_p(2)))
   vel%y(ijk(1)  ,ijk(2)  ,ijk(3)+1) = vel%y(ijk(1)  ,ijk(2)  ,ijk(3)+1) + ((w(5)*v_p(2)))
   vel%y(ijk(1)+1,ijk(2)  ,ijk(3)+1) = vel%y(ijk(1)+1,ijk(2)  ,ijk(3)+1) + ((w(6)*v_p(2)))
   vel%y(ijk(1)  ,ijk(2)+1,ijk(3)+1) = vel%y(ijk(1)  ,ijk(2)+1,ijk(3)+1) + ((w(7)*v_p(2)))
   vel%y(ijk(1)+1,ijk(2)+1,ijk(3)+1) = vel%y(ijk(1)+1,ijk(2)+1,ijk(3)+1) + ((w(8)*v_p(2)))

   vel%z(ijk(1)  ,ijk(2)  ,ijk(3)  ) = vel%z(ijk(1)  ,ijk(2)  ,ijk(3)  ) + ((w(1)*v_p(3)))
   vel%z(ijk(1)+1,ijk(2)  ,ijk(3)  ) = vel%z(ijk(1)+1,ijk(2)  ,ijk(3)  ) + ((w(2)*v_p(3)))
   vel%z(ijk(1)  ,ijk(2)+1,ijk(3)  ) = vel%z(ijk(1)  ,ijk(2)+1,ijk(3)  ) + ((w(3)*v_p(3)))
   vel%z(ijk(1)+1,ijk(2)+1,ijk(3)  ) = vel%z(ijk(1)+1,ijk(2)+1,ijk(3)  ) + ((w(4)*v_p(3)))
   vel%z(ijk(1)  ,ijk(2)  ,ijk(3)+1) = vel%z(ijk(1)  ,ijk(2)  ,ijk(3)+1) + ((w(5)*v_p(3)))
   vel%z(ijk(1)+1,ijk(2)  ,ijk(3)+1) = vel%z(ijk(1)+1,ijk(2)  ,ijk(3)+1) + ((w(6)*v_p(3)))
   vel%z(ijk(1)  ,ijk(2)+1,ijk(3)+1) = vel%z(ijk(1)  ,ijk(2)+1,ijk(3)+1) + ((w(7)*v_p(3)))
   vel%z(ijk(1)+1,ijk(2)+1,ijk(3)+1) = vel%z(ijk(1)+1,ijk(2)+1,ijk(3)+1) + ((w(8)*v_p(3)))

  enddo

  __GETTIME(7,2)!--Timer stop
  call cond_limit_func(dn,vel,infompi,nc1)

  __WRT_DEBUG_OUT("Dmtsp3")
 end subroutine Dmtsp3

 !************************************ END DMTSP3 **********************************

 !******************************* EMTSP3 ********************************************
  subroutine Emtsp3(particule,dn,vel,n1,n2, &
      &             gstep,s_min, & 
      &             infompi,&
      &             nc1)

  use defs_mpitype,only     : mpitype

  integer,intent(in) :: n1,n2
  integer,intent(in) :: nc1(3)
  real(dp),intent(in) :: gstep(3),s_min(3)
  type(mpitype),intent(in) :: infompi
  type(particletype),intent(in) :: particule(:)
  type(arr3Dtype),intent(inout) :: vel
  real(dp),intent(inout) :: dn(:,:,:)

  integer :: nn,ijk(3)
  real(dp) :: sqp,smp,q2m
  real(dp) :: w(nb_voisins)
  real(dp),dimension(3) :: gstep_inv,s_m,s_f,s_a,v_p

  __WRT_DEBUG_IN("Emtsp3")
  __GETTIME(12,1)!--Timer start

  !--Species partial moment arrays = 0
  dn = zero
  vel%x = zero;  vel%y = zero;  vel%z = zero

  gstep_inv = one/gstep

  do nn = n1,n2
   !--On collecte la composante X,Y,Z de la vitesse
   v_p = particule(nn)%vel
   smp = particule(nn)%mass !--On collecte la masse de la particule
   sqp = particule(nn)%char !--On collecte la charge de la particule
   q2m = (sqp**two)/smp
   
   !--Relative position in cell s_f=(xf,yf,zf) center of the particule
   s_m = one +(particule(nn)%pos-s_min)*gstep_inv

   ijk = int(s_m)

#ifdef HAVE_DEBUG   
   if(any(ijk>nc1)) then
    print *,"ERROR Emtsp3",infompi%me,nn
    stop
   endif
#endif

   s_f = s_m-real(ijk,dp)

   !--Sequence of indices of B at cell corners
   !   (i,j,k)

   !--Trilinear weights
   s_a = one-s_f

   w(1) = s_a(1)*s_a(2)*s_a(3)*q2m
   w(2) = s_f(1)*s_a(2)*s_a(3)*q2m
   w(3) = s_a(1)*s_f(2)*s_a(3)*q2m
   w(4) = s_f(1)*s_f(2)*s_a(3)*q2m
   w(5) = s_a(1)*s_a(2)*s_f(3)*q2m
   w(6) = s_f(1)*s_a(2)*s_f(3)*q2m
   w(7) = s_a(1)*s_f(2)*s_f(3)*q2m
   w(8) = s_f(1)*s_f(2)*s_f(3)*q2m

   !--Number density (in DN)
   dn(ijk(1)  ,ijk(2)  ,ijk(3)  ) = dn(ijk(1)  ,ijk(2)  ,ijk(3)  ) + w(1)
   dn(ijk(1)+1,ijk(2)  ,ijk(3)  ) = dn(ijk(1)+1,ijk(2)  ,ijk(3)  ) + w(2)
   dn(ijk(1)  ,ijk(2)+1,ijk(3)  ) = dn(ijk(1)  ,ijk(2)+1,ijk(3)  ) + w(3)
   dn(ijk(1)+1,ijk(2)+1,ijk(3)  ) = dn(ijk(1)+1,ijk(2)+1,ijk(3)  ) + w(4)
   dn(ijk(1)  ,ijk(2)  ,ijk(3)+1) = dn(ijk(1)  ,ijk(2)  ,ijk(3)+1) + w(5)
   dn(ijk(1)+1,ijk(2)  ,ijk(3)+1) = dn(ijk(1)+1,ijk(2)  ,ijk(3)+1) + w(6)
   dn(ijk(1)  ,ijk(2)+1,ijk(3)+1) = dn(ijk(1)  ,ijk(2)+1,ijk(3)+1) + w(7)
   dn(ijk(1)+1,ijk(2)+1,ijk(3)+1) = dn(ijk(1)+1,ijk(2)+1,ijk(3)+1) + w(8)

   !--vx-collection  (in UX)
   vel%x(ijk(1)  ,ijk(2)  ,ijk(3)  ) = vel%x(ijk(1)  ,ijk(2)  ,ijk(3)  ) + ((w(1)*v_p(1)))
   vel%x(ijk(1)+1,ijk(2)  ,ijk(3)  ) = vel%x(ijk(1)+1,ijk(2)  ,ijk(3)  ) + ((w(2)*v_p(1)))
   vel%x(ijk(1)  ,ijk(2)+1,ijk(3)  ) = vel%x(ijk(1)  ,ijk(2)+1,ijk(3)  ) + ((w(3)*v_p(1)))
   vel%x(ijk(1)+1,ijk(2)+1,ijk(3)  ) = vel%x(ijk(1)+1,ijk(2)+1,ijk(3)  ) + ((w(4)*v_p(1)))
   vel%x(ijk(1)  ,ijk(2)  ,ijk(3)+1) = vel%x(ijk(1)  ,ijk(2)  ,ijk(3)+1) + ((w(5)*v_p(1)))
   vel%x(ijk(1)+1,ijk(2)  ,ijk(3)+1) = vel%x(ijk(1)+1,ijk(2)  ,ijk(3)+1) + ((w(6)*v_p(1)))
   vel%x(ijk(1)  ,ijk(2)+1,ijk(3)+1) = vel%x(ijk(1)  ,ijk(2)+1,ijk(3)+1) + ((w(7)*v_p(1)))
   vel%x(ijk(1)+1,ijk(2)+1,ijk(3)+1) = vel%x(ijk(1)+1,ijk(2)+1,ijk(3)+1) + ((w(8)*v_p(1)))

   vel%y(ijk(1)  ,ijk(2)  ,ijk(3)  ) = vel%y(ijk(1)  ,ijk(2)  ,ijk(3)  ) + ((w(1)*v_p(2)))
   vel%y(ijk(1)+1,ijk(2)  ,ijk(3)  ) = vel%y(ijk(1)+1,ijk(2)  ,ijk(3)  ) + ((w(2)*v_p(2)))
   vel%y(ijk(1)  ,ijk(2)+1,ijk(3)  ) = vel%y(ijk(1)  ,ijk(2)+1,ijk(3)  ) + ((w(3)*v_p(2)))
   vel%y(ijk(1)+1,ijk(2)+1,ijk(3)  ) = vel%y(ijk(1)+1,ijk(2)+1,ijk(3)  ) + ((w(4)*v_p(2)))
   vel%y(ijk(1)  ,ijk(2)  ,ijk(3)+1) = vel%y(ijk(1)  ,ijk(2)  ,ijk(3)+1) + ((w(5)*v_p(2)))
   vel%y(ijk(1)+1,ijk(2)  ,ijk(3)+1) = vel%y(ijk(1)+1,ijk(2)  ,ijk(3)+1) + ((w(6)*v_p(2)))
   vel%y(ijk(1)  ,ijk(2)+1,ijk(3)+1) = vel%y(ijk(1)  ,ijk(2)+1,ijk(3)+1) + ((w(7)*v_p(2)))
   vel%y(ijk(1)+1,ijk(2)+1,ijk(3)+1) = vel%y(ijk(1)+1,ijk(2)+1,ijk(3)+1) + ((w(8)*v_p(2)))

   vel%z(ijk(1)  ,ijk(2)  ,ijk(3)  ) = vel%z(ijk(1)  ,ijk(2)  ,ijk(3)  ) + ((w(1)*v_p(3)))
   vel%z(ijk(1)+1,ijk(2)  ,ijk(3)  ) = vel%z(ijk(1)+1,ijk(2)  ,ijk(3)  ) + ((w(2)*v_p(3)))
   vel%z(ijk(1)  ,ijk(2)+1,ijk(3)  ) = vel%z(ijk(1)  ,ijk(2)+1,ijk(3)  ) + ((w(3)*v_p(3)))
   vel%z(ijk(1)+1,ijk(2)+1,ijk(3)  ) = vel%z(ijk(1)+1,ijk(2)+1,ijk(3)  ) + ((w(4)*v_p(3)))
   vel%z(ijk(1)  ,ijk(2)  ,ijk(3)+1) = vel%z(ijk(1)  ,ijk(2)  ,ijk(3)+1) + ((w(5)*v_p(3)))
   vel%z(ijk(1)+1,ijk(2)  ,ijk(3)+1) = vel%z(ijk(1)+1,ijk(2)  ,ijk(3)+1) + ((w(6)*v_p(3)))
   vel%z(ijk(1)  ,ijk(2)+1,ijk(3)+1) = vel%z(ijk(1)  ,ijk(2)+1,ijk(3)+1) + ((w(7)*v_p(3)))
   vel%z(ijk(1)+1,ijk(2)+1,ijk(3)+1) = vel%z(ijk(1)+1,ijk(2)+1,ijk(3)+1) + ((w(8)*v_p(3)))

  enddo

  __GETTIME(12,2)!--Timer stop
  call cond_limit_func(dn,vel,infompi,nc1)

  __WRT_DEBUG_OUT("Emtsp3")
 end subroutine Emtsp3
 !************************************ END EMTSP3 **********************************
  !********************************** FMTSP3 *****************************************
 subroutine Fmtsp3(particule,dna,dn_e_pl,dn_e_incdt,vela,n1,n2,gstep,s_min,&
      &            infompi,nc1)
  !--Calcul de la densite de charge (rho) au temps 0    -> dna
  !--et des composantes du courant (Ji) au temps 0      -> vela%x,vela%y,vela%z

  use defs_mpitype,only     : mpitype

  integer,intent(in) :: n1,n2
  integer,intent(in) :: nc1(3)
  real(dp),intent(in) :: gstep(3),s_min(3)
  type(mpitype),intent(in) :: infompi
  type(particletype),intent(in) :: particule(:)
  real(dp),intent(inout) :: dna(:,:,:),dn_e_pl(:,:,:),dn_e_incdt(:,:,:)
  type(arr3Dtype),intent(inout) ::  vela

  integer :: nn,ijk(3)
  real(dp) :: sqp
  real(dp) :: w(nb_voisins)
  real(dp),dimension(3) :: s_f,s_a,v_p,s_m,gstep_inv

  __WRT_DEBUG_IN("Fmtsp3")
  __GETTIME(68,1)!--Timer start

  !--Species partial moment arrays = 0
  dna = zero
  dn_e_pl = zero
  dn_e_incdt = zero
  vela%x = zero;  vela%y = zero;  vela%z = zero

  gstep_inv = one/gstep

  do nn = n1,n2
   !--On collecte la composante X,Y,Z de la vitesse
   v_p = particule(nn)%vel
   sqp = particule(nn)%char !--On collecte la charge de la particule
   
   !--Relative position in cell s_f=(xf,yf,zf) center of the particule
   s_m = one + (particule(nn)%pos-s_min)*gstep_inv
   
   ijk = int(s_m)

#ifdef HAVE_DEBUG   
   if(any(ijk>nc1)) then
    print *,"ERROR Fmtsp3",infompi%me,nn
    print *,"Particule%vel ",particule(nn)%vel
    print *,"Particule%pos ",particule(nn)%pos
    print *,"Particule%mass ",particule(nn)%mass
    print *,"Particule%char ",particule(nn)%char
    print *,"Particule%dir ",particule(nn)%dir
    print *,"ijk ",ijk
    print *,"nc1 ",nc1
    stop
   endif
#endif
   
   s_f = s_m-real(ijk,dp)
   
   !--Sequence of indices of B at cell corners
   !--Trilinear weight
   s_a = one-s_f

   w(1) = s_a(1)*s_a(2)*s_a(3)*sqp
   w(2) = s_f(1)*s_a(2)*s_a(3)*sqp
   w(3) = s_a(1)*s_f(2)*s_a(3)*sqp
   w(4) = s_f(1)*s_f(2)*s_a(3)*sqp
   w(5) = s_a(1)*s_a(2)*s_f(3)*sqp
   w(6) = s_f(1)*s_a(2)*s_f(3)*sqp
   w(7) = s_a(1)*s_f(2)*s_f(3)*sqp
   w(8) = s_f(1)*s_f(2)*s_f(3)*sqp
   
   if ((particule(nn)%orig == zero).and.(particule(nn)%exc == zero)) then
   
     dn_e_incdt(ijk(1)  ,ijk(2)  ,ijk(3)  ) = dn_e_incdt(ijk(1)  ,ijk(2)  ,ijk(3)  ) + w(1)
     dn_e_incdt(ijk(1)+1,ijk(2)  ,ijk(3)  ) = dn_e_incdt(ijk(1)+1,ijk(2)  ,ijk(3)  ) + w(2)
     dn_e_incdt(ijk(1)  ,ijk(2)+1,ijk(3)  ) = dn_e_incdt(ijk(1)  ,ijk(2)+1,ijk(3)  ) + w(3)
     dn_e_incdt(ijk(1)+1,ijk(2)+1,ijk(3)  ) = dn_e_incdt(ijk(1)+1,ijk(2)+1,ijk(3)  ) + w(4)
     dn_e_incdt(ijk(1)  ,ijk(2)  ,ijk(3)+1) = dn_e_incdt(ijk(1)  ,ijk(2)  ,ijk(3)+1) + w(5)
     dn_e_incdt(ijk(1)+1,ijk(2)  ,ijk(3)+1) = dn_e_incdt(ijk(1)+1,ijk(2)  ,ijk(3)+1) + w(6)
     dn_e_incdt(ijk(1)  ,ijk(2)+1,ijk(3)+1) = dn_e_incdt(ijk(1)  ,ijk(2)+1,ijk(3)+1) + w(7)
     dn_e_incdt(ijk(1)+1,ijk(2)+1,ijk(3)+1) = dn_e_incdt(ijk(1)+1,ijk(2)+1,ijk(3)+1) + w(8)
 
   else
     dn_e_pl(ijk(1)  ,ijk(2)  ,ijk(3)  ) = dn_e_pl(ijk(1)  ,ijk(2)  ,ijk(3)  ) + w(1)
     dn_e_pl(ijk(1)+1,ijk(2)  ,ijk(3)  ) = dn_e_pl(ijk(1)+1,ijk(2)  ,ijk(3)  ) + w(2)
     dn_e_pl(ijk(1)  ,ijk(2)+1,ijk(3)  ) = dn_e_pl(ijk(1)  ,ijk(2)+1,ijk(3)  ) + w(3)
     dn_e_pl(ijk(1)+1,ijk(2)+1,ijk(3)  ) = dn_e_pl(ijk(1)+1,ijk(2)+1,ijk(3)  ) + w(4)
     dn_e_pl(ijk(1)  ,ijk(2)  ,ijk(3)+1) = dn_e_pl(ijk(1)  ,ijk(2)  ,ijk(3)+1) + w(5)
     dn_e_pl(ijk(1)+1,ijk(2)  ,ijk(3)+1) = dn_e_pl(ijk(1)+1,ijk(2)  ,ijk(3)+1) + w(6)
     dn_e_pl(ijk(1)  ,ijk(2)+1,ijk(3)+1) = dn_e_pl(ijk(1)  ,ijk(2)+1,ijk(3)+1) + w(7)
     dn_e_pl(ijk(1)+1,ijk(2)+1,ijk(3)+1) = dn_e_pl(ijk(1)+1,ijk(2)+1,ijk(3)+1) + w(8)
  endif
  
   !--vx-collection  (in UX)
   vela%x(ijk(1)  ,ijk(2)  ,ijk(3)  ) = vela%x(ijk(1)  ,ijk(2)  ,ijk(3)  ) + ((w(1)*v_p(1)))
   vela%x(ijk(1)+1,ijk(2)  ,ijk(3)  ) = vela%x(ijk(1)+1,ijk(2)  ,ijk(3)  ) + ((w(2)*v_p(1)))
   vela%x(ijk(1)  ,ijk(2)+1,ijk(3)  ) = vela%x(ijk(1)  ,ijk(2)+1,ijk(3)  ) + ((w(3)*v_p(1)))
   vela%x(ijk(1)+1,ijk(2)+1,ijk(3)  ) = vela%x(ijk(1)+1,ijk(2)+1,ijk(3)  ) + ((w(4)*v_p(1)))
   vela%x(ijk(1)  ,ijk(2)  ,ijk(3)+1) = vela%x(ijk(1)  ,ijk(2)  ,ijk(3)+1) + ((w(5)*v_p(1)))
   vela%x(ijk(1)+1,ijk(2)  ,ijk(3)+1) = vela%x(ijk(1)+1,ijk(2)  ,ijk(3)+1) + ((w(6)*v_p(1)))
   vela%x(ijk(1)  ,ijk(2)+1,ijk(3)+1) = vela%x(ijk(1)  ,ijk(2)+1,ijk(3)+1) + ((w(7)*v_p(1)))
   vela%x(ijk(1)+1,ijk(2)+1,ijk(3)+1) = vela%x(ijk(1)+1,ijk(2)+1,ijk(3)+1) + ((w(8)*v_p(1)))

   vela%y(ijk(1)  ,ijk(2)  ,ijk(3)  ) = vela%y(ijk(1)  ,ijk(2)  ,ijk(3)  ) + ((w(1)*v_p(2)))
   vela%y(ijk(1)+1,ijk(2)  ,ijk(3)  ) = vela%y(ijk(1)+1,ijk(2)  ,ijk(3)  ) + ((w(2)*v_p(2)))
   vela%y(ijk(1)  ,ijk(2)+1,ijk(3)  ) = vela%y(ijk(1)  ,ijk(2)+1,ijk(3)  ) + ((w(3)*v_p(2)))
   vela%y(ijk(1)+1,ijk(2)+1,ijk(3)  ) = vela%y(ijk(1)+1,ijk(2)+1,ijk(3)  ) + ((w(4)*v_p(2)))
   vela%y(ijk(1)  ,ijk(2)  ,ijk(3)+1) = vela%y(ijk(1)  ,ijk(2)  ,ijk(3)+1) + ((w(5)*v_p(2)))
   vela%y(ijk(1)+1,ijk(2)  ,ijk(3)+1) = vela%y(ijk(1)+1,ijk(2)  ,ijk(3)+1) + ((w(6)*v_p(2)))
   vela%y(ijk(1)  ,ijk(2)+1,ijk(3)+1) = vela%y(ijk(1)  ,ijk(2)+1,ijk(3)+1) + ((w(7)*v_p(2)))
   vela%y(ijk(1)+1,ijk(2)+1,ijk(3)+1) = vela%y(ijk(1)+1,ijk(2)+1,ijk(3)+1) + ((w(8)*v_p(2)))

   vela%z(ijk(1)  ,ijk(2)  ,ijk(3)  ) = vela%z(ijk(1)  ,ijk(2)  ,ijk(3)  ) + ((w(1)*v_p(3)))
   vela%z(ijk(1)+1,ijk(2)  ,ijk(3)  ) = vela%z(ijk(1)+1,ijk(2)  ,ijk(3)  ) + ((w(2)*v_p(3)))
   vela%z(ijk(1)  ,ijk(2)+1,ijk(3)  ) = vela%z(ijk(1)  ,ijk(2)+1,ijk(3)  ) + ((w(3)*v_p(3)))
   vela%z(ijk(1)+1,ijk(2)+1,ijk(3)  ) = vela%z(ijk(1)+1,ijk(2)+1,ijk(3)  ) + ((w(4)*v_p(3)))
   vela%z(ijk(1)  ,ijk(2)  ,ijk(3)+1) = vela%z(ijk(1)  ,ijk(2)  ,ijk(3)+1) + ((w(5)*v_p(3)))
   vela%z(ijk(1)+1,ijk(2)  ,ijk(3)+1) = vela%z(ijk(1)+1,ijk(2)  ,ijk(3)+1) + ((w(6)*v_p(3)))
   vela%z(ijk(1)  ,ijk(2)+1,ijk(3)+1) = vela%z(ijk(1)  ,ijk(2)+1,ijk(3)+1) + ((w(7)*v_p(3)))
   vela%z(ijk(1)+1,ijk(2)+1,ijk(3)+1) = vela%z(ijk(1)+1,ijk(2)+1,ijk(3)+1) + ((w(8)*v_p(3)))

  enddo
  dna = dn_e_pl+dn_e_incdt
  
  __GETTIME(68,2)!--Timer stop

  !--Boundary conditions
  call cond_limit_func(dna,vela,infompi,nc1)
  call cond_limit_func(dn_e_pl,infompi,nc1)
  call cond_limit_func(dn_e_incdt,infompi,nc1)
  
  __WRT_DEBUG_OUT("Fmtsp3")
 end subroutine Fmtsp3
 !**************************** FMTSP3 ***********************************************
end module part_moment
