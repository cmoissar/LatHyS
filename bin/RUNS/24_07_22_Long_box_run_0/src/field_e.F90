!!=============================================================
!!=============================================================
!!module: m_ecalc
!! NAME
!!  m_ecalc (MMancini,RModolo)
!!
!! FUNCTION
!!  Contains the routine for the Efield computation
!! NOTE
module field_e

 use defs_basis
 use defs_arr3Dtype
 use m_writeout
 use m_timing,only         : time_get

#include "q-p_common.h"

 implicit none
 private

 public ::         &
      ecalc3,      &!--Compute Electric field
      MEfield       !--Average of Electric on cells

contains
 !!#####################################################################

 !!=============================================================
 !!routine: m_ecalc/ecalc3
 !!
 !! FUNCTION
 !!  Compute Electric field
 !! IN 
 !!  idisp
 !!  dmin
 !!  rmu0
 !!  resis=resistivity coefficient
 !!  infompi(mpitype)=mpi informations
 !!  nc1(3)=grid size -1
 !!  gstep(3)=Grid step
 !!  e_conv(3),e1_conv(3)=Electric Field at x=0,x=1
 !!  Bfield(arr3Dtype)=Magntic Field
 !!  vel(:,:,:)=Velocity Field
 !!  dn(:,:,:)=Density Field
 !!  pe(:,:,:)=Electronic Pressure Field
 !!  resistivity(:,:,:)=Resistivity Field
 !!
 !! OUT
 !! SIDE EFFECT
 !!  Efield(arr3Dtype)=Electric Field
 !! NOTES
 !!  Electric field calculation (3-D). It is calculated on an interlaced grid,
 !!  size (nc1(1)+1,nc1(2)+1,nc1(3)+1), from B,dn,Ji,pe on the (nc1(1),nc1(2),nc1(3)) grid.
 subroutine ecalc3(Efield,Bfield,vel,&
      &            dn,pe,resistivity,&
      &            e_conv,e1_conv,gstep,&
      &            dmin,rmu0,resis,&
      &            idisp,nc1,&
      &            infompi)

  use defs_mpitype,only     : mpitype
  use field_cond_limit,only       : cond_limit_func
  use defs_parametre,only : dn_lim_inf_conduc,planetname
  use defs_grid,only :ncm
  use defs_variable,only : dn_e_pl,dn_e_incdt,s_min,s_min_loc,s_max,s_max_loc

  integer,intent(in) :: idisp
  integer,intent(in) :: nc1(3)
  real(dp),intent(in) :: dmin,rmu0,resis
  real(dp),intent(in) :: gstep(3) !--Grid step
  real(dp),intent(in) :: e_conv(3),e1_conv(3)
  type(mpitype),intent(in) :: infompi
  type(arr3Dtype),intent(inout) :: Efield
  type(arr3Dtype),intent(in) :: Bfield,vel
  real(dp),dimension(:,:,:),intent(in) :: dn,pe
  real(dp),dimension(:,:,:),intent(inout) :: resistivity
  !--resistivity is intent(inout) because is not used 
  !  if used put it intent(in)

  type(arr3Dtype)::Tfield
  integer :: ii,jj,kk,size_smooth_E
  real(dp) :: rmu0i,dni,resisloc,dni_lim,dni_pl,dni_inc,resis_m,x,y,z,weight
  real(dp),dimension(3) :: e_local,b_i,curl_b,resis_v
  real(dp),dimension(3) :: gstep_invq,g_p,u_i,u_b,r_b
  real(dp),dimension(8) :: pe_local
  real(dp),dimension(8,3) :: b_eight,selector
  ! character(len=500)::msg

  __WRT_DEBUG_IN("ecalc3")
  !__GETTIME(51,1)!--Timer start

  !--gstep_invq=1/(4*dx), 1/(4*dy) and 1/(4*dz) for differential operators
  gstep_invq = quarter/gstep
  dni_lim=1./dn_lim_inf_conduc
  rmu0i = one/rmu0

  !--Defining selector to mix Bfield composent
  selector(:,1) = gstep_invq(1)*real((/ -1, 1,-1, 1,-1, 1,-1, 1 /),dp)  
  selector(:,2) = gstep_invq(2)*real((/ -1,-1, 1, 1,-1,-1, 1, 1 /),dp)
  selector(:,3) = gstep_invq(3)*real((/ -1,-1,-1,-1, 1, 1, 1, 1 /),dp)   

  !--Include dispersion if IDISP = 1
  if(idisp == 1)then

   do kk = 2,nc1(3)
    do jj = 2,nc1(2)
     do ii = 2,nc1(1)

      !--Charge density at cell centre
      !__GETTIME(52,1)!--Timer start dni in ecalc3
!      dni_pl = eighth*(&
!           & dn_e_pl(ii-1,jj-1,kk-1) + dn_e_pl(ii  ,jj-1,kk-1)+ &
!           & dn_e_pl(ii-1,jj  ,kk-1) + dn_e_pl(ii  ,jj  ,kk-1)+  &
!           & dn_e_pl(ii-1,jj-1,kk  ) + dn_e_pl(ii  ,jj-1,kk  )+ &
!           & dn_e_pl(ii-1,jj  ,kk  ) + dn_e_pl(ii  ,jj  ,kk  ))
!      dni_inc = eighth*(&
!           & dn_e_incdt(ii-1,jj-1,kk-1) + dn_e_incdt(ii  ,jj-1,kk-1)+ &
!           & dn_e_incdt(ii-1,jj  ,kk-1) + dn_e_incdt(ii  ,jj  ,kk-1)+  &
!           & dn_e_incdt(ii-1,jj-1,kk  ) + dn_e_incdt(ii  ,jj-1,kk  )+ &
!           & dn_e_incdt(ii-1,jj  ,kk  ) + dn_e_incdt(ii  ,jj  ,kk  ))
!      dni=dni_pl+dni_inc
!      if (dni.gt.0._dp)then
!		resis_m=(dni_inc/30._dp+dni_pl)/dni
!       else
                resis_m=1._dp
!	endif
      dni = eighth*(&
           & dn(ii-1,jj-1,kk-1) + dn(ii  ,jj-1,kk-1)+ &
           & dn(ii-1,jj  ,kk-1) + dn(ii  ,jj  ,kk-1)+  &
           & dn(ii-1,jj-1,kk  ) + dn(ii  ,jj-1,kk  )+ &
           & dn(ii-1,jj  ,kk  ) + dn(ii  ,jj  ,kk  ))
       !__GETTIME(52,2)!--Timer stop dni in ecalc3

      !--E index
      if(dni>=dmin)then
       !--Inverse density
       dni = one/dni

       !__GETTIME(53,1)!--Timer start b_eight in ecalc3
       !--Recording the twice used part of Bfield
       b_eight(:,1) = (/ Bfield%x(ii-1,jj-1,kk-1),Bfield%x(ii  ,jj-1,kk-1), &
            &            Bfield%x(ii-1,jj  ,kk-1),Bfield%x(ii  ,jj  ,kk-1), &
            &            Bfield%x(ii-1,jj-1,kk  ),Bfield%x(ii  ,jj-1,kk  ), &
            &            Bfield%x(ii-1,jj  ,kk  ),Bfield%x(ii  ,jj  ,kk  ) /)
       b_eight(:,2) = (/ Bfield%y(ii-1,jj-1,kk-1),Bfield%y(ii  ,jj-1,kk-1), &
            &            Bfield%y(ii-1,jj  ,kk-1),Bfield%y(ii  ,jj  ,kk-1), &
            &            Bfield%y(ii-1,jj-1,kk  ),Bfield%y(ii  ,jj-1,kk  ), &
            &            Bfield%y(ii-1,jj  ,kk  ),Bfield%y(ii  ,jj  ,kk  ) /)
       b_eight(:,3) = (/ Bfield%z(ii-1,jj-1,kk-1),Bfield%z(ii  ,jj-1,kk-1), &
            &            Bfield%z(ii-1,jj  ,kk-1),Bfield%z(ii  ,jj  ,kk-1), &
            &            Bfield%z(ii-1,jj-1,kk  ),Bfield%z(ii  ,jj-1,kk  ), &
            &            Bfield%z(ii-1,jj  ,kk  ),Bfield%z(ii  ,jj  ,kk  ) /)
       !__GETTIME(53,2)!--Timer stop b_eight in ecalc3

       !--<B>,<Ji> at cell centre
       !__GETTIME(54,1)!--Timer start JxB in ecalc3
       b_i(:) = eighth*sum(b_eight,dim=1)

       u_i(1) = eighth*&
            &         ( vel%x(ii-1,jj-1,kk-1) + vel%x(ii  ,jj-1,kk-1) &
            &         + vel%x(ii-1,jj  ,kk-1) + vel%x(ii  ,jj  ,kk-1) &
            &         + vel%x(ii-1,jj-1,kk  ) + vel%x(ii  ,jj-1,kk  ) &
            &         + vel%x(ii-1,jj  ,kk  ) + vel%x(ii  ,jj  ,kk  ) )
       u_i(2) = eighth*&
            &         ( vel%y(ii-1,jj-1,kk-1) + vel%y(ii  ,jj-1,kk-1) &
            &         + vel%y(ii-1,jj  ,kk-1) + vel%y(ii  ,jj  ,kk-1) &
            &         + vel%y(ii-1,jj-1,kk  ) + vel%y(ii  ,jj-1,kk  ) &
            &         + vel%y(ii-1,jj  ,kk  ) + vel%y(ii  ,jj  ,kk  ) )
       u_i(3) = eighth*&
            &         ( vel%z(ii-1,jj-1,kk-1) + vel%z(ii  ,jj-1,kk-1) &
            &         + vel%z(ii-1,jj  ,kk-1) + vel%z(ii  ,jj  ,kk-1) &
            &         + vel%z(ii-1,jj-1,kk  ) + vel%z(ii  ,jj-1,kk  ) &
            &         + vel%z(ii-1,jj  ,kk  ) + vel%z(ii  ,jj  ,kk  ))


       !--<Ji> X <B>
       !u_b = u_i.cross.b_i
       u_b = (/u_i(2)*b_i(3)-u_i(3)*b_i(2),&
            &  u_i(3)*b_i(1)-u_i(1)*b_i(3),&
            &  u_i(1)*b_i(2)-u_i(2)*b_i(1) /) 
       !__GETTIME(54,2)!--Timer stop JxB in ecalc3

       !--grad(PE)
       !__GETTIME(55,1)!--Timer start gred_PE in ecalc3
       !--Recording pe_local
       pe_local = (/ pe(ii-1,jj-1,kk-1),pe(ii  ,jj-1,kk-1), &
            &        pe(ii-1,jj  ,kk-1),pe(ii  ,jj  ,kk-1), &
            &        pe(ii-1,jj-1,kk  ),pe(ii  ,jj-1,kk  ),&
            &        pe(ii-1,jj  ,kk  ),pe(ii  ,jj  ,kk  )/) 
       g_p = matmul(pe_local,selector)
       !__GETTIME(55,2)!--Timer stop gred_PE in ecalc3

       !--No dispersion if IDISP = 0 
       !--E = - (U X B + GRAD(PE) ) 
       !--wihtout dispersion no resistivity
       e_local = -dni*(u_b + g_p)

       !--curl B0
       !__GETTIME(56,1)!--Timer start curl_B in ecalc3
       ! curl_b.x = dbzdy - dbydz
       ! curl_b.y = dbxdz - dbzdx
       ! curl_b.z = dbydx - dbxdy
       curl_b(1) = sum(selector(:,2)*b_eight(:,3)-selector(:,3)*b_eight(:,2))
       curl_b(2) = sum(selector(:,3)*b_eight(:,1)-selector(:,1)*b_eight(:,3)) 
       curl_b(3) = sum(selector(:,1)*b_eight(:,2)-selector(:,2)*b_eight(:,1))
       !__GETTIME(56,2)!--Timer stop curl_B in ecalc3

       !__GETTIME(57,1)!--Timer start e_local in ecalc3
       !--(curl B) x B
       !r_b = curl_b .cross. b_i
       r_b = (/curl_b(2)*b_i(3)-curl_b(3)*b_i(2),&
            &  curl_b(3)*b_i(1)-curl_b(1)*b_i(3),& 
            &  curl_b(1)*b_i(2)-curl_b(2)*b_i(1) /) 

       !--resistivity actually is not used
       ! resisloc = eighth*(resistivity(ii-1,jj-1,kk-1) + resistivity(ii  ,jj-1,kk-1) &
       !      &           + resistivity(ii-1,jj  ,kk-1) + resistivity(ii  ,jj  ,kk-1) &
       !      &           + resistivity(ii-1,jj-1,kk  ) + resistivity(ii  ,jj-1,kk  ) &
       !      &           + resistivity(ii-1,jj  ,kk  ) + resistivity(ii  ,jj  ,kk  ))

       resisloc = resis*resis_m!*100._dp
       
       !--E = [curl(B) x B /(MU0) - (U X B) - GRAD(PE) + RESIS (curl B)/ MUO0 ]/N
       e_local = e_local + dni*rmu0i*r_b+resisloc*(curl_b-b_i*dot_product(b_i,curl_b)/dot_product(b_i,b_i))


       if (dni.lt.dni_lim) e_local = -dni*g_p+ dni*rmu0i*r_b+resisloc*(curl_b-b_i*dot_product(b_i,curl_b)/dot_product(b_i,b_i))


       if(any(abs(e_local) >40.0_dp)) then
        e_local(1) = quarter*(Efield%x(ii,jj-1,kk  ) + Efield%x(ii,jj+1,kk  )&
             &              + Efield%x(ii,jj  ,kk-1) + Efield%x(ii,jj  ,kk+1))                                                   
        e_local(2) = quarter*(Efield%y(ii,jj-1,kk  ) + Efield%y(ii,jj+1,kk  )&
             &              + Efield%y(ii,jj  ,kk-1) + Efield%y(ii,jj  ,kk+1))        
        e_local(3) = quarter*(Efield%z(ii,jj-1,kk  ) + Efield%z(ii,jj+1,kk  )&
             &              + Efield%z(ii,jj  ,kk-1) + Efield%z(ii,jj  ,kk+1))                                                   
       end if
       !__GETTIME(57,2)!--Timer start e_local in ecalc3

      else
       e_local = zero
      end if

      Efield%x(ii,jj,kk) = e_local(1)
      Efield%y(ii,jj,kk) = e_local(2)
      Efield%z(ii,jj,kk) = e_local(3)

     enddo
    enddo
   enddo
  else
   do kk = 2,nc1(3)
    do jj = 2,nc1(2)
     do ii = 2,nc1(1)

      !--Charge density at cell centre
      !__GETTIME(52,1)!--Timer start dni in ecalc3
      dni = eighth*(&
           & dn(ii-1,jj-1,kk-1) + dn(ii  ,jj-1,kk-1)+ &
           & dn(ii-1,jj  ,kk-1) + dn(ii  ,jj  ,kk-1)+  &
           & dn(ii-1,jj-1,kk  ) + dn(ii  ,jj-1,kk  )+ &
           & dn(ii-1,jj  ,kk  ) + dn(ii  ,jj  ,kk  ))
      !__GETTIME(52,2)!--Timer stop dni in ecalc3

      !--E index
      if(dni>=dmin)then
       !--Inverse density
       dni = one/dni

       !__GETTIME(53,1)!--Timer start b_eight in ecalc3
       !--Recording the twice used part of Bfield
       b_eight(:,1) = (/ Bfield%x(ii-1,jj-1,kk-1),Bfield%x(ii  ,jj-1,kk-1), &
            &            Bfield%x(ii-1,jj  ,kk-1),Bfield%x(ii  ,jj  ,kk-1), &
            &            Bfield%x(ii-1,jj-1,kk  ),Bfield%x(ii  ,jj-1,kk  ), &
            &            Bfield%x(ii-1,jj  ,kk  ),Bfield%x(ii  ,jj  ,kk  ) /)
       b_eight(:,2) = (/ Bfield%y(ii-1,jj-1,kk-1),Bfield%y(ii  ,jj-1,kk-1), &
            &            Bfield%y(ii-1,jj  ,kk-1),Bfield%y(ii  ,jj  ,kk-1), &
            &            Bfield%y(ii-1,jj-1,kk  ),Bfield%y(ii  ,jj-1,kk  ), &
            &            Bfield%y(ii-1,jj  ,kk  ),Bfield%y(ii  ,jj  ,kk  ) /)
       b_eight(:,3) = (/ Bfield%z(ii-1,jj-1,kk-1),Bfield%z(ii  ,jj-1,kk-1), &
            &            Bfield%z(ii-1,jj  ,kk-1),Bfield%z(ii  ,jj  ,kk-1), &
            &            Bfield%z(ii-1,jj-1,kk  ),Bfield%z(ii  ,jj-1,kk  ), &
            &            Bfield%z(ii-1,jj  ,kk  ),Bfield%z(ii  ,jj  ,kk  ) /)
       !__GETTIME(53,2)!--Timer stop b_eight in ecalc3

       !--<B>,<Ji> at cell centre
       !__GETTIME(54,1)!--Timer start JxB in ecalc3

       b_i(:) = eighth*sum(b_eight,dim=1)

       u_i(1) = eighth*&
            &         ( vel%x(ii-1,jj-1,kk-1) + vel%x(ii  ,jj-1,kk-1) &
            &         + vel%x(ii-1,jj  ,kk-1) + vel%x(ii  ,jj  ,kk-1) &
            &         + vel%x(ii-1,jj-1,kk  ) + vel%x(ii  ,jj-1,kk  ) &
            &         + vel%x(ii-1,jj  ,kk  ) + vel%x(ii  ,jj  ,kk  ) )
       u_i(2) = eighth*&
            &         ( vel%y(ii-1,jj-1,kk-1) + vel%y(ii  ,jj-1,kk-1) &
            &         + vel%y(ii-1,jj  ,kk-1) + vel%y(ii  ,jj  ,kk-1) &
            &         + vel%y(ii-1,jj-1,kk  ) + vel%y(ii  ,jj-1,kk  ) &
            &         + vel%y(ii-1,jj  ,kk  ) + vel%y(ii  ,jj  ,kk  ) )
       u_i(3) = eighth*&
            &         ( vel%z(ii-1,jj-1,kk-1) + vel%z(ii  ,jj-1,kk-1) &
            &         + vel%z(ii-1,jj  ,kk-1) + vel%z(ii  ,jj  ,kk-1) &
            &         + vel%z(ii-1,jj-1,kk  ) + vel%z(ii  ,jj-1,kk  ) &
            &         + vel%z(ii-1,jj  ,kk  ) + vel%z(ii  ,jj  ,kk  ))

       !--<Ji> X <B>
       !u_b = u_i.cross.b_i
       u_b = (/u_i(2)*b_i(3)-u_i(3)*b_i(2),&
            &  u_i(3)*b_i(1)-u_i(1)*b_i(3),&
            &  u_i(1)*b_i(2)-u_i(2)*b_i(1) /) 
       !__GETTIME(54,2)!--Timer stop JxB in ecalc3

       !--grad(PE)
       !__GETTIME(55,1)!--Timer start gred_PE in ecalc3
       !--Recording pe_local
       pe_local = (/ pe(ii-1,jj-1,kk-1),pe(ii  ,jj-1,kk-1), &
            &        pe(ii-1,jj  ,kk-1),pe(ii  ,jj  ,kk-1), &
            &        pe(ii-1,jj-1,kk  ),pe(ii  ,jj-1,kk  ),&
            &        pe(ii-1,jj  ,kk  ),pe(ii  ,jj  ,kk  )/) 
       g_p = matmul(pe_local,selector)
       !__GETTIME(55,2)!--Timer stop gred_PE in ecalc3

       !--No dispersion if IDISP = 0 
       !--E = - (U X B + GRAD(PE) ) / N
       !--wihtout dispersion no resistivity
       e_local = -dni*(u_b + g_p)


      else
       e_local = zero
      end if

      Efield%x(ii,jj,kk) = e_local(1)
      Efield%y(ii,jj,kk) = e_local(2)
      Efield%z(ii,jj,kk) = e_local(3)

     enddo
    enddo
   enddo
  endif

  ! !--Smoothing
Tfield=Efield
  do kk = 2,nc1(3)
   do jj = 2,nc1(2)
     do ii = 2,nc1(1)
            dni=(2.*Efield%x(ii,jj,kk)-Efield%x(ii-1,jj,kk)-Efield%x(ii+1,jj,kk))**2+&
    &           (2.*Efield%y(ii,jj,kk)-Efield%y(ii,jj-1,kk)-Efield%y(ii,jj+1,kk))**2+&
    &           (2.*Efield%z(ii,jj,kk)-Efield%z(ii,jj,kk-1)-Efield%z(ii,jj,kk+1))**2
            dni=sqrt(dni)
        if (dni.gt.10._dp) then
        Tfield%x(ii,jj,kk) = sum(Efield%x(max(ii-1,1):min(ii+1,ncm(1)),max(jj-1,1):min(jj+1,ncm(2)),max(kk-1,1):min(kk+1,ncm(3))))&
             &      /(max(1,size(Efield%x(max(ii-1,1):min(ii+1,ncm(1)),max(jj-1,1):min(jj+1,ncm(2)),max(kk-1,1):min(kk+1,ncm(3))))))
        Tfield%y(ii,jj,kk) = sum(Efield%y(max(ii-1,1):min(ii+1,ncm(1)),max(jj-1,1):min(jj+1,ncm(2)),max(kk-1,1):min(kk+1,ncm(3))))&
             &      /(max(1,size(Efield%y(max(ii-1,1):min(ii+1,ncm(1)),max(jj-1,1):min(jj+1,ncm(2)),max(kk-1,1):min(kk+1,ncm(3))))))
        Tfield%z(ii,jj,kk) = sum(Efield%z(max(ii-1,1):min(ii+1,ncm(1)),max(jj-1,1):min(jj+1,ncm(2)),max(kk-1,1):min(kk+1,ncm(3))))&
             &      /(max(1,size(Efield%z(max(ii-1,1):min(ii+1,ncm(1)),max(jj-1,1):min(jj+1,ncm(2)),max(kk-1,1):min(kk+1,ncm(3))))))
       end if
      enddo
     enddo
    enddo
     
   
  Efield=Tfield
  ! Add smoothing close to simulation box boundaries to dump waves
  if (trim(planetname) == 'ganymede') then
  
  Tfield = Efield
  size_smooth_E = 5*gstep(1)
  weight = 3._dp/4._dp
   
  do kk = 2,nc1(3)
    do jj = 2,nc1(2)
      do ii = 2,nc1(1)
        ! check if the grid point is in the smoothing region
        x = (ii-1)*gstep(1)+s_min_loc(1)
        y = (jj-1)*gstep(2)+s_min_loc(2)
        z = (kk-1)*gstep(3)+s_min_loc(3)
       ! if (((x >= s_max(1)-size_smooth_E).or.&
       !    &	((y<=s_min(2)+size_smooth_E).or.(y>=s_max(2)-size_smooth_E))).or.&
       !    &	((z<=s_min(3)+size_smooth_E).or.(z>=s_max(3)-size_smooth_E))) then
       if (((y<=s_min(2)+size_smooth_E).or.(y>=s_max(2)-size_smooth_E)).or. &
           ((z<=s_min(3)+size_smooth_E).or.(z>=s_max(3)-size_smooth_E))) then

        Tfield%x(ii,jj,kk) = weight*sum(Efield%x(max(ii-1,1):min(ii+1,ncm(1)),&
        & max(jj-1,1):min(jj+1,ncm(2)),max(kk-1,1):min(kk+1,ncm(3))))&
        &      /(max(1,size(Efield%x(max(ii-1,1):min(ii+1,ncm(1)),&
        & max(jj-1,1):min(jj+1,ncm(2)),max(kk-1,1):min(kk+1,ncm(3)))))) + &
             &          (1.-weight)*Efield%x(ii,jj,kk)
        Tfield%y(ii,jj,kk) = weight*sum(Efield%y(max(ii-1,1):min(ii+1,ncm(1)),&
        & max(jj-1,1):min(jj+1,ncm(2)),max(kk-1,1):min(kk+1,ncm(3))))&
        &      /(max(1,size(Efield%y(max(ii-1,1):min(ii+1,ncm(1)),&
        & max(jj-1,1):min(jj+1,ncm(2)),max(kk-1,1):min(kk+1,ncm(3))))))+ &
             &          (1.-weight)*Efield%y(ii,jj,kk)              
        Tfield%z(ii,jj,kk) = weight*sum(Efield%z(max(ii-1,1):min(ii+1,ncm(1)),&
        & max(jj-1,1):min(jj+1,ncm(2)),max(kk-1,1):min(kk+1,ncm(3))))&
        &      /(max(1,size(Efield%z(max(ii-1,1):min(ii+1,ncm(1)),&
        & max(jj-1,1):min(jj+1,ncm(2)),max(kk-1,1):min(kk+1,ncm(3))))))  + &   
        &          (1.-weight)*Efield%z(ii,jj,kk)            
       
!           if  (((y<=s_min(2)+size_smooth_E).or.(y>=s_max(2)-size_smooth_E)).or.&
!           &    ((z<=s_min(3)+size_smooth_E).or.(z>=s_max(3)-size_smooth_E))) then
               Tfield%x(ii,jj,kk) = (Tfield%x(ii,jj,kk)-e_conv(1))*0.5_dp+e_conv(1)
               Tfield%y(ii,jj,kk) = (Tfield%y(ii,jj,kk)-e_conv(2))*0.5_dp+e_conv(2)
               Tfield%z(ii,jj,kk) = (Tfield%z(ii,jj,kk)-e_conv(3))*0.5_dp+e_conv(3) 
!           endif
           
           ! we smooth following X direction
!           Tfield%x(ii,jj,kk) = quarter*Tfield%x(ii-1,jj,kk) + half*Tfield%x(ii,jj,kk) + &
!           	&	quarter*Tfield%x(ii+1,jj,kk)
!	   Tfield%y(ii,jj,kk) = quarter*Tfield%y(ii-1,jj,kk) + half*Tfield%y(ii,jj,kk) + &
!           	&	quarter*Tfield%y(ii+1,jj,kk)
!           Tfield%z(ii,jj,kk) = quarter*Tfield%z(ii-1,jj,kk) + half*Tfield%z(ii,jj,kk) + &
!           	&	quarter*Tfield%z(ii+1,jj,kk)
!	   ! we smooth following Y direction
!	   Tfield%x(ii,jj,kk) = quarter*Tfield%x(ii,jj-1,kk) + half*Tfield%x(ii,jj,kk) + &
!	   	&	quarter*Tfield%x(ii,jj+1,kk)
!     	   Tfield%y(ii,jj,kk) = quarter*Tfield%y(ii,jj-1,kk) + half*Tfield%y(ii,jj,kk) + &
!	   	&	quarter*Tfield%y(ii,jj+1,kk)         	
!           Tfield%z(ii,jj,kk) = quarter*Tfield%z(ii,jj-1,kk) + half*Tfield%z(ii,jj,kk) + &
!	   	&	quarter*Tfield%z(ii,jj+1,kk)
!	   ! we smooth following Z direction
!	   Tfield%x(ii,jj,kk) = quarter*Tfield%x(ii,jj,kk-1) + half*Tfield%x(ii,jj,kk) + &
!	   	   	&	quarter*Tfield%x(ii,jj,kk+1)
!	   Tfield%y(ii,jj,kk) = quarter*Tfield%y(ii,jj,kk-1) + half*Tfield%y(ii,jj,kk) + &
!	   	   	&	quarter*Tfield%y(ii,jj,kk+1)         	
!	   Tfield%z(ii,jj,kk) = quarter*Tfield%z(ii,jj,kk-1) + half*Tfield%z(ii,jj,kk) + &
!	   		&	quarter*Tfield%z(ii,jj,kk+1)
        endif ! we are in the smoothing region
      enddo
    enddo 
  enddo
  Efield = Tfield
  endif   ! end of smoothing for Ganymede
    
    


  !--Boundaries (periodic)
  !__GETTIME(58,1)!--Timer start cond_limit in ecalc3
  call cond_limit_func(Efield,infompi,nc1+1,e_conv,e1_conv)

  !__GETTIME(58,2)!--Timer stop cond_limit in ecalc3

  !__GETTIME(51,2)!--Timer stop




  __WRT_DEBUG_OUT("ecalc3")

 end subroutine ecalc3

 !!=============================================================
 !!routine: m_ecalc/ecalc3
 !!
 !! FUNCTION
 !!  Average on any cell (so on the grid of the Magnetic Field)
 !!  of the Electric Field
 !! IN 
 !!  nc1(3)=grid size -1
 !! OUT
 !! SIDE EFFECT
 !!  Efield(arr3Dtype)=Electric Field
 !! NOTES
 !!  Electric field calculation (3-D). It is calculated on an interlaced grid,
 !!  size (nc1(1)+1,nc1(2)+1,nc1(3)+1), from B,dn,Ji,pe on the (nc1(1),nc1(2),nc1(3)) grid.
 subroutine MEfield(Efield,nc1)

  integer,intent(in) :: nc1(3)
  type(arr3Dtype),intent(inout) :: Efield

  __WRT_DEBUG_IN("MEfield")

  !--Upbound of Efield components are nc1+1
  Efield%x(1:nc1(1),1:nc1(2),1:nc1(3)) = eighth*(&
       &     Efield%x(1:nc1(1),1:nc1(2),1:nc1(3)) &
       &   + Efield%x(2:      ,1:nc1(2),1:nc1(3)) &
       &   + Efield%x(2:      ,2:      ,1:nc1(3)) &
       &   + Efield%x(1:nc1(1),2:      ,1:nc1(3)) &
       &   + Efield%x(1:nc1(1),1:nc1(2),2:      ) &    
       &   + Efield%x(2:      ,1:nc1(2),2:      ) &
       &   + Efield%x(2:      ,2:      ,2:      ) &
       &   + Efield%x(1:nc1(1),2:      ,2:      ))

  Efield%y(1:nc1(1),1:nc1(2),1:nc1(3)) = eighth*(&
       &     Efield%y(1:nc1(1),1:nc1(2),1:nc1(3)) &
       &   + Efield%y(2:      ,1:nc1(2),1:nc1(3)) &
       &   + Efield%y(2:      ,2:      ,1:nc1(3)) &
       &   + Efield%y(1:nc1(1),2:      ,1:nc1(3)) &
       &   + Efield%y(1:nc1(1),1:nc1(2),2:      ) &
       &   + Efield%y(2:      ,1:nc1(2),2:      ) &
       &   + Efield%y(2:      ,2:      ,2:      ) &
       &   + Efield%y(1:nc1(1),2:      ,2:      ))

  Efield%z(1:nc1(1),1:nc1(2),1:nc1(3)) = eighth*(&
       &     Efield%z(1:nc1(1),1:nc1(2),1:nc1(3)) &
       &   + Efield%z(2:      ,1:nc1(2),1:nc1(3)) &
       &   + Efield%z(2:      ,2:      ,1:nc1(3)) &
       &   + Efield%z(1:nc1(1),2:      ,1:nc1(3)) &
       &   + Efield%z(1:nc1(1),1:nc1(2),2:      ) &
       &   + Efield%z(2:      ,1:nc1(2),2:      ) &
       &   + Efield%z(2:      ,2:      ,2:      ) &
       &   + Efield%z(1:nc1(1),2:      ,2:      ))

  __WRT_DEBUG_OUT("MEfield")
 end subroutine MEfield

end module field_e
