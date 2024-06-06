!=============================================================
!=============================================================
module m_current_computation
 use defs_basis
 use defs_variable
 use defs_grid
 use defs_arr3Dtype
 use defs_parametre, only : fildat,gstep,ns,planetname
 use m_writeout
 use netcdf
 use defs_basic_cdf
 use diag_wrt_common_cdf
#include "q-p_common.h"

 implicit none

 private :: &
    read_field_cdf, &
    read_dim_field_cdf, &
    header_2Dcut_J_VOTABLE, &
    header_2Dcut_JxB_VOTABLE, &
    extract_XY, &
    getstrf

 public :: &
    compute_current

contains

subroutine compute_current(run_name)
   ! this routine computes the currents
   character(len=*),intent(in) :: run_name
   !Local Variable
   character(len=64) :: varname1,varname2,varname3,varname4
   character(len=4)  :: prefix,prefix_output,prefixb
   character(len=25) :: filename
   real(dp),dimension(:,:,:),allocatable :: Ax, Ay, Az
   real(dp),dimension(:,:,:),allocatable :: Jx, Jy, Jz,Jtot
   real(dp),dimension(:,:,:),allocatable :: JcrossB_x,JcrossB_y,JcrossB_z,JcrossB_tot
   real(dp),dimension(:,:,:),allocatable :: GradPmag_x,GradPmag_y,GradPmag_z,Pmag,GradPmag_tot
      integer :: ncm_tot(3)
   real(dp)  :: centr(3),radius,gstep(3),c_wpi,gstep_invq(3),B0
   real(dp) :: bxi,byi,bzi
   integer :: ival,jval,kval,ii,jj,kk
    character(len=8) :: planet
    real(dp) :: sgn,nrm_J
    character(len=5) :: coord

   prefix = "Magw"
   varname1 = "Bx"
   varname2 = "By"
   varname3 = "Bz"
   filename = prefix//'_'//trim(run_name)//'.nc'


  call read_dim_field_cdf(filename,varname1,ncm_tot,radius,centr,gstep,planet,c_wpi,B0)

   allocate(Ax(ncm_tot(1),ncm_tot(2),ncm_tot(3)))
   allocate(Ay(ncm_tot(1),ncm_tot(2),ncm_tot(3)))
   allocate(Az(ncm_tot(1),ncm_tot(2),ncm_tot(3)))
   !allocate(Amod(ncm_tot(1),ncm_tot(2),ncm_tot(3)))
 Ax = zero;      Ay = zero;     Az = zero;!
    allocate(Pmag(ncm_tot(1),ncm_tot(2),ncm_tot(3)))
    Pmag = zero


  call read_field_cdf(filename,varname1,varname2,varname3,Ax,Ay,Az)
  !*******
  ! renormalizing in simulation unit
  sgn = -1
  coord = 'MSO'
  nrm_J = 1./(4*pi*1.e-7*c_wpi*1000) ! unit of J nA/m2
Ax = -Ax
Ay = -Ay

   Pmag = (Ax**2+Ay**2+Az**2)



    gstep_invq = (quarter)/(gstep)
        allocate(Jx(ncm_tot(1),ncm_tot(2),ncm_tot(3)));   Jx(:,:,:) = zero
        allocate(Jy(ncm_tot(1),ncm_tot(2),ncm_tot(3)));   Jy(:,:,:) = zero
        allocate(Jz(ncm_tot(1),ncm_tot(2),ncm_tot(3)));   Jz(:,:,:) = zero
        allocate(Jtot(ncm_tot(1),ncm_tot(2),ncm_tot(3)));  Jtot(:,:,:) = zero


        Jx(1:ncm_tot(1)-1,1:ncm_tot(2)-1,1:ncm_tot(3)-1) = gstep_invq(2)*(&
               &   - Az(1:ncm_tot(1)-1,1:ncm_tot(2)-1,1:ncm_tot(3)-1) &
               &   - Az(2:ncm_tot(1)  ,1:ncm_tot(2)-1,1:ncm_tot(3)-1) &
               &   + Az(1:ncm_tot(1)-1,2:ncm_tot(2)  ,1:ncm_tot(3)-1) &
               &   + Az(2:ncm_tot(1)  ,2:ncm_tot(2)  ,1:ncm_tot(3)-1) &
               &   - Az(1:ncm_tot(1)-1,1:ncm_tot(2)-1,2:ncm_tot(3)  ) &
               &   - Az(2:ncm_tot(1)  ,1:ncm_tot(2)-1,2:ncm_tot(3)  ) &
               &   + Az(1:ncm_tot(1)-1,2:ncm_tot(2)  ,2:ncm_tot(3)  ) &
               &   + Az(2:ncm_tot(1)  ,2:ncm_tot(2)  ,2:ncm_tot(3)  ))&
               -gstep_invq(3)*(&
               &   - Ay(1:ncm_tot(1)-1,1:ncm_tot(2)-1,1:ncm_tot(3)-1) &
               &   - Ay(2:ncm_tot(1)  ,1:ncm_tot(2)-1,1:ncm_tot(3)-1) &
               &   - Ay(1:ncm_tot(1)-1,2:ncm_tot(2)  ,1:ncm_tot(3)-1) &
               &   - Ay(2:ncm_tot(1)  ,2:ncm_tot(2)  ,1:ncm_tot(3)-1) &
               &   + Ay(1:ncm_tot(1)-1,1:ncm_tot(2)-1,2:ncm_tot(3)  ) &
               &   + Ay(2:ncm_tot(1)  ,1:ncm_tot(2)-1,2:ncm_tot(3)  ) &
               &   + Ay(1:ncm_tot(1)-1,2:ncm_tot(2)  ,2:ncm_tot(3)  ) &
               &   + Ay(2:ncm_tot(1)  ,2:ncm_tot(2)  ,2:ncm_tot(3)  ))

         Jy(1:ncm_tot(1)-1,1:ncm_tot(2)-1,1:ncm_tot(3)-1) = gstep_invq(3)*(&
               &   - Ax(1:ncm_tot(1)-1,1:ncm_tot(2)-1,1:ncm_tot(3)-1) &
               &   - Ax(2:ncm_tot(1)  ,1:ncm_tot(2)-1,1:ncm_tot(3)-1) &
               &   - Ax(1:ncm_tot(1)-1,2:ncm_tot(2)  ,1:ncm_tot(3)-1) &
               &   - Ax(2:ncm_tot(1)  ,2:ncm_tot(2)  ,1:ncm_tot(3)-1) &
               &   + Ax(1:ncm_tot(1)-1,1:ncm_tot(2)-1,2:ncm_tot(3)  ) &
               &   + Ax(2:ncm_tot(1)  ,1:ncm_tot(2)-1,2:ncm_tot(3)  ) &
               &   + Ax(1:ncm_tot(1)-1,2:ncm_tot(2)  ,2:ncm_tot(3)  ) &
               &   + Ax(2:ncm_tot(1)  ,2:ncm_tot(2)  ,2:ncm_tot(3)  )) &
               - gstep_invq(1)*(&
               &   - Az(1:ncm_tot(1)-1,1:ncm_tot(2)-1,1:ncm_tot(3)-1) &
               &   + Az(2:ncm_tot(1)  ,1:ncm_tot(2)-1,1:ncm_tot(3)-1) &
               &   - Az(1:ncm_tot(1)-1,2:ncm_tot(2)  ,1:ncm_tot(3)-1) &
               &   + Az(2:ncm_tot(1)  ,2:ncm_tot(2)  ,1:ncm_tot(3)-1) &
               &   - Az(1:ncm_tot(1)-1,1:ncm_tot(2)-1,2:ncm_tot(3)  ) &
               &   + Az(2:ncm_tot(1)  ,1:ncm_tot(2)-1,2:ncm_tot(3)  ) &
               &   - Az(1:ncm_tot(1)-1,2:ncm_tot(2)  ,2:ncm_tot(3)  ) &
               &   + Az(2:ncm_tot(1)  ,2:ncm_tot(2)  ,2:ncm_tot(3)  ))

          Jz(1:ncm_tot(1)-1,1:ncm_tot(2)-1,1:ncm_tot(3)-1) = gstep_invq(1)*(&
               &   - Ay(1:ncm_tot(1)-1,1:ncm_tot(2)-1,1:ncm_tot(3)-1) &
               &   + Ay(2:ncm_tot(1)  ,1:ncm_tot(2)-1,1:ncm_tot(3)-1) &
               &   - Ay(1:ncm_tot(1)-1,2:ncm_tot(2)  ,1:ncm_tot(3)-1) &
               &   + Ay(2:ncm_tot(1)  ,2:ncm_tot(2)  ,1:ncm_tot(3)-1) &
               &   - Ay(1:ncm_tot(1)-1,1:ncm_tot(2)-1,2:ncm_tot(3)  ) &
               &   + Ay(2:ncm_tot(1)  ,1:ncm_tot(2)-1,2:ncm_tot(3)  ) &
               &   - Ay(1:ncm_tot(1)-1,2:ncm_tot(2)  ,2:ncm_tot(3)  ) &
               &   + Ay(2:ncm_tot(1)  ,2:ncm_tot(2)  ,2:ncm_tot(3)  ))&
               - gstep_invq(2)*(&
               &   - Ax(1:ncm_tot(1)-1,1:ncm_tot(2)-1,1:ncm_tot(3)-1) &
               &   - Ax(2:ncm_tot(1)  ,1:ncm_tot(2)-1,1:ncm_tot(3)-1) &
               &   + Ax(1:ncm_tot(1)-1,2:ncm_tot(2)  ,1:ncm_tot(3)-1) &
               &   + Ax(2:ncm_tot(1)  ,2:ncm_tot(2)  ,1:ncm_tot(3)-1) &
               &   - Ax(1:ncm_tot(1)-1,1:ncm_tot(2)-1,2:ncm_tot(3)  ) &
               &   - Ax(2:ncm_tot(1)  ,1:ncm_tot(2)-1,2:ncm_tot(3)  ) &
               &   + Ax(1:ncm_tot(1)-1,2:ncm_tot(2)  ,2:ncm_tot(3)  ) &
               &   + Ax(2:ncm_tot(1)  ,2:ncm_tot(2)  ,2:ncm_tot(3)  ))
      Jtot = sqrt(Jx**2+Jy**2+Jz**2)

      write(*,*) 'Current calculation .... done '

        allocate(JcrossB_x(ncm_tot(1),ncm_tot(2),ncm_tot(3)));   JcrossB_x(:,:,:) = zero
        allocate(JcrossB_y(ncm_tot(1),ncm_tot(2),ncm_tot(3)));   JcrossB_y(:,:,:) = zero
        allocate(JcrossB_z(ncm_tot(1),ncm_tot(2),ncm_tot(3)));   JcrossB_z(:,:,:) = zero
        allocate(JcrossB_tot(ncm_tot(1),ncm_tot(2),ncm_tot(3)));  JcrossB_tot(:,:,:) = zero

   do ii=1, ncm_tot(1)-1
     do jj=1,ncm_tot(2)-1
       do kk=1,ncm_tot(3)-1
          bxi = eighth*(Ax(ii,jj,kk)+Ax(ii+1,jj,kk)+Ax(ii,jj+1,kk)+Ax(ii+1,jj+1,kk)+ &
                Ax(ii,jj,kk+1)+Ax(ii+1,jj,kk+1)+Ax(ii,jj+1,kk+1)+Ax(ii+1,jj+1,kk+1))
          byi = eighth*(Ay(ii,jj,kk)+Ay(ii+1,jj,kk)+Ay(ii,jj+1,kk)+Ay(ii+1,jj+1,kk)+ &
                Ay(ii,jj,kk+1)+Ay(ii+1,jj,kk+1)+Ay(ii,jj+1,kk+1)+Ay(ii+1,jj+1,kk+1))
          bzi = eighth*(Az(ii,jj,kk)+Az(ii+1,jj,kk)+Az(ii,jj+1,kk)+Az(ii+1,jj+1,kk)+ &
                Az(ii,jj,kk+1)+Az(ii+1,jj,kk+1)+Az(ii,jj+1,kk+1)+Az(ii+1,jj+1,kk+1))

          JcrossB_x(ii,jj,kk) =  Jy(ii,jj,kk)*bzi-Jz(ii,jj,kk)*byi
          JcrossB_y(ii,jj,kk) = -Jx(ii,jj,kk)*bzi+Jz(ii,jj,kk)*bzi
          JcrossB_z(ii,jj,kk) =  Jx(ii,jj,kk)*byi-Jy(ii,jj,kk)*bxi
       enddo
     enddo
   enddo

   JcrossB_tot = sqrt(JcrossB_x**2+JcrossB_y**2+JcrossB_z**2)


write(*,*) 'JxB calculation .... done '


!! compute grdap(Pmag)
        allocate(GradPmag_x(ncm_tot(1),ncm_tot(2),ncm_tot(3)));   GradPmag_x(:,:,:) = zero
        allocate(GradPmag_y(ncm_tot(1),ncm_tot(2),ncm_tot(3)));   GradPmag_y(:,:,:) = zero
        allocate(GradPmag_z(ncm_tot(1),ncm_tot(2),ncm_tot(3)));   GradPmag_z(:,:,:) = zero
        allocate(GradPmag_tot(ncm_tot(1),ncm_tot(2),ncm_tot(3)));   GradPmag_tot(:,:,:) = zero
       GradPmag_x(1:ncm_tot(1)-1,1:ncm_tot(2)-1,1:ncm_tot(3)-1) = gstep_invq(1)*(&
               &   - Pmag(1:ncm_tot(1)-1,1:ncm_tot(2)-1,1:ncm_tot(3)-1) &
               &   + Pmag(2:ncm_tot(1)  ,1:ncm_tot(2)-1,1:ncm_tot(3)-1) &
               &   - Pmag(1:ncm_tot(1)-1,2:ncm_tot(2)  ,1:ncm_tot(3)-1) &
               &   + Pmag(2:ncm_tot(1)  ,2:ncm_tot(2)  ,1:ncm_tot(3)-1) &
               &   - Pmag(1:ncm_tot(1)-1,1:ncm_tot(2)-1,2:ncm_tot(3)  ) &
               &   + Pmag(2:ncm_tot(1)  ,1:ncm_tot(2)-1,2:ncm_tot(3)  ) &
               &   - Pmag(1:ncm_tot(1)-1,2:ncm_tot(2)  ,2:ncm_tot(3)  ) &
               &   + Pmag(2:ncm_tot(1)  ,2:ncm_tot(2)  ,2:ncm_tot(3)  ))

       GradPmag_y(1:ncm_tot(1)-1,1:ncm_tot(2)-1,1:ncm_tot(3)-1) = gstep_invq(2)*(&
               &   - Pmag(1:ncm_tot(1)-1,1:ncm_tot(2)-1,1:ncm_tot(3)-1) &
               &   - Pmag(2:ncm_tot(1)  ,1:ncm_tot(2)-1,1:ncm_tot(3)-1) &
               &   + Pmag(1:ncm_tot(1)-1,2:ncm_tot(2)  ,1:ncm_tot(3)-1) &
               &   + Pmag(2:ncm_tot(1)  ,2:ncm_tot(2)  ,1:ncm_tot(3)-1) &
               &   - Pmag(1:ncm_tot(1)-1,1:ncm_tot(2)-1,2:ncm_tot(3)  ) &
               &   - Pmag(2:ncm_tot(1)  ,1:ncm_tot(2)-1,2:ncm_tot(3)  ) &
               &   + Pmag(1:ncm_tot(1)-1,2:ncm_tot(2)  ,2:ncm_tot(3)  ) &
               &   + Pmag(2:ncm_tot(1)  ,2:ncm_tot(2)  ,2:ncm_tot(3)  ))

        GradPmag_z(1:ncm_tot(1)-1,1:ncm_tot(2)-1,1:ncm_tot(3)-1) = gstep_invq(3)*(&
               &   - Pmag(1:ncm_tot(1)-1,1:ncm_tot(2)-1,1:ncm_tot(3)-1) &
               &   - Pmag(2:ncm_tot(1)  ,1:ncm_tot(2)-1,1:ncm_tot(3)-1) &
               &   - Pmag(1:ncm_tot(1)-1,2:ncm_tot(2)  ,1:ncm_tot(3)-1) &
               &   - Pmag(2:ncm_tot(1)  ,2:ncm_tot(2)  ,1:ncm_tot(3)-1) &
               &   + Pmag(1:ncm_tot(1)-1,1:ncm_tot(2)-1,2:ncm_tot(3)  ) &
               &   + Pmag(2:ncm_tot(1)  ,1:ncm_tot(2)-1,2:ncm_tot(3)  ) &
               &   + Pmag(1:ncm_tot(1)-1,2:ncm_tot(2)  ,2:ncm_tot(3)  ) &
               &   + Pmag(2:ncm_tot(1)  ,2:ncm_tot(2)  ,2:ncm_tot(3)  ))
write(*,*) 'Grad(Pmag) calculation .... done '


!! normalisation
Jx = sgn*nrm_J*Jx
Jy = sgn*nrm_J*Jy
Jz =     nrm_J*Jz
Jtot =   nrm_J*Jtot
!!
JcrossB_x = sgn*nrm_J*JcrossB_x*1.e-9
JcrossB_y = sgn*nrm_J*JcrossB_y*1.e-9
JcrossB_z =     nrm_J*JcrossB_z*1.e-9
JcrossB_tot = nrm_J*JcrossB_tot*1.e-9
!!
GradPmag_x = sgn*GradPmag_x*1.e-3/(c_wpi*1000.*2*mu0)
GradPmag_y = sgn*GradPmag_y*1.e-3/(c_wpi*1000.*2*mu0)
GradPmag_z =     GradPmag_z*1.e-3/(c_wpi*1000.*2*mu0)
GradPmag_tot = sqrt(GradPmag_x**2+GradPmag_y**2+GradPmag_z**2)

   ! extract plane cuts in XML
           varname1 = "Jcur_"
           ival = int(centr(3)/gstep(1))
           call extract_XY(Jtot,Jx,Jy,Jz,ival,centr,gstep,ncm_tot,&
           & varname1,run_name,radius,prefix,planet,nrm_J,sgn,coord)
           varname1 = "JxB_"
           ival = int(centr(3)/gstep(1))
           call extract_XY(JcrossB_tot,JcrossB_x,JcrossB_y,JcrossB_z,ival,centr,gstep,ncm_tot,varname1,run_name, &
           & radius,prefix,planet,nrm_J,sgn,coord)
         !  ival = int(centr(2)/gstep(2))
         !  call extract_XZ(Jtot,Jx,Jy,Jz,ival,centr,gstep,ncm_tot-1,varname,run_name,radius,prefix,planet,nrm_J,sgn,coord)
         !  ival = int(centr(1)/gstep(1))
         !  call extract_YZ(Jtot,Jx,Jy,Jz,ival,centr,gstep,ncm_tot-1,varname,run_name,radius,prefix,planet,nrm_J,sgn,coord)

 ! save cubes
 prefix_output = "Jcur"
 varname1 = "Jx";   varname2 = "Jy";   varname3 = "Jz"; varname4 = "Jtot"
 call save_field_cdf(prefix_output,run_name,Jx,Jy,Jz,Jtot,varname1,varname2,varname3,varname4,&
        & ncm_tot,centr,radius,gstep,planet,c_wpi)
 write(*,*) "Saving Current file .... done"

 prefix_output = "GPma"
 varname1 = "GPmagx";   varname2 = "GPmagy";   varname3 = "GPmagz"
 varname4 ="GPmagtot"
 call save_field_cdf(prefix_output,run_name,GradPmag_x,GradPmag_y,GradPmag_z,GradPmag_tot,varname1,varname2,&
   varname3,varname4,ncm_tot,centr,radius,gstep,planet,c_wpi)
 write(*,*) "Saving Grad(Pmag) file .... done"

 deallocate(Jx,Jy,Jz,Jtot,JcrossB_x,JcrossB_y,JcrossB_z,JcrossB_tot,Pmag,GradPmag_x,GradPmag_y,GradPmag_z)


  end subroutine compute_current

 !!############################################################
     !! diag_Ele/read_field_cdf
     !! This routines reads the global output
     !! simulation files. Concern files are :
     !!		- Magw
     !!		- Elew
     !!		- Velw
subroutine read_field_cdf(run_name,varname1,varname2,varname3,Ax,Ay,Az)

      character(len=*),intent(in) :: run_name
      character(len=*),intent(in) :: varname1,varname2,varname3
      real(dp),dimension(:,:,:),intent(inout) :: Ax,Ay,Az

      !--Local variables
      character(len=64) :: varname
      integer  :: stId,ncid
      logical :: file_e
      character(len=50) :: write_name
      character(len=64) :: filename

 __WRT_DEBUG_IN("read_field_cdf")

      !--Inquire if the file exists
      inquire( file=trim(run_name), exist=file_e )
      !--No file: return
      if(.not.(file_e)) then
       call wrtout(6,"File: "//trim(run_name)//" does not exits","PERS")
       return
      endif

      !--Open NetCDF file
      stId = nf90_open(trim(run_name), nf90_nowrite, ncid)
      call test_cdf(stId)


      call get_simple_variable_cdf(ncid,varname1 ,Ax(:,:,:) )
      call get_simple_variable_cdf(ncid,varname2 ,Ay(:,:,:) )
      call get_simple_variable_cdf(ncid,varname3 ,Az(:,:,:) )



      stId = nf90_close(ncid);  call test_cdf(stId)
      write(*,*) 'close file ',run_name
 __WRT_DEBUG_OUT("read_field_cdf")
    end subroutine read_field_cdf

 !!##############################################################
   !! IMPEX_WS/read_dim_field_cdf
   !! This routines reads a simulation file
   !! and extract the dimension of a given variabe
   !! associated to a tag
   subroutine read_dim_field_cdf(cube_name,var_name,ncm_tot,radius,centr,gstep,planet,c_wpi,B0)

   character(len=*),intent(in) :: cube_name
   character(len=*),intent(in) :: var_name
   integer,intent(inout)         :: ncm_tot(3)
   real(dp),intent(out)        :: radius,centr(3),gstep(3),c_wpi,B0
   character(len=8),intent(out) :: planet

   !--Others
   integer  :: stId,ncid
   logical :: file_e
   character(len=50) :: write_name
   integer :: varfieldId, numDims,var_id
   integer,dimension(nf90_max_var_dims) :: DimfieldId

__WRT_DEBUG_IN("read_dim_field_cdf")

   !--Inquire if the file exists
   inquire( file=trim(cube_name), exist=file_e )
   !--No file: return
   if(.not.(file_e)) then
     call wrtout(6,"File: "//trim(cube_name)//" does not exits","PERS")
      return
   endif

   !--Open NetCDF file
   stId = nf90_open(trim(cube_name), nf90_nowrite, ncid)
   call test_cdf(stId)


   !--Get number of point for fields in X,Y,Z
   StId = nf90_inq_varId(ncid,var_name,varfieldId)
   StId = nf90_Inquire_Variable(ncid,varfieldId,dimids = DimfieldId)
   StId = nf90_Inquire_Dimension(ncid,DimfieldId(1), len=ncm_tot(1))
   StId = nf90_Inquire_Dimension(ncid,DimfieldId(2), len=ncm_tot(2))
   StId = nf90_Inquire_Dimension(ncid,DimfieldId(3), len=ncm_tot(3))

    !--Get the grid step
    call get_simple_variable_cdf(ncid,"gstep",gstep)

   !--Get the obstacle radius
   call get_simple_variable_cdf(ncid,"r_planet",radius)

   !--Get the obstacle position
   call get_simple_variable_cdf(ncid,"s_centr",centr)

    !--Get the refernce length (ion inertial length)
    call get_simple_variable_cdf(ncid,"phys_length",c_wpi)
    call get_simple_variable_cdf(ncid,"phys_mag",B0)

   !--Get Planetname
   stId = nf90_inq_varid(ncid,"planetname",var_id)
   stId = nf90_get_var(ncid,var_id,planet)
  ! write(*,*) 'Planet :',planet


    stId = nf90_close(ncid);  call test_cdf(stId)
  !  write(*,*) 'close file ',cube_name

__WRT_DEBUG_OUT("read_dim_field_cdf")
   end subroutine read_dim_field_cdf

   !!##############################################################################
     !! m_IMPEX_cut/extract_XY
     !! This routine extracts data in XY plane
     !! and creates a new ascii file
     subroutine extract_XY(A,Ax,Ay,Az, cut_val,centr,gstep,ncm_tot,varname,run_name,radius,prefix,planetname,nrm,sgn,coord)
     character(len=*),intent(in) :: run_name,planetname
     character(len=*),intent(in) :: varname
     character(len=*),intent(in) :: prefix,coord
     real(dp),dimension(:,:,:),intent(in) :: A,Ax,Ay,Az
     real(dp),intent(in) :: radius,centr(3),gstep(3),nrm,sgn
     integer,intent(in)  :: ncm_tot(3),cut_val
     ! Local variable
     real(dp),dimension(:,:),allocatable :: A_XY,Ax_XY,Ay_XY,Az_XY
     character(len=50) :: write_name
     character(len=2) :: plane
     integer :: iunit,i,j
     character(len=20) :: tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9

     allocate(A_XY(ncm_tot(1),ncm_tot(2)));     A_XY(:,:) = zero
     allocate(Ax_XY(ncm_tot(1),ncm_tot(2)));    Ax_XY(:,:) = zero
     allocate(Ay_XY(ncm_tot(1),ncm_tot(2)));    Ay_XY(:,:) = zero
     allocate(Az_XY(ncm_tot(1),ncm_tot(2)));    Az_XY(:,:) = zero
     ! extract information for Ax
     A_XY(:,:) =  A(:,:,cut_val)
     Ax_XY(:,:) =  Ax(:,:,cut_val)
     Ay_XY(:,:) =  Ay(:,:,cut_val)
     Az_XY(:,:) =  Az(:,:,cut_val)

     !-- filename of XY extracted dat
     write_name = trim(varname)//'_'//'XY'//'_'//trim(run_name(1:len_trim(run_name)-7))//".xml"
     call wrtout(6," ======= Creation of file : "//trim(write_name),"PERS")
     plane = 'XY'

     print *,'sgn,radius',sgn,radius

     iunit = 1
     open(UNIT = iunit, FILE = write_name, FORM = 'FORMATTED',  &
         ACTION = 'write',STATUS = 'UNKNOWN')

     if (trim(varname)=='Jcur_') &
     & call header_2Dcut_J_VOTABLE(iunit,prefix,plane,cut_val,planetname,coord)
     if (trim(varname)=='JxB_') &
     & call header_2Dcut_JxB_VOTABLE(iunit,prefix,plane,cut_val,planetname,coord)


     do i=1,ncm_tot(1)
       do j = 1,ncm_tot(2)
            write(iunit,'(a)') '<TR>'
    call getstrf(sgn*(float(i)*gstep(1)-centr(1))/radius,tmp1,'(f9.3)')
    call getstrf(sgn*(float(j)*gstep(2)-centr(2))/radius,tmp2,'(f9.3)')
    call getstrf((float(cut_val)*gstep(3)-centr(3))/radius,tmp3,'(f9.3)')
    call  getstrf(A_XY(i,j),tmp4,'(e12.3)') ; call getstrf(Ax_XY(i,j),tmp5,'(e12.3)')
    call getstrf(Ay_XY(i,j),tmp6,'(e12.3)') ; call getstrf(Az_XY(i,j),tmp7,'(e12.3)')
            write(iunit,'(15a)') &
              '<TD>',trim(adjustl(tmp1)),'</TD><TD>',trim(adjustl(tmp2)),'</TD><TD>',trim(adjustl(tmp3)), &
           '</TD><TD>',trim(adjustl(tmp4)),'</TD><TD>',trim(adjustl(tmp5)),'</TD><TD>',trim(adjustl(tmp6)), &
           '</TD><TD>',trim(adjustl(tmp7)),'</TD>'
   !         write(iunit,'(a,f9.3,a,f9.3,a,f9.3,a,4(f12.3,a))') &
   !           '<TD>',sgn*(float(i)*gstep(1)-centr(1))/radius,'</TD><TD>',sgn*(float(j)*gstep(2)-centr(2))/radius, &
   !           '</TD><TD>',(float(cut_val)*gstep(3)-centr(3))/radius,'</TD><TD>', A_XY(i,j),'</TD><TD>',&
   !           Ax_XY(i,j),'</TD><TD>',Ay_XY(i,j),'</TD><TD>',Az_XY(i,j),'</TD>'
            write(iunit,'(a)') '</TR>'
       enddo
     enddo
     ! finalize the file
     write(iunit,'(a)') '</TABLEDATA>'
     write(iunit,'(a)') '</DATA>'
     write(iunit,'(a)') '</TABLE>'
     write(iunit,'(a)') '</RESSOURCE>'
     write(iunit,'(a)') '</VOTABLE>'

     close(iunit)


     deallocate(A_XY,Ax_XY,Ay_XY,Az_XY)
  end subroutine extract_XY

  !!##########################################################
  !! m_IMPEX_cut/header_2Dcut_J_IMPEX
  !! This routine writes the header of the xml
  !! file for the 2D cut IMPEX diag
  subroutine header_2Dcut_J_VOTABLE(iunit,prefix,plane,cutval,planetname,coord)
  integer, intent(in) :: iunit,cutval
  character(len=*),intent(in) :: prefix,plane,planetname,coord

    write(iunit,'(a)') '<?xml version="1.0"?>'
    write(iunit,'(a)') '<VOTABLE version="1.2" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" '
    write(iunit,'(a)')  'xmlns="http://www.ivoa.net/xml/VOTable/v1.2" xmlns:stc="http://www.ivoa.net/xml/STC/v1.30" >'
    write(iunit,'(a)')  '<RESSOURCE name="IMPEx 2D cut">'
    write(iunit,'(a)')  '<TABLE name="results">'


     write(iunit,'(a)')  '<GROUP ID="PosFrame" ref="MSO">'
     write(iunit,'(a)') '<PARAM name="Radius" datatype="float" arraysize="*" value="3393." unit="km" ucd="phys.size"/>'
     write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C1" ref="col1"/>'
     write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C2" ref="col2"/>'
     write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C3" ref="col3"/>'
     write(iunit,'(a)') '</GROUP>'
     write(iunit,'(a)') '<GROUP ID="FieldFrame" ref="MSO">'
     write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C1" ref="col5"/>'
     write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C2" ref="col6"/>'
     write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C3" ref="col7"/>'
     write(iunit,'(a)') '</GROUP>'




      write(iunit,'(3a)') '<DESCRIPTION> Current components </DESCRIPTION>'
      write(iunit,'(3a)') '<FIELD name="X" ID="col1" ucd="pos.cartesian.x" ', &
      '  utype="stc:AstroCoords.Position3D.Value3.C1" datatype="float" width="9"/>'
      write(iunit,'(3a)') '<FIELD name="Y" ID="col2" ucd="pos.cartesian.y" ', &
      '  utype="stc:AstroCoords.Position3D.Value3.C2" datatype="float" width="9"/>'
      write(iunit,'(3a)') '<FIELD name="Z" ID="col3" ucd="pos.cartesian.z" ', &
      '  utype="stc:AstroCoords.Position3D.Value3.C3" datatype="float" width="9"/>'
      write(iunit,'(3a)') '<FIELD name="Jtot" ID="col4" ucd="phys.current" ',&
      '  utype="" datatype="float" width="12" unit="nA.m-2"/>'
      write(iunit,'(3a)') '<FIELD name="Jx" ID="col5" ucd="phys.current" ',&
      '  utype="" datatype="float" width="12" unit="nA.m-2"/>'
      write(iunit,'(3a)') '<FIELD name="Jy" ID="col6" ucd="phys.current" ',&
      '  utype="" datatype="float" width="12" unit="nA.m-2"/>'
      write(iunit,'(3a)') '<FIELD name="Jz" ID="col7" ucd="phys.current"',&
      '  utype="" datatype="float" width="12" unit="nA.m-2"/>'


   write(iunit,'(a)') '<DATA>'
    write(iunit,'(a)') '<TABLEDATA>'


end subroutine header_2Dcut_J_VOTABLE

 !!##########################################################
  !! m_IMPEX_cut/header_2Dcut_J_IMPEX
  !! This routine writes the header of the xml
  !! file for the 2D cut IMPEX diag
  subroutine header_2Dcut_JxB_VOTABLE(iunit,prefix,plane,cutval,planetname,coord)
  integer, intent(in) :: iunit,cutval
  character(len=*),intent(in) :: prefix,plane,planetname,coord

    write(iunit,'(a)') '<?xml version="1.0"?>'
    write(iunit,'(a)') '<VOTABLE version="1.2" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" '
    write(iunit,'(a)')  'xmlns="http://www.ivoa.net/xml/VOTable/v1.2" xmlns:stc="http://www.ivoa.net/xml/STC/v1.30" >'
    write(iunit,'(a)')  '<RESSOURCE name="IMPEx 2D cut">'
    write(iunit,'(a)')  '<TABLE name="results">'


     write(iunit,'(a)')  '<GROUP ID="PosFrame" ref="MSO">'
     write(iunit,'(a)') '<PARAM name="Radius" datatype="float" arraysize="*" value="3393." unit="km" ucd="phys.size"/>'
     write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C1" ref="col1"/>'
     write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C2" ref="col2"/>'
     write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C3" ref="col3"/>'
     write(iunit,'(a)') '</GROUP>'
     write(iunit,'(a)') '<GROUP ID="FieldFrame" ref="MSO">'
     write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C1" ref="col5"/>'
     write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C2" ref="col6"/>'
     write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C3" ref="col7"/>'
     write(iunit,'(a)') '</GROUP>'




      write(iunit,'(3a)') '<DESCRIPTION> Current components </DESCRIPTION>'
      write(iunit,'(3a)') '<FIELD name="X" ID="col1" ucd="pos.cartesian.x" ', &
      '  utype="stc:AstroCoords.Position3D.Value3.C1" datatype="float" width="9"/>'
      write(iunit,'(3a)') '<FIELD name="Y" ID="col2" ucd="pos.cartesian.y" ', &
      '  utype="stc:AstroCoords.Position3D.Value3.C2" datatype="float" width="9"/>'
      write(iunit,'(3a)') '<FIELD name="Z" ID="col3" ucd="pos.cartesian.z" ', &
      '  utype="stc:AstroCoords.Position3D.Value3.C3" datatype="float" width="9"/>'
      write(iunit,'(3a)') '<FIELD name="JxB_tot" ID="col4" ucd="phys.current" ',&
      '  utype="" datatype="float" width="12" unit="nA.m-3"/>'
      write(iunit,'(3a)') '<FIELD name="JxB_x" ID="col5" ucd="phys.current" ',&
      '  utype="" datatype="float" width="12" unit="nN.m-3"/>'
      write(iunit,'(3a)') '<FIELD name="JxB_y" ID="col6" ucd="phys.current" ',&
      '  utype="" datatype="float" width="12" unit="nN.m-3"/>'
      write(iunit,'(3a)') '<FIELD name="JxB_z" ID="col7" ucd="phys.current"',&
      '  utype="" datatype="float" width="12" unit="nN.m-3"/>'


   write(iunit,'(a)') '<DATA>'
    write(iunit,'(a)') '<TABLEDATA>'


end subroutine header_2Dcut_JxB_VOTABLE

subroutine getstrf(inp,outp,frmt)
 real(dp),intent(in) :: inp
 character(len=*),intent(in) :: frmt
 character(len=20),intent(out) :: outp
 write(outp,frmt) inp
end subroutine getstrf



!!############################################################
!!=============================================================
!!subroutine: diag_Elec/save_field_cdf
!! FUNCTION
!!  Create the name of the file containing fields information
!! INPUT
!!  filwrt=suffix containgt data
!! OUTPUT
!!  file :  file where fields will be recorded


subroutine save_field_cdf(prefix_output,run_name,Ax,Ay,Az,Atot,varname1,varname2,varname3,varname4,&
        & ncm_tot,centr,radius,gstep,planet,c_wpi)
character(len=*), intent(in) :: prefix_output, run_name, varname1,varname2,varname3,varname4
real(dp),dimension(:,:,:),intent(in) :: Ax, Ay, Az, Atot
character(len=30) :: write_name
character(len=*),intent(in) :: planet
real(dp),intent(in) :: c_wpi
integer :: ncid, stId,ii
integer,dimension(3),intent(in) :: ncm_tot
real(dp),dimension(3),intent(in) :: centr,gstep
real(dp),intent(in) :: radius
integer :: dimid(7), varid(75)


  !--Creation du fichier de diagnostic de champ a acces sequentiel
      write_name = prefix_output//'_'//trim(run_name)//".nc"
      call wrtout(6," ======= Creation of file : "//trim(write_name),"PERS")

      !--Open the file for write Particles
      stId = nf90_create(trim(write_name),nf90_64BIT_OFFSET, ncid)
      call test_cdf(stId)

    call set_global_attribute_cdf(ncid,"Fields")

    !--This is needed by m_wrt_common_cdf to use global number of points
      !and not relative to processus
  !--Define the dimensions that will define the size of the file.
  stId = nf90_def_dim(ncid, "dim_scalar"      , 1         , dimid(1))
  call test_cdf(stId)
  stId = nf90_def_dim(ncid, "bidimensional"   , 2         , dimid(2))
  call test_cdf(stId)
  stId = nf90_def_dim(ncid, "space_dimension" , 3         , dimid(3))
  call test_cdf(stId)
  stId = nf90_def_dim(ncid, "size_x"          , ncm_tot(1)    , dimid(4))
  call test_cdf(stId)
  stId = nf90_def_dim(ncid, "size_y"          , ncm_tot(2)    , dimid(5))
  call test_cdf(stId)
  stId = nf90_def_dim(ncid, "size_z"          , ncm_tot(3)    , dimid(6))
  call test_cdf(stId)
  stId = nf90_def_dim(ncid, "dimplanet",20, dimid(7))
  call test_cdf(stId)

      ii = 1
  !    call common_def_var_cdf(ncid,varid,dimid,ii)
  !    write(*,*)  "defining variable dimension .... ok"

    stId = nf90_def_var(ncid, "gstep",QP_NF90_DP,dimid(3), varid(ii))
  call test_cdf(stId); ii = ii+1
      stId = nf90_def_var(ncid, "s_centr",QP_NF90_DP,dimid(3), varid(ii))
  call test_cdf(stId); ii = ii+1
    stId = nf90_def_var(ncid, "r_planet",QP_NF90_DP,dimid(3), varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "planetname",nf90_char,dimid(7), varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "phys_length",QP_NF90_DP,dimid(1), varid(ii))
  call test_cdf(stId); ii = ii+1


      stId = nf90_def_var(ncid, "s_min",   QP_NF90_DP,dimid(3), varid(ii))
      call test_cdf(stId); ii = ii+1
      stId = nf90_def_var(ncid, "s_max",   QP_NF90_DP,dimid(3), varid(ii))
      call test_cdf(stId); ii = ii+1
      write(*,*)  "defining box dimension .... ok"

 	!do ish=1,n_spe
                    stId = nf90_def_var(ncid,varname1, QP_NF90_DP,dimid(4:6),varid(ii))
            call test_cdf(stId); ii = ii+1
                    stId = nf90_def_var(ncid,varname2, QP_NF90_DP,dimid(4:6),varid(ii))
            call test_cdf(stId); ii = ii+1
                    stId = nf90_def_var(ncid,varname3, QP_NF90_DP,dimid(4:6),varid(ii))
            call test_cdf(stId); ii = ii+1
                    stId = nf90_def_var(ncid,varname4, QP_NF90_DP,dimid(4:6),varid(ii))
            call test_cdf(stId); ii = ii+1

  	!enddo
       ii = ii-1

       write(*,*)  "swithcing to writing mode"
      !--Switch to write mode
        stId = nf90_enddef(ncid); call test_cdf(stId)

        ii = 1
        !--Write common variables into the file
  !      call common_put_var_cdf(ncid,varid,ii)
    stId = nf90_put_var(ncid, varid(ii), gstep)
  call test_cdf(stId); ii = ii+1
    stId = nf90_put_var(ncid, varid(ii), centr)
  call test_cdf(stId); ii = ii+1
      stId = nf90_put_var(ncid, varid(ii), radius)
  call test_cdf(stId); ii = ii+1
        stId = nf90_put_var(ncid, varid(ii), planet)
  call test_cdf(stId); ii = ii+1
        stId = nf90_put_var(ncid, varid(ii), c_wpi)
  call test_cdf(stId); ii = ii+1

        stId = nf90_put_var(ncid, varid(ii), s_min)
        call test_cdf(stId); ii = ii+1
        stId = nf90_put_var(ncid, varid(ii), s_max)
        call test_cdf(stId); ii = ii+1

 	!do ish=1,n_spe
            stId = nf90_put_var(ncid, varid(ii), Ax(:,:,:))
            call test_cdf(stId); ii = ii+1
            stId = nf90_put_var(ncid, varid(ii), Ay(:,:,:))
            call test_cdf(stId); ii = ii+1
            stId = nf90_put_var(ncid, varid(ii), Az(:,:,:))
            call test_cdf(stId); ii = ii+1
            stId = nf90_put_var(ncid, varid(ii), Atot(:,:,:))
            call test_cdf(stId); ii = ii+1

  	!enddo
         ii = ii-1

        !--Close the file
         stId = nf90_close(ncid); call test_cdf(stId)

   end subroutine save_field_cdf

end module m_current_computation
