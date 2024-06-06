!!===================================================
!!===================================================
module m_IMPEX_spectro_orbit


 use defs_basis
 use defs_variable
 use defs_grid
 use defs_arr3Dtype
 use defs_particletype
 use defs_species
 use defs_parametre,only : fildat
 use m_writeout
 use m_VO
 use ImpexTreeXML_generator
#ifdef HAVE_NETCDF 
 use netcdf
 use defs_basic_cdf
 use diag_wrt_common_cdf
#endif
#include "q-p_common.h"

 implicit none
 private 

#ifdef HAVE_NETCDF
 private ::		  &
      create_file_name,	  &
      save_spectro_orbit
      !extract_XY,	  &
      !extract_XZ
      !create_file_diag,   &
      !read_moment_species_cdf,&
      !read_field_cdf_1_array

 public ::                &
      !read_part_cdf,	&
      extract_spectro
contains
 !********************************************************************
 ! Auteur				:	 RModolo, SHess
 ! Date					:	 05/04/12
 ! Institution				:	LATMOS/CNRS/IPSL
 ! Derniere modification		:	10/04/12	
 ! Resume	
 ! 
 !********************************************************************
 
 !!##################################################################
 !! m_IMPEX_spectro_orbit/extract_spectro
 !! This routine reads the particles 'p3' files
 !! and extract ion spectrogram along 
 !! S/C virtual orbit
 subroutine extract_spectro(traj_name,run_name,X0_traj,Y0_traj,Z0_traj, &
    		YYYY,MM,DD,hh,minut,ss,angle,spacecraft)
 character(len=50),intent(in) :: traj_name
 character(len=20),intent(in) :: run_name
 character(len=*), intent(in) :: spacecraft
 real(dp),dimension(:),intent(inout) :: X0_traj,Y0_traj,Z0_traj
 integer,dimension(:),intent(in) :: YYYY,MM,DD,hh,minut
 real(dp),dimension(:),intent(in) :: ss   
 real(dp),intent(in) :: angle
! Local variables
 character(len=40) :: filename
 character(len=8) :: planet
 integer :: rank,StId,ncid,var_id
 real(dp) :: centr(3),gstep(3),radius,Vspeed,c_wpi,box_size(6),orient(3)
 integer :: npm,nb_procs,n,nptot
 type(particletype),dimension(:),allocatable :: particule
 ! Spectrogram energy table
  integer :: nb_ene=94	! number of energy channel
  real(dp),dimension(nb_ene) :: tab_E 
  real(dp) :: r_cyl,r_cyl_km
  real(dp),dimension(:,:),allocatable :: Ion_spectra
  integer :: nb_traj,t,i,min_t,compt,index_E,sgn
    character(len=5) :: coord 
  real(dp) :: min_dist_to_traj,traj(3),distance
  real(dp) :: mass,weight,speed_part,energy
  character(len=600)::numdatid 
   real(dp),dimension(size(Y0_traj)) :: X_traj,Y_traj,Z_traj
  
 !  X_traj=X0_traj
 !  Y_traj=Y0_traj*cos((angle+90)/180.*pi)-Z0_traj*sin((angle+90)/180.*pi)
 !  Z_traj=Z0_traj*cos((angle+90)/180.*pi)+Y0_traj*sin((angle+90)/180.*pi)
 X_traj = X0_traj
 Y_traj = Y0_traj
 Z_traj = Z0_traj
 
 
   print *,' Max & Min location of the trajectory after rotation'
   print *, 'angle',angle
   print *, 'X MSO',maxval(X_traj),minval(X_traj)
   print *, 'Y MSO',maxval(Y_traj),minval(Y_traj)
   print *, 'Z MSO',maxval(Z_traj),minval(Z_traj)
      
   !-- Energy table (comes from Mars-Express) in eV
    tab_E = (/1., 2. , 2.8, 3.8, 4.9, 6., 7.4, 9., 10.6, 12., 14., 15.7, 17.9, & 
    	20., 22.7, 25.4,  28.4 , 31.6, 35., 39., 42., 47., 52., 57., 63. , 69., &
    	76., 83., 93., 104., 116., 127., 138., 150., 172., 184., 207., 218., 241., &
    	263., 286., 309., 331., 365., 400., 443., 467., 513., 558., 603., 660., &
    	716., 773., 841., 921., 1000., 1079., 1181., 1283., 1385., 1510., 1646., &
    	1782., 1941., 2099., 2280., 2484., 2700., 2927., 3176., 3754., 4071., & 
    	4423., 4808., 5216., 5669., 6156., 6689., 7256., 7879., 8959., 9295., &
    	10100., 10962., 11903., 12923.,14033., 15246., 16549., 17977., 19518., &
    	21196., 23009., 24992./)   
    	
   !-- Get the length of the trajectory file
   nb_traj = size(X_traj,1)
   
   !-- Allocate Ion_spectra array
   allocate(Ion_spectra(nb_traj,nb_ene));	Ion_spectra = zero
   
   compt = 0
   
    rank=0
    call create_file_name(filename,run_name,rank)
    print *,'Reading file ',filename
 
    stId = nf90_open(trim(filename), nf90_nowrite, ncid)
    call test_cdf(stId) 
    
 
    !--Get the maximum number of particle per process
     call get_simple_variable_cdf(ncid,"npm",npm)  

     !--Get the grid step
     call get_simple_variable_cdf(ncid,"gstep",gstep)
 
    !--Get the obstacle radius
    call get_simple_variable_cdf(ncid,"r_planet",radius)
    
    !--Get the obstacle position
    call get_simple_variable_cdf(ncid,"s_centr",centr)  
    
    !--Get the number of processus which generated the file
     call get_simple_variable_cdf(ncid,"nproc",nb_procs)   
     
    !--Get the refernce speed (Alfven speed)
     call get_simple_variable_cdf(ncid,"phys_speed",Vspeed)  
   ! print *,'Alfven speed',Vspeed 
    
    !--Get the refernce length (ion inertial length)
    call get_simple_variable_cdf(ncid,"phys_length",c_wpi)  
   ! print *,'Ion inertia length',c_wpi

    stId = nf90_inq_varid(ncid,"planetname",var_id)
    stId = nf90_get_var(ncid,var_id,planet)

     stId = nf90_close(ncid);  call test_cdf(stId)
    write(*,*) 'close file ',filename   

   ! coordinate transformation for each object
   select case (trim(planet))
   case ("mars","mercury")
     sgn = -1
     coord='MSO'
   case ("ganymede")
     sgn=+1
     coord = 'GPhiO'
   case("titan")
     sgn = +1
     coord='TIIS'
     case default 
     planet="mars    "
     sgn = -1
     coord = 'MSO'
    end select
   
   !--Convert MSO trajectory in simulation unit and coordinate
   X_traj = (sgn*X_traj*radius + centr(1))
   Y_traj = (sgn*Y_traj*radius + centr(2))
   Z_traj =  (Z_traj*radius + centr(3))
   
   
    !-- Radius of cylinder over which the spectrogram is constructed
    r_cyl_km = 1000.	! in km
    r_cyl = r_cyl_km/(c_wpi)
    
    
    !-- Allocate particle arrays
    allocate(particule(npm))	
    
    
    write(*,*) ' '
    write(*,*) '---------------------------'
    write(*,*) ' Ion Spectrogram diagnostic'
    write(*,*) '---------------------------'
    write(*,*)
    
   print *, 'X simu',maxval(X_traj),minval(X_traj)
   print *, 'Y simu',maxval(Y_traj),minval(Y_traj)
   print *, 'Z simu',maxval(Z_traj),minval(Z_traj)
    
   
    !-- Open the diferent "p3_" files and get the 
    !-- partiles information
    do n=0,nb_procs-1
      !-- Set to zero all particle values
      call set_zero_particle(particule,1,npm)
       !-- prepare the name of the file to read
       call create_file_name(filename,run_name,n) 
       
       
       !-- Open the file
       stId = nf90_open(trim(filename), nf90_nowrite, ncid)
       call test_cdf(stId)        
       print *,'Treating file :',trim(filename)
      
      !--Get the number of particle in the current process file
       call get_simple_variable_cdf(ncid,"nptot",nptot)  
      
       call get_var_particle_cdf(particule,nptot,ncid)
       
       stId = nf90_close(ncid);  call test_cdf(stId)
       
!       print *,' Max & Min position for x',maxval(particule(:nptot)%pos(1)),minval(particule(:nptot)%pos(1))
!       print *,' Max & Min position for y',maxval(particule(:nptot)%pos(2)),minval(particule(:nptot)%pos(2))
!       print *,' Max & Min position for z',maxval(particule(:nptot)%pos(3)),minval(particule(:nptot)%pos(3))
       
       !-- For a given position in the trajectory we collect
       !-- all particles in a tube centred on S/C orbit
        do t=1,nb_traj
  !      print *,'time ',t
        traj(1) = X_traj(t); traj(2) = Y_traj(t);	traj(3) = Z_traj(t)
        !-- we sweep on all particles
          do i = 1,nptot
          distance = sqrt(dot_product(particule(i)%pos-traj,particule(i)%pos-traj))
          !-- If the particle is in the cylinder
          if (distance <= r_cyl) then
            compt = compt + 1
            mass = particule(i)%mass/particule(i)%char
            if (mass==1._dp) then ! H+ composition
            if ((particule(i)%orig == 0).and.(particule(i)%exc == 0)) then
            weight = particule(i)%char
            if (mass == 2) weight = weight/2.
            speed_part = dot_product(particule(i)%vel,particule(i)%vel)
            energy = 0.5*amu_pmass*mass*speed_part*(Vspeed*1.e3)**2/e_Cb
            index_E = minloc(abs(tab_E(:)-energy),dim=1)
        !    print *,'Particle with energy [eV] index',energy, index_E,mass
            Ion_spectra(t,index_E) = Ion_spectra(t,index_E) + weight
            endif ! on comosition
           endif
          endif          
          enddo
        enddo
       
!       do i = 1,nptot
!       !-- For a given particle, try to find the minimum distance to the trajectory
!       !-- trajectory is in the simulation coordinate system and in simulation unit
!        traj(1) = X_traj(1); traj(2) = Y_traj(1);	traj(3) = Z_traj(1)
!        min_dist_to_traj = sqrt(dot_product(particule(i)%pos-centr,traj))
!        min_t = 1
!        do t = 2,nb_traj
!          traj(1) = X_traj(t); traj(2) = Y_traj(t);	traj(3) = Z_traj(t)
!          distance = sqrt(dot_product(particule(i)%pos-centr,traj))
!          if (distance < min_dist_to_traj) then
!            min_t = t
!            min_dist_to_traj = distance
!          endif  
!        enddo
!        !-- If the particle is in the cylinder
!        if (min_dist_to_traj <= r_cyl) then
!          compt = compt + 1
!          mass = particule(i)%mass/particule(i)%char
!          weight = particule(i)%char
!          if (mass == 2) weight = weight/2.
!          speed_part = dot_product(particule(i)%vel,particule(i)%vel)
!          energy = 0.5*amu_pmass*mass*speed_part*(Vspeed*1.e3)**2/e_Cb
!          index_E = minloc(abs(tab_E(:)-energy),dim=1)
!          print *,'Particle with energy [eV] index',energy, index_E,mass
!          Ion_spectra(min_t,index_E) = Ion_spectra(min_t,index_E) + weight
!        endif
!      enddo  
       
    enddo
    
    !-- Number of particle in the cylinder around S/C orbit
    print *,'# of particle counted ',compt

   
   call save_spectro_orbit(traj_name,run_name,YYYY,MM,DD,hh,minut,ss, &
   		Ion_spectra,tab_E,nb_ene,planet,coord,angle)    
    
    
   box_size=(/ -centr(1),ncm_tot(1)*gstep(1)-centr(1),-centr(2),ncm_tot(2)*gstep(2)-centr(2),-centr(3),ncm_tot(3)*gstep(3)-centr(3) /)
   box_size=box_size*c_wpi
   !--Write tree.xml metadata
   call write_numdat_XML(fildat,planet,'Spectro',.false.,.false.,.false.,.true.,.false.,trim(spacecraft),angle,&
	YYYY,MM,DD,hh,minut,ss,numdatid,box_size(1),box_size(2),box_size(3),box_size(4),box_size(5),box_size(6),orient)
   call write_granule_TS_XML(traj_name,'Spectro',run_name,YYYY,MM,DD,hh,minut,ss,numdatid,angle)


   
    deallocate(particule,Ion_spectra)
    		
 end subroutine extract_spectro
 
  !!=============================================================
  !!subroutine: m_IMPEX_spectro_orbit/create_file_name
  !! FUNCTION 
  !!  Create the name of the file containing particles information
  !! INPUT
  !!  filwrt=suffix containgt data
  !!  me=processus identity
  !! OUTPUT
  !!  name_file=name of the file where particles will be recorded
  subroutine create_file_name(name_file,filwrt,me)
 
   integer,intent(in) :: me
   character(len=40),intent(out) :: name_file 
   character(len=*),intent(in) :: filwrt 
   
   write(name_file,'(a3,i4.4,a1,a)')"p3_",me,'_',trim(filwrt)
 #ifdef HAVE_NETCDF
   name_file = trim(name_file)//".nc"
 #endif
 
 end subroutine create_file_name
 
    !!###################################################
    !! m_IMPEX_diag/save_spectro_orbit
    !! This routines dumps into a file
    !! the extracted spectrogram along a trajectory
    subroutine save_spectro_orbit(traj_name,run_name,YYYY,MM,DD,hh,minut,ss, &
   		Ion_spectra,tab_E,nb_ene,planetname,coord,angle)
    character(len=*),intent(in) ::traj_name,run_name,planetname,coord
    integer,dimension(:),intent(in) :: YYYY,MM,DD,hh,minut
    real(dp),dimension(:),intent(in) :: ss,tab_E
    integer,intent(in) :: nb_ene
    real(dp),intent(in) :: Ion_spectra(:,:)
    real(dp),intent(in) :: angle
    ! Local variables
    character(len=50) :: write_name,prefix
    character(len=23)::votime
    integer :: iunit,i,nb
    character(len=5) ::msg_angle

    write(msg_angle,'(f5.1)')angle
   
    prefix = 'Spectro'
    !--Create file name for orbital value 
    write(write_name,'(a3,a1,a,a1,a,a,a,a)')trim(prefix),"_",trim(run_name(1:len_trim(run_name)-7)),"_",trim(traj_name(1:len_trim(traj_name)-4)),"_",trim(adjustl(msg_angle)),".xml"
    write(*,*) 'Saving file'  ,write_name     
    
    !-- Get length of the array
    nb = size(YYYY,1)

    ! in VOTale format
    iunit = 1
     open(UNIT = iunit,FILE = write_name, FORM = 'FORMATTED', STATUS = 'UNKNOWN', &
    	ACTION = 'WRITE')
    call header_spectro_VOTABLE(iunit,prefix,planetname,coord)
    
     do i=1,nb
          write(iunit,'(a)') '<TR>' 
    call get_VOTIME(YYYY(i),MM(i),DD(i),hh(i),minut(i),ss(i),votime)
          write(iunit,'(3a,94(f11.2),a)') '<TD>',votime,'</TD> <TD>',Ion_spectra(i,:),'</TD>'
          write(iunit,'(a)') '</TR>'	
    enddo
  ! finalize the file
  write(iunit,'(a)') '</TABLEDATA>'
  write(iunit,'(a)') '</DATA>'
  write(iunit,'(a)') '</TABLE>'
  write(iunit,'(a)') '</RESSOURCE>'
  write(iunit,'(a)') '</VOTABLE>'      
    close(iunit)

    
!    iunit = 1
!    open(UNIT = iunit,FILE = write_name, FORM = 'FORMATTED', STATUS = 'UNKNOWN', &
!    	ACTION = 'WRITE')
!    do i=1,nb
!      write(iunit,'(i4,4(a1,i2),a1,f8.3,a1,94(f11.2))') YYYY(i)," ",MM(i)," ",DD(i)," ",hh(i)," ",minut(i), &
!      	" ",ss(i)," ",Ion_spectra(i,:)
!    enddo	
!    close(iunit)

    end subroutine save_spectro_orbit 
     !!######################################################
    !! m_IMPEX_diag/header_orbit_VOTABLE
    !! This routine write the header of the 
    !! VOTABLE file
    subroutine header_spectro_VOTABLE(iunit,prefix,planetname,coord)
    character(len=*),intent(in) :: prefix,planetname,coord
    integer,intent(in) :: iunit
  write(iunit,'(a)') '<?xml version="1.0"?>'
  write(iunit,'(a)') '<VOTABLE version="1.2" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" '
  write(iunit,'(a)')  'xmlns="http://www.ivoa.net/xml/VOTable/v1.2" xmlns:stc="http://www.ivoa.net/xml/STC/v1.30" >'  
  write(iunit,'(a)')  '<RESSOURCE name="IMPEx fly-through">'
  write(iunit,'(a)')  '<TABLE name="results">'
  
  select case (trim(planetname))
  case("mars","mercury")
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
  case("ganymede")
   write(iunit,'(a)')  '<GROUP ID="PosFrame" ref="GPHIO">'
   write(iunit,'(a)') '<PARAM name="Radius" datatype="float" arraysize="*" value="2634." unit="km" ucd="phys.size"/>'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C1" ref="col1"/>'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C2" ref="col2"/>'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C3" ref="col3"/>'
   write(iunit,'(a)') '</GROUP>'
   write(iunit,'(a)') '<GROUP ID="FieldFrame" ref="MSO">'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C1" ref="col5"/>'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C2" ref="col6"/>'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C3" ref="col7"/>'
   write(iunit,'(a)') '</GROUP>'  
  case("titan")
   write(iunit,'(a)')  '<GROUP ID="PosFrame" ref="MSO">'
   write(iunit,'(a)') '<PARAM name="Radius" datatype="float" arraysize="*" value="2575." unit="km" ucd="phys.size"/>'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C1" ref="col1"/>'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C2" ref="col2"/>'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C3" ref="col3"/>'
   write(iunit,'(a)') '</GROUP>'
   write(iunit,'(a)') '<GROUP ID="FieldFrame" ref="MSO">'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C1" ref="col5"/>'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C2" ref="col6"/>'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C3" ref="col7"/>'
   write(iunit,'(a)') '</GROUP>'   
  end select 
  
 
    write(iunit,'(3a)') '<DESCRIPTION> Magnetic field components </DESCRIPTION>'
    write(iunit,'(3a)') '<FIELD name="Time" ID="col1" ucd="time" ref="UTC" ', &
    	'  utype="stc:AstroCoords.Time.TimeInstant.ISOTime" datatype="integer" width="4" />' 
    write(iunit,'(3a)') '<FIELD name="Counts" ID="col2" ucd="" ref=" ',trim(coord),&
        ' " utype="" datatype="float" width="10" unit="" arraysize="94"/>'        

  
  write(iunit,'(a)') '<DATA>'
  write(iunit,'(a)') '<TABLEDATA>'  
    end subroutine header_spectro_VOTABLE

 #endif
 
 end module m_IMPEX_spectro_orbit
