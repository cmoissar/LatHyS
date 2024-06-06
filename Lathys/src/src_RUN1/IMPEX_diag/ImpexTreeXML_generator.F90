!!===================================================
!!===================================================
module IMPEXTreeXML_generator


 use defs_basis
 use defs_variable
 use defs_grid
 use defs_arr3Dtype
 use defs_particletype
 use defs_species
 use m_VO
 use m_writeout
 use m_logo,only: take_date

#ifdef HAVE_NETCDF 
 use netcdf
 use defs_basic_cdf
 use diag_wrt_common_cdf
#endif
#include "q-p_common.h"

 implicit none
 private 

 public ::                    &
      write_numdat_XML,       &
      write_granule_TS_XML,   &  
      write_granule_3D_XML,   &  
      write_granule_2D_XML

contains
 !********************************************************************
 ! Auteur				:	 RModolo, SHess
 ! Date					:	 05/04/12
 ! Institution				:	LATMOS/CNRS/IPSL
 ! Derniere modification		:	10/04/12	
 ! Resume	
 ! 
 !********************************************************************

subroutine write_numdat_XML(fildat,planetname,datatype,ts,c3d,c2d,sp,ln,spcraft,angle,YYYY,MM,DD,hh,minut,ss,numdatid,xminmso,xmaxmso,yminmso,ymaxmso,zminmso,zmaxmso,orient)
  character(len=8),intent(in) :: planetname
  character(len=*),intent(in) ::datatype,spcraft,fildat
  character(len=600),intent(out) :: numdatid
  real(dp),intent(in) :: angle,xminmso,xmaxmso,yminmso,ymaxmso,zminmso,zmaxmso
  logical,intent(in) :: ts,sp,ln,c3d,c2d
  integer,dimension(:),intent(in) :: YYYY,MM,DD,hh,minut
  real(dp),dimension(:),intent(in) :: ss,orient(3)

  character(len=10)::chem,atm,ppmn,ppcs,ptype
  character(len=200)::string_name,string_name2,string_name4,string_name5
  character(len=600)::string,write_name
  character(len=10) :: date
  character(len=5) ::msg,msg_angle
  integer ::tmp_unit,valeur1(8),sz
  character(len=8) ::planetnam
  character(len=16) ::PID
  
  select case (datatype)
    case ("Opl")
      chem="O  " ; atm="8" ; ppmn="16" ; ppcs="1";ptype="Ion";PID="Planetary_O+    ";
    case ("O2p")
      chem="O2 " ; atm=" " ; ppmn="32" ; ppcs="1";ptype="Ion";PID="Planetary_O2+   ";
    case ("CO2")
      chem="CO2" ; atm=" " ; ppmn="44" ; ppcs="1";ptype="Ion";PID="Planetary_CO2+  ";
    case ("Hpl")
      chem="H  " ; atm="1" ; ppmn="1" ; ppcs="1";ptype="Proton";PID="Planetary_H+    ";
    case ("Hsw")
      chem="H  " ; atm="1" ; ppmn="1" ; ppcs="1";ptype="Proton";PID="SolarWind_H+    ";
    case ("Hes")
      chem="He " ; atm="2" ; ppmn="4" ; ppcs="2";ptype="Ion";PID="SolarWind_He+   ";
    case default
    PID = trim(datatype)
end select

    if (angle /= -1) then
    write(msg_angle,'(f5.1)')angle
    else
    msg_angle =""
    endif

planetnam=planetname
print *,planetnam,ichar(planetnam(1:1))
if (ichar(planetnam(1:1)).gt.97) planetnam(1:1)=char(ichar(planetnam(1:1))-32)
print *,planetnam,ichar(planetnam(1:1))

write(string_name2, "(3a)") trim(adjustl(fildat)),".",trim(adjustl(datatype))
if (ts.or.sp.or.c2d) write(string_name2, "(3a)") trim(adjustl(string_name2)),"_",trim(adjustl(spcraft)) ! for c2D use spcraft for plane XY,YZ,XZ
if (ts.or.sp) write(string_name2, "(3a)") trim(adjustl(string_name2)),"_",trim(adjustl(msg_angle))
if (ts.or.sp) write(string_name5, "(6a)") trim(adjustl(PID))," @",trim(adjustl(spcraft))," (",trim(adjustl(msg_angle)),")"
if (c3d) write(string_name5, "(2a)") trim(adjustl(PID)),"/3D"
if (c3d) write(string_name2, "(2a)") trim(adjustl(string_name2)),"_3D"
if (c2d) write(string_name5, "(3a)") trim(adjustl(PID)),"_",trim(adjustl(spcraft))
if (c2d) write(string_name2, "(2a)") trim(adjustl(string_name2)),"_2D"
if (ln) write(string_name5, "(2a)") trim(adjustl(PID)),"/FL"
if (ln) write(string_name2, "(2a)") trim(adjustl(string_name2)),"_FL"

string_name4=string_name2

    write(write_name,'(2a)') trim(string_name4),".numdata.xml"
    write(*,*) 'Saving file'  ,trim(write_name)

  valeur1 = take_date()


  write(msg,'(i4)')valeur1(1)
  write(date,'(a2,a2,a1,i2.2,a1,i2.2)')&
       & '20',trim(adjustl(msg(3:))),'-',&!--Year
       & valeur1(2),'-',&   !--Month
       & valeur1(3)   !--Day

    ! in VOTale format
    tmp_unit = 1001
     open(UNIT = tmp_unit,FILE = write_name, FORM = 'FORMATTED', STATUS = 'UNKNOWN', &
    	ACTION = 'WRITE')

write(tmp_unit,'(A)') "<?xml version='1.0' encoding='UTF-8'?>"
write(tmp_unit,'(A)') '<Spase  xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"'
write(tmp_unit,'(A)') '      xmlns="http://www.impex.latmos.ipsl.fr"'
write(tmp_unit,'(A)') '      xsi:schemaLocation="http://hess.page.latmos.ipsl.fr/impex/impex-0_9.xsd">'
write(tmp_unit,'(A)') "<Version>2.2.2</Version>"
write(tmp_unit,'(A)') "<NumericalOutput>"
write(string_name, "(3a)")  trim(planetnam),"_",trim(adjustl(fildat))!SimRun
write(string_name2, "(3a)") trim(adjustl(string_name)),"/",trim(adjustl(datatype))
if (c2d) write(string_name2, "(2a)") trim(adjustl(string_name2)),"/2D"
if (ts.or.sp.or.c2d) write(string_name2, "(3a)") trim(adjustl(string_name2)),"/",trim(adjustl(spcraft))
if (ts.or.sp) write(string_name2, "(3a)") trim(adjustl(string_name2)),"/",trim(adjustl(msg_angle))
if (c3d) write(string_name2, "(2a)") trim(adjustl(string_name2)),"/3D"
if (ln) write(string_name2, "(2a)") trim(adjustl(string_name2)),"/FL"
numdatid="impex://LATMOS/Hybrid/"//trim(adjustl(string_name2))

write(string, "(3a)") "<ResourceID>",trim(adjustl(numdatid)),"</ResourceID>"
write(tmp_unit,'(A)') trim(adjustl(string))
write(tmp_unit,'(A)') "<ResourceHeader>"
write(string,'(3a)') "<ResourceName>",trim(string_name5),"</ResourceName>"
write(tmp_unit,'(A)') trim(adjustl(string))
write(string,'(3a)') "<ReleaseDate>",trim(adjustl(date)),"T00:00:00.000</ReleaseDate>"
write(tmp_unit,'(A)') trim(adjustl(string))
write(tmp_unit,'(A)') "<Description />"
write(tmp_unit,'(A)') "<Contact>"
write(tmp_unit,'(A)') "<PersonID>LATMOS</PersonID>"
write(tmp_unit,'(A)') "<Role>DataProducer</Role>"
write(tmp_unit,'(A)') "</Contact>"
write(tmp_unit,'(A)') "</ResourceHeader>"
write(tmp_unit,'(A)') "<AccessInformation>"
write(tmp_unit,'(A)') "<RepositoryID>impex://LATMOS</RepositoryID>"
write(tmp_unit,'(A)') "<AccessURL>"
write(tmp_unit,'(A)') "<URL>ftp://latmos.ipsl.fr/IMPEX/</URL>"
write(tmp_unit,'(A)') "</AccessURL>"
if (c3d) then 
	write(tmp_unit,'(A)') "<Format>NetCDF</Format>" 
else 
	write(tmp_unit,'(A)') "<Format>XML</Format>"
endif
write(tmp_unit,'(A)') "</AccessInformation>" 
  select case (datatype)
  case ("Mag")
    write(tmp_unit,'(A)') "<MeasurementType>MagneticField</MeasurementType>"
  case ("Ele")
    write(tmp_unit,'(A)') "<MeasurementType>ElectricField</MeasurementType>"
  case ("The")
    write(tmp_unit,'(A)') "<MeasurementType>ThermalPlasma</MeasurementType>"
  case ("Vel")
    write(tmp_unit,'(A)') "<MeasurementType>ThermalPlasma</MeasurementType>"
  case ("Den")
    write(tmp_unit,'(A)') "<MeasurementType>ThermalPlasma</MeasurementType>"
  case ("Jcu")
    write(tmp_unit,'(A)') "<MeasurementType>Current</MeasurementType>"
  case default
    if (sp) then 
	write(tmp_unit,'(A)') "<MeasurementType>Spectrum</MeasurementType>"
    else
	write(tmp_unit,'(A)') "<MeasurementType>IonComposition</MeasurementType>"
    endif
  end select
if (ts.or.sp) then 
 write(tmp_unit,'(A)') "<TemporalDescription>"
 write(tmp_unit,'(A)') "<TimeSpan>"
 sz=size(YYYY)
 write(string,'(a,i4,5(a,i2.2),a)') "<StartDate>",YYYY(1),"-",MM(1),"-",DD(1),"T",hh(1),":",minut(1),":",int(ss(1)),".000</StartDate>"
 write(tmp_unit,'(A)') trim(adjustl(string))
 write(string,'(a,i4,5(a,i2.2),a)') "<StopDate>",YYYY(sz),"-",MM(sz),"-",DD(sz),"T",hh(sz),":",minut(sz),":",int(ss(sz)),".000</StopDate>"
 write(tmp_unit,'(A)') trim(adjustl(string))
 write(tmp_unit,'(A)') "</TimeSpan>"
 !<Cadence>PT1800S</Cadence>
 write(tmp_unit,'(A)') "</TemporalDescription>"
endif
if (c3d.or.c2d.or.ln) then 
 write(tmp_unit,'(A)') "<SpatialDescription>"
 if (c3d) write(tmp_unit,'(A)') "<Dimension>3</Dimension>"
 if (c2d) write(tmp_unit,'(A)') "<Dimension>2</Dimension>"
 if (ln) write(tmp_unit,'(A)') "<Dimension>1</Dimension>"
write(tmp_unit,'(A)') "<CoordinateSystem>"
write(tmp_unit,'(A)') "<CoordinateRepresentation>Cartesian</CoordinateRepresentation>"
write(tmp_unit,'(A)') "<CoordinateSystemName>MSO</CoordinateSystemName>"
write(tmp_unit,'(A)') "</CoordinateSystem>"
write(tmp_unit,'(A)') "<Units>km</Units>"
write(tmp_unit,'(A)') "<UnitsConversion>1000 > m</UnitsConversion>"
if (c3d) then
write(string,'(a,3(f10.1,a),a)') "<RegionBegin>",xminmso," ",yminmso," ",zminmso," ","</RegionBegin>"
write(tmp_unit,'(A)') trim(adjustl(string))
write(string,'(a,3(f10.1,a),a)') "<RegionEnd>",xmaxmso," ",ymaxmso," ",zmaxmso," ","</RegionEnd>"
write(tmp_unit,'(A)') trim(adjustl(string))
endif
if (c2d) then 
write(string,'(a,3(f5.2,a),a)') "<PlaneNormalVector>",orient(1)," ",orient(2)," ",orient(3)," ","</PlaneNormalVector>"
write(tmp_unit,'(A)') trim(adjustl(string))
write(string,'(a,3(f5.2,a),a)') "<PlanePoint>0 0 0</PlanePoint>"
write(tmp_unit,'(A)') trim(adjustl(string))
endif
write(tmp_unit,'(A)') "</SpatialDescription>"
endif
write(string,'(3a)') "<SimulatedRegion>",trim(planetnam),"</SimulatedRegion>"
write(tmp_unit,'(A)') trim(adjustl(string))
write(string,'(3a)') "<InputResourceID>impex://LATMOS/Hybrid/",trim(adjustl(string_name)),"/SimRun</InputResourceID>"
write(tmp_unit,'(A)') trim(adjustl(string))



if (ts.or.c3d.or.c2d) then 
  select case (datatype)
  case ("Jcu")
	write(tmp_unit,'(A)') "<Parameter>"
	write(tmp_unit,'(A)') "<Name>TotalCurrent</Name>"
	write(tmp_unit,'(A)') "<ParameterKey>Jtot</ParameterKey>"
	write(tmp_unit,'(A)') "<Units>mA.km-2</Units>"
	write(tmp_unit,'(A)') "<Field>"
	write(tmp_unit,'(A)') "<Qualifier>Total</Qualifier>"
	write(tmp_unit,'(A)') "<FieldQuantity>Current</FieldQuantity>"
	write(tmp_unit,'(A)') "</Field>"
	write(tmp_unit,'(A)') "</Parameter>"
	write(tmp_unit,'(A)') "<Parameter>"
	write(tmp_unit,'(A)') "<Name>Current</Name>"
	write(tmp_unit,'(A)') "<ParameterKey>Jx Jy Jz</ParameterKey>"
	write(tmp_unit,'(A)') "<Units>mA.km-2</Units>"
	write(tmp_unit,'(A)') "<UnitsConverion>1.e-3 > A.m^-2 </UnitsConversion>"
	write(tmp_unit,'(A)') "<Field>"
	write(tmp_unit,'(A)') "<Qualifier>Vector</Qualifier>"
	write(tmp_unit,'(A)') "<FieldQuantity>Current</FieldQuantity>"
	write(tmp_unit,'(A)') "</Field>"
	write(tmp_unit,'(A)') "</Parameter>"
  case ("Mag")
        if (ts.or.c2d) then
	write(tmp_unit,'(A)') "<Parameter>"
	write(tmp_unit,'(A)') "<Name>TotalMagneticField</Name>"
	write(tmp_unit,'(A)') "<ParameterKey>Btot</ParameterKey>"
	write(tmp_unit,'(A)') "<Units>nT</Units>"
	write(tmp_unit,'(A)') "<UnitsConversion> 1.e-9 > T </UnitsConversion>"
	write(tmp_unit,'(A)') "<Field>"
	write(tmp_unit,'(A)') "<Qualifier>Total</Qualifier>"
	write(tmp_unit,'(A)') "<FieldQuantity>Magnetic</FieldQuantity>"
	write(tmp_unit,'(A)') "</Field>"
	write(tmp_unit,'(A)') "</Parameter>"
	endif
	write(tmp_unit,'(A)') "<Parameter>"
	write(tmp_unit,'(A)') "<Name>MagneticField</Name>"
	write(tmp_unit,'(A)') "<ParameterKey>Bx By Bz</ParameterKey>"
	write(tmp_unit,'(A)') "<Units>nT</Units>"
	write(tmp_unit,'(A)') "<UnitsConversion> 1.e-9 > T </UnitsConversion>"	
	write(tmp_unit,'(A)') "<Field>"
	write(tmp_unit,'(A)') "<Qualifier>Vector</Qualifier>"
	write(tmp_unit,'(A)') "<FieldQuantity>Magnetic</FieldQuantity>"
	write(tmp_unit,'(A)') "</Field>"
	write(tmp_unit,'(A)') "</Parameter>"
  case ("Ele")
        if (ts.or.c2d) then    
	write(tmp_unit,'(A)') "<Parameter>"
	write(tmp_unit,'(A)') "<Name>TotalElectricField</Name>"
	write(tmp_unit,'(A)') "<ParameterKey>Etot</ParameterKey>"
	write(tmp_unit,'(A)') "<Units>mV.m-1</Units>"
	write(tmp_unit,'(A)') "<UnitsConversion> 1.e-3 > V.m^-1 </UnitsConversion>"	
	write(tmp_unit,'(A)') "<Field>"
	write(tmp_unit,'(A)') "<Qualifier>Total</Qualifier>"
	write(tmp_unit,'(A)') "<FieldQuantity>Electric</FieldQuantity>"
	write(tmp_unit,'(A)') "</Field>"
	write(tmp_unit,'(A)') "</Parameter>"
	endif
	write(tmp_unit,'(A)') "<Parameter>"
	write(tmp_unit,'(A)') "<Name>ElectricField</Name>"
	write(tmp_unit,'(A)') "<ParameterKey>Ex Ey Ez</ParameterKey>"
	write(tmp_unit,'(A)') "<Units>mV.m-1</Units>"
	write(tmp_unit,'(A)') "<UnitsConversion> 1.e-3 > V.m^-1 </UnitsConversion>"		
	write(tmp_unit,'(A)') "<Field>"
	write(tmp_unit,'(A)') "<Qualifier>Vector</Qualifier>"
	write(tmp_unit,'(A)') "<FieldQuantity>Electric</FieldQuantity>"
	write(tmp_unit,'(A)') "</Field>"
	write(tmp_unit,'(A)') "</Parameter>"
  case ("The")
	write(tmp_unit,'(A)') "<Parameter>"
	write(tmp_unit,'(A)') "<Name>n_e</Name>"
	write(tmp_unit,'(A)') "<ParameterKey>Density</ParameterKey>"
	write(tmp_unit,'(A)') "<Units>cm^-3</Units>"
	write(tmp_unit,'(A)') "<UnitsConversion> 1.e6 > m^-3 </UnitsConversion>"		
	write(tmp_unit,'(A)') "<Particle>"
	write(tmp_unit,'(A)') "<ParticleType>Electron</ParticleType>"
	write(tmp_unit,'(A)') "<Qualifier>Total</Qualifier>"
	write(tmp_unit,'(A)') "<ParticleQuantity>NumberDensity</ParticleQuantity>"
	write(tmp_unit,'(A)') "</Particle>"
	write(tmp_unit,'(A)') "</Parameter>"
        if (ts.or.c2d) then
	write(tmp_unit,'(A)') "<Parameter>"
	write(tmp_unit,'(A)') "<Name>|U|</Name>"
	write(tmp_unit,'(A)') "<ParameterKey>Utot</ParameterKey>"
	write(tmp_unit,'(A)') "<Units> km.s^-1 </Units>"		
        write(tmp_unit,'(A)') "<UnitsConversion> 1.e3 > m.s^-1 </UnitsConversion>"			
	write(tmp_unit,'(A)') "<Particle>"
	write(tmp_unit,'(A)') "<ParticleType>Ion</ParticleType>"
	write(tmp_unit,'(A)') "<Qualifier>Total</Qualifier>"
	write(tmp_unit,'(A)') "<ParticleQuantity>Velocity</ParticleQuantity>"
	write(tmp_unit,'(A)') "</Particle>"
	write(tmp_unit,'(A)') "</Parameter>"
	endif
	write(tmp_unit,'(A)') "<Parameter>"
	write(tmp_unit,'(A)') "<Name>U</Name>"
	write(tmp_unit,'(A)') "<ParameterKey>Ux Uy Uz</ParameterKey>"
	write(tmp_unit,'(A)') "<Units>km.s^-1</Units>"
	write(tmp_unit,'(A)') "<UnitsConversion> 1.e3 > m.s^-1 </UnitsConversion>"	
	write(tmp_unit,'(A)') "<Particle>"
	write(tmp_unit,'(A)') "<ParticleType>Ion</ParticleType>"
	write(tmp_unit,'(A)') "<Qualifier>Vector</Qualifier>"
	write(tmp_unit,'(A)') "<ParticleQuantity>Velocity</ParticleQuantity>"
	write(tmp_unit,'(A)') "</Particle>"
	write(tmp_unit,'(A)') "</Parameter>"
	write(tmp_unit,'(A)') "<Parameter>"
	write(tmp_unit,'(A)') "<Name>T_e</Name>"
	write(tmp_unit,'(A)') "<ParameterKey>Temperature</ParameterKey>"
	write(tmp_unit,'(A)') "<Units>eV</Units>"
	write(tmp_unit,'(A)') "<Particle>"
	write(tmp_unit,'(A)') "<ParticleType>Electron</ParticleType>"
	write(tmp_unit,'(A)') "<ParticleQuantity>Temperature</ParticleQuantity>"
	write(tmp_unit,'(A)') "</Particle>"
	write(tmp_unit,'(A)') "</Parameter>"

  case default
	write(tmp_unit,'(A)') "<Parameter>"
        write(tmp_unit,'(A)') "<Name>n_i</Name>"
	write(tmp_unit,'(A)') "<ParameterKey>Density</ParameterKey>"
	write(tmp_unit,'(A)') "<Units>cm^-3</Units>"
	write(tmp_unit,'(A)') "<UnitsConversion> 1.e6 > m^-3 </UnitsConversion>"		
	write(tmp_unit,'(A)') "<Particle>"
	write(string,'(3a)') "<PopulationID>",trim(PID),"</PopulationID>"
	write(tmp_unit,'(A)') trim(adjustl(string))
	write(string,'(3a)') "<ParticleType>",trim(ptype),"</ParticleType>"
	write(tmp_unit,'(A)') trim(adjustl(string))
	write(tmp_unit,'(A)') "<ParticleQuantity>NumberDensity</ParticleQuantity>"
	write(string,'(3a)') "<ChemicalFormula>",trim(chem),"</ChemicalFormula>"
	write(tmp_unit,'(A)') trim(adjustl(string))
	if (trim(atm).ne." ") then
	write(string,'(3a)') "<AtomicNumber>",trim(atm),"</AtomicNumber>"	
 	write(tmp_unit,'(A)') trim(adjustl(string))
	endif
	write(string,'(3a)') "<PopulationMassNumber>",trim(ppmn),"</PopulationMassNumber>"
	write(tmp_unit,'(A)') trim(adjustl(string))
	write(string,'(3a)') "<PopulationChargeState>",trim(ppcs),"</PopulationChargeState>"
	write(tmp_unit,'(A)') trim(adjustl(string))
	write(tmp_unit,'(A)') "</Particle>"
	write(tmp_unit,'(A)') "</Parameter>"
	if (ts.or.c2d) then
	write(tmp_unit,'(A)') "<Parameter>"
        write(tmp_unit,'(A)') "<Name>|U_i|</Name>"
	write(tmp_unit,'(A)') "<ParameterKey>Utot</ParameterKey>"
	write(tmp_unit,'(A)') "<Units>km.s^-1</Units>"
	write(tmp_unit,'(A)') "<Particle>"
	write(string,'(3a)') "<PopulationID>",trim(PID),"</PopulationID>"
	write(tmp_unit,'(A)') trim(adjustl(string))
	write(string,'(3a)') "<ParticleType>",trim(ptype),"</ParticleType>"
	write(tmp_unit,'(A)') trim(adjustl(string))
	write(tmp_unit,'(A)') "<Qualifier>Total</Qualifier>"
	write(tmp_unit,'(A)') "<ParticleQuantity>Velocity</ParticleQuantity>"
	write(string,'(3a)') "<ChemicalFormula>",trim(chem),"</ChemicalFormula>"
	write(tmp_unit,'(A)') trim(adjustl(string))
	if (trim(atm).ne." ") then
	write(string,'(3a)') "<AtomicNumber>",trim(atm),"</AtomicNumber>"	
 	write(tmp_unit,'(A)') trim(adjustl(string))
	endif
	write(string,'(3a)') "<PopulationMassNumber>",trim(ppmn),"</PopulationMassNumber>"
	write(tmp_unit,'(A)') trim(adjustl(string))
	write(string,'(3a)') "<PopulationChargeState>",trim(ppcs),"</PopulationChargeState>"
	write(tmp_unit,'(A)') trim(adjustl(string))
	write(tmp_unit,'(A)') "</Particle>"
	write(tmp_unit,'(A)') "</Parameter>"
        endif
	write(tmp_unit,'(A)') "<Parameter>"
	write(tmp_unit,'(A)') "<Name>U_i</Name>"
	write(tmp_unit,'(A)') "<ParameterKey>Ux Uy Uz</ParameterKey>"
	write(tmp_unit,'(A)') "<Units>km.s^-1</Units>"
        write(tmp_unit,'(A)') "<UnitsConversion> 1.e3 > m.s^-1 </UnitsConversion>"		
	write(tmp_unit,'(A)') "<Particle>"
	write(string,'(3a)') "<PopulationID>",trim(PID),"</PopulationID>"
	write(tmp_unit,'(A)') trim(adjustl(string))
	write(string,'(3a)') "<ParticleType>",trim(ptype),"</ParticleType>"
	write(tmp_unit,'(A)') trim(adjustl(string))
	write(tmp_unit,'(A)') "<Qualifier>Vector</Qualifier>"
	write(tmp_unit,'(A)') "<ParticleQuantity>Velocity</ParticleQuantity>"
	write(string,'(3a)') "<ChemicalFormula>",trim(chem),"</ChemicalFormula>"
	write(tmp_unit,'(A)') trim(adjustl(string))
	if (trim(atm).ne." ") then
	write(string,'(3a)') "<AtomicNumber>",trim(atm),"</AtomicNumber>"	
 	write(tmp_unit,'(A)') trim(adjustl(string))
	endif
	write(string,'(3a)') "<PopulationMassNumber>",trim(ppmn),"</PopulationMassNumber>"
	write(tmp_unit,'(A)') trim(adjustl(string))
	write(string,'(3a)') "<PopulationChargeState>",trim(ppcs),"</PopulationChargeState>"
	write(tmp_unit,'(A)') trim(adjustl(string))
	write(tmp_unit,'(A)') "</Particle>"
	write(tmp_unit,'(A)') "</Parameter>"
	write(tmp_unit,'(A)') "<Parameter>"
	write(tmp_unit,'(A)') "<Name>T_i</Name>"
	write(tmp_unit,'(A)') "<ParameterKey>Temperature</ParameterKey>"
	write(tmp_unit,'(A)') "<Units>eV</Units>"
	write(tmp_unit,'(A)') "<Particle>"
	write(string,'(3a)') "<PopulationID>",trim(PID),"</PopulationID>"
	write(tmp_unit,'(A)') trim(adjustl(string))
	write(string,'(3a)') "<ParticleType>",trim(ptype),"</ParticleType>"
	write(tmp_unit,'(A)') trim(adjustl(string))
	write(tmp_unit,'(A)') "<ParticleQuantity>Temperature</ParticleQuantity>"
	write(string,'(3a)') "<ChemicalFormula>",trim(chem),"</ChemicalFormula>"
	write(tmp_unit,'(A)') trim(adjustl(string))
	if (trim(atm).ne."") then
        write(string,'(3a)') "<AtomicNumber>",trim(atm),"</AtomicNumber>"	
 	write(tmp_unit,'(A)') trim(adjustl(string))
	endif
	write(string,'(3a)') "<PopulationMassNumber>",trim(ppmn),"</PopulationMassNumber>"
	write(tmp_unit,'(A)') trim(adjustl(string))
	write(string,'(3a)') "<PopulationChargeState>",trim(ppcs),"</PopulationChargeState>"
	write(tmp_unit,'(A)') trim(adjustl(string))
	write(tmp_unit,'(A)') "</Particle>"
	write(tmp_unit,'(A)') "</Parameter>"
  end select
endif
if (sp) then 
	write(tmp_unit,'(A)') "<Parameter>"
	write(tmp_unit,'(A)') "<Name>EnergySpectrum</Name>"
	write(tmp_unit,'(A)') "<ParameterKey>Spectrum</ParameterKey>"
	write(tmp_unit,'(A)') "<Particle>"
	write(tmp_unit,'(A)') "<ParticleType>Ion</ParticleType>"
	write(tmp_unit,'(A)') "<ParticleQuantity>NumberFlux</ParticleQuantity>"
	write(tmp_unit,'(A)') "</Particle>"
	write(tmp_unit,'(A)') "</Parameter>"
endif

if (c3d) write(tmp_unit,'(A)') "<SimulationProduct>3DCubes</SimulationProduct>"
if (c2d) write(tmp_unit,'(A)') "<SimulationProduct>2DCuts</SimulationProduct>"
if (ts) write(tmp_unit,'(A)') "<SimulationProduct>TimeSeries</SimulationProduct>"
if (sp) write(tmp_unit,'(A)') "<SimulationProduct>Spectra</SimulationProduct>"
if (ln) write(tmp_unit,'(A)') "<SimulationProduct>Lines</SimulationProduct>"
if (sp.or.ts) then
write(tmp_unit,'(A)') "<Property>"
	write(tmp_unit,'(A)') "<Name>Spacecraft</Name>"
	write(tmp_unit,'(A)') "<PropertyQuantity>Platform</PropertyQuantity>"
	write(string,'(3a)') "<PropertyValue>",trim(adjustl(spcraft)),"</PropertyValue>" 
	write(tmp_unit,'(A)') trim(adjustl(string))
	write(tmp_unit,'(A)') "</Property>" 
	write(tmp_unit,'(A)') "<Property>" 
	write(tmp_unit,'(A)') "<Name>IMF Clock Angle</Name>" 
	write(tmp_unit,'(A)') "<PropertyQuantity>IMFClockAngle</PropertyQuantity>" 
	write(tmp_unit,'(A)') "<Units>degrees</Units>" 
        write(string,'(3a)') "<PropertyValue>",trim(adjustl(msg_angle)),"</PropertyValue>" 
	write(tmp_unit,'(A)') trim(adjustl(string))
	write(tmp_unit,'(A)') "</Property>"
endif
	
write(tmp_unit,'(A)') "</NumericalOutput>"

write(tmp_unit,'(A)') "</Spase>"

 close(tmp_unit)

end subroutine write_numdat_XML


subroutine write_granule_TS_XML(traj_name,prefix,run_name,YYYY,MM,DD,hh,minut,ss,numdatid,angle)
  character(len=*),intent(in) ::traj_name,prefix,run_name,numdatid
  integer,dimension(:),intent(in) :: YYYY,MM,DD,hh,minut
  real(dp),dimension(:),intent(in) :: ss
  real(dp),intent(in) :: angle

  character(len=200)::write_name,write_name2
  character(len=600)::string
  character(len=10) :: date
  character(len=4) ::msg
  integer ::tmp_unit,valeur1(8),sz
  character(len=5) ::msg_angle

     if (angle /= -1) then
     write(msg_angle,'(f5.1)')angle
     else
     msg_angle =""
    endif

    write(write_name2,'(a3,a1,a,a1,a,a,a)')trim(prefix),"_",trim(run_name(1:len_trim(run_name)-7)),"_",trim(traj_name(1:len_trim(traj_name)-4)),"_",trim(adjustl(msg_angle))

    write(write_name,'(2a)')trim(write_name2),".granule.xml"
    write(*,*) 'Saving file'  ,trim(write_name)

  valeur1 = take_date()

  write(msg,'(i4)')valeur1(1)
  write(date,'(a2,a2,a1,i2.2,a1,i2.2)')&
       & '20',trim(adjustl(msg(3:))),'-',&!--Year
       & valeur1(2),'-',&   !--Month
       & valeur1(3)   !-Day

    ! in VOTale format
    tmp_unit = 1002
     open(UNIT = tmp_unit,FILE = write_name, FORM = 'FORMATTED', STATUS = 'UNKNOWN', &
    	ACTION = 'WRITE')

write(tmp_unit,'(A)') "<?xml version='1.0' encoding='UTF-8'?>"
write(tmp_unit,'(A)') '<Spase  xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"'
write(tmp_unit,'(A)') '      xmlns="http://www.impex.latmos.ipsl.fr"'
write(tmp_unit,'(A)') '      xsi:schemaLocation="http://hess.page.latmos.ipsl.fr/impex/impex-0_6.xsd">'
write(tmp_unit,'(A)') "<Version>2.2.2</Version>"
write(tmp_unit,'(A)') "<Granule>"
write(string,'(5a)')  "<ResourceID>",trim(adjustl(numdatid)),"/",trim(adjustl(write_name2)),"</ResourceID>"
write(tmp_unit,'(A)') trim(adjustl(string))
write(string,'(3a)') "<ReleaseDate>",trim(adjustl(date)),"T00:00:00.000</ReleaseDate>"
write(tmp_unit,'(A)') trim(adjustl(string))
write(string,'(3a)') "<ParentID>",trim(adjustl(numdatid)),"</ParentID>"
write(tmp_unit,'(A)') trim(adjustl(string))
sz=size(YYYY)
write(string,'(a,i4,5(a,i2.2),a)') "<StartDate>",YYYY(1),"-",MM(1),"-",DD(1),"T",hh(1),":",minut(1),":",int(ss(1)),".000</StartDate>"
write(tmp_unit,'(A)') trim(adjustl(string))
write(string,'(a,i4,5(a,i2.2),a)') "<StopDate>",YYYY(sz),"-",MM(sz),"-",DD(sz),"T",hh(sz),":",minut(sz),":",int(ss(sz)),".000</StopDate>"
write(tmp_unit,'(A)') trim(adjustl(string))
write(tmp_unit,'(A)') "<Source>"
write(tmp_unit,'(A)') "<SourceType>Data</SourceType>"
write(string,'(3a)') "<URL>",trim(adjustl(write_name2)),".xml</URL>"
write(tmp_unit,'(A)') trim(adjustl(string))
write(tmp_unit,'(A)') "</Source>"
write(tmp_unit,'(A)') "</Granule>"
write(tmp_unit,'(A)') "</Spase>"
 close(tmp_unit)
end subroutine write_granule_TS_XML

subroutine write_granule_3D_XML(prefix,run_name,xminmso,xmaxmso,yminmso,ymaxmso,zminmso,zmaxmso,numdatid)
  character(len=*),intent(in) ::prefix,run_name,numdatid
  real(dp),intent(in) :: xminmso,xmaxmso,yminmso,ymaxmso,zminmso,zmaxmso

  character(len=200)::write_name,write_name2
  character(len=600)::string
  character(len=10) :: date
  character(len=4) ::msg
  integer ::tmp_unit,valeur1(8),sz

    write(write_name2,'(3a)')trim(prefix),"_",trim(run_name)

    write(write_name,'(2a)')trim(adjustl(write_name2(1:len_trim(write_name2)-7))),"_3DCube.granule.xml"
    write(*,*) 'Saving file'  ,write_name     

  valeur1 = take_date()

  write(msg,'(i4)')valeur1(1)
  write(date,'(a2,a2,a1,i2.2,a1,i2.2)')&
       & '20',trim(adjustl(msg(3:))),'-',&!--Year
       & valeur1(2),'-',&   !--Month
       & valeur1(3)   !-Day

    ! in VOTale format
    tmp_unit = 1002
     open(UNIT = tmp_unit,FILE = write_name, FORM = 'FORMATTED', STATUS = 'UNKNOWN', &
    	ACTION = 'WRITE')

write(tmp_unit,'(A)') "<?xml version='1.0' encoding='UTF-8'?>"
write(tmp_unit,'(A)') '<Spase  xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"'
write(tmp_unit,'(A)') '      xmlns="http://www.impex.latmos.ipsl.fr"'
write(tmp_unit,'(A)') '      xsi:schemaLocation="http://hess.page.latmos.ipsl.fr/impex/impex-0_6.xsd">'
write(tmp_unit,'(A)') "<Version>2.2.2</Version>"
write(tmp_unit,'(A)') "<Granule>"
write(string,'(5a)')  "<ResourceID>",trim(adjustl(numdatid)),"/",trim(adjustl(write_name2)),"</ResourceID>"
write(tmp_unit,'(A)') trim(adjustl(string))
write(string,'(3a)') "<ReleaseDate>",trim(adjustl(date)),"T00:00:00.000</ReleaseDate>"
write(tmp_unit,'(A)') trim(adjustl(string))
write(string,'(3a)') "<ParentID>",trim(adjustl(numdatid)),"</ParentID>"
write(tmp_unit,'(A)') trim(adjustl(string))
write(string,'(a,3(f10.1,a),a)') "<RegionBegin>",xminmso," ",yminmso," ",zminmso," ","</RegionBegin>"
write(tmp_unit,'(A)') trim(adjustl(string))
write(string,'(a,3(f10.1,a),a)') "<RegionEnd>",xmaxmso," ",ymaxmso," ",zmaxmso," ","</RegionEnd>"
write(tmp_unit,'(A)') trim(adjustl(string))
write(tmp_unit,'(A)') "<Source>"
write(tmp_unit,'(A)') "<SourceType>Data</SourceType>"
write(string,'(3a)') "<URL>",trim(adjustl(write_name2)),".nc</URL>"
write(tmp_unit,'(A)') trim(adjustl(string))
write(tmp_unit,'(A)') "</Source>"
write(tmp_unit,'(A)') "</Granule>"
write(tmp_unit,'(A)') "</Spase>"
 close(tmp_unit)
end subroutine write_granule_3D_XML

subroutine write_granule_2D_XML(prefix,run_name,xminmso,xmaxmso,yminmso,ymaxmso,zminmso,zmaxmso,numdatid)
  character(len=*),intent(in) ::prefix,run_name,numdatid
  real(dp),intent(in) :: xminmso,xmaxmso,yminmso,ymaxmso,zminmso,zmaxmso

  character(len=200)::write_name,write_name2
  character(len=600)::string
  character(len=10) :: date
  character(len=4) ::msg
  integer ::tmp_unit,valeur1(8),sz

    write(write_name2,'(3a)')trim(prefix),"_",trim(run_name(1:len_trim(run_name)-7))

    write(write_name,'(2a)')trim(adjustl(write_name2)),"_2DCut.granule.xml"
    write(*,*) 'Saving file'  ,write_name     

  valeur1 = take_date()

  write(msg,'(i4)')valeur1(1)
  write(date,'(a2,a2,a1,i2.2,a1,i2.2)')&
       & '20',trim(adjustl(msg(3:))),'-',&!--Year
       & valeur1(2),'-',&   !--Month
       & valeur1(3)   !-Day

    ! in VOTale format
    tmp_unit = 1002
     open(UNIT = tmp_unit,FILE = write_name, FORM = 'FORMATTED', STATUS = 'UNKNOWN', &
    	ACTION = 'WRITE')

write(tmp_unit,'(A)') "<?xml version='1.0' encoding='UTF-8'?>"
write(tmp_unit,'(A)') '<Spase  xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"'
write(tmp_unit,'(A)') '      xmlns="http://www.impex.latmos.ipsl.fr"'
write(tmp_unit,'(A)') '      xsi:schemaLocation="http://hess.page.latmos.ipsl.fr/impex/impex-0_6.xsd">'
write(tmp_unit,'(A)') "<Version>2.2.2</Version>"
write(tmp_unit,'(A)') "<Granule>"
write(string,'(5a)')  "<ResourceID>",trim(adjustl(numdatid)),"/",trim(adjustl(write_name2)),"</ResourceID>"
write(tmp_unit,'(A)') trim(adjustl(string))
write(string,'(3a)') "<ReleaseDate>",trim(adjustl(date)),"T00:00:00.000</ReleaseDate>"
write(tmp_unit,'(A)') trim(adjustl(string))
write(string,'(3a)') "<ParentID>",trim(adjustl(numdatid)),"</ParentID>"
write(tmp_unit,'(A)') trim(adjustl(string))
write(string,'(a,3(f10.1,a),a)') "<RegionBegin>",xminmso," ",yminmso," ",zminmso," ","</RegionBegin>"
write(tmp_unit,'(A)') trim(adjustl(string))
write(string,'(a,3(f10.1,a),a)') "<RegionEnd>",xmaxmso," ",ymaxmso," ",zmaxmso," ","</RegionEnd>"
write(tmp_unit,'(A)') trim(adjustl(string))
write(tmp_unit,'(A)') "<Source>"
write(tmp_unit,'(A)') "<SourceType>Data</SourceType>"
write(string,'(3a)') "<URL>",trim(adjustl(write_name2)),".xml</URL>"
write(tmp_unit,'(A)') trim(adjustl(string))
write(tmp_unit,'(A)') "</Source>"
write(tmp_unit,'(A)') "</Granule>"
write(tmp_unit,'(A)') "</Spase>"
 close(tmp_unit)
end subroutine write_granule_2D_XML

end module IMPEXTreeXML_generator

