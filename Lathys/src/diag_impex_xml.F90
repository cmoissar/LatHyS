!!=============================================================
!!=============================================================
module diag_impex_xml

 use defs_basis
 use defs_parametre
 use defs_variable
 use defs_grid
 use environment
 use m_logo,only: take_date

#include "q-p_common.h"

 implicit none
 private

 public ::                       &
      write_impex_xml

contains
 !!############################################################

 subroutine write_impex_xml
  
  character(len=200)::string_name,string1
  character(len=600)::string
  character(len=10) :: date
  character(len=4) ::msg
  real(dp) :: valid_min(3),valid_max(3),grid_cell(3),box_size(3)
  real(dp) :: time_duration,step_duration
  integer ::tmp_unit,valeur1(8),i
  real(dp) :: rphi,rpsi,mag_val(3)
  real(dp) ::spe1_charge,spe1_mass,dens_sw0

tmp_unit=1001


  valeur1 = take_date()

  write(msg,'(i4)')valeur1(1)
  write(date,'(a2,a2,a1,i2.2,a1,i2.2)')&
       & '20',trim(adjustl(msg(3:))),'-',&!--Year
       & valeur1(2),'-',&   !--Month
       & valeur1(3)   !--Day

 open(UNIT   = tmp_unit, &
       FILE   = fildat//'.input.xml', &
       ACTION = 'write', &
       STATUS = 'unknown', &
       ACCESS = 'sequential')
write(tmp_unit,'(A)') "<?xml version='1.0' encoding='UTF-8'?>"
write(tmp_unit,'(A)') '<Spase  xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"'
write(tmp_unit,'(A)') '      xmlns="http://www.impex.latmos.ipsl.fr"'
write(tmp_unit,'(A)') '      xsi:schemaLocation="http://impex.latmos.ipsl.fr/doc/impex+spase_latest.xsd">'
write(tmp_unit,'(A)') "<Version>2.2.2</Version>"
write(tmp_unit,'(A)') "<SimulationRun>"
write(string_name, "(3a)")  trim(adjustl(planetname)),"_",trim(adjustl(fildat))
write(string, "(3a)") "<ResourceID>impex://LATMOS/Hybrid/",trim(adjustl(string_name)),"/SimRun</ResourceID>"
write(tmp_unit,'(A)') trim(adjustl(string))
write(tmp_unit,'(A)') "<Description />"
write(tmp_unit,'(A)') "<ResourceHeader>"
write(string,'(3a)') "<ResourceName>LatHyS_",trim(adjustl(string_name)),"</ResourceName>"
write(tmp_unit,'(A)') trim(adjustl(string))
write(string,'(3a)') "<ReleaseDate>",trim(adjustl(date)),"T00:00:00.000</ReleaseDate>"
write(tmp_unit,'(A)') trim(adjustl(string))
write(tmp_unit,'(A)') "<Contact>"
write(tmp_unit,'(A)') "<PersonID>LATMOS</PersonID>"
write(tmp_unit,'(A)') "<Role>DataProducer</Role>"
write(tmp_unit,'(A)') "</Contact>"
write(tmp_unit,'(A)') "</ResourceHeader>"
write(tmp_unit,'(A)') "<Model>"
write(tmp_unit,'(A)') "<ModelID>impex://LATMOS/Hybrid</ModelID>"
write(tmp_unit,'(A)') "</Model>"
!write(tmp_unit,'(A)') "<SimulationType>Hybrid</SimulationType>"
write(tmp_unit,'(A)') "<TemporalDependence>No</TemporalDependence>"
!write(tmp_unit,'(A)') "<Language>FORTRAN2003</Language>"
write(string,'(3a)') "<ObservedRegion>",trim(adjustl(planetname)),"</ObservedRegion>"
write(tmp_unit,'(A)') trim(adjustl(string))

time_duration=nhm*dt*Spe%ref%inv_gyro
step_duration=dt*Spe%ref%inv_gyro
write(tmp_unit,'(A)') "<SimulationTime>"
write(string1,'(f10.3)') time_duration
write(string,'(3a)') "<Duration>PT",trim(adjustl(string1)),"S</Duration>"
write(tmp_unit,'(A)') trim(adjustl(string))
write(tmp_unit,'(A)') "<TimeStart>00:00:00</TimeStart>"
write(string1,'(f10.3)') step_duration
write(string,'(3a)') "<TimeStep>PT",trim(adjustl(string1)),"S</TimeStep>"
write(tmp_unit,'(A)') trim(adjustl(string))
write(tmp_unit,'(A)') "</SimulationTime>"

valid_min=-(Spe%P%centr)*Spe%ref%c_omegapi
valid_max=(nc_tot*gstep-Spe%P%centr)*Spe%ref%c_omegapi
box_size=nc_tot*gstep*Spe%ref%c_omegapi
grid_cell=gstep*Spe%ref%c_omegapi
write(tmp_unit,'(A)') "<SimulationDomain>"
write(tmp_unit,'(A)') "<CoordinateSystem>"
write(tmp_unit,'(A)') "<CoordinateRepresentation>Cartesian</CoordinateRepresentation>"
select case(trim(adjustl(planetname)))
  case("mars")
        write(tmp_unit,'(A)') "<CoordinateSystemName>MSO</CoordinateSystemName>"
  case default
        write(tmp_unit,'(A)') "<CoordinateSystemName></CoordinateSystemName>"
  end select!
write(tmp_unit,'(A)') "</CoordinateSystem>"
write(tmp_unit,'(A)') "<SpatialDimension>3</SpatialDimension>"
write(tmp_unit,'(A)') "<VelocityDimension>3</VelocityDimension>"
write(tmp_unit,'(A)') "<FieldDimension>3</FieldDimension>"
write(tmp_unit,'(A)') "<Units>km</Units>"
write(tmp_unit,'(A)') "<UnitsConversion> 1000 > m </UnitsConversion>"
write(tmp_unit,'(A)') "<AxesLabel>X Y Z</AxesLabel>"
write(string1,'(3(f10.1))') box_size(1),box_size(2),box_size(3)
write(string,'(3a)') "<BoxSize>",trim(adjustl(string1)),"</BoxSize>"
write(tmp_unit,'(A)') trim(adjustl(string))
write(string1,'(3(f10.1))') valid_min(1),valid_min(2),valid_min(3)
write(string,'(3a)') "<ValidMin>",trim(adjustl(string1)),"</ValidMin>"
write(tmp_unit,'(A)') trim(adjustl(string))
write(string1,'(3(f10.1))') valid_max(1),valid_max(2),valid_max(3)
write(string,'(3a)') "<ValidMax>",trim(adjustl(string1)),"</ValidMax>"
write(tmp_unit,'(A)') trim(adjustl(string))
write(tmp_unit,'(A)') "<GridStructure>Constant</GridStructure>"
write(string1,'(3(f8.1))') grid_cell(1),grid_cell(2),grid_cell(3)
write(string,'(3a)') "<GridCellSize>",trim(adjustl(string1)),"</GridCellSize>"
write(tmp_unit,'(A)') trim(adjustl(string))
write(tmp_unit,'(A)') "<Symmetry>Axial</Symmetry>"
write(tmp_unit,'(A)') "<BoundaryConditions>"
write(tmp_unit,'(A)') "<ParticleBoundary>"
write(tmp_unit,'(A)') "<FrontWall> absorbing </FrontWall>"
write(tmp_unit,'(A)') "<BackWall> absorbing </BackWall>"
write(tmp_unit,'(A)') "<SideWall> absorbing </SideWall>"
write(tmp_unit,'(A)') "<Obstacle> absorbing </Obstacle>"
write(tmp_unit,'(A)') "</ParticleBoundary>"
write(tmp_unit,'(A)') "<FieldBoundary>"
write(tmp_unit,'(A)') "<FrontWall> IMF </FrontWall>"
write(tmp_unit,'(A)') "<BackWall> Neuman zero-gradient </BackWall>"
write(tmp_unit,'(A)') "<SideWall> periodic </SideWall>"
write(tmp_unit,'(A)') "<Obstacle> absorbing </Obstacle>"
write(tmp_unit,'(A)') "</FieldBoundary>"
write(tmp_unit,'(A)') "</BoundaryConditions>"
write(tmp_unit,'(A)') "</SimulationDomain>"


write(tmp_unit,'(A)') "<RegionParameter>"
write(string,'(3a)') "<SimulatedRegion>",trim(adjustl(planetname)),"</SimulatedRegion>"
write(tmp_unit,'(A)') trim(adjustl(string))
write(string1,'(f8.2)') Spe%P%radius*sPE%ref%c_omegapi
write(string,'(3a)') '<Radius Units="km">',trim(adjustl(string1)),'</Radius>'
write(tmp_unit,'(A)') trim(adjustl(string))
write(string,'(3a)') "<SubLatitude>",trim(adjustl(string1)),"</SubLatitude>"
write(string1,'(f8.2)') Spe%P%sslat
write(tmp_unit,'(A)') trim(adjustl(string))
write(string1,'(f8.2)') Spe%P%ssl
write(string,'(3a)') "<SubLongitude>",trim(adjustl(string1)),"</SubLongitude>"
write(tmp_unit,'(A)') trim(adjustl(string))
write(tmp_unit,'(A)') "<Property>"
write(tmp_unit,'(A)') "<Name>Exobase Altitude</Name>"
write(tmp_unit,'(A)') "<Description>Altitude of the exobase</Description>"
write(tmp_unit,'(A)') "<PropertyQuantity>Positional</PropertyQuantity>"
write(tmp_unit,'(A)') "<Units>km</Units>"
write(string1,'(f8.2)') Spe%P%r_lim*Spe%ref%c_omegapi
write(string,'(3a)') "<PropertyValue>",trim(adjustl(string1)),"</PropertyValue>"
write(tmp_unit,'(A)') trim(adjustl(string))
write(tmp_unit,'(A)') "</Property>"
write(tmp_unit,'(A)') "<Property>"
write(tmp_unit,'(A)') "<Name>Ionopause Altitude</Name>"
write(tmp_unit,'(A)') "<Description>Altitude of the ionopause in the model"
write(tmp_unit,'(A)') "That is the altitude below which ion densities are computed analytically</Description>"
write(tmp_unit,'(A)') "<PropertyQuantity>Positional</PropertyQuantity>"
write(string1,'(f8.2)') Spe%P%r_iono*Spe%ref%c_omegapi
write(string,'(3a)') "<PropertyValue>",trim(adjustl(string1)),"</PropertyValue>"
write(tmp_unit,'(A)') trim(adjustl(string))
write(tmp_unit,'(A)') "</Property>"
write(tmp_unit,'(A)') "</RegionParameter>"

 rphi = phi*deg_to_rad
 rpsi = psi*deg_to_rad
 mag_val=(/cos(rphi)*sin(rpsi),sin(rphi)*sin(rpsi),cos(rpsi)/)
 mag_val=mag_val*Spe%ref%mag*1E9
!+++++++++ IMF ++++++++++++
write(tmp_unit,'(A)') "<InputField>"
write(tmp_unit,'(A)') "<Name>IMF</Name>"
write(tmp_unit,'(A)') "<Description>Interplanetary Magnetic Field</Description>"
write(tmp_unit,'(A)') "<SimulatedRegion>Heliosphere</SimulatedRegion>"
write(tmp_unit,'(A)') "<FieldQuantity>Magnetic</FieldQuantity>"
write(tmp_unit,'(A)') "<Units>nT</Units>"
write(tmp_unit,'(A)') "<InputLabel>Bx By Bz</InputLabel>"
write(string1,'(3(f8.2))') mag_val(1),mag_val(2),mag_val(3)
write(string,'(3a)') "<InputValue>",trim(adjustl(string1)),"</InputValue>"
write(tmp_unit,'(A)') trim(adjustl(string))
write(string1,'(f8.2)') Spe%ref%mag*1E9
write(string,'(3a)') "<ValidMin>",trim(adjustl(string1)),"</ValidMin>"
write(tmp_unit,'(A)') trim(adjustl(string))
write(string,'(3a)') "<ValidMax>",trim(adjustl(string1)),"</ValidMax>"
write(tmp_unit,'(A)') trim(adjustl(string))
write(tmp_unit,'(A)') "</InputField>"

!++++++++ SW population +++++
dens_sw0=Spe%ref%density/Spe%S(1)%percent
!1) e-
write(tmp_unit,'(A)') '<InputPopulation>'
write(tmp_unit,'(A)') '<Name>Solar Wind electrons</Name>'
write(tmp_unit,'(A)') '<Description>Solar wind electron fluid</Description>'
write(tmp_unit,'(A)') '<SimulatedRegion>Heliosphere</SimulatedRegion>'
write(tmp_unit,'(A)') '<SimulatedRegion>Incident</SimulatedRegion>'
write(tmp_unit,'(A)') '<PopulationType>electron</PopulationType>'
write(tmp_unit,'(A)') '<PopulationMassNumber>0</PopulationMassNumber>'
write(tmp_unit,'(A)') '<PopulationChargeState>-1</PopulationChargeState>' 
write(string1,'((es10.2E2))') dens_sw0*1E-6
write(string,'(3a)')  '<PopulationDensity Units="cm^-3">',trim(adjustl(string1)),'</PopulationDensity>'
write(tmp_unit,'(A)') trim(adjustl(string))
write(string1,'((f8.2))') Spe%betae*(Spe%ref%mag)**2/dens_sw0/(2.*4.*pi*1e-7*1.6E-19)
write(string,'(3a)')  '<PopulationTemperature Units="eV" UnitsConversion="11605 > K">',& 
&       trim(adjustl(string1)),'</PopulationTemperature>'
write(tmp_unit,'(A)') trim(adjustl(string))
write(string1,'((f8.2))') Spe%ref%alfvenspeed*Spe%S(1)%vxs
write(string,'(3a)') '<PopulationFlowSpeed Units="km/s">',trim(adjustl(string1)),'</PopulationFlowSpeed>'
write(tmp_unit,'(A)') trim(adjustl(string))
if (mod(ipe,2) == 0) then
write(tmp_unit,'(A)') '<Profile>Adiabatic</Profile>'
else
write(tmp_unit,'(A)') '<Profile>Isothermal</Profile>'
endif
write(tmp_unit,'(A)') '<Distribution>Maxwellian</Distribution>'
write(tmp_unit,'(A)') '</InputPopulation>'

spe1_mass=1._dp
spe1_charge=1._dp

do i=1,ns
write(tmp_unit,'(A)') '<InputPopulation>'
write(string,'(3a)') '<Name>Solar Wind ',trim(adjustl(Spe%S(i)%name)),'</Name>'
write(tmp_unit,'(A)') trim(adjustl(string))
write(string,'(3a)') '<Description>Solar wind ',trim(adjustl(Spe%S(i)%name)),'</Description>'
write(tmp_unit,'(A)') trim(adjustl(string))
write(tmp_unit,'(A)') '<ObservedRegion>Heliosphere</ObservedRegion>'
write(tmp_unit,'(A)') '<PopulationType>ion</PopulationType>'
write(string1,'((i3))') int(Spe%S(i)%rmass*spe1_mass)
write(string,'(3a)')  '<PopulationMassNumber>',trim(adjustl(string1)),'</PopulationMassNumber>'
write(tmp_unit,'(A)') trim(adjustl(string))
write(string1,'((i3))') int(Spe%S(i)%rcharge*spe1_charge)
write(string,'(3a)')  '<PopulationChargeState>',trim(adjustl(string1)),'</PopulationChargeState>'
write(tmp_unit,'(A)') trim(adjustl(string))
write(string1,'((es10.2E2))') dens_sw0*Spe%S(i)%percent*1E-6
write(string,'(3a)')  '<PopulationDensity Units="cm^-3">',trim(adjustl(string1)),'</PopulationDensity>'
write(tmp_unit,'(A)') trim(adjustl(string))
write(string1,'(f8.2)') Spe%S(i)%betas*(Spe%ref%mag)**2/dens_sw0/Spe%S(i)%percent/(2.*4.*pi*1e-7*1.6E-19)
write(string,'(3a)')  '<PopulationTemperature Units="eV" UnitsConversion="11605 > K">',&
&       trim(adjustl(string1)),'</PopulationTemperature>'
write(tmp_unit,'(A)') trim(adjustl(string))
write(string1,'(f8.2)') Spe%ref%alfvenspeed*Spe%S(i)%vxs
write(string,'(3a)') '<PopulationFlowSpeed Units="km/s">',trim(adjustl(string1)),'</PopulationFlowSpeed>'
write(tmp_unit,'(A)') trim(adjustl(string))
write(tmp_unit,'(A)') '<Distribution>Maxwellian</Distribution>'
write(tmp_unit,'(A)') '</InputPopulation>'
enddo
!--------SW------

!++++ Exo-Iono +++
write(tmp_unit,'(A)') '<InputPopulation>'
write(tmp_unit,'(A)') '<Name>Ionospheric electrons</Name>'
write(tmp_unit,'(A)') '<Description>Electron fluid in the planet ionosphere</Description>'
write(string,'(3a)') "<SimulatedRegion>",trim(adjustl(planetname)),"</SimulatedRegion>"
write(tmp_unit,'(A)') trim(adjustl(string))
write(tmp_unit,'(A)') '<PopulationType>electron</PopulationType>'
write(tmp_unit,'(A)') '<PopulationMassNumber>0</PopulationMassNumber>'
write(tmp_unit,'(A)') '<PopulationChargeState>-1</PopulationChargeState>' 
write(string1,'((f8.2))') Spe%betae*(Spe%ref%mag)**2/dens_sw0/(2.*4.*pi*1e-7*1.6E-19)*Spe%tempe_ratio
write(string,'(3a)')  '<PopulationTemperature Units="eV" UnitsConversion="11605 > K">',&
        & trim(adjustl(string1)),'</PopulationTemperature>'
write(tmp_unit,'(A)') trim(adjustl(string))
write(tmp_unit,'(A)') '<PopulationFlowSpeed Units="km/s">0</PopulationFlowSpeed>'
 select case(ipe)
  case(0)
        write(string1,'((f8.2))') Spe%ref%density*1E-6
        write(string,'(3a)')  '<Profile> Temperature for a density of ',trim(adjustl(string1)),'cm^-3,'
        write(tmp_unit,'(A)') trim(adjustl(string))
        write(tmp_unit,'(A)') 'Adiabatic profile</Profile>'
  case(1)
        write(tmp_unit,'(A)') '<Profile>Isothermal profile</Profile>'
  case(2)
        write(string1,'((f8.2))') Spe%ref%density*1E-5
        write(string,'(3a)')  '<Profile> Temperature for a density of ',trim(adjustl(string1)),'cm^-3,'
        write(tmp_unit,'(A)') trim(adjustl(string))
        write(string1,'((f8.2))') Spe%ref%density*3.98*1E-6
        write(string,'(3a)')  'Adiabatic profile below ',trim(adjustl(string1)),'cm^-3, hydrostatic above.</Profile>'
        write(tmp_unit,'(A)') trim(adjustl(string))
  case(3)
        write(string1,'((f8.2))') Spe%ref%density*1E-5
        write(string,'(3a)')  '<Profile> Temperature for a density of ',trim(adjustl(string1)),'cm^-3,'
        write(tmp_unit,'(A)') trim(adjustl(string))
        write(string,'(A)')  'Isothermal profile for lower densities, hydrostatic for higher.</Profile>'
  case default
    stop
  end select
 write(tmp_unit,'(A)') '</InputPopulation>'

do i=1,atmosphere%n_species
write(tmp_unit,'(A)') '<InputPopulation>'
if (atmosphere%species(i)%charge.eq.0) then 
write(string,'(3a)') '<Name>Exospheric ',trim(adjustl(atmosphere%species(i)%name)),'</Name>'
write(tmp_unit,'(A)') trim(adjustl(string))
write(string,'(3a)') '<Description>Exospheric ',trim(adjustl(atmosphere%species(i)%name)),'</Description>'
write(tmp_unit,'(A)') trim(adjustl(string))
write(string,'(3a)') "<SimulatedRegion>",trim(adjustl(planetname)),"</SimulatedRegion>"
write(tmp_unit,'(A)') trim(adjustl(string))
write(tmp_unit,'(A)') '<PopulationType>molecule</PopulationType>'
else
write(string,'(3a)') '<Name>Ionospheric ',trim(adjustl(atmosphere%species(i)%name)),'</Name>'
write(tmp_unit,'(A)') trim(adjustl(string))
write(string,'(3a)') '<Description>Ionospheric ',trim(adjustl(atmosphere%species(i)%name)),'</Description>'
write(tmp_unit,'(A)') trim(adjustl(string))
write(string,'(3a)') "<SimulatedRegion>",trim(adjustl(planetname)),"</SimulatedRegion>"
write(tmp_unit,'(A)') trim(adjustl(string))
write(tmp_unit,'(A)') '<PopulationType>ion</PopulationType>'
endif
write(string1,'((i3))') int(atmosphere%species(i)%mass*spe1_mass)
write(string,'(3a)')  '<PopulationMassNumber>',trim(adjustl(string1)),'</PopulationMassNumber>'
write(tmp_unit,'(A)') trim(adjustl(string))
write(string1,'((i3))') int(atmosphere%species(i)%charge*spe1_charge)
write(string,'(3a)')  '<PopulationChargeState >',trim(adjustl(string1)),'</PopulationChargeState>'
write(tmp_unit,'(A)') trim(adjustl(string))
write(tmp_unit,'(A)') '<PopulationTemperature Units="eV" UnitsConversion="11605 > K">0</PopulationTemperature>'
write(tmp_unit,'(A)') '<PopulationFlowSpeed Units="km/s">0</PopulationFlowSpeed>'
if (LEN_TRIM(adjustl(atmosphere%species(i)%description)).ne.0) &
        & write(tmp_unit,'(A)') trim(adjustl(atmosphere%species(i)%description))
write(tmp_unit,'(A)') '<Distribution>Dirac</Distribution>'
write(tmp_unit,'(A)') '</InputPopulation>'
enddo
!--------Exo-Iono------
!+++++ photoprod +++
do i=1,atmosphere%n_pp
write(tmp_unit,'(A)') '<InputProcess>'
write(string,'(5a)') '<Name>',&
        & trim(adjustl(atmosphere%photo_reactions(i)%mother%name)),&
        & ' + photon  = ',trim(adjustl(atmosphere%photo_reactions(i)%daughter%name)),&
        & ' + e-</Name>'
write(tmp_unit,'(A)') trim(adjustl(string))
write(tmp_unit,'(A)') '<ProcessType>PhotoIonization</ProcessType>'
write(tmp_unit,'(A)') '<Units>s-1</Units>'
write(string1,'((es10.2E2))')  atmosphere%photo_reactions(i)%frequency
write(string,'(3a)') '<ProcessCoefficient>',trim(adjustl(string1)),'</ProcessCoefficient>'
write(tmp_unit,'(A)') trim(adjustl(string))
write(tmp_unit,'(A)') '<ProcessCoeffType>Frequency</ProcessCoeffType>'
if (LEN_TRIM(adjustl(atmosphere%photo_reactions(i)%description)).ne.0) &
        & write(tmp_unit,'(A)') trim(adjustl(atmosphere%photo_reactions(i)%description))
write(tmp_unit,'(A)') '</InputProcess>'
enddo
!+++++ chex +++
do i=1,atmosphere%n_exc
write(tmp_unit,'(A)') '<InputProcess>'
write(string,'(5a)') '<Name>',trim(adjustl(atmosphere%exc_reactions(i)%ion%name)),&
        & '  =  ',trim(adjustl(atmosphere%exc_reactions(i)%neutral%name)),&
        & ' Charge Exchange</Name>'
write(tmp_unit,'(A)') trim(adjustl(string))
write(tmp_unit,'(A)') '<ProcessType>ChargeExchange</ProcessType>'
write(tmp_unit,'(A)') '<Units>cm-2</Units>'
write(string1,'((es10.2E2))')  atmosphere%exc_reactions(i)%cross_section/(dt*Spe%ref%Alfvenspeed*Spe%ref%inv_gyro*1.e9)
write(string,'(3a)') '<ProcessCoefficient>',trim(adjustl(string1)),'</ProcessCoefficient>'
write(tmp_unit,'(A)') trim(adjustl(string))
write(tmp_unit,'(A)') '<ProcessCoeffType>CrossSection</ProcessCoeffType>'
if (LEN_TRIM(adjustl(atmosphere%exc_reactions(i)%description)).ne.0) &
        & write(tmp_unit,'(A)') trim(adjustl(atmosphere%exc_reactions(i)%description))
write(tmp_unit,'(A)') '</InputProcess>'
enddo
!+++++ ei +++
do i=1,atmosphere%n_ei
write(tmp_unit,'(A)') '<InputProcess>'
write(string,'(5a)') '<Name>',trim(adjustl(atmosphere%ei_reactions(i)%neutral%name)),&
        & ' + e-  = ',trim(adjustl(atmosphere%ei_reactions(i)%ion%name)),&
        & ' + 2e-</Name>'
write(tmp_unit,'(A)') trim(adjustl(string))
write(tmp_unit,'(A)') '<ProcessType>ElectronImpact</ProcessType>'
if (LEN_TRIM(adjustl(atmosphere%ei_reactions(i)%description)).ne.0) &
        & write(tmp_unit,'(A)') trim(adjustl(atmosphere%ei_reactions(i)%description))
write(tmp_unit,'(A)') '</InputProcess>'
enddo


!++++++++++++++++++++AUTRE+++++++++++++++++++++++++++


!++++++++++++++++++++++++++++++++++++++++++++++++++++

write(tmp_unit,'(A)') '<InputParameter>'
write(tmp_unit,'(A)') '<Name>Derived Parameters</Name>'
write(tmp_unit,'(A)') '<ParameterQuantity> Other </ParameterQuantity>'
write(tmp_unit,'(A)') '<Property>'
write(tmp_unit,'(A)') '<Name>Solar UV Flux @ 10.7</Name>'
write(tmp_unit,'(A)') '<PropertyQuantity>SolarUVFlux</PropertyQuantity>'
!write(tmp_unit,'(A)') '<Units>km/s</Units>'
write(string1,'(f8.2)')  atmosphere%F107
write(string,'(3a)') '<PropertyValue>',trim(adjustl(string1)),'</PropertyValue>'
write(tmp_unit,'(A)') trim(adjustl(string))
write(tmp_unit,'(A)') '</Property>' 
write(tmp_unit,'(A)') '<Property>'
write(tmp_unit,'(A)') '<Name>Solar UV Flux @ 10.7 Average</Name>'
write(tmp_unit,'(A)') '<PropertyQuantity>ActivityIndex</PropertyQuantity>'
!write(tmp_unit,'(A)') '<Units>km/s</Units>'
write(string1,'(f8.2)')  atmosphere%F107_Avg
write(string,'(3a)') '<PropertyValue>',trim(adjustl(string1)),'</PropertyValue>'
write(tmp_unit,'(A)') trim(adjustl(string))
write(tmp_unit,'(A)') '</Property>' 
write(tmp_unit,'(A)') '<Property>' 
write(tmp_unit,'(A)') '<Name>Alfven Speed</Name>' 
write(tmp_unit,'(A)') '<PropertyQuantity>AlfvenVelocity</PropertyQuantity>'
write(tmp_unit,'(A)') '<Units>km/s</Units>'
write(string1,'(f8.2)')  Spe%ref%alfvenspeed
write(string,'(3a)') '<PropertyValue>',trim(adjustl(string1)),'</PropertyValue>'
write(tmp_unit,'(A)') trim(adjustl(string))
write(tmp_unit,'(A)') '</Property>' 
write(tmp_unit,'(A)') '<Property>'
write(tmp_unit,'(A)') '<Name>Alfven Mach Number</Name>'
write(tmp_unit,'(A)') '<PropertyQuantity>AlfvenMachNumber</PropertyQuantity>'
write(string1,'(f8.2)')  Spe%S(1)%vxs
write(string,'(3a)') '<PropertyValue>',trim(adjustl(string1)),'</PropertyValue>'
write(tmp_unit,'(A)') trim(adjustl(string))
write(tmp_unit,'(A)') '</Property>' 
write(tmp_unit,'(A)') '<Property>'
write(tmp_unit,'(A)') '<Name>Ion Inertial Length</Name>'
write(tmp_unit,'(A)') '<PropertyQuantity>Other</PropertyQuantity>'
write(tmp_unit,'(A)') '<Units>km</Units>'
write(string1,'(f8.2)')  Spe%ref%c_omegapi
write(string,'(3a)') '<PropertyValue>',trim(adjustl(string1)),'</PropertyValue>'
write(tmp_unit,'(A)') trim(adjustl(string))
write(tmp_unit,'(A)') '</Property>'
write(tmp_unit,'(A)') '</InputParameter>'
write(tmp_unit,'(A)') "</SimulationRun>"
write(tmp_unit,'(A)') "</Spase>"
 close(tmp_unit)

 end subroutine write_impex_xml

end module diag_impex_xml
