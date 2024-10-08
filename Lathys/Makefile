############################ -*- Mode: Makefile -*- ###########################
## Makefile --- Gestion Particules + Pickup avec une vitesse dirige suiavnt Z
## 
## Auteurs         : Ronan Modolo <Ronan.Modolo@cetp.ipsl.fr>
##                   Marco Mancini <man74cio@gmail.com>
##                   Sebastien Hess <sebastien.hess@latmos.ipsl.fr>
## Cr le         : 02 sep 2010 16:08

## Dern. mod. par  : RModolo <LATMOS/IPSL/UVSQ>
## Dern. mod. le   : 31 mai 2013 13:38
###############################################################################
#/opt/open-mpi/openmpi-1.4.2-ifort/bin/mpif90
#mpif90

COMPILER=mpiifort

IFORTFLAG=-O3 \
	-DHAVE_TIMING -DHAVE_NETCDF  \
	-stand f08 -DIntelFortran \
        #-DHAVE_NO_PLANET \
        #-DHAVE_FREQUENT_B_OUTPUT \
	#-pg -g -check all,noarg_temp_created  -traceback -WB  \
	# -opt-report 3 -opt-report-phase=pgo \
	# -parallel -par-report -vec-report1 \
	# -mpe=mpilog\
	# -openmp\
	# -DHAVE_WAVE_TEST\
	# -DHAVE_NETCDF\
	# -DHAVE_DOUBLE_PRECISION\
	# -DHAVE_MAXSIZEWRITE\
	# -DHAVE_NO_PLANET\
	# -heap-arrays 1024 -save -zero \

GFORTFLAG=-O3\
	-DHAVE_TIMING \
	-std=f2008 -w \
	 # -g -pg -Wall  -Wuninitialized -pedantic\
	 # -fbounds-check  -fbacktrace -DHAVE_DEBUG \
	# -Wconversion -Wno-tabs\
	# -fno-automatic \
	# -fcheck-array-temporaries  -Wconversion -Wno-tabs\
	# -Warray-temporaries\
	# -DHAVE_DOUBLE_PRECISION \
	# -DHAVE_MAXSIZEWRITE

PGFORTFLAG_WS= -fast -Mipa=fast -fpic -Bstatic_pgi -DNOTHAVE_MPI\
        -DHAVE_TIMING -DHAVE_NETCDF -DPCGIFortran -tp istanbul #-DHAVE_DEBUG\

PGFORTFLAG= -fast -Mipa=fast\
        -DHAVE_TIMING -DHAVE_NETCDF -DPGIFortran  #-DHAVE_DEBUG\


FFLAGS = $(IFORTFLAG)
FC = $(COMPILER)
FC_WS = mpif90
CDF_DIR = $(NETCDF_ROOT)
INC_CDF = $(CDF_DIR)/include/
LIB_CDF = $(CDF_DIR)/lib/

#Directories:
HERE := $(CURDIR)
SOURCES := $(HERE)/src
OBJ_DIR := $(HERE)/objects

#Dependencies
VPATH = $(SOURCES):$(OBJ_DIR)

OBJS_WS=defs_basis.o\
        defs_mpitype.o\
        defs_parametre.o\
        m_writeout.o\
        defs_basic_cdf.o\
        defs_arr3Dtype.o\
        defs_counts_types.o\
        defs_species.o\
        defs_particletype.o\
        defs_tregister.o\
        defs_diag_type.o\
        defs_variable.o\
        defs_grid.o\
        m_logo.o\
        diag_wrt_common_cdf.o\


OBJ_TEST = defs_basis.o\
	defs_parametre.o\
	defs_basic_cdf.o\
	getTimeInfo.o\

#OBJS_0= most basic modules
OBJS_0 = defs_basis.o\
	 defs_mpitype.o\
 	 defs_parametre.o\
	 m_writeout.o\

#OBJS_1= definitions bloc
OBJS_1 = \
	defs_counts_types.o\
	defs_diag_type.o\
	defs_basic_cdf.o\
	defs_species.o\
	defs_particletype.o\
	defs_tregister.o\
        defs_arr3Dtype.o\
        defs_variable.o\
	defs_grid.o\
	defs_atmospheretype.o\

#OBJS_2= miscellaneous bloc (except m_writeout (0) and m_restart (8))
OBJS_2 = \
	m_logo.o\
	m_timing.o\
	m_cmdline.o\
	m_rand_gen.o\
	particle_fluxes.o\
	field_add_waves.o\
	time_variation.o\

#Hack : particle_fluxes.o should be in particle bloc, but time_variation.o need it

#OBJS_3= atmosphere bloc
OBJS_3 =\
	atm_photoproduction.o\
	atm_ionosphere.o\
	atm_charge_exchange.o\
	atm_sections_efficaces.o\
	atm_magnetic_fields.o\
	atm_external_atmosphere.o\

#OBJS_4= environment bloc
OBJS_4 =\
	env_mars_FSU90.o\
	env_mars.o\
	env_mars_try.o\
	env_mercure.o\
	env_moon.o\
	env_ganymede.o\
	env_titan.o\
	env_earth.o\
	env_venus.o\
        environment.o\

#OBJS_5= field bloc
OBJS_5 = \
        field_cond_limit.o\
        field_lissage.o\
	field_pe.o\
	field_b.o\
	field_e.o\
        field.o\

#OBJS_6= particle bloc
OBJS_6 = \
	diag_wrt_common_cdf.o\
	diag_flux_part_imp.o\
	m_distribution_function.o\
	particle_sort.o\
        part_moment.o\
        particle_init.o\
	particle_creation.o\
	particle_com.o\
        particle.o\

#OBJS_7= diagnosis bloc
OBJS_7 = \
        diag_energy.o\
	diag_fields.o\
	diag_iono.o\
	diag_prod.o\
	diag_particles.o\
	diag_tm_results.o\
	diag_moment_species.o\
	diagnostique.o\

#OBJS_8= main bloc
OBJS_8 = \
	initialisation.o\
        m_restart.o\
        time_schedule.o\
        diag_impex_xml.o\
	hyb_3d.o

OBJS_diag = $(OBJS_0) $(OBJS_1) $(OBJS_2) $(OBJS_3) $(OBJS_4) $(OBJS_5) $(OBJS_6) $(OBJS_7) \
	80_mpi_enseble/m_verify_reassemble_mpi.o\
	80_mpi_enseble/m_split_cdf.o\
	80_mpi_enseble/m_merge_global_cdf.o\
	80_mpi_enseble/m_reasseble_unform.o\
	80_mpi_enseble/MPI_reassemble_file.o
	
OBJS_diag_IMPEX = $(OBJS_WS) \
        IMPEX_diag/m_VO.o\
        IMPEX_diag/ImpexTreeXML_generator.o\
	IMPEX_diag/m_IMPEX_cut.o\
	IMPEX_diag/m_IMPEX_moment_orbit.o\
	IMPEX_diag/m_IMPEX_spectro_orbit.o\
	IMPEX_diag/mpi_IMPEX_diag.o
	
OBJS_web_IMPEX = $(OBJS_WS) \
	IMPEX_WS/IMPEX_WS_InterpoleFields.o
	
OBJS_webfieldline_IMPEX = $(OBJS_WS) \
	IMPEX_FL/IMPEX_WS_getFieldLine.o	

OBJS_webSpectra_IMPEX = $(OBJS_WS) \
	IMPEX_Spectra/m_IMPEX_spectro_orbit.o\
	IMPEX_Spectra/IMPEX_getSpectra.o
	
OBJS_diag_sput = $(OBJS_0) $(OBJS_1) $(OBJS_2) $(OBJS_3) $(OBJS_4) $(OBJS_5) $(OBJS_6) $(OBJS_7) \
	sputtering_diag/photoprod_sputt.o\
	sputtering_diag/CX_sputt.o\
	sputtering_diag/eimpact_sputt.o\
	sputtering_diag/sputtering_input.o
	
OBJS_diag_ionflux = $(OBJS_0) $(OBJS_1) $(OBJS_2) $(OBJS_3) $(OBJS_4) $(OBJS_5) $(OBJS_6) $(OBJS_7) \
	Ion_flux/m_flux_computation.o\
	Ion_flux/Ion_flux_diag.o

OBJECTS = $(OBJS_0) $(OBJS_1) $(OBJS_2) $(OBJS_3) $(OBJS_4) $(OBJS_5) $(OBJS_6) $(OBJS_7) $(OBJS_8)

#%.o %.mod: %.F90 
#	cd $(OBJ_DIR) ;\
#	$(FC) $(FFLAGS) -I$(INC_CDF) -I$(SOURCES) -c  $< -o  $@ 

%.o %.mod: %.F90
	cd $(OBJ_DIR) ;\
	$(FC) $(FFLAGS) -I$(INC_CDF) -I$(SOURCES) -c  $< -o  $@ 



quiet_plasma: compo_0 compo_1 compo_2 compo_3 compo_4
	cd $(OBJ_DIR) ;\
	$(FC)  $(FFLAGS) -o $@ $(OBJECTS) -L$(LIB_CDF) -lnetcdf  -lnetcdff
	@mv $(OBJ_DIR)/quiet_plasma $(HERE)/. 2> /dev/null

compo_0: $(OBJS_0) $(OBJS_1) $(OBJS_2)

compo_1: $(OBJS_3) $(OBJS_4)

compo_2: $(OBJS_5) $(OBJS_6)

compo_3: $(OBJS_7)

compo_4: $(OBJS_8)

test: $(OBJ_TEST)
	cd $(OBJ_DIR) ;\
	$(FC_WS) $(PGFORTFLAG_WS) -o $@ $(OBJ_TEST) -L$(LIB_CDF) -lnetcdf -lnetcdff
	@mv $(OBJ_DIR)/test $(HERE)/. 2> /dev/null

diag: $(OBJS_diag)
	cd $(OBJ_DIR) ;\
	$(FC) $(FFLAGS) -o $@ $(OBJS_diag) -L$(LIB_CDF) -lnetcdf -lnetcdff
	@mv $(OBJ_DIR)/diag $(HERE)/. 2> /dev/null

diag_IMPEX: $(OBJS_diag_IMPEX)
	cd $(OBJ_DIR) ;\
	$(FC) $(PGFORTFLAG) -o $@ $(OBJS_diag_IMPEX) -L$(LIB_CDF) -lnetcdf -lnetcdff
	@mv $(OBJ_DIR)/diag_IMPEX $(HERE)/. 2> /dev/null	

diag_SPUTTERING: $(OBJS_diag_sput)
	cd $(OBJ_DIR) ;\
	$(FC) $(FFLAGS) -o $@ $(OBJS_diag_sput) -L$(LIB_CDF) -lnetcdf -lnetcdff
	@mv $(OBJ_DIR)/diag_SPUTTERING $(HERE)/. 2> /dev/null	

interpole_fields: $(OBJS_web_IMPEX)
	cd $(OBJ_DIR) ;\
	$(FC) $(FFLAGS) -o $@ $(OBJS_web_IMPEX) -L$(LIB_CDF) -lnetcdf -lnetcdff
	@mv $(OBJ_DIR)/interpole_fields $(HERE)/. 2> /dev/null	

getFieldLine: $(OBJS_webfieldline_IMPEX)
	cd $(OBJ_DIR) ;\
	$(FC_WS) $(PGFORTFLAG_WS) -o $@ $(OBJS_webfieldline_IMPEX) -L$(LIB_CDF) -lnetcdf -lnetcdff
	@mv $(OBJ_DIR)/getFieldLine $(HERE)/. 2> /dev/null

interpole_spectra: $(OBJS_webSpectra_IMPEX)
	cd $(OBJ_DIR) ;\
	$(FC_WS) $(PGFORTFLAG_WS) -o $@ $(OBJS_webSpectra_IMPEX) -L$(LIB_CDF) -lnetcdf -lnetcdff
	@mv $(OBJ_DIR)/interpole_spectra $(HERE)/. 2> /dev/null	
	
ion_flux: $(OBJS_diag_ionflux)
	cd $(OBJ_DIR) ;\
	$(FC) $(FFLAGS) -o $@ $(OBJS_diag_ionflux) -L$(LIB_CDF) -lnetcdf -lnetcdff
	@mv $(OBJ_DIR)/ion_flux $(HERE)/. 2> /dev/null	

all:    quiet_plasma\
	diag\
	diag_IMPEX\
	diag_SPUTTERING\
	interpole_fields\
	getFieldLine\
	ion_flux

clean:
	cd $(OBJ_DIR); rm -f  *.mod  *.o *.oo 
	cd $(OBJ_DIR)/80_*; rm -f  *.o *.oo
	cd $(OBJ_DIR)/IMPEX_*; rm -f  *.o *.oo	
	cd $(MOD_DIR); rm -f  *.mod  
	rm -f quiet_plasma diag diag_IMPEX diag_SPUTTERING interpole_fields getFieldLine interpole_spectra
	@bash ./tools/simple_configure 2>/dev/null 

cleanall: clean 
