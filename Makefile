FC = gfortran
FCFLAGS = -g -c -fdefault-real-8 -fPIC -fno-second-underscore -fbacktrace -fno-align-commons -fbounds-check -std=legacy
LDFLAGS =

MODDIR := .mod
ifneq ($(MODDIR),)
  $(shell test -d $(MODDIR) || mkdir -p $(MODDIR))
  FCFLAGS+= -J $(MODDIR)
endif

SRCS =  MOM_error_handler.F90\
	MOM_document.F90\
        pkg/MOM6/src/framework/MOM_string_functions.F90\
        MOM_file_parser.F90\
        pkg/MOM6/src/framework/MOM_unit_scaling.F90\
        MOM_io.F90 pkg/MOM6/src/core/MOM_verticalGrid.F90 pkg/MOM6/src/framework/MOM_unit_scaling.F90\
	MOM_get_input.F90\
	pkg/pyRemapping/remap.f90\
        pkg/MOM6/src/ALE/polynomial_functions.F90\
	pkg/MOM6/src/ALE/regrid_solvers.F90 pkg/MOM6/src/ALE/regrid_edge_values.F90\
	pkg/MOM6/src/ALE/regrid_consts.F90 pkg/MOM6/src/ALE/P1M_functions.F90\
	pkg/MOM6/src/ALE/regrid_consts.F90 pkg/MOM6/src/ALE/PCM_functions.F90\
	pkg/MOM6/src/ALE/regrid_consts.F90 pkg/MOM6/src/ALE/P3M_functions.F90\
	pkg/MOM6/src/ALE/regrid_consts.F90 pkg/MOM6/src/ALE/PLM_functions.F90\
	pkg/MOM6/src/ALE/regrid_consts.F90 pkg/MOM6/src/ALE/PPM_functions.F90\
	pkg/MOM6/src/ALE/regrid_consts.F90 pkg/MOM6/src/ALE/PQM_functions.F90\
	pkg/MOM6/src/ALE/regrid_interp.F90 pkg/MOM6/src/ALE/MOM_remapping.F90\
	pkg/MOM6/src/ALE/coord_zlike.F90 pkg/MOM6/src/ALE/coord_sigma.F90\
	MOM_EOS.F90 MOM_grid.F90\
	MOM_regridding.F90

OBJECTS = $(SRCS:.F90=.o)
TARGET = libRegrid.a


$(TARGET): $(OBJECTS)
	rm -f $@
	ar cr $@ $^
	python setup.py config_fc --f90flags="-g -c -fdefault-real-8 -fPIC -fno-second-underscore -fbacktrace -fno-align-commons -fbounds-check" --fcompiler=gfortran build
	python setup.py install

%.o: %.F90
	$(FC) $(FCFLAGS) -I pkg/MOM6/config_src/memory/dynamic_symmetric -I pkg/MOM6/src/framework -c $< -o $@


test: $(TARGET)
	python test.py

clean:
	 rm -rf build *.o *.mod *.MOD .mod libRegrid.a pkg/MOM6/src/framework/*.o pkg/MOM6/src/core/*.o pkg/MOM6/src/ALE/*.o pkg/pyRemapping/*.o
