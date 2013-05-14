#
# Makefile for libpfasst.
#

LIBPFASST ?= $(PWD)

include $(LIBPFASST)/Makefile.defs

#
# libpfasst
#

$(shell $(PY) mk/version.py src/pf_version.f90)

VPATHS = $(LIBPFASST)/src

FSRC = $(wildcard src/*.f90)
CSRC = $(wildcard ls src/*.c)
OBJ  = $(subst src,build,$(FSRC:.f90=.o) $(CSRC:.c=.o))

build/libpfasst.a: $(OBJ)
	$(AR) build/libpfasst.a $(OBJ)

include $(LIBPFASST)/Makefile.rules

rossby_adjustment: examples/rossby_adjustment/lib/libspatialdiscretization.a examples/rossby_adjustment/*.f90
		$(FC) $(FFLAGS) -I$(CURDIR)/examples/rossby_adjustment/inc -o $@ $^ -Lbuild -lpfasst -L$(CURDIR)/examples/rossby_adjustment/lib -lspatialdiscretization -lserialio $(LDFLAGS)

boussinesq_2d: examples/boussinesq_2d/lib/libspatialdiscretization.a examples/boussinesq_2d/*.f90
		$(FC) $(FFLAGS) -I$(CURDIR)/examples/boussinesq_2d/inc -o $@ $^ -Lbuild -lpfasst -L$(CURDIR)/examples/boussinesq_2d/lib -lspatialdiscretization -lserialio $(LDFLAGS)

.PHONY: clean

include Makefile.external

#
# dependencies
#

build/pf_utils.o:       build/pf_dtype.o
build/pf_timer.o:       build/pf_dtype.o
build/pf_cycle.o:       build/pf_dtype.o
build/pf_explicit.o:    build/pf_timer.o
build/pf_hooks.o:       build/pf_timer.o
build/pf_imex.o:        build/pf_timer.o
build/pf_implicit.o:    build/pf_timer.o
build/pf_mpi.o:         build/pf_timer.o
build/pf_pthreads.o:    build/pf_timer.o
build/pf_logger.o:      build/pf_hooks.o
build/pf_restrict.o:    build/pf_utils.o build/pf_timer.o build/pf_hooks.o
build/pf_interpolate.o: build/pf_restrict.o build/pf_hooks.o
build/pf_parallel.o:    build/pf_interpolate.o build/pf_hooks.o build/pf_cycle.o
build/sdc_quadrature.o: build/pf_dtype.o build/sdc_poly.o
build/pf_quadrature.o:  build/sdc_quadrature.o
build/pf_pfasst.o:      build/pf_utils.o build/pf_quadrature.o

build/pfasst.o:         build/pf_parallel.o build/pf_pfasst.o \
                        build/pf_implicit.o build/pf_explicit.o build/pf_imex.o \
                        build/pf_mpi.o build/pf_pthreads.o build/pf_cpthreads.o \
                        build/pf_version.o build/pf_logger.o build/pf_cycle.o
