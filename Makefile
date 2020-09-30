#####################################################
###                                               ###
###      Makefile for the core_accr code          ###
###           James Cadman (2020)                 ###
###                                               ###
#####################################################

# Compiler variables:
FC     = gfortran

# To compile with OpenMP flags
# FFLAGS = -O3 -frecord-marker=4 -fdefault-real-8 -Wunused -fbounds-check -fopenmp

# To compile without OpenMP
FFLAGS = -O3 -frecord-marker=4 -fdefault-real-8 -Wunused -fbounds-check

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

SOURCESAF90 = calc_fg.f90 disc_module.f90 everything_mod.f90 setup.f90 \
				planet_form_new.f90

OBJECTSA    = $(SOURCESAF90:.f90=.o)

# Create executable files:
build: core_accr

core_accr:  $(OBJECTSA)
	$(FC) $(FFLAGS) -o $@ $(OBJECTSA)

# Clean statements:
clean:
	\rm *.o *.mod core_accr

# End Makefile
