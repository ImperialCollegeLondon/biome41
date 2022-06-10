##
## Makefile for BIOME4 series models
## Uses netCDF v3.x I/O
## Jed O. Kaplan, 19 October 1999
##

################################################################
## Edit these three to indicate the path for the netcdf include
## file 'netcdf.h', the name of the netcdf library file, and the
## path to that library file.
################################################################
NETCDFINCDIR = /apps/netcdf/4.4.4-fortran/include/
NETCDFLIB    = -lnetcdff
NETCDFLIBDIR = /apps/netcdf/4.4.4-fortran/lib/

INCDIRS = -I$(NETCDFINCDIR)
LIBDIRS = -L$(NETCDFLIBDIR)
LIBS    = $(NETCDFLIB)

################################################################
## If you want to use another compiler instead of the
## the GNU g77 fortran compiler, change value for compile in the
## following line.
################################################################
FC = gfortran

####################
## Can add a -g here
####################
#OTHERFLAGS = -g

################################################################
## You should not have to edit anything below this line        #
################################################################

MODELOBJS = biome4main.o biome4setup.o biome4driver.o biome4.o

FFLAGS = $(OTHERFLAGS) -Wall $(INCDIRS)

################################################################

all::	model

model:	$(MODELOBJS)
	$(FC) -c f90getopt.F90
	$(FC) f90getopt.F90 -o biome4 $(MODELOBJS) $(INCDIRS) $(LIBDIRS) $(LIBS)

clean::
	-rm *.o
