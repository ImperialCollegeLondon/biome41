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
NETCDFINCDIR = /opt/netcdf/include
NETCDFLIB    = -lnetcdf
NETCDFLIBDIR = /opt/netcdf/lib

INCDIRS = -I$(NETCDFINCDIR)
LIBDIRS = -L$(NETCDFLIBDIR)
LIBS    = $(NETCDFLIB)

################################################################
## If you want to use another compiler instead of the
## the GNU g77 fortran compiler, change value for compile in the
## following line. 
################################################################
FC = g77

####################
## Can add a -g here
####################
#OTHERFLAGS = -g

################################################################
## You should not have to edit anything below this line        #
################################################################

MODELOBJS = biome4main.o biome4setup.o biome4driver.o biome4.o

FFLAGS = $(OTHERFLAGS) -fno-silent -Wall $(INCDIRS)

################################################################

all::	model

model:	$(MODELOBJS)
	$(FC) -o biome4 $(MODELOBJS) $(INCDIRS) $(LIBDIRS) $(LIBS)

clean::	
	-rm *.o
