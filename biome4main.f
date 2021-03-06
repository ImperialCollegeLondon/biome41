c---------------------------------------------------------------------------
c    The BIOME4-system:       biome4main.f   1.0b1   22.10.99
c
c       Copyright (c) 1999 by Jed O. Kaplan
c
c       See COPYING file for copying and redistribution conditions.
c
c       This program is free software; you can redistribute it and/or modify
c       it under the terms of the GNU General Public License as published by
c       the Free Software Foundation; version 2 of the License.
c
c       This program is distributed in the hope that it will be useful,
c       but WITHOUT ANY WARRANTY; without even the implied warranty of
c       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c       GNU General Public License for more details.
c
c       Contact info: jkaplan@bgc-jena.mpg.de
c---------------------------------------------------------------------------
c
c                       B I O M E 4 M A I N . F
c
c---------------------------------------------------------------------------
c
c      This program is the main part of the BIOME4 system.  It calls the
c      setup and driver subroutines and contains a couple of utility
c      subroutines used by this and other parts of the system.
c
c      Author:       Jed O. Kaplan
c      Date:       22 October 1999
c      Version:       v1.0b1
c      Revised:
c
c------------------------------------------------------------------------

      program biome4main

      use f90getopt

      implicit none
      include 'netcdf.inc'

      ! Local declarations for command line arguments
      character(120) infile,outfile,settings
      integer stderr
      integer status
      real co2
      
      logical infile_not_set,outfile_not_set
      logical settings_not_set,co2_not_set
      integer limits(4)
      integer inputid,outputid
      integer vartypes(100),location(100),list(100)
      integer noutvars

      real globalparms(4)

      ! Declaration of longopts only (optional)
      ! ----------------------------------
      ! option_s derived type:
      !   1st value = long option name (character array, max. 80)
      !   2nd value = if option has value (boolean)
      !   3rd value = short option name (single character), same as in getopt()
      ! option_s is not needed if you just use short options
      type(option_s) :: opts(5)
      opts(1) = option_s("help",    .false.,  "h")
      opts(2) = option_s("infile",  .true.,  "i")
      opts(3) = option_s("outfile", .true.,  "o")
      opts(4) = option_s("co2",     .true.,  "c")
      opts(5) = option_s("settings",.true.,  "s")

      stderr = 0
      infile_not_set = .true.
      outfile_not_set = .true.
      settings_not_set = .true.
      co2_not_set = .true.

      ! If no options were committed
      ! ----------------------------
      if (command_argument_count() .eq. 0) then
          write(stderr,*) "ERROR:  Use options -h or --help for details"
          stop
      end if

      ! Processing options
      ! ------------------------
      ! Process short options one by one and long options if specified in option_s
      !
      ! getopt(optstr, longopt):
      !  - optstr = character of short option character without a space
      !             ":" after a character says that this option requires a value
      !  - opts   = longopts, if specified in option_s (optional)
      do
            select case(getopt("i:o:h", opts))
            case(char(0))
                  exit
            case("i") ! option --infile
                  infile = trim(optarg)
                  infile_not_set = .false.
            case("o") ! option --outfile
                  outfile = trim(optarg)
                  outfile_not_set = .false.
            case("s") ! option --settings
                  settings = trim(optarg)
                  settings_not_set = .false.
            case("c") ! option --co2
                  if (isnum(trim(optarg)) > 0) then ! Check for number in "optarg"
                      read(optarg,*) co2 ! Convert character string to double precision
                      co2_not_set = .false.
                  else
                      write(stderr,*) "ERROR: -c/--co2 is not a number."
                      stop
                  end if
            case("h") ! help output
                  print*, "Usage: biome4 -i input.nc -o output.nc \"
                  print*, "              -s settings.txt -c 410 "
                  print*, "Options:"
                  print*, "  -h  --help      Print this help screen"
                  print*, "Required arguments:"
                  print*, "  -i  --input     Input data file"
                  print*, "  -o  --output    Output data file"
                  print*, "  -s  --settings  Global settings file"
                  print*, "  -c  --co2       CO2 (ppm)"
                  stop
            end select
      end do
      
      ! Check for missing settings.
      if (infile_not_set .or. outfile_not_set .or. 
     >    settings_not_set .or. co2_not_set) then
           write(stderr,*) "ERROR: missing required settings"
           stop
      end if

c-------------------------------------

      call biome4setup(infile,outfile,settings,co2,
     >inputid,outputid,limits,globalparms,noutvars,
     >list,location,vartypes)

      call biome4driver(inputid,outputid,limits,
     >globalparms,noutvars,list,location,vartypes)

c-------------------------------------
c     Close files and clean up

5     status=nf_close(inputid)
      if (status.ne.nf_noerr) call handle_err(status)

      status=nf_close(outputid)
      if (status.ne.nf_noerr) call handle_err(status)

      end

c---------------------------------------------------------------------------
      subroutine handle_err(status)

c     NOTE always add this subroutine to programs
c     that use the NetCDF libraries. It handles the
c     errors and you will be much happier you did.

      implicit none
      include 'netcdf.inc'

      integer status

      if (status.ne.nf_noerr) then
       print *,nf_strerror(status)
       stop
      end if

      return
      end

c-------------------------------------
      integer function length(string)

c     Returns the length of a character string excluding trailing blanks

      implicit none

      character*(*) string

      length=len(string)

      do while (string(length:length).eq.' ')
       length=length-1
      end do

      end

c-------------------------------------
      subroutine timestring(timestr)

c     Returns a character string with the current time and date

      implicit none

      character*30 timestr

      timestr=ctime(time())

      return
      end

c-------------------------------------
