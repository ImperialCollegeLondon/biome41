c---------------------------------------------------------------------------
c    The BIOME4-system:	biome4main.f	1.0b1	22.10.99
c
c	Copyright (c) 1999 by Jed O. Kaplan
c     
c	See COPYING file for copying and redistribution conditions.
c
c	This program is free software; you can redistribute it and/or modify
c	it under the terms of the GNU General Public License as published by
c	the Free Software Foundation; version 2 of the License.
c
c	This program is distributed in the hope that it will be useful,
c	but WITHOUT ANY WARRANTY; without even the implied warranty of
c	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c	GNU General Public License for more details.
c
c	Contact info: jkaplan@bgc-jena.mpg.de
c---------------------------------------------------------------------------
c
c			B I O M E 4 M A I N . F
c
c---------------------------------------------------------------------------
c
c      This program is the main part of the BIOME4 system.  It calls the
c      setup and driver subroutines and contains a couple of utility 
c      subroutines used by this and other parts of the system.
c
c      Author:	Jed O. Kaplan
c      Date:	22 October 1999
c      Version:	v1.0b1
c      Revised:
c
c------------------------------------------------------------------------

      program biome4main
      
      use f90getopt

      implicit none
      include 'netcdf.inc'

      integer status

      integer limits(4)
      integer inputid,outputid
      integer vartypes(100),location(100),list(100)
      integer noutvars

      real globalparms(4)

c-------------------------------------

      call biome4setup(inputid,outputid,limits,
     >globalparms,noutvars,list,location,vartypes)

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
