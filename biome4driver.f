c---------------------------------------------------------------------------
c    The BIOME4-system:       biome4driver.f  1.0b1    22.10.99
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
c                      B I O M E 4 D R I V E R . F
c
c---------------------------------------------------------------------------
c
c      This subroutine is the driver to the BIOME4 model.  It fills the input
c      array and calls the biome4 subroutine which is the central part of the
c      model.  This subroutine also handles the output array and writes the
c      user-selected output variables to the netCDF output file.
c
c      Author:  Jed O. Kaplan
c      Date:    22 October 1999
c      Version: v1.0b1
c      Revised: 18.11.99 by JOK
c       Divided incoming sun value by 10. (line 175)because in 10 minute
c       dataset the input variable is x10.  This should be set as a flag.
c      Revised: 06.12.99 by KS
c       dimensions setting for small outputfiles (lon,lat)
c       varaibls: xo,yo,thislono  in lines 133-136
c      Revised: 14.12.99 by KS
c       change offset 29 to 28
c       varaibls: vars_out(k+28)  in line 188
c      Revised: 08.01.00 by JK
c       revised section on writing output variables to accommodate three
c       dimensional variable output (ie. montlhy values).
c       Also changed the type and dimension checking to query the output
c       dataset.  This may be more time-consuming but it is failsafe.
c------------------------------------------------------------------------

      subroutine biome4driver(inputid,outputid,limits,
     >           globalparms,noutvars,list,location,vartypes)

      implicit none
      include 'netcdf.inc'

c-------------------------------
c     variables

      logical diagmode

      integer*2 temp(12),prec(12),sun(12),tmin
      integer*2 shortval,ice,barren

      integer status
      integer k
      integer limits(4)
      integer inputid,outputid
      integer thispixel(2),thispixelo(2),thispixthr(3)
      integer getpixel(3),moncount(3),layercount(3)
      integer x,y,xo,yo
      integer maxx,maxy,minx,miny
      integer vartypes(100),location(100),list(100)
      integer var,noutvars
      integer xspan,yspan
      integer jobsize
      integer count,mark
      integer dot,dots
      integer ndims,i

      real globalparms(4)
      real p,co2
      real water
      real vars_out(50)
      real outputdata(500)
      real thislon,thislat
      real whc(2),perc(2)
      real realval
      real arrval(12)

      parameter (ice=28,barren=27)

      data moncount /1,1,12/
      data layercount /1,1,2/

1     format(I3,A1,$)
2     format(A1,$)

      p=globalparms(1)
      co2=globalparms(2)
      water=globalparms(3)

      if (globalparms(4).eq.1) then
       diagmode=.true.
      else
       diagmode=.false.
      end if

c-------------------------------

      minx=limits(1)
      maxx=limits(2)
      miny=limits(3)
      maxy=limits(4)

      xspan=maxx-minx
      yspan=miny-maxy

      jobsize=xspan*yspan

      write(*,*),'jobsize is',jobsize,' pixels'
      write(*,*)

      dot=nint(jobsize/100.)
      count=0
      mark=0
      dots=0

      if (.not.diagmode) write(*,1)mark,'%'

      do y=maxy,miny
       do x=minx,maxx

        status=nf_get_var1_real(inputid,1,x,thislon)
        if (status.ne.nf_noerr) call handle_err(status)

        status=nf_get_var1_real(inputid,2,y,thislat)
        if (status.ne.nf_noerr) call handle_err(status)

        thispixel(1)=x
        thispixel(2)=y

        getpixel(1)=x
        getpixel(2)=y
        getpixel(3)=1


        xo=x-minx+1
        yo=y-maxy+1
        thispixelo(1)=xo
        thispixelo(2)=yo

        thispixthr(1)=xo
        thispixthr(2)=yo
        thispixthr(3)=1

        status=nf_get_vara_int2(inputid,5,getpixel,moncount,temp)
        if (status.ne.nf_noerr) call handle_err(status)

        status=nf_get_vara_int2(inputid,6,getpixel,moncount,prec)
        if (status.ne.nf_noerr) call handle_err(status)

        status=nf_get_vara_int2(inputid,7,getpixel,moncount,sun)
        if (status.ne.nf_noerr) call handle_err(status)

        status=nf_get_vara_real(inputid,8,getpixel,layercount,whc)
        if (status.ne.nf_noerr) call handle_err(status)

        status=nf_get_vara_real(inputid,9,getpixel,layercount,perc)
        if (status.ne.nf_noerr) call handle_err(status)

        status=nf_get_var1_int2(inputid,10,thispixel,tmin)
        if (status.ne.nf_noerr) call handle_err(status)

c-------------------------------

        if (whc(1).ne.water)then

         if (whc(1).eq.-4.) then

          status=nf_put_var1_real(outputid,1,xo,thislon)
          if (status.ne.nf_noerr) call handle_err(status)

          status=nf_put_var1_real(outputid,2,yo,thislat)
          if (status.ne.nf_noerr) call handle_err(status)

          status=nf_put_var1_int2(outputid,5,thispixelo,ice)
          if (status.ne.nf_noerr) call handle_err(status)

         else if (whc(1).eq.-1.) then

          status=nf_put_var1_real(outputid,1,xo,thislon)
          if (status.ne.nf_noerr) call handle_err(status)

          status=nf_put_var1_real(outputid,2,yo,thislat)
          if (status.ne.nf_noerr) call handle_err(status)

          status=nf_put_var1_int2(outputid,5,thispixelo,barren)
          if (status.ne.nf_noerr) call handle_err(status)

         else

          vars_out(1)=thislat
          vars_out(2)=co2
          vars_out(3)=p
          vars_out(4)=real(tmin)/10.

          do k=1,12
           vars_out(k+4)=real(temp(k))/10.
           vars_out(k+16)=real(prec(k))
           vars_out(k+28)=real(sun(k))/10. !29  -->  28 in old file
          end do

          vars_out(41)=perc(1)       !perc top=k1
          vars_out(42)=perc(2)       !perc bottom=k2
          vars_out(43)=whc(1)        !whc top=k5
          vars_out(44)=whc(2)        !whc bottom=k6

          vars_out(46)=globalparms(4) !diagnostic mode indicator 1=on other=off
          vars_out(49)=thislon

c-------------------------------

          call biome4(vars_out,outputdata)

c-------------------------------

          status=nf_put_var1_real(outputid,1,xo,thislon)
          if (status.ne.nf_noerr) call handle_err(status)

          status=nf_put_var1_real(outputid,2,yo,thislat)
          if (status.ne.nf_noerr) call handle_err(status)

          do var=1,noutvars

           status=nf_inq_varndims(outputid,var+4,ndims)
           if (status.ne.nf_noerr) call handle_err(status)

           if (ndims.gt.2) then

           do i=1,12
            arrval(i)=outputdata(location(list(var))+(i-1))
           end do

            status=nf_put_vara_real
     >       (outputid,var+4,thispixthr,moncount,arrval)
            if (status.ne.nf_noerr) call handle_err(status)

           else

           if (vartypes(list(var)).eq.nf_short) then
            shortval=nint(outputdata(location(list(var))))
            status=nf_put_var1_int2(outputid,var+4,thispixelo,shortval)
            if (status.ne.nf_noerr) call handle_err(status)
           else
            realval=outputdata(location(list(var)))
            status=nf_put_var1_real(outputid,var+4,thispixelo,realval)
            if (status.ne.nf_noerr) call handle_err(status)
           end if

           end if

          end do

         end if

         else

         status=nf_put_var1_real(outputid,1,xo,thislon)
         if (status.ne.nf_noerr) call handle_err(status)

         status=nf_put_var1_real(outputid,2,yo,thislat)
         if (status.ne.nf_noerr) call handle_err(status)

        end if

c-------------------------------
c       here is the dot counter

        if (.not.diagmode) then
         count=count+1
         if (count.eq.dot) then
          dots=dots+1
          if (dots.eq.10) then
          mark=mark+10
           write(*,1)mark,'%'
           dots=0
          else
           write(*,2)'.'
          end if
          count=0
         end if
        end if

c-------------------------------

       end do
      end do

      return

      end
