c---------------------------------------------------------------------------
c    The BIOME4-system: biome4setup.f   1.0b1   22.10.99
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
c                     B I O M E 4 S E T U P . F
c
c---------------------------------------------------------------------------
c
c      This subroutine is the performs the setup to initialize a run of the
c      BIOME4 model.  The subroutine opens the netCDF input file and generates
c      a netCDF output file with the appropriate number of dimensions.
c      It also initializes the output variables selected by the user in the
c      biome4options file.  Biome4setup and the biome4driver subroutines expect
c      the input datafile to be of a certain standard; the list of dimensions
c      and variables is below.
c
c-------------------------------------
c     Standard for netCDF input data files (generated with ncdump):
c
c netcdf <<filename>> {
c dimensions:
c       lon = <<arbitrary size>> ;
c       lat = <<arbitrary size>> ;
c       time = 12 ;
c       soil_layer = 2 ;
c variables:
c       float lon(lon) ;
c              lon:long_name = "longitude" ;
c              lon:units = "degrees_east" ;
c              lon:missing_value = -9999.f ;
c       float lat(lat) ;
c              lat:long_name = "latitude" ;
c              lat:units = "degrees_north" ;
c              lat:missing_value = -9999.f ;
c       int time(time) ;
c              time:long_name = "time" ;
c              time:units = "month" ;
c              time:missing_value = -9999 ;
c       int soil_layer(soil_layer) ;
c              soil_layer:long_name = "soil layer, 0-30cm, 30cm-bottom" ;
c              soil_layer:units = "layer" ;
c              soil_layer:missing_value = -9999 ;
c       short temp(time, lat, lon) ;
c              temp:long_name = "monthly mean temperature" ;
c              temp:units = "degC" ;
c              temp:missing_value = -9999s ;
c       short prec(time, lat, lon) ;
c              prec:long_name = "monthly total precipitation" ;
c              prec:units = "mm" ;
c              prec:missing_value = -9999s ;
c       short sun(time, lat, lon) ;
c              sun:long_name = "mean monthly percent of possible sunshine" ;
c              sun:units = "percent" ;
c              sun:missing_value = -9999s ;
c       float whc(soil_layer, lat, lon) ;
c              whc:long_name = "soil water holding capacity" ;
c              whc:units = "mm/m" ;
c              whc:missing_value = -9999.f ;
c       float perc(soil_layer, lat, lon) ;
c              perc:long_name = "soil water percolation index" ;
c              perc:units = "mm/hr" ;
c              perc:missing_value = -9999.f ;
c       short tmin(lat, lon) ;
c              tmin:long_name = "annual absolute mimimum temperature" ;
c              tmin:units = "degC" ;
c              tmin:missing_value = -9999s ;
c
c // global attributes:
c              :title = "<<dataset title>>" ;
c }
c
c-------------------------------------
c
c      Author:       Jed O. Kaplan
c      Date:       22 October 1999
c      Version:       v1.0b1
c      Revised: 18.11.99 by JOK
c       Moved opening input data file and the box calculation to the
c       beginning of the program so that the output file is dynamically
c       sized to the run size.
c       Also changed the dimensions in the outfile to the calculated size.
c      Revised: 10.01.00 by JOK
c       Added two new global attributes to the output file: CO2 concentration
c       and resolution on the grid.  The grid resolution is set automatically
c       at the moment but should be automatic.
c      Revised: 12.01.00 by JOK
c       Added path and filename of the driving data input file.
c       Corrected output file creation so that output is placed in directory
c       specified in the options file and not simply in the current directory.
c
c------------------------------------------------------------------------
c     Program code begins here:

      subroutine biome4setup(infile,outfile,inputid,outputid,limits,
     >globalparms,noutvars,list,location,vartypes)

      implicit none
      include 'netcdf.inc'

c     variables

      character*15 outname
      character*20 dimname,varnames(100),var_units(100),var_unit
      character*20 varname,dimunits
      character*30 timestr
      character*60 var_longnames(100)
      character*60 dim_longname,var_longname
      character*65 inputpath,outputpath
      character*80 optionsfile,attributefile
      character*80 output_title,input_title
      character*120 header,infile,outfile

      integer lonsize,latsize
      integer inputid,outputid,status
      integer iopt
      integer length,pos,noutvars
      integer list(100)
      integer i
      integer dimsize,dimtype
      integer vardims(100),vardim
      integer dimids(100,4),dimid(4)
      integer vartypes(100),vartype
      integer location(100)
      integer whcid
      integer typenum(6),num
      integer varid,ice
      parameter (ice=28)

      integer limits(4)

      real lonlatbox(4)
      real var_missing
      real globalparms(4)
      real p,co2
      real gridres
      real water

      parameter (p=1E5)

c-------------------------------------

      data typenum
     > /nf_byte,nf_char,nf_short,nf_int,nf_float,nf_double/

c-------------------------------------
c     Read in the user run options

      optionsfile='biome4options'
      open(99,file=optionsfile,status='old')

      read(99,*)inputpath
      read(99,*)outputpath

      read(99,*) co2

      read(99,*) globalparms(4) !diagnostic mode option

      i=0
15    i=i+1
      read(99,*)list(i)
      if (list(i).ne.99) then
       goto 15
      else
       continue
      end if

      noutvars=i-1

      read(99,*) iopt

      if(iopt.eq.2) then
       lonlatbox(1)=-180.0
       lonlatbox(2)=180.0
       lonlatbox(3)=-60.0
       lonlatbox(4)=90.0
      else if(iopt.eq.5)then
       read(99,*) lonlatbox(1)
       read(99,*) lonlatbox(2)
       read(99,*) lonlatbox(3)
       read(99,*) lonlatbox(4)
      end if

c     Read the output attributes file
      print*,'reading attributes file'
      attributefile='biome4outvars'
      open(97,file=attributefile,status='old')

      do i=1,17
       read(97,*)header
      end do

      i=0
5     i=i+1
       read(97,*,end=10)varnames(i),vardims(i),
     >           dimids(i,1),dimids(i,2),dimids(i,3),dimids(i,4),
     >           num,var_longnames(i),var_units(i),location(i)
       vartypes(i)=typenum(num)
      go to 5

10    continue

c-------------------------------------
c     Open the netCDF input file

      print*,'opening input file', infile

      status=nf_open(infile, nf_nowrite,inputid)
      if (status.ne.nf_noerr) call handle_err(status)

      input_title=infile

c-------------------
c      get the boundaries of the box
c      find the x and y values for the min and max lat and lon

      status=nf_inq_dimlen(inputid,1,lonsize)
      if (status.ne.nf_noerr) call handle_err(status)

      status=nf_inq_dimlen(inputid,2,latsize)
      if (status.ne.nf_noerr) call handle_err(status)

      call box(inputid,lonsize,latsize,lonlatbox,limits)

c      print*,lonlatbox
c      print*,limits
c      stop

c-------------------------------------
c     Create the output file

      print*,'creating output file',outfile

      status=nf_create(outfile,nf_clobber,outputid)
      if (status.ne.nf_noerr) call handle_err(status)

      call timestring(timestr)

      output_title=
     >'BIOME4 output file, generated '//timestr(1:length(timestr))

      status=nf_put_att_text
     > (outputid,nf_global,'title',length(output_title),
     >  output_title(1:length(output_title)))
      if (status.ne.nf_noerr) call handle_err(status)

      status=nf_put_att_text
     > (outputid,nf_global,'input_file',length(input_title),
     >  input_title(1:length(input_title)))
      if (status.ne.nf_noerr) call handle_err(status)


      status=nf_put_att_real
     >(outputid,nf_global,'pCO2',nf_float,1,co2)
      if (status.ne.nf_noerr) call handle_err(status)

      gridres=30.

      status=nf_put_att_real
     >(outputid,nf_global,'resolution',nf_float,1,gridres)
      if (status.ne.nf_noerr) call handle_err(status)

c-----------------
c     define lon, lat and time dimensions and variables

      dimname='lon'
      dimsize=1+limits(2)-limits(1)
      dim_longname='longitude'
      dimtype=nf_float
      dimunits='degrees_east'

      call definedim
     >(outputid,dimname,dimsize,dim_longname,dimtype,dimunits)

      dimname='lat'
      dimsize=1+limits(3)-limits(4)
      dim_longname='latitude'
      dimtype=nf_float
      dimunits='degrees_north'

      call definedim
     >(outputid,dimname,dimsize,dim_longname,dimtype,dimunits)

      dimname='time'
      dimsize=12
      dim_longname='time'
      dimtype=nf_int
      dimunits='month'

      call definedim
     >(outputid,dimname,dimsize,dim_longname,dimtype,dimunits)

      dimname='layer'
      dimsize=2
      dim_longname='soil layer'
      dimtype=nf_int
      dimunits='layer'

      call definedim
     >(outputid,dimname,dimsize,dim_longname,dimtype,dimunits)

c-----------------
c     define the variables to be output
c     a library of common outputs is stored in the file biome4outvars

      do pos=1,noutvars

       varname=varnames(list(pos))
       vardim=vardims(list(pos))

       do i=1,4
        dimid(i)=dimids(list(pos),i)
       end do

       vartype=vartypes(list(pos))
       var_longname=var_longnames(list(pos))
       var_unit=var_units(list(pos))

       var_missing=-9999.

       call makevar(outputid,varname,vartype,vardim,dimid,
     >              var_longname,var_unit,var_missing)

      end do

c-------------------------------------
c     The output file is defined.  Take it out of define mode so it
c     is ready for writing.

      status=nf_enddef(outputid)
      if (status.ne.nf_noerr) call handle_err(status)

c     fill the time and layer dimension variables with data

      varid=3
      do i=1,12
       status=nf_put_var1_int(outputid,varid,i,i)
       if (status.ne.nf_noerr) call handle_err(status)
      end do

      varid=4
      do i=1,2
       status=nf_put_var1_int(outputid,varid,i,i)
       if (status.ne.nf_noerr) call handle_err(status)
      end do

c-------------------------------------
c     Find the flag for water

      status=nf_inq_varid(inputid,'whc',whcid)
      if (status.ne.nf_noerr) call handle_err(status)

      status=nf_get_att_real(inputid,whcid,'missing_value',water)
      if (status.ne.nf_noerr) call handle_err(status)

      globalparms(1)=p
      globalparms(2)=co2
      globalparms(3)=water

      print*,'setup complete'

c-------------------
c      the program is set up. return

       return

       end

c---------------------------------------------------------------------
      subroutine definedim
     >(outputid,dimname,dimsize,dim_longname,dimtype,dimunits)

      implicit none
      include 'netcdf.inc'

      integer status

      character*20 dimunits,dimname
      character*60 dim_longname

      integer outputid,dimsize,dimid
      integer dimtype,dimids(4)

      status=nf_def_dim(outputid,dimname,dimsize,dimid)
      if (status.ne.nf_noerr) call handle_err(status)

      dimids(1)=dimid

      call makevar(outputid,dimname,dimtype,1,dimids,
     >             dim_longname,dimunits,-9999.)

      return
      end

c-------------------------------------
       subroutine makevar(outputid,varname,vartype,vardim,dimid,
     >                    var_longname,var_unit,var_missing)

      implicit none
      include 'netcdf.inc'

      integer outputid,status,varid
      integer vardim
      integer vartype,i
      integer dimid(4),ldimid(vardim)
      integer length
      integer imissing
      integer*2 i2missing

      real var_missing

      character*60 var_longname
      character*20 varname,var_unit

      imissing=nint(var_missing)
      i2missing=nint(var_missing)

      do i=1,vardim
       ldimid(i)=dimid(i)
      end do

      status=nf_def_var(outputid,varname,vartype,vardim,ldimid,varid)
      if (status.ne.nf_noerr) call handle_err(status)

      status=nf_put_att_text
     >  (outputid,varid,'long_name',length(var_longname),var_longname)
      if (status.ne.nf_noerr) call handle_err(status)

      status=nf_put_att_text
     >  (outputid,varid,'units',length(var_unit),var_unit)
      if (status.ne.nf_noerr) call handle_err(status)

      if (vartype.eq.nf_byte) then
       status=nf_put_att_int1
     > (outputid,varid,'missing_value',vartype,1,var_missing)
       if (status.ne.nf_noerr) call handle_err(status)
       status=nf_put_att_int1
     > (outputid,varid,'_FillValue',vartype,1,var_missing)
       if (status.ne.nf_noerr) call handle_err(status)

      else if (vartype.eq.nf_char) then
c       status=nf_put_att_text
c     >  (outputid,varid,'missing_value',vartype,length(missing),missing)
c       if (status.ne.nf_noerr) call handle_err(status)

      else if (vartype.eq.nf_short) then
       status=nf_put_att_int2
     >  (outputid,varid,'missing_value',vartype,1,i2missing)
       if (status.ne.nf_noerr) call handle_err(status)
       status=nf_put_att_int2
     >  (outputid,varid,'_FillValue',vartype,1,i2missing)
       if (status.ne.nf_noerr) call handle_err(status)

      else if (vartype.eq.nf_int) then
       status=nf_put_att_int
     >  (outputid,varid,'missing_value',vartype,1,imissing)
       if (status.ne.nf_noerr) call handle_err(status)
       status=nf_put_att_int
     >  (outputid,varid,'_FillValue',vartype,1,imissing)
       if (status.ne.nf_noerr) call handle_err(status)

      else if (vartype.eq.nf_float) then
       status=nf_put_att_real
     >  (outputid,varid,'missing_value',vartype,1,var_missing)
       if (status.ne.nf_noerr) call handle_err(status)
       status=nf_put_att_real
     >  (outputid,varid,'_FillValue',vartype,1,var_missing)
       if (status.ne.nf_noerr) call handle_err(status)

      else if (vartype.eq.nf_double) then
       status=nf_put_att_double
     >  (outputid,varid,'missing_value',vartype,1,var_missing)
       if (status.ne.nf_noerr) call handle_err(status)
       status=nf_put_att_double
     >  (outputid,varid,'_FillValue',vartype,1,var_missing)
       if (status.ne.nf_noerr) call handle_err(status)

      end if

      return
      end

c-------------------------------------

      subroutine box(inputid,lonsize,latsize,lonlatbox,limits)
c     finds the x,y index limits of the selected lat,lon box

      implicit none
      include 'netcdf.inc'

      integer inputid
      integer lonsize,latsize
      integer limits(4)
      integer minx,maxx,miny,maxy
      integer status,i

      real lonlatbox(4),milo,malo,mila,mala
      real lon(lonsize),lat(latsize)

c----------------------

      milo=lonlatbox(1)
      malo=lonlatbox(2)
      mila=lonlatbox(3)
      mala=lonlatbox(4)

      status=nf_get_var_real(inputid,1,lon)
      if (status.ne.nf_noerr) call handle_err(status)

      status=nf_get_var_real(inputid,2,lat)
      if (status.ne.nf_noerr) call handle_err(status)

      do i=1,lonsize
       if (lon(i).le.milo) minx=i
       if (lon(i).le.malo) maxx=i
      end do

      maxy=1
      do i=1,latsize
       if (lat(i).ge.mila) miny=i
       if (lat(i).ge.mala) maxy=i
      end do

      limits(1)=minx
      limits(2)=maxx
      limits(3)=miny
      limits(4)=maxy

      return
      end


c-------------------------------------
