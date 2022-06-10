# BIOME 4 model

This repository is being used to update the original Biome 4 code created by Jed 
Kaplan with the following aims:

* ensure it works with more recent releases of the netcdf library,
* maintain an updated Makefile for use with gfortran
* clear out some of the compilation warnings (mostly the trivial whitespace ones at the moment!)
* adding a more user friendly command line interface via [`f90getopt`](https://github.com/haniibrahim/f90getopt)

## New interface

The program now sets the input file, output file, global settings and CO2 
concentration from the command line

    Usage: biome4 -i input.nc -o output.nc \
                  -s settings.txt -c 410 
    Options:
      -h  --help      Print this help screen
    Required arguments:
      -i  --input     Input data file
      -o  --output    Output data file
      -s  --settings  Global settings file
      -c  --co2       CO2 (ppm)

## Original use

The `biome4` executable expected to find _two_ files in the directory from 
which it was run:

* `biome4outvars` set the netcdf naming and formatting for possible output
  variables.

* `biome4options` set some fairly static options, such as the output variable
  choices, but also set things that would change with _every_ run, including:
     
    * the CO2 concentration
    * the path to the input data (which then _had_ to be called `inputdata.nc`
    * the output directory

So if you had a directory with the following structure.

    /rds/general/home/dorme/biome4_2012/inputdata.nc
    /rds/general/home/dorme/biome4_2012/biome4outvars
    /rds/general/home/dorme/biome4_2012/biome4options
    /rds/general/home/dorme/biome4_2012/out/

The first three lines of biome4options could then be:

    /rds/general/home/dorme/biome4_2012/
    /rds/general/home/dorme/biome4_2012/out/
    400
 
And you could then:

    cd /rds/general/home/dorme/biome4_2012/
    /path/to/the/compile/biome4
    
All of the options at run time were then read from those two files.

This is really inflexible - for multiple runs, all the input files would have
to have the same file name and each run would need a separate option file.

## Other thoughts

* The latitude and longitude limits are either set as a window (which subsets the provided data) or
  fixed at -180:180 and -60:90. This doesn't work when the input file itself is a subset, so maybe
  read those dimensions from the input file.
