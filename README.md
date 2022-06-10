# BIOME 4 model

This repository is being used to update the original Biome 4 code created by Jed 
Kaplan with the following aims:

* ensure it works with more recent releases of the netcdf library,
* maintain an updated Makefile for use with gfortran
* clear out some of the compilation warnings (mostly the trivial whitespace ones at the moment!)
* potentially add a more user friendly command line interface via [`f90getopt`](https://github.com/haniibrahim/f90getopt)

## Original use

The `biome4` executable expects to find _two_ files in the directory from which it is run:

* biome4outvars sets which variables are saved in the outputs

* biome4options sets a range of things including the CO2 concentration, something again about which variables are outputted but most critically the first two lines are the paths to the input data and output location. The input data _has_ to be at that input path and has to be called inputdata.nc. 

So if you had a directory with the following structure.

    /rds/general/home/dorme/biome4_2012/inputdata.nc
    /rds/general/home/dorme/biome4_2012/biome4outvars
    /rds/general/home/dorme/biome4_2012/biome4options
    /rds/general/home/dorme/biome4_2012/out/

The first two lines of biome4options should then be:

    /rds/general/home/dorme/biome4_2012/
    /rds/general/home/dorme/biome4_2012/out/

And you could then:

    cd /rds/general/home/dorme/biome4_2012/
    /path/to/the/compile/biome4
    
All of the options at run time are then read from those two files.

## Proposed UI

The vague plan is to try and get `biome4` to work with command line arguments, along the lines of:

    biome4 --input myinputfile.nc --output path/outfile.nc --options biome4options --outvars biome4outvars

It is conceivable that those last two options might be merged into a single config file?

## Other thoughts

* The latitude and longitude limits are either set as a window (which subsets the provided data) or
  fixed at -180:180 and -60:90. This doesn't work when the input file itself is a subset, so maybe
  read those dimensions from the input file.
