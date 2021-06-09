# lpjg2nc

This is a first test git repository for a new conversion tool for classical LPJ-GUESS output (.out-files) to netCDF4. 

Reguirements:

* Python 3 
  * netCDF4


Load conda environment:
 - ?> source activate env.name

How to use:
* Single file
  - python lpjg2nc.py <in-file.out> <out-file.nc>

* Run multiple jobs on compute-nodes
  - sbatch nodejob_lpjg2nc.sh <in-file.out> <out-file.nc>
