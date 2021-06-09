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


Installation guide:

  Download a copy:
    $ git clone https://github.com/nierad/lpjg2nc.git

  Load Anaconda and create an environment:
    e.g. $ module load Anaconda (use "$ module spider Anaconda" to search for installations)
    conda create --name py36 python=3.6

  Activate environment:
    $ conda activate py36 (on some machines it is $ source activate py36)

  Now install missing packages:
    $ conda install netCDF4

  That should do!

Running it:

  Activate your environment
    $ conda activate py36

  Run jobs 
    $ python lpjg2nc.py <in-file.out> <out-file.nc>

  if you're done
    $ conda deactivate

