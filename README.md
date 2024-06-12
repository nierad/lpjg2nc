# lpjg2nc

A new conversion tool for classical LPJ-GUESS output (.out-files) to netCDF4. As opposed to out2nc it will create a lat-lon-grid which will be quick-viewable with standard viewers like ncview or panoply. It currently reads all classical LPJ-GUESS outputs (daily, yearly, and both, column and row monthly formats. 

<h2>Available Scripts</h2>

* lpjg2nc.py: This will simply convert an LPJ-GUESS out-file into a netCDF4-file, rechunked for quick viewing.
* nodejob_lpjg2nc.sh: A script to run lpjg2nc.py on a node. Take it as a template. 

<h2>Reguirements</h2>

* Python 3 
* netCDF4

<h2>Running it </h2>

Activate your environment
 ```
conda activate py36 (or "source activate py36" on some machines, e.g. Tetralith)
```
Run jobs: 
```
python lpjg2nc.py <in-file.out> <out-file.nc>
```
  Run job on compute-nodes
```
sbatch nodejob_lpjg2nc.sh <in-file.out> <out-file.nc>
```
When you're done
```
conda deactivate
```

<h2> Installation guide </h2>

Download a copy:
```
git clone https://github.com/nierad/lpjg2nc.git
cd lpjg2nc
```

Load Anaconda:
  e.g.
```
module load Anaconda (use "module spider Anaconda" to search for installations)
```
...and create an environment (py36 is the name here , but can be anything...)
  **Note:** Conda-environments will use the folder ```.conda``` in your HOME-dir. You might want to link it to some other place, as over time it can become quite large. Any conda-environment will be stored there.
```
conda create --name py36 python=3.6
```
Activate environment:
```
conda activate py36 (on some machines it is $ source activate py36)
```
Now install missing packages:
```
conda install netCDF4
```
Deactivate environment:
```
conda deactivate
```
That should do!

