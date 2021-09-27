## module load Anaconda3 first
# use environment py36: conda  activate py36
# might have to ?> conda install netCDF4
# Then rechunk
# nccopy -c time/20,latitude/36,longitude/72 in.nc out1.nc
# then 
# nccopy -c time/1,latitude/180,longitude/360 out1.nc out.nc
#
# Stats:
# Monthly file, 6 Variables:
# original out-file: 6.1 GB
# Time to process nc-file: ~9min -> 1.8 GB
# first nccopy 5:20 min -> 8.6 GB
# Second nccopy: 1:10 s -> 8.6 GB

import sys
import os
import numpy as np
import netCDF4 as nc

ifile=sys.argv[1]
ofile=sys.argv[2]
slabsize=1
if len(sys.argv)>3:
    slabsize = sys.argv[3]

#Check file size and time period
lcnt=0
gcnt=1
with open(ifile,"r") as fi:
    header=fi.readline()
    lo,la,y=fi.readline().split()[:3]
    sy=np.int(y)

nyears = 1
ey=sy
print(header)
# Check time-range and output-frequency
idx_u = os.path.basename(os.path.splitext(ifile)[0]).index("_")
freq  = os.path.basename(os.path.splitext(ifile)[0])[idx_u+1:]

if freq == "monthly":
    if header.split()[3]=="Jan" and header.split()[14]=="Dec":
        freq     = "monthly_col"
        tdimsize = nyears*12
        tunit    = "months since 01 "+str(sy)
        barup    = 100
    else:
        freq = "monthly_row"
        tdimsize = nyears*12
        tunit    = "months since 01 "+str(sy)
        #    ncells   = (lcnt-1)/(12*nyears)
        barup    = 100
elif freq == "yearly":
        tdimsize = nyears
        tunit    = "yearss since "+str(sy)
#        ncells   = (lcnt-1)/nyears
        barup    = 1000
elif freq == "daily":
    tdimsize = nyears*365
    tunit    = "days since 01 January "+str(sy)
#    ncells   = (lcnt-1)/(365*nyears)
    barup    = 1
else:
    sys.exit(-1)

nocells = 25799 # T255=25799, T159=10407

print("\n Reading "+str(nyears)+" years of "+freq+" data for "+str(np.int(nocells))+" cells from "+ifile)

# Get number of columns to skip
skipcol=4
if freq == "yearly" or freq == "monthly_col":
    skipcol=3

if freq == "monthly_col":
    ovars = np.array([os.path.basename(os.path.splitext(ifile)[0])[0:idx_u]])
else:
    ovars = header.split()[skipcol:]

print(" Variables: ",end="")
for v in ovars:
    print(" "+v,end="")
print("\n")

# create netCDF4 dataset
mons = np.array(range(12))+1
if freq == "yearly":
    tsteps = np.array(range(sy,ey+1,1))
elif freq == "monthly_col" or freq == "monthly_row":
    tsteps = np.zeros([tdimsize])
    idx=0
    for y in range(sy,ey+1,1):
        for m in range(12):
            tsteps[idx] = y*100+m+1
            idx+=1
else:
    tsteps = np.zeros([tdimsize])
    idx=0
    for y in range(sy,ey+1,1):
        for m in range(365):
            tsteps[idx] = y*1000+m+1
            idx+=1
        
nodata_value=-1.e16

slabsize=1
# create netCDF4 dataset
#open  dataset
ds = nc.Dataset(ofile,'w',format='NETCDF4')

# define netCDF dimensions
time        = ds.createDimension('time' ,tdimsize)
ncells      = ds.createDimension('ncells',nocells)

# create dim-variables
time      = ds.createVariable('time' ,'i2',('time',))
#CRMlatitude  = ds.createVariable('latitude' ,'f4',('latitude' ,))
#CRMlongitude = ds.createVariable('longitude','f4',('longitude',))

time.units      = tunit

#create variables (test for int v float!!!)
ovl=[]
ovl.append(ds.createVariable("latitude", 'f4',('ncells',),fill_value=nodata_value,))
ovl.append(ds.createVariable("longitude",'f4',('ncells',),fill_value=nodata_value,))
ovl[0].units = "degrees North" 
ovl[1].units = "degrees East" 

idx=2
for v in ovars:
    ovl.append(ds.createVariable(ovars[idx-2],'f4',('time','ncells',),fill_value=nodata_value,))# chunksizes=(1,nocells,),))
        
    ovl[idx].units = "dummy"
    #monthly_ba.units        = "% of total annual burned area"
    #monthly_ba.long_name    = "Monthly burned area climatology" 
    #ovl[idx].grid_mapping = 'latitude longitude'
    #ovl[idx].coordinates     = 'longitude latitude'
    idx+=1

#write coordinates
time[:]      = tsteps[:]
lons = np.zeros(nocells)
lats = np.zeros(nocells)
print(freq)
# Now read the whole file if input is annual
if freq == "yearly" or freq == "monthly_col":
    # Now read whole data in chunks of slabsize cells
    ovtmp=np.zeros([nocells,tdimsize,len(ovars)])
    ovtmp[:,:,:] = np.nan
    
    with open(ifile,"r") as fi:    
        header=fi.readline()
        tcnt=0
        ccnt=0
        for line in fi:
            # check for new coordinate
            if freq == "yearly":
               ovtmp[ccnt,tcnt,:] = np.array(line.split()[skipcol:])
            else:
               ovtmp[ccnt,(12*tcnt):(12*(tcnt+1)),0] = np.array(line.split()[skipcol:])

            if ccnt < nocells: 
               lons[ccnt] = np.float(line.split()[0])
               lats[ccnt] = np.float(line.split()[1])

            ccnt+=1
            if (ccnt % tdimsize*barup) == 0:
                    print("\r {:3.0f}% processed".format(tcnt*100/(tdimsize)),end="")
            if ccnt % nocells == 0:
                tcnt += 1

    ovl[0][:] = lats[:]
    ovl[1][:] = lons[:]
    for x in range(len(ovl)-2):
        for t in range(tdimsize):
            ovl[x+2][t,:] = ovtmp[:,t,x]
            

else:
    # Now read whole data in chunks of slabsize tdimsize
    ovtmp=np.zeros([tdimsize,nocells,len(ovars)])
    with open(ifile,"r") as fi:    
        header=fi.readline()
        tcnt=0
        ccnt=0
        for line in fi:
            # check for new coordinate
            if tcnt % tdimsize == 0:
                tcnt = 0

            ovtmp[tcnt,ccnt,:] = np.array(line.split()[skipcol:])
            
    
            if (tcnt+1) % tdimsize == 0:
                ovl[0][ccnt] = np.float(line.split()[0])
                ovl[1][ccnt] = np.float(line.split()[1])
                ccnt+=1


            if tcnt % tdimsize*barup == 0:
                print("\r {:3.0f}% processed".format(ccnt*100/(nocells)),end="")
            tcnt+=1

    for x in range(2,len(ovl),1):
        ovl[x][:,:] = ovtmp[:,:,x-2]


print("\r {:3.0f}% processed".format(100))
print("\n Output written to: "+ofile)
ds.close()
