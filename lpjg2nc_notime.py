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

def latlon2ixjy (LAT,LON,nx,ny,mtype='array'):

    ex = 0
    chk = np.array(LAT,dtype=float,ndmin=1)
    if np.any(chk<-90.) or np.any(chk>90.):
        print("Excess LAT in latlon2ixjy")
        ex = 1
    chk = np.array(LON,dtype=float,ndmin=1)
    if np.any(chk<-180.) or np.any(chk>180.):
        print("Excess LON in latlon2ixjy")
        ex = 1
    if ex == 1:
        sys.exit(-1)

    valx = np.array(np.rint(((LON- 0.25) + 180.) / 360. * nx))
    valy = np.array(np.rint(((LAT- 0.25) +  90.) / 180. * ny))

    XY = np.copy(valx)    
    XY = np.append([XY], [valy],axis=0).reshape(2,-1)

    if mtype == 'array':    
        XY = np.fix(XY).astype(int)
        return XY
    elif mtype == 'mask':
        MASK = np.zeros((nx,ny))
        MASK[np.fix(XY[0]).astype(int),np.fix(XY[1]).astype(int)] = 1        
        return MASK.T
    elif mtype == "single":
        return int(XY[0,0]),int(XY[1,0])
    else:
        print("Wrong mtype in latlon2ixjy. Either 'array' or 'mask'!")
        sys.exit(-1)

ifile=sys.argv[1]
ofile=sys.argv[2]
slabsize=1
if len(sys.argv)>3:
    slabsize = sys.argv[3]

#Check file size and time period
lcnt=0
gcnt=1

sy = 0
ey = 0
has_time=1
with open(ifile,"r") as fi:
    header=fi.readline()
    try:
        header.index("Year")
    except:
        has_time=0
    lo,la,y=fi.readline().split()[:3]
    if has_time:
        sy=int(y)
    else:
        sy = 0

    lcnt+=2
    tflag=1
    for line in fi:
        if tflag:
            lox,lax,y=line.split()[:3]
            if lox != lo or lax != la:
                tflag=0
                break
            else:
                gcnt+=1
                if has_time:
                    ey=int(y)
        lcnt+=1

nyears = ey-sy+1
print(header)
# Check time-range and output-frequency
if gcnt/nyears == 1:
    freq     = "yearly"
    tdimsize = nyears
    tunit    = "yearss since "+str(sy)
    ncells   = (lcnt-1)/nyears
    barup    = 1000
    if has_time:
        if header.split()[3]=="Jan" and header.split()[14]=="Dec":
            freq     = "monthly_col"
            tdimsize = nyears*12
            tunit    = "months since 01 "+str(sy)
            ncells   = (lcnt-1)/nyears
            barup    = 100
elif gcnt/nyears == 12:
    freq = "monthly_row"
    tdimsize = nyears*12
    tunit    = "months since 01 "+str(sy)
    ncells   = (lcnt-1)/(12*nyears)
    barup    = 100
elif gcnt/nyears == 365:
    freq = "daily"
    tdimsize = nyears*365
    tunit    = "days since 01 January "+str(sy)
    ncells   = (lcnt-1)/(365*nyears)
    barup    = 1
else:
    sys.exit(-1)

ncells = 59191

print("\n Reading "+str(nyears)+" years of "+freq+" data for "+str(int(ncells))+" cells from "+ifile)

# Get number of columns to skip
skipcol=4
if freq == "yearly" or freq == "monthly_col":
    skipcol=3
if has_time == 0:
    skipcol=2
    
if freq == "monthly_col":
    ovars = np.array([os.path.basename(os.path.splitext(ifile)[0])[0:]])
else:
    ovars = header.split()[skipcol:]

print(" Variables: ",end="")
for v in ovars:
    print(" "+v,end="")
print("\n")

# create netCDF4 dataset
nlat=360
nlon=720

lons = np.array(range(nlon))*0.5-179.75
lats = np.array(range(nlat))*0.5- 89.75
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
latitude    = ds.createDimension('latitude',nlat)
longitude   = ds.createDimension('longitude',nlon)

# create dim-variables
time      = ds.createVariable('time' ,'i2',('time',))
latitude  = ds.createVariable('latitude' ,'f4',('latitude' ,))
longitude = ds.createVariable('longitude','f4',('longitude',))

time.units      = tunit
latitude.units  = "degrees North" 
longitude.units = "degrees East" 

#create variables (test for int v float!!!)
ovl=[]
idx=0
for v in ovars:
    if freq == "yearly" or freq == "monthly_col":
        ovl.append(ds.createVariable(ovars[idx],'f4',('time','latitude','longitude',),fill_value=nodata_value, chunksizes=(1,360,720,),))
    else:
        ovl.append(ds.createVariable(ovars[idx],'f4',('time','latitude','longitude',),fill_value=nodata_value, chunksizes=(tdimsize,1,1,),))
        
    ovl[idx].units = "dummy"
    #monthly_ba.units        = "% of total annual burned area"
    #monthly_ba.long_name    = "Monthly burned area climatology" 
    ovl[idx].grid_mapping = 'latitude longitude'
    ovl[idx].coordinates     = 'longitude latitude'
    idx+=1

#write coordinates
longitude[:] = lons[:]
latitude[:]  = lats[:]
time[:]      = tsteps[:]

# Now read the whole file if input is annual
if freq == "yearly" or freq == "monthly_col":
    # Now read whole data in chunks of slabsize cells
    ovtmp=np.zeros([nlat,nlon,tdimsize,len(ovars)])
    ovtmp[:,:,:,:] = np.nan
    with open(ifile,"r") as fi:    
        header=fi.readline()
        tcnt=0
        cnt =0
        ccnt=0
        for line in fi:
            # check for new coordinate
            tborder = tdimsize
            if freq == "monthly_col":
                tborder = tdimsize / 12
            if cnt % tborder == 0:
                lox,lay = latlon2ixjy(float(line.split()[1]),float(line.split()[0]),nlon,nlat,mtype="single")
                cnt = 0
        
            if freq == "yearly":
                ovtmp[lay,lox,cnt,:] = np.array(line.split()[skipcol:])
            else:
                ovtmp[lay,lox,(12*cnt):(12*(cnt+1)),0] = line.split()[skipcol:]
            cnt+=1
            tcnt+=1
            if (tcnt) % (tdimsize*barup) == 0:
                    print("\r {:3.0f}% processed".format(tcnt*100/(tdimsize*ncells)),end="")

    for x in range(len(ovl)):
        for t in range(tdimsize):
            ovl[x][t,:,:] = ovtmp[:,:,t,x]
            
else:
    # Now read whole data in chunks of slabsize tdimsize
    ovtmp=np.zeros([tdimsize,len(ovars)])
    with open(ifile,"r") as fi:    
        header=fi.readline()
        tcnt=0
        cnt =0
    
        for line in fi:
            # check for new coordinate
            if cnt % tdimsize == 0:
                #print(float(line.split()[1]),float(line.split()[0]))
                lox, lay = latlon2ixjy (float(line.split()[1]),float(line.split()[0]),nlon,nlat,mtype="single")
                ovtmp[:,:] = np.nan
                cnt = 0

            ovtmp[cnt,:] = np.array(line.split()[skipcol:])
            
    
            if (cnt+1) % tdimsize == 0:
                for x in range(len(ovl)):
                    ovl[x][:,lay,lox] = ovtmp[:,x]

            if (tcnt) % (tdimsize*barup) == 0:
                print("\r {:3.0f}% processed".format(tcnt*100/(tdimsize*ncells)),end="")
            cnt+=1
            tcnt+=1

print("\r {:3.0f}% processed".format(100))
print("\n Output written to: "+ofile)
ds.close()
