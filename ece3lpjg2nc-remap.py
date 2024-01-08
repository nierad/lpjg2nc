import logging
import datetime
import dateutil.parser

import netCDF4
import numpy as np
import pandas as pd


#====
import sys
import os
import gzip


# Logger object
log = logging.getLogger(__name__)

# Experiment name
exp_name_ = None

# Table root
table_root_ = None

# Reference date i.e. the start date given by user as a command line parameter
ref_date_ = None

# lpjg_path_ is the directory where the data files (LPJG .out-files) are located
lpjg_path_ = None

# ncpath_ is the tmp directory where the temporary netcdf files will be placed
ncpath_ = None
ncpath_created_ = False

target_grid_ = "T255"
gridfile_ = os.path.join(os.path.dirname(__file__), ".", "ingrid_T255_unstructured.txt")

# weights file and remap command to use
# If weights_file_ is empty then remap_command_ will be used and remapping weights will be computed each time the script is called.
# If weights_file_ is not empty then weights_command will be used to generate the weights_file_ if it is not found.
# This speeds up considerably the remapping command.
#weights_file_ = "weights_remapycon"
weights_file_ = ""
weights_command_ = "genycon"
remap_command_ = "remapycon"

# target grid used for the remapping, n80 for T159 or n128 for T255, which is set in coords()
target_res_ = ""

# list of requested entries for the land use axis
landuse_requested_ = []

# the cmor prefix (e.g. CMIP6) is currently needed to treat the possible requests for land use types, but  might be
# unnecessary in the future depending on how much of the request will be handled in already when writing the model
# output
cmor_prefix_ = None

_months = ['jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec']

# various things extracted from Michael.Mischurow out2nc tool: ec_earth.py
grids = {
    80: [18, 25, 36, 40, 45, 54, 60, 64, 72, 72, 80, 90, 96, 100, 108, 120, 120, 128, 135, 144, 144, 150, 160, 160, 180,
         180, 180, 192, 192, 200, 200, 216, 216, 216, 225, 225, 240, 240, 240, 256, 256, 256, 256, 288, 288, 288, 288,
         288, 288, 288, 288, 288, 300, 300, 300, 300, 320, 320, 320, 320, 320, 320, 320, 320, 320, 320, 320, 320, 320,
         320, 320, 320, 320, 320, 320, 320, 320, 320, 320, 320],
    128: [18, 25, 36, 40, 45, 50, 60, 64, 72, 72, 80, 90, 90, 100, 108, 120, 120, 125, 128, 144, 144, 150, 160, 160,
          180, 180, 180, 192, 192, 200, 216, 216, 216, 225, 240, 240, 240, 250, 250, 256, 270, 270, 288, 288, 288, 300,
          300, 320, 320, 320, 320, 324, 360, 360, 360, 360, 360, 360, 360, 375, 375, 375, 375, 384, 384, 400, 400, 400,
          400, 405, 432, 432, 432, 432, 432, 432, 432, 450, 450, 450, 450, 450, 480, 480, 480, 480, 480, 480, 480, 480,
          480, 480, 486, 486, 486, 500, 500, 500, 500, 500, 500, 500, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512,
          512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512],
}
grids = {i: j + j[::-1] for i, j in grids.items()}


def rnd(x, digits=3):
    return round(x, digits)


def coords(df, root, meta):
    global gridfile_, target_grid_, target_res_
    # deg is 128 in N128
    # common deg: 32, 48, 80, 128, 160, 200, 256, 320, 400, 512, 640
    # correspondence to spectral truncation:
    # t159 = n80; t255 = n128; t319 = n160; t639 = n320; t1279 = n640
    # i.e. t(2*X -1) = nX
    # number of longitudes in the regular grid: deg * 4
    # At deg >= 319 polar correction might have to be applied (see Courtier and Naughton, 1994)
    deg = 128
    target_res_ = "n128"
    grid_size = len(df)
    if grid_size == 10407:
        target_grid_ = "T159"
        target_res_ = "n80"
        deg = 80
        #gridfile_ = os.path.join(os.path.dirname(__file__), "resources/lpjg-grid-content", "ingrid_T159_unstructured.txt")
        gridfile_ = os.path.join(os.path.dirname(__file__), "", "ingrid_T159_unstructured.txt")
    elif grid_size != 25799:
        log.error("Current grid with %i cells is not supported!", grid_size)
        exit(-1)

    lons = [lon for num in grids[deg] for lon in np.linspace(0, 360, num, False)]
    x, w = np.polynomial.legendre.leggauss(deg * 2)
    lats = np.arcsin(x) * 180 / -np.pi
    lats = [lats[i] for i, n in enumerate(grids[deg]) for _ in range(n)]

    if 'i' not in root.dimensions:
        root.createDimension('i', len(lons))
        root.createDimension('j', 1)

        #latitude = root.createVariable('lat', 'f4', ('j', 'i'))
        #latitude.standard_name = 'latitude'
        #latitude.long_name = 'latitude coordinate'
        #latitude.units = 'degrees_north'
        ## latitude.bounds = 'lat_vertices'
        #latitude[:] = lats

        #longitude = root.createVariable('lon', 'f4', ('j', 'i'))
        #longitude.standard_name = 'longitude'
        #longitude.long_name = 'longitude coordinate'
        #longitude.units = 'degrees_east'
        ## longitude.bounds = 'lon_vertices'
        #longitude[:] = lons

    #run_lons = [rnd(i) for i in (df.index.levels[0].values + 360.0) % 360.0]
    #run_lats = [rnd(i) for i in df.index.levels[1]]
    run_lons = df.index.levels[0].values
    run_lats = df.index.levels[1].values

    df.index.set_levels([run_lons, run_lats], inplace=True)
    df = df.reindex([(rnd(i), rnd(j)) for i, j in zip(lons, lats)], fill_value=meta['missing'])
    return df, ('j', 'i')


# TODO: if LPJG data that has been run on the regular grid is also used, the corresponding coords function
# from Michael Mischurow's regular.py should be added here (possibly with some modifications)


# Executes the processing loop.
# used the nemo2cmor.py execute as template
def get_lpj_freq(frequency):
    if frequency == "yr" or frequency == "yrPt":
        return "yearly"
    if frequency == "mon":
        return "monthly"
    if frequency == "day":
        return "daily"
    return None


# checks that the time resolution in the .out data file matches the requested frequency
def check_time_resolution(lpjgfile, freq):
    with open(lpjgfile) as f:
        header = next(f).lower().split()
    if freq == "mon":
        return 'mth' in header or header[-12:] == _months
    elif freq == "day":
        return 'day' in header
    elif freq.startswith("yr"):
        # to find out if it is yearly data have to check that it is neither monthly nor daily
        if 'mth' in header or header[-12:] == _months or 'day' in header:
            return False
        else:
            return True
    else:
        return False  # LPJ-Guess only supports yearly, monthly or daily time resolutions


# Returns first and last year present in the .out data file
def find_timespan(lpjgfile):
    df = pd.read_csv(lpjgfile, delim_whitespace=True, usecols=['Year'], dtype=np.int32)

    firstyr = df['Year'].min()
    lastyr = df['Year'].max()

    return firstyr, lastyr


# Divides the .out file by year to a set of temporary files This approach was chosen to avoid problems with trying to
#  keep huge amounts of data in the memory as the original .out-files can have even >100 years worth of daily data
def divide_years(lpjgfile, firstyr, lastyr, outname):
    files = {}
    filenames = []
    with open(lpjgfile) as f:
        header = next(f)
        # create the yearly files and write header to each
        for yr in range(firstyr, lastyr + 1):
            fname = os.path.join(ncpath_, outname + "_" + str(yr) + ".out")
            filenames.append(fname)
            files[yr] = open(fname, 'w')
            files[yr].write(header)

        # assign the data lines in the .out-file to correct yearly file
        for line in f:
            yr = int(line.split()[2])
            if yr < firstyr:
                continue
            files[yr].write(line)

    for yr in files:
        files[yr].close()

    return filenames


def run_exit(command):

    print(command)
    ret = os.WEXITSTATUS(os.system(command))
    #ret = 0
    if ret!=0:
        print("Error " + str(ret) + " when running command " + command)
        sys.exit(ret)

# this function builds upon a combination of _get and save_nc functions from the out2nc.py tool originally by Michael
#  Mischurow
def create_lpjg_netcdf(freq, inputfile, outname, outfile, outdims):
    global ncpath_, gridfile_, target_res_

    # checks for additional dimensions besides lon,lat&time (for those dimensions where the dimension actually exists
    #  in lpjg data)
    is_land_use = "landUse" in outdims.split()
    is_veg_type = "vegtype" in outdims.split()
    is_sdepth = "sdepth" in outdims.split()
    is_extra = "extra" in outdims.split()
    extras = None

    # assigns a flag to handle two different possible monthly LPJ-Guess formats
    months_as_cols = False
    if freq == "mon":
        with open(inputfile) as f:
            header = next(f).lower().split()
            months_as_cols = header[-12:] == _months

    if freq == "mon" and not months_as_cols:
        idx_col = [0, 1, 2, 3]
    elif freq == "day":
        idx_col = [0, 1]
    else:
        idx_col = [0, 1, 2]

    #df = pd.read_csv(inputfile, delim_whitespace=True, index_col=idx_col, dtype=np.float64, compression='infer')
    df = pd.read_csv(inputfile, delim_whitespace=True, dtype=np.float64, compression='infer')
    df.rename(columns=lambda x: x.lower(), inplace=True)

    # round lon and lat to make sure duplicates are removed, this can happen with e.g. the ssp126 lu files with values of -78.3333 instead of -78.33333
    df['lon'] = rnd( (df['lon'] + 360.0) % 360.0 )
    df['lat'] = rnd( df['lat'] )

    # set index after having rounded lon and lat
    idx_names = [df.columns[i] for i in idx_col ]
    df.set_index(idx_names, inplace=True, drop=True)

    if is_land_use:
        # NOTE: The following treatment of landuse types is likely to change depending on how the lut data requests
        # will be treated when creating the .out-files
        if not landuse_requested_:  # no specific request for land use types, pick all types present in the .out-file
            landuse_types = list(df.columns.values)
        elif cmor_prefix_ == "CMIP6":
            # NOTE: the land use files in the .out-files should match the CMIP6 requested ones (in content if not in
            # name) for future CMIP6 runs this is just a temporary placeholder solution for testing purposes!
            landuse_types = ['psl', 'pst', 'crp', 'urb']
        else:
            # for now skip the variable entirely if there is not exact matches in the .out-file for all the requested
            #  landuse types (except for CMIP6-case of course)
            colnames = list(df.columns.values)
            for lut in landuse_requested_:
                if lut not in colnames:
                    return None
            landuse_types = landuse_requested_

        df_list = []
        for lut in range(len(landuse_types)):
            colname = landuse_types[lut]
            df_list.append(get_lpjg_datacolumn(df, freq, colname, months_as_cols))

    elif is_veg_type:
        pfts = list(df.columns.values)
        df_list = []
        for p in range(len(pfts)):
            df_list.append(get_lpjg_datacolumn(df, freq, pfts[p], months_as_cols))

    elif is_sdepth:
        depths = list(df.columns.values)
        if "year" in depths:
            if freq == "mon" or freq == "day":
                depths = list(depths[2:])
            else:
                depths = list(depths[1:])
        df_list = []
        for sd in range(len(depths)):
            df_list.append(get_lpjg_datacolumn(df, freq, depths[sd], months_as_cols))

    elif is_extra:
        extras = list(df.columns.values)
        df_list = []
        for p in range(len(extras)):
            print("reading column "+extras[p])
            df_list.append(get_lpjg_datacolumn(df, freq, extras[p], months_as_cols))

    else:  # regular variable
        colname = ""
        if not months_as_cols:
            if "total" not in list(df.columns.values):
                return None
            else:
                colname = "total"
        df = get_lpjg_datacolumn(df, freq, colname, months_as_cols)
        df_list = [df]

    if freq.startswith("yr"):
        str_year=str(int(df_list[0].columns[0]))
    else:
        str_year=str(int(df_list[0].columns[0][1]))

    #log.info( "Creating lpjg netcdf file for variable " + outname + " for year " + str_year )
    print( "Creating lpjg netcdf file " + outname + " for year " + str_year )

    #ncfile = os.path.join(ncpath_, outname + "_" + freq + "_" + str_year + ".nc")
    ncfile = outfile
    # Note that ncfile could be named anything, it will be deleted later and the cmorization takes care of proper
    # naming conventions for the final file

    # temporary netcdf file name (will be removed after remapping is done)
    #temp_ncfile = os.path.join(ncpath_, 'LPJGtemp.nc')
    temp_ncfile = os.path.join(ncpath_, outname+'_LPJGtemp')
    root = netCDF4.Dataset(temp_ncfile, 'w')  # now format is NETCDF4

    root.createDimension('time', None)
    timev = root.createVariable('time', 'f4', ('time',))
    refyear = int(ref_date_.year)
    if freq == "mon":
        curyear, tres = int(df_list[0].columns[0][1]), 'month'
        t_since_fyear = (curyear - refyear) * 12
    elif freq == "day":
        curyear, tres = int(df_list[0].columns[0][1]), 'day'
        t_since_fyear = (date(curyear, 1, 1) - date(refyear, 1, 1)).days
    else:
        curyear, tres = int(df_list[0].columns[0]), 'year'
        t_since_fyear = curyear - refyear
    timev[:] = np.arange(t_since_fyear, t_since_fyear + df_list[0].shape[1])
    timev.units = '{}s since {}-01-01'.format(tres, refyear)
    timev.calendar = "proleptic_gregorian"

    meta = {"missing": 1.e+20}  # the missing/fill value could/should be taken from the target header info if available
    # and does not need to be in a meta dict since coords only needs the fillvalue anyway, but do it like this (i.e.
    # out2nc-style) for now

    N_dfs = len(df_list)
    df_normalised = []
    dimensions = []

    print("Converting to the global gaussian grid")

    for l in range(N_dfs):
        # TODO: if different LPJG grids possible you need an if-check here to choose which function is called
        df_out, dimensions = coords(df_list[l], root, meta)
        df_normalised.append(df_out)

    if N_dfs == 1:
        dimensions = 'time', dimensions[0], dimensions[1]

        variable = root.createVariable(outname, 'f4', dimensions, zlib=True,
                                       shuffle=False, complevel=5, fill_value=meta['missing'])
        if outname == "tsl":
            variable[:] = df_normalised[0].values.T  # TODO: see out2nc for what to do here if you have the LPJG regular grid
        else:
            dumvar = df_normalised[0].values.T
            #variable[:] = np.where(dumvar < 1.e+20, dumvar, meta['missing'])   # TODO: see out2nc for what to do here if you have the LPJG regular grid
            variable[:] = np.where(dumvar < 1.e+20, dumvar, 0.)   # TODO: see out2nc for what to do here if you have the LPJG regular grid

    elif is_extra:
        dimensions = 'time', dimensions[0], dimensions[1]
        for l in range(N_dfs):
            print("creating variable "+extras[l])
            variable = root.createVariable(extras[l], 'f4', dimensions, zlib=True,
                                           shuffle=False, complevel=2, fill_value=meta['missing'])
            if outname == "tsl":
                variable[:] = df_normalised[l].values.T  # TODO: see out2nc for what to do here if you have the LPJG regular grid
            else:
                dumvar = df_normalised[l].values.T
                #variable[:] = np.where(dumvar < meta['missing'], dumvar, meta['missing'])   # TODO: see out2nc for what to do here if you have the LPJG regular grid
                variable[:] = np.where(dumvar < meta['missing'], dumvar, 0.)   # TODO: see out2nc for what to do here if you have the LPJG regular grid

    else:
        root.createDimension('fourthdim', N_dfs)

        dimensions = 'time', 'fourthdim', dimensions[0], dimensions[1]
        variable = root.createVariable(outname, 'f4', dimensions, zlib=True,
                                       shuffle=False, complevel=5, fill_value=meta['missing'])
        for l in range(N_dfs):
            if outname == "tsl":
                variable[:, l, :, :] = df_normalised[l].values.T  # TODO: see out2nc for what to do here if you have the LPJG regular grid
            else:
                dumvar = df_normalised[l].values.T
                #variable[:, l, :, :] = np.where(dumvar < 1.e+20, dumvar, meta['missing'])   # TODO: see out2nc for what to do here if you have the LPJG regular grid
                variable[:, l, :, :] = np.where(dumvar < 1.e+20, dumvar, 0.)   # TODO: see out2nc for what to do here if you have the LPJG regular grid

    root.sync()
    root.close()

    # do the remapping
#    cdo = Cdo()

#    if target_grid_ == "T159":
#        cdo.remapycon('n80', input="-setgrid," + gridfile_ + " " + temp_ncfile,
#                      output=ncfile)
#    else:
#        cdo.remapycon('n128',input="-setgrid," + gridfile_ + " " + temp_ncfile,
#                      output=ncfile)  # TODO: add remapping for possible other grids

    if weights_file_ != "":
        if not os.path.exists(weights_file_):
            run_exit("cdo -O -L " + weights_command_ + "," + target_res_ + " -setgrid," + gridfile_ + " " + temp_ncfile + " " + weights_file_)
        run_exit("cdo -O -L -f nc4c -z zip_2 -remap," + target_res_ + "," + weights_file_ + " -setgrid," + gridfile_ + " " + temp_ncfile + " " +ncfile)
    else:
        run_exit("cdo -O -L -f nc4c -z zip_2 -" + remap_command_ + "," + target_res_ + " -setgrid," + gridfile_ + " " + temp_ncfile + " " +ncfile)

    os.remove(temp_ncfile)

    return ncfile


# Extracts single column from the .out-file
def get_lpjg_datacolumn(df, freq, colname, months_as_cols):
    if freq == "day":
        # create a single time column so that extra days won't be added to 
        # the time axis (if there are both leap and non-leap years)
        # Time axis needs to be modified on first call 
        if "year" in list(df.columns.values):
            df['timecolumn'] = df['year'] + 0.001 * df['day']
            df.drop(columns=['year', 'day'], inplace=True)
            df.set_index('timecolumn', append=True, inplace=True)
        df = df[[colname]]
        df = df.unstack()
         
    elif freq.startswith("yr"):
        df = df.pop(colname)
        df = df.unstack()
    elif freq == "mon":
        if months_as_cols:
            df = df.unstack()
            df = df.reindex(sorted(df.columns, key=(lambda x: (x[1], _months.index(x[0])))),
                            axis=1, copy=False)
        else:
            df = df.pop(colname)
            df = df.unstack().unstack()
            df = df.reindex(sorted(df.columns, key=(lambda x: (x[1], x[0]))),
                            axis=1, copy=False)
    return df



# setup default values for ece2cmor code
outdims = "longitude latitude time"
ncpath_ = ""
ncpath_created_ = ""
ref_date_ = datetime.datetime.combine(dateutil.parser.parse("1850-01-01"), datetime.datetime.min.time())

ifile=sys.argv[1]
ofile=sys.argv[2]

print(ifile)
print(ofile)

#Check file size and time period
lcnt=0
gcnt=1

#with open(ifile,"r") as fi:
#    header=fi.readline()
#    lo,la,y=fi.readline().split()[:3]
#    sy=np.int(y)
try:
    if ifile.endswith(".gz"):
        fi=gzip.open(ifile,"r")
    else:
        fi=open(ifile,"r")
except Error:
    print("Error opening ifile "+ifile)
    sys.exit(1)

header=fi.readline()
lo,la,y=fi.readline().split()[:3]
sy=np.int(y)

fi.close()

nyears = 1
ey=sy
# Check time-range and output-frequency
try:
    idx_u = os.path.basename(os.path.splitext(ifile)[0]).index("_")
    freq  = os.path.basename(os.path.splitext(ifile)[0])[idx_u+1:]
except ValueError:
    freq = "yearly"
if freq != "monthly" or freq != "yearly":
    freq = "yearly"

if freq == "monthly":
    if header.split()[3]=="Jan" and header.split()[14]=="Dec":
        freq     = "monthly_col"
        tdimsize = nyears*12
        tunit    = "months since 01 "+str(sy)
        barup    = 100
        cmorfreq = "mon"
    else:
        freq = "monthly_row"
        tdimsize = nyears*12
        tunit    = "months since 01 "+str(sy)
        #    ncells   = (lcnt-1)/(12*nyears)
        cmorfreq = "mon"
        barup    = 100
        #print("error monthly_row not supported")
        #sys.exit(-1)
        outdims="longitude latitude extra time"
elif freq == "yearly":
    if header.split()[3]=="Total":
        freq     = "yearly_col"
    else:
        freq     = "yearly_row"
        outdims="longitude latitude extra time"
    tdimsize = nyears
    tunit    = "yearss since "+str(sy)
    #ncells   = (lcnt-1)/nyears
    barup    = 1000
    cmorfreq = "yr"
elif freq == "daily":
    tdimsize = nyears*365
    tunit    = "days since 01 January "+str(sy)
#    ncells   = (lcnt-1)/(365*nyears)
    barup    = 1
    cmorfreq = "day"
else:
    sys.exit(-1)

#print("\n Reading "+str(nyears)+" years of "+freq+" data for "+str(np.int(nocells))+" cells from "+ifile)

# Get number of columns to skip
skipcol=4

if freq == "monthly_col":
    skipcol=3
    #ovars = np.array([os.path.basename(os.path.splitext(ifile)[0])[0:idx_u]])
    outname = ovars[0]
else:
    #ovars = header.split()[skipcol:]
    outname = ofile

#print("ovars:")
#print(ovars)

ncfile = create_lpjg_netcdf(cmorfreq, ifile, outname, ofile, outdims)
