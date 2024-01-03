from .utils import utils, Point
import numpy as np
import pandas as pd
import netCDF4 as nc
import os
from datetime import datetime, timedelta
import mikeio
import requests
import time
from mikecore.DfsFile import eumUnit, eumItem
#from mikeio.eum import EUMType, ItemInfo
from mikeio import ItemInfo, EUMType, EUMUnit

class HYCOM:
    def __init__(self, params) -> None:
        """
        This class downloads data from HYCOM model with 3-hourly intervals
        params is a dictionary with the following keys:
            start_date (datetime): The start of the time range to download the data. Make sure that the hour is set to a multiplier of 3 (0, 3, 6, ...)
            end_date (datetime): The end of the time range to download the data.
            variables (list): It can include the following elements: surf_el, water_u, water_v, water_temp, and salinity
            location (list with 2 / 4 elements): [lon, lat], or [North, West, East, South]
            out_folder (raw string): The main directory to store the data
            reset_log (True / False): Flag for reseting the log file
        """
        self.start_date = params["start_date"]
        self.end_date = params["end_date"]
        self.variables = params["variables"]
        self.point = None
        self.bbox = None
        if len(params["location"]) == 2:
            self.location = Point(lon=params["location"][0], lat=params["location"][1])
        elif len(params["location"]) == 4:
            self.location = params["location"]                      # [North, West, East, South]
            self.north = self.location[0]
            self.west = self.location[1]
            self.east = self.location[2]
            self.south = self.location[3]
        self.currd = params["out_folder"]
        utils.mkch(path=self.currd)
        self.reset_log = params["reset_log"]
        self.dates = [self.start_date]
        while self.dates[-1] < self.end_date:
            self.dates.append(self.dates[-1] + timedelta(hours=3))


    def download(self):
        """
        Download the data for the specified information
        """
        utils.create_log(self.currd, self.reset_log)
        start = time.time()
        for date in self.dates:
            self._download_date(date=date) 
        end = time.time()
        run_time = end - start
        print("The data downloaded in " + str(round(run_time, 2)) + " seconds!")
        return
    
    def _download_date(self, date):
        url = self._create_url(date=date)
        fname = datetime.strftime(date, "%Y%m%dT%H")
        timeout = 60
        if self._check_downloaded(date):
            return
        try:
            start_vars = time.time()
            r = requests.get(url, timeout=timeout)
            if r.status_code == 200:
                utils.mkch(" - ".join(self.variables))
                utils.mkch(str(date.year))
                utils.mkch(str(date.month))
                with open(fname+".nc", "wb") as nc:
                    nc.write(r.content)
                end_vars = time.time()
                print(" - ".join(self.variables)+ " - " + fname + " - " + str(round(end_vars - start_vars, 2)) + " seconds")
            else:
                utils.write_log("-".join(self.variables) + " " + fname)
        except:
            utils.write_log("-".join(self.variables) + " " + fname)
        utils.mkch(self.currd)
        return
    
    def _create_url(self, date):
        model, experiment = self._find_experiment_model(date)
        
        if isinstance(self.location, Point):
            if experiment == "expt_53.X":
                base_url = "https://ncss.hycom.org/thredds/ncss/{model}/{experiment}/data/{year}?{vars}&latitude={lat}&longitude={lon}&time={year}-{month}-{day}T{hour}%3A00%3A00Z&vertCoord=&accept=netcdf4"
            else:
                base_url = "https://ncss.hycom.org/thredds/ncss/{model}/{experiment}?{vars}&latitude={lat}&longitude={lon}&time={year}-{month}-{day}T{hour}%3A00%3A00Z&vertCoord=&accept=netcdf4"
        elif isinstance(self.location, list):
            if experiment == "expt_53.X":
                base_url = "https://ncss.hycom.org/thredds/ncss/{model}/{experiment}/data/{year}?{vars}&north={north}&west={west}&east={east}&south={south}&disableProjSubset=on&horizStride=1&time={year}-{month}-{day}T{hour}%3A00%3A00Z&vertCoord=&addLatLon=true&accept=netcdf4"
            else:
                base_url = "https://ncss.hycom.org/thredds/ncss/{model}/{experiment}?{vars}&north={north}&west={west}&east={east}&south={south}&disableProjSubset=on&horizStride=1&time={year}-{month}-{day}T{hour}%3A00%3A00Z&vertCoord=&addLatLon=true&accept=netcdf4"
        
        vars = "&".join(["var="+v for v in self.variables])
        year = str(date.year).zfill(4)
        month = str(date.month).zfill(2)
        day = str(date.day).zfill(2)
        hour = str(date.hour).zfill(2)

        if isinstance(self.location, Point):
            url = base_url.format(model=model, experiment=experiment, vars=vars, year=year, month=month, day=day, hour=hour, lat=str(self.location.lat), lon=str(self.location.lon))
        elif isinstance(self.location, list):
            url = base_url.format(model=model, experiment=experiment, vars=vars, year=year, month=month, day=day, hour=hour, north=self.north, west=self.west, east=self.east, south=self.south)
        return url
    
    def _find_experiment_model(self, date):
        model = "GLBv0.08"
        if date >= datetime(1994, 1, 1, 0, 0, 0) and date <= datetime(2015, 12, 31, 0, 0, 0):
            experiment = "expt_53.X"
        elif date > datetime(2015, 12, 31, 0, 0, 0) and date <= datetime(2016, 4, 30, 0, 0, 0):
            experiment =  "expt_56.3"
        elif date > datetime(2016, 4, 30, 0, 0, 0) and date <= datetime(2017, 1, 31, 0, 0, 0):
            experiment =  "expt_57.2"
        elif date > datetime(2017, 1, 31, 0, 0, 0) and date <= datetime(2017, 5, 31, 0, 0, 0):
            experiment =  "expt_92.8"
        elif date > datetime(2017, 5, 31, 0, 0, 0) and date <= datetime(2017, 9, 30, 0, 0, 0):
            experiment =  "expt_57.7"
        elif date > datetime(2017, 9, 30, 0, 0, 0) and date <= datetime(2017, 12, 31, 0, 0, 0):
            experiment =  "expt_92.9"
        elif date > datetime(2017, 12, 31, 0, 0, 0) and date <= datetime(2020, 2, 18, 0, 0, 0):
            experiment =  "expt_93.0"
        elif date > datetime(2020, 2, 18, 0, 0, 0):
            model = "GLBy0.08"
            experiment =  "expt_93.0"
        else:
            raise ValueError(datetime.strftime(date, "%Y-%m-%d %H:%M:%S") + " is outside the scope of this class!")
        
        return model, experiment
    
    def _check_downloaded(self, date):
        fname = datetime.strftime(date, "%Y%m%dT%H.nc")
        utils.mkch(self.currd)
        utils.mkch(" - ".join(self.variables))
        utils.mkch(str(date.year))
        utils.mkch(str(date.month))
        if fname in os.listdir(os.getcwd()):
            utils.mkch(self.currd)
            return True
        utils.mkch(self.currd)
        return False
    
    def check_missings(self, logged=True):
        steps = []
        date = self.start_date
        while date < self.end_date:
            steps.append(date.strftime("%Y%m%dT%H.nc"))
            date = date + timedelta(hours=3)
        for root, dirs, files in os.walk(self.currd):
            for f in files:
                if f in steps:
                    steps.remove(f)
        
        if logged:
            utils.create_log(self.currd, self.reset_log)
            for elem in steps:
                utils.write_log("-".join(self.variables) + " " + elem.split(".")[0])
        
        output = [datetime.strptime(f, "%Y%m%dT%H.nc") for f in steps]
        return output
    
    def download_missings(self):
        missings = open("log.txt", "r").readlines()
        utils.create_log(self.currd, True)
        for l in missings:
            date = datetime.strptime(l.split(" ")[1], "%Y%m%dT%H\n")
            self._download_date(date=date)

    def to_dfs(self):
        if isinstance(self.location, Point):
            self._to_dfs0()
        else:
            self._to_dfs2()

    def _to_dfs0(self):
        utils.mkch(self.currd)
        utils.create_log(path=self.currd, reset=True, log_name="corrupted files.txt")

        data = {key: None for key in self.variables}
        year = -1
        for date in self.dates:
            utils.mkch(" - ".join(self.variables))
            utils.mkch(str(date.year))
            utils.mkch(str(date.month))
            fname = date.strftime("%Y%m%dT%H.nc")
            loaded_data = self._read_nc_point(fname=fname)
            for var in self.variables:
                try:
                    data[var] = np.append(data[var], loaded_data[var])
                except:
                    data[var] = loaded_data[var]
            if date.year > year:
                year = date.year
                print("Data for {year} is loading!".format(year=str(year)))
            utils.mkch(self.currd) 
        depths = self._read_depth_point()
        dfs = self._create_dfs_point(data, depths)
        itemInfos = self._create_ItemInfo_point(depths)
        utils.mkch(self.currd)
        for var in self.variables:
            dfs[var].to_dfs0(var+".dfs0", items=itemInfos[var])
        


    def _read_nc_point(self, fname):
        loaded = False
        try:
            ds = nc.Dataset(fname)
            loaded = True
        except:
            utils.mkch(self.currd) 
            utils.write_log(error=fname, log_name="corrupted files.txt")
            loaded_data = {}
            for var in self.variables:
                if var == "surf_el":
                    loaded_data[var] = np.nan
                else:
                    loaded_data[var] = np.repeat(np.nan, 40) 
            return loaded_data
        
        if loaded:
            loaded_data = {}
            for var in self.variables:
                if var == "surf_el":
                    data = ds[var][:].data[0]
                    if np.abs(data) > 1000.0:
                        data = np.nan
                    loaded_data[var] = data
                else:
                    loaded_data[var] = ds[var][0, 0, :].data
            utils.mkch(self.currd) 
            return loaded_data
        
    def _read_depth_point(self):
        depths = None
        i = 0
        utils.mkch(self.currd)
        while depths == None:
            utils.mkch(" - ".join(self.variables))
            utils.mkch(str(self.dates[i].year))
            utils.mkch(str(self.dates[i].month))
            fname = self.dates[i].strftime("%Y%m%dT%H.nc")
            try:
                ds = nc.Dataset(fname)
                depths = ds["depth"][0, 0, :].data
                break
            except:
                i = i + 1
        utils.mkch(self.currd)
        return depths
    
    def _create_dfs_point(self, loaded_data, depths):
        scale_factor = 0.001
        dfs = {}
        for var in self.variables:
            if var == "surf_el":
                df = pd.DataFrame(data=loaded_data[var], index=self.dates)
            else:
                df = pd.DataFrame(data=loaded_data[var], index=self.dates, columns=depths)
                df.replace(-30000, np.nan, inplace=True)
            df.sort_index(inplace=True)
            df = df * scale_factor
            dfs[var] = df
        return dfs
    
    def _create_ItemInfo_point(self, depths):
        itemInfos = {}
        if "surf_el" in self.variables:
            itemInfos["surf_el"] = [ItemInfo("Water Level", EUMType.Water_Level, EUMUnit.meter)]
        if "water_u" in self.variables:
            itemInfos["water_u"] = [ItemInfo("Current U "+str(d)+"m", EUMType.u_velocity_component, EUMUnit.meter_per_sec) for d in depths]
        if "water_v" in self.variables:
            itemInfos["water_v"] = [ItemInfo("Current V "+str(d)+"m", EUMType.v_velocity_component, EUMUnit.meter_per_sec) for d in depths]
        if "water_temp" in self.variables:
            itemInfos["water_temp"] = [ItemInfo("Water Temp  "+str(d)+"m", EUMType.Temperature, EUMUnit.degree_Celsius) for d in depths]
        if "salinity" in self.variables:
            itemInfos["salinity"] = [ItemInfo("Salinity "+str(d)+"m", EUMType.Salinity, EUMUnit.PSU) for d in depths]
        return itemInfos
    
    def _to_dfs2(self):
        utils.mkch(self.currd)
        utils.create_log(path=self.currd, reset=True, log_name="corrupted files.txt")

        
        for date in self.dates:
            utils.mkch(" - ".join(self.variables))
            utils.mkch(str(date.year))
            utils.mkch(str(date.month))
            fname = date.strftime("%Y%m%dT%H.nc")
            lon, lat, loaded = self._get_lon_lat(fname)
            if loaded:
                break

        utils.mkch(self.currd)

        
        year = self.dates[0].year
        while year <= self.dates[-1].year: 
            print("Data for {year} is loading!".format(year=str(year)))
            data = {key: None for key in self.variables}
            dates = [d for d in self.dates if d.year == year]
            for date in dates:
                utils.mkch(" - ".join(self.variables))
                utils.mkch(str(date.year))
                utils.mkch(str(date.month))
                fname = date.strftime("%Y%m%dT%H.nc")
                loaded_data = self._read_nc_bbox(fname=fname, lon=lon, lat=lat)
                for var in self.variables:
                    try:
                        data[var] = np.append(data[var], loaded_data[var], axis=0)
                    except:
                        data[var] = loaded_data[var] 
                utils.mkch(self.currd) 
            depths = self._read_depth_bbox()
            das = self._create_das_bbox(data, depths, lon, lat, dates)
            utils.mkch(self.currd)
            for var in self.variables:
                mds = mikeio.Dataset(das[var])
                utils.mkch(var)
                mds.to_dfs(str(year)+".dfs2")
                utils.mkch(self.currd)
            year = year + 1

    def _get_lon_lat(self, fname):
        loaded = False
        try:
            ds = nc.Dataset(fname)
            loaded = True
        except:
            return np.nan, np.nan, loaded
        if loaded:
            lon = ds["lon"][:].data
            if lon[int(len(lon)/2)] > self.east:
                lon = lon - 360
            lat = ds["lat"][:].data
            return lon, lat, True


    def _read_nc_bbox(self, fname, lon, lat):
        loaded = False
        try:
            ds = nc.Dataset(fname)
            loaded = True
        except:
            utils.mkch(self.currd) 
            utils.write_log(error=fname, log_name="corrupted files.txt")
            loaded_data = {}
            for var in self.variables:
                if var == "surf_el":
                    nan_array = np.empty([1, len(lat), len(lon)])
                    nan_array[:] = np.nan
                    loaded_data[var] = nan_array
                else:
                    nan_array = np.empty([1, 40, len(lat), len(lon)])
                    nan_array[:] = np.nan
                    loaded_data[var] = nan_array
            return loaded_data
        
        if loaded:
            loaded_data = {}
            for var in self.variables:
                if var == "surf_el":
                    data = ds[var][:, :, :].data
                    data[np.abs(data) > 1000.0] = np.nan
                    loaded_data[var] = data
                else:
                    data = ds[var][:, :, :, :].data
                    data[np.abs(data) > 1000.0] = np.nan
                    loaded_data[var] = data
            utils.mkch(self.currd) 
            return loaded_data

    def _read_depth_bbox(self):
        depths = None
        i = 0
        utils.mkch(self.currd)
        while depths == None:
            utils.mkch(" - ".join(self.variables))
            utils.mkch(str(self.dates[i].year))
            utils.mkch(str(self.dates[i].month))
            fname = self.dates[i].strftime("%Y%m%dT%H.nc")
            try:
                ds = nc.Dataset(fname)
                depths = ds["depth"][:].data
                break
            except:
                i = i + 1
        utils.mkch(self.currd)
        return depths
    
    def _create_das_bbox(self, loaded_data, depths, lon, lat, dates):
        lon = np.round(lon, 2)
        lat = np.round(lat, 2)
        das = {}
        ItemInfos = self._create_ItemInfo_bbox(depths)
        geometry = mikeio.Grid2D(x=lon, y=lat, projection="LONG/LAT")
        for var in self.variables:
            if var == "surf_el":
                das[var] = [mikeio.DataArray(time=pd.DatetimeIndex(dates), data=loaded_data[var], geometry=geometry, item=ItemInfos[var])]
            else:
                das[var] = [mikeio.DataArray(time=pd.DatetimeIndex(dates), data=loaded_data[var][:, i, :, :], geometry=geometry, item=ItemInfos[var][i]) for i in range(len(depths))]
        return das

    def _create_ItemInfo_bbox(self, depths):
        itemInfos = {}
        if "surf_el" in self.variables:
            itemInfos["surf_el"] = ItemInfo("Water Level", EUMType.Water_Level, EUMUnit.meter)
        if "water_u" in self.variables:
            itemInfos["water_u"] = [ItemInfo("Current U at " + str(d) + "m", EUMType.u_velocity_component, EUMUnit.meter_per_sec) for d in depths]
        if "water_v" in self.variables:
            itemInfos["water_v"] = [ItemInfo("Current V at " + str(d) + "m", EUMType.v_velocity_component, EUMUnit.meter_per_sec) for d in depths]
        if "water_temp" in self.variables:
            itemInfos["water_temp"] = [ItemInfo("Water Temp at " + str(d) + "m", EUMType.Temperature, EUMUnit.degree_Celsius) for d in depths]
        if "salinity" in self.variables:
            itemInfos["salinity"] = [ItemInfo("Salinity at " + str(d) + "m", EUMType.Salinity, EUMUnit.PSU) for d in depths]
        return itemInfos
    
    def dfs22dfs0(self, point):
        if len(self.location) == 2:
            raise SyntaxError("dfs22dfs0 is only available for the cases with bounding box!")
            return
        utils.mkch(self.currd)
        for var in self.variables:
            try:
                os.chdir(var)
                os.chdir(self.currd)
            except:
                raise SyntaxError("You should call to_dfs function first before being able to extract point data")
                return
        os.chdir(self.currd)
        if isinstance(point, list):
            point = Point(lon=point[0], lat=point[1])
        for var in self.variables:
            os.chdir(var)
            for year in range(self.start_date.year, self.end_date.year+1):
                tmp = mikeio.read(str(year) + ".dfs2").interp(x=point.lon, y=point.lat)
                try:
                    df = mikeio.Dataset.concat([df, tmp])
                except:
                    df = tmp
            os.chdir(self.currd)
            df.to_dfs(var+"_"+str(point)+".dfs0")
            del df
