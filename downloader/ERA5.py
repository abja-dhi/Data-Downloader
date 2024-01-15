from .utils import utils, Point
from p_tools import uv2spddir
import numpy as np
import pandas as pd
import netCDF4 as nc
from datetime import datetime, timedelta
import mikeio
import os
import time
import cdsapi
from mikeio import ItemInfo, EUMType, EUMUnit

class ERA5:
    def __init__(self, params):
        self.start_year = params["start_year"]
        self.end_year = params["end_year"]
        self.variables = params["variables"]
        self.items = self._items_definition()
        self.north = params["bbox"][0]
        self.west = params["bbox"][1]
        self.east = params["bbox"][2]
        self.south = params["bbox"][3]
        self.currd = params["out_folder"]
        self.currd_nc = os.path.join(self.currd, "nc")
        self.currd_dfs = os.path.join(self.currd, "dfs")
        
        utils.mkch(self.currd)

    def _create_dict(self):
        self.dict = {'product_type': 'reanalysis',
                     'format': 'netcdf',
                     'variable': self.variables,
                     'month': [str(m).zfill(2) for m in range(1, 13)],
                     'day': [str(d).zfill(2) for d in range(1, 32)],
                     'time': [datetime.strftime(datetime(2000, 1, 1, t, 0, 0), "%H:%M") for t in range(24)],
                     'area': [self.north, self.west, self.south, self.east]}
    
    def _download(self, year):
        self._create_dict()
        self.dict['year'] = str(year)
        c = cdsapi.Client()
        c.retrieve('reanalysis-era5-single-levels', self.dict, str(year) + '.nc')

    def download(self):
        utils.mkch(self.currd_nc)
        start = time.time()
        for y in range(self.start_year, self.end_year+1):
            start_y = time.time()
            self._download(y)
            end_y = time.time()
            print(str(y) + " data downloaded in " + str(round(end_y - start_y, 2)) + " seconds!")
        end = time.time()
        print("Downloading finished in " + str(round(end - start, 2)) + " seconds!")

    def nc2dfs(self):
        start = time.time()
        files = [f for f in os.listdir(self.currd_nc) if ".nc" in f]
        
        for fname in files:
            self._nc2dfs(fname)
            
        end = time.time()
        print("All files converted to dfs2 in " + str(round(end-start, 2)) + " seconds!")

    def _create_item_info(self, item):
        if item == "u10":
            return ItemInfo("Wind U", EUMType.Wind_Velocity, EUMUnit.meter_per_sec)
        elif item == "v10":
            return ItemInfo("Wind V", EUMType.Wind_Velocity, EUMUnit.meter_per_sec)
        elif item == "t2m":
            return ItemInfo("2m Air Temperature", EUMType.Temperature, EUMUnit.degree_Kelvin)
        elif item == "msl":
            return ItemInfo("Mean Sea Level Pressure", EUMType.Pressure, EUMUnit.pascal)
        elif item == "siconc":
            return ItemInfo("Sea Ice Area Fraction", EUMType.Snow_Cover_Percentage)
        elif item == "sst":
            return ItemInfo("Sea Surface Temperature", EUMType.Temperature, EUMUnit.degree_Kelvin)
        elif item == "sp":
            return ItemInfo("Surface Pressure", EUMType.Pressure, EUMUnit.pascal)
        elif item == "d2m":
            return ItemInfo("2m Dewpoint Temperature", EUMType.Temperature, EUMUnit.degree_Kelvin)
        elif item == "blh":
            return ItemInfo("Boundary Layer Height", EUMType.Boundary_Layer_Thickness, EUMUnit.meter)
        elif item == "chnk":
            return ItemInfo("Charnock", EUMType.Charnock_constant)
        elif item == "zust":
            return ItemInfo("Friction Velocity", EUMType.Velocity, EUMUnit.meter_per_sec)
        elif item == "msdwlwrf":
            return ItemInfo("Mean Surface Downward Long Wave Radiation Flux", EUMType.Radiation_intensity, EUMUnit.watt_per_meter_pow_2)
        elif item == "msdwswrf":
            return ItemInfo("Mean Surface Downward Short Wave Radiation Flux", EUMType.Radiation_intensity, EUMUnit.watt_per_meter_pow_2)
                            

    def _nc2dfs(self, fname):
        start_y = time.time()
        data = self._parse_nc(fname)
        times = pd.DatetimeIndex(data["time"])
        geometry = mikeio.Grid2D(x=data["lon"], y=data["lat"], projection="LONG/LAT")
        das = [mikeio.DataArray(data=data[key], time=times, geometry=geometry, item=self._create_item_info(self.items[key])) for key in self.items.keys()]
        mds = mikeio.Dataset(das)
        utils.mkch(self.currd_dfs)
        mds.to_dfs(fname.replace(".nc", ".dfs2"))
        end_y = time.time()
        print(fname.replace(".nc", ".dfs2") + " converted in " + str(round(end_y - start_y, 2)) + " seconds!")

    def _parse_nc(self, fname):
        utils.mkch(self.currd_nc)
        ds = nc.Dataset(fname)
        data = {key: np.flip(ds[self.items[key]][:, :, :].data, axis=1) for key in self.items.keys()}
        for key in self.items.keys():
            data[key][data[key] == ds[self.items[key]]._FillValue] = np.nan
        data["time"] = [datetime(1900, 1, 1, 0, 0, 0) + timedelta(hours=int(t)) for t in ds["time"][:].data]
        data["lon"] = ds["longitude"][:].data
        data["lat"] = np.flip(ds["latitude"][:].data)
        return data

    def _items_definition(self):
        items = {"10m_u_component_of_wind": "u10",
                 "10m_v_component_of_wind": "v10",
                 "2m_temperature": "t2m",
                 "mean_sea_level_pressure": "msl",
                 "sea_ice_cover": "siconc",
                 "sea_surface_temperature": "sst",
                 "surface_pressure": "sp",
                 "2m_dewpoint_temperature": "d2m",
                 "boundary_layer_height": "blh",
                 "charnock": "chnk",
                 "friction_velocity": "zust",
                 "mean_surface_downward_long_wave_radiation_flux": "msdwlwrf",
                 "mean_surface_downward_short_wave_radiation_flux": "msdwswrf"
                 }
        return {key: items[key] for key in self.variables}
            

    def extract_point(self, point, name="", interp=False):
        if name == "":
            name = str(point)
        utils.mkch(self.currd_dfs)
        files = [f for f in os.listdir(os.getcwd()) if ".dfs2" in f]
        if interp:
            df = mikeio.read(files[0]).interp(x=point.lon, y=point.lat, n_nearest=4)
        else:
            df = mikeio.read(files[0]).sel(x=point.lon, y=point.lat)
        print(files[0].split(".")[0] + " extracted!")
        for f in files[1:]:
            if interp:
                tmp = mikeio.read(f).interp(x=point.lon, y=point.lat, n_nearest=4)
            else:
                tmp = mikeio.read(f).sel(x=point.lon, y=point.lat)
            df = mikeio.Dataset.concat([df, tmp])
            print(f.split(".")[0] + " extracted!")
        utils.mkch(self.currd)
        df.to_dfs("ERA5_"+ name + ".dfs0")
        
    
