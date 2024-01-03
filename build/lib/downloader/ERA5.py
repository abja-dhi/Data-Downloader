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
            
        # Concatenate all years 
        utils.mkch(self.currd_dfs)
        files = [f for f in os.listdir(os.getcwd()) if ".dfs2" in f]
        df = mikeio.read(files[0])
        for f in files[1:]:
            df = mikeio.Dataset.concat([df, mikeio.read(f)])
        utils.mkch(self.currd)
        df.to_dfs("ERA5.dfs2")
        end = time.time()
        print("All files converted to dfs2 in " + str(round(end-start, 2)) + " seconds!")

    def _nc2dfs(self, fname):
        start_y = time.time()
        times, u10, v10, lons, lats = self._parse_nc(fname)
        times = pd.DatetimeIndex(times)
        geometry = mikeio.Grid2D(x=lons, y=lats, projection="LONG/LAT")
        u_da = mikeio.DataArray(data=u10,time=times, geometry=geometry, item=ItemInfo("Wind U", EUMType.Wind_Velocity, EUMUnit.meter_per_sec))
        v_da = mikeio.DataArray(data=v10,time=times, geometry=geometry, item=ItemInfo("Wind V", EUMType.Wind_Velocity, EUMUnit.meter_per_sec))
        mds = mikeio.Dataset([u_da, v_da])
        utils.mkch(self.currd_dfs)
        mds.to_dfs(fname.replace(".nc", ".dfs2"))
        end_y = time.time()
        print(fname.replace(".nc", ".dfs2") + " converted in " + str(round(end_y - start_y, 2)) + " seconds!")

    def _parse_nc(self, fname):
        utils.mkch(self.currd_nc)
        ds = nc.Dataset(fname)
        u10 = np.flip(ds["u10"][:, :, :].data, axis=1)
        v10 = np.flip(ds["v10"][:, :, :].data, axis=1)
        ttt = ds["time"][:].data
        times = []
        for t in ttt:
            times.append(datetime(1900, 1, 1, 0, 0, 0) + timedelta(hours=int(t)))
        lons = ds["longitude"][:].data
        lats = np.flip(ds["latitude"][:].data)

        return times, u10, v10, lons, lats

    def extract_point(self, point, name="", interp=False, vars="uv"):
        if name == "":
            name = str(point)
        utils.mkch(self.currd)
        df = mikeio.read("ERA5.dfs2")
        if interp:
            df_extract = df.interp(x=point.lon, y=point.lat, n_nearest=4)
        else:
            df_extract = df.sel(x=point.lon, y=point.lat)
        if vars == "uv":
            df_extract.to_dfs("ERA5_"+ name + ".dfs0")
        elif vars == "spddir":
            data, items = self._uv2spddir(df_extract)
            data.to_dfs0("ERA5_"+ name + ".dfs0", items=items)
        
    def _uv2spddir_from_mikeio(self, df):
        u = df[0].values
        v = df[1].values
        spd, dir = uv2spddir(u, v, 'from')
        times = df.time
        data = pd.DataFrame({"WS": spd, "WD": dir}, index=times)
        ws_item = ItemInfo("Wind Speed", EUMType.Wind_speed, EUMUnit.meter_per_sec)
        wd_item = ItemInfo("Wind Direction", EUMType.Wind_Direction, EUMUnit.degree)
        items = [ws_item, wd_item]
        return data, items
    
    def _uv2spddir_from_dataframe(self, df):
        df["Wind Speed"], df["Wind Direction"] = uv2spddir(df["Wind U"], df["Wind V"], 'from')
        df.drop(["Wind U", "Wind V"], axis=1, inplace=True)
        ws_item = ItemInfo("Wind Speed", EUMType.Wind_speed, EUMUnit.meter_per_sec)
        wd_item = ItemInfo("Wind Direction", EUMType.Wind_Direction, EUMUnit.degree)
        items = [ws_item, wd_item]
        return df, items

    

    def get_averaging(self, fname="", vars='uv', window=120):
        utils.mkch(self.currd)
        if "ERA5" not in fname:
            name = "ERA5_" + fname + ".dfs0"
        df = mikeio.read(name).to_dataframe()
        
        if "spddir" in vars:
            pass
        df_averaged = df.rolling(window=timedelta(minutes=window), center=True).mean()
        U_item = ItemInfo("U velocity", EUMType.Wind_Velocity, EUMUnit.meter_per_sec)
        V_item = ItemInfo("V velocity", EUMType.Wind_Velocity, EUMUnit.meter_per_sec)
        df_averaged.to_dfs0("ERA5_"+ fname + "_" + str(int(window/60)) + "h_averaged.dfs0", items=[U_item, V_item])
        

        
