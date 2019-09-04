###############################################################################
# imports
###############################################################################
from __future__ import division
import numpy as np
#import h5py
import gdal
from osgeo import osr
#import os
from pyproj import Proj, transform

#gdal.AllRegister()
###############################################################################
# read data
###############################################################################

def Retrieve(HDFFPN,extractD,dstFPN):
    '''
    '''
    print ('open SMAP hdf',HDFFPN)
    '''
    with h5py.File(HDFFPN, mode="r") as f:
        datagrid = '/%(hdffolder)s/%(hdfgrid)s' %{'hdffolder':extractD['hdffolder'],'hdfgrid':extractD['hdfgrid']}

        data = f[datagrid][:]
        #units = f[name].attrs['units']
        #longname = f[name].attrs['long_name']
        _FillValue = f[datagrid].attrs['_FillValue']
        valid_max = f[datagrid].attrs['valid_max']
        valid_min = f[datagrid].attrs['valid_min']
        invalid = np.logical_or(data > valid_max,
                                data < valid_min)
        invalid = np.logical_or(invalid, data == _FillValue)
        data[invalid] = np.nan
        data =  np.ma.masked_where(np.isnan(data), data)
        
        # get geodata
        lonGrid = '/%(llfolder)s/%(longrid)s' %{'llfolder':extractD['lonlatfolder'],'longrid':extractD['longrid']}
        latGrid = '/%(llfolder)s/%(latgrid)s' %{'llfolder':extractD['lonlatfolder'],'latgrid':extractD['latgrid']}
        
        #print ('lonGrid',lonGrid)
        #print ('latGrid',latGrid)
        
        longitude = f[lonGrid][:] #x
        latitude = f[latGrid][:]   #y
        #print ('longitude',longitude)
        #print ('llatitude',latitude)
        #The null (__FillValues) of lon and lat = -9999
        latitude_masked = np.ma.masked_where(latitude==-9999, latitude)
        longitude_masked = np.ma.masked_where(longitude==-9999, longitude)

        # reproject coordinates
        inProj = Proj(init='epsg:4326')
        #EASE GRID 2 = epsg:6933 DOES NOT WORK
        outProj = Proj('+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs')


        longitude,latitude = transform(inProj, outProj, longitude_masked, latitude_masked)
        
        # set up geotransform
        num_cols = data.shape[1]
        num_rows = data.shape[0]
        
        xmin = longitude[longitude!=-9999].min()
        xmax = longitude[longitude!=-9999].max()
        ymin = latitude[latitude!=-9999].min()
        ymax = latitude[latitude!=-9999].max()
    
        xres = (xmax-xmin)/num_cols
        yres = (ymax-ymin)/num_rows
        
        geotransform = (xmin, xres, 0, ymax, 0, -yres)
        
        #print (int(num_cols), int(num_rows))
        # create geotiff
        driver = gdal.GetDriverByName("Gtiff")
        if extractD['celltype'].lower() == 'float32':
            raster = driver.Create(dstFPN,
                                   int(num_cols), int(num_rows),
                                   1, gdal.GDT_Float32)
        else:
            print (extractD['celltype'])
            PLEASADD

        # set geotransform and projection
        
        srs = osr.SpatialReference()
        #srs.ImportFromEPSG(6933)

        #EPSG 6933 for global data, only works with version 4.8 or later for cylindrical PROJ.4
        #wkt_text = r'PROJCS["unnamed",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],TOWGS84[0,0,0,0,0,0,0],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9108"]],AUTHORITY["EPSG","4326"]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["standard_parallel_1",30],PARAMETER["central_meridian",0],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["Meter",1],AUTHORITY["epsg","6933"]]'
        #srs.ImportFromWkt(wkt_text)
        #+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs
        srs.ImportFromProj4('+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs')
        
        #print ('projprint',srs.ExportToWkt())
        raster.SetGeoTransform(geotransform)
        raster.SetProjection(srs.ExportToWkt())
        
        # write data array to raster
        data.data[np.isnan(data.data)]=-9999
        #print (data.data)

        raster.GetRasterBand(1).WriteArray(data.data)
        raster.GetRasterBand(1).SetNoDataValue(-9999)
        raster.FlushCache()
        raster = None
        '''

