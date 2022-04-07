import ee
from osgeo import gdal, gdal_array
import numpy as np
import matplotlib
from datetime import datetime

from landsat_mosaic import *
from biomass_model import biomassModel


# Get bounds for our Area of Intrest in Brazil
brazil = ee.FeatureCollection("USDOS/LSIB/2017").filter(ee.Filter.eq('COUNTRY_NA', 'Brazil'));

#BBox takes: west south east north
Regions = { "Xingu" : [30, ee.Geometry.BBox( -53.8, -12.1, -53.3, -11.7 )],  # Xingu River
       "South" : [30, ee.Geometry.BBox( -53.6, -29.7, -53.1, -29.3 )],   # Parque Estadual Quarta Colonia (Southern Tip)
       "Amazon": [30, ee.Geometry.BBox( -63.9,  -4.1, -63.4,  -3.7 )],   #Deep Amazon
       "Bahia" : [30, ee.Geometry.BBox( -42.3, -11.6, -41.8, -11.2 )],   #Near Irece, Bahia
       "Full"  : [300,ee.Geometry.BBox( -75  , -34  , -34  ,   6   )]
   }
   
# Predictors always include 'biome' and 'forestage'
# Available: {'elevation', 'aspect', 'slope', 'hillshade', 'GPP', 'Ec', 'Es', 'Ei', 'B1','B2','B3','B4','B5','B7','NDVI','NBR','NDWI', 'EVI'}
predictorList = {'elevation', 'aspect', 'B3', 'B4', 'B5', 'B7'}
targetYear = 2025
yearOffset = 5
runName = "Brazil{0}_{1}".format(targetYear, datetime.now().strftime("%b%d-%H%M"))
modelResult, mapResult, predictors, samples = biomassModel(targetYear, yearOffset, brazil, predictorList)


# Clip bounds and rescale
modelResult = modelResult.set('importance', ee.Dictionary(modelResult.get('importance')).map(lambda k,i: ee.Number(i).divide(1e6).toInt()))
predictors = predictors.clip(brazil)
mapResult = mapResult.clip(brazil)


# Export samples used to make model
ee.batch.Export.table.toCloudStorage(
        collection=samples,
        description='{0}_Samples_Locations'.format(runName),
        bucket='earthshot-trial-eeoutputs'
        ).start()


# Export Map results
for regionName, aoi in Regions.items():
    ee.batch.Export.image.toCloudStorage(
           image=predictors.toInt16(),
           description="{0}_{1}_{2}m_{3}_inputs".format(runName, regionName, aoi[0], targetYear),
           bucket="earthshot-trial-eeoutputs",
           scale=aoi[0],
           maxPixels=1e10,
           region=aoi[1],
           formatOptions= {'cloudOptimized': True}
       ).start()

    ee.batch.Export.image.toCloudStorage(
           image=mapResult.multiply(100).toInt16(),
           description="{0}_{1}_{2}m_{3}_outputs".format(runName, regionName, aoi[0], targetYear),
           bucket="earthshot-trial-eeoutputs",
           scale=aoi[0],
           maxPixels=1e10,
           region=aoi[1],
           formatOptions= {'cloudOptimized': True}
       ).start()
   
# Print errors from model to console
print(runName, modelResult.getInfo())
