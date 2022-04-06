import ee
from osgeo import gdal, gdal_array
import numpy as np
import matplotlib
from datetime import datetime

from landsat_mosaic import *


"""
ee.Initialize()
forestAge = ee.Image("projects/earthshot-trial/assets/forest_age_00_southamerica").toFloat().divide(10)
abiomass = ee.ImageCollection("NASA/ORNL/biomass_carbon_density/v1").first()
biome = ee.Image("OpenLandMap/PNV/PNV_BIOME-TYPE_BIOME00K_C/v01")
brazil = ee.FeatureCollection("USDOS/LSIB/2017").filter(ee.Filter.eq('COUNTRY_NA', 'Brazil'));

biome = biome.clip(brazil)
forestAge = forestAge.clip(brazil)
abiomass = abiomass.clip(brazil)

print(forestAge.getInfo())
print(abiomass.getInfo())

"""

def getPredictors(targetYear, predictorList, useBrazil300=False):

    # Cast as set
    predictorList = set(predictorList)

    # Static Predictors fixed to 2010
    # Biomes
    predictors = ee.Image("OpenLandMap/PNV/PNV_BIOME-TYPE_BIOME00K_C/v01").select([0], ['biome']).unmask()
    
    # Forest Age
    forestAge = (ee.Image("projects/earthshot-trial/assets/forest_age_00_southamerica").unmask()
                    .toFloat().divide(10).select([0], ['forestage']))
    oldgrowth = forestAge.eq(300).select([0],['oldgrowth'])
    forestAge = forestAge.subtract(oldgrowth.multiply(200)).add(targetYear-2010)
    predictors = predictors.addBands(oldgrowth).addBands(forestAge)
    
    """
    predictors = predictors.addBands(
            # Advance forest age in non-oldgrowth areas
            forestAge.add(forestAge.lt(300).multiply(targetYear-2010))
        )
    """
    
    
    # Model is based on data from 5 years Previous:
    targetYear = targetYear - 5 


    # Elevation & Terrain
    terrainPredictors = predictorList.intersection({'elevation', 'aspect', 'slope', 'hillshade'})
    if terrainPredictors:
        predictors = predictors.addBands(
            ee.Terrain.products(ee.Image('USGS/GMTED2010')
                .select(['be75'], ['elevation']))
                .select(list(terrainPredictors))
        )

    #High-frequency predictors
    # EVI
    tightDateFilter = ee.Filter.date('{0}-01-01'.format(targetYear), '{0}-03-01'.format(targetYear))
    if 'EVI' in predictorList:
        predictors = predictors.addBands(
            ee.ImageCollection('Oxford/MAP/EVI_5km_Monthly')
                .filter(tightDateFilter).first()
                .select([0], ['EVI'])
        )

    # Evapotranspiration
    evapPredictors = predictorList.intersection({'GPP', 'Ec', 'Es', 'Ei'})
    if evapPredictors:
        predictors = predictors.addBands(
            ee.ImageCollection("CAS/IGSNRR/PML/V2")
                .filter(tightDateFilter).first()
                .select(list(evapPredictors))
        )

    # Median Landsat composites
    medLandsatPredictors = predictorList.intersection({'B1','B2','B3','B4','B5','B6','B7','NDVI','NBR','NDWI'})
    if medLandsatPredictors:
            
        # if we're exporting at 300 in Brazil, there's some pre-cooked maps
        if useBrazil300 and targetYear in [2005,2010,2015,2020]:
            SR = ee.Image.loadGeoTIFF('gs://earthshot-trial-eeoutputs/MedianComposite_Brazil_300m_{0}-{1}.tif'.format(targetYear-1, targetYear))
        else:
            SR = getMedianMosaic('{0}-01-01'.format(targetYear-1), '{0}-12-31'.format(targetYear))
        
        # Add normalized differences (only actually get calculated if they're used)
        if predictorList.intersection({'NDVI','NBR','NDWI'}):
            
            nd = ((SR.normalizedDifference(['B4','B3']).select([0],['NDVI']))
                    .addBands(SR.normalizedDifference(['B4','B7']).select([0],['NBR']))
                    .addBands(SR.normalizedDifference(['B2','B5']).select([0],['NDWI']))
                ).multiply(10000).toInt16()
            SR = SR.addBands(nd)
            
        predictors = predictors.addBands(SR.select(list(medLandsatPredictors)))

    return predictors
    
    
def getSample(inputs, aoi, size, classBand, scale, split):
    
    # Sample 50 pixels per biome within Brazil at a scale of 300 m
    sample = inputs.stratifiedSample(
                  numPoints= 100,
                  classBand= classBand,
                  region= aoi,
                  scale= scale,
                  geometries= True
                ).merge( 
                inputs.sample(
                  numPixels= 2000,
                  region= aoi,
                  scale= scale,
                  geometries= True
                ))

    sample = sample.randomColumn()

    # Add a random value field to the sample and use it to approximately split 80%
    # of the features into a training set and 20% into a validation set.
    training = sample.filter(ee.Filter.lte('random', split))
    validation = sample.filter(ee.Filter.gt('random', split))
    
    return training, validation
    
    
def calcError(result, target):
    # This join is obnoxious. There's gotta be a simpler way to do this...
    
    # Join these two feature collections
    joined = ee.Join.inner('first','second').apply(result, target, ee.Filter.equals(leftField='system:index',rightField='system:index'))
    v = ee.Number(0)
    error = joined.iterate(lambda p,v:
                ee.Number(v).add(ee.Number(ee.Feature(p.get('first')).get('agb')).subtract(
                       ee.Number(ee.Feature(p.get('second')).get('agb')))
                        .pow(2))
                , 0)
                  
    # Get mean-squared error
    error = ee.Number(error).divide(result.size())
    valError =  error.sqrt()
    
    # Calc R^2
    variance = target.reduceColumns('variance', ['agb']).get('variance')
    r2 = ee.Number(1).subtract(error.divide(variance))
    
    return valError, r2
    
    
def biomassModel(targetYear, aoi, predictorList):
    # Aboveground Biomass
    abiomass = ee.ImageCollection("NASA/ORNL/biomass_carbon_density/v1").first().select('agb').unmask()
    
    #Predictors
    predictors = getPredictors(2010, predictorList, useBrazil300=True)
    trainingSample, valSample = getSample(abiomass.addBands(predictors), aoi, size=1000, classBand='biome', scale=300, split=0.5)
    
    # Train a 100-tree random forest classifier from the training sample.
    model = ee.Classifier.smileRandomForest(100).setOutputMode('regression').train(
                      features= trainingSample,
                      classProperty= 'agb',
                      inputProperties= predictors.bandNames()
                    )

    # Get a confusion matrix and overall accuracy for the validation sample.
    valResult = valSample.classify(model).select(['classification'], ['agb'])
    valTarget = valSample.select('agb')
    valError, valR2 = calcError(valResult, valTarget)
    modelResult =  model.explain().set('ValError', valError).set('R2', valR2);

    # Classify the reflectance image from the trained classifier.
    predictors = getPredictors(targetYear, predictorList)
    predictedBiomass = predictors.classify(model)
    errorMap = predictedBiomass.subtract(abiomass)
    
    mapResult = predictedBiomass.select([0],['predicted']).addBands(abiomass).addBands(errorMap.select([0],['error']))
    
    return modelResult, mapResult, predictors, trainingSample.merge(valSample)
          
          
########################################
          
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
targetYear = 2015
runName = "Brazil{0}_{1}".format(targetYear, datetime.now().strftime("%b%d-%H%M"))
modelResult, mapResult, predictors, samples = biomassModel(targetYear, brazil, predictorList)

# Clip bounds and rescale
modelResult = modelResult.set('importance', ee.Dictionary(modelResult.get('importance')).map(lambda k,i: ee.Number(i).divide(1e6).toInt()))
predictors = predictors.clip(brazil)
mapResult = mapResult.clip(brazil)
"""
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
            description="{0}_{1}_30m_{2}_inputs".format(runName, regionName, targetYear),
            bucket="earthshot-trial-eeoutputs",
            scale=aoi[0],
            maxPixels=1e10,
            region=aoi[1],
            formatOptions= {'cloudOptimized': True}
        ).start()
    
    ee.batch.Export.image.toCloudStorage(
            image=mapResult.multiply(100).toInt16(),
            description="{0}_{1}_30m_{2}_outputs".format(runName, regionName, targetYear),
            bucket="earthshot-trial-eeoutputs",
            scale=aoi[0],
            maxPixels=1e10,
            region=aoi[1],
            formatOptions= {'cloudOptimized': True}
        ).start()
        
# Print errors from model to console
print(runName, modelResult.getInfo())

"""