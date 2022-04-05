# Build median mosaic.
# Joe Hughes 2022-APR-02
#
# Adapted from LandTrendr.js:
# * @author Zhiqiang Yang (USDA Forest Service)
# * @author Robert Kennedy (Oregon State University)
# * @author Justin Braaten (Oregon State University)
# * @author Joe Hughes (Oregon State University)
# * @author Peter Clary (Oregon State University)
#

import ee
ee.Initialize()


# Apply a cloud, shadow, snow, and/or water mask to an SR image from its QA Band
def maskSRcollection(img, sensor, maskThese):
    maskOptions = {'cloud':1<<3, 'shadow':1<<4, 'snow':1<<5, 'water':1<<7}
    
    # Select the bands and rename similar bands back to the TM names
    if sensor in ['LC08', 'LC09']:
        dat = img.select(['SR_B2','SR_B3','SR_B4','SR_B5','SR_B6','SR_B7'],['B1', 'B2', 'B3', 'B4', 'B5', 'B7']).unmask() #OLI
    else:
        dat = img.select(['SR_B1','SR_B2','SR_B3','SR_B4','SR_B5','SR_B7'],['B1', 'B2', 'B3', 'B4', 'B5', 'B7']).unmask()  #TM/ETM
    
    # Mask specified obstructions
    if maskThese: 
        mask = 0
        for m in maskThese:
            mask += maskOptions[m]
        dat = dat.mask(img.select('QA_PIXEL').bitwiseAnd(mask).eq(0))
        
    return dat
    
# Collect all of the masked imagery between dates and within the aoi for a sensor
def getSRcollection(startDate, endDate, aoi, sensor, maskThese):
    maskThese = [m.lower() for m in maskThese]
    
    # get a landsat collection for given year, day range, and sensor
    srCollection = (ee.ImageCollection('LANDSAT/'+ sensor + '/C02/T1_L2')
                     .filterDate(startDate, endDate)
                     .filterBounds(aoi)
                    )
                     
    # extract appropriate bands and mask away clouds, shadows, etc
    srCollection = srCollection.map(lambda img:
                        maskSRcollection(img, sensor, maskThese)
                    )

    return srCollection

    
# For each year, build a medoid mosaic and add it to the collection
def getMedianMosaic(startDate, endDate, aoi, sensors=['LT05', 'LE07', 'LC08', 'LC09'], maskThese=['cloud','shadow','snow']):
    
    # Build a combined collection of SR imagery from each sensor for the dates and aoi specified
    c = ee.ImageCollection([ee.Image([0,0,0,0,0,0]).mask(ee.Image(0))])
    for sensor in sensors:
        c = c.merge(getSRcollection(startDate, endDate, aoi, sensor, maskThese))
    
    # Mosaic by selecting the unmasked pixel with the smallest distance across
    # all bands from the median pixel of each band. Simply choosing the raw  
    # median will give band values for a pixel selected from multiple images
    median = ee.ImageCollection(c).median()
    diff = c.map(lambda img:
                        ee.Image(img).subtract(median).abs()
                        .reduce('sum').addBands(img)
                ).reduce(ee.Reducer.min(4))
    img = diff.select([1,2,3,4,5,6], ['B1', 'B2', 'B3', 'B4', 'B5', 'B7'])
    
    # Rescale SR data to uint16 range
    img = img.multiply(0.0000275).add(-0.2).multiply(10000).toUint16()
    
    # Set date to a nominal value halfway between the start and end dates
    midDate = ee.Date((ee.Date(startDate).millis().add(ee.Date(endDate).millis())).divide(2) )
    return img.set('system:time_start', midDate)
    
    
# Create a collection of median mosaics for a list of start and end dates
def medianMosiacCollection(dates, aoi, sensors=['LT05', 'LE07', 'LC08', 'LC09'], maskThese=['cloud','shadow','snow']):
    return ee.ImageCollection(
                ee.List(dates).map(lambda d: 
                    getMedianMosaic(ee.List(d).get(0), ee.List(d).get(1), aoi, sensors, maskThese)
                )
            )
            
# Get bounds for our Area of Intrest in Brazil
brazil = ee.FeatureCollection("USDOS/LSIB/2017").filter(ee.Filter.eq('COUNTRY_NA', 'Brazil'));
BrazilBox = ee.Geometry.BBox(-75, -34, -34, 6)

# Export a median image
task = ee.batch.Export.image.toCloudStorage(
                            image=getMedianMosaic('2009-01-01', '20010-12-31', brazil).select('B5','B4','B3'),
                            description="MedianBrazil_L5789_2009-2010",
                            bucket="earthshot-trial-eeoutputs",
                            scale=1000,
                            maxPixels=2e8,
                            region=BrazilBox
                            )
task.start()