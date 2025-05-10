from flask import Flask, redirect, render_template, request, url_for
import ee
import folium
from datetime import datetime, timedelta
import sqlite3

service_account = 'afflux@afflux-1923.iam.gserviceaccount.com'
credentials = ee.ServiceAccountCredentials(service_account, 'afflux.json')
ee.Initialize(credentials)

app = Flask(__name__)
db_name = 'login.db'

app.jinja_env.globals['url_for'] = url_for

def earth_engine(loc, a_date_srt, e_date_srt):
    # Authenticate to ee Engine using a service account
    
    
    admin2 = ee.FeatureCollection("FAO/GAUL_SIMPLIFIED_500m/2015/level2")
    s2 = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")

    def pre_date(fetched_date_str):
        # Convert the fetched date to a datetime object
        fetched_date = datetime.strptime(fetched_date_str, "%Y-%m-%d").date()

        # Example timedelta object
        time_difference = timedelta(days=29)

        # Subtract the timedelta from the fetched date
        result_date = fetched_date - time_difference
        result = (str(result_date))
        return result

    before_start= pre_date(a_date_srt)
    before_end= pre_date(e_date_srt)

    # Now set the same parameters for AFTER the flood.
    after_start= a_date_srt
    after_end= e_date_srt

    polarization = "VH" 
    pass_direction = "DESCENDING" 
    difference_threshold = 1.25 
    # MOdification done for taking the results of DUBAI
    #dubai = admin2.filter(ee.Filter.eq('ADM0_NAME', 'United Arab Emirates'))
    kerala = admin2.filter(ee.Filter.eq('ADM2_NAME', loc))
    geometry = kerala.geometry()
    
    #rename selected geometry feature 
    aoi = ee.FeatureCollection(geometry)

    #location for folium
    centroid = geometry.centroid()
    coordinates = centroid.coordinates()
    long = ee.Number(coordinates.get(0)).format('%.2f')
    lat = ee.Number(coordinates.get(1)).format('%.2f')
    longitude = float(long.getInfo())
    latitude = float(lat.getInfo() )
    collection = ee.ImageCollection('COPERNICUS/S1_GRD')\
        .filter(ee.Filter.eq('instrumentMode','IW'))\
        .filter(ee.Filter.listContains('transmitterReceiverPolarisation', polarization))\
        .filter(ee.Filter.eq('orbitProperties_pass',pass_direction)) \
        .filter(ee.Filter.eq('resolution_meters',10))\
        .filter(ee.Filter.bounds(aoi))\
        .select(polarization)
    
    before_collection = collection.filterDate(before_start, before_end)
    after_collection = collection.filterDate(after_start,after_end)

    def dates(imgcol):
        range = imgcol.reduceColumns(ee.Reducer.minMax(), ["system:time_start"])
        printed = ee.String('from ')\
            .cat(ee.Date(range.get('min')).format('YYYY-MM-dd'))\
            .cat(' to ')\
            .cat(ee.Date(range.get('max')).format('YYYY-MM-dd'))
        return printed
    
    before_count = before_collection.size()
    print(ee.String('Tiles selected: Before Flood ').cat('(').cat(before_count).cat(')'),dates(before_collection), before_collection)
    
    # print dates of after images to console
    after_count = before_collection.size()
    print(ee.String('Tiles selected: After Flood ').cat('(').cat(after_count).cat(')'),dates(after_collection), after_collection)

    # Create a mosaic of selected tiles and clip to study area
    before = before_collection.mosaic().clip(aoi)
    after = after_collection.mosaic().clip(aoi)

    # Apply reduce the radar speckle by smoothing  
    smoothing_radius = 50
    before_filtered = before.focal_mean(smoothing_radius, 'circle', 'meters')
    after_filtered = after.focal_mean(smoothing_radius, 'circle', 'meters')

    # //------------------------------- FLOOD EXTENT CALCULATION -------------------------------//

    # Calculate the difference between the before and after images
    difference = after_filtered.divide(before_filtered)

    # Apply the predefined difference-threshold and create the flood extent mask 
    threshold = difference_threshold
    difference_binary = difference.gt(threshold)

    # Refine flood result using additional datasets
        
    # Include JRC layer on surface water seasonality to mask flood pixels from areas
    # of "permanent" water (where there is water > 10 months of the year)
    swater = ee.Image('JRC/GSW1_0/GlobalSurfaceWater').select('seasonality')
    swater_mask = swater.gte(10).updateMask(swater.gte(10))
        
    # Flooded layer where perennial water bodies (water > 10 mo/yr) is assigned a 0 value
    flooded_mask = difference_binary.where(swater_mask, 0)
    # final flooded area without pixels in perennial waterbodies
    flooded = flooded_mask.updateMask(flooded_mask)
        
    # Compute connectivity of pixels to eliminate those connected to 8 or fewer neighbours
    # This operation reduces noise of the flood extent product 
    connections = flooded.connectedPixelCount()    
    flooded = flooded.updateMask(connections.gte(8))
        
    # Mask out areas with more than 5 percent slope using a Digital Elevation Model 
    DEM = ee.Image('WWF/HydroSHEDS/03VFDEM')
    terrain = ee.Algorithms.Terrain(DEM)
    slope = terrain.select('slope')
    flooded = flooded.updateMask(slope.lt(5))

    # Calculate flood extent area
    # Create a raster layer containing the area information of each pixel 
    flood_pixelarea = flooded.select(polarization).multiply(ee.Image.pixelArea())

    # Sum the areas of flooded pixels
    # default is set to 'bestEffort: true' in order to reduce computation time, for a more 
    # accurate result set bestEffort to false and increase 'maxPixels'. 
    flood_stats = flood_pixelarea.reduceRegion(
        reducer=ee.Reducer.sum(),  # Example reducer: Sum the pixel values
        geometry=aoi,
        scale=10,
        bestEffort=True
    )

    # Convert the flood extent to hectares (area calculations are originally given in meters)  
    flood_area_ha = flood_stats.getNumber(polarization).divide(10000).round()

    population_count = ee.Image('JRC/GHSL/P2016/POP_GPW_GLOBE_V1/2015').clip(aoi)

    # Calculate the amount of exposed population get GHSL projection
    GHSLprojection = population_count.projection()

    # Reproject flood layer to GHSL scale
    flooded_res1 = flooded.reproject(crs=GHSLprojection)

    # Create a raster showing exposed population only using the resampled flood layer
    population_exposed = population_count.updateMask(flooded_res1).updateMask(population_count)

    # Sum pixel values of exposed population raster 
    stats = population_exposed.reduceRegion(
        reducer=ee.Reducer.sum(),
        geometry=aoi,
        scale=250,
        maxPixels=1e9
    )

    # get number of exposed people as an integer
    number_pp_exposed = stats.getNumber('population_count').round()

    '''//----------------------------- Affected agricultural land ----------------------------//
    // using MODIS Land Cover Type Yearly Global 500m
    // filter image collection by the most up-to-date MODIS Land Cover product '''

    LC = ee.ImageCollection('MODIS/006/MCD12Q1') \
        .filterDate('2014-01-01', after_end) \
        .sort('system:index', False) \
        .select("LC_Type1") \
        .first() \
        .clip(aoi)

    # Extract only cropland pixels using the classes cropland (>60%) and Cropland/Natural 
    # Vegetation Mosaics: mosaics of small-scale cultivation 40-60% incl. natural vegetation
    cropmask = LC.eq(12).Or(LC.eq(14))
    cropland = LC.updateMask(cropmask)
    
    # get MODIS projection
    MODISprojection = LC.projection()

    # Reproject flood layer to MODIS scale
    flooded_res = flooded.reproject(crs= MODISprojection)

    # Calculate affected cropland using the resampled flood layer
    cropland_affected = flooded_res.updateMask(cropland)

    # get pixel area of affected cropland layer
    crop_pixelarea = cropland_affected.multiply(ee.Image.pixelArea())  # calculate the area of each pixel 

    # sum pixels of affected cropland layer
    crop_stats = crop_pixelarea.reduceRegion(
        reducer=ee.Reducer.sum(),
        geometry=aoi,
        scale=500,
        maxPixels=1e9
    )

    # convert area to hectares
    crop_area_ha = crop_stats.getNumber(polarization).divide(10000).round()

    '''//-------------------------------- Affected urban area ------------------------------//

    // Using the same MODIS Land Cover Product 
    // Filter urban areas'''
    urbanmask = LC.eq(13)
    urban = LC.updateMask(urbanmask)

    # Calculate affected urban areas using the resampled flood layer
    urban_affected = urban.mask(flooded_res).updateMask(urban)

    # get pixel area of affected urban layer
    urban_pixelarea = urban_affected.multiply(ee.Image.pixelArea())  # calculate the area of each pixel 

    # sum pixels of affected cropland layer
    urban_stats = urban_pixelarea.reduceRegion(
        reducer=ee.Reducer.sum(),
        geometry=aoi,
        scale=500,
        bestEffort=True
    )

    # convert area to hectares
    urban_area_ha = urban_stats.getNumber('LC_Type1').divide(10000).round()

    '''Map view'''
    center = [latitude,longitude]
    m = folium.Map(location=center, zoom_start=9)
    # Displaying image
    # Add the after_filtered layer
    folium.TileLayer(
        tiles=after_filtered.getMapId({'min': -25, 'max': 0})['tile_fetcher'].url_format,
        attr='Google Earth Engine',
        overlay=True,
        name='after_filtered'
    ).add_to(m)

    # Add the before_filtered layer
    folium.TileLayer(
        tiles=before_filtered.getMapId({'min': -25, 'max': 0})['tile_fetcher'].url_format,
        attr='Google Earth Engine',
        overlay=True,
        name='before_filtered',
        show=False  # Checkbox turned off initially
    ).add_to(m)

    # Add the difference layer
    folium.TileLayer(
        tiles=difference.getMapId({'min': 0, 'max': 2})['tile_fetcher'].url_format,
        attr='Google Earth Engine',
        overlay=True,
        name='difference',
        show=False  # Checkbox turned off initially
    ).add_to(m)

    # Add the population_count layer
    folium.TileLayer(
        tiles=population_count.getMapId({'min': 0, 'max': 200,'palette':['060606','337663','337663','ffffff']})['tile_fetcher'].url_format,
        attr='Google Earth Engine',
        overlay=True,
        name='population Density',
        show=False  # Checkbox turned off initially
    ).add_to(m)

    # Add the flooded layer
    folium.TileLayer(
        tiles=flooded.getMapId({  'min': 0.0, 'max': 100.0,'palette': ['0000ff']})['tile_fetcher'].url_format,
        attr='Google Earth Engine',
        overlay=True,
        name='flooded'
    ).add_to(m)

    # Add the population_exposed layer
    folium.TileLayer(
        tiles=population_exposed.getMapId({  'min': 0, 'max': 200.0, 'palette': ['yellow', 'orange', 'red']})['tile_fetcher'].url_format,
        attr='Google Earth Engine',
        overlay=True,
        name='population exposed'
    ).add_to(m)

    # Add the MODIS Land Cover layer
    folium.TileLayer(
        tiles=LC.getMapId({  'min': 1.0, 'max': 17.0, 'palette': ['05450a', '086a10', '54a708', '78d203', '009900', 'c6b044', 'dcd159',
    'dade48', 'fbff13', 'b6ff05', '27ff87', 'c24f44', 'a5a5a5', 'ff6d4c',
    '69fff8', 'f9ffa4', '1c0dff']})['tile_fetcher'].url_format,
        attr='Google Earth Engine',
        overlay=True,
        name='Land Cover',
        show=False  # Checkbox turned off initially
    ).add_to(m)

    # Add the Cropland layer
    folium.TileLayer(
        tiles=cropland.getMapId({  'min': 0, 'max': 14.0, 'palette': ['30b21c']})['tile_fetcher'].url_format,
        attr='Google Earth Engine',
        overlay=True,
        name='cropland',
        show=False  # Checkbox turned off initially
    ).add_to(m)

    # Add the Urban layer
    folium.TileLayer(
        tiles=urban.getMapId({  'min': 0, 'max': 13.0, 'palette': ['grey']})['tile_fetcher'].url_format,
        attr='Google Earth Engine',
        overlay=True,
        name='urban',
        show=False  # Checkbox turned off initially
    ).add_to(m)

    # Add the Affected urban layer
    folium.TileLayer(
        tiles=urban_affected.getMapId({  'min': 0, 'max': 13.0, 'palette': ['grey']})['tile_fetcher'].url_format,
        attr='Google Earth Engine',
        overlay=True,
        name='urban affected'
    ).add_to(m)
    

    # Add the layer control
    folium.LayerControl().add_to(m)
    


    after_c=dates(after_collection).getInfo()
    f_hect=flood_area_ha.getInfo()
    population=number_pp_exposed.getInfo()
    MODIS_date = ee.String(LC.get('system:index')).slice(0, 4)
    MODIS_d=MODIS_date.getInfo()
    crop_ha=crop_area_ha.getInfo()
    urban_a=urban_area_ha.getInfo()
    # Render the map


    return m,after_start,after_end,after_c,f_hect,population,MODIS_d,crop_ha,urban_a

# Loading page
@app.route('/')
def loading():
    return render_template('index.html')
# Home page
@app.route('/home')
def home():
    return render_template('home.html')
#login
@app.route('/signlog', methods=['GET','POST'])
def signlog():
    conn = sqlite3.connect(db_name)
    if request.method == 'POST':
        if 'signup' in request.form:
            # Get form data for sign-up
            name = request.form['name']
            email=request.form['email']
            password = request.form['password']

            # Perform sign-up logic (e.g., store user details in a database)
            # Connect to the database
            #conn = sqlite3.connect(db_name)
            cursor = conn.cursor()

            # Execute queries to store or retrieve data from the database
            # For example, to insert user data into a table:
            cursor.execute('INSERT INTO login (name, email, password) VALUES (?, ?, ?)', (name, email, password))
            # Commit the changes and close the connection
            conn.commit()
            conn.close()
            return render_template('signlog.html', show_login=True)  # Redirect to the login page after sign-up
        elif 'login' in request.form:
            # Get form data for login
            username = request.form['email']
            password = request.form['password']
            cursor = conn.cursor()

            # Execute a query to retrieve the user details from the table
            cursor.execute('SELECT * FROM login WHERE email=? AND password=?', (username, password))
            user = cursor.fetchone()

            conn.close()
            # Perform login logic (e.g., verify user credentials)
            if user:
                return render_template('home.html')  # Return a success message if login is successful
            else:
                error = "Incorrect password & usernmae. Please try again."
                return render_template('signlog.html', show_login=True,error=error)
    else:   
        return render_template('signlog.html')

@app.route('/process_form', methods=["POST", "GET"])
def process_form():
    loc = request.form['place']
    a_date_srt = request.form['startdate']
    e_date_srt = request.form['enddate']

    m,ab,ae,af_c,f_hect,population,MODIS_d,crop_ha,urban_a = earth_engine(loc,
                         a_date_srt, e_date_srt)

    # Do something with the form data
    # For example, you can pass the data to a template
    
    
    return render_template('split.html', folium_map=m._repr_html_(),after_start=ab,
                           after_end=ae,
                           after_c=af_c,
                           f_hect=f_hect,
                           population=population,
                           MODIS_d=MODIS_d,
                           crop_ha=crop_ha,
                           urban_a=urban_a)
                                                                                                                                       
@app.route('/split')
def split():
    return render_template('split.html')

if __name__ == '__main__':
    app.run(debug=True)
