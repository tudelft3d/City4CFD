# Importing the required libraries
import osmnx as ox
import numpy as np
import os
import geopandas as gpd
from shapely.geometry import Point
from geopy.geocoders import Nominatim
from pyproj import CRS
# 
# USER INPUT PARAMETERS
#
Hmax = 230                            # Maximum height of the building the region of interest
rbuildings = 1200                     # Radius for building polygons
rpolygons = 4000                      # Radius for non-building polygons (only active when Hmax < 0)
outdir = "data"                       # Define the output directory where polygons are downloaded 
#
# Initialize geocoder
#
geolocator = Nominatim(user_agent="osm_downloader")
#
# Define the file read function
#
def parse_city_line(line):
    parts = line.strip().split(",")
    parts = [p.strip() for p in parts]

    city = parts[0]
    country = parts[1] if len(parts) > 1 else ""
    lat = parts[2] if len(parts) > 2 and parts[2] else None
    lon = parts[3] if len(parts) > 3 and parts[3] else None
    crs_code = parts[4] if len(parts) > 4 and parts[4] else "EPSG:4326"

    # Convert lat/lon to float if available
    try:
        lat = float(lat) if lat else None
        lon = float(lon) if lon else None
    except ValueError:
        lat, lon = None, None

    # Convert CRS
    try:
        crs = CRS.from_user_input(crs_code)
    except:
        print(f"Invalid CRS for {city}, defaulting to EPSG:4326")
        crs = CRS.from_epsg(4326)

    return city, country, lat, lon, crs
#
# Ensure data directory exists
#
os.makedirs(outdir, exist_ok=True)
#
# Load city names, coordinates, and CRS
#
cities = []
with open("cities.txt", "r") as file:
    for line in file:
        city, country, lat, lon, crs = parse_city_line(line)
        cities.append((f"{city}, {country}", lat, lon, crs))
#
# Radius in meters for downloading the OSM data
# First element is for buildings and second element is for other features
#
# Setup the radius based on COST72 guidelines
if (Hmax > 0):
    radii = [rbuildings, 15*Hmax + 2*rbuildings]
else:
    radii = [rbuildings, rpolygons]
#
# Define OSM tags
#
tags = {
    "buildings": {
         "building": True,
         "height": True,
         "building:levels": True
     },
    "vegetation": {
        "landuse": ["forest", "grass", "meadow", "orchard"],
        "leisure": ["park", "nature_reserve", "garden", "dog_park"],
        "natural": ["grassland", "wood"]
    },
    "water": {
        "natural": ["water", "sea"],
        "waterway": True,
    },
    "ocean": {
        "natural": ["coastline"]
    }
}
#
# Process each city
#
for city, lat, lon, crs in cities:
    print(f"Processing: {city} (CRS: {crs.to_string()})")

    # Use geocoder if coordinates are missing
    if lat is None or lon is None:
        location = geolocator.geocode(city)
        if not location:
            print(f"Could not find coordinates for {city}. Skipping...")
            continue
        lat, lon = location.latitude, location.longitude

    print(f"Coordinates for {city}: ({lat}, {lon})")

    # Create circular clipping area in specified CRS
    if(crs.to_string() != "EPSG:4326"):
        city_center = gpd.GeoDataFrame(geometry=[Point(lon, lat)], crs="EPSG:4326").to_crs(crs)
    else:
        city_center = gpd.GeoDataFrame(geometry=[Point(lon, lat)], crs=crs)
    buffer_distance = radii[1]

    clipping_area = city_center.buffer(buffer_distance)

    print(f"City4CFD: {city} | {city_center})")

    # Download and process OSM features
    for category, osm_tags in tags.items():
        radius_meters = radii[0] if category == "buildings" else radii[1]

        print(f"\nDownloading {category} within {radius_meters}m of {city}...")

        try:
            gdf = ox.features_from_point((lat, lon), tags=osm_tags, dist=radius_meters)
            if gdf.empty:
                print(f"No {category} data found for {city}.")
                continue

            gdf = gdf.to_crs(crs)
            clipped_gdf = gpd.clip(gdf, clipping_area)

            if clipped_gdf.empty:
                print(f"No {category} data remains after clipping for {city}.")
                continue

            filename = f"{outdir}/{city.replace(',', '').replace(' ', '_').lower()}_{category}.geojson"
            clipped_gdf.to_file(filename, driver="GeoJSON")
            print(f"Saved clipped {category} data to {filename}")

            # Calculate and print the extent with 10 unit offset in each direction
            bounds = clipped_gdf.total_bounds  # returns (xmin, ymin, xmax, ymax)
            offset_bounds = (bounds[0] - 10, bounds[1] - 10, bounds[2] + 10, bounds[3] + 10)
            print(f"Extended bounds for {category} (xmin-10 ymin-10 xmax+10 ymax+10): {bounds[0]-10} {bounds[1]-10} {bounds[2]+10} {bounds[3]+10}")

        except Exception as e:
            print(f"Error downloading {category} for {city}: {e}")
