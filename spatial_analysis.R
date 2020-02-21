# Load the sf package
library(sf)

# Read in the trees shapefile
trees <- st_read("trees.shp")

# Read in the neighborhood shapefile
neighborhoods <- st_read("neighborhoods.shp")

# Read in the parks shapefile
parks <- st_read("parks.shp")

# View the first few trees
head(trees)


# Load the raster package
library(raster)

# Read in the tree canopy single-band raster
canopy <- raster("canopy.tif")

# Read in the manhattan Landsat image multi-band raster
manhattan <- brick("manhattan.tif")

# Get the class for the new objects
class(canopy)
class(manhattan)

# Identify how many layers each object has
nlayers(canopy)
nlayers(manhattan)

# Load the sf package
library(sf)

# ... and the dplyr package
library(dplyr)

# Read in the trees shapefile
trees <- st_read("trees.shp")

# Use filter() to limit to honey locust trees
honeylocust <- trees %>% filter(species == "honeylocust")

# Count the number of rows
nrow(honeylocust)

# Limit to tree_id and boroname variables
honeylocust_lim <- honeylocust %>% select(tree_id, boroname) 

# Use head() to look at the first few records
head(honeylocust_lim)

# Create a standard, non-spatial data frame with one column
df <- data.frame(a = 1:3)

# Add a list column to your data frame
df$b <- list(1:4, 1:5, 1:10)

# Look at your data frame with head
head(df)

# Convert your data frame to a tibble and print on console
as_tibble(df)

# Pull out the third observation from both columns individually
df$a[3]
df$b[3]

# Read in the parks shapefile
parks <- st_read("parks.shp")

# Compute the areas of the parks
areas <- st_area(parks)

# Create a quick histogram of the areas using hist
hist(areas, xlim = c(0, 200000), breaks = 1000)

# Filter to parks greater than 30000 (square meters)
big_parks <- parks %>% filter(unclass(areas) > 30000)

# Plot just the geometry of big_parks
plot(st_geometry(big_parks))

# Plot the parks object using all defaults
plot(parks)

# Plot just the acres attribute of the parks data
plot(parks["acres"])

# Create a new object of just the parks geometry
parks_geo <- st_geometry(parks)

# Plot the geometry of the parks data
plot(parks_geo)

# Load the raster package
library(raster)

# Read in the rasters
canopy <- raster("canopy.tif")
manhattan <- brick("manhattan.tif")

# Get the extent of the canopy object
extent(canopy)

# Get the CRS of the manhattan object
crs(manhattan)

# Determine the number of grid cells in both raster objects
ncell(manhattan)
ncell(canopy)

# Check if the data is in memory
inMemory(canopy)

# Use head() to peek at the first few records
head(canopy)

# Use getValues() to read the values into a vector
vals <- getValues(canopy)

# Use hist() to create a histogram of the values
hist(vals)

# Plot the canopy raster (single raster)
plot(canopy)

# Plot the manhattan raster (as a single image for each layer)
plot(manhattan)

# Plot the manhattan raster as an image
plotRGB(manhattan)