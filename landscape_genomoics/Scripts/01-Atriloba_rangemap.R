## ------------------------------------------------
## 0. Packages (install them once manually if needed)
## ------------------------------------------------
# install.packages(c("rgbif", "sf", "dplyr", "lwgeom"))

library(rgbif)
library(sf)
library(dplyr)
library(lwgeom)

## ------------------------------------------------
## 0b. Paths + run configuration
## ------------------------------------------------
data_dir <- path.expand("~/Desktop/LandscapeGenomicsPipeline/")
out_dir <- file.path(data_dir, "outputs", "01-rangemap")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## Choose ONE:
use_gbif_snapshot <- TRUE

## If `use_gbif_snapshot == TRUE`, set your GBIF download key (snapshot)
## DOI: 10.15468/dl.kwmm6k
download_key <- "0042383-260226173443078"

## ------------------------------------------------
## 1. (Optional) GBIF download – skip if you already have a key
## ------------------------------------------------
## If you want to re-download from GBIF instead of using a fixed key,
## uncomment this block and comment out the hard-coded key section below.

# if (Sys.getenv("GBIF_USER") == "" ||
#     Sys.getenv("GBIF_PWD")  == "" ||
#     Sys.getenv("GBIF_EMAIL") == "") {
#   stop("Set GBIF_USER, GBIF_PWD, and GBIF_EMAIL in your .Renviron before running.")
# }
#
if (!use_gbif_snapshot) {
  sp <- name_backbone(name = "Asimina triloba", kingdom = "Plantae")
  species_key <- sp$speciesKey

  dw <- occ_download(
    pred("taxonKey", species_key),
    pred("hasCoordinate", TRUE),
    pred("country", "US"),
    pred_not(pred("basisOfRecord", "FOSSIL_SPECIMEN")),
    pred_not(pred("hasGeospatialIssue", TRUE)),
    pred_gte("year", 1950),   # >= 1950
    pred_lte("year", 2025),   # <= 2025 (optional, can omit)
    format = "SIMPLE_CSV"
  )

  occ_download_wait(dw$key)
  download_key <- dw$key
}

## ------------------------------------------------
## 2. Use an existing GBIF download (fixed snapshot)
## ------------------------------------------------

if (Sys.getenv("GBIF_USER") == "" ||
    Sys.getenv("GBIF_PWD")  == "" ||
    Sys.getenv("GBIF_EMAIL") == "") {
  stop("Set GBIF_USER, GBIF_PWD, and GBIF_EMAIL in your .Renviron before running.")
}

d <- occ_download_get(download_key, overwrite = TRUE) %>%
  occ_download_import()

# Clean the data
d_clean <- d |>
  # 1) drop obvious coordinate problems
  filter(
    !is.na(decimalLongitude),
    !is.na(decimalLatitude),
    decimalLongitude >= -100, decimalLongitude <= -70,   # rough East/Central US box
    decimalLatitude  >=  25,  decimalLatitude  <=  45
  ) |>
  # 2) remove flagged geospatial issues
  # drop rows with GEOSPATIAL issues (issues is a string; NA means no issues recorded)
  filter(
    is.na(issue) | !grepl("GEOSPATIAL", issue)
  ) |>
  # 3) remove very imprecise coordinates
  filter(
    is.na(coordinateUncertaintyInMeters) |
      coordinateUncertaintyInMeters <= 20000    # keep <= 20 km
  ) |>
  # 4) keep reasonable basisOfRecord
  filter(
    basisOfRecord %in% c("PRESERVED_SPECIMEN", "HUMAN_OBSERVATION", "OBSERVATION")
  ) |>
  # 5) drop old records if you want a “modern” range
  filter(
    is.na(year) | year >= 1950
  ) |>
  # 6) remove duplicated coordinates
  distinct(decimalLongitude, decimalLatitude, .keep_all = TRUE)

## ------------------------------------------------
## 3. Convert to sf and intersect with ecoregions
## ------------------------------------------------

occ_sf <- st_as_sf(
  d_clean,
  coords = c("decimalLongitude", "decimalLatitude"),
  crs = 4326,
  remove = FALSE
)

ecoregions <- st_read("~/Desktop/LandscapeGenomicsPipeline/InputData/us_eco_l3_state_boundaries/us_eco_l3_state_boundaries.shp",
                      quiet = TRUE) |>
  st_transform(4326) |>
  st_make_valid()

# Spatial join: keep only points inside ecoregions
sf_use_s2(FALSE)  # avoid s2 topology issues for this step
occ_sf_ec <- st_join(occ_sf, ecoregions, left = FALSE)

## ------------------------------------------------
## 4. Build a buffered hull range polygon
## ------------------------------------------------

# Union all points into one geometry
occ_union <- st_union(occ_sf_ec)

# Concave hull (nicer outline than convex hull); tune ratio as needed
asimina_hull <- st_concave_hull(occ_union, ratio = 0.15)

# Buffer the hull by ~50 km to generalize the range
asimina_range_50km <- st_transform(asimina_hull, 3857) |>
  st_buffer(dist = 50 * 1000) |>          # 50 km
  st_transform(4326)

# Constrain to ecoregions (intersection)
asimina_range_ecol <- st_intersection(ecoregions, asimina_range_50km) |>
  st_make_valid()

# after asimina_range_ecol is created
asimina_range_ecol <- st_make_valid(asimina_range_ecol)

# dissolve internal boundaries into a single polygon
asimina_range_simple <- st_union(asimina_range_ecol)

# (optionally cast back to MULTIPOLYGON for consistency)
asimina_range_simple <- st_cast(asimina_range_simple, "MULTIPOLYGON")

## ------------------------------------------------
## 5. Save range polygon as shapefile
## ------------------------------------------------
range_gpkg <- file.path(out_dir, "Asimina_triloba_range_ecol.gpkg")
st_write(
  asimina_range_simple,
  range_gpkg,
  layer = "range",
  delete_dsn = TRUE
)

## Optional: also write a shapefile (more fragile than GeoPackage)
range_shp <- file.path(out_dir, "Asimina_triloba_range_ecol.shp")
st_write(
  asimina_range_simple,
  range_shp,
  delete_dsn = TRUE
)
