library(sfarrow)
library(tmap)
library(ncdf4)
library(dplyr)

setwd("D:/FloodA5/test/square")
mesh <- sfarrow::st_read_parquet("square_mesh18_standard.parquet")
nc_file <- nc_open("square_rainpoint_40hr.h5")


tmap_mode("plot")
tm_shape(mesh) +
  tm_borders(
  )

data_20hr <- data.frame(id = ncvar_get(nc_file, "mesh/cell_ids"), 
                        depth = ncvar_get(nc_file, "frames/000021/water_depth"),
                        saturation = ncvar_get(nc_file, "frames/000021/saturation"),
                        velocity = ncvar_get(nc_file, "frames/000021/velocity"),
                        volume = ncvar_get(nc_file, "frames/000021/volume"))

data_40hr <- data.frame(id = ncvar_get(nc_file, "mesh/cell_ids"), 
                   depth = ncvar_get(nc_file, "frames/000041/water_depth"),
                   saturation = ncvar_get(nc_file, "frames/000041/saturation"),
                   velocity = ncvar_get(nc_file, "frames/000041/velocity"),
                   volume = ncvar_get(nc_file, "frames/000041/volume"))

nc_close(nc_file)


names(data_20hr) <- paste0(names(data_20hr), "_20hr")
names(data_40hr) <- paste0(names(data_40hr), "_40hr")


mesh <- dplyr::left_join(mesh, data_20hr, by = c("cell_id" = "id_20hr"))
mesh <- dplyr::left_join(mesh, data_40hr, by = c("cell_id" = "id_40hr"))
mesh <- mesh |> dplyr::select(-neighbours)

sf::st_crs(mesh) <- "EPSG:4326"
sf::st_write(mesh, "square_mesh18_standard_output.gpkg")


tm_shape(mesh) + tm_borders("grey") +
  tm_shape(dplyr::filter(mesh, saturation_40hr == 1)) +
    tm_polygons(
      fill = "depth_40hr", 
      fill.scale = tm_scale_continuous(values = "-viridis"), # Use 'values' instead of 'palette'
      fill.legend = tm_legend(title = "vals")
    )




# Ensure it's in a projected CRS to calculate true area and perimeter
mesh_utm <- sf::st_transform(mesh, crs = 32630) |>
  dplyr::filter(saturation_40hr == 1) |>
  sf::st_union() |>
  sf::st_buffer(100) |> sf::st_buffer(-100) |>
  sf::st_as_sf() |> sf::st_set_geometry("x")

# get the Polsby-Popper score (circularity index)
mesh_utm |>
  dplyr::mutate(
    area = sf::st_area(mesh_utm),
    perimeter = sf::st_perimeter(mesh_utm),
    circularity = (4 * pi * area) / (perimeter^2)
  )




