library(sfarrow)
library(tmap)
library(ncdf4)
library(dplyr)

fscore <- function(preddepth, truedepth, threshold = 0.001) {
  
  preddepth[is.na(preddepth)] <- 0
  truedepth[is.na(truedepth)] <- 0
  
  pred <- preddepth * 0
  pred[preddepth >= threshold] <- 1
  obs <- truedepth * 0
  obs[truedepth >= threshold] <- 1
  
  true_pred  <- as.numeric(pred == 1 & obs == 1)
  overpred   <- as.numeric(pred == 1 & obs == 0)
  underpred  <- as.numeric(pred == 0 & obs == 1)
  
  tp_count <- sum(true_pred)
  op_count <- sum(overpred)
  up_count <- sum(underpred)
  
  flood_metric <- tp_count / (tp_count + op_count + up_count)
  
  return(flood_metric)
}
rsq <- function(preddepth, truedepth) {
  
  preddepth[is.na(preddepth)] <- 0
  truedepth[is.na(truedepth)] <- 0
  
  model <- lm(preddepth ~ truedepth)
  r2 <- summary(model)$r.squared
  
  return(r2)
}
rmse <- function(preddpeth, truedepth) {
  
  preddepth[is.na(preddepth)] <- 0
  truedepth[is.na(truedepth)] <- 0
  
  squared_errors <- (truedepth - preddepth)^2
  mean_sq_error  <- mean(squared_errors)
  rmse     <- sqrt(mean_sq_error)
  return(rmse)
  
}


setwd("D:/FloodA5/test/carlisle")
mesh18 <- sfarrow::st_read_parquet("carlisle_mesh18_standard.parquet")

pts <- mesh18 |> sf::st_drop_geometry() |> 
  sf::st_as_sf(coords = c("center_lon", "center_lat"), 
               crs = 4326)
ptsvect <- pts |> sf::st_transform(crs = 27700) |> terra::vect()


nc_file <- nc_open("carlisle_standard_res18.h5")

lisflood_dir <- "F:/OneDrive - University of Canterbury/LISFLOOD-FP_TestCases/Carlisle/carlisle_5m"

lisflood_files <- list.files(lisflood_dir, pattern = "*.wd$")

hrs <- 1:120
metrics <- {}

for(i in hrs){
  message(i)
  
  lis_data <- terra::rast(file.path(lisflood_dir, lisflood_files[hrs[i+1]]))
  names(lis_data) <- "lis_depth"
  lis_data[is.na(lis_data)] <- 0
  
  hrstr <- sprintf("frames/%.6d/", i+1)
  data <- data.frame(id = ncvar_get(nc_file, "mesh/cell_ids"), 
                     depth = ncvar_get(nc_file, paste0(hrstr, "water_depth")),
                     saturation = ncvar_get(nc_file, paste0(hrstr, "saturation")))
  
  
  joineddata <- dplyr::left_join(pts, data, by = c("cell_id" = "id"))
  
  x <- terra::extract(lis_data, ptsvect)
  joineddata$lis_depth <- x$lis_depth
  
  metrics <- rbind(metrics, "Hours" = i,
                   "F" = fscore(joineddata$depth, joineddata$lis_depth),
                   "R2" = rsq(joineddata$depth, joineddata$lis_depth),
                   "RMSE" = rmse(joineddata$depth, joineddata$lis_depth))
}









library(ggplot2)

ggplot(data = joineddata, mapping = aes(x = lis_depth, y = depth)) + 
  geom_point()+
  coord_equal()




 

