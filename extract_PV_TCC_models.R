###########################LOAD RASTERS & EXTRACT PIXEL VALUES#######################
library(pacman)
library(rgee)
library(googledrive)
library(mapview)
library(sp)
library(raster)
#install.packages('terra')
library(terra)

ee$Initialize(project='ee-avatitsi')
####################################################################################

##load my dataset from GEE
#Composite Sentinel-2 (median)
composite_s2_path<-"projects/ee-avatitsi/assets/S2_2020_median" #from asset details: image ID (#S2_median)
composite_s2 <- ee$Image(composite_s2_path)
composite_s2_info <- composite_s2$getInfo()
print(composite_s2_info)
#vis_params_s2 <- list(min = 0.02, max = 0.3, bands = c("B8", "B4", "B3"))
#s2<- Map$addLayer(composite_s2,vis_params_s2, name='composite_s2')

#Composite indices
composite_VI_path<-"projects/ee-avatitsi/assets/VIs_composite2_2020" #from asset details: image ID (#VIs_composite)
composite_VI <- ee$Image(composite_VI_path)
composite_VI_info <- composite_VI$getInfo()
print(composite_VI_info)
#vis_params_VI <- list(min = -1, max = 1, bands = c("NDVI"))
#VI<- Map$addLayer(composite_VI,vis_params_VI, name='composite_s2')

#elevation_slope_aspect
elevation_slope_aspect_path<-"projects/ee-avatitsi/assets/elevation_slope_aspect" #from asset details: image ID
elevation_slope_aspect <- ee$Image(elevation_slope_aspect_path)
elevation_slope_aspect_info <- elevation_slope_aspect$getInfo()
print(elevation_slope_aspect_info)

#forest height
forest_h_path<-"projects/ee-avatitsi/assets/Forest_height_2019_NAFR_gr_laea" #from asset details: image ID
forest_h <- ee$Image(forest_h_path)
forest_h_info <- forest_h$getInfo()
print(forest_h_info)

#lai
lai <- composite_VI$expression(
  expression = '3.618 * evi - 0.0118',
  opt_map = list(evi = composite_VI$select('EVI'))  # Selecting the 'EVI' band
)
#vis_params_VI <- list(min = -1, max = 3)
#VI<- Map$addLayer(lai,vis_params_VI, name='lai')

#display multiple layers
#s2+VI

#maes
maes_path<-"projects/ee-avatitsi/assets/life_maes_reclass2"
maes<-ee$Image(maes_path)
maes_info<-maes$getInfo()
print(maes_info)

##reproject and resample
composite_s2_laea<-composite_s2$reproject(crs = 'EPSG:3035',scale=20)
composite_s2_laea_info <- composite_s2_laea$projection()$getInfo()

new_crs <- composite_s2_laea$projection() #$getInfo()
new_scale <- new_crs$nominalScale() #$getInfo()

composite_VI_laea<-composite_VI$reproject(crs = new_crs, scale = new_scale)
composite_VI_laea_info <- composite_VI_laea$projection()$getInfo()

elevation_slope_aspect_laea_resample<-elevation_slope_aspect$reproject(crs = new_crs, scale = new_scale)$resample()
elevation_slope_aspect_laea_resample_info<-elevation_slope_aspect_laea_resample$projection()$getInfo()

forest_h_laea_resample<-forest_h$reproject(crs = new_crs, scale = new_scale)$resample()
forest_h_laea_resample_info<-forest_h_laea_resample$projection()$getInfo()

lai_laea<-lai$reproject(crs = new_crs,scale = new_scale)
lai_laea_info <- lai_laea$projection()$getInfo()

maes_laea<-maes$reproject(crs = new_crs,scale = new_scale)
maes_laea_info <- maes_laea$projection()$getInfo()

##clip and write raster
geometry<-ee$Geometry$Polygon(
  list(
    list(
      c(22.56855837991539,37.59321695598985),
      c(23.84022586038414,37.59321695598985),
      c(23.84022586038414,38.28205157052546),
      c(22.56855837991539,38.28205157052546),
      c(22.56855837991539,37.59321695598985)
    )
  )
)
s2_clip<-maes_laea$clip(geometry) #composite_s2_laea, composite_VI_laea, elevation_slope_aspect_laea_resample, lai_laea, forest_h_laea_resample, maes_laea
s2<-ee_as_rast(image= s2_clip, maxPixels = 1e13)
##############################################################
#extract pixel values
library(foreign)
library(shapefiles)
library(sf)

#read shapefile (laea)
points30x30_laea_shp <- st_read("E:/new/katerina/PhD/FINAL_POINTS_TCC/POINTS/TCC_ALL_etrs_CORINE_LIFEMAP_30x30.shp")
head(points30x30_laea_shp)
print(nrow(points30x30_laea_shp))

##batch process
# Function to process a batch of polygons
process_batch <- function(batch_polygons, composite, scale) {
  PV <- ee_extract(
    x = composite,
    y = batch_polygons,
    scale = new_scale,
    fun = ee$Reducer$mean(),
    sf = TRUE
  )
  return(PV)
}

# Initialize an empty list to store results
all_results <- list()
batch_size <- 2000
# Process polygons in batches
num_batches <- ceiling(nrow(points30x30_laea_shp) / batch_size)
for (i in 1:num_batches) {
  start_idx <- (i - 1) * batch_size + 1
  end_idx <- min(i * batch_size, nrow(points30x30_laea_shp))
  batch_polygons <- points30x30_laea_shp[start_idx:end_idx, ]
  
  # Process the current batch
  batch_result <- process_batch(batch_polygons, forest_h_laea_resample, scale = 20)#composite_s2_laea,composite_VI_laea,composite_VI_laea elevation_slope_aspect_laea_resample, lai_laea, forest_h_laea_resample
  
  # Append the result to the list
  all_results[[i]] <- batch_result
}

# Combine all batch results into a single data frame
combined_results <- do.call(rbind, all_results)

# Define the output file path
output_file <- "E:/new/katerina/PhD/TCC/pixel_values_rgee/laea2020_20_corrected/PV2020_foresth_laea.csv"

# Check if the file already exists
if (file.exists(output_file)) {
  stop("The file '", output_file, "' already exists. Set a different filename or delete the existing file.")
} else {
  # Save the combined results to a CSV file
  write.csv(combined_results, output_file, na = "NA", row.names = FALSE)
}
####################################Model development##############################################
##Load dataset
pv_all<- read.csv("E:/new/katerina/PhD/TCC/pixel_values_rgee/laea2020_20_corrected/PV2020.csv",header=TRUE) #If a header row exists then, the header should be set TRUE else header should set to FALSE.
head(pv_all)
###############Load dataset as data.table######################
library(data.table)
pv_all_dt <-fread("E:/new/katerina/PhD/TCC/pixel_values_rgee/laea2020_20_corrected/PV2020.csv",header=TRUE )
summary(pv_all_dt)
head(pv_all_dt)
anyNA(pv_all_dt)
pv_all_sub_dt<-pv_all_dt[,c(4:30)]
#h2o.describe(pv_all_sub_dt)
pv_all_sub_dt$pososto <- as.numeric(pv_all_sub_dt$pososto)
###############h2o install, load, intitiate##########################
##initiate h2o
library(h2o)
h2o.init(nthreads = -1, max_mem_size="12G")

##split dataset
library (tidyr)
library(dplyr)
library(rsample)
library(recipes)

###################split dataset####################################
pv_all_sub_dt_h2o<-h2o.uploadFile("E:/new/katerina/PhD/TCC/pixel_values_rgee/laea2020_20_corrected/PV2020.csv")
pv<-as.h2o(pv_all_sub_dt)
seed <- 1234
splits <- h2o.splitFrame(pv_all_sub_dt_h2o, seed = seed, ratios = c(0.7))

train <- splits[[1]]
test <- splits[[2]]
#################################################################
# For k_fold strategy we need to provide fold column
train$fold <- h2o.kfold_column(train, nfolds = 10, seed = seed)

categorical_columns <- c("MAES_Name")  # Categorical columns for target encoding #"MAES_Name",
target_column <- "pososto"           # Dependent variable (numerical)

# Apply Target Encoding on the categorical variable
encoded_train <- h2o.targetencoder(
  x = categorical_columns,
  y = target_column,
  training_frame = train,
  blending = TRUE,           
  inflection_point = 15,     
  smoothing = 10,            
  seed=seed
)

# New target encoded train and test sets
transformed_train <- h2o.transform(encoded_train, train, as_training=TRUE)
transformed_train_sub<-transformed_train[ , c(1,6:30,32,33)]
transformed_test <- h2o.transform(encoded_train, test, noise=0)
transformed_test_sub<-transformed_test[ , c(1,6:30,32)]

##
#h2o.exportFile(transformed_test, path = "E:/new/katerina/PhD/TCC/uncertainty_maps/transformed_test.csv", force = TRUE)
#h2o.exportFile(transformed_train, path = "E:/new/katerina/PhD/TCC/uncertainty_maps/transformed_train.csv", force = TRUE)
#y <- target_column #"pososto"
#x <- setdiff(names(train[ , c(1:27)]), y) #c(1:61)]
##
nfolds <- 10
##############train models###############################
# Train & Cross-validate a GBM
features_with_te <- setdiff(names(transformed_train_sub), c(target_column, categorical_columns))

my_gbm <- h2o.gbm(x = features_with_te,
                  y = target_column,
                  training_frame = transformed_train_sub,
                  distribution = "gaussian",
                  ntrees = 90,
                  max_depth = 6,
                  min_rows = 10,
                  learn_rate = 0.1,
                  sample_rate = 0.7,
                  col_sample_rate=0.4,
                  nbins=20,
                  fold_column = "fold",
                  keep_cross_validation_predictions = TRUE,
                  seed = 1)

# Stacked results train gbm
MSE_gbm_t<-h2o.performance(my_gbm, newdata = transformed_train_sub)@metrics$MSE
RMSE_gbm_t<-h2o.performance(my_gbm, newdata = transformed_train_sub)@metrics$RMSE
R2_gbm_t<-h2o.performance(my_gbm, newdata = transformed_train_sub)@metrics$r2
MRD_gbm_t<-h2o.performance(my_gbm, newdata = transformed_train_sub)@metrics$mean_residual_deviance
RMSLE_gbm_t<-h2o.performance(my_gbm, newdata = transformed_train_sub)@metrics$rmsle
MAE_gbm_t<-h2o.performance(my_gbm, newdata = transformed_train_sub)@metrics$mae

results_gbm_t<-c(MSE_gbm_t,RMSE_gbm_t, MAE_gbm_t,RMSLE_gbm_t,MRD_gbm_t,R2_gbm_t)
columnnames_gbm_t<- c("MSE",	"RMSE",	"MAE",	"RMSLE","MRD", "R2")
RESULT_MAT_gbm_t<-matrix(results_gbm_t, ncol = 6)
colnames(RESULT_MAT_gbm_t)<-columnnames_gbm_t
RESULT_MAT_gbm_t

# Stacked results test gbm
MSE_gbm<-h2o.performance(my_gbm, newdata = transformed_test_sub)@metrics$MSE
RMSE_gbm<-h2o.performance(my_gbm, newdata = transformed_test_sub)@metrics$RMSE
R2_gbm<-h2o.performance(my_gbm, newdata = transformed_test_sub)@metrics$r2
MRD_gbm<-h2o.performance(my_gbm, newdata = transformed_test_sub)@metrics$mean_residual_deviance
RMSLE_gbm<-h2o.performance(my_gbm, newdata = transformed_test_sub)@metrics$rmsle
MAE_gbm<-h2o.performance(my_gbm, newdata = transformed_test_sub)@metrics$mae

results_gbm<-c(MSE_gbm,RMSE_gbm, MAE_gbm,RMSLE_gbm,MRD_gbm,R2_gbm)
columnnames_gbm<- c("MSE",	"RMSE",	"MAE",	"RMSLE","MRD", "R2")
RESULT_MAT_gbm<-matrix(results_gbm, ncol = 6)
colnames(RESULT_MAT_gbm)<-columnnames_gbm
RESULT_MAT_gbm

# Train & Cross-validate a RF
my_rf <- h2o.randomForest(x = features_with_te,
                          y = target_column,
                          training_frame = transformed_train_sub,
                          ntrees = 600,
                          fold_column = "fold",
                          categorical_encoding = 'enum',
                          keep_cross_validation_predictions = TRUE,
                          seed = 1)


# Stacked results train RF
MSE_rf_t<-h2o.performance(my_rf, newdata = transformed_train_sub)@metrics$MSE
RMSE_rf_t<-h2o.performance(my_rf, newdata = transformed_train_sub)@metrics$RMSE
MAE_rf_t<-h2o.performance(my_rf, newdata = transformed_train_sub)@metrics$mae
RMSLE_rf_t<-h2o.performance(my_rf, newdata = transformed_train_sub)@metrics$rmsle
MRD_rf_t<-h2o.performance(my_rf, newdata = transformed_train_sub)@metrics$mean_residual_deviance
R2_rf_t<-h2o.performance(my_rf, newdata = transformed_train_sub)@metrics$r2

results_rf_t<-c(MSE_rf_t,RMSE_rf_t, MAE_rf_t,RMSLE_rf_t,MRD_rf_t,R2_rf_t)
columnnames_rf_t<- c("MSE",	"RMSE",	"MAE",	"RMSLE","MRD", "R2")
RESULT_MAT_rf_t<-matrix(results_rf_t, ncol = 6)
colnames(RESULT_MAT_rf_t)<-columnnames_rf_t
RESULT_MAT_rf_t

# Stacked results test RF
MSE_rf<-h2o.performance(my_rf, newdata = transformed_test_sub)@metrics$MSE
RMSE_rf<-h2o.performance(my_rf, newdata = transformed_test_sub)@metrics$RMSE
R2_rf<-h2o.performance(my_rf, newdata = transformed_test_sub)@metrics$r2
MRD_rf<-h2o.performance(my_rf, newdata = transformed_test_sub)@metrics$mean_residual_deviance
RMSLE_rf<-h2o.performance(my_rf, newdata = transformed_test_sub)@metrics$rmsle
MAE_rf<-h2o.performance(my_rf, newdata = transformed_test_sub)@metrics$mae

results_rf<-c(MSE_rf,RMSE_rf, MAE_rf,RMSLE_rf,MRD_rf,R2_rf)
columnnames_rf<- c("MSE",	"RMSE",	"MAE",	"RMSLE","MRD", "R2")
RESULT_MAT_rf<-matrix(results_rf, ncol = 6)
colnames(RESULT_MAT_rf)<-columnnames_rf
RESULT_MAT_rf

# Train & Cross-validate a glm
my_glm <- h2o.glm(family = "gaussian",
                  x = features_with_te,
                  y = target_column,
                  training_frame = transformed_train_sub,
                  lambda = 0,
                  alpha=1,
                  fold_column = "fold",
                  keep_cross_validation_predictions = TRUE,
                  seed=1 )

# Stacked results train glm
MSE_glm_t<-h2o.performance(my_glm, newdata = transformed_train_sub)@metrics$MSE
RMSE_glm_t<-h2o.performance(my_glm, newdata = transformed_train_sub)@metrics$RMSE
R2_glm_t<-h2o.performance(my_glm, newdata = transformed_train_sub)@metrics$r2
MRD_glm_t<-h2o.performance(my_glm, newdata = transformed_train_sub)@metrics$mean_residual_deviance
RMSLE_glm_t<-h2o.performance(my_glm, newdata = transformed_train_sub)@metrics$rmsle
MAE_glm_t<-h2o.performance(my_glm, newdata = transformed_train_sub)@metrics$mae

results_glm_t<-c(MSE_glm_t,RMSE_glm_t, MAE_glm_t,RMSLE_glm_t,MRD_glm_t,R2_glm_t)
columnnames_glm_t<- c("MSE",	"RMSE",	"MAE",	"RMSLE","MRD", "R2")
RESULT_MAT_glm_t<-matrix(results_glm_t, ncol = 6)
colnames(RESULT_MAT_glm_t)<-columnnames_gbm_t
RESULT_MAT_glm_t

# Stacked results testglm
MSE_glm<-h2o.performance(my_glm, newdata = transformed_test_sub)@metrics$MSE
RMSE_glm<-h2o.performance(my_glm, newdata = transformed_test_sub)@metrics$RMSE
R2_glm<-h2o.performance(my_glm, newdata = transformed_test_sub)@metrics$r2
MRD_glm<-h2o.performance(my_glm, newdata = transformed_test_sub)@metrics$mean_residual_deviance
RMSLE_glm<-h2o.performance(my_glm, newdata = transformed_test_sub)@metrics$rmsle
MAE_glm<-h2o.performance(my_glm, newdata = transformed_test_sub)@metrics$mae

results_glm<-c(MSE_glm,RMSE_glm, MAE_glm,RMSLE_glm,MRD_glm,R2_glm)
columnnames_glm<- c("MSE",	"RMSE",	"MAE",	"RMSLE","MRD", "R2")
RESULT_MAT_glm<-matrix(results_glm, ncol = 6)
colnames(RESULT_MAT_glm)<-columnnames_glm
RESULT_MAT_glm


# Train & Cross-validate a deep learning

my_dl <- h2o.deeplearning(x =features_with_te,
                          y = target_column,
                          distribution = "gaussian",
                          hidden = c(200),
                          epochs = 500,
                          train_samples_per_iteration = -2,
                          reproducible = TRUE,
                          activation = "Tanh",
                          single_node_mode = FALSE,
                          #balance_classes = FALSE,#default
                          force_load_balance = FALSE,
                          seed = 1,
                          fold_column = "fold",
                          #nfolds = nfolds,
                          keep_cross_validation_predictions = TRUE,
                          score_training_samples = 0,
                          score_validation_samples = 0,
                          training_frame = transformed_train_sub,
                          stopping_rounds = 10)

# Stacked results train dl
MSE_dl_t<-h2o.performance(my_dl, newdata = transformed_train_sub)@metrics$MSE
RMSE_dl_t<-h2o.performance(my_dl, newdata = transformed_train_sub)@metrics$RMSE
MAE_dl_t<-h2o.performance(my_dl, newdata = transformed_train_sub)@metrics$mae
RMSLE_dl_t<-h2o.performance(my_dl, newdata = transformed_train_sub)@metrics$rmsle
MRD_dl_t<-h2o.performance(my_dl, newdata = transformed_train_sub)@metrics$mean_residual_deviance
R2_dl_t<-h2o.performance(my_dl, newdata = transformed_train_sub)@metrics$r2

results_dl_t<-c(MSE_dl_t,RMSE_dl_t, MAE_dl_t,RMSLE_dl_t,MRD_dl_t,R2_dl_t)
columnnames_dl_t<- c("MSE",	"RMSE",	"MAE",	"RMSLE","MRD", "R2")
RESULT_MAT_dl_t<-matrix(results_dl_t, ncol = 6)
colnames(RESULT_MAT_dl_t)<-columnnames_dl_t
RESULT_MAT_dl_t

# Stacked results test dl
MSE_dl<-h2o.performance(my_dl, newdata = transformed_test_sub)@metrics$MSE
RMSE_dl<-h2o.performance(my_dl, newdata = transformed_test_sub)@metrics$RMSE
R2_dl<-h2o.performance(my_dl, newdata = transformed_test_sub)@metrics$r2
MRD_dl<-h2o.performance(my_dl, newdata = transformed_test_sub)@metrics$mean_residual_deviance
RMSLE_dl<-h2o.performance(my_dl, newdata = transformed_test_sub)@metrics$rmsle
MAE_dl<-h2o.performance(my_dl, newdata = transformed_test_sub)@metrics$mae

results_dl<-c(MSE_dl,RMSE_dl, MAE_dl,RMSLE_dl,MRD_dl,R2_dl)
columnnames_dl<- c("MSE",	"RMSE",	"MAE",	"RMSLE","MRD", "R2")
RESULT_MAT_dl<-matrix(results_dl, ncol = 6)
colnames(RESULT_MAT_dl)<-columnnames_dl
RESULT_MAT_dl

# Train a stacked ensemble
ensemble <- h2o.stackedEnsemble(x = features_with_te,
                                y = target_column,
                                training_frame = transformed_train_sub,
                                seed = 123,
                                metalearner_algorithm="AUTO",
                                base_models = list(my_gbm, my_glm, my_rf, my_dl)) #,my_glm

# Eval ensemble performance on a test set
perf <- h2o.performance(ensemble, newdata = transformed_test_sub)

# Stacked results train ensemble
MSE_ens_t<-h2o.performance(ensemble, newdata = transformed_train_sub)@metrics$MSE
RMSE_ens_t<-h2o.performance(ensemble, newdata = transformed_train_sub)@metrics$RMSE
MAE_ens_t<-h2o.performance(ensemble, newdata = transformed_train_sub)@metrics$mae
RMSLE_ens_t<-h2o.performance(ensemble, newdata = transformed_train_sub)@metrics$rmsle
MRD_ens_t<-h2o.performance(ensemble, newdata = transformed_train_sub)@metrics$mean_residual_deviance
R2_ens_t<-h2o.performance(ensemble, newdata = transformed_train_sub)@metrics$r2

results_ens_t<-c(MSE_ens_t,RMSE_ens_t, MAE_ens_t,RMSLE_ens_t,MRD_ens_t,R2_ens_t)
columnnames_ens_t<- c("MSE",	"RMSE",	"MAE",	"RMSLE","MRD", "R2")
RESULT_MAT_ens_t<-matrix(results_ens_t, ncol = 6)
colnames(RESULT_MAT_ens_t)<-columnnames_ens_t
RESULT_MAT_ens_t

# Stacked results test ensemble
MSE_ens<-h2o.performance(ensemble, newdata = transformed_test_sub)@metrics$MSE
RMSE_ens<-h2o.performance(ensemble, newdata = transformed_test_sub)@metrics$RMSE
R2_ens<-h2o.performance(ensemble, newdata = transformed_test_sub)@metrics$r2
MRD_ens<-h2o.performance(ensemble, newdata = transformed_test_sub)@metrics$mean_residual_deviance
RMSLE_ens<-h2o.performance(ensemble, newdata = transformed_test_sub)@metrics$rmsle
MAE_ens<-h2o.performance(ensemble, newdata = transformed_test_sub)@metrics$mae

results_ens<-c(MSE_ens,RMSE_ens, MAE_ens,RMSLE_ens,MRD_ens,R2_ens)
columnnames_ens<- c("MSE",	"RMSE",	"MAE",	"RMSLE","MRD", "R2")
RESULT_MAT_ens<-matrix(results_ens, ncol = 6)
colnames(RESULT_MAT_ens)<-columnnames_ens
RESULT_MAT_ens
####################explain model#############
exm_gbm <- h2o.explain(my_gbm, transformed_test_sub)
exm_rf <- h2o.explain(my_rf, transformed_test_sub)
exm_glm <- h2o.explain(my_glm, transformed_test_sub)
exm_dl <- h2o.explain(my_dl, transformed_test_sub)

###############################################

# save the model
setwd('E:/new/katerina/PhD/TCC/models_tc/model_tc_04022026')
getwd()
model_path <- h2o.saveModel(object = ensemble, path = getwd(), force = TRUE)
print(model_path)

# load the model
saved_model <- h2o.loadModel(model_path)

# download the model built above to your local machine
my_local_model <- h2o.download_model(ensemble, path = "E:/new/katerina/PhD/TCC/models_tc/model_tc_04022026")

# upload the model that you just downloded above to the H2O cluster
uploaded_model <- h2o.upload_model("E:/new/katerina/PhD/TCC/models_tc/model_tc_04022026/StackedEnsemble_model_R_1770194253674_6")
#################################################

##STACK MODELS FOR FOREST SPECIES
pv_all<- read.csv("E:/new/katerina/PhD/TCC/pixel_values_rgee/laea2020_20_corrected/PV2020_tempCon.csv",header=TRUE)
head(pv_all)

library(data.table)
#load data using fread (data.table)
pv_all_dt <-fread("E:/new/katerina/PhD/TCC/pixel_values_rgee/laea2020_20_corrected/PV2020_tempCon.csv",header=TRUE )
pv_all_sub_dt<-pv_all_dt[,c(4:29)]#,c(4:27)
pv_all_sub_dt$pososto <- as.numeric(pv_all_sub_dt$pososto)
#################################################################
#split dataset
set.seed(123)
split <- initial_split(pv_all_sub_dt, prop=0.7,strata = "pososto")
train <- training(split)
test <- testing(split) 

###
blueprint_tc <- recipe(pososto ~ ., data = pv_all_sub_dt)
summary(blueprint_tc)
head(blueprint_tc)

## Create training & test sets for h2o
train_h2o <- prep(blueprint_tc, training = train, retain = TRUE) %>%
  juice() %>%
  as.h2o()

test_h2o <- prep(blueprint_tc, training =  train) %>%
  bake(new_data = test) %>%
  as.h2o()

y <- "pososto"
x <- setdiff(names(train[ , c(1:26)]), y)#c(1:23,24)

#h2o.exportFile(test_h2o, path = "E:/new/katerina/PhD/TCC/uncertainty_maps/PV2020_tempCon_test.csv", force = TRUE)
#h2o.exportFile(train_h2o, path = "E:/new/katerina/PhD/TCC/uncertainty_maps/PV2020_tempCon.csv_train.csv", force = TRUE)

######################################################################
nfolds <- 10

# search criteria for grid search
search_criteria <- list(strategy = "RandomDiscrete", max_models = 50, seed = 1, max_runtime_secs = 900, stopping_tolerance = 0.001,stopping_rounds = 15)
#########################################################################
##hyperparameters
# GBM hyperparameters
gbm_params <- list(learn_rate = seq(0.1, 1, 0.1),
                   max_depth = seq(2, 10, 1),
                   sample_rate = seq(0.5, 1.0, 0.1),
                   min_rows=seq(1,10,1),
                   ntrees=c(20,50,100),
                   col_sample_rate = seq(0.1, 1.0, 0.1))

# RF hyperparameters
rf_params <- list(
  ntrees = c(500,550,600),
  max_depth = seq(5,15,5),
  min_rows = c(1,10),
  mtries = c(-1,2), 
  stopping_tolerance = 0.001
)

#DL hyperparameters
dl_params<-list(
  activation = c("Rectifier", "Maxout", "Tanh", "RectifierWithDropout", "MaxoutWithDropout", "TanhWithDropout"), 
  hidden = list(c(64), c(64, 32), c(5128, 64, 32)),
  epochs = c(10,31, 100, 500),
  rho = c(0.9, 0.95, 0.99),
  epsilon = c(1e-10, 1e-8, 1e-6),
  momentum_start = c(0, 0.5),
  momentum_stable = c(0.99, 0.5, 0),
  input_dropout_ratio = c(0, 0.1, 0.2),
  max_w2 = c(10, 100, 1000, 3.4e+38)
)

#GLM hyperparameters
hyper_params <- list( lambda = c(1, 0.5, 0.1, 0.01, 0.001, 0.0001, 0.00001, 0))
#######################################################################
# Train and validate a random grid of GBMs
gbm_grid <- h2o.grid("gbm", x = x, y = y,
                     grid_id = "gbm_grid",
                     training_frame = train_h2o,
                     hyper_params = gbm_params,
                     nfolds = nfolds,
                     keep_cross_validation_predictions = TRUE,
                     search_criteria = search_criteria)

gbm_gridperf <- h2o.getGrid(grid_id = "gbm_grid",
                            sort_by = "r2",
                            decreasing = TRUE)


best_gbm <- h2o.getModel(gbm_gridperf@model_ids[[1]])


best_gbm_perf <- h2o.performance(model = best_gbm,
                                 newdata = test_h2o)
h2o.r2(best_gbm_perf)


print(best_gbm@model[["model_summary"]])

################################################################
# Train and validate a random grid of RF
rf_grid<- h2o.grid("randomForest", x = x, y = y,
                   grid_id = "rf_grid",
                   training_frame = train_h2o,
                   hyper_params = rf_params,
                   nfolds = nfolds,
                   keep_cross_validation_predictions = TRUE,
                   search_criteria = search_criteria)

rf_gridperf <- h2o.getGrid(grid_id = "rf_grid",
                           sort_by = "r2",
                           decreasing = TRUE)


best_rf<- h2o.getModel(rf_gridperf@model_ids[[7]])


best_rf_perf<- h2o.performance(model = best_rf,
                               newdata = test_h2o)
h2o.r2(best_rf_perf)


print(best_rf@model[["model_summary"]])


h2o.varimp(best_rf)
################################################################################
# Train and validate a random grid of glm
glm_grid <- h2o.grid(algorithm = "glm",
                     family = "gaussian",
                     x = x,
                     y = y,
                     training_frame = train_h2o,
                     grid_id = "glm_grid",
                     #lambda = 0,
                     #alpha=1,
                     nfolds = nfolds,
                     keep_cross_validation_predictions = TRUE,
                     hyper_params = hyper_params,
                     search_criteria = search_criteria,
                     #compute_p_values = TRUE,
                     seed=1 )

glm_gridperf <- h2o.getGrid(grid_id = "glm_grid",
                            sort_by = "r2",
                            decreasing = TRUE)


best_glm <- h2o.getModel(glm_gridperf@model_ids[[2]])


best_glm_perf <- h2o.performance(model = best_glm,
                                 newdata = test_h2o)
h2o.r2(best_glm_perf)


print(best_glm@model[["model_summary"]])

####################################################################
# Train & Cross-validate a deep learning

dl_grid <- h2o.grid(algorithm = "deeplearning", 
                    x = x,
                    y = y,
                    #weights_column = weights,
                    grid_id = "dl_grid",
                    training_frame = train_h2o,
                    #validation_frame = valid,
                    nfolds = nfolds,
                    keep_cross_validation_predictions = TRUE,
                    fold_assignment = "Random",
                    hyper_params = dl_params,
                    search_criteria = search_criteria,
                    seed = 42
)

dl_gridperf <- h2o.getGrid(grid_id = "dl_grid",
                           sort_by = "r2",
                           decreasing = TRUE)


best_dl <- h2o.getModel(dl_gridperf@model_ids[[1]])
best_dl <- h2o.getModel("dl_grid_model_23")


best_dl_perf <- h2o.performance(model = best_dl,
                                newdata = test_h2o)
h2o.r2(best_dl_perf)

print(best_dl@model[["model_summary"]])
############################################################################
gbm_models <- h2o.getGrid(grid_id = gbm_grid@grid_id, sort_by = "r2", decreasing = TRUE)@model_ids
glm_models <- h2o.getGrid(grid_id = glm_grid@grid_id, sort_by = "r2", decreasing = TRUE)@model_ids
dl_models <- h2o.getGrid(grid_id = dl_grid@grid_id, sort_by = "r2", decreasing = TRUE)@model_ids
rf_models <- h2o.getGrid(grid_id = rf_grid@grid_id, sort_by = "r2", decreasing = TRUE)@model_ids

# Train a stacked ensemble
ensemble <- h2o.stackedEnsemble(x = x,
                                y = y,
                                training_frame = train_h2o,
                                seed = 123,
                                metalearner_algorithm="AUTO",
                                base_models = list(best_gbm,best_glm,best_dl,best_rf))

# Eval ensemble performance on a test set
perf <- h2o.performance(ensemble, newdata = test_h2o)
################################################################################
# Stacked results train ensemble
MSE_ens_t<-h2o.performance(ensemble, newdata = train_h2o)@metrics$MSE
RMSE_ens_t<-h2o.performance(ensemble, newdata = train_h2o)@metrics$RMSE
MAE_ens_t<-h2o.performance(ensemble, newdata = train_h2o)@metrics$mae
RMSLE_ens_t<-h2o.performance(ensemble, newdata = train_h2o)@metrics$rmsle
MRD_ens_t<-h2o.performance(ensemble, newdata = train_h2o)@metrics$mean_residual_deviance
R2_ens_t<-h2o.performance(ensemble, newdata = train_h2o)@metrics$r2

results_ens_t<-c(MSE_ens_t,RMSE_ens_t, MAE_ens_t,RMSLE_ens_t,MRD_ens_t,R2_ens_t)
columnnames_ens_t<- c("MSE",	"RMSE",	"MAE",	"RMSLE","MRD", "R2")
RESULT_MAT_ens_t<-matrix(results_ens_t, ncol = 6)
colnames(RESULT_MAT_ens_t)<-columnnames_ens_t
RESULT_MAT_ens_t

# Stacked results test ensemble
MSE_ens<-h2o.performance(ensemble, newdata = test_h2o)@metrics$MSE
RMSE_ens<-h2o.performance(ensemble, newdata = test_h2o)@metrics$RMSE
R2_ens<-h2o.performance(ensemble, newdata = test_h2o)@metrics$r2
MRD_ens<-h2o.performance(ensemble, newdata = test_h2o)@metrics$mean_residual_deviance
RMSLE_ens<-h2o.performance(ensemble, newdata = test_h2o)@metrics$rmsle
MAE_ens<-h2o.performance(ensemble, newdata = test_h2o)@metrics$mae

results_ens<-c(MSE_ens,RMSE_ens, MAE_ens,RMSLE_ens,MRD_ens,R2_ens)
columnnames_ens<- c("MSE",	"RMSE",	"MAE",	"RMSLE","MRD", "R2")
RESULT_MAT_ens<-matrix(results_ens, ncol = 6)
colnames(RESULT_MAT_ens)<-columnnames_ens
RESULT_MAT_ens
##############################################
# save the species-specific model
setwd('E:/new/katerina/PhD/TCC/models_tc/model_tc_04022026')
getwd()
model_path <- h2o.saveModel(object = ensemble, path = getwd(), force = TRUE)
print(model_path)

# load the model
saved_model <- h2o.loadModel(model_path)

# download the model built above to your local machine
my_local_model <- h2o.download_model(ensemble, path = "E:/new/katerina/PhD/TCC/models_tc/model_tc_04022026")

# upload the model that you just downloded above to the H2O cluster
#uploaded_model <- h2o.upload_model("E:/new/katerina/PhD/TCC/models_tc/model_tc_04022026/StackedEnsemble_model_R_1770194253674_7")#TempCon
#uploaded_model <- h2o.upload_model("E:/new/katerina/PhD/TCC/models_tc/model_tc_04022026/StackedEnsemble_model_R_1770198156105_1")#medCon
#uploaded_model <- h2o.upload_model("E:/new/katerina/PhD/TCC/models_tc/model_tc_04022026/StackedEnsemble_model_R_1770199800149_1")#medDec
#uploaded_model <- h2o.upload_model("E:/new/katerina/PhD/TCC/models_tc/model_tc_04022026/StackedEnsemble_model_R_1770201514776_1")#medScl
#uploaded_model <- h2o.upload_model("E:/new/katerina/PhD/TCC/models_tc/model_tc_04022026/StackedEnsemble_model_R_1770203236340_1")#mixed
#uploaded_model <- h2o.upload_model("E:/new/katerina/PhD/TCC/models_tc/model_tc_04022026/StackedEnsemble_model_R_1770204918231_1")#TempDec
###############################################
###uncertainty maps###########
##predictions-observations
pred_y = predict(ensemble, newdata=transformed_test_sub)
test_y=transformed_test_sub[,27]
pred_y <- as.vector(pred_y)
test_y <- as.vector(test_y)

##export csv
pred_test<-data.frame(Predictions = pred_y, Actual = test_y)  # Should return TRUE)
write.csv(pred_test, "E:/new/katerina/PhD/TCC/uncertainty_maps/pred_obs_test.csv", row.names = FALSE)
test_vector <- as.vector(test)
write.csv(test_vector, "E:/new/katerina/PhD/TCC/uncertainty_maps/test.csv", row.names = FALSE)

#polynomial model
data <- read.csv("E:/new/katerina/PhD/TCC/uncertainty_maps/rmse.csv")
data
sub_data<-data[,c(1,3)]

degree<-2
pol_model<-lm(sub_data$RMSE_ABS~ poly(sub_data$MEAN_PRED,degree,raw=TRUE))
summary(pol_model)
coefficients(pol_model)

#plot
plot(sub_data$MEAN_PRED,sub_data$RMSE_ABS, main = "RMSE uncertainty (relative)", xlab = "TCC (%)", ylab = "RMSE (%)", pch = 19, col = "dark green",cex = 1.1)
x_seq <- seq(min(sub_data$MEAN_PRED), max(sub_data$MEAN_PRED), length.out = 10)
#lines(x_seq,sub_data$RMSE_ABS, col = "red", lwd = 1.5, lty='dashed')
smooth_spline <- smooth.spline(sub_data$MEAN_PRED, sub_data$RMSE_ABS)
lines(smooth_spline, col = "red", lwd = 1.5, lty = 'dashed')

#uncertainty maps
compo_tcc<-raster("E:/new/katerina/PhD/TCC/OUTPUTS/predictions_corrected/tc_stack_all_clip_laea.tif")

calculate_rmse <- function(x) {
  0.03895028* x^2 - 5.51101488 * x + 214.48720665
}

rmse_map <- calc(compo_tcc, calculate_rmse)

calculate_rmse <- function(x) {
  y = -0.005009337*x^2 +0.556172771*x +8.702767499
}

writeRaster(rmse_map, "E:/new/katerina/PhD/TCC/TCC_text_plots/plots/rmse_uncertainty/rmse_absolute_map.tif", format="GTiff", overwrite=TRUE)
###################################################
# Load libraries
library(h2o)
library(ggplot2)
library(dplyr)
library(shapviz)

# Start H2O
h2o.init()
h2o.init(max_mem_size = "12G")
###################################################################################
#Load saved model

ensemble <- h2o.loadModel("E:/new/katerina/PhD/TCC/models_tc/model_tc_04022026/StackedEnsemble_model_R_1770194253674_6")
#uploaded_model <- h2o.upload_model("E:/new/katerina/PhD/TCC/models_tc/model_tc_04022026/StackedEnsemble_model_R_1770194253674_7")#TempCon
#uploaded_model <- h2o.upload_model("E:/new/katerina/PhD/TCC/models_tc/model_tc_04022026/StackedEnsemble_model_R_1770198156105_1")#medCon
#uploaded_model <- h2o.upload_model("E:/new/katerina/PhD/TCC/models_tc/model_tc_04022026/StackedEnsemble_model_R_1770199800149_1")#medDec
#uploaded_model <- h2o.upload_model("E:/new/katerina/PhD/TCC/models_tc/model_tc_04022026/StackedEnsemble_model_R_1770201514776_1")#medScl
#uploaded_model <- h2o.upload_model("E:/new/katerina/PhD/TCC/models_tc/model_tc_04022026/StackedEnsemble_model_R_1770203236340_1")#mixed
#uploaded_model <- h2o.upload_model("E:/new/katerina/PhD/TCC/models_tc/model_tc_04022026/StackedEnsemble_model_R_1770204918231_1")#TempDec


################################
#load test set
test<- read.csv("E:/new/katerina/PhD/TCC/uncertainty_maps/transformed_test.csv")
test_sub<-test[ , c(1,6:30,32)]
test_sub_df <- as.data.frame(test_sub)

# Identify response and features 
y <- "pososto" 
x <- setdiff(names(test_sub_df), y)

###############################
set.seed(142)
# create background dataset

train<- read.csv("E:/new/katerina/PhD/TCC/uncertainty_maps/transformed_train.csv")
train_sub<-train[ , c(1,6:30,32,33)]
train_sub_df <- as.data.frame(train_sub)
background_hf <-as.h2o(train_sub_df[sample(1:nrow(train_sub_df), 50), ])

#######################

test_global_df <- test_sub_df %>%
  mutate(canopy_bin = ntile(pososto, 6)) %>%  
  group_by(pososto) %>%
  slice_sample(n = 500) %>% 
  ungroup() %>%
  select(-pososto)
########################

test_global_hf <- as.h2o(test_global_df)


contrib_global <- h2o.predict_contributions(
  ensemble, 
  newdata = test_global_hf,
  background_frame = background_hf
)

contrib_g_df <- as.data.frame(contrib_global)
shap_matrix_g <- as.matrix(contrib_g_df[, setdiff(colnames(contrib_g_df), "BiasTerm")])
baseline_g <- as.numeric(contrib_g_df$BiasTerm[1])

shp_global <- shapviz(
  object = shap_matrix_g, 
  X = test_global_df, 
  baseline = baseline_g
)


sv_importance(
  shp_global, 
  kind = "beeswarm", 
  max_display = 5,  
  alpha = 0.3,       
  size = 1.5         
) + 
  scale_color_gradient(low = "blue", high = "red")+ # Μπλε σε Κόκκινο
  coord_cartesian(xlim = c(-20, 20))+theme_minimal()

p<-sv_importance(
  shp_global, kind = "bar", 
  max_display = 5 )+
  coord_cartesian(xlim = c(0, 8))

# Make bars thinner + set color
p$layers[[1]]$aes_params$fill <- "#2C7BB6"
p$layers[[1]]$aes_params$width <- 0.5

p +
  theme_minimal(base_size = 14) +   
  theme(
    #panel.grid.major.y = element_blank(), 
    #panel.grid.minor = element_blank(),
    axis.title.y = element_blank(),       
    plot.title = element_text(face = "bold")
  ) +
  labs(
    x = "Mean |SHAP value|"
  )+
  geom_text(
    aes(label = round(after_stat(x), 2)),
    color = "#2C7BB6",
    stat = "identity",
    hjust = -0.2,   
    size = 4
  ) +
  theme(
    axis.line.x = element_line(color = "lightgrey", linewidth = 0.7),
    axis.ticks.x = element_line(color = "lightgrey"),
    axis.text.x = element_text(color = "black"),
    axis.title.x = element_text(color = "black"),
    axis.line.y = element_line(color = "grey", linewidth = 0.7),
    axis.ticks.y = element_line(color = "lightgrey")
  )

