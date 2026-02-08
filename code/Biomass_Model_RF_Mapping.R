###########################################################################################
#
#  this script construct random forest based biomass model and apply to canopy hegiht
#  maps to produce plant aboveground biomass maps
#
#    --- Last updated:  2025.03.15 By Daryl Yang <yangd@ornl.gov>
###########################################################################################

#******************** close all devices and delete all variables *************************#
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files
dlm <- .Platform$file.sep # <--- What is the platform specific delimiter?
#*****************************************************************************************#

#****************************** load required libraries **********************************#
### install and load required R packages
list.of.packages <- c("ggplot2", "randomForest", "caTools", "ggpmisc", "terra", 
                      "foreach", "doParallel")  
# check for dependencies and install if needed
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
# load libraries
invisible(lapply(list.of.packages, library, character.only = TRUE))
#*****************************************************************************************#

#************************************ user parameters ************************************#
# define an output directory to store outputs
out.dir <- "***"
# create output directory if not exist
if (! file.exists(out.dir)) dir.create(out.dir,recursive=TRUE)
# creat an temporary to store files temporarily generated during the course of processing
#temp.dir <- file.path(paste0(out.dir, "/", 'temporary'))
#if (! file.exists(temp.dir)) dir.create(temp.dir,recursive=TRUE)

# define how many models you want to build?
nmodel <- 20

# define the output resolution of biomass map: recommend to be at least 10 times the
# chm image resolution
reso <- 1 # m
#*****************************************************************************************#

#*************************************** load data ***************************************#
# define directory to biomass training files
biomss.database.path <- "***"
biomass.dir <- list.files(biomss.database.path, pattern = 'biomass_data_cleaned.csv',
                       full.names = TRUE)
data.org <- read.csv(biomass.dir)
str(data.org)
data.org <- na.omit(data.org)

# load in canopy height model map
chm.path <- "***"
chm.dir <- list.files(chm.path, pattern = 'UAS_SfM.tif$',
                      full.names = TRUE, recursive = TRUE)
chm.rst <- terra::rast(chm.dir)
chm.rst[chm.rst < 0] <- NA
#*****************************************************************************************#

#******************************* model test & evaluation *********************************#
print('Testing Model with leave one out cross valuation (LOOCV)')
# extract data for random forest regression
rf.data <- data.org[, c(15:20)]
# test model performance and uncertainty using LOOCV
rf.out <- c()
for (i in 1:nrow(rf.data))
{
  # hold one sample out for validation
  rf.data.holdout <- rf.data[i, ]
  # select remaining samples as training
  rf.data.train <- rf.data[-i,]
  
  # build multiple model to quantify uncertainty 
  # by randomly selecting 90% of the samples
  holdout.pred <- c()
  for (j in 1:nmodel)
  {
    cat(j)
    # select a portion of training data for each model train 
    train.split <- sample.split(rf.data.train, SplitRatio = 0.9) 
    data.train <- subset(rf.data.train, train.split == "TRUE") 
    
    # build a random forest model
    biomass.model <- randomForest(ABG_kg_m2 ~ ., data = data.train,
                                  ntree = 100, mtry = 2, importance = TRUE)
    # apply model to the holdout validation data
    rf.model.val <- predict(biomass.model, newdata = rf.data.holdout)
    holdout.pred <- cbind(holdout.pred, rf.model.val)
  }
  
  # calculate predicted mean and standard deviation from all models 
  holdout.pred <- data.frame(holdout.pred)
  # mean and standard deviation
  val.mean <- apply(holdout.pred, 1, FUN = mean)
  val.sd <- apply(holdout.pred, 1, FUN = sd)
  # store predicted mean and standard deviation
  pred.matr <- cbind(rf.data.holdout$ABG_kg_m2, val.mean, val.sd)
  pred.matr <- data.frame(pred.matr)
  names(pred.matr) <- c('truth', 'predicted', 'unc')
  rf.out <- rbind(rf.out, pred.matr)
}


#### plot model LOOCV validation result
plt.data <- rf.out
# calculate absolute prediction error
mae <- sum(abs(plt.data$predicted-plt.data$truth)/nrow(plt.data))
mae <- format(round(mae, 3), nsmall = 3)

formula <- y ~ x
ggplot(data = plt.data, aes(x = truth, y = predicted)) + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed", size=1.5) +
  geom_errorbar(aes(ymin = predicted-unc, ymax = predicted+unc), color = 'grey') + 
  geom_smooth(method = 'lm', formula = formula, col = 'black') +
  geom_point(size = 1, shape = 21, color = 'black', fill='black') + 
  ylim(c(0, 6)) + xlim(c(0, 6)) +
  labs(x = 'Ground AGB (kg/m2)', y = 'Predicted AGB (kg/m2)') +
  theme(legend.position = 'none') +
  theme(axis.text = element_text(size=12), axis.title=element_text(size=13)) +
  stat_poly_eq(aes(label = paste(..eq.label.., sep = "~~~")), 
               label.x = 0.65, label.y = 0.33,
               eq.with.lhs = "italic(hat(y))~`=`~",
               eq.x.rhs = "~italic(x)",
               formula = formula, parse = TRUE, size = 4, hjust = 0) +
  stat_poly_eq(aes(label = paste(..rr.label.., sep = "~~~")), 
               label.x = 0.65, label.y = 0.25,
               formula = formula, parse = TRUE, size = 4, hjust = 0) +
  annotate('text', x= 4.0, y = 0.9, label = paste0('MAE = ', mae), size = 4, hjust = 0) +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

png.name <- paste0(out.dir, "/",'model_eva.pdf')
ggsave(png.name, plot = last_plot(), width = 12, height = 12, units = 'cm')
#*****************************************************************************************#

#************************************** apply model **************************************#
print(paste0('Applying Model to: ', basename(chm.dir)))
# calculate window based height stats for applying biomass model
window.size <- round(reso/xres(chm.rst))
grid.min <- aggregate(chm.rst, window.size, min, na.rm = TRUE)
grid.max <- aggregate(chm.rst, window.size, max, na.rm = TRUE)
grid.mean <- aggregate(chm.rst, window.size, mean, na.rm = TRUE)
grid.median <- aggregate(chm.rst, window.size, median, na.rm = TRUE)
grid.90percent <- aggregate(chm.rst, window.size, 
                            fun=function(i) quantile(i, probs=0.9,na.rm = TRUE),
                            cores = 6)
# stack window stats to multiple layer raster
grid.stack <- c(grid.min, grid.max, grid.mean, grid.median, grid.90percent) 
rm(grid.min, grid.max, grid.mean, grid.median, grid.90percent)
# convert raster to dataframe for apply random forest model
stack.dataframe <- as.data.frame(grid.stack, xy=TRUE)
coords <- stack.dataframe[, c(1,2)]
predictors <- stack.dataframe[, -c(1,2)]
names(predictors) <- names(data.org[, c(16:20)])

rm(grid.stack)
### apply model on application data
#Setup backend to use many processors
totalCores = detectCores()
#Leave one core to avoid overload your computer
cluster <- makeCluster(10) #totalCores[1]-1
registerDoParallel(cluster)
result <- foreach(i=1:100) %dopar% {
  train.split <- caTools::sample.split(rf.data.train, SplitRatio = 0.9)
  data.train <- subset(rf.data.train, train.split == "TRUE")
  biomass.model <- randomForest::randomForest(ABG_kg_m2 ~ ., data = data.train,
                                              ntree = 90, mtry = 2, importance = TRUE)
  biomass.pred <- predict(biomass.model, newdata = predictors)
}
pred.df <- data.frame(result)

### calculate mean biomass prediction
pred.mean <- apply(pred.df, 1, FUN = mean, na.rm = TRUE)
# convert biomass dataframe to raster
pred.mean <- cbind(coords, pred.mean)
biomass.pred.rst <- terra::rast(pred.mean, type = 'xyz')
terra::crs(biomass.pred.rst) <- terra::crs(chm.rst)
biomass.pred.rst[biomass.pred.rst == 0] <- NA
basename <- basename(chm.dir)
outname <- paste0(out.dir, '/', 'AGB_', reso, 'm_', basename)
writeRaster(biomass.pred.rst, outname, overwrite=TRUE)

### calculate biomass prediction uncertainty
pred.unc <- apply(pred.df, 1, FUN = sd, na.rm = TRUE)
# convert biomass dataframe to raster
pred.unc <- cbind(coords, pred.unc)
biomass.unc.rst <- terra::rast(pred.unc, type = 'xyz')
terra::crs(biomass.unc.rst) <- terra::crs(chm.rst)
basename <- basename(chm.dir)
outname <- paste0(out.dir, '/', 'AGB_UNC_', reso, 'm_', basename)
writeRaster(biomass.unc.rst, outname, overwrite=TRUE)
#*****************************************************************************************#













