####Skaff et al. Thermal thresholds incrase the vulnerability of coastal Los Angeles to temperature-linked increases in West Nile virus transmission####

#####Loading packages and functions####
library(lubridate)
library(randomForest)
library(dplyr)
library(forestFloor)
library(ROSE)
require(pROC)
library(googledrive)
library(parallel)
library(foreach)
library(ggplot2)
library(sf)
library(plyr)
library(readr)
library(stringr)
library(stargazer)
library(cowplot)
library(colorspace)
library(sp)
library(caret)
library(RColorBrewer)
library(plotly)
library(fitdistrplus)
library(tidyr)
library(mgcv)
library(scales)
library(directlabels)
library(zoo)
library(raster)
library(mediation)
library(scoring)
library(wesanderson)
library(standardize)
library(smoothr)
library(rgdal)
library(leaflet)
library(mapview)
library(medflex)
library(MASS)
library(jpeg)
library(extrafont)
library(rsq)

####Downloading data####
proj_system<-"+init=epsg:3310"
la_metro<-load_sf_from_googledrive("la_metro_shapefile")
la_metro_proj<-spTransform(as_Spatial(la_metro), CRS(proj_system))

#Calsurv mosquito trapping data
calsurv_summed_pools_all_predsLA<- load_csv_from_googledrive("data/calsurv")


####Preparing data####

#for quinquefasciatus abundnace
calsurv_quin<-calsurv_summed_pools_all_predsLA[calsurv_summed_pools_all_predsLA$species=="Culex quinquefasciatus",]

#labeling with the section of LA each mosq point falls in
mosq_points<-SpatialPoints(calsurv_quin[,c("longitude","latitude")], 
                           proj4string =CRS("+proj=longlat +datum=NAD83 +no_defs") )



#making wnv pres variable
calsurv_quin[which(calsurv_quin$num_wnv_pos ==0),"wnv_pres"] <- 0
calsurv_quin[which(calsurv_quin$num_wnv_pos >0),"wnv_pres"] <- 1



###adding value of 1 for NA num_trap
calsurv_quin[which(is.na(calsurv_quin$num_trap)),"num_trap"]<-1

calsurv_quin[which(is.na(calsurv_quin$trap_nights)),"trap_nights"]<-1



####generating temporal predictors####
calsurv_quin$year_num<-year(calsurv_quin$date)

calsurv_quin$month_num<-month(calsurv_quin$date)

calsurv_quin$week_num<-week(calsurv_quin$date)

calsurv_quin$days_since_2000<-as.numeric(difftime(calsurv_quin$date,"2000-01-01",units = "days"))

#removing data with NA, data before 2003, and data with more than 1 trap night or number of traps operated
calsurv_quin<-calsurv_quin[!is.na(calsurv_quin$wnv_pres) & calsurv_quin$year_num>2003 & calsurv_quin$trap_nights==1 &calsurv_quin$num_trap==1 ,]


#code trap types
calsurv_quin[is.na(calsurv_quin$trap_type), "trap_type"]<-"Unknown trap type"
calsurv_quin$trap_type<-as.factor(calsurv_quin$trap_type)

##removing some trap types with very few observations
calsurv_quin<-calsurv_quin[!(calsurv_quin$trap_type %in% c("OTHER", "USDS", "Unknown trap type","REST","BGSENT")),]

#one hot encoding traps
oneHotMdl=function(x) {
  trap = factor(x)
  model.matrix(~trap+0)
}

trap_hotencode<-oneHotMdl(as.factor(as.character(calsurv_quin$trap_type)))
calsurv_quin<-cbind(calsurv_quin,trap_hotencode)


#removing some agencies with very few observations
calsurv_quin<-calsurv_quin[!(calsurv_quin$agency_code %in% c("PASA", "NWST", "SGVA")),]

#one hot encoding agency types
calsurv_quin[is.na(calsurv_quin$agency_code), "agency_code"]<-"Unknown"
calsurv_quin$agency_code<-as.factor(calsurv_quin$agency_code)

oneHotMdl=function(x) {
  agency = factor(x)
  model.matrix(~agency+0)
}

agency_hotencode<-oneHotMdl(calsurv_quin$agency_code)
calsurv_quin<-cbind(calsurv_quin,agency_hotencode)


#removing all virus detection methods except RTPCR
calsurv_quin<-calsurv_quin[calsurv_quin$RTPCR_test_pools>0,]

#removing observations where RTPCR may have differed from RAMP 
calsurv_quin<-calsurv_quin[!(calsurv_quin$RAMP_test_pools>0 & calsurv_quin$num_wnv_pos==1),]


##calculating integrated column soil moisture at different lags
calsurv_quin$soil_moist_integrated_week_lag1<-rowSums(calsurv_quin[,c("soil_moist_level_1_week_lag1","soil_moist_level_2_week_lag1","soil_moist_level_3_week_lag1")])

calsurv_quin$soil_moist_integrated_week_lag2<-rowSums(calsurv_quin[,c("soil_moist_level_1_week_lag2","soil_moist_level_2_week_lag2","soil_moist_level_3_week_lag2")])

calsurv_quin$soil_moist_integrated_week_lag3<-rowSums(calsurv_quin[,c("soil_moist_level_1_week_lag3","soil_moist_level_2_week_lag3","soil_moist_level_3_week_lag3")])

calsurv_quin$soil_moist_integrated_month_lag1<-rowSums(calsurv_quin[,c("soil_moist_level_1_month_lag1","soil_moist_level_2_month_lag1","soil_moist_level_3_month_lag1")])

calsurv_quin$soil_moist_integrated_month_lag2<-rowSums(calsurv_quin[,c("soil_moist_level_1_month_lag2","soil_moist_level_2_month_lag2","soil_moist_level_3_month_lag2")])

calsurv_quin$soil_moist_integrated_month_lag3<-rowSums(calsurv_quin[,c("soil_moist_level_1_month_lag3","soil_moist_level_2_month_lag3","soil_moist_level_3_month_lag3")])

calsurv_quin$soil_moist_integrated_season_lag1<-rowSums(calsurv_quin[,c("soil_moist_level_1_season_lag1","soil_moist_level_2_season_lag1","soil_moist_level_3_season_lag1")])

calsurv_quin$soil_moist_integrated_season_lag2<-rowSums(calsurv_quin[,c("soil_moist_level_1_season_lag2","soil_moist_level_2_season_lag2","soil_moist_level_3_season_lag2")])

calsurv_quin$soil_moist_integrated_season_lag3<-rowSums(calsurv_quin[,c("soil_moist_level_1_season_lag3","soil_moist_level_2_season_lag3","soil_moist_level_3_season_lag3")])


####Generating a dataset to use for random forest####

#randomForest for WNV Presence
calsurv_quin_rf<-calsurv_quin
calsurv_quin_rf$wnv_pres <- as.factor(calsurv_quin_rf$wnv_pres)


#removing covariates with all zeros
calsurv_quin_rf<-calsurv_quin_rf[,colSums(calsurv_quin_rf != 0, na.rm=T) > 0]

#removing unused predictor variables
calsurv_quin_rf1<-calsurv_quin_rf[,c(-grep("clc_pct", colnames(calsurv_quin_rf)),
                                     -grep("slope", colnames(calsurv_quin_rf)),
                                     -grep("wet", colnames(calsurv_quin_rf)),
                                     -grep("canal", colnames(calsurv_quin_rf)),
                                     -grep("catch", colnames(calsurv_quin_rf)),
                                     -grep("maintenance", colnames(calsurv_quin_rf)),
                                     -grep("Spreading", colnames(calsurv_quin_rf)),
                                     -grep("nat_drain", colnames(calsurv_quin_rf)),
                                     -grep("pumpstation", colnames(calsurv_quin_rf)),
                                     -grep("lag7", colnames(calsurv_quin_rf)),
                                     -grep("lag6", colnames(calsurv_quin_rf)),
                                     -grep("lag5", colnames(calsurv_quin_rf)),
                                     -grep("degree_days_week_lag4" , colnames(calsurv_quin_rf)),
                                     -grep("degree_days_week_lag3" , colnames(calsurv_quin_rf)),
                                     -grep("_50", colnames(calsurv_quin_rf)),
                                     -grep("culvert", colnames(calsurv_quin_rf)),
                                     -grep("degree_days", colnames(calsurv_quin_rf)),
                                     -grep("soil_moist_level", colnames(calsurv_quin_rf)),
                                     -grep("lag4" , colnames(calsurv_quin_rf))
                                     )]

#Adding 2 important wetland variables and integrated seasonal soil mositure aggregations to the dataset used for RF predictions
calsurv_quin_rf1<-data.frame(calsurv_quin_rf1,calsurv_quin_rf[,c(
  grep("wetpalus_area", colnames(calsurv_quin_rf)),
  grep("wetriverine_area", colnames(calsurv_quin_rf)) )])

#removing unnecessary buffers from dataset
calsurv_quin_rf1<-calsurv_quin_rf1[,-grep("_50", colnames(calsurv_quin_rf1))]

calsurv_quin_rf<-calsurv_quin_rf1

#setting seed and creating validation and tuning datasets
set.seed(1)
train.index <- createDataPartition(calsurv_quin_rf$wnv_pres,times=1, p = .2, list = F)


##creating variables and shapefiles with unique mosquito sampling sites
unique_sites<-calsurv_quin_rf[-train.index,c("site.ID","latitude","longitude")] %>% 
  dplyr::group_by(site.ID) %>% 
  dplyr::summarize(latitude=first(latitude), longitude=first(longitude), obs_num=n())

unique_sites_sp<-SpatialPointsDataFrame(unique_sites[,c("longitude","latitude")],data=unique_sites[,c("site.ID","obs_num")], proj4string =CRS("+proj=longlat +datum=NAD83 +no_defs") )

unique_sites_sp_trans<-spTransform(unique_sites_sp, CRSobj = CRS(proj_system) )



#removing unnecessary covariates. don't need both grvd and co2 trap types since they're just opposites of each other
calsurv_quin_rf_cov<- calsurv_quin_rf[,c(7,8,15,19,26:99,101:103,106,107:126)]

#converting trap types to factor
calsurv_quin_rf_cov[,"trapGRVD"]<-as.factor(calsurv_quin_rf_cov[,"trapGRVD"])



####Tuning random forest in parallel####


no_cores <- detectCores() - 1

cl <- makeCluster(no_cores)
doParallel::registerDoParallel(cl)


#grid search values
samp_prop<-expand.grid(wnv=c(.25,.5,.75,1), nownv=c(.25,.5,.75,1), trees=c(500), mtry=c(103))
#samp_prop<-expand.grid(wnv=c(.25), nownv=c(.25),  trees=c(500))

tuning_data<-data.frame()
system.time(
  for (i in 1:nrow(samp_prop)){
quin_rf_wnvpres <-  randomForest(
  calsurv_quin_rf_cov[train.index,],
  calsurv_quin_rf$wnv_pres[train.index],
  strata=calsurv_quin_rf$wnv_pres[train.index],
  sampsize = c(ceiling(as.numeric(samp_prop$nownv[i]*table(calsurv_quin_rf$wnv_pres[train.index])[2])),ceiling(as.numeric(samp_prop$wnv[i]*table(calsurv_quin_rf$wnv_pres[train.index])[2]))),
  ntree=samp_prop$trees[i],
  importance=T, 
  keep.inbag = T,
  mtry=samp_prop$mtry[i]
  #mtry = length(calsurv_quin_rf_cov[train.index,]) 
  )

    
    rf.roc<-roc(data.frame(calsurv_quin_rf)$wnv_pres[train.index],quin_rf_wnvpres$votes[,2])
    AUC<-auc(rf.roc)
    quin_rf_wnvpres_pred <- predict(quin_rf_wnvpres)
    conf<-caret::confusionMatrix(quin_rf_wnvpres_pred, calsurv_quin_rf$wnv_pres[train.index], positive="1")
    
    tuning_data<-rbind(tuning_data,data.frame(
      mtry=samp_prop$mtry[i],

      wnv_prop=samp_prop$wnv[i],ntree=samp_prop$trees[i],nownv_prop=samp_prop$nownv[i], auc=AUC, sens=conf$byClass[1], spec=conf$byClass[2], pos_pred=conf$byClass[3],neg_pred=conf$byClass[4], acc=conf$overall[1], balacc=conf$byClass[11], brier=mean(brierscore(calsurv_quin_rf$wnv_pres[train.index]~ predict( quin_rf_wnvpres, type="prob")), na.rm=T)))
    print(i)
  }
  
)

print("Order by smallest difference between sens and spec")
print(tuning_data[order(abs(tuning_data$sens-tuning_data$spec)),])

stopCluster(cl)


####Testing model on validation set -- Random OOB CV####

##testing model on validation set
#randomized prediction on test set
validation_data<-data.frame()
no_cores <- detectCores() - 1
ntree<-c(71,71,71,71,72,72,72) #500 tree


cl <- makeCluster(no_cores)
doParallel::registerDoParallel(cl)
system.time(
  quin_rf_wnvpres_testrand <- foreach(ntree=ntree, 
                                      .combine=randomForest::combine, .packages='randomForest') %dopar% 
    randomForest(calsurv_quin_rf_cov[-train.index,],
                 calsurv_quin_rf$wnv_pres[-train.index],
                 strata=calsurv_quin_rf$wnv_pres[-train.index],
                 sampsize = c(ceiling(as.numeric(.25*table(calsurv_quin_rf$wnv_pres[-train.index])[2])),
                              ceiling(as.numeric(.25*table(calsurv_quin_rf$wnv_pres[-train.index])[2]))),
                 ntree=ntree, 
                 importance=T, 
                 keep.inbag = T,
                 mtry = length(calsurv_quin_rf_cov[-train.index,]) )
  
)

quin_rf_wnvpres_pred_rand_prob<-predict(quin_rf_wnvpres_testrand, type="prob")

rf.roc<-roc(data.frame(calsurv_quin_rf)$wnv_pres[-train.index],quin_rf_wnvpres_testrand$votes[,2])
AUC<-auc(rf.roc)
quin_rf_wnvpres_pred_rand <- predict(quin_rf_wnvpres_testrand)
conf<-caret::confusionMatrix(quin_rf_wnvpres_pred_rand, calsurv_quin_rf$wnv_pres[-train.index], positive="1")

validation_data<-plyr::rbind.fill(validation_data,data.frame(model="random_cv",
                                                             mtry=length(calsurv_quin_rf_cov[-train.index,]),
                                                             trees=sum(ntree),
                                                             auc=AUC, 
                                                             sens=conf$byClass[1], 
                                                             spec=conf$byClass[2], 
                                                             pos_pred=conf$byClass[3],
                                                             neg_pred=conf$byClass[4], 
                                                             acc=conf$overall[1], 
                                                             balacc=conf$byClass[11]))


stopCluster(cl)
print(validation_data)


####Testing model on validation set#### 
#--CV Blocked on Year


###Leave one out year
LOO_year_preds<-data.frame()

cl <- makeCluster(no_cores)
doParallel::registerDoParallel(cl)

unique_years<-unique(calsurv_quin_rf[-train.index,"year_num"])
ntree<-c(72,72,72,72,72,72,68)#500trees



system.time(
  for (i in 1:length(unique_years)){
    
    df<-calsurv_quin_rf[-train.index,]
    #selecting rows in the test year that have site IDs in non-test years so that they can be trained on
    df_indx<-which(!(calsurv_quin_rf[-train.index,]$year_num %in% data.frame(unique_years)[i,]))
    cov_df<-calsurv_quin_rf_cov[-train.index,]
    cov_df_train<-cov_df[df_indx,]
    cov_df_test<-cov_df[-df_indx,]
    resp_df<-calsurv_quin_rf[-train.index,]$wnv_pres
    resp_df_train<-resp_df[df_indx]
    resp_df_test<-resp_df[-df_indx]
    
    
    
    quin_rf_wnvpres_testyear <- foreach(ntree=ntree, 
                                        .combine=randomForest::combine, .packages='randomForest') %dopar%  
      randomForest(cov_df_train,resp_df_train,
                   strata=resp_df_train,
                   sampsize = c(ceiling(as.numeric(.25*table(resp_df_train)[2])),
                                ceiling(as.numeric(.25*table(resp_df_train)[2]))),
                   ntree=ntree, 
                   importance=T, 
                   keep.inbag = T,
                   mtry = length(calsurv_quin_rf_cov[-train.index,]) )
    
    quin_rf_wnvpres_pred_rand_vote <- predict(quin_rf_wnvpres_testyear,cov_df_test ,type="prob")
    quin_rf_wnvpres_pred_rand <- predict(quin_rf_wnvpres_testyear, cov_df_test)
    
    LOO_year_preds<-rbind.fill(LOO_year_preds,
                               data.frame(year=unique_years[i], 
                                          obs=resp_df_test, 
                                          votes=quin_rf_wnvpres_pred_rand_vote[,2], 
                                          preds=quin_rf_wnvpres_pred_rand ))
    
    
    print(i)
  }
)
stopCluster(cl)

#CALCULATE AUC AND CONFusion matrix


rf.roc<-roc(LOO_year_preds$obs,LOO_year_preds$votes)
AUC<-auc(rf.roc)
conf<-caret::confusionMatrix(LOO_year_preds$preds, LOO_year_preds$obs, positive="1")

validation_data<-plyr::rbind.fill(validation_data,
                                  data.frame(model="year_cv",
                                             trees=sum(ntree),
                                             mtry=length(calsurv_quin_rf_cov[-train.index,]) ,
                                             auc=AUC, sens=conf$byClass[1], 
                                             spec=conf$byClass[2], 
                                             pos_pred=conf$byClass[3],
                                             neg_pred=conf$byClass[4], 
                                             acc=conf$overall[1], 
                                             balacc=conf$byClass[11]))

print(validation_data)


###Leave one out month
LOO_month_preds<-data.frame()

cl <- makeCluster(no_cores)
doParallel::registerDoParallel(cl)

unique_months<-as.Date(unique(cut(calsurv_quin_rf[-train.index,"date"], "months")))
ntree<-c(29,29,29,29,29,29,26) #200 tree



system.time(
  for (i in 1:length(unique_months)){
    
    df<-calsurv_quin_rf[-train.index,]
    #selecting rows in the test month that have site IDs in non-test months so that they can be trained on
    df_indx<-which(!(as.Date(cut(calsurv_quin_rf[-train.index,]$date, "months")) %in% data.frame(unique_months)[i,]))
    cov_df<-calsurv_quin_rf_cov[-train.index,]
    cov_df_train<-cov_df[df_indx,]
    cov_df_test<-cov_df[-df_indx,]
    resp_df<-calsurv_quin_rf[-train.index,]$wnv_pres
    resp_df_train<-resp_df[df_indx]
    resp_df_test<-resp_df[-df_indx]
    
    
    quin_rf_wnvpres_testmonth <- foreach(ntree=ntree,
                                         .combine=randomForest::combine, .packages='randomForest') %dopar% 
      randomForest(cov_df_train,
                   resp_df_train,
                   strata=resp_df_train,
                   sampsize = c(ceiling(as.numeric(.25*table(resp_df_train)[2])),
                                ceiling(as.numeric(.25*table(resp_df_train)[2]))),
                   ntree=ntree, 
                   importance=T, 
                   keep.inbag = T,
                   mtry = length(calsurv_quin_rf_cov[-train.index,]) )
    
    quin_rf_wnvpres_pred_rand_vote <- predict(quin_rf_wnvpres_testmonth,cov_df_test ,type="prob")
    quin_rf_wnvpres_pred_rand <- predict(quin_rf_wnvpres_testmonth, cov_df_test)
    
    LOO_month_preds<-rbind.fill(LOO_month_preds,data.frame(month=unique_months[i], obs=resp_df_test, votes=quin_rf_wnvpres_pred_rand_vote[,2], preds=quin_rf_wnvpres_pred_rand ))
    
    
    print(i)
  }
)
stopCluster(cl)

#CALCULATE AUC AND CONFusion matrix


rf.roc<-roc(LOO_month_preds$obs,LOO_month_preds$votes)
AUC<-auc(rf.roc)
conf<-caret::confusionMatrix(LOO_month_preds$preds, LOO_month_preds$obs, positive="1")

validation_data<-plyr::rbind.fill(validation_data,
                                  data.frame(model="month_cv",
                                             trees=sum(ntree),
                                             mtry=length(calsurv_quin_rf_cov[-train.index,]) ,
                                             auc=AUC, sens=conf$byClass[1], 
                                             spec=conf$byClass[2], 
                                             pos_pred=conf$byClass[3],
                                             neg_pred=conf$byClass[4], 
                                             acc=conf$overall[1], 
                                             balacc=conf$byClass[11]))

print(validation_data)




####Testing model on validation set#### 
#--CV Blocked on Site

#going to do blocked cross validation for surveillance sites
LOO_site_preds<-data.frame()


#creating distance matrix
dist_mat<-data.frame(as.matrix(dist(unique_sites_sp_trans@coords, diag=T, upper=F)))
colnames(dist_mat)<- unique_sites_sp$site.ID
rownames(dist_mat)<- unique_sites_sp$site.ID


#subsetting distance matrix to only include sites with more than 10 measurements in order to be able to get a robust measure of error (AUC) for the site

num_samples<-calsurv_quin_rf[-train.index,] %>% 
  group_by(site.ID) %>% 
  dplyr::summarize(num_samples=length(unique(date)))
num_samples_10<-data.frame(num_samples[num_samples$num_samples>=10,"site.ID"])[,1]

dist_mat<-dist_mat[,colnames(dist_mat) %in% num_samples_10]

n_sites<-length(num_samples_10)
site_samp<-sample(1:length(dist_mat), size=n_sites, replace=F)

ntree<-c(29,29,29,29,29,29,26) #200 tree

cl <- makeCluster(no_cores)
doParallel::registerDoParallel(cl)

system.time(
  for (i in 1:length(site_samp)){
    
    far_sites<-rownames(dist_mat[which(dist_mat[site_samp[i]]>=10000),])
    df_indx<-which((calsurv_quin_rf[-train.index,]$site.ID %in% far_sites))
    test_indx<-which((calsurv_quin_rf[-train.index,]$site.ID %in% colnames(dist_mat[site_samp[i]])))
    cov_df<-calsurv_quin_rf_cov[-train.index,]
    cov_df_train<-cov_df[df_indx,]
    cov_df_test<-cov_df[test_indx,]
    resp_df<-calsurv_quin_rf[-train.index,]$wnv_pres
    resp_df_train<-resp_df[df_indx]
    resp_df_test<-resp_df[test_indx]
    
    
    quin_rf_wnvpres_testsite <- foreach(ntree=ntree, 
                                        .combine=randomForest::combine, .packages='randomForest') %dopar%  
      randomForest(cov_df_train,resp_df_train,
                   strata=resp_df_train,sampsize = c(ceiling(as.numeric(.25*table(resp_df_train)[2])),
                                                     ceiling(as.numeric(.25*table(resp_df_train)[2]))),
                   ntree=ntree, 
                   importance=T, 
                   keep.inbag = T,
                   mtry = length(cov_df_train))
    
    quin_rf_wnvpres_pred_rand_vote <- predict(quin_rf_wnvpres_testsite,cov_df_test ,type="prob")
    
    quin_rf_wnvpres_pred_rand <- predict(quin_rf_wnvpres_testsite, cov_df_test)
    
    LOO_site_preds<-rbind.fill(LOO_site_preds,
                               data.frame(site=colnames(dist_mat[site_samp[i]]), 
                                          obs=resp_df_test, 
                                          votes=quin_rf_wnvpres_pred_rand_vote[,2], 
                                          preds=quin_rf_wnvpres_pred_rand ))
    
    print(i)  
  }
)

stopCluster(cl)

#CALCULATE AUC AND CONFusion matrix
rf.roc<-roc(LOO_site_preds$obs,LOO_site_preds$votes)
AUC<-auc(rf.roc)
conf<-caret::confusionMatrix(LOO_site_preds$preds, LOO_site_preds$obs, positive="1")

validation_data<-plyr::rbind.fill(validation_data[1:2,],data.frame(model="site_cv",trees=200,mtry=length(cov_df_train),auc=AUC, sens=conf$byClass[1], spec=conf$byClass[2], pos_pred=conf$byClass[3],neg_pred=conf$byClass[4], acc=conf$overall[1], balacc=conf$byClass[11]))




####Forest Floor model on validation set -- generates feature contributions/marginal effects####

ff_quin = forestFloor(
  rf.fit = quin_rf_wnvpres_testrand,       # mandatory
  X=calsurv_quin_rf_cov[-train.index,],
  y=calsurv_quin_rf$wnv_pres[-train.index],
  calc_np = FALSE,    # TRUE or FALSE both works, makes no difference
  binary_reg = T #takes no effect here when rfo$type="regression"
)


#combining with original dataset
facet_data<-data.frame(data.frame(calsurv_quin_rf[-train.index,]),ff_quin$FCmatrix)

facet_data$predicted<-predict(quin_rf_wnvpres_testrand, type="prob")[,2]


####Discretization of temperature####

facet_data$tmean_month_lag1_colors<-cut(facet_data$tmean_month_lag1, c(-1,21,22.7,40))

facet_data$tmin_month_lag1_colors<-cut(facet_data$tmin_month_lag1, c(0,14.8,15.4,30))

my_brewer = rev(brewer.pal(n = 11, "Spectral")[c(2,4,10)])
my_brewer2 = rev(brewer.pal(n = 11, "Spectral")[c(2,4,10)])




####Creating regions -- Clustering surveillance points on map based on feature contribution####


allfiles<-list.files("data/TopoWx", full.names=T)#topowx temperature rasters

max_files<-allfiles[grep( "max", allfiles)]
min_files<-allfiles[grep( "min", allfiles)]



begin<-Sys.time()
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)


clusterEvalQ(cl, 
             {library("dplyr")
               library("tidyr")
               library("raster")
               library("SearchTrees")
               library("stringr")
             } )


clusterExport(cl, c( "la_metro"), envir=environment())


topowx_data<-parLapply(cl=cl, as.list(allfiles), function(x){
  temp<-brick(x)
  temp_la<-crop(temp, extent(la_metro))
}
)


#all months
#subsetting to 2006 and for June-october
topowx_data_sum<-list()
year_indx<-c(7:17,24:34)
for (i in 1:length(year_indx)){
  topowx_data_sum[[i]]<- topowx_data[[year_indx[i]]][[which(month(as.Date(
    gsub("\\.","-",gsub("X","",topowx_data[[1]]@data@names)))) %in% c(6:10))]]
  
}

max_normals<-mean(stack(topowx_data_sum[1:11]))
min_normals<-mean(stack(topowx_data_sum[12:22]))

mean_normals<-(max_normals+min_normals)/2

meannorm_la<-mask(mean_normals, la_metro)



## kmeans classification 
v <- getValues(meannorm_la)
i <- which(!is.na(v))
v <- na.omit(v)



E <- kmeans(v, 3, iter.max = 1000, nstart = 100)

E$cluster[E$cluster==rownames(E$centers)[which(E$centers==min(E$centers))]]<-4
E$cluster[E$cluster==rownames(E$centers)[which(E$centers!=max(E$centers) & E$centers!=min(E$centers))]]<-5
E$cluster[E$cluster==rownames(E$centers)[which(E$centers==max(E$centers))]]<-6


E$cluster<-E$cluster-3

kmeans_raster <- raster(meannorm_la)


kmeans_raster[i] <- E$cluster
plot(kmeans_raster)

range(v[E$cluster==1])
range(v[E$cluster==2])
range(v[E$cluster==3])



poly_region<-rasterToPolygons(kmeans_raster, dissolve=T)



#adding region to analysis dataset
region_rast<-kmeans_raster


site_region<-raster::extract(x=region_rast, y=unique_sites_sp)

#for sites with NA finding closest raster pixel
extract_NA <- apply(X = coordinates(unique_sites_sp[is.na(site_region),]), 
                    MARGIN = 1, 
                    FUN = function(xy) values(region_rast)
                    [which.min(replace(distanceFromPoints(region_rast, xy), is.na(region_rast), NA))])

site_region[is.na(site_region)]<-extract_NA


#adding region to facet data
facet_data<-left_join(facet_data, data.frame(cbind(site.ID=unique_sites_sp$site.ID,cluster=site_region)), by="site.ID")




####Figure 1####


lancet_colors<-c("#00468BFF","#42B540FF", "#ED0000FF" )

poly_region$layer<-factor(poly_region$layer,levels =c("1","2","3")  ,labels=c("Coastal", "Central", "Inland"))

map_regions<-mapview(poly_region, map.types="CartoDB.Positron", trim=T, col.regions=lancet_colors, alpha.regions=.6,na.color="transparent", lwd=.005, legend = TRUE, layer.name="Zones")+mapview(la_metro,  alpha.regions=.001, map.types="CartoDB.Positron", legend=F)

####Table S1####

print(validation_data)



#### Figure S1####
#gini importance
importance_data<-data.frame(randomForest::importance(quin_rf_wnvpres_testrand)[,3:4])
importance_data<-importance_data[rev(order(importance_data$MeanDecreaseGini)),]
importance_data$var_names<-rownames(importance_data)

importance_data$var_names<- factor(importance_data$var_names, levels = rev(importance_data$var_names))


import_gini<-ggplot() + 
  geom_bar(data=importance_data[1:10,], aes(y=MeanDecreaseGini, x=var_names), stat="identity" ) + 
  coord_flip() + 
  scale_x_discrete("", labels=rev(c("Mean Temp-1 month", "# mosq in pool", "Min Temp-1 month", 
                                    "Mean Temp-1 quarter", "# mosq captured","Mean Temp-2 month","Longitude",
                                    "Precip-2 quarter" ,"Min Temp-3 week","Min Temp-2 month"))) + 
  scale_y_continuous("Mean Decrease in Gini Index")+
  theme_bw() +guides(fill=F) + 
  theme(panel.grid = element_blank(), 
        panel.border = element_blank(), 
        axis.ticks.y=element_blank(),
        axis.line.x = element_line(colour = "black"), 
        aspect.ratio = 6/4)



#Variable importance measured by mean decrease in model accuracy when the variable is permuted from the model.

importance_data<-data.frame(randomForest::importance(quin_rf_wnvpres_testrand)[,3:4])
importance_data<-importance_data[rev(order(importance_data$MeanDecreaseAccuracy)),]
importance_data$var_names<-rownames(importance_data)

importance_data$var_names<- factor(importance_data$var_names, levels = rev(importance_data$var_names))

import_acc<-ggplot() + 
  geom_bar(data=importance_data[1:10,], aes(y=MeanDecreaseAccuracy, x=var_names), stat="identity" ) + 
  coord_flip() +
  scale_x_discrete("", labels=rev(c( "# mosq in pool", "Mean Temp-1 month", "Mean Temp-1 quarter",
                                     "Max Temp-1 week","Longitude", "# mosq captured","Diurnal-3 month" ,
                                     "Mean Temp-2 month","Min Temp-2 week","Forest cover-1000m"))) + 
  scale_y_continuous("Mean Decrease in Model Accuracy")+
  theme_bw() +guides(fill=F)+ 
  theme(panel.grid = element_blank(), 
        panel.border = element_blank(), 
        axis.ticks.y=element_blank(),
        axis.line.x = element_line(colour = "black"), 
        aspect.ratio = 6/4)



Supp_Fig1<-plot_grid(import_gini, import_acc, ncol = 2,labels = "auto")



####Figure 2####
#Feature contribution plots showing the marginal effects of minimum temperature and degree days at 1 month lags. 
#figure 2a
dd_fc<-ggplot(data=facet_data) + 
  geom_point(aes(x=tmean_month_lag1, y=tmean_month_lag1.1), alpha=.1, shape=1, color="black")  + 
  scale_x_continuous("Monthly mean temp.(°C)") + 
  scale_y_continuous("Δ Cx. infection probability", limits=c(-.56,.46)) + 
  scale_color_manual("",values=(my_brewer2), labels=c("Inhibitory", "Transitional", "Favorable"))+ 
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"),legend.position = c(.2,.8),legend.text=element_text(size=16)) + 
  geom_smooth(data=facet_data[facet_data$tmean_month_lag1<28,],aes(x=tmean_month_lag1, y=tmean_month_lag1.1, color=tmean_month_lag1_colors), method = "gam", formula=y~s(x,bs="cs", k=100),se = F, size=1.2) + 
  guides(colour = guide_legend(override.aes = list(size=5)))

#determining maximal effect of transitional range
smooth_data<-ggplot_build(dd_fc)$data[[2]]
(smooth_data[smooth_data$x==21.00033,"y"]*-1)+smooth_data[smooth_data$x==22.70017,"y"]

#summing the feature contribution of all temperature predictors
facet_data$temp_sum<-rowSums(ff_quin$FCmatrix[,c(grep("tmean",colnames(ff_quin$FCmatrix)),
                                                 grep("tmin",colnames(ff_quin$FCmatrix)),
                                                 grep("tmax",colnames(ff_quin$FCmatrix)),
                                                 grep("diurnal",colnames(ff_quin$FCmatrix)))])



#generating colors for the most important non-tmean contribution

#extracting all temperatue feature contributions
fc_temp<-ff_quin$FCmatrix[,c(grep("tmean",colnames(ff_quin$FCmatrix)),
                             grep("tmin",colnames(ff_quin$FCmatrix)),
                             grep("tmax",colnames(ff_quin$FCmatrix)),
                             grep("diurnal",colnames(ff_quin$FCmatrix)))]



#figure 2b
temp_fc<-ggplot(data=facet_data) +
  #geom_point(aes(x=tmean_month_lag1, y=temp_sum, color=diurnal_color, alpha=tmax_color)) +
  geom_point(aes(x=tmean_month_lag1, y=temp_sum), alpha=.1, shape=1, color="black") +
  scale_x_continuous("Monthly mean temp. (°C)") + 
  scale_y_continuous("Δ Cx. infection probability", limits=c(-.56,.46))+
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"), legend.position = "top",legend.text=element_text(size=16))+ 
  geom_smooth(data=facet_data[facet_data$tmean_month_lag1<28,],aes(x=tmean_month_lag1, y=temp_sum), method = "gam", formula=y~s(x,bs="cr", k=45), color="red",se = T, size=1.2)+
  #scale_color_manual("Monthly Diurnal Var.", values = c("grey", "red"), labels=c("<18°C",">18°C")) + 
  #scale_alpha_manual("Monthly Diurnal Var.", values = c(.4, .8), labels=c("<18°C",">18°C")) + 
  guides(color = guide_legend(override.aes = list(alpha=1,size=3)),shape = guide_legend(override.aes = list(alpha=1,size=2))) + 
  scale_shape_manual("",breaks=c(1,-1),values=c(2,6),labels=c("Increase","Decrease") ) +
  geom_vline(xintercept = c(21,22.7), linetype="dashed", color="grey10")
#+ ggtitle("Marginal effects - all temp. predictors")


temp_fc_noleg<-temp_fc+theme(legend.position = c(.2,.8))
temp_fc_colleg<-get_legend(temp_fc+guides(shape=F)+theme(legend.justification = "center"))
figure2<-plot_grid(dd_fc,temp_fc_noleg,ncol=1, labels=c("a","b",NA), rel_heights = c(1.2,1.2))


###Figure S2####

facet_data[facet_data$diurnal_month_lag1<=(18),"diurnal_color"]<-"black"
facet_data[facet_data$diurnal_month_lag1> (18),"diurnal_color"]<-"red"


#Figure S2a
temp_fc_diurnal<-ggplot(data=facet_data) +
  geom_point(aes(x=tmean_month_lag1, y=temp_sum, color=diurnal_color, alpha=diurnal_color)) +
  #geom_point(aes(x=tmean_month_lag1, y=temp_sum), alpha=.1, shape=1, color="black") +
  scale_x_continuous("Monthly mean temp. (°C)") + 
  scale_y_continuous("Δ Cx. infection probability", limits=c(-.56,.46))+
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"), legend.position = "top",legend.text=element_text(size=16))+ 
  geom_smooth(data=facet_data[facet_data$tmean_month_lag1<28,],aes(x=tmean_month_lag1, y=temp_sum), method = "gam", formula=y~s(x,bs="cr", k=45), color="red",se = T, size=1.2)+
  scale_color_manual("Monthly Diurnal Var.", values=c("black", "red"),breaks=c("black", "red"), labels=c("<18°C",">18°C")) + 
  scale_alpha_manual("Monthly Diurnal Var.", values = c(.1, .5), labels=c("<18°C",">18°C")) + 
  guides(color = guide_legend(override.aes = list(alpha=1,size=3)),shape = guide_legend(override.aes = list(alpha=1,size=2))) + 
  scale_shape_manual("",breaks=c(1,-1),values=c(2,6),labels=c("Increase","Decrease") ) +
  geom_vline(xintercept = c(21,22.7), linetype="dashed", color="grey10")
#+ ggtitle("Marginal effects - all temp. predictors")



#Figure S2b
diurnal_fc<-ggplot(data=facet_data) + 
  geom_point(aes(x=diurnal_month_lag1, y=diurnal_month_lag1.1+diurnal_season_lag1.1), alpha=.1, shape=1, color="black")  + 
  scale_x_continuous("Monthly diurnal variation (°C)") + 
scale_y_continuous("Δ Cx. infection probability") + 
#scale_color_manual("",values=(my_brewer2), labels=c("Inhibitory", "Transitional", "Favorable"))+ 
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"),legend.position = c(.2,.8),legend.text=element_text(size=16)) + 
  geom_smooth(data=facet_data[facet_data$diurnal_month_lag1<20,],aes(x=diurnal_month_lag1, y=diurnal_month_lag1.1+diurnal_season_lag1.1), method = "gam", formula=y~s(x,bs="cs", k=10),se = F, size=1.2) + 
  guides(colour = guide_legend(override.aes = list(size=5)))
#+ ggtitle("Marginal effects - monthly mean temp.")


figureS2<-plot_grid(diurnal_fc,temp_fc_diurnal,ncol=1, labels=c("a","b"), rel_heights = c(1.2,1.2))

ggsave("Supp_fig_diurnal_7_25.png",figureS2, path="Personal_Folders/Nick/Results", width=6,height=10, units = "in", dpi = 300)



####Figure S3####
#Maps of LA temperature 

#pal<-choose_palette()
pal <- function (n, h = c(-92, 43), c = 100, l = c(51, 68), power = 0.711864406779661, 
                 fixup = TRUE, gamma = NULL, alpha = 1, ...) 
{
  if (!is.null(gamma)) 
    warning("'gamma' is deprecated and has no effect")
  if (n < 1L) 
    return(character(0L))
  h <- rep(h, length.out = 2L)
  c <- c[1L]
  l <- rep(l, length.out = 2L)
  power <- rep(power, length.out = 2L)
  rval <- seq(1, -1, length = n)
  rval <- hex(polarLUV(L = l[2L] - diff(l) * abs(rval)^power[2L], 
                       C = c * abs(rval)^power[1L], H = ifelse(rval > 0, h[1L], 
                                                               h[2L])), fixup = fixup, ...)
  if (!missing(alpha)) {
    alpha <- pmax(pmin(alpha, 1), 0)
    alpha <- format(as.hexmode(round(alpha * 255 + 1e-04)), 
                    width = 2L, upper.case = TRUE)
    rval <- paste(rval, alpha, sep = "")
  }
  return(rval)
}





month_labeller <- function(variable,value){
  return(month.abb[value])
}


#Figure S3a
dd_fc_map<- ggplot(data=facet_data[facet_data$month_num %in% c(6:10),]) + geom_jitter(aes(x=longitude, y=latitude, color=tmean_month_lag1.1), width=.01, height=.01, alpha=.1)+ geom_sf(data=la_metro, fill="transparent") + theme_bw()+ facet_wrap(~month_num, ncol = 1, labeller = month_labeller) + theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),legend.position = "right") + scale_color_distiller("Marginal Effect\nMean Temp\n(1 month) [a]",palette="Spectral")

#Figure S3b
tempsum_fc_map<- ggplot(data=facet_data[facet_data$month_num %in% c(6:10),]) + geom_jitter(aes(x=longitude, y=latitude, color=max_fc$fc_group[facet_data$month_num %in% c(6:10)]), width=.01, height=.01)+ geom_sf(data=la_metro, fill="transparent") + theme_bw()+ facet_wrap(~month_num, ncol = 1, labeller = month_labeller) + theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),legend.position = "right") + scale_color_brewer("Marginal Effect\nAll Temp [b]",palette="Spectral")

#Figure S3c
predicted_map<- ggplot(data=facet_data[facet_data$month_num %in% c(6:10),]) + geom_jitter(aes(x=longitude, y=latitude, color=predicted), width=.01, height=.01, alpha=.1)+ geom_sf(data=la_metro, fill="transparent") + theme_bw()+ facet_wrap(~month_num, ncol = 1, labeller = month_labeller) + theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),legend.position = "right") + scale_color_distiller("Predicted Prob.\nWNV Presence [c]",palette="Spectral")


dd_fc_map_noleg<-dd_fc_map+guides(color=F)
dd_fc_map_leg<-  get_legend(dd_fc_map)
tempsum_fc_map_noleg<-tempsum_fc_map+guides(color=F)
tempsum_fc_map_leg<-  get_legend(tempsum_fc_map)
predicted_fc_map_noleg<-predicted_map+guides(color=F)
predicted_map_leg<-  get_legend(predicted_map)


temp_legends<-plot_grid(dd_fc_map_leg,tempsum_fc_map_leg,predicted_map_leg, ncol=1)

temp_maps<-plot_grid(dd_fc_map_noleg, tempsum_fc_map_noleg,predicted_fc_map_noleg,ncol=3, labels="auto")

figure3<-plot_grid(temp_maps, temp_legends, ncol=2, rel_widths = c(1,.3))






####Calculating projected temperature changes under rcp4.5####
load("data/midD_Nick_Tmean_2040to2069_minus_1996to2025_RCP45_050219.rda")



# #for each combination of first 2 raster dimensions, take mean over the 3rd (selected months)
jul_oct_ave<-apply(midD[,,c("Jul", "Oct"),"multi-model-mean"], c(2,1), mean)

jul_oct_rast<-raster::flip(raster(jul_oct_ave, xmn=-118.8715, xmx=-116.829, ymn=33.37322 , ymx=34.33791 ), direction=2)

projection(jul_oct_rast)<-CRS("+proj=longlat +datum=NAD83 +no_defs")

jul_oct_rast_mask<-raster::mask(jul_oct_rast, mask = la_metro)



###August, sept
aug_sep_ave<-apply(midD[,,c("Aug", "Sep"),"multi-model-mean"], c(2,1), mean)

aug_sep_rast<-raster::flip(raster(aug_sep_ave, xmn=-118.8715, xmx=-116.829, ymn=33.37322 , ymx=34.33791 ), direction=2)

projection(aug_sep_rast)<-CRS("+proj=longlat +datum=NAD83 +no_defs")

aug_sep_rast_mask<-raster::mask(aug_sep_rast, mask = la_metro)






##extracting increases in temp by region

aug_sep_regional<-raster::extract( aug_sep_rast_mask,spTransform(poly_region, CRSobj = CRS(proj4string( aug_sep_rast_mask))), fun=mean, na.rm=T)
aug_sep_regional

jul_oct_regional<-raster::extract( jul_oct_rast_mask,spTransform(poly_region, CRSobj = CRS(proj4string( jul_oct_rast_mask))), fun=mean, na.rm=T)
jul_oct_regional


projected_increases_df<-data.frame(region=c(rep(c("Coastal", "Central", "Inland"), 7)  ), 
                                   month=c("6","6","6","7","7","7","8","8","8","9","9","9","10","10","10",
                                           "aug_sep","aug_sep", "aug_sep", "jul_oct","jul_oct","jul_oct"), 
                                   increase=c(jun_regional[,1],jul_regional[,1],aug_regional[,1],
                                              sep_regional[,1], oct_regional[,1], aug_sep_regional[,1], 
                                              jul_oct_regional[,1]))




####Figure S5####


facet_data$month_fac<-factor(facet_data$month_num, levels =7:10, labels = c("July", "August", "September", "October"))

facet_data$cluster_fac<-factor(facet_data$cluster, levels = 1:3, labels = c("Coastal", "Central", "Inland"))


for (i in 1:nrow(facet_data)){
  if(facet_data$month_num[i] %in% 7:10){
  facet_data[i,"tmean_warmed"]<-facet_data$tmean_month_lag1[i]+projected_increases_df[projected_increases_df$region==facet_data$cluster_fac[i] & projected_increases_df$month==facet_data$month_num[i],"increase"]
  print(i)
  }
  else(next)
}



  transcritical_hist<-ggplot() + 
  geom_histogram(data=facet_data[facet_data$month_num %in% c(7:10), ], 
                 aes(x=tmean_month_lag1), color="white", bins=30, alpha=.5) +
    geom_histogram(data=facet_data[facet_data$month_num %in% c(7:10), ], 
                   aes(x=tmean_warmed), color="white",fill="red", bins=30, alpha=.4) +
  geom_vline(xintercept = c(21,22.7), color="black", size=1, linetype="dashed") +
  facet_grid(cluster_fac~month_fac, scales="free_y" )+
  scale_x_continuous("Monthly mean temp. (°C)" )+
  scale_y_continuous("Number of observations",expand=c(0,0))+
  theme_bw()




##determining % of observations in the favorable range by example years

cluster1_temps<-facet_data$tmean_month_lag1[facet_data$year_num==2011 & facet_data$month_num %in% c(7:10) & facet_data$cluster==1]

length(cluster1_temps[cluster1_temps>22.7])/length(cluster1_temps)


cluster2_temps<-facet_data$tmean_month_lag1[facet_data$year_num==2011 & facet_data$month_num %in% c(7:10) & facet_data$cluster==2]

length(cluster2_temps[cluster2_temps>22.7])/length(cluster2_temps)


cluster3_temps<-facet_data$tmean_month_lag1[facet_data$year_num==2011 & facet_data$month_num %in% c(7:10) & facet_data$cluster==3]

length(cluster3_temps[cluster3_temps>22.7])/length(cluster3_temps)



cluster1_temps<-facet_data$tmean_month_lag1[facet_data$year_num==2014 & facet_data$month_num %in% c(7:10) & facet_data$cluster==1]

length(cluster1_temps[cluster1_temps>22.7])/length(cluster1_temps)


cluster2_temps<-facet_data$tmean_month_lag1[facet_data$year_num==2014 & facet_data$month_num %in% c(7:10) & facet_data$cluster==2]

length(cluster2_temps[cluster2_temps>22.7])/length(cluster2_temps)


cluster3_temps<-facet_data$tmean_month_lag1[facet_data$year_num==2014 & facet_data$month_num %in% c(7:10) & facet_data$cluster==3]

length(cluster3_temps[cluster3_temps>22.7])/length(cluster3_temps)





####Figure S6####

pal <- colorNumeric(
  palette = "Spectral",
  domain = c(1,2.02),
  na.color = "transparent",
  reverse = T
)


#Figure S6a
aug_sep_rast_map_5_24<-leaflet(la_metro)%>% 
  addProviderTiles("CartoDB.Positron")%>%
  addRasterImage(aug_sep_rast, colors = pal, opacity = 0.7) %>%
  addLegend( position="bottomright", pal=pal, na.label = NA, title ="Δ Tmean (ᵒC)<br>Aug-Sept<br>2040-2069", values =aug_sep_rast@data@values, bins=c(1,1.2,1.4,1.6,1.8,2) ) %>%
  addPolygons( fill = NA, color="black", weight = 3)


#Figure S6b
jul_oct_rast_map_5_24<-leaflet(la_metro)%>% 
  addProviderTiles("CartoDB.Positron")%>%
  addRasterImage(jul_oct_rast, colors = pal, opacity = 0.7) %>%
  addLegend( position="bottomright", pal=pal, na.label = NA, title ="Δ Tmean (ᵒC)<br>July/Oct<br>2040-2069", values =aug_sep_rast@data@values, bins=c(1,1.2,1.4,1.6,1.8,2)) %>% #leaving this same as august and sept. so that the legend is the same between plots
  addPolygons( fill = NA, color="black", weight = 3)






####Figure 3####

#Figure 3c
scaleFUN <- function(x) sprintf("%.2f", x)

#generating weekly aggregated smooth line
weekly_data_2011_2014<-facet_data[which(!is.na(facet_data$cluster_fac) & facet_data$year_num %in% c(2011,2014)),]%>%
  group_by(cluster, date=as.Date(cut(date, "week"))) %>% 
  dplyr::summarize(predicted=mean(predicted), 
                   tmean_month_lag1=mean(tmean_month_lag1), 
                   tmean_month_lag1.1=mean(tmean_month_lag1.1))

#Figure 3c
tmean_fc_2011_2014<-ggplot(data=facet_data[which(!is.na(facet_data$cluster_fac)& facet_data$year_num %in% c(2011,2014)),]) +
  geom_point( aes(x=yday(date), y=tmean_month_lag1.1, color=as.factor(cluster)), alpha=.05, size=2,   stroke=1) +
  geom_line( data=weekly_data_2011_2014[year(weekly_data_2011_2014$date)!=2013,],aes(x=yday(date), y=tmean_month_lag1.1, group=cluster, color=as.factor(cluster)), size=.8) +
  facet_wrap(~year(date), ncol = 1) +
  geom_vline(aes(xintercept=c(151)), linetype="dotted")+
  geom_vline(aes(xintercept=c(305)), linetype="dotted")+
  theme_bw()+
  theme(strip.background=element_rect(fill="white"), strip.text = element_text(face = "bold", size=11))+
  scale_x_continuous("Month",breaks=c(1,31,60,91,121,152,182,213,244,274,305, 335), labels=c(1:12), expand = c(0.01, 0.01))+
  scale_y_continuous("Δ Cx. infection probability", labels=scaleFUN, limits=c(-.35,.25)) +
  scale_color_manual("LA Region", values = c("#00468BFF","#42B540FF", "#ED0000FF" ))+
  scale_fill_distiller("LA Region Average \nMonthly Mean Temp. (°C)",palette="Spectral")+
  guides(color=guide_legend(override.aes = list(alpha=1, fill=NA)))


#Figure 3b
#mean temp values not contributions
tmean_2011_2014<-ggplot(data=facet_data[which(!is.na(facet_data$cluster)& facet_data$year_num %in% c(2011,2014)),]) +
  geom_point( aes(x=yday(date), y=tmean_month_lag1, color=as.factor(cluster)), alpha=.035, size=2,   stroke=1) +
  geom_line( data=weekly_data_2011_2014[year(weekly_data_2011_2014$date)!=2013,],aes(x=yday(date), y=tmean_month_lag1, group=cluster, color=as.factor(cluster)), size=.8) +
  geom_rect(aes(ymin=21,ymax=22.3,xmin=-Inf,xmax=Inf),alpha=0.01,fill="#FDAE61")+
  facet_wrap(~year(date), ncol=1) +
  geom_vline(aes(xintercept=c(151)), linetype="dotted") +
  geom_vline(aes(xintercept=c(305)), linetype="dotted") +
  theme_bw()+
  theme(strip.background=element_rect(fill="white"), strip.text = element_text(face = "bold", size=11))+
  scale_x_continuous("Month", breaks=c(1,31,60,91,121,152,182,213,244,274,305, 335), labels=c(1:12), expand = c(0.01, 0.01))+
  scale_y_continuous("Monthly mean temp. (°C)") +
  scale_color_manual("LA Region", values = c("#00468BFF","#42B540FF", "#ED0000FF" ))+
  guides(color=guide_legend(override.aes = list(alpha=1, fill=NA)))





tmean_fc_2011_2014_noleg<-tmean_fc_2011_2014+guides(color=F, fill=F)
tmean_fc_2011_2014_leg<-get_legend(tmean_fc_2011_2014)
tmean_2011_2014_noleg<-tmean_2011_2014+guides(color=F)



#merging just the yearly comparison not overall timeseries yet
fig3_2011_2014_noseries<-plot_grid(tmean_2011_2014_noleg,tmean_fc_2011_2014_noleg,tmean_fc_2011_2014_leg, nrow = 1, labels=c("a","b",NA), rel_widths = c(1,1.1,.4))




#calculating monthly data, including regional anomalies


tmean_week<-facet_data[facet_data$year_num >2005,] %>% 
  group_by(week=week(date), cluster) %>% 
  dplyr::summarize(tmean_ref=mean(tmean_month_lag1))

monthly_anomaly_intensity<-facet_data[facet_data$year_num >2005 & facet_data$month_num %in% 1:12,] %>% 
  group_by(cluster, date=cut(date, "week"), week=week(date))  %>% 
  dplyr::summarise(tmean_1=mean(tmean_month_lag1), 
                   mean_predicted=mean(predicted))%>%
  left_join(.,tmean_week, by=c("week","cluster")) %>% 
  group_by(cluster, month=month(date), year=year(date))%>% 
  dplyr::summarise(tmean_mean=mean(tmean_1), 
                   mean_prob_pred=mean(mean_predicted),
                   tmean_anom=mean(tmean_1)-mean(tmean_ref) )

monthly_anomaly_intensity$cluster_fac<-factor(monthly_anomaly_intensity$cluster, levels=c(1,2,3), 
                                              labels = c("Coastal", "Central", "Inland"))

dates<-data.frame(date=seq(from=as.Date("2006_01_01","%Y_%m_%d"), to=as.Date("2016_12_31","%Y_%m_%d"), by = 1))
dates$month<-month(dates$date)
dates$year<-year(dates$date)



#average monthly anomalies for the whole region
monthly_anomaly_intensity_agg<-monthly_anomaly_intensity %>% 
  group_by(month, year)%>% 
  dplyr::summarise(tmean_anom=mean(tmean_anom))

la_anom<-full_join(dates, monthly_anomaly_intensity_agg, by=c("month","year")) 


weekly_data<-facet_data[which(!is.na(facet_data$cluster_fac)),] %>% 
  group_by(cluster_fac, date=as.Date(cut(date, "week"))) %>% 
  dplyr::summarize(predicted=mean(predicted), 
                   tmean=mean(tmean_month_lag1),
                   tmean_fc=mean(tmean_month_lag1.1), 
                   alltemp_fc=mean(temp_sum) )


#Figure 3a
predicted_fc_fullseries<-ggplot(data=weekly_data)+
 geom_rect(data=la_anom[!is.na(la_anom$tmean_anom),],aes(xmin=date+.001,xmax=date+30,ymin=-.3,ymax=.3, 
                                                         fill=tmean_anom), alpha=.03)+ 
  geom_line( aes(x=date, y=tmean_fc, group=cluster_fac, color=as.factor(cluster_fac)), size=.75) +
  theme_bw()+
  scale_x_date("",breaks=seq(from=as.Date("2006_01_01","%Y_%m_%d"), to=as.Date("2016_01_01","%Y_%m_%d"), 
                             by = "year") ,labels=c(2006:2016),expand = c(0, 0))+
  scale_y_continuous("Δ Cx. infection probability") +
  scale_color_manual("", values = c("#00468BFF","#42B540FF", "#ED0000FF" )) +
  scale_fill_distiller("LA Region Monthly\nMean Temp. Anomaly (°C)",palette="Spectral")+
  guides(color=guide_legend(override.aes = list(alpha=1, fill=NA)))+
  theme(legend.position="top",legend.key.height=unit(.7,"line"))



predicted_fc_fullseries_noleg<-predicted_fc_fullseries+guides(color=F, fill=F)
predicted_fc_fullseries_leg<-get_legend(predicted_fc_fullseries)




fig3<-plot_grid(predicted_fc_fullseries,
                              plot_grid(tmean_2007_2015_noleg,tmean_fc_2007_2015_noleg, 
                                        ncol=2, labels=c("b","c"), 
                                        rel_widths = c(1,1.09)),
                              rel_heights = c(1,1.3) ,
                              nrow=2, 
                              labels=c("a", NA))



####Figure S4####

#Figure S4a
tmean_fullseries<-ggplot(data=weekly_data)+
  geom_rect(aes(ymin=21,ymax=22.3,xmin=as.Date("2006-01-01")+1,xmax=as.Date("2016-12-31")-1),
            alpha=0.008,fill="#FDAE61")+
  geom_line( aes(x=date, y=tmean, group=cluster_fac, color=as.factor(cluster_fac)), size=.75) +
  theme_bw()+
  scale_x_date("",breaks=seq(from=as.Date("2006_01_01","%Y_%m_%d"), to=as.Date("2016_01_01","%Y_%m_%d"), 
                             by = "year") ,labels=c(2006:2016),expand = c(0, 0))+
  scale_y_continuous("Monthly mean temp. (°C)") +
  scale_color_manual("", values = c("#00468BFF","#42B540FF", "#ED0000FF" )) +
  guides(color=guide_legend(override.aes = list(alpha=1, fill=NA)))+
  theme(legend.position="top",legend.key.height=unit(.7,"line"))

#Figure S4b
predicted_fullseries<-ggplot(data=weekly_data)+
  geom_line( aes(x=date, y=predicted, group=cluster_fac, color=as.factor(cluster_fac)), size=.75) +
  theme_bw()+
  scale_x_date("",breaks=seq(from=as.Date("2006_01_01","%Y_%m_%d"), to=as.Date("2016_01_01","%Y_%m_%d"), by = "year") ,labels=c(2006:2016),expand = c(0, 0))+
  scale_y_continuous("Cx. infection probability", limits=c(0,.8)) +
  scale_color_manual("", values = c("#00468BFF","#42B540FF", "#ED0000FF" )) +
  guides(color=guide_legend(override.aes = list(alpha=1, fill=NA)))+
  theme(legend.position="top",legend.key.height=unit(.7,"line"))




pred_legend_fullseries<-get_legend(predicted_fullseries)
predicted_fullseries_noleg<-predicted_fullseries+guides(color=F)

figure_S4<-plot_grid(pred_legend_fullseries, tmean_fullseries+guides(color=F),predicted_fullseries_noleg, nrow=3, rel_heights = c(.2,1,1), labels = c(NA, "a","b"))






#monthly human cases
human_cases_monthly<-read.csv("data/human_la_aggregated_monthly.csv", header=T)
monthly_human_anom<-full_join(monthly_anomaly_intensity, human_cases_monthly, by=c("cluster","month","year"))
monthly_human_anom[is.na(monthly_human_anom$incidence), c("incidence", "cases")]<-0


monthly_human_anom<-data.frame(monthly_human_anom[!is.na(monthly_human_anom$tmean_mean),])
monthly_human_anom$cluster<-as.factor(monthly_human_anom$cluster)


model.Y<-glm.nb(cases~ tmean_mean*cluster+mean_prob_pred+ offset(log((pop2010))),data=data.frame(monthly_human_anom[monthly_human_anom$year>2005 & monthly_human_anom$month %in% 8:9 ,]), control = glm.control(maxit = 50))



###Intermediary models for Aug/Sept
model.M<-lm(mean_prob_pred~tmean_mean*relevel(cluster, ref=1),data=data.frame(monthly_human_anom[monthly_human_anom$year>2005& monthly_human_anom$month %in% 8:9,])) 

#same model just changing he reference level to a different climate zone
model.M_cent<-lm(mean_prob_pred~tmean_mean*relevel(cluster, ref=2),data=data.frame(monthly_human_anom[monthly_human_anom$year>2005& monthly_human_anom$month %in% 8:9,])) 

model.M_inland<-lm(mean_prob_pred~tmean_mean*relevel(cluster, ref=3),data=data.frame(monthly_human_anom[monthly_human_anom$year>2005& monthly_human_anom$month %in% 8:9,])) 


###Intermediary models for July/October
model.M_JO<-lm(mean_prob_pred~tmean_mean,data=data.frame(monthly_human_anom[monthly_human_anom$year>2005& monthly_human_anom$month %in% c(7,10),])) 

#same model just changing he reference level to a different climate zone
model.M_JO_cent<-lm(mean_prob_pred~tmean_mean*relevel(cluster, ref=2),data=data.frame(monthly_human_anom[monthly_human_anom$year>2005& monthly_human_anom$month %in% c(7,10),])) 

model.M_JO_inland<-lm(mean_prob_pred~tmean_mean*relevel(cluster, ref=3),data=data.frame(monthly_human_anom[monthly_human_anom$year>2005& monthly_human_anom$month %in% c(7,10),])) 



####Mediation analysis####
#mediation analysis for august + september
#coastal_ref
medflex_glm<-glm(cases~ tmean_mean+mean_prob_pred+tmean_mean*cluster+ offset(log((pop2010))), data = data.frame(monthly_human_anom[monthly_human_anom$year>2005 & monthly_human_anom$month %in% c(8:9), c("tmean_mean", "mean_prob_pred", "mean_abund", "cluster",  "pop2010", "lag_cases", "lag_mean_prob_pred", "cases", "ni_cases")]), family=negative.binomial(model.Y$theta))

expData <- neImpute(medflex_glm)

neMod1 <- neModel(cases~ (tmean_mean0+tmean_mean1)*cluster+ offset(log((pop2010))),  family =  negative.binomial(model.Y$theta), expData = expData, se="robust")


coef<-neMod1$neModelFit$coefficients
confint<-confint(neMod1)

#have to re-run the model with different reference level to derive confidence intervals for indirect effects in different regions
#central ref

monthly_human_anom_cent <- within(monthly_human_anom, cluster <- relevel(cluster, ref = "2"))

medflex_glm_cent<-glm(cases~ (tmean_mean+mean_prob_pred)+tmean_mean*cluster+ offset(log((pop2010))),  family = negative.binomial(model.Y$theta), data = data.frame(monthly_human_anom_cent[monthly_human_anom_cent $year>2005 & monthly_human_anom_cent $month %in% c(8:9), c("tmean_mean", "mean_prob_pred", "mean_abund", "cluster","cluster_fac",  "pop2010", "lag_cases", "lag_mean_prob_pred", "cases")]))

expData_cent <- neImpute(medflex_glm_cent)

neMod1_cent <- neModel(cases~ (tmean_mean0+tmean_mean1)*cluster+ offset(log((pop2010))),  family =  negative.binomial(.2), expData = expData_cent , se="robust")


coef_cent<-neMod1_cent$neModelFit$coefficients
confint_cent<-confint(neMod1_cent)


#inland ref

monthly_human_anom_inland <- within(monthly_human_anom, cluster <- relevel(cluster, ref = "3"))

medflex_glm_inland<-glm(cases~ (tmean_mean+mean_prob_pred)+tmean_mean*cluster+ offset(log((pop2010))),  family = negative.binomial(model.Y$theta), data = data.frame(monthly_human_anom_inland[monthly_human_anom_inland $year>2005 & monthly_human_anom_inland $month %in% c(8:9), c("tmean_mean", "mean_prob_pred", "mean_abund", "cluster","cluster_fac",  "pop2010", "lag_cases", "lag_mean_prob_pred", "cases", "ni_cases")]))

expData_inland <- neImpute(medflex_glm_inland)

neMod1_inland <- neModel(cases~ (tmean_mean0+tmean_mean1)*cluster+ offset(log((pop2010))),  family =  negative.binomial(.2), expData = expData_inland , se="robust")


coef_inland<-neMod1_inland$neModelFit$coefficients
confint_inland<-confint(neMod1_inland)


#estimating the change in incidence rate ratio for each zone due to mediation effects
coastal_RR<-exp(coef["tmean_mean1"])
coastal_RR_lci<-exp(confint["tmean_mean1",1])
coastal_RR_uci<-exp(confint["tmean_mean1",2])


central_RR<-exp(coef_cent["tmean_mean1"])
central_RR_lci<-exp(confint_cent["tmean_mean1",1])
central_RR_uci<-exp(confint_cent["tmean_mean1",2])


inland_RR<-exp(coef_inland["tmean_mean1"])
inland_RR_lci<-exp(confint_inland["tmean_mean1",1])
inland_RR_uci<-exp(confint_inland["tmean_mean1",2])





#mediation analysis for july and october
theta_glm<-glm.nb(cases~ (tmean_mean+mean_prob_pred)+tmean_mean*cluster+ offset(log((pop2010))), data = data.frame(monthly_human_anom[monthly_human_anom$year>2005 & monthly_human_anom$month %in% c(7,10), c("tmean_mean", "mean_prob_pred", "mean_abund", "cluster",  "pop2010", "lag_cases", "lag_mean_prob_pred", "cases")]))

#coastal_ref
medflex_glm_JO<-glm(cases~ (tmean_mean+mean_prob_pred)+tmean_mean*cluster+ offset(log((pop2010))),  family = negative.binomial(theta_glm$theta), data = data.frame(monthly_human_anom[monthly_human_anom$year>2005 & monthly_human_anom$month %in% c(7,10), c("tmean_mean", "mean_prob_pred", "mean_abund", "cluster",  "pop2010", "lag_cases", "lag_mean_prob_pred", "cases", "ni_cases")]))

expData_JO <- neImpute(medflex_glm_JO)

neMod1_JO <- neModel(cases~ (tmean_mean0+tmean_mean1)*factor(cluster)+ offset(log((pop2010))),  family =  negative.binomial(.2), expData = expData_JO, se="robust")



coef_JO<-neMod1_JO$neModelFit$coefficients
confint_JO<-confint(neMod1_JO)

#have to re-run the model with different reference level to derive confidence intervals for indirect effects in different regions
#central ref

monthly_human_anom_cent_JO <- within(monthly_human_anom, cluster <- relevel(cluster, ref = "2"))

medflex_glm_JO_cent<-glm(cases~ (tmean_mean+mean_prob_pred)+tmean_mean*cluster+ offset(log((pop2010))),  family = negative.binomial(theta_glm$theta), data = data.frame(monthly_human_anom_cent_JO[monthly_human_anom_cent_JO $year>2005 & monthly_human_anom_cent_JO $month %in% c(7,10), c("tmean_mean", "mean_prob_pred", "mean_abund", "cluster","cluster_fac",  "pop2010", "lag_cases", "lag_mean_prob_pred", "cases")]))

expData_cent_JO <- neImpute(medflex_glm_JO_cent)

neMod1_cent_JO <- neModel(cases~ (tmean_mean0+tmean_mean1)*cluster+ offset(log((pop2010))),  family =  negative.binomial(.2), expData = expData_cent_JO , se="robust")

coef_cent_JO<-neMod1_cent_JO$neModelFit$coefficients
confint_cent_JO<-confint(neMod1_cent_JO)


#inland ref


monthly_human_anom_inland_JO <- within(monthly_human_anom, cluster <- relevel(cluster, ref = "3"))

medflex_glm_JO_inland<-glm(cases~ (tmean_mean+mean_prob_pred)+tmean_mean*cluster+ offset(log((pop2010))),  family = negative.binomial(theta_glm$theta), data = data.frame(monthly_human_anom_inland_JO[monthly_human_anom_inland_JO $year>2005 & monthly_human_anom_inland_JO $month %in% c(7,10), c("tmean_mean", "mean_prob_pred", "mean_abund", "cluster","cluster_fac",  "pop2010", "lag_cases", "lag_mean_prob_pred", "cases")]))

expData_inland_JO <- neImpute(medflex_glm_JO_inland)

neMod1_inland_JO <- neModel(cases~ (tmean_mean0+tmean_mean1)*cluster+ offset(log((pop2010))),  family =  negative.binomial(.2), expData = expData_inland_JO , se="robust")

coef_inland_JO<-neMod1_inland_JO$neModelFit$coefficients
confint_inland_JO<-confint(neMod1_inland_JO)


#estimating the change in incidence rate ratio for each zone due to mediation effects
coastal_RR_JO<-exp(coef_JO["tmean_mean1"])
coastal_RR_JO_lci<-exp(confint_JO["tmean_mean1",1])
coastal_RR_JO_uci<-exp(confint_JO["tmean_mean1",2])


central_RR_JO<-exp(coef_cent_JO["tmean_mean1"])
central_RR_JO_lci<-exp(confint_cent_JO["tmean_mean1",1])
central_RR_JO_uci<-exp(confint_cent_JO["tmean_mean1",2])


inland_RR_JO<-exp(coef_inland_JO["tmean_mean1"])
inland_RR_JO_lci<-exp(confint_inland_JO["tmean_mean1",1])
inland_RR_JO_uci<-exp(confint_inland_JO["tmean_mean1",2])



#####Figure 4####

###base image for mediation diagram, from danslane's GitHub page
my_image<-readJPEG("data/MediationPlotter-master/MediationPlots.jpg")

specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))




#Extract terms for plot from mediation object#

a1_1<-specify_decimal(model.M[["coefficients"]][[ "tmean_mean"]]*100,2) #multiplying by 100 to make it infected mosquitoes per 100
b1_1<-specify_decimal((medflex_glm[["coefficients"]][[ "mean_prob_pred"]]),2)

c1_1<-specify_decimal(exp(coef["tmean_mean0"]),2)

a1_1_se<-summary(model.M)$coefficients[["tmean_mean","Std. Error"]]
b1_1_se<-summary(medflex_glm)$coefficients[["mean_prob_pred","Std. Error"]]




a1_1_lci<-specify_decimal(as.numeric(a1_1)-(1.96*(a1_1_se*100)),2)
a1_1_uci<-specify_decimal(as.numeric(a1_1)+(1.96*(a1_1_se*100)),2)
b1_1_lci<-specify_decimal((as.numeric(b1_1)-(1.96*b1_1_se))/100,2)
b1_1_uci<-specify_decimal((as.numeric(b1_1)+(1.96*b1_1_se))/100,2)
b1_1<-specify_decimal(as.numeric(b1_1)/100,2) #recoding it to scale to percentage after it's used for calculation

c1_1_lci<-specify_decimal(as.numeric(exp(confint["tmean_mean0",1])),2)
c1_1_uci<-specify_decimal(as.numeric(exp(confint["tmean_mean0",2])),2)


a1_2<-specify_decimal(model.M_cent[["coefficients"]][[ "tmean_mean"]]*100,2)
c1_2<-specify_decimal(exp(coef_cent["tmean_mean0"]),2)

a1_2_se<-summary(model.M_cent)$coefficients[["tmean_mean","Std. Error"]]



a1_2_lci<-specify_decimal(as.numeric(a1_2)-(1.96*(a1_1_se*100)),2)
a1_2_uci<-specify_decimal(as.numeric(a1_2)+(1.96*(a1_1_se*100)),2)
c1_2_lci<-specify_decimal(as.numeric(exp(confint_cent["tmean_mean0",1])),2)
c1_2_uci<-specify_decimal(as.numeric(exp(confint_cent["tmean_mean0",2])),2)



a1_3<-specify_decimal(model.M_inland[["coefficients"]][[ "tmean_mean"]]*100,2)
c1_3<-specify_decimal(exp(coef_inland["tmean_mean0"]),2)

a1_3_se<-summary(model.M_inland)$coefficients[["tmean_mean","Std. Error"]]



a1_3_lci<-specify_decimal(as.numeric(a1_3)-(1.96*(a1_1_se*100)),2)
a1_3_uci<-specify_decimal(as.numeric(a1_3)+(1.96*(a1_1_se*100)),2)
c1_3_lci<-specify_decimal(as.numeric(exp(confint_inland["tmean_mean0",1])),2)
c1_3_uci<-specify_decimal(as.numeric(exp(confint_inland["tmean_mean0",2])),2)


a<-paste("Marginal effect of T on pᵢ\n(linear scale)\n",paste("Coastal: ",a1_1," [",a1_1_lci,",",a1_1_uci,"]",sep=""),paste("Central: ",a1_2," [",a1_2_lci,",",a1_2_uci,"]",sep=""),paste("Inland: ",a1_3," [",a1_3_lci,",",a1_3_uci,"]",sep=""),sep="\n")


b<-paste("Marginal effect of pᵢ on Hᵢ\n(logarithmic scale)\n",paste(b1_1," [",b1_1_lci,",",b1_1_uci,"]",sep=""),sep="\n")

c<-paste("Direct Effect\nIRR associated with effect of T on Hᵢ\n",paste("Coastal: ",c1_1," [",c1_1_lci,",",c1_1_uci,"]a",sep=""),paste("Central: ",c1_2," [",c1_2_lci,",",c1_2_uci,"]a",sep=""),paste("Inland: ",c1_3," [",c1_3_lci,",",c1_3_uci,"]a",sep=""),sep="\n")



#plotting DAG for August and September
png(filename="Results/Fig7_DAG_aug_sep_7_25.png", width = 8.5,height=6.5,units = "in", res = 300)
par(mar = c(0,0,0,0),family="Arial") # set zero margins on all 4 sides
plot(x = NULL, y = NULL, xlim = c(0,1584), ylim = c(0,891), pch = '',
     xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', xaxs = 'i', yaxs = 'i',
     bty = 'n') # plot empty figure
rasterImage(my_image, xleft = 0, ybottom = 0, xright = 1584, ytop = 891) # plot jpeg

#Plot coefficients
text(1584/2,90, c,cex=1)
text(200,500,a,cex=1)
text(1330,500,b,cex=1)

#Plot variable names

text(240,185,"Monthly Mean\nTemperature\n(T, °C)",cex=1.2)
text(797,719,"Cx. Infection\nProbability\n(pᵢ, i/100)",cex=1.2)
text(1335,185,"Human WNV\nIncidence\n(Hᵢ, cases/100,000)",cex=1.2)

#center text
PE <- paste("Mediation Effect\nIRR of T mediated through pᵢ\n",
            paste("Coastal = ",specify_decimal(exp(coef["tmean_mean1"]),2)," [", specify_decimal(exp(confint["tmean_mean1",1]) ,2),",", specify_decimal(exp(confint["tmean_mean1",2]) ,2),"]a",sep=""),
            paste("Central = ",specify_decimal(exp(coef_cent["tmean_mean1"]),2)," [", specify_decimal(exp(confint_cent["tmean_mean1",1]) ,2),",", specify_decimal(exp(confint_cent["tmean_mean1",2]) ,2),"]ab",sep=""),
            paste("Inland = ",specify_decimal(exp(coef_inland["tmean_mean1"]),2)," [", specify_decimal(exp(confint_inland["tmean_mean1",1]) ,2),",", specify_decimal(exp(confint_inland["tmean_mean1",2]) ,2),"]b",sep=""),
            sep="\n")

text((1584/2),(891/2),PE,font=2, cex=1.4)

#Add main title
text(1584/2,845,paste("August-September",  paste("\nR²=",specify_decimal(rsq.sse(medflex_glm, adj = T), 2)),sep=""), font=2,cex=1.5 )

#Add letter
text(10,845,"a",pos=4,  cex=1.5, font=2)

dev.off()








#plotting dag for july , october

a1_1<-specify_decimal(model.M_JO[["coefficients"]][[ "tmean_mean"]]*100,2) #multiplying by 100 to make it infected mosquitoes per 100
b1_1<-specify_decimal((medflex_glm_JO[["coefficients"]][[ "mean_prob_pred"]]),2)
c1_1<-specify_decimal(exp(coef_JO["tmean_mean0"]),2)

a1_1_se<-summary(model.M_JO)$coefficients[["tmean_mean","Std. Error"]]
b1_1_se<-summary(medflex_glm_JO)$coefficients[["mean_prob_pred","Std. Error"]]



a1_1_lci<-specify_decimal(as.numeric(a1_1)-(1.96*(a1_1_se*100)),2)
a1_1_uci<-specify_decimal(as.numeric(a1_1)+(1.96*(a1_1_se*100)),2)
b1_1_lci<-specify_decimal((as.numeric(b1_1)-(1.96*b1_1_se))/100,2)
b1_1_uci<-specify_decimal((as.numeric(b1_1)+(1.96*b1_1_se))/100,2)
b1_1<-specify_decimal(as.numeric(b1_1)/100,2) #recoding it to scale to percent after it's used for calculation
c1_1_lci<-specify_decimal(as.numeric(exp(confint_JO["tmean_mean0",1])),2)
c1_1_uci<-specify_decimal(as.numeric(exp(confint_JO["tmean_mean0",2])),2)


a1_2<-specify_decimal(model.M_JO_cent[["coefficients"]][[ "tmean_mean"]]*100,2)
c1_2<-specify_decimal(exp(coef_cent_JO["tmean_mean0"]),2)

a1_2_se<-summary(model.M_JO_cent)$coefficients[["tmean_mean","Std. Error"]]




a1_2_lci<-specify_decimal(as.numeric(a1_2)-(1.96*(a1_1_se*100)),2)
a1_2_uci<-specify_decimal(as.numeric(a1_2)+(1.96*(a1_1_se*100)),2)
c1_2_lci<-specify_decimal(as.numeric(exp(confint_cent_JO["tmean_mean0",1])),2)
c1_2_uci<-specify_decimal(as.numeric(exp(confint_cent_JO["tmean_mean0",2])),2)



a1_3<-specify_decimal(model.M_JO_inland[["coefficients"]][[ "tmean_mean"]]*100,2)
c1_3<-specify_decimal(exp(coef_inland_JO["tmean_mean0"]),2)

a1_3_se<-summary(model.M_JO_inland)$coefficients[["tmean_mean","Std. Error"]]


a1_3_lci<-specify_decimal(as.numeric(a1_3)-(1.96*(a1_1_se*100)),2)
a1_3_uci<-specify_decimal(as.numeric(a1_3)+(1.96*(a1_1_se*100)),2)
c1_3_lci<-specify_decimal(as.numeric(exp(confint_inland_JO["tmean_mean0",1])),2)
c1_3_uci<-specify_decimal(as.numeric(exp(confint_inland_JO["tmean_mean0",2])),2)




a<-paste("Marginal effect of T on pᵢ\n(linear scale)\n",paste("Coastal: ",a1_1," [",a1_1_lci,",",a1_1_uci,"]",sep=""),paste("Central: ",a1_2," [",a1_2_lci,",",a1_2_uci,"]",sep=""),paste("Inland: ",a1_3," [",a1_3_lci,",",a1_3_uci,"]",sep=""),sep="\n")


b<-paste("Marginal effect of pᵢ on Hᵢ\n(logarithmic scale)\n",paste(b1_1," [",b1_1_lci,",",b1_1_uci,"]",sep=""),sep="\n")

c<-paste("Direct Effect\nIRR associated with effect of T on Hᵢ\n",paste("Coastal: ",c1_1," [",c1_1_lci,",",c1_1_uci,"]a",sep=""),paste("Central: ",c1_2," [",c1_2_lci,",",c1_2_uci,"]a",sep=""),paste("Inland: ",c1_3," [",c1_3_lci,",",c1_3_uci,"]a",sep=""),sep="\n")




#plotting july october DAG
png(filename="Results/Fig7_DAG_july_oct_7_25.png", width = 8.5,height=6.5,units = "in", res = 300)
par(mar = c(0,0,0,0),family="Arial") # set zero margins on all 4 sides
plot(x = NULL, y = NULL, xlim = c(0,1584), ylim = c(0,891), pch = '',
     xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', xaxs = 'i', yaxs = 'i',
     bty = 'n') # plot empty figure
rasterImage(my_image, xleft = 0, ybottom = 0, xright = 1584, ytop = 891) # plot jpeg


#Plot coefficients
text(1584/2,90, c,cex=1)
text(200,500,a,cex=1)
text(1330,500,b,cex=1)

#Plot variable names
text(240,185,"Monthly Mean\nTemperature\n(T, °C)",cex=1.2)
text(797,719,"Cx. Infection\nProbability\n(pᵢ, i/100)",cex=1.2)
text(1335,185,"Human WNV\nIncidence\n(Hᵢ, cases/100,000)",cex=1.2)


#center text
PE <- paste("Mediation Effect\nIRR of T mediated through pᵢ\n",
            paste("Coastal = ",specify_decimal(exp(coef_JO["tmean_mean1"]),2)," [", specify_decimal(exp(confint_JO["tmean_mean1",1]) ,2),",", specify_decimal(exp(confint_JO["tmean_mean1",2]) ,2),"]a",sep=""),
            paste("Central = ",specify_decimal(exp(coef_cent_JO["tmean_mean1"]),2)," [", specify_decimal(exp(confint_cent_JO["tmean_mean1",1]) ,2),",", specify_decimal(exp(confint_cent_JO["tmean_mean1",2]) ,2),"]a",sep=""),
            paste("Inland = ",specify_decimal(exp(coef_inland_JO["tmean_mean1"]),2)," [", specify_decimal(exp(confint_inland_JO["tmean_mean1",1]) ,2),",", specify_decimal(exp(confint_inland_JO["tmean_mean1",2]) ,2),"]a",sep=""),sep="\n")

text((1584/2),(891/2),PE,font=2, cex=1.4)

#Add main title
text(1584/2,845,paste("July/October",  paste("\nR²=",specify_decimal(rsq.sse(medflex_glm_JO, adj = T), 2)),sep=""), font=2,cex=1.5 )

#Add letter
text(10,845,"b",pos=4,  cex=1.5, font=2)

dev.off()





####Figure 5####
CT_pop<-read.csv("data/ct_pop.csv", header=T)
CT_pop<-unique(CT_pop[,c("GEOID", "tract_pop2010")])


human_CT<-read.csv("data/WNV_census_tract_studyarea_4_12_18.csv", header=T)

human_CT_la<-human_CT[human_CT$study_area=="la_metro",]
human_CT_la$GEOID<-as.character(human_CT_la$GEOID)
human_CT_la<-human_CT_la[human_CT_la$agent %in% c("WNF", "WNND"),]

human_CT_la<-human_CT_la[!is.na(as.Date(human_CT_la$onsetdt, format="%m/%d/%Y")) | !is.na(as.Date(human_CT_la$hospadmitdt, format="%m/%d/%Y")),]

human_CT_la$date<-as.Date(human_CT_la$onsetdt, format="%m/%d/%Y")

#converting the 11 cases with only a hospital admit date to an onset date by subtracting 4 days, which is the mean difference between them.
human_CT_la$date[is.na(human_CT_la$date)]<-as.Date(human_CT_la$hospadmitdt[is.na(human_CT_la$date)], format="%m/%d/%Y")-4


CT<-load_rgdal_from_googledrive("LA_censustract_shapefile")
CT_la<-crop(CT, extent(la_metro))

##jittering census tracts with 20km random displacement
CT_la@data$latitude<-jitter(coordinates(gCentroid(CT_la, byid = T))[,2], amount=.02)
CT_la@data$longitude<-jitter(coordinates(gCentroid(CT_la, byid = T))[,1], amount=.02)

CT_la@data$GEO_ID<-gsub("1400000US0", "",CT_la@data$GEO_ID )
CT_human<-full_join(human_CT_la, CT_la@data, by=c("GEOID"="GEO_ID"))

CT_human$GEOID<-as.numeric(CT_human$GEOID)
CT_human<-left_join(CT_human, CT_pop, by="GEOID")
CT_human$tract_pop2010<-as.numeric(CT_human$tract_pop2010)


#using TopoWx climate data

#all months
#subsetting to after 2006 and for july_sept
topowx_data_sum8_9<-list()
year_indx<-c(12,15,29,32)
for (i in 1:length(year_indx)){
  topowx_data_sum8_9[[i]]<- topowx_data[[year_indx[i]]][[which(month(as.Date(gsub("\\.","-",gsub("X","",topowx_data[[1]]@data@names)))) %in% c(9))]]
  
}
mean_2011_8_9<-mask(mean(stack(topowx_data_sum8_9[c(1,3)])), la_metro)
mean_2014_8_9<-mask(mean(stack(topowx_data_sum8_9[c(2,4)])), la_metro)




mean_2011_8_9_dff<-mean_2011_8_9
mean_2014_8_9_dff<-mean_2014_8_9


#coding different transitions
#no inhibitory to inhibitory found
mean_2011_8_9_dff[mean_2011_8_9<21& mean_2014_8_9<21]<-0
mean_2011_8_9_dff[mean_2011_8_9<21& mean_2014_8_9>21 &mean_2014_8_9<22.7 ]<-1
mean_2011_8_9_dff[mean_2011_8_9>=21& mean_2011_8_9<=22.7 &mean_2014_8_9>=22.7 ]<-2
mean_2011_8_9_dff[mean_2011_8_9<21& mean_2014_8_9>22.7]<-3
mean_2011_8_9_dff[mean_2011_8_9>22.7& mean_2014_8_9>22.7]<-4


#calculating the number of human cases within each climate scenario for 2011 and 2014
inhib_trans_cases<-length(
  over(SpatialPointsDataFrame(coords = CT_human[CT_human$year %in% c(2011)& month(CT_human$date) %in% c(9) ,c("longitude", "latitude")],proj4string=CRS(proj4string(mean_2011_8_9_dff)), data = CT_human[CT_human$year %in% c(2011)& month(CT_human$date) %in% c(9) ,]),
       change_2011_2014_8_9[change_2011_2014_8_9$layer=="Inhibitory→Transitional",])[which(!is.na(over(SpatialPointsDataFrame(coords = CT_human[CT_human$year %in% c(2011)& month(CT_human$date) %in% c(9) ,c("longitude", "latitude")],proj4string=CRS(proj4string(mean_2011_8_9_dff)), data = CT_human[CT_human$year %in% c(2011)& month(CT_human$date) %in% c(9) ,]),change_2011_2014_8_9[change_2011_2014_8_9$layer=="Inhibitory→Transitional",]))),])

trans_fav_cases<-length(
  over(SpatialPointsDataFrame(coords = CT_human[CT_human$year %in% c(2011)& month(CT_human$date) %in% c(9) ,c("longitude", "latitude")],proj4string=CRS(proj4string(mean_2011_8_9_dff)), data = CT_human[CT_human$year %in% c(2011)& month(CT_human$date) %in% c(9) ,]),
       change_2011_2014_8_9[change_2011_2014_8_9$layer=="Transitional→Favorable",])[which(!is.na(over(SpatialPointsDataFrame(coords = CT_human[CT_human$year %in% c(2011)& month(CT_human$date) %in% c(9) ,c("longitude", "latitude")],proj4string=CRS(proj4string(mean_2011_8_9_dff)), data = CT_human[CT_human$year %in% c(2011)& month(CT_human$date) %in% c(9) ,]),change_2011_2014_8_9[change_2011_2014_8_9$layer=="Transitional→Favorable",]))),])


inhib_fav_cases<-length(
  over(SpatialPointsDataFrame(coords = CT_human[CT_human$year %in% c(2011)& month(CT_human$date) %in% c(9) ,c("longitude", "latitude")],proj4string=CRS(proj4string(mean_2011_8_9_dff)), data = CT_human[CT_human$year %in% c(2011)& month(CT_human$date) %in% c(9) ,]),
       change_2011_2014_8_9[change_2011_2014_8_9$layer=="Inhibitory→Favorable",])[which(!is.na(over(SpatialPointsDataFrame(coords = CT_human[CT_human$year %in% c(2011)& month(CT_human$date) %in% c(9) ,c("longitude", "latitude")],proj4string=CRS(proj4string(mean_2011_8_9_dff)), data = CT_human[CT_human$year %in% c(2011)& month(CT_human$date) %in% c(9) ,]),change_2011_2014_8_9[change_2011_2014_8_9$layer=="Inhibitory→Favorable",]))),])

fav_fav_cases<-length(
  over(SpatialPointsDataFrame(coords = CT_human[CT_human$year %in% c(2011)& month(CT_human$date) %in% c(9) ,c("longitude", "latitude")],proj4string=CRS(proj4string(mean_2011_8_9_dff)), data = CT_human[CT_human$year %in% c(2011)& month(CT_human$date) %in% c(9) ,]),
       change_2011_2014_8_9[change_2011_2014_8_9$layer=="Favorable→Favorable",])[which(!is.na(over(SpatialPointsDataFrame(coords = CT_human[CT_human$year %in% c(2011)& month(CT_human$date) %in% c(9) ,c("longitude", "latitude")],proj4string=CRS(proj4string(mean_2011_8_9_dff)), data = CT_human[CT_human$year %in% c(2011)& month(CT_human$date) %in% c(9) ,]),change_2011_2014_8_9[change_2011_2014_8_9$layer=="Favorable→Favorable",]))),])


inhib_trans_cases_2014<-length(
  over(SpatialPointsDataFrame(coords = CT_human[CT_human$year %in% c(2014)& month(CT_human$date) %in% c(9) ,c("longitude", "latitude")],proj4string=CRS(proj4string(mean_2014_8_9_dff)), data = CT_human[CT_human$year %in% c(2014)& month(CT_human$date) %in% c(9) ,]),
       change_2011_2014_8_9[change_2011_2014_8_9$layer=="Inhibitory→Transitional",])[which(!is.na(over(SpatialPointsDataFrame(coords = CT_human[CT_human$year %in% c(2014)& month(CT_human$date) %in% c(9) ,c("longitude", "latitude")],proj4string=CRS(proj4string(mean_2014_8_9_dff)), data = CT_human[CT_human$year %in% c(2014)& month(CT_human$date) %in% c(9) ,]),change_2011_2014_8_9[change_2011_2014_8_9$layer=="Inhibitory→Transitional",]))),])

trans_fav_cases_2014<-length(
  over(SpatialPointsDataFrame(coords = CT_human[CT_human$year %in% c(2014)& month(CT_human$date) %in% c(9) ,c("longitude", "latitude")],proj4string=CRS(proj4string(mean_2014_8_9_dff)), data = CT_human[CT_human$year %in% c(2014)& month(CT_human$date) %in% c(9) ,]),
       change_2011_2014_8_9[change_2011_2014_8_9$layer=="Transitional→Favorable",])[which(!is.na(over(SpatialPointsDataFrame(coords = CT_human[CT_human$year %in% c(2014)& month(CT_human$date) %in% c(9) ,c("longitude", "latitude")],proj4string=CRS(proj4string(mean_2014_8_9_dff)), data = CT_human[CT_human$year %in% c(2014)& month(CT_human$date) %in% c(9) ,]),change_2011_2014_8_9[change_2011_2014_8_9$layer=="Transitional→Favorable",]))),])


inhib_fav_cases_2014<-length(
  over(SpatialPointsDataFrame(coords = CT_human[CT_human$year %in% c(2014)& month(CT_human$date) %in% c(9) ,c("longitude", "latitude")],proj4string=CRS(proj4string(mean_2014_8_9_dff)), data = CT_human[CT_human$year %in% c(2014)& month(CT_human$date) %in% c(9) ,]),
       change_2011_2014_8_9[change_2011_2014_8_9$layer=="Inhibitory→Favorable",])[which(!is.na(over(SpatialPointsDataFrame(coords = CT_human[CT_human$year %in% c(2014)& month(CT_human$date) %in% c(9) ,c("longitude", "latitude")],proj4string=CRS(proj4string(mean_2014_8_9_dff)), data = CT_human[CT_human$year %in% c(2014)& month(CT_human$date) %in% c(9) ,]),change_2011_2014_8_9[change_2011_2014_8_9$layer=="Inhibitory→Favorable",]))),])

fav_fav_cases_2014<-length(
  over(SpatialPointsDataFrame(coords = CT_human[CT_human$year %in% c(2014)& month(CT_human$date) %in% c(9) ,c("longitude", "latitude")],proj4string=CRS(proj4string(mean_2014_8_9_dff)), data = CT_human[CT_human$year %in% c(2014)& month(CT_human$date) %in% c(9) ,]),
       change_2011_2014_8_9[change_2011_2014_8_9$layer=="Favorable→Favorable",])[which(!is.na(over(SpatialPointsDataFrame(coords = CT_human[CT_human$year %in% c(2014)& month(CT_human$date) %in% c(9) ,c("longitude", "latitude")],proj4string=CRS(proj4string(mean_2014_8_9_dff)), data = CT_human[CT_human$year %in% c(2014)& month(CT_human$date) %in% c(9) ,]),change_2011_2014_8_9[change_2011_2014_8_9$layer=="Favorable→Favorable",]))),])


mean_2011_8_9_dff_0<-as.factor(mean_2011_8_9_dff)

change_2011_2014_8_9<-rasterToPolygons(mean_2011_8_9_dff_0,  n=4, na.rm=TRUE, digits=12, dissolve=TRUE)
change_2011_2014_8_9$layer<-factor(change_2011_2014_8_9$layer,levels =c("1","2","3","4")  ,labels=c("Inhibitory→Transitional","Transitional→Favorable","Inhibitory→Favorable", "Favorable→Favorable"))

mapviewOptions(legend.pos="bottomright")



map_2011_2014_8_9_rand<-mapview(change_2011_2014_8_9, map.types="CartoDB.Positron", trim=T,zcol="layer" ,col.regions=c("green","purple", "blue","transparent" ), alpha.regions=.4, legend=TRUE,layer.name="Temperature Shift (2011→2014)", lwd=.00001 )+
  mapview(la_metro,  alpha.regions=.001, map.types="CartoDB.Positron", legend=FALSE) + 
  mapview(spsample(change_2011_2014_8_9[change_2011_2014_8_9$layer=="Transitional→Favorable",], n=trans_fav_cases, type="clustered"), cex=4, lwd=1.3, map.types="CartoDB.Positron",legend = FALSE, alpha.regions=.8, layer.name="Human WNV Cases", col.regions="red")+ 
  mapview(spsample(change_2011_2014_8_9[change_2011_2014_8_9$layer=="Favorable→Favorable",], n=fav_fav_cases, type="clustered"), cex=4, lwd=1.3, map.types="CartoDB.Positron",legend = FALSE, alpha.regions=.8, layer.name="Human WNV Cases", col.regions="red")+
  mapview(spsample(change_2011_2014_8_9[change_2011_2014_8_9$layer=="Inhibitory→Transitional",], n=inhib_trans_cases_2014, type="clustered"), cex=4, lwd=1.3, map.types="CartoDB.Positron",legend = FALSE, alpha.regions=.8, layer.name="Human WNV Cases", col.regions="yellow") + 
  mapview(spsample(change_2011_2014_8_9[change_2011_2014_8_9$layer=="Transitional→Favorable",], n=trans_fav_cases_2014, type="clustered"), cex=4, lwd=1.3, map.types="CartoDB.Positron",legend = FALSE, alpha.regions=.8, layer.name="Human WNV Cases", col.regions="yellow")+ 
  mapview(spsample(change_2011_2014_8_9[change_2011_2014_8_9$layer=="Inhibitory→Favorable",], n=inhib_fav_cases_2014, type="clustered"), cex=4, lwd=1.3, map.types="CartoDB.Positron",legend = FALSE, alpha.regions=.8, layer.name="Human WNV Cases", col.regions="yellow")+ 
  mapview(spsample(change_2011_2014_8_9[change_2011_2014_8_9$layer=="Favorable→Favorable",], n=fav_fav_cases_2014, type="clustered"), cex=4, lwd=1.3, map.types="CartoDB.Positron",legend = FALSE, alpha.regions=.8, layer.name="Human WNV Cases", col.regions="yellow")


