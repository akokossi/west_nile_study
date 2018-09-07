

#===============STARTING DATA PREPARATION FOR RANDOM FOREST MODEL====================
#' @description run_rf_modelrun random forest model
#' @param df - data frame containing training and validation sets
#' @export
run_rf_model<-function(df)
{

  data1<-df
  # Rescaling the data for modeling:

  data1$LATITUDE<-(data1$LATITUDE-mean(data1$LATITUDE))/sd(data1$LATITUDE)
  data1$LONGITUDE<-(data1$LONGITUDE-mean(data1$LONGITUDE))/sd(data1$LONGITUDE)
  data1$seasonindex<-(data1$seasonindex-mean(data1$seasonindex))/sd(data1$seasonindex)


  ### Transform categorical variables into binaries:


  dims<-as.character(unique(as.character(data1$SPECIES)))
  k1<-1:length(dims)
  specieslist<-paste("SPECIES_",k1,sep="")
  specformula<-as.formula(paste("RESULT ~seasonindex+LATITUDE+LONGITUDE+",paste(specieslist,collapse="+",sep=""),sep=""))
  dd<-data.frame(varname=specieslist)
  for (i in 1:length(dims)){
    data1[,paste('SPECIES','_', i,sep="")]<-ifelse(as.character(data1[,"SPECIES"])==as.character(dims[i]),1,0)
    dd[i,"vardesc"]<-as.character(dims[i])
  }


  # defining holdout
  data1$holdout<-as.numeric(ifelse(as.numeric(data1$WEEK+data1$SEASON_YEAR*100)>=201824,1,0))

  data1$BLOCK<-as.factor(data1$BLOCK)
  data1$RESULT<-as.factor(data1$RESULT)
  data1$WEEK<-as.factor(data1$WEEK)
  data1$SPECIES<-as.factor(data1$SPECIES)

  df.train <- data1[data1$holdout==0 && data1$t_nfreq>data1$nfreq & data1$t_nfreq>1000,]
  df.test <- data1[data1$holdout==1 ,]

  df.train$BLOCK<-as.factor(df.train$BLOCK)
  df.train$RESULT<-as.factor(df.train$RESULT)
  df.train$WEEK<-as.factor(df.train$WEEK)
  df.train$SPECIES<-as.factor(df.train$SPECIES)

  df.test$BLOCK<-as.factor(df.test$BLOCK)
  df.test$RESULT<-as.factor(df.test$RESULT)
  df.test$WEEK<-as.factor(df.test$WEEK)
  df.test$SPECIES<-as.factor(df.test$SPECIES)


  # ==============RANDOM FOREST Training function

  rf<-randomForest( specformula, data = df.train, ntree =length(unique(df.train$BLOCK)),weights=t_nfreq)
  return(rf)
}
#=========================================
#' @description plot_my_time_series plot time series
#' @param df - data frame containing and validation training sets
#' @export
plot_my_time_series<-function(df)
{
  data1<-df

  mapdata<-sqldf("select BLOCK,LONGITUDE,LATITUDE,avg(NUMBER_OF_MOSQUITOES) as NMOSQ from data1
                 group by BLOCK,LONGITUDE,LATITUDE")

  ### MISSING VALUES HANDLING IN SEASONALITY ANALYSIS: set number of mosquitoes to zero for cases where they are missing
  seasonalitydata<-as.data.frame(sqldf("
                                       select RESULT,
                                       YEAR_WEEK,
                                       SUM(ifnull(NUMBER_OF_MOSQUITOES,0)) as NMOSQ from data1
                                       group by RESULT,YEAR_WEEK"))
  reshaped_seas<-reshape(seasonalitydata, idvar = "YEAR_WEEK", timevar = "RESULT", direction = "wide")
  reshaped_seas[is.na(reshaped_seas)]<-0
  reshaped_seas$positive_PCT<- round(100*reshaped_seas$NMOSQ.positive/(reshaped_seas$NMOSQ.positive+reshaped_seas$NMOSQ.negative),4)

  ggplot(reshaped_seas, aes(YEAR_WEEK, positive_PCT)) + geom_line() + xlab("") + ylab("weekly Views")
}


#=========================================
#' @description run_mixed_model run mixed effect model on data set
#' @param df - data frame containing training set
#' @export
run_mixed_model<-function(df)
{
  # convert category variables into factors in preparation for regression
  data1<-df
  #data1<-data1saved

  # Rescaling the data for modeling:

  data1$LATITUDE<-(data1$LATITUDE-mean(data1$LATITUDE))/sd(data1$LATITUDE)
  data1$LONGITUDE<-(data1$LONGITUDE-mean(data1$LONGITUDE))/sd(data1$LONGITUDE)
  data1$seasonindex<-(data1$seasonindex-mean(data1$seasonindex))/sd(data1$seasonindex)

  # defining holdout
  data1$holdout<-as.numeric(ifelse(as.numeric(data1$WEEK+data1$SEASON_YEAR*100)>=201824,1,0))

  data1$BLOCK<-as.factor(data1$BLOCK)
  data1$RESULT<-as.factor(data1$RESULT)
  data1$WEEK<-as.factor(data1$WEEK)
  data1$SPECIES<-as.factor(data1$SPECIES)



  df.train <- data1[data1$holdout==0 && data1$t_nfreq>data1$nfreq,]              #get training set
  df.test <- data1[data1$holdout==1,]

  df.train$BLOCK<-as.factor(df.train$BLOCK)
  df.train$RESULT<-as.factor(df.train$RESULT)
  df.train$WEEK<-as.factor(df.train$WEEK)
  df.train$SPECIES<-as.factor(df.train$SPECIES)

  df.test$BLOCK<-as.factor(df.test$BLOCK)
  df.test$RESULT<-as.factor(df.test$RESULT)
  df.test$WEEK<-as.factor(df.test$WEEK)
  df.test$SPECIES<-as.factor(df.test$SPECIES)


  m <- glmer(RESULT ~seasonindex +LONGITUDE+LATITUDE+SPECIES+(1|SPECIES) , data = df.train, family =
               binomial(link=logit),
             control=glmerControl(optimizer = "Nelder_Mead", tolPwrss = 1e-10,optCtrl=list(maxfun=2e5) ),
             weights=t_nfreq)



  # print the mod results without correlations among fixed effects
  print(m, corr = FALSE)
  return(m)
}
################################################################################################

#===============================================
#' @description create Google Map for the Chicago land area and add locations studied
#' @param df - data frame containing the 2 columns of interest
#' @export
googl_map<-function(df)
{
  data1<-df
  mapdata<-sqldf("select BLOCK,LONGITUDE,LATITUDE,avg(NUMBER_OF_MOSQUITOES) as NMOSQ from data1
                 group by BLOCK,LONGITUDE,LATITUDE")

  counties <- map_data("county")
  il_county <- subset(counties, region == "illinois")

  states <- map_data("state")
  il_df <- subset(states, region == "illinois")


  mapdata2<-mapdata
  names(mapdata2)<-tolower(names(mapdata2))

  il_county2<-il_county

  il_county2$longitude<- round(as.numeric(il_county2$long),5)
  il_county2$latitude<- round(as.numeric(il_county2$lat),5)

  mapdata2$longitude<-  round(as.numeric(mapdata2$longitude),5)
  mapdata2$latitude<-  round(as.numeric(mapdata2$latitude),5)

  myplotcords<-mapdata2[,c("longitude","latitude","nmosq","block")]
  towns_fortify <- fortify(myplotcords, region="BLOCK")


  sites <-myplotcords
  # State boundaries from the maps package. The fill option must be TRUE.
  bounds <- map('state', c('illinois'), fill=TRUE, plot=FALSE)
  # A custom icon.
  icons <- awesomeIcons(
    icon = 'disc',
    iconColor = 'black',
    library = 'ion', # Options are 'glyphicon', 'fa', 'ion'.
    markerColor = 'blue',
    squareMarker = TRUE
  )

  map <- leaflet(data = sites) %>%

    addProviderTiles("CartoDB.Positron", group = "Map") %>%
    addProviderTiles("Esri.WorldImagery", group = "Satellite") %>%
    addProviderTiles("Esri.WorldShadedRelief", group = "Relief") %>%
    addMarkers(~longitude, ~latitude, label = ~block, group = "Sites") %>%
    addPolygons(data=bounds, group="States", weight=2, fillOpacity = 0) %>%
    addScaleBar(position = "bottomleft") %>%
    addLayersControl(
      baseGroups = c("Map", "Satellite", "Relief"),
      overlayGroups = c("Sites", "States"),
      options = layersControlOptions(collapsed = FALSE)
    )
  invisible(print(map))
}
###### UTILITY FUNCTION DEFINITIONS
#' @description create cross tabulation given to factor columns
#' @param df - data frame containing the 2 columns of interest
#' @param col1 - first column in cross tab
#' @param col2 - second column in cross tab
#' @export
view_freq<-function(df,col1,col2)
{
  return(table(df[,col1],df[,col2]))
}
