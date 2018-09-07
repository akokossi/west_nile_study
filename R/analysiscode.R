
#' @import lme4
#' @import sqldf
#' @import randomForest
#' @import ggplot2
#' @import stringr
#' @import dplyr
#' @import ggmap
#' @import maps
#' @import mapdata
#' @import leaflet
#' @import maps
#' @import htmlwidgets
#' @import caret
source("file:///C:/Users/penof/OneDrive/Documents/myprojects/RProjects/westnileproject/R/mufuncs.R")
##### GLOBAL PARAMETERS
data_path<<-"file:///C:/Users/penof/OneDrive/Documents/my reading/Interviews/ATT/new/west-nile-virus-wnv-mosquito-test-results.csv"


#==============================================


data1<-read.csv(data_path,check.names=FALSE)
names(data1)<-gsub(" ","_",names(data1))
data1$YEAR_WEEK<- as.numeric(data1$SEASON_YEAR)*100+as.numeric(data1$WEEK)

freqdata<-sqldf("select a1.BLOCK,a1.RESULT, nfreq,t_nfreq from
                (select BLOCK,RESULT,SUM(NUMBER_OF_MOSQUITOES) as nfreq from data1
                group by BLOCK,RESULT) a1
                inner join (select BLOCK,SUM(NUMBER_OF_MOSQUITOES) as t_nfreq from data1
                group by BLOCK) a2 on a1.BLOCK=a2.BLOCK")

# Manual imputation of longitute and latitude using neighborhood average or using O'Aire airport coordinates
# since most of the missing coordinates is O'Hare airport
imputX<-sqldf("select BLOCK,AVG(case when LONGITUDE is null then  -87.836723 else LONGITUDE   end) as avgLong,AVG(case when LATITUDE is null then  41.977226 else LATITUDE end ) as avglat from data1
              group by BLOCK " )

x1<-sqldf("select YEAR_WEEK,AVG(NUMBER_OF_MOSQUITOES) as avgfreq from data1
          group by YEAR_WEEK " )
aggx1<-mean(x1$avgfreq)

# Creating seasonality lift data from yearly by week data on number of mosquitoes
# longitude and latitude are used in the model
data1<-sqldf(paste("select
                   a1.SEASON_YEAR,
                   a1.WEEK, a1.YEAR_WEEK,
                   TEST_ID,
                   a1.BLOCK,
                   TRAP,
                   TRAP_TYPE,
                   TEST_DATE,
                   NUMBER_OF_MOSQUITOES,
                   a1.RESULT,
                   SPECIES,
                   (case when LATITUDE is null then avglat else LATITUDE end) as LATITUDE,
                   (case when LONGITUDE is null then avglong else LONGITUDE end) as LONGITUDE,
                   nfreq,t_nfreq,avgfreq/",aggx1," as seasonindex from data1 a1
                   inner join freqdata a2 on a1.BLOCK=a2.BLOCK and a1.RESULT=a2.RESULT
                   inner join x1 a3 on a1.YEAR_WEEK=a3.YEAR_WEEK
                   inner join imputX a4 on a1.BLOCK=a4.BLOCK
                   ",sep=""))

data1saved<-data1


check_missing<-as.data.frame(sapply(data1, function(x) sum(is.na(x) | is.null(x) | as.character(x)=="" )))

names(check_missing)<-"missing_values"
check_missing

#===================================


### MISSING VALUES HANDLING: set number of mosquitoes to zero for cases where they are missing
cross_by_species<-view_freq(data1,"SPECIES","RESULT")
view_freq(data1,"SPECIES","RESULT")

#==== Google Map
googl_map(data1)
#
#==================PREDICT Generalized Mixed Model FULL DATASET===============================
data1<-data1saved
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


m<-run_mixed_model(data1)
data1$probability <- predict(m, data1, re.form=NA, type="response")
data1$linfitted <- predict(m, data1, re.form=NA, type="link")
data1$pred<-ifelse(data1$probability>0.00211,"positive","negative")   #0.00276
table(data1$RESULT,data1$pred  )


#====================PREDICT Generalized Mixed Model Validation DATASET======

df.test$probability <- predict(m, df.test, re.form=NA, type="response")
df.test$linfitted <- predict(m, df.test, re.form=NA, type="link")
df.test$pred<-ifelse(df.test$probability>0.00211,"positive","negative")   #0.00276
table(df.test$RESULT,df.test$pred  )


# ==============Random forest prediction and confusion matrix for the full dataset
data1<-data1saved
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



rf<-run_rf_model(data1)
data1$predicted.response <- predict(rf , data1)

table(data1$predicted.response,data1$RESULT)

#============FULL prediction and   Confusion matrix validation========
df.test$predicted.response <- predict(rf , df.test)

table(df.test$predicted.response,df.test$RESULT)

#==============VARIABLE IMPORTANCE====

importance(rf)
varImp(rf)



