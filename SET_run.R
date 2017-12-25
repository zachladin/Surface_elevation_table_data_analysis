########################################################################################################
##SET Data Analyses - Run Script
#Zachary Ladin (zach@udel.edu)

########################################################################################################
#Clear the Global Environment.
rm(list=ls())

########################################################################################################
#Set the working directory.
setwd("~/Dropbox (ZachTeam)/ZachGreg (1)/SET/SET_analysis_final_2")

########################################################################################################
#If needed, install required packages.
# install.packages("knitr")
# install.packages("stringr")
# install.packages("lme4")
# install.packages("lmerTest")
# install.packages("plyr")
# install.packages("reshape")
# install.packages("RColorBrewer")
# install.packages("colorRamps")
# install.packages("ggplot2")
# install.packages("grid")
# install.packages("data.table")
# install.packages("lubridate")
# install.packages("ggmap")
# install.packages("rjags")
# install.packages("runjags")


########################################################################################################

#Load required packages and functions in source files. If packages are not installed, use install.packages("knitr", "lme4", . . .).

message("Loading packages.")
library(knitr)
library(stringr)
library(lme4)
library(lmerTest)
library(plyr)
library(reshape)
library(RColorBrewer)
library(colorRamps)
library(ggplot2)
library(grid)
library(data.table)
library(lubridate)
library(ggmap)
library(rjags)
library(runjags)
library(scales)
library(spatstat)
library(ggsn)


#Load functions
message("Loading source files.")
source(paste(getwd(),"R_source","summaryFunction.R",sep="/"))
source(paste(getwd(), "R_source","summaryDelta.R",sep="/"))
source(paste(getwd(), "R_source","formatSETdata.R",sep="/"))
source(paste(getwd(), "R_source", "getSummariesAndPlot.R",sep="/"))
source(paste(getwd(), "R_source", "getSlope.R",sep="/"))
source(paste(getwd(), "R_source", "getSlopeObserver.R",sep="/"))
source(paste(getwd(), "R_source", "getSlopeBayesNoObs.R",sep="/"))
source(paste(getwd(), "R_source", "getSlopeBayesObserver.R",sep="/"))
source(paste(getwd(), "R_source", "obsEffect.R",sep="/"))
source(paste(getwd(), "R_source", "obsEffectBayes.R",sep="/"))
source(paste(getwd(), "R_source", "writeJAGSmodels.R", sep="/"))
source(paste(getwd(), "R_source", "runJAGSmodel.R",sep="/"))
source(paste(getwd(), "R_source", "getMaps.R",sep="/"))
source(paste(getwd(), "R_source", "getMapsBayes.R",sep="/"))
source(paste(getwd(), "R_source", "runBayesYearVisitObs.R",sep="/"))
source(paste(getwd(), "R_source", "runBayesYearVisitNoObs.R",sep="/"))
source(paste(getwd(), "R_source", "SummaryDataTable1.R",sep="/"))
########################################################################################################
#Read in raw data from file.
rawData<-read.csv(file=paste(getwd(),"Data" ,"Raw_SET_data","SET_All_Raw_Data_Combined_2016.csv" , sep="/"), header=TRUE)
head(rawData)

#Create table of unique SET stations
SETnames<-unique(rawData[,c("Unit_Code","Site_Name","Plot_Name")])
length(SETnames[,1])


#test names to get them straight
SETnames1<-unique(rawData$Plot_Name)

coords<-read.csv(paste(getwd(), "Data","SETcoords.csv",sep="/"),header=TRUE)
SETnames2<-unique(coords$Plot_Name)

setdiff(SETnames1, SETnames2)
#Save SET names as .csv file.
write.csv(SETnames, file=paste(getwd(),"Data","SETloc.csv",sep="/"),row.names=FALSE)



#data checking
PH.data.check<-subset(rawData, Plot_Name=="PH2-4 (Shallow)")
PH.data.formatted.check<-subset(rawData, Plot_Name=="PH2-4 (Shallow)")

########################################################################################################
#Fomat SET data (may take over 5 min to run, for full dataset). However, you only need to do this one time!)
message("Formatting SET data and computing the change in height (mm) of marsh elevation (make take several minutes).")
#Use system.time() wrapper around tasks to see how long run times are.
system.time(formatSETdata(dataIn=rawData))

########################################################################################################
#Read in formatted data.
data<-read.csv(paste(getwd(), "Data", "SET.delta.melt.csv",sep="/"),header=TRUE)
head(data)

#unique list of refuges
unique(sort(as.character(data$Unit_Code)))
length(unique(sort(as.character(data$Plot_Name))))

#get range of delta values
range(data$value)


summary.all<-summaryDelta(new.dataIn=data)
range(summary.all$mean)
summary.all<-data.frame(Unit_Code="All",summary.all)

########################################################################################################
#trim whitespace around Unit_Code
data$Unit_Code<-trimws(data$Unit_Code)

#make a table of NWRs included in data, and # of SET stations (Table 1.)

#remove "NWR" from RefugeName
data$RefugeName2<-gsub("NWR","",data$RefugeName)

data.table.1<-unique(data[,c("State","RefugeName2","Unit_Code")])

#remove NAs
data.table.1<-data.table.1[complete.cases(data.table.1),]

#get freqency of SET stations per refuge
sub.data<-unique(data[,c("State","RefugeName2","Unit_Code","Site_Name","Plot_Name")])
count.stations<-table(sub.data$Unit_Code, sub.data$Plot_Name)
count.stations.total<-rowSums(count.stations)
count.stations.total.1<-data.frame(Unit_Code=names(count.stations.total), nSET=as.data.frame(count.stations.total))

#get years sampled (start and end), and get number of observers
sub.data.years<-unique(data[,c("Unit_Code","Year","Last_Name")])

#get refugeList
refugeList<-unique(sort(as.character(sub.data.years$Unit_Code)))

#loop over refugeList to get range of years per refuge
range.years<-list()
for(i in 1:length(refugeList)){
  #subset data for 1 refuge
  new.data<-subset(sub.data.years, Unit_Code==refugeList[i])
  #get Unit_Code
  refugeName<-as.character(unique(new.data$Unit_Code))
  #get starting year
  startYear<-min(new.data$Year)
  #get ending year
  endYear<-max(new.data$Year)
  #get total years
  nYears<-endYear-startYear
  #get n unique observers (data recorders)
  nObs<-length(unique(new.data$Last_Name))

  range.years.1<-data.frame(Unit_Code=refugeName, startYear=startYear, endYear=endYear, nYears=nYears,nObs=nObs)
  range.years<-rbind(range.years, range.years.1)

}


#######################################################################################################
#get plotList
# plotList<-unique(sort(as.character(data$Plot_Name)))
# 
# n.Obs.out<-list()
# for(i in 1:length(plotList)){
#   #subset data for 1 refuge
#   new.data<-subset(data, Plot_Name==plotList[i])
#   #refugeName
#   refugeName<-as.character(unique(new.data$Unit_Code))
#   #get Unit_Code
#   plotName<-as.character(unique(new.data$Plot_Name))
#   #get starting year
#   startYear<-min(new.data$Year)
#   #get ending year
#   endYear<-max(new.data$Year)
#   #get total years
#   nYears<-endYear-startYear
#   #get n unique observers (data recorders)
#   nObs<-length(unique(new.data$Last_Name))
# 
#   n.Obs.1<-data.frame(Unit_Code=refugeName, startYear=startYear, endYear=endYear, nYears=nYears,nObs=nObs)
# 
#   n.Obs.out<-rbind(n.Obs.out, n.Obs.1)
# }
# n.Obs.out.unique<-unique(n.Obs.out)
#######################################################################################################


#combine with SummaryTable.1
data.table.2<-na.omit(merge(data.table.1, count.stations.total.1, by="Unit_Code",all.x=TRUE))

#now combine with range.years
data.table.3<-merge(data.table.2, range.years, by="Unit_Code")


#Now sort table by Latitude (reorder Unit_Code factor)
#create data.frame with NWRs ordered from North to South
RefugeOrder<-data.frame(Unit_Code=c("MEC","RHC","PKR","SPT","JHC","OYS","WRT","EBF","BMH","PMH", "PMH_shallow","ESV","BKB"), orderNum=seq(1:13))

#rename column headers
colnames(data.table.3)<-c("Unit_Code", "State","Refuge_Name",
                            "nSET","Start_Year","End_Year","nYears","nObservers")

#reorder columns
data.table.4<-data.table.3[,c("Refuge_Name","Unit_Code","State","nSET","nObservers","nYears","Start_Year","End_Year")]

#create integer column of order for NWRs
data.table.order<-merge(data.table.4, RefugeOrder, by=c("Unit_Code"),all.x=TRUE)
data.table.order<-data.table.order[order(data.table.order$orderNum),]
data.table.5<-data.table.order
data.table.5$orderNum<-NULL

#remove row.names from table
row.names(data.table.5)<-NULL

#create new Results folders to save table
dir.create(paste(getwd(),"Results",sep="/"))
dir.create(paste(getwd(),"Results","Tables",sep="/"))

#save Table as .csv file
write.csv(data.table.5,file=paste(getwd(),"Results","Tables","Table_1.csv",sep="/"),row.names=FALSE)

###################################################################
#Use getSummariesAndPlot function to first summarize data and then generate plots at multiple scales.
system.time(getSummariesAndPlot(dataIn=data))

###################################################################
#slopeFunction loop to average slopes from 4 positions and compile by stations
###################################################################
#get station level slopes (mean and SE) - no Observer effects
slopes<-getSlope(dataIn=data)
head(slopes)
#########################################################################
#get station level slopes (mean and SE) - with Observer effects
slopes.obs<-getSlopeObserver(dataIn=data)
head(slopes.obs)

########################################################################
#make table of Refuge-level slopes (Table 2)

#read in estimates without observer effects
refuge.noObs<-read.csv(paste(getwd(),"Results","Tables","All_Refuge-level_slopes.csv",sep="/"))
refuge.noObs.1<-refuge.noObs[,c("State","Unit_Code","n","mean","SE", "CV")]
colnames(refuge.noObs.1)<-c("State","Unit_Code","nSET","Mean","SE","CV")


#read in estimates with observer effects
refuge.Obs<-read.csv(paste(getwd(),"Results","Tables","All_Refuge-level_slopes.Obs.csv",sep="/"))
refuge.Obs.1<-refuge.Obs[,c("Unit_Code","mean","SE","CV")]
colnames(refuge.Obs.1)<-c("Unit_Code","Mean.obs","SE.obs","CV.obs")

#combine into one table
Freq.table<-merge(refuge.noObs.1, refuge.Obs.1, by="Unit_Code")

#add # of unique observers
range.years.obs<-range.years[,c("Unit_Code","nObs")]

Freq.table.2<-merge(Freq.table, range.years.obs, by="Unit_Code")

#reorder table columns
Freq.table.3<-Freq.table.2[,c("Unit_Code","State","nSET","nObs","Mean","SE","CV","Mean.obs", "SE.obs","CV.obs")]

#now sort by Latitude
#create integer column of order for NWRs
Freq.table.3.order<-merge(Freq.table.3, RefugeOrder, by=c("Unit_Code"),all.x=TRUE)
Freq.table.3.order<-Freq.table.3.order[order(Freq.table.3.order$orderNum),]
Freq.table.4<-Freq.table.3.order
Freq.table.4$orderNum<-NULL

#save Table as .csv file
write.csv(Freq.table.4,file=paste(getwd(),"Results","Tables","Table2_Frequentist_NWR_means.csv",sep="/"),row.names=FALSE)

########################################################################
#make table of Station-level slopes - No Observer effect (Table 3)

#read in estimates without observer effects
station.noObs<-read.csv(paste(getwd(),"Results","Tables","All_Station-level_slopes.csv",sep="/"))
names(station.noObs)
station.noObs.1<-station.noObs[,c("State","Unit_Code","Plot_Name","Lat","Long","mean","SE", "CV")]
colnames(station.noObs.1)<-c("State","Unit_Code","SET_Name","Lat","Long","Mean","SE","CV")

Table3.noObs<-station.noObs.1
head(Table3.noObs)
#add Obs model SET station-level estimates

#read in estimates with observer effects
station.Obs<-read.csv(paste(getwd(),"Results","Tables","All_Station-level_slopes.Obs.csv",sep="/"))
station.Obs.1<-station.Obs[,c("Unit_Code","Plot_Name","mean","SE","CV")]
colnames(station.Obs.1)<-c("Unit_Code","SET_Name","Mean.obs","SE.obs","CV.obs")

Table3.Obs<-station.Obs.1
head(Table3.Obs)

Table3.out<-merge(Table3.noObs, Table3.Obs, by="SET_Name")
head(Table3.out)
Table3.out$Unit_Code.y<-NULL

colnames(Table3.out)[3]<-"Unit_Code"
head(Table3.out)
#now sort table by Latitude to Longitude
#get order of NWR from .csv file (Unit_Order.csv)
RefugeOrder<-data.frame(Unit_Code=c("MEC","RHC","PKR","SPT","JHC","OYS","WRT","EBF","BMH","PMH", "PMH_shallow","ESV","BKB"), orderNum=seq(1:13))

Table5.out<-merge(Table3.out, RefugeOrder, by="Unit_Code")
head(Table5.out)
Table5.out$orderNum<-as.integer(Table5.out$orderNum)
Table6.out<-Table5.out[ order(Table5.out$orderNum,Table5.out$SET_Name,Table5.out$Lat),]
head(Table6.out)

#remove ordering columns
Table6.out$orderNum<-NULL
Table6.out$Lat<-NULL
Table6.out$Long<-NULL

#save Table as .csv file
write.csv(Table6.out,file=paste(getwd(),"Results","Tables","Table3_Frequentist_Station_means.csv",sep="/"),row.names=FALSE)

########################################################################
#Run obsEffect function to get difference between normalized slopes for models with and without Observer effects.
#read in data
slopes<-read.csv(paste(getwd(),"Results","Tables","All_Station-level_slopes.csv",sep="/"),header=TRUE)
slopes.obs<-read.csv(paste(getwd(),"Results","Tables","All_Station-level_slopes.Obs.csv",sep="/"),header=TRUE)

#Remove NAs from slopes and slopes.obs (Stations where less than 2 visits occurred)
slopes.1<-slopes[is.na(slopes$mean)==FALSE,]
slopes.obs.1<-slopes.obs[is.na(slopes.obs$mean)==FALSE,]

#Run obsEffect function to plot comparisons between models with and without Observer covariate.
obsEffect(dataIn_1=slopes.1,dataIn_2=slopes.obs.1)

#########################################################################
#make a map (radius = magnitude, color=trend direction)

#read in slope data
slopes<-read.csv(paste(getwd(),"Results","Tables","All_Station-level_slopes.csv",sep="/"),header=TRUE)
slopes$Model<-"No_Observer"

slopes.obs<-read.csv(paste(getwd(),"Results","Tables","All_Station-level_slopes.Obs.csv",sep="/"),header=TRUE)
slopes.obs$Model<-"Observer"

slopes.obs<-slopes.obs[is.na(slopes.obs$mean)==FALSE,]

#########################################################################
#Use getMaps function to plot SET trends and save maps.
getMaps(dataIn=slopes.obs)

#########################################################################

#build JAGS models with writeJAGSmodels function
writeJAGSmodels()

#########################################################################
#test getSlopeBayes functions

#first noObs
getSlopeBayesNoObs(dataIn=data)

#then with Observer covariate
getSlopeBayesObserver(dataIn=data)






#########################################################################
########################################################################
#make figure of comparison of CV among Frequentist and Bayesian noObs and Obs (Figure 2 and 3)

#frequentist
f.noObs<-read.csv(paste(getwd(),"Results","Tables","All_Station-level_slopes.csv",sep="/"))
f.noObs$Model<-"Frequentist_no_Observer"
f.noObs<-f.noObs[,c("Unit_Code","Site_Name","Plot_Name","mean","SD","Model")]
colnames(f.noObs)[4]<-"Mean"

f.Obs<-read.csv(paste(getwd(),"Results","Tables","All_Station-level_slopes.Obs.csv",sep="/"))
f.Obs$Model<-"Frequentist_Observer"
f.Obs<-f.Obs[,c("Unit_Code","Site_Name","Plot_Name","mean","SD","Model")]
colnames(f.Obs)[4]<-"Mean"

#Bayesian
b.noObs<-read.csv(paste(getwd(),"Results","Results_Bayes","Tables","All_Bayes_Station-level_slopes.csv",sep="/"))
b.noObs$Model<-"Bayesian_no_Observer"
b.noObs<-b.noObs[,c("Unit_Code","Site_Name","Plot_Name","Mean","SD","Model")]

b.Obs<-read.csv(paste(getwd(),"Results","Results_Bayes","Tables","All_Bayes_Station-level_slopes.Obs.csv",sep="/"))
b.Obs$Model<-"Bayesian_Observer"
b.Obs<-b.Obs[,c("Unit_Code","Site_Name","Plot_Name","Mean","SD","Model")]

#stack data
data.compare<-rbind(f.noObs, f.Obs, b.noObs, b.Obs)
head(data.compare)

#create Unit.Model column
data.compare$Unit.Model<-paste(data.compare$Unit_Code, data.compare$Model, sep=".")
#get summary of Means
data.summary<-summaryFunction(dataIn=data.compare, response="Mean",factor="Unit.Model")
colnames(data.summary)[3]<-"Mean"

#get relative standard deviation
data.summary$Rsd<-data.summary$SD/abs(data.summary$Mean)
head(data.summary)

Unit.and.Model<-read.table(text=as.character(data.summary$Unit.Model),sep=".",colClasses="character")
unique(data.compare$Model)

data.summary.2<-cbind(Unit.and.Model, data.summary)
colnames(data.summary.2)[1:2]<-c("Unit_Code","Model")
head(data.summary.2)

data.summary.2$Unit_Code <- factor(data.summary.2$Unit_Code,
                                   levels = c("MEC","RHC","PKR","SPT","JHC","OYS","WRT","EBF","BMH","PMH", "PMH_shallow","ESV","BKB"))

#now Plot
#now plot unit results
compare.plot<-ggplot(data.summary.2)+
  #coord_flip()+
  aes(x=Unit_Code, y=Rsd)+
  geom_bar(stat="identity",position="dodge",aes(fill=Model))+
  ylim(0,8)+
  theme(panel.background=element_rect(fill='white',color="black"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.line=element_line(color="black"))+
  theme(panel.background=element_rect(color="black"))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=12, color="black"),
        axis.title.x = element_text(size=13, hjust=0.5, vjust=1.9),
        axis.title.y = element_text(angle = 90, vjust=1.2, size=13))+
  labs(x="SET stations", y="Percent relative SD")
compare.plot
compare.plot.out<-compare.plot+facet_wrap(~Unit_Code,scales="free_x")

ggsave(compare.plot.out, file="Percent_relative_SD_NWRs_figure.png",path=paste(getwd(),"Results","Figures",sep="/"))
print(compare.plot.out)
###########################################################################
data.compare.2<-na.omit(data.compare)

data.summary.2$Unit_Code <- factor(data.summary.2$Unit_Code,
                         levels = c("MEC","RHC","PKR","SPT","JHC","OYS","WRT","EBF","BMH","PMH", "PMH_shallow","ESV","BKB"))

#now plot unit results
compare.mean.plot<-ggplot(data.summary.2)+
  #coord_flip()+
  aes(x=Unit_Code, y=Mean, ymin=Mean-SD, ymax=Mean+SD,group=Model)+
  geom_errorbar(width=0.1, size=0.5, position=position_dodge(width=0.95),aes(color=Model))+
  geom_bar(stat="identity",position="dodge",aes(fill=Model))+
  #ylim(-100,100)+
  theme(panel.background=element_rect(fill='white',color="black"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.line=element_line(color="black"))+
  theme(panel.background=element_rect(color="black"))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=12, color="black"),
        axis.title.x = element_text(size=13, hjust=0.5, vjust=1.9),
        axis.title.y = element_text(angle = 90, vjust=1.2, size=13))+
  labs(x="SET stations", y="Change in SET height (mm/year)")
compare.mean.plot
compare.mean.plot.out<-compare.mean.plot+facet_wrap(~Unit_Code,scales="free_x")

ggsave(compare.mean.plot.out, file="Mean_change_in_elevation_NWRs_figure.png",path=paste(getwd(),"Results","Figures",sep="/"))

print(compare.mean.plot.out)


#########################################################################
#########################################################################
########################################################################
#make table of Bayes-estimated Refuge-level slopes (Table 4)

#read in estimates without observer effects
refuge.noObs<-read.csv(paste(getwd(),"Results","Results_Bayes","Tables","All_Bayes_Refuge-level_slopes.csv",sep="/"))
refuge.noObs.1<-refuge.noObs[,c("State","Unit_Code","n","mean","SE")]
colnames(refuge.noObs.1)<-c("State","Unit_Code","nSET","Mean","SE")


#read in estimates with observer effects
refuge.Obs<-read.csv(paste(getwd(),"Results","Results_Bayes","Tables","All_Bayes_Refuge-level_slopes.Obs.csv",sep="/"))
refuge.Obs.1<-refuge.Obs[,c("Unit_Code","mean","SE")]
colnames(refuge.Obs.1)<-c("Unit_Code","Mean.obs","SE.obs")

#combine into one table
Freq.table<-merge(refuge.noObs.1, refuge.Obs.1, by="Unit_Code")


#get years sampled (start and end), and get number of observers
sub.data.years<-unique(data[,c("Unit_Code","Year","Last_Name")])

#get refugeList
refugeList<-unique(sort(as.character(sub.data.years$Unit_Code)))

#loop over refugeList to get range of years per refuge
range.years<-list()
for(i in 1:length(refugeList)){
  #subset data for 1 refuge
  new.data<-subset(sub.data.years, Unit_Code==refugeList[i])
  #get Unit_Code
  refugeName<-as.character(unique(new.data$Unit_Code))
  #get starting year
  startYear<-min(new.data$Year)
  #get ending year
  endYear<-max(new.data$Year)
  #get total years
  nYears<-endYear-startYear
  #get n unique observers (data recorders)
  nObs<-length(unique(new.data$Last_Name))

  range.years.1<-data.frame(Unit_Code=refugeName, startYear=startYear, endYear=endYear, nYears=nYears,nObs=nObs)
  range.years<-rbind(range.years, range.years.1)

}

#add # of unique observers
range.years.obs<-range.years[,c("Unit_Code","nObs")]

Freq.table.2<-merge(Freq.table, range.years.obs, by="Unit_Code")

#reorder table columns
Freq.table.3<-Freq.table.2[,c("Unit_Code","State","nSET","nObs","Mean","SE","Mean.obs", "SE.obs")]

#now sort by Latitude
#create integer column of order for NWRs
Freq.table.3.order<-merge(Freq.table.3, RefugeOrder, by=c("Unit_Code"),all.x=TRUE)
Freq.table.3.order<-Freq.table.3.order[order(Freq.table.3.order$orderNum),]
Freq.table.4<-Freq.table.3.order
Freq.table.4$orderNum<-NULL

#save Table as .csv file
write.csv(Freq.table.4,file=paste(getwd(),"Results","Results_Bayes","Tables","Table4_Bayesian_NWR_means.csv",sep="/"),row.names=FALSE)

########################################################################
#make table of Station-level slopes - No Observer effect (Table 5)

#read in estimates without observer effects
station.noObs<-read.csv(paste(getwd(),"Results","Results_Bayes","Tables","All_Bayes_Station-level_slopes.csv",sep="/"))
names(station.noObs)
station.noObs.1<-station.noObs[,c("State","Unit_Code","Plot_Name","Lat","Long","Mean","Time.series.SE")]
colnames(station.noObs.1)<-c("State","Unit_Code","SET_Name","Lat","Long","Mean","SE")

Table3.noObs<-station.noObs.1
head(Table3.noObs)
#add Obs model SET station-level estimates

#read in estimates with observer effects
station.Obs<-read.csv(paste(getwd(),"Results","Results_Bayes","Tables","All_Bayes_Station-level_slopes.Obs.csv",sep="/"))
station.Obs.1<-station.Obs[,c("Unit_Code","Plot_Name","Mean","Time.series.SE")]
colnames(station.Obs.1)<-c("Unit_Code","SET_Name","Mean.obs","SE.obs")

Table3.Obs<-station.Obs.1
head(Table3.Obs)

Table3.out<-merge(Table3.noObs, Table3.Obs, by="SET_Name")
head(Table3.out)
Table3.out$Unit_Code.y<-NULL

colnames(Table3.out)[3]<-"Unit_Code"
head(Table3.out)
#now sort table by Latitude to Longitude
#get order of NWR from .csv file (Unit_Order.csv)
RefugeOrder<-data.frame(Unit_Code=c("MEC","RHC","PKR","SPT","JHC","OYS","WRT","EBF","BMH","PMH", "PMH_shallow","ESV","BKB"), orderNum=seq(1:13))

Table5.out<-merge(Table3.out, RefugeOrder, by="Unit_Code")
head(Table5.out)
Table5.out$orderNum<-as.integer(Table5.out$orderNum)
Table6.out<-Table5.out[ order(Table5.out$orderNum,Table5.out$SET_Name,Table5.out$Lat),]
head(Table6.out)

#remove ordering columns
Table6.out$orderNum<-NULL
Table6.out$Lat<-NULL
Table6.out$Long<-NULL

#save Table as .csv file
write.csv(Table6.out,file=paste(getwd(),"Results","Results_Bayes","Tables","Table5_Bayesian_Station_means.csv",sep="/"),row.names=FALSE)
########################################################################


#########################################################################
#Now for Bayes estimates, make a map (radius = magnitude, color=trend direction)

#read in slope data
slopesB<-read.csv(paste(getwd(),"Results","Results_Bayes","Tables","All_Bayes_Station-level_slopes.csv",sep="/"),header=TRUE)
slopesB<-slopesB[is.na(slopesB$Mean)==FALSE,]

# slopesB$Model<-"No_Observer"
# slopesB<-slopesB[is.na(slopesB$Mean)==FALSE,]

slopes.obsB<-read.csv(paste(getwd(),"Results","Results_Bayes","Tables","All_Bayes_Station-level_slopes.Obs.csv",sep="/"),header=TRUE)
slopes.obsB$Model<-"Observer"
slopes.obsB<-slopes.obsB[is.na(slopes.obsB$Mean)==FALSE,]

#Run obsEffect function to plot comparisons between models with and without Observer covariate.
obsEffectBayes(dataIn_1=slopesB,dataIn_2=slopes.obsB)

#########################################################################
#plot means of Bayes with and without Observer covariate




#########################################################################
#Use getMaps function to plot SET trends and save maps.

slopes.obsB<-read.csv(paste(getwd(),"Results","Results_Bayes","Tables","All_Bayes_Station-level_slopes.Obs.csv",sep="/"),header=TRUE)
slopes.obsB$Model<-"Observer"

getMapsBayes(dataIn=slopes.obsB)
head(slopes.obsB)
#########################################################################
#Appendix B
set.data<-read.csv(paste(getwd(),"Data","SET.delta.melt.csv",sep="/"),header=TRUE)

#full refugeList
refugeList<-c("MEC","RHC","PKR","SPT","JHC","OYS","WRT","EBF","BMH","PMH", "PMH_shallow","ESV","BKB")

for (i in 1:length(refugeList)){
  refuge.data<-subset(set.data, Unit_Code==refugeList[1])

  smiUnitList<-sort(unique(as.character(refuge.data$Site_Name)))

  site.delta.SET<-list()
  for(j in 1:length(smiUnitList)){
    new.data<-subset(refuge.data, Site_Name==smiUnitList[1])
    siteName<-unique(as.character(new.data$Site_Name))
    refugeName<-unique(as.character(new.data$Unit_Code))
    print(refugeName)
    stateName<-unique(as.character(new.data$State))


    #plot regional estimates
    nStations=length(unique(new.data$Plot_Name))


    new.summary<-summaryDelta(new.dataIn=new.data)
    new.summary<-data.frame(State=stateName,Unit_Code=refugeName,Site_Name=siteName, new.summary)

    #now plot unit results
    SETplot.site<-ggplot(new.summary)+
      aes(x=year.visit, y=mean,group=1)+
      geom_smooth(method=lm,fullrange=FALSE, linetype=1, color="white")+
      geom_errorbar(ymax=new.summary$mean+new.summary$SE, ymin=new.summary$mean-new.summary$SE, width=0.2, size=0.6, color=I("grey30"))+
      geom_line(color="deepskyblue",linetype=1, size=1)+
      geom_point(size=2, color="deepskyblue")+
      theme(panel.background=element_rect(fill='white',color="black"))+
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
      theme(axis.line=element_line(color="black"))+
      theme(panel.background=element_rect(color="black"))+
      theme(axis.text.x = element_text(angle = 55, hjust = 1.2,vjust=1.2, size=10, color="black"),
            axis.text.y = element_text(size=12, color="black"),
            axis.title.x = element_text(size=13, hjust=0.5, vjust=1.9),
            axis.title.y = element_text(angle = 90, vjust=1.2, size=13))+
      ylim(min(new.summary$mean-new.summary$SE)-30,max(new.summary$mean+new.summary$SE)+30)+
      #xlim(200.8,2014.2)
      labs(x="Year (visits)", y="Change in elevation (mm)")
    #scale_x_discrete(limits=c(2009,2014),breaks=c(2009,2010,2011,2012,2013,2014,2015), "Year")+
    #scale_y_continuous(limits=c(-5,35),breaks=c(-5,0,5,10,15,20,25,30,35), "Change in elevation (mm)")
    plot.5<-SETplot.site+
      ggtitle(paste(siteName," SET data (n = ", nStations,")", "in ",refugeName,".",sep=""))

    print(plot.5)

    cat("\n\n\\pagebreak\n")

    #save figure
    # myFilepath<-paste(getwd(), "Results","Data_Summary_Results","SMI_unit_Summary",refugeName,siteName,sep="/")
    # ggsave(plot.5, filename=paste(siteName, "_delta_SET_year_visit.pdf",sep=""),path=myFilepath, width=9,height=6.5, limitsize=FALSE)


    # write.csv(site.delta.SET, file=paste(getwd(), "Results","Data_Summary_Results", paste("All_SMI_units", "_delta_SET_year_visit.csv",sep=""),sep="/"),row.names=FALSE)
  }
}
########################################################################################################################
#Appendix C
#read
set.data<-read.csv(paste(getwd(),"Data","SET.delta.melt.csv",sep="/"),header=TRUE)

#full refugeList
refugeList<-c("MEC","RHC","PKR","SPT","JHC","OYS","WRT","EBF","BMH","PMH", "PMH_shallow","ESV","BKB")

for (i in 1:length(refugeList)){
  refuge.data<-subset(set.data, Unit_Code==refugeList[i])

  smiUnitList<-sort(unique(as.character(refuge.data$Site_Name)))

  for(j in 1:length(smiUnitList)){
    new.data<-subset(refuge.data, Site_Name==smiUnitList[j])

    stationList<-sort(unique(as.character(new.data$Plot_Name)))

    for(k in 1:length(stationList)){

      new.data<-subset(set.data, Plot_Name==stationList[k])
      stationName<-unique(as.character(new.data$Plot_Name))
      siteName<-unique(as.character(new.data$Site_Name))
      refugeName<-unique(as.character(new.data$Unit_Code))

      new.summary<-summaryDelta(new.dataIn=new.data)
      new.summary<-data.frame(State=stateName,Unit_Code=refugeName,Site_Name=siteName,Plot_Name=stationName, new.summary)

      #now plot unit results
      SETplot.station<-ggplot(new.summary)+
        aes(x=year.visit, y=mean,group=1)+
        geom_smooth(method=lm,fullrange=FALSE, linetype=1, color="white")+
        geom_errorbar(ymax=new.summary$mean+new.summary$SE, ymin=new.summary$mean-new.summary$SE, width=0.2, size=0.6, color=I("grey30"))+
        geom_line(color="deepskyblue",linetype=1, size=1)+
        geom_point(size=2, color="deepskyblue")+
        theme(panel.background=element_rect(fill='white',color="black"))+
        theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
        theme(axis.line=element_line(color="black"))+
        theme(panel.background=element_rect(color="black"))+
        theme(axis.text.x = element_text(angle = 55, hjust = 1.2,vjust=1.2, size=10, color="black"),
              axis.text.y = element_text(size=12, color="black"),
              axis.title.x = element_text(size=13, hjust=0.5, vjust=1.9),
              axis.title.y = element_text(angle = 90, vjust=1.2, size=13))+
        ylim(min(new.summary$mean-new.summary$SE)-30,max(new.summary$mean+new.summary$SE)+30)+
        #xlim(200.8,2014.2)
        labs(x="Year (visits)", y="Change in elevation (mm)")
      #scale_x_discrete(limits=c(2009,2014),breaks=c(2009,2010,2011,2012,2013,2014,2015), "Year")+
      #scale_y_continuous(limits=c(-5,35),breaks=c(-5,0,5,10,15,20,25,30,35), "Change in elevation (mm)")
      plot.6<-SETplot.station+
        ggtitle(paste(stationName," SET data (n = ", nStations,")", " in ",refugeName,".",sep=""))

      print(plot.6)


###################################################################################################################################
##SET Appendix - Refuge Summaries (based on DE Atlas unit summaries)
      #read in required data
      
      #Read in raw data from file.
      rawData<-read.csv(file=paste(getwd(),"Data" ,"Raw_SET_data","SET_All_Raw_Data_Combined_2016.csv" , sep="/"), header=TRUE)
      #Read in formatted data.
      data<-read.csv(paste(getwd(), "Data", "SET.delta.melt.csv",sep="/"),header=TRUE)
      #read in slope data
      Table3<-read.csv(paste(getwd(),"Results","Tables","Table3_Frequentist_Station_means.csv",sep="/"),header=TRUE)
      
      #read in data melt
      
      #remove NAs
      Table3.1<-na.omit(Table3)
      
      #redefine data
      sub.data.1<-Table3.1
      
      #remove NAs
      sub.data<-sub.data.1[is.na(sub.data.1$Mean)==FALSE,]
      
      #Add Trend column
      sub.data$Trend<-ifelse(sub.data$Mean<0,"Negative",
                             ifelse(sub.data$Mean>0,"Positive","Zero"))
      sub.data$Trend<-as.factor(as.character(sub.data$Trend))
      
      #add refuge Names 
      refugeNames<-na.omit(unique(data[,c("Unit_Code","RefugeName")]))
      
      sub.data.merge<-merge(sub.data, refugeNames, by=c("Unit_Code"))
      
      #remove "NWR" from RefugeName
      sub.data.merge$RefugeName<-gsub("NWR","",sub.data.merge$RefugeName)
      sub.data.merge$RefugeName<-trimws(sub.data.merge$RefugeName)
      sub.data.merge$SET_Name<-as.character(sub.data.merge$SET_Name)
      sub.data.merge$SET_Name<-gsub("*","",sub.data.merge$SET_Name, fixed=TRUE)
      
      #add Lat/Long coords
      latLong<-na.omit(unique(data[,c("Plot_Name","Lat","Long")]))
      colnames(latLong)<-c("SET_Name","Latitude","Longitude")
      
      #merge with trend data
      data.merge<-merge(sub.data.merge, latLong, by=c("SET_Name"),all.x=TRUE)
      
      data.merge$Longitude<-as.numeric(as.character(data.merge$Longitude))
      data.merge$Latitude<-as.numeric(as.character(data.merge$Latitude))
      
      #change column header names of data
      colnames(data.merge)<-c("SET_Name","Unit_Code","State","Mean","SE","CV","Mean.obs","SE.obs","CV.obs","Trend","RefugeName",
                              "Latitude","Longitude")
      ############################################################################################################################
      #get refugeList
      refugeList<-c("MEC","RHC","PKR","SPT","JHC","OYS","WRT","EBF","BMH","PMH", "PMH_shallow","ESV","BKB")
      
      #refugeList<-c("EBF")
      
      #now loop
      for(i in 1:length(refugeList)){
        
        new.data<-subset(data.merge, Unit_Code==refugeList[i])
        
        refugeName<-unique(as.character(new.data$RefugeName))
        
        #get number of species at each point, and for park
        nSETs<-length(unique(new.data$SET_Name))
        
        setList<-paste(unique(as.character(new.data$SET_Name)), collapse=", ")
        
        #get dates
        sub.all.data<-subset(data, Unit_Code==refugeList[i])
        startYear<-min(sub.all.data$Year)
        endYear<-max(sub.all.data$Year)
        
        #get min and max trends
        minTrend<-round(min(new.data$Mean.obs),2)
        maxTrend<-round(max(new.data$Mean.obs),2)
        
        #get average and SE
        meanTrend<-round(mean(new.data$Mean.obs), 2)
        seTrend<-round(sd(new.data$Mean.obs)/sqrt(nSETs),2)
        
        maxLon=max(new.data$Longitude,na.rm=TRUE)*0.9992
        minLon=min(new.data$Longitude,na.rm=TRUE)*1.0008
        minLat=min(new.data$Latitude,na.rm=TRUE)*0.9992
        maxLat=max(new.data$Latitude,na.rm=TRUE)*1.0008
        
        bbox.1<-c(left=minLon, bottom=minLat,right=maxLon,top=maxLat)
        
        #get base map (requires internet connection)
        map.1 <- get_map(location = bbox.1, maptype="terrain")
        #ggmap(map.1)
        
        #get true map coords
        map.minLon=as.numeric(attr(map.1, "bb")[2])
        map.maxLon=as.numeric(attr(map.1, "bb")[4])
        map.minLat=as.numeric(attr(map.1, "bb")[1])
        map.maxLat=as.numeric(attr(map.1, "bb")[3])
        
        #######################################################################
        #add text for header
        
        par(family="Times",ps=14)
        plot.new()
        text(0.5,0.5, paste(refugeName," NWR ", "Summary", sep=""))
        
        #add summary paragraph
        cat("Surface elevation table (SET) data was collected within ", refugeName," NWR ", " at the following ", nSETs, " SET stations ","(",setList,")", " between ",startYear, " and ", endYear, ".", " Trends ranged from ", minTrend, " to ", maxTrend,".", " Overall, the average refuge-level trend was ",meanTrend," \u00B1 ", seTrend," (mean"," \u00B1 ","SE).",sep="")
        
        #add 2 blank lines
        cat("\\newline")
        cat("\\newline")
        
        #######################################################################
        color.df<-data.frame(Trend=c("Positive","Negative"),
                             TrendOrder=c(1,2),
                             myColors=c("royalblue2","red"))
        
        new.data.2<-merge(new.data, color.df, by=c("Trend"),all.x=TRUE) 
        new.data.3<-new.data.2[order(new.data.2$TrendOrder,decreasing=FALSE),]
        
        myColors<-paste(unique(new.data.3$myColors))
        
        #scalebar placement 
        scalebar.x<-map.minLon+(abs(map.maxLon - map.minLon)/2)
        scalebar.y<-(abs(map.maxLat - map.minLat)/10)+map.minLat
        
        scalebar.dist<-as.numeric(ifelse(abs(map.maxLon - map.minLon) < 0.03, 0.5,1))
        
        #adjust factor levels
        new.data.3$Trend<-factor(new.data.3$Trend, levels=c("Positive","Negative"))
        
        #put points on map
        map.all<-ggmap(map.1)+
          geom_point(data=new.data.3, stat="identity",aes(x=Longitude, y=Latitude, color=Trend,size=abs(Mean.obs)),alpha=0.7)+
          # geom_label_repel(data=new.data.2, stat="identity", aes(x=Longitude, y=Latitude,label=as.character(SpeciesCount)),
          #           label.size=0.15,color="black",alpha=0.8,max.iter = 100)+
          scale_color_manual(values=myColors)+
          # scale_color_gradientn(colors=colorFunc(20),limits=c(0,190))+
          guides(color=guide_legend(order=1),size=guide_legend(title = "Magnitude",order=2))+
          labs(x="Longitude",y="Latitude")+
          #ggtitle(refugeName)+
          theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
          scalebar(anchor=c(x=scalebar.x, y=scalebar.y), x.min=map.minLon, x.max=map.maxLon, y.min=map.minLat, y.max=map.maxLat,
                   dist=scalebar.dist, dd2km= TRUE, model='WGS84',st.dist=0.02, st.size=4,height=0.02)
        #map.all
        
        
        
        #add Species Richness
        plot.new()
        text(0.5,0.5,paste("Trends of changes in Salt marsh Elevation (mm)"))
        
        #add figure heading
        cat(paste("Figure 1. ","Map showing locations, trend directionality, and trend magnitude of SET stations within ",refugeName, ".", " Radii of positve (blue) and negative (red) points show the magnitude of linear trends.",sep=""))
        
        print(map.all)
        
        # #add pagebreak
        cat("\\pagebreak")
        
        #add Species Density
        plot.new()
        text(0.5,0.5,paste(" "))
        
        
        cat(paste("Table 1. ", "Trends of surface elevation table (SET) stations within ",refugeName,"."," Means and SEs for trends are shown for both linear models with and without observers.",sep=""))
        
        new.data.4<-new.data.3[,c("SET_Name","Mean","SE","CV","Mean.obs","SE.obs","CV.obs")]
        row.names(new.data.4)<-NULL
        
        
        #kable(new.data.4)
        print(xtable(new.data.4,caption=NULL),comment=FALSE, floating=TRUE)
        
        cat("\n\\clearpage\n")
        
        plot.new()
        text(0.5,0.5,paste("Summary of Raw Data - Surface Elevation change (mm)"))
        
        new.data<-subset(data, Unit_Code==refugeList[i])
        
        #print(unique(as.character(new.data$RefugeName)))
        
        refugeName<-unique(as.character(new.data$RefugeName))
        stateName<-unique(as.character(new.data$State))
        
        new.summary<-summaryDelta(new.dataIn=new.data)
        new.summary<-data.frame(State=stateName,Unit_Code=refugeName[1], new.summary)
        
        #plot regional estimates
        nStations=length(unique(new.data$Plot_Name))
        
        #now plot unit results
        NWRplot<-ggplot(new.summary)+
          aes(x=year.visit, y=mean,group=1)+
          geom_smooth(method=lm,fullrange=FALSE, linetype=1, color="white")+
          geom_errorbar(ymax=new.summary$mean+new.summary$SE, ymin=new.summary$mean-new.summary$SE, width=0.2, size=0.6, color=I("grey50"))+
          geom_line(color="deepskyblue3",linetype=1, size=1)+
          geom_point(size=1.5, color="deepskyblue3")+
          theme(panel.background=element_rect(fill='white',color="black"))+
          theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
          theme(axis.line=element_line(color="black"))+
          theme(panel.background=element_rect(color="black"))+
          theme(axis.text.x = element_text(angle = 55, hjust = 1.2,vjust=1.2, size=10, color="black"),
                axis.text.y = element_text(size=12, color="black"),
                axis.title.x = element_text(size=13, hjust=0.5, vjust=1.9),
                axis.title.y = element_text(angle = 90, vjust=1.2, size=13))+
          #ylim(min(new.summary$mean-new.summary$SE)-30,max(new.summary$mean+new.summary$SE)+30)+
          #xlim(200.8,2014.2)
          labs(x="Year (visits)", y="Elevation change (mm)")+
          #scale_x_discrete(limits=c(2009,2014),breaks=c(2009,2010,2011,2012,2013,2014,2015), "Year")+
          scale_y_continuous(expand=c(0.5,0.5))
        
        
        #add figure caption
        cat(paste("Figure 2. ", "Elevation change (mm) within ",refugeName,"."," Means and SEs for each year (and within year visits) are shown along with overall linear trend.",sep=""))
        
        print(NWRplot)
        ###########################################################################################
        #add pagebreak
        cat("\n\\clearpage\n")
        
        # plot.new()
        # text(0.5,0.5,paste("Summary of Raw Data (Site-level estimates)"))
        
        new.data$Site_Name<-as.factor(as.character(new.data$Site_Name))
        
        smiUnitList<-sort(unique(as.character(new.data$Site_Name)))
        
        summary.out<-list()
        for(j in 1:length(smiUnitList)){
          sub.data<-subset(new.data, Site_Name==smiUnitList[j])
          
          siteName<-unique(as.character(sub.data$Site_Name))
          new.summary<-summaryDelta(new.dataIn=sub.data)
          new.summary<-data.frame(State=stateName,Unit_Code=refugeName[1],Site_Name=siteName, new.summary)
          
          
          #plot regional estimates
          nStations=length(unique(sub.data$Plot_Name))
          
          
          summary.out<-rbind(summary.out, new.summary)
        }
        
        
        #now plot unit results
        SETplot.site<-ggplot(summary.out)+
          aes(x=year.visit, y=mean,group=1)+
          geom_smooth(method=lm,fullrange=FALSE, linetype=1, color="white")+
          geom_errorbar(ymax=summary.out$mean+summary.out$SE, ymin=summary.out$mean-summary.out$SE, width=0.2, size=0.6, color=I("grey30"))+
          geom_line(color="chartreuse4",linetype=1, size=1)+
          geom_point(size=2, color="chartreuse4")+
          theme(panel.background=element_rect(fill='white',color="black"))+
          theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
          theme(axis.line=element_line(color="black"))+
          theme(panel.background=element_rect(color="black"))+
          theme(axis.text.x = element_text(angle = 55, hjust = 1.2,vjust=1.2, size=10, color="black"),
                axis.text.y = element_text(size=12, color="black"),
                axis.title.x = element_text(size=13, hjust=0.5, vjust=1.9),
                axis.title.y = element_text(angle = 90, vjust=1.2, size=13))+
          #ylim(min(summary.out$mean-summary.out$SE)-20,max(summary.out$mean+summary.out$SE)+20)+
          #xlim(200.8,2014.2)
          labs(x="Year (visits)", y="Evelation change (mm)")+
          #scale_x_discrete(limits=c(2009,2014),breaks=c(2009,2010,2011,2012,2013,2014,2015), "Year")+
          scale_y_continuous(expand=c(0.5,0.5))
        SitePlot<-SETplot.site+facet_wrap(~Site_Name, scales="free",ncol=3)
        
        #add figure caption
        cat(paste("Figure 3. ", "Summary of Site-level elevation change (mm) data within ",refugeName,"."," Means and SEs (gray bars) for each year (and within year visits) are shown along with overall linear trends. Note potentially different scales for both x and y axes.",sep=""))
        
        print(SitePlot)
        
        
        cat("\n\\clearpage\n")
        
        ##########################################################################################
        # plot.new()
        # text(0.5,0.5,paste("Summary of Raw Data (SET station-level estimates)"))
        
        new.data$Site_Name<-as.factor(as.character(new.data$Site_Name))
        
        smiUnitList<-sort(unique(as.character(new.data$Site_Name)))
        
        
        site.out<-list()
        for(j in 1:length(smiUnitList)){
          sub.data<-subset(new.data, Site_Name==smiUnitList[j])
          
          stationList<-sort(unique(as.character(sub.data$Plot_Name)))
          
          summary.out<-list()
          for(k in 1:length(stationList)){
            
            station.data<-subset(sub.data, Plot_Name==stationList[k])
            stationName<-unique(as.character(station.data$Plot_Name))
            siteName<-unique(as.character(station.data$Site_Name))
            refugeName<-unique(as.character(station.data$RefugeName))
            
            new.summary<-summaryDelta(new.dataIn=station.data)
            new.summary<-data.frame(State=stateName,Unit_Code=refugeName,Site_Name=siteName,Plot_Name=stationName, new.summary)
            summary.out<-rbind(summary.out, new.summary)
          }
          
          site.out<-rbind(site.out, summary.out)
        }
        
        colNum<-ifelse(length(unique(site.out$Plot_Name)) < 6, 3, 4)
        
        #now plot unit results
        SETplot.station<-ggplot(site.out)+
          aes(x=year.visit, y=mean,group=1)+
          geom_smooth(method=lm,fullrange=FALSE, linetype=1, color="white")+
          geom_errorbar(ymax=site.out$mean+site.out$SE, ymin=site.out$mean-site.out$SE, width=0.2, size=0.6, color=I("grey30"))+
          geom_line(color="brown3",linetype=1, size=1)+
          geom_point(size=2, color="brown3")+
          theme(panel.background=element_rect(fill='white',color="black"))+
          theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
          theme(axis.line=element_line(color="black"))+
          theme(panel.background=element_rect(color="black"))+
          theme(axis.text.x = element_text(angle = 55, hjust = 1.2,vjust=1.2, size=9, color="black"),
                axis.text.y = element_text(size=10, color="black"),
                axis.title.x = element_text(size=13, hjust=0.5, vjust=1.9),
                axis.title.y = element_text(angle = 90, vjust=1.2, size=13))+
          #ylim(min(site.out$mean-site.out$SE)-50,max(site.out$mean+site.out$SE)+50)+
          #xlim(200.8,2014.2)
          labs(x="Year (visits)", y="Elevation change (mm)")+
          #scale_x_discrete(limits=c(2009,2014),breaks=c(2009,2010,2011,2012,2013,2014,2015), "Year")+
          scale_y_continuous(expand=c(1,1))
        StationPlot<-SETplot.station+facet_wrap(~Plot_Name, scales="free",ncol=colNum)
        
        #add figure caption
        cat(paste("Figure 4. ", "Summary of SET station-level elevation change (mm) data within ",refugeName,"."," Means and SEs (gray bars) for each year (and within year visits) are shown along with overall linear trends. Note potentially different scales for both x and y axes.",sep=""))
        
        print(StationPlot)
        
        ##########################################################################################
        
        cat("\n\\clearpage\n")
      }
   
      
      
      
      
      ################################################################################################
      
      
      #not sure why, but needed to fill in Refuge Names for NAs for some reason
      
      #load data
      all.data <- read.csv(paste(getwd(),"SET.delta.melt.csv",sep="/"),header=TRUE)
      
      
      refugeName.df<-unique(na.omit(all.data[,c("RefugeName","Unit_Code")]))
      refugeName.df<-refugeName.df[refugeName.df$Unit_Code !="PMH_shallow",]
      
      
      #merge with all.data
      all.data$Unit_Code<-as.character(all.data$Unit_Code)
      all.data$Unit_Code<-gsub("_shallow","",all.data$Unit_Code)
      
      all.data.new<-merge(all.data,refugeName.df, by="Unit_Code",all.x=TRUE)
      all.data.new$RefugeName.x<-NULL
      colnames(all.data.new)[colnames(all.data.new) == "RefugeName.y"] <- "RefugeName"
      
      #get refuge Names
      all.data.new$RefugeName<-as.character(all.data.new$RefugeName)
      all.data.new$RefugeName<-gsub(" Shallow","",all.data.new$RefugeName)
      #remove "NWR"
      all.data.new$RefugeName<-gsub(" NWR","",all.data.new$RefugeName)
      unique(all.data.new$RefugeName)
      
      #get refuge List
      refugeList<-sort(as.character(unique(all.data.new$RefugeName)))
      
      #write as .csv
      write.csv(all.data.new, file=paste(getwd(),"set.data.csv",sep="/"),row.names=FALSE)
      