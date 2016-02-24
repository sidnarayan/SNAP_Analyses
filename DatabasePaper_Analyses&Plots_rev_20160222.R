library(ggplot2)
library(relaimpo)
library(rgdal)
library(ggmap)
library(dplyr)
library(tidyr)
library(magrittr)
library(maps)
library(mapdata)
library(rworldmap)
require(gridExtra)
library(ggrepel)# for plotting multiple plots 

# Data Input 
alldata_field<-read.table("C:/Rfiles/SNAP/WaveReductionData_20160222.txt", header = TRUE, sep="\t",skip=1)
alldata_projects<-read.table("C:/Rfiles/SNAP/Projectdata_20160223.txt", header = TRUE, sep="\t",skip=1)

fielddata<-read.table("C:/Rfiles/SNAP/WaveReductionData_20160222.txt", header = TRUE, sep="\t",skip=5)
projectdata<-read.table("C:/Rfiles/SNAP/Projectdata_20160223.txt", header = TRUE,skip=5)

#Dataframe Transformations and Calculations 
fielddata<-fielddata %>%
  mutate(
    hvbyh=hv/h,
    Hbyh=Hs.i/h,
    HbyB=Hs.i/B,
    HabitatGrp=ifelse(Habitat=="CRC"|Habitat=="CRF"|Habitat=="CRW","CR",as.character(Habitat))
  )
projectdata<-projectdata %>%
  mutate(
    Reported.Cost.Mean=ifelse(Reported.Cost.V2==0.0,Reported.Cost.V1,((Reported.Cost.V1+Reported.Cost.V2)/2)),
    Habitat.Cost.Reported.Total=ifelse(Cost.Unit=="Total",Reported.Cost.Mean,
                                       ifelse(Cost.Unit=="Ha",Reported.Cost.Mean*(Area/10000),
                                              Reported.Cost.Mean*Area)),
    Habitat.Cost.USD2014=Habitat.Cost.Reported.Total*CPI.2014*USD.2014,
    Habitat.Cost.Per.SquareMeter=Habitat.Cost.USD2014 / Area,
    HabitatName=ifelse(Habitat=="CR","Coral Reef",ifelse(Habitat=="OR","Oyster Reef",ifelse(Habitat=="SM","Salt Marsh","Mangrove")))
  )

# Variable Definitions
wvdata<-subset(fielddata,fielddata$Hs.i!=-999)
reefdata<-subset(wvdata, wvdata$HabitatGrp=="CR")
vegdata<-subset(wvdata, wvdata$HabitatGrp!="CR")
itvegdata<-subset(vegdata, vegdata$Habitat!="SG") # much better than grep() for subsetting! Note: Can use use %in% for multiple strings
mgdata<-subset(vegdata, vegdata$Habitat=="MG")
smdata<-subset(vegdata, vegdata$Habitat=="SM")
sgdata<-subset(vegdata, vegdata$Habitat=="SG")
mgdata_submerge<-subset(mgdata, mgdata$hvbyh<=1)
mgdata_emerge<-subset(mgdata, mgdata$hvbyh>1)
smdata_submerge<-subset(smdata, smdata$hvbyh<1)
smdata_emerge<-subset(smdata, smdata$hvbyh>1)
reefdata_noXoutlier<-subset(reefdata, reefdata$BbyL<50)
mgdata_noXoutlier<-subset(mgdata, mgdata$BbyL<50)

#Estimating median widths
reefdata_Border<-reefdata[order(reefdata$B),]
reef_B_median<-median(reefdata_Border$B,na.rm=TRUE)
mgdata_Border<-mgdata[order(mgdata$B),]
mg_B_median<-median(mgdata_Border$B,na.rm=TRUE)
smdata_Border<-smdata[order(smdata$B),]
sm_B_median<-median(smdata_Border$B,na.rm=TRUE)
sgdata_Border<-sgdata[order(sgdata$B),]
sg_B_median<-median(sgdata_Border$B,na.rm=TRUE)

################ CREATE AND WRITE OUT SUPPLEMENTARY INFO TABLES ######################
#Note: except for first two, all tables to be written after generating in code below.

suppltable_waves<-fielddata[c("ID","First_Author","Year","Country","Habitat","Latitude","Longitude","Measurement.Type","n","Variable","Depth_Control","h","Tp","B","hv","Hs.i","Hs.iERR","Hs.t","Hs.tErr","R","Source_Text")]
suppltable_costs<-projectdata[c("ID","First_Author","Year","Country","Habitat","Latitude","Longitude","Width","Length","Area","Objective","Exposure","Habitat.Cost.Reported.Total","Cost.Year","Currency","CPI.2014","USD.2014","Habitat.Cost.USD2014","Habitat.Cost.Per.SquareMeter","Reported.Benefit.Type")]
write.table(suppltable_waves, "C:/Rfiles/SNAP/wavereductionmeasurements.txt", sep="\t",row.names=F)
write.table(suppltable_costs, "C:/Rfiles/SNAP/projectcostsandbenefits.txt", sep="\t",row.names=F)

write.table(nearbyprojects, "C:/Rfiles/SNAP/nearbyprojects.txt", sep="\t",row.names=F)
write.table(struct_comparison, "C:/Rfiles/SNAP/StructureToProjectCostRatios.txt", sep="\t",row.names=F)
write.table(meanRRvals, "C:/Rfiles/SNAP/MeanRRVals.txt", sep="\t",row.names=F)
write.table(struceffcostvals, "C:/Rfiles/SNAP/StruceffcostVals.txt", sep="\t",row.names=F)
write.table(plotdata, "C:/Rfiles/SNAP/Plotdata.txt", sep="\t",row.names=F)
write.table(benprojectdata, "C:/Rfiles/SNAP/beneficialprojects.txt", sep="\t",row.names=F)


###############  Random Effects Model for Wave Height Reduction Meta-Analyses #####################
#calculated using Eq 4.19, Eq 12.6, Eq 11.2 in Borenstein, 2009.

meta_analyses <- function(habitat) {
  habitat %>%
    mutate(
      RR=Hs.t/Hs.i,
      Yi=log(RR), #Effect size
      Vyi=(Hs.iErr^2)/(n*Hs.i^2)+(Hs.tErr^2)/(n*Hs.t^2), # Eq 4.19 for Swithin from Borenstein
      Wi=1/Vyi,
      Wi2=Wi^2,
      Wi.Yi=Wi*Yi,
      Wi.Yi2=Wi*Yi^2,
      df=nrow(reefdata)-1,
      Q=sum(Wi.Yi2)-((sum(Wi.Yi))^2/sum(Wi)),
      C=sum(Wi)-sum(Wi2)/sum(Wi),
      T2=(Q-df)/C,
      Vyi_star=Vyi+T2,
      Wi_star=1/Vyi_star,
      Wi_star.Yi=Wi_star*Yi,
      r=(1-RR)/B,
      M_star=sum(Wi_star.Yi)/sum(Wi_star),
      Vm_star=1/sum(Wi_star),
      Sem_star=sqrt(Vm_star),
      Z_star=M_star/Sem_star,
      MeanRR=100*(1-exp(M_star)),
      LLm=100*(1-exp(M_star-1.96*Sem_star)),
      ULm=100*(1-exp(M_star+1.96*Sem_star)),
      MeanRR_energy=(1-(1-MeanRR/100)^2)*100,
      LLm_energy=(1-(1-LLm/100)^2)*100,
      ULm_energy=(1-(1-ULm/100)^2)*100
    )
}

meta_reefs<-meta_analyses(reefdata)
meta_marshes<-meta_analyses(smdata)
meta_mangroves<-meta_analyses(mgdata)
meta_grasses<-meta_analyses(sgdata)

meanRRvals=data.frame(matrix(NA,ncol=8,nrow=4))
colnames(meanRRvals)<-c("Habitat","Mean Wave Height Reduction %","95% CI Lower", "95% CI Upper","Mean Wave Energy Reduction %","95% CI Lower (Energy)","95% CI Upper (Energy)","Average Measured Wave Height (m)")
meanRRvals[1:4,1]=c("Coral Reefs","Salt_Marshes","Mangroves","Seagrass/Kelp Beds")
meanRRvals[1:4,2]=c(mean(meta_reefs$MeanRR),mean(meta_marshes$MeanRR),mean(meta_mangroves$MeanRR),mean(meta_grasses$MeanRR))
meanRRvals[1:4,3]=c(mean(meta_reefs$ULm),mean(meta_marshes$ULm),mean(meta_mangroves$ULm),mean(meta_grasses$ULm))
meanRRvals[1:4,4]=c(mean(meta_reefs$LLm),mean(meta_marshes$LLm),mean(meta_mangroves$LLm),mean(meta_grasses$LLm))
meanRRvals[1:4,5]=c(mean(meta_reefs$MeanRR_energy),mean(meta_marshes$MeanRR_energy),mean(meta_mangroves$MeanRR_energy),mean(meta_grasses$MeanRR_energy))
meanRRvals[1:4,6]=c(mean(meta_reefs$ULm_energy),mean(meta_marshes$ULm_energy),mean(meta_mangroves$ULm_energy),mean(meta_grasses$ULm_energy))
meanRRvals[1:4,7]=c(mean(meta_reefs$LLm_energy),mean(meta_marshes$LLm_energy),mean(meta_mangroves$LLm_energy),mean(meta_grasses$LLm_energy))
meanRRvals[1:4,8]=c(mean(reefdata$Hs.i),mean(smdata$Hs.i),mean(mgdata$Hs.i),mean(sgdata$Hs.i))

############# PROJECT DATA CALCULATIONS #################

projdata<-subset(projectdata, projectdata$Latitude!="NaN")
benprojectdata<-subset(projectdata, Reported.Benefit.Type!="NaN")
creationprojectdata<-subset(projectdata,projectdata$Habitat.Cost.Per.SquareMeter!="NaN")

creationprojectdata<-group_by(creationprojectdata, HabitatName)
meanvals<-creationprojectdata %>%
  mutate(
    MeanUnitCost=mean(Habitat.Cost.Per.SquareMeter),
    Cost95Upper=MeanUnitCost+sd(Habitat.Cost.Per.SquareMeter)*1.96/sqrt(length(HabitatName)),
    Cost95Lower=MeanUnitCost-sd(Habitat.Cost.Per.SquareMeter)*1.96/sqrt(length(HabitatName)),
    MeanArea=mean(Area),
    Area95Upper=MeanArea+sd(Area)*1.96/sqrt(length(HabitatName)),
    Area95Lower=MeanArea-sd(Area)*1.96/sqrt(length(HabitatName))
    )
creationprojectdata<-ungroup(creationprojectdata)
meanvals<-meanvals %>%
  dplyr::select(HabitatName, MeanUnitCost, Cost95Upper, Cost95Lower, MeanArea, Area95Upper, Area95Lower) %>%
  distinct(HabitatName)

############ Replacement Cost Ratio Calculations ######################

# OBTAINING NEARBY PROJECT SITES FOR MAPPING THAT ARE < BbyL DECIMAL DEGREES FROM A FIELD EXPERIMENT SITE (0.46=50km; 0.09=10km; 0.01=1km).
for(i in 1:length(benprojectdata$ID))
{
  for (j in 1:length(fielddata$ID))
  {
    if((abs(fielddata$Latitude[j]-benprojectdata$Latitude[i])<=0.46)&&(abs(fielddata$Longitude[j]-benprojectdata$Longitude[i])<=0.46)&&(fielddata$R[j]>0))
    {
      benprojectdata$NearbyExp[i]="Yes"
      break
    }
    
    
    else
      benprojectdata$NearbyExp[i]="No"
    
  }
}
# OBTAINING NBD CONDITIONS FROM MEASUREMENTS AND NEARBY PROJECTS AND CALCULATING REPLACEMENT COST RATIOS
nearbyprojects<-data.frame(matrix(NA,ncol=23,nrow=length(fielddata$ID)))
colnames(nearbyprojects)<-c("Proj_ID","Exp_ID","Proj_Study","Exp_Study","Habitat","Country","Proj_lat","Proj_lon","Exp_lat","Exp_lon","Proj_area","Proj_width","Exp_width","Exp_depth","Exp_Hs","Exp_Tp","Exp_R","Exp_R_per_m","Exp_r","Proj_cost_m2","Proj_B_type","Proj_Rpercent","Proj_cost_per_m2")
k=0
for(i in 1:length(projdata$ID)) {
  for (j in 1:length(fielddata$ID)) {
    deltaLat<-abs(fielddata$Latitude[j]-projdata$Latitude[i])* pi/180
    deltaLon<-abs(fielddata$Longitude[j]-projdata$Longitude[i])* pi/180
    x=deltaLon*cos((fielddata$Latitude[j]*pi/180+projdata$Latitude[i]*pi/180)/2)
    y=deltaLat
    distance=6371*sqrt(x*x+y*y)
    if((distance<=50)&&(fielddata$R[j]>0))
    {
      k=k+1
      nearbyprojects$Proj_ID[k]<-projdata$ID[i]
      nearbyprojects$Exp_ID[k]<-fielddata$ID[j]
      nearbyprojects$Proj_Study[k]<-as.character(projdata$First_Author[i])
      nearbyprojects$Exp_Study[k]<-as.character(fielddata$First_Author[j])
      nearbyprojects$Habitat[k]<-as.character(projdata$HabitatName[i])
      nearbyprojects$Country[k]<-as.character(projdata$Country[i])
      nearbyprojects$Proj_lat[k]<-as.numeric(projdata$Latitude[i])
      nearbyprojects$Proj_lon[k]<-as.numeric(projdata$Longitude[i])
      nearbyprojects$Exp_lat[k]<-as.numeric(fielddata$Latitude[j])
      nearbyprojects$Exp_lon[k]<-as.numeric(fielddata$Longitude[j])
      nearbyprojects$Proj_area[k]<-as.numeric(projdata$Area[i])
      nearbyprojects$Proj_width[k]<-as.numeric(projdata$Width[i])
      nearbyprojects$Exp_width[k]<-as.numeric(fielddata$B[j])
      nearbyprojects$Exp_depth[k]<-as.numeric(fielddata$h[j])
      nearbyprojects$Exp_Hs[k]<-as.numeric(fielddata$Hs.i[j])
      nearbyprojects$Exp_Tp[k]<-as.numeric(fielddata$Tp[j])
      nearbyprojects$Exp_R[k]<-as.numeric(fielddata$R[j])
      nearbyprojects$Exp_R_per_m[k]<-as.numeric(fielddata$R[j]/fielddata$B[j])
      # nearbyprojects$Exp_r[k]<-as.numeric(fielddata$r[j]) # This r = R/B = r-actual/Hi
      nearbyprojects$Proj_cost_m2[k]<-as.numeric(projdata$Habitat.Cost.Per.SquareMeter[i])
      nearbyprojects$Proj_B_type[k]<-as.character(projdata$Reported.Benefit.Type[i])
    }
}}
# nearbyprojects<-subset(nearbyprojects,nearbyprojects$Proj_ID!="NA")
nearbyprojects<-nearbyprojects[order(nearbyprojects$Proj_ID, -nearbyprojects$Exp_R),]
# nearbyprojects<-nearbyprojects[!(duplicated(nearbyprojects$Proj_ID)&duplicated(nearbyprojects$Exp_lat)&duplicated(nearbyprojects$Exp_lon)),] #In-built function!

#CALCULATING COST AND WAVE REDUCTION OF PROJECTS USING NEARBY FIELD MEASUREMENTS
nearbyprojects<-filter(nearbyprojects,Proj_width>1&Proj_width!="NA")
nearbyprojects<-filter(nearbyprojects,Exp_Hs!=1)
nearbyprojects<-nearbyprojects %>%
  mutate(
    NbD_Kt=ifelse((1-Exp_r*Proj_width)>0.05,1-Exp_r*Proj_width,0.05),
    NbD_RPercent=(1-NbD_Kt^2)*100,
    Proj_Rpercent=ifelse((Exp_R*100*Proj_width/Exp_width)<100,Exp_R*100*Proj_width/Exp_width,100),
    Proj_cost_per_m=Proj_cost_m2*Proj_width
    )
nearbyprojects<-filter(nearbyprojects, Proj_cost_per_m!="NA")

#Calculating dimensions and cost of comparable breakwater for wave reduction
#PROCEDURE AND ASSUMPTIONS
#Costs and wave reduction per m length of coastline are calculated for a hypothetical nature-based defense project 
#based on Wave reduction from field measurements at the site and costs from a nearby project site (<10 km).
#Habitat width is assumed to be the same as the width for the nearby project.
#Total Wave reduction by a project is calculated as project width times reduction per m from the field measurement 
#Multiple measurements at the same location are treated as separate projects if there is variation in at least one
#of the following variables: Hs, project width, depth, rate of reduction or cost per m.
#The comparable structure is assumed to be a submerged breakwater (for wave reduction only)
# Structure freeboard (Rc) and crest width (Bc) are calculated to obtain a Kt roughly similar to the habitat Kt.
# Rc and Bc are calculated using standard engineeering formulae for the desired Kt from van der Meer et al. (2005).


## Replacement Cost Ratios and Structure Cost Curves
slope=1/100 #bottom slope
repstrucvol=36.288 #m3 per m of cross-section from NACCS Appendix T, Fig II-10
repstruccost=27600 #US$ per mfrom NACCS Appendix T, Table II-10
unitcost_NACCS=repstruccost/repstrucvol #cost per m3 per m for USA from NACCS Report Appendix T, Section II.4
struc_slope=1.5 #from NACCS Appendix T, Fig II-10

struct_comparison<-nearbyprojects[c("Proj_ID","Exp_ID","Country","Habitat","Proj_width","Exp_depth","Exp_Hs","Exp_Tp","Exp_R_per_m","Proj_Rpercent","Proj_cost_per_m")]
colnames(struct_comparison)<-c("Proj_ID","Exp_ID","Country","Habitat","Proj_width","Proj_depth","Hs","Tp","Proj_R_per_m","Proj_R_Total","Proj_cost_per_m")
struct_comparison<-struct_comparison %>% 
  mutate(
    Struc_depth=ifelse(Proj_depth!="NaN",Proj_depth,1),
    Proj_Kt=ifelse((1-Proj_R_per_m*Proj_width)>0,1-Proj_R_per_m*Proj_width,0.05),
    SSP=tan(slope)/sqrt(Hs*2*pi/(9.81*Tp*Tp)),
    Struc_Rc=as.numeric(c(-0.075,-0.021,-0.032,-0.55,-0.055,-0.12,-0.012,-0.012,-0.023,-0.014,-0.016,-0.021,-0.032)), #based on formulae for desired Kt
    Struc_Bc=as.numeric(c(2,4,4,4,6,2,2,2,6,6,2,4,4)), #based on formulae for desired Kt
    Struc_Kt=ifelse(Struc_Bc/Hs>12,(Struc_Kt=-0.35*Struc_Rc/Hs+0.51*(1-exp(-0.41*SSP))/((Struc_Bc/Hs)^0.65)),(Struc_Kt=-0.4*Struc_Rc/Hs+0.64*(1-exp(-0.41*SSP))/((Struc_Bc/Hs)^0.31))),
    Struc_Ht=Struc_depth+Struc_Rc,
    Struc_vol=(Struc_Ht*(Struc_Bc+Struc_Bc+Struc_Ht*struc_slope)/2),
    Struc_cost_per_m=ifelse((Country=="Vietnam"|Country=="China"),Struc_vol*unitcost_NACCS/10,Struc_vol*unitcost_NACCS),
    RCratio=Struc_cost_per_m/Proj_cost_per_m,
    Struc_max_effective_unitcost_per_m=Proj_cost_per_m/Struc_vol,
    BeneficialProj_cost_per_m=Proj_cost_per_m/RCratio,
    Kt_Family=as.numeric(ifelse(Proj_Kt==0.05,100,ifelse(Proj_Kt<0.25,75,ifelse(Proj_Kt<0.5,50,25))))
    )
struct_comparison<-struct_comparison %>% 
  group_by(Proj_width) %>%
  mutate(
    WidthAvg_RCratio=mean(RCratio)
    ) %>%
  ungroup() # Remember to ungroup!

Kt_All<-struct_comparison %>%
  distinct(Proj_depth,Proj_width) %>%
  select(Proj_ID, Country, Habitat, Proj_width, Proj_depth, Hs, Struc_Bc, Proj_cost_per_m,Proj_R_Total, Struc_max_effective_unitcost_per_m, Kt_Family) %>%
  arrange(Proj_depth)

Kt_Struc<-Kt_All %>%
  distinct(Proj_depth) %>%
  mutate(
    #Kt = 0.05, ## Rc Estimated for given depth and Hs, for Kt=0.05
    Rc0.05=as.numeric(c(-0.01,-0.015,-0.01,-0.012,-0.015,-0.015)), 
        Ht0.05=Proj_depth+Rc0.05,
    Vol0.05=(Ht0.05*(Struc_Bc+Struc_Bc+Ht0.05*struc_slope)/2),
    Struc_cost_0.05=Vol0.05*unitcost_NACCS,
    #Kt = 0.25,  ## Rc Estimated for given depth and Hs, for Kt=0.25
    Rc0.25=as.numeric(c(-0.14,-0.105,-0.16,-0.072,-0.115,-0.08)), 
    Ht0.25=Proj_depth+Rc0.25,
    Vol0.25=(Ht0.25*(Struc_Bc+Struc_Bc+Ht0.25*struc_slope)/2),
    Struc_cost_0.25=Vol0.25*unitcost_NACCS,
    #Kt = 0.5, ## Rc Estimated for given depth and Hs, for Kt=0.5
    Rc0.5=as.numeric(c(-0.3,-0.215,-0.35,-0.145,-0.24,-0.165)),
    Ht0.5=Proj_depth+Rc0.5,
    Vol0.5=(Ht0.5*(Struc_Bc+Struc_Bc+Ht0.5*struc_slope)/2),
    Struc_cost_0.5=Vol0.5*unitcost_NACCS,
    #Kt = 0.75, ## Rc Estimated for given depth and Hs, for Kt=0.75
    Rc0.75=as.numeric(c(-0.465,-0.33,-0.535,-0.22,-0.365,-0.247)),
    Ht0.75=Proj_depth+Rc0.75,
    Vol0.75=(Ht0.75*(Struc_Bc+Struc_Bc+Ht0.75*struc_slope)/2),
    Struc_cost_0.75=Vol0.75*unitcost_NACCS,
    #Kt = 0.05 & Hs=0.1 & Bc=2,  ## Rc Estimated for all depths, for Kt=0.05, Hs = 0.1
    Rc0.05_Hs0.1_Bc2=as.numeric(c(-0.012,-0.012,-0.012,-0.012,-0.012,-0.012)), 
    Ht0.05_Hs0.1_Bc2=Proj_depth+Rc0.05_Hs0.1_Bc2,
    Vol0.05_Hs0.1_Bc2=(Ht0.05_Hs0.1_Bc2*(2+2+Ht0.05_Hs0.1_Bc2*struc_slope)/2),
    StrucCost0.05_Hs0.1_Bc2=Vol0.05_Hs0.1_Bc2*unitcost_NACCS,
    #Kt = 0.75 & Hs=0.1 & Bc=2,  ## Rc Estimated for all depths, for Kt=0.75, Hs = 0.1
    Rc0.75_Hs0.1_Bc2=as.numeric(c(-0.213,-0.213,-0.213,-0.213,-0.213,-0.213)), 
    Ht0.75_Hs0.1_Bc2=Proj_depth+Rc0.75_Hs0.1_Bc2,
    Vol0.75_Hs0.1_Bc2=(Ht0.75_Hs0.1_Bc2*(2+2+Ht0.75_Hs0.1_Bc2*struc_slope)/2),
    StrucCost0.75_Hs0.1_Bc2=Vol0.75_Hs0.1_Bc2*unitcost_NACCS,
    #Kt = 0.05 & Hs=0.2 & Bc=2,  ## Rc Estimated for all depths, for Kt=0.05, Hs = 0.2
    Rc0.05_Hs0.2_Bc2=as.numeric(c(-0.015,-0.015,-0.015,-0.015,-0.015,-0.015)), 
    Ht0.05_Hs0.2_Bc2=Proj_depth+Rc0.05_Hs0.2_Bc2,
    Vol0.05_Hs0.2_Bc2=(Ht0.05_Hs0.2_Bc2*(2+2+Ht0.05_Hs0.2_Bc2*struc_slope)/2),
    StrucCost0.05_Hs0.2_Bc2=Vol0.05_Hs0.2_Bc2*unitcost_NACCS,
    #Kt = 0.75 & Hs=0.2 & Bc=2,  ## Rc Estimated for all depths, for Kt=0.75, Hs = 0.2
    Rc0.75_Hs0.2_Bc2=as.numeric(c(-0.39,-0.39,-0.39,-0.39,-0.39,-0.39)), 
    Ht0.75_Hs0.2_Bc2=Proj_depth+Rc0.75_Hs0.2_Bc2,
    Vol0.75_Hs0.2_Bc2=(Ht0.75_Hs0.2_Bc2*(2+2+Ht0.75_Hs0.2_Bc2*struc_slope)/2),
    StrucCost0.75_Hs0.2_Bc2=Vol0.75_Hs0.2_Bc2*unitcost_NACCS,
    #Kt = 0.05 & Hs=0.3 & Bc=2,  ## Rc Estimated for all depths, for Kt=0.05, Hs = 0.3
    Rc0.05_Hs0.3_Bc2=as.numeric(c(-0.01,-0.01,-0.01,-0.01,-0.01,-0.01)), 
    Ht0.05_Hs0.3_Bc2=Proj_depth+Rc0.05_Hs0.3_Bc2,
    Vol0.05_Hs0.3_Bc2=(Ht0.05_Hs0.3_Bc2*(2+2+Ht0.05_Hs0.3_Bc2*struc_slope)/2),
    StrucCost0.05_Hs0.3_Bc2=Vol0.05_Hs0.3_Bc2*unitcost_NACCS,
    #Kt = 0.75 & Hs=0.3 & Bc=2,  ## Rc Estimated for all depths, for Kt=0.75, Hs = 0.3
    Rc0.75_Hs0.3_Bc2=as.numeric(c(-0.53,-0.53,-0.53,-0.53,-0.53,-0.53)), 
    Ht0.75_Hs0.3_Bc2=Proj_depth+Rc0.75_Hs0.3_Bc2,
    Vol0.75_Hs0.3_Bc2=(Ht0.75_Hs0.3_Bc2*(2+2+Ht0.75_Hs0.3_Bc2*struc_slope)/2),
    StrucCost0.75_Hs0.3_Bc2=Vol0.75_Hs0.3_Bc2*unitcost_NACCS
    )
Kt_All<-Kt_All %>%
  arrange(-Kt_Family) %>%
  mutate(
    Kt_Family=as.character(Kt_Family)
    )
Kt_Mangrove<-filter(Kt_All,Habitat=="Mangrove")
Kt_SaltMarsh<-filter(Kt_All,Habitat!="Mangrove")

mg_struc_maxunitcost_mean<-mean(Kt_Mangrove$Struc_max_effective_unitcost_per_m)
sm_struc_maxunitcost_mean<-mean(Kt_SaltMarsh$Struc_max_effective_unitcost_per_m)
mg_struc_maxunitcost_sd<-sd(Kt_Mangrove$Struc_max_effective_unitcost_per_m)
sm_struc_maxunitcost_sd<-sd(Kt_SaltMarsh$Struc_max_effective_unitcost_per_m)
mg_struc_maxunitcost_margin<-qt(0.975,length(Kt_Mangrove)-1)*mg_struc_maxunitcost_sd/sqrt(length(Kt_Mangrove))
sm_struc_maxunitcost_margin<-qt(0.975,length(Kt_SaltMarsh)-1)*sm_struc_maxunitcost_sd/sqrt(length(Kt_SaltMarsh))

struceffcostvals=data.frame(matrix(NA,ncol=4,nrow=2))
colnames(struceffcostvals)<-c("Habitat","Structural Replacement value (per m3 per m)","95% CI Upper", "95% CI Lower")
struceffcostvals[1:2,1]=c("Mangrove","Salt-Marsh")
struceffcostvals[1:2,2]=c(mg_struc_maxunitcost_mean,sm_struc_maxunitcost_mean)
struceffcostvals[1:2,3]=c(mg_struc_maxunitcost_mean-mg_struc_maxunitcost_margin, sm_struc_maxunitcost_mean-sm_struc_maxunitcost_margin)
struceffcostvals[1:2,4]=c(mg_struc_maxunitcost_mean+mg_struc_maxunitcost_margin, sm_struc_maxunitcost_mean+sm_struc_maxunitcost_margin)
struceffcostvals

##########Cost Curve Plots ###################
costcurveplot_All<-ggplot()+scale_linetype_identity()+
  geom_line(data=Kt_Struc, aes(x=Kt_Struc$Proj_depth,y=Kt_Struc$StrucCost0.05_Hs0.2_Bc2,linetype="solid"),show.legend=TRUE)+
  geom_line(data=Kt_Struc, aes(x=Kt_Struc$Proj_depth,y=Kt_Struc$StrucCost0.75_Hs0.2_Bc2,linetype="1F"),show.legend=TRUE)+
  geom_line(data=Kt_Struc, aes(x=Kt_Struc$Proj_depth,y=Kt_Struc$StrucCost0.05_Hs0.2_Bc2/10,linetype="longdash"),show.legend=TRUE)+
  geom_line(data=Kt_Struc, aes(x=Kt_Struc$Proj_depth,y=Kt_Struc$StrucCost0.75_Hs0.2_Bc2/10,linetype="dotted"),show.legend=TRUE)+
  geom_point(data=Kt_All, position="jitter",aes(Kt_All$Proj_depth,Kt_All$Proj_cost_per_m, colour=Kt_All$Habitat,size=Kt_All$Proj_width),shape=1)+
  geom_text(data=Kt_All[c(3,5,11:13),],aes(Proj_depth,Proj_cost_per_m,label=paste(round(Proj_R_Total,digits=0),"%",sep="")),size=3,vjust=1)+
  scale_linetype_manual("Structure wave reduction %",values=c(5,3,2,1),labels=c("solid"="100(Europe/USA)","1F"="25(Europe/USA)","longdash"="100(Vietnam)","dotted"="25(Vietnam)"))+
  scale_size_continuous("Project Width",range=c(3,15),breaks=c(500,1500,2500,3000))+
  scale_colour_manual("Habitat",values=c("green","violet"),labels=c("Mangrove","Salt-Marsh"))+
  theme_bw()+labs(title="Replacement Structure Cost Curves",x="Depth (m)",y="Cost per m coastline (US $)")+
  guides(linetype=guide_legend(order=1))+
  guides(colour=guide_legend(order=2, override.aes=list(linetype=c(0,0)),fill="NA"))+
  guides(size=guide_legend(order=4,override.aes=list(linetype=c(0,0,0))))
costcurveplot_All

# costcurveplot<-function(costfactor, plotdata, textdata,breakdata,coord_xlim,coord_ylim,Title,override_colour_aes,override_size_aes)
# {
#   ggplot()+scale_linetype_identity()+
#     geom_line(data=Kt_Struc, aes(x=Kt_Struc$Proj_depth,y=costfactor*Kt_Struc$StrucCost0.05_Hs0.2_Bc2,linetype="solid"),show.legend=TRUE)+
#     geom_line(data=Kt_Struc, aes(x=Kt_Struc$Proj_depth,y=costfactor*Kt_Struc$StrucCost0.75_Hs0.2_Bc2*costfactor,linetype="1F"),show.legend=TRUE)+
#     geom_line(data=Kt_Struc, aes(x=Kt_Struc$Proj_depth,y=costfactor*Kt_Struc$StrucCost0.05_Hs0.3_Bc2*costfactor,linetype="longdash"),show.legend=TRUE)+
#     geom_line(data=Kt_Struc, aes(x=Kt_Struc$Proj_depth,y=costfactor*Kt_Struc$StrucCost0.75_Hs0.3_Bc2*costfactor,linetype="dotted"),show.legend=TRUE)+
#     geom_line(data=Kt_Struc, aes(x=Kt_Struc$Proj_depth,y=costfactor*Kt_Struc$StrucCost0.05_Hs0.1_Bc2*costfactor,linetype="longdash"),show.legend=TRUE)+
#     geom_line(data=Kt_Struc, aes(x=Kt_Struc$Proj_depth,y=costfactor*Kt_Struc$StrucCost0.75_Hs0.1_Bc2*costfactor,linetype="dotted"),show.legend=TRUE)+  
#     geom_point(data=plotdata, aes(Proj_depth,Proj_cost_per_m, size=Proj_width),shape=1)+
#     geom_text(data=textdata,aes(Proj_depth,Proj_cost_per_m,label=paste(round(Proj_R_Total,digits=0),"%",sep="")),size=3,hjust=1,vjust=0)+
#     scale_linetype_manual("Structure wave reduction %",values=c(5,3,2,1),labels=c("solid"="100(Hs=0.2)","1F"="25(Hs=0.2)","longdash"="100(Hs=0.3/0.1)","dotted"="25(Hs=0.3/0.1)"))+
#     scale_size_continuous("Project Width",range=c(5,11),breaks=breakdata)+
#     coord_cartesian(xlim=coord_xlim,ylim=coord_ylim)+
#     theme_bw()+labs(title=Title,x="Depth (m)",y="Cost per m coastline (US $)")+
#     guides(linetype=guide_legend(order=1))+
#     guides(colour=guide_legend(order=2, override.aes=override_colour_aes,fill="NA"))+
#     guides(size=guide_legend(order=3,override.aes=override_size_aes))
# }
# 
# costcurveplot_Mangrove<-costcurveplot(costfactor=0.1,plotdata=Kt_Mangrove,textdata=subset(Kt_Mangrove[c(1:2,4,6:7),]),breakdata=c(800,1200,1500,1600),
#                                       coord_xlim=c(1,1.75),coord_ylim=c(0,1000),Title="Mangrove Habitat and Replacement Structure Costs",
#                                       override_colour_aes=list(linetype=c(0)),override_size_aes=list(linetype=c(0,0,0)))
# 
# costcurveplot_SaltMarsh<-costcurveplot(costfactor=1,plotdata=Kt_SaltMarsh,textdata=subset(Kt_SaltMarsh[4:6,]),breakdata=c(100,400,800,1600,2800,3000),
#                                        coord_xlim=c(1,4),coord_ylim=c(0,15000),Title="Salt-Marsh Habitat and Replacement Structure Costs",
#                                        override_colour_aes=list(linetype=c(0)),override_size_aes=list(linetype=c(0,0,0,0,0)))

costcurveplot_SaltMarsh<-ggplot()+scale_linetype_identity()+
  geom_line(data=Kt_Struc, aes(x=Kt_Struc$Proj_depth,y=Kt_Struc$StrucCost0.05_Hs0.2_Bc2,linetype="solid"),show.legend=TRUE)+
  geom_line(data=Kt_Struc, aes(x=Kt_Struc$Proj_depth, y=Kt_Struc$StrucCost0.75_Hs0.2_Bc2,linetype="1F"),show.legend=TRUE)+
  geom_line(data=Kt_Struc, aes(x=Kt_Struc$Proj_depth,y=Kt_Struc$StrucCost0.05_Hs0.3_Bc2,linetype="longdash"),show.legend=TRUE)+
  geom_line(data=Kt_Struc, aes(x=Kt_Struc$Proj_depth, y=Kt_Struc$StrucCost0.75_Hs0.3_Bc2,linetype="dotted"),show.legend=TRUE)+
  geom_line(data=Kt_Struc, aes(x=Kt_Struc$Proj_depth,y=Kt_Struc$StrucCost0.05_Hs0.1_Bc2,linetype="longdash"),show.legend=TRUE)+
  geom_line(data=Kt_Struc, aes(x=Kt_Struc$Proj_depth, y=Kt_Struc$StrucCost0.75_Hs0.1_Bc2,linetype="dotted"),show.legend=TRUE)+
  geom_point(data=Kt_SaltMarsh, aes(Kt_SaltMarsh$Proj_depth,Kt_SaltMarsh$Proj_cost_per_m, size=Kt_SaltMarsh$Proj_width),shape=1)+
  geom_text_repel(data=Kt_SaltMarsh,aes(Proj_depth,Proj_cost_per_m,label=paste(round(Proj_R_Total,digits=0),"%",sep="")),size=3)+
  scale_linetype_manual("Structure wave reduction %",values=c(5,3,2,1),labels=c("solid"="100(Hs=0.2)","1F"="25(Hs=0.2)","longdash"="100(Hs=0.3/0.1)","dotted"="25(Hs=0.3/0.1)"))+
  scale_size_continuous("Project Width",range=c(5,11),breaks=c(100,400,800,1600,2800,3000))+
  coord_cartesian(xlim=c(0,4),ylim=c(0,15000))+
  theme_bw()+labs(title="Salt-Marsh Habitat and Replacement Structure Costs",x="Depth (m)",y="Cost per m coastline (US $)")+
  guides(linetype=guide_legend(order=1))+
  guides(colour=guide_legend(order=2, override.aes=list(linetype=c(0)),fill="NA"))+
  guides(size=guide_legend(order=3,override.aes=list(linetype=c(0,0,0,0,0))))
costcurveplot_SaltMarsh

costcurveplot_Mangrove<-ggplot()+scale_linetype_identity()+
  geom_line(data=Kt_Struc, aes(x=Kt_Struc$Proj_depth,y=Kt_Struc$StrucCost0.05_Hs0.2_Bc2/10,linetype="solid"),show.legend=TRUE)+
  geom_line(data=Kt_Struc, aes(x=Kt_Struc$Proj_depth, y=Kt_Struc$StrucCost0.75_Hs0.2_Bc2/10,linetype="1F"),show.legend=TRUE)+
  geom_line(data=Kt_Struc, aes(x=Kt_Struc$Proj_depth,y=Kt_Struc$StrucCost0.05_Hs0.3_Bc2/10,linetype="longdash"),show.legend=TRUE)+
  geom_line(data=Kt_Struc, aes(x=Kt_Struc$Proj_depth, y=Kt_Struc$StrucCost0.75_Hs0.3_Bc2/10,linetype="dotted"),show.legend=TRUE)+
  geom_line(data=Kt_Struc, aes(x=Kt_Struc$Proj_depth,y=Kt_Struc$StrucCost0.05_Hs0.1_Bc2/10,linetype="longdash"),show.legend=TRUE)+
  geom_line(data=Kt_Struc, aes(x=Kt_Struc$Proj_depth, y=Kt_Struc$StrucCost0.75_Hs0.1_Bc2/10,linetype="dotted"),show.legend=TRUE)+  
  geom_point(data=Kt_Mangrove, aes(Kt_Mangrove$Proj_depth,Kt_Mangrove$Proj_cost_per_m, size=Kt_Mangrove$Proj_width),shape=1)+
  geom_text_repel(data=Kt_Mangrove,aes(Proj_depth,Proj_cost_per_m,label=paste(round(Proj_R_Total,digits=0),"%",sep="")),size=3)+
  scale_linetype_manual("Structure wave reduction %",values=c(5,3,2,1),labels=c("solid"="100(Hs=0.2)","1F"="25(Hs=0.2)","longdash"="100(Hs=0.3/0.1)","dotted"="25(Hs=0.3/0.1)"))+
  scale_size_continuous("Project Width",range=c(5,11),breaks=c(800,1200,1500,1600))+
  coord_cartesian(xlim=c(0.5,2),ylim=c(0,1000))+
  theme_bw()+labs(title="Mangrove Habitat and Replacement Structure Costs",x="Depth (m)",y="Cost per m coastline (US $)")+
  guides(linetype=guide_legend(order=1))+
  guides(colour=guide_legend(order=2, override.aes=list(linetype=c(0)),fill="NA"))+
  guides(size=guide_legend(order=3,override.aes=list(linetype=c(0,0,0))))
costcurveplot_Mangrove

grid.arrange(costcurveplot_Mangrove, costcurveplot_SaltMarsh,nrow=2)

#######################  WAVE ATTENUATION GLOBAL PLOTS ########################
## THIS FIGURE IS NOW GENERATED IN ARCGIS USING ARCGIS BASEMAP LAYERS.

plotdata<-fielddata
# (using http://stackoverflow.com/questions/16028659/plots-on-a-map-using-ggplot2)
map<-ggplot()
world_map<-map_data(map="world")
world_map<-subset(world_map, world_map$region!="Antarctica")
p<-ggplot()+coord_fixed()
base_world<-p+geom_path(data=world_map, aes(x=long, y=lat, group=group),alpha=0.8)
##CAUTION: TOO MANY +S WHEN COPY-PASTING FROM R FILE RESULTS IN 'INVALID UNARY OPERATOR' ERROR
map<-base_world+geom_point(aes(plotdata$Longitude,plotdata$Latitude, colour=plotdata$Habitat, size=plotdata$R*100), shape=1, alpha=0.8)+
  geom_text(size=2, label="Habitat")+labs(x="",y="Latitude", title="Field Measurements of Wave Reduction in Coastal Habitats")+
  scale_colour_manual("Habitat Type", labels=c("Reef Crest", "Reef Flat", "Whole Reef","Mangrove","Seagrass/Kelp","Salt Marsh"), values=c("red","orange","yellow2","green","blue","violet"))+
  theme_bw()+scale_size_continuous("Wave Reduction %               ") # scale_** changes legend titles
map

map_benefits<-base_world+geom_point(data=benprojectdata, aes(x=Longitude, y=Latitude, colour=benprojectdata$Habitat, shape=benprojectdata$NearbyExp), size=4,alpha=0.8)+
  labs(x="Longitude", y="Latitude", title="Nature-based Defense Projects Reporting Coastal Protection Benefits")+
  scale_shape_manual("Nearby Field Measurement", values=c(1, 3), labels=c("No","Yes"))+
  theme_bw()+scale_colour_manual("Habitat Type", labels=c("Coral Reef","Mangrove","Salt Marsh"), values=c("orange","green","violet"))
map_benefits

grid.arrange(map, map_benefits, nrow=2)

################### TREND ANALYSES PLOTS ###############################

## Hs Reduction Versus Width Plots
#Reefs
nls_reef_B<-nls(R~a*log(B)+b,reefdata,start=list(a=0.1,b=0.01))
nls_R2_reef_B<-1-sum(residuals(nls_reef_B)^2)/sum((reefdata$R - mean(reefdata$R))^2)
CRBplot<-ggplot(data=reefdata,group=Habitat)+geom_point(aes(x=B,y=R,shape=Habitat),size=3)+ylim(-0.2,1.5)+
  labs(title="Reefs",x="",y="Wave Reduction %")+scale_shape_manual("Reef Environment",values=c(2,4,3),labels=c("Crest","Flat","Whole Reef"))+
  theme_bw()+theme(legend.position=c(1,0), legend.justification=c(1,0))+
  geom_smooth(data=reefdata,aes(x=B,y=R),method="nls",formula=y~I(a*log(x)+b),start=list(a=0.1,b=0.01),se=FALSE,linetype=2,colour="black")+
  scale_y_continuous(labels=function(x)x*100)+annotate('text',x=2800,y=0.8,label=paste("R^2==",format(nls_R2_reef_B, digits=4)),parse=TRUE)
CRBplot

#Others
Width_plot<-function(plotdata,Title,Shapedata)
{
  ggplot(data=plotdata)+geom_point(aes(x=B,y=R),shape=Shapedata)+
    labs(title=Title,x="Habitat Width (m)",y="Wave Reduction %")+
    theme_bw()+scale_y_continuous(labels=function(x)x*100)
}
MGBplot<-Width_plot(mgdata,"Mangroves",0)
SMBplot<-Width_plot(smdata,"Salt-Marshes",5)
SGBplot<-Width_plot(sgdata,"Seagrass/Kelp",8)

grid.arrange(CRBplot,MGBplot,SMBplot,SGBplot,ncol=2)

## Hs Reduction Versus Incoming Hs Plots
#Variable Transformations For Absolute and Percentage Reduction Extents vs Incoming Hs (Excluding All Hs.i = 1)
wvdata<-subset(fielddata,fielddata$Hs.i!=1)

#Absolute Hs Plots
# Reefs
reef_Hs.i_lm_Abs<-lm((reefdata$Hs.i - reefdata$Hs.t)~reefdata$Hs.i)
CRHsplot_Abs<-ggplot(data=reefdata,group=Habitat)+geom_point(aes(x=Hs.i,y=(Hs.i-Hs.t),shape=Habitat))+
  labs(title="Reefs",x="Incoming Wave Height, Hs.i (m)",y="Absolute Wave Attenuation (Hs.i - Hs.t) (m)")+scale_shape_manual("Reef Environment",values=c(2,4,3),labels=c("Crest","Flat","Whole Reef"))+
  theme_bw()+theme(legend.position=c(0,1), legend.justification=c(0,1))+
  geom_smooth(data=reefdata,aes(x=Hs.i,y=Hs.i-Hs.t),method=lm,formula=y~x,se=FALSE,linetype=2,colour="black")+
  annotate('text',x=1.05,y=1.7,label=paste("R^2==",format(summary(reef_Hs.i_lm_Abs)$adj.r.squared, digits=4)),parse=TRUE)+
  annotate('text',x=1.65,y=1.67,label=paste("slope==",format(summary(reef_Hs.i_lm_Abs)$coefficients[2,1], digits=4)),parse=TRUE)

#Others
Hs.i_abs_plot<-function(plotdata,Title,Shapedata,r2_x,r2_y,slope_x,slope_y,lmdata)
{
  ggplot(data=plotdata)+geom_point(aes(x=Hs.i,y=(Hs.i-Hs.t)),shape=Shapedata)+
    labs(title=Title,x="Incoming Wave Height, Hs.i (m)",y="Absolute Wave Attenuation (Hs.i - Hs.t) (m)")+
    theme_bw()+geom_smooth(data=plotdata,aes(x=Hs.i,y=Hs.i-Hs.t),method=lm,formula=y~x,se=FALSE,linetype=2,colour="black")+
    annotate('text',x=r2_x,y=r2_y,label=paste("R^2==",format(summary(mg_Hs.i_lm_Abs)$adj.r.squared, digits=4)),parse=TRUE)+
    annotate('text',x=slope_x,y=slope_y,label=paste("slope==",format(summary(lmdata)$coefficients[2,1], digits=4)),parse=TRUE)
  
}
mg_Hs.i_lm_Abs<-lm((mgdata$Hs.i - mgdata$Hs.t)~mgdata$Hs.i)
sm_Hs.i_lm_Abs<-lm((smdata$Hs.i - smdata$Hs.t)~smdata$Hs.i)
sg_Hs.i_lm_Abs<-lm((sgdata$Hs.i - sgdata$Hs.t)~sgdata$Hs.i)
MGHsplot_Abs<-Hs.i_abs_plot(mgdata,"Mangroves",0,0.12,0.31,0.24,0.305,mg_Hs.i_lm_Abs)
SMHsplot_Abs<-Hs.i_abs_plot(smdata,"Salt-Marshes",5,0.14,0.62,0.44,0.61,sm_Hs.i_lm_Abs)
SGHsplot_Abs<-Hs.i_abs_plot(sgdata,"Seagrass/Kelp",8,0.15,0.46,0.45,0.45,sg_Hs.i_lm_Abs)

# Relative Hs Plots
#Reefs
reef_Hs.i_lm<-lm(R~Hs.i, reefdata)
nls_reef_H<-nls(R~a*log(Hs.i)+b,reefdata,start=list(a=0.1,b=0.01))
nls_R2_reef_H<-1-sum(residuals(nls_reef_H)^2)/sum((reefdata$R - mean(reefdata$R))^2)
CRHsplot<-ggplot(data=reefdata,group=Habitat)+geom_point(aes(x=Hs.i,y=R,shape=Habitat))+
  labs(title="Reefs",x="Incoming Wave Height, Hs.i (m)",y="Wave Reduction %")+scale_shape_manual("Reef Environment",values=c(2,4,3),labels=c("Crest","Flat","Whole Reef"))+
  theme_bw()+theme(legend.position=c(1,0), legend.justification=c(1,0))+
  geom_smooth(data=reefdata,aes(x=Hs.i,y=R),method="nls",formula=y~I(a*log(x)+b),se=FALSE,start=list(a=0.1,b=0.01),linetype=2,colour="black")

#Others
Hs.i_rel_plot<-function(plotdata,Title,Shapedata)
{
  ggplot(data=plotdata)+geom_point(aes(x=Hs.i,y=R),shape=Shapedata)+
    labs(title=Title,x="Incoming Wave Height, Hs.i (m)",y="Wave Reduction %")+
    theme_bw()+scale_y_continuous(labels=function(x)x*100)
}
MGHsplot<-Hs.i_rel_plot(mgdata,"Mangroves",0)
SMHsplot<-Hs.i_rel_plot(smdata,"Salt-Marshes",5)
SGHsplot<-Hs.i_rel_plot(sgdata,"Salt-Marshes",8)

grid.arrange(CRHsplot_Abs,MGHsplot_Abs,SMHsplot_Abs,SGHsplot_Abs,ncol=2)
grid.arrange(CRHsplot,MGHsplot,SMHsplot,SGHsplot,ncol=2)

## DIMENSIONLESS PARAMETER PLOTS BY HABITAT

#B/L
wvdata<-subset(fielddata,fielddata$Hs.i!=-999)
#Reefs
reefdata<-wvdata[grep("CR",wvdata$Habitat),] # grep() subsets data frame by rows based on a char variable
MVRdata_R_reef_all<-reefdata[c("Habitat","Hs.i","Hs.t","R","BbyL")]
MVRdata_R_reef_all<-MVRdata_R_reef_all[complete.cases(MVRdata_R_reef_all),]
MVRdata_R_reef_out<-subset(MVRdata_R_reef_all, MVRdata_R_reef_all$BbyL>50)
MVRdata_R_reef<-reefdata_noXoutlier[c("Habitat","Hs.i","Hs.t","R","BbyL")]
MVRdata_R_reef<-MVRdata_R_reef[complete.cases(MVRdata_R_reef),]
nls_reef_BbyL<-nls(R~a*log(BbyL)+b,MVRdata_R_reef,start=list(a=0.1,b=0.01))
nls_R2_reef_BbyL<-1-sum(residuals(nls_reef_BbyL)^2)/sum(MVRdata_R_reef$R - mean(MVRdata_R_reef$R)^2)
sub<-ggplot(data=MVRdata_R_reef_all, group=Habitat)+geom_point(aes(x=BbyL,y=R,shape=Habitat))+
  geom_rect(data=MVRdata_R_reef_all$BbyL, xmin=0, xmax=20,ymin=-0.2, ymax=1,fill="light blue",hvbyh=0.7)+
  theme_bw()+theme(legend.position="none")+labs(x="B/L",y="Wave Reduction %")
sub$layers<-rev(sub$layers)
REEFBbyLplot<-ggplot(data=MVRdata_R_reef, group=Habitat)+geom_point(aes(x=BbyL,y=R,shape=Habitat))+
  labs(title="Reefs",x="Relative Width, B/L",y="Wave Reduction %")+scale_shape_manual("Reef Environment",values=c(2,4,3),labels=c("Crest","Flat","Whole Reef"))+
  theme_bw()+theme(legend.position=c(1,0),legend.justification=c(1,0))+
  geom_smooth(data=MVRdata_R_reef,aes(x=BbyL,y=R),method=nls,formula=y~I(a*log(x)+b),start=list(a=0.1,b=0.01),se=FALSE,linetype=2,colour="black")+
  annotate('text',x=8,y=0.6,label=paste("R^2==",format(nls_R2_reef_BbyL, digits=4)),parse=TRUE)+
  annotation_custom(ggplotGrob(sub),xmin=3,xmax=10,ymin=-0.2,ymax=0.3)

grid.arrange(REEFHbyhplot,REEFBbyLplot,ncol=2)

#Others
BbyLplot<-function(plotdata,Shapedata,Title)
{
  ggplot(data=plotdata)+geom_point(aes(x=BbyL,y=R),shape=Shapedata)+
    labs(title=Title,x="Relative Width, B/L",y="Wave Attenuation %")+theme_bw()
}

MGBbyLplot<-BbyLplot(mgdata,0,"Mangroves")
SMBbyLplot<-BbyLplot(smdata,5,"Salt-Marshes")
SGBbyLplot<-BbyLplot(sgdata,8,"Seagrass/Kelp")

grid.arrange(REEFBbyLplot, MGBbyLplot,SMBbyLplot,SGBbyLplot,ncol=2)

#H/h
wvdata<-subset(fielddata,fielddata$Hs.i!=1)
#Reefs
reefdata<-wvdata[grep("CR",wvdata$Habitat),] # grep() subsets data frame by rows based on a char variable
MVRdata_R_REEF<-reefdata[c("Habitat","Hs.i","Hs.t","R", "Hbyh")]
MVRdata_R_REEF<-MVRdata_R_REEF[complete.cases(MVRdata_R_REEF),]
nls_reef_Hbyh<-nls(R~a*log(Hbyh)+b,MVRdata_R_REEF,start=list(a=0.1,b=0.01))
nls_R2_reef_Hbyh<-1-sum(residuals(nls_reef_Hbyh)^2)/sum(MVRdata_R_REEF$R - mean(MVRdata_R_REEF$R)^2)
REEFHbyhplot<-ggplot(data=MVRdata_R_REEF, group=Habitat)+geom_point(aes(x=Hbyh,y=R,shape=Habitat))+
  labs(title="Reefs",x="Relative wave height, H/h",y="Wave Attenuation %")+scale_shape_manual("Reef Environment",values=c(2,4,3),labels=c("Crest","Flat","Whole Reef"))+
  theme_bw()+theme(legend.position=c(1,0),legend.justification=c(1,0))+
  geom_smooth(data=MVRdata_R_REEF,aes(x=Hbyh,y=R),method=nls,formula=y~I(a*log(x)+b),se=FALSE,start=list(a=0.1,b=0.01), linetype=2,colour="black")+
  annotate('text',x=1.5,y=0.6,label=paste("R^2==",format(nls_R2_reef_Hbyh, digits=4)),parse=TRUE)+
geom_vline(xintercept=0.78, linetype="dotted", color="red")

#Hbyh #Others
Hbyhplot<-function(plotdata,Shapedata,Title)
{
  ggplot(data=plotdata)+geom_point(aes(x=Hbyh,y=R),shape=Shapedata)+
    labs(title=Title,x="Relative wave height, H/h",y="Wave Attenuation %")+theme_bw()
}

MGHbyhplot<-Hbyhplot(mgdata,0,"Mangroves")
SMHbyhplot<-Hbyhplot(smdata,5,"Salt-Marshes")
SGHbyhplot<-Hbyhplot(sgdata,8,"Seagrass/Kelp")

grid.arrange(REEFHbyhplot, MGHbyhplot,SMHbyhplot,SGHbyhplot,ncol=2)

#hv/h #MANGROVES
MVRdata_R_MG_sub<-mgdata_submerge[c("Habitat","Hs.i","Hs.t","R","hvbyh")]
MVRdata_R_MG_sub<-MVRdata_R_MG_sub[complete.cases(MVRdata_R_MG_sub),]
MVRdata_R_MG_out<-subset(MVRdata_R_MG_sub, MVRdata_R_MG_sub$hvbyh<0.5)
MVRdata_R_MG_sub<-subset(MVRdata_R_MG_sub, MVRdata_R_MG_sub$hvbyh>0.5) ## Caution! Removing outlier!
MVRdata_R_MG_glm_sub<-lm(R~hvbyh, MVRdata_R_MG_sub)
summary(MVRdata_R_MG_glm_sub) 
MGhvbyhplot_submerge<-ggplot(data=MVRdata_R_MG_sub)+geom_point(aes(x=hvbyh,y=R),shape=0)+
  labs(title="Mangroves",x="Relative vegetation height, hv/h",y="Wave Attenuation %")+
  theme_bw()+geom_smooth(data=MVRdata_R_MG_sub,aes(x=hvbyh,y=R),method=lm,formula=y~x,se=FALSE,linetype=2,colour="black")+
  annotate('text',x=0.7,y=0.36,label=paste("R^2==",format(summary(MVRdata_R_MG_glm_sub)$adj.r.squared, digits=4)),parse=TRUE)
MGhvbyhplot_submerge
MVRdata_R_MG_emg<-mgdata_emerge[c("Habitat","Hs.i","Hs.t","R","hvbyh")]
MVRdata_R_MG_emg<-MVRdata_R_MG_emg[complete.cases(MVRdata_R_MG_emg),]
MGhvbyhplot_both<-MGhvbyhplot_submerge+geom_point(data=MVRdata_R_MG_emg,aes(x=hvbyh,y=R),shape=0)+
  geom_vline(xintercept=1.0,linetype="dotted",color="red")+
  geom_point(data=MVRdata_R_MG_out,aes(x=hvbyh,y=R),shape=0) #Add back outlier point
MGhvbyhplot_both

#hv/h #MARSHES
MVRdata_R_SM_all<-smdata[c("Habitat","Hs.i","Hs.t","R","hvbyh")]
MVRdata_R_SM_all<-MVRdata_R_SM_all[complete.cases(MVRdata_R_SM_all),]
MVRdata_R_SM_sub<-smdata_submerge[c("Habitat","Hs.i","Hs.t","R","hvbyh")]
MVRdata_R_SM_sub<-MVRdata_R_SM_sub[complete.cases(MVRdata_R_SM_sub),]
MVRdata_R_SM_out<-subset(MVRdata_R_SM_sub, MVRdata_R_SM_sub$hvbyh<0.1)
MVRdata_R_SM_sub<-subset(MVRdata_R_SM_sub, MVRdata_R_SM_sub$hvbyh>0.1) ## Caution! Removing outlier!
MVRdata_R_SM_glm_sub<-lm(R~hvbyh, MVRdata_R_SM_sub)
summary(MVRdata_R_SM_glm_sub) 
SMhvbyhplot_submerge<-ggplot(data=MVRdata_R_SM_sub)+geom_point(aes(x=hvbyh,y=R),shape=8)+
  labs(title="Salt-Marshes",x="Relative vegetation height, hv/h",y="Wave Attenuation %")+
  theme_bw()+geom_smooth(data=MVRdata_R_SM_sub,aes(x=hvbyh,y=R),method=lm,formula=y~x,se=FALSE,linetype=2,colour="black")+
  annotate('text',x=0.7,y=0.36,label=paste("R^2==",format(summary(MVRdata_R_SM_glm_sub)$adj.r.squared, digits=4)),parse=TRUE)
SMhvbyhplot_submerge
MVRdata_R_SM_emg<-smdata_emerge[c("Habitat","Hs.i","Hs.t","R","hvbyh")]
MVRdata_R_SM_emg<-MVRdata_R_SM_emg[complete.cases(MVRdata_R_SM_emg),]
SMhvbyhplot_both<-SMhvbyhplot_submerge+geom_point(data=MVRdata_R_SM_emg,aes(x=hvbyh,y=R),shape=8)+
  geom_vline(xintercept=1,linetype="dotted",color="red")+
  geom_point(data=MVRdata_R_SM_out,aes(x=hvbyh,y=R),shape=8) #Add back outlier point
SMhvbyhplot_both

#hv/h SEAGRASS/KELP
MVRdata_R_SG<-sgdata[c("Habitat","Hs.i","Hs.t","R","hvbyh")]
MVRdata_R_SG<-MVRdata_R_SG[complete.cases(MVRdata_R_SG),]
MVRdata_R_SG_glm<-lm(R~hvbyh, MVRdata_R_SG)
summary(MVRdata_R_SG_glm) 
SGhvbyhplot<-ggplot(data=MVRdata_R_SG)+geom_point(aes(x=hvbyh,y=R),shape=0)+
  labs(title="Seagrass/Kelp",x="Relative vegetation height, hv/h",y="Wave Attenuation %")+
  theme_bw()+theme(legend.position=c(1,0),legend.justification=c(0,1))+
  geom_smooth(data=MVRdata_R_SG,aes(x=hvbyh,y=R),method=lm,formula=y~x,se=FALSE,linetype=2,colour="black")+
  annotate('text',x=0.7,y=0.36,label=paste("R^2==",format(summary(MVRdata_R_SG_glm)$adj.r.squared, digits=4)),parse=TRUE)
SGhvbyhplot

grid.arrange(MGhvbyhplot,SMhvbyhplot,nrow=2)

## THE END ##