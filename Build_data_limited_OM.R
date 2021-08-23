
#-----------------------------------------
#Data limited approach to OM construction
#-----------------------------------------

#devtools::install_github("Blue-Matter/openMSE")
#install.packages("openMSE") #Only run once

library(openMSE)


#-----------------------------------------
#Step 1. Create file system for OM objects 
#-----------------------------------------

#Initialize a new OM
#Creates blank set of files in the working directory, initializes a blank object. Overwrite=FALSE prevents overwriting an existing OM
#To get started, We will create an OM that contains pre-populated generic options for: Fleet, Obs, & Imp 
#Note that we will leave the Stock object blank for now, allowing us to focus on specifying eel life history parameters

#Create file system for OM - note that Stock is not specified. Slots will be left blank, we will populate later
OMinit(name = "eel_data_limited_OM", Generic_IncE, Precise_Unbiased, Perfect_Imp, overwrite=FALSE)


#-------------------------------------------------------------------
#Step 2. Populate parameter values for Stock object & Fleet object
#-------------------------------------------------------------------

#The easiest way to do this is to enter the parameter values as ranges in the Excel spreadsheets that were created in step 1.

#As you enter parameter values in excel, enter notes in the eel_OM.rmd file. Later, we will use R markdown to auto-create documentation for OM

#Parameter choices and corresponding notes have already been made. This is a good opportunity to review the preliminary choices that were made by the analyst

#Explore fishlife for steepness
library(FishLife)
fish_out<-Plot_taxa(Search_species(Genus="Ophichthus", Species="remiger", add_ancestors=TRUE)$match_taxonomy, mfrow=c(3,2))


#-------------------------------------------------
#Step 3. Load the OM from excel into R
#-------------------------------------------------

eel_data_limited_OM <- XL2OM('eel_data_limited_OM')


#--------------------------------------------------
#Step 4. Tuning to depletion of VULNERABLE biomass
#--------------------------------------------------

#Objective: Tune Depletion (in units spawning biomass), such that 95% of simulations have vulnerable biomass between 0.6 and 0.76
eel_data_limited_OM@D<-c(0.6, 0.7)
eel_data_limited_OM@nsim<-200
Hist<-Simulate(eel_data_limited_OM)
X<-rowSums(Hist@TSdata$VBiomass[, eel_data_limited_OM@nyears,]) /  Hist@Ref$ReferencePoints$VB0
Y<-rowSums(Hist@TSdata$SBiomass[, eel_data_limited_OM@nyears,]) /  Hist@Ref$ReferencePoints$SSB0
cbind(X,Y)
hist(X)
quantile(X, probs=c(0.025, 0.975))

#Updated saved OM
OMinit(name = "eel_data_limited_OM", eel_data_limited_OM, overwrite=TRUE)

#Generate plots of the historical dynamics
plot(Hist, html=FALSE)

#Generate plot for trap time series used in tuning
png(file=paste0(getwd(), "/images/Relative effort.png"), width=6, height=4, units="in", res=200, bg="white",pointsize=12)  
plot(eel_data_limited_OM@EffYears, eel_data_limited_OM@EffUpper, type = "l", las=1, xlab = "Year", ylab = "Relative fishing effort")
dev.off()

#--------------------------------------------------------------------
#Step 4. Use R markdown to auto-create html documentation for the OM
#---------------------------------------------------------------------

OMdoc('eel_data_limited_OM')


#-----------------------------------------------------------------
#Step 5. Specify a preliminary set of management procedures (MPs)
#-----------------------------------------------------------------

eel_data_limited_OM@interval<-5
  
################
#Note that some our MPs rely on specified quantities used assessment, will set these. 
#openMSE has default assumptions for specifying these values, if they are not set by the user
#These quantities are guesses and need to be reviewed
eel_data_limited_OM@reps<-100
Data <- new('Data')
Data@Mort<-median(eel_data_limited_OM@M)
Data@CV_Mort<-0.1  
Data@steep<-median(eel_data_limited_OM@h)
Data@CV_steep<-0.3  #These quantities are guesses and need to be reviewed
Data@vbK<-median(eel_data_limited_OM@K)
Data@CV_vbK<-0.1
Data@vbLinf<-median(eel_data_limited_OM@Linf)
Data@CV_vbLinf<-0.1
Data@vbt0<-median(eel_data_limited_OM@t0)
Data@MaxAge<-eel_data_limited_OM@maxage
Data@L50<-median(eel_data_limited_OM@L50)
Data@CV_L50<-0.1
Data@wla<-eel_data_limited_OM@a
Data@wlb<-eel_data_limited_OM@b
eel_data_limited_OM@cpars$Data <- Data

#################
#Create plots of priors 
png(file=paste0(getwd(), "/images/R_prior.png"), width=6, height=6, units="in", res=200, bg="white",pointsize=12)  
reps<-10000
Mvec <- MSEtool::trlnorm(reps, Data@Mort, Data@CV_Mort)
Kvec <- MSEtool::trlnorm(reps, Data@vbK, Data@CV_vbK)
Linfvec = MSEtool::trlnorm(reps, Data@vbLinf, Data@CV_vbLinf)
t0vec <- rep(Data@vbt0, reps)
t0vec[!is.finite(t0vec)] <- 0
hvec <- MSEtool::sample_steepness2(reps, Data@steep, Data@CV_steep)
rsamp <- getr(x=1, Data, Mvec, Kvec, Linfvec, t0vec, hvec, 
              maxage = Data@MaxAge, r_reps = reps)
par(mfrow=c(2,2))
hist(rsamp, las=1,prob=TRUE, xlab="r", main="r prior from M & h")
lines(density(rsamp), col="red", lwd=2)
rprior<-rnorm(reps, 0.9, 0.1)
hist(rprior, las=1, prob=TRUE, xlab="r", main="r prior IMARPE")
lines(density(rprior), col="red", lwd=2)
hist(Mvec, prob=TRUE, las=1, xlab="M", main="M prior")
lines(density(Mvec), col="red", lwd=2)
hist(hvec, las=1, prob=TRUE, xlab="h", main="h prior")
lines(density(hvec), col="red", lwd=2)
dev.off()

###############
#Create the data-rich MPs for surplus production models
SP_SS1 <- make_MP(SP_SS, HCR_MSY, fix_dep=TRUE, fix_n=TRUE, use_r_prior=TRUE,  SR_type = "BH")
SP_SS2 <- make_MP(SP_SS, HCR_MSY, fix_dep=TRUE, fix_n=TRUE, use_r_prior=TRUE, start = list(r_prior = c(0.9, 0.10)))


#################
#Combine MPs into a vector, adding reference MPs and data-limited MPs for comparison
MPs<-c("SP_SS1", "SP_SS2", "Lratio_BHI", "Islope1", "LstepCC3", "FMSYref")

#-------------------------------------------------------
#Step 6.Run MSE
#------------------------------------------------------

#Set sims, set projection time period
eel_data_limited_OM@nsim<-48
eel_data_limited_OM@proyears<-100

#Run MSE
MSE <- runMSE(OM=eel_data_limited_OM, MPs=MPs, checkMPs = FALSE, parallel = TRUE)

#Save the MSE
saveRDS(MSE, file="MSE.rds")

#Load the MSE (instead of re-running)
MSE<-readRDS("MSE.rds")



#---------------------------------------------------------------------------
#Step 7. Evaluate the results, including application of performance metrics
#---------------------------------------------------------------------------

#Check available performance metrics
avail("PM")

#----------------------------------------
# PNOF Probability of not overfishing 
# Probabilyt F < FMSY across years
# This is a DLMtool built in function

PNOF(MSE) #Across all years
PNOF(MSE, Yrs=-10) #Across last 10 years
PNOF(MSE, Yrs=10) #Across first 10 years

#----------------------------------------
# P50 Probability that spawning biomass > 0.5SBmsy
# Over all projection years
# This is a DLMtool built in function

P50(MSE) #Across all years
P50(MSE, Yrs=-10)  #Across last 10 years
P50(MSE, Yrs=10) #Across first 10 years

#-------------------------------------------------------------------------------
#LTY   Long term yield
#Probability Yield > 0.5 FMSY Yield in the last 10 years of the projection
#This is a  DLMtool built in function 

LTY(MSE)

#-------------------------------------------------------------------------------
#STY   Short term yield
#Probability Yield > 0.5 FMSY Yield in the first 10 years of the projection
#This is a  DLMtool built in function 

STY(MSE)

#-------------------------------------------------------------------------------
#Interannual variation in yield
#Probability that inter-annual change in yield < 20%, across years
#This is a  DLMtool built in function 

AAVY(MSE) #Across all years
AAVY(MSE, Yrs=-10) #Across last 10 years
AAVY(MSE, Yrs=10) #Across first 10 years

#------------------------------------------------------------------------------
#Mertic Relatuve change in catches
#Performance: Mean proportion of years that metric exceeded 1 (i.e., Ref = 1)
#------------------------------------------------------------------------------

CREL <- function(MSEobj = NULL, Ref = 0.8, Yrs = NULL) {
  Yrs <- ChkYrs(Yrs, MSEobj)
  PMobj <- new("PMobj")
  PMobj@Name <- paste0("Prob Yield > 0.8 current (Years ",  Yrs[1], " - ", Yrs[2], ")")
  PMobj@Caption<- PMobj@Name
  PMobj@Ref <- Ref
  
  #Calculate a metric for each sim, MP, and year
  RefC<-apply(MSEobj@CB_hist[,(MSEobj@nyears-4):MSEobj@nyears], MARGIN=1, FUN=mean, na.rm=TRUE)
  PMobj@Stat <- array(NA, dim=c(MSEobj@nsim, MSEobj@nMPs, MSEobj@proyears))
  for(i in 1:MSEobj@nsim){
    PMobj@Stat[i,,]<-MSEobj@Catch[i,,] / RefC[i]
  }
  PMobj@Stat <- PMobj@Stat[,,Yrs[1]:Yrs[2]]
  PMobj@Prob <- calcProb(PMobj@Stat > PMobj@Ref, MSEobj)
  PMobj@Mean <- calcMean(PMobj@Prob)
  PMobj@MPs <- MSEobj@MPs
  PMobj
}

class(CREL) <- "PM"

CREL(MSE) #All years
CREL(MSE, Yrs=-10) #Last 10 years
CREL(MSE, Yrs=10) #First 10 years

#----------------------------------------------------------------------------------
#Metric: Relative change in spawning biomass
#Performance: Mean proportion of years that metric exceeded 1 (i.e., Ref = 1)

BREL <- function(MSEobj = NULL, Ref = 0.8, Yrs = NULL) {
  Yrs <- ChkYrs(Yrs, MSEobj)
  PMobj <- new("PMobj")
  PMobj@Name <- paste0("Prob SSB > 0.8 current (Years ",  Yrs[1], " - ", Yrs[2], ")")
  PMobj@Caption<- PMobj@Name
  PMobj@Ref <- Ref
  
  #Calculate a metric for each sim, MP, and year
  RefB<-apply(MSEobj@SSB_hist[,(MSEobj@nyears-4):MSEobj@nyears], MARGIN=1, FUN=mean, na.rm=TRUE)
  PMobj@Stat <- array(NA, dim=c(MSEobj@nsim, MSEobj@nMPs, MSEobj@proyears))
  for(i in 1:MSEobj@nsim){
    PMobj@Stat[i,,]<-MSEobj@SSB[i,,] / RefB[i]
  }
  PMobj@Stat <- PMobj@Stat[,,Yrs[1]:Yrs[2]]
  PMobj@Prob <- calcProb(PMobj@Stat > PMobj@Ref, MSEobj)
  PMobj@Mean <- calcMean(PMobj@Prob)
  PMobj@MPs <- MSEobj@MPs
  PMobj
}

class(BREL) <- "PM"

BREL(MSE) #All years
BREL(MSE, Yrs=-10) #Last 10 years
BREL(MSE, Yrs=10) #First 10 years

#----------------------
#Plots
#----------------------

#---------------
#Time series
png(file=paste0(getwd(), "/images/Time_series_MSYseeking.png"), width=8, height=6, units="in", res=200, bg="white",pointsize=12)  
Pplot2(MSE, MPs=MSE@MPs[1:3], traj="quant",  YVar = c("F_FMSY", "SSB_SSBMSY", "Yield"), quants=c(0.25, 0.75),  oneIt = FALSE)
dev.off()

png(file=paste0(getwd(), "/images/Time_series_stableBio.png"), width=8, height=6, units="in", res=200, bg="white",pointsize=12)  
Pplot2(MSE, MPs=MSE@MPs[4:6], traj="quant",  YVar = c("F_FMSY", "SSB_SSBMSY", "Yield"), quants=c(0.25, 0.75),  oneIt = FALSE)
dev.off()
#------------------------------------------
#Stock status - last 5 years of projections
#png(file=paste0(getwd(), "/images/Kobe.png"), width=7, height=5, units="in", res=200, bg="white",pointsize=12)  
#Kplot(MSE)
#dev.off()

#png(file=paste0(getwd(), "/images/Cplot.png"), width=6, height=6, units="in", res=200, bg="white",pointsize=12)  
#Cplot(MSE, lastYrs = 10)
#dev.off()

#----------------
#Trade-off plots

#Long term
png(file=paste0(getwd(), "/images/Long_term_prob.png"), width=8, height=8, units="in", res=200, bg="white",pointsize=12)  
TradePlot(MSE, PMlist=list("P50", "PNOF", "LTY", "AAVY", "BREL", "CREL"), Yrs=list(P50=-10, PNOF=-10, AAVY=-10, CREL=-10, BREL=-10), Satisficed=FALSE, legend=FALSE, Lims = 0)
dev.off()

#Short term
png(file=paste0(getwd(), "/images/Short_term_prob.png"), width=8, height=8, units="in", res=200, bg="white",pointsize=12)  
TradePlot(MSE, PMlist=list("P50", "PNOF", "STY", "AAVY", "BREL", "CREL"), Yrs=list(P50=10, PNOF=10, AAVY=10, CREL=10, BREL=10), Satisficed=FALSE, legend=FALSE,  Lims = 0)
dev.off()

#All years
png(file=paste0(getwd(), "/images/All_years_prob.png"), width=8, height=8, units="in", res=200, bg="white",pointsize=12)  
TradePlot(MSE, PMlist=list("P50", "PNOF", "STY", "AAVY", "BREL", "CREL"), Satisficed=FALSE, legend=FALSE, Lims = 0)
dev.off()
