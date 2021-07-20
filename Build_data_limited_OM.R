
#install.packages("openMSE") #Only run once
library(openMSE)

#-----------------------------------------
#Data limited approach to OM construction
#-----------------------------------------

#-----------------------------------------
#Step 1. Create file system for OM objects 
#-----------------------------------------

#Initialize a new OM
#Creates blank set of files in the working directory, initializes a blank object. Overwrite=FALSE prevents overwriting an existing OM
#To get started, We will create an OM that contains pre-populated generic options for: Fleet, Obs, & Imp 
#Note that we will leave the Stock object blank for now, allowing us to focus on specifying eel life history parameters

#See Fleet options
avail('Fleet')
plot(Generic_IncE, Albacore) 

#See Obs options
avail('Obs')
plot(Precise_Unbiased, Albacore) 

#See Imp options
avail('Imp')
plot(Perfect_Imp, Albacore) 

#Create file system for OM - note that Stock is not specified. Slots will be left blank, we will populate later
OMinit(name = "eel_data_limited_OM", Generic_IncE, Precise_Unbiased, Perfect_Imp, overwrite=FALSE)


#--------------------------------------------------
#Step 2. Populate parameter values for Stock object
#--------------------------------------------------

#The easiest way to do this is to enter the parameter values as ranges in the Excel spreadsheets that were created in step 1.

#As you enter parameter values in excel, enter notes in the eel_OM.rmd file. Later, we will use R markdown to auto-create documentation for OM

#Parameter choices and corresponding notes have already been made. This is good opportunity to review the preliminary choices that were made by the analyst

#Explore fishlife for steepness
#Life history for Pacific chub
library(FishLife)
fish_out<-Plot_taxa(Search_species(Genus="Ophichthus", Species="remiger", add_ancestors=TRUE)$match_taxonomy, mfrow=c(3,2))


#-------------------------------------------------
#Step 3. Load the OM from excel into R
#-------------------------------------------------

eel_data_limited_OM <- XL2OM('eel_data_limited_OM')


#--------------------------------------------------------------------
#Step 4. Use R markdown to auto-create html documentation for the OM
#---------------------------------------------------------------------

OMdoc('eel_data_limited_OM')



