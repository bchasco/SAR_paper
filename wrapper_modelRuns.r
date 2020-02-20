library(TMB)
library(TMBhelper)
# rm(list=ls())
rootDir <- "C:/Chasco/PROJECTS/SAR_PAPER/" #root directory
setwd(rootDir)
eStartYr <- 2000 #start year of the environmental data
eLastYr <- 2015 #end year of the environmental data
startYr <- 2000 #start year of the pit data
lastYr <- 2015 #end year of the pit data
minJ <- 100
maxJ <- 190
ny <- 75 #number of projection years
setSeed <- TRUE
runSensitivity <- FALSE #Do you want to run a sensitivity analysis
numberOfMarineVariables <- 1 #number of marine variables
#Estimated (1 is estimate, 0 is don't estimate)
rearType <- c("W")
re_j <- 1 #Day effect
re_t <- 0 #Year effect
re_jt <- 1 #Day X Year effect
fixed_mar <- 1 #mar/env surv parameters
mean_s <- 1 #Mean survival parameters
cov_pars <- 0 #covariance parameters

calibration_flag <- 0 #run calibration (0= NO, 1= Yes)
calibration_files <- c("calibration.out",
                       "ch1_Historical_2020rerun_wild_sep90_bonarrival.out",
                       "ch1_Historical_2020rerun_wild_sep50_bonarrival.out",
                       "ch1_Historical_2020rerun_wild_bonarrival.out")
calibration_file <- calibration_files[1]

reCompile <- TRUE
getSD <- TRUE

retro <- 0 # 1 = do the retrospective, 0 = don't do the retrospective
#Do you want to create SAR projections for all of the models you just ran
nsim <- 1 #number of simulations
SensitivitySimulation <- 1 #Just pick a simulation number it doesn't matter which one
createSARProjections <- FALSE #Do you want to also do the projections
simRCP <- c("stationary", "8.5")
simQuants <- c("50")
if(setSeed){
  set.seed(100)
}

AIC <- 100000

#marine variables

marVars <- c(
  'ersstArc.win'
  ,'ersstWAcoast.sum'
  ,'ersstArc.spr'
  ,'cui.spr'
)
load("envData.rData")
marVars <- c(names(envdata)[2:13],names(envdata)[18:42])

# #Freshwater variables
fwVars <- c(
  'F.S2.ParrSu'
  # ,'F.S3.ParrF'
  # ,'T.S2.ParrSu','T.S4.ParrW'
  # ,'aprmayjunetempLGR','aprmayjuneflowLGR'
  # ,'aprmayjunBON.scrollTEMP','aprmayjunFLOW.obs'
)
#Combine the marine and freshwater variables
myVars <-c(marVars)
#'Tranport' or 'In-river' - subset the pit data
myTrans <- c("In-river")


saveOutput <- TRUE
# AICoutput <- as.data.frame(matrix(NA,100000,10))
# names(AICoutput) <- c("AIC",
                      # "vars",
                      # "rear",
                      # "varNames",
                      # "re_j",
                      # "re_t",
                      # "re_jt",
                      # "nvar",
                      # "gr",
                      # "cor")
if(reCompile){
  try(dyn.unload("integrated2.dll"))
  compile("integrated2.cpp")
}
dyn.load("integrated2.dll")

icnt <- 0
jcnt <- 0
for(nn in 2:2){
  numberOfMarineVariables <- max(nn,1)
  if(nn==0){
    fixed_mar <- 0
  }
  if(nn>0){
    fixed_mar <- 1
  }
  #Get all of the model combinations you are interested in
  myTmpMarVars <- combn(1:length(marVars),numberOfMarineVariables)
  #Total number of models
  nmod <- ncol(myTmpMarVars)

  if(nn==0) mods <- 1:1
  if(nn>0) mods <- 1:nmod
  
  for(jj in 1:1){
    for(tt in 0:0){
      for(jt in 1:1){
        re_j <- jj
        re_t <- tt
        re_jt <- jt
        
        for(mm in 1:1){
          for(rr in c("H")){
            jcnt <- jcnt + 1
            rear <- rr
            tmpMarVars <- c(11,24)#c(myTmpMarVars[,mm]) #The temporary marine variables
            
            source("create_DataAndPars.r")
            if(numberOfMarineVariables>1){
              if(abs(cor(subData[tmpMarVars])[1,2])>0.7){
                next
              }
            }
            
            icnt <- icnt + 1
            print(paste(icnt,jcnt))
            source("create_MapAndObj.r")
            # AICoutput[icnt,1] <- thisAIC
            # if(nmod>0)
            #   AICoutput[icnt,2] <- paste(tmpMarVars,collapse = "_")
            # if(nmod==0)
            #   AICoutput[icnt,2] <- NA
            # 
            # AICoutput[icnt,3] <- rear
            # 
            # if(nn==0)
            #   AICoutput[icnt,4] <- NA
            # else
            #   AICoutput[icnt,4] <- paste(names(subData)[tmpMarVars],collapse = "_")
            # 
            # AICoutput[icnt,5] <- re_j
            # AICoutput[icnt,6] <- re_t
            # AICoutput[icnt,7] <- re_jt
            # AICoutput[icnt,8] <- nn
            # AICoutput[icnt,9] <- gr
            # if(numberOfMarineVariables>1){
            #   AICoutput[icnt,10] <- abs(cor(subData[tmpMarVars])[1,2])
            # }

            file<-paste0("out_",rear,"_bestDailyMod.rData",collapse = "")
            if(saveOutput){
              save(out,
                   obj,
                   gr,
                   thisAIC,
                   SD,
                   data,
                   rep,
                   minJ,
                   maxJ,
                   subData, #This is the pit data
                   startYr, #Initial year
                   lastYr, #Last year
                   marVars, #marine variables
                   fwVars, #Freshwater variables
                   myVars, #Marine and freshwater variables
                   tmpMarVars, #Subset of marine variables for a particular model
                   pitTotal, #Raw pit total data
                   env_mu,
                   env_sc,
                   env_mu_2000_2015,
                   env_sc_2000_2015,
                   tmpenv, #All of the environmental data
                   rearType, #The names of the transportion and rearing combinations in the model
                   file=
                     file)  
            }
            
            if(thisAIC<AIC){
              AIC <- thisAIC
            }#enf if
          }#end rear
        }#end env. model
      }#end jt
    }#end tt
  }#end jj
}#end nvars

# save(AICoutput, file="AICoutput.rData")
