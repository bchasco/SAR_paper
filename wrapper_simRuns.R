library(TMB)
library(TMBhelper)
rm(list=ls())
rootDir <- "C:/NOAA/PROJECTS/SAR_PAPER/" #root directory
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
re_j <- 0 #Day effect
re_t <- 0 #Year effect
re_jt <- 0 #Day X Year effect
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
getSD <- FALSE

retro <- 0 # 1 = do the retrospective, 0 = don't do the retrospective
#Do you want to create SAR projections for all of the models you just ran
nsim <- 1 #number of simulations
SensitivitySimulation <- 1 #Just pick a simulation number it doesn't matter which one
createSARProjections <- FALSE #Do you want to also do the projections

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
)
#Combine the marine and freshwater variables
myVars <-c(marVars)
#'Tranport' or 'In-river' - subset the pit data
myTrans <- c("In-river")

if(reCompile){
  try(dyn.unload("integrated2_sim.dll"))
  compile("integrated2_sim.cpp")
}
dyn.load("integrated2_sim.dll")
# 
icnt <- 0
jcnt <- 0

simCor_j <- 0.5
simCor_jt1 <- 0.5
simCor_jt2 <- 0.5
sim_n <- 1

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
        
        for(mm in 1:1){ #mods
          for(rr in c("W")){
            jcnt <- jcnt + 1
            rear <- rr
            tmpMarVars <- c(15,36)# W = 15,36, H = 11,24, c(myTmpMarVars[,mm]) #The temporary marine variables
            
            source("create_DataAndPars_sim.r")
            if(numberOfMarineVariables>1){
              if(abs(cor(subData[tmpMarVars])[1,2])>0.7){
                next
              }
            }
            file<-paste0("out_",rear,"_NoMar_Mod.rData",collapse = "")
            
            icnt <- icnt + 1
            print(paste(icnt,jcnt))
            source("create_MapAndObj_sim.r")
          }#end rear
        }#end env. model
      }#end jt
    }#end tt
  }#end jj
}#end nvars

# save(AICoutput, file="AICoutput.rData")
nsim <- 500
nex <- 3
nset <- 3

experiments <- matrix(0,3,3)
experiments[1,] <- c(0.1,0.5,0.9) #rho_j
experiments[2,] <- c(0.1,0.5,0.9) #rho_jt
experiments[3,] <- c(0.5,1,5) #sample size

simMat <- matrix(NA,nex*nset*nsim,10)

icnt <- 1

for(ex in 3:3){
  data$simCor_j <- 1/(1+exp(-rep$frho_j))
  data$simCor_jt1 <- atan(rep$frho1_jt)*2/3.154
  data$simCor_jt2 <- atan(rep$frho2_jt)*2/3.154
  sim_n <- 1
  for(ss in 2:2){
    source("create_DataAndPars_sim.r")
    
    if(ex==1){ #correlation for j
      data$simCor_j <- experiments[1,ss]
    }
    if(ex==2){ #correlation for jt
      data$simCor_jt1 <- experiments[2,ss]
      data$simCor_jt2 <- experiments[2,ss]
    }
    if(ex==3){ #sample size n and k
      data$s_n <- round(data$s_n * experiments[3,ss])
      data$s_k <- round(data$s_k * experiments[3,ss])
    }
    
    source("create_MapAndObj_sim.r")
    
    for(i in 1:nsim){
      set.seed(10*i)
      
      sim <- obj$simulate(complete=TRUE)
      sim$s_k <- sim$sim_k
      
      print(sum(sim$s_k))
      print(paste(ex,experiments[ex,ss],sum(sim$s_n)))
      
      if(sum(is.na(sim$sim_k))==0){
        obj_sim <- MakeADFun(data = sim,
                             parameters = parameters,
                             map = myMap,
                             random=c("eps_j"
                                      ,"eps_t"
                                      ,"eps_jt"
                                      ,"eps_x"
                                      ,"frho_Rx"
                             ),
                             silent = TRUE,
                             bias.correct=TRUE,
                             DLL = "integrated2_sim")
        
        
        sim_out <- nlminb(obj_sim$par,obj_sim$fn,obj_sim$gr)
        sim_rep <- obj_sim$report()
        print(max(obj_sim$gr()))
        SD <- sdreport(obj_sim)
        print(SD)
        #The true value depends on the experiment
        if(ex==1){
          cor_j <- sim$simCor_j  
          cor1_jt <- atan(rep$frho1_jt)*2/3.154  
          cor2_jt <- atan(rep$frho2_jt)*2/3.154  
        }
        if(ex==2){
          cor_j <- 1/(1+exp(-rep$frho_j))  
          cor1_jt <- sim$simCor_jt1  
          cor2_jt <- sim$simCor_jt2  
        }
        if(ex==3){
          cor_j <- 1/(1+exp(-rep$frho_j))  
          cor1_jt <- atan(rep$frho1_jt)*2/3.154  
          cor2_jt <- atan(rep$frho2_jt)*2/3.154  
        }
        
        simMat[icnt,] <- c(ex,
                              experiments[ex,ss],
                              (exp(sim_rep$mu_s)-exp(rep$mu_s))/exp(rep$mu_s)*100,
                        (plogis(sim_rep$frho_j) - data$simCor_j)/data$simCor_j*100,
                        # (atan(sim_rep$frho_t)*2/3.154-atan(rep$frho_t)*2/3.154)/(atan(rep$frho_t)*2/3.154)*100,
                        (atan(sim_rep$frho1_jt)*2/3.154 - data$simCor_jt1)/data$simCor_jt1*100,
                        (atan(sim_rep$frho2_jt)*2/3.154- data$simCor_jt2)/data$simCor_jt2*100,
                        (exp(sim_rep$fpsi_j) - exp(rep$fpsi_j))/exp(rep$fpsi_j)*100,
                        # (exp(sim_rep$fpsi_t)-exp(rep$fpsi_t))/exp(rep$fpsi_t)*100,
                        (exp(sim_rep$fpsi_jt)-exp(rep$fpsi_jt))/exp(rep$fpsi_jt)*100,
                        c((exp(sim_rep$beta_mar)-exp(rep$beta_mar))/exp(rep$beta_mar))*100)
        icnt <- icnt + 1
        
      }
    }
  }
}

#Melt the output
library(reshape2)
simMat <- as.data.frame(simMat)
names(simMat) <- c("ex","ss","mu_s","rho_j","rho1_jt","rho2_jt","psi_j","psi_jt","beta1","beta2")
print(simMat)
simMelt <- melt(simMat,
                # measure.vars = names(simMat),
                id.vars = c("ex","ss"),
                variable.name="pars")

#Save the output
save(file="CorrelationAndSampleSizeExperiment.rData",simMelt)

