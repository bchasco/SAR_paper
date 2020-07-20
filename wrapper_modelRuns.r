library(TMB)
library(TMBhelper)
# rm(list=ls())
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
rearType <- c("H")
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


saveOutput <- FALSE
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
simCor_j <- 0.9
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
          for(rr in c("H")){
            jcnt <- jcnt + 1
            rear <- rr
            tmpMarVars <- c(11,24)# W = 15,36, H = 11,24, c(myTmpMarVars[,mm]) #The temporary marine variables
            
            source("create_DataAndPars.r")
            if(numberOfMarineVariables>1){
              if(abs(cor(subData[tmpMarVars])[1,2])>0.7){
                next
              }
            }
            file<-paste0("out_",rear,"_NoMar_Mod.rData",collapse = "")
            
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
# nsim <- 500
# simMat <- matrix(NA,nsim,10)
# for(i in 1:nsim){
#   set.seed(10*i)
# 
#   sim <- obj$simulate(complete=TRUE)
#   sim$s_k <- sim$sim_k
#   sim$s_n <- sim$s_n * 1
#   if(sum(is.na(sim$sim_k))==0){
#     obj_sim <- MakeADFun(data = sim,
#                          parameters = parameters,
#                          map = myMap,
#                          random=c("eps_j"
#                                   ,"eps_t"
#                                   ,"eps_jt"
#                                   ,"eps_x"
#                                   ,"frho_Rx"
#                          ),
#                          silent = TRUE,
#                          bias.correct=TRUE,
#                          DLL = "integrated2")
#     
#     
#     sim_out <- nlminb(obj_sim$par,obj_sim$fn,obj_sim$gr)
#     
#     sim_rep <- obj_sim$report()
#     
#     simMat[i,] <- c(exp(sim_rep$mu_s),
#                     plogis(sim_rep$frho_j),
#                     atan(sim_rep$frho_t)*2/3.154,
#                     atan(sim_rep$frho1_jt)*2/3.154,
#                     atan(sim_rep$frho2_jt)*2/3.154,
#                     exp(sim_rep$fpsi_j),
#                     exp(sim_rep$fpsi_t),
#                     exp(sim_rep$fpsi_jt),
#                     sim_rep$beta_mar)
#     
#       # unlist(sim_out$par),sum(sim$s_k))
# 
#     par(mfrow=c(1,1))
#     data.frame(n=sim$s_n,k=sim$sim_k,j=sim$j)
#     df <- data.frame(n=sim$s_n,k=sim$sim_k,j=sim$j)
#     ag <- aggregate(list(k=df$k,n=df$n),by=list(df$j),sum)
#     par(mfrow=c(2,2))
#     
#     plot(ag$Group.1,ag$k/ag$n,
#          cex=1.3,
#          ylab="Survival",
#          col="darkgrey",
#          pch=16)
#     lines(plogis(sim_rep$mu_s+sim_rep$eps_j))
#     
#     plot(sim$s_eps_j)
#     lines(sim_rep$eps_j)
#     
#     matplot(sim_rep$eps_jt[1,,c(1,10)], 
#             type="l",
#             ylab="dayXyear deviates",
#             lty=1,
#             lwd=2)
#     matpoints(sim$s_eps_jt[,c(1,10)], 
#               pch=16,
#               cex=1.2)
#     
#     
#     print(paste(i,sum(sim$sim_k>0),sum(sim$sim_k)))
#   }
# }
# 
# 
# library(reshape2)
# simMat <- as.data.frame(simMat)
# names(simMat) <- c("mu_s","rho_j","rho_t","rho1_jt","rho2_jt","psi_j","psi_t","psi_jt","beta1","beta2")
# myLines <- data.frame(pars=names(simMat),
#                       value=c(exp(rep$mu_s),
#                               simCor_j,#plogis(rep$frho_j),
#                               atan(rep$frho_t)*2/3.154,
#                               atan(rep$frho1_jt)*2/3.154,
#                               atan(rep$frho2_jt)*2/3.154,
#                               exp(rep$fpsi_j),
#                               exp(rep$fpsi_t),
#                               exp(rep$fpsi_jt),
#                               rep$beta_mar))
# simMat <- as.data.frame(t((t(simMat)-myLines$value)/myLines$value))
# simMelt <- melt(simMat,
#                 measure.vars = names(simMat),
#                 variable.name="pars")
# 
# for(ii in names(simMat)){
#   if(!is.na(mean(simMelt$value[simMelt$pars==ii]))){
#     dat_q <- quantile(simMelt$value[simMelt$pars==ii], probs=c(0.025,0.975))
#     simMelt <- simMelt[!(simMelt$pars==ii & (simMelt$value<dat_q[1] | simMelt$value>=dat_q[2])),]
#   }
# }
# 
# library(ggplot2)
# 
# 
# p <- ggplot(simMelt, aes(factor(pars),value)) +
#   facet_wrap(~pars, scales = "free") +
#   geom_violin(aes(factor(pars)))+
#   geom_hline(aes(yintercept=0), colour="blue")
# 
# print(p)
