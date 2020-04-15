library(TMB)
library(TMBhelper)
library(lme4)

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
fixed_mar <- 1 #mar/env surv parameters
mean_s <- 1 #Mean survival parameters
cov_pars <- 0 #covariance parameters

calibration_flag <- 0 #run calibration (0= NO, 1= Yes)
calibration_files <- c("calibration.out",
                       "ch1_Historical_2020rerun_wild_sep90_bonarrival.out",
                       "ch1_Historical_2020rerun_wild_sep50_bonarrival.out",
                       "ch1_Historical_2020rerun_wild_bonarrival.out")
calibration_file <- calibration_files[1]


retro <- 0 # 1 = do the retrospective, 0 = don't do the retrospective, DEPRECATED
#Do you want to create SAR projections for all of the models you just ran
SensitivitySimulation <- 1 #Just pick a simulation number it doesn't matter which one, DEPRECATED
createSARProjections <- FALSE #Do you want to also do the projections, DEPRECATED

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

#Recompile if necessary
compile("integrated3_sim.cpp")
dyn.load("integrated3_sim.dll")

simCor_j <- 0.5
simCor_jt1 <- 0.5
simCor_jt2 <- 0.5


nk_dim <- 1
getSD <- FALSE


#Run the best fit model
re_j <- 1
re_t <- 0
re_jt <- 1
rear <- "W"
tmpMarVars <- c(15,36) #lowest AIC model

best_eps_j <- rep(0,91)
use_best_eps_j <- 0 #Don't use best eps_j to start with

#Create the objects and estimate the parameters for the lowest AIC model
source("create_DataAndPars_sim3.r")
source("create_MapAndObj_sim3.r")
out <- nlminb(obj$par,obj$fn,obj$gr)

#REport the results
rep_true <- obj$report()
#Save the original samples size
original_s_k <- data$s_k
original_s_n <- data$s_n


#Form of the lme4 model
form1<-as.formula(paste("resMatrix ~ var1 + var2 + julian + julian2 + day*yr + (1|yr)"))

betas <- matrix(0,length(data$yr),2)
for(i in 1:length(data$yr)){
  betas[i,1] <- subData[data$yr[i]+1, tmpMarVars[1]]
  betas[i,2] <- subData[data$yr[i]+1, tmpMarVars[2]]
}


best_eps_j <- rep_true$eps_j[,1]

sim_n <- 300
mu_lme <- as.data.frame(matrix(NA,sim_n,3))
mu_TMB <- as.data.frame(matrix(NA,sim_n,3))
sd_lme <- as.data.frame(matrix(NA,sim_n,3))
sd_TMB <- as.data.frame(matrix(NA,sim_n,3))

names(mu_lme) <- c("mu", "beta[SST]","beta[CUI]")
names(mu_TMB) <- c("mu", "beta[SST]","beta[CUI]")
names(sd_lme) <- c("mu", "beta[SST]","beta[CUI]")
names(sd_TMB) <- c("mu", "beta[SST]","beta[CUI]")

# set.seed(800)
for(ii in 1:sim_n){
  
  #Flags
  use_best_eps_j <- 0 #use the eps_j from the MLE model (1=yes, 0=randomly generate)
  re_j <- 1 #Estimate daily random effects
  
  source("create_DataAndPars_sim3.r")
  source("create_MapAndObj_sim3.r")

  #Simulated data set
  sim <- obj$simulate(complete=TRUE)
  
  #Simulated number of fish surviving
  resMatrix <- as.matrix(cbind(sim$sim_k,sim$s_n))
  
  #Create the dataframe for the lme4
  dfSim <- data.frame(yr=sim$yr,
                   julian = sim$j,
                   julian2 = (sim$j-mean(sim$j))^2,
                   var1=betas[,1],
                   var2=betas[,2])
  #Standardize the julian date to make the estimation easier
  dfSim$julian <- (dfSim$julian - mean(dfSim$julian))/sd((dfSim$julian))
  dfSim$julian2 <- (dfSim$julian2 - mean(dfSim$julian2))/sd((dfSim$julian2))
  #Add the resMatrix to the data.frame
  dfSim$resMatrix <- resMatrix
  
  #update the data and the parameters
  data$s_k <- sim$sim_k
  parameters$eps_j <- parameters$eps_j*0
  # #Create the object
  sim_obj <- MakeADFun(data = data,
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
                       DLL = "integrated3_sim")
  out_sim <- nlminb(sim_obj$par,sim_obj$fn,sim_obj$gr)   #Run the model
  rep_sim <- sim_obj$report()   #Report the results

  dfMelt <- matrix(0,nrow = sum(dfSim$resMatrix[,2]), ncol = ncol(dfSim))
  df_1 <- as.data.frame(cbind(rep(dfSim$yr,dfSim$resMatrix[,1]),
        rep(dfSim$julian,dfSim$resMatrix[,1]),
        rep(dfSim$julian2,dfSim$resMatrix[,1]),
        rep(dfSim$var1,dfSim$resMatrix[,1]),
        rep(dfSim$var2,dfSim$resMatrix[,1])))
  df_1$surv <- 1
  df_0 <- as.data.frame(cbind(rep(dfSim$yr,(dfSim$resMatrix[,2]-dfSim$resMatrix[,1])),
                rep(dfSim$julian,(dfSim$resMatrix[,2]-dfSim$resMatrix[,1])),
                rep(dfSim$julian2,(dfSim$resMatrix[,2]-dfSim$resMatrix[,1])),
                rep(dfSim$var1,(dfSim$resMatrix[,2]-dfSim$resMatrix[,1])),
                rep(dfSim$var2,(dfSim$resMatrix[,2])-dfSim$resMatrix[,1])))
  df_0$surv <- 0
  names(df_0) <- c(names(dfSim)[1:5],"surv")
  names(df_1) <- c(names(dfSim)[1:5],"surv")
  dfMelt <- as.data.frame(rbind(df_0,df_1))
  mod1 <- glm(surv~var1+var2+julian+julian2+yr*julian, 
                # control=glmerControl(optimizer = "bobyqa"),
                # REML=FALSE,
                data=dfMelt, 
               family="binomial")
  
  #Get the standard for the TMB model
  SD <- sdreport(sim_obj)

  print(max(sim_obj$gr()))
  print(SD)
  
  #Save the parameter and standard error estimates
  mu_lme[ii,] <- c(plogis(coef(mod1)[1]),coef(mod1)[2:3])
  mu_TMB[ii,] <- c(plogis(SD$value[1]),SD$value[2:3])
  sd_lme[ii,] <- sqrt(diag(vcov(mod1)))[1:3]
  sd_TMB[ii,] <- SD$sd

  mod1_pred <- plogis(predict(mod1))
  
  # #Mized effect model using LME4 
  # df_obs <- data.frame(out=rep("obs",length(data$s_k)),yr=data$yr+2000, j=data$j+minJ, val=data$s_k/data$s_n)
  # df_MV <- data.frame(out=rep("MV",length(data$s_k)),yr=data$yr+2000, j=data$j+minJ, val=(rep_sim$s_hat))
  # # df_re <- data.frame(out=rep("re",length(data$s_k)),yr=data$yr+2000, j=data$j+minJ, val=mod1_pred)
  # df_fixed <- data.frame(out=rep("fixed",nrow(dfMelt)),yr=dfMelt$yr+2000, j=rep(data$j,times=dfSim$resMatrix[,2])+minJ, val=mod1_pred)
  # 
  # 
  # df_MV_mu <- data.frame(out=rep("MV_mean",(max(data$j)+1)*(max(data$yr)+1)),
  #                        yr=rep(sort(unique(data$yr)+2000),each=length(data$best_eps_j)),
  #                        j=rep(minJ:maxJ,times=length(unique(data$yr))),
  #                        val=plogis(rep_sim$mu_s+c(t(t(matrix(sim$s_eps_j,length(sim$s_eps_j),length(rep_sim$eMar2)))+rep_sim$eMar2))))
  # dfAll <- rbind(df_obs,
  #             df_MV,
  #             df_fixed,
  #             df_MV_mu)
  print(ii)
  # source("fig_quadraticRealization_ggplot.r")
}

mu_lme$Model <- "glm"
mu_TMB$Model <- "TMB"
mu_lme$stat <- "A"
mu_TMB$stat <- "A"

sd_lme$Model <- "glm"
sd_TMB$Model <- "TMB"
sd_lme$stat <- "B"
sd_TMB$stat <- "B"

df <- rbind(mu_lme,mu_TMB,sd_lme,sd_TMB)

df <- melt(df,
                # measure.vars = names(simMat),
                id.vars = c("Model","stat"),
                variable.name="pars")

df <- na.omit(df)
df <- df[!(df$pars=="beta[SST]" & df$stat=="B" & df$value>100),]
df_par <- df
df$value[df$pars=="mu" & df$stat=="A"] <- (df$value[df$pars=="mu" & df$stat=="A"] - plogis(rep_true$mu_s))/plogis(rep$mu_s)*100
df$value[df$pars=="beta[CUI]" & df$stat=="A"] <- (df$value[df$pars=="beta[CUI]" & df$stat=="A"] - (rep_true$beta_mar[2,1]))/(rep_true$beta_mar[2,1])*100
df$value[df$pars=="beta[SST]" & df$stat=="A"] <- (df$value[df$pars=="beta[SST]" & df$stat=="A"] - (rep_true$beta_mar[1,1]))/(rep_true$beta_mar[1,1])*100


save(file="simQuadractic.rData",df, df_par)
