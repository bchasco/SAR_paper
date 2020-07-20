library(pROC)
pred <- rep(rep$s_hat,data$s_n)
obs <- NA
for(i in 1:length(data$s_n)){
  obs <- c(obs,rep("surv",data$s_k[i]),rep("dead",data$s_n[i]-data$s_k[i]))
}
obs <- na.omit(obs)
roc(obs,pred)
sum(data$s_n)
  