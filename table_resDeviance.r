hmods <- c("out_H_Null_Mod.rData",
           "out_H_Mar_Mod.rData", 
           "out_H_Day_Mod.rData", 
           "out_H_DayYear_Mod.rData", 
           "out_H_NoMar_Mod.rData", 
           "out_H_bestMod.rData")

wmods <- c("out_W_Null_Mod.rData",
           "out_W_Mar_Mod.rData", 
           "out_W_Day_Mod.rData", 
           "out_W_DayYear_Mod.rData", 
           "out_W_NoMar_Mod.rData", 
           "out_W_bestMod.rData")

df <- data.frame(Hatchery=NA,Wild=NA,H_percent=NA,W_percent=NA)
icnt <- 1
for(i in hmods){
  load(i)
  # df[icnt,1] <- sum((data$s_k/data$s_n - rep$s_hat)^2 * data$s_n) 
  df[icnt,1] <- out$objective
  df[icnt,3] <- 1-df[icnt,1]/df[1,1]
  icnt <- icnt + 1
}
icnt <- 1
for(i in wmods){
  load(i)
  # df[icnt,2] <- sum((data$s_k/data$s_n - rep$s_hat)^2 * data$s_n) 
  df[icnt,2] <- out$objective
  df[icnt,4] <- 1-df[icnt,2]/df[1,2]
  icnt <- icnt + 1
}
# df$H_percent <- (-df$Hatchery+max(df$Hatchery))/max(df$Hatchery)
# df$W_percent <- (-df$Wild+max(df$Wild))/max(df$Wild)

row.names(df) <- c("Null", "Marine", "Day", 
                   "Day/year","Day + day/year", 
                   "Marine covariates + day + day/year")

write.csv(df,"table_resDeviance.csv", sep="\t")