library(TMBhelper)
Hf <- "out_H_bestMod.rData"
Hf2 <- "out_H_bestDailyMod.rData"
Wf <- "out_W_bestMod.rData"
files <- c(Hf,Wf,Hf2)
load(files[1])
par <- c(exp(rep$mu_s)
         ,NA
         ,atan(rep$frho2_jt)*2/3.154
         ,atan(rep$frho1_jt)*2/3.154
         ,NA
         ,exp(rep$fpsi_jt)
         ,rep$beta_mar)
parSD_pos <- c(exp(rep$mu_s + 1.96 * sqrt(diag(SD$cov.fixed)["mu_s"]))
                  ,NA
                  ,atan(rep$frho2_jt + 1.96 * sqrt(diag(SD$cov.fixed)["frho2_jt"]))*2/3.154
               ,atan(rep$frho1_jt + 1.96 * sqrt(diag(SD$cov.fixed)["frho1_jt"]))*2/3.154
               ,NA
                  ,exp(rep$fpsi_jt + 1.96 * sqrt(diag(SD$cov.fixed)["fpsi_jt"]))
                  ,c(rep$beta_mar + 1.96 * sqrt(diag(SD$cov.fixed)[names(diag(SD$cov.fixed))=="beta_mar"])))
parSD_neg <- c(exp(rep$mu_s - 1.96 * sqrt(diag(SD$cov.fixed)["mu_s"]))
               ,NA
               ,atan(rep$frho2_jt - 1.96 * sqrt(diag(SD$cov.fixed)["frho2_jt"]))*2/3.154
               ,atan(rep$frho1_jt - 1.96 * sqrt(diag(SD$cov.fixed)["frho1_jt"]))*2/3.154
               ,NA
               ,exp(rep$fpsi_jt - 1.96 * sqrt(diag(SD$cov.fixed)["fpsi_jt"]))
               ,c(rep$beta_mar - 1.96 * sqrt(diag(SD$cov.fixed)[names(diag(SD$cov.fixed))=="beta_mar"])))

load(files[2])
par_w <- c(exp(rep$mu_s)
         ,atan(rep$frho_j)*2/3.154
         ,atan(rep$frho2_jt)*2/3.154
         ,atan(rep$frho1_jt)*2/3.154
         ,exp(rep$fpsi_j)
         ,exp(rep$fpsi_jt)
         ,rep$beta_mar)

parSD_pos_w <- c(exp(rep$mu_s + 1.96 * sqrt(diag(SD$cov.fixed)["mu_s"]))
               ,atan(rep$frho_j + 1.96 * sqrt(diag(SD$cov.fixed)["frho_j"]))*2/3.154
               ,atan(rep$frho2_jt + 1.96 * sqrt(diag(SD$cov.fixed)["frho2_jt"]))*2/3.154
               ,atan(rep$frho1_jt + 1.96 * sqrt(diag(SD$cov.fixed)["frho1_jt"]))*2/3.154
               ,exp(rep$fpsi_j + 1.96 * sqrt(diag(SD$cov.fixed)["fpsi_j"]))
               ,exp(rep$fpsi_jt + 1.96 * sqrt(diag(SD$cov.fixed)["fpsi_jt"]))
               ,c(rep$beta_mar + 1.96 * sqrt(diag(SD$cov.fixed)[names(diag(SD$cov.fixed))=="beta_mar"])))

parSD_neg_w <- c(exp(rep$mu_s - 1.96 * sqrt(diag(SD$cov.fixed)["mu_s"]))
               ,atan(rep$frho_j - 1.96 * sqrt(diag(SD$cov.fixed)["frho_j"]))*2/3.154
               ,atan(rep$frho2_jt - 1.96 * sqrt(diag(SD$cov.fixed)["frho2_jt"]))*2/3.154
               ,atan(rep$frho1_jt - 1.96 * sqrt(diag(SD$cov.fixed)["frho1_jt"]))*2/3.154
               ,exp(rep$fpsi_j - 1.96 * sqrt(diag(SD$cov.fixed)["fpsi_j"]))
               ,exp(rep$fpsi_jt - 1.96 * sqrt(diag(SD$cov.fixed)["fpsi_jt"]))
               ,c(rep$beta_mar - 1.96 * sqrt(diag(SD$cov.fixed)[names(diag(SD$cov.fixed))=="beta_mar"])))

load(files[3])
par_h2 <- c(exp(rep$mu_s)
           ,1/(1+exp(rep$frho_j))
           ,atan(rep$frho2_jt)*2/3.154
           ,atan(rep$frho1_jt)*2/3.154
           ,exp(rep$fpsi_j)
           ,exp(rep$fpsi_jt)
           ,rep$beta_mar)

parSD_pos_h2 <- c(exp(rep$mu_s + 1.96 * sqrt(diag(SD$cov.fixed)["mu_s"]))
                 ,1/(1+exp(rep$frho_j - 1.96 * sqrt(diag(SD$cov.fixed)["frho_j"])))
                 ,atan(rep$frho2_jt + 1.96 * sqrt(diag(SD$cov.fixed)["frho2_jt"]))*2/3.154
                 ,atan(rep$frho1_jt + 1.96 * sqrt(diag(SD$cov.fixed)["frho1_jt"]))*2/3.154
                 ,exp(rep$fpsi_j + 1.96 * sqrt(diag(SD$cov.fixed)["fpsi_j"]))
                 ,exp(rep$fpsi_jt + 1.96 * sqrt(diag(SD$cov.fixed)["fpsi_jt"]))
                 ,c(rep$beta_mar + 1.96 * sqrt(diag(SD$cov.fixed)[names(diag(SD$cov.fixed))=="beta_mar"])))

parSD_neg_h2 <- c(exp(rep$mu_s - 1.96 * sqrt(diag(SD$cov.fixed)["mu_s"]))
                 ,1/(1+exp(rep$frho_j + 1.96 * sqrt(diag(SD$cov.fixed)["frho_j"])))
                 ,atan(rep$frho2_jt - 1.96 * sqrt(diag(SD$cov.fixed)["frho2_jt"]))*2/3.154
                 ,atan(rep$frho1_jt - 1.96 * sqrt(diag(SD$cov.fixed)["frho1_jt"]))*2/3.154
                 ,exp(rep$fpsi_j - 1.96 * sqrt(diag(SD$cov.fixed)["fpsi_j"]))
                 ,exp(rep$fpsi_jt - 1.96 * sqrt(diag(SD$cov.fixed)["fpsi_jt"]))
                 ,c(rep$beta_mar - 1.96 * sqrt(diag(SD$cov.fixed)[names(diag(SD$cov.fixed))=="beta_mar"])))

df <- data.frame(Hatchery=paste0(round(par,3)," ( ",round(parSD_neg,3),", ",round(parSD_pos,3)," )"),
                 Hatchery=paste0(round(par_h2,3)," ( ",round(parSD_neg_h2,3),", ",round(parSD_pos_h2,3)," )"),
                 Wild=paste0(round(par_w,3)," ( ",round(parSD_neg_w,3),", ",round(parSD_pos_w,3)," )"))

row.names(df) <- c("Mean annual survival"
                    ,"Correlation of day effect"
                    ,"Correlation of day in day/year effect"
                    ,"Correlation of year in day/year effect"
                    ,"Process error for day effect"
                    ,"Process error for day/year effect"
                    ,"Effect of first marine covariate"
                    ,"Effect of second marine covariate")
# names(par) <- c("Hatchery", "Wild")
write.csv(df,file="table_bestFitMods.csv", sep = "\t")
     