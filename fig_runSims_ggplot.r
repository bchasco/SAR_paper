#Plot the output
load("CorrelationAndSampleSizeExperiment.rData")
simMelt$Experiment <- "Manipulation"
simMelt$Experiment[simMelt$ex==3 & simMelt$ss==1] <- "theta[mle], gamma[mle]"
# simMelt <- as.data.frame(simMelt[simMelt$value<200 & simMelt$value>(-200),])
simMelt$pars <- as.character(simMelt$pars)

#Parameter expressions
simMelt$pars[simMelt$pars=="mu_s"] <- "mu"
simMelt$pars[simMelt$pars=="beta1"] <- "beta[CUI.spr]"
simMelt$pars[simMelt$pars=="beta2"] <- "beta[PDO.sum]"
simMelt$pars[simMelt$pars=="psi_j"] <- "psi[j]"
simMelt$pars[simMelt$pars=="psi_jt"] <- "phi[]"
simMelt$pars[simMelt$pars=="rho_j"] <- "rho[j]"
simMelt$pars[simMelt$pars=="rho2_jt"] <- "tau^(t)"
simMelt$pars[simMelt$pars=="rho1_jt"] <- "tau^(j)"
#Experimental expressions
simMelt$ex[simMelt$ex==3] <- "n[jt]"
simMelt$ex[simMelt$ex==2] <- "tau^(t)"
simMelt$ex[simMelt$ex==1] <- "rho[j]"

#Load the plotting libraries and plot
library(ggplot2)
library(wesanderson)
library(viridis)
p <- ggplot(simMelt, aes(factor(ss),value, fill=Experiment)) +
  facet_grid(vars(factor(pars)),vars(ex), scales = "free", labeller = label_parsed) +
  geom_violin(trim=FALSE,aes(factor(ss)))+
  ylab("Percent difference")+
  xlab("Trial")+
  theme_bw() +
  scale_fill_manual(values=plasma(n=2,alpha=0.3), 
                    labels=expression(paste(bold(theta)[paste(h,",",e)],",",bold(gamma)[paste(h,",",e)]),paste(bold(theta)[mle],",",bold(gamma)[mle])))+
  theme(strip.text.x = element_text(size = 14, colour = "black"))+
  theme(strip.text.y = element_text(size = 14, colour = "black"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        legend.text=element_text(size=c(14)),
        legend.title = element_blank(),
        legend.text.align = 0,
        legend.box.margin=margin(0,20,0,0))+
  theme(legend.key = element_rect(size = 6),
              legend.key.height = unit(1, "cm"),
              legend.key.width = unit(1, "cm"),
        legend.key.size = unit(30,"cm"))+
    geom_hline(aes(yintercept=0), 
             colour="blue",
             size=0.1)
print(p)
