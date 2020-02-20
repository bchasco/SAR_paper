library(ggplot2)
rm(list=ls())

png(paste0("fig_envEffect_ggplot.png"), width=6, height=6, units = "in", res=600)
wd <- "C:/Chasco/PROJECTS/SAR_PAPER/" #root directory
setwd(wd)
Hf <- "out_H_bestMod.rData"
Wf <- "out_W_bestMod.rData"
load(Hf)
load(Wf)
h_env <- c(11,24)
w_env <- c(15,36)

load("envData.rData")
marVars <- c(names(envdata)[2:13],names(envdata)[18:42])
h_names <- names(subData[h_env])
w_names <- names(subData[w_env])

dev <- seq(-1.96,1.96,0.05)
nd <- length(dev)
beta <- c(plogis(rep$mu_s+cbind(dev)%*%rep$beta_mar[1,1]),
                 plogis(rep$mu_s+cbind(dev)%*%rep$beta_mar[2,1]))
beta <- (beta - plogis(rep$mu_s))/plogis(rep$mu_s)*100
df <- data.frame(Origin = rep("Wild", nd),
                       dev = dev,
                 beta = rep(w_names,each=nd),
                       surv=beta)

load(Hf)
dev <- seq(-1.96,1.96,0.1)
nd <- length(dev)
beta <- c(plogis(rep$mu_s+cbind(dev)%*%rep$beta_mar[1,1]),
          plogis(rep$mu_s+cbind(dev)%*%rep$beta_mar[2,1]))
beta <- (beta - plogis(rep$mu_s))/plogis(rep$mu_s)*100
df <- rbind(df,data.frame(Origin = rep("Hatchery", nd),
                 dev = dev,
                 beta = rep(h_names,each=nd),
                 surv=beta))



p <- ggplot(data=df[df$Origin=="Wild",], aes(x=dev, 
                         y = surv, 
                         group=beta,
                         colour=beta))+
  geom_line(size = 0.5,alpha = 1)+ 
  theme(panel.background = element_blank(),
        axis.title=element_text(size=8)) + 
  ylab("% change in survival") + 
  xlab("")+
  theme(legend.position=c(.5,.75),
        legend.title = element_blank())+
  ylim(-100,200) + 
  theme(axis.text = element_text(size=rel(1.5))) + 
  theme(axis.title = element_text(size=rel(1.5))) +
  theme(legend.position=c(.5,.75),
        legend.title = element_blank(),
        legend.text = element_text(size=rel(1.5)))+
  theme(axis.line = element_line(colour = "black", size=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

p2 <- ggplot(data=df[df$Origin=="Hatchery",], aes(x=dev, 
                                             y = surv, 
                                             group=beta,
                                             colour=beta))+
  geom_line(size = 0.5,alpha = 1)+ 
  theme(panel.background = element_blank()) + 
  ylab("% change in survival") + 
  xlab("Environmental deviate")+
  theme(legend.position=c(.5,.75),
        legend.title = element_blank(),
        axis.title=element_text(size=8))+
  ylim(-100,200) + 
  theme(axis.text = element_text(size=rel(1.5))) + 
  theme(axis.title = element_text(size=rel(1.5))) +
  theme(legend.position=c(.5,.75),
        legend.title = element_blank(),
        legend.text = element_text(size=rel(1.5)))+
  theme(axis.line = element_line(colour = "black", size=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

  
print(p)
print(p2)
              
o3 <- cowplot::plot_grid(p,p2, nrow = 2)
print(o3)

ggsave(plot = o3, file = "EnvironmentalEffect.png", 
       type = "cairo-png",  bg = "transparent",
       width = 8, height = 8, units = "cm", dpi = 300)

dev.off()
