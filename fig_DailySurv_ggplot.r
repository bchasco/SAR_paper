library(ggplot2)
rm(list=ls())

# tiff(paste0("fig_DailySurv_ggplot.tiff"), width=800, height=600)
wd <- "C:/Chasco/PROJECTS/SAR_PAPER/" #root directory
setwd(wd)

Hf <- "out_H_bestDailyMod.rData"
Wf <- "out_W_bestMod.rData"
load(Wf)
nj <- length(minJ:maxJ)
df <- rbind(data.frame(Origin = rep("Wild", nj),
                       day = minJ:maxJ,
                       surv=SD$value[names(SD$value)=="s_j"],
                       sd = SD$sd[names(SD$value)=="s_j"]))
df_n_w <- aggregate(list(total=data$s_n), by=list(day = data$j), sum)
df_n_w$total <- df_n_w$total/sum(df_n_w$total) 
df_n_w$Origin <- "Wild"


load(Hf)
df <- rbind(df,data.frame(Origin = rep("Hatchery", nj),
                       day = minJ:maxJ,
                       surv=SD$value[names(SD$value)=="s_j"],
                       sd = SD$sd[names(SD$value)=="s_j"]))

df_n_h <- aggregate(list(total=data$s_n), by=list(day = data$j), sum)
df_n_h$total <- df_n_h$total/sum(df_n_h$total) 
df_n_h$Origin <- "Hatchery"

df_n <- as.data.frame(rbind(df_n_w,
              df_n_h))

scale_fill_surv <- function(...){
  ggplot2:::manual_scale(
    'fill', 
    values = setNames(c('red', 'blue'), c("Hatchery","Wild")), 
    ...
  )
}

scale_color_surv <- function(...){
  ggplot2:::manual_scale(
    'color', 
    values = setNames(c('red', 'blue'), c("Hatchery","Wild")), 
    ...
  )
}

p <- ggplot(data=df, aes(x=as.Date(day, origin = as.Date("2018-01-01")), 
                         y = plogis(surv), 
                         group=factor(as.factor(Origin), levels=c("Wild","Hatchery")),
                         colour=Origin)) +
  geom_ribbon(aes(ymin=plogis(surv - 1.96 * sd),
                  ymax=plogis(surv + 1.96 * sd),
                  fill=Origin),
              alpha = 0.25,
              colour=NA)+
  scale_fill_surv() +
  scale_color_surv() +
  geom_line(size = 0.5,alpha = 1) +
  xlim(minJ,maxJ)+
  scale_x_date(date_labels = "%b") +
  theme(panel.background = element_blank()) +
  ylab("Marine survival") +
  xlab("")+
  theme(axis.text = element_text(size=rel(0.5))) + 
  theme(axis.title = element_text(size=rel(0.5))) +
  theme(legend.position=c(.8,.75),
        legend.title = element_blank(),
        legend.text = element_text(size=rel(0.5)))+
  theme(axis.line = element_line(colour = "black", size=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
# print(p)
p2 <- ggplot(data=df_n, aes(x=as.Date(day+minJ,origin = as.Date("2018-01-01")),
                            y=total, 
                            group=factor(as.factor(Origin), levels=c("Wild","Hatchery"))))+
  geom_ribbon(aes(ymin=rep(0,dim(df_n)[1]),
                  ymax=total,
                  fill=Origin),
              alpha=0.25)+
  scale_fill_surv() +
  xlim(minJ,maxJ) +
  scale_x_date(date_labels = "%b") +
  xlab("Bonneville passage date")+
  ylab("Migration rate") +
  theme(panel.background = element_blank()) +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size=rel(.5))) + 
  theme(axis.title = element_text(size=rel(.5))) +
  theme(legend.position="none",
        legend.title = element_blank(),
        legend.text = element_text(size=rel(0.5)))+
  theme(axis.line = element_line(colour = "black", size=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 


  
print(p)
print(p2)
              
o3 <- cowplot::plot_grid(p,p2, nrow = 2)

print(o3)

ggsave(plot = o3, file = "fig_DailySurvival_ggplot.png", 
       type = "cairo-png",  bg = "transparent",
       width = 6, height = 6, units = "cm", dpi = 400)

# dev.off()
