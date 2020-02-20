library(ggplot2)
rm(list=ls())

tiff(paste0("fig_AnnualSurv_ggplot.tiff"), width=800, height=600)
wd <- "C:/NOAA/PROJECTS/COLLEAGUES/BURKE/SAR_PAPER/"
setwd(wd)

Hf <- "out_H_bestMod.rData"
Wf <- "out_W_bestMod.rData"
load(Wf)
nt <- length(startYr:lastYr)
nj <- length(minJ:maxJ)
df <- rbind(data.frame(Origin = rep("Wild", nt),
                       yr = startYr:lastYr,
                       surv=SD$value[names(SD$value)=="s_t"],
                       sd = SD$sd[names(SD$value)=="s_t"]))
df_n <- aggregate(list(total=data$s_n, surv=data$s_k), by=list(yr = data$yr + startYr - data$yShift), sum)


load(Hf)
df <- rbind(df,
            data.frame(Origin = rep("Hatchery", nt),
                 yr = startYr:lastYr,
                 surv=SD$value[names(SD$value)=="s_t"],
                 sd = SD$sd[names(SD$value)=="s_t"]))

df_n <- rbind(df_n,
              aggregate(list(total=data$s_n, surv=data$s_k), by=list(yr = data$yr + startYr - data$yShift), sum))
df_n$Origin <- rep(c("Wild","Hatchery"),each=nt)



# myYears <- c(2002,2003,2009,2010,2012,2013)
# myYears <- startYr:lastYr
# df_n <- df_n[df_n$yr%in%myYears,]
df_n$ylim <- aggregate(df$surv+1.96*df$sd,by=list(yr=df$yr, Origin=df$Origin),max)$x

ylim.all <- max(plogis(df_n$ylim),df_n$surv/df_n$total)
# df <- df[df$yr%in%myYears,]


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

p <- ggplot(data=df, aes(x=yr, y = plogis(surv), 
                         group=factor(Origin),
                         colour=Origin)) +
  geom_line(size = 0.75,alpha = 1) + 
  geom_ribbon(aes(ymin=plogis(surv - 1.96 * sd),
                  ymax=plogis(surv + 1.96 * sd),
                  fill=Origin),
              alpha = 0.15,
              colour=NA)+
  # facet_wrap(~Origin, nrow=2)+
  scale_fill_surv() +
  scale_color_surv() +
  ylim(0,ylim.all) +
  xlim(2000,2015) +
  theme(panel.background = element_blank()) + 
  theme_bw()+
  ylab("Marine survival (ocean entry year)") + 
  xlab("Migration year")

p <- p + 
  geom_point(data=df_n, 
             aes(x=yr, y = surv/total, 
                 group=as.factor(Origin)),
             shape = 16,
             size = 4) +
  theme(axis.text = element_text(size=rel(2))) + 
  theme(axis.title = element_text(size=rel(2))) +
  theme(legend.position=c(.8,.75),
        legend.title = element_blank(),
        legend.text = element_text(size=rel(2))) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

print(p)

dev.off()
