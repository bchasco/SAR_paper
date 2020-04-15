library(ggplot2)
# rm(list=ls())

ylim <- 0.05
tiff("fig_DayxYearSurvival_Alt_ggplot.tiff", width=600, height=800)
wd <- "C:/NOAA/PROJECTS/SAR_PAPER/" #root directory
setwd(wd)
Hf <- "out_H_bestMod.rData"
Wf <- "out_W_bestMod.rData"
rearType <- "H"
load(Hf)
nt <- length(startYr:lastYr)
nj <- length(minJ:maxJ)
Hdf <- data.frame(origin2 = rep("Hatchery", length(startYr:lastYr)*length(minJ:maxJ)),
                                yr = rep(startYr:lastYr, each = length(minJ:maxJ)),
                                j = rep(minJ:maxJ, length(startYr:lastYr)),
                                surv=c(t(matrix(SD$value[names(SD$value)=="s_tj"],nt,nj))),
                                sd = c(t(matrix(SD$sd[names(SD$value)=="s_tj"],nt,nj))))

ag_h <- aggregate(list(total=data$s_n, surv=data$s_k), by=list(yr = data$yr + startYr), sum)
ag_h$origin2 <- "H"

timing <- data.frame(origin2 = rep("Hatchery", length(data$yr[(data$s_k/data$s_n)<ylim])),
                      yr = data$yr[(data$s_k/data$s_n)<ylim] + startYr,
                     j = data$j[(data$s_k/data$s_n)<ylim]+minJ,
                     s_n = data$s_n[(data$s_k/data$s_n)<ylim],
                     s_k = data$s_k[(data$s_k/data$s_n)<ylim])

load(Wf)
Wdf <- data.frame(origin2 = rep("Wild", length(startYr:lastYr)*length(minJ:maxJ)),
                  yr = rep(startYr:lastYr, each = length(minJ:maxJ)),
                  j = rep(minJ:maxJ, length(startYr:lastYr)),
                  surv=c(t(matrix(SD$value[names(SD$value)=="s_tj"],nt,nj))),
                  sd = c(t(matrix(SD$sd[names(SD$value)=="s_tj"],nt,nj))))

ag_w <- aggregate(list(total=data$s_n, surv=data$s_k), by=list(yr = data$yr + startYr), sum)
ag_w$origin2 <- "Wild"


df <- rbind(Hdf,Wdf)

# myYears <- c(2002,2003,2009,2010,2012,2013)
myYears <- startYr:lastYr
df <- df[df$yr%in%myYears,]


scale_fill_surv <- function(...){
  ggplot2:::manual_scale(
    'fill', 
    values = setNames(c('red', 'blue'), c("H","W")), 
    ...
  )
}

scale_color_surv <- function(...){
  ggplot2:::manual_scale(
    'color', 
    values = setNames(c('red', 'blue'), c("H","W")), 
    ...
  )
}

p <- ggplot(data=df, aes(x=j, 
                         y = plogis(surv), factor(yr))) +
  geom_line(size = 0.2,alpha = 0.5, aes(group=yr)) + 
  scale_x_date(date_labels = "%b") +
  facet_wrap(~origin2,scales = "free_y", nrow = 2) +
  scale_color_identity()+
  xlim(110,160) +
  # ylim(0,0.3) +
  theme(panel.background = element_blank()) +
  # ylim(0,0.2) +
  ylab("Marine survival") +
  xlab("Calendar day") +
  theme(axis.line = element_line(colour = "black", size=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

print(p)
dev.off()
