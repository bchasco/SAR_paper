library(ggplot2)
# rm(list=ls())

ylim <- 0.05
tiff("fig_DayxYearSurvival_ggplot.tiff", width=600, height=800)
wd <- "C:/Chasco/PROJECTS/SAR_PAPER/" #root directory
setwd(wd)
Hf <- "out_H_bestMod.rData"
Wf <- "out_W_bestMod.rData"
rearType <- "H"
load(Hf)
nt <- length(startYr:lastYr)
nj <- length(minJ:maxJ)
Hdf <- data.frame(origin = rep("H", length(startYr:lastYr)*length(minJ:maxJ)),
                                yr = rep(startYr:lastYr, each = length(minJ:maxJ)),
                                j = rep(minJ:maxJ, length(startYr:lastYr)),
                                surv=c(t(matrix(SD$value[names(SD$value)=="s_tj"],nt,nj))),
                                sd = c(t(matrix(SD$sd[names(SD$value)=="s_tj"],nt,nj))))

ag_h <- aggregate(list(total=data$s_n, surv=data$s_k), by=list(yr = data$yr + startYr), sum)
ag_h$origin <- "H"

timing <- data.frame(origin = rep("H", length(data$yr[(data$s_k/data$s_n)<ylim])),
                      yr = data$yr[(data$s_k/data$s_n)<ylim] + startYr,
                     j = data$j[(data$s_k/data$s_n)<ylim]+minJ,
                     s_n = data$s_n[(data$s_k/data$s_n)<ylim],
                     s_k = data$s_k[(data$s_k/data$s_n)<ylim])

load(Wf)
Wdf <- data.frame(origin = rep("W", length(startYr:lastYr)*length(minJ:maxJ)),
                  yr = rep(startYr:lastYr, each = length(minJ:maxJ)),
                  j = rep(minJ:maxJ, length(startYr:lastYr)),
                  surv=c(t(matrix(SD$value[names(SD$value)=="s_tj"],nt,nj))),
                  sd = c(t(matrix(SD$sd[names(SD$value)=="s_tj"],nt,nj))))

ag_w <- aggregate(list(total=data$s_n, surv=data$s_k), by=list(yr = data$yr + startYr), sum)
ag_w$origin <- "W"


timing <- rbind(timing,
                data.frame(origin = rep("W", length(data$yr[(data$s_k/data$s_n)<ylim])),
                           yr = data$yr[(data$s_k/data$s_n)<ylim] + startYr,
                     j = data$j[(data$s_k/data$s_n)<ylim]+minJ,
                     s_n = data$s_n[(data$s_k/data$s_n)<ylim],
                      s_k = data$s_k[(data$s_k/data$s_n)<ylim]))
df_n <- rbind(ag_h,ag_w)
for(i in unique(df_n$yr)){
  df_n$ylim[df_n$yr==i] <- max(plogis(df$surv[df$yr==i]+2*df$sd[df$yr==i] ))
}
df <- rbind(Hdf,Wdf)

# myYears <- c(2002,2003,2009,2010,2012,2013)
myYears <- startYr:lastYr
df_n <- df_n[df_n$yr%in%myYears,]
df <- df[df$yr%in%myYears,]
timing <- timing[timing$yr%in%myYears,]


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

p <- ggplot(data=df, aes(x=as.Date(j, origin = as.Date("2018-01-01")), y = plogis(surv), factor(origin))) +
  geom_line(size = 0.2,alpha = 0.5, aes(colour=factor(origin))) + 
  geom_ribbon(aes(ymin=plogis(surv - 1.96 * sd), 
                  ymax=plogis(surv + 1.96 * sd), 
                  # colour = factor(origin),
                  fill=factor(origin)), 
              alpha = 0.15)  + 
  scale_fill_surv() +
  scale_color_surv() +
  scale_x_date(date_labels = "%b") + 
  facet_wrap(~yr,scales = "free_y") +
  xlim(110,160) + 
  # ylim(0,0.3) + 
  theme(panel.background = element_blank()) + 
  # ylim(0,0.2) +
  ylab("Marine survival (ocean entry year)") + 
  xlab("Calendar day") +
  theme(axis.line = element_line(colour = "black", size=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 


p <- p + 
  geom_point(data=timing, 
             aes(x=as.Date(j, origin = as.Date("2018-01-01")), y = s_k/s_n, 
                 colour = factor(origin)),
             shape = 16,
             size = 2.5, alpha = 0.25) +
  scale_fill_surv() +
  scale_color_surv() +
  scale_x_date(date_labels = "%b") + 
  facet_wrap(~yr, scales = "free_y") +
  theme(legend.position = "none")
  # scale_y_continuous(
  #   "PIT observations per day", 
  #   sec.axis = sec_axis(~ . * 0.0001, name = "Survival")
  # )
df_n$label <- paste(df_n$origin," = ", df_n$surv,"/",df_n$total)
# p <- p + 
#   geom_text(data = df_n[df_n$origin=="W",], 
#             aes(label = label, x=as.Date(150, origin = as.Date("2018-01-01")), y = 0.1))
geom.text.size = 4
theme.size = (17/5) * geom.text.size

p <- p + 
  geom_text(data = df_n[df_n$origin=="W",], 
            aes(label = label, x=as.Date(150, origin = as.Date("2018-01-01")),y=ylim), size=geom.text.size)
p <- p + 
  geom_text(data = df_n[df_n$origin=="H",], 
            aes(label = label, x=as.Date(150, origin = as.Date("2018-01-01")),y=ylim*0.9), size=geom.text.size)

print(p)
dev.off()
