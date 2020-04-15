library(ggplot2)
library(ggpubr)

load(file="simQuadractic.rData")

label_panel <- aggregate(list(value=df$value),by=list(stat=df$stat),max)
label_panel$label <- LETTERS[1:nrow(label_panel)]
label_panel$Model <- "TMB"

med_points <- aggregate(list(value=df$value),by=list(pars=df$pars,stat=df$stat),median)
line_0 <- aggregate(list(value=df$value),by=list(pars=df$pars,stat=df$stat,Model=df$Model),median)
line_0 <- line_0[line_0$stat=="A",]
line_0$line_0 <- 0

pt_se <- aggregate(list(value=df$value),by=list(pars=df$pars,stat=df$stat,Model=df$Model),median)
pt_se <- pt_se[pt_se$stat=="B",]

#https://stackoverflow.com/questions/51228076/ggplot-split-violin-plot-with-horizontal-mean-lines
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, draw_group = function(self, data, ..., draw_quantiles = NULL){
  data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
  grp <- data[1,'group']
  newdata <- plyr::arrange(transform(data, x = if(grp%%2==1) xminv else xmaxv), if(grp%%2==1) y else -y)
  newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
  newdata[c(1,nrow(newdata)-1,nrow(newdata)), 'x'] <- round(newdata[1, 'x']) 
  if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
    stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 
                                              1))
    quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
    aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
    aesthetics$alpha <- rep(1, nrow(quantiles))
    both <- cbind(quantiles, aesthetics)
    quantile_grob <- GeomPath$draw_panel(both, ...)
    ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
  }
  else {
    ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
  }
})

geom_split_violin <- function (mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, position = position, show.legend = show.legend, inherit.aes = inherit.aes, params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

#Here's the code for t

group_name = c(expression(mu),
               expression(beta[CUI.spr]),
               expression(beta[PDO.sum]))

p <- ggplot(df[df$stat=="A",],aes(pars,value, fill=Model))+
  geom_split_violin(trim=TRUE) +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", 
               width = 0.25,
               position = position_dodge(width = .25)
  ) +
  ylab("Percent difference")+
  xlab("Trial")+
  theme_bw() +
  xlab("")+
  scale_fill_manual(values=plasma(n=3,alpha=0.3))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))+
  scale_x_discrete(labels=group_name)+
  geom_hline(data=line_0,
             aes(yintercept=line_0),
             colour="blue",
             size=1) +
  theme(text = element_text(size=16))
  
  
p2 <- ggplot(df[df$stat=="B",],aes(pars,value, fill=Model))+
  geom_split_violin(trim=TRUE) +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", 
               width = 0.25,
               position = position_dodge(width = .25)
  ) +
  ylab("Standard error")+
  xlab("Trial")+
  theme_bw() +
  xlab("")+
  scale_x_discrete(labels=group_name)+
  scale_fill_manual(values=plasma(n=3,alpha=0.3))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))+
  theme(text = element_text(size=16))
  

figure <- ggarrange(p,p2,
                    labels = c("A", "B"),
                    nrow = 2,
                    label.x = c(1,1),
                    label.y = c(1,1),
                    hjust = 5,
                    vjust = 3,
                    # fig.lab.pos = "top.right",
                    common.legend=TRUE,
                    legend="right")

figure


