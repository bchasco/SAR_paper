library(ggplot2)

load("H.rData")
load("W.rData")

df <- data.frame(Origin=c(rep("Hatchery",nrow(H)*2),rep("Wild",nrow(W)*2)),
                        var=c(H$firstVar,H$secondVar,W$firstVar,W$secondVar),
                        wt=c(exp(-H$deltaAIC),exp(-H$deltaAIC),exp(-W$deltaAIC),exp(-W$deltaAIC)))
df$var <- as.character(df$var)
df$var[is.na(df$var)] <- "blank"
df <- aggregate(df$wt,by=list(var=df$var,Origin=df$Origin),sum)
df$x[df$Origin=="Wild"] <- df$x[df$Origin=="Wild"]/sum(df$x[df$Origin=="Wild"])
df$x[df$Origin=="Hatchery"] <- df$x[df$Origin=="Hatchery"]/sum(df$x[df$Origin=="Hatchery"])

scale_fill_surv <- function(...){
  ggplot2:::manual_scale(
    'fill', 
    values = setNames(c('red', 'blue'), c("Hatchery","Wild")), 
    ...
  )
}

p<- ggplot(data=df, aes(x=as.factor(var),y=x,  fill=Origin), alpha=0.5) +
  facet_grid(.~Origin)+
  coord_flip() +
  scale_fill_surv() +
  scale_alpha(range = c(0.5)) + 
  xlab("Environmental variables") +
  ylab("Aggregated AIC wt") +
  geom_bar(stat="identity", position = position_dodge(), alpha=0.5) +
  theme_bw() +
  theme(legend.position = "none")
  theme(axis.text.x = element_text(size=0.2),
        axis.text.y = element_text(size=0.2))

print(p)

ggsave(plot = p, file = "fig_EnvironmentalVariableWt_ggplot.png", 
       type = "cairo-png",  bg = "transparent",
       width = 8, height = 8, units = "in", dpi = 400)

