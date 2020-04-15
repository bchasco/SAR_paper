#Set seed to set.seed(800)

p <- ggplot(dfAll[dfAll$val<0.4,],aes(x=as.Date(j, origin = as.Date("2018-01-01")),y=val))+
  scale_x_date(date_labels = "%b") +
  facet_wrap(~yr,scale="free")+
  xlim(as.Date(minJ+5, origin = as.Date("2018-01-01")),as.Date(maxJ+5, origin = as.Date("2018-01-01")))+
  # ylim(0,0.25)+
  geom_point(data=dfAll[dfAll$out=="obs" & dfAll$val<0.4,],
             aes(as.Date(j, origin = as.Date("2018-01-01")),val), color=rgb(0.5,0.5,0.5,0.5))+
  geom_line(data=dfAll[dfAll$out=="MV" & dfAll$val<0.4,],
            aes(as.Date(j, origin = as.Date("2018-01-01")),val), colour="blue", size=0.9)+
  geom_line(data=dfAll[dfAll$out=="fixed" & dfAll$val<0.4,],
            aes(as.Date(j, origin = as.Date("2018-01-01")),val), colour="red", size=0.9) +
  theme_bw() +
  ylab("Marine survival")+
  xlab("Calendar day")+
  # scale_fill_manual(values=plasma(n=3,alpha=0.3))+
  theme(strip.text.x = element_text(size = 10, colour = "black"))+
  theme(strip.text.y = element_text(size = 10, colour = "black"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))
  #   # geom_line(data=df[df$out=="MV_mean",],aes(j,val), colour="darkgreen", size=0.9)
print(p)
