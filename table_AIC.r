load("H.Rdata")
load("W.Rdata")
W <- data.frame(W[,1:11],firstVar=W$firstVar,secondVar=W$secondVar)

df <- rbind(H[H$deltaAIC<4,],
            W[W$deltaAIC<4,])
write.csv(df,"table_AICoutput.csv")
