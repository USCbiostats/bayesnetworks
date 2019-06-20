library(pcalg) 
mydir <- "/Users/dthomas/Duncan/Genetics/Stat Gen P01/Project 3/Resubmission"
setwd(mydir)
dat <- read.table("P3 simulation 2.dat")

dat2 = dat[,c(2:13,24:33,44:53)]
pc.fit2 = pc(list(C=cor(dat2), n=nrow(dat2)), labels=paste(1:32), indepTest = gaussCItest, alpha=0.01)
plot(pc.fit2,main="fit32")

dat4a = dat[,c(2:13,24:33,44:53,64:73)]
pc.fit4a = pc(list(C=cor(dat4a), n=nrow(dat4a)), labels=paste(1:42), indepTest = gaussCItest, alpha=0.01)
plot(pc.fit4a,main="fit42")

dat4 = dat[,c(2:13,24:38,44:53,64:73)]
pc.fit4 = pc(list(C=cor(dat4), n=nrow(dat4)), labels=paste(1:47), indepTest = gaussCItest, alpha=0.01)
plot(pc.fit4,main="fit47")

dat5 = dat[,c(2:13,24:33,44:53,64:83)]
pc.fit5 = pc(list(C=cor(dat5), n=nrow(dat5)), labels=paste(1:52), indepTest = gaussCItest, alpha=0.01)
plot(pc.fit5,main="fit52")

dat5a = dat[,c(2:13,24:38,44:53,64:83)]
pc.fit5a = pc(list(C=cor(dat5a), n=nrow(dat5a)), labels=paste(1:57), indepTest = gaussCItest, alpha=0.01)
plot(pc.fit5a,main="fit57")

dat6 = dat[,c(2:13,24:33,44:53,64:83)]
pc.fit6 = pc(list(C=cor(dat6), n=nrow(dat6)), labels=paste(1:62), indepTest = gaussCItest, alpha=0.01)
plot(pc.fit6,main="fit62")

dat6a = dat[,c(2:13,24:38,44:53,64:83)]
pc.fit6a = pc(list(C=cor(dat6a), n=nrow(dat6a)), labels=paste(1:67), indepTest = gaussCItest, alpha=0.01)
plot(pc.fit6a,main="fit67")

dat7 = dat[,c(2:73)]
pc.fit7 = pc(list(C=cor(dat7), n=nrow(dat7)), labels=paste(1:72), indepTest = gaussCItest, alpha=0.01)
plot(pc.fit7,main="fit72")

dat7 = dat[,c(2:78)]
pc.fit7 = pc(list(C=cor(dat7), n=nrow(dat7)), labels=paste(1:77), indepTest = gaussCItest, alpha=0.01)
plot(pc.fit7,main="fit77")

dat8 = dat[,c(2:83)]
pc.fit8 = pc(list(C=cor(dat8), n=nrow(dat8)), labels=paste(1:82), indepTest = gaussCItest, alpha=0.01)
plot(pc.fit8,main="fit82")


dat2a = dat[,c(3:13,24:33,44:53,64:73)]
pc.fit2a = pc(list(C=cor(dat2a), n=nrow(dat2a)), labels=paste(1:41), indepTest = gaussCItest, alpha=0.01)
plot(pc.fit2a,main="fit41")
