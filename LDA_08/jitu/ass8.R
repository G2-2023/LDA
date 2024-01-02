########################################
########################################
# MATHDAT WITH MIXTURES
########################################
install.packages("R2jags")
library(R2jags)
getwd()

###############################################################
# DATA
###############################################################
dat2 <- read.table("/cloud/project/LDA_08/jitu/autism_short.dat", header=T, 
                   na.strings="NA", dec=".",strip.white = T, 
                   col.names = c("id","group","y1","y2","y3","y4","y5"))

head(dat2)

# restructure data into short format
N <- length(levels(as.factor(dat2$id))) # no of persons
dim(dat2)
#Todo from here
dat1 <- data.frame(matrix(NA,N,7))
for(i in 1:N){#i<-1
  idat <- dat2[dat2$id==i,]
  dat1[i,1:6] <- idat$dep
  dat1[i,7] <- idat$order[1]
}

colnames(dat1) <- c(paste0("dep",1:6),"treat")
head(dat1)

N <- dim(dat1)[1]

# illustrate the data
plot(dat1[,paste0("dep",1:6)])

# or use a Spagettiplot
colo <- c("black","red")
plot(1:6,apply(dat1[dat1$treat==1,paste0("dep",1:6)],2,mean,na.rm=T),type="l",lwd=2,ylim=c(0,6),ylab="dep",xlab="time",col=colo[2])
par(new=T)
plot(1:6,apply(dat1[dat1$treat==0,paste0("dep",1:6)],2,mean,na.rm=T),type="l",lwd=2,ylim=c(0,6),ylab="dep",xlab="time",col=colo[1])
for(i in 1:N){
  par(new=T)
  plot(1:6,jitter(as.numeric(dat1[i,paste0("dep",1:6)])),type="l",ylim=c(0,6),ylab="",xlab="",col=colo[dat1$treat[i]+1],lty=2)
}


###############################################################
# PART 2: THIS IS OUR BASELINE MODEL IN JAGS
# USE CODE from lgm_jags_bock.txt
###############################################################

y <- as.matrix(dat1[,paste0("dep",1:6)])
data.jags <- list(y=y,x=dat1$treat+1,N=N,psi0=diag(3))
params <- c("beta0","sigmay","sigmaeta")

###############################################################
# RUN MODELS
###############################################################
model1 <- jags.parallel(data=data.jags, 
                        parameters.to.save=params,
                        n.iter=4000, n.chains=3,n.thin=1,n.burnin = 2000,
                        model.file="lgm_jags_bock.txt")

est1 <- model1$BUGSoutput$summary
round(est1,2)

#############################################################
# traceplots
samps <- as.mcmc(model1)
#xyplot, densityplot
pdf(file="xyplot.pdf",pointsize=6)
plot(samps)
dev.off()

# Rhat statistic over time (iterations)
pdf("gelmanplot.pdf")
par(new=F)
for (v in 1:nvar(samps)){
  gelman.plot(samps[,v],ylim=c(0,3),col=c("black",NA))
  par(new=T)
}
abline(h=1.1)
dev.off()
###############################################################

# alternatively here is the correctly specified model
# that includes the treatment and the same code as the initial model
# from exercise 6 (model1.piece)

model2 <- jags.parallel(data=data.jags, 
                        parameters.to.save=params,
                        n.iter=4000, n.chains=3,n.thin=1,n.burnin = 2000,
                        model.file="lgm_jags_bock2.txt")

est2 <- model2$BUGSoutput$summary
round(est2,2)

# works totally fine.
# we now constrain the relevant mean and covariance structures.
# (including residual variances)
model3 <- jags.parallel(data=data.jags, 
                        parameters.to.save=params,
                        n.iter=4000, n.chains=3,n.thin=1,n.burnin = 2000,
                        model.file="lgm_jags_bock3.txt")

est3 <- model3$BUGSoutput$summary
round(est3,2)





###############################################################
# PART 3: THIS IS OUR MIXTURE MODEL IN JAGS
# We first use a naive model, then one relating to our model3,
# the last model is specific application of piecewise models
# USE CODE from gmm_jags_bock.txt and expand to the different versions
###############################################################

params <- c("beta0","deltabeta","sigmay","sigmaeta","alpha","prob")

###############################################################
# check convergence with 20k its
###############################################################
model4 <- jags.parallel(data=data.jags, 
                        parameters.to.save=params,
                        n.iter=20000, n.chains=3,n.thin=1,n.burnin = 1,
                        model.file="gmm_jags_bock.txt")

est4 <- model4$BUGSoutput$summary
round(est4,2)

###############################################################
# traceplots
samps <- as.mcmc(model4)
#xyplot, densityplot
pdf(file="xyplot.pdf",pointsize=6)
plot(samps)
dev.off()

# Rhat statistic over time (iterations)
pdf("gelmanplot.pdf")
par(new=F)
for (v in 1:nvar(samps)){
  gelman.plot(samps[,v],ylim=c(0,3),col=c("black",NA))
  par(new=T)
}
abline(h=1.1)
dev.off()
###############################################################

###############################################################
# check model estimates with 10k burnin, 10k additional its
###############################################################
model5 <- jags.parallel(data=data.jags, 
                        parameters.to.save=params,
                        n.iter=10000, n.chains=3,n.thin=1,n.burnin = 5000,
                        model.file="gmm_jags_bock.txt")

est5 <- model5$BUGSoutput$summary
round(est5,2)

# crude plot of the results
plot(0,0,col="white",xlim=c(0,5),ylim=c(0,6),xlab="time",ylab="dep")
# class #1
abline(est5["beta0[1,1]","mean"],est5["beta0[2,1]","mean"])
abline(est5["beta0[1,1]","mean"]+sqrt(est5["sigmaeta[1,1,1]","mean"]),
       est5["beta0[2,1]","mean"]+sqrt(est5["sigmaeta[2,2,1]","mean"]),lty=2)
abline(est5["beta0[1,1]","mean"]-sqrt(est5["sigmaeta[1,1,1]","mean"]),
       est5["beta0[2,1]","mean"]-sqrt(est5["sigmaeta[2,2,1]","mean"]),lty=2)
# class #2
abline(est5["beta0[1,2]","mean"],est5["beta0[2,2]","mean"],col="red")
abline(est5["beta0[1,2]","mean"]+sqrt(est5["sigmaeta[1,1,2]","mean"]),
       est5["beta0[2,2]","mean"]+sqrt(est5["sigmaeta[2,2,2]","mean"]),lty=2,col="red")
abline(est5["beta0[1,2]","mean"]-sqrt(est5["sigmaeta[1,1,2]","mean"]),
       est5["beta0[2,2]","mean"]-sqrt(est5["sigmaeta[2,2,2]","mean"]),lty=2,col="red")


par(new=T)
plot(1:6,apply(dat1[dat1$treat==1,paste0("dep",1:6)],2,mean,na.rm=T),type="l",lwd=2,ylim=c(0,6),ylab="dep",xlab="time",col=colo[2],lty=3)
par(new=T)
plot(1:6,apply(dat1[dat1$treat==0,paste0("dep",1:6)],2,mean,na.rm=T),type="l",lwd=2,ylim=c(0,6),ylab="dep",xlab="time",col=colo[1],lty=3)

legend("bottomleft",c("mean","+/-1sd","observed"),lty=1:3)
# this somewhat resembles the grouping of the data

###############################################################
# include a partly informed model
# The researcher knows that persons switched between time point 3 and 4.
# She also assumes that the changes are parallel
# She only forgot which person was in which group (e.g., anonymous data)
# -> This model is fairly identical to lgm_jags_bock3.txt above
###############################################################
model6 <- jags.parallel(data=data.jags, 
                        parameters.to.save=params,
                        n.iter=20000, n.chains=3,n.thin=1,n.burnin = 1,
                        model.file="gmm_jags_bock3.txt")

est6 <- model6$BUGSoutput$summary
round(est6,2)

###############################################################
# traceplots
samps <- as.mcmc(model6)
#xyplot, densityplot
pdf(file="xyplot.pdf",pointsize=6)
plot(samps)
dev.off()

# Rhat statistic over time (iterations)
pdf("gelmanplot.pdf")
par(new=F)
for (v in 1:nvar(samps)){
  gelman.plot(samps[,v],ylim=c(0,3),col=c("black",NA))
  par(new=T)
}
abline(h=1.1)
dev.off()
###############################################################




###############################################################
# include class specific residual variances
###############################################################
model7 <- jags.parallel(data=data.jags, 
                        parameters.to.save=params,
                        n.iter=20000, n.chains=3,n.thin=1,n.burnin = 10000,
                        model.file="gmm_jags_bock3.txt")

est7 <- model7$BUGSoutput$summary
round(est7,2)
# compare to 
round(est2,2)
# looks rather similar


###############################################################
# include no class specific residual variances
###############################################################
model8 <- jags.parallel(data=data.jags, 
                        parameters.to.save=params,
                        n.iter=20000, n.chains=3,n.thin=1,n.burnin = 10000,
                        model.file="gmm_jags_bock2.txt")

est8 <- model8$BUGSoutput$summary
round(est8,2)
# compare to 
round(est3,2)
# looks rather similar



###############################################################
# extract class-memberships
###############################################################
params2 <- ("C")
model7.c <- jags.parallel(data=data.jags, 
                        parameters.to.save=params2,
                        n.iter=20000, n.chains=3,n.thin=1,n.burnin = 10000,
                        model.file="gmm_jags_bock3.txt")

model8.c <- jags.parallel(data=data.jags, 
                          parameters.to.save=params2,
                          n.iter=20000, n.chains=3,n.thin=1,n.burnin = 10000,
                          model.file="gmm_jags_bock2.txt")

# estimated class membership (median)
C.est7 <- model7.c$BUGSoutput$summary[-(N+1),"50%"]
C.est8 <- model8.c$BUGSoutput$summary[-(N+1),"50%"]

# comparison to observed grouping
table(C.est7,dat1$treat) # group-specific variances
table(C.est8,dat1$treat) # no group-specific variances
table(C.est7,C.est8)


colo <- c("black","blue")
# average trend
plot(1:6,apply(dat1[C.est7==1,paste0("dep",1:6)],2,mean,na.rm=T),type="l",lwd=2,ylim=c(0,6),ylab="dep",xlab="time")
par(new=T)
plot(1:6,apply(dat1[C.est7==2,paste0("dep",1:6)],2,mean,na.rm=T),type="l",lwd=2,ylim=c(0,6),ylab="",xlab="time",col="blue")
# individual slopes
for(i in 1:N){
  par(new=T)
  plot(1:6,jitter(as.numeric(dat1[i,paste0("dep",1:6)])),type="l",ylim=c(0,6),ylab="",xlab="",col=colo[C.est7[i]],lty=2)
}

# or version 8
# average trend
plot(1:6,apply(dat1[C.est8==1,paste0("dep",1:6)],2,mean,na.rm=T),type="l",lwd=2,ylim=c(0,6),ylab="dep",xlab="time")
par(new=T)
plot(1:6,apply(dat1[C.est8==2,paste0("dep",1:6)],2,mean,na.rm=T),type="l",lwd=2,ylim=c(0,6),ylab="",xlab="time",col="blue")
# individual slopes
for(i in 1:N){
  par(new=T)
  plot(1:6,jitter(as.numeric(dat1[i,paste0("dep",1:6)])),type="l",ylim=c(0,6),ylab="",xlab="",col=colo[C.est8[i]],lty=2)
}

# or separate for each observed group
par(mfrow=c(1,2))

plot(1:6,apply(dat1[C.est8==1,paste0("dep",1:6)],2,mean,na.rm=T),type="l",lwd=2,ylim=c(0,6),ylab="dep",xlab="time")
par(new=T)
plot(1:6,apply(dat1[dat1$treat==1,paste0("dep",1:6)],2,mean,na.rm=T),type="l",lwd=2,ylim=c(0,6),ylab="",xlab="",col="red")
# individual slopes
for(i in 1:N){
  if(dat1$treat[i]==1){
    par(new=T)
    plot(1:6,jitter(as.numeric(dat1[i,paste0("dep",1:6)])),type="l",ylim=c(0,6),ylab="",xlab="",col=colo[C.est8[i]],lty=2)
  }
}

plot(1:6,apply(dat1[C.est8==2,paste0("dep",1:6)],2,mean,na.rm=T),type="l",lwd=2,ylim=c(0,6),ylab="dep",xlab="time")
par(new=T)
plot(1:6,apply(dat1[dat1$treat==0,paste0("dep",1:6)],2,mean,na.rm=T),type="l",lwd=2,ylim=c(0,6),ylab="",xlab="",col="red")
# individual slopes
for(i in 1:N){
  if(dat1$treat[i]==0){
    par(new=T)
    plot(1:6,jitter(as.numeric(dat1[i,paste0("dep",1:6)])),type="l",ylim=c(0,6),ylab="",xlab="",col=colo[C.est8[i]],lty=2)
  }
}






