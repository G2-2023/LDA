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

# JAGS model definition
modelString <- "
  model {
    for(i in 1:N){
      for(j in 1:6){
        y[i,j] ~ dnorm(mu[i,j], tau)
        mu[i,j] <- alpha[class[i]] + beta[class[i]] * (j-1)
      }
      class[i] ~ dcat(pi[])
    }

    for(k in 1:2){
      alpha[k] ~ dnorm(0, .0001)
      beta[k] ~ dnorm(0, .0001)I(0,) # Censored normal distribution for class 2
    }
    
    tau ~ dgamma(0.01, 0.01)
    pi[1] <- 1 - p
    pi[2] <- p
    p ~ dbeta(1,1)
  }
"

# Initial values
initValues <- function(){
  list(alpha=c(0,0), beta=c(0,0), class=rep(1, N), p=0.5) # class is now a numeric vector
}

# Data list for JAGS
dataList <- list(
  y = as.matrix(dat1[,1:6]),
  N = N
)

# Parameters to monitor
params <- c("alpha", "beta", "tau", "class", "p")

# Run the model
library(rjags)
jagsModel <- jags.model(textConnection(modelString), data = dataList, inits = initValues, n.chains = 3)
update(jagsModel, 1000) # Burn-in
samples <- coda.samples(jagsModel, variable.names = params, n.iter = 5000)

# Check results
print(samples)
