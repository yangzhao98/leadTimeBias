## Required information
# 1) n.scr   = number of participants at screening
# 2) can.scr = number of cancers detected at screening
# 3) py      = number of person-years contributed in the yth year after screening
# 4) J       = incidence of cancers for unscreened population
# 5) o       = number of interval cancers per month
# 6) time    = time since screening
# 7) age     = age at screening test
# 8) N       = number of time intervals

# sink("model.txt")
model <- "
  model {
    can.scr ~ dbin(p.scr, n.scr)
    p.scr <- (S*J*(exp(-lambda*age)-exp(-J*age))/(J-lambda)/(exp(-J*age)+(J*(exp(-lambda*age)-exp(J*age))/(J-lambda))))
    for (i in 1:N) {
      o[i] ~ dpois(theta.obs[i])
      theta.new[i] <- py[i]*J*(1-exp(-lambda*(time[i]-0.5)))
      theta.und[i] <- can.scr*exp(logit.S)*(exp(-lambda*(time[i]-1))-exp(-lambda*time[i]))
      theta.obs[i] <- theta.new[i] + theta.und[i]
    }
    logit.S ~ dunif(-4,4)
    S <- 1/(1+exp(logit.S))
    lambda ~ dgamma(0.001,0.001)I(0.05,5)
    mst <- 1/lambda
  }
"
# sink()

## Total

inits <- function() { list(logit.S = -2, lambda = 0.5)}
params <- c("S", "lambda", "mst", "p.scr")

# S <- 0.85
# logit.S <- log(1-S)/S
# logit.S

library(R2OpenBUGS)
# library(R2WinBUGS)
library(MASS)
library(coda)
library(tidyverse)
library(runjags)

## Total
data <- list(n.scr = 71307,
             can.scr = 131,
             py = c(67824, 62667, 39853),
             J = 1730/2078934,
             o = c(44, 29, 22),
             N = 3,
             time = c(1,2,3),
             age = 60)

## Using runjags package for Bayesian analysis
leadtime.sim.jags <- run.jags(data = data, inits = inits,
                              model = model,
                              monitor = params,
                              n.chains = 3,
                              sample = 1e4,
                              burnin = 5e3,
                              # silent.jags = TRUE,
                              adapt = 1e3)

out <- summary(leadtime.sim.jags)[,c("Mean", "SD", "Lower95", "Median", "Upper95", "psrf")]
pdf(file="FigS3_ConvergencePlot.pdf",width=10,height=6,onefile=FALSE)
plot(leadtime.sim.jags, var = c("S", "lambda", "mst"),plot.type=c("trace"))
dev.off()
knitr::kable(out, caption = "Estimates based posterior distribution")
gelman.plot(leadtime.sim.jags)

leadtime.result <- data.frame(as.mcmc(leadtime.sim.jags))
head(leadtime.result)

leadtime <- function(OS, Status, lambda) {
  if (Status == 1) {
    s <- (1 - exp(-lambda*OS) - lambda*OS*exp(-lambda*OS))/(lambda*(1 - exp(-lambda*OS)))
  } else {
    s <- (1 - exp(-lambda*OS))/lambda
  }
  return(s)
}