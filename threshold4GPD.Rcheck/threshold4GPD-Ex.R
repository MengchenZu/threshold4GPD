pkgname <- "threshold4GPD"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('threshold4GPD')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("Data_Randomization")
### * Data_Randomization

flush(stderr()); flush(stdout())

### Name: Data_Randomization
### Title: Data_Randomization
### Aliases: Data_Randomization

### ** Examples

x <- rnorm(100)
x <- sort(x)
x <- round(x*10)/10
Data_Randomization(x)




cleanEx()
nameEx("Expl.diag")
### * Expl.diag

flush(stderr()); flush(stdout())

### Name: Expl.diag
### Title: producing the exponential diagnostic plots
### Aliases: Expl.diag

### ** Examples

library(mvtnorm)
xbvn<-rmvnorm(6000, sigma=matrix(c(1,0.7,0.7,1),2,2))

# Transform margins to exponential
xbvn.exp<- -log(1-pnorm(xbvn))
Expl.diag(apply(xbvn.exp,1,min), k=30, q1=0, param="Rate")
Expl.diag(apply(xbvn.exp,1,min), k=30, q1=0, param="InvRate")




cleanEx()
nameEx("NHPP.diag")
### * NHPP.diag

flush(stderr()); flush(stdout())

### Name: NHPP.diag
### Title: producing the Non-Homogeneous Poisson Process diagnostic plots
### Aliases: NHPP.diag

### ** Examples

## insert an easy example for test run
set.seed(1)
xnorm<-abs(rnorm(5000))
thresholds_xnorm <- Threshold_Generator_Uniform(xnorm, wetCuttingPoint=0.1)
ChangedThresholds_xnorm <- Threshold_Randomization(xnorm,thresholds_xnorm,seed = 1)
nhpp <-NHPP.diag(xnorm, u= ChangedThresholds_xnorm, plots=c("LRT", "WN", "PS"))

##insert an example for test run
set.seed(1)
x <- rnorm(100000)
x <- x[x > quantile(x, probs = 0.9)]
x <- sort(x)
x <- Data_Randomization(x,seed = 1)
x <- sort(x)
thresholds_x <- Threshold_Generator_Uniform(x, wetCuttingPoint=0.1)
ChangedThresholds_x <- Threshold_Randomization(x,thresholds_x,seed = 1)
nhpp <-NHPP.diag(x, u= ChangedThresholds_x,
                 M=365, nbs=1000, alpha=0.05, plots=c("LRT", "WN", "PS"), UseQuantiles=FALSE,
                 cex.lab=1.5, cex.axis=1.4,cex.main=2,mgp=c(4.2,1,0))

#########################
# View different plots for easy example:
# likelihood ratio plot only
NHPP.diag(xnorm, u= ChangedThresholds_xnorm, plots=c("LRT"))

# white noise plot only
NHPP.diag(xnorm, u= ChangedThresholds_xnorm, plots=c("WN"))

# parameter stability plot only
nhpp <-NHPP.diag(xnorm, u= ChangedThresholds_xnorm, plots=c("PS"))

# likelihood ratio test for assistance plot only
NHPP.diag(xnorm, u= ChangedThresholds_xnorm, plots=c("LRTEST"))

# Kolmogorov-Smirnov goodness-of-fit test for assistance plot only
NHPP.diag(xnorm, u= ChangedThresholds_xnorm, plots=c("KS"))

# Anderson-Darling test for assistance plot only
NHPP.diag(xnorm, u= ChangedThresholds_xnorm, plots=c("AD"))

# a showoff of all the plots
NHPP.diag(xnorm, u= ChangedThresholds_xnorm, plots=c("LRT", "WN", "PS", "LRTEST", "KS", "AD"))




cleanEx()
nameEx("Threshold_Generator_Density")
### * Threshold_Generator_Density

flush(stderr()); flush(stdout())

### Name: Threshold_Generator_Density
### Title: Thershold_Generator_Density
### Aliases: Threshold_Generator_Density

### ** Examples

set.seed(1)
xnorm<-abs(rnorm(5000))
Threshold_Generator_Density(xnorm)




cleanEx()
nameEx("Threshold_Generator_Uniform")
### * Threshold_Generator_Uniform

flush(stderr()); flush(stdout())

### Name: Threshold_Generator_Uniform
### Title: Thershold_Generator_Uniform
### Aliases: Threshold_Generator_Uniform

### ** Examples

set.seed(1)
xnorm<-abs(rnorm(5000))
Threshold_Generator_Uniform(xnorm)




cleanEx()
nameEx("Threshold_Randomization")
### * Threshold_Randomization

flush(stderr()); flush(stdout())

### Name: Threshold_Randomization
### Title: Threshold_Randomization
### Aliases: Threshold_Randomization

### ** Examples

set.seed(1)
xnorm<-abs(rnorm(5000))
thresholds_xnorm <- Threshold_Generator_Uniform(xnorm)
thresholds_xnorm
Threshold_Randomization(xnorm,thresholds_xnorm,seed = 1)




### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
