##automate the threshold selection based on non-homogeneous point process
##Warren Jin, Created in Jan 2016, refined on 12 Feb 2016

#######################################################################################################
# These following functions to create plots and tests in article
#
# "Exploiting structure of maximum likelihood estimators for extreme value threshold selection"
#  Jennifer Wadsworth from University of Cambridge
#  Refined by Warren Jin from CSIRO
#######################################################################################################
#
#



#############################################################################################################
# NHPP.diag
#############################################################################################################
#
# Details:
#
# Function to produce diagnostic plots and test statistics for the NHPP model
#
# Arguments:
#
# x - data vector
# u - optional; vector of candidate thresholds
# k - number of thresholds to consider (if u unspecified)
# q1 -  lowest quantile for the threshold sequence
# q2 -  upper quantile limit for the threshold sequence (q2 itself is not used as a threshold, but rather the uppermost
#       threshold will be at the (q2-1/k) quantile)
# par - starting values for the optimization
# M - number of superpositions or "blocks" / "years" the process corresponds to (can affect the optimization)
# nbs - number of simulations used to assess the null distribution of the LRT, and produce the p-value
# alpha - significance level of the LRT
# plots - which plots to produce; "LRT"= likelihood ratio test, "WN" = white noise, "PS" = parameter stability
# UseQuantiles - logical; use quantiles as the thresholds in the plot?
# pmar - vector of length 4 giving the arguments for the plot margins in par(mar=c(*,*,*,*))
# head - number of thresholds inserted between the first threshold and the second threshold temporary to reduce the
#       variance at the begining of likelihood ratio statistic test
# MLE_model - which model of maximum likelihood estimation is used;
#       "median"= use the median data of 7 maximum likelihood estimation (not accurate enough but suitable for most data vector)
#       "min"= use the min data of 7 maximum likelihood estimation (most accurate but not suitable for some small data vector)
#       "second_min"= use the second min data of 7 maximum likelihood estimation (almost same as "min")
# use_Data_Randomization - whether or not use Data_Randomization to randomize the given data
# use_Threshold_Randomization - whether or not use Threshold_Randomization to randomize the thresholds
# DEBUGing - default as false, set to true if you want to debug
# ... - other parameter
#
# Value: (results:)
#
# Plots the requested diagnostics. Also returns a list with components:
#
# MLEall - MLEs from all thresholds
# Cov.xi - Joint asymptotic covariance matrix for xi
# WN - values of the white noise process
# LRT - values of the LT test statistic vs threshold
# pval - p-value of the LR test
# k - final number of thresholds used
# thresh - threshold selected by LR procedure
# mle.u - MLE from selected threshold

#' @title producing the Non-Homogeneous Poisson Process diagnostic plots
#'
#' @description Function to produce diagnostic plots and test statistics for the NHPP model
#' @details Function to produce diagnostic plots and test statistics for the NHPP model
#' @param x - data vector
#' @param u - optional; vector of candidate thresholds
#' @param k - number of thresholds to consider (if u unspecified)
#' @param q1 -  lowest quantile for the threshold sequence
#' @param q2 -  upper quantile limit for the threshold sequence (q2 itself is not used as a threshold, but rather the uppermost
#'       threshold will be at the (q2-1/k) quantile)
#' @param par - starting values for the optimization
#' @param M - number of superpositions or "blocks" / "years" the process corresponds to (can affect the optimization)
#' @param nbs - number of simulations used to assess the null distribution of the LRT, and produce the p-value
#' @param alpha - significance level of the LRT
#' @param plots - which plots to produce;
#'       "LRT"= likelihood ratio
#'       "WN" = white noise
#'       "PS" = parameter stability
#'       "LRTEST" = likelihood ratio test for assistance
#'       "KS" = Kolmogorov-Smirnov goodness-of-fit test for assistance
#'       "AD" = Anderson-Darling test for assistance
#' @param UseQuantiles - logical; use quantiles as the thresholds in the plot?
#' @param pmar - vector of length 4 giving the arguments for the plot margins in par(mar=c(*,*,*,*))
#' @param head - number of thresholds inserted between the first threshold and the second threshold temporary to reduce the
#'       variance at the begining of likelihood ratio statistic test
#' @param MLE_model - which model of maximum likelihood estimation is used;
#'       "median"= use the median data of 7 maximum likelihood estimation (not accurate enough but suitable for most data vector)
#'       "min"= use the min data of 7 maximum likelihood estimation (most accurate but not suitable for some small data vector)
#'       "second_min"= use the second min data of 7 maximum likelihood estimation (almost same as "min")
#' @param use_Data_Randomization - whether or not use Data_Randomization to randomize the given data
#' @param use_Threshold_Randomization - whether or not use Threshold_Randomization to randomize the thresholds
#' @param DEBUGing - default as false, set to true if you want to debug
#' @param ... - other parameter
#'
#' @import ismev
#' @import mgcv
#' @import ADGofTest
#' @import dgof
#'
#' @return MLEall - MLEs from all thresholds
#' @return Cov.xi - Joint asymptotic covariance matrix for xi
#' @return WN - values of the white noise process
#' @return LRT - values of the LT test statistic vs threshold
#' @return pval - p-value of the LR test
#' @return k - final number of thresholds used
#' @return thresh - threshold selected by LR procedure
#' @return mle.u - MLE from selected threshold
#'
#' @export
#'
#' @examples
#' ## insert an easy example for test run
#' set.seed(1)
#' xnorm<-abs(rnorm(5000))
#' thresholds_xnorm <- Threshold_Generator_Uniform(xnorm, wetCuttingPoint=0.1)
#' ChangedThresholds_xnorm <- Threshold_Randomization(xnorm,thresholds_xnorm,seed = 1)
#' nhpp <-NHPP.diag(xnorm, u= ChangedThresholds_xnorm, plots=c("LRT", "WN", "PS"))
#'
#' ##insert an example for test run
#' set.seed(1)
#' x <- rnorm(100000)
#' x <- x[x > quantile(x, probs = 0.9)]
#' x <- sort(x)
#' x <- Data_Randomization(x,seed = 1)
#' x <- sort(x)
#' thresholds_x <- Threshold_Generator_Uniform(x, wetCuttingPoint=0.1)
#' ChangedThresholds_x <- Threshold_Randomization(x,thresholds_x,seed = 1)
#' nhpp <-NHPP.diag(x, u= ChangedThresholds_x,
#'                  M=365, nbs=1000, alpha=0.05, plots=c("LRT", "WN", "PS"), UseQuantiles=FALSE,
#'                  cex.lab=1.5, cex.axis=1.4,cex.main=2,mgp=c(4.2,1,0))
#'
#' #########################
#' # View different plots for easy example:
#' # likelihood ratio plot only
#' NHPP.diag(xnorm, u= ChangedThresholds_xnorm, plots=c("LRT"))
#'
#' # white noise plot only
#' NHPP.diag(xnorm, u= ChangedThresholds_xnorm, plots=c("WN"))
#'
#' # parameter stability plot only
#' nhpp <-NHPP.diag(xnorm, u= ChangedThresholds_xnorm, plots=c("PS"))
#'
#' # likelihood ratio test for assistance plot only
#' NHPP.diag(xnorm, u= ChangedThresholds_xnorm, plots=c("LRTEST"))
#'
#' # Kolmogorov-Smirnov goodness-of-fit test for assistance plot only
#' NHPP.diag(xnorm, u= ChangedThresholds_xnorm, plots=c("KS"))
#'
#' # Anderson-Darling test for assistance plot only
#' NHPP.diag(xnorm, u= ChangedThresholds_xnorm, plots=c("AD"))
#'
#' # a showoff of all the plots
#' NHPP.diag(xnorm, u= ChangedThresholds_xnorm, plots=c("LRT", "WN", "PS", "LRTEST", "KS", "AD"))
#'

NHPP.diag<-function(x, u= NULL, k, q1=NULL, q2=1, par=NULL, M=NULL, nbs=1000, alpha=0.05,
                    plots=c("LRT", "WN", "PS", "LRTEST", "KS", "AD"), UseQuantiles=TRUE,
                    DEBUGing=FALSE,pmar=c(5.5,7,3,3),head=5,MLE_model = NULL,
                    use_Data_Randomization=FALSE, use_Threshold_Randomization=FALSE,...)
{
  if(use_Data_Randomization) x <- Data_Randomization(x)

  if(is.null(u)) u <- Threshold_Generator_Density(x)

  if(use_Threshold_Randomization) u <- Threshold_Randomization(x,u)

  if(is.null(MLE_model)) {
    if(length(u) > 40) MLE_model = "min"
    else if(length(u) > 30) MLE_model = "second_min"
    else MLE_model = "median"
  }

  u1 <- u[1]
  u2 <- u[2]
  add_length = 0
  # Added by Mengchen, for "robust head LR statistic"
  if(head >0) {
    temp <- x[x > u[1] & x < u[2]]
    temp <- unique(temp)
    u1 <- u[1]
    u2 <- u[2]
    threshold_start <- temp[1]
    threshold_end <- temp[length(temp)]
    if(length(temp) > 400) {
      add_thresholds <- seq(threshold_start , threshold_end , len = head+2)
      add_thresholds <- add_thresholds[c(-1,-length(add_thresholds))]
    } else if(length(temp)>100) {
      add_thresholds <- temp[temp = seq(1,length(temp),by=100)]
      add_thresholds <- add_thresholds[c(-1,-length(add_thresholds))]
    } else {
      add_thresholds <- NULL
    }
    u<- sort(c(u,add_thresholds))
  }
  # Ended of Edit

  ##Warren Jin inserts on 25 Jan 2016 the following code to sort u from low to high
  u <- unique(sort(u))
  k<-length(u);

  ## Added by Mengchen on 13 July 2016, for "find out lowest quantile for the threshold sequence"
  if(is.null(q1))
  {
    de <- density(x)
    q1 <- de$x[which.max(de$y)]
    q1 <- ecdf(x)(q1)
  }
  ## End of Added

  if(is.null(M))
  {
    M<-length(x[x>min(u)])/3
  }

  if(is.null(par))
  {
    require(ismev)
    ppf<-pp.fit(xdat=x,thresh=min(u),npy=length(x)/M,show=F)
    par<-ppf$mle
  }

  par(mfrow=c(length(plots), 1))
  par(mar=pmar)

  par(las=1)

  J1<-Joint.MLE.NHPP(x=x, u=u, k=k, q1=q1, q2=q2, par=par, M=M, MLE_model = MLE_model)
  warn<-any(eigen(J1$Cov.xi)$val<=.Machine$double.xmin*10^10)
  ##Warren Jin inserts on 25 Jan 2016 the following code to avoid duplicated thresholds and save time, especially for large k
  uOriginal=uOriginalIn=u
  thresholdsLenMax<- round(length(uOriginalIn)*3)
  thresholdsLenMin <- round(length(uOriginalIn)*0.5)
  while(warn){
    if(thresholdsLenMax >= length(uOriginalIn)*2+1  || length(u) >= thresholdsLenMin ){
      ###to remove a threshold u element which has shortest distance with its neighbors
      #to remove threshold elements which have lowest difference among the two consecutive diagonal element of J1$Cov.xi
      ind2Rm <- NULL
      diag4cxi<- diag(J1$Cov.xi)
      #Diagonal elements are assumed (or should) increase slightly
      while((any(diff(c(0,diag4cxi))< .Machine$double.xmin*10^10))){
        #prod(diff(c(0,diag4cxi))) < 0){
        if(DEBUGing && sum(diff(c(0,diag4cxi))<0) >= k/2) browser();
        ##choose the one with the maximum downward change to remove
        (newInd<- which.min(diff(c(0,diag4cxi))))
        (actualNewInd <- newInd + sum(ind2Rm<=newInd))
        ind2Rm <- c(ind2Rm,actualNewInd)
        diag4cxi <- diag4cxi[-(newInd)]
        u <- u[-(newInd)]
        k <- k-1;
      }
      if(DEBUGing && length(ind2Rm)<1) browser();
    }else{
      #regenerate new thresholds list
      k<- thresholdsLenMax <- thresholdsLenMax+1
      print(paste(' thresholdsLenMax becomes',thresholdsLenMax))
      u<- quantile(x,seq(q1+(q2-q1)/3,q2,length.out=k))
      u <- as.numeric(u[-length(u)])
      u <- sort(unique(c(u,uOriginalIn)))
      k=length(u)
      uOriginal<- u;
    }

    J1<-Joint.MLE.NHPP(x=x, u=u, q1=q1, q2=q2, par=par,M=M, MLE_model = MLE_model)
    warn<-any(eigen(J1$Cov.xi)$val <=.Machine$double.xmin*10^10)
  }

  #Add by Mengchen, for showing deleted thresholds, using "robust head LR statistic"
  u_deleted <- uOriginal[!uOriginal %in% u]
  remain_length <- length(u[u<=u1|u>=u2])
  u_deleted <- u_deleted[u_deleted<=u1|u_deleted>=u2]
  if(length(u_deleted)!=0) print(sprintf(' Note: %d thresholds used after %d thresholds (%s) have been removed',
                                         remain_length,length(u_deleted),
                                         paste(u_deleted,collapse=', ')))
  #End of Added

  while(any(eigen(J1$Cov.xi)$val<=.Machine$double.xmin*10^10))
  {
    k<-k-1
    ##Warren Jin inserts the code to avoid duplicated thresholds and save time, especially for large k
    thresholdsTemp<- quantile(x,seq(q1,q2,len=k+1))
    if(length(unique( thresholdsTemp)) <= k) next;
    ##End of insertion

    J1<-Joint.MLE.NHPP(x=x,k=k, q1=q1, q2=q2, par=par,M=M, MLE_model = MLE_model)
  }
  if(warn){warning(paste("Estimated covariance matrix for xi not positive definite for initial k. Final k:", k))}

  if(length(u[u<=u1|u>=u2])<= 3) {
    print("too less thresholds.")
    return(NULL)
  } else {
    ##normalised xi star
    wn<-C1(k)%*%J1$mle[,3]/sqrt(diag(C1(k)%*%J1$Cov.xi%*%t(C1(k))))
    nl<-norm.LRT(x=wn,u=u[-c(1,k+1)])

    # Added by Mengchen, for "robust head LR statistic"
    if(head > 0) {
      add_length <- length(u[u>u1&u<u2])
      if(add_length != 0) {
        u <- u[-c(seq(2,add_length+1,by=1))]
        wn <- wn[-c(seq(2,add_length+1,by=1))]
        nl <- nl[-c(seq(2,add_length+1,by=1)-1),]
      }
    }
    # Ended of Edit

    p.value.list <- NULL
    for(i in 1:5) {
      nlt<-NULL
      for(j in 1:nbs)
      {
        nlt[j]<-max(norm.LRT(x=rnorm(k-1),u[-c(1,k+1)])[,2])
      }

      p.value.list<-c(p.value.list,length(nlt[nlt>max(nl[,2])])/nbs)
    }
    pval <- median(p.value.list)

    if(pval<alpha){ustar<-nl[nl[,2]==max(nl[,2]),1]
    } else{ustar<-min(u)}
    ustar4print<-as.numeric(sprintf('%.3f',ustar))
    ind<-u[-(k+1)]==ustar
    theta.hat<-J1$mle[ind,]

    # Added by Mengchen on 17 Jul 2016
    plotslength <- length(plots)
    if(plotslength == 0 | plotslength == 1) { par(mfrow=c(1,1))
    } else if(plotslength == 2) { par(mfrow=c(2,1))
    } else if(plotslength == 3) { par(mfrow=c(3,1))
    } else if(plotslength == 4) { par(mfrow=c(2,2))
    } else {par(mfrow=c(3,2))}
    # End of Added

    if(is.element("LRT",plots)){
      plot(u[-c(k+1)], c(rep(NA,2),nl[,2]),xlab="Threshold",ylab="LR statistic", main=paste("p-value:", pval))
    }

    if(is.element("WN",plots)){
      plot(u[-c(k+1)], c(NA,wn),xlab="Threshold",ylab="White noise")
      abline(h=0,col=2)
      abline(v=ustar,col=4)
      text(ustar,mean(c(NA,wn),na.rm=TRUE),bquote(u==.(ustar4print)),pos=4,col='blue')
    }

    if(is.element("PS",plots)){
      if(add_length != 0) {
        TradCI<-cbind(J1$mle[,3]-qnorm(0.975)*sqrt(diag(J1$Cov.xi)),J1$mle[,3]+qnorm(0.975)*sqrt(diag(J1$Cov.xi)))
        plot(u[-(k+1)],J1$mle[,3][-c(seq(2,add_length+1,by=1))],ylim=c(min(TradCI[,1]),max(TradCI[,2])), xlab="Threshold", ylab=expression(hat(xi)))
        lines(u[-(k+1)],TradCI[,1][-c(seq(2,add_length+1,by=1))], lty=2)
        lines(u[-(k+1)],TradCI[,2][-c(seq(2,add_length+1,by=1))], lty=2)
        abline(v=ustar,col=4)
        text(ustar,mean(J1$mle[,3][-c(seq(2,add_length+1,by=1))]),bquote(u==.(ustar4print)),pos=4,col='blue')
      } else {
        TradCI<-cbind(J1$mle[,3]-qnorm(0.975)*sqrt(diag(J1$Cov.xi)),J1$mle[,3]+qnorm(0.975)*sqrt(diag(J1$Cov.xi)))
        plot(u[-(k+1)],J1$mle[,3],ylim=c(min(TradCI[,1]),max(TradCI[,2])), xlab="Threshold", ylab=expression(hat(xi)))
        lines(u[-(k+1)],TradCI[,1], lty=2)
        lines(u[-(k+1)],TradCI[,2], lty=2)
        abline(v=ustar,col=4)
        text(ustar,mean(J1$mle[,3]),bquote(u==.(ustar4print)),pos=4,col='blue')
      }
    }

    # Added by Mengchen on 15 Jul 2016
    if(is.element("LRTEST",plots)) {
      pvallist <- NULL
      for(f in 1:length(nl[,2]))
      {
        pvallist[f] <- length(nlt[nlt>nl[f,2]])/nbs
      }
      LR_Forward <- K.ForwardStop(pvallist)
      LR_Strong <- K.StrongStop(pvallist)
      #par(mfrow=c(1,1))
      plot(u,c(NA,NA,pvallist),xlab="Threshold",ylab="LR TEST",type="b",col="blue")
      abline(h=0.05,col="red")

      #plot(u,c(NA,NA,LR_Forward$forward.list),xlab="Threshold",ylab="LR TEST",type="b",col="blue")
      #abline(h=0.05,col="red")
      #plot(u,c(NA,NA,LR_Strong$strong.list),xlab="Threshold",ylab="LR TEST",type="b",col="blue")
      #abline(h=0.05,col="red")
    }

    if(is.element("KS",plots)) {
      require(dgof)
      kstest <- 1
      for(i in 1:length(wn)) {
        kstest <- c(kstest,ks.test(wn[i:length(wn)],pnorm,0,1)$p.value)
      }
      kstest <- kstest[-1]
      KS_Forward <- K.ForwardStop(kstest)
      KS_Strong <- K.StrongStop(kstest)
      #par(mfrow=c(1,1))
      plot(u,c(NA,kstest),xlab="Threshold",ylab="KS TEST",type="b",col="blue")
      abline(h=0.05,col="red")

      #plot(u,c(NA,KS_Forward$forward.list),xlab="Threshold",ylab="LR TEST",type="b",col="blue")
      #abline(h=0.05,col="red")
      #plot(u,c(NA,KS_Strong$strong.list),xlab="Threshold",ylab="LR TEST",type="b",col="blue")
      #abline(h=0.05,col="red")
    }

    if(is.element("AD",plots)) {
      require(ADGofTest)
      adtest <- 1
      for(i in 1:length(wn))
      {
        adtest <- c(adtest,ad.test(wn[i:length(wn)],pnorm,0,1)$p.value)
      }
      adtest <- adtest[-1]
      AD_Forward <- K.ForwardStop(adtest)
      AD_Strong <- K.StrongStop(adtest)
      #par(mfrow=c(1,1))
      plot(u,c(NA,adtest),xlab="Threshold",ylab="AD TEST",type="b",col="blue")
      abline(h=0.05,col="red")

      #plot(u,c(NA,AD_Forward$forward.list),xlab="Threshold",ylab="LR TEST",type="b",col="blue")
      #abline(h=0.05,col="red")
      #plot(u,c(NA,AD_Strong$strong.list),xlab="Threshold",ylab="LR TEST",type="b",col="blue")
      #abline(h=0.05,col="red")
    }
    # End of Added

    return(list(MLEall=J1$mle, Cov.xi=J1$Cov.xi, WN=wn, LRT = nl, pval = pval, k=k, thresh=ustar, mle.u=theta.hat))
  }
}



#############################################################################################################
# Expl.diag
#############################################################################################################
#
# Details:
#
# Function to produce diagnostic plots and test statistics for the rate / inverse rate parameter of the Exponential model
#
# Arguments:
#
# x - data vector
# u - optional; vector of candidate thresholds
# k - number of thresholds to consider (if u unspecified)
# q1 -  lowest quantile for the threshold sequence
# q2 -  upper quantile limit for the threshold sequence (q2 itself is not used as a threshold, but rather the uppermost
#       threshold will be at the (q2-1/k) quantile)
# nbs - number of simulations used to assess the null distribution of the LRT, and produce the p-value
# alpha - significance level of the LRT
# plots - which plots to produce; "LRT"= likelihood ratio test, "WN" = white noise, "PS" = parameter stability
# UseQuantiles - logical; use quantiles as the thresholds in the plot?
# param - character specifying "InvRate" or "Rate" for either inverse rate parameter / rate parameter, respectively
# pmar - vector of length 4 giving the arguments for the plot margins in par(mar=c(*,*,*,*))
# ... - other parameter
#
# Value: (results:)
#
# Plots the requested diagnostics. Also returns a list with components:
#
# MLE - MLEs from all thresholds
# Cov - Joint asymptotic covariance matrix for xi
# WN - values of the white noise process
# LRT - values of the LT test statistic vs threshold
# pval - p-value of the LR test
# k - final number of thresholds used
# thresh - threshold selected by the LR procedure
# mle.u - MLE from selected threshold

#' @title producing the exponential diagnostic plots
#' @description Function to produce diagnostic plots and test statistics for the rate / inverse rate parameter of the Exponential model
#'
#' @param x - data vector
#' @param u - optional; vector of candidate thresholds
#' @param k - number of thresholds to consider (if u unspecified)
#' @param q1 -  lowest quantile for the threshold sequence
#' @param q2 -  upper quantile limit for the threshold sequence (q2 itself is not used as a threshold, but rather the uppermost
#       threshold will be at the (q2-1/k) quantile)
#' @param nbs - number of simulations used to assess the null distribution of the LRT, and produce the p-value
#' @param alpha - significance level of the LRT
#' @param plots - which plots to produce; "LRT"= likelihood ratio test, "WN" = white noise, "PS" = parameter stability
#' @param UseQuantiles - logical; use quantiles as the thresholds in the plot?
#' @param param - character specifying "InvRate" or "Rate" for either inverse rate parameter / rate parameter, respectively
#' @param pmar - vector of length 4 giving the arguments for the plot margins in par(mar=c(*,*,*,*))
#' @param ... - other parameter
#'
#' @return MLE - MLEs from all thresholds
#' @return Cov - Joint asymptotic covariance matrix for xi
#' @return WN - values of the white noise process
#' @return LRT - values of the LT test statistic vs threshold
#' @return pval - p-value of the LR test
#' @return k - final number of thresholds used
#' @return thresh - threshold selected by the LR procedure
#' @return mle.u - MLE from selected threshold
#'
#' @export
#'
#' @examples
#' library(mvtnorm)
#' xbvn<-rmvnorm(6000, sigma=matrix(c(1,0.7,0.7,1),2,2))
#'
#' # Transform margins to exponential
#' xbvn.exp<- -log(1-pnorm(xbvn))
#' Expl.diag(apply(xbvn.exp,1,min), k=30, q1=0, param="Rate")
#' Expl.diag(apply(xbvn.exp,1,min), k=30, q1=0, param="InvRate")
#'

Expl.diag<-function(x, u= NULL, k, q1, q2=1, nbs=1000, alpha=0.05, plots=c("LRT", "WN", "PS"), UseQuantiles=TRUE, param="InvRate",pmar=c(5.5,7,3,3),...)
{
  unull<-is.null(u)
  if(!unull){k<-length(u)}

  par(mfrow=c(length(plots), 1))
  par(mar=pmar)

  par(las=1)

  J1<-Joint.MLE.Expl(x=x,u=u,k=k, q1=q1, q2=q2, param=param)
  warn<-any(eigen(J1$Cov)$val<=.Machine$double.xmin*10^10)
  if(!unull&&warn){stop("Estimated covariance matrix for eta not positive definite: try different thresholds")}

  while(any(eigen(J1$Cov)$val<=.Machine$double.xmin*10^10))
  {
    k<-k-1
    J1<-Joint.MLE.Expl(x=x,k=k, q1=q1, q2=q2, param=param)
  }
  if(warn){warning(paste("Estimated covariance matrix for 1/eta not positive definite for initial k. Final k:", k))}

  if(unull){u<-quantile(x,seq(q1,1,len=k+1))}

  wn<-C1(k)%*%J1$mle/sqrt(diag(C1(k)%*%J1$Cov%*%t(C1(k))))
  nl<-norm.LRT(x=wn,u=u[-c(1,k+1)])

  nlt<-NULL
  for(j in 1:nbs)
  {
    nlt[j]<-max(norm.LRT(x=rnorm(k-1),u[-c(1,k+1)])[,2])
  }

  pval<-length(nlt[nlt>max(nl[,2])])/nbs

  if(pval<alpha){ustar<-nl[nl[,2]==max(nl[,2]),1]}
  else{ustar<-min(u)}
  ind<-u[-(k+1)]==ustar
  theta.hat<-J1$mle[ind]

  if(unull){qs<- seq(q1, q2, len=k+1)[-(k+1)]}

  if(is.element("LRT",plots)){
    if(unull&&UseQuantiles)
    {plot(qs, c(NA,NA,nl[,2]),xlab="Quantile",ylab="LR statistic", main=paste("p-value:", pval),...)}
    else{plot(u[-c(k+1)], c(NA,NA,nl[,2]),xlab="Threshold",ylab="LR statistic", main=paste("p-value:", pval),...)}
  }

  if(is.element("WN",plots)){
    if(unull&&UseQuantiles){plot(qs, c(NA,wn),xlab="Quantile",ylab="White noise",...)
      abline(h=0,col=2)
      abline(v=mean(x<=ustar),col=4)
    }
    else{plot(u[-c(k+1)], c(NA,wn),xlab="Threshold",ylab="White noise",...)
      abline(h=0,col=2)
      abline(v=ustar,col=4)}
  }

  if(is.element("PS",plots)){
    TradCI<-cbind(J1$mle-qnorm(0.975)*sqrt(diag(J1$Cov)),J1$mle+qnorm(0.975)*sqrt(diag(J1$Cov)))
    if(UseQuantiles)
    {
      if(param=="InvRate"){plot(qs,J1$mle,ylim=c(min(TradCI[,1]),max(TradCI[,2])), xlab="Quantile", ylab=expression(hat(eta)),...)}
      else if(param=="Rate"){plot(qs,J1$mle,ylim=c(min(TradCI[,1]),max(TradCI[,2])), xlab="Quantile", ylab=expression(hat(theta)),...)}
      lines(qs,TradCI[,1], lty=2)
      lines(qs,TradCI[,2], lty=2)
      abline(v=mean(x<=ustar),col=4)
    }
    else{
      if(param=="InvRate"){plot(u[-(k+1)],J1$mle,ylim=c(min(TradCI[,1]),max(TradCI[,2])), xlab="Threshold", ylab=expression(hat(eta)),...)}
      else if(param=="InvRate"){plot(u[-(k+1)],J1$mle,ylim=c(min(TradCI[,1]),max(TradCI[,2])), xlab="Threshold", ylab=expression(hat(theta)),...)}
      lines(u[-(k+1)],TradCI[,1], lty=2)
      lines(u[-(k+1)],TradCI[,2], lty=2)
      abline(v=ustar,col=4)
    }
  }

  return(list(MLE=J1$mle, Cov=J1$Cov, WN=wn, LRT = nl, pval = pval, k=k, thresh=ustar, mle.u=theta.hat))
}



#####################################################################################
# nhpp.nll
########################

# Details:

# Negative log likelihood for the NHPP model, to minimize for MLEs

# Arguments:

# theta - parameter vector (mu, sigma, xi)
# x - data vector
# u - threshold
# M - number of superpositions or "blocks" / "years" the process corresponds to (affects estimation of mu, sigma,
# but these can be changed post-hoc to correspond to any number)


# Value: (results:)

# negative log-likelihood value

#' @title Non-Homogeneous Poisson Process negative log-likelihood
#' @description Negative log likelihood for the NHPP model, to minimize for MLEs
#'
#' @param theta - parameter vector (mu, sigma, xi)
#' @param x - data vector
#' @param u - threshold
#' @param M - number of superpositions or "blocks" / "years" the process corresponds to (affects estimation of mu, sigma,
#'  but these can be changed post-hoc to correspond to any number)
#'
#' @return negative log-likelihood value
#'

nhpp.nll<-function(theta,x,u,M)
{
  x<-x[x>u]
  mu<-theta[1]; sig<-theta[2]; xi<-theta[3]
  Nu<-length(x)

  if(sig<=0|| any(1+xi*(x-mu)/sig < 0)){return(10e10)}

  else{
    if(abs(xi)>1e-10)
    {
      nll<- (1/xi+1)*sum(log(1+xi*(x-mu)/sig)) + Nu*log(sig) + M*(1+xi*(u-mu)/sig)^(-1/xi)
    }
    else{
      nll<-Nu*log(sig)+sum((x-mu)/sig) + M*exp(-(u-mu)/sig)
    }
    return(nll)
  }
}



#####################################################################################
# Joint.MLE.NHPP
#####################################################################################
#
# Details:
#
# Calculates the MLEs of the parameters (mu,sigma,xi), and joint asymptotic covariance matrix of these MLEs
# over a range of thresholds as supplied by the user.
#
# Arguments:
#
# x - vector of data
# u - vector of thresholds. If not supplied, then k thresholds between quantiles (q1, q2) will be used
# k - number of thresholds to consider if u not supplied
# q1 -  lowest quantile for the threshold sequence
# q2 -  upper quantile limit for the threshold sequence (q2 itself is not used as a threshold, but rather the uppermost
#       threshold will be at the (q2-1/k) quantile)
# par - starting values for the optimization
# M - number of superpositions or "blocks" / "years" the process corresponds to (affects estimation of mu, sigma,
# but these can be changed post-hoc to correspond to any number)
# MLE_model - which model of maximum likelihood estimation is used;
#       "median"= use the median data of 7 maximum likelihood estimation (not accurate enough but suitable for most data vector)
#       "min"= use the min data of 7 maximum likelihood estimation (most accurate but not suitable for some small data vector)
#       "second_min"= use the second min data of 7 maximum likelihood estimation (almost same as "min")
#
# Value: (returns:)
#
# mle - matrix of MLEs above the supplied thresholds; columns are (mu, sigma, xi)
# Cov.all - joint asymptotic covariance matrix of all MLEs
# Cov.mu - joint asymptotic covariance matrix of MLEs for mu
# Cov.sig - joint asymptotic covariance matrix of MLEs for sig
# Cov.xi - joint asymptotic covariance matrix of MLEs for xi

#' @title joint maximum likelihood estimation Non-Homogeneous Poisson Process
#' @description the MLEs of the parameters (mu,sigma,xi), and joint asymptotic covariance matrix of these MLEs
#' over a range of thresholds as supplied by the user.
#'
#' @param x - vector of data
#' @param u - vector of thresholds. If not supplied, then k thresholds between quantiles (q1, q2) will be used
#' @param k - number of thresholds to consider if u not supplied
#' @param q1 - lowest quantile for the threshold sequence
#' @param q2 -  upper quantile limit for the threshold sequence (q2 itself is not used as a threshold, but rather the uppermost
#'       threshold will be at the (q2-1/k) quantile)
#' @param par - starting values for the optimization
#' @param M - number of superpositions or "blocks" / "years" the process corresponds to (affects estimation of mu, sigma,
#'  but these can be changed post-hoc to correspond to any number)
#' @param MLE_model - which model of maximum likelihood estimation is used;
#'       "median"= use the median data of 7 maximum likelihood estimation (not accurate enough but suitable for most data vector)
#'       "min"= use the min data of 7 maximum likelihood estimation (most accurate but not suitable for some small data vector)
#'       "second_min"= use the second min data of 7 maximum likelihood estimation (almost same as "min")
#'
#' @return mle - matrix of MLEs above the supplied thresholds; columns are (mu, sigma, xi)
#' @return Cov.all - joint asymptotic covariance matrix of all MLEs
#' @return Cov.mu - joint asymptotic covariance matrix of MLEs for mu
#' @return Cov.sig - joint asymptotic covariance matrix of MLEs for sig
#' @return Cov.xi - joint asymptotic covariance matrix of MLEs for xi
#'

Joint.MLE.NHPP<-function(x, u=NULL, k, q1, q2=1, par, M, MLE_model = "median")
{
  x <- sort(x)
  if(!is.null(u))
  {
    k<-length(u)
    x<-x[x>u[1]]
    # add threshold above all data, to "close" final region
    u<-c(u,max(x)+1)
  }
  else{
    u<-quantile(x,seq(q1,q2,len=k+1))
  }

  I<-Iinv<-list()
  thetahat<-matrix(NA,ncol=3,nrow=k)

  for(i in 1:k)
  {
    if(u[i]> x[5]) {
      u_list1 <- u[i]
      u_list2 <- u[i]
      test_number <- 3
      for(j in 1:test_number) {
        temp <- max(x[x<u_list1[j]])
        u_list1 <- c(u_list1,temp)
      }
      for(j in 1:test_number) {
        temp <- min(x[x>u_list2[j]])
        u_list2 <- c(u_list2,temp)
      }
      u_list <- sort(unique(c(u_list1,u_list2)))
      opt <- matrix(0,1,3)
      mu <- 1
      sig <- 1
      xi <- 1
      if(i == 1) {upper_number = test_number + 1
      } else upper_number = 2*test_number+1
      for(j in 1:upper_number) {
        temp <- optim(nhpp.nll, par=par, x=x, u=u_list[j],M=M, hessian=F)$par
        opt <- rbind(opt,temp)
        test <- temp[1]
        mu <- c(mu,test)
        test <- temp[2]
        sig <- c(sig,test)
        test <- temp[3]
        xi <- c(xi,test)
      }
      opt <- opt[-1,]
      mu <- mu[-1]
      sig <- sig[-1]
      xi <- xi[-1]
      opt_median <- c(median(mu), median(sig), median(xi))
      thetahat[i,]<-opt_median

      I_list <- list()
      Iinv_list <- list()
      for(j in 1:upper_number)
      {
        #It used to be this one
        #if(opt[j,3]>-0.5)
        #change to avoid all zero in I_list[[j]], sometimes avoid "the integral is probably divergent"
        if(opt[j,3]>-0.5&abs(opt[j,3])>1e-05)
        {
          I_list[[j]] <- E.Info.Mat(theta=opt[j,], u=u_list[j], M=M)$EIM
          Iinv_list[[j]] <- solve(I_list[[j]])
        } else{I_list[[j]]<-matrix(0,3,3)}
      }
      I_temp1 <- NULL
      I_temp2 <- NULL
      I_temp3 <- NULL
      I_temp4 <- NULL
      I_temp5 <- NULL
      I_temp6 <- NULL
      I_temp7 <- NULL
      I_temp8 <- NULL
      I_temp9 <- NULL
      for(j in 1:upper_number) {
        I_temp1 <- c(I_temp1, I_list[[j]][1,1])
        I_temp2 <- c(I_temp2, I_list[[j]][1,2])
        I_temp3 <- c(I_temp3, I_list[[j]][1,3])
        I_temp4 <- c(I_temp4, I_list[[j]][2,1])
        I_temp5 <- c(I_temp5, I_list[[j]][2,2])
        I_temp6 <- c(I_temp6, I_list[[j]][2,3])
        I_temp7 <- c(I_temp7, I_list[[j]][3,1])
        I_temp8 <- c(I_temp8, I_list[[j]][3,2])
        I_temp9 <- c(I_temp9, I_list[[j]][3,3])
      }
      #Iinv_median <- matrix(0,3,3)
      if(MLE_model == "median") {
        I_median <- matrix(0,3,3)
        I_median[1,1] <- median(I_temp1)
        I_median[1,2] <- median(I_temp2)
        I_median[1,3] <- median(I_temp3)
        I_median[2,1] <- median(I_temp4)
        I_median[2,2] <- median(I_temp5)
        I_median[2,3] <- median(I_temp6)
        I_median[3,1] <- median(I_temp7)
        I_median[3,2] <- median(I_temp8)
        I_median[3,3] <- median(I_temp9)
        I[[i]] <- I_median
      } else if(MLE_model == "min") {
        I_min <- matrix(0,3,3)
        I_min[1,1] <- min(I_temp1)
        I_min[1,2] <- min(I_temp2)
        I_min[1,3] <- min(I_temp3)
        I_min[2,1] <- min(I_temp4)
        I_min[2,2] <- min(I_temp5)
        I_min[2,3] <- min(I_temp6)
        I_min[3,1] <- min(I_temp7)
        I_min[3,2] <- min(I_temp8)
        I_min[3,3] <- min(I_temp9)
        I[[i]] <- I_min
      } else if(MLE_model == "second_min") {
        I_second_min <- matrix(0,3,3)
        I_second_min[1,1] <- min(I_temp1[I_temp1>min(I_temp1)])
        I_second_min[1,2] <- min(I_temp2[I_temp2>min(I_temp2)])
        I_second_min[1,3] <- min(I_temp3[I_temp3>min(I_temp3)])
        I_second_min[2,1] <- min(I_temp4[I_temp4>min(I_temp4)])
        I_second_min[2,2] <- min(I_temp5[I_temp5>min(I_temp5)])
        I_second_min[2,3] <- min(I_temp6[I_temp6>min(I_temp6)])
        I_second_min[3,1] <- min(I_temp7[I_temp7>min(I_temp7)])
        I_second_min[3,2] <- min(I_temp8[I_temp8>min(I_temp8)])
        I_second_min[3,3] <- min(I_temp9[I_temp9>min(I_temp9)])
        I[[i]] <- I_second_min
      }
      if(max(I[[i]]) == 0 & min(I[[i]]) == 0) {Iinv[[i]] <- matrix(0,3,3)
      } else Iinv[[i]]<-solve(I[[i]])
    }
    else {
      opt<-optim(nhpp.nll, par=par, x=x, u=u[i],M=M, hessian=F)
      thetahat[i,]<-opt$par

      ###
      # Deal with xi<-0.5
      if(thetahat[i,3]>-0.5)
        ###
      {
        I[[i]]<- E.Info.Mat(theta=opt$par, u=u[i], M=M)$EIM
        Iinv[[i]]<-solve(I[[i]])
      }
      else{I[[i]]<-Iinv[[i]]<-matrix(0,3,3)}
    }
  }

  Wcov<-list()
  Wcov1<-NULL
  for(i in 1:k)
  {
    Wcov[[i]]<-matrix(0,3,3)
    for(j in 1:k)
    {
      Wcov[[i]]<-cbind(Wcov[[i]],Iinv[[min(i,j)]])
    }
    Wcov1<-rbind(Wcov1,Wcov[[i]])
  }
  Wcov1<-Wcov1[,-c(1:3)]

  CovT<-Wcov1

  Cov.mu<-CovT[seq(1,3*k,by=3),seq(1,3*k,by=3)]
  Cov.sig<-CovT[seq(2,3*k,by=3),seq(2,3*k,by=3)]
  Cov.xi<-CovT[seq(3,3*k,by=3),seq(3,3*k,by=3)]

  return(list(mle=thetahat, Cov.all=CovT, Cov.mu=Cov.mu, Cov.sig = Cov.sig, Cov.xi=Cov.xi))
}



#######################################################################################################
# Joint.MLE.Expl
#######################################################################################################
#
# Details:
#
# Calculates the MLEs of the rate parameter, and joint asymptotic covariance matrix of these MLEs
# over a range of thresholds as supplied by the user.
#
# Arguments:
#
# x - vector of data
# u - vector of thresholds. If not supplied, then k thresholds between quantiles (q1, q2) will be used
# k - number of thresholds to consider if u not supplied
# q1 -  lowest quantile for the threshold sequence
# q2 -  upper quantile limit for the threshold sequence (q2 itself is not used as a threshold, but rather the uppermost
#       threshold will be at the (q2-1/k) quantile)
# param - character specifying "InvRate" or "Rate" for either inverse rate parameter / rate parameter, respectively
#
# Value: (returns:)
#
# mle - vector of MLEs above the supplied thresholds
# cov - joint asymptotic covariance matrix of these MLEs

#' @title Joint.MLE.Expl
#' @description Calculates the MLEs of the rate parameter, and joint asymptotic covariance matrix of these MLEs
#' over a range of thresholds as supplied by the user.
#'
#' @param x - vector of data
#' @param u - vector of thresholds. If not supplied, then k thresholds between quantiles (q1, q2) will be used
#' @param k - number of thresholds to consider if u not supplied
#' @param q1 - lowest quantile for the threshold sequence
#' @param q2 -  upper quantile limit for the threshold sequence (q2 itself is not used as a threshold, but rather the uppermost
#'       threshold will be at the (q2-1/k) quantile)
#' @param param - character specifying "InvRate" or "Rate" for either inverse rate parameter / rate parameter, respectively
#'
#' @return mle - vector of MLEs above the supplied thresholds
#' @return cov - joint asymptotic covariance matrix of these MLEs
#'

Joint.MLE.Expl<-function(x, u=NULL, k, q1, q2=1, param)
{
  if(!is.element(param, c("InvRate","Rate"))){stop("param should be one of InvRate or Rate")}
  if(!is.null(u))
  {
    k<-length(u)
    x<-x[x>u[1]]
    # add threshold above all data, to "close" final region
    u<-c(u,max(x)+1)
  }
  else{u<-quantile(x,seq(q1,q2,len=k+1))}

  I<-n<-m<-thetahat<-NULL

  for(i in 1:k)
  {
    if(param=="InvRate")
    {
      thetahat[i]<-mean(x[x>=u[i]]- u[i])
    }
    else if(param=="Rate")
    {
      thetahat[i]<-1/mean(x[x>=u[i]]- u[i])
    }
    n[i]<-length(x[x>=u[i]&x<=u[i+1]])
  }
  for(i in 1:k)
  {
    m[i]<-sum(n[i:k])
    I[i]<- 1/thetahat[i]^2
  }

  Tcov<-matrix(0,k,k)
  for(i in 1:k)
  {
    for(j in 1:k)
    {
      Tcov[i,j]<-1/(I[min(i,j)]*m[min(i,j)])
    }
  }
  CovT<-Tcov
  return(list(mle=thetahat, Cov=CovT))
}



###################################################################################
# norm.LRT
###################################################################################
#
# Details:
#
# Evaluates the likelihood ratio statistics for testing white noise
#
# Arguments:
#
# x - vector of white noise process (WNP, usually normalized estimates of xi or the exponential rate paramter 1/eta)
# u - vector of thresholds that are associated to the WNP

#' @title normal distribution likelihood ratio statistics
#' @description Evaluates the likelihood ratio statistics for testing white noise
#'
#' @param x - vector of white noise process (WNP, usually normalized estimates of xi or the exponential rate paramter 1/eta)
#' @param u - vector of thresholds that are associated to the WNP
#'
#' @return v
#' @return likelihood ratio
#'

norm.LRT<-function(x,u)
{
  l<-length(u)
  v<-u[-c(1)] # means two or more obs available for std dev calculation

  #Edit by Mengchen on 22 May 2016
  for(i in 1:length(x))
  {
    if(is.na(x[i]) || is.nan(x[i])) x[i] <- 0
  }
  #End of Edit

  lr<-NULL
  for(i in 1:length(v))
  {
    n1<-length(x[u<=v[i]])
    num<-nll.norm(theta=c(mean(x[u<=v[i]]),sd(x[u<=v[i]])*sqrt((n1-1)/n1)),x=x[u<=v[i]])
    den<-nll.norm(theta=c(0,1),x=x[u<=v[i]])
    lr[i]<--2*(num-den)
  }
  return(cbind(v,lr))
}



###################################################################################
# nll.norm - negative log liklihood for the normal distribution
###################################################################################

nll.norm<-function(theta,x){
  if(theta[2]<0){return(10e10)}
  else{
    return(-sum(dnorm(x,mean=theta[1],sd=theta[2],log=T)))
  }
}



###################################################################################
# C1
###################################################################################
#
# Details:
#
# Produces "Contrast matrix" with (1,-1) elements running down the two diagonals
#
# Arguments:
#
# k - number of columns (k-1 = number of rows)
#
# Value: (returns:)
#
# (k-1) x k contrast matrix

C1<-function(k)
{
  C<-matrix(0,k-1,k)
  for(i in 1:(k-1))
  {
    C[i,i]<-1
    C[i,i+1]<--1
  }
  return(C)
}


###################################################################################
# Functions for the expected information matrix used in Joint.MLE.NHPP
###################################################################################
#
# Details:
#
# Functions with names of form "d2ldmu2" are derivatives (with expectation over the random NUMBER of points
# already incorporated). Functions with names of form "i.d2ldmu2", are in a form ready for integration; this integration
# yields the expectation over the random LOCATIONS of the points.
#
# Arguments:
#
# x - dummy variable over which to integrate
# mu - location parameter
# sig - scale parameter
# xi - shape parameter
# u - threshold above which NHPP model assumed
# M - number of superpositions or "blocks" / "years" the process corresponds to

d2ldmu2<-function(x, mu, sig, xi, u, M)
{
  c1<-(1+(xi/sig)*(u-mu))^(-1/xi)

  p1<-M*c1*(xi*(1+xi)/sig^2)*(1+(xi/sig)*(x-mu))^(-2)
  p2<- -M*((1+xi)/sig^2)*(1+(xi/sig)*(u-mu))^(-1/xi-2)
  return(p1+p2)
}

i.d2ldmu2<-function(x, mu, sig, xi, u, M)
{
  d2ldmu2(x=x,mu=mu,sig=sig,xi=xi,u=u,M=M)* (1+(xi/sig)*(u-mu))^(1/xi) * (1/sig)*(1+(xi/sig)*(x-mu))^(-1/xi-1)
}



d2ldmudsig<-function(x, mu, sig, xi, u, M)
{
  c1<-(1+(xi/sig)*(u-mu))^(-1/xi)

  p1<-M*c1*(xi*(1+xi)/sig^3)*(x-mu)*(1+(xi/sig)*(x-mu))^(-2)
  p2<-M*c1*(-(1+xi)/sig^2)*(1+(xi/sig)*(x-mu))^(-1)
  p3<- M*(1/sig^2)*(1+(xi/sig)*(u-mu))^(-1/xi-1)
  p4<- -M*((1+xi)/sig^3)*(u-mu)*(1+(xi/sig)*(u-mu))^(-1/xi-2)
  return(p1+p2+p3+p4)
}

i.d2ldmudsig<-function(x, mu, sig, xi, u, M)
{
  d2ldmudsig(x=x,mu=mu,sig=sig,xi=xi,u=u, M=M)* (1+(xi/sig)*(u-mu))^(1/xi) * (1/sig)*(1+(xi/sig)*(x-mu))^(-1/xi-1)
}

#########

d2ldmudxi<-function(x, mu, sig, xi, u, M)
{
  c1<-(1+(xi/sig)*(u-mu))^(-1/xi)

  p1<-M*c1* (1/sig)*(1+(xi/sig)*(x-mu))^(-1)
  p2<-M*c1*(-(1+xi)/sig^2)*(x-mu)*(1+(xi/sig)*(x-mu))^(-2)
  p3.1<-(1/sig)*((-1/xi^2)*log(1+(xi/sig)*(u-mu)) + (1/sig)*(1/xi+1)*(u-mu)*(1+(xi/sig)*(u-mu))^(-1))
  p3<-M*p3.1*(1+(xi/sig)*(u-mu))^(-1/xi-1)

  return(p1+p2+p3)
}

i.d2ldmudxi<-function(x, mu, sig, xi, u, M)
{
  d2ldmudxi(x=x,mu=mu,sig=sig,xi=xi,u=u, M=M)* (1+(xi/sig)*(u-mu))^(1/xi) * (1/sig)*(1+(xi/sig)*(x-mu))^(-1/xi-1)
}



##########


d2ldsig2<-function(x, mu, sig, xi, u, M)
{
  c1<-(1+(xi/sig)*(u-mu))^(-1/xi)

  p1<-M*c1*(-2*(1+xi)/sig^3)*(x-mu)*(1+(xi/sig)*(x-mu))^(-1)
  p2<-M*c1*(xi*(1+xi)/sig^4)*((x-mu)^2)*(1+(xi/sig)*(x-mu))^(-2)
  p3<- M*(2/sig^3)*(u-mu)*(1+(xi/sig)*(u-mu))^(-1/xi-1)
  p4<- -M*((1+xi)/sig^4)*((u-mu)^2)*(1+(xi/sig)*(u-mu))^(-1/xi-2)
  p5<-M*c1*1/sig^2
  return(p1+p2+p3+p4+p5)
}

i.d2ldsig2<-function(x, mu, sig, xi, u, M)
{
  d2ldsig2(x=x,mu=mu,sig=sig,xi=xi,u=u, M=M)* (1+(xi/sig)*(u-mu))^(1/xi) * (1/sig)*(1+(xi/sig)*(x-mu))^(-1/xi-1)
}



#########

d2ldsigdxi<-function(x, mu, sig, xi, u, M)
{
  c1<-(1+(xi/sig)*(u-mu))^(-1/xi)

  p1<- M*c1*((x-mu)/sig^2)*(1+(xi/sig)*(x-mu))^(-1)
  p2<-M*c1*(-(1+xi)/sig^3)*((x-mu)^2)*(1+(xi/sig)*(x-mu))^(-2)
  p3.1<-((u-mu)/sig^2)*((-1/xi^2)*log(1+(xi/sig)*(u-mu)) + (1/sig)*(1/xi+1)*(u-mu)*(1+(xi/sig)*(u-mu))^(-1))
  p3<-M*p3.1*(1+(xi/sig)*(u-mu))^(-1/xi-1)

  return(p1+p2+p3)
}


i.d2ldsigdxi<-function(x, mu, sig, xi, u, M)
{
  d2ldsigdxi(x=x,mu=mu,sig=sig,xi=xi,u=u, M=M)* (1+(xi/sig)*(u-mu))^(1/xi) * (1/sig)*(1+(xi/sig)*(x-mu))^(-1/xi-1)
}


#########

d2ldxi2<-function(x, mu, sig, xi, u, M)
{
  c1<-(1+(xi/sig)*(u-mu))^(-1/xi)

  p1<- M*c1*(-2/xi^3)*log(1+(xi/sig)*(x-mu))#
  p2<- M*c1*2*((x-mu)/(sig*xi^2))*(1+(xi/sig)*(x-mu))^(-1)#
  p3<-M*c1*((1/xi+1)/sig^2)*((x-mu)^2)*(1+(xi/sig)*(x-mu))^(-2)#

  p4.1<-((-1/xi^2)*log(1+(xi/sig)* (u-mu)) + (1/sig)*(1/xi)*(u-mu)*(1+(xi/sig)*(u-mu))^(-1))#
  p4<- -(p4.1^2)*M*(1+(xi/sig)*(u-mu))^(-1/xi)#

  p5.1<-((2/xi^3)*log(1+(xi/sig)*(u-mu)) - (2/sig)*(1/xi^2)*(u-mu)*(1+(xi/sig)*(u-mu))^(-1)- (1/sig^2)*(1/xi)*((u-mu)^2)*(1+(xi/sig)*(u-mu))^(-2))#
  p5<- (p5.1)*M*(1+(xi/sig)*(u-mu))^(-1/xi)#

  return(p1+p2+p3+p4+p5)
}

i.d2ldxi2<-function(x, mu, sig, xi, u, M)
{
  d2ldxi2(x=x,mu=mu,sig=sig,xi=xi,u=u, M=M)* (1+(xi/sig)*(u-mu))^(1/xi) * (1/sig)*(1+(xi/sig)*(x-mu))^(-1/xi-1)
}

## version of d2ldxi2 and i.d2ldxi2 with limit as xi \to 0. (This derivative has most trouble with small absolute
## values of xi.) Also possible to take such limits in other derivatives, but not implemented as they are generally
## less problematic.


d2ldxi2.xi0<-function(x, mu, sig, xi, u, M)
{
  c1<-exp(-(u-mu)/sig)
  q1<-((x-mu)/sig)
  q2<-((u-mu)/sig)

  p1<--((2/3)*q1^3-q1^2)*M*c1
  p2<--(-(2/3)*q2^3+(1/4)*q2^4)*M*c1

  return(p1+p2)
}

i.d2ldxi2.xi0<-function(x, mu, sig, xi, u, M)
{
  d2ldxi2.xi0(x=x,mu=mu,sig=sig,xi=xi,u=u, M=M)* (1/sig)*exp(-(x-u)/sig)
}



#############################################################################################
# E.Info.Mat
#############################################################################################
#
# Details:
#
# Calculates the numerically-integrated expected information matrix for an NHPP with specified parameters
#
# Arguments:
#
# theta - vector of parameters (mu,sigma, xi)
# u - threshold for NHPP
# M - number of superpositions or "blocks" / "years" the process corresponds to
#
# Value: (returns:)
#
# EIM - expected information matrix
# Errors - vector of errors from the numerical integration of the 6 unique components

#' @title E.Info.Mat
#' @description Calculates the numerically-integrated expected information matrix for an NHPP with specified parameters
#'
#' @param theta - vector of parameters (mu,sigma, xi)
#' @param u - threshold for NHPP
#' @param M - number of superpositions or "blocks" / "years" the process corresponds to
#'
#' @return EIM - expected information matrix
#' @return Errors - vector of errors from the numerical integration of the 6 unique components
#'

E.Info.Mat<-function(theta, u, M)
{
  if(theta[3]<0)
  {
    up<-theta[1] - theta[2]/theta[3]
  }
  else{up<-Inf}

  Ed2ldmu2<-integrate(i.d2ldmu2,low=u,upper=up,mu=theta[1], sig=theta[2],
                      xi=theta[3], u=u,M=M, abs.tol=0)
  Ed2ldmudsig<-integrate(i.d2ldmudsig,low=u,upper=up,mu=theta[1], sig=theta[2],
                         xi=theta[3], u=u,M=M, abs.tol=0)
  Ed2ldmudxi<-integrate(i.d2ldmudxi,low=u,upper=up,mu=theta[1], sig=theta[2],
                        xi=theta[3], u=u,M=M, abs.tol=0)
  Ed2ldsig2<-integrate(i.d2ldsig2,low=u,upper=up,mu=theta[1], sig=theta[2],
                       xi=theta[3], u=u,M=M, abs.tol=0)
  Ed2ldsigdxi<-integrate(i.d2ldsigdxi,low=u,upper=up,mu=theta[1], sig=theta[2],
                         xi=theta[3], u=u,M=M, abs.tol=0)
  if(abs(theta[3])>0.0001)
  {
    Ed2ldxi2<-integrate(i.d2ldxi2,low=u,upper=up,mu=theta[1], sig=theta[2],
                        xi=theta[3], u=u,M=M, abs.tol=0)
  }
  else{
    Ed2ldxi2<-integrate(i.d2ldxi2.xi0,low=u,upper=up,mu=theta[1], sig=theta[2],
                        xi=theta[3], u=u,M=M, abs.tol=0)
  }

  Errors<-c(Ed2ldmu2$abs.err,Ed2ldmudsig$abs.err,Ed2ldmudxi$abs.err,Ed2ldmudsig$abs.err,Ed2ldsig2$abs.err,
            Ed2ldsigdxi$abs.err,Ed2ldmudxi$abs.err,Ed2ldsigdxi$abs.err,Ed2ldxi2$abs.err)

  EIM<--matrix(c(Ed2ldmu2$val,Ed2ldmudsig$val,Ed2ldmudxi$val,Ed2ldmudsig$val,Ed2ldsig2$val,Ed2ldsigdxi$val,
                 Ed2ldmudxi$val,Ed2ldsigdxi$val,Ed2ldxi2$val), byrow=T, nrow=3)

  return(list(EIM=EIM,Errors=Errors))
}



#############################################################################################
# K.ForwardStop
#############################################################################################
#
# Details:
#
# calculate the cutoff k based on ForwardStop
#
# Arguments:
#
# pvalues - the p-values for each threshold in the threshold list
# alpha - the pre-specified significance level, usually take 0.01 or 0.05
#
# Value: (returns:)
#
# forward - if the given null hypothesis should be rejected. Take 0 as no rejection and 1 as make rejection
#

#' @title K.ForwardStop
#' @description calculate the cutoff k based on ForwardStop
#'
#' @param pvalues - the p-values for each threshold in the threshold list
#' @param alpha - the pre-specified significance level, usually take 0.01 or 0.05
#'
#' @return forward - if the given null hypothesis should be rejected. Take 0 as no rejection and 1 as make rejection
#' @export
#'

K.ForwardStop<-function(pvalues, alpha = 0.05)
{
  pvalues <- rev(pvalues)

  plength <- length(pvalues)
  Y <- 1
  for(i in 1:plength) {
    singleY <- -log(1-pvalues[i])
    Y <- c(Y, singleY)
  }
  Y <- Y[-1]

  forward.list <- 1
  for(i in 1:plength) {
    singleForward.list <- sum(Y[1:i])/i
    forward.list <- c(forward.list, singleForward.list)
  }
  forward.list <- forward.list[-1]

  for(f in 1:length(forward.list)) {
    if(forward.list[f] >= 1) forward.list[f] = 1
    if(forward.list[f] <= 0) forward.list[f] = 0
  }

  forward.list2 <- 1
  for(i in 1:plength) {
    singleForward.list2 <- forward.list[i] <= alpha
    forward.list2 <- c(forward.list2, singleForward.list2)
  }
  forward.list2 <- forward.list2[-1]

  forward <- max(forward.list2)
  returnlist <- list()
  returnlist$forward <- forward
  returnlist$forward.list <- forward.list
  return(returnlist)
}



#############################################################################################
# K.StrongStop
#############################################################################################
#
# Details:
#
# calculate the cutoff k based on StrongStop
#
# Arguments:
#
# pvalues - the p-values for each threshold in the threshold list
# alpha - the pre-specified significance level, usually take 0.01 or 0.05
#
# Value: (returns:)
#
# strong - if the given null hypothesis should be rejected. Take 0 as no rejection and 1 as make rejection
#

#' @title K.StrongStop
#' @description calculate the cutoff k based on StrongStop
#'
#' @param pvalues - the p-values for each threshold in the threshold list
#' @param alpha - the pre-specified significance level, usually take 0.01 or 0.05
#'
#' @return strong - if the given null hypothesis should be rejected. Take 0 as no rejection and 1 as make rejection
#' @export
#'

K.StrongStop<-function(pvalues, alpha = 0.05)
{
  pvalues <- rev(pvalues)

  plength <- length(pvalues)
  Y <- 1
  for(i in 1:plength) {
    singleY <- -log(pvalues[i])
    Y <- c(Y, singleY)
  }
  Y <- Y[-1]

  Z <- 1
  for(i in 1:plength) {
    singleZ <- sum(Y[i:plength]/seq(i,plength,1))
    Z <- c(Z, singleZ)
  }
  Z <- Z[-1]

  q <- 1
  for(i in 1:plength) {
    singleq <- exp(-Z[i]) * plength / i
    q <- c(q, singleq)
  }
  q <- q[-1]

  strong.list <- 1
  for(i in 1:plength) {
    singleStrong.list <- (q[i] <= alpha)
    strong.list <- c(strong.list, singleStrong.list)
  }
  strong.list <- strong.list[-1]

  strong <- max(strong.list)

  returnlist <- list()
  returnlist$strong <- strong
  returnlist$strong.list <- q
  return(returnlist)
}



#############################################################################################
# Data_Randomization
#############################################################################################
#
# Details:
#
# add the randomization of data to reduce the limitation of minimum measured value
#
# Arguments:
#
# x - data vector
# random_digit - the precision of data vector
# round_digit - the precision of return randomization data vector
# seed - presetting random seed
#
# Values: (Returns:)
#
# RandomData - randomization data vector
#

#' @title Data_Randomization
#' @description add the randomization of data to reduce the limitation of minimum measured value
#'
#' @param x - data vector
#' @param random_digit - the precision of data vector
#' @param round_digit - the precision of return randomization data vector
#' @param seed - presetting random seed
#'
#' @return RandomData - randomization data vector
#' @export
#'
#' @examples
#' x <- rnorm(100)
#' x <- sort(x)
#' x <- round(x*10)/10
#' Data_Randomization(x)
#'

Data_Randomization <- function(x, random_digit = 2, round_digit = 4, seed = 1) {
  set.seed(seed)
  Random <- runif(length(x), -5/(10^random_digit), 4.9999/(10^random_digit))
  RandomData <- round((x+Random) * (10^round_digit)) / (10^round_digit)
  RandomData <- sort(RandomData)
  return(RandomData)
}



#############################################################################################
# Threshold_Randomization
#############################################################################################
#
# Details:
#
# add the randomization of thresholds to increase the randomization of thresholds
#
# Arguments:
#
# x - data vector
# thresholds - thresholds vector
# seed - presetting random seed
#
# Values: (Returns:)
#
# thresholds - randomization thresholds vector
#

#' @title Threshold_Randomization
#' @description add the randomization of thresholds to increase the randomization of thresholds
#'
#' @param x - data vector
#' @param thresholds - thresholds vector
#' @param seed - presetting random seed
#'
#' @return thresholds - randomization thresholds vector
#' @export
#'
#' @examples
#' set.seed(1)
#' xnorm<-abs(rnorm(5000))
#' thresholds_xnorm <- Threshold_Generator_Uniform(xnorm)
#' thresholds_xnorm
#' Threshold_Randomization(xnorm,thresholds_xnorm,seed = 1)
#'


Threshold_Randomization <- function(x, thresholds, seed = 1) {
  set.seed(seed)
  for(i in 1:length(thresholds)) {
    if(thresholds[i] <= min(x)) Random <- thresholds[i]
    else Random <- runif(1,max(x[x<=thresholds[i]]),thresholds[i])
    thresholds[i] <- Random
  }
  return(thresholds)
}



#############################################################################################
# Thershold_Generator_Uniform
#############################################################################################
#
# Details:
#
# generate a candidate thresholds list for given data set based on uniform distance among the distribution
#
# Arguments:
#
# x - data vector
# wetCuttingPoint - the minimum precision which can be considered as "rain"
# ks - number of candidate thresholds
#
# Values: (Returns:)
#
# thresholds - candidate thresholds list
#

#' @title Thershold_Generator_Uniform
#' @description generate a candidate thresholds list for given data set based on uniform distance among the distribution
#'
#' @param x - data vector
#' @param wetCuttingPoint - the minimum precision which can be considered as "rain"
#' @param ks - number of candidate thresholds
#'
#' @return thresholds - candidate thresholds list
#' @export
#'
#' @examples
#' set.seed(1)
#' xnorm<-abs(rnorm(5000))
#' Threshold_Generator_Uniform(xnorm)
#'

Threshold_Generator_Uniform<-function(x, wetCuttingPoint=0.1, ks = NULL) {
  set_number <- 8

  if(is.null(ks)) ks <- KvsN(x, wetCuttingPoint)

  q1 <- min(x)
  q2 <- (length(x) - set_number) / length(x)

  threshold_q1 <- q1
  threshold_q2 <- quantile(x, probs = q2)

  thresholds <- seq(threshold_q1 , threshold_q2 , len = ks + 1)

  return(thresholds)
}



#############################################################################################
# Thershold_Generator_Density
#############################################################################################
#
# Details:
#
# generate a candidate thresholds list for given data set based on variant distance among the distribution density
#
# Arguments:
#
# x - data vector
# wetCuttingPoint - the minimum precision which can be considered as "rain"
# ks - number of candidate thresholds
# q - presetting value to control the density (usually set between 1 and 2)
#     (the bigger q is, the more intensive thresholds in tail part)
#
# Values: (Returns:)
#
# thresholds - candidate thresholds list
#

#' @title Thershold_Generator_Density
#' @description generate a candidate thresholds list for given data set based on variant distance among the distribution density
#'
#' @param x - data vector
#' @param wetCuttingPoint - the minimum precision which can be considered as "rain"
#' @param ks - number of candidate thresholds
#' @param q - presetting value to control the density (usually set between 1 and 2)
#'     (the bigger q is, the more intensive thresholds in tail part)
#'
#' @return thresholds - candidate thresholds list
#' @export
#'
#' @examples
#' set.seed(1)
#' xnorm<-abs(rnorm(5000))
#' Threshold_Generator_Density(xnorm)
#'

Threshold_Generator_Density <-function(x, wetCuttingPoint=0.1, ks = NULL, q = NULL){
  set_number <- 8
  set_unique_number <- 5

  if(is.null(ks)) ks <- KvsN(x, wetCuttingPoint)

  de <- density(x)

  q1 <- min(x)
  q2 <- (length(x) - set_number) / length(x)
  threshold_q1 <- q1
  threshold_q2 <- quantile(x, probs = q2)

  thresholds <- seq(threshold_q1 , threshold_q2 , len = ks + 1)

  threshold_middle <- 1
  for(i in 1:(length(thresholds)-1)) {
    threshold_middle <- c(threshold_middle,((thresholds[i]+thresholds[i+1])/2))
  }
  threshold_middle <- threshold_middle[-1]

  interval_density <- 1
  for(i in 1:length(threshold_middle)) {
    interval_density <- c(interval_density,de$y[which.min(abs(de$x - threshold_middle[i]))])
  }
  interval_density <- interval_density[-1]

  if(is.null(q)) { q <- 1}

  balanced_interval_density <- 1
  for(i in 1:length(interval_density)) {
    balanced_interval_density <- c(balanced_interval_density,interval_density[i] ^ q)
  }
  balanced_interval_density <- balanced_interval_density[-1]

  sum_balanced_interval_density <- sum(balanced_interval_density)
  new_interval_prob <- 1
  for(i in 1:length(interval_density)) {
    new_interval_prob <- c(new_interval_prob, sum(balanced_interval_density[1:i])/sum_balanced_interval_density)
  }
  new_interval_prob <- new_interval_prob[-1]

  changed_new_interval_prob <- 1
  for(i in 1:length(new_interval_prob)) {
    changed_new_interval_prob <- c(changed_new_interval_prob, new_interval_prob[i]*(q2-q1)+q1)
  }
  changed_new_interval_prob <- changed_new_interval_prob[-1]

  new_thresholds <- quantile(x, probs = changed_new_interval_prob)

  new_thresholds <- c(threshold_q1, new_thresholds)
  new_data <- x[x >= threshold_q1]
  new_data <- new_data[new_data <= threshold_q2]
  for(i in length(new_thresholds):2) {
    regress <- 0

    record <- new_data[new_data <= new_thresholds[i]]
    record <- record[record >= new_thresholds[i-1]]
    unique_record <- unique(record)

    record2 <- new_data[new_data <= new_thresholds[i]]
    record2 <- record2[record2 >= threshold_q1]
    unique_record2 <- unique(record2)

    if(length(record2) < set_number || length(unique_record2) < set_unique_number) {
      new_thresholds <- new_thresholds[-i]
    }else if(length(record) < set_number || length(unique_record) < set_unique_number) {
      regress <- 1
      while(regress) {
        new_thresholds[i-1] <- max(new_data[new_data < new_thresholds[i-1]])
        record3 <- new_data[new_data <= new_thresholds[i]]
        record3 <- record3[record3 >= new_thresholds[i-1]]
        unique_record3 <- unique(record3)

        if(length(record3) >= set_number && length(unique_record3) >= set_unique_number) {
          regress <- FALSE
        }
      }
    }
  }

  return(new_thresholds)
}

#############################################################################################
# KvsN
#############################################################################################
#
# Details:
#
# Function to determine an appropriate number of thresholds
#
# Arguments:
#
# x - data vector
# wetCuttingPoint - the minimum precision which can be considered as "rain"
#
# Values: (Returns:)
#
# ks - an appropriate number of thresholds
#
#' @title KvsN
#'
#' @description Function to determine an appropriate number of thresholds
#'
#' @param x - data vector
#' @param wetCuttingPoint - the minimum precision which can be considered as "rain"
#'
#' @return ks - an appropriate number of thresholds
#'

###Function to determine an appropriate number of thresholds
KvsN = function(x, wetCuttingPoint = 0.1){
  KagainstObs<-c(10,12,15,23,30,33,37,40,50)
  (nObs <- 1000*(1:length(KagainstObs)))
  lmKvsN<-lm(KagainstObs ~ nObs)
  (Kpredicted=ceiling(predict(lmKvsN,newdata=data.frame(nObs=sum(x>wetCuttingPoint)))[1]))
  (ks <- ceiling(min(max(Kpredicted,9),51)))

  return(ks)
}

