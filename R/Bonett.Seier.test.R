#' @title Bonett-Seier Test for Equality of Variability Measures
#'
#' @description Bonett-Seier t-test for comparing variability measures of two independent samples.
#' @export
#' @param x A vector specifying the first sample.
#' @param y A vector specifying the second sample.
#' @param alternative A character string specifying the alternative hypothesis; "two.sided" for two-tailed, "less" for left-tailed, and "greater" for right-tailed alternative hypothesis. The default option is "two.sided".
#' @param alpha A numeric specifying the significance level. The default option is 0.05.
#'
#' @return 
#' \item{Statistic}{The test statistic.}
#' \item{p.value}{The p-value.}
#' \item{Estimate}{The ratio of variability measures.}
#' \item{Lower.CI}{The lower bound of the confidence interval.}
#' \item{Upper.CI}{The upper bound of the confidence interval.}
#'
#' @examples 
#'  set.seed(123)
#'  x<-runif(10)
#'  y<-runif(15)
#'  Bonett.Seier.test(x,y,"two.sided",0.05)
#'
#' @importFrom stats median pnorm qnorm sd

Bonett.Seier.test<-function(x,y,alternative=c("two.sided", "less", "greater"),alpha=0.05)
{
 eta.x<-median(x, na.rm=TRUE)
 eta.y<-median(y, na.rm=TRUE)
 tau.x<-mean(abs(x-eta.x), na.rm=TRUE)
 tau.y<-mean(abs(y-eta.y), na.rm=TRUE)
 n.x<-sum(!is.na(x))
 n.y<-sum(!is.na(y))
 c.x<-n.x/(n.x-1)
 c.y<-n.y/(n.y-1)
 mu.x<-mean(x, na.rm=TRUE)
 mu.y<-mean(y, na.rm=TRUE)
 delta.x<-(mu.x-eta.x)/tau.x
 delta.y<-(mu.y-eta.y)/tau.y
 sigma.x<-sd(x, na.rm=TRUE)
 sigma.y<-sd(y, na.rm=TRUE)
 gamma.x<-sigma.x^2/tau.x^2
 gamma.y<-sigma.y^2/tau.y^2
 var.ln.tau.x<-(delta.x^2 + gamma.x - 1)/n.x
 var.ln.tau.y<-(delta.y^2 + gamma.y - 1)/n.y
 
 est<-(c.x*tau.x)/(c.y*tau.y)
 se<-exp(qnorm(alpha/2, lower.tail=FALSE)*sqrt(var.ln.tau.x+var.ln.tau.y))
 conf.upper<-est*se
 conf.lower<-est/se
 
 q<-log(1/est)/sqrt(var.ln.tau.x+var.ln.tau.y)
 
 if(alternative=="less")
 {
  p.value<-pnorm(q,lower.tail=TRUE) 
 }
 else if(alternative=="greater")
 {
  p.value<-pnorm(q,lower.tail=FALSE) 
 } 
 else #if(alternative=="two.sided")
 {
  p.value<-2*pnorm(abs(q),lower.tail=FALSE)   	
 }
 
 return(list(Statistic=q, p.value=p.value, Estimate=est, Lower.CI=conf.lower, Upper.CI=conf.upper))
} 
