#' @title Range-Preserving Two-Sample T-Test for Equality of Means
#'
#' @description Range-preserving two-sample t-test for comparing means.
#' @export
#' @param x A vector specifying the first sample.
#' @param y A vector specifying the second sample.
#' @param transformation A character string specifying the transformation on the sample means; "none" for real-valued data, "log" for positive-valued data, and "logit" for data in a fixed interval. The default option is "none".
#' @param alternative A character string specifying the alternative hypothesis; "two.sided" for two-tailed, "less" for left-tailed, and "greater" for right-tailed alternative hypothesis. The default option is "two.sided".
#' @param paired A boolean specifying whether or not the samples are paired or independent. The default option is FALSE.
#' @param alpha A numeric specifying the significance level. The default option is 0.05.
#' @param plot.ci A boolean specifying whether or not the deviation between the two samples should be plotted. The default option is TRUE.
#' @param plot.ici A boolean specifying whether or not the inferential confidence intervals for the two samples should be plotted. The default option is TRUE.
#' @param ici.interval A vector with two elements specifying the range of values for the plot of inferential confidence intervals. The default option is NULL, in which case, an appropriate range is automatically specified.
#' @param add.individual.ci A boolean specifying whether or not the confidence intervals with the confidence level implied by \code{alpha} should also be plotted. The default option is TRUE.
#' @param by.ici A numeric specifying te scales in the plot for ici.interval. The default option is 0.5.
#' @param xlab.ici A character string specifying the x axis labels used for the inferential confidence intervals and confidence intervals. The default option is "".
#' @param rounds A numeric specifying the number of decimal places to be rounded. The default option is 3.
#' @param rounds.plot A numeric specifying the number of decimal places to be rounded. The default option is NULL, in which case it is set equal to \code{rounds}.
#' @param ci.interval A vector with two elements specifying the range of values for the plot of confidence intervals. The default option is NULL, in which case, an appropriate range is automatically specified.
#' @param by.ci A numeric specifying te scales in the plot for ci.interval. The default option is 0.5.
#' @param name.x A character string specifying the label for the x variable. The default option is NULL, in which case, it is set to "x".
#' @param name.y A character string specifying the label for the y variable. The default option is NULL, in which case, it is set to "y".
#' @param pool.mean A boolean specifying whether or not the sample means should be pooled for the degrees of freedom. The default option is FALSE.
#' @param logit.interval A vector with two elements specifying the lower and upper bound if logit transformation is used. The default option is NULL.
#' @param rho A numeric specifying the correlation coefficient between the two samples. The default option is NULL, in which case, the sample correlation coefficient is estimated automatically.
#' @param tol A numeric specifying the cut-off value for a numerical zero. The default option is 1e-7.
#'
#' @return 
#' \item{CI}{A data frame displaying the effect size estimate, lower and upper bound of the effect size confidence interval, and the degrees of freedom.}
#' \item{ICI}{A data frame displaying the sample mean for both x and y, lower and upper bounds of the inferential confidence intervals, the degrees of freedom of the inferential confidence intervals, the inferential significance level, and the significance level.}
#' \item{Statistic}{A data frame displaying the test statistic and p-value of the generalized two-sample t-test.}
#' \item{Ind.CI}{A data frame displaying the sample mean for both x and y, lower and upper bounds of the individual confidence intervals, the degrees of freedom of the individual confidence intervals, the sample standard deviations, and the significance level.}
#' \item{Effect.Sizes}{A data frame displaying Cohen's d and either log ratio (for transformation="log") or log odds (for transformation="logit").}
#'
#' @examples 
#'  set.seed(123)
#'  x<-runif(10)
#'  y<-runif(15)
#'  two.sample.mean(x,y,"logit","two.sided",paired=FALSE, ici.interval=c(0,1),
#'   by.ici=0.2, logit.interval=c(0,1), rounds=2, name.x="xvar", name.y="yvar")
#'
#' @importFrom grDevices dev.new
#' @importFrom graphics abline axis box plot rect text title
#' @importFrom stats cor na.omit pt qt sd uniroot var

two.sample.mean<-function(x,y,transformation=c("none","log","logit"),alternative=c("two.sided", "less", "greater"),
paired=FALSE,alpha=0.05,plot.ci=TRUE,plot.ici=TRUE,ici.interval=NULL,add.individual.ci=TRUE, by.ici=0.5, 
xlab.ici="", rounds=3, rounds.plot=NULL, ci.interval=NULL, by.ci=0.5, name.x=NULL, name.y=NULL, 
pool.mean=FALSE,logit.interval=NULL, rho=NULL, tol=1e-7)
{
 transformation <- match.arg(transformation)
 alternative <- match.arg(alternative)
 
 mu.x <- mean(x, na.rm=TRUE)
 mu.y <- mean(y, na.rm=TRUE)
 v.x <- var(x, na.rm=TRUE)
 v.y <- var(y, na.rm=TRUE)
 s.x <- sd(x, na.rm=TRUE)
 s.y <- sd(y, na.rm=TRUE)
 n.x <- sum(!is.na(x))
 n.y <- sum(!is.na(y))
 conf.level <- 1 - alpha
  
 if(is.null(logit.interval) && transformation=="logit")
 {
  stop("For logit transformation, logit.interval needs to be specified.")
 }
 
 if(pool.mean==TRUE)
 {
  weight <- (n.x/v.x)/(n.x/v.x + n.y/v.y)
  mu.pool <- weight*mu.x + (1-weight)*mu.y
 }
 
 if(transformation=="none")
 {
  g.mu.x <- mu.x
  g.mu.y <- mu.y
  g.mu.p.x <- 1
  g.mu.p.y <- 1
 }
 
 if(transformation=="log")
 {
  g.mu.x <- log(mu.x)
  g.mu.y <- log(mu.y)
  if(pool.mean==FALSE)
  {  
   g.mu.p.x <- 1/mu.x
   g.mu.p.y <- 1/mu.y
  }
 if(pool.mean==TRUE)
 {
   g.mu.p.x <- 1/mu.pool
   g.mu.p.y <- g.mu.p.x
  }    
 }
 
 if(transformation=="logit")
 {
  a <- min(logit.interval)
  b <- max(logit.interval)
  g.mu.x <- log(mu.x-a) - log(b-mu.x)
  g.mu.y <- log(mu.y-a) - log(b-mu.y)
  if(pool.mean==FALSE)
  {
   g.mu.p.x <- 1/(mu.x-a) + 1/(b-mu.x)
   g.mu.p.y <- 1/(mu.y-a) + 1/(b-mu.y)
  }
  if(pool.mean==TRUE)
  {
   g.mu.p.x <- 1/(mu.pool-a) + 1/(b-mu.pool)
   g.mu.p.y <- g.mu.p.x 
  } 
 }
 
 ### Effect size calculation
 xn <- na.omit(x)
 yn <- na.omit(y)
 cohen.m <- (mean(xn) - mean(yn))/sqrt((var(xn)*(length(xn)-1)+var(yn)*(length(yn)-1))/(length(xn)+length(yn)-2))
 logratio.m <- NA
 if(transformation=="log")
 {
  logratio.m <- log(mean(xn)/mean(yn))
 }
 logodds.m <- NA
 if(transformation=="logit")
 {
  logodds.m <- log(((mean(xn)-a)/(b-mean(xn)))/((mean(yn)-a)/(b-mean(yn))))
 } 
  
 ess <- c(cohen.m,logratio.m,logodds.m)
 names(ess) <- c("Cohen.d","Log.Ratio","Log.Odds")
  
 if(paired==FALSE)
 {
  statistic <- (g.mu.x - g.mu.y)/sqrt(g.mu.p.x^2*v.x/n.x + g.mu.p.y^2*v.y/n.y)
  df.numerator <- (g.mu.p.x^2*v.x/n.x + g.mu.p.y^2*v.y/n.y)^2
  df.denominator <- (g.mu.p.x^2*v.x/n.x)^2/(n.x-1) + (g.mu.p.y^2*v.y/n.y)^2/(n.y-1)
  df <- df.numerator/df.denominator
  
  ### general confidence interval for the difference ###
  estimate <- g.mu.x - g.mu.y
  halfwidth.ci <- qt(alpha/2, df=df, lower.tail=FALSE)*sqrt(g.mu.p.x^2*v.x/n.x + g.mu.p.y^2*v.y/n.y)
  Lower <- estimate - halfwidth.ci
  Upper <- estimate + halfwidth.ci
 } 
 else ###paired==TRUE
 {
  if(n.x != n.y)
  {
   stop("For paired t-test, the two sample sizes must be equal.")
  }
  n <- n.x
  if(is.null(rho))
  {
   rho.xy <- cor(x,y, use="pairwise.complete.obs")
  }
  else
  {
   rho.xy <- rho
  } 
  statistic <- sqrt(n)*(g.mu.x - g.mu.y)/sqrt(g.mu.p.x^2*v.x - 2*g.mu.p.x*g.mu.p.y*rho.xy*s.x*s.y + g.mu.p.y^2*v.y)
  df <- n-1
  
  ### general confidence interval for the difference ###
  estimate <- g.mu.x - g.mu.y
  halfwidth.ci <- qt(alpha/2, df=df, lower.tail=FALSE)*sqrt(g.mu.p.x^2*v.x - 2*g.mu.p.x*g.mu.p.y*rho.xy*s.x*s.y + g.mu.p.y^2*v.y)/sqrt(n)
  Lower <- estimate - halfwidth.ci
  Upper <- estimate + halfwidth.ci  
 }
 
 ### finding alpha.star ###
 if(paired==FALSE)
 {
  trhs.init <- qt(alpha/2, df=df, lower.tail=FALSE)*sqrt(g.mu.p.x^2*s.x^2/n.x + g.mu.p.y^2*s.y^2/n.y)/(g.mu.p.x*s.x/sqrt(n.x) + g.mu.p.y*s.y/sqrt(n.y))
  alpha.star.x <- 2*pt(trhs.init, df=n.x-1, lower.tail=FALSE)
  alpha.star.y <- 2*pt(trhs.init, df=n.y-1, lower.tail=FALSE)
  if(n.x==n.y)
  {
   alpha.star <- alpha.star.x
  }
  else ### when n.x is not equal to n.y
  {
   alpha.star.min <- min(alpha.star.x, alpha.star.y)
   alpha.star.max <- max(alpha.star.x, alpha.star.y)
   
   crit.fun <- function(alpha.star)
   {
    a1 <- qt(alpha/2, df=df, lower.tail=FALSE)
    b1 <- qt(alpha.star/2, df=n.x-1, lower.tail=FALSE)*g.mu.p.x*s.x/sqrt(n.x) + qt(alpha.star/2, df=n.y-1, lower.tail=FALSE)*g.mu.p.y*s.y/sqrt(n.y)
    b2 <- sqrt(g.mu.p.x^2*s.x^2/n.x + g.mu.p.y^2*s.y^2/n.y)
    return(a1 - b1/b2)
   }
   
   alpha.star <- uniroot(crit.fun, interval=c(alpha.star.min,alpha.star.max),tol=tol)$root 
  }
 }
 else ### paired==TRUE
 {
  trhs <- qt(alpha/2, df=df, lower.tail=FALSE)*sqrt(g.mu.p.x^2*s.x^2 - 2*g.mu.p.x*g.mu.p.y*rho.xy*s.x*s.y + g.mu.p.y^2*s.y^2)/(g.mu.p.x*s.x + g.mu.p.y*s.y)
  alpha.star <- 2*pt(trhs, df=n-1, lower.tail=FALSE)
 }
 
 if(alternative=="two.sided")
 {
  p.value <- 2*pt(abs(statistic), df=df, lower.tail=FALSE)
  if(transformation=="none")
  {
   Upper.x.transform <- mu.x + qt(alpha.star/2, df=n.x-1, lower.tail=FALSE)*s.x/sqrt(n.x)
   Lower.x.transform <- mu.x - qt(alpha.star/2, df=n.x-1, lower.tail=FALSE)*s.x/sqrt(n.x)
   Upper.y.transform <- mu.y + qt(alpha.star/2, df=n.y-1, lower.tail=FALSE)*s.y/sqrt(n.y)
   Lower.y.transform <- mu.y - qt(alpha.star/2, df=n.y-1, lower.tail=FALSE)*s.y/sqrt(n.y)
   
   Upper.x.ind.transform <- mu.x + qt(alpha/2, df=n.x-1, lower.tail=FALSE)*s.x/sqrt(n.x)
   Lower.x.ind.transform <- mu.x - qt(alpha/2, df=n.x-1, lower.tail=FALSE)*s.x/sqrt(n.x)
   Upper.y.ind.transform <- mu.y + qt(alpha/2, df=n.y-1, lower.tail=FALSE)*s.y/sqrt(n.y)
   Lower.y.ind.transform <- mu.y - qt(alpha/2, df=n.y-1, lower.tail=FALSE)*s.y/sqrt(n.y)   
   
   Lower.transform <- Lower
   Upper.transform <- Upper
   Estimate.transform <- estimate
  } 
  if(transformation=="log")
  {
   Upper.x.transform <- mu.x*exp(qt(alpha.star/2, df=n.x-1, lower.tail=FALSE)*s.x/(mu.x*sqrt(n.x)))
   Lower.x.transform <- mu.x*exp(-qt(alpha.star/2, df=n.x-1, lower.tail=FALSE)*s.x/(mu.x*sqrt(n.x)))
   Upper.y.transform <- mu.y*exp(qt(alpha.star/2, df=n.y-1, lower.tail=FALSE)*s.y/(mu.y*sqrt(n.y)))
   Lower.y.transform <- mu.y*exp(-qt(alpha.star/2, df=n.y-1, lower.tail=FALSE)*s.y/(mu.y*sqrt(n.y)))
   
   Upper.x.ind.transform <- mu.x*exp(qt(alpha/2, df=n.x-1, lower.tail=FALSE)*s.x/(mu.x*sqrt(n.x)))
   Lower.x.ind.transform <- mu.x*exp(-qt(alpha/2, df=n.x-1, lower.tail=FALSE)*s.x/(mu.x*sqrt(n.x)))
   Upper.y.ind.transform <- mu.y*exp(qt(alpha/2, df=n.y-1, lower.tail=FALSE)*s.y/(mu.y*sqrt(n.y)))
   Lower.y.ind.transform <- mu.y*exp(-qt(alpha/2, df=n.y-1, lower.tail=FALSE)*s.y/(mu.y*sqrt(n.y)))   
   
   Lower.transform <- exp(Lower)
   Upper.transform <- exp(Upper)
   Estimate.transform <- exp(estimate)  
  } 
  if(transformation=="logit")
  {
   exp.part.upper.x <- exp(qt(alpha.star/2, df=n.x-1, lower.tail=FALSE)*(b-a)*s.x/((mu.x-a)*(b-mu.x)*sqrt(n.x)))
   exp.part.lower.x <- exp(-qt(alpha.star/2, df=n.x-1, lower.tail=FALSE)*(b-a)*s.x/((mu.x-a)*(b-mu.x)*sqrt(n.x)))   
   exp.part.upper.y <- exp(qt(alpha.star/2, df=n.y-1, lower.tail=FALSE)*(b-a)*s.y/((mu.y-a)*(b-mu.y)*sqrt(n.y)))
   exp.part.lower.y <- exp(-qt(alpha.star/2, df=n.y-1, lower.tail=FALSE)*(b-a)*s.y/((mu.y-a)*(b-mu.y)*sqrt(n.y)))      
   Upper.x.transform <- (a+b*((mu.x-a)/(b-mu.x))*exp.part.upper.x)/(1+((mu.x-a)/(b-mu.x))*exp.part.upper.x)
   Lower.x.transform <- (a+b*((mu.x-a)/(b-mu.x))*exp.part.lower.x)/(1+((mu.x-a)/(b-mu.x))*exp.part.lower.x)
   Upper.y.transform <- (a+b*((mu.y-a)/(b-mu.y))*exp.part.upper.y)/(1+((mu.y-a)/(b-mu.y))*exp.part.upper.y)
   Lower.y.transform <- (a+b*((mu.y-a)/(b-mu.y))*exp.part.lower.y)/(1+((mu.y-a)/(b-mu.y))*exp.part.lower.y)
   
   exp.part.upper.x.ind <- exp(qt(alpha/2, df=n.x-1, lower.tail=FALSE)*(b-a)*s.x/((mu.x-a)*(b-mu.x)*sqrt(n.x)))
   exp.part.lower.x.ind <- exp(-qt(alpha/2, df=n.x-1, lower.tail=FALSE)*(b-a)*s.x/((mu.x-a)*(b-mu.x)*sqrt(n.x)))   
   exp.part.upper.y.ind <- exp(qt(alpha/2, df=n.y-1, lower.tail=FALSE)*(b-a)*s.y/((mu.y-a)*(b-mu.y)*sqrt(n.y)))
   exp.part.lower.y.ind <- exp(-qt(alpha/2, df=n.y-1, lower.tail=FALSE)*(b-a)*s.y/((mu.y-a)*(b-mu.y)*sqrt(n.y)))      
   Upper.x.ind.transform <- (a+b*((mu.x-a)/(b-mu.x))*exp.part.upper.x.ind)/(1+((mu.x-a)/(b-mu.x))*exp.part.upper.x.ind)
   Lower.x.ind.transform <- (a+b*((mu.x-a)/(b-mu.x))*exp.part.lower.x.ind)/(1+((mu.x-a)/(b-mu.x))*exp.part.lower.x.ind)
   Upper.y.ind.transform <- (a+b*((mu.y-a)/(b-mu.y))*exp.part.upper.y.ind)/(1+((mu.y-a)/(b-mu.y))*exp.part.upper.y.ind)
   Lower.y.ind.transform <- (a+b*((mu.y-a)/(b-mu.y))*exp.part.lower.y.ind)/(1+((mu.y-a)/(b-mu.y))*exp.part.lower.y.ind)   
   
   Lower.transform <- exp(Lower)
   Upper.transform <- exp(Upper)
   Estimate.transform <- exp(estimate)    
  }     
 }
 if(alternative=="less")
 {
  p.value <- pt(statistic, df=df, lower.tail=TRUE)
 }
 if(alternative=="greater")
 {
  p.value <- pt(statistic, df=df, lower.tail=FALSE)
 }
 
 if(is.null(name.x)) name.x <- "x"
 if(is.null(name.y)) name.y <- "y" 
 
 if(transformation=="none") 
 {
  if(xlab.ici=="")
  {
   ci.string <- "Difference"
  }
  else
  {
   ci.string <- paste("Difference in ", xlab.ici, sep="")
  } 
  conname <- paste(name.x, "-", name.y, sep="")
 } 
 if(transformation=="log") 
 {
  if(xlab.ici=="")
  { 	
   ci.string <- "Ratio"  
  }
  else
  {
   ci.string <- paste("Ratio for ", xlab.ici, sep="")   
  } 
  conname <- paste(name.x, "/", name.y, sep="")
 } 
 if(transformation=="logit")
 { 
  if(xlab.ici=="")
  {  	
  ci.string <- "Odds Ratio" 
  }
  else
  {
   ci.string <- paste("Odds Ratio for ", xlab.ici, sep="")   
  }      
  conname <- paste("odds.", name.x, "/", "odds.", name.y, sep="")
 } 
 
 if(is.null(rounds.plot))
 {
  if(rounds < Inf)
  {
   rounds.plot <- rounds
  }
  else
  {
   rounds.plot <- 8
  }
 }
 if(rounds.plot == Inf)
 {
  rounds.plot <- 8
 }
 
 if(plot.ci == TRUE)
 {
  dev.new(width=5,height=5)
  text.ci <- paste(paste(conf.level * 100, "%", sep=""), "Confidence Interval")
  if(is.null(ci.interval))
  {
   Width <- max(Upper.transform - Lower.transform)
   Wide.hwidth <- 1.1*Width/2
   Midpoints <- c((Upper.transform + Lower.transform)/2)	
   upper.limit.ci <- ceiling(max(Midpoints + Wide.hwidth))
   lower.limit.ci <- floor(min(Midpoints - Wide.hwidth))
  }  
  else
  {
   upper.limit.ci <- max(ci.interval)
   lower.limit.ci <- min(ci.interval)
  }
  
  plot(Estimate.transform, 1, axes = FALSE, 
  xlab = ci.string, ylab = "", xlim = c(lower.limit.ci, upper.limit.ci), ylim = c(0,2), type="n")
  if(transformation=="none")
  {
   abline(v = 0, lty = 3, lwd = 2)        	
  }
  else
  {
   abline(v = 1, lty = 3, lwd = 2)           	
  }
  
  bw.ci <- 0.3  
  p.ci <- 0.4
  cex.ci <- 1.2
  rect(Estimate.transform, (1-bw.ci), Estimate.transform, (1+bw.ci), density=20, lwd = 2)
  rect(Lower.transform, (1-bw.ci), Upper.transform, (1+bw.ci), density=20, lwd = 2)
  text(Lower.transform, (1+p.ci), sprintf(paste("%.", rounds.plot, "f", sep=""), Lower.transform), cex=cex.ci)
  text(Upper.transform, (1+p.ci), sprintf(paste("%.", rounds.plot, "f", sep=""), Upper.transform), cex=cex.ci)     
 
  axis(1, at = seq(lower.limit.ci, upper.limit.ci, by=by.ci))
  axis(2, at = 1, labels = conname, las=0)
  box()
  if(transformation=="logit")
  {
   title(main = c(text.ci, ci.string, paste("Transformation: ", transformation, " with min=", a, " and max=", b, sep=""), 
   paste("DF: ", round(df, rounds), ", p-value: ", sprintf(paste("%.", rounds.plot, "f", sep=""), p.value), sep=""))) 
  }
  else
  {
   title(main = c(text.ci, ci.string, paste("Transformation: ", transformation, sep=""), 
   paste("DF: ", round(df, rounds), ", p-value: ", sprintf(paste("%.", rounds.plot, "f", sep=""), p.value), sep=""))) 
  }  
 }

 if(plot.ici == TRUE)
 {
  dev.new(width=5,height=5)
  text.ici <- paste("Inferential Confidence Intervals with alpha=",alpha,sep="")    
   
  if(is.null(ici.interval))
  {
   Width <- max(c(Upper.x.transform - Lower.x.transform), c(Upper.y.transform - Lower.y.transform), 
   c(Upper.x.ind.transform - Lower.x.ind.transform), c(Upper.y.ind.transform - Lower.y.ind.transform))
   Wide.hwidth <- 1.1*Width/2
   Midpoints.x <- c((Upper.x.transform + Lower.x.transform)/2)
   Midpoints.y <- c((Upper.y.transform + Lower.y.transform)/2)
   Midpoints.x.ind <- c((Upper.x.ind.transform + Lower.x.ind.transform)/2)
   Midpoints.y.ind <- c((Upper.y.ind.transform + Lower.y.ind.transform)/2)   
   max.mid <- max(Midpoints.x, Midpoints.y, Midpoints.x.ind, Midpoints.y.ind)
   min.mid <- min(Midpoints.x, Midpoints.y, Midpoints.x.ind, Midpoints.y.ind)
   upper.limit.ici <- ceiling(max.mid + Wide.hwidth)
   lower.limit.ici <- floor(min.mid - Wide.hwidth)
   ici.interval <- c(lower.limit.ici, upper.limit.ici)
  }  
  else
  {
   upper.limit.ici <- max(ici.interval)
   lower.limit.ici <- min(ici.interval)
  }
     
  plot(c(mu.x, mu.y), 1:2, axes = FALSE, 
  xlab = xlab.ici, ylab = "", xlim = c(lower.limit.ici, upper.limit.ici), ylim = c(0,3), type="n")

  bw.ici <- 0.3
  p.ici <- 0.4  
  cex.ici <- 1
  
  rect(mu.x, (1-bw.ici), mu.x, (1+bw.ici), density=20, lwd = 2)
  rect(Lower.x.transform, (1-bw.ici), Upper.x.transform, (1+bw.ici), density=20, lwd = 2)
  text(Lower.x.transform, (1+p.ici), sprintf(paste("%.", rounds.plot, "f", sep=""), Lower.x.transform), cex=cex.ici)
  text(Upper.x.transform, (1+p.ici), sprintf(paste("%.", rounds.plot, "f", sep=""), Upper.x.transform), cex=cex.ici)     
   
  if(add.individual.ci == TRUE)
  {
   rect(Lower.x.ind.transform, (1-bw.ici), Upper.x.ind.transform, (1+bw.ici), density=0, lwd = 2, lty=2)   	
  }    
   
  rect(mu.y, (2-bw.ici), mu.y, (2+bw.ici), density=20, lwd = 2)
  rect(Lower.y.transform, (2-bw.ici), Upper.y.transform, (2+bw.ici), density=20, lwd = 2)
  text(Lower.y.transform, (2+p.ici), sprintf(paste("%.", rounds.plot, "f", sep=""), Lower.y.transform), cex=cex.ici)
  text(Upper.y.transform, (2+p.ici), sprintf(paste("%.", rounds.plot, "f", sep=""), Upper.y.transform), cex=cex.ici)     
   
  if(add.individual.ci == TRUE)
  {
   rect(Lower.y.ind.transform, (2-bw.ici), Upper.y.ind.transform, (2+bw.ici), density=0, lwd = 2, lty=2)   	
  }
  
  axis(1, at = seq(lower.limit.ici, upper.limit.ici, by=by.ici))
  axis(2, at = 1:2, labels = c(name.x,name.y), las=0)
  box()
  if(transformation=="logit")
  { 
   title(main = c(text.ici, paste("Transformation: ", transformation, " with min=", a, " and max=", b, sep=""), 
   paste("DF(", name.x,"): ", (n.x-1), ", DF(", name.y,"): ", (n.y-1)), paste("alpha*: ", sprintf(paste("%.", rounds.plot, "f", sep=""), alpha.star), 
   ", p-value: ", sprintf(paste("%.", rounds.plot, "f", sep=""), p.value), sep="")))
  }
  else
  { 
   title(main = c(text.ici, paste("Transformation: ", transformation, sep=""), 
   paste("DF(", name.x,"): ", (n.x-1), ", DF(", name.y,"): ", (n.y-1)), paste("alpha*: ", sprintf(paste("%.", rounds.plot, "f", sep=""), alpha.star), 
   ", p-value: ", sprintf(paste("%.", rounds.plot, "f", sep=""), p.value), sep="")))
  }   
 }   

 ICI <- data.frame(
 Estimate.x = round(mu.x, rounds), 
 Lower.x = round(Lower.x.transform, rounds), 
 Upper.x = round(Upper.x.transform, rounds),
 df.x = n.x-1,
 Estimate.y = round(mu.y, rounds),
 Lower.y = round(Lower.y.transform, rounds), 
 Upper.y = round(Upper.y.transform, rounds), 
 df.y = n.y-1,
 alpha.star = round(alpha.star, rounds),
 alpha = alpha)
 row.names(ICI) <- conname  
 colnames(ICI) <- c(paste("Estimate.",name.x,sep=""), paste("Lower.",name.x,sep=""), paste("Upper.",name.x,sep=""), paste("df.",name.x,sep=""), 
 paste("Estimate.",name.y,sep=""), paste("Lower.",name.y,sep=""), paste("Upper.",name.y,sep=""), paste("df.",name.y,sep=""), "alpha*", "alpha") 
   
 CI <- data.frame(
 Estimate = round(Estimate.transform, rounds), 
 Lower = round(Lower.transform, rounds), 
 Upper = round(Upper.transform, rounds), 
 df = round(df, rounds))
 row.names(CI) <- conname

 Stats.p.value <- data.frame(
 Statistic = round(statistic, rounds), 
 p.value = round(p.value, rounds)) 
 row.names(Stats.p.value) <- conname 
  	
 Ind.CI <- data.frame(
 Estimate.x = round(mu.x, rounds), 
 Lower.x = round(Lower.x.ind.transform, rounds), 
 Upper.x = round(Upper.x.ind.transform, rounds),
 df.x = n.x-1,
 sd.x = round(s.x, rounds),
 Estimate.y = round(mu.y, rounds),
 Lower.y = round(Lower.y.ind.transform, rounds), 
 Upper.y = round(Upper.y.ind.transform, rounds), 
 df.y = n.y-1,
 sd.y = round(s.y, rounds), 
 alpha = alpha)
 row.names(Ind.CI) <- conname  
 colnames(Ind.CI) <- c(paste("Estimate.",name.x,sep=""), paste("Lower.",name.x,sep=""), paste("Upper.",name.x,sep=""), paste("df.",name.x,sep=""), 
 paste("sd.",name.x,sep=""), 
 paste("Estimate.",name.y,sep=""), paste("Lower.",name.y,sep=""), paste("Upper.",name.y,sep=""), paste("df.",name.y,sep=""), 
 paste("sd.",name.y,sep=""), "alpha")   

 Analysis <- list(CI = CI, ICI = ICI, Statistic = Stats.p.value, Ind.CI = Ind.CI, Effect.Sizes = ess)  
 return(Analysis)
}

