#' @title Two-Sample T-Test for Equality of Variability Mesures
#'
#' @description Two-sample t-test for comparing variability measures.
#' @export
#' @param x A vector specifying the first sample.
#' @param y A vector specifying the second sample.
#' @param alternative A character string specifying the alternative hypothesis; "two.sided" for two-tailed, "less" for left-tailed, and "greater" for right-tailed alternative hypothesis. The default option is "two.sided".
#' @param scale.option A character string specifying the transformation on each observation; "Levene.med" for Levene's transformation (absolute difference from the sample median), "Levene.med.0" for Levene's transformation with zero removal for odd sample(s), "Levene.med.00" for Levene's transformation with zero removal for both odd and even sample(s), and "Variance" for squared difference from the sample mean. The default option is "Levene.med.0".
#' @param scale.adj A boolean specifying whether or not any constant should be multiplied to each transformed observation. This is useful to obtain unbiased estimates. The default option is TRUE.
#' @param scale.00 A scale to be applied for an even sample when "Levene.med.00" is chosen. The default option is 2, but the square root of 2 is another viable option.
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
#' @param tol A numeric specifying the cut-off value for a numerical zero. The default option is 1e-7.
#'
#' @return 
#' \item{CI}{A data frame displaying the effect size estimate, lower and upper bound of the effect size confidence interval, and the degrees of freedom.}
#' \item{ICI}{A data frame displaying the sample mean for both x and y, lower and upper bounds of the inferential confidence intervals, the degrees of freedom of the inferential confidence intervals, the inferential significance level, and the significance level.}
#' \item{Statistic}{A data frame displaying the test statistic and p-value of the generalized two-sample t-test.}
#' \item{Ind.CI}{A data frame displaying the sample mean for both x and y, lower and upper bounds of the individual confidence intervals, the degrees of freedom of the individual confidence intervals, the sample standard deviations, and the significance level.}
#' \item{Effect.Sizes}{A data frame displaying Cohen's d and log ratio for the mean and variability measure comparisons.}
#'
#' @examples
#'  set.seed(123) 
#'  x<-rexp(10)
#'  y<-rexp(15) 
#'  two.sample.var(x,y,alternative="two.sided",scale.option="Levene.Med.0",scale.adj=TRUE,paired=FALSE)
#'
#' @importFrom stats na.omit var

two.sample.var<-function(x,y,alternative=c("two.sided", "less", "greater"),
scale.option=c("Levene.Med","Levene.Med.0","Levene.Med.00","Variance"), scale.adj=TRUE, scale.00=2,
paired=FALSE,alpha=0.05,plot.ci=TRUE,plot.ici=TRUE,ici.interval=NULL,add.individual.ci=TRUE, by.ici=0.5, 
xlab.ici="", rounds=3, rounds.plot=NULL, ci.interval=NULL, by.ci=0.5, name.x=NULL, name.y=NULL, 
pool.mean=FALSE,logit.interval=NULL, tol=1e-7)
{
 scale.data.xy<-data.scale.xy(x,y,scale.option=scale.option,scale.adj=scale.adj,scale.00=scale.00,paired=paired)
 new.x<-scale.data.xy$trans.x
 new.y<-scale.data.xy$trans.y
 orig.rho<-scale.data.xy$rho
 xn <- na.omit(x)
 yn <- na.omit(y)
 cohen.m <- (mean(xn) - mean(yn))/sqrt((var(xn)*(length(xn)-1)+var(yn)*(length(yn)-1))/(length(xn)+length(yn)-2))
 logratio.m <- NA
 if(min(xn) > 0 && min(yn) > 0)
 {
  logratio.m <- log(mean(xn)/mean(yn))
 }
 newx <- na.omit(new.x)
 newy <- na.omit(new.y)
 cohen.v <- (mean(newx) - mean(newy))/sqrt((var(newx)*(length(newx)-1)+var(newy)*(length(newy)-1))/(length(newx)+length(newy)-2))
 logratio.v <- log(mean(newx)/mean(newy))
 
 results<-two.sample.mean(new.x,new.y,transformation="log",alternative=alternative,
 paired=paired,alpha=alpha,plot.ci=plot.ci,plot.ici=plot.ici,ici.interval=ici.interval,add.individual.ci=add.individual.ci, 
 by.ici=by.ici, xlab.ici=xlab.ici, rounds=rounds, rounds.plot=rounds.plot, ci.interval=ci.interval, by.ci=by.ci, 
 name.x=name.x, name.y=name.y, pool.mean=pool.mean,logit.interval=logit.interval, rho=orig.rho, tol=tol)
 ess <- c(cohen.m,logratio.m,cohen.v,logratio.v)
 names(ess) <- c("Cohen.d.Mean","Log.Ratio.Mean","Cohen.d.Var","Log.Ratio.Var")
 Analysis <- list(CI = results$CI, ICI = results$ICI, Statistic = results$Statistic, Ind.CI = results$Ind.CI, 
 Effect.Sizes = ess) 
 return(Analysis)
}

