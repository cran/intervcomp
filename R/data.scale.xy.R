#' @title Data Transformation for Testing Equality of Variability Measures
#'
#' @description Data transformation for testing equality of variability measures (mean absolute deviation from median or variance).
#' @export
#' @param x A vector specifying the first sample.
#' @param y A vector specifying the second sample.
#' @param scale.option A character string specifying the transformation on each observation; "Levene.med" for Levene's transformation (absolute difference from the sample median), "Levene.med.0" for Levene's transformation with zero removal for odd sample(s), "Levene.med.00" for Levene's transformation with zero removal for both odd and even sample(s), and "Variance" for squared difference from the sample mean. The default option is "Levene.med.0".
#' @param scale.adj A boolean specifying whether or not any constant should be multiplied to each transformed observation. This is useful to obtain unbiased estimates. The default option is TRUE.
#' @param scale.00 A scale to be applied for an even sample when "Levene.med.00" is chosen. The default option is 2, but the square root of 2 is another viable option.
#' @param paired A boolean specifying whether or not the samples are paired or independent. The default option is FALSE.
#'
#' @return 
#' \item{trans.x}{Transformed first sample.}
#' \item{trans.y}{Transformed second sample.}
#' \item{rho}{Sample correlation coefficient. For independent samples, it returns a NULL.}
#'
#' @examples 
#'  set.seed(123)
#'  x<-runif(10)
#'  y<-runif(15)
#'  data.scale.xy(x, y, scale.option="Levene.Med.0", scale.adj=TRUE, paired=FALSE)
#'
#' @importFrom stats cor median

data.scale.xy<-function(x, y, scale.option=c("Levene.Med.0","Levene.Med","Levene.Med.00","Variance"), scale.adj=TRUE, scale.00=2, paired=FALSE)
{
 scale.option<-match.arg(scale.option)	
 rho<-NULL
 if(paired==TRUE) ###Removing NAs
 {
  xy<-rbind(x,y)
  na.xy<-apply(xy, 2, is.na)
  na.sum<-colSums(na.xy)
  if(sum(na.sum) > 0)
  {
   na.cols<-which(na.sum > 0)
   if((dim(xy)[2] - length(na.cols)) < 3)
   {
   	stop("not enough data.")
   }
   xy.new<-xy[,-na.cols]
   x<-xy.new[1,]
   y<-xy.new[2,]
   rho<-cor(x,y)
  } 
 }
 
 n.x<-sum(!is.na(x))
 n.y<-sum(!is.na(y))
 
 m.x<-floor(n.x/2)*2+1
 m.y<-floor(n.y/2)*2+1
 
 if(scale.option=="Variance")
 {
  mean.x<-mean(x, na.rm=TRUE)
  mean.y<-mean(y, na.rm=TRUE)
  trans.x<-(x-mean.x)^2
  trans.y<-(y-mean.y)^2
  if(scale.adj==TRUE)
  {
   trans.x<-1/(1-1/n.x)*trans.x
   trans.y<-1/(1-1/n.y)*trans.y
  }
 }
 
 if(scale.option=="Levene.Med")
 {
  median.x<-median(x, na.rm=TRUE)
  median.y<-median(y, na.rm=TRUE)
  trans.x<-abs(x-median.x)
  trans.y<-abs(y-median.y)
  if(scale.adj==TRUE)
  {
   trans.x<-((m.x/(m.x-1))^(5/6))*trans.x
   trans.y<-((m.y/(m.y-1))^(5/6))*trans.y   
  }
 }
 
 if(scale.option=="Levene.Med.0")
 {
  if(n.x%%2==0) ### x has an even sample size
  {
   median.x<-median(x, na.rm=TRUE)  
   trans.x<-abs(x-median.x)
   if(scale.adj==TRUE)
   {
    trans.x<-((m.x/(m.x-1))^(5/6))*trans.x   	       
   }
  }
  
  if(n.y%%2==0) ### y has an even sample size
  {
   median.y<-median(y, na.rm=TRUE)  
   trans.y<-abs(y-median.y)
   if(scale.adj==TRUE)
   {
    trans.y<-((m.y/(m.y-1))^(5/6))*trans.y   	       
   }
  }
   
  if(n.x%%2==1)  ### x has an odd sample size
  {
   median.x<-median(x, na.rm=TRUE)  
   trans.x<-abs(x-median.x)
   if(scale.adj==TRUE)
   {
    trans.x<-(1-1/n.x)*((m.x/(m.x-1))^(5/6))*trans.x   	       
   }
   zero.x<-which.min(trans.x)[1]
   #trans.x<-trans.x[-zero.x]
   trans.x[zero.x]<-NA
  }
   
  if(n.y%%2==1)  ### y has an odd sample size
  {
   median.y<-median(y, na.rm=TRUE)  
   trans.y<-abs(y-median.y)
   if(scale.adj==TRUE)
   {
    trans.y<-(1-1/n.y)*((m.y/(m.y-1))^(5/6))*trans.y   	       
   }
   zero.y<-which.min(trans.y)[1]
   #trans.y<-trans.y[-zero.y]
   trans.y[zero.y]<-NA
  }   
 }   
 
 if(scale.option=="Levene.Med.00")
 {
  if(n.x%%2==0) ### x has an even sample size
  {
   median.x<-median(x, na.rm=TRUE)  
   trans.x<-abs(x-median.x)
   if(scale.adj==TRUE)
   {
    trans.x<-(1-1/n.x)*((m.x/(m.x-1))^(5/6))*trans.x   	       
   }
   zero.mult.x<-order(trans.x)[c(1,2)]
   order.x<-order(x[zero.mult.x])
   zero.x<-zero.mult.x[order.x[1]]
   mult.x<-zero.mult.x[order.x[2]]  
   trans.x[mult.x]<-scale.00*trans.x[mult.x] 
   #trans.x<-trans.x[-zero.x] 
   trans.x[zero.x]<-NA       
  }
  
  if(n.y%%2==0) ### y has an even sample size
  {
   median.y<-median(y, na.rm=TRUE)  
   trans.y<-abs(y-median.y)
   if(scale.adj==TRUE)
   {
    trans.y<-(1-1/n.y)*((m.y/(m.y-1))^(5/6))*trans.y   	       
   }
   zero.mult.y<-order(trans.y)[c(1,2)]
   order.y<-order(y[zero.mult.y])
   zero.y<-zero.mult.y[order.y[1]]
   mult.y<-zero.mult.y[order.y[2]]  
   trans.y[mult.y]<-scale.00*trans.y[mult.y] 
   #trans.y<-trans.y[-zero.y]
   trans.y[zero.y]<-NA             
  }
   
  if(n.x%%2==1)  ### x has an odd sample size
  {
   median.x<-median(x, na.rm=TRUE)  
   trans.x<-abs(x-median.x)
   if(scale.adj==TRUE)
   {
    trans.x<-(1-1/n.x)*((m.x/(m.x-1))^(5/6))*trans.x   	       
   }
   zero.x<-which.min(trans.x)[1]
   #trans.x<-trans.x[-zero.x]
   trans.x[zero.x]<-NA
  }
   
  if(n.y%%2==1)  ### y has an odd sample size
  {
   median.y<-median(y, na.rm=TRUE)  
   trans.y<-abs(y-median.y)
   if(scale.adj==TRUE)
   {
    trans.y<-(1-1/n.y)*((m.y/(m.y-1))^(5/6))*trans.y   	       
   }
   zero.y<-which.min(trans.y)[1]
   #trans.y<-trans.y[-zero.y]
   trans.y[zero.y]<-NA   
  }   
 }     
 
 return(list(trans.x=trans.x, trans.y=trans.y, rho=rho))
}


