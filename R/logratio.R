#' @title Log Ratio Analysis for the Implicit Association Test (IAT)
#'
#' @description Log ratio analysis for the IAT.
#' @export
#' @param rt A vector specifying all the reaction times.
#' @param subject A vector specifying the subject IDs for the \code{rt} vector.
#' @param block_type A vector specifying the block type for the \code{rt} vector. There should be two blocks.
#' @param trial_num A vector specifying the trial number for each observation in the \code{rt} vector.
#' @param group A data frame with two columns specifying the subject IDs and corresponding group. There should be two groups.
#' @param block_order A character string specifying the order of the two groups.
#' There are two options: "original" puts the first group in the numerator of the ratio, and
#' "reverse" puts the second group in the numerator of the ratio.
#' @param min_limit A numeric specifying the lower limit for the reaction times to be included. The default option is 400.
#' @param max_limit A numeric specifying the upper limit for the reaction times to be included. The default option is 10000.
#' @param rt_min A numeric specifying the minimum time required for reaction. The default option is 200.
#' @param correctvec A vector specifying whether or not the response is correct (0 for incorrect, 1 for correct). 
#' The default option is NULL, in which case, all the responses are assumed to be correct.
#' @param trace A boolean specifying whether or not the progress should be displayed on the screen. The default option is FALSE.
#'
#' @return 
#' \item{scores}{A list containing the IAT scores using the log ratio analysis.}
#'
#' @examples
#' data(reactiontimes)
#' data(grouping)
#' rt <- reactiontimes$rt
#' subject <- reactiontimes$subId
#' block_type <- reactiontimes$block_type
#' trial_num <- reactiontimes$trial_num
#' block_order <- "reverse"
#' \donttest{
#' results <- logratio(rt=rt, subject=subject, block_type=block_type, trial_num=trial_num,
#' group=grouping, block_order=block_order, trace=TRUE)
#' femaleRatioLog<-results$`0`
#' maleRatioLog<-results$`1`
#' two.sample.var(femaleRatioLog,maleRatioLog,alternative="two.sided",
#' scale.option="Levene.Med.0",scale.adj=TRUE,paired=FALSE)
#' Bonett.Seier.test(femaleRatioLog,maleRatioLog,alternative="two.sided")
#'}
#'
#' @importFrom stats median

logratio <- function(rt, subject, block_type, trial_num, group, block_order=c("original","reverse"), 
min_limit=400, max_limit=10000, rt_min=200, correctvec=NULL, trace=FALSE)
{
 total <- length(rt)
 block_order <- match.arg(block_order)
 
 if(length(subject)!=total)
 {
  stop("Length of subject vector is not equal to length of rt")
 }
 if(length(block_type)!=total)
 {
  stop("Length of block_type vector is not equal to length of rt")
 }
 if(length(trial_num)!=total)
 {
  stop("Length of trial_num vector is not equal to length of rt")
 }
 sub_ids <- unique(subject)
 if(length(sub_ids)!=dim(group)[1])
 {
  stop("Length of unique subject IDs is not equal to length of group")
 } 
 
 if(is.null(correctvec))
 {
  correctvec <- rep(1,total)
 } 
 
 block_type01 <- unique(block_type)
 if(length(block_type01) != 2)
 {
  stop("Length of unique groups is not equal to 2") 	
 } 
 numSubs <- length(sub_ids)
 group01 <- unique(group[,2])
 if(length(group01) != 2)
 {
  stop("Length of unique groups is not equal to 2") 	
 }
 else
 {
  if(block_order=="original")
  {
   firstgroup <- group01[1]
   secondgroup <- group01[2]
  }
  else #block_order=="reverse"
  {
   firstgroup <- group01[2]
   secondgroup <- group01[1]
  }  
  message("The order of grouping is: ", firstgroup, " (numerator), ",
  secondgroup, " (denominator).",sep="")
 }
 num0 <- length(which(group[,2]==group01[1]))
 num1 <- length(which(group[,2]==group01[2]))
 subject0 <- group[which(group[,2]==group01[1]),1]
 subject1 <- group[which(group[,2]==group01[2]),1]
 
 allRatio<-rep(NA, numSubs)
 i_all<-1
 
 ### gender IAT
 Ratio_0<-rep(NA, num0)
 i_0<-1
 Ratio_1<-rep(NA, num1)
 i_1<-1
 
 # trial numbers to analyze
 totTrials<-sort(unique(trial_num))
 
 trialRatio <- rep(NA, length(totTrials))
 i_trial<-1

 # generate IAT scores
 subNum<-1 
 
 for(id in sub_ids)
 {
  if(trace==TRUE)
  {
   message(paste("Processing", subNum, "of", numSubs, "subjects"))
  }
  
  i_trial <- 1
  for(trialNum in totTrials)
  {
   sRow<-which((subject==id)*(trial_num==trialNum)*(block_type==block_type01[1])==1)
   gRow<-which((subject==id)*(trial_num==trialNum)*(block_type==block_type01[2])==1)
   
   includes<-(rt[sRow]>=min_limit)*(rt[sRow]<=max_limit)
   includeg<-(rt[gRow]>= min_limit)*(rt[gRow]<= max_limit)
   include<-includes*includeg*correctvec[sRow]*correctvec[gRow]
   if(include==1)
   {
   	if(block_order=="original")
   	{
     trialRatio[i_trial] <- (rt[sRow] - rt_min)/(rt[gRow] - rt_min)   		
   	}
   	else #block_order == "reverse"
   	{
     trialRatio[i_trial] <- (rt[gRow] - rt_min)/(rt[sRow] - rt_min)    		
   	}
   }
   else
   {
    trialRatio[i_trial] <- NA
   }
   i_trial <- i_trial + 1
  }
  
  # take median of all trials for subject
  ratioMedian <- median(trialRatio,na.rm=TRUE)
  
  # add to total group
  allRatio[i_all] <- ratioMedian
  i_all <- i_all+1
  
  # check if group0
  if(is.element(id, subject0))
  {
   Ratio_0[i_0] <- ratioMedian
   i_0 <- i_0 + 1
  }
  
  # check if group1
  if(is.element(id, subject1))
  {
   Ratio_1[i_1] <- ratioMedian
   i_1 <- i_1 + 1
  }
  subNum <- subNum + 1
 }
 if(trace==TRUE)
 {
  message(paste("Done"))
 } 
 
 # take log of the ratios and analyze
 LogRatio_0 <- log(Ratio_0)
 LogRatio_1 <- log(Ratio_1)
 
 scores <- list(LogRatio_0, LogRatio_1)
 names(scores) <- as.character(group01)
 
 return(scores)
} 
