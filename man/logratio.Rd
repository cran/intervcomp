% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logratio.R
\name{logratio}
\alias{logratio}
\title{Log Ratio Analysis for the Implicit Association Test (IAT)}
\usage{
logratio(rt, subject, block_type, trial_num, group,
  block_order = c("original", "reverse"), min_limit = 400,
  max_limit = 10000, rt_min = 200, correctvec = NULL,
  trace = FALSE)
}
\arguments{
\item{rt}{A vector specifying all the reaction times.}

\item{subject}{A vector specifying the subject IDs for the \code{rt} vector.}

\item{block_type}{A vector specifying the block type for the \code{rt} vector. There should be two blocks.}

\item{trial_num}{A vector specifying the trial number for each observation in the \code{rt} vector.}

\item{group}{A data frame with two columns specifying the subject IDs and corresponding group. There should be two groups.}

\item{block_order}{A character string specifying the order of the two groups.
There are two options: "original" puts the first group in the numerator of the ratio, and
"reverse" puts the second group in the numerator of the ratio.}

\item{min_limit}{A numeric specifying the lower limit for the reaction times to be included. The default option is 400.}

\item{max_limit}{A numeric specifying the upper limit for the reaction times to be included. The default option is 10000.}

\item{rt_min}{A numeric specifying the minimum time required for reaction. The default option is 200.}

\item{correctvec}{A vector specifying whether or not the response is correct (0 for incorrect, 1 for correct). 
The default option is NULL, in which case, all the responses are assumed to be correct.}

\item{trace}{A boolean specifying whether or not the progress should be displayed on the screen. The default option is FALSE.}
}
\value{
\item{scores}{A list containing the IAT scores using the log ratio analysis.}
}
\description{
Log ratio analysis for the IAT.
}
\examples{
data(reactiontimes)
data(grouping)
rt <- reactiontimes$rt
subject <- reactiontimes$subId
block_type <- reactiontimes$block_type
trial_num <- reactiontimes$trial_num
block_order <- "reverse"
\donttest{
results <- logratio(rt=rt, subject=subject, block_type=block_type, trial_num=trial_num,
group=grouping, block_order=block_order, trace=TRUE)
femaleRatioLog<-results$`0`
maleRatioLog<-results$`1`
two.sample.var(femaleRatioLog,maleRatioLog,alternative="two.sided",
scale.option="Levene.Med.0",scale.adj=TRUE,paired=FALSE)
Bonett.Seier.test(femaleRatioLog,maleRatioLog,alternative="two.sided")
}

}
