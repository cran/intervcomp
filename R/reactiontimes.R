#' @title Reaction Time (RT) Data for the Implicit Association Test
#'
#' @name reactiontimes
#'
#' @description Each subject recorded one hundred RT observations for each of the
#' two pairings (congruent vs. incongrunet). That is, in total, each subject has
#' two hundred RT observations in total. The variables are as follows:
#'
#' \itemize{
#'   \item subId: Subject IDs
#'   \item stim_pairing: Pairings of the stimulus.
#'   \item block_type: Block type (gaygd or strgd).
#'   \item trial_num: Trial number.
#'   \item rt: Reaction time.
#'   \item correct: Whether or not the response was correct (1=Correct, 0=Incorrect).
#' }
#'
#' @docType data
#'
#' @usage data(reactiontimes)
#'
#' @format A data frame with 16600 rows and 6 columns.
#'
#' @keywords datasets
#'
#' @references Lemm, K.M. (2006).
#' Positive associations among interpersonal contact, motivation, and implicit
#' and explicit attitudes toward gay men. Journal of Homosexuality 51:2, 79-99.
#'
#' @examples
#' data(reactiontimes)
NULL
