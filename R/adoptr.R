#' Adaptive Optimal Two-Stage Designs
#'
#' The \pkg{adoptr} package provides functionality to explore custom optimal
#' two-stage designs for one- or two-arm superiority tests. More than two arms
#' can be compared via chi-squared tests or ANOVA.
#' For more details on the theoretical background see
#' <doi:10.1002/sim.8291> and <doi:10.18637/jss.v098.i09>.
#' \pkg{adoptr} makes heavy use of the S4 class system.
#' A good place to start learning about it can be found
#' \href{http://adv-r.had.co.nz/OO-essentials.html}{here}.
#'
#'
#'
#' @section Quickstart:
#'
#' For a sample workflow and a quick demo of the capabilities, see
#' \href{https://imbi-heidelberg.github.io/adoptr/articles/adoptr.html}{here}.
#'
#' A more detailed description of the background and the usage of \pkg{adoptr}
#' can be found \href{https://imbi-heidelberg.github.io/adoptr/articles/adoptr_jss.html}{here}
#' or here <doi:10.18637/jss.v098.i09> .
#'
#' A variety of examples is presented in the validation report hosted
#' \href{https://imbi-heidelberg.github.io/adoptr-validation-report/}{here}.
#'
#'
#'
#' @section Designs:
#'
#' \pkg{adoptr} currently supports \code{\link{TwoStageDesign}},
#' \code{\link{GroupSequentialDesign-class}}, and \code{\link{OneStageDesign-class}}.
#'
#'
#'
#' @section Data distributions:
#'
#' The implemented data distributions are \code{\link{Normal}}, \code{\link{Binomial}},
#' \code{\link{Student}}, \code{\link{Survival}}, \code{\link{ChiSquared}} (including
#' \code{\link{Pearson2xK}} and \code{\link{ZSquared}}) and \code{\link{ANOVA}}.
#'
#'
#'
#'
#' @section Priors:
#'
#' Both \code{\link{ContinuousPrior}} and \code{\link{PointMassPrior}} are
#' supported for the single parameter of a \code{\link{DataDistribution}}.
#'
#'
#'
#' @section Scores:
#'
#' See \code{\link{Scores}} for information on the basic system of representing
#' scores.
#' Available scores are \code{\link{ConditionalPower}},
#' \code{\link{ConditionalSampleSize}}, \code{\link{Power}}, and
#' \code{\link{ExpectedSampleSize}}..
#'
#' @import methods
#' @docType package
#' @name adoptr
NULL
