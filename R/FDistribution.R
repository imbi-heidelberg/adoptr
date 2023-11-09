#' F-Distribution
#'
#' Implements the F-distribution used for an ANOVA or for the comparison of the fit of two
#' nested regression models. In both cases, the test statistic follows a F-distribution.
#' \code{NestedModel} is used to compare the fit of two regression models, where one model contains
#' the independent variables of the smaller model as a subset. Then, one can use ANOVA to determine
#' whether more variance can be explained by adding more independent variables.
#' In the class \code{ANOVA}, the number of independent variables of the smaller model is set to \eqn{1}
#' in order to match the degrees of freedom and we obtain a one-way ANOVA.
#'
#' @slot p_inner number of parameters in smaller model
#' @slot p_outer number of parameters in bigger model
#'
#' @template DataDistributionTemplate
#'
#' @include DataDistribution.R
#'
#' @rdname NestedModels-class
#' @exportClass NestedModels
setClass("NestedModels", representation(
    p_inner  = "numeric",
    p_outer =  "numeric"),
    contains = "DataDistribution")


#' @include DataDistribution.R
#'
#' @rdname ANOVA-class
#' @exportClass ANOVA
setClass("ANOVA", contains = "NestedModels")



#' @param p_inner number of independent variables in smaller model
#' @param p_outer number of independent variables in bigger model
#'
#' @examples
#' model <- NestedModels(2, 4)
#'
#' @seealso see \code{\link{probability_density_function}} and
#'    \code{\link{cumulative_distribution_function}} to evaluate the pdf
#'    and the cdf, respectively.
#'
#' @rdname NestedModels-class
#' @export
NestedModels <- function(p_inner, p_outer) {
    if (any(p_inner < 0, p_outer < 0))
        stop("Number of parameters may not be less than 0.")
    new("NestedModels", p_inner = p_inner, p_outer = p_outer)
}


#' Analysis of Variance
#'
#' ANOVA is used to test whether there is a significant difference between the means of groups.
#' The sample size which \code{adoptr} returns is the total sample size,
#'
#' @param n_groups number of groups to be compared
#'
#' @examples
#' model <- ANOVA(3L)
#'
#' @seealso see \code{\link{probability_density_function}} and
#'    \code{\link{cumulative_distribution_function}} to evaluate the pdf
#'    and the cdf, respectively.
#'
#' @rdname ANOVA-class
#' @export
ANOVA <- function(n_groups) {
    new("ANOVA", p_outer = n_groups, p_inner = 1L)
}


#' @export
get_tau_ANOVA <- function(means, common_sd) {
    mean((means - mean(means))^2) / common_sd^2
}


#' @examples
#' probability_density_function(Binomial(.2, FALSE), 1, 50, .3)
#'
#' @details If the distribution is \code{\link{Binomial}},
#'   \ifelse{html}{\out{theta}}{\eqn{theta}} denotes the rate difference between
#'   intervention and control group.
#'   Then, the mean is assumed to be
#'   \ifelse{html}{\out{&radic; n  theta}}{\eqn{\sqrt{n} theta}}.
#'
#' @rdname probability_density_function
#' @export
setMethod("probability_density_function", signature("NestedModels", "numeric", "numeric", "numeric"),
          function(dist, x, n, theta, ...) {
              return(stats::df(x, df1 = dist@p_outer - dist@p_inner, df2 = n - dist@p_outer,  ncp = n * theta))
          })


#' @examples
#' cumulative_distribution_function(Binomial(.1, TRUE), 1, 50, .3)
#'
#' @details If the distribution is \code{\link{Binomial}},
#'   \ifelse{html}{\out{theta}}{\eqn{theta}} denotes the rate difference between
#'   intervention and control group.
#'   Then, the mean is assumed to be
#'   \ifelse{html}{\out{&radic; n  theta}}{\eqn{\sqrt{n} theta}}.
#'
#' @rdname cumulative_distribution_function
#' @export
setMethod("cumulative_distribution_function", signature("NestedModels", "numeric", "numeric", "numeric"),
          function(dist, x, n, theta, ...) {
              ret <- x
              ret[x==Inf] <- 1
              ret[x==-Inf] <- 0
              ret[!is.infinite(x)] <- stats::pf(x, df1 = dist@p_outer - dist@p_inner, df2 = n - dist@p_outer,  ncp = n * theta)
              #cat('df1: ', dist@p_outer - dist@p_inner, ' df2: ', n - dist@p_outer, ' ncp: ', n * theta, ' n: ', n, '\n')
              return(ret)
        })


#' @param probs vector of probabilities
#' @rdname BinomialDataDistribution-class
#' @export
setMethod("quantile", signature("NestedModels"),
          function(x, probs, n, theta, ...) { # must be x to conform with generic
              return(stats::qf(probs, df1 = x@p_outer - x@p_inner, df2 = n - x@p_outer,  ncp = n * theta))
          })



#' @details Note that \code{simulate} for class \code{Binomial} simulates the
#'    normal approximation of the test statistic.
#'
#' @rdname BinomialDataDistribution-class
#'
#' @param object object of class \code{Binomial}
#' @param nsim number of simulation runs
#' @param seed random seed
#'
#' @export
setMethod("simulate", signature("NestedModels", "numeric"),
          function(object, nsim, n, theta, seed = NULL, ...) {
              if (!is.null(seed)) set.seed(seed)
              return(stats::rf(nsim, df1 = x@p_outer - x@p_inner, df2 = n - x@p_outer,  ncp = n * theta))
})

setMethod("print", signature('NestedModels'), function(x, ...) {
    glue::glue(
        "{class(x)[1]}<p_inner={x@p_inner}, p_outer={x@p_outer}>"
    )
})

setMethod("print", signature('ANOVA'), function(x, ...) {
    glue::glue(
        "{class(x)[1]}<n_groups={x@p_inner}>"
    )
})









