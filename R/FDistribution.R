#' Chi-Squared data distribution
#'
#' Viel Spaß, Nico.
#' Implements the normal approximation for a test on rates.
#' The reponse rate in the control group,
#' \ifelse{html}{\out{r<sub>C</sub>}}{\eqn{r_C}}, has to be specified by
#' \code{rate_control}.
#' The null hypothesis is:
#' \ifelse{html}{\out{r<sub>E</sub> &le; r<sub>C</sub>}}{\eqn{r_E <= r_C}},
#' where \ifelse{html}{\out{r<sub>E</sub>}}{\eqn{r_E}} denotes the response rate
#' in the invervention group.
#' It is tested against the alternative
#' \ifelse{html}{\out{r<sub>E</sub> > r<sub>C</sub>}}{\eqn{r_E > r_C}}.
#' The test statistic is given as
#' \ifelse{html}{\out{X<sub>1</sub> = &radic;n (r<sub>E</sub> - r<sub>C</sub>) / &radic;(2  r<sub>0</sub> (1-r<sub>0</sub>))}}{\eqn{X_1 = \sqrt{n}(r_E - r_C) / \sqrt{2 r_0 (1- r_0)}}},
#' where \ifelse{html}{\out{r<sub>0</sub>}}{\eqn{r_0}} denotes the mean between
#' \ifelse{html}{\out{r<sub>E</sub>}}{\eqn{r_E}} and
#' \ifelse{html}{\out{r<sub>C</sub>}}{\eqn{r_C}} in the two-armed case,
#' and \ifelse{html}{\out{r<sub>E</sub>}}{\eqn{r_E}} in the one-armed case.#'
#' All priors have to be defined for the rate difference
#' \ifelse{html}{\out{r<sub>E</sub> - r<sub>C</sub>}}{\eqn{r_E - r_C}}.
#'
#' @slot p1 cf. parameter 'df'
#' @slot p2 cf. parameter 'df'
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


setClass("ANOVA", contains = "NestedModels")




#' @param rate_control assumed response rate in control group
#' @param two_armed logical indicating if a two-armed trial is regarded
#'
#' @examples
#' datadist <- Binomial(rate_control = 0.2, two_armed = FALSE)
#'
#' @seealso see \code{\link{probability_density_function}} and
#'    \code{\link{cumulative_distribution_function}} to evaluate the pdf
#'    and the cdf, respectively.
#'
#' @rdname NestedModels-class
#' @export
NestedModels <- function(p_inner, p_outer) {
    if (any(p_inner < 0, p_outer < 0))
        stop("Numbe of parameters may not be less than 0.")
    new("NestedModels", p_inner = p_inner, p_outer = p_outer)
}


#' @param rate_control assumed response rate in control group
#' @param two_armed logical indicating if a two-armed trial is regarded
#'
#' @examples
#' datadist <- Binomial(rate_control = 0.2, two_armed = FALSE)
#'
#' @seealso see \code{\link{probability_density_function}} and
#'    \code{\link{cumulative_distribution_function}} to evaluate the pdf
#'    and the cdf, respectively.
#'
#' @rdname ChiSquaredDataDistribution-class
#' @export
ANOVA <- function(n_groups) {
    new("ANOVA", p_inner = n_groups, p_outer = 0L)
}



get_tau_ANOVA <- function(means, common_sd) {
    mean((means - mean(means))^2)/common_sd^2
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
              return(stats::df(x, df1 = dist@p_inner - dist@p_outer, df2 = n-dist@p_inner,  ncp = n * theta))
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
              ret[x==-Inf] <- 1
              ret[x==-Inf] <- 0
              ret[!is.infinite(x)] <- stats::pf(x, df1 = dist@p_inner - dist@p_outer, df2 = n-dist@p_inner,  ncp = n * theta)
              return(ret)
        })


#' @param probs vector of probabilities
#' @rdname BinomialDataDistribution-class
#' @export
setMethod("quantile", signature("NestedModels"),
          function(x, probs, n, theta, ...) { # must be x to conform with generic
              return(stats::qf(probs, df1 = x@p_inner - x@p_outer, df2 = n-x@p_inner,  ncp = n * theta))
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
              return(stats::rf(nsim, df1 = x@p_inner - x@p_outer, df2 = n-x@p_inner,  ncp = n * theta))
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









