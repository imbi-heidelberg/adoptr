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
#' @slot df1 cf. parameter 'df'
#' @slot df2 cf. parameter 'df'
#'
#' @template DataDistributionTemplate
#'
#' @include DataDistribution.R
#'
#' @rdname FDataDistribution-class
#' @exportClass FDistribution
setClass("FDistribution", representation(
    df1  = "numeric",
    df2 =  "numeric"),
    contains = "DataDistribution")


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
#' @rdname FDataDistribution-class
#' @export
FDistribution <- function(df1, df2) {
    if (any(df1 < 0, df2 < 0))
        stop("The degrees of freedom may not be less than 0.")
    new("FDistribution", df1 = df1, df2 = df2)
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
#' @rdname ChiSquaredlDataDistribution-class
#' @export
Pearson2xK <- function(n_groups) {
    new("ChiSquared", df = n_groups-1L, multiplier = 1/n_groups)
}

#' @example
#' H1 <- PointMassPrior(get_ncp_Pearson2xK(c(.3, .25, .4)), 1)
get_ncp_Pearson2xK <- function(p_vector) {
    n_groups <- length(p_vector)
    mean_p <- mean(p_vector)
    deltas <- p_vector - mean_p
    tau <- (sum(deltas^2)/n_groups - mean(deltas)^2) / (mean_p * (1-mean_p))
    tau
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
setMethod("probability_density_function", signature("ChiSquared", "numeric", "numeric", "numeric"),
          function(dist, x, n, theta, ...) {
              return(stats::dchisq(x, df = dist@df, ncp = n / dist@multiplier * theta))
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
setMethod("cumulative_distribution_function", signature("ChiSquared", "numeric", "numeric", "numeric"),
          function(dist, x, n, theta, ...) {
              return(stats::pchisq(x, df = dist@df, ncp = n / dist@multiplier * theta))
        })


#' @param probs vector of probabilities
#' @rdname BinomialDataDistribution-class
#' @export
setMethod("quantile", signature("ChiSquared"),
          function(x, probs, n, theta, ...) { # must be x to conform with generic
              return(stats::qchisq(probs, df = dist@df, ncp = n / dist@multiplier * theta))
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
setMethod("simulate", signature("ChiSquared", "numeric"),
          function(object, nsim, n, theta, seed = NULL, ...) {
              if (!is.null(seed)) set.seed(seed)
              return(stats::rchisq(nsim, df = dist@df, ncp = n / dist@multiplier * theta))
})

setMethod("print", signature('ChiSquared'), function(x, ...) {
    glue::glue(
        "{class(x)[1]}<df={x@df}>"
    )
})



