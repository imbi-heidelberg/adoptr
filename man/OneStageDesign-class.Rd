% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/OneStageDesign.R
\docType{class}
\name{OneStageDesign-class}
\alias{OneStageDesign-class}
\alias{OneStageDesign}
\alias{OneStageDesign,numeric-method}
\alias{TwoStageDesign,OneStageDesign-method}
\alias{TwoStageDesign,OneStageDesignSurvival-method}
\alias{plot,OneStageDesign-method}
\title{One-stage designs}
\usage{
OneStageDesign(n, ...)

\S4method{OneStageDesign}{numeric}(n, c, event_rate)

\S4method{TwoStageDesign}{OneStageDesign}(n1, event_rate, order = 5L, eps = 0.01, ...)

\S4method{TwoStageDesign}{OneStageDesignSurvival}(n1, order = 5L, eps = 0.01, ...)

\S4method{plot}{OneStageDesign}(x, y, ...)
}
\arguments{
\item{n}{sample size (stage-one sample size)}

\item{...}{further optional arguments}

\item{c}{rejection boundary (\ifelse{html}{\out{c = c<sub>1</sub><sup>f</sup> = c<sub>1</sub><sup>e</sup>}}{\eqn{c = c_1^f = c_1^e}})}

\item{event_rate}{probability that a subject in either group will eventually have an event,
only needs to be specified for time-to-event endpoints.}

\item{n1}{\code{OneStageDesign} object to convert, overloaded from
\code{\link{TwoStageDesign}}}

\item{order}{integer >= 2, default is 5; order of Gaussian quadrature
integration rule to use for new TwoStageDesign.}

\item{eps}{numeric > 0, default = .01; the single critical value c must be
split in a continuation interval [c1f, c1e]; this is given by c +/- eps.}

\item{x}{design to plot}

\item{y}{not used}
}
\description{
\code{OneStageDesign} implements a one-stage design as special case of
a two-stage design, i.e. as sub-class of \code{\link{TwoStageDesign}}.
This is possible by defining \ifelse{html}{\out{n<sub>2</sub> = 0}}{\eqn{n_2=0}},
\ifelse{html}{\out{c = c<sub>1</sub><sup>f</sup> = c<sub>1</sub><sup>e</sup>}}{\eqn{c = c_1^f = c_1^e}},
\ifelse{html}{\out{c<sub>2</sub>(x<sub>1</sub>) = ifelse(x<sub>1</sub> < c, Inf, -Inf)}}{\eqn{c_2(x_1) = ifelse(x_1 < c, Inf, -Inf)}}.
No integration pivots etc are required (set to \code{NaN}).
}
\details{
Note that the default \code{\link{plot,TwoStageDesign-method}} method
is not supported for \code{OneStageDesign} objects.
}
\examples{
design <- OneStageDesign(30, 1.96)
summary(design)
design_twostage <- TwoStageDesign(design)
summary(design_twostage)
design_survival <- OneStageDesign(30,1.96,0.7)

TwoStageDesign(design_survival)
}
\seealso{
\code{\link{TwoStageDesign}}, \code{\link{GroupSequentialDesign-class}}
}
