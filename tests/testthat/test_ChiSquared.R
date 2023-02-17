# Two-sided, one group test (Z^2 ~ Chi-Squared)
H_0 <- PointMassPrior(0, 1)
H_1 <- PointMassPrior(get_tau_ZSquared(.4), 1)
datadist <- ZSquared(FALSE)
ess <- ExpectedSampleSize(datadist, H_1)
power <- Power(datadist, H_1)
toer  <- Power(datadist, H_0)
initial_design <- get_initial_design(
    theta = sqrt(H_1@theta),
    alpha = .025,
    beta  = .2,
    type_design  = "two-stage",
    dist  = Normal(FALSE),
    order = 7L
)
initial_design@c1f <- .2
initial_design@c1e <- initial_design@c1e^2*.8
initial_design@c2_pivots <- initial_design@c2_pivots^2*.8
evaluate(ess, initial_design)
evaluate(power, initial_design)
evaluate(toer, initial_design)

opt_res <- minimize(
    ess,
    subject_to(
        power >= 0.8,
        toer  <= .05
    ),
    initial_design,
    opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-06, maxeval = 10000)
)
opt_d <- opt_res$design
plot(opt_d)

evaluate(ess, opt_d)
evaluate(power, opt_d)
evaluate(toer, opt_d)




## Normal, for reference
H_0 <- PointMassPrior(0, 1)
H_1 <- PointMassPrior(.4, 1)
datadist <- Normal(FALSE)
ess <- ExpectedSampleSize(datadist, H_1)
power <- Power(datadist, H_1)
toer  <- Power(datadist, H_0)
initial_design <- get_initial_design(
    theta = H_1@theta,
    alpha = .025,
    beta  = .2,
    type_design  = "two-stage",
    dist  = datadist,
    order = 7L
)
evaluate(ess, initial_design)
evaluate(power, initial_design)
evaluate(toer, initial_design)

opt_res <- minimize(
    ess,
    subject_to(
        power >= 0.8,
        toer  <= .025
    ),
    initial_design
)
opt_d2 <- opt_res$design
plot(opt_d2)

evaluate(ess, opt_d2)
evaluate(power, opt_d2)
evaluate(toer, opt_d2)


## Some reality checks
calc_ncp <- function(p_vector){
    n_groups <- length(p_vector)
    mean_p <- mean(p_vector)
    deltas <- p_vector - mean_p
    (sum(deltas^2)/n_groups - mean(deltas)^2) / (mean_p * (1-mean_p))
}

ni <- 1000
nsim <- 1e5
p1 <- .25
p2 <- .45
p3 <- .3
Chis <- replicate(nsim,
                  {
                      X1 <- rbinom(1L, ni, p1)
                      X2 <- rbinom(1L, ni, p2)
                      X3 <- rbinom(1L, ni, p3)
                      tab <-
                          matrix(c(X1, X2, X3, ni - X1, ni - X2, ni - X3), byrow = TRUE, ncol = 3)
                      chisq.test(tab)$statistic
                  })

plot(density(Chis))
ref_x <- seq(30, 170, .01)
ref_y <- dchisq(ref_x, df = (2L-1L)*(3L-1L), ncp =  3*ni*calc_ncp(c(p1, p2, p3)) )
lines(ref_x, ref_y)









