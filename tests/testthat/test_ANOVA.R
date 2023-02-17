# ANOVA with 2 groups is basically a two-sided t-test
H_0 <- PointMassPrior(0, 1)
H_1 <- PointMassPrior(get_tau_ANOVA(c(0, 0.4), 1), 1)

datadist <- ANOVA(2)
ess <- ExpectedSampleSize(datadist, H_1)
power <- Power(datadist, H_1)
toer  <- Power(datadist, H_0)
initial_design <- get_initial_design(
    theta = sqrt(H_1@theta),
    alpha = .025,
    beta  = .2,
    type_design  = "two-stage",
    dist  = Student(TRUE),
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
