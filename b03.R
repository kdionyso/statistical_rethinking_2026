install.packages(c("coda", "mvtnorm", "devtools", "loo", "dagitty", "shape"))
devtools::install_github("rmcelreath/rethinking", force = TRUE)
library(rethinking)
data(bangladesh)
d <- bangladesh
str(d)

dat <- list(
    C = d$use.contraception,
    D = as.integer(d$district),
    U = d$urban #,
    #    A = standardize(d$age.centered),
    #    K = d$living.children
)

mCDU = ulam(
    alist(
        C ~ bernoulli(p),
        logit(p) <- a[D] + b[D] * U,
        transpars > vector[61]:a <<- abar[1] + v[, 1],
        transpars > vector[61]:b <<- abar[2] + v[, 2],
        transpars > matrix[61, 2]:v <- compose_noncentered(sigma, L_Rho, Z),
        matrix[2, 61]:Z ~ normal(0, 1),
        vector[2]:abar ~ normal(0, 1),
        vector[2]:sigma ~ exponential(1),
        cholesky_factor_corr[2]:L_Rho ~ lkj_corr_cholesky(4),
        # Correlation
        gq > matrix[2, 2]:Rho <<- Chol_to_Corr(L_Rho)
    ),
    data = dat,
    cores = 4,
    chains = 4,
    iter = 2000,
    log_lik = TRUE,
    cmdstan = TRUE
)

mCDUcov_nc <- ulam(
    alist(
        C ~ bernoulli(p),
        logit(p) <- a[D] + b[D] * U,
        # define effects using other parameters
        # this is the non-centered Cholesky machine
        transpars > vector[61]:a <<- abar[1] + v[, 1],
        transpars > vector[61]:b <<- abar[2] + v[, 2],
        transpars > matrix[61, 2]:v <-
            compose_noncentered(sigma, L_Rho, Z),
        # priors - note that none have parameters inside them
        # that is what makes them non-centered
        matrix[2, 61]:Z ~ normal(0, 1),
        vector[2]:abar ~ normal(0, 1),
        cholesky_factor_corr[2]:L_Rho ~ lkj_corr_cholesky(4),
        vector[2]:sigma ~ exponential(1),
        # convert Cholesky to Corr matrix
        gq > matrix[2, 2]:Rho <<- Chol_to_Corr(L_Rho)
    ),
    data = dat,
    chains = 4,
    cores = 4
)
t <- precis(mCDU, depth = 2)

# covariance - non-centered
post <- extract.samples(mCDU)
par(mfrow = c(1, 1))
s <- mCDU@cstanfit$summary()
# m <- cbind(as.data.frame(s)$rhat, as.data.frame(s)$ess_bulk)
plot(
    s$ess_bulk,
    s$rhat,
    ylab = "Rhat",
    xlab = "Effective Sample Size",
    main = "Rhat vs Effective Sample Size"
)
abline(h = 1.01, col = "red", lty = 2)
abline(v = 200, col = "blue", lty = 2)

par(mfrow = c(1, 1))

# Show HMC energy diagnostics
diagnostics <- mCDU@cstanfit$sampler_diagnostics(inc_warmup = FALSE)
hmc_energy <- diagnostics[,, grepl("energy", dimnames(diagnostics)$variable)] # mCDU@cstanfit$sampler_diagnostics()
min_max <- c(min(hmc_energy), max(hmc_energy))
dens_min_max <- c(
    min(density(as.matrix(hmc_energy))$y),
    max(density(as.matrix(hmc_energy))$y)
)
#plot(NULL, xlim=min_max, ylim=dens_min_max, ylab="density", xlab="Energy")

plot(density(hmc_energy), lwd = 2)
dens(hmc_energy, add = TRUE, col = 2)

install.packages("ggplot2")
install.packages("bayesplot")
library(ggplot2)
library(bayesplot)
color_scheme_set("red")
np <- nuts_params(mCDU@cstanfit)
mcmc_nuts_energy(np) + ggtitle("NUTS Energy Diagnostic")
# m_diagnostics <-as.matrix(diagnostics[grepl("energy",dimnames(diagnostics)$variable)])
# dens(
#     m_diagnostics,
# #    main = "HMC energy"
#     xlim=c(1200,1400),
#     ylim=c(0,0.007)
# )

dashboard(mCDU)


# posterior rho
par(mfrow = c(1, 2))
post <- extract.samples(mCDU)
dens(
    post$Rho[, 1, 2],
    xlim = c(-1, 1),
    lwd = 3,
    col = 2,
    xlab = "posterior correlation a,b"
)
abline(v = 0, lty = 2, lwd = 0.5)
prior_rho <- rlkjcorr(1e4, 2, eta = 4)
dens(prior_rho[, 1, 2], lwd = 2, lty = 2, add = TRUE)


# posterior MVN of a,b
install.packages("ellipse")
library(ellipse)
plot(NULL, xlim = c(-2, 1), ylim = c(-1, 2), xlab = "a", ylab = "b")
abline(v = 0, lty = 2, lwd = 0.5)
abline(h = 0, lty = 2, lwd = 0.5)
SIGMA <- cov(cbind(apply(post$a, 2, mean), apply(post$b, 2, mean)))
MU <- apply(post$abar, 2, mean)
for (l in seq(from = 0.25, to = 0.95, len = 5)) {
    el <- ellipse(SIGMA, centre = MU, level = l)
    #lines( (el) , col=2 , lwd=3 )
    polygon((el), col = col.alpha(2, 0.25), border = NA)
}

Y <- rmvnorm(6, c(0, 0), sigma = SIGMA)
points(Y, lwd = 4, col = "white")
points(Y, lwd = 3, col = 2)

plot(NULL, xlim = c(-2.5, 2.5), ylim = c(-2.5, 2.5), xlab = "a", ylab = "b")
for (i in 1:10) {
    RHO <- rlkjcorr(1, 2, eta = 4)
    s <- rexp(1, 1)
    tau <- rexp(1, 1)
    SIGMA <- diag(c(s, tau)) %*% RHO %*% diag(c(s, tau))
    el <- ellipse(SIGMA, centre = MU, level = 0.89)
    lines((el), col = col.alpha(2, 0.5), lwd = 3)
    #polygon( (el) , col=col.alpha(2,0.25) , border=NA )
}

SIGMA <- rlkjcorr(1e4, 2, eta = 4)
dens(SIGMA[, 1, 2], lwd = 3, col = 2, xlab = "correlation")


# Probability of urban vs rural

library(httpgd)
library(languageserver)
library(vscDebugger)

post <- extract.samples(mCDUcov_nc)
mu_rural <- inv_logit(post$a)
mu_urban <- inv_logit(post$b + post$a)
par(mfrow = c(1, 1))
plot(
    NULL,
    xlim = c(0.1, 0.6),
    ylim = c(0.2, 0.7),
    xlab = "Rural",
    ylab = "Urban",
    cex = 1
)
points(
    apply(mu_rural, 2, mean),
    apply(mu_urban, 2, mean),
    lwd = 2,
    main = "Probability of Contraception Use: Urban vs Rural"
)
abline(h = 0.5, lty = 2, col = "gray", lw = 2)
abline(v = 0.5, lty = 2, col = "gray", lw = 2)
