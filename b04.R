install.packages(c("coda", "mvtnorm", "devtools", "loo", "dagitty", "shape"))
devtools::install_github("rmcelreath/rethinking", force = TRUE)
library(rethinking)

str(d)
install.packages("ggplot2")
install.packages("bayesplot")
library(ggplot2)
library(bayesplot)


# Fake data

library(rethinking)

set.seed(8672)

sim_dat <- function(
    N_groups = 30,
    N_id = 200,
    a0 = (-2),
    bZY = (-1),
    bXY = 0
) {
    g <- sample(1:N_groups, size = N_id, replace = TRUE) # sample into groups
    Ug <- rnorm(N_groups, 1.5) # group confounds
    X <- rnorm(N_id, Ug[g]) # individual varying trait
    Z <- rnorm(N_groups) # group varying trait (observed)
    Y <- rbern(N_id, p = inv_logit(a0 + bXY * X + Ug[g] + bZY * Z[g]))
    dat <- list(
        g = g,
        Ug = Ug,
        X = X,
        Z = Z,
        Y = Y,
        Ng = N_groups,
        N_id = N_id,
        a0 = a0,
        bZY = bZY,
        bXY = bXY
    )
    return(dat)
}

dat <- sim_dat(N_groups = 30, N_id = 2000, bZY = 1, bXY = 0)

# dat$groups <-
# Plot
plot(dat$g, dat$X, xlab = "group", ylab = "X", lwd = 2, col = 2, cex = 3)
# Fixed effects model

datN <- list(
    X = dat$X,
    Z = dat$Z,
    Y = dat$Y,
    N = length(dat$X)
    # Ng = dat$Ng
)
dat2 <- list(
    g = dat$g,
    X = dat$X,
    Z = dat$Z,
    Y = dat$Y,
    Ng = dat$Ng
)

mN <- ulam(
    alist(
        Y ~ bernoulli(p),
        logit(p) <- a + bXY * X + bZY * Z[g], # +Ug[g] +,
        # Ug ~ normal(0, 1),
        a ~ normal(0, 10),
        bXY ~ normal(0, 1),
        bZY ~ normal(0, 1)
    ),
    data = dat2,
    cores = 4,
    chains = 4,
    log_lik = TRUE,
    cmdstan = TRUE
)

mFE <- ulam(
    alist(
        Y ~ bernoulli(p),
        logit(p) <- a[g] + bXY * X + bZY * Z[g], # +Ug[g] +,
        # Ug ~ normal(0, 1),
        a[g] ~ normal(0, 10),
        bXY ~ normal(0, 1),
        bZY ~ normal(0, 1)
    ),
    data = dat2,
    cores = 4,
    chains = 4,
    log_lik = TRUE,
    cmdstan = TRUE
)
precis(mFE, depth = 2)
# par(mfrow = c(1, 2))
postN <- extract.samples(mN)
post <- extract.samples(mFE)
# dens(post$bXY, main = "bXY", ylim = c(0, 12), xlim = c(-1, 1))
# dens(postN$bXY, add = TRUE, col = 2, lt = 2)
# dens(post$bZY, main = "bZY", ylim = c(0, 7), xlim = c(-1, 3))
# dens(postN$bZY, add = TRUE, col = 2, lt = 2)

# Multi-level partial pooling model

mPP <- ulam(
    alist(
        Y ~ bernoulli(p),
        logit(p) <- a[g] + bXY * X + bZY * Z[g], # +Ug[g] +,
        # Ug ~ normal(0, 1),
        transpars > vector[Ng]:a <<- abar + z * tau,
        abar ~ normal(0, 1),
        z[g] ~ normal(0, 1),
        tau ~ exponential(1),
        bXY ~ normal(0, 1),
        bZY ~ normal(0, 1)
    ),
    data = dat2,
    cores = 4,
    chains = 4,
    log_lik = TRUE,
    cmdstan = TRUE,
    sample = TRUE
)
postPP <- extract.samples(mPP)


# Mine new
par(mfrow = c(1, 1))
print(min(dat2$X))
print(max(dat2$X))
plot(NULL, xlab = "range", ylab = "Density", xlim = c(-4, 6), ylim = c(0, 0.7))
N_gr <- 5
for (i in 1:N_gr) {
    print(min(dat2$X[dat2$g == i]))
    print(max(dat2$X[dat2$g == i]))
    dens(dat2$X[dat2$g == i], add = TRUE, lw = i / N_gr * 2)
}

X <- rnorm(1000, 0, 1)
bux <- rexp(1000, 1) #rnorm(1000,5,1)
muu <- rnorm(1000, 0, 1)
musigma <- rexp(1000, 1)
U <- sapply(1:N_groups, function(q) {
    return(rnorm(1000, muu, musigma))
})
X_pred <- sapply(1:N_groups, function(q) return(aX + bux * U[, q]))

plot(NULL, xlab = "range", ylab = "Density", xlim = c(-6, 6), ylim = c(-6, 6))
# for (i in 1:100){
plot(X_pred[, 1] ~ U[, 1])
# }
plot(NULL, xlab = "range", ylab = "Density", xlim = c(-6, 6), ylim = c(0, 0.7))
for (i in 1:N_gr) {
    dens(dat2$X[dat2$g == i], add = TRUE, col = "black", lw = i / N_gr * 2)
    dens(X_pred[, i], add = TRUE, lw = i / N_gr * 2, col = "blue")
}

# Latent Mundlak machine (partial pooling with modeled confound)

mMM <- ulam(
    alist(
        Y ~ bernoulli(p),
        logit(p) <- a[g] + bXY * X + bZY * Z[g] + buy * Ug[g],
        # Ug ~ normal(0, 1),
        transpars > vector[Ng]:a <<- abar + z * tau,
        X ~ normal(mu, sigma),
        mu <- aX + bux * Ug[g],
        vector[Ng]:Ug ~ normal(0, 1),
        sigma ~ exponential(1),
        aX ~ normal(0, 1),
        z[g] ~ normal(0, 1),
        tau ~ exponential(1),
        bXY ~ normal(0, 1),
        bZY ~ normal(0, 1),
        bux ~ normal(0, 1),
        buy ~ normal(0, 1),
        abar ~ normal(0, 1)
    ),
    data = dat2,
    cores = 4,
    chains = 4,
    log_lik = TRUE,
    cmdstan = TRUE,
    sample = TRUE,
    refresh = 0
)
# Zero divergences, but some Rhats large and some variables below 200 ess.
dashboard(mMM)
# results from diagnose
# The following parameters had rank-normalized split R-hat greater than 1.01:
#  Ug[1], Ug[2], Ug[3], Ug[4], Ug[5], Ug[6], Ug[7], Ug[8], Ug[9], Ug[10], Ug[11], Ug[12], Ug[13], Ug[14], Ug[15], Ug[16], Ug[17], Ug[18], Ug[19], Ug[20], Ug[21], Ug[22], Ug[23], Ug[24], Ug[25], Ug[26], Ug[27], Ug[28], Ug[29], Ug[30], aX, bux, buy, abar, a[1], a[4], a[5], a[6], a[8], a[11], a[13], a[15], a[16], a[18], a[21], a[22], a[23], a[29]
#Such high values indicate incomplete mixing and biased estimation.
#You should consider regularizing your model with additional prior information or a more effective parameterization.
mMM@cstanfit$cmdstan_diagnose()

# Richard Mc Elreath
mMMRMCE <- ulam(
    alist(
        Y ~ bernoulli(p),
        logit(p) <- a[g] + bXY * X + bZY * Z[g] + buy * Ug[g],
        # Ug ~ normal(0, 1),
        transpars > vector[Ng]:a <<- abar + z * tau,
        X ~ normal(mu, sigma),
        mu <- aX + bux * Ug[g],
        vector[Ng]:Ug ~ normal(0, 1),
        sigma ~ exponential(1),
        aX ~ normal(0, 1),
        z[g] ~ normal(0, 1),
        tau ~ exponential(1),
        bXY ~ normal(0, 1),
        bZY ~ normal(0, 1),
        bux ~ exponential(1), #normal(0, 1),
        buy ~ normal(0, 1),
        abar ~ normal(0, 1)
    ),
    data = dat2,
    cores = 4,
    chains = 4,
    log_lik = TRUE,
    cmdstan = TRUE,
    sample = TRUE,
    refresh = 0
)

# Better Rhat (below 1.01), not divergences, no below 200
dashboard(mMMRMCE)
# results from diagnose
# No regularization required, all good!
mMMRMCE@cstanfit$cmdstan_diagnose()

N_groups <- max(unique(dat$g))
dat3 <- dat2
dat3$Xbar <- sapply(1:N_groups, function(j) mean(dat3$X[dat3$g == j]))
# Richard Mc Elreath - bars
mMMRMCEbar <- ußlam(
    alist(
        Y ~ bernoulli(p),
        logit(p) <- a[g] + bXY * X + bZY * Z[g] + buy * Xbar[g],
        # Ug ~ normal(0, 1),
        transpars > vector[Ng]:a <<- abar + z * tau,
        aX ~ normal(0, 1),
        z[g] ~ normal(0, 1),
        tau ~ exponential(1),
        bXY ~ normal(0, 1),
        bZY ~ normal(0, 1),
        # bux ~ exponential(1), #normal(0, 1),
        buy ~ normal(0, 1),
        abar ~ normal(0, 1)
    ),
    data = dat3,
    cores = 4,
    chains = 4,
    log_lik = TRUE,
    cmdstan = TRUE,
    sample = TRUE,
    refresh = 0
)

# Better Rhat (below 1.01), not divergences, no below 200
dashboard(mMMRMCEbar)
# results from diagnose
# No regularization required, all good!
mMMRMCEbar@cstanfit$cmdstan_diagnose()

# Half-normals for exponentials
mMMHN <- ulam(
    alist(
        Y ~ bernoulli(p),
        logit(p) <- a[g] + bXY * X + bZY * Z[g] + buy * Ug[g],
        # Ug ~ normal(0, 1),
        transpars > vector[Ng]:a <<- abar + z * tau,
        X ~ normal(mu, sigma),
        mu <- aX + bux * Ug[g],
        vector[Ng]:Ug ~ normal(0, 1),
        sigma ~ half_normal(0, 1),
        aX ~ normal(0, 1),
        z[g] ~ normal(0, 1),
        tau ~ half_normal(0, 1),
        bXY ~ normal(0, 1),
        bZY ~ normal(0, 1),
        bux ~ normal(0, 1), #normal(0, 1),
        buy ~ normal(0, 1),
        abar ~ normal(0, 1)
    ),
    data = dat2,
    cores = 4,
    chains = 4,
    log_lik = TRUE,
    cmdstan = TRUE,
    sample = TRUE,
    refresh = 0
)

# Rhat (still above 1.01), no divergences, some variables below 200 ess_bulk
dashboard(mMMHN)
# results from diagnose
# The following parameters had rank-normalized split R-hat greater than 1.01:
#  Ug[1], Ug[2], Ug[3], Ug[4], Ug[5], Ug[6], Ug[7], Ug[8], Ug[9], Ug[10], Ug[11], Ug[12], Ug[13], Ug[14], Ug[15], Ug[16], Ug[17], Ug[18], Ug[19], Ug[20], Ug[21], Ug[22], Ug[23], Ug[24], Ug[25], Ug[26], Ug[27], Ug[28], Ug[29], Ug[30], aX, bux, buy, abar, a[1], a[3], a[4], a[5], a[7], a[8], a[9], a[11], a[12], a[13], a[14], a[15], a[16], a[18], a[20], a[21], a[22], a[23], a[24], a[25], a[26], a[27], a[28], a[29], a[30]
# Such high values indicate incomplete mixing and biased estimation.
# You should consider regularizing your model with additional prior information or a more effective parameterization.
mMMHN@cstanfit$cmdstan_diagnose()


# Half-normals for exponentials
mMMHNRMCE <- ulam(
    alist(
        Y ~ bernoulli(p),
        logit(p) <- a[g] + bXY * X + bZY * Z[g] + buy * Ug[g],
        # Ug ~ normal(0, 1),
        transpars > vector[Ng]:a <<- abar + z * tau,
        X ~ normal(mu, sigma),
        mu <- aX + bux * Ug[g],
        vector[Ng]:Ug ~ normal(0, 1),
        sigma ~ half_normal(0, 1),
        aX ~ normal(0, 1),
        z[g] ~ normal(0, 1),
        tau ~ half_normal(0, 1),
        bXY ~ normal(0, 1),
        bZY ~ normal(0, 1),
        bux ~ half_normal(0, 1), #normal(0, 1),
        buy ~ normal(0, 1),
        abar ~ normal(0, 1)
    ),
    data = dat2,
    cores = 4,
    chains = 4,
    log_lik = TRUE,
    cmdstan = TRUE,
    sample = TRUE,
    refresh = 0
)

# Rhat (still above 1.01), butt none above 1.04, no divergences, some variables below 200 ess_bulk
dashboard(mMMHNRMCE)
# results from diagnose
# Ug[1], Ug[2], Ug[3], Ug[4], Ug[5], Ug[6], Ug[7], Ug[8], Ug[9], Ug[10], Ug[11], Ug[12], Ug[13], Ug[14], Ug[15], Ug[16], Ug[17], Ug[18], Ug[19], Ug[20], Ug[21], Ug[22], Ug[23], Ug[24], Ug[25], Ug[26], Ug[27], Ug[28], Ug[29], Ug[30], aX, bux, buy, abar, a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9], a[10], a[11], a[12], a[13], a[14], a[15], a[16], a[17], a[18], a[19], a[20], a[21], a[22], a[23], a[24], a[25], a[26], a[27], a[28], a[29], a[30]
# Such high values indicate incomplete mixing and biased estimation.
#You should consider regularizing your model with additional prior information or a more effective parameterization.
mMMHNRMCE@cstanfit$cmdstan_diagnose()

# Half-normals for exponentials
mMMHNRMCEexp <- ulam(
    alist(
        Y ~ bernoulli(p),
        logit(p) <- a[g] + bXY * X + bZY * Z[g] + buy * Ug[g],
        # Ug ~ normal(0, 1),
        transpars > vector[Ng]:a <<- abar + z * tau,
        X ~ normal(mu, sigma),
        mu <- aX + bux * Ug[g],
        vector[Ng]:Ug ~ normal(0, 1),
        sigma ~ half_normal(0, 1),
        aX ~ normal(0, 1),
        z[g] ~ normal(0, 1),
        tau ~ half_normal(0, 1),
        bXY ~ normal(0, 1),
        bZY ~ normal(0, 1),
        bux ~ exponential(1), #normal(0, 1),
        buy ~ normal(0, 1),
        abar ~ normal(0, 1)
    ),
    data = dat2,
    cores = 4,
    chains = 4,
    log_lik = TRUE,
    cmdstan = TRUE,
    sample = TRUE,
    refresh = 0
)

# Rhat (below or equal to 1.01), no divergences, no variables below 200 ess_bulk
dashboard(mMMHNRMCEexp)
# results from diagnose
# No problems detected, all good!
mMMHNRMCEexp@cstanfit$cmdstan_diagnose()

# Half-normals for exponentials
mMMHNML <- ulam(
    alist(
        Y ~ bernoulli(p),
        logit(p) <- a[g] + bXY * X + bZY * Z[g] + buy * Ug[g],
        # Ug ~ normal(0, 1),
        transpars > vector[Ng]:a <<- abar + z * tau,
        X ~ normal(mu, sigma),
        mu <- aX + bux * Ug[g],
        vector[Ng]:Ug ~ normal(0, 1),
        sigma ~ half_normal(0, 1),
        aX ~ normal(0, 1),
        z[g] ~ normal(0, 1),
        tau ~ half_normal(0, 1),
        bXY ~ normal(0, 1),
        bZY ~ normal(0, 1),
        bux ~ normal(buxbar, sigmabuxbar), #normal(0, 1),
        buy ~ normal(0, 1),
        abar ~ normal(0, 1),
        buxbar ~ normal(0, 1),
        sigmabuxbar ~ half_normal(0, 1)
    ),
    data = dat2,
    cores = 4,
    chains = 4,
    log_lik = TRUE,
    cmdstan = TRUE,
    sample = TRUE,
    refresh = 0
)

# Rhat (above 1.01), 20 divergences, some variables below 200 ess_bulk
dashboard(mMMHNML)
# results from diagnose
# The following parameters had rank-normalized split R-hat greater than 1.01:
#  Ug[1], Ug[2], Ug[3], Ug[4], Ug[5], Ug[6], Ug[7], Ug[8], Ug[9], Ug[10], Ug[11], Ug[12], Ug[13], Ug[14], Ug[15], Ug[16], Ug[17], Ug[18], Ug[19], Ug[20], Ug[21], Ug[22], Ug[23], Ug[24], Ug[25], Ug[26], Ug[27], Ug[28], Ug[29], Ug[30], aX, bux, buy, abar, buxbar, a[1], a[2], a[3], a[5], a[6], a[7], a[8], a[9], a[10], a[11], a[12], a[13], a[14], a[15], a[16], a[17], a[18], a[19], a[20], a[21], a[22], a[23], a[24], a[25], a[26], a[27], a[28], a[29], a[30]
#Such high values indicate incomplete mixing and biased estimation.
# You should consider regularizing your model with additional prior information or a more effective parameterization.
mMMHNML@cstanfit$cmdstan_diagnose()


# Half-normals for exponentials
mMMHNMLNC <- ulam(
    alist(
        Y ~ bernoulli(p),
        logit(p) <- a[g] + bXY * X + bZY * Z[g] + buy * Ug[g],
        # Ug ~ normal(0, 1),
        transpars > vector[Ng]:a <<- abar + z * tau,
        vector[N]:X ~ normal(mu, sigma),
        mu <- aX[1] + bux[1] * Ug[g],
        transpars > vector[1]:aX <<- abarux[1] + v[, 1],
        transpars > vector[1]:bux <<- abarux[2] + v[, 2],
        transpars > matrix[1, 2]:v <- compose_noncentered(
            sigma_bux,
            L_rho_bux,
            Z_bux
        ),
        matrix[2, 1]:Z_bux ~ normal(0, 1),
        vector[2]:abarux ~ normal(0, 1),
        cholesky_factor_corr[2]:L_rho_bux ~ lkj_corr_cholesky(4),
        vector[2]:sigma_bux ~ exponential(1),

        vector[Ng]:Ug ~ normal(0, 1),
        sigma ~ exponential(1),
        # aX ~ normal(0, 1),
        z[g] ~ normal(0, 1),
        tau ~ exponential(1),
        bXY ~ normal(0, 1),
        bZY ~ normal(0, 1),
        buy ~ normal(0, 1),
        abar ~ normal(0, 1),
        gq > matrix[2, 2]:Rho_bux <<- Chol_to_Corr(L_rho_bux)
    ),
    data = dat3,
    cores = 4,
    chains = 4,
    log_lik = TRUE,
    cmdstan = TRUE,
    sample = TRUE,
    refresh = 0,
    control = list(adapt_delta = 0.99, max_treedepth = 13)
)

# Rhat (above 1.01), 20 divergences, some variables below 200 ess_bulk
dashboard(mMMHNMLNC)
# results from diagnose
# The following parameters had rank-normalized split R-hat greater than 1.01:
#  Ug[1], Ug[2], Ug[3], Ug[4], Ug[5], Ug[6], Ug[7], Ug[8], Ug[9], Ug[10], Ug[11], Ug[12], Ug[13], Ug[14], Ug[15], Ug[16], Ug[17], Ug[18], Ug[19], Ug[20], Ug[21], Ug[22], Ug[23], Ug[24], Ug[25], Ug[26], Ug[27], Ug[28], Ug[29], Ug[30], aX, bux, buy, abar, buxbar, a[1], a[2], a[3], a[5], a[6], a[7], a[8], a[9], a[10], a[11], a[12], a[13], a[14], a[15], a[16], a[17], a[18], a[19], a[20], a[21], a[22], a[23], a[24], a[25], a[26], a[27], a[28], a[29], a[30]
#Such high values indicate incomplete mixing and biased estimation.
# You should consider regularizing your model with additional prior information or a more effective parameterization.
mMMHNMLNC@cstanfit$cmdstan_diagnose()


mMMpriors <- ulam(
    alist(
        Y ~ bernoulli(p),
        logit(p) <- a[g] + bXY * X + bZY * Z[g] + buy * Ug[g],
        # Ug ~ normal(0, 1),
        transpars > vector[Ng]:a <<- abar + z * tau,
        X ~ normal(mu, sigma),
        mu <- aX + bux * Ug[g],
        vector[Ng]:Ug ~ normal(0, 1),
        sigma ~ exponential(1),
        aX ~ normal(1, 1),
        z[g] ~ normal(0, 1),
        tau ~ exponential(1),
        bXY ~ normal(0, 1),
        bZY ~ normal(0, 1),
        bux ~ normal(buxbar, buxsigma),
        buxbar ~ normal(1, 1),
        buxsigma ~ exponential(1),
        buy ~ normal(0, 1),
        abar ~ normal(0, 1)
    ),
    data = dat2,
    cores = 4,
    chains = 4,
    log_lik = TRUE,
    cmdstan = TRUE,
    sample = TRUE,
    refresh = 0,
    # warmup = 1,
    # iter = 1,
    # sample_prior = TRUE
)

# Prior Predictive
sapply()
# Zero divergences, but some Rhats large and some variables below 200 ess.
dashboard(mMMpriors)
mMMpriors@cstanfit$cmdstan_diagnose()

sim_priors <- function(
    inp_g,
    inp_Ng
) {
    n_samples <- 1000
    tau <- rexp(n_samples, 1)
    bXY <- rnorm(n_samples, 0, 1)
    # bZY <- rnorm(n_samples, 0, 1)

    buxbar <- rnorm(n_samples, 1, 1)
    buxsigma <- rexp(n_samples, 1)

    bux <- rexp(n_samples, 1) #rnorm(n_samples, 0,1)#buxbar, buxsigma)

    buy <- rnorm(n_samples, 0, 1)
    abar <- rnorm(n_samples, 0, 1)

    sigma <- rexp(n_samples, 1)
    aX <- rep(0, n_samples) #rnorm(n_samples, 0, 1)
    z <- sapply(1:inp_Ng, function(q) {
        return(rnorm(n_samples, 0, 1))
    }) #rnorm(n_samples, 0, 1)

    Ug <- sapply(1:inp_Ng, function(q) {
        return(rnorm(n_samples, 0, 1))
    })

    mu <- aX +
        bux *
            sapply(1:inp_Ng, function(q) {
                return(Ug[, q])
            })[, inp_g]

    pred_X <- rnorm(n_samples, mu, sigma)

    a <- abar + z * tau
    p <- inv_logit(
        a[, inp_g] + bXY * pred_X + buy * Ug[, inp_g] #+ bZY * inp_Z[inp_g]
    )
    Y <- rbern(p)

    # Ug ~ normal(0, 1),

    return(list(
        Y = Y,
        p = p,
        a = a,
        mu = mu,
        pred_X = pred_X,
        Ug = Ug,
        abar = abar,
        z = z,
        tau = tau,
        mu = mu,
        sigma = sigma,
        aX = aX,
        bux = bux,
        buy = buy,
        buxsigma = buxsigma,
        buxbar = buxbar,
        bXY = bXY
    ))
}
par(mfrow = c(1, 1))
plot(NULL, xlab = "X_orig", ylab = "Density", xlim = c(-4, 6), ylim = c(0, 0.3))
means <- sapply(dat$g, function(q) {
    return(sim_priors(inp_g = dat$g[i], inp_Ng = dat$Ng)$mu)
})

for (i in 1:dat$Ng) {
    dens(dat$X[dat$g == i], col = "black", lw = i / dat$Ng * 2, add = TRUE)
    means <- sim_priors(inp_g = dat$g[dat$g == i], inp_Ng = dat$Ng)$mu
    dens(means, add = TRUE, col = "blue", lw = i / dat$Ng * 2)
}
# for (i in 1:length(dat$X)) {
#     sim_d_priors <- sim_priors(inp_g=dat$g[i],inp_Ng=dat$Ng)
#     points(dat$X[i],mean(sim_d_priors$mu))
#     pic <- PI(sim_d_priors$mu)
#     print(pic)
#     # lines(c(i,), c(i,))
#     # lines(c(i,), c(i,))
# }

# for (i in 1:length(dat$X)) {
#     sim_d_priors <- sim_priors(inp_g=dat$g[i],inp_Ng=dat$Ng)
#     points(dat$X[i],mean(sim_d_priors$mu))
#     pic <- PI(sim_d_priors$mu)
#     print(pic)
#     # lines(c(i,), c(i,))
#     # lines(c(i,), c(i,))
# }
# mMM2b <- ulam(
#     alist(
#         Y ~ bernoulli(p),
#         logit(p) <- a[g] + bXY * X + bZY * Z[g] + buy * Ug[g],
#         # Ug ~ normal(0, 1),
#         transpars > vector[Ng]:a <<- abar + z * tau,
#         X ~ normal(mu, sigma),
#         mu <- aX + bux * Ug[g],
#         vector[Ng]:Ug ~ normal(0, 1),
#         bux ~ normal(0, 1),
#         aX ~ normal(0, 1),
#         z[g] ~ normal(0, 1),
#         sigma ~ half_normal(0, 1),
#         # bubar ~ normal(0, 1),
#         tau ~ half_normal(0, 1),
#         bXY ~ normal(0, 1),
#         bZY ~ normal(0, 1),
#         buy ~ normal(0, 1),
#         abar ~ normal(0, 1)
#     ),
#     data = dat2,
#     cores = 4,
#     chains = 4,
#     log_lik = TRUE,
#     cmdstan = TRUE,
#     sample = TRUE,
#     refresh = 0 #,
#     # iter = 3000
# )
# dashboard(mMM2b)

# dat3 <- list(
#     g = dat$g,
#     X = dat$X,
#     Z = dat$Z,
#     Y = dat$Y,
#     Ng = dat$Ng,
#     N = length(dat$X)
# )

# mMM3 <- ulam(
#     alist(
#         Y ~ bernoulli(p),
#         logit(p) <- a[g] + bXY * X + bZY * Z[g] + buy * Ug[g],
#         # Ug ~ normal(0, 1),
#         transpars > vector[Ng]:a <<- abar + z * tau,

#         X ~ normal(mu, sigma),
#         mu <- au[1] + bu[1] * Ug[g],
#         vector[Ng]:Ug ~ normal(0, 1),
#         transpars > vector[1]:au <<- abaru[1] + vu[, 1],
#         transpars > vector[1]:bu <<- abaru[2] + vu[, 2],
#         transpars > matrix[1, 2]:vu <- compose_noncentered(
#             sigma_u,
#             L_Rho_u,
#             Z_u
#         ),
#         matrix[2, 1]:Z_u ~ normal(0, 1),
#         vector[2]:abaru ~ normal(0, 1),
#         vector[2]:sigma_u ~ exponential(1),
#         cholesky_factor_corr[2]:L_Rho_u ~ lkj_corr_cholesky(4),
#         sigma ~ exponential(1),
#         # aX ~ normal(0, 1),
#         z[g] ~ normal(0, 1),
#         tau ~ exponential(1),
#         bXY ~ normal(0, 1),
#         bZY ~ normal(0, 1),
#         # bux ~ normal(mubux, sigmabux), #normal(0, 1),
#         # transpars > vector[1]:bux <<- amubar + zmubux * taubux,
#         # bux ~ normal(mubux, taubux),
#         buy ~ normal(0, 1),
#         abar ~ normal(0, 1),
#         # taubux ~ exponential(1),
#         # vector[1]:zmubux ~ normal(0, 1),
#         # amubar ~ normal(0, 1),
#         # Correlation
#         gq > matrix[2, 2]:Rho_u <<- Chol_to_Corr(L_Rho_u)
#     ),
#     data = dat2,
#     cores = 4,
#     chains = 4,
#     log_lik = TRUE,
#     cmdstan = TRUE,
#     sample = TRUE,
#     control = list(max_treedepth = 15, adapt_delta = 0.99)
#     # iter = 2000
#     # refresh = 0
# )
#         dashboard(mMM3)
# }

mMM4 <- ulam(
    alist(
        Y ~ bernoulli(p),
        logit(p) <- a[g] + bXY * X + bZY * Z[g] + buy * Ug[g],
        # Ug ~ normal(0, 1),
        transpars > vector[Ng]:a <<- abar + z * tau,

        X ~ normal(mu, sigma),
        mu <- au[1] + bu[1] * Ug[g],
        vector[Ng]:Ug ~ normal(0, 1),
        transpars > vector[1]:au <<- abaru[1] + vu[, 1],
        transpars > vector[1]:bu <<- abaru[2] + vu[, 2],
        transpars > matrix[1, 2]:vu <- compose_noncentered(
            sigma_u,
            L_Rho_u,
            Z_u
        ),
        matrix[2, 1]:Z_u ~ normal(0, 1),
        vector[2]:abaru ~ normal(0, 1),
        vector[2]:sigma_u ~ half_normal(1),
        cholesky_factor_corr[2]:L_Rho_u ~ lkj_corr_cholesky(2),
        sigma ~ exponential(1),
        # aX ~ normal(0, 1),
        z[g] ~ normal(0, 1),
        tau ~ exponential(1),
        bXY ~ normal(0, 1),
        bZY ~ normal(0, 1),
        # bux ~ normal(mubux, sigmabux), #normal(0, 1),
        # transpars > vector[1]:bux <<- amubar + zmubux * taubux,
        # bux ~ normal(mubux, taubux),
        buy ~ normal(0, 1),
        abar ~ normal(0, 1),
        # taubux ~ exponential(1),
        # vector[1]:zmubux ~ normal(0, 1),
        # amubar ~ normal(0, 1),
        # Correlation
        gq > matrix[2, 2]:Rho_u <<- Chol_to_Corr(L_Rho_u)
    ),
    data = dat2,
    cores = 4,
    chains = 4,
    log_lik = TRUE,
    cmdstan = TRUE,
    sample = TRUE #,
    # iter = 5000 #,
    # control = list(max_treedepth = 12, adapt_delta = 0.99)
    # iter = 2000
    # refresh = 0
)
dashboard(mMM4)

postMM <- extract.samples(mMM)
postMM2 <- extract.samples(mMM2)
postMM3 <- extract.samples(mMM3)
postMM2b <- extract.samples(mMM2b)

par(mfrow = c(2, 2))
dens(post$bXY, main = "bXY", ylim = c(0, 12), xlim = c(-0.3, 1))
dens(postN$bXY, add = TRUE, col = 2, lt = 2)
dens(postPP$bXY, add = TRUE, col = 4, lt = 2)
dens(postMM$bXY, add = TRUE, col = 6, lt = 2, lw = 3)
dens(postMM2$bXY, add = TRUE, col = 8, lt = 2, lw = 3)
dens(postMM3$bXY, add = TRUE, col = 25, lt = 2, lw = 3)
dens(post$bZY, main = "bZY", ylim = c(0, 7), xlim = c(-1, 3))
dens(postN$bZY, add = TRUE, col = 2, lt = 2)
dens(postPP$bZY, add = TRUE, col = 4, lt = 2)
dens(postMM$bZY, add = TRUE, col = 6, lt = 2, lw = 3)
dens(postMM2$bZY, add = TRUE, col = 8, lt = 2, lw = 3)
dens(postMM3$bZY, add = TRUE, col = 25, lt = 2, lw = 3)
dens(postMM2b$bZY, add = TRUE, col = 35, lt = 2, lw = 3)
dens(postMM$bux, col = 6)
dens(postMM2$bux, add = TRUE, col = 8, lw = 3)
dens(postMM3$bux, add = TRUE, col = 25, lw = 3)
dens(postMM2b$bux, add = TRUE, col = 35, lt = 2, lw = 3)

apply(postMM3$bux, 2, mean)
apply(postMM2b$bux, 2, mean)
compare(mMM, mMM2, mMM3, func = "PSIS")
