install.packages(c("coda", "mvtnorm", "devtools", "loo", "dagitty", "shape"))
install.packages("devtools")
devtools::install_github("IRkernel/IRkernel")
IRkernel::installspec()
devtools::install_github("rmcelreath/rethinking", force = TRUE)
library(rethinking)
data(bangladesh)
d <- bangladesh
str(d)

dat <- list(
    C = d$use.contraception,
    D = as.integer(d$district),
    U = as.integer(d$urban),
    A = standardize(d$age.centered),
    K = d$living.children, # use ordered data here
    alpha = rep(2, max(d$living.children) - 1),
    MAXK = max(d$living.children),
    MAXD = max(d$district)
)

# Original model
mPPO <- ulam(
    alist(
        C ~ bernoulli(p),
        logit(p) <- a[D] +
            bK[D] *
                sum(delta_j[1:K]) +
            bA * A +
            bU[D] * U,
        bA ~ normal(0, 0.5),
        transpars > vector[MAXD]:a <<- abar[1] + v[, 1],
        transpars > vector[MAXD]:bK <<- abar[2] + v[, 2],
        transpars > vector[MAXD]:bU <<- abar[3] + v[, 3],

        transpars > matrix[MAXD, 3]:v <- compose_noncentered(sigma, L_Rho, Z),
        vector[3]:abar ~ normal(0, 1),
        matrix[3, MAXD]:Z ~ normal(0, 1),
        vector[3]:sigma ~ exponential(1),
        cholesky_factor_corr[3]:L_Rho ~ lkj_corr_cholesky(4),
        vector[MAXK]:delta_j <<- append_row(0, delta),
        simplex[MAXK - 1]:delta ~ dirichlet(alpha),
        gq > matrix[3, 3]:Rho <<- Chol_to_Corr(L_Rho)
    ),
    data = dat,
    chains = 4,
    cores = 4,
    cmdstan = TRUE,
    log_lik = TRUE
)

dashboard(mPPO)
mPPO@cstanfit$cmdstan_diagnose()
precis(mPPO, depth = 2)

# Fixed effects model with confounding variable
mFE <- ulam(
    alist(
        C ~ bernoulli(p),
        logit(p) <- a[D] +
            bK[D] *
                sum(delta_j[1:K]) +
            bA * A +
            bU[D] * U,
        bA ~ normal(0, 0.5),
        transpars > vector[MAXD]:a <<- abar + z_a * sigma_a,
        transpars > vector[MAXD]:bK <<- bkbar + z_bk * sigma_bk,
        transpars > vector[MAXD]:bU <<- bubar + z_bu * sigma_bu,

        abar ~ normal(0, 1),
        bkbar ~ normal(0, 1),
        bubar ~ normal(0, 1),
        vector[MAXD]:z_a ~ normal(0, 1),
        vector[MAXD]:z_bk ~ normal(0, 1),
        vector[MAXD]:z_bu ~ normal(0, 1),
        c(sigma_a, sigma_bk, sigma_bu) ~ exponential(1),

        vector[MAXK]:delta_j <<- append_row(0, delta),
        simplex[MAXK - 1]:delta ~ dirichlet(alpha)
    ),
    data = dat,
    chains = 4,
    cores = 4,
    cmdstan = TRUE,
    log_lik = TRUE
)

dashboard(mFE)
mFE@cstanfit$cmdstan_diagnose()
precis(mFE, depth = 2)
precis(mFE, depth = 2, pars = "delta_j")

par(mfrow = c(1, 1))
postMFE <- extract.samples(mFE)
dens(postMFE$sigma_a, xlab = "std", cex = 4)
dens(postMFE$sigma_bu, add = TRUE, col = 2, lw = 2)
dens(postMFE$sigma_bk, add = TRUE, col = 4, lw = 2)
curve(dexp(x, 1), add = TRUE, lty = 2)

par(mfrow = c(1, 1))
plot(
    apply(postMFE$bK, 2, mean),
    apply(postMFE$a, 2, mean),
    col = 2,
    xlab = "bK[D]",
    ylab = "a[D]",
    cex = 2,
    lwd = 4
)
# Multi-level model with confounding variable

# Mundlak model with confounding variable

# Trial

# Latent Mundlak model with confounding variable
dat <- list(
    C = d$use.contraception,
    D = as.integer(d$district),
    U = as.integer(d$urban),
    A = standardize(d$age.centered),
    K = d$living.children, # use ordered data here
    alpha = rep(2, max(d$living.children) - 1),
    MAXK = max(d$living.children),
    MAXD = max(d$district)
)
mLMMCF <- ulam(
    alist(
        C ~ bernoulli(p),
        logit(p) <- a[D] +
            bA * A +
            bU * U +
            bK * sum(delta_j[1:K]) +
            bUG * UG[D],

        # parameters for top model
        transpars > vector[MAXD]:a <<- abar + z * tau,
        abar ~ normal(0, 1),
        vector[MAXD]:z ~ normal(0, 1),
        tau ~ exponential(1),

        bK ~ normal(0, 1),
        bA ~ normal(0, 0.5),
        bU ~ normal(0, 1),
        bUG ~ normal(0, 1),

        # U model
        U ~ bernoulli(pU),
        logit(pU) <- aU + bUD * UG[D],
        aU ~ normal(0, 1),
        bUD ~ exponential(1),

        # K model
        K ~ ordered_logistic(phi, cutpoints),
        phi <- aKK + bKK * UG[D],
        bKK ~ exponential(1),
        aKK ~ normal(0, 1),
        cutpoints ~ normal(0, 1),

        # UG model
        vector[MAXD]:UG ~ normal(0, 1),

        # the rest
        vector[MAXK]:delta_j <<- append_row(0, delta),
        simplex[MAXK - 1]:delta ~ dirichlet(alpha) #,
        # # gq > matrix[2, 2]:Rho <<- Chol_to_Corr(L_Rho)
        #bK * sum(delta(1:K)) +
    ),
    data = dat,
    chains = 4,
    cores = 4,
    cmdstan = TRUE,
    log_lik = TRUE,
    sample = TRUE
)

dashboard(mLMMCF)
mLMMCF@cstanfit$cmdstan_diagnose()
precis(mLMMCF, depth = 2)

compare(mPPO, mLMMCF)
