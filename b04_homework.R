install.packages(c("coda", "mvtnorm", "devtools", "loo", "dagitty", "shape"))
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
    #    A = standardize(d$age.centered),
    #    K = d$living.children
)

# Original model
mFE <- ulam(
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
        # transpars > vector[MAXD]:bA <<- abar[4] + v[, 4],

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

dashboard(mFE)
mFE@cstanfit$cmdstan_diagnose()
precis(mFE, depth = 2)

# Fixed effects model with confounding variable
mFECF <- ulam(
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
        # transpars > vector[MAXD]:bA <<- abar[4] + v[, 4],

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

dashboard(mFE)
mFE@cstanfit$cmdstan_diagnose()
precis(mFE, depth = 2)


# Multi-level model with confounding variable
mMLCF <- ulam(
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
        # transpars > vector[MAXD]:bA <<- abar[4] + v[, 4],

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

dashboard(mMLCF)
mMLCF@cstanfit$cmdstan_diagnose()
precis(mMLCF, depth = 2)

# Mundlak model with confounding variable
mMMCF <- ulam(
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
        # transpars > vector[MAXD]:bA <<- abar[4] + v[, 4],

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

dashboard(mMMCF)
mMMCF@cstanfit$cmdstan_diagnose()
precis(mMMCF, depth = 2)


# Mundlak model with confounding variable
dat <- list(
    C = d$use.contraception,
    D = as.integer(d$district),
    U = as.integer(d$urban),
    A = standardize(d$age.centered),
    # K = d$living.children, # use ordered data here
    # alpha = rep(2, max(d$living.children) - 1),
    # MAXK = max(d$living.children),
    MAXD = max(d$district)
    #    A = standardize(d$age.centered),
    #    K = d$living.children
)
mLMMCF <- ulam(
    alist(
        C ~ bernoulli(p),
        logit(p) <- a[D] + bA * A + bU * U + bUG * UG[D],

        # parameters for top model
        transpars > vector[MAXD]:a <<- abar + z * tau,
        abar ~ normal(0, 1),
        vector[MAXD]:z ~ normal(0, 1),
        tau ~ exponential(1),

        # bK ~ normal(0, 1),
        bA ~ normal(0, 0.5),
        bU ~ normal(0, 1),
        bUG ~ normal(0, 1),

        # U model
        U ~ bernoulli(pU),
        logit(pU) <- aU + bUD * UG[D],
        aU ~ normal(0, 1),
        bUD ~ exponential(1),

        # # K model
        # K ~ dordlogit(phi, cutpoints),
        # phi <- aKK + bKK * UG[D],
        # bKK ~ exponential(1),
        # aKK ~ normal(0, 1),
        # cutpoints ~ normal(2.5, 3),

        # UG model
        vector[MAXD]:UG ~ normal(0, 1) #,

        # # the rest
        # vector[K]:delta_j <<- append_row(0, delta),
        # simplex[K - 1]:delta ~ dirichlet(alpha) #,
        # # gq > matrix[2, 2]:Rho <<- Chol_to_Corr(L_Rho)
        #bK * sum(delta(1:K)) +
    ),
    data = dat,
    chains = 4,
    cores = 4,
    cmdstan = TRUE,
    log_lik = TRUE
)

dashboard(mLMMCF)
mLMMCF@cstanfit$cmdstan_diagnose()
precis(mLMMCF, depth = 2)
