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
datFE <- list(
    C = d$use.contraception,
    D = as.integer(d$district),
    U = as.integer(d$urban),
    A = standardize(d$age.centered),
    K = d$living.children, # use ordered data here
    alpha = rep(2, max(d$living.children) - 1),
    MAXK = max(d$living.children),
    MAXD = max(d$district),
    Kbar = sapply(1:max(d$district), function(q) {
        return(
            if (length(d$living.children[d$district == q]) > 0) {
                return(
                    sum(d$living.children[d$district == q]) /
                        (sum(as.integer(d$district == 1)) *
                            max(d$living.children))
                )
            } else {
                return(0)
            }
        )
    }),
    Ubar = sapply(1:max(d$district), function(q) {
        return(
            if (length(d$urban[d$district == q]) > 0) {
                return(mean(d$urban[d$district == q]))
            } else {
                return(0)
            }
        )
    })
)
mFE <- ulam(
    alist(
        C ~ bernoulli(p),
        logit(p) <- (abar + z_a * sigma_a) +
            (bkbar + z_bk * sigma_bk) *
                sum(delta_j[1:K]) -
            (bkbar + z_bk * sigma_bk) * Kbar[D] +
            bA * A +
            (bubar + z_bu * sigma_bu) * U -
            (bubar + z_bu * sigma_bu) * Ubar[D],
        bA ~ normal(0, 0.5),
        # a <- ,
        # bK <- ,
        # bU <- ,

        abar ~ normal(0, 1),
        bkbar ~ normal(0, 1),
        bubar ~ normal(0, 1),
        z_a ~ normal(0, 1),
        z_bk ~ normal(0, 1),
        z_bu ~ normal(0, 1),
        c(sigma_a, sigma_bk, sigma_bu) ~ exponential(1),

        vector[MAXK]:delta_j <<- append_row(0, delta),
        simplex[MAXK - 1]:delta ~ dirichlet(alpha)
    ),
    data = datFE,
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
# Seems like they are correlated anyway
plot(
    with(postMFE, (bkbar + z_bk * sigma_bk)),
    with(postMFE, (abar + z_a * sigma_a)),
    col = 2,
    xlab = "bK",
    ylab = "a",
    cex = 2,
    lwd = 4
)
# Multi-level model with confounding variable
datML <- list(
    C = d$use.contraception,
    D = as.integer(d$district),
    U = as.integer(d$urban),
    A = standardize(d$age.centered),
    K = d$living.children, # use ordered data here
    alpha = rep(2, max(d$living.children) - 1),
    MAXK = max(d$living.children),
    MAXD = max(d$district)
)
mML <- ulam(
    alist(
        C ~ bernoulli(p),
        logit(p) <- a[D] +
            bK[D] *
                sum(delta_j[1:K]) -
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

dashboard(mML)
mML@cstanfit$cmdstan_diagnose()
precis(mML, depth = 2)
precis(mML, depth = 2, pars = "delta_j")

par(mfrow = c(1, 1))
postmML <- extract.samples(mML)
dens(postmML$sigma_a, xlab = "std", cex = 4)
dens(postmML$sigma_bu, add = TRUE, col = 2, lw = 2)
dens(postmML$sigma_bk, add = TRUE, col = 4, lw = 2)
curve(dexp(x, 1), add = TRUE, lty = 2)

par(mfrow = c(1, 1))
# Seems like they are correlated anyway
plot(
    apply(postmML$bK, 2, mean),
    apply(postmML$a, 2, mean),
    col = 2,
    xlab = "bK[D]",
    ylab = "a[D]",
    cex = 2,
    lwd = 4
)

# Correlated, try non-centered multilevel model

mML_nc_pp_corr <- ulam(
    alist(
        C ~ bernoulli(p),
        logit(p) <- a[D] +
            bK[D] *
                sum(delta_j[1:K]) -
            bA * A +
            bU[D] * U,
        bA ~ normal(0, 0.5),

        transpars > vector[MAXD]:a <<- abar[1] + v[, 1],
        transpars > vector[MAXD]:bK <<- abar[2] + v[, 2],
        transpars > vector[MAXD]:bU <<- abar[3] + v[, 3],

        vector[3]:abar ~ normal(0, 1),

        transpars > matrix[MAXD, 3]:v <- compose_noncentered(sigma, L_Rho, Z),
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

dashboard(mML_nc_pp_corr)
mML_nc_pp_corr@cstanfit$cmdstan_diagnose()
precis(mML_nc_pp_corr, depth = 2)
precis(mML_nc_pp_corr, depth = 2, pars = "delta_j")

par(mfrow = c(1, 1))
postmML_nc_pp_corr <- extract.samples(mML_nc_pp_corr)
dens(postmML_nc_pp_corr$sigma[, 1], xlab = "std", cex = 4)
dens(postmML_nc_pp_corr$sigma[, 2], add = TRUE, col = 2, lw = 2)
dens(postmML_nc_pp_corr$sigma[, 3], add = TRUE, col = 4, lw = 2)
curve(dexp(x, 1), add = TRUE, lty = 2)

par(mfrow = c(3, 1))
# Seems like they are correlated anyway
plot(
    apply(postmML_nc_pp_corr$bK, 2, mean),
    apply(postmML_nc_pp_corr$a, 2, mean),
    col = 2,
    xlab = "bK[D]",
    ylab = "a[D]",
    cex = 2,
    lwd = 4
)
plot(
    apply(postmML_nc_pp_corr$bU, 2, mean),
    apply(postmML_nc_pp_corr$a, 2, mean),
    col = 2,
    xlab = "bU[D]",
    ylab = "a[D]",
    cex = 2,
    lwd = 4
)

plot(
    apply(postmML_nc_pp_corr$bU, 2, mean),
    apply(postmML_nc_pp_corr$bK, 2, mean),
    col = 2,
    xlab = "bU[D]",
    ylab = "bK[D]",
    cex = 2,
    lwd = 4
)

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
postLMMCF <- extract.samples(mLMMCF)

compare(mPPO, mLMMCF)

### Plotting

par(mfrow = c(5, 1))
plot(
    NULL,
    xlim = c(0, max(dat$D) * 0.1),
    ylim = c(0, 1),
    xlab = "district",
    ylab = "proportion[C]",
    main = "mPPO"
)
res <- rbern(
    1:length(dat$D),
    prob = apply(inv_logit(mPPO$p), 2, mean)
)
prop_pred <- sapply(1:max(dat$D), function(q) {
    return(sum(res[dat$D == q]) / sum(as.integer(dat$D == q)))
})
prop <- sapply(1:max(dat$D), function(q) {
    return(sum(dat$C[dat$D == q]) / sum(as.integer(dat$D == q)))
})
for (i in 1:length(dat$D)) {
    points(
        i * 0.1,
        prop[i],
        col = "black",
        cex = 2,
        lwd = 4
    )
    points(
        i * 0.1,
        prop_pred[i],

        col = 4,
        cex = 2,
        lwd = 4
    )
    lines(
        c(i * 0.1, i * 0.1),
        c(prop[i], prop_pred[i]),
        col = 4,
        lwd = 4,
        lt = 2
    )
}

plot(
    NULL,
    xlim = c(0, max(dat$D) * 0.1),
    ylim = c(0, 1),
    xlab = "district",
    ylab = "proportion[C]",
    main = "mFE"
)
res <- rbern(
    1:length(dat$D),
    prob = apply(inv_logit(postMFE$p), 2, mean)
)
prop_pred <- sapply(1:max(dat$D), function(q) {
    return(sum(res[dat$D == q]) / sum(as.integer(dat$D == q)))
})
prop <- sapply(1:max(dat$D), function(q) {
    return(sum(dat$C[dat$D == q]) / sum(as.integer(dat$D == q)))
})
for (i in 1:length(dat$D)) {
    points(
        i * 0.1,
        prop[i],
        col = "black",
        cex = 2,
        lwd = 4
    )
    points(
        i * 0.1,
        prop_pred[i],

        col = 4,
        cex = 2,
        lwd = 4
    )
    lines(
        c(i * 0.1, i * 0.1),
        c(prop[i], prop_pred[i]),
        col = 4,
        lwd = 4,
        lt = 2
    )
}

plot(
    NULL,
    xlim = c(0, max(dat$D) * 0.1),
    ylim = c(0, 1),
    xlab = "district",
    ylab = "proportion[C]",
    main = "mML"
)
res <- rbern(
    1:length(dat$D),
    prob = apply(inv_logit(postmML$p), 2, mean)
)
prop_pred <- sapply(1:max(dat$D), function(q) {
    return(sum(res[dat$D == q]) / sum(as.integer(dat$D == q)))
})
prop <- sapply(1:max(dat$D), function(q) {
    return(sum(dat$C[dat$D == q]) / sum(as.integer(dat$D == q)))
})
for (i in 1:length(dat$D)) {
    points(
        i * 0.1,
        prop[i],
        col = "black",
        cex = 2,
        lwd = 4
    )
    points(
        i * 0.1,
        prop_pred[i],

        col = 4,
        cex = 2,
        lwd = 4
    )
    lines(
        c(i * 0.1, i * 0.1),
        c(prop[i], prop_pred[i]),
        col = 4,
        lwd = 4,
        lt = 2
    )
}


plot(
    NULL,
    xlim = c(0, max(dat$D) * 0.1),
    ylim = c(0, 1),
    xlab = "district",
    ylab = "proportion[C]",
    main = "ML_nc_pp_corr"
)
res <- rbern(
    1:length(dat$D),
    prob = apply(inv_logit(postmML_nc_pp_corr$p), 2, mean)
)
prop_pred <- sapply(1:max(dat$D), function(q) {
    return(sum(res[dat$D == q]) / sum(as.integer(dat$D == q)))
})
prop <- sapply(1:max(dat$D), function(q) {
    return(sum(dat$C[dat$D == q]) / sum(as.integer(dat$D == q)))
})
for (i in 1:length(dat$D)) {
    points(
        i * 0.1,
        prop[i],
        col = "black",
        cex = 2,
        lwd = 4
    )
    points(
        i * 0.1,
        prop_pred[i],

        col = 4,
        cex = 2,
        lwd = 4
    )
    lines(
        c(i * 0.1, i * 0.1),
        c(prop[i], prop_pred[i]),
        col = 4,
        lwd = 4,
        lt = 2
    )
}


plot(
    NULL,
    xlim = c(0, max(dat$D) * 0.1),
    ylim = c(0, 1),
    xlab = "district",
    ylab = "proportion[C]",
    main = "mLMMCF"
)
res <- rbern(
    1:length(dat$D),
    prob = apply(inv_logit(postLMMCF$p), 2, mean)
)
prop_pred <- sapply(1:max(dat$D), function(q) {
    return(sum(res[dat$D == q]) / sum(as.integer(dat$D == q)))
})
prop <- sapply(1:max(dat$D), function(q) {
    return(sum(dat$C[dat$D == q]) / sum(as.integer(dat$D == q)))
})
for (i in 1:length(dat$D)) {
    points(
        i * 0.1,
        prop[i],
        col = "black",
        cex = 2,
        lwd = 4
    )
    points(
        i * 0.1,
        prop_pred[i],

        col = 4,
        cex = 2,
        lwd = 4
    )
    lines(
        c(i * 0.1, i * 0.1),
        c(prop[i], prop_pred[i]),
        col = 4,
        lwd = 4,
        lt = 2
    )
}

compare(mPPO, mFE, mML, mML_nc_pp_corr, mLMMCF)
mPPO@cstanfit$cmdstan_diagnose() # chains ok # TODO: some issues with parameters
kap <- as.data.frame(precis(mPPO, depth = 3))
lkap <- kap[,
    grepl("rhat", colnames(kap)) |
        grepl("", colnames(kap)) |
        grepl("ess_bulk", colnames(kap))
]
lkap[(lkap$rhat > 1.01) | (lkap$ess_bulk < 200), ]
mFE@cstanfit$cmdstan_diagnose() # chains have divergences 0.2%
kap <- as.data.frame(precis(mFE, depth = 3))
lkap <- kap[,
    grepl("rhat", colnames(kap)) |
        grepl("", colnames(kap)) |
        grepl("ess_bulk", colnames(kap))
]
lkap[(lkap$rhat > 1.01) | (lkap$ess_bulk < 200), ]
mML@cstanfit$cmdstan_diagnose() # chains ok, sigma_bu rhat ~=1.014
kap <- as.data.frame(precis(mML, depth = 3))
lkap <- kap[,
    grepl("rhat", colnames(kap)) |
        grepl("", colnames(kap)) |
        grepl("ess_bulk", colnames(kap))
]
lkap[(lkap$rhat > 1.01) | (lkap$ess_bulk < 200), ]
mML_nc_pp_corr@cstanfit$cmdstan_diagnose() # chains with 1.02>rhat>1.01 # sigma[2], bK[7], v[30,2] # TODO: Something weird with parameters
kap <- as.data.frame(precis(mML_nc_pp_corr, depth = 3))
lkap <- kap[,
    grepl("rhat", colnames(kap)) |
        grepl("", colnames(kap)) |
        grepl("ess_bulk", colnames(kap))
]
lkap[(lkap$rhat > 1.01) | (lkap$ess_bulk < 200), ]
mLMMCF@cstanfit$cmdstan_diagnose() # # chains with 1.02>rhat>1.01 aU
kap <- as.data.frame(precis(mLMMCF, depth = 3))
lkap <- kap[,
    grepl("rhat", colnames(kap)) |
        grepl("", colnames(kap)) |
        grepl("ess_bulk", colnames(kap))
]
lkap[(lkap$rhat > 1.01) | (lkap$ess_bulk < 200), ]
