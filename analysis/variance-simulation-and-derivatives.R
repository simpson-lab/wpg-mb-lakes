## This is needed but not sure you already load it
library('tidyr')
library('ggplot2')
## feel free to separate out code from functions into how you want it for your
## repo

##' Simulate from the posterior distribution of the variance of a Gamma LS GAM
##'  using Gaussian approximation to the posterior
##'
##' @param model the fitted GAM
##' @param data the new data locations you want to get variance for
##' @param nsims the number of posterior draws wanted - low by default to avoid
##'   excessive computation, but needs to be 10,000+ for quantile based
##'   intervals
##' @param unconditional logical; use the smoothness selection corrected version
##'   of the Bayesian covariance matrix of the model?
`gammals_var` <- function(model, data, nsims = 100,
                          unconditional  = FALSE, ...) {
    ## Simulate variance from posterior
    sim <- sim_gammals_var(model = model, data = data,
                           nsims = nsims, unconditional = unconditional)
    ## process results into a tibble
    colnames(sim) <- paste0("sim", seq_len(nsims))
    tbl <- as_tibble(sim) %>%
        bind_cols(data) %>%
        tibble::rowid_to_column(var = "row")
    tbl <- pivot_longer(tbl, cols = matches("^sim"),
                        names_to = "simulation",
                        values_to = "variance")
    tbl
}

##' Simulate from the posterior distribution of the derivative of the variance of
##'  a Gamma LS GAM using Gaussian approximation to the posterior. Derivatives
##'  are computed using left finite differences only. Other arguments as per
##'  above.
##'
##' @param var character; the variable to shift along for derivatives
##' @param eps numeric; the value to shift data by for finite differences
`gammals_var_deriv` <- function(model, data, var, nsims = 100,
                                unconditional  = FALSE, eps = 1e-07, ...) {
    ## simulate predict variance for data
    var1 <- sim_gammals_var(model = model, data = data,
                            nsims = nsims, unconditional = unconditional)
    data2 <- data # copy
    ## shift the variable of interest by eps
    data2[[var]] <- data2[[var]] + eps
    ## predict for shifted data
    var2 <- sim_gammals_var(model = model, data = data2,
                            nsims = nsims, unconditional = unconditional)

    ## compute finite differences - think I have this the right way round
    sim_d <- var2 - var1
    ## process into a tibble
    colnames(sim_d) <- paste0("sim", seq_len(nsims))
    tbl <- as_tibble(sim_d) %>%
        bind_cols(data) %>%
        tibble::rowid_to_column(var = "row")
    tbl <- pivot_longer(tbl, cols = matches("^sim"),
                        names_to = "simulation",
                        values_to = "derivative")
    tbl
}

##' The internal workhorse does all the cool stuff 
`sim_gammals_var` <- function(model, data, nsims = 100, unconditional = FALSE,
                              ...) {
    ## prediction matrix
    Xp <- predict(model, newdata = data, type = 'lpmatrix')
    ## Bayesian covariance matrix
    Vb <- vcov(model, unconditional = unconditional)
    ## which coefs go with the scale linear predictor
    scale_take <- grepl('^s\\.1', colnames(Xp)) |
        colnames(Xp) %in% c('(Intercept).1')
    ## simply later code so form the compliment to select mean linear predictor
    mu_take <- !scale_take
    ## subset Xp matrix into mean and scale parts
    Xp_mu <- Xp[, mu_take, drop = FALSE]
    Xp_scale <- Xp[, scale_take, drop = FALSE]
    ## subset Bayesian VCOV matrix into mean and scale parts
    Vb_mu <- Vb[mu_take, mu_take, drop = FALSE]
    Vb_scale <- Vb[scale_take, scale_take, drop = FALSE]
    ## model parameters
    coefs <- coef(model)
    ## split into mean and scale parts
    coefs_mu <- coefs[mu_take]
    coefs_scale <- coefs[scale_take]

    ## Simulate from posterior using Gaussian approximation
    betas_mu <- mvnfast::rmvn(n = nsims,
                              mu = coefs_mu,
                              sigma = Vb_mu)
    betas_scale <- mvnfast::rmvn(n = nsims,
                                 mu = coefs_scale,
                                 sigma = Vb_scale)

    ## Predict for mean
    fit_mu <- Xp_mu %*% t(betas_mu) ## predict on internal scale
    ilink_mu <- inv_link(model, parameter = "location") ## link function
    fit_mu <- ilink_mu(fit_mu) # apply g-1() this is just identity so redundant
    ## This model is parameterised in terms of the log-mean, so we still need to
    ## transform to the actual data scale using exp()
    fit_mu <- exp(fit_mu)
    
    ## Predict for scale
    fit_scale <- Xp_scale %*% t(betas_scale) # predict on internal scale
    ilink_scale <- inv_link(model, parameter = "scale") # does this work? Yes
    fit_scale <- ilink_scale(fit_scale) # apply g-1()
    ## fit scale even after usig inverse link is log for scale parameter, so
    ## more back transforming
    fit_scale <- exp(fit_scale)

    ## variance is mean * scale
    fit_var_draws <- fit_mu * fit_scale
    ## return
    fit_var_draws
}

## simulations from posterior of variance of a Gamma LS model
var_sim <- gammals_var(m.gammals, data = newd, nsims = 100)

## plot those simulations
ggplot(var_sim, aes(x = year, y = variance, group = simulation)) +
    geom_line(alpha = 0.2) +
    facet_grid(pigment ~ lake, scales = "free_y")

## derivatives of Gamma LS variance using Gaussian approximation
## note these are all being done on the variance scale itself ad only via
## left finite differences
var_d <- gammals_var_deriv(m.gammals, data = newd, nsims = 100, var = "year",
                           eps = 0.0001)

## plot derivatives
ggplot(var_d, aes(x = year, y = derivative, group = simulation)) +
    geom_line(alpha = 0.2) +
    facet_grid(pigment ~ lake, scales = "free_y")

## Sanity check code to confirm this is doing the right thing
## Xp <- predict(m.gammals, newdata = newd, type = 'lpmatrix')
## scale_take <- grepl('^s\\.1', colnames(Xp)) | colnames(Xp) %in% c('(Intercept).1')
## beta_mu <- coef(m.gammals)[!scale_take]
## beta_scale <- coef(m.gammals)[scale_take]
## fit_mu <- Xp[, !scale_take] %*% beta_mu
## fit_scale <- Xp[, scale_take] %*% beta_scale

## ilink_mu <- inv_link(m.gammals, parameter = "location")
## ilink_scale <- inv_link(m.gammals, parameter = "scale")

## fit_mu <- exp(ilink_mu(fit_mu))
## fit_scale <- exp(ilink_scale(fit_scale))

## fit_var <- fit_mu * fit_scale

## head(fit_var)
