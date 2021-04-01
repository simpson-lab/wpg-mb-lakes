EXAMPLES <- FALSE # should examples be run?

##' Simulate from the posterior distribution of the mean of a Gamma LS GAM
##'  using Gaussian approximation to the posterior
##'
##' @param model         the fitted GAM
##' @param data          the new data locations you want to get variance for
##' @param nsims         the number of posterior draws wanted - low by default
##'                      to avoid excessive computation, but needs to be
##'                      10,000+ for quantile-based intervals
##' @param unconditional logical; use the smoothness selection corrected version
##'                      of the Bayesian covariance matrix of the model?
`gammals_mean` <- function(model, data, nsims = 100,
                         unconditional  = FALSE, ...) {
    ## Simulate variance from posterior
    sim <- sim_gammals_mean(model = model, data = data,
                            nsims = nsims, unconditional = unconditional)
    ## process results into a tibble
    colnames(sim) <- paste0("sim", seq_len(nsims))
    tbl <- as_tibble(sim) %>%
        bind_cols(data) %>%
        tibble::rowid_to_column(var = "row")
    tbl <- pivot_longer(tbl,
                        cols = matches("^sim"),
                        names_to = "simulation",
                        values_to = "mean")
    tbl
}

##' Simulate from the posterior distribution of the variance of a Gamma LS GAM
##'  using Gaussian approximation to the posterior
##'
##' @param model         the fitted GAM
##' @param data          the new data locations you want to get variance for
##' @param nsims         the number of posterior draws wanted - low by default
##'                      to avoid excessive computation, but needs to be
##'                      10,000+ for quantile-based intervals
##' @param unconditional logical; use the smoothness selection corrected version
##'                      of the Bayesian covariance matrix of the model?
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
    tbl <- pivot_longer(tbl,
                        cols = matches("^sim"),
                        names_to = "simulation",
                        values_to = "variance")
    tbl
}

##' Simulate from the posterior distribution of the derivative of the mean
##' of a location -scale Gamma GAM using Gaussian approximation to the
##' posterior. Derivatives are computed using left finite differences only.
##' Other arguments as per above.
##'
##' @param var character; the variable to shift along for derivatives
##' @param eps numeric; the value to shift data by for finite differences
`gammals_mean_deriv` <- function(model, data, var, nsims = 100,
                                unconditional  = FALSE, eps = 1e-07, ...) {
    
    ## f'(x) = (f(x + eps) - f(x))/eps as eps --> 0
    
    ## prediction matrix
    Xp1 <- predict(model, newdata = data, type = 'lpmatrix')
    ## model parameters
    coefs <- coef(model)
    ## Bayesian covariance matrix
    Vb <- vcov(model, unconditional = unconditional)
    ## which coefs go with the theta linear predictor
    theta_take <- grepl('^s\\.1', colnames(Xp1)) |
        colnames(Xp1) %in% c('(Intercept).1')
    mu_take <- !theta_take
    ## inverse link functions
    ilink_mu <- inv_link(model, parameter = "location") # mu inv link function

    ## Simulate from posterior using Gaussian approximation
    betas <- mvnfast::rmvn(n = nsims,
                           mu = coefs,
                           sigma = Vb)
    
    ## predict variance for data
    mu1 <- est_gammals_mean(betas, Xp1, mu_take, ilink_mu)
    data2 <- data # copy
    ## shift the variable of interest by eps
    data2[[var]] <- data2[[var]] + eps
    ## predict for shifted data
    ## prediction matrix
    Xp2 <- predict(model, newdata = data2, type = 'lpmatrix')
    mu2 <- est_gammals_mean(betas, Xp2, mu_take, ilink_mu)

    ## compute finite differences
    sim_d <- (mu2 - mu1) / eps
    
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

##' Simulate from the posterior distribution of the derivative of the variance
##' of a location -scale Gamma GAM using Gaussian approximation to the
##' posterior. Derivatives are computed using left finite differences only.
##' Other arguments as per above.
##'
##' @param var character; the variable to shift along for derivatives
##' @param eps numeric; the value to shift data by for finite differences
`gammals_var_deriv` <- function(model, data, var, nsims = 100,
                                unconditional  = FALSE, eps = 1e-07, ...) {
    
    ## f'(x) = (f(x + eps) - f(x))/eps as eps --> 0
    
    ## prediction matrix
    Xp1 <- predict(model, newdata = data, type = 'lpmatrix')
    ## model parameters
    coefs <- coef(model)
    ## Bayesian covariance matrix
    Vb <- vcov(model, unconditional = unconditional)
    ## which coefs go with the theta linear predictor
    theta_take <- grepl('^s\\.1', colnames(Xp1)) |
        colnames(Xp1) %in% c('(Intercept).1')
    mu_take <- !theta_take
    ## inverse link functions
    ilink_mu <- inv_link(model, parameter = "location") # mu inv link function
    ilink_theta <- inv_link(model, parameter = "scale") # theta inverse link

    ## Simulate from posterior using Gaussian approximation
    betas <- mvnfast::rmvn(n = nsims,
                           mu = coefs,
                           sigma = Vb)
    
    ## predict variance for data
    var1 <- est_gammals_var(betas, Xp1, mu_take, theta_take,
                            ilink_mu, ilink_theta)
    data2 <- data # copy
    ## shift the variable of interest by eps
    data2[[var]] <- data2[[var]] + eps
    ## predict for shifted data
    ## prediction matrix
    Xp2 <- predict(model, newdata = data2, type = 'lpmatrix')
    var2 <- est_gammals_var(betas, Xp2, mu_take, theta_take,
                            ilink_mu, ilink_theta)

    ## compute finite differences
    sim_d <- (var2 - var1) / eps
    
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

`est_gammals_mean` <- function(betas, Xp, mu_take, ilink_mu) {
    ## subset Xp matrix into mean part
    Xp_mu <- Xp[, mu_take, drop = FALSE]
    
    ## Predict for mean
    fit_mu <- Xp_mu %*% t(betas[, mu_take, drop = FALSE]) # predict on internal scale
    fit_mu <- ilink_mu(fit_mu) # apply g-1() this is just identity so redundant
    ## This model is parameterised in terms of the log-mean, so we still need to
    ## transform to the actual data scale using exp()
    fit_mu <- exp(fit_mu)
    fit_mu
}


`est_gammals_var` <- function(betas, Xp, mu_take, theta_take,
                              ilink_mu, ilink_theta) {
    ## subset Xp matrix into mean and scale parts
    Xp_mu <- Xp[, mu_take, drop = FALSE]
    Xp_theta <- Xp[, theta_take, drop = FALSE]

    ## Predict for mean
    fit_mu <- Xp_mu %*% t(betas[, mu_take, drop = FALSE]) # predict on internal scale
    fit_mu <- ilink_mu(fit_mu) # apply g-1() this is just identity so redundant
    ## This model is parameterised in terms of the log-mean, so we still need to
    ## transform to the actual data scale using exp()
    fit_mu <- exp(fit_mu)
    
    ## Predict for theta
    fit_theta <- Xp_theta %*%
        t(betas[, theta_take, drop = FALSE]) # predict on internal scale
    fit_theta <- ilink_theta(fit_theta) # apply g-1()
    ## fit theta even after using inverse link is log for theta parameter, so
    ## more back transforming
    fit_theta <- exp(fit_theta)

    ## variance is mu * s where s is the scale in sense of rgamma, not theta
    ## From ?rgamma Var(y) = shape * scale^2 = (1 / theta) * (mu * theta)^2
    ## From ?gammals Var(y) = mu * scale = mu * s
    ##   where scale = s = mu * theta. Hence from ?gammals we arrive finally
    ##   at: Var(y) = mu * s = mu * (mu * theta)
    fit_var_draws <- fit_mu * (fit_mu * fit_theta)
    ## return
    fit_var_draws
}

##' The internal workhorse does all the cool stuff 
`sim_gammals_mean` <- function(model, data, nsims = 100, unconditional = FALSE,
                               ...) {
    ## prediction matrix
    Xp <- predict(model, newdata = data, type = 'lpmatrix')
    ## model parameters
    coefs <- coef(model)
    ## Bayesian covariance matrix
    Vb <- vcov(model, unconditional = unconditional)
    ## which coefs go with the theta linear predictor
    theta_take <- grepl('^s\\.1', colnames(Xp)) |
        colnames(Xp) %in% c('(Intercept).1')

    ## Simulate from posterior using Gaussian approximation
    betas <- mvnfast::rmvn(n = nsims,
                           mu = coefs,
                           sigma = Vb)
    
    ## simplify later code so form the compliment to select mean
    ## linear predictor
    mu_take <- !theta_take
    ## subset Xp matrix into mean and theta parts
    Xp_mu <- Xp[, mu_take, drop = FALSE]

    ## Predict for mean
    fit_mu <- Xp_mu %*% t(betas[, mu_take, drop = FALSE]) # predict on internal scale
    ilink_mu <- inv_link(model, parameter = "location") # link function
    fit_mu <- ilink_mu(fit_mu) # apply g-1() this is just identity so redundant
    ## This model is parameterised in terms of the log-mean, so we still need to
    ## transform to the actual data scale using exp()
    fit_mu <- exp(fit_mu)
    fit_mu
}

##' The internal workhorse does all the cool stuff 
`sim_gammals_var` <- function(model, data, nsims = 100, unconditional = FALSE,
                              ...) {
    ## prediction matrix
    Xp <- predict(model, newdata = data, type = 'lpmatrix')
    ## model parameters
    coefs <- coef(model)
    ## Bayesian covariance matrix
    Vb <- vcov(model, unconditional = unconditional)
    ## which coefs go with the theta linear predictor
    theta_take <- grepl('^s\\.1', colnames(Xp)) |
        colnames(Xp) %in% c('(Intercept).1')

    ## Simulate from posterior using Gaussian approximation
    betas <- mvnfast::rmvn(n = nsims,
                           mu = coefs,
                           sigma = Vb)
    
    ## simplify later code so form the compliment to select mean
    ## linear predictor
    mu_take <- !theta_take
    ## subset Xp matrix into mean and theta parts
    Xp_mu <- Xp[, mu_take, drop = FALSE]
    Xp_theta <- Xp[, theta_take, drop = FALSE]

    ## Predict for mean
    fit_mu <- Xp_mu %*% t(betas[, mu_take, drop = FALSE]) # predict on internal scale
    ilink_mu <- inv_link(model, parameter = "location") # link function
    fit_mu <- ilink_mu(fit_mu) # apply g-1() this is just identity so redundant
    ## This model is parameterised in terms of the log-mean, so we still need to
    ## transform to the actual data scale using exp()
    fit_mu <- exp(fit_mu)
    
    ## Predict for theta
    fit_theta <- Xp_theta %*%
        t(betas[, theta_take, drop = FALSE]) # predict on internal scale
    ilink_theta <- inv_link(model, parameter = "scale") # theta inverse link
    fit_theta <- ilink_theta(fit_theta) # apply g-1()
    ## fit scale even after using inverse link is log for theta parameter, so
    ## more back transforming
    fit_theta <- exp(fit_theta)

    ## variance is mu * s where s is the scale in sense of rgamma, not theta
    ## From ?rgamma Var(y) = shape * scale^2 = (1 / theta) * (mu * theta)^2
    ## From ?gammals Var(y) = mu * scale = mu * s
    ##   where scale = s = mu * theta. Hence from ?gammals we arrive finally
    ##   at: Var(y) = mu * s = mu * (mu * theta)
    fit_var_draws <- fit_mu * (fit_mu * fit_theta)
    ## return
    fit_var_draws
}

## Examples -----------------------------------------------------------

if(EXAMPLES) {
    
    library('ggplot2')
    library('readr')
    library('dplyr')
    library('tidyr')
    library('readxl')
    library('gratia')
    ## number of simulations
    K <- 25
    
    ## Need to create `newd` for the prediction locations
    m.gammals <- read_rds('models/lakes-gammals-fs.rds')
    
    ## import data and model (see models.R for info on data pre-processing) ----
    wpg <- read_xlsx('data/wpg/Final Core 1 Summary data for Report.xlsx') %>%
        select(SECTION_NO, YEAR, DIATOX, ALLO, PHEO_B, PHEO_A, CHL_A, CANTH, B_CAR,
               DEPTH_CM) %>%
        mutate(thickness = DEPTH_CM - lag(DEPTH_CM),
               interval = lag(YEAR) - YEAR) %>%
        rename(sample = SECTION_NO)
    colnames(wpg) <- tolower(colnames(wpg))
    wpg <- wpg[-(1:2), ]
    wpg <- select(wpg, -chl_a, -pheo_a) %>%
        pivot_longer(cols = -c('sample', 'depth_cm', 'thickness', 'interval', 'year'),
                     names_to = 'pigment',
                     values_to = 'conc') %>%
        mutate(pigment = factor(pigment),
               lake = 'Lake Winnipeg') %>%
        filter(!is.na(conc))
    
    mb <- read_xlsx('data/mb/Manitoba pigs isotope Core 1 April 2014.xlsx') %>%
        select(SAMPLE, MID_DEPTH_CM, YEAR, ALLOX, DIATOX, CANTH, PHEO_B, BCAROT, PHEO_A,
               CHLA) %>%
        rename(depth_cm = MID_DEPTH_CM,
               allo = ALLOX,
               b_car = BCAROT) %>%
        mutate(thickness = depth_cm - lag(depth_cm),
               interval = lag(YEAR) - YEAR)
    colnames(mb) <- tolower(colnames(mb))
    mb <- mb[-(1:2), ] %>%
        select(-pheo_a, -chla) %>%
        pivot_longer(cols = -c('sample', 'depth_cm', 'thickness', 'interval', 'year'),
                     names_to = 'pigment',
                     values_to = 'conc') %>%
        mutate(pigment = factor(pigment),
               lake = 'Lake Manitoba') %>%
        filter(!is.na(conc))
    
    ## combine lake datasets
    lakes <- rbind(wpg, mb) %>%
        mutate(lake_pigment = interaction(lake, pigment, drop = TRUE),
               lake = factor(lake),
               pigment.expr = case_when(pigment == 'allo' ~ 'Alloxanthin',
                                        pigment == 'b_car' ~ 'beta-carotene',
                                        pigment == 'canth' ~ 'Canthaxanthin',
                                        pigment == 'diatox' ~ 'Diatoxanthin',
                                        pigment == 'pheo_b' ~ 'Pheophytin~b'),
               lake.expr = case_when(lake == 'Lake Manitoba' ~ 'Lake~Manitoba',
                                     lake == 'Lake Winnipeg' ~ 'Lake~Winnipeg'),
               pigment.expr = factor(pigment.expr, levels = c('Diatoxanthin',
                                                              'Alloxanthin',
                                                              'Pheophytin~b',
                                                              'Canthaxanthin',
                                                              'beta-carotene'))) %>%
        filter(year >= 1800)
    
    ## create new data for regularly-spaced predictions
    group_by(lakes, lake) %>% summarise(interval = mean(interval))
    newd <- with(lakes,
                 expand_grid(year = seq(min(year), max(year), by = 1),
                             pigment = unique(lakes$pigment),
                             lake = unique(lakes$lake))) %>%
        mutate(interval = case_when(lake == 'Lake Manitoba' ~ 4.16,
                                    lake == 'Lake Winnipeg' ~ 2.70,
                                    TRUE ~ NA_real_),
               lake_pigment = interaction(lake, pigment, drop = TRUE))

    ## simulate for mean
    set.seed(1)
    mu_sim <- gammals_mean(m.gammals, data = newd, nsims = K)
    
    ## plot those simulations
    ggplot(mu_sim, aes(x = year, y = mean, group = simulation)) +
        geom_line(alpha = 0.2) +
        facet_wrap(~ pigment + lake, ncol = 2, scales = "free_y")
    
    ## simulations from posterior of variance of a Gamma LS model
    set.seed(1)
    var_sim <- gammals_var(m.gammals, data = newd, nsims = K)
    
    ## plot those simulations
    ggplot(var_sim, aes(x = year, y = variance, group = simulation)) +
        geom_line(alpha = 0.2) +
        facet_wrap(~ pigment + lake, ncol = 2, scales = "free_y")
    
    ## derivatives of Gamma LS mean using Gaussian approximation
    ## note these are all being done on the mean scale itself and only via
    ## left finite differences
    set.seed(1)
    mu_d <- gammals_mean_deriv(m.gammals, data = newd, nsims = K,
                                var = "year", eps = 0.00001)
    
    ## plot derivatives
    ggplot(mu_d, aes(x = year, y = derivative, group = simulation)) +
        geom_line(alpha = 0.2) +
        facet_wrap(~ pigment + lake, ncol = 2, scales = "free_y")
    
    ## derivatives of Gamma LS variance using Gaussian approximation
    ## note these are all being done on the variance scale itself ad only via
    ## left finite differences
    var_d <- gammals_var_deriv(m.gammals, data = newd, nsims = K,
                               var = "year", eps = 0.00001)
    
    ## plot derivatives
    ggplot(var_d, aes(x = year, y = derivative, group = simulation)) +
        geom_line(alpha = 0.2) +
        facet_wrap(~ pigment + lake, ncol = 2, scales = "free_y")
    
    ## check curves for K == 1
    if(K == 1) {
        cowplot::plot_grid(ggplot(var_sim, aes(year, variance)) +
                               facet_grid(lake ~ pigment, scales = 'free_y') +
                               geom_line(),
                           ggplot(var_d, aes(year, derivative)) +
                               facet_grid(lake ~ pigment, scales = 'free_y') +
                               geom_hline(yintercept = 0, color = 'red') +
                               geom_line(),
                           ncol = 1)
        
    }
    
    ## Sanity check code to confirm this is doing the right thing
    ## Xp <- predict(m.gammals, newdata = newd, type = 'lpmatrix')
    ## scale_take <- grepl('^s\\.1', colnames(Xp)) |
    ##     colnames(Xp) %in% c('(Intercept).1')
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
    
}
