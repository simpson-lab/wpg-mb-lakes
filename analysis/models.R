library('readr')     # for data reading
library('readxl')    # for data reading
library('dplyr')     # for data wrangling
library('tidyr')     # for data wrangling
library('mgcv')      # for modelling
library('ggplot2')   # for plotting
library('cowplot')   # for plotting
library('gratia')    # for plotting
library('extrafont') # for plot fonts
loadfonts(device = 'win', quiet = TRUE)
theme_set(theme_bw(base_family = 'Times New Roman'))

##################################### change to 1e4 ######################################
K <- 100 # number of simulations

# pigments:
# diatox (diatoxanthin  = diatoms)
# allo   (alloxanthin   = cryptophytes)
# pheo_b (pheophytin b  = chlorophytes)
# canth  (canthaxanthin = potentially N2 fixing cyanos)
# b_car  (beta-carotene = all phytoplankton)

# import Lake Winnipeg data ----
# assume interval between mean dating of two slices = interval within 1st slice
wpg <- read_xlsx('data/wpg/Final Core 1 Summary data for Report.xlsx') %>%
  select(SECTION_NO, YEAR, DIATOX, ALLO, PHEO_B, PHEO_A, CHL_A, CANTH, B_CAR,
         DEPTH_CM) %>%
  mutate(thickness = DEPTH_CM - lag(DEPTH_CM),
         interval = lag(YEAR) - YEAR) %>%
  rename(sample = SECTION_NO)

# rename columns to lower case names
colnames(wpg) <- tolower(colnames(wpg))

# low degradation in early years
plot(chl_a / pheo_a ~ year, wpg)

# remove first two rows due to low degradation
wpg <- wpg[-(1:2), ]

wpg <- select(wpg, -chl_a, -pheo_a)

# change table orientation to long, remove NAs
wpg <- pivot_longer(wpg,
                    cols = -c('sample', 'depth_cm', 'thickness', 'interval', 'year'),
                    names_to = 'pigment',
                    values_to = 'conc') %>%
  mutate(pigment = factor(pigment),
         lake = 'Lake Winnipeg') %>%
  filter(!is.na(conc))

# import Lake Manitoba data ----
mb <- read_xlsx('data/mb/Manitoba pigs isotope Core 1 April 2014.xlsx') %>%
  select(SAMPLE, MID_DEPTH_CM, YEAR, ALLOX, DIATOX, CANTH, PHEO_B, BCAROT, PHEO_A,
         CHLA) %>%
  rename(depth_cm = MID_DEPTH_CM,
         allo = ALLOX,
         b_car = BCAROT) %>%
  mutate(thickness = depth_cm - lag(depth_cm),
         interval = lag(YEAR) - YEAR)

# rename columns to lower case names
colnames(mb) <- tolower(colnames(mb))

# no sudden changes in degradation, remove top two for consistency
plot(chla/pheo_a ~ year, mb)

# remove first two rows due to low degradation
mb <- mb[-(1:2), ] %>%
  select(-pheo_a, -chla)

# change table orientation to long, remove NAs
mb <- pivot_longer(mb,
                   cols = -c('sample', 'depth_cm', 'thickness', 'interval', 'year'),
                   names_to = 'pigment',
                   values_to = 'conc') %>%
  mutate(pigment = factor(pigment),
         lake = 'Lake Manitoba') %>%
  filter(!is.na(conc))

# combine lake datasets
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
                                                        'beta-carotene')))


# plot data ----
# mb starts in 1801, wpg starts much earlier
ggplot(lakes, aes(year, conc)) +
  facet_wrap(pigment ~ lake, scales = 'free_y', ncol = 2) +
  geom_vline(xintercept = 1800, color = 'red') +
  geom_point() +
  scale_color_gradient(low = 'darkorange', high = 'blue') +
  labs(x = NULL, y = 'Concentration')

# remove data prior to 1800
lakes <- filter(lakes, year >= 1800)

# modelling ----
# read in model
m.gammals <- read_rds('models/lakes-gammals-fs.rds')

# gamma location-scale model (20 seconds)
tictoc::tic()
m.gammals <- gam(list(conc ~
                      # formula for mean
                        s(year, lake_pigment, k = 15, bs = 'fs', xt = list(bs = 'cr')),
                      # formula for scale
                      ~ s(interval, k = 10) + # number of years /sample
                        s(year, lake_pigment, k = 10, bs = 'fs')),
                 family = gammals(),
                 data = lakes,
                 method = 'REML',
                 control = gam.control(nthreads = 4))
tictoc::toc()
#saveRDS(m.gammals, 'models/lakes-gammals-fs.rds')

# check models
lakes <- mutate(lakes,
                e = resid(m.gammals),
                estar = abs(conc - m.gammals$fitted.values[, 1]))
appraise(m.gammals)
plot(m.gammals, pages = 1, scale = 0)

cowplot::plot_grid(ggplot(lakes, aes(year, e, group = pigment)) +
                     facet_grid(lake ~ .) +
                     geom_point(),
                   ggplot(lakes, aes(interval, e, group = pigment)) +
                     facet_grid(lake ~ .) +
                     geom_point(),
                   ggplot(lakes, aes(log(conc), e, group = pigment)) +
                     facet_grid(pigment ~ lake, scale = 'free_x') +
                     geom_point(),
                   ggplot(lakes, aes(e, group = pigment)) +
                     facet_grid(lake ~ .) +
                     geom_density())

# create new data for regularly-spaced predictions
group_by(lakes, lake) %>% summarise(interval = mean(interval))
newd <- with(lakes,
             expand_grid(year = seq(min(year), max(year), by = 1),
                         pigment = unique(lakes$pigment),
                         lake = unique(lakes$lake))) %>%
  mutate(interval = case_when(lake == 'Lake Manitoba' ~ 4.16,
                              lake == 'Lake Winnipeg' ~ 2.70,
                              TRUE ~ NA_real_),
         lake_pigment = interaction(lake, pigment, drop = TRUE))

set.seed(1)
slopes.mu <- derivatives(m.gammals, newdata = newd, term = 's(year,lake_pigment)',
                         type = 'central', interval = 'simultaneous')
slopes.shape <- derivatives(m.gammals, newdata = newd, term = 's.1(year,lake_pigment)',
                            interval = 'simultaneous')

# predictions
pred <-
  predict(m.gammals, newd, se.fit = TRUE, type = 'link') %>%
  as.data.frame() %>%
  as_tibble() %>%
  rename(log.mean = fit.1,
         shape = fit.2,
         se.log.mean = se.fit.1,
         se.shape = se.fit.2) %>%
  mutate(mean = exp(log.mean),
         mean.lwr = exp(log.mean - 1.96 * se.log.mean),
         mean.upr = exp(log.mean + 1.96 * se.log.mean),
         s2 = mean * shape,
         signif.mu = case_when(slopes.mu$lower > 0 ~ 'red',
                               slopes.mu$upper < 0 ~ 'blue',
                               TRUE ~ 'transparent'),
         signif.shape = case_when(slopes.shape$lower > 0 ~ 'red',
                                  slopes.shape$upper < 0 ~ 'blue',
                                  TRUE ~ 'transparent'),
         signif.s2 = case_when(signif.mu == 'incr' & signif.shape == 'incr' ~ 'red',
                                   signif.mu == 'decr' & signif.shape == 'decr' ~ 'blue',
                                   TRUE ~ 'transparent')) %>%
  bind_cols(newd) %>%
  mutate(pigment.expr = case_when(pigment == 'allo' ~ 'Alloxanthin',
                                  pigment == 'b_car' ~ 'beta-carotene',
                                  pigment == 'canth' ~ 'Canthaxanthin',
                                  pigment == 'diatox' ~ 'Diatoxanthin',
                                  pigment == 'pheo_b' ~ 'Pheophytin~b'),
         pigment.expr = factor(pigment.expr, levels = c('Diatoxanthin',
                                                        'Alloxanthin',
                                                        'Pheophytin~b',
                                                        'Canthaxanthin',
                                                        'beta-carotene')),
         lake.expr = case_when(lake == 'Lake Manitoba' ~ 'Lake~Manitoba',
                               lake == 'Lake Winnipeg' ~ 'Lake~Winnipeg'))

################# change to 1e4 predictions ########################
# sims <- smooth_samples(m.gammals, n = K, newdata = newd, seed = 1, ncores = 4)
# pred.mean <-
#   filter(sims, !grepl('s.1', smooth)) %>%
#   select(-by_variable,-smooth, -term)
# pred.mean <- pivot_wider(pred.mean, names_from = .x2, values_from = value)
# pred.mean <- pivot_longer(pred.mean, cols = unique(lakes$lake_pigment))
# pred.mean <- filter(pred.mean, !is.na(value))
# pred.mean <-  mutate(pred.mean, sim = exp(`NA` + value)) %>%
#   group_by(.x1, name) %>%
#   summarize(mean.median = quantile(sim, 0.5),
#             mean.lwr = quantile(sim, 0.025),
#             mean.upr = quantile(sim, 0.975)) %>%
#   mutate(year = .x1,
#          lake = substr(name, start = 1, regexpr('\\.', name) - 1),
#          pigment = substr(name, start = regexpr('\\.', name) + 1, stop = nchar(name)),
#          pigment.expr = case_when(pigment == 'allo' ~ 'Alloxanthin',
#                                   pigment == 'b_car' ~ 'beta-carotene',
#                                   pigment == 'canth' ~ 'Canthaxanthin',
#                                   pigment == 'diatox' ~ 'Diatoxanthin',
#                                   pigment == 'pheo_b' ~ 'Pheophytin~b'),
#          pigment.expr = factor(pigment.expr, levels = c('Diatoxanthin',
#                                                         'Alloxanthin',
#                                                         'Pheophytin~b',
#                                                         'Canthaxanthin',
#                                                         'beta-carotene')),
#          lake.expr = case_when(lake == 'Lake Manitoba' ~ 'Lake~Manitoba',
#                                lake == 'Lake Winnipeg' ~ 'Lake~Winnipeg'))
# 
# mutate(sigma2 = )
# group_by(smooth, term, .x1, .x2) %>%
#   summarize(median = median(value),
#             lwr = quantile(value, 0.025),
#             upr = quantile(value, 0.975))
# #predict(m.gammals, newd, se.fit = TRUE, type = 'link') %>%
# as.data.frame() %>%
#   as_tibble() %>%
#   rename(log.mean = fit.1,
#          shape = fit.2,
#          se.log.mean = se.fit.1,
#          se.shape = se.fit.2) %>%
#   mutate(mean = exp(log.mean),
#          mean.lwr = exp(log.mean - 1.96 * se.log.mean),
#          mean.upr = exp(log.mean + 1.96 * se.log.mean),
#          s2 = mean * shape) %>%
#   bind_cols(newd)

# plot predictions
p.mean.mb <- ggplot() +
  facet_grid(pigment.expr ~ lake.expr, scales = 'free', labeller = label_parsed) +
  geom_ribbon(aes(year, ymin = mean.lwr, ymax = mean.upr),
              filter(pred, lake == 'Lake Manitoba'), alpha = 0.3) +
  geom_point(aes(year, conc), filter(lakes, lake == 'Lake Manitoba'), alpha = 0.3) +
  geom_line(aes(year, mean, color = signif.mu), filter(pred, lake == 'Lake Manitoba'),
            lwd = 2) +
  geom_line(aes(year, mean), filter(pred, lake == 'Lake Manitoba')) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = NULL, y = 'Mean concentration') +
  scale_color_manual(values = c('blue', 'red', 'transparent')) +
  theme(legend.position = 'none')
p.mean.wpg <- ggplot() +
  facet_grid(pigment.expr ~ lake.expr, scales = 'free', labeller = label_parsed) +
  geom_ribbon(aes(year, ymin = mean.lwr, ymax = mean.upr),
              filter(pred, lake == 'Lake Winnipeg'), alpha = 0.3) +
  geom_point(aes(year, conc), filter(lakes, lake == 'Lake Winnipeg'), alpha = 0.3) +
  geom_line(aes(year, mean, color = signif.mu), filter(pred, lake == 'Lake Winnipeg'),
            lwd = 2) +
  geom_line(aes(year, mean), filter(pred, lake == 'Lake Winnipeg')) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = NULL, y = 'Mean concentration') +
  scale_color_manual(values = c('blue', 'red', 'transparent')) +
  theme(legend.position = 'none')

p.var.mb <- ggplot() +
  facet_grid(pigment.expr ~ lake.expr, scales = 'free', labeller = label_parsed) +
  geom_line(aes(year, s2, color = signif.s2), filter(pred, lake == 'Lake Manitoba'),
            lwd = 2) +
  geom_line(aes(year, s2), filter(pred, lake == 'Lake Manitoba')) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = NULL, y = 'Concentration variance') +
  scale_color_manual(values = c('blue', 'red', 'transparent')) +
  theme(legend.position = 'none')
p.var.wpg <- ggplot() +
  facet_grid(pigment.expr ~ lake.expr, scales = 'free', labeller = label_parsed) +
  geom_line(aes(year, s2, color = signif.s2), filter(pred, lake == 'Lake Winnipeg'),
            lwd = 2) +
  geom_line(aes(year, s2), filter(pred, lake == 'Lake Winnipeg')) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = NULL, y = 'Concentration variance') +
  scale_color_manual(values = c('blue', 'red', 'transparent')) +
  theme(legend.position = 'none')

p <- plot_grid(p.mean.mb, p.mean.wpg, NULL, p.var.mb, p.var.wpg,
               nrow = 1, rel_widths = c(1, .925, .05, 1, .925))
ggsave('figures/mean-variance-predictions.pdf', plot = p, width = 12, height = 5.6,
       dpi = 300, device = cairo_pdf)

head(predict(m.gammals, se.fit = FALSE, type = 'link', newdata = newd) %>%
       as.data.frame() %>%
       mutate())


# posterior samples for the mean, but the variance?
k <- 1e4
sims <- simulate(m.gammals, nsim = k, newdata = newd) %>%
  as.data.frame() %>%
  bind_cols(newd) %>%
  pivot_longer(-c('year', 'pigment', 'lake', 'lake_pigment'),
               names_to = 'sim') %>%
  group_by(year, pigment, lake) %>%
  summarise(lwr = quantile(value, probs = 0.025),
            upr = quantile(value, probs = 0.975),
            median = median(value),
            mu = mean(value),
            s2 = var(value))

# posterior mean
ggplot(sims) +
  facet_wrap(pigment ~ lake, ncol = 2, scales = 'free') +
  geom_ribbon(aes(year, ymin = lwr, ymax = upr), alpha = 0.1) +
  geom_line(aes(year, median))

# variance in the posterior mean
ggplot(sims) +
  facet_wrap(pigment ~ lake, ncol = 2, scales = 'free') +
  geom_line(aes(year, s2))

# derivatives plot ####
slopes.mu <-
  mutate(slopes.mu,
         name = as.character(fs_var),
         lake = substr(name, start = 1, regexpr('\\.', name) - 1),
         pigment = substr(name, regexpr('\\.', name) + 1, stop = nchar(name)),
         pigment.expr = case_when(pigment == 'allo' ~ 'Alloxanthin',
                                  pigment == 'b_car' ~ 'beta-carotene',
                                  pigment == 'canth' ~ 'Canthaxanthin',
                                  pigment == 'diatox' ~ 'Diatoxanthin',
                                  pigment == 'pheo_b' ~ 'Pheophytin~b'),
         pigment.expr = factor(pigment.expr, levels = c('Diatoxanthin',
                                                        'Alloxanthin',
                                                        'Pheophytin~b',
                                                        'Canthaxanthin',
                                                        'beta-carotene')),
         lake.expr = case_when(lake == 'Lake Manitoba' ~ 'Lake~Manitoba',
                               lake == 'Lake Winnipeg' ~ 'Lake~Winnipeg'))
slopes.shape <-
  mutate(slopes.shape,
         name = as.character(fs_var),
         lake = substr(name, start = 1, regexpr('\\.', name) - 1),
         pigment = substr(name, regexpr('\\.', name) + 1, stop = nchar(name)),
         pigment.expr = case_when(pigment == 'allo' ~ 'Alloxanthin',
                                  pigment == 'b_car' ~ 'beta-carotene',
                                  pigment == 'canth' ~ 'Canthaxanthin',
                                  pigment == 'diatox' ~ 'Diatoxanthin',
                                  pigment == 'pheo_b' ~ 'Pheophytin~b'),
         pigment.expr = factor(pigment.expr, levels = c('Diatoxanthin',
                                                        'Alloxanthin',
                                                        'Pheophytin~b',
                                                        'Canthaxanthin',
                                                        'beta-carotene')),
         lake.expr = case_when(lake == 'Lake Manitoba' ~ 'Lake~Manitoba',
                               lake == 'Lake Winnipeg' ~ 'Lake~Winnipeg'))
ggplot() +
  facet_grid(pigment.expr ~ lake.expr, labeller = label_parsed) +
  geom_hline(yintercept = 0) +
  geom_ribbon(aes(data, ymin = lower, ymax = upper), slopes.mu, alpha = 0.3,
              fill = 'forestgreen') +
  geom_line(aes(data, derivative), slopes.mu, color = 'forestgreen') +
  geom_ribbon(aes(data, ymin = lower, ymax = upper), slopes.shape, alpha = 0.3,
              fill = 'goldenrod') +
  geom_line(aes(data, derivative), slopes.shape, color = 'goldenrod') +
  labs(x = NULL, y = 'Slope')
ggsave('figures/derivatives.png', height = 8, width = 8)
