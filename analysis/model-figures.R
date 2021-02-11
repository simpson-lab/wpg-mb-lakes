library('readr')  # for reading rds files
library('readxl') # for data reading
library('dplyr')  # for data wrangling
library('tidyr')  # for data wrangling
library('gratia') # for plotting
source('analysis/default-figure-styling.R')

# import data and model (see models.R for info on data pre-processing) ----
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
                                                        'beta-carotene'))) %>%
  filter(year >= 1800)

# read in model
m.gammals <- read_rds('models/lakes-gammals-fs.rds')

# predictions ----
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

# find derivatives
set.seed(1)
slopes.mu <- derivatives(m.gammals, newdata = newd, term = 's(year,lake_pigment)',
                         type = 'central', interval = 'simultaneous') %>%
  rename(year = data,
         lake_pigment = fs_var,
         lower.mean = lower,
         upper.mean = upper) %>%
  select(year, lake_pigment, derivative, se, crit, lower.mean, upper.mean)

slopes.shape <- derivatives(m.gammals, newdata = newd, term = 's.1(year,lake_pigment)',
                            interval = 'simultaneous') %>%
  rename(year = data,
         lake_pigment = fs_var,
         lower.shape = lower,
         upper.shape = upper) %>%
  select(year, lake_pigment, derivative, se, crit, lower.shape, upper.shape)

# predictions
pred <-
  predict(m.gammals, newd, se.fit = TRUE, type = 'link') %>%
  as.data.frame() %>%
  as_tibble() %>%
  rename(log.mean = fit.1,
         shape = fit.2,
         se.log.mean = se.fit.1,
         se.shape = se.fit.2) %>%
  bind_cols(newd) %>%
  left_join(slopes.mu, by = c('year', 'lake_pigment')) %>%
  left_join(slopes.shape, by = c('year', 'lake_pigment')) %>%
  mutate(mean = exp(log.mean),
         mean.lwr = exp(log.mean - 1.96 * se.log.mean),
         mean.upr = exp(log.mean + 1.96 * se.log.mean),
         shape.lwr = shape - 1.96 * se.shape,
         shape.upr = shape + 1.96 * se.shape,
         s2 = mean * shape,
         signif.mu = lower.mean > 0 | upper.mean < 0,
         signif.shape = lower.shape > 0 | upper.shape < 0,
         signif.s2 = TRUE) %>%
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
                               lake == 'Lake Winnipeg' ~ 'Lake~Winnipeg')) %>%
  arrange(lake_pigment) %>%
  mutate(segm.mu.bool = signif.mu != lag(signif.mu) | lake_pigment != lag(lake_pigment),
         segm.mu = 1,
         segm.shape.bool = signif.shape != lag(signif.shape) |
           lake_pigment != lag(lake_pigment),
         segm.shape = 1,
         segm.s2.bool = signif.s2 != lag(signif.s2) | lake_pigment != lag(lake_pigment),
         segm.s2 = 1)

## different groups for each significant segment
# mean
for(i in 2:nrow(pred)) {
  if(pred$segm.mu.bool[i]) {
    pred$segm.mu[i] <- pred$segm.mu[i - 1] + 1
  } else {
    pred$segm.mu[i] <- pred$segm.mu[i - 1]
  }
}

# shape
for(i in 2:nrow(pred)) {
  if(pred$segm.shape.bool[i]) {
    pred$segm.shape[i] <- pred$segm.shape[i - 1] + 1
  } else {
    pred$segm.shape[i] <- pred$segm.shape[i - 1]
  }
}

# variance
for(i in 2:nrow(pred)) {
  if(pred$segm.s2.bool[i]) {
    pred$segm.s2[i] <- pred$segm.s2[i - 1] + 1
  } else {
    pred$segm.s2[i] <- pred$segm.s2[i - 1]
  }
}

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

# plot predictions ####
# mean
p.mean.mb <- ggplot() +
  facet_grid(pigment.expr ~ lake.expr, scales = 'free', labeller = label_parsed) +
  geom_ribbon(aes(year, ymin = mean.lwr, ymax = mean.upr),
              filter(pred, lake == 'Lake Manitoba'), alpha = 0.3) +
  
  # lines
  geom_line(aes(year, mean, group = segm.mu), color = 'red', # mean highlight
            filter(pred, lake == 'Lake Manitoba', signif.mu), lwd = 2) +
  geom_line(aes(year, mean), filter(pred, lake == 'Lake Manitoba')) +
  
  # datapoints
  geom_point(aes(year, conc), filter(lakes, lake == 'Lake Manitoba'), alpha = 0.3) +
  
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = 'Year C.E.', y = expression(Mean~concentration~(nmol~g^{-1}~C))) +
  theme(strip.text.y = element_blank())

p.mean.wpg <- ggplot() +
  facet_grid(pigment.expr ~ lake.expr, scales = 'free', labeller = label_parsed) +
  geom_ribbon(aes(year, ymin = mean.lwr, ymax = mean.upr),
              filter(pred, lake == 'Lake Winnipeg'), alpha = 0.3) +
  # lines
  geom_line(aes(year, mean, group = segm.mu), color = 'red', # mean highlight
            filter(pred, lake == 'Lake Winnipeg', signif.mu), lwd = 2) +
  geom_line(aes(year, mean), filter(pred, lake == 'Lake Winnipeg')) +
  
  # datapoints
  geom_point(aes(year, conc), filter(lakes, lake == 'Lake Winnipeg'), alpha = 0.3) +
  
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = NULL, y = NULL) +
  theme(strip.text.y = element_blank())

# mean with variance highlight
p.mean.mb2 <- ggplot() +
  facet_grid(pigment.expr ~ lake.expr, scales = 'free', labeller = label_parsed) +
  geom_ribbon(aes(year, ymin = mean.lwr, ymax = mean.upr),
              filter(pred, lake == 'Lake Manitoba'), alpha = 0.3) +
  
  # lines
  # geom_line(aes(year, mean, group = segm.s2), color = '#5E8BDE', # variance highlight
  #           filter(pred, lake == 'Lake Manitoba', signif.), lwd = 1.5,  alpha = 0.5) +
  geom_line(aes(year, mean, group = segm.mu), color = '#5E8BDE', # TEMP variance highlight
            filter(pred, lake == 'Lake Manitoba', signif.mu), lwd = 5) +
  geom_line(aes(year, mean, group = segm.mu), color = 'red', # mean highlight
            filter(pred, lake == 'Lake Manitoba', signif.mu), lwd = 2) +
  geom_line(aes(year, mean), filter(pred, lake == 'Lake Manitoba')) +
  
  # datapoints
  geom_point(aes(year, conc), filter(lakes, lake == 'Lake Manitoba'), alpha = 0.3) +
  
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = 'Year C.E.', y = expression(Mean~concentration~(nmol~g^{-1}~C)))

p.mean.wpg2 <- ggplot() +
  facet_grid(pigment.expr ~ lake.expr, scales = 'free', labeller = label_parsed) +
  geom_ribbon(aes(year, ymin = mean.lwr, ymax = mean.upr),
              filter(pred, lake == 'Lake Winnipeg'), alpha = 0.3) +
  # lines
  # geom_line(aes(year, mean, group = segm.s2), color = '#5E8BDE', # variance highlight
  #           filter(pred, lake == 'Lake Winnipeg', signif.), lwd = 1.5,  alpha = 0.5) +
  geom_line(aes(year, mean, group = segm.mu), color = '#5E8BDE', # TEMP variance highlight
            filter(pred, lake == 'Lake Winnipeg', signif.mu), lwd = 5) +
  geom_line(aes(year, mean, group = segm.mu), color = 'red', # mean highlight
            filter(pred, lake == 'Lake Winnipeg', signif.mu), lwd = 2) +
  geom_line(aes(year, mean), filter(pred, lake == 'Lake Winnipeg')) +
  
  # datapoints
  geom_point(aes(year, conc), filter(lakes, lake == 'Lake Winnipeg'), alpha = 0.3) +
  
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = NULL, y = NULL)

p.means.nolabs <-
  plot_grid(p.mean.mb2 + theme(axis.title.x = element_blank(),
                               strip.text.y = element_text(color = 'transparent')),
                     p.mean.wpg2,
                     nrow = 1) %>%
  plot_grid(get_plot_component(p.mean.mb, pattern = 'xlab-b'),
            nrow = 2,
            rel_heights = c(0.95, 0.05))
#p2pdf('mean-predictions-nolabs.pdf', p.means.nolabs, scale = 2, y.plots = 1.25)

p.means <- plot_grid(get_plot_component(p.mean.mb, pattern = 'ylab-l'),
                            p.mean.mb2 + theme(axis.title = element_blank(),
                                               strip.text.y = element_text(color = 'transparent')),
                            NULL,
                            p.mean.wpg2,
                            rel_widths = c(0.15, 1, 0, 1),
                            nrow = 1) %>%
  plot_grid(get_plot_component(p.mean.mb, pattern = 'xlab-b'),
            nrow = 2,
            rel_heights = c(0.95, 0.05)) + 
  draw_text(LABELS[1:10],
            x = sort(rep(c(.07, 0.54), 5)),
            y = rep(seq(.95, by = -0.178, length.out = 5), 2),
            family = 'serif')

#p2pdf('mean-predictions.pdf', p.means, scale = 2, y.plots = 1.25)

# shape
p.shape.mb <- ggplot() +
  facet_grid(pigment.expr ~ lake.expr, labeller = label_parsed) +
  geom_line(aes(year, shape, group = segm.shape),
            filter(pred, lake == 'Lake Manitoba', signif.shape), lwd = 2, color = 'grey30') +
  geom_line(aes(year, shape), filter(pred, lake == 'Lake Manitoba')) +
  scale_y_continuous(limits = c(0, 5)) +
  labs(x = NULL, y = expression(Shape~parameter~(nmol~g^{-1}~C))) +
  theme(strip.background.y = element_blank(),
        strip.text.y = element_blank())
p.shape.wpg <- ggplot() +
  facet_grid(pigment.expr ~ lake.expr, labeller = label_parsed) +
  geom_line(aes(year, shape, group = segm.shape),
            filter(pred, lake == 'Lake Winnipeg', signif.shape), lwd = 2, color = 'grey30') +
  geom_line(aes(year, shape), filter(pred, lake == 'Lake Winnipeg')) +
  scale_y_continuous(limits = c(0, 5)) +
  labs(x = NULL, y = NULL) +
  theme(strip.background.y = element_blank(),
        strip.text.y = element_blank())

# variance
p.var.mb <- ggplot() +
  facet_grid(pigment.expr ~ lake.expr, scales = 'free', labeller = label_parsed) +
  geom_line(aes(year, s2), filter(pred, lake == 'Lake Manitoba', signif.s2),
            color = '#5E8BDE', lwd = 2) +
  geom_line(aes(year, s2), filter(pred, lake == 'Lake Manitoba')) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = NULL, y = expression(paste(Concentration~variance~(nmol^2~g^{-2}~C)))) +
  theme(strip.background.y = element_blank(),
        strip.text.y = element_blank())
p.var.wpg <- ggplot() +
  facet_grid(pigment.expr ~ lake.expr, scales = 'free', labeller = label_parsed) +
  geom_line(aes(year, s2), filter(pred, lake == 'Lake Winnipeg', signif.s2),
            color = '#5E8BDE', lwd = 2) +
  geom_line(aes(year, s2), filter(pred, lake == 'Lake Winnipeg')) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = NULL, y = NULL)

# full figure
p.full <- plot_grid(plot_grid(p.mean.mb + xlab(NULL), p.mean.wpg, NULL,
                              p.shape.mb, p.shape.wpg, NULL,
                              p.var.mb, p.var.wpg,
                              nrow = 1, rel_widths = c(1, .95, .1, 1, .95, .1, 1, 1)),
                    get_plot_component(p.mean.mb, pattern = 'xlab-b'),
                    nrow = 2,
                    rel_heights = c(0.95, 0.05))
#p2pdf('mean-shape-variance-predictions.pdf', p.full, scale = 2, x.plots = 2)

head(predict(m.gammals, se.fit = FALSE, type = 'link', newdata = newd) %>%
       as.data.frame() %>%
       mutate())

# posterior samples for the mean, but the variance?
K <- 1e4
sims <- simulate(m.gammals, nsim = K, newdata = newd) %>%
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
         name = as.character(lake_pigment),
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
         name = as.character(lake_pigment),
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
p.deriv <- 
  ggplot() +
  facet_grid(pigment.expr ~ lake.expr, labeller = label_parsed) +
  geom_hline(yintercept = 0) +
  geom_ribbon(aes(year, ymin = lower.mean, ymax = upper.mean), slopes.mu, alpha = 0.2,
              fill = 'forestgreen') +
  geom_line(aes(year, derivative), slopes.mu, color = 'forestgreen') +
  geom_ribbon(aes(year, ymin = lower.shape, ymax = upper.shape), slopes.shape, alpha = 0.2,
              fill = 'goldenrod') +
  geom_line(aes(year, derivative), slopes.shape, color = 'goldenrod') +
  labs(x = 'Year C. E.', y = 'Slope')
# p2pdf('derivatives.pdf', p.deriv, scale = 2.5)
