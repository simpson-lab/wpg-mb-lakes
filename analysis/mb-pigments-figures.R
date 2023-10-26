library('readxl') # for data reading
library('dplyr')  # for data wrangling
library('tidyr')  # for data wrangling
library('gratia') # for plotting
source('analysis/default-figure-styling.R')
source('analysis/variance-simulation-and-derivatives.R')

# import data and model ----
mb <- bind_rows(read_xlsx('data/mb/Manitoba pigs isotope Core 1 April 2014.xlsx') %>%
                  mutate(core = 'Core~1'),
                read_xlsx('data/mb/Manitoba pigs isotope Core 2 April 2014.xlsx') %>%
                  mutate(core = 'Core~2',
                         CHLA = CHL_A),
                read_xlsx('data/mb/Manitoba pigs isotope Core 3 April 2014.xlsx') %>%
                  mutate(core = 'Core~3'),
                read_xlsx('data/mb/Manitoba pigs isotope Core 4 April 2014.xlsx') %>%
                  mutate(core = 'Core~4')) %>%
  select(core, YEAR, ALLOX, DIATOX, CANTH, PHEO_B, BCAROT, CHL_PHEO) %>%
  rename(allo = ALLOX,
         b_car = BCAROT) %>%
  mutate(interval = lag(YEAR) - YEAR) %>%
  filter(interval > 0) %>%
  rename_with(tolower)

## plot small facet of chla/pheoa
chlpheo <- 
  ggplot(mb, aes(x = year, y = chl_pheo))+ 
  geom_line(linewidth = 0.75)+ 
  geom_point(size = 1.5)+
  facet_wrap(~ core, scales = 'free_y', labeller = label_parsed)+
  expand_limits(y = 0) +
  labs(y = expression(Chlorophyll~a/Pheophytin~a), x = 'Year C.E.')+ 
  theme_bw()+
  theme(strip.background = element_blank())+
  theme(strip.text = element_text(size = 11))
chlpheo

p2pdf('chlpheo.pdf', chlpheo, scale = 2, y.plots=1.25)

mb <- mb %>%
  select(! chl_pheo) %>%
  pivot_longer(cols = -c('core', 'interval', 'year'),
               names_to = 'pigment',
               values_to = 'conc') %>%
  mutate(pigment = factor(pigment)) %>%
  filter(!is.na(conc)) %>%
  mutate(core_pigment = interaction(core, pigment, sep = '_'),
         pigment.expr = case_when(pigment == 'allo' ~ 'Cryptophytes',
                                  pigment == 'b_car' ~ 'Total~production',
                                  pigment == 'canth' ~ 'Cyanobacteria',
                                  pigment == 'diatox' ~ 'Diatoms',
                                  pigment == 'pheo_b' ~ 'Chlorophytes'),
         pigment.expr = factor(pigment.expr, levels = c('Diatoms',
                                                        'Cryptophytes',
                                                        'Chlorophytes',
                                                        'Cyanobacteria',
                                                        'Total~production')))

# import the pigment model
m <- readRDS('models/mb-pigments-gammals.rds')

# create predictions ----
newd <- expand_grid(year = seq(1801, 2010, length.out = 400),
                    core = unique(mb$core),
                    pigment = unique(mb$pigment)) %>%
  mutate(core_pigment = interaction(core, pigment, sep = '_')) %>%
  left_join(group_by(mb, core) %>% summarize(interval = mean(interval)),
            by = 'core')

# find derivatives
set.seed(1)

# splitting into two because memory demand is too high with a single value
newd.12 <- filter(newd, core %in% c('Core~1', 'Core~2'))
newd.34 <- filter(newd, core %in% c('Core~3', 'Core~4'))

slopes.mu.12 <-
  gammals_mean_deriv(model = m, data = newd.12, var = 'year', nsims = 1e4) %>%
  group_by(year, core_pigment) %>%
  summarize(mu.deriv = median(derivative),
            lwr.mu.deriv = quantile(derivative, probs = 0.025),
            upr.mu.deriv = quantile(derivative, probs = 0.975),
            .groups = 'drop')

slopes.mu.34 <-
  gammals_mean_deriv(model = m, data = newd.34, var = 'year', nsims = 1e4) %>%
  group_by(year, core_pigment) %>%
  summarize(mu.deriv = median(derivative),
            lwr.mu.deriv = quantile(derivative, probs = 0.025),
            upr.mu.deriv = quantile(derivative, probs = 0.975),
            .groups = 'drop')
slopes.mu <- bind_rows(slopes.mu.12, slopes.mu.34)

slopes.s2.12 <-
  gammals_var_deriv(model = m, data = newd.12, var = 'year', nsims = 1e4) %>%
  group_by(year, core_pigment) %>%
  summarize(s2.deriv = median(derivative),
            lwr.s2.deriv = quantile(derivative, probs = 0.025),
            upr.s2.deriv = quantile(derivative, probs = 0.975),
            .groups = 'drop')
slopes.s2.34 <-
  gammals_var_deriv(model = m, data = newd.34, var = 'year', nsims = 1e4) %>%
  group_by(year, core_pigment) %>%
  summarize(s2.deriv = median(derivative),
            lwr.s2.deriv = quantile(derivative, probs = 0.025),
            upr.s2.deriv = quantile(derivative, probs = 0.975),
            .groups = 'drop')
slopes.s2 <- bind_rows(slopes.s2.12, slopes.s2.34)

rm(slopes.mu.12, slopes.mu.34, slopes.s2.12, slopes.s2.34)

# predictions
mu.12 <- gammals_mean(model = m, data = newd.12, nsims = 1e4) %>%
  group_by(year, core_pigment) %>%
  summarize(mu = median(mean),
            lwr.mu = quantile(mean, probs = 0.025),
            upr.mu = quantile(mean, probs = 0.975),
            .groups = 'drop')
mu.34 <- gammals_mean(model = m, data = newd.34, nsims = 1e4) %>%
  group_by(year, core_pigment) %>%
  summarize(mu = median(mean),
            lwr.mu = quantile(mean, probs = 0.025),
            upr.mu = quantile(mean, probs = 0.975),
            .groups = 'drop')
mu <- bind_rows(mu.12, mu.34)
rm(mu.12, mu.34)

s2.12 <- gammals_var(model = m, data = newd.12, nsims = 1e4) %>%
  group_by(year, core_pigment) %>%
  summarize(s2 = median(variance),
            lwr.s2 = quantile(variance, probs = 0.025),
            upr.s2 = quantile(variance, probs = 0.975),
            .groups = 'drop') %>%
  mutate(sd = sqrt(s2),
         lwr.sd = sqrt(lwr.s2),
         upr.sd = sqrt(upr.s2))
s2.34 <- gammals_var(model = m, data = newd.34, nsims = 1e4) %>%
  group_by(year, core_pigment) %>%
  summarize(s2 = median(variance),
            lwr.s2 = quantile(variance, probs = 0.025),
            upr.s2 = quantile(variance, probs = 0.975),
            .groups = 'drop') %>%
  mutate(sd = sqrt(s2),
         lwr.sd = sqrt(lwr.s2),
         upr.sd = sqrt(upr.s2))
s2 <- bind_rows(s2.12, s2.34)
rm(s2.12, s2.34)

pred <-
  newd %>%
  left_join(mu, by = c('year', 'core_pigment')) %>%
  left_join(s2, by = c('year', 'core_pigment')) %>%
  left_join(slopes.mu, by = c('year', 'core_pigment')) %>%
  left_join(slopes.s2, by = c('year', 'core_pigment')) %>%
  mutate(signif.mu = lwr.mu.deriv > 0 | upr.mu.deriv < 0,
         signif.s2 = lwr.s2.deriv > 0 | upr.s2.deriv < 0) %>%
  mutate(pigment.expr = case_when(pigment == 'allo' ~ 'Cryptophytes',
                                  pigment == 'b_car' ~ 'Total~production',
                                  pigment == 'canth' ~ 'Cyanobacteria',
                                  pigment == 'diatox' ~ 'Diatoms',
                                  pigment == 'pheo_b' ~ 'Chlorophytes'),
         pigment.expr = factor(pigment.expr, levels = c('Diatoms',
                                                        'Cryptophytes',
                                                        'Chlorophytes',
                                                        'Cyanobacteria',
                                                        'Total~production'))) %>%
  arrange(core_pigment) %>%
  mutate(segm.mu.bool = signif.mu != lag(signif.mu) |
           core_pigment != lag(core_pigment),
         segm.mu = 1,
         segm.s2.bool = signif.s2 != lag(signif.s2) |
           core_pigment != lag(core_pigment),
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

# variance
for(i in 2:nrow(pred)) {
  if(pred$segm.s2.bool[i]) {
    pred$segm.s2[i] <- pred$segm.s2[i - 1] + 1
  } else {
    pred$segm.s2[i] <- pred$segm.s2[i - 1]
  }
}


# create plots ----
p.mu <-
  mutate(pred, upr.mu = if_else(upr.mu <= 300, upr.mu, 300)) %>%
  ggplot() +
  facet_grid(pigment.expr ~ core, scales = 'free_y', labeller = label_parsed) +
  geom_point(aes(year, conc), mb, alpha = 0.3) +
  geom_ribbon(aes(year, ymin = lwr.mu, ymax = upr.mu), alpha = 0.3) +
  geom_line(aes(year, mu, group = segm.mu), filter(pred, signif.mu),
            color = 'red', lwd = 2) +
  geom_line(aes(year, mu)) +
  scale_x_continuous('Year C.E.', breaks=seq(1800, 2000, by = 10),
                     labels = c(1800, rep('', 4), 1850, rep('', 4),
                                1900, rep('', 4), 1950, rep('', 4), 2000))+
  scale_color_brewer(type = 'qual', palette = 6) +
  scale_fill_brewer(type = 'qual', palette = 6) +
  scale_y_continuous(expression(Concentration~(nmol~g^{-1}~C)),
                     limits = c(0, NA)) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size=11),
        axis.title = element_text(size=11))
p.mu

# truncate CIs to a max of 3000
p.s2 <-
  mutate(pred,
         s2 = if_else(s2 <= 3000, s2, NA_real_),
         upr.s2 = if_else(upr.s2 <= 3000, upr.s2, 3000),
         upr.s2 = if_else(pigment.expr == 'Total~production' &
                            upr.s2 > 500, 500, upr.s2)) %>%
  ggplot() +
  facet_grid(pigment.expr ~ core, scales = 'free_y', labeller = label_parsed) +
  geom_ribbon(aes(year, ymin = lwr.s2, ymax = upr.s2), alpha = 0.3) +
  geom_line(aes(year, s2, group = segm.s2), filter(pred, signif.s2),
            color = '#5E8BDE', lwd = 2) +
  geom_line(aes(year, s2)) +
  scale_color_brewer(type = 'qual', palette = 6) +
  scale_fill_brewer(type = 'qual', palette = 6) +
  ylim(c(0, NA)) +
  labs(x = 'Year C.E.',
       y = expression(paste(Concentration~variance~(nmol^2~g^{-2}~C)))) +
  theme(strip.text = element_text(size=12),
        axis.title = element_text(size=12))

# p2pdf('mb-pigments.pdf', p.mu, width = 5, height = 3.5, scale = 2)
# p2pdf('mb-pigments-var.pdf', p.s2, width = 5, height = 3.5, scale = 2)
p.full <- plot_grid(p.mu, p.s2, ncol = 1, labels = c('a.', 'b.'))
# p2pdf('mb-pigments-mean-variance.pdf', p.full, width = 5, height = 7, scale=2)

# figure of akinetes from Core 1
akin <-
  read_xlsx('data/mb/Manitoba pigs isotope Core 1 April 2014.xlsx') %>%
  mutate(pigment.expr = 'Akinetes') %>%
  ggplot(aes(YEAR, CYANOS)) +
  facet_grid(pigment.expr ~ ., scales = 'free_y', labeller = label_parsed) +
  geom_point(alpha = 0.3) +
  scale_x_continuous('Year C.E.', breaks=seq(1800, 2000, by = 10),
                     labels = c(1800, rep('', 4), 1850, rep('', 4),
                                1900, rep('', 4), 1950, rep('', 4), 2000))+
  scale_color_brewer(type = 'qual', palette = 6) +
  scale_fill_brewer(type = 'qual', palette = 6) +
  scale_y_continuous(expression(atop(Concentration,
                                     x10^{3}~(g^{-1}~wet))),
                     limits = c(0, NA)) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size=11),
        axis.title = element_text(size=11))
akin

# p2pdf('mb-akinetes.pdf', akin, width = 5, height = 4)
