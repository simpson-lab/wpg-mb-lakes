library('readr')  # for reading rds files
library('readxl') # for data reading
library('dplyr')  # for data wrangling
library('tidyr')  # for data wrangling
library('gratia') # for plotting
source('analysis/default-figure-styling.R')
source('analysis/variance-simulation-and-derivatives.R')

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
  select(SAMPLE, MID_DEPTH_CM, YEAR, ALLOX, DIATOX, CANTH, PHEO_B, BCAROT,
         PHEO_A, CHLA) %>%
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
m.gammals <- readRDS('models/lakes-gammals-fs.rds')

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
slopes.mu <-
  gammals_mean_deriv(model = m.gammals, data = newd, var = 'year',nsims =1e4)%>%
  group_by(year, lake_pigment) %>%
  summarize(mu.deriv = median(derivative),
            lwr.mu.deriv = quantile(derivative, probs = 0.025),
            upr.mu.deriv = quantile(derivative, probs = 0.975),
            .groups = 'drop')

slopes.s2 <-
  gammals_var_deriv(model = m.gammals, data = newd, var = 'year',nsims = 1e4)%>%
  group_by(year, lake_pigment) %>%
  summarize(s2.deriv = median(derivative),
            lwr.s2.deriv = quantile(derivative, probs = 0.025),
            upr.s2.deriv = quantile(derivative, probs = 0.975),
            .groups = 'drop')

# predictions
mu <- gammals_mean(model = m.gammals, data = newd, nsims = 1e4) %>%
  group_by(year, lake_pigment) %>%
  summarize(mu = median(mean),
            lwr.mu = quantile(mean, probs = 0.025),
            upr.mu = quantile(mean, probs = 0.975),
            .groups = 'drop')
s2 <- gammals_var(model = m.gammals, data = newd, nsims = 1e4) %>%
  group_by(year, lake_pigment) %>%
  summarize(s2 = median(variance),
            lwr.s2 = quantile(variance, probs = 0.025),
            upr.s2 = quantile(variance, probs = 0.975),
            .groups = 'drop') %>%
  mutate(sd = sqrt(s2),
         lwr.sd = sqrt(lwr.s2),
         upr.sd = sqrt(upr.s2))

pred <-
  newd %>%
  left_join(mu, by = c('year', 'lake_pigment')) %>%
  left_join(s2, by = c('year', 'lake_pigment')) %>%
  left_join(slopes.mu, by = c('year', 'lake_pigment')) %>%
  left_join(slopes.s2, by = c('year', 'lake_pigment')) %>%
  mutate(signif.mu = lwr.mu.deriv > 0 | upr.mu.deriv < 0,
         signif.s2 = lwr.s2.deriv > 0 | upr.s2.deriv < 0) %>%
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
  mutate(segm.mu.bool = signif.mu != lag(signif.mu) |
           lake_pigment != lag(lake_pigment),
         segm.mu = 1,
         segm.s2.bool = signif.s2 != lag(signif.s2) |
           lake_pigment != lag(lake_pigment),
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

# plot predictions ####
# mean
p.mu.mb <- ggplot() +
  facet_grid(pigment.expr ~ lake.expr, scales = 'free',
             labeller = label_parsed) +
  geom_ribbon(aes(year, ymin = lwr.mu, ymax = upr.mu),
              filter(pred, lake == 'Lake Manitoba'), alpha = 0.3) +
  
  # lines
  geom_line(aes(year, mu, group = segm.mu), color = 'red', # mean highlight
            filter(pred, lake == 'Lake Manitoba', signif.mu), lwd = 2) +
  geom_line(aes(year, mu), filter(pred, lake == 'Lake Manitoba')) +
  
  # datapoints
  geom_point(aes(year,conc), filter(lakes, lake=='Lake Manitoba'),alpha = 0.3) +
  
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = 'Year C.E.', y = expression(Mean~concentration~(nmol~g^{-1}~C))) +
  theme(strip.text.y = element_blank())

p.mu.wpg <- ggplot() +
  facet_grid(pigment.expr ~ lake.expr, scales = 'free', labeller=label_parsed) +
  geom_ribbon(aes(year, ymin = lwr.mu, ymax = upr.mu),
              filter(pred, lake == 'Lake Winnipeg'), alpha = 0.3) +
  # lines
  geom_line(aes(year, mu, group = segm.mu), color = 'red', # mean highlight
            filter(pred, lake == 'Lake Winnipeg', signif.mu), lwd = 2) +
  geom_line(aes(year, mu), filter(pred, lake == 'Lake Winnipeg')) +
  
  # datapoints
  geom_point(aes(year, conc), filter(lakes, lake=='Lake Winnipeg'), alpha=0.3) +
  
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = NULL, y = NULL) +
  theme(strip.text.y = element_blank())

# mean with variance highlight
p.mu.mb2 <- ggplot() +
  facet_grid(pigment.expr ~ lake.expr, scales = 'free', labeller=label_parsed) +
  geom_ribbon(aes(year, ymin = lwr.mu, ymax = upr.mu),
              filter(pred, lake == 'Lake Manitoba'), alpha = 0.3) +
  
  # lines
  geom_line(aes(year, mu, group = segm.s2), color = '#5E8BDE', # s2 highlight
            filter(pred, lake == 'Lake Manitoba', signif.s2), lwd = 2.5) +
  geom_line(aes(year, mu, group = segm.mu), color = 'red', # mean highlight
            filter(pred, lake == 'Lake Manitoba', signif.mu), lwd = 2) +
  geom_line(aes(year, mu), filter(pred, lake == 'Lake Manitoba')) +
  
  # datapoints
  geom_point(aes(year, conc), filter(lakes, lake == 'Lake Manitoba'),alpha=0.3)+
  
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = 'Year C.E.', y = expression(Mean~concentration~(nmol~g^{-1}~C)))

p.mu.wpg2 <- ggplot() +
  facet_grid(pigment.expr ~ lake.expr, scales = 'free', labeller = label_parsed) +
  geom_ribbon(aes(year, ymin = lwr.mu, ymax = upr.mu),
              filter(pred, lake == 'Lake Winnipeg'), alpha = 0.3) +
  # lines
  geom_line(aes(year, mu, group = segm.s2), color = '#5E8BDE', # s2 highlight
            filter(pred, lake == 'Lake Winnipeg', signif.s2), lwd = 2.5) +
  geom_line(aes(year, mu, group = segm.mu), color = 'red', # mean highlight
            filter(pred, lake == 'Lake Winnipeg', signif.mu), lwd = 2) +
  geom_line(aes(year, mu), filter(pred, lake == 'Lake Winnipeg')) +
  
  # datapoints
  geom_point(aes(year, conc), filter(lakes, lake == 'Lake Winnipeg'), alpha = 0.3) +
  
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = NULL, y = NULL)

p.mus.nolabs <-
  plot_grid(p.mu.mb2 + theme(axis.title.x = element_blank(),
                               strip.text.y = element_text(color = 'transparent')),
            p.mu.wpg2,
            nrow = 1) %>%
  plot_grid(get_plot_component(p.mu.mb, pattern = 'xlab-b'),
            nrow = 2,
            rel_heights = c(0.95, 0.05))
#p2pdf('mean-predictions-nolabs.pdf', p.mus.nolabs, scale = 2, y.plots = 1.25)

p.mus <- plot_grid(get_plot_component(p.mu.mb, pattern = 'ylab-l'),
                     p.mu.mb2 + theme(axis.title = element_blank(),
                                        strip.text.y = element_text(color = 'transparent')),
                     NULL,
                     p.mu.wpg2,
                     rel_widths = c(0.15, 1, 0, 1),
                     nrow = 1) %>%
  plot_grid(get_plot_component(p.mu.mb, pattern = 'xlab-b'),
            nrow = 2,
            rel_heights = c(0.95, 0.05)) + 
  draw_text(LABELS[1:10],
            x = sort(rep(c(.07, 0.54), 5)),
            y = rep(seq(.95, by = -0.178, length.out = 5), 2),
            family = 'serif')

#p2pdf('mean-predictions.pdf', p.mus, scale = 2, y.plots = 1.25)

# variance
p.var.mb <- ggplot() +
  facet_grid(pigment.expr ~ lake.expr, scales = 'free', labeller=label_parsed) +
  geom_ribbon(aes(year, ymin = lwr.s2, ymax = upr.s2),
              filter(pred, lake == 'Lake Manitoba'), alpha = 0.3) +
  geom_line(aes(year, s2, group = segm.s2),
            filter(pred, lake == 'Lake Manitoba', signif.s2),
            color = '#5E8BDE', lwd = 2) +
  geom_line(aes(year, s2), filter(pred, lake == 'Lake Manitoba')) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = NULL,
       y = expression(paste(Concentration~variance~(nmol^2~g^{-2}~C)))) +
  theme(strip.background.y = element_blank(),
        strip.text.y = element_blank())

p.var.wpg <- ggplot() +
  facet_grid(pigment.expr ~ lake.expr, scales = 'free', labeller=label_parsed) +
  geom_ribbon(aes(year, ymin = lwr.s2, ymax = upr.s2),
              filter(pred, lake == 'Lake Winnipeg'), alpha = 0.3) +
  geom_line(aes(year, s2, group = segm.s2),
            filter(pred, lake == 'Lake Winnipeg', signif.s2),
            color = '#5E8BDE', lwd = 2) +
  geom_line(aes(year, s2), filter(pred, lake == 'Lake Winnipeg')) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = NULL, y = NULL)

# full figure
p.full <- plot_grid(plot_grid(p.mu.mb + xlab(NULL), p.mu.wpg, NULL,
                              p.var.mb, p.var.wpg,
                              nrow = 1,
                              rel_widths = c(1, .95, .1, 1, .95)),
                    get_plot_component(p.mu.mb, pattern = 'xlab-b'),
                    nrow = 2,
                    rel_heights = c(0.95, 0.05))
#p2pdf('mean-variance-predictions.pdf', p.full, scale = 2, x.plots = 4/3)
