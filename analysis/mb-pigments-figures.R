library('readxl') # for data reading
library('dplyr')  # for data wrangling
library('tidyr')  # for data wrangling
library('gratia') # for plotting
source('analysis/default-figure-styling.R')

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
  select(core, YEAR, ALLOX, DIATOX, CANTH, PHEO_B, BCAROT, PHEO_A, CHLA) %>%
  rename(allo = ALLOX,
         b_car = BCAROT) %>%
  mutate(interval = lag(YEAR) - YEAR) %>%
  filter(interval > 0)
colnames(mb) <- tolower(colnames(mb))
mb <- mb %>%
  select(-pheo_a, -chla) %>%
  pivot_longer(cols = -c('core', 'interval', 'year'),
               names_to = 'pigment',
               values_to = 'conc') %>%
  mutate(pigment = factor(pigment)) %>%
  filter(!is.na(conc)) %>%
  mutate(core_pigment = interaction(core, pigment, sep = '_'),
         pigment.expr = case_when(pigment == 'allo' ~ 'Alloxanthin',
                                  pigment == 'b_car' ~ 'beta-carotene',
                                  pigment == 'canth' ~ 'Canthaxanthin',
                                  pigment == 'diatox' ~ 'Diatoxanthin',
                                  pigment == 'pheo_b' ~ 'Pheophytin~b'),
         pigment.expr = factor(pigment.expr, levels = c('Diatoxanthin',
                                                        'Alloxanthin',
                                                        'Pheophytin~b',
                                                        'Canthaxanthin',
                                                        'beta-carotene')))

m <- readr::read_rds('models/mb-pigments-gammals.rds')

# create predictions ----
pred <- expand_grid(year = seq(1801, 2010, length.out = 400),
                    core = unique(mb$core),
                    pigment = unique(mb$pigment)) %>%
  mutate(core_pigment = interaction(core, pigment, sep = '_')) %>%
  left_join(group_by(mb, core) %>% summarize(interval = mean(interval)),
            by = 'core')
pred <- cbind(pred,
              predict(m, newdata = pred, type = 'link', se.fit = TRUE)) %>%
  as_tibble() %>%
  rename(fit.mu = fit.1,
         fit.shape = fit.2,
         se.mu = se.fit.1,
         se.shape = se.fit.2) %>%
  mutate(mu = exp(fit.mu),
         mu.lwr = exp(fit.mu - 1.96 * se.mu),
         mu.upr = exp(fit.mu + 1.96 * se.mu),
         var = mu * fit.shape,
         pigment.expr = case_when(pigment == 'allo' ~ 'Alloxanthin',
                                  pigment == 'b_car' ~ 'beta-carotene',
                                  pigment == 'canth' ~ 'Canthaxanthin',
                                  pigment == 'diatox' ~ 'Diatoxanthin',
                                  pigment == 'pheo_b' ~ 'Pheophytin~b'),
         pigment.expr = factor(pigment.expr, levels = c('Diatoxanthin',
                                                        'Alloxanthin',
                                                        'Pheophytin~b',
                                                        'Canthaxanthin',
                                                        'beta-carotene')))

# create plots ----
p.mean <- ggplot() +
  facet_grid(pigment.expr ~ core, scales = 'free_y', labeller = label_parsed) +
  geom_point(aes(year, conc), mb, alpha = 0.3) +
  geom_ribbon(aes(year, ymin = mu.lwr, ymax = mu.upr, fill = core), pred, alpha = 0.3) +
  geom_line(aes(year, mu, color = core), pred) +
  theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 6) +
  scale_fill_brewer(type = 'qual', palette = 6) +
  ylim(c(0, NA)) +
  labs(x = 'Year C.E.', y = expression(Mean~concentration~(nmol~g^{-1}~C)))
p.shape <- ggplot() +
  facet_grid(pigment.expr ~ core, scales = 'free_y', labeller = label_parsed) +
  #geom_ribbon(aes(year, ymin = shape.lwr, ymax = shape.upr, fill = core), pred, alpha = 0.3) +
  geom_line(aes(year, fit.shape, color = core), pred) +
  theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 6) +
  scale_fill_brewer(type = 'qual', palette = 6) +
  ylim(c(0, NA)) +
  labs(x = 'Year C.E.', y = expression(Shape~parameter~(nmol~g^{-1}~C)))
p.var <- ggplot() +
  facet_grid(pigment.expr ~ core, scales = 'free_y', labeller = label_parsed) +
  #geom_ribbon(aes(year, ymin = var.lwr, ymax = var.upr, fill = core), pred, alpha = 0.3) +
  geom_line(aes(year, mu, color = core), pred) +
  theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 6) +
  scale_fill_brewer(type = 'qual', palette = 6) +
  ylim(c(0, NA)) +
  labs(x = 'Year C.E.', y = expression(paste(Concentration~variance~(nmol^2~g^{-2}~C))))

# p2pdf('mb-pigments.pdf', p.mean, width = 5, height = 3.5, scale = 2)
p.full <- plot_grid(p.mean, p.var, ncol = 1, labels = c('a.', 'b.'))
# p2pdf('mb-pigments-mean-variance.pdf', p.full, width = 5, height = 7, scale = 2)
