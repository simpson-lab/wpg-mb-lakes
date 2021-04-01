library('readxl')    # for data reading
library('dplyr')     # for data wrangling
library('tidyr')     # for data wrangling
library('ggplot2')   # for plotting
library('mgcv')      # for modelling
library('gratia')    # for plotting
theme_set(theme_bw())
source('analysis/variance-simulation-and-derivatives.R')
SAMPLING.DATE <- lubridate::decimal_date(as.POSIXlt('2014-04-01'))

mb <-
  bind_rows(read_xlsx('data/mb/Manitoba pigs isotope Core 1 April 2014.xlsx') %>%
              mutate(core = 'Core~1',
                     interval = c(SAMPLING.DATE, YEAR[-length(YEAR)]) - YEAR,
                     weight = interval / mean(interval)),
            read_xlsx('data/mb/Manitoba pigs isotope Core 2 April 2014.xlsx') %>%
              mutate(core = 'Core~2',
                     CHLA = CHL_A,
                     interval = c(SAMPLING.DATE, YEAR[-length(YEAR)]) - YEAR,
                     weight = interval / mean(interval)),
            read_xlsx('data/mb/Manitoba pigs isotope Core 3 April 2014.xlsx') %>%
              mutate(core = 'Core~3',
                     interval = c(SAMPLING.DATE, YEAR[-length(YEAR)]) - YEAR,
                     weight = interval / mean(interval)),
            read_xlsx('data/mb/Manitoba pigs isotope Core 4 April 2014.xlsx') %>%
              mutate(core = 'Core~4',
                     interval = c(SAMPLING.DATE, YEAR[-length(YEAR)]) - YEAR,
                     weight = interval / mean(interval))) %>%
  select(core, YEAR, ALLOX, DIATOX, CANTH, PHEO_B, BCAROT, PHEO_A, CHLA,
         interval, weight) %>%
  rename(allo = ALLOX,
         b_car = BCAROT)

colnames(mb) <- tolower(colnames(mb))

# check weights
ggplot(mb) +
  facet_grid(core ~ ., labeller = label_parsed) +
  geom_point(aes(year, weight))

# no sudden changes in degradation
ggplot(mb) +
  facet_grid(core ~ ., labeller = label_parsed) +
  geom_point(aes(year, chla/pheo_a))

mb <- mb %>%
  select(-pheo_a, -chla) %>%
  pivot_longer(cols = -c('core', 'interval', 'year', 'weight'),
               names_to = 'pigment',
               values_to = 'conc') %>%
  mutate(pigment = factor(pigment)) %>%
  filter(!is.na(conc)) %>%
  mutate(core_pigment = interaction(core, pigment, sep = '_'))

# similar smoothness between pigments, some different trends between cores
ggplot(mb, aes(year, conc, color = core, fill = core)) +
  facet_grid(pigment ~ ., scales = 'free_y') +
  geom_point(alpha = 0.3) +
  geom_smooth(alpha = 0.3, formula = y ~ s(x), method = 'gam') +
  scale_color_brewer(type = 'qual', palette = 6) +
  scale_fill_brewer(type = 'qual', palette = 6)

# no zeros => can use gamma location-scale model
sum(mb$conc == 0)

m <- readRDS('models/mb-pigments-gammals.rds')

# less than a minute on 4 cores
m <- gam(list(conc ~ s(year, core_pigment, bs = 'fs', k = 10), # mean formula
              ~ s(year, core_pigment, bs = 'fs', k = 10) +     # scale formula
                s(interval, bs = 'cr', k = 5)),
         family = gammals(),
         data = mb,
         method = 'REML',
         weights = weight,
         control = gam.control(nthreads = 4))
#saveRDS(m, 'models/mb-pigments-gammals.rds')

appraise(m) # diagnostics ok
draw(m)

mb <- mutate(mb,
             e = resid(m),
             estar = abs(conc - m$fitted.values[, 1]))

# for some: low concentrations overestimated, high underestimated
cowplot::plot_grid(ggplot(mb, aes(year, e, group = pigment)) +
                     facet_grid(core ~ .) +
                     geom_point(),
                   ggplot(mb, aes(interval, e, group = pigment)) +
                     facet_grid(core ~ .) +
                     geom_point(),
                   ggplot(mb, aes(log(conc), e, group = pigment)) +
                     facet_grid(pigment ~ core, scale = 'free_x') +
                     geom_point(),
                   ggplot(mb, aes(e, group = pigment)) +
                     facet_grid(core ~ .) +
                     geom_density())

newd <- with(mb,
             expand_grid(year = seq(min(year), max(year), by = 1),
                         pigment = unique(mb$pigment),
                         core = unique(mb$core))) %>%
  left_join(group_by(mb, core) %>% summarize(interval = mean(interval)),
            by = 'core') %>%
  mutate(core_pigment = interaction(core, pigment, drop = TRUE, sep = '_'))

sims <-
  gammals_var(model = m, data = newd, nsims = 100) %>%
  group_by(year, core, pigment) %>%
  summarize(variance = mean(variance), .groups = 'drop')

ggplot() +
  facet_grid(pigment ~ core, scales = 'free_y',labeller = label_both)+
  geom_point(aes(year, estar), mb, alpha = 0.3) +
  geom_line(aes(year, sqrt(variance)), sims, color = 'red')
