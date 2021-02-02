library('readxl')    # for data reading
library('dplyr')     # for data wrangling
library('tidyr')     # for data wrangling
library('ggplot2')   # for plotting
library('mgcv')      # for modelling
library('gratia')    # for plotting
theme_set(theme_bw())

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

# no sudden changes in degradation
ggplot(mb) +
  facet_grid(core ~ .) +
  geom_point(aes(year, chla/pheo_a))

mb <- mb %>%
  select(-pheo_a, -chla) %>%
  pivot_longer(cols = -c('core', 'interval', 'year'),
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

m <- readr::read_rds('models/mb-cores-gammals.rds')

# ~ 1.25 minutes on 4 cores
m <- gam(list(conc ~ s(year, core_pigment, bs = 'fs', k = 10), # mean formula
              ~ s(year, core_pigment, bs = 'fs', k = 10) +     # scale formula
                s(interval, bs = 'cr', k = 5)),
         family = gammals(),
         data = mb,
         method = 'REML',
         control = gam.control(nthreads = 4))
#saveRDS(m, 'models/mb-pigments-gammals.rds')

appraise(m) # diagnostics ok
draw(m)
