library('readr')     # for data reading
library('readxl')    # for data reading
library('dplyr')     # for data wrangling
library('tidyr')     # for data wrangling
library('mgcv')      # for modelling
library('gratia')    # for plotting
library('ggplot2')   # for plotting
theme_set(theme_bw())
source('analysis/variance-simulation-and-derivatives.R') # for estimates of V(Y)


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
                    cols = -c('sample', 'depth_cm', 'thickness', 'interval',
                              'year'),
                    names_to = 'pigment',
                    values_to = 'conc') %>%
  mutate(pigment = factor(pigment),
         lake = 'Lake Winnipeg',
         weight = interval / mean(interval)) %>%
  filter(!is.na(conc))

# import Lake Manitoba data ----
mb <- read_xlsx('data/mb/Manitoba pigs isotope Core 1 April 2014.xlsx') %>%
  select(SAMPLE, MID_DEPTH_CM, YEAR, ALLOX, DIATOX, CANTH, PHEO_B, BCAROT,
         PHEO_A, CHLA) %>%
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
                   cols = -c('sample', 'depth_cm', 'thickness', 'interval',
                             'year'),
                   names_to = 'pigment',
                   values_to = 'conc') %>%
  mutate(pigment = factor(pigment),
         lake = 'Lake Manitoba',
         weight = interval / mean(interval)) %>%
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
  geom_point(aes(size = weight)) +
  scale_color_gradient(low = 'darkorange', high = 'blue') +
  scale_size('Weight', range = c(0.75, 1.75)) +
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
                        s(year,lake_pigment, k=15, bs='fs', xt = list(bs='cr')),
                      # formula for scale
                      ~ s(interval, k = 10) + # number of years /sample
                        s(year, lake_pigment, k = 10, bs = 'fs')),
                 family = gammals(),
                 data = lakes,
                 weights = weight,
                 method = 'REML',
                 control = gam.control(nthreads = 4))
tictoc::toc()
#saveRDS(m.gammals, 'models/lakes-gammals-fs.rds')

# check models
lakes <- mutate(lakes,
                e = resid(m.gammals),
                estar = abs(conc - m.gammals$fitted.values[, 1]))
appraise(m.gammals)
draw(m.gammals)

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

newd <- with(lakes,
             expand_grid(year = seq(min(year), max(year), by = 1),
                         pigment = unique(lakes$pigment),
                         lake = unique(lakes$lake))) %>%
  mutate(interval = case_when(lake == 'Lake Manitoba' ~ 4.16,
                              lake == 'Lake Winnipeg' ~ 2.70,
                              TRUE ~ NA_real_),
         lake_pigment = interaction(lake, pigment, drop = TRUE))

sims <-
  gammals_var(model = m.gammals, data = newd, nsims = 100) %>%
  group_by(year, lake, pigment) %>%
  summarize(variance = mean(variance), .groups = 'drop')

ggplot() +
  facet_wrap(lake ~ pigment, ncol = 2, scales = 'free_y', dir = 'h',
             labeller = label_both) +
  geom_point(aes(year, estar), lakes, alpha = 0.3) +
  geom_line(aes(year, sqrt(variance)), sims, color = 'red')
