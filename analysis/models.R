library('readr')   # for data reading
library('readxl')  # for data reading
library('dplyr')   # for data wrangling
library('tidyr')   # for data wrangling
library('mgcv')    # for modelling
library('ggplot2') # for plotting
library('gratia')  # for plotting
theme_set(theme_bw())

# pigments:
# diatox (diatoxanthin  = diatoms)
# allo   (alloxanthin   = cryptophytes)
# pheo_b (pheophytin b  = chlorophytes)
# canth  (canthaxanthin = potentially N2 fixing cyanos)
# b_car  (beta-carotene = all phytoplankton)

# import Lake Winnipeg data ----
# assume interval between mean dating of two slices = interval within 1st slice
wpg <- read_xlsx('data/wpg/Final Core 1 Summary data for Report.xlsx') %>%
  select(SECTION_NO, YEAR, DIATOX, ALLO, PHEO_B, PHEO_A, CHL_A, CANTH, B_CAR, DEPTH_CM) %>%
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

# add thickness and interval of first row
mb <- mb[-(1:2), ] %>%
  select(-pheo_a, -chla)

# change table orientation to long, remove NAs
mb <- pivot_longer(mb,
                   cols = -c('sample', 'depth_cm', 'thickness', 'interval',
                             'year'),
                   names_to = 'pigment',
                   values_to = 'conc') %>%
  mutate(pigment = factor(pigment),
         lake = 'Lake Manitoba') %>%
  filter(!is.na(conc))

# combine lake datasets
lakes <- rbind(wpg, mb) %>%
  mutate(lake_pigment = interaction(lake, pigment, drop = TRUE),
         lake = factor(lake))

# plot data ----
# mb starts in 1801, wpg starts much earlier
ggplot(lakes, aes(year, conc)) +
  facet_wrap(pigment ~ lake, scales = 'free_y', ncol = 2) +
  geom_vline(xintercept = 1800, color = 'red') +
  geom_point() +
  scale_color_gradient(low = 'darkorange', high = 'blue')

# remove data prior to 1800
lakes <- filter(lakes, year >= 1800)

# modelling ----
# read in model
m.gammals <- read_rds('models/lakes-gammals-fs.rds') # gamma location-scale model
m.twlss <- read_rds('models/lakes-twlss-fs.rds')     # tweedie location-scale model

# gamma location-scale model
tictoc::tic()
m.gammals <- gam(list(conc ~ # formula for mean
                        s(year, k = 20) +
                        s(year, lake_pigment, k = 10, bs = 'fs', xt = list(bs = 'cr')),
                      # formula for scale
                      ~ s(interval, k = 10) + # number of years /sample
                        s(year, lake_pigment, k = 10, bs = 'fs')),
                 family = gammals(),
                 data = lakes,
                 method = 'REML',
                 control = gam.control(nthreads = 4))
tictoc::toc()
#saveRDS(m.gammals, 'models/lakes-gammals-fs.rds')

# tweedie location-scale-shape model
# tictoc::tic()
# m.twlss <- gam(list(conc ~ # formula for mean
#                       s(year, k = 20) +
#                       s(year, lake_pigment, k = 10, bs = 'fs', xt = list(bs = 'cr')),
#                     # formula for power
#                     ~ 1, # constant
#                     # formula for scale
#                     ~ s(interval, k = 10) + # number of years /sample
#                       s(year, lake, k = 10, bs = 'fs')),
#                family = twlss(),
#                data = lakes,
#                method = 'REML',
#                control = gam.control(nthreads = 4, maxit = 1e4, trace = FALSE))
# tictoc::toc()
#saveRDS(m.twlss, 'models/lakes-twlss-fs.rds')

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
newd <- with(lakes,
             expand_grid(year = seq(min(year), max(year), by = 1),
                         pigment = unique(lakes$pigment),
                         lake = unique(lakes$lake),
                         interval = 1)) %>%
  mutate(lake_pigment = interaction(lake, pigment, drop = TRUE))

# create predictions
pred <-
  predict(m.gammals, newd, se.fit = TRUE) %>%
  as.data.frame() %>%
  as_tibble() %>%
  rename(logmean = fit.1,
         sd = fit.2,
         se.logmean = se.fit.1,
         se.sd = se.fit.2) %>%
  mutate(mean = exp(logmean),
         mean.lwr = exp(logmean - 1.96 * se.logmean),
         mean.upr = exp(logmean + 1.96 * se.logmean),
         s2 = sd^2) %>%
  bind_cols(newd) %>%
  filter(lake == 'Lake Winnipeg' | (lake == 'Lake Manitoba' & year > 1800)) %>%
  mutate(pigment.expr = case_when(pigment == 'allo' ~ 'Alloxanthin',
                                  pigment == 'b_car' ~ 'beta-carotene',
                                  pigment == 'canth' ~ 'Canthaxanthin',
                                  pigment == 'diatox' ~ 'Diatoxanthin',
                                  pigment == 'pheo_b' ~ 'Pheophytin~b'),
         lake.expr = case_when(lake == 'Lake Manitoba' ~ 'Lake~Manitoba',
                               lake == 'Lake Winnipeg' ~ 'Lake~Winnipeg'))

lakes <-  mutate(lakes,
                 pigment.expr = case_when(pigment == 'allo' ~ 'Alloxanthin',
                                          pigment == 'b_car' ~ 'beta-carotene',
                                          pigment == 'canth' ~ 'Canthaxanthin',
                                          pigment == 'diatox' ~ 'Diatoxanthin',
                                          pigment == 'pheo_b' ~ 'Pheophytin~b'),
                 lake.expr = case_when(lake == 'Lake Manitoba' ~ 'Lake~Manitoba',
                                       lake == 'Lake Winnipeg' ~ 'Lake~Winnipeg'))

# plot predictions
ggplot() +
  facet_wrap(pigment.expr ~ lake.expr, scales = 'free', labeller = label_parsed,
             ncol = 2) +
  geom_ribbon(aes(year, ymin = mean.lwr, ymax = mean.upr), pred, alpha = 0.3) +
  geom_point(aes(year, conc), lakes, alpha = 0.3) +
  geom_line(aes(year, mean), pred) +
  labs(x = NULL, y = 'Mean concentration')

ggplot() +
  facet_grid(pigment ~ lake, scales = 'free') +
#  geom_point(aes(year, estar), lakes, alpha = 0.3) +
  geom_line(aes(year, sd, group = lake_pigment), pred) +
  labs(x = NULL, y = 'Concentration variance')

