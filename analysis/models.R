library('readr')   # for data reading
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

# import Lake Winnipeg data
wpg <- read_csv('data/wpg/core-1-summary-data-final.csv',
                col_type = paste(rep('n', 59), collapse = '')) %>%
  select(YEAR, DIATOX, ALLO, PHEO_B, CANTH, B_CAR)

# rename columns to lower case names
colnames(wpg) <- tolower(colnames(wpg))

# change table orientation to long, remove NAs
wpg <- pivot_longer(wpg, cols = -year, names_to = 'pigment',
                    values_to = 'conc') %>%
  mutate(pigment = factor(pigment),
         lake = 'Lake Winnipeg') %>%
  filter(!is.na(conc))

ggplot(wpg, aes(year, conc)) +
  facet_grid(pigment ~ .) +
  geom_point(alpha = 0.5) +
  geom_smooth(aes(color = pigment), se = FALSE, method = 'gam',
              formula = y ~ s(x, k = 10, bs = 'ad'))

# read in models
m.wpg <- readRDS('models/wpg-twlss-fs.rds')

####################################################
#### add lake effect in mean and scale formulae ####
####################################################
m.wpg <- gam(list(conc ~ # formula for mean
                    s(year, k = 20, bs = 'ad') +
                    s(year, pigment, k = 10, bs = 'fs', xt = list(bs = 'cr')),
                  # formula for power
                  ~ 1,
                  # formula for scale
                  ~ 1),
             family = twlss(),
             data = wpg, # change to lakes
             method = 'REML',
             control = gam.control(nthreads = 4, trace = TRUE))
#saveRDS(m.wpg, 'models/wpg-twlss-fs.rds')

# check models
appraise(m.wpg)
draw(m.wpg, scales = 'free', residuals = TRUE)
