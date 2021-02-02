library('readxl') # for data reading
library('dplyr')  # for data wrangling
library('tidyr')  # for data tidying
source('analysis/default-figure-styling.R')

COLNAMES <- read_xls('data/mb/Lake Manitoba dating Core 1.xls',
                     sheet = 'age calculation CRS B',
                     range = 'A11:AD13',
                     col_names = FALSE) %>%
  apply(MARGIN = 2,
        FUN = function(x) paste0('col.', paste(x[!is.na(x)], collapse = ' ')))
t(t(COLNAMES))

read.core.data <- function(core) {
  read_xls(paste0('data/mb/Lake Manitoba dating Core ', core, '.xls'),
           sheet = 'age calculation CRS B',
           range = 'A15:AD29',
           col_names = FALSE) %>%
    transmute(core = paste('Core', core), # fix column names
              depth = ...3,
              pb.210 = ...8,
              se.pb.210 = ...9,
              cs.137 = ...13 / 60, # from decays/minute/g to Bq/g
              se.cs.137 = ...14 / 60, # from decays/minute/g to Bq/g
              year = ...30,
              se.year = ...27, # uncorrelated age SE, ...28 is correlated
              pb.210.lwr = pb.210 - se.pb.210,  # add 1 SE ranges
              pb.210.upr = pb.210 + se.pb.210,
              cs.137.lwr = cs.137 - se.cs.137,
              cs.137.upr = cs.137 + se.cs.137,
              year.lwr = year - se.year,
              year.upr = year + se.year)
}

# single-core plot ----
CS.COEF <- 1e4
YEAR.A <- 1800
YEAR.B <- .2
core1 <- read.core.data(1) %>%
  mutate(pb.year = (year - YEAR.A) / YEAR.B,
         pb.year.lwr = (year - YEAR.A - se.year) / YEAR.B,
         pb.year.upr = (year - YEAR.A + se.year) / YEAR.B)

core1.tidy <-
  select(core1, depth, pb.210, cs.137) %>%
  pivot_longer(-'depth', values_to = 'est') %>%
  left_join(core1 %>%
              select(depth, se.pb.210, se.cs.137) %>%
              pivot_longer(-'depth', values_to = 'se') %>%
              mutate(name = substr(name, 4, nchar(name))),
            by = c('depth', 'name')) %>%
  mutate(slope = case_when(name == 'pb.210' ~ 1,
                           name == 'cs.137' ~ CS.COEF,
                           name == 'year' ~ YEAR.B),
         intercept = case_when(name == 'pb.210' ~ 0,
                               name == 'cs.137' ~ 0,
                               name == 'year' ~ -1800),
         lwr = (est - se) * slope + intercept,
         upr = (est + se) * slope + intercept,
         est = est * slope + intercept)

p.triple <-
  ggplot(filter(core1.tidy, name != 'year')) +
  
  # year line
  geom_line(aes(pb.year, depth), core1) +
  geom_ribbon(aes(xmin = pb.year.lwr, xmax = pb.year.upr, y = depth), core1, alpha = 0.3) +
  
  # points with SE
  geom_errorbar(aes(xmin = lwr, xmax = upr, y = depth,
                    color = factor(name, levels = c('pb.210','cs.137')))) +
  geom_point(aes(est, depth, shape =  factor(name, levels = c('pb.210','cs.137')),
                 color = factor(name, levels = c('pb.210','cs.137')))) +
  
  # other
  labs(x = expression(atop(''^{210}~Pb~activity~(Bq~g^{-1}~dry~mass),
                           ''^{137}~Cs~activity~(10^{-4}~Bq~g^{-1}~dry~mass))),
       y = 'Depth (cm)') +
  scale_y_reverse() +
  scale_x_continuous(sec.axis = dup_axis(~ . * YEAR.B + YEAR.A, name = 'Year C. E.')) +
  scale_shape_manual(NULL, values = c(1, 19)) +
  scale_color_brewer(NULL, type = 'qual', palette = 6) +
  theme(legend.position = c(0.3, -0.13), legend.text = element_blank(),
        legend.key.width = unit(0.05, 'in')); p.triple

# p2pdf('pb-cs-year.pdf', p.triple, scale = 2, height = 2)

# multi-panel appendix figure ----
# name rapair is ok
cores <- bind_rows(read.core.data(1),
                   read.core.data(2),
                   read.core.data(3),
                   read.core.data(4)) %>%
  mutate(pb.year = (year - YEAR.A) / YEAR.B,
         pb.year.lwr = (year - YEAR.A - se.year) / YEAR.B,
         pb.year.upr = (year - YEAR.A + se.year) / YEAR.B)

cores.tidy <-
  select(cores, core, depth, pb.210, cs.137) %>%
  pivot_longer(-c('depth', 'core'), values_to = 'est') %>%
  left_join(cores %>%
              select(core, depth, se.pb.210, se.cs.137) %>%
              pivot_longer(-c('depth', 'core'), values_to = 'se') %>%
              mutate(name = substr(name, 4, nchar(name))),
            by = c('core', 'depth', 'name')) %>%
  mutate(slope = case_when(name == 'pb.210' ~ 1,
                           name == 'cs.137' ~ CS.COEF,
                           name == 'year' ~ YEAR.B),
         intercept = case_when(name == 'pb.210' ~ 0,
                               name == 'cs.137' ~ 0,
                               name == 'year' ~ -1800),
         lwr = (est - se) * slope + intercept,
         upr = (est + se) * slope + intercept,
         est = est * slope + intercept)

p.dating <- 
  plot_grid(ggplot(filter(cores.tidy, name != 'year')) +
              facet_grid(. ~ core) +
              geom_errorbar(aes(xmin = lwr, xmax = upr, y = depth,
                                color = factor(name, levels = c('pb.210','cs.137')))) +
              geom_point(aes(est, depth,
                             shape = factor(name, levels = c('pb.210','cs.137')),
                             color = factor(name, levels = c('pb.210','cs.137')))) +
              labs(x = expression(atop(''^{210}~Pb~activity~(Bq~g^{-1}~dry~mass),
                                       ''^{137}~Cs~activity~(10^{-4}~Bq~g^{-1}~dry~mass))),
                   y = 'Depth (cm)') +
              scale_y_reverse() +
              scale_shape_manual(NULL, values = c(1, 19)) +
              scale_color_brewer(NULL, type = 'qual', palette = 6) +
              theme(legend.position = c(0.37, -0.13), legend.text = element_blank(),
                    legend.key.width = unit(0.05, 'in')),
            ggplot(cores) +
              facet_grid(. ~ core) +
              geom_ribbon(aes(xmin = year.lwr, xmax = year.upr, y = depth), alpha = 0.2) +
              geom_line(aes(year, depth)) +
              labs(x = 'Estimated year C. E.',
                   y = 'Depth (cm)') +
              scale_y_reverse(),
            labels = LABELS,
            label_fontfamily = 'serif',
            ncol = 1,
            rel_heights = c(1.1, 1)); p.dating
# p2pdf('dating-appendix.pdf', p.dating, scale = 3)
