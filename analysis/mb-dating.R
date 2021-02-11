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
core1 <- read.core.data(1)

# fixed coefficients for secondary axis
YEAR.A <- 1740
PB.COEF <- 0.25
CS.COEF <- 1e4 * PB.COEF

core1.tidy <-
  select(core1, depth, pb.210, cs.137) %>%
  pivot_longer(-'depth', values_to = 'est') %>%
  left_join(core1 %>%
              select(depth, se.pb.210, se.cs.137) %>%
              pivot_longer(-'depth', values_to = 'se') %>%
              mutate(name = substr(name, nchar('se.') + 1, nchar(name))),
            by = c('depth', 'name')) %>%
  mutate(slope = case_when(name == 'pb.210' ~ PB.COEF,
                           name == 'cs.137' ~ CS.COEF),
         intercept = case_when(name == 'pb.210' ~ YEAR.A,
                               name == 'cs.137' ~ YEAR.A),
         lwr = (est - se) * slope + intercept,
         upr = (est + se) * slope + intercept,
         est = est * slope + intercept)

isotope.lab <-
  expression(atop(''^{210}~Pb~activity~(Bq~g^{-1}~dry~mass),
                  ''^{137}~Cs~activity~(10^{-4}~Bq~g^{-1}~dry~mass)))
p.triple <-
  ggplot(core1.tidy) +
  
  # year line
  geom_line(aes(year, depth), core1) +
  geom_ribbon(aes(xmin = year.lwr, xmax = year.upr, y = depth),
              core1, alpha = 0.3) +
  
  # points with SE
  geom_errorbarh(aes(xmin = lwr, xmax = upr, y = depth,
                     color = factor(name, levels = c('pb.210','cs.137'))),
                 height = 0) +
  geom_point(aes(est, depth,
                 shape = factor(name, levels = c('pb.210','cs.137')),
                 color = factor(name, levels = c('pb.210','cs.137')))) +
  
  # other
  labs(x = 'Year C. E.',
       y = 'Depth (cm)') +
  scale_y_reverse(name = 'Depth (cm)',
                  breaks = seq(0, 32.5, by = 2.5),
                  labels = c(0, '', 5, '', 10, '', 15, '', 20, '',
                             25, '', 30, ''),
                  limits = c(33, 0)) +
  scale_x_continuous(sec.axis =
                       dup_axis(~ (. - YEAR.A) / PB.COEF,
                                name = isotope.lab,
                                breaks = 0:10 * 100,
                                labels = c(0, '', 200, '', 400, '',
                                           600, '', 800, '', 1000)),
                     breaks = seq(1750, 2000, by = 25),
                     labels = c(1750, '', 1800, '', 1850, '',
                                1900, '', 1950, '', 2000),
                     minor_breaks = seq(1750, 2000, by = 25)) +
  scale_shape_manual(NULL, values = c(19, 1)) +
  scale_color_brewer(NULL, type = 'qual', palette = 6) +
  theme(legend.position = c(0.2, 1.1),
        legend.text = element_blank(),
        legend.key.height = unit(0.275, 'in'),
        legend.spacing.x = unit(-1, 'in'))

#p2pdf('pb-cs-year.pdf', p.triple, scale = 1.5, height = 4)

# multi-panel appendix figure ----
# name rapair is ok
cores <- bind_rows(read.core.data(1),
                   read.core.data(2),
                   read.core.data(3),
                   read.core.data(4))
# y axis for depth
depth.scale <- scale_y_reverse(name = 'Depth (cm)',
                               breaks = seq(0, 37.5, by = 2.5),
                               labels = c(0, '', 5, '', 10, '', 15, '', 20, '',
                                          25, '', 30, '', 35, ''))

cores.tidy <-
  select(cores, core, depth, pb.210, cs.137) %>%
  pivot_longer(-c('depth', 'core'), values_to = 'est') %>%
  left_join(cores %>%
              select(core, depth, se.pb.210, se.cs.137) %>%
              pivot_longer(-c('depth', 'core'), values_to = 'se') %>%
              mutate(name = substr(name, 4, nchar(name))),
            by = c('core', 'depth', 'name')) %>%
  mutate(slope = case_when(name == 'pb.210' ~ 1,
                           name == 'cs.137' ~ 1e4),# don't scale by Pb
         lwr = (est - se) * slope,
         upr = (est + se) * slope,
         est = est * slope)

p.dating <- 
  plot_grid(ggplot(cores.tidy) +
              facet_grid(. ~ core) +
              geom_errorbarh(aes(xmin = lwr, xmax = upr, y = depth,
                                 color = factor(name,
                                                levels = c('pb.210','cs.137'))),
                             height = 0) +
              geom_point(aes(est, depth,
                             shape = factor(name, levels = c('pb.210','cs.137')),
                             color = factor(name, levels = c('pb.210','cs.137')))) +
              depth.scale +
              scale_x_continuous(name = isotope.lab,
                                 breaks = 0:11 * 125,
                                 labels = c(0, '', 250, '', 500, '', 750, '',
                                            1000, '', 1250, '')) +
              scale_shape_manual(NULL, values = c(19, 1)) +
              scale_color_brewer(NULL, type = 'qual', palette = 6) +
              theme(legend.position = c(0.37, -0.3),
                    legend.text = element_blank(),
                    legend.key.height = unit(0.265, 'in'),
                    legend.spacing.y = unit(-1, 'in')),
            ggplot(cores) +
              facet_grid(. ~ core) +
              geom_ribbon(aes(xmin = year.lwr, xmax = year.upr, y = depth),
                          alpha = 0.2) +
              geom_line(aes(year, depth)) +
              scale_x_continuous(name = 'Estimated year C. E.',
                                 breaks = seq(1750, 2000, by = 25),
                                 labels = c(1750, '', 1800, '', 1850, '',
                                            1900, '', 1950, '', 2000),
                                 minor_breaks = seq(1750, 2000, by = 25)) +
              depth.scale,
            labels = LABELS,
            label_fontfamily = 'serif',
            ncol = 1,
            rel_heights = c(1.1, 1))
# p2pdf('dating-appendix.pdf', p.dating, scale = 3)
