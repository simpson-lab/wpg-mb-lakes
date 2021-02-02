library('readxl') # for data reading
library('dplyr')  # for data wrangling
library('tidyr')  # for data tidying
source('analysis/default-figure-styling.R')

# import data
mb <- bind_rows(read_xlsx('data/mb/Manitoba pigs isotope Core 1 April 2014.xlsx') %>%
                  select(YEAR, D15N, PERCENTN, D13C, PERCENTC) %>%
                  rename(PERCENT_N = PERCENTN, PERCENT_C = PERCENTC) %>%
                  mutate(CORE = 'Manitoba 1'),
                read_xlsx('data/mb/Manitoba pigs isotope Core 2 April 2014.xlsx') %>%
                  select(CORE, YEAR, D15N, PERCENT_N, D13C, PERCENT_C),
                read_xlsx('data/mb/Manitoba pigs isotope Core 3 April 2014.xlsx') %>%
                  select(CORE, YEAR, D15N, PERCENT_N, D13C, PERCENT_C),
                read_xlsx('data/mb/Manitoba pigs isotope Core 4 April 2014.xlsx') %>%
                  select(CORE, YEAR, D15N, PERCENT_N, D13C, PERCENT_C)) %>%
  mutate(CORE = paste0('Core~', substr(CORE, nchar(CORE), nchar(CORE)))) %>%
  pivot_longer(-c('YEAR', 'CORE')) %>%
  mutate(lab = case_when(name == 'D15N' ~ 'delta^{15}~N~\'\211\'',
                         name == 'PERCENT_N' ~ 'N~\'\045\'',
                         name == 'D13C' ~ 'delta^{13}~C~\'\211\'',
                         name == 'PERCENT_C' ~ 'C~\'\045\'') %>%
           factor(levels = c('delta^{15}~N~\'\211\'', 'delta^{13}~C~\'\211\'',
                             'N~\'\045\'', 'C~\'\045\'')))

# create blank data for D15N plot
b <- tibble(x = 1800,
            y = c(4, 8),
            lab = unique(mb$lab)[1])

p.isotopes <-
  ggplot(filter(mb, CORE == 'Core~1')) +
  facet_grid(lab ~ ., scales = 'free_y', labeller = label_parsed, switch = 'y') +
  geom_point(aes(YEAR, value), alpha = 0.5) +
  geom_blank(aes(x, y), b) + # add blank data to change ylims
  labs(x = 'Year C.E.', y = NULL) +
  theme(strip.background = element_blank(), strip.placement = 'outside')
#p2pdf('isotopes.pdf', p.isotopes, x.plots = 1.5, y.plots = 2, device = 'pdf')

p.isotopes.all <-
  ggplot(mb) +
  facet_grid(lab ~ CORE, scales = 'free_y', labeller = label_parsed, switch = 'y') +
  geom_point(aes(YEAR, value), alpha = 0.5) +
  geom_blank(aes(x, y), b) + # add blank data to change ylims
  labs(x = 'Year C.E.', y = NULL) +
  theme(strip.background = element_blank(), strip.placement = 'outside')
#p2pdf('isotopes-1-4.pdf', p.isotopes.all, x.plots = 3, y.plots = 2, device = 'pdf')
