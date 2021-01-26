library('readxl')    # for data reading
library('dplyr')     # for data wrangling
source('analysis/default-figure-styling.R')

mb <- read_xlsx('data/mb/Manitoba pigs isotope Core 1 April 2014.xlsx') %>%
  select(YEAR, D15N, PERCENTN, D13C, PERCENTC)

p.isotopes <- plot_grid(ggplot(mb) +
                          geom_point(aes(YEAR, D15N), alpha = 0.5) +
                          labs(x = 'Year CE', y = expression(delta^{15}~N)),
                        ggplot(mb) +
                          geom_point(aes(YEAR, PERCENTN), alpha = 0.5) +
                          labs(x = 'Year CE', y = expression('%'~N~g^{-1}~dry~mass)),
                        ggplot(mb) +
                          geom_point(aes(YEAR, D13C), alpha = 0.5) +
                          labs(x = 'Year CE', y = expression(delta^{13}~C)),
                        ggplot(mb) +
                          geom_point(aes(YEAR, PERCENTC), alpha = 0.5) +
                          labs(x = 'Year CE', y = expression('%'~C~g^{-1}~dry~mass)),
                        ncol = 2,
                        byrow = FALSE,
                        labels = LABELS,
                        label_fontfamily = 'serif')

# p2pdf('isotopes.pdf', p.isotopes, x.plots = 2, y.plots = 1.5)

# all four cores
mb.all <- bind_rows(read_xlsx('data/mb/Manitoba pigs isotope Core 1 April 2014.xlsx') %>%
                      select(YEAR, D15N, PERCENTN, D13C, PERCENTC) %>%
                      rename(PERCENT_N = PERCENTN, PERCENT_C = PERCENTC) %>%
                      mutate(CORE = 'Manitoba 1'),
                    read_xlsx('data/mb/Manitoba pigs isotope Core 2 April 2014.xlsx') %>%
                      select(CORE, YEAR, D15N, PERCENT_N, D13C, PERCENT_C),
                    read_xlsx('data/mb/Manitoba pigs isotope Core 3 April 2014.xlsx') %>%
                      select(CORE, YEAR, D15N, PERCENT_N, D13C, PERCENT_C),
                    read_xlsx('data/mb/Manitoba pigs isotope Core 4 April 2014.xlsx') %>%
                      select(CORE, YEAR, D15N, PERCENT_N, D13C, PERCENT_C)) %>%
  mutate(CORE = paste0('Core ', substr(CORE, nchar(CORE), nchar(CORE))))

p.isotopes.all <-
  plot_grid(ggplot(mb.all) +
              facet_grid(CORE ~ .) +
              geom_point(aes(YEAR, D15N), alpha = 0.5) +
              labs(x = 'Year CE', y = expression(delta^{15}~N)),
            ggplot(mb.all) +
              facet_grid(CORE ~ .) +
              geom_point(aes(YEAR, PERCENT_N), alpha = 0.5) +
              labs(x = 'Year CE', y = expression('%'~N~g^{-1}~dry~mass)),
            ggplot(mb.all) +
              facet_grid(CORE ~ .) +
              geom_point(aes(YEAR, D13C), alpha = 0.5) +
              labs(x = 'Year CE', y = expression(delta^{13}~C)),
            ggplot(mb.all) +
              facet_grid(CORE ~ .) +
              geom_point(aes(YEAR, PERCENT_C), alpha = 0.5) +
              labs(x = 'Year CE', y = expression('%'~C~g^{-1}~dry~mass)),
            ncol = 2,
            byrow = FALSE,
            labels = LABELS,
            label_fontfamily = 'serif')

# p2pdf('isotopes-1-4.pdf', p.isotopes.all, y.plots = 2.5, x.plots = 2)
