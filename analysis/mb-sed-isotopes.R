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
                        ncol = 2, byrow = FALSE)

p2pdf('isotopes.pdf', p.isotopes, x.plots = 2, y.plots = 2)
