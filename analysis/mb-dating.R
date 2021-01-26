library('readxl')    # for data reading
library('dplyr')     # for data wrangling
library('tidyr')     # for data wrangling
library('ggplot2')   # for plotting
library('cowplot')   # for plotting
library('extrafont') # for plot fonts
loadfonts(device = 'win', quiet = TRUE)
theme_set(theme_bw(base_family = 'serif', base_line_size = 0.05))

COLNAMES <- read_xls('data/mb/Lake Manitoba dating Core 1.xls',
                     sheet = 'age calculation CRS B',
                     range = 'A11:AD13',
                     col_names = FALSE) %>%
  apply(MARGIN = 2,
        FUN = function(x) paste0('col.', paste(x[!is.na(x)], collapse = ' ')))
cat(COLNAMES, sep = '\n')

read.core.data <- function(core) {
  read_xls(paste0('data/mb/Lake Manitoba dating Core ', core, '.xls'),
           sheet = 'age calculation CRS B',
           range = 'A15:AD29',
           col_names = FALSE) %>%
    transmute(core = paste('Core', core), # fix column names
              depth = ...3,
              pb.210 = ...8,
              se = ...9,
              year = ...30,
              pb.210.lwr = pb.210 - se,  # add 1 SE ranges
              pb.210.upr = pb.210 + se)
}

# name rapair is ok
cores <- bind_rows(read.core.data(1),
                   read.core.data(2),
                   read.core.data(3),
                   read.core.data(4))

p.dating <- 
  plot_grid(ggplot(cores) +
              facet_wrap(. ~ core) +
              geom_point(aes(pb.210, depth)) +
              geom_errorbar(aes(xmin = pb.210.lwr, xmax = pb.210.upr, y = depth)) +
              labs(x = expression(paste(''^{210}~Pb~activity~(Bq~g^{-1}~dry~mass))),
                   y = 'Mid depth (cm)') +
              scale_y_reverse(),
            ggplot(cores) +
              facet_wrap(. ~ core) +
              geom_line(aes(year, depth)) +
              labs(x = 'Estimated year CE',
                   y = 'Mid depth (cm)') +
              scale_y_reverse(),
            labels = c('a.', 'b.'),
            label_fontfamily = 'serif',
            ncol = 1)

ggsave('figures/210Pb-dating.pdf', p.dating, width = 6, height = 8, dpi = 300,
       device = cairo_pdf)
