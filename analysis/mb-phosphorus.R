library('readxl') # for data reading
library('dplyr')  # for data wrangling
library('tidyr')  # for data tidying
source('analysis/default-figure-styling.R')

phosphorus <- read_xlsx('data/mb/Manitoba sed P Core 1 April 2014.xlsx') %>%
  pivot_longer(cols = -c('CORE', 'DEPTH_CM', 'YEAR'), names_to = 'TYPE',
               values_to = 'p') %>%
  mutate(type = case_when(TYPE == 'TP_MG_GDRY' ~ 'Total',
                          TYPE == 'APATITE_P' ~ 'Apatite',
                          TYPE == 'NONAPATITE_P' ~ 'Non-apatite',
                          TYPE == 'EXCHANGEABLE' ~ 'Exchangeable',
                          TYPE == 'ORGANIC' ~ 'Organic') %>%
           factor(levels = c('Total', 'Apatite', 'Non-apatite', 'Exchangeable',
                             'Organic')))

foo <- function(x, d = phosphorus, facet = FALSE) {
  p <- 
    ggplot(filter(d, type == x)) +
    geom_point(aes(YEAR, p), alpha = 0.3) +
    labs(x = 'Year C.E.', y = expression(Phosphorus~(mg~g^{-1}~dry)),
         subtitle = x) +
    theme(axis.title = element_blank(), plot.subtitle = element_text(face = 'bold'))
  
  if(facet) {
    p + facet_grid(core ~ .)
  } else {
    p
  }
}

p1 <- plot_grid(get_plot_component(foo('Total') +
                                     theme(axis.title = NULL),
                                   pattern = 'ylab-l'),
                plot_grid(plot_grid(foo('Total'),
                                    foo('Apatite'),
                                    foo('Non-apatite'),
                                    foo('Exchangeable'),
                                    foo('Organic'),
                                    labels = LABELS,
                                    nrow = 3,
                                    label_fontfamily = 'serif'),
                          get_plot_component(foo('Total') +
                                               theme(axis.title = NULL),
                                             pattern = 'xlab-b'),
                          ncol = 1,
                          rel_heights = c(1, 0.05)),
                rel_widths = c(0.05, 1))

# p2pdf('phosphorus-1.pdf', p1, scale = 2)

# add core 2
phosphorus2 <-
  bind_rows(mutate(phosphorus, core = 'Core 1'),
            read_xlsx('data/mb/Manitoba sed P Core 2 April 2014.xlsx') %>%
              pivot_longer(cols = -c('CORE', 'DEPTH_CM', 'YEAR'), names_to = 'TYPE',
                           values_to = 'p') %>%
              mutate(type = case_when(TYPE == 'TP_MG_GDRY' ~ 'Total',
                                      TYPE == 'APATITE_P' ~ 'Apatite',
                                      TYPE == 'NONAPATITE_P' ~ 'Non-apatite',
                                      TYPE == 'EXCHANGEABLE' ~ 'Exchangeable',
                                      TYPE == 'ORGANIC' ~ 'Organic') %>%
                       factor(levels = c('Total', 'Apatite', 'Non-apatite',
                                         'Exchangeable', 'Organic')),
                     core = 'Core 2'))

p2 <- plot_grid(get_plot_component(foo('Total') +
                                     theme(axis.title = NULL),
                                   pattern = 'ylab-l'),
                plot_grid(plot_grid(foo('Total', phosphorus2, facet = TRUE),
                                    foo('Apatite', phosphorus2, facet = TRUE),
                                    foo('Non-apatite', phosphorus2, facet = TRUE),
                                    foo('Exchangeable', phosphorus2, facet = TRUE),
                                    foo('Organic', phosphorus2, facet = TRUE),
                                    labels = LABELS,
                                    nrow = 3,
                                    label_fontfamily = 'serif'),
                          get_plot_component(foo('Total') +
                                               theme(axis.title = NULL),
                                             pattern = 'xlab-b'),
                          ncol = 1,
                          rel_heights = c(1, 0.05)),
                rel_widths = c(0.05, 1))

# p2pdf('phosphorus-2.pdf', p2, y.plots = 2, scale = 2)
