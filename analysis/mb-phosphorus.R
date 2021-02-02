library('readxl') # for data reading
library('dplyr')  # for data wrangling
library('tidyr')  # for data tidying
source('analysis/default-figure-styling.R')

phosphorus <- read_xlsx('data/mb/Manitoba sed P Core 1 April 2014.xlsx') %>%
  pivot_longer(cols = -c('CORE', 'DEPTH_CM', 'YEAR'), names_to = 'TYPE',
               values_to = 'p') %>%
  mutate(core = 'Core 1') %>%
  bind_rows(read_xlsx('data/mb/Manitoba sed P Core 2 April 2014.xlsx') %>%
              pivot_longer(cols = -c('CORE', 'DEPTH_CM', 'YEAR'), names_to = 'TYPE',
                           values_to = 'p') %>%
              mutate(core = 'Core 2')) %>%
  mutate(type = case_when(TYPE == 'TP_MG_GDRY' ~ 'Total',
                          TYPE == 'APATITE_P' ~ 'Apatite',
                          TYPE == 'NONAPATITE_P' ~ 'Non-apatite',
                          TYPE == 'EXCHANGEABLE' ~ 'Exchangeable',
                          TYPE == 'ORGANIC' ~ 'Organic') %>%
           factor(levels = c('Total', 'Apatite', 'Non-apatite',
                             'Exchangeable', 'Organic'))) %>%
  rename(year = YEAR) %>%
  filter(!is.na(p))

p1 <- ggplot(phosphorus, aes(year, p, color = core, fill = core, shape = core)) +
  facet_grid(type ~ ., scales = 'free_y') +
  geom_point() +
  geom_smooth(formula = y ~ s(x), method = 'gam', alpha = 0.3, level = 0.95) +
  labs(x = 'Year C.E.', y = expression(Phosphorus~(mg~g^{-1}~dry))) +
  ylim(c(0, NA)) +
  scale_y_continuous(breaks = scales::breaks_pretty(n = 3)) +
  scale_color_brewer(NULL, type = 'qual', palette = 6) +
  scale_fill_brewer(NULL, type = 'qual', palette = 6) +
  scale_shape_discrete(NULL) +
  theme(legend.position = 'top')
p2pdf('mb-phosphorus.pdf', p1, scale = 2)
# p2pdf('phosphorus-2.pdf', p2, y.plots = 2, scale = 2)
